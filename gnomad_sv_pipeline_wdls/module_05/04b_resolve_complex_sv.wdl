# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


#Resolve complex SV for a single chromosome
workflow resolve_complex_sv {
  File vcf
  File vcf_idx
  String prefix
  String contig
  Int max_shards_per_chrom
  Int min_variants_per_shard
  File cytobands
  File cytobands_idx
  File mei_bed
  File discfile_list
  File discfile_idx_list
  File pe_blacklist
  File pe_blacklist_idx
  File svc_acct_key
  File rf_cutoffs
  String inv_only

  Array[String] discfiles = read_lines(discfile_list)
  Array[File] discfile_idxs = read_lines(discfile_idx_list)

  # Get SR count cutoff from RF metrics to use in single-ender rescan procedure
  call get_se_cutoff {
    input:
      rf_cutoffs=rf_cutoffs
  }

  #Shard vcf for complex resolution
  #Note: as of Nov 2, 2018, return lists of variant IDs for each shard. This should
  # dramatically improve sharding speed
  call shard_vcf {
    input:
      vcf=vcf,
      vcf_idx=vcf_idx,
      max_shards=max_shards_per_chrom,
      min_variants_per_shard=min_variants_per_shard,
      prefix="${prefix}.${contig}",
      inv_only=inv_only
  }

  #Scatter over shards and resolve variants per shard
  scatter ( VID_list in shard_vcf.VID_lists ) {

    #Prep files for svtk resolve using remote tabixing-enabled pysam
    call resolve_prep {
      input:
        vcf=vcf,
        VIDs_list=VID_list,
        chrom=contig,
        discfiles=discfiles,
        discfile_idxs=discfile_idxs,
        svc_acct_key=svc_acct_key
    }

    #Run svtk resolve
    call svtk_resolve as resolve_cpx_per_shard {
      input:
        full_vcf=resolve_prep.subsetted_vcf,
        noref_vcf=resolve_prep.noref_vcf,
        chrom=contig,
        cytobands=cytobands,
        cytobands_idx=cytobands_idx,
        mei_bed=mei_bed,
        pe_blacklist=pe_blacklist,
        pe_blacklist_idx=pe_blacklist_idx,
        se_pe_cutoff=get_se_cutoff.median_PE_cutoff,
        noref_vids=resolve_prep.noref_vids,
        merged_discfile=resolve_prep.merged_discfile,
        merged_discfile_idx=resolve_prep.merged_discfile_idx
    }

    #Add unresolved variants back into resolved VCF
    call restore_unresolved_cnv as restore_unresolved_cnv_per_shard {
      input:
        resolved_vcf=resolve_cpx_per_shard.rs_vcf,
        unresolved_vcf=resolve_cpx_per_shard.un_vcf,
        chrom=contig
    }
  }

  #Merge across shards
  call concat_vcfs as concat_resolved_per_shard {
    input:
      vcfs=restore_unresolved_cnv_per_shard.res,
      vcftype="resolved",
      prefix=prefix
  }

  output {
    File resolved_vcf_merged = concat_resolved_per_shard.concat_vcf
    File resolved_vcf_merged_idx = concat_resolved_per_shard.concat_vcf_idx
  }
}


#Get SE cutoff: first quartile of PE cutoff from SR random forest across all batches
#Defaults to 4 if first quartile < 4
task get_se_cutoff {
  File rf_cutoffs

  command <<<
    mkdir rf_cutoff_files/
    cat ${rf_cutoffs} | gsutil cp -I rf_cutoff_files/
    while read file; do
      /opt/sv-pipeline/04_variant_resolution/scripts/convert_poisson_p.py \
      $( awk -F '\t' '{if ( $5=="PE_log_pval") print $2 }' $file | head -n1 )
    done < <( find rf_cutoff_files/ -name "*cutoffs" ) | \
    Rscript -e "cat(max(c(4,floor(quantile(as.numeric(scan('stdin',quiet=T)),probs=0.25)))),sep='\n')" > \
    median_cutoff.txt
  >>>

  output {
    Int median_PE_cutoff = read_tsv("median_cutoff.txt")[0][0]
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:40e6a9e956302d32137d0c4a2f33779d23a70ea603d88a8ec03653090bb70107"
    preemptible: 1
  }
}


#Split VCF into chunks for parallelized CPX resolution
task shard_vcf {
  File vcf
  File vcf_idx
  Int max_shards
  Int min_variants_per_shard
  String prefix
  String inv_only

  command <<<
    if [ ${inv_only} == "TRUE" ]; then
      /opt/sv-pipeline/04_variant_resolution/scripts/shardVCF_preResolveCPX_invOnly_part1.sh \
      -L ${min_variants_per_shard} \
      -S ${max_shards} \
      -P ${prefix} \
      -T ${vcf_idx} \
      ${vcf}
    else
      /opt/sv-pipeline/04_variant_resolution/scripts/shardVCF_preResolveCPX_part1.sh \
      -L ${min_variants_per_shard} \
      -S ${max_shards} \
      -P ${prefix} \
      -T ${vcf_idx} \
      ${vcf}
    fi
  >>>

  output {
    Array[File] VID_lists = glob("*.VIDs.list")
  }
  
  runtime {
    preemptible: 1
    docker: "talkowski/sv-pipeline@sha256:ca76dffed573c9c792b8362e594bae23b830045d7ce8585bc90202bdcfe4e9a4"
    disks: "local-disk 500 SSD"
  }
}

#Prep files for svtk resolve
task resolve_prep {
  File vcf
  File VIDs_list
  String chrom
  Array[String] discfiles
  Array[File] discfile_idxs
  File svc_acct_key

  command <<<
    #First, subset VCF to variants of interest
    zcat ${vcf} | sed -n '1,1000p' | fgrep "#" > header.vcf
    zcat ${vcf} | fgrep -v "#" | fgrep -wf ${VIDs_list} | cat header.vcf - | bgzip -c \
    > input.vcf.gz
    #Second, extract all-ref variants from VCF. These break svtk resolve with
    # remote tabixing enabled
    svtk vcf2bed input.vcf.gz input.bed
    fgrep -v "#" input.bed | awk -v FS="\t" '{ if ($6!="") print $4 }' > noref.VIDs.list
    cat <( zcat input.vcf.gz | fgrep "#" ) <( zcat input.vcf.gz | fgrep -wf noref.VIDs.list ) \
      | vcf-sort | bgzip -c \
      > noref.vcf.gz
    #Third, use remote tabix to pull down the discfile chunks within ±2kb of all
    # INVERSION breakpoints, and bgzip / tabix
    echo "LOCALIZING ALL DISCFILE SHARDS..."
    fgrep -v "#" input.bed | fgrep INV | awk -v OFS="\t" -v buffer=2000 \
    '{ print $1, $2-buffer, $2+buffer"\n"$1, $3-buffer, $3+buffer }' \
    | awk -v OFS="\t" '{ if ($2<1) $2=1; print $1, $2, $3 }' \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | bedtools merge -i - \
    > regions_to_tabix.bed
    while read gs_path; do
      gsutil signurl -d 24h ${svc_acct_key} $gs_path | sed '1d' | cut -f 4;
    done < ${write_tsv(discfiles)} > signed_URLs.list
    paste signed_URLs.list ${write_tsv(discfile_idxs)} \
    | awk -v OFS="\t" '{ print $1, $2, "disc"NR"shard" }' > discfiles.list;
    while read url idx slice; do
      echo "REMOTE TABIXING $slice..."
      svtk remote-tabix -R regions_to_tabix.bed "$url" "$idx" \
      | awk '{ if ($1==$4 && $3==$6) print }' \
      | bgzip -c > $slice.txt.gz
    done < discfiles.list;
    #Fourth, merge PE files and add one artificial pair corresponding to the chromosome of interest
    #This makes it so that svtk doesn't break downstream
    zcat disc*shard.txt.gz \
    | cat - \
          <( echo -e "${chrom}\t1\t+\t${chrom}\t2\t+\tDUMMY_SAMPLE_IGNORE" ) \
          <( echo -e "chr${chrom}\t1\t+\t${chrom}\t2\t+\tDUMMY_SAMPLE_IGNORE" ) \
    | sort -Vk1,1 -k2,2n -k5,5n -k7,7 | bgzip -c > discfile.txt.gz;
    rm disc*shard.txt.gz;
    tabix -s 1 -b 2 -e 2 -f discfile.txt.gz;
  >>>

  output {
    File subsetted_vcf = "input.vcf.gz"
    File noref_vcf = "noref.vcf.gz"
    File noref_vids = "noref.VIDs.list"
    File merged_discfile = "discfile.txt.gz"
    File merged_discfile_idx = "discfile.txt.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline-remote-pysam@sha256:76966a6d63b99d98e74eed9b8efb42dd40bcb125c31d22fda679beee801b732d"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
  }
}

#Resolve complex SV
task svtk_resolve {
  File full_vcf
  File noref_vcf
  String chrom
  File cytobands
  File cytobands_idx
  File mei_bed
  File pe_blacklist
  File pe_blacklist_idx
  Int se_pe_cutoff
  File noref_vids
  File merged_discfile
  File merged_discfile_idx

  command <<<
    #Run svtk resolve on variants after all-ref exclusion
    svtk resolve \
      ${noref_vcf} \
      all_batches.resolved.${chrom}.shard.vcf \
      -p AllBatches_CPX_${chrom} \
      -u all_batches.unresolved.${chrom}.shard.vcf \
      --discfile ${merged_discfile} \
      --mei-bed ${mei_bed} \
      --cytobands ${cytobands} \
      --min-rescan-pe-support ${se_pe_cutoff} \
      -x ${pe_blacklist};
    #Add all-ref variants back into resolved VCF
    #Note: requires modifying the INFO field with sed & awk given pysam C bug
    zcat ${full_vcf} | fgrep -v "#" | fgrep -wvf ${noref_vids} \
      | sed -e 's/;MEMBERS=[^\t]*\t/\t/g' \
      | awk -v OFS="\t" '{ $8=$8";MEMBERS="$3; print }' \
      | cat all_batches.resolved.${chrom}.shard.vcf - \
      | vcf-sort \
      > all_batches.resolved.${chrom}.shard.vcf2
    mv all_batches.resolved.${chrom}.shard.vcf2 \
      all_batches.resolved.${chrom}.shard.vcf; 
    #Sanity check for failed svtk jobs
    if [ $( fgrep -v "#" all_batches.resolved.${chrom}.shard.vcf | wc -l ) -eq 0 ] && \
       [ $( fgrep -v "#" all_batches.unresolved.${chrom}.shard.vcf | wc -l ) -eq 0 ]; then
      print "ERROR: BOTH OUTPUT VCFS EMPTY"
      exit 0
    fi
    bgzip -f all_batches.resolved.${chrom}.shard.vcf; 
    bgzip -f all_batches.unresolved.${chrom}.shard.vcf
  >>>

  output {
    File rs_vcf = "all_batches.resolved.${chrom}.shard.vcf.gz"
    File un_vcf = "all_batches.unresolved.${chrom}.shard.vcf.gz"
  }  

  runtime {
    docker: "talkowski/sv-pipeline@sha256:24c91e6eadac380ed4d1ce284e6cedfcc6ac082645e2d1e9e04079795c1d717a"
    preemptible: 1
    memory: "16 GB"
    disks: "local-disk 30 SSD"
  }
}

#Restore unresolved CNVs to resolved VCF
task restore_unresolved_cnv {
  File resolved_vcf
  File unresolved_vcf
  String chrom

  command <<<
    #Add unresolved CNVs to resolved VCF and wipe unresolved status
    zcat ${unresolved_vcf} \
      | fgrep -e "<DEL>" -e "<DUP>" -e "SVTYPE=DEL" -e "SVTYPE=DUP" -e "SVTYPE=CNV" -e "SVTYPE=MCNV" \
      | sed -r -e 's/;EVENT=[^;]*;/;/' -e 's/;UNRESOLVED[^;]*;/;/g' \
      | sed -r -e 's/;UNRESOLVED_TYPE[^;]*;/;/g' -e 's/;UNRESOLVED_TYPE[^\t]*\t/\t/g' \
      | cat <( zcat ${resolved_vcf} ) - \
      > all_batches.resolved_plus_cnv.${chrom}.vcf

    #Add other unresolved variants & retain unresolved status (except for inversion single enders)
    zcat ${unresolved_vcf} \
      | fgrep -v -e "<DEL>" -e "<DUP>" -e "SVTYPE=DEL" -e "SVTYPE=DUP" -e "SVTYPE=CNV" -e "SVTYPE=MCNV" \
      | fgrep -v "INVERSION_SINGLE_ENDER" \
      | fgrep -v "#" \
      >> all_batches.resolved_plus_cnv.${chrom}.vcf

    #Add inversion single enders as SVTYPE=BND
    zcat ${unresolved_vcf} \
      | fgrep -v -e "<DEL>" -e "<DUP>" -e "SVTYPE=DEL" -e "SVTYPE=DUP" -e "SVTYPE=CNV" -e "SVTYPE=MCNV" \
      | fgrep "INVERSION_SINGLE_ENDER" \
      | fgrep -v "#" \
      | sed -e 's/SVTYPE=INV/SVTYPE=BND/g' \
      >> all_batches.resolved_plus_cnv.${chrom}.vcf

    #Sort, clean, and compress
    cat all_batches.resolved_plus_cnv.${chrom}.vcf \
      | vcf-sort -c \
      | /opt/sv-pipeline/04_variant_resolution/scripts/postCPX_cleanup.py \
        /dev/stdin /dev/stdout \
      | bgzip -c \
      > all_batches.resolved_plus_cnv.${chrom}.vcf.gz
  >>>

  output {
    File res = "all_batches.resolved_plus_cnv.${chrom}.vcf.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:40e6a9e956302d32137d0c4a2f33779d23a70ea603d88a8ec03653090bb70107"
    preemptible: 1
    disks: "local-disk 100 SSD"
  }
}

#Combine multiple VCFs
task concat_vcfs {
  Array[File] vcfs
  String vcftype
  String prefix

  command <<<
    vcf-concat ${sep=' ' vcfs} | vcf-sort -c | bgzip -c > ${prefix}.${vcftype}.vcf.gz
    tabix -p vcf -f ${prefix}.${vcftype}.vcf.gz
  >>>

  output {
    File concat_vcf = "${prefix}.${vcftype}.vcf.gz"
    File concat_vcf_idx = "${prefix}.${vcftype}.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:40e6a9e956302d32137d0c4a2f33779d23a70ea603d88a8ec03653090bb70107"
    preemptible: 1
    disks: "local-disk 250 SSD"
  }
}

