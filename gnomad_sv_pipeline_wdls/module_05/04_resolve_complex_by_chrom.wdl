workflow resolve_complex_by_chrom {
  File vcf
  File vcf_idx
  String contig
  Int max_shards
  Int min_variants_per_shard
  File cytobands
  File cytobands_idx
  File mei_bed
  File discfile_list
  File discfile_idx_list
  File pe_blacklist
  File pe_blacklist_idx
  File svc_acct_key
  Int se_pe_cutoff
  
  Array[String] discfiles = read_lines(discfile_list)
  Array[File] discfile_idxs = read_lines(discfile_idx_list)

  call shard_vcf {
    input:
      vcf=vcf,
      vcf_idx=vcf_idx,
      max_shards=max_shards,
      min_per_shard=min_variants_per_shard,
      prefix=contig
  }

  scatter (vcf_shard in shard_vcf.vcf_shards) {
    call svtk_resolve as resolve_cpx_per_shard {
      input:
        vcf=vcf_shard,
        chrom=contig,
        cytobands=cytobands,
        cytobands_idx=cytobands_idx,
        discfiles=discfiles,
        discfile_idxs=discfile_idxs,
        mei_bed=mei_bed,
        pe_blacklist=pe_blacklist,
        pe_blacklist_idx=pe_blacklist_idx,
        svc_acct_key=svc_acct_key,
        se_pe_cutoff=se_pe_cutoff
    }

    call restore_unresolved_cnv as restore_unresolved_cnv_per_shard {
      input:
        resolved_vcf=resolve_cpx_per_shard.rs_vcf,
        unresolved_vcf=resolve_cpx_per_shard.un_vcf,
        chrom=contig
    }
  }

  call concat_vcfs as concat_resolved_per_shard {
    input:
      vcfs=restore_unresolved_cnv_per_shard.res,
      vcftype="resolved"
  }

  # call concat_vcfs as concat_unresolved_per_shard {
  #   input:
  #     vcfs=restore_unresolved_cnv_per_shard.unres,
  #     vcftype="unresolved"
  # }

  output {
  File res_vcf = concat_resolved_per_shard.concat_vcf
  # File unres_vcf = concat_unresolved_per_shard.concat_vcf
  }
}

#Split VCF into chunks for parallelized CPX resolution
task shard_vcf {
  File vcf
  File vcf_idx
  Int max_shards
  Int min_per_shard
  String prefix

  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/shardVCF_preResolveCPX.sh \
    -L ${min_per_shard} \
    -S ${max_shards} \
    -P ${prefix} \
    -T ${vcf_idx} \
    ${vcf}
  >>>

  output {
    Array[File] vcf_shards = glob("${prefix}.shard_*.vcf.gz")
  }
  
  runtime {
    preemptible: 3
    docker: "talkowski/sv-pipeline@sha256:b359f2cb0c9d5f5a55eb4c41fd362f4e574bf3f8f0f395a2907837571b367ee0"
    disks: "local-disk 1000 SSD"
  }
}

#Resolve complex SV
task svtk_resolve {
  File vcf
  String chrom
  File cytobands
  File cytobands_idx
  File mei_bed
  Array[String] discfiles
  Array[File] discfile_idxs
  File svc_acct_key
  File pe_blacklist
  File pe_blacklist_idx
  Int se_pe_cutoff

  command <<<
    #First, extract all-ref variants from VCF. These break svtk resolve with
    # remote tabixing enabled
    svtk vcf2bed ${vcf} input.bed
    awk -v FS="\t" '{ if ($6!="") print $4 }' input.bed > noref.VIDs.list
    cat <( zcat ${vcf} | fgrep "#" ) <( zcat ${vcf} | fgrep -wf noref.VIDs.list ) \
      | vcf-sort | bgzip -c \
      > noref.vcf.gz
    #Run svtk resolve on variants after all-ref exclusion
    while read gs_path; do
      gsutil signurl -d 24h ${svc_acct_key} $gs_path | sed '1d' | cut -f 4;
    done < ${write_tsv(discfiles)} > signed_URLs.list;
    paste signed_URLs.list ${write_tsv(discfile_idxs)} > discfiles.list
    svtk resolve \
      noref.vcf.gz \
      all_batches.resolved.${chrom}.shard.vcf \
      -p AllBatches_CPX_${chrom} \
      -u all_batches.unresolved.${chrom}.shard.vcf \
      --discfile-list discfiles.list \
      --mei-bed ${mei_bed} \
      --cytobands ${cytobands} \
      --min-rescan-pe-support ${se_pe_cutoff} \
      -x ${pe_blacklist}
    #Add all-ref variants back into resolved VCF
    #Note: requires modifying the INFO field with sed & awk given pysam C bug
    zcat ${vcf} | fgrep -v "#" | fgrep -wvf noref.VIDs.list \
      | sed -e 's/;MEMBERS=[^\t]*\t/\t/g' \
      | awk -v OFS="\t" '{ $8=$8";MEMBERS="$3; print }' \
      | cat all_batches.resolved.${chrom}.shard.vcf - \
      | vcf-sort \
      > all_batches.resolved.${chrom}.shard.vcf2
    mv all_batches.resolved.${chrom}.shard.vcf2 \
      all_batches.resolved.${chrom}.shard.vcf
  >>>

  output {
    File rs_vcf = "all_batches.resolved.${chrom}.shard.vcf"
    File un_vcf = "all_batches.unresolved.${chrom}.shard.vcf"
  }  

  runtime {
    docker: "talkowski/sv-pipeline-remote-pysam@sha256:0c21137179665254ca0d9ebe4d21251ae2ff6679337fd9b3e9d6e6ab808db6a8"
    preemptible: 3
    memory: "8 GB"
    disks: "local-disk 250 SSD"
  }
}

#Restore unresolved CNVs to resolved VCF
task restore_unresolved_cnv {
  File resolved_vcf
  File unresolved_vcf
  String chrom

  command <<<
    #Add unresolved CNVs to resolved VCF and wipe unresolved status
    fgrep -e "<DEL>" -e "<DUP>" -e "SVTYPE=DEL" -e "SVTYPE=DUP" -e "SVTYPE=CNV" -e "SVTYPE=MCNV" ${unresolved_vcf} \
      | sed -r -e 's/;EVENT=[^;]*;/;/' -e 's/;UNRESOLVED[^;]*;/;/g' \
      | sed -r -e 's/;UNRESOLVED_TYPE[^;]*;/;/g' -e 's/;UNRESOLVED_TYPE[^\t]*\t/\t/g' \
      | cat ${resolved_vcf} - \
      > all_batches.resolved_plus_cnv.${chrom}.vcf

    #Add other unresolved variants & retain unresolved status (except for inversion single enders)
    fgrep -v -e "<DEL>" -e "<DUP>" -e "SVTYPE=DEL" -e "SVTYPE=DUP" -e "SVTYPE=CNV" -e "SVTYPE=MCNV" ${unresolved_vcf} \
      | fgrep -v "INVERSION_SINGLE_ENDER" \
      | fgrep -v "#" \
      >> all_batches.resolved_plus_cnv.${chrom}.vcf

    #Add inversion single enders as SVTYPE=BND
    fgrep -v -e "<DEL>" -e "<DUP>" -e "SVTYPE=DEL" -e "SVTYPE=DUP" -e "SVTYPE=CNV" -e "SVTYPE=MCNV" ${unresolved_vcf} \
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
    docker: "talkowski/sv-pipeline@sha256:b359f2cb0c9d5f5a55eb4c41fd362f4e574bf3f8f0f395a2907837571b367ee0"
    preemptible: 1
    disks: "local-disk 1000 SSD"
  }
}

#Combine multiple VCFs
task concat_vcfs {
  Array[File] vcfs
  String vcftype

  command <<<
    vcf-concat ${sep=' ' vcfs} | vcf-sort -c | bgzip -c > all_batches.${vcftype}.vcf.gz
    tabix -p vcf -f all_batches.${vcftype}.vcf.gz
  >>>

  output {
    File concat_vcf = "all_batches.${vcftype}.vcf.gz"
    File concat_vcf_idx = "all_batches.${vcftype}.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:b359f2cb0c9d5f5a55eb4c41fd362f4e574bf3f8f0f395a2907837571b367ee0"
    preemptible: 1
    disks: "local-disk 1000 SSD"
  }
}
