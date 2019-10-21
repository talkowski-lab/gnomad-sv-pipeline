# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# This is an analysis WDL that wraps three steps in the Talkowski SV pipeline:
#  1) minGQ optimization
#  2) minGQ filter application
#  3) post-filter VCF QC

# This is the second build of this workflow, which enumerates many more fine-grained
#  minGQ filtering conditions, and may not be optimized for small cohorts with fewer
#  variants


import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:minGQ_roc_opt_subworkflow/versions/1/plain-WDL/descriptor" as roc_opt_sub
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:master_SV_VCF_QC/versions/81/plain-WDL/descriptor" as QC
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:compute_simple_AFs_singleChrom/versions/14/plain-WDL/descriptor" as calcAF


workflow minGQ_filter_workflow_v2 {
  File vcf
  File vcf_idx
  String prefix
  File contiglist
  File trios_famfile
  String optimize_minSizes
  String optimize_maxSizes
  String optimize_minFreqs
  String optimize_maxFreqs
  String optimize_includeSVTYPEs
  String optimize_includeFILTERs
  String optimize_excludeFILTERs
  String optimize_includeEV
  String optimize_excludeEV
  Int optimize_maxSVperTrio
  Float ROC_maxFDR_PCRMINUS
  Float ROC_maxFDR_PCRPLUS
  Int ROC_minGQ
  Int ROC_maxGQ
  Int ROC_stepGQ
  Int ROC_shards
  Int min_SV_per_proband_per_condition
  Float max_noCallRate
  Int global_minGQ
  String ref_build
  File Sanders_2015_tarball
  File Collins_2017_tarball
  File Werling_2018_tarball
  File PCRPLUS_samples_list

  Array[Array[String]] contigs = read_tsv(contiglist)


  # Get list of PCRMINUS samples
  call get_sample_lists {
    input:
      vcf=vcf,
      vcf_idx=vcf_idx,
      PCRPLUS_samples_list=PCRPLUS_samples_list,
      prefix=prefix
  }

  # Shard VCF per-chromosome and add AF annotation
  scatter ( contig in contigs ) {
    call calcAF.getAFs_singleChrom as getAFs {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        contig=contig[0],
        sv_per_shard=1000,
        prefix=prefix
    }
    #Split VCF into PCR+ and PCR-
    call split_PCR_vcf {
      input:
        vcf=getAFs.vcf_wAFs,
        prefix="${prefix}.${contig[0]}",
        PCRPLUS_samples_list=PCRPLUS_samples_list
    }
    # Annotate PCR-specific AFs
    call calcAF.getAFs_singleChrom as getAFs_byPCR {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        contig=contig[0],
        sv_per_shard=1000,
        prefix=prefix,
        sample_pop_assignments=get_sample_lists.sample_PCR_labels
    }
    # Gather table of AC/AN/AF for PCRPLUS and PCRMINUS samples
    call get_AF_tables {
      input:
        vcf=getAFs_byPCR.vcf_wAFs,
        vcf_idx=getAFs_byPCR.vcf_wAFs_idx,
        prefix="${prefix}.${contig[0]}"
    }
  }


  # Make master table of AC/AN/AF for all variants
  call combine_roc_opt_results as cat_AF_table_PCRPLUS {
    input:
      shards=get_AF_tables.PCRPLUS_AF_table,
      outfile="${prefix}.PCRPLUS.AF_preMinGQ.txt"
  }
  call combine_roc_opt_results as cat_AF_table_PCRMINUS {
    input:
      shards=get_AF_tables.PCRMINUS_AF_table,
      outfile="${prefix}.PCRMINUS.AF_preMinGQ.txt"
  }


  # Get one table of relevant variants & stats per family
  ###PCRPLUS
  call split_famfile as split_famfile_PCRPLUS {
    input:
      vcf=split_PCR_vcf.PCRPLUS_vcf[0],
      vcf_idx=split_PCR_vcf.PCRPLUS_vcf_idx[0],
      famfile=trios_famfile,
      fams_per_shard=1,
      prefix="${prefix}.PCRPLUS"
  }
  scatter ( fam in split_famfile_PCRPLUS.famfile_shards ) {
    call collect_trio_SVdat as collect_trio_SVdat_PCRPLUS {
      input:
        vcf_shards=split_PCR_vcf.PCRPLUS_vcf,
        famfile=fam
    }
  }
  call gather_trio_dat as gather_trio_dat_PCRPLUS {
    input:
      files=collect_trio_SVdat_PCRPLUS.trio_SVdat,
      prefix="${prefix}.PCRPLUS"
  }
  ###PCRMINUS
  call split_famfile as split_famfile_PCRMINUS {
    input:
      vcf=split_PCR_vcf.PCRMINUS_vcf[0],
      vcf_idx=split_PCR_vcf.PCRMINUS_vcf_idx[0],
      famfile=trios_famfile,
      fams_per_shard=1,
      prefix="${prefix}.PCRMINUS"
  }
  scatter ( fam in split_famfile_PCRMINUS.famfile_shards ) {
    call collect_trio_SVdat as collect_trio_SVdat_PCRMINUS {
      input:
        vcf_shards=split_PCR_vcf.PCRMINUS_vcf,
        famfile=fam
    }
  }
  call gather_trio_dat as gather_trio_dat_PCRMINUS {
    input:
      files=collect_trio_SVdat_PCRMINUS.trio_SVdat,
      prefix="${prefix}.PCRMINUS"
  }


  # Get table of all conditions to evaluate
  call enumerate_conditions {
    input:
      prefix=prefix,
      condition_shards=ROC_shards,
      optimize_minSizes=optimize_minSizes,
      optimize_maxSizes=optimize_maxSizes,
      optimize_minFreqs=optimize_minFreqs,
      optimize_maxFreqs=optimize_maxFreqs,
      optimize_includeSVTYPEs=optimize_includeSVTYPEs,
      optimize_includeFILTERs=optimize_includeFILTERs,
      optimize_excludeFILTERs=optimize_excludeFILTERs,
      optimize_includeEV=optimize_includeEV,
      optimize_excludeEV=optimize_excludeEV
  }


  # Scatter over each shard of conditions and send the trio data for ROC optimization
  scatter ( shard in enumerate_conditions.minGQ_conditions_table_noHeader_shards ) {
    ### PCRPLUS
    call roc_opt_sub.minGQ_roc_opt_subworkflow as roc_opt_PCRPLUS {
      input:
        trio_tarball=gather_trio_dat_PCRPLUS.tarball,
        prefix="${prefix}.PCRPLUS",
        trios_list=split_famfile_PCRPLUS.cleaned_trios_famfile,
        conditions_table=shard,
        maxSVperTrio=optimize_maxSVperTrio,
        ROC_maxFDR=ROC_maxFDR_PCRPLUS,
        ROC_minGQ=ROC_minGQ,
        ROC_maxGQ=ROC_maxGQ,
        ROC_stepGQ=ROC_stepGQ,
        min_SV_per_proband_per_condition=min_SV_per_proband_per_condition
    }
    ### PCRMINUS
    call roc_opt_sub.minGQ_roc_opt_subworkflow as roc_opt_PCRMINUS {
      input:
        trio_tarball=gather_trio_dat_PCRMINUS.tarball,
        prefix="${prefix}.PCRMINUS",
        trios_list=split_famfile_PCRMINUS.cleaned_trios_famfile,
        conditions_table=shard,
        maxSVperTrio=optimize_maxSVperTrio,
        ROC_maxFDR=ROC_maxFDR_PCRMINUS,
        ROC_minGQ=ROC_minGQ,
        ROC_maxGQ=ROC_maxGQ,
        ROC_stepGQ=ROC_stepGQ,
        min_SV_per_proband_per_condition=min_SV_per_proband_per_condition
    }
  }


  # Merge ROC results to build minGQ filtering lookup tree
  ###PCRPLUS
  call combine_roc_opt_results as combine_roc_optimal_PCRPLUS {
    input:
      shards=roc_opt_PCRPLUS.ROC_optimal_merged,
      outfile="${prefix}.PCRPLUS.minGQ_condition_opts.txt"
  }
  call combine_roc_opt_results as combine_roc_stats_PCRPLUS {
    input:
      shards=roc_opt_PCRPLUS.distrib_stats_merged,
      outfile="${prefix}.minGQ_condition_distrib_stats.txt"
  }
  ###PCRMINUS
  call combine_roc_opt_results as combine_roc_optimal_PCRMINUS {
    input:
      shards=roc_opt_PCRMINUS.ROC_optimal_merged,
      outfile="${prefix}.PCRMINUS.minGQ_condition_opts.txt"
  }
  call combine_roc_opt_results as combine_roc_stats_PCRMINUS {
    input:
      shards=roc_opt_PCRMINUS.distrib_stats_merged,
      outfile="${prefix}.minGQ_condition_distrib_stats.txt"
  }


  # Create final minGQ filtering tree
  ###PCRPLUS
  call build_filter_tree as build_tree_PCRPLUS {
    input:
      conditions_table=enumerate_conditions.minGQ_conditions_table,
      condition_optimizations=combine_roc_optimal_PCRPLUS.merged_file,
      condition_distrib_stats=combine_roc_stats_PCRPLUS.merged_file,
      prefix="${prefix}.PCRPLUS"
  }
  ###PCRMINUS
  call build_filter_tree as build_tree_PCRMINUS {
    input:
      conditions_table=enumerate_conditions.minGQ_conditions_table,
      condition_optimizations=combine_roc_optimal_PCRMINUS.merged_file,
      condition_distrib_stats=combine_roc_stats_PCRMINUS.merged_file,
      prefix="${prefix}.PCRMINUS"
  }


  # Apply filter per chromosome
  ###PCRPLUS
  scatter ( vcf_shard in split_PCR_vcf.PCRPLUS_vcf ) {
    call apply_minGQ_filter as apply_filter_PCRPLUS {
      input:
        vcf=vcf_shard,
        minGQ_lookup_table=build_tree_PCRPLUS.filter_lookup_table,
        prefix="${prefix}.PCRPLUS",
        PCR_status="PCRPLUS",
        maxNCR=max_noCallRate,
        global_minGQ=global_minGQ
    }
  }
  ###PCRMINUS
  scatter ( vcf_shard in split_PCR_vcf.PCRMINUS_vcf ) {
    call apply_minGQ_filter as apply_filter_PCRMINUS {
      input:
        vcf=vcf_shard,
        minGQ_lookup_table=build_tree_PCRMINUS.filter_lookup_table,
        prefix="${prefix}.PCRMINUS",
        PCR_status="PCRMINUS",
        maxNCR=max_noCallRate,
        global_minGQ=global_minGQ
    }
  }


  # Merge filtered VCFs by PCR status & across chromosomes
  scatter (i in range(length(contigs))) {
    call merge_PCR_VCFs {
      input:
        PCRPLUS_vcf=apply_filter_PCRPLUS.filtered_vcf[i],
        PCRMINUS_vcf=apply_filter_PCRMINUS.filtered_vcf[i],
        prefix=prefix
    }
  }
  call combine_vcfs {
    input:
      vcfs=merge_PCR_VCFs.merged_vcf,
      prefix=prefix
  }


  # Run QC on filtered VCF
  call QC.master_vcf_qc as filtered_VCF_QC {
    input:
      vcf=combine_vcfs.vcf,
      vcf_idx=combine_vcfs.vcf_idx,
      famfile=trios_famfile,
      ref_build=ref_build,
      prefix="${prefix}_minGQ_filtered_VCF",
      sv_per_shard=10000,
      samples_per_shard=50,
      Sanders_2015_tarball=Sanders_2015_tarball,
      Collins_2017_tarball=Collins_2017_tarball,
      Werling_2018_tarball=Werling_2018_tarball,
      contiglist=contiglist
  }


  # Final output
  output {
    File filtered_VCF = combine_vcfs.vcf
    File filtered_VCF_idx = combine_vcfs.vcf_idx
    File filtered_VCF_QC_output = filtered_VCF_QC.sv_vcf_qc_output
    File AF_table_preMinGQ_PCRPLUS = cat_AF_table_PCRPLUS.merged_file
    File AF_table_preMinGQ_PCRMINUS = cat_AF_table_PCRMINUS.merged_file
    # File minGQ_filter_lookup_table = build_filter_tree.filter_lookup_table
    # File minGQ_ordered_tree_hierarchy = build_filter_tree.ordered_tree_hierarchy
  }
}


# Get lists of PCRPLUS and PCRMINUS samples present in input VCF
task get_sample_lists {
  File vcf
  File vcf_idx
  File PCRPLUS_samples_list
  String prefix

  command <<<
    set -euo pipefail
    tabix -H ${vcf} | fgrep -v "##" | cut -f10- | sed 's/\t/\n/g' > all_samples.list
    fgrep -wf ${PCRPLUS_samples_list} all_samples.list > "${prefix}.PCRPLUS.samples.list" || true
    fgrep -wvf ${PCRPLUS_samples_list} all_samples.list > "${prefix}.PCRMINUS.samples.list" || true
    cat \
      <( awk -v OFS="\t" '{ print $1, "PCRPLUS" }' "${prefix}.PCRPLUS.samples.list" || true ) \
      <( awk -v OFS="\t" '{ print $1, "PCRMINUS" }' "${prefix}.PCRMINUS.samples.list" || true ) \
    > "${prefix}.PCR_status_assignments.txt"
  >>>

  output {
    File updated_PCRPLUS_samples_list = "${prefix}.PCRPLUS.samples.list"
    File updated_PCRMINUS_samples_list = "${prefix}.PCRMINUS.samples.list"
    File sample_PCR_labels = "${prefix}.PCR_status_assignments.txt"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:193d18c26100fdd603c569346722513f5796685e990ec3abcaeb4be887062a1a"
    disks: "local-disk 50 HDD"
    preemptible: 1
    maxRetries: 1
  }
}


# Split a VCF into two parts, corresponding to PCR+ and PCR-
task split_PCR_vcf {
  File vcf
  String prefix
  File PCRPLUS_samples_list

  command <<<
    set -euo pipefail
    #Get index of PCR+ samples
    PCRPLUS_idxs=$( zcat ${vcf} | sed -n '1,500p' | fgrep "#" | fgrep -v "##" \
                    | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print NR, $1 }' \
                    | fgrep -wf ${PCRPLUS_samples_list} | cut -f1 | paste -s -d, )
    #Get PCR+ VCF
    zcat ${vcf} \
    | cut -f1-9,"$PCRPLUS_idxs" \
    | bgzip -c \
    > "${prefix}.PCRPLUS.vcf.gz"
    tabix -f "${prefix}.PCRPLUS.vcf.gz"
    #Get PCR- VCF
    zcat ${vcf} \
    | cut --complement -f"$PCRPLUS_idxs" \
    | bgzip -c \
    > "${prefix}.PCRMINUS.vcf.gz"
    tabix -f "${prefix}.PCRMINUS.vcf.gz"
  >>>

  output {
    File PCRPLUS_vcf = "${prefix}.PCRPLUS.vcf.gz"
    File PCRPLUS_vcf_idx = "${prefix}.PCRPLUS.vcf.gz.tbi"
    File PCRMINUS_vcf = "${prefix}.PCRMINUS.vcf.gz"
    File PCRMINUS_vcf_idx = "${prefix}.PCRMINUS.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:193d18c26100fdd603c569346722513f5796685e990ec3abcaeb4be887062a1a"
    disks: "local-disk 50 HDD"
    preemptible: 1
    maxRetries: 1
  }
}


# Get a simple table with ID/AC/AN/AF per variant, prior to minGQ
task get_AF_tables {
  File vcf
  File vcf_idx
  String prefix

  command <<<
    set -euo pipefail
    for PCR in PCRPLUS PCRMINUS; do
      svtk vcf2bed --no-header --no-samples -i "$PCR"_AC -i "$PCR"_AN ${vcf} "$PCR".bed
      awk -v OFS="\t" '{ print $4, $6, $7 }' "$PCR".bed > ${prefix}."$PCR".AF_preMinGQ.txt
    done
  >>>

  output {
    File PCRPLUS_AF_table = "${prefix}.PCRPLUS.AF_preMinGQ.txt"
    File PCRMINUS_AF_table = "${prefix}.PCRMINUS.AF_preMinGQ.txt"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:193d18c26100fdd603c569346722513f5796685e990ec3abcaeb4be887062a1a"
    disks: "local-disk 50 HDD"
    preemptible: 1
    maxRetries: 1
  }
}


# Shard a trio famfile to keep only trios that are all represented in the vcf header
task split_famfile {
  File vcf
  File vcf_idx
  File famfile
  String prefix
  
  Int fams_per_shard

  command <<<
    #Get list of sample IDs & column numbers from VCF header
    tabix -H ${vcf} | fgrep -v "##" | sed 's/\t/\n/g' \
      | awk -v OFS="\t" '{ print $1, NR }' > vcf_header_columns.txt
    #Iterate over families & subset VCF
    while read famID pro fa mo prosex pheno; do
      pro_idx=$( awk -v ID=$pro '{ if ($1==ID) print $2 }' vcf_header_columns.txt )
      fa_idx=$( awk -v ID=$fa '{ if ($1==ID) print $2 }' vcf_header_columns.txt )
      mo_idx=$( awk -v ID=$mo '{ if ($1==ID) print $2 }' vcf_header_columns.txt )
      if ! [ -z $pro_idx ] && ! [ -z $fa_idx ] && ! [ -z $mo_idx ]; then
        fgrep -w "$famID" ${famfile}
      fi
    done < ${famfile} \
    > "${prefix}.cleaned_trios.fam"
    split -l ${fams_per_shard} --numeric-suffixes=00001 -a 5 ${prefix}.cleaned_trios.fam famfile_shard_
  >>>

  output {
    File cleaned_trios_famfile = "${prefix}.cleaned_trios.fam"
    Array[File] famfile_shards = glob("famfile_shard_*")
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:193d18c26100fdd603c569346722513f5796685e990ec3abcaeb4be887062a1a"
    disks: "local-disk 30 HDD"
    preemptible: 1
    maxRetries: 1
  }
}


# Collect a single table of all relevant variants for a single family
task collect_trio_SVdat {
  Array[File] vcf_shards
  File famfile

  command <<<
    for wrapper in 1; do
      #Write header
      echo -e "#famID\tVID\tSVLEN\tAF\tSVTYPE\tFILTER\tpro_EV\tpro_AC\tfa_AC\tmo_AC\tpro_GQ\tfa_GQ\tmo_GQ"
      #Iterate over list of VCF shards
      while read vcf; do
        #Get list of sample IDs & column numbers from VCF header
        zcat "$vcf" | head -n1000 | fgrep "#" | fgrep -v "##" | sed 's/\t/\n/g' \
          | awk -v OFS="\t" '{ print $1, NR }' > vcf_header_columns.txt
        #Iterate over families & subset VCF
        while read famID pro fa mo prosex pheno; do
          pro_idx=$( awk -v ID=$pro '{ if ($1==ID) print $2 }' vcf_header_columns.txt )
          fa_idx=$( awk -v ID=$fa '{ if ($1==ID) print $2 }' vcf_header_columns.txt )
          mo_idx=$( awk -v ID=$mo '{ if ($1==ID) print $2 }' vcf_header_columns.txt )
          if ! [ -z $pro_idx ] && ! [ -z $fa_idx ] && ! [ -z $mo_idx ]; then
            #Subset vcf to only 
            zcat "$vcf" | cut -f1-9,"$pro_idx","$fa_idx","$mo_idx" \
            | grep -e '\#\|[0-1]\/1\|MULTIALLELIC' \
            | bgzip -c > $famID.vcf.gz
            #Get list of CNVs in proband that are ≥5kb have ≥50% coverage in either parent
            svtk vcf2bed -i SVTYPE --no-header $famID.vcf.gz stdout \
            | awk -v OFS="\t" '{ if ($NF ~ /DEL|DUP|MCNV/) print $1, $2, $3, $4, $NF, $6 }' \
            > $famID.CNVs.bed
            fgrep -w $pro $famID.CNVs.bed \
            | awk -v OFS="\t" '{ if ($3-$2>=5000 && $5!="MCNV") print $1, $2, $3, $4, $5 }' \
            > $pro.CNVs.gt5kb.bed
            fgrep -w $fa $famID.CNVs.bed > $fa.CNVs.bed
            fgrep -w $mo $famID.CNVs.bed > $mo.CNVs.bed
            #Deletions
            awk -v OFS="\t" '{ if ($NF=="DEL") print $0, "1" }' $pro.CNVs.gt5kb.bed \
            | bedtools coverage -a - \
              -b <( awk '{ if ($5 ~ /DEL|MCNV/) print $0 }' $fa.CNVs.bed ) \
            | awk -v OFS="\t" '{ if ($NF>=0.5) $NF=1; else $NF=0; print $1, $2, $3, $4, $5, $6, $NF }' \
            | bedtools coverage -a - \
              -b <( awk '{ if ($5 ~ /DEL|MCNV/) print $0 }' $mo.CNVs.bed ) \
            | awk -v OFS="\t" '{ if ($NF>=0.5) $NF=1; else $NF=0; print $4, $6, $7, $NF }' \
            > $famID.RD_genotype_update.txt
            #Duplications
            awk -v OFS="\t" '{ if ($NF=="DUP") print $0, "1" }' $pro.CNVs.gt5kb.bed \
            | bedtools coverage -a - \
              -b <( awk '{ if ($5 ~ /DUP|MCNV/) print $0 }' $fa.CNVs.bed ) \
            | awk -v OFS="\t" '{ if ($NF>=0.5) $NF=1; else $NF=0; print $1, $2, $3, $4, $5, $6, $NF }' \
            | bedtools coverage -a - \
              -b <( awk '{ if ($5 ~ /DUP|MCNV/) print $0 }' $mo.CNVs.bed ) \
            | awk -v OFS="\t" '{ if ($NF>=0.5) $NF=1; else $NF=0; print $4, $6, $7, $NF }' \
            >> $famID.RD_genotype_update.txt
            #Get variant stats
            /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/gather_trio_genos.py \
              --ac-adj $famID.RD_genotype_update.txt \
              --no-header \
              $famID.vcf.gz stdout "$pro" "$fa" "$mo" \
            | awk -v famID="$famID" -v OFS="\t" '{ print famID, $0 }'
          fi
        done < ${famfile}
      done < ${write_lines(vcf_shards)}
    done | bgzip -c > "trio_variant_info.txt.gz"
  >>>

  output {
    File trio_SVdat = "trio_variant_info.txt.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:193d18c26100fdd603c569346722513f5796685e990ec3abcaeb4be887062a1a"
    preemptible: 1
    maxRetries: 2
    memory: "4 GB"
    disks: "local-disk 50 HDD"
  }
}


# Gather all trio SV data into a single tarball (helps with Cromwell file localization)
task gather_trio_dat {
  Array[File] files
  String prefix

  command <<<
    tar -czvf ${prefix}.tar.gz -T ${write_lines(files)}
  >>>

  output {
    File tarball = "${prefix}.tar.gz"
  }

  runtime {
    preemptible: 1
    maxRetries: 1
    docker : "talkowski/sv-pipeline@sha256:193d18c26100fdd603c569346722513f5796685e990ec3abcaeb4be887062a1a"
    disks: "local-disk 250 SSD"
    memory: "4 GB"
  }  
}


# Enumerate all minGQ conditions to test
task enumerate_conditions {
  String prefix
  Int condition_shards
  String optimize_minSizes
  String optimize_maxSizes
  String optimize_minFreqs
  String optimize_maxFreqs
  String optimize_includeSVTYPEs
  String optimize_includeFILTERs
  String optimize_excludeFILTERs
  String optimize_includeEV
  String optimize_excludeEV

  command <<<
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/create_minGQ_tranches_table.R \
      --min.sizes "${optimize_minSizes}" \
      --max.sizes "${optimize_maxSizes}" \
      --min.freqs "${optimize_minFreqs}" \
      --max.freqs "${optimize_maxFreqs}" \
      --svtype.include "${optimize_includeSVTYPEs}" \
      --filter.include "${optimize_includeFILTERs}" \
      --filter.exclude "${optimize_excludeFILTERs}" \
      --ev.include "${optimize_includeEV}" \
      --ev.exclude "${optimize_excludeEV}" \
      "${prefix}.minGQ_conditions.txt"
    fgrep -v "#" "${prefix}.minGQ_conditions.txt" \
      > "${prefix}.minGQ_conditions.noHeader.txt"
    /opt/sv-pipeline/04_variant_resolution/scripts/evenSplitter.R \
      -S ${condition_shards} \
      "${prefix}.minGQ_conditions.noHeader.txt" \
      "${prefix}.minGQ_conditions.noHeader.shard"
  >>>

  output {
    File minGQ_conditions_table = "${prefix}.minGQ_conditions.txt"
    File minGQ_conditions_table_noHeader = "${prefix}.minGQ_conditions.noHeader.txt"
    Array[File] minGQ_conditions_table_noHeader_shards = glob("${prefix}.minGQ_conditions.noHeader.shard*")
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:193d18c26100fdd603c569346722513f5796685e990ec3abcaeb4be887062a1a"
    preemptible: 1
    maxRetries: 1
  }
}


# Merge ROC optimal cutoffs or stats
task combine_roc_opt_results {
  Array[File] shards
  String outfile

  command <<<
    cat ${sep=" " shards} | fgrep -v "#" | sort -Vk1,1 > ${outfile}
  >>>

  output {
    File merged_file = "${outfile}"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:4900cae92f1f8bc98c54f89444a00e134ac4c86ca55543e2646f024270a29a69"
    preemptible: 1
    maxRetries: 1
    memory: "4 GB"
    disks: "local-disk 50 HDD"
  }
}


# Build final minGQ filtering tree
task build_filter_tree {
  File conditions_table
  File condition_optimizations
  File condition_distrib_stats
  String prefix

  command <<<
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/create_minGQ_lookup_table.R \
      "${conditions_table}" \
      "${condition_distrib_stats}" \
      "${condition_optimizations}" \
      "${prefix}.minGQ.ordered_tree_hierarchy.txt" \
      "${prefix}.minGQ.filter_lookup_table.txt"
  >>>


  output {
    File ordered_tree_hierarchy = "${prefix}.minGQ.ordered_tree_hierarchy.txt"
    File filter_lookup_table = "${prefix}.minGQ.filter_lookup_table.txt"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:193d18c26100fdd603c569346722513f5796685e990ec3abcaeb4be887062a1a"
    preemptible: 1
    maxRetries: 1
  }
}


# Apply minGQ filter to VCF
task apply_minGQ_filter {
  File vcf
  File minGQ_lookup_table
  String prefix
  String PCR_status
  Float maxNCR
  Int global_minGQ

  command <<<
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/apply_minGQ_filter.py \
      --minGQ "${global_minGQ}" \
      --maxNCR "${maxNCR}" \
      --cleanAFinfo \
      --prefix "${PCR_status}" \
      "${vcf}" \
      "${minGQ_lookup_table}" \
      stdout \
    | fgrep -v "##INFO=<ID=AN," \
    | fgrep -v "##INFO=<ID=AC," \
    | fgrep -v "##INFO=<ID=AF," \
    | fgrep -v "##INFO=<ID=N_BI_GENOS," \
    | fgrep -v "##INFO=<ID=N_HOMREF," \
    | fgrep -v "##INFO=<ID=N_HET," \
    | fgrep -v "##INFO=<ID=N_HOMALT," \
    | fgrep -v "##INFO=<ID=FREQ_HOMREF," \
    | fgrep -v "##INFO=<ID=FREQ_HET," \
    | fgrep -v "##INFO=<ID=FREQ_HOMALT," \
    | bgzip -c \
    > "${prefix}.minGQ_filtered.vcf.gz"
  >>>

  output {
    File filtered_vcf = "${prefix}.minGQ_filtered.vcf.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:193d18c26100fdd603c569346722513f5796685e990ec3abcaeb4be887062a1a"
    disks: "local-disk 20 SSD"
    memory: "4 GB"
    preemptible: 1
    maxRetries: 1
  }
}


# Merge PCRPLUS and PCRMINUS VCFs for a single chromosome
task merge_PCR_VCFs {
  File PCRPLUS_vcf
  File PCRMINUS_vcf
  String prefix

  command <<<
    #Sanitize FILTER columns
    zcat "${PCRPLUS_vcf}" | cut -f7 | grep -ve '^#' | sed '1d' > PCRPLUS_filters.txt
    zcat "${PCRMINUS_vcf}" | cut -f7 | grep -ve '^#' | sed '1d' > PCRMINUS_filters.txt
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/merge_filter_columns.py \
      PCRPLUS_filters.txt \
      PCRMINUS_filters.txt \
      merged_filters.txt
    #Write new VCF header
    zcat "${PCRPLUS_vcf}" | sed -n '1,1000p' | grep -e '^##' > "${prefix}.minGQ_filtered.vcf"
    zcat "${PCRMINUS_vcf}" | sed -n '1,1000p' | grep -e '^##' | fgrep "NOCALL_RATE" >> "${prefix}.minGQ_filtered.vcf"
    #Column-wise merger
    paste \
      <( zcat "${PCRPLUS_vcf}" | grep -ve '^##' | cut -f1-6 ) \
      <( cat <( echo -e "FILTER" ) merged_filters.txt ) \
      <( zcat "${PCRPLUS_vcf}" | grep -ve '^##' | cut -f8- ) \
      <( zcat "${PCRMINUS_vcf}" | grep -ve '^##' | cut -f10- ) \
    >> "${prefix}.minGQ_filtered.vcf"
    /opt/sv-pipeline/scripts/drop_empty_records.py \
      "${prefix}.minGQ_filtered.vcf" \
      "${prefix}.minGQ_filtered.no_blanks.vcf"
    #Bgzip & tabix
    bgzip -f "${prefix}.minGQ_filtered.no_blanks.vcf"
  >>>

  output {
    File merged_vcf = "${prefix}.minGQ_filtered.no_blanks.vcf.gz"
  }

  runtime {
    preemptible: 1
    maxRetries: 1
    docker: "talkowski/sv-pipeline@sha256:193d18c26100fdd603c569346722513f5796685e990ec3abcaeb4be887062a1a"
    disks: "local-disk 250 SSD"
    memory: "4 GB"
  }
}


# Merge per-chromosome VCF shards
task combine_vcfs {
  Array[File] vcfs
  String prefix
  
  command {
    vcf-concat ${sep=" " vcfs} | vcf-sort | bgzip -c > "${prefix}.minGQ_filtered.vcf.gz";
    tabix -p vcf "${prefix}.minGQ_filtered.vcf.gz"
  }

  runtime {
    preemptible: 0
    maxRetries: 1
    docker: "talkowski/sv-pipeline@sha256:193d18c26100fdd603c569346722513f5796685e990ec3abcaeb4be887062a1a"
    disks: "local-disk 250 SSD"
    memory: "4 GB"
  }

  output {
    File vcf="${prefix}.minGQ_filtered.vcf.gz"
    File vcf_idx="${prefix}.minGQ_filtered.vcf.gz.tbi"
  }
}