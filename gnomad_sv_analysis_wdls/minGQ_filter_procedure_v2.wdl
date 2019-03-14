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


import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:master_SV_VCF_QC/versions/73/plain-WDL/descriptor" as QC
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:compute_simple_AFs_singleChrom/versions/4/plain-WDL/descriptor" as calcAF


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
  Int min_SV_per_proband_per_condition
  Float max_noCallRate
  Int global_minGQ
  String ref_build
  File Sanders_2015_tarball
  File Collins_2017_tarball
  File Werling_2018_tarball
  File PCRPLUS_samples_list

  Array[Array[String]] contigs=read_tsv(contiglist)


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


  # Get table of all conditions to evaluate
  call enumerate_conditions {
    input:
      prefix=prefix,
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


  # Scatter over each condition and send the trio data for ROC optimization
  Array[Array[String]] conditions = read_tsv(enumerate_conditions.minGQ_conditions_table_noHeader)
  scatter ( condition in conditions ) {
    # Subset variants to condition of interest & merge across trios
    # Also computes median & Q1/Q3 variants per proband
    # If median > min_SV_per_proband_per_condition, also runs ROC
    ###PCRPLUS
    call filter_merge_variants_withROC as ROC_PCRPLUS {
      input:
        trio_SVdat=collect_trio_SVdat_PCRPLUS.trio_SVdat,
        prefix="${prefix}.PCRPLUS",
        trios_list=split_famfile_PCRPLUS.cleaned_trios_famfile,
        condition_id=condition[0],
        minSVLEN=condition[1],
        maxSVLEN=condition[2],
        minAF=condition[3],
        maxAF=condition[4],
        includeSVTYPE=condition[5],
        excludeSVTYPE=condition[6],
        includeFILTER=condition[7],
        excludeFILTER=condition[8],
        includeEV=condition[9],
        excludeEV=condition[10],
        maxSVperTrio=optimize_maxSVperTrio,
        ROC_maxFDR=ROC_maxFDR_PCRPLUS,
        ROC_minGQ=ROC_minGQ,
        ROC_maxGQ=ROC_maxGQ,
        ROC_stepGQ=ROC_stepGQ,
        min_SV_per_proband_per_condition=min_SV_per_proband_per_condition
    }
    ###PCRMINUS
    call filter_merge_variants_withROC as ROC_PCRMINUS {
      input:
        trio_SVdat=collect_trio_SVdat_PCRMINUS.trio_SVdat,
        prefix="${prefix}.PCRMINUS",
        trios_list=split_famfile_PCRMINUS.cleaned_trios_famfile,
        condition_id=condition[0],
        minSVLEN=condition[1],
        maxSVLEN=condition[2],
        minAF=condition[3],
        maxAF=condition[4],
        includeSVTYPE=condition[5],
        excludeSVTYPE=condition[6],
        includeFILTER=condition[7],
        excludeFILTER=condition[8],
        includeEV=condition[9],
        excludeEV=condition[10],
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
  call combine_roc_opt_results as combine_roc_PCRPLUS {
    input:
      condition_optimizations=ROC_PCRPLUS.ROC_optimal,
      condition_distrib_stats=ROC_PCRPLUS.distrib_stats,
      prefix="${prefix}.PCRPLUS"
  }
  ###PCRMINUS
  call combine_roc_opt_results as combine_roc_PCRMINUS {
    input:
      condition_optimizations=ROC_PCRMINUS.ROC_optimal,
      condition_distrib_stats=ROC_PCRMINUS.distrib_stats,
      prefix="${prefix}.PCRMINUS"
  }


  # Create final minGQ filtering tree
  ###PCRPLUS
  call build_filter_tree as build_tree_PCRPLUS {
    input:
      conditions_table=enumerate_conditions.minGQ_conditions_table,
      condition_optimizations=combine_roc_PCRPLUS.combined_optimizations,
      condition_distrib_stats=combine_roc_PCRPLUS.combined_distrib_stats,
      prefix="${prefix}.PCRPLUS"
  }
  ###PCRMINUS
  call build_filter_tree as build_tree_PCRMINUS {
    input:
      conditions_table=enumerate_conditions.minGQ_conditions_table,
      condition_optimizations=combine_roc_PCRMINUS.combined_optimizations,
      condition_distrib_stats=combine_roc_PCRMINUS.combined_distrib_stats,
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
    # File minGQ_filter_lookup_table = build_filter_tree.filter_lookup_table
    # File minGQ_ordered_tree_hierarchy = build_filter_tree.ordered_tree_hierarchy
  }
}


# Split a VCF into two parts, corresponding to PCR+ and PCR-
task split_PCR_vcf {
  File vcf
  String prefix
  File PCRPLUS_samples_list

  command <<<
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
    docker: "talkowski/sv-pipeline@sha256:7c08627780b347f52ba49220c289a602984090a019107a1cac5a4158623bf462"
    disks: "local-disk 50 HDD"
    preemptible: 1
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
    docker: "talkowski/sv-pipeline@sha256:7c08627780b347f52ba49220c289a602984090a019107a1cac5a4158623bf462"
    disks: "local-disk 30 HDD"
    preemptible: 1
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
            #Get variant stats
            zcat "$vcf" | cut -f1-9,"$pro_idx","$fa_idx","$mo_idx" \
              | fgrep -v "MULTIALLELIC" \
              | /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/gather_trio_genos.py \
                --no-header stdin stdout "$pro" "$fa" "$mo" \
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
    docker: "talkowski/sv-pipeline@sha256:7c08627780b347f52ba49220c289a602984090a019107a1cac5a4158623bf462"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 50 HDD"
  }
}


# Enumerate all minGQ conditions to test
task enumerate_conditions {
  String prefix
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
  >>>

  output {
    File minGQ_conditions_table = "${prefix}.minGQ_conditions.txt"
    File minGQ_conditions_table_noHeader = "${prefix}.minGQ_conditions.noHeader.txt"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:7c08627780b347f52ba49220c289a602984090a019107a1cac5a4158623bf462"
    preemptible: 1
  }
}


# Subset variants to meet a given set of conditions, merge across trios, 
# and run ROC if condition has enough variants per sample
task filter_merge_variants_withROC {
  Array[File] trio_SVdat
  String prefix
  File trios_list
  String condition_id
  String minSVLEN
  String maxSVLEN
  String minAF
  String maxAF
  String includeSVTYPE
  String excludeSVTYPE
  String includeFILTER
  String excludeFILTER
  String includeEV
  String excludeEV
  Int maxSVperTrio
  Float ROC_maxFDR
  Int ROC_minGQ
  Int ROC_maxGQ
  Int ROC_stepGQ
  Int min_SV_per_proband_per_condition


  command <<<
    #Iterate over families and process them one at a time
    while read famdat; do
      #Subset to variants matching condition
      /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/subset_minGQ_trio_data.R \
        --min.size "${minSVLEN}" \
        --max.size "${maxSVLEN}" \
        --min.freq "${minAF}" \
        --max.freq "${maxAF}" \
        --svtype.include "${includeSVTYPE}" \
        --svtype.exclude "${excludeSVTYPE}" \
        --filter.include "${includeFILTER}" \
        --filter.exclude "${excludeFILTER}" \
        --ev.include "${includeEV}" \
        --ev.exclude "${excludeEV}" \
        --max.variants "${maxSVperTrio}" \
        "$famdat" /dev/stdout
    done < ${write_lines(trio_SVdat)} \
    | gzip -c \
    > "${prefix}.${condition_id}.trio_variants.txt.gz"
    #Compute median # of filtered calls per trio
    if [ $( zcat "${prefix}.${condition_id}.trio_variants.txt.gz" | wc -l ) -gt 0 ]; then
      /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/helper_median_counts_per_trio.R \
        --ID "${condition_id}" \
        "${prefix}.${condition_id}.trio_variants.txt.gz" \
        "${trios_list}" \
        "${prefix}.${condition_id}.perTrio_distrib_stats.txt"
      med=$( fgrep -v "#" "${prefix}.${condition_id}.perTrio_distrib_stats.txt" | cut -f2 )
    else
      echo -e "#condition\thetsPerProband_median\thetsPerProband_Q1\thetsPerProband_Q2\n${condition_id}\t0\t0\t0" \
      > "${prefix}.${condition_id}.perTrio_distrib_stats.txt"
      med=0
    fi
    #Run ROC if enough variants per proband
    echo -e "FINISHED FILTERING. FOUND $med MEDIAN QUALIFYING VARIANTS PER CHILD."
    if [ "$med" -gt ${min_SV_per_proband_per_condition} ]; then
      echo -e "STARTING ROC OPTIMIZATION."
      /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/optimize_minGQ_ROC_v2.R \
        --prefix "${condition_id}" \
        --fdr "${ROC_maxFDR}" \
        --minGQ "${ROC_minGQ}" \
        --maxGQ "${ROC_maxGQ}" \
        --step "${ROC_stepGQ}" \
        "${prefix}.${condition_id}.trio_variants.txt.gz" \
        "${trios_list}" \
        "./"
      gzip "${condition_id}.minGQ_ROC.data.txt"
    else
      echo -e "TOO FEW VARIANTS FOR ROC OPTIMIZATION."
      touch "${condition_id}.minGQ_ROC.data.txt.gz"
      touch "${condition_id}.minGQ_ROC.optimal.txt"
      touch "${condition_id}.minGQ_ROC.plot.pdf"
    fi
  >>>

  output {
    File distrib_stats = "${prefix}.${condition_id}.perTrio_distrib_stats.txt"
    File ROC_data = "${condition_id}.minGQ_ROC.data.txt.gz"
    File ROC_optimal = "${condition_id}.minGQ_ROC.optimal.txt"
    File ROC_plot = "${condition_id}.minGQ_ROC.plot.pdf"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:7c08627780b347f52ba49220c289a602984090a019107a1cac5a4158623bf462"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 50 HDD"
  }
}


# Merge all ROC optimal cutoffs into single file for tree reconstruction
task combine_roc_opt_results {
  Array[File] condition_optimizations
  Array[File] condition_distrib_stats
  String prefix

  command <<<
    find / -name "*.minGQ_ROC.optimal.txt" \
    | xargs -I {} cat {} | fgrep -v "#" | sort -Vk1,1 \
    > "${prefix}.minGQ_condition_opts.txt"
    find / -name "*.perTrio_distrib_stats.txt" \
    | xargs -I {} cat {} | fgrep -v "#" | sort -Vk1,1 \
    > "${prefix}.minGQ_condition_distrib_stats.txt"
  >>>

  output {
    File combined_optimizations = "${prefix}.minGQ_condition_opts.txt"
    File combined_distrib_stats = "${prefix}.minGQ_condition_distrib_stats.txt"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:7c08627780b347f52ba49220c289a602984090a019107a1cac5a4158623bf462"
    preemptible: 1
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
    docker: "talkowski/sv-pipeline@sha256:7c08627780b347f52ba49220c289a602984090a019107a1cac5a4158623bf462"
    preemptible: 1
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
    docker: "talkowski/sv-pipeline@sha256:52c97ee16d3b5c62872bd43bbae75c2b2547d52558e449c808c476e4838f8f7b"
    disks: "local-disk 20 SSD"
    memory: "4 GB"
    preemptible: 1
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
    docker : "talkowski/sv-pipeline@sha256:7c08627780b347f52ba49220c289a602984090a019107a1cac5a4158623bf462"
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
    docker : "talkowski/sv-pipeline@sha256:7c08627780b347f52ba49220c289a602984090a019107a1cac5a4158623bf462"
    disks: "local-disk 250 SSD"
    memory: "4 GB"
  }

  output {
    File vcf="${prefix}.minGQ_filtered.vcf.gz"
    File vcf_idx="${prefix}.minGQ_filtered.vcf.gz.tbi"
  }
}