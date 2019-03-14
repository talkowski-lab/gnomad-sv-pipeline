# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


#This is an analysis WDL that wraps three steps in the Talkowski SV pipeline:
# 1) minGQ optimization
# 2) minGQ filter application
# 3) post-filter VCF QC


import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:optimize_GQ_filter/versions/25/plain-WDL/descriptor" as opt
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:apply_GQ_filter/versions/4/plain-WDL/descriptor" as filt
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:master_SV_VCF_QC/versions/73/plain-WDL/descriptor" as vcf_qc


workflow optimize_GQ_filter {
  File vcf
  File vcf_idx
  String prefix
  File contiglist
  File trios_famfile
  Int fams_per_shard
  Int size_boundary
  Int minGQ_global
  Float max_noCallRate
  Float max_fdr_small_lowSR
  Float max_fdr_small_highSR
  Float max_fdr_large_lowSR
  Float max_fdr_large_highSR
  Int minGQ_ROC
  Int maxGQ_ROC
  Int GQstepsize_ROC
  String ref_build
  File Sanders_2015_tarball
  File Collins_2017_tarball
  File Werling_2018_tarball


  #Optimize minGQ
  call opt.optimize_GQ_filter as ROC_opt {
    input:
      vcf=vcf,
      vcf_idx=vcf_idx,
      prefix=prefix,
      famfile=trios_famfile,
      fams_per_shard=fams_per_shard,
      size_boundary=size_boundary,
      max_fdr_small_lowSR=max_fdr_small_lowSR,
      max_fdr_small_highSR=max_fdr_small_highSR,
      max_fdr_large_lowSR=max_fdr_large_lowSR,
      max_fdr_large_highSR=max_fdr_large_highSR,
      minGQ=minGQ_ROC,
      maxGQ=maxGQ_ROC,
      GQstepsize=GQstepsize_ROC,
      contiglist=contiglist
  }


  #Apply minGQ filter to small lowSR variants
  call filt.apply_GQ_filter as small_lowSR_filt {
    input:
      vcf=ROC_opt.small_lowSR_vcf_subset,
      vcf_idx=ROC_opt.small_lowSR_vcf_subset_idx,
      prefix="${prefix}.small_lowSR",
      contiglist=contiglist,
      minGQ_global=minGQ_global,
      minGQ_perSVTYPE_table=ROC_opt.small_lowSR_minGQ_ROC_table,
      max_noCallRate=max_noCallRate
  }


  #Apply minGQ filter to small highSR variants
  call filt.apply_GQ_filter as small_highSR_filt {
    input:
      vcf=ROC_opt.small_highSR_vcf_subset,
      vcf_idx=ROC_opt.small_highSR_vcf_subset_idx,
      prefix="${prefix}.small_highSR",
      contiglist=contiglist,
      minGQ_global=minGQ_global,
      minGQ_perSVTYPE_table=ROC_opt.small_highSR_minGQ_ROC_table,
      max_noCallRate=max_noCallRate
  }


  #Apply minGQ filter to large lowSR variants
  call filt.apply_GQ_filter as large_lowSR_filt {
    input:
      vcf=ROC_opt.large_lowSR_vcf_subset,
      vcf_idx=ROC_opt.large_lowSR_vcf_subset_idx,
      prefix="${prefix}.large_lowSR",
      contiglist=contiglist,
      minGQ_global=minGQ_global,
      minGQ_perSVTYPE_table=ROC_opt.large_lowSR_minGQ_ROC_table,
      max_noCallRate=max_noCallRate
  }


  #Apply minGQ filter to large highSR variants
  call filt.apply_GQ_filter as large_highSR_filt {
    input:
      vcf=ROC_opt.large_highSR_vcf_subset,
      vcf_idx=ROC_opt.large_highSR_vcf_subset_idx,
      prefix="${prefix}.large_highSR",
      contiglist=contiglist,
      minGQ_global=minGQ_global,
      minGQ_perSVTYPE_table=ROC_opt.large_highSR_minGQ_ROC_table,
      max_noCallRate=max_noCallRate
  }


  #Merge filtered VCFs
  call concat_vcfs {
    input:
      small_lowSR_vcf=small_lowSR_filt.filtered_vcf,
      small_highSR_vcf=small_highSR_filt.filtered_vcf,
      large_lowSR_vcf=large_lowSR_filt.filtered_vcf,
      large_highSR_vcf=large_highSR_filt.filtered_vcf,
      outfile_prefix="${prefix}.miGQ_filtered"
  }


  #Merge ROC optimization plots & tables
  call gather_ROC_opt_summary_dat {
    input:
      prefix=prefix,
      small_lowSR_plot=ROC_opt.small_lowSR_minGQ_ROC_plot,
      small_lowSR_table=ROC_opt.small_lowSR_minGQ_ROC_table,
      small_highSR_plot=ROC_opt.small_highSR_minGQ_ROC_plot,
      small_highSR_table=ROC_opt.small_highSR_minGQ_ROC_table,
      large_lowSR_plot=ROC_opt.large_lowSR_minGQ_ROC_plot,
      large_lowSR_table=ROC_opt.large_lowSR_minGQ_ROC_table,
      large_highSR_plot=ROC_opt.large_highSR_minGQ_ROC_plot,
      large_highSR_table=ROC_opt.large_highSR_minGQ_ROC_table
  }


  #Run QC on filtered & merged VCF
  call vcf_qc.master_vcf_qc as filtered_vcf_qc {
    input:
      vcf=concat_vcfs.concat_vcf,
      vcf_idx=concat_vcfs.concat_vcf_idx,
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


  #Final outputs
  output {
    File filtered_vcf = concat_vcfs.concat_vcf
    File filtered_vcf_idx = concat_vcfs.concat_vcf_idx
    File filtered_vcf_QC = filtered_vcf_qc.sv_vcf_qc_output
    File minGQ_optimization_summary_data = gather_ROC_opt_summary_dat.tarball
  }
}


#Task to combine filtered VCFs
task concat_vcfs {
  File small_lowSR_vcf
  File small_highSR_vcf
  File large_lowSR_vcf
  File large_highSR_vcf
  String outfile_prefix

  command <<<
    vcf-concat ${small_lowSR_vcf} ${small_highSR_vcf} ${large_lowSR_vcf} ${large_highSR_vcf} \
    | vcf-sort -c | bgzip -c > ${outfile_prefix}.vcf.gz; 
    tabix -p vcf -f "${outfile_prefix}.vcf.gz"
  >>>

  output {
    File concat_vcf = "${outfile_prefix}.vcf.gz"
    File concat_vcf_idx = "${outfile_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:d3844f6c7c26da55e679c9c521882d54dbecd169f884f09f05a12d6565bf6063"
    preemptible: 1
    memory: "8 GB"
    disks: "local-disk 1000 SSD"
  }
}


#Task to collect ROC optimization summary data
task gather_ROC_opt_summary_dat {
  String prefix
  File small_lowSR_plot
  File small_lowSR_table
  File small_highSR_plot
  File small_highSR_table
  File large_lowSR_plot
  File large_lowSR_table
  File large_highSR_plot
  File large_highSR_table

  command <<<
    mkdir "${prefix}_minGQ_ROC_optimization_summary_data"
    mv ${small_lowSR_plot} "${prefix}_minGQ_ROC_optimization_summary_data"/
    mv ${small_lowSR_table} "${prefix}_minGQ_ROC_optimization_summary_data"/
    mv ${small_highSR_plot} "${prefix}_minGQ_ROC_optimization_summary_data"/
    mv ${small_highSR_table} "${prefix}_minGQ_ROC_optimization_summary_data"/
    mv ${large_lowSR_plot} "${prefix}_minGQ_ROC_optimization_summary_data"/
    mv ${large_lowSR_table} "${prefix}_minGQ_ROC_optimization_summary_data"/
    mv ${large_highSR_plot} "${prefix}_minGQ_ROC_optimization_summary_data"/
    mv ${large_highSR_table} "${prefix}_minGQ_ROC_optimization_summary_data"/
    tar -czvf "${prefix}_minGQ_ROC_optimization_summary_data.tar.gz" \
      "${prefix}_minGQ_ROC_optimization_summary_data"
  >>>

  output {
    File tarball = "${prefix}_minGQ_ROC_optimization_summary_data.tar.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:d3844f6c7c26da55e679c9c521882d54dbecd169f884f09f05a12d6565bf6063"
    preemptible: 1
    disks: "local-disk 30 SSD"
  }
}
