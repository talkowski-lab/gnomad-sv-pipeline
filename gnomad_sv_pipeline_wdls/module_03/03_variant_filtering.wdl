import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:03_filter_vcf/versions/10/plain-WDL/descriptor" as filter_vcf
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:03_filter_outliers/versions/6/plain-WDL/descriptor" as filter_outliers
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:master_SV_VCF_QC/versions/47/plain-WDL/descriptor" as vcf_qc

workflow variant_filtering {
  File evidence_metrics
  File manta_vcf
  File delly_vcf
  File melt_vcf
  File depth_vcf
  String batch
  File famfile
  File trios_famfile
  String ref_build
  File Sanders_2015_tarball
  File Collins_2017_tarball
  File Werling_2018_tarball
  Int outlier_cutoff_nIQR
  Array[String] samples

  call adjudicate_SV {
    input:
      metrics=evidence_metrics,
      batch=batch
  }

  call filter_vcf.RF_filter_vcf as filter_delly {
    input:
      vcf=delly_vcf,
      metrics=evidence_metrics,
      scores=adjudicate_SV.scores,
      cutoffs=adjudicate_SV.cutoffs,
      prefix="${batch}.delly"
  }
  
  call filter_vcf.RF_filter_vcf as filter_manta {
    input:
      vcf=manta_vcf,
      metrics=evidence_metrics,
      scores=adjudicate_SV.scores,
      cutoffs=adjudicate_SV.cutoffs,
      prefix="${batch}.manta"
  }
  
  call filter_vcf.RF_filter_vcf as filter_melt {
    input:
      vcf=melt_vcf,
      metrics=evidence_metrics,
      scores=adjudicate_SV.scores,
      cutoffs=adjudicate_SV.cutoffs,
      prefix="${batch}.melt"
  }
  
  call filter_vcf.RF_filter_vcf as filter_depth {
    input:
      vcf=depth_vcf,
      metrics=evidence_metrics,
      scores=adjudicate_SV.scores,
      cutoffs=adjudicate_SV.cutoffs,
      prefix="${batch}.depth"
  }

  call filter_outliers.filter_outlier_samples as outlier_exclusion {
    input:
      delly_vcf=filter_delly.filtered_vcf,
      manta_vcf=filter_manta.filtered_vcf,
      melt_vcf=filter_melt.filtered_vcf,
      depth_vcf=filter_depth.filtered_vcf,
      N_IQR_cutoff=outlier_cutoff_nIQR,
      batch=batch,
      samples=samples
  }

  call update_famfile {
    input:
      famfile=famfile,
      excluded_samples=outlier_exclusion.outlier_samples_excluded,
      batch=batch
  }

  call merge_pesr_vcfs {
    input:
      delly_vcf=outlier_exclusion.delly_vcf_noOutliers,
      manta_vcf=outlier_exclusion.manta_vcf_noOutliers,
      melt_vcf=outlier_exclusion.melt_vcf_noOutliers,
      batch=batch
  }

  call vcf_qc.master_vcf_qc as delly_qc {
    input:
      vcf=outlier_exclusion.delly_vcf_noOutliers,
      famfile=trios_famfile,
      ref_build=ref_build,
      prefix="${batch}.delly_03_filtered_vcf",
      sv_per_shard=10000,
      samples_per_shard=100,
      Sanders_2015_tarball=Sanders_2015_tarball,
      Collins_2017_tarball=Collins_2017_tarball,
      Werling_2018_tarball=Werling_2018_tarball
  }

  call vcf_qc.master_vcf_qc as manta_qc {
    input:
      vcf=outlier_exclusion.manta_vcf_noOutliers,
      famfile=trios_famfile,
      ref_build=ref_build,
      prefix="${batch}.manta_03_filtered_vcf",
      sv_per_shard=10000,
      samples_per_shard=100,
      Sanders_2015_tarball=Sanders_2015_tarball,
      Collins_2017_tarball=Collins_2017_tarball,
      Werling_2018_tarball=Werling_2018_tarball
  }

  call vcf_qc.master_vcf_qc as melt_qc {
    input:
      vcf=outlier_exclusion.melt_vcf_noOutliers,
      famfile=trios_famfile,
      ref_build=ref_build,
      prefix="${batch}.melt_03_filtered_vcf",
      sv_per_shard=10000,
      samples_per_shard=100,
      Sanders_2015_tarball=Sanders_2015_tarball,
      Collins_2017_tarball=Collins_2017_tarball,
      Werling_2018_tarball=Werling_2018_tarball
  }

  call vcf_qc.master_vcf_qc as pesr_qc {
    input:
      vcf=merge_pesr_vcfs.merged_pesr_vcf,
      famfile=trios_famfile,
      ref_build=ref_build,
      prefix="${batch}.pesr_merged_03_filtered_vcf",
      sv_per_shard=10000,
      samples_per_shard=100,
      Sanders_2015_tarball=Sanders_2015_tarball,
      Collins_2017_tarball=Collins_2017_tarball,
      Werling_2018_tarball=Werling_2018_tarball
  }

  call vcf_qc.master_vcf_qc as depth_qc {
    input:
      vcf=outlier_exclusion.depth_vcf_noOutliers,
      famfile=trios_famfile,
      ref_build=ref_build,
      prefix="${batch}.depth_03_filtered_vcf",
      sv_per_shard=10000,
      samples_per_shard=100,
      Sanders_2015_tarball=Sanders_2015_tarball,
      Collins_2017_tarball=Collins_2017_tarball,
      Werling_2018_tarball=Werling_2018_tarball
  }

  output {
    File filtered_delly_vcf = outlier_exclusion.delly_vcf_noOutliers
    File filtered_delly_vcf_qc = delly_qc.sv_vcf_qc_output
    File filtered_manta_vcf = outlier_exclusion.manta_vcf_noOutliers
    File filtered_manta_vcf_qc = manta_qc.sv_vcf_qc_output
    File filtered_melt_vcf = outlier_exclusion.melt_vcf_noOutliers
    File filtered_melt_vcf_qc = melt_qc.sv_vcf_qc_output
    File filtered_pesr_vcf = merge_pesr_vcfs.merged_pesr_vcf
    File filtered_pesr_vcf_qc = pesr_qc.sv_vcf_qc_output
    File filtered_depth_vcf = outlier_exclusion.depth_vcf_noOutliers
    File filtered_depth_vcf_qc = depth_qc.sv_vcf_qc_output
    File scores = adjudicate_SV.scores
    File cutoffs = adjudicate_SV.cutoffs
    File RF_intermediate_files = adjudicate_SV.RF_intermediate_files
    File outlier_samples_excluded = outlier_exclusion.outlier_samples_excluded
    File batch_samples_postOutlierExclusion = outlier_exclusion.filtered_batch_samples_list
    File famfile_postOutlierExclusion = update_famfile.filtered_famfile
  }
}


# Adjudicate SV
task adjudicate_SV {
  File metrics
  String batch

  command <<<
    svtk adjudicate ${metrics} ${batch}.scores ${batch}.cutoffs
    mkdir ${batch}.RF_intermediate_files
    mv *_trainable.txt ${batch}.RF_intermediate_files/
    mv *_testable.txt ${batch}.RF_intermediate_files/
    tar -czvf ${batch}.RF_intermediate_files.tar.gz ${batch}.RF_intermediate_files
  >>>

  output {
    File scores = "${batch}.scores"
    File cutoffs = "${batch}.cutoffs"
    File RF_intermediate_files = "${batch}.RF_intermediate_files.tar.gz"
  }

  runtime {
      docker: "talkowski/sv-pipeline@sha256:7e7e6163d6ac0fc5781eb99ee5a7eec4db37506f48d00f5063b96123f9ca5024"
      memory: "50 GB"
      disks: "local-disk 100 SSD"
      preemptible: 3
  }
}


# Exclude outlier samples from famfile
task update_famfile {
  File famfile
  File excluded_samples
  String batch

  command <<<
    fgrep -wvf ${excluded_samples} ${famfile} > \
      ${batch}.outlier_samples_removed.fam
  >>>

  output {
    File filtered_famfile = "${batch}.outlier_samples_removed.fam"
  }

  runtime {
      docker: "talkowski/sv-pipeline@sha256:7e7e6163d6ac0fc5781eb99ee5a7eec4db37506f48d00f5063b96123f9ca5024"
      preemptible: 3
  }
}


# Merge PESR VCFs
task merge_pesr_vcfs {
  File delly_vcf
  File manta_vcf
  File melt_vcf
  String batch

  command <<<
    vcf-concat ${delly_vcf} ${manta_vcf} ${melt_vcf} \
      | vcf-sort -c \
      | bgzip -c > \
      ${batch}.filtered_pesr_merged.vcf.gz
  >>>

  output {
    File merged_pesr_vcf = "${batch}.filtered_pesr_merged.vcf.gz"
  }

  runtime {
      docker: "talkowski/sv-pipeline@sha256:78c53ed494b7a86c91c98aac15b2c122882fb8636617fe233d72fdd08f933e66"
      disks: "local-disk 100 HDD"
      preemptible: 3
  }
}
