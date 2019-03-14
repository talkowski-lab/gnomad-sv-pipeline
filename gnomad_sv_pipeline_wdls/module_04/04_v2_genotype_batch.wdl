import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04_v2_genotype_pesr_part1/versions/12/plain-WDL/descriptor" as genotype_pesr_part1
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04_v2_genotype_pesr_part2/versions/43/plain-WDL/descriptor" as genotype_pesr_part2
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04_v2_genotype_depth_part1/versions/7/plain-WDL/descriptor" as genotype_depth_part1
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04_v2_genotype_depth_part2/versions/14/plain-WDL/descriptor" as genotype_depth_part2
# import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:master_SV_VCF_QC/versions/49/plain-WDL/descriptor" as vcf_qc

workflow genotype {
  File batch_pesr_vcf
  File batch_depth_vcf
  File cohort_pesr_vcf
  File cohort_depth_vcf
  String batch
  Int n_per_split
  String ref_build
  File trios_famfile
  File Sanders_2015_tarball
  File Werling_2018_tarball
  File Collins_2017_tarball

  String coveragefile     # batch coverage file
  File coveragefile_idx
  File medianfile         # batch median file
  File famfile            # batch famfile
  File svc_acct_key
  File rf_cutoffs         # Random forest cutoffs
  File seed_cutoffs
  File samples_list   # List of samples in batch
  Int n_RD_genotype_bins  # number of RdTest bins
  String discfile
  File discfile_idx
  File pesr_blacklist
  String splitfile
  File splitfile_idx
  String reference_build  #hg19 or hg38

  # Convert samples_list to array
  Array[String] samples = read_lines(samples_list)

  call add_batch_samples as add_batch_samples_pesr {
    input:
      batch_vcf=batch_pesr_vcf,
      cohort_vcf=cohort_pesr_vcf,
      prefix="${batch}.pesr"
  }
  
  call add_batch_samples as add_batch_samples_depth {
    input:
      batch_vcf=batch_depth_vcf,
      cohort_vcf=cohort_depth_vcf,
      prefix="${batch}.depth"
  }

  call genotype_pesr_part1.genotype_pesr_part1 {
    input:
      batch_vcf=batch_pesr_vcf,
      batch=batch,
      coveragefile=coveragefile,
      coveragefile_idx=coveragefile_idx,
      medianfile=medianfile,
      famfile=famfile,
      svc_acct_key=svc_acct_key,
      rf_cutoffs=rf_cutoffs,
      seed_cutoffs=seed_cutoffs,
      samples=samples,
      n_RD_genotype_bins=n_RD_genotype_bins,
      discfile=discfile,
      discfile_idx=discfile_idx,
      splitfile=splitfile,
      splitfile_idx=splitfile_idx,
      pesr_blacklist=pesr_blacklist,
      n_per_RD_split=n_per_split,
      n_per_PE_split=n_per_split,
      n_per_SR_split=n_per_split,
      reference_build=reference_build
  }

  call genotype_pesr_part2.genotype_pesr_part2 {
    input:
      cohort_vcf=add_batch_samples_pesr.updated_vcf,
      RD_pesr_sepcutoff=genotype_pesr_part1.RD_pesr_sepcutoff,
      RD_depth_sepcutoff=genotype_pesr_part1.RD_depth_sepcutoff,
      PE_metrics=genotype_pesr_part1.PE_metrics,
      SR_metrics=genotype_pesr_part1.SR_metrics,
      n_per_split=n_per_split,
      n_RdTest_bins=n_RD_genotype_bins,
      batch=batch,
      medianfile=medianfile,
      famfile=famfile,
      svc_acct_key=svc_acct_key,
      samples=samples,
      coveragefile=coveragefile,
      coveragefile_idx=coveragefile_idx,
      discfile=discfile,
      discfile_idx=discfile_idx,
      splitfile=splitfile,
      splitfile_idx=splitfile_idx
  }

  # call vcf_qc.master_vcf_qc as pesr_qc {
  #   input:
  #     vcf=genotype_pesr_part2.genotyped_vcf,
  #     famfile=trios_famfile,
  #     ref_build=ref_build,
  #     prefix="${batch}_genotyped_pesr_vcf",
  #     sv_per_shard=10000,
  #     samples_per_shard=100,
  #     Sanders_2015_tarball=Sanders_2015_tarball,
  #     Werling_2018_tarball=Werling_2018_tarball,
  #     Collins_2017_tarball=Collins_2017_tarball
  # }

  call genotype_depth_part1.genotype_depth_part1 {
    input:
      batch_vcf=batch_depth_vcf,
      batch=batch,
      coveragefile=coveragefile,
      coveragefile_idx=coveragefile_idx,
      medianfile=medianfile,
      famfile=famfile,
      svc_acct_key=svc_acct_key,
      rf_cutoffs=rf_cutoffs,
      seed_cutoffs=seed_cutoffs,
      samples=samples,
      n_RD_genotype_bins=n_RD_genotype_bins,
      n_per_RD_split=n_per_split,
      reference_build=reference_build
  }

  call genotype_depth_part2.genotype_depth_part2 {
    input:
      cohort_vcf=add_batch_samples_depth.updated_vcf,
      RD_pesr_sepcutoff=genotype_depth_part1.RD_pesr_sepcutoff,
      RD_depth_sepcutoff=genotype_depth_part1.RD_depth_sepcutoff,
      n_per_split=n_per_split,
      n_RdTest_bins=n_RD_genotype_bins,
      batch=batch,
      medianfile=medianfile,
      famfile=famfile,
      svc_acct_key=svc_acct_key,
      samples=samples,
      coveragefile=coveragefile,
      coveragefile_idx=coveragefile_idx
  }

  #Index vcfs
  call index_vcf as index_pesr {
    input:
      vcf=genotype_pesr_part2.genotyped_vcf,
      prefix="${batch}.pesr"
  }
  call index_vcf as index_depth {
    input:
      vcf=genotype_depth_part2.genotyped_vcf,
      prefix="${batch}.depth"
  }

  # call vcf_qc.master_vcf_qc as depth_qc {
  #   input:
  #     vcf=genotype_depth_part2.genotyped_vcf,
  #     famfile=trios_famfile,
  #     ref_build=ref_build,
  #     prefix="${batch}_genotyped_depth_vcf",
  #     sv_per_shard=10000,
  #     samples_per_shard=100,
  #     Sanders_2015_tarball=Sanders_2015_tarball,
  #     Werling_2018_tarball=Werling_2018_tarball,
  #     Collins_2017_tarball=Collins_2017_tarball
  # }

  output {
    File genotyped_pesr_vcf = genotype_pesr_part2.genotyped_vcf
    File genotyped_pesr_vcf_idx = index_pesr.vcf_idx
    # File genotyped_pesr_vcf_qc = pesr_qc.sv_vcf_qc_output
    File genotyped_depth_vcf = genotype_depth_part2.genotyped_vcf
    File genotyped_depth_vcf_idx = index_depth.vcf_idx
    # File genotyped_depth_vcf_qc = depth_qc.sv_vcf_qc_output
    File sr_background_fail = genotype_pesr_part2.background_fail
    File sr_bothside_pass = genotype_pesr_part2.bothside_pass
    File pesr_genotyping_RD_pesr_sepcutoff = genotype_pesr_part1.RD_pesr_sepcutoff
    File pesr_genotyping_RD_depth_sepcutoff = genotype_pesr_part1.RD_depth_sepcutoff
    File depth_genotyping_RD_pesr_sepcutoff = genotype_depth_part1.RD_pesr_sepcutoff
    File depth_genotying_RD_depth_sepcutoff = genotype_depth_part1.RD_depth_sepcutoff
    # File pesr_merged_pe_counts = genotype_pesr_part2.pe_counts
    # File pesr_merged_sr_counts = genotype_pesr_part2.sr_counts
  }
}


task add_batch_samples {
  File batch_vcf
  File cohort_vcf
  String prefix
  
  command <<<
  	/opt/sv-pipeline/04_variant_resolution/scripts/add_batch_samples.py ${batch_vcf} ${cohort_vcf} ${prefix}.vcf;
    bgzip ${prefix}.vcf
  >>>
  
  output {
  	File updated_vcf = "${prefix}.vcf.gz"
  }
  
  runtime {
    docker: 'talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c'
    preemptible: 3
  }
}


task integrate_pesr_GQ {
  File vcf
  File RD_melted_genotypes
  File RD_vargq
  File PE_genotypes
  File PE_vargq
  File SR_genotypes
  File SR_vargq
  
  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/IntegrateGQ.sh \
      ${vcf} \
      ${RD_melted_genotypes} \
      ${RD_vargq} \
      ${PE_genotypes} \
      ${PE_vargq} \
      ${SR_genotypes} \
      ${SR_vargq}
  >>>
  
  output {
    File genotypes = "genotype.indiv.txt.gz"
    File varGQ = "genotype.variant.txt.gz"
  }
    
  runtime {
    docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
    preemptible: 3
    memory: "16 GB"
    disks: "local-disk 200 SSD"
  }
}


task integrate_depth_GQ {
  File vcf
  File RD_melted_genotypes
  File RD_vargq
  
  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/IntegrateGQ_depthonly.sh \
      ${vcf} \
      ${RD_melted_genotypes} \
      ${RD_vargq}
  >>>
  
  output {
    File genotypes = "genotype.indiv.depth.txt.gz"
    File varGQ = "genotype.variant.depth.txt.gz"
  }
    
  runtime {
    docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
    preemptible: 3
    memory: "16 GB"
    disks: "local-disk 50 SSD"
  }
}


task add_genotypes {
  File vcf
  File genotypes
  File varGQ
  String prefix
  
  command <<<
  	echo ""
    /opt/sv-pipeline/04_variant_resolution/scripts/add_genotypes.py \
        ${vcf} \
        ${genotypes} \
        ${varGQ} \
        ${prefix}.genotyped.vcf;
    vcf-sort -c ${prefix}.genotyped.vcf | bgzip -c > ${prefix}.genotyped.vcf.gz
  >>>
  
  output {
      File genotyped_vcf = "${prefix}.genotyped.vcf.gz"
  }
  
  runtime {
    docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
    preemptible: 3
    memory: "16 GB"
    disks: "local-disk 50 SSD"
  }
}


task index_vcf {
  File vcf
  String prefix

  command <<<
    tabix -f -p vcf ${vcf}
    mv $( find `pwd` -name "*.tbi" ) ${prefix}.genotyped.vcf.gz.tbi
  >>>

  output {
    File vcf_idx = "${prefix}.genotyped.vcf.gz.tbi"
  }
  
  runtime {
    docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
    preemptible: 3
    disks: "local-disk 50 SSD"
  }
}