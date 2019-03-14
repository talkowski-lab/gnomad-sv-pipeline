import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04_v2_RD_genotyping_train/versions/7/plain-WDL/descriptor" as RD_genotype_train
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04_v2_PE_genotyping_train/versions/6/plain-WDL/descriptor" as PE_genotype_train
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04_v2_SR_genotyping_train/versions/6/plain-WDL/descriptor" as SR_genotype_train

workflow genotype_pesr_part1 {
  File batch_vcf
  String batch
  String coveragefile     # batch coverage file
  File coveragefile_idx
  File medianfile         # batch median file
  File famfile            # batch famfile
  File svc_acct_key
  File rf_cutoffs         # Random forest cutoffs
  File seed_cutoffs
  Array[String] samples   # List of samples in batch
  Int n_RD_genotype_bins  # number of RdTest bins
  Int n_per_RD_split      # number of variants per RdTest split
  Int n_per_PE_split
  String discfile
  File discfile_idx
  File pesr_blacklist
  String splitfile
  Int n_per_SR_split
  File splitfile_idx
  String reference_build  #hg19 or hg38

  call RD_genotype_train.RD_genotype_train {
    input:
      vcf=batch_vcf,
      coveragefile=coveragefile,
      coveragefile_idx=coveragefile_idx,
      medianfile=medianfile,
      famfile=famfile,
      svc_acct_key=svc_acct_key,
      rf_cutoffs=rf_cutoffs,
      seed_cutoffs=seed_cutoffs,
      samples=samples,
      prefix=batch,
      n_bins=n_RD_genotype_bins,
      n_per_split=n_per_RD_split,
      reference_build=reference_build
  }

  call PE_genotype_train.PE_genotype_train {
    input:
      batch_vcf=batch_vcf,
      discfile=discfile,
      n_per_split=n_per_PE_split,
      medianfile=medianfile,
      discfile_idx=discfile_idx,
      svc_acct_key=svc_acct_key,
      samples=samples,
      batch_ID=batch,
      RF_cutoffs=rf_cutoffs,
      RD_genotypes=RD_genotype_train.genotypes,
      RD_melted_genotypes=RD_genotype_train.melted_genotypes,
      blacklist=pesr_blacklist
  }

  call SR_genotype_train.SR_genotype_train {
    input:
      batch_vcf=batch_vcf,
      splitfile=splitfile,
      n_per_split=n_per_SR_split,
      medianfile=medianfile,
      splitfile_idx=splitfile_idx,
      svc_acct_key=svc_acct_key,
      samples=samples,
      batch_ID=batch,
      RF_cutoffs=rf_cutoffs,
      RD_melted_genotypes=RD_genotype_train.melted_genotypes,
      PE_train=PE_genotype_train.PE_train,
      PE_genotypes=PE_genotype_train.PE_genotypes
  }

  output {
    File RD_pesr_sepcutoff = RD_genotype_train.pesr_sepcutoff
    File RD_depth_sepcutoff = RD_genotype_train.depth_sepcutoff
    File PE_metrics = PE_genotype_train.PE_metrics
    File SR_metrics = SR_genotype_train.SR_metrics
  }
}