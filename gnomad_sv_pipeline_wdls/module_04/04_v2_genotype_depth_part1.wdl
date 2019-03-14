import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04_v2_RD_genotyping_train/versions/7/plain-WDL/descriptor" as RD_genotype_train

workflow genotype_depth_part1 {
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

  output {
    File RD_pesr_sepcutoff = RD_genotype_train.pesr_sepcutoff
    File RD_depth_sepcutoff = RD_genotype_train.depth_sepcutoff
  }
}