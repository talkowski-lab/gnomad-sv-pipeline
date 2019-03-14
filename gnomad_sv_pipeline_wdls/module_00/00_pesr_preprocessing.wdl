import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:00_pesr_processing_single_algorithm/versions/4/plain-WDL/descriptor" as pp

workflow preprocess_pesr {
  String sample         # Sample ID
  File manta_vcf        # Manta VCF
  File delly_vcf        # Delly VCF
  File melt_vcf         # Melt VCF
  File contigs          # .fai file of whitelisted contigs
  Int min_svsize        # Minimum SV length to include

  call pp.preprocess_algorithm as process_manta {
    input:
      vcf=manta_vcf,
      contigs=contigs,
      min_svsize=min_svsize,
      algorithm="manta",
      sample=sample
  }

  call pp.preprocess_algorithm as process_delly {
    input:
      vcf=delly_vcf,
      contigs=contigs,
      min_svsize=min_svsize,
      algorithm="delly",
      sample=sample
  }
  
  call pp.preprocess_algorithm as process_melt {
    input:
      vcf=melt_vcf,
      contigs=contigs,
      min_svsize=min_svsize,
      algorithm="melt",
      sample=sample
  }

  output {
    File std_manta_vcf = process_manta.std_vcf
    File std_delly_vcf = process_delly.std_vcf
    File std_melt_vcf = process_melt.std_vcf
  }
}