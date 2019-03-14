import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:05_CleanVCF/versions/93/plain-WDL/descriptor" as CleanVCF_chr
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:master_SV_VCF_QC/versions/75/plain-WDL/descriptor" as QC

workflow CleanVCF_scatter{

  File vcf
  File chrlist
  File backgroundlist
  File famfile
  Int max_shards_per_chrom_step1
  Int min_records_per_shard_step1
  Int samples_per_step2_shard
  File trio_famfile
  String ref_build
  String prefix
  File Sanders_2015_tarball
  File Collins_2017_tarball
  File Werling_2018_tarball
  File? outlier_samples_list


  Array[Array[String]] chrs=read_tsv(chrlist)

    
  scatter ( chr in chrs ){
  	call CleanVCF_chr.CleanVCF {
      input:
      	vcf=vcf,
        Chr=chr[0],
        backgroundlist=backgroundlist,
        famfile=famfile,
        prefix=prefix,
        max_shards_per_chrom_step1=max_shards_per_chrom_step1,
        min_records_per_shard_step1=min_records_per_shard_step1,
        samples_per_step2_shard=samples_per_step2_shard,
        outlier_samples_list=outlier_samples_list
    }
  }

  call combine {
    input:
      vcfs=CleanVCF.out
  }

  call QC.master_vcf_qc as QC_all {
    input:
      vcf=combine.out,
      vcf_idx=combine.idx,
      famfile=trio_famfile,
      ref_build=ref_build,
      prefix="${prefix}_cleanedVCF",
      sv_per_shard=10000,
      samples_per_shard=50,
      Sanders_2015_tarball=Sanders_2015_tarball,
      Collins_2017_tarball=Collins_2017_tarball,
      Werling_2018_tarball=Werling_2018_tarball,
      contiglist=chrlist
  }

  # call subset_pass {
  #   input:
  #     vcf=combine.out,
  #     prefix=prefix
  # }

  # call QC.master_vcf_qc as QC_pass {
  #   input:
  #     vcf=subset_pass.filtered_vcf,
  #     famfile=trio_famfile,
  #     ref_build=ref_build,
  #     prefix="${prefix}_cleanedVCF_filterPass",
  #     sv_per_shard=10000,
  #     samples_per_shard=100,
  #     Sanders_2015_tarball=Sanders_2015_tarball,
  #     Collins_2017_tarball=Collins_2017_tarball,
  #     Werling_2018_tarball=Werling_2018_tarball
  # }
    
  # call subset_fail {
  #   input:
  #     vcf=combine.out,
  #     prefix=prefix
  # }

  # call QC.master_vcf_qc as QC_fail {
  #   input:
  #     vcf=subset_fail.filtered_vcf,
  #     famfile=trio_famfile,
  #     ref_build=ref_build,
  #     prefix="${prefix}_cleanedVCF_filterFail",
  #     sv_per_shard=10000,
  #     samples_per_shard=100,
  #     Sanders_2015_tarball=Sanders_2015_tarball,
  #     Collins_2017_tarball=Collins_2017_tarball,
  #     Werling_2018_tarball=Werling_2018_tarball
  # }
    
	output {
  	File cleaned_vcf = combine.out
    File cleaned_vcf_idx = combine.idx
    File all_variants_QC = QC_all.sv_vcf_qc_output
    # File passing_variants_QC = QC_pass.sv_vcf_qc_output
    # File failing_variants_QC = QC_fail.sv_vcf_qc_output
  }
}


# Merge per-chromosome VCF shards
task combine {

	Array[File] vcfs
  String prefix
  
  command {
    vcf-concat ${sep=" " vcfs} | vcf-sort | bgzip -c > ${prefix}.cleanedvcf.vcf.gz;
    tabix -p vcf ${prefix}.cleanedvcf.vcf.gz
  }

  runtime {
    preemptible: 1
    docker : "talkowski/sv-pipeline@sha256:facb963613f57bf6c70072c9356241e3ffe47c5d0550beaf9b21f805315846b0"
    disks: "local-disk 500 SSD"
    memory: "8 GB"
  }

  output {
    File out="${prefix}.cleanedvcf.vcf.gz"
    File idx="${prefix}.cleanedvcf.vcf.gz.tbi"
  }
}


# Task to sunset variants with VCF FILTER = PASS | MULTIALLELIC
task subset_pass {
  File vcf
  String prefix

  command <<<
    zcat ${vcf} \
      | awk -v FS="\t" -v OFS="\t" \
        '{ if ($1~"#" || $7=="PASS" || $7=="MULTIALLELIC") print $0 }' \
      | vcf-sort \
      | bgzip -c \
      > ${prefix}.passing_variants.vcf.gz
  >>>

  runtime {
    docker: "talkowski/sv-pipeline@sha256:facb963613f57bf6c70072c9356241e3ffe47c5d0550beaf9b21f805315846b0"
    preemptible: 1
    disks: "local-disk 500 SSD"
  }

  output {
    File filtered_vcf = "${prefix}.passing_variants.vcf.gz"
  }
}


# Task to sunset variants with VCF FILTER != PASS | MULTIALLELIC
task subset_fail {
  File vcf
  String prefix

  command <<<
    zcat ${vcf} \
      | awk -v FS="\t" -v OFS="\t" \
        '{ if ($1~"#" || ($7!="PASS" && $7!="MULTIALLELIC") ) print $0 }' \
      | vcf-sort \
      | bgzip -c \
      > ${prefix}.failing_variants.vcf.gz
  >>>

  runtime {
    docker: "talkowski/sv-pipeline@sha256:facb963613f57bf6c70072c9356241e3ffe47c5d0550beaf9b21f805315846b0"
    preemptible: 1
    disks: "local-disk 500 SSD"
  }

  output {
    File filtered_vcf = "${prefix}.failing_variants.vcf.gz"
  }
}
