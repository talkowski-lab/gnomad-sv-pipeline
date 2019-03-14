# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# Helper WDL to parallelize collection of Mendelian violation rate data for 
# the Talkowski lab SV pipeline

workflow mvr_colection_helper {
  File vcf
  File vcf_idx
  String contig
  String prefix
  File trios_famfile
  File PCRPLUS_samples_list
  Int sv_per_shard

  call shard_vcf {
    input:
      vcf=vcf,
      vcf_idx=vcf_idx,
      contig=contig,
      sv_per_shard=sv_per_shard
  }

  scatter ( shard in shard_vcf.shard_vcfs ){
    call gather_MVR_data {
    input:
      vcf=shard,
      prefix="${prefix}.${contig}",
      famfile=trios_famfile,
      PCRPLUS_samples_list=PCRPLUS_samples_list
    }
  }

  output {
    Array[File] mvr_data = gather_MVR_data.MVR_data
  }
}


# Shard VCF into fixed size chunks
task shard_vcf {
  File vcf
  File vcf_idx
  String contig
  Int sv_per_shard

  command {
    #Tabix chromosome of interest
    tabix -h ${vcf} ${contig} | bgzip -c > ${contig}.vcf.gz
    #Then shard VCF
    /opt/sv-pipeline/scripts/shard_VCF.sh \
      ${contig}.vcf.gz \
      ${sv_per_shard} \
      "vcf.shard."
  }

  output {
    Array[File] shard_vcfs = glob("vcf.shard.*.vcf.gz")
  }
  
  runtime {
    preemptible: 1
    docker: "talkowski/sv-pipeline@sha256:07160ad5fad8b8b9faa60a64caf9990e374a47fa63e8f2160d3645f5e4545c48"
    memory: "4 GB"
    disks: "local-disk 250 SSD"
  }
}


# Subset compute all data needed for downstream filter determination
task gather_MVR_data {
  File vcf
  String prefix
  File famfile
  File PCRPLUS_samples_list

  command <<<
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/count_mendelian_violations.py \
      ${vcf} ${famfile} ${PCRPLUS_samples_list} "${prefix}.MVR_data.txt"
    gzip -f "${prefix}.MVR_data.txt"
  >>>

  output {
    File MVR_data = "${prefix}.MVR_data.txt.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:58b67cb4e4edf285b89250d2ebab72e17c0247e3bf6891c2c2fcda646b2a6cf4"
    preemptible: 1
    disks: "local-disk 20 HDD"
    memory: "4 GB"
  }
}
