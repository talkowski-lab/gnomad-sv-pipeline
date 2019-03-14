# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


#This is an analysis WDL to apply a per-sample GQ cutoff to all variants in an SV VCF 

workflow apply_GQ_filter {
  File vcf
  File vcf_idx
  String prefix
  File contiglist
  Int minGQ_global
  File minGQ_perSVTYPE_table
  Float max_noCallRate

  Array[Array[String]] contigs=read_tsv(contiglist)

  #Split vcf per chromosome
  scatter ( contig in contigs ) {
    #Subset vcf to contig
    call shard_vcf {
        input:
          vcf=vcf,
          vcf_idx=vcf_idx,
          prefix=prefix,
          contig=contig[0]
      }

    #Apply minGQ filter
    call filter_GQ {
      input:
        vcf=shard_vcf.shard,
        prefix="${prefix}.${contig[0]}",
        minGQ_global=minGQ_global,
        minGQ_perSVTYPE_table=minGQ_perSVTYPE_table,
        max_noCallRate=max_noCallRate
    }
  }

  #Merge sharded VCFs
  call combine {
    input:
      vcfs=filter_GQ.filtered_vcf,
      prefix=prefix
  }

  output {
    File filtered_vcf = combine.out
    File filtered_vcf_idx = combine.idx
  }
}


# Shard VCF per chromosome
task shard_vcf {
  File vcf
  File vcf_idx
  String prefix
  String contig

  command {
    tabix -h ${vcf} ${contig} | bgzip -c > "${prefix}.${contig}.vcf.gz"
  }

  output {
    File shard = "${prefix}.${contig}.vcf.gz"
  }
  
  runtime {
    preemptible: 1
    docker: "talkowski/sv-pipeline@sha256:6bcf2b506fc66b13f5aa5e99ccf19e01891aec963b147b09b59e6510116f1adc"
    memory: "4 GB"
    disks: "local-disk 275 SSD"
  }
}


# Apply minGQ filter 
task filter_GQ {
  File vcf
  String prefix
  Int minGQ_global
  File minGQ_perSVTYPE_table
  Float max_noCallRate

  command {
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/apply_minGQ_filter.py \
      --dropEmpties \
      -m ${minGQ_global} \
      --maxNCR ${max_noCallRate} \
      -t ${minGQ_perSVTYPE_table} \
      ${vcf} \
      "${prefix}.filtered.vcf"
    bgzip -f "${prefix}.filtered.vcf"
  }

  output {
    File filtered_vcf = "${prefix}.filtered.vcf.gz"
  }
  
  runtime {
    preemptible: 1
    docker: "talkowski/sv-pipeline@sha256:6bcf2b506fc66b13f5aa5e99ccf19e01891aec963b147b09b59e6510116f1adc"
    memory: "4 GB"
    disks: "local-disk 30 SSD"
  }
}


# Merge VCF shards
task combine {
  Array[File] vcfs
  String prefix
  
  command {
    vcf-concat ${sep=" " vcfs} | vcf-sort | bgzip -c > "${prefix}.minGQ_filtered.vcf.gz";
    tabix -p vcf "${prefix}.minGQ_filtered.vcf.gz"
  }

  runtime {
    preemptible: 1
    docker : "talkowski/sv-pipeline@sha256:6bcf2b506fc66b13f5aa5e99ccf19e01891aec963b147b09b59e6510116f1adc"
    disks: "local-disk 500 SSD"
    memory: "4 GB"
  }

  output {
    File out="${prefix}.minGQ_filtered.vcf.gz"
    File idx="${prefix}.minGQ_filtered.vcf.gz.tbi"
  }
}
