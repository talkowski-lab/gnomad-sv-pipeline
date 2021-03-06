# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Helper workflow to calculate basic AF statistics for a single chromosome on an input VCF

workflow getAFs_singleChrom {
	File vcf
  File vcf_idx
  String contig
  Int sv_per_shard
	String prefix
  File? sample_pop_assignments  #Two-column file with sample ID & pop assignment. "." for pop will ignore sample
  File? famfile                 #Used for M/F AF calculations
  String? drop_empty_records


  # Tabix to chromosome of interest, and shard input VCF for stats collection
  call shard_vcf {
    input:
      vcf=vcf,
      vcf_idx=vcf_idx,
      contig=contig,
      sv_per_shard=sv_per_shard
  }

  # Scatter over VCF shards
  scatter ( shard in shard_vcf.shard_vcfs ) {
    # Collect AF summary stats
  	call compute_shard_AFs {
  		input:
        vcf=shard,
        prefix="${prefix}.${contig}",
        sample_pop_assignments=sample_pop_assignments,
        famfile=famfile
      }
  	}

  # Merge shards into single VCF
  call combine_sharded_vcfs {
    input:
      vcfs=compute_shard_AFs.shard_wAFs,
      prefix="${prefix}.${contig}",
      drop_empty_records=drop_empty_records
  }

  # Final output
  output {
    File vcf_wAFs = combine_sharded_vcfs.vcf_out
    File vcf_wAFs_idx = combine_sharded_vcfs.vcf_out_idx
  }
}


# Shard VCF into fixed size chunks
task shard_vcf {
  File vcf
  File vcf_idx
  String contig
  Int sv_per_shard

  command {
    set -euo pipefail
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
    maxRetries: 1
    docker: "talkowski/sv-pipeline@sha256:193d18c26100fdd603c569346722513f5796685e990ec3abcaeb4be887062a1a"
    memory: "4 GB"
    disks: "local-disk 250 SSD"
  }
}


# Subset a vcf to a single chromosome, and add global AF information (no subpop)
task compute_shard_AFs {
  File vcf
  String prefix
  File? sample_pop_assignments
  File? famfile


  command <<<
    set -euo pipefail
    optionals=" "
    if [ ${default="SKIP" sample_pop_assignments} != "SKIP" ]; then
      optionals="$( echo "$optionals" ) -p ${sample_pop_assignments}"
    fi
    if [ ${default="SKIP" famfile} != "SKIP" ]; then
      optionals="$( echo "$optionals" ) -f ${famfile}"
    fi
    echo -e "OPTIONALS INTERPRETED AS: $optionals"
    echo -e "NOW RUNNING: /opt/sv-pipeline/05_annotation/scripts/compute_AFs.py $( echo "$optionals" ) ${vcf} stdout"
    #Tabix chromosome of interest & compute AN, AC, and AF
    /opt/sv-pipeline/05_annotation/scripts/compute_AFs.py $optionals "${vcf}" stdout \
    | bgzip -c \
    > "${prefix}.wAFs.vcf.gz"
  >>>

  output {
    File shard_wAFs = "${prefix}.wAFs.vcf.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:193d18c26100fdd603c569346722513f5796685e990ec3abcaeb4be887062a1a"
    preemptible: 1
    maxRetries: 1
    memory: "4 GB"
    disks: "local-disk 20 SSD"
  }
}


# Merge VCF shards & drop records with zero remaining non-ref alleles
task combine_sharded_vcfs {
  Array[File] vcfs
  String prefix
  String? drop_empty_records
  
  command {
    set -euo pipefail
    vcf-concat ${sep=" " vcfs} \
    | vcf-sort \
    > merged.vcf
    if [ ${default="TRUE" drop_empty_records} == "TRUE" ]; then
      /opt/sv-pipeline/05_annotation/scripts/prune_allref_records.py \
        merged.vcf stdout \
      | bgzip -c \
      > "${prefix}.wAFs.vcf.gz"
    else
      cat merged.vcf | bgzip -c > "${prefix}.wAFs.vcf.gz"
    fi
    tabix -p vcf "${prefix}.wAFs.vcf.gz"
  }

  runtime {
    preemptible: 1
    maxRetries: 1
    docker: "talkowski/sv-pipeline@sha256:193d18c26100fdd603c569346722513f5796685e990ec3abcaeb4be887062a1a"
    disks: "local-disk 250 SSD"
    bootDiskSizeGb: 30
    memory: "4 GB"
  }

  output {
    File vcf_out = "${prefix}.wAFs.vcf.gz"
    File vcf_out_idx = "${prefix}.wAFs.vcf.gz.tbi"
  }
}

