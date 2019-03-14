# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


#This is a helper WDL that runs svtk vcf2bed parallelized across many shards for
# a single chromosome

workflow sharded_vcf2bed {
  File vcf
  File vcf_idx
  String contig
  Int sv_per_shard
  String prefix

  # Tabix to chromosome of interest, and shard input VCF for stats collection
  call shard_vcf {
    input:
      vcf=vcf,
      vcf_idx=vcf_idx,
      contig=contig,
      sv_per_shard=sv_per_shard
  }

  # Scatter over VCF shards
  scatter (shard in shard_vcf.shard_vcfs) {
    # Run vcf2bed
    call vcf2bed_sub {
      input:
        vcf=shard,
        prefix="${prefix}.shard"
      }
    }

  # Merge vcf2bed_sub outputs
  call merge_vcf2bed_sub {
    input:
      vcf2bed_sub_shards=vcf2bed_sub.vcf2bed_sub_out,
      prefix=prefix
  }

  output {
    File vcf2bed_out=merge_vcf2bed_sub.merged_vcf2bed_out
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
    docker: "talkowski/sv-pipeline@sha256:ec7e6f578ba2a8796399fc6f0f9864ec2d34a4921c769a8a54bbbf5254337a8b"
    memory: "4 GB"
    disks: "local-disk 270 SSD"
  }
}


# Run vcf2bed_sub on an input vcf
task vcf2bed_sub {
  File vcf
  String prefix

  command {
    svtk vcf2bed \
      --info ALL \
      --include-filters \
      --no-samples \
      ${vcf} \
      stdout \
      | bgzip -c \
      > "${prefix}.vcf2bed.bed.gz"
  }

  output {
    File vcf2bed_sub_out = "${prefix}.vcf2bed.bed.gz"
  }
  
  runtime {
    preemptible: 1
    docker: "talkowski/sv-pipeline@sha256:ec7e6f578ba2a8796399fc6f0f9864ec2d34a4921c769a8a54bbbf5254337a8b"
    memory: "4 GB"
    disks: "local-disk 25 SSD"
  }
}

# Merge vcf2bed_sub shards
task merge_vcf2bed_sub {
  Array[File] vcf2bed_sub_shards
  String prefix

  command <<<
    zcat ${vcf2bed_sub_shards[0]} | sed -n '1p' > header.txt
    zcat ${sep=' ' vcf2bed_sub_shards} | fgrep -v "#" \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | cat header.txt - \
      | bgzip -c \
      > "${prefix}.vcf2bed_sub.bed.gz"
  >>>

  output {
    File merged_vcf2bed_out = "${prefix}.vcf2bed_sub.bed.gz"
  }
  
  runtime {
    preemptible: 1
    docker: "talkowski/sv-pipeline@sha256:ec7e6f578ba2a8796399fc6f0f9864ec2d34a4921c769a8a54bbbf5254337a8b"
    memory: "4 GB"
    disks: "local-disk 200 SSD"
  }
}


