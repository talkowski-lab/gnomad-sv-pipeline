# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License

# Workflow to shard a filtered vcf & run vcfcluster (sub-sub-sub workflow)
workflow sharded_cluster {
  File vcf
  Int dist
  Float frac
  Int max_shards
  Int min_per_shard
  String prefix
  String contig
  String svtype
  Float sample_overlap
  String do_blacklist
  File blacklist
  File blacklist_idx
  Int svsize
  Array[String] svtypes

  #New as of November 2, 2018: perform sharding and return list of variant IDs
  # for each shard, rather than VCF shards themselves, which should dramatically
  # improve speed of sharding task (previously took 1-6 hours for 14k samples in
  # gnomAD v2)
  call shard_vcf {
    input:
      vcf=vcf,
      dist=dist,
      frac=frac,
      max_shards=max_shards,
      min_per_shard=min_per_shard,
      prefix="${prefix}.${contig}.${svtype}"
  }

  #Run vcfcluster per shard
  scatter ( VIDs_list in shard_vcf.VID_list_shards ) {
    call vcfcluster {
      input:
        vcf=vcf,
        VIDs=VIDs_list,
        prefix="${prefix}.${contig}.${svtype}",
        dist=dist,
        frac=frac,
        sample_overlap=sample_overlap,
        do_blacklist=do_blacklist,
        blacklist=blacklist,
        blacklist_idx=blacklist_idx,
        svsize=svsize,
        svtypes=svtypes
    }
  }

  #Merge shards per svtype
  call concat_vcfs as concat_shards {
    input:
      vcfs=vcfcluster.clustered_vcf,
      prefix="${prefix}.${contig}.${svtype}"
  }

  #Output
  output {
    File clustered_vcf = concat_shards.concat_vcf
  }
}


#Intelligently shard a VCF for parallelized clustering
task shard_vcf {
  File vcf
  Int dist
  Float frac
  Int max_shards
  Int min_per_shard
  String prefix

  command <<<
    tabix -f -p vcf ${vcf}
    /opt/sv-pipeline/04_variant_resolution/scripts/shardVCF_preClustering_part1.sh \
    -D ${dist} \
    -R ${frac} \
    -L ${min_per_shard} \
    -S ${max_shards} \
    -P ${prefix} \
    ${vcf}
  >>>

  output {
    Array[File] VID_list_shards = glob("*.VIDs.list")
  }
  
  runtime {
    preemptible: 1
    docker: "talkowski/sv-pipeline@sha256:ca76dffed573c9c792b8362e594bae23b830045d7ce8585bc90202bdcfe4e9a4"
    disks: "local-disk 250 SSD"
  }
}


#Run svtk vcfcluster
task vcfcluster {
  File vcf
  File VIDs
  String prefix
  Int dist
  Float frac
  Float sample_overlap
  String do_blacklist
  File blacklist
  File blacklist_idx
  Int svsize
  Array[String] svtypes
  
  command <<<
    #Generate random hash for tag
    prefix="${prefix}.$( cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 16 | head -n 1 | tr '[:lower:]' '[:upper:]' )"
    #Prep vcf
    zcat ${vcf} | sed -n '1,1000p' | fgrep "#" > header.vcf
    zcat ${vcf} | fgrep -v "#" | fgrep -wf ${VIDs} | cat header.vcf - | bgzip -c \
    > input.vcf.gz
    #Run clustering
    echo "input.vcf.gz" > unclustered_vcfs.list;
    if [ ${do_blacklist} == "YES" ]; then
      svtk vcfcluster unclustered_vcfs.list ${prefix}.vcf \
        -d ${dist} \
        -f ${frac} \
        -x ${blacklist} \
        -z ${svsize} \
        -p ${prefix} \
        -t ${sep=',' svtypes} \
        -o ${sample_overlap} \
        --preserve-ids \
        --preserve-genotypes \
        --preserve-header
      else
        svtk vcfcluster unclustered_vcfs.list ${prefix}.vcf \
        -d ${dist} \
        -f ${frac} \
        -z ${svsize} \
        -p ${prefix} \
        -t ${sep=',' svtypes} \
        -o ${sample_overlap} \
        --preserve-ids \
        --preserve-genotypes \
        --preserve-header
      fi
    bgzip -f ${prefix}.vcf
  >>>

  output {
    File clustered_vcf = "${prefix}.vcf.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:47339fc9bed97946a6c3aa374619b1059ee971ed3728a9d5982b2d7ec82edcf8"
    preemptible: 1
    memory: "8 GB"
    disks: "local-disk 20 SSD"
  }
}


#Merge multiple vcfs
task concat_vcfs {
  Array[File] vcfs
  String prefix

  command <<<
    vcf-concat ${sep=' ' vcfs} | vcf-sort -c | bgzip -c > ${prefix}.vcf.gz
    tabix -f -p vcf ${prefix}.vcf.gz
  >>>

  output {
    File concat_vcf = "${prefix}.vcf.gz"
    File concat_vcf_idx = "${prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:b359f2cb0c9d5f5a55eb4c41fd362f4e574bf3f8f0f395a2907837571b367ee0"
    preemptible: 1
    disks: "local-disk 500 SSD"
  }
}