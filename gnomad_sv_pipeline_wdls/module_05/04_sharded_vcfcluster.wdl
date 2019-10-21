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
    set -eu -o pipefail

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
    docker: "talkowski/sv-pipeline@sha256:703a19f84f498989ba8ffde110a3462cfecfbd7ade1084a151fac5fff742c266"
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
    set -eu -o pipefail

    # Don't generate random characters for vcf name, it produces problems with caching on cromwell
    # You *could* pass a seed like so:
    # INPUT_HASH=$(tr -d '/+' < <(openssl enc -a -aes-256-ctr -pass pass:"$SEED" -nosalt </dev/zero 2>/dev/null) | head -c16)
    # But if you hash filtered input vcf, you accomplish the same goal of avoiding similar-named files in the loop,
    # without introducing randomness:
    INPUT_HASH=$(md5sum ${vcf} | awk '{print $1}')
    # concat prefix and hash to create unique vcf name:
    VCF_NAME="${prefix}-$INPUT_HASH"

    #Prep vcf
    zcat ${vcf} | sed -n '1,1000p' | fgrep "#" > header.vcf
    zcat ${vcf} | fgrep -v "#" | fgrep -wf ${VIDs} | cat header.vcf - | bgzip -c \
    > input.vcf.gz
    #Run clustering
    echo "input.vcf.gz" > unclustered_vcfs.list;
    if [ ${do_blacklist} == "YES" ]; then
      svtk vcfcluster unclustered_vcfs.list $VCF_NAME.vcf \
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
        svtk vcfcluster unclustered_vcfs.list $VCF_NAME.vcf \
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
    bgzip -f $VCF_NAME.vcf
  >>>

  output {
    # need to use glob since cromwell will not be aware of the value of INPUT hash
    File clustered_vcf = glob("${prefix}*.vcf.gz")[0]
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:703a19f84f498989ba8ffde110a3462cfecfbd7ade1084a151fac5fff742c266"
    preemptible: 1
    maxRetries: 1
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
    docker: "talkowski/sv-pipeline@sha256:703a19f84f498989ba8ffde110a3462cfecfbd7ade1084a151fac5fff742c266"
    preemptible: 1
    maxRetries: 1
    disks: "local-disk 500 SSD"
  }
}