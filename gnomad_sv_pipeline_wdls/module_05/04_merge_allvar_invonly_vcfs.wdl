# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# Workflow to parallelize same-bp overlap filter per chromosome
workflow same_bp_filter {
  String vcf
  File vcf_idx
  String prefix
  File contiglist
  File svc_acct_key
  File bothside_pass
  File background_fail

  Array[Array[String]] contigs = read_tsv(contiglist)

  #Run same-bp overlap filter, scattered by chromosome
  scatter (contig in contigs) {

    #Remote tabix each vcf & join into a single vcf
    call subset_vcf {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        contig=contig[0],
        prefix=prefix,
        svc_acct_key=svc_acct_key
    }

    #Run same-bp overlap filter per chromosome
    call bp_overlap_filter {
      input:
        vcf=subset_vcf.subsetted_vcf,
        prefix="${prefix}.${contig[0]}",
        bothside_pass=bothside_pass,
        background_fail=background_fail
    }
  }

  #Merge filtered vcfs across chromosomes
  call concat_vcfs {
    input:
      vcfs=bp_overlap_filter.bp_filtered_vcf,
      prefix="${prefix}.non_redundant"
  }

  output {
    File filtered_vcf = concat_vcfs.concat_vcf
    File filtered_vcf_idx = concat_vcfs.concat_vcf_idx
  }
}


#Remote tabix a single chromosome per VCFs
task subset_vcf {
  String vcf
  File vcf_idx
  String contig
  String prefix
  File svc_acct_key

  command <<<
    #Remote tabix to chromosome of interest
    url=$( gsutil signurl -d 24h ${svc_acct_key} "$vcf" | sed '1d' | cut -f 4 );
    echo $url;
    svtk remote_tabix --header "$url" ${vcf_idx} "${contig}:0-300000000" > "${prefix}.${contig}.vcf"
    bgzip -f "${prefix}.${contig}.vcf"
    tabix -p vcf -f "${prefix}.${contig}.vcf.gz"
  >>>

  output {
    File subsetted_vcf = "${prefix}.${contig}.vcf.gz"
    File subsetted_vcf_idx = "${prefix}.${contig}.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline-remote-pysam@sha256:9fd37fb64e28e54d53172dd30d68c36f0815f21af465381dac281d53755edd86"
    preemptible: 1
    disks: "local-disk 50 SSD"
  }
}


# Run Harrison's overlapping breakpoint filter prior to complex resolution
task bp_overlap_filter {
  File vcf
  String prefix
  File bothside_pass
  File background_fail

  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/overlapbpchange.sh \
    ${vcf} \
    ${background_fail} \
    ${bothside_pass};
    mv non_redundant.vcf.gz "${prefix}.non_redundant.vcf.gz"
  >>>

  output {
    File bp_filtered_vcf = "${prefix}.non_redundant.vcf.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:b359f2cb0c9d5f5a55eb4c41fd362f4e574bf3f8f0f395a2907837571b367ee0"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 250 SSD"
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
    disks: "local-disk 1000 SSD"
  }
}