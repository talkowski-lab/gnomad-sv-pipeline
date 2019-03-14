# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# Workflow to parallelize vcf clustering per chromosome
workflow pesr_depth_overlap {
  String pesr_vcf
  File pesr_vcf_idx
  String depth_vcf
  File depth_vcf_idx
  File contigs
  Array[String] samples
  File svc_acct_key

  Array[Array[String]] contiglist = read_tsv(contigs)

  scatter (contig in contiglist) {
    call subset_vcf as subset_pesr_vcf {
      input:
        vcf=pesr_vcf,
        vcf_idx=pesr_vcf_idx,
        contig=contig[0],
        prefix="all_batches.pesr",
        svc_acct_key=svc_acct_key
    }
    
    call subset_vcf as subset_depth_vcf {
      input:
        vcf=depth_vcf,
        vcf_idx=depth_vcf_idx,
        contig=contig[0],
        prefix="all_batches.depth",
        svc_acct_key=svc_acct_key
    }
    
    call merge_pesr_depth {
      input:
        pesr_vcf=subset_pesr_vcf.subsetted_vcf,
        depth_vcf=subset_depth_vcf.subsetted_vcf,
        contig=contig[0]
    }
  }

  call concat_vcfs {
    input:
      vcfs=merge_pesr_depth.merged_vcf,
      prefix="all_batches.pesr_depth"
  }

  output {
    File merged_vcf = concat_vcfs.concat_vcf
    File merged_vcf_idx = concat_vcfs.concat_vcf_idx
  }
}

task subset_vcf {
  String vcf
  File vcf_idx
  String contig
  String prefix
  File svc_acct_key

  command <<<
    # tabix -p vcf ${vcf};
    # tabix -h ${vcf} ${contig} | bgzip -c > ${prefix}.${contig}.vcf.gz
    url=$( gsutil signurl -d 24h ${svc_acct_key} ${vcf} | sed '1d' | cut -f 4 );
    echo $url;
    svtk remote_tabix --header "$url" "${vcf_idx}" "$contig" \
      | bgzip -c \
      > "${prefix}.${contig}.vcf.gz"
  >>>

  output {
    File subsetted_vcf = "${prefix}.${contig}.vcf.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline-remote-pysam@sha256:0c21137179665254ca0d9ebe4d21251ae2ff6679337fd9b3e9d6e6ab808db6a8"
    preemptible: 3
    disks: "local-disk 100 SSD"
  }
}

task merge_pesr_depth {
  File pesr_vcf
  File depth_vcf
  String contig
  
  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/PESR_RD_merge_wrapper.sh \
      ${pesr_vcf} \
      ${depth_vcf} \
      ${contig} \
      all_batches.pesr_depth.${contig}.vcf.gz
  >>>

  output {
    File merged_vcf = "all_batches.pesr_depth.${contig}.vcf.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:3f9d99b8154dff67eb33b0da0a4358ac149461d65f819e7eb64958953d478900"
    preemptible: 1
    memory: "16 GB"
    disks: "local-disk 500 SSD"
  }
}

task concat_vcfs {
  Array[File] vcfs
  String prefix

  command <<<
    vcf-concat ${sep=' ' vcfs} | vcf-sort -c | bgzip -c > ${prefix}.vcf.gz;
    tabix -p vcf -f ${prefix}.vcf.gz
  >>>

  output {
    File concat_vcf = "${prefix}.vcf.gz"
    File concat_vcf_idx = "${prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:b359f2cb0c9d5f5a55eb4c41fd362f4e574bf3f8f0f395a2907837571b367ee0"
    preemptible: 1
    memory: "8 GB"
    disks: "local-disk 5000 SSD"
  }
}