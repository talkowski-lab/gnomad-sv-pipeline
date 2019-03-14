import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:06_annotate_per_chrom/versions/14/plain-WDL/descriptor" as annotate_by_chrom

# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# Workflow to parallelize VCF annotation by chromosome
workflow parallelized_annotation {
  String vcf
  File vcf_idx
  String prefix
  File contiglist
  File protein_coding_gtf
  # File antisense_gtf
  File lincRNA_gtf
  # File processed_transcript_gtf
  # File pseudogene_gtf
  File promoter_bed
  File noncoding_bed
  File svc_acct_key

  Array[Array[String]] contigs = read_tsv(contiglist)

  #Annotate, scattered by chromosome
  scatter (contig in contigs) {

    #Remote tabix each chromosome
    call subset_vcf {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        contig=contig[0],
        prefix="${prefix}.${contig[0]}",
        svc_acct_key=svc_acct_key
    }

    #Annotate per chromosome
    call annotate_by_chrom.annotate as annotate {
      input:
        vcf=subset_vcf.subsetted_vcf,
        prefix="${prefix}.${contig[0]}",
        protein_coding_gtf=protein_coding_gtf,
        lincRNA_gtf=lincRNA_gtf,
        promoter_bed=promoter_bed,
        noncoding_bed=noncoding_bed
    }
  }

  #Merge integrated vcfs across chromosomes
  call concat_vcfs {
    input:
      vcfs=annotate.annotated_vcf,
      prefix="${prefix}.annotated"
  }

  output {
    File annotated_vcf = concat_vcfs.concat_vcf
    File annotated_vcf_idx = concat_vcfs.concat_vcf_idx
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
    url=$( gsutil signurl -d 24h ${svc_acct_key} ${vcf} | sed '1d' | cut -f 4 );
    echo "$url";
    svtk remote_tabix --header "$url" ${vcf_idx} "${contig}:0-300000000" \
      | bgzip -c > "${prefix}.${contig}.vcf.gz"
    tabix -p vcf -f "${prefix}.${contig}.vcf.gz"
  >>>

  output {
    File subsetted_vcf = "${prefix}.${contig}.vcf.gz"
    File subsetted_vcf_idx = "${prefix}.${contig}.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline-remote-pysam@sha256:13da9601b97e08ce2abb1aca494551dc7c09920e46dcca11768cd6aff3db37e5"
    preemptible: 1
    disks: "local-disk 50 SSD"
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
    docker: "talkowski/sv-pipeline@sha256:17553d54115b8a5f51dc2691dd891e7b9991bb5b7365105ae478cd4a92938a30"
    preemptible: 1
    disks: "local-disk 1000 SSD"
  }
}