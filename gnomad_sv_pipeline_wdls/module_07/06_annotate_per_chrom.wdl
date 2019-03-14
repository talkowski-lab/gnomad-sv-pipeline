# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# Workflow to annotate output of cleanVCF
workflow annotate {
  File vcf
  String prefix
  File protein_coding_gtf
  # File antisense_gtf
  File lincRNA_gtf
  # File processed_transcript_gtf
  # File pseudogene_gtf
  File promoter_bed
  File noncoding_bed

  call annotate_coding as annotate_protein_coding {
    input:
      vcf=vcf,
      gtf=protein_coding_gtf,
      prefix=prefix,
      gene_set="protein_coding"
  }

  # call annotate_coding as annotate_antisense {
  #   input:
  #     vcf=vcf,
  #     gtf=antisense_gtf,
  #     prefix=prefix,
  #     gene_set="antisense"
  # }

  call annotate_coding as annotate_lincRNA {
    input:
      vcf=vcf,
      gtf=lincRNA_gtf,
      prefix=prefix,
      gene_set="lincRNA"
  }

  # call annotate_coding as annotate_processed_transcript {
  #   input:
  #     vcf=vcf,
  #     gtf=processed_transcript_gtf,
  #     prefix=prefix,
  #     gene_set="processed_transcript"
  # }

  # call annotate_coding as annotate_pseudogene {
  #   input:
  #     vcf=vcf,
  #     gtf=pseudogene_gtf,
  #     prefix=prefix,
  #     gene_set="pseudogene"
  # }

  call annotate_noncoding as annotate_promoter {
    input:
      vcf=vcf,
      bed=promoter_bed,
      prefix=prefix,
      noncoding_set="promoter"
  }

  call annotate_noncoding as annotate_noncoding_elements {
    input:
      vcf=vcf,
      bed=noncoding_bed,
      prefix=prefix,
      noncoding_set="noncoding"
  }
  
  call merge_annotations {
  	input:
      vcf=vcf,
      protein_coding_vcf=annotate_protein_coding.annotated_vcf,
      lincRNA_vcf=annotate_lincRNA.annotated_vcf,
      promoter_vcf=annotate_promoter.annotated_vcf,
      noncoding_vcf=annotate_noncoding_elements.annotated_vcf,
      prefix=prefix
  }
  
  output {
  	File annotated_vcf = merge_annotations.annotated_vcf
  }
}

task annotate_coding {
  File vcf
  File gtf
  String prefix
  String gene_set

  command <<<
    svtk annotate \
      --gencode ${gtf} \
      ${vcf} \
      ${prefix}.${gene_set}.vcf
    orig=$( zcat ${vcf} | cut -f1 | fgrep -v "#" | wc -l )
    new=$( cut -f1 ${prefix}.${gene_set}.vcf | fgrep -v "#" | wc -l )
    if [ "$new" -ne "$orig" ]; then
      echo "ANNOTATED VCF DOES NOT HAVE THE SAME NUMBER OF RECORDS AS INPUT VCF ($new vs $orig)"
      exit 1
    fi
    bgzip -f ${prefix}.${gene_set}.vcf
  >>>

  output {
    File annotated_vcf = "${prefix}.${gene_set}.vcf.gz"
  }

  runtime {
      preemptible: 1
      disks: "local-disk 50 SSD"
      memory: "4 GB"
      docker: "talkowski/sv-pipeline@sha256:9d172112828a0c06dc43685b189fe67551f89d1406d34d351565e9531d5a5860"
  }
}

task annotate_noncoding {
  File vcf
  File bed
  String prefix
  String noncoding_set
  
  command <<<
    svtk annotate \
      --noncoding ${bed} \
      ${vcf} \
      ${prefix}.${noncoding_set}.vcf
    orig=$( zcat ${vcf} | cut -f1 | fgrep -v "#" | wc -l )
    new=$( cut -f1 ${prefix}.${noncoding_set}.vcf | fgrep -v "#" | wc -l )
    if [ "$new" -ne "$orig" ]; then
      echo "ANNOTATED VCF DOES NOT HAVE THE SAME NUMBER OF RECORDS AS INPUT VCF ($new vs $orig)"
      exit 1
    fi
    bgzip -f ${prefix}.${noncoding_set}.vcf
  >>>

  output {
    File annotated_vcf = "${prefix}.${noncoding_set}.vcf.gz"
  }

  runtime {
      preemptible: 1
      disks: "local-disk 50 SSD"
      memory: "4 GB"
      docker: "talkowski/sv-pipeline@sha256:9d172112828a0c06dc43685b189fe67551f89d1406d34d351565e9531d5a5860"
  }
}

task merge_annotations {
  File vcf
  File protein_coding_vcf
  # File antisense_vcf
  File lincRNA_vcf
  # File processed_transcript_vcf
  # File pseudogene_vcf
  File promoter_vcf
  File noncoding_vcf
  String prefix
  
  command <<<
    /opt/sv-pipeline/05_annotation/scripts/merge_annotations.py \
      ${vcf} \
      ${protein_coding_vcf} \
      ${lincRNA_vcf} \
      ${promoter_vcf} \
      ${noncoding_vcf} \
      ${prefix}.annotated.vcf
    bgzip ${prefix}.annotated.vcf
    orig=$( zcat ${vcf} | cut -f1 | fgrep -v "#" | wc -l )
    new=$( zcat ${prefix}.annotated.vcf.gz | cut -f1 | fgrep -v "#" | wc -l )
    if [ "$new" -ne "$orig" ]; then
      echo "ANNOTATED VCF DOES NOT HAVE THE SAME NUMBER OF RECORDS AS INPUT VCF ($new vs $orig)"
      exit 1
    fi
  >>>
  
  output {
    File annotated_vcf = "${prefix}.annotated.vcf.gz"
  }

  runtime {
      preemptible: 1
      disks: "local-disk 250 SSD"
      memory: "8 GB"
      docker: "talkowski/sv-pipeline@sha256:9d172112828a0c06dc43685b189fe67551f89d1406d34d351565e9531d5a5860"
  }
}