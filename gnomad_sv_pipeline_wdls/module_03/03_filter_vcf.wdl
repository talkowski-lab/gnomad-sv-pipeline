workflow RF_filter_vcf {
  File vcf
  File metrics
  File scores
  File cutoffs
  String prefix

  call filter_vcf {
    input:
      vcf=vcf,
      scores=scores,
      prefix=prefix
  }

  call rewrite_SR_coords {
    input:
      vcf=filter_vcf.filtered_vcf,
      metrics=metrics,
      cutoffs=cutoffs,
      prefix=prefix
  }
  
  call annotate_RF_evidence {
    input:
      vcf=rewrite_SR_coords.corrected_vcf,
      scores=scores,
      prefix=prefix
  }

  output {
    File filtered_vcf = annotate_RF_evidence.annotated_vcf
  }
}

task filter_vcf {
  File vcf
  File scores
  String prefix

  command <<<
    cat \
        <(sed -e '1d' ${scores} | fgrep -e DEL -e DUP | awk '($3>=0.5)' | cut -f1 | fgrep -w -f - <(zcat ${vcf})) \
        <(sed -e '1d' ${scores} | fgrep -e INV -e BND -e INS | awk '($3>=0.5)' | cut -f1 | fgrep -w -f - <(zcat ${vcf}) | sed -e 's/SVTYPE=DEL/SVTYPE=BND/' -e 's/SVTYPE=DUP/SVTYPE=BND/' -e 's/<DEL>/<BND>/' -e 's/<DUP>/<BND>/') \
      | cat <(sed -n -e '/^#/p' <(zcat ${vcf})) - \
      | vcf-sort -c \
      | bgzip -c \
      > ${prefix}.filtered.vcf.gz
  >>>

  output {
    File filtered_vcf = "${prefix}.filtered.vcf.gz"
  }

  runtime {
      docker: "talkowski/sv-pipeline@sha256:7e7e6163d6ac0fc5781eb99ee5a7eec4db37506f48d00f5063b96123f9ca5024"
      preemptible: 3
  }
}

task rewrite_SR_coords {
  File vcf
  File metrics
  File cutoffs
  String prefix

  command <<<
    set -o pipefail;
    /opt/sv-pipeline/03_variant_filtering/scripts/rewrite_SR_coords.py ${vcf} ${metrics} ${cutoffs} stdout \
      | vcf-sort -c \
      | bgzip -c \
      > ${prefix}.corrected_coords.vcf.gz
  >>>

  output {
    File corrected_vcf = "${prefix}.corrected_coords.vcf.gz"
  }

  runtime {
      docker: "talkowski/sv-pipeline@sha256:7e7e6163d6ac0fc5781eb99ee5a7eec4db37506f48d00f5063b96123f9ca5024"
      memory: "10 GB"
      preemptible: 3
  }
}

task annotate_RF_evidence {
  File vcf
  File scores
  String prefix
  
  command <<<
    /opt/sv-pipeline/03_variant_filtering/scripts/annotate_RF_evidence.py ${vcf} ${scores} ${prefix}.with_evidence.vcf;
    bgzip ${prefix}.with_evidence.vcf
  >>>
  
  output {
    File annotated_vcf = "${prefix}.with_evidence.vcf.gz"
  }
  
  runtime {
    docker: "talkowski/sv-pipeline@sha256:7e7e6163d6ac0fc5781eb99ee5a7eec4db37506f48d00f5063b96123f9ca5024"
    preemptible: 3
  }
}