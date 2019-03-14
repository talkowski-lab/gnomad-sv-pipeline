workflow preprocess_depth {
  Array[File] beds
  String batch

  call concat_batch as preprocess_DELs {
    input:
      beds=beds,
      batch=batch,
      svtype="DEL"
  }

  call concat_batch as preprocess_DUPs {
    input:
      beds=beds,
      batch=batch,
      svtype="DUP"
  }

  output {
    File del_bed = preprocess_DELs.bed
    File dup_bed = preprocess_DUPs.bed
    File del_bed_idx = preprocess_DELs.bed_idx
    File dup_bed_idx = preprocess_DUPs.bed_idx
  }
}

task concat_batch {
  Array[File] beds
  String svtype
  String batch

  command <<<
    zcat ${sep=' ' beds} \
      | sed -e '/^#chr/d' -e 's/cn.MOPS/cnmops/g' \
      | awk -v svtype=${svtype} '($6==svtype)' \
      | sort -k1,1V -k2,2n \
      | awk -v OFS="\t" -v svtype=${svtype} -v batch=${batch} '{$4=batch"_"svtype"_"NR; print}' \
      | cat <(echo -e "#chr\tstart\tend\tname\tsample\tsvtype\tsources") - \
      | bgzip -c \
      > ${batch}.${svtype}.bed.gz;
  tabix -p bed ${batch}.${svtype}.bed.gz
  >>>

  output {
    File bed="${batch}.${svtype}.bed.gz"
    File bed_idx="${batch}.${svtype}.bed.gz.tbi"
  }
  
  runtime {
    docker: "talkowski/sv-pipeline"
    preemptible: 3
  }
}