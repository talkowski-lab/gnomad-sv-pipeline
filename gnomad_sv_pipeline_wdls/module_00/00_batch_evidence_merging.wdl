workflow evidence_merging {
  Array[File] PE_files
  Array[File] PE_indexes
  Array[File] SR_files
  Array[File] SR_indexes
  Array[File] bincov_files
  Array[File] bincov_indexes
  Array[File] BAF_files
  Array[File] BAF_indexes
  Array[String] samples
  String batch
  File inclusion_bed

  call merge_PESR_files as merge_PE_files {
    input:
      files=PE_files,
      indexes=PE_indexes,
      batch=batch,
      evidence="PE",
      inclusion_bed=inclusion_bed
  }
  
  call merge_PESR_files as merge_SR_files {
    input:
      files=SR_files,
      indexes=SR_indexes,
      batch=batch,
      evidence="SR",
      inclusion_bed=inclusion_bed
  }

  call make_bincov_matrix {
    input:
      samples=samples,
      filepaths=bincov_files,
      batch=batch,
  }

  call merge_PESR_files as merge_BAF_files {
    input:
      files=BAF_files,
      indexes=BAF_indexes,
      batch=batch,
      evidence="BAF",
      inclusion_bed=inclusion_bed
  }

  output {
    File merged_PE = merge_PE_files.merged
    File merged_PE_idx = merge_PE_files.merged_idx
    File merged_SR = merge_SR_files.merged
    File merged_SR_idx = merge_SR_files.merged_idx
    File merged_bincov = make_bincov_matrix.bincov_matrix
    File merged_bincov_idx = make_bincov_matrix.bincov_matrix_idx
    File merged_BAF = merge_BAF_files.merged
    File merged_BAF_idx = merge_BAF_files.merged_idx
  }
}

task merge_PESR_files {
  Array[File] files
  Array[File] indexes
  String batch
  String evidence
  File inclusion_bed

  command <<<
    tmpdir=$(mktemp -d);
    cmd="sort -m -k1,1V -k2,2n -T $tmpdir";
    while read file; do
      cmd="$cmd <( tabix -h -R ${inclusion_bed} $file )"
    done < ${write_tsv(files)};
    echo "$cmd"
    eval "$cmd" | bgzip -c > ${batch}.${evidence}.txt.gz;
    tabix -f -s1 -b 2 -e 2 ${batch}.${evidence}.txt.gz
  >>>

  output {
    File merged = "${batch}.${evidence}.txt.gz"
    File merged_idx = "${batch}.${evidence}.txt.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline-remote-pysam"
    memory: "8 GB"
    disks: "local-disk 5000 HDD"
  }
}

task make_bincov_matrix {
  Array[String] samples
  Array[File] filepaths
  String batch

  command <<<
    paste ${write_tsv(samples)} ${write_tsv(filepaths)} > samples.key;
    makeMatrix.sh -z -N -o ${batch}.bincov.bed.gz samples.key
  >>>

  output {
    File bincov_matrix = "${batch}.bincov.bed.gz"
    File bincov_matrix_idx = "${batch}.bincov.bed.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline-remote-pysam"
    disks: "local-disk 1000 HDD"
  }
}
