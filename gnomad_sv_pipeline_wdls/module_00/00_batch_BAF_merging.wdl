workflow evidence_merging {
  Array[File] BAF_files
  Array[File] BAF_indexes
  Array[String] samples
  String batch
  File inclusion_bed

  call merge_PESR_files as merge_BAF_files {
    input:
      files=BAF_files,
      indexes=BAF_indexes,
      batch=batch,
      evidence="BAF",
      inclusion_bed=inclusion_bed
  }

  output {
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
