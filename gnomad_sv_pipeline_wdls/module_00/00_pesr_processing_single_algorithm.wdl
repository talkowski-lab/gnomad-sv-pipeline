workflow preprocess_algorithm {
  File vcf
  File contigs
  String sample
  String algorithm
  Int min_svsize

  call standardize_vcf {
    input: 
      raw_vcf=vcf,
      algorithm=algorithm,
      group=sample,
      contigs=contigs,
      min_svsize=min_svsize
  }

  call sort_vcf {
    input: 
      unsorted_vcf=standardize_vcf.std_vcf,
      algorithm=algorithm,
      group=sample
  }

  output {
    File std_vcf = sort_vcf.sorted_vcf
  }
}

task standardize_vcf {
  File raw_vcf
  File contigs
  Int min_svsize
  String algorithm
  String group

  command {
    svtk standardize --prefix ${algorithm}_${group} --contigs ${contigs} --min-size ${min_svsize} ${raw_vcf} ${algorithm}.${group}.vcf ${algorithm}
  }

  output { 
    File std_vcf="${algorithm}.${group}.vcf"
    String group_="${group}"
  }
  
  runtime {
    docker: "talkowski/sv-pipeline"
  }
}

task sort_vcf {
  File unsorted_vcf
  String algorithm
  String group
 
  command {
    vcf-sort -c ${unsorted_vcf} | bgzip -c > ${algorithm}.${group}.vcf.gz;
    tabix -p vcf ${algorithm}.${group}.vcf.gz
  }
  
  output {
    File sorted_vcf="${algorithm}.${group}.vcf.gz"
  }
  
  runtime {
    docker: "talkowski/sv-pipeline"
  }
}