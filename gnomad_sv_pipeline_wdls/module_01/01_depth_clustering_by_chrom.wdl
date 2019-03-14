workflow bedcluster_by_chrom {
  String batch
  String svtype
  File bed
  File contigs
  Float frac
  String flags
  
  Array[Array[String]] contiglist = read_tsv(contigs)

  scatter (contig in contiglist) {
    call bedcluster {
      input:
        batch=batch,
        svtype=svtype,
        chrom=contig[0],
        bed=bed,
        frac=frac,
        flags=flags
    }
  }

  call concat_beds {
    input:
      batch=batch,
      svtype=svtype,
      beds=bedcluster.clustered_bed
  }

  output {
    File clustered_bed = concat_beds.merged_bed
  }
}

task bedcluster {
  String batch
  String svtype
  String chrom
  File bed
  
  Float frac
  String flags

  command {
    tabix -p bed ${bed};
    svtk bedcluster ${bed} -r ${chrom} \
      -p ${batch}_depth_${svtype}_${chrom} \
      -f ${frac} \
      ${flags} \
      > ${batch}.${svtype}.${chrom}.bed
  }
  
  output {
    File clustered_bed="${batch}.${svtype}.${chrom}.bed"
  } 
  
  runtime {
    docker: "talkowski/sv-pipeline@sha256:a89824ac34b915f605d09bcf57516bc76d950bd762ad5c1f336d421be917be55"
    preemptible: 3
  }
}

task concat_beds {
  String batch
  String svtype
  Array[File] beds

  command <<<
    awk 'FNR==1 && NR!=1 { while (/^#chrom/) getline; } 1 {print}' ${sep=' ' beds} > ${batch}.${svtype}.bed
  >>>
  
  output {
    File merged_bed = "${batch}.${svtype}.bed"
  }
  
  runtime {
    docker: "talkowski/sv-pipeline@sha256:a89824ac34b915f605d09bcf57516bc76d950bd762ad5c1f336d421be917be55"
    preemptible: 3
  }
}