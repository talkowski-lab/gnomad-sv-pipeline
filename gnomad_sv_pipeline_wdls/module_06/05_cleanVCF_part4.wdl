workflow Clean4{

  Array[File] RD_CN_revises
  File normal_revise_vcf
  
  scatter ( RD_CN_revise in  RD_CN_revises ){
    call cleanvcf4 {
      input:
        RD_CN_revise=RD_CN_revise,
        normal_revise_vcf=normal_revise_vcf,
    }
  }

  call combine as combine_revised {
    input:
      shards=cleanvcf4.out,
      outfile="revise.vcf.lines.txt.gz"
  }

  call combine as combine_multi_IDs {
    input:
      shards=cleanvcf4.multi_IDs,
      outfile="multi.geno.ids.txt.gz"
  }

  output {
      File out=combine_revised.out
      File multi_IDs=combine_multi_IDs.out
  }
}

task cleanvcf4 {
  File RD_CN_revise
  File normal_revise_vcf

  command {
      bash /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part4.sh ${RD_CN_revise} ${normal_revise_vcf} 
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:703a19f84f498989ba8ffde110a3462cfecfbd7ade1084a151fac5fff742c266"
    disks: "local-disk 200 SSD"
    memory: "16 GB"
  }

  output {
    File out="revise.vcf.lines.txt.gz"
    File multi_IDs="multi.geno.ids.txt.gz"
  }
}

task combine {
  Array[File] shards
  String outfile

  command {
      zcat ${sep=" " shards} | bgzip -c > ${outfile}
  }

  runtime {
    preemptible: 1
    docker: "talkowski/sv-pipeline@sha256:703a19f84f498989ba8ffde110a3462cfecfbd7ade1084a151fac5fff742c266"
    disks: "local-disk 250 SSD"
    memory: "8 GB"
  }

  output {
    File out="${outfile}"
  }
}
