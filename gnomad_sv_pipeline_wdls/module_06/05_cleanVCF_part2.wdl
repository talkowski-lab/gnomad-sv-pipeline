workflow Clean {

	Array[File] whitelists
  File normal_revise_vcf
  File multi_cnvs
  File vcftools_idx

	scatter ( white in whitelists ){
		call cleanvcf2{
      input:
        normal_revise_vcf=normal_revise_vcf,
        whitelist=white,
        multi_cnvs=multi_cnvs,
        vcftools_idx=vcftools_idx
  	}
	}

  call combine{
    input:
      shards=cleanvcf2.out
  }
 	
  output {
    File out=combine.out
  }
}


task cleanvcf2 {
  
  File normal_revise_vcf
  File whitelist
  File multi_cnvs
  File vcftools_idx

  command {
      bash /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part2.sh ${normal_revise_vcf} ${whitelist} ${multi_cnvs} "output.txt"
  }

  runtime {
    preemptible: 1
    docker: "talkowski/sv-pipeline@sha256:703a19f84f498989ba8ffde110a3462cfecfbd7ade1084a151fac5fff742c266"
    disks: "local-disk 250 SSD"
    bootDiskSizeGb: 30
    memory: "32 GB"
  }

  output {
    File out="output.txt"
  }
}


task combine {
	
  Array[File] shards

  command {
    cat ${sep=" " shards} > output.txt
  }
  
  runtime {
    preemptible: 1
    docker: "talkowski/sv-pipeline@sha256:703a19f84f498989ba8ffde110a3462cfecfbd7ade1084a151fac5fff742c266"
    disks: "local-disk 200 SSD"
    memory: "4 GB"
  }
  
  output {
		File out="output.txt"
  }
}
