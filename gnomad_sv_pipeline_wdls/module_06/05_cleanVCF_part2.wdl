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
    docker : "talkowski/sv-pipeline@sha256:a21921b3517e9a10439188b48561f876bcb19bdf56bf21c44264b5bed09b0851"
    disks: "local-disk 250 SSD"
    bootDiskSizeGb: 30
    memory: "16 GB"
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
    docker : "talkowski/sv-pipeline@sha256:facb963613f57bf6c70072c9356241e3ffe47c5d0550beaf9b21f805315846b0"
    disks: "local-disk 200 SSD"
    memory: "4 GB"
  }
  
  output {
		File out="output.txt"
  }
}
