workflow make_cohort_VCFs {
  File pesr_vcfs_list
  File depth_vcfs_list
  
  call merge_vcfs as merge_pesr_vcfs {
    input:
      vcfs_list=pesr_vcfs_list,
      prefix="all_batches.pesr"
  }
  
  call merge_vcfs as merge_depth_vcfs {
    input:
      vcfs_list=depth_vcfs_list,
      prefix="all_batches.depth"
  }

  output {
  	File cohort_pesr_vcf = merge_pesr_vcfs.merged_vcf
    File cohort_depth_vcf = merge_depth_vcfs.merged_vcf
  }
}

task merge_vcfs {
  File vcfs_list
  String prefix

  command {
    /opt/sv-pipeline/04_variant_resolution/scripts/merge_vcfs.sh ${vcfs_list} ${prefix}
  }

  output {
    File merged_vcf = "${prefix}.vcf.gz"
  }
  
  runtime {
    docker: "talkowski/sv-pipeline@sha256:aaf0b5fa587fbe4f4d137532a4c1be292f9ea104422494e1a7d8ac7a5d8459e6"
    preemptible: 3
  }
}