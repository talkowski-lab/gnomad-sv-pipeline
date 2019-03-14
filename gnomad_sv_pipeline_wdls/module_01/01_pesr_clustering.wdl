import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:01_pesr_clustering_single_algorithm/versions/15/plain-WDL/descriptor" as single

workflow cluster_pesr {
  Array[File] manta_vcfs
  Array[File] delly_vcfs
  Array[File] melt_vcfs
  File contigs
  String batch
  File trios_famfile
  String ref_build
  File Sanders_2015_tarball
  File Collins_2017_tarball
  File Werling_2018_tarball
  
  Int dist
  Float frac
  File blacklist
  Int svsize
  String flags

  call single.cluster_pesr_algorithm as cluster_manta {
    input:
      vcfs=manta_vcfs,
      batch=batch,
      algorithm="manta",
      contigs=contigs,
      dist=dist,
      frac=frac,
      blacklist=blacklist,
      svsize=svsize,
      flags=flags,
      svtypes="DEL,DUP,INV,BND,INS",
      famfile=famfile,
      ref_build=ref_build,
      Sanders_2015_tarball=Sanders_2015_tarball,
      Werling_2018_tarball=Werling_2018_tarball
  }

  call single.cluster_pesr_algorithm as cluster_delly {
    input:
      vcfs=delly_vcfs,
      batch=batch,
      algorithm="delly",
      contigs=contigs,
      dist=dist,
      frac=frac,
      blacklist=blacklist,
      svsize=svsize,
      flags=flags,
      svtypes="DEL,DUP,INV,BND",
      famfile=famfile,
      ref_build=ref_build,
      Sanders_2015_tarball=Sanders_2015_tarball,
      Werling_2018_tarball=Werling_2018_tarball
  }
  
  call single.cluster_pesr_algorithm as cluster_melt {
    input:
      vcfs=melt_vcfs,
      batch=batch,
      algorithm="melt",
      contigs=contigs,
      dist=dist,
      frac=frac,
      blacklist=blacklist,
      svsize=svsize,
      flags=flags,
      svtypes="INS",
      famfile=famfile,
      ref_build=ref_build,
      Sanders_2015_tarball=Sanders_2015_tarball,
      Werling_2018_tarball=Werling_2018_tarball
  }

  # call merge_vcf_qc {
  #   input:
  #     manta_vcf_qc=cluster_manta.clustered_vcf_qc,
  #     delly_vcf_qc=cluster_delly.clustered_vcf_qc,
  #     melt_vcf_qc=cluster_melt.clustered_vcf_qc
  # }

  output {
    File manta_vcf = cluster_manta.clustered_vcf
    File manta_vcf_qc = cluster_manta.clustered_vcf_qc
    File delly_vcf = cluster_delly.clustered_vcf
    File delly_vcf_qc = cluster_delly.clustered_vcf_qc
    File melt_vcf = cluster_melt.clustered_vcf
    File melt_vcf_qc = cluster_melt.clustered_vcf_qc
    # File merged_vcf_qc = merge_vcf_qc.merged_qc
  }
}

# task merge_vcf_qc {
#   File manta_vcf_qc
#   File delly_vcf_qc
#   File melt_vcf_qc

#   command <<<
#     mkdir merged_pesr_clustering_vcf_qc/
#     mv ${manta_vcf_qc} merged_pesr_clustering_vcf_qc/
#     mv ${delly_vcf_qc} merged_pesr_clustering_vcf_qc/
#     mv ${melt_vcf_qc} merged_pesr_clustering_vcf_qc/
#     tar -czvf merged_pesr_clustering_vcf_qc.tar.gz \
#       merged_pesr_clustering_vcf_qc
#   >>>

#   output {
#     File merged_qc = "merged_pesr_clustering_vcf_qc.tar.gz"
#   }

#     runtime {
#     docker: "talkowski/sv-pipeline"
#     preemptible: 3
#   }
# }