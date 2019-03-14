# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# This is an analysis WDL for downstream processing of Talkowski SV pipeline callests
# that determines categories of SV with high Mendelian violation rates based on
# parent-child trio analyses, and tags those sites as LOW_QUALITY in the
# VCF FILTER field

# QC is performed on the final VCF separated by LOW_QUALITY and non-LOW_QUALITY


import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:master_SV_VCF_QC/versions/73/plain-WDL/descriptor" as QC
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:MVR_collection_helper/versions/5/plain-WDL/descriptor" as collect


workflow assign_lowQuality_sites {
  File vcf
  File vcf_idx
  String prefix
  File contiglist
  File trios_famfile
  File PCRPLUS_samples_list
  Int sv_per_shard

  Array[Array[String]] contigs=read_tsv(contiglist)

  # Shard VCF per-chromosome and collect MVR data
  scatter ( contig in contigs ) {
    call collect.mvr_colection_helper as gather_MVR_data_perChrom {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        contig=contig[0],
        prefix=prefix,
        trios_famfile=trios_famfile,
        PCRPLUS_samples_list=PCRPLUS_samples_list,
        sv_per_shard=sv_per_shard
    }
    call combine_MVR_data as combine_MVR_data_perChrom {
      input:
        MVR_data=gather_MVR_data_perChrom.mvr_data,
        prefix=prefix
    }
  }

  # Merge MVR data
  call combine_MVR_data as combine_MVR_data_crossChrom {
    input:
      MVR_data=combine_MVR_data_perChrom.merged_data,
      prefix=prefix
  }


  # Final outputs
  output {
    File merged_MVR_data = combine_MVR_data_crossChrom.merged_data
  }
}


# Combine MVR data from per-chromosome shards
task combine_MVR_data {
  Array[File] MVR_data
  String prefix

  command <<<
    zcat ${MVR_data[0]} | sed -n '1p' > "${prefix}.merged_MVR_data.txt"
    zcat ${sep=' ' MVR_data} | fgrep -v "#" >> "${prefix}.merged_MVR_data.txt"
    gzip -f "${prefix}.merged_MVR_data.txt"
  >>>

  output {
    File merged_data = "${prefix}.merged_MVR_data.txt.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:07160ad5fad8b8b9faa60a64caf9990e374a47fa63e8f2160d3645f5e4545c48"
    preemptible: 1
    disks: "local-disk 30 HDD"
  }  
}

