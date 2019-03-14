# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


#This is an analysis WDL to help parameterize an optimal per-sample GQ cutoff 
# based on minimizing de novo rate while retaining as many resolved variants 
# as possible based on a ROC curve

import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:minGQ_ROC_helper/versions/5/plain-WDL/descriptor" as ROC_helper
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:sharded_vcf2bed/versions/5/plain-WDL/descriptor" as parallel_v2b

workflow optimize_GQ_filter {
  File vcf
  File vcf_idx
  File famfile
  String prefix
  Int fams_per_shard
  Int size_boundary
  Float max_fdr_small_lowSR
  Float max_fdr_small_highSR
  Float max_fdr_large_lowSR
  Float max_fdr_large_highSR
  Int minGQ
  Int maxGQ
  Int GQstepsize
  File contiglist

  Array[Array[String]] contigs = read_tsv(contiglist)


  #Split famfile into even chunks
  call split_famfile {
    input:
      famfile=famfile,
      vcf=vcf,
      fams_per_shard=fams_per_shard
  }


  #Run parallelized svtk vcf2bed
  scatter ( contig in contigs ) {
    call parallel_v2b.sharded_vcf2bed as parallel_vcf2bed {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        contig=contig[0],
        sv_per_shard=10000,
        prefix="${prefix}.${contig[0]}"
    }
  }
  call merge_vcf2bed {
    input:
      vcf2bed_shards=parallel_vcf2bed.vcf2bed_out,
      prefix=prefix
  }


  #Optimization: small & lowSR
  call filter_vcf as filter_small_lowSR {
    input:
      vcf=vcf,
      bed=merge_vcf2bed.merged_vcf2bed_out,
      prefix="${prefix}.small_lowSR",
      min_size=-2,
      max_size=size_boundary,
      highSR="FALSE"
  }
  call ROC_helper.optimize_ROC as ROC_small_lowSR {
    input:
      vcf=filter_small_lowSR.filtered_vcf,
      famfiles=split_famfile.famfile_shards,
      prefix="${prefix}.small_lowSR",
      max_fdr=max_fdr_small_lowSR,
      minGQ=minGQ,
      maxGQ=maxGQ,
      GQstepsize=GQstepsize
  }


  #Optimization: small & highSR
  call filter_vcf as filter_small_highSR {
    input:
      vcf=vcf,
      bed=merge_vcf2bed.merged_vcf2bed_out,
      prefix="${prefix}.small_highSR",
      min_size=-2,
      max_size=size_boundary,
      highSR="TRUE"
  }
  call ROC_helper.optimize_ROC as ROC_small_highSR {
    input:
      vcf=filter_small_highSR.filtered_vcf,
      famfiles=split_famfile.famfile_shards,
      prefix="${prefix}.small_highSR",
      max_fdr=max_fdr_small_highSR,
      minGQ=minGQ,
      maxGQ=maxGQ,
      GQstepsize=GQstepsize
  }


  #Optimization: large & lowSR
  call filter_vcf as filter_large_lowSR {
    input:
      vcf=vcf,
      bed=merge_vcf2bed.merged_vcf2bed_out,
      prefix="${prefix}.large_lowSR",
      min_size=size_boundary,
      max_size=400000000,
      highSR="FALSE"
  }
  call ROC_helper.optimize_ROC as ROC_large_lowSR {
    input:
      vcf=filter_large_lowSR.filtered_vcf,
      famfiles=split_famfile.famfile_shards,
      prefix="${prefix}.large_lowSR",
      max_fdr=max_fdr_large_lowSR,
      minGQ=minGQ,
      maxGQ=maxGQ,
      GQstepsize=GQstepsize
  }


  #Optimization: large & highSR
  call filter_vcf as filter_large_highSR {
    input:
      vcf=vcf,
      bed=merge_vcf2bed.merged_vcf2bed_out,
      prefix="${prefix}.large_highSR",
      min_size=size_boundary,
      max_size=400000000,
      highSR="TRUE"
  }
  call ROC_helper.optimize_ROC as ROC_large_highSR {
    input:
      vcf=filter_large_highSR.filtered_vcf,
      famfiles=split_famfile.famfile_shards,
      prefix="${prefix}.large_highSR",
      max_fdr=max_fdr_large_highSR,
      minGQ=minGQ,
      maxGQ=maxGQ,
      GQstepsize=GQstepsize
  }


  output {
    File small_lowSR_minGQ_ROC_plot = ROC_small_lowSR.minGQ_ROC_plot
    File small_lowSR_minGQ_ROC_table = ROC_small_lowSR.minGQ_ROC_table
    File small_lowSR_vcf_subset = filter_small_lowSR.filtered_vcf
    File small_lowSR_vcf_subset_idx = filter_small_lowSR.filtered_vcf_idx
    File small_highSR_minGQ_ROC_plot = ROC_small_highSR.minGQ_ROC_plot
    File small_highSR_minGQ_ROC_table = ROC_small_highSR.minGQ_ROC_table
    File small_highSR_vcf_subset = filter_small_highSR.filtered_vcf
    File small_highSR_vcf_subset_idx = filter_small_highSR.filtered_vcf_idx
    File large_lowSR_minGQ_ROC_plot = ROC_large_lowSR.minGQ_ROC_plot
    File large_lowSR_minGQ_ROC_table = ROC_large_lowSR.minGQ_ROC_table
    File large_lowSR_vcf_subset = filter_large_lowSR.filtered_vcf
    File large_lowSR_vcf_subset_idx = filter_large_lowSR.filtered_vcf_idx
    File large_highSR_minGQ_ROC_plot = ROC_large_highSR.minGQ_ROC_plot
    File large_highSR_minGQ_ROC_table = ROC_large_highSR.minGQ_ROC_table
    File large_highSR_vcf_subset = filter_large_highSR.filtered_vcf
    File large_highSR_vcf_subset_idx = filter_large_highSR.filtered_vcf_idx
  }
}


#Shard famfile
task split_famfile {
  File famfile
  File vcf
  Int fams_per_shard

  command <<<
    #Get list of sample IDs & column numbers from VCF header
    zcat ${vcf} | head -n1000 | fgrep "#" | fgrep -v "##" | sed 's/\t/\n/g' \
      | awk -v OFS="\t" '{ print $1, NR }' > vcf_header_columns.txt
    #Iterate over families & subset VCF
    while read famID pro fa mo prosex pheno; do
      pro_idx=$( awk -v ID=$pro '{ if ($1==ID) print $2 }' vcf_header_columns.txt )
      fa_idx=$( awk -v ID=$fa '{ if ($1==ID) print $2 }' vcf_header_columns.txt )
      mo_idx=$( awk -v ID=$mo '{ if ($1==ID) print $2 }' vcf_header_columns.txt )
      if ! [ -z $pro_idx ] && ! [ -z $fa_idx ] && ! [ -z $mo_idx ]; then
        fgrep -w "$famID" ${famfile}
      fi
    done < ${famfile} \
      | split -l ${fams_per_shard} --numeric-suffixes=001 -a 3 /dev/stdin famfile_shard_
  >>>

  output {
    Array[File] famfile_shards = glob("famfile_shard_*")
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:d3844f6c7c26da55e679c9c521882d54dbecd169f884f09f05a12d6565bf6063"
    disks: "local-disk 250 HDD"
    preemptible: 1
  }
}


#Merge vcf2bed shards
task merge_vcf2bed {
  Array[File] vcf2bed_shards
  String prefix

  command <<<
    zcat ${vcf2bed_shards[0]} | sed -n '1p' > header.txt
    zcat ${sep=' ' vcf2bed_shards} | fgrep -v "#" \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | cat header.txt - \
      | bgzip -c \
      > "${prefix}.vcf2bed.bed.gz"
  >>>

  output {
    File merged_vcf2bed_out = "${prefix}.vcf2bed.bed.gz"
  }
  
  runtime {
    preemptible: 1
    docker: "talkowski/sv-pipeline@sha256:d3844f6c7c26da55e679c9c521882d54dbecd169f884f09f05a12d6565bf6063"
    memory: "4 GB"
    disks: "local-disk 200 SSD"
  }
}


#Filter VCF based on size and SR background rate
task filter_vcf {
  File vcf
  File bed
  String prefix
  Int min_size
  Int max_size
  String highSR

  command <<<
    #Get list of variants to keep
    size_idx=$( zcat ${bed} | sed -n '1p' | sed 's/\t/\n/g' \
                  | awk '{ if ($1=="SVLEN") print NR }' )
    if [ ${highSR} == "TRUE" ]; then
      zcat ${bed} \
        | fgrep -v "#" \
        | fgrep -w "HIGH_SR_BACKGROUND" \
        | awk -v OFS="\t" -v FS="\t" \
          -v size_idx="$size_idx" \
          -v min_size=${min_size} \
          -v max_size=${max_size} \
          '{ if ($(size_idx) >= min_size && $(size_idx) < max_size) print $4 }' \
        > VIDs_to_keep.list
    else
      zcat ${bed} \
        | fgrep -v "#" \
        | fgrep -wv "HIGH_SR_BACKGROUND" \
        | awk -v OFS="\t" -v FS="\t" \
          -v size_idx="$size_idx" \
          -v min_size=${min_size} \
          -v max_size=${max_size} \
          '{ if ($(size_idx) >= min_size && $(size_idx) < max_size) print $4 }' \
        > VIDs_to_keep.list
    fi

    #Filter vcf
    zcat ${vcf} | sed -n '1,2000p' | fgrep "#" > header.vcf
    zcat ${vcf} | fgrep -wf VIDs_to_keep.list | cat header.vcf - | bgzip -c \
    > "${prefix}.filtered.vcf.gz"
    tabix -p vcf "${prefix}.filtered.vcf.gz"
  >>>

  output {
    File filtered_vcf = "${prefix}.filtered.vcf.gz"
    File filtered_vcf_idx = "${prefix}.filtered.vcf.gz.tbi"
  }
  
  runtime {
    preemptible: 1
    docker: "talkowski/sv-pipeline@sha256:d3844f6c7c26da55e679c9c521882d54dbecd169f884f09f05a12d6565bf6063"
    memory: "4 GB"
    disks: "local-disk 200 SSD"
  }
}

