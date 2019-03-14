# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# Workflow to perform final sample pruning & compute all relevant AF statistics
# for a VCF from the Talkowski SV pipeline

import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:compute_simple_AFs_singleChrom/versions/10/plain-WDL/descriptor" as calcAF

workflow prune_and_add_vafs {
  File vcf
  File vcf_idx
  String prefix
  File? sample_pop_assignments  #Two-column file with sample ID & pop assignment. "." for pop will ignore sample
  File? prune_list              #List of samples to be excluded from the output vcf
  File? famfile                 #Used for M/F AF calculations
  Int sv_per_shard
  File contiglist

  Array[Array[String]] contigs=read_tsv(contiglist)


  #Iterate over chromosomes
  scatter (contig in contigs) {
    #Prune VCF
    call prune_vcf {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        contig=contig[0],
        prune_list=prune_list,
        prefix=prefix
    }
    #Compute AC, AN, and AF per population & sex combination
    call calcAF.getAFs_singleChrom as getAFs {
      input:
        vcf=prune_vcf.pruned_vcf,
        vcf_idx=prune_vcf.pruned_vcf_idx,
        contig=contig[0],
        sv_per_shard=sv_per_shard,
        prefix=prefix,
        sample_pop_assignments=sample_pop_assignments,
        famfile=famfile
    }
  }


  #Merge pruned VCFs with allele info
  call concat_vcfs {
    input:
      vcfs=getAFs.vcf_wAFs,
      outfile_prefix="${prefix}.pruned_wAFs"
  }

  output {
    File output_vcf = concat_vcfs.concat_vcf
    File output_vcf_idx = concat_vcfs.concat_vcf_idx
  }
}


#Shard vcf into single chromosome shards & drop pruned samples
task prune_vcf {
  File vcf
  File vcf_idx
  String contig
  File? prune_list
  String prefix

  command <<<
    #Tabix chromosome of interest
    tabix -h ${vcf} ${contig} | bgzip -c > ${contig}.vcf.gz
    #Get column indexes corresponding to samples to drop, if any exist
    if [ "${default="SKIP" prune_list}" != "SKIP" ]; then
      dropidx=$( zcat ${contig}.vcf.gz | sed -n '1,500p' | fgrep "#" | fgrep -v "##" \
                 | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print NR, $1 }' \
                 | fgrep -wf ${prune_list} | cut -f1 | paste -s -d, )
      zcat ${contig}.vcf.gz \
      | cut --complement -f"$dropidx" \
      | bgzip -c \
      > "${prefix}.${contig}.pruned.vcf.gz"
    else
      cp "${contig}.vcf.gz" "${prefix}.${contig}.pruned.vcf.gz"
    fi
    tabix -f "${prefix}.${contig}.pruned.vcf.gz"
  >>>

  output {
    File pruned_vcf = "${prefix}.${contig}.pruned.vcf.gz"
    File pruned_vcf_idx = "${prefix}.${contig}.pruned.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:831595263ca60288fa8512602d1a3f1fcc23c3f31a6a8f0db2e597138b5e3d36"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 250 SSD"
  }
}


#General task to combine multiple VCFs
task concat_vcfs {
  Array[File] vcfs
  String outfile_prefix

  command <<<
    vcf-concat ${sep=' ' vcfs} | vcf-sort -c | bgzip -c > ${outfile_prefix}.vcf.gz; 
    tabix -p vcf -f "${outfile_prefix}.vcf.gz"
  >>>

  output {
    File concat_vcf = "${outfile_prefix}.vcf.gz"
    File concat_vcf_idx = "${outfile_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:7f1e8ae2c7ce7779fc3de75a67841fde21bbd5be657911abd26e8551aba9e8a5"
    preemptible: 1
    memory: "16 GB"
    disks: "local-disk 250 SSD"
  }
}


