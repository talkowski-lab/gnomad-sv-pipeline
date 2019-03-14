import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04_genotype_CPX_CNVs_perBatch/versions/19/plain-WDL/descriptor" as rd_gt_perbatch

# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# Workflow to perform depth-based genotyping for a single vcf shard scattered 
# across batches on predicted CPX CNVs from 04b
workflow genotype_CPX_CNVs {
  File vcf
  File gt_input_files
  Int n_per_split_small
  Int n_per_split_large
  Int n_RdTest_bins
  File svc_acct_key
  String prefix
  File famfile
  String contig

  Array[Array[String]] gt_input_array = read_tsv(gt_input_files)

  # Convert VCF to bed of CPX CNV intervals
  call get_cpx_cnv_intervals {
    input:
      vcf=vcf,
      prefix="${prefix}.${contig}"
  }

  # Scatter over each batch (row) in gt_input_files and run depth genotyping
  scatter (gt_inputs in gt_input_array) {
    call rd_gt_perbatch.genotype_CPX_CNVs_perBatch as gt_batch {
      input:
        cpx_bed=get_cpx_cnv_intervals.CPX_CNV_BED,
        batch=gt_inputs[0],
        coveragefile=gt_inputs[1],
        coveragefile_idx=gt_inputs[2],
        RD_depth_sepcutoff=gt_inputs[3],
        sampleslist=gt_inputs[4],
        famfile=gt_inputs[5],
        medianfile=gt_inputs[6],
        n_per_split_small=n_per_split_small,
        n_per_split_large=n_per_split_large,
        n_RdTest_bins=n_RdTest_bins,
        svc_acct_key=svc_acct_key
    }
  }

  # Merge melted genotypes across all batches
  call merge_melted_gts {
    input:
      melted_gts=gt_batch.genotypes,
      prefix="${prefix}.${contig}"
  }

  # Parse genotyping results
  call parse_gts {
    input:
      vcf=vcf,
      intervals=get_cpx_cnv_intervals.CPX_CNV_BED,
      genotypes=merge_melted_gts.merged_genotypes,
      prefix="${prefix}.${contig}",
      famfile=famfile,
      contig=contig
  }

  # Final output
  output {
    File cpx_depth_gt_resolved_vcf = parse_gts.cpx_depth_gt_resolved_vcf
    File reclassification_table = parse_gts.reclassification_table
    File interval_genotype_counts_table = parse_gts.gt_counts_table
  }
}


# Get CNV intervals from complex SV for depth genotyping
task get_cpx_cnv_intervals {
  File vcf
  String prefix

  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/gather_cpx_intervals_for_rd_gt.sh \
      ${vcf} \
      ${prefix}.complex_CNV_intervals_to_test.bed.gz
  >>>

  output {
    File CPX_CNV_BED = "${prefix}.complex_CNV_intervals_to_test.bed.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:c2af5febc8967dff0b7a10cd764b292f43029ffd119e40832cef3fcbc3df1c1f"
    preemptible: 1
    memory: "8 GB"
    disks: "local-disk 100 HDD"
  }
}


# Merge output from per-batch genotyping
task merge_melted_gts {
  Array[File] melted_gts
  String prefix

  command <<<
    while read file; do
      zcat $file
    done < ${write_tsv(melted_gts)} \
      | sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V \
      | bgzip -c \
      > ${prefix}.CPX_intervals.merged_rd_genos.bed.gz
  >>>

  output {
    File merged_genotypes = "${prefix}.CPX_intervals.merged_rd_genos.bed.gz"
  }

    runtime {
    docker: "talkowski/sv-pipeline@sha256:c2af5febc8967dff0b7a10cd764b292f43029ffd119e40832cef3fcbc3df1c1f"
    preemptible: 1
    disks: "local-disk 100 HDD"
  }
}


# Parse genotyping results
task parse_gts {
  File vcf
  File intervals
  File genotypes
  File famfile
  String prefix
  String contig

  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/process_posthoc_cpx_depth_regenotyping.sh \
      -R ${prefix}.CPXregenotyping_reclassification_table.${contig}.txt \
      -G ${prefix}.CPXregenotyping_raw_genotype_counts_table.${contig}.txt \
      ${vcf} \
      ${intervals} \
      ${genotypes} \
      ${famfile} \
      ${prefix}.postCPXregenotyping.${contig}.vcf.gz
  >>>

  output {
    File cpx_depth_gt_resolved_vcf = "${prefix}.postCPXregenotyping.${contig}.vcf.gz"
    File reclassification_table = "${prefix}.CPXregenotyping_reclassification_table.${contig}.txt"
    File gt_counts_table = "${prefix}.CPXregenotyping_raw_genotype_counts_table.${contig}.txt"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:c2af5febc8967dff0b7a10cd764b292f43029ffd119e40832cef3fcbc3df1c1f"
    preemptible: 1
    disks: "local-disk 100 HDD"
  }
}


# Combine multiple VCFs
task concat_vcfs {
  Array[File] vcfs
  String prefix

  command <<<
    vcf-concat ${sep=' ' vcfs} | vcf-sort -c | bgzip -c > ${prefix}.mod04b_final.vcf.gz
  >>>

  output {
    File concat_vcf = "${prefix}.mod04b_final.vcf.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:c2af5febc8967dff0b7a10cd764b292f43029ffd119e40832cef3fcbc3df1c1f"
    preemptible: 1
    disks: "local-disk 300 HDD"
  }
}
