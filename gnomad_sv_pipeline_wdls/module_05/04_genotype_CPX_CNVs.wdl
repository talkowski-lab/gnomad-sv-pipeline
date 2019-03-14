import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04_genotype_CPX_CNVs_perBatch/versions/18/plain-WDL/descriptor" as rd_gt_perbatch

# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# Workflow to perform depth-based genotyping scattered across batches
# on predicted CPX CNVs from 04b

workflow genotype_CPX_CNVs {
  File vcf
  File gt_input_files
  Int n_per_split_small
  Int n_per_split_large
  Int n_RdTest_bins
  File svc_acct_key
  String prefix
  File famfile
  File contigs

  Array[Array[String]] gt_input_array = read_tsv(gt_input_files)
  Array[Array[String]] contiglist = read_tsv(contigs)

  # Convert VCF to bed of CPX CNV intervals
  call get_cpx_cnv_intervals {
    input:
      vcf=vcf,
      prefix=prefix
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
      prefix=prefix
  }

  # Scatter over contigs and parse genotyping results
  scatter (contig in contiglist) {
    # Prep files per contig
    call prep_chrom_files {
      input:
        vcf=vcf,
        intervals=get_cpx_cnv_intervals.CPX_CNV_BED,
        genotypes=merge_melted_gts.merged_genotypes,
        prefix=prefix,
        contig=contig[0]
    }

    # Parse genotyping results
    call parse_gts {
      input:
        vcf=prep_chrom_files.vcf_shard,
        intervals=prep_chrom_files.intervals_shard,
        genotypes=prep_chrom_files.genotypes_shard,
        prefix=prefix,
        famfile=famfile,
        contig=contig[0]
    }
  }

  # Merge VCF shards across contigs
  call concat_vcfs {
    input:
      vcfs=parse_gts.cpx_depth_gt_resolved_vcf,
      prefix=prefix
  }

  # Final output
  output {
    File cpx_depth_gt_resolved_vcf = concat_vcfs.concat_vcf
    # File CPX_interval_genotypes = merge_melted_gts.merged_genotypes
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
    docker: "talkowski/sv-pipeline@sha256:5762d98ff5c88c0f2380b7de39f309316165b70880998022602432579da42344"
    preemptible: 3
    disks: "local-disk 250 SSD"
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
    docker: "talkowski/sv-pipeline@sha256:5762d98ff5c88c0f2380b7de39f309316165b70880998022602432579da42344"
    preemptible: 3
    disks: "local-disk 250 SSD"
  }
}


# Prep per-chromosome files for gt result parsing
task prep_chrom_files {
  File vcf
  File intervals
  File genotypes
  String prefix
  String contig

  command <<<
    tabix -f ${vcf};
    tabix -f ${intervals};
    tabix -f ${genotypes};
    zcat ${intervals} \
      | head -n1 \
      > ${prefix}.complex_CNV_intervals_to_test.${contig}.bed;
    tabix ${vcf} ${contig} \
      | cut -f3 \
      | fgrep -wf - <( zcat ${intervals} ) \
      | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
      >> ${prefix}.complex_CNV_intervals_to_test.${contig}.bed;
    bgzip -f ${prefix}.complex_CNV_intervals_to_test.${contig}.bed;
    zcat ${genotypes} \
      | head -n1 \
      > ${prefix}.CPX_intervals.merged_rd_genos.${contig}.bed;
    zcat ${prefix}.complex_CNV_intervals_to_test.${contig}.bed.gz \
      | fgrep -v "#" | cut -f4 \
      | fgrep -wf - <( zcat ${genotypes} ) \
      >> ${prefix}.CPX_intervals.merged_rd_genos.${contig}.bed;
    bgzip -f ${prefix}.CPX_intervals.merged_rd_genos.${contig}.bed;
    tabix -h ${vcf} ${contig} | vcf-sort | bgzip -c \
      > ${prefix}.sharded_vcf.${contig}.vcf.gz
  >>>

  output {
    File vcf_shard = "${prefix}.sharded_vcf.${contig}.vcf.gz"
    File intervals_shard = "${prefix}.complex_CNV_intervals_to_test.${contig}.bed.gz"
    File genotypes_shard = "${prefix}.CPX_intervals.merged_rd_genos.${contig}.bed.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:5762d98ff5c88c0f2380b7de39f309316165b70880998022602432579da42344"
    preemptible: 3
    disks: "local-disk 250 SSD"
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
      ${vcf} \
      ${intervals} \
      ${genotypes} \
      ${famfile} \
      ${prefix}.postCPXregenotyping.${contig}.vcf.gz
  >>>

  output {
    File cpx_depth_gt_resolved_vcf = "${prefix}.postCPXregenotyping.${contig}.vcf.gz"
    File reclassification_table = "${prefix}.CPXregenotyping_reclassification_table.${contig}.txt"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:7c59e97cf19f9d2b203b5e1d2e77722dc61c49770a73c560805d82d3e622d053"
    preemptible: 3
    disks: "local-disk 250 SSD"
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
    docker: "talkowski/sv-pipeline@sha256:5762d98ff5c88c0f2380b7de39f309316165b70880998022602432579da42344"
    preemptible: 3
    disks: "local-disk 250 SSD"
  }
}
