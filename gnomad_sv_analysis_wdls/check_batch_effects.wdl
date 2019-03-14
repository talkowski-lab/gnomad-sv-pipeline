# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# This is an analysis WDL to perform pairwise comparisons of batches in the 
# Talkowski lab SV pipeline, and mark sites that appear batch-specific


import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:prune_and_add_vfs/versions/13/plain-WDL/descriptor" as calcAF
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:gather_batch_effects_helper/versions/8/plain-WDL/descriptor" as helper


workflow check_batch_effects {
	File vcf
  File vcf_idx
  File sample_batch_assignments
  File batches_list
  File sample_pop_assignments
  File probands_list
  File PCRPLUS_samples_list
  File famfile
  File contiglist
  Int variants_per_shard
  String prefix

  Array[Array[String]] batches = read_tsv(batches_list)

  # Shard VCF per batch, compute pops-specific AFs, and convert to table of VID & AF stats
  scatter ( batch in batches ) {
    # Get list of samples to include & exclude per batch
    call get_batch_samples_list {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        batch=batch[0],
        sample_batch_assignments=sample_batch_assignments,
        probands_list=probands_list
    }
    # Prune VCF to samples 
    call calcAF.prune_and_add_vafs as getAFs {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        prefix=batch[0],
        sample_pop_assignments=sample_pop_assignments,
        prune_list=get_batch_samples_list.exclude_samples_list,
        famfile=famfile,
        sv_per_shard=25000,
        contiglist=contiglist
    }
    # Get minimal table of AF data per batch
    call get_freq_table {
      input:
        vcf=getAFs.output_vcf,
        prefix=batch[0]
    }
  }

  # Merge frequency results per batch into a single table of all variants with AF data across batches
  call merge_freq_tables {
    input:
      tables=get_freq_table.freq_data,
      batches_list=batches_list,
      prefix=prefix
  }

  # Generate matrix of correlation coefficients for all batches, by population & SVTYPE
  Array[String] populations = ["AFR", "ASN", "EUR", "HSP"]
  scatter ( pop in populations ) {
    call make_correlation_matrices {
      input:
        freq_table=merge_freq_tables.merged_table,
        pop=pop,
        batches_list=batches_list,
        prefix=prefix
    }
  }

  # Make list of nonredundant pairs of batches to be evaluated
  call make_batch_pairs_list {
    input:
      batches_list=batches_list,
      prefix=prefix
  }
  Array[Array[String]] batch_pairs = read_tsv(make_batch_pairs_list.batch_pairs_list)

  # Compute AF stats per pair of batches & determine variants with batch effects
  scatter ( pair in batch_pairs ) {
    call helper.check_batch_effects as check_batch_effects {
      input:
        freq_table=merge_freq_tables.merged_table,
        batch1=pair[0],
        batch2=pair[1],
        prefix=prefix,
        variants_per_shard=variants_per_shard
    }
  }
  # Collect results from pairwise batch effect detection
  call merge_variant_failure_lists as merge_pairwise_checks {
    input:
      fail_variant_lists=check_batch_effects.batch_effect_variants,
      prefix="${prefix}.pairwise_comparisons"
  }

  # Perform one-vs-all comparison of AFs per batch to find batch-specific sites
  scatter ( batch in batches ) {
    call helper.check_batch_effects as one_vs_all_comparison {
      input:
        freq_table=merge_freq_tables.merged_table,
        batch1=batch[0],
        batch2="ALL_OTHERS",
        prefix=prefix,
        variants_per_shard=variants_per_shard
    }
  }
  # Collect results from pairwise batch effect detection
  call merge_variant_failure_lists as merge_one_vs_all_checks {
    input:
      fail_variant_lists=one_vs_all_comparison.batch_effect_variants,
      prefix="${prefix}.one_vs_all_comparisons"
  }

  # Distill final table of variants to be reclassified
  call make_reclassification_table {
    input:
      freq_table=merge_freq_tables.merged_table,
      pairwise_fails=merge_pairwise_checks.fails_per_variant_all,
      pairwise_PCRMINUS_fails=merge_pairwise_checks.fails_per_variant_PCRMINUS_to_PCRMINUS,
      onevsall_fails=merge_one_vs_all_checks.fails_per_variant_all,
      onevsall_PCRMINUS_fails=merge_one_vs_all_checks.fails_per_variant_PCRMINUS_to_PCRMINUS,
      prefix=prefix
  }

  # Apply batch effect labels
  call apply_batch_effect_labels {
    input:
      vcf=vcf,
      vcf_idx=vcf_idx,
      reclassification_table=make_reclassification_table.reclassification_table,
      PCRPLUS_samples_list=PCRPLUS_samples_list,
      prefix=prefix
  }

  output {
    File labeled_vcf = apply_batch_effect_labels.labeled_vcf
    File labeled_vcf_idx = apply_batch_effect_labels.labeled_vcf_idx
  }
}


# Get list of samples to include & exclude per batch
# Always exclude probands from all batches
task get_batch_samples_list {
  File vcf
  File vcf_idx
  String batch
  File sample_batch_assignments
  File probands_list

  command <<<
    # Get list of all samples present in VCF header
    tabix -H ${vcf} | fgrep -v "##" | cut -f10- | sed 's/\t/\n/g' | sort -Vk1,1 \
    > all_samples.list
    # Get list of samples in batch
    fgrep -w ${batch} ${sample_batch_assignments} | cut -f1 \
    | fgrep -wf - all_samples.list \
    | fgrep -wvf ${probands_list} \
    > "${batch}.samples.list"
    # Get list of samples not in batch
    fgrep -wv ${batch} ${sample_batch_assignments} | cut -f1 \
    cat - ${probands_list} | sort -Vk1,1 | uniq \
    | fgrep -wf - all_samples.list \
    > "${batch}.exclude_samples.list"
  >>>

  output {
    File include_samples_list = "${batch}.samples.list"
    File exclude_samples_list = "${batch}.exclude_samples.list"
  }

  runtime {
    preemptible: 1
    docker : "talkowski/sv-pipeline@sha256:58b67cb4e4edf285b89250d2ebab72e17c0247e3bf6891c2c2fcda646b2a6cf4"
    disks: "local-disk 50 SSD"
  }
}


# Run vcf2bed and subset to just include VID, SVTYPE, SVLEN, _AC, and _AN
task get_freq_table {
  File vcf
  String prefix

  command <<<
    #Run vcf2bed
    svtk vcf2bed \
      --info ALL \
      --no-samples \
      ${vcf} "${prefix}.vcf2bed.bed"
    #Cut to necessary columns
    idxs=$( sed -n '1p' "${prefix}.vcf2bed.bed" \
            | sed 's/\t/\n/g' \
            | awk -v OFS="\t" '{ print $1, NR }' \
            | grep -e 'name\|SVLEN\|SVTYPE\|_AC\|_AN' \
            | fgrep -v "OTH" \
            | cut -f2 \
            | paste -s -d\, )
    cut -f"$idxs" "${prefix}.vcf2bed.bed" \
    | sed 's/^name/\#VID/g' \
    | gzip -c \
    > "${prefix}.frequencies.preclean.txt.gz"
    #Clean frequencies
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/clean_frequencies_table.R \
      "${prefix}.frequencies.preclean.txt.gz" \
      "${prefix}.frequencies.txt"
  >>>

  output {
    File freq_data = "${prefix}.frequencies.txt.gz"
    # File freqs_preclean = "${prefix}.frequencies.preclean.txt.gz"
  }

  runtime {
    docker : "talkowski/sv-pipeline@sha256:5a5bb420b823b129a15ca49ac2e6d3862984c5dd1b8ec91a49ee1ba803685d83"
    disks: "local-disk 50 HDD"
    memory: "4 GB"
    preemptible: 1
  }
}


# Combine frequency data across batches
task merge_freq_tables {
  Array[File] tables
  File batches_list
  String prefix

  command <<<
    #Get list of batch IDs and batch table paths
    while read batch; do
      echo "$batch"
      find / -name "$batch.frequencies.txt.gz"
    done < ${batches_list} | paste - - \
    > input.list
    #Make sure all input files have the same number of lines
    nlines=$( while read batch file; do
                zcat "$file" | wc -l
              done < input.list | sort | uniq | wc -l )
    if [ "$nlines" -gt 1 ]; then
      echo "AT LEAST ONE INPUT FILE HAS A DIFFERENT NUMBER OF LINES"
      exit 0
    fi
    #Prep files for paste joining
    echo "PREPPING FILES FOR MERGING"
    while read batch file; do
        #Header
        zcat "$file" | sed -n '1p' | cut -f1-3
        #Body
        zcat "$file" | sed '1d' \
        | sort -Vk1,1 \
        | cut -f1-3
    done < <( sed -n '1p' input.list ) \
    > header.txt
    while read batch file; do
      for wrapper in 1; do
        #Header
        zcat "$file" | sed -n '1p' \
        | cut -f4- | sed 's/\t/\n/g' \
        | awk -v batch="$batch" '{ print $1"."batch }' \
        | paste -s
        #Body
        zcat "$file" | sed '1d' \
        | sort -Vk1,1 \
        | cut -f4-
      done > "$batch.prepped.txt" 
    done < input.list
    #Join files with simple paste
    paste \
      header.txt \
      $( awk -v ORS=" " '{ print $1".prepped.txt" }' input.list ) \
    | gzip -c \
    > "${prefix}.merged_AF_table.txt.gz"
    # #Make master table (OLD IMPLEMENTATION -- TOO MEMORY INEFFICIENT)
    # /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/merge_batch_freq_tables.R \
    #   input.list \
    #   "${prefix}.merged_AF_table.txt"
    # gzip "${prefix}.merged_AF_table.txt"
  >>>

  output {
    File merged_table = "${prefix}.merged_AF_table.txt.gz"
  }

  runtime {
    docker : "talkowski/sv-pipeline@sha256:02e52b8a3f158ee2a3e2c299385a92d3f75baf6aae2a810e79b39d343887f8da"
    disks: "local-disk 100 HDD"
    memory: "16 GB"
    preemptible: 1
  }
}


# Calculate & plot cross-batch correlation coefficient matrixes
task make_correlation_matrices {
  File freq_table
  String pop
  File batches_list
  String prefix

  command <<<
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/correlate_batches_singlePop.R \
      ${batches_list} \
      ${freq_table} \
      "${pop}" \
      "${prefix}.${pop}"
  >>>

  output {
    Array[File] corr_matrixes = glob("${prefix}.${pop}.*.R2_matrix.txt")
    Array[File] heatmaps = glob("${prefix}.${pop}.*heatmap*.pdf")
    Array[File] dotplots = glob("${prefix}.${pop}.*perBatch_R2_sina_plot.pdf")
  }

  runtime {
    docker : "talkowski/sv-pipeline@sha256:f4a86948e52c171cc7e58034e3d8bc4819a0d6932925fd4874bca6a91baa5c3c"
    disks: "local-disk 50 HDD"
    memory: "8 GB"
    preemptible: 1
  }
}


# Generate list of all pairs of batches to be compared
task make_batch_pairs_list {
  File batches_list
  String prefix

  command <<<
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/make_batch_pairs_list.R \
      ${batches_list} \
      "${prefix}.nonredundant_batch_pairs.txt"
  >>>

  output {
    File batch_pairs_list = "${prefix}.nonredundant_batch_pairs.txt"
  }

  runtime {
    docker : "talkowski/sv-pipeline@sha256:e06c6701a50be7513e3eae7397d3a5d9290fe0dd9a457e574de2f58f5d337b69"
    preemptible: 1
  }
}


# Merge lists of batch effect checks and count total number of times each variant failed
task merge_variant_failure_lists {
  Array[File] fail_variant_lists
  String prefix

  command <<<
    #Write list of paths to all batch effect variant lists
    while read path; do
      echo "$path"
    done < ${write_lines(fail_variant_lists)} \
    > fail_variant_lists.paths.txt
    #Get master list of PCR+ to PCR+ failures
    fgrep "PCRPLUS" fail_variant_lists.paths.txt \
    | fgrep -v "PCRMINUS" \
    | xargs -I {} cat {} \
    | sort -Vk1,1 | uniq -c \
    | awk -v OFS="\t" '{ print $2, $1 }' \
    > "${prefix}.PCRPLUS_to_PCRPLUS.failures.txt"
    #Get master list of PCR- to PCR- failures
    fgrep -v "PCRPLUS" fail_variant_lists.paths.txt \
    | fgrep "PCRMINUS" \
    | xargs -I {} cat {} \
    | sort -Vk1,1 | uniq -c \
    | awk -v OFS="\t" '{ print $2, $1 }' \
    > "${prefix}.PCRMINUS_to_PCRMINUS.failures.txt"
    #Get master list of PCR+ to PCR- failures
    fgrep "PCRPLUS" fail_variant_lists.paths.txt \
    | fgrep "PCRMINUS" \
    | xargs -I {} cat {} \
    | sort -Vk1,1 | uniq -c \
    | awk -v OFS="\t" '{ print $2, $1 }' \
    > "${prefix}.PCRPLUS_to_PCRMINUS.failures.txt"
    #Get master list of all possible failures
    cat fail_variant_lists.paths.txt \
    | xargs -I {} cat {} \
    | sort -Vk1,1 | uniq -c \
    | awk -v OFS="\t" '{ print $2, $1 }' \
    > "${prefix}.all.failures.txt"
  >>>

  output {
    File fails_per_variant_PCRPLUS_to_PCRPLUS = "${prefix}.PCRPLUS_to_PCRPLUS.failures.txt"
    File fails_per_variant_PCRMINUS_to_PCRMINUS = "${prefix}.PCRMINUS_to_PCRMINUS.failures.txt"
    File fails_per_variant_PCRPLUS_to_PCRMINUS = "${prefix}.PCRPLUS_to_PCRMINUS.failures.txt"
    File fails_per_variant_all = "${prefix}.all.failures.txt"
  }

  runtime {
    docker : "talkowski/sv-pipeline@sha256:e06c6701a50be7513e3eae7397d3a5d9290fe0dd9a457e574de2f58f5d337b69"
    preemptible: 1
  }
}


# Consolidate all batch effect check results into a single table with reclassification per variant
task make_reclassification_table {
  File freq_table
  File pairwise_fails
  File pairwise_PCRMINUS_fails
  File onevsall_fails
  File onevsall_PCRMINUS_fails
  String prefix

  command <<<
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/make_batch_effect_reclassification_table.R \
      ${freq_table} \
      ${pairwise_fails} \
      ${pairwise_PCRMINUS_fails} \
      ${onevsall_fails} \
      ${onevsall_PCRMINUS_fails} \
      "${prefix}.batch_effect_reclassification_table.txt"
  >>>

  output {
    File reclassification_table = "${prefix}.batch_effect_reclassification_table.txt"
  }

  runtime {
    docker : "talkowski/sv-pipeline@sha256:a182bde10e0540e915ec3c0ef2c89b281e392c962ccc50f7a1c4503c20698ae4"
    preemptible: 1
    memory: "8 GB"
  }
}


# Apply batch effect labels to VCF
task apply_batch_effect_labels {
  File vcf
  File vcf_idx
  File reclassification_table
  File PCRPLUS_samples_list
  String prefix

  command <<<
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/label_batch_effects.py \
    ${vcf} \
    ${reclassification_table} \
    ${PCRPLUS_samples_list} \
    stdout \
    | bgzip -c \
    > "${prefix}.batch_effects_labeled.vcf.gz"
    tabix -p vcf -f "${prefix}.batch_effects_labeled.vcf.gz"
  >>>

  output {
    File labeled_vcf = "${prefix}.batch_effects_labeled.vcf.gz"
    File labeled_vcf_idx = "${prefix}.batch_effects_labeled.vcf.gz.tbi"
  }

  runtime {
    docker : "talkowski/sv-pipeline@sha256:a182bde10e0540e915ec3c0ef2c89b281e392c962ccc50f7a1c4503c20698ae4"
    preemptible: 0
    disks: "local-disk 50 HDD"
    memory: "4 GB"
  }
}


