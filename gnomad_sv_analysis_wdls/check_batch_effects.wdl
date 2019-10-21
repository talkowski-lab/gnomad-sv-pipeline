# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# This is an analysis WDL to perform pairwise comparisons of batches in the 
# Talkowski lab SV pipeline, and mark sites that appear batch-specific


import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:prune_and_add_vfs/versions/17/plain-WDL/descriptor" as calcAF
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:gather_batch_effects_helper/versions/10/plain-WDL/descriptor" as helper


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
  File AF_PCRPLUS_preMinGQ
  File AF_PCRMINUS_preMinGQ


  Array[Array[String]] batches = read_tsv(batches_list)
  Array[Array[String]] contigs = read_tsv(contiglist)

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
        contiglist=contiglist,
        drop_empty_records="FALSE"
    }
    # Get minimal table of AF data per batch, split by ancestry
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
  call merge_freq_tables as merge_freq_tables_allPops {
    input:
      tables=get_freq_table.freq_data_allPops,
      batches_list=batches_list,
      prefix=prefix
  }

  # Compare frequencies before and after minGQ, and generate list of variants
  # that are significantly different between the steps
  call compare_freqs_prePost_minGQ {
    input:
      AF_PCRPLUS_preMinGQ=AF_PCRPLUS_preMinGQ,
      AF_PCRMINUS_preMinGQ=AF_PCRMINUS_preMinGQ,
      AF_postMinGQ_table=merge_freq_tables_allPops.merged_table,
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
  scatter ( contig in contigs ) {
    call apply_batch_effect_labels as apply_labels_perContig {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        contig=contig[0],
        reclassification_table=make_reclassification_table.reclassification_table,
        PCRPLUS_samples_list=PCRPLUS_samples_list,
        minGQ_prePost_PCRPLUS_fails=compare_freqs_prePost_minGQ.PCRPLUS_fails,
        minGQ_prePost_PCRMINUS_fails=compare_freqs_prePost_minGQ.PCRMINUS_fails,
        prefix="${prefix}.${contig[0]}"
    }
  }
  call concat_vcfs as merge_labeled_vcfs {
    input:
      vcfs=apply_labels_perContig.labeled_vcf,
      outfile_prefix="${prefix}.batch_effects_labeled_merged"
  }

  output {
    File labeled_vcf = merge_labeled_vcfs.concat_vcf
    File labeled_vcf_idx = merge_labeled_vcfs.concat_vcf_idx
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
    set -euo pipefail
    # Get list of all samples present in VCF header
    tabix -H ${vcf} | fgrep -v "##" | cut -f10- | sed 's/\t/\n/g' | sort -Vk1,1 \
    > all_samples.list
    # Get list of samples in batch
    fgrep -w ${batch} ${sample_batch_assignments} | cut -f1 \
    | fgrep -wf - all_samples.list \
    | fgrep -wvf ${probands_list} \
    > "${batch}.samples.list" || true
    # Get list of samples not in batch
    fgrep -wv ${batch} ${sample_batch_assignments} | cut -f1 \
    cat - ${probands_list} | sort -Vk1,1 | uniq \
    | fgrep -wf - all_samples.list \
    > "${batch}.exclude_samples.list" || true
  >>>

  output {
    File include_samples_list = "${batch}.samples.list"
    File exclude_samples_list = "${batch}.exclude_samples.list"
  }

  runtime {
    preemptible: 1
    maxRetries: 1
    docker : "talkowski/sv-pipeline@sha256:aef8156983cec6ac6a91fa6461b197a63835e5487fc9523ec857f947cfac660e"
    disks: "local-disk 50 SSD"
  }
}


# Run vcf2bed and subset to just include VID, SVTYPE, SVLEN, _AC, and _AN
task get_freq_table {
  File vcf
  String prefix

  command <<<
    set -euo pipefail
    #Run vcf2bed
    svtk vcf2bed \
      --info ALL \
      --no-samples \
      ${vcf} "${prefix}.vcf2bed.bed"
    ### Create table of freqs by ancestry
    #Cut to necessary columns
    idxs=$( sed -n '1p' "${prefix}.vcf2bed.bed" \
            | sed 's/\t/\n/g' \
            | awk -v OFS="\t" '{ print $1, NR }' \
            | grep -e 'name\|SVLEN\|SVTYPE\|_AC\|_AN' \
            | fgrep -v "OTH" \
            | cut -f2 \
            | paste -s -d\, || true )
    cut -f"$idxs" "${prefix}.vcf2bed.bed" \
    | sed 's/^name/\#VID/g' \
    | gzip -c \
    > "${prefix}.frequencies.preclean.txt.gz"
    #Clean frequencies
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/clean_frequencies_table.R \
      "${prefix}.frequencies.preclean.txt.gz" \
      "${prefix}.frequencies.txt"
    ### Create table of freqs, irrespective of ancestry
    #Cut to necessary columns
    idxs=$( sed -n '1p' "${prefix}.vcf2bed.bed" \
            | sed 's/\t/\n/g' \
            | awk -v OFS="\t" '{ if ($1=="name" || $1=="SVLEN" || $1=="SVTYPE" || $1=="AC" || $1=="AN") print NR }' \
            | paste -s -d\, || true )
    cut -f"$idxs" "${prefix}.vcf2bed.bed" \
    | sed 's/^name/\#VID/g' \
    | gzip -c \
    > "${prefix}.frequencies.allPops.txt.gz"
  >>>

  output {
    File freq_data = "${prefix}.frequencies.txt.gz"
    File freq_data_allPops = "${prefix}.frequencies.allPops.txt.gz"
    # File freqs_preclean = "${prefix}.frequencies.preclean.txt.gz"
  }

  runtime {
    docker : "talkowski/sv-pipeline@sha256:aef8156983cec6ac6a91fa6461b197a63835e5487fc9523ec857f947cfac660e"
    disks: "local-disk 50 HDD"
    memory: "4 GB"
    preemptible: 1
    maxRetries: 1
  }
}


# Combine frequency data across batches
task merge_freq_tables {
  Array[File] tables
  File batches_list
  String prefix

  command <<<
    set -euo pipefail
    #Get list of batch IDs and batch table paths
    while read batch; do
      echo "$batch"
      find / -name "$batch.frequencies*txt.gz"
    done < ${batches_list} | paste - - \
    > input.list
    #Make sure all input files have the same number of lines
    while read batch file; do
      zcat "$file" | wc -l
    done < input.list > nlines.list
    nlines=$( cat nlines.list | sort | uniq | wc -l )
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
    docker : "talkowski/sv-pipeline@sha256:aef8156983cec6ac6a91fa6461b197a63835e5487fc9523ec857f947cfac660e"
    disks: "local-disk 100 HDD"
    memory: "16 GB"
    preemptible: 1
    maxRetries: 1
  }
}


# Compare 
task compare_freqs_prePost_minGQ {
  File AF_PCRPLUS_preMinGQ
  File AF_PCRMINUS_preMinGQ
  File AF_postMinGQ_table
  String prefix

  command <<<
    set -euo pipefail
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/compare_freqs_pre_post_minGQ.R \
      ${AF_PCRPLUS_preMinGQ} \
      ${AF_PCRMINUS_preMinGQ} \
      ${AF_postMinGQ_table} \
      ./ \
      "${prefix}."
  >>>

  output {
    File PCRPLUS_fails = "${prefix}.PCRPLUS_minGQ_AF_prePost_fails.VIDs.list"
    File PCRMINUS_fails = "${prefix}.PCRMINUS_minGQ_AF_prePost_fails.VIDs.list"
    File minGQ_prePost_comparison_data = "${prefix}.minGQ_AF_prePost_comparison.data.txt.gz"
    File minGQ_prePost_comparison_plot = "${prefix}.minGQ_AF_prePost_comparison.plot.png"
  }

  runtime {
    docker : "talkowski/sv-pipeline@sha256:6552363ba2a13147940084658eff35171242c3e71442bc6936a2b5e53b3501ce"
    disks: "local-disk 30 HDD"
    memory: "8 GB"
    preemptible: 1
    maxRetries: 1
  }
}


# Calculate & plot cross-batch correlation coefficient matrixes
task make_correlation_matrices {
  File freq_table
  String pop
  File batches_list
  String prefix

  command <<<
    set -euo pipefail
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
    docker : "talkowski/sv-pipeline@sha256:aef8156983cec6ac6a91fa6461b197a63835e5487fc9523ec857f947cfac660e"
    disks: "local-disk 50 HDD"
    memory: "8 GB"
    preemptible: 1
    maxRetries: 1
  }
}


# Generate list of all pairs of batches to be compared
task make_batch_pairs_list {
  File batches_list
  String prefix

  command <<<
    set -euo pipefail
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/make_batch_pairs_list.R \
      ${batches_list} \
      "${prefix}.nonredundant_batch_pairs.txt"
  >>>

  output {
    File batch_pairs_list = "${prefix}.nonredundant_batch_pairs.txt"
  }

  runtime {
    docker : "talkowski/sv-pipeline@sha256:aef8156983cec6ac6a91fa6461b197a63835e5487fc9523ec857f947cfac660e"
    preemptible: 1
    maxRetries: 1
  }
}


# Merge lists of batch effect checks and count total number of times each variant failed
task merge_variant_failure_lists {
  Array[File] fail_variant_lists
  String prefix

  command <<<
    set -euo pipefail
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
    > "${prefix}.PCRPLUS_to_PCRPLUS.failures.txt" || true
    #Get master list of PCR- to PCR- failures
    fgrep -v "PCRPLUS" fail_variant_lists.paths.txt \
    | fgrep "PCRMINUS" \
    | xargs -I {} cat {} \
    | sort -Vk1,1 | uniq -c \
    | awk -v OFS="\t" '{ print $2, $1 }' \
    > "${prefix}.PCRMINUS_to_PCRMINUS.failures.txt" || true
    #Get master list of PCR+ to PCR- failures
    fgrep "PCRPLUS" fail_variant_lists.paths.txt \
    | fgrep "PCRMINUS" \
    | xargs -I {} cat {} \
    | sort -Vk1,1 | uniq -c \
    | awk -v OFS="\t" '{ print $2, $1 }' \
    > "${prefix}.PCRPLUS_to_PCRMINUS.failures.txt" || true
    #Get master list of all possible failures
    cat fail_variant_lists.paths.txt \
    | xargs -I {} cat {} \
    | sort -Vk1,1 | uniq -c \
    | awk -v OFS="\t" '{ print $2, $1 }' \
    > "${prefix}.all.failures.txt" || true
  >>>

  output {
    File fails_per_variant_PCRPLUS_to_PCRPLUS = "${prefix}.PCRPLUS_to_PCRPLUS.failures.txt"
    File fails_per_variant_PCRMINUS_to_PCRMINUS = "${prefix}.PCRMINUS_to_PCRMINUS.failures.txt"
    File fails_per_variant_PCRPLUS_to_PCRMINUS = "${prefix}.PCRPLUS_to_PCRMINUS.failures.txt"
    File fails_per_variant_all = "${prefix}.all.failures.txt"
  }

  runtime {
    docker : "talkowski/sv-pipeline@sha256:aef8156983cec6ac6a91fa6461b197a63835e5487fc9523ec857f947cfac660e"
    preemptible: 1
    maxRetries: 1
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
    set -euo pipefail
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
    docker : "talkowski/sv-pipeline@sha256:aef8156983cec6ac6a91fa6461b197a63835e5487fc9523ec857f947cfac660e"
    preemptible: 1
    maxRetries: 1
    memory: "8 GB"
  }
}


# Apply batch effect labels to VCF
task apply_batch_effect_labels {
  File vcf
  File vcf_idx
  String contig
  File reclassification_table
  File PCRPLUS_samples_list
  File minGQ_prePost_PCRPLUS_fails
  File minGQ_prePost_PCRMINUS_fails
  String prefix

  command <<<
    set -euo pipefail
    tabix -h ${vcf} ${contig} \
    | /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/label_batch_effects.py \
        --unstable-af-pcrplus ${minGQ_prePost_PCRPLUS_fails} \
        --unstable-af-pcrminus ${minGQ_prePost_PCRMINUS_fails} \
        stdin \
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
    docker : "talkowski/sv-pipeline@sha256:4a428d92c780b70f0fc57ca2024c2edc19aefd0c6e08dc9f01515bd4f7804a4d"
    preemptible: 1
    maxRetries: 1
    disks: "local-disk 50 HDD"
    memory: "4 GB"
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
    docker: "talkowski/sv-pipeline@sha256:4a428d92c780b70f0fc57ca2024c2edc19aefd0c6e08dc9f01515bd4f7804a4d"
    preemptible: 0
    maxRetries: 1
    memory: "16 GB"
    disks: "local-disk 250 SSD"
  }
}


