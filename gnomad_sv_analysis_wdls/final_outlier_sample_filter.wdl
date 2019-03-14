import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:gather_attribute_paths_multiSampleSet/versions/7/plain-WDL/descriptor" as getAttribute

# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License

# This is an analysis WDL to identify & filter outliers from VCFs 
# after minGQ filtering at the end of the Talkowski SV pipeline

# Treats PCR+ and PCR- samples separately

workflow filter_outlier_samples {
  File vcf
  File vcf_idx
  File PCRPLUS_samples_list
  Int N_IQR_cutoff_PCRPLUS
  Int N_IQR_cutoff_PCRMINUS
  String prefix
  File autosomes_list

  Array[Array[String]] contigs=read_tsv(autosomes_list)

  # Write original list of unfiltered samples and split by PCR status
  call write_samples_list {
    input:
      vcf=vcf,
      vcf_idx=vcf_idx,
      PCRPLUS_samples_list=PCRPLUS_samples_list,
      prefix=prefix
  }

  # Get count of biallelic autosomal variants per sample
  scatter ( contig in contigs ) {
    call count_svtypes {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        prefix=prefix,
        contig=contig[0]
    }
  }
  call combine_counts {
    input:
      svcounts=count_svtypes.sv_counts,
      prefix=prefix
  }

  # Get outliers
  call identify_outliers as identify_PCRPLUS_outliers {
    input:
      svcounts=combine_counts.summed_svcounts,
      N_IQR_cutoff=N_IQR_cutoff_PCRPLUS,
      samples_list=write_samples_list.plus_samples_list,
      prefix="${prefix}.PCRPLUS"
  }
  call identify_outliers as identify_PCRMINUS_outliers {
    input:
      svcounts=combine_counts.summed_svcounts,
      N_IQR_cutoff=N_IQR_cutoff_PCRMINUS,
      samples_list=write_samples_list.minus_samples_list,
      prefix="${prefix}.PCRMINUS"
  }

  # Exclude outliers from vcf
  call exclude_outliers {
    input:
      vcf=vcf,
      vcf_idx=vcf_idx,
      plus_outliers_list=identify_PCRPLUS_outliers.outliers_list,
      minus_outliers_list=identify_PCRMINUS_outliers.outliers_list,
      outfile="${prefix}.outliers_removed.vcf.gz",
      prefix=prefix
  }

  # Write new list of samples without outliers
  call filter_sample_list {
    input:
      original_samples_list=write_samples_list.samples_list,
      outlier_samples=exclude_outliers.merged_outliers_list,
      prefix=prefix
  }

  # Final outputs
  output {
    File vcf_noOutliers = exclude_outliers.vcf_no_outliers
    File vcf_noOutliers_idx = exclude_outliers.vcf_no_outliers_idx
    File nooutliers_samples_list = filter_sample_list.filtered_samples_list
    File excluded_samples_list = exclude_outliers.merged_outliers_list
    File svcounts_per_sample_data = combine_counts.summed_svcounts
    File svcounts_per_sample_plots_PCRPLUS = identify_PCRPLUS_outliers.svcount_distrib_plots
    File svcounts_per_sample_plots_PCRMINUS = identify_PCRMINUS_outliers.svcount_distrib_plots
  }
}


# Write original list of samples
task write_samples_list {
  File vcf
  File vcf_idx
  File PCRPLUS_samples_list
  String prefix

  command <<<
    tabix -H ${vcf} | fgrep -v "##" \
    | cut -f10- | sed 's/\t/\n/g' > "${prefix}.samples.list"
    fgrep -wf ${PCRPLUS_samples_list} "${prefix}.samples.list" \
    > "${prefix}.PCRPLUS.samples.list"
    fgrep -wvf ${PCRPLUS_samples_list} "${prefix}.samples.list" \
    > "${prefix}.PCRMINUS.samples.list"
  >>>

  output {
    File samples_list = "${prefix}.samples.list"
    File plus_samples_list = "${prefix}.PCRPLUS.samples.list"
    File minus_samples_list = "${prefix}.PCRMINUS.samples.list"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:ef7584fc2cd354567d98b7f0d8ba4c83ac79f73a8c337ebaf765f8ff008c274c"
    preemptible: 1
    disks: "local-disk 50 HDD"
  }
}


# Count biallelic SV per sample for a single chromosome
task count_svtypes {
  File vcf
  File vcf_idx
  String prefix
  String contig

  command <<<
    tabix --print-header "${vcf}" "${contig}" \
    | fgrep -v "MULTIALLELIC" \
    | fgrep -v "PESR_GT_OVERDISPERSION" \
    | svtk count-svtypes --no-header stdin \
    | awk -v OFS="\t" -v chr="${contig}" '{ print $0, chr }' \
    > "${prefix}.${contig}.svcounts.txt"
  >>>

  output {
    File sv_counts = "${prefix}.${contig}.svcounts.txt"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:d9dad81e3e62423b4488b598268c9d6a657c882fef8af5f0bc46a1440afeb0c9"
    preemptible: 1
    disks: "local-disk 50 HDD"
  }
}


# Combine SV count files across chromosomes
task combine_counts {
  Array[File] svcounts
  String prefix

  command <<<
    while read file; do
      cat "$file"
    done < ${write_lines(svcounts)} \
    > merged_svcounts.txt
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/sum_svcounts_perSample.R \
      merged_svcounts.txt \
      "${prefix}.summed_svcounts_per_sample.txt"
  >>>


  output {
    File summed_svcounts = "${prefix}.summed_svcounts_per_sample.txt"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:58bb6f91431f2dcf17f5f796e9a7bb50d27c58fa259471837489a0d4db7ae6aa"
    preemptible: 1
    disks: "local-disk 30 HDD"
    memory: "4 GB"
  }
}


# Identify the list of outlier samples & generate distribution plots
task identify_outliers {
  File svcounts
  Int N_IQR_cutoff
  File samples_list
  String prefix

  command <<<
    # Subset input data to specified samples
    sed -n '1p' ${svcounts} > filtered_counts.input.txt
    sed '1d' ${svcounts} | fgrep -wf ${samples_list} >> filtered_counts.input.txt
    # Return list of samples exceeding cutoff for at least one sv class
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/determine_svcount_outliers.R  \
      -p "${prefix}" \
      -I "${N_IQR_cutoff}" \
      filtered_counts.input.txt \
      "${prefix}_svcount_outlier_plots/"
    cat "${prefix}_svcount_outlier_plots/${prefix}.SV_count_outlier_samples.txt" \
      | fgrep -v "#" | cut -f1 | sort | uniq \
      > "${prefix}.SV_count_outliers.samples.list"
    tar -cvzf "${prefix}_svcount_outlier_plots.tar.gz" "${prefix}_svcount_outlier_plots/"
  >>>

  output {
    File outliers_list = "${prefix}.SV_count_outliers.samples.list"
    File svcount_distrib_plots = "${prefix}_svcount_outlier_plots.tar.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:3acc88238e71346ea1e8fce673a3b613e54b8fa3f2748cee90f3f8fc7c39fc65"
    preemptible: 1
    disks: "local-disk 20 HDD"
    memory: "4 GB"
  }
}


# Exclude outliers from VCF
task exclude_outliers {
  File vcf
  File vcf_idx
  File plus_outliers_list
  File minus_outliers_list
  String outfile
  String prefix

  command <<<
    cat ${plus_outliers_list} ${minus_outliers_list} \
      | sort -Vk1,1 | uniq \
      > "${prefix}.SV_count_outliers.samples.list"
    tabix -H ${vcf} | fgrep -v "##" | \
      sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' | \
      fgrep -wf "${prefix}.SV_count_outliers.samples.list" | cut -f2 > \
      indexes_to_exclude.txt
    if [ $( cat indexes_to_exclude.txt | wc -l ) -gt 0 ]; then
      zcat ${vcf} | \
      cut --complement -f$( cat indexes_to_exclude.txt | paste -s -d, ) | \
      bgzip -c \
      > "${prefix}.subsetted_preEmptyRemoval.vcf.gz"
      /opt/sv-pipeline/scripts/drop_empty_records.py \
        "${prefix}.subsetted_preEmptyRemoval.vcf.gz" \
        stdout | \
      bgzip -c > ${outfile}
    else
      cp ${vcf} ${outfile}
    fi
    tabix -p vcf -f "${outfile}"
  >>>

  output {
    File merged_outliers_list = "${prefix}.SV_count_outliers.samples.list"
    File vcf_no_outliers = "${outfile}"
    File vcf_no_outliers_idx = "${outfile}.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:d3a7eabbd8a6c79cfa0426dac8fbe6527a445d0a104e270eecbbd4ad66cb74fe"
    preemptible: 1
    disks: "local-disk 100 HDD"
  }
}


# Write new list of samples per prefix after outlier filtering
task filter_sample_list {
  File original_samples_list
  File outlier_samples
  String prefix

  command <<<
    fgrep -wvf ${outlier_samples} ${original_samples_list} > \
    ${prefix}.outliers_excluded.samples.list
  >>>

  output {
    File filtered_samples_list = "${prefix}.outliers_excluded.samples.list"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:3acc88238e71346ea1e8fce673a3b613e54b8fa3f2748cee90f3f8fc7c39fc65"
    preemptible: 1
  }
}
