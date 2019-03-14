import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:gather_attribute_paths_multiSampleSet/versions/7/plain-WDL/descriptor" as getAttribute

# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# Workflow to identify & filter outliers from VCFs after module 03 (random forest)
workflow filter_outlier_samples {
	File delly_vcf
  File manta_vcf
  File melt_vcf
  File depth_vcf
  Int N_IQR_cutoff
  String batch
  Array[String] samples

  # Write original list of unfiltered samples
  call write_samples_list {
    input:
      samples=samples
  }

  # Get delly outliers
  call identify_outliers as get_delly_outliers {
    input:
      vcf=delly_vcf,
      N_IQR_cutoff=N_IQR_cutoff,
      outfile="delly_outliers.txt"
  }

  # Get manta outliers
  call identify_outliers as get_manta_outliers {
    input:
      vcf=manta_vcf,
      N_IQR_cutoff=N_IQR_cutoff,
      outfile="manta_outliers.txt"
  }

  # Get melt outliers
  call identify_outliers as get_melt_outliers {
    input:
      vcf=melt_vcf,
      N_IQR_cutoff=N_IQR_cutoff,
      outfile="melt_outliers.txt"
  }

  # Get depth outliers
  call identify_outliers as get_depth_outliers {
    input:
      vcf=depth_vcf,
      N_IQR_cutoff=N_IQR_cutoff,
      outfile="depth_outliers.txt"
  }

  # Merge list of outliers
  call cat_outliers {
    input:
      delly_outliers=get_delly_outliers.outliers_list,
      manta_outliers=get_manta_outliers.outliers_list,
      melt_outliers=get_melt_outliers.outliers_list,
      depth_outliers=get_depth_outliers.outliers_list,
      batch=batch
  }

  # Exclude outliers from delly
  call exclude_outliers as exclude_outliers_delly {
    input:
      vcf=delly_vcf,
      outliers_list=cat_outliers.outliers_list,
      outfile="${batch}.delly.outliers_removed.vcf.gz"
  }

  # Exclude outliers from manta
  call exclude_outliers as exclude_outliers_manta {
    input:
      vcf=manta_vcf,
      outliers_list=cat_outliers.outliers_list,
      outfile="${batch}.manta.outliers_removed.vcf.gz"
  }

  # Exclude outliers from melt
  call exclude_outliers as exclude_outliers_melt {
    input:
      vcf=melt_vcf,
      outliers_list=cat_outliers.outliers_list,
      outfile="${batch}.melt.outliers_removed.vcf.gz"
  }

  # Exclude outliers from depth
  call exclude_outliers as exclude_outliers_depth {
    input:
      vcf=depth_vcf,
      outliers_list=cat_outliers.outliers_list,
      outfile="${batch}.depth.outliers_removed.vcf.gz"
  }

  # Write new list of samples without outliers
  call filter_sample_list {
    input:
      original_samples_list=write_samples_list.samples_list,
      outlier_samples=cat_outliers.outliers_list,
      batch=batch
  }

  # Final outputs
  output {
    File delly_vcf_noOutliers = exclude_outliers_delly.vcf_no_outliers
    File manta_vcf_noOutliers = exclude_outliers_manta.vcf_no_outliers
    File melt_vcf_noOutliers = exclude_outliers_melt.vcf_no_outliers
    File depth_vcf_noOutliers = exclude_outliers_depth.vcf_no_outliers
    File filtered_batch_samples_list = filter_sample_list.filtered_samples_list
    File outlier_samples_excluded = cat_outliers.outliers_list
  }

}


# Write original list of samples
task write_samples_list {
  Array[String] samples

  command <<<
    cat ${write_tsv(samples)} > samples.list
  >>>

  output {
    File samples_list = "samples.list"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:7e7e6163d6ac0fc5781eb99ee5a7eec4db37506f48d00f5063b96123f9ca5024"
    preemptible: 3
  }
}


# Identify the list of outlier samples from a single VCF
task identify_outliers {
  File vcf
  Int N_IQR_cutoff
  String outfile

  command <<<
    # Count sv per class per sample
    svtk count-svtypes ${vcf} svcounts.txt

    # Return list of samples exceeding cutoff for at least one sv class
    /opt/sv-pipeline/03_variant_filtering/scripts/get_outliers_from_svcounts.helper.R \
      svcounts.txt \
      ${N_IQR_cutoff} \
      ${outfile}
  >>>

  output {
    File outliers_list = "${outfile}"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:7e7e6163d6ac0fc5781eb99ee5a7eec4db37506f48d00f5063b96123f9ca5024"
    preemptible: 3
    disks: "local-disk 100 HDD"
    memory: "8 GB"
  }
}


# Merge outlier sample lists across algorithms
task cat_outliers {
  File delly_outliers
  File manta_outliers
  File melt_outliers
  File depth_outliers
  String batch

  command <<<
    cat ${delly_outliers} ${manta_outliers} ${melt_outliers} ${depth_outliers} \
      | sort | uniq > ${batch}.post03_outliers.samples.list
  >>>

  output {
    File outliers_list = "${batch}.post03_outliers.samples.list"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:7e7e6163d6ac0fc5781eb99ee5a7eec4db37506f48d00f5063b96123f9ca5024"
    preemptible: 3
  }
}


# Exclude outliers from VCF
task exclude_outliers {
  File vcf
  File outliers_list
  String outfile

  command <<<
    zcat ${vcf} | fgrep "#" | fgrep -v "##" | \
      sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' | \
      fgrep -wf ${outliers_list} | cut -f2 > \
      indexes_to_exclude.txt
    if [ $( cat indexes_to_exclude.txt | wc -l ) -gt 0 ]; then
      zcat ${vcf} | \
      cut --complement -f$( cat indexes_to_exclude.txt | paste -s -d, ) | \
      vcftools --mac 1 --vcf - --recode --recode-INFO-all --stdout | \
      bgzip -c > ${outfile}
    else
      cp ${vcf} ${outfile}
    fi
  >>>

  output {
    File vcf_no_outliers = "${outfile}"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:7e7e6163d6ac0fc5781eb99ee5a7eec4db37506f48d00f5063b96123f9ca5024"
    preemptible: 3
    disks: "local-disk 100 HDD"
  }
}


# Write new list of samples per batch after outlier filtering
task filter_sample_list {
  File original_samples_list
  File outlier_samples
  String batch

  command <<<
    fgrep -wvf ${outlier_samples} ${original_samples_list} > \
    ${batch}.post03_outliers_excluded.samples.list
  >>>

  output {
    File filtered_samples_list = "${batch}.post03_outliers_excluded.samples.list"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:7e7e6163d6ac0fc5781eb99ee5a7eec4db37506f48d00f5063b96123f9ca5024"
    preemptible: 3
  }
}
