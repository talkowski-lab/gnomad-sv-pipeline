import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:gather_attribute_paths_multiSampleSet/versions/7/plain-WDL/descriptor" as getAttribute

# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# Workflow to preprocess all files needed for batch integration in module 04b
workflow preprocess_04b_files {

  File sample_set_list
  File svcActKeyJson
  String workspaceProject
  String workspaceName
  Float min_sr_background_fail_batches

  # Get cohort_genotyped_pesr_vcf_list
  call getAttribute.gather_attribute_paths_multiSampleSet as get_pesr_vcf_list {
    input:
      sample_set_list=sample_set_list,
      Attribute="genotyped_pesr_vcf",
      svcActKeyJson=svcActKeyJson,
      workspaceProject=workspaceProject,
      workspaceName=workspaceName
  }

  # Get cohort_genotyped_pesr_vcf_idx_list
  call getAttribute.gather_attribute_paths_multiSampleSet as get_pesr_vcf_idx_list {
    input:
      sample_set_list=sample_set_list,
      Attribute="genotyped_pesr_vcf_idx",
      svcActKeyJson=svcActKeyJson,
      workspaceProject=workspaceProject,
      workspaceName=workspaceName
  }

  # Get cohort_genotyped_depth_vcf_list
  call getAttribute.gather_attribute_paths_multiSampleSet as get_depth_vcf_list {
    input:
      sample_set_list=sample_set_list,
      Attribute="genotyped_depth_vcf",
      svcActKeyJson=svcActKeyJson,
      workspaceProject=workspaceProject,
      workspaceName=workspaceName
  }

  # Get cohort_genotyped_depth_vcf_idx_list
  call getAttribute.gather_attribute_paths_multiSampleSet as get_depth_vcf_idx_list {
    input:
      sample_set_list=sample_set_list,
      Attribute="genotyped_depth_vcf_idx",
      svcActKeyJson=svcActKeyJson,
      workspaceProject=workspaceProject,
      workspaceName=workspaceName
  }

  # Get cohort_discfile_list
  call getAttribute.gather_attribute_paths_multiSampleSet as get_discfile_list {
    input:
      sample_set_list=sample_set_list,
      Attribute="PE_file",
      svcActKeyJson=svcActKeyJson,
      workspaceProject=workspaceProject,
      workspaceName=workspaceName
  }

  # Get cohort_discfile_idx_list
  call getAttribute.gather_attribute_paths_multiSampleSet as get_discfile_idx_list {
    input:
      sample_set_list=sample_set_list,
      Attribute="PE_file_idx",
      svcActKeyJson=svcActKeyJson,
      workspaceProject=workspaceProject,
      workspaceName=workspaceName
  }

  # Get cohort_bincov_list
  call getAttribute.gather_attribute_paths_multiSampleSet as get_bincov_list {
    input:
      sample_set_list=sample_set_list,
      Attribute="bincov_file",
      svcActKeyJson=svcActKeyJson,
      workspaceProject=workspaceProject,
      workspaceName=workspaceName
  }

  # Get cohort_bincov_idx_list
  call getAttribute.gather_attribute_paths_multiSampleSet as get_bincov_idx_list {
    input:
      sample_set_list=sample_set_list,
      Attribute="bincov_file_idx",
      svcActKeyJson=svcActKeyJson,
      workspaceProject=workspaceProject,
      workspaceName=workspaceName
  }

  # Get cohort_sr_genotyping_bothside_pass_list
  call getAttribute.gather_attribute_paths_multiSampleSet as get_bothside_pass_list {
    input:
      sample_set_list=sample_set_list,
      Attribute="genotyping_sr_bothside_pass",
      svcActKeyJson=svcActKeyJson,
      workspaceProject=workspaceProject,
      workspaceName=workspaceName
  }
  call clean_bothside_pass {
    input:
      gs_list=get_bothside_pass_list.attribute_list
  }

  # Get cohort_sr_genotyping_background_fail_list
  call getAttribute.gather_attribute_paths_multiSampleSet as get_background_fail_list {
    input:
      sample_set_list=sample_set_list,
      Attribute="genotyping_sr_background_fail",
      svcActKeyJson=svcActKeyJson,
      workspaceProject=workspaceProject,
      workspaceName=workspaceName
  }
  call clean_background_fail {
    input:
      gs_list=get_background_fail_list.attribute_list,
      min_sr_background_fail_batches=min_sr_background_fail_batches
  }

  # Create cohort famfile
  call getAttribute.gather_attribute_paths_multiSampleSet as get_famfile_list {
    input:
      sample_set_list=sample_set_list,
      Attribute="famfile_after_outlier_exclusion",
      svcActKeyJson=svcActKeyJson,
      workspaceProject=workspaceProject,
      workspaceName=workspaceName
  }
  Array[File] famfile_list_array = read_lines(get_famfile_list.attribute_list)
  call cat_files as clean_famfile_list {
    input:
      files=famfile_list_array,
      outfile="merged_famfile.fam"
  }

  # Get cohort random forest cutoffs
  call getAttribute.gather_attribute_paths_multiSampleSet as get_RF_cutoffs_list {
    input:
      sample_set_list=sample_set_list,
      Attribute="RF_cutoffs",
      svcActKeyJson=svcActKeyJson,
      workspaceProject=workspaceProject,
      workspaceName=workspaceName
  }

  # Get cohort depth_genotying_RD_depth_sepcutoff
  call getAttribute.gather_attribute_paths_multiSampleSet as get_depth_gt_RD_depth_sep_list {
    input:
      sample_set_list=sample_set_list,
      Attribute="depth_genotying_RD_depth_sepcutoff",
      svcActKeyJson=svcActKeyJson,
      workspaceProject=workspaceProject,
      workspaceName=workspaceName
  }

  # Get cohort bincov median file list
  call getAttribute.gather_attribute_paths_multiSampleSet as get_medianfile_list {
    input:
      sample_set_list=sample_set_list,
      Attribute="medianfile",
      svcActKeyJson=svcActKeyJson,
      workspaceProject=workspaceProject,
      workspaceName=workspaceName
  }

  # Get sample lists after outlier exclusion
  call getAttribute.gather_attribute_paths_multiSampleSet as get_samples_list {
    input:
      sample_set_list=sample_set_list,
      Attribute="sample_list_after_outlier_exclusion",
      svcActKeyJson=svcActKeyJson,
      workspaceProject=workspaceProject,
      workspaceName=workspaceName
  }

  # Outputs
  output {
    File genotyped_pesr_vcf_list = get_pesr_vcf_list.attribute_list
    File genotyped_pesr_vcf_idx_list = get_pesr_vcf_idx_list.attribute_list
    File genotyped_depth_vcf_list = get_depth_vcf_list.attribute_list
    File genotyped_depth_vcf_idx_list = get_depth_vcf_idx_list.attribute_list
    File discfile_list = get_discfile_list.attribute_list
    File discfile_idx_list = get_discfile_idx_list.attribute_list
    File bincov_list = get_bincov_list.attribute_list
    File bincov_idx_list = get_bincov_idx_list.attribute_list
    File sr_bothside_pass = clean_bothside_pass.bothside_pass
    File sr_background_fail = clean_background_fail.background_fail
    File famfile_list = get_famfile_list.attribute_list
    File famfile = clean_famfile_list.merged_file
    File RF_cutoffs = get_RF_cutoffs_list.attribute_list
    File depth_gt_RD_sep = get_depth_gt_RD_depth_sep_list.attribute_list
    File medianfile_list = get_medianfile_list.attribute_list
    File samples_list = get_samples_list.attribute_list
  }
}


# Clean sr bothside pass list
task clean_bothside_pass {
  File gs_list

  command <<<
    mkdir lists/
    while read gspath; do
      gsutil cp $gspath lists/
    done < ${gs_list}
    cat lists/* | sort | uniq -c | \
    awk -v nbatch=$( find lists/ -name "*.txt" | wc -l ) -v OFS="\t" \
    '{ print $1/nbatch, $2 }' > cohort_sr_genotyping_bothside_pass_list.txt
  >>>

  output {
    File bothside_pass = "cohort_sr_genotyping_bothside_pass_list.txt"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:26fbe9ea59a4043466277062905532178767b48a214a0bc5782a5a01f242d158"
    preemptible: 3
  }
}


# Clean sr background fail list
task clean_background_fail {
  File gs_list
  Float min_sr_background_fail_batches

  command <<<
    mkdir lists/
    while read gspath; do
      gsutil cp $gspath lists/
    done < ${gs_list}
    cat lists/* | sort | uniq -c | \
    awk -v nbatch=$( find lists/ -name "*.txt" | wc -l ) -v min=${min_sr_background_fail_batches} -v OFS="\t" \
    '{ if ($1/nbatch>=min) print $2 }' > cohort_sr_genotyping_background_fail_list.txt
  >>>

  output {
    File background_fail = "cohort_sr_genotyping_background_fail_list.txt"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:26fbe9ea59a4043466277062905532178767b48a214a0bc5782a5a01f242d158"
    preemptible: 3
  }
}


# Cat files from array
task cat_files {
  Array[File] files
  String outfile

  command <<<
    cat ${sep=' ' files} > ${outfile}
  >>>

  output {
    File merged_file = "${outfile}"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:26fbe9ea59a4043466277062905532178767b48a214a0bc5782a5a01f242d158"
    preemptible: 3
  }
}
