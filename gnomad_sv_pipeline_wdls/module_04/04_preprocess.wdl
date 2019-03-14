import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:gather_attribute_paths_multiSampleSet/versions/7/plain-WDL/descriptor" as getAttribute

# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# Workflow to preprocess all files needed for per-batch genotyping in module 04a
workflow preprocess_04a_files {

  File sample_set_list
  File svcActKeyJson
  String workspaceProject
  String workspaceName

  # Get cohort_filtered_pesr_vcf_list
  call getAttribute.gather_attribute_paths_multiSampleSet as get_pesr_vcf_list {
    input:
      sample_set_list=sample_set_list,
      Attribute="filtered_pesr_vcf",
      svcActKeyJson=svcActKeyJson,
      workspaceProject=workspaceProject,
      workspaceName=workspaceName
  }

  # Get cohort_filtered_depth_vcf_list
  call getAttribute.gather_attribute_paths_multiSampleSet as get_depth_vcf_list {
    input:
      sample_set_list=sample_set_list,
      Attribute="filtered_depth_vcf",
      svcActKeyJson=svcActKeyJson,
      workspaceProject=workspaceProject,
      workspaceName=workspaceName
  }

  # Outputs
  output {
    File filtered_pesr_vcf_list = get_pesr_vcf_list.attribute_list
    File filtered_depth_vcf_list = get_depth_vcf_list.attribute_list
  }
}