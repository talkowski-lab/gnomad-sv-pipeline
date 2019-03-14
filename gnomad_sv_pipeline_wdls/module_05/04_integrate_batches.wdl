import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04_vcfcluster_by_chrom/versions/57/plain-WDL/descriptor" as vcfcluster_by_chrom
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04_pesr_depth_overlap/versions/26/plain-WDL/descriptor" as pesr_depth_overlap
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04_bp_overlap_filter_by_chrom/versions/2/plain-WDL/descriptor" as bp_overlap_filter_by_chrom
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04_resolve_complex_SV/versions/101/plain-WDL/descriptor" as resolve_complex_sv
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04_integrate_resolved_vcfs/versions/2/plain-WDL/descriptor" as integrate_resolved_vcfs
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04_genotype_CPX_CNVs/versions/30/plain-WDL/descriptor" as cpx_gt
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:master_SV_VCF_QC/versions/49/plain-WDL/descriptor" as vcf_qc

workflow batch_integration {
  File pesr_vcf_list
  File pesr_vcf_idx_list
  File depth_vcf_list
  File depth_vcf_idx_list
  File blacklist
  File blacklist_idx
  Array[String] samples
  File contigs
  Int max_shards_per_chrom
  Int min_variants_per_shard
  File cytobands
  File cytobands_idx
  File discfile_list
  File discfile_idx_list
  File bincov_list
  File bincov_idx_list
  File mei_bed
  File pe_blacklist
  File pe_blacklist_idx
  String prefix
  File svc_acct_key
  File famfile
  File trios_famfile
  String ref_build
  File Sanders_2015_tarball
  File Collins_2017_tarball
  File Werling_2018_tarball
  File sr_background_fail
  File sr_bothend_pass
  File rf_cutoffs
  File batches_list
  File depth_gt_rd_sep_list
  File medianfile_list
  File famfile_list
  File sampleslist_list

  # Array[String] pesr_vcfs = read_lines(pesr_vcf_list)
  # Array[File] pesr_vcf_idxs = read_lines(pesr_vcf_idx_list)
  # Array[String] depth_vcfs = read_lines(depth_vcf_list)
  # Array[File] depth_vcf_idxs = read_lines(depth_vcf_idx_list)

  # Cluster PESR VCFs across batches
  call vcfcluster_by_chrom.vcfcluster_by_chrom as cluster_pesr {
    input:
      vcf_list=pesr_vcf_list,
      vcf_idx_list=pesr_vcf_idx_list,
      batches_list=batches_list,
      prefix="AllBatches_pesr",
      dist=300,
      frac=0.1,
      sample_overlap=0.5,
      do_blacklist="YES",
      blacklist=blacklist,
      blacklist_idx=blacklist_idx,
      svsize=50,
      svtypes=["DEL","DUP","INV","BND","INS"],
      contiglist=contigs,
      max_shards_per_chrom_svtype=100,
      min_variants_per_shard_per_chrom_svtype=100,
      svc_acct_key=svc_acct_key
  }

  # Update SR background fail & bothside pass files
  call update_sr_list as update_background_fail {
    input:
      vcf=cluster_pesr.clustered_vcf,
      original_list=sr_background_fail,
      outfile="sr_background_fail.updated.txt"
  }

  call update_sr_list as update_bothside_pass {
    input:
      vcf=cluster_pesr.clustered_vcf,
      original_list=sr_bothend_pass,
      outfile="sr_bothside_pass.updated.txt"
  }

  # Run QC on clustered PESR VCFs 
  call vcf_qc.master_vcf_qc as pesr_vcf_qc {
    input:
      vcf=cluster_pesr.clustered_vcf,
      famfile=trios_famfile,
      ref_build=ref_build,
      prefix="${prefix}_PESR_VCF",
      sv_per_shard=10000,
      samples_per_shard=100,
      Sanders_2015_tarball=Sanders_2015_tarball,
      Collins_2017_tarball=Collins_2017_tarball,
      Werling_2018_tarball=Werling_2018_tarball
  }

  # Cluster RD VCFs across batches
  call vcfcluster_by_chrom.vcfcluster_by_chrom as cluster_depth {
    input:
      vcf_list=depth_vcf_list,
      vcf_idx_list=depth_vcf_idx_list,
      batches_list=batches_list,
      prefix="AllBatches_depth",
      dist=500000,
      frac=0.5,
      sample_overlap=0.5,
      do_blacklist="NO",
      blacklist=blacklist,
      blacklist_idx=blacklist_idx,
      svsize=5000,
      svtypes=["DEL","DUP"],
      contiglist=contigs,
      max_shards_per_chrom_svtype=100,
      min_variants_per_shard_per_chrom_svtype=100,
      svc_acct_key=svc_acct_key
  }

  # Run QC on clustered RD VCFs 
  call vcf_qc.master_vcf_qc as rd_vcf_qc {
    input:
      vcf=cluster_depth.clustered_vcf,
      famfile=trios_famfile,
      ref_build=ref_build,
      prefix="${prefix}_RD_VCF",
      sv_per_shard=10000,
      samples_per_shard=100,
      Sanders_2015_tarball=Sanders_2015_tarball,
      Collins_2017_tarball=Collins_2017_tarball,
      Werling_2018_tarball=Werling_2018_tarball
  }

  # Scatter over chromosomes and merge clustered PESR & RD VCFs
  call pesr_depth_overlap.pesr_depth_overlap {
    input:
      pesr_vcf=cluster_pesr.clustered_vcf,
      pesr_vcf_idx=cluster_pesr.clustered_vcf_idx,
      depth_vcf=cluster_depth.clustered_vcf,
      depth_vcf_idx=cluster_depth.clustered_vcf_idx,
      samples=samples,
      contigs=contigs,
      svc_acct_key=svc_acct_key
  }

  # Update SR background fail & bothside pass files
  call update_sr_list as update_background_fail_second {
    input:
      vcf=pesr_depth_overlap.merged_vcf,
      original_list=update_background_fail.updated_list,
      outfile="sr_background_fail.updated2.txt"
  }

  call update_sr_list as update_bothside_pass_second {
    input:
      vcf=pesr_depth_overlap.merged_vcf,
      original_list=update_bothside_pass.updated_list,
      outfile="sr_bothside_pass.updated2.txt"
  }
  
  # QC merged PESR/RD VCF
  call vcf_qc.master_vcf_qc as pesr_rd_vcf_qc {
    input:
      vcf=pesr_depth_overlap.merged_vcf,
      famfile=trios_famfile,
      ref_build=ref_build,
      prefix="${prefix}_PESR_RD_VCF",
      sv_per_shard=10000,
      samples_per_shard=100,
      Sanders_2015_tarball=Sanders_2015_tarball,
      Collins_2017_tarball=Collins_2017_tarball,
      Werling_2018_tarball=Werling_2018_tarball
  }

  # Subset inversions from PESR+RD VCF
  call subset_inversions {
    input:
      vcf=pesr_depth_overlap.merged_vcf,
      prefix=prefix
  }

  # Resolve variants - inversions only
  call resolve_complex_sv.resolve_complex_sv as resolve_cpx_inv {
    input:
      vcf=subset_inversions.inversion_vcf,
      contigs=contigs,
      max_shards_per_chrom=1,
      min_variants_per_shard=1,
      cytobands=cytobands,
      cytobands_idx=cytobands_idx,
      discfile_list=discfile_list,
      discfile_idx_list=discfile_idx_list,
      mei_bed=mei_bed,
      pe_blacklist=pe_blacklist,
      pe_blacklist_idx=pe_blacklist_idx,
      svc_acct_key=svc_acct_key,
      rf_cutoffs=rf_cutoffs
  }

  # Filter merged VCF for variant calls with identical breakpoints
  call bp_overlap_filter_by_chrom.same_bp_filter as bp_overlap_filter_by_chrom {
    input:
      vcf=pesr_depth_overlap.merged_vcf,
      vcf_idx=pesr_depth_overlap.merged_vcf_idx,
      prefix=prefix,
      contiglist=contigs,
      svc_acct_key=svc_acct_key,
      bothside_pass=update_bothside_pass_second.updated_list,
      background_fail=update_background_fail_second.updated_list
  }

  # Resolve variants - full VCF
  call resolve_complex_sv.resolve_complex_sv as resolve_cpx_all {
    input:
      vcf=bp_overlap_filter_by_chrom.filtered_vcf,
      contigs=contigs,
      max_shards_per_chrom=max_shards_per_chrom,
      min_variants_per_shard=min_variants_per_shard,
      cytobands=cytobands,
      cytobands_idx=cytobands_idx,
      discfile_list=discfile_list,
      discfile_idx_list=discfile_idx_list,
      mei_bed=mei_bed,
      pe_blacklist=pe_blacklist,
      pe_blacklist_idx=pe_blacklist_idx,
      svc_acct_key=svc_acct_key,
      rf_cutoffs=rf_cutoffs
  }

  # Integrate inv-only and all-variant resolved VCFs
  call integrate_resolved_vcfs.integrate_invonly_allvars as integrate_resolved {
    input:
      inv_res_vcf=resolve_cpx_inv.resolved_vcf_merged,
      all_res_vcf=resolve_cpx_all.resolved_vcf_merged,
      inv_res_vcf_idx=resolve_cpx_inv.resolved_vcf_merged_idx,
      all_res_vcf_idx=resolve_cpx_all.resolved_vcf_merged_idx,
      prefix=prefix,
      background_fail=update_background_fail_second.updated_list,
      bothside_pass=update_bothside_pass_second.updated_list,
      contiglist=contigs,
      svc_acct_key=svc_acct_key
  }

  # Apply consistent variant naming scheme to integrated VCF
  call rename {
    input:
      vcf=integrate_resolved.integrated_vcf,
      prefix=prefix
  }
  
  # Update SR background fail & bothside pass files
  call update_sr_list as update_background_fail_third {
    input:
      vcf=rename.renamed_vcf,
      original_list=update_background_fail_second.updated_list,
      outfile="sr_background_fail.updated3.txt"
  }

  call update_sr_list as update_bothside_pass_third {
    input:
      vcf=rename.renamed_vcf,
      original_list=update_bothside_pass_second.updated_list,
      outfile="sr_bothside_pass.updated3.txt"
  }

  # Prep input file for depth genotyping of complex intervals
  call make_cpx_cnv_input_file {
    input:
      batches_list=batches_list,
      depth_gt_rd_sep_list=depth_gt_rd_sep_list,
      bincov_list=bincov_list,
      bincov_idx_list=bincov_idx_list,
      famfile_list=famfile_list,
      medianfile_list=medianfile_list,
      sampleslist_list=sampleslist_list,
      prefix=prefix
  }

  # Depth-based genotyping of complex intervals
  call cpx_gt.genotype_CPX_CNVs as genotype_CPX_CNVs {
    input:
      vcf=rename.renamed_vcf,
      gt_input_files=make_cpx_cnv_input_file.CPX_CNV_gt_input,
      n_per_split_small=2500,
      n_per_split_large=250,
      n_RdTest_bins=100000,
      svc_acct_key=svc_acct_key,
      prefix=prefix,
      contigs=contigs,
      famfile=famfile
  }

  # QC VCF after complex resolution
  call vcf_qc.master_vcf_qc as postCPX_vcf_qc {
    input:
      vcf=genotype_CPX_CNVs.cpx_depth_gt_resolved_vcf,
      famfile=trios_famfile,
      ref_build=ref_build,
      prefix="${prefix}_postCPX_VCF",
      sv_per_shard=10000,
      samples_per_shard=100,
      Sanders_2015_tarball=Sanders_2015_tarball,
      Collins_2017_tarball=Collins_2017_tarball,
      Werling_2018_tarball=Werling_2018_tarball
  }

  # # Merge various VCF QC outputs into single tarball
  # call merge_vcf_qc {
  #   input:
  #     pesr_vcf_qc_tar=pesr_vcf_qc.sv_vcf_qc_output,
  #     rd_vcf_qc_tar=rd_vcf_qc.sv_vcf_qc_output,
  #     pesr_rd_vcf_qc_tar=pesr_rd_vcf_qc.sv_vcf_qc_output,
  #     postCPX_vcf_qc_tar=postCPX_vcf_qc.sv_vcf_qc_output,
  #     prefix=prefix
  # }

  # Final outputs
  output {
    File final_04_vcf = genotype_CPX_CNVs.cpx_depth_gt_resolved_vcf
    File final_04_vcf_qc = postCPX_vcf_qc.sv_vcf_qc_output
    # File CPX_CNVs_to_test = get_cpx_cnv_intervals.CPX_CNV_BED
    File pesr_vcf = cluster_pesr.clustered_vcf
    File depth_vcf = cluster_depth.clustered_vcf
    File pesr_depth_vcf = pesr_depth_overlap.merged_vcf
    # File complex_resolved = resolve_cpx_all.resolved_vcf_merged
    # File complex_unresolved = resolve_cpx_all.unresolved_vcf_merged
    # File merged_vcf_qc = merge_vcf_qc.merged_vcf_qc
    File updated_bothside_pass = update_bothside_pass_third.updated_list
    File updated_background_fail = update_background_fail_third.updated_list
    # File cpx_genotypes = genotype_CPX_CNVs.CPX_interval_genotypes
  }
}


# Update either SR bothside_pass or background_fail files
task update_sr_list {
  File vcf
  File original_list
  String outfile

  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/trackpesr_ID.sh \
      ${vcf} \
      ${original_list} \
      ${outfile}
  >>>

  output {
    File updated_list = "${outfile}"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:b359f2cb0c9d5f5a55eb4c41fd362f4e574bf3f8f0f395a2907837571b367ee0"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
  }
}


# Subset inversions from PESR VCF
task subset_inversions {
  File vcf
  String prefix

  command <<<
    #Write header
    zcat ${vcf} | fgrep "#" > ${prefix}.inversions_only.vcf;
    #Get inversions
    zcat ${vcf} | fgrep -v "#" | fgrep "SVTYPE=INV" \
      >> ${prefix}.inversions_only.vcf;
    #Sort & compress output
    cat ${prefix}.inversions_only.vcf \
      | vcf-sort \
      | bgzip -c \
      >> ${prefix}.inversions_only.vcf.gz
  >>>

  output {
    File inversion_vcf = "${prefix}.inversions_only.vcf.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:b359f2cb0c9d5f5a55eb4c41fd362f4e574bf3f8f0f395a2907837571b367ee0"
    preemptible: 1
    disks: "local-disk 500 SSD"
  }
}


# Rename variants in VCF
task rename {
  File vcf
  String prefix

  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/rename.py \
      --prefix ${prefix} ${vcf} - \
      | bgzip -c > "${prefix}.04_renamed.vcf.gz"
  >>>

  output {
    File renamed_vcf = "${prefix}.04_renamed.vcf.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:b359f2cb0c9d5f5a55eb4c41fd362f4e574bf3f8f0f395a2907837571b367ee0"
    preemptible: 1
    memory: "8 GB"
    disks: "local-disk 1000 SSD"
  }
}


# Create input file for per-batch genotyping of predicted CPX CNV intervals
task make_cpx_cnv_input_file {
  File batches_list
  File bincov_list
  File bincov_idx_list
  File depth_gt_rd_sep_list
  File famfile_list
  File medianfile_list
  File sampleslist_list
  String prefix

  command <<<
    paste ${batches_list} ${bincov_list} ${bincov_idx_list} ${depth_gt_rd_sep_list} \
    ${sampleslist_list} ${famfile_list} ${medianfile_list} > \
    ${prefix}.cpx_cnv_genotyping_input.txt
  >>>

  output {
    File CPX_CNV_gt_input = "${prefix}.cpx_cnv_genotyping_input.txt"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:b359f2cb0c9d5f5a55eb4c41fd362f4e574bf3f8f0f395a2907837571b367ee0"
    preemptible: 1
    memory: "8 GB"
    disks: "local-disk 500 SSD"
  }
}


# Merge VCF QC outputs
task merge_vcf_qc {
  File pesr_vcf_qc_tar
  File rd_vcf_qc_tar
  File pesr_rd_vcf_qc_tar
  File postCPX_vcf_qc_tar
  String prefix

  command <<<
    mkdir "${prefix}_merged_vcf_qc_04_integrate_batches/"
    cp ${pesr_vcf_qc_tar} ${rd_vcf_qc_tar} \
      ${pesr_rd_vcf_qc_tar} ${postCPX_vcf_qc_tar} \
      ${prefix}_merged_vcf_qc_04_integrate_batches/
    tar -czvf ${prefix}_merged_vcf_qc_04_integrate_batches.tar.gz \
      ${prefix}_merged_vcf_qc_04_integrate_batches
  >>>

  output {
    File merged_vcf_qc = "${prefix}_merged_vcf_qc_04_integrate_batches.tar.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:b359f2cb0c9d5f5a55eb4c41fd362f4e574bf3f8f0f395a2907837571b367ee0"
    preemptible: 1
    disks: "local-disk 500 SSD"
  }
}
