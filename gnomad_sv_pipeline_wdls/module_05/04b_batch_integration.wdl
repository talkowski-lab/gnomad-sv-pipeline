# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


#This is the third rebuild of module 04b from Sept 2018, which has been modified
# to scatter once by chromosome at the beginning of the workflow, run all steps
# in 04b, and gather across chromosomes once at the end.


#Imports:
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04b_vcfcluster_single_chrom/versions/11/plain-WDL/descriptor" as clust
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04b_resolve_complex_sv/versions/28/plain-WDL/descriptor" as resolve
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04b_scatter_CPX_genotyping/versions/12/plain-WDL/descriptor" as cpx_gt
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:master_SV_VCF_QC/versions/75/plain-WDL/descriptor" as vcf_qc


# Master workflow to parallelize batch integration per chromosome
workflow integrate_batches {
  #Inputs
  File pesr_vcf_list
  File pesr_vcf_idx_list
  File depth_vcf_list
  File depth_vcf_idx_list
  File blacklist
  File blacklist_idx
  Array[String] samples
  File contiglist
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

  Array[Array[String]] contigs = read_tsv(contiglist)

  #Prep input file for depth genotyping of complex intervals
  call make_cpx_cnv_input_file {
    input:
      batches_list=batches_list,
      depth_gt_rd_sep_list=depth_gt_rd_sep_list,
      bincov_list=bincov_list,
      bincov_idx_list=bincov_idx_list,
      famfile_list=famfile_list,
      medianfile_list=medianfile_list,
      sampleslist_list=sampleslist_list,
      prefix="${prefix}"
  }

  #Scatter per chromosome
  scatter ( contig in contigs ) {

    #Subset PESR VCFs to single chromosome & cluster
    #Note: also subsets bothside_pass and background_fail files to variants 
    #present on chromosome of interest
    call clust.vcfcluster_single_chrom as cluster_pesr {
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
        contig=contig[0],
        max_shards_per_chrom_svtype=100,
        min_variants_per_shard_per_chrom_svtype=100,
        svc_acct_key=svc_acct_key,
        subset_sr_lists="TRUE",
        bothside_pass=sr_bothend_pass,
        background_fail=sr_background_fail
    }

    #Subset RD VCFs to single chromosome & cluster
    call clust.vcfcluster_single_chrom as cluster_depth {
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
        contig=contig[0],
        max_shards_per_chrom_svtype=100,
        min_variants_per_shard_per_chrom_svtype=100,
        svc_acct_key=svc_acct_key,
        subset_sr_lists="FALSE",
        bothside_pass=sr_bothend_pass,
        background_fail=sr_background_fail
    }

    #Update SR background fail & bothside pass files (1)
    call update_sr_list as update_background_fail {
      input:
        vcf=cluster_pesr.clustered_vcf,
        original_list=cluster_pesr.filtered_bothside_pass,
        outfile="sr_background_fail.${contig[0]}.updated.txt"
    }
    call update_sr_list as update_bothside_pass {
      input:
        vcf=cluster_pesr.clustered_vcf,
        original_list=cluster_pesr.filtered_background_fail,
        outfile="sr_bothside_pass.${contig[0]}.updated.txt"
    }

    #Merge PESR & RD VCFs
    call merge_pesr_depth {
      input:
        pesr_vcf=cluster_pesr.clustered_vcf,
        depth_vcf=cluster_depth.clustered_vcf,
        contig=contig[0]
    }

    #Update SR background fail & bothside pass files (2)
    call update_sr_list as update_background_fail_second {
      input:
        vcf=merge_pesr_depth.merged_vcf,
        original_list=update_background_fail.updated_list,
        outfile="sr_background_fail.${contig[0]}.updated2.txt"
    }
    call update_sr_list as update_bothside_pass_second {
      input:
        vcf=merge_pesr_depth.merged_vcf,
        original_list=update_bothside_pass.updated_list,
        outfile="sr_bothside_pass.${contig[0]}.updated2.txt"
    }

    #Subset inversions from PESR+RD VCF
    call subset_inversions {
      input:
        vcf=merge_pesr_depth.merged_vcf,
        prefix="${prefix}.${contig[0]}"
    }

    #Resolve inversion-only VCF
    call resolve.resolve_complex_sv as resolve_cpx_inv {
      input:
        vcf=subset_inversions.inversion_vcf,
        vcf_idx=subset_inversions.inversion_vcf_idx,
        prefix="${prefix}.inv_only",
        contig=contig[0],
        max_shards_per_chrom=max_shards_per_chrom,
        min_variants_per_shard=100,
        cytobands=cytobands,
        cytobands_idx=cytobands_idx,
        discfile_list=discfile_list,
        discfile_idx_list=discfile_idx_list,
        mei_bed=mei_bed,
        pe_blacklist=pe_blacklist,
        pe_blacklist_idx=pe_blacklist_idx,
        svc_acct_key=svc_acct_key,
        rf_cutoffs=rf_cutoffs,
        inv_only="TRUE"
    }

    #Run same-bp overlap filter on full vcf
    call bp_overlap_filter {
      input:
        vcf=merge_pesr_depth.merged_vcf,
        prefix="${prefix}.${contig[0]}",
        bothside_pass=update_bothside_pass_second.updated_list,
        background_fail=update_background_fail_second.updated_list
    }

    #Resolve all-variants VCF after same-bp overlap filter
    call resolve.resolve_complex_sv as resolve_cpx_all {
      input:
        vcf=bp_overlap_filter.bp_filtered_vcf,
        vcf_idx=bp_overlap_filter.bp_filtered_vcf_idx,
        prefix="${prefix}.all_variants",
        contig=contig[0],
        max_shards_per_chrom=max_shards_per_chrom,
        min_variants_per_shard=100,
        cytobands=cytobands,
        cytobands_idx=cytobands_idx,
        discfile_list=discfile_list,
        discfile_idx_list=discfile_idx_list,
        mei_bed=mei_bed,
        pe_blacklist=pe_blacklist,
        pe_blacklist_idx=pe_blacklist_idx,
        svc_acct_key=svc_acct_key,
        rf_cutoffs=rf_cutoffs,
        inv_only="FALSE"
    }

    #Integrate inv-only and all-variants resolved VCFs
    call integrate_resolved_vcfs {
      input:
        inv_res_vcf=resolve_cpx_inv.resolved_vcf_merged,
        all_res_vcf=resolve_cpx_all.resolved_vcf_merged,
        prefix="${prefix}.resolved.${contig[0]}"
    }

    #Apply consistent variant naming scheme to integrated VCF
    call rename {
      input:
        vcf=integrate_resolved_vcfs.integrated_vcf,
        prefix="${prefix}.${contig[0]}"
    }
    
    #Update SR background fail & bothside pass files
    call update_sr_list as update_background_fail_third {
      input:
        vcf=rename.renamed_vcf,
        original_list=update_background_fail_second.updated_list,
        outfile="sr_background_fail.${contig[0]}.updated3.txt"
    }
    call update_sr_list as update_bothside_pass_third {
      input:
        vcf=rename.renamed_vcf,
        original_list=update_bothside_pass_second.updated_list,
        outfile="sr_bothside_pass.${contig[0]}.updated3.txt"
    }

    #Depth-based genotyping of complex intervals
    call cpx_gt.scatter_CPX_genotyping as genotype_CPX_CNVs {
      input:
        vcf=rename.renamed_vcf,
        n_master_vcf_shards=200,
        n_master_min_vars_per_vcf_shard=5000,
        gt_input_files=make_cpx_cnv_input_file.CPX_CNV_gt_input,
        n_per_split_small=2500,
        n_per_split_large=250,
        n_RdTest_bins=100000,
        svc_acct_key=svc_acct_key,
        prefix=prefix,
        contig=contig[0],
        famfile=famfile
    }
  }

  #Merge PESR+RD VCFs for midpoint QC
  call concat_vcfs as concat_midpoint_vcfs {
    input:
      vcfs=merge_pesr_depth.merged_vcf,
      outfile_prefix="${prefix}.pesr_rd_merged"
  }

  #Run midpoint QC on merged PESR+RD VCF across all chromosomes
  call vcf_qc.master_vcf_qc as midpoint_qc {
    input:
      vcf=concat_midpoint_vcfs.concat_vcf,
      vcf_idx=concat_midpoint_vcfs.concat_vcf_idx,
      famfile=trios_famfile,
      ref_build=ref_build,
      prefix="${prefix}_pesr_rd_merged_VCF",
      sv_per_shard=10000,
      samples_per_shard=100,
      Sanders_2015_tarball=Sanders_2015_tarball,
      Collins_2017_tarball=Collins_2017_tarball,
      Werling_2018_tarball=Werling_2018_tarball,
      contiglist=contiglist
  }

  #Merge final resolved vcfs
  call concat_vcfs as concat_final_vcfs {
    input:
      vcfs=genotype_CPX_CNVs.cpx_depth_gt_resolved_vcf,
      outfile_prefix="${prefix}.resolved_regenotyped"
  }

  #Run final QC on resolved VCF across all chromosomes
  call vcf_qc.master_vcf_qc as final_qc {
    input:
      vcf=concat_final_vcfs.concat_vcf,
      vcf_idx=concat_final_vcfs.concat_vcf_idx,
      famfile=trios_famfile,
      ref_build=ref_build,
      prefix="${prefix}_resolved_VCF",
      sv_per_shard=10000,
      samples_per_shard=100,
      Sanders_2015_tarball=Sanders_2015_tarball,
      Collins_2017_tarball=Collins_2017_tarball,
      Werling_2018_tarball=Werling_2018_tarball,
      contiglist=contiglist
  }

  #Merge SR background fail & bothside pass files
  call cat_vid_lists as cat_background_fail {
    input:
      vid_lists=update_background_fail_third.updated_list,
      outfile="${prefix}.sr_background_fail.txt"
  }
  call cat_vid_lists as cat_bothside_pass {
    input:
      vid_lists=update_bothside_pass_third.updated_list,
      outfile="${prefix}.sr_bothside_pass.txt"
  }

  #Final outputs
  output {
    File final_04b_vcf = concat_final_vcfs.concat_vcf
    File final_04b_vcf_idx = concat_final_vcfs.concat_vcf_idx
    File final_04b_vcf_qc = final_qc.sv_vcf_qc_output
    File updated_sr_background_fail = cat_background_fail.merged_file
    File updated_sr_bothside_pass = cat_bothside_pass.merged_file
  }
}


#Create input file for per-batch genotyping of predicted CPX CNV intervals
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
    docker: "talkowski/sv-pipeline@sha256:86855b7f43894e2d3a9f61c0e4bf3782746a9282daad60559f58ece978680c9b"
    preemptible: 1
    disks: "local-disk 10 SSD"
  }
}


#Update either SR bothside_pass or background_fail files
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
    docker: "talkowski/sv-pipeline@sha256:86855b7f43894e2d3a9f61c0e4bf3782746a9282daad60559f58ece978680c9b"
    preemptible: 1
    disks: "local-disk 50 SSD"
  }
}


#Merge PESR + RD VCFs
task merge_pesr_depth {
  File pesr_vcf
  File depth_vcf
  String contig
  
  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/PESR_RD_merge_wrapper.sh \
      ${pesr_vcf} \
      ${depth_vcf} \
      ${contig} \
      all_batches.pesr_depth.${contig}.vcf.gz;
    tabix -p vcf -f all_batches.pesr_depth.${contig}.vcf.gz
  >>>

  output {
    File merged_vcf = "all_batches.pesr_depth.${contig}.vcf.gz"
    File merged_vcf_idx = "all_batches.pesr_depth.${contig}.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:3f9d99b8154dff67eb33b0da0a4358ac149461d65f819e7eb64958953d478900"
    preemptible: 1
    memory: "32 GB"
    disks: "local-disk 500 SSD"
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
      >> ${prefix}.inversions_only.vcf.gz;
    #Index output
    tabix -p vcf -f "${prefix}.inversions_only.vcf.gz"
  >>>

  output {
    File inversion_vcf = "${prefix}.inversions_only.vcf.gz"
    File inversion_vcf_idx = "${prefix}.inversions_only.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:40e6a9e956302d32137d0c4a2f33779d23a70ea603d88a8ec03653090bb70107"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 300 SSD"
  }
}


#Run Harrison's overlapping breakpoint filter prior to complex resolution
task bp_overlap_filter {
  File vcf
  String prefix
  File bothside_pass
  File background_fail

  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/overlapbpchange.sh \
    ${vcf} \
    ${background_fail} \
    ${bothside_pass};
    mv non_redundant.vcf.gz "${prefix}.non_redundant.vcf.gz";
    tabix -p vcf -f "${prefix}.non_redundant.vcf.gz"
  >>>

  output {
    File bp_filtered_vcf = "${prefix}.non_redundant.vcf.gz"
    File bp_filtered_vcf_idx = "${prefix}.non_redundant.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:40e6a9e956302d32137d0c4a2f33779d23a70ea603d88a8ec03653090bb70107"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 300 SSD"
  }
}


#Merge inversion-only and all-variant cpx-resolved outputs
task integrate_resolved_vcfs {
  File inv_res_vcf
  File all_res_vcf
  String prefix

  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/Complex_Inversion_Integration.sh \
      ${inv_res_vcf} \
      ${all_res_vcf} \
      ${prefix}.integrated_resolved.vcf.gz;
      tabix -p vcf -f "${prefix}.integrated_resolved.vcf.gz"
  >>>

  output {
    File integrated_vcf = "${prefix}.integrated_resolved.vcf.gz"
    File integrated_vcf_idx = "${prefix}.integrated_resolved.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:40e6a9e956302d32137d0c4a2f33779d23a70ea603d88a8ec03653090bb70107"
    preemptible: 1
    memory: "8 GB"
    disks: "local-disk 300 SSD"
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
    docker: "talkowski/sv-pipeline@sha256:40e6a9e956302d32137d0c4a2f33779d23a70ea603d88a8ec03653090bb70107"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
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
    docker: "talkowski/sv-pipeline@sha256:86855b7f43894e2d3a9f61c0e4bf3782746a9282daad60559f58ece978680c9b"
    preemptible: 1
    memory: "8 GB"
    disks: "local-disk 1000 SSD"
  }
}


#Cat & sort VID lists
task cat_vid_lists {
  Array[File] vid_lists
  String outfile

  command <<<
    cat ${sep=' ' vid_lists} | sort -Vk1,1 > ${outfile}
  >>>

  output {
    File merged_file = "${outfile}"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:86855b7f43894e2d3a9f61c0e4bf3782746a9282daad60559f58ece978680c9b"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 20 SSD"
  }  
}
