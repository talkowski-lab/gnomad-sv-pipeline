# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


#This is an analysis WDL to perform PCA-based ancestry assignment and infer
# relatedness from a SV VCF output by the Talkowski Lab SV pipeline


workflow get_ancestries_and_relatedness {
  File vcf
  File vcf_idx
  File sample_pop_training_labels
  File sample_batch_assignments
  File PCRPLUS_samples_list
  File famfile
  File autosomes_list
  File trios_famfile
  Float min_nonnull_gt_frac_forPCA
  Float min_nonnull_gt_frac_forKING
  Float min_maf_forPCA
  Float min_maf_forKING
  Int min_qual
  String prefix

  Array[Array[String]] contigs = read_tsv(autosomes_list)

  #Filter VCF to common autosomal PASS variants
  scatter ( contig in contigs ) {
    call filter_vcf as prePCA_vcf_filter {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        contig=contig[0],
        min_maf=min_maf_forPCA,
        min_qual=min_qual,
        min_nonnull_gt_frac=min_nonnull_gt_frac_forPCA,
        prefix="${prefix}.${contig[0]}.filtered.common"
    }
  }
  call concat_vcfs as concat_prePCA_vcf {
    input:
      vcfs=prePCA_vcf_filter.filtered_vcf,
      outfile_prefix="${prefix}.filtered.common"
  }

  #Compute LD-pruned GRM from filtered VCF of common variants only (for PCA)
  call createGRM {
    input:
      vcf=concat_prePCA_vcf.concat_vcf,
      max_LD=0.2,
      LD_prune_distance=1000,
      prefix=prefix
  }

  #Perform PCA and assign ancestry labels from GRM of common variants
  call PCA_assign_pops {
    input:
      GRM=createGRM.GRM,
      sample_pop_training_labels=sample_pop_training_labels,
      sample_batch_assignments=sample_batch_assignments,
      PCRPLUS_samples_list=PCRPLUS_samples_list,
      prefix=prefix
  }

  #Filter VCF to all autosomal PASS variants (with looser AF restriction)
  scatter ( contig in contigs ) {
    call filter_vcf as preKING_vcf_filter {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        contig=contig[0],
        min_maf=min_maf_forKING,
        min_qual=min_qual,
        min_nonnull_gt_frac=min_nonnull_gt_frac_forKING,
        prefix="${prefix}.${contig[0]}.filtered.looser"
    }
  }
  call concat_vcfs as concat_preKING_vcf {
    input:
      vcfs=preKING_vcf_filter.filtered_vcf,
      outfile_prefix="${prefix}.filtered.looser"
  }

  #Make PLINK BED, FAM, and BIM files from filtered VCF with no AF restrictions
  call make_plink_files {
    input:
      vcf=concat_preKING_vcf.concat_vcf,
      prefix=prefix
  }

  #Infer relatedness with KING
  call run_king {
    input:
      plink_bed=make_plink_files.plink_bed,
      plink_fam=make_plink_files.plink_fam,
      plink_bim=make_plink_files.plink_bim,
      prefix=prefix
  }
  call infer_relateds {
    input:
      king_metrics_filtered=run_king.king_metrics_filtered,
      king_metrics_all=run_king.king_metrics_all,
      trios_famfile=trios_famfile,
      n_unrelated_pairs=10000,
      prefix=prefix
  }

  output {
    File sample_pop_assignments = PCA_assign_pops.new_population_assignments
    File PCA_loadings = PCA_assign_pops.PCA_loadings
    File related_samples_to_prune = infer_relateds.related_samples_to_prune
  }
}


# Generate VCFs of autosomal PASS variants
task filter_vcf {
  File vcf
  File vcf_idx
  String contig
  Float min_maf
  Float min_nonnull_gt_frac
  Float min_qual
  String prefix

  command <<<
    tabix -h ${vcf} ${contig} \
    | vcftools \
      --vcf - \
      --chr ${contig} \
      --keep-filtered PASS \
      --remove-INFO PCRPLUS_DEPLETED \
      --remove-INFO UNSTABLE_AF_PCRPLUS \
      --remove-INFO VARIABLE_ACROSS_BATCHES \
      --max-missing ${min_nonnull_gt_frac} \
      --maf ${min_maf} \
      --minQ ${min_qual} \
      --recode \
      --recode-INFO-all \
      --stdout \
    | bgzip -c \
    > "${prefix}.vcf.gz"
  >>>

  output {
    File filtered_vcf = "${prefix}.vcf.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:50875d61110ab2139f657c7fa1311fa9be9c7df11fe10c3ab08ecfd8e06e5250"
    preemptible: 1
    maxRetries: 1
    memory: "4 GB"
    disks: "local-disk 250 HDD"
  }
}


# Combine multiple VCFs without sorting
task concat_vcfs {
  Array[File] vcfs
  String outfile_prefix

  command <<<
    vcf-concat ${sep=' ' vcfs} | bgzip -c > ${outfile_prefix}.vcf.gz; 
    tabix -p vcf -f "${outfile_prefix}.vcf.gz"
  >>>

  output {
    File concat_vcf = "${outfile_prefix}.vcf.gz"
    File concat_vcf_idx = "${outfile_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:50875d61110ab2139f657c7fa1311fa9be9c7df11fe10c3ab08ecfd8e06e5250"
    preemptible: 0
    maxRetries: 1
    memory: "8 GB"
    disks: "local-disk 250 HDD"
  }
}


# Create genetic relatedness matrix (GRM) of allele dosage
task createGRM {
  File vcf
  Float max_LD
  Float LD_prune_distance
  String prefix

  command <<<
    # Prune input VCF on LD
    bcftools +prune \
      --max-LD ${max_LD} \
      --window ${LD_prune_distance}kb \
      -Oz \
      -o "${prefix}.LD_pruned.vcf.gz" \
      ${vcf}
    # Then compute GRM
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/vcf2grm.py \
      "${prefix}.LD_pruned.vcf.gz" stdout \
      | bgzip -c \
      > "${prefix}.grm.txt.gz"
  >>>

  output {
    File GRM = "${prefix}.grm.txt.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:50875d61110ab2139f657c7fa1311fa9be9c7df11fe10c3ab08ecfd8e06e5250"
    preemptible: 1
    maxRetries: 1
    memory: "4 GB"
    disks: "local-disk 50 HDD"
  }
}


# Run PCA on GRM, use SVM to classify populations, and label samples
task PCA_assign_pops {
  File GRM
  File sample_pop_training_labels
  File sample_batch_assignments
  File PCRPLUS_samples_list
  String prefix

  command <<<
    mkdir "${prefix}_PCA_plots"
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/runPCA_labelAncestries.R \
      -p "./${prefix}_PCA_plots/${prefix}" \
      --batchAssignments ${sample_batch_assignments} \
      --PCRPLUSsamples ${PCRPLUS_samples_list} \
      --confidence 0.80 \
      ${GRM} \
      ${sample_pop_training_labels} \
      /opt/sv-pipeline/ref/gnomAD_population_colors.txt
    cp ./${prefix}_PCA_plots/${prefix}.PCA_loadings.txt.gz \
      ${prefix}.PCA_loadings.txt.gz
    cp ./${prefix}_PCA_plots/${prefix}.new_population_labels.txt.gz \
      ${prefix}.new_population_labels.txt.gz
    zcat ${prefix}.new_population_labels.txt.gz \
    | fgrep -v "#" \
    | cut -f1-2 \
    > ${prefix}.PCA_population_assignments.txt
    tar -czvf "${prefix}_PCA_plots.tar.gz" \
      "./${prefix}_PCA_plots"
  >>>

  output {
    File pca_plots_tarball = "${prefix}_PCA_plots.tar.gz"
    File new_population_assignments = "${prefix}.PCA_population_assignments.txt"
    File PCA_loadings = "${prefix}.PCA_loadings.txt.gz"
    File SVM_classifications = "${prefix}.new_population_labels.txt.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:e644bde07541047bd1e89dfa38dca76e33a29d9ac66940e321207447a0537003"
    preemptible: 1
    maxRetries: 1
    memory: "32 GB"
    bootDiskSizeGb: 30
    disks: "local-disk 100 HDD"
  }
}


# Make plink files required for KING
task make_plink_files {
  File vcf
  String prefix

  command <<<
    plink2 --vcf ${vcf} --make-bed --out "${prefix}.plink_input"
  >>>

  output {
    File plink_bed = "${prefix}.plink_input.bed"
    File plink_fam = "${prefix}.plink_input.fam"
    File plink_bim = "${prefix}.plink_input.bim"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:50875d61110ab2139f657c7fa1311fa9be9c7df11fe10c3ab08ecfd8e06e5250"
    preemptible: 1
    maxRetries: 1
    memory: "4 GB"
    disks: "local-disk 30 HDD"
  }
}


# Run KING
task run_king {
  File plink_bed
  File plink_fam
  File plink_bim
  String prefix

  command <<<
    set -euo pipefail
    #Run KING
    king -b ${plink_bed} --fam ${plink_fam} --bim ${plink_bim} --related --degree 2
    #Clean output for all samples
    cut -f2-4,7-14 king.kin | gzip -c > ${prefix}.king_metrics.txt.gz
    # #Filter on HetConc > 0.2, IBS0 < 0.006, and Kinship > 0.1
    #Assumes all pairs not meeting these criteria are, at best, distant relatives
    zcat ${prefix}.king_metrics.txt.gz \
    | sed -n '1p' | sed 's/\t/\n/g' \
    | awk -v OFS="\t" '{ print $1, NR }' \
    > header_indexes.txt
    HetConc_idx=$( fgrep -w HetConc header_indexes.txt | cut -f2 )
    HetHet_idx=$( fgrep -w HetHet header_indexes.txt | cut -f2 )
    IBS0_idx=$( fgrep -w IBS0 header_indexes.txt | cut -f2 )
    Kinship_idx=$( fgrep Kinship header_indexes.txt | cut -f2 )
    zcat ${prefix}.king_metrics.txt.gz \
    | awk -v OFS="\t" \
      -v HetConc="$HetConc_idx" \
      -v HetHet="$HetHet_idx" \
      -v IBS0="$IBS0_idx" \
      -v Kinship="$Kinship_idx" \
      '{ if ($1=="ID1" || ($(HetConc)>0.2 && $(IBS0)<0.006 $(Kinship)>0.1)) print $0 }' \
    | gzip -c \
    > ${prefix}.king_metrics.candidate_pairs_subsetted.txt.gz
  >>>

  output {
    File king_metrics_all = "${prefix}.king_metrics.txt.gz"
    File king_metrics_filtered = "${prefix}.king_metrics.candidate_pairs_subsetted.txt.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:ee412e454f2183d270698844a5e10149b97d919592a521b5ecbba65225ded9e4"
    preemptible: 1
    maxRetries: 1
    memory: "32 GB"
    disks: "local-disk 50 HDD"
  }
}


# Infer relatedness between samples from KING output
task infer_relateds {
  File king_metrics_filtered
  File king_metrics_all
  File trios_famfile
  Int n_unrelated_pairs
  String prefix

  command <<<
    set -euo pipefail
    # Select random pairs of samples not known to have any relatives to spike into SVM training
    awk -v FS="\t" -v OFS="\n" '{ print $2, $3, $4 }' ${trios_famfile} \
    > samples_in_families.txt
    zcat ${king_metrics_all} \
    | cut -f1,2 \
    | fgrep -wvf samples_in_families.txt \
    | shuf --random-source=${king_metrics_all} \
    | head -n ${n_unrelated_pairs} \
    > unrelated_sample_pairs.supplement.txt || true
    #Parse KING results and determine which samples should be pruned
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/parse_KING_results.R \
      -p ${prefix} \
      -u unrelated_sample_pairs.supplement.txt \
      ${king_metrics_filtered} \
      ${trios_famfile}
    #Tarball all king plots
    mkdir ${prefix}_KING_plots
    mv *.pdf ${prefix}_KING_plots/
    tar -czvf ${prefix}_KING_plots.tar.gz ${prefix}_KING_plots
  >>>

  output {
    File related_samples_to_prune = "${prefix}.related_samles_to_prune.txt"
    File inferred_relationships = "${prefix}.inferred_pairwise_relationships.txt.gz"
    File KING_inference_accuracy = "${prefix}.relationship_inference_accuracy.txt"
    File KING_plots = "${prefix}_KING_plots.tar.gz"
    File unrelated_sample_supplement = "unrelated_sample_pairs.supplement.txt"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:ee412e454f2183d270698844a5e10149b97d919592a521b5ecbba65225ded9e4"
    preemptible: 1
    maxRetries: 1
    memory: "32 GB"
    disks: "local-disk 50 HDD"
  }
}

