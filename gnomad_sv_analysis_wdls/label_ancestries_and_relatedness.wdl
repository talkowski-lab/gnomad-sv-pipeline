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
  String prefix

  Array[Array[String]] contigs = read_tsv(autosomes_list)

  #Filter VCF to common autosomal PASS variants
  scatter ( contig in contigs ) {
    call filter_vcf as prePCA_vcf_filter {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        contig=contig[0],
        min_maf=0.01,
        prefix="${prefix}.${contig[0]}.filtered.common"
    }
  }
  call concat_vcfs as concat_prePCA_vcf {
    input:
      vcfs=prePCA_vcf_filter.filtered_vcf,
      outfile_prefix="${prefix}.filtered.common"
  }

  #Compute GRM from filtered VCF of common variants only (for PCA)
  call createGRM {
    input:
      vcf=concat_prePCA_vcf.concat_vcf,
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
        min_maf=0.00015,
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
  call infer_relateds {
    input:
      plink_bed=make_plink_files.plink_bed,
      plink_fam=make_plink_files.plink_fam,
      plink_bim=make_plink_files.plink_bim,
      trios_famfile=trios_famfile,
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
  String prefix

  command <<<
    tabix -h ${vcf} ${contig} \
    | vcftools \
      --vcf - \
      --chr ${contig} \
      --keep-filtered PASS \
      --remove-INFO PCRPLUS_DEPLETED \
      --remove-INFO PESR_GT_OVERDISPERSION \
      --max-missing 0.99 \
      --maf ${min_maf} \
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
    docker: "talkowski/sv-pipeline@sha256:25eebd3a2dfeaaf94fbfa85b30ecbbeeafb92633edcefa93b178c053317fcd8b"
    preemptible: 1
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
    docker: "talkowski/sv-pipeline@sha256:25eebd3a2dfeaaf94fbfa85b30ecbbeeafb92633edcefa93b178c053317fcd8b"
    preemptible: 0
    memory: "8 GB"
    disks: "local-disk 250 HDD"
  }
}


# Create genetic relatedness matrix (GRM) of allele dosage
task createGRM {
  File vcf
  String prefix

  command <<<
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/vcf2grm.py \
      ${vcf} stdout \
      | bgzip -c \
      > "${prefix}.grm.txt.gz"
  >>>

  output {
    File GRM = "${prefix}.grm.txt.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:25eebd3a2dfeaaf94fbfa85b30ecbbeeafb92633edcefa93b178c053317fcd8b"
    preemptible: 1
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
      --confidence 0.85 \
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
    docker: "talkowski/sv-pipeline@sha256:c2af5febc8967dff0b7a10cd764b292f43029ffd119e40832cef3fcbc3df1c1f"
    preemptible: 1
    memory: "8 GB"
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
    docker: "talkowski/sv-pipeline@sha256:c2af5febc8967dff0b7a10cd764b292f43029ffd119e40832cef3fcbc3df1c1f"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 30 HDD"
  }
}


# Infer relatedness between samples with KING
task infer_relateds {
  File plink_bed
  File plink_fam
  File plink_bim
  File trios_famfile
  String prefix

  command <<<
    #Run KING
    king -b ${plink_bed} --fam ${plink_fam} --bim ${plink_bim} --related --degree 2
    #Clean output for all samples
    cut -f2-4,7-14 king.kin | gzip -c > ${prefix}.king_metrics.txt.gz
    #Filter on HetConc > 0.25, IBD1Seg > 0.15, and Kinship > 0.05
    #Assumes all pairs not meeting these criteria are, at best, distant relatives
    zcat ${prefix}.king_metrics.txt.gz \
    | sed -n '1p' | sed 's/\t/\n/g' \
    | awk -v OFS="\t" '{ print $1, NR }' \
    > header_indexes.txt
    HetConc_idx=$( fgrep HetConc header_indexes.txt | cut -f2 )
    IBD1Seg_idx=$( fgrep IBD1Seg header_indexes.txt | cut -f2 )
    Kinship_idx=$( fgrep Kinship header_indexes.txt | cut -f2 )
    zcat ${prefix}.king_metrics.txt.gz \
    | awk -v OFS="\t" \
      -v HetConc="$HetConc_idx" \
      -v IBD1Seg="$IBD1Seg_idx" \
      -v Kinship="$Kinship_idx" \
      '{ if ($1=="ID1" || ($(HetConc)>0.25 && $(IBD1Seg)>0.15 && $(Kinship)>0.05)) print $0 }' \
    | gzip -c \
    > ${prefix}.king_metrics.candidate_pairs_subsetted.txt.gz
    #Parse KING results and determine which samples should be pruned
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/parse_KING_results.R \
      -p ${prefix} \
      ${prefix}.king_metrics.candidate_pairs_subsetted.txt.gz \
      ${trios_famfile}
  >>>

  output {
    File related_samples_to_prune = "${prefix}.related_samles_to_prune.txt"
    File inferred_relationships = "${prefix}.inferred_pairwise_relationships.txt.gz"
    File KING_inference_accuracy = "${prefix}.relationship_inference_accuracy.txt"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:64d98e5d6f85d34ba14a32aa615d69cfe303df9c952bdfeb86ca0d7737856e41"
    preemptible: 1
    memory: "32 GB"
    disks: "local-disk 50 HDD"
  }
}

