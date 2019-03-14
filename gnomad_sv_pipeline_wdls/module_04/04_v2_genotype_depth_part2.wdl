workflow genotype_depth_part2 {
  File cohort_vcf
  File RD_pesr_sepcutoff
  File RD_depth_sepcutoff
  Int n_per_split
  Int n_RdTest_bins
  String batch

  File medianfile
  File famfile
  File svc_acct_key
  Array[String] samples

  String coveragefile
  File coveragefile_idx
  
  call split_variants {
    input:
      vcf=cohort_vcf,
      n_per_split=n_per_split
  }

  scatter (lt5kb_bed in split_variants.lt5kb_beds) {
    call RdTest_genotype as RD_genotype_lt5kb {
      input:
        bed=lt5kb_bed,
        coveragefile=coveragefile,
        coveragefile_idx=coveragefile_idx,
        svc_acct_key=svc_acct_key,
        medianfile=medianfile,
        famfile=famfile,
        samples=samples,
        gt_cutoffs=RD_pesr_sepcutoff,
        n_bins=n_RdTest_bins,
        prefix=basename(lt5kb_bed, ".bed")
    }

    call make_subset_vcf as make_subset_vcf_lt5kb {
      input:
        vcf=cohort_vcf,
        bed=lt5kb_bed
    }

    call integrate_depth_GQ as integrate_GQ_lt5kb {
      input:
        vcf=make_subset_vcf_lt5kb.subset_vcf,
        RD_melted_genotypes=RD_genotype_lt5kb.melted_genotypes,
        RD_vargq=RD_genotype_lt5kb.varGQ,
    }

    call add_genotypes as add_genotypes_lt5kb {
      input:
        vcf=make_subset_vcf_lt5kb.subset_vcf,
        genotypes=integrate_GQ_lt5kb.genotypes,
        varGQ=integrate_GQ_lt5kb.varGQ,
        prefix=basename(lt5kb_bed, ".bed")
    }
  }
  
  scatter (gt5kb_bed in split_variants.gt5kb_beds) {
    call RdTest_genotype as RD_genotype_gt5kb {
      input:
        bed=gt5kb_bed,
        coveragefile=coveragefile,
        coveragefile_idx=coveragefile_idx,
        svc_acct_key=svc_acct_key,
        medianfile=medianfile,
        famfile=famfile,
        samples=samples,
        gt_cutoffs=RD_depth_sepcutoff,
        n_bins=n_RdTest_bins,
        prefix=basename(gt5kb_bed)
    }

    call make_subset_vcf as make_subset_vcf_gt5kb {
      input:
        vcf=cohort_vcf,
        bed=gt5kb_bed
    }

    call integrate_depth_GQ as integrate_GQ_gt5kb {
      input:
        vcf=make_subset_vcf_gt5kb.subset_vcf,
        RD_melted_genotypes=RD_genotype_gt5kb.melted_genotypes,
        RD_vargq=RD_genotype_gt5kb.varGQ,
    }

    call add_genotypes as add_genotypes_gt5kb {
      input:
        vcf=make_subset_vcf_gt5kb.subset_vcf,
        genotypes=integrate_GQ_gt5kb.genotypes,
        varGQ=integrate_GQ_gt5kb.varGQ,
        prefix=basename(gt5kb_bed, ".bed")
    }
  }

  call concat_genotyped_vcfs {
    input:
      lt5kb_vcfs=add_genotypes_lt5kb.genotyped_vcf,
      gt5kb_vcfs=add_genotypes_gt5kb.genotyped_vcf,
      batch=batch
  }

  output {
    File genotyped_vcf = concat_genotyped_vcfs.genotyped_vcf
  }
}

task split_variants {
  File vcf
  Int n_per_split

  command <<<
    svtk vcf2bed ${vcf} stdout \
      | awk -v OFS="\t" '(($5=="DEL" || $5=="DUP") && $3-$2>=5000) {print $1, $2, $3, $4, $6, $5}' \
      | split -l ${n_per_split} -a 6 - gt5kb.;
    svtk vcf2bed ${vcf} stdout \
      | awk -v OFS="\t" '(($5=="DEL" || $5=="DUP") && $3-$2<5000) {print $1, $2, $3, $4, $6, $5}' \
      | split -l ${n_per_split} -a 6 - lt5kb.;
  >>>

  output {
    Array[File] lt5kb_beds = glob("lt5kb.*")
    Array[File] gt5kb_beds = glob("gt5kb.*")
  }
    
  runtime {
      preemptible: 3
      docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
  }
}

task make_subset_vcf {
  File vcf
  File bed
  String prefix = basename(bed, ".bed")

  command <<<
    zcat ${vcf} | fgrep -e "#" > ${prefix}.vcf;
    zcat ${vcf} | fgrep -w -f <(cut -f4 ${bed}) >> ${prefix}.vcf;
    bgzip ${prefix}.vcf
  >>>

  output {
    File subset_vcf = "${prefix}.vcf.gz"
  }
    
  runtime {
      preemptible: 3
      docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
  }
}

task RdTest_genotype {
  File bed
  String coveragefile
  File medianfile
  File svc_acct_key
  File coveragefile_idx
  File famfile
  Array[String] samples
  File gt_cutoffs
  Int n_bins
  String prefix

  command <<<
    /opt/RdTest/localize_bincov.sh ${bed} ${coveragefile} ${coveragefile_idx} ${svc_acct_key};
    Rscript /opt/RdTest/RdTest.R \
      -b ${bed} \
      -c local_coverage.bed.gz \
      -m ${medianfile} \
      -f ${famfile} \
      -n ${prefix} \
      -w ${write_tsv(samples)} \
      -i ${n_bins} \
      -r ${gt_cutoffs} \
      -y /opt/RdTest/bin_exclude.bed.gz \
      -g TRUE;
    /opt/sv-pipeline/04_variant_resolution/scripts/merge_RdTest_genotypes.py ${prefix}.geno ${prefix}.gq rd.geno.cnv.bed;
    sort -k1,1V -k2,2n rd.geno.cnv.bed | uniq | bgzip -c > rd.geno.cnv.bed.gz
  >>>

  output {
    File genotypes = "${prefix}.geno"
    File copy_states = "${prefix}.median_geno"
    File metrics = "${prefix}.metrics"
    File gq = "${prefix}.gq"
    File varGQ = "${prefix}.vargq"
    File melted_genotypes = "rd.geno.cnv.bed.gz"
  }

  runtime {
      preemptible: 3
      docker: "talkowski/sv-pipeline-rdtest@sha256:764635fce650adac449b013058388a55653e8c7e6c075452a80f6e2a104754cd"
      disks: "local-disk 40 SSD"
  }
}

task integrate_depth_GQ {
  File vcf
  File RD_melted_genotypes
  File RD_vargq
  
  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/IntegrateGQ_depthonly.sh \
      ${vcf} \
      ${RD_melted_genotypes} \
      ${RD_vargq}
  >>>
  
  output {
    File genotypes = "genotype.indiv.depth.txt.gz"
    File varGQ = "genotype.variant.depth.txt.gz"
  }
    
  runtime {
    docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
    preemptible: 3
    memory: "16 GB"
    disks: "local-disk 200 SSD"
  }
}

task add_genotypes {
  File vcf
  File genotypes
  File varGQ
  String prefix
  
  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/add_genotypes.py \
        ${vcf} \
        ${genotypes} \
        ${varGQ} \
        ${prefix}.genotyped.vcf;
    vcf-sort -c ${prefix}.genotyped.vcf | bgzip -c > ${prefix}.genotyped.vcf.gz
  >>>
  
  output {
      File genotyped_vcf = "${prefix}.genotyped.vcf.gz"
  }
  
  runtime {
    docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
    preemptible: 3
    memory: "16 GB"
  }
}

task concat_genotyped_vcfs {
  Array[File] lt5kb_vcfs
  Array[File] gt5kb_vcfs
  String batch

  command <<<
    vcf-concat ${sep=' ' lt5kb_vcfs} ${sep=' ' gt5kb_vcfs} \
      | vcf-sort -c \
      | bgzip -c > ${batch}.depth.vcf.gz
  >>>

  output {
    File genotyped_vcf = "${batch}.depth.vcf.gz"
  }
  
  runtime {
    docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
    preemptible: 3
    memory: "16 GB"
  }
}