workflow genotype_pesr_part2 {
  File cohort_vcf
  File RD_pesr_sepcutoff
  File RD_depth_sepcutoff
  File PE_metrics
  File SR_metrics
  Int n_per_split
  Int n_RdTest_bins
  String batch

  File medianfile
  File famfile
  File svc_acct_key
  Array[String] samples

  String coveragefile
  File coveragefile_idx
  
  String discfile
  File discfile_idx

  String splitfile
  File splitfile_idx

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

    call count_pe as count_pe_lt5kb {
      input:
        vcf=make_subset_vcf_lt5kb.subset_vcf,
        discfile=discfile,
        discfile_idx=discfile_idx,
        medianfile=medianfile,
        svc_acct_key=svc_acct_key,
        samples=samples
    }

    call genotype_PE_part2 as genotype_PE_part2_lt5kb {
      input:
        PE_counts=count_pe_lt5kb.pe_counts,
        PE_metrics=PE_metrics
    }

    call count_sr as count_sr_lt5kb {
      input:
        vcf=make_subset_vcf_lt5kb.subset_vcf,
        splitfile=splitfile,
        splitfile_idx=splitfile_idx,
        medianfile=medianfile,
        svc_acct_key=svc_acct_key,
        samples=samples
    }

    call genotype_SR_part2 as genotype_SR_part2_lt5kb {
      input:
        vcf=make_subset_vcf_lt5kb.subset_vcf,
        SR_counts=count_sr_lt5kb.sr_counts,
        SR_sum=count_sr_lt5kb.sr_sum,
        SR_metrics=SR_metrics
    }

    call integrate_GQ as integrate_GQ_lt5kb {
      input:
        vcf=make_subset_vcf_lt5kb.subset_vcf,
        RD_melted_genotypes=RD_genotype_lt5kb.melted_genotypes,
        RD_vargq=RD_genotype_lt5kb.varGQ,
        PE_genotypes=genotype_PE_part2_lt5kb.genotypes,
        PE_vargq=genotype_PE_part2_lt5kb.varGQ,
        SR_genotypes=genotype_SR_part2_lt5kb.genotypes,
        SR_vargq=genotype_SR_part2_lt5kb.varGQ
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

    call count_pe as count_pe_gt5kb {
      input:
        vcf=make_subset_vcf_gt5kb.subset_vcf,
        discfile=discfile,
        discfile_idx=discfile_idx,
        medianfile=medianfile,
        svc_acct_key=svc_acct_key,
        samples=samples
    }

    call genotype_PE_part2 as genotype_PE_part2_gt5kb {
      input:
        PE_counts=count_pe_gt5kb.pe_counts,
        PE_metrics=PE_metrics
    }

    call count_sr as count_sr_gt5kb {
      input:
        vcf=make_subset_vcf_gt5kb.subset_vcf,
        splitfile=splitfile,
        splitfile_idx=splitfile_idx,
        medianfile=medianfile,
        svc_acct_key=svc_acct_key,
        samples=samples
    }

    call genotype_SR_part2 as genotype_SR_part2_gt5kb {
      input:
        vcf=make_subset_vcf_gt5kb.subset_vcf,
        SR_counts=count_sr_gt5kb.sr_counts,
        SR_sum=count_sr_gt5kb.sr_sum,
        SR_metrics=SR_metrics
    }

    call integrate_GQ as integrate_GQ_gt5kb {
      input:
        vcf=make_subset_vcf_gt5kb.subset_vcf,
        RD_melted_genotypes=RD_genotype_gt5kb.melted_genotypes,
        RD_vargq=RD_genotype_gt5kb.varGQ,
        PE_genotypes=genotype_PE_part2_gt5kb.genotypes,
        PE_vargq=genotype_PE_part2_gt5kb.varGQ,
        SR_genotypes=genotype_SR_part2_gt5kb.genotypes,
        SR_vargq=genotype_SR_part2_gt5kb.varGQ
    }

    call add_genotypes as add_genotypes_gt5kb {
      input:
        vcf=make_subset_vcf_gt5kb.subset_vcf,
        genotypes=integrate_GQ_gt5kb.genotypes,
        varGQ=integrate_GQ_gt5kb.varGQ,
        prefix=basename(gt5kb_bed, ".bed")
    }
  }

  scatter (bca_bed in split_variants.bca_beds) {
    call make_subset_vcf as make_subset_vcf_bca {
      input:
        vcf=cohort_vcf,
        bed=bca_bed
    }

    call count_pe as count_pe_bca {
      input:
        vcf=make_subset_vcf_bca.subset_vcf,
        discfile=discfile,
        discfile_idx=discfile_idx,
        medianfile=medianfile,
        svc_acct_key=svc_acct_key,
        samples=samples
    }

    call genotype_PE_part2 as genotype_PE_part2_bca {
      input:
        PE_counts=count_pe_bca.pe_counts,
        PE_metrics=PE_metrics
    }

    call count_sr as count_sr_bca {
      input:
        vcf=make_subset_vcf_bca.subset_vcf,
        splitfile=splitfile,
        splitfile_idx=splitfile_idx,
        medianfile=medianfile,
        svc_acct_key=svc_acct_key,
        samples=samples
    }

    call genotype_SR_part2 as genotype_SR_part2_bca {
      input:
        vcf=make_subset_vcf_bca.subset_vcf,
        SR_counts=count_sr_bca.sr_counts,
        SR_sum=count_sr_bca.sr_sum,
        SR_metrics=SR_metrics
    }

    call integrate_pesr_GQ as integrate_GQ_bca {
      input:
        vcf=make_subset_vcf_bca.subset_vcf,
        PE_genotypes=genotype_PE_part2_bca.genotypes,
        PE_vargq=genotype_PE_part2_bca.varGQ,
        SR_genotypes=genotype_SR_part2_bca.genotypes,
        SR_vargq=genotype_SR_part2_bca.varGQ
    }

    call add_genotypes as add_genotypes_bca {
      input:
        vcf=make_subset_vcf_bca.subset_vcf,
        genotypes=integrate_GQ_bca.genotypes,
        varGQ=integrate_GQ_bca.varGQ,
        prefix=basename(bca_bed, ".bed")
    }
  }

  call concat_genotyped_vcfs {
    input:
      lt5kb_vcfs=add_genotypes_lt5kb.genotyped_vcf,
      gt5kb_vcfs=add_genotypes_gt5kb.genotyped_vcf,
      bca_vcfs=add_genotypes_bca.genotyped_vcf,
      batch=batch
  }

  call triple_stream_cat as cat_background_fail {
    input:
      files_a=genotype_SR_part2_lt5kb.background_fail,
      files_b=genotype_SR_part2_gt5kb.background_fail,
      files_c=genotype_SR_part2_bca.background_fail,
      outfile="${batch}.genotype_SR_part2_background_fail.txt"
  }

  call triple_stream_cat as cat_bothside_pass {
    input:
      files_a=genotype_SR_part2_lt5kb.bothside_pass,
      files_b=genotype_SR_part2_gt5kb.bothside_pass,
      files_c=genotype_SR_part2_bca.bothside_pass,
      outfile="${batch}.genotype_SR_part2_bothside_pass.txt"
  }

  # call merge_pe_counts {
  #   input:
  #     lt5kb_files=count_pe_lt5kb.pe_counts,
  #     gt5kb_files=count_pe_gt5kb.pe_counts,
  #     bca_files=count_pe_bca.pe_counts,
  #     batch=batch
  # }

  # call merge_sr_counts {
  #   input:
  #     lt5kb_files=count_sr_lt5kb.sr_counts,
  #     gt5kb_files=count_sr_gt5kb.sr_counts,
  #     bca_files=count_sr_bca.sr_counts,
  #     batch=batch
  # }

  output {
    File genotyped_vcf = concat_genotyped_vcfs.genotyped_vcf
    File background_fail = cat_background_fail.merged_file
    File bothside_pass = cat_bothside_pass.merged_file
    # File pe_counts = merge_pe_counts.merged_pe_counts
    # File sr_counts = merge_sr_counts.merged_sr_counts
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
    svtk vcf2bed ${vcf} stdout \
      | awk -v OFS="\t" '($5!="DEL" && $5!="DUP") {print $1, $2, $3, $4, $6, $5}' \
      | split -l ${n_per_split} -a 6 - bca.;
  >>>

  output {
    Array[File] lt5kb_beds = glob("lt5kb.*")
    Array[File] gt5kb_beds = glob("gt5kb.*")
    Array[File] bca_beds = glob("bca.*")
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
    sort -k1,1V -k2,2n rd.geno.cnv.bed | uniq | bgzip -c > rd.geno.cnv.bed.gz;
    echo "test"
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

task count_pe {
  File vcf
  String discfile
  File discfile_idx
  File medianfile
  File svc_acct_key
  Array[String] samples

  String prefix = basename(vcf, ".vcf")

  command <<<
    url=$(gsutil signurl -d 24h ${svc_acct_key} ${discfile} | sed '1d' | cut -f 4);
    svtk vcf2bed --split-bnd --no-header ${vcf} test.bed;
    awk -v OFS="\t" -v window=5000 '{if ($2-window>0){print $1,$2-window,$2+window}else{print $1,0,$2+window}}' test.bed  >> region.bed;
    awk -v OFS="\t" -v window=5000 '{if ($3-window>0){print $1,$3-window,$3+window}else{print $1,0,$3+window}}' test.bed  >> region.bed;
    sort -k1,1 -k2,2n region.bed > region.sorted.bed;
    bedtools merge -i region.sorted.bed > region.merged.bed;
    svtk remote_tabix "$url" ${discfile_idx} -R region.merged.bed | bgzip -c > PE.txt.gz;
    tabix -b 2 -e 2 PE.txt.gz;
    svtk count-pe --index PE.txt.gz.tbi -s ${write_tsv(samples)} --medianfile ${medianfile} ${vcf} PE.txt.gz ${prefix}.pe_counts.txt;
    gzip ${prefix}.pe_counts.txt
  >>>

  output {
    File pe_counts = "${prefix}.pe_counts.txt.gz"
  }

  runtime {
      preemptible: 3
      docker: "talkowski/sv-pipeline-remote-pysam@sha256:41a84644c1f7d339813c1176fdd6d42ed1ac770e430b053975d47da6e99f5f26"
  }
}

task genotype_PE_part2 {
  File PE_counts
  File PE_metrics
  
  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/PE_genotype.opt_part2.sh \
      ${PE_counts} \
      ${PE_metrics}
  >>>

  output {
    File genotypes = "pe.geno.withquality.txt.gz"
    File varGQ = "pe.variant.quality.final.txt.gz"
  }

  runtime {
      preemptible: 3
      docker: "talkowski/sv-pipeline-rdtest@sha256:764635fce650adac449b013058388a55653e8c7e6c075452a80f6e2a104754cd"
      disks: "local-disk 50 SSD"
  }
}

task count_sr {
  File vcf
  String splitfile
  File splitfile_idx
  File medianfile
  File svc_acct_key
  Array[String] samples

  String prefix = basename(vcf, ".vcf")

  command <<<
    url=$(gsutil signurl -d 24h ${svc_acct_key} ${splitfile} | sed '1d' | cut -f 4);
    svtk vcf2bed --split-bnd --no-header ${vcf} test.bed;
    awk -v OFS="\t" '{if ($2-250>0){print $1,$2-250,$2+250}else{print $1,0,$2+250}}' test.bed  >> region.bed;
    awk -v OFS="\t" '{if ($3-250>0){print $1,$3-250,$3+250}else{print $1,0,$3+250}}' test.bed  >> region.bed;
    sort -k1,1 -k2,2n region.bed > region.sorted.bed;
    bedtools merge -i region.sorted.bed > region.merged.bed;
    svtk remote_tabix "$url" ${splitfile_idx} -R region.merged.bed | bgzip -c > SR.txt.gz;
    tabix -b 2 -e 2 SR.txt.gz;
    svtk count-sr --index SR.txt.gz.tbi -s ${write_tsv(samples)} --medianfile ${medianfile} ${vcf} SR.txt.gz ${prefix}.sr_counts.txt;
    /opt/sv-pipeline/04_variant_resolution/scripts/sum_SR.sh ${prefix}.sr_counts.txt ${prefix}.sr_sum.txt.gz;
    gzip ${prefix}.sr_counts.txt
  >>>

  output {
    File sr_counts = "${prefix}.sr_counts.txt.gz"
    File sr_sum = "${prefix}.sr_sum.txt.gz"
  }

  runtime {
      preemptible: 3
      docker: "talkowski/sv-pipeline-remote-pysam@sha256:41a84644c1f7d339813c1176fdd6d42ed1ac770e430b053975d47da6e99f5f26"
  }
}

task genotype_SR_part2 {
  File vcf
  File SR_counts
  File SR_sum
  File SR_metrics
  
  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/SR_genotype.opt_part2.sh \
      ${vcf} \
      ${SR_counts} \
      ${SR_sum} \
      ${SR_metrics}
  >>>

  output {
    File genotypes = "sr.geno.withquality.txt.gz"
    File varGQ = "sr.variant.quality.final.txt.gz"
    File background_fail = "background.variant.fail.txt"
    File bothside_pass = "both.pass.txt"
  }

  runtime {
      preemptible: 3
      docker: "talkowski/sv-pipeline-rdtest@sha256:764635fce650adac449b013058388a55653e8c7e6c075452a80f6e2a104754cd"
      disks: "local-disk 60 SSD"
      memory: "16 GB"
  }
}

task integrate_GQ {
  File vcf
  File RD_melted_genotypes
  File RD_vargq
  File PE_genotypes
  File PE_vargq
  File SR_genotypes
  File SR_vargq
  
  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/IntegrateGQ.sh \
      ${vcf} \
      ${RD_melted_genotypes} \
      ${RD_vargq} \
      ${PE_genotypes} \
      ${PE_vargq} \
      ${SR_genotypes} \
      ${SR_vargq}
  >>>
  
  output {
    File genotypes = "genotype.indiv.txt.gz"
    File varGQ = "genotype.variant.txt.gz"
  }
    
  runtime {
    docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
    preemptible: 3
    memory: "16 GB"
    disks: "local-disk 200 SSD"
  }
}

task integrate_pesr_GQ {
  File vcf
  File PE_genotypes
  File PE_vargq
  File SR_genotypes
  File SR_vargq
  
  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/IntegrateGQ_PESR.sh \
      ${vcf} \
      ${PE_genotypes} \
      ${PE_vargq} \
      ${SR_genotypes} \
      ${SR_vargq};
    echo "test"
  >>>
  
  output {
    File genotypes = "genotype.indiv.txt.gz"
    File varGQ = "genotype.variant.txt.gz"
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
  Array[File] bca_vcfs
  String batch

  command <<<
    vcf-concat ${sep=' ' lt5kb_vcfs} ${sep=' ' gt5kb_vcfs} ${sep=' ' bca_vcfs} \
      | vcf-sort -c \
      | bgzip -c > ${batch}.pesr.vcf.gz
  >>>

  output {
    File genotyped_vcf = "${batch}.pesr.vcf.gz"
  }
  
  runtime {
    docker: "talkowski/sv-pipeline@sha256:239467dd98790097ed34637bdb5e32ffb6bc617bc295035792f037e12c5e5e7c"
    preemptible: 3
    memory: "32 GB"
    disks: "local-disk 200 HDD"
  }
}

task triple_stream_cat {
  Array[File] files_a
  Array[File] files_b
  Array[File] files_c
  String outfile

  command <<<
    cat ${sep=' ' files_a} ${sep=' ' files_b} ${sep=' ' files_c} | sort -Vk1,1 | uniq > ${outfile}
  >>>

  output {
    File merged_file = "${outfile}"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
    preemptible: 3
  }
}

task merge_pe_counts {
  Array[File] lt5kb_files
  Array[File] gt5kb_files
  Array[File] bca_files
  String batch

  command <<<
    echo -e "name\tsample\tcount" > ${batch}.04a_merged_pe_counts.txt
    while read file; do
      zcat $file | sed '1d' | awk -v OFS="\t" '{ if ($NF>0) print $0 }'
    done < <( cat ${sep='\n' lt5kb_files} ${sep='\n' gt5kb_files} ${sep='\n' bca_files} ) >> ${batch}.04a_merged_pe_counts.txt
    gzip -f ${batch}.04a_merged_pe_counts.txt
  >>>

  output {
    File merged_pe_counts = "${batch}.04a_merged_pe_counts.txt.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
    preemptible: 3
    disks: "local-disk 200 HDD"
  }
}

task merge_sr_counts {
  Array[File] lt5kb_files
  Array[File] gt5kb_files
  Array[File] bca_files
  String batch

  command <<<
    echo -e "name\tcoord\tsample\tcount" > ${batch}.04a_merged_sr_counts.txt
    while read file; do
      zcat $file | sed '1d' | awk -v OFS="\t" '{ if ($NF>0) print $0 }'
    done < <( cat ${sep='\n' lt5kb_files} ${sep='\n' gt5kb_files} ${sep='\n' bca_files} ) >> ${batch}.04a_merged_sr_counts.txt
    gzip -f ${batch}.04a_merged_sr_counts.txt
  >>>

  output {
    File merged_sr_counts = "${batch}.04a_merged_sr_counts.txt.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
    preemptible: 3
    disks: "local-disk 200 HDD"
  }
}

