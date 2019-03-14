workflow RD_genotype_train {
  File vcf                # VCF to genotype
  String coveragefile     # batch coverage file
  File coveragefile_idx
  File medianfile         # batch median file
  File famfile            # batch famfile
  File svc_acct_key
  File rf_cutoffs         # Random forest cutoffs
  File seed_cutoffs
  Array[String] samples   # List of samples in batch
  String prefix           # batch/algorithm ID to use in output files
  Int n_bins              # number of RdTest bins
  Int n_per_split         # number of variants per RdTest split
  String reference_build  #hg19 or hg38

  call make_training_bed {
    input:
      sample_ID=samples[0],
      reference_build=reference_build
  }
 
  call RdTest_genotype as genotype_train {
    input:
      bed=make_training_bed.bed,
      coveragefile=coveragefile,
      medianfile=medianfile,
      coveragefile_idx=coveragefile_idx,
      svc_acct_key=svc_acct_key,
      famfile=famfile,
      samples=samples,
      gt_cutoffs=seed_cutoffs,
      n_bins=n_bins,
      prefix="train"
  }
 
  call generate_cutoff {
    input:
      copy_states=genotype_train.copy_states,
      max_copystate=4,
      prefix="train"
  }

  call update_cutoff {
    input:
      rf_cutoffs=rf_cutoffs,
      gt_cutoffs=generate_cutoff.cutoffs
  }

  call split_variants {
    input:
      vcf=vcf,
      n_per_split=n_per_split
  }

  scatter (pesr_bed in split_variants.pesr_beds) {
    call RdTest_genotype as genotype_pesr {
      input:
        bed=pesr_bed,
        coveragefile=coveragefile,
        coveragefile_idx=coveragefile_idx,
        svc_acct_key=svc_acct_key,
        medianfile=medianfile,
        famfile=famfile,
        samples=samples,
        gt_cutoffs=update_cutoff.pesr_sepcutoff,
        n_bins=n_bins,
        prefix=basename(pesr_bed)
    }
  }
  
  scatter (gt5kb_bed in split_variants.gt5kb_beds) {
    call RdTest_genotype as genotype_gt5kb {
      input:
        bed=gt5kb_bed,
        coveragefile=coveragefile,
        coveragefile_idx=coveragefile_idx,
        svc_acct_key=svc_acct_key,
        medianfile=medianfile,
        famfile=famfile,
        samples=samples,
        gt_cutoffs=update_cutoff.depth_sepcutoff,
        n_bins=n_bins,
        prefix=basename(gt5kb_bed)
    }
  }

  call merge_genotype_results {
    input:
      pesr_genotypes=genotype_pesr.genotypes,
      gt5kb_genotypes=genotype_gt5kb.genotypes,
      pesr_GQ=genotype_pesr.gq,
      gt5kb_GQ=genotype_gt5kb.gq,
      pesr_varGQ=genotype_pesr.varGQ,
      gt5kb_varGQ=genotype_gt5kb.varGQ
  }

  output {
    File genotypes = merge_genotype_results.genotypes
    File melted_genotypes = merge_genotype_results.melted_genotypes
    File GQ = merge_genotype_results.GQ
    File varGQ = merge_genotype_results.varGQ
    File pesr_sepcutoff = update_cutoff.pesr_sepcutoff
    File depth_sepcutoff = update_cutoff.depth_sepcutoff
  }
}

task make_training_bed {
  String sample_ID
  String reference_build

  command <<<
    if [ ${reference_build} == "hg19" ]; then
      awk -v OFS="\t" -v sample="${sample_ID}" '{$5=sample; print $1, $2, $3, $4, $5, $6}' /opt/RdTest/1kg.train.loci.bed > train.bed
    else
      awk -v OFS="\t" -v sample="${sample_ID}" '{$5=sample; print $1, $2, $3, $4, $5, $6}' /opt/RdTest/train_hg38_reviewed_final.bed > train.bed
    fi
  >>>

  output {
    File bed = "train.bed"
  }

  runtime {
      preemptible: 3
      docker: "talkowski/sv-pipeline-rdtest@sha256:764635fce650adac449b013058388a55653e8c7e6c075452a80f6e2a104754cd"
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
  >>>

  output {
    File genotypes = "${prefix}.geno"
    File copy_states = "${prefix}.median_geno"
    File metrics = "${prefix}.metrics"
    File gq = "${prefix}.gq"
    File varGQ = "${prefix}.vargq"
  }

  runtime {
      preemptible: 3
      docker: "talkowski/sv-pipeline-rdtest@sha256:764635fce650adac449b013058388a55653e8c7e6c075452a80f6e2a104754cd"
      disks: "local-disk 40 SSD"
  }
}

task generate_cutoff {
  File copy_states
  Int max_copystate
  String prefix

  command {
    Rscript /opt/RdTest/generate_cutoff.R ${copy_states} ${max_copystate} ${prefix}.cutoffs
  }

  output {
    File cutoffs = "${prefix}.cutoffs"
  }

  runtime {
      preemptible: 3
      docker: "talkowski/sv-pipeline-rdtest@sha256:764635fce650adac449b013058388a55653e8c7e6c075452a80f6e2a104754cd"
  }
}

task update_cutoff {
  File rf_cutoffs
  File gt_cutoffs 

  command <<<
    sep=$( awk -F'\t' '{if ($1=="PESR" && $6==1000 && $5=="RD_Median_Separation") print $2}' ${rf_cutoffs});
    cat ${gt_cutoffs} | awk -v var=$sep '{if ($1=="1" && $4>1-var) $4=1-var; else if ($1=="2" && $4<1+var) $4=1+var; print}' | tr ' ' '\t' > pesr_sepcutoff.txt;
    sep=$( awk -F'\t' '{if ($1=="Depth" && $5=="RD_Median_Separation") print $2}' ${rf_cutoffs} | sort -nr | head -n 1);
    cat ${gt_cutoffs} | awk -v var=$sep '{if ($1=="1" && $4>1-var) $4=1-var; else if ($1=="2" && $4<1+var) $4=1+var; print}'|tr ' ' '\t' > depth_sepcutoff.txt;
  >>> 
  
  output {
    File pesr_sepcutoff = "pesr_sepcutoff.txt"
    File depth_sepcutoff = "depth_sepcutoff.txt"
  }

  runtime {
      preemptible: 3
      docker: "talkowski/sv-pipeline-rdtest@sha256:764635fce650adac449b013058388a55653e8c7e6c075452a80f6e2a104754cd"
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
      | split -l ${n_per_split} -a 6 - pesr.;
  >>>

  output {
    Array[File] pesr_beds = glob("pesr.*")
    Array[File] gt5kb_beds = glob("gt5kb.*")
  }
    
  runtime {
      preemptible: 3
      docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
  }
}

task merge_genotype_results {
  Array[File] pesr_genotypes
  Array[File] gt5kb_genotypes
  Array[File] pesr_GQ
  Array[File] gt5kb_GQ
  Array[File] pesr_varGQ
  Array[File] gt5kb_varGQ

  command <<<
    cat ${sep=' ' pesr_genotypes} ${sep=' ' gt5kb_genotypes} | awk '!_[$0]++' > rd.geno.all;
    cat ${sep=' ' pesr_GQ} ${sep=' ' gt5kb_GQ} | awk '!_[$0]++' > rd.GQ.all;
    cat ${sep=' ' pesr_varGQ} ${sep=' ' gt5kb_varGQ} | awk '!_[$0]++' > rd.varGQ.all;
    
    /opt/sv-pipeline/04_variant_resolution/scripts/merge_RdTest_genotypes.py rd.geno.all rd.GQ.all rd.geno.cnv.bed;
    sort -k1,1V -k2,2n rd.geno.cnv.bed | uniq | bgzip -c > rd.geno.cnv.bed.gz
  >>>

  output {
    File genotypes = "rd.geno.all"
    File GQ = "rd.GQ.all"
    File melted_genotypes = "rd.geno.cnv.bed.gz"
    File varGQ = "rd.varGQ.all"
  }
    
  runtime {
      preemptible: 0
      docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
      memory: "32 GB"
      disks: "local-disk 50 SSD"
  }
}