workflow PE_genotype_train {
  File batch_vcf    # variants from just the batch in question
  String discfile
  Int n_per_split
  File medianfile
  File discfile_idx
  File svc_acct_key
  Array[String] samples
  String batch_ID
  File RF_cutoffs
  File RD_genotypes
  File RD_melted_genotypes
  File blacklist

  call vcf2bed as make_batch_bed {
    input:
      vcf=batch_vcf,
      prefix=batch_ID
  }
  
  call split_vcf as split_batch_vcf {
    input:
      vcf=batch_vcf,
      n_per_split=n_per_split
  }

  scatter (vcf in split_batch_vcf.vcfs) {
    call count_pe as count_batch_pe {
      input:
        vcf=vcf,
        discfile=discfile,
        discfile_idx=discfile_idx,
        medianfile=medianfile,
        svc_acct_key=svc_acct_key,
        samples=samples
    }
  }

  call merge_pe_counts {
    input:
      count_list=count_batch_pe.pe_counts
  }

  call genotype_PE_part1 {
    input:
      bed=make_batch_bed.bed,
      RF_cutoffs=RF_cutoffs,
      PE_counts=merge_pe_counts.counts,
      RD_genotypes=RD_genotypes,
      RD_melted_genotypes=RD_melted_genotypes,
      blacklist=blacklist
  }

  output {
    File PE_genotypes = genotype_PE_part1.genotypes
    File PE_varGQ = genotype_PE_part1.varGQ
    File PE_metrics = genotype_PE_part1.PE_metrics
    File PE_train = genotype_PE_part1.PE_train
  }
}

task vcf2bed {
  File vcf
  String prefix

  command {
    svtk vcf2bed ${vcf} -i ALGORITHMS ${prefix}.bed
  }

  output {
    File bed = "${prefix}.bed"
  }

  runtime {
      preemptible: 3
      docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
  }
}

task split_vcf {
  File vcf
  Int n_per_split

  command <<<
    if [[ ${vcf} == *.gz ]] ; then
      zcat ${vcf} | sed -n -e '/^#/p' > header.vcf;
      zcat ${vcf} | sed -e '/^#/d' | split -l ${n_per_split} - pe;
    else
      sed -n -e '/^#/p' ${vcf} > header.vcf;
      sed -e '/^#/d' ${vcf} | split -l ${n_per_split} - pe;
    fi
    for f in pe*; do cat header.vcf $f > $f.vcf; done
  >>>

  output {
    Array[File] vcfs = glob("pe*.vcf")
  }

  runtime {
      preemptible: 3
      docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
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

task merge_pe_counts {
  Array[File] count_list

  command {
    zcat ${sep=' ' count_list} | fgrep -v -e "name" | gzip -c > pe_counts.txt.gz
  }

  output {
    File counts = "pe_counts.txt.gz"
  }

  runtime {
      preemptible: 3
      docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
      disks: "local-disk 50 SSD"
  }
}

task genotype_PE_part1 {
  File bed
  File RF_cutoffs
  File PE_counts
  File RD_genotypes
  File RD_melted_genotypes
  File blacklist
  
  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/PE_genotype.sh \
      ${bed} \
      ${PE_counts} \
      ${RD_genotypes} \
      ${RD_melted_genotypes} \
      ${RF_cutoffs} \
      ${blacklist} \
      /opt/RdTest/generate_cutoff_PE.R 
  >>>

  output {
    File PE_train = "pe.train.include.txt"
    File PE_metrics = "pe_metric_file.txt"
    File genotypes = "pe.geno.withquality.txt.gz"
    File varGQ = "pe.variant.quality.final.txt.gz"
  }

  runtime {
      preemptible: 0
      docker: "talkowski/sv-pipeline-rdtest@sha256:764635fce650adac449b013058388a55653e8c7e6c075452a80f6e2a104754cd"
      disks: "local-disk 50 SSD"
  }
}