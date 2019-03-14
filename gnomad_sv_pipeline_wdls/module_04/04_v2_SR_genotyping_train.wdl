workflow SR_genotype_train {
  File batch_vcf
  String splitfile
  Int n_per_split
  File medianfile
  File splitfile_idx
  File svc_acct_key
  Array[String] samples
  String batch_ID
  File RF_cutoffs
  File RD_melted_genotypes
  File PE_train
  File PE_genotypes

  call split_vcf as split_batch_vcf {
    input:
      vcf=batch_vcf,
      n_per_split=n_per_split
  }

  scatter (vcf in split_batch_vcf.vcfs) {
    call count_sr as count_batch_sr {
      input:
        vcf=vcf,
        splitfile=splitfile,
        splitfile_idx=splitfile_idx,
        medianfile=medianfile,
        svc_acct_key=svc_acct_key,
        samples=samples
    }
  }

  call merge_sr_counts {
    input:
      count_list=count_batch_sr.sr_counts,
      sum_list=count_batch_sr.sr_sum
  }

  call genotype_SR_part1 {
    input:
      vcf=batch_vcf,
      RF_cutoffs=RF_cutoffs,
      SR_counts=merge_sr_counts.counts,
      SR_sum=merge_sr_counts.sum,
      RD_melted_genotypes=RD_melted_genotypes,
      PE_train=PE_train,
      samples=samples,
      PE_genotypes=PE_genotypes
  }

  output {
    File SR_metrics = genotype_SR_part1.SR_metrics
  }
}

task split_vcf {
  File vcf
  Int n_per_split

  command <<<
      if [[ ${vcf} == *.gz ]] ; then
      echo "gzipped";
      zcat ${vcf} | sed -n -e '/^#/p' > header.vcf;
      zcat ${vcf} | sed -e '/^#/d' | split -l ${n_per_split} - sr;
    else
      echo "plaintext";
      sed -n -e '/^#/p' ${vcf} > header.vcf;
      sed -e '/^#/d' ${vcf} | split -l ${n_per_split} - sr;
    fi
    for f in sr*; do cat header.vcf $f | bgzip -c > $f.vcf.gz; done
  >>>

  output {
    Array[File] vcfs = glob("sr*.vcf.gz")
  }

  runtime {
      preemptible: 3
      docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
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

task merge_sr_counts {
  Array[File] count_list
  Array[File] sum_list

  command {
    zcat ${sep=' ' count_list} | fgrep -v -e "name" | gzip -c > sr_counts.txt.gz;
    cat ${sep=' ' sum_list} > sr_sum.txt.gz
  }

  output {
    File counts = "sr_counts.txt.gz"
    File sum = "sr_sum.txt.gz"
  }

  runtime {
      preemptible: 3
      docker: "talkowski/sv-pipeline@sha256:e5c7ce65c2e0c851261679b62095a13f42d0e4b4fef70b1d0183f2767e4ec53c"
      disks: "local-disk 60 SSD"
  }
}

task genotype_SR_part1 {
  File vcf
  File SR_counts
  File SR_sum
  File RD_melted_genotypes
  File RF_cutoffs
  Array[String] samples
  File PE_train
  File PE_genotypes
  
  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/SR_genotype.opt_part1.sh \
      ${vcf} \
      ${SR_counts} \
      ${SR_sum} \
      ${RD_melted_genotypes} \
      ${RF_cutoffs} \
      ${write_tsv(samples)} \
      ${PE_train} \
      ${PE_genotypes}
  >>>

  output {
    File SR_metrics = "sr_metric_file.txt"
  }

  runtime {
      preemptible: 0
      docker: "talkowski/sv-pipeline-rdtest@sha256:764635fce650adac449b013058388a55653e8c7e6c075452a80f6e2a104754cd"
      disks: "local-disk 60 SSD"
      memory: "16 GB"
  }
}