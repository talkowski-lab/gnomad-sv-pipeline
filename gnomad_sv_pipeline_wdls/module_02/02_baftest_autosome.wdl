# Run baftest on a single autosome, parallelizing across a fixed split size
workflow baftest_autosome {
  File vcf                      # Input VCF
  String baf_metrics            # Matrix of BAF statistics
  File baf_metrics_idx          # Tabix index of BAF matrix
  Array[String] samples			# list of samples in batch
  File svc_acct_key				# Service account json
  String batch                  # Batch ID
  String algorithm              # Algorithm ID
  String chrom                  # Chromosome being processed
  Int split_size                # Number of lines in each baftest split

  # Compute the length of the suffix needed to accomodate all splits
  call compute_suffix_len {
    input:
      vcf=vcf,
      chrom=chrom,
      split_size=split_size
  }

  # Split the VCF into smaller chunks
  call split_vcf {
    input:
      vcf=vcf,
      batch=batch,
      algorithm=algorithm,
      chrom=chrom,
      split_size=split_size,
      suffix_len=compute_suffix_len.len
  }

  # Run baftest on each split
  scatter (split in split_vcf.split_beds) {
    # Run baftest
    call baftest {
      input:
        bed=split,
        prefix=basename(split),
        baf_metrics=baf_metrics,
        baf_metrics_idx=baf_metrics_idx,
        samples=samples,
        batch=batch,
        svc_acct_key=svc_acct_key
    }
  }

  # Merge splits into single file
  call merge_splits {
    input:
      stats=baftest.stats,
      prefix="${batch}.${algorithm}.${chrom}"
  }

  output {
    File stats = merge_splits.merged_stats
  }
}

# Compute the length of the suffix necessary to accommodate all splits
task compute_suffix_len {
  File vcf
  String chrom
  Int split_size

  command <<<
    tabix -p vcf ${vcf};
    python3 <<CODE
    import numpy as np
    import pysam
    vcf = pysam.VariantFile('${vcf}')
    for i, record in enumerate(vcf.fetch('${chrom}')):
      continue
    n_records = i + 1
    n_splits = int(np.ceil(n_records / ${split_size}))
    suffix_len = max(int(np.ceil(np.log10(n_splits))), 1)
    print(suffix_len)
    CODE
  >>>

  output {
    Int len = read_int(stdout())
  }
  
  runtime {
  	preemptible: 3
  	docker: "talkowski/sv-pipeline-remote-pysam"
  }
}

# Split VCF into fixed size chunks
task split_vcf {
  File vcf
  String batch
  String algorithm
  String chrom
  
  Int split_size
  Int suffix_len

  command <<<
    tabix -p vcf ${vcf};
    tabix -h ${vcf} ${chrom} \
      | svtk vcf2bed --no-header stdin stdout \
      | fgrep -e "DEL" -e "DUP" \
      | awk -v OFS="\t" '{print $1, $2, $3, $4, $6, $5}' \
      | awk '($3-$2>=10000 && $3-$2<10000000)' \
      | split -a ${suffix_len} -d -l 300 - ${batch}.${algorithm}.split.gt10kb.
    tabix -h ${vcf} ${chrom} \
      | svtk vcf2bed --no-header stdin stdout \
      | fgrep -e "DEL" -e "DUP" \
      | awk -v OFS="\t" '{print $1, $2, $3, $4, $6, $5}' \
      | awk '($3-$2<10000)' \
      | sort -k1,1V -k2,2n \
      | split -a ${suffix_len} -d -l ${split_size} - ${batch}.${algorithm}.split.
  >>>

  output {
    Array[File] split_beds = glob("${batch}.${algorithm}.split.*")
  }
  
  runtime {
  	preemptible: 3
  	docker: "talkowski/sv-pipeline-remote-pysam"
  }
}

# Run baftest
task baftest {
  File bed
  String baf_metrics
  File baf_metrics_idx
  Array[String] samples
  File svc_acct_key
  String prefix
  String batch
 
  command <<<
    echo -e "sample\tgroup\tbatch" > batch.key;
    awk -v batch=${batch} -v OFS="\t" '{print $1, $1, batch}' ${write_tsv(samples)} >> batch.key;
  	url=$(gsutil signurl -d 24h ${svc_acct_key} ${baf_metrics} | sed '1d' | cut -f 4);
    start=$(cut -f2 ${bed} | sort -k1,1n | head -n1);
    end=$(cut -f3 ${bed} | sort -k1,1n | tail -n1);
    chrom=$(cut -f1 ${bed} | head -n1);
    svtk remote_tabix "$url" ${baf_metrics_idx} "$chrom":"$start"-"$end" | bgzip -c > local_baf.bed.gz;
    tabix -b2 local_baf.bed.gz;
    svtk baf-test ${bed} local_baf.bed.gz --batch batch.key > ${prefix}.metrics
  >>>
  
  output {
    File stats = "${prefix}.metrics"
  }
  
  runtime {
  	preemptible: 3
  	memory: "10 GB"
    disks: "local-disk 50 SSD"
  	docker: "talkowski/sv-pipeline-remote-pysam"
  }
}

# Merge split baftest results into single file
task merge_splits {
  Array[File] stats
  String prefix

  command <<<
    echo -n "chrom start end name samples svtype delstat snp_ratio " > ${prefix}.stats;
    echo -n "del_loglik dupstat KS_stat KS_pval total_case_snps " >> ${prefix}.stats;
    echo -n "total_snps n_nonROH_cases n_samples mean_control_snps " >> ${prefix}.stats;
    echo -n "n_nonROH_controls n_controls" >> ${prefix}.stats;
    sed -i -e 's/ /\t/g' ${prefix}.stats;
    while read split; do
      cat $split;
    done < ${write_tsv(stats)} >> ${prefix}.stats
  >>>

  output {
    File merged_stats = "${prefix}.stats"
  }
  
  runtime {
  	preemptible: 3
  	docker: "talkowski/sv-pipeline-remote-pysam"
  }
}