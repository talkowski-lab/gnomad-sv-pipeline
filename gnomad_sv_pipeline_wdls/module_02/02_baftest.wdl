import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:02_baftest_autosome/versions/12/plain-WDL/descriptor" as auto

# Parallelize baftest on a single VCF across chromosomes
workflow baftest_by_chrom {
  File vcf                      # Input VCF
  String baf_metrics            # Matrix of BAF statistics
  File baf_metrics_idx          # Tabix index of BAF matrix
  File autosome_contigs         # Autosomes .fai
  File svc_acct_key				# Service account json
  Array[String] samples			# List of samples in batch
  String batch                  # Batch ID
  String algorithm              # Algorithm ID
  Int split_size                # Number of lines in each baftest split

  Array[Array[String]] autosomes = read_tsv(autosome_contigs) 

  # Run baftest on each autosome
  scatter (autosome in autosomes) {
    call auto.baftest_autosome {
      input:
        vcf=vcf,
        baf_metrics=baf_metrics,
        baf_metrics_idx=baf_metrics_idx,
        batch=batch,
        algorithm=algorithm,
        chrom=autosome[0],
        split_size=split_size,
        samples=samples,
        svc_acct_key=svc_acct_key
    }
  }
  
  # Combine baftest results into single file
  call merge_baftest {
    input:
      autosomes=baftest_autosome.stats,
      prefix="${batch}.${algorithm}"
  }

  output {
    File baftest = merge_baftest.merged_stats
  }
}

# Combine per-chromosome baftest results into single table
task merge_baftest {
  Array[File] autosomes
  String prefix

  command <<<
    while read split; do
      sed -e '1d' $split;
    done < ${write_tsv(autosomes)} | cat <(head -n1 ${autosomes[0]}) - > ${prefix}.stats
  >>>

  output {
    File merged_stats = "${prefix}.stats"
  }
  
  runtime {
  	preemptible: 3
  	docker: "talkowski/sv-pipeline-remote-pysam"
  }
}