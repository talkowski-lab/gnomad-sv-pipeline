import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:02_rdtest_autosome/versions/12/plain-WDL/descriptor" as auto
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:02_rdtest_allosome/versions/7/plain-WDL/descriptor" as allo

# Parallelize rdtest on a single VCF across chromosomes
workflow rdtest_by_chrom {
  File vcf                      # Input VCF
  String coveragefile           # Bincov matrix
  File coveragefile_idx         # Tabix index of bincov matrix
  File medianfile               # Median coverage of each sample
  File famfile                  # Batch fam file 
  File autosome_contigs         # Autosomes .fai
  File allosome_contigs         # Allosomes .fai
  File svc_acct_key				# Service account json
  String batch                  # Batch ID
  String algorithm              # Algorithm ID
  Int split_size                # Number of lines in each rdtest split

  Array[Array[String]] autosomes = read_tsv(autosome_contigs) 
  Array[Array[String]] allosomes = read_tsv(allosome_contigs) 

  # Run rdtest on each autosome
  scatter (autosome in autosomes) {
    call auto.rdtest_autosome {
      input:
        vcf=vcf,
        coveragefile=coveragefile,
        coveragefile_idx=coveragefile_idx,
        medianfile=medianfile,
        famfile=famfile,
        batch=batch,
        algorithm=algorithm,
        chrom=autosome[0],
        split_size=split_size,
        svc_acct_key=svc_acct_key
    }
  }
  
  # Run rdtest on each allosome
  scatter (allosome in allosomes) {
    call allo.rdtest_allosome {
      input:
        vcf=vcf,
        coveragefile=coveragefile,
        coveragefile_idx=coveragefile_idx,
        medianfile=medianfile,
        famfile=famfile,
        batch=batch,
        algorithm=algorithm,
        chrom=allosome[0],
        split_size=split_size,
        svc_acct_key=svc_acct_key
    }
  }

  # Combine rdtest results into single file
  call merge_rdtest {
    input:
      autosomes=rdtest_autosome.stats,
      allosomes=rdtest_allosome.stats,
      prefix="${batch}.${algorithm}"
  }

  output {
    File rdtest = merge_rdtest.merged_stats
  }
}

# Combine per-chromosome rdtest results into single table
task merge_rdtest {
  Array[File] autosomes
  Array[File] allosomes
  String prefix

  command <<<
    cat ${write_tsv(autosomes)} ${write_tsv(allosomes)} > splits.list;
    while read split; do
      sed -e '1d' $split;
    done < splits.list | cat <(head -n1 ${autosomes[0]}) - > ${prefix}.stats
  >>>

  output {
    File merged_stats = "${prefix}.stats"
  }
  
  runtime {
  	preemptible: 3
  	docker: "talkowski/sv-pipeline-remote-pysam"
  }
}