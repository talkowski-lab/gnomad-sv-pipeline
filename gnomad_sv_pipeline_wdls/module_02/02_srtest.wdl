import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:02_srtest_autosome/versions/12/plain-WDL/descriptor" as auto
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:02_srtest_allosome/versions/11/plain-WDL/descriptor" as allo

# Parallelize srtest on a single VCF across chromosomes
workflow srtest_by_chrom {
  File vcf                      # Input VCF
  String splitfile              # Split read file
  String medianfile             # Medianfile
  File splitfile_idx            # Tabix index of split read file
  File famfile                  # Batch fam file 
  File autosome_contigs         # Autosomes .fai
  File allosome_contigs         # Allosomes .fai
  File svc_acct_key             # Service account json
  String batch                  # Batch ID
  String algorithm              # Algorithm ID
  Int split_size                # Number of lines in each srtest split

  Array[Array[String]] autosomes = read_tsv(autosome_contigs) 
  Array[Array[String]] allosomes = read_tsv(allosome_contigs) 

  # Run srtest on each autosome
  scatter (autosome in autosomes) {
    call auto.srtest_autosome {
      input:
        vcf=vcf,
        splitfile=splitfile,
        medianfile=medianfile,
        splitfile_idx=splitfile_idx,
        batch=batch,
        algorithm=algorithm,
        chrom=autosome[0],
        split_size=split_size,
        svc_acct_key=svc_acct_key
    }
  }
  
  # Run srtest on each allosome
  scatter (allosome in allosomes) {
    call allo.srtest_allosome {
      input:
        vcf=vcf,
        splitfile=splitfile,
        medianfile=medianfile,
        splitfile_idx=splitfile_idx,
        famfile=famfile,
        batch=batch,
        algorithm=algorithm,
        chrom=allosome[0],
        split_size=split_size,
        svc_acct_key=svc_acct_key
    }
  }

  # Combine srtest results into single file
  call merge_srtest {
    input:
      autosomes=srtest_autosome.stats,
      allosomes=srtest_allosome.stats,
      prefix="${batch}.${algorithm}"
  }

  output {
    File srtest = merge_srtest.merged_stats
  }
}

# Combine per-chromosome srtest results into single table
task merge_srtest {
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
    disks: "local-disk 100 SSD"
  }
}