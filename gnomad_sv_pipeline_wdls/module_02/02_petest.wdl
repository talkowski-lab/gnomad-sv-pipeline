import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:02_petest_autosome/versions/14/plain-WDL/descriptor" as auto
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:02_petest_allosome/versions/10/plain-WDL/descriptor" as allo

# Parallelize petest on a single VCF across chromosomes
workflow petest_by_chrom {
  File vcf                      # Input VCF
  String discfile                 # Discordant pair file
  String medianfile				# Medianfile
  File discfile_idx             # Tabix index of discordant pair file
  File famfile                  # Batch fam file 
  File autosome_contigs         # Autosomes .fai
  File allosome_contigs         # Allosomes .fai
  File svc_acct_key
  String batch                  # Batch ID
  String algorithm              # Algorithm ID
  Int split_size                # Number of lines in each petest split

  Array[Array[String]] autosomes = read_tsv(autosome_contigs) 
  Array[Array[String]] allosomes = read_tsv(allosome_contigs) 

  # Run petest on each autosome
  scatter (autosome in autosomes) {
    call auto.petest_autosome {
      input:
        vcf=vcf,
        discfile=discfile,
        medianfile=medianfile,
        discfile_idx=discfile_idx,
        batch=batch,
        algorithm=algorithm,
        chrom=autosome[0],
        split_size=split_size,
        svc_acct_key=svc_acct_key
    }
  }
  
  # Run petest on each allosome
  scatter (allosome in allosomes) {
    call allo.petest_allosome {
      input:
        vcf=vcf,
        discfile=discfile,
        medianfile=medianfile,
        discfile_idx=discfile_idx,
        famfile=famfile,
        batch=batch,
        algorithm=algorithm,
        chrom=allosome[0],
        split_size=split_size,
        svc_acct_key=svc_acct_key
    }
  }

  # Combine petest results into single file
  call merge_petest {
    input:
      autosomes=petest_autosome.stats,
      allosomes=petest_allosome.stats,
      prefix="${batch}.${algorithm}"
  }

  output {
    File petest = merge_petest.merged_stats
  }
}

# Combine per-chromosome petest results into single table
task merge_petest {
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
  	docker: "talkowski/sv-pipeline"
  }
}