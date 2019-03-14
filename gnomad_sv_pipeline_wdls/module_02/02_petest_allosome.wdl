# Run petest on a single allosome, parallelizing across a fixed split size
workflow petest_allosome {
  File vcf                      # Input VCF
  String discfile               # Discordant pair file
  File medianfile             # Median file
  File discfile_idx             # Tabix index of discordant pair file
  File famfile                  # Batch fam file (used to segregate by sex)
  File svc_acct_key             # Service account key
  String batch                  # Batch ID
  String algorithm              # Algorithm ID
  String chrom                  # Chromosome being processed
  Int split_size                # Number of lines in each petest split

  # Compute the length of the suffix needed to accomodate all splits
  call compute_suffix_len {
    input:
      vcf=vcf,
      chrom=chrom,
      split_size=split_size
  }

  # Make list of males for petest whitelisting
  call choose_sex as choose_males {
    input:
      famfile=famfile,
      sex="1"
  }

  # Make list of females for petest whitelisting
  call choose_sex as choose_females {
    input:
      famfile=famfile,
      sex="2"
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

  # Run petest on each split
  scatter (split in split_vcf.split_vcfs) {
    # Add VCF header to split
    call reheader_split {
      input:
        vcf=vcf,
        split=split
    }

    # Run petest on males
    call petest as petest_males {
      input:
        vcf=reheader_split.split_w_header,
        prefix=basename(split),
        medianfile=medianfile,
        discfile=discfile,
        discfile_idx=discfile_idx,
        svc_acct_key=svc_acct_key,
        whitelist=choose_males.samples
    }
    
    # Run petest on females
    call petest as petest_females {
      input:
        vcf=reheader_split.split_w_header,
        prefix=basename(split),
        medianfile=medianfile,
        discfile=discfile,
        discfile_idx=discfile_idx,
        svc_acct_key=svc_acct_key,
        whitelist=choose_females.samples
    }

    # Combine male and female test results
    call merge_allosomes {
      input:
        male_petest=petest_males.stats,
        female_petest=petest_females.stats,
        chrom=chrom
    }
  }

  # Merge splits into single file
  call merge_splits {
    input:
      stats=merge_allosomes.merged_petest,
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

# Select all samples of a given sex from a fam file
task choose_sex {
  File famfile
  String sex

  command <<<
    awk -v sex=${sex} '($5==sex) {print $2}' ${famfile} > ${sex}.list
  >>>

  output {
    File samples = "${sex}.list"
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

  command {
    tabix -p vcf ${vcf};
    tabix ${vcf} ${chrom} | sort -R | split -a ${suffix_len} -d -l ${split_size} - ${batch}.${algorithm}.split.
  }

  output {
    Array[File] split_vcfs = glob("${batch}.${algorithm}.split.*")
  }
      
  runtime {
  	preemptible: 3
    docker: "talkowski/sv-pipeline-remote-pysam"
  }
}

# Restore VCF header to split files
task reheader_split {
  File vcf
  File split

  command {
    cat <(zcat ${vcf} | sed -n -e '/^#/p') ${split} | bgzip -c > ${basename(split)}.vcf.gz
  }

  output {
    File split_w_header = "${basename(split)}.vcf.gz"
  }
      
  runtime {
  	preemptible: 3
    docker: "talkowski/sv-pipeline-remote-pysam"
  }
}

# Run petest, whitelisting by sex
task petest {
  File vcf
  String discfile
  File medianfile
  File discfile_idx
  File svc_acct_key
  File whitelist
  String prefix
 
  command { 
    url=$(gsutil signurl -d 24h ${svc_acct_key} ${discfile} | sed '1d' | cut -f 4);
    echo $url;
    svtk pe-test -o 1000 --index ${discfile_idx} --medianfile ${medianfile} --samples ${whitelist} ${vcf} "$url" ${prefix}.stats
  }
  
  output {
    File stats = "${prefix}.stats"
  }
      
  runtime {
  	preemptible: 3
    docker: "talkowski/sv-pipeline-remote-pysam"
  }
}

# Combine male and female test results
# On chrX, use female test results unless the variant appears only in males
# On chrY, use male test results
task merge_allosomes {
  File male_petest
  File female_petest
  String chrom

  command <<<
    python3 <<CODE
    import pandas as pd
    males = pd.read_table("${male_petest}")
    females = pd.read_table("${female_petest}")
    if "${chrom}" == 'Y' or "${chrom}" == 'chrY':
      males.to_csv("${basename(male_petest)}.merged.csv", sep='\t', index=False, na_rep='NA')
    else:
      male_only = females.log_pval.isnull()
      females.loc[male_only] = males
      females.to_csv("${basename(male_petest)}.merged.csv", sep='\t', index=False, na_rep='NA')
    CODE
  >>>

  output {
    File merged_petest = "${basename(male_petest)}.merged.csv"
  }
      
  runtime {
  	preemptible: 3
    docker: "talkowski/sv-pipeline-remote-pysam"
  }
}

# Merge split petest results into single file
task merge_splits {
  Array[File] stats
  String prefix

  command <<<
    while read split; do
      sed -e '1d' $split;
    done < ${write_tsv(stats)} | cat <(head -n1 ${stats[0]}) - > ${prefix}.stats
  >>>

  output {
    File merged_stats = "${prefix}.stats"
  }
      
  runtime {
  	preemptible: 3
    docker: "talkowski/sv-pipeline-remote-pysam"
  }
}