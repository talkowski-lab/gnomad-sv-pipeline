# Run rdtest on a single allosome, parallelizing across a fixed split size
workflow rdtest_allosome {
  File vcf                      # Input VCF
  String coveragefile             # Bincov matrix
  File coveragefile_idx         # Tabix index of bincov matrix
  File medianfile               # Median coverage of each sample
  File famfile                  # Batch fam file 
  File svc_acct_key				# Service account json
  String batch                  # Batch ID
  String algorithm              # Algorithm ID
  String chrom                  # Chromosome being processed
  Int split_size                # Number of lines in each rdtest split

  # Compute the length of the suffix needed to accomodate all splits
  call compute_suffix_len {
    input:
      vcf=vcf,
      chrom=chrom,
      split_size=split_size
  }

  # Make list of males for rdtest whitelisting
  call choose_sex as choose_males {
    input:
      famfile=famfile,
      sex="1"
  }

  # Make list of females for rdtest whitelisting
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

  # Run rdtest on each split
  scatter (split in split_vcf.split_beds) {
    # Run rdtest on males
    call rdtest as rdtest_males {
      input:
        bed=split,
        prefix=basename(split),
        coveragefile=coveragefile,
        coveragefile_idx=coveragefile_idx,
        medianfile=medianfile,
        famfile=famfile,
        whitelist=choose_males.samples,
        svc_acct_key=svc_acct_key
    }
    
    # Run rdtest on females
    call rdtest as rdtest_females {
      input:
        bed=split,
        prefix=basename(split),
        coveragefile=coveragefile,
        coveragefile_idx=coveragefile_idx,
        medianfile=medianfile,
        famfile=famfile,
        whitelist=choose_females.samples,
        svc_acct_key=svc_acct_key 
    }

    # Combine male and female test results
    call merge_allosomes {
      input:
        male_rdtest=rdtest_males.stats,
        female_rdtest=rdtest_females.stats,
        chrom=chrom
    }
  }

  # Merge splits into single file
  call merge_splits {
    input:
      stats=merge_allosomes.merged_rdtest,
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

  command <<<
    tabix -p vcf ${vcf};
    tabix -h ${vcf} ${chrom} \
      | svtk vcf2bed --no-header stdin stdout \
      | fgrep -e "DEL" -e "DUP" \
      | awk -v OFS="\t" '{print $1, $2, $3, $4, $6, $5}' \
      | awk '($3-$2>=10000)' \
      > ${batch}.${algorithm}.split.gt10kb;
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

# Run rdtest, whitelisting by sex
task rdtest {
  File bed
  String coveragefile
  File coveragefile_idx
  File medianfile
  File famfile
  File whitelist
  File svc_acct_key
  String prefix
 
  command <<<
  	url=$(gsutil signurl -d 24h ${svc_acct_key} ${coveragefile} | sed '1d' | cut -f 4);
    start=$(cut -f2 ${bed} | sort -k1,1n | head -n1);
    end=$(cut -f3 ${bed} | sort -k1,1n | tail -n1);
    chrom=$(cut -f1 ${bed} | head -n1);
    svtk remote_tabix --header "$url" ${coveragefile_idx} "$chrom":"$start"-"$end" |sed 's/Chr/chr/g'|sed 's/Start/start/g'|sed 's/End/end/' | bgzip -c > local_coverage.bed.gz;
    tabix -p bed local_coverage.bed.gz;
    Rscript /opt/RdTest/RdTest.R \
      -b ${bed} \
      -n ${prefix} \
      -c local_coverage.bed.gz \
      -m ${medianfile} \
      -f ${famfile} \
      -w ${whitelist}
  >>>
  
  output {
    File stats = "${prefix}.metrics"
  }
  
  runtime {
  	preemptible: 3
  	docker: "talkowski/sv-pipeline-rdtest"
  }
}

# Combine male and female test results
# On chrX, use female test results unless the variant appears only in males
# On chrY, use male test results
task merge_allosomes {
  File male_rdtest
  File female_rdtest
  String chrom

  command <<<
    python3 <<CODE
    import pandas as pd
    males = pd.read_table("${male_rdtest}")
    females = pd.read_table("${female_rdtest}")
    if "${chrom}" == 'Y' or "${chrom}" == 'chrY':
      males.to_csv("${basename(male_rdtest)}.merged.csv", sep='\t', index=False, na_rep='NA')
    else:
      male_only = females.P == 'No_samples_for_analysis'
      females.loc[male_only] = males
      females.to_csv("${basename(male_rdtest)}.merged.csv", sep='\t', index=False, na_rep='NA')
    CODE
  >>>

  output {
    File merged_rdtest = "${basename(male_rdtest)}.merged.csv"
  }
  
  runtime {
  	preemptible: 3
  	docker: "talkowski/sv-pipeline-remote-pysam"
  }
}

# Merge split rdtest results into single file
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