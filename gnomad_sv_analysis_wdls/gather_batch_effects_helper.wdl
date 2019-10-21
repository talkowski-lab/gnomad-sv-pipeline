# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# This is an analysis WDL to perform pairwise comparisons of batches in the 
# Talkowski lab SV pipeline, and mark sites that appear batch-specific


workflow check_batch_effects {
	File freq_table
  String batch1
  String batch2
  String prefix
  Int variants_per_shard

  # Shard frequency table
  call shard_table {
    input:
      freq_table=freq_table,
      variants_per_shard=variants_per_shard
  }

  # Scatter over shards and compute AF correlations for each variant
  scatter ( shard in shard_table.shards ) {
    call compare_batches {
      input:
        freq_table=shard,
        batch1=batch1,
        batch2=batch2,
        prefix=prefix
    }
  }

  # Combine shards, perform bonferroni correction to determine significant batch effects, and plot AF correlation scatter
  call combine_shards {
    input:
      freq_tables=compare_batches.results,
      batch1=batch1,
      batch2=batch2,
      prefix=prefix
  }

  # Outputs
  output {
    File comparison_table = combine_shards.merged_table
    File batch_effect_variants = combine_shards.batch_effect_variants
    File scatterplots_tarball = combine_shards.correlation_scatterplots_tarball
  }
}


# Shard a frequency table into an even number of evenly sized shards
task shard_table {
  File freq_table
  Int variants_per_shard

  command <<<
    set -euo pipefail
    #Split variant lines
    zcat ${freq_table} | sed '1d' | \
    split -l ${variants_per_shard} --numeric-suffixes=00001 -a 5 /dev/stdin freq_table_shard_ || true
    #Add header & gzip each shard
    zcat ${freq_table} | sed -n '1p' > header.txt
    maxshard=$( find / -name "freq_table_shard_*" | awk -v FS="_" '{ print $NF }' \
                | sort -Vrk1,1 | sed -n '1p' || true )
    for i in $( seq -w 00001 "$maxshard" ); do
      cat header.txt "freq_table_shard_$i" \
      | gzip -c \
      > "freq_table_shard_$i.txt.gz" || true
    done
  >>>

  output {
    Array[File] shards = glob("freq_table_shard_*.txt.gz")
  }

  runtime {
    docker : "talkowski/sv-pipeline@sha256:aef8156983cec6ac6a91fa6461b197a63835e5487fc9523ec857f947cfac660e"
    preemptible: 1
    maxRetries: 1
  }
}


# Compare AF stats per variant between a pair of batches
task compare_batches {
  File freq_table
  String batch1
  String batch2
  String prefix

  command <<<
    set -euo pipefail
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/find_batch_effects.shard_helper.R \
      ${freq_table} \
      "${batch1}" \
      "${batch2}" \
      "${prefix}.${batch1}_vs_${batch2}.results.txt"
    gzip "${prefix}.${batch1}_vs_${batch2}.results.txt"
  >>>

  output {
    File results = "${prefix}.${batch1}_vs_${batch2}.results.txt.gz"
  }

  runtime {
    docker : "talkowski/sv-pipeline@sha256:aef8156983cec6ac6a91fa6461b197a63835e5487fc9523ec857f947cfac660e"
    memory: "4 GB"
    preemptible: 1
    maxRetries: 1
  }
}


# Merge sharded comparison results and perform analysis for batch effects
task combine_shards {
  Array[File] freq_tables
  String batch1
  String batch2
  String prefix
  
  command <<<
    set -euo pipefail
    #Write header
    zcat ${freq_tables[0]} | sed -n '1p' > header.txt || true
    #Iterate over files and cat
    while read file; do
      zcat "$file" | sed '1d'
    done < ${write_lines(freq_tables)} \
    | cat header.txt - \
    | gzip -c \
    > "${prefix}.${batch1}_vs_${batch2}.AF_comparison_table.txt.gz" || true
    #Analyze
    mkdir "${batch1}_vs_${batch2}"
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/find_batch_effects.R \
      "${prefix}.${batch1}_vs_${batch2}.AF_comparison_table.txt.gz" \
      "${batch1}" \
      "${batch2}" \
      "${batch1}_vs_${batch2}/${prefix}"
    gzip -f "${batch1}_vs_${batch2}/${prefix}.${batch1}_vs_${batch2}.freq_table_wBonferroni.txt"
    tar -czvf "${batch1}_vs_${batch2}.tar.gz" \
      "${batch1}_vs_${batch2}"
  >>>

  output {
    File merged_table = "${batch1}_vs_${batch2}/${prefix}.${batch1}_vs_${batch2}.freq_table_wBonferroni.txt.gz"
    File batch_effect_variants = "${batch1}_vs_${batch2}/${prefix}.${batch1}_vs_${batch2}.batch_effect_variants.txt"
    File correlation_scatterplots_tarball = "${batch1}_vs_${batch2}.tar.gz"
  }

  runtime {
    docker : "talkowski/sv-pipeline@sha256:aef8156983cec6ac6a91fa6461b197a63835e5487fc9523ec857f947cfac660e"
    memory: "4 GB"
    preemptible: 1
    maxRetries: 1
  }
}

