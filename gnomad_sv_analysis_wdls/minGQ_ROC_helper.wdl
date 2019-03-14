# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


#This is a helper WDL that performs minGQ ROC optimization for a filtered VCF
#See parent WDL: optimize_GQ_filter

workflow optimize_ROC {
  File vcf
  Array[File] famfiles
  String prefix
  String max_fdr
  Int minGQ
  Int maxGQ
  Int GQstepsize

  #Get de novo stats for each famfile shard
  scatter( famfile in famfiles ){
    call gather_denovo_stats as gather_stats {
      input:
        vcf=vcf,
        famfile=famfile,
        minGQ=minGQ,
        maxGQ=maxGQ,
        GQstepsize=GQstepsize
    }

    call cat_denovo_stats as cat_stats_pershard {
      input:
        dn_stats=gather_stats.dn_stats_glob,
        prefix="${prefix}.pershard"
    }
  }

  #Merge de novo stats files across all shards
  call cat_denovo_stats as cat_stats_allshards {
    input:
      dn_stats=cat_stats_pershard.merged_stats,
      prefix="${prefix}"
  }

  #Run ROC analysis
  call ROC_optimization {
    input:
      merged_stats=cat_stats_allshards.merged_stats,
      prefix="${prefix}",
      max_fdr=max_fdr
  }

  output {
    File minGQ_ROC_plot = ROC_optimization.tarball
    File minGQ_ROC_table = ROC_optimization.minGQ_table
  }
}


#Gather de novo stats for a set of trios from a vcf
task gather_denovo_stats {
  File vcf
  File famfile
  Int minGQ
  Int maxGQ
  Int GQstepsize

  command <<<
    #Get list of sample IDs & column numbers from VCF header
    zcat ${vcf} | head -n1000 | fgrep "#" | fgrep -v "##" | sed 's/\t/\n/g' \
      | awk -v OFS="\t" '{ print $1, NR }' > vcf_header_columns.txt
    #Iterate over families & subset VCF
    while read famID pro fa mo prosex pheno; do
      pro_idx=$( awk -v ID=$pro '{ if ($1==ID) print $2 }' vcf_header_columns.txt )
      fa_idx=$( awk -v ID=$fa '{ if ($1==ID) print $2 }' vcf_header_columns.txt )
      mo_idx=$( awk -v ID=$mo '{ if ($1==ID) print $2 }' vcf_header_columns.txt )
      if ! [ -z $pro_idx ] && ! [ -z $fa_idx ] && ! [ -z $mo_idx ]; then
        echo -e "ANALYZING $famID, which contains $pro, $fa, and $mo"
        #Get variant stats
        zcat ${vcf} | cut -f1-9,"$pro_idx","$fa_idx","$mo_idx" \
          | fgrep -v "MULTIALLELIC" \
          | /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/gather_trio_genos.py \
            stdin stdout "$pro" "$fa" "$mo" \
          | gzip -c \
          > "$famID".trio_variant_info.txt.gz
        #Titrate GQs & count de novos
        /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/compute_denovo_stats.R \
          --famID "$famID" \
          --minGQ ${minGQ} \
          --maxGQ ${maxGQ} \
          --step ${GQstepsize} \
          "$famID".trio_variant_info.txt.gz \
          "$famID".trio_denovo_summary.txt
        gzip -f "$famID".trio_denovo_summary.txt
      fi
    done < ${famfile}
  >>>

  output {
    Array[File] dn_stats_glob = glob("*.trio_denovo_summary.txt.gz")
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:d3844f6c7c26da55e679c9c521882d54dbecd169f884f09f05a12d6565bf6063"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 250 HDD"
  }
}


#Combine de novo stats into one long melted table
task cat_denovo_stats {
  Array[File] dn_stats
  String prefix

  command <<<
    zcat ${dn_stats[0]} | sed -n '1p' > ${prefix}.merged_denovo_stats.txt
    while read statsfile; do
      zcat "$statsfile" | sed '1d'
    done < ${write_lines(dn_stats)} \
      | sort -nk1,1 \
      >> ${prefix}.merged_denovo_stats.txt
    gzip -f ${prefix}.merged_denovo_stats.txt
  >>>

  output {
    File merged_stats = "${prefix}.merged_denovo_stats.txt.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:d3844f6c7c26da55e679c9c521882d54dbecd169f884f09f05a12d6565bf6063"
    preemptible: 1
  }
}


#Run de novo ROC analysis
task ROC_optimization {
  File merged_stats
  String prefix
  String max_fdr

  command <<<
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/optimize_GQ_ROC.R \
      --prefix ${prefix} \
      --fdr ${max_fdr} \
      -S /opt/sv-pipeline/ref/vcf_qc_refs/SV_colors.txt \
      ${merged_stats} \
      ./${prefix}_minGQ_ROC_results/
      cp ./${prefix}_minGQ_ROC_results/${prefix}.minGQ_ROC_results.txt \
      ${prefix}.minGQ_ROC_results.txt
    tar -czvf ${prefix}_minGQ_ROC_results.tar.gz ./${prefix}_minGQ_ROC_results
  >>>

  output {
    File tarball = "${prefix}_minGQ_ROC_results.tar.gz"
    File minGQ_table = "${prefix}.minGQ_ROC_results.txt"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:d3844f6c7c26da55e679c9c521882d54dbecd169f884f09f05a12d6565bf6063"
    preemptible: 1
    memory: "8 GB"
    disks: "local-disk 20 HDD"
  }
}