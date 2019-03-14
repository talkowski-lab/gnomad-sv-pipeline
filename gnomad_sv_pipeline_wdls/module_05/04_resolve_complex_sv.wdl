import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04_resolve_complex_by_chrom/versions/63/plain-WDL/descriptor" as resolve_complex_by_chrom

workflow resolve_complex_sv {
  File vcf
  File contigs
  Int max_shards_per_chrom
  Int min_variants_per_shard
  File cytobands
  File cytobands_idx
  File mei_bed
  File discfile_list
  File discfile_idx_list
  File pe_blacklist
  File pe_blacklist_idx
  File svc_acct_key
  File rf_cutoffs

  Array[Array[String]] contiglist = read_tsv(contigs)

  # Get SR count cutoff from RF metrics to use in single-ender rescan procedure
  call get_se_cutoff {
    input:
      rf_cutoffs=rf_cutoffs
  }


  scatter (contig in contiglist) {
    call subset_vcf {
      input:
        vcf=vcf,
        chrom=contig[0]
    }

    call resolve_complex_by_chrom.resolve_complex_by_chrom as resolve_perChrom {
      input:
        vcf=subset_vcf.single_chrom,
        vcf_idx=subset_vcf.single_chrom_idx,
        contig=contig[0],
        max_shards=max_shards_per_chrom,
        min_variants_per_shard=min_variants_per_shard,
        cytobands=cytobands,
        cytobands_idx=cytobands_idx,
        discfile_list=discfile_list,
        discfile_idx_list=discfile_idx_list,
        mei_bed=mei_bed,
        pe_blacklist=pe_blacklist,
        pe_blacklist_idx=pe_blacklist_idx,
        svc_acct_key=svc_acct_key,
        se_pe_cutoff=get_se_cutoff.median_PE_cutoff
    }
  }

  call resolve_complex_by_chrom.concat_vcfs as concat_resolved {
    input:
      vcfs=resolve_perChrom.res_vcf,
      vcftype="resolved"
  }

  # call resolve_complex_by_chrom.concat_vcfs as concat_unresolved {
  #   input:
  #     vcfs=resolve_perChrom.unres_vcf,
  #     vcftype="unresolved"
  # }

  output {
    File resolved_vcf_merged = concat_resolved.concat_vcf
    File resolved_vcf_merged_idx = concat_resolved.concat_vcf_idx
    # File unresolved_vcf_merged = concat_unresolved.concat_vcf
  }
}

#Subset VCF per chromosome
task subset_vcf {
  File vcf
  String chrom

  String prefix = basename(vcf, ".vcf.gz")
  
  command <<<
    tabix -p vcf ${vcf};
    tabix --print-header ${vcf} ${chrom} | bgzip -c > ${prefix}.${chrom}.vcf.gz
    tabix -f ${prefix}.${chrom}.vcf.gz
  >>>

  output {
    File single_chrom = "${prefix}.${chrom}.vcf.gz"
    File single_chrom_idx = "${prefix}.${chrom}.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:96d07aa2c7c3e8bd12f2621a0644a5a8fca99f922926922724497ad2aad9364d"
    preemptible: 3
    disks: "local-disk 1000 SSD"
  }
}

# Get SE cutoff
task get_se_cutoff {
  File rf_cutoffs

  command <<<
    mkdir rf_cutoff_files/
    cat ${rf_cutoffs} | gsutil cp -I rf_cutoff_files/
    while read file; do
      /opt/sv-pipeline/04_variant_resolution/scripts/convert_poisson_p.py \
      $( awk -F '\t' '{if ( $5=="PE_log_pval") print $2 }' $file | head -n1 )
    done < <( find rf_cutoff_files/ -name "*cutoffs" ) | \
    Rscript -e "cat(floor(median(scan('stdin',quiet=T))),sep='\n')" > \
    median_cutoff.txt
  >>>

  output {
    Int median_PE_cutoff = read_tsv("median_cutoff.txt")[0][0]
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:96d07aa2c7c3e8bd12f2621a0644a5a8fca99f922926922724497ad2aad9364d"
    preemptible: 3
  }
}

# 
##Combine multiple VCFs
#task concat_vcfs {
#  Array[File] vcfs
#
#  command <<<
#    vcf-concat ${sep=' ' vcfs} | vcf-sort -c | bgzip -c > all_batches.vcf.gz
#  >>>
#
#  output {
#    File concat_vcf = "all_batches.vcf.gz"
#  }
#
#  runtime {
#    docker: "talkowski/sv-pipeline@sha256:b0455d30df2fbdbd4649466d968cada0a44d02a7159d94982308b629dd1aef78"
#    preemptible: 3
#  }
#}
