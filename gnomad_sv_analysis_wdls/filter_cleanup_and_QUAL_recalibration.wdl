# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# This is an analysis WDL to perform FILTER cleanup and recalibrate
# variant QUAL scores at the end of the Talkowski SV pipeline


workflow filter_cleanup_qual_recalibration {
  File vcf
  File vcf_idx
  File PCRPLUS_samples_list
  File famfile
  Float min_callrate_global
  Float min_callrate_smallDels
  File contiglist
  String prefix

  Array[Array[String]] contigs = read_tsv(contiglist)

  scatter ( contig in contigs ) {
    call cleanup {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        contig=contig[0],
        PCRPLUS_samples_list=PCRPLUS_samples_list,
        famfile=famfile,
        min_callrate_global=min_callrate_global,
        min_callrate_smallDels=min_callrate_smallDels,
        prefix=prefix
    }
  }

  call concat_vcfs {
    input:
      vcfs=cleanup.out_vcf,
      outfile_prefix="${prefix}.cleaned_filters_qual_recalibrated"
  }

    output {
      File cleaned_vcf = concat_vcfs.concat_vcf
      File cleaned_vcf_idx = concat_vcfs.concat_vcf_idx
    }
}


# Applies filters & cleanup to VCF for a single chromosome
task cleanup {
  File vcf
  File vcf_idx
  String contig
  File PCRPLUS_samples_list
  File famfile
  Float min_callrate_global
  Float min_callrate_smallDels
  String prefix

  command <<<
    set -euo pipefail
    #Subset to chromosome of interest
    tabix -h ${vcf} ${contig} | bgzip -c > input.vcf.gz
    #Get list of PCR- samples
    tabix -H ${vcf} | fgrep -v "##" | cut -f10- | sed 's/\t/\n/g' \
    > all.samples.list
    fgrep -wvf ${PCRPLUS_samples_list} all.samples.list \
    > pcrminus.samples.list
    #Restrict famfiles
    while read ptn; do fgrep -w $ptn ${famfile}; done < all.samples.list > revised.fam
    fgrep -wf pcrminus.samples.list revised.fam > revised.pcrminus.fam
    #Compute fraction of missing genotypes per variant
    zcat input.vcf.gz \
    | awk '{ if ($7 !~ /MULTIALLELIC/) print $0 }' \
    | bgzip -c \
    > input.noMCNV.vcf.gz
    plink2 \
      --missing variant-only \
      --max-alleles 2 \
      --keep-fam revised.pcrminus.fam \
      --fam revised.fam \
      --vcf input.noMCNV.vcf.gz
    fgrep -v "#" plink2.vmiss \
    | awk -v OFS="\t" '{ print $2, 1-$NF }' \
    > callrates.txt
    #Clean up VCF
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/filter_cleanup_and_QUAL_recalibration.py \
      --callrate-table callrates.txt \
      --min-callrate-global ${min_callrate_global} \
      --min-callrate-smallDels ${min_callrate_smallDels} \
      input.vcf.gz \
      stdout \
    | bgzip -c \
    > "${prefix}.${contig}.cleaned_filters_qual_recalibrated.vcf.gz"
    # tabix -p vcf -f "${prefix}.cleaned_filters_qual_recalibrated.vcf.gz"
  >>>

  output {
    File out_vcf = "${prefix}.${contig}.cleaned_filters_qual_recalibrated.vcf.gz"
    # File out_vcf_idx = "${prefix}.cleaned_filters_qual_recalibrated.vcf.gz.tbi"
  }

  runtime {
    docker : "talkowski/sv-pipeline@sha256:4587376100d71d66fb864740f95e0cc5f343bb1fe6e892f5b8116c789c38333f"
    preemptible: 1
    maxRetries: 0
    disks: "local-disk 50 HDD"
    memory: "4 GB"
  }
}


#General task to combine and sort multiple VCFs
task concat_vcfs {
  Array[File] vcfs
  String outfile_prefix

  command <<<
    set -euo pipefail
    vcf-concat ${sep=' ' vcfs} | bgzip -c > ${outfile_prefix}.vcf.gz; 
    tabix -p vcf -f "${outfile_prefix}.vcf.gz"
  >>>

  output {
    File concat_vcf = "${outfile_prefix}.vcf.gz"
    File concat_vcf_idx = "${outfile_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:4587376100d71d66fb864740f95e0cc5f343bb1fe6e892f5b8116c789c38333f"
    preemptible: 0
    maxRetries: 1
    memory: "4 GB"
    bootDiskSizeGb: 30
    disks: "local-disk 250 HDD"
  }
}

