# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# This is an analysis WDL to perform FILTER cleanup and recalibrate
# variant QUAL scores at the end of the Talkowski SV pipeline


workflow filter_cleanup_qual_recalibration {
  File vcf
  File vcf_idx
  File contiglist
  File PCRPLUS_samples_list
  String prefix

  Array[Array[String]] contigs = read_tsv(contiglist)

  scatter ( contig in contigs ) {
    call cleanup {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        contig=contig[0],
        PCRPLUS_samples_list=PCRPLUS_samples_list,
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
  String prefix

  command <<<
    tabix -h ${vcf} ${contig} \
    | bgzip -c \
    > input.vcf.gz
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/filter_cleanup_and_QUAL_recalibration.py \
      input.vcf.gz \
      ${PCRPLUS_samples_list} \
      stdout \
    | bgzip -c \
    > "${prefix}.cleaned_filters_qual_recalibrated.vcf.gz"
    # tabix -p vcf -f "${prefix}.cleaned_filters_qual_recalibrated.vcf.gz"
  >>>

  output {
    File out_vcf = "${prefix}.cleaned_filters_qual_recalibrated.vcf.gz"
    # File out_vcf_idx = "${prefix}.cleaned_filters_qual_recalibrated.vcf.gz.tbi"
  }

  runtime {
    docker : "talkowski/sv-pipeline@sha256:112c926b786933690781a93a5c0153fc56568fe5f2e03722d067d5fd6c2bf046"
    preemptible: 0
    disks: "local-disk 50 HDD"
    memory: "4 GB"
  }
}


#General task to combine multiple VCFs
#NO SORT
task concat_vcfs {
  Array[File] vcfs
  String outfile_prefix

  command <<<
    vcf-concat ${sep=' ' vcfs} | bgzip -c > ${outfile_prefix}.vcf.gz; 
    tabix -p vcf -f "${outfile_prefix}.vcf.gz"
  >>>

  output {
    File concat_vcf = "${outfile_prefix}.vcf.gz"
    File concat_vcf_idx = "${outfile_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:5a5bb420b823b129a15ca49ac2e6d3862984c5dd1b8ec91a49ee1ba803685d83"
    preemptible: 0
    memory: "8 GB"
    disks: "local-disk 250 HDD"
  }
}

