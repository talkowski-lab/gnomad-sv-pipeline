# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# Workflow to parallelize integration of all-variant and inv-only svtk resolve results per chromosome
workflow integrate_invonly_allvars {
  String inv_res_vcf
  String all_res_vcf
  File inv_res_vcf_idx
  File all_res_vcf_idx
  String prefix
  File contiglist
  File svc_acct_key
  File bothside_pass
  File background_fail

  Array[Array[String]] contigs = read_tsv(contiglist)

  #Merge, scattered by chromosome
  scatter (contig in contigs) {

    #Remote tabix each vcf
    call subset_vcf as subset_inv {
      input:
        vcf=inv_res_vcf,
        vcf_idx=inv_res_vcf_idx,
        contig=contig[0],
        prefix="${prefix}.inv_only.${contig[0]}",
        svc_acct_key=svc_acct_key
    }
    call subset_vcf as subset_all {
      input:
        vcf=all_res_vcf,
        vcf_idx=all_res_vcf_idx,
        contig=contig[0],
        prefix="${prefix}.all_variants.${contig[0]}",
        svc_acct_key=svc_acct_key
    }

    #Run integration per chromosome
    call integrate_resolved_vcfs {
      input:
        inv_res_vcf=subset_inv.subsetted_vcf,
        all_res_vcf=subset_all.subsetted_vcf,
        prefix="${prefix}.resolved.${contig[0]}"
    }
  }

  #Merge integrated vcfs across chromosomes
  call concat_vcfs {
    input:
      vcfs=integrate_resolved_vcfs.integrated_vcf,
      prefix="${prefix}.resolved"
  }

  output {
    File integrated_vcf = concat_vcfs.concat_vcf
    File integrated_vcf_idx = concat_vcfs.concat_vcf_idx
  }
}


#Remote tabix a single chromosome per VCFs
task subset_vcf {
  String vcf
  File vcf_idx
  String contig
  String prefix
  File svc_acct_key

  command <<<
    #Remote tabix to chromosome of interest
    url=$( gsutil signurl -d 24h ${svc_acct_key} "$vcf" | sed '1d' | cut -f 4 );
    echo $url;
    svtk remote_tabix --header "$url" ${vcf_idx} "${contig}:0-300000000" > "${prefix}.${contig}.vcf"
    bgzip -f "${prefix}.${contig}.vcf"
    tabix -p vcf -f "${prefix}.${contig}.vcf.gz"
  >>>

  output {
    File subsetted_vcf = "${prefix}.${contig}.vcf.gz"
    File subsetted_vcf_idx = "${prefix}.${contig}.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline-remote-pysam@sha256:9fd37fb64e28e54d53172dd30d68c36f0815f21af465381dac281d53755edd86"
    preemptible: 1
    disks: "local-disk 50 SSD"
  }
}


# Merge inversion-only and all-variant cpx-resolved outputs
task integrate_resolved_vcfs {
  File inv_res_vcf
  File all_res_vcf
  String prefix

  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/Complex_Inversion_Integration.sh \
      ${inv_res_vcf} \
      ${all_res_vcf} \
      ${prefix}.integrated_resolved.vcf.gz
  >>>

  output {
    File integrated_vcf = "${prefix}.integrated_resolved.vcf.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:b359f2cb0c9d5f5a55eb4c41fd362f4e574bf3f8f0f395a2907837571b367ee0"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 250 SSD"
  }
}


#Merge multiple vcfs
task concat_vcfs {
  Array[File] vcfs
  String prefix

  command <<<
    vcf-concat ${sep=' ' vcfs} | vcf-sort -c | bgzip -c > ${prefix}.vcf.gz
    tabix -f -p vcf ${prefix}.vcf.gz
  >>>

  output {
    File concat_vcf = "${prefix}.vcf.gz"
    File concat_vcf_idx = "${prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:b359f2cb0c9d5f5a55eb4c41fd362f4e574bf3f8f0f395a2907837571b367ee0"
    preemptible: 1
    disks: "local-disk 1000 SSD"
  }
}