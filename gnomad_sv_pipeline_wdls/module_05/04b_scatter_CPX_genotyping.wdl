import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04b_genotype_CPX_CNVs/versions/12/plain-WDL/descriptor" as cpx_gt

# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# Workflow to perform depth-based genotyping for a single vcf shard scattered 
# across batches on predicted CPX CNVs from 04b
workflow scatter_CPX_genotyping {
  File vcf
  Int n_master_vcf_shards
  Int n_master_min_vars_per_vcf_shard
  File gt_input_files
  Int n_per_split_small
  Int n_per_split_large
  Int n_RdTest_bins
  File svc_acct_key
  String prefix
  File famfile
  String contig

  # Shard VCF into even slices
  call shard_vcf {
    input:
      vcf=vcf,
      prefix="${prefix}.${contig}",
      n_shards=n_master_vcf_shards,
      min_vars_per_shard=n_master_min_vars_per_vcf_shard
  }

  # Scatter genotyping over shards
  scatter ( shard in shard_vcf.vcf_shards ) {
    # Run genotyping
    call cpx_gt.genotype_CPX_CNVs as genotype_shard {
      input:
        vcf=shard,
        gt_input_files=gt_input_files,
        n_per_split_large=n_per_split_large,
        n_per_split_small=n_per_split_small,
        n_RdTest_bins=n_RdTest_bins,
        svc_acct_key=svc_acct_key,
        prefix=prefix,
        famfile=famfile,
        contig=contig
    }
  }

  # Merge VCF shards
  call concat_vcfs {
    input:
      vcfs=genotype_shard.cpx_depth_gt_resolved_vcf,
      outfile_prefix="${prefix}.${contig}.resolved"
  }

  # Output merged VCF
  output {
    File cpx_depth_gt_resolved_vcf = concat_vcfs.concat_vcf
    File cpx_depth_gt_resolved_vcf_idx = concat_vcfs.concat_vcf_idx
  }
 }




#Shard a vcf into even chunks
task shard_vcf {
  File vcf
  String prefix
  Int n_shards
  Int min_vars_per_shard

  command <<<
    zcat ${vcf} | sed -n '1,1000p' | fgrep "#" > header.vcf;
    nrecords=$( zcat ${vcf} | cut -f1 | fgrep -v "#" | wc -l );
    rec_per_shard=$( echo "$(( $nrecords / ${n_shards} ))" | cut -f1 -d\. );
    if [ $rec_per_shard -lt ${min_vars_per_shard} ]; then
      rec_per_shard=${min_vars_per_shard}
    fi;
    zcat ${vcf} | fgrep -v "#" \
    | split -l $rec_per_shard --numeric-suffixes=001 -a 3 /dev/stdin vcf_records_ ;
    max_suf=$( find `pwd` -name "vcf_records_*" | awk -v FS="_" '{ print $NF }' | sort -nrk1,1 | sed -n '1p' )
    for i in $( seq -w 001 "$max_suf" ); do
      cat header.vcf vcf_records_"$i" | bgzip -c > ${prefix}.shard_"$i".vcf.gz
      rm vcf_records_"$i"
    done
  >>>

  output {
    Array[File] vcf_shards = glob("${prefix}.shard_*.vcf.gz")
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:c2af5febc8967dff0b7a10cd764b292f43029ffd119e40832cef3fcbc3df1c1f"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 500 HDD"
  }
}


#General task to combine multiple VCFs
task concat_vcfs {
  Array[File] vcfs
  String outfile_prefix

  command <<<
    vcf-concat ${sep=' ' vcfs} | vcf-sort -c | bgzip -c > ${outfile_prefix}.vcf.gz; 
    tabix -p vcf -f "${outfile_prefix}.vcf.gz"
  >>>

  output {
    File concat_vcf = "${outfile_prefix}.vcf.gz"
    File concat_vcf_idx = "${outfile_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:c2af5febc8967dff0b7a10cd764b292f43029ffd119e40832cef3fcbc3df1c1f"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 500 HDD"
  }
}


