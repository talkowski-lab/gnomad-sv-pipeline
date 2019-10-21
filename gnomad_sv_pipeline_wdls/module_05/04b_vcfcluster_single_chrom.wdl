import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:04_vcfcluster_tasks_per_chrom/versions/17/plain-WDL/descriptor" as vcfcluster_tasks

# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License


# Workflow to run parallelized vcf clustering for a single chromosome
workflow vcfcluster_single_chrom {
  File vcf_list
  File vcf_idx_list
  String prefix
  Int dist
  Float frac
  Float sample_overlap
  String do_blacklist
  File blacklist
  File blacklist_idx
  File batches_list
  Int svsize
  Array[String] svtypes
  String contig
  Int max_shards_per_chrom_svtype
  Int min_variants_per_shard_per_chrom_svtype
  File svc_acct_key
  String subset_sr_lists
  File bothside_pass
  File background_fail

    #Remote tabix each vcf & join into a single vcf
    call join_vcfs {
      input:
        vcf_list=vcf_list,
        vcf_idx_list=vcf_idx_list,
        batches_list=batches_list,
        contig=contig,
        prefix=prefix,
        svc_acct_key=svc_acct_key
    }

    #Run vcfcluster per chromosome
    call vcfcluster_tasks.cluster_single_chrom as cluster {
      input:
        vcf=join_vcfs.joined_vcf,
        vcf_idx=join_vcfs.joined_vcf_idx,
        contig=contig,
        prefix=prefix,
        max_shards=max_shards_per_chrom_svtype,
        min_per_shard=min_variants_per_shard_per_chrom_svtype,
        dist=dist,
        frac=frac,
        sample_overlap=sample_overlap,
        do_blacklist=do_blacklist,
        blacklist=blacklist,
        blacklist_idx=blacklist_idx,
        svsize=svsize,
        svtypes=svtypes
    }

    #Subset bothside_pass & background_fail to chromosome of interest
    call subset_variant_list as subset_bothside_pass {
      input:
        vid_list=bothside_pass,
        vcf=join_vcfs.joined_vcf,
        prefix=prefix,
        contig=contig,
        run=subset_sr_lists
    }
    call subset_variant_list as subset_background_fail {
      input:
        vid_list=background_fail,
        vcf=join_vcfs.joined_vcf,
        prefix=prefix,
        contig=contig,
        run=subset_sr_lists
    }

  output {
    File clustered_vcf = cluster.clustered_vcf
    File clustered_vcf_idx = cluster.clustered_vcf_idx
    File filtered_bothside_pass = subset_bothside_pass.filtered_vid_list
    File filtered_background_fail = subset_background_fail.filtered_vid_list
  }
}


#Task to remote tabix a single chromosome for all VCFs, then merge row-wise
task join_vcfs {
  File vcf_list
  File vcf_idx_list
  File batches_list
  String contig
  String prefix
  File svc_acct_key

  command <<<
    set -euo pipefail
    #Remote tabix all vcfs to chromosome of interest
    echo "REMOTE TABIXING VCFs"
    while read batch vcf_path vcf_idx_path; do
      gsutil cp "$vcf_idx_path" "$batch.vcf.gz.idx";
      url=$( gsutil signurl -d 24h ${svc_acct_key} "$vcf_path" | sed '1d' | cut -f 4 );
      echo $url;
      svtk remote_tabix --header "$url" "$batch.vcf.gz.idx" "${contig}:0-300000000" |sed "s/AN=[0-9]*;//g"|sed "s/AC=[0-9]*;//g"> "$batch.${contig}.vcf"
      bgzip -f "$batch.${contig}.vcf"
    done < <( paste ${batches_list} ${vcf_list} ${vcf_idx_list} )
    find `pwd` -name "*.vcf.gz"
    while read batch; do
      find `pwd` -name "$batch.${contig}.vcf.gz"
    done < ${batches_list} > subsetted_vcfs.list
    #Sanity check to make sure all subsetted VCFs have same number of records
    while read vcf; do
      zcat "$vcf" | cut -f1 | fgrep -v "#" | wc -l
    done < subsetted_vcfs.list \
    > records_per_vcf.txt
    if [ $( cat records_per_vcf.txt | sort | uniq | wc -l ) -gt 1 ]; then
      echo -e "ERROR: INCONSISTENT NUMBER OF RECORDS PER VCF DETECTED"
      cat records_per_vcf.txt
      exit 0
    fi
    # echo "CONTENTS OF subsetted_vcfs.list:"
    # while read vcf; do
    #   zcat $vcf | fgrep -v "#" | head -n100 | wc -l
    # done < subsetted_vcfs.list
    #Join vcfs
    /opt/sv-pipeline/04_variant_resolution/scripts/join_vcfs_paste_implementation.sh \
      subsetted_vcfs.list \
      "${prefix}.joined";
    echo "FINISHED join_vcfs_parallel_implementation.sh; RESULTS:"
    find `pwd` -name "${prefix}.joined.vcf*"
    /opt/sv-pipeline/04_variant_resolution/scripts/make_concordant_multiallelic_alts.py \
      $( find `pwd` -name "${prefix}.joined.vcf.gz" ) \
      subsetted_vcfs.list \
      ${prefix}.unclustered.vcf;
    #bgzip ${prefix}.unclustered.vcf
    # zcat ${prefix}.joined.vcf.gz \
    cat ${prefix}.unclustered.vcf \
      | sed -e 's/:RD,PE,SR/:7/g' \
      | sed -e 's/:PE,SR/:6/g' \
      | sed -e 's/:RD,SR/:5/g' \
      | sed -e 's/:RD,PE/:3/g' \
      | sed -e 's/:PE\t/:2\t/g' -e 's/:SR\t/:4\t/g' -e 's/:RD\t/:1\t/g' \
      | sed -e 's/ID=EV,Number=.,Type=String/ID=EV,Number=1,Type=Integer/g' \
      | bgzip -c > ${prefix}.unclustered.vcf.gz
    tabix -f -p vcf "${prefix}.unclustered.vcf.gz"
  >>>

  output {
    File joined_vcf = "${prefix}.unclustered.vcf.gz"
    File joined_vcf_idx = "${prefix}.unclustered.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline-remote-pysam@sha256:e3268dac103ee611214d0fdc44ab0485d4dc1f5f795512cf6bb1f4cfc3da1dc8"
    preemptible: 1
    memory: "8 GB"
    disks: "local-disk 1000 SSD"
  }
}


#Subset a single variant list
task subset_variant_list {
  File vid_list
  File vcf
  String prefix
  String contig
  String run

  command <<<
    set -e
    if [ ${run} == "TRUE" ]; then
      #Get list of variants present in VCF
      zcat ${vcf} | cut -f1-3 | fgrep -v "#" | cut -f3 > valid_vids.list;
      #Restrict input list to valid VIDs
      fgrep -wf valid_vids.list ${vid_list} > "${prefix}.${contig}.VIDs.list"
    else
      echo "" > "${prefix}.${contig}.VIDs.list"
    fi
  >>>

  output {
    File filtered_vid_list = "${prefix}.${contig}.VIDs.list"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:703a19f84f498989ba8ffde110a3462cfecfbd7ade1084a151fac5fff742c266"
    preemptible: 1
    disks: "local-disk 30 SSD"
  }
}
