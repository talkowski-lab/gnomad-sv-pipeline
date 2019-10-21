import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:05_CleanVCF2/versions/35/plain-WDL/descriptor" as Clean2
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:05_CleanVCF3/versions/23/plain-WDL/descriptor" as clean4

workflow CleanVCF {

  File vcf
  String Chr
  File backgroundlist
  File bothsides_pass_list
  File famfile
  String prefix
  Int max_shards_per_chrom_step1
  Int min_records_per_shard_step1
  Int samples_per_step2_shard
  File? outlier_samples_list

  call subsetvcf {
    input:
      vcf=vcf,
      chr=Chr
  }

  call shard_vcf {
    input:
      vcf=subsetvcf.splitvcf,
      prefix="${prefix}.${Chr}",
      n_shards=max_shards_per_chrom_step1,
      min_vars_per_shard=min_records_per_shard_step1
  }

  scatter ( vcf_shard in shard_vcf.vcf_shards ) {
    call cleanvcf1a {
      input:
        VCF=vcf_shard,
        Backgroundlist=backgroundlist,
        bothsides_pass_list=bothsides_pass_list,
        Famfile=famfile
    }
  }

  call combine_step1_outputs {
    input:
      intermediate_vcfs=cleanvcf1a.intermediate_vcf,
      sexchr_revises=cleanvcf1a.Sex,
      prefix=prefix
  }

  call cleanvcf1b {
    input:
      intermediate_vcf=combine_step1_outputs.merged_intermediate_vcf
  }

  call split {
    input:
      whitelist=cleanvcf1a.Whitelist[0],
      samples_per_step2_shard=samples_per_step2_shard
  }

  call Clean2.Clean {
    input:
    	whitelists=split.whitelists,
      normal_revise_vcf=cleanvcf1b.Normal,
      multi_cnvs=cleanvcf1b.Multi,
      vcftools_idx=cleanvcf1b.vcftools_idx
  }

  call cleanvcf3 {
    input:
  	 RD_CN_revise=Clean.out
  }

  call clean4.Clean4 {
    input:
    	RD_CN_revises=cleanvcf3.shards,
      normal_revise_vcf=cleanvcf1b.Normal,
  }

  call cleanvcf5 {
    input:
    	revise_vcf_lines=Clean4.out,
      normal_revise_vcf=cleanvcf1b.Normal,
      famfile=famfile,
      sexchr_revise=combine_step1_outputs.merged_sexchr_revise,
      multi_IDs=Clean4.multi_IDs,
      outlier_samples_list=outlier_samples_list
  }

  call drop_redundant_CNVs {
    input:
      vcf=cleanvcf5.Polished,
		  Chr=Chr
  }

  call stitch_fragmented_CNVs {
    input:
      vcf=drop_redundant_CNVs.cleaned_vcf_shard,
      chr=Chr,
      prefix=prefix
  }

  call final_cleanup {
    input:
      vcf=stitch_fragmented_CNVs.stitched_vcf_shard,
      chr=Chr,
      prefix=prefix
  }
  
  output {
  	File out=final_cleanup.final_cleaned_shard
  }
}


task subsetvcf {
	
  File vcf
  String chr

  command {
    	tabix -p vcf ${vcf};
      tabix -h ${vcf} ${chr} | bgzip -c > test.${chr}.vcf.gz
  }

  runtime {
  	preemptible: 1
    docker: "talkowski/sv-pipeline@sha256:703a19f84f498989ba8ffde110a3462cfecfbd7ade1084a151fac5fff742c266"
    disks: "local-disk 500 SSD"
    memory: "8 GB"
  }

  output {
  	File splitvcf="test.${chr}.vcf.gz"
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
    docker: "talkowski/sv-pipeline@sha256:703a19f84f498989ba8ffde110a3462cfecfbd7ade1084a151fac5fff742c266"
    preemptible: 1
    maxRetries: 1
    memory: "4 GB"
    disks: "local-disk 500 SSD"
  }
}


#CleanVCF 1a is sharded
task cleanvcf1a {

	File VCF
  File Backgroundlist
  File bothsides_pass_list
  File Famfile

  command {
  	bash /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part1.sh \
           ${VCF} \
           ${Backgroundlist} \
           ${Famfile}
    /opt/sv-pipeline/04_variant_resolution/scripts/add_bothsides_support_filter.py \
      --bgzip \
      --outfile int.w_bothsides.vcf.gz \
      int.vcf.gz \
      ${bothsides_pass_list}
  }
  
  runtime {
    preemptible: 1
    docker: "talkowski/sv-pipeline@sha256:358b979a14d316639fa8b0a9b35edae40ce08470ac4eda8ed4abce9673b42e85"
    disks: "local-disk 30 HDD"
    memory: "4 GB"
    maxRetries: 1
  }

  output {
    # File Combined="all.combined.bed.gz"         #Moved to 2
    # File Multi="multi.cnvs.txt"                 #Moved to 1b
    # File Normal="normal.revise.vcf.gz"          #Moved to 1b
    File Whitelist="whitelist.txt"
    File Sex="sexchr.revise.txt"
    File intermediate_vcf="int.w_bothsides.vcf.gz"
    # File FullVar="fullvar.afternormal.list.txt" #No longer needed
  }
}


task combine_step1_outputs {
  Array[File] intermediate_vcfs
  Array[File] sexchr_revises
  String prefix

  command <<<
    vcf-concat ${sep=' ' intermediate_vcfs} | vcf-sort | bgzip -c \
      > ${prefix}.cleanVCF_step1.intermediate_vcf.merged.vcf.gz;
    cat ${sep=' ' sexchr_revises} > ${prefix}.cleanVCF_step1.sexchr_revise.merged.txt
  >>>

  output {
    File merged_intermediate_vcf = "${prefix}.cleanVCF_step1.intermediate_vcf.merged.vcf.gz"
    File merged_sexchr_revise = "${prefix}.cleanVCF_step1.sexchr_revise.merged.txt"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:703a19f84f498989ba8ffde110a3462cfecfbd7ade1084a151fac5fff742c266"
    preemptible: 1
    maxRetries: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
  }
}


task cleanvcf1b {

  File intermediate_vcf

  command {
    bash /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part1b.sh ${intermediate_vcf}
  }
  
  runtime {
    preemptible: 1
    maxRetries: 1
    docker: "talkowski/sv-pipeline@sha256:703a19f84f498989ba8ffde110a3462cfecfbd7ade1084a151fac5fff742c266"
    disks: "local-disk 100 SSD"
    memory: "4 GB"
  }

  output {
    File Multi="multi.cnvs.txt"
    File Normal="normal.revise.vcf.gz"
    File vcftools_idx = "normal.revise.vcf.gz.csi"
    # File int_afternormalfix = "int.afternormalfix.bed.gz" #Moved to step 2
  }
}


task split {

  File whitelist
  Int samples_per_step2_shard

  command {
      mkdir output;
      split -l ${samples_per_step2_shard} ${whitelist} output/whiteblack.
  }

  runtime {
    preemptible: 1
    maxRetries: 1
    docker: "talkowski/sv-pipeline@sha256:703a19f84f498989ba8ffde110a3462cfecfbd7ade1084a151fac5fff742c266"
    disks: "local-disk 100 SSD"
    memory: "4 GB"
  }

  output {
   Array[File] whitelists=glob("output/*")
  }
}


task cleanvcf3{

	File RD_CN_revise

  command {
    bash /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part3.sh ${RD_CN_revise} 
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:703a19f84f498989ba8ffde110a3462cfecfbd7ade1084a151fac5fff742c266"
    disks: "local-disk 100 SSD"
    memory: "8 GB"
  }
  
  output {
	 Array[File] shards = glob("shards/*")
  }
}


task cleanvcf5 {

  File revise_vcf_lines
  File normal_revise_vcf
  File famfile
  File sexchr_revise
  File multi_IDs
  File? outlier_samples_list
  
  command {
    if [ ${default="SKIP" outlier_samples_list} == "SKIP" ]; then
      echo "" > outliers.txt
    else
      cat ${outlier_samples_list} > outliers.txt
    fi
    touch polished.vcf.gz
    bash /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part5.sh ${revise_vcf_lines} ${normal_revise_vcf} ${famfile} ${sexchr_revise} ${multi_IDs} outliers.txt
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:703a19f84f498989ba8ffde110a3462cfecfbd7ade1084a151fac5fff742c266"
    disks: "local-disk 100 SSD"
    memory: "8 GB"
    maxRetries: 1
  }
  
  output {
    File Polished="polished.vcf.gz"
  }
}


task drop_redundant_CNVs {

  File vcf
  String Chr

  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/resolve_CPX_CNV_redundancies.sh \
      ${vcf} \
      ${Chr}.shard.no_CNV_redundancies.vcf.gz
  >>>

  output {
    File cleaned_vcf_shard = "${Chr}.shard.no_CNV_redundancies.vcf.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:703a19f84f498989ba8ffde110a3462cfecfbd7ade1084a151fac5fff742c266"
    preemptible: 0
    maxRetries: 1
    disks: "local-disk 100 SSD"
  }
}


# Stitch fragmented RD-only calls found in 100% of the same samples
task stitch_fragmented_CNVs {

  File vcf
  String chr
  String prefix

  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/stitch_fragmented_CNVs.sh \
      ${vcf} \
      "${chr}.shard.fragmented_CNVs_stitched1.vcf.gz"
    /opt/sv-pipeline/04_variant_resolution/scripts/stitch_fragmented_CNVs.sh \
      "${chr}.shard.fragmented_CNVs_stitched1.vcf.gz" \
      "${chr}.shard.fragmented_CNVs_stitched.vcf.gz"
  >>>

  output {
    File stitched_vcf_shard = "${chr}.shard.fragmented_CNVs_stitched.vcf.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:703a19f84f498989ba8ffde110a3462cfecfbd7ade1084a151fac5fff742c266"
    preemptible: 0
    maxRetries: 1
    disks: "local-disk 100 SSD"
  }
}


# Final VCF cleanup
task final_cleanup {
  File vcf
  String chr
  String prefix

  command <<<
    /opt/sv-pipeline/04_variant_resolution/scripts/rename_after_vcfcluster.py \
      --chrom ${chr} \
      --prefix ${prefix} \
      ${vcf} stdout \
    | fgrep -v "##INFO=<ID=HIGH_SR_BACKGROUND" \
    | /opt/sv-pipeline/04_variant_resolution/scripts/sanitize_filter_field.py stdin stdout \
    | fgrep -v "##INFO=<ID=MEMBERS,Number=.,Type=String," \
    | bgzip -c > "${prefix}.${chr}.final_cleanup.vcf.gz"
  >>>

  output {
    File final_cleaned_shard = "${prefix}.${chr}.final_cleanup.vcf.gz"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:703a19f84f498989ba8ffde110a3462cfecfbd7ade1084a151fac5fff742c266"
    preemptible: 0
    maxRetries: 1
    disks: "local-disk 100 SSD"
  }
}
