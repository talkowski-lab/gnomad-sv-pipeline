import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:master_SV_VCF_QC/versions/47/plain-WDL/descriptor" as vcf_qc

workflow cluster_pesr_algorithm {
  Array[File] vcfs
  File contigs
  String batch
  String algorithm
  File famfile
  File trios_famfile
  String ref_build
  File Sanders_2015_tarball
  File collins_2017_tarball
  File Werling_2018_tarball

  # VCFcluster parameters
  Int dist
  Float frac
  File blacklist
  Int svsize
  String svtypes
  String flags

  Array[Array[String]] contiglist = read_tsv(contigs)

  scatter (contig in contiglist) {
    call vcfcluster {
      input:
        vcfs=vcfs,
        batch=batch,
        algorithm=algorithm,
        chrom=contig[0],
        dist=dist,
        frac=frac,
        blacklist=blacklist,
        svsize=svsize,
        flags=flags,
        svtypes=svtypes
    }
  }

  call concat_vcfs {
    input:
      vcfs=vcfcluster.clustered_vcf,
      batch=batch,
      algorithm=algorithm
  }

  call vcf_qc.master_vcf_qc as vcf_qc {
    input:
      vcf=concat_vcfs.vcf,
      famfile=trios_famfile,
      ref_build=ref_build,
      prefix="${batch}_clustered_${algorithm}_vcf",
      sv_per_shard=10000,
      samples_per_shard=100,
      Sanders_2015_tarball=Sanders_2015_tarball,
      Collins_2017_tarball=Collins_2017_tarball,
      Werling_2018_tarball=Werling_2018_tarball
  }

  output {
    File clustered_vcf = concat_vcfs.vcf
    File clustered_vcf_qc = vcf_qc.sv_vcf_qc_output
  }
}

task vcfcluster {
  Array[File] vcfs
  String batch
  String algorithm
  String chrom

  # VCFcluster parameters
  Int dist
  Float frac
  File blacklist
  Int svsize
  String svtypes
  String flags

  command <<<
    for f in ${sep=' ' vcfs}; do tabix -p vcf -f $f; done;
    tabix -p bed ${blacklist};

    svtk vcfcluster ${write_tsv(vcfs)} stdout \
        -r ${chrom} \
        -p ${batch}_${algorithm}_${chrom} \
        -d ${dist} \
        -f ${frac} \
        -x ${blacklist} \
        -z ${svsize} \
        -t ${svtypes} \
        ${flags} \
      | vcf-sort -c \
      | bgzip -c > ${batch}.${algorithm}.${chrom}.vcf.gz
  >>>

  output {
    File clustered_vcf="${batch}.${algorithm}.${chrom}.vcf.gz"
  }
  
  runtime {
    docker: "talkowski/sv-pipeline"
    disks: "local-disk 300 HDD"
    preemptible: 3
  }
}

task concat_vcfs {
  Array[File] vcfs
  String batch
  String algorithm

  command {
    vcf-concat ${sep=' ' vcfs} | vcf-sort -c | bgzip -c > ${batch}.${algorithm}.vcf.gz;
    tabix -p vcf ${batch}.${algorithm}.vcf.gz;
  }

  output {
    File vcf="${batch}.${algorithm}.vcf.gz"
    File idx="${batch}.${algorithm}.vcf.gz.tbi"
  }
  
  runtime {
    docker: "talkowski/sv-pipeline"
    disks: "local-disk 300 HDD"
    preemptible: 3
  }
}