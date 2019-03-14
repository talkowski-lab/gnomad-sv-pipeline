import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:01_depth_clustering_by_chrom/versions/4/plain-WDL/descriptor" as dibc
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:master_SV_VCF_QC/versions/47/plain-WDL/descriptor" as vcf_qc

workflow cluster_depth {
  File del_bed
  File dup_bed
  File contigs
  Float frac
  String flags
  String batch
  File famfile
  File trios_famfile
  String ref_build
  File Sanders_2015_tarball
  File Collins_2017_tarball
  File Werling_2018_tarball

  call dibc.bedcluster_by_chrom as cluster_DELs {
    input:
      batch=batch,
      svtype="DEL",
      bed=del_bed,
      contigs=contigs,
      frac=frac,
      flags=flags
  }

  call dibc.bedcluster_by_chrom as cluster_DUPs {
    input:
      batch=batch,
      svtype="DUP",
      bed=dup_bed,
      contigs=contigs,
      frac=frac,
      flags=flags
  }

  call make_rdtest_bed {
    input:
      dels=cluster_DELs.clustered_bed,
      dups=cluster_DUPs.clustered_bed,
      batch=batch,
  }

  call make_depth_vcf {
    input:
      bed=make_rdtest_bed.bed,
      batch=batch,
      contigs=contigs
  }

  call vcf_qc.master_vcf_qc as vcf_qc {
    input:
      vcf=make_depth_vcf.vcf,
      famfile=trios_famfile,
      ref_build=ref_build,
      prefix="${batch}_clustered_depth_vcf",
      sv_per_shard=10000,
      samples_per_shard=100,
      Sanders_2015_tarball=Sanders_2015_tarball,
      Collins_2017_tarball=Collins_2017_tarball,
      Werling_2018_tarball=Werling_2018_tarball
  }

  output {
    File clustered_vcf = make_depth_vcf.vcf
    File clustered_vcf_qc = vcf_qc.sv_vcf_qc_output
  }
}

task make_rdtest_bed { 
  File dels
  File dups
  File script
  String batch

  command <<<
    cat \
        <(python3 ${script} ${dels} | sed '1d') \
        <(python3 ${script} ${dups} | sed '1d') \
      | sort -k1,1V -k2,2n \
      | cat <(echo -e "#chrom start end name samples svtype" | sed -e 's/ /\t/g') - \
      > ${batch}.depth.bed;
  >>>
  
  output {
    File bed = "${batch}.depth.bed"
  }
  
  runtime {
    docker: "talkowski/sv-pipeline@sha256:a89824ac34b915f605d09bcf57516bc76d950bd762ad5c1f336d421be917be55"
    preemptible: 3
  }
}

task make_depth_vcf {
  File bed
  File contigs
  String batch
  
  command <<<
    cut -f5 ${bed} | sed -e '1d' -e 's/,/\n/g' | sort -u > samples.list;
    svtk rdtest2vcf --contigs ${contigs} ${bed} samples.list ${batch}.depth.vcf.gz;
  >>>

  output {
    File vcf = "${batch}.depth.vcf.gz"
  }
  
  runtime {
    docker: "talkowski/sv-pipeline@sha256:a89824ac34b915f605d09bcf57516bc76d950bd762ad5c1f336d421be917be55"
    preemptible: 3
  }
}