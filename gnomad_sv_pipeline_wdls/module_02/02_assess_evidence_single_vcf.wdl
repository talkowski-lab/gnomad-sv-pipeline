import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:02_petest/versions/10/plain-WDL/descriptor" as petest
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:02_srtest/versions/12/plain-WDL/descriptor" as srtest
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:02_rdtest/versions/7/plain-WDL/descriptor" as rdtest
import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:02_baftest/versions/10/plain-WDL/descriptor" as baftest

# Assess PE, SR, RD, and BAF evidence for every variant in a VCF
workflow assess_evidence {
  File vcf                      # Input VCF
  File discfile                 # Discordant pair file
  File discfile_idx             # Tabix index of discordant pair file
  File splitfile                # Split read file
  File splitfile_idx            # Tabix index of split read file
  File coveragefile             # Bincov matrix
  File coveragefile_idx         # Tabix index of bincov matrix
  File medianfile               # Median coverage of each sample
  File baf_metrics              # Matrix of BAF statistics
  File baf_metrics_idx          # Tabix index of BAF matrix
  File famfile                  # Batch fam file 
  File autosome_contigs         # Autosomes .fai
  File allosome_contigs         # Allosomes .fai
  File rmsk                     # Repeatmasker track
  File segdups                  # Seg dups track
  String batch                  # Batch ID
  String algorithm              # Algorithm ID
  Int PE_split_size             # Number of lines in each petest split
  Int SR_split_size             # Number of lines in each srtest split
  Int RD_split_size             # Number of lines in each rdtest split
  Int BAF_split_size            # Number of lines in each baftest split
  File svc_acct_key
  Array[String] samples

  if (algorithm == "depth") {
    call rdtest.rdtest_by_chrom as rdtest_depth {
      input:
        vcf=vcf,
        svc_acct_key=svc_acct_key,
        coveragefile=coveragefile,
        coveragefile_idx=coveragefile_idx,
        medianfile=medianfile,
        famfile=famfile,
        autosome_contigs=autosome_contigs,  
        allosome_contigs=allosome_contigs,
        batch=batch,
        algorithm=algorithm,
        split_size=RD_split_size
    }
    call baftest.baftest_by_chrom as baftest_depth {
      input:
        vcf=vcf,
        samples=samples,
        svc_acct_key=svc_acct_key,
        baf_metrics=baf_metrics,
        baf_metrics_idx=baf_metrics_idx,
        autosome_contigs=autosome_contigs,
        batch=batch,
        algorithm=algorithm,
        split_size=BAF_split_size
    }
    
    call aggregate_depth {
      input:
        vcf=vcf,
        rdtest=rdtest_depth.rdtest,
        baftest=baftest_depth.baftest,
        segdups=segdups,
        rmsk=rmsk
    }
  }
  if (algorithm == "melt") {
    call srtest.srtest_by_chrom as srtest_melt {
      input:
        vcf=vcf,
        svc_acct_key=svc_acct_key,
        splitfile=splitfile,
        medianfile=medianfile,
        splitfile_idx=splitfile_idx,
        famfile=famfile,
        autosome_contigs=autosome_contigs,
        allosome_contigs=allosome_contigs,
        batch=batch,
        algorithm=algorithm,
        split_size=SR_split_size
    }
    
    call aggregate_melt {
      input:
        vcf=vcf,
        srtest=srtest_melt.srtest,
        segdups=segdups,
        rmsk=rmsk
    }
  }
  if ((algorithm != "melt") && (algorithm !="depth")) {
    call petest.petest_by_chrom {
      input:
        vcf=vcf,
        svc_acct_key=svc_acct_key,
        discfile=discfile,
        medianfile=medianfile,
        discfile_idx=discfile_idx,
        famfile=famfile,
        autosome_contigs=autosome_contigs,
        allosome_contigs=allosome_contigs,
        batch=batch,
        algorithm=algorithm,
        split_size=PE_split_size
    }

    call srtest.srtest_by_chrom {
      input:
        vcf=vcf,
        svc_acct_key=svc_acct_key,
        splitfile=splitfile,
        medianfile=medianfile,
        splitfile_idx=splitfile_idx,
        famfile=famfile,
        autosome_contigs=autosome_contigs,
        allosome_contigs=allosome_contigs,
        batch=batch,
        algorithm=algorithm,
        split_size=SR_split_size
    }
    call rdtest.rdtest_by_chrom as rdtest_pesr {
      input:
        vcf=vcf,
        svc_acct_key=svc_acct_key,
        coveragefile=coveragefile,
        coveragefile_idx=coveragefile_idx,
        medianfile=medianfile,
        famfile=famfile,
        autosome_contigs=autosome_contigs,  
        allosome_contigs=allosome_contigs,
        batch=batch,
        algorithm=algorithm,
        split_size=RD_split_size
    }
    call baftest.baftest_by_chrom as baftest_pesr {
      input:
        vcf=vcf,
        samples=samples,
		svc_acct_key=svc_acct_key,
        baf_metrics=baf_metrics,
        baf_metrics_idx=baf_metrics_idx,
        autosome_contigs=autosome_contigs,
        batch=batch,
        algorithm=algorithm,
        split_size=BAF_split_size
    }
    
    call aggregate_pesr {
      input:
        vcf=vcf,
        petest=petest_by_chrom.petest,
        srtest=srtest_by_chrom.srtest,
        rdtest=rdtest_pesr.rdtest,
        baftest=baftest_pesr.baftest,
        segdups=segdups,
        rmsk=rmsk
    }
  }

  Array[File?] outs = [aggregate_pesr.metrics, aggregate_depth.metrics,aggregate_melt.metrics,]

  output {
    File metrics = select_first(outs)
  }
}

task aggregate_depth {
  File vcf
  File rdtest
  File baftest
  File segdups
  File rmsk

  command <<<
    /opt/sv-pipeline/02_evidence_assessment/02e_metric_aggregation/scripts/aggregate.py -v ${vcf} \
      -r ${rdtest} -b ${baftest} \
      --segdups ${segdups} --rmsk ${rmsk} \
      depth.metrics
  >>>

  output {
    File metrics = "depth.metrics"
  }
  runtime {
  	preemptible: 3
    docker: "talkowski/sv-pipeline-remote-pysam"
    memory: "8 GB"
    disks: "local-disk 100 HDD"
  }
}

task aggregate_pesr {
  File vcf
  File petest
  File srtest
  File rdtest
  File baftest
  File segdups
  File rmsk

  command <<<
    /opt/sv-pipeline/02_evidence_assessment/02e_metric_aggregation/scripts/aggregate.py -v ${vcf} \
      -p ${petest} -s ${srtest} -r ${rdtest} -b ${baftest} \
      --segdups ${segdups} --rmsk ${rmsk} \
      pesr.metrics
  >>>

  output {
    File metrics = "pesr.metrics"
  }
  runtime {
  	preemptible: 3
    docker: "talkowski/sv-pipeline-remote-pysam"
    memory: "20 GB"
    disks: "local-disk 100 HDD"
  }
}
task aggregate_melt {
  File vcf
  File srtest
  File segdups
  File rmsk

  command <<<
    /opt/sv-pipeline/02_evidence_assessment/02e_metric_aggregation/scripts/aggregate.py -v ${vcf} \
      -s ${srtest}\
      --segdups ${segdups} --rmsk ${rmsk} \
      melt.metrics
  >>>

  output {
    File metrics = "melt.metrics"
  }
  runtime {
  	preemptible: 3
    docker: "talkowski/sv-pipeline-remote-pysam"
    memory: "8 GB"
    disks: "local-disk 100 HDD"
  }
}