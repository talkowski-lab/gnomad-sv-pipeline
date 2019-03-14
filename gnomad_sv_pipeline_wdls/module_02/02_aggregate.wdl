import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:02_assess_evidence_single_vcf/versions/31/plain-WDL/descriptor"  as assess
workflow assess_evidence_batch {
  File mantavcf                      # Input VCF
  File meltvcf                      # Input VCF
  File dellyvcf                      # Input VCF
  File depthvcf                      # Input VCF
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
  Int PE_split_size             # Number of lines in each petest split
  Int SR_split_size             # Number of lines in each srtest split
  Int RD_split_size             # Number of lines in each rdtest split
  Int BAF_split_size            # Number of lines in each baftest split
  File svc_acct_key
  Array[String] samples
    call assess.assess_evidence as assessmanta{input:
       vcf=mantavcf,
       samples=samples,
       svc_acct_key=svc_acct_key,
       discfile=discfile,
       discfile_idx=discfile_idx,
       splitfile=splitfile,
       splitfile_idx=splitfile_idx,
       coveragefile=coveragefile,
       coveragefile_idx=coveragefile_idx,
       medianfile=medianfile,
       baf_metrics=baf_metrics,
       baf_metrics_idx=baf_metrics_idx,
       famfile=famfile,
       autosome_contigs=autosome_contigs,
       allosome_contigs=allosome_contigs,
       rmsk=rmsk,
       segdups=segdups,
       batch=batch,
       algorithm="manta",
       PE_split_size=PE_split_size,
       SR_split_size=SR_split_size,
       RD_split_size=RD_split_size,
       BAF_split_size=BAF_split_size,
    }
    call assess.assess_evidence as assessmelt{input:
       vcf=meltvcf,
       samples=samples,
       svc_acct_key=svc_acct_key,
       discfile=discfile,
       discfile_idx=discfile_idx,
       splitfile=splitfile,
       splitfile_idx=splitfile_idx,
       coveragefile=coveragefile,
       coveragefile_idx=coveragefile_idx,
       medianfile=medianfile,
       baf_metrics=baf_metrics,
       baf_metrics_idx=baf_metrics_idx,
       famfile=famfile,
       autosome_contigs=autosome_contigs,
       allosome_contigs=allosome_contigs,
       rmsk=rmsk,
       segdups=segdups,
       batch=batch,
       algorithm="melt",
       PE_split_size=PE_split_size,
       SR_split_size=SR_split_size,
       RD_split_size=RD_split_size,
       BAF_split_size=BAF_split_size,}
    call assess.assess_evidence as assessdelly{input:
       vcf=dellyvcf,
       samples=samples,
       svc_acct_key=svc_acct_key,
       discfile=discfile,
       discfile_idx=discfile_idx,
       splitfile=splitfile,
       splitfile_idx=splitfile_idx,
       coveragefile=coveragefile,
       coveragefile_idx=coveragefile_idx,
       medianfile=medianfile,
       baf_metrics=baf_metrics,
       baf_metrics_idx=baf_metrics_idx,
       famfile=famfile,
       autosome_contigs=autosome_contigs,
       allosome_contigs=allosome_contigs,
       rmsk=rmsk,
       segdups=segdups,
       batch=batch,
       algorithm="delly",
       PE_split_size=PE_split_size,
       SR_split_size=SR_split_size,
       RD_split_size=RD_split_size,
       BAF_split_size=BAF_split_size,}
    call assess.assess_evidence as assessdepth{input:
       vcf=depthvcf,
       samples=samples,
       svc_acct_key=svc_acct_key,
       discfile=discfile,
       discfile_idx=discfile_idx,
       splitfile=splitfile,
       splitfile_idx=splitfile_idx,
       coveragefile=coveragefile,
       coveragefile_idx=coveragefile_idx,
       medianfile=medianfile,
       baf_metrics=baf_metrics,
       baf_metrics_idx=baf_metrics_idx,
       famfile=famfile,
       autosome_contigs=autosome_contigs,
       allosome_contigs=allosome_contigs,
       rmsk=rmsk,
       segdups=segdups,
       batch=batch,
       algorithm="depth",
       PE_split_size=PE_split_size,
       SR_split_size=SR_split_size,
       RD_split_size=RD_split_size,
       BAF_split_size=BAF_split_size,}
    call aggregate_metric{input:
        batch=batch,mantametric=assessmanta.metrics,dellymetric=assessdelly.metrics,meltmetric=assessmelt.metrics,depthmetric=assessdepth.metrics}
    output{
        File metrics=aggregate_metric.metrics
    }
}
task aggregate_metric{
    String batch
    File mantametric
    File dellymetric
    File depthmetric
    File meltmetric
    command <<<
        python3 <<CODE
        import pandas as pd
        metrics = ["${mantametric}","${dellymetric}","${depthmetric}","${meltmetric}"]
        dfs=[]
        for df in metrics:
            dfs.append(pd.read_table(df))
        df = pd.concat(dfs)
        df.to_csv("${batch}.metrics", index=False, sep='\t')
        CODE
        >>>
    output{
        File metrics="${batch}.metrics"
    }
    runtime {
  	preemptible: 3
      docker: "talkowski/sv-pipeline-remote-pysam"
      memory: "20 GB"
      disks: "local-disk 100 HDD"
    }
}