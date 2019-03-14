#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Collect per-individual summary data for formal gnomAD analysis


###Set master parameters
options(stringsAsFactors=F,scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Read & clean vcf2bed input for basic counts
read.vcf2bed.forCounts <- function(vcf2bed.in,pass.only=F){
  #Read data
  dat <- read.table(vcf2bed.in,comment.char="",header=T,sep="\t")
  colnames(dat)[1] <- "chrom"
  #Restrict to sites with at least one observed alternative allele
  dat <- dat[union(which(dat$AC>0),grep("MULTIALLELIC",dat$FILTER,fixed=T)),]
  #Subset to passing variants for final analysis, if optioned
  if(pass.only==T){
    dat <- dat[which(dat$FILTER=="PASS" | dat$FILTER=="MULTIALLELIC"),]
  }
  #Drop columns not being used (to save memory)
  cols.to.drop <- c("start","end","ALGORITHMS","CHR2","CPX_INTERVALS",
                    "END","EVIDENCE","SOURCE","STRANDS","UNRESOLVED_TYPE",
                    colnames(dat)[grep("N_",colnames(dat),fixed=T)],
                    colnames(dat)[grep("PROTEIN_CODING__",colnames(dat),fixed=T)],
                    colnames(dat)[grep("LINCRNA__",colnames(dat),fixed=T)],
                    colnames(dat)[grep("NONCODING_",colnames(dat),fixed=T)])
  dat <- dat[,-which(colnames(dat) %in% cols.to.drop)]
  #Convert numeric columns
  numeric.columns <- sort(unique(c(grep("FREQ",colnames(dat),fixed=T),
                                   grep("AN",colnames(dat),fixed=T),
                                   grep("AC",colnames(dat),fixed=T),
                                   grep("AF",colnames(dat),fixed=T))))
  dat[,numeric.columns] <- apply(dat[,numeric.columns],2,as.numeric)
  return(dat)
}
#Read & clean vcf2bed input for functional data
read.vcf2bed.forFunc <- function(vcf2bed.in,pass.only=F){
  #Read data
  dat <- read.table(vcf2bed.in,comment.char="",header=T,sep="\t")
  colnames(dat)[1] <- "chrom"
  #Subset to passing variants for final analysis, if optioned
  if(pass.only==T){
    dat <- dat[which(dat$FILTER=="PASS" | dat$FILTER=="MULTIALLELIC"),]
  }
  #Drop columns not being used (to save memory)
  cols.to.drop <- c("start","end","ALGORITHMS","CHR2","CPX_INTERVALS",
                    "END","EVIDENCE","SOURCE","STRANDS","UNRESOLVED_TYPE",
                    colnames(dat)[grep("AFR_",colnames(dat),fixed=T)],
                    colnames(dat)[grep("ASN_",colnames(dat),fixed=T)],
                    colnames(dat)[grep("EUR_",colnames(dat),fixed=T)],
                    colnames(dat)[grep("HSP_",colnames(dat),fixed=T)],
                    colnames(dat)[grep("OTH_",colnames(dat),fixed=T)],
                    colnames(dat)[grep("LINCRNA__",colnames(dat),fixed=T)],
                    colnames(dat)[grep("NONCODING_",colnames(dat),fixed=T)])
  dat <- dat[,-which(colnames(dat) %in% cols.to.drop)]
  #Convert numeric columns
  numeric.columns <- sort(unique(c(grep("FREQ",colnames(dat),fixed=T),
                                   grep("AN",colnames(dat),fixed=T),
                                   grep("AC",colnames(dat),fixed=T),
                                   grep("AF",colnames(dat),fixed=T))))
  numeric.columns <- setdiff(numeric.columns,grep("SPAN",colnames(dat),fixed=T))
  dat[,numeric.columns] <- apply(dat[,numeric.columns],2,as.numeric)
  return(dat)
}
#Read a single sample's VID/genotype list and parse genoytpes
readSampleData <- function(ID,PERSAMPLEDIR){
  #Load sample/VID/genotype list, if it exists
  VID.list.in <- paste(PERSAMPLEDIR,"/",ID,".VIDs_genotypes.txt.gz",sep="")
  if(file.exists(VID.list.in)){
    VID.list <- read.table(VID.list.in,
                           header=F,colClasses=c("character","character","NULL"))
    colnames(VID.list) <- c("name","geno")
    VID.list <- VID.list[which(!(VID.list$geno %in% c("./.","./2","0/0"))),]
    #Parse genotypes
    mCNV.genos.idx <- grep("./",VID.list$geno,fixed=T)
    parsed.bi.genos <- sapply(VID.list$geno[-mCNV.genos.idx],function(gt){
      sum(as.numeric(unlist(strsplit(gt,split="/"))))
    })
    parsed.mCNV.genos <- sapply(VID.list$geno[mCNV.genos.idx],function(gt){
      as.numeric(unlist(strsplit(gt,split="/"))[2])-2
    })
    VID.list$geno[-mCNV.genos.idx] <- parsed.bi.genos
    VID.list$geno[mCNV.genos.idx] <- parsed.mCNV.genos
    colnames(VID.list)[2] <- "sample.AC"
  }else{
    VID.list <- NULL
  }
  return(VID.list)
}
#Summarize a single sample's data into a vector of counts
summarizeSampleVariants.forCounts <- function(ID,PERSAMPLEDIR,pop,dat,alleles=F){
  #Read data
  VID.list <- readSampleData(ID,PERSAMPLEDIR)
  if(!is.null(VID.list)){
    #Join sample VID list with master variant data
    sample.dat <- merge(VID.list,dat,on="name",sort=F)
    rm(VID.list)
    sample.dat$sample.AC <- as.numeric(sample.dat$sample.AC)
    #Convert MCNV to MCNV_DEL and MCNV_DUP based on sample.AC
    sample.dat$SVTYPE[which(sample.dat$SVTYPE=="MCNV" & sample.dat$sample.AC<0)] <- "MCNV_DEL"
    sample.dat$SVTYPE[which(sample.dat$SVTYPE=="MCNV" & sample.dat$sample.AC>0)] <- "MCNV_DUP"
    #Count variants per SVTYPE
    SVTYPES.for.perSample.counts <- c("DEL","MCNV_DEL","DUP","MCNV_DUP","INS","INV","CPX","BND","CTX")
    all.counts <- sapply(SVTYPES.for.perSample.counts,function(svtype){
      if(alleles==T){
        sum(abs(sample.dat$sample.AC[which(sample.dat$SVTYPE==svtype)]))
      }else{
        length(which(sample.dat$SVTYPE==svtype & sample.dat$sample.AC!=0))
      }
    })
    names(all.counts) <- paste("all",SVTYPES.for.perSample.counts,sep=".")
    all.counts <- c("all.SITES"=sum(all.counts),
                    "all.VARIANTS"=sum(all.counts[which(names(all.counts)!="all.BND")]),
                    all.counts)
    #Count rare variants per SVTYPE
    rare.counts <- sapply(SVTYPES.for.perSample.counts,function(svtype){
      if(alleles==T){
        sum(abs(sample.dat$sample.AC[which(sample.dat$SVTYPE==svtype & sample.dat$AF<0.01)]))
      }else{
        length(which(sample.dat$SVTYPE==svtype & sample.dat$sample.AC!=0 & sample.dat$AF<0.01))
      }
    })
    names(rare.counts) <- paste("rare",SVTYPES.for.perSample.counts,sep=".")
    rare.counts <- c("rare.SITES"=sum(rare.counts),
                    "rare.VARIANTS"=sum(rare.counts[which(names(rare.counts)!="rare.BND")]),
                    rare.counts)
    #Count singleton variants per SVTYPE
    singleton.counts <- sapply(SVTYPES.for.perSample.counts,function(svtype){
      if(alleles==T){
        sum(abs(sample.dat$sample.AC[which(sample.dat$SVTYPE==svtype & sample.dat$AC==1)]))
      }else{
        length(which(sample.dat$SVTYPE==svtype & sample.dat$sample.AC!=0 & sample.dat$AC==1))
      }
    })
    names(singleton.counts) <- paste("singleton",SVTYPES.for.perSample.counts,sep=".")
    singleton.counts <- c("singleton.SITES"=sum(singleton.counts),
                    "singleton.VARIANTS"=sum(singleton.counts[which(names(singleton.counts)!="singleton.BND")]),
                    singleton.counts)
    #Count large (>1Mb) variants per SVTYPE
    large.counts <- sapply(SVTYPES.for.perSample.counts,function(svtype){
      if(alleles==T){
        sum(abs(sample.dat$sample.AC[which(sample.dat$SVTYPE==svtype & sample.dat$SVLEN>=1000000)]))
      }else{
        length(which(sample.dat$SVTYPE==svtype & sample.dat$SVLEN>=1000000))
      }
    })
    names(large.counts) <- paste("large",SVTYPES.for.perSample.counts,sep=".")
    large.counts <- c("large.SITES"=sum(large.counts),
                    "large.VARIANTS"=sum(large.counts[which(names(large.counts)!="large.BND")]),
                    large.counts)
    #Count homozygous variants per SVTYPE
    homozygous.counts <- sapply(SVTYPES.for.perSample.counts,function(svtype){
      if(alleles==T){
        sum(abs(sample.dat$sample.AC[which(sample.dat$SVTYPE==svtype & abs(sample.dat$sample.AC)>1)]))
      }else{
        length(which(sample.dat$SVTYPE==svtype & abs(sample.dat$sample.AC)>1))
      }
    })
    names(homozygous.counts) <- paste("homozygous",SVTYPES.for.perSample.counts,sep=".")
    homozygous.counts <- c("homozygous.SITES"=sum(homozygous.counts),
                    "homozygous.VARIANTS"=sum(homozygous.counts[which(names(homozygous.counts)!="homozygous.BND")]),
                    homozygous.counts)
    #Prep output vector of counts & sample ID/pop
    out.vect <- c("sample"=ID,"pop"=pop,all.counts,rare.counts,singleton.counts,
             large.counts,homozygous.counts)
  }else{
    out.vect <- NULL
  }
  return(out.vect)
}
#Create a data frame of summary counts per sample
summarizeAllSamples.forCounts <- function(pop.assignments,PERSAMPLEDIR,dat,alleles=F){
  #Iterate over all samples
  res <- do.call("rbind", lapply(split(pop.assignments,seq(nrow(pop.assignments))),function(vals){
    ID <- as.character(vals[1])
    pop <- as.character(vals[2])
    summarizeSampleVariants.forCounts(ID,PERSAMPLEDIR,pop,dat,alleles=alleles)
  }))
  #Format data frame
  res <- as.data.frame(res)
  res[,-c(1:2)] <- apply(res[,-c(1:2)],2,as.numeric)
  return(res)
}
#Summarize counts of complex SV subtypes per sample
summarizeSampleComplexCounts <- function(ID,PERSAMPLEDIR,pop,dat,alleles=F){
  #Read data
  VID.list <- readSampleData(ID,PERSAMPLEDIR)
  if(!is.null(VID.list)){
    #Join sample VID list with master variant data for CPX only
    sample.dat <- merge(VID.list,dat[which(dat$SVTYPE=="CPX"),],on="name",sort=F)
    rm(VID.list)
    sample.dat$sample.AC <- as.numeric(sample.dat$sample.AC)
    #Count variants per CPXTYPE
    CPXTYPES.for.perSample.counts <- sort(unique(dat$CPX_TYPE))
    all.counts <- sapply(CPXTYPES.for.perSample.counts,function(cpxtype){
      if(alleles==T){
        sum(abs(sample.dat$sample.AC[which(sample.dat$CPX_TYPE==cpxtype)]))
      }else{
        length(which(sample.dat$CPX_TYPE==cpxtype & sample.dat$sample.AC!=0))
      }
    })
    names(all.counts) <- CPXTYPES.for.perSample.counts
    all.counts <- c("ALL_CPX"=sum(all.counts),
                    all.counts)
    #Prep output vector of counts & sample ID/pop
    out.vect <- c("sample"=ID,"pop"=pop,all.counts)
  }else{
    out.vect <- NULL
  }
  return(out.vect)
}
#Create a data frame of complex subtype counts per sample
summarizeAllSamples.forComplexCounts <- function(pop.assignments,PERSAMPLEDIR,dat,alleles=F){
  #Iterate over all samples
  res <- do.call("rbind", lapply(split(pop.assignments,seq(nrow(pop.assignments))),function(vals){
    ID <- as.character(vals[1])
    pop <- as.character(vals[2])
    summarizeSampleComplexCounts(ID,PERSAMPLEDIR,pop,dat,alleles=alleles)
  }))
  #Format data frame
  res <- as.data.frame(res)
  res[,-c(1:2)] <- apply(res[,-c(1:2)],2,as.numeric)
  return(res)
}
#Summarize a single sample's data into a vector of functional summary data
summarizeSampleVariants.forFunc <- function(ID,PERSAMPLEDIR,pop,dat,alleles=F){
  #Read data
  VID.list <- readSampleData(ID,PERSAMPLEDIR)
  if(!is.null(VID.list)){
    #Join sample VID list with master variant data
    sample.dat <- merge(VID.list,dat,on="name",sort=F)
    rm(VID.list)
    sample.dat$sample.AC <- as.numeric(sample.dat$sample.AC)
    #Convert MCNV to MCNV_DEL and MCNV_DUP based on sample.AC
    sample.dat$SVTYPE[which(sample.dat$SVTYPE=="MCNV" & sample.dat$sample.AC<0)] <- "MCNV_DEL"
    sample.dat$PROTEIN_CODING__MCNV_LOSS <- NA
    if(length(which(sample.dat$SVTYPE=="MCNV_DEL"))>0){
      sample.dat$PROTEIN_CODING__MCNV_LOSS[which(sample.dat$SVTYPE=="MCNV_DEL")] <- sample.dat$PROTEIN_CODING__MSV_EXON_OVR[which(sample.dat$SVTYPE=="MCNV_DEL")]
    }
    sample.dat$SVTYPE[which(sample.dat$SVTYPE=="MCNV" & sample.dat$sample.AC>0)] <- "MCNV_DUP"
    sample.dat$PROTEIN_CODING__MCNV_GAIN <- NA
    if(length(which(sample.dat$SVTYPE=="MCNV_DUP"))>0){
      sample.dat$PROTEIN_CODING__MCNV_GAIN[which(sample.dat$SVTYPE=="MCNV_DUP")] <- sample.dat$PROTEIN_CODING__MSV_EXON_OVR[which(sample.dat$SVTYPE=="MCNV_DUP")]
    }
    
    #Count variants per functional class
    func.for.perSample.counts <- c("PROTEIN_CODING__LOF","PROTEIN_CODING__DUP_LOF",
                                   "PROTEIN_CODING__COPY_GAIN","PROTEIN_CODING__DUP_PARTIAL",
                                   "PROTEIN_CODING__MCNV_LOSS","PROTEIN_CODING__MCNV_GAIN",
                                   "PROTEIN_CODING__INV_SPAN","PROTEIN_CODING__INTRONIC",
                                   "PROTEIN_CODING__PROMOTER","PROTEIN_CODING__UTR")
    #Count all variants per functional consequence
    all.variants <- sapply(func.for.perSample.counts,function(func){
      length(which(!is.na(sample.dat[which(colnames(sample.dat)==func)]) & sample.dat$sample.AC!=0))
    })
    names(all.variants) <- paste("all",func.for.perSample.counts,sep=".")
    #Count all unique genes per functional consequence
    all.genes <- sapply(func.for.perSample.counts,function(func){
      length(unique(as.character(unlist(strsplit(sample.dat[,which(colnames(sample.dat)==func)][which(!is.na(sample.dat[which(colnames(sample.dat)==func)]) & sample.dat$sample.AC!=0)],split=",")))))
    })
    names(all.genes) <- paste("all",func.for.perSample.counts,"genes",sep=".")
    #Count rare variants per functional consequence
    rare.variants <- sapply(func.for.perSample.counts,function(func){
      length(which(!is.na(sample.dat[which(colnames(sample.dat)==func)]) & sample.dat$sample.AC!=0 & sample.dat$AF<0.01))
    })
    names(rare.variants) <- paste("rare",func.for.perSample.counts,sep=".")
    #Count rare unique genes per functional consequence
    rare.genes <- sapply(func.for.perSample.counts,function(func){
      length(unique(as.character(unlist(strsplit(sample.dat[,which(colnames(sample.dat)==func)][which(!is.na(sample.dat[which(colnames(sample.dat)==func)]) & sample.dat$sample.AC!=0 & sample.dat$AF<0.01)],split=",")))))
    })
    names(rare.genes) <- paste("rare",func.for.perSample.counts,"genes",sep=".")
    #Count singleton variants per functional consequence
    singleton.variants <- sapply(func.for.perSample.counts,function(func){
      length(which(!is.na(sample.dat[which(colnames(sample.dat)==func)]) & sample.dat$sample.AC!=0 & sample.dat$AC==1))
    })
    names(singleton.variants) <- paste("singleton",func.for.perSample.counts,sep=".")
    #Count singleton unique genes per functional consequence
    singleton.genes <- sapply(func.for.perSample.counts,function(func){
      length(unique(as.character(unlist(strsplit(sample.dat[,which(colnames(sample.dat)==func)][which(!is.na(sample.dat[which(colnames(sample.dat)==func)]) & sample.dat$sample.AC!=0 & sample.dat$AC==1)],split=",")))))
    })
    names(singleton.genes) <- paste("singleton",func.for.perSample.counts,"genes",sep=".")
    #Count homozygous variants per functional consequence
    homozygous.variants <- sapply(func.for.perSample.counts,function(func){
      length(which(!is.na(sample.dat[which(colnames(sample.dat)==func)]) & abs(sample.dat$sample.AC)>1))
    })
    names(homozygous.variants) <- paste("homozygous",func.for.perSample.counts,sep=".")
    #Count homozygous unique genes per functional consequence
    homozygous.genes <- sapply(func.for.perSample.counts,function(func){
      length(unique(as.character(unlist(strsplit(sample.dat[,which(colnames(sample.dat)==func)][which(!is.na(sample.dat[which(colnames(sample.dat)==func)]) & abs(sample.dat$sample.AC)>1)],split=",")))))
    })
    names(homozygous.genes) <- paste("homozygous",func.for.perSample.counts,"genes",sep=".")
    #Count LoF variants by SVTYPE
    svtypes.for.lof <- c("DEL","DUP","MCNV_DEL","INS","INV","CPX")
    lof.variants.bytype <- sapply(svtypes.for.lof,function(svtype){
      length(which(sample.dat$SVTYPE==svtype & (!is.na(sample.dat$PROTEIN_CODING__LOF) | !is.na(sample.dat$PROTEIN_CODING__MCNV_LOSS))))
    })
    names(lof.variants.bytype) <- paste("lof",svtypes.for.lof,"variants",sep=".")
    #Count genes with LoF variants by SVTYPE
    lof.genes.bytype <- sapply(svtypes.for.lof,function(svtype){
      keep.idx <- which(sample.dat$SVTYPE==svtype & (!is.na(sample.dat$PROTEIN_CODING__LOF) | !is.na(sample.dat$PROTEIN_CODING__MCNV_LOSS)))
      length(sort(unique(as.character(unlist(strsplit(c(sample.dat$PROTEIN_CODING__LOF[keep.idx],
                                                        sample.dat$PROTEIN_CODING__MCNV_LOSS[keep.idx]),split=","))))))
    })
    names(lof.genes.bytype) <- paste("lof",svtypes.for.lof,"genes",sep=".")
    #Prep output vector of counts & sample ID/pop
    out.vect <- c("sample"=ID,"pop"=pop,all.variants,all.genes,rare.variants,rare.genes,
                  singleton.variants,singleton.genes,homozygous.variants,homozygous.genes,
                  lof.variants.bytype,lof.genes.bytype)
  }else{
    out.vect <- NULL
  }
  return(out.vect)
}
#Create a data frame of summary counts per sample
summarizeAllSamples.forFunc <- function(pop.assignments,PERSAMPLEDIR,dat){
  #Iterate over all samples
  res <- do.call("rbind", lapply(split(pop.assignments,seq(nrow(pop.assignments))),function(vals){
    ID <- as.character(vals[1])
    pop <- as.character(vals[2])
    summarizeSampleVariants.forFunc(ID,PERSAMPLEDIR,pop,dat)
  }))
  #Format data frame
  res <- as.data.frame(res)
  res[,-c(1:2)] <- apply(res[,-c(1:2)],2,as.numeric)
  return(res)
}


################
###RSCRIPT BLOCK
################
require(optparse,quietly=T)
###List of command-line options
option_list <- list(
  make_option(c("--prefix"), type="character", default="gnomAD_v2_SV",
              help="prefix used for naming outfiles [default %default]",
              metavar="character")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog VCF2BED PERSAMPLEDIR SAMPLE_POP_ASSIGNMENTS OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
vcf2bed.in <- args$args[1]
PERSAMPLEDIR <- args$args[2]
pop.assignments.in <- args$args[3]
OUTDIR <- args$args[4]
prefix <- opts$prefix

# #Dev parameters (local)
# vcf2bed.in <- "~/scratch/gnomAD_v2_SV_MASTER.vcf2bed.bed.gz"
# PERSAMPLEDIR <- "~/scratch/gnomAD_v2_SV_MASTER_perSample_VIDs_merged/"
# #pop.assignments.in <- "~/scratch/gnomAD_v2_SV_MASTER.PCA_population_assignments.txt"
# pop.assignments.in <- "~/scratch/gnomAD_v2_SV_MASTER.pop_assignment_shard_00001"
# OUTDIR <- "~/scratch/gnomAD_v2_test_plots/"
# prefix <- "gnomAD_v2_SV_MASTER"

###Create output directory, if necessary
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

###Process input data
pop.assignments <- read.table(pop.assignments.in,header=F)
colnames(pop.assignments) <- c("sample","pop")
pop.assignments$pop[which(pop.assignments$pop==".")] <- "OTH"

###Gather site counts
cat("NOW SUMMARIZING SV COUNTS PER SAMPLE\n")
dat.forCounts <- read.vcf2bed.forCounts(vcf2bed.in)
summary.forCounts <- summarizeAllSamples.forCounts(pop.assignments,PERSAMPLEDIR,dat.forCounts)
write.table(summary.forCounts,paste(OUTDIR,"/",prefix,".sitesPerSample_summary.txt",sep=""),
            col.names=T,row.names=F,sep="\t",quote=F)

###Gather allele counts
cat("NOW SUMMARIZING SV COUNTS PER SAMPLE\n")
summary.forAlleleCounts <- summarizeAllSamples.forCounts(pop.assignments,PERSAMPLEDIR,dat.forCounts,alleles=T)
write.table(summary.forAlleleCounts,paste(OUTDIR,"/",prefix,".allelesPerSample_summary.txt",sep=""),
            col.names=T,row.names=F,sep="\t",quote=F)

###Gather count of complex subtypes per sample
cat("NOW SUMMARIZING COMPLEX SV COUNTS PER SAMPLE\n")
summary.forComplex <- summarizeAllSamples.forComplexCounts(pop.assignments,PERSAMPLEDIR,dat.forCounts)
write.table(summary.forComplex,paste(OUTDIR,"/",prefix,".complexCountsPerSample_summary.txt",sep=""),
            col.names=T,row.names=F,sep="\t",quote=F)

###Gather functional counts
cat("NOW SUMMARIZING FUNCTIONAL DATA PER SAMPLE\n")
dat.forFunc <- read.vcf2bed.forFunc(vcf2bed.in)
summary.forFunc <- summarizeAllSamples.forFunc(pop.assignments,PERSAMPLEDIR,dat.forFunc)
write.table(summary.forFunc,paste(OUTDIR,"/",prefix,".functionalPerSample_summary.txt",sep=""),
            col.names=T,row.names=F,sep="\t",quote=F)


