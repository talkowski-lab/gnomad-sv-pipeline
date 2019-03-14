#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis helper script

# Summarize stats from VCF downsampling routine


###Set master parameters
options(stringsAsFactors=F,scipen=1000)
svtypes <- c("DEL","DUP","INS","INV","CPX","BND")
allpops <- c("AFR","EAS","EUR","AMR","SAS")


###Helper functions
#Summarize median # of singletons per sample
medianSingletonsPerSample <- function(singletons,svtypes,pops,PCRMINUS.samples){
  dtmp <- merge(singletons,pops,sort=F,by="sample")
  dtmp <- dtmp[which(dtmp$sample %in% PCRMINUS.samples),]
  res <- lapply(allpops,function(pop){
    meds <- sapply(svtypes,function(svtype){
      median(dtmp$count[which(dtmp$svtype==svtype & dtmp$pop==pop)],na.rm=T)
    })
    out.v <- c("samples"=length(unique(as.character(dtmp$sample[which(dtmp$pop==pop)]))),meds)
    names(out.v) <- paste(pop,names(out.v),sep=".")
    return(out.v)
  })
  return(t(as.data.frame(unlist(res))))
}

#Count total number of sites discovered by svtype
count.sites <- function(bed,svtypes){
  counts <- as.data.frame(t(sapply(svtypes,function(svtype){
    length(which(bed$SVTYPE==svtype))
  })))
  return(counts)
}

#Process SNV gene data
read.snvdata <- function(SNVdata.in,gene.metadata.in){
  #Read & clean
  snv.data <- read.table(SNVdata.in,header=T,comment.char="")
  metadata <- read.table(gene.metadata.in,header=T)
  merged <- merge(x=snv.data,y=metadata,by="gene",sort=F)
  merged <- merged[which(!(merged$chrom %in% c("chrX","chrY"))),]
  #Assign oe deciles
  merged$mis_oe_dec <- ceiling(10*rank(merged$mis_oe)/(nrow(merged)+1))
  merged$ptv_oe_dec <- ceiling(10*rank(merged$ptv_oe)/(nrow(merged)+1))
  #Assign oe percentiles
  merged$mis_oe_cent <- ceiling(100*rank(merged$mis_oe)/(nrow(merged)+1))
  merged$ptv_oe_cent <- ceiling(100*rank(merged$ptv_oe)/(nrow(merged)+1))
  #Return formatted data
  return(merged)
}

#Organize table of genes from SV data
getSVbed <- function(bed,genes,prefix=NULL){
  #LoF
  lof.any <- as.character(unlist(strsplit(bed$PROTEIN_CODING__LOF[which(!is.na(bed$PROTEIN_CODING__LOF))],split=",")))
  lof.del <- as.character(unlist(strsplit(bed$PROTEIN_CODING__LOF[which(!is.na(bed$PROTEIN_CODING__LOF)
                                                                        & bed$SVTYPE=="DEL")],split=",")))
  lof.other <- as.character(unlist(strsplit(bed$PROTEIN_CODING__LOF[which(!is.na(bed$PROTEIN_CODING__LOF)
                                                                          & bed$SVTYPE!="DEL")],split=",")))
  #Dup GC
  cg.dup <- as.character(unlist(strsplit(bed$PROTEIN_CODING__COPY_GAIN[which(!is.na(bed$PROTEIN_CODING__COPY_GAIN))],split=",")))
  #Dup LoF
  plof.dup <- as.character(unlist(strsplit(bed$PROTEIN_CODING__DUP_LOF[which(!is.na(bed$PROTEIN_CODING__DUP_LOF))],split=",")))
  #Inversion span
  inv.span <- as.character(unlist(strsplit(as.character(bed$PROTEIN_CODING__INV_SPAN[which(!is.na(bed$PROTEIN_CODING__INV_SPAN))]),split=",")))
  #Collect vector per gene
  res <- as.data.frame(t(sapply(genes,function(gene){
    g.lof.any <- length(which(lof.any==gene))
    g.lof.del <- length(which(lof.del==gene))
    g.lof.other <- length(which(lof.other==gene))
    g.cg.dup <- length(which(cg.dup==gene))
    g.plof.dup <- length(which(plof.dup==gene))
    g.inv.span <- length(which(inv.span==gene))
    g.out <- as.integer(c(g.lof.any,g.lof.del,g.lof.other,g.cg.dup,g.plof.dup,g.inv.span))
    g.out[which(is.na(g.out))] <- 0
    g.out <- c(gene,g.out)
    return(g.out)
  })))
  colnames(res) <- c("gene","lof.any","lof.del","lof.other","cg","plof","inv")
  if(!is.null(prefix)){
    colnames(res)[-1] <- paste(prefix,colnames(res)[-1],sep=".")
  }
  rownames(res) <- 1:nrow(res)
  return(res)
}


###Read command-line arguments
args <- commandArgs(trailingOnly=T)
bed.in <- as.character(args[1])
singletons.in <- as.character(args[2])
pop.assignments.in <- as.character(args[3])
PCRMINUS.samples.in <- as.character(args[4])
SNVdata.in <- as.character(args[5])
gene_metadata.in <- as.character(args[6])
OUTPREFIX <- as.character(args[7])

# bed.in <- "~/scratch/downsampled.vcf2bed.bed.gz"
# singletons.in <- "~/scratch/singletons_per_sample.txt"
# pop.assignments.in <- "~/scratch/gnomAD_v2_SV.sample_population_assignments.rough_initial_simplified.Nov19_2018_updated_Laurent_PCA.txt"
# PCRMINUS.samples.in <- "~/scratch/PCRMINUS_samples.list"
# SNVdata.in <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/gnomAD_v2.1_canonical_constraint.condensed.txt.gz"
# gene_metadata.in <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/Gencode.v19.autosomal_canonical_gene_metadata.txt.gz"
# OUTPREFIX <- "~/scratch/downsample_test"


###Read data
bed <- read.table(bed.in,header=T,comment.char="",sep="\t")
bed <- bed[which(bed$FILTER=="PASS" | bed$FILTER=="MULTIALLELIC"),]
singletons <- read.table(singletons.in,header=T,sep="\t")
pops <- read.table(pop.assignments.in,header=F,sep="\t")
colnames(pops) <- c("sample","pop")
PCRMINUS.samples <- as.character(read.table(PCRMINUS.samples.in,header=F)[,1])
gene.data <- read.snvdata(SNVdata.in,gene_metadata.in)
genes <- sort(unique(as.character(gene.data$gene)))


###Run analysis
#Analysis 1: get median count of singletons per sample
median.singletons <- medianSingletonsPerSample(singletons,svtypes,pops,PCRMINUS.samples)

#Analysis 2: get total count of sites in VCF by svtype
sites.count <- count.sites(bed,svtypes)

#Analysis 3: get total counts of variants per gene by predicted effect
sv.pergene <- getSVbed(bed,genes)


###Write out results
write.table(median.singletons,paste(OUTPREFIX,"median_singletons_per_sample.txt",sep="."),
            col.names=T,row.names=F,quote=F,sep="\t")
write.table(sites.count,paste(OUTPREFIX,"total_sv_sites.txt",sep="."),
            col.names=T,row.names=F,quote=F,sep="\t")
write.table(sv.pergene,paste(OUTPREFIX,"sv_per_gene.txt",sep="."),
            col.names=T,row.names=F,quote=F,sep="\t")

