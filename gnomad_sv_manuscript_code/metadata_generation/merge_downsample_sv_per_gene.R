#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis helper script

# Get average SV per gene from downsampling experiments


###Set master parameters & load libraries
options(stringsAsFactors=F,scipen=1000)
effects <- c("lof.any","lof.del","lof.other","cg","plof","inv")


###Helper functions
#Import list of sv-per-gene counts, and take mean across downsampling points
process.counts <- function(counts.list.in,seed.table.in){
  #Read seeds
  down.sizes <- as.integer(read.table(seed.table.in,header=F,sep="\t")[,1])
  #Read tables as list
  counts.list <- as.character(read.table(counts.list.in,header=F,sep="\t")[,1])
  counts <- lapply(counts.list,function(l){
    read.table(l,header=T,sep="\t")
  })
  #Sanity check to make sure all counts have same number of rows
  if(length(unique(unlist(lapply(counts,nrow)))) > 1){
    stop("Some sv_per_gene.txt tables don't have the same number of lines.")
  }
  #Take mean counts per gene per downsample size
  down.points <- sort(unique(down.sizes))
  res <- do.call("cbind", lapply(down.points,function(n){
    idxs <- which(down.sizes==n)
    merged <- do.call("rbind", counts[idxs])
    meaned <- as.data.frame(t(sapply(sort(unique(as.character(merged$gene))),function(gene){
      as.numeric(apply(merged[which(merged$gene==gene),-1],2,mean,na.rm=T))
    })))
    colnames(meaned) <- paste(effects,n,sep=".")
    return(meaned)
  }))
  #Add gene name & return output data frame
  res <- data.frame("gene"=rownames(res),res)
  rownames(res) <- NULL
  return(res)
}


###Read command-line arguments
args <- commandArgs(trailingOnly=T)
counts.list.in <- as.character(args[1])
seed.table.in <- as.character(args[2])
OUTFILE <- as.character(args[3])

# #Dev parameters (local)
# counts.list.in <- "~/scratch/sv_per_gene.input.list"
# seed.table.in <- "~/scratch/tmp_seeds_input.txt"
# OUTFILE <- "~/scratch/merged_sv_per_gene.test.txt"


###Process data & write out
res <- process.counts(counts.list.in,seed.table.in)
write.table(res,OUTFILE,col.names=T,row.names=F,sep="\t",quote=F)

