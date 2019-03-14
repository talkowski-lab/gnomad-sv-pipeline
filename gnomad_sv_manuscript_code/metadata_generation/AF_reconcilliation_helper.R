#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Helper script to sanitize AFs for 1kG comparisons


###Set master parameters
options(stringsAsFactors=F,scipen=1000)

###Read command line args
args <- commandArgs(trailingOnly=TRUE)

###Read data
dat <- read.table(args[1],header=T,comment.char="",sep="\t")

###Clean AFs
AFs.to.clean <- grep(",",dat$AF,fixed=T)
cleaned.AFs <- as.vector(as.numeric(sapply(dat$AF[AFs.to.clean],function(s){
  return(sum(round(as.numeric(unlist(strsplit(s,split=","))),10)[-2]))
})))
cleaned.AFs[which(cleaned.AFs>1)] <- 1
dat$AF[AFs.to.clean] <- cleaned.AFs

###Write data
colnames(dat)[1] <- "#chr"
write.table(dat,args[1],col.names=T,row.names=F,sep="\t",quote=F)
