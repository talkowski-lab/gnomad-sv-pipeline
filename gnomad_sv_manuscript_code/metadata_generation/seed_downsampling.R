#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis helper script

# Create list of unique downsample seeds


###Set master parameters
options(stringsAsFactors=F,scipen=1000)
sample.sizes <- c(1,2,3,4,5,6,7,8,9,10,25,50,75,100,250,500,750,1000,2500,5000,7500,10000)
sample.sizes <- sample.sizes[which(sample.sizes<=10000)]

###Read command-line arguments
args <- commandArgs(trailingOnly=T)
N <- as.numeric(args[1])
# master.seed <- as.numeric(args[2])
OUTFILE <- as.character(args[2])

###Create data frame of random seeds & sample sizes
s <- as.numeric(sapply(sample.sizes,rep,times=N))
set.seed(123456789)
r <- ceiling(runif(length(s), 0, 10^12))
out <- data.frame("sample.size"=s,
                  "random.seed"=r)
write.table(out,OUTFILE,col.names=F,row.names=F,sep="\t",quote=F)
