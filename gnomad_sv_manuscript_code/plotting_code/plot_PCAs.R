#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Plot 3D PCA and supplementary PCA grids for formal gnomAD analysis


###Set master parameters
options(stringsAsFactors=F,scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Read input data and consolidate labels
import.data <- function(PCA.in,pop.assignments.in,PCRPLUS.samples.in,
                        batch.assignments.in,pops,keep.n.pcs=20){
  #Read raw data
  pca <- read.table(PCA.in,header=T,comment.char="",sep="\t")[,1:(keep.n.pcs+1)]
  colnames(pca)[1] <- "sample"
  pop.assignments <- read.table(pop.assignments.in,header=F,comment.char="",sep="\t")
  colnames(pop.assignments) <- c("sample","pop")
  plus.samples.list <- read.table(PCRPLUS.samples.in,header=F)[,1]
  pcr.assignments <- data.frame("sample"=pca$sample,"pcr"="PCRMINUS")
  pcr.assignments$pcr[which(pcr.assignments$sample %in% plus.samples.list)] <- "PCRPLUS"
  batch.assignments <- read.table(batch.assignments.in,header=F,comment.char="",sep="\t")
  colnames(batch.assignments) <- c("sample","batch")
  #Merge data
  merged <- merge(x=pop.assignments,y=pcr.assignments,by="sample",sort=F,all=F)
  merged <- merge(x=merged,y=batch.assignments,by="sample",sort=F,all=F)
  merged <- merge(x=merged,y=pca,by="sample",sort=F,all=F)
  return(merged)
}


#####################
###PLOTTING FUNCTIONS
#####################
#3D scatterplot of PC1-3, colored by ancestry
scatter.3d <- function(dat,pops){
  #Get ancestry colors
  pop.cols <- as.character(sapply(dat$sample,function(ID){pops$color[which(pops$pop==as.character(dat$pop[which(dat$sample==ID)]))]}))
  # #Helper function to project small dots onto the background panels
  # panelfirst <- function(pmat){
  #   x <- dat$PC1
  #   y <- dat$PC2
  #   z <- dat$PC3
  #   XY <- trans3D(x, y, z = rep(min(z), length(z)), pmat = pmat)
  #   scatter2D(XY$x, XY$y, colvar = NULL, pch=19,
  #             cex = 0.05, add = TRUE, col=adjustcolor(pop.cols,alpha=0.2))
  #   
  #   YZ <- trans3D(x = rep(max(x), length(x)), y, z, pmat = pmat)
  #   scatter2D(YZ$x, YZ$y, colvar = NULL, pch=19,
  #             cex = 0.05, add = TRUE, col=adjustcolor(pop.cols,alpha=0.2))
  # }
  #Prep 3D plot
  par(mar=c(2,0,0,0))
  scatter3D(x=dat$PC1,y=dat$PC2,z=dat$PC4,pch=19,cex=0.35,colvar=NULL,
            col=adjustcolor(pop.cols,alpha=0.4),theta=305,phi=25,bty="b2",
            xlab="PC 1",
            ylab="PC 2",
            zlab="PC 3",
            ticktype="detailed")
}

#2D scatterplot, PCs 1v2
scatter.2d <- function(dat,pops){
  #Get ancestry colors
  pop.cols <- as.character(sapply(dat$sample,function(ID){pops$color[which(pops$pop==as.character(dat$pop[which(dat$sample==ID)]))]}))
  #Prep plot area
  par(mar=c(2.5,2.5,0.5,0.5),bty="n")
  plot(x=range(dat$PC1,na.rm=T),
       y=range(dat$PC2,na.rm=T),
       type="n",xaxt="n",xlab="",yaxt="n",ylab="")
  #Add points
  points(dat$PC1,dat$PC2,pch=21,cex=0.2,col=pop.cols,bg=adjustcolor(pop.cols,alpha=0.2),lwd=0.5)
  #Add axes
  axis(1,at=seq(-100,100,5),labels=NA,tck=-0.03)
  axis(2,at=seq(-100,100,5),labels=NA,tck=-0.03)
  sapply(seq(-100,100,5),function(i){
    axis(1,at=i,tick=F,labels=i,line=-0.8,cex.axis=0.7)
    axis(2,at=i,tick=F,labels=i,line=-0.6,cex.axis=0.7,las=2)
  })
  mtext(1,line=1.25,text="Principal Component 1",cex=1.05)
  mtext(2,line=1.5,text="Principal Component 2",cex=1.05)
}


################
###RSCRIPT BLOCK
################
require(optparse,quietly=T)
require(plot3D,quietly=T)
require(umap,quietly=T)
###List of command-line options
option_list <- list(
  make_option(c("--prefix"), type="character", default="gnomAD_v2_SV",
              help="prefix used for naming outfiles [default %default]",
              metavar="character"),
  make_option(c("-P", "--populations"), type="character", default=NULL,
              help="tab-delimited file specifying populations and HEX colors [default %default]",
              metavar="character")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog PCA pop_assignments PCRPLUS_samples batch_assignments OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 5){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
PCA.in <- args$args[1]
pop.assignments.in <- args$args[2]
PCRPLUS.samples.in <- args$args[3]
batch.assignments.in <- args$args[4]
OUTDIR <- args$args[5]
prefix <- opts$prefix
pops.file <- opts$populations

# #Dev parameters (local)
# PCA.in <- "~/scratch/gnomAD_v2_SV_MASTER.PCA_loadings.txt.gz"
# pop.assignments.in <- "~/scratch/gnomAD_v2_SV_MASTER.updated_pop_assignments.txt"
# PCRPLUS.samples.in <- "~/scratch/gnomAD_v2_SV_PCRPLUS.samples.list"
# batch.assignments.in <- "~/scratch/gnomAD_v2_SV.sample_batch_assignments.txt"
# OUTDIR <- "~/scratch/gnomAD_v2_test_plots/"
# prefix <- "gnomAD_v2_SV"
# pops.file <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/gnomAD_population_colors.txt"

###Create output directory, if necessary
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

###Sets populations & colors
if(!is.null(pops.file)){
  pops <- read.table(pops.file,sep="\t",header=T,comment.char="",check.names=F)
}

###Process input data
cat("NOW LOADING AND CLEANING DATA\n")
dat <- import.data(PCA.in,pop.assignments.in,PCRPLUS.samples.in,
                   batch.assignments.in,pops)

#Plot 3D scatterplot of PCs 1, 2, and 4
png(paste(OUTDIR,"/",prefix,".PCA_by_population.3d.png",sep=""),
    height=350*4,width=350*4,res=350)
scatter.3d(dat=dat,pops=pops)
dev.off()

#Plot 2D scatterplot of PCs 1 v 2
png(paste(OUTDIR,"/",prefix,".PCA_by_population.2d.png",sep=""),
    height=350*2.6,width=350*2.6,res=400)
scatter.2d(dat=dat,pops=pops)
dev.off()

#Plot UMAP of top 4 PCs
dat.umap <- umap(dat[,-c(1:4)][,1:4])
pop.cols <- as.character(sapply(dat$sample,function(ID){pops$color[which(pops$pop==as.character(dat$pop[which(dat$sample==ID)]))]}))
png(paste(OUTDIR,"/",prefix,".PCA_UMAP_by_population.png",sep=""),
    height=350*4,width=350*4,res=350)
par(mar=rep(0.25,4),bty="n")
plot(dat.umap$layout,col=pop.cols,pch=19,cex=0.25,
     xaxt="n",yaxt="n",xlab="",ylab="")
dev.off()
