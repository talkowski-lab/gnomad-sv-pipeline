#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Plot 3D PCA and supplementary PCA grids for formal gnomAD analysis


###Set master parameters
options(stringsAsFactors=F, scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Read input data and consolidate labels
import.data <- function(PCA.in, pop.assignments.in, PCRPLUS.samples.in, 
                        batch.assignments.in, famfile.in, pops, keep.n.pcs=20){
  #Read raw data
  pca <- read.table(PCA.in, header=T, comment.char="", sep="\t")[, 1:(keep.n.pcs+1)]
  colnames(pca)[1] <- "sample"
  pop.assignments <- read.table(pop.assignments.in, header=F, comment.char="", sep="\t")
  colnames(pop.assignments) <- c("sample", "pop")
  plus.samples.list <- read.table(PCRPLUS.samples.in, header=F)[, 1]
  pcr.assignments <- data.frame("sample"=pca$sample, "pcr"="PCRMINUS")
  pcr.assignments$pcr[which(pcr.assignments$sample %in% plus.samples.list)] <- "PCRPLUS"
  batch.assignments <- read.table(batch.assignments.in, header=F, comment.char="", sep="\t")
  colnames(batch.assignments) <- c("sample", "batch")
  famfile <- read.table(famfile.in, header=F, comment.char="", sep="\t")[, c(2, 5)]
  colnames(famfile) <- c("sample", "sexcode")
  famfile$sex <- NA
  famfile$sex[which(famfile$sexcode==1)] <- "MALE"
  famfile$sex[which(famfile$sexcode==2)] <- "FEMALE"
  famfile <- famfile[, c(1, 3)]
  #Merge data
  merged <- merge(x=pop.assignments, y=pcr.assignments, by="sample", sort=F, all=F)
  merged <- merge(x=merged, y=batch.assignments, by="sample", sort=F, all=F)
  merged <- merge(x=merged, y=famfile, by="sample", sort=F, all=F)
  merged <- merge(x=merged, y=pca, by="sample", sort=F, all=F)
  return(merged)
}


#####################
###PLOTTING FUNCTIONS
#####################
#3D scatterplot of PC1-3, colored by ancestry
scatter.3d <- function(dat, pops, theta=305, phi=25){
  #Get ancestry colors
  pop.cols <- as.character(sapply(dat$sample, function(ID){pops$color[which(pops$pop==as.character(dat$pop[which(dat$sample==ID)]))]}))
  # #Helper function to project small dots onto the background panels
  # panelfirst <- function(pmat){
  #   x <- dat$PC1
  #   y <- dat$PC2
  #   z <- dat$PC3
  #   XY <- trans3D(x, y, z = rep(min(z), length(z)), pmat = pmat)
  #   scatter2D(XY$x, XY$y, colvar = NULL, pch=19, 
  #             cex = 0.05, add = TRUE, col=adjustcolor(pop.cols, alpha=0.2))
  #   
  #   YZ <- trans3D(x = rep(max(x), length(x)), y, z, pmat = pmat)
  #   scatter2D(YZ$x, YZ$y, colvar = NULL, pch=19, 
  #             cex = 0.05, add = TRUE, col=adjustcolor(pop.cols, alpha=0.2))
  # }
  #Prep 3D plot
  par(mar=c(2, 0, 0, 0))
  scatter3D(x=dat$PC1, y=dat$PC2, z=dat$PC3, pch=19, cex=0.35, colvar=NULL, 
            col=adjustcolor(pop.cols, alpha=0.4), theta=theta, phi=phi, bty="b2", 
            xlab="PC 1", 
            ylab="PC 2", 
            zlab="PC 3", 
            ticktype="detailed")
}

#2D scatterplot, PCs 1v2
scatter.2d <- function(dat, pops, x.idx=1, y.idx=2, x.axis=T, y.axis=T, lab.cex=1){
  #Get ancestry colors
  pop.cols <- as.character(sapply(dat$sample, function(ID){pops$color[which(pops$pop==as.character(dat$pop[which(dat$sample==ID)]))]}))
  #Prep plot area
  par(mar=c(2.5, 2.5, 0.5, 0.5), bty="n")
  x.vals <- dat[, which(colnames(dat) == paste("PC", x.idx, sep=""))]
  y.vals <- dat[, which(colnames(dat) == paste("PC", y.idx, sep=""))]
  plot(x=range(x.vals, na.rm=T), 
       y=range(y.vals, na.rm=T), 
       type="n", xaxt="n", xlab="", yaxt="n", ylab="")
  #Add points
  points(x.vals, y.vals, pch=21, cex=0.2, col=pop.cols, 
         bg=adjustcolor(pop.cols, alpha=0.2), lwd=0.5)
  #Add axes
  if(x.axis==T){
    axis(1, at=seq(-100, 100, 10), labels=NA, tck=-0.03)
    sapply(seq(-100, 100, 10), function(i){
      axis(1, at=i, tick=F, labels=i, line=-0.8, cex.axis=0.7)
    })
    mtext(1, line=1.25, text=paste("Principal Component", x.idx), cex=1.05*lab.cex)
  }
  if(y.axis==T){
    axis(2, at=seq(-100, 100, 10), labels=NA, tck=-0.03)
    sapply(seq(-100, 100, 10), function(i){
      axis(2, at=i, tick=F, labels=i, line=-0.5, cex.axis=0.7, las=2)
    })
    mtext(2, line=1.5, text=paste("Principal Component", y.idx), cex=1.05*lab.cex)
  }
}




################
###RSCRIPT BLOCK
################
require(optparse, quietly=T)
require(plot3D, quietly=T)
require(umap, quietly=T)
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
args <- parse_args(OptionParser(usage="%prog PCA pop_assignments PCRPLUS_samples batch_assignments famfile OUTDIR", 
                                option_list=option_list), 
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 6){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
PCA.in <- args$args[1]
pop.assignments.in <- args$args[2]
PCRPLUS.samples.in <- args$args[3]
batch.assignments.in <- args$args[4]
famfile.in <- args$args[5]
OUTDIR <- args$args[6]
prefix <- opts$prefix
pops.file <- opts$populations

# #Dev parameters (local)
# PCA.in <- "~/scratch/gnomAD_v2_SV_MASTER.PCA_loadings.txt.gz"
# pop.assignments.in <- "~/scratch/gnomAD_v2_SV_MASTER.updated_pop_assignments.txt"
# PCRPLUS.samples.in <- "~/scratch/gnomAD_v2_SV_PCRPLUS.samples.list"
# batch.assignments.in <- "~/scratch/gnomAD_v2_SV.sample_batch_assignments.txt"
# famfile.in <- "~/scratch/merged_famfile.fam"
# OUTDIR <- "~/scratch/gnomAD_v2_test_plots/"
# prefix <- "gnomAD_v2_SV"
# pops.file <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/gnomAD_population_colors.txt"

###Create output directory, if necessary
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

###Sets populations & colors
if(!is.null(pops.file)){
  pops <- read.table(pops.file, sep="\t", header=T, comment.char="", check.names=F)
}

###Process input data
cat("NOW LOADING AND CLEANING DATA\n")
dat <- import.data(PCA.in, pop.assignments.in, PCRPLUS.samples.in, 
                   batch.assignments.in, famfile.in, pops)

###Set colors
pop.cols <- as.character(sapply(dat$sample, 
                                function(ID){pops$color[which(pops$pop==as.character(dat$pop[which(dat$sample==ID)]))]}))
sex.cols <- as.character(sapply(dat$sex, 
                                function(sex){
                                  if(is.na(sex)){
                                    pops$color[which(pops$pop=="OTH")]
                                  }else{
                                    if(sex=="MALE"){
                                      "#00B0CF"
                                    }else if(sex=="FEMALE"){
                                      "#F064A5"
                                    }else{pops$color[which(pops$pop=="OTH")]
                                    }
                                  }
                                }))
pcr.cols <- as.character(sapply(dat$pcr, 
                                function(pcr){
                                  if(pcr=="PCRPLUS"){
                                    "#1F6AAB"
                                  }else if(pcr=="PCRMINUS"){
                                    "#E01925"
                                  }else{pops$color[which(pops$pop=="OTH")]
                                  }
                                }))


#Plot 3D scatterplot of PCs 1-3
png(paste(OUTDIR, "/", prefix, ".PCA_by_population.3d.png", sep=""), 
    height=350*4, width=350*4, res=350)
scatter.3d(dat=dat, pops=pops, theta=220, phi=30)
dev.off()

#Plot 2D scatterplot of PCs 1 v 2
png(paste(OUTDIR, "/", prefix, ".PCA_by_population.2d.png", sep=""), 
    height=350*2.6, width=350*2.6, res=400)
scatter.2d(dat=dat, pops=pops)
dev.off()

#Plot UMAP of top 3 PCs
dat.umap <- umap(dat[, -c(1:5)][, 1:3])
png(paste(OUTDIR, "/", prefix, ".PCA_UMAP_by_population.png", sep=""), 
    height=350*4, width=350*4, res=350)
par(mar=rep(0.25, 4), bty="n")
plot(dat.umap$layout, col=pop.cols, pch=19, cex=0.25, 
     xaxt="n", yaxt="n", xlab="", ylab="")
dev.off()

# #Plot top 10 PCs, colored by population, batch, sex, and PCR status
# #By pop
# png(paste(OUTDIR, "/", prefix, ".PCA_grid.col_by_pop.png", sep=""), 
#     height=350*8, width=350*8, res=350)
# par(oma=rep(1,4))
# pairs(dat[, 6:15], col=pop.cols, cex=0.15, xaxt="n", yaxt="n")
# dev.off()
# #By sex
# png(paste(OUTDIR, "/", prefix, ".PCA_grid.col_by_sex.png", sep=""), 
#     height=350*8, width=350*8, res=350)
# par(oma=rep(1,4))
# pairs(dat[, 6:15], col=sex.cols, cex=0.15, xaxt="n", yaxt="n")
# dev.off()
# #By PCR
# png(paste(OUTDIR, "/", prefix, ".PCA_grid.col_by_pcr.png", sep=""), 
#     height=350*8, width=350*8, res=350)
# par(oma=rep(1,4))
# pairs(dat[, 6:15], col=pcr.cols, cex=0.15, xaxt="n", yaxt="n")
# dev.off()

#Plot grid of top 3 PCs, colored by population, for supplement
# png(paste(OUTDIR, "/", prefix, ".topThreePCs_by_pop.png", sep=""),
#     height=350*6, width=350*6, res=350)
# par(oma=rep(1,4))
# pairs(dat[, 6:8], col=pop.cols, cex=0.15)
# dev.off()

#Plot strip of top 3 PCs, colored by population, for supplement
png(paste(OUTDIR, "/", prefix, ".topThreePCs_by_pop.strip.png", sep=""),
    height=350*2.6, width=3*350*2.6, res=400)
par(mfrow=c(1, 3))
scatter.2d(dat=dat, pops=pops, lab.cex=0.7)
mtext(3, text="PC1 x PC2", font=2, line=-0.75, cex=0.8)
scatter.2d(dat=dat, pops=pops, x.idx=1, y.idx=3, lab.cex=0.7)
mtext(3, text="PC1 x PC3", font=2, line=-0.75, cex=0.8)
scatter.2d(dat=dat, pops=pops, x.idx=2, y.idx=3, lab.cex=0.7)
mtext(3, text="PC2 x PC3", font=2, line=-0.75, cex=0.8)
dev.off()




