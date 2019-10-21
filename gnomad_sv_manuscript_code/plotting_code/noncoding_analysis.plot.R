#!/usr/bin/env Rscript

# Copyright (c) 2019 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Plot APS-based noncoding CWAS results


### Set master parameters
options(stringsAsFactors=F, scipen=1000)


####################
### HELPER FUNCTIONS
####################
# Read & clean vcf2bed input
read.vcf2bed <- function(vcf2bed.in, all.aps.in=NULL, sdsr.cov.in=NULL, 
                         phastCons.in=NULL, exon.hits.in=NULL, exclude.allosomes=T){
  #Read data
  dat <- read.table(vcf2bed.in, comment.char="", header=T, sep="\t")
  colnames(dat)[1] <- "chrom"
  if(exclude.allosomes==T){
    dat <- dat[which(dat$chrom %in% 1:22), ]
  }
  #Restrict to sites with at least one observed alternative allele
  dat <- dat[union(which(dat$AC>0), grep("MULTIALLELIC", dat$FILTER, fixed=T)), ]
  #Drop columns not being used (to save memory)
  cols.to.drop <- c("start", "end", "CHR2", "CPX_INTERVALS", 
                    "END", "SOURCE", "STRANDS", "UNRESOLVED_TYPE", 
                    "LINCRNA__LOF", "LINCRNA__DUP_LOF", "LINCRNA__COPY_GAIN", 
                    "LINCRNA__DUP_PARTIAL", "LINCRNA__MSV_EXON_OVR", 
                    "LINCRNA__INTRONIC", "LINCRNA__INV_SPAN", "LINCRNA__UTR")
  dat <- dat[, -which(colnames(dat) %in% cols.to.drop)]
  #Convert numeric columns
  numeric.columns <- sort(unique(c(grep("FREQ", colnames(dat), fixed=T), 
                                   grep("AN", colnames(dat), fixed=T), 
                                   grep("AC", colnames(dat), fixed=T), 
                                   grep("AF", colnames(dat), fixed=T))))
  numeric.columns <- setdiff(numeric.columns, grep("SPAN", colnames(dat), fixed=T))
  dat[, numeric.columns] <- apply(dat[, numeric.columns], 2, as.numeric)
  #Read & add APS, if optioned
  if(!is.null(all.aps.in)){
    aps <- read.table(all.aps.in, header=T, sep="\t", comment.char="")
    colnames(aps)[1] <- "VID"
    dat <- merge(dat, aps, by.x="name", by.y="VID", all.x=T, all.y=F, sort=F)
  }
  # Read & add SD/SR coverage, if optioned
  if(!is.null(sdsr.cov.in)){
    sdsr.cov <- read.table(sdsr.cov.in, header=F, sep="\t")
    colnames(sdsr.cov) <- c("name", "SDSR_COV")
    dat <- merge(dat, sdsr.cov, by="name", all.x=T, all.y=F, sort=F)
  }
  # Read & add phastCons data, if optioned
  if(!is.null(phastCons.in)){
    phastCons <- read.table(phastCons.in, header=T, sep="\t", comment.char="")
    phastCons[, -1] <- apply(phastCons[, -1], 2, as.numeric)
    colnames(phastCons)[1] <- "name"
    dat <- merge(dat, phastCons, by="name", all.x=T, all.y=F, sort=F)
  }
  # Read & add # of exons overlapped, if optioned
  if(!is.null(exon.hits.in)){
    exons <- read.table(exon.hits.in, header=F, sep="\t")
    colnames(exons) <- c("name", "EXONS_OVERLAPPED")
    dat <- merge(dat, exons, by="name", all.x=T, all.y=F, sort=F)
  }
  return(dat)
}

# Load APS test stats
import.stats <- function(stats.in, min.sv=10){
  dat <- read.table(stats.in, header=T, comment.char="", sep="\t")
  colnames(dat)[1] <- "svi"
  dat <- dat[which(dat$N_SV >= min.sv), ]
  return(dat)
}

# Process SNV gene data (for constrained gene info)
read.snvdata <- function(SNVdata.in, gene.metadata.in){
  # Read & clean
  snv.data <- read.table(SNVdata.in, header=T, comment.char="")
  metadata <- read.table(gene.metadata.in, header=T)
  merged <- merge(x=snv.data, y=metadata, by="gene", sort=F)
  merged <- merged[which(!(merged$chrom %in% c("chrX", "chrY"))), ]
  # Assign oe deciles
  merged$mis_oe_dec <- ceiling(10*rank(merged$mis_oe)/(nrow(merged)+1))
  merged$ptv_oe_dec <- ceiling(10*rank(merged$ptv_oe)/(nrow(merged)+1))
  # Assign oe to 40 bins
  merged$mis_oe_binrank <- ceiling(40*rank(merged$mis_oe)/(nrow(merged)+1))
  merged$ptv_oe_binrank <- ceiling(40*rank(merged$ptv_oe)/(nrow(merged)+1))
  # Assign oe percentiles
  merged$mis_oe_cent <- ceiling(100*rank(merged$mis_oe)/(nrow(merged)+1))
  merged$ptv_oe_cent <- ceiling(100*rank(merged$ptv_oe)/(nrow(merged)+1))
  # Return formatted data
  return(merged)
}


# Calculate APS mean & 95% CI with bootstrapping
calc.aps <- function(v.aps, boot.n=100, conf=0.95){
  v.aps <- v.aps[which(!is.na(v.aps))]
  helper.getFracSingle <- function(v.aps, indices){sum(v.aps[indices])/length(v.aps[indices])}
  point.est <- helper.getFracSingle(v.aps, indices=1:length(v.aps))
  calc.ci <- function(v.aps, n, conf){
    set.seed(0)
    boot.obj <- boot(data=v.aps, statistic=helper.getFracSingle, R=n)
    ci <- boot.ci(boot.obj, conf=conf, type="basic")$basic[4:5]
    return(ci)
  }
  ci <- calc.ci(v.aps, n=boot.n, conf=conf)
  return(c(point.est, ci))
}

# Gather values for distance-based analyses
get.dist.vals <- function(dat){
  annos <- sort(unique(dat$anno))
  annos <- annos[which(!(annos %in% c("ANY", "NONE")))]
  dists <- c("IN_GENE", "10KB", "100KB", "1MB")
  res <- lapply(c("DEL", "DUP"), function(CNV){
    as.data.frame(sapply(annos, function(anno){
      as.numeric(sapply(dists, function(dist){
        dat$APS[which(dat$svi==paste("PARTIAL", CNV, sep="_")
                      & dat$anno==anno
                      & dat$cons=="ANY"
                      & dat$gset=="ANY"
                      & dat$grel==dist)]
      }))
    }))
  })
  names(res) <- c("DEL", "DUP")
  return(res)
}

# Gather values for subpanel dotplot for a single annotation category
get.blowout.vals <- function(dat, anno){
  vals <- list(cbind(rbind(dat[which(dat$svi=="PARTIAL_DEL" 
                                     & dat$anno==anno & dat$cons=="ANY" 
                                     & dat$gset=="ANY" & dat$grel=="ANY"), 7:10],
                           dat[which(dat$svi=="PARTIAL_DUP" 
                                     & dat$anno==anno & dat$cons=="ANY" 
                                     & dat$gset=="ANY" & dat$grel=="ANY"), 7:10]),
                     "color"=c(del.color, dup.color),
                     "alpha"=c(1, 1)),
               cbind(rbind(dat[which(dat$svi=="FULL_DEL" 
                                     & dat$anno==anno & dat$cons=="ANY" 
                                     & dat$gset=="ANY" & dat$grel=="ANY"), 7:10],
                           dat[which(dat$svi=="FULL_DUP" 
                                     & dat$anno==anno & dat$cons=="ANY" 
                                     & dat$gset=="ANY" & dat$grel=="ANY"), 7:10]),
                     "color"=c(del.color, dup.color),
                     "alpha"=c(1, 1)),
               cbind(rbind(dat[which(dat$svi=="PARTIAL_DEL" 
                                     & dat$anno==anno & dat$cons=="UNCONSERVED" 
                                     & dat$gset=="ANY" & dat$grel=="ANY"), 7:10],
                           dat[which(dat$svi=="PARTIAL_DUP" 
                                     & dat$anno==anno & dat$cons=="UNCONSERVED" 
                                     & dat$gset=="ANY" & dat$grel=="ANY"), 7:10]),
                     "color"=c(del.color, dup.color),
                     "alpha"=c(1, 1)),
               cbind(rbind(dat[which(dat$svi=="PARTIAL_DEL" 
                                     & dat$anno==anno & dat$cons=="CONSERVED" 
                                     & dat$gset=="ANY" & dat$grel=="ANY"), 7:10],
                           dat[which(dat$svi=="PARTIAL_DUP" 
                                     & dat$anno==anno & dat$cons=="CONSERVED" 
                                     & dat$gset=="ANY" & dat$grel=="ANY"), 7:10]),
                     "color"=c(del.color, dup.color),
                     "alpha"=c(1, 1)),
               cbind(rbind(dat[which(dat$svi=="PARTIAL_DEL" 
                                     & dat$anno==anno & dat$cons=="ANY" 
                                     & dat$gset=="unconstrained" & dat$grel=="1MB"), 7:10],
                           dat[which(dat$svi=="PARTIAL_DUP" 
                                     & dat$anno==anno & dat$cons=="ANY" 
                                     & dat$gset=="unconstrained" & dat$grel=="1MB"), 7:10]),
                     "color"=c(del.color, dup.color),
                     "alpha"=c(1, 1)),
               cbind(rbind(dat[which(dat$svi=="PARTIAL_DEL" 
                                     & dat$anno==anno & dat$cons=="ANY" 
                                     & dat$gset=="constrained" & dat$grel=="1MB"), 7:10],
                           dat[which(dat$svi=="PARTIAL_DUP" 
                                     & dat$anno==anno & dat$cons=="ANY" 
                                     & dat$gset=="constrained" & dat$grel=="1MB"), 7:10]),
                     "color"=c(del.color, dup.color),
                     "alpha"=c(1, 1)))
  names(vals) <- c("Partially overlaps element",
                   "Overlaps entire element",
                   "Unconserved elements",
                   "Conserved elements",
                   "Near (\u00B11Mb) unconstrained genes",
                   "Near (\u00B11Mb) constrained genes")
  return(vals)
}

# Compute APS over various phastCons intervals
aps.by.phastcons <- function(bed, measure="sum", n.bins=100){
  # Get values
  measure.idx <- which(colnames(bed)==paste("phastCons", measure, sep="_"))
  vals <- as.numeric(bed[, measure.idx])
  vals <- vals[which(!is.na(vals))]
  bins <- quantile(vals, probs=seq(0, 1, by=1/n.bins))
  res <- t(sapply(1:(length(bins)-1), function(i){
    min.val <- bins[i]
    max.val <- bins[i+1]
    calc.aps(bed$APS[which(bed[, measure.idx]>=min.val
                           & bed[, measure.idx]<=max.val)])
  }))
  colnames(res) <- c("APS", "lower_CI", "upper_CI")
  return(res)
}


#################
### PLOTTING CODE
#################
# Helper function to plot a single point estimate & CI
add.point <- function(at, val, ci, horiz=T, color="black", alpha=1, ci.lwd=1, ci.alpha=1){
  if(horiz==T){
    segments(x0=min(ci), x1=max(ci), y0=at, y1=at, col="white", lend="round", lwd=ci.lwd)
    segments(x0=min(ci), x1=max(ci), y0=at, y1=at, 
             col=adjustcolor(color, ci.alpha), lend="round", lwd=2)
    points(x=val, y=at, pch=19, col="white")
    points(x=val, y=at, pch=21, col=color, bg=adjustcolor(color, alpha))
  }else{
    segments(y0=min(ci), y1=max(ci), x0=at, x1=at, col="white", lend="round", lwd=ci.lwd)
    segments(y0=min(ci), y1=max(ci), x0=at, x1=at, 
             col=adjustcolor(color, ci.alpha), lend="round", lwd=ci.lwd)
    points(y=val, x=at, pch=19, col="white")
    points(y=val, x=at, pch=21, col=color, bg=adjustcolor(color, alpha))
  }
}

# Helper function to plot a series of point estimates & CIs
plot.dots <- function(values, sig.cutoff=0.05, parmar=c(0.5, 8, 3, 0.5),
                      bg.lines=list(c(0, "gray70", 2)), 
                      ci.lwd=1.5, horiz=F){
  # Values should be a list of data frames, where each data frame has columns:
  # 1) APS (point estimate)
  # 2) lower_CI
  # 3) upper_CI
  # 4) p
  # 5) color
  # 6) alpha
  # Names of values should be desired category labels
  
  # Get plot values
  n.rows <- length(values)
  aps.range <- as.numeric(range(unlist(lapply(values, function(df){if(!is.na(df)){range(df[, 1], na.rm=T)}}))))
  
  # Prep plot area
  par(bty="n", mar=parmar)
  y.at <- seq(-1, 1, 0.1)
  if(horiz==T){
    plot(y=c(min(aps.range)-0.05, max(aps.range)+0.05), 
         x=c(0.5, n.rows-0.5), type="n",
         xaxt="n", xlab="", yaxt="n", ylab="")
    if(!is.null(bg.lines)){
      lapply(bg.lines, function(v){
        abline(h=v[1], col=v[2], lty=as.numeric(v[3]))
      })
    }
    axis(2, at=y.at, labels=NA, tck=-0.03)
    axis(2, at=y.at, tick=F, line=-0.7, cex.axis=0.85, las=2)
    mtext(2, line=1.6, text="APS")
  }else{
    plot(x=c(min(aps.range)-0.05, max(aps.range)+0.05), 
         y=c(-0.5, -n.rows+0.5), type="n",
         xaxt="n", xlab="", yaxt="n", ylab="")
    if(!is.null(bg.lines)){
      lapply(bg.lines, function(v){
        abline(v=v[1], col=v[2], lty=as.numeric(v[3]))
      })
    }
    axis(3, at=y.at, labels=NA, tck=-0.03)
    axis(3, at=y.at, tick=F, line=-0.7, cex.axis=0.85)
    mtext(3, line=1.1, text="APS")
  }
  
  # Add points & labels
  sapply(1:n.rows, function(i){
    df <- values[[i]]
    if(!is.na(df)){
      n.subrows <- nrow(df)
      if(n.subrows > 0){
        buffer <- 1/(n.subrows+3)
        sapply(n.subrows:1, function(k){
          if(horiz==T){
            add.point(at=(i-1)+((k+1)*buffer), 
                      val=as.numeric(df[k, 1]),
                      ci=as.numeric(df[k, 2:3]),
                      color=df[k, 5], alpha=df[k, 6], 
                      ci.lwd=ci.lwd, horiz=F)
          }else{
            add.point(at=(-i+1)-((k+1)*buffer), 
                      val=as.numeric(df[k, 1]),
                      ci=as.numeric(df[k, 2:3]),
                      color=df[k, 5], alpha=df[k, 6], 
                      ci.lwd=ci.lwd, horiz=T)
          }
        })
      }
      if(horiz==F){
        axis(2, at=-i+c(0.15, 0.85), tck=0, col="gray85", labels=NA)
        axis(2, at=-i+0.5, tick=F, las=2, cex.axis=0.85, line=-0.9,
             labels=names(values)[i])
      }else{
        axis(1, at=i-c(0.15, 0.5), tck=0, col="gray85", labels=NA)
        text(x=i, y=par("usr")[3]-(0.05*(par("usr")[4]-par("usr")[3])), 
             srt=45, cex=0.85, xpd=T, pos=2, labels=names(values)[i])
      }
    }
  })
}

### Plot swarms of partial- vs whole-element CNVs
plot.paired_cnvs <- function(dat, mode="svi"){
  # Get plot data
  annos <- sort(unique(dat$anno))
  annos <- annos[which(!(annos %in% c("ANY", "NONE")))]
  plot.dat <- do.call("cbind", lapply(c("DEL", "DUP"), function(CNV){
    if(mode=="svi"){
      do.call("cbind", lapply(c("PARTIAL", "FULL"), function(ovr){
        sapply(annos, function(anno){
          dat$APS[which(dat$anno==anno & dat$svi==paste(ovr, CNV, sep="_")
                        & dat$cons=="ANY" & dat$gset=="ANY" & dat$grel=="ANY")]
        })
      }))
    }else if(mode=="cons"){
      do.call("cbind", lapply(c("UNCONSERVED", "CONSERVED"), function(cons){
        sapply(annos, function(anno){
          d <- dat$APS[which(dat$anno==anno & dat$svi==paste("PARTIAL", CNV, sep="_")
                             & dat$cons==cons & dat$gset=="ANY" & dat$grel=="ANY")]
          if(length(d)==0){
            d <- NA
          }
          return(d)
        })
      }))
    }
  }))
  colnames(plot.dat) <- c("DEL.partial", "DEL.full",
                          "DUP.partial", "DUP.full")
  ylims <- range(plot.dat, na.rm=T)
  colors <- c(rep(del.color, 2), rep(dup.color, 2))
  
  # Prep plot area
  par(bty="n", mar=c(2, 2.5, 1.25, 0.25))
  plot(x=c(0, 4), y=c(min(ylims), max(ylims)+0.03), type="n",
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i")
  abline(h=0, lty=2, col="gray50")
  sapply(1:4, function(x){
    axis(1, at=x-c(0.8, 0.2), tck=0, labels=NA)
    if(mode=="svi"){
      axis(1, at=x-0.5, labels=rep(c("Partial", "Full"), 2)[x],
           line=-1, cex.axis=0.9, tick=F)
    }else if(mode=="cons"){
      axis(1, at=x-0.5, labels=rep(c("Uncons.", "Cons."), 2)[x],
           line=-1, cex.axis=0.9, tick=F)
    }
  })
  
  axis(1, at=c(0.2, 1.8), tck=0, labels=NA, line=1)
  axis(1, at=c(2.2, 3.8), tck=0, labels=NA, line=1)
  axis(1, at=c(1, 3), tick=F, line=0.1, labels=c("DEL", "DUP"))
  axis(2, at=axTicks(2), tck=-0.025, labels=NA)
  axis(2, at=axTicks(2), tick=F, las=2, cex.axis=0.85, line=-0.7)
  mtext(2, text="APS", line=1.5)
  
  # Add violins
  sapply(1:4, function(i){
    mean.i <- mean(plot.dat[, i], na.rm=T)
    beeswarm(plot.dat[, i], add=T, at=i-0.5, pch=19, col=colors[i], cex=0.8)
    segments(x0=i-0.3, x1=i-0.7, y0=mean.i, y1=mean.i, lwd=2, lend="round")
    # text(x=i-0.5, y=mean.i, pos=3, cex=0.8, labels=formatC(round(mean.i, 2), digits=2), font=2)
  })
  
  # Add p-values
  lapply(list(c(1, 2), c(3, 4)), function(k){
    p <- format(t.test(plot.dat[, k[1]], plot.dat[, k[2]],
                       alternative="less", paired=T, na.rm=T)$p.value, scientific=T)
    p.parts <- as.numeric(unlist(strsplit(p, split="e")))
    axis(3, at=k-0.5, tck=0.03, labels=NA)
    axis(3, at=min(k), line=-0.9, tick=F, cex.axis=0.9,
         labels=bquote(italic("P") == .(formatC(round(p.parts[1], 1), digits=2)) * "x10"^{.(formatC(round(p.parts[2], 0), digits=0))}))
  })
}

# Plot phastCons correlation with APS
plot.phastCons.dots <- function(bed, measure="sum", n.bins=100, 
                                color, title=NULL, ylims=NULL, 
                                xlabel="phastCons Sum Percentile", 
                                ax.labels=TRUE, cex.labels=0.8, cex.ax.titles=1, 
                                tck=NULL, yline=-0.4, parmar=c(2, 3.25, 1.75, 2)){
  require(zoo, quietly=T)
  # Compute binned APS
  means <- aps.by.phastcons(bed, measure=measure, n.bins=n.bins)[, 1]
  d <- (1:length(means))
  
  #Prep plot area
  if(is.null(ylims)){
    ylims <- range(means, na.rm=T)
  }
  par(mar=parmar)
  plot(x=c(0, length(means)), y=ylims, type="n", 
       xaxt="n", yaxt="n", xlab="", ylab="")
  # axis(1, at=seq(0, 100, 10), labels=NA)
  sapply(seq(0, 100, 25), function(p){
    axis(1, at=p, tick=F, cex.axis=cex.labels, line=-0.9, 
         labels=bquote(.(p)^'th'))
  })
  axis(2, at=round(axTicks(2), 2), labels=NA, tck=tck)
  if(ax.labels==T){
    mtext(1, line=1, text=xlabel, cex=cex.ax.titles)
    mtext(2, line=2.15, text="APS", cex=cex.ax.titles)
  }
  axis(2, at=round(axTicks(2), 2), line=yline, cex.axis=cex.labels, tick=F, las=2)
  mtext(3, text=title, line=0.2, cex=0.9)
  abline(h=0, lty=2, col="gray70")
  
  #Add points & rolling mean
  points(x=d-0.5, y=means, pch=21, col=color, cex=0.4)
  points(x=d-0.5, y=rollapply(means, 21, mean, na.rm=T, partial=T), 
         lwd=2, type="l", col=color)
  
  #Add correlation coefficient
  cor.res <- cor.test(x=d-0.5, y=means, method="spearman")
  text(x=par("usr")[1]-(0.025*(par("usr")[2]-par("usr")[1])), 
       y=par("usr")[3]+(0.885*(par("usr")[4]-par("usr")[3])), 
       pos=4, cex=cex.labels, 
       labels=bquote(rho == .(format(round(cor.res$estimate, 2), nsmall=2))))
  # cor.res <- cor.test(x=d-0.5, y=means, method="pearson")
  cor.p <- cor.res$p.value
  if(cor.p>10^-100){
    cor.p <- format(cor.p, scientific=T)
    cor.p.base <- format(as.numeric(strsplit(cor.p, split="e")[[1]][1]), nsmall=2)
    cor.p.exp <- format(as.numeric(strsplit(cor.p, split="e")[[1]][2]), nsmall=0)
    text(x=par("usr")[2]+(0.025*(par("usr")[2]-par("usr")[1])), 
         y=par("usr")[3]+(0.075*(par("usr")[4]-par("usr")[3])), 
         pos=2, cex=cex.labels, 
         labels=bquote(italic(P) == .(format(round(as.numeric(cor.p.base), 2), nsmall=2))*"x"*10^.(cor.p.exp)))
  }else{
    text(x=par("usr")[2]+(0.025*(par("usr")[2]-par("usr")[1])), 
         y=par("usr")[3]+(0.075*(par("usr")[4]-par("usr")[3])), 
         pos=2, cex=cex.labels, 
         labels=bquote(italic(P) < 10^-100))
  }
  box()
}


#################
### RSCRIPT BLOCK
#################
require(optparse, quietly=T)
require(boot, quietly=T)
require(beeswarm, quietly=T)
### List of command-line options
option_list <- list(
  make_option(c("--prefix"), type="character", default="gnomAD_v2_SV", 
              help="prefix used for naming outfiles [default %default]", 
              metavar="character"), 
  make_option(c("-S", "--svtypes"), type="character", default=NULL, 
              help="tab-delimited file specifying SV types and HEX colors [default %default]", 
              metavar="character")
)

### Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog APS_STATS OUTDIR", 
                                option_list=option_list), 
                   positional_arguments=TRUE)
opts <- args$options

### Checks for appropriate positional arguments
if(length(args$args) != 9){
  stop("Incorrect number of required positional arguments\n")
}

### Writes args & opts to vars
vcf2bed.in <- args$args[1]
all.aps.in <- args$args[2]
sdsr.cov.in <- args$args[3]
phastCons.in <- args$args[4]
exon.hits.in <- args$args[5]
SNVdata.in <- args$args[6]
gene_metadata.in <- args$args[7]
stats.in <- args$args[8]
OUTDIR <- args$args[9]
prefix <- opts$prefix
svtypes.file <- opts$svtypes

# # Dev parameters (local)
# vcf2bed.in <- "~/scratch/gnomAD-SV_v2_rev1.vcf2bed.bed.gz"
# all.aps.in <- "~/scratch/gnomAD-SV_v2_rev1.APS_stats.txt"
# sdsr.cov.in <- "~/scratch/gnomAD-SV_v2_rev1.variant_SD_SR_coverage.txt.gz"
# phastCons.in <- "~/scratch/gnomAD-SV_v2_rev1.phastCons_stats.gz"
# exon.hits.in <- "~/scratch/gnomAD-SV_v2_rev1.exons_overlapped_per_SV.txt.gz"
# SNVdata.in <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/gnomAD_v2.1_canonical_constraint.condensed.txt.gz"
# gene_metadata.in <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/Gencode.v19.autosomal_canonical_gene_metadata.txt.gz"
# stats.in <- "~/scratch/gnomAD-SV_v2_rev1.merged_cwas_aps_stats.txt"
# OUTDIR <- "~/scratch/gnomAD_v2_test_plots/"
# prefix <- "gnomAD-SV_v2_rev1"
# svtypes.file <- "~/Desktop/Collins/Talkowski/code/sv-pipeline/ref/vcf_qc_refs/SV_colors.txt"

### Create output directory, if necessary
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

### Sets sv types & colors
if(!is.null(svtypes.file)){
  svtypes <- read.table(svtypes.file, sep="\t", header=F, comment.char="", check.names=F)
  svtypes <- as.data.frame(apply(svtypes, 2, as.character))
  colnames(svtypes) <- c("svtype", "color")
}else{
  require(RColorBrewer, quietly=T)
  svtypes.v <- unique(dat$SVTYPE)
  svtypes.c <- brewer.pal(length(svtypes.v), "Dark2")
  svtypes <- data.frame("svtype"=svtypes.v, 
                        "color"=svtypes.c)
}
del.color <- svtypes$color[which(svtypes$svtype=="DEL")]
dup.color <- svtypes$color[which(svtypes$svtype=="DUP")]
mcnv.color <- svtypes$color[which(svtypes$svtype=="MCNV")]

### Set other parameters
anno.labels <- c("ANY"="Any annotation",
                 "CHMM_S12"="ChromHMM bivalent enhancers",
                 "CHMM_S13"="ChromHMM polycomb repressed",
                 "CHMM_S6"="ChromHMM genic enhancers",
                 "CHMM_S7"="ChromHMM enhancers",
                 "DHS"="DNAseI hypersensitive sites",
                 "ENH"="Enhancer Atlas predictions",
                 "HAR"="Human accelerated regions",
                 "LOOP"="Chromatin loop boundaries",
                 "NONE"="No annotations",
                 "RECOMB"="Recombination hotspots",
                 "SUPER"="Predicted super enhancers",
                 "TBR"="TAD boundaries",
                 "TFBS"="TF binding sites",
                 "UCNE"="Ultraconserved elements",
                 "VISTA"="VISTA validated enhancers")

### Process input data
bed.all <- read.vcf2bed(vcf2bed.in, all.aps.in, sdsr.cov.in, phastCons.in, exon.hits.in)
bed <- bed.all[which(bed.all$FILTER=="PASS" & bed.all$SDSR_COV < 0.3), ]
dat <- import.stats(stats.in, min.sv=10)
gene.data <- read.snvdata(SNVdata.in, gene_metadata.in)
constrained.genes <- gene.data$gene[which(gene.data$pli>0.9)]

### Compute baseline APS from full vcf2bed
noncoding.idxs <- which(bed$EXONS_OVERLAPPED==0)
intergenic.idxs <- which(is.na(bed$PROTEIN_CODING__LOF)
                         & is.na(bed$PROTEIN_CODING__COPY_GAIN)
                         & is.na(bed$PROTEIN_CODING__DUP_LOF)
                         & is.na(bed$PROTEIN_CODING__DUP_PARTIAL)
                         & is.na(bed$PROTEIN_CODING__INTRONIC)
                         & is.na(bed$PROTEIN_CODING__PROMOTER)
                         & is.na(bed$PROTEIN_CODING__UTR)
                         & bed$PROTEIN_CODING__INTERGENIC=="True")
constrained.lof.idxs <- which(sapply(bed$PROTEIN_CODING__LOF, function(genes){
  if(is.na(genes)){
    return(F)
  }else{
    if(any(sapply(strsplit(genes, split=",", fixed=T), function(g){g %in% constrained.genes}))){
      return(T)
    }else{
      return(F)
    }
  }
}))
constrained.ied.idxs <- which(sapply(bed$PROTEIN_CODING__DUP_LOF, function(genes){
  if(is.na(genes)){
    return(F)
  }else{
    if(any(sapply(strsplit(genes, split=",", fixed=T), function(g){g %in% constrained.genes}))){
      return(T)
    }else{
      return(F)
    }
  }
}))
constrained.lof_ied.idxs <- unique(c(constrained.lof.idxs, constrained.ied.idxs))
baseline.aps <- t(data.frame("all"=calc.aps(bed$APS[which(bed$SDSR_COV<0.3)]),
                             "all_DEL"=calc.aps(bed$APS[which(bed$SDSR_COV<0.3 & bed$SVTYPE=="DEL")]),
                             "all_DUP"=calc.aps(bed$APS[which(bed$SDSR_COV<0.3 & bed$SVTYPE=="DUP")]),
                             "intergenic"=calc.aps(bed$APS[intersect(intergenic.idxs, which(bed$SDSR_COV<0.3))]),
                             "intergenic_DEL"=calc.aps(bed$APS[intersect(intergenic.idxs, which(bed$SDSR_COV<0.3 & bed$SVTYPE=="DEL"))]),
                             "intergenic_DUP"=calc.aps(bed$APS[intersect(intergenic.idxs, which(bed$SDSR_COV<0.3 & bed$SVTYPE=="DUP"))]),
                             "coding_DEL"=calc.aps(bed$APS[which(bed$SDSR_COV<0.3 & bed$SVTYPE=="DEL" & (!is.na(bed$PROTEIN_CODING__LOF) | !is.na(bed$PROTEIN_CODING__DUP_LOF)))]),
                             "coding_DUP"=calc.aps(bed$APS[which(bed$SDSR_COV<0.3 & bed$SVTYPE=="DUP" & (!is.na(bed$PROTEIN_CODING__LOF) | !is.na(bed$PROTEIN_CODING__DUP_LOF)))]),
                             "constrained_DEL"=calc.aps(bed$APS[intersect(constrained.lof_ied.idxs, which(bed$SDSR_COV<0.3 & bed$SVTYPE=="DEL"))]),
                             "constrained_DUP"=calc.aps(bed$APS[intersect(constrained.lof_ied.idxs, which(bed$SDSR_COV<0.3 & bed$SVTYPE=="DUP"))])))
colnames(baseline.aps) <- c("APS", "lower_CI", "upper_CI")

### Plot master panel
# Gather data
main.values.bench <- list("All intergenic CNVs"=cbind(baseline.aps[which(rownames(baseline.aps) %in% c("intergenic_DEL", "intergenic_DUP")), ],
                                                      "p"=c(NA, NA),
                                                      "color"=c(del.color, dup.color),
                                                      "alpha"=c(0.15, 0.15)),
                          "Protein-altering (pLoF & IED) CNVs"=cbind(baseline.aps[which(rownames(baseline.aps) %in% c("coding_DEL", "coding_DUP")), ],
                                                                     "p"=c(NA, NA),
                                                                     "color"=c(del.color, dup.color),
                                                                     "alpha"=c(1, 1)),
                          "Protein-altering (constrained genes)"=cbind(baseline.aps[which(rownames(baseline.aps) %in% c("constrained_DEL", "constrained_DUP")), ],
                                                                       "p"=c(NA, NA),
                                                                       "color"=c(del.color, dup.color),
                                                                       "alpha"=c(1, 1)))
annos <- unique(dat$anno)
annos.max <- sapply(annos, function(anno){
  max(dat$APS[intersect(grep("PARTIAL", dat$svi, fixed=T),
                        which(dat$anno==anno & dat$cons=="ANY" & dat$gset=="ANY" & dat$grel=="ANY"))])
})
annos.alpha <- sapply(annos, function(anno){
  ps <- dat$p[intersect(grep("PARTIAL", dat$svi, fixed=T),
                        which(dat$anno==anno & dat$cons=="ANY" & dat$gset=="ANY" & dat$grel=="ANY"))]
  sapply(ps, function(p){if(p<0.05/(2*length(annos))){1}else{0.15}})
})
main.values.annos <- lapply(annos[order(-annos.max)], function(anno){
  cbind(rbind(dat[which(dat$svi=="PARTIAL_DEL" 
                        & dat$anno==anno & dat$cons=="ANY" 
                        & dat$gset=="ANY" & dat$grel=="ANY"), 7:10],
              dat[which(dat$svi=="PARTIAL_DUP" 
                        & dat$anno==anno & dat$cons=="ANY" 
                        & dat$gset=="ANY" & dat$grel=="ANY"), 7:10]),
        "color"=c(del.color, dup.color),
        "alpha"=annos.alpha[, which(colnames(annos.alpha)==anno)])
})
names(main.values.annos) <- sapply(annos[order(-annos.max)], function(anno){anno.labels[which(names(anno.labels)==anno)]})
main.values <- c(main.values.bench,
                 NA, 
                 main.values.annos[which(names(main.values.annos) == "No annotations")],
                 main.values.annos[which(names(main.values.annos) == "Any annotation")],
                 NA,
                 main.values.annos[which(!(names(main.values.annos) %in% c("No annotations", "Any annotation")))])
# Vertical
pdf(paste(OUTDIR, "/", prefix, ".noncoding_APS.main_categories.long.pdf", sep=""),
    height=5.5, width=4)
plot.dots(main.values, parmar=c(0.5, 11.5, 2, 0.5),
          bg.lines=list(c(0, "gray50", 2),
                        c(baseline.aps[which(rownames(baseline.aps)=="coding_DEL"), 1], del.color, 2),
                        c(baseline.aps[which(rownames(baseline.aps)=="coding_DUP"), 1], dup.color, 2)))
dev.off()
# Horizontal
pdf(paste(OUTDIR, "/", prefix, ".noncoding_APS.main_categories.horiz.pdf", sep=""),
    height=3, width=7)
plot.dots(main.values, parmar=c(9, 7, 0.5, 0.5), horiz=T,
          bg.lines=list(c(0, "gray50", 2),
                        c(baseline.aps[which(rownames(baseline.aps)=="coding_DEL"), 1], del.color, 2),
                        c(baseline.aps[which(rownames(baseline.aps)=="coding_DUP"), 1], dup.color, 2)))
dev.off()


# ### Plot blowout subpanels (no longer used)
# pdf(paste(OUTDIR, "/", prefix, ".noncoding_APS.main_categories.blowouts_singlePanel.pdf", sep=""),
#     height=6, width=4)
# plot.dots(c(main.values.bench,
#             NA, NA,
#             get.blowout.vals(dat, "ANY"),
#             NA, NA,
#             get.blowout.vals(dat, "DHS"),
#             NA, NA,
#             get.blowout.vals(dat, "ENH")),
#           parmar=c(0.5, 11.5, 2, 0.5),
#           bg.lines=list(c(0, "black", 1),
#                         c(baseline.aps[which(rownames(baseline.aps)=="coding_DEL"), 1], del.color, 2),
#                         c(baseline.aps[which(rownames(baseline.aps)=="coding_DUP"), 1], dup.color, 2)))
# dev.off()

### Plot comparison of partial vs whole element CNVs
pdf(paste(OUTDIR, "/", prefix, ".noncoding_APS.whole_vs_partial.APS.pdf", sep=""),
    height=2.4, width=2.4)
plot.paired_cnvs(dat, mode="svi")
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".noncoding_APS.conserved_vs_unconserved.APS.pdf", sep=""),
    height=2.2, width=2.8)
plot.paired_cnvs(dat, mode="cons")
dev.off()

### Plot comparison of APS vs distance from nearest gene


### Plot correlation of phastCons vs APS
parmar <- c(2.25, 3, 0.5, 0.75)
pdf(paste(OUTDIR, "/", prefix, ".noncoding_APS.phastCons_corr.DEL.pdf", sep=""), 
    height=1.6, width=2)
plot.phastCons.dots(bed[intersect(noncoding.idxs, which(bed$SVTYPE=="DEL")), ], 
                    ylims=c(-0.05, 0.1), color=del.color, parmar=parmar, n.bins=100, tck=-0.05, 
                    cex.labels=0.75, cex.ax.titles=0.9, xlabel="phastCons Percentile")
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".noncoding_APS.phastCons_corr.DUP.pdf", sep=""), 
    height=1.6, width=2)
plot.phastCons.dots(bed[intersect(noncoding.idxs, which(bed$SVTYPE=="DUP")), ], 
                    ylims=c(-0.05, 0.1), color=dup.color, parmar=parmar, n.bins=100, tck=-0.05, 
                    cex.labels=0.75, cex.ax.titles=0.9, xlabel="phastCons Percentile")
dev.off()
