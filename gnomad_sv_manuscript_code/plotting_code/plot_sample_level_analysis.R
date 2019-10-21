#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Plot per-individual summary data (both sites and functional effects) for formal gnomAD analysis


###Set master parameters
options(stringsAsFactors=F, scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Gather waterfall data from an input data frame
getWaterfallData <- function(dat){
  #Get ancestry data
  dat$pop[which(dat$pop=="SAS")] <- "OTH"
  ancestries <- sort(unique(as.character(dat$pop)))
  ancestries <- c(ancestries[which(ancestries!="OTH")], 
                  ancestries[which(ancestries=="OTH")])
  #Get global & ancestry medians and means
  global.medians <- apply(dat[, -c(1:2)], 2, median, na.rm=T)
  ancestry.medians <- lapply(ancestries, function(anc){
    apply(dat[which(dat$pop==anc), -c(1:2)], 2, median, na.rm=T)
  })
  names(ancestry.medians) <- ancestries
  global.means <- apply(dat[, -c(1:2)], 2, mean, na.rm=T)
  ancestry.means <- lapply(ancestries, function(anc){
    apply(dat[which(dat$pop==anc), -c(1:2)], 2, mean, na.rm=T)
  })
  names(ancestry.means) <- ancestries
  #Get list of data frames, one df per ancestry, sorted by total N per sample
  ancestry.perSamp <- lapply(ancestries, function(anc){
    subdat <- dat[which(dat$pop==anc), -c(1:2)]
    sums <- apply(subdat, 1, sum, na.rm=T)
    subdat[order(-sums), ]
  })
  names(ancestry.perSamp) <- ancestries
  #Get number of samples per ancestry
  nsamp.perAncestry <- unlist(lapply(ancestry.perSamp, nrow))
  #Prep output list of values
  out.list <- list("global.medians"=global.medians, 
                   "global.means"=global.means, 
                   "ancestries"=ancestries, 
                   "nsamp.perAncestry"=nsamp.perAncestry, 
                   "ancestry.medians"=ancestry.medians, 
                   "ancestry.means"=ancestry.means, 
                   "ancestry.perSamp"=ancestry.perSamp)
  return(out.list)
}
#Plot waterfall data
plotWaterfall <- function(wdat, colors, ylabel=NULL, y.at.nmax=5, titles=F, 
                          MCNV=T, cat.labels=NULL, cex.ylabs=0.8, cex.ytitle=0.6,
                          median.adj=-0.025){
  #Set graphical parameters
  MCNV.dens <- 60
  
  #Drop extreme outliers for visualization purposes
  wdat$ancestry.perSamp <- lapply(wdat$ancestry.perSamp, function(subdat){
    totals <- apply(subdat, 1, sum, na.rm=T)
    q3 <- quantile(totals, 0.75)
    iqr <- IQR(totals)
    outliers <- which(totals>q3+(3*iqr))
    print(length(outliers)/nrow(subdat))
    if(length(outliers)>0){
      subdat <- subdat[-outliers, ]
    }
    return(subdat)
  })
  
  #Prep plot area & scaling
  n.samps <- sum(wdat$nsamp.perAncestry)
  x.buffer.small <- 0.02*n.samps
  x.buffer.large <- 0.05*n.samps
  y.max <- 1.05*max(unlist(lapply(wdat$ancestry.perSamp, function(d){apply(d, 1, sum, na.rm=T)})))
  n.ancestries <- length(wdat$ancestries)
  xlims <- c(0, n.samps+(n.ancestries*(x.buffer.large+(1.5*x.buffer.small)))+((2*x.buffer.large)+x.buffer.small))
  plot(x=xlims, y=c(0, y.max), type="n", 
       xaxt="n", xlab="", xaxs="i", yaxt="n", ylab="", yaxs="i")
  
  #Plot master median of all samples
  rect(xleft=xlims[2]-(x.buffer.large+x.buffer.small), 
       xright=xlims[2]-x.buffer.small, 
       ybottom=-par("usr")[4], ytop=sum(wdat$global.medians), 
       col=NA, lwd=2)
  # segments(x0=xlims[2]-x.buffer.small, x1=xlims[1]+x.buffer.large, 
  #          y0=sum(wdat$global.medians), 
  #          y1=sum(wdat$global.medians), 
  #          col="gray50")
  rect(xleft=xlims[2]-(x.buffer.large+x.buffer.small), 
       xright=xlims[2]-x.buffer.small, 
       ybottom=c(0, cumsum(wdat$global.medians)[-length(wdat$global.medians)]), 
       ytop=cumsum(wdat$global.medians), 
       bty="n", border=NA, col=colors)
  if(MCNV==T){
    rect(xleft=xlims[2]-(x.buffer.large+x.buffer.small), 
         xright=xlims[2]-x.buffer.small, 
         ybottom=c(0, cumsum(wdat$global.medians)[-length(wdat$global.medians)])[grep("MCNV", names(wdat$global.medians), fixed=T)], 
         ytop=cumsum(wdat$global.medians)[grep("MCNV", names(wdat$global.medians), fixed=T)], 
         bty="n", border=NA, lend="butt", lwd=0.5, col="white", density=MCNV.dens)
  }
  text(x=xlims[2]-((0.5*x.buffer.large)+x.buffer.small), 
       y=sum(wdat$global.medians)+(median.adj*(par("usr")[4]-par("usr")[3])), 
       font=2, pos=3, cex=0.7, xpd=T, 
       labels=prettyNum(round(sum(wdat$global.medians), 0), big.mark=","))
  #Add title, if optioned
  if(titles==T){
    axis(3, at=c(xlims[2]-(x.buffer.large+x.buffer.small), xlims[2]-x.buffer.small), 
         tck=0, labels=NA, line=0.25)
    axis(3, at=mean(c(xlims[2]-(x.buffer.large+x.buffer.small), xlims[2]-x.buffer.small)), 
         tick=F, labels="Median", cex.axis=0.7, xpd=T, padj=0, font=2, line=-0.7)
  }
  #Add connector to right axis
  right.axis.at <- seq(par("usr")[3], par("usr")[4], by=(par("usr")[4]-par("usr")[3])/length(wdat$global.medians[which(wdat$global.medians>0)]))
  right.axis.at <- right.axis.at-(0.5*(right.axis.at[2]-right.axis.at[1]))
  right.axis.at <- right.axis.at[-1]
  sapply(1:length(wdat$global.medians[which(wdat$global.medians>0)]), function(i){
    axis(4, at=right.axis.at[i], labels=NA, lwd=0.75, lend="butt")
    axis(4, at=right.axis.at[i], tick=F, cex.axis=0.7, las=2, line=-0.4, 
         labels=paste(prettyNum(round(wdat$global.medians[which(wdat$global.medians>0)][i], 0), big.mark=","), 
                      cat.labels[which(wdat$global.medians>0)][i], sep=" "), 
         col.axis=colors[which(wdat$global.medians>0)][i])
    axis(4, at=right.axis.at[i], tick=F, cex.axis=0.7, las=2, line=-0.4, 
         labels=prettyNum(round(wdat$global.medians[which(wdat$global.medians>0)][i], 0), big.mark=","), 
         col.axis="black")
  })
  segments(x0=rep(xlims[2]-x.buffer.small, length(wdat$global.medians[which(wdat$global.medians>0)])), 
           x1=rep(par("usr")[2], length(wdat$global.medians[which(wdat$global.medians>0)])), 
           y0=as.numeric((c(0, cumsum(wdat$global.medians)[-length(wdat$global.medians)])+cumsum(wdat$global.medians))/2)[which(wdat$global.medians>0)], 
           y1=right.axis.at, 
           lwd=0.75, lend="butt")
  
  #Iterate over ancestries & plot
  sapply(1:length(wdat$ancestries), function(i){
    #Determine starting coordinate
    if(i>1){
      left.pos <- sum(wdat$nsamp.perAncestry[1:(i-1)])
    }else{
      left.pos <- 0
    }
    left.pos <- left.pos+(x.buffer.large*i)+(1.5*x.buffer.small*i)
    #Plot individual sample bars
    sapply(1:nrow(wdat$ancestry.perSamp[[i]]), function(k){
      rect(xleft=left.pos+k-1, xright=left.pos+k, 
           ybottom=c(0, cumsum(as.numeric(wdat$ancestry.perSamp[[i]][k, ]))[-ncol(wdat$ancestry.perSamp[[i]])]), 
           ytop=cumsum(as.numeric(wdat$ancestry.perSamp[[i]][k, ])), 
           bty="n", border=NA, col=adjustcolor(colors, alpha=0.5))
      if(MCNV==T){
        rect(xleft=left.pos+k-1, xright=left.pos+k, 
             ybottom=c(0, cumsum(as.numeric(wdat$ancestry.perSamp[[i]][k, ]))[-ncol(wdat$ancestry.perSamp[[i]])])[grep("MCNV", names(wdat$ancestry.medians[[i]]), fixed=T)], 
             ytop=cumsum(as.numeric(wdat$ancestry.perSamp[[i]][k, ]))[grep("MCNV", names(wdat$ancestry.medians[[i]]), fixed=T)], 
             bty="n", border=NA, lend="butt", lwd=0.5, col="white", density=MCNV.dens)
      }
    })
    #Plot ancestry median
    segments(x0=left.pos-(1.5*x.buffer.small), x1=left.pos+wdat$nsamp.perAncestry[i], 
             y0=sum(wdat$ancestry.medians[[i]]), y1=sum(wdat$ancestry.medians[[i]]), 
             col="gray40", lwd=0.75)
    rect(xleft=left.pos-(1.5*x.buffer.small), xright=left.pos-(0.5*x.buffer.small), 
         ybottom=-par("usr")[4], ytop=sum(wdat$ancestry.medians[[i]]), 
         col=NA, lwd=2)
    rect(xleft=left.pos-(1.5*x.buffer.small), xright=left.pos-(0.5*x.buffer.small), 
         ybottom=c(0, cumsum(wdat$ancestry.medians[[i]])[-length(wdat$ancestry.medians[[i]])]), 
         ytop=cumsum(wdat$ancestry.medians[[i]]), 
         bty="n", border=NA, col=colors)
    if(MCNV==T){
      rect(xleft=left.pos-(1.5*x.buffer.small), xright=left.pos-(0.5*x.buffer.small), 
           ybottom=c(0, cumsum(wdat$ancestry.medians[[i]])[-length(wdat$ancestry.medians[[i]])])[grep("MCNV", names(wdat$ancestry.medians[[i]]), fixed=T)], 
           ytop=cumsum(wdat$ancestry.medians[[i]])[grep("MCNV", names(wdat$ancestry.medians[[i]]), fixed=T)], 
           bty="n", border=NA, lend="butt", lwd=0.5, col="white", density=MCNV.dens)
    }
    text(x=left.pos-(1.5*x.buffer.small), 
         y=sum(wdat$ancestry.medians[[i]])+(median.adj*(par("usr")[4]-par("usr")[3])), 
         pos=3, cex=0.7, xpd=T, 
         labels=prettyNum(round(median(apply(wdat$ancestry.perSamp[[i]], 1, sum)), 0), big.mark=","))
    #Add title, if optioned
    if(titles==T){
      axis(3, at=c(left.pos-(1.5*x.buffer.small), left.pos+wdat$nsamp.perAncestry[i]), tck=0, labels=NA, line=0.25)
      # axis(3, at=mean(c(left.pos-(1.5*x.buffer.small), left.pos+wdat$nsamp.perAncestry[i])), 
      #      tick=F, labels=paste(pops$name[which(pops$pop==wdat$ancestries[i])], 
      #                          "\n(n=", prettyNum(wdat$nsamp.perAncestry[i], big.mark=","), ")", sep=""), 
      #      cex.axis=0.7, xpd=T, line=-0.8, padj=0)
      axis(3, at=mean(c(left.pos-(1.5*x.buffer.small), left.pos+wdat$nsamp.perAncestry[i])), 
           tick=F, labels=pops$pop[which(pops$pop==wdat$ancestries[i])], 
           cex.axis=0.7, xpd=T, line=-0.7, padj=0)
    }
  })
  #Clean up
  axis(1, at=c(par("usr")[1], par("usr")[2]), labels=NA, tck=0)
  if(length(axTicks(2))>y.at.nmax){
    y.at <- axTicks(2)[seq(1, length(axTicks(2)), 2)]
    y.at <- c(y.at, y.at[length(y.at)]+(y.at[2]-y.at[1]))
  }else{
    y.at <- axTicks(2)
  }
  axis(2, at=y.at, labels=NA)
  axis(2, at=y.at, tick=F, labels=prettyNum(y.at, big.mark=","), 
       line=-0.4, cex.axis=cex.ylabs, las=2)
  mtext(2, line=2.5, text=ylabel, cex=cex.ytitle)
}

#Generate jitter residuals for a vector of values based on their density
sina.jitter <- function(vals){
  d <- density(vals)
  dv <- approx(d$x, d$y, xout=vals)
  dv <- dv$y/max(dv$y)
  dv.j <- sapply(1:length(vals), function(i){
    jitter(0, amount=dv[i])
  })
  return(dv.j)
}

#Generate sina points to add to existing plot
sina.plot <- function(vals, y.at, color, width=0.1, horiz=T, cex=0.25){
  j <- (width*sina.jitter(vals))+y.at
  if(horiz==T){
    points(x=vals, y=j, pch=21, cex=cex, col=color, bg=adjustcolor(color, alpha=0.3))
    segments(x0=median(vals), x1=median(vals), 
             y0=y.at-width, y1=y.at+width, lwd=3, lend="butt")
  }else{
    points(x=j, y=vals, pch=21, cex=cex, col=color, bg=adjustcolor(color, alpha=0.3))
    segments(x0=y.at-width, x1=y.at+width, y0=median(vals), y1=median(vals), lwd=3, lend="butt")
  }
}

#Plot single panel vertical sina plot for a given svtype, by ancestry
plot.sina.singleSvtypeByPop <- function(wdat, svtype, ordered.pops, title=NULL){
  #Get plot values
  plot.vals <- lapply(ordered.pops, function(pop){
      subdf <- wdat$ancestry.perSamp[[which(names(wdat$ancestry.perSamp)==pop)]]
      if(svtype!="ALL"){
        subdf[, grep(paste(".", svtype, sep=""), colnames(subdf), fixed=T)]
      }else{
        apply(subdf, 1, sum, na.rm=T)
      }
  })
  if(is.null(title)){
    title <- svtype
  }
  #Prep plot area
  par(mar=c(0.5, 2, 1.5, 0.5), bty="n")
  plot(x=c(0, length(plot.vals)), y=range(unlist(plot.vals), na.rm=T), 
       type="n", xlab="", xaxt="n", ylab="", yaxt="n")  
  #Add jitter & color by pop
  sapply(1:length(ordered.pops), function(i){
    pop <- ordered.pops[i]
    pop.col <- pops$color[which(pops$pop==pop)]
    sina.plot(vals=plot.vals[[i]], y.at=i-0.5, color=pop.col, 
              width=0.4, horiz=F, cex=0.1)
  })
  #Add axes
  axis(3, at=c(0, length(plot.vals)), tck=0, labels=NA)
  mtext(3, line=0.2, text=title, cex=0.7)
  y.ax.at <- axTicks(2)
  if(length(y.ax.at)>5){
    y.ax.at <- y.ax.at[seq(1, length(y.ax.at), 2)]
  }
  axis(2, at=y.ax.at, labels=NA)
  if(max(y.ax.at)>1000){
    axis(2, at=y.ax.at, tick=F, 
         labels=paste(round(y.ax.at/1000, 2), "k", sep=""), 
         las=2, cex.axis=0.8, line=-0.4)
  }else{
    axis(2, at=y.ax.at, tick=F, labels=prettyNum(y.ax.at, big.mark=","), 
         las=2, cex.axis=0.8, line=-0.4)
  }
}

#Wrapper to plot horizontal grid of counts of SV per population by SVTYPE
wrapper.gridSinaSvtypeByPop <- function(wdat){
  par(mfrow=c(1, 8))
  sapply(c("ALL", "DEL", "DUP", "MCNV_DEL", "MCNV_DUP", "INS", "INV", "CPX"), function(svtype){
    if(svtype=="ALL"){
      title <- "All SVs"
    }else if(svtype=="MCNV_DEL"){
      title <- "MCNV (Loss)"
    }else if(svtype=="MCNV_DUP"){
      title <- "MCNV (Gain)"
    }else{
      title <- svtype
    }
    plot.sina.singleSvtypeByPop(wdat=wdat, svtype=svtype, 
                                ordered.pops=pops$pop, 
                                title=title)
  })
}


#Plot single panel vertical sina plot for a set of functional category, restricted by population
plot.sina.MultiFunc <- function(wdat, funcs, func.labels, func.cols, 
                                pop="ALL", ymax=NULL, boxes=F, vline=NULL, 
                                pt.labels=F){
  # #DEV:
  # wdat <- all.func.wdat.genes
  # funcs <- c("PROTEIN_CODING__LOF.genes", 
  #            "PROTEIN_CODING__DUP_LOF.genes", 
  #            "PROTEIN_CODING__COPY_GAIN.genes", 
  #            "PROTEIN_CODING__MCNV_LOSS.genes", 
  #            "PROTEIN_CODING__MCNV_GAIN.genes")
  # func.labels <- c("pLoF", "IED", "CG", "MCNV\n(pLoF)", "MCNV\n(CG)")
  # func.cols <- c(svtypes$color[which(svtypes$svtype=="DEL")], 
  #                svtypes$color[which(svtypes$svtype=="MCNV")], 
  #                svtypes$color[which(svtypes$svtype=="DUP")], 
  #                svtypes$color[which(svtypes$svtype=="DEL")], 
  #                svtypes$color[which(svtypes$svtype=="DUP")])
  # pop <- "EUR"
  #Get plot values
  if(pop=="ALL"){
    d.tmp <- do.call("rbind", wdat$ancestry.perSamp)
    plot.dat <- as.data.frame(do.call("cbind", lapply(funcs, function(func){
      d.tmp[, grep(func, colnames(d.tmp), fixed=T)]
    })))
  }else{
    d.tmp <- wdat$ancestry.perSamp[[which(names(wdat$ancestry.perSamp)==pop)]]
    plot.dat <- as.data.frame(do.call("cbind", lapply(funcs, function(func){
      d.tmp[, grep(func, colnames(d.tmp), fixed=T)]
    })))
  }
  colnames(plot.dat) <- funcs
  if(is.null(ymax)){
    ymax <- max(apply(plot.dat, 2, quantile, probs=0.999))
  }
  #Prep plot area
  par(mar=c(2.5, 2.75, 0.5, 0.5), bty="n")
  plot(x=c(0.1, ncol(plot.dat)+0.25), y=c(0, ymax), 
       type="n", xlab="", xaxt="n", ylab="", yaxt="n", yaxs="i")
  if(!is.null(vline)){
    abline(v=vline, lty=2, col="gray80")
  }
  # axis(1, at=c(-100, 100))
  #Add jitter & color by class
  sapply(1:ncol(plot.dat), function(i){
    if(boxes==T){
     boxplot(plot.dat[, i], at=i-0.5, col=func.cols[i], 
             outline=F, add=T, lty=1, staplewex=0, 
             xaxt="n", yaxt="n", xlab="", ylab="")
      if(pt.labels==T){
        text(x=i-0.625, y=median(plot.dat[, i]), pos=4, cex=0.65, xpd=T, 
             labels=prettyNum(median(plot.dat[, i]), big.mark=","))
      }
    }else{
      sina.plot(vals=plot.dat[, i], y.at=i-0.5, color=func.cols[i], 
                width=0.3, horiz=F, cex=0.025)
      if(pt.labels==T){
        text(x=i-0.525, y=median(plot.dat[, i]), pos=4, cex=0.65, xpd=T, 
             labels=prettyNum(median(plot.dat[, i]), big.mark=","))
      }
    }
    axis(1, at=i-0.5, las=2, labels=func.labels[i], srt=2, tick=F, cex.axis=0.85, line=-0.8)
  })
  #Add axes
  y.ax.at <- axTicks(2)
  if(length(y.ax.at)>5){
    y.ax.at <- y.ax.at[seq(1, length(y.ax.at), 2)]
  }
  axis(2, at=c(-100, 10000))
  axis(2, at=y.ax.at, labels=NA)
  if(max(y.ax.at)>1000){
    axis(2, at=y.ax.at, tick=F, 
         labels=paste(round(y.ax.at/1000, 2), "k", sep=""), 
         las=2, cex.axis=0.8, line=-0.4)
  }else{
    axis(2, at=y.ax.at, tick=F, labels=prettyNum(y.ax.at, big.mark=","), 
         las=2, cex.axis=0.8, line=-0.4)
  }
  if(pop!="ALL"){
    mtext(2, line=1.9, text=paste("Genes per", pop, "Genome", sep=" "))
  }else{
    mtext(2, line=1.9, text="Genes per Genome")
  }
}


################
###RSCRIPT BLOCK
################
require(optparse, quietly=T)
require(zoo, quietly=T)
###List of command-line options
option_list <- list(
  make_option(c("--prefix"), type="character", default="gnomAD_v2_SV", 
              help="prefix used for naming outfiles [default %default]", 
              metavar="character"), 
  make_option(c("-S", "--svtypes"), type="character", default=NULL, 
              help="tab-delimited file specifying SV types and HEX colors [default %default]", 
              metavar="character"), 
  make_option(c("-P", "--populations"), type="character", default=NULL, 
              help="tab-delimited file specifying populations and HEX colors [default %default]", 
              metavar="character")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog SITES_SUMMARY FUNCTIONAL_SUMMARY OUTDIR", 
                                option_list=option_list), 
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 5){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
sites.in <- args$args[1]
alleles.in <- args$args[2]
functional.in <- args$args[3]
PCRMINUS.samples.list.in <- args$args[4]
OUTDIR <- args$args[5]
prefix <- opts$prefix
svtypes.file <- opts$svtypes
pops.file <- opts$populations

# #Dev parameters (local)
# sites.in <- "~/scratch/gnomAD-SV_v2_rev1.sites.perSample_summary.txt.gz"
# alleles.in <- "~/scratch/gnomAD-SV_v2_rev1.alleles.perSample_summary.txt.gz"
# functional.in <- "~/scratch/gnomAD-SV_v2_rev1.functional.perSample_summary.txt.gz"
# PCRMINUS.samples.list.in <- "~/scratch/PCRMINUS_samples.list"
# OUTDIR <- "~/scratch/gnomAD_v2_test_plots/"
# prefix <- "gnomAD-SV_v2_rev1"
# svtypes.file <- "~/Desktop/Collins/Talkowski/code/sv-pipeline/ref/vcf_qc_refs/SV_colors.txt"
# pops.file <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/gnomAD_population_colors.txt"

###Create output directory, if necessary
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

###Process input data
cat("NOW LOADING DATA\n")
sites <- read.table(sites.in, header=T, sep="\t")
alleles <- read.table(alleles.in, header=T, sep="\t")
functional <- read.table(functional.in, header=T, sep="\t")
PCRMINUS.samples <- as.character(read.table(PCRMINUS.samples.list.in, header=F)[, 1])

###Restrict to PCRMINUS samples
sites <- sites[which(sites$sample %in% PCRMINUS.samples), ]
alleles <- alleles[which(alleles$sample %in% PCRMINUS.samples), ]
functional <- functional[which(functional$sample %in% PCRMINUS.samples), ]

###Sets sv types & colors
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

###Sets populations & colors
if(!is.null(pops.file)){
  pops <- read.table(pops.file, sep="\t", header=T, comment.char="", check.names=F)
}


###Get count of sites for AFR vs non-AFR
median(sites$all.VARIANTS[which(sites$pop=="AFR")])
median(sites$all.VARIANTS[which(sites$pop!="AFR")])
###Get count of homozygous sites for EAS vs non-EAS
median(sites$homozygous.VARIANTS[which(sites$pop=="EAS")])
median(sites$homozygous.VARIANTS[which(sites$pop!="EAS")])


###Waterfall plots of variant counts per sample
site.colors <- c(svtypes$color[which(svtypes$svtype=="DEL")], 
                 svtypes$color[which(svtypes$svtype=="DEL")], 
                 svtypes$color[which(svtypes$svtype=="DUP")], 
                 svtypes$color[which(svtypes$svtype=="DUP")], 
                 svtypes$color[which(svtypes$svtype=="INS")], 
                 svtypes$color[which(svtypes$svtype=="INV")], 
                 svtypes$color[which(svtypes$svtype=="CPX")])
site.labels <- c("DEL", "MCNV (Loss)", "DUP", "MCNV (Gain)", 
                 "INS", "INV", "CPX")
# site.colors <- sapply(site.colors, adjustcolor, alpha=0.75)
all.sites.wdat <- getWaterfallData(dat=sites[, which(colnames(sites) %in% c("sample", "pop", 
                                                                           "all.DEL", "all.MCNV_DEL", 
                                                                           "all.DUP", "all.MCNV_DUP", 
                                                                           "all.INS", "all.INV", "all.CPX"))])
homozygous.sites.wdat <- getWaterfallData(dat=sites[, which(colnames(sites) %in% c("sample", "pop", 
                                                                                  "homozygous.DEL", "homozygous.MCNV_DEL", 
                                                                                  "homozygous.DUP", "homozygous.MCNV_DUP", 
                                                                                  "homozygous.INS", "homozygous.INV", "homozygous.CPX"))])
rare.sites.wdat <- getWaterfallData(dat=sites[, which(colnames(sites) %in% c("sample", "pop", 
                                                                            "rare.DEL", "rare.MCNV_DEL", 
                                                                            "rare.DUP", "rare.MCNV_DUP", 
                                                                            "rare.INS", "rare.INV", "rare.CPX"))])
singleton.sites.wdat <- getWaterfallData(dat=sites[, which(colnames(sites) %in% c("sample", "pop", 
                                                                                 "singleton.DEL", "singleton.MCNV_DEL", 
                                                                                 "singleton.DUP", "singleton.MCNV_DUP", 
                                                                                 "singleton.INS", "singleton.INV", "singleton.CPX"))])
#Generate waterfall plots
png(paste(OUTDIR, "/", prefix, ".sites_perSample.waterfall.png", sep=""), 
    width=9*300, height=4*300, res=400)
layout(matrix(1:4, nrow=4), heights=c(3, 2, 1.5, 1.5))
par(mar=c(0.5, 3.5, 1.75, 5.5), bty="n")
plotWaterfall(wdat=all.sites.wdat, colors=site.colors, 
              ylabel="All", titles=T, cat.labels=site.labels)
par(mar=c(0.5, 3.5, 0.5, 5.5), bty="n")
plotWaterfall(wdat=homozygous.sites.wdat, colors=site.colors, 
              ylabel="Homozygous", y.at.nmax=4, cat.labels=site.labels)
plotWaterfall(wdat=rare.sites.wdat, colors=site.colors, 
              ylabel="Rare", y.at.nmax=3, cat.labels=site.labels)
plotWaterfall(wdat=singleton.sites.wdat, colors=site.colors, 
              ylabel="Singleton", y.at.nmax=3, cat.labels=site.labels)
dev.off()

#Generate smaller plot of just counts per sample for main figure
png(paste(OUTDIR, "/", prefix, ".sites_perSample.waterfall.small.png", sep=""), 
    width=10.5*300, height=1.5*300, res=400)
par(mar=c(0.5, 3.5, 1.1, 5.5), bty="n")
plotWaterfall(wdat=all.sites.wdat, colors=site.colors, 
              ylabel="SVs per Genome", titles=T, cat.labels=site.labels, 
              y.at.nmax=10, cex.ylabs=0.6, median.adj=-0.06, cex.ytitle=0.8)
dev.off()

#Generate sina plot of counts per class by ancestry
png(paste(OUTDIR, "/", prefix, ".sites_perSample.sina.small.png", sep=""), 
    width=8*300, height=1.1*300, res=400)
wrapper.gridSinaSvtypeByPop(wdat=all.sites.wdat)
dev.off()

# ###Waterfall plots of allele counts per sample
# allele.labels <- site.labels
# allele.colors <- site.colors
# # site.colors <- sapply(site.colors, adjustcolor, alpha=0.75)
# all.alleles.wdat <- getWaterfallData(dat=alleles[, which(colnames(alleles) %in% c("sample", "pop", 
#                                                                            "all.DEL", "all.MCNV_DEL", 
#                                                                            "all.DUP", "all.MCNV_DUP", 
#                                                                            "all.INS", "all.INV", "all.CPX"))])
# homozygous.alleles.wdat <- getWaterfallData(dat=alleles[, which(colnames(alleles) %in% c("sample", "pop", 
#                                                                                   "homozygous.DEL", "homozygous.MCNV_DEL", 
#                                                                                   "homozygous.DUP", "homozygous.MCNV_DUP", 
#                                                                                   "homozygous.INS", "homozygous.INV", "homozygous.CPX"))])
# rare.alleles.wdat <- getWaterfallData(dat=alleles[, which(colnames(alleles) %in% c("sample", "pop", 
#                                                                             "rare.DEL", "rare.MCNV_DEL", 
#                                                                             "rare.DUP", "rare.MCNV_DUP", 
#                                                                             "rare.INS", "rare.INV", "rare.CPX"))])
# singleton.alleles.wdat <- getWaterfallData(dat=alleles[, which(colnames(alleles) %in% c("sample", "pop", 
#                                                                                  "singleton.DEL", "singleton.MCNV_DEL", 
#                                                                                  "singleton.DUP", "singleton.MCNV_DUP", 
#                                                                                  "singleton.INS", "singleton.INV", "singleton.CPX"))])
# #Generate plot
# png(paste(OUTDIR, "/", prefix, ".alleles_perSample.waterfall.png", sep=""), 
#     width=9*300, height=4*300, res=400)
# layout(matrix(1:4, nrow=4), heights=c(3, 2, 1.5, 1.5))
# par(mar=c(0.5, 3.5, 1.75, 5.5), bty="n")
# plotWaterfall(wdat=all.alleles.wdat, colors=allele.colors, 
#               ylabel="All", titles=T, cat.labels=allele.labels)
# par(mar=c(0.5, 3.5, 0.5, 5.5), bty="n")
# plotWaterfall(wdat=homozygous.alleles.wdat, colors=allele.colors, 
#               ylabel="Homozygous", y.at.nmax=4, cat.labels=allele.labels)
# plotWaterfall(wdat=rare.alleles.wdat, colors=allele.colors, 
#               ylabel="Rare", y.at.nmax=3, cat.labels=allele.labels)
# plotWaterfall(wdat=singleton.alleles.wdat, colors=allele.colors, 
#               ylabel="Singleton", y.at.nmax=3, cat.labels=allele.labels)
# dev.off()



###Waterfall plots of functional counts per sample
func.colors <- c(svtypes$color[which(svtypes$svtype=="DEL")], 
                 svtypes$color[which(svtypes$svtype=="MCNV")], 
                 svtypes$color[which(svtypes$svtype=="DUP")], 
                 svtypes$color[which(svtypes$svtype=="DEL")], 
                 svtypes$color[which(svtypes$svtype=="DUP")])
func.labels <- c("pLoF", "IED", "CG", "pLoF (MCNV)", "CG (MCNV)")
# all.func.wdat <- getWaterfallData(dat=functional[, which(colnames(functional) %in% c("sample", "pop", 
#                                                                                     "all.PROTEIN_CODING__LOF", 
#                                                                                     "all.PROTEIN_CODING__MCNV_LOSS", 
#                                                                                     "all.PROTEIN_CODING__COPY_GAIN", 
#                                                                                     "all.PROTEIN_CODING__MCNV_GAIN"))])
all.func.wdat.genes <- getWaterfallData(dat=functional[, which(colnames(functional) %in% c("sample", "pop", 
                                                                                          "all.PROTEIN_CODING__LOF.genes", 
                                                                                          "all.PROTEIN_CODING__DUP_LOF.genes", 
                                                                                          "all.PROTEIN_CODING__COPY_GAIN.genes", 
                                                                                          "all.PROTEIN_CODING__MCNV_LOSS.genes", 
                                                                                          "all.PROTEIN_CODING__MCNV_GAIN.genes"))])
# homozygous.func.wdat <- getWaterfallData(dat=functional[, which(colnames(functional) %in% c("sample", "pop", 
#                                                                                            "homozygous.PROTEIN_CODING__LOF", 
#                                                                                            "homozygous.PROTEIN_CODING__MCNV_LOSS", 
#                                                                                            "homozygous.PROTEIN_CODING__COPY_GAIN", 
#                                                                                            "homozygous.PROTEIN_CODING__MCNV_GAIN"))])
homozygous.func.wdat.genes <- getWaterfallData(dat=functional[, which(colnames(functional) %in% c("sample", "pop", 
                                                                                                 "homozygous.PROTEIN_CODING__LOF.genes", 
                                                                                                 "homozygous.PROTEIN_CODING__DUP_LOF.genes", 
                                                                                                 "homozygous.PROTEIN_CODING__COPY_GAIN.genes", 
                                                                                                 "homozygous.PROTEIN_CODING__MCNV_LOSS.genes", 
                                                                                                 "homozygous.PROTEIN_CODING__MCNV_GAIN.genes"))])
# rare.func.wdat <- getWaterfallData(dat=functional[, which(colnames(functional) %in% c("sample", "pop", 
#                                                                                      "rare.PROTEIN_CODING__LOF", 
#                                                                                      "rare.PROTEIN_CODING__MCNV_LOSS", 
#                                                                                      "rare.PROTEIN_CODING__COPY_GAIN", 
#                                                                                      "rare.PROTEIN_CODING__MCNV_GAIN"))])
rare.func.wdat.genes <- getWaterfallData(dat=functional[, which(colnames(functional) %in% c("sample", "pop", 
                                                                                           "rare.PROTEIN_CODING__LOF.genes", 
                                                                                           "rare.PROTEIN_CODING__DUP_LOF.genes", 
                                                                                           "rare.PROTEIN_CODING__COPY_GAIN.genes", 
                                                                                           "rare.PROTEIN_CODING__MCNV_LOSS.genes", 
                                                                                           "rare.PROTEIN_CODING__MCNV_GAIN.genes"))])
singleton.func.wdat.genes <- getWaterfallData(dat=functional[, which(colnames(functional) %in% c("sample", "pop", 
                                                                                           "singleton.PROTEIN_CODING__LOF.genes", 
                                                                                           "singleton.PROTEIN_CODING__DUP_LOF.genes", 
                                                                                           "singleton.PROTEIN_CODING__COPY_GAIN.genes", 
                                                                                           "singleton.PROTEIN_CODING__MCNV_LOSS.genes", 
                                                                                           "singleton.PROTEIN_CODING__MCNV_GAIN.genes"))])
lof.func.wdat <- getWaterfallData(dat=functional[, which(colnames(functional) %in% c("sample", "pop", 
                                                                                    "lof.DEL.variants", 
                                                                                    "lof.MCNV_DEL.variants", 
                                                                                    "lof.DUP.variants", 
                                                                                    "lof.INS.variants", 
                                                                                    "lof.INV.variants", 
                                                                                    "lof.CPX.variants"))])
lof.colors <- c(svtypes$color[which(svtypes$svtype=="DEL")], 
                svtypes$color[which(svtypes$svtype=="DUP")], 
                svtypes$color[which(svtypes$svtype=="MCNV")], 
                svtypes$color[which(svtypes$svtype=="INS")], 
                svtypes$color[which(svtypes$svtype=="INV")], 
                svtypes$color[which(svtypes$svtype=="CPX")])
lof.labels <- c("DEL", "DUP", "MCNV", "INS", "INV", "CPX")
#Generate plot
png(paste(OUTDIR, "/", prefix, ".coding_perSample.waterfall.png", sep=""), 
    width=9*300, height=4.25*300, res=400)
layout(matrix(1:4, nrow=4), heights=c(2, 1, 0.75, 1.5))
par(mar=c(0.5, 3.5, 1.75, 5.5), bty="n")
plotWaterfall(wdat=all.func.wdat.genes, colors=func.colors, 
              ylabel="All", y.at.nmax=10, MCNV=T, titles=T, 
              cat.labels=func.labels)
par(mar=c(0.5, 3.5, 0.5, 5.5), bty="n")
plotWaterfall(wdat=homozygous.func.wdat.genes, colors=func.colors, 
              ylabel="Homozygous", y.at.nmax=10, MCNV=T, 
              cat.labels=func.labels)
plotWaterfall(wdat=rare.func.wdat.genes, colors=func.colors, 
              ylabel="Rare", y.at.nmax=4, MCNV=T, 
              cat.labels=func.labels)
plotWaterfall(wdat=lof.func.wdat, colors=lof.colors, 
              ylabel="pLoF SVs", y.at.nmax=4, MCNV=F, 
              cat.labels=lof.labels)
dev.off()


#Generate smaller plots of just counts per sample for main figure
func.colors.mod <- c(svtypes$color[which(svtypes$svtype=="DEL")], 
                 svtypes$color[which(svtypes$svtype=="MCNV")], 
                 svtypes$color[which(svtypes$svtype=="DUP")], 
                 svtypes$color[which(svtypes$svtype=="DEL")], 
                 svtypes$color[which(svtypes$svtype=="DUP")])
func.labels <- c("pLoF", "IED", "CG", "pLoF (MCNV)", "CG (MCNV)")
png(paste(OUTDIR, "/", prefix, ".genes_perSample.waterfall.small_long.png", sep=""), 
    width=10*300, height=1.75*300, res=400)
par(mar=c(0.5, 3.5, 1.75, 5.5), bty="n")
plotWaterfall(wdat=all.func.wdat.genes, colors=func.colors.mod, 
              ylabel="Genes per Genome", titles=T, cat.labels=func.labels)
dev.off()
png(paste(OUTDIR, "/", prefix, ".genes_perSample.rare.waterfall.small_long.png", sep=""), 
    width=10*300, height=1.75*300, res=400)
par(mar=c(0.5, 3.5, 1.75, 5.5), bty="n")
plotWaterfall(wdat=rare.func.wdat.genes, colors=func.colors.mod, 
              ylabel="Genes per Genome", titles=T, cat.labels=func.labels)
dev.off()


#Sina plots of counts of genes disrupted per genome by functional class
#All SV
png(paste(OUTDIR, "/", prefix, ".genes_perSample.all_sv.sina.png", sep=""), 
    width=3*300, height=2.9*300, res=400)
plot.sina.MultiFunc(wdat=all.func.wdat.genes, 
                    funcs=c("PROTEIN_CODING__LOF.genes", 
                            "PROTEIN_CODING__COPY_GAIN.genes", 
                            "PROTEIN_CODING__DUP_LOF.genes", 
                            "PROTEIN_CODING__MCNV_LOSS.genes", 
                            "PROTEIN_CODING__MCNV_GAIN.genes"), 
                    func.labels=c("pLoF", "CG", "IED", "pLoF", "CG"), 
                    func.cols=c(svtypes$color[which(svtypes$svtype=="DEL")], 
                                svtypes$color[which(svtypes$svtype=="DUP")], 
                                svtypes$color[which(svtypes$svtype=="MCNV")], 
                                svtypes$color[which(svtypes$svtype=="DEL")], 
                                svtypes$color[which(svtypes$svtype=="DUP")]), 
                    pop="ALL", boxes=F, ymax=275, 
                    vline=3)
axis(1, at=3, col.ticks="gray80", labels=NA, tck=-0.2)
points(x=0.5, y=122.4, xpd=T, pch=23, cex=1, lwd=1.7, bg="white")
dev.off()
#Rare SV
png(paste(OUTDIR, "/", prefix, ".genes_perSample.rare_sv.sina.png", sep=""), 
    width=2*300, height=2.9*300, res=400)
plot.sina.MultiFunc(wdat=rare.func.wdat.genes, 
                    funcs=c("PROTEIN_CODING__LOF.genes", 
                            "PROTEIN_CODING__COPY_GAIN.genes",
                            "PROTEIN_CODING__DUP_LOF.genes"), 
                    func.labels=c("pLoF", "CG", "IED"), 
                    func.cols=c(svtypes$color[which(svtypes$svtype=="DEL")], 
                                svtypes$color[which(svtypes$svtype=="DUP")], 
                                svtypes$color[which(svtypes$svtype=="MCNV")]), 
                    pop="ALL", ymax=18, boxes=T)
points(x=0.5, y=16.3, xpd=T, pch=23, cex=1, lwd=1.7, bg="white")
# axis(2, at=18, lwd=2, labels=NA)
dev.off()


###Write out tables of global and per-population medians
#Table of sites
sites.table.out <- rbind(all.sites.wdat$global.medians, 
                         do.call("rbind", all.sites.wdat$ancestry.medians), 
                         homozygous.sites.wdat$global.medians, 
                         do.call("rbind", homozygous.sites.wdat$ancestry.medians), 
                         rare.sites.wdat$global.medians, 
                         do.call("rbind", rare.sites.wdat$ancestry.medians), 
                         singleton.sites.wdat$global.medians, 
                         do.call("rbind", singleton.sites.wdat$ancestry.medians))
colnames(sites.table.out) <- gsub("all.", "", colnames(sites.table.out), fixed=T)
sites.table.out <- cbind(apply(sites.table.out, 1, sum), sites.table.out)
colnames(sites.table.out)[1] <- "ALL"
rownames(sites.table.out) <- NULL
pop.col <- as.character(sapply(c("all", "homozygous", "rare", "singleton"), function(prefix){paste(prefix, c("ALL", all.sites.wdat$ancestries), sep=".")}))
pop.N <- rep(c(sum(all.sites.wdat$nsamp.perAncestry), 
           all.sites.wdat$nsamp.perAncestry), 4)
sites.table.out <- data.frame("pop"=pop.col, 
                              "pop.N"=pop.N, 
                              sites.table.out)
write.table(sites.table.out, paste(OUTDIR, "/", prefix, ".median_sites_per_sample_by_population.txt", sep=""), 
            col.names=T, row.names=F, sep="\t", quote=F)

#Table of functional data
func.table.out <- rbind(all.func.wdat.genes$global.medians, 
                         do.call("rbind", all.func.wdat.genes$ancestry.medians), 
                         homozygous.func.wdat.genes$global.medians, 
                         do.call("rbind", homozygous.func.wdat.genes$ancestry.medians), 
                         rare.func.wdat.genes$global.medians, 
                         do.call("rbind", rare.func.wdat.genes$ancestry.medians))
colnames(func.table.out) <- gsub("all.", "", colnames(func.table.out), fixed=T)
func.table.out <- cbind(apply(func.table.out, 1, sum), func.table.out)
colnames(func.table.out)[1] <- "ALL"
rownames(func.table.out) <- NULL
pop.col <- as.character(sapply(c("all", "homozygous", "rare"), function(prefix){paste(prefix, c("ALL", all.func.wdat.genes$ancestries), sep=".")}))
pop.N <- rep(c(sum(all.func.wdat.genes$nsamp.perAncestry), 
               all.func.wdat.genes$nsamp.perAncestry), 3)
func.table.out <- data.frame("pop"=pop.col, 
                              "pop.N"=pop.N, 
                              func.table.out)
write.table(func.table.out, paste(OUTDIR, "/", prefix, ".median_genes_per_sample_by_population.txt", sep=""), 
            col.names=T, row.names=F, sep="\t", quote=F)
