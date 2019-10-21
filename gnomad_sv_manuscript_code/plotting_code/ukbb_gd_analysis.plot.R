#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Plot genomic disorder frequency comparison between UKBB and gnomAD


###Set master parameters
options(stringsAsFactors=F, scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Read and postprocess input data
import.data <- function(data.in){
  #Read data
  dat <- read.table(data.in, header=F, sep="\t")
  colnames(dat) <- c("chr", "start", "end", "ID", "CNV", "UKBB.N", "UKBB.k", "gnomAD.N", "gnomAD.k", 
                     "Coe_CASE.N", "Coe_CASE.k", "Coe_CTRL.N", "Coe_CTRL.k")
  #Calculate mean & 95% CI
  binom.stats <- function(N, k){
    d <- binom.test(x=k, n=N)
    return(as.numeric(c(d$estimate, d$conf.int)))
  }
  UKBB.stats <- do.call("rbind", lapply(1:nrow(dat), function(i){
    binom.stats(N=dat$UKBB.N[i], k=dat$UKBB.k[i])
  }))
  colnames(UKBB.stats) <- c("UKBB.freq", "UKBB.lower", "UKBB.upper")
  gnomAD.stats <- do.call("rbind", lapply(1:nrow(dat), function(i){
    binom.stats(N=dat$gnomAD.N[i], k=dat$gnomAD.k[i])
  }))
  colnames(gnomAD.stats) <- c("gnomAD.freq", "gnomAD.lower", "gnomAD.upper")
  combined.stats <- do.call("rbind", lapply(1:nrow(dat), function(i){
    binom.stats(N=dat$gnomAD.N[i]+dat$UKBB.N[i], k=dat$gnomAD.k[i]+dat$UKBB.k[i])
  }))
  colnames(combined.stats) <- c("combined.freq", "combined.lower", "combined.upper")
  Coe_CTRL.stats <- do.call("rbind", lapply(1:nrow(dat), function(i){
    binom.stats(N=dat$Coe_CTRL.N[i], k=dat$Coe_CTRL.k[i])
  }))
  colnames(Coe_CTRL.stats) <- c("Coe_CTRL.freq", "Coe_CTRL.lower", "Coe_CTRL.upper")
  Coe_CASE.stats <- do.call("rbind", lapply(1:nrow(dat), function(i){
    binom.stats(N=dat$Coe_CASE.N[i], k=dat$Coe_CASE.k[i])
  }))
  colnames(Coe_CASE.stats) <- c("Coe_CASE.freq", "Coe_CASE.lower", "Coe_CASE.upper")
  dat <- cbind(dat, UKBB.stats, gnomAD.stats, combined.stats, Coe_CASE.stats)
  freqmax.stats <- do.call("rbind", lapply(1:nrow(dat), function(i){
    freqmax.idx <- head(which(c(dat$UKBB.freq[i], dat$gnomAD.freq[i])==max(c(dat$UKBB.freq[i], dat$gnomAD.freq[i]))), 1)
    c(dat[i, which(colnames(dat)==paste(c("UKBB", "gnomAD")[freqmax.idx], "freq", sep="."))], 
      dat[i, which(colnames(dat)==paste(c("UKBB", "gnomAD")[freqmax.idx], "lower", sep="."))], 
      dat[i, which(colnames(dat)==paste(c("UKBB", "gnomAD")[freqmax.idx], "upper", sep="."))])
  }))
  colnames(freqmax.stats) <- c("freqmax.freq", "freqmax.lower", "freqmax.upper")
  dat <- cbind(dat, freqmax.stats)
  #Calculate binomial p-value that gnomAD is different from UKBB
  dat$binom.p <- sapply(1:nrow(dat), function(i){
    binom.test(x=dat$gnomAD.k[i], n=dat$gnomAD.N[i], p=dat$UKBB.freq[i])$p.value
  })
  #Calculate Fisher's exact test for cases vs controls
  get.fisher.stats <- function(dat, control.N, control.k, i){
    fisher <- fisher.test(matrix(c(control.N-control.k, control.k, 
                                   dat$Coe_CASE.N[i]-dat$Coe_CASE.k[i], dat$Coe_CASE.k[i]), 
                                 nrow=2, byrow=F))
    return(as.numeric(c(fisher$estimate, fisher$conf.int, fisher$p.value)))
  }
  #Calculate Fisher's exact test for Coe_CASE vs UKBB
  DD.stats.UKBB <- do.call("rbind", lapply(1:nrow(dat), function(i){
    get.fisher.stats(dat=dat, 
                     control.N=dat$UKBB.N[i], 
                     control.k=dat$UKBB.k[i], 
                     i=i)
  }))
  colnames(DD.stats.UKBB) <- c("DD.UKBB.OR", "DD.UKBB.lower", 
                               "DD.UKBB.upper", "DD.UKBB.fisher.p")
  #Calculate Fisher's exact test for Coe_CASE vs gnomAD
  DD.stats.gnomAD <- do.call("rbind", lapply(1:nrow(dat), function(i){
    get.fisher.stats(dat=dat, 
                     control.N=dat$gnomAD.N[i], 
                     control.k=dat$gnomAD.k[i], 
                     i=i)
  }))
  colnames(DD.stats.gnomAD) <- c("DD.gnomAD.OR", "DD.gnomAD.lower", 
                                 "DD.gnomAD.upper", "DD.gnomAD.fisher.p")
  #Calculate Fisher's exact test for Coe_CASE vs combined controls
  DD.stats.combined <- do.call("rbind", lapply(1:nrow(dat), function(i){
    combined.N <- dat$UKBB.N[i]+dat$gnomAD.N[i]
    combined.k <- dat$UKBB.k[i]+dat$gnomAD.k[i]
    get.fisher.stats(dat=dat, 
                     control.N=combined.N, 
                     control.k=combined.k, 
                     i=i)
  }))
  colnames(DD.stats.combined) <- c("DD.combined.OR", "DD.combined.lower", 
                                   "DD.combined.upper", "DD.combined.fisher.p")
  #Calculate Fisher's exact test for Coe_CASE vs max of controls
  DD.stats.freqmax <- do.call("rbind", lapply(1:nrow(dat), function(i){
    freq.max <- head(which(c(dat$UKBB.freq[i], dat$gnomAD.freq[i])==max(c(dat$UKBB.freq[i], dat$gnomAD.freq[i]))), 1)
    freqmax.control.N <- c(dat$UKBB.N[i], dat$gnomAD.N[i])[freq.max]
    freqmax.control.k <- c(dat$UKBB.k[i], dat$gnomAD.k[i])[freq.max]
    get.fisher.stats(dat=dat, 
                     control.N=freqmax.control.N, 
                     control.k=freqmax.control.k, 
                     i=i)
  }))
  colnames(DD.stats.freqmax) <- c("DD.freqmax.OR", "DD.freqmax.lower", 
                                  "DD.freqmax.upper", "DD.freqmax.fisher.p")
  #Calculate Fisher's exact test for Coe_CASE vs Coe_CTRL
  DD.stats.CoeCoe <- do.call("rbind", lapply(1:nrow(dat), function(i){
    get.fisher.stats(dat=dat, 
                     control.N=dat$Coe_CTRL.N[i], 
                     control.k=dat$Coe_CTRL.k[i], 
                     i=i)
  }))
  colnames(DD.stats.CoeCoe) <- c("DD.CoeCoe.OR", "DD.CoeCoe.lower", 
                                 "DD.CoeCoe.upper", "DD.CoeCoe.fisher.p")
  dat <- cbind(dat, DD.stats.UKBB, DD.stats.gnomAD, DD.stats.combined, DD.stats.freqmax, DD.stats.CoeCoe)
  return(dat)
}
#Scatterplot of UKBB vs gnomAD carrier freqs
gd.scatter.freq <- function(dat, svtypes, ax.lims=c(0, 0.0015), 
                            axes=T, highlight=NULL, pt.labels=F, ax.by=NULL,
                            plot.fit=F, plot.ci=T){
  #Prep plot values
  if(is.null(ax.lims)){
    ax.lims <- c(0, max(c(dat$UKBB.freq, dat$gnomAD.freq)))
  }
  pt.colors <- sapply(dat$CNV, function(CNV){
    svtypes$color[which(svtypes$svtype==CNV)]
  })
  #Prep plot area
  if(axes==T){
    par(mar=c(3.1, 3.1, 0.75, 0.75))
    pt.cex <- 1
  }else{
    par(mar=c(1.25, 1.25, 0.3, 0.3))
    pt.cex <- 1
  }
  ax.lims <- c(min(ax.lims), 1.1*max(ax.lims))
  plot(x=ax.lims, y=ax.lims, type="n", 
       xaxt="n", xlab="", yaxt="n", ylab="")
  # abline(0, 1, col="gray50", lty=2)
  if(!is.null(highlight)){
    rect(xleft=highlight[1], xright=highlight[2], 
         ybottom=highlight[1], ytop=highlight[2], 
         bty="n", border=NA, col=adjustcolor("yellow", alpha=0.3))
  }
  #Add axes
  if(axes==T){
    if(is.null(ax.by)){
      tick.at <- axTicks(1)
      ax.at <- axTicks(1)
    }else{
      tick.at <- seq(0, 2*max(ax.lims), by=ax.by/2)
      ax.at <- seq(0, 2*max(ax.lims), by=ax.by)
    }
    axis(1, at=tick.at, labels=NA, tck=-0.02)
    axis(1, at=ax.at, tick=F, labels=paste(round(100*ax.at, 3), "%", sep=""), 
         cex.axis=0.8, line=-0.9)
    mtext(1, line=1, text="UKBB CNV Freq.")
    axis(2, at=tick.at, labels=NA, tck=-0.02)
    axis(2, at=ax.at, tick=F, labels=paste(round(100*ax.at, 3), "%", sep=""), 
         cex.axis=0.8, line=-0.75, las=2)
    mtext(2, line=2.25, text="gnomAD CNV Freq.")
  }else{
    ax.at <- seq(0, 1, by=0.01)
    axis(1, at=ax.at, labels=NA, tck=-0.015)
    axis(1, at=ax.at, tick=F, labels=paste(round(100*ax.at, 2), "%", sep=""), 
         cex.axis=0.8, line=-1.1)
    axis(2, at=ax.at, labels=NA, tck=-0.015)
    axis(2, at=ax.at, tick=F, labels=paste(round(100*ax.at, 2), "%", sep=""), 
         cex.axis=0.8, line=-0.85, las=2)
  }
  #Add confidence intervals
  if(plot.ci==T){
    segments(x0=dat$UKBB.lower, x1=dat$UKBB.upper, 
             y0=dat$gnomAD.freq, y1=dat$gnomAD.freq, 
             col=adjustcolor(pt.colors, alpha=1/4), lwd=1)
    segments(x0=dat$UKBB.freq, x1=dat$UKBB.freq, 
             y0=dat$gnomAD.lower, y1=dat$gnomAD.upper, 
             col=adjustcolor(pt.colors, alpha=1/4), lwd=1)
  }
  #Add linear trend
  if(plot.fit==T){
    fit <- lm(gnomAD.freq ~ UKBB.freq, data=dat)
    abline(fit, lwd=2, col="gray50", lty=1)
  }
  #Add points
  points(x=dat$UKBB.freq, y=dat$gnomAD.freq, pch=21, cex=pt.cex, bg=pt.colors)
  if(pt.labels==T){
    text(x=dat$UKBB.freq, y=dat$gnomAD.freq, labels=dat$ID, cex=0.7, srt=30)
  }
  #Cleanup
  if(axes==T){
    box()
  }else{
    box(lwd=2)
  }
}
#Swarmplot of GD size for different vs not different freqs between UKBB and gnomAD
plot.sizeSwarm <- function(dat, svtypes){
  #Get plot data
  sizes.diff <- log10(dat$end-dat$start)[which(dat$binom.p<0.05/nrow(dat))]
  sizes.nodiff <- log10(dat$end-dat$start)[which(dat$binom.p>=0.05/nrow(dat))]
  yrange <- range(c(sizes.diff, sizes.nodiff), na.rm=T)
  yrange[2] <- yrange[2]+(0.2*(yrange[2]-yrange[1]))
  pt.colors <- sapply(dat$CNV, function(CNV){
    svtypes$color[which(svtypes$svtype==CNV)]
  })
  pt.colors.diff <- pt.colors[which(dat$binom.p<0.05/nrow(dat))]
  pt.colors.nodiff <- pt.colors[which(dat$binom.p>=0.05/nrow(dat))]
  #Prep plot area
  par(mar=c(1.75, 2.4, 0.25, 0.25), bty="n")
  plot(x=c(0, 2), y=yrange, type="n", 
       xaxt="n", yaxt="n", xlab="", ylab="")
  #Add swarms
  beeswarm(x=sizes.nodiff, add=T, at=0.5, 
           pwbg=pt.colors.nodiff, pch=21, cex=0.75, lwd=0.5, 
           corral="random", corralWidth=1)
  beeswarm(x=sizes.diff, add=T, at=1.5, 
           pwbg=pt.colors.diff, pch=21, cex=0.75, lwd=0.5, 
           corral="random", corralWidth=1)
  #Add medians
  segments(x0=(0:1)+0.2, x1=(1:2)-0.2, 
           y0=c(median(sizes.nodiff), median(sizes.diff)), 
           y1=c(median(sizes.nodiff), median(sizes.diff)), 
           lwd=3)
  #Add difference bar
  nodiff.max <- max(sizes.nodiff)
  diff.max <- max(sizes.diff)
  max.max <- max(c(diff.max, nodiff.max))
  spacer <- 0.05*(par("usr")[4]-par("usr")[3])
  segments(x0=c(0.5, 0.5, 1.5), x1=c(0.5, 1.5, 1.5), 
           y0=c(nodiff.max+spacer, 
                max.max+(1.5*spacer), 
                diff.max+spacer), 
           y1=rep(max.max+(1.5*spacer), 3))
  pval <- wilcox.test(sizes.nodiff, sizes.diff, alternative="greater")$p.value
  text(x=1, y=max.max+(2.5*spacer), 
       labels=bquote(italic(P) == .(round(pval, 3))))
  #Add axes
  axis(1, at=0.5, cex.axis=0.85, tick=F, line=-0.2, labels="Same\nFreqs")
  axis(1, at=1.5, cex.axis=0.85, tick=F, line=-0.2, labels="Diff.\nFreqs")
  axis(2, at=log10(as.numeric(unlist(sapply(1:10, function(x){(1:9)*(10^x)})))), 
       labels=NA, tck=-0.025)
  axis(2, at=1:10, tck=-0.05, labels=NA)
  axis(2, at=4:7, line=-0.7, tick=F, las=2, cex.axis=0.75, labels=c("10kb", "100kb", "1Mb", "10Mb"))
  mtext(2, line=1.5, text="GD Size")
}
#Swarmplot of GD freqs across populations
plot.freqByPop <- function(bypop, pops){
  #Get plot data
  plot.dat <- bypop[, -c(1:5)]
  pwcol <- rep(NA, times=nrow(bypop))
  pwcol[which(bypop$ID=="2q13dupNPHP1")] <- "black"
  yrange <- range(plot.dat)
  #Prep plot area
  par(mar=c(1.75, 2.5, 0.25, 0.25), bty="n")
  plot(x=c(0, nrow(pops)), y=yrange, type="n", 
       xaxt="n", yaxt="n", xlab="", ylab="")
  #Add swarms
  beeswarm(x=plot.dat, add=T, at=(1:nrow(pops))-0.5, 
           bg=pops$color, pwcol=rep(pwcol, nrow(pops)), 
           pch=21, cex=0.8, lwd=1.2, 
           corral="random", corralWidth=0.75)
  #Add axes
  axis(1, at=(1:nrow(pops))-0.5, cex.axis=0.85, tick=F,
       las=2, line=-1, labels=pops$pop)
  axis(2, at=c(0, 2*max(axTicks(2))), labels=NA, tck=0)
  axis(2, at=axTicks(2), labels=NA, tck=-0.025)
  axis(2, at=axTicks(2), line=-0.8, tick=F, las=2, cex.axis=0.7, 
       labels=paste(format(round(100*axTicks(2), 1), 1), "%", sep=""))
  mtext(2, line=1.6, text="CNV Frequency")
}
#Get scatterplot of control carrier freq vs DD odds ratio
gd.scatter.or_vs_freq <- function(dat, svtypes, comparison, comparison.freq=NULL, comparison.OR=NULL, 
                                  pt.labels=F, xlabel=NULL, title=NULL, cex.axlabs=1, 
                                  ylabel="Odds Ratio in Dev. Disorders", plot.ci=T, plot.fit=T){
  #Prep plot values
  if(is.null(comparison.freq)){
    comparison.freq <- comparison
  }
  if(is.null(comparison.OR)){
    comparison.OR <- comparison
  }
  freq <- log10(dat[, which(colnames(dat)==paste(comparison.freq, "freq", sep="."))])
  freq.lower <- log10(dat[, which(colnames(dat)==paste(comparison.freq, "lower", sep="."))])
  freq.upper <- log10(dat[, which(colnames(dat)==paste(comparison.freq, "upper", sep="."))])
  OR.mean <- log2(dat[, which(colnames(dat)==paste("DD", comparison.OR, "OR", sep="."))])
  OR.lower <- log2(dat[, which(colnames(dat)==paste("DD", comparison.OR, "lower", sep="."))])
  OR.upper <- log2(dat[, which(colnames(dat)==paste("DD", comparison.OR, "upper", sep="."))])
  #Fit lm to data prior to cleaning
  lm.fit <- lm(OR.mean[which(!is.infinite((OR.mean)) & freq<0 & !is.infinite(freq))] ~ freq[which(!is.infinite((OR.mean)) & freq<0 & !is.infinite(freq))])
  cor.fit <- cor.test(OR.mean[which(!is.infinite((OR.mean)) & freq<0 & !is.infinite(freq))], 
                      freq[which(!is.infinite((OR.mean)) & freq<0 & !is.infinite(freq))])
  print(cor.fit)
  #Clean plot values
  if(any(is.infinite(freq[which(freq<0)]))){
    some.infs.x <- T
    old.freq.max <- max(freq[which(!is.infinite(freq))])
    freq.range <- range(freq[which(!is.infinite(freq))])
    freq[which(is.infinite(freq))] <- max(freq.range)+(0.075*(freq.range[2]-freq.range[1]))
  }else{
    some.infs.x <- F
  }
  if(any(is.infinite(OR.mean[which(OR.mean>0)]))){
    some.infs.y <- T
    old.OR.max <- max(OR.mean[which(!is.infinite(OR.mean))])
    OR.range <- range(OR.mean[which(!is.infinite(OR.mean))])
    OR.mean[which(is.infinite(OR.mean))] <- max(OR.range)+(0.075*(OR.range[2]-OR.range[1]))
  }else{
    some.infs.y <- F
  }
  freq[which(is.infinite(freq) & freq<0)] <- min(freq[which(!is.infinite(freq))])
  OR.mean[which(is.infinite(OR.mean) & OR.mean<0)] <- min(OR.mean[which(!is.infinite(OR.mean))])
  freq.range <- range(freq)
  OR.range <- range(OR.mean)
  pt.colors <- sapply(dat$CNV, function(CNV){
    svtypes$color[which(svtypes$svtype==CNV)]
  })
  #Prep plot area
  if(is.null(title)){
    par(mar=c(3.25, 3.25, 1, 1))
  }else{
    par(mar=c(3.25, 3.25, 2.1, 2.1))
  }
  plot(x=freq.range, y=c(OR.range[1], OR.range[2]+1), type="n", 
       xaxt="n", xlab="", yaxt="n", ylab="")
  if(plot.fit==T){
    abline(lm.fit, lwd=2, col="gray50", lty=1)
  }
  #Add axes
  if(some.infs.x==T){
    x.at <- seq(floor(min(freq.range)), floor(max(old.freq.max)), 1)
    x.inf.at <- max(freq.range)
  }else{
    x.at <- seq(floor(min(freq.range)), ceiling(max(freq.range)), 1)
    x.inf.at <- NULL
  }
  if(some.infs.y==T){
    y.tick.at <- seq(floor(min(OR.range)), floor(max(old.OR.max)), 1)
    y.labels.at <- seq(floor(min(OR.range)), floor(max(old.OR.max)), 2)
    y.inf.at <- max(OR.range)
  }else{
    y.tick.at <- seq(floor(min(OR.range)), ceiling(max(OR.range)), 1)
    y.labels.at <- seq(floor(min(OR.range)), ceiling(max(OR.range)), 2)
    y.inf.at <- NULL
  }
  axis(1, at=-log10(as.numeric(sapply(x.at, function(i){(1:9)*10^-i}))), labels=NA, tck=-0.01)
  axis(1, at=x.at, labels=NA, tck=-0.02)
  sapply(x.at, function(x){
    axis(1, at=x, tick=F, labels=paste(round(100*(10^x), 5), "%", sep=""), 
         cex.axis=0.8, line=-0.9)
  })
  if(some.infs.x==T){
    axis(1, at=x.inf.at, labels=NA, tck=-0.02)
    axis(1, at=x.inf.at, tick=F, labels="0%", cex.axis=0.8, line=-0.9)
  }
  mtext(1, line=1, text=xlabel, cex=cex.axlabs)
  axis(2, at=y.tick.at, labels=NA, tck=-0.02)
  axis(2, at=y.labels.at, tick=F, labels=round(2^y.labels.at, 3), 
       cex.axis=0.8, line=-0.7, las=2)
  if(some.infs.y==T){
    axis(2, at=y.inf.at, labels=NA, tck=-0.02)
    axis(2, at=y.inf.at, tick=F, labels="Inf.", cex.axis=0.8, las=2, line=-0.75)
  }
  mtext(2, line=1.65, text=ylabel, cex=cex.axlabs)
  #Add 95% CIs
  if(plot.ci==T){
    segments(x0=freq.lower, x1=freq.upper, 
             y0=OR.mean, y1=OR.mean, 
             col=adjustcolor(pt.colors, alpha=1/4), lwd=1.5)
    segments(x0=freq, x1=freq, 
             y0=OR.lower, y1=OR.upper, 
             col=adjustcolor(pt.colors, alpha=1/4), lwd=1.5)
  }
  #Add axis breaks
  if(some.infs.x==T){
    rect(xleft=old.freq.max+(0.015*(old.freq.max-min(freq.range))), 
         xright=par("usr")[2],
         ybottom=par("usr")[3], 
         ytop=par("usr")[4],
         bty="n", border=NA, col="white")
    axis.break(1, breakpos=old.freq.max+(0.03*(old.freq.max-min(freq.range))))
    axis.break(3, breakpos=old.freq.max+(0.03*(old.freq.max-min(freq.range))))
    # rect(xleft=old.freq.max+(0.015*(old.freq.max-min(freq.range))), 
    #      xright=old.freq.max+(0.045*(old.freq.max-min(freq.range))),
    #      ybottom=par("usr")[3], ytop=par("usr")[4],
    #      bty="n", border=NA, col="gray90")
    # abline(v=c(old.freq.max+(0.015*(old.freq.max-min(freq.range))),
    #            old.freq.max+(0.045*(old.freq.max-min(freq.range)))),
    #        lty=2, col="gray80")
    # abline(v=old.freq.max+(0.03*(old.freq.max-min(freq.range))),
    #        lty=2, col="gray80")
  }
  if(some.infs.y==T){
    rect(xleft=par("usr")[1], 
         xright=par("usr")[2],
         ybottom=old.OR.max+(0.015*(old.OR.max-min(OR.range))),
         ytop=par("usr")[4],
         bty="n", border=NA, col="white")
    axis.break(2, breakpos=old.OR.max+(0.03*(old.OR.max-min(OR.range))))
    axis.break(4, breakpos=old.OR.max+(0.03*(old.OR.max-min(OR.range))))
    # rect(xleft=par("usr")[1], 
    #      xright=par("usr")[2],
    #      ybottom=old.OR.max+(0.015*(old.OR.max-min(OR.range))),
    #      ytop=par("usr")[4],
    #      bty="n", border=NA, col="gray90", density=15)
    # rect(xleft=par("usr")[1], xright=par("usr")[2],
    #      ybottom=old.OR.max+(0.015*(old.OR.max-min(OR.range))),
    #      ytop=old.OR.max+(0.045*(old.OR.max-min(OR.range))),
    #      bty="n", border=NA, col="gray90")
    # abline(h=c(old.OR.max+(0.015*(old.OR.max-min(OR.range))),
    #            old.OR.max+(0.045*(old.OR.max-min(OR.range)))),
    #        lty=2, col="gray80")
    # abline(h=old.OR.max+(0.03*(old.OR.max-min(OR.range))),
    #        lty=2, col="gray80")
  }
  mtext(3, line=1, text=title, font=2)
  r2 <- cor.fit$estimate^2
  cor.p <- format(cor.fit$p.value, scientific=T)
  cor.p.base <- format(as.numeric(strsplit(cor.p, split="e")[[1]][1]), nsmall=2)
  cor.p.exp <- format(as.numeric(strsplit(cor.p, split="e")[[1]][2]), nsmall=0)
  mtext(3, line=-0.1, cex=0.95*cex.axlabs, 
        text=bquote(R^2 == .(format(round(r2, 2), nsmall=2)) * ";" ~ italic(P) == .(format(round(as.numeric(cor.p.base), 2), nsmall=2))*"x"*10^.(cor.p.exp)))
  
  #Add points
  points(x=freq, y=OR.mean, pch=21, cex=cex.axlabs, bg=pt.colors)
  if(pt.labels==T){
    text(x=freq, y=OR.mean, labels=dat$ID, cex=0.7, srt=30)
  }
  #Cleanup
  box()
}

#Plot ordered ORs for all 54 GDs
plot.ordered.ors <- function(dat, svtypes, comparison, log2=F){
  #Prep plot values
  OR.mean <- dat[, which(colnames(dat)==paste("DD", comparison, "OR", sep="."))]
  OR.lower <- dat[, which(colnames(dat)==paste("DD", comparison, "lower", sep="."))]
  OR.upper <- dat[, which(colnames(dat)==paste("DD", comparison, "upper", sep="."))]
  if(log2==T){
    OR.mean <- log2(OR.mean)
    OR.lower <- log2(OR.lower)
    OR.upper <- log2(OR.upper)
  }
  pt.colors <- sapply(dat$CNV, function(CNV){
    svtypes$color[which(svtypes$svtype==CNV)]
  })
  pt.order <- order(-OR.mean)
  #Prep plot area
  par(mar=c(0.75, 2.25, 0.25, 2.75), bty="n")
  plot(x=c(0, length(OR.mean)), y=range(OR.mean[which(!is.infinite(OR.mean))]), type="n", 
       xaxt="n", yaxt="n", xlab="", ylab="")
  # if(log2==T){
  #   abline(h=log2(1), col="gray80")
  #   # abline(h=log2(10), lty=2)
  # }else{
  #   abline(h=1, col="gray80")
  #   # abline(h=10, lty=2)
  # }
  OR.quartiles <- quantile(OR.mean, probs=seq(0, 1, 0.25))
  OR.quartiles[which(is.infinite(OR.quartiles) & OR.quartiles<0)] <- par("usr")[3]
  carriersByQuartile <- sapply(1:4, function(i){
    sum(dat$combined.freq[which(OR.mean>OR.quartiles[i] & OR.mean<=OR.quartiles[i+1])])
  })
  segments(x0=1:length(OR.mean), x1=1:length(OR.mean), 
           y0=OR.lower[pt.order], y1=OR.upper[pt.order], 
           col=adjustcolor(pt.colors[pt.order], alpha=1/4), lwd=1.5)
  par(xpd=T)
  sapply(length(OR.mean):1, function(i){
    points(x=i, y=OR.mean[pt.order[i]], pch=21, bg=pt.colors[pt.order[i]], lwd=0.75, cex=0.9)
  })
  par(xpd=F)
  abline(h=OR.quartiles[c(2, 4:5)], lty=2)
  axis(1, at=c(-100, 100), labels=NA)
  mtext(1, text="GDs Ordered by Revised OR", cex=0.75, line=-0.2)
  if(log2==T){
    axis(2, at=c(-100, 0:10), labels=NA, tck=-0.02)
    axis(2, at=seq(0, 10, 2), labels=2^seq(0, 10, 2), tick=F, las=2, cex.axis=0.8, line=-0.8)
  }else{
    axis(2, at=axTicks(2)[-1], tick=F, las=2, cex.axis=0.8, line=-0.8)
    axis(2, at=c(-100, axTicks(2)[-1]), labels=NA, tck=-0.02)
  }
  mtext(2, text="Revised OR", line=1)
  axis(4, at=OR.quartiles[1:2], tck=0.075, labels=NA)
  axis(4, at=OR.quartiles[4:5], tck=0.075, labels=NA)
  axis(4, at=OR.quartiles[1:2], tck=-0.075, labels=NA)
  axis(4, at=OR.quartiles[4:5], tck=-0.075, labels=NA)
  axis(4, at=(OR.quartiles[c(1, 4)]+OR.quartiles[c(2, 5)])/2, tick=F, line=-0.9, 
       labels=paste(round(100*carriersByQuartile[c(1, 4)], 2), "%", sep=""), las=2, cex=0.8)
  print(dat$ID[pt.order])
}


################
###RSCRIPT BLOCK
################
require(optparse, quietly=T)
require(plotrix, quietly=T)
require(beeswarm, quietly=T)
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
args <- parse_args(OptionParser(usage="%prog DATA BYPOP OUTDIR", 
                                option_list=option_list), 
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
data.in <- args$args[1]
data_bypop.in <- args$args[2]
OUTDIR <- args$args[3]
prefix <- opts$prefix
svtypes.file <- opts$svtypes
pops.file <- opts$populations

# #Dev parameters (local)
# data.in <- "~/scratch/gnomAD-SV_v2_rev1.UKBB_GD_comparison_results.txt"
# data_bypop.in <- "~/scratch/gnomAD-SV_v2_rev1.GD_carrier_rates_byPop.bed"
# OUTDIR <- "~/scratch/gnomAD_v2_test_plots/"
# svtypes.file <- "~/Desktop/Collins/Talkowski/code/sv-pipeline/ref/vcf_qc_refs/SV_colors.txt"
# prefix <- "gnomAD-SV_v2_rev1"
# pops.file <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/gnomAD_population_colors.txt"

###Create output directory, if necessary
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

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


###Load & process data
dat <- import.data(data.in)
bypop <- read.table(data_bypop.in, header=T, sep="\t", comment.char="")
bypop <- bypop[, which(colnames(bypop) != "OTH")]

###Print list of GDs above Bonferroni-corrected significance
print(dat[which(dat$binom.p<0.05/nrow(dat)), ])

###Print correlation statistics between gnomAD and UKBB
cor(dat$UKBB.freq, dat$gnomAD.freq)^2
format(cor.test(dat$UKBB.freq, dat$gnomAD.freq)$p.value, scientific=T)

# ###Swarmplot of GD sizes where UKBB and gnomAD do and do not differ
# pdf(paste(OUTDIR, "/", prefix, ".UKBB_vs_gnomAD.size_swarm.pdf", sep=""), 
#     height=2.25, width=1.5)
# plot.sizeSwarm(dat, svtypes)
# dev.off()

###Swarmplot of GD freqs by pop
pdf(paste(OUTDIR, "/", prefix, ".GD_freq_by_pop.swarm.pdf", sep=""), 
    height=2.25, width=1.75)
plot.freqByPop(bypop, pops[which(pops$pop != "OTH"), ])
dev.off()

###Plot large scatter for zoomed panel
pdf(paste(OUTDIR, "/", prefix, ".UKBB_GD_freq_comparison.zoomed.pdf", sep=""), 
    height=2.75, width=2.75)
gd.scatter.freq(dat, svtypes, ax.lims=c(0, 0.001), ax.by=0.0005, plot.fit=T, plot.ci=T)
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".UKBB_GD_freq_comparison.zoomed_smaller.pdf", sep=""), 
    height=2.25, width=2.25)
gd.scatter.freq(dat, svtypes, ax.lims=c(0, 0.001), ax.by=0.0005, plot.fit=T, plot.ci=T)
dev.off()


###Plot smaller scatter for zoomed out panel
pdf(paste(OUTDIR, "/", prefix, ".UKBB_GD_freq_comparison.full.pdf", sep=""), 
    height=1.5, width=1.5)
gd.scatter.freq(dat, svtypes, ax.lims=NULL, axes=F, plot.fit=T, plot.ci=T)
dev.off()

###Full-size plot with all points labeled -- for viewing only
pdf(paste(OUTDIR, "/", prefix, ".UKBB_GD_freq_comparison.full.with_labels.pdf", sep=""), 
    height=5, width=5)
gd.scatter.freq(dat, svtypes, ax.lims=NULL, axes=T, pt.labels=T, plot.ci=F)
dev.off()

###Zoomed-in plot with all points labeled -- for viewing only
pdf(paste(OUTDIR, "/", prefix, ".UKBB_GD_freq_comparison.zoomed.with_labels.pdf", sep=""), 
    height=5, width=5)
gd.scatter.freq(dat, svtypes, axes=T, pt.labels=T, ax.lims=c(0, 0.001), plot.ci=F)
dev.off()

###Zoomed-in plot with all points labeled -- for viewing only
pdf(paste(OUTDIR, "/", prefix, ".UKBB_GD_freq_comparison.super_zoomed.with_labels.pdf", sep=""), 
    height=5, width=5)
gd.scatter.freq(dat, svtypes, ax.lims=c(0, 0.0008), axes=T, pt.labels=T, plot.ci=F)
dev.off()

###Four-panel grid of freq vs OR in DDs
pdf(paste(OUTDIR, "/", prefix, ".UKBB_GD_DD_OR_vs_freq.fourpanel.pdf", sep=""), 
    height=6, width=6)
par(mfrow=c(2, 2))
gd.scatter.or_vs_freq(dat=dat, svtypes=svtypes, 
                      comparison="UKBB", pt.labels=F, 
                      comparison.OR="CoeCoe", 
                      xlabel="UKBB CNV Freq.", 
                      ylabel="CNV Odds Ratio (DDs)", 
                      title="UKBB Only")
gd.scatter.or_vs_freq(dat=dat, svtypes=svtypes, 
                      comparison="gnomAD", pt.labels=F, 
                      comparison.OR="CoeCoe", 
                      xlabel="gnomAD CNV Freq.", 
                      ylabel="CNV Odds Ratio (DDs)", 
                      title="gnomAD Only")
gd.scatter.or_vs_freq(dat=dat, svtypes=svtypes, 
                      comparison="combined", pt.labels=F, 
                      comparison.OR="CoeCoe", 
                      xlabel="UKBB + gnomAD CNV Freq.", 
                      ylabel="CNV Odds Ratio (DDs)", 
                      title="UKBB + gnomAD")
gd.scatter.or_vs_freq(dat=dat, svtypes=svtypes, 
                      comparison="freqmax", pt.labels=F, 
                      comparison.OR="CoeCoe", 
                      xlabel="Max CNV Freq.", 
                      ylabel="CNV Odds Ratio (DDs)", 
                      title="Max Freq. (UKBB or gnomAD)")
dev.off()

###Single plot of freq vs OR in DDs
pdf(paste(OUTDIR, "/", prefix, ".UKBB_GD_DD_OR_vs_freq.UKBB_plus_gnomAD.pdf", sep=""), 
    height=2.4, width=2.4)
gd.scatter.or_vs_freq(dat=dat, svtypes=svtypes, 
                      comparison="combined", pt.labels=F, 
                      comparison.OR="CoeCoe", 
                      xlabel="gnomAD+UKBB CNV Freq.", 
                      cex.axlabs=0.9, 
                      ylabel="Reported Odds Ratio (OR)",
                      plot.ci=T)
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".UKBB_GD_DD_OR_vs_freq.UKBB_gnomAD_freqmax.pdf", sep=""), 
    height=2.4, width=2.4)
gd.scatter.or_vs_freq(dat=dat, svtypes=svtypes, 
                      comparison="freqmax", pt.labels=F, 
                      comparison.OR="CoeCoe", 
                      xlabel="Max Freq. (gnomAD, UKBB)", 
                      cex.axlabs=0.9, 
                      ylabel="Reported Odds Ratio (OR)",
                      plot.ci=T)
dev.off()

###Plot ordered GD ORs
pdf(paste(OUTDIR, "/", prefix, ".UKBB_GD_ordered_DD_ORs.UKBB_plus_gnomAD.pdf", sep=""), 
    height=1.25, width=2.4)
plot.ordered.ors(dat=dat, svtypes=svtypes, comparison="combined", log2=T)
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".UKBB_GD_ordered_DD_ORs.UKBB_plus_gnomAD.long.pdf", sep=""), 
    height=1.25, width=3.4)
plot.ordered.ors(dat=dat, svtypes=svtypes, comparison="combined", log2=T)
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".UKBB_GD_ordered_DD_ORs.UKBB_gnomAD_freqmax.pdf", sep=""), 
    height=1.25, width=2.4)
plot.ordered.ors(dat=dat, svtypes=svtypes, comparison="freqmax", log2=T)
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".UKBB_GD_ordered_DD_ORs.UKBB_gnomAD_freqmax.long.pdf", sep=""), 
    height=1.25, width=3.4)
plot.ordered.ors(dat=dat, svtypes=svtypes, comparison="freqmax", log2=T)
dev.off()


###Write data out
write.table(dat, paste(OUTDIR, "/", prefix, ".UKBB_GD_comparison_data.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)
