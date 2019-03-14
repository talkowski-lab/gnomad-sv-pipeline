#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Plot chromosome-wide SV density for formal gnomAD analysis


###Set master parameters
options(stringsAsFactors=F,scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Read density pileups
loadPileup <- function(dat.in,centromeres.in,window=11,smooth=T){
  #Read data
  dat <- read.table(dat.in,header=T,comment.char="")
  colnames(dat)[1] <- "chr"
  cents <- read.table(centromeres.in,header=T,comment.char="")
  colnames(cents)[1] <- "chr"
  #Exclude bins overlapping centromeres
  drop.rows <- unlist(sapply(unique(dat$chr),function(k){
    cent.pos <- as.numeric(cents[which(cents$chr==k),2:3])
    which(dat$chr==k & dat$start<cent.pos[2] & dat$end>cent.pos[1])
  }))
  if(length(drop.rows)>0){
    dat <- dat[-drop.rows,]
  }
  #Sum counts
  dat$ALL <- apply(dat[,-c(1:3)],1,sum,na.rm=T)
  #Smooth counts per chromosome
  require(zoo,quietly=T)
  if(smooth==T){
    for(k in unique(dat$chr)){
      dat[which(dat$chr==k),-c(1:3)] <- apply(dat[which(dat$chr==k),-c(1:3)],2,function(vals){
        rollapply(as.numeric(vals),width=window,FUN=mean,partial=T)
      })
    }
  }
  #Assign normalized distance from centromere
  #Range: -1 = at pter, 0 = in centromere, +1 = at qter
  dat$cdist.norm <- as.numeric(unlist(sapply(unique(dat$chr),function(k){
    midpoints <- (dat$end[which(dat$chr==k)]+dat$end[which(dat$chr==k)])/2
    c.mid <- mean(as.numeric(cents[which(cents$chr==k),2:3]))
    c.dists <- midpoints-c.mid
    plen <- c.mid
    qlen <- max(midpoints)-c.mid
    c.dists.n <- c.dists
    c.dists.n[which(c.dists<0)] <- c.dists.n[which(c.dists<0)]/plen
    c.dists.n[which(c.dists>0)] <- c.dists.n[which(c.dists>0)]/qlen
    c.dists.n[which(c.dists.n< -1)] <- -1
    c.dists.n[which(c.dists.n>1)] <- 1
    return(c.dists.n)
  })))
  return(dat)
}
#Gather meta-chromosome average for a single SVTYPE
metaAverage <- function(dat,SVTYPE="ALL",n.bins=250){
  p.bins <- seq(-1,0,by=1/n.bins)
  q.bins <- seq(0,1,by=1/n.bins)
  col.idx <- which(colnames(dat)==SVTYPE)
  p.means <- sapply(1:(length(p.bins)-1),function(i){
    mean(dat[which(dat$cdist.norm>=p.bins[i] & dat$cdist.norm<p.bins[i+1]),col.idx],na.rm=T)
  })
  q.means <- sapply(1:(length(q.bins)-1),function(i){
    mean(dat[which(dat$cdist.norm>q.bins[i] & dat$cdist.norm<=q.bins[i+1]),col.idx],na.rm=T)
  })
  means <- c(p.means,q.means)
  means.norm <- means/mean(means,na.rm=T)
  out.df <- data.frame("norm.pos"=c(p.bins[-length(p.bins)],q.bins[-1]),
                       "mean"=means,"mean.norm"=means.norm)
  return(out.df)
}
#Split densities into terminal, interstitial, and pericentromeric bins
calc.meanByContext <- function(dat,meta.svtypes,ter.buf=0.05,cen.buf=0.05){
  #Helper function to calculate mean, 95% CI, and p-value that the true mean isn't 1
  get.ci <- function(vals){
    vals <- vals[which(!is.na(vals))]
    k <- 1.96*(sd(vals,na.rm=T)/sqrt(length(vals)))
    p.less <- t.test(vals,mu=1,alternative="less")$p.value
    p.greater <- t.test(vals,mu=1,alternative="greater")$p.value
    return(c(log2(c(mean(vals)-k,mean(vals),mean(vals)+k)),p.less,p.greater))
  }
  res <- lapply(meta.svtypes,function(svtype){
    #Mean-normalize all values
    vals <- as.numeric(dat[,which(colnames(dat)==svtype)])
    vals <- vals/mean(vals,na.rm=T)
    #Calculate stats
    ter.idx <- which(dat$cdist.norm<=-1+ter.buf | dat$cdist.norm>=1-ter.buf)
    cen.idx <- which(dat$cdist.norm>=-cen.buf & dat$cdist.norm<=ter.buf)
    int.idx <- which(!(1:nrow(dat) %in% c(ter.idx,cen.idx)))
    ter.stats <- get.ci(vals[ter.idx])
    int.stats <- get.ci(vals[int.idx])
    cen.stats <- get.ci(vals[cen.idx])
    return(data.frame("ter"=ter.stats,
                      "int"=int.stats,
                      "cen"=cen.stats))
  })
  names(res) <- meta.svtypes
  return(res)
}


#####################
###PLOTTING FUNCTIONS
#####################
#Master function to plot raw coverage from all four samples per contig
generateCovPlotsPerChrom <- function(mat,
                                     colors,
                                     contigs.top,
                                     contigs.middle,
                                     contigs.bottom,
                                     labels.on.top=F,
                                     fill=T,norm=F){
  #Normalize data, if optioned
  if(norm==T){
    mat <- lapply(mat,function(chr.mat){
      chr.mat[,-c(1:3)] <- apply(as.data.frame(chr.mat[,-c(1:3)]),2,function(vals){
        vals <- vals/mean(vals,na.rm=T)
        return(vals)
      })
      return(chr.mat)
    })
  }
  
  #Set spacer
  spacer <- 30000000
  
  #Determine total length to be plotted
  chr.lengths <- as.numeric(unlist(lapply(mat,function(chr.mat){
    return(max(chr.mat[,3]))
  })))
  #Add spacer distance between contigs
  chr.lengths[1:22] <- chr.lengths[1:22]+spacer
  
  #Mini helper function to plot a row of summary coverage values
  plotMiniCovSummary <- function(contigs,ymax=NULL){
    #Prep plot area
    if(labels.on.top==T){
      par(mar=c(0.2,2,1.5,1),bty="n",bg="white")
      lpos=3
    }else{
      par(mar=c(1.5,2,0.2,1),bty="n",bg="white")
      lpos=1
    }
    if(is.null(ymax)){
      ymax <- as.integer(ceiling(quantile(unlist(lapply(mat,function(chr.mat){chr.mat[,4]})),probs=0.995)))
    }
    plot(x=c(-0.01*sum(chr.lengths[contigs]),sum(chr.lengths[contigs])-spacer),
         y=c(-0.1*ymax,1.1*ymax),type="n",
         xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i")
    #Add background shading rectangle
    rect(xleft=par("usr")[1],xright=par("usr")[2],
         ybottom=par("usr")[3],ytop=par("usr")[4],
         border="gray99",col="gray99")
    box(which="plot",col="white",lwd=3)
    #Add contig positions & labels
    sapply(1:length(contigs),function(i){
      axis(lpos,at=c(c(0,cumsum(chr.lengths[contigs]))[i],
                     c(0,cumsum(chr.lengths[contigs])-spacer)[i+1]),
           labels=NA,col="gray30",tck=0,line=0.1)
      axis(lpos,at=mean(c(c(0,cumsum(chr.lengths[contigs]))[i],
                          c(0,cumsum(chr.lengths[contigs])-spacer)[i+1])),
           tick=F,labels=paste("chr",contigs[i],sep=""),line=-0.5)
    })
    #Add y axis & gridlines
    y.at <- axTicks(2)[seq(1,length(axTicks(2)),by=2)]
    # y.at <- c(y.at,max(y.at)+(y.at[2]-y.at[1]))
    axis(2,at=axTicks(2),labels=NA,tck=-0.05,col="gray40")
    axis(2,at=y.at,labels=NA,tck=-0.1)
    axis(2,at=y.at,tick=F,las=2,cex.axis=1,line=-0.4)
    # mtext(2,line=2.25,text="Copy Number")
    abline(h=axTicks(2),col="gray80",lwd=0.7,lty=3)
    # abline(h=y.at,col="gray80")
    abline(h=0)
    #Add coverage values
    sapply(1:length(contigs),function(i){
      sapply(1:(ncol(mat[[1]])-3),function(s){
        # cov.vals <- rollmean(cov[[s]][[contigs[i]]][,4],k=11,na.pad=T)
        vals <- mat[[contigs[i]]][,s+3]
        plot.vals <- smooth.spline(x=mat[[contigs[i]]][,2]+c(0,cumsum(chr.lengths[contigs]))[i],
                                   y=vals,spar=0.32*mean(c(1,rep(chr.lengths[contigs[1]]/chr.lengths[contigs[i]],4))))
        if(fill==T){
          polygon(x=c(plot.vals$x,rev(plot.vals$x)),
                  y=c(plot.vals$y,rep(0,times=length(plot.vals$y))),
                  border=NA,col=adjustcolor(colors[s],alpha=0.7))
        }
        points(plot.vals$x,plot.vals$y,type="l",lwd=1.25,col=colors[s])
      })
    })
    #Add cleanup rectangles
    rect(xleft=par("usr")[1],xright=par("usr")[2],
         ybottom=c(par("usr")[3],ymax),ytop=c(0,par("usr")[4]),
         border="white",lty=1,lwd=1,col="white")
    abline(h=par("usr")[3:4],lwd=1,col="white")
    # abline(h=c(0,ymax),col="gray80")
    box(lwd=2,col="white")
    rect(xleft=(cumsum(chr.lengths[contigs])-spacer),
         xright=cumsum(chr.lengths[contigs]),
         ybottom=par("usr")[3],ytop=par("usr")[4],
         border="white",lty=0,lwd=0,col="white")
    rect(xleft=c(par("usr")[1],tail((cumsum(chr.lengths[contigs])-spacer),1)),
         xright=c(0,par("usr")[2]),
         ybottom=par("usr")[3],ytop=par("usr")[4],
         border="white",lty=0,lwd=0,col="white")
    #Add N-mask rectangles
    sapply(1:length(contigs),function(i){
      rect(xleft=Nmasks[which(Nmasks[,1]==contigs[i]),2]+c(0,cumsum(chr.lengths[contigs]))[i],
           xright=Nmasks[which(Nmasks[,1]==contigs[i]),3]+c(0,cumsum(chr.lengths[contigs]))[i],
           ybottom=0,ytop=ymax,
           border=NA,col="gray80")
      sapply(1:nrow(Nmasks[which(Nmasks[,1]==contigs[i]),]),function(k){
        if((Nmasks[which(Nmasks[,1]==contigs[i]),3]-Nmasks[which(Nmasks[,1]==contigs[i]),2])[k]>10000000){
          text(x=mean(c(Nmasks[which(Nmasks[,1]==contigs[i]),2][k]+c(0,cumsum(chr.lengths[contigs]))[i],
                        Nmasks[which(Nmasks[,1]==contigs[i]),3][k]+c(0,cumsum(chr.lengths[contigs]))[i])),
               y=ymax/2,labels="N",col="gray72",font=3,cex=1.2)
        }
      })
    })
  }
  
  #####Plot stacked panels
  par(mfrow=c(3,1))
  plotMiniCovSummary(contigs.top)
  plotMiniCovSummary(contigs.middle)
  plotMiniCovSummary(contigs.bottom)
}
###Plot meta-chromosome density for a single SVTYPE
plot.metaDist <- function(meta.dat,color,fill=T,norm=F,xlabel=NULL){
  #Clean meta dat
  meta.dat <- meta.dat[which(!is.na(meta.dat$mean)),]
  if(norm==T){
    meta.dat$mean <- meta.dat$mean.norm
  }
  plot.vals <- rollapply(data=meta.dat$mean.norm,width=6,mean,partial=T)
  #Set parameters
  ymax <- 1.075*max(max(plot.vals,na.rm=T),
              2*mean(plot.vals,na.rm=T))
  lpos <- 1
  #Prep plot area
  par(bty="n",bg="white")
  plot(x=1.015*range(meta.dat$norm.pos),y=c(-0.05*ymax,ymax),type="n",
       xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i")
  #Add background shading rectangle
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=par("usr")[3],ytop=par("usr")[4],
       border="gray99",col="gray99")
  box(which="plot",col="white",lwd=3)
  #Add contig positions & labels
  axis(lpos,at=par("usr")[1:2],
       labels=NA,col="gray30",tck=0,line=0.1)
  if(is.null(xlabel)){
    axis(lpos,at=mean(par("usr")[1:2]),
         tick=F,labels="Meta-chromosome",line=-0.6)
  }else{
    axis(lpos,at=mean(par("usr")[1:2]),
         tick=F,labels=xlabel,line=-0.6,cex.axis=0.8)
  }
  #Add y axis & gridlines
  y.at <- axTicks(2)
  y.at <- c(y.at,max(y.at)+(y.at[2]-y.at[1]))
  axis(2,at=y.at,labels=NA,tck=-0.05,col="gray40")
  # axis(2,at=y.at,labels=NA,tck=-0.1)
  axis(2,at=y.at,tick=F,las=2,cex.axis=1,line=-0.4,labels=round(y.at,2))
  # mtext(2,line=2.25,text="Copy Number")
  abline(h=y.at,col="gray80",lwd=0.7,lty=3)
  # abline(h=y.at,col="gray80")
  # abline(v=0,col="#963231")
  if(norm==T){
    abline(h=1,col="gray50")
  }
  #Add coverage values
  # plot.vals <- smooth.spline(x=meta.dat$norm.pos,y=meta.dat$mean,spar=0.2)
  # plot.vals <- data.frame("x"=meta.dat$norm.pos,
                          # "y"=meta.dat$mean)
  if(fill==T){
    polygon(x=c(meta.dat$norm.pos,rev(meta.dat$norm.pos)),
            y=c(plot.vals,rep(0,times=length(plot.vals))),
            border=NA,col=adjustcolor(color,alpha=0.7))
  }
  points(meta.dat$norm.pos,plot.vals,type="l",lwd=1.25,col=color)
  #Cleanup
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=par("usr")[3],ytop=0,
       bty="n",border=NA,col="white")
  abline(h=0)
}
###Plot point estimates and CIs for ter/int/cen averages by class
plot.metaByContext <- function(dat,meta.svtypes,colors){
  #Get point estimates and CIs
  plot.vals <- calc.meanByContext(dat,meta.svtypes)
  #Prep plot area
  ylims <- c(1.1*min(as.numeric(unlist(lapply(plot.vals,range,na.rm=T)))),
             1.1*max(as.numeric(unlist(lapply(plot.vals,range,na.rm=T)))))
  par(bty="n",mar=c(0.25,2.8,1.25,0.25))
  plot(x=c(0,length(plot.vals)-0.4),y=ylims,type="n",
       xaxt="n",yaxt="n",xlab="",ylab="")
  abline(h=0,col="gray50")
  #Add points per svtype
  sapply(1:length(plot.vals),function(i){
    segments(x0=i-c(1,0.8,0.6),
             x1=i-c(1,0.8,0.6),
             y0=as.numeric(plot.vals[[i]][1,]),
             y1=as.numeric(plot.vals[[i]][3,]),
             col=colors[i],lwd=2)
    points(x=i-c(1,0.8,0.6),
           y=as.numeric(plot.vals[[i]][2,]),
           pch=19,col=colors[i],cex=1.25)
    text(x=i-c(1,0.8,0.6),
         y=as.numeric(plot.vals[[i]][2,]),
         labels=c("T","I","C"),cex=0.6,font=2,col="white")
    #Add category label
    axis(3,at=i-c(1,0.6),tck=0,labels=NA,line=0.1)
    axis(3,at=i-0.8,tick=F,line=-0.9,cex.axis=0.75,labels=meta.svtypes[i],col.axis=colors[i])
    #Add p-values
    par(xpd=T)
    sapply(1:3,function(k){
      if(plot.vals[[i]][4,k]<0.05/(3*length(plot.vals))){
        text(x=(i-c(1,0.8,0.6))[k],
             y=plot.vals[[i]][2,k],
             pos=1,labels="*")
      }
      if(plot.vals[[i]][5,k]<0.05/(3*length(plot.vals))){
        text(x=(i-c(1,0.8,0.6))[k],
             y=plot.vals[[i]][2,k]-(0.04*(par("usr")[4]-par("usr")[3])),
             pos=3,labels="*")
      }
    })
    par(xpd=F)
  })

  #Clean up
  axis(2,at=axTicks(2),labels=NA,tck=-0.03)
  sapply(axTicks(2),function(y){
    axis(2,at=y,tick=F,las=2,cex.axis=0.9,line=-0.5,
         labels=bquote("2"^.(y)))
  })
  mtext(2,line=1.75,text="SV Fold-Enrichment")
}


################
###RSCRIPT BLOCK
################
require(optparse,quietly=T)
require(zoo,quietly=T)
###List of command-line options
option_list <- list(
  make_option(c("--prefix"), type="character", default="gnomAD_v2_SV",
              help="prefix used for naming outfiles [default %default]",
              metavar="character"),
  make_option(c("-S", "--svtypes"), type="character", default=NULL,
              help="tab-delimited file specifying SV types and HEX colors [default %default]",
              metavar="character")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog PILEUP CENTROMERES NMASK OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
dat.in <- args$args[1]
centromeres.in <- args$args[2]
nmask.in <- args$args[3]
OUTDIR <- args$args[4]
prefix <- opts$prefix
svtypes.file <- opts$svtypes

# #Dev parameters (local)
# dat.in <- "~/scratch/gnomAD_v2_SV_MASTER.binned_counts_bySVTYPE.autosomes.bed.gz"
# centromeres.in <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/h37.centromeres.bed.gz"
# nmask.in <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/h37.large_Nmasked_blocks.bed"
# OUTDIR <- "~/scratch/gnomAD_v2_test_plots/"
# prefix <- "gnomAD_v2_SV_MASTER"
# svtypes.file <- "~/Desktop/Collins/Talkowski/code/sv-pipeline/ref/vcf_qc_refs/SV_colors.txt"

###Create output directory, if necessary
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

###Process input data
cat("NOW LOADING DATA\n")
dat.raw <- loadPileup(dat.in,centromeres.in,smooth=F)
dat <- loadPileup(dat.in,centromeres.in)

###Sets sv types & colors
if(!is.null(svtypes.file)){
  svtypes <- read.table(svtypes.file,sep="\t",header=F,comment.char="",check.names=F)
  svtypes <- as.data.frame(apply(svtypes,2,as.character))
  colnames(svtypes) <- c("svtype","color")
}else{
  require(RColorBrewer,quietly=T)
  svtypes.v <- unique(dat$SVTYPE)
  svtypes.c <- brewer.pal(length(svtypes.v),"Dark2")
  svtypes <- data.frame("svtype"=svtypes.v,
                        "color"=svtypes.c)
}

###Set analysis parameters & read helper files
meta.svtypes <- c("ALL","DEL","DUP","MCNV","INS","INV","CPX")
Nmasks <- read.table(nmask.in,header=F)
names(Nmasks) <- c("chr","start","end")

###Gather meta-chromosome averages per percentile
meta.means <- lapply(meta.svtypes,function(svtype){
  m <- metaAverage(dat=dat,SVTYPE=svtype,n.bins=500)
})
names(meta.means) <- meta.svtypes

###Generate genome-wide panel of idiograms - ALL SV
png(paste(OUTDIR,"/",prefix,".SV_density_idiograms.png",sep=""),
    height=800,width=2600,res=300,bg="white")
generateCovPlotsPerChrom(mat=lapply(1:22,function(k){return(as.data.frame(dat[which(dat$chr==k),c(1:3,which(colnames(dat)=="ALL"))]))}),
                         colors="gray10",
                         contigs.top=1:5,
                         contigs.middle=6:12,
                         contigs.bottom=13:22,
                         fill=T)
dev.off()
png(paste(OUTDIR,"/",prefix,".SV_density_idiograms.shorter.png",sep=""),
    height=600,width=(2/3)*2600,res=300,bg="white")
generateCovPlotsPerChrom(mat=lapply(1:22,function(k){return(as.data.frame(dat[which(dat$chr==k),c(1:3,which(colnames(dat)=="ALL"))]))}),
                         colors="gray10",
                         contigs.top=1:5,
                         contigs.middle=6:12,
                         contigs.bottom=13:22,
                         fill=T)
dev.off()

###Generate "meta-chromosome" plots per SVTYPE
#Large plot: all svtypes
png(paste(OUTDIR,"/",prefix,".SV_density_metaIdiograms.All_SV.png",sep=""),
    height=750,width=2.5*(2600/7),res=300,bg="white")
par(mar=c(1.5,2,0.2,1))
plot.metaDist(meta.dat=meta.means[[which(names(meta.means)=="ALL")]],
              color="gray10",norm=T,
              xlabel="Meta-chromosome (Mean of Autosomes)")
dev.off()
#Smaller panels: individual svtypes
png(paste(OUTDIR,"/",prefix,".SV_density_metaIdiograms.bySVTYPE.png",sep=""),
    height=650,width=3*(2600/7),res=300,bg="white")
par(mfrow=c(2,length(meta.svtypes[-1])/2),
    mar=c(1.5,2,1.5,1))
lapply(meta.svtypes[-1],function(svtype){
  svt.col <- svtypes$color[which(svtypes$svtype==svtype)]
  plot.metaDist(meta.dat=meta.means[[which(names(meta.means)==svtype)]],
                color=svt.col,norm=T)
})
dev.off()

###Plot meta means by ter/int/cen
pdf(paste(OUTDIR,"/",prefix,".SV_density_meta.byChromContext.dotplots.pdf",sep=""),
    height=1.9,width=3.25)
plot.metaByContext(dat,meta.svtypes,
                   colors=c("gray10",
            svtypes$color[which(svtypes$svtype=="DEL")],
            svtypes$color[which(svtypes$svtype=="DUP")],
            svtypes$color[which(svtypes$svtype=="MCNV")],
            svtypes$color[which(svtypes$svtype=="INS")],
            svtypes$color[which(svtypes$svtype=="INV")],
            svtypes$color[which(svtypes$svtype=="CPX")]))
dev.off()
###Generate genome-wide panel of idiograms - ALL SV
png(paste(OUTDIR,"/",prefix,".SV_density_idiograms_bySVTYPE.png",sep=""),
    height=800,width=2600,res=300,bg="white")
generateCovPlotsPerChrom(mat=lapply(1:22,function(k){return(as.data.frame(dat[which(dat$chr==k),c(1:3,which(colnames(dat) %in% c("DEL","DUP","INS")))]))}),
                         colors=c(svtypes$color[which(svtypes$svtype=="DEL")],
                                  svtypes$color[which(svtypes$svtype=="DUP")],
                                  svtypes$color[which(svtypes$svtype=="INS")]),
                         contigs.top=1:5,
                         contigs.middle=6:12,
                         contigs.bottom=13:22,
                         fill=F,norm=T)
dev.off()

###Plot grid of boxplots by SVTYPE
# par(mfrow=c(2,2),bty="n")
# sapply(c("DEL","DUP","INS","INV+CPX"),function(svtype){
#   if(svtype=="INV+CPX"){
#     plot.dat <- data.frame("chr"=dat.raw$chr,
#                            "sv"=dat.raw$INV+dat.raw$CPX)
#     color <- svtypes$color[which(svtypes$svtype=="CPX")]
#   }else{
#     plot.dat <- data.frame("chr"=dat.raw$chr,
#                            "sv"=dat.raw[,which(colnames(dat)==svtype)])
#     color <- svtypes$color[which(svtypes$svtype==svtype)]
#   }
#   boxplot(sv ~ chr,data=plot.dat,outline=F,lty=1,staplewex=0,col=color,
#           xaxt="n",yaxt="n")
# })


