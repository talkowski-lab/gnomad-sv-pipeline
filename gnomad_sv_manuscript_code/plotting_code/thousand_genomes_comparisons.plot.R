#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Plot 1000 Genomes comparisons for formal gnomAD analysis


###Set master parameters
options(stringsAsFactors=F,scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Read & compute all overlaps for an input file list
import.overlaps <- function(dat.in,gpop){
  #Read list of file paths
  paths.list <- as.character(read.table(dat.in,header=F)[,1])
  tpops <- c("ALL","AFR","AMR","EAS","EUR")
  #Iterate over 1kG pops
  res <- lapply(tpops,function(tpop){
    #gnomAD sensitivity vs 1kGs
    path <- paths.list[grep(paste("gnomAD_",gpop,"_vs_1000G_Sudmant_",tpop,sep=""),paths.list,fixed=T)]
    sens.d <- read.table(path,header=T,comment.char="",sep="\t")
    colnames(sens.d)[1] <- "chr"
    sens.d[,ncol(sens.d)] <- "NO_OVR" #overwrite overlap method 3 to restrict overlaps
    sens.AF.pairs <- get.AFpairs(sens.d)
    sens.svtypeOverlap <- calc.svtypeOverlap(sens.d)
    sens.sizeOverlap <- calc.sizeOverlap(sens.d)
    sens.freqOverlap <- calc.freqOverlap(sens.d)
    #1kGs sensitivity vs gnomAD
    path <- paths.list[grep(paste("1000G_Sudmant_",tpop,"_vs_gnomAD_",gpop,sep=""),paths.list,fixed=T)]
    spec.d <- read.table(path,header=T,comment.char="",sep="\t")
    colnames(spec.d)[1] <- "chr"
    colnames(spec.d)[which(colnames(spec.d)=="SVTYPE")] <- "svtype"
    spec.d[,ncol(spec.d)] <- "NO_OVR" #overwrite overlap method 3 to restrict overlaps
    spec.d$length <- spec.d$end-spec.d$start
    spec.AF.pairs <- get.AFpairs(spec.d)
    spec.svtypeOverlap <- calc.svtypeOverlap(spec.d)
    spec.sizeOverlap <- calc.sizeOverlap(spec.d)
    spec.freqOverlap <- calc.freqOverlap(spec.d)
    #Format & return list of values
    out <- list("sens"=list("AF.pairs"=sens.AF.pairs,
                            "svtypeOverlap"=sens.svtypeOverlap,
                            "sizeOverlap"=sens.sizeOverlap,
                            "freqOverlap"=sens.freqOverlap),
                "spec"=list("AF.pairs"=spec.AF.pairs,
                            "svtypeOverlap"=spec.svtypeOverlap,
                            "sizeOverlap"=spec.sizeOverlap,
                            "freqOverlap"=spec.freqOverlap))
    return(out)
  })
  names(res) <- tpops
  return(res)
}
#Get closest AF between overlap conditions for each site
get.AFpairs <- function(dat){
  AF.pairs <- do.call("rbind", apply(dat[,which(colnames(dat) %in% c("AF","ovr1a","ovr1b",
                                                                     "ovr2a","ovr2b"))],
                                     1,function(vals){
                                       t.AF <- as.numeric(vals[1])
                                       g.AFs <- as.numeric(vals[-1][which(vals[-1] != "NO_OVR")])
                                       if(length(g.AFs)>0){
                                         AF.dist <- abs(t.AF-g.AFs)
                                         best <- head(g.AFs[which(AF.dist==min(AF.dist))],1)
                                         return(c("t"=t.AF,"g"=best))
                                       }else{
                                         return(NULL)
                                       }
                                     }))
}
#Allele frequency correlation plot between datasets
plot.AFcorr <- function(AF.pairs,color="gray50",t.pop=NULL,g.pop=NULL,
                        xaxis=T,yaxis=T,min.AF=0.01){
  #Prepare plot
  AF.pairs <- as.data.frame(AF.pairs[,c(2,1)])
  AF.pairs <- AF.pairs[which(AF.pairs[,1]>=min.AF | AF.pairs[,2]>=min.AF),]
  logscale.all <- log10(as.numeric(unlist(sapply(c(0:9),function(i){(1:9)*(10^i)}))))
  logscale.major <- 0:9
  major.labels <- sapply(logscale.major,function(i){expression(paste(i^"th"))})
  par(mar=c(3.2,1.5,1.5,3.2))
  plot(x=log10(c(0.01,1)),
       y=log10(c(0.01,1)),
       type="n",xaxt="n",yaxt="n",xlab="",ylab="")
  abline(a=0,b=1,col="gray10")
  if(xaxis==T){
    axis(1,at=-logscale.all,labels=NA,tck=-0.015,lwd=0.7)
    axis(1,at=-logscale.major,labels=NA,tck=-0.03,lwd=1.1)
    mtext(1,text=bquote("log"[10] ~ "(gnomAD" ~ .(g.pop) ~ "AF)"),line=2)
  }
  if(yaxis==T){
    axis(4,at=-logscale.all,labels=NA,tck=-0.015,lwd=0.7)
    axis(4,at=-logscale.major,labels=NA,tck=-0.03,lwd=1.1)
    mtext(4,text=bquote("log"[10] ~ "(1000G" ~ .(t.pop) ~ "AF)"),line=2)
    sapply(-logscale.major,function(i){
      # axis(1,at=i,labels=bquote('10'^.(i)),tick=F,line=-0.6,cex.axis=0.8)
      # axis(4,at=i,labels=bquote('10'^.(i)),tick=F,line=-0.4,cex.axis=0.8,las=2)
      axis(1,at=i,labels=i,tick=F,line=-0.6,cex.axis=0.8)
      axis(4,at=i,labels=i,tick=F,line=-0.4,cex.axis=0.8,las=2)
    })
  }
  
  #Add points
  pt.cex <- 0.4
  alpha <- 0.25
  points(x=log10(AF.pairs[,1]),y=log10(AF.pairs[,2]),
         pch=19,cex=pt.cex,lwd=0,
         col=adjustcolor(color,alpha=alpha))
  
  #Add stats
  # abline(lm(AF.pairs[,1] ~ AF.pairs[,2]),col="gray10",lty=2)
  AB.cor <- format(round(cor(AF.pairs[,1],AF.pairs[,2])^2,3),nsmall=3)
  text(x=par("usr")[1],y=par("usr")[4]-(0.085*(par("usr")[4]-par("usr")[3])),
       labels=bquote(italic(R)^2 == .(AB.cor)),cex=1.4,pos=4)
}
#Gather overlap by size deciles
calc.sizeOverlap <- function(dat){
  d <- quantile(dat$length,probs=seq(0,1,0.1))
  ovr <- as.data.frame(t(sapply(1:(length(d)-1),function(i){
    members <- which(dat$length>=d[i] & dat$length<d[i+1])
    denom <- length(members)
    strict <- length(intersect(members,which(dat$ovr1a != "NO_OVR" | dat$ovr2a != "NO_OVR")))
    loose <- length(intersect(members,which(dat$ovr1a != "NO_OVR" | dat$ovr1b != "NO_OVR" | dat$ovr2a != "NO_OVR" | dat$ovr2b != "NO_OVR" | dat$ovr3 != "NO_OVR")))
    return(c(strict/denom,loose/denom))
  })))
  colnames(ovr) <- c("strict","loose")
  return(ovr)
}
#Gather overlap by frequency deciles
calc.freqOverlap <- function(dat){
  singles <- which(dat$AF==min(dat$AF,na.rm=T))
  singles.strict <- length(intersect(singles,which(dat$ovr1a != "NO_OVR" | dat$ovr2a != "NO_OVR")))
  singles.loose <- length(intersect(singles,which(dat$ovr1a != "NO_OVR" | dat$ovr1b != "NO_OVR" | dat$ovr2a != "NO_OVR" | dat$ovr2b != "NO_OVR" | dat$ovr3 != "NO_OVR")))
  singles.ovr <- data.frame("strict"=singles.strict/length(singles),
                            "loose"=singles.loose/length(singles))
  d <- quantile(dat$AF[-singles],probs=seq(0,1,0.2))
  d.ovr <- as.data.frame(t(sapply(1:(length(d)-1),function(i){
    members <- which(dat$AF>=d[i] & dat$AF<d[i+1])
    denom <- length(members)
    strict <- length(intersect(members,which(dat$ovr1a != "NO_OVR" | dat$ovr2a != "NO_OVR")))
    loose <- length(intersect(members,which(dat$ovr1a != "NO_OVR" | dat$ovr1b != "NO_OVR" | dat$ovr2a != "NO_OVR" | dat$ovr2b != "NO_OVR" | dat$ovr3 != "NO_OVR")))
    return(c(strict/denom,loose/denom))
  })))
  colnames(d.ovr) <- c("strict","loose")
  ovr <- rbind(singles.ovr,d.ovr)
  return(ovr)
}
#Gather overlap by frequency deciles
calc.svtypeOverlap <- function(dat){
  #Get overlap per SVTYPE
  svt <- as.character(unique(dat$svtype))
  svt.order <- c("DEL","DUP","MCNV","INS","INV","CPX","BND")
  svt <- svt.order[which(svt.order %in% svt)]
  ovr <- as.data.frame(t(sapply(svt,function(s){
    members <- which(dat$svtype==s)
    denom <- length(members)
    strict <- length(intersect(members,which(dat$ovr1a != "NO_OVR" | dat$ovr2a != "NO_OVR")))
    loose <- length(intersect(members,which(dat$ovr1a != "NO_OVR" | dat$ovr1b != "NO_OVR" | dat$ovr2a != "NO_OVR" | dat$ovr2b != "NO_OVR" | dat$ovr3 != "NO_OVR")))
    return(c(strict/denom,loose/denom))
  })))
  colnames(ovr) <- c("strict","loose")
  #Get overlap for all variants and append to top of ovr
  denom.all <- nrow(dat)
  strict.all <- length(which(dat$ovr1a != "NO_OVR" | dat$ovr2a != "NO_OVR"))
  loose.all <- length(which(dat$ovr1a != "NO_OVR" | dat$ovr1b != "NO_OVR" | dat$ovr2a != "NO_OVR" | dat$ovr2b != "NO_OVR" | dat$ovr3 != "NO_OVR"))
  ovr <- rbind("ALL"=c(strict.all/denom.all,loose.all/denom.all),ovr)
  return(ovr)
}
#Generic function to plot decile bars with two-tone shading
plot.bars <- function(ovr,color,deciles=T,
                      add.singletons=F,
                      xlabel=NULL,ylabel="SV Captured",
                      xaxis=T,yaxis=T){
  par(bty="n",mar=c(2.5,5.5,1.5,0.5))
  plot(x=c(0,nrow(ovr)),y=c(0,1),type="n",
       xaxt="n",xlab="",yaxt="n",ylab="",yaxs="i")
  if(xaxis==T){
    if(deciles==T){
      if(add.singletons==T){
        axis(1,at=0.5,labels="S",tick=F,line=-0.7,cex.axis=0.8)
        sapply(2:nrow(ovr),function(i){
          if(i==2){
            axis(1,at=i-0.5,labels=bquote(1^'st'),
                 tick=F,line=-0.7,cex.axis=0.8)
          }else if(i==3){
            axis(1,at=i-0.5,labels=bquote(2^'nd'),
                 tick=F,line=-0.7,cex.axis=0.8)
          }else if(i==4){
            axis(1,at=i-0.5,labels=bquote(3^'rd'),
                 tick=F,line=-0.7,cex.axis=0.8)
          }else{
            axis(1,at=i-0.5,labels=bquote(.(i-1)^'th'),
                 tick=F,line=-0.7,cex.axis=0.8)
          }
        })
      }else{
        sapply(1:nrow(ovr),function(i){
          if(i==1){
            axis(1,at=i-0.5,labels=bquote(1^'st'),
                 tick=F,line=-0.7,cex.axis=0.8)
          }else if(i==2){
            axis(1,at=i-0.5,labels=bquote(2^'nd'),
                 tick=F,line=-0.7,cex.axis=0.8)
          }else if(i==3){
            axis(1,at=i-0.5,labels=bquote(3^'rd'),
                 tick=F,line=-0.7,cex.axis=0.8)
          }else{
            axis(1,at=i-0.5,labels=bquote(.(i)^'th'),
                 tick=F,line=-0.7,cex.axis=0.8)
          }
        })
      }
    }else{
      axis(1,at=(1:nrow(ovr))-0.5,labels=rownames(ovr),
           tick=F,las=2,line=-0.8,cex.axis=0.9)
    }
    mtext(1,line=1.5,text=xlabel)
  }
  if(yaxis==T){
    axis(2,at=seq(0,1,0.2),labels=NA)
    axis(2,at=seq(0,1,0.2),tick=F,line=-0.4,las=2,cex.axis=0.8,
         labels=paste(seq(0,100,20),"%",sep=""))
    mtext(2,line=2.5,text=ylabel)
  }
  rect(xleft=(1:nrow(ovr))-0.8,xright=(1:nrow(ovr))-0.2,
       ybottom=0,ytop=1,col="gray98",border="gray90")
  segments(x0=par("usr")[1],x1=nrow(ovr)-0.2,
           y0=0.5,y1=0.5,col="gray60")
  rect(xleft=(1:nrow(ovr))-0.8,xright=(1:nrow(ovr))-0.2,
       ybottom=0,ytop=ovr$strict,
       col=color,border=NA,bty="n")
  rect(xleft=(1:nrow(ovr))-0.8,xright=(1:nrow(ovr))-0.2,
       ybottom=ovr$strict,ytop=ovr$loose,
       col="white",border=NA,bty="n")
  rect(xleft=(1:nrow(ovr))-0.8,xright=(1:nrow(ovr))-0.2,
       ybottom=ovr$strict,ytop=ovr$loose,
       col=adjustcolor(color,alpha=0.5),border=NA,bty="n")
  rect(xleft=(1:nrow(ovr))-0.8,xright=(1:nrow(ovr))-0.2,
       ybottom=0,ytop=ovr$loose,col=NA)
  na.comps <- which(is.na(ovr$strict) | is.na(ovr$loose))
  if(length(na.comps)>0){
    rect(xleft=na.comps-0.8,xright=na.comps-0.2,
         ybottom=0,ytop=1,col="white")
    rect(xleft=na.comps-0.8,xright=na.comps-0.2,
         ybottom=0,ytop=1,col="gray80",density=15)
    text(x=na.comps-0.5,y=0.5,cex=0.7,labels="N/A",col="gray10",srt=90)
  }
}


################
###RSCRIPT BLOCK
################
require(optparse,quietly=T)
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
args <- parse_args(OptionParser(usage="%prog ALL AFR EAS EUR AMR OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 6){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
ALL.in <- args$args[1]
AFR.in <- args$args[2]
EAS.in <- args$args[3]
EUR.in <- args$args[4]
AMR.in <- args$args[5]
OUTDIR <- args$args[6]
prefix <- opts$prefix
svtypes.file <- opts$svtypes
pops.file <- opts$populations

# #Dev parameters (local)
# ALL.in <- "~/scratch/1kg_overlaps/1kG_ALL_comparisons.list"
# AFR.in <- "~/scratch/1kg_overlaps/1kG_AFR_comparisons.list"
# EAS.in <- "~/scratch/1kg_overlaps/1kG_EAS_comparisons.list"
# EUR.in <- "~/scratch/1kg_overlaps/1kG_EUR_comparisons.list"
# AMR.in <- "~/scratch/1kg_overlaps/1kG_AMR_comparisons.list"
# svtypes.file <- "~/Desktop/Collins/Talkowski/code/sv-pipeline/ref/vcf_qc_refs/SV_colors.txt"
# OUTDIR <- "~/scratch/gnomAD_v2_test_plots/"
# prefix <- "gnomAD_v2_SV_MASTER"
# pops.file <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/gnomAD_population_colors.txt"

###Create output directory, if necessary
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

###Process input data
cat("NOW LOADING AND CLEANING DATA\n")
gpops <- c("ALL","AFR","AMR","EAS","EUR")
tpops <- c("ALL","AFR","AMR","EAS","EUR")
dat.in.l <- list(ALL.in,AFR.in,AMR.in,EAS.in,EUR.in)
dat <- lapply(1:length(dat.in.l),function(i){
  import.overlaps(dat.in=unlist(dat.in.l[i]),gpop=gpops[i])
})
names(dat) <- gpops

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

###Sets populations & colors
if(!is.null(pops.file)){
  pops <- read.table(pops.file,sep="\t",header=T,comment.char="",check.names=F)
}

###Plot population strips of barplots (sens + spec) and AF correlations
sapply(gpops,function(gpop){
  #Assign variables
  if(gpop=="AMR"){
    tpop <- "AMR"
  }else{
    tpop <- gpop
  }
  if(gpop=="ALL"){
    p.color <- "gray20"
  }else{
    p.color <- pops$color[which(pops$pop==gpop)]
  }
  plot.data <- dat[[which(names(dat)==gpop)]]
  plot.data.sens <- plot.data[[which(names(plot.data)==tpop)]]$sens
  plot.data.spec <- plot.data[[which(names(plot.data)==tpop)]]$spec
  
  #Barplots of gnomAD sensitivity vs 1kG as truth
  png(paste(OUTDIR,"/",prefix,".gnomAD_",gpop,"_vs_1000Genomes_",tpop,".png",sep=""),
      height=6*350,width=3*400,res=400)
  par(mfrow=c(3,1))
  #overlap by SVTYPE barplots
  plot.bars(ovr=plot.data.sens$svtypeOverlap,color=p.color,
            deciles=F,ylabel="1000 Genomes SV\nFound in gnomAD")
  text(x=0.5,y=plot.data.sens$svtypeOverlap[1,2],pos=3,font=2,
       labels=paste(round(100*plot.data.sens$svtypeOverlap[1,2],0),"%",sep=""))
  #overlap by size barplots
  plot.bars(ovr=plot.data.sens$sizeOverlap,color=p.color,
            xlabel="SV Size Decile",ylabel="1000 Genomes SV\nFound in gnomAD")
  #overlap by AF barplots
  plot.bars(ovr=plot.data.sens$freqOverlap,color=p.color,
            add.singletons=T,
            xlabel="AF Quintile",ylabel="1000 Genomes SV\nFound in gnomAD")
  dev.off()
  
  #Barplots of 1kG sensitivity vs gnomAD as truth
  png(paste(OUTDIR,"/",prefix,".1000Genomes_",tpop,"_vs_gnomAD_",gpop,".png",sep=""),
      height=6*350,width=3*400,res=400)
  par(mfrow=c(3,1))
  #overlap by SVTYPE barplots
  plot.bars(ovr=plot.data.spec$svtypeOverlap,color=p.color,
            deciles=F,ylabel="gnomAD SV Found\nin 1000 Genomes")
  text(x=0.5,y=plot.data.spec$svtypeOverlap[1,2],pos=3,font=2,
       labels=paste(round(100*plot.data.spec$svtypeOverlap[1,2],0),"%",sep=""))
  #overlap by size barplots
  plot.bars(ovr=plot.data.spec$sizeOverlap,color=p.color,
            xlabel="SV Size Decile",ylabel="gnomAD SV Found\nin 1000 Genomes")
  #overlap by AF barplots
  plot.bars(ovr=plot.data.spec$freqOverlap,color=p.color,
            add.singletons=T,
            xlabel="AF Quintile",ylabel="gnomAD SV Found\nin 1000 Genomes")
  dev.off()
  
  #AF correlation scatterplots
  png(paste(OUTDIR,"/",prefix,".gnomAD_",gpop,"_AF_correlations_vs_1000Genomes.png",sep=""),
      height=3*length(tpops)*300,width=3*300,res=400)
  par(mfrow=c(length(tpops),1))
  sapply(tpops,function(tpop){
    if(tpop==gpop | (tpop=="AMR" & gpop=="AMR")){
      bwd <- 3
      if(tpop=="ALL"){
       psc.color <- "gray20"
      }else{
        psc.color <- p.color
      }
    }else{
      bwd <- 0.01
      psc.color <- "gray65"
    }
    plot.AFcorr(plot.data[[which(names(plot.data)==tpop)]]$sens$AF.pairs,
                color=psc.color,t.pop=tpop,g.pop=gpop)
    box(lwd=bwd)
  })
  dev.off()
})

#Make single AF correlation plot of ALL vs ALL for ED figure
png(paste(OUTDIR,"/",prefix,".gnomAD_ALL_AF_correlations_vs_1000Genomes_ALL.singlepanel.png",sep=""),
    height=5*300,width=5*300,res=400)
plot.AFcorr(dat$ALL$ALL$sens$AF.pairs,
            color="gray20",t.pop="ALL",g.pop="ALL")
dev.off()



