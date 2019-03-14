#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Plot de novo rates by SVTYPE for formal gnomAD analysis


###Set master parameters
options(stringsAsFactors=F,scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Read & clean de novo rate data
importDeNovoRates <- function(dat.in,svtypes.forAnalysis,PCRMINUS.samples){
  dat <- read.table(dat.in,header=T,comment.char="")
  colnames(dat)[1] <- "proband"
  filts <- c("PASS","FAIL")
  #Get unique variant categories
  cats <- unique(unlist(lapply(strsplit(colnames(dat)[grep("PASS_",colnames(dat),fixed=T)],split="PASS_"),function(l){l[[2]]})))
  #Consolidate by SVTYPE groups
  for(cat in cats){
    for(filt in filts){
      #Get indexes
      cnv.idxs <- sapply(c("DEL","DUP"),function(svtype){
        grep(paste(svtype,filt,cat,sep="_"),colnames(dat),fixed=T)
      })
      noncnv.idxs <- sapply(c("INS","INV","CPX","CTX"),function(svtype){
        grep(paste(svtype,filt,cat,sep="_"),colnames(dat),fixed=T)
      })
      #Create new columns
      cnv.newvals <- apply(dat[,cnv.idxs],1,sum)
      noncnv.newvals <- apply(dat[,noncnv.idxs],1,sum)
      dat <- cbind(dat,
                   cnv.newvals,
                   noncnv.newvals)
      colnames(dat)[c(ncol(dat)-1,ncol(dat))] <- c(paste("CNV",filt,cat,sep="_"),
                                                   paste("BCA",filt,cat,sep="_"))
    }
  }
  #Drop unnecessary columns
  cols.to.drop <- sort(unique(as.vector(sapply(c("DEL","DUP","INS","INV","CPX","CTX"),function(svtype){
    grep(paste(svtype,"_",sep=""),colnames(dat),fixed=T)
  }))))
  if(length(cols.to.drop)>0){
    dat <- dat[,-cols.to.drop]
  }
  #Calculate various summary rates
  for(svtype in c("ALL","CNV","BCA","BND")){
    for(filt in filts){
      #Mendelian variants: MENDELIAN_CHILD_HET, MENDELIAN_CHILD_HOM, UNTRANSMITTED_HET
      mendelian.cats <- c("MENDELIAN_CHILD_HET","MENDELIAN_CHILD_HOM","UNTRANSMITTED_HET")
      mendelian.idxs <- sapply(mendelian.cats,function(cat){
        grep(paste(svtype,filt,cat,sep="_"),colnames(dat),fixed=T)
      })
      mendelian.newvals <- apply(dat[,mendelian.idxs],1,sum)
      #Mendelian violations: APPARENT_DE_NOVO_HET, APPARENT_DE_NOVO_HOM, UNTRANSMITTED_HOM
      violation.cats <- c("APPARENT_DE_NOVO_HET","APPARENT_DE_NOVO_HOM","UNTRANSMITTED_HOM")
      violation.idxs <- sapply(violation.cats,function(cat){
        grep(paste(svtype,filt,cat,sep="_"),colnames(dat),fixed=T)
      })
      violation.newvals <- apply(dat[,violation.idxs],1,sum)
      #Mendelian violation rate
      mvr.denom.newvals <- (mendelian.newvals+violation.newvals)
      mvr.newvals <- violation.newvals/mvr.denom.newvals
      #Het de novo rate
      het_mendel <- dat[,grep(paste(svtype,filt,"MENDELIAN_CHILD_HET",sep="_"),colnames(dat),fixed=T)]
      het_denovo <- dat[,grep(paste(svtype,filt,"APPARENT_DE_NOVO_HET",sep="_"),colnames(dat),fixed=T)]
      het_dnr.denom.newvals <- (het_mendel+het_denovo)
      het_dnr.newvals <- het_denovo/het_dnr.denom.newvals
      #Hom genotype FDR
      hom_mendel <- dat[,grep(paste(svtype,filt,"MENDELIAN_CHILD_HOM",sep="_"),colnames(dat),fixed=T)]
      hom_denovo <- dat[,grep(paste(svtype,filt,"APPARENT_DE_NOVO_HOM",sep="_"),colnames(dat),fixed=T)]
      hom_fdr.denom.newvals <- (hom_mendel+hom_denovo)
      hom_fdr.newvals <- hom_denovo/hom_fdr.denom.newvals
      #Het genotype FNR
      parent_hom_untrans <- dat[,grep(paste(svtype,filt,"UNTRANSMITTED_HOM",sep="_"),colnames(dat),fixed=T)]
      het_fnr.denom.newvals <- (het_mendel+parent_hom_untrans)
      het_fnr.newvals <- parent_hom_untrans/het_fnr.denom.newvals
      #Add columns
      dat <- cbind(dat,
                   mvr.denom.newvals,
                   mvr.newvals,
                   het_dnr.denom.newvals,
                   het_dnr.newvals,
                   hom_fdr.denom.newvals,
                   hom_fdr.newvals,
                   het_fnr.denom.newvals,
                   het_fnr.newvals)
      new.colnames <- as.vector(sapply(c("MVR","HET_DNR","SPONT_HOMOZYGOTE_RATE","HET_FNR"),function(cat){
        c(paste(svtype,filt,cat,"SV_COUNT",sep="_"),
          paste(svtype,filt,cat,sep="_"))
      }))
      colnames(dat)[(ncol(dat)-7):ncol(dat)] <- new.colnames
    }
  }
  return(dat)
}
#Generate jitter residuals for a vector of values based on their density
sina.jitter <- function(vals){
  d <- density(vals)
  dv <- approx(d$x,d$y,xout=vals)
  dv <- dv$y/max(dv$y)
  dv.j <- sapply(1:length(vals),function(i){
    jitter(0,amount=dv[i])
  })
  return(dv.j)
}
#Generate sina points to add to existing plot
sina.plot <- function(vals,y.at,color,width=0.1){
  j <- (width*sina.jitter(vals))+y.at
  points(x=vals,y=j,pch=21,cex=0.25,col=color,bg=adjustcolor(color,alpha=0.3))
  segments(x0=median(vals),x1=median(vals),
           y0=y.at-width,y1=y.at+width,lwd=3)
  # points(x=median(vals),y=y.at,pch=18)
}
#Plot jitter of de novo rates by SVTYPE
plot.denovorates <- function(dat,sv.groups,sv.groups.labels,filters,filter.labels,filter.cols){
  #Prep plot area
  dnrs <- dat[,intersect(grep("_HET_DNR",colnames(dat),fixed=T),
                         grep("_HET_DNR_SV_COUNT",colnames(dat),fixed=T,invert=T))]
  counts <- dat[,grep("_HET_DNR_SV_COUNT",colnames(dat),fixed=T)]
  plot(x=c(0,max(dnrs,na.rm=T)),y=c(-0.1,-3.5),type="n",
       xaxt="n",xlab="",yaxt="n",ylab="")
  # abline(v=0.05,col="gray80",lty=2)
  # abline(v=0,col="gray80")
  x.ticks <- seq(0,max(axTicks(1)),by=0.05)
  axis(3,at=x.ticks,labels=NA)
  sapply(x.ticks,function(x){
    axis(3,at=x,tick=F,line=-0.5,cex.axis=0.8,
         labels=paste(round(100*x,1),"%",sep=""))
  })
  mtext(3,line=1.35,text=expression("Apparent Het." ~ italic("De Novo") ~"Rate"))
  #Add sina plots
  sapply(1:3,function(i){
    sapply(2:3,function(k){
      rates <- dnrs[,which(colnames(dnrs)==paste(sv.groups[i],filters[k],"HET_DNR",sep="_"))]
      totals <- counts[,which(colnames(counts)==paste(sv.groups[i],filters[k],"HET_DNR_SV_COUNT",sep="_"))]
      #Sina plot
      sina.plot(vals=rates,
                color=filter.cols[k],
                y.at=(-i+1)-(1.075-((4:1)[k]/4)))
      #Count & median DNR label
      med.count <- median(totals)
      med.dnr <- median(rates)
      axis(2,at=(-i+1)-(1.075-((4:1)[k]/4)),tick=F,line=-0.8,las=2,cex.axis=0.8,col.axis=filter.cols[k],
           labels=paste(prettyNum(round(med.count,0),big.mark=",")," / Child (",round(100*med.dnr,1),"%)",sep=""))
    })
    #Group label
    axis(2,at=-i+0.95,tick=F,line=-0.8,las=2,cex.axis=0.9,font=2,
         labels=sv.groups.labels[i])
  })
  #Add BNDs on their own
  rates <- dnrs$BND_FAIL_HET_DNR
  totals <- counts$BND_FAIL_HET_DNR_SV_COUNT
  #Sina plot
  sina.plot(vals=rates,
            color=filter.cols[3],
            y.at=(-4+1)-(1.075-((4:1)[2]/4)))
  #Count & median DNR label
  med.count <- median(totals)
  med.dnr <- median(rates)
  axis(2,at=(-4+1)-(1.075-((4:1)[2]/4)),tick=F,line=-0.8,las=2,cex.axis=0.8,col.axis=filter.cols[3],
       labels=paste(prettyNum(round(med.count,0),big.mark=",")," / Child (",round(100*med.dnr,1),"%)",sep=""))
  axis(2,at=-4+0.95,tick=F,line=-0.8,las=2,cex.axis=0.9,font=2,
       labels=sv.groups.labels[4])
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
              metavar="character")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog DENOVO_DATA OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
dat.in <- args$args[1]
PCRMINUS.samples.list.in <- args$args[2]
OUTDIR <- args$args[3]
prefix <- opts$prefix
svtypes.file <- opts$svtypes

# #Dev parameters (local)
# dat.in <- "~/scratch/gnomAD_v2_SV_MASTER.PCRMINUS_trios_deNovoStats.perSample_summary.txt.gz"
# PCRMINUS.samples.list.in <- "~/scratch/gnomAD_v2_SV_MASTER_wRelateds.PCRMINUS.samples.list"
# OUTDIR <- "~/scratch/gnomAD_v2_test_plots/"
# prefix <- "gnomAD_v2_SV_MASTER"
# svtypes.file <- "~/Desktop/Collins/Talkowski/code/sv-pipeline/ref/vcf_qc_refs/SV_colors.txt"

###Create output directory, if necessary
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

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

###Set analysis parameters
filters <- c("ANY","PASS","FAIL")
filter.labels <- c("All SVs","PASS-only","Non-PASS")
sv.groups <- c("ALL","CNV","BCA","BND")
sv.groups.labels <- c("All SVs","DEL/DUP","INS/INV/CPX/CTX","BND")
# filter.cols <- c("any"="gray60",
#                 "pass"="#4DAC26",
#                 "fail"="#D01C8B")
filter.cols <- c("ANY"="gray60",
                 "PASS"="#4DAC26",
                 "FAIL"="#AC26A1")
# svtypes.forAnalysis <- c("ALL","DEL","DUP","INS","INV","CPX","BND")
# svtype.colors <- c("gray50",
#                    svtypes$color[which(svtypes$svtype=="DEL")],
#                    svtypes$color[which(svtypes$svtype=="DUP")],
#                    svtypes$color[which(svtypes$svtype=="INS")],
#                    svtypes$color[which(svtypes$svtype=="INV")],
#                    svtypes$color[which(svtypes$svtype=="CPX")])


###Process input data
cat("NOW LOADING DATA\n")
PCRMINUS.samples <- as.character(read.table(PCRMINUS.samples.list.in,header=F)[,1])
dat <- importDeNovoRates(dat.in,svtypes.forAnalysis,PCRMINUS.samples)

###Plot de novo rates
png(paste(OUTDIR,"/",prefix,".het_deNovo_rates.png",sep=""),
    width=6*300,height=5*300,res=400)
par(mar=c(0.25,7,3,0.25),bty="n")
plot.denovorates(dat,sv.groups,sv.groups.labels,filters,filter.labels,filter.cols)
dev.off()

###Print sums of variants per trio (for supplementary QC table)
apply(dat[,intersect(grep("ALL_PASS",colnames(dat),fixed=T),
               grep("COUNT",colnames(dat),fixed=T))],
      2,median)
apply(dat[,tail(intersect(grep("ALL_PASS",colnames(dat),fixed=T),
                     grep("COUNT",colnames(dat),fixed=T,invert=T)),4)],
      2,median)

  