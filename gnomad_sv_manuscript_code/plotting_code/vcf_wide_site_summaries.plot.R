#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Generate VCF-wide distribution plots for formal gnomAD analysis


###Set master parameters
options(stringsAsFactors=F,scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Read & clean vcf2bed input
read.vcf2bed <- function(vcf2bed.in){
  #Read data
  dat <- read.table(vcf2bed.in,comment.char="",header=T,sep="\t")
  colnames(dat)[1] <- "chrom"
  #Restrict to sites with at least one observed alternative allele
  dat <- dat[union(which(dat$AC>0),grep("MULTIALLELIC",dat$FILTER,fixed=T)),]
  #Drop columns not being used (to save memory)
  cols.to.drop <- c("start","end","ALGORITHMS","CHR2","CPX_INTERVALS",
                    "END","EVIDENCE","SOURCE","STRANDS","UNRESOLVED_TYPE",
                    "LINCRNA__LOF","LINCRNA__DUP_LOF","LINCRNA__COPY_GAIN",
                    "LINCRNA__DUP_PARTIAL","LINCRNA__MSV_EXON_OVR",
                    "LINCRNA__INTRONIC","LINCRNA__INV_SPAN","LINCRNA__UTR")
  dat <- dat[,-which(colnames(dat) %in% cols.to.drop)]
  #Convert numeric columns
  numeric.columns <- sort(unique(c(grep("FREQ",colnames(dat),fixed=T),
                                   grep("AN",colnames(dat),fixed=T),
                                   grep("AC",colnames(dat),fixed=T),
                                   grep("AF",colnames(dat),fixed=T))))
  numeric.columns <- setdiff(numeric.columns,grep("SPAN",colnames(dat),fixed=T))
  dat[,numeric.columns] <- apply(dat[,numeric.columns],2,as.numeric)
  return(dat)
}

#Plot samples in study vs prior SV studies
plot.sampleSize <- function(priors,pops,pop.assignments){
  ordered.pops <- c("OTH","SAS","EUR","EAS","AMR","AFR")
  #Format input data
  plot.dat <- as.data.frame(t(sapply(priors$Study,function(s){
    s.dat <- sapply(ordered.pops,function(p){
      priors[which(priors$Study==s),which(colnames(priors)==paste(p,".pop",sep=""))]
    })
  })))
  g.dat <- sapply(ordered.pops,function(p){
    length(which(pop.assignments$pop==p))
  })
  # g.dat[which(names(g.dat)=="OTH")] <- g.dat[which(names(g.dat)=="OTH")]+(14891-sum(g.dat))
  plot.dat <- rbind(plot.dat,"gnomAD-SV"=g.dat)
  plot.dat <- apply(plot.dat,2,as.numeric)
  rownames(plot.dat) <- c(priors$Study,"gnomAD-SV")
  plot.dat <- plot.dat[order(-apply(plot.dat,1,sum)),]
  #Prep plot values
  plot.colors <-sapply(ordered.pops,function(p){
    if(p=="SAS"){
      pops$color[which(pops$pop=="OTH")]
    }else{
      pops$color[which(pops$pop==p)]
    }
  })
  #Plot
  par(bty="n")
  plot(x=c(0,1.15*15000),y=c(0,-nrow(plot.dat)),type="n",
       xaxt="n",xlab="",xaxs="i",yaxt="n",ylab="")
  x.at <- axTicks(1)
  if(length(x.at)>5){
    x.at <- x.at[seq(1,length(x.at),2)]
    x.at <- c(x.at,x.at[length(x.at)]+(x.at[2]-x.at[1]))
  }
  axis(1,at=x.at,labels=NA,tck=-0.03)
  sapply(x.at,function(k){
    axis(1,at=k,tick=F,line=-0.9,cex.axis=0.75,
         labels=paste(round(k/1000,0),"k",sep=""))
  })
  mtext(1,text="Sample Size",line=1,cex=0.85)
  sapply(1:nrow(plot.dat),function(i){
    rect(ybottom=-i+0.2,ytop=-i+0.8,
         xleft=c(0,cumsum(as.numeric(plot.dat[i,-ncol(plot.dat)]))),
         xright=cumsum(as.numeric(plot.dat[i,])),
         border=NA,bty="n",col=plot.colors)
    rect(ybottom=-i+0.2,ytop=-i+0.8,
         xleft=0,xright=sum(as.numeric(plot.dat[i,])),
         col=NA)
    par(xpd=T)
    text(x=sum(as.numeric(plot.dat[i,]))-(0.05*(par("usr")[2]-par("usr")[1])),
         y=-i+0.5,pos=4,cex=0.6,
         labels=prettyNum(sum(as.numeric(plot.dat[i,])),big.mark=","))
    par(xpd=F)
    if(rownames(plot.dat)[i]=="gnomAD-SV"){
      axis(2,at=-i+(2/3),tick=F,line=-0.8,cex.axis=0.75,las=2,
           labels="gnomAD-SV")
      axis(2,at=-i+(1/3),tick=F,line=-0.8,cex.axis=0.5,las=2,
           labels="This study")
    }else{
      year <- paste("(",priors$Year[which(priors$Study==rownames(plot.dat)[i])],")",sep="")
      axis(2,at=-i+0.5,tick=F,line=-0.8,cex.axis=0.75,las=2,
           labels=priors$Cohort[which(priors$Study==rownames(plot.dat)[i])])
      # if(rownames(plot.dat)[i]=="HehirKwa"){
      #   axis(2,at=-i+(1/3),tick=F,line=-0.8,cex.axis=0.6,las=2,
      #        labels=bquote("Hehir-Kwa" ~ italic("et al.") ~  .(year)))
      # }else{
      #   axis(2,at=-i+(1/3),tick=F,line=-0.8,cex.axis=0.6,las=2,
      #        labels=bquote(.(priors$Author[which(priors$Study==rownames(plot.dat)[i])]) ~ italic("et al.") ~ .(year)))
      # }
    }
  })
}

#Plot total SV sites discovered per study vs prior SV studies
plot.SVcountsByStudy <- function(priors,dat,svtypes){
  #Format input data
  plot.dat <- as.data.frame(t(sapply(priors$Study,function(s){
    s.dat <- sapply(c("DEL","DUP","MCNV","INS","INV","CPX","OTH"),function(p){
      if(p!="CPX"){
        priors[which(priors$Study==s),which(colnames(priors)==p)]
      }else{
        return(0)
      }
    })
  })))
  g.dat <- sapply(c("DEL","DUP","MCNV","INS","INV","CPX"),function(p){
    length(which(dat$SVTYPE==p))
  })
  g.dat <- c(g.dat,"OTH"=length(which(!(dat$SVTYPE %in% c("DEL","DUP","MCNV","INS","INV","CPX")))))
  plot.dat <- rbind(plot.dat,"gnomAD-SV"=g.dat)
  plot.dat <- apply(plot.dat,2,as.numeric)
  rownames(plot.dat) <- c(priors$Study,"gnomAD-SV")
  # plot.dat <- plot.dat[order(-apply(plot.dat,1,sum)),]
  plot.dat <- plot.dat[match(c("gnomAD-SV","Sudmant","HehirKwa","Chiang"),rownames(plot.dat)),]
  #Prep plot values
  plot.colors <-sapply(c("DEL","DUP","MCNV","INS","INV","CPX","OTH"),function(p){
    svtypes$color[which(svtypes$svtype==p)]
  })
  #Plot
  par(bty="n")
  plot(x=c(0,1.15*max(apply(plot.dat,1,sum))),y=c(0,-nrow(plot.dat)),type="n",
       xaxt="n",xlab="",xaxs="i",yaxt="n",ylab="")
  # x.at <- axTicks(1)
  x.at <- seq(0,1000000,250000)
  if(length(x.at)>5){
    x.at <- x.at[seq(1,length(x.at),2)]
    x.at <- c(x.at,x.at[length(x.at)]+(x.at[2]-x.at[1]))
  }
  axis(1,at=x.at,labels=NA,tck=-0.03)
  sapply(x.at,function(k){
    axis(1,at=k,tick=F,line=-0.9,cex.axis=0.75,
         labels=paste(round(k/1000,0),"k",sep=""))
  })
  mtext(1,text="SV Sites Discovered",line=1,cex=0.85)
  sapply(1:nrow(plot.dat),function(i){
    rect(ybottom=-i+0.2,ytop=-i+0.8,
         xleft=c(0,cumsum(as.numeric(plot.dat[i,-ncol(plot.dat)]))),
         xright=cumsum(as.numeric(plot.dat[i,])),
         border=NA,bty="n",col=plot.colors)
    rect(ybottom=-i+0.2,ytop=-i+0.8,
         xleft=0,xright=sum(as.numeric(plot.dat[i,])),
         col=NA)
    par(xpd=T)
    text(x=sum(as.numeric(plot.dat[i,]))-(0.05*(par("usr")[2]-par("usr")[1])),
         y=-i+0.5,pos=4,cex=0.6,
         labels=prettyNum(sum(as.numeric(plot.dat[i,])),big.mark=","))
    par(xpd=F)
    if(rownames(plot.dat)[i]=="gnomAD-SV"){
      axis(2,at=-i+(2/3),tick=F,line=-0.8,cex.axis=0.75,las=2,
           labels="gnomAD-SV")
      axis(2,at=-i+(1/3),tick=F,line=-0.8,cex.axis=0.5,las=2,
           labels="This study")
    }else{
      year <- paste("(",priors$Year[which(priors$Study==rownames(plot.dat)[i])],")",sep="")
      axis(2,at=-i+0.5,tick=F,line=-0.8,cex.axis=0.75,las=2,
           labels=priors$Cohort[which(priors$Study==rownames(plot.dat)[i])])
      # if(rownames(plot.dat)[i]=="HehirKwa"){
      #   axis(2,at=-i+(1/3),tick=F,line=-0.8,cex.axis=0.6,las=2,
      #        labels=bquote("Hehir-Kwa" ~ italic("et al.") ~  .(year)))
      # }else{
      #   axis(2,at=-i+(1/3),tick=F,line=-0.8,cex.axis=0.6,las=2,
      #        labels=bquote(.(priors$Author[which(priors$Study==rownames(plot.dat)[i])]) ~ italic("et al.") ~ .(year)))
      # }
    }
  })
}

#Plot simple log-scaled bars of SV by count
plot.totalCounts <- function(dat,svtypes,thousandG=F){
  #Gather data
  counts <- lapply(svtypes$svtype[which(svtypes$svtype!="OTH")],function(svtype){
    return(log10(length(which(dat$SVTYPE==svtype))))
  })
  names(counts) <- svtypes$svtype[which(svtypes$svtype!="OTH")]
  counts$OTH <- log10(length(which(!(dat$SVTYPE %in% svtypes$svtype[which(svtypes$svtype!="OTH")]))))
  counts <- lapply(counts,function(val){if(is.infinite(val)){val <- 0}; return(val)})
  counts <- unlist(counts)
  sudmant.counts <- log10(c(42279,6025,2929,16631+168,786,NA,NA))
  names(sudmant.counts) <- names(counts)
  #Plot
  log.minor <- log10(as.numeric(sapply(0:10,function(i){(1:8)*(10^i)})))
  log.major <- 1:9
  par(mar=c(3,3.5,3,0.5),bty="n")
  plot(x=c(0,2*nrow(svtypes)),y=c(0,1.02*max(counts)),type="n",
       xaxt="n",xlab="",yaxt="n",ylab="",yaxs="i")
  if(thousandG==T){
    rect(xleft=seq(1,2*length(counts),2)-0.8,xright=seq(1,2*length(counts),2),
         ybottom=0,ytop=sudmant.counts,
         lwd=0.5,col="gray70")
    rect(xleft=seq(2,2*length(counts),2)-1,xright=seq(2,2*length(counts),2)-0.2,
         ybottom=0,ytop=counts,
         lwd=0.5,col=svtypes$color)
  }else{
    rect(xleft=seq(1,2*length(counts),2)-0.8,xright=seq(2,2*length(counts),2)-0.2,
         ybottom=0,ytop=counts,
         lwd=0.5,col=svtypes$color)
  }
  segments(x0=seq(1,2*length(counts),2)-0.8,
           x1=seq(2,2*length(counts),2)-0.2,
           y0=par("usr")[3],y1=par("usr")[3],lwd=2)
  axis(1,at=seq(1,2*length(counts),2),tick=F,line=-0.9,
       labels=svtypes$svtype,cex.axis=0.8,las=2)
  axis(2,at=log.minor,labels=NA,tck=-0.02,lwd=0.9)
  axis(2,at=log.major,labels=NA,tck=-0.04,lwd=1.1)
  axis(2,at=1:6,tick=F,line=-0.5,cex.axis=0.8,las=2,
       labels=c("10","100","1k","10k","100k","1M"))
  mtext(2,text="SV Sites Discovered",line=2)
  sapply(1:length(counts),function(i){
    if(thousandG==T){
      if(is.na(sudmant.counts[i])){
        sudmant.y <- 0
        sudmant.lab <- 0
      }else{
        sudmant.y <- sudmant.counts[i]
        sudmant.lab <- 10^sudmant.counts[i]
      }
      text(x=2*i-1.4,y=sudmant.y+(0.025*(par("usr")[4]-par("usr")[3])),
           labels=prettyNum(sudmant.lab,big.mark=","),
           srt=90,adj=0,xpd=T,cex=0.7,col="gray60")
      text(x=2*i-0.55,y=counts[i]+(0.025*(par("usr")[4]-par("usr")[3])),
           labels=prettyNum(10^counts[i],big.mark=","),
           srt=90,adj=0,xpd=T,cex=0.7)
    }else{
      text(x=2*i-1,y=counts[i]+(0.025*(par("usr")[4]-par("usr")[3])),
           labels=prettyNum(10^counts[i],big.mark=","),
           srt=90,adj=0,xpd=T,cex=0.7)
    }
  })
}

#Build summary table of counts, sizes, frequencies
build.table <- function(dat.wrelateds,dat.all,dat){
  svtypes.fortable <- c("DEL","DUP","MCNV","INS","INV","CPX","CTX","BND")
  #Total count of sites
  sites.count <- sapply(svtypes.fortable,function(svtype){
    length(which(dat.wrelateds$SVTYPE==svtype))
  })
  sites.count <- c(sum(sites.count),sites.count)
  #Pct of PASS sites
  pass.count <- sapply(svtypes.fortable,function(svtype){
    length(which(dat.wrelateds$SVTYPE==svtype & (dat.wrelateds$FILTER=="PASS" | dat.wrelateds$FILTER=="MULTIALLELIC")))
  })
  pass.count <- c(sum(pass.count),pass.count)
  pct.pass <- pass.count/sites.count
  #Count of PASS sites in unrelated genomes
  pass.count.unrelated <- sapply(svtypes.fortable,function(svtype){
    length(which(dat.all$SVTYPE==svtype & (dat.all$FILTER=="PASS" | dat.all$FILTER=="MULTIALLELIC")))
  })
  pass.count.unrelated <- c(sum(pass.count.unrelated),pass.count.unrelated)
  #Count of variants by freq bin
  count.byfreq <- t(sapply(svtypes.fortable,function(svtype){
    singletons <- length(which(dat$SVTYPE==svtype & dat$AC==1))
    rare <- length(which(dat$SVTYPE==svtype & dat$AC>1 & dat$AF<0.01))
    common <- length(which(dat$SVTYPE==svtype & (dat$AF>=0.01 | dat$FILTER=="MULTIALLELIC")))
    c(singletons,rare,common)
  }))
  colnames(count.byfreq) <- c("Singleton","Rare","Common")
  count.byfreq <- rbind(c(length(which(dat$AC==1)),
                          length(which(dat$AC>1 & dat$AF<0.01)),
                          length(which(dat$AF>=0.01))),
                        count.byfreq)
  #Median size
  med.size <- sapply(svtypes.fortable,function(svtype){
    median(dat$SVLEN[which(dat$SVTYPE==svtype)])
  })
  med.size <- c(median(dat$SVLEN,na.rm=T),med.size)
  med.size[which(med.size<1)] <- NA
  #Format output table
  table.out <- cbind(data.frame("Abbrev."=c("ALL",svtypes.fortable),
                          "Total Sites Discovered"=sites.count,
                          "Pct.Pass"=pct.pass,
                          "Passing Sites"=pass.count.unrelated),
                     count.byfreq,
                     data.frame("Med.Size"=med.size))
  return(table.out)
}

#Plot rolling mean size distribution lines
plot.sizes <- function(dat,svtypes,step=0.02,xlims=c(50,10000000)){
  #Iterate over SVTYPES and get log-scaled count of SV by class per bin
  require(zoo,quietly=T)
  size.bins <- seq(log10(50),log10(50000000),step)
  plot.dat <- as.data.frame(sapply(svtypes$svtype[which(svtypes$svtype!="OTH")],function(svtype){
    res <- sapply(1:length(size.bins),function(i){
      nrow(dat[which(dat$SVLEN>=10^size.bins[i] & dat$SVLEN<10^(size.bins[i]+step) & dat$SVTYPE==svtype),])
    })
    return(res)
  }))
  plot.dat <- apply(plot.dat,2,function(vals){rollapply(vals,sum,width=5,partial=T)})
  plot.dat <- apply(plot.dat,2,log10)
  plot.dat[which(is.infinite(plot.dat))] <- 0
  plot.dat <- as.data.frame(plot.dat)
  
  #Prep plot area
  par(bty="n",mar=c(2,3.5,0.5,1))
  L1.peak <- log10(6000)
  SVA.peak <- log10(1200)
  ALU.peak <- log10(280)
  log.minor <- log10(as.numeric(sapply(0:10,function(i){(1:8)*(10^i)})))
  log.mid <- log10(as.numeric(sapply(0:10,function(i){c(1,5)*(10^i)})))
  log.major <- 1:9
  plot(x=log10(xlims),y=c(0,1.05*max(plot.dat)),type="n",
       xaxt="n",xaxs="i",xlab="",yaxt="n",ylab="",yaxs="i")
  # abline(v=c(L1.peak,SVA.peak,ALU.peak),lty=2,col="gray80")
  # sapply(c(L1.peak,SVA.peak,ALU.peak),function(pos){
  #   axis(3,at=pos,labels=NA)
  # })
  # abline(h=log.minor,v=log.minor,col="gray97",lwd=0.5)
  # abline(h=log.mid,v=log.mid,col="gray94",lwd=0.75)
  # abline(h=log.mid,v=log.major,col="gray91")
  # rect(xleft=c(L1.peak,SVA.peak,ALU.peak)-(0.01*(par("usr")[2]-par("usr")[1])),
  #      xright=c(L1.peak,SVA.peak,ALU.peak)+(0.01*(par("usr")[2]-par("usr")[1])),
  #      ybottom=par("usr")[3],ytop=par("usr")[4],
  #      border=NA,bty="n",col=adjustcolor("gray85",alpha=0.5))
  
  #Plot lines per SV class
  sapply(1:ncol(plot.dat),function(i){
    points(x=size.bins,y=plot.dat[,i],lwd=2.5,col=svtypes$color[i],type="l")
  })
  
  #Add axes
  axis(1,at=log.minor,labels=NA,tck=-0.0175,lwd=0.9)
  # axis(1,at=log.mid,labels=NA,tck=-0.02)
  axis(1,at=log.major,labels=NA,tck=-0.035,lwd=1.1)
  sapply(1:6,function(i){
    axis(1,at=i+1,tick=F,line=-0.9,cex.axis=0.65,
         labels=c("100bp","1kb","10kb","100kb","1Mb","10Mb")[i])
  })
  
  mtext(1,text="SV Size",line=1)
  axis(2,at=log.minor,labels=NA,tck=-0.0175,lwd=0.9)
  # axis(2,at=log.mid,labels=NA,tck=-0.02)
  axis(2,at=log.major,labels=NA,tck=-0.035,lwd=1.1)
  axis(2,at=1:6,tick=F,line=-0.6,cex.axis=0.8,las=2,
       labels=c("10","100","1k","10k","100k","1M"))
  mtext(2,text="SV Discovered",line=2)
  
  #Legend
  # legend("topright",pch=19,col=svtypes$color,
  #        legend=svtypes$svtype[which(svtypes$svtype!="OTH")],
  #        border=NA,bty="n",cex=0.8)
}

#Plot distribution of allele counts by SVTYPE
plot.freqByType <- function(dat,ymin=NULL,axlabel.cex=1){
  #Gather data
  exclude <- grep("MULTIALLELIC",dat$FILTER,fixed=T)
  if(length(exclude)>0){
    dat <- dat[-exclude,]
  }
  dat <- dat[which(dat$SVTYPE!="MCNV" & !(dat$chrom %in% c("X","Y"))),]
  AC.cutoffs <- c(1:9,seq(10,90,10),seq(100,900,100),seq(1000,25000,1000))
  AF.dat <- as.data.frame(sapply(c("INS","DUP","DEL","INV","CPX"),function(svtype){
    sapply(AC.cutoffs,function(AC){
      length(which(dat$AC[which(dat$SVTYPE==svtype)]<=AC))/length(which(dat$SVTYPE==svtype))
    })
  }))
  # AF.dat$OTH <- as.numeric(sapply(AC.cutoffs,function(AC){
  #   return(length(which(dat$AC[which(!(dat$SVTYPE %in% svtypes$svtype[which(svtypes$svtype!="OTH")]))]<=AC))/length(which(!(dat$SVTYPE %in% svtypes$svtype[which(svtypes$svtype!="OTH")]))))
  # }))
  #Plot
  if(is.null(ymin)){
    ymin <- log10(floor(100*min(AF.dat,na.rm=T))/100)
  }else{
    ymin <- log10(ymin)
  }
  AF.dat <- as.data.frame(apply(AF.dat,2,log10))
  xrange <- c(-0.25,max(log10(AC.cutoffs)))
  common.threshold <- min(as.numeric(dat$AC[which(as.numeric(dat$AF)>=0.01 & as.numeric(dat$AN)==max(dat$AN,na.rm=T))]),na.rm=T)
  par(mar=c(2,3.25,0.5,0.5),bty="n")
  plot(x=xrange,y=c(ymin,log10(1)),type="n",
       xaxt="n",xlab="",yaxt="n",ylab="")
  # rect(xleft=-0.1,xright=0.1,ybottom=par("usr")[3],ytop=par("usr")[4],
  #      border=NA,bty="n",col="gray90")
  # abline(v=log10(common.threshold),col="gray80")
  segments(x0=log10(common.threshold),x1=log10(common.threshold),
           y0=par("usr")[3],y1=log10(1),col="gray80")
  text(x=log10(common.threshold)+(0.025*(par("usr")[2]-par("usr")[1])),
       y=par("usr")[3]+(0.08*(par("usr")[4]-par("usr")[3])),
       labels="Rare\n(AF<1%)",pos=2,cex=0.7)
  text(x=log10(common.threshold)-(0.025*(par("usr")[2]-par("usr")[1])),
       y=par("usr")[3]+(0.08*(par("usr")[4]-par("usr")[3])),
       labels="Common\n(AF>1%)",pos=4,cex=0.7)
  sapply(c("INS","DUP","DEL","INV","CPX"),function(svtype){
    points(x=log10(AC.cutoffs),y=AF.dat[,which(colnames(AF.dat)==svtype)],
           col=svtypes$color[which(svtypes$svtype==svtype)],type="l",lwd=3)
  })
  sapply(c("INS","DUP","DEL","INV","CPX"),function(svtype){
    # points(x=log10(AC.cutoffs)[1],y=AF.dat[,which(colnames(AF.dat)==svtype)][1],
    #        col="black",pch="-",cex=1.4)
    # points(x=log10(AC.cutoffs)[1],y=AF.dat[,which(colnames(AF.dat)==svtype)][1],
    #        col=svtypes$color[which(svtypes$svtype==svtype)],pch="-",cex=1.2)
    rect(xleft=log10(AC.cutoffs)[1]-(0.03*(par("usr")[2]-par("usr")[1])),
         xright=log10(AC.cutoffs)[1]+(0.03*(par("usr")[2]-par("usr")[1])),
         ybottom=AF.dat[,which(colnames(AF.dat)==svtype)][1]-(0.01*(par("usr")[4]-par("usr")[3])),
         ytop=AF.dat[,which(colnames(AF.dat)==svtype)][1]+(0.01*(par("usr")[4]-par("usr")[3])),
         col=svtypes$color[which(svtypes$svtype==svtype)])
  })
  logscale.all <- log10(as.numeric(unlist(sapply(c(0:9),function(i){(1:9)*(10^i)}))))
  logscale.major <- 0:9
  axis(1,at=logscale.all,labels=NA,tck=-0.015,lwd=0.7)
  axis(1,at=logscale.major,labels=NA,tck=-0.03,lwd=1.1)
  sapply(1:5,function(i){
    axis(1,at=i-1,tick=F,line=-0.9,cex.axis=0.8,
         labels=c("1","10","100","1k","10k")[i])
  })
  mtext(1,text="Allele Count",line=1,cex=axlabel.cex)
  logscale.pct.all <- log10((1:100)/100)
  logscale.pct.major <- log10(seq(10,100,10)/100)
  axis(2,at=logscale.pct.all,labels=NA,tck=-0.015,lwd=0.9)
  axis(2,at=logscale.pct.major,labels=NA,tck=-0.03,lwd=1.1)
  axis(2,at=logscale.pct.major,tick=F,line=-0.5,cex.axis=0.8,las=2,
       labels=paste(seq(10,100,10),"%",sep=""))
  mtext(2,text="Fraction of SV",line=2.25,cex=axlabel.cex)
}

#Get proportion of singletons and 95% CI for a set of variants
calc.fracSingletons.singleClass <- function(dat,boot.n=100,conf=0.95){
  ACs <- dat$AC
  ACs <- ACs[which(!is.na(ACs) & ACs>0)]
  helper.getFracSingle <- function(ACs,indices){length(which(ACs[indices]==1))/length(ACs[indices])}
  point.est <- helper.getFracSingle(ACs,indices=1:length(ACs))
  calc.ci <- function(ACs,n,conf){
    set.seed(0)
    boot.obj <- boot(data=ACs,statistic=helper.getFracSingle,R=n)
    ci <- boot.ci(boot.obj,conf=conf,type="basic")$basic[4:5]
    return(ci)
  }
  ci <- calc.ci(ACs,n=boot.n,conf=conf)
  return(c(point.est,ci))
}

#Gather fraction of singletons for a single SV class in bins by size
calc.fracSingletons.singleClass.binned <- function(dat,sd_sr_cov,svtype,
                                                   max.cov=1/3,include.ins=F,
                                                   d=c(0,1000,10000,100000,1000000000)){
  if(svtype!="ALL"){
    tmpdat <- dat[which(dat$SVTYPE==svtype & dat$FILTER=="PASS"),]
  }else{
    tmpdat <- dat[which(dat$FILTER=="PASS"),]
  }
  if(include.ins==F){
    tmpdat <- tmpdat[which(tmpdat$SVTYPE!="INS"),]
  }
  tmpdat <- tmpdat[which(tmpdat$name %in% sd_sr_cov[which(sd_sr_cov[,2]<=max.cov),1]),]
  if(is.null(d)){
    d <- quantile(tmpdat$SVLEN,probs=seq(0,1,0.1))
  }
  d.singles <- t(sapply(2:length(d),function(i){
    calc.fracSingletons.singleClass(dat=tmpdat[which(tmpdat$SVLEN>=d[i-1] & tmpdat$SVLEN<d[i]),])
  }))
  return(d.singles)
}

#Plot fraction of singletons by size for a single SV class
plot.fracSingletons.bySize <- function(dat,sd_sr_cov,svtype,
                                       max.cov=0.3,include.ins=T,
                                       d=c(0,1000,10000,100000,1000000,1000000000),
                                       d.labels=c("<1kb","1-10kb","10-100kb","100kb-1Mb",">1Mb"),
                                       color="black",title=NULL,mar=c(3,2.75,0.5,0.5),ylab.cex=1.2,
                                       ylims=NULL){
  #Get plot data
  p.dat <- calc.fracSingletons.singleClass.binned(dat,sd_sr_cov,svtype,
                                                  max.cov,include.ins,d)
  exclude <- which(apply(p.dat,1,function(vals){any(vals==0)}))
  if(length(exclude)>0){
    p.dat[exclude,1:3] <- NA
  }
  p.dat[which(p.dat>1)] <- 1
  p.dat[which(p.dat<0)] <- 0
  if(is.null(ylims)){
    ylims <- range(p.dat,na.rm=T)
  }
  #Prep plot area
  par(mar=mar,bty="n")
  plot(x=c(0.25,length(d)-1.25),y=ylims,type="n",
       xlab="",xaxt="n",ylab="",yaxt="n")
  #Plot dots and CIs
  segments(x0=(1:(length(d)-1))-0.5,
           x1=(1:(length(d)-1))-0.5,
           y0=p.dat[,2],y1=p.dat[,3],
           lend="round",lwd=2,col=color)
  points(x=(2:length(d))-1.5,y=p.dat[,1],pch=19,col=color)
  #Add axes
  text(x=(1:(length(d)-1))+0.1,y=par("usr")[3]-(0.025*(par("usr")[4]-par("usr")[3])),
       xpd=T,srt=45,labels=d.labels,cex=0.85,pos=2)
  axis(2,at=seq(0,1,0.1),labels=NA,tck=-0.03)
  axis(2,at=seq(0,1,0.1),tick=F,labels=paste(seq(0,100,10),"%",sep=""),
       las=2,cex.axis=0.7,line=-0.7)
  mtext(2,line=1.7,text="Singleton Proportion",cex=ylab.cex)
  mtext(3,line=0,text=paste("     ",title,sep=""),cex=1.2*ylab.cex)
}


# #Plot distribution of allele counts by SVLEN decile
# plot.freqBySize <- function(dat){
#   #Gather data
#   AC.cutoffs <- c(1:9,seq(10,90,10),seq(100,900,100),seq(1000,25000,1000))
#   AF.dat.sub <- dat[which(dat$SVTYPE %in% c("DEL","DUP","INV","CPX")),]
#   size.cutoffs <- quantile(AF.dat.sub$SVLEN,probs=(0:5)/5)
# 
#   # AF.dat.sub <- dat[which(!(dat$SVTYPE %in% c("BND","CTX"))),]
#   # AF.dat.sub <- AF.dat.sub[-grep("INS:ME",AF.dat.sub$svtype,fixed=T),]
#   AF.dat <- as.data.frame(sapply(1:5,function(i){
#     sapply(AC.cutoffs,function(AC){
#       length(which(AF.dat.sub$AC[which(AF.dat.sub$SVLEN>size.cutoffs[i] & AF.dat.sub$SVLEN<=size.cutoffs[i+1])]<=AC))/length(which(AF.dat.sub$SVLEN>size.cutoffs[i] & AF.dat.sub$SVLEN<=size.cutoffs[i+1]))
#     })
#   }))
#   #Plot
#   ymin <- log10(floor(100*min(AF.dat,na.rm=T))/100)
#   AF.dat <- as.data.frame(apply(AF.dat,2,log10))
#   xrange <- c(-0.25,max(log10(AC.cutoffs)))
#   common.threshold <- min(dat$AC[which(dat$AF>=0.01 & dat$AN==max(dat$AN))])
#   col.pal <- colorRampPalette(c("#440154","#365C8C","#25A584","#FDE725"))(5)
#   par(mar=c(2.5,3.25,0.5,0.5),bty="n")
#   plot(x=xrange,y=c(ymin,log10(1)),type="n",
#        xaxt="n",xlab="",yaxt="n",ylab="")
#   # rect(xleft=-0.1,xright=0.1,ybottom=par("usr")[3],ytop=par("usr")[4],
#   #      border=NA,bty="n",col="gray90")
#   # abline(v=log10(common.threshold),col="gray80")
#   segments(x0=log10(common.threshold),x1=log10(common.threshold),
#            y0=par("usr")[3],y1=log10(1),col="gray80")
#   text(x=log10(common.threshold),y=par("usr")[3]+(0.08*(par("usr")[4]-par("usr")[3])),
#        labels="Rare\n(AF<1%)",pos=2,cex=0.7)
#   text(x=log10(common.threshold),y=par("usr")[3]+(0.08*(par("usr")[4]-par("usr")[3])),
#        labels="Common\n(AF>1%)",pos=4,cex=0.7)
#   sapply(1:5,function(i){
#     points(x=log10(AC.cutoffs),y=AF.dat[,i],
#            col=col.pal[i],type="l",lwd=3)
#   })
#   sapply(c("OTH","INS","DUP","DEL","INV","CPX"),function(svtype){
#     points(x=log10(AC.cutoffs)[1],y=AF.dat[,which(colnames(AF.dat)==svtype)][1],
#            bg=svtypes$color[which(svtypes$svtype==svtype)],pch=21,cex=1.2)
#   })
#   logscale.all <- log10(as.numeric(unlist(sapply(c(0:9),function(i){(1:9)*(10^i)}))))
#   logscale.major <- 0:9
#   axis(1,at=logscale.all,labels=NA,tck=-0.015,lwd=0.7)
#   axis(1,at=logscale.major,labels=NA,tck=-0.03,lwd=1.1)
#   axis(1,at=c(0:4),tick=F,line=-0.7,cex.axis=0.8,
#        labels=c("1","10","100","1k","10k"))
#   mtext(1,text="Maximum Allele Count",line=1.3)
#   logscale.pct.all <- log10((1:100)/100)
#   logscale.pct.major <- log10(seq(10,100,10)/100)
#   axis(2,at=logscale.pct.all,labels=NA,tck=-0.015,lwd=0.9)
#   axis(2,at=logscale.pct.major,labels=NA,tck=-0.03,lwd=1.1)
#   axis(2,at=logscale.pct.major,tick=F,line=-0.5,cex.axis=0.8,las=2,
#        labels=paste(seq(10,100,10),"%",sep=""))
#   mtext(2,text="Fraction of SV Sites",line=2.25)
# }

#Gather Hardy-Weinberg data
makeHWEmat <- function(dat,pop=NULL){
  sub.dat <- dat[which(!(dat$chrom %in% c("X","Y"))),]
  cols.to.exclude <- grep("MULTIALLELIC",sub.dat$FILTER,fixed=T)
  if(length(cols.to.exclude)>0){
    sub.dat <- sub.dat[-cols.to.exclude,]
  }
  # sub.dat <- sub.dat[-grep("PESR_GT_OVERDISPERSION",sub.dat$FILTER,fixed=T),]
  # sub.dat <- sub.dat[-grep("UNRESOLVED",sub.dat$FILTER,fixed=T),]
  if(!is.null(pop)){
    n.genos.idx <- which(colnames(sub.dat)==paste(pop,"_N_BI_GENOS",sep=""))
    n.homref.idx <- which(colnames(sub.dat)==paste(pop,"_N_HOMREF",sep=""))
    n.het.idx <- which(colnames(sub.dat)==paste(pop,"_N_HET",sep=""))
    n.homalt.idx <- which(colnames(sub.dat)==paste(pop,"_N_HOMALT",sep=""))
  }else{
    n.genos.idx <- which(colnames(sub.dat)=="N_BI_GENOS")
    n.homref.idx <- which(colnames(sub.dat)=="N_HOMREF")
    n.het.idx <- which(colnames(sub.dat)=="N_HET")
    n.homalt.idx <- which(colnames(sub.dat)=="N_HOMALT")
  }
  sub.dat <- sub.dat[which(sub.dat[,n.genos.idx]>0),]
  sub.dat <- sub.dat[which(apply(sub.dat[,c(n.het.idx,n.homalt.idx)],1,sum)>0),]
  HWE.mat <- data.frame("AA"=as.numeric(sub.dat[,n.homref.idx]),
                        "AB"=as.numeric(sub.dat[,n.het.idx]),
                        "BB"=as.numeric(sub.dat[,n.homalt.idx]))
  HWE.mat <- HWE.mat[complete.cases(HWE.mat),]
  return(HWE.mat)
}
#Hardy-Weinberg ternary plot
plot.HWE <- function(dat,pop=NULL,title=NULL,full.legend=F,lab.cex=1){
  require(HardyWeinberg,quietly=T)
  #Gather HW p-values & colors
  HWE.mat <- makeHWEmat(dat=dat,pop=pop)
  HW.p <- HWChisqStats(X=HWE.mat,x.linked=F,pvalues=T)
  HW.cols <- rep("#4DAC26",times=length(HW.p))
  HW.cols[which(HW.p<0.05)] <- "#81F850"
  HW.cols[which(HW.p<0.05/length(HW.p))] <- "#AC26A1"
  
  #Generate HW plot frame
  par(mar=c(1,1,1,1),bty="n")
  plot(x=1.15*c(-1/sqrt(3),1/sqrt(3)),y=c(-0.15,1.15),type="n",
       xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")
  segments(x0=c(-1/sqrt(3),0,1/sqrt(3)),
           x1=c(0,1/sqrt(3),-1/sqrt(3)),
           y0=c(0,1,0),y1=c(1,0,0))
  HWTernaryPlot(X=HWE.mat,n=max(HWE.mat,na.rm=T),newframe=F,
                vbounds=F,mafbounds=F,
                region=1,vertexlab=NA,
                alpha=0.05,
                curvecols=c("#4DAC26","#81F850",NA,NA),pch=NA)
  
  #Add axes
  text(x=c(-1/sqrt(3),1/sqrt(3)),y=0,labels=c("0/0","1/1"),
       pos=1,cex=0.8,xpd=T,font=2)
  text(x=0,y=1,labels="0/1",pos=3,cex=0.8,xpd=T,font=2)
  
  #Finish HW plot
  HWTernaryPlot(X=HWE.mat,n=max(HWE.mat,na.rm=T),newframe=F,
                vbounds=F,mafbounds=F,
                region=1,vertexlab=NA,
                alpha=0.03/nrow(HWE.mat),
                curvecols=c("#4DAC26","#AC26A1",NA,NA),
                pch=21,cex=0.3,signifcolour=F,markercol=NA,
                markerbgcol=adjustcolor(HW.cols,alpha=0.25))
  segments(x0=c(-1/sqrt(3),0,1/sqrt(3)),
           x1=c(0,1/sqrt(3),-1/sqrt(3)),
           y0=c(0,1,0),y1=c(1,0,0))
  
  #Add legend
  n.pass <- length(which(HW.p>=0.05))
  print(paste("PASS: ",n.pass/length(HW.p),sep=""))
  n.nom <- length(which(HW.p<0.05 & HW.p>=0.05/nrow(HWE.mat)))
  print(paste("NOMINAL FAILS: ",n.nom/length(HW.p),sep=""))
  n.bonf <- length(which(HW.p<0.05/nrow(HWE.mat)))
  print(paste("BONFERRONI FAILS: ",n.bonf/length(HW.p),sep=""))
  legend("right",pch=19,col=c("#4DAC26","#81F850","#AC26A1"),pt.cex=1.3,
         legend=c(paste(round(100*(n.pass/nrow(HWE.mat)),0),"%",sep=""),
                  paste(round(100*(n.nom/nrow(HWE.mat)),0),"%",sep=""),
                  paste(round(100*(n.bonf/nrow(HWE.mat)),0),"%",sep="")),
         bty="n",bg=NA,cex=0.7)
  text(x=par("usr")[2],y=par("usr")[4]-(0.2*(par("usr")[4]-par("usr")[3])),pos=2,cex=0.7,
       labels=paste(title,"\n \n ",sep=""),font=2)
  text(x=par("usr")[2],y=par("usr")[4]-(0.2*(par("usr")[4]-par("usr")[3])),pos=2,cex=0.7,
       labels=paste(" \n",prettyNum(max(apply(HWE.mat,1,sum),na.rm=T),big.mark=",")," Samples\n ",sep=""))
  text(x=par("usr")[2],y=par("usr")[4]-(0.2*(par("usr")[4]-par("usr")[3])),pos=2,cex=0.7,
       labels=paste(" \n \n",prettyNum(nrow(HWE.mat),big.mark=",")," SV",sep=""))
}

#Allele frequency correlation plot between populations
plot.crossPopCorr <- function(dat,pops,popA,popB){
  #Subset to biallelic, autosomal sites with non-zero AF in both populations
  AN.A.idx <- which(colnames(dat)==paste(popA,"_AN",sep=""))
  AC.A.idx <- which(colnames(dat)==paste(popA,"_AC",sep=""))
  AF.A.idx <- which(colnames(dat)==paste(popA,"_AF",sep=""))
  AN.B.idx <- which(colnames(dat)==paste(popB,"_AN",sep=""))
  AC.B.idx <- which(colnames(dat)==paste(popB,"_AC",sep=""))
  AF.B.idx <- which(colnames(dat)==paste(popB,"_AF",sep=""))
  subdat <- dat[which(dat[,AF.A.idx]>0 & dat[,AF.B.idx]>0 & !(dat$chrom %in% c("X","Y"))),]
  rows.to.drop <- sort(unique(c(grep("MULTIALLELIC",subdat$FILTER,fixed=T),
                                grep("PESR_GT_OVERDISPERSION",subdat$FILTER,fixed=T),
                                grep("UNRESOLVED",subdat$FILTER,fixed=T))))
  if(length(rows.to.drop)>0){
    subdat <- subdat[-rows.to.drop,]
  }
  
  #Calculate p-values with chi-square test
  pvals <- sapply(1:nrow(subdat),function(i){
    AN.A <- subdat[i,AN.A.idx]
    AC.A <- subdat[i,AC.A.idx]
    AN.B <- subdat[i,AN.B.idx]
    AC.B <- subdat[i,AC.B.idx]
    return(chisq.test(data.frame(c(AN.A-AC.A,AC.A),c(AN.B-AC.B,AC.B)))$p.value)
  })
  pvals.bonf <- p.adjust(pvals,method="bonferroni")
  
  #Prepare plot
  logscale.all <- log10(as.numeric(unlist(sapply(c(0:9),function(i){(1:9)*(10^i)}))))
  logscale.major <- 0:9
  major.labels <- sapply(logscale.major,function(i){expression(paste(i^"th"))})
  par(mar=c(2.6,2.6,1.5,1.5))
  plot(x=log10(c(min(subdat[,AF.A.idx]),1)),
       y=log10(c(min(subdat[,AF.B.idx]),1)),
       type="n",xaxt="n",yaxt="n",xlab="",ylab="")
  axis(1,at=-logscale.all,labels=NA,tck=-0.015,lwd=0.7)
  axis(1,at=-logscale.major,labels=NA,tck=-0.03,lwd=1.1)
  mtext(1,text=paste(popA," AF",sep=""),line=1.25)
  axis(2,at=-logscale.all,labels=NA,tck=-0.015,lwd=0.7)
  axis(2,at=-logscale.major,labels=NA,tck=-0.03,lwd=1.1)
  mtext(2,text=paste(popB," AF",sep=""),line=1.6)
  sapply(-logscale.major,function(i){
    axis(1,at=i,labels=bquote('10'^.(i)),tick=F,line=-0.7,cex.axis=0.8)
    axis(2,at=i,labels=bquote('10'^.(i)),tick=F,line=-0.6,cex.axis=0.8,las=2)
  })
  
  #Add points
  pt.cex <- 0.3
  alpha <- 0.1
  col.A <- pops$color[which(pops$pop==popA)]
  col.B <- pops$color[which(pops$pop==popB)]
  points(x=log10(subdat[which(pvals.bonf<=0.05 & subdat[,AF.A.idx]>subdat[,AF.B.idx]),AF.A.idx]),
         y=log10(subdat[which(pvals.bonf<=0.05 & subdat[,AF.A.idx]>subdat[,AF.B.idx]),AF.B.idx]),
         pch=19,cex=pt.cex,lwd=0,
         col=adjustcolor(col.A,alpha=alpha))
  points(x=log10(subdat[which(pvals.bonf<=0.05 & subdat[,AF.B.idx]>subdat[,AF.A.idx]),AF.A.idx]),
         y=log10(subdat[which(pvals.bonf<=0.05 & subdat[,AF.B.idx]>subdat[,AF.A.idx]),AF.B.idx]),
         pch=19,cex=pt.cex,lwd=0,
         col=adjustcolor(col.B,alpha=alpha))
  points(x=log10(subdat[which(pvals.bonf>0.05),AF.A.idx]),
         y=log10(subdat[which(pvals.bonf>0.05),AF.B.idx]),
         pch=19,cex=pt.cex,lwd=0,
         col=adjustcolor("gray50",alpha=alpha))
  
  #Add stats
  # abline(lm(subdat[,AF.B.idx] ~ subdat[,AF.A.idx]),lty=2,col="gray50")
  frac.A <- length(which(pvals.bonf<=0.05 & subdat[,AF.A.idx]>subdat[,AF.B.idx]))/nrow(subdat)
  text(x=par("usr")[2]+(0.03*(par("usr")[2]-par("usr")[1])),
       y=par("usr")[3]+(0.05*(par("usr")[4]-par("usr")[3])),
       pos=2,cex=0.8,col=col.A,
       labels=paste(round(100*frac.A,digits=1),"%",sep=""))
  frac.B <- length(which(pvals.bonf<=0.05 & subdat[,AF.B.idx]>subdat[,AF.A.idx]))/nrow(subdat)
  text(x=par("usr")[1]-(0.03*(par("usr")[2]-par("usr")[1])),
       y=par("usr")[4]-(0.05*(par("usr")[4]-par("usr")[3])),
       pos=4,cex=0.8,col=col.B,
       labels=paste(round(100*frac.B,digits=1),"%",sep=""))
  AB.cor <- format(round(cor(subdat[,AF.A.idx],subdat[,AF.B.idx])^2,3),nsmall=3)
  mtext(3,line=0.1,text=bquote(italic(R)^2 == .(AB.cor)))
}

###Perform Watterson estimator analysis of mutation rate
#Individual function to get mutation rate estimate for a single pop & svtype
getMu <- function(dat,pop=NULL,svtype=NULL,Ne=10000){
  #Format variables
  if(!is.null(pop)){
    pop <- paste(pop,"_",sep="")
  }
  #Get number of chromosomes assessed
  pop.AN.idx <- which(colnames(dat)==paste(pop,"N_BI_GENOS",sep=""))
  n <- 2*max(dat[,pop.AN.idx],na.rm=T)
  #Get number of autosomal biallelic segregating sites
  pop.AF.idx <- which(colnames(dat)==paste(pop,"AF",sep=""))
  if(!is.null(svtype)){
    K <- length(which(dat[,pop.AF.idx]>0 & dat$SVTYPE==svtype & !(dat$chrom %in% c("X","Y")) & dat$SVTYPE != "MCNV"))
  }else{
    K <- length(which(dat[,pop.AF.idx]>0 & !(dat$chrom %in% c("X","Y")) & dat$SVTYPE != "MCNV"))
  }
  #Get n-1th harmonic number
  harmsum <- sum(sapply(1:(n-1),function(k){1/k}))
  #Get Watterson estimator
  theta.hat <- K/harmsum
  #Solve for mutation rate
  mu <- theta.hat/(4*Ne)
  return(mu)
}
#Wrapper function to get mutation rate estimates for all svtypes & populations
getAllMus <- function(dat,Ne){
  mu.svtypes <- c("DEL","DUP","INS","INV","CPX")
  #Get mutation rates across classes
  mu.AllClasses <- sapply(mu.svtypes,function(svtype){
    getMu(dat,pop=NULL,svtype=svtype,Ne=Ne)
  })
  #Get mutation rate of all SV across populations
  mu.AllPops <- sapply(pops$pop,function(pop){
      getMu(dat,pop=pop,svtype=NULL,Ne=Ne)
  })
  #Get mutation rate by svtype & population
  mu.ClassByPop <- sapply(pops$pop,function(pop){
    sapply(mu.svtypes,function(svtype){
      getMu(dat,pop=pop,svtype=svtype,Ne=Ne)
    })
  })
  #Get mean by class across pops
  mu.PopMeanByClass <- apply(mu.ClassByPop,1,function(vals){
    lower <- as.numeric(t.test(vals)$conf.int[1])
    mean <- as.numeric(t.test(vals)$estimate)
    upper <- as.numeric(t.test(vals)$conf.int[2])
    return(c(lower,mean,upper))
  })
  mu.PopMeanByClass <- cbind(apply(mu.PopMeanByClass,1,sum),mu.PopMeanByClass)
  colnames(mu.PopMeanByClass)[1] <- "ALL"
  rownames(mu.PopMeanByClass) <- c("CI.lowerBound","mean","CI.upperBound")
  return(list("muByClass"=mu.AllClasses,
              "muByPop"=mu.AllPops,
              "muByClassAndPop"=mu.ClassByPop,
              "mu.Means"=mu.PopMeanByClass))
}
#Plot mutation rates
plotMus <- function(dat,Ne){
  mu.dat <- getAllMus(dat,Ne)
  plot.dat <- mu.dat$mu.Means
  mu.svtypes <- c("DEL","DUP","INS","INV","CPX")
  mu.cols <- c("gray30",sapply(mu.svtypes,function(svtype){
    svtypes$color[which(svtypes$svtype==svtype)]
  }))
  werling.rates <- c(166,80,47,37,0,1)/1038
  #Prep plot area
  par(mar=c(1.5,3.5,0.5,0.5),bty="n")
  plot(x=c(-0.025,ncol(plot.dat))+0.2,y=c(0,max(plot.dat)),type="n",
       xaxt="n",yaxt="n",xlab="",ylab="")
  axis(1,at=c(-4,10),tck=0,labels=NA,lwd=1.5)
  axis(2,at=c(-1,1),tck=0,labels=NA,lwd=1.5)
  axis(1,at=(1:ncol(plot.dat))-0.5,tick=F,line=-0.8,
       labels=c("All SV",mu.svtypes))
  axis(2,labels=NA)
  axis(2,tick=F,las=2,cex.axis=0.8,line=-0.4)
  mtext(2,line=2,text="Mutation Rate (SVs per generation)")
  #Add points & CIs
  text.buf <- 0.04*(par("usr")[4]-par("usr")[3])
  sapply(1:ncol(plot.dat),function(i){
    segments(x0=i-0.5,x1=i-0.5,
             y0=plot.dat[1,i],y1=plot.dat[3,i],
             lend="round",lwd=2,col=mu.cols[i])
    points(x=i-0.5,y=plot.dat[2,i],pch=19,cex=1.5,col=mu.cols[i])
    par(xpd=T)
    # text(x=i-c(0.6,0.45,0.6),y=plot.dat[,i]+c(-text.buf,0,text.buf),
    #      cex=c(0.6,0.7,0.6),labels=format(round(plot.dat[,i],4),nsmall=4),pos=4,
    #      col=mu.cols[i],font=c(3,2,3))
    text(x=i-0.45,y=plot.dat[2,i],
         cex=1,labels=format(round(plot.dat[2,i],3),nsmall=3),pos=4,
         col=mu.cols[i],font=2)
    par(xpd=F)
    points(x=i-0.5,y=werling.rates[i],pch=23,lwd=1.5)
  })
  legend("topright",border=NA,bty="n",
         legend=c(as.expression(bquote(mu ~ "from Watterson" ~ hat(theta[italic(W)]) ~ "in gnomAD, " ~ N[e] == .(prettyNum(Ne,big.mark=",")))),
                  expression("Rate of validated" ~ italic("de novo") ~ "SV from 519 quartets")),
         pch=c(19,23),pt.lwd=1.5,lwd=c(2,NA),lty=c(1,NA),cex=0.8,pt.cex=c(1.25,1))
}


###################
###COMPLEX SV TABLE
###################
#Master function to make complex SV table
masterCPXtable <- function(dat,
                           cpx.types=c("delINV","INVdel","delINVdel",
                                       "dupINV","INVdup","dupINVdup",
                                       "delINVdup","dupINVdel",
                                       "piDUP_FR","piDUP_RF",
                                       "dDUP","dDUP_iDEL","INS_iDEL")){
  #Set parameters
  filler.color <- "gray50"
  panels <- c("CPX Class","Involved SV\nSignatures","Rearrangement Schematic","Variants",
              "SV per\nGenome","Median\nSize","Size\nDistribution","Mean AF")
  n.panels <- length(panels)
  panel.widths=c(0.75,0.55,2,0.5,0.5,0.5,1,0.5)
  panel.heights=c(1,rep(1,times=length(cpx.types)+1),0.75)

  #Empty header text box with single value
  headerTextBox <- function(text,cex=1){
    plot(x=c(0,1),y=c(0,1),type="n",
         xaxt="n",xlab="",xaxs="i",
         yaxt="n",ylab="",yaxs="i")
    text(x=0.5,y=par("usr")[3],pos=3,labels=text,cex=cex,font=2)
    axis(1,at=c(par("usr")[1],par("usr")[2]),tck=0,labels=NA)
  }
  
  #Empty text box with single value
  textBox <- function(text,cex=1,pos=2){
    plot(x=c(0,1),y=c(0,1),type="n",
         xaxt="n",xlab="",xaxs="i",
         yaxt="n",ylab="",yaxs="i")
    text(x=par("usr")[2],y=0.5,pos=pos,labels=text,cex=cex)
  }
  
  #Blank box with nothing
  blankBox <- function(){
    plot(x=c(0,1),y=c(0,1),type="n",
         xaxt="n",xlab="",xaxs="i",
         yaxt="n",ylab="",yaxs="i")
  }
  
  #Boxes for SV types involved in CPX class
  cpx.signature <- function(cpx.type){
    plot(x=c(0,4),y=c(0,1),type="n",
         xaxt="n",xlab="",xaxs="i",
         yaxt="n",ylab="",yaxs="i")
    rect(xleft=(0:3)+0.2,xright=(0:3)+0.8,
         ybottom=0.2,ytop=0.8,col="gray95",border="gray90")
    if(length(grep("DEL",cpx.type,ignore.case=T))>0){
      rect(xleft=0.2,xright=0.8,ybottom=0.2,ytop=0.8,
           col=svtypes$color[which(svtypes$svtype=="DEL")])
    }
    if(length(grep("DUP",cpx.type,ignore.case=T))>0){
      rect(xleft=1.2,xright=1.8,ybottom=0.2,ytop=0.8,
           col=svtypes$color[which(svtypes$svtype=="DUP")])
    }
    if(length(c(grep("INS",cpx.type,ignore.case=T),
                grep("dDUP",cpx.type,ignore.case=T)))>0){
      rect(xleft=2.2,xright=2.8,ybottom=0.2,ytop=0.8,
           col=svtypes$color[which(svtypes$svtype=="INS")])
    }
    if(length(c(grep("INV",cpx.type,ignore.case=T),
                grep("piDUP",cpx.type,ignore.case=T)))>0){
      rect(xleft=3.2,xright=3.8,ybottom=0.2,ytop=0.8,
           col=svtypes$color[which(svtypes$svtype=="INV")])
    }
  }
  
  #Simple barplot
  cpx.simpleBar <- function(value,max.value=1){
    plot(x=c(0,max.value),y=c(0,1),type="n",
         xaxt="n",xlab="",xaxs="i",
         yaxt="n",ylab="",yaxs="i")
    rect(xleft=0,xright=value,
         ybottom=0.25,ytop=0.75,
         col=filler.color)
  }
  
  #Simple size density plot
  cpx.SVLEN.dens <- function(dat,cpx.type){
    require(zoo,quietly=T)
    if(cpx.type=="ALL"){
      sizes <- log10(dat$SVLEN[which(dat$SVTYPE=="CPX")])
    }else{
      sizes <- log10(dat$SVLEN[which(dat$SVTYPE=="CPX" & dat$CPX_TYPE==cpx.type)])
    }
    h <- hist(sizes,breaks=seq(0,10,0.05),plot=F)$counts
    h <- rollapply(h,11,mean,partial=T)
    h <- h/max(h,na.rm=T)
    out.df <- data.frame("min.size"=10^seq(0,9.95,0.05),
                         "dens"=h)
    plot(x=log10(c(50,10000000)),y=c(0,1),type="n",
         xaxt="n",xlab="",xaxs="i",
         yaxt="n",ylab="",yaxs="i")
    polygon(x=log10(c(out.df$min.size,rev(out.df$min.size))),
            y=c(out.df$dens,rep(0,times=nrow(out.df))),
            col=filler.color)
    axis(1,at=c(par("usr")[1],par("usr")[2]),tck=0,labels=NA)
    abline(v=median(sizes),col=svtypes$color[which(svtypes$svtype=="CPX")],lwd=2)
  }

  #Pct axis
  pctAxis <- function(max.value=1){
    plot(x=c(0,max.value),y=c(0,1),type="n",
         xaxt="n",xlab="",xaxs="i",
         yaxt="n",ylab="",yaxs="i")
    axis(3,at=axTicks(3),labels=NA,tck=0.125)
    axis(3,at=axTicks(3),tick=F,line=-2.25,cex.axis=0.7,
         labels=paste(round(100*axTicks(3),1),"%",sep=""))
    # mtext(3,text="Pct.",line=-2.5,cex=0.7)
  }
  
  #Size axis
  sizeAxis <- function(){
    plot(x=log10(c(50,10000000)),y=c(0,1),type="n",
         xaxt="n",xlab="",xaxs="i",
         yaxt="n",ylab="",yaxs="i")
    axis(3,at=as.numeric(unlist(sapply(1:7,function(x){log10((1:9)*(10^x))}))),labels=NA,tck=0.1,lwd=0.7)
    axis(3,at=2:7,labels=NA,tck=0.2)
    labels <- c("10bp","100bp","1kb","10kb","100kb","1Mb","10Mb")
    sapply(2:7,function(i){
      axis(3,at=i,tick=F,labels=labels[i],line=-2.25,cex.axis=0.7)
    })
    # mtext(3,text="SV Size",line=-2.5,cex=0.7)
  }
  
  #Prep figure layout
  layout(matrix(1:(n.panels*(length(cpx.types)+3)),byrow=T,nrow=length(cpx.types)+3),
         widths=panel.widths,heights=panel.heights)
  
  #Headers
  par(bty="n",mar=c(0.1,1,0.1,1),xpd=T)
  sapply(panels,function(title){
    headerTextBox(title)
  })
  
  #Add top row for summary of all CPX
  par(bty="n",mar=c(0.5,0.5,0.5,0.5),xpd=F)
  textBox(text="All CPX SV")
  blankBox()
  blankBox()
  textBox(text=prettyNum(length(which(dat$SVTYPE=="CPX")),big.mark=","))
  textBox(text=prettyNum(round(sum(c(dat$N_HET[which(dat$SVTYPE=="CPX")],
                                     dat$N_HOMALT[which(dat$SVTYPE=="CPX")]))/max(dat$N_BI_GENOS[which(dat$SVTYPE=="CPX")]),
                               1),big.mark=","))
  textBox(text=paste(format(round(median(dat$SVLEN[which(dat$SVTYPE=="CPX")])/1000,1),nsmall=1),"kb",sep=""))
  cpx.SVLEN.dens(dat,cpx.type="ALL")
  textBox(text=paste(format(round(100*mean(dat$AF[which(dat$SVTYPE=="CPX")],na.rm=T),2),nsmall=2),"%",sep=""))
  
  #Iterate over cpx types - one row per type
  par(bty="n",mar=c(0.5,0.5,0.5,0.5),xpd=F)
  sapply(1:length(cpx.types),function(i){
    cpx.type <- cpx.types[i]
    #CPX class (text)
    textBox(text=cpx.type)
    #Signatures
    cpx.signature(cpx.type)
    #Placeholder for rearrangement schematic (done in illustrator)
    blankBox()
    #Count of variants
    textBox(text=prettyNum(length(which(dat$SVTYPE=="CPX" & dat$CPX_TYPE==cpx.type)),big.mark=","))
    # #Pct of all CPX
    # cpx.simpleBar(value=length(which(dat$SVTYPE=="CPX" & dat$CPX_TYPE==cpx.type))/length(which(dat$SVTYPE=="CPX")),
    #               max.value=max(table(dat$CPX_TYPE[which(dat$SVTYPE=="CPX")])/length(which(dat$SVTYPE=="CPX"))))
    #Mean # of variants per genome
    textBox(text=prettyNum(round(sum(c(dat$N_HET[which(dat$SVTYPE=="CPX" & dat$CPX_TYPE==cpx.type)],
                       dat$N_HOMALT[which(dat$SVTYPE=="CPX" & dat$CPX_TYPE==cpx.type)]))/max(dat$N_BI_GENOS[which(dat$SVTYPE=="CPX")]),
                       1),big.mark=","))
    #Median size (text)
    textBox(text=paste(format(round(median(dat$SVLEN[which(dat$SVTYPE=="CPX" & dat$CPX_TYPE==cpx.type)])/1000,1),nsmall=1),"kb",sep=""))
    #Size distribution
    cpx.SVLEN.dens(dat,cpx.type)
    #Mean AF (text)
    textBox(text=paste(format(round(100*mean(dat$AF[which(dat$SVTYPE=="CPX" & dat$CPX_TYPE==cpx.type)],na.rm=T),2),nsmall=2),"%",sep=""))
  })
  
  #Axes, if necessary
  par(bty="n",mar=c(0.1,0.5,0.5,0.1),xpd=T)
  sapply(rep(0,4),function(i){blankBox()})
  # pctAxis(max.value=max(table(dat$CPX_TYPE[which(dat$SVTYPE=="CPX")])/length(which(dat$SVTYPE=="CPX"))))
  sapply(rep(0,2),function(i){blankBox()})
  sizeAxis()
}

################################
###MERGED SIMPLE + COMPLEX TABLE
################################
#Master function to make merged canonical + complex SV table
masterSVtableFigure <- function(dat.wrelateds,dat.all,dat,svtypes,med.sitesPerSample){
  #Set parameters
  filler.color <- "gray50"
  simple.rows <- c("All SV","DEL","DUP","MCNV","INS","INV","CTX","CPX")
  complex.rows <- c("delINV\nINVdel","delINVdel",
                    "dupINV\nINVdup","dupINVdup",
                    "delINVdup\ndupINVdel",
                    "piDUP (FR)\npiDUP (RF)",
                    "dDUP","dDUP-iDEL\nINS-iDEL")
  cpx.groups <- list(c("delINV","INVdel"),"delINVdel",
                     c("dupINV","INVdup"),"dupINVdup",
                     c("delINVdup","dupINVdel"),
                     c("piDUP_FR","piDUP_RF"),
                     "dDUP",c("dDUP_iDEL","INS_iDEL"))
  #Order complex groups by size
  new.cpx.order <- order(-unlist(lapply(cpx.groups,function(classes){
    median(dat$SVLEN[which(dat$SVTYPE=="CPX" & dat$CPX_TYPE %in% classes)])
  })))
  complex.rows <- complex.rows[new.cpx.order]
  cpx.groups <- cpx.groups[new.cpx.order]
  sv.rows <- c(simple.rows,complex.rows,"BND")
  panels <- c("","SV Class","Abbrev.","Mutational\nSignature","","Ref. Allele\nStructure","Alt. Allele\nStructure(s)",
              "Resolved\nVariants","SV per\nGenome","SV Size","Proportion\nSingletons","")
  n.panels <- length(panels)
  panel.widths=c(0.2,0.7,0.5,0.4,0.2,0.65,1,0.5,0.5,0.6,0.5,0.2)
  panel.heights=c(0.9,rep(1,times=length(sv.rows)),0.75)
  
  #Empty header text box with single value
  headerTextBox <- function(text,cex=1){
    plot(x=c(0,1),y=c(0,1),type="n",
         xaxt="n",xlab="",xaxs="i",
         yaxt="n",ylab="",yaxs="i")
    if(text!=""){
      text(x=0.5,y=par("usr")[3],pos=3,labels=text,cex=cex,font=2)
      axis(1,at=c(par("usr")[1],par("usr")[2]),tck=0,labels=NA)
    }
  }
  
  #Jewel corresponding to SV class
  jewel <- function(svtype,svtypes){
    plot(x=c(0,1),y=c(0,1),type="n",
         xaxt="n",xlab="",xaxs="i",
         yaxt="n",ylab="",yaxs="i")
    if(svtype=="ALL"){
      color <- "gray30"
    }else if(svtype %in% c("CTX","BND")){
      color <- svtypes$color[which(svtypes$svtype=="OTH")]
    }else{
      color <- svtypes$color[which(svtypes$svtype==svtype)]
    }
    points(x=0.5,y=0.5,pch=21,bg=color,cex=2.5)
  }
  
  #Empty text box with single value
  textBox <- function(text,cex=1,pos=2,x=1,font=1,color="black"){
    plot(x=c(0,1),y=c(0,1),type="n",
         xaxt="n",xlab="",xaxs="i",
         yaxt="n",ylab="",yaxs="i")
    text(x=x,y=0.5,pos=pos,labels=text,cex=cex,font=font,col=color,xpd=T)
  }
  
  #Blank box with nothing
  blankBox <- function(){
    plot(x=c(0,1),y=c(0,1),type="n",
         xaxt="n",xlab="",xaxs="i",
         yaxt="n",ylab="",yaxs="i")
  }
  
  #Get text corresponding to each sv class
  getClassText <- function(svtype){
    if(svtype=="ALL"){
      text <- "All SV"
    }else if(svtype=="DEL"){
      text <- "Deletion"
    }else if(svtype=="DUP"){
      text <- "Duplication"
    }else if(svtype=="MCNV"){
      text <- "Multiallelic CNV"
    }else if(svtype=="INS"){
      text <- "Insertion"
    }else if(svtype=="INV"){
      text <- "Inversion"
    }else if(svtype=="CTX"){
      text <- "Reciprocal\nTranslocation"
    }else if(svtype=="BND"){
      text <- "Breakend\n(Unresolved)"
    }else if(svtype=="CPX"){
      text <- "All Complex SV"
    }
    return(text)
  }
  
  #Get text corresponding to each complex group
  getCPXSubclassText <- function(cpx.subclasses){
    if("delINV" %in% cpx.subclasses){
      text <- "Deletion-Flanked\nInversion"
    }else if("delINVdel" %in% cpx.subclasses){
      text <- "Paired-Deletion\nInversion"
    }else if("dupINV" %in% cpx.subclasses){
      text <- "Dup.-Flanked\nInversion"
    }else if("dupINVdup" %in% cpx.subclasses){
      text <- "Paired-Dup.\nInversion"
    }else if("delINVdup" %in% cpx.subclasses){
      text <- "Paired-Del./Dup.\nInversion"
    }else if("piDUP_FR" %in% cpx.subclasses){
      text <- "Palindromic\nInverted Dup."
    }else if("dDUP" %in% cpx.subclasses){
      text <- "Dispersed\nDuplication"
    }else if("dDUP_iDEL" %in% cpx.subclasses){
      text <- "Insertion with\nIns. Site Del."
    }else{
      text <- "N/A"
    }
    return(text)
  }
  
  #Boxes for SV types involved in CPX class
  signatures <- function(svtype="CPX",cpx.subclasses=NULL){
    cpx.type <- paste(cpx.subclasses,collapse="")
    plot(x=c(0,4),y=c(0,1),type="n",
         xaxt="n",xlab="",xaxs="i",
         yaxt="n",ylab="",yaxs="i")
    rect(xleft=(0:3)+0.2,xright=(0:3)+0.8,
         ybottom=0.2,ytop=0.8,col="gray95",border="gray90")
    if(length(grep("DEL",svtype,ignore.case=T))>0 
       | svtype=="MCNV"
       | length(grep("DEL",cpx.type,ignore.case=T))>0){
      rect(xleft=0.2,xright=0.8,ybottom=0.2,ytop=0.8,
           col=svtypes$color[which(svtypes$svtype=="DEL")])
    }
    if(length(grep("DUP",svtype,ignore.case=T))>0 
       | svtype=="MCNV"
       | length(grep("DUP",cpx.type,ignore.case=T))>0){
      rect(xleft=1.2,xright=1.8,ybottom=0.2,ytop=0.8,
           col=svtypes$color[which(svtypes$svtype=="DUP")])
    }
    if(length(grep("INS",svtype,ignore.case=T))>0 
       | length(c(grep("INS",cpx.type,ignore.case=T),
                  grep("dDUP",cpx.type,ignore.case=T)))>0){
      rect(xleft=2.2,xright=2.8,ybottom=0.2,ytop=0.8,
           col=svtypes$color[which(svtypes$svtype=="INS")])
    }
    if(length(grep("INV",svtype,ignore.case=T))>0 
       | length(c(grep("INV",cpx.type,ignore.case=T),
                  grep("piDUP",cpx.type,ignore.case=T)))>0){
      rect(xleft=3.2,xright=3.8,ybottom=0.2,ytop=0.8,
           col=svtypes$color[which(svtypes$svtype=="INV")])
    }
  }
  
  #Get count of sites, pct pass, and count of final variants
  getCounts <- function(dat.wrelateds,dat,svtype,cpx.subclasses=NULL){
    if(svtype=="ALL"){
      sites <- nrow(dat.wrelateds)
      # sites.pass <- length(which(dat.wrelateds$FILTER %in% c("PASS","MULTIALLELIC")))
      final.sv <- nrow(dat)
    }else if(svtype=="CPX" & !is.null(cpx.subclasses)){
      sites <- length(which(dat.wrelateds$SVTYPE==svtype & dat.wrelateds$CPX_TYPE %in% cpx.subclasses))
      final.sv <- length(which(dat$SVTYPE==svtype & dat$CPX_TYPE %in% cpx.subclasses))
    }else{
      sites <- length(which(dat.wrelateds$SVTYPE==svtype))
      # sites.pass <- length(which(dat.wrelateds$SVTYPE==svtype & dat.wrelateds$FILTER %in% c("PASS","MULTIALLELIC")))
      final.sv <- length(which(dat$SVTYPE==svtype))
    }
    # pct.pass <- sites.pass/sites
    pct.pass <- final.sv/sites
    return(c(sites,pct.pass,final.sv))
  }
  
  #Get count of median sites per sample
  getSitesPerSample <- function(med.sitesPerSample,svtype){
    if(svtype %in% c("CTX","BND")){
      x <- 0
    }else if(svtype=="MCNV"){
      x <- sum(med.sitesPerSample[1,grep("MCNV_",colnames(med.sitesPerSample),fixed=T)])
    }else{
      x <- sum(med.sitesPerSample[1,which(colnames(med.sitesPerSample)==svtype)])
    }
    return(x)
  }
  
  #Simple barplot
  cpx.simpleBar <- function(value,max.value=1){
    plot(x=c(0,max.value),y=c(0,1),type="n",
         xaxt="n",xlab="",xaxs="i",
         yaxt="n",ylab="",yaxs="i")
    rect(xleft=0,xright=value,
         ybottom=0.25,ytop=0.75,
         col=filler.color)
  }
  
  #Simple size density plot
  SVLEN.dens <- function(dat,svtype,cpx.subclasses=NULL,font=1){
    require(zoo,quietly=T)
    if(svtype=="ALL"){
      sizes <- log10(dat$SVLEN)
      line.col <- "gray30"
    }else if(svtype=="CPX"){
      line.col <- svtypes$color[which(svtypes$svtype=="CPX")]
      if(is.null(cpx.subclasses)){
        sizes <- log10(dat$SVLEN[which(dat$SVTYPE=="CPX")])
      }else{
        sizes <- log10(dat$SVLEN[which(dat$SVTYPE=="CPX" & dat$CPX_TYPE %in% cpx.subclasses)])
      }
    }else{
      line.col <- svtypes$color[which(svtypes$svtype==svtype)]
      sizes <- log10(dat$SVLEN[which(dat$SVTYPE==svtype)])
    }
    # h <- hist(sizes,breaks=seq(0,10,0.05),plot=F)$counts
    # h <- rollapply(h,11,mean,partial=T)
    # h <- h/max(h,na.rm=T)
    # out.df <- data.frame("min.size"=10^seq(0,9.95,0.05),
    #                      "dens"=h)
    plot(x=log10(c(50,10000000)),y=c(0,1),type="n",
         xaxt="n",xlab="",xaxs="i",
         yaxt="n",ylab="",yaxs="i")
    rect(xleft=log10(50),xright=log10(10000000),ybottom=0.1,ytop=0.6,col="white",border="gray85")
    segments(x0=0:10,x1=0:10,y0=0.1,y1=0.6,col="gray90")
    # segments(x0=log10(c(50,10000000)),x1=log10(c(50,10000000)),y0=0.1,y1=0.6,col="black",lwd=1.5)
    # abline(v=0:10,col="gray90")
    # polygon(x=log10(c(out.df$min.size,rev(out.df$min.size))),
    #         y=c(out.df$dens,rep(0,times=nrow(out.df))),
    #         col=filler.color)
    # axis(1,at=c(par("usr")[1],par("usr")[2]),tck=0,labels=NA)
    # abline(v=median(sizes,na.rm=T),col=line.col,lwd=3,lend="round")
    require(vioplot,quietly=T)
    viosizes <- as.numeric(sizes[which(!is.na(sizes) & !is.infinite(sizes) & !is.nan(sizes))])
    if(length(viosizes)>0){
      vioplot(viosizes,horizontal=T,add=T,col=line.col,at=0.35,drawRect=F,h=0.25,wex=0.5)
      segments(x0=median(sizes,na.rm=T),
               x1=median(sizes,na.rm=T),
               y0=0.15,y1=0.55,col="black",lwd=2,lend="round")
      text(x=median(sizes,na.rm=T),y=0.5,pos=3,cex=1,font=font,xpd=T,
           labels=paste(format(round(median((10^sizes)/1000,na.rm=T),1),n.small=1),"kb",sep=""))
    }
  }
  
  #Calculate proportion of singletons for a single SV class
  calc.fracSingletons.singleClass <- function(dat,svtype,cpx.subclasses=NULL,boot.n=100,conf=0.95){
    if(svtype=="ALL"){
      ACs <- dat$AC
    }else if(svtype=="CPX" & !is.null(cpx.subclasses)){
      ACs <- dat[which(dat$SVTYPE==svtype & dat$CPX_TYPE %in% cpx.subclasses),]$AC
    }else{
      ACs <- dat[which(dat$SVTYPE==svtype),]$AC
    }
    ACs <- ACs[which(!is.na(ACs) & ACs>0)]
    helper.getFracSingle <- function(ACs,indices){length(which(ACs[indices]==1))/length(ACs[indices])}
    point.est <- helper.getFracSingle(ACs,indices=1:length(ACs))
    calc.ci <- function(ACs,n,conf){
      set.seed(0)
      boot.obj <- boot(data=ACs,statistic=helper.getFracSingle,R=n)
      ci <- boot.ci(boot.obj,conf=conf,type="basic")$basic[4:5]
      return(ci)
    }
    ci <- calc.ci(ACs,n=boot.n,conf=conf)
    return(c(point.est,ci))
  }
  
  #Plot point estimate & 95% CI for proportion singletons
  propSingles <- function(dat,svtype,cpx.subclasses=NULL,font=1){
    if(svtype=="ALL"){
      point.col <- "gray30"
    }else if(svtype=="CPX"){
      point.col <- svtypes$color[which(svtypes$svtype=="CPX")]
    }else if(svtype %in% c("CTX","BND")){
      point.col <- svtypes$color[which(svtypes$svtype=="OTH")]
    }else{
      point.col <- svtypes$color[which(svtypes$svtype==svtype)]
    }
    s.dat <- calc.fracSingletons.singleClass(dat,svtype=svtype,cpx.subclasses=cpx.subclasses)
    plot(x=c(0.3,0.9),y=c(0,1),type="n",
         xaxt="n",xlab="",xaxs="i",
         yaxt="n",ylab="",yaxs="i")
    # abline(v=seq(0,1,0.25),col="gray90")
    # abline(v=c(0,1))
    rect(xleft=0.3,xright=0.9,ybottom=0.1,ytop=0.6,col="white",border="gray85")
    segments(x0=seq(0.3,0.9,0.1),x1=seq(0.3,0.9,0.1),
             y0=0.1,y1=0.6,col="gray90")
    # segments(x0=c(0.3,0.9),x1=c(0.3,0.9),
    #          y0=0.1,y1=0.6,col="black",lwd=1.5)
    segments(x0=s.dat[2],x1=s.dat[3],y0=0.35,y1=0.35,lend="round",lwd=1.5)
    rect(xleft=s.dat[1]-0.015,xright=s.dat[1]+0.015,
         ybottom=0.15,ytop=0.55,col=point.col,lwd=1)
    # segments(x0=s.dat[2],x1=s.dat[3],y0=0.35,y1=0.35,lend="round",lwd=1.5)
    # points(x=s.dat[1],y=0.35,pch=21,bg=point.col,cex=1.5)
    text(x=s.dat[1],y=0.5,pos=3,cex=1,font=font,xpd=T,
         labels=paste(format(round(100*s.dat[1],1),n.small=1),"%",sep=""))
  }
  
  #Pct axis
  pctAxis <- function(range=c(0.3,0.9),ticks.at=seq(0.3,0.9,0.1),labels.at=seq(0.3,0.9,0.2)){
    plot(x=range,y=c(0,1),type="n",
         xaxt="n",xlab="",xaxs="i",
         yaxt="n",ylab="",yaxs="i")
    axis(3,at=ticks.at,labels=NA,tck=0.15)
    sapply(labels.at,function(x){
      axis(3,at=x,tick=F,cex.axis=0.75,las=2,hadj=1,line=-1.6,
           labels=paste(round(100*x,1),"%",sep=""))
    })
    # mtext(3,text="Pct.",line=-2.5,cex=0.7)
  }
  
  #Size axis
  sizeAxis <- function(){
    plot(x=log10(c(50,10000000)),y=c(0,1),type="n",
         xaxt="n",xlab="",xaxs="i",
         yaxt="n",ylab="",yaxs="i")
    axis(3,at=as.numeric(unlist(sapply(1:7,function(x){log10((1:9)*(10^x))}))),labels=NA,tck=0.075,lwd=0.5)
    axis(3,at=1:8,labels=NA,tck=0.15)
    labels <- c("10bp","100bp","1kb","10kb","100kb","1Mb","10Mb")
    sapply(seq(2,7,by=1),function(i){
      axis(3,at=i,tick=F,labels=labels[i],cex.axis=0.75,las=2,hadj=1,line=-1.6)
    })
    # mtext(3,text="SV Size",line=-2.5,cex=0.7)
  }
  
  #Prep figure layout
  layout(matrix(1:(n.panels*(length(sv.rows)+2)),byrow=T,nrow=length(sv.rows)+2),
         widths=panel.widths,heights=panel.heights)
  
  #Headers
  par(bty="n",mar=c(0.1,1,0.1,1),xpd=T)
  sapply(panels,function(title){
    headerTextBox(title)
  })
  
  #Add top row for summary of all SV
  par(bty="n",mar=c(0.4,0.5,0.4,0.5),xpd=F)
  jewel("ALL",svtypes)
  textBox(text=getClassText("ALL"),font=2)
  textBox(text="ALL",font=2)
  textBox("Varies",color="gray50",pos=NULL,x=0.5,font=3)
  blankBox()
  textBox("Varies",color="gray50",pos=NULL,x=0.5,font=3)
  textBox("Varies",color="gray50",pos=NULL,x=0.5,font=3)
  # textBox(text=prettyNum(getCounts(dat=dat,dat.wrelateds=dat.wrelateds,svtype="ALL")[1],big.mark=","),font=2)
  # textBox(paste(format(round(100*getCounts(dat=dat,dat.wrelateds=dat.wrelateds,svtype="ALL")[2],1),n.small=1),"%",sep=""),font=2)
  textBox(text=paste(prettyNum(getCounts(dat=dat,dat.wrelateds=dat.wrelateds,svtype="ALL")[3],big.mark=","),"*",sep=""),font=2)
  # textBox(text=prettyNum(round(sum(c(dat$N_HET,dat$N_HOMALT),na.rm=T)/max(dat$N_BI_GENOS,na.rm=T),
  #                              1),big.mark=","))
  textBox(prettyNum(getSitesPerSample(med.sitesPerSample,svtype="ALL"),big.mark=","),font=2)
  SVLEN.dens(dat,svtype="ALL",cpx.subclasses=NULL,font=2)
  # textBox(text=paste(format(round(median(dat$SVLEN)/1000,1),nsmall=1),"kb",sep=""),font=2)
  propSingles(dat,svtype="ALL",font=2)
  blankBox()
  
  #Iterate over canonical classes
  sapply(c("DEL","DUP","MCNV","INS","INV","CTX"),function(svtype){
    #Jewel
    jewel(svtype,svtypes)
    #SV class (text)
    textBox(getClassText(svtype))
    #Abbreviation (text)
    textBox(text=svtype)
    #Signatures
    if(svtype %in% c("CTX","BND")){
      textBox("N/A",color="gray50",pos=NULL,x=0.5)
    }else{
      signatures(svtype)
    }
    blankBox()
    #Placeholder for rearrangement schematic (done in illustrator)
    blankBox()
    blankBox()
    # textBox(text=prettyNum(getCounts(dat=dat,dat.wrelateds=dat.wrelateds,svtype=svtype)[1],big.mark=","))
    # textBox(paste(format(round(100*getCounts(dat=dat,dat.wrelateds=dat.wrelateds,svtype=svtype)[2],1),n.small=1),"%",sep=""))
    textBox(text=prettyNum(getCounts(dat=dat,dat.wrelateds=dat.wrelateds,svtype=svtype)[3],big.mark=","))
    #Median # of variants per genome
    textBox(prettyNum(getSitesPerSample(med.sitesPerSample,svtype=svtype),big.mark=","))
    #Median size (text) & size distribution
    if(svtype %in% c("CTX","BND")){
      textBox(text="N/A",color="gray50",pos=NULL,x=0.5)
      # textBox(text="N/A")
    }else{
      SVLEN.dens(dat,svtype=svtype,cpx.subclasses=NULL)
      # textBox(text=paste(format(round(median(dat$SVLEN[which(dat$SVTYPE==svtype)])/1000,1),nsmall=1),"kb",sep=""))
    }
    if(svtype %in% c("MCNV","BND")){
      textBox("N/A",color="gray50",pos=NULL,x=0.5)
    }else{
      propSingles(dat,svtype=svtype)
    }
    blankBox()
  })
  
  #One summary row for all CPX SV
  jewel("CPX",svtypes)
  textBox(getClassText("CPX"),font=2)
  textBox(text="CPX",font=2)
  textBox(text="Varies",color="gray50",pos=NULL,x=0.5,font=3)
  blankBox()
  textBox(text="Varies",color="gray50",pos=NULL,x=0.5,font=3)
  textBox(text="Varies",color="gray50",pos=NULL,x=0.5,font=3)
  # textBox(text=prettyNum(getCounts(dat=dat,dat.wrelateds=dat.wrelateds,svtype="CPX")[1],big.mark=","),font=2)
  # textBox(paste(format(round(100*getCounts(dat=dat,dat.wrelateds=dat.wrelateds,svtype="CPX")[2],1),n.small=1),"%",sep=""),font=2)
  textBox(text=prettyNum(getCounts(dat=dat,dat.wrelateds=dat.wrelateds,svtype="CPX")[3],big.mark=","),font=2)
  textBox(prettyNum(getSitesPerSample(med.sitesPerSample,svtype="CPX"),big.mark=","),font=2)
  SVLEN.dens(dat,svtype="CPX",cpx.subclasses=NULL,font=2)
  # textBox(text=paste(format(round(median(dat$SVLEN[which(dat$SVTYPE=="CPX")])/1000,1),nsmall=1),"kb",sep=""))
  propSingles(dat,svtype="CPX",font=2)
  blankBox()
  
  #Iterate over cpx types - one row per group of types
  sapply(1:length(cpx.groups),function(i){
    #Get parameters
    cpx.group <- complex.rows[i]
    cpx.subclasses <- cpx.groups[[i]]
    #Jewel
    jewel(svtype="CPX",svtypes)
    #Category name
    textBox(text=getCPXSubclassText(cpx.subclasses),cex=1)
    #Abbrev.
    if(length(cpx.subclasses)>1){
      textBox(text=cpx.group,cex=1)
    }else{
      textBox(text=cpx.group)
    }
    #Mutation signature
    signatures(svtype="CPX",cpx.subclasses=cpx.subclasses)
    blankBox()
    #Spacer for alt allele schematics
    sapply(rep(0,2),function(i){blankBox()})
    #Counts of sites
    # textBox(text=prettyNum(getCounts(dat=dat,dat.wrelateds=dat.wrelateds,
    #                                  svtype="CPX",cpx.subclasses=cpx.subclasses)[1],big.mark=","))
    # textBox(paste(format(round(100*getCounts(dat=dat,dat.wrelateds=dat.wrelateds,
    #                                          svtype="CPX",cpx.subclasses=cpx.subclasses)[2],1),n.small=1),"%",sep=""))
    textBox(text=prettyNum(getCounts(dat=dat,dat.wrelateds=dat.wrelateds,
                                     svtype="CPX",cpx.subclasses=cpx.subclasses)[3],big.mark=","))
    #TODO: collect per-sample counts of each complex type
    textBox("TBD")
    #Size distribution
    SVLEN.dens(dat,svtype="CPX",cpx.subclasses=cpx.subclasses)
    # textBox(text=paste(format(round(median(dat$SVLEN[which(dat$SVTYPE=="CPX" & dat$CPX_TYPE %in% cpx.subclasses)])/1000,1),nsmall=1),"kb",sep=""))
    propSingles(dat,svtype="CPX",cpx.subclasses=cpx.subclasses)
    #Final spacer
    blankBox()
  })
  
  #Add final row for BNDs  
  #Jewel
  jewel(svtype="BND",svtypes)
  #SV class (text)
  textBox(getClassText("BND"))
  #Abbreviation (text)
  textBox(text="BND")
  #Signatures
  textBox("N/A",color="gray50",pos=NULL,x=0.5)
  blankBox()
  #Placeholder for rearrangement schematic (done in illustrator)
  blankBox()
  blankBox()
  # textBox(text=prettyNum(getCounts(dat=dat,dat.wrelateds=dat.wrelateds,svtype="BND")[1],big.mark=","))
  # textBox(paste(format(round(100*getCounts(dat=dat,dat.wrelateds=dat.wrelateds,svtype="BND")[2],1),n.small=1),"%",sep=""))
  textBox(text=paste(prettyNum(getCounts(dat=dat.all,dat.wrelateds=dat.wrelateds,svtype="BND")[3],big.mark=","),"*",sep=""))
  #Median # of variants per genome
  # textBox(prettyNum(getSitesPerSample(med.sitesPerSample,svtype="BND"),big.mark=","))
  textBox("N/A",color="gray50",font=3)
  #Median size (text) & size distribution
  textBox(text="N/A",color="gray50",pos=NULL,x=0.5)
  textBox("N/A",color="gray50",pos=NULL,x=0.5)
  blankBox()
  
  #Axes, if necessary
  par(bty="n",mar=c(0.1,0.5,0.1,0.5),xpd=T)
  sapply(rep(0,9),function(i){blankBox()})
  # pctAxis(max.value=max(table(dat$CPX_TYPE[which(dat$SVTYPE=="CPX")])/length(which(dat$SVTYPE=="CPX"))))
  # sapply(rep(0,1),function(i){blankBox()})
  sizeAxis()
  pctAxis()
}

#Doubleton analysis table
doubleton.analysis <- function(dat){
  #Identify doubletons
  doubletons <- dat[which(dat$N_HET==2 & dat$AC==2),]
  nhet.pop.idxs <- which(colnames(doubletons) %in% paste(pops$pop[which(pops$pop!="OTH")],"_N_HET",sep=""))
  doubletons.samepop <- doubletons[which(apply(doubletons[,nhet.pop.idxs],1,max)==2),]
  doubletons.crosspop <- doubletons[which(apply(doubletons[,nhet.pop.idxs],1,max)==1),]
  #Iterate over svtypes and compute fraction of samepop vs crosspop doubletons
  sapply(c("DEL","DUP","INS","INV",""),function(svtype){
    all <- length(which(doubletons$SVTYPE==svtype))
    samepop <- length(which(doubletons.samepop$SVTYPE==svtype))
    crosspop <- length(which(doubletons.crosspop$SVTYPE==svtype))
    frac.cross <- crosspop/all
    return(c(all,samepop,crosspop,frac.cross))
  })
}

#Get table of gross chromosomal abnormalities
gather.bca.table <- function(dat,max.AF=0.01){
  res <- as.data.frame(t(sapply(c("DEL","DUP","INV","CTX","CPX"),function(svtype){
    idxs <- which(dat$SVTYPE==svtype 
                  & !(dat$chrom %in% c("X","Y"))
                  & dat$AF<max.AF 
                  & (dat$SVLEN>=1000000 | dat$SVTYPE=="CTX"))
    # if(svtype=="CPX"){
    #   idxs <- unique(c(idxs,
    #                    which(dat$SVTYPE==svtype 
    #                          & !(dat$chrom %in% c("X","Y"))
    #                          & dat$AF<max.AF 
    #                          & dat$CPX_TYPE=="CCR")))
    # }
    variants <- length(idxs)
    carriers <- sum(dat$N_HET[idxs])
    #One CCR is actually a complex translocation, and should be counted
    if(svtype=="CPX"){carriers <- carriers+1; variants <- variants+1}
    # if(svtype=="CTX"){carriers <- 17; variants <- 15}
    denom <- max(dat$N_BI_GENOS[idxs],na.rm=T)
    rate <- carriers/denom
    binom.ci <- binom.test(carriers,denom)$conf.int
    c(variants,carriers,rate,binom.ci)
  })))
  colnames(res) <- c("variants","carriers","mean","lower","upper")
  res <- rbind(apply(res,2,sum),res)
  rownames(res)[1] <- "ALL"
  return(res)
}

#Plot basic metadata for gross chromosomal abnormalities
plot.chrom.abnormalities <- function(dat,max.AF=0.01,
                                     category.labels=c("Deletion",
                                                       "Duplication",
                                                       "Inversion",
                                                       "Reciprocal\nTranslocation",
                                                       "Complex SV")){
  #Get plotting data
  plot.dat <- gather.bca.table(dat=dat,max.AF=max.AF)
  plot.dat <- plot.dat[-1,]
  xmax <- max(plot.dat[,-c(1:2)])
  point.colors <- as.character(sapply(rownames(plot.dat),function(svtype){
    if(svtype=="CTX"){
      svtype <- "OTH"
    }
    svtypes$color[which(svtypes$svtype==svtype)]
  }))
  #Prep plot area
  layout(matrix(c(1:4),nrow=1,byrow=T),
         widths=c(4,1.5,3,5))
  #Plot class label
  par(mar=c(0.3,0.3,3,0.3),bty="n",xpd=T)
  plot(x=c(0,1),y=c(0,-nrow(plot.dat)),type="n",
       xlab="",ylab="",xaxt="n",yaxt="n")
  text(x=1,y=(-1:-nrow(plot.dat))+0.5,pos=2,
       labels=category.labels)
  mtext(3,line=1.2,text="Rare SV")
  mtext(3,line=0,text=(bquote("" >= "1Mb")))
  axis(3,at=c(0,1),labels=NA,tck=0)
  #Plot number of variants
  par(xpd=T)
  plot(x=c(0,1),y=c(0,-nrow(plot.dat)),type="n",
       xlab="",ylab="",xaxt="n",yaxt="n")
  text(x=0.5,y=(-1:-nrow(plot.dat))+0.5,cex=1,
       labels=plot.dat$variants)
  mtext(3,line=0,text="SV")
  axis(3,at=c(0,1),labels=NA,tck=0)
  #Plot number of carriers
  plot(x=c(0,1),y=c(0,-nrow(plot.dat)),type="n",
       xlab="",ylab="",xaxt="n",yaxt="n")
  text(x=0.5,y=(-1:-nrow(plot.dat))+0.5,cex=0.95,
       labels=paste(plot.dat$carriers,"/",prettyNum(max(dat$N_BI_GENOS,na.rm=T),big.mark=",")))
  mtext(3,line=0,text="Samples",cex=0.95)
  axis(3,at=c(0,1),labels=NA,tck=0)
  #Plot points & CIs
  par(mar=c(0.3,0.4,3,0.3),bty="n",xpd=F)
  plot(x=c(0,xmax),y=c(0,-nrow(plot.dat)),type="n",
       xlab="",ylab="",xaxt="n",yaxt="n")
  abline(v=axTicks(3),col="gray90")
  abline(v=0,col="gray50")
  segments(x0=plot.dat$lower,x1=plot.dat$upper,
           y0=(-1:-nrow(plot.dat))+0.5,
           y1=(-1:-nrow(plot.dat))+0.5,
           lend="round",lwd=2)
  points(x=plot.dat$mean,
         y=(-1:-nrow(plot.dat))+0.5,
         pch=21,cex=1.5,bg=point.colors)
  par(xpd=T)
  text(x=plot.dat$mean,y=(-1:-nrow(plot.dat))+0.5,
       pos=3,labels=paste(format(round(100*plot.dat$mean,2),nsmall=2),"%",sep=""))
  par(xpd=F)
  axis(3,at=axTicks(3),labels=NA,tck=-0.025)
  sapply(axTicks(3)[seq(1,length(axTicks(3)),2)],function(x){
    axis(3,at=x,tick=F,line=-0.8,cex.axis=0.9,
         labels=paste(round(100*x,1),"%",sep=""))
  })
  mtext(3,line=1.1,text="Samples (Pct.)")
  # #Plot schematic label
  # plot(x=c(0,1),y=c(0,-nrow(plot.dat)),type="n",
  #      xlab="",ylab="",xaxt="n",yaxt="n")
  # mtext(3,line=-0.9,text="Schematic")
}


################
###RSCRIPT BLOCK
################
require(optparse,quietly=T)
require(boot,quietly=T)
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
args <- parse_args(OptionParser(usage="%prog VCF2BED PRIORSTUDIES POPASSIGNMENTS OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 7){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
vcf2bed.in <- args$args[1]
vcf2bed.wrelateds.in <- args$args[2]
priorstudies.in <- args$args[3]
popassignments.in <- args$args[4]
sd_sr_coverage.in <- args$args[5]
medians.perSample.in <- args$args[6]
OUTDIR <- args$args[7]
prefix <- opts$prefix
svtypes.file <- opts$svtypes
pops.file <- opts$populations

# #Dev parameters (local)
# vcf2bed.in <- "~/scratch/gnomAD_v2_SV_MASTER.vcf2bed.bed.gz"
# vcf2bed.wrelateds.in <- "~/scratch/gnomAD_v2_SV_MASTER_wRelateds.vcf2bed.bed.gz"
# priorstudies.in <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/prior_SV_study_stats.txt"
# popassignments.in <- "~/scratch/gnomAD_v2_SV_MASTER.updated_pop_assignments.txt"
# sd_sr_coverage.in <- "~/scratch/gnomAD_v2_SV_MASTER.variant_SD_SR_coverage.txt.gz"
# medians.perSample.in <- "~/scratch/gnomAD_v2_SV_MASTER.median_sites_per_sample_by_population.txt"
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
dat.wrelateds <- read.vcf2bed(vcf2bed.wrelateds.in)
dat.all <- read.vcf2bed(vcf2bed.in)
dat <- dat.all[which(dat.all$FILTER=="PASS" | dat.all$FILTER=="MULTIALLELIC"),]
priors <- read.table(priorstudies.in,header=T,sep="\t")
priors <- priors[-which(priors$Study=="Mills"),]
pop.assignments <- read.table(popassignments.in,header=F)
colnames(pop.assignments) <- c("sample","pop")
pop.assignments$pop[which(pop.assignments$pop==".")] <- "OTH"
sd_sr_cov <- read.table(sd_sr_coverage.in,header=F,sep="\t")
med.sitesPerSample <- read.table(medians.perSample.in,header=T,sep="\t")

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

###Plots comparison of total N and total sites discovered vs other studies
cat("NOW STARTING BARPLOTS VS PRIOR STUDIES\n")
pdf(paste(OUTDIR,"/",prefix,".sample_size_vs_prior_studies.pdf",sep=""),
    height=1.5,width=1.75)
par(mar=c(2,3.9,0.1,1.2))
plot.sampleSize(priors,pops,pop.assignments)
dev.off()
pdf(paste(OUTDIR,"/",prefix,".sv_sites_vs_prior_studies.pdf",sep=""),
    height=1.5,width=1.9)
par(mar=c(2,3.9,0.1,1.65))
plot.SVcountsByStudy(priors,dat.wrelateds,svtypes)
dev.off()

###Plot simple bars of counts
cat("NOW STARTING COUNT OF SV BY TYPE\n")
pdf(paste(OUTDIR,"/",prefix,".site_counts_by_type.pdf",sep=""),
    height=2.75,width=2.25)
plot.totalCounts(dat=dat,svtypes=svtypes,thousandG=F)
dev.off()
pdf(paste(OUTDIR,"/",prefix,".site_counts_by_type.all_sites_unfiltered.pdf",sep=""),
    height=2.75,width=2.25)
plot.totalCounts(dat=dat.all,svtypes=svtypes,thousandG=F)
dev.off()

###Build table of counts, sizes, singletons, etc
callset.summary.table <- build.table(dat.wrelateds,dat.all,dat)
write.table(callset.summary.table,
            paste(OUTDIR,"/",prefix,".callset_summary_table.txt",sep=""),
            col.names=F,row.names=F,sep="\t",quote=F)

###Plot size distribution
cat("NOW STARTING SIZE DISTRIBUTIONS\n")
pdf(paste(OUTDIR,"/",prefix,".size_distribution_by_type.pdf",sep=""),
    height=2,width=2.9)
plot.sizes(dat=dat,svtypes=svtypes)
dev.off()

###Plot frequency distribution by SVTYPE
cat("NOW STARTING FREQUENCY DISTRIBUTIONS\n")
pdf(paste(OUTDIR,"/",prefix,".site_frequency_distributions_bySVtype.pdf",sep=""),
    height=2,width=2.3)
plot.freqByType(dat=dat)
dev.off()
#Plot three-panel of SV with and without predicted coding effects
pdf(paste(OUTDIR,"/",prefix,".site_frequency_distributions_bySVtype.by_genomic_context.pdf",sep=""),
    height=1.75,width=5.5)
par(mfrow=c(1,3))
plot.freqByType(dat=dat[which(!is.na(dat$PROTEIN_CODING__LOF)
                              | !is.na(dat$PROTEIN_CODING__COPY_GAIN)
                              | !is.na(dat$PROTEIN_CODING__DUP_LOF)),],
                ymin=0.38,axlabel.cex=0.8)
plot.freqByType(dat=dat[which(is.na(dat$PROTEIN_CODING__LOF)
                              & is.na(dat$PROTEIN_CODING__COPY_GAIN)
                              & is.na(dat$PROTEIN_CODING__DUP_LOF)
                              & (!is.na(dat$PROTEIN_CODING__INTRONIC))
                                 | !is.na(dat$PROTEIN_CODING__DUP_PARTIAL)),],
                ymin=0.38,axlabel.cex=0.8)
plot.freqByType(dat=dat[which(dat$PROTEIN_CODING__INTERGENIC=="True"),],
                ymin=0.38,axlabel.cex=0.8)
dev.off()

###Plot proportion of singletons as a function of SV size
pdf(paste(OUTDIR,"/",prefix,".prop_singletons_vs_size.pdf",sep=""),
    height=2.3,width=2.25)
plot.fracSingletons.bySize(dat=dat,sd_sr_cov=sd_sr_cov,
                           svtype="ALL",max.cov=0.3,include.ins=T,
                           d=c(-1,1000,10000,100000,1000000,1000000000),
                           d.labels=c("<1kb","1-10kb","10-100kb","100kb-1Mb",">1Mb"),
                           mar=c(4,2.75,0.5,1),ylab.cex=1)
mtext(1,line=3,text="SV Size",cex=1.2)
dev.off()
#Plot grid for all sv classes
pdf(paste(OUTDIR,"/",prefix,".prop_singletons_vs_size.by_svtype.pdf",sep=""),
    height=4,width=5)
par(mfrow=c(2,3))
sapply(c("ALL","DEL","DUP","INS","INV","CPX"),function(svtype){
  if(svtype!="ALL"){
    color <- svtypes$color[which(svtypes$svtype==svtype)]
    title <- svtype
  }else{
    color <- "black"
    title <- "All SV"
  }
  plot.fracSingletons.bySize(dat=dat,sd_sr_cov=sd_sr_cov,
                             svtype=svtype,max.cov=0.3,include.ins=T,
                             d=c(-1,1000,10000,100000,1000000,1000000000),
                             d.labels=c("<1kb","1-10kb","10-100kb","100kb-1Mb",">1Mb"),
                             color=color,title=title,mar=c(3.5,2.75,2.5,1),ylab.cex=0.7,
                             ylims=NULL)
})
dev.off()

###Plot HWE ternary graphs
cat("NOW STARTING HARDY-WEINBERG EQUILIBRIUM CALCULATIONS\n")
png(paste(OUTDIR,"/",prefix,".HWE_all_samples.png",sep=""),
    height=1000,width=1000,res=400)
plot.HWE(dat=dat,pop=NULL,title="All Samples")
dev.off()
if(!is.null(pops.file)){
  sapply(1:nrow(pops),function(i){
    pop <- pops$pop[i]
    if(pop=="AFR"){
      pop.name <- "African"
    }else{
      pop.name <- pops$name[i]
    }
    print(pop)
    png(paste(OUTDIR,"/",prefix,".HardyWeinberg_ternary_plot.",pop,"_samples.png",sep=""),
        height=1000,width=1000,res=400)
    plot.HWE(dat=dat,pop=pop,title=pop.name)
    dev.off()
  })
}
# #HWE by SVTYPE
# if(!is.null(svtypes.file)){
#   sapply(c("DEL","DUP","INS","INV","CPX"),function(svtype){
#     png(paste(OUTDIR,"/",prefix,".HardyWeinberg_ternary_plot.all_samples.",svtype,".png",sep=""),
#         height=1000,width=1000,res=400)
#     plot.HWE(dat=dat[which(dat$SVTYPE==svtype),],pop=NULL,title="All Samples")
#     dev.off()
#   })
# }
# #HWE for LoF variants only
# png(paste(OUTDIR,"/",prefix,".HWE_all_samples.LOF.png",sep=""),
#     height=1000,width=1000,res=400)
# plot.HWE(dat=dat[which(!is.na(dat$PROTEIN_CODING__LOF)),],
#          pop=NULL,title="All Samples")
# dev.off()
# if(!is.null(pops.file)){
#   sapply(1:nrow(pops),function(i){
#     pop <- pops$pop[i]
#     pop.name <- pops$name[i]
#     png(paste(OUTDIR,"/",prefix,".HardyWeinberg_ternary_plot.",pop,"_samples.LOF.png",sep=""),
#         height=1000,width=1000,res=400)
#     plot.HWE(dat=dat[which(!is.na(dat$PROTEIN_CODING__LOF)),],
#              pop=pop,title=pop.name)
#     dev.off()
#   })
# }
# #HWE for intronic variants only
# png(paste(OUTDIR,"/",prefix,".HWE_all_samples.intergenic.png",sep=""),
#     height=1000,width=1000,res=400)
# plot.HWE(dat=dat[which(!is.na(dat$PROTEIN_CODING__INTERGENIC)),],
#          pop=NULL,title="All Samples")
# dev.off()
# if(!is.null(pops.file)){
#   sapply(1:nrow(pops),function(i){
#     pop <- pops$pop[i]
#     pop.name <- pops$name[i]
#     png(paste(OUTDIR,"/",prefix,".HardyWeinberg_ternary_plot.",pop,"_samples.intergenic.png",sep=""),
#         height=1000,width=1000,res=400)
#     plot.HWE(dat=dat[which(!is.na(dat$PROTEIN_CODING__INTERGENIC)),],
#              pop=pop,title=pop.name)
#     dev.off()
#   })
# }

###Plot cross-population allele frequency correlations
cat("NOW STARTING CROSS-POPULATION ALLELE FREQUENCY CORRELATIONS\n")
sapply(1:nrow(pops),function(a){
  sapply(2:nrow(pops),function(b){
    if(b>a){
      popA <- pops$pop[a]
      popB <- pops$pop[b]
      if(popA != "OTH" & popB != "OTH"){
        png(paste(OUTDIR,"/",prefix,".frequency_correlation.",popA,"_vs_",popB,".png",sep=""),
            height=1000,width=1000,res=400)
        plot.crossPopCorr(dat,pops,popA=popA,popB=popB)
        dev.off()
      }
    }
  })
})

###PLOT COMPLEX SV BREAKDOWN TABLE
cat("NOW STARTING COMPLEX SV SUBTYPES TABLE\n")
pdf(paste(OUTDIR,"/",prefix,".complex_sv_subtypes.pdf",sep=""),
    height=5.25,width=8)
masterCPXtable(dat)
dev.off()

###PLOT MERGED SIMPLE + COMPLEX SV BREAKDOWN TABLE
cat("NOW STARTING SIMPLE + COMPLEX SV SUBTYPES TABLE\n")
pdf(paste(OUTDIR,"/",prefix,".simple_and_complex_sv_classes.pdf",sep=""),
    height=7.6,width=8.75)
masterSVtableFigure(dat.wrelateds=dat.wrelateds,dat.all=dat.all,dat=dat,
                    svtypes=svtypes,med.sitesPerSample=med.sitesPerSample)
dev.off()

###PLOT MUTATION RATE ESTIMATES
cat("NOW STARTING MUTATION RATE ANALYSIS\n")
pdf(paste(OUTDIR,"/",prefix,".mutation_rate_estimates.pdf",sep=""),
    height=2.9,width=5)
plotMus(dat,Ne=10000)
dev.off()

###Chromosomal abnormalities
pdf(paste(OUTDIR,"/",prefix,".chrom_abnormalities.pdf",sep=""),
    height=2,width=3)
plot.chrom.abnormalities(dat=dat)
dev.off()

###Wilcoxon test of allele frequencies between SV classes
format(wilcox.test(dat$AF[which(dat$SVTYPE=="DEL")],
            dat$AF[which(dat$SVTYPE!="DEL" & dat$SVTYPE!="CPX")],
            alternative="less")$p.value,
       scientific=T)
format(wilcox.test(dat$AF[which(dat$SVTYPE=="CPX")],
            dat$AF[which(dat$SVTYPE!="CPX")],
            alternative="less")$p.value,
       scientific=T)
