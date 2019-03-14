#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Generate plots corresponding to all VCF downsampling analyses


###Set master parameters
options(stringsAsFactors=F,scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Import table of total sites vs sample size
import.sites <- function(totalsites.in,
                         extrapolate=c(10000,100000,1000000,1000000000),
                         include.bnd=F){
  dat <- read.table(totalsites.in,header=T,sep="\t")
  if(include.bnd==F){
    dat <- dat[,which(colnames(dat)!="BND")]
  }
  dat$ALL <- apply(dat[,-1],1,sum,na.rm=T)
  #Compute medians
  n.points <- sort(unique(dat$N_samples))
  medians <- do.call("rbind", lapply(n.points,function(n){
    apply(dat[which(dat$N_samples==n),-1],2,median,na.rm=T)
  }))
  medians <- data.frame("N_samples"=n.points,medians)
  #Fit linear model to each SVTYPE and extrapolate
  #Note: fit on sample sizes >100 to avoid noise at low downsamples
  extrapolations <- sapply(2:ncol(dat),function(i){
    if(any(dat[,i]>0)){
      p.dat <- data.frame("n"=log10(medians[which(medians$N_samples>=100),1]),
                          "sv"=log10(medians[which(medians$N_samples>=100),i]))
      fit <- lm(sv ~ n, data=p.dat, na.action="na.omit")
      predict(fit,data.frame("n"=log10(extrapolate)))
    }else{
      rep(NA,times=length(extrapolate))
    }
  })
  colnames(extrapolations) <- colnames(dat)[-1]
  extrapolations <- data.frame("N_samples"=extrapolate,
                               extrapolations)
  return(list("medians"=medians,
              "raw.data"=dat,
              "extrapolations"=extrapolations))
}
#Import table of singletons per sample by ancestry
import.singletons <- function(singletons.in,seed.table.in,pops,include.bnd=F){
  dat <- read.table(singletons.in,header=T,sep="\t")
  if(include.bnd==F){
    dat <- dat[,-grep("BND",colnames(dat))]
  }
  dat$N <- as.integer(read.table(seed.table.in,header=F)[,1])
  pop.dat <- lapply(pops$pop,function(pop){
    col.idxs <- grep(pop,colnames(dat),fixed=T)
    if(length(col.idxs)>0){
      subdat <- cbind("N"=dat$N,dat[,col.idxs])
      subdat <- subdat[which(subdat[,1]>0 & subdat[,2]>0),]
      subdat <- subdat[order(subdat[,1]),]
      colnames(subdat) <- gsub(paste(pop,".",sep=""),"",colnames(subdat),fixed=T)
      subdat$ALL <- apply(subdat[,-c(1:2)],1,sum,na.rm=T)
      return(subdat)
    }
  })
  names(pop.dat) <- pops$pop
  return(pop.dat)
}
#Import and format gene data
import.genes <- function(svpergene.in,SNVdata.in){
  a <- read.table(svpergene.in,header=T,sep="\t",comment.char="")
  b <- read.table(SNVdata.in,header=T,sep="\t",comment.char="")
  b$oe_pct <- round(99*rank(b$ptv_oe)/nrow(b),0)+1
  gene.dat <- merge(a,b,by="gene",sort=F)
  gene.dat$constrained <- 0
  # gene.dat$constrained[which(gene.dat$pli>0.9)] <- 1
  gene.dat$constrained[which(gene.dat$oe_pct<=15)] <- 1
  gene.dat$unconstrained <- 0
  # gene.dat$unconstrained[which(gene.dat$pli<0.1)] <- 1
  gene.dat$unconstrained[which(gene.dat$oe_pct>85)] <- 1
  return(gene.dat)
}
#Helper function to get mean and bootstrapped 95% CI for a vector of values
boot.mean.ci <- function(vals,boot.n=100,conf=0.95){
  helper.mean <- function(vals,indices){mean(vals[indices])}
  point.est <- helper.mean(vals,indices=1:length(vals))
  calc.ci <- function(vals,n,conf){
    set.seed(0)
    boot.obj <- boot(data=vals,statistic=helper.mean,R=n)
    ci <- boot.ci(boot.obj,conf=conf,type="basic")$basic[4:5]
    return(ci)
  }
  ci <- calc.ci(vals,n=boot.n,conf=conf)
  return(c(point.est,ci))
}
#Gather average count of functional sv per gene vs sample size for a single functional class
get.gene.trend <- function(gene.dat,func,seed.table.in,invert.uncons=F){
  seeds <- read.table(seed.table.in,header=F,sep="\t",comment.char="")
  Ns <- sort(unique(seeds[,1]))
  all.genes <- as.data.frame(t(sapply(Ns,function(n){
    c(n,boot.mean.ci(gene.dat[,which(colnames(gene.dat)==paste(func,n,sep="."))]))
  })))
  colnames(all.genes) <- c("N","mean","lower","upper")
  constrained.genes <- as.data.frame(t(sapply(Ns,function(n){
    c(n,boot.mean.ci(gene.dat[which(gene.dat$constrained==1),
                              which(colnames(gene.dat)==paste(func,n,sep="."))]))
  })))
  colnames(constrained.genes) <- c("N","mean","lower","upper")
  if(invert.uncons==T){
    unconstrained.genes <- as.data.frame(t(sapply(Ns,function(n){
      c(n,boot.mean.ci(gene.dat[which(gene.dat$constrained!=1),
                                which(colnames(gene.dat)==paste(func,n,sep="."))]))
    }))) 
  }else{
    unconstrained.genes <- as.data.frame(t(sapply(Ns,function(n){
      c(n,boot.mean.ci(gene.dat[which(gene.dat$unconstrained==1),
                                which(colnames(gene.dat)==paste(func,n,sep="."))]))
    })))
  }
  colnames(unconstrained.genes) <- c("N","mean","lower","upper")
  return(list("all"=all.genes,
              "constrained"=constrained.genes,
              "unconstrained"=unconstrained.genes))
}

#Gather pct of genes with at least one functional SV vs sample size for a single functional class
get.pct.genes.bySampSize <- function(gene.dat,func,seed.table.in){
  seeds <- read.table(seed.table.in,header=F,sep="\t",comment.char="")
  Ns <- sort(unique(seeds[,1]))
  all.genes <- as.data.frame(t(sapply(Ns,function(n){
    vals <- gene.dat[,which(colnames(gene.dat)==paste(func,n,sep="."))]
    vals[which(vals>1)] <- 1
    c(n,boot.mean.ci(vals))
  })))
  colnames(all.genes) <- c("N","mean","lower","upper")
  constrained.genes <- as.data.frame(t(sapply(Ns,function(n){
    vals <- gene.dat[which(gene.dat$constrained==1),
                     which(colnames(gene.dat)==paste(func,n,sep="."))]
    vals[which(vals>1)] <- 1
    c(n,boot.mean.ci(vals))
  })))
  colnames(constrained.genes) <- c("N","mean","lower","upper")
  unconstrained.genes <- as.data.frame(t(sapply(Ns,function(n){
    vals <- gene.dat[which(gene.dat$unconstrained==1),
                     which(colnames(gene.dat)==paste(func,n,sep="."))]
    vals[which(vals>1)] <- 1
    c(n,boot.mean.ci(vals))
  })))
  colnames(unconstrained.genes) <- c("N","mean","lower","upper")
  return(list("all"=all.genes,
              "constrained"=constrained.genes,
              "unconstrained"=unconstrained.genes))
}


#####################
###PLOTTING FUNCTIONS
#####################
#Plot total SV sites vs sample size
plot.totalSites <- function(sv.sites,svtypes,empirical.counts,max.N){
  #Prep plot area
  par(bty="n",mar=c(2.25,2.75,0.5,1))
  xlims <- log10(c(10,max.N))
  ylims <- c(1,max(log10(max(c(empirical.counts,max(sv.sites$raw.data[,-1]))))))
  plot(x=xlims,y=ylims,type="n",
       xaxt="n",xlab="",yaxt="n",ylab="")
  
  #Add points per SVTYPE
  sapply(c("ALL",svtypes$svtype),function(svtype){
    idx <- which(colnames(sv.sites$medians)==svtype)
    if(length(idx)>0){
      if(svtype=="ALL"){
        color <- "gray20"
      }else{
        color <- svtypes$color[which(svtypes$svtype==svtype)]
      }
      emp <- empirical.counts[which(names(empirical.counts)==svtype)]
      #Smoothed spline of actual downsampled VCFs
      #Points of actual downsampled VCF medians
      # points(x=log10(sv.sites$medians[,1]),
      #        y=log10(sv.sites$medians[,idx]),
      #        cex=0.5,col=color,pch=19,lwd=0.75)
      lines(smooth.spline(x=log10(c(sv.sites$medians[,1],max.N)),
                          y=log10(c(sv.sites$medians[,idx],emp))),
            col=color,lwd=3)
      points(x=log10(max.N),y=log10(emp),
             lwd=1,cex=0.9,col="black",pch=23,bg=color,xpd=T)
      
      # #Projection
      # segments(x0=log10(sv.sites$extrapolations[1,1]),
      #          x1=log10(sv.sites$extrapolations[nrow(sv.sites$extrapolations),1]),
      #          y0=sv.sites$extrapolations[1,idx],
      #          y1=sv.sites$extrapolations[nrow(sv.sites$extrapolations),idx],
      #          lty=3,col=color,lwd=2)
      # points(x=log10(sv.sites$extrapolations[which(sv.sites$extrapolations[,1] %in% extrapolate),1]),
      #        y=sv.sites$extrapolations[which(sv.sites$extrapolations[,1] %in% extrapolate),idx],
      #        lwd=1.5,cex=0.7,col=color,pch=23,bg="white")
    }
  })
  
  #Add axes
  log.minor <- log10(as.numeric(sapply(0:10,function(i){(1:8)*(10^i)})))
  log.mid <- log10(as.numeric(sapply(0:10,function(i){c(1,5)*(10^i)})))
  log.major <- 1:9
  axis(1,at=log.minor,labels=NA,tck=-0.0175,lwd=0.9)
  # axis(1,at=log.mid,labels=NA,tck=-0.02)
  axis(1,at=log.major,labels=NA,tck=-0.035,lwd=1.1)
  sapply(0:6,function(i){
    axis(1,at=i,tick=F,line=-0.9,cex.axis=0.8,
         labels=c("1","10","100","1k","10k","100k","1M")[i+1])
  })
  mtext(1,text="Sample Size",line=1)
  axis(2,at=log.minor,labels=NA,tck=-0.0175,lwd=0.9)
  # axis(2,at=log.mid,labels=NA,tck=-0.02)
  axis(2,at=log.major,labels=NA,tck=-0.035,lwd=1.1)
  axis(2,at=1:6,tick=F,line=-0.7,cex.axis=0.8,las=2,
       labels=c("10","100","1k","10k","100k","1M"))  
  mtext(2,text="SV Discovered",line=1.8)
}
#Plot singletons per sample vs sample size
plot.singletons <- function(singletons,pops,med.sitesPerSample,max.N,ymax=NULL){
  #Prep plot area
  par(bty="n",mar=c(2.4,2.65,0.3,0.3))
  xlims <- c(0.5,log10(1000*ceiling((max(do.call("rbind",singletons)[,2]))/1000)))
  if(is.null(ymax)){
    ymax <- log10(max(do.call("rbind",singletons)[,-1]))
  }
  ylims <- c(0,log10(ymax))
  plot(x=c(min(xlims),max(xlims)+0.5),y=ylims,type="n",
       xaxt="n",xlab="",yaxt="n",ylab="")
  rect(xleft=c(0.5,max(xlims)+0.25),
       xright=c(0.75,max(xlims)+0.5),
       ybottom=par("usr")[3],ytop=par("usr")[4],
       border="gray90",col=NA,lwd=2)
  
  #Add axes
  log.minor <- log10(as.numeric(sapply(0:10,function(i){(1:8)*(10^i)})))
  log.mid <- log10(as.numeric(sapply(0:10,function(i){c(1,5)*(10^i)})))
  log.mid <- log.mid
  log.major <- 0:9
  log.major <- log.major
  axis(1,at=log.minor[which(log.minor<=max(xlims) & log.minor>=1)],labels=NA,tck=-0.0175,lwd=0.9)
  axis(1,at=log.minor[which(log.minor<=max(xlims) & log.minor>=1)],labels=NA,tck=0,lwd=1.1)
  axis(1,at=log.mid[which(log.mid<=max(xlims) & log.mid>=1)],labels=NA,tck=-0.025)
  axis(1,at=log.major[which(log.major<=max(xlims) & log.major>=1)],labels=NA,tck=-0.035,lwd=1.1)
  sapply(log.major[which(log.major<=max(xlims))],function(i){
      axis(1,at=i,tick=F,line=-0.7,cex.axis=0.8,
           labels=c("1","10","100","1k","10k","100k","1M")[i+1])
  })
  axis(1,at=0.625,tck=-0.035,lwd=1.1,labels=NA)
  axis(1,at=0.625,tick=F,line=-0.7,cex.axis=0.8,labels=1)
  axis(1,at=max(xlims)+0.375,tick=F,line=-1.65,labels="Pop.\nMax",cex.axis=0.7,padj=1)
  mtext(1,text="Population Reference Size",line=1.3)
  axis(2,at=log.minor,labels=NA,tck=-0.0175,lwd=0.9)
  # axis(2,at=log.mid,labels=NA,tck=-0.02)
  axis(2,at=log.major,labels=NA,tck=-0.035,lwd=1.1)
  axis(2,at=0:6,tick=F,line=-0.5,cex.axis=0.8,las=2,
       labels=c("1","10","100","1k","10k","100k","1M"))  
  mtext(2,text="Singleton SV per Sample",line=1.7)
  
  #Add points per population
  sapply(rev(pops$pop),function(pop){
    p.dat <- as.data.frame(singletons[which(names(singletons)==pop)])
    p.dat <- p.dat[which(p.dat[,2]>=10),]
    colnames(p.dat) <- gsub(paste(pop,".",sep=""),"",colnames(p.dat),fixed=T)
    color <- pops$color[which(pops$pop==pop)]
    #Gather median of each downsample point
    p.dat.meds <- t(sapply(unique(p.dat$N),function(N){
      c(median(p.dat$samples[which(p.dat$N==N)]),
        median(p.dat$ALL[which(p.dat$N==N)]))
    }))
    x <- log10(p.dat.meds[,1])
    y <- log10(p.dat.meds[,2])
    # x.med <- unique(sort(p.dat$samples))
    # y.med <- sapply(x.med,function(n){
    #   median(p.dat$ALL[which(p.dat$samples==n)],na.rm=T)
    # })
    # x.med <- log10(x.med)
    # y.med <- log10(y.med)
    # lines(smooth.spline(x=x,y=y),
    #       col=color,lwd=3)
    points(x=x,y=y,col=color,lwd=3,type="l")
    # abline(lm(y ~ x),lwd=3,col=color)
    #Points of actual downsampled VCFs
    # points(x=x.med,y=y.med,cex=0.65,col=color,pch=19,lwd=0.75)
    #Points of n=1 and popmax
    segments(x0=c(0.5,max(xlims)+0.25),
             x1=c(0.75,max(xlims)+0.5),
             col=color,lwd=3,lend="butt",
             y0=c(log10(med.sitesPerSample[which(med.sitesPerSample[,1]==paste("all",pop,sep=".")),]$ALL),
                  log10(med.sitesPerSample[which(med.sitesPerSample[,1]==paste("singleton",pop,sep=".")),]$ALL)),
             y1=c(log10(med.sitesPerSample[which(med.sitesPerSample[,1]==paste("all",pop,sep=".")),]$ALL),
                  log10(med.sitesPerSample[which(med.sitesPerSample[,1]==paste("singleton",pop,sep=".")),]$ALL)))
  })
}
#Plot curves for mean # of functional SV per gene
plot.svpergene <- function(gene.dat,func,seed.table.in,color,max.N,
                           xlabel=NULL,ylabel=NULL,mar=c(3,3,0.3,0.3),legend=F){
  #Gather plotting data
  plot.dat <- get.gene.trend(gene.dat=gene.dat,func=func,seed.table.in=seed.table.in)
  #Fit log-log models
  fit.loglog <- function(vals,N){
    ll.df <- data.frame("N"=N,"vals"=vals)
    ll.df <- ll.df[which(ll.df$N>=1000),]
    fit <- lm(log10(vals) ~ log10(N),data=ll.df)
    newdat <- data.frame("N"=sort(c(10^seq(0,7,0.025),max.N)))
    newdat$vals <- 10^predict.lm(fit,newdata=newdat,type="response")
    return(newdat)
  }
  fits <- lapply(plot.dat,function(df){
    mean.fit <- fit.loglog(vals=df$mean,N=df$N)
    lower.fit <- fit.loglog(vals=df$lower,N=df$N)
    upper.fit <- fit.loglog(vals=df$upper,N=df$N)
    data.frame("N"=mean.fit$N,
               "mean"=mean.fit$vals,
               "lower"=lower.fit$vals,
               "upper"=upper.fit$vals)
  })
  #Prepare plot area
  par(mar=mar,bty="n")
  plot(x=c(0,max.N),y=c(0,max(unlist(lapply(fits,function(df){df$upper[which(df$N<=max.N)]})))),
       type="n",xlab="",xaxt="n",ylab="",yaxt="n")
  axis(1,at=c(axTicks(1),10*max(axTicks(1))),labels=NA,tck=-0.03)
  sapply(axTicks(1),function(x){
    axis(1,at=x,tick=F,labels=paste(round(x/1000,0),"k",sep=""),line=-0.85,cex.axis=0.85)
  })
  mtext(1,line=1.1,text=xlabel)
  axis(2,at=c(axTicks(2),10*max(axTicks(2))),labels=NA,tck=-0.03)
  axis(2,at=axTicks(2),tick=F,labels=round(axTicks(2),2),line=-0.7,cex.axis=0.85,las=2)
  mtext(2,line=1.4,text=ylabel)
  #Plot curves & CIs
  points(fits$all$N,fits$all$mean,type="l",lwd=2,col=color)
  polygon(x=c(fits$all$N,rev(fits$all$N)),
          y=c(fits$all$lower,rev(fits$all$upper)),
          border=NA,bty="n",col=adjustcolor(color,alpha=0.15))
  points(fits$constrained$N,fits$constrained$mean,type="l",lwd=2,col=color,lty=5)
  polygon(x=c(fits$constrained$N,rev(fits$constrained$N)),
          y=c(fits$constrained$lower,rev(fits$constrained$upper)),
          border=NA,bty="n",col=adjustcolor(color,alpha=0.15))
  points(fits$unconstrained$N,fits$unconstrained$mean,type="l",lwd=2,col=color,lty=3)
  polygon(x=c(fits$unconstrained$N,rev(fits$unconstrained$N)),
          y=c(fits$unconstrained$lower,rev(fits$unconstrained$upper)),
          border=NA,bty="n",col=adjustcolor(color,alpha=0.15))
  #Legend
  if(legend==T){
    legend("topleft",lwd=2,lty=c(1,5,3),col=color,
           legend=c("All Genes","Constrained (pLI>0.9)","Unconstrained (pLI<0.1)"),
           cex=0.7,border=NA,bty="n")
  }
}
#Plot curves for fraction of genes with at least one functional vs sample size
plot.fracNonZeroGenes <- function(gene.dat,func,seed.table.in,color,max.N,
                           xlabel=NULL,ylabel=NULL,mar=c(2,4,0.3,0.3),legend=F){
  #Gather plotting data
  plot.dat <- get.pct.genes.bySampSize(gene.dat=gene.dat,func=func,seed.table.in=seed.table.in)
  #Fit log-log models
  fit.loglog <- function(vals,N){
    ll.df <- data.frame("N"=N,"vals"=vals)
    ll.df <- ll.df[which(ll.df$N>=1000),]
    fit <- lm(log10(vals) ~ log10(N),data=ll.df)
    newdat <- data.frame("N"=sort(c(10^seq(0,7,0.025),max.N)))
    newdat$vals <- 10^predict.lm(fit,newdata=newdat,type="response")
    return(newdat)
  }
  fits <- lapply(plot.dat,function(df){
    mean.fit <- fit.loglog(vals=df$mean,N=df$N)
    lower.fit <- fit.loglog(vals=df$lower,N=df$N)
    upper.fit <- fit.loglog(vals=df$upper,N=df$N)
    data.frame("N"=mean.fit$N,
               "mean"=mean.fit$vals,
               "lower"=lower.fit$vals,
               "upper"=upper.fit$vals)
  })
  #Prepare plot area
  par(mar=mar,bty="n")
  plot(x=c(0,max.N),y=c(0,max(unlist(lapply(fits,function(df){df$upper[which(df$N<=max.N)]})))),
       type="n",xlab="",xaxt="n",ylab="",yaxt="n")
  axis(1,at=c(axTicks(1),10*max(axTicks(1))),labels=NA,tck=-0.03)
  sapply(axTicks(1)[seq(1,length(axTicks(1)),2)],function(x){
    axis(1,at=x,tick=F,labels=paste(round(x/1000,0),"k",sep=""),line=-0.85,cex.axis=0.75)
  })
  mtext(1,line=1,text=xlabel)
  axis(2,at=c(axTicks(2),10*max(axTicks(2))),labels=NA,tck=-0.03)
  axis(2,at=axTicks(2),tick=F,labels=paste(round(100*axTicks(2),1),"%",sep=""),
       line=-0.7,cex.axis=0.75,las=2)
  mtext(2,line=1.8,text=ylabel)
  #Plot curves & CIs
  points(fits$all$N,fits$all$mean,type="l",lwd=2,col=color)
  polygon(x=c(fits$all$N,rev(fits$all$N)),
          y=c(fits$all$lower,rev(fits$all$upper)),
          border=NA,bty="n",col=adjustcolor(color,alpha=0.15))
  points(fits$constrained$N,fits$constrained$mean,type="l",lwd=2,col=color,lty=5)
  polygon(x=c(fits$constrained$N,rev(fits$constrained$N)),
          y=c(fits$constrained$lower,rev(fits$constrained$upper)),
          border=NA,bty="n",col=adjustcolor(color,alpha=0.15))
  points(fits$unconstrained$N,fits$unconstrained$mean,type="l",lwd=2,col=color,lty=3)
  polygon(x=c(fits$unconstrained$N,rev(fits$unconstrained$N)),
          y=c(fits$unconstrained$lower,rev(fits$unconstrained$upper)),
          border=NA,bty="n",col=adjustcolor(color,alpha=0.15))
  #Legend
  if(legend==T){
    legend("topleft",lwd=2,lty=c(1,5,3),col=color,
           legend=c("All Genes","Constrained (pLI>0.9)","Unconstrained (pLI<0.1)"),
           cex=0.7,border=NA,bty="n")
  }
}


###############################
###SV constraint power analyses
###############################
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
sina.plot <- function(vals,y.at,color,width=0.1,horiz=T,cex=0.25,median=T){
  j <- (width*sina.jitter(vals))+y.at
  if(horiz==T){
    points(x=vals,y=j,pch=21,cex=cex,col=color,bg=adjustcolor(color,alpha=0.3))
    if(median==T){
      segments(x0=median(vals),x1=median(vals),
               y0=y.at-width,y1=y.at+width,lwd=2)
    }else{
      segments(x0=mean(vals),x1=mean(vals),
               y0=y.at-width,y1=y.at+width,lwd=2)
    }
  }else{
    points(x=j,y=vals,pch=21,cex=cex,col=color,bg=adjustcolor(color,alpha=0.3))
    if(median==T){
      segments(x0=y.at-width,x1=y.at+width,y0=median(vals),y1=median(vals),lwd=1.5)
    }else{
      segments(x0=y.at-width,x1=y.at+width,y0=mean(vals),y1=mean(vals),lwd=1.5)
    }
  }
}

#Sina of SNV pLoF ObsExp for constrained vs not constrained genes
snvlof.obsexp.sina <- function(gene.dat,colors=c(svtypes$color[which(svtypes$svtype=="DEL")],
                                                 "grey70")){
  #Get plot data
  cons.oe <- gene.dat$ptv_oe[which(gene.dat$constrained==1)]
  uncons.oe <- gene.dat$ptv_oe[which(gene.dat$constrained==0)]
  uncons.oe <- uncons.oe[which(!is.na(uncons.oe))]
  ylims <- c(0,max(c(cons.oe,uncons.oe)))
  #Prep plot area
  par(mar=c(3.1,2.75,0.5,2.5),bty="n")
  plot(x=c(0,2),y=ylims,type="n",
       xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i")
  abline(h=1,lty=2)
  #Add sinas
  sina.plot(cons.oe,y.at=0.5,color=colors[1],width=0.35,horiz=F,median=F,cex=0.04)
  sina.plot(uncons.oe,y.at=1.5,color=colors[2],width=0.35,horiz=F,median=F,cex=0.04)
  #Add group labels
  axis(1,at=0.5,line=-1,tick=F,labels=bquote("pLI" > 0.9),font=2)
  axis(1,at=1.5,line=-1,tick=F,labels=bquote("pLI" <= 0.9),font=2)
  axis(1,at=0.5,line=-0.1,tick=F,cex.axis=0.9,
       labels=bquote(mu == .(format(round(mean(cons.oe),2),nsmall=2))))
  axis(1,at=0.5,line=0.6,tick=F,cex.axis=0.9,
       labels=bquote(sigma == .(format(round(sd(cons.oe),2),nsmall=2))))
  axis(1,at=1.5,line=-0.1,tick=F,cex.axis=0.9,
       labels=bquote(mu == .(format(round(mean(uncons.oe),2),nsmall=2))))
  axis(1,at=1.5,line=0.6,tick=F,cex.axis=0.9,
       labels=bquote(sigma == .(format(round(sd(uncons.oe),2),nsmall=2))))
  #Add axes
  axis(2,at=axTicks(2),tck=-0.02,labels=NA)
  axis(2,at=axTicks(2),tick=F,line=-0.6,labels=axTicks(2),las=2,cex.axis=0.8)
  mtext(2,line=1.5,text="SNV pLoF (Obs:Exp)")
  #Annotate difference between groups
  axis(4,at=c(mean(cons.oe),mean(uncons.oe)),lwd=2,tck=0.03,labels=NA)
  axis(4,at=mean(c(mean(cons.oe),mean(uncons.oe))),tick=F,las=2,line=-0.8,
       labels=paste(100*(1-round(mean(cons.oe)/mean(uncons.oe),3)),"%",sep=""))
}

#Predict average number of rare pLoF or CG per gene 
pergene.extrapolations <- function(gene.dat){
  #Helper to fit log-log models
  fit.loglog <- function(vals,N){
    ll.df <- data.frame("N"=N,"vals"=vals)
    ll.df <- ll.df[which(ll.df$N>=1000),]
    fit <- lm(log10(vals) ~ log10(N),data=ll.df)
    newdat <- data.frame("N"=sort(c(10^seq(0,8,0.025),max.N)))
    newdat$vals <- 10^predict.lm(fit,newdata=newdat,type="response")
    return(newdat)
  }
  #Fit models & project data
  projected.data <- lapply(c("lof.any","cg"),function(func){
    #Gather downsampling data
    down.dat <- get.gene.trend(gene.dat=gene.dat,func=func,seed.table.in=seed.table.in,invert.uncons=T)
    #Fit model
    fits <- lapply(down.dat,function(df){
      mean.fit <- fit.loglog(vals=df$mean,N=df$N)
      lower.fit <- fit.loglog(vals=df$lower,N=df$N)
      upper.fit <- fit.loglog(vals=df$upper,N=df$N)
      data.frame("N"=mean.fit$N,
                 "mean"=mean.fit$vals,
                 "lower"=lower.fit$vals,
                 "upper"=upper.fit$vals)
    })
    #Return just data frame of N & mean
    data.frame("N"=fits$unconstrained$N,"mean"=fits$unconstrained$mean)
  })
  projected.data <- merge(projected.data[[1]],projected.data[[2]],
                          by="N",suffixes=c(".lof",".cg"))
}

#Plot extrapolated number of SV per gene
plot.extrapolations <- function(projected.data,xmax=50000000){
  #Prep plot area
  if(is.null(xmax)){
    xmax <- max(projected.data$N)
  }
  ymax <- max(projected.data[which(projected.data$N<=xmax),-1])
  par(mar=c(2.5,3,0.3,0.3),bty="n")
  plot(x=c(0,xmax),y=c(0,ymax),
       type="n",xlab="",xaxt="n",ylab="",yaxt="n")
  #Plot lines
  points(x=projected.data$N,y=projected.data$mean.lof,
       type="l",lwd=2,col=svtypes$color[which(svtypes$svtype=="DEL")])
  points(x=projected.data$N,y=projected.data$mean.cg,
         type="l",lwd=2,col=svtypes$color[which(svtypes$svtype=="DUP")])
  #Add axes
  axis(1,at=axTicks(1),labels=NA,tck=-0.03)
  axis(1,at=axTicks(1),labels=axTicks(1)/1000000,tick=F,line=-0.75,cex.axis=0.75)
  mtext(1,text="Sample Size (Millions)",line=1)
  axis(2,at=axTicks(2),labels=NA,tck=-0.03)
  axis(2,at=axTicks(2),labels=axTicks(2),tick=F,line=-0.6,cex.axis=0.75,las=2)
  mtext(2,text="Mean SV per Gene",line=1.25)
}

#Simulate observed SV at a given sample size
simulate.sv <- function(gene.dat,mu,cv=0.36,theta=1.47,
                        constrained.obs.frac=0.218,
                        seed=123456789){
  #Simulate "true" count of expected SV per gene
  set.seed(seed)
  # #Normal
  # exp.true <- rnorm(n=nrow(gene.dat),
  #                     mean=mu,
  #                     sd=cv*mu)
  # exp.true[which(exp.true<0)] <- 0
  # #Poisson
  # exp.true <- rpois(n=nrow(gene.dat),
  #                   lambda=mu)
  # #Geometric
  # exp.true <- rgeom(n=nrow(gene.dat),prob=1/mu)
  #Negative binomial
  exp.true <- rnegbin(n=nrow(gene.dat),mu=mu,theta=theta)
  #Simulate labels
  sim.labels <- c(rep(1,length(which(gene.dat$constrained==1))),
                  rep(0,length(which(gene.dat$constrained!=1))))
  #Adjust observed values for constrained genes according to constrained.obs.frac
  obs <- exp.true
  obs[which(sim.labels==1)] <- constrained.obs.frac*obs[which(sim.labels==1)]
  obs <- round(obs)
  #Prepare simulated output
  sim.dat <- data.frame("obs"=obs,"exp.true"=exp.true,"constrained"=sim.labels)
  return(sim.dat)
}

#Inject noise into expected SV per gene
inject.noise <- function(exp.true,r2=0.5){
  #Get values
  n <- length(exp.true)
  rho <- sqrt(r2)
  theta <- acos(rho)
  x1 <- exp.true
  x2 <- rnorm(n,2,0.5)
  X <- cbind(x1,x2)
  Xctr <- scale(X,center=T,scale=F)
  Id <- diag(n)
  Q <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))
  P <- tcrossprod(Q)
  x2o <- (Id-P) %*% Xctr[,2]
  Xc2 <- cbind(Xctr[,1],x2o)
  Y <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))
  x <- Y[ , 2] + (1 / tan(theta)) * Y[,1]
  fit <- lm(x1 ~ x)
  newvals <- predict.lm(fit,newdata=data.frame("x1"=x1))
  newvals[which(newvals<0)] <- 0
  return(newvals)
}
#Alternative implementation (much faster)
inject.noise.2 <- function(exp.true,r2=0.5,theta=1.47,seed=123456789){
  y <- exp.true
  set.seed(seed)
  x <- rnbinom(length(y),size=theta,mu=mean(exp.true))
  rho <- sqrt(r2)
  y.perp <- residuals(glm.nb(x ~ y))
  y.prescale <- rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
  fit <- lm(exp.true ~ y.prescale)
  newvals <- predict.lm(fit,newdata=data.frame("y.prescale"=y.prescale))
  newvals[which(newvals<0)] <- 0
  return(newvals)
}

#Estimate fraction of constrained genes captured under model spec at a single sample size
simulate.power <- function(gene.dat,mu,model.r2,
                           cv=0.36,theta=1.47,constrained.obs.frac=0.218,
                           constrained.min.p=0.001,
                           seed=123456789){
  #Old simulation-based implementation:
  # #Simulate SV data per gene
  # sim.dat <- simulate.sv(gene.dat=gene.dat,
  #                        mu=mu,cv=cv,theta=theta,
  #                        constrained.obs.frac=constrained.obs.frac,
  #                        seed=seed)
  # sim.dat$exp.noised <- inject.noise.2(sim.dat$exp.true,r2=model.r2)
  # #Calculate p-value from chisquared test given expected number of variants
  # testable <- sim.dat[which(sim.dat$exp.noised>=10),]
  # p <- pchisq(q=testable$obs,df=testable$exp.noised)
  # #Return fraction of constrained genes with p > constrained.min.p
  # length(which(p<constrained.min.p & testable$constrained==1))/length(which(sim.dat$constrained==1))
  #New simpler implementation:
  #one-sample t-test
  # power <- as.numeric(pwr.t.test(n=mu,d=(1-constrained.obs.frac)*model.r2,
  #            sig.level=0.1,type="one.sample",alternative="greater")$power)
  #one-sample proportion test
  if(mu>1){
    power <- as.numeric(pwr.p.test(h=(1-constrained.obs.frac)*model.r2,
                                   n=mu,
                                   sig.level=0.1,
                                   alternative="greater")$power)
  }else{
    power <- NA
  }
  if(is.nan(power)){power <- NA}
  return(power)
}

#Calculate power for all models and all sample sizes
simulate.power.all <- function(gene.dat,
                               projected.data,
                               model.r2s=seq(0.1,0.9,0.1)){
  #pLoF
  lof.power <- as.data.frame(t(sapply(projected.data[,2],
                function(mu){
    sapply(model.r2s,function(model.r2){
      simulate.power(gene.dat=gene.dat,
                     mu=mu,
                     model.r2=model.r2)
    })
  })))
  lof.power <- cbind(projected.data$N,lof.power)
  colnames(lof.power) <- c("N",paste("lof",model.r2s,sep="."))
  #CG
  cg.power <- as.data.frame(t(sapply(projected.data[,3],
                                      function(mu){
                                        sapply(model.r2s,function(model.r2){
                                          simulate.power(gene.dat=gene.dat,
                                                         mu=mu,
                                                         model.r2=model.r2)
                                        })
                                      })))
  cg.power <- cbind(projected.data$N,cg.power)
  colnames(cg.power) <- c("N",paste("cg",model.r2s,sep="."))
  return(list("lof"=lof.power,
              "cg"=cg.power))
}

#Plot power curves
plot.power <- function(power.dat,xlims=NULL){
  #Get plot values
  if(is.null(xlims)){
    xmin <- min(power.dat$N[which(!is.na(power.dat[,10]))])
    xmax <- max(power.dat$N)
    xlims <- c(xmin,xmax)
  }
  col.pal <- colorRampPalette(c("#440154","#365C8C","#25A584","#FDE725"))(ncol(power.dat)-1)
  #Prep plot area
  par(mar=c(2.5,3,0.25,1),bty="n")
  plot(x=log10(xlims),y=c(0,1),type="n",
       xaxt="n",yaxt="n",xlab="",ylab="")
  abline(h=0.5,lty=2)
  #Plot power curves
  sapply(2:ncol(power.dat),function(i){
    points(x=log10(power.dat$N),y=power.dat[,i],
           type="l",lwd=2,col=col.pal[i-1])
  })
  #Add axes
  axis(1,at=log10(as.numeric(unlist(sapply(3:9,function(i){(1:9)*(10^i)})))),
       labels=NA,tck=-0.015)
  axis(1,at=3:8,labels=NA,tck=-0.03)
  axis(1,at=3:8,line=-0.85,cex.axis=0.75,tick=F,
       labels=c("1k","10k","100k","1M","10M","100M"))
  mtext(1,line=1,text="Sample Size")
  axis(2,at=seq(0,1,0.2),labels=NA,tck=-0.03)
  axis(2,at=seq(0,1,0.2),tick=F,line=-0.6,cex.axis=0.75,las=2,
       labels=paste(seq(0,100,20),"%",sep=""))
  mtext(2,line=1.75,text=bquote("Power (1-"*beta*")"))
}


################
###RSCRIPT BLOCK
################
require(optparse,quietly=T)
require(boot,quietly=T)
require(pwr,quietly=T)
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
args <- parse_args(OptionParser(usage="%prog singletons_per_sample total_sv_sites sv_per_gene OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 8){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
vcf2bed.in <- args$args[1]
singletons.in <- args$args[2]
totalsites.in <- args$args[3]
svpergene.in <- args$args[4]
SNVdata.in <- args$args[5]
seed.table.in <- args$args[6]
medians.perSample.in <- args$args[7]
OUTDIR <- args$args[8]
prefix <- opts$prefix
svtypes.file <- opts$svtypes
pops.file <- opts$populations

# #Dev parameters (local)
# vcf2bed.in <- "~/scratch/gnomAD_v2_SV_MASTER.vcf2bed.bed.gz"
# singletons.in <- "~/scratch/gnomAD_v2_SV_MASTER.downsampling.singletons_per_sample.txt.gz"
# totalsites.in <- "~/scratch/gnomAD_v2_SV_MASTER.downsampling.total_sv_sites.txt.gz"
# svpergene.in <- "~/scratch/gnomAD_v2_SV_MASTER.downsampling.sv_per_gene.txt.gz"
# SNVdata.in <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/gnomAD_v2.1_canonical_constraint.condensed.txt.gz"
# seed.table.in <- "~/scratch/gnomAD_v2_SV_MASTER.downsampling_seeds.txt"
# medians.perSample.in <- "~/scratch/gnomAD_v2_SV_MASTER.median_sites_per_sample_by_population.txt"
# OUTDIR <- "~/scratch/gnomAD_v2_test_plots/"
# svtypes.file <- "~/Desktop/Collins/Talkowski/code/sv-pipeline/ref/vcf_qc_refs/SV_colors.txt"
# prefix <- "gnomAD_v2_SV_MASTER"
# pops.file <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/gnomAD_population_colors.txt"

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

###Sets populations & colors
if(!is.null(pops.file)){
  pops <- read.table(pops.file,sep="\t",header=T,comment.char="",check.names=F)
}

###Process input data
cat("NOW LOADING AND CLEANING DATA\n")
sv.sites <- import.sites(totalsites.in)
singletons <- import.singletons(singletons.in,seed.table.in,pops)
vcf2bed <- read.table(vcf2bed.in,header=T,comment.char="")
vcf2bed <- vcf2bed[which(vcf2bed$FILTER %in% c("PASS","MULTIALLELIC")),]
empirical.counts <- sapply(svtypes$svtype,function(svtype){length(which(vcf2bed$SVTYPE==svtype))})
empirical.counts <- c(sum(empirical.counts),empirical.counts)
names(empirical.counts) <- c("ALL",svtypes$svtype)
max.N <- max(vcf2bed$N_BI_GENOS,na.rm=T)
rm(vcf2bed)
med.sitesPerSample <- read.table(medians.perSample.in,header=T,sep="\t")
gene.dat <- import.genes(svpergene.in,SNVdata.in)

###Plot total sites discovered as a function of sample size
pdf(paste(OUTDIR,"/",prefix,".total_sites_discovered_bySampleSize.pdf",sep=""),
    height=2,width=1.8)
plot.totalSites(sv.sites,svtypes,empirical.counts,max.N)
dev.off()

###Plot singletons per sample as a function of sample size
pdf(paste(OUTDIR,"/",prefix,".singletons_per_sample_bySampleSize.pdf",sep=""),
    height=2.5,width=2.5)
plot.singletons(singletons,pops[which(!(pops$pop %in% c("OTH","SAS"))),],
                med.sitesPerSample=med.sitesPerSample,
                max.N=max.N,ymax=10000)
dev.off()

###Plot panels of functional SV per gene
pdf(paste(OUTDIR,"/",prefix,".lof_cg_per_gene_bySampleSize.pdf",sep=""),
    height=3,width=2.25)
par(mfrow=c(2,1))
plot.svpergene(gene.dat=gene.dat,func="lof.any",
               seed.table.in=seed.table.in,
               color=svtypes$color[which(svtypes$svtype=="DEL")],
               max.N=max.N,xlabel="",ylabel="pLoF SV / Gene",
               mar=c(1.5,3,1.8,0.3))
plot.svpergene(gene.dat=gene.dat,func="cg",
               seed.table.in=seed.table.in,
               color=svtypes$color[which(svtypes$svtype=="DUP")],
               max.N=max.N,xlabel="gnomAD Sample Size",ylabel="CG SV / Gene")
dev.off()

###Plot panels of non-zero genes
pdf(paste(OUTDIR,"/",prefix,".lof_cg_nonZero_genes_bySampleSize.pdf",sep=""),
    height=2.5,width=2.25)
par(mfrow=c(2,1))
plot.fracNonZeroGenes(gene.dat=gene.dat,func="lof.any",
               seed.table.in=seed.table.in,
               color=svtypes$color[which(svtypes$svtype=="DEL")],
               max.N=max.N,xlabel="",ylabel="Genes with\nat least 1 SV",
               mar=c(1.5,4,0.8,0.3))
text(x=par("usr")[1],y=par("usr")[4]-(0.15*(par("usr")[4]-par("usr")[3])),
     pos=4,labels="pLoF",cex=0.9)
plot.fracNonZeroGenes(gene.dat=gene.dat,func="cg",
               seed.table.in=seed.table.in,
               color=svtypes$color[which(svtypes$svtype=="DUP")],
               max.N=max.N,xlabel="Sample Size",ylabel="Genes with\nat least 1 SV")
text(x=par("usr")[1],y=par("usr")[4]-(0.15*(par("usr")[4]-par("usr")[3])),
     pos=4,labels="CG",cex=0.9)
dev.off()


###SV constraint analyses
#Plot of empirical SNV constraint data
png(paste(OUTDIR,"/",prefix,".snv_lof_constraint_obsexp.png",sep=""),
    height=1100,width=1000,res=400)
snvlof.obsexp.sina(gene.dat)
dev.off()
#Project mean number of SV per gene (pLI≤0.9)
projected.data <- pergene.extrapolations(gene.dat)
#Plot projections
pdf(paste(OUTDIR,"/",prefix,".projected_lof_cg_per_gene_vs_sample_size.pdf",sep=""),
    height=2.25,width=2.75)
plot.extrapolations(projected.data=projected.data,xmax=50000000)
dev.off()
#Estimate power vs sample size
power.data <- simulate.power.all(gene.dat=gene.dat,
                                projected.data=projected.data)
#Plot power curves
pdf(paste(OUTDIR,"/",prefix,".predicted_sv_constraint_power.pdf",sep=""),
    height=2.25,width=5.25)
par(mfrow=c(1,2))
plot.power(power.dat=power.data$lof,
           xlims=c(100000,50000000))
plot.power(power.dat=power.data$cg,
           xlims=c(100000,50000000))
dev.off()

