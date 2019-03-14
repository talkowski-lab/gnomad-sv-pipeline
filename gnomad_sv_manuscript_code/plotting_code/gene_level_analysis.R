#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Perform gene-level analyses for formal gnomAD analysis


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
#Process SNV gene data
read.snvdata <- function(SNVdata.in,gene.metadata.in){
  #Read & clean
  snv.data <- read.table(SNVdata.in,header=T,comment.char="")
  metadata <- read.table(gene.metadata.in,header=T)
  merged <- merge(x=snv.data,y=metadata,by="gene",sort=F)
  merged <- merged[which(!(merged$chrom %in% c("chrX","chrY"))),]
  #Assign oe deciles
  merged$mis_oe_dec <- ceiling(10*rank(merged$mis_oe)/(nrow(merged)+1))
  merged$ptv_oe_dec <- ceiling(10*rank(merged$ptv_oe)/(nrow(merged)+1))
  #Assign oe to 40 bins
  merged$mis_oe_binrank <- ceiling(40*rank(merged$mis_oe)/(nrow(merged)+1))
  merged$ptv_oe_binrank <- ceiling(40*rank(merged$ptv_oe)/(nrow(merged)+1))
  #Assign oe percentiles
  merged$mis_oe_cent <- ceiling(100*rank(merged$mis_oe)/(nrow(merged)+1))
  merged$ptv_oe_cent <- ceiling(100*rank(merged$ptv_oe)/(nrow(merged)+1))
  #Return formatted data
  return(merged)
}
#Read external gene list and restrict to autosomal genes considered in gnomAD
importGenelist <- function(genelist.dir,filename,gene.data){
  genes <- read.table(paste(genelist.dir,"/",filename,sep=""),header=F)
  genes <- as.character(genes[,1])
  genes <- genes[which(genes %in% gene.data$gene)]
  return(genes)
}

#Count fraction of SV with functional effects by SVTYPE
count.fracLOFSV <- function(dat){
  subdat <- dat[-unique(c(grep("MULTIALLELIC",dat$FILTER,fixed=T),
                          grep("UNRESOLVED",dat$FILTER,fixed=T))),]
  subdat <- dat[which(!(dat$chrom %in% c("X","Y"))),]
  #All SV - pLoF
  all.lof <- length(which(!is.na(subdat$PROTEIN_CODING__LOF)))/nrow(subdat)
  #Deletions - pLoF
  del.lof <- length(which(subdat$SVTYPE=="DEL" & !is.na(subdat$PROTEIN_CODING__LOF)))/length(which(subdat$SVTYPE=="DEL"))
  #Insertions - pLoF
  ins.lof <- length(which(subdat$SVTYPE=="INS" & !is.na(subdat$PROTEIN_CODING__LOF)))/length(which(subdat$SVTYPE=="INS"))
  #Inversions - pLoF
  inv.lof <- length(which(subdat$SVTYPE=="INV" & !is.na(subdat$PROTEIN_CODING__LOF)))/length(which(subdat$SVTYPE=="INV"))
  #Complex - pLoF
  cpx.lof <- length(which(subdat$SVTYPE=="CPX" & !is.na(subdat$PROTEIN_CODING__LOF)))/length(which(subdat$SVTYPE=="CPX"))
  #DUP/CPX - CG
  dup.cg <- length(which(subdat$SVTYPE %in% c("DUP","CPX") & !is.na(subdat$PROTEIN_CODING__COPY_GAIN)))/length(which(subdat$SVTYPE %in% c("DUP","CPX")))
  #Duplications - pLoF
  dup.lof <- length(which(subdat$SVTYPE=="DUP" & !is.na(subdat$PROTEIN_CODING__DUP_LOF)))/length(which(subdat$SVTYPE=="DUP"))
  # #Duplications - partial gene
  # dup.partial <- length(which(subdat$SVTYPE=="DUP" & !is.na(subdat$PROTEIN_CODING__DUP_PARTIAL)))/length(which(subdat$SVTYPE=="DUP"))
  #Inversions + complex - span
  invcpx.span <- length(which(subdat$SVTYPE %in% c("INV","CPX") & !is.na(subdat$PROTEIN_CODING__INV_SPAN)))/length(which(subdat$SVTYPE %in% c("INV","CPX")))
  #Intronic, no coding effect
  all.intronic <- length(which(is.na(subdat$PROTEIN_CODING__LOF) 
                               & is.na(subdat$PROTEIN_CODING__DUP_LOF)
                               & is.na(subdat$PROTEIN_CODING__COPY_GAIN)
                               & is.na(subdat$PROTEIN_CODING__DUP_PARTIAL)
                               & is.na(subdat$PROTEIN_CODING__MSV_EXON_OVR)
                               & !is.na(subdat$PROTEIN_CODING__INTRONIC)))/nrow(subdat)
  #Promoter, no coding effect
  all.promoter <- length(which(is.na(subdat$PROTEIN_CODING__LOF) 
                               & is.na(subdat$PROTEIN_CODING__DUP_LOF)
                               & is.na(subdat$PROTEIN_CODING__COPY_GAIN)
                               & is.na(subdat$PROTEIN_CODING__DUP_PARTIAL)
                               & is.na(subdat$PROTEIN_CODING__MSV_EXON_OVR)
                               & !is.na(subdat$PROTEIN_CODING__PROMOTER)))/nrow(subdat)
  #Format output
  v.out <- c("all.lof"=all.lof,"del.lof"=del.lof,"ins.lof"=ins.lof,"inv.lof"=inv.lof,"cpx.lof"=cpx.lof,
             "dup.cg"=dup.cg,"dup.lof"=dup.lof,"invcpx.span"=invcpx.span,"all.intronic"=all.intronic,"all.promoter"=all.promoter)
  return(v.out)
}
#Barplots of SV fractions by functional category
barplot.fracByEffect <- function(frac.mat,base.colors,category.labels){
  #Prep plot area
  ymax <- 1.05*max(frac.mat,na.rm=T)
  par(mar=c(2.5,3,1,0.5),bty="n")
  plot(x=c(0,nrow(frac.mat)),y=c(0,ymax),type="n",
       xlab="",xaxt="n",ylab="",yaxt="n",yaxs="i")
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=par("usr")[3],ytop=par("usr")[4],
       col="white",border=NA,bty="n")
  if(length(axTicks(2))>6){
    y.at <- axTicks(2)[seq(1,length(axTicks(2)),2)]
    y.at <- c(y.at,y.at[length(y.at)]+(y.at[2]-y.at[1]))
  }else{
    y.at <- axTicks(2)
  }
  axis(2,at=y.at,labels=NA,tck=-0.06)
  axis(2,at=y.at,tick=F,cex.axis=0.8,las=2,line=-0.6,
       labels=paste(round(100*y.at,0),"%",sep=""))
  mtext(2,line=1.9,text="Pct. of Sites")
  
  #Add bars
  sapply(1:nrow(frac.mat),function(i){
    rect(xleft=i-0.85,xright=i-0.75,ybottom=0,ytop=frac.mat[i,1],col="black")
    rect(xleft=i-c(0.75,0.55,0.35),xright=i-c(0.55,0.35,0.15),
         ybottom=0,ytop=frac.mat[i,2:4],
         col=c(adjustcolor(base.colors[i],alpha=0.15),
               adjustcolor(base.colors[i],alpha=0.5),
               base.colors[i]))
    axis(1,at=i-0.5,tick=F,labels=category.labels[i],cex.axis=0.8,padj=1,line=-1.5)
    text(x=i-0.8,y=frac.mat[i,1]-(0.05*(par("usr")[4]-par("usr")[3])),pos=3,xpd=T,
         labels=paste(round(100*frac.mat[i,1],1),"%",sep=""),cex=0.85)
  })
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
#Count fraction of singletons per functional class
count.fracSingletons <- function(dat){
  subdat <- dat[-unique(c(grep("MULTIALLELIC",dat$FILTER,fixed=T),
                          grep("UNRESOLVED",dat$FILTER,fixed=T))),]
  subdat <- dat[which(!(dat$chrom %in% c("X","Y"))),]
  lof.all.idx <- which(!is.na(subdat$PROTEIN_CODING__LOF))
  lof.del.idx <- which(!is.na(subdat$PROTEIN_CODING__LOF) & subdat$SVTYPE=="DEL")
  lof.invcpx.idx <- which(!is.na(subdat$PROTEIN_CODING__LOF) & subdat$SVTYPE %in% c("INV","CPX"))
  lof.ins.idx <- which(!is.na(subdat$PROTEIN_CODING__LOF) & subdat$SVTYPE=="INS")
  # lof.dup.idx <- which(!is.na(subdat$PROTEIN_CODING__LOF) & subdat$SVTYPE=="DUP")
  cg.idx <- which(!is.na(subdat$PROTEIN_CODING__COPY_GAIN))
  dup.lof.idx <- which(!is.na(subdat$PROTEIN_CODING__DUP_LOF))
  inv.span.idx <- which(!is.na(subdat$PROTEIN_CODING__INV_SPAN))
  intronic.idx <- which(is.na(subdat$PROTEIN_CODING__LOF) 
                        & is.na(subdat$PROTEIN_CODING__DUP_LOF)
                        & is.na(subdat$PROTEIN_CODING__COPY_GAIN)
                        & is.na(subdat$PROTEIN_CODING__DUP_PARTIAL)
                        & is.na(subdat$PROTEIN_CODING__MSV_EXON_OVR)
                        & !is.na(subdat$PROTEIN_CODING__INTRONIC))
  promoter.idx <- which(is.na(subdat$PROTEIN_CODING__LOF) 
                        & is.na(subdat$PROTEIN_CODING__DUP_LOF)
                        & is.na(subdat$PROTEIN_CODING__COPY_GAIN)
                        & is.na(subdat$PROTEIN_CODING__DUP_PARTIAL)
                        & is.na(subdat$PROTEIN_CODING__MSV_EXON_OVR)
                        & is.na(subdat$PROTEIN_CODING__INTRONIC)
                        & !is.na(subdat$PROTEIN_CODING__PROMOTER))
  intergenic.idx <- which(is.na(subdat$PROTEIN_CODING__LOF) 
                        & is.na(subdat$PROTEIN_CODING__DUP_LOF)
                        & is.na(subdat$PROTEIN_CODING__COPY_GAIN)
                        & is.na(subdat$PROTEIN_CODING__DUP_PARTIAL)
                        & is.na(subdat$PROTEIN_CODING__MSV_EXON_OVR)
                        & is.na(subdat$PROTEIN_CODING__INTRONIC)
                        & is.na(subdat$PROTEIN_CODING__PROMOTER)
                        & subdat$PROTEIN_CODING__INTERGENIC=="True")
  all.sv.idx <- 1:nrow(subdat)
  fracSingles <- do.call("rbind", lapply(list(lof.all.idx,lof.del.idx,lof.ins.idx,lof.invcpx.idx,
                            cg.idx,dup.lof.idx,inv.span.idx,
                            intronic.idx,promoter.idx,intergenic.idx,
                            all.sv.idx),function(idx){
                              c(length(idx),
                                calc.fracSingletons.singleClass(dat=subdat[idx,]))
                            }))
  colnames(fracSingles) <- c("N_SV","fracSingletons","lower95CI","upper95CI")
  #Format output
  fracSingles <- cbind("Annotation"=c("LOF.ANY","LOF.DEL","LOF.INS","LOF.INV_CPX",
                                      "CG","DUP_LOF","INV_SPAN",
                                      "INTRONIC","PROMOTER","INTERGENIC","ALL.SV"),
                       fracSingles)
  return(fracSingles)
}
#Count fraction of singletons per SVTYPE, split by those with and without coding effects
count.fracSingletons.coding_vs_noncoding <- function(dat){
  subdat <- dat[-unique(c(grep("MULTIALLELIC",dat$FILTER,fixed=T),
                          grep("UNRESOLVED",dat$FILTER,fixed=T))),]
  subdat <- dat[which(!(dat$chrom %in% c("X","Y"))),]
  sv.with.coding.effects <- which(!is.na(subdat$PROTEIN_CODING__LOF)
                                  | !is.na(subdat$PROTEIN_CODING__DUP_LOF)
                                  | !is.na(subdat$PROTEIN_CODING__COPY_GAIN)
                                  | !is.na(subdat$PROTEIN_CODING__MSV_EXON_OVR))
  intergenic.sv <- setdiff(which(subdat$PROTEIN_CODING__INTERGENIC=="True"),
                           sv.with.coding.effects)
  all.sv.idx <- 1:nrow(subdat)
  del.lof <- intersect(which(subdat$SVTYPE=="DEL"),sv.with.coding.effects)
  del.noncoding <- intersect(which(subdat$SVTYPE=="DEL"),intergenic.sv)
  dup.ied <- intersect(which(subdat$SVTYPE=="DUP" & !is.na(subdat$PROTEIN_CODING__DUP_LOF)),
                       sv.with.coding.effects)
  dup.cg <- intersect(which(subdat$SVTYPE=="DUP" & !is.na(subdat$PROTEIN_CODING__COPY_GAIN)),
                       sv.with.coding.effects)
  dup.noncoding <- intersect(which(subdat$SVTYPE=="DUP"),intergenic.sv)
  ins.lof <- intersect(which(subdat$SVTYPE=="INS"),sv.with.coding.effects)
  ins.noncoding <- intersect(which(subdat$SVTYPE=="INS"),intergenic.sv)
  inv.lof <- intersect(which(subdat$SVTYPE=="INV"),sv.with.coding.effects)
  inv.noncoding <- intersect(which(subdat$SVTYPE=="INV"),intergenic.sv)
  cpx.lof <- intersect(which(subdat$SVTYPE=="CPX" & !is.na(subdat$PROTEIN_CODING__LOF)),
                       sv.with.coding.effects)
  cpx.cg <- intersect(which(subdat$SVTYPE=="CPX" & !is.na(subdat$PROTEIN_CODING__COPY_GAIN)),
                       sv.with.coding.effects)
  cpx.noncoding <- intersect(which(subdat$SVTYPE=="CPX"),intergenic.sv)
  fracSingles <- do.call("rbind", lapply(list(all.sv.idx,del.lof,del.noncoding,
                                              dup.ied,dup.cg,dup.noncoding,
                                              ins.lof,ins.noncoding,inv.lof,inv.noncoding,
                                              cpx.lof,cpx.cg,cpx.noncoding),function(idx){
                                                c(length(idx),
                                                  calc.fracSingletons.singleClass(dat=subdat[idx,]))
                                              }))
  colnames(fracSingles) <- c("N_SV","fracSingletons","lower95CI","upper95CI")
  #Format output
  fracSingles <- cbind("Annotation"=c("ALL.SV","DEL.LOF","DEL.NON",
                                      "DUP.IED","DUP.CG","DUP.NON",
                                      "INS.LOF","INS.NON","INV.LOF","INV.NON",
                                      "CPX.LOF","CPX.CG","CPX.NON"),
                       fracSingles)
  return(fracSingles)
}
#Dotplot of fraction of singletons per functional annotation by class
dotplot.fracSingletons <- function(fracSingletons,horiz=F,add.pt.labs=F,add.N=T){
  plot.dat <- as.data.frame(apply(fracSingletons[,-1],2,as.numeric))
  rownames(plot.dat) <- fracSingletons[,1]
  class.colors <- c(rep(svtypes$color[which(svtypes$svtype=="DEL")],4),
                    svtypes$color[which(svtypes$svtype=="DUP")],
                    svtypes$color[which(svtypes$svtype=="MCNV")],
                    svtypes$color[which(svtypes$svtype=="INV")],
                    rep("grey35",2),"grey70","black")
  class.labels <- c("pLoF (All)","pLoF (DEL)","pLoF (INS)","pLoF (INV & CPX)",
                    "CG","IED","Whole-Gene INV",
                    "Intronic","Promoter","Intergenic","All Autosomal SV")
  new.order <- rev(c(which(rownames(plot.dat)=="ALL.SV"),
                     which(rownames(plot.dat)=="LOF.ANY"),
                     intersect(order(-plot.dat$fracSingletons),
                               which(!(rownames(plot.dat) %in% c("ALL.SV","LOF.ANY"))))))
  if(horiz==T){
    new.order <- rev(new.order)
  }
  plot.dat <- plot.dat[new.order,]
  class.colors <- class.colors[new.order]
  class.labels <- class.labels[new.order]
  if(horiz==T){
    par(mar=c(3.5,2.5,0.25,0.25),bty="n")
    plot(y=range(plot.dat[,-1]),
         x=c(0.4,nrow(plot.dat)-0.1),
         type="n",xaxt="n",xlab="",yaxt="n",ylab="")
    abline(h=plot.dat$fracSingletons[which(rownames(plot.dat)=="ALL.SV")],lty=2)
    segments(y0=plot.dat[,3],y1=plot.dat[,4],
             x0=(1:nrow(plot.dat))-0.5,
             x1=(1:nrow(plot.dat))-0.5,
             lwd=2,lend="round",
             col=class.colors)
    points(y=plot.dat[,2],x=c(1:nrow(plot.dat))-0.5,pch=19,
           col=class.colors)
    if(add.pt.labs==T){
      text(x=c(2:nrow(plot.dat))-0.65,y=plot.dat[-1,2],
           col=class.colors[-1],cex=0.5,pos=4,xpd=T,
           labels=paste(format(round(100*plot.dat[-1,2],1),nsmall=1),"%",sep=""))
      text(x=1-0.65,y=plot.dat[1,2]+0.01,
           col=class.colors[1],cex=0.5,pos=4,xpd=T,
           labels=paste(format(round(100*plot.dat[1,2],1),nsmall=1),"%",sep=""))
    }
    par(xpd=T)
    sapply(1:nrow(plot.dat),function(i){
      if(add.N==T){
        text(x=i-0.15,y=par("usr")[3]-(0.03*(par("usr")[4]-par("usr")[3])),
             labels=class.labels[i],srt=40,pos=2,cex=0.7)
        text(x=i+0.175,y=par("usr")[3]-(0.05*(par("usr")[4]-par("usr")[3])),srt=40,col="gray40",cex=0.5,pos=2,
             labels=paste("N=",prettyNum(plot.dat[i,1],big.mark=",")," SVs",sep=""))
      }else{
        text(x=i-0.08,y=par("usr")[3]-(0.03*(par("usr")[4]-par("usr")[3])),
             labels=class.labels[i],srt=40,pos=2,cex=0.7)
      }
    })
    par(xpd=F)
    axis(2,at=seq(0.4,0.7,0.1),labels=NA,tck=-0.04)
    sapply(seq(0.4,0.7,0.1),function(x){
      axis(2,at=x,line=-0.7,cex.axis=0.65,tick=F,las=2,
           labels=paste(round(x*100,0),"%",sep=""))
    })
    mtext(2,text="Singleton Proportion",line=1.5,cex=0.8)
  }else{
    par(mar=c(0.5,6,2.25,0.5),bty="n")
    plot(x=range(plot.dat[,-1]),
         y=c(0.25,nrow(plot.dat)-0.25),
         type="n",xaxt="n",xlab="",yaxt="n",ylab="")
    abline(v=plot.dat$fracSingletons[which(rownames(plot.dat)=="ALL.SV")],lty=2)
    segments(x0=plot.dat[,3],x1=plot.dat[,4],
             y0=(1:nrow(plot.dat))-0.5,
             y1=(1:nrow(plot.dat))-0.5,
             lwd=2,lend="round",
             col=class.colors)
    points(x=plot.dat[,2],y=c(1:nrow(plot.dat))-0.5,pch=19,
           col=class.colors)
    axis(2,at=(1:nrow(plot.dat))-0.3,tick=F,line=-0.8,las=2,cex.axis=0.8,
         labels=class.labels)
    sapply(1:nrow(plot.dat),function(i){
      axis(2,at=i-0.7,tick=F,line=-0.8,las=2,cex.axis=0.65,
           labels=paste("N=",prettyNum(plot.dat[i,1],big.mark=",")," SVs",sep=""),
           col.axis=class.colors[i])
    })
    axis(3,at=seq(0.4,0.7,0.1),labels=NA,tck=-0.04)
    sapply(seq(0.4,0.7,0.1),function(x){
      axis(3,at=x,line=-0.7,cex.axis=0.65,tick=F,
           labels=paste(round(x*100,0),"%",sep=""))
    })
    mtext(3,text="Singleton Proportion",line=1.15,cex=0.9)
  }
}
#Dotplot of fraction of singletons per functional annotation by class
dotplot.fracSingletons.coding_vs_non <- function(fracSingletons,add.N=T){
  plot.dat <- as.data.frame(apply(fracSingletons.coding_vs_non[,-1],2,as.numeric))
  rownames(plot.dat) <- fracSingletons.coding_vs_non[,1]
  class.colors <- c("black",
                    svtypes$color[which(svtypes$svtype=="DEL")],
                    adjustcolor(svtypes$color[which(svtypes$svtype=="DEL")],alpha=0.3),
                    rep(svtypes$color[which(svtypes$svtype=="DUP")],2),
                    adjustcolor(svtypes$color[which(svtypes$svtype=="DUP")],alpha=0.3),
                    svtypes$color[which(svtypes$svtype=="INS")],
                    adjustcolor(svtypes$color[which(svtypes$svtype=="INS")],alpha=0.3),
                    svtypes$color[which(svtypes$svtype=="INV")],
                    adjustcolor(svtypes$color[which(svtypes$svtype=="INV")],alpha=0.3),
                    rep(svtypes$color[which(svtypes$svtype=="CPX")],2),
                    adjustcolor(svtypes$color[which(svtypes$svtype=="CPX")],alpha=0.3))
  label.colors <- c("black",
                    rep(svtypes$color[which(svtypes$svtype=="DEL")],2),
                    rep(svtypes$color[which(svtypes$svtype=="DUP")],3),
                    rep(svtypes$color[which(svtypes$svtype=="INS")],2),
                    rep(svtypes$color[which(svtypes$svtype=="INV")],2),
                    rep(svtypes$color[which(svtypes$svtype=="CPX")],3))
    class.labels <- c("All Autosomal SV","DEL (pLoF)","DEL (Intergenic)",
                    "DUP (IED)","DUP (CG)","DUP (Intergenic)",
                    "INS (pLoF)","INS (Intergenic)","INV (pLoF)","INV (Intergenic)",
                    "CPX (pLoF)","CPX (CG)","CPX (Intergenic)")
    par(mar=c(3.5,2.25,0.25,0.25),bty="n")
    plot(y=range(plot.dat[,-1]),
         x=c(0.25,nrow(plot.dat)-0.4),
         type="n",xaxt="n",xlab="",yaxt="n",ylab="")
    abline(h=plot.dat$fracSingletons[which(rownames(plot.dat)=="ALL.SV")],lty=2)
    segments(y0=plot.dat[,3],y1=plot.dat[,4],
             x0=(1:nrow(plot.dat))-0.5,
             x1=(1:nrow(plot.dat))-0.5,
             lwd=2,lend="round",
             col="white")
    segments(y0=plot.dat[,3],y1=plot.dat[,4],
             x0=(1:nrow(plot.dat))-0.5,
             x1=(1:nrow(plot.dat))-0.5,
             lwd=2,lend="round",
             col=class.colors)
    points(y=plot.dat[,2],x=c(1:nrow(plot.dat))-0.5,pch=19,
           col="white")
    points(y=plot.dat[,2],x=c(1:nrow(plot.dat))-0.5,pch=19,
           col=class.colors)
    # text(x=c(1:nrow(plot.dat))-0.625,y=plot.dat[,2],
    #      col=label.colors,cex=0.5,pos=4,xpd=T,
    #      labels=paste(format(round(100*plot.dat[,2],1),nsmall=1),"%",sep=""))
    par(xpd=T)
    sapply(1:nrow(plot.dat),function(i){
      if(add.N==T){
        text(x=i-0.15,y=par("usr")[3]-(0.03*(par("usr")[4]-par("usr")[3])),
             labels=class.labels[i],srt=40,pos=2,cex=0.7)
        text(x=i+0.175,y=par("usr")[3]-(0.05*(par("usr")[4]-par("usr")[3])),srt=40,col="gray40",cex=0.5,pos=2,
             labels=paste("N=",prettyNum(plot.dat[i,1],big.mark=",")," SV",sep=""))
      }else{
        text(x=i-0.08,y=par("usr")[3]-(0.03*(par("usr")[4]-par("usr")[3])),
             labels=class.labels[i],srt=40,pos=2,cex=0.7)
      }
    })
    par(xpd=F)
    axis(2,at=seq(0.4,1,0.1),labels=NA,tck=-0.04)
    sapply(seq(0.4,1,0.1),function(x){
      axis(2,at=x,line=-0.7,cex.axis=0.65,tick=F,las=2,
           labels=paste(round(x*100,0),"%",sep=""))
    })
    mtext(2,text="Singleton Proportion",line=1.5,cex=0.8)
}


#Organize table of genes from SV data
getSVdat <- function(dat,genes,prefix=NULL){
  #pLoF
  lof.any <- as.character(unlist(strsplit(dat$PROTEIN_CODING__LOF[which(!is.na(dat$PROTEIN_CODING__LOF))],split=",")))
  lof.del <- as.character(unlist(strsplit(dat$PROTEIN_CODING__LOF[which(!is.na(dat$PROTEIN_CODING__LOF)
                                                                        & dat$SVTYPE=="DEL")],split=",")))
  lof.other <- as.character(unlist(strsplit(dat$PROTEIN_CODING__LOF[which(!is.na(dat$PROTEIN_CODING__LOF)
                                                                          & dat$SVTYPE!="DEL")],split=",")))
  #Dup GC
  cg.dup <- as.character(unlist(strsplit(dat$PROTEIN_CODING__COPY_GAIN[which(!is.na(dat$PROTEIN_CODING__COPY_GAIN))],split=",")))
  #Dup pLoF
  plof.dup <- as.character(unlist(strsplit(dat$PROTEIN_CODING__DUP_LOF[which(!is.na(dat$PROTEIN_CODING__DUP_LOF))],split=",")))
  #Inversion span
  inv.span <- as.character(unlist(strsplit(dat$PROTEIN_CODING__INV_SPAN[which(!is.na(dat$PROTEIN_CODING__INV_SPAN))],split=",")))
  #Collect vector per gene
  res <- as.data.frame(t(sapply(genes,function(gene){
    g.lof.any <- length(which(lof.any==gene))
    g.lof.del <- length(which(lof.del==gene))
    g.lof.other <- length(which(lof.other==gene))
    g.cg.dup <- length(which(cg.dup==gene))
    g.plof.dup <- length(which(plof.dup==gene))
    g.inv.span <- length(which(inv.span==gene))
    g.out <- as.integer(c(g.lof.any,g.lof.del,g.lof.other,g.cg.dup,g.plof.dup,g.inv.span))
    g.out[which(is.na(g.out))] <- 0
    g.out <- c(gene,g.out)
    return(g.out)
  })))
  colnames(res) <- c("gene","lof.any","lof.del","lof.other","cg","plof","inv")
  if(!is.null(prefix)){
    colnames(res)[-1] <- paste(prefix,colnames(res)[-1],sep=".")
  }
  rownames(res) <- 1:nrow(res)
  return(res)
}
#Get merged table of SNV stats and SV stats by frequency bin
getSVdat.all <- function(dat,snv.data){
  #Gather data
  genes <- sort(unique(as.character(snv.data$gene)))
  all.counts <- getSVdat(dat,genes,prefix="all")
  common.counts <- getSVdat(dat[which(dat$AF>=0.01),],genes,prefix="common")
  rare.counts <- getSVdat(dat[which(dat$AF<0.01 & dat$AC>1),],genes,prefix="rare")
  singleton.counts <- getSVdat(dat[which(dat$AC==1),],genes,prefix="singleton")
  #Merge data
  merged <- merge(x=snv.data,y=all.counts,by="gene",all.x=T,sort=F)
  merged <- merge(x=merged,y=common.counts,by="gene",all.x=T,sort=F)
  merged <- merge(x=merged,y=rare.counts,by="gene",all.x=T,sort=F)
  merged <- merge(x=merged,y=singleton.counts,by="gene",all.x=T,sort=F)
  merged[,-c(which(colnames(merged) %in% c("gene","chrom")))] <- apply(merged[,-c(which(colnames(merged) %in% c("gene","chrom")))],2,as.numeric)
  return(merged)
}

#Stacked barplot of % of genes with at least one SV per category
barplot.genesFracByEffect <- function(gene.data){
  ngenes <- nrow(gene.data)
  #Get cumulative fractions per type
  plot.dat <- sapply(c("lof.any","lof.del","lof.other","plof","cg","inv"),function(effect){
    idxs <- grep(paste(".",effect,sep=""),colnames(gene.data),fixed=T)
    sing.hit <- length(which(gene.data[,idxs[4]]>0))/ngenes
    singOrRare.hit <- length(which(gene.data[,idxs[4]]>0 | gene.data[,idxs[3]]>0))/ngenes
    singOrRareOrCommon.hit <- length(which(gene.data[,idxs[4]]>0 | gene.data[,idxs[3]]>0 | gene.data[,idxs[2]]>0))/ngenes
    return(c(sing.hit,singOrRare.hit,singOrRareOrCommon.hit))
  })
  
  #Prepare plot area
  par(bty="n",mar=c(0.1,5,2.5,0.75))
  plot(x=c(0,1),y=c(0,-6),type="n",
       xaxt="n",xlab="",yaxt="n",ylab="")
  axis(2,at=-(0:5)-0.5,tick=F,las=2,cex.axis=0.8,line=-1.2,
       labels=c("pLoF (All)","pLoF (DEL Only)","pLoF (Not DEL)",
                "IED","CG","Inverted Gene"))
  axis(3,at=seq(0,1,0.2),labels=NA,tck=-0.04)
  sapply(seq(0,1,0.2),function(l){
    axis(3,at=l,tick=F,line=-0.6,cex.axis=0.8,
         labels=paste(100*l,"%",sep=""))
  })
  mtext(3,line=1.4,text="Pct. of All Autosomal Genes")
  
  #Plot bars
  bar.base.cols <- c(svtypes$color[which(svtypes$svtype=="DEL")],
                     svtypes$color[which(svtypes$svtype=="DEL")],
                     svtypes$color[which(svtypes$svtype=="DEL")],
                     svtypes$color[which(svtypes$svtype=="MCNV")],
                     svtypes$color[which(svtypes$svtype=="DUP")],
                     svtypes$color[which(svtypes$svtype=="INV")])
  sapply(1:ncol(plot.dat),function(i){
    rect(xleft=0,xright=1,ybottom=-i+0.8,ytop=-i+0.2,
         col="white",border="gray85")
    rect(xleft=c(0,plot.dat[,i]),xright=c(plot.dat[,i],1),
         ybottom=-i+0.2,ytop=-i+0.8,border=NA,bty="n",
         col=c(bar.base.cols[i],
               adjustcolor(bar.base.cols[i],alpha=0.5),
               adjustcolor(bar.base.cols[i],alpha=0.15),
               "gray98"))
    rect(xleft=0,xright=max(plot.dat[,i]),
         ybottom=-i+0.2,ytop=-i+0.8,col=NA)
    text(x=max(plot.dat[,i]),y=-i+0.45,cex=0.7,pos=4,
         labels=paste(format(round(100*max(plot.dat[,i]),1),nsmall=1),"%",sep=""))
  })
}

#Mini histogram of # of SV per gene by functional effect
minihist.byFunc <- function(gene.data,func,title=NULL,base.col,x.max=5){
  #Get hist data
  counts <- sapply(0:x.max,function(k){
    length(which(gene.data[,which(colnames(gene.data) == paste("all",func,sep="."))]==k))
  })
  counts <- c(counts,length(which(gene.data[,which(colnames(gene.data) == paste("all",func,sep="."))]>x.max)))
  counts <- counts/1000
  #Plot hist
  par(mar=c(2,3,1.5,0.5),bty="n")
  plot(x=c(0,x.max+2),y=c(0,max(counts)),type="n",
       xaxt="n",xlab="",yaxt="n",ylab="")
  rect(xleft=(0:(x.max+1))+0.1,xright=c(0:(x.max+1))+0.9,
       ybottom=0,ytop=counts,col=base.col)
  sapply((0:x.max),function(k){
    axis(1,at=k+0.5,labels=k,cex.axis=0.8,tick=F,line=-1.3)
  })
  axis(1,at=x.max+1.5,labels=paste(">",x.max,sep=""),cex.axis=0.8,tick=F,line=-1.3)
  
  mtext(1,text="SV per Gene",line=0.75,cex=0.85)
  axis(2,at=axTicks(2),labels=NA)
  axis(2,at=axTicks(2),tick=F,labels=paste(axTicks(2),"k",sep=""),
       line=-0.4,las=2,cex.axis=0.8)
  mtext(2,text="Genes",line=1.75,cex=0.9)
  mtext(3,text=title)
}


##############################
###COUNT PER GENE CORRELATIONS
##############################
#Run permutation test to gather expected correlation between two functional classes
permute.perGeneCounts <- function(gene.data,x.category,y.category,times=1000,max.ax=3){
  #Get actual values
  x <- gene.data[,which(colnames(gene.data)==x.category)]
  y <- gene.data[,which(colnames(gene.data)==y.category)]
  #Get shuffled values
  shuf.res <- lapply(1:times,function(i){
    x.shuf <- sample(x,size=length(x),replace=F)
    y.shuf <- sample(y,size=length(y),replace=F)
    shuf.dat <- sapply(0:(max.ax+1),function(x.ct){
      sapply(0:(max.ax+1),function(y.ct){
        if(x.ct>max.ax & y.ct>max.ax){
          length(which(x.shuf>max.ax &
                         y.shuf>max.ax))
        }else if(x.ct>max.ax & y.ct<=max.ax){
          length(which(x.shuf>max.ax &
                         y.shuf==y.ct))
        }else if(x.ct<=max.ax & y.ct>max.ax){
          length(which(x.shuf==x.ct &
                         y.shuf>max.ax))
        }else{
          length(which(x.shuf==x.ct &
                         y.shuf==y.ct))
        }
      })
    })
  })
  #Get means of shuffled values
  shuf.dat.means <- sapply(0:(max.ax+1),function(x.ct){
    sapply(0:(max.ax+1),function(y.ct){
      mean(sapply(shuf.res,function(df){
        return(df[x.ct+1,y.ct+1])
      }))
    })
  })
  return(shuf.dat.means)
}
#Generic heatmap of counts per gene for two functional classes
plot.perGeneHeatmaps <- function(gene.data,x.category,y.category,max.ax=3,max.col=0.2,
                                 x.label=NULL,y.label=NULL){
  #Get plot data
  x.idx <- which(colnames(gene.data)==x.category)
  y.idx <- which(colnames(gene.data)==y.category)
  plot.dat <- sapply(0:(max.ax+1),function(x.ct){
    sapply(0:(max.ax+1),function(y.ct){
      if(x.ct>max.ax & y.ct>max.ax){
        length(which(gene.data[,x.idx]>max.ax &
                       gene.data[,y.idx]>max.ax))
      }else if(x.ct>max.ax & y.ct<=max.ax){
        length(which(gene.data[,x.idx]>max.ax &
                       gene.data[,y.idx]==y.ct))
      }else if(x.ct<=max.ax & y.ct>max.ax){
        length(which(gene.data[,x.idx]==x.ct &
                       gene.data[,y.idx]>max.ax))
      }else{
        length(which(gene.data[,x.idx]==x.ct &
                       gene.data[,y.idx]==y.ct))
      }
    })
  })
  plot.dat <- plot.dat/sum(plot.dat)
  plot.dat[which(plot.dat>max.col)] <- max.col
  plot.dat <- floor(1000*plot.dat)
  #Prep plot area
  par(mar=c(3,3,0.5,0.5))
  plot(x=c(0,max.ax+2),y=c(0,max.ax+2),type="n",
       xlab="",ylab="",xaxt="n",yaxt="n",xaxs="i",yaxs="i")
  mtext(1,line=1.5,text=x.label,cex=0.9)
  axis(1,at=seq(0,max.ax,1)+0.5,tick=F,line=-0.8,labels=seq(0,max.ax,1))
  axis(1,at=max.ax+1.5,tick=F,line=-0.8,
       labels=paste(">",max.ax,sep=""))
  mtext(2,line=1.5,text=y.label,cex=0.9)
  axis(2,at=seq(0,max.ax,1)+0.5,tick=F,line=-0.8,labels=seq(0,max.ax,1),las=2)
  axis(2,at=max.ax+1.5,tick=F,line=-0.8,las=2,
       labels=paste(">",max.ax,sep=""))
  #Add heatmap
  col.pal <- colorRampPalette(c("#440154","#365C8C","#25A584","#FDE725"))(length(seq(0,ceiling(1000*max.col),1)))
  sapply(0:(max.ax+1),function(x){
    sapply(0:(max.ax+1),function(y){
      rect(xleft=x,xright=x+1,
           ybottom=y,ytop=y+1,
           border=NA,col=col.pal[plot.dat[y+1,x+1]+1])
    })
  })
}
#Key for per gene heatmaps
plot.perGeneHeatmapsKey <- function(max.col=0.2){
  col.pal <- colorRampPalette(c("#440154","#365C8C","#25A584","#FDE725"))(length(seq(0,ceiling(1000*max.col),1)))
  par(mar=rep(0.1,4))
  plot(x=c(0,1),y=c(0,length(col.pal)),type="n",
       xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab="",ylab="")
  rect(xleft=0,xright=1,ybottom=(1:length(col.pal))-1,ytop=1:length(col.pal),
       border=col.pal,col=col.pal)
  box()
}
#Scaled dotplots of counts per gene for two functional classes vs permuted expectation
plot.perGeneDotplots <- function(gene.data,x.category,y.category,max.ax=3,max.size=5000,
                                 x.label=NULL,y.label=NULL,perm.times=1000){
  #Get empirical data
  x.idx <- which(colnames(gene.data)==x.category)
  y.idx <- which(colnames(gene.data)==y.category)
  plot.dat <- sapply(0:(max.ax+1),function(x.ct){
    sapply(0:(max.ax+1),function(y.ct){
      if(x.ct>max.ax & y.ct>max.ax){
        length(which(gene.data[,x.idx]>max.ax &
                       gene.data[,y.idx]>max.ax))
      }else if(x.ct>max.ax & y.ct<=max.ax){
        length(which(gene.data[,x.idx]>max.ax &
                       gene.data[,y.idx]==y.ct))
      }else if(x.ct<=max.ax & y.ct>max.ax){
        length(which(gene.data[,x.idx]==x.ct &
                       gene.data[,y.idx]>max.ax))
      }else{
        length(which(gene.data[,x.idx]==x.ct &
                       gene.data[,y.idx]==y.ct))
      }
    })
  })
  plot.dat[which(plot.dat>max.size)] <- max.size
  plot.dat <- sqrt(plot.dat)/(2*sqrt(max.size))
  #Get permuted expectations
  plot.exp <- permute.perGeneCounts(gene.data,x.category,y.category,times=perm.times,max.ax)
  plot.exp[which(plot.exp>max.size)] <- max.size
  plot.exp <- sqrt(plot.exp)/(2*sqrt(max.size))
  #Prep plot area
  par(mar=c(3,3,0.5,0.5),bty="n")
  plot(x=c(0,max.ax+2),y=c(0,max.ax+2),type="n",
       xlab="",ylab="",xaxt="n",yaxt="n")
  abline(h=par("usr")[1],v=par("usr")[3])
  mtext(1,line=1.5,text=x.label,cex=0.9)
  axis(1,at=seq(0,max.ax,1)+0.5,tick=F,line=-0.8,labels=seq(0,max.ax,1))
  axis(1,at=max.ax+1.5,tick=F,line=-0.8,
       labels=paste(">",max.ax,sep=""))
  mtext(2,line=1.5,text=y.label,cex=0.9)
  axis(2,at=seq(0,max.ax,1)+0.5,tick=F,line=-0.8,labels=seq(0,max.ax,1),las=2)
  axis(2,at=max.ax+1.5,tick=F,line=-0.8,las=2,
       labels=paste(">",max.ax,sep=""))
  #Add circles
  sapply(0:(max.ax+1),function(x){
    sapply(0:(max.ax+1),function(y){
      #Grey shaded circle for expectation
      draw.circle(x=x+0.5,y=y+0.5,radius=plot.exp[x+1,y+1],
                  border=NA,col="gray80")
      #Black open circle for empirical observation
      draw.circle(x=x+0.5,y=y+0.5,radius=plot.dat[y+1,x+1],
                  col=NA,lwd=2)
    })
  })
}
#Scaled dotplot key
plot.perGeneDotplotKey <- function(max.size=5000,marks=c(10,100,1000,5000)){
  marks <- sqrt(marks)/(2*sqrt(max.size))
  par(mar=rep(0.1,4),bty="n")
  plot(x=c(0,1),y=c(0,length(marks)),type="n",
       xaxt="n",yaxt="n",xlab="",ylab="")
  sapply(1:length(marks),function(i){
    draw.circle(x=0.5,y=(i/1.5)-0.5,
                radius=marks[i],
                lwd=2,col="gray80")
  })
}


##################################
###CODE FOR CONSTRAINT COMPARISONS
##################################
#Prep covariates matrics for regression analysis
prepCovariates <- function(gene.data,y){
  cov <- data.frame("gene"=gene.data$gene,
                    "y.SV_count"=gene.data[which(colnames(gene.data)==y)],
                    "gene_length"=scale(log10(as.numeric(gene.data$gene_length)),center=T,scale=T),
                    "exon_count"=scale(log10(as.numeric(gene.data$exon_count)),center=T,scale=T),
                    # "exon_min"=scale(log10(as.numeric(gene.data$exon_min)),center=T,scale=T),
                    # "exon_max"=scale(log10(as.numeric(gene.data$exon_min)),center=T,scale=T),
                    # "exon_mean"=scale(log10(as.numeric(gene.data$exon_mean)),center=T,scale=T),
                    "exon_median"=scale(log10(as.numeric(gene.data$exon_median)),center=T,scale=T),
                    # "exon_harm.mean"=scale(log10(as.numeric(gene.data$exon_harm.mean)),center=T,scale=T),
                    "exon_sum"=scale(log10(as.numeric(gene.data$exon_mean)),center=T,scale=T),
                    "intron_count"=scale(log10(as.numeric(gene.data$intron_count)),center=T,scale=T),
                    # "intron_min"=scale(log10(as.numeric(gene.data$intron_min)),center=T,scale=T),
                    # "intron_max"=scale(log10(as.numeric(gene.data$intron_min)),center=T,scale=T),
                    # "intron_mean"=scale(log10(as.numeric(gene.data$intron_mean)),center=T,scale=T),
                    "intron_median"=scale(log10(as.numeric(gene.data$intron_median)),center=T,scale=T),
                    # "intron_harm.mean"=scale(log10(as.numeric(gene.data$intron_harm.mean)),center=T,scale=T),
                    "intron_sum"=scale(log10(as.numeric(gene.data$intron_mean)),center=T,scale=T),
                    "segdup"=round(gene.data$segdup,0))
  colnames(cov)[2] <- "y.SV_count"
  chrom.dummy.mat <- as.data.frame(sapply(unique(gene.data$chrom),function(chr){
    v <- rep(0,times=nrow(gene.data))
    v[which(gene.data$chrom==chr)] <- 1
    return(v)
  }))
  cov <- cbind(cov,chrom.dummy.mat)
  # rownames(cov) <- gene.data$gene
  return(cov)
}

#Fit model for predicting # of rare functional SV
fitConstraintModel <- function(cov,train.genes=which(gene.data$ptv_oe_dec>=5 & gene.data$ptv_oe_dec<=9)){
  #Fit Poisson model
  # outlier.gene.cutoff <- quantile(cov$y.SV_count,0.99)
  # outlier.gene.cutoff <- 3
  # train.genes <- intersect(train.genes,which(cov$y.SV_count<=outlier.gene.cutoff))
  # glm.fit <- glm(y.SV_count ~ ., data=cov[train.genes,-1], family="poisson")
  glm.fit <- glm.nb(y.SV_count ~ ., data=cov[train.genes,-1])
  #Evaluate variance explained
  # af <- anova(glm.fit)
  # af$PctExp <- 100*af$`Sum Sq`/sum(af$`Sum Sq`)
  # pct.expl <- sum(af$PctExp[-length(af$PctExp)])
  #Apply model to all genes
  glm.fit.vals <- predict.glm(glm.fit,newdata=cov[,-(1:2)],type="response")
  #Return output data frame
  res <- data.frame("gene"=cov$gene,
                    "SV_count_raw"=cov$y.SV_count,
                    "SV_count_glm_exp"=glm.fit.vals)
  return(res)
}

#Make summed dotplots
plotSummedDots <- function(deciles,SV_count_obs,SV_count_exp,
                           color,title=NULL,ymax=NULL,
                           xlabel="SNV pLoF Constraint Percentile",
                           ax.labels=TRUE,cex.labels=0.8,
                           tck=NULL,yline=-0.4,conf.int=F,parmar=c(2,3.25,1.75,2)){
  require(zoo,quietly=T)
  d <- sort(unique(as.numeric(deciles)))
  #Compute summed obs/exp
  means <- sapply(d,function(i){
    exp.sum <- sum(SV_count_exp[which(deciles==i)])
    obs.sum <- sum(SV_count_obs[which(deciles==i)])
    return(obs.sum/exp.sum)
  })
  means[which(means<0)] <- NA
  #Compute CI
  if(conf.int==T){
    cis <- t(sapply(d,function(i){
      oe.vals <- SV_count_obs[which(deciles==i)]/SV_count_exp[which(deciles==i)]
      get.mean <- function(vals,indices){mean(vals[indices])}
      set.seed(i)
      boot.obj <- boot(data=oe.vals,statistic=get.mean,R=1000)
      ci <- boot.ci(boot.obj,conf=0.9,type="basic")$basic[4:5]
      return(ci)
    }))
  }
  
  #Prep plot area
  if(is.null(ymax)){
    ymax <- max(means,na.rm=T)
  }
  par(mar=parmar)
  plot(x=c(0,length(d)),y=c(0,ymax),type="n",
       xaxt="n",yaxt="n",xlab="",ylab="")
  # axis(1,at=seq(0,100,10),labels=NA)
  sapply(seq(0,100,20),function(p){
    axis(1,at=p,tick=F,cex.axis=cex.labels,line=-0.9,
         labels=bquote(.(p)^'th'))
  })
  # text(x=(1:10)-0.75,y=par('usr')[3]-(0.115*(par("usr")[4]-par("usr")[3])),
  #      labels=paste(seq(0,90,10),"-",seq(10,100,10),"%",sep=""),
  #      cex=0.7,srt=45)
  axis(2,at=axTicks(2),labels=NA,tck=tck)
  if(ax.labels==T){
    mtext(1,line=0.9,text=xlabel,cex=cex.labels)
    mtext(2,line=2.15,text="Rare SV Obs/Exp",cex=cex.labels)
  }
  axis(2,at=axTicks(2),line=yline,cex.axis=cex.labels,
       labels=paste(100*round(axTicks(2),2),"%",sep=""),tick=F,las=2)
  mtext(3,text=title,line=0.2,cex=0.9)
  abline(h=1,lwd=2,col="gray80")
  
  #Add points & rolling mean/ci
  # points(x=d-0.5,y=means,col=color,type="l")
  if(conf.int==T){
    polygon(x=c(d-0.5,rev(d-0.5)),
            y=c(rollapply(cis[,1],21,mean,na.rm=T,partial=T),
                rev(rollapply(cis[,2],21,mean,na.rm=T,partial=T))),
            border=NA,bty="n",col=adjustcolor(color,alpha=0.1))
  }
  points(x=d-0.5,y=means,pch=21,col=color,cex=0.4)
  points(x=d-0.5,y=rollapply(means,21,mean,na.rm=T,partial=T),
         lwd=2,type="l",col=color)
  
  #Add correlation coefficient
  cor.res <- cor.test(x=d-0.5,y=means,method="spearman")
  text(x=par("usr")[1]-(0.025*(par("usr")[2]-par("usr")[1])),
       y=par("usr")[3]+(0.885*(par("usr")[4]-par("usr")[3])),
       pos=4,cex=cex.labels,
       labels=bquote(rho == .(format(round(cor.res$estimate,2),nsmall=2))))
  # cor.res <- cor.test(x=d-0.5,y=means,method="pearson")
  cor.p <- cor.res$p.value
  if(cor.p>10^-100){
    cor.p <- format(cor.p,scientific=T)
    cor.p.base <- format(as.numeric(strsplit(cor.p,split="e")[[1]][1]),nsmall=2)
    cor.p.exp <- format(as.numeric(strsplit(cor.p,split="e")[[1]][2]),nsmall=0)
    text(x=par("usr")[2]+(0.025*(par("usr")[2]-par("usr")[1])),
         y=par("usr")[3]+(0.075*(par("usr")[4]-par("usr")[3])),
         pos=2,cex=cex.labels,
         labels=bquote(italic(P) == .(format(round(as.numeric(cor.p.base),2),nsmall=2))*"x"*10^.(cor.p.exp)))
  }else{
    text(x=par("usr")[2]+(0.025*(par("usr")[2]-par("usr")[1])),
         y=par("usr")[3]+(0.075*(par("usr")[4]-par("usr")[3])),
         pos=2,cex=cex.labels,
         labels=bquote(italic(P) < 10^-100))
  }
  box()
}

#Wrapper for all constraint analyses for a single class of SV
constraintModelWrapper <- function(gene.data,y,snv.class="ptv",color,title,
                                   ymax=NULL,return=FALSE,ax.labels=TRUE,
                                   conf.int=F,cex.labels=0.8,
                                   tck=NULL,yline=-0.4,parmar=c(2,3.25,1.75,2)){
  if(snv.class=="ptv"){
    # train.genes <- which(gene.data$pli<=0.1)
    train.genes <- which(gene.data$ptv_oe_dec>=5 & gene.data$ptv_oe_dec<=10)
    deciles <- gene.data$ptv_oe_cent
    xlabel <- "SNV pLoF Constraint Percentile"
  }else{
    train.genes <- which(gene.data$mis_oe_dec>=5 & gene.data$mis_oe_dec<=10)
    deciles <- gene.data$mis_oe_cent
    xlabel <- "Missense Constraint Percentile"
  }
  cov <- prepCovariates(gene.data=gene.data,y=y)
  res <- fitConstraintModel(cov=cov,train.genes=train.genes)
  print(paste("VARIANCE EXPLAINED: ",100*(cor(res[,2],res[,3],use="complete.obs")^2),"%",sep=""))
  plotSummedDots(deciles=deciles,
                 SV_count_obs=res$SV_count_raw,
                 SV_count_exp=res$SV_count_glm_exp,
                 color=color,title=title,ymax,
                 xlabel=xlabel,
                 ax.labels=ax.labels,
                 cex.labels=cex.labels,
                 conf.int=conf.int,
                 parmar=parmar,
                 tck=tck,yline=yline)
  if(return==T){
    return(res)
  }
}

# #Plot per-gene SV-only constraint data (boxplot)
# 
# #Wrapper for SV-only constraint analysis (per-gene significance)
# perGeneSVConstraintWrapper <- function(gene.data,y,color,return=F){
#   #Gather data
#   train.genes <- which(gene.data$ptv_oe_dec>=5 & gene.data$ptv_oe<10)
#   cov <- prepCovariates(gene.data=gene.data,y=y)
#   res <- fitConstraintModel(cov=cov,train.genes=train.genes)
#   res$oe <- res[,2]/res[,3]
#   res$p <- apply(res[,2:3],1,function(vals){
#     obs <- as.numeric(vals[1])
#     exp <- as.numeric(vals[2])
#     ppois(q=obs,lambda=exp)
#   })
#   res$q <- p.adjust(res$p,method="fdr")
#   res$bonf <- p.adjust(res$p,method="bonf")
#   #Plot data
#   
# }

#Mini histogram of # of SV per gene by constraint category
minihist.byConstraint <- function(gene.data,min.pli=0,max.pli=1,
                                  title=NULL,x.max=10){
  #Get hist data
  counts <- sapply(0:x.max,function(k){
    length(which(gene.data[,which(colnames(gene.data) == "rare.lof.any")]==k
                 & gene.data$pli>=min.pli & gene.data$pli<=max.pli))
  })
  counts <- c(counts,length(which(gene.data[,which(colnames(gene.data) == "rare.lof.any")]>x.max
                                  & gene.data$pli>=min.pli & gene.data$pli<=max.pli)))
  counts <- counts/1000
  #Plot hist
  base.col <- "gray70"
  par(mar=c(2,3,2.5,0.5),bty="n")
  plot(x=c(0,x.max+2),y=c(0,max(counts)),type="n",
       xaxt="n",xlab="",yaxt="n",ylab="")
  rect(xleft=(0:(x.max+1))+0.1,xright=c(0:(x.max+1))+0.9,
       ybottom=0,ytop=counts,col=base.col)
  sapply((0:x.max),function(k){
    axis(1,at=k+0.5,labels=k,cex.axis=0.8,tick=F,line=-1.3)
  })
  axis(1,at=x.max+1.5,labels=paste("> ",x.max,sep=""),cex.axis=0.8,tick=F,line=-1.3)
  mtext(1,text="Rare pLoF SV per Gene",line=0.75,cex=0.75)
  axis(2,at=axTicks(2),labels=NA)
  axis(2,at=axTicks(2),tick=F,labels=paste(axTicks(2),"k",sep=""),
       line=-0.4,las=2,cex.axis=0.8)
  mtext(2,text="Genes",line=1.75,cex=0.75)
  mtext(3,text=title,cex=0.8)
  text(x=mean(c(par("usr")[1],par("usr")[2])),
       y=par("usr")[4],pos=1,cex=0.7,
       labels=paste("n = ",prettyNum(length(which(gene.data$pli>=min.pli & gene.data$pli<=max.pli)),
                                     big.mark=","),
                    " genes",sep=""))
  text(x=mean(c(par("usr")[1],par("usr")[2])),
       y=0.9*par("usr")[4],pos=1,cex=0.7,
       labels=paste("~",format(round(mean(gene.data$rare.lof.any[which(gene.data$pli>=min.pli & gene.data$pli<=max.pli)]),2),nsmall=2),
                    " pLoF SV / gene",sep=""))
}


##########################
###HUMAN KNOCKOUT ANALYSIS
##########################
#Gather all data for KO analysis
getKOdata <- function(dat,gene.data){
  #Subset to biallelic, autosomal, resolved pLoF SV
  subdat <- dat[which(!(dat$chrom %in% c("X","Y")) & !is.na(dat$PROTEIN_CODING__LOF)),]
  rows.to.exclude <- grep("MULTIALLELIC",subdat$FILTER,fixed=T)
  if(length(rows.to.exclude)>0){
    subdat <- subdat[-rows.to.exclude,]
  }
  
  #Helper function to return total number of nonredundant genes and total number of variants
  countKO <- function(subdat,gene.data){
    n.variants <- nrow(subdat)
    f.variants <- n.variants/nrow(dat[which(!(dat$chrom %in% c("X","Y"))),])
    genes <- unique(sort(unlist(strsplit(subdat$PROTEIN_CODING__LOF,split=","))))
    genes <- genes[which(genes %in% gene.data$gene)]
    n.genes <- length(genes)
    f.genes <- n.genes/nrow(gene.data[which(!(gene.data$chrom %in% c("chrX","chrY"))),])
    return(c(n.variants,f.variants,n.genes,f.genes))
  }
  
  #Gather counts
  a <- countKO(subdat,gene.data)
  subdat <- subdat[which(subdat$N_HOMALT>0),]
  b <- countKO(subdat,gene.data)
  subdat <- subdat[which(subdat$POPMAX_AF<0.01),]
  c <- countKO(subdat,gene.data)
  subdat <- subdat[which(subdat$N_HOMALT>1),]
  d <- countKO(subdat,gene.data)
  
  #Prep & return output data frame
  out.res <- t(data.frame(a,b,c,d))
  colnames(out.res) <- c("n.variants","f.variants","n.genes","f.genes")
  rownames(out.res) <- c("All pLoF SV","pLoF SV with > 0 Homozygotes",
                         "Rare (AF < 1%) pLoF SV with > 0 Homozygotes",
                         "Rare (AF < 1%) pLoF SV with > 1 Homozygote")
  return(out.res)
}
#Human KO plot
plotKO <- function(KO.data){
  #Helper function for generic barplot
  genericBarplot <- function(vals,label.format="count"){
    par(mar=c(0.5,0.5,2.5,1.75))
    vals <- rev(vals)
    plot(x=c(0,1.1*max(vals)),y=c(0,length(vals)),type="n",
         xlab="",xaxt="n",ylab="",yaxt="n")
    rect(xleft=0,xright=vals,ybottom=(1:length(vals))-0.8,ytop=(1:length(vals))-0.2,
         col=rev(c(adjustcolor(svtypes$color[which(svtypes$svtype=="DEL")],alpha=0.4),
                   adjustcolor(svtypes$color[which(svtypes$svtype=="DEL")],alpha=0.6),
                   adjustcolor(svtypes$color[which(svtypes$svtype=="DEL")],alpha=0.8),
                   svtypes$color[which(svtypes$svtype=="DEL")])),
         lwd=0.75)
    abline(v=0,lwd=0.8)
    if(length(axTicks(3))>5){
      x.at <- axTicks(3)[seq(1,length(axTicks(3)),2)]
      x.at <- c(x.at,x.at[length(x.at)]+(x.at[2]-x.at[1]))
    }else{
      x.at <- axTicks(3)
    }
    if(label.format=="count"){
      text(x=vals-(0.05*(par("usr")[2]-par("usr")[1])),y=(1:length(vals))-0.55,pos=4,cex=0.7,
           labels=prettyNum(vals,big.mark=","),xpd=T)
      axis(3,at=x.at,labels=NA,tck=-0.05,lwd=0.75)
      axis(3,at=x.at,tick=F,cex.axis=0.7,line=-0.6,
           labels=paste(format(round(x.at/1000,1),nsmall=1),"k",sep=""))
    }else{
      text(x=vals-(0.05*(par("usr")[2]-par("usr")[1])),y=(1:length(vals))-0.55,pos=4,cex=0.7,
           labels=paste(format(round(100*vals,1),nsmall=1),"%",sep=""),xpd=T)
      axis(3,at=x.at,labels=NA,tck=-0.05,lwd=0.75)
      axis(3,at=x.at,tick=F,cex.axis=0.7,line=-0.6,
           labels=paste(round(100*x.at,1),"%",sep=""))
    }
  }
  
  #Prep plot layout
  layout(matrix(1:5,nrow=1,byrow=T),
         widths=c(2,1,1,1,1))
  #Panel 1: labels
  par(mar=c(0.5,0.5,2.5,1.75),bty="n")
  plot(x=c(0,1),y=c(0,nrow(KO.data)),type="n",
       xlab="",xaxt="n",ylab="",yaxt="n")
  text(x=1.2*par("usr")[2],y=(1:nrow(KO.data))-0.55,pos=2,cex=0.75,
       labels=rev(rownames(KO.data)),xpd=T)
  #Panel 2: count of variants
  genericBarplot(KO.data[,1],label.format="count")
  mtext(3,line=1.25,text="SV (Count)",cex=0.7)
  #Panel 3: pct of variants
  genericBarplot(KO.data[,2],label.format="pct")
  mtext(3,line=1.25,text="SV (Pct.)",cex=0.7)
  #Panel 4: count of genes
  genericBarplot(KO.data[,3],label.format="count")
  mtext(3,line=1.25,text="Genes (Count)",cex=0.7)
  #Panel 5: pct of genes
  genericBarplot(KO.data[,4],label.format="pct")
  mtext(3,line=1.25,text="Genes (Pct.)",cex=0.7)
}
#Gather genes from 3-way comparison of existing KO studies
intersect.KOstudies <- function(SV.kos,narasimhan.kos,saleheen.kos){
  #Get intersections
  a.set <- SV.kos
  b.set <- narasimhan.kos
  c.set <- saleheen.kos
  a <- a.set[which(!(a.set %in% b.set) & !(a.set %in% c.set))]
  b <- b.set[which(!(b.set %in% a.set) & !(b.set %in% c.set))]
  c <- c.set[which(!(c.set %in% a.set) & !(c.set %in% b.set))]
  ab <- a.set[which((a.set %in% b.set) & !(a.set %in% c.set))]  
  ac <- a.set[which(!(a.set %in% b.set) & (a.set %in% c.set))]
  bc <- b.set[which(!(b.set %in% a.set) & (b.set %in% c.set))]
  abc <- a.set[which((a.set %in% b.set) & (a.set %in% c.set))]
  #Return list
  return(list("sv"=a,"nara"=b,"sale"=c,
              "sv_nara"=ab,"sv_sale"=ac,"nara_sale"=bc,
              "sv_nara_sale"=abc))
}
#Upset plot of three-way comparison vs existing KO studies
plot.KOstudyComparison <- function(ko.intersections){
  #Prep layout
  layout(matrix(1:2,nrow=2),heights=c(3,1))
  #Top panel: upset bars
  h <- rev(sort(unlist(lapply(ko.intersections,length))))
  par(mar=c(0.1,12,0.5,0.5),bty="n")
  plot(x=c(0,length(ko.intersections)),y=c(0,1.05*max(h)),type="n",
       xlab="",xaxt="n",ylab="",yaxt="n")
  col.bars <- rep("gray70",times=length(h))
  col.bars[grep("sv",names(h))] <- svtypes$color[which(svtypes$svtype=="DEL")]
  rect(xleft=(1:length(h))-0.8,xright=(1:length(h))-0.2,
       ybottom=0,ytop=h,col=col.bars)
  text(x=(1:length(h))-0.5,y=h-(0.025*(par("usr")[4]-par("usr")[3])),pos=3,cex=0.7,
       labels=prettyNum(h,big.mark=","))
  axis(2,at=axTicks(2),labels=NA)
  axis(2,at=axTicks(2),tick=F,line=-0.4,cex.axis=0.7,las=2,
       labels=prettyNum(axTicks(2),big.mark=","))
  mtext(2,line=2,text="Genes with\nHomozygous pLoF")
  
  #Bottom panel: points & study labels
  par(mar=c(0.1,12,0.1,0.5),bty="n")
  plot(x=c(0,length(h)),y=c(0,3),type="n",
       xlab="",xaxt="n",ylab="",yaxt="n")
  point.cex <- 1.25
  sapply(1:3,function(y){
    points(x=(1:length(h))-0.5,y=rep(y-0.5,length(h)),
           pch=19,col="gray95",cex=point.cex)
  })
  studies <- rev(c("sv","nara","sale"))
  seg.lwd <- 2
  segments(x0=((1:length(h))-0.5)[grep("sv_nara",names(h))],
           x1=((1:length(h))-0.5)[grep("sv_nara",names(h))],
           y0=rep(2-0.5,length(h))[grep("sv_nara",names(h))],
           y1=rep(3-0.5,length(h))[grep("sv_nara",names(h))],
           lwd=seg.lwd)
  segments(x0=((1:length(h))-0.5)[grep("nara_sale",names(h))],
           x1=((1:length(h))-0.5)[grep("nara_sale",names(h))],
           y0=rep(1-0.5,length(h))[grep("nara_sale",names(h))],
           y1=rep(2-0.5,length(h))[grep("nara_sale",names(h))],
           lwd=seg.lwd)
  segments(x0=((1:length(h))-0.5)[intersect(grep("sv",names(h)),grep("sale",names(h)))],
           x1=((1:length(h))-0.5)[intersect(grep("sv",names(h)),grep("sale",names(h)))],
           y0=rep(1-0.5,length(h))[intersect(grep("sv",names(h)),grep("sale",names(h)))],
           y1=rep(3-0.5,length(h))[intersect(grep("sv",names(h)),grep("sale",names(h)))],
           lwd=seg.lwd)
  sapply(1:3,function(i){
    if(studies[i]=="sv"){
      pt.col <- svtypes$color[which(svtypes$svtype=="DEL")]
    }else{
      pt.col <- "black"
    }
    points(x=((1:length(h))-0.5)[grep(studies[i],names(h))],
           y=rep(i-0.5,length(h))[grep(studies[i],names(h))],
           pch=21,cex=point.cex,bg=pt.col)
  })
  axis(2,at=(1:3)-0.5,tick=F,line=-0.9,las=2,cex.axis=0.75,
       labels=c(expression("Saleheen"~italic("et al., Nature")~"(2017)"),
                expression("Narasimhan"~italic("et al., Science")~"(2016)"),
                "gnomAD-SV (this study)"))
}
#Plot boxplots of PTV or missense obs/exp across a list of gene sets
boxplots.oe <- function(genesets.list,gene.data,snv.category="ptv",base.colors,ymax=2){
  #Gather oe data
  oe.idx <- which(colnames(gene.data)==paste(snv.category,"_oe",sep=""))
  oe.vals <- lapply(genesets.list,function(genes){
    oe.vals <- gene.data[which(gene.data$gene %in% genes),oe.idx]
    oe.vals <- oe.vals[which(!is.na(oe.vals))]
    return(oe.vals)
  })
  
  #Prep plot area
  if(snv.category=="ptv"){
    label="pLoF SNV\nObs/Exp"
  }else{
    label="Missense\nObs/Exp"
  }
  par(mar=c(0.5,4.5,0.5,0.5),bty="n")
  plot(x=c(0.75,length(oe.vals)+0.25),y=c(0,ymax),type="n",
       xlab="",xaxt="n",ylab="",yaxt="n")
  abline(h=median(oe.vals[[1]]),lwd=2,col=base.colors[1])
  boxplot(oe.vals,add=T,outline=F,staplewex=0,lty=1,col=base.colors,
          xaxt="n",yaxt="n")
  axis(2,at=seq(0,2,0.5),labels=NA)
  axis(2,at=seq(0,2,0.5),tick=F,las=2,cex.axis=0.7,line=-0.4,
       labels=paste(seq(0,200,50),"%",sep=""))
  mtext(2,line=2.25,text=label,cex=1.2)
}
boxplots.size <- function(genesets.list,gene.data,base.colors){
  #Gather oe data
  size.idx <- which(colnames(gene.data)=="exon_sum")
  size.vals <- lapply(genesets.list,function(genes){
    size.vals <- gene.data[which(gene.data$gene %in% genes),size.idx]
    size.vals <- size.vals[which(!is.na(size.vals))]
    return(size.vals)
  })
  
  #Prep plot area
  label <- "Gene CDS\nLength (kb)"
  ymax <- 1.02*max(boxplot(size.vals,plot=F)$stats,na.rm=T)
  par(mar=c(0.5,4.5,0.5,0.5),bty="n")
  plot(x=c(0.75,length(size.vals)+0.25),y=c(0,ymax),type="n",
       xlab="",xaxt="n",ylab="",yaxt="n")
  abline(h=median(size.vals[[1]]),lwd=2,col=base.colors[1])
  boxplot(size.vals,add=T,outline=F,staplewex=0,lty=1,col=base.colors,
          xaxt="n",yaxt="n")
  axis(2,at=axTicks(2),labels=NA)
  axis(2,at=axTicks(2),tick=F,las=2,cex.axis=0.7,line=-0.4,
       labels=round(axTicks(2)/1000,0))
  mtext(2,line=2.25,text=label,cex=1.2)
}
getGeneListFraction <- function(comparison.set,genesets.list){
  #Helper: compute fraction of genes in list X that are present in list Y
  XinY <- function(X,Y){
    length(which(X %in% Y))/length(X)
  }
  vals <- unlist(lapply(genesets.list,XinY,Y=comparison.set))
  #Normalize to all genes (assumes all genes is first list
  vals <- vals/vals[1]
  return(vals)
}
getGeneListFractionMultiple <- function(comparison.sets.list,genesets.list,
                                        comparison.set.names,geneset.names){
  res <- t(do.call("rbind", lapply(comparison.sets.list,getGeneListFraction,genesets.list=genesets.list)))
  colnames(res) <- comparison.set.names
  rownames(res) <- geneset.names
  return(res)
}
#Generate dotplots of KO gene set enrichments
dotplots.KOenrichments <- function(KO.enrichment.data,colors){
  #Round down enrichments to max of 2
  KO.enrichment.data[which(KO.enrichment.data>2)] <- 2
  #Prep plot area
  par(mar=c(0.5,7,2.5,1),bty="n",xpd=F)
  plot(x=c(0,2),y=c(0,-ncol(KO.enrichment.data)),type="n",
       xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")
  #Add bars
  buffer <- 0.2
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=-(1:ncol(KO.enrichment.data))+buffer,
       ytop=-(1:ncol(KO.enrichment.data))+1-buffer,
       border="gray90",col="gray95")
  abline(v=1)
  sapply(1:ncol(KO.enrichment.data),function(i){
    top <- -i+1-buffer; bottom <- -i+buffer
    breaks <- rev(seq(bottom,top,abs(top-bottom)/nrow(KO.enrichment.data)))
    # rect(xleft=0,xright=KO.enrichment.data[,i],
    #      ybottom=breaks[-length(breaks)],ytop=breaks[-1],
    #      col=colors)
    points(x=KO.enrichment.data[,i],
           y=(breaks[-length(breaks)]+breaks[-1])/2,
           bg=colors,pch=21)
    axis(2,at=-i+0.5,tick=F,line=-0.8,las=2,cex.axis=0.8,
         labels=colnames(KO.enrichment.data)[i])
    axis(4,at=c(top,bottom),labels=NA,tck=0.02)
    axis(2,at=c(top,bottom),labels=NA,tck=0.02)
  })
  axis(3,at=seq(0,2,0.5),labels=NA)
  axis(3,at=seq(0,2,0.5),tick=F,labels=paste(seq(0,200,50),"%",sep=""),
       cex.axis=0.75,line=-0.5)
  mtext(3,line=1.4,text="Gene Set Enrichment",cex=1.2)
}


############################
###RARE LOF CARRIER ANALYSIS
############################
#Gather fraction of individuals from a given population with a rare pLoF in a given gene list
countRareLoFCarrierRate <- function(pop="ALL",genelist,dat,maxAF=0.001,mode="ALL"){
  #Get index of all variants with pLoF of at least one gene in list
  lof.in.genelist <- which(unlist(lapply(strsplit(dat$PROTEIN_CODING__LOF,split=","),function(genes){any(genes %in% genelist)})))
  subdat <- dat[lof.in.genelist,]
  #Set filtering indexes
  if(pop=="ALL"){
    prefix <- ""
  }else{
    prefix <- paste(pop,"_",sep="")
  }
  AF.idx <- which(colnames(subdat)==paste(prefix,"AF",sep=""))
  genos.idx <- which(colnames(subdat)==paste(prefix,"N_BI_GENOS",sep=""))
  het.idx <- which(colnames(subdat)==paste(prefix,"N_HET",sep=""))
  homalt.idx <- which(colnames(subdat)==paste(prefix,"N_HOMALT",sep=""))
  #Count non-ref individuals by SVTYPE
  svtypes.forCarrierAnalysis <- c("DEL","DUP","INS","INV","CPX")
  fracs <- sapply(svtypes.forCarrierAnalysis,function(svtype){
    hits <- subdat[which(subdat[,AF.idx]<maxAF & subdat$SVTYPE==svtype),]
    if(mode=="ALL"){
      sum(hits[,het.idx]+hits[,homalt.idx])/max(hits[,genos.idx],na.rm=T)
    }else if(mode=="HET"){
      sum(hits[,het.idx])/max(hits[,genos.idx],na.rm=T)
    }else if(mode=="HOM"){
      sum(hits[,homalt.idx])/max(hits[,genos.idx],na.rm=T)
    }else{
      stop("mode must be one of HOM, HET, ALL")
    }
  })
  return(fracs)
}
#Plot rare pLoF carrier rates for a list of gene lists
plot.rareLoFCarrierByList <- function(pop="ALL",genelists.list,dat,genelist.labels=NULL,barlabels=T,small=F,xmax=NULL,xlabel=NULL){
  res <- do.call("rbind", lapply(genelists.list,countRareLoFCarrierRate,pop=pop,dat=dat))
  #Prep plot area
  plot.colors <- c(svtypes$color[which(svtypes$svtype=="DEL")],
                   svtypes$color[which(svtypes$svtype=="DUP")],
                   svtypes$color[which(svtypes$svtype=="INS")],
                   svtypes$color[which(svtypes$svtype=="INV")],
                   svtypes$color[which(svtypes$svtype=="CPX")])
  par(bty="n")
  if(is.null(xmax)){
    xmax <- max(apply(res,1,sum))
  }
  plot(x=c(0,1.1*xmax),y=c(0,-nrow(res)),type="n",
       xaxt="n",xlab="",xaxs="i",yaxt="n",ylab="",yaxs="i")
  #Plot bars
  sapply(1:nrow(res),function(i){
    rect(xleft=c(0,cumsum(res[i,])[-ncol(res)]),
         xright=cumsum(res[i,]),
         ybottom=-i+0.2,ytop=-i+0.8,
         bty="n",border=NA,col=plot.colors)
    rect(xleft=0,xright=sum(res[i,]),
         ybottom=-i+0.2,ytop=-i+0.8,
         col=NA)
    if(barlabels==T){
      if(small==T){
        text(x=sum(res[i,]),y=-i+0.5,cex=0.6,
             pos=4,labels=paste(round(100*sum(res[i,]),1),"%",sep=""))
      }else{
        text(x=sum(res[i,]),y=-i+0.5,cex=0.8,
             pos=4,labels=paste(round(100*sum(res[i,]),1),"%",sep=""))
      }
    }
    if(!is.null(genelist.labels)){
      if(small==T){
        axis(2,at=-i+0.5,tick=F,line=-0.8,
             labels=genelist.labels[i],cex.axis=0.7,las=2)
      }else{
        axis(2,at=-i+0.5,tick=F,line=-0.8,
             labels=genelist.labels[i],cex.axis=0.85,las=2)
      }
    }
  })
  #Clean up
  if(small!=T){
    axis(1,at=axTicks(1),labels=NA)
    axis(1,at=axTicks(1),tick=F,line=-0.5,cex.axis=0.7,
         labels=paste(round(100*axTicks(1),1),"%",sep=""))
    if(is.null(xlabel)){
      xlabel <- "Very Rare (AF<0.1%) SV pLoF Carrier Rate"
    }
    mtext(1,line=1.5,text=xlabel)
  }
}

#Plot rare pLoF carrier rates for a list of gene lists, split by population
plot.rareLoFCarrierByListAndPop <- function(pops,dat,genelists.list,barlabels=T,xmax=NULL,xlabel=NULL){
  res <- lapply(genelists.list,function(g){
    sapply(pops,function(pop){
      countRareLoFCarrierRate(pop=pop,genelist=g,dat=dat)
    })
  })
  #Prep plot area
  plot.colors <- c(svtypes$color[which(svtypes$svtype=="DEL")],
                   svtypes$color[which(svtypes$svtype=="DUP")],
                   svtypes$color[which(svtypes$svtype=="INS")],
                   svtypes$color[which(svtypes$svtype=="INV")],
                   svtypes$color[which(svtypes$svtype=="CPX")])
  if(is.null(xmax)){
    xmax <- max(unlist(lapply(res,function(df){apply(df,2,sum)})))
  }
  par(bty="n",mfrow=c(length(res),1),mar=c(1.5,2,1.5,1))
  sapply(1:length(res),function(i){
    plot(x=c(0,1.125*xmax),y=c(0,-ncol(res[[i]])),type="n",
         xaxt="n",xlab="",xaxs="i",yaxt="n",ylab="",yaxs="i")
    #Plot bars
    plot.df <- res[[i]]
    sapply(1:ncol(plot.df),function(k){
      rect(xleft=c(0,cumsum(plot.df[,k])[-nrow(plot.df)]),
           xright=cumsum(plot.df[,k]),
           ybottom=-k+0.2,ytop=-k+0.8,
           bty="n",border=NA,col=plot.colors)
      rect(xleft=0,xright=sum(plot.df[,k]),
           ybottom=-k+0.2,ytop=-k+0.8,
           col=NA)
      if(barlabels==T){
          text(x=sum(plot.df[,k]),y=-k+0.5,cex=0.9,
               pos=4,labels=paste(round(100*sum(plot.df[,k]),1),"%",sep=""))
      }
      axis(2,at=-k+0.5,tick=F,labels=pops[k],cex.axis=0.9,las=2,line=-0.8)
    })
    axis(1,at=axTicks(1),labels=NA,tck=-0.04)
    axis(1,at=axTicks(1),tick=F,line=-0.75,cex.axis=0.8,
         labels=paste(round(100*axTicks(1),1),"%",sep=""))
  })
}

#Plot vertical barplot for rare pLoF carrier rates for a list of gene lists
plot.rareLoFCarrierByList.vertical <- function(pop="ALL",genelists.list,dat,modes,genelist.labels=NULL,
                                                 ymax=NULL,y.break=NULL,ylabel=NULL){
  res <- do.call("rbind", lapply(1:length(genelists.list),function(i){
    countRareLoFCarrierRate(genelist=genelists.list[[i]],pop=pop,dat=dat,mode=modes[i])
  }))
  #Prep plot area
  plot.colors <- c(svtypes$color[which(svtypes$svtype=="DEL")],
                   svtypes$color[which(svtypes$svtype=="DUP")],
                   svtypes$color[which(svtypes$svtype=="INS")],
                   svtypes$color[which(svtypes$svtype=="INV")],
                   svtypes$color[which(svtypes$svtype=="CPX")])
  par(bty="n")
  if(is.null(ymax)){
    ymax <- max(apply(res,1,sum))
  }
  if(!is.null(y.break)){
    ymax <- ymax-(y.break[2]-y.break[1])
  }
  par(mar=c(5.75,3.75,1,0.5),bty="n",xpd=F)
  plot(x=c(0,nrow(res))+1,y=c(0,ymax),type="n",
       xaxt="n",xlab="",xaxs="i",yaxt="n",ylab="",yaxs="i")
  par(xpd=T)
  rect(xleft=1,xright=par("usr")[2],
       ybottom=par("usr")[3]-(0.06*(par("usr")[4]-par("usr")[3])),
       ytop=par("usr")[3]-(0.01*(par("usr")[4]-par("usr")[3])),
       bty="n",border=NA,col="gray90")
  par(xpd=F)
  #Plot bars
  sapply(1:nrow(res),function(i){
    rect(ybottom=c(0,cumsum(res[i,])[-ncol(res)]),
         ytop=cumsum(res[i,]),
         xleft=i+0.2,xright=i+0.8,
         bty="n",border=NA,col=plot.colors)
    rect(ybottom=0,ytop=sum(res[i,]),
         xleft=i+0.2,xright=i+0.8,
         col=NA)
    text(y=sum(res[i,])-(0.025*(par("usr")[4]-par("usr")[3])),
         x=i+0.5,cex=0.75,
         pos=3,labels=paste(round(100*sum(res[i,]),1),"%",sep=""))
    if(!is.null(genelist.labels)){
      par(xpd=T)
      text(x=i+0.5,y=par("usr")[3]-(0.035*(par("usr")[4]-par("usr")[3])),
           labels=prettyNum(length(genelists.list[[i]]),big.mark=","),
           cex=0.8)
      text(x=i+1,y=par("usr")[3]-(0.09*(par("usr")[4]-par("usr")[3])),
           labels=genelist.labels[i],srt=40,pos=2,cex=0.9)
      par(xpd=F)
    }
  })
  #Clean up bottom panel
  axis(2,at=c(par("usr")[3],ymax),labels=NA,tck=0)
  if(!is.null(y.break)){
    axis(2,at=axTicks(2)[which(axTicks(2)<y.break[1])],labels=NA,tck=-0.025)
    axis(2,at=axTicks(2)[which(axTicks(2)<y.break[1])],tick=F,line=-0.6,cex.axis=0.9,las=2,
         labels=paste(round(100*axTicks(2)[which(axTicks(2)<y.break[1])],1),"%",sep=""))
  }
  mtext(2,line=2.2,text=ylabel,cex=1.4)
  #Add top panel, if necessary
  if(!is.null(y.break)){
    #Clear all data above future axis break
    rect(xleft=par("usr")[1]+0.2,xright=par("usr")[2]-0.2,
         ybottom=y.break[1],
         ytop=par("usr")[4],
         col="white",border="white",lwd=4)
    #Add stacked rectangles as required
    sapply(which(apply(res,1,sum)>y.break[1]),function(i){
      revised.vals <- cumsum(res[i,])-(y.break[2]-y.break[1])
      rect(xleft=i+0.2,xright=i+0.8,
           ybottom=c(y.break[1],revised.vals[-length(revised.vals)]),
           ytop=revised.vals,
           bty="n",border=NA,col=plot.colors)
      rect(xleft=i+0.2,xright=i+0.8,
           ybottom=0,ytop=max(revised.vals),
           col=NA)
      #Add bar labels
      text(y=max(revised.vals)-(0.025*(par("usr")[4]-par("usr")[3])),
           x=i+0.5,cex=0.75,
           pos=3,labels=paste(round(100*sum(res[i,]),1),"%",sep=""))
    })
    #Add top axis
    axis(2,at=axTicks(2)[which(axTicks(2)>y.break[1])],labels=NA,tck=-0.025)
    axis(2,at=axTicks(2)[which(axTicks(2)>y.break[1])],tick=F,line=-0.6,cex.axis=0.9,las=2,
         labels=paste(round(100*(axTicks(2)[which(axTicks(2)>y.break[1])]+(y.break[2]-y.break[1])),1),"%",sep=""))
    #Add axis break
    rect(xleft=par("usr")[1],xright=par("usr")[2],
         ybottom=y.break[1]-(0.01*(par("usr")[4]-par("usr")[3])),
         ytop=y.break[1]+(0.01*(par("usr")[4]-par("usr")[3])),
         col="white",bty="n",border=NA)
    axis.break(2,breakpos=y.break[1])
    # abline(h=c(y.break[1]-(0.01*(par("usr")[4]-par("usr")[3])),
    #            y.break[1]+(0.01*(par("usr")[4]-par("usr")[3]))),
    #        lty=2,lwd=0.75)
  }
}



################
###RSCRIPT BLOCK
################
require(optparse,quietly=T)
require(plotrix,quietly=T)
require(boot,quietly=T)
require(MASS,quietly=T)
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
args <- parse_args(OptionParser(usage="%prog VCF2BED SNV_DATA GENE_METADATA GENELIST_DIRECTORY OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 5){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
vcf2bed.in <- args$args[1]
SNVdata.in <- args$args[2]
gene_metadata.in <- args$args[3]
genelist.dir <- args$args[4]
OUTDIR <- args$args[5]
prefix <- opts$prefix
svtypes.file <- opts$svtypes
pops.file <- opts$populations

# #Dev parameters (local)
# vcf2bed.in <- "~/scratch/gnomAD_v2_SV_MASTER.vcf2bed.bed.gz"
# SNVdata.in <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/gnomAD_v2.1_canonical_constraint.condensed.txt.gz"
# gene_metadata.in <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/Gencode.v19.autosomal_canonical_gene_metadata.txt.gz"
# genelist.dir <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/genelists/"
# svtypes.file <- "~/Desktop/Collins/Talkowski/code/sv-pipeline/ref/vcf_qc_refs/SV_colors.txt"
# OUTDIR <- "~/scratch/gnomAD_v2_test_plots/"
# prefix <- "gnomAD_v2_SV_MASTER"
# pops.file <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/gnomAD_population_colors.txt"

###Create output directory, if necessary
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

###Process input data
dat.all <- read.vcf2bed(vcf2bed.in)
dat <- dat.all[which(dat.all$FILTER=="PASS" | dat.all$FILTER=="MULTIALLELIC"),]
snv.data <- read.snvdata(SNVdata.in,gene_metadata.in)
genes <- sort(unique(as.character(snv.data$gene)))

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

###Get fraction of SV class with functional consequences by frequency bin
all.fracs <- count.fracLOFSV(dat[which(!(dat$chrom %in% c("X","Y"))),])
common.fracs <- count.fracLOFSV(dat[which(dat$AF>=0.01 & !(dat$chrom %in% c("X","Y"))),])
rare.fracs <- count.fracLOFSV(dat[which(dat$AF<0.01 & dat$AC>1 & !(dat$chrom %in% c("X","Y"))),])
singleton.fracs <- count.fracLOFSV(dat[which(dat$AC==1 & !(dat$chrom %in% c("X","Y"))),])
merged.fracs <- data.frame("all"=all.fracs,
                           "common"=common.fracs,
                           "rare"=rare.fracs,
                           "singleton"=singleton.fracs)
pdf(paste(OUTDIR,"/",prefix,".site_fraction_by_functional_effct_barplots.pdf",sep=""),
    height=2.8,width=4.75)
par(mfrow=c(2,1))
barplot.fracByEffect(frac.mat=merged.fracs[1:5,],
                     base.colors=c("gray30",
                                   svtypes$color[which(svtypes$svtype=="DEL")],
                                   svtypes$color[which(svtypes$svtype=="INS")],
                                   svtypes$color[which(svtypes$svtype=="INV")],
                                   svtypes$color[which(svtypes$svtype=="CPX")]),
                     category.labels=c("pLoF\n(All SV)","pLoF\n(DEL)","pLoF\n(INS)",
                                       "pLoF\n(INV)","pLoF\n(CPX)"))
barplot.fracByEffect(frac.mat=merged.fracs[-c(1:5),],
                     base.colors=c(svtypes$color[which(svtypes$svtype=="MCNV")],
                                   svtypes$color[which(svtypes$svtype=="DUP")],
                                   svtypes$color[which(svtypes$svtype=="INV")],
                                   "gray30","gray30"),
                     category.labels=c("IED\n(DUP)","CG\n(DUP+CPX)",
                                       "Inverted\nGene\n(INV+CPX)","Intronic\n(All SV)","Promoter\n(All SV)"))
dev.off()

###Overlap SV & SNV data into master gene table
gene.data <- getSVdat.all(dat=dat,snv.data=snv.data)
write.table(gene.data,paste(OUTDIR,"/",prefix,".gene_SV_data.txt",sep=""),
            col.names=T,row.names=F,sep="\t",quote=F)
constrained.genes <- gene.data$gene[which(gene.data$pli>0.9)]
unconstrained.genes <- gene.data$gene[which(gene.data$pli<0.1)]

###Get proportion of singletons per functional class
fracSingletons <- count.fracSingletons(dat)
write.table(fracSingletons,paste(OUTDIR,"/",prefix,".singleton_fraction_by_coding_annotation.txt",sep=""),
            col.names=T,row.names=F,sep="\t",quote=F)
#Dotplot of singleton fraction
pdf(paste(OUTDIR,"/",prefix,".fraction_singletons_by_coding_annotation.pdf",sep=""),
    height=4.5,width=2.6)
dotplot.fracSingletons(fracSingletons)
dev.off()
pdf(paste(OUTDIR,"/",prefix,".fraction_singletons_by_coding_annotation.horizontal.pdf",sep=""),
    height=2.25,width=4.3)
dotplot.fracSingletons(fracSingletons,horiz=T,add.pt.labs=F,add.N=F)
dev.off()

###Get proportion of singletons by SVTYPE for those with coding effects vs intergenic
fracSingletons.coding_vs_non <- count.fracSingletons.coding_vs_noncoding(dat)
pdf(paste(OUTDIR,"/",prefix,".fraction_singletons_by_SVTYPE_matched_coding_vs_non.pdf",sep=""),
    height=2.25,width=4.9)
dotplot.fracSingletons.coding_vs_non(fracSingletons.coding_vs_non,add.N=F)
dev.off()

###Barplot of pct of genes with each functional class
pdf(paste(OUTDIR,"/",prefix,".pct_genes_with_functional_SV.pdf",sep=""),
    height=1.75,width=4.25)
barplot.genesFracByEffect(gene.data)
dev.off()

###Histogram quadrant of SV per gene by functional class
pdf(paste(OUTDIR,"/",prefix,".functional_SV_per_gene_hist.pdf",sep=""),
    height=2.5,width=3.25)
par(mfrow=c(2,2))
minihist.byFunc(gene.data=gene.data,func="lof.any",title="pLoF",
                base.col=svtypes$color[which(svtypes$svtype =="DEL")])
minihist.byFunc(gene.data=gene.data,func="plof",title="IED",
                base.col=svtypes$color[which(svtypes$svtype =="MCNV")])
minihist.byFunc(gene.data=gene.data,func="cg",title="CG",
                base.col=svtypes$color[which(svtypes$svtype =="DUP")])
minihist.byFunc(gene.data=gene.data,func="inv",title="Inverted",
                base.col=svtypes$color[which(svtypes$svtype =="INV")])
dev.off()
pdf(paste(OUTDIR,"/",prefix,".functional_SV_per_gene_hist.long.pdf",sep=""),
    height=0.5*2.5,width=6.25)
par(mfrow=c(1,4))
minihist.byFunc(gene.data=gene.data,func="lof.any",title="pLoF",
                base.col=svtypes$color[which(svtypes$svtype =="DEL")])
minihist.byFunc(gene.data=gene.data,func="plof",title="IED",
                base.col=svtypes$color[which(svtypes$svtype =="MCNV")])
minihist.byFunc(gene.data=gene.data,func="cg",title="CG",
                base.col=svtypes$color[which(svtypes$svtype =="DUP")])
minihist.byFunc(gene.data=gene.data,func="inv",title="Inverted",
                base.col=svtypes$color[which(svtypes$svtype =="INV")])
dev.off()
pdf(paste(OUTDIR,"/",prefix,".functional_SV_per_gene_hist.tall_noInv.pdf",sep=""),
    height=4,width=1.5)
par(mfrow=c(3,1))
minihist.byFunc(gene.data=gene.data,func="lof.any",title="pLoF",
                base.col=svtypes$color[which(svtypes$svtype =="DEL")])
minihist.byFunc(gene.data=gene.data,func="plof",title="IED",
                base.col=svtypes$color[which(svtypes$svtype =="MCNV")])
minihist.byFunc(gene.data=gene.data,func="cg",title="CG",
                base.col=svtypes$color[which(svtypes$svtype =="DUP")])
dev.off()

###Heatmaps of count of SV per gene between two functional classes
pdf(paste(OUTDIR,"/",prefix,".functional_SV_per_gene_heatmaps.pdf",sep=""),
    height=2.5,width=7.5)
par(mfrow=c(1,3))
plot.perGeneHeatmaps(gene.data,x.category="all.lof.any",y.category="all.plof",
                     x.label="pLoF per Gene",y.label="IED per Gene")
plot.perGeneHeatmaps(gene.data,x.category="all.lof.any",y.category="all.cg",
                     x.label="pLoF per Gene",y.label="CG per Gene")
plot.perGeneHeatmaps(gene.data,x.category="all.plof",y.category="all.cg",
                     x.label="CG per Gene",y.label="IED per Gene")
dev.off()
pdf(paste(OUTDIR,"/",prefix,".functional_SV_per_gene_heatmaps.key.pdf",sep=""),
    height=2,width=0.25)
plot.perGeneHeatmapsKey()
dev.off()

###Circle plots of count of SV per gene between two functional classes vs permuted expectation
pdf(paste(OUTDIR,"/",prefix,".functional_SV_per_gene_dotplots_vs_permutation.pdf",sep=""),
    height=2.5,width=7.5)
par(mfrow=c(1,3))
plot.perGeneDotplots(gene.data,x.category="all.lof.any",y.category="all.plof",
                     x.label="pLoF per Gene",y.label="IED per Gene")
plot.perGeneDotplots(gene.data,x.category="all.lof.any",y.category="all.cg",
                     x.label="pLoF per Gene",y.label="CG per Gene")
plot.perGeneDotplots(gene.data,x.category="all.plof",y.category="all.cg",
                     x.label="IED per Gene",y.label="CG per Gene")
dev.off()
pdf(paste(OUTDIR,"/",prefix,".functional_SV_per_gene_dotplots_vs_permutation.key.pdf",sep=""),
    height=2.5,width=0.5)
plot.perGeneDotplotKey()
dev.off()

###Constraint comparisons
pdf(paste(OUTDIR,"/",prefix,".snv_constraint_comparison_by_SV_class.SV_lof_and_dup.pdf",sep=""),
    height=2*1.75,width=5)
par(mfrow=c(2,2))
constraintModelWrapper(gene.data=gene.data,
                       y="rare.lof.any",
                       color=svtypes$color[which(svtypes$svtype=="DEL")],
                       title="pLoF SV",
                       ymax=1.5)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.cg",
                       color=svtypes$color[which(svtypes$svtype=="DUP")],
                       title="Gene Duplications",
                       ymax=1.5)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.lof.any",
                       snv.class="mis",
                       color=svtypes$color[which(svtypes$svtype=="DEL")],
                       title="pLoF SV",
                       ymax=1.5)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.cg",
                       snv.class="mis",
                       color=svtypes$color[which(svtypes$svtype=="DUP")],
                       title="Gene Duplications",
                       ymax=1.5)
dev.off()
pdf(paste(OUTDIR,"/",prefix,".snv_constraint_comparison_by_SV_class.SV_duplof_and_inv.pdf",sep=""),
    height=2*1.75,width=5)
par(mfrow=c(2,2))
constraintModelWrapper(gene.data=gene.data,
                       y="rare.plof",
                       color=svtypes$color[which(svtypes$svtype=="MCNV")],
                       title="Int. Exon Dup.",
                       ymax=2)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.inv",
                       color=svtypes$color[which(svtypes$svtype=="INV")],
                       title="Gene Inversions",
                       ymax=2)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.plof",
                       snv.class="mis",
                       color=svtypes$color[which(svtypes$svtype=="MCNV")],
                       title="Int. Exon Dup",
                       ymax=2)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.inv",
                       snv.class="mis",
                       color=svtypes$color[which(svtypes$svtype=="INV")],
                       title="Gene Inversions",
                       ymax=2)
dev.off()
pdf(paste(OUTDIR,"/",prefix,".snv_constraint_comparison_by_SV_class.SV_lof_by_subsets.pdf",sep=""),
    height=2*1.5,width=6)
par(mfrow=c(2,3))
constraintModelWrapper(gene.data=gene.data,
                       y="rare.lof.any",
                       color=svtypes$color[which(svtypes$svtype=="DEL")],
                       title="All pLoF SV",
                       ymax=1.5,cex.labels=0.7)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.lof.del",
                       color=svtypes$color[which(svtypes$svtype=="DEL")],
                       title="pLoF Deletions",
                       ymax=1.5,cex.labels=0.7)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.lof.other",
                       color=svtypes$color[which(svtypes$svtype=="DEL")],
                       title="Other pLoF SV",
                       ymax=1.5,cex.labels=0.7)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.lof.any",
                       color=svtypes$color[which(svtypes$svtype=="DEL")],
                       title="All pLoF SV",
                       snv.class="mis",
                       ymax=1.5,cex.labels=0.7)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.lof.del",
                       color=svtypes$color[which(svtypes$svtype=="DEL")],
                       title="pLoF Deletions",
                       snv.class="mis",
                       ymax=1.5,cex.labels=0.7)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.lof.other",
                       color=svtypes$color[which(svtypes$svtype=="DEL")],
                       title="Other pLoF SV",
                       snv.class="mis",
                       ymax=1.5,cex.labels=0.7)
dev.off()
pdf(paste(OUTDIR,"/",prefix,".snv_constraint_comparison_by_SV_class.all_panels.pdf",sep=""),
    height=2.25,width=6)
par(mfrow=c(2,4))
parmar <- c(1.5,3.25,1.5,0.75)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.lof.any",
                       color=svtypes$color[which(svtypes$svtype=="DEL")],
                       title="",ymax=1.5,ax.labels=F,cex.labels=0.65,
                       parmar=parmar)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.plof",
                       color=svtypes$color[which(svtypes$svtype=="MCNV")],
                       title="",ymax=1.5,ax.labels=F,cex.labels=0.65,
                       parmar=parmar)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.cg",
                       color=svtypes$color[which(svtypes$svtype=="DUP")],
                       title="",ymax=1.5,ax.labels=F,cex.labels=0.65,
                       parmar=parmar)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.inv",
                       color=svtypes$color[which(svtypes$svtype=="INV")],
                       title="",ymax=1.5,ax.labels=F,cex.labels=0.65,
                       parmar=parmar)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.lof.any",
                       color=svtypes$color[which(svtypes$svtype=="DEL")],
                       title="",ymax=1.5,ax.labels=F,
                       snv.class="mis",cex.labels=0.65,
                       parmar=parmar)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.plof",
                       color=svtypes$color[which(svtypes$svtype=="MCNV")],
                       title="",ymax=1.5,ax.labels=F,
                       snv.class="mis",cex.labels=0.65,
                       parmar=parmar)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.cg",
                       color=svtypes$color[which(svtypes$svtype=="DUP")],
                       title="",ymax=1.5,ax.labels=F,
                       snv.class="mis",cex.labels=0.65,
                       parmar=parmar)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.inv",
                       color=svtypes$color[which(svtypes$svtype=="INV")],
                       title="",ymax=1.5,ax.labels=F,
                       snv.class="mis",cex.labels=0.65,
                       parmar=parmar)
dev.off()
pdf(paste(OUTDIR,"/",prefix,".snv_constraint_comparison_by_SV_class.all_panels.lof_constraint_only.pdf",sep=""),
    height=0.45*2.25,width=5.75)
par(mfrow=c(1,4))
parmar <- c(1.25,2,0.5,0.75)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.lof.any",
                       color=svtypes$color[which(svtypes$svtype=="DEL")],
                       title="",ymax=1.5,ax.labels=F,cex.labels=0.85,
                       tck=-0.02,yline=-0.85,parmar=parmar)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.plof",
                       color=svtypes$color[which(svtypes$svtype=="MCNV")],
                       title="",ymax=1.5,ax.labels=F,cex.labels=0.85,
                       tck=-0.02,yline=-0.85,parmar=parmar)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.cg",
                       color=svtypes$color[which(svtypes$svtype=="DUP")],
                       title="",ymax=1.5,ax.labels=F,cex.labels=0.85,
                       tck=-0.02,yline=-0.85,parmar=parmar)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.inv",
                       color=svtypes$color[which(svtypes$svtype=="INV")],
                       title="",ymax=1.5,ax.labels=F,cex.labels=0.85,
                       tck=-0.02,yline=-0.85,parmar=parmar)
dev.off()
pdf(paste(OUTDIR,"/",prefix,".snv_constraint_comparison_by_SV_class.all_panels.missense_constraint_only.pdf",sep=""),
    height=0.45*2.25,width=5.75)
par(mfrow=c(1,4))
parmar <- c(1.25,2,0.5,0.75)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.lof.any",
                       snv.class="mis",
                       color=svtypes$color[which(svtypes$svtype=="DEL")],
                       title="",ymax=1.5,ax.labels=F,cex.labels=0.85,
                       tck=-0.02,yline=-0.85,parmar=parmar)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.plof",
                       snv.class="mis",
                       color=svtypes$color[which(svtypes$svtype=="MCNV")],
                       title="",ymax=1.5,ax.labels=F,cex.labels=0.85,
                       tck=-0.02,yline=-0.85,parmar=parmar)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.cg",
                       snv.class="mis",
                       color=svtypes$color[which(svtypes$svtype=="DUP")],
                       title="",ymax=1.5,ax.labels=F,cex.labels=0.85,
                       tck=-0.02,yline=-0.85,parmar=parmar)
constraintModelWrapper(gene.data=gene.data,
                       y="rare.inv",
                       snv.class="mis",
                       color=svtypes$color[which(svtypes$svtype=="INV")],
                       title="",ymax=1.5,ax.labels=F,cex.labels=0.85,
                       tck=-0.02,yline=-0.85,parmar=parmar)
dev.off()

#Build table of gene symbol, observed rare pLoF deletion, and expected rare pLoF deletion (for Konrad/pLoF paper)
deletion.constraint.data <- constraintModelWrapper(gene.data=gene.data,
                                                   y="rare.lof.del",
                                                   color=svtypes$color[which(svtypes$svtype=="DEL")],
                                                   title="pLoF Deletions",
                                                   ymax=1.5,return=T)
colnames(deletion.constraint.data) <- c("#gene","observed_rare_biallelic_LoF_deletions","expected_rare_biallelic_LoF_deletions")
write.table(deletion.constraint.data,
            paste(OUTDIR,"/",prefix,".rare_biallelic_LoF_deletions_per_gene.obs_exp.txt",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")

#Trio of mini histograms with # of pLoF SV per gene by constrained classification
pdf(paste(OUTDIR,"/",prefix,".rare_LoF_SV_per_gene_by_constraint_classification.pdf",sep=""),
    height=1.8,width=5.5)
par(mfrow=c(1,3))
minihist.byConstraint(gene.data=gene.data,min.pli=0.9,max.pli=1,
                      title="Constrained\n(pLI > 0.9)",x.max=5)
minihist.byConstraint(gene.data=gene.data,min.pli=0.1,max.pli=0.9,
                      title="Indeterminate\n(0.1 < pLI < 0.9)",x.max=5)
text(x=mean(c(3,floor(par("usr")[2]))),y=0.2*par("usr")[4],
     labels=paste(prettyNum(length(which(gene.data$pli>0.1 & gene.data$pli<0.9 & gene.data$rare.lof.any>2)),big.mark=","),
                  " genes with\nambiguous constraint\nand > 2 rare pLoF SV",sep=""),cex=0.8,pos=3)
minihist.byConstraint(gene.data=gene.data,min.pli=0,max.pli=0.1,
                      title="Not Constrained\n(pLI < 0.1)",x.max=5)
dev.off()

#Human KO analysis part 1: SV data only
KO.data <- getKOdata(dat=dat,gene.data=gene.data)
pdf(paste(OUTDIR,"/",prefix,".homozygous_LoF_breakdown.pdf",sep=""),
    height=1,width=5.5)
plotKO(KO.data=KO.data)
dev.off()

#Human KO analysis part 2: comparison to other KO studies
narasimhan.kos <- importGenelist(genelist.dir,
                                 "human_KOs.Narasimhan_2016.genes.list",
                                 gene.data)
saleheen.kos <- importGenelist(genelist.dir,
                               "human_KOs.Saleheen_2017.genes.list",
                               gene.data)
# SV.kos <- sort(unique(as.character(unlist(strsplit(dat$PROTEIN_CODING__LOF[which(dat$N_HOMALT>0 & !is.na(dat$PROTEIN_CODING__LOF) & !(dat$chrom %in% c("X","Y")))],split=",")))))
SV.kos <- sort(unique(as.character(unlist(strsplit(dat$PROTEIN_CODING__LOF[which(dat$N_HOMALT>0 & dat$POPMAX_AF<0.01 & !is.na(dat$PROTEIN_CODING__LOF) & !(dat$chrom %in% c("X","Y")))],split=",")))))
SV.kos <- SV.kos[which(SV.kos %in% gene.data$gene)]
ko.intersections <- intersect.KOstudies(SV.kos,narasimhan.kos,saleheen.kos)
pdf(paste(OUTDIR,"/",prefix,".human_knockout_cross_study_comparison.pdf",sep=""),
    height=2.4,width=4)
plot.KOstudyComparison(ko.intersections)
dev.off()

#Human KO analysis part 3: comparison vs gnomAD SNV data
base.colors <- c("gray70","gray30","white",NA,"#FFF2E1","#DAC691",
                 svtypes$color[which(svtypes$svtype=="DEL")])
pdf(paste(OUTDIR,"/",prefix,".human_knockout_ptv_oe_boxplots.pdf",sep=""),
    height=1.25,width=2)
boxplots.oe(genesets.list=list(gene.data$gene,
                               gene.data$gene[which(gene.data$pli>0.9)],
                               gene.data$gene[which(gene.data$pli<0.1)],
                               NA,
                               saleheen.kos,
                               narasimhan.kos,
                               SV.kos),
            gene.data=gene.data,snv.category="ptv",base.colors=base.colors,ymax=2)
dev.off()
pdf(paste(OUTDIR,"/",prefix,".human_knockout_mis_oe_boxplots.pdf",sep=""),
    height=1.25,width=2)
boxplots.oe(genesets.list=list(gene.data$gene,
                               gene.data$gene[which(gene.data$pli>0.9)],
                               gene.data$gene[which(gene.data$pli<0.1)],
                               NA,
                               saleheen.kos,
                               narasimhan.kos,
                               SV.kos),
            gene.data=gene.data,snv.category="mis",base.colors=base.colors,ymax=1.5)
dev.off()
pdf(paste(OUTDIR,"/",prefix,".human_knockout_size_boxplots.pdf",sep=""),
    height=1.25,width=2)
boxplots.size(genesets.list=list(gene.data$gene,
                                 gene.data$gene[which(gene.data$pli>0.9)],
                                 gene.data$gene[which(gene.data$pli<0.1)],
                                 NA,
                                 saleheen.kos,
                                 narasimhan.kos,
                                 SV.kos),
              gene.data=gene.data,base.colors=base.colors)
dev.off()

#Human KO analysis part 4: enrichment vs external datasets
AD.genes <- importGenelist(genelist.dir,
                           "Autosomal_dominant_disease.genes.list",
                           gene.data)
AR.genes <- importGenelist(genelist.dir,
                           "Autosomal_recessive_disease.genes.list",
                           gene.data)
FDA.genes <- importGenelist(genelist.dir,
                            "FDA_drug_targets.genes.list",
                            gene.data)
essential.genes <- importGenelist(genelist.dir,
                                  "Culture_essential.genes.list",
                                  gene.data)
olfactory.genes <- importGenelist(genelist.dir,
                                  "Olfactory_receptors.genes.list",
                                  gene.data)
KO.enrichment.data <- getGeneListFractionMultiple(comparison.sets.list=list(essential.genes,olfactory.genes,AD.genes,AR.genes,FDA.genes),
                                                  comparison.set.names=c("Cell Culture\nEssential Genes",
                                                                         "Olfactory\nReceptors",
                                                                         "Autosomal Dominant\nDisease Genes",
                                                                         "Autosomal Recessive\nDisease Genes",
                                                                         "FDA-Approved\nDrug Targets"),
                                                  genesets.list=list(gene.data$gene,
                                                                     gene.data$gene[which(gene.data$pli>0.9)],
                                                                     gene.data$gene[which(gene.data$pli<0.1)],
                                                                     saleheen.kos,
                                                                     narasimhan.kos,
                                                                     SV.kos),
                                                  geneset.names=c("All genes","pLI>0.9","pLI<0.1",
                                                                  "Saleheen","Narasimhan","gnomAD-SV"))
pdf(paste(OUTDIR,"/",prefix,".human_knockout_geneset_enrichments.pdf",sep=""),
    height=2.5,width=3.5)
dotplots.KOenrichments(KO.enrichment.data=KO.enrichment.data,
                       colors=c("gray70","gray30","white","#FFF2E1","#DAC691",
                                svtypes$color[which(svtypes$svtype=="DEL")]))
dev.off()

#Carrier rates for rare pLoF SV in various medically relevant gene lists
ACMG.genes <- importGenelist(genelist.dir,
                             "ACMG_59.genes.list",
                             gene.data)
DDD.genes <- importGenelist(genelist.dir,
                            "DDD_2017.genes.list",
                            gene.data)
clingenHC.genes <- importGenelist(genelist.dir,
                                  "ClinGen_haploinsufficient_high_confidence.genes.list",
                                  gene.data)
clingenMC.genes <- importGenelist(genelist.dir,
                                  "ClinGen_haploinsufficient_medium_confidence.genes.list",
                                  gene.data)
clingenLC.genes <- importGenelist(genelist.dir,
                                  "ClinGen_haploinsufficient_low_confidence.genes.list",
                                  gene.data)
dominant_DD.genes <- importGenelist(genelist.dir,
                                    "DDG2P_confirmed_LoF_dominant_DD.genes.list",
                                    gene.data)
recessive_DD.genes <- importGenelist(genelist.dir,
                                    "DDG2P_confirmed_LoF_recessive_DD.genes.list",
                                    gene.data)
#Main panel - cohort-wide
pdf(paste(OUTDIR,"/",prefix,".rare_diagnostic_lof_carrierRates_ALL.pdf",sep=""),
    height=3.5,width=4.25)
par(mar=c(2.5,6.5,0.5,1),xpd=T)
plot.rareLoFCarrierByList(pop="ALL",dat=dat,
                          genelists.list=list(ACMG.genes,DDD.genes,clingenHC.genes,
                                              clingenMC.genes,clingenLC.genes),
                          genelist.labels=c("ACMG 59\nReportable",
                                            "DDD 76",
                                            "ClinGen\nHaploinsufficient\n(High Conf.)",
                                            "ClinGen\nHaploinsufficient\n(Med. Conf.)",
                                            "ClinGen\nHaploinsufficient\n(Low Conf.)"),
                          xlabel="Very Rare (AF<0.1%) pLoF SV Carrier Rate")
dev.off()
#Smaller panel - cohort-wide
pdf(paste(OUTDIR,"/",prefix,".rare_diagnostic_lof_carrierRates_ALL.smaller.pdf",sep=""),
    height=2.5,width=3.25)
par(mar=c(2.5,4.5,0.5,1.25),xpd=T)
plot.rareLoFCarrierByList(pop="ALL",dat=dat,
                          genelists.list=list(ACMG.genes,DDD.genes,clingenHC.genes,
                                              clingenMC.genes,clingenLC.genes),
                          genelist.labels=c("ACMG 59\nReportable",
                                            "DDD 76",
                                            "ClinGen HI\n(High-Conf.)",
                                            "ClinGen HI\n(Med.-Conf.)",
                                            "ClinGen HI\n(Low-Conf.)"),
                          xlabel="Very Rare pLoF SV Carrier Rate")
dev.off()

#Smaller panels - individual populations
sapply(c("AFR","EAS","EUR","AMR"),function(pop){
  pdf(paste(OUTDIR,"/",prefix,".rare_diagnostic_lof_carrierRates_",pop,".pdf",sep=""),
      height=0.75,width=2.5)
  par(mar=c(0.5,4,0.5,1.25),xpd=T)
  plot.rareLoFCarrierByList(pop=pop,dat=dat,
                            genelists.list=list(ACMG.genes,DDD.genes,clingenHC.genes,
                                                clingenMC.genes,clingenLC.genes),
                            genelist.labels=c("ACMG",
                                              "DDD",
                                              "ClinGen High",
                                              "ClinGen Med.",
                                              "ClinGen\ Low."),
                            small=T,xmax=0.02)
  dev.off()
})

#Blowout panels - individual populations
pdf(paste(OUTDIR,"/",prefix,".rare_diagnostic_lof_carrierRates.ACMG_DDD_clinGenLC.CrossPop.pdf",sep=""),
    height=2.5,width=1.9)
par(xpd=T)
plot.rareLoFCarrierByListAndPop(pops=c("AFR","EAS","EUR","AMR"),
                                dat=dat,
                                genelists.list=list(ACMG.genes,DDD.genes,clingenLC.genes))
par(xpd=F)
dev.off()

#Vertical panel
pdf(paste(OUTDIR,"/",prefix,".rare_diagnostic_lof_carrierRates.vertical.pdf",sep=""),
    height=3.75,width=2.5)
plot.rareLoFCarrierByList.vertical(pop="ALL",dat=dat,
                                   genelists.list=list(ACMG.genes,dominant_DD.genes,
                                                       clingenHC.genes,
                                                       recessive_DD.genes,recessive_DD.genes),
                                   modes=c(rep("ALL",3),"HET","HOM"),
                                   genelist.labels=c("ACMG Reportable",
                                                     "Dominant DD",
                                                     "ClinGen Haploinsuff.",
                                                     "Recessive DD",
                                                     "Recessive DD"),
                                   y.break=c(0.015,0.085),ymax=0.10,
                                   ylabel="Samples with pLoF SV")
dev.off()

#Vertical panel, no homozygous recessive DD
pdf(paste(OUTDIR,"/",prefix,".rare_diagnostic_lof_carrierRates.vertical.noHomDD.pdf",sep=""),
    height=3.75,width=2.5)
plot.rareLoFCarrierByList.vertical(pop="ALL",dat=dat,
                                   genelists.list=list(ACMG.genes,dominant_DD.genes,
                                                       clingenHC.genes,
                                                       recessive_DD.genes),
                                   modes=c(rep("ALL",3),"HET"),
                                   genelist.labels=c("ACMG Reportable",
                                                     "Dominant DD",
                                                     "ClinGen Haploinsuff.",
                                                     "Recessive DD"),
                                   y.break=c(0.015,0.085),ymax=0.10,
                                   ylabel="Samples with pLoF SV")
dev.off()
