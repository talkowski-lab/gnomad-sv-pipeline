#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Perform tertiary exploration of singleton rates for various supplementary figures


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

#Get indexes corresponding to SV with a functional consequence stratified by PTV O:E
get.func.indexes.byConstraint <- function(dat,func="LOF",sd_sr.max=1.1,bins=NULL){
  #Instantiate constraint bins (if null)
  if(is.null(bins)){
    bins <- data.frame(seq(0,90,10),seq(10,100,10))
  }
  #Subset to all variants meeting func criteria
  func.idx <- which(colnames(dat)==paste("PROTEIN_CODING__",func,sep=""))
  subdat <- dat[which(!is.na(dat[,func.idx]) 
                      & !(dat$chrom %in% c("X","Y"))
                      & dat$SD_SR_OVERLAP<sd_sr.max),]
  multiGene <- grep(",",subdat[,func.idx],fixed=T)
  if(length(multiGene)>0){
    subdat <- subdat[-multiGene,]
  }
  if(func %in% c("INV_SPAN","DUP_PARTIAL","INTRONIC","UTR","PROMOTER","NEAREST_TSS")){
    subdat <- subdat[which(is.na(subdat$PROTEIN_CODING__LOF)
                           & is.na(subdat$PROTEIN_CODING__DUP_LOF)
                           & is.na(subdat$PROTEIN_CODING__COPY_GAIN)),]
  }
  if(func=="NEAREST_TSS"){
    subdat <- subdat[which(subdat$PROTEIN_CODING__INTERGENIC=="True"
                           & is.na(subdat$PROTEIN_CODING__PROMOTER)
                           & is.na(subdat$PROTEIN_CODING__INTRONIC)
                           & is.na(subdat$PROTEIN_CODING__UTR)),]
  }
  genes.per.sv <- strsplit(subdat[,func.idx],split=",")
  min.oe.per.sv <- unlist(lapply(genes.per.sv,function(genes){
    gene.matches <- snv.data$ptv_oe_cent[which(snv.data$gene %in% genes)]
    if(length(gene.matches)==0){
      return(NA)
    }else{
      return(min(gene.matches))
    }
  }))
  subdat <- subdat[which(!is.na(min.oe.per.sv)),]
  min.oe.per.sv <- min.oe.per.sv[which(!is.na(min.oe.per.sv))]
  #Iterate over constraint deciles & return indexes
  sv.idxs <- apply(bins,1,function(bin){
    hits <- which(min.oe.per.sv>as.numeric(bin[1])
                  & min.oe.per.sv<=as.numeric(bin[2]))
    which(dat$name %in% subdat$name[hits])
  })
  #Add first row for all genes meeting functional criteria
  # and spacer row before getting to constraint
  sv.idxs <- c(list(which(dat$name %in% subdat$name)),
               list(NA),
               sv.idxs)
  return(sv.idxs)
}

#Get proportion of singletons and 95% CI for a set of variants
calc.fracSingletons.single <- function(dat,indexes,boot.n=100,conf=0.95){
  if(is.na(indexes)){
    point.est <- NA
    ci <- c(NA,NA)
    n.SV <- 0
  }else{
    ACs <- dat$AC[indexes]
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
    n.SV <- length(indexes)
  }
  return(c(point.est,ci,n.SV))
}

#Dotplot of fraction of singletons per functional annotation by class
dotplot.fracSingletons.general <- function(singleton.dat,class.colors,class.labels,
                                           add.counts=F,add.sizes=F,parmar=c(2.25,2.25,0.25,0.25),
                                           cex.xlabels=0.9,dashed.lines=NULL,ylims=NULL,
                                           ytck=-0.04,cex.ylabels=0.65,line.ylabel=1.5){
  singleton.dat <- apply(singleton.dat,2,as.numeric)
  par(mar=parmar,bty="n")
  if(is.null(ylims)){
    ylims <- range(as.numeric(singleton.dat[,1]),na.rm=T)
    ylims <- c(ylims[1]-(0.2*(ylims[2]-ylims[1])),
               ylims[2]+(0.2*(ylims[2]-ylims[1])))
  }
  plot(y=ylims,
       x=c(0.25,nrow(singleton.dat)-0.4),
       type="n",xaxt="n",xlab="",yaxt="n",ylab="")
  if(is.null(dashed.lines)){
    segments(x0=0.5,x1=nrow(singleton.dat)-0.5,
             y0=singleton.dat[1,1],y1=singleton.dat[1,1],
             lty=2)
  }else{
    lapply(dashed.lines,function(line.coords){
      segments(x0=as.numeric(line.coords[1])+0.5,
               x1=as.numeric(line.coords[2])-0.5,
               y0=as.numeric(line.coords[3]),
               y1=as.numeric(line.coords[3]),
               lty=2,col=line.coords[4])
    })
  }
  segments(y0=singleton.dat[,2],
           y1=singleton.dat[,3],
           x0=(1:nrow(singleton.dat))-0.5,
           x1=(1:nrow(singleton.dat))-0.5,
           lwd=2,lend="round",
           col=class.colors)
  points(y=singleton.dat[,1],
         x=c(1:nrow(singleton.dat))-0.5,
         pch=19,col=class.colors)
  if(add.counts==T){
    sapply(1:nrow(singleton.dat),function(i){
      if(!is.na(singleton.dat[i,4])){
        if(singleton.dat[i,4]>=1000){
          count.lab <- paste(round(singleton.dat[i,4]/1000,0),"k",sep="")
        }else{
          count.lab <- prettyNum(singleton.dat[i,4],big.mark=",")
        }
        class.labels[i] <- paste(class.labels[i],"\n",count.lab,sep="")
      }
    })
  }
  sapply(1:nrow(singleton.dat),function(i){
    axis(1,at=i-0.5,labels=class.labels[i],tick=F,line=-1.6,cex.axis=cex.xlabels,padj=1)
  })
  axis(2,at=axTicks(2),labels=NA,tck=ytck)
  sapply(axTicks(2),function(x){
    axis(2,at=x,line=-0.7,cex.axis=cex.ylabels,tick=F,las=2,
         labels=paste(round(x*100,0),"%",sep=""))
  })
  mtext(2,text="Singleton Proportion",line=line.ylabel,cex=0.8)
}

#Helper function to add linear fit to constraint decile analyses
helper.add.lm <- function(singleton.data,color,at=c(2.5,11.5)){
  lm.dat <- data.frame("d"=1:10,"p"=singleton.data[-c(1:2),1])
  y.ends <- predict(lm(p ~ d, data=lm.dat),newdata=data.frame("d"=c(1,10)))
  segments(x0=at[1],x1=at[2],
           y0=y.ends[1],y1=y.ends[2],
           lwd=4,col=adjustcolor(color,alpha=0.33))
}




################
###RSCRIPT BLOCK
################
require(optparse,quietly=T)
require(plotrix,quietly=T)
require(boot,quietly=T)
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
args <- parse_args(OptionParser(usage="%prog VCF2BED EXON_OVR WHOLEGENE_OVR INV_CENTROMERES INV_RECOMB ALL_RECOMB SD_SR_COV OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 10){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
vcf2bed.in <- args$args[1]
SNVdata.in <- args$args[2]
gene_metadata.in <- args$args[3]
exon_overlap.in <- args$args[4]
wholegene_overlap.in <- args$args[5]
inv_centromeres.in <- args$args[6]
inv_recomb.in <- args$args[7]
all_recomb.in <- args$args[8]
sd_sr_overlap.in <- args$args[9]
OUTDIR <- args$args[10]
prefix <- opts$prefix
svtypes.file <- opts$svtypes

# #Dev parameters (local)
# vcf2bed.in <- "~/scratch/gnomAD_v2_SV_MASTER.vcf2bed.bed.gz"
# SNVdata.in <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/gnomAD_v2.1_canonical_constraint.condensed.txt.gz"
# gene_metadata.in <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/Gencode.v19.autosomal_canonical_gene_metadata.txt.gz"
# exon_overlap.in <- "~/scratch/gnomAD_v2_SV_MASTER.exons_overlapped_per_SV.txt.gz"
# wholegene_overlap.in <- "~/scratch/gnomAD_v2_SV_MASTER.whole_genes_overlapped_per_SV.txt.gz"
# inv_centromeres.in <- "~/scratch/gnomAD_v2_SV_MASTER.centromeres_per_inversion.txt"
# inv_recomb.in <- "~/scratch/gnomAD_v2_SV_MASTER.recombination_hotspots_per_inversion.txt"
# all_recomb.in <- "~/scratch/gnomAD_v2_SV_MASTER.recombination_hotspots.all_sv.txt"
# sd_sr_overlap.in <- "~/scratch/gnomAD_v2_SV_MASTER.variant_SD_SR_coverage.txt.gz"
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
#Add columns for exon/wholegene/centromere/recomb overlap
exon_overlap <- read.table(exon_overlap.in,header=F,sep="\t")
colnames(exon_overlap) <- c("name","N_EXONS")
dat.all <- merge(dat.all,exon_overlap,on="name",sort=F,all.x=T)
rm(exon_overlap)
wholegene_overlap <- read.table(wholegene_overlap.in,header=F,sep="\t")
colnames(wholegene_overlap) <- c("name","N_WHOLEGENES")
dat.all <- merge(dat.all,wholegene_overlap,on="name",sort=F,all.x=T)
rm(wholegene_overlap)
inv_centromeres <- read.table(inv_centromeres.in,header=F,sep="\t")
colnames(inv_centromeres) <- c("name","N_CENTROMERES_SPANNED_INV")
dat.all <- merge(dat.all,inv_centromeres,on="name",sort=F,all.x=T)
rm(inv_centromeres)
inv_recomb <- read.table(inv_recomb.in,header=F,sep="\t")
colnames(inv_recomb) <- c("name","N_RECOMB_SPANNED_INV")
dat.all <- merge(dat.all,inv_recomb,on="name",sort=F,all.x=T)
rm(inv_recomb)
all_recomb <- read.table(all_recomb.in,header=F,sep="\t")
colnames(all_recomb) <- c("name","N_RECOMB_SPANNED_ALL")
dat.all <- merge(dat.all,all_recomb,on="name",sort=F,all.x=T)
rm(all_recomb)
sd_sr_overlap <- read.table(sd_sr_overlap.in,header=F,sep="\t")
colnames(sd_sr_overlap) <- c("name","SD_SR_OVERLAP")
dat.all <- merge(dat.all,sd_sr_overlap,on="name",sort=F,all.x=T)
rm(sd_sr_overlap)
#Filter to final analysis variants
dat <- dat.all[which(dat.all$FILTER=="PASS" | dat.all$FILTER=="MULTIALLELIC"),]
#Import SNV data
snv.data <- read.snvdata(SNVdata.in,gene_metadata.in)

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


#######################################################################
###Analysis of LOF deletions by count of exons and count of whole genes
#######################################################################
lof.dels.by.exon.count.idxs <- lapply(1:8,function(x){
  which(dat$SVTYPE=="DEL" 
        & dat$N_EXONS==x
        & dat$SD_SR_OVERLAP<0.3
        & !is.na(dat$PROTEIN_CODING__LOF)
        & !(dat$chrom %in% c("X","Y")))
})
lof.dels.by.exon.count.fracSingles <- do.call("rbind",lapply(lof.dels.by.exon.count.idxs,
                                                             calc.fracSingletons.single,dat=dat))

lof.dels.by.exon.count.fracSingles <- rbind(calc.fracSingletons.single(dat=dat,
                                                                       indexes=which(dat$SVTYPE=="DEL" 
                                                                                     & dat$SD_SR_OVERLAP<0.3
                                                                                     & !(dat$chrom %in% c("X","Y")))),
                                            c(NA,NA,NA,""),
                                            calc.fracSingletons.single(dat=dat,
                                                                       indexes=which(dat$SVTYPE=="DEL" 
                                                                                     & dat$SD_SR_OVERLAP<0.3
                                                                                     & dat$N_EXONS==0
                                                                                     & is.na(dat$PROTEIN_CODING__LOF)
                                                                                     & !(dat$chrom %in% c("X","Y")))),
                                            lof.dels.by.exon.count.fracSingles,
                                            calc.fracSingletons.single(dat=dat,
                                                                       indexes=which(dat$SVTYPE=="DEL" 
                                                                                     & dat$SD_SR_OVERLAP<0.3
                                                                                     & dat$N_EXONS>8
                                                                                     & !is.na(dat$PROTEIN_CODING__LOF)
                                                                                     & !(dat$chrom %in% c("X","Y")))))


###Analysis of LOF deletions by count of whole-genes
lof.dels.by.wholegene.count.idxs <- lapply(1:2,function(x){
  which(dat$SVTYPE=="DEL" 
        & dat$N_WHOLEGENES==x
        & dat$SD_SR_OVERLAP<0.3
        & !is.na(dat$PROTEIN_CODING__LOF)
        & !(dat$chrom %in% c("X","Y")))
})
lof.dels.by.wholegene.count.fracSingles <- do.call("rbind",lapply(lof.dels.by.wholegene.count.idxs,
                                                                  calc.fracSingletons.single,dat=dat))
lof.dels.by.wholegene.count.fracSingles <- rbind(calc.fracSingletons.single(dat=dat,
                                                                            indexes=which(dat$SVTYPE=="DEL" 
                                                                                          & dat$N_WHOLEGENES==0
                                                                                          & dat$SD_SR_OVERLAP<0.3
                                                                                          & is.na(dat$PROTEIN_CODING__LOF)
                                                                                          & !(dat$chrom %in% c("X","Y")))),
                                                 lof.dels.by.wholegene.count.fracSingles,
                                                 calc.fracSingletons.single(dat=dat,
                                                                            indexes=which(dat$SVTYPE=="DEL" 
                                                                                          & dat$N_WHOLEGENES>2
                                                                                          & dat$SD_SR_OVERLAP<0.3
                                                                                          & !is.na(dat$PROTEIN_CODING__LOF)
                                                                                          & !(dat$chrom %in% c("X","Y")))))


#Generate plot for both exon-level and gene-level plots
del.plot.dat <- rbind(lof.dels.by.exon.count.fracSingles,
                      c(NA,NA,NA,""),
                      lof.dels.by.wholegene.count.fracSingles)
pdf(paste(OUTDIR,"/",prefix,".frac_singletons.DEL_by_exon_count.pdf",sep=""),
    height=2,width=4.75)
dotplot.fracSingletons.general(singleton.dat=del.plot.dat,
                               class.colors=rep(svtypes$color[which(svtypes$svtype=="DEL")],
                                                length(del.plot.dat)),
                               class.labels=c("All\nDEL","",0:8,">8","",0:2,">2"),
                               add.counts=F)
# axis(1,at=13.8,labels=">10",tick=F,line=-0.9,cex.axis=0.9)
dev.off()


#############################
###Analysis of INV by context
#############################
inv.bycontext.idx <- list(which((dat$SVTYPE=="INV" | dat$SVTYPE=="CPX")),
                          which((dat$SVTYPE=="INV" | dat$SVTYPE=="CPX")
                                & is.na(dat$PROTEIN_CODING__LOF)
                                & is.na(dat$PROTEIN_CODING__COPY_GAIN)
                                & is.na(dat$PROTEIN_CODING__DUP_LOF)),
                          which((dat$SVTYPE=="INV" | dat$SVTYPE=="CPX")
                                & !is.na(dat$PROTEIN_CODING__LOF)),
                          which((dat$SVTYPE=="INV" | dat$SVTYPE=="CPX")
                                & !is.na(dat$PROTEIN_CODING__COPY_GAIN)),
                          which((dat$SVTYPE=="INV" | dat$SVTYPE=="CPX")
                                & dat$N_RECOMB_SPANNED_INV==0
                                & is.na(dat$PROTEIN_CODING__LOF)
                                & is.na(dat$PROTEIN_CODING__COPY_GAIN)),
                          which((dat$SVTYPE=="INV" | dat$SVTYPE=="CPX")
                                & dat$N_RECOMB_SPANNED_INV==1
                                & is.na(dat$PROTEIN_CODING__LOF)
                                & is.na(dat$PROTEIN_CODING__COPY_GAIN)),
                          which((dat$SVTYPE=="INV" | dat$SVTYPE=="CPX")
                                & dat$N_RECOMB_SPANNED_INV==2
                                & is.na(dat$PROTEIN_CODING__LOF)
                                & is.na(dat$PROTEIN_CODING__COPY_GAIN)),
                          which((dat$SVTYPE=="INV" | dat$SVTYPE=="CPX")
                                & dat$N_RECOMB_SPANNED_INV>2
                                & is.na(dat$PROTEIN_CODING__LOF)
                                & is.na(dat$PROTEIN_CODING__COPY_GAIN)))
inv.bycontext.singleton.data <- do.call("rbind", lapply(inv.bycontext.idx,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
inv.bycontext.singleton.data <- rbind(inv.bycontext.singleton.data[1,],
                                      c(NA,NA,NA,NA),
                                      inv.bycontext.singleton.data[3:4,],
                                      c(NA,NA,NA,NA),
                                      inv.bycontext.singleton.data[5:8,])
pdf(paste(OUTDIR,"/",prefix,".frac_singletons.inv_by_context.pdf",sep=""),
    height=2,width=2.4)
dotplot.fracSingletons.general(singleton.dat=inv.bycontext.singleton.data,
                               class.colors=rep(svtypes$color[which(svtypes$svtype=="INV")],
                                                length(inv.bycontext.singleton.data)),
                               class.labels=c("Any\nAny","",
                                              "LoF\nAny","CG\nAny","",
                                              "No\n0","No\n1","No\n2","No\n>2"),
                               add.counts=T,cex.xlabels=0.65,
                               parmar=c(2.25,2.25,0.25,0.25))
dev.off()


####################################################################
###Analyses of SV by functional context stratified by PTV o:e DECILE
####################################################################
#pLOF
lof.indexes.byConstraint <- get.func.indexes.byConstraint(dat,func="LOF")
lof.byconstraint.singleton.data <- do.call("rbind", lapply(lof.indexes.byConstraint,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#IED
ied.indexes.byConstraint <- get.func.indexes.byConstraint(dat,func="DUP_LOF")
ied.byconstraint.singleton.data <- do.call("rbind", lapply(ied.indexes.byConstraint,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#CG
cg.indexes.byConstraint <- get.func.indexes.byConstraint(dat,func="COPY_GAIN")
cg.byconstraint.singleton.data <- do.call("rbind", lapply(cg.indexes.byConstraint,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#Whole-gene inversion
wginv.indexes.byConstraint <- get.func.indexes.byConstraint(dat,func="INV_SPAN")
wginv.byconstraint.singleton.data <- do.call("rbind", lapply(wginv.indexes.byConstraint,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#Partial gene duplication
partialdup.indexes.byConstraint <- get.func.indexes.byConstraint(dat,func="DUP_PARTIAL")
partialdup.byconstraint.singleton.data <- do.call("rbind", lapply(partialdup.indexes.byConstraint,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#Intronic
intronic.indexes.byConstraint <- get.func.indexes.byConstraint(dat,func="INTRONIC")
intronic.byconstraint.singleton.data <- do.call("rbind", lapply(intronic.indexes.byConstraint,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#Promoter
promoter.indexes.byConstraint <- get.func.indexes.byConstraint(dat,func="PROMOTER")
promoter.byconstraint.singleton.data <- do.call("rbind", lapply(promoter.indexes.byConstraint,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#UTR
utr.indexes.byConstraint <- get.func.indexes.byConstraint(dat,func="UTR")
utr.byconstraint.singleton.data <- do.call("rbind", lapply(utr.indexes.byConstraint,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#Nearest TSS
nearest.indexes.byConstraint <- get.func.indexes.byConstraint(dat,func="NEAREST_TSS")
nearest.byconstraint.singleton.data <- do.call("rbind", lapply(nearest.indexes.byConstraint,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#Generate plot
byconstraint.parmar <- c(5,3,2,0.75)
byconstraint.cex.xlabels <- 1
byconstraint.ylims <- c(0.25,0.8)
byconstraint.cex.ylabels <- 0.8
byconstraint.line.ylabel <- 2
byconstraint.ytck <- -0.02
pdf(paste(OUTDIR,"/",prefix,".frac_singletons_by_constraint.pdf",sep=""),
    height=7,width=8)
par(mfrow=c(3,3))
#pLoF
dotplot.fracSingletons.general(singleton.dat=lof.byconstraint.singleton.data,
                               class.colors=rep(svtypes$color[which(svtypes$svtype=="DEL")],12),
                               class.labels=c("ALL","",1:10),
                               add.counts=F,cex.xlabels=byconstraint.cex.xlabels,
                               parmar=byconstraint.parmar,ylims=byconstraint.ylims,
                               ytck=byconstraint.ytck,cex.ylabels=byconstraint.cex.ylabels,
                               line.ylabel=byconstraint.line.ylabel)
mtext(3,line=0.1,text="pLoF",font=2)
helper.add.lm(singleton.data=lof.byconstraint.singleton.data,
              color=svtypes$color[which(svtypes$svtype=="DEL")])
#IED
dotplot.fracSingletons.general(singleton.dat=ied.byconstraint.singleton.data,
                               class.colors=rep(svtypes$color[which(svtypes$svtype=="MCNV")],12),
                               class.labels=c("ALL","",1:10),
                               add.counts=F,cex.xlabels=byconstraint.cex.xlabels,
                               parmar=byconstraint.parmar,ylims=byconstraint.ylims,
                               ytck=byconstraint.ytck,cex.ylabels=byconstraint.cex.ylabels,
                               line.ylabel=byconstraint.line.ylabel)
mtext(3,line=0.1,text="IED",font=2)
helper.add.lm(singleton.data=ied.byconstraint.singleton.data,
              color=svtypes$color[which(svtypes$svtype=="MCNV")])
#CG
dotplot.fracSingletons.general(singleton.dat=cg.byconstraint.singleton.data,
                               class.colors=rep(svtypes$color[which(svtypes$svtype=="DUP")],12),
                               class.labels=c("ALL","",1:10),
                               add.counts=F,cex.xlabels=byconstraint.cex.xlabels,
                               parmar=byconstraint.parmar,ylims=byconstraint.ylims,
                               ytck=byconstraint.ytck,cex.ylabels=byconstraint.cex.ylabels,
                               line.ylabel=byconstraint.line.ylabel)
mtext(3,line=0.1,text="CG",font=2)
helper.add.lm(singleton.data=cg.byconstraint.singleton.data,
              color=svtypes$color[which(svtypes$svtype=="DUP")])
#Whole-gene inversion
dotplot.fracSingletons.general(singleton.dat=wginv.byconstraint.singleton.data,
                               class.colors=rep(svtypes$color[which(svtypes$svtype=="INV")],12),
                               class.labels=c("ALL","",1:10),
                               add.counts=F,cex.xlabels=byconstraint.cex.xlabels,
                               parmar=byconstraint.parmar,ylims=byconstraint.ylims,
                               ytck=byconstraint.ytck,cex.ylabels=byconstraint.cex.ylabels,
                               line.ylabel=byconstraint.line.ylabel)
mtext(3,line=1.1,text="Whole-Gene INV",font=2)
mtext(3,line=0.1,text="(No Coding Effect)",font=3,cex=0.8)
helper.add.lm(singleton.data=wginv.byconstraint.singleton.data,
              color=svtypes$color[which(svtypes$svtype=="INV")])
#Partial gene duplication
dotplot.fracSingletons.general(singleton.dat=partialdup.byconstraint.singleton.data,
                               class.colors=rep(svtypes$color[which(svtypes$svtype=="DUP")],12),
                               class.labels=c("ALL","",1:10),
                               add.counts=F,cex.xlabels=byconstraint.cex.xlabels,
                               parmar=byconstraint.parmar,ylims=byconstraint.ylims,
                               ytck=byconstraint.ytck,cex.ylabels=byconstraint.cex.ylabels,
                               line.ylabel=byconstraint.line.ylabel)
mtext(3,line=1.1,text="Partial CG",font=2)
mtext(3,line=0.1,text="(No Coding Effect)",font=3,cex=0.8)
helper.add.lm(singleton.data=partialdup.byconstraint.singleton.data,
              color=svtypes$color[which(svtypes$svtype=="DUP")])
#Intronic
dotplot.fracSingletons.general(singleton.dat=intronic.byconstraint.singleton.data,
                               class.colors=rep("grey35",12),
                               class.labels=c("ALL","",1:10),
                               add.counts=F,cex.xlabels=byconstraint.cex.xlabels,
                               parmar=byconstraint.parmar,ylims=byconstraint.ylims,
                               ytck=byconstraint.ytck,cex.ylabels=byconstraint.cex.ylabels,
                               line.ylabel=byconstraint.line.ylabel)
mtext(3,line=1.1,text="Intronic",font=2)
mtext(3,line=0.1,text="(No Coding Effect)",font=3,cex=0.8)
helper.add.lm(singleton.data=intronic.byconstraint.singleton.data,
              color="gray35")
#Promoter
dotplot.fracSingletons.general(singleton.dat=promoter.byconstraint.singleton.data,
                               class.colors=rep("grey35",12),
                               class.labels=c("ALL","",1:10),
                               add.counts=F,cex.xlabels=byconstraint.cex.xlabels,
                               parmar=byconstraint.parmar,ylims=byconstraint.ylims,
                               ytck=byconstraint.ytck,cex.ylabels=byconstraint.cex.ylabels,
                               line.ylabel=byconstraint.line.ylabel)
mtext(3,line=1.1,text="Promoter",font=2)
mtext(3,line=0.1,text="(No Coding Effect)",font=3,cex=0.8)
helper.add.lm(singleton.data=promoter.byconstraint.singleton.data,
              color="gray35")
#UTR
dotplot.fracSingletons.general(singleton.dat=utr.byconstraint.singleton.data,
                               class.colors=rep("grey35",12),
                               class.labels=c("ALL","",1:10),
                               add.counts=F,cex.xlabels=byconstraint.cex.xlabels,
                               parmar=byconstraint.parmar,ylims=byconstraint.ylims,
                               ytck=byconstraint.ytck,cex.ylabels=byconstraint.cex.ylabels,
                               line.ylabel=byconstraint.line.ylabel)
mtext(3,line=1.1,text="UTR",font=2)
mtext(3,line=0.1,text="(No Coding Effect)",font=3,cex=0.8)
helper.add.lm(singleton.data=utr.byconstraint.singleton.data,
              color="gray35")
#Nearest Gene
dotplot.fracSingletons.general(singleton.dat=nearest.byconstraint.singleton.data,
                               class.colors=rep("grey35",12),
                               class.labels=c("ALL","",1:10),
                               add.counts=F,cex.xlabels=byconstraint.cex.xlabels,
                               parmar=byconstraint.parmar,ylims=byconstraint.ylims,
                               ytck=byconstraint.ytck,cex.ylabels=byconstraint.cex.ylabels,
                               line.ylabel=byconstraint.line.ylabel)
mtext(3,line=1.1,text="Nearest Gene",font=2)
mtext(3,line=0.1,text="(Strictly Intergenic)",font=3,cex=0.8)
helper.add.lm(singleton.data=nearest.byconstraint.singleton.data,
              color="gray35")
dev.off()


############################################################################
###Analyses of SV by functional context stratified by PTV o:e CLASSIFICATION
############################################################################
#pLOF
lof.indexes.byConstraintClass <- get.func.indexes.byConstraint(dat,func="LOF",bins=data.frame(c(0,15,85),c(15,85,100)))[3:5]
lof.byConstraintClass.singleton.data <- do.call("rbind", lapply(lof.indexes.byConstraintClass,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#IED
ied.indexes.byConstraintClass <- get.func.indexes.byConstraint(dat,func="DUP_LOF",bins=data.frame(c(0,15,85),c(15,85,100)))[3:5]
ied.byConstraintClass.singleton.data <- do.call("rbind", lapply(ied.indexes.byConstraintClass,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#CG
cg.indexes.byConstraintClass <- get.func.indexes.byConstraint(dat,func="COPY_GAIN",bins=data.frame(c(0,15,85),c(15,85,100)))[3:5]
cg.byConstraintClass.singleton.data <- do.call("rbind", lapply(cg.indexes.byConstraintClass,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#Whole-gene inversion
wginv.indexes.byConstraintClass <- get.func.indexes.byConstraint(dat,func="INV_SPAN",bins=data.frame(c(0,15,85),c(15,85,100)))[3:5]
wginv.byConstraintClass.singleton.data <- do.call("rbind", lapply(wginv.indexes.byConstraintClass,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#Partial gene duplication
partialdup.indexes.byConstraintClass <- get.func.indexes.byConstraint(dat,func="DUP_PARTIAL",bins=data.frame(c(0,15,85),c(15,85,100)))[3:5]
partialdup.byConstraintClass.singleton.data <- do.call("rbind", lapply(partialdup.indexes.byConstraintClass,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#Intronic
intronic.indexes.byConstraintClass <- get.func.indexes.byConstraint(dat,func="INTRONIC",bins=data.frame(c(0,15,85),c(15,85,100)))[3:5]
intronic.byConstraintClass.singleton.data <- do.call("rbind", lapply(intronic.indexes.byConstraintClass,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#Promoter
promoter.indexes.byConstraintClass <- get.func.indexes.byConstraint(dat,func="PROMOTER",bins=data.frame(c(0,15,85),c(15,85,100)))[3:5]
promoter.byConstraintClass.singleton.data <- do.call("rbind", lapply(promoter.indexes.byConstraintClass,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#UTR
utr.indexes.byConstraintClass <- get.func.indexes.byConstraint(dat,func="UTR",bins=data.frame(c(0,15,85),c(15,85,100)))[3:5]
utr.byConstraintClass.singleton.data <- do.call("rbind", lapply(utr.indexes.byConstraintClass,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#Nearest TSS
nearest.indexes.byConstraintClass <- get.func.indexes.byConstraint(dat,func="NEAREST_TSS",bins=data.frame(c(0,15,85),c(15,85,100)))[3:5]
nearest.byConstraintClass.singleton.data <- do.call("rbind", lapply(nearest.indexes.byConstraintClass,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#Merge data
singleton.data.byConstraintClass <- rbind(lof.byConstraintClass.singleton.data,
                                          c(NA,NA,NA,NA),
                                          ied.byConstraintClass.singleton.data,
                                          c(NA,NA,NA,NA),
                                          cg.byConstraintClass.singleton.data,
                                          c(NA,NA,NA,NA),
                                          partialdup.byConstraintClass.singleton.data,
                                          c(NA,NA,NA,NA),
                                          intronic.byConstraintClass.singleton.data,
                                          c(NA,NA,NA,NA),
                                          promoter.byConstraintClass.singleton.data,
                                          c(NA,NA,NA,NA),
                                          utr.byConstraintClass.singleton.data,
                                          c(NA,NA,NA,NA),
                                          nearest.byConstraintClass.singleton.data)
#Generate plot
pdf(paste(OUTDIR,"/",prefix,".frac_singletons_by_constraint_category.pdf",sep=""),
    height=2.25,width=8)
dotplot.fracSingletons.general(singleton.dat=singleton.data.byConstraintClass,
                               class.colors=c(rep(svtypes$color[which(svtypes=="DEL")],4),
                                              rep(svtypes$color[which(svtypes=="MCNV")],4),
                                              rep(svtypes$color[which(svtypes=="DUP")],4),
                                              rep(svtypes$color[which(svtypes=="DUP")],4),
                                              rep("grey35",16)),
                               class.labels=rep(c("T","M","B",NA),9),
                               add.counts=T,cex.xlabels=0.8,
                               parmar=c(1.5,2.25,0.25,0.25),dashed.lines=F,
                               ytck=-0.02)
dev.off()




#####################################################################
###Analyses of SV by functional context stratified by PTV o:e SEXTILE
#####################################################################
#pLOF
lof.indexes.byConstraintSextile <- get.func.indexes.byConstraint(dat,func="LOF",bins=data.frame(seq(0,(5/6)*100,(1/6)*100),seq((1/6)*100,100,(1/6)*100)))[-c(1:2)]
lof.byConstraintSextile.singleton.data <- do.call("rbind", lapply(lof.indexes.byConstraintSextile,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#IED
ied.indexes.byConstraintSextile <- get.func.indexes.byConstraint(dat,func="DUP_LOF",bins=data.frame(seq(0,(5/6)*100,(1/6)*100),seq((1/6)*100,100,(1/6)*100)))[-c(1:2)]
ied.byConstraintSextile.singleton.data <- do.call("rbind", lapply(ied.indexes.byConstraintSextile,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#CG
cg.indexes.byConstraintSextile <- get.func.indexes.byConstraint(dat,func="COPY_GAIN",bins=data.frame(seq(0,(5/6)*100,(1/6)*100),seq((1/6)*100,100,(1/6)*100)))[-c(1:2)]
cg.byConstraintSextile.singleton.data <- do.call("rbind", lapply(cg.indexes.byConstraintSextile,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#Whole-gene inversion
wginv.indexes.byConstraintSextile <- get.func.indexes.byConstraint(dat,func="INV_SPAN",bins=data.frame(seq(0,(5/6)*100,(1/6)*100),seq((1/6)*100,100,(1/6)*100)))[-c(1:2)]
wginv.byConstraintSextile.singleton.data <- do.call("rbind", lapply(wginv.indexes.byConstraintSextile,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#Partial gene duplication
partialdup.indexes.byConstraintSextile <- get.func.indexes.byConstraint(dat,func="DUP_PARTIAL",bins=data.frame(seq(0,(5/6)*100,(1/6)*100),seq((1/6)*100,100,(1/6)*100)))[-c(1:2)]
partialdup.byConstraintSextile.singleton.data <- do.call("rbind", lapply(partialdup.indexes.byConstraintSextile,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#Intronic
intronic.indexes.byConstraintSextile <- get.func.indexes.byConstraint(dat,func="INTRONIC",bins=data.frame(seq(0,(5/6)*100,(1/6)*100),seq((1/6)*100,100,(1/6)*100)))[-c(1:2)]
intronic.byConstraintSextile.singleton.data <- do.call("rbind", lapply(intronic.indexes.byConstraintSextile,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#Promoter
promoter.indexes.byConstraintSextile <- get.func.indexes.byConstraint(dat,func="PROMOTER",bins=data.frame(seq(0,(5/6)*100,(1/6)*100),seq((1/6)*100,100,(1/6)*100)))[-c(1:2)]
promoter.byConstraintSextile.singleton.data <- do.call("rbind", lapply(promoter.indexes.byConstraintSextile,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#UTR
utr.indexes.byConstraintSextile <- get.func.indexes.byConstraint(dat,func="UTR",bins=data.frame(seq(0,(5/6)*100,(1/6)*100),seq((1/6)*100,100,(1/6)*100)))[-c(1:2)]
utr.byConstraintSextile.singleton.data <- do.call("rbind", lapply(utr.indexes.byConstraintSextile,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#Nearest TSS
nearest.indexes.byConstraintSextile <- get.func.indexes.byConstraint(dat,func="NEAREST_TSS",bins=data.frame(seq(0,(5/6)*100,(1/6)*100),seq((1/6)*100,100,(1/6)*100)))[-c(1:2)]
nearest.byConstraintSextile.singleton.data <- do.call("rbind", lapply(nearest.indexes.byConstraintSextile,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
#All SV
autosomal.singleton.data <- calc.fracSingletons.single(dat,indexes=which(!(dat$chrom %in% c("X","Y"))))
#Merge data
singleton.data.byConstraintSextile <- rbind(autosomal.singleton.data,
                                            c(NA,NA,NA,NA),
                                            lof.byConstraintSextile.singleton.data,
                                          c(NA,NA,NA,NA),
                                          ied.byConstraintSextile.singleton.data,
                                          c(NA,NA,NA,NA),
                                          cg.byConstraintSextile.singleton.data,
                                          c(NA,NA,NA,NA),
                                          partialdup.byConstraintSextile.singleton.data,
                                          c(NA,NA,NA,NA),
                                          intronic.byConstraintSextile.singleton.data,
                                          c(NA,NA,NA,NA),
                                          promoter.byConstraintSextile.singleton.data,
                                          c(NA,NA,NA,NA),
                                          utr.byConstraintSextile.singleton.data,
                                          c(NA,NA,NA,NA),
                                          nearest.byConstraintSextile.singleton.data)
#Generate plot
pdf(paste(OUTDIR,"/",prefix,".frac_singletons_by_constraint_sextile.pdf",sep=""),
    height=1.05*2,width=1.05*8.75)
dotplot.fracSingletons.general(singleton.dat=singleton.data.byConstraintSextile,
                               class.colors=c("black",NA,
                                              rep(svtypes$color[which(svtypes=="DEL")],7),
                                              rep(svtypes$color[which(svtypes=="MCNV")],7),
                                              rep(svtypes$color[which(svtypes=="DUP")],7),
                                              rep(svtypes$color[which(svtypes=="DUP")],7),
                                              rep("grey35",28)),
                               class.labels=c(NA,NA,rep(c(1:6,NA),8)),
                               add.counts=T,cex.xlabels=0.75,
                               parmar=c(1.5,2.25,1.75,0.25),
                               ytck=-0.02,ylims=c(0.3,0.7))
axis(3,tick=F,at=seq(0,7*7,by=7)+5,line=mean(c(-0.2,-0.85)),
     labels=c("pLoF","IED","CG",rep("",5)),
     cex.axis=0.85)
axis(3,tick=F,at=seq(0,7*7,by=7)+5,line=-0.2,
     labels=c(rep("",3),"Partial CG","Intronic","Promoter","UTR","Nearest Gene"),
     cex.axis=0.85)
axis(3,tick=F,at=seq(0,7*7,by=7)+5,line=-0.85,
     labels=c(rep("",3),rep("(No Coding Effect)",4),"(Strictly Intergenic)"),
     cex.axis=0.6,col.axis="grey30")
sapply(seq(0,7*7,by=7)+2,function(i){
  axis(3,labels=NA,at=c(i+0.5,i+5.5),tck=0,line=0.1)
})
dev.off()


#############################################
###Analysis of CNVs at recombination hotspots
#############################################
cnv.recomb.idx <- list(which(dat$SVTYPE=="DEL"),
                       which(dat$SVTYPE=="DEL"
                             & !is.na(dat$PROTEIN_CODING__LOF)),
                       which(dat$SVTYPE=="DEL"
                             & dat$N_RECOMB_SPANNED_ALL==0
                             & is.na(dat$PROTEIN_CODING__LOF)),
                       which(dat$SVTYPE=="DEL"
                             & dat$N_RECOMB_SPANNED_ALL==1
                             & is.na(dat$PROTEIN_CODING__LOF)),
                       which(dat$SVTYPE=="DEL"
                             & dat$N_RECOMB_SPANNED_ALL==2
                             & is.na(dat$PROTEIN_CODING__LOF)),
                       which(dat$SVTYPE=="DEL"
                             & dat$N_RECOMB_SPANNED_ALL>2
                             & is.na(dat$PROTEIN_CODING__LOF)),
                       which(dat$SVTYPE=="DUP"),
                       which(dat$SVTYPE=="DUP"
                             & !is.na(dat$PROTEIN_CODING__DUP_LOF)),
                       which(dat$SVTYPE=="DUP"
                             & !is.na(dat$PROTEIN_CODING__COPY_GAIN)),
                       which(dat$SVTYPE=="DUP"
                             & dat$N_RECOMB_SPANNED_ALL==0
                             & is.na(dat$PROTEIN_CODING__DUP_LOF)
                             & is.na(dat$PROTEIN_CODING__COPY_GAIN)),
                       which(dat$SVTYPE=="DUP"
                             & dat$N_RECOMB_SPANNED_ALL==1
                             & is.na(dat$PROTEIN_CODING__DUP_LOF)
                             & is.na(dat$PROTEIN_CODING__COPY_GAIN)),
                       which(dat$SVTYPE=="DUP"
                             & dat$N_RECOMB_SPANNED_ALL==2
                             & is.na(dat$PROTEIN_CODING__DUP_LOF)
                             & is.na(dat$PROTEIN_CODING__COPY_GAIN)),
                       which(dat$SVTYPE=="DUP"
                             & dat$N_RECOMB_SPANNED_ALL>2
                             & is.na(dat$PROTEIN_CODING__DUP_LOF)
                             & is.na(dat$PROTEIN_CODING__COPY_GAIN)))
cnv.recomb.singleton.data <- do.call("rbind", lapply(cnv.recomb.idx,function(indexes){
  calc.fracSingletons.single(dat=dat,indexes=indexes)
}))
del.recomb.singleton.plot.data <- rbind(cnv.recomb.singleton.data[1,],
                                        c(NA,NA,NA,NA),
                                        cnv.recomb.singleton.data[2,],
                                        c(NA,NA,NA,NA),
                                        cnv.recomb.singleton.data[3:6,])
pdf(paste(OUTDIR,"/",prefix,".frac_singletons.del_by_recomb_context.pdf",sep=""),
    height=2.25,width=2.25)
dotplot.fracSingletons.general(singleton.dat=del.recomb.singleton.plot.data,
                               class.colors=c(rep(svtypes$color[which(svtypes$svtype=="DEL")],8)),
                               class.labels=c("Any\nAny","",
                                              "LoF\nAny","",
                                              "No\n0","No\n1","No\n2","No\n>2"),
                               add.counts=T,cex.xlabels=0.65,
                               parmar=c(2.5,2.25,0.25,0.25))
dev.off()
dup.recomb.singleton.plot.data <- rbind(cnv.recomb.singleton.data[7,],
                                        c(NA,NA,NA,NA),
                                        cnv.recomb.singleton.data[8:9,],
                                        c(NA,NA,NA,NA),
                                        cnv.recomb.singleton.data[10:13,])
pdf(paste(OUTDIR,"/",prefix,".frac_singletons.dup_by_recomb_context.pdf",sep=""),
    height=2,width=2.4)
dotplot.fracSingletons.general(singleton.dat=dup.recomb.singleton.plot.data,
                               class.colors=c(rep(svtypes$color[which(svtypes$svtype=="DUP")],10)),
                               class.labels=c("Any\nAny","",
                                              "IED\nAny","CG\nAny","",
                                              "No\n0","No\n1","No\n2","No\n>2"),
                               add.counts=T,cex.xlabels=0.65,
                               parmar=c(2.25,2.25,0.25,0.25))
dev.off()


############################################################
###Joint plot of INV, DEL, and DUP vs recombination hotspots
############################################################
combined.recomb.plot.data <- rbind(inv.bycontext.singleton.data,
                                   rep(NA,4),
                                   del.recomb.singleton.plot.data,
                                   rep(NA,4),
                                   dup.recomb.singleton.plot.data)
pdf(paste(OUTDIR,"/",prefix,".frac_singletons.inv_del_dup_by_recomb_context.pdf",sep=""),
    height=2.5,width=5.6)
dotplot.fracSingletons.general(singleton.dat=combined.recomb.plot.data,
                               class.colors=c(rep(svtypes$color[which(svtypes$svtype=="INV")],
                                                  nrow(inv.bycontext.singleton.data)+1),
                                              rep(svtypes$color[which(svtypes$svtype=="DEL")],
                                                  nrow(del.recomb.singleton.plot.data)+1),
                                              rep(svtypes$color[which(svtypes$svtype=="DUP")],
                                                  nrow(dup.recomb.singleton.plot.data))),
                               class.labels=c("Any\nAny","",
                                              "LoF\nAny","CG\nAny","",
                                              "No\n0","No\n1","No\n2","No\n>2",
                                              "",
                                              "Any\nAny","",
                                              "LoF\nAny","",
                                              "No\n0","No\n1","No\n2","No\n>2",
                                              "",
                                              "Any\nAny","",
                                              "IED\nAny","CG\nAny","",
                                              "No\n0","No\n1","No\n2","No\n>2"),
                               add.counts=F,cex.xlabels=0.55,
                               parmar=c(1.75,2.25,0.25,0.25),
                               dashed.lines=list(c(0,nrow(inv.bycontext.singleton.data),
                                                   inv.bycontext.singleton.data[1,1],
                                                   svtypes$color[which(svtypes$svtype=="INV")]),
                                                 c(nrow(inv.bycontext.singleton.data)+1,
                                                   nrow(inv.bycontext.singleton.data)+1+nrow(del.recomb.singleton.plot.data),
                                                   del.recomb.singleton.plot.data[1,1],
                                                   svtypes$color[which(svtypes$svtype=="DEL")]),
                                                 c(nrow(inv.bycontext.singleton.data)+1+nrow(del.recomb.singleton.plot.data)+1,
                                                   nrow(inv.bycontext.singleton.data)+1+nrow(del.recomb.singleton.plot.data)+1+nrow(dup.recomb.singleton.plot.data),
                                                   dup.recomb.singleton.plot.data[1,1],
                                                   svtypes$color[which(svtypes$svtype=="DUP")])))
dev.off()

