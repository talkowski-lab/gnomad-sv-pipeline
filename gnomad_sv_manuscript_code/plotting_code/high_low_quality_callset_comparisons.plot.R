#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Compare a few key results at different quality thresholds


### Set master parameters
options(stringsAsFactors=F, scipen=1000)


####################
### HELPER FUNCTIONS
####################
# Read & clean vcf2bed input
read.vcf2bed <- function(vcf2bed.in, sd_sr_cov.in=NULL, quals.in=NULL, aps.in=NULL){
  # Read data
  dat <- read.table(vcf2bed.in, comment.char="", header=T, sep="\t")
  colnames(dat)[1] <- "chrom"
  # Restrict to sites with at least one observed alternative allele
  dat <- dat[union(which(dat$AC>0), grep("MULTIALLELIC", dat$FILTER, fixed=T)), ]
  # Drop columns not being used (to save memory)
  cols.to.drop <- c("start", "end", "CHR2", "CPX_INTERVALS", 
                    "END", "SOURCE", "STRANDS", "UNRESOLVED_TYPE", 
                    "LINCRNA__LOF", "LINCRNA__DUP_LOF", "LINCRNA__COPY_GAIN", 
                    "LINCRNA__DUP_PARTIAL", "LINCRNA__MSV_EXON_OVR", 
                    "LINCRNA__INTRONIC", "LINCRNA__INV_SPAN", "LINCRNA__UTR")
  dat <- dat[, -which(colnames(dat) %in% cols.to.drop)]
  # Convert numeric columns
  numeric.columns <- sort(unique(c(grep("FREQ", colnames(dat), fixed=T), 
                                   grep("AN", colnames(dat), fixed=T), 
                                   grep("AC", colnames(dat), fixed=T), 
                                   grep("AF", colnames(dat), fixed=T))))
  numeric.columns <- setdiff(numeric.columns, grep("SPAN", colnames(dat), fixed=T))
  dat[, numeric.columns] <- apply(dat[, numeric.columns], 2, as.numeric)
  # Read & add SDSR cov, if optioned
  if(!is.null(sd_sr_cov.in)){
    sdsr_cov <- read.table(sd_sr_cov.in, header=F, sep="\t")
    colnames(sdsr_cov) <- c("name", "SDSR_COV")
    dat <- merge(dat, sdsr_cov, by="name", all.x=T, all.y=F, sort=F)
  }
  # Read & add QUAL, if optioned
  if(!is.null(quals.in)){
    quals <- read.table(quals.in, header=T, sep="\t", comment.char="")
    colnames(quals)[1] <- "name"
    dat <- merge(dat, quals, by="name", all.x=T, all.y=F, sort=F)
  }
  # Read & add APS, if optioned
  if(!is.null(aps.in)){
    aps <- read.table(aps.in, header=T, sep="\t", comment.char="")
    colnames(aps)[1] <- "VID"
    dat <- merge(dat, aps, by.x="name", by.y="VID", all.x=T, all.y=F, sort=F)
    dat$APS[which(dat$chrom %in% c("X", "Y"))] <- NA
  }
  return(dat)
}

# Process SNV gene data
read.snvdata <- function(SNVdata.in, gene.metadata.in){
  #Read & clean
  snv.data <- read.table(SNVdata.in, header=T, comment.char="")
  metadata <- read.table(gene.metadata.in, header=T)
  merged <- merge(x=snv.data, y=metadata, by="gene", sort=F)
  merged <- merged[which(!(merged$chrom %in% c("chrX", "chrY"))), ]
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

# Organize table of genes from SV data
getSVdat <- function(dat, genes, prefix=NULL){
  #pLoF
  lof.any <- as.character(unlist(strsplit(dat$PROTEIN_CODING__LOF[which(!is.na(dat$PROTEIN_CODING__LOF))], split=",")))
  lof.del <- as.character(unlist(strsplit(dat$PROTEIN_CODING__LOF[which(!is.na(dat$PROTEIN_CODING__LOF)
                                                                        & dat$SVTYPE=="DEL")], split=",")))
  lof.other <- as.character(unlist(strsplit(dat$PROTEIN_CODING__LOF[which(!is.na(dat$PROTEIN_CODING__LOF)
                                                                          & dat$SVTYPE!="DEL")], split=",")))
  #Dup GC
  cg.dup <- as.character(unlist(strsplit(dat$PROTEIN_CODING__COPY_GAIN[which(!is.na(dat$PROTEIN_CODING__COPY_GAIN))], split=",")))
  #Dup pLoF
  plof.dup <- as.character(unlist(strsplit(dat$PROTEIN_CODING__DUP_LOF[which(!is.na(dat$PROTEIN_CODING__DUP_LOF))], split=",")))
  #Inversion span
  inv.span <- as.character(unlist(strsplit(dat$PROTEIN_CODING__INV_SPAN[which(!is.na(dat$PROTEIN_CODING__INV_SPAN))], split=",")))
  #Collect vector per gene
  res <- as.data.frame(t(sapply(genes, function(gene){
    g.lof.any <- length(which(lof.any==gene))
    g.lof.del <- length(which(lof.del==gene))
    g.lof.other <- length(which(lof.other==gene))
    g.cg.dup <- length(which(cg.dup==gene))
    g.plof.dup <- length(which(plof.dup==gene))
    g.inv.span <- length(which(inv.span==gene))
    g.out <- as.integer(c(g.lof.any, g.lof.del, g.lof.other, g.cg.dup, g.plof.dup, g.inv.span))
    g.out[which(is.na(g.out))] <- 0
    g.out <- c(gene, g.out)
    return(g.out)
  })))
  colnames(res) <- c("gene", "lof.any", "lof.del", "lof.other", "cg", "plof", "inv")
  if(!is.null(prefix)){
    colnames(res)[-1] <- paste(prefix, colnames(res)[-1], sep=".")
  }
  rownames(res) <- 1:nrow(res)
  return(res)
}

# Get merged table of SNV stats and SV stats by frequency bin
getSVdat.all <- function(dat, snv.data, include.cpx=T, require.SR=F){
  #Gather data
  genes <- sort(unique(as.character(snv.data$gene)))
  if(include.cpx==F){
    dat <- dat[which(dat$SVTYPE != "CPX"), ]
  }
  if(require.SR==T){
    dat <- dat[grep("SR", dat$EVIDENCE), ]
  }
  all.counts <- getSVdat(dat, genes, prefix="all")
  common.counts <- getSVdat(dat[which(dat$AF>=0.01), ], genes, prefix="common")
  rare.counts <- getSVdat(dat[which(dat$AF<0.01 & dat$AC>1), ], genes, prefix="rare")
  singleton.counts <- getSVdat(dat[which(dat$AC==1), ], genes, prefix="singleton")
  #Merge data
  merged <- merge(x=snv.data, y=all.counts, by="gene", all.x=T, sort=F)
  merged <- merge(x=merged, y=common.counts, by="gene", all.x=T, sort=F)
  merged <- merge(x=merged, y=rare.counts, by="gene", all.x=T, sort=F)
  merged <- merge(x=merged, y=singleton.counts, by="gene", all.x=T, sort=F)
  merged[, -c(which(colnames(merged) %in% c("gene", "chrom")))] <- apply(merged[, -c(which(colnames(merged) %in% c("gene", "chrom")))], 2, as.numeric)
  return(merged)
}

# Gather Hardy-Weinberg data
makeHWEmat <- function(dat, pop=NULL){
  sub.dat <- dat[which(!(dat$chrom %in% c("X", "Y"))), ]
  cols.to.exclude <- grep("MULTIALLELIC", sub.dat$FILTER, fixed=T)
  if(length(cols.to.exclude)>0){
    sub.dat <- sub.dat[-cols.to.exclude, ]
  }
  if(!is.null(pop)){
    n.genos.idx <- which(colnames(sub.dat)==paste(pop, "_N_BI_GENOS", sep=""))
    n.homref.idx <- which(colnames(sub.dat)==paste(pop, "_N_HOMREF", sep=""))
    n.het.idx <- which(colnames(sub.dat)==paste(pop, "_N_HET", sep=""))
    n.homalt.idx <- which(colnames(sub.dat)==paste(pop, "_N_HOMALT", sep=""))
  }else{
    n.genos.idx <- which(colnames(sub.dat)=="N_BI_GENOS")
    n.homref.idx <- which(colnames(sub.dat)=="N_HOMREF")
    n.het.idx <- which(colnames(sub.dat)=="N_HET")
    n.homalt.idx <- which(colnames(sub.dat)=="N_HOMALT")
  }
  sub.dat <- sub.dat[which(sub.dat[, n.genos.idx]>0), ]
  sub.dat <- sub.dat[which(apply(sub.dat[, c(n.het.idx, n.homalt.idx)], 1, sum)>0), ]
  HWE.mat <- data.frame("AA"=as.numeric(sub.dat[, n.homref.idx]), 
                        "AB"=as.numeric(sub.dat[, n.het.idx]), 
                        "BB"=as.numeric(sub.dat[, n.homalt.idx]))
  HWE.mat <- HWE.mat[complete.cases(HWE.mat), ]
  return(HWE.mat)
}

#Individual function to get mutation rate estimate for a single pop & svtype
getMu <- function(dat, pop=NULL, svtype=NULL, Ne=10000){
  # Format variables
  if(!is.null(pop)){
    pop <- paste(pop, "_", sep="")
  }
  # Get number of chromosomes assessed
  pop.AN.idx <- which(colnames(dat)==paste(pop, "N_BI_GENOS", sep=""))
  n <- 2*max(dat[, pop.AN.idx], na.rm=T)
  # Get number of autosomal biallelic segregating sites
  pop.AF.idx <- which(colnames(dat)==paste(pop, "AF", sep=""))
  if(!is.null(svtype)){
    K <- length(which(dat[, pop.AF.idx]>0 & dat$SVTYPE==svtype & !(dat$chrom %in% c("X", "Y")) & dat$SVTYPE != "MCNV"))
  }else{
    K <- length(which(dat[, pop.AF.idx]>0 & !(dat$chrom %in% c("X", "Y")) & dat$SVTYPE != "MCNV"))
  }
  # Get n-1th harmonic number
  harmsum <- sum(sapply(1:(n-1), function(k){1/k}))
  # Get Watterson estimator
  theta.hat <- K/harmsum
  # Solve for mutation rate
  mu <- theta.hat/(4*Ne)
  return(mu)
}

# Wrapper function to get mutation rate estimates for all svtypes & populations
getAllMus <- function(dat, Ne){
  mu.svtypes <- c("DEL", "DUP", "INS", "INV", "CPX")
  # Get mutation rates across classes
  mu.AllClasses <- sapply(mu.svtypes, function(svtype){
    getMu(dat, pop=NULL, svtype=svtype, Ne=Ne)
  })
  # Get mutation rate of all SV across populations
  mu.AllPops <- sapply(pops$pop, function(pop){
    getMu(dat, pop=pop, svtype=NULL, Ne=Ne)
  })
  # Get mutation rate by svtype & population
  mu.ClassByPop <- sapply(pops$pop, function(pop){
    sapply(mu.svtypes, function(svtype){
      getMu(dat, pop=pop, svtype=svtype, Ne=Ne)
    })
  })
  # Get mean by class across pops
  mu.PopMeanByClass <- apply(mu.ClassByPop, 1, function(vals){
    lower <- as.numeric(t.test(vals)$conf.int[1])
    mean <- as.numeric(t.test(vals)$estimate)
    upper <- as.numeric(t.test(vals)$conf.int[2])
    return(c(lower, mean, upper))
  })
  mu.PopMeanByClass <- cbind(apply(mu.PopMeanByClass, 1, sum), mu.PopMeanByClass)
  colnames(mu.PopMeanByClass)[1] <- "ALL"
  rownames(mu.PopMeanByClass) <- c("CI.lowerBound", "mean", "CI.upperBound")
  return(list("muByClass"=mu.AllClasses, 
              "muByPop"=mu.AllPops, 
              "muByClassAndPop"=mu.ClassByPop, 
              "mu.Means"=mu.PopMeanByClass))
}

#Prep covariates matrics for regression analysis
prepCovariates <- function(gene.data, y){
  cov <- data.frame("gene"=gene.data$gene, 
                    "y.SV_count"=gene.data[which(colnames(gene.data)==y)], 
                    "gene_length"=scale(log10(as.numeric(gene.data$gene_length)), center=T, scale=T), 
                    "exon_count"=scale(log10(as.numeric(gene.data$exon_count)), center=T, scale=T), 
                    "exon_median"=scale(log10(as.numeric(gene.data$exon_median)), center=T, scale=T), 
                    "exon_sum"=scale(log10(as.numeric(gene.data$exon_mean)), center=T, scale=T), 
                    "intron_count"=scale(log10(as.numeric(gene.data$intron_count)), center=T, scale=T), 
                    "intron_median"=scale(log10(as.numeric(gene.data$intron_median)), center=T, scale=T), 
                    "intron_sum"=scale(log10(as.numeric(gene.data$intron_mean)), center=T, scale=T), 
                    "segdup"=round(gene.data$segdup, 0))
  colnames(cov)[2] <- "y.SV_count"
  chrom.dummy.mat <- as.data.frame(sapply(unique(gene.data$chrom), function(chr){
    v <- rep(0, times=nrow(gene.data))
    v[which(gene.data$chrom==chr)] <- 1
    return(v)
  }))
  cov <- cbind(cov, chrom.dummy.mat)
  # rownames(cov) <- gene.data$gene
  return(cov)
}

#Fit model for predicting # of rare functional SV
fitConstraintModel <- function(cov, train.genes=which(gene.data$ptv_oe_dec>=5 & gene.data$ptv_oe_dec<=9)){
  #Fit negative binomial model
  glm.fit <- glm.nb(y.SV_count ~ ., data=cov[train.genes, -1])
  #Apply model to all genes
  glm.fit.vals <- predict.glm(glm.fit, newdata=cov[, -(1:2)], type="response")
  #Return output data frame
  res <- data.frame("gene"=cov$gene, 
                    "SV_count_raw"=cov$y.SV_count, 
                    "SV_count_glm_exp"=glm.fit.vals)
  return(res)
}

# Read external gene list and restrict to autosomal genes considered in gnomAD
importGenelist <- function(genelist.dir, filename, gene.data){
  genes <- read.table(paste(genelist.dir, "/", filename, sep=""), header=F)
  genes <- as.character(genes[, 1])
  genes <- genes[which(genes %in% gene.data$gene)]
  return(genes)
}

# Gather fraction of individuals from a given population with a rare pLoF in a given gene list
countRareLoFCarrierRate <- function(pop="ALL", genelist, dat, maxAF=0.001, mode="ALL"){
  # Get index of all variants with pLoF of at least one gene in list
  lof.in.genelist <- which(unlist(lapply(strsplit(dat$PROTEIN_CODING__LOF, split=","), function(genes){any(genes %in% genelist)})))
  subdat <- dat[lof.in.genelist, ]
  # Set filtering indexes
  if(pop=="ALL"){
    prefix <- ""
  }else{
    prefix <- paste(pop, "_", sep="")
  }
  AF.idx <- which(colnames(subdat)==paste(prefix, "AF", sep=""))
  genos.idx <- which(colnames(subdat)==paste(prefix, "N_BI_GENOS", sep=""))
  het.idx <- which(colnames(subdat)==paste(prefix, "N_HET", sep=""))
  homalt.idx <- which(colnames(subdat)==paste(prefix, "N_HOMALT", sep=""))
  # Count non-ref individuals by SVTYPE
  svtypes.forCarrierAnalysis <- c("DEL", "DUP", "INS", "INV", "CPX")
  fracs <- sapply(svtypes.forCarrierAnalysis, function(svtype){
    hits <- subdat[which(subdat[, AF.idx]<maxAF & subdat$SVTYPE==svtype), ]
    if(mode=="ALL"){
      sum(hits[, het.idx]+hits[, homalt.idx])/max(hits[, genos.idx], na.rm=T)
    }else if(mode=="HET"){
      sum(hits[, het.idx])/max(hits[, genos.idx], na.rm=T)
    }else if(mode=="HOM"){
      sum(hits[, homalt.idx])/max(hits[, genos.idx], na.rm=T)
    }else{
      stop("mode must be one of HOM, HET, ALL")
    }
  })
  return(fracs)
}

# Get table of gross chromosomal abnormalities
gather.bca.table <- function(dat, max.AF=0.01){
  res <- as.data.frame(t(sapply(c("DEL", "DUP", "INV", "CTX", "CPX"), function(svtype){
    idxs <- which(dat$SVTYPE==svtype 
                  & !(dat$chrom %in% c("X", "Y"))
                  & dat$AF<max.AF 
                  & (dat$SVLEN>=1000000 | dat$SVTYPE=="CTX"))
    variants <- length(idxs)
    carriers <- sum(dat$N_HET[idxs])
    # One CCR is actually a complex translocation,  and should be counted
    if(svtype=="CPX"){carriers <- carriers+1; variants <- variants+1}
    denom <- max(dat$N_BI_GENOS[idxs], na.rm=T)
    rate <- carriers/denom
    binom.ci <- binom.test(carriers, denom)$conf.int
    c(variants, carriers, rate, binom.ci)
  })))
  colnames(res) <- c("variants", "carriers", "mean", "lower", "upper")
  res <- rbind(apply(res, 2, sum), res)
  rownames(res)[1] <- "ALL"
  return(res)
}



######################
### PLOTTING FUNCTIONS
######################
# Plot simple log-scaled bars of SV by count
plot.totalCounts <- function(dat, svtypes, thousandG=F, ymax=NULL){
  # Gather data
  counts <- lapply(svtypes$svtype[which(svtypes$svtype!="OTH")], function(svtype){
    return(log10(length(which(dat$SVTYPE==svtype))))
  })
  names(counts) <- svtypes$svtype[which(svtypes$svtype!="OTH")]
  # counts$OTH <- log10(length(which(!(dat$SVTYPE %in% svtypes$svtype[which(svtypes$svtype!="OTH")]))))
  counts <- lapply(counts, function(val){if(is.infinite(val)){val <- 0}; return(val)})
  counts <- unlist(counts)
  sudmant.counts <- log10(c(42279, 6025, 2929, 16631+168, 786, NA, NA))
  names(sudmant.counts) <- names(counts)
  
  # Plot
  log.minor <- log10(as.numeric(sapply(0:10, function(i){(1:8)*(10^i)})))
  log.major <- 1:9
  if(is.null(ymax)){
    ymax <- max(counts, na.rm=T)
  }
  par(mar=c(2.5, 3, 2, 0.5), bty="n")
  plot(x=c(0, 2*nrow(svtypes)), y=c(0, 1.02*ymax), type="n", 
       xaxt="n", xlab="", yaxt="n", ylab="", yaxs="i")
  if(thousandG==T){
    rect(xleft=seq(1, 2*length(counts), 2)-0.8, xright=seq(1, 2*length(counts), 2), 
         ybottom=0, ytop=sudmant.counts, 
         lwd=0.5, col="gray70")
    rect(xleft=seq(2, 2*length(counts), 2)-1, xright=seq(2, 2*length(counts), 2)-0.2, 
         ybottom=0, ytop=counts, 
         lwd=0.5, col=svtypes$color)
  }else{
    rect(xleft=seq(1, 2*length(counts), 2)-0.8, xright=seq(2, 2*length(counts), 2)-0.2, 
         ybottom=0, ytop=counts, 
         lwd=0.5, col=svtypes$color)
  }
  segments(x0=seq(1, 2*length(counts), 2)-0.8, 
           x1=seq(2, 2*length(counts), 2)-0.2, 
           y0=par("usr")[3], y1=par("usr")[3], lwd=2)
  axis(1, at=seq(1, 2*length(counts), 2), tick=F, line=-0.9, 
       labels=svtypes$svtype[which(svtypes$svtype!="OTH")], cex.axis=0.8, las=2)
  axis(2, at=log.minor, labels=NA, tck=-0.02, lwd=0.9)
  axis(2, at=log.major, labels=NA, tck=-0.04, lwd=1.1)
  axis(2, at=1:6, tick=F, line=-0.5, cex.axis=0.8, las=2, 
       labels=c("10", "100", "1k", "10k", "100k", "1M"))
  mtext(2, text="SVs Discovered", line=2)
  sapply(1:length(counts), function(i){
    if(thousandG==T){
      if(is.na(sudmant.counts[i])){
        sudmant.y <- 0
        sudmant.lab <- 0
      }else{
        sudmant.y <- sudmant.counts[i]
        sudmant.lab <- 10^sudmant.counts[i]
      }
      text(x=2*i-1.4, y=sudmant.y+(0.025*(par("usr")[4]-par("usr")[3])), 
           labels=prettyNum(sudmant.lab, big.mark=","), 
           srt=90, adj=0, xpd=T, cex=0.7, col="gray60")
      text(x=2*i-0.55, y=counts[i]+(0.025*(par("usr")[4]-par("usr")[3])), 
           labels=prettyNum(10^counts[i], big.mark=","), 
           srt=90, adj=0, xpd=T, cex=0.7)
    }else{
      text(x=2*i-1, y=counts[i]+(0.025*(par("usr")[4]-par("usr")[3])), 
           labels=prettyNum(10^counts[i], big.mark=","), 
           srt=90, adj=0, xpd=T, cex=0.7)
    }
  })
}

# Plot rolling mean size distribution lines
plot.sizes <- function(dat, svtypes, step=0.02, xlims=c(50, 10000000), ymax=NULL){
  # Iterate over SVTYPES and get log-scaled count of SV by class per bin
  size.bins <- seq(log10(50), log10(50000000), step)
  plot.dat <- as.data.frame(sapply(svtypes$svtype[which(svtypes$svtype!="OTH")], function(svtype){
    res <- sapply(1:length(size.bins), function(i){
      nrow(dat[which(dat$SVLEN>=10^size.bins[i] & dat$SVLEN<10^(size.bins[i]+step) & dat$SVTYPE==svtype), ])
    })
    return(res)
  }))
  plot.dat <- apply(plot.dat, 2, function(vals){rollapply(vals, sum, width=5, partial=T)})
  plot.dat <- apply(plot.dat, 2, log10)
  plot.dat[which(is.infinite(plot.dat))] <- 0
  plot.dat <- as.data.frame(plot.dat)
  
  # Prep plot area
  par(bty="n", mar=c(2, 3.5, 0.5, 1))
  L1.peak <- log10(6000)
  SVA.peak <- log10(1200)
  ALU.peak <- log10(280)
  log.minor <- log10(as.numeric(sapply(0:10, function(i){(1:8)*(10^i)})))
  log.mid <- log10(as.numeric(sapply(0:10, function(i){c(1, 5)*(10^i)})))
  log.major <- 1:9
  if(is.null(ymax)){
    ymax <- max(plot.dat, na.rm=T)
  }
  plot(x=log10(xlims), y=c(0, 1.05*ymax), type="n", 
       xaxt="n", xaxs="i", xlab="", yaxt="n", ylab="", yaxs="i")

  # Plot lines per SV class
  sapply(1:ncol(plot.dat), function(i){
    points(x=size.bins, y=plot.dat[, i], lwd=2.5, col=svtypes$color[i], type="l")
  })
  
  # Add axes
  axis(1, at=log.minor, labels=NA, tck=-0.0175, lwd=0.9)
  axis(1, at=log.major, labels=NA, tck=-0.035, lwd=1.1)
  sapply(1:6, function(i){
    axis(1, at=i+1, tick=F, line=-0.9, cex.axis=0.65, 
         labels=c("100bp", "1kb", "10kb", "100kb", "1Mb", "10Mb")[i])
  })
  
  mtext(1, text="SV Size", line=1)
  axis(2, at=log.minor, labels=NA, tck=-0.0175, lwd=0.9)
  axis(2, at=log.major, labels=NA, tck=-0.035, lwd=1.1)
  axis(2, at=1:6, tick=F, line=-0.6, cex.axis=0.8, las=2, 
       labels=c("10", "100", "1k", "10k", "100k", "1M"))
  mtext(2, text="SV Discovered", line=2)
}

#Plot distribution of allele counts by SVTYPE
plot.freqByType <- function(dat, ymin=NULL, axlabel.cex=1){
  # Gather data
  exclude <- grep("MULTIALLELIC", dat$FILTER, fixed=T)
  if(length(exclude)>0){
    dat <- dat[-exclude, ]
  }
  dat <- dat[which(dat$SVTYPE!="MCNV" & !(dat$chrom %in% c("X", "Y"))), ]
  AC.cutoffs <- c(1:9, seq(10, 90, 10), seq(100, 900, 100), seq(1000, 25000, 1000))
  AF.dat <- as.data.frame(sapply(c("INS", "DUP", "DEL", "INV", "CPX"), function(svtype){
    sapply(AC.cutoffs, function(AC){
      length(which(dat$AC[which(dat$SVTYPE==svtype)]<=AC))/length(which(dat$SVTYPE==svtype))
    })
  }))
  
  # Plot
  if(is.null(ymin)){
    ymin <- log10(floor(100*min(AF.dat, na.rm=T))/100)
  }else{
    ymin <- log10(ymin)
  }
  AF.dat <- as.data.frame(apply(AF.dat, 2, log10))
  xrange <- c(-0.25, max(log10(AC.cutoffs)))
  common.threshold <- min(as.numeric(dat$AC[which(as.numeric(dat$AF)>=0.01 & as.numeric(dat$AN)==max(dat$AN, na.rm=T))]), na.rm=T)
  par(mar=c(2, 3.25, 0.5, 0.5), bty="n")
  plot(x=xrange, y=c(ymin, log10(1)), type="n", 
       xaxt="n", xlab="", yaxt="n", ylab="")
  segments(x0=log10(common.threshold), x1=log10(common.threshold), 
           y0=par("usr")[3], y1=log10(1), col="gray80")
  text(x=log10(common.threshold)+(0.025*(par("usr")[2]-par("usr")[1])), 
       y=par("usr")[3]+(0.08*(par("usr")[4]-par("usr")[3])), 
       labels="Rare\n(AF<1%)", pos=2, cex=0.7)
  text(x=log10(common.threshold)-(0.025*(par("usr")[2]-par("usr")[1])), 
       y=par("usr")[3]+(0.08*(par("usr")[4]-par("usr")[3])), 
       labels="Common\n(AF>1%)", pos=4, cex=0.7)
  sapply(c("INS", "DUP", "DEL", "INV", "CPX"), function(svtype){
    points(x=log10(AC.cutoffs), y=AF.dat[, which(colnames(AF.dat)==svtype)], 
           col=svtypes$color[which(svtypes$svtype==svtype)], type="l", lwd=3)
  })
  sapply(c("INS", "DUP", "DEL", "INV", "CPX"), function(svtype){
    rect(xleft=log10(AC.cutoffs)[1]-(0.03*(par("usr")[2]-par("usr")[1])), 
         xright=log10(AC.cutoffs)[1]+(0.03*(par("usr")[2]-par("usr")[1])), 
         ybottom=AF.dat[, which(colnames(AF.dat)==svtype)][1]-(0.01*(par("usr")[4]-par("usr")[3])), 
         ytop=AF.dat[, which(colnames(AF.dat)==svtype)][1]+(0.01*(par("usr")[4]-par("usr")[3])), 
         col=svtypes$color[which(svtypes$svtype==svtype)])
  })
  logscale.all <- log10(as.numeric(unlist(sapply(c(0:9), function(i){(1:9)*(10^i)}))))
  logscale.major <- 0:9
  axis(1, at=logscale.all, labels=NA, tck=-0.015, lwd=0.7)
  axis(1, at=logscale.major, labels=NA, tck=-0.03, lwd=1.1)
  sapply(1:5, function(i){
    axis(1, at=i-1, tick=F, line=-0.9, cex.axis=0.8, 
         labels=c("1", "10", "100", "1k", "10k")[i])
  })
  mtext(1, text="Allele Count", line=1, cex=axlabel.cex)
  logscale.pct.all <- log10((1:100)/100)
  logscale.pct.major <- log10(seq(10, 100, 10)/100)
  axis(2, at=logscale.pct.all, labels=NA, tck=-0.015, lwd=0.9)
  axis(2, at=logscale.pct.major, labels=NA, tck=-0.03, lwd=1.1)
  axis(2, at=logscale.pct.major, tick=F, line=-0.5, cex.axis=0.8, las=2, 
       labels=paste(seq(10, 100, 10), "%", sep=""))
  mtext(2, text="Fraction of SVs", line=2.25, cex=axlabel.cex)
}

# Hardy-Weinberg ternary plot
plot.HWE <- function(dat, pop=NULL, title=NULL, full.legend=F, lab.cex=1){
  # Gather HW p-values & colors
  HWE.mat <- makeHWEmat(dat=dat, pop=pop)
  HW.p <- HWChisqStats(X=HWE.mat, x.linked=F, pvalues=T)
  HW.cols <- rep("#4DAC26", times=length(HW.p))
  HW.cols[which(HW.p<0.05)] <- "#81F850"
  HW.cols[which(HW.p<0.05/length(HW.p))] <- "#AC26A1"
  
  # Generate HW plot frame
  par(mar=c(1, 1, 1, 1), bty="n")
  plot(x=1.15*c(-1/sqrt(3), 1/sqrt(3)), y=c(-0.15, 1.15), type="n", 
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  segments(x0=c(-1/sqrt(3), 0, 1/sqrt(3)), 
           x1=c(0, 1/sqrt(3), -1/sqrt(3)), 
           y0=c(0, 1, 0), y1=c(1, 0, 0))
  HWTernaryPlot(X=HWE.mat, n=max(HWE.mat, na.rm=T), newframe=F, 
                vbounds=F, mafbounds=F, 
                region=1, vertexlab=NA, 
                alpha=0.05, 
                curvecols=c("#4DAC26", "#81F850", NA, NA), pch=NA)
  
  # Add axes
  text(x=c(-1/sqrt(3), 1/sqrt(3)), y=0, labels=c("0/0", "1/1"), 
       pos=1, cex=0.8, xpd=T, font=2)
  text(x=0, y=1, labels="0/1", pos=3, cex=0.8, xpd=T, font=2)
  
  # Finish HW plot
  HWTernaryPlot(X=HWE.mat, n=max(HWE.mat, na.rm=T), newframe=F, 
                vbounds=F, mafbounds=F, 
                region=1, vertexlab=NA, 
                alpha=0.03/nrow(HWE.mat), 
                curvecols=c("#4DAC26", "#AC26A1", NA, NA), 
                pch=21, cex=0.3, signifcolour=F, markercol=NA, 
                markerbgcol=adjustcolor(HW.cols, alpha=0.25))
  segments(x0=c(-1/sqrt(3), 0, 1/sqrt(3)), 
           x1=c(0, 1/sqrt(3), -1/sqrt(3)), 
           y0=c(0, 1, 0), y1=c(1, 0, 0))
  
  # Add legend
  n.pass <- length(which(HW.p>=0.05))
  print(paste("PASS: ", n.pass/length(HW.p), sep=""))
  n.nom <- length(which(HW.p<0.05 & HW.p>=0.05/nrow(HWE.mat)))
  print(paste("NOMINAL FAILS: ", n.nom/length(HW.p), sep=""))
  n.bonf <- length(which(HW.p<0.05/nrow(HWE.mat)))
  print(paste("BONFERRONI FAILS: ", n.bonf/length(HW.p), sep=""))
  legend("topright", pch=19, col=c("#4DAC26", "#81F850", "#AC26A1"), pt.cex=1.3, 
         legend=c(paste(round(100*(n.pass/nrow(HWE.mat)), 0), "%", sep=""), 
                  paste(round(100*(n.nom/nrow(HWE.mat)), 0), "%", sep=""), 
                  paste(round(100*(n.bonf/nrow(HWE.mat)), 0), "%", sep="")), 
         bty="n", bg=NA, cex=0.7)
}

# Plot mutation rates
plotMus <- function(dat, Ne, werling=T, legend=T, ymax=NULL){
  mu.dat <- getAllMus(dat, Ne)
  plot.dat <- mu.dat$mu.Means
  mu.svtypes <- c("DEL", "DUP", "INS", "INV", "CPX")
  mu.cols <- c("gray30", sapply(mu.svtypes, function(svtype){
    svtypes$color[which(svtypes$svtype==svtype)]
  }))
  werling.rates <- c(166, 80, 47, 37, 0, 1)/1038
  
  # Prep plot area
  if(is.null(ymax)){
    ymax <- max(plot.dat, na.rm=T)
  }
  par(mar=c(1.2, 3, 0.5, 0.5), bty="n")
  plot(x=c(-0.025, ncol(plot.dat))+0.2, y=c(0, ymax), type="n", 
       xaxt="n", yaxt="n", xlab="", ylab="")
  axis(1, at=c(-4, 10), tck=0, labels=NA, lwd=1.5)
  axis(2, at=c(-1, 1), tck=0, labels=NA, lwd=1.5)
  # axis(1, at=(1:ncol(plot.dat))-0.5, tick=F, line=-0.8, 
  #      labels=c("All SV", mu.svtypes))
  mtext(1, line=0.1, text="SV Class")
  axis(2, labels=NA, tck=-0.025)
  axis(2, tick=F, las=2, cex.axis=0.8, line=-0.6)
  mtext(2, line=2, text="Mutation Rate")
  
  # Add points & CIs
  text.buf <- 0.04*(par("usr")[4]-par("usr")[3])
  sapply(1:ncol(plot.dat), function(i){
    segments(x0=i-0.5, x1=i-0.5, 
             y0=plot.dat[1, i], y1=plot.dat[3, i], 
             lend="round", lwd=2, col=mu.cols[i])
    points(x=i-0.5, y=plot.dat[2, i], pch=19, cex=1.1, col=mu.cols[i])
    par(xpd=T)
    if(i==3){
      text(x=i-0.65, y=plot.dat[3, i], srt=55, 
           cex=0.7, labels=format(round(plot.dat[2, i], 3), nsmall=3), pos=4, 
           col=mu.cols[i])
    }else{
      text(x=i-0.45, y=plot.dat[2, i], srt=55, 
           cex=0.7, labels=format(round(plot.dat[2, i], 3), nsmall=3), pos=4, 
           col=mu.cols[i])
    }
    par(xpd=F)
    if(werling==T){
      points(x=i-0.5, y=werling.rates[i], pch=23, lwd=1.5)
    }
  })
  if(legend==T){
    if(werling==T){
      legend("topright", border=NA, bty="n", 
             legend=c(as.expression(bquote(mu ~ "from Watterson" ~ hat(theta[italic(W)]) ~ "in gnomAD,  " ~ N[e] == .(prettyNum(Ne, big.mark=",")))), 
                      expression("Rate of validated" ~ italic("de novo") ~ "SV from 519 quartets")), 
             pch=c(19, 23), pt.lwd=1.5, lwd=c(2, NA), lty=c(1, NA), cex=0.8, pt.cex=c(1.25, 1))
    }else{
      legend("topright", border=NA, bty="n", 
             legend=as.expression(bquote(mu ~ "from Watterson" ~ hat(theta[italic(W)]) ~ "in gnomAD,  " ~ N[e] == .(prettyNum(Ne, big.mark=",")))), 
             pch=19, pt.lwd=1.5, lwd=2, lty=1, cex=0.8, pt.cex=1.25)
    }
  }
}

# Make summed dotplots for constraint comparisons
plotSummedDots <- function(deciles, SV_count_obs, SV_count_exp, 
                           color, title=NULL, ymax=NULL, 
                           xlabel="SNV pLoF Constraint Percentile", 
                           ax.labels=TRUE, cex.labels=0.8, 
                           tck=NULL, yline=-0.4, conf.int=F, parmar=c(2, 3.25, 1.75, 2)){
  require(zoo, quietly=T)
  d <- sort(unique(as.numeric(deciles)))
  # Compute summed obs/exp
  means <- sapply(d, function(i){
    exp.sum <- sum(SV_count_exp[which(deciles==i)])
    obs.sum <- sum(SV_count_obs[which(deciles==i)])
    return(obs.sum/exp.sum)
  })
  means[which(means<0)] <- NA
  # Compute CI
  if(conf.int==T){
    cis <- t(sapply(d, function(i){
      oe.vals <- SV_count_obs[which(deciles==i)]/SV_count_exp[which(deciles==i)]
      get.mean <- function(vals, indices){mean(vals[indices])}
      set.seed(i)
      boot.obj <- boot(data=oe.vals, statistic=get.mean, R=1000)
      ci <- boot.ci(boot.obj, conf=0.9, type="basic")$basic[4:5]
      return(ci)
    }))
  }
  
  # Prep plot area
  if(is.null(ymax)){
    ymax <- max(means, na.rm=T)
  }
  par(mar=parmar)
  plot(x=c(0, length(d)), y=c(0, ymax), type="n", 
       xaxt="n", yaxt="n", xlab="", ylab="")
  sapply(seq(0, 100, 20), function(p){
    axis(1, at=p, tick=F, cex.axis=0.85*cex.labels, line=-0.9, 
         labels=bquote(.(p)^'th'))
  })
  axis(2, at=axTicks(2), labels=NA, tck=tck)
  if(ax.labels==T){
    mtext(1, line=0.9, text=xlabel, cex=cex.labels)
    mtext(2, line=1.75, text="Rare SV Obs/Exp", cex=cex.labels)
  }
  axis(2, at=axTicks(2), line=yline, cex.axis=0.85*cex.labels, 
       labels=paste(100*round(axTicks(2), 2), "%", sep=""), tick=F, las=2)
  mtext(3, text=title, line=0.2, cex=0.9)
  abline(h=1, lwd=2, col="gray80")
  
  # Add points & rolling mean/ci
  if(conf.int==T){
    polygon(x=c(d-0.5, rev(d-0.5)), 
            y=c(rollapply(cis[, 1], 21, mean, na.rm=T, partial=T), 
                rev(rollapply(cis[, 2], 21, mean, na.rm=T, partial=T))), 
            border=NA, bty="n", col=adjustcolor(color, alpha=0.1))
  }
  points(x=d-0.5, y=means, pch=21, col=color, cex=0.4)
  points(x=d-0.5, y=rollapply(means, 21, mean, na.rm=T, partial=T), 
         lwd=2, type="l", col=color)
  
  #Add correlation coefficient
  cor.res <- cor.test(x=d-0.5, y=means, method="spearman")
  text(x=par("usr")[1]-(0.025*(par("usr")[2]-par("usr")[1])), 
       y=par("usr")[3]+(0.885*(par("usr")[4]-par("usr")[3])), 
       pos=4, cex=0.8*cex.labels, 
       labels=bquote(rho == .(format(round(cor.res$estimate, 2), nsmall=2))))
  cor.p <- cor.res$p.value
  if(cor.p>10^-100){
    cor.p <- format(cor.p, scientific=T)
    cor.p.base <- format(as.numeric(strsplit(cor.p, split="e")[[1]][1]), nsmall=2)
    cor.p.exp <- format(as.numeric(strsplit(cor.p, split="e")[[1]][2]), nsmall=0)
    text(x=par("usr")[2]+(0.025*(par("usr")[2]-par("usr")[1])), 
         y=par("usr")[3]+(0.075*(par("usr")[4]-par("usr")[3])), 
         pos=2, cex=0.8*cex.labels, 
         labels=bquote(italic(P) == .(format(round(as.numeric(cor.p.base), 2), nsmall=2))*"x"*10^.(cor.p.exp)))
  }else{
    text(x=par("usr")[2]+(0.025*(par("usr")[2]-par("usr")[1])), 
         y=par("usr")[3]+(0.075*(par("usr")[4]-par("usr")[3])), 
         pos=2, cex=0.8*cex.labels, 
         labels=bquote(italic(P) < 10^-100))
  }
  box()
}

# Wrapper for all constraint analyses for a single class of SV
constraintModelWrapper <- function(gene.data, y, snv.class="ptv", color, title, 
                                   ymax=NULL, return=FALSE, ax.labels=TRUE, 
                                   conf.int=F, cex.labels=0.8, 
                                   tck=NULL, yline=-0.4, parmar=c(2, 3.25, 1.75, 2)){
  if(snv.class=="ptv"){
    train.genes <- which(gene.data$ptv_oe_dec>=5 & gene.data$ptv_oe_dec<=10)
    deciles <- gene.data$ptv_oe_cent
    xlabel <- "SNV pLoF Constraint Pct."
  }else{
    train.genes <- which(gene.data$mis_oe_dec>=5 & gene.data$mis_oe_dec<=10)
    deciles <- gene.data$mis_oe_cent
    xlabel <- "Missense Constraint Pct."
  }
  cov <- prepCovariates(gene.data=gene.data, y=y)
  res <- fitConstraintModel(cov=cov, train.genes=train.genes)
  print(paste("VARIANCE EXPLAINED: ", 100*(cor(res[, 2], res[, 3], use="complete.obs")^2), "%", sep=""))
  plotSummedDots(deciles=deciles, 
                 SV_count_obs=res$SV_count_raw, 
                 SV_count_exp=res$SV_count_glm_exp, 
                 color=color, title=title, ymax, 
                 xlabel=xlabel, 
                 ax.labels=ax.labels, 
                 cex.labels=cex.labels, 
                 conf.int=conf.int, 
                 parmar=parmar, 
                 tck=tck, yline=yline)
  if(return==T){
    return(res)
  }
}

# Plot vertical barplot for rare pLoF carrier rates for a list of gene lists
plot.rareLoFCarrierByList.vertical <- function(pop="ALL", genelists.list, dat, modes, genelist.labels=NULL, 
                                               ymax=NULL, y.break=NULL, ylabel=NULL, ax.lab.cex=0.75,
                                               parmar=c(5.75, 3.75, 1, 0.5)){
  res <- do.call("rbind", lapply(1:length(genelists.list), function(i){
    countRareLoFCarrierRate(genelist=genelists.list[[i]], pop=pop, dat=dat, mode=modes[i])
  }))
  # Prep plot area
  plot.colors <- c(svtypes$color[which(svtypes$svtype=="DEL")], 
                   svtypes$color[which(svtypes$svtype=="DUP")], 
                   svtypes$color[which(svtypes$svtype=="INS")], 
                   svtypes$color[which(svtypes$svtype=="INV")], 
                   svtypes$color[which(svtypes$svtype=="CPX")])
  par(bty="n")
  if(is.null(ymax)){
    ymax <- max(apply(res, 1, sum))
  }
  if(!is.null(y.break)){
    ymax <- ymax-(y.break[2]-y.break[1])
  }
  par(mar=parmar, bty="n", xpd=F)
  plot(x=c(0, nrow(res))+1, y=c(0, ymax), type="n", 
       xaxt="n", xlab="", xaxs="i", yaxt="n", ylab="", yaxs="i")
  par(xpd=T)
  rect(xleft=1, xright=par("usr")[2], 
       ybottom=par("usr")[3]-(0.06*(par("usr")[4]-par("usr")[3])), 
       ytop=par("usr")[3]-(0.01*(par("usr")[4]-par("usr")[3])), 
       bty="n", border=NA, col="gray90")
  par(xpd=F)
  # Plot bars
  sapply(1:nrow(res), function(i){
    rect(ybottom=c(0, cumsum(res[i, ])[-ncol(res)]), 
         ytop=cumsum(res[i, ]), 
         xleft=i+0.2, xright=i+0.8, 
         bty="n", border=NA, col=plot.colors)
    rect(ybottom=0, ytop=sum(res[i, ]), 
         xleft=i+0.2, xright=i+0.8, 
         col=NA)
    if(sum(res[i, ])<max(y.break)){
      text(y=sum(res[i, ])-(0.025*(par("usr")[4]-par("usr")[3])), 
           x=i+0.5, cex=0.75*ax.lab.cex, xpd=T,
           pos=3, labels=paste(round(100*sum(res[i, ]), 2), "%", sep=""))
    }
    if(!is.null(genelist.labels)){
      par(xpd=T)
      text(x=i+0.5, y=par("usr")[3]-(0.035*(par("usr")[4]-par("usr")[3])), 
           labels=prettyNum(length(genelists.list[[i]]), big.mark=","), 
           cex=0.8*ax.lab.cex)
      text(x=i+1, y=par("usr")[3]-(0.09*(par("usr")[4]-par("usr")[3])), 
           labels=genelist.labels[i], srt=40, pos=2, cex=0.9*ax.lab.cex)
      par(xpd=F)
    }
  })
  # Clean up bottom panel
  axis(2, at=c(par("usr")[3], ymax), labels=NA, tck=0)
  if(!is.null(y.break)){
    axis(2, at=axTicks(2)[which(axTicks(2)<y.break[1])], labels=NA, tck=-0.025)
    axis(2, at=axTicks(2)[which(axTicks(2)<y.break[1])], tick=F, line=-0.6, cex.axis=0.9, las=2, 
         labels=paste(round(100*axTicks(2)[which(axTicks(2)<y.break[1])], 1), "%", sep=""))
  }
  mtext(2, line=1.9, text=ylabel, cex=1.2*ax.lab.cex)
  # Add top panel, if necessary
  if(!is.null(y.break)){
    # Clear all data above future axis break
    rect(xleft=par("usr")[1]+0.2, xright=par("usr")[2]-0.2, 
         ybottom=y.break[1], 
         ytop=par("usr")[4], 
         col="white", border="white", lwd=4)
    # Add stacked rectangles as required
    sapply(which(apply(res, 1, sum)>y.break[1]), function(i){
      revised.vals <- cumsum(res[i, ])-(y.break[2]-y.break[1])
      rect(xleft=i+0.2, xright=i+0.8, 
           ybottom=c(y.break[1], revised.vals[-length(revised.vals)]), 
           ytop=revised.vals, 
           bty="n", border=NA, col=plot.colors)
      rect(xleft=i+0.2, xright=i+0.8, 
           ybottom=0, ytop=max(revised.vals), 
           col=NA)
      # Add bar labels
      text(y=max(revised.vals)-(0.025*(par("usr")[4]-par("usr")[3])), 
           x=i+0.5, cex=0.75*ax.lab.cex, xpd=T,
           pos=3, labels=paste(round(100*sum(res[i, ]), 2), "%", sep=""))
    })
    # Add top axis
    axis(2, at=axTicks(2)[which(axTicks(2)>y.break[1])], labels=NA, tck=-0.025)
    axis(2, at=axTicks(2)[which(axTicks(2)>y.break[1])], tick=F, line=-0.6, cex.axis=0.9, las=2, 
         labels=paste(round(100*(axTicks(2)[which(axTicks(2)>y.break[1])]+(y.break[2]-y.break[1])), 1), "%", sep=""))
    # Add axis break
    rect(xleft=par("usr")[1], xright=par("usr")[2], 
         ybottom=y.break[1]-(0.01*(par("usr")[4]-par("usr")[3])), 
         ytop=y.break[1]+(0.01*(par("usr")[4]-par("usr")[3])), 
         col="white", bty="n", border=NA)
    axis.break(2, breakpos=y.break[1], brw=0.06)
  }
}

#Plot carrier rate for gross chromosomal abnormalities
plot.chrom.abnormalities <- function(dat, max.AF=0.01, xmax=NULL,
                                     category.labels=c("DEL", 
                                                       "DUP", 
                                                       "INV", 
                                                       "CTX", 
                                                       "CPX"),
                                     cex.toplabels=0.9){
  # Get plotting data
  plot.dat <- gather.bca.table(dat=dat, max.AF=max.AF)
  plot.dat <- plot.dat[-1, ]
  if(is.null(xmax)){
    xmax <- max(plot.dat[, -c(1:2)])
  }
  point.colors <- as.character(sapply(rownames(plot.dat), function(svtype){
    if(svtype=="CTX"){
      svtype <- "OTH"
    }
    svtypes$color[which(svtypes$svtype==svtype)]
  }))
  
  # Prep plot area
  par(mar=c(0.3, 2.5, 2, 0.3), bty="n", xpd=F)
  plot(x=c(-0.05*xmax, xmax), y=c(0, -nrow(plot.dat)), type="n", 
       xlab="", ylab="", xaxt="n", yaxt="n")
  abline(v=axTicks(3), col="gray90")
  abline(v=0, col="gray50")
  
  # Plot points & CIs
  segments(x0=plot.dat$lower, x1=plot.dat$upper, 
           y0=(-1:-nrow(plot.dat))+0.5, 
           y1=(-1:-nrow(plot.dat))+0.5, 
           lend="round", lwd=2)
  points(x=plot.dat$mean, 
         y=(-1:-nrow(plot.dat))+0.5, 
         pch=21, cex=1.1, bg=point.colors)
  par(xpd=T)
  pct.lab.pos <- sapply(plot.dat$mean, function(x){if(x<0.01){4}else{2}})
  pct.lab.coord <- sapply(1:nrow(plot.dat), function(i){if(plot.dat$mean[i]<0.01){plot.dat$upper[i]-0.001}else{plot.dat$lower[i]+0.001}})
  text(x=pct.lab.coord, y=(-1:-nrow(plot.dat))+0.45, cex=0.75, 
       pos=pct.lab.pos, labels=paste(format(round(100*plot.dat$mean, 2), nsmall=2), "%", sep=""))
  par(xpd=F)
  axis(2, at=(-1:-nrow(plot.dat))+0.5, line=-0.9, tick=F, labels=category.labels, las=2, cex.axis=0.9)
  axis(3, at=axTicks(3), labels=NA, tck=-0.025)
  sapply(axTicks(3)[seq(1, length(axTicks(3)), 2)], function(x){
    axis(3, at=x, tick=F, line=-0.8, cex.axis=0.9, 
         labels=paste(round(100*x, 1), "%", sep=""))
  })
  mtext(3, line=1.1, text="Samples (Pct.)", cex=cex.toplabels)
}


#################
### RSCRIPT BLOCK
#################
require(optparse, quietly=T)
require(zoo, quietly=T)
require(HardyWeinberg, quietly=T)
require(MASS, quietly=T)
require(plotrix, quietly=T)
### List of command-line options
option_list <- list(
  make_option(c("--prefix"),  type="character",  default="gnomAD_v2_SV", 
              help="prefix used for naming outfiles [default %default]", 
              metavar="character"), 
  make_option(c("-S",  "--svtypes"),  type="character",  default=NULL, 
              help="tab-delimited file specifying SV types and HEX colors [default %default]", 
              metavar="character"), 
  make_option(c("-P",  "--populations"),  type="character",  default=NULL, 
              help="tab-delimited file specifying populations and HEX colors [default %default]", 
              metavar="character")
)


### Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog VCF2BED SD_SR_COV QUALS APS_STATS SNV_DATA GENE_METADATA GENELIST_DIR OUTDIR", 
                                option_list=option_list), 
                   positional_arguments=TRUE)
opts <- args$options


### Checks for appropriate positional arguments
if(length(args$args) != 8){
  stop("Incorrect number of required positional arguments\n")
}


### Writes args & opts to vars
vcf2bed.in <- args$args[1]
sd_sr_cov.in <- args$args[2]
quals.in <- args$args[3]
aps.in <- args$args[4]
SNVdata.in <- args$args[5]
gene_metadata.in <- args$args[6]
genelist.dir <- args$args[7]
OUTDIR <- args$args[8]
prefix <- opts$prefix
svtypes.file <- opts$svtypes
pops.file <- opts$populations


# # Dev parameters (local)
# vcf2bed.in <- "~/scratch/gnomAD-SV_v2_rev1.vcf2bed.bed.gz"
# sd_sr_cov.in <- "~/scratch/gnomAD-SV_v2_rev1.variant_SD_SR_coverage.txt.gz"
# quals.in <- "~/scratch/gnomAD-SV_v2_rev1.QUAL_per_SV.txt.gz"
# aps.in <- "~/scratch/gnomAD-SV_v2_rev1.APS_stats.txt"
# SNVdata.in <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/gnomAD_v2.1_canonical_constraint.condensed.txt.gz"
# gene_metadata.in <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/Gencode.v19.autosomal_canonical_gene_metadata.txt.gz"
# genelist.dir <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/genelists/"
# OUTDIR <- "~/scratch/gnomAD_v2_test_plots/"
# prefix <- "gnomAD-SV_v2_rev1"
# svtypes.file <- "~/Desktop/Collins/Talkowski/code/sv-pipeline/ref/vcf_qc_refs/SV_colors.txt"
# pops.file <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/gnomAD_population_colors.txt"


### Create output directory,  if necessary
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}


### Process input data
cat("NOW LOADING AND CLEANING DATA\n")
dat.lq <- read.vcf2bed(vcf2bed.in, sd_sr_cov.in, quals.in, aps.in)
dat.lq <- dat.lq[which(dat.lq$SVTYPE!="BND"), ]
dat.mq <- dat.lq[which(dat.lq$FILTER=="PASS" | dat.lq$FILTER=="MULTIALLELIC"), ]
dat.hq <- dat.mq[which(dat.mq$QUAL>500), ]
snv.data <- read.snvdata(SNVdata.in, gene_metadata.in)
genes <- sort(unique(as.character(snv.data$gene)))
gene.data.lq <- getSVdat.all(dat=dat.lq, snv.data=snv.data)
gene.data.mq <- getSVdat.all(dat=dat.mq, snv.data=snv.data)
gene.data.hq <- getSVdat.all(dat=dat.hq, snv.data=snv.data)


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


### Sets populations & colors
if(!is.null(pops.file)){
  pops <- read.table(pops.file, sep="\t", header=T, comment.char="", check.names=F)
}


### Import selected gene lists
ACMG.genes <- importGenelist(genelist.dir, "ACMG_59.genes.list", gene.data.lq)
DDD.genes <- importGenelist(genelist.dir, "DDD_2017.genes.list", gene.data.lq)
clingenHC.genes <- importGenelist(genelist.dir, 
                                  "ClinGen_haploinsufficient_high_confidence.genes.list", 
                                  gene.data.lq)
dominant_DD.genes <- importGenelist(genelist.dir, 
                                    "DDG2P_confirmed_LoF_dominant_DD.genes.list", 
                                    gene.data.lq)
recessive_DD.genes <- importGenelist(genelist.dir, 
                                     "DDG2P_confirmed_LoF_recessive_DD.genes.list", 
                                     gene.data.lq)



### Plot simple bars of counts
pdf(paste(OUTDIR, "/", prefix, ".site_counts_by_type.lq_callset.pdf", sep=""), 
    height=2.25, width=2.25)
plot.totalCounts(dat=dat.lq, svtypes=svtypes, thousandG=F, ymax=log10(200000))
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".site_counts_by_type.mq_callset.pdf", sep=""), 
    height=2.25, width=2.25)
plot.totalCounts(dat=dat.mq, svtypes=svtypes, thousandG=F, ymax=log10(200000))
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".site_counts_by_type.hq_callset.pdf", sep=""), 
    height=2.25, width=2.25)
plot.totalCounts(dat=dat.hq, svtypes=svtypes, thousandG=F, ymax=log10(200000))
dev.off()


### Plot size distributions
pdf(paste(OUTDIR, "/", prefix, ".size_distribution_by_type.lq_callset.pdf", sep=""), 
    height=2, width=2.9)
plot.sizes(dat=dat.lq, svtypes=svtypes, ymax=5)
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".size_distribution_by_type.mq_callset.pdf", sep=""), 
    height=2, width=2.9)
plot.sizes(dat=dat.mq, svtypes=svtypes, ymax=5)
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".size_distribution_by_type.hq_callset.pdf", sep=""), 
    height=2, width=2.9)
plot.sizes(dat=dat.hq, svtypes=svtypes, ymax=5)
dev.off()


### Plot frequency distributions by SVTYPE
pdf(paste(OUTDIR, "/", prefix, ".site_frequency_distributions_bySVtype.lq_callset.pdf", sep=""), 
    height=2, width=2.3)
plot.freqByType(dat=dat.lq, ymin=0.40)
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".site_frequency_distributions_bySVtype.mq_callset.pdf", sep=""), 
    height=2, width=2.3)
plot.freqByType(dat=dat.mq, ymin=0.40)
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".site_frequency_distributions_bySVtype.hq_callset.pdf", sep=""), 
    height=2, width=2.3)
plot.freqByType(dat=dat.hq, ymin=0.40)
dev.off()


### Plot Hardy-Weinberg Equilibrium
png(paste(OUTDIR, "/", prefix, ".HWE_all_samples.lq_callset.png", sep=""), 
    height=800, width=800, res=400)
plot.HWE(dat=dat.lq, pop=NULL)
dev.off()
png(paste(OUTDIR, "/", prefix, ".HWE_all_samples.mq_callset.png", sep=""), 
    height=800, width=800, res=400)
plot.HWE(dat=dat.mq, pop=NULL)
dev.off()
png(paste(OUTDIR, "/", prefix, ".HWE_all_samples.hq_callset.png", sep=""), 
    height=800, width=800, res=400)
plot.HWE(dat=dat.hq, pop=NULL)
dev.off()


### Plot mutation rate estimates
pdf(paste(OUTDIR, "/", prefix, ".mutation_rate_estimates.lq_callset.pdf", sep=""), 
    height=2, width=2.5)
plotMus(dat.lq, Ne=10000, werling=F, legend=F, ymax=0.5)
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".mutation_rate_estimates.mq_callset.pdf", sep=""), 
    height=2, width=2.5)
plotMus(dat.mq, Ne=10000, werling=F, legend=F, ymax=0.5)
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".mutation_rate_estimates.hq_callset.pdf", sep=""), 
    height=2, width=2.5)
plotMus(dat.hq, Ne=10000, werling=F, legend=F, ymax=0.5)
dev.off()


### Plot constraint comparisons for pLoF and CG
parmar <- c(2.25, 3, 0.5, 0.75)
pdf(paste(OUTDIR, "/", prefix, ".snv_constraint_comparison.pLoF.lq_callset.pdf", sep=""),
    height=1.5, width=2)
constraintModelWrapper(gene.data=gene.data.lq,
                       y="rare.lof.any",
                       color=svtypes$color[which(svtypes$svtype=="DEL")],
                       title="", ymax=1.5, ax.labels=T, cex.labels=0.85,
                       tck=-0.02, yline=-0.85, parmar=parmar)
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".snv_constraint_comparison.pLoF.mq_callset.pdf", sep=""),
    height=1.5, width=2)
constraintModelWrapper(gene.data=gene.data.mq,
                       y="rare.lof.any",
                       color=svtypes$color[which(svtypes$svtype=="DEL")],
                       title="", ymax=1.5, ax.labels=T, cex.labels=0.85,
                       tck=-0.02, yline=-0.85, parmar=parmar)
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".snv_constraint_comparison.pLoF.hq_callset.pdf", sep=""),
    height=1.5, width=2)
constraintModelWrapper(gene.data=gene.data.hq,
                       y="rare.lof.any",
                       color=svtypes$color[which(svtypes$svtype=="DEL")],
                       title="", ymax=1.5, ax.labels=T, cex.labels=0.85,
                       tck=-0.02, yline=-0.85, parmar=parmar)
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".snv_constraint_comparison.CG.lq_callset.pdf", sep=""),
    height=1.5, width=2)
constraintModelWrapper(gene.data=gene.data.lq,
                       y="rare.cg",
                       color=svtypes$color[which(svtypes$svtype=="DUP")],
                       title="", ymax=1.5, ax.labels=T, cex.labels=0.85,
                       tck=-0.02, yline=-0.85, parmar=parmar)
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".snv_constraint_comparison.CG.mq_callset.pdf", sep=""),
    height=1.5, width=2)
constraintModelWrapper(gene.data=gene.data.mq,
                       y="rare.cg",
                       color=svtypes$color[which(svtypes$svtype=="DUP")],
                       title="", ymax=1.5, ax.labels=T, cex.labels=0.85,
                       tck=-0.02, yline=-0.85, parmar=parmar)
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".snv_constraint_comparison.CG.hq_callset.pdf", sep=""),
    height=1.5, width=2)
constraintModelWrapper(gene.data=gene.data.hq,
                       y="rare.cg",
                       color=svtypes$color[which(svtypes$svtype=="DUP")],
                       title="", ymax=1.5, ax.labels=T, cex.labels=0.85,
                       tck=-0.02, yline=-0.85, parmar=parmar)
dev.off()


### Plot rare pLoF carrier rates in selected gene lists
pdf(paste(OUTDIR, "/", prefix, ".rare_medical_lof_carrierRates.lq_callset.pdf", sep=""), 
    height=2.25, width=2)
plot.rareLoFCarrierByList.vertical(pop="ALL", dat=dat.lq, 
                                   genelists.list=list(ACMG.genes, dominant_DD.genes, 
                                                       clingenHC.genes, 
                                                       recessive_DD.genes), 
                                   modes=c(rep("ALL", 3), "HET"), 
                                   genelist.labels=c("Incidental Report.", 
                                                     "Dominant DD", 
                                                     "ClinGen Haploinsuff.", 
                                                     "Recessive DD"), 
                                   y.break=c(0.02, 0.05), ymax=0.09, 
                                   ylabel="Samples with pLoF SV",
                                   parmar=c(4.25, 3.5, 1, 0.5))
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".rare_medical_lof_carrierRates.mq_callset.pdf", sep=""), 
    height=2.25, width=2)
plot.rareLoFCarrierByList.vertical(pop="ALL", dat=dat.mq, 
                                   genelists.list=list(ACMG.genes, dominant_DD.genes, 
                                                       clingenHC.genes, 
                                                       recessive_DD.genes), 
                                   modes=c(rep("ALL", 3), "HET"), 
                                   genelist.labels=c("Incidental Report.", 
                                                     "Dominant DD", 
                                                     "ClinGen Haploinsuff.", 
                                                     "Recessive DD"), 
                                   y.break=c(0.02, 0.05), ymax=0.09, 
                                   ylabel="Samples with pLoF SV",
                                   parmar=c(4.25, 3.5, 1, 0.5))
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".rare_medical_lof_carrierRates.hq_callset.pdf", sep=""), 
    height=2.25, width=2)
plot.rareLoFCarrierByList.vertical(pop="ALL", dat=dat.hq, 
                                   genelists.list=list(ACMG.genes, dominant_DD.genes, 
                                                       clingenHC.genes, 
                                                       recessive_DD.genes), 
                                   modes=c(rep("ALL", 3), "HET"), 
                                   genelist.labels=c("Incidental Report.", 
                                                     "Dominant DD", 
                                                     "ClinGen Haploinsuff.", 
                                                     "Recessive DD"), 
                                   y.break=c(0.02, 0.05), ymax=0.09, 
                                   ylabel="Samples with pLoF SV",
                                   parmar=c(4.25, 3.5, 1, 0.5))
dev.off()



### Plot chromosomal abnormality carrier rate
pdf(paste(OUTDIR, "/", prefix, ".chrom_abnormalities.lq_callset.pdf", sep=""), 
    height=1.5, width=2)
plot.chrom.abnormalities(dat.lq, xmax=0.02)
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".chrom_abnormalities.mq_callset.pdf", sep=""), 
    height=1.5, width=2)
plot.chrom.abnormalities(dat.mq, xmax=0.02)
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".chrom_abnormalities.hq_callset.pdf", sep=""), 
    height=1.5, width=2)
plot.chrom.abnormalities(dat.hq, xmax=0.02)
dev.off()



