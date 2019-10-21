#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# SNV/SV LD analyses for gnomAD-SV v2 rev1


### Set master parameters
options(stringsAsFactors=F, scipen=1000)


####################
### HELPER FUNCTIONS
####################
# Import & merge all data
import.data <- function(vcf2bed.in, eur_ld.in, afr_ld.in, 
                        gwas_hits.in, sd_sr_coverage.in, quals.in, 
                        minAF=0.01, minSNV=2000, exclude.allosomes=T){
  # Read data
  dat <- read.table(vcf2bed.in, comment.char="", header=T, sep="\t")
  colnames(dat)[1] <- "chrom"
  # Restrict to biallelic PASS sites
  dat <- dat[which(dat$FILTER=="PASS"), ]
  # Restrict to autosomes, if optioned
  if(exclude.allosomes == T){
    dat <- dat[which(dat$chrom %in% c(1:22)), ]
  }
  # Drop columns not being used (to save memory)
  cols.to.drop <- c("svtype", "CHR2", "CPX_INTERVALS", 
                    "CPX_TYPE", "END2", "END", "POS2",
                    "SOURCE", "STRANDS", "UNRESOLVED_TYPE",
                    colnames(dat)[grep("AMR_", colnames(dat), fixed=T)],
                    colnames(dat)[grep("EAS_", colnames(dat), fixed=T)],
                    colnames(dat)[grep("OTH_", colnames(dat), fixed=T)])
  dat <- dat[, -which(colnames(dat) %in% cols.to.drop)]
  # Convert numeric columns
  numeric.columns <- sort(unique(c(grep("FREQ", colnames(dat), fixed=T), 
                                   grep("AN", colnames(dat), fixed=T), 
                                   grep("AC", colnames(dat), fixed=T), 
                                   grep("AF", colnames(dat), fixed=T))))
  numeric.columns <- setdiff(numeric.columns, grep("SPAN", colnames(dat), fixed=T))
  dat[, numeric.columns] <- apply(dat[, numeric.columns], 2, as.numeric)
  # Read & add SD/SR coverage and QUAL score
  sd_sr_cov <- read.table(sd_sr_coverage.in, header=F, comment.char="")
  colnames(sd_sr_cov) <- c("name", "SDSR_COV")
  dat <- merge(dat, sd_sr_cov, all.x=T, all.y=F, by="name", sort=F)
  quals <- read.table(quals.in, header=T, comment.char="")
  colnames(quals)[1] <- c("name")
  dat <- merge(dat, quals, all.x=T, all.y=F, by="name", sort=F)
  # Read & add EUR LD results
  eur.ld <- read.table(eur_ld.in, header=T, sep="\t", comment.char="")
  eur.ld <- eur.ld[, which(colnames(eur.ld) %in% c("VID", "N_SNV_COMPARED", "MEAN_R2", "MEDIAN_R2", "MAX_R2",
                                                   "N_SNV_WITH_R2_0.8", "DIST_TO_NEAREST_R2_0.8",
                                                   "N_SNV_WITH_R2_0.2", "DIST_TO_NEAREST_R2_0.2",
                                                   "SNV_WITH_R2_0.8"))]
  colnames(eur.ld) <- c("name", "EUR_SNV_COMPARED", 
                        "EUR_MEAN_R2", "EUR_MEDIAN_R2", "EUR_MAX_R2", 
                        "EUR_N_SNV_R2_0.8", "EUR_DIST_R2_0.8", 
                        "EUR_N_SNV_R2_0.2", "EUR_DIST_R2_0.2",
                        "EUR_SNV_R2_0.8")
  eur.ld <- eur.ld[which(eur.ld$EUR_SNV_COMPARED >= minSNV), ]
  dat <- merge(dat, eur.ld, by="name", all.x=T, sort=F)
  # Read & add AFR LD results
  afr.ld <- read.table(afr_ld.in, header=T, sep="\t", comment.char="")
  afr.ld <- afr.ld[, which(colnames(afr.ld) %in% c("VID", "N_SNV_COMPARED", "MEAN_R2", "MEDIAN_R2", "MAX_R2",
                                                   "N_SNV_WITH_R2_0.8", "DIST_TO_NEAREST_R2_0.8",
                                                   "N_SNV_WITH_R2_0.2", "DIST_TO_NEAREST_R2_0.2",
                                                   "SNV_WITH_R2_0.8"))]
  colnames(afr.ld) <- c("name", "AFR_SNV_COMPARED", 
                        "AFR_MEAN_R2", "AFR_MEDIAN_R2", "AFR_MAX_R2", 
                        "AFR_N_SNV_R2_0.8", "AFR_DIST_R2_0.8", 
                        "AFR_N_SNV_R2_0.2", "AFR_DIST_R2_0.2",
                        "AFR_SNV_R2_0.8")
  afr.ld <- afr.ld[which(afr.ld$AFR_SNV_COMPARED >= minSNV), ]
  dat <- merge(dat, afr.ld, by="name", all.x=T, sort=F)
  # Remove variants not present in either EUR or AFR LD results
  dat <- dat[which(!is.na(dat$EUR_MAX_R2) | !is.na(dat$AFR_MAX_R2)), ]
  # Keep variants with either EUR or AFR AF >= minAF
  dat <- dat[which(dat$EUR_AF >= minAF | dat$AFR_AF >= minAF), ]
  # Annotate variants with GWAS catalog overlap
  gwas_hits <- read.table(gwas_hits.in, header=F)
  gwas_hits$IN_GWAS <- 1
  colnames(gwas_hits)[1] <- "name"
  dat <- merge(dat, gwas_hits, by="name", sort=F, all.x=T, all.y=F)
  return(dat)
}

# General function to plot peak R2 for subsets of variants, split by EUR/AFR
plot.rsq <- function(dat, indexes, pop.split=T, colors=NULL, h=0.05, bg=T){
  # Get plot data
  if(pop.split==F){
    plot.dat <- lapply(indexes, function(idxs){
      vals <- apply(dat[idxs, which(colnames(dat) %in% c("AFR_MAX_R2", "EUR_MAX_R2"))], 1, max, na.rm=T)
      return(vals[which(!is.na(vals))])
    })
  }else{
    plot.dat <- lapply(indexes, function(idxs){
      plot.dat.afr <- dat$AFR_MAX_R2[idxs]
      plot.dat.afr <- plot.dat.afr[which(!is.na(plot.dat.afr))]
      plot.dat.eur <- dat$EUR_MAX_R2[idxs]
      plot.dat.eur <- plot.dat.eur[which(!is.na(plot.dat.eur))]
      return(list(plot.dat.afr, plot.dat.eur))
    })
  }
  # Prep plot area
  par(mar=c(2.25, 3.5, 0.5, 0.5), bty="n")
  plot(x=c(0, length(indexes)), y=c(0, 1), type="n",
       xaxt="n", yaxt="n", xlab="", ylab="")
  if(bg==T){
    abline(h=c(0.8, 0.2), lty=2)
  }
  axis(1, at=par("usr")[1:2], tck=0, labels=NA)
  axis(2, at=seq(0, 1, 0.2), labels=NA, tck=-0.025)
  axis(2, at=seq(0, 1, 0.2), tick=F, las=2, line=-0.5)
  if(pop.split==F){
    mtext(2, text=bquote("Max" ~ R^{2} ~ "(AFR, EUR)"), line=1.75)
  }else{
    mtext(2, text=bquote("Max" ~ R^{2}), line=1.75)
  }
  # Add violins
  iqr.lwd <- 1
  med.cex <- 1
  if(pop.split==F){
    sapply(1:length(indexes), function(i){
      vioplot(plot.dat[[i]], add=T, at=i-0.5, col=colors[i], h=h, drawRect=F)
      segments(x0=i-0.5, x1=i-0.5, 
               y0=quantile(plot.dat[[i]], probs=0.25, na.rm=T),
               y1=quantile(plot.dat[[i]], probs=0.75, na.rm=T),
               lwd=iqr.lwd, lend="butt")
      points(x=i-0.5, y=median(plot.dat[[i]], na.rm=T), pch=19, cex=med.cex)
    })
  }else{
    sapply(1:length(indexes), function(i){
      # AFR
      vioplot(plot.dat[[i]][[1]], add=T, at=i-0.65, h=h, drawRect=F, wex=0.35,
              col=pops$color[which(pops$pop=="AFR")])
      segments(x0=i-0.65, x1=i-0.65, 
               y0=quantile(plot.dat[[i]][[1]], probs=0.25, na.rm=T),
               y1=quantile(plot.dat[[i]][[1]], probs=0.75, na.rm=T),
               lwd=iqr.lwd, lend="butt")
      points(x=i-0.65, y=median(plot.dat[[i]][[1]], na.rm=T), pch=19, cex=med.cex)
      # EUR
      vioplot(plot.dat[[i]][[2]], add=T, at=i-0.35, h=h, drawRect=F, wex=0.35,
              col=pops$color[which(pops$pop=="EUR")])
      segments(x0=i-0.35, x1=i-0.35, 
               y0=quantile(plot.dat[[i]][[2]], probs=0.25, na.rm=T),
               y1=quantile(plot.dat[[i]][[2]], probs=0.75, na.rm=T),
               lwd=iqr.lwd, lend="butt")
      points(x=i-0.35, y=median(plot.dat[[i]][[2]], na.rm=T), pch=19, cex=med.cex)
    })
  }
}

# Get average LD data by quality vs SV types
get.ldByQual <- function(dat, maxSDSR=0.30, require.sr=T, qual.bin.step=50, mode="cumulative"){
  qual.bins <- seq(0, 1000, qual.bin.step)
  if(require.sr==T){
    dat <- dat[grep("SR", dat$EVIDENCE, fixed=T), ]
  }
  dat.byQual <- lapply(c("ALL", "DEL", "DUP", "INS", "INVCPX"), function(svtype){
    if(svtype=="ALL"){
      if(mode=="cumulative"){
        kept <- sapply(qual.bins[-length(qual.bins)], function(minq){
          median(apply(dat[which(dat$SDSR_COV<=maxSDSR & dat$QUAL>=minq), 
                           which(colnames(dat) %in% c("AFR_MAX_R2", "EUR_MAX_R2"))], 
                       1, max, na.rm=T))
        })
        filtered <- sapply(qual.bins[-length(qual.bins)], function(minq){
          median(apply(dat[which(dat$SDSR_COV<=maxSDSR & dat$QUAL<minq), 
                           which(colnames(dat) %in% c("AFR_MAX_R2", "EUR_MAX_R2"))], 
                       1, max, na.rm=T))
        })
        return(list("kept"=kept, "filtered"=filtered))
      }else{
        sapply(2:length(qual.bins), function(i){
          median(apply(dat[which(dat$SDSR_COV<=maxSDSR & dat$QUAL>=qual.bins[i-1] & dat$QUAL<qual.bins[i]), 
                           which(colnames(dat) %in% c("AFR_MAX_R2", "EUR_MAX_R2"))], 
                       1, max, na.rm=T))
        })
      }
    }else if(svtype=="INVCPX"){
      if(mode=="cumulative"){
        kept <- sapply(qual.bins[-length(qual.bins)], function(minq){
          median(apply(dat[which(dat$SDSR_COV<=maxSDSR & dat$QUAL>=minq & dat$SVTYPE %in% c("INV", "CPX")), 
                           which(colnames(dat) %in% c("AFR_MAX_R2", "EUR_MAX_R2"))], 
                       1, max, na.rm=T))
        })
        filtered <- sapply(qual.bins[-length(qual.bins)], function(minq){
          median(apply(dat[which(dat$SDSR_COV<=maxSDSR & dat$QUAL<minq & dat$SVTYPE %in% c("INV", "CPX")), 
                           which(colnames(dat) %in% c("AFR_MAX_R2", "EUR_MAX_R2"))], 
                       1, max, na.rm=T))
        })
        return(list("kept"=kept, "filtered"=filtered))
      }else{
        sapply(2:length(qual.bins), function(i){
          median(apply(dat[which(dat$SDSR_COV<=maxSDSR & dat$QUAL>=qual.bins[i-1] & dat$QUAL<qual.bins[i] & dat$SVTYPE %in% c("INV", "CPX")), 
                           which(colnames(dat) %in% c("AFR_MAX_R2", "EUR_MAX_R2"))], 
                       1, max, na.rm=T))
        })
      }
    }else{
      if(mode=="cumulative"){
        kept <- sapply(qual.bins[-length(qual.bins)], function(minq){
          median(apply(dat[which(dat$SDSR_COV<=maxSDSR & dat$QUAL>=minq & dat$SVTYPE==svtype), 
                           which(colnames(dat) %in% c("AFR_MAX_R2", "EUR_MAX_R2"))], 
                       1, max, na.rm=T))
        })
        filtered <- sapply(qual.bins[-length(qual.bins)], function(minq){
          median(apply(dat[which(dat$SDSR_COV<=maxSDSR & dat$QUAL<minq & dat$SVTYPE==svtype), 
                           which(colnames(dat) %in% c("AFR_MAX_R2", "EUR_MAX_R2"))], 
                       1, max, na.rm=T))
        })
        return(list("kept"=kept, "filtered"=filtered))
      }else{
        sapply(2:length(qual.bins), function(i){
          median(apply(dat[which(dat$SDSR_COV<=maxSDSR & dat$QUAL>=qual.bins[i-1] & dat$QUAL<qual.bins[i] & dat$SVTYPE==svtype), 
                           which(colnames(dat) %in% c("AFR_MAX_R2", "EUR_MAX_R2"))], 
                       1, max, na.rm=T))
        })
      }
    }
  })
  names(dat.byQual) <- c("ALL", "DEL", "DUP", "INS", "INVCPX")
  return(dat.byQual)
}

# Plot LD as a function of QUAL, for a single SV type
plot.ldVsQual <- function(dat.byQual, svtype=NULL, qual.bin.step=50,
                              color="black", mode="cumulative"){
  # Get plot data
  qual.bins <- seq(0, 1000, qual.bin.step)
  if(is.null(svtype)){
    svtype <- "ALL"
  }
  if(mode=="cumulative"){
    keep.vals <- dat.byQual[[which(names(dat.byQual) == svtype)]]$kept
    filtered.vals <- dat.byQual[[which(names(dat.byQual) == svtype)]]$filtered
    vals <- c(keep.vals, filtered.vals)
    colors <- c("PASS"="#4DAC26", "FAIL"="#AC26A1")
  }else if(mode=="marginal"){
    vals <- dat.byQual[[which(names(dat.byQual) == svtype)]]
  }else{
    stop(paste("Specified mode (", mode, ") is not recognized", sep=""))
  }
  ymin <- max(c(0, min(vals, na.rm=T)-0.05))
  ymax <- min(c(1, max(vals, na.rm=T)+0.05))
  yrange <- c(ymin, ymax)
  # Prep plot area
  par(mar=c(2.25, 4.1, 0.5, 0.5))
  plot(x=c(0, 1000), y=yrange, type="n",
       xaxt="n", xlab="", yaxt="n", ylab="")
  # Add plot background
  if(mode=="cumulative"){
    rmean.keep.vals <- rollapply(keep.vals, 5, mean, partial=T, na.rm=T)
    rmean.filtered.vals <- rollapply(filtered.vals, 5, mean, partial=T, na.rm=T)
    points(x=(qual.bins[-length(qual.bins)]+qual.bins[-1])/2,
           y=rmean.filtered.vals, type="l", lwd=3, col=adjustcolor(colors[2], alpha=0.3))
    points(x=(qual.bins[-length(qual.bins)]+qual.bins[-1])/2,
           y=rmean.keep.vals, type="l", lwd=3, col=adjustcolor(colors[1], alpha=0.3))
  }else{
    abline(v=qual.bins, col="gray90", lwd=0.75)
  }
  # Add points
  if(mode=="cumulative"){
    points(x=(qual.bins[-length(qual.bins)]+qual.bins[-1])/2, 
           y=filtered.vals, pch=20, col=colors[2], cex=0.85)
    points(x=(qual.bins[-length(qual.bins)]+qual.bins[-1])/2, 
           y=keep.vals, pch=20, col=colors[1], cex=0.85)
  }else{
    points(x=(qual.bins[-length(qual.bins)]+qual.bins[-1])/2, 
           y=vals, pch=15, col=color, cex=0.85)
  }
  # Add axes
  y.lab <- bquote("Max" ~ R^{2} ~ "(SV, SNV)")
  if(mode=="cumulative"){
    x.lab <- "Min. QUAL Score"
  }else{
    x.lab <- "Binned QUAL Score"
  }
  x.at <- seq(0, 1000, 250)
  axis(1, at=x.at, labels=NA, tck=-0.03)
  sapply(x.at, function(x){
    axis(1, at=x, tick=F, line=-0.8, cex.axis=0.85)
  })
  mtext(1, text=x.lab, line=1.15)
  axis(2, at=axTicks(2), labels=NA, tck=-0.03)
  sapply(axTicks(2), function(y){
    axis(2, at=y, tick=F, las=2, line=-0.55, cex.axis=0.85)
  })
  mtext(2, text=y.lab, line=1.5)
  if(mode=="marginal"){
    mtext(2, text="Marginal", line=3)
  }
  box(bty="o")
}

# Barplot of GWAS tagged variants by consequence
gwas.hits.barplot <- function(dat){
  # Get data
  gwas.tagged.idxs <- which(dat$SDSR_COV<=0.3 
                            & (dat$EUR_MAX_R2>=0.8 | dat$AFR_MAX_R2>=0.8)
                            & !is.na(dat$IN_GWAS))
  lof.idxs <- intersect(gwas.tagged.idxs, which(!is.na(dat$PROTEIN_CODING__LOF)))
  partialdup.idxs <- setdiff(intersect(gwas.tagged.idxs, which(!is.na(dat$PROTEIN_CODING__DUP_PARTIAL))), lof.idxs)
  utr.idxs <- setdiff(intersect(gwas.tagged.idxs, which(!is.na(dat$PROTEIN_CODING__UTR))), c(lof.idxs, partialdup.idxs))
  promoter.idxs <- setdiff(intersect(gwas.tagged.idxs, which(!is.na(dat$PROTEIN_CODING__PROMOTER))), c(lof.idxs, partialdup.idxs, utr.idxs))
  intronic.idxs <- setdiff(intersect(gwas.tagged.idxs, which(!is.na(dat$PROTEIN_CODING__INTRONIC))), 
                           c(lof.idxs, partialdup.idxs, utr.idxs, promoter.idxs))
  intergenic.idxs <- setdiff(intersect(gwas.tagged.idxs, which(dat$PROTEIN_CODING__INTERGENIC=="True")), 
                           c(lof.idxs, partialdup.idxs, utr.idxs, promoter.idxs, intronic.idxs))
  counts <- data.frame("csq"=c("pLoF", "Partial CG", "UTR", "Promoter", "Intronic", "Intergenic"),
                       "count"=c(length(lof.idxs), length(partialdup.idxs), length(utr.idxs),
                                 length(promoter.idxs), length(intronic.idxs), length(intergenic.idxs)),
                       "color"=c(svtypes$color[which(svtypes$svtype=="DEL")],
                                 svtypes$color[which(svtypes$svtype=="DUP")],
                                 "gray25", "gray50", "gray75", "white"))
  # Prep plot area
  par(mar=c(0.75, 3, 0.75, 7.5), bty="n")
  plot(x=c(0, 1), y=c(0, sum(counts$count)), type="n",
       xaxt="n", yaxt="n", xlab="", ylab="")
  axis(2, at=axTicks(2), labels=NA, tck=-0.06)
  axis(2, at=axTicks(2), tick=F, las=2, line=-0.6, cex.axis=0.8)
  mtext(2, line=1.8, text="SVs in Strong LD with GWAS Loci")
  # Add rectangles
  rect(xleft=0, xright=0.75, bty="n", border=NA,
       ybottom=c(0, cumsum(counts$count[-nrow(counts)])),
       ytop=cumsum(counts$count),
       col=counts$color)
  rect(xleft=0, xright=0.75, ybottom=0, ytop=sum(counts$count), col=NA)
  # Add points and connectors
  pt.at <- seq(0, sum(counts$count), by=sum(counts$count)/(nrow(counts)+1))
  pt.at <- pt.at[-c(1, length(pt.at))]
  # sapply(1:nrow(counts), function(i){
  #   mdpt <- mean(c(0, cumsum(counts$count))[c(i, i+1)])
  #   segments(x0=c(0.75, 0.8, 0.9), x1=c(0.8, 0.9, 0.95),
  #            y0=c(mdpt, mdpt, pt.at[i]),
  #            y1=c(mdpt, pt.at[i], pt.at[i]))
  # })
  points(x=rep(0.95, nrow(counts)), y=pt.at, pch=22, bg=counts$color, cex=1.5)
  axis(4, at=pt.at, tick=F, line=-0.8, 
       labels=paste(counts$count, counts$csq, "SVs"), las=2)
}


# Gather functional enrichments for GWAS hits
get.gwas.enrichment <- function(dat, eur.only=T){
  # Get baseline data
  if(eur.only==T){
    tagged.all.idxs <- which(dat$SDSR_COV<=0.3 & dat$EUR_MAX_R2>=0.8)
  }else{
    tagged.all.idxs <- which(dat$SDSR_COV<=0.3 & (dat$EUR_MAX_R2>=0.8 | dat$AFR_MAX_R2>=0.8))
  }
  n.tagged.all <- length(tagged.all.idxs)
  tagged.gwas.idxs <- intersect(tagged.all.idxs, which(!is.na(dat$IN_GWAS)))
  n.tagged.gwas <- length(tagged.gwas.idxs)
  tagged.notgwas.idxs <- intersect(tagged.all.idxs, which(is.na(dat$IN_GWAS)))
  n.tagged.notgwas <- length(tagged.notgwas.idxs)
  # Calculate enrichment per annotation
  annotations <- c("LOF", "DUP_LOF", "COPY_GAIN", "PROMOTER", "UTR",
                   "DUP_PARTIAL", "INTRONIC", "INTERGENIC")
  annotations <- paste("PROTEIN_CODING__", annotations, sep="")
  sapply(annotations, function(anno){
    # Get counts
    col.idx <- which(colnames(dat)==anno)
    if(anno=="PROTEIN_CODING__INTERGENIC"){
      hits.idxs <- which(dat[, col.idx]=="True")
    }else{
      hits.idxs <- which(!is.na(dat[, col.idx]))
    }
    tagged.hits.idxs <- intersect(tagged.all.idxs, hits.idxs)
    tagged.gwas.hits.idxs <- intersect(tagged.gwas.idxs, tagged.hits.idxs)
    n.tagged.gwas.hits <- length(tagged.gwas.hits.idxs)
    tagged.notgwas.hits.idxs <- intersect(tagged.notgwas.idxs, tagged.hits.idxs)
    n.tagged.notgwas.hits <- length(tagged.notgwas.hits.idxs)
    # Calculate enrichment
    f.res <- fisher.test(matrix(c(n.tagged.notgwas - n.tagged.notgwas.hits, n.tagged.notgwas.hits,
                         n.tagged.gwas - n.tagged.gwas.hits, n.tagged.gwas.hits),
           nrow=2, byrow=F))
    c(n.tagged.gwas.hits, f.res$estimate, f.res$conf.int, f.res$p.value)
  })
}

# Dotplots of GWAS functional enrichments
plot.gwas.enrich <- function(gwas.enrich){
  # Set plot parameters
  func.metadata <- data.frame("class"=paste("PROTEIN_CODING__",
                                            c("LOF", "DUP_LOF", "COPY_GAIN", "PROMOTER", "UTR",
                                              "DUP_PARTIAL", "INTRONIC", "INTERGENIC"),
                                            sep=""),
                              "label"=c("pLoF", "CG", "IED", "Promoter", "UTR", 
                                        "Partial CG", "Intronic", "Intergenic"),
                              "color"=c(svtypes$color[which(svtypes$svtype=="DEL")],
                                        svtypes$color[which(svtypes$svtype=="MCNV")],
                                        svtypes$color[which(svtypes$svtype=="DUP")],
                                        rep("gray35", 2),
                                        svtypes$color[which(svtypes$svtype=="DUP")],
                                        "gray35", "gray60"))
  cols.to.keep <- which(apply(gwas.enrich, 2, function(col){!any(is.infinite(col))}))
  plot.dat <- gwas.enrich[-1, cols.to.keep]
  rownames(plot.dat) <- c("OR", "lower_CI", "upper_CI", "p")
  plot.dat <- plot.dat[, order(-plot.dat[1, ])]
  
  # Set plot area
  par(bty="n", mar=c(0.5, 4, 2.5, 0.5))
  plot(x=c(0, 1.5*max(plot.dat[1, ])), y=c(-0.4, -ncol(plot.dat)), type="n",
       xaxt="n", yaxt="n", xlab="", ylab="")
  abline(v=0)
  abline(v=1, lty=2)
  axis(3, at=axTicks(3), tck=-0.025, labels=NA)
  axis(3, at=axTicks(3), tick=F, line=-0.7, cex.axis=0.85)
  mtext(3, line=1.1, text="Odds Ratio")
  
  # Add points & CIs with labels
  sapply(1:ncol(plot.dat), function(i){
    class <- colnames(plot.dat)[i]
    color.i <- func.metadata$color[which(func.metadata$class==class)]
    label.i <- func.metadata$label[which(func.metadata$class==class)]
    segments(x0=plot.dat[2, i], x1=plot.dat[3, i],
             y0=-i+0.25, y1=-i+0.25, lend="round", col=color.i)
    points(x=plot.dat[1, i], y=-i+0.25, pch=19, col=color.i)
    axis(2, at=-i+0.25, tick=F, line=-0.9, labels=label.i, las=2,
         col.axis=color.i)
    pval <- plot.dat[4, i]
    if(pval<0.05/ncol(plot.dat)){
      sig <- "***"
    }else if(pval<0.05){
      sig <- "*"
    }else{
      sig <- ""
    }
    text(x=plot.dat[1, i], y=-i-0.1, pos=3, labels=sig, cex=1.25)
  })
}


#################
### RSCRIPT BLOCK
#################
require(optparse, quietly=T)
require(vioplot, quietly=T)
require(zoo, quietly=T)
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
              metavar="character"), 
  make_option(c("-m",  "--minAF"),  type="numeric",  default=0.01, 
              help="minimum SV AF to include in analysis [default %default]", 
              metavar="numeric")
)

### Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog VCF2BED EUR_LD_DATA AFR_LD_DATA SD_SR_COV QUALS OUTDIR", 
                                option_list=option_list), 
                   positional_arguments=TRUE)
opts <- args$options

### Checks for appropriate positional arguments
if(length(args$args) != 7){
  stop("Incorrect number of required positional arguments\n")
}

### Writes args & opts to vars
vcf2bed.in <- args$args[1]
eur_ld.in <- args$args[2]
afr_ld.in <- args$args[3]
gwas_hits.in <- args$args[4]
sd_sr_coverage.in <- args$args[5]
quals.in <- args$args[6]
OUTDIR <- args$args[7]
prefix <- opts$prefix
svtypes.file <- opts$svtypes
pops.file <- opts$populations
minAF <- opts$minAF

# #Dev parameters (local)
# vcf2bed.in <- "~/scratch/gnomAD-SV_v2_rev1.vcf2bed.bed.gz"
# eur_ld.in <- "~/scratch/gnomAD.EUR.SNV_SV.LD.results.bed.gz"
# afr_ld.in <- "~/scratch/gnomAD.AFR.SNV_SV.LD.results.bed.gz"
# gwas_hits.in <- "~/scratch/gnomAD-SV_v2_rev1.SVs_in_LD_to_UKBB_or_GWAS_cat.vids.list"
# sd_sr_coverage.in <- "~/scratch/gnomAD-SV_v2_rev1.variant_SD_SR_coverage.txt.gz"
# quals.in <- "~/scratch/gnomAD-SV_v2_rev1.QUAL_per_SV.txt.gz"
# OUTDIR <- "~/scratch/gnomAD_v2_test_plots/"
# svtypes.file <- "~/Desktop/Collins/Talkowski/code/sv-pipeline/ref/vcf_qc_refs/SV_colors.txt"
# pops.file <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/gnomAD_population_colors.txt"
# prefix <- "gnomAD-SV_v2_rev1"
# minAF <- 0.01


### Create output directory,  if necessary
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

### Process input data
dat <- import.data(vcf2bed.in, eur_ld.in, afr_ld.in, gwas_hits.in, 
                   sd_sr_coverage.in, quals.in,
                   minAF=minAF, minSNV=2000,
                   exclude.allosomes=T)

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


### Format supplementary table of SVs in LD with SNVs
tab.dat <- dat[which(dat$SDSR_COV<0.3 & (dat$AFR_N_SNV_R2_0.8>0 | dat$EUR_N_SNV_R2_0.8>0)), ]
tab.out <- data.frame("Chrom"=tab.dat$chrom, "SV Start"=tab.dat$start, "SV End"=tab.dat$end,
                      "SV Class"=tab.dat$SVTYPE, "SV ID"=tab.dat$name,
                      "AFR AF"=tab.dat$AFR_AF, "EUR AF"=tab.dat$EUR_AF,
                      "AFR Max R2 to SNV or indel"=tab.dat$AFR_MAX_R2,
                      "EUR Max R2 to SNV or indel"=tab.dat$EUR_MAX_R2,
                      "AFR SNV or indel in Strong LD"=tab.dat$AFR_SNV_R2_0.8, 
                      "EUR SNV or indel in Strong LD"=tab.dat$EUR_SNV_R2_0.8)
tab.out <- tab.out[order(tab.out$Chrom, tab.out$SV.Start, tab.out$SV.End), ]
colnames(tab.out) <- sapply(colnames(tab.out), function(str){gsub(".", " ", str, fixed=T)})
write.table(tab.out, paste(OUTDIR, "/", prefix, ".SV_SNV_strongLD_table.txt", sep=""),
            quote=F, sep="\t", col.names=T, row.names=F)
dat <- dat[, -which(colnames(dat) %in% c("AFR_SNV_R2_0.8", "EUR_SNV_R2_0.8"))]


### Plot LD by SDSR coverage bin
idx.by.SDSR <- lapply(seq(0, 0.8, 0.2), function(sdsrmin){
  which(dat$SDSR_COV>=sdsrmin & dat$SDSR_COV<=(sdsrmin+0.2))
})
# Both pops
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD.bySDSRcov.pdf", sep=""), height=2.5, width=4.5)
plot.rsq(dat, idx.by.SDSR, pop.split=F, colors=rep("gray50", 5))
sapply(1:5, function(i){
  axis(1, at=i-0.5, line=-0.8, tick=F, cex.axis=0.8,
       labels=paste(20*(i-1), "-", 20*i, "%", sep=""))
})
mtext(1, line=1.25, text="SD/SR Coverage")
dev.off()
# Split by pop
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD.bySDSRcov.popsplit.pdf", sep=""), height=2.5, width=4.5)
plot.rsq(dat, idx.by.SDSR, pop.split=T)
sapply(1:5, function(i){
  axis(1, at=i-0.5, line=-0.8, tick=F, cex.axis=0.8,
       labels=paste(20*(i-1), "-", 20*i, "%", sep=""))
})
mtext(1, line=1.25, text="SD/SR Coverage")
dev.off()


### Plot LD by SV type
idx.by.svtype <- lapply(c("ALL", "DEL", "DUP", "INS", "INVCPX"), function(svtype){
  if(svtype=="INVCPX"){
    which(dat$SVTYPE %in% c("INV", "CPX"))
  }else if(svtype=="ALL"){
    1:nrow(dat)
  }else{
    which(dat$SVTYPE==svtype)
  }
})
colors.by.svtype <- c("gray50", sapply(c("DEL", "DUP", "INS", "INV"), function(svtype){
  svtypes$color[which(svtypes$svtype==svtype)]
}))
# Both pops
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD.bySVtype.pdf", sep=""), height=2.5, width=4.5)
plot.rsq(dat, idx.by.svtype, pop.split=F, colors=colors.by.svtype)
sapply(1:5, function(i){
  axis(1, at=i-0.5, line=-0.8, tick=F, 
       labels=c("All SVs", "DEL", "DUP", "INS", "INV & CPX")[i])
})
mtext(1, line=1.25, text="SV Class")
dev.off()
# Split by pop
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD.bySVtype.popsplit.pdf", sep=""), height=2.5, width=4.5)
plot.rsq(dat, idx.by.svtype, pop.split=T)
sapply(1:5, function(i){
  axis(1, at=i-0.5, line=-0.8, tick=F, 
       labels=c("All SVs", "DEL", "DUP", "INS", "INV & CPX")[i])
})
mtext(1, line=1.25, text="SV Class")
dev.off()


### Plot LD by SV type, clean regions only
idx.by.svtype.clean <- lapply(c("ALL", "DEL", "DUP", "INS", "INVCPX"), function(svtype){
  if(svtype=="INVCPX"){
    which(dat$SVTYPE %in% c("INV", "CPX") & dat$SDSR_COV<0.30)
  }else if(svtype=="ALL"){
    which(dat$SDSR_COV<0.30)
  }else{
    which(dat$SVTYPE==svtype & dat$SDSR_COV<0.30)
  }
})
# Both pops
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD.bySVtype.clean_only.pdf", sep=""), height=2.5, width=4.5)
plot.rsq(dat, idx.by.svtype.clean, pop.split=F, colors=colors.by.svtype)
sapply(1:5, function(i){
  axis(1, at=i-0.5, line=-0.8, tick=F, 
       labels=c("All SVs", "DEL", "DUP", "INS", "INV & CPX")[i])
})
mtext(1, line=1.25, text="SV Class")
dev.off()
# Split by pop
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD.bySVtype.popsplit.clean_only.pdf", sep=""), height=2.5, width=4.5)
plot.rsq(dat, idx.by.svtype.clean, pop.split=T)
sapply(1:5, function(i){
  axis(1, at=i-0.5, line=-0.8, tick=F, 
       labels=c("All SVs", "DEL", "DUP", "INS", "INV & CPX")[i])
})
mtext(1, line=1.25, text="SV Class")
dev.off()


### PLot LD by QUAL score bin
idx.by.qual <- lapply(seq(0, 900, 100), function(qmin){
  which(dat$QUAL>qmin & dat$QUAL <= qmin+100)
})
# Both pops
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD.byQUAL.pdf", sep=""), height=2.5, width=6)
plot.rsq(dat, idx.by.qual, pop.split=F, colors=rep("gray50", 10))
sapply(1:10, function(i){
  axis(1, at=i-0.5, line=-0.8, tick=F, cex.axis=0.8, 
       labels=paste(i-1, "-", i, sep=""))
})
mtext(1, line=1.25, text="QUAL Score (x100)")
dev.off()
# Split by pop
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD.byQUAL.popsplit.pdf", sep=""), height=2.5, width=6)
plot.rsq(dat, idx.by.qual, pop.split=T)
sapply(1:10, function(i){
  axis(1, at=i-0.5, line=-0.8, tick=F, cex.axis=0.8, 
       labels=paste(i-1, "-", i, sep=""))
})
mtext(1, line=1.25, text="QUAL Score (x100)")
dev.off()


### PLot LD by QUAL score bin, clean regions only
idx.by.qual.clean <- lapply(seq(0, 900, 100), function(qmin){
  which(dat$QUAL>qmin & dat$QUAL <= qmin+100 & dat$SDSR_COV<0.30)
})
# Both pops
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD.byQUAL.clean_only.pdf", sep=""), height=2.5, width=6)
plot.rsq(dat, idx.by.qual.clean, pop.split=F, colors=rep("gray50", 10))
sapply(1:10, function(i){
  axis(1, at=i-0.5, line=-0.8, tick=F, cex.axis=0.8, 
       labels=paste(i-1, "-", i, sep=""))
})
mtext(1, line=1.25, text="QUAL Score (x100)")
dev.off()
# Split by pop
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD.byQUAL.popsplit.clean_only.pdf", sep=""), height=2.5, width=6)
plot.rsq(dat, idx.by.qual.clean, pop.split=T)
sapply(1:10, function(i){
  axis(1, at=i-0.5, line=-0.8, tick=F, cex.axis=0.8, 
       labels=paste(i-1, "-", i, sep=""))
})
mtext(1, line=1.25, text="QUAL Score (x100)")
dev.off()


### Plot LD by size
idx.by.size <- lapply(list(c(-2, 500), c(500, 1000), c(1000, 5000), 
                           c(5000, 50000), c(50000, 300000000)), function(sizelims){
  which(dat$SVLEN>=sizelims[1] & dat$SVLEN<sizelims[2])
})
# Both pops
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD.bysize.pdf", sep=""), height=2.5, width=4.5)
plot.rsq(dat, idx.by.size, pop.split=F, colors=rep("gray50", 5))
sapply(1:5, function(i){
  axis(1, at=i-0.5, line=-0.8, tick=F, 
       labels=c("0-0.5", "0.5-1", "1-5", "5-50", ">50")[i])
})
mtext(1, line=1.25, text="SV Size (kb)")
dev.off()
# Both pops
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD.bysize.popsplit.pdf", sep=""), height=2.5, width=4.5)
plot.rsq(dat, idx.by.size, pop.split=T)
sapply(1:5, function(i){
  axis(1, at=i-0.5, line=-0.8, tick=F, 
       labels=c("0-0.5", "0.5-1", "1-5", "5-50", ">50")[i])
})
mtext(1, line=1.25, text="SV Size (kb)")
dev.off()


### Plot LD by size, clean regions only
idx.by.size.clean <- lapply(list(c(-2, 500), c(500, 1000), c(1000, 5000), 
                           c(5000, 50000), c(50000, 300000000)), function(sizelims){
                             which(dat$SVLEN>=sizelims[1] & dat$SVLEN<sizelims[2] & dat$SDSR_COV<0.30)
                           })
# Both pops
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD.bysize.clean_only.pdf", sep=""), height=2.5, width=4.5)
plot.rsq(dat, idx.by.size.clean, pop.split=F, colors=rep("gray50", 5))
sapply(1:5, function(i){
  axis(1, at=i-0.5, line=-0.8, tick=F, 
       labels=c("0-0.5", "0.5-1", "1-5", "5-50", ">50")[i])
})
mtext(1, line=1.25, text="SV Size (kb)")
dev.off()
# Both pops
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD.bysize.popsplit.clean_only.pdf", sep=""), height=2.5, width=4.5)
plot.rsq(dat, idx.by.size.clean, pop.split=T)
sapply(1:5, function(i){
  axis(1, at=i-0.5, line=-0.8, tick=F, 
       labels=c("0-0.5", "0.5-1", "1-5", "5-50", ">50")[i])
})
mtext(1, line=1.25, text="SV Size (kb)")
dev.off()


### Plot LD by AF
idx.by.af <- lapply(list(c(0, 0.02), c(0.02, 0.05), c(0.05, 0.1), 
                         c(0.1, 0.5), c(0.5, 1)), function(aflims){
                           which(dat$AF>aflims[1] & dat$AF<=aflims[2])
                         })
# Both pops
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD.byAF.pdf", sep=""), height=2.5, width=4.5)
plot.rsq(dat, idx.by.af, pop.split=F, colors=rep("gray50", 5))
sapply(1:5, function(i){
  axis(1, at=i-0.5, line=-0.8, tick=F, 
       labels=c("<2%", "2-5%", "5-10%", "10-50%", ">50%")[i])
})
mtext(1, line=1.25, text="Global SV AF")
dev.off()
# Both pops
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD.byAF.popsplit.pdf", sep=""), height=2.5, width=4.5)
plot.rsq(dat, idx.by.af, pop.split=T)
sapply(1:5, function(i){
  axis(1, at=i-0.5, line=-0.8, tick=F, 
       labels=c("<2%", "2-5%", "5-10%", "10-50%", ">50%")[i])
})
mtext(1, line=1.25, text="Global SV AF")
dev.off()


### Plot LD by AF, clean regions only
idx.by.af.clean <- lapply(list(c(0, 0.02), c(0.02, 0.05), c(0.05, 0.1), 
                         c(0.1, 0.5), c(0.5, 1)), function(aflims){
                           which(dat$AF>aflims[1] & dat$AF<=aflims[2] & dat$SDSR_COV<0.30)
                         })
# Both pops
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD.byAF.clean_only.pdf", sep=""), height=2.5, width=4.5)
plot.rsq(dat, idx.by.af.clean, pop.split=F, colors=rep("gray50", 5))
sapply(1:5, function(i){
  axis(1, at=i-0.5, line=-0.8, tick=F, 
       labels=c("<2%", "2-5%", "5-10%", "10-50%", ">50%")[i])
})
mtext(1, line=1.25, text="Global SV AF")
dev.off()
# Both pops
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD.byAF.popsplit.clean_only.pdf", sep=""), height=2.5, width=4.5)
plot.rsq(dat, idx.by.af.clean, pop.split=T)
sapply(1:5, function(i){
  axis(1, at=i-0.5, line=-0.8, tick=F, 
       labels=c("<2%", "2-5%", "5-10%", "10-50%", ">50%")[i])
})
mtext(1, line=1.25, text="Global SV AF")
dev.off()


### Plot LD by QUAL for QUAL calibration supplementary figure, in clean regions only
dat.byQual.cumulative <- get.ldByQual(dat, maxSDSR=0.30, require.sr=T, qual.bin.step=50, mode="cumulative")
dat.byQual.marginal <- get.ldByQual(dat, maxSDSR=0.30, require.sr=T, qual.bin.step=50, mode="marginal")
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD_byQUAL.all.cumulative.pdf", sep=""),
    height=2, width=2.25)
plot.ldVsQual(dat.byQual.cumulative,  svtype="ALL", mode="cumulative")
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD_byQUAL.all.marginal.pdf", sep=""),
    height=2, width=2.25)
plot.ldVsQual(dat.byQual.marginal,  svtype="ALL", mode="marginal")
dev.off()
sapply(c("DEL", "DUP", "INS", "INVCPX"), function(svtype){
  pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD_byQUAL.", svtype, ".cumulative.pdf", sep=""),
      height=2, width=2.25)
  plot.ldVsQual(dat.byQual.cumulative,  svtype=svtype, mode="cumulative")
  dev.off()
  pdf(paste(OUTDIR, "/", prefix, ".SNV_SV_LD_byQUAL.", svtype, ".marginal.pdf", sep=""),
      height=2, width=2.25)
  if(svtype=="INVCPX"){
    plot.ldVsQual(dat.byQual.marginal,  svtype=svtype, mode="marginal",
                      color=svtypes$color[which(svtypes$svtype=="CPX")])
  }else{
    plot.ldVsQual(dat.byQual.marginal,  svtype=svtype, mode="marginal",
                      color=svtypes$color[which(svtypes$svtype==svtype)])
  }
  dev.off()
})

### Calculate fraction of SVs with at least one strong tagged SNP
# Clean regions only
length(which(apply(dat[which(dat$SDSR_COV<=0.3), which(colnames(dat) %in% c("AFR_MAX_R2", "EUR_MAX_R2"))], 1, max, na.rm=T)>=0.2))/nrow(dat)
length(which(apply(dat[which(dat$SDSR_COV<=0.3), which(colnames(dat) %in% c("AFR_MAX_R2", "EUR_MAX_R2"))], 1, max, na.rm=T)>=0.8))/nrow(dat)
#Raw count
length(which(apply(dat[which(dat$SDSR_COV<=0.3), which(colnames(dat) %in% c("AFR_MAX_R2", "EUR_MAX_R2"))], 1, max, na.rm=T)>=0.2))
length(which(apply(dat[which(dat$SDSR_COV<=0.3), which(colnames(dat) %in% c("AFR_MAX_R2", "EUR_MAX_R2"))], 1, max, na.rm=T)>=0.8))

### Calculate median maximum correlation coefficient in clean regions
# All SVs
vals <- apply(dat[which(dat$SDSR_COV<=0.3), which(colnames(dat) %in% c("AFR_MAX_R2", "EUR_MAX_R2"))], 1, max, na.rm=T)
c(length(vals), median(vals, na.rm=T))
# By SV type (DEL/DUP/INS)
sapply(c("DEL", "DUP", "INS"), function(svtype){
  vals <- apply(dat[which(dat$SDSR_COV<=0.3 & dat$SVTYPE==svtype), which(colnames(dat) %in% c("AFR_MAX_R2", "EUR_MAX_R2"))], 1, max, na.rm=T)
  c(length(vals), median(vals, na.rm=T))
})
# INV & CPX
vals <- apply(dat[which(dat$SDSR_COV<=0.3 & dat$SVTYPE %in% c("INV", "CPX")), which(colnames(dat) %in% c("AFR_MAX_R2", "EUR_MAX_R2"))], 1, max, na.rm=T)
c(length(vals), median(vals, na.rm=T))
# All SVs, EUR vs AFR
median(dat$EUR_MAX_R2[which(dat$SDSR_COV<=0.3)], na.rm=T)
median(dat$AFR_MAX_R2[which(dat$SDSR_COV<=0.3)], na.rm=T)
format(t.test(dat$AFR_MAX_R2[which(dat$SDSR_COV<=0.3)], 
                   dat$EUR_MAX_R2[which(dat$SDSR_COV<=0.3)], 
                   alternative="less")$p.value, scientific=T)


### Count number of GWAS-linked SVs
length(which(!is.na(dat$IN_GWAS) & dat$SDSR_COV<=0.3))


### Barplot of GWAS catalog breakdown
pdf(paste(OUTDIR, "/", prefix, ".GWAS_tagged_SVs.func_consequence.barplot.pdf", sep=""),
    height=3, width=3)
gwas.hits.barplot(dat)
dev.off()

### Functional enrichments of GWAS catalog hits
gwas.enrich <- get.gwas.enrichment(dat, eur.only=T)
pdf(paste(OUTDIR, "/", prefix, ".GWAS_tagged_SVs.func_enrich.dotplot.pdf", sep=""),
    height=2.5, width=2.5)
plot.gwas.enrich(gwas.enrich)
dev.off()

# ### Write list of candidate variants for analysis
# hq.candidates <- dat$name[intersect(grep("RD", dat$EVIDENCE, fixed=T),
#                                     intersect(grep("SR", dat$EVIDENCE, fixed=T),
#                                               which((!is.na(dat$PROTEIN_CODING__PROMOTER)
#                                                      | !is.na(dat$PROTEIN_CODING__LOF)
#                                                      | !is.na(dat$PROTEIN_CODING__UTR))
#                                                     & dat$SVTYPE=="DEL" & dat$EUR_MAX_R2>0.8
#                                                     & !is.na(dat$IN_GWAS)
#                                                     & dat$SDSR_COV<0.05)))]
# write.table(hq.candidates, 
#             paste(OUTDIR, "/", prefix, ".highQuality_GWAS_deletions_forEvaluation.txt", sep=""),
#             sep="\t", col.names=F, quote=F, row.names=F)


