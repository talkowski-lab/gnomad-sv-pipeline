#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Plot de novo rates by SVTYPE for formal gnomAD analysis


###Set master parameters
options(stringsAsFactors=F, scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Read & clean de novo rate data
importDeNovoRates <- function(dat.in, svtypes.forAnalysis, PCRMINUS.samples){
  dat <- read.table(dat.in, header=T, comment.char="")
  colnames(dat)[1] <- "proband"
  filts <- c("PASS", "FAIL")
  #Get unique variant categories
  cats <- unique(unlist(lapply(strsplit(colnames(dat)[grep("PASS_", colnames(dat), fixed=T)], split="PASS_"), function(l){l[[2]]})))
  #Consolidate by SVTYPE groups
  for(cat in cats){
    for(filt in filts){
      #Get indexes
      cnv.idxs <- sapply(c("DEL", "DUP"), function(svtype){
        grep(paste(svtype, filt, cat, sep="_"), colnames(dat), fixed=T)
      })
      noncnv.idxs <- sapply(c("INS", "INV", "CPX", "CTX"), function(svtype){
        grep(paste(svtype, filt, cat, sep="_"), colnames(dat), fixed=T)
      })
      #Create new columns
      cnv.newvals <- apply(dat[, cnv.idxs], 1, sum)
      noncnv.newvals <- apply(dat[, noncnv.idxs], 1, sum)
      dat <- cbind(dat, 
                   cnv.newvals, 
                   noncnv.newvals)
      colnames(dat)[c(ncol(dat)-1, ncol(dat))] <- c(paste("CNV", filt, cat, sep="_"), 
                                                   paste("BCA", filt, cat, sep="_"))
    }
  }
  #Drop unnecessary columns
  cols.to.drop <- sort(unique(as.vector(sapply(c("DEL", "DUP", "INS", "INV", "CPX", "CTX"), function(svtype){
    grep(paste(svtype, "_", sep=""), colnames(dat), fixed=T)
  }))))
  if(length(cols.to.drop)>0){
    dat <- dat[, -cols.to.drop]
  }
  #Calculate various summary rates
  for(svtype in c("ALL", "CNV", "BCA", "BND")){
    for(filt in filts){
      #Mendelian variants: MENDELIAN_CHILD_HET, MENDELIAN_CHILD_HOM, UNTRANSMITTED_HET
      mendelian.cats <- c("MENDELIAN_CHILD_HET", "MENDELIAN_CHILD_HOM", "UNTRANSMITTED_HET")
      mendelian.idxs <- sapply(mendelian.cats, function(cat){
        grep(paste(svtype, filt, cat, sep="_"), colnames(dat), fixed=T)
      })
      mendelian.newvals <- apply(dat[, mendelian.idxs], 1, sum)
      #Mendelian violations: APPARENT_DE_NOVO_HET, APPARENT_DE_NOVO_HOM, UNTRANSMITTED_HOM
      violation.cats <- c("APPARENT_DE_NOVO_HET", "APPARENT_DE_NOVO_HOM", "UNTRANSMITTED_HOM")
      violation.idxs <- sapply(violation.cats, function(cat){
        grep(paste(svtype, filt, cat, sep="_"), colnames(dat), fixed=T)
      })
      violation.newvals <- apply(dat[, violation.idxs], 1, sum)
      #Mendelian violation rate
      mvr.denom.newvals <- (mendelian.newvals+violation.newvals)
      mvr.newvals <- violation.newvals/mvr.denom.newvals
      #Het de novo rate
      het_mendel <- dat[, grep(paste(svtype, filt, "MENDELIAN_CHILD_HET", sep="_"), colnames(dat), fixed=T)]
      het_denovo <- dat[, grep(paste(svtype, filt, "APPARENT_DE_NOVO_HET", sep="_"), colnames(dat), fixed=T)]
      het_dnr.denom.newvals <- (het_mendel+het_denovo)
      het_dnr.newvals <- het_denovo/het_dnr.denom.newvals
      #Hom genotype FDR
      hom_mendel <- dat[, grep(paste(svtype, filt, "MENDELIAN_CHILD_HOM", sep="_"), colnames(dat), fixed=T)]
      hom_denovo <- dat[, grep(paste(svtype, filt, "APPARENT_DE_NOVO_HOM", sep="_"), colnames(dat), fixed=T)]
      hom_fdr.denom.newvals <- (hom_mendel+hom_denovo)
      hom_fdr.newvals <- hom_denovo/hom_fdr.denom.newvals
      #Het genotype FNR
      parent_hom_untrans <- dat[, grep(paste(svtype, filt, "UNTRANSMITTED_HOM", sep="_"), colnames(dat), fixed=T)]
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
      new.colnames <- as.vector(sapply(c("MVR", "HET_DNR", "SPONT_HOMOZYGOTE_RATE", "HET_FNR"), function(cat){
        c(paste(svtype, filt, cat, "SV_COUNT", sep="_"), 
          paste(svtype, filt, cat, sep="_"))
      }))
      colnames(dat)[(ncol(dat)-7):ncol(dat)] <- new.colnames
    }
  }
  return(dat)
}

# Read & clean de novo data, binned by QUAL score
import.denovoByQual <- function(byqual.in, mode="marginal"){
  df <- read.table(byqual.in, header=T, sep="\t", comment.char="")[, -1]
  res <- lapply(c("ALL", "DEL", "DUP", "INS", "INV_CPX"), function(svtype){
    if(svtype != "INV_CPX"){
      dn.idx <- intersect(grep(paste(svtype, "_", sep=""), colnames(df)),
                          grep("_DE_NOVO_", colnames(df)))
      inh.idx <- intersect(grep(paste(svtype, "_", sep=""), colnames(df)),
                           grep("_MENDELIAN_", colnames(df)))
      minQuals <- unlist(lapply(strsplit(colnames(df)[dn.idx], split="_"), function(v){as.numeric(v[[2]])}))
      dn.k.tab <- df[, dn.idx]
      inh.k.tab <- df[, inh.idx]
      all.k.tab <- dn.k.tab + inh.k.tab
    }else{
      dn.idx.inv <- intersect(grep("INV_", colnames(df)),
                          grep("_DE_NOVO_", colnames(df)))
      dn.idx.cpx <- intersect(grep("CPX_", colnames(df)),
                              grep("_DE_NOVO_", colnames(df)))
      inh.idx.inv <- intersect(grep("INV_", colnames(df)),
                           grep("_MENDELIAN_", colnames(df)))
      inh.idx.cpx <- intersect(grep("CPX_", colnames(df)),
                               grep("_MENDELIAN_", colnames(df)))
      minQuals <- unlist(lapply(strsplit(colnames(df)[dn.idx.inv], split="_"), function(v){as.numeric(v[[2]])}))
      dn.k.tab <- df[, dn.idx.inv] + df[, dn.idx.cpx]
      inh.k.tab <- df[, inh.idx.inv] + df[, inh.idx.cpx]
      all.k.tab <- dn.k.tab + inh.k.tab
    }
    if(mode=="marginal"){
      dnr.tab <- dn.k.tab / all.k.tab
      dnr.v <- as.numeric(apply(dnr.tab, 2, mean, na.rm=T))
      names(dnr.v) <- minQuals
      return(dnr.v)
    }else if(mode=="cumulative"){
      #Upper
      dn.ksum.kept.tab <- as.data.frame(sapply(1:ncol(dn.k.tab), function(i){
        apply(as.data.frame(dn.k.tab[, i:ncol(dn.k.tab)]), 1, sum)
      }))
      inh.ksum.kept.tab <- as.data.frame(sapply(1:ncol(inh.k.tab), function(i){
        apply(as.data.frame(inh.k.tab[, i:ncol(inh.k.tab)]), 1, sum)
      }))
      all.ksum.kept.tab <- dn.ksum.kept.tab + inh.ksum.kept.tab
      dnr.sum.kept.tab <- dn.ksum.kept.tab / all.ksum.kept.tab
      dnr.kept.v <- as.numeric(apply(dnr.sum.kept.tab, 2, mean, na.rm=T))
      #Lower
      dn.ksum.filtered.tab <- as.data.frame(sapply(2:ncol(dn.k.tab), function(i){
        apply(as.data.frame(dn.k.tab[, 1:i]), 1, sum)
      }))
      dn.ksum.filtered.tab <- cbind("V0"=rep(0, times=nrow(dn.ksum.filtered.tab)),
            dn.ksum.filtered.tab)
      inh.ksum.filtered.tab <- as.data.frame(sapply(2:ncol(inh.k.tab), function(i){
        apply(as.data.frame(inh.k.tab[, 1:i]), 1, sum)
      }))
      inh.ksum.filtered.tab <- cbind("V0"=rep(0, times=nrow(inh.ksum.filtered.tab)),
                                     inh.ksum.filtered.tab)
      all.ksum.filtered.tab <- dn.ksum.filtered.tab + inh.ksum.filtered.tab
      dnr.sum.filtered.tab <- dn.ksum.filtered.tab / all.ksum.filtered.tab
      dnr.filtered.v <- as.numeric(apply(dnr.sum.filtered.tab, 2, mean, na.rm=T))
      return(list("kept"=dnr.kept.v, "filtered"=dnr.filtered.v))
    }else{
      stop(paste("Specified mode (", mode, ") is not recognized", sep=""))
    }
  })
  names(res) <- c("ALL", "DEL", "DUP", "INS", "INV_CPX")
  return(res)
}

#Generate jitter residuals for a vector of values based on their density
sina.jitter <- function(vals){
  d <- density(vals)
  dv <- approx(d$x, d$y, xout=vals)
  dv <- dv$y/max(dv$y)
  dv.j <- sapply(1:length(vals), function(i){
    jitter(0, amount=dv[i])
  })
  return(dv.j)
}
#Generate sina points to add to existing plot
sina.plot <- function(vals, y.at, color, width=0.1){
  j <- (width*sina.jitter(vals))+y.at
  points(x=vals, y=j, pch=21, cex=0.25, col=color, bg=adjustcolor(color, alpha=0.3))
  segments(x0=median(vals), x1=median(vals), 
           y0=y.at-width, y1=y.at+width, lwd=3)
  # points(x=median(vals), y=y.at, pch=18)
}
#Plot jitter of de novo rates by SVTYPE
plot.denovorates <- function(dat, sv.groups, sv.groups.labels, filters, filter.labels, filter.cols){
  #Prep plot area
  dnrs <- dat[, intersect(grep("_HET_DNR", colnames(dat), fixed=T), 
                         grep("_HET_DNR_SV_COUNT", colnames(dat), fixed=T, invert=T))]
  counts <- dat[, grep("_HET_DNR_SV_COUNT", colnames(dat), fixed=T)]
  plot(x=c(0, max(dnrs, na.rm=T)), y=c(-0.1, -3.5), type="n", 
       xaxt="n", xlab="", yaxt="n", ylab="")
  # abline(v=0.05, col="gray80", lty=2)
  abline(v=0)
  x.ticks <- seq(0, max(axTicks(1)), by=0.05)
  axis(3, at=x.ticks, labels=NA)
  sapply(x.ticks, function(x){
    axis(3, at=x, tick=F, line=-0.5, cex.axis=0.8, 
         labels=paste(round(100*x, 1), "%", sep=""))
  })
  mtext(3, line=1.35, text=expression("Het." ~ italic("De Novo") ~"Rate"))
  #Add sina plots
  sapply(1:3, function(i){
    sapply(2:3, function(k){
      rates <- dnrs[, which(colnames(dnrs)==paste(sv.groups[i], filters[k], "HET_DNR", sep="_"))]
      totals <- counts[, which(colnames(counts)==paste(sv.groups[i], filters[k], "HET_DNR_SV_COUNT", sep="_"))]
      #Sina plot
      sina.plot(vals=rates, 
                color=filter.cols[k], 
                y.at=(-i+1)-(1.075-((4:1)[k]/4)))
      #Count & median DNR label
      med.count <- median(totals)
      med.dnr <- median(rates)
      axis(2, at=(-i+1)-(1.075-((4:1)[k]/4)), tick=F, line=-0.8, las=2, cex.axis=0.8, col.axis=filter.cols[k], 
           labels=paste(prettyNum(round(med.count, 0), big.mark=","), " / Child (", round(100*med.dnr, 1), "%)", sep=""))
    })
    #Group label
    axis(2, at=-i+0.95, tick=F, line=-0.8, las=2, cex.axis=0.9, font=2, 
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
  axis(2, at=(-4+1)-(1.075-((4:1)[2]/4)), tick=F, line=-0.8, las=2, cex.axis=0.8, col.axis=filter.cols[3], 
       labels=paste(prettyNum(round(med.count, 0), big.mark=","), " / Child (", round(100*med.dnr, 1), "%)", sep=""))
  axis(2, at=-4+0.95, tick=F, line=-0.8, las=2, cex.axis=0.9, font=2, 
       labels=sv.groups.labels[4])
}

# Plot de novo rate as a function of QUAL, for a single SV type
plot.denovoVsQual <- function(dat.byQual, svtype=NULL, qual.bin.step=50,
                           color="black", mode="cumulative"){
  # Get plot data
  qual.bins <- seq(0, 1000, qual.bin.step)
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
  ymax <- max(vals, na.rm=T)
  # Prep plot area
  par(mar=c(2.25, 4.1, 0.5, 0.5))
  plot(x=c(0, 1000), y=c(ymax+0.025, 0), ylim=c(ymax+0.025, 0), type="n",
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
  y.lab <- bquote("Het." ~ italic("De Novo") ~ "Rate")
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
    axis(2, at=y, tick=F, las=2, line=-0.7, cex.axis=0.85, 
         labels=paste(100*y, "%", sep=""))
  })
  mtext(2, text=y.lab, line=2)
  if(mode=="marginal"){
    mtext(2, text="Marginal", line=3)
  }
  box(bty="o")
}


################
###RSCRIPT BLOCK
################
require(optparse, quietly=T)
require(zoo, quietly=T)
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
args <- parse_args(OptionParser(usage="%prog DENOVO_DATA PCRMINUS_SAMPLES DENOVO_BY_QUAL OUTDIR", 
                                option_list=option_list), 
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
dat.in <- args$args[1]
PCRMINUS.samples.list.in <- args$args[2]
byqual.in <- args$args[3]
OUTDIR <- args$args[4]
prefix <- opts$prefix
svtypes.file <- opts$svtypes

# #Dev parameters (local)
# dat.in <- "~/scratch/gnomAD-SV_v2_rev1.PCRMINUS_trios_deNovoStats.perSample_summary.txt.gz"
# PCRMINUS.samples.list.in <- "~/scratch/gnomAD-SV_v2_rev1_wRelateds.PCRMINUS.samples.list"
# byqual.in <- "~/scratch/gnomAD-SV_v2_rev1.PCRMINUS_trios_deNovoStats_binnedByQual.perSample_summary.txt.gz"
# # byqual.in <- "~/scratch/dnv_by_qual_gs/test_dn_byQual.txt"
# OUTDIR <- "~/scratch/gnomAD_v2_test_plots/"
# prefix <- "gnomAD-SV_v2_rev1"
# svtypes.file <- "~/Desktop/Collins/Talkowski/code/sv-pipeline/ref/vcf_qc_refs/SV_colors.txt"

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

###Set analysis parameters
filters <- c("ANY", "PASS", "FAIL")
filter.labels <- c("All SVs", "PASS-only", "Non-PASS")
sv.groups <- c("ALL", "CNV", "BCA", "BND")
sv.groups.labels <- c("All SVs", "DEL/DUP", "INS/INV/CPX/CTX", "BND")
# filter.cols <- c("any"="gray60", 
#                 "pass"="#4DAC26", 
#                 "fail"="#D01C8B")
filter.cols <- c("ANY"="gray60", 
                 "PASS"="#4DAC26", 
                 "FAIL"="#AC26A1")
# svtypes.forAnalysis <- c("ALL", "DEL", "DUP", "INS", "INV", "CPX", "BND")
# svtype.colors <- c("gray50", 
#                    svtypes$color[which(svtypes$svtype=="DEL")], 
#                    svtypes$color[which(svtypes$svtype=="DUP")], 
#                    svtypes$color[which(svtypes$svtype=="INS")], 
#                    svtypes$color[which(svtypes$svtype=="INV")], 
#                    svtypes$color[which(svtypes$svtype=="CPX")])


###Process input data
cat("NOW LOADING DATA\n")
PCRMINUS.samples <- as.character(read.table(PCRMINUS.samples.list.in, header=F)[, 1])
dat <- importDeNovoRates(dat.in, svtypes.forAnalysis, PCRMINUS.samples)
dat.byQual.marginal <- import.denovoByQual(byqual.in, mode="marginal")
dat.byQual.cumulative <- import.denovoByQual(byqual.in, mode="cumulative")

###Plot de novo rates
png(paste(OUTDIR, "/", prefix, ".het_deNovo_rates.png", sep=""), 
    width=6*300, height=5*300, res=400)
par(mar=c(0.25, 7, 3, 0.25), bty="n")
plot.denovorates(dat, sv.groups, sv.groups.labels, filters, filter.labels, filter.cols)
dev.off()
png(paste(OUTDIR, "/", prefix, ".het_deNovo_rates.narrow.png", sep=""), 
    width=4.7*300, height=4.4*300, res=400)
par(mar=c(0.25, 6.5, 2.25, 0.25), bty="n")
plot.denovorates(dat, sv.groups, sv.groups.labels, filters, filter.labels, filter.cols)
dev.off()

###Print sums of variants per trio (for supplementary QC table)
apply(dat[, intersect(grep("ALL_PASS", colnames(dat), fixed=T), 
               grep("COUNT", colnames(dat), fixed=T))], 
      2, median)
apply(dat[, tail(intersect(grep("ALL_PASS", colnames(dat), fixed=T), 
                     grep("COUNT", colnames(dat), fixed=T, invert=T)), 4)], 
      2, median)


### Plot de novo rates vs QUAL scores
pdf(paste(OUTDIR, "/", prefix, ".deNovoRate_byQUAL.all.cumulative.pdf", sep=""),
    height=2, width=2.25)
plot.denovoVsQual(dat.byQual.cumulative,  svtype="ALL", mode="cumulative")
dev.off()
pdf(paste(OUTDIR, "/", prefix, ".deNovoRate_byQUAL.all.marginal.pdf", sep=""),
    height=2, width=2.25)
plot.denovoVsQual(dat.byQual.marginal,  svtype="ALL", mode="marginal")
dev.off()
sapply(c("DEL", "DUP", "INS", "INV_CPX"), function(svtype){
  pdf(paste(OUTDIR, "/", prefix, ".deNovoRate_byQUAL.", svtype, ".cumulative.pdf", sep=""),
      height=2, width=2.25)
  plot.denovoVsQual(dat.byQual.cumulative,  svtype=svtype, mode="cumulative")
  dev.off()
  pdf(paste(OUTDIR, "/", prefix, ".deNovoRate_byQUAL.", svtype, ".marginal.pdf", sep=""),
      height=2, width=2.25)
  if(svtype=="INV_CPX"){
    plot.denovoVsQual(dat.byQual.marginal,  svtype=svtype, mode="marginal",
                      color=svtypes$color[which(svtypes$svtype=="CPX")])
  }else{
    plot.denovoVsQual(dat.byQual.marginal,  svtype=svtype, mode="marginal",
                      color=svtypes$color[which(svtypes$svtype==svtype)])
  }
  dev.off()
})


  