#! /usr/bin/env Rscript

#CNView-Multicopy: a visualization and annotation tool for multiallelic copy number variation from whole-genome sequencing

# Copyright (c) 2016 Ryan Collins <rcollins@chgr.mgh.harvard.edu>
# Distributed under terms of the MIT license.

# NOTE: Requires bedtools executable to be on your environment $PATH

#Main plotting function, called later by Rscript (see bottom of script)
CNView_multi <- function(chr,start,end,          #region to be plotted
                         covmatrix,              #Absolute path of coverage matrix. Header with sample IDs required
                         smoothing=0.05,         #Smoothing factor. [0,1]
                         window=NULL,            #window with which to pad interval; NULL defaults to ±20%
                         yscale=NULL,            #vector of copy states to be represented on y axis
                         normSamp=NULL,          #number of bins to randomly select for genome-wide normalization. NULL=number of bins in plot interval
                         tranches=16,            #number of tranches for copy-state histogram windows
                         UCSCtracks=c("Gene",    #append UCSC sequence context information; choose either NULL
                                      "SegDup",  #or up to three among "Gene", "Gap", "RepMask", "blacklist", and "SegDup"
                                      "Gap"),
                         genesymbols=TRUE,       #print gene symbols below UCSC gene body annotations
                         gcex=1,                 #global scaling for all fonts
                         title=NULL,             #option to add custom title. Overrides default
                         legend=T,               #logical option to plot legend
                         output=NULL,            #path to output as pdf. If NULL, will plot to active device
                         noUnix=FALSE,           #logical option to specify a non-unix OS (i.e. no awk, needed to read data)
                         quiet=FALSE){           #logical option to disable verbose output
  ##Sanity check input##
  if(!(is.numeric(c(start,end)) & end > start)){
    stop("INPUT ERROR: Improper input coordinates")}
  if(!(is.numeric(smoothing)) | smoothing<0 | smoothing>1){
    stop('INPUT ERROR: smoothing parameter must be float ~ [0,1]')}
  suppressWarnings(if(!(is.na(highlight)) & !(is.null(highlight))){
    if(length(highlightcol)!=length(highlight)){
      stop("INPUT ERROR: highlightcol must be same length as intervals to highlight")}})
  if(noUnix==T){
    warning('noUnix parameter specified as TRUE; operation speed will be substantially slower')
  }
  
  ##Set preferences##
  options(scipen=1000, #disables scientific notation
          warn=-1) #disables warnings
  
  ##Interleaver Helper FX##
  interleave <- function(v1,v2){
    ord1 <- 2*(1:length(v1))-1
    ord2 <- 2*(1:length(v2))
    c(v1,v2)[order(c(ord1,ord2))]
  }
  
  ##Random File Sampling FX Direct From File##
  ##Credit: Martin Morgan @
  #http://stackoverflow.com/questions/15532810/reading-40-gb-csv-file-into-r-using-bigmemory/18282037#18282037
  fsample <-function(fname, n, seed, header=FALSE, ..., reader = read.csv){
      set.seed(seed)
      con <- file(fname, open="r")
      hdr <- if (header) {
        readLines(con, 1L)
      } else character()
      buf <- readLines(con, n)
      n_tot <- length(buf)
      repeat {
        txt <- readLines(con, n)
        if ((n_txt <- length(txt)) == 0L)
          break
        n_tot <- n_tot + n_txt
        n_keep <- rbinom(1, n_txt, n_txt / n_tot)
        if (n_keep == 0L)
          next
        keep <- sample(n_txt, n_keep)
        drop <- sample(n, n_keep)
        buf[drop] <- txt[keep]
      }
      reader(textConnection(c(hdr, buf)), header=header, ...)
    }
  
  ##Loads required packages; installs if necessary##
  if("RMySQL" %in% rownames(installed.packages()) == FALSE)
  {warning("RMySQL package not installed.  Attempting to install from CRAN...")
    install.packages("RMySQL",repos="http://cran.rstudio.com/")}
  suppressPackageStartupMessages(library(RMySQL))
  if("colorspace" %in% rownames(installed.packages()) == FALSE)
  {warning("colorspace package not installed.  Attempting to install from CRAN...")
    install.packages("colorspace",repos="http://cran.rstudio.com/")}
  suppressPackageStartupMessages(library(colorspace))
  
  ##Parameter cleanup##
  if(is.null(window)){
    window <- round(0.2*(end-start),0)
  }
  if(!(is.null(UCSCtracks))){
    UCSCtracks <- rev(UCSCtracks)
  }
  
  ##Subset & Load Plotting Values##
  if(quiet==F){cat("Filtering & loading coverage matrix...")}
  if(noUnix==TRUE){
    cov <- read.table(covmatrix,header=T,sep="\t")
    cov <- cov[which(cov[,1]==chr & cov[,2]<=end & cov[,3]>=start),]
  }else{
    subcovmatrix <- tempfile()
    system(paste("head -n1 ",covmatrix," > ",subcovmatrix,sep=""))
    system(paste("awk -v OFS=\"\t\" '{ if ($1==\"",chr,"\" && $2<=",end+window," && $3>=",start-window,") print $0 }' ",covmatrix," >> ",
                 subcovmatrix,sep=""))
    cov <- read.table(subcovmatrix,header=T,sep="\t")
  }
  if(quiet==F){cat(" Complete\n")}
  
  ##Calculate number of bins to sample for normalization##
  if(is.null(normSamp)){
    normSamp <- nrow(cov)
  }
  
  ##Sample & Load Normalization Values
  if(quiet==F){cat("Sampling & loading normalization matrix...")}
  normcov <- fsample(covmatrix,normSamp,sample(1:1000,1),sep="\t",header=T)
  if(quiet==F){cat(" Complete\n")}
  
  ##Calculate Median Bin Coverage Per Library##
  covMeds <- apply(normcov[which(apply(normcov[,-c(1:3)],1,median)>0),
                           -c(1:3)],
                   2,median)
  
  ##Normalize Coverage Matrix by Library Medians##
  for(i in 4:ncol(cov)){
    cov[,i] <- 2*cov[,i]/covMeds[i-3]
  }
  
  ##Rewrite Bins as NA Coverage Where Third Quartile == 0##
  cov[apply(cov[,-c(1:3)],1,function(vals){
    return(quantile(vals,0.75))
  })==0,-c(1:3)] <- NA
  
  ##Subset Coverage Matrix Into N Tranches & Estimate Copy State##
  tranche_size <- ceiling(nrow(cov)/tranches)
  trancheCov <- as.data.frame(t(sapply(1:tranches,function(i){
    tCov <- cov[seq((tranche_size*(i-1))+1,tranche_size*i),]
    #Report Median Coverage Per Sample In Tranche#
    tCovNorm <- apply(tCov[,-c(1:3)],2,function(vals){return(median(vals,na.rm=T))})
    #Rough Estimate of Copy State by Rounding to Nearest Integer#
    tCovEst <- round(tCovNorm,0)
    if(any(!(is.na(tCov[,-c(1:3)])))){
      #K-Means Clustering Where K=Number of Unique Predicted Copy States#
      tCovClust <- kmeans(tCovNorm,centers=unique(tCovEst))
      #Calculate Cluster Adjustment#
      cAdj <- sapply(as.vector(tCovClust$centers),function(val){
        return(round(val,0)-val)
      })
      #Apply Cluster Adjustment#
      for(i in 1:length(cAdj)){
        tCovNorm[which(tCovClust$cluster==i)] <- tCovNorm[which(tCovClust$cluster==i)]+cAdj[i]
      }
      tCovEst <- tCovClust$cluster
      for(i in 1:length(cAdj)){
        tCovEst[which(tCovClust$cluster==i)] <- round(mean(tCovNorm[which(tCovClust$cluster==i)]),0)
      }
    }
    return(list(tCovNorm,tCovEst))
  })))
  cleanCov <- t(as.data.frame(trancheCov[[1]]))
  rownames(cleanCov) <- 1:tranches
  cleanCN <- t(as.data.frame(trancheCov[[2]]))
  rownames(cleanCN) <- 1:tranches

  ##Plot All Libraries -- DEV##
  par(mfrow=c(2,1))
  ymax=max(cleanCN,na.rm=T)+1
  CNcolors <- c("black",rainbow_hcl(ymax-1,c=75,l=70))
  plot(lowess(cov[,2],cov[,4],f=0.25*smoothing),type="l",ylim=c(0,ymax),
       col=adjustcolor("black",alpha=0.3),
       xaxs="i",yaxs="i",
       panel.first=c(rect(xleft=cov[seq(1,nrow(cov),by=tranche_size),2],
                          xright=cov[seq(tranche_size,nrow(cov)+tranche_size,by=tranche_size),2],
                          ybottom=par("usr")[3],
                          ytop=par("usr")[4],
                     col=c("gray95","gray87"),
                     border="gray50"),
                     abline(h=2,lwd=2),
                     abline(h=0:ymax,lty=2,col="gray50")))
  for(i in 5:ncol(cov)){
    points(lowess(cov[,2],cov[,i],f=0.25*smoothing),type="l",
           col=adjustcolor("black",alpha=0.3))
  }
  
  ##Plot All Tranche Coverages -- DEV##
  plot(x=c(1,tranches),
       y=c(0,ymax),
       type="n",
       xaxs="i",yaxs="i",
       panel.first=c(rect(xleft=0:(tranches-1),
                          xright=1:tranches,
                          ybottom=par("usr")[3],
                          ytop=par("usr")[4],
                          col=c("gray95","gray87"),
                          border="gray50"),
                     abline(h=2,lwd=2),
                     abline(h=0:ymax,lty=2,col="gray50")))
  for(i in 1:ncol(cleanCov)){
    points(x=1:tranches-(sample(30:70,tranches)/100),
           y=cleanCov[,i],
           pch=19,cex=0.7,
           col=adjustcolor(CNcolors[cleanCN[,i]+1],alpha=0.5))
  }
  
  ##Disconnect from UCSC, if necessary##
  if(exists("UCSC")){
    dbDisconnect(UCSC)
  }
  
  #Remove temporary file##
  if(exists(subcovmatrix)==T){
    if(file.exists(subcovmatrix)){
      system(paste("rm ",subcovmatrix,sep=""))
    }    
  }
  if(exists(subnormmatrix)==T){
    if(file.exists(subnormmatrix)){
      system(paste("rm ",subnormmatrix,sep=""))
    }
  }
  
  ##Finish up##
  if(quiet==F){cat(paste("\n** FINISHED ON ",date()," **\n\n",sep=""))}
}