#! /usr/bin/env Rscript

#CNView: Code to plot normalized coverage for CNV visualization from WGS data

# Copyright (c) 2016 Ryan Collins <rcollins@chgr.mgh.harvard.edu>
# Distributed under terms of the MIT license.

#Main plotting function, called later by Rscript (see bottom of script)
CNView <- function(chr,start,end,            #region to be plotted
                   sampleID,                 #Character vector of IDs of samples to plot
                   covmatrix,                #Absolute path of coverage matrix. Header with sample IDs required
                   compression="optimize",   #compression factor for rebinning, if desired
                   highlight=NA,             #list of coordinate pairs; intervals to highlight; defaults to query interval; NULL disables
                   highlightcol="gold",      #vector of colors to shade each highlighted interval
                   window=0,                 #distance to append to both sides of input interval for viewing
                   yscale="optimize",        #vector of values to be represented on y axis
                   normDist=5000000,        #distance outside region to normalize (both sides). Must either be int or "genome"
                   UCSCtracks=c("Gene",      #append UCSC sequence context information; choose either NULL
                                "SegDup",    #or up to three among "Gene", "Gap", "RepMask", "blacklist", and "SegDup"
                                "Gap"),
                   title=NULL,               #option to add custom title. Overrides default
                   legend=T,                 #logical option to plot legend
                   output=NULL,              #path to output as pdf. If NULL, will plot to active device
                   plot=TRUE,                #logical option to disable plot step; mandatory for output!=NULL
                   returnData=FALSE,         #logical option to return df of all normalized coverage values
                   quiet=FALSE){             #logical option to disable verbose output
  
  ##Sanity check input##
  if(!(is.numeric(c(start,end)) & end > start)){
    stop("INPUT ERROR: Improper input coordinates")}
  if(!(is.character(sampleID))){
    stop("INPUT ERROR: Improper sampleID")}
  if(compression!="optimize"){
    if(!(compression >= 1 & all.equal(compression,as.integer(compression)))){
      stop('INPUT ERROR: compression parameter must either be "optomize" or positive integer > 1')}}
  suppressWarnings(if(!(is.na(highlight)) & !(is.null(highlight))){
    if(length(highlightcol)!=length(highlight)){
      stop("INPUT ERROR: highlightcol must be same length as intervals to highlight")}})
  if(normDist!="genome"){
    if(window>normDist | !(all.equal(window,as.integer(window)))){
      stop('INPUT ERROR: window must be a positive, whole number less than normDist')}
  }
  if(!(is.null(UCSCtracks)) & length(as.vector(UCSCtracks)) > 3){
    stop('INPUT ERROR: UCSCtracks must be a vector of no more than three track names (see documentation for options)')}
  if(plot==F && !(is.null(output))){
    stop('INPUT ERROR: plot must be TRUE for output other than NULL')}
  if(file.exists(covmatrix)==F){
    stop('INPUT ERROR: coverage matrix file not found')
  }
  
  ##Set preferences##
  options(scipen=1000, #disables scientific notation
          warn=-1) #disables warnings
  
  ##Counts number of samples to plot##
  nsamp=length(sampleID)
  
  ##Interleaver Helper FX##
  interleave <- function(v1,v2){
    ord1 <- 2*(1:length(v1))-1
    ord2 <- 2*(1:length(v2))
    c(v1,v2)[order(c(ord1,ord2))]
  }
  
  ##Loads required packages; installs if necessary##
  if("RMySQL" %in% rownames(installed.packages()) == FALSE)
  {warning("RMySQL package not installed.  Attempting to install from CRAN...")
   install.packages("RMySQL",repos="http://cran.rstudio.com/")}
  suppressPackageStartupMessages(library(RMySQL))
  if("plyr" %in% rownames(installed.packages()) == FALSE)
  {warning("plyr package not installed.  Attempting to install from CRAN...")
   install.packages("plyr",repos="http://cran.rstudio.com/")}
  suppressPackageStartupMessages(library(plyr))
  if("MASS" %in% rownames(installed.packages()) == FALSE)
  {warning("MASS package not installed.  Attempting to install from CRAN...")
   install.packages("MASS",repos="http://cran.rstudio.com")}
  suppressPackageStartupMessages(library(MASS))
  
  ##Connects to UCSC, if necessary##
  if(!(is.null(UCSCtracks))){
    if(quiet==F){cat("Attempting to connect to UCSC Genome Browser...")}
    UCSC <- dbConnect(MySQL(),
                      user='genome',
                      dbname='hg19',
                      host='genome-mysql.cse.ucsc.edu')
    if(quiet==F){cat(" Complete\n")}
  }
  
  ##Parameter cleanup##
  if(!(is.null(highlight))){
    if(is.na(highlight)){
      highlight <- list(c(start,end))
      highlightcol <- "gold"
    }
  }
  if(!(is.null(UCSCtracks))){
    UCSCtracks <- rev(UCSCtracks)
  }
  
  ##Subset & Load Plotting Values##
  if(quiet==F){cat("Filtering & loading coverage matrix...")}
  if(normDist!="genome"){
    subcovmatrix <- tempfile()
    system(paste("head -n1 ",covmatrix," > ",subcovmatrix,sep=""))
    system(paste("awk -v OFS=\"\t\" '{ if ($1==\"",chr,"\" && $2<=",end+normDist," && $3>=",start-normDist,") print $0 }' ",covmatrix," >> ",
                 subcovmatrix,sep=""))
  }else{
    subcovmatrix <- covmatrix
  }
  cov <- read.table(subcovmatrix,header=T,sep="\t")
  if(quiet==F){cat(" Complete\n")}
  
  ##Rebin Helper FX##
  rebin <- function(df,compression){
    Chr <- df[1,1]
    Start <- df[1,2]
    End <- df[compression,3]
    for(i in 2:(floor(nrow(df)/compression))) {
      Chr <- c(Chr,as.character(df[((i-1)*compression)+1,1]))
      Start <- c(Start,as.integer(df[((i-1)*compression)+1,2]))
      End <- c(End,as.integer(df[i*compression,3]))
    }
    newvals <- apply(df[,4:ncol(df)],2,
                     function(vals,compression){
                       newcol <- sum(vals[1:compression])
                       for(i in 2:(floor(length(vals)/compression))) {
                         newcol <- c(newcol,as.integer(sum(vals[(((i-1)*compression)+1):(i*compression)])))
                       }
                       return(newcol)
                     },compression)
    return(as.data.frame(cbind(Chr,Start,End,newvals)))
  }
  
  ##Rebins values##
  if(quiet==F){cat("Compressing coverage matrix [")}
  obinsize <- cov[1,3]-cov[1,2]
  if(compression=="optimize"){
    compression <- round((end-start+(2*window))/120000)
  }
  binsize <- compression*obinsize
  if(quiet==F){cat(paste(prettyNum(binsize,big.mark=",")," bp bins]... ",sep=""))}
  if(compression>1){
    res <- rebin(cov,compression)    
  } else {
    res <- cov
  }
  if(quiet==F){cat(" Complete\n")}
  
  ##Scale each col within that sample by median##
  if(quiet==F){cat("Performing intra-sample normalization...")}
  res[,4:ncol(res)] <- apply(res[,4:ncol(res)],2,
                             function(vals){
                               nvals <- as.numeric(vals)/median(as.numeric(vals))
                               nvals[is.infinite(nvals)] <- NA
                               return(nvals)
                             })
  if(quiet==F){cat(" Complete\n")}
  
  ##Normalize each row across all samples##
  if(quiet==F){cat("Performing inter-sample normalization...")}
  colnames(res)[1:3] <- c("Chr","Start","End")
  names <- colnames(res)
  oncol <- ncol(res)
  res[,4:oncol] <- data.frame(t(apply(res[,4:oncol],1,function(vals){
    return(scale(as.numeric(vals)))
  })))
  res[is.na(res)] <- 0
  res$mean <- apply(res[,4:oncol],1,mean)
  res$sd <- apply(res[,4:oncol],1,sd)
  res$median <- apply(res[,4:oncol],1,median)
  res$mad <- apply(res[,4:oncol],1,mad)
  if(quiet==F){cat(" Complete\n")}
  
  ##Subset View Window##
  plotSet <- as.data.frame(apply(res[which(as.integer(as.character(res$Start)) <= end+window & 
                                             as.integer(as.character(res$End)) >= start-window & 
                                             as.character(res$Chr) == as.character(chr)),],
                                 2,function(col){
                                   return(as.numeric(as.character(col)))
                                 }))
  plotSet[is.na(plotSet)] <- 0
  
  ##Get Sample Indexes##
  sampIdx <- as.vector(sapply(as.vector(sampleID),function(val){grep(val,colnames(plotSet),ignore.case=T)}))
  
  ##Output Options##  
  if(plot==T){
    if(!(is.null(output))){
      if(quiet==F){cat(paste("Plotting samples to ",output,"...",sep=""))}
      pdf(output,width=10.5,height=(5+(2.5*nsamp)))
    } else {
      if(quiet==F){cat("Plotting samples to screen...")}
    }
    par(mar=c(0,0,0,0),
        oma=c(5,4,4,2))
    if(!(is.null(UCSCtracks))){
      layout(matrix(c(1:(nsamp+1)),byrow=T),heights=c(rep(4,nsamp),1))
    } else {
      par(mfrow=c(nsamp,1))
    }
    
    ##Plot per sample##
    for(k in 1:nsamp){
      
      ##Generate Colors##
      if(pt(plotSet[1,sampIdx[k]],df=(ncol(plotSet)-8)) >= 1-(0.05/(ncol(plotSet)-7))){
        colval <- "blue"
      }else if(pt(plotSet[1,sampIdx[k]],df=(ncol(plotSet)-8)) <= 0.05/(ncol(plotSet)-7)){
        colval <- "red"
      }else{
        colval <- "gray40"
      }
      for(i in 2:nrow(plotSet)){
        if(pt(plotSet[i,sampIdx[k]],df=(ncol(plotSet)-8)) >= 1-(0.05/(ncol(plotSet)-7))){
          colval <- c(colval,"blue")
        } else if(pt(plotSet[i,sampIdx[k]],df=(ncol(plotSet)-8)) <= 0.05/(ncol(plotSet)-7)){
          colval <- c(colval,"red")
        } else {
          colval <- c(colval,"gray40")
        }
      }
      colval <- interleave(colval,colval)
      for(i in 2:length(colval)){
        if(colval[i]=="blue"){
          colval[i-1] <- "blue"
        }else if(colval[i]=="red"){
          colval[i-1] <- "red"
        }
      }
      
      ####Plot main####
      plot(as.numeric(plotSet$Start),
           as.numeric(plotSet[,sampIdx[k]]),
           type="n",xlim=c(max(start-window,0),
                           max(plotSet$Start)),
           if(yscale=="optimize"){
             ylim=c(min(0,min(plotSet[,sampIdx[k]]))-1,
                    max(0,max(plotSet[,sampIdx[k]]))+1)
           }else{
             ylim=yscale
           },
           xlab="",xaxt="n",ylab="",xaxs="i",yaxt="n",
           panel.first=c(abline(h=seq(round_any(par("usr")[3],2),
                                      round_any(par("usr")[4],2),by=2),
                                col="gray80"),
                         polygon(y=c(2*plotSet$mad,
                                     rev(-2*plotSet$mad)),
                                 x=c(as.numeric(as.character(plotSet$Start)),
                                     rev(as.numeric(as.character(plotSet$Start)))),
                                 col="gray85",border="gray60"),
                         polygon(y=c(plotSet$mad,
                                     rev(-plotSet$mad)),
                                 x=c(as.numeric(as.character(plotSet$Start)),
                                     rev(as.numeric(as.character(plotSet$Start)))),
                                 col="gray75",border="gray55"),
                         points(as.integer(as.character(plotSet$Start)),
                                as.numeric(as.character(plotSet$median)),
                                lty=2,type="l"),
                         if(!(is.null(highlight))){
                           for(i in 1:length(highlight)){
                             rect(highlight[[i]][1],
                                  par("usr")[3],
                                  highlight[[i]][2],
                                  par("usr")[4],
                                  col=adjustcolor(highlightcol[i],alpha=0.2),
                                  border=NA)
                             abline(v=highlight[[i]][1],
                                    lty=3,lwd=2)
                             abline(v=highlight[[i]][2],
                                    lty=3,lwd=2)
                           }
                         },
                         segments(x0=interleave(plotSet[seq(1,(nrow(plotSet)-1)),2],
                                                plotSet[seq(2,(nrow(plotSet))),2]),
                                  y0=interleave(plotSet[seq(1,(nrow(plotSet)-1)),sampIdx[k]],
                                                plotSet[seq(1,(nrow(plotSet)-1)),sampIdx[k]]),
                                  x1=interleave(plotSet[seq(2,(nrow(plotSet))),2],
                                                plotSet[seq(2,(nrow(plotSet))),2]),
                                  y1=interleave(plotSet[seq(1,(nrow(plotSet)-1)),sampIdx[k]],
                                                plotSet[seq(2,(nrow(plotSet))),sampIdx[k]]),
                                  lwd=3,
                                  col=colval)))
      #Y Axis
      axis(2,at=seq(round_any(par("usr")[3],2),
                    round_any(par("usr")[4],2),by=2),
           las=2)
      mtext(paste("Norm. Depth t Score",sep=""),
            side=2,line=2,cex=0.7)
      #Print Sample ID if >1 sample
      if(nsamp>1){
        text(x=mean(par("usr")[1:2]),y=par("usr")[4],labels=names(plotSet)[sampIdx[k]],cex=1.2,pos=1,font=2)
      }
      
      ##Legend & title if first sample
      if(k==1){
        #Title
        if(!(is.null(title))){
          mtext(text=title,outer=T,side=3,line=1,font=2,cex=1.6)
        }else{
          if(nsamp==1){
            mtext(text=paste("Normalized Sequencing Depth of ",sampleID,sep=""),
                  outer=T,side=3,line=1,font=2,cex=1.6)
          }else{
            mtext(text=paste("Normalized Sequencing Depth of ",nsamp," Samples",sep=""),
                  outer=T,side=3,line=1,font=2,cex=1.6)
          }
        }
        mtext(text=paste("chr",chr," : ",prettyNum(max((start-window),0),big.mark=",")," - ",
                         prettyNum(end+window,big.mark=","),sep=""),
              outer=T,side=3,line=0)
        #Legend
        if(legend==T){
          if(max(plotSet[,sampIdx])+min(plotSet[,sampIdx]) >= 0){
            legend("topright",
                   legend=c(paste("p(Dup) < Bonferroni (df=",ncol(plotSet)-8,")",sep=""),
                            paste("p(Del) < Bonferroni (df=",ncol(plotSet)-8,")",sep=""),
                            "Median t Score",
                            "1 * MAD",
                            "2 * MAD"),
                   pch=c(NA,NA,NA,15,15),pt.cex=c(1,1,1,1.5,1.5),
                   lty=c(1,1,2,NA,NA),lwd=c(4,4,1,NA,NA),
                   col=c("blue","red","black","gray54","lightgray"),
                   bg="white",cex=0.8)
            text(x=par("usr")[1],
                 y=0.95*par("usr")[4],
                 labels=paste(prettyNum(binsize,big.mark=",")," bp Bins",sep=""),
                 font=4,pos=4)
          } else if(max(plotSet[,sampIdx])+min(plotSet[,sampIdx]) < 0){
            legend("bottomright",
                   legend=c(paste("p(Dup) < Bonferroni (df=",ncol(plotSet)-8,")",sep=""),
                            paste("p(Del) < Bonferroni (df=",ncol(plotSet)-8,")",sep=""),
                            "Median t Score",
                            "1 * MAD",
                            "2 * MAD"),
                   pch=c(NA,NA,NA,15,15),pt.cex=c(1,1,1,1.5,1.5),
                   lty=c(1,1,2,NA,NA),lwd=c(4,4,1,NA,NA),
                   col=c("blue","red","black","gray54","lightgray"),
                   bg="white",cex=0.7)
            text(x=par("usr")[1],
                 y=0.95*par("usr")[3],
                 labels=paste(prettyNum(binsize,big.mark=",")," bp Bins",sep=""),
                 font=4,pos=4)
          }
        }
      }
    }
    if(quiet==F){cat(" Complete\n")}
    
    ##UCSC plot##
    if(!(is.null(UCSCtracks))){
      if(quiet==F){cat("Appending UCSC tracks...")}
      plot(plotSet$Start,
           c(1,2,rep(3,nrow(plotSet)-2)),
           ylim=c(0,3),
           type="n",xaxt="n",yaxt="n",
           ylab="",xlab="",xaxs="i")
      axis(1,at=seq(min(plotSet$Start),
                    max(plotSet$Start),
                    by=(max(plotSet$Start)-min(plotSet$Start))/8),
           labels=prettyNum(seq(min(plotSet$Start),
                                max(plotSet$Start),
                                by=(max(plotSet$Start)-min(plotSet$Start))/8),
                            big.mark=","),
           cex.axis=0.9)
      if("Gap" %in% UCSCtracks){
        gaps <- dbGetQuery(UCSC,paste("SELECT chromStart, chromEnd, type FROM gap WHERE `chrom` = 'chr",chr,"' ",
                                      "AND `chromStart` <= ",end+window," ",
                                      "AND `chromEnd` >= ",start-window,sep=""))
        if(nrow(gaps) > 0){
          for(i in 1:nrow(gaps)){
            rect(xleft=gaps$chromStart[i],
                 ybottom=grep("Gap",UCSCtracks)-.9,
                 xright=gaps$chromEnd[i],
                 ytop=grep("Gap",UCSCtracks)-.1,
                 col="cadetblue3",border=NA)
          }
        }
      }
      if("SegDup" %in% UCSCtracks){
        segDups <- dbGetQuery(UCSC,paste("SELECT chromStart, chromEnd FROM genomicSuperDups WHERE `chrom` = 'chr",chr,"' ",
                                         "AND `chromStart` <= ",end+window," ",
                                         "AND `chromEnd` >= ",start-window,sep=""))
        if(nrow(segDups) > 0){
          for(i in 1:nrow(segDups)){
            rect(xleft=segDups$chromStart[i],
                 ybottom=grep("SegDup",UCSCtracks)-.9,
                 xright=segDups$chromEnd[i],
                 ytop=grep("SegDup",UCSCtracks)-.1,
                 border=NA,col="darkorange")
          }
        }
      }
      if("RepMask" %in% UCSCtracks){
        repeatMasker <- dbGetQuery(UCSC,paste("SELECT genoStart, genoEnd, repClass FROM rmsk WHERE `genoName` = 'chr",chr,"' ",
                                              "AND `genoStart` <= ",end+window," ",
                                              "AND `genoEnd` >= ",start-window,sep=""))
        if(nrow(repeatMasker) > 0){
          for(i in 1:nrow(repeatMasker)){
            rect(xleft=repeatMasker$genoStart[i],
                 ybottom=grep("RepMask",UCSCtracks)-.9,
                 xright=repeatMasker$genoEnd[i],
                 ytop=grep("RepMask",UCSCtracks)-.1,
                 border=NA,col="lightsteelblue4")
          }
        }
      }
      if("Gene" %in% UCSCtracks){
        genes <- unique(dbGetQuery(UCSC,paste("SELECT txStart, txEnd, name2, strand, exonStarts, exonEnds FROM refGene WHERE `chrom` = 'chr",chr,"' ",
                                              "AND `txStart` <= ",end+window," ",
                                              "AND `txEnd` >= ",start-window,sep="")))
        if(nrow(genes) > 0){
          for(i in 1:nrow(genes)){
            rect(xleft=genes$txStart[i],
                 ybottom=grep("Gene",UCSCtracks)-.6,
                 xright=genes$txEnd[i],
                 ytop=grep("Gene",UCSCtracks)-.1,
                 border=NA,col="lightgreen")
            segments(x0=genes$txStart[i],
                     x1=genes$txEnd[i],
                     y0=grep("Gene",UCSCtracks)-0.35,
                     y1=grep("Gene",UCSCtracks)-0.35,
                     col="darkgreen")
          }
          exons <- unique(data.frame("start"=as.character(unlist(sapply(genes$exonStarts,
                                                                        function(string){return(strsplit(string,split=","))}))),
                                     "end"=as.character(unlist(sapply(genes$exonEnds,
                                                                      function(string){return(strsplit(string,split=","))})))))
          for(i in 1:nrow(exons)){
            rect(xleft=exons[i,1],
                 xright=exons[i,2],
                 ybottom=grep("Gene",UCSCtracks)-.55,
                 ytop=grep("Gene",UCSCtracks)-.15,
                 border=NA,col="darkgreen")
          }
          for(i in unique(genes$name2)){
            text(x=(max(c(min(genes[which(genes$name2==i),1]),par("usr")[1]))+min(c(max(genes[which(genes$name2==i),2]),par("usr")[2])))/2,
                 y=grep("Gene",UCSCtracks)-.8,
                 labels=i,
                 cex=0.75,font=4)
          }
        }
      }
      axis(2,at=c(seq(0.5,length(UCSCtracks)-0.5)),
           labels=UCSCtracks,
           cex.axis=0.6,las=1,tick=F)
      if(quiet==F){cat(" Complete\n")}
    }
    
    ##Adds X axis##
    axis(1,at=seq(min(plotSet$Start),
                  max(plotSet$Start),
                  by=(max(plotSet$Start)-min(plotSet$Start))/8),
         labels=prettyNum(seq(min(plotSet$Start),
                              max(plotSet$Start),
                              by=(max(plotSet$Start)-min(plotSet$Start))/8),
                          big.mark=","),
         cex.axis=0.9)
    mtext(paste("chr",chr," Coordinate (bp)",sep=""),
          side=1,outer=T,line=2)
    
    ##Close Output##
    if(!(is.null(output))){
      dev.off()
    }
  }
  
  ##Disconnect from UCSC, if necessary##
  if(exists("UCSC")){
    dbDisconnect(UCSC)
  }
  
  ##Return Values as df if specified#
  if(returnData==T){
    if(quiet==F){cat("Returning normalized values\n")}
    if(quiet==F){cat(paste("\n** FINISHED ON ",date()," **\n\n",sep=""))}
    return(plotSet)
  }
  
  #Remove temporary file##
  if(normDist!="genome"){
    system(paste("rm ",subcovmatrix,sep=""))
  }
  
  ##Finish up##
  if(quiet==F){cat(paste("\n** FINISHED ON ",date()," **\n\n",sep=""))}
}

####################################################
####Rscript functionality for command line usage####
####################################################

#Disables factor default
options(stringsAsFactors=F)

#installs optparse if not already installed
if("optparse" %in% rownames(installed.packages()) == FALSE)
{install.packages("optparse",repos="http://cran.rstudio.com")}
suppressWarnings(suppressPackageStartupMessages(library(optparse)))

#list of Rscript options
option_list <- list(
  make_option(c("-c", "--compression"), type="integer", default=NULL,
              help="compression scalar for rebinning, if desired [default '%default']", 
              metavar="integer"),
  make_option(c("-i","--highlight"), type="character", default=NA,
              help="tab-delimited list of coordinate pairs for intervals to highlight and color as third column; NULL disables highlighting [default %default]",
              metavar="character"),
  make_option(c("-w","--window"), type="integer", default=NULL,
              help="distance to append to both sides of input interval for viewing [default 60% of plot interval]",
              metavar="integer"),
  make_option(c("--ymin"), type="integer", default=NULL,
              help="minimum value for y axis [default %default]", 
              metavar="integer"),
  make_option(c("--ymax"), type="integer", default=NULL,
              help="maximum value for y axis [default %default]", 
              metavar="integer"),
  make_option(c("-n","--normDist"), type="integer", default=5000000,
              help="distance outside region to use for normalization (both sides) [default %default]",
              metavar="integer"),
  make_option(c("-t","--title"), type="character", default=NULL,
              help="custom title for plot [default %default]",
              metavar="character"),
  make_option(c("-q","--quiet"), action="store_true", default=FALSE,
              help="disable verbose output [default %default]"),
  make_option(c("-l", "--nolegend"), action="store_false", default=TRUE,
              help="disable legend on plot [default %default]")
)

#Get command-line arguments & options
parser <- OptionParser(usage="%prog [options] chr start end samples.list covmatrix.bed outfile",
                       option_list=option_list,add_help_option=T)
args <- parse_args(parser,positional_arguments=TRUE)
opts <- args$options

#checks for appropriate positional arguments
if(length(args$args) != 6) {
  print_help(parser)
  stop("Incorrect number of required positional arguments")
}

#Prints startup message
if(opts$quiet==F){cat("+-------------------+\n| CNView Visualizer |\n|     (c) 2016      |\n+-------------------+\n")}

#Processes arguments & cleans options
if(file.exists(args$args[4])==T){
  if(opts$quiet==F){cat(paste("Reading sample IDs from ",args$args[4],"\n",sep=""))}
  samps <- read.table(args$args[4])[,1]
}else{
  if(opts$quiet==F){cat(paste("Sample ID file '",args$args[4],"' not found, assuming single sample ID provided\n",sep=""))}
  samps <- args$args[4]
}
suppressWarnings(if(!(is.na(opts$highlight)) & !(is.null(opts$highlight))){
  hightable <- read.table(opts$highlight,sep="\t",header=F)[,1:3]
  highopt <- as.list(as.data.frame(t(hightable[,1:2])))
  highcolopt <- as.vector(hightable[,3])
}else{
  highopt <- opts$highlight
  highcolopt <- "gold"
})
if(!(is.null(opts$ymin)) & !(is.null(opts$ymax))){
  yscale=c(opts$ymin,opts$ymax)
}else{
  yscale="optimize"
}
if(is.null(opts$compression)){
  compression <- "optimize"
}else{
  compression <- opts$compression
}
if(is.null(opts$window)){
  window=round(0.6*(as.numeric(args$args[3])-as.numeric(args$args[2])),0)
}else{
  window=round(opts$window,0)
}

#Runs CNView function
if(opts$quiet==T){
  sink("/dev/null")
}
CNView(chr=args$args[1],
       start=as.numeric(args$args[2]),
       end=as.numeric(args$args[3]),
       sampleID=as.vector(samps),
       covmatrix=args$args[5],
       compression=compression,
       highlight=highopt,
       highlightcol=highcolopt,
       window=window,
       yscale=yscale,
       normDist=opts$normDist,
       UCSCtracks=c("Gene","SegDup","Gap"),
       title=opts$title,
       output=args$args[6],
       plot=T,
       returnData=F,
       quiet=opts$quiet)