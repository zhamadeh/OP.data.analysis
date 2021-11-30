collectBreaksAllFiles <- function(datapath="DATA/rdata/"){
  files <- datapath
  
  
  breaks.all.files <- list()
  breaksConfInt.all.files <- list()
  summaryBreaks <- list()
  n=1
  for (file in files) {
    message("Reading ... " , basename(file), " ... ",round(  (n/length(files))*100  ,  digits = 1  ) , "%"  )
    n=n+1
    data <- get(load(file))[c('breaks', 'confint','ID')]
    data$breaks$ID <- data$ID
    summaryBreaks[[basename(file)]] <- breakpointR::summarizeBreaks(data)
    breakpoints <- data$breaks
    breaks.confint <- data$confint
    if (length(breakpoints)) {
      suppressWarnings( breaks.all.files[[file]] <- breakpoints ) #TODO check if this can be done without warnings
    }  
    if (length(breaks.confint)) {
      suppressWarnings( breaksConfInt.all.files[[file]] <- breaks.confint ) 
    }  
  }
  return(breaks.all.files)
}


#Use p-values to call hotspots for breakpoints

breakpointHotspotter <- function(breaks.all.files){
  
  #Take in all breakpoints and calcualte p-values
  
  names(breaks.all.files) <- NULL
  gr.list=breaks.all.files
  bw=1000000
  pval=1e-8
  names(gr.list) <- NULL
  gr <- do.call(c, gr.list)
  gr <- GenomicRanges::sort(gr)
  
  ## Iterate over chromosomes and calculate p-values
  pranges.list <- GenomicRanges::GRangesList()
  hotspots=data.frame()
  count=1
  
  df = data.frame(chr=c(),midpoint=c(),density=c(),pvalue=c(),null_midpoints=c(),null_density=c())
  
  
  for (chrom in seqlevels(gr)) {
    grc <- gr[seqnames(gr)==chrom]
    if (length(grc)>1) {
      midpoints <- (start(grc)+end(grc))/2
      kde <- stats::density(midpoints,bw=bw,kernel='gaussian')
      # Random distribution of genomic events
      kde.densities <- numeric()
      
      for (i1 in seq_len(100)) {
        midpoints.r <- round(stats::runif(length(midpoints),1,seqlengths(gr)[chrom]))
        kde.r <- stats::density(midpoints.r,bw=bw,kernel='gaussian')
        kde.densities <- c(kde.densities, kde.r$y)
      }
      # Use ecdf to calculate p-values 
      p <- 1-stats::ecdf(kde.densities)(kde$y)
      pvalues <- data.frame(chromosome=chrom,start=kde$x,pvalue=p)
      # Make GRanges
      pvalues$end <- pvalues$start
      pvalues$chromosome <- factor(pvalues$chromosome, levels=seqlevels(gr))
      pvalues <- as(pvalues,'GRanges')
      seqlevels(pvalues) <- seqlevels(gr)
      suppressWarnings(
        seqlengths(pvalues) <- seqlengths(gr)[names(seqlengths(pvalues))]
      )
      # Resize from pointsize to bandwidth
      suppressWarnings(
        pvalues <- GenomicRanges::resize(pvalues, width=bw, fix='center')
      )
      pvalues <- trim(pvalues)
      ## Find regions where p-value is below specification
      mask <- pvalues$pvalue <= pval
      rle.pvals <- rle(mask)
      rle.pvals$values <- cumsum(rle.pvals$values+1)
      pvalues$group <- inverse.rle(rle.pvals)
      if (length(which(mask))>0) {
        pvalues.split <- split(pvalues[mask],pvalues$group[mask])
        pranges <- unlist(endoapply(pvalues.split, function(x) { y <- x[1]; end(y) <- end(x)[length(x)]; y$pvalue <- min(x$pvalue); return(y) }))
        pranges$group <- NULL
        pranges$num.events <- GenomicRanges::countOverlaps(pranges,grc)
        pranges.list[[chrom]] <- pranges
      }
      for (el in 1:length(pranges)){
        tmp = as.data.frame(gr[queryHits(findOverlaps(gr,pranges[el],type = "any"))])
        tmp$count=count
        hotspots <- rbind(tmp,hotspots)
        count=count+1
      }
      
      
      
      mp = kde$x
      ds = kde$y
      pv = pvalues$pvalue
      exp = data.frame(mp=mp,ds=ds)
      exp$type="BREAKPOINTS"
      exp$p = pv
      
      null_mp = kde.r$x
      null_ds = kde.r$y
      null = data.frame(mp=null_mp,ds=null_ds)
      null$type = "NULL"
      null$p = pv
      
      tmp = rbind(exp,null)
      tmp$chr = chrom
      
      #tmp2 = data.frame(pvalues=pv)
      #tmp2$chr = chrom
      
      df <- rbind(tmp,df)
      #df_p <- rbind(tmp2,df_p)
    }
    count=count+1
  }
  pranges <- unlist(pranges.list, use.names=FALSE)
  names(pranges) <- NULL
  
  write.table(df,"INPUT/07.densityPvalueSummary.txt",col.names = T,row.names = F, quote = F,sep="\t")
  
  return(hotspots)
}




 

breaks.all.files = collectBreaksAllFiles(list.files("../NA878/data_filtered/",full.names = T))
hotspots = breakpointHotspotter(breaks.all.files)


hotspots$count<- as.factor(hotspots$count)


# Create file structure for each inversion: 
## reads.bed files for all libraries involved
## breakpoints.bed for breakpoint coordinates
## chr_breakpoints.pdf for chromosome-specific ideogram plotting
## genotype.txt for summary of genotype results

# One summary file for all inversions, each row:
## chr, mean(start), mean(end), mean(width), # of libraries, % BLM, % RECQL5, % BLM-RECQL5, % WT 

#breakpointR::plotBreakpointsPerChr(files2plot = "../StraVa/DATA/rdata/",chromosomes = "chr12")

savingAndPrinting <- function(hotspots,hotpath="HOTSPOT_EVENTS",readspath="../StraVa/DATA/browserfiles/"){
  
  # Directory for creating file structure
  
  if (!file.exists(hotpath) ) { dir.create(hotpath)}
  
  # Count corresponds to a unique inversion so convert to factor for easier iterating
  hotspots$count=as.factor(hotspots$count)
  # Also easier for dealing with IDs as factor
  hotspots$ID <- as.factor(hotspots$ID)
  
  # Initialize empty dataframe for summary
  
  summary <- data.frame(chr=c(),start= c(),end=c(),count=c(), width=c(),n=c(),perc=c())
  
  
  numOfLibs = length(quality_libraries)
  
  # Iterate through each level of count (unique inversion) and 1) filter 2) summarize 3) save 4) plot
  hotspots$count<-as.factor(hotspots$count)
  
  level=levels(hotspots$count)[1]
  for (level in levels(hotspots$count)){
    message("Level: ",level)
    
    # 1) FILTER
    tmp <- dplyr::filter(hotspots, count==level)
    chr = as.character(tmp[1,1]) # Set chromosome
    
    tmp$ID <- droplevels(tmp$ID) # Drop unused levles
    
    
    datapath=paste0(hotpath,"/",chr,"-",level,"/")
    if (!file.exists(datapath) ) { dir.create(datapath)}
    readsdatapath=paste0(hotpath,"/",chr,"-",level,"/reads/")
    if (!file.exists(readsdatapath) ) { dir.create(readsdatapath)}
    
    numOfLibsInvolved=length(levels(droplevels(tmp$ID)))
    
    
  
    
    if (length(tmp$ID)>12){
      files2transfer=paste0(readspath,tmp$ID[1:12],"_reads.bed.gz")
    } else {
      files2transfer=paste0(readspath,tmp$ID,"_reads.bed.gz")
    }
    
    row <- data.frame(chr=chr,start= mean(tmp$start),end=mean(tmp$end),count=tmp$count[1], width=mean(tmp$width),n=numOfLibsInvolved)
    
    summary <- rbind(summary,row)
    
    file.copy(files2transfer,readsdatapath)
    files2plot=paste0("../StraVa/DATA/rdata/",levels(tmp$ID),".RData")
    
    export(select(tmp,c(seqnames,start,end)),paste0(datapath,"breakpoints.bed"),format = "bed")

    plotBreakpointsPerChr(files2plot,plotspath = datapath,chromosomes = c(chr))
    
  }
  return(summary)
}
  
summary <- savingAndPrinting(hotspots,hotpath="HOTSPOT_EVENTS/")
  
summaryHot = data.frame()

for (level in levels(hotspots$count)){
  tmp <- dplyr::filter(hotspots, count==level)
  row = data.frame(chr = tmp$seqnames[1],start=min(tmp$start),end=max(tmp$end))
  summaryHot <- rbind(summaryHot,row)
}
write.table(summaryHot,"hotspotsSummaryNA878.bed",quote = F,sep = "\t",row.names = F,col.names = F)




list.dirs("HOTSPOT_EVENTS/",recursive = F)

for (folder in list.dirs("HOTSPOT_EVENTS/",recursive = F)){
  files = list.files(paste0(folder,"/reads/"))
}

