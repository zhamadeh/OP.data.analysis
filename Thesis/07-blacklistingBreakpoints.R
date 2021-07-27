#######################################################
####         Workflow Part 7: Blacklisting          ###
#######################################################

frequencyFilterBreakpoints <- function(summaryBreaks.df, blacklist="Input/00.centromeres2.txt"){
  
  summaryBreaks<-read.table("Input/05.breakpoints.metrics.quality.gene.txt",header = T)
  summaryBreaks.df<- GRanges(summaryBreaks)
  
  centromeres <- read.table(blacklist,header=F) #%>% select(-c(V4))
  centromeres <-centromeres %>% dplyr::rename("seqnames"=V1,"start"=V2,"end"=V3)
  centroGRange <- GRanges(centromeres)
  
  suppressWarnings(breakGRanges <-summaryBreaks.df[-queryHits(findOverlaps(summaryBreaks.df, centroGRange, type="any")),])
  breakGRanges$freq = 0
  
  
  for (row in 1:length(breakGRanges)){
    tmp = breakGRanges[row,]
    breakGRanges[row,]$freq = countOverlaps(tmp,breakGRanges, maxgap = 1000000)
  }
  
  #breakGRanges$allele_freq=breakGRanges$freq/(length(levels(as.factor(summaryBreaks$library))))
  
  inversions <- breakGRanges[breakGRanges$freq>20,]
  breaks <- breakGRanges[breakGRanges$freq<=20,]
  breaks <- as.data.frame(breaks)
  #breaks <- filter(breaks,width<1000000)
  
  write.table(breaks,"Input/06.sisterChromatidExchanges.txt",row.names = F,col.names = T,quote = F,sep="\t")
  
  #breaks$filenames <- tools::file_path_sans_ext(breaks$filenames)
  breaks$library <- as.factor(breaks$library)
  levels(breaks$library) <- droplevels(breaks$library)
  
  counter=0
  for (f in levels(breaks$library)){
    counter=counter+1
    
    file = paste0("Output/bpr/data/",f,".RData")
    message("Working on ",basename(file), " ... ", (counter/length(levels(breaks$library)))*100,"%")
    tmp <- get(load(file))
    
    seqinfo <- tmp$breaks@seqinfo
    seqnameLevels <- levels(tmp$breaks@seqnames)
    
    
    tmp$breaks <-GRanges(filter(breaks,library==tmp$ID) %>% select(c(seqnames,start,end,width,strand,genoT,deltaW)))
    tmp$breaks@seqinfo <- seqinfo
    levels(tmp$breaks@seqnames) <- seqnameLevels
    
    tmp$confint
    tmp$confint<-GRanges(filter(breaks,library==tmp$ID) %>% select(c(seqnames,CI.start,CI.end,width,strand,genoT,deltaW)))
    tmp$confint@seqinfo <- seqinfo
    levels(tmp$confint@seqnames) <- seqnameLevels
    
    save(tmp, file=paste0(cleanDatapath,"/",basename(file)))
  }
  
}
