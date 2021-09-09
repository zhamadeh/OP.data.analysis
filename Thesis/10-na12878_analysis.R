##########################################################
######         Workflow Part 10: NA878 cells       #####
##########################################################



######################## Collect Data ######################

# Datasets
na878_breakpoints = collectBreaksAllFiles(datapath = "../NA878/data/") # 59K
quality = read.table("../NA878/quality.lansdorp_na12878.txt",header=F)

# List of good libraries
quality_good = filter(quality, V2 > 0.7)$V1 

# Merge list of good libraries
na878_breakpoints_quality = merge(na878_breakpoints, quality, by.x="library",by.y="V1") # 28k 
na878_breakpoints_quality_good = filter(na878_breakpoints_quality, V2 > 0.7) # 28k 




######################## Filtering data ######################

# 1) Remove homozygous events
na878_breakpoints_quality_good_SCE = filter(na878_breakpoints_quality_good,genoT!= "cc-ww" | genoT!= "cc-ww") # 26k 

# 2) Filter out events that are too close to each other (2Mb)
na878_breakpoints_quality_good_SCE$library <- as.factor(na878_breakpoints_quality_good_SCE$library)
filteredSCEs = data.frame()
n=1
for (level in levels(na878_breakpoints_quality_good_SCE$library)){
  message(  round(     (      (n/length(levels(na878_breakpoints_quality_good_SCE$library)))   *100        )     ,digits = 2)   ,"% ... complete"   )
  tmp = filter(na878_breakpoints_quality_good_SCE, library==level)
  tmp$seqnames <- droplevels(tmp$seqnames)
  
  level2=levels(tmp$seqnames)[2]
  for (level2 in levels(tmp$seqnames)){
    
    tmp2 = filter(tmp, seqnames==level2)
    tmp3 <- GRanges(tmp2)
    
    overlaps = countOverlaps(tmp3, tmp3, type="any",maxgap = 2065000)
    overlap_rows = which(overlaps %in% c(1))
    
    tmp4 = tmp2[overlap_rows,]

    filteredSCEs <- rbind(filteredSCEs,tmp4)
    
  }
  n=n+1
}

# 3) Filter out events that are too close to each other (2Mb)
na878_breakpoints_quality_good_SCE_filtered = filteredSCEs # 18 k 
na878_breakpoints_quality_good_SCE_filtered_10kb <- filter(na878_breakpoints_quality_good_SCE_filtered, width < 10000) # 3k 

na878_breakpoints_quality_good_SCE_filtered_10kb$library<-droplevels(na878_breakpoints_quality_good_SCE_filtered_10kb$library)
length(levels(na878_breakpoints_quality_good_SCE_filtered_10kb$library))
length(levels(na878_breakpoints_quality_good_SCE_filtered$library))


######################## Exporting Data ######################

# Write to enrichment folder for enrichment analysis
write.table(select(na878_breakpoints_quality_good_SCE_filtered_10kb,c(seqnames, start,end)), "Thesis/Enrichment/GenomePermute/na878_sces.bed", col.names = F, row.names = F, quote = F, sep="\t")




######################## Plotting Data ######################
# Take putative SCEs and reprint breaksPlot with just them
n=1
for (file in quality_good){
  
  message(  round(     (      (n/length(quality_good))   *100        )     ,digits = 2)   ,"% ... complete"   )
  
  filename = paste0("../NA878/data/",file,".RData")
  tmp = get(load(filename))
  
  
  seqinfo <- tmp$breaks@seqinfo
  seqnameLevels <- levels(tmp$breaks@seqnames)
  
  tmp$breaks <-GRanges(filter(na878_breakpoints_quality_good_SCE_filtered,library==tmp$ID) %>% select(c(seqnames,start,end,width,strand,genoT,deltaW)))
  tmp$breaks@seqinfo <- seqinfo
  levels(tmp$breaks@seqnames) <- seqnameLevels
  
  
  tmp$confint = tmp$confint[queryHits(findOverlaps(tmp$confint, tmp$breaks, type="any")),]
  #$confint<-GRanges(filter(breaks,library==tmp$ID) %>% select(c(seqnames,CI.start,CI.end,width,strand,genoT,deltaW)))
  tmp$confint@seqinfo <- seqinfo
  levels(tmp$confint@seqnames) <- seqnameLevels
  
  save(tmp, file=paste0("../NA878/data_filtered/",file,"_filtered.RData"))
  
  n=n+1
}

breakpointR::plotBreakpoints(files2plot = paste0("../NA878/data_filtered/",quality_good,"_filtered.RData"),file = "breaksPlot_na12878_good70_filtered.pdf")
