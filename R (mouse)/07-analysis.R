##########################################################
######         Workflow Part 10: NA878 cells       #####
##########################################################
library(dplyr)
library(GenomicRanges)
library(tidyr)
library(ggplot2)
library(breakpointR)


######################## Collect Data ######################

# Datasets

mouse = read.table("INPUT/SUMMARY/MOUSE/06.gene.full.txt",header = T)
# Merge list of good libraries
mouse$Library <- as.factor(mouse$Library)
length(levels(mouse$Library))

######################## Filtering data ######################

# 1) Remove homozygous events
mouse = filter(mouse,genoT!= "cc-ww" | genoT!= "cc-ww") # 26k 

mouse.df = data.frame()
n=1
for (level in levels(mouse$Library)){
  message(  round(     (      (n/length(levels(mouse$Library)))   *100        )     ,digits = 2)   ,"% ... complete"   )
  tmp = filter(mouse, Library==level)
  tmp$seqnames <- droplevels(tmp$seqnames)
  
  level2=levels(tmp$seqnames)[1]
  for (level2 in levels(tmp$seqnames)){
    
    tmp2 = filter(tmp, seqnames==level2)
    tmp3 <- GRanges(tmp2)
    
    overlaps = countOverlaps(tmp3, tmp3, type="any",maxgap = 2500000)
    overlap_rows = which(overlaps %in% c(1))
    
    tmp4 = tmp2[overlap_rows,]

    mouse.df <- rbind(mouse.df,tmp4)
    
  }
  n=n+1
}
mouse.df$Library<-droplevels(mouse.df$Library)
length(levels(mouse.df$Library))

# 4) Filter out events that are too close to each other (2Mb)
mouse.df_10kb <- filter(mouse.df, width < 150000) # 3k 
mouse.df_10kb$Library<-droplevels(mouse.df_10kb$Library)
                                      

######################## Exporting Data ######################

# Write to enrichment folder for enrichment analysis
write.table(select(mouse.df_10kb,c(seqnames, start,end)), "R (mouse)/Enrichment/sces.bed", col.names = F, row.names = F, quote = F, sep="\t")


######################## Plotting Data ######################
# Take putative SCEs and reprint breaksPlot with just them

cleanRData <- function(data=na878_breakpoints_good_hetero_SCE_CM_bl){
  n=1
  for (file in levels(data$Library)){
    
    message(  round(     (      (n/length(levels(data$Library)))   *100        )     ,digits = 2)   ,"% ... complete"   )
    filename = paste0("INPUT/DATA/BREAKPOINTR/MOUSE/",file,".RData")
    tmp = get(load(filename))
    
    seqinfo <- tmp$breaks@seqinfo
    seqnameLevels <- levels(tmp$breaks@seqnames)
    
    tmp$breaks <-GRanges(filter(data,Library==tmp$ID) %>% select(c(seqnames,start,end,width,strand,genoT,deltaW)))
    tmp$breaks@seqinfo <- seqinfo
    levels(tmp$breaks@seqnames) <- seqnameLevels
    
    tmp$confint = tmp$confint[queryHits(findOverlaps(tmp$confint, tmp$breaks, type="any")),]
    #$confint<-GRanges(filter(breaks,Library==tmp$ID) %>% select(c(seqnames,CI.start,CI.end,width,strand,genoT,deltaW)))
    tmp$confint@seqinfo <- seqinfo
    levels(tmp$confint@seqnames) <- seqnameLevels
    
    save(tmp, file=paste0("INPUT/DATA/BREAKPOINTR/MOUSE_BL/",file,"filtered.RData"))
    
    n=n+1
  }
  
}

cleanRData(mouse.df_10kb)
breakpointR::plotBreakpoints(files2plot = list.files("INPUT/DATA/BREAKPOINTR/MOUSE_BL/",full.names = T),file = "mouse_breaksPlot_bl.pdf")


