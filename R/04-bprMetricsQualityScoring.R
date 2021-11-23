#######################################################
####     Workflow Part 4: Quality Scoring           ###
#######################################################

library(stringr)
library(dplyr)
library("ggplot2")
# Load ggplot2 package

qualityMetrics <- function(){
  #load data
  df = read.table("INPUT/2021/02.metrics_bam_summ.txt",header=T)
  merge=df
  merge$background = NA
  merge$coverage = NA
  
  merge <- filter(merge,merge$Quality!="NA")
  
  for (i in 1:nrow(merge)){
  	filename = paste0("INPUT/DATA/BREAKPOINTR/2021-BLM-RECQL5/",merge[i,1],".processed.bam.RData")
  	tmp = get(load(filename))$lib.metrics[1][[1]]
  	tmp2 = get(load(filename))$lib.metrics[2][[1]]
  	merge[i,]$background = tmp
  	merge[i,]$coverage = tmp2
  }
  
  merge$quality=NA
  
  merge$quality[merge$coverage<25]="poor"
  merge$quality[merge$coverage>=25]="okay"
  merge$quality[merge$coverage>=50]="acceptable"
  merge$quality[merge$coverage>=100]="good"
  merge$quality[merge$coverage>=175]="very_good"

  merge$quality[merge$background>=0.079]="okay"
  merge$quality[merge$evenness.mean>=0.21]="okay"
  merge$quality[merge$spikiness>=0.172]="okay"
  
  write.table(merge,"INPUT/2021/04.library.quality_2021.txt")
}


