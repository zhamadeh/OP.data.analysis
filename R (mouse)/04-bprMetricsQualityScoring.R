#######################################################
####     Workflow Part 4: Quality Scoring           ###
#######################################################

library(stringr)
library(dplyr)
library("ggplot2")
# Load ggplot2 package

qualityMetrics <- function(){
  #load data
  df = read.table("INPUT/SUMMARY/MOUSE/metrics_details.txt",header=T)
  merge=df
  
  merge$quality=NA
  
  merge$quality[merge$coverage<25]="poor"
  merge$quality[merge$coverage>=25]="okay"
  merge$quality[merge$coverage>=50]="acceptable"
  merge$quality[merge$coverage>=100]="good"
  merge$quality[merge$coverage>=175]="very_good"

  merge$quality[merge$background>=0.079]="okay"
  merge$quality[merge$evenness.mean>=0.21]="okay"
  merge$quality[merge$spikiness>=0.172]="okay"
  
  write.table(merge,"INPUT/SUMMARY/MOUSE//05.library.quality.txt")
}


