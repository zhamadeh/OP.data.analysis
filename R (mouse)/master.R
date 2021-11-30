#######################################################
####                  Master Script                 ###
#######################################################


#Packages
install.packages("BiocManager")
BiocManager::install("breakpointR")
library(breakpointR)
for (file in list.files("Thesis/",full.names = T)){
  source(file)
}

#### Running individual functions from each script ####

# Part 1
collectLibraryStats("../../Data/BAM/All_good_libraries/ALL_GOOD_BAM_FILES/")

# Part 2
runBreakpointR("../../Data/BAM/All_good_libraries/ALL_GOOD_BAM_FILES/")

# Part 3
breakpoints <- collectBreaksAllFiles(datapath = "INPUT/DATA/BREAKPOINTR/MOUSE/")
breakpoints = breakpoints %>% separate(library,c("library","b"),"[.]") %>% select(-c(b))
write.table(breakpoints,"INPUT/SUMMARY/MOUSE/03.library.breakpoints.txt",col.names = T,row.names = F,quote = F,sep = "\t")

s# Part 4
bammetrics <- read.table("INPUT/SUMMARY/MOUSE/01.library.metrics_mouse.txt",header=T)
metrics <- read.table("INPUT/SUMMARY/MOUSE/02.metrics.bam.summary.txt",header=T)

#do one of these to make filenames the same
bammetrics = bammetrics %>% separate(file,c("Library","b"),"[.]") %>% select(-c(b))
metrics$Library <- paste0(metrics$Library,".processed.bam")

metrics_bam_summ = merge(metrics,bammetrics,by.x="Library",by.y="file")
write.table(metrics_bam_summ,"INPUT/SUMMARY/MOUSE/04.metrics_full.txt",row.names = F,col.names = T,sep = "\t",quote = F)
qualityMetrics()

# Part 5
plottingPairs()

# Part 6
qualityFilterLibraries()

# Part 7
frequencyFilterBreakpoints(blacklist="Input/00.centromeres2.txt")

# Part 8
rePlottingBPR()

# Part 9 will change everytime depending on data
sces <- read.table("Input/06.sisterChromatidExchanges.txt",header=T)

