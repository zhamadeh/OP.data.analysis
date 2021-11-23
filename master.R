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
breakpoints <- collectBreaksAllFiles()
breakpoints = breakpoints %>% separate(Library,c("Library","b"),"[.]") %>% select(-c(b,library))
write.table(breakpoints,"INPUT/2021/03.library.breakpoints.txt",col.names = T,row.names = F,quote = F,sep = "\t")

s# Part 4
bammetrics <- read.table("INPUT/2021/01.library.metrics_2021.txt",header=T)
metrics <- read.table("INPUT/2021/00.metrics.summary_2021.txt",header=T)
bammetrics = bammetrics %>% separate(file,c("Library","b"),"[.]") %>% select(-c(b))
metrics_bam_summ = merge(metrics,bammetrics,by="Library")
write.table(metrics_bam_summ,"INPUT/2021/02.metrics_bam_summ.txt",row.names = F,col.names = T,sep = "\t",quote = F)
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

