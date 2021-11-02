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
write.table(breakpoints,"Input/03.library.breakpoints.txt",col.names = T,row.names = F,quote = F,sep = "\t")

# Part 4
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

