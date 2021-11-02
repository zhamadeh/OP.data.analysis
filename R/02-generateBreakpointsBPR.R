#######################################################
####     Workflow Part 2: Generate breakpoints      ###
#######################################################

### Run breakpointR and generate rdata files ###

runBreakpointR <- function(dir){
  breakpointr(inputfolder=dir, outputfolder="Output/bpr/", pairedEndReads=TRUE, numCPU=2,windowsize=175,binMethod="reads",peakTh=0.3875,min.mapq=7.75,trim=6.5,background=0.15)
}
