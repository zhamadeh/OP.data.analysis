#######################################################
####         Workflow Part 8: Re-plotting BPR       ###
#######################################################

rePlottingBPR <- function(){
  breakpointR::plotBreakpoints(files2plot = list.files("Output/bpr/data.blacklisted/",full.names = T),file = "Output/breaksPlot_bl.pdf")
}