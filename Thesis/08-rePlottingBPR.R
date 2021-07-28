#######################################################
####         Workflow Part 8: Re-plotting BPR       ###
#######################################################

breakpointR::plotBreakpoints(files2plot = list.files("Output/bpr/data.blacklisted/",full.names = T),file = "Output/breaksPlot_bl.pdf")


ggplot(breaks)+geom_bar()