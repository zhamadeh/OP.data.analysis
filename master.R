library(breakpointR)
source("R/readCountPlotting.R")
source("R/collectLibraryMetrics.R")
source("R/plottingPairs.R")

args = commandArgs(trailingOnly=TRUE)

### Collect library metrics from BAM files ### 

collectLibraryStats(folder = args[1])

### Run breakpointR and generate rdata files ###

breakpointr(inputfolder=args[1], outputfolder="Output/bpr/", pairedEndReads=TRUE, numCPU=2,windowsize=175,binMethod="reads",peakTh=0.3875,min.mapq=7.75,trim=6.5,background=0.15)

### Quality metrics assembly

qualityMetrics()

### Plot pairs

plottingPairs()

### Run breakpointR and generate rdata files ###

readPlotting(rdata="Output/bpr/data/",plot.dir = "Output/bpr/plots/",cluster.metrics="merge.metrics.background.quality.txt",features=c("coverage","background","spikiness","evenness.mean","good"),numOfLibs=10)




### PCA analysis ### 

merge <- read.table("merge.quality.metrics.complete.txt",header=T)
merge.pca =  prcomp(merge[,c(3:7)],center = T,scale. = T)

fviz_pca_ind(merge.pca, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = merge$quality, 
             col.ind = "black", 
             palette = c("#b1c926", "#32a852", "#c98d26","red","#0acca5"), 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Library Quality")+
  theme(text = element_text(size=18)) +
  save("Output/pca.png")
