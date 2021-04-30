#library(M3C)
#install.packages("factoextra")
#library(factoextra)

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
  theme(text = element_text(size=18))
