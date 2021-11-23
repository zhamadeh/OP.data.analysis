#######################################################
####     Workflow Part 5: Plot Pairs      ###
#######################################################

library(stringr)
library(dplyr)
library("ggplot2")
# Load ggplot2 package

plottingPairs <- function(){
  
  merge=read.table("INPUT/2021/04.library.quality_2021.txt",header=T)
  merge$quality<-as.factor(merge$quality)
  merge$quality
  my_cols=c("#c98d26","#32a852","red")
  
  panel.cor <- function(x, y){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y), digits=2)
    txt <- paste0("R = ", r)
    cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
  }
  # Customize upper panel
  upper.panel<-function(x, y){
    points(x,y, pch = 19, col = my_cols[merge$quality])
  }
  
  panel.hist <- function(x, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
  }
  
  # Create the plots
  pairs(merge[, c("coverage","background","spikiness","Quality","evenness.mean", "Coverage_after_1Gb_sequencing")],
        lower.panel = panel.cor,
        upper.panel = upper.panel,
        diag.panel = panel.hist) +
    ggsave("2021.pairs.libraryQuality.png")
  
  
}
