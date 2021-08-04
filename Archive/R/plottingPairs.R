library(stringr)
library(dplyr)
library("ggplot2")
# Load ggplot2 package

qualityMetrics <- function(){
  #load data
  df = read.table("Input/all.bam.new.metrics.txt",header=T)
  df=df[- grep("na", df$file),]
  df=df[- grep("Undetermined", df$file),]
  merge=df
  merge$background = NA
  
  for (i in 1:nrow(merge)){
  	filename = paste0("Output/bpr/data/",merge[i,1],".RData")
  	tmp = get(load(filename))$lib.metrics[1][[1]]
  	tmp2 = get(load(filename))$lib.metrics[2][[1]]
  	merge[i,]$background = tmp
  	merge[i,]$coverage = tmp2
  }
  
  merge$quality=NA
  
  merge$quality[merge$coverage<25]="poor"
  merge$quality[merge$coverage>=25]="okay"
  merge$quality[merge$coverage>=50]="acceptable"
  merge$quality[merge$coverage>=100]="good"
  merge$quality[merge$coverage>=175]="very_good"

  merge$quality[merge$background>=0.079]="okay"
  merge$quality[merge$evenness.mean>=0.21]="okay"
  merge$quality[merge$spikiness>=0.172]="okay"
  
  write.table(merge,"merge.metrics.background.quality.txt")
}

plottingPairs <- function(){
  
  merge=read.table("merge.metrics.background.csv",header=T)
  merge=merge[- grep("na", merge$file),]
  
  my_cols=c("#b1c926", "#32a852", "#c98d26","red","#0acca5")

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
  pairs(merge[, c("spikiness","evenness.med","evenness.mean","coverage","background")],
  	  lower.panel = panel.cor,
  	  upper.panel = upper.panel,
  	  diag.panel = panel.hist) +
    ggsave("Output/pairs.all.png")
  
  
  
  
  
  brdu = dplyr::filter(merge,brdu=="yes")
  good = dplyr::filter(merge,quality %in% c("good","very_good"))
  merge= rbind(brdu,good)
  
  
  my_cols=c("green","red")
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
    points(x,y, pch = 19, col = my_cols[merge$brdu])
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
  pairs(merge[, c("spikiness","evenness.med","evenness.mean","coverage","background")],
        lower.panel = panel.cor,
        upper.panel = upper.panel,
        diag.panel = panel.hist)+
    ggsave("Output/pairs.brdu.png")
}

