library(stringr)
library(dplyr)
library("ggplot2")
# Load ggplot2 package


#load data
df = read.table("Input/all.bam.new.metrics2.txt",header=T)
quality <- read.csv("../library.quality.annotated.2021.csv",sep="\t")
quality$quality<-as.factor(quality$quality)
levels(quality$quality) = c( "acceptable" ,"good", "okay", "poor", "really_good")

interesting.libraries = quality[quality$cool=="***",]$filename


qt = select(quality,c(filename,quality))
merge = merge(qt,df,by.x="filename",by.y="file")
merge$background = NA

for (i in 1:nrow(merge)){
	filename = paste0("../StraVa/DATA/rdata/",merge[i,1],".RData")
	tmp = get(load(filename))$lib.metrics[1][[1]]
	tmp2 = get(load(filename))$lib.metrics[2][[1]]
	merge[i,]$background = tmp
	merge[i,]$coverage = tmp2
}





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
	  diag.panel = panel.hist)


