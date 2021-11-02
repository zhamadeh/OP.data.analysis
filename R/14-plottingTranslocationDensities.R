#######################################################
####                  Packages                     ###
#######################################################

library(breakpointR)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringdist)
library(tibble)
library(rtrack )

#######################################################
####                  Data                     ###
#######################################################

hotspots=read.table("hotspots_zh_data.txt",header = T)

chrA = filter(hotspots,count=="25")
chrB = filter(hotspots,count=="59")

shared = intersect(chrA$ID,chrB$ID)
chrB_inter = filter(chrA,ID %in% shared)
chrA_inter = filter(chrA,ID %in% shared)
chrB_inter = filter(chrB,ID %in% shared)
chrA_inter$midpoints = chrA_inter$start + (chrA_inter$width/2)

x= chrA_inter$midpoints
x=filter(chrA_inter,midpoints > 130777373,midpoints < 130851519)$midpoints
x=filter(chrB_inter,midpoints >23270202,midpoints < 23310774)$midpoints

chrA_inter
min(chrA_inter$start)
max(chrA_inter$end)

chrA_inter=filter(chrA_inter,start > 130780000,end < 130820000)


length = max(chrA_inter$end)-min(chrA_inter$start)
chrA_inter$new_start=chrA_inter$start-min(chrA_inter$start)
chrA_inter$new_end=chrA_inter$end-min(chrA_inter$start)

vec = c(rep(0,1,length ))
row=1
for (row in 1:nrow(chrA_inter)){
  tmp = (chrA_inter[row,])
  #tmp
  weight = 1/tmp$width
  #weight
  #length(vec)
  #tmp$width
  #length(vec[(tmp$new_start):(tmp$new_end)])
  #sum(vec[(tmp$new_start):(tmp$new_end)])
  vec[(tmp$new_start):(tmp$new_end)]=vec[tmp$new_start:tmp$new_end]+weight
}

t = as.data.frame(vec)
x=vec
t$position = 1:length(vec)
t$coordinate=t$position+min(chrA_inter$start)
ggplot(t)+geom_density(aes(coordinate))
# Create a histogram
hist(x, freq = FALSE, main = "Histogram and density")
# Calculate density
dx <- density(x)
# Add density
lines(dx, lwd = 2, col = "red")
# Plot the density without histogram
plot(dx, lwd = 2, col = "blue",
     main = "Density",xlab="CHROMOSOME 9",cex.axis=1.2,cex.lab=1.5,yaxt="n",font.lab=2)


# Add the data-poins with noise in the X-axis
rug(jitter(x))
abline(v=(mean(x)-(1*(sd(x)))),col="red",lwd=3,lty=2)
abline(v=(mean(x)+(1*(sd(x)))),col="red",lwd=3,lty=2)

mean(x)-(1*(sd(x)))
mean(x)+(1*(sd(x)))