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
library(rtracklayer )

#######################################################
####                  Data                     ###
#######################################################

#write.table(hotspots,"hotspots.txt",quote = F,row.names = F,col.names = T,sep = "\t")
hotspots=read.table("hotspots_zh_data.txt",header = T)
hotspots=read.table("hotspots.txt",header = T)


#######################################################
####                  Genotype                    ###
#######################################################

gene = data.frame(library=list.files("SCEs/BLM-RECQL5/Data-BreakpointR/data/"))
gene$gene  <- "gene"
gene$ID2 = gene$library
suppressWarnings(gene <- gene %>% separate(ID2, c("a","b","c","d","e","f"), "[_-]+"))
gene$gene <- "gene"
for (row in 1:nrow(gene)){
  
  for (letter in c("a","b","c","d","e","f")){
    if (is.na(gene[row,letter])!=T){
      if (gene[row,letter]=="WT" | gene[row,letter]=="wt"){
        gene[row,"gene"]="WT"
      }
      else if (gene[row,letter]=="blm" | gene[row,letter]=="BLM" ) {
        gene[row,"gene"]="BLM"
      }
      
      else if (gene[row,letter]=="RECQL5" | gene[row,letter]=="recql5" | gene[row,letter]=="RECQ5" | gene[row,letter]=="recq5" ) {
        if (gene[row,"gene"]=="BLM"){
          gene[row,"gene"]="BLM/RECQL5"
        }
        else{
          gene[row,"gene"]="RECQL5"
        }
      }
    }
  }
  
}
gene= select(gene,-c(a,b,c,d,e,f))
gene = filter(gene,gene=="BLM/RECQL5")

breaks.all.files = collectBreaksAllFiles(paste0("../StraVa/DATA/rdata/",gene$library))
hotspots = breakpointHotspotter(breaks.all.files)


#######################################################
####     Matrix of count signatures (n x m)       ###
#######################################################

hotspots$count<-as.factor(hotspots$count)
hotspots$ID<-as.factor(hotspots$ID)
df = data.frame(matrix(0,ncol=length(levels(hotspots$count)),nrow=length(levels(hotspots$ID))))
hotspots$count = as.factor(paste0(hotspots$count,"-",hotspots$seqnames))
colnames(df) = levels(hotspots$count)
rownames(df) = levels(hotspots$ID)

level=levels(hotspots$count)[1]

for (level in levels(hotspots$count)){
  tmp = filter(hotspots, count==level)
  #print(level)
  for (row in 1:length(df[,level])){
    #print(row)
    if (row.names(df[row,]) %in% tmp$ID  ){ df[row,level] = 1}
  }
}


#######################################################
####  Matrix comparing hotspot signatures (n x n)    ###
#######################################################


#levels(hotspots$seqnames) <- factor(c("chr1" ,"chr2","chr3",  "chr4",  "chr5",  "chr7" , "chr8" , "chr9" , "chr10", "chr11", "chr12" ,"chr13", "chr14" ,"chr15" ,"chr16", "chr17" ,"chr19", "chr20" ,"chr21", "chr22" , "chrX" ))



matrix=data.frame(matrix(0,ncol=length(levels(hotspots$count)),nrow=length(levels(hotspots$count))))
colnames(matrix) = levels(hotspots$count)
rownames(matrix) = levels(hotspots$count)


#for (col in 1:ncol(df)){
#  tmp = df[,col]
#  for (i in 1:ncol(df)){
#    matrix[col,i]=1- ( sum(stringdist(df[,i],tmp,method="lv"))/length(tmp))
#  }}
col=colnames(df)[1]
tmp=df[,"38"]
i=colnames(df)[2]
for (col in colnames(df)){
  tmp = df[,col]
  
  
  for (i in colnames(df)){
    vec= df[,i] + tmp
    matrix[col,i]=length(vec[vec=="2"])/(sum(vec)/1)
    #matrix[col,i] = 1- ( sum(stringdist(df[,i],tmp,method="lv"))/ (sum(df[,i]) + sum(tmp))        )
  }
}




dt2 <- matrix %>%
  tibble::rownames_to_column() %>%
  gather(colname, value, -rowname)
dt2$rowname <- as.factor(dt2$rowname)
dt2$colname <- as.factor(dt2$colname)


levels(dt2$rowname) = factor(c( "chr1-1"  ,"chr2-3"  , "chr2-5" ,  "chr4-11" , "chr4-7"  , "chr4-9"  , "chr7-13",  "chr8-15" , "chr9-17" ,
                                "chr10-19" ,"chr10-21", "chr12-23", "chr13-25" ,"chr13-26" ,"chr14-28", "chr15-30" ,"chr16-32" ,"chr17-34","chr17-36" ,"chr19-38" ,"chr19-39" ,"chr20-41", "chr20-42" ,"chr20-43" ,"chr21-45" ,"chr22-47" ,"chrX-49" , "chrX-50" ))
levels(dt2$colname) = factor(c("chr1-1"  ,"chr2-3"  , "chr2-5" ,  "chr4-11" , "chr4-7"  , "chr4-9"  , "chr7-13",  "chr8-15" , "chr9-17" ,
                               "chr10-19" ,"chr10-21", "chr12-23", "chr13-25" ,"chr13-26" ,"chr14-28", "chr15-30" ,"chr16-32" ,"chr17-34","chr17-36" ,"chr19-38" ,"chr19-39" ,"chr20-41", "chr20-42" ,"chr20-43" ,"chr21-45" ,"chr22-47" ,"chrX-49" , "chrX-50" ))

  
  
  
levels(dt2$rowname) = factor(c(  "chr1-1" ,"chr1-2"  ,"chr2-4"  , "chr2-5"  , "chr3-7" ,  "chr4-9" ,  "chr5-11" , "chr5-13" , "chr7-15" , "chr7-16" ,
                                 "chr7-17" , "chr8-19" , "chr8-20",  "chr8-21"  ,"chr9-23" , "chr9-24" , "chr9-25", "chr10-27",  "chr11-29" ,"chr11-30" , "chr12-32" ,"chr13-34", "chr14-36", "chr15-38", "chr15-39", "chr15-40", "chr15-41", "chr16-43", "chr16-44", "chr16-45", "chr17-47" ,"chr17-49", "chr19-51" ,"chr19-52" ,"chr20-54", "chr21-56", "chr22-58","chr22-59",  "chrX-61" , "chrX-62" ))
levels(dt2$colname) = factor(c(  "chr1-1" ,"chr1-2"  ,"chr2-4"  , "chr2-5"  , "chr3-7" ,  "chr4-9" ,  "chr5-11" , "chr5-13" , "chr7-15" , "chr7-16" ,
                                 "chr7-17" , "chr8-19" , "chr8-20",  "chr8-21"  ,"chr9-23" , "chr9-24" , "chr9-25","chr10-27",  "chr11-29" ,"chr11-30" , "chr12-32" ,"chr13-34", "chr14-36", "chr15-38", "chr15-39", "chr15-40",  "chr15-41", "chr16-43", "chr16-44", "chr16-45", "chr17-47" ,"chr17-49", "chr19-51" ,"chr19-52" ,"chr20-54", "chr21-56", "chr22-58", "chr22-59",  "chrX-61" , "chrX-62" ))


levels(dt2$rowname) = factor(c(  "1-chr1" ,  "2-chr1"  ,"4-chr2"  ,"5-chr2" ,"7-chr3"  , "8-chr3","10-chr4",  "12-chr5" , "13-chr5",  "15-chr6" , "17-chr7" , "19-chr8" ,  "20-chr8" , "21-chr8" , "23-chr9" , "24-chr9" , "26-chr9" 
                                ,"27-chr9" , "29-chr9" , "30-chr9" , "32-chr12", "34-chr13" ,"36-chr13" ,"38-chr15", "40-chr16", "41-chr16", "42-chr16" ,"44-chr17"
                                ,"45-chr17" ,"47-chr17" ,"48-chr17",  "50-chr19" ,"52-chr20", "54-chr21" ,"55-chr21", "57-chr22", "59-chrX"  ,  "60-chrX" ))
levels(dt2$colname) =factor(c(  "1-chr1" ,  "2-chr1"  ,"4-chr2"  ,"5-chr2" ,"7-chr3"  , "8-chr3","10-chr4",  "12-chr5" , "13-chr5",  "15-chr6" , "17-chr7" , "19-chr8" ,  "20-chr8" , "21-chr8" , "23-chr9" , "24-chr9" , "26-chr9" 
                                ,"27-chr9" , "29-chr9" , "30-chr9" , "32-chr12", "34-chr13" ,"36-chr13" ,"38-chr15", "40-chr16", "41-chr16", "42-chr16" ,"44-chr17"
                                ,"45-chr17" ,"47-chr17" ,"48-chr17",  "50-chr19" ,"52-chr20", "54-chr21" ,"55-chr21", "57-chr22", "59-chrX"  ,  "60-chrX" ))



ggplot(dt2, aes(x = rowname, y = colname, fill = value)) +
  geom_tile() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradient2(
    low = "#213859",
    mid= "#213859",
    high = "#eb4034",
    midpoint = 0.35
  ) #+ ggsave("TRANSLOCATIONS/zh_317_good.png")
ggplot(dt4, aes(x = as.numeric(chrA), y = as.numeric(chrB), fill = value)) +
  geom_tile() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradient2(
    low = "#213859",
    mid= "#213859",
    high = "#eb4034",
    midpoint = 0.35
  )

dt3 = dt2 %>% separate(rowname,c("chrA","hotspot1"),"-")  
dt4 = dt3 %>% separate(colname,c("chrB","hotspot2"),"-")

#dt2 = dplyr::filter(dt2, dt2$value > 0.75)
translocations = data.frame()
for (row in 1:nrow(dt4)){
  if (dt4$value[row] > 0.4){
      if (dt4[row,2]!= dt4[row,4]){
      print(dt4[row,])
      row = dt4[row,]
      translocations=rbind(translocations,row)
      }
  }
}

translocations=translocations[c(),]
A="51"
B="40"
row=7
for (row in 1:nrow(translocations)){
  
  A=translocations[row,1]
  B=translocations[row,3]
  tmp <- filter(hotspots,seqnames %in% c(A,B))
  chrA = as.character(tmp[1,1])
  chrB = as.character(tail(tmp$seqnames,n=1))
  transpath = paste0("TRANSLOCATIONS/",chrA,"-",chrB,"/")
  ctrlpath=paste0(transpath,"Control/")
  if (!file.exists(transpath) ) { dir.create(transpath)}
  if (!file.exists(ctrlpath) ) { dir.create(ctrlpath)}
  chrAdf = filter(tmp, seqnames==A)
  chrBdf = filter(tmp, seqnames==B)
  
  #plotBreakpointsPerChr(files2plot = paste0("../NA878/data_filtered/",chrAdf$ID,"filtered.RData"),plotspath = transpath,chromosomes = chrA)
  #plotBreakpointsPerChr(files2plot = paste0("../NA878/data_filtered/",chrBdf$ID,"filtered.RData"),plotspath = transpath,chromosomes = chrB)
  
  if (length(intersect(chrAdf$ID,chrBdf$ID))<15){
    plotBreakpointsPerChr(files2plot = paste0("../StraVa/DATA/rdata/",intersect(chrAdf$ID,chrBdf$ID),".RData"),plotspath = transpath,chromosomes = chrA)
    plotBreakpointsPerChr(files2plot = paste0("../StraVa/DATA/rdata/",intersect(chrAdf$ID,chrBdf$ID),".RData"),plotspath = transpath,chromosomes = chrB)
  } else {
    plotBreakpointsPerChr(files2plot = paste0("../StraVa/DATA/rdata/",intersect(chrAdf$ID,chrBdf$ID)[1:10],".RData"),plotspath = transpath,chromosomes = chrA)
    plotBreakpointsPerChr(files2plot = paste0("../StraVa/DATA/rdata/",intersect(chrAdf$ID,chrBdf$ID)[1:10],".RData"),plotspath = transpath,chromosomes = chrB)
  }
  plotBreakpointsPerChr(files2plot = paste0("../StraVa/DATA/rdata/",setdiff( hotspots$ID,intersect(chrAdf$ID,chrBdf$ID))[1:10],".RData"),plotspath = ctrlpath,chromosomes = chrA)
  plotBreakpointsPerChr(files2plot = paste0("../StraVa/DATA/rdata/",setdiff( hotspots$ID,intersect(chrAdf$ID,chrBdf$ID))[1:10],".RData"),plotspath = ctrlpath,chromosomes = chrB)
  

}


x=filter(chrA_inter,midpoints > 130777373,midpoints < 130851519)$midpoints
x=filter(chrB_inter,midpoints >23270202,midpoints < 23310774)$midpoints
#x=chrB_inter$midpoints

# Create a histogram
hist(x, freq = FALSE, main = "Histogram and density")

# Calculate density
dx <- density(x)

# Add density
lines(dx, lwd = 2, col = "red")

# Plot the density without histogram
plot(dx, lwd = 2, col = "blue",
     main = "Density",xlab="CHROMOSOME 9",cex.axis=1.2,cex.lab=1.5,yaxt="n",font.lab=2)

#axis(1)

# Add the data-poins with noise in the X-axis
rug(jitter(x))
abline(v=(mean(x)-(1*(sd(x)))),col="red",lwd=3,lty=2)
abline(v=(mean(x)+(1*(sd(x)))),col="red",lwd=3,lty=2)
mean(x)-(1*(sd(x)))
mean(x)+(1*(sd(x)))
