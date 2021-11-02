
######################################################################
####   Packages   ####
######################################################################

library(rtracklayer)
library(tidyverse)
library(dplyr)
library(plyr)
library(tidyr)
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(GenomicRanges))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(GenomicAlignments))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(DNAcopy))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(doParallel))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(tiydyverse))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(ggplot2))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(cowplot))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(caTools))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(gplots))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(ReorderCluster))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(BiocManager))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(BiocManager::install("AneuFinder"))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(AneuFinder))))
######################################################################
####   Assemble datasets   ####
######################################################################

cutoff=3000
file = "binsize_5e+05_stepsize_5e+05_StrandSeq_CNV.bed.gz"

copyNumberFinder <- function(cutoff, file){
  
  CNV<- import(file)
  cnvPerCell<-data.frame()
  cnvPerCellSummary=data.frame()
  
  for (i in 1:length(CNV)){
    
    ## ASSIGN GENOTYPE BASED OFF FILE NAME
    file=strsplit(CNV[i]@listData[[1]]@trackLine@description,split = " ")[[1]][4]
    
    message("Reading file: ",file," ... ",round((i/length(CNV))*100,2),"%")
    
    ## RETRIEVE FIRST LIBRARY AND REMOVE SMALL SEGMENTS < 500Kb and HIGH PLOIDY STATES TO ASSIGN  PLOIDY
    tmp <- as.data.frame(CNV[i])
    tmp$name <- as.factor(tmp$name)
    tmp=tmp[tmp$width>100000,]
    tmp=tmp[tmp$seqnames!="chrY",]
    #HIGH PLOIDY STATES
    tmp <- tmp[tmp$name!="20-somy" & tmp$name!="19-somy" &tmp$name!="18-somy" &tmp$name!="17-somy" & tmp$name!="16-somy" &tmp$name!="15-somy" &tmp$name!="14-somy"&tmp$name!="13-somy" &tmp$name!="12-somy"&tmp$name!="11-somy"&tmp$name!="10-somy"&tmp$name!="zero-inflation",]
    tmp$name <-droplevels(tmp$name)
    
    
    ## BUILD PLOIDY TABLE SUMMARY TO CLASSIFY NATIVE STATE AND GAINS & LOSSES
    ploidyTable <- tmp %>% group_by(name) %>% dplyr::summarize(n(),sum=sum(width))
    ploidy = as.numeric(strsplit(as.character(ploidyTable[which.max(ploidyTable$sum),]$name),split = "[-]")[[1]][1])
    if (ploidy==0){
      ploidy=1
    } else if (ploidy >2){
      next
    }
    
    
    ## REMOVE SMALL SEGMENTS  AFFTER HAVING  ASSIGNED PLOIDY  STATE
    tmp=tmp[tmp$width>cutoff,]
    ploidyTable <- tmp %>% group_by(name) %>% dplyr::summarize(n(),sum=sum(width))
    
    
    ## TURN PLOIDY INTO NUMERIC AND CLASSIFY  GAINS/LOSSES
    ploidyTable=  separate(ploidyTable,col = name,sep="-",into = c("name"))
    ploidyTable$name<- as.numeric(ploidyTable$name)
    ploidyTable =  ploidyTable[ploidyTable$name!=ploidy,]
    ploidyTable$type =  ifelse(ploidyTable$name < ploidy, "loss", ifelse((ploidyTable$name > ploidy), "gain", "unclear"))
    ## REPEAT FOR WHOLE DF
    tmp=  separate(tmp,col = name,sep="-",into = c("name"))
    tmp$name<- as.numeric(tmp$name)
    tmp =  tmp[tmp$name!=ploidy,]
    tmp$type =  ifelse(tmp$name < ploidy, "loss", ifelse((tmp$name > ploidy), "gain", "unclear"))
    
    
    ## QUANTIFY GAINS/LOSSES
    totalGain= sum(filter(ploidyTable,type=="gain")$sum)
    totalGainSeg = sum(filter(ploidyTable,type=="gain")$`n()`)
    totalLoss= sum(filter(ploidyTable,type=="loss")$sum)
    totalLossSeg = sum(filter(ploidyTable,type=="loss")$`n()`)
    
    
    #SUMMARY STATS
    row <- data.frame(ID=id,ploidy=ploidy,gainSeg=totalGainSeg,totalGain=totalGain,lossSeg=totalLossSeg,totalLoss=totalLoss,file=file)
    cnvPerCellSummary <- rbind(row,cnvPerCellSummary)
    
    #SUMMARIZE INDIVIDUAL EVENTS
    df = tmp %>% select(c(seqnames,start,end,width,name,type))
    if (nrow(df)>0){
      df$file=file
      df$ploidy=ploidy
      df$gene=id
      cnvPerCell<- rbind(df,cnvPerCell)
    }
  }
  
  write.table(cnvPerCell,file = paste0("CNVs/",cutoff/10000,"cnvPerCell.txt"),quote = F,row.names = F,col.names = T,sep = "\t")
  write.table(cnvPerCellSummary,paste0("CNVs/",cutoff/10000,"cnvPerCellSummary.txt"),quote = F,row.names = F,col.names = T,sep = "\t")
  
}

copyNumberFinder(30000000,file)
copyNumberFinder(20000000,"../CoNoVariants/Input/Thesis/binsize_1e+05_stepsize_1e+05_StrandSeq_CNV.bed.gz")
copyNumberFinder(10000000,"../CoNoVariants/Input/Thesis/binsize_1e+05_stepsize_1e+05_StrandSeq_CNV.bed.gz")
copyNumberFinder(7500000,"../CoNoVariants/Input/Thesis/binsize_1e+05_stepsize_1e+05_StrandSeq_CNV.bed.gz")
copyNumberFinder(5000000,"../CoNoVariants/Input/Thesis/binsize_1e+05_stepsize_1e+05_StrandSeq_CNV.bed.gz")
copyNumberFinder(2500000,"../CoNoVariants/Input/Thesis/binsize_1e+05_stepsize_1e+05_StrandSeq_CNV.bed.gz")
copyNumberFinder(1000000,"../CoNoVariants/Input/Thesis/binsize_1e+05_stepsize_1e+05_StrandSeq_CNV.bed.gz")



for (cutoff in c(3000,2000,1000,750,500,250,100)){
  print(cutoff)
  cnvPerCell <- read.table(paste0("CNVs/",cutoff,"cnvPerCell.txt"),header=T)
  cnvPerCellSummary <- read.table(paste0("CNVs/",cutoff,"cnvPerCellSummary.txt"),header=T)
  
  cnvPerCell$type<-as.factor(cnvPerCell$type)
  gains <- dplyr::filter(cnvPerCell,cnvPerCell$type=="gain")
  gains$file<-as.factor(gains$file)
  
  cnvpath=paste0("CNVs/CNV_",cutoff,"kb/")
  if (!file.exists(cnvpath) ) { dir.create(cnvpath)}
  
  for (bamfile in levels(droplevels(gains$file))){
    print(bamfile)
    filepath=paste0(cnvpath,"/",bamfile,"/")
    if (!file.exists(filepath) ) { dir.create(filepath)}
    
    breaks = dplyr::filter(gains,gains$file==bamfile)
    breakpointSegments = data.frame(seqnames=c(),start=c(),end=c())
    
    for (rows in 1:nrow(breaks)){
      row = breaks[rows,1:4]
      rowA = data.frame(seqnames=as.character(row$seqnames),start=as.numeric(row$start),end=(as.numeric(row$start)+1))
      rowB = data.frame(seqnames=as.character(row$seqnames),start=as.numeric(row$end),end=(as.numeric(row$end)+1))
      
      breakpointSegments =  rbind(breakpointSegments,rowA)
      breakpointSegments =  rbind(breakpointSegments,rowB)
    }
    
    tmp = get(load(paste0("../NA878/data/",bamfile,".RData")))
    
    seqinfo <- tmp$breaks@seqinfo
    seqnameLevels <- levels(tmp$breaks@seqnames)
    tmp$breaks = GRanges(breakpointSegments)
    tmp$breaks@seqinfo <- seqinfo
    levels(tmp$breaks@seqnames) <- seqnameLevels
    
    tmp$confint = tmp$breaks
    tmp$confint@seqinfo <- seqinfo
    levels(tmp$confint@seqnames) <- seqnameLevels
    
    save(tmp, file=paste0(filepath,bamfile,"_bl.RData"))
    #chr=levels(droplevels(as.factor(breaks$seqnames)))[1]
    
    
    for (chr in levels(droplevels(as.factor(breaks$seqnames)))){
      breakpointR::plotBreakpointsPerChr(files2plot = paste0(filepath,bamfile,"_bl.RData"), plotspath = filepath,chromosomes = chr)
    }
  }
  
}




