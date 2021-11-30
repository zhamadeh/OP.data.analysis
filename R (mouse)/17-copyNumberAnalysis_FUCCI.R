
######################################################################
####   Packages   ####
######################################################################

library(rtracklayer)
library(tidyverse)
library(dplyr)
library(plyr)
library(tidyr)
require(GenomicRanges)
require(GenomicAlignments)

######################################################################
####   Assemble datasets   ####
######################################################################


copyNumberFinder(30000000,"CNVs/FUCCI/Data-AneuFinder/BLM-RECQL5/dnacopy/binsize_1e+05_stepsize_1e+05_StrandSeq_CNV.bed.gz")

copyNumberFinder <- function(cutoff, file){

  CNV<- import(file)
  cnvPerCell<-data.frame()
  cnvPerCellSummary=data.frame()
  
  for (i in 1:length(CNV)){
    
    ## ASSIGN GENOTYPE BASED OFF FILE NAME
    ID <- strsplit(strsplit(CNV[i]@listData[[1]]@trackLine@description,split = " ")[[1]][4],split = "[-_.]")[[1]]
    file=strsplit(CNV[i]@listData[[1]]@trackLine@description,split = " ")[[1]][4]
    
    message("Reading file: ",file," ... ",round((i/length(CNV))*100,2),"%")
    
    ## GENOTYPE BRUTE FORCE
    for (j in ID){
      if ("blm" %in% tolower(ID) ){
        if ("recq5" %in% tolower(ID) ||"recql5" %in% tolower(ID) ){
          id <- "BLM/RECQL5"
        } else {
          id <- "BLM"
        }
      } else if ("recq5" %in% tolower(ID) ||"recql5" %in% tolower(ID) ){
        if (! "blm" %in% tolower(ID) ){
          id <- "RECQL5"
        }
      } else {id <- "WT" }
    }
    
    
    ## RETRIEVE FIRST LIBRARY AND REMOVE SMALL SEGMENTS < 500Kb and HIGH PLOIDY STATES TO ASSIGN  PLOIDY
    tmp <- as.data.frame(CNV[i])
    tmp$name <- as.factor(tmp$name)
    tmp=tmp[tmp$width>500000,]
    tmp=tmp[tmp$seqnames!="chrY",]
    #HIGH PLOIDY STATES
    tmp <- tmp[tmp$name!="20-somy" & tmp$name!="19-somy" &tmp$name!="18-somy" &tmp$name!="17-somy" & tmp$name!="16-somy" &tmp$name!="15-somy" &tmp$name!="14-somy"&tmp$name!="13-somy" &tmp$name!="12-somy"&tmp$name!="11-somy"&tmp$name!="10-somy"&tmp$name!="zero-inflation",]
    tmp$name <-droplevels(tmp$name)
    
    
    ## BUILD PLOIDY TABLE SUMMARY TO CLASSIFY NATIVE STATE AND GAINS & LOSSES
    tmp$width <- as.numeric(tmp$width)
    ploidyTable <- tmp %>% group_by(name) %>% dplyr::summarize(n(),sum=sum(width))
    ploidy = as.numeric(strsplit(as.character(ploidyTable[which.max(ploidyTable$sum),]$name),split = "[-]")[[1]][1])
    if (ploidy==0){
      ploidy=1
    } else if (ploidy >2){
      next
    }
    
    tmp
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
  
  write.table(cnvPerCell,file = paste0("CNVs/FUCCI/",cutoff/1000000,"Mb_cnvPerCell.txt"),quote = F,row.names = F,col.names = T,sep = "\t")
  write.table(cnvPerCellSummary,paste0("CNVs/FUCCI/",cutoff/1000000,"Mb_cnvPerCellSummary.txt"),quote = F,row.names = F,col.names = T,sep = "\t")
  
}



######################################################################
               ####   Data   ####
######################################################################


cnvPerCell <- read.table(paste0("CNVs/FUCCI/30Mb_cnvPerCell.txt"),header=T)
cnvPerCellSummary <- read.table(paste0("CNVs/FUCCI/30Mb_cnvPerCellSummary.txt"),header=T)
cnvPerCell$type<-as.factor(cnvPerCell$type)



######################################################################
                ####   Duplications   ####
######################################################################


gains <- dplyr::filter(cnvPerCell,cnvPerCell$type=="gain")
gains$file<-as.factor(gains$file)
  
cnvpath=paste0("CNVs/FUCCI/Duplications/")
if (!file.exists(cnvpath) ) { dir.create(cnvpath)}


for (chr in levels(droplevels(gains$seqnames)) ){
    filepath=paste0(cnvpath,"/",chr,"/")
    if (!file.exists(filepath) ) { dir.create(filepath)}
    
    breaks = dplyr::filter(gains,seqnames==chr)
    
    breakpointR::plotBreakpointsPerChr(files2plot = paste0("SCEs/BLM-RECQL5/Data-BreakpointR/data/",levels(droplevels(breaks$file)),".RData"), plotspath = filepath,chromosomes = chr)
   
    karpath=paste0(filepath,"/karyogram/")
    if (!file.exists(karpath) ) { dir.create(karpath)}
    
    for (file in list.files("CNVs/FUCCI/Data-AneuFinder/BLM-RECQL5/dnacopy/RData/",full.names = T)){
      
      filename = strsplit(basename(file),"bam")[[1]][1]
      
      if (paste0(filename,"bam") %in% levels(droplevels(breaks$file))){
        plot(file,type="karyogram",plot.breakpoints=F,both.strands=F) +ggsave(paste0(karpath,basename(file),".png"))
      }
    }
    
}

######################################################################
                  ####   Deletions  ####
######################################################################


deletions <- dplyr::filter(cnvPerCell,cnvPerCell$type=="loss")
deletions$file<-as.factor(deletions$file)

cnvpath=paste0("CNVs/FUCCI/Deletions/")
if (!file.exists(cnvpath) ) { dir.create(cnvpath)}


for (chr in levels(droplevels(deletions$seqnames)) ){
  filepath=paste0(cnvpath,"/",chr,"/")
  if (!file.exists(filepath) ) { dir.create(filepath)}
  
  breaks = dplyr::filter(deletions,seqnames==chr)
  
  breakpointR::plotBreakpointsPerChr(files2plot = paste0("SCEs/BLM-RECQL5/Data-BreakpointR/data/",levels(droplevels(breaks$file)),".RData"), plotspath = filepath,chromosomes = chr)
  
  karpath=paste0(filepath,"/karyogram/")
  if (!file.exists(karpath) ) { dir.create(karpath)}
  
  for (file in list.files("CNVs/FUCCI/Data-AneuFinder/BLM-RECQL5/dnacopy/RData/",full.names = T)){
    
    filename = strsplit(basename(file),"bam")[[1]][1]
    
    if (paste0(filename,"bam") %in% levels(droplevels(breaks$file))){
      plot(file,type="karyogram",plot.breakpoints=F,both.strands=F) +ggsave(paste0(karpath,basename(file),".png"))
    }
  }
  
}

######################################################################
          ####   Manual curation of CNVs  ####
######################################################################
keep = data.frame()
remove=data.frame()

#### CHROMOSOME 2 ####
two = gains[(gains$seqnames=="chr2"),]
keep = rbind(keep,two[(two$seqnames=="chr2") & (two$file=="blm-recq5-2-single-f16-r5-c5_S572_.trimmed.mdup.bam"),])
remove =rbind(remove,(two %>% anti_join(keep)) )

#### CHROMOSOME 3 ####
three = gains[(gains$seqnames=="chr3"),]
keep = rbind(keep,three[(three$seqnames=="chr3") & (three$file=="blm-recq5-1-single-f15-r6-c3_S528_.trimmed.mdup.bam"),])
remove =rbind(remove,(three %>% anti_join(keep)) )

#### CHROMOSOME 4 ####
remove =rbind(remove, gains[(gains$seqnames=="chr4") ,])

#### CHROMOSOME 5 ####
five = gains[(gains$seqnames=="chr5"),]
keep = rbind(keep,five[(five$seqnames=="chr5") & (five$file=="blm-recq5-1-single-f14-r6-c6_S482_.trimmed.mdup.bam" | five$file=="blm_recq5-2-single-c16-w36_.trimmed.mdup.bam" ),])
remove =rbind(remove, (five %>% anti_join(keep))  )

#### CHROMOSOME 6 ####
six = gains[(gains$seqnames=="chr6"),]
keep = rbind(keep,six[(six$seqnames=="chr6") & (six$file=="blm_recq5-2-single-c16-w36_.trimmed.mdup.bam" | six$file=="blm-recq5-2-single-f16-r2-c4_S550_.trimmed.mdup.bam" ),])
remove =rbind(remove, (six %>% anti_join(keep))  )

#### CHROMOSOME 7 ####
seven = gains[(gains$seqnames=="chr7"),]
keep = rbind(keep,seven[(seven$seqnames=="chr7") & (seven$file=="blm-1-single-c10-w33_.trimmed.mdup.bam" | seven$file=="blm-recq5-1-single-f13-r7-c1_S435_.trimmed.mdup.bam" ),])
remove =rbind(remove, (seven %>% anti_join(keep))  )

#### CHROMOSOME 8 ####
eight = gains[(gains$seqnames=="chr8"),]
remove = rbind(remove,eight[(eight$seqnames=="chr8") & (eight$file=="recq5-2-single-f11-r5-c5_S327_.trimmed.mdup.bam" ),])
keep =rbind(keep, (eight %>% anti_join(remove))  )

#### CHROMOSOME 9 ####
nine = gains[(gains$seqnames=="chr9"),]
keep = rbind(keep,nine[(nine$seqnames=="chr9") & (nine$file=="blm-1-single-f7-r6-c2_S135_.trimmed.mdup.bam" | nine$file=="blm-recq5-2-single-f8-r7-c2_S191_.trimmed.mdup.bam" ),])
remove =rbind(remove, (nine %>% anti_join(keep))  )

#### CHROMOSOME 10 ####
ten = gains[(gains$seqnames=="chr10"),]
remove =rbind(remove, ten )

#### CHROMOSOME 11 ####
eleven = gains[(gains$seqnames=="chr11"),]
keep = rbind(keep,eleven[(eleven$seqnames=="chr11") & (eleven$file=="blm-1-single-f7-r6-c2_S135_.trimmed.mdup.bam" | eleven$file=="blm-recq5-2-single-f8-r7-c2_S191_.trimmed.mdup.bam" ),])
remove = rbind(remove, filter(eleven, file %in% c("zeid_RECQ5_UV300s_hoechst3_cluster15_well27_S517_.trimmed.mdup.bam","zeid_RECQ5_UV300s_hoechst3_cluster15_well31_S521_.trimmed.mdup.bam" )))
keep =rbind(keep, (eleven %>% anti_join(remove))  )

#### CHROMOSOME 12 ####
twelve = gains[(gains$seqnames=="chr12"),]
keep = rbind(keep, filter(twelve, file %in% c("recq5-2-single-f11-r5-c5_S327_.trimmed.mdup.bam")))
remove =rbind(remove, (twelve %>% anti_join(keep))  )

#### CHROMOSOME 13 ####
thirteen = gains[(gains$seqnames=="chr13"),]
keep = rbind(keep,thirteen[(thirteen$seqnames=="chr13") & (thirteen$file=="blm-1-single-c10-w33_.trimmed.mdup.bam" | thirteen$file=="blm_recq5-2-single-c15-w21_.trimmed.mdup.bam" ),])
remove =rbind(remove, (thirteen %>% anti_join(keep))  )

#### CHROMOSOME 14 ####
fourteen = gains[(gains$seqnames=="chr14"),]
remove =rbind(remove, fourteen )

#### CHROMOSOME 15 ####
fifteen = gains[(gains$seqnames=="chr15"),]
remove =rbind(remove, fifteen )

#### CHROMOSOME 16 ####
sixteen = gains[(gains$seqnames=="chr16"),]
keep = rbind(keep,sixteen[(sixteen$seqnames=="chr16") & (sixteen$file=="blm-recq5-2-single-f16-r4-c2_S562_.trimmed.mdup.bam"  ),])
remove =rbind(remove, (sixteen %>% anti_join(keep))  )

#### CHROMOSOME 17 ####
seventeen = gains[(gains$seqnames=="chr17"),]
remove =rbind(remove, seventeen )

#### CHROMOSOME 18 ####
eighteen = gains[(gains$seqnames=="chr18"),]
keep = rbind(keep,eighteen[(eighteen$seqnames=="chr18") & (eighteen$file=="zeid_RECQ5_UV300s_hoechst3_cluster15_well13_S503_.trimmed.mdup.bam" | eighteen$file=="blm-1-single-c5-w42-3d6_S142_.trimmed.mdup.bam" ),])
remove =rbind(remove, (eighteen %>% anti_join(keep))  )

#### CHROMOSOME 19 ####
nineteen = gains[(gains$seqnames=="chr19"),]
keep = rbind(keep,nineteen[(nineteen$seqnames=="chr19") & (nineteen$file=="blm-recq5-1-single-c11-w1-3d6_S395_.trimmed.mdup.bam" ),])
remove =rbind(remove, (nineteen %>% anti_join(keep))  )

#### CHROMOSOME 20 ####
twenty = gains[(gains$seqnames=="chr20"),]
remove =rbind(remove,twenty )

#### CHROMOSOME 21 ####
twentyone = gains[(gains$seqnames=="chr21"),]
remove =rbind(remove, twentyone)

#### CHROMOSOME 22 ####
twentytwo = gains[(gains$seqnames=="chr22"),]
remove =rbind(remove, twentytwo)

#### CHROMOSOME X ####
x = gains[(gains$seqnames=="chrX"),]
keep = rbind(keep,x)


curatedCNVs = keep
write.table(curatedCNVs,"CNVs/FUCCI/30Mb_curatedCallsetDuplications.txt",row.names = F,col.names = T,quote = F,sep = "\t")
write.table(remove,"CNVs/FUCCI/30Mb_curatedCallsetDuplications_FALSE_POSITIVES.txt",row.names = F,col.names = T,quote = F,sep = "\t")


######################################################################
####   Manual curation of CNVs  ####
######################################################################
keep = data.frame()
remove=data.frame()


#### CHROMOSOME 1 ####
one = deletions[(deletions$seqnames=="chr1"),]
keep = rbind(keep,one[(one$seqnames=="chr1") & (one$file=="blm_recq5-2-single-c16-w16_.trimmed.mdup.bam"),])
remove =rbind(remove,(one %>% anti_join(keep)) )

#### CHROMOSOME 2 ####
two = deletions[(deletions$seqnames=="chr2"),]
remove =rbind(remove, two )

#### CHROMOSOME 3 ####
three = deletions[(deletions$seqnames=="chr3"),]
remove =rbind(remove,three )

#### CHROMOSOME 4 ####
four = deletions[(deletions$seqnames=="chr4"),]
keep = rbind(keep,four[(four$seqnames=="chr4") & (four$file=="recq5-2-single-c6-w12-3d6_S161_.trimmed.mdup.bam"),])
remove =rbind(remove,(four %>% anti_join(keep)) )

#### CHROMOSOME 5 ####
five = deletions[(deletions$seqnames=="chr5"),]
keep = rbind(keep,five[(five$seqnames=="chr5") & (five$file=="blm-1-single-f7-r1-c1_S99_.trimmed.mdup.bam" ),])
remove =rbind(remove, (five %>% anti_join(keep))  )

#### CHROMOSOME 6 ####
six = deletions[(deletions$seqnames=="chr6"),]
remove =rbind(remove,six )

#### CHROMOSOME 7 ####
seven = deletions[(deletions$seqnames=="chr7"),]
keep = rbind(keep,seven)

#### CHROMOSOME 8 ####
eight = deletions[(deletions$seqnames=="chr8"),]
remove = rbind(remove,filter(eight, file %in% c("blm-1-single-c5-w24-3d6_S124_.trimmed.mdup.bam"," blm-recq5-2-single-f12-r1-c2_S345_.trimmed.mdup.bam","blm-recq5-2-single-c4-w19-3d6_S85_.trimmed.mdup.bam",
                                                "zeid_BLM_UV300s_hoechst3_cluster14_well33_S474_.trimmed.mdup.bam","blm-recq5-2-single-f8-r6-c1_S183_.trimmed.mdup.bam") ) )
keep =rbind(keep, (eight %>% anti_join(remove))  )

#### CHROMOSOME 9 ####
nine = deletions[(deletions$seqnames=="chr9"),]
remove =rbind(remove, nine   )

#### CHROMOSOME 10 ####
ten = deletions[(deletions$seqnames=="chr10"),]
keep =rbind(keep, ten )

#### CHROMOSOME 11 ####
eleven = deletions[(deletions$seqnames=="chr11"),]
keep =rbind(keep, eleven )

#### CHROMOSOME 12 ####
twelve = deletions[(deletions$seqnames=="chr12"),]
keep = rbind(keep, twelve   )

#### CHROMOSOME 13 ####
thirteen = deletions[(deletions$seqnames=="chr13"),]
keep =rbind(keep, thirteen )

#### CHROMOSOME 14 ####
fourteen = deletions[(deletions$seqnames=="chr14"),]
remove =rbind(remove, fourteen )

#### CHROMOSOME 15 ####
fifteen = deletions[(deletions$seqnames=="chr15"),]
remove =rbind(remove, fifteen )

#### CHROMOSOME 16 ####
sixteen = deletions[(deletions$seqnames=="chr16"),]
keep = rbind(keep,sixteen[(sixteen$seqnames=="chr16") & (sixteen$file=="blm-recq5-2-single-f16-r1-c3_S542_.trimmed.mdup.bam"  ),])
remove =rbind(remove, (sixteen %>% anti_join(keep))  )

#### CHROMOSOME 17 ####
seventeen = deletions[(deletions$seqnames=="chr17"),]
keep =rbind(keep, seventeen )

#### CHROMOSOME 18 ####
eighteen = deletions[(deletions$seqnames=="chr18"),]
remove =rbind(remove, eighteen   )

#### CHROMOSOME 19 ####
nineteen = deletions[(deletions$seqnames=="chr19"),]
remove =rbind(remove, nineteen )

#### CHROMOSOME 20 ####
twenty = deletions[(deletions$seqnames=="chr20"),]
remove =rbind(remove,twenty )

#### CHROMOSOME 21 ####
twentyone = deletions[(deletions$seqnames=="chr21"),]
remove =rbind(remove, twentyone)

#### CHROMOSOME 22 ####
twentytwo = deletions[(deletions$seqnames=="chr22"),]
remove =rbind(remove, twentytwo)

#### CHROMOSOME X ####
x = deletions[(deletions$seqnames=="chrX"),]
remove = rbind(remove,x)


curatedDeletions = keep
curatedDeletionFalsePositives = remove
write.table(curatedDeletions,"CNVs/FUCCI/30Mb_curatedCallsetDeletions.txt",row.names = F,col.names = T,quote = F,sep = "\t")
write.table(curatedDeletionFalsePositives,"CNVs/FUCCI/30Mb_curatedCallsetDeletions_FALSE_POSITIVES.txt",row.names = F,col.names = T,quote = F,sep = "\t")


