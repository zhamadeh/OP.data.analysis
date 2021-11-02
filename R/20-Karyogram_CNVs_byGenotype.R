
cnvPerCell=curatedDeletions
#cnvPerCell <- read.table(paste0("CNVs/3000cnvPerCell.txt"),header=T)
cnvPerCell$gene<-as.factor(cnvPerCell$gene)


data = read.table("../hg38_chrLengths") 
data$V1<-paste0("chr",data$V1)
data=select(data,-c(V2))
colnames(data)=c("chr","size")
data$chr<-as.factor(data$chr)
data$chr = factor(data$chr,levels=c("chr1" ,"chr2" ,"chr3" , "chr4"  ,"chr5" , "chr6"  ,"chr7"  ,"chr8",  "chr9" ,
                                    "chr10", "chr11" ,"chr12", "chr13" ,"chr14" ,"chr15" ,"chr16" ,"chr17", "chr18" ,"chr19" , 
                                    "chr20", "chr21","chr22" , "chrX" , "chrY" ))

data$size=data$size/1000000


## CFS

cfs = read.table("Anxilliary/Fragile_sites/hg38_cfs.bed")
cfs$V1=gsub("chrx","chrX",cfs$V1)
cfs$V3=cfs$V3/1000000
cfs$V2=cfs$V2/1000000
cfs$V1<-as.factor(cfs$V1)
cfs$V1 = factor(cfs$V1,levels=c("chr1" ,"chr2" ,"chr3" , "chr4"  ,"chr5" , "chr6"  ,"chr7"  ,"chr8",  "chr9" ,
                                "chr10", "chr11" ,"chr12", "chr13" ,"chr14" ,"chr15" ,"chr16" ,"chr17", "chr18" ,"chr19" , 
                                "chr20", "chr21","chr22" , "chrX" , "chrY" ))




g="RECQL5"
  
  SNP=dplyr::filter(cnvPerCell,gene==g)
  SNP=SNP %>% select(c(seqnames,start,end))
  
  colnames(SNP)=c("seqnames","start","end")
  SNP$width = SNP$end-SNP$start
  
  
 
  
  
  
  #SNP = select(SNP,c(seqnames,start,end,width))
  colnames(SNP)= c( "chr" ,"start"  ,  "end"  ,"width" )
  
  SNP$chr<-as.factor(SNP$chr)
  SNP$chr = factor(SNP$chr,levels=c("chr1" ,"chr2" ,"chr3" , "chr4"  ,"chr5" , "chr6"  ,"chr7"  ,"chr8",  "chr9" ,
                                      "chr10", "chr11" ,"chr12", "chr13" ,"chr14" ,"chr15" ,"chr16" ,"chr17", "chr18" ,"chr19" , 
                                      "chr20", "chr21","chr22" , "chrX" , "chrY" ))
  SNP$start=SNP$start/1000000
  SNP$end = SNP$end/1000000
  
  
  

  
  # plotting
  
    ggplot() +
      geom_segment(data = data,
                   aes(x = chr, xend = chr, y = 0, yend = size),
                   lineend = "round", color = "lightgrey", size = 5) +
      geom_rect(data = SNP,
                   aes(xmin = as.integer(chr) - 0.25, xmax = as.integer(chr) + 0.25,
                       ymin = start, ymax = end)
                   ,fill="#1d705d",size = 0.25,alpha=0.2) +
      geom_rect(data = cfs, aes(xmin = as.integer(V1) - 0.4, xmax = as.integer(V1) - 0.3,
                                ymin = V2, ymax = V3) ,fill="black",
                size = 0.25) +
      theme_classic() +
      theme(text = element_text(size=18),axis.line=element_blank())+
      labs(x="",y="CHROMOSOME POSITION (Mb)")+
      ggsave("OUTPUT/Karyograms/Deletions_RECQL5_curated_karyogram.png")
 
    
    g="BLM/RECQL5"
    
    SNP=dplyr::filter(cnvPerCell,gene==g)
    SNP=SNP %>% select(c(seqnames,start,end))
    colnames(SNP)=c("seqnames","start","end")
    SNP$width = SNP$end-SNP$start
    
    
    
    
    
    
    #SNP = select(SNP,c(seqnames,start,end,width))
    colnames(SNP)= c( "chr" ,"start"  ,  "end"  ,"width" )
    
    SNP$chr<-as.factor(SNP$chr)
    SNP$chr = factor(SNP$chr,levels=c("chr1" ,"chr2" ,"chr3" , "chr4"  ,"chr5" , "chr6"  ,"chr7"  ,"chr8",  "chr9" ,
                                      "chr10", "chr11" ,"chr12", "chr13" ,"chr14" ,"chr15" ,"chr16" ,"chr17", "chr18" ,"chr19" , 
                                      "chr20", "chr21","chr22" , "chrX" , "chrY" ))
    SNP$start=SNP$start/1000000
    SNP$end = SNP$end/1000000
    
    
    
    
    
    # plotting
    
    ggplot() +
      geom_segment(data = data,
                   aes(x = chr, xend = chr, y = 0, yend = size),
                   lineend = "round", color = "lightgrey", size = 5) +
      geom_rect(data = SNP,
                aes(xmin = as.integer(chr) - 0.25, xmax = as.integer(chr) + 0.25,
                    ymin = start, ymax = end)
                ,fill="#1d705d",size = 0.25,alpha=0.2) +
      geom_rect(data = cfs, aes(xmin = as.integer(V1) - 0.4, xmax = as.integer(V1) - 0.3,
                                ymin = V2, ymax = V3) ,fill="black",
                size = 0.25) +
      theme_classic() +
      theme(text = element_text(size=18),axis.line=element_blank())+
      labs(x="",y="CHROMOSOME POSITION (Mb)")+
      ggsave("OUTPUT/Karyograms/Deletions_BLM-RECQL5_curated_karyogram.png")
    
    
    
    
    
    g="BLM"
    
    SNP=dplyr::filter(cnvPerCell,gene==g)
    SNP=SNP %>% select(c(seqnames,start,end))
    colnames(SNP)=c("seqnames","start","end")
    SNP$width = SNP$end-SNP$start
    
    
    
    
    
    
    #SNP = select(SNP,c(seqnames,start,end,width))
    colnames(SNP)= c( "chr" ,"start"  ,  "end"  ,"width" )
    
    SNP$chr<-as.factor(SNP$chr)
    SNP$chr = factor(SNP$chr,levels=c("chr1" ,"chr2" ,"chr3" , "chr4"  ,"chr5" , "chr6"  ,"chr7"  ,"chr8",  "chr9" ,
                                      "chr10", "chr11" ,"chr12", "chr13" ,"chr14" ,"chr15" ,"chr16" ,"chr17", "chr18" ,"chr19" , 
                                      "chr20", "chr21","chr22" , "chrX" , "chrY" ))
    SNP$start=SNP$start/1000000
    SNP$end = SNP$end/1000000
    
    
    
    
    
    # plotting
    
    ggplot() +
      geom_segment(data = data,
                   aes(x = chr, xend = chr, y = 0, yend = size),
                   lineend = "round", color = "lightgrey", size = 5) +
      geom_rect(data = SNP,
                aes(xmin = as.integer(chr) - 0.25, xmax = as.integer(chr) + 0.25,
                    ymin = start, ymax = end)
                ,fill="#1d705d",size = 0.25,alpha=0.2) +
      geom_rect(data = cfs, aes(xmin = as.integer(V1) - 0.4, xmax = as.integer(V1) - 0.3,
                                ymin = V2, ymax = V3) ,fill="black",
                size = 0.25) +
      theme_classic() +
      theme(text = element_text(size=18),axis.line=element_blank())+
      labs(x="",y="CHROMOSOME POSITION (Mb)")+
      ggsave("OUTPUT/Karyograms/Deletions_BLM_curated_karyogram.png")
    



