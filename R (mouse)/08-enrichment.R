
chromosomeLengths <- read.table("INPUT/SUMMARY/MOUSE/chromosomeLength_names",header=T)
chromsomeNames = select(chromosomeLengths,c(Chromosome,Accession))
mouse.df = merge(mouse.df,chromsomeNames,by.x="seqnames",by.y="Accession") 
mouse.df = mouse.df %>% select(-c(seqnames))

blm = filter(mouse.df,gene=="BLM")
wt = filter(mouse.df,gene=="WT")

write.table(select(blm,c(Chromosome, start,end)),"R (mouse)/Enrichment/blm.sces.bed",quote=F,row.names = F,col.names = F,sep = "\t")
write.table(select(wt,c(Chromosome, start,end)),"R (mouse)/Enrichment/wt.sces.bed",quote=F,row.names = F,col.names = F,sep = "\t")

chrLen = select(chromosomeLengths,c(Chromosome,Length))
write.table(chrLen,"R (mouse)/Enrichment/GRCh37_chrlen.txt",quote=F,row.names = F,col.names = F,sep = "\t")
write.table(data.frame(),"R (mouse)/Enrichment/GRCh37_gaps.bed",quote=F,row.names = F,col.names = F,sep = "\t")

mouse_g4= read.table("Supplemental_File_S2.csv",sep=",",header=T)
mouse_g4<-as.data.frame(mouse_g4)
mouse_g4=filter(mouse_g4,mm10 == "mm10")
mouse_g4=select(mouse_g4,c(mm10.chr,mm10.start,mm10.end))
mouse_g4$mm10.chr=gsub("chr","",mouse_g4$mm10.chr)

write.table(mouse_g4,"R (mouse)/Enrichment/g4s.bed",quote=F,row.names = F,col.names = F,sep = "\t")


#Loading the data
blm=read.table("R (mouse)/Enrichment/blm.sces.bed_g4s.bed_1000_permutations_0bp_flank.txt",header=T)
blm$model="blm"
blm$enrichment = blm$Overlaps/mean(blm$Overlaps[2:nrow(blm)])
wt=read.table("R (mouse)/Enrichment/wt.sces.bed_g4s.bed_1000_permutations_0bp_flank.txt",header=T)
wt$model="wt"
wt$enrichment = wt$Overlaps/mean(wt$Overlaps[2:nrow(wt)])

#binding all 4 dataframes together
true <- rbind(blm[1,],wt[1,])
bind<- rbind(blm[2:nrow(blm),],wt[2:nrow(wt),])

pval = data.frame(model=c("blm","wt"),p = c("<0.001","0.051"),ast = c("***",""))

#plotting
ggplot(bind)+geom_violin(aes(model,enrichment,fill=model),lwd=0.9) +
  geom_boxplot(width=0.05,aes(model,enrichment),lwd=0.9) +
  #scale_fill_manual(values=c("#53b0db","#53b0db","#53b0db","#6ceb70"))+
  theme_classic(base_size = 19) +
  geom_point(data=true,aes(model,enrichment),fill="red",colour="black",pch=21, size=5)+
  theme(legend.position = "none")  +
  labs(x="Cells",y="Enrichment") +
  geom_text(pval, mapping=aes(x=model,y=1.32,label=p))+
  geom_text(pval, mapping=aes(x=model,y=1.29,label=ast))
  #ggsave(paste0("Plots/enrichmentPlot",args[1],".png"))
