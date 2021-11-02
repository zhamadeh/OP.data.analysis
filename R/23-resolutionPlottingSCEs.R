
na = read.table("SCEs/NA878/na878_sces.bed",header=F)
na$width=na$V3-na$V2
colnames(na)=c("seqnames","start","end","width")
na$cell="NA12878"
blm=read.table("INPUT/06.sisterChromatidExchanges.txt",header=T) %>% select(c(seqnames,start,end,width))
blm$cell="HAP1"
merge=rbind(blm,na)

ggplot(merge) + stat_ecdf(aes(width,color=cell),size=0.4) +
  scale_x_log10() +
  theme_classic() +
  ylab("SCEs Mapped (%)") +
  xlab("Resolution") +
  annotation_logticks(sides = "b") +
  theme(text = element_text(size=15),legend.title = element_blank(),
        legend.position=c(0.3,0.5))+
  geom_density(aes(width),size=1)#+
  #scale_x_discrete(breaks=sc)
  #scale_x_discrete(labels=c("<10 bp","<100 bp", "<1 Kbp","<10 Kbp","<100 Kbp","<1 Mbp","<10 Mbp"),
                   #breaks = c(10,100,1000,10000,100000,1000000,10000000))
