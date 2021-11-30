#######################################################
####         Workflow Part 8: Re-plotting BPR       ###
#######################################################
library(tidyverse)
library(plyr)
library(ggpubr)



breakpointSummary <- read.table("INPUT/SUMMARY/MOUSE/06.gene.full.txt",header=T) #%>% select(-c(CI.start,CI.end,genoT))
breakpointSummary = mouse.df
breakpointSummary$library=breakpointSummary$Library
breakpointSummary$seqnames = paste0('chr',breakpointSummary$Chromosome)




########### PLOT #1: SCEs/CHR/LIB vs LENGTH ###############
###########################################################

breakpointSummary$gene <- as.factor(breakpointSummary$gene)
breakpointSummary$library <- as.factor(breakpointSummary$library)

#sces per chroomosome
suppressMessages(all <- as.data.frame(breakpointSummary %>% group_by(seqnames) %>% dplyr::summarize(ALL=n())))
suppressMessages(b <- as.data.frame(breakpointSummary %>% filter(gene=="BLM")  %>% group_by(seqnames) %>% dplyr::summarize(BLM=n())))
suppressMessages(w <- as.data.frame(breakpointSummary %>% filter(gene=="WT")  %>% group_by(seqnames)%>% dplyr::summarize(WT=n())))
byChr <- merge(b,w,by="seqnames")
byChr <- merge(byChr,all,by="seqnames",all=T)
byChr[is.na(byChr)] <- 0
write.table(byChr,"OUTPUT/SummaryTables/MOUSE/perChrom.txt",quote=F,row.names = F,col.names = T,sep="\t")

breakpointSummary$Library=breakpointSummary$library
numOfLibsPerGene <- data.frame(gene=character(),n=numeric())
b <- (as.data.frame(breakpointSummary %>% filter(gene=="BLM")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
numOfLibsPerGene= add_row(numOfLibsPerGene, gene="BLM",n=as.numeric(length(levels(b$Library))))
w <- (as.data.frame(breakpointSummary %>% filter(gene=="WT")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
numOfLibsPerGene= add_row(numOfLibsPerGene, gene="WT",n=as.numeric(length(levels(w$Library))))
numOfLibsPerGene= add_row(numOfLibsPerGene, gene="ALL",n=as.numeric(length(levels(breakpointSummary$Library))))
write.table(numOfLibsPerGene,"OUTPUT/SummaryTables/MOUSE/numOfLibsPerGene.txt",quote=F,row.names = F,col.names = F,sep="\t")


numOfLibsPerGene=read.table("OUTPUT/SummaryTables/MOUSE/numOfLibsPerGene.txt",header=F)
lengths <- read.table("OUTPUT/SummaryTables/MOUSE/chrlen.txt",header=F)
lengths$V1<-paste0("chr",lengths$V1)
colnames(lengths)=c("Chromosome","Length")
byChr<- read.table("OUTPUT/SummaryTables/MOUSE/perChrom.txt",header=T)

for (i in c("BLM","WT","ALL")){
  print(byChr[,i])
  for (j in 1:nrow(numOfLibsPerGene)){
    print(j)
    if (numOfLibsPerGene[j,1]==i){
      #print(numOfLibsPerGene[j,2])
      byChr[,i] = byChr[,i]/ numOfLibsPerGene[j,2]
    }
  }
}

scePerChrPerGeneVsLength <- merge(lengths,byChr,by.x="Chromosome",by.y="seqnames")
tidy <- gather(scePerChrPerGeneVsLength, gene,sce_per_chr,BLM:ALL)
all <- filter(tidy,gene=="ALL")
tidy <- filter(tidy,gene!="ALL")
tidy$Length<- tidy$Length/1e+6

ggplot(tidy) +  geom_smooth(aes(Length,sce_per_chr,group=gene,color=gene),se=F,method="lm", na.rm=TRUE,size=2)+
                   theme_classic(base_size = 25) +
                   ylab("SCEs/Chr/Lib")+
                   xlab("Chrom Length")+
  theme(legend.position = c(0.25, 0.8),legend.title = element_blank())+
  scale_colour_grey()+
  geom_point(aes(Length,sce_per_chr,group=gene),color="blue", na.rm=TRUE)+
  ggsave("OUTPUT/Plots/MOUSE/chrLengthvsSCEs.png")

tidy$Length<-as.numeric(tidy$Length)
tidy$gene<-as.factor(tidy$gene)
tidy <- select(tidy,-Chromosome)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

tgc=summarySE(tidy, measurevar="sce_per_chr", groupvars=c("Length"))
tgc$Length=as.numeric(tgc$Length)

ggplot(tgc, aes(x=Length, y=sce_per_chr)) + 
  geom_errorbar(aes(ymin=sce_per_chr-se, ymax=sce_per_chr+se),width=10,size=1,color="grey") +
  geom_smooth(method="lm",se=F,size=3,color="black") +
  #geom_crossbar(aes(ymin=sce_per_chr-se, ymax=sce_per_chr+se))+
  theme_classic(base_size = 25) +
  ylab("SCEs / Chromosome")+
  xlab("Chromosome Length")+
  geom_point(size=4,color="#4075a8")+
  stat_cor(method = "pearson",p.accuracy = 0.001,size=7)

ggsave("Output/Plots/scePerChrPerGeneVsLength.png")






################## Plot #2 ##############################

suppressMessages(test<-as.data.frame(breakpointSummary %>%
                                       group_by(Library) %>%
                                       dplyr::summarize(n())))
test$gene <- "gene"
test$library <- test$Library
suppressWarnings(test <- test %>% separate(Library, c("a","b","c","d","e","f"), "[_-]+"))

test[test$b=="c38",]$gene = "BLM"
test[test$b=="wt",]$gene = "WT"


test <- select(test,-c(a,b,c,d,e,f))

test <- select(test,c("n()","gene"))
colnames(test)=c("sces","gene")




p=data.frame(gene=c("BLM","WT"),p=c("<0.001",""),ast = c("***",""))

ggplot(test) + geom_jitter(aes(gene,sces, color=gene))+ geom_boxplot(aes(gene,sces),width=0.1,coef = 5) +
  theme_classic()+
  theme(text=element_text(size=15),legend.position = "none")+
  labs(x="Knockout Cells",y="SCEs Per Library") +
  geom_text(p,mapping=aes(gene,85,label=p))+
  geom_text(p,mapping=aes(gene,80,label=ast))+
  ggsave("OUTPUT/Plots/MOUSE/SCEperGene.png")



my_comparisons <- list( c("WT", "RECQL5"), c("WT", "BLM/RECQL5"), c("WT", "BLM") )
ggboxplot(test, x = "gene", y = "sces",
          fill="#4075a8",  add = "jitter",width=0.1, add.params = list(color = "gene"),
		  xlab="Gene",ylab="SCEs/library") +
	stat_compare_means(comparisons = my_comparisons,label.y = c(31, 34, 37)) +
	stat_compare_means(label = "p.signif", method = "t.test",
				   ref.group = "WT")+
  scale_colour_grey()



breakpointSummary$width=breakpointSummary$end- breakpointSummary$start

suppressMessages(suppressWarnings(ggplot(breakpointSummary) + geom_smooth(aes(Reads_per_Mb, width,color=gene),se=F) +
                                    scale_y_log10() +
                                    theme_classic() +
                                    theme(text=element_text(size=15))+
                                    geom_hline(yintercept=10000, linetype="dashed", color = "red") +
                                    ggsave("OUTPUT/Plots/MOUSE/resolutionVsDepth.png")))


sce_summary=data.frame(gene=character(),SCE=numeric(),mean_resolution=numeric(),median_resolution=numeric())
b<- filter(breakpointSummary, gene=="BLM")
sce_summary<- add_row(sce_summary,gene="BLM",SCE=nrow(b),mean_resolution=mean(b$width),median_resolution=median(b$width))

w<- filter(breakpointSummary, gene=="WT")
sce_summary<- add_row(sce_summary,gene="WT",SCE=nrow(w),mean_resolution=mean(w$width),median_resolution=median(w$width))
write.table(sce_summary,"OUTPUT//SummaryTables/MOUSE/SCE_summary.txt",quote=F,row.names = F,col.names = T,sep="\t")




#plot here
ggplot(breakpointSummary) + stat_ecdf(aes(width,color=gene)) +
                                    scale_x_log10() +
                                    theme_classic() +
                                    ylab("SCEs Mapped (%)") +
                                    xlab("Resolution") +
                                    annotation_logticks(sides = "b") +
                                    theme(text = element_text(size=15))+
                                    geom_density(aes(width),size=1.1)+
                                    geom_vline(xintercept=median(breakpointSummary$width), linetype="dashed", color = "red") +
                                    geom_text(aes(x=5000, label=paste0("Median\n",median(width)," bp"), y=0.8))
                                    ggsave("OUTPUT/Plots/MOUSE/breakpointResolution.png")
                                    