#This script is for plotting out enrichment test results
suppressMessages(suppressWarnings(library(ggplot2)))


dataframe = read.table("Thesis/Enrichment/GenomePermute/na878_sces.bed_GRCh37_G4_1-3.bed_1000_permutations_0bp_flank.txt",header=T)
dataframe = read.table("Thesis/Enrichment/GenomePermute/na878_sces.bed_exons.bed_1000_permutations_0bp_flank.txt",header=T)
#dataframe = read.table("Thesis/Enrichment/GenomePermute/na878_sces.bed_GRCh37_ens75_genes.bed_1000_permutations_0bp_flank.txt",header=T)
dataframe = read.table("Thesis/Enrichment/GenomePermute/na878_sces.bed_G4_K.bed_1000_permutations_0bp_flank.txt",header=T)
dataframe = read.table("Thesis/Enrichment/GenomePermute/na878_sces.bed_introngs.bed_1000_permutations_0bp_flank.txt",header=T)
#dataframe = read.table("Thesis/Enrichment/GenomePermute/na878_sces.bed_hg38_cfs.bed_1000_permutations_0bp_flank.txt",header=T)


#binding all 4 dataframes together
true <- dataframe[1,]
bind<- dataframe[2:1001,]
mean = median(bind$Overlaps)
true$enrichment = true$Overlaps/mean
bind$enrichment = bind$Overlaps/mean


#plotting
ggplot(bind)+geom_density(aes(enrichment)) +
	theme_classic(base_size = 19) +
  geom_vline(mapping = aes(xintercept=true$enrichment[1]))+
	theme(legend.position = "none")
        

gaps <- read.table("Thesis/Enrichment/GenomePermute/exons2.bed") %>% select(c(V1,V2,V3))
 gaps$V1<-as.character(gaps$V1)
#gaps=separate(gaps,V1,into=c("chr","start","end") ,"[:-]")
gaps$V1=gsub("x","X",gaps$V1)

write.table(gaps,"Thesis/Enrichment/GenomePermute/exons2.bed",sep = "\t",quote=F,row.names = F,col.names = F)
