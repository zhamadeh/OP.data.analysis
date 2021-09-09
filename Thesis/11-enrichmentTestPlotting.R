#This script is for plotting out enrichment test results
suppressMessages(suppressWarnings(library(tidyverse)))


dataframe = read.table("Thesis/Enrichment/GenomePermute/na878_sces.bed_GRCh37_G4_1-3.bed_1000_permutations_0bp_flank.txt",header=T)
dataframe = read.table("Thesis/Enrichment/GenomePermute/na878_sces.bed_exons.bed_1000_permutations_0bp_flank.txt",header=T)
dataframe = read.table("Thesis/Enrichment/GenomePermute/na878_sces.bed",header=T)
dataframe = read.table("Thesis/Enrichment/GenomePermute/na878_sces.bed_G4_K_plus.txt_1000_permutations_0bp_flank.txt",header=T)
dataframe = read.table("Thesis/Enrichment/GenomePermute/na878_sces.bed_introngs.bed_1000_permutations_0bp_flank.txt",header=T)


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



