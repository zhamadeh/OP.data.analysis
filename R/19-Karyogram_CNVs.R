
######################################################################
                ####   Plotting Karyogram   ####
######################################################################

################ Step 1: load genomic range data ####################

# Needs to be in dataframe format, can be points or ranges
SNP=select(curatedDeletions,c(seqnames,start,end))
#cnvPerCell <- read.table(paste0("CNVs/FUCCI/30Mb_cnvPerCell.txt"),header=T)
#SNP=filter(cnvPerCell,type=='gain') %>% select(c(seqnames,start,end))

# Give bed file column names for first three rows
colnames(SNP)=c("chr","start","end")

# Convert chrom to factor
SNP$chr<-as.factor(SNP$chr)

# Reorder levels for printing 
SNP$chr = factor(SNP$chr,levels=c("chr1" ,"chr2" ,"chr3" , "chr4"  ,"chr5" , "chr6"  ,"chr7"  ,"chr8",  "chr9" ,
                                  "chr10", "chr11" ,"chr12", "chr13" ,"chr14" ,"chr15" ,"chr16" ,"chr17", "chr18" ,"chr19" , 
                                  "chr20", "chr21","chr22" , "chrX" , "chrY" ))

# Divide by 1Mb to simplify axis in plotting
SNP$start=SNP$start/1000000
SNP$end = SNP$end/1000000


################ Step 2: load Genome data (lengths) ####################

# Load chromosome lengths
data = read.table("Anxilliary/hg38_chrLengths") 

# Make sure they have "chr" nomenclature
data$V1<-paste0("chr",data$V1)

# Get rid of 'start' column you only need the lengths
data=select(data,-c(V2))

# Rename columns to chromosome and size
colnames(data)=c("chr","size")

# Change to factor
data$chr<-as.factor(data$chr)

# Reorder levels for plotting
data$chr = factor(data$chr,levels=c("chr1" ,"chr2" ,"chr3" , "chr4"  ,"chr5" , "chr6"  ,"chr7"  ,"chr8",  "chr9" ,
                     "chr10", "chr11" ,"chr12", "chr13" ,"chr14" ,"chr15" ,"chr16" ,"chr17", "chr18" ,"chr19" , 
                     "chr20", "chr21","chr22" , "chrX" , "chrY" ))

# Divide by 1Mb to clean up axis
data$size=data$size/1000000


################ Step 3: load secondary dataset (eg. CFSs) ####################

# Load data
cfs = read.table("Anxilliary/Fragile_sites/hg38_cfs.bed")

# Change lowercase 'x'
cfs$V1=gsub("chrx","chrX",cfs$V1)

# Divide by 1Mb for axis
cfs$V3=cfs$V3/1000000
cfs$V2=cfs$V2/1000000

# Change to factor and reorder levels
cfs$V1<-as.factor(cfs$V1)
cfs$V1 = factor(cfs$V1,levels=c("chr1" ,"chr2" ,"chr3" , "chr4"  ,"chr5" , "chr6"  ,"chr7"  ,"chr8",  "chr9" ,
                                  "chr10", "chr11" ,"chr12", "chr13" ,"chr14" ,"chr15" ,"chr16" ,"chr17", "chr18" ,"chr19" , 
                                  "chr20", "chr21","chr22" , "chrX" , "chrY" ))



################ Step 4: Plotting  ####################

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
  ggsave("OUTPUT/Karyograms/karyogram_curatedDeletions.png")

