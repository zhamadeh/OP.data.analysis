#######################################################
####     Workflow Part 6: Complete Dataset          ###
#######################################################

qualityFilterLibraries <- function(){

	metrics <- read.table("INPUT/SUMMARY/MOUSE/05.library.quality.txt",header  =T) 
	metrics$Library <- as.factor(metrics$Library)
	#metrics%>% group_by(quality)%>% summarize(n())
	
	breakpoints<- read.table("INPUT/SUMMARY/MOUSE/03.library.breakpoints.txt",header=T)
  breakpoints$Library<-as.factor(breakpoints$library)
  
  #breakpoints$Library<-paste0(breakpoints$Library,".processed.bam")
  breakpoints_full = merge(breakpoints,metrics,by="Library")
  
  #write.table(breakpoints_full,"INPUT/05.breakpoints.complete.txt",quote = F,col.names = T,row.names = F,sep = "\t")
  
  #breakpoints_full$library <- breakpoints_full$Library
	suppressWarnings(breakpoints_full <- breakpoints_full %>% separate(library, c("a","b","c","d","e","f"), "[_-]+"))
	breakpoints_full$gene <- "gene"
	
	breakpoints_full[breakpoints_full$b=="c38",]$gene = "BLM"
	breakpoints_full[breakpoints_full$b=="wt",]$gene = "WT"
	
	
	breakpoints_full <- select(breakpoints_full,-c(a,b,c,d,e,f))
	write.table(breakpoints_full,"INPUT/SUMMARY/MOUSE/06.gene.full.txt",quote = F,col.names = T,row.names = F,sep = "\t")

}
