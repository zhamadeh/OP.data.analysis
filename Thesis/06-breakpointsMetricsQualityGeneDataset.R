#######################################################
####     Workflow Part 6: Complete Dataset          ###
#######################################################

qualityFilterLibraries <- function(){

	metrics <- read.table("Input/02.library.quality.txt",header  =T) 
	metrics$file <- as.factor(metrics$file)
	
	breakpoints<- read.table("Input/03.library.breakpoints.txt",header=T)
  breakpoints$ID<-as.factor(breakpoints$ID)
  
  breakpoints_full = merge(breakpoints,metrics,by.y="file",by.x="ID")
  
  write.table(breakpoints_full,"Input/04.breakpoints.metrics.quality.txt",quote = F,col.names = T,row.names = F,sep = "\t")
  
  breakpoints_full$library <- breakpoints_full$ID
	suppressWarnings(breakpoints_full <- breakpoints_full %>% separate(ID, c("a","b","c","d","e","f"), "[_-]+"))
	breakpoints_full$gene <- "gene"
	for (row in 1:nrow(breakpoints_full)){

		for (letter in c("a","b","c","d","e","f")){
			if (is.na(breakpoints_full[row,letter])!=T){
				if (breakpoints_full[row,letter]=="WT" | breakpoints_full[row,letter]=="wt"){
					breakpoints_full[row,"gene"]="WT"
				}
				else if (breakpoints_full[row,letter]=="blm" | breakpoints_full[row,letter]=="BLM" ) {
					breakpoints_full[row,"gene"]="BLM"
				}

				else if (breakpoints_full[row,letter]=="RECQL5" | breakpoints_full[row,letter]=="recql5" | breakpoints_full[row,letter]=="RECQ5" | breakpoints_full[row,letter]=="recq5" ) {
					if (breakpoints_full[row,"gene"]=="BLM"){
						breakpoints_full[row,"gene"]="BLM/RECQL5"
					}
					else{
						breakpoints_full[row,"gene"]="RECQL5"
					}
				}
			}
		}

	}

	breakpoints_full <- select(breakpoints_full,-c(a,b,c,d,e,f))
	write.table(breakpoints_full,"Input/05.breakpoints.metrics.quality.gene.txt",quote = F,col.names = T,row.names = F,sep = "\t")

}
