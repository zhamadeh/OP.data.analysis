#######################################################
####     Workflow Part 6: Complete Dataset          ###
#######################################################

qualityFilterLibraries <- function(){

	metrics <- read.table("INPUT/2021/04.library.quality_2021.txt",header  =T) 
	metrics$Library <- as.factor(metrics$Library)
	#metrics%>% group_by(quality)%>% summarize(n())
	
	breakpoints<- read.table("INPUT/2021/03.library.breakpoints.txt",header=T)
  breakpoints$library<-as.factor(breakpoints$Library)
  
  breakpoints_full = merge(breakpoints,metrics,by="Library")
  
  #write.table(breakpoints_full,"INPUT/05.breakpoints.complete.txt",quote = F,col.names = T,row.names = F,sep = "\t")
  
  #breakpoints_full$library <- breakpoints_full$Library
	suppressWarnings(breakpoints_full <- breakpoints_full %>% separate(library, c("a","b","c","d","e","f"), "[_-]+"))
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
			  else if (breakpoints_full[row,letter]=="RECQL1"| breakpoints_full[row,letter]=="recql1"){
			    breakpoints_full[row,"gene"]="RECQL1"
			  }
			  else if (breakpoints_full[row,letter]=="RTEL"| breakpoints_full[row,letter]=="rtel"){
			    breakpoints_full[row,"gene"]="RTEL1"
			  }
			}
		}

	}

	breakpoints_full <- select(breakpoints_full,-c(a,b,c,d,e,f))
	write.table(breakpoints_full,"INPUT/2021/05.breakpoints.metrics.quality.gene_2021.txt",quote = F,col.names = T,row.names = F,sep = "\t")

}
