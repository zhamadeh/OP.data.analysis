#######################################################
####     Workflow Part 3: Collect breakpoints       ###
#######################################################

## Collect breakpoints from RData files
collectBreaksAllFiles <- function(datapath="../StraVa/DATA/rdata/"){
  files <- list.files(datapath, pattern=".RData$", full.names=TRUE)
  
  
  breaks.all.files <- data.frame()
  n=1
  for (file in files) {
    message("Reading ... " , basename(file), " ... ",round(  (n/length(files))*100  ,  digits = 1  ) , "%"  )
    n=n+1
    data <- get(load(file))[c('breaks', 'confint','ID')]
    data$breaks$ID <- data$ID
    breakpoints <- as.data.frame(data$breaks)

    if (length(breakpoints)) {
      breaks.all.files <- rbind(breakpoints,breaks.all.files)
    }  
  }
  return(breaks.all.files)
}

