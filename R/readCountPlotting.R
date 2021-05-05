

readPlotting <- function(rdata="Output/bpr/data/",plot.dir = "Output/bpr/plots/",cluster.metrics="merge.metrics.background.txt",features=c("coverage","background","spikiness","evenness"),numOfLibs=20){
    
    cluster.metrics = read.table(cluster.metrics,header=T)
    
    for (feature in features){
        
        plots <- list()
        
        
        if (feature=="coverage"){
            df = cluster.metrics[order(cluster.metrics[,feature] ,decreasing = T),][1:numOfLibs,]
        } else {df = cluster.metrics[order(cluster.metrics[,feature] ,decreasing = F),][1:numOfLibs,]}
        
        file=df$file[1]
        
        for (file in df$file){
            f =file
            file=paste0(rdata,file,".RData")

            data <- get(load(file))#[[1]]
            
            filename <- data$ID
            ptm <- startTimedMessage("Plotting ", filename, " ...")
            
            bamfile <- data$ID
            reads <- data$fragments
            chroms2plot <- GenomeInfoDb::seqlevels(reads)
            breaks <- data$breaks
            counts <- data$counts
            lib.metrics <- data$lib.metrics
            lib.metrics["spikiness"]= cluster.metrics[cluster.metrics$file==f,]$spikiness
            lib.metrics["evenness"]= cluster.metrics[cluster.metrics$file==f,]$evenness.med
            lib.metrics <- round(lib.metrics, digits = 5)
            lib.metrics <- paste(names(lib.metrics), lib.metrics, sep = '=')
            lib.metrics <- paste(lib.metrics, collapse = "  |  ")
            
            #Skip chromosomes shorter then 5-times of the bin size 200kb
            if (any(seqlengths(reads) < 200000*5)) {
                message(" Skipping short chromosomes/contigs!")
                keep.chroms <- names(seqlengths(reads)[seqlengths(reads) >= 200000*5])
                reads <- GenomeInfoDb::keepSeqlevels(reads, keep.chroms, pruning.mode = 'coarse')
            }
            
            binned.data <- unlist(GenomicRanges::tileGenome(seqlengths(reads), tilewidth = 200000))
            
            #counts overlaps between bins and our reads
            Watsonreads <- GenomicRanges::countOverlaps(binned.data, reads[strand(reads)=='-'])
            Crickreads <- GenomicRanges::countOverlaps(binned.data, reads[strand(reads)=='+'])
            bothreads <- Watsonreads + Crickreads
            
            mcols(binned.data)$bothreads <- bothreads
            mcols(binned.data)$Watsonreads <- Watsonreads
            mcols(binned.data)$Crickreads <- Crickreads
            
            #transform bin coordinates of each chromosome into genomewide coordinates (cumulative sum of bin coordintes)
            cum.seqlengths <- cumsum(as.numeric(seqlengths(binned.data)))
            cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
            #get positions of ends of each chromosome to plot lones between the chromosomes
            if (length(cum.seqlengths) > 1) {
                chr.lines <- data.frame( y=cum.seqlengths[-length(cum.seqlengths)] )
            } else {
                chr.lines <- data.frame( y=0 )
            }
            #get positions of each chromosomes names
            chr.label.pos <- round( cum.seqlengths.0 + (0.5 * seqlengths(binned.data) ) )
            names(chr.label.pos) <- gsub("chr", "", names(chr.label.pos)) #line to add to exclude chr
            
            #transform chromosome based coordinates into genomewide coordinates
            trans.reads <- transCoord(binned.data)
            trans.breaks <- transCoord(breaks)
            trans.counts <- transCoord(counts)
            
            dfplot.reads <- as.data.frame(trans.reads)
            dfplot.breaks <- as.data.frame(trans.breaks)
            dfplot.counts <- as.data.frame(trans.counts)
            
            my_theme <- theme(
                legend.position="none",
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank()
            )
            
            
            ### PLOT READS
            
            #get midpoint values for each genomic bin
            dfplot.reads$midpoint <- dfplot.reads$start.genome + ( (dfplot.reads$end.genome - dfplot.reads$start.genome) %/% 2 )
            
            #filter bins with extremely high amount of reads
            Crickreads.outlier <- stats::quantile(dfplot.reads$Crickreads, 0.999)
            Watsonreads.outlier <- stats::quantile(dfplot.reads$Watsonreads, 0.999)
            #set outlier bins to the limit
            dfplot.reads$Crickreads[dfplot.reads$Crickreads >= Crickreads.outlier] <- Crickreads.outlier
            dfplot.reads$Watsonreads[dfplot.reads$Watsonreads >= Watsonreads.outlier] <- Watsonreads.outlier
            
            #construct ggplot
            dfplot.reads$mCrickreads <- -dfplot.reads$Crickreads
            ggplt1 <-
                ggplot(dfplot.reads) +
                geom_linerange(aes_string(ymin=0, ymax='mCrickreads', x='midpoint'), color="paleturquoise4", size=0.2) +
                geom_linerange(aes_string(ymin=0, ymax='Watsonreads', x='midpoint'), color="sandybrown", size=0.2 ) +
                geom_linerange(data=chr.lines, aes_string(ymin=-Inf, ymax=Inf, x='y'), col='black') + xlab(NULL) +
                ylab("Read counts") +
                xlab("Chromosomes") +
                scale_x_continuous(breaks=chr.label.pos, labels=names(chr.label.pos), expand = c(0,0)) +
                #theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
                my_theme +
                ggtitle(bamfile, subtitle = lib.metrics)
            
            #if (length(breaks)) {
            #	p <- suppressWarnings( cowplot::plot_grid(ggplt1, ggplt2, ggplt3, ncol=1, align="v", rel_heights = c(3,3,2)) )
            #} else {
            #	p <- suppressWarnings( cowplot::plot_grid(ggplt1, ncol=1, align="v", rel_heights = 3) )
            #}
            p=ggplt1
            plots[[length(plots)+1]] <- p
            stopTimedMessage(ptm)
        }
        
        printFile=paste0("Output/bpr/plots/",feature,".plot.pdf")
        message("Printing to PDF ",printFile)
            
        grDevices::pdf(printFile, width=max(10, length(chroms2plot)), height=2)
        bquiet = lapply(plots, print)
        d <- grDevices::dev.off()
        
        
    }
    
}






