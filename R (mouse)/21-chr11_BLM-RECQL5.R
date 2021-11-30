

chr11 = filter(cnvPerCell,seqnames=="chr11")


breakpointR::plotBreakpointsPerChr(files2plot = paste0("SCEs/Data-BreakpointR/data/",levels(droplevels(chr11$file)),".RData"),plotspath = "CNVs/",chromosomes = "chr11")
breakpointR::plotBreakpointsPerChr(files2plot = paste0("SCEs/Data-BreakpointR/data/",levels(droplevels(chr11$file)),".RData"),plotspath = "CNVs/",chromosomes = "chr12")


cnvpath=paste0("CNVs/chr11Duplication/cat")
if (!file.exists(cnvpath) ) { dir.create(cnvpath)}

export(chr11,"chr11_breakpoints.bed",format = "bed")
  
file.copy( paste0("SCEs/Data-BreakpointR/browserfiles/",levels(droplevels(chr11$file)),"_reads.bed.gz"),cnvpath)
  
plot(file,type="karyogram",plot.breakpoints=F,both.strands=F )
