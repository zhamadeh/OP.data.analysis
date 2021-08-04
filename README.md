# One-pot Strand-seq Data analysis

### Short description
This is the master repository for all of Zeid's data analysis used in his PhD thesis. It has 9 parts each of which are described below to summarize and visualize findings of reccurrent and non-recurrent structural variants in single RecQ helicase deficient cells. 

### Steps of Data analysis

##### 1) Collect Initial Library Metrics
This involves collecting sequencing quality metrics such as spikiness, coverage and complexity from raw BAM files and storing them in `Input/01.library.metrics.txt`

##### 2) Generate breakpoint coordinates
This involves running breakpointR and generating RData files stored in `Output/bpr/data`

##### 3) Collect breakpoints from RData files
This involves collecting all breakpoints from above RData files into one file stored in `Input/03.library.breakpoints.txt`

##### 4) Colect metrics with quality scoring
This involves adding new BPR-specfic metrics to BAM-file metrics and doing automated threshold-based quality scoring and storing results in `Input/02.library.quality.txt`

##### 5) Plotting pairs of quality
This is just to plot Pairs plot for quality scoring

##### 6) Merge breakpoints with metrics, quality and genotype metadata
This merges breakpoints with all metadata files and assigns genotype too based off file name

##### 7) Blacklisting repetitive regions
Removing blaclisted regions stored in `Input/00.centromeres2.txt` and re-storing RData files in `Output/bpr/bl_data/`

##### 8) Re-plotting breakpointR plots
To convey efficiency of blacklisting

##### 9) Plotting summary
A few plots for summarizing SCE analsis



