    library("dplyr") ## for filtering and selecting rows
    library("plyr")  ## for renmaing factors
    library("reshape2") ##  for melting dataframe

Disclaimer
----------

If you are viewing this page then hopefully you want to access
some really large data files containing raw RNA transcript counts and
estimates of transcripts per million. These files and this analysis will take up considerable space and time. 

Kallisto Gather
---------------

To obtain the data analyzed in this markdown file, I ran the [kallisto program](../UNIXworkflow/04_04_kallisto.md) the Stampede Cluster at the Texas Advacned Computing Facility. It runs really fast! The data are exported as abunance files in a subdirectory for every sample. 

The kallisto output gives you read counts for sample in an abundance
file for every single sample. This portion of the code goes through and
finds each samples’ abundance.tsv file, extracts the data, and combines
it all into a dataframe. The `counts` file is unnormalized, but the
`tpm` is the data after being normalized by transcripts per million.
This script was developed with assistance from Anna Batthenhouse and
Dennis Whylie.

Rather than examine unique transcripts, my analyses will focus on
gene-level exprrssion. I use some string splitting to take the very long
transcript identifying and create a `geneids` file that has all the
database identifiers for each transcript. Then, I’ll save the dount
data.

    setwd("../data/GSE99765_Dissociation/")
    ## this will create lists of all the samples
    kallistoDirs = dir(".")
    kallistoDirs = kallistoDirs[!grepl("\\.(R|py|pl|sh|xlsx?|txt|tsv|csv|org|md|obo|png|jpg|pdf)$",
    kallistoDirs, ignore.case=TRUE)]

    kallistoFiles = paste0(kallistoDirs, "/abundance.tsv")
    names(kallistoFiles) = kallistoDirs
    if(file.exists(kallistoFiles))
      kallistoData = lapply(
      kallistoFiles,
      read.table,
      sep = "\t",
      row.names = 1,
      header = TRUE
    )

    ## this for loop uses the reduce function to make two data frame with counts or tpm from all the samples. note, only counts are used

    ids = Reduce(f=union, x=lapply(kallistoData, rownames))
    if (all(sapply(kallistoData, function(x) {all(rownames(x)==ids)}))) {
      count = data.frame(
        id = ids,
        sapply(kallistoData, function(x) {x$est_counts}),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
      tpm = data.frame(
        id = ids,
        sapply(kallistoData, function(x) {x$tpm}),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    }

    ## make a dataframe with the parts of the gene id as columns
    geneids <- count[c(1)] 
    geneids$ENSMUST <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 1)
    geneids$ENSMUSG <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 2)
    geneids$OTTMUSG <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 3)
    geneids$OTTMUST <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 4)
    geneids$transcript <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 5)
    geneids$gene <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 6)
    geneids$length <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 7)
    geneids$structure1 <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 8)
    geneids$structure2 <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 9)
    geneids$structure3 <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 10)
    geneids$transcript_lenght <- as.factor(paste(geneids$transcript, geneids$length, sep="_"))


    # tpm to tpmbygene - note: not used for subsequent analyses
    tpmbygene <-  full_join(geneids, tpm) # merge tpm and genids

    ## Joining, by = "id"

    tpmbygene <- tpmbygene[-c(1:6,8:12)]   ## keep gene name and tpm for samples)
    tpmbygene <- melt(tpmbygene, id=c("gene")) ## lenghten 
    tpmbygene <- dcast(tpmbygene, gene ~ variable, value.var= "value", fun.aggregate=sum) #then widen by sum
    row.names(tpmbygene) <- tpmbygene$gene ## make gene the row name
    tpmbygene[1] <- NULL ## make gene the row name
    tpmbygene <- round(tpmbygene) #round all value to nearest 1s place


    # counts to prkcz_counts - used just for pkmz
    prkcz_counts <- full_join(geneids, count) # merge tpm and genids

    ## Joining, by = "id"

    prkcz_counts <- prkcz_counts[-c(1,2,4,5,12)]   ## keep gene name and tpm for samples)
    prkcz_counts <- dplyr::filter(prkcz_counts, grepl('Prkcz', gene)) %>%
                             arrange(transcript)
    prkcz_counts

    ##                 ENSMUSG transcript  gene length structure1   structure2
    ## 1 ENSMUSG00000029053.16  Prkcz-001 Prkcz   2641 UTR5:1-139 CDS:140-1918
    ## 2 ENSMUSG00000029053.16  Prkcz-002 Prkcz   4083 UTR5:1-538 CDS:539-1768
    ## 3 ENSMUSG00000029053.16  Prkcz-003 Prkcz    694 UTR5:1-118  CDS:119-409
    ## 4 ENSMUSG00000029053.16  Prkcz-004 Prkcz    586 UTR5:1-447  CDS:448-586
    ## 5 ENSMUSG00000029053.16  Prkcz-006 Prkcz    837 UTR5:1-473  CDS:474-837
    ##       structure3 100-CA1-1  100-CA1-2 100-CA1-3 100-CA3-1 100-CA3-4
    ## 1 UTR3:1919-2641   62.9547  127.49200   52.5939  18.69530         0
    ## 2 UTR3:1769-4083  394.0450 1119.62000  353.4060 319.41900       421
    ## 3   UTR3:410-694    0.0000    0.00000    0.0000   0.00000         0
    ## 4           <NA>    0.0000    0.00000    0.0000   0.00000         0
    ## 5           <NA>    0.0000    7.88499    0.0000   3.88531         0
    ##   100-DG-2 100-DG-3 101-CA1-1   101-CA1-2 101-CA1-3 101-CA3-1 101-CA3-4
    ## 1        0  35.8385   33.9907 0.00000e+00   15.8552  90.73590   23.7489
    ## 2      238 959.8530 1235.3400 9.99991e+00   29.6528 381.43800  194.2510
    ## 3        0   0.0000    0.0000 0.00000e+00    0.0000   0.00000    0.0000
    ## 4        0   0.0000    0.0000 9.03355e-05    5.4920   0.00000    0.0000
    ## 5        0  13.3085   12.6655 0.00000e+00    0.0000   7.82569    0.0000
    ##   101-DG-3 101-DG-4
    ## 1        0        0
    ## 2        7       98
    ## 3        0        0
    ## 4        0        0
    ## 5        0        0

    prkcz_counts <- prkcz_counts[-c(1,3:7)]   ## keep gene name and tpm for samples)
    row.names(prkcz_counts) <- prkcz_counts$transcript ## make gene the row name
    prkcz_counts[1] <- NULL ## make gene the row name
    prkcz_counts <- round(prkcz_counts) #round all value to nearest 1s place

    # count to countbygene
    countbygene <- full_join(geneids, count) # merge count and genids

    ## Joining, by = "id"

    countbygene <- countbygene[-c(1:6,8:12)]   ## rkeep gene name and counts for samples)
    countbygene <- melt(countbygene, id=c("gene")) ## lenghten 
    countbygene  <- dcast(countbygene, gene ~ variable, value.var= "value", fun.aggregate=sum) #then widen by sum
    row.names(countbygene) <- countbygene$gene ## make gene the row name
    countbygene[1] <- NULL ## make gene the row name
    countbygene <- round(countbygene) #round all value to nearest 1s place

    write.csv(geneids, "../geneids.csv", row.names=F)
    write.csv(prkcz_counts, "../GSE99765_Dissociation_prkcz_counts.csv", row.names=T)
    write.csv(countbygene, "../GSE99765_DissociationCountData.csv", row.names=T)

Session Info
------------

    sessionInfo()

    ## R version 3.5.0 (2018-04-23)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS  10.14
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] bindrcpp_0.2.2 reshape2_1.4.3 plyr_1.8.4     dplyr_0.7.6   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.18     knitr_1.20       bindr_0.1.1      magrittr_1.5    
    ##  [5] tidyselect_0.2.4 R6_2.2.2         rlang_0.2.2      stringr_1.3.1   
    ##  [9] tools_3.5.0      htmltools_0.3.6  yaml_2.2.0       assertthat_0.2.0
    ## [13] rprojroot_1.3-2  digest_0.6.17    tibble_1.4.2     crayon_1.3.4    
    ## [17] purrr_0.2.5      glue_1.3.0       evaluate_0.11    rmarkdown_1.10  
    ## [21] stringi_1.2.4    compiler_3.5.0   pillar_1.3.0     backports_1.1.2 
    ## [25] pkgconfig_2.0.2
