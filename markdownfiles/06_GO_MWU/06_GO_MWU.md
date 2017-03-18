GO\_MWU uses continuous measure of significance (such as fold-change or
-log(p-value) ) to identify GO categories that are significantly
enriches with either up- or down-regulated genes. The advantage - no
need to impose arbitrary significance cutoff.

If the measure is binary (0 or 1) the script will perform a typical "GO
enrichment" analysis based Fisher's exact test: it will show GO
categories over-represented among the genes that have 1 as their
measure.

On the plot, different fonts are used to indicate significance and color
indicates enrichment with either up (red) or down (blue) regulated
genes. No colors are shown for binary measure analysis.

The tree on the plot is hierarchical clustering of GO categories based
on shared genes. Categories with no branch length between them are
subsets of each other.

The fraction next to GO category name indicates the fracton of "good"
genes in it; "good" genes being the ones exceeding the arbitrary
absValue cutoff (option in gomwuPlot). For Fisher's based test, specify
absValue=0.5. This value does not affect statistics and is used for
plotting only.

Stretch the plot manually to match tree to text

Mikhail V. Matz, UT Austin, February 2015; <matz@utexas.edu>

################################################################ 

NOTES: This program drains memory and creates some very large
intermediate files, especially for the biological process catagory.

First, I run the stats from the command line to make sure its working.
Once I've generated the temp files, I comment out then stats portions
and recreate the plots by kniting the rmd file.

    library(ape)

    ## Warning: package 'ape' was built under R version 3.3.2

    source("gomwu.functions.R")

    # set output file for figures 
    knitr::opts_chunk$set(fig.path = '../../figures/06_GO_MMU/')

From Experitment 1: Dissociation Test Molecular Function (MF)
-------------------------------------------------------------

    # input files
    input="01_dissociation_GOpvals.csv" 
    goAnnotations="goAnnotations.tab" 
    goDatabase="go.obo" 
    goDivision="MF" # either MF, or BP, or CC

    # Calculating stats
    #gomwuStats(input, goDatabase, goAnnotations, #goDivision,
    #   perlPath="perl", 
    #   largest=0.1,  
    #   smallest=5,   
    #   clusterCutHeight=0.25)  

    # Data viz
    gomwuPlot(input,goAnnotations,goDivision,
        absValue=-log(0.05,10),  
        level1=0.001, 
        level2=0.0005, 
        level3=0.0001, 
        txtsize=1.4,    
        treeHeight=0.5, 
      colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") 
    )

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

![](../../figures/06_GO_MMU/01_dissociationMF-1.png)

    ## GO terms dispayed:  15 
    ## "Good genes" accounted for:  228 out of 706 ( 32% )

From Experitment 1: Dissociation Test Cellular Component (CC)
-------------------------------------------------------------

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

![](../../figures/06_GO_MMU/01_dissociationCC-1.png)

    ## GO terms dispayed:  23 
    ## "Good genes" accounted for:  429 out of 751 ( 57% )

From Experitment 1: Biological Process (BP)
-------------------------------------------

From Experitment 2: Stress Molecular Function (MF)
--------------------------------------------------

From Experitment 2: Stress Cellular Component (CC)
--------------------------------------------------

From Experitment 3: Cognition Cellular Component (CC)
-----------------------------------------------------

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

![](../../figures/06_GO_MMU/03_behaviorCC-1.png)

    ## GO terms dispayed:  20 
    ## "Good genes" accounted for:  555 out of 1298 ( 43% )

From Experitment 3: Cognition Molecular Function (MF)
-----------------------------------------------------

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

![](../../figures/06_GO_MMU/03_behaviorMF-1.png)

    ## GO terms dispayed:  29 
    ## "Good genes" accounted for:  541 out of 1236 ( 44% )

Now for Presence/Absence GO analysis
------------------------------------

-   First, I exported the list of overlapping genes from the previous
    05\_combo.Rmd file
-   Now, I'll use some joining majic to create a file with all the genes
    in the study and a 0 or 1 based on whether or not they were in the
    list of genes

<!-- -->

    library(dplyr)

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    #list of all genes. I only care about column 1
    allgenes <- read.csv("01_dissociation_GOpvals.csv", header = T)
    allgenes <- select(allgenes, gene)

    # first, identify genes differentially expressed by region in all experiments
    intersection <- read.csv("05_combo_intersection.csv", header=T)
    intersection$PresAbs <- 1
    intersection <- full_join(allgenes, intersection)

    ## Joining, by = "gene"

    ## Warning in full_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
    ## factors with different levels, coercing to character vector

    intersection <- intersection %>%
       mutate(PresAbs = replace(PresAbs,is.na(PresAbs),0))
    str(intersection)

    ## 'data.frame':    16272 obs. of  2 variables:
    ##  $ gene   : chr  "0610007P14Rik" "0610009B22Rik" "0610009L18Rik" "0610009O20Rik" ...
    ##  $ PresAbs: num  0 0 0 0 0 0 0 0 0 0 ...

    write.csv(intersection, "./05_combo_intersection_allgenes.csv", row.names = F)

The intersection: Cellular component
------------------------------------

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

![](../../figures/06_GO_MMU/05_intersectionCC-1.png)

    ## GO terms dispayed:  6 
    ## "Good genes" accounted for:  42 out of 111 ( 38% )

The intersection: Molecular Function
------------------------------------

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

![](../../figures/06_GO_MMU/05_intersectionMF-1.png)

    ## GO terms dispayed:  2 
    ## "Good genes" accounted for:  17 out of 107 ( 16% )
