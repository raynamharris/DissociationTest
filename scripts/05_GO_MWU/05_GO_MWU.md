GO\_MWU uses continuous measure of significance (such as fold-change or
-log(p-value) ) to identify GO categories that are significantly
enriches with either up- or down-regulated genes. The advantage - no
need to impose arbitrary significance cutoff.

If the measure is binary (0 or 1) the script will perform a typical “GO
enrichment” analysis based Fisher’s exact test: it will show GO
categories over-represented among the genes that have 1 as their
measure.

On the plot, different fonts are used to indicate significance and color
indicates enrichment with either up (red) or down (blue) regulated
genes. No colors are shown for binary measure analysis.

The tree on the plot is hierarchical clustering of GO categories based
on shared genes. Categories with no branch length between them are
subsets of each other.

The fraction next to GO category name indicates the fracton of “good”
genes in it; “good” genes being the ones exceeding the arbitrary
absValue cutoff (option in gomwuPlot). For Fisher’s based test, specify
absValue=0.5. This value does not affect statistics and is used for
plotting only.

Stretch the plot manually to match tree to text

Mikhail V. Matz, UT Austin, February 2015; <matz@utexas.edu>

################################################################ 

NOTES: This program drains memory and creates some very large
intermediate files, especially for the biological process catagory.

First, I run the stats from the command line to make sure its working.
Once I’ve generated the temp files, I comment out then stats portions
and recreate the plots by kniting the rmd file.

    library(ape)
    library(dplyr)

    ## Warning: package 'dplyr' was built under R version 3.5.1

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    source("gomwu.functions.R")

    # set output file for figures 
    knitr::opts_chunk$set(fig.path = '../../figures/05_GO_MMU/',
                         eval=F, include=F)

Dissociation vs Homogenization Molecular Function (MF)
------------------------------------------------------

    # input files
    input="01_dissociation_GOpvals.csv" 
    goAnnotations="goAnnotations.tab" 
    goDatabase="go.obo" 
    goDivision="MF" # either MF, or BP, or CC

    # Calculating stats
    #gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.1, smallest=5,clusterCutHeight=0.25)  

    # fewer catagories
    gomwuPlot(input,goAnnotations,goDivision,
        absValue=-log(0.05,10),  
        level1=0.001, 
        level2=0.001, 
        level3=0.0001, 
        txtsize=1,    
        treeHeight=0.5, 
      #colors=c("#d9d9d9","#525252","#d9d9d9","#525252")
        #colors=c("blue","green","blue","green") 
        colors=c("#969696","#000000","#bdbdbd","#252525") 
    )

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

![](../../figures/05_GO_MMU/MF-1.png)

    ## GO terms dispayed:  10 
    ## "Good genes" accounted for:  176 out of 1025 ( 17% )

    ##                                                                                                   pval
    ## 11/66 hydrogen ion transmembrane transporter activity                                     6.418384e-04
    ## 10/36 oxidoreductase activity, acting on NAD(P)H, quinone or similar compound as acceptor 3.221611e-04
    ## 50/596 oxidoreductase activity                                                            1.423643e-08
    ## 12/62 oxidoreductase activity, acting on NAD(P)H                                          4.228964e-04
    ## 32/433 ligase activity                                                                    3.221611e-04
    ## 19/245 ligase activity, forming carbon-nitrogen bonds                                     1.951264e-05
    ## 8/128 helicase activity                                                                   6.418384e-04
    ## 15/55 rRNA binding                                                                        7.701595e-05
    ## 42/88 structural constituent of ribosome                                                  1.000308e-15
    ## 74/325 structural molecule activity                                                       7.176782e-06
    ##                                                                                           direction
    ## 11/66 hydrogen ion transmembrane transporter activity                                             1
    ## 10/36 oxidoreductase activity, acting on NAD(P)H, quinone or similar compound as acceptor         1
    ## 50/596 oxidoreductase activity                                                                    1
    ## 12/62 oxidoreductase activity, acting on NAD(P)H                                                  1
    ## 32/433 ligase activity                                                                            0
    ## 19/245 ligase activity, forming carbon-nitrogen bonds                                             0
    ## 8/128 helicase activity                                                                           0
    ## 15/55 rRNA binding                                                                                1
    ## 42/88 structural constituent of ribosome                                                          1
    ## 74/325 structural molecule activity                                                               1
    ##                                                                                             color
    ## 11/66 hydrogen ion transmembrane transporter activity                                     #000000
    ## 10/36 oxidoreductase activity, acting on NAD(P)H, quinone or similar compound as acceptor #000000
    ## 50/596 oxidoreductase activity                                                            #000000
    ## 12/62 oxidoreductase activity, acting on NAD(P)H                                          #000000
    ## 32/433 ligase activity                                                                    #969696
    ## 19/245 ligase activity, forming carbon-nitrogen bonds                                     #969696
    ## 8/128 helicase activity                                                                   #969696
    ## 15/55 rRNA binding                                                                        #000000
    ## 42/88 structural constituent of ribosome                                                  #000000
    ## 74/325 structural molecule activity                                                       #000000
    ##                                                                                           enriched
    ## 11/66 hydrogen ion transmembrane transporter activity                                           NA
    ## 10/36 oxidoreductase activity, acting on NAD(P)H, quinone or similar compound as acceptor       NA
    ## 50/596 oxidoreductase activity                                                                  NA
    ## 12/62 oxidoreductase activity, acting on NAD(P)H                                                NA
    ## 32/433 ligase activity                                                                          NA
    ## 19/245 ligase activity, forming carbon-nitrogen bonds                                           NA
    ## 8/128 helicase activity                                                                         NA
    ## 15/55 rRNA binding                                                                              NA
    ## 42/88 structural constituent of ribosome                                                        NA
    ## 74/325 structural molecule activity                                                             NA

Dissociation vs Homogenization Molecular Function (MF)
------------------------------------------------------

    # input files
    input="01_dissociation_GOpvals.csv" 
    goAnnotations="goAnnotations.tab" 
    goDatabase="go.obo" 
    goDivision="CC" # either MF, or BP, or CC

    # Calculating stats
    #gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.1, smallest=5,clusterCutHeight=0.25)  

    # Data viz

    gomwuPlot(input,goAnnotations,goDivision,
        absValue=-log(0.05,10),  
        level1=0.00000001, 
        level2=0.00000001, 
        level3=0.000000001, 
        txtsize=1,    
        treeHeight=0.5, 
      #colors=c("#d9d9d9","#525252","#d9d9d9","#525252")
        #colors=c("blue","green","blue","green") 
        colors=c("#969696","#000000","#bdbdbd","#252525") 
    )

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

![](../../figures/05_GO_MMU/CC-1.png)

    ## GO terms dispayed:  14 
    ## "Good genes" accounted for:  475 out of 1084 ( 44% )

    ##                                                           pval direction
    ## 16/238 nuclear body                               6.306685e-09         0
    ## 38/529 nucleoplasm part                           6.346047e-14         0
    ## 60/208 ribosome                                   1.039230e-15         1
    ## 35/106 large ribosomal subunit                    5.200361e-12         1
    ## 33/53 cytosolic large ribosomal subunit           1.001255e-15         1
    ## 63/194 cytosolic part                             3.569336e-13         1
    ## 22/68 small ribosomal subunit                     6.201553e-11         1
    ## 20/41 cytosolic small ribosomal subunit           7.892550e-14         1
    ## 174/1126 extracellular membrane-bounded organelle 1.962229e-11         1
    ## 115/628 anchoring junction                        3.605100e-09         1
    ## 87/356 cell-substrate junction                    1.071864e-11         1
    ## 83/522 synapse part                               2.416116e-09         0
    ## 101/985 extracellular region                      4.122054e-15         1
    ## 122/957 extracellular space                       3.355731e-14         1
    ##                                                     color enriched
    ## 16/238 nuclear body                               #969696       NA
    ## 38/529 nucleoplasm part                           #969696       NA
    ## 60/208 ribosome                                   #000000       NA
    ## 35/106 large ribosomal subunit                    #000000       NA
    ## 33/53 cytosolic large ribosomal subunit           #000000       NA
    ## 63/194 cytosolic part                             #000000       NA
    ## 22/68 small ribosomal subunit                     #000000       NA
    ## 20/41 cytosolic small ribosomal subunit           #000000       NA
    ## 174/1126 extracellular membrane-bounded organelle #000000       NA
    ## 115/628 anchoring junction                        #000000       NA
    ## 87/356 cell-substrate junction                    #000000       NA
    ## 83/522 synapse part                               #969696       NA
    ## 101/985 extracellular region                      #000000       NA
    ## 122/957 extracellular space                       #000000       NA

table of significance, top 10 ish Go terms
------------------------------------------
