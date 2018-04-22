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
    library(dplyr)

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

    # many catagories
    gomwuPlot(input,goAnnotations,goDivision,
        absValue=-log(0.05,10),  
        level1=0.05, 
        level2=0.05, 
        level3=0.001, 
        txtsize=1.4,    
        treeHeight=0.5, 
      #colors=c("#d9d9d9","#525252","#d9d9d9","#525252")
        #colors=c("blue","green","blue","green") 
        colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") 
    )

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

![](../../figures/05_GO_MMU/MF-1.png)

    ## GO terms dispayed:  58 
    ## "Good genes" accounted for:  533 out of 1025 ( 52% )

    ##                                                                                                     pval
    ## 2/20 threonine-type peptidase activity                                                      2.661592e-02
    ## 2/8 GABA receptor activity                                                                  2.897855e-02
    ## 2/86 ATP-dependent helicase activity                                                        4.877209e-03
    ## 2/58 RNA-dependent ATPase activity                                                          1.287400e-02
    ## 8/128 helicase activity                                                                     6.418384e-04
    ## 9/187 transferase activity, transferring one-carbon groups                                  1.000226e-02
    ## 0/23 tRNA methyltransferase activity                                                        4.367171e-02
    ## 6/122 S-adenosylmethionine-dependent methyltransferase activity                             8.876988e-03
    ## 3/15 histone methyltransferase activity (H3-K4 specific)                                    4.315074e-03
    ## 6/76 N-methyltransferase activity                                                           2.320277e-02
    ## 4/55 histone methyltransferase activity                                                     1.202857e-02
    ## 2/7 histone acetyltransferase activity (H4-K16 specific)                                    4.367171e-02
    ## 11/51 histone acetyltransferase activity                                                    1.533767e-03
    ## 16/93 acetyltransferase activity                                                            2.320277e-02
    ## 43/392 protein binding transcription factor activity                                        4.149728e-02
    ## 16/135 transcription corepressor activity                                                   4.182971e-02
    ## 43/453 chromatin binding                                                                    6.136555e-03
    ## 3/44 methylated histone residue binding                                                     8.954459e-03
    ## 18/153 histone binding                                                                      4.367171e-02
    ## 32/433 ligase activity                                                                      3.221611e-04
    ## 19/245 ligase activity, forming carbon-nitrogen bonds                                       1.951264e-05
    ## 4/52 small conjugating protein-specific protease activity                                   4.742286e-02
    ## 7/66 ubiquitinyl hydrolase activity                                                         4.877209e-03
    ## 1/20 ankyrin binding                                                                        3.904201e-02
    ## 20/144 voltage-gated ion channel activity                                                   1.202857e-02
    ## 3/5 voltage-gated calcium channel activity involved in cardiac muscle cell action potential 3.935120e-02
    ## 19/104 calcium ion transmembrane transporter activity                                       2.320277e-02
    ## 2/7 calcium:cation antiporter activity                                                      3.904201e-02
    ## 5/51 electron carrier activity                                                              3.403735e-02
    ## 9/150 coenzyme binding                                                                      4.204888e-02
    ## 12/62 oxidoreductase activity, acting on NAD(P)H                                            4.228964e-04
    ## 10/36 oxidoreductase activity, acting on NAD(P)H, quinone or similar compound as acceptor   3.221611e-04
    ## 8/51 antioxidant activity                                                                   5.390440e-03
    ## 50/596 oxidoreductase activity                                                              1.423643e-08
    ## 7/32 oxidoreductase activity, acting on peroxide as acceptor                                4.927224e-03
    ## 5/17 glutathione peroxidase activity                                                        2.089963e-02
    ## 42/302 GTPase binding                                                                       3.824039e-02
    ## 11/66 hydrogen ion transmembrane transporter activity                                       6.418384e-04
    ## 3/25 ATPase activity, coupled to transmembrane movement of ions, rotational mechanism       2.129751e-02
    ## 3/5 IgG binding                                                                             3.004682e-02
    ## 5/22 intramolecular oxidoreductase activity                                                 1.886267e-02
    ## 2/7 GPI-linked ephrin receptor activity                                                     4.204888e-02
    ## 73/669 kinase activity                                                                      3.935120e-02
    ## 42/403 protein serine/threonine kinase activity                                             2.787177e-02
    ## 17/72 integrin binding                                                                      1.897340e-03
    ## 40/230 cell adhesion molecule binding                                                       2.661592e-02
    ## 4/13 peptide antigen binding                                                                1.585990e-02
    ## 5/10 TAP binding                                                                            1.657911e-02
    ## 28/195 cytokine receptor binding                                                            4.367171e-02
    ## 146/1091 receptor binding                                                                   1.533767e-03
    ## 16/108 cytokine activity                                                                    6.136555e-03
    ## 20/149 glycosaminoglycan binding                                                            1.435711e-02
    ## 23/117 growth factor binding                                                                8.797936e-03
    ## 13/75 cytokine binding                                                                      7.276519e-03
    ## 4/11 BH domain binding                                                                      5.040794e-03
    ## 15/55 rRNA binding                                                                          7.701595e-05
    ## 42/88 structural constituent of ribosome                                                    1.000308e-15
    ## 74/325 structural molecule activity                                                         7.176782e-06
    ##                                                                                             direction
    ## 2/20 threonine-type peptidase activity                                                              1
    ## 2/8 GABA receptor activity                                                                          0
    ## 2/86 ATP-dependent helicase activity                                                                0
    ## 2/58 RNA-dependent ATPase activity                                                                  0
    ## 8/128 helicase activity                                                                             0
    ## 9/187 transferase activity, transferring one-carbon groups                                          0
    ## 0/23 tRNA methyltransferase activity                                                                0
    ## 6/122 S-adenosylmethionine-dependent methyltransferase activity                                     0
    ## 3/15 histone methyltransferase activity (H3-K4 specific)                                            0
    ## 6/76 N-methyltransferase activity                                                                   0
    ## 4/55 histone methyltransferase activity                                                             0
    ## 2/7 histone acetyltransferase activity (H4-K16 specific)                                            0
    ## 11/51 histone acetyltransferase activity                                                            0
    ## 16/93 acetyltransferase activity                                                                    0
    ## 43/392 protein binding transcription factor activity                                                0
    ## 16/135 transcription corepressor activity                                                           0
    ## 43/453 chromatin binding                                                                            0
    ## 3/44 methylated histone residue binding                                                             0
    ## 18/153 histone binding                                                                              0
    ## 32/433 ligase activity                                                                              0
    ## 19/245 ligase activity, forming carbon-nitrogen bonds                                               0
    ## 4/52 small conjugating protein-specific protease activity                                           0
    ## 7/66 ubiquitinyl hydrolase activity                                                                 0
    ## 1/20 ankyrin binding                                                                                0
    ## 20/144 voltage-gated ion channel activity                                                           0
    ## 3/5 voltage-gated calcium channel activity involved in cardiac muscle cell action potential         0
    ## 19/104 calcium ion transmembrane transporter activity                                               0
    ## 2/7 calcium:cation antiporter activity                                                              0
    ## 5/51 electron carrier activity                                                                      1
    ## 9/150 coenzyme binding                                                                              1
    ## 12/62 oxidoreductase activity, acting on NAD(P)H                                                    1
    ## 10/36 oxidoreductase activity, acting on NAD(P)H, quinone or similar compound as acceptor           1
    ## 8/51 antioxidant activity                                                                           1
    ## 50/596 oxidoreductase activity                                                                      1
    ## 7/32 oxidoreductase activity, acting on peroxide as acceptor                                        1
    ## 5/17 glutathione peroxidase activity                                                                1
    ## 42/302 GTPase binding                                                                               0
    ## 11/66 hydrogen ion transmembrane transporter activity                                               1
    ## 3/25 ATPase activity, coupled to transmembrane movement of ions, rotational mechanism               1
    ## 3/5 IgG binding                                                                                     1
    ## 5/22 intramolecular oxidoreductase activity                                                         1
    ## 2/7 GPI-linked ephrin receptor activity                                                             0
    ## 73/669 kinase activity                                                                              0
    ## 42/403 protein serine/threonine kinase activity                                                     0
    ## 17/72 integrin binding                                                                              1
    ## 40/230 cell adhesion molecule binding                                                               1
    ## 4/13 peptide antigen binding                                                                        1
    ## 5/10 TAP binding                                                                                    1
    ## 28/195 cytokine receptor binding                                                                    1
    ## 146/1091 receptor binding                                                                           1
    ## 16/108 cytokine activity                                                                            1
    ## 20/149 glycosaminoglycan binding                                                                    1
    ## 23/117 growth factor binding                                                                        1
    ## 13/75 cytokine binding                                                                              1
    ## 4/11 BH domain binding                                                                              1
    ## 15/55 rRNA binding                                                                                  1
    ## 42/88 structural constituent of ribosome                                                            1
    ## 74/325 structural molecule activity                                                                 1
    ##                                                                                                   color
    ## 2/20 threonine-type peptidase activity                                                       firebrick1
    ## 2/8 GABA receptor activity                                                                  dodgerblue2
    ## 2/86 ATP-dependent helicase activity                                                        dodgerblue2
    ## 2/58 RNA-dependent ATPase activity                                                          dodgerblue2
    ## 8/128 helicase activity                                                                     dodgerblue2
    ## 9/187 transferase activity, transferring one-carbon groups                                  dodgerblue2
    ## 0/23 tRNA methyltransferase activity                                                        dodgerblue2
    ## 6/122 S-adenosylmethionine-dependent methyltransferase activity                             dodgerblue2
    ## 3/15 histone methyltransferase activity (H3-K4 specific)                                    dodgerblue2
    ## 6/76 N-methyltransferase activity                                                           dodgerblue2
    ## 4/55 histone methyltransferase activity                                                     dodgerblue2
    ## 2/7 histone acetyltransferase activity (H4-K16 specific)                                    dodgerblue2
    ## 11/51 histone acetyltransferase activity                                                    dodgerblue2
    ## 16/93 acetyltransferase activity                                                            dodgerblue2
    ## 43/392 protein binding transcription factor activity                                        dodgerblue2
    ## 16/135 transcription corepressor activity                                                   dodgerblue2
    ## 43/453 chromatin binding                                                                    dodgerblue2
    ## 3/44 methylated histone residue binding                                                     dodgerblue2
    ## 18/153 histone binding                                                                      dodgerblue2
    ## 32/433 ligase activity                                                                      dodgerblue2
    ## 19/245 ligase activity, forming carbon-nitrogen bonds                                       dodgerblue2
    ## 4/52 small conjugating protein-specific protease activity                                   dodgerblue2
    ## 7/66 ubiquitinyl hydrolase activity                                                         dodgerblue2
    ## 1/20 ankyrin binding                                                                        dodgerblue2
    ## 20/144 voltage-gated ion channel activity                                                   dodgerblue2
    ## 3/5 voltage-gated calcium channel activity involved in cardiac muscle cell action potential dodgerblue2
    ## 19/104 calcium ion transmembrane transporter activity                                       dodgerblue2
    ## 2/7 calcium:cation antiporter activity                                                      dodgerblue2
    ## 5/51 electron carrier activity                                                               firebrick1
    ## 9/150 coenzyme binding                                                                       firebrick1
    ## 12/62 oxidoreductase activity, acting on NAD(P)H                                             firebrick1
    ## 10/36 oxidoreductase activity, acting on NAD(P)H, quinone or similar compound as acceptor    firebrick1
    ## 8/51 antioxidant activity                                                                    firebrick1
    ## 50/596 oxidoreductase activity                                                               firebrick1
    ## 7/32 oxidoreductase activity, acting on peroxide as acceptor                                 firebrick1
    ## 5/17 glutathione peroxidase activity                                                         firebrick1
    ## 42/302 GTPase binding                                                                       dodgerblue2
    ## 11/66 hydrogen ion transmembrane transporter activity                                        firebrick1
    ## 3/25 ATPase activity, coupled to transmembrane movement of ions, rotational mechanism        firebrick1
    ## 3/5 IgG binding                                                                              firebrick1
    ## 5/22 intramolecular oxidoreductase activity                                                  firebrick1
    ## 2/7 GPI-linked ephrin receptor activity                                                     dodgerblue2
    ## 73/669 kinase activity                                                                      dodgerblue2
    ## 42/403 protein serine/threonine kinase activity                                             dodgerblue2
    ## 17/72 integrin binding                                                                       firebrick1
    ## 40/230 cell adhesion molecule binding                                                        firebrick1
    ## 4/13 peptide antigen binding                                                                 firebrick1
    ## 5/10 TAP binding                                                                             firebrick1
    ## 28/195 cytokine receptor binding                                                             firebrick1
    ## 146/1091 receptor binding                                                                    firebrick1
    ## 16/108 cytokine activity                                                                     firebrick1
    ## 20/149 glycosaminoglycan binding                                                             firebrick1
    ## 23/117 growth factor binding                                                                 firebrick1
    ## 13/75 cytokine binding                                                                       firebrick1
    ## 4/11 BH domain binding                                                                       firebrick1
    ## 15/55 rRNA binding                                                                           firebrick1
    ## 42/88 structural constituent of ribosome                                                     firebrick1
    ## 74/325 structural molecule activity                                                          firebrick1
    ##                                                                                             enriched
    ## 2/20 threonine-type peptidase activity                                                          DISS
    ## 2/8 GABA receptor activity                                                                      HOMO
    ## 2/86 ATP-dependent helicase activity                                                            HOMO
    ## 2/58 RNA-dependent ATPase activity                                                              HOMO
    ## 8/128 helicase activity                                                                         HOMO
    ## 9/187 transferase activity, transferring one-carbon groups                                      HOMO
    ## 0/23 tRNA methyltransferase activity                                                            HOMO
    ## 6/122 S-adenosylmethionine-dependent methyltransferase activity                                 HOMO
    ## 3/15 histone methyltransferase activity (H3-K4 specific)                                        HOMO
    ## 6/76 N-methyltransferase activity                                                               HOMO
    ## 4/55 histone methyltransferase activity                                                         HOMO
    ## 2/7 histone acetyltransferase activity (H4-K16 specific)                                        HOMO
    ## 11/51 histone acetyltransferase activity                                                        HOMO
    ## 16/93 acetyltransferase activity                                                                HOMO
    ## 43/392 protein binding transcription factor activity                                            HOMO
    ## 16/135 transcription corepressor activity                                                       HOMO
    ## 43/453 chromatin binding                                                                        HOMO
    ## 3/44 methylated histone residue binding                                                         HOMO
    ## 18/153 histone binding                                                                          HOMO
    ## 32/433 ligase activity                                                                          HOMO
    ## 19/245 ligase activity, forming carbon-nitrogen bonds                                           HOMO
    ## 4/52 small conjugating protein-specific protease activity                                       HOMO
    ## 7/66 ubiquitinyl hydrolase activity                                                             HOMO
    ## 1/20 ankyrin binding                                                                            HOMO
    ## 20/144 voltage-gated ion channel activity                                                       HOMO
    ## 3/5 voltage-gated calcium channel activity involved in cardiac muscle cell action potential     HOMO
    ## 19/104 calcium ion transmembrane transporter activity                                           HOMO
    ## 2/7 calcium:cation antiporter activity                                                          HOMO
    ## 5/51 electron carrier activity                                                                  DISS
    ## 9/150 coenzyme binding                                                                          DISS
    ## 12/62 oxidoreductase activity, acting on NAD(P)H                                                DISS
    ## 10/36 oxidoreductase activity, acting on NAD(P)H, quinone or similar compound as acceptor       DISS
    ## 8/51 antioxidant activity                                                                       DISS
    ## 50/596 oxidoreductase activity                                                                  DISS
    ## 7/32 oxidoreductase activity, acting on peroxide as acceptor                                    DISS
    ## 5/17 glutathione peroxidase activity                                                            DISS
    ## 42/302 GTPase binding                                                                           HOMO
    ## 11/66 hydrogen ion transmembrane transporter activity                                           DISS
    ## 3/25 ATPase activity, coupled to transmembrane movement of ions, rotational mechanism           DISS
    ## 3/5 IgG binding                                                                                 DISS
    ## 5/22 intramolecular oxidoreductase activity                                                     DISS
    ## 2/7 GPI-linked ephrin receptor activity                                                         HOMO
    ## 73/669 kinase activity                                                                          HOMO
    ## 42/403 protein serine/threonine kinase activity                                                 HOMO
    ## 17/72 integrin binding                                                                          DISS
    ## 40/230 cell adhesion molecule binding                                                           DISS
    ## 4/13 peptide antigen binding                                                                    DISS
    ## 5/10 TAP binding                                                                                DISS
    ## 28/195 cytokine receptor binding                                                                DISS
    ## 146/1091 receptor binding                                                                       DISS
    ## 16/108 cytokine activity                                                                        DISS
    ## 20/149 glycosaminoglycan binding                                                                DISS
    ## 23/117 growth factor binding                                                                    DISS
    ## 13/75 cytokine binding                                                                          DISS
    ## 4/11 BH domain binding                                                                          DISS
    ## 15/55 rRNA binding                                                                              DISS
    ## 42/88 structural constituent of ribosome                                                        DISS
    ## 74/325 structural molecule activity                                                             DISS

    # fewer catagories
    gomwuPlot(input,goAnnotations,goDivision,
        absValue=-log(0.05,10),  
        level1=0.001, 
        level2=0.001, 
        level3=0.0001, 
        txtsize=1.4,    
        treeHeight=0.5, 
      #colors=c("#d9d9d9","#525252","#d9d9d9","#525252")
        #colors=c("blue","green","blue","green") 
        colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") 
    )

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

![](../../figures/05_GO_MMU/MF-2.png)

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
    ##                                                                                                 color
    ## 11/66 hydrogen ion transmembrane transporter activity                                      firebrick1
    ## 10/36 oxidoreductase activity, acting on NAD(P)H, quinone or similar compound as acceptor  firebrick1
    ## 50/596 oxidoreductase activity                                                             firebrick1
    ## 12/62 oxidoreductase activity, acting on NAD(P)H                                           firebrick1
    ## 32/433 ligase activity                                                                    dodgerblue2
    ## 19/245 ligase activity, forming carbon-nitrogen bonds                                     dodgerblue2
    ## 8/128 helicase activity                                                                   dodgerblue2
    ## 15/55 rRNA binding                                                                         firebrick1
    ## 42/88 structural constituent of ribosome                                                   firebrick1
    ## 74/325 structural molecule activity                                                        firebrick1
    ##                                                                                           enriched
    ## 11/66 hydrogen ion transmembrane transporter activity                                         DISS
    ## 10/36 oxidoreductase activity, acting on NAD(P)H, quinone or similar compound as acceptor     DISS
    ## 50/596 oxidoreductase activity                                                                DISS
    ## 12/62 oxidoreductase activity, acting on NAD(P)H                                              DISS
    ## 32/433 ligase activity                                                                        HOMO
    ## 19/245 ligase activity, forming carbon-nitrogen bonds                                         HOMO
    ## 8/128 helicase activity                                                                       HOMO
    ## 15/55 rRNA binding                                                                            DISS
    ## 42/88 structural constituent of ribosome                                                      DISS
    ## 74/325 structural molecule activity                                                           DISS

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
        level1=0.05, 
        level2=0.05, 
        level3=0.00000001, 
        txtsize=1.4,    
        treeHeight=0.5, 
      #colors=c("#d9d9d9","#525252","#d9d9d9","#525252")
        #colors=c("blue","green","blue","green") 
        colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") 
    )

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

![](../../figures/05_GO_MMU/CC-1.png)

    ## GO terms dispayed:  103 
    ## "Good genes" accounted for:  923 out of 1084 ( 85% )

    ##                                                                                        pval
    ## 2/20 proteasome core complex                                                   3.870397e-02
    ## 61/675 nucleolus                                                               3.081924e-02
    ## 60/208 ribosome                                                                1.039230e-15
    ## 35/106 large ribosomal subunit                                                 5.200361e-12
    ## 33/53 cytosolic large ribosomal subunit                                        1.001255e-15
    ## 22/68 small ribosomal subunit                                                  6.201553e-11
    ## 20/41 cytosolic small ribosomal subunit                                        7.892550e-14
    ## 11/153 ubiquitin ligase complex                                                6.526940e-03
    ## 4/15 protein serine/threonine phosphatase complex                              2.243001e-02
    ## 20/303 transferase complex                                                     3.850809e-05
    ## 5/63 protein acetyltransferase complex                                         4.616268e-02
    ## 55/645 catalytic complex                                                       1.250819e-03
    ## 9/74 methyltransferase complex                                                 7.691701e-05
    ## 1/10 Set1C/COMPASS complex                                                     2.243001e-02
    ## 2/59 protein kinase complex                                                    6.526940e-03
    ## 1/8 IkappaB kinase complex                                                     3.025213e-02
    ## 63/194 cytosolic part                                                          3.569336e-13
    ## 6/112 chromosome, telomeric region                                             3.426720e-02
    ## 36/360 nuclear chromosome part                                                 1.495433e-02
    ## 1/14 Sin3-type complex                                                         2.211665e-02
    ## 1/18 centromeric heterochromatin                                               2.211665e-02
    ## 31/262 chromatin                                                               6.526940e-03
    ## 3/71 heterochromatin                                                           2.795009e-02
    ## 18/331 chromosome                                                              1.912849e-05
    ## 4/76 chromosome, centromeric region                                            2.243001e-02
    ## 53/592 chromosomal part                                                        4.222348e-05
    ## 5/46 transcription elongation factor complex                                   2.211665e-02
    ## 4/71 PML body                                                                  1.250819e-03
    ## 38/529 nucleoplasm part                                                        6.346047e-14
    ## 2/14 mRNA cleavage and polyadenylation specificity factor complex              3.821024e-02
    ## 2/19 mRNA cleavage factor complex                                              4.019587e-02
    ## 6/99 nuclear speck                                                             3.013573e-02
    ## 1/34 Cajal body                                                                2.677524e-02
    ## 16/238 nuclear body                                                            6.306685e-09
    ## 122/957 extracellular space                                                    3.355731e-14
    ## 101/985 extracellular region                                                   4.122054e-15
    ## 23/247 proteinaceous extracellular matrix                                      9.952163e-04
    ## 51/384 extracellular matrix                                                    1.246720e-08
    ## 113/1046 Golgi apparatus                                                       6.883147e-04
    ## 0/9 Golgi transport complex                                                    1.733040e-02
    ## 39/272 side of membrane                                                        5.641465e-06
    ## 0/11 immunoglobulin complex                                                    3.106132e-02
    ## 4/42 cell periphery                                                            2.211665e-02
    ## 61/381 plasma membrane region                                                  3.426720e-02
    ## 4/17 MHC protein complex                                                       6.241488e-03
    ## 4/10 MHC class I protein complex                                               4.550434e-03
    ## 219/1505 membrane-bounded vesicle                                              5.698732e-08
    ## 174/1126 extracellular membrane-bounded organelle                              1.962229e-11
    ## 8/30 integrin complex                                                          2.050602e-02
    ## 89/507 cell surface                                                            1.392771e-08
    ## 5/12 phagocytic cup                                                            2.715347e-03
    ## 88/590 cell-cell junction                                                      2.795009e-02
    ## 115/628 anchoring junction                                                     3.605100e-09
    ## 87/356 cell-substrate junction                                                 1.071864e-11
    ## 2/44 spindle microtubule                                                       4.811212e-02
    ## 8/102 ciliary basal body                                                       7.799447e-05
    ## 1/46 spindle pole                                                              3.046713e-02
    ## 151/1308 cytoskeletal part                                                     1.250819e-03
    ## 0/14 gamma-tubulin complex                                                     2.939908e-02
    ## 40/516 microtubule organizing center                                           1.987097e-03
    ## 138/1238 cytoskeleton                                                          5.794676e-03
    ## 6/92 centriole                                                                 6.241488e-03
    ## 7/130 microtubule organizing center part                                       3.285735e-05
    ## 0/10 mRNA cap binding complex                                                  3.081924e-02
    ## 5/19 sarcoplasmic reticulum membrane                                           5.133588e-03
    ## 44/250 transporter complex                                                     2.618733e-03
    ## 25/129 cation channel complex                                                  2.795009e-02
    ## 23/161 site of polarized growth                                                4.992248e-02
    ## 2/8 ribbon synapse                                                             3.271411e-02
    ## 8/49 neuron projection membrane                                                3.000094e-02
    ## 123/864 cell projection part                                                   3.113335e-04
    ## 40/257 axon part                                                               4.690764e-02
    ## 6/12 axon initial segment                                                      2.831850e-02
    ## 0/16 inhibitory synapse                                                        4.726500e-02
    ## 3/15 asymmetric synapse                                                        2.795009e-02
    ## 69/506 synapse                                                                 1.362744e-03
    ## 39/237 synaptic membrane                                                       1.246720e-08
    ## 42/210 postsynaptic density                                                    5.704602e-07
    ## 83/522 synapse part                                                            2.416116e-09
    ## 25/130 neuron spine                                                            2.931824e-03
    ## 11/58 dendritic shaft                                                          5.189127e-04
    ## 68/413 dendrite                                                                2.751683e-05
    ## 189/1239 neuron part                                                           8.395302e-03
    ## 6/76 pore complex                                                              7.680869e-03
    ## 2/18 nuclear periphery                                                         4.811212e-02
    ## 51/347 vacuole                                                                 2.243001e-02
    ## 0/7 TORC1 complex                                                              1.194071e-02
    ## 28/175 myelin sheath                                                           5.826920e-07
    ## 10/87 oxidoreductase complex                                                   1.109797e-03
    ## 24/170 mitochondrial membrane part                                             1.086768e-05
    ## 6/42 respiratory chain complex I                                               3.967990e-06
    ## 67/725 mitochondrial part                                                      5.369859e-07
    ## 6/14 respiratory chain complex IV                                              1.226445e-02
    ## 8/26 cytochrome complex                                                        3.113335e-04
    ## 51/495 mitochondrial membrane                                                  3.722259e-07
    ## 131/1373 organelle membrane                                                    2.211665e-02
    ## 39/386 organelle inner membrane                                                1.537315e-08
    ## 2/11 respiratory chain complex III                                             4.992248e-02
    ## 2/19 proton-transporting two-sector ATPase complex, proton-transporting domain 1.450611e-02
    ## 2/10 proton-transporting ATP synthase complex, coupling factor F(o)            2.243001e-02
    ## 5/29 proton-transporting two-sector ATPase complex                             9.841256e-05
    ## 125/1424 mitochondrion                                                         3.907507e-04
    ## 4/16 mitochondrial proton-transporting ATP synthase complex                    9.889536e-05
    ##                                                                                direction
    ## 2/20 proteasome core complex                                                           1
    ## 61/675 nucleolus                                                                       0
    ## 60/208 ribosome                                                                        1
    ## 35/106 large ribosomal subunit                                                         1
    ## 33/53 cytosolic large ribosomal subunit                                                1
    ## 22/68 small ribosomal subunit                                                          1
    ## 20/41 cytosolic small ribosomal subunit                                                1
    ## 11/153 ubiquitin ligase complex                                                        0
    ## 4/15 protein serine/threonine phosphatase complex                                      0
    ## 20/303 transferase complex                                                             0
    ## 5/63 protein acetyltransferase complex                                                 0
    ## 55/645 catalytic complex                                                               0
    ## 9/74 methyltransferase complex                                                         0
    ## 1/10 Set1C/COMPASS complex                                                             0
    ## 2/59 protein kinase complex                                                            0
    ## 1/8 IkappaB kinase complex                                                             0
    ## 63/194 cytosolic part                                                                  1
    ## 6/112 chromosome, telomeric region                                                     0
    ## 36/360 nuclear chromosome part                                                         0
    ## 1/14 Sin3-type complex                                                                 0
    ## 1/18 centromeric heterochromatin                                                       0
    ## 31/262 chromatin                                                                       0
    ## 3/71 heterochromatin                                                                   0
    ## 18/331 chromosome                                                                      0
    ## 4/76 chromosome, centromeric region                                                    0
    ## 53/592 chromosomal part                                                                0
    ## 5/46 transcription elongation factor complex                                           0
    ## 4/71 PML body                                                                          0
    ## 38/529 nucleoplasm part                                                                0
    ## 2/14 mRNA cleavage and polyadenylation specificity factor complex                      0
    ## 2/19 mRNA cleavage factor complex                                                      0
    ## 6/99 nuclear speck                                                                     0
    ## 1/34 Cajal body                                                                        0
    ## 16/238 nuclear body                                                                    0
    ## 122/957 extracellular space                                                            1
    ## 101/985 extracellular region                                                           1
    ## 23/247 proteinaceous extracellular matrix                                              1
    ## 51/384 extracellular matrix                                                            1
    ## 113/1046 Golgi apparatus                                                               0
    ## 0/9 Golgi transport complex                                                            0
    ## 39/272 side of membrane                                                                1
    ## 0/11 immunoglobulin complex                                                            1
    ## 4/42 cell periphery                                                                    1
    ## 61/381 plasma membrane region                                                          1
    ## 4/17 MHC protein complex                                                               1
    ## 4/10 MHC class I protein complex                                                       1
    ## 219/1505 membrane-bounded vesicle                                                      1
    ## 174/1126 extracellular membrane-bounded organelle                                      1
    ## 8/30 integrin complex                                                                  1
    ## 89/507 cell surface                                                                    1
    ## 5/12 phagocytic cup                                                                    1
    ## 88/590 cell-cell junction                                                              1
    ## 115/628 anchoring junction                                                             1
    ## 87/356 cell-substrate junction                                                         1
    ## 2/44 spindle microtubule                                                               0
    ## 8/102 ciliary basal body                                                               0
    ## 1/46 spindle pole                                                                      0
    ## 151/1308 cytoskeletal part                                                             0
    ## 0/14 gamma-tubulin complex                                                             0
    ## 40/516 microtubule organizing center                                                   0
    ## 138/1238 cytoskeleton                                                                  0
    ## 6/92 centriole                                                                         0
    ## 7/130 microtubule organizing center part                                               0
    ## 0/10 mRNA cap binding complex                                                          0
    ## 5/19 sarcoplasmic reticulum membrane                                                   0
    ## 44/250 transporter complex                                                             0
    ## 25/129 cation channel complex                                                          0
    ## 23/161 site of polarized growth                                                        0
    ## 2/8 ribbon synapse                                                                     0
    ## 8/49 neuron projection membrane                                                        0
    ## 123/864 cell projection part                                                           0
    ## 40/257 axon part                                                                       0
    ## 6/12 axon initial segment                                                              0
    ## 0/16 inhibitory synapse                                                                0
    ## 3/15 asymmetric synapse                                                                0
    ## 69/506 synapse                                                                         0
    ## 39/237 synaptic membrane                                                               0
    ## 42/210 postsynaptic density                                                            0
    ## 83/522 synapse part                                                                    0
    ## 25/130 neuron spine                                                                    0
    ## 11/58 dendritic shaft                                                                  0
    ## 68/413 dendrite                                                                        0
    ## 189/1239 neuron part                                                                   0
    ## 6/76 pore complex                                                                      0
    ## 2/18 nuclear periphery                                                                 0
    ## 51/347 vacuole                                                                         1
    ## 0/7 TORC1 complex                                                                      0
    ## 28/175 myelin sheath                                                                   1
    ## 10/87 oxidoreductase complex                                                           1
    ## 24/170 mitochondrial membrane part                                                     1
    ## 6/42 respiratory chain complex I                                                       1
    ## 67/725 mitochondrial part                                                              1
    ## 6/14 respiratory chain complex IV                                                      1
    ## 8/26 cytochrome complex                                                                1
    ## 51/495 mitochondrial membrane                                                          1
    ## 131/1373 organelle membrane                                                            1
    ## 39/386 organelle inner membrane                                                        1
    ## 2/11 respiratory chain complex III                                                     1
    ## 2/19 proton-transporting two-sector ATPase complex, proton-transporting domain         1
    ## 2/10 proton-transporting ATP synthase complex, coupling factor F(o)                    1
    ## 5/29 proton-transporting two-sector ATPase complex                                     1
    ## 125/1424 mitochondrion                                                                 1
    ## 4/16 mitochondrial proton-transporting ATP synthase complex                            1
    ##                                                                                      color
    ## 2/20 proteasome core complex                                                    firebrick1
    ## 61/675 nucleolus                                                               dodgerblue2
    ## 60/208 ribosome                                                                 firebrick1
    ## 35/106 large ribosomal subunit                                                  firebrick1
    ## 33/53 cytosolic large ribosomal subunit                                         firebrick1
    ## 22/68 small ribosomal subunit                                                   firebrick1
    ## 20/41 cytosolic small ribosomal subunit                                         firebrick1
    ## 11/153 ubiquitin ligase complex                                                dodgerblue2
    ## 4/15 protein serine/threonine phosphatase complex                              dodgerblue2
    ## 20/303 transferase complex                                                     dodgerblue2
    ## 5/63 protein acetyltransferase complex                                         dodgerblue2
    ## 55/645 catalytic complex                                                       dodgerblue2
    ## 9/74 methyltransferase complex                                                 dodgerblue2
    ## 1/10 Set1C/COMPASS complex                                                     dodgerblue2
    ## 2/59 protein kinase complex                                                    dodgerblue2
    ## 1/8 IkappaB kinase complex                                                     dodgerblue2
    ## 63/194 cytosolic part                                                           firebrick1
    ## 6/112 chromosome, telomeric region                                             dodgerblue2
    ## 36/360 nuclear chromosome part                                                 dodgerblue2
    ## 1/14 Sin3-type complex                                                         dodgerblue2
    ## 1/18 centromeric heterochromatin                                               dodgerblue2
    ## 31/262 chromatin                                                               dodgerblue2
    ## 3/71 heterochromatin                                                           dodgerblue2
    ## 18/331 chromosome                                                              dodgerblue2
    ## 4/76 chromosome, centromeric region                                            dodgerblue2
    ## 53/592 chromosomal part                                                        dodgerblue2
    ## 5/46 transcription elongation factor complex                                   dodgerblue2
    ## 4/71 PML body                                                                  dodgerblue2
    ## 38/529 nucleoplasm part                                                        dodgerblue2
    ## 2/14 mRNA cleavage and polyadenylation specificity factor complex              dodgerblue2
    ## 2/19 mRNA cleavage factor complex                                              dodgerblue2
    ## 6/99 nuclear speck                                                             dodgerblue2
    ## 1/34 Cajal body                                                                dodgerblue2
    ## 16/238 nuclear body                                                            dodgerblue2
    ## 122/957 extracellular space                                                     firebrick1
    ## 101/985 extracellular region                                                    firebrick1
    ## 23/247 proteinaceous extracellular matrix                                       firebrick1
    ## 51/384 extracellular matrix                                                     firebrick1
    ## 113/1046 Golgi apparatus                                                       dodgerblue2
    ## 0/9 Golgi transport complex                                                    dodgerblue2
    ## 39/272 side of membrane                                                         firebrick1
    ## 0/11 immunoglobulin complex                                                     firebrick1
    ## 4/42 cell periphery                                                             firebrick1
    ## 61/381 plasma membrane region                                                   firebrick1
    ## 4/17 MHC protein complex                                                        firebrick1
    ## 4/10 MHC class I protein complex                                                firebrick1
    ## 219/1505 membrane-bounded vesicle                                               firebrick1
    ## 174/1126 extracellular membrane-bounded organelle                               firebrick1
    ## 8/30 integrin complex                                                           firebrick1
    ## 89/507 cell surface                                                             firebrick1
    ## 5/12 phagocytic cup                                                             firebrick1
    ## 88/590 cell-cell junction                                                       firebrick1
    ## 115/628 anchoring junction                                                      firebrick1
    ## 87/356 cell-substrate junction                                                  firebrick1
    ## 2/44 spindle microtubule                                                       dodgerblue2
    ## 8/102 ciliary basal body                                                       dodgerblue2
    ## 1/46 spindle pole                                                              dodgerblue2
    ## 151/1308 cytoskeletal part                                                     dodgerblue2
    ## 0/14 gamma-tubulin complex                                                     dodgerblue2
    ## 40/516 microtubule organizing center                                           dodgerblue2
    ## 138/1238 cytoskeleton                                                          dodgerblue2
    ## 6/92 centriole                                                                 dodgerblue2
    ## 7/130 microtubule organizing center part                                       dodgerblue2
    ## 0/10 mRNA cap binding complex                                                  dodgerblue2
    ## 5/19 sarcoplasmic reticulum membrane                                           dodgerblue2
    ## 44/250 transporter complex                                                     dodgerblue2
    ## 25/129 cation channel complex                                                  dodgerblue2
    ## 23/161 site of polarized growth                                                dodgerblue2
    ## 2/8 ribbon synapse                                                             dodgerblue2
    ## 8/49 neuron projection membrane                                                dodgerblue2
    ## 123/864 cell projection part                                                   dodgerblue2
    ## 40/257 axon part                                                               dodgerblue2
    ## 6/12 axon initial segment                                                      dodgerblue2
    ## 0/16 inhibitory synapse                                                        dodgerblue2
    ## 3/15 asymmetric synapse                                                        dodgerblue2
    ## 69/506 synapse                                                                 dodgerblue2
    ## 39/237 synaptic membrane                                                       dodgerblue2
    ## 42/210 postsynaptic density                                                    dodgerblue2
    ## 83/522 synapse part                                                            dodgerblue2
    ## 25/130 neuron spine                                                            dodgerblue2
    ## 11/58 dendritic shaft                                                          dodgerblue2
    ## 68/413 dendrite                                                                dodgerblue2
    ## 189/1239 neuron part                                                           dodgerblue2
    ## 6/76 pore complex                                                              dodgerblue2
    ## 2/18 nuclear periphery                                                         dodgerblue2
    ## 51/347 vacuole                                                                  firebrick1
    ## 0/7 TORC1 complex                                                              dodgerblue2
    ## 28/175 myelin sheath                                                            firebrick1
    ## 10/87 oxidoreductase complex                                                    firebrick1
    ## 24/170 mitochondrial membrane part                                              firebrick1
    ## 6/42 respiratory chain complex I                                                firebrick1
    ## 67/725 mitochondrial part                                                       firebrick1
    ## 6/14 respiratory chain complex IV                                               firebrick1
    ## 8/26 cytochrome complex                                                         firebrick1
    ## 51/495 mitochondrial membrane                                                   firebrick1
    ## 131/1373 organelle membrane                                                     firebrick1
    ## 39/386 organelle inner membrane                                                 firebrick1
    ## 2/11 respiratory chain complex III                                              firebrick1
    ## 2/19 proton-transporting two-sector ATPase complex, proton-transporting domain  firebrick1
    ## 2/10 proton-transporting ATP synthase complex, coupling factor F(o)             firebrick1
    ## 5/29 proton-transporting two-sector ATPase complex                              firebrick1
    ## 125/1424 mitochondrion                                                          firebrick1
    ## 4/16 mitochondrial proton-transporting ATP synthase complex                     firebrick1
    ##                                                                                enriched
    ## 2/20 proteasome core complex                                                       DISS
    ## 61/675 nucleolus                                                                   HOMO
    ## 60/208 ribosome                                                                    DISS
    ## 35/106 large ribosomal subunit                                                     DISS
    ## 33/53 cytosolic large ribosomal subunit                                            DISS
    ## 22/68 small ribosomal subunit                                                      DISS
    ## 20/41 cytosolic small ribosomal subunit                                            DISS
    ## 11/153 ubiquitin ligase complex                                                    HOMO
    ## 4/15 protein serine/threonine phosphatase complex                                  HOMO
    ## 20/303 transferase complex                                                         HOMO
    ## 5/63 protein acetyltransferase complex                                             HOMO
    ## 55/645 catalytic complex                                                           HOMO
    ## 9/74 methyltransferase complex                                                     HOMO
    ## 1/10 Set1C/COMPASS complex                                                         HOMO
    ## 2/59 protein kinase complex                                                        HOMO
    ## 1/8 IkappaB kinase complex                                                         HOMO
    ## 63/194 cytosolic part                                                              DISS
    ## 6/112 chromosome, telomeric region                                                 HOMO
    ## 36/360 nuclear chromosome part                                                     HOMO
    ## 1/14 Sin3-type complex                                                             HOMO
    ## 1/18 centromeric heterochromatin                                                   HOMO
    ## 31/262 chromatin                                                                   HOMO
    ## 3/71 heterochromatin                                                               HOMO
    ## 18/331 chromosome                                                                  HOMO
    ## 4/76 chromosome, centromeric region                                                HOMO
    ## 53/592 chromosomal part                                                            HOMO
    ## 5/46 transcription elongation factor complex                                       HOMO
    ## 4/71 PML body                                                                      HOMO
    ## 38/529 nucleoplasm part                                                            HOMO
    ## 2/14 mRNA cleavage and polyadenylation specificity factor complex                  HOMO
    ## 2/19 mRNA cleavage factor complex                                                  HOMO
    ## 6/99 nuclear speck                                                                 HOMO
    ## 1/34 Cajal body                                                                    HOMO
    ## 16/238 nuclear body                                                                HOMO
    ## 122/957 extracellular space                                                        DISS
    ## 101/985 extracellular region                                                       DISS
    ## 23/247 proteinaceous extracellular matrix                                          DISS
    ## 51/384 extracellular matrix                                                        DISS
    ## 113/1046 Golgi apparatus                                                           HOMO
    ## 0/9 Golgi transport complex                                                        HOMO
    ## 39/272 side of membrane                                                            DISS
    ## 0/11 immunoglobulin complex                                                        DISS
    ## 4/42 cell periphery                                                                DISS
    ## 61/381 plasma membrane region                                                      DISS
    ## 4/17 MHC protein complex                                                           DISS
    ## 4/10 MHC class I protein complex                                                   DISS
    ## 219/1505 membrane-bounded vesicle                                                  DISS
    ## 174/1126 extracellular membrane-bounded organelle                                  DISS
    ## 8/30 integrin complex                                                              DISS
    ## 89/507 cell surface                                                                DISS
    ## 5/12 phagocytic cup                                                                DISS
    ## 88/590 cell-cell junction                                                          DISS
    ## 115/628 anchoring junction                                                         DISS
    ## 87/356 cell-substrate junction                                                     DISS
    ## 2/44 spindle microtubule                                                           HOMO
    ## 8/102 ciliary basal body                                                           HOMO
    ## 1/46 spindle pole                                                                  HOMO
    ## 151/1308 cytoskeletal part                                                         HOMO
    ## 0/14 gamma-tubulin complex                                                         HOMO
    ## 40/516 microtubule organizing center                                               HOMO
    ## 138/1238 cytoskeleton                                                              HOMO
    ## 6/92 centriole                                                                     HOMO
    ## 7/130 microtubule organizing center part                                           HOMO
    ## 0/10 mRNA cap binding complex                                                      HOMO
    ## 5/19 sarcoplasmic reticulum membrane                                               HOMO
    ## 44/250 transporter complex                                                         HOMO
    ## 25/129 cation channel complex                                                      HOMO
    ## 23/161 site of polarized growth                                                    HOMO
    ## 2/8 ribbon synapse                                                                 HOMO
    ## 8/49 neuron projection membrane                                                    HOMO
    ## 123/864 cell projection part                                                       HOMO
    ## 40/257 axon part                                                                   HOMO
    ## 6/12 axon initial segment                                                          HOMO
    ## 0/16 inhibitory synapse                                                            HOMO
    ## 3/15 asymmetric synapse                                                            HOMO
    ## 69/506 synapse                                                                     HOMO
    ## 39/237 synaptic membrane                                                           HOMO
    ## 42/210 postsynaptic density                                                        HOMO
    ## 83/522 synapse part                                                                HOMO
    ## 25/130 neuron spine                                                                HOMO
    ## 11/58 dendritic shaft                                                              HOMO
    ## 68/413 dendrite                                                                    HOMO
    ## 189/1239 neuron part                                                               HOMO
    ## 6/76 pore complex                                                                  HOMO
    ## 2/18 nuclear periphery                                                             HOMO
    ## 51/347 vacuole                                                                     DISS
    ## 0/7 TORC1 complex                                                                  HOMO
    ## 28/175 myelin sheath                                                               DISS
    ## 10/87 oxidoreductase complex                                                       DISS
    ## 24/170 mitochondrial membrane part                                                 DISS
    ## 6/42 respiratory chain complex I                                                   DISS
    ## 67/725 mitochondrial part                                                          DISS
    ## 6/14 respiratory chain complex IV                                                  DISS
    ## 8/26 cytochrome complex                                                            DISS
    ## 51/495 mitochondrial membrane                                                      DISS
    ## 131/1373 organelle membrane                                                        DISS
    ## 39/386 organelle inner membrane                                                    DISS
    ## 2/11 respiratory chain complex III                                                 DISS
    ## 2/19 proton-transporting two-sector ATPase complex, proton-transporting domain     DISS
    ## 2/10 proton-transporting ATP synthase complex, coupling factor F(o)                DISS
    ## 5/29 proton-transporting two-sector ATPase complex                                 DISS
    ## 125/1424 mitochondrion                                                             DISS
    ## 4/16 mitochondrial proton-transporting ATP synthase complex                        DISS

    gomwuPlot(input,goAnnotations,goDivision,
        absValue=-log(0.05,10),  
        level1=0.00000001, 
        level2=0.00000001, 
        level3=0.000000001, 
        txtsize=1.4,    
        treeHeight=0.5, 
      #colors=c("#d9d9d9","#525252","#d9d9d9","#525252")
        #colors=c("blue","green","blue","green") 
        colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") 
    )

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

![](../../figures/05_GO_MMU/CC-2.png)

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
    ##                                                         color enriched
    ## 16/238 nuclear body                               dodgerblue2     HOMO
    ## 38/529 nucleoplasm part                           dodgerblue2     HOMO
    ## 60/208 ribosome                                    firebrick1     DISS
    ## 35/106 large ribosomal subunit                     firebrick1     DISS
    ## 33/53 cytosolic large ribosomal subunit            firebrick1     DISS
    ## 63/194 cytosolic part                              firebrick1     DISS
    ## 22/68 small ribosomal subunit                      firebrick1     DISS
    ## 20/41 cytosolic small ribosomal subunit            firebrick1     DISS
    ## 174/1126 extracellular membrane-bounded organelle  firebrick1     DISS
    ## 115/628 anchoring junction                         firebrick1     DISS
    ## 87/356 cell-substrate junction                     firebrick1     DISS
    ## 83/522 synapse part                               dodgerblue2     HOMO
    ## 101/985 extracellular region                       firebrick1     DISS
    ## 122/957 extracellular space                        firebrick1     DISS

table of significance, top 10 ish Go terms
------------------------------------------
