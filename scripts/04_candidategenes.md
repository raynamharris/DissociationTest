    library(readxl)
    library(ggplot2)
    library(cowplot)
    library(dplyr)
    library(ggrepel)
    library(plyr)

    # set output file for figures and general chunk settings
    knitr::opts_chunk$set(fig.path = '../figures/04_candidategenes/', echo = T, message = F, warning=F)

Comparison with other candidate gen lists
=========================================

1.  Sanes and Licthman 1999 analysis
2.  Cho et al 2015
3.  Cembrowksi

Sanes and licthman
------------------

Import file with information about differentially expressed genes (aka
`dissociation`) and the list of genes found in the reference
transcriptome (aka `geneids`).

    dissociation <- read.csv("../results/volcanoTreatment.csv", header = T, row.names = 1)
    geneids <- read.csv("../data/geneids.csv", header = T)

    # Filter out non-significant genes, sort by p-value, rename column, round to 2 decimal places. 

    DEGs <- dissociation %>%
      filter(direction != "none") %>%
      arrange((padj))

    # Compare to Molecules implicated in hippocampal LTP
    # create list of candidate genes Sanes and Lichtman 1999 

    sanesLichtman <- c("Gria1", "Gria2", 
           "Grm1", "Grm4", "Grm5", "Grm7",
           "Grin1", "Grin2a", "Grin2d", 
           "Th", "Drd1",
           "Adrb1", "Adrb2", "Adrb3",
           "Adra1a", "Adra1b", "Adra1d", 
           "Adra2a", "Adra2b", "Adra2c",
           "Oprm1", "Oprk1", "Oprd1",
           "Chrm1", "Chrm2", "Chrm3", "Chrm4", "Chrm5",
           "Chrna1", "Chrna7", "Chrna3", 
           "Chrnb1", "Chrnb2", "Chrnb3",
           "Gabra1", "Gabra2",  "Gabra3", "Gabra5", "Gabra6",
           "Gabrb1", "Gabrb2", "Gabrb3",  "Gabrr1",  "Gabbr1",
           "Cnr1", "Pnoc", "Oprl1",
           "Htr1a", "Htr1b", "Htr1f",
           "Htr2a", "Htr2c", "Htr2b",
           "Htr3a", "Htr3b", "Htr5a", "Htr5b",
           "Htr7", "Htr6", "Htr4", 
           "Edn1", "Egf", "Fgf2",
           "Nrg1", "Nrg2", "Nrg3",
           "Erbb4", "Ngf", "Bdnf", "Ntrk2",
           "Nos1", "Nos2", "Nos3",
           "Il1b", # interleukin 1 beta
           "Inhba", "Calm1", "Calm2", "Calm3",
           "Nrgn", "Calb1", "Calb2", "Gap43", "S100b",
           "Cacna1c", "Cacna1d", "Cacna1s", "Cacna1f",
           "Cacna1b", "Cacna1a", "Cacna1e",
           "Cnga2", "Syp", "Napa",
           "Vamp1", "Vamp2", "Vamp3", "Vamp4", "Vamp5", "Vamp8",
           "Rab3a", "Stx1b", "Syn1", "Snap25",
           "Dlg4", "Rara", "Rarb", "Creb1",
           "Egr1", "Egr2", "Epha5", "Efna5", "Ncam1", "Ncam2",
           "Cdh1", "Cdh2", "Thy1", "Icam5",
           "L1cam", "Ptn", 
            "Itga1", "Itga10","Itga11", "Itga2",  "Itga2b",
           "Itga3",  "Itga4", "Itga5",  "Itga6",  "Itga7",  
           "Itga8",  "Itga9",  "Itgad",  "Itgae", "Itgal",  
           "Itgam",  "Itgav",  "Itgax",  "Itgb1",  "Itgb1bp1",
           "Itgb2",  "Itgb2l", "Itgb3",  "Itgb3bp",  "Itgb4",
           "Itgb5",  "Itgb6", "Itgb7",  "Itgb8",  "Itgbl1",
           "Cd47", "Tnc",
           "Mapk1", "Mapk10", "Mapk11", "Mapk12", "Mapk14", 
           "Mapk3", "Mapk4", "Mapk6", "Mapk7", "Mapk8", "Mapk9",
           "Src", "Fyn", 
           "Prkcg", "Prkg1", "Prkcz", 
           "Capn1", "Capn10", "Capn11", "Capn12", "Capn13",
           "Capn15", "Capn2", "Capn3", "Capn5", "Capn6", 
           "Capn7", "Capn8", "Capn9",
           "Cast", "Serpine2", "Plat", "Plg", "Ube3a",
           "Pla2g10", "Pla2g12a", "Pla2g12b", "Pla2g15",
           "Pla2g16", "Pla2g1b", "Pla2g2a", "Pla2g2c", "Pla2g2d",
           "Pla2g2e", "Pla2g2f", "Pla2g3",  "Pla2g4a", "Pla2g4b",
           "Pla2g4e", "Pla2g4f", "Pla2g5", "Pla2g6", "Pla2g7",
           "Plcb1", "Plcb2", "Plcb3", "Plcb4",
           "Plcg1", "Plcg2",
           "Parp1", "Ppp3ca", "Ppp3cb", "Ppp3cc",
           "Phpt1", "Ache",
           "Adcy1", "Adcy10", "Adcy2", "Adcy3", "Adcy4",
           "Adcy5", "Adcy6", "Adcy7", "Adcy8", "Adcy9",
           "Gucy1a2", "Gucy1a3", "Gucy1b2", "Gucy1b3",
           "Gucy2c", "Gucy2d", "Gucy2e", "Gucy2g",
           "Sptan1", "Sptbn1", "Gfap", "Stmn4",
           "Ccr7", "Mas1",
           "Homer1", "Homer2", "Homer3" )

    sanesLichtman[order(sanesLichtman)] # print list alphabetically

    ##   [1] "Ache"     "Adcy1"    "Adcy10"   "Adcy2"    "Adcy3"    "Adcy4"   
    ##   [7] "Adcy5"    "Adcy6"    "Adcy7"    "Adcy8"    "Adcy9"    "Adra1a"  
    ##  [13] "Adra1b"   "Adra1d"   "Adra2a"   "Adra2b"   "Adra2c"   "Adrb1"   
    ##  [19] "Adrb2"    "Adrb3"    "Bdnf"     "Cacna1a"  "Cacna1b"  "Cacna1c" 
    ##  [25] "Cacna1d"  "Cacna1e"  "Cacna1f"  "Cacna1s"  "Calb1"    "Calb2"   
    ##  [31] "Calm1"    "Calm2"    "Calm3"    "Capn1"    "Capn10"   "Capn11"  
    ##  [37] "Capn12"   "Capn13"   "Capn15"   "Capn2"    "Capn3"    "Capn5"   
    ##  [43] "Capn6"    "Capn7"    "Capn8"    "Capn9"    "Cast"     "Ccr7"    
    ##  [49] "Cd47"     "Cdh1"     "Cdh2"     "Chrm1"    "Chrm2"    "Chrm3"   
    ##  [55] "Chrm4"    "Chrm5"    "Chrna1"   "Chrna3"   "Chrna7"   "Chrnb1"  
    ##  [61] "Chrnb2"   "Chrnb3"   "Cnga2"    "Cnr1"     "Creb1"    "Dlg4"    
    ##  [67] "Drd1"     "Edn1"     "Efna5"    "Egf"      "Egr1"     "Egr2"    
    ##  [73] "Epha5"    "Erbb4"    "Fgf2"     "Fyn"      "Gabbr1"   "Gabra1"  
    ##  [79] "Gabra2"   "Gabra3"   "Gabra5"   "Gabra6"   "Gabrb1"   "Gabrb2"  
    ##  [85] "Gabrb3"   "Gabrr1"   "Gap43"    "Gfap"     "Gria1"    "Gria2"   
    ##  [91] "Grin1"    "Grin2a"   "Grin2d"   "Grm1"     "Grm4"     "Grm5"    
    ##  [97] "Grm7"     "Gucy1a2"  "Gucy1a3"  "Gucy1b2"  "Gucy1b3"  "Gucy2c"  
    ## [103] "Gucy2d"   "Gucy2e"   "Gucy2g"   "Homer1"   "Homer2"   "Homer3"  
    ## [109] "Htr1a"    "Htr1b"    "Htr1f"    "Htr2a"    "Htr2b"    "Htr2c"   
    ## [115] "Htr3a"    "Htr3b"    "Htr4"     "Htr5a"    "Htr5b"    "Htr6"    
    ## [121] "Htr7"     "Icam5"    "Il1b"     "Inhba"    "Itga1"    "Itga10"  
    ## [127] "Itga11"   "Itga2"    "Itga2b"   "Itga3"    "Itga4"    "Itga5"   
    ## [133] "Itga6"    "Itga7"    "Itga8"    "Itga9"    "Itgad"    "Itgae"   
    ## [139] "Itgal"    "Itgam"    "Itgav"    "Itgax"    "Itgb1"    "Itgb1bp1"
    ## [145] "Itgb2"    "Itgb2l"   "Itgb3"    "Itgb3bp"  "Itgb4"    "Itgb5"   
    ## [151] "Itgb6"    "Itgb7"    "Itgb8"    "Itgbl1"   "L1cam"    "Mapk1"   
    ## [157] "Mapk10"   "Mapk11"   "Mapk12"   "Mapk14"   "Mapk3"    "Mapk4"   
    ## [163] "Mapk6"    "Mapk7"    "Mapk8"    "Mapk9"    "Mas1"     "Napa"    
    ## [169] "Ncam1"    "Ncam2"    "Ngf"      "Nos1"     "Nos2"     "Nos3"    
    ## [175] "Nrg1"     "Nrg2"     "Nrg3"     "Nrgn"     "Ntrk2"    "Oprd1"   
    ## [181] "Oprk1"    "Oprl1"    "Oprm1"    "Parp1"    "Phpt1"    "Pla2g10" 
    ## [187] "Pla2g12a" "Pla2g12b" "Pla2g15"  "Pla2g16"  "Pla2g1b"  "Pla2g2a" 
    ## [193] "Pla2g2c"  "Pla2g2d"  "Pla2g2e"  "Pla2g2f"  "Pla2g3"   "Pla2g4a" 
    ## [199] "Pla2g4b"  "Pla2g4e"  "Pla2g4f"  "Pla2g5"   "Pla2g6"   "Pla2g7"  
    ## [205] "Plat"     "Plcb1"    "Plcb2"    "Plcb3"    "Plcb4"    "Plcg1"   
    ## [211] "Plcg2"    "Plg"      "Pnoc"     "Ppp3ca"   "Ppp3cb"   "Ppp3cc"  
    ## [217] "Prkcg"    "Prkcz"    "Prkg1"    "Ptn"      "Rab3a"    "Rara"    
    ## [223] "Rarb"     "S100b"    "Serpine2" "Snap25"   "Sptan1"   "Sptbn1"  
    ## [229] "Src"      "Stmn4"    "Stx1b"    "Syn1"     "Syp"      "Th"      
    ## [235] "Thy1"     "Tnc"      "Ube3a"    "Vamp1"    "Vamp2"    "Vamp3"   
    ## [241] "Vamp4"    "Vamp5"    "Vamp8"

    # confirm that all all Sanes and Lichtman genes are in the reference transcriptome
    sanesLichtman_reference <- geneids %>%
      dplyr::filter(gene %in% sanesLichtman) %>%
      dplyr::select(gene) %>%
      distinct() %>%
      arrange(gene) %>%
      droplevels()
    sanesLichtman_reference <- sanesLichtman_reference[,c(1)]
    str(sanesLichtman_reference)

    ##  Factor w/ 243 levels "Ache","Adcy1",..: 1 2 3 4 5 6 7 8 9 10 ...

    # identify which of the Sanes and Lichtman genes are present in my samples
    sanesLichtman_present <- dissociation %>%
      dplyr::filter(gene %in% sanesLichtman) %>%
      droplevels()
    sanesLichtman_present <- sanesLichtman_present[,c(1)]
    sanesLichtman_present

    ##   [1] Ache     Adcy1    Adcy2    Adcy3    Adcy5    Adcy6    Adcy7   
    ##   [8] Adcy8    Adcy9    Adra1a   Adra1d   Adra2a   Adra2b   Adra2c  
    ##  [15] Adrb1    Bdnf     Cacna1a  Cacna1b  Cacna1c  Cacna1d  Cacna1e 
    ##  [22] Calb1    Calb2    Calm1    Calm2    Calm3    Capn1    Capn10  
    ##  [29] Capn15   Capn2    Capn3    Capn5    Capn7    Cd47     Cdh2    
    ##  [36] Chrm1    Chrm2    Chrm3    Chrm4    Chrna7   Chrnb2   Cnr1    
    ##  [43] Creb1    Dlg4     Egr1     Egr2     Epha5    Erbb4    Fgf2    
    ##  [50] Fyn      Gabbr1   Gabra1   Gabra2   Gabra3   Gabra5   Gabrb1  
    ##  [57] Gabrb2   Gabrb3   Gap43    Gfap     Gria1    Gria2    Grin1   
    ##  [64] Grin2a   Grin2d   Grm1     Grm4     Grm5     Grm7     Gucy1a2 
    ##  [71] Gucy1a3  Gucy1b3  Gucy2e   Gucy2g   Homer1   Homer2   Homer3  
    ##  [78] Htr1a    Htr1b    Htr2a    Htr2c    Htr3a    Htr4     Htr5a   
    ##  [85] Htr5b    Htr6     Htr7     Icam5    Il1b     Itga1    Itga10  
    ##  [92] Itga11   Itga2b   Itga3    Itga4    Itga5    Itga6    Itga7   
    ##  [99] Itga8    Itga9    Itgam    Itgav    Itgb1    Itgb1bp1 Itgb2   
    ## [106] Itgb4    Itgb5    Itgb6    Itgb8    Itgbl1   L1cam    Mapk1   
    ## [113] Mapk10   Mapk11   Mapk12   Mapk14   Mapk3    Mapk4    Mapk6   
    ## [120] Mapk7    Mapk8    Mapk9    Mas1     Napa     Ncam1    Ncam2   
    ## [127] Nos1     Nos3     Nrg1     Nrg2     Nrg3     Nrgn     Ntrk2   
    ## [134] Oprd1    Oprl1    Parp1    Phpt1    Pla2g12a Pla2g15  Pla2g16 
    ## [141] Pla2g3   Pla2g4e  Pla2g6   Pla2g7   Plat     Plcb1    Plcb3   
    ## [148] Plcb4    Plcg1    Plcg2    Pnoc     Ppp3ca   Ppp3cb   Ppp3cc  
    ## [155] Prkcg    Prkcz    Prkg1    Ptn      Rab3a    Rara     S100b   
    ## [162] Serpine2 Snap25   Sptan1   Sptbn1   Src      Stmn4    Stx1b   
    ## [169] Syn1     Syp      Th       Thy1     Ube3a    Vamp1    Vamp2   
    ## [176] Vamp3    Vamp4    Vamp8   
    ## 178 Levels: Ache Adcy1 Adcy2 Adcy3 Adcy5 Adcy6 Adcy7 Adcy8 Adcy9 ... Vamp8

    # identify whichof the Sanes and Lichtman genes are differentially expressed in this analysis
    sanesLichtman_DEGs <- DEGs %>%
      dplyr::filter(gene %in% sanesLichtman) %>%
      arrange(gene)
    sanesLichtman_DEGs

    ##     gene   pvalue       lfc        padj direction
    ## 1 Grin2a 1.570140 -1.659562 0.026906675      HOMO
    ## 2   Il1b 1.523350  2.405914 0.029967447      DISS
    ## 3  Itga5 1.795090  3.054466 0.016029141      DISS
    ## 4  Itgam 1.304370  1.746838 0.049616993      DISS
    ## 5  Itgb5 1.695894  1.978733 0.020142152      DISS
    ## 6  Mapk3 2.008779  1.606075 0.009799875      DISS

    # Grin2a NMDAR 2A
    # Il1b   Interleukin 1 Beta
    # Itga5  Integrin Subunit Alpha 5
    # Itgam  Integrin Subunit Alpha M
    # Itgb5  Integrin Subunit Beta 5
    # Mapk3  MAP Kinase 3

    # percent DEGs in the sanes lichtman list
    round(6/243*100,2)

    ## [1] 2.47

    round(6/178*100,2)

    ## [1] 3.37

Cho et al 2015 anlaysis
-----------------------

Cho et al 2015 used RNA sequencing to quantify transcript levels in the
mouse hippocampus after contextual fear conditioning. The Cho dataset
provides a snapshot of gene expression changes associated with
hippocampal learning and memory 30 min and 4 h after an experiment. The
Cho data are available at
<http://science.sciencemag.org/content/suppl/2015/09/30/350.6256.82.DC1>
and <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72064>. The
file `../data/aac7368-Cho.SM.Table.S2.xls` of differentially expressed
genes was used as a representative dataset for learning and memory
associated gene expression patterns.

In this analysis, I compared the Cho et al differentially expressed
genes (DEGs) to my experimental results (referred to as the Harris data)
to identify the interaction between genes that are differentially
expressed following fear condition and following a chemical
manipulation.

This analysis prints the list of genes that are differentially expressed
in both experiments. The images show that only a few (red dots) of genes
that respond to chemical dissociation are up-regulated or down-regulated
following fear-conditioning.

The number of differential gene expression in the harris data set and
then the Cho data sets at 4 h and 30 min.

    # read data and set factors
    dissociation <- read.csv("../results/volcanoTreatment.csv", header = T, row.names = 1)
    dissociation$direction <- factor(dissociation$direction, c("HOMO", "none", "DISS"))

In the Cho et al. data, there are 9 differentially expressed genes after
4 hours with LFC &gt; 1. But there are 35 with lfc &gt; 0.25. I will go
with those since that is about the cuttoff used in the Cho paper.

    S2 <- as.data.frame(readxl::read_excel("../data/aac7368-Cho.SM.Table.S2.xls", skip = 1 ))

    fourhoursRNA <- rename(S2, c(`RNA fold change (4 h/control), log2` ="lfc", 
                       `p-value (4 h)` = "pvalue",
                       `Gene Symbol` = "gene"))
    fourhoursRNA <- wrangleCho(fourhoursRNA)

    volcanoplot1 <- plotvolcano(fourhoursRNA, fourhoursRNA$lfc, fourhoursRNA$log10p, 
                                plottitle = "Cho DEGs - 4 h") +
                    scale_color_manual(values = c("none" = "grey",
                                    "fear-conditioned" = "#018571",
                                    "control" = "#a6611a"))
    volcanoplot1

![](../figures/04_candidategenes/fourhours-1.png)

    Chofourhoursgenes <- fourhoursRNA %>%
      filter(direction != "none")
    Chofourhoursgenes

    ##             gene        lfc    log10p       pvalue        direction
    ## 1          Npas4  0.4155533  1.848876 1.416199e-02 fear-conditioned
    ## 2            Fos -0.6719830  6.157898 6.951882e-07          control
    ## 3          Vat1l -0.4320860  5.258212 5.518076e-06          control
    ## 4          Capsl -0.4687226  3.060348 8.702668e-04          control
    ## 5         Col8a2 -0.7632370  7.260695 5.486624e-08          control
    ## 6        Wfikkn2 -0.6316636  2.824650 1.497442e-03          control
    ## 7          Lama5 -0.5116774  3.934160 1.163697e-04          control
    ## 8         Tmem72 -1.6548897  7.715675 1.924530e-08          control
    ## 9           Pon3 -0.5167425  5.414085 3.854033e-06          control
    ## 10       Sostdc1 -0.8069617 14.441635 3.617140e-15          control
    ## 11          Prlr -0.6768509  8.810605 1.546661e-09          control
    ## 12        Sema3b -0.5417905  3.582509 2.615114e-04          control
    ## 13          Otx2 -0.6570360  5.026343 9.411464e-06          control
    ## 14        Pcdhb2  0.5429030  2.381734 4.152078e-03 fear-conditioned
    ## 15          Epn3 -0.9624654  4.451133 3.538889e-05          control
    ## 16       Ccdc135 -0.9159923  5.133405 7.355206e-06          control
    ## 17          Aqp1 -1.3821718 12.836310 1.457775e-13          control
    ## 18         Sulf1 -0.4668857 11.231934 5.862272e-12          control
    ## 19         Cldn2 -1.0033125  8.527052 2.971308e-09          control
    ## 20        Calml4 -1.0874741  5.187103 6.499753e-06          control
    ## 21         Xrcc6 -0.4408808  2.515905 3.048558e-03          control
    ## 22         Rsph1 -0.7219688  1.560257 2.752597e-02          control
    ## 23       Tinagl1 -0.4192463  1.924291 1.190444e-02          control
    ## 24         Clic6 -1.0178515  7.843026 1.435405e-08          control
    ## 25         Kcne2 -1.1696214 12.652911 2.223763e-13          control
    ## 26        Slc4a5 -1.2353303  6.039947 9.121213e-07          control
    ## 27          Mid1 -0.4123224  1.655534 2.210376e-02          control
    ## 28        Gpr133 -0.7087897  2.246372 5.670590e-03          control
    ## 29         Serhl -0.6882239  2.088772 8.151321e-03          control
    ## 30       Efcab10 -0.5056588  1.694729 2.019625e-02          control
    ## 31            F8  0.4094271  1.653764 2.219402e-02 fear-conditioned
    ## 32           Mia -0.4531486  4.740758 1.816527e-05          control
    ## 33           Ace -0.4773821  3.898966 1.261927e-04          control
    ## 34            F5 -0.7884346  6.388467 4.088205e-07          control
    ## 35         Kcnj2  0.5207135  2.697591 2.006360e-03 fear-conditioned
    ## 36         Sgms2 -0.6627646  2.154366 7.008651e-03          control
    ## 37         Folr1 -1.1021230 12.302898 4.978535e-13          control
    ## 38       Slc24a5 -0.5351884  6.606556 2.474252e-07          control
    ## 39        Slc4a2 -0.4088131  7.232030 5.860972e-08          control
    ## 40           Ttr -1.3565859  5.464636 3.430551e-06          control
    ## 41         Enpp2 -0.5505385  4.584654 2.602234e-05          control
    ## 42         Mtcp1  0.4068933  2.223980 5.970628e-03 fear-conditioned
    ## 43            Kl -0.5682186 11.045384 9.007737e-12          control
    ## 44           Lbp -0.4462483  8.124935 7.500061e-09          control
    ## 45 1500015O10Rik -0.8228514 13.169578 6.767410e-14          control
    ## 46        Kcnj13 -0.5225731  7.976999 1.054388e-08          control
    ## 47        Sowahc -0.4104584  1.544910 2.851609e-02          control
    ## 48          Mrc1  0.4931112  1.347372 4.493946e-02 fear-conditioned
    ## 49       Slc2a12 -0.4052413  5.427146 3.739847e-06          control
    ## 50         Abca4 -0.9014487  9.452195 3.530243e-10          control
    ## 51         Mdfic -0.7068837  1.741788 1.812224e-02          control
    ## 52          Msx1 -0.8377806  4.101077 7.923613e-05          control
    ## 53        Atp2a3 -0.4401118  2.139664 7.249959e-03          control
    ## 54      Itpripl1 -0.5435444  2.733488 1.847193e-03          control

    Cho4hrCandidates <- Chofourhoursgenes$gene

    Cho4hrCandidates

    ##  [1] "Npas4"         "Fos"           "Vat1l"         "Capsl"        
    ##  [5] "Col8a2"        "Wfikkn2"       "Lama5"         "Tmem72"       
    ##  [9] "Pon3"          "Sostdc1"       "Prlr"          "Sema3b"       
    ## [13] "Otx2"          "Pcdhb2"        "Epn3"          "Ccdc135"      
    ## [17] "Aqp1"          "Sulf1"         "Cldn2"         "Calml4"       
    ## [21] "Xrcc6"         "Rsph1"         "Tinagl1"       "Clic6"        
    ## [25] "Kcne2"         "Slc4a5"        "Mid1"          "Gpr133"       
    ## [29] "Serhl"         "Efcab10"       "F8"            "Mia"          
    ## [33] "Ace"           "F5"            "Kcnj2"         "Sgms2"        
    ## [37] "Folr1"         "Slc24a5"       "Slc4a2"        "Ttr"          
    ## [41] "Enpp2"         "Mtcp1"         "Kl"            "Lbp"          
    ## [45] "1500015O10Rik" "Kcnj13"        "Sowahc"        "Mrc1"         
    ## [49] "Slc2a12"       "Abca4"         "Mdfic"         "Msx1"         
    ## [53] "Atp2a3"        "Itpripl1"

    Cho4hrCandidatesGrepl <- c("Fos|Col8a|Wfikkn2|Lama|Tmem|Pon|Sostdc|Prlr|Sema|Otx|Pcdhb|Epn|Ccdc|Aqp|Cldn|Calml|Rsph|Clic|Kcne|Slc4|Gpr|Serhl|Efcab|F5|Kcnj|Sgms|Folr|Slc24a|Ttr|Enpp|Kl|1500015O10Rik|Kcnj|Abca|Mdfic|Msx|Itpripl")


    # see if cho 4 hour DEGs are in my DEG list
    Cho_DEGs <- DEGs %>%
      dplyr::filter(gene %in% Cho4hrCandidates) %>%
      arrange(gene)
    Cho_DEGs

    ##    gene   pvalue      lfc       padj direction
    ## 1 Enpp2 1.310489 1.768966 0.04892278      DISS

    # see if related cho 4 hour DEGs are in my DEG list
    Cho_DEGs <- DEGs %>%
      dplyr::filter(grepl(Cho4hrCandidatesGrepl, gene)) %>%
      arrange(gene)
    Cho_DEGs

    ##       gene   pvalue       lfc         padj direction
    ## 1   Cldn11 3.979408  3.139729 0.0001048558      DISS
    ## 2    Clic4 2.742689  2.008972 0.0018084688      DISS
    ## 3    Enpp2 1.310489  1.768966 0.0489227752      DISS
    ## 4     Fosb 1.478562  1.585131 0.0332229625      DISS
    ## 5    Gpr34 1.310489  2.079302 0.0489227752      DISS
    ## 6    Gpr84 1.728492  2.858254 0.0186856271      DISS
    ## 7    Kcnj6 1.359987 -2.392068 0.0436529130      HOMO
    ## 8   Sema5a 2.526452  3.251563 0.0029754201      DISS
    ## 9  Tmem119 3.022088  2.543998 0.0009504111      DISS
    ## 10 Tmem170 1.532416  2.526577 0.0293483740      DISS
    ## 11 Tmem88b 4.169077  3.141197 0.0000677521      DISS

At there Fosb is the only overlap.

    S2 <- as.data.frame(readxl::read_excel("../data/aac7368-Cho.SM.Table.S2.xls", skip = 1 ))

    thirtyminRNA <- rename(S2, c(`RNA fold change (30 min/control), log2` ="lfc", 
                       `p-value (30 min)` = "pvalue",
                       `Gene Symbol` = "gene"))
    thirtyminRNA <- wrangleCho(thirtyminRNA)

    volcanoplot3 <- plotvolcano(thirtyminRNA, thirtyminRNA$lfc, thirtyminRNA$log10p, 
                                plottitle = "Cho DEGs - 30 min")  +
                    scale_color_manual(values = c("none" = "grey",
                                    "fear-conditioned" = "#018571",
                                    "control" = "#a6611a"))
    volcanoplot3

![](../figures/04_candidategenes/thirtymin-1.png)

    Cho30mingenes <- thirtyminRNA %>%
      filter(direction != "none")
    Cho30mingenes

    ##             gene        lfc    log10p       pvalue        direction
    ## 1          Npas4  0.6768515 14.827528 1.487551e-15 fear-conditioned
    ## 2            Fos  0.8341890 26.810336 1.547619e-27 fear-conditioned
    ## 3           Junb  0.4452631  9.499079 3.168989e-10 fear-conditioned
    ## 4            Arc  0.4472752 18.742271 1.810211e-19 fear-conditioned
    ## 5           Btg2  0.4155126  2.254590 5.564291e-03 fear-conditioned
    ## 6          Arl4d  0.4605700  2.976364 1.055933e-03 fear-conditioned
    ## 7          Capsl -0.4169848  2.884578 1.304435e-03          control
    ## 8         Col8a2 -0.4218772  5.566678 2.712201e-06          control
    ## 9          Acta2 -0.4334486  1.353017 4.435918e-02          control
    ## 10        Tmem72 -0.4520730  7.531901 2.938323e-08          control
    ## 11         Ackr3  0.4480152  2.289752 5.131548e-03 fear-conditioned
    ## 12          Prlr -0.4379552  5.407868 3.909598e-06          control
    ## 13         Cyr61  0.5767257  4.346678 4.501134e-05 fear-conditioned
    ## 14          Epn3 -0.7182916  4.936213 1.158209e-05          control
    ## 15          Aqp1 -0.5998278  6.961992 1.091461e-07          control
    ## 16         Cldn2 -0.4630600 10.258710 5.511754e-11          control
    ## 17        Calml4 -0.5080975  5.666329 2.156111e-06          control
    ## 18        Col9a3 -0.4282895  2.348704 4.480190e-03          control
    ## 19       Tinagl1 -0.4740280  2.305269 4.951430e-03          control
    ## 20         Clic6 -0.4399658  2.641955 2.280581e-03          control
    ## 21          Nqo1 -0.4859055  2.570114 2.690831e-03          control
    ## 22          Aass -0.4322532  1.557089 2.772750e-02          control
    ## 23         Hspg2 -0.4206140  1.534223 2.922650e-02          control
    ## 24         Nr4a1  0.4076412 17.734942 1.841018e-18 fear-conditioned
    ## 25            F5 -0.5154348  3.029245 9.348772e-04          control
    ## 26         Ltc4s -0.7486669  5.066732 8.575678e-06          control
    ## 27         Folr1 -0.5512630  6.888991 1.291246e-07          control
    ## 28           Ttr -0.4649020  1.921448 1.198262e-02          control
    ## 29    St6galnac2 -0.6133743  3.244462 5.695576e-04          control
    ## 30            Kl -0.4362041  5.450433 3.544599e-06          control
    ## 31 1500015O10Rik -0.6183895  2.611249 2.447658e-03          control
    ## 32       Slc2a12 -0.4092759  4.637800 2.302503e-05          control
    ## 33          Fosb  0.4801394 12.831460 1.474143e-13 fear-conditioned
    ## 34          Msx1 -0.5967465  3.049871 8.915150e-04          control

    Cho30minCandidates <- Cho30mingenes$gene

    Cho30minCandidatesGrepl <- c("Npas|Fos|Cyr|Epn|Aqp|Calml|F5|Ltc|Folr|St6galnac|1500015O10Rik|Msx1")

    # see if cho 4 hour DEGs are in my DEG list
    Cho_DEGs <- DEGs %>%
      dplyr::filter(gene %in% Cho30minCandidates) %>%
      arrange(gene)
    Cho_DEGs

    ##   gene   pvalue      lfc       padj direction
    ## 1 Fosb 1.478562 1.585131 0.03322296      DISS

    # see if related cho 4 hour DEGs are in my DEG list
    Cho_DEGs <- DEGs %>%
      dplyr::filter(grepl(Cho30minCandidatesGrepl, gene)) %>%
      arrange(gene)
    Cho_DEGs

    ##   gene   pvalue      lfc       padj direction
    ## 1 Fosb 1.478562 1.585131 0.03322296      DISS

Fos B and Jun overlap at 10 min

    tenmin <- rename(S2, c(`RNA fold change (10 min/control), log2` ="lfc", 
                       `p-value (10 min)` = "pvalue",
                       `Gene Symbol` = "gene"))
    tenmin <- wrangleCho(tenmin)

    Cho10min <- tenmin %>%
      filter(direction != "none")
    Cho10min

    ##          gene        lfc    log10p       pvalue        direction
    ## 1       Npas4  2.1166706 41.169699 6.765522e-42 fear-conditioned
    ## 2         Fos  0.8290431 15.870962 1.345977e-16 fear-conditioned
    ## 3        Junb  0.5444082 12.027934 9.377044e-13 fear-conditioned
    ## 4         Arc  0.4609445 10.312106 4.874092e-11 fear-conditioned
    ## 5       Ecel1 -0.4406157  5.414690 3.848661e-06          control
    ## 6        Ier2  0.6448592  5.103180 7.885335e-06 fear-conditioned
    ## 7        Pygl -0.4104363  2.390089 4.072964e-03          control
    ## 8  Rpl31-ps12  0.4906620  2.284587 5.192937e-03 fear-conditioned
    ## 9      Entpd2  0.4276484  1.440618 3.625620e-02 fear-conditioned
    ## 10      Cerkl -0.7314787  1.303379 4.973031e-02          control

    Cho10minCandidates <- Cho10min$gene
    Cho10minCandidates

    ##  [1] "Npas4"      "Fos"        "Junb"       "Arc"        "Ecel1"     
    ##  [6] "Ier2"       "Pygl"       "Rpl31-ps12" "Entpd2"     "Cerkl"

    Cho10minCandidatesGrepl <- c("Npas|Fos|Jun|Arc|Ecel|Ier|Pygl|Rpl31-ps|Entp|Cerk")

    # see if cho 4 hour DEGs are in my DEG list
    Cho_DEGs <- DEGs %>%
      dplyr::filter(gene %in% Cho10minCandidates) %>%
      arrange(gene)
    Cho_DEGs

    ## [1] gene      pvalue    lfc       padj      direction
    ## <0 rows> (or 0-length row.names)

    # see if related cho 4 hour DEGs are in my DEG list
    Cho_DEGs <- DEGs %>%
      dplyr::filter(grepl(Cho10minCandidatesGrepl, gene)) %>%
      arrange(gene)
    Cho_DEGs

    ##   gene   pvalue      lfc         padj direction
    ## 1 Fosb 1.478562 1.585131 0.0332229625      DISS
    ## 2  Jun 3.495075 1.596999 0.0003198344      DISS

Jun overlaps at 5 min

    fivemin <- rename(S2, c(`RNA fold change (5 min/control), log2` ="lfc", 
                       `p-value (5 min)` = "pvalue",
                       `Gene Symbol` = "gene"))
    fivemin <- wrangleCho(fivemin)

    Cho5min <- fivemin %>%
      filter(direction != "none")
    Cho5min

    ##       gene        lfc   log10p       pvalue        direction
    ## 1     Junb  0.4407774 2.038749 9.146427e-03 fear-conditioned
    ## 2 Slc26a10  0.7345409 1.680363 2.087552e-02 fear-conditioned
    ## 3     Epn3 -0.4790191 1.531283 2.942502e-02          control
    ## 4     Nme2  0.6842354 4.292430 5.100002e-05 fear-conditioned
    ## 5      Prl  0.9715696 2.623884 2.377477e-03 fear-conditioned
    ## 6    Serhl -0.6814371 2.528999 2.958016e-03          control

    Cho5minCandidates <- Cho5min$gene
    Cho5minCandidates

    ## [1] "Junb"     "Slc26a10" "Epn3"     "Nme2"     "Prl"      "Serhl"

    Cho5minCandidatesGrepl <- c( "Jun|Slc26a|Epn|Nme|Prl|Serhl")

    # see if cho 4 hour DEGs are in my DEG list
    Cho_DEGs <- DEGs %>%
      dplyr::filter(gene %in% Cho5minCandidates) %>%
      arrange(gene)
    Cho_DEGs

    ## [1] gene      pvalue    lfc       padj      direction
    ## <0 rows> (or 0-length row.names)

    # see if related cho 4 hour DEGs are in my DEG list
    Cho_DEGs <- DEGs %>%
      dplyr::filter(grepl(Cho5minCandidatesGrepl, gene)) %>%
      arrange(gene)
    Cho_DEGs

    ##   gene   pvalue      lfc         padj direction
    ## 1  Jun 3.495075 1.596999 0.0003198344      DISS

Cembrowksi comparsion…. incomplete
----------------------------------

    # table downloaded from biojupies http://amp.pharm.mssm.edu/biojupies/notebook/zGZJvFUQF
    cembrowski <- read.table("../data/BioJupiesCembrowskiSignature.txt", header = 1)
    head(cembrowski)

    ##   gene_symbol     logFC    AveExpr         t      P.Value    adj.P.Val
    ## 1       PROX1 13.876409  3.6431261 10.811670 2.600771e-09 7.296508e-07
    ## 2      IGFBP5 10.655410  4.0444132  6.024944 1.062618e-05 4.829865e-04
    ## 3       C1QL2 10.009691  1.0010774 14.536638 2.117321e-11 2.088064e-08
    ## 4         DSP  9.995924  0.8950099 10.403023 4.755595e-09 1.199737e-06
    ## 5      STXBP6  9.664500  3.8723598 11.969999 5.138066e-10 2.039186e-07
    ## 6      MFSD2A  9.315131 -0.6966044 13.533442 6.935708e-11 4.614397e-08
    ##           B
    ## 1 11.551588
    ## 2  3.492611
    ## 3 14.128875
    ## 4 10.111786
    ## 5 13.103359
    ## 6 12.302552

    cembrowski <- rename(cembrowski, c(`logFC` ="lfc", 
                       `P.Value` = "pvalue",
                       `gene_symbol` = "gene"))


    wranglecembroski <- function(df){
      data <- df  
      data$log10p <- -log10(data$pvalue) 
      data <- data %>% select(gene, lfc, log10p, pvalue) %>%
        mutate(direction = ifelse(data$lfc > 1 & data$log10p > 1.3, 
                            yes = "DG", 
                            no = ifelse(data$lfc < -1 & data$log10p > 1.3, 
                                        yes = "CA1", 
                                        no = "none")))
      data$direction <- as.factor(data$direction)
      data$direction <- factor(data$direction, c("CA1", "none", "DG"))
      return(data)
    }


    cembrowski <- wranglecembroski(cembrowski)
    head(cembrowski)

    ##     gene       lfc    log10p       pvalue direction
    ## 1  PROX1 13.876409  8.584898 2.600771e-09        DG
    ## 2 IGFBP5 10.655410  4.973623 1.062618e-05        DG
    ## 3  C1QL2 10.009691 10.674213 2.117321e-11        DG
    ## 4    DSP  9.995924  8.322795 4.755595e-09        DG
    ## 5 STXBP6  9.664500  9.289200 5.138066e-10        DG
    ## 6 MFSD2A  9.315131 10.158909 6.935708e-11        DG

    volcanoplot1b <- plotvolcano(cembrowski, fourhoursRNA$lfc, fourhoursRNA$log10p, 
                                plottitle = "Cembrowski") +
      scale_color_manual(values = c("none" = "grey",
                                    "DG" = "#018571",
                                    "CA1" = "#a6611a"))
    volcanoplot1b

![](../figures/04_candidategenes/cembrowksi-1.png)
