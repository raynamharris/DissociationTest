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

Cho et al. data at 4 hours

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

Cho et al. data at 30 min

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

    ### Plotting their data with my differential exprssion

    overlap30min <- overlap(thirtyminRNA) # no DEGs
    overlap4h <- overlap(fourhoursRNA)
    summary(overlap4h$color)

    ##           absent             none          control fear-conditioned 
    ##             1566            11403                9                0 
    ##             HOMO             DISS 
    ##               11              138

    fivemin <- rename(S2, c(`RNA fold change (5 min/control), log2` ="lfc", 
                       `p-value (5 min)` = "pvalue",
                       `Gene Symbol` = "gene"))
    fivemin <- wrangleCho(fivemin)

    tenmin <- rename(S2, c(`RNA fold change (10 min/control), log2` ="lfc", 
                       `p-value (10 min)` = "pvalue",
                       `Gene Symbol` = "gene"))
    tenmin <- wrangleCho(tenmin)

    overlaptenmin <- overlap(tenmin) # no DEGs
    overlapfivemin <- overlap(fivemin)
    summary(tenmin$color)

    ## Length  Class   Mode 
    ##      0   NULL   NULL

    summary(fivemin$color)

    ## Length  Class   Mode 
    ##      0   NULL   NULL

    overlap4h %>%
      filter(Cho == "control",
             Harris =="HOMO")  %>%
      select(gene, lfc.x, log10p) %>%
      arrange(gene)

    ## [1] gene   lfc.x  log10p
    ## <0 rows> (or 0-length row.names)

    overlap4h %>%
      filter(Cho == "control",
             Harris =="DISS")  %>%
      select(gene, lfc.x, log10p) %>%
      arrange(gene)

    ## [1] gene   lfc.x  log10p
    ## <0 rows> (or 0-length row.names)

    overlap4h %>%
      filter(Cho == "control")  %>%
      select(gene, lfc.x, log10p) %>%
      arrange(gene)

    ##     gene     lfc.x    log10p
    ## 1   Aqp1 -1.382172 12.836310
    ## 2 Calml4 -1.087474  5.187103
    ## 3  Cldn2 -1.003313  8.527052
    ## 4  Clic6 -1.017852  7.843026
    ## 5  Folr1 -1.102123 12.302898
    ## 6  Kcne2 -1.169621 12.652911
    ## 7 Slc4a5 -1.235330  6.039947
    ## 8 Tmem72 -1.654890  7.715675
    ## 9    Ttr -1.356586  5.464636

    suzyvolcano1 <- suzyvolcano(overlap30min, overlap30min$lfc.x,overlap30min$log10p, 
                                plottitle = NULL)
    suzyvolcano1

![](../figures/04_candidategenes/overlap-1.png)

    suzyvolcano1 <- suzyvolcano1 + 
        annotate("text", label = "338", x = -0.9, y = 2, size = 3, color = "black") + 
        annotate("text", label = "261", x = 0.9, y = 2, size = 3, color = "black")


    suzyvolcano2 <- suzyvolcano(overlap4h, overlap4h$lfc.x, overlap4h$log10p, 
                                plottitle = NULL)

    suzyvolcano2

![](../figures/04_candidategenes/overlap-2.png)

    suzyvolcano2 <- suzyvolcano2 + 
        annotate("text", label = "435", x = -0.9, y = 2, size = 3, color = "black") + 
        annotate("text", label = "593", x = 0.9, y = 2, size = 3, color = "black")


    # side by side plots 
    figure4 <- plot_grid(suzyvolcano1, suzyvolcano2, nrow = 1, labels = c('A', 'B'))

    figure4

![](../figures/04_candidategenes/overlap-3.png)

    ggsave(
      "../figures/figure4.png",
      figure4,
      width = 6,
      height = 4,
      dpi = 1200
    )

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
