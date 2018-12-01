Import file with information about differentially expressed genes.

    dissociation <- read.csv("../results/01_dissociation_volcanoTreatment.csv", header = T, row.names = 1)
    tail(dissociation)

    ##         gene     pvalue        lfc      padj color
    ## 16704   Zxdb 0.05155582 -0.3295400 0.8880638  none
    ## 16705   Zxdc 0.03591434  0.3481733 0.9206311  none
    ## 16706 Zyg11b 0.74313214 -0.9573600 0.1806624  none
    ## 16707    Zyx 0.47664005  1.2908846 0.3337029  none
    ## 16708  Zzef1 0.02879537  0.1606400 0.9358465  none
    ## 16709   Zzz3 0.11026785  0.4241942 0.7757685  none

Filter out non-significant genes, sort by p-value, rename column, round
to 2 decimal places.

    DEGs <- dissociation %>%
      filter(color != "none") %>%
      arrange((padj))
    DEGs$pvalue <- NULL # drop log pvalue columns
    names(DEGs)[4] <- "upregulated in"
    DEGs$lfc <- signif(DEGs$lfc, digits = 2)
    DEGs$padj <- signif(DEGs$padj, digits = 3)

    write.csv(DEGs, file = "../results/SuppTable1.csv", row.names = F)

Compare to Molecules implicated in hippocampal LTP
--------------------------------------------------

    # import list of all genes in mouse reference transcriptome used in this analysis
    geneids <- read.csv("../data/geneids.csv", header = T)

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
           "Src", "Fyn")
    sanesLichtman[order(sanesLichtman)] # print list alphabetically

    ##   [1] "Adra1a"   "Adra1b"   "Adra1d"   "Adra2a"   "Adra2b"   "Adra2c"  
    ##   [7] "Adrb1"    "Adrb2"    "Adrb3"    "Bdnf"     "Cacna1a"  "Cacna1b" 
    ##  [13] "Cacna1c"  "Cacna1d"  "Cacna1e"  "Cacna1f"  "Cacna1s"  "Calb1"   
    ##  [19] "Calb2"    "Calm1"    "Calm2"    "Calm3"    "Cd47"     "Cdh1"    
    ##  [25] "Cdh2"     "Chrm1"    "Chrm2"    "Chrm3"    "Chrm4"    "Chrm5"   
    ##  [31] "Chrna1"   "Chrna3"   "Chrna7"   "Chrnb1"   "Chrnb2"   "Chrnb3"  
    ##  [37] "Cnga2"    "Cnr1"     "Creb1"    "Dlg4"     "Drd1"     "Edn1"    
    ##  [43] "Efna5"    "Egf"      "Egr1"     "Egr2"     "Epha5"    "Erbb4"   
    ##  [49] "Fgf2"     "Fyn"      "Gabbr1"   "Gabra1"   "Gabra2"   "Gabra3"  
    ##  [55] "Gabra5"   "Gabra6"   "Gabrb1"   "Gabrb2"   "Gabrb3"   "Gabrr1"  
    ##  [61] "Gap43"    "Gria1"    "Gria2"    "Grin1"    "Grin2a"   "Grin2d"  
    ##  [67] "Grm1"     "Grm4"     "Grm5"     "Grm7"     "Htr1a"    "Htr1b"   
    ##  [73] "Htr1f"    "Htr2a"    "Htr2b"    "Htr2c"    "Htr3a"    "Htr3b"   
    ##  [79] "Htr4"     "Htr5a"    "Htr5b"    "Htr6"     "Htr7"     "Icam5"   
    ##  [85] "Il1b"     "Inhba"    "Itga1"    "Itga10"   "Itga11"   "Itga2"   
    ##  [91] "Itga2b"   "Itga3"    "Itga4"    "Itga5"    "Itga6"    "Itga7"   
    ##  [97] "Itga8"    "Itga9"    "Itgad"    "Itgae"    "Itgal"    "Itgam"   
    ## [103] "Itgav"    "Itgax"    "Itgb1"    "Itgb1bp1" "Itgb2"    "Itgb2l"  
    ## [109] "Itgb3"    "Itgb3bp"  "Itgb4"    "Itgb5"    "Itgb6"    "Itgb7"   
    ## [115] "Itgb8"    "Itgbl1"   "L1cam"    "Mapk1"    "Mapk10"   "Mapk11"  
    ## [121] "Mapk12"   "Mapk14"   "Mapk3"    "Mapk4"    "Mapk6"    "Mapk7"   
    ## [127] "Mapk8"    "Mapk9"    "Napa"     "Ncam1"    "Ncam2"    "Ngf"     
    ## [133] "Nos1"     "Nos2"     "Nos3"     "Nrg1"     "Nrg2"     "Nrg3"    
    ## [139] "Nrgn"     "Ntrk2"    "Oprd1"    "Oprk1"    "Oprl1"    "Oprm1"   
    ## [145] "Pnoc"     "Ptn"      "Rab3a"    "Rara"     "Rarb"     "S100b"   
    ## [151] "Snap25"   "Src"      "Stx1b"    "Syn1"     "Syp"      "Th"      
    ## [157] "Thy1"     "Tnc"      "Vamp1"    "Vamp2"    "Vamp3"    "Vamp4"   
    ## [163] "Vamp5"    "Vamp8"

    # confirm that all all Sanes and Lichtman genes are in the reference transcriptome
    sanesLichtman_reference <- geneids %>%
      dplyr::filter(gene %in% sanesLichtman) %>%
      dplyr::select(gene) %>%
      distinct() %>%
      arrange(gene) %>%
      droplevels()
    sanesLichtman_reference <- sanesLichtman_reference[,c(1)]
    sanesLichtman_reference

    ##   [1] Adra1a   Adra1b   Adra1d   Adra2a   Adra2b   Adra2c   Adrb1   
    ##   [8] Adrb2    Adrb3    Bdnf     Cacna1a  Cacna1b  Cacna1c  Cacna1d 
    ##  [15] Cacna1e  Cacna1f  Cacna1s  Calb1    Calb2    Calm1    Calm2   
    ##  [22] Calm3    Cd47     Cdh1     Cdh2     Chrm1    Chrm2    Chrm3   
    ##  [29] Chrm4    Chrm5    Chrna1   Chrna3   Chrna7   Chrnb1   Chrnb2  
    ##  [36] Chrnb3   Cnga2    Cnr1     Creb1    Dlg4     Drd1     Edn1    
    ##  [43] Efna5    Egf      Egr1     Egr2     Epha5    Erbb4    Fgf2    
    ##  [50] Fyn      Gabbr1   Gabra1   Gabra2   Gabra3   Gabra5   Gabra6  
    ##  [57] Gabrb1   Gabrb2   Gabrb3   Gabrr1   Gap43    Gria1    Gria2   
    ##  [64] Grin1    Grin2a   Grin2d   Grm1     Grm4     Grm5     Grm7    
    ##  [71] Htr1a    Htr1b    Htr1f    Htr2a    Htr2b    Htr2c    Htr3a   
    ##  [78] Htr3b    Htr4     Htr5a    Htr5b    Htr6     Htr7     Icam5   
    ##  [85] Il1b     Inhba    Itga1    Itga10   Itga11   Itga2    Itga2b  
    ##  [92] Itga3    Itga4    Itga5    Itga6    Itga7    Itga8    Itga9   
    ##  [99] Itgad    Itgae    Itgal    Itgam    Itgav    Itgax    Itgb1   
    ## [106] Itgb1bp1 Itgb2    Itgb2l   Itgb3    Itgb3bp  Itgb4    Itgb5   
    ## [113] Itgb6    Itgb7    Itgb8    Itgbl1   L1cam    Mapk1    Mapk10  
    ## [120] Mapk11   Mapk12   Mapk14   Mapk3    Mapk4    Mapk6    Mapk7   
    ## [127] Mapk8    Mapk9    Napa     Ncam1    Ncam2    Ngf      Nos1    
    ## [134] Nos2     Nos3     Nrg1     Nrg2     Nrg3     Nrgn     Ntrk2   
    ## [141] Oprd1    Oprk1    Oprl1    Oprm1    Pnoc     Ptn      Rab3a   
    ## [148] Rara     Rarb     S100b    Snap25   Src      Stx1b    Syn1    
    ## [155] Syp      Th       Thy1     Tnc      Vamp1    Vamp2    Vamp3   
    ## [162] Vamp4    Vamp5    Vamp8   
    ## 164 Levels: Adra1a Adra1b Adra1d Adra2a Adra2b Adra2c Adrb1 Adrb2 ... Vamp8

    # identify which of the Sanes and Lichtman genes are present in my samples
    sanesLichtman_present <- dissociation %>%
      dplyr::filter(gene %in% sanesLichtman) %>%
      droplevels()
    sanesLichtman_present <- sanesLichtman_present[,c(1)]
    sanesLichtman_present

    ##   [1] Adra1a   Adra1d   Adra2a   Adra2b   Adra2c   Adrb1    Bdnf    
    ##   [8] Cacna1a  Cacna1b  Cacna1c  Cacna1d  Cacna1e  Calb1    Calb2   
    ##  [15] Calm1    Calm2    Calm3    Cd47     Cdh2     Chrm1    Chrm2   
    ##  [22] Chrm3    Chrm4    Chrna7   Chrnb2   Cnr1     Creb1    Dlg4    
    ##  [29] Egr1     Egr2     Epha5    Erbb4    Fgf2     Fyn      Gabbr1  
    ##  [36] Gabra1   Gabra2   Gabra3   Gabra5   Gabrb1   Gabrb2   Gabrb3  
    ##  [43] Gap43    Gria1    Gria2    Grin1    Grin2a   Grin2d   Grm1    
    ##  [50] Grm4     Grm5     Grm7     Htr1a    Htr1b    Htr2a    Htr2c   
    ##  [57] Htr3a    Htr4     Htr5a    Htr5b    Htr6     Htr7     Icam5   
    ##  [64] Il1b     Itga1    Itga10   Itga11   Itga2b   Itga3    Itga4   
    ##  [71] Itga5    Itga6    Itga7    Itga8    Itga9    Itgam    Itgav   
    ##  [78] Itgb1    Itgb1bp1 Itgb2    Itgb4    Itgb5    Itgb6    Itgb8   
    ##  [85] Itgbl1   L1cam    Mapk1    Mapk10   Mapk11   Mapk12   Mapk14  
    ##  [92] Mapk3    Mapk4    Mapk6    Mapk7    Mapk8    Mapk9    Napa    
    ##  [99] Ncam1    Ncam2    Nos1     Nos3     Nrg1     Nrg2     Nrg3    
    ## [106] Nrgn     Ntrk2    Oprd1    Oprl1    Pnoc     Ptn      Rab3a   
    ## [113] Rara     S100b    Snap25   Src      Stx1b    Syn1     Syp     
    ## [120] Th       Thy1     Vamp1    Vamp2    Vamp3    Vamp4    Vamp8   
    ## 126 Levels: Adra1a Adra1d Adra2a Adra2b Adra2c Adrb1 Bdnf ... Vamp8

    # identify whichof the Sanes and Lichtman genes are differentially expressed in this analysis
    sanesLichtman_DEGs <- DEGs %>%
      dplyr::filter(gene %in% sanesLichtman) %>%
      arrange(gene)
    sanesLichtman_DEGs

    ##     gene  lfc   padj upregulated in
    ## 1 Grin2a -1.7 0.0269           HOMO
    ## 2   Il1b  2.4 0.0300           DISS
    ## 3  Itga5  3.1 0.0160           DISS
    ## 4  Itgam  1.7 0.0496           DISS
    ## 5  Itgb5  2.0 0.0201           DISS
    ## 6  Mapk3  1.6 0.0098           DISS

    candidates <- dissociation %>%
      dplyr::filter(grepl('Mapk', gene)) %>%
      dplyr::select(gene) %>%
      distinct() %>%
      arrange(gene) %>%
      droplevels()
    candidates <- candidates[,c(1)]
    candidates

    ##  [1] Mapk1     Mapk10    Mapk11    Mapk12    Mapk14    Mapk1ip1  Mapk1ip1l
    ##  [8] Mapk3     Mapk4     Mapk6     Mapk7     Mapk8     Mapk8ip1  Mapk8ip2 
    ## [15] Mapk8ip3  Mapk9     Mapkap1   Mapkapk2  Mapkapk3  Mapkapk5  Mapkbp1  
    ## 21 Levels: Mapk1 Mapk10 Mapk11 Mapk12 Mapk14 Mapk1ip1 Mapk1ip1l ... Mapkbp1
