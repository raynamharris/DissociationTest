Import file with information about differentially expressed genes (aka
`dissociation`) and the list of genes found in the reference
transcriptome (aka `geneids`).

    dissociation <- read.csv("../results/01_dissociation_volcanoTreatment.csv", header = T, row.names = 1)
    geneids <- read.csv("../data/geneids.csv", header = T)

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
    sanesLichtman_reference

    ##   [1] Ache     Adcy1    Adcy10   Adcy2    Adcy3    Adcy4    Adcy5   
    ##   [8] Adcy6    Adcy7    Adcy8    Adcy9    Adra1a   Adra1b   Adra1d  
    ##  [15] Adra2a   Adra2b   Adra2c   Adrb1    Adrb2    Adrb3    Bdnf    
    ##  [22] Cacna1a  Cacna1b  Cacna1c  Cacna1d  Cacna1e  Cacna1f  Cacna1s 
    ##  [29] Calb1    Calb2    Calm1    Calm2    Calm3    Capn1    Capn10  
    ##  [36] Capn11   Capn12   Capn13   Capn15   Capn2    Capn3    Capn5   
    ##  [43] Capn6    Capn7    Capn8    Capn9    Cast     Ccr7     Cd47    
    ##  [50] Cdh1     Cdh2     Chrm1    Chrm2    Chrm3    Chrm4    Chrm5   
    ##  [57] Chrna1   Chrna3   Chrna7   Chrnb1   Chrnb2   Chrnb3   Cnga2   
    ##  [64] Cnr1     Creb1    Dlg4     Drd1     Edn1     Efna5    Egf     
    ##  [71] Egr1     Egr2     Epha5    Erbb4    Fgf2     Fyn      Gabbr1  
    ##  [78] Gabra1   Gabra2   Gabra3   Gabra5   Gabra6   Gabrb1   Gabrb2  
    ##  [85] Gabrb3   Gabrr1   Gap43    Gfap     Gria1    Gria2    Grin1   
    ##  [92] Grin2a   Grin2d   Grm1     Grm4     Grm5     Grm7     Gucy1a2 
    ##  [99] Gucy1a3  Gucy1b2  Gucy1b3  Gucy2c   Gucy2d   Gucy2e   Gucy2g  
    ## [106] Homer1   Homer2   Homer3   Htr1a    Htr1b    Htr1f    Htr2a   
    ## [113] Htr2b    Htr2c    Htr3a    Htr3b    Htr4     Htr5a    Htr5b   
    ## [120] Htr6     Htr7     Icam5    Il1b     Inhba    Itga1    Itga10  
    ## [127] Itga11   Itga2    Itga2b   Itga3    Itga4    Itga5    Itga6   
    ## [134] Itga7    Itga8    Itga9    Itgad    Itgae    Itgal    Itgam   
    ## [141] Itgav    Itgax    Itgb1    Itgb1bp1 Itgb2    Itgb2l   Itgb3   
    ## [148] Itgb3bp  Itgb4    Itgb5    Itgb6    Itgb7    Itgb8    Itgbl1  
    ## [155] L1cam    Mapk1    Mapk10   Mapk11   Mapk12   Mapk14   Mapk3   
    ## [162] Mapk4    Mapk6    Mapk7    Mapk8    Mapk9    Mas1     Napa    
    ## [169] Ncam1    Ncam2    Ngf      Nos1     Nos2     Nos3     Nrg1    
    ## [176] Nrg2     Nrg3     Nrgn     Ntrk2    Oprd1    Oprk1    Oprl1   
    ## [183] Oprm1    Parp1    Phpt1    Pla2g10  Pla2g12a Pla2g12b Pla2g15 
    ## [190] Pla2g16  Pla2g1b  Pla2g2a  Pla2g2c  Pla2g2d  Pla2g2e  Pla2g2f 
    ## [197] Pla2g3   Pla2g4a  Pla2g4b  Pla2g4e  Pla2g4f  Pla2g5   Pla2g6  
    ## [204] Pla2g7   Plat     Plcb1    Plcb2    Plcb3    Plcb4    Plcg1   
    ## [211] Plcg2    Plg      Pnoc     Ppp3ca   Ppp3cb   Ppp3cc   Prkcg   
    ## [218] Prkcz    Prkg1    Ptn      Rab3a    Rara     Rarb     S100b   
    ## [225] Serpine2 Snap25   Sptan1   Sptbn1   Src      Stmn4    Stx1b   
    ## [232] Syn1     Syp      Th       Thy1     Tnc      Ube3a    Vamp1   
    ## [239] Vamp2    Vamp3    Vamp4    Vamp5    Vamp8   
    ## 243 Levels: Ache Adcy1 Adcy10 Adcy2 Adcy3 Adcy4 Adcy5 Adcy6 Adcy7 ... Vamp8

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

    # percent DEGs in the sanes lichtman list
    round(6/243*100,2)

    ## [1] 2.47

    round(6/178*100,2)

    ## [1] 3.37

    # Grin2a NMDAR 2A
    # Il1b   Interleukin 1 Beta
    # Itga5  Integrin Subunit Alpha 5
    # Itgam  Integrin Subunit Alpha M
    # Itgb5  Integrin Subunit Beta 5
    # Mapk3  MAP Kinase 3
