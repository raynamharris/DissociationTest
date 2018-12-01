Import file with information about differentially expressed genes.

    dissociation <- read.csv("../results/01_dissociation_volcanoTreatment.csv", header = T, row.names = 1)
    head(dissociation)

    ##            gene     pvalue        lfc      padj color
    ## 1 0610007P14Rik 0.06147827 -0.5281158 0.8680040  none
    ## 2 0610009B22Rik 0.01810658  0.2891023 0.9591652  none
    ## 3 0610009L18Rik 0.01886522 -0.3916516 0.9574912  none
    ## 4 0610009O20Rik 0.03752116  0.2571160 0.9172312  none
    ## 5 0610010F05Rik 0.03988827 -0.2712169 0.9122455  none
    ## 6 0610010K14Rik 0.02680878  0.2984152 0.9401372  none

    str(dissociation)

    ## 'data.frame':    12157 obs. of  5 variables:
    ##  $ gene  : Factor w/ 12157 levels "0610007P14Rik",..: 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ pvalue: num  0.0615 0.0181 0.0189 0.0375 0.0399 ...
    ##  $ lfc   : num  -0.528 0.289 -0.392 0.257 -0.271 ...
    ##  $ padj  : num  0.868 0.959 0.957 0.917 0.912 ...
    ##  $ color : Factor w/ 3 levels "DISS","HOMO",..: 3 3 3 3 3 3 3 3 3 3 ...

Filter out non-significant genes, sort by p-value, rename column, round
to 2 decimal places.

    DEGes <- dissociation %>%
      filter(color != "none") %>%
      arrange((padj))
    DEGes$pvalue <- NULL # drop log pvalue columns
    names(DEGes)[4] <- "upregulated in"
    DEGes$lfc <- signif(DEGes$lfc, digits = 2)
    DEGes$padj <- signif(DEGes$padj, digits = 3)

    write.csv(DEGes, file = "../results/SuppTable1.csv", row.names = F)

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
                       "Egr1", "Egr2")
    sanesLichtman[order(sanesLichtman)] # print list alphabetically

    ##   [1] "Adra1a"  "Adra1b"  "Adra1d"  "Adra2a"  "Adra2b"  "Adra2c"  "Adrb1"  
    ##   [8] "Adrb2"   "Adrb3"   "Bdnf"    "Cacna1a" "Cacna1b" "Cacna1c" "Cacna1d"
    ##  [15] "Cacna1e" "Cacna1f" "Cacna1s" "Calb1"   "Calb2"   "Calm1"   "Calm2"  
    ##  [22] "Calm3"   "Chrm1"   "Chrm2"   "Chrm3"   "Chrm4"   "Chrm5"   "Chrna1" 
    ##  [29] "Chrna3"  "Chrna7"  "Chrnb1"  "Chrnb2"  "Chrnb3"  "Cnga2"   "Cnr1"   
    ##  [36] "Creb1"   "Dlg4"    "Drd1"    "Edn1"    "Egf"     "Egr1"    "Egr2"   
    ##  [43] "Erbb4"   "Fgf2"    "Gabbr1"  "Gabra1"  "Gabra2"  "Gabra3"  "Gabra5" 
    ##  [50] "Gabra6"  "Gabrb1"  "Gabrb2"  "Gabrb3"  "Gabrr1"  "Gap43"   "Gria1"  
    ##  [57] "Gria2"   "Grin1"   "Grin2a"  "Grin2d"  "Grm1"    "Grm4"    "Grm5"   
    ##  [64] "Grm7"    "Htr1a"   "Htr1b"   "Htr1f"   "Htr2a"   "Htr2b"   "Htr2c"  
    ##  [71] "Htr3a"   "Htr3b"   "Htr4"    "Htr5a"   "Htr5b"   "Htr6"    "Htr7"   
    ##  [78] "Il1b"    "Inhba"   "Napa"    "Ngf"     "Nos1"    "Nos2"    "Nos3"   
    ##  [85] "Nrg1"    "Nrg2"    "Nrg3"    "Nrgn"    "Ntrk2"   "Oprd1"   "Oprk1"  
    ##  [92] "Oprl1"   "Oprm1"   "Pnoc"    "Rab3a"   "Rara"    "Rarb"    "S100b"  
    ##  [99] "Snap25"  "Stx1b"   "Syn1"    "Syp"     "Th"      "Vamp1"   "Vamp2"  
    ## [106] "Vamp3"   "Vamp4"   "Vamp5"   "Vamp8"

    # confirm that all all Sanes and Lichtman genes are in the reference transcriptome
    sanesLichtman_reference <- geneids %>%
      dplyr::filter(gene %in% sanesLichtman) %>%
      dplyr::select(gene) %>%
      distinct() %>%
      arrange(gene) %>%
      droplevels()
    sanesLichtman_reference <- sanesLichtman_reference[,c(1)]
    sanesLichtman_reference

    ##   [1] Adra1a  Adra1b  Adra1d  Adra2a  Adra2b  Adra2c  Adrb1   Adrb2  
    ##   [9] Adrb3   Bdnf    Cacna1a Cacna1b Cacna1c Cacna1d Cacna1e Cacna1f
    ##  [17] Cacna1s Calb1   Calb2   Calm1   Calm2   Calm3   Chrm1   Chrm2  
    ##  [25] Chrm3   Chrm4   Chrm5   Chrna1  Chrna3  Chrna7  Chrnb1  Chrnb2 
    ##  [33] Chrnb3  Cnga2   Cnr1    Creb1   Dlg4    Drd1    Edn1    Egf    
    ##  [41] Egr1    Egr2    Erbb4   Fgf2    Gabbr1  Gabra1  Gabra2  Gabra3 
    ##  [49] Gabra5  Gabra6  Gabrb1  Gabrb2  Gabrb3  Gabrr1  Gap43   Gria1  
    ##  [57] Gria2   Grin1   Grin2a  Grin2d  Grm1    Grm4    Grm5    Grm7   
    ##  [65] Htr1a   Htr1b   Htr1f   Htr2a   Htr2b   Htr2c   Htr3a   Htr3b  
    ##  [73] Htr4    Htr5a   Htr5b   Htr6    Htr7    Il1b    Inhba   Napa   
    ##  [81] Ngf     Nos1    Nos2    Nos3    Nrg1    Nrg2    Nrg3    Nrgn   
    ##  [89] Ntrk2   Oprd1   Oprk1   Oprl1   Oprm1   Pnoc    Rab3a   Rara   
    ##  [97] Rarb    S100b   Snap25  Stx1b   Syn1    Syp     Th      Vamp1  
    ## [105] Vamp2   Vamp3   Vamp4   Vamp5   Vamp8  
    ## 109 Levels: Adra1a Adra1b Adra1d Adra2a Adra2b Adra2c Adrb1 Adrb2 ... Vamp8

    # identify which of the Sanes and Lichtman genes are present in my samples
    sanesLichtman_present <- dissociation %>%
      dplyr::filter(gene %in% sanesLichtman) %>%
      droplevels()
    sanesLichtman_present <- sanesLichtman_present[,c(1)]
    sanesLichtman_present

    ##  [1] Adra1a  Adra1d  Adra2a  Adra2b  Adra2c  Adrb1   Bdnf    Cacna1a
    ##  [9] Cacna1b Cacna1c Cacna1d Cacna1e Calb1   Calb2   Calm1   Calm2  
    ## [17] Calm3   Chrm1   Chrm2   Chrm3   Chrm4   Chrna7  Chrnb2  Cnr1   
    ## [25] Creb1   Dlg4    Egr1    Egr2    Erbb4   Fgf2    Gabbr1  Gabra1 
    ## [33] Gabra2  Gabra3  Gabra5  Gabrb1  Gabrb2  Gabrb3  Gap43   Gria1  
    ## [41] Gria2   Grin1   Grin2a  Grin2d  Grm1    Grm4    Grm5    Grm7   
    ## [49] Htr1a   Htr1b   Htr2a   Htr2c   Htr3a   Htr4    Htr5a   Htr5b  
    ## [57] Htr6    Htr7    Il1b    Napa    Nos1    Nos3    Nrg1    Nrg2   
    ## [65] Nrg3    Nrgn    Ntrk2   Oprd1   Oprl1   Pnoc    Rab3a   Rara   
    ## [73] S100b   Snap25  Stx1b   Syn1    Syp     Th      Vamp1   Vamp2  
    ## [81] Vamp3   Vamp4   Vamp8  
    ## 83 Levels: Adra1a Adra1d Adra2a Adra2b Adra2c Adrb1 Bdnf ... Vamp8

    # identify whichof the Sanes and Lichtman genes are differentially expressed in this analysis
    sanesLichtman_DEGs <- DEGes %>%
      dplyr::filter(gene %in% sanesLichtman)
    sanesLichtman_DEGs

    ##     gene  lfc   padj upregulated in
    ## 1 Grin2a -1.7 0.0269           HOMO
    ## 2   Il1b  2.4 0.0300           DISS
