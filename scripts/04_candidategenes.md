    library(readxl)
    library(ggplot2)
    library(cowplot)
    library(dplyr)
    library(plyr)
    # note: won't run if plyr loaded before dplyr

    # set output file for figures and general chunk settings
    knitr::opts_chunk$set(fig.path = '../figures/04_candidategenes/', echo = T, message = F, warning=F)

Comparison with other candidate gene lists
==========================================

1.  Sanes and Licthman 1999
2.  Cho et al 2015
3.  Cahoy et al 2008

Sanes and licthman
------------------

Import file with information about differentially expressed genes (aka
`dissociation`) and the list of genes found in the reference
transcriptome (aka `geneids`).

    dissociation <- read.csv("../results/volcanoTreatment.csv", header = T, row.names = 1)
    dissociation$lfc <- round(dissociation$lfc,2)
    dissociation$padj <- formatC(dissociation$padj, format = "e", digits = 2)

    geneids <- read.csv("../data/geneids.csv", header = T)

    # Filter out non-significant genes, sort by p-value, rename column, round to 2 decimal places. 

    DEGs <- dissociation %>%
      dplyr::filter(direction != "neither") %>%
      arrange((padj))

    # Compare to Molecules implicated in hippocampal LTP from sanes

    supptable2 <- read.csv("../data/SanesLichtman.csv", check.names = F)
    head(supptable2, 10)

    ##                   Sanes & Lichtman Molecules
    ## 1                        GLUTAMATE RECEPTORS
    ## 2                               GluR1; GluR2
    ## 3             mGluR1; mGluR4; mGluR5; mGluR7
    ## 4             NMDA NR2A; NMDA NR2D; NMDA NR1
    ## 5                    OTHER NEUROTRANSMITTERS
    ## 6  norepinephrine and b-adrenergic receptors
    ## 7       adenosine and adenosine 2A receptors
    ## 8         dopamine and D1 dopamine receptors
    ## 9              mu and delta opioid receptors
    ## 10                   acetylcholine receptors
    ##                                Related Transcripts
    ## 1                                                 
    ## 2                                    Gria1; Gria2 
    ## 3                          Grm1; Grm4; Grm5; Grm7 
    ## 4                           Grin1; Grin2a; Grin2d 
    ## 5                                                 
    ## 6                             Adrb1; Adrb2; Adrb3 
    ## 7   Adra1a; Adra1b; Adra1d; Adra2a; Adra2b; Adra2c
    ## 8                                        Th; Drd1 
    ## 9                                     Oprm1; Oprd1
    ## 10  Chrna1; Chrna7; Chrna3; Chrnb1; Chrnb2; Chrnb3

    # create list of candidate genes Sanes and Lichtman 1999 

    sanesLichtman <- c("Gria1", "Gria2", 
           "Grm1", "Grm4", "Grm5", "Grm7",
           "Grin1", "Grin2a", "Grin2d", 
           "Th", "Drd1",
           "Adrb1", "Adrb2", "Adrb3",
           "Adra2a", "Adra2b", "Adra2c",
           "Oprm1", "Oprd1",
           "Chrm1", "Chrm2", "Chrm3", "Chrm4", "Chrm5",
           "Chrna1", "Chrna7", "Chrna3", 
           "Chrnb1", "Chrnb2", "Chrnb3",
           "Gabra1", "Gabra2",  "Gabra3", "Gabra5", "Gabra6",
           "Gabrb1", "Gabrb2", "Gabrb3",  
           "Gabrr1",  "Gabbr1",
           "Cnr1", "Cnr2",
           "Pnoc", "Oprl1",
           "Htr1a", "Htr1b", "Htr1f",
           "Htr2a", "Htr2c", "Htr2b",
           "Htr3a", "Htr3b", "Htr5a", "Htr5b",
           "Htr7", "Htr6", "Htr4", 
           "Edn1", "Egf", "Fgf2",
           "Nrg1", "Nrg2", "Nrg3",
           "Erbb4", "Ngf", "Bdnf", "Ntrk2",
           "Nos1", "Nos3",
           "Il1b",
           "Inhba", "Calm1", "Calm2", "Calm3",
           "Nrgn", "Calb1", "Calb2", "Gap43", "S100b",
           "Cacna1c", "Cacna1d", "Cacna1s", "Cacna1f",
           "Cacna1b", "Cacna1a", "Cacna1e",
           "Cnga2", "Syp", "Napa",
           "Vamp1", "Vamp2", "Vamp3", "Vamp4", "Vamp5", "Vamp8",
           "Rab3a", "Stx1b", "Syn1", "Snap25",
           "Dlg4", "Rarb", "Creb1",
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
           "Itpka", "Itpkb", "Itpkc", 
           "Mapk1", "Mapk10", "Mapk11", "Mapk12", "Mapk14", 
           "Mapk3", "Mapk4", "Mapk6", "Mapk7", "Mapk8", "Mapk9",
           "Src", "Fyn", 
           "Prkacb", "Prkar1b",
           "Prkcg", "Prkg1", "Prkcz", 
           "Camk1",  "Camk2",  "Camk4",
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
           "Adcy1", 
           "Gucy1a2", "Gucy1a3", "Gucy1b2", "Gucy1b3",
           "Gucy2c", "Gucy2d", "Gucy2e", "Gucy2g",
           "Sptan1", "Sptbn1", "Gfap", "Stmn4",
           "Ccr7", "Mas1",
           "Homer1", "Homer2", "Homer3" )

    #sanesLichtman[order(sanesLichtman)] # print list alphabetically

    # confirm that all all Sanes and Lichtman genes are in the reference transcriptome
    sanesLichtman_reference <- geneids %>%
      dplyr::filter(gene %in% sanesLichtman) %>%
      dplyr::select(gene) %>%
      distinct() %>%
      arrange(gene) %>%
      droplevels()
    sanesLichtman_reference <- sanesLichtman_reference[,c(1)]
    str(sanesLichtman_reference)

    ##  Factor w/ 236 levels "Ache","Adcy1",..: 1 2 3 4 5 6 7 8 9 10 ...

    # identify which of the Sanes and Lichtman genes are present in my samples
    sanesLichtman_present <- dissociation %>%
      dplyr::filter(gene %in% sanesLichtman) %>%
      droplevels()
    sanesLichtman_present <- sanesLichtman_present[,c(1)]
    str(sanesLichtman_present)

    ##  Factor w/ 175 levels "Ache","Adcy1",..: 1 2 3 4 5 6 7 8 9 10 ...

    # identify whichof the Sanes and Lichtman genes are differentially expressed in this analysis
    sanesLichtman_DEGs <- DEGs %>%
      dplyr::filter(gene %in% sanesLichtman) %>%
        dplyr::filter(direction != "neither") %>%
      arrange(gene)
    sanesLichtman_DEGs

    ##       gene   pvalue   lfc     padj direction
    ## 1  Cacna1e 1.094199 -1.28 8.05e-02      HOMO
    ## 2   Gabrb1 1.987167 -1.07 1.03e-02      HOMO
    ## 3   Grin2a 1.570140 -1.66 2.69e-02      HOMO
    ## 4     Il1b 1.523350  2.41 3.00e-02      DISS
    ## 5    Itga5 1.795090  3.05 1.60e-02      DISS
    ## 6    Itgam 1.304370  1.75 4.96e-02      DISS
    ## 7    Itgb4 1.200986  2.93 6.30e-02      DISS
    ## 8    Itgb5 1.695894  1.98 2.01e-02      DISS
    ## 9    Itpkb 1.033912  1.53 9.25e-02      DISS
    ## 10   Mapk3 2.008779  1.61 9.80e-03      DISS

    # what are these genes?

    # Cacna1e Calcium Voltage-Gated Channel Subunit Alpha1 
    # Gabrb1  Gamma-Aminobutyric Acid Type A Receptor Beta1 Subunit
    # Grin2a  NMDAR 2A
    # Il1b    Interleukin 1 Beta
    # Itga5   Integrin Subunit Alpha 5
    # Itgam   Integrin Subunit Alpha M
    # Itgb4   Integrin Subunit Beta 4
    # Itgb5   Integrin Subunit Beta 5
    # Itpkb   Inositol 1,4,5-Trisphosphate 3-Kinase B 
    # Mapk3   MAP Kinase 3

    # percent DEGs in the sanes lichtman list
    round(11/236*100,2)

    ## [1] 4.66

    # percent DEGs in the sanes lichtman list AND present
    round(11/175*100,2)

    ## [1] 6.29

Cho et al 2015 anlaysis
-----------------------

Cho et al 2015 used RNA sequencing to quantify transcript levels in the
mouse hippocampus after contextual fear conditioning. The Cho dataset
provides a snapshot of gene expression changes associated with
hippocampal learning and memory 30 min and 4 h after an experiment. The
Cho data are available at
<a href="http://science.sciencemag.org/content/suppl/2015/09/30/350.6256.82.DC1" class="uri">http://science.sciencemag.org/content/suppl/2015/09/30/350.6256.82.DC1</a>
and
<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72064" class="uri">https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72064</a>.
The file `../data/aac7368-Cho.SM.Table.S2.xls` of differentially
expressed genes was used as a representative dataset for learning and
memory associated gene expression patterns.

In this analysis, I compared the Cho et al differentially expressed
genes (DEGs) to my experimental results (referred to as the Harris data)
to identify the interaction between genes that are differentially
expressed following fear condition and following a chemical
manipulation.

This analysis prints the list of genes that are differentially expressed
in both experiments. The images show that only a few (red dots) of genes
that respond to chemical dissociation are up-regulated or down-regulated
following fear-conditioning.

The Cho data is not tidy It is very wide, with fold change and pvalue
scores for gene expression at four different timepoints. I want to
subset this one dataframe into four smaller data frames, one for each
timepoint, each with column headings: gene, lfc, and pvalue.

    # read supplemental data from cho et al
    S2 <- as.data.frame(readxl::read_excel("../data/aac7368-Cho.SM.Table.S2.xls", skip = 1 ))
    S2 <- rename(S2, c(`Gene Symbol` = "gene")) # rename gene column
    names(S2)

    ##  [1] "Refseq accession"                      
    ##  [2] "gene"                                  
    ##  [3] "RPF fold change (5 min/control), log2" 
    ##  [4] "RPF fold change (10 min/control), log2"
    ##  [5] "RPF fold change (30 min/control), log2"
    ##  [6] "RPF fold change (4 h/control), log2"   
    ##  [7] "RNA fold change (5 min/control), log2" 
    ##  [8] "RNA fold change (10 min/control), log2"
    ##  [9] "RNA fold change (30 min/control), log2"
    ## [10] "RNA fold change (4 h/control), log2"   
    ## [11] "p-value (5 min)"                       
    ## [12] "p-value (10 min)"                      
    ## [13] "p-value (30 min)"                      
    ## [14] "p-value (4 h)"                         
    ## [15] "FDR (5 min)"                           
    ## [16] "FDR (10 min)"                          
    ## [17] "FDR (30 min)"                          
    ## [18] "FDR (4 h)"                             
    ## [19] "Description"

So I wrote a few functions to:

1.  `subset_df`: subset the columns of interest and reanme them
2.  `determineChoDEGs`: catagorize the genes as differentially expressed
    or not
3.  `createDEGlist`: create a list of just DEG gene names
4.  `comparetoDISS`: crossreference the list of DEG names with the the
    dissociated DEGs

<!-- -->

    # subset data
    subset_df <- function(df, lfc_col, pval_col){
      df <- df %>% select(gene, lfc_col, pval_col)
      colnames(df) <- c("gene","lfc","pvalue")
      return(df)
    }

    # identify DEGs in cho dataset. Call anything with pvalue <0.05 and lfc greater than 0.3
    determineChoDEGs <- function(df){
      data <- df  
      data$log10p <- (-log10(data$pvalue)) 
      data <- data %>% 
        dplyr::mutate(direction = ifelse(data$lfc > 0.3 & data$pvalue < 0.05, 
                            yes = "fear-conditioned", 
                            no = ifelse(data$lfc < -0.3 & data$pvalue < 0.05, 
                                        yes = "control", 
                                        no = "neither")))
      data$direction <- as.factor(data$direction)
      data$direction <- factor(data$direction, c("control", "neither", "fear-conditioned"))
      return(data)
    }

    # extract just the list of Cho DEGs
    createDEGlist <- function(df){
      DEGs_df <- df %>%
        dplyr::filter(direction != "neither") %>%
        dplyr::arrange(gene)
      DEGs_list <- DEGs_df$gene
      return(DEGs_list)
    }

    # see which Cho DEGs are also differentailly expressed in the dissociation treatment
    comparetoDISS <- function(listofDEGs){
      comparison <- DEGs %>%
        dplyr::filter(gene %in% listofDEGs) %>%
        dplyr::filter(direction != "neither") %>%
        dplyr::arrange(gene)
      return(comparison)
    }

### four hour DEGs

In the Cho et al. data, there are 9 differentially expressed genes after
4 hours with LFC &gt; 1. But there are 35 with lfc &gt; 0.25. I will go
with those since that is about the cuttoff used in the Cho paper.

    fourhours_df <- subset_df(S2, "RNA fold change (4 h/control), log2", "p-value (4 h)")
    fourhours_df <- determineChoDEGs(fourhours_df)
    fourhourDEGs <- createDEGlist(fourhours_df)
    fourhourcomparison <- comparetoDISS(fourhourDEGs)
    fourhourcomparison

    ##    gene   pvalue  lfc     padj direction
    ## 1 Enpp2 1.310489 1.77 4.89e-02      DISS
    ## 2   Fn1 1.570140 1.77 2.69e-02      DISS

    # what genes are they?

    # Enpp2 - Ectonucleotide Pyrophosphatase/Phosphodiesterase 2
    # Fn1 - Fibronectin 1 - cell adhesion

### thirty minute DEGs

    # create 30 min df and rename column file and column headings for easy analysis
    thirtymin <- subset_df(S2, "RNA fold change (30 min/control), log2", "p-value (30 min)")

    thirtymin_df <- determineChoDEGs(thirtymin)
    thirtyminDEGs <- createDEGlist(thirtymin_df)
    thirtymincomparison <- comparetoDISS(thirtyminDEGs)
    thirtymincomparison

    ##    gene   pvalue  lfc     padj direction
    ## 1  Btg2 1.375269 1.39 4.21e-02      DISS
    ## 2 Enpp2 1.310489 1.77 4.89e-02      DISS
    ## 3  Fosb 1.478562 1.59 3.32e-02      DISS
    ## 4  Junb 1.235639 1.03 5.81e-02      DISS

    # what genes are they?

    # Btg2 - BTG Anti-Proliferation Factor 2
    # Enpp2 - Ectonucleotide Pyrophosphatase/Phosphodiesterase 2
    # Fosb - Fos Proto-Oncogene, AP-1 Transcription Factor Subunit
    # Junb - Jun Proto-Oncogene, AP-1 Transcription Factor Subunit

### ten minute DEGs

    # create 10 min df and rename column file and column headings for easy analysis
    tenmin <- subset_df(S2, "RNA fold change (10 min/control), log2", "p-value (10 min)")

    tenmin_df <- determineChoDEGs(tenmin)
    tenminDEGs <- createDEGlist(tenmin_df)
    tenmincomparison <- comparetoDISS(tenminDEGs)
    tenmincomparison

    ##   gene   pvalue  lfc     padj direction
    ## 1 Btg2 1.375269 1.39 4.21e-02      DISS
    ## 2 Ier2 1.240069 1.33 5.75e-02      DISS
    ## 3 Junb 1.235639 1.03 5.81e-02      DISS

    # what genes are they?

    # Btg2 - BTG Anti-Proliferation Factor 2
    # Junb - Jun Proto-Oncogene, AP-1 Transcription Factor Subunit
    # Ier2 - Immediate Early Response 2

### five minute DEGs

    # create 5 min df and rename column file and column headings for easy analysis
    fivemin <- subset_df(S2, "RNA fold change (5 min/control), log2", "p-value (5 min)")

    fivemin_df <- determineChoDEGs(fivemin)
    fiveminDEGs <- createDEGlist(fivemin_df)
    fivemincomparison <- comparetoDISS(fiveminDEGs)
    fivemincomparison

    ##   gene   pvalue  lfc     padj direction
    ## 1 Junb 1.235639 1.03 5.81e-02      DISS

    # what genes is it?

    # Junb - Jun Proto-Oncogene, AP-1 Transcription Factor Subunit

Cahoy et al 2008
----------------

The Cahoy et al 2008 paper
<a href="http://www.jneurosci.org/content/jneuro/28/1/264.full.pdf" class="uri">http://www.jneurosci.org/content/jneuro/28/1/264.full.pdf</a>
validates a bunch of marker genes using RNAseq and IHC. I took the
information in supplementary table 1
(<a href="http://www.jneurosci.org/content/jneuro/suppl/2008/01/03/28.1.264.DC1/JN-RM-4178_Supplemental_Data.pdf" class="uri">http://www.jneurosci.org/content/jneuro/suppl/2008/01/03/28.1.264.DC1/JN-RM-4178_Supplemental_Data.pdf</a>)
of their paper to create lists of cell-type specific marker genes.

    # conversion of some marker names to gene names
    # GLT-1 = Slc1a2
    # Connexin 30  = GJB6
    # Aquaporin 4 =  Aqp4
    # Connexin 47 = GJC2
    # Neurofilament = NEFL, NEFH, NEFM
    # Synaptotagmin I = SYT1
    # KCC2 = SLC12A5

    # list of marker genes
    astrocyte_markers <- c("Slc1a2", "Gfap", "Gjb6", "Fgfr3", "Aqp4", 
                   "Aldoc")
    oligodendrocyte_markers <- c("Gjc2", "Sox10", "Mag", "Mog", "Mbp",
                         "Cspg4", "Pdgfra", "Ugt8a", "Gal3st1", "Mobp", "Mal")
    microglia_markers <- c("Cd68", "Ptprc", "Tnf")
    neuron_markers <- c("Nefl", "Nefh", "Nefm", 
                "Gabra1", "Syt1", "Slc12a5", "Snap25",
                "Kcnq2", "Sv2b")

    # make data frames of genes expression results for markers 
    marker_expression <- function(celltype, markers){
        df <- dissociation %>%
        dplyr::filter(gene %in% c(markers)) %>%
        dplyr::mutate(marker = celltype) %>%
        droplevels()
        print(df)
    }

    astrocyte <- marker_expression("astrocyte", astrocyte_markers)

    ##     gene      pvalue   lfc     padj direction    marker
    ## 1  Aldoc 0.458931608  0.93 3.48e-01   neither astrocyte
    ## 2   Aqp4 0.015725804 -0.15 9.64e-01   neither astrocyte
    ## 3  Fgfr3 0.021833547  0.16 9.51e-01   neither astrocyte
    ## 4   Gfap 0.211045576  0.71 6.15e-01   neither astrocyte
    ## 5   Gjb6 0.097233543  0.65 7.99e-01   neither astrocyte
    ## 6 Slc1a2 0.005610261  0.04 9.87e-01   neither astrocyte

    oligodendrocyte <- marker_expression("oligodendrocyte", oligodendrocyte_markers)

    ##       gene    pvalue  lfc     padj direction          marker
    ## 1    Cspg4 0.8478768 1.38 1.42e-01   neither oligodendrocyte
    ## 2  Gal3st1 0.2889825 2.10 5.14e-01   neither oligodendrocyte
    ## 3     Gjc2 1.0179365 2.39 9.60e-02      DISS oligodendrocyte
    ## 4      Mag 4.3490754 3.31 4.48e-05      DISS oligodendrocyte
    ## 5      Mal 3.6347966 3.20 2.32e-04      DISS oligodendrocyte
    ## 6      Mbp 2.0955285 1.95 8.03e-03      DISS oligodendrocyte
    ## 7     Mobp 3.3551113 2.60 4.41e-04      DISS oligodendrocyte
    ## 8      Mog 1.6436839 2.48 2.27e-02      DISS oligodendrocyte
    ## 9   Pdgfra 0.6213240 1.24 2.39e-01   neither oligodendrocyte
    ## 10   Sox10 1.2400691 2.16 5.75e-02      DISS oligodendrocyte
    ## 11   Ugt8a 0.7590945 1.68 1.74e-01   neither oligodendrocyte

    microglia <- marker_expression("microglia", microglia_markers)

    ##    gene    pvalue  lfc     padj direction    marker
    ## 1  Cd68 1.0403953 2.35 9.11e-02      DISS microglia
    ## 2 Ptprc 0.1376531 1.37 7.28e-01   neither microglia
    ## 3   Tnf 1.6554860 2.40 2.21e-02      DISS microglia

    neuron <- marker_expression("neuron", neuron_markers)

    ##      gene      pvalue   lfc     padj direction marker
    ## 1  Gabra1 0.851423539 -1.05 1.41e-01   neither neuron
    ## 2   Kcnq2 0.183099176 -0.41 6.56e-01   neither neuron
    ## 3    Nefh 0.126491732  0.59 7.47e-01   neither neuron
    ## 4    Nefl 0.100423320  0.30 7.94e-01   neither neuron
    ## 5    Nefm 0.148997007 -0.37 7.10e-01   neither neuron
    ## 6 Slc12a5 0.573153720 -0.87 2.67e-01   neither neuron
    ## 7  Snap25 0.087332131  0.37 8.18e-01   neither neuron
    ## 8    Sv2b 0.002359648 -0.07 9.95e-01   neither neuron
    ## 9    Syt1 0.118374975 -0.33 7.61e-01   neither neuron

    # combine four small data frames into one and sort
    marker_df <- rbind.data.frame(astrocyte, oligodendrocyte, 
                     microglia, neuron)
    marker_df <- marker_df %>% select(marker, gene, lfc, padj, direction) 
    marker_df$direction <- as.character(marker_df$direction)
    marker_df <- arrange(marker_df, marker, direction)
    marker_df 

    ##             marker    gene   lfc     padj direction
    ## 1        astrocyte   Aldoc  0.93 3.48e-01   neither
    ## 2        astrocyte    Aqp4 -0.15 9.64e-01   neither
    ## 3        astrocyte   Fgfr3  0.16 9.51e-01   neither
    ## 4        astrocyte    Gfap  0.71 6.15e-01   neither
    ## 5        astrocyte    Gjb6  0.65 7.99e-01   neither
    ## 6        astrocyte  Slc1a2  0.04 9.87e-01   neither
    ## 7        microglia    Cd68  2.35 9.11e-02      DISS
    ## 8        microglia     Tnf  2.40 2.21e-02      DISS
    ## 9        microglia   Ptprc  1.37 7.28e-01   neither
    ## 10          neuron  Gabra1 -1.05 1.41e-01   neither
    ## 11          neuron   Kcnq2 -0.41 6.56e-01   neither
    ## 12          neuron    Nefh  0.59 7.47e-01   neither
    ## 13          neuron    Nefl  0.30 7.94e-01   neither
    ## 14          neuron    Nefm -0.37 7.10e-01   neither
    ## 15          neuron Slc12a5 -0.87 2.67e-01   neither
    ## 16          neuron  Snap25  0.37 8.18e-01   neither
    ## 17          neuron    Sv2b -0.07 9.95e-01   neither
    ## 18          neuron    Syt1 -0.33 7.61e-01   neither
    ## 19 oligodendrocyte    Gjc2  2.39 9.60e-02      DISS
    ## 20 oligodendrocyte     Mag  3.31 4.48e-05      DISS
    ## 21 oligodendrocyte     Mal  3.20 2.32e-04      DISS
    ## 22 oligodendrocyte     Mbp  1.95 8.03e-03      DISS
    ## 23 oligodendrocyte    Mobp  2.60 4.41e-04      DISS
    ## 24 oligodendrocyte     Mog  2.48 2.27e-02      DISS
    ## 25 oligodendrocyte   Sox10  2.16 5.75e-02      DISS
    ## 26 oligodendrocyte   Cspg4  1.38 1.42e-01   neither
    ## 27 oligodendrocyte Gal3st1  2.10 5.14e-01   neither
    ## 28 oligodendrocyte  Pdgfra  1.24 2.39e-01   neither
    ## 29 oligodendrocyte   Ugt8a  1.68 1.74e-01   neither

    write.csv(marker_df, "../results/markergenes.csv", col.names = FALSE)
