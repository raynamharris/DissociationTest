Comparing candidate genes lists
===============================

The purpose of this script is to determine if there is any overlap in
the genes identified as differentially expressed in my study (aka
Dissociation test) and other experiments. I compare my data to Sanes and
Lichtman 1999 and Cho et a.l 2015 to see if these genes are implicated
in long-term potention or fear-condition, respectively. I compare it to
Cahoy et al. 2008 to see if they overlap cell-type specific marker
genes.

Data sets used in these comparisons.

1.  My data: Dissociation Test
2.  Sanes and Lichtman 1999
3.  Cho et al 2015
4.  Cahoy et al 2008

My data: Dissociation Test
--------------------------

Import file with information about differentially expressed genes (aka
`dissociation`) and the list of genes found in the reference
transcriptome (aka `geneids`). Then filter out the non-differentially
expressed genes to create `DEGs`. In this table “direction” refers to
wether genes were upregulated in the dissocated (DISS) or homogenized
(HOMO) samples.

    dissociation <- read.csv("../results/volcanoTreatment.csv", header = T, row.names = 1)
    dissociation$lfc <- round(dissociation$lfc,2)
    dissociation$padj <- formatC(dissociation$padj, format = "e", digits = 2)
    dissociation$gene <- str_to_upper(dissociation$gene)

    geneids <- read.csv("../data/geneids.csv", header = T)
    geneids$gene <- str_to_upper(geneids$gene)

    # Filter out non-significant genes, sort by p-value, rename column, round to 2 decimal places. 

    DEGs <- dissociation %>%
      dplyr::filter(direction != "neither") %>%
      arrange((padj))
    kable(head(DEGs))

<table>
<thead>
<tr>
<th style="text-align:left;">
gene
</th>
<th style="text-align:right;">
pvalue
</th>
<th style="text-align:right;">
lfc
</th>
<th style="text-align:left;">
padj
</th>
<th style="text-align:left;">
direction
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
GABRB1
</td>
<td style="text-align:right;">
1.987167
</td>
<td style="text-align:right;">
-1.07
</td>
<td style="text-align:left;">
1.03e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
CLDN11
</td>
<td style="text-align:right;">
3.979408
</td>
<td style="text-align:right;">
3.14
</td>
<td style="text-align:left;">
1.05e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RGS7BP
</td>
<td style="text-align:right;">
1.942452
</td>
<td style="text-align:right;">
-1.36
</td>
<td style="text-align:left;">
1.14e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
TREM2
</td>
<td style="text-align:right;">
1.942452
</td>
<td style="text-align:right;">
3.25
</td>
<td style="text-align:left;">
1.14e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SYNGR2
</td>
<td style="text-align:right;">
1.922096
</td>
<td style="text-align:right;">
3.35
</td>
<td style="text-align:left;">
1.20e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SLC2A5
</td>
<td style="text-align:right;">
2.911672
</td>
<td style="text-align:right;">
4.08
</td>
<td style="text-align:left;">
1.23e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
</tbody>
</table>

Sanes and Licthman 1999
-----------------------

Sanes and Lichtman 1999
<a href="https://www.nature.com/articles/nn0799_597" class="uri">https://www.nature.com/articles/nn0799_597</a>
is a review paper that discusses a bunch of genes that had been
implicated in long-term potentiation (LTP). They have this one gigantic
table of protein names, organized into categories (e.g. calcium
channels, enzymes, glutamate receptors). I obtained the gene names for
as many of these molecules as I could, and put those genes into a list
called `sanesLichtman`.

    # Compare to Molecules implicated in hippocampal LTP from sanes

    supptable2 <- read.csv("../data/SanesLichtman.csv", check.names = F)
    headsupp <- head(supptable2, 10)

    headsupp %>%
      kable() 

<table>
<thead>
<tr>
<th style="text-align:left;">
Sanes & Lichtman Molecules
</th>
<th style="text-align:left;">
Related Transcripts
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
GLUTAMATE RECEPTORS
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
GluR1; GluR2
</td>
<td style="text-align:left;">
Gria1; Gria2
</td>
</tr>
<tr>
<td style="text-align:left;">
mGluR1; mGluR4; mGluR5; mGluR7
</td>
<td style="text-align:left;">
Grm1; Grm4; Grm5; Grm7
</td>
</tr>
<tr>
<td style="text-align:left;">
NMDA NR2A; NMDA NR2D; NMDA NR1
</td>
<td style="text-align:left;">
Grin1; Grin2a; Grin2d
</td>
</tr>
<tr>
<td style="text-align:left;">
OTHER NEUROTRANSMITTERS
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
norepinephrine and b-adrenergic receptors
</td>
<td style="text-align:left;">
Adrb1; Adrb2; Adrb3
</td>
</tr>
<tr>
<td style="text-align:left;">
adenosine and adenosine 2A receptors
</td>
<td style="text-align:left;">
Adra1a; Adra1b; Adra1d; Adra2a; Adra2b; Adra2c
</td>
</tr>
<tr>
<td style="text-align:left;">
dopamine and D1 dopamine receptors
</td>
<td style="text-align:left;">
Th; Drd1
</td>
</tr>
<tr>
<td style="text-align:left;">
mu and delta opioid receptors
</td>
<td style="text-align:left;">
Oprm1; Oprd1
</td>
</tr>
<tr>
<td style="text-align:left;">
acetylcholine receptors
</td>
<td style="text-align:left;">
Chrna1; Chrna7; Chrna3; Chrnb1; Chrnb2; Chrnb3
</td>
</tr>
</tbody>
</table>

    # create list of candidate genes Sanes and Lichtman 1999 
    # (note, list creation not show, only the alphabetized version is shown)

    sanesLichtman[order(sanesLichtman)] # print list alphabetically

    ##   [1] "ACHE"     "ADCY1"    "ADRA2A"   "ADRA2B"   "ADRA2C"   "ADRB1"   
    ##   [7] "ADRB2"    "ADRB3"    "BDNF"     "CACNA1A"  "CACNA1B"  "CACNA1C" 
    ##  [13] "CACNA1D"  "CACNA1E"  "CACNA1F"  "CACNA1S"  "CALB1"    "CALB2"   
    ##  [19] "CALM1"    "CALM2"    "CALM3"    "CAMK1"    "CAMK2"    "CAMK4"   
    ##  [25] "CAPN1"    "CAPN10"   "CAPN11"   "CAPN12"   "CAPN13"   "CAPN15"  
    ##  [31] "CAPN2"    "CAPN3"    "CAPN5"    "CAPN6"    "CAPN7"    "CAPN8"   
    ##  [37] "CAPN9"    "CAST"     "CCR7"     "CD47"     "CDH1"     "CDH2"    
    ##  [43] "CHRM1"    "CHRM2"    "CHRM3"    "CHRM4"    "CHRM5"    "CHRNA1"  
    ##  [49] "CHRNA3"   "CHRNA7"   "CHRNB1"   "CHRNB2"   "CHRNB3"   "CNGA2"   
    ##  [55] "CNR1"     "CNR2"     "CREB1"    "DLG4"     "DRD1"     "EDN1"    
    ##  [61] "EFNA5"    "EGF"      "EGR1"     "EGR2"     "EPHA5"    "ERBB4"   
    ##  [67] "FGF2"     "FYN"      "GABBR1"   "GABRA1"   "GABRA2"   "GABRA3"  
    ##  [73] "GABRA5"   "GABRA6"   "GABRB1"   "GABRB2"   "GABRB3"   "GABRR1"  
    ##  [79] "GAP43"    "GFAP"     "GRIA1"    "GRIA2"    "GRIN1"    "GRIN2A"  
    ##  [85] "GRIN2D"   "GRM1"     "GRM4"     "GRM5"     "GRM7"     "GUCY1A2" 
    ##  [91] "GUCY1A3"  "GUCY1B2"  "GUCY1B3"  "GUCY2C"   "GUCY2D"   "GUCY2E"  
    ##  [97] "GUCY2G"   "HOMER1"   "HOMER2"   "HOMER3"   "HTR1A"    "HTR1B"   
    ## [103] "HTR1F"    "HTR2A"    "HTR2B"    "HTR2C"    "HTR3A"    "HTR3B"   
    ## [109] "HTR4"     "HTR5A"    "HTR5B"    "HTR6"     "HTR7"     "ICAM5"   
    ## [115] "IL1B"     "INHBA"    "ITGA1"    "ITGA10"   "ITGA11"   "ITGA2"   
    ## [121] "ITGA2B"   "ITGA3"    "ITGA4"    "ITGA5"    "ITGA6"    "ITGA7"   
    ## [127] "ITGA8"    "ITGA9"    "ITGAD"    "ITGAE"    "ITGAL"    "ITGAM"   
    ## [133] "ITGAV"    "ITGAX"    "ITGB1"    "ITGB1BP1" "ITGB2"    "ITGB2L"  
    ## [139] "ITGB3"    "ITGB3BP"  "ITGB4"    "ITGB5"    "ITGB6"    "ITGB7"   
    ## [145] "ITGB8"    "ITGBL1"   "ITPKA"    "ITPKB"    "ITPKC"    "L1CAM"   
    ## [151] "MAPK1"    "MAPK10"   "MAPK11"   "MAPK12"   "MAPK14"   "MAPK3"   
    ## [157] "MAPK4"    "MAPK6"    "MAPK7"    "MAPK8"    "MAPK9"    "MAS1"    
    ## [163] "NAPA"     "NCAM1"    "NCAM2"    "NGF"      "NOS1"     "NOS3"    
    ## [169] "NRG1"     "NRG2"     "NRG3"     "NRGN"     "NTRK2"    "OPRD1"   
    ## [175] "OPRL1"    "OPRM1"    "PARP1"    "PHPT1"    "PLA2G10"  "PLA2G12A"
    ## [181] "PLA2G12B" "PLA2G15"  "PLA2G16"  "PLA2G1B"  "PLA2G2A"  "PLA2G2C" 
    ## [187] "PLA2G2D"  "PLA2G2E"  "PLA2G2F"  "PLA2G3"   "PLA2G4A"  "PLA2G4B" 
    ## [193] "PLA2G4E"  "PLA2G4F"  "PLA2G5"   "PLA2G6"   "PLA2G7"   "PLAT"    
    ## [199] "PLCB1"    "PLCB2"    "PLCB3"    "PLCB4"    "PLCG1"    "PLCG2"   
    ## [205] "PLG"      "PNOC"     "PPP3CA"   "PPP3CB"   "PPP3CC"   "PRKACB"  
    ## [211] "PRKAR1B"  "PRKCG"    "PRKCZ"    "PRKG1"    "PTN"      "RAB3A"   
    ## [217] "RARB"     "S100B"    "SERPINE2" "SNAP25"   "SPTAN1"   "SPTBN1"  
    ## [223] "SRC"      "STMN4"    "STX1B"    "SYN1"     "SYP"      "TH"      
    ## [229] "THY1"     "TNC"      "UBE3A"    "VAMP1"    "VAMP2"    "VAMP3"   
    ## [235] "VAMP4"    "VAMP5"    "VAMP8"

After checking to see how many of the Sanes Lichtman genes are in my
reference transcriptome and then in my dataset, I determine which ones
are differentially expressed.

    # confirm that all all Sanes and Lichtman genes are in the reference transcriptome
    sanesLichtman_reference <- geneids %>%
      dplyr::filter(gene %in% sanesLichtman) %>%
      dplyr::select(gene) %>%
      distinct() %>%
      arrange(gene) %>%
      droplevels()
    sanesLichtman_reference <- sanesLichtman_reference[,c(1)]
    str(sanesLichtman_reference)

    ##  chr [1:236] "ACHE" "ADCY1" "ADRA2A" "ADRA2B" "ADRA2C" "ADRB1" "ADRB2" ...

    # identify which of the Sanes and Lichtman genes are present in my samples
    sanesLichtman_present <- dissociation %>%
      dplyr::filter(gene %in% sanesLichtman) %>%
      droplevels()
    sanesLichtman_present <- sanesLichtman_present[,c(1)]
    str(sanesLichtman_present)

    ##  chr [1:175] "ACHE" "ADCY1" "ADRA2A" "ADRA2B" "ADRA2C" "ADRB1" "BDNF" ...

    # identify whichof the Sanes and Lichtman genes are differentially expressed in this analysis
    sanesLichtman_DEGs <- DEGs %>%
      dplyr::filter(gene %in% sanesLichtman) %>%
        dplyr::filter(direction != "neither") %>%
      arrange(gene)
    kable(sanesLichtman_DEGs)

<table>
<thead>
<tr>
<th style="text-align:left;">
gene
</th>
<th style="text-align:right;">
pvalue
</th>
<th style="text-align:right;">
lfc
</th>
<th style="text-align:left;">
padj
</th>
<th style="text-align:left;">
direction
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
CACNA1E
</td>
<td style="text-align:right;">
1.094199
</td>
<td style="text-align:right;">
-1.28
</td>
<td style="text-align:left;">
8.05e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
GABRB1
</td>
<td style="text-align:right;">
1.987167
</td>
<td style="text-align:right;">
-1.07
</td>
<td style="text-align:left;">
1.03e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
GRIN2A
</td>
<td style="text-align:right;">
1.570140
</td>
<td style="text-align:right;">
-1.66
</td>
<td style="text-align:left;">
2.69e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
IL1B
</td>
<td style="text-align:right;">
1.523350
</td>
<td style="text-align:right;">
2.41
</td>
<td style="text-align:left;">
3.00e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ITGA5
</td>
<td style="text-align:right;">
1.795090
</td>
<td style="text-align:right;">
3.05
</td>
<td style="text-align:left;">
1.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ITGAM
</td>
<td style="text-align:right;">
1.304370
</td>
<td style="text-align:right;">
1.75
</td>
<td style="text-align:left;">
4.96e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ITGB4
</td>
<td style="text-align:right;">
1.200986
</td>
<td style="text-align:right;">
2.93
</td>
<td style="text-align:left;">
6.30e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ITGB5
</td>
<td style="text-align:right;">
1.695894
</td>
<td style="text-align:right;">
1.98
</td>
<td style="text-align:left;">
2.01e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ITPKB
</td>
<td style="text-align:right;">
1.033913
</td>
<td style="text-align:right;">
1.53
</td>
<td style="text-align:left;">
9.25e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MAPK3
</td>
<td style="text-align:right;">
2.008780
</td>
<td style="text-align:right;">
1.61
</td>
<td style="text-align:left;">
9.80e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
</tbody>
</table>

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

Cho et al 2015
--------------

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
genes (DEGs) to my experimental results to identify the interaction
between genes that are differentially expressed following fear condition
and following a chemical manipulation. *Note: log fold change in the Cho
dataset is very, very small. Only two genes have a log fold change
greater than one, and only about 10 have a log fold change less than
one. So, I use a liberal cutoff of +/- 0.3 for fold change differences.*

The Cho data is not tidy It is very wide, with fold change and p-value
scores for gene expression at four different time points. I want to
subset this one data frame into four smaller data frames, one for each
time point, each with column headings: gene, lfc, and p-value.

    # read supplemental data from cho et al
    S2 <- as.data.frame(readxl::read_excel("../data/aac7368-Cho.SM.Table.S2.xls", skip = 1 ))
    S2 <- rename(S2, c(`Gene Symbol` = "gene")) # rename gene column
    S2$gene <- str_to_upper(S2$gene)
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

1.  `subset_df`: subset the columns of interest and rename them
2.  `determineChoDEGs`: categorize the genes as differentially expressed
    or not
3.  `createDEGlist`: create a list of just DEG gene names
4.  `comparetoDISS`: crossreference the list of DEG names with the
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
with those since that is about the cutoff used in the Cho paper.

    fourhours_df <- subset_df(S2, "RNA fold change (4 h/control), log2", "p-value (4 h)")
    fourhours_df <- determineChoDEGs(fourhours_df)
    fourhourDEGs <- createDEGlist(fourhours_df)
    fourhourcomparison <- comparetoDISS(fourhourDEGs)
    kable(fourhourcomparison)

<table>
<thead>
<tr>
<th style="text-align:left;">
gene
</th>
<th style="text-align:right;">
pvalue
</th>
<th style="text-align:right;">
lfc
</th>
<th style="text-align:left;">
padj
</th>
<th style="text-align:left;">
direction
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
ENPP2
</td>
<td style="text-align:right;">
1.310489
</td>
<td style="text-align:right;">
1.77
</td>
<td style="text-align:left;">
4.89e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
FN1
</td>
<td style="text-align:right;">
1.570140
</td>
<td style="text-align:right;">
1.77
</td>
<td style="text-align:left;">
2.69e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
</tbody>
</table>

    # what genes are they?

    # Enpp2 - Ectonucleotide Pyrophosphatase/Phosphodiesterase 2
    # Fn1 - Fibronectin 1 - cell adhesion

### thirty minute DEGs

    # create 30 min df and rename column file and column headings for easy analysis
    thirtymin <- subset_df(S2, "RNA fold change (30 min/control), log2", "p-value (30 min)")

    thirtymin_df <- determineChoDEGs(thirtymin)
    thirtyminDEGs <- createDEGlist(thirtymin_df)
    thirtymincomparison <- comparetoDISS(thirtyminDEGs)
    kable(thirtymincomparison)

<table>
<thead>
<tr>
<th style="text-align:left;">
gene
</th>
<th style="text-align:right;">
pvalue
</th>
<th style="text-align:right;">
lfc
</th>
<th style="text-align:left;">
padj
</th>
<th style="text-align:left;">
direction
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BTG2
</td>
<td style="text-align:right;">
1.375269
</td>
<td style="text-align:right;">
1.39
</td>
<td style="text-align:left;">
4.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ENPP2
</td>
<td style="text-align:right;">
1.310489
</td>
<td style="text-align:right;">
1.77
</td>
<td style="text-align:left;">
4.89e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
FOSB
</td>
<td style="text-align:right;">
1.478562
</td>
<td style="text-align:right;">
1.59
</td>
<td style="text-align:left;">
3.32e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
JUNB
</td>
<td style="text-align:right;">
1.235639
</td>
<td style="text-align:right;">
1.03
</td>
<td style="text-align:left;">
5.81e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
</tbody>
</table>

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
    kable(tenmincomparison)

<table>
<thead>
<tr>
<th style="text-align:left;">
gene
</th>
<th style="text-align:right;">
pvalue
</th>
<th style="text-align:right;">
lfc
</th>
<th style="text-align:left;">
padj
</th>
<th style="text-align:left;">
direction
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BTG2
</td>
<td style="text-align:right;">
1.375269
</td>
<td style="text-align:right;">
1.39
</td>
<td style="text-align:left;">
4.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
IER2
</td>
<td style="text-align:right;">
1.240069
</td>
<td style="text-align:right;">
1.33
</td>
<td style="text-align:left;">
5.75e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
JUNB
</td>
<td style="text-align:right;">
1.235639
</td>
<td style="text-align:right;">
1.03
</td>
<td style="text-align:left;">
5.81e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
</tbody>
</table>

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
    kable(fivemincomparison)

<table>
<thead>
<tr>
<th style="text-align:left;">
gene
</th>
<th style="text-align:right;">
pvalue
</th>
<th style="text-align:right;">
lfc
</th>
<th style="text-align:left;">
padj
</th>
<th style="text-align:left;">
direction
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
JUNB
</td>
<td style="text-align:right;">
1.235639
</td>
<td style="text-align:right;">
1.03
</td>
<td style="text-align:left;">
5.81e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
</tbody>
</table>

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

    ## [1] gene      pvalue    lfc       padj      direction marker   
    ## <0 rows> (or 0-length row.names)

    oligodendrocyte <- marker_expression("oligodendrocyte", oligodendrocyte_markers)

    ## [1] gene      pvalue    lfc       padj      direction marker   
    ## <0 rows> (or 0-length row.names)

    microglia <- marker_expression("microglia", microglia_markers)

    ## [1] gene      pvalue    lfc       padj      direction marker   
    ## <0 rows> (or 0-length row.names)

    neuron <- marker_expression("neuron", neuron_markers)

    ## [1] gene      pvalue    lfc       padj      direction marker   
    ## <0 rows> (or 0-length row.names)

    # combine four small data frames into one and sort
    marker_df <- rbind.data.frame(astrocyte, oligodendrocyte, 
                     microglia, neuron)
    marker_df <- marker_df %>% select(marker, gene, lfc, padj, direction) 
    marker_df$direction <- as.character(marker_df$direction)
    marker_df <- arrange(marker_df, marker, direction)
    kable(marker_df) 

<table>
<thead>
<tr>
<th style="text-align:left;">
marker
</th>
<th style="text-align:left;">
gene
</th>
<th style="text-align:right;">
lfc
</th>
<th style="text-align:left;">
padj
</th>
<th style="text-align:left;">
direction
</th>
</tr>
</thead>
<tbody>
<tr>
</tr>
</tbody>
</table>

    write.csv(marker_df, "../results/markergenes.csv")
