    library(readxl)
    library(ggplot2)
    library(cowplot)
    library(dplyr)
    library(plyr)
    library(stringr) # for sentence to upper case conversion
    library(kableExtra) # for pretty tables
    # note: won't run if plyr loaded before dplyr

    # set output file for figures and general chunk settings
    knitr::opts_chunk$set(fig.path = '../figures/04_candidategenes/', echo = T, message = F, warning=F)

Gene expression datasets or candidate gene lists in this analysis
-----------------------------------------------------------------

1.  My comparisons of homogenized and dissociated tissues
2.  Sanes and Licthman 1999
3.  Cho et al 2015
4.  Cahoy et al 2008

### My data: Dissociation Test

Import file with information about differentially expressed genes (aka
`dissociation`) and the list of genes found in the reference
transcriptome (aka `geneids`). Then filter out the non-differentially
expressed genes to create `DEGs`.

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
    kable(DEGs)

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
<tr>
<td style="text-align:left;">
IRF8
</td>
<td style="text-align:right;">
1.897994
</td>
<td style="text-align:right;">
2.77
</td>
<td style="text-align:left;">
1.26e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
JUND
</td>
<td style="text-align:right;">
1.897994
</td>
<td style="text-align:right;">
1.06
</td>
<td style="text-align:left;">
1.26e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
PLEKHB1
</td>
<td style="text-align:right;">
1.897994
</td>
<td style="text-align:right;">
1.49
</td>
<td style="text-align:left;">
1.26e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
GRIN2B
</td>
<td style="text-align:right;">
1.878071
</td>
<td style="text-align:right;">
-1.66
</td>
<td style="text-align:left;">
1.32e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
SORL1
</td>
<td style="text-align:right;">
1.871085
</td>
<td style="text-align:right;">
-1.41
</td>
<td style="text-align:left;">
1.35e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
FCRLS
</td>
<td style="text-align:right;">
2.854734
</td>
<td style="text-align:right;">
2.67
</td>
<td style="text-align:left;">
1.40e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CST3
</td>
<td style="text-align:right;">
1.837860
</td>
<td style="text-align:right;">
1.53
</td>
<td style="text-align:left;">
1.45e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RPL28
</td>
<td style="text-align:right;">
1.809330
</td>
<td style="text-align:right;">
1.17
</td>
<td style="text-align:left;">
1.55e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
UNC93B1
</td>
<td style="text-align:right;">
1.809330
</td>
<td style="text-align:right;">
2.28
</td>
<td style="text-align:left;">
1.55e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
LGMN
</td>
<td style="text-align:right;">
2.803848
</td>
<td style="text-align:right;">
1.84
</td>
<td style="text-align:left;">
1.57e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
D6WSU163E
</td>
<td style="text-align:right;">
1.798157
</td>
<td style="text-align:right;">
2.32
</td>
<td style="text-align:left;">
1.59e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SSTR1
</td>
<td style="text-align:right;">
1.798157
</td>
<td style="text-align:right;">
4.26
</td>
<td style="text-align:left;">
1.59e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
BIRC6
</td>
<td style="text-align:right;">
1.795090
</td>
<td style="text-align:right;">
-1.32
</td>
<td style="text-align:left;">
1.60e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
CDKL5
</td>
<td style="text-align:right;">
1.795090
</td>
<td style="text-align:right;">
-1.67
</td>
<td style="text-align:left;">
1.60e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
GM10401
</td>
<td style="text-align:right;">
1.795090
</td>
<td style="text-align:right;">
7.62
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
LY86
</td>
<td style="text-align:right;">
1.795090
</td>
<td style="text-align:right;">
2.81
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
OLIG1
</td>
<td style="text-align:right;">
1.795090
</td>
<td style="text-align:right;">
1.67
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
ZFP36
</td>
<td style="text-align:right;">
1.795090
</td>
<td style="text-align:right;">
1.79
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
RALGAPA1
</td>
<td style="text-align:right;">
1.787807
</td>
<td style="text-align:right;">
-1.19
</td>
<td style="text-align:left;">
1.63e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
RPS8
</td>
<td style="text-align:right;">
1.770256
</td>
<td style="text-align:right;">
1.04
</td>
<td style="text-align:left;">
1.70e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CD9
</td>
<td style="text-align:right;">
1.765684
</td>
<td style="text-align:right;">
2.36
</td>
<td style="text-align:left;">
1.72e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RPH3A
</td>
<td style="text-align:right;">
1.765684
</td>
<td style="text-align:right;">
2.26
</td>
<td style="text-align:left;">
1.72e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SLC2A1
</td>
<td style="text-align:right;">
1.752242
</td>
<td style="text-align:right;">
1.70
</td>
<td style="text-align:left;">
1.77e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CLIC4
</td>
<td style="text-align:right;">
2.742689
</td>
<td style="text-align:right;">
2.01
</td>
<td style="text-align:left;">
1.81e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
BIN2
</td>
<td style="text-align:right;">
1.738056
</td>
<td style="text-align:right;">
2.49
</td>
<td style="text-align:left;">
1.83e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
GPR84
</td>
<td style="text-align:right;">
1.728492
</td>
<td style="text-align:right;">
2.86
</td>
<td style="text-align:left;">
1.87e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
H2-K1
</td>
<td style="text-align:right;">
1.728492
</td>
<td style="text-align:right;">
2.09
</td>
<td style="text-align:left;">
1.87e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RC3H2
</td>
<td style="text-align:right;">
1.728492
</td>
<td style="text-align:right;">
-1.97
</td>
<td style="text-align:left;">
1.87e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
SLCO2B1
</td>
<td style="text-align:right;">
1.728492
</td>
<td style="text-align:right;">
2.17
</td>
<td style="text-align:left;">
1.87e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
TGFBR2
</td>
<td style="text-align:right;">
1.728492
</td>
<td style="text-align:right;">
2.44
</td>
<td style="text-align:left;">
1.87e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CDH9
</td>
<td style="text-align:right;">
1.700501
</td>
<td style="text-align:right;">
3.65
</td>
<td style="text-align:left;">
1.99e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CTSD
</td>
<td style="text-align:right;">
1.700501
</td>
<td style="text-align:right;">
1.26
</td>
<td style="text-align:left;">
1.99e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
PTGDS
</td>
<td style="text-align:right;">
1.700501
</td>
<td style="text-align:right;">
2.53
</td>
<td style="text-align:left;">
1.99e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
EPHA6
</td>
<td style="text-align:right;">
1.697706
</td>
<td style="text-align:right;">
-1.70
</td>
<td style="text-align:left;">
2.01e-02
</td>
<td style="text-align:left;">
HOMO
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
MT1
</td>
<td style="text-align:right;">
1.697706
</td>
<td style="text-align:right;">
1.42
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
INPP5D
</td>
<td style="text-align:right;">
1.690780
</td>
<td style="text-align:right;">
2.84
</td>
<td style="text-align:left;">
2.04e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ICAM1
</td>
<td style="text-align:right;">
1.685250
</td>
<td style="text-align:right;">
2.87
</td>
<td style="text-align:left;">
2.06e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
KCNMA1
</td>
<td style="text-align:right;">
1.685250
</td>
<td style="text-align:right;">
-1.92
</td>
<td style="text-align:left;">
2.06e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
MYRF
</td>
<td style="text-align:right;">
1.655486
</td>
<td style="text-align:right;">
2.21
</td>
<td style="text-align:left;">
2.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
NOS1AP
</td>
<td style="text-align:right;">
1.655486
</td>
<td style="text-align:right;">
-1.23
</td>
<td style="text-align:left;">
2.21e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
SIGLECH
</td>
<td style="text-align:right;">
1.654668
</td>
<td style="text-align:right;">
2.11
</td>
<td style="text-align:left;">
2.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SPARC
</td>
<td style="text-align:right;">
1.655486
</td>
<td style="text-align:right;">
1.79
</td>
<td style="text-align:left;">
2.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
TNF
</td>
<td style="text-align:right;">
1.655486
</td>
<td style="text-align:right;">
2.40
</td>
<td style="text-align:left;">
2.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
LAPTM5
</td>
<td style="text-align:right;">
3.649884
</td>
<td style="text-align:right;">
2.31
</td>
<td style="text-align:left;">
2.24e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
PLAU
</td>
<td style="text-align:right;">
3.649884
</td>
<td style="text-align:right;">
3.90
</td>
<td style="text-align:left;">
2.24e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
PLP1
</td>
<td style="text-align:right;">
3.649884
</td>
<td style="text-align:right;">
2.71
</td>
<td style="text-align:left;">
2.24e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MT-CO3
</td>
<td style="text-align:right;">
1.645408
</td>
<td style="text-align:right;">
1.40
</td>
<td style="text-align:left;">
2.26e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SERPINB1A
</td>
<td style="text-align:right;">
1.645917
</td>
<td style="text-align:right;">
4.36
</td>
<td style="text-align:left;">
2.26e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MOG
</td>
<td style="text-align:right;">
1.643684
</td>
<td style="text-align:right;">
2.48
</td>
<td style="text-align:left;">
2.27e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MAL
</td>
<td style="text-align:right;">
3.634797
</td>
<td style="text-align:right;">
3.20
</td>
<td style="text-align:left;">
2.32e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
PLEK
</td>
<td style="text-align:right;">
3.634797
</td>
<td style="text-align:right;">
2.50
</td>
<td style="text-align:left;">
2.32e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
STAT3
</td>
<td style="text-align:right;">
1.618314
</td>
<td style="text-align:right;">
1.56
</td>
<td style="text-align:left;">
2.41e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
TLR7
</td>
<td style="text-align:right;">
1.615572
</td>
<td style="text-align:right;">
3.89
</td>
<td style="text-align:left;">
2.42e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
PLPPR4
</td>
<td style="text-align:right;">
1.603397
</td>
<td style="text-align:right;">
-1.25
</td>
<td style="text-align:left;">
2.49e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
CELSR2
</td>
<td style="text-align:right;">
1.595174
</td>
<td style="text-align:right;">
-1.13
</td>
<td style="text-align:left;">
2.54e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
LGALS9
</td>
<td style="text-align:right;">
1.581335
</td>
<td style="text-align:right;">
2.74
</td>
<td style="text-align:left;">
2.62e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
TLN1
</td>
<td style="text-align:right;">
1.581335
</td>
<td style="text-align:right;">
1.34
</td>
<td style="text-align:left;">
2.62e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CYBA
</td>
<td style="text-align:right;">
1.579646
</td>
<td style="text-align:right;">
3.80
</td>
<td style="text-align:left;">
2.63e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ERMN
</td>
<td style="text-align:right;">
1.570140
</td>
<td style="text-align:right;">
2.53
</td>
<td style="text-align:left;">
2.69e-02
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
MT-CO1
</td>
<td style="text-align:right;">
1.570140
</td>
<td style="text-align:right;">
1.24
</td>
<td style="text-align:left;">
2.69e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SIPA1
</td>
<td style="text-align:right;">
1.570140
</td>
<td style="text-align:right;">
2.71
</td>
<td style="text-align:left;">
2.69e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
LRRC7
</td>
<td style="text-align:right;">
1.563559
</td>
<td style="text-align:right;">
-1.78
</td>
<td style="text-align:left;">
2.73e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
VASP
</td>
<td style="text-align:right;">
1.563559
</td>
<td style="text-align:right;">
2.31
</td>
<td style="text-align:left;">
2.73e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SHROOM3
</td>
<td style="text-align:right;">
1.553443
</td>
<td style="text-align:right;">
6.23
</td>
<td style="text-align:left;">
2.80e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MT-CYTB
</td>
<td style="text-align:right;">
1.544161
</td>
<td style="text-align:right;">
1.27
</td>
<td style="text-align:left;">
2.86e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
QDPR
</td>
<td style="text-align:right;">
1.544161
</td>
<td style="text-align:right;">
1.55
</td>
<td style="text-align:left;">
2.86e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SLC12A2
</td>
<td style="text-align:right;">
1.544161
</td>
<td style="text-align:right;">
2.07
</td>
<td style="text-align:left;">
2.86e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
TMEM170
</td>
<td style="text-align:right;">
1.532416
</td>
<td style="text-align:right;">
2.53
</td>
<td style="text-align:left;">
2.93e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SEMA5A
</td>
<td style="text-align:right;">
2.526452
</td>
<td style="text-align:right;">
3.25
</td>
<td style="text-align:left;">
2.98e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
APOD
</td>
<td style="text-align:right;">
1.523350
</td>
<td style="text-align:right;">
2.05
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
CADM2
</td>
<td style="text-align:right;">
1.523350
</td>
<td style="text-align:right;">
-1.21
</td>
<td style="text-align:left;">
3.00e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
CCL4
</td>
<td style="text-align:right;">
1.523350
</td>
<td style="text-align:right;">
1.94
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
ERDR1
</td>
<td style="text-align:right;">
1.523350
</td>
<td style="text-align:right;">
1.33
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
PLEKHG3
</td>
<td style="text-align:right;">
1.523350
</td>
<td style="text-align:right;">
2.65
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
SEPT4
</td>
<td style="text-align:right;">
1.523350
</td>
<td style="text-align:right;">
1.84
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
C5AR1
</td>
<td style="text-align:right;">
2.521755
</td>
<td style="text-align:right;">
3.66
</td>
<td style="text-align:left;">
3.01e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CMTM5
</td>
<td style="text-align:right;">
2.521755
</td>
<td style="text-align:right;">
2.66
</td>
<td style="text-align:left;">
3.01e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CTCFL
</td>
<td style="text-align:right;">
2.516209
</td>
<td style="text-align:right;">
3.08
</td>
<td style="text-align:left;">
3.05e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CRYAB
</td>
<td style="text-align:right;">
2.513033
</td>
<td style="text-align:right;">
2.81
</td>
<td style="text-align:left;">
3.07e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CD82
</td>
<td style="text-align:right;">
1.494433
</td>
<td style="text-align:right;">
2.18
</td>
<td style="text-align:left;">
3.20e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
JUN
</td>
<td style="text-align:right;">
3.495075
</td>
<td style="text-align:right;">
1.60
</td>
<td style="text-align:left;">
3.20e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CSPG5
</td>
<td style="text-align:right;">
1.493249
</td>
<td style="text-align:right;">
1.28
</td>
<td style="text-align:left;">
3.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MAGI2
</td>
<td style="text-align:right;">
1.493249
</td>
<td style="text-align:right;">
-1.27
</td>
<td style="text-align:left;">
3.21e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
UCP2
</td>
<td style="text-align:right;">
1.493249
</td>
<td style="text-align:right;">
2.49
</td>
<td style="text-align:left;">
3.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RPL26
</td>
<td style="text-align:right;">
2.492461
</td>
<td style="text-align:right;">
1.27
</td>
<td style="text-align:left;">
3.22e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SMOC2
</td>
<td style="text-align:right;">
1.485930
</td>
<td style="text-align:right;">
3.81
</td>
<td style="text-align:left;">
3.27e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MT-ND3
</td>
<td style="text-align:right;">
1.483492
</td>
<td style="text-align:right;">
1.70
</td>
<td style="text-align:left;">
3.28e-02
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
MFNG
</td>
<td style="text-align:right;">
1.478562
</td>
<td style="text-align:right;">
4.05
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
PCDHGA8
</td>
<td style="text-align:right;">
1.478562
</td>
<td style="text-align:right;">
1.42
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
EMP2
</td>
<td style="text-align:right;">
1.477974
</td>
<td style="text-align:right;">
2.35
</td>
<td style="text-align:left;">
3.33e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CD37
</td>
<td style="text-align:right;">
1.474133
</td>
<td style="text-align:right;">
4.14
</td>
<td style="text-align:left;">
3.36e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SCRG1
</td>
<td style="text-align:right;">
1.473390
</td>
<td style="text-align:right;">
2.63
</td>
<td style="text-align:left;">
3.36e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
GAB1
</td>
<td style="text-align:right;">
2.462986
</td>
<td style="text-align:right;">
2.90
</td>
<td style="text-align:left;">
3.44e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
PRR18
</td>
<td style="text-align:right;">
1.447684
</td>
<td style="text-align:right;">
2.29
</td>
<td style="text-align:left;">
3.57e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MT-ND1
</td>
<td style="text-align:right;">
1.443602
</td>
<td style="text-align:right;">
1.21
</td>
<td style="text-align:left;">
3.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RASGRP1
</td>
<td style="text-align:right;">
1.422330
</td>
<td style="text-align:right;">
-1.08
</td>
<td style="text-align:left;">
3.78e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
ROCK2
</td>
<td style="text-align:right;">
1.411287
</td>
<td style="text-align:right;">
-1.05
</td>
<td style="text-align:left;">
3.88e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
ASAP1
</td>
<td style="text-align:right;">
1.396225
</td>
<td style="text-align:right;">
-1.88
</td>
<td style="text-align:left;">
4.02e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
CD274
</td>
<td style="text-align:right;">
1.385171
</td>
<td style="text-align:right;">
4.24
</td>
<td style="text-align:left;">
4.12e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RPL32
</td>
<td style="text-align:right;">
1.385171
</td>
<td style="text-align:right;">
1.33
</td>
<td style="text-align:left;">
4.12e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SYNPR
</td>
<td style="text-align:right;">
1.385171
</td>
<td style="text-align:right;">
2.31
</td>
<td style="text-align:left;">
4.12e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CACNA2D1
</td>
<td style="text-align:right;">
1.380369
</td>
<td style="text-align:right;">
-1.22
</td>
<td style="text-align:left;">
4.17e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
POLD1
</td>
<td style="text-align:right;">
1.380369
</td>
<td style="text-align:right;">
2.23
</td>
<td style="text-align:left;">
4.17e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
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
CRTC2
</td>
<td style="text-align:right;">
1.359374
</td>
<td style="text-align:right;">
1.78
</td>
<td style="text-align:left;">
4.37e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
GM42878
</td>
<td style="text-align:right;">
1.359987
</td>
<td style="text-align:right;">
-3.80
</td>
<td style="text-align:left;">
4.37e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
KCNJ6
</td>
<td style="text-align:right;">
1.359987
</td>
<td style="text-align:right;">
-2.39
</td>
<td style="text-align:left;">
4.37e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
C1QC
</td>
<td style="text-align:right;">
2.356667
</td>
<td style="text-align:right;">
2.60
</td>
<td style="text-align:left;">
4.40e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MOBP
</td>
<td style="text-align:right;">
3.355111
</td>
<td style="text-align:right;">
2.60
</td>
<td style="text-align:left;">
4.41e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CNP
</td>
<td style="text-align:right;">
4.349075
</td>
<td style="text-align:right;">
2.45
</td>
<td style="text-align:left;">
4.48e-05
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
IL1A
</td>
<td style="text-align:right;">
4.349075
</td>
<td style="text-align:right;">
3.06
</td>
<td style="text-align:left;">
4.48e-05
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MAG
</td>
<td style="text-align:right;">
4.349075
</td>
<td style="text-align:right;">
3.31
</td>
<td style="text-align:left;">
4.48e-05
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MSN
</td>
<td style="text-align:right;">
1.342431
</td>
<td style="text-align:right;">
1.75
</td>
<td style="text-align:left;">
4.55e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MT-ND4
</td>
<td style="text-align:right;">
1.342431
</td>
<td style="text-align:right;">
1.15
</td>
<td style="text-align:left;">
4.55e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RTN2
</td>
<td style="text-align:right;">
1.342431
</td>
<td style="text-align:right;">
1.47
</td>
<td style="text-align:left;">
4.55e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MYC
</td>
<td style="text-align:right;">
1.340703
</td>
<td style="text-align:right;">
2.35
</td>
<td style="text-align:left;">
4.56e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
FA2H
</td>
<td style="text-align:right;">
2.336664
</td>
<td style="text-align:right;">
3.30
</td>
<td style="text-align:left;">
4.61e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
BCAS1
</td>
<td style="text-align:right;">
1.328015
</td>
<td style="text-align:right;">
1.70
</td>
<td style="text-align:left;">
4.70e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ADCY9
</td>
<td style="text-align:right;">
1.323975
</td>
<td style="text-align:right;">
-1.42
</td>
<td style="text-align:left;">
4.74e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
GATM
</td>
<td style="text-align:right;">
1.323975
</td>
<td style="text-align:right;">
2.12
</td>
<td style="text-align:left;">
4.74e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
GJB1
</td>
<td style="text-align:right;">
1.323975
</td>
<td style="text-align:right;">
4.76
</td>
<td style="text-align:left;">
4.74e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
PROS1
</td>
<td style="text-align:right;">
1.323975
</td>
<td style="text-align:right;">
2.29
</td>
<td style="text-align:left;">
4.74e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
PSAT1
</td>
<td style="text-align:right;">
1.323975
</td>
<td style="text-align:right;">
1.31
</td>
<td style="text-align:left;">
4.74e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RAB3IL1
</td>
<td style="text-align:right;">
1.320833
</td>
<td style="text-align:right;">
2.76
</td>
<td style="text-align:left;">
4.78e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
PCSK1N
</td>
<td style="text-align:right;">
3.319088
</td>
<td style="text-align:right;">
1.30
</td>
<td style="text-align:left;">
4.80e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SOCS3
</td>
<td style="text-align:right;">
3.319088
</td>
<td style="text-align:right;">
2.26
</td>
<td style="text-align:left;">
4.80e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CD14
</td>
<td style="text-align:right;">
4.311691
</td>
<td style="text-align:right;">
3.38
</td>
<td style="text-align:left;">
4.88e-05
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ARL4C
</td>
<td style="text-align:right;">
1.310489
</td>
<td style="text-align:right;">
1.79
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
DHRS3
</td>
<td style="text-align:right;">
1.310489
</td>
<td style="text-align:right;">
2.95
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
GPR34
</td>
<td style="text-align:right;">
1.310489
</td>
<td style="text-align:right;">
2.08
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
NTNG1
</td>
<td style="text-align:right;">
1.310489
</td>
<td style="text-align:right;">
3.08
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
P2RY12
</td>
<td style="text-align:right;">
1.310489
</td>
<td style="text-align:right;">
1.74
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
TLR2
</td>
<td style="text-align:right;">
1.310489
</td>
<td style="text-align:right;">
2.65
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
CYP27A1
</td>
<td style="text-align:right;">
1.308202
</td>
<td style="text-align:right;">
3.58
</td>
<td style="text-align:left;">
4.92e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
GNA12
</td>
<td style="text-align:right;">
1.308128
</td>
<td style="text-align:right;">
1.30
</td>
<td style="text-align:left;">
4.92e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ASAP3
</td>
<td style="text-align:right;">
1.304370
</td>
<td style="text-align:right;">
3.89
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
LCP1
</td>
<td style="text-align:right;">
1.304370
</td>
<td style="text-align:right;">
2.73
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
CLEC2L
</td>
<td style="text-align:right;">
2.299919
</td>
<td style="text-align:right;">
1.94
</td>
<td style="text-align:left;">
5.01e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
FRMD8
</td>
<td style="text-align:right;">
1.298154
</td>
<td style="text-align:right;">
2.79
</td>
<td style="text-align:left;">
5.03e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
NCF1
</td>
<td style="text-align:right;">
1.295307
</td>
<td style="text-align:right;">
1.80
</td>
<td style="text-align:left;">
5.07e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SPRY2
</td>
<td style="text-align:right;">
2.294344
</td>
<td style="text-align:right;">
1.79
</td>
<td style="text-align:left;">
5.08e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MCAM
</td>
<td style="text-align:right;">
2.290432
</td>
<td style="text-align:right;">
2.86
</td>
<td style="text-align:left;">
5.12e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
DNAJB2
</td>
<td style="text-align:right;">
1.288835
</td>
<td style="text-align:right;">
1.52
</td>
<td style="text-align:left;">
5.14e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MT-ATP6
</td>
<td style="text-align:right;">
1.288599
</td>
<td style="text-align:right;">
1.25
</td>
<td style="text-align:left;">
5.15e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MT-CO2
</td>
<td style="text-align:right;">
1.287949
</td>
<td style="text-align:right;">
1.36
</td>
<td style="text-align:left;">
5.15e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
OLFML3
</td>
<td style="text-align:right;">
2.285323
</td>
<td style="text-align:right;">
2.09
</td>
<td style="text-align:left;">
5.18e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
WASF2
</td>
<td style="text-align:right;">
1.284085
</td>
<td style="text-align:right;">
1.67
</td>
<td style="text-align:left;">
5.20e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
TRF
</td>
<td style="text-align:right;">
6.274558
</td>
<td style="text-align:right;">
2.72
</td>
<td style="text-align:left;">
5.31e-07
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CCDC3
</td>
<td style="text-align:right;">
1.273796
</td>
<td style="text-align:right;">
5.38
</td>
<td style="text-align:left;">
5.32e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
FAT3
</td>
<td style="text-align:right;">
1.267581
</td>
<td style="text-align:right;">
-2.18
</td>
<td style="text-align:left;">
5.40e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
SMOX
</td>
<td style="text-align:right;">
1.267581
</td>
<td style="text-align:right;">
1.58
</td>
<td style="text-align:left;">
5.40e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
TAOK1
</td>
<td style="text-align:right;">
1.267581
</td>
<td style="text-align:right;">
-1.31
</td>
<td style="text-align:left;">
5.40e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
TMEM63A
</td>
<td style="text-align:right;">
1.267581
</td>
<td style="text-align:right;">
2.46
</td>
<td style="text-align:left;">
5.40e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
VSIR
</td>
<td style="text-align:right;">
1.267581
</td>
<td style="text-align:right;">
2.06
</td>
<td style="text-align:left;">
5.40e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
GPR17
</td>
<td style="text-align:right;">
1.263140
</td>
<td style="text-align:right;">
2.60
</td>
<td style="text-align:left;">
5.46e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
PLLP
</td>
<td style="text-align:right;">
3.262969
</td>
<td style="text-align:right;">
2.76
</td>
<td style="text-align:left;">
5.46e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ADGRL3
</td>
<td style="text-align:right;">
1.244328
</td>
<td style="text-align:right;">
-1.35
</td>
<td style="text-align:left;">
5.70e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
CMTM6
</td>
<td style="text-align:right;">
1.240069
</td>
<td style="text-align:right;">
1.70
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
DOCK8
</td>
<td style="text-align:right;">
1.240069
</td>
<td style="text-align:right;">
3.07
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
GP1BB
</td>
<td style="text-align:right;">
1.240069
</td>
<td style="text-align:right;">
3.10
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
IL16
</td>
<td style="text-align:right;">
1.240069
</td>
<td style="text-align:right;">
3.32
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
RFTN1
</td>
<td style="text-align:right;">
1.240069
</td>
<td style="text-align:right;">
3.55
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
SOX10
</td>
<td style="text-align:right;">
1.240069
</td>
<td style="text-align:right;">
2.16
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
SEPP1
</td>
<td style="text-align:right;">
2.238045
</td>
<td style="text-align:right;">
1.54
</td>
<td style="text-align:left;">
5.78e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
NDRG1
</td>
<td style="text-align:right;">
3.236922
</td>
<td style="text-align:right;">
2.87
</td>
<td style="text-align:left;">
5.80e-04
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
<tr>
<td style="text-align:left;">
D630003M21RIK
</td>
<td style="text-align:right;">
1.229676
</td>
<td style="text-align:right;">
3.12
</td>
<td style="text-align:left;">
5.89e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
KCNK6
</td>
<td style="text-align:right;">
1.214729
</td>
<td style="text-align:right;">
2.84
</td>
<td style="text-align:left;">
6.10e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SLC4A8
</td>
<td style="text-align:right;">
1.214729
</td>
<td style="text-align:right;">
-1.04
</td>
<td style="text-align:left;">
6.10e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
RELA
</td>
<td style="text-align:right;">
1.208984
</td>
<td style="text-align:right;">
1.52
</td>
<td style="text-align:left;">
6.18e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MT-ND4L
</td>
<td style="text-align:right;">
1.206898
</td>
<td style="text-align:right;">
1.06
</td>
<td style="text-align:left;">
6.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
PTMA
</td>
<td style="text-align:right;">
2.207107
</td>
<td style="text-align:right;">
1.19
</td>
<td style="text-align:left;">
6.21e-03
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
ZFP521
</td>
<td style="text-align:right;">
1.200986
</td>
<td style="text-align:right;">
2.72
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
ILKAP
</td>
<td style="text-align:right;">
1.199962
</td>
<td style="text-align:right;">
1.61
</td>
<td style="text-align:left;">
6.31e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
NCKAP1L
</td>
<td style="text-align:right;">
1.199962
</td>
<td style="text-align:right;">
2.11
</td>
<td style="text-align:left;">
6.31e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MPEG1
</td>
<td style="text-align:right;">
4.189362
</td>
<td style="text-align:right;">
2.42
</td>
<td style="text-align:left;">
6.47e-05
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CCRL2
</td>
<td style="text-align:right;">
1.188123
</td>
<td style="text-align:right;">
3.01
</td>
<td style="text-align:left;">
6.48e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
GJC3
</td>
<td style="text-align:right;">
1.188123
</td>
<td style="text-align:right;">
1.85
</td>
<td style="text-align:left;">
6.48e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CSRNP3
</td>
<td style="text-align:right;">
2.183588
</td>
<td style="text-align:right;">
-1.46
</td>
<td style="text-align:left;">
6.55e-03
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
CSRNP1
</td>
<td style="text-align:right;">
1.178821
</td>
<td style="text-align:right;">
1.52
</td>
<td style="text-align:left;">
6.62e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
KLF4
</td>
<td style="text-align:right;">
1.178821
</td>
<td style="text-align:right;">
1.74
</td>
<td style="text-align:left;">
6.62e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
TBC1D16
</td>
<td style="text-align:right;">
1.176886
</td>
<td style="text-align:right;">
1.48
</td>
<td style="text-align:left;">
6.65e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
WSCD1
</td>
<td style="text-align:right;">
1.174797
</td>
<td style="text-align:right;">
1.50
</td>
<td style="text-align:left;">
6.69e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ATXN1
</td>
<td style="text-align:right;">
1.171603
</td>
<td style="text-align:right;">
-1.52
</td>
<td style="text-align:left;">
6.74e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
TMEM88B
</td>
<td style="text-align:right;">
4.169077
</td>
<td style="text-align:right;">
3.14
</td>
<td style="text-align:left;">
6.78e-05
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ARFGEF3
</td>
<td style="text-align:right;">
1.165241
</td>
<td style="text-align:right;">
-1.52
</td>
<td style="text-align:left;">
6.84e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
NFKBIA
</td>
<td style="text-align:right;">
4.165023
</td>
<td style="text-align:right;">
2.10
</td>
<td style="text-align:left;">
6.84e-05
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ARHGAP32
</td>
<td style="text-align:right;">
1.164014
</td>
<td style="text-align:right;">
-1.28
</td>
<td style="text-align:left;">
6.85e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
ATF5
</td>
<td style="text-align:right;">
1.164014
</td>
<td style="text-align:right;">
1.58
</td>
<td style="text-align:left;">
6.85e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
DUSP10
</td>
<td style="text-align:right;">
1.164014
</td>
<td style="text-align:right;">
2.49
</td>
<td style="text-align:left;">
6.85e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
PRDM5
</td>
<td style="text-align:right;">
1.164014
</td>
<td style="text-align:right;">
3.59
</td>
<td style="text-align:left;">
6.85e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
PHLDB1
</td>
<td style="text-align:right;">
3.164286
</td>
<td style="text-align:right;">
2.27
</td>
<td style="text-align:left;">
6.85e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ARPC1B
</td>
<td style="text-align:right;">
1.161877
</td>
<td style="text-align:right;">
2.05
</td>
<td style="text-align:left;">
6.89e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MT-ND2
</td>
<td style="text-align:right;">
1.158319
</td>
<td style="text-align:right;">
1.19
</td>
<td style="text-align:left;">
6.95e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
PDE4D
</td>
<td style="text-align:right;">
1.158319
</td>
<td style="text-align:right;">
-1.75
</td>
<td style="text-align:left;">
6.95e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
C1QB
</td>
<td style="text-align:right;">
5.150771
</td>
<td style="text-align:right;">
2.28
</td>
<td style="text-align:left;">
7.07e-06
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
GRAP
</td>
<td style="text-align:right;">
1.150260
</td>
<td style="text-align:right;">
2.96
</td>
<td style="text-align:left;">
7.08e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RPL18
</td>
<td style="text-align:right;">
1.150260
</td>
<td style="text-align:right;">
1.04
</td>
<td style="text-align:left;">
7.08e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SUV420H2
</td>
<td style="text-align:right;">
1.150150
</td>
<td style="text-align:right;">
2.25
</td>
<td style="text-align:left;">
7.08e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
TAF1
</td>
<td style="text-align:right;">
1.145942
</td>
<td style="text-align:right;">
-1.21
</td>
<td style="text-align:left;">
7.15e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
IL6RA
</td>
<td style="text-align:right;">
1.142011
</td>
<td style="text-align:right;">
2.21
</td>
<td style="text-align:left;">
7.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
NPR1
</td>
<td style="text-align:right;">
1.139976
</td>
<td style="text-align:right;">
3.34
</td>
<td style="text-align:left;">
7.24e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CLOCK
</td>
<td style="text-align:right;">
1.139462
</td>
<td style="text-align:right;">
-1.24
</td>
<td style="text-align:left;">
7.25e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
SMTNL2
</td>
<td style="text-align:right;">
1.139462
</td>
<td style="text-align:right;">
5.83
</td>
<td style="text-align:left;">
7.25e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MT-ATP8
</td>
<td style="text-align:right;">
1.132693
</td>
<td style="text-align:right;">
1.11
</td>
<td style="text-align:left;">
7.37e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
PIEZO1
</td>
<td style="text-align:right;">
1.132693
</td>
<td style="text-align:right;">
3.06
</td>
<td style="text-align:left;">
7.37e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
PIM1
</td>
<td style="text-align:right;">
1.132318
</td>
<td style="text-align:right;">
2.41
</td>
<td style="text-align:left;">
7.37e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RPL23A
</td>
<td style="text-align:right;">
1.132693
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:left;">
7.37e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RPL31
</td>
<td style="text-align:right;">
1.132693
</td>
<td style="text-align:right;">
1.09
</td>
<td style="text-align:left;">
7.37e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ATP2B1
</td>
<td style="text-align:right;">
1.130993
</td>
<td style="text-align:right;">
-1.21
</td>
<td style="text-align:left;">
7.40e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
CD63
</td>
<td style="text-align:right;">
1.130993
</td>
<td style="text-align:right;">
1.58
</td>
<td style="text-align:left;">
7.40e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CSF1
</td>
<td style="text-align:right;">
1.130993
</td>
<td style="text-align:right;">
1.88
</td>
<td style="text-align:left;">
7.40e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
PRR36
</td>
<td style="text-align:right;">
1.130993
</td>
<td style="text-align:right;">
-1.76
</td>
<td style="text-align:left;">
7.40e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
RPLP1
</td>
<td style="text-align:right;">
1.130993
</td>
<td style="text-align:right;">
1.13
</td>
<td style="text-align:left;">
7.40e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MYO1F
</td>
<td style="text-align:right;">
1.126676
</td>
<td style="text-align:right;">
4.78
</td>
<td style="text-align:left;">
7.47e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
INSIG1
</td>
<td style="text-align:right;">
1.124876
</td>
<td style="text-align:right;">
1.61
</td>
<td style="text-align:left;">
7.50e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ITPK1
</td>
<td style="text-align:right;">
1.120462
</td>
<td style="text-align:right;">
1.53
</td>
<td style="text-align:left;">
7.58e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CHADL
</td>
<td style="text-align:right;">
1.113943
</td>
<td style="text-align:right;">
1.82
</td>
<td style="text-align:left;">
7.69e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
LHFPL2
</td>
<td style="text-align:right;">
1.114186
</td>
<td style="text-align:right;">
2.15
</td>
<td style="text-align:left;">
7.69e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
NLRP3
</td>
<td style="text-align:right;">
1.113943
</td>
<td style="text-align:right;">
2.20
</td>
<td style="text-align:left;">
7.69e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SH3GL1
</td>
<td style="text-align:right;">
1.113943
</td>
<td style="text-align:right;">
1.56
</td>
<td style="text-align:left;">
7.69e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CIDEA
</td>
<td style="text-align:right;">
1.105338
</td>
<td style="text-align:right;">
5.46
</td>
<td style="text-align:left;">
7.85e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RACK1
</td>
<td style="text-align:right;">
1.105338
</td>
<td style="text-align:right;">
1.03
</td>
<td style="text-align:left;">
7.85e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SEMA4G
</td>
<td style="text-align:right;">
1.099518
</td>
<td style="text-align:right;">
2.30
</td>
<td style="text-align:left;">
7.95e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RPL23
</td>
<td style="text-align:right;">
1.095896
</td>
<td style="text-align:right;">
1.23
</td>
<td style="text-align:left;">
8.02e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
DUSP1
</td>
<td style="text-align:right;">
2.095528
</td>
<td style="text-align:right;">
1.35
</td>
<td style="text-align:left;">
8.03e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
LAG3
</td>
<td style="text-align:right;">
2.095528
</td>
<td style="text-align:right;">
3.32
</td>
<td style="text-align:left;">
8.03e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MBP
</td>
<td style="text-align:right;">
2.095528
</td>
<td style="text-align:right;">
1.95
</td>
<td style="text-align:left;">
8.03e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ADAM17
</td>
<td style="text-align:right;">
1.094199
</td>
<td style="text-align:right;">
2.27
</td>
<td style="text-align:left;">
8.05e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
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
ATF3
</td>
<td style="text-align:right;">
2.091979
</td>
<td style="text-align:right;">
2.17
</td>
<td style="text-align:left;">
8.09e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
C1QL1
</td>
<td style="text-align:right;">
2.091979
</td>
<td style="text-align:right;">
3.10
</td>
<td style="text-align:left;">
8.09e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
HEXB
</td>
<td style="text-align:right;">
6.091588
</td>
<td style="text-align:right;">
2.35
</td>
<td style="text-align:left;">
8.10e-07
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CPM
</td>
<td style="text-align:right;">
1.088817
</td>
<td style="text-align:right;">
2.85
</td>
<td style="text-align:left;">
8.15e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
TYROBP
</td>
<td style="text-align:right;">
3.085288
</td>
<td style="text-align:right;">
3.34
</td>
<td style="text-align:left;">
8.22e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
COX4I1
</td>
<td style="text-align:right;">
1.082307
</td>
<td style="text-align:right;">
1.04
</td>
<td style="text-align:left;">
8.27e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
NECTIN1
</td>
<td style="text-align:right;">
1.082307
</td>
<td style="text-align:right;">
1.84
</td>
<td style="text-align:left;">
8.27e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
OSTF1
</td>
<td style="text-align:right;">
1.082307
</td>
<td style="text-align:right;">
2.36
</td>
<td style="text-align:left;">
8.27e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
GNAQ
</td>
<td style="text-align:right;">
1.081900
</td>
<td style="text-align:right;">
-1.18
</td>
<td style="text-align:left;">
8.28e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
MAML2
</td>
<td style="text-align:right;">
1.081900
</td>
<td style="text-align:right;">
2.19
</td>
<td style="text-align:left;">
8.28e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SYT14
</td>
<td style="text-align:right;">
1.078212
</td>
<td style="text-align:right;">
-1.50
</td>
<td style="text-align:left;">
8.35e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
AP1S1
</td>
<td style="text-align:right;">
1.076883
</td>
<td style="text-align:right;">
1.07
</td>
<td style="text-align:left;">
8.38e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
FAM124A
</td>
<td style="text-align:right;">
1.076883
</td>
<td style="text-align:right;">
2.36
</td>
<td style="text-align:left;">
8.38e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
HR
</td>
<td style="text-align:right;">
1.076883
</td>
<td style="text-align:right;">
2.45
</td>
<td style="text-align:left;">
8.38e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
NFATC1
</td>
<td style="text-align:right;">
1.076883
</td>
<td style="text-align:right;">
3.56
</td>
<td style="text-align:left;">
8.38e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RPS6KB2
</td>
<td style="text-align:right;">
1.076883
</td>
<td style="text-align:right;">
1.65
</td>
<td style="text-align:left;">
8.38e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
STK32B
</td>
<td style="text-align:right;">
1.076292
</td>
<td style="text-align:right;">
5.43
</td>
<td style="text-align:left;">
8.39e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
UNC5B
</td>
<td style="text-align:right;">
1.069620
</td>
<td style="text-align:right;">
1.96
</td>
<td style="text-align:left;">
8.52e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CDKN1A
</td>
<td style="text-align:right;">
1.068459
</td>
<td style="text-align:right;">
2.51
</td>
<td style="text-align:left;">
8.54e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ADGRG1
</td>
<td style="text-align:right;">
2.064602
</td>
<td style="text-align:right;">
1.86
</td>
<td style="text-align:left;">
8.62e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
FAM102B
</td>
<td style="text-align:right;">
1.057831
</td>
<td style="text-align:right;">
1.53
</td>
<td style="text-align:left;">
8.75e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
GAREM
</td>
<td style="text-align:right;">
1.051773
</td>
<td style="text-align:right;">
-1.39
</td>
<td style="text-align:left;">
8.88e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
PDE4C
</td>
<td style="text-align:right;">
1.045171
</td>
<td style="text-align:right;">
-4.00
</td>
<td style="text-align:left;">
9.01e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
SH3D19
</td>
<td style="text-align:right;">
1.042289
</td>
<td style="text-align:right;">
2.43
</td>
<td style="text-align:left;">
9.07e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
TSPAN2
</td>
<td style="text-align:right;">
1.042289
</td>
<td style="text-align:right;">
1.93
</td>
<td style="text-align:left;">
9.07e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CSRP1
</td>
<td style="text-align:right;">
2.042615
</td>
<td style="text-align:right;">
1.70
</td>
<td style="text-align:left;">
9.07e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
LTBP3
</td>
<td style="text-align:right;">
2.042615
</td>
<td style="text-align:right;">
2.29
</td>
<td style="text-align:left;">
9.07e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MCL1
</td>
<td style="text-align:right;">
2.042615
</td>
<td style="text-align:right;">
1.41
</td>
<td style="text-align:left;">
9.07e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RPL29
</td>
<td style="text-align:right;">
2.042615
</td>
<td style="text-align:right;">
1.16
</td>
<td style="text-align:left;">
9.07e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RPL36A
</td>
<td style="text-align:right;">
2.042615
</td>
<td style="text-align:right;">
1.57
</td>
<td style="text-align:left;">
9.07e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
STXBP5L
</td>
<td style="text-align:right;">
2.042615
</td>
<td style="text-align:right;">
-1.94
</td>
<td style="text-align:left;">
9.07e-03
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
TNFAIP6
</td>
<td style="text-align:right;">
1.041573
</td>
<td style="text-align:right;">
2.76
</td>
<td style="text-align:left;">
9.09e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CD68
</td>
<td style="text-align:right;">
1.040395
</td>
<td style="text-align:right;">
2.35
</td>
<td style="text-align:left;">
9.11e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
MMP15
</td>
<td style="text-align:right;">
1.036560
</td>
<td style="text-align:right;">
1.94
</td>
<td style="text-align:left;">
9.19e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
LAMB2
</td>
<td style="text-align:right;">
1.035656
</td>
<td style="text-align:right;">
2.70
</td>
<td style="text-align:left;">
9.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SELPLG
</td>
<td style="text-align:right;">
6.035454
</td>
<td style="text-align:right;">
2.97
</td>
<td style="text-align:left;">
9.22e-07
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
IL10RA
</td>
<td style="text-align:right;">
1.032086
</td>
<td style="text-align:right;">
2.49
</td>
<td style="text-align:left;">
9.29e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SLC15A3
</td>
<td style="text-align:right;">
4.030623
</td>
<td style="text-align:right;">
3.45
</td>
<td style="text-align:left;">
9.32e-05
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RPS27A
</td>
<td style="text-align:right;">
1.029974
</td>
<td style="text-align:right;">
1.13
</td>
<td style="text-align:left;">
9.33e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RPL35A
</td>
<td style="text-align:right;">
1.024753
</td>
<td style="text-align:right;">
1.47
</td>
<td style="text-align:left;">
9.45e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CX3CR1
</td>
<td style="text-align:right;">
3.022089
</td>
<td style="text-align:right;">
2.10
</td>
<td style="text-align:left;">
9.50e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
TMEM119
</td>
<td style="text-align:right;">
3.022089
</td>
<td style="text-align:right;">
2.54
</td>
<td style="text-align:left;">
9.50e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
TBX1
</td>
<td style="text-align:right;">
1.019646
</td>
<td style="text-align:right;">
5.68
</td>
<td style="text-align:left;">
9.56e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CSF1R
</td>
<td style="text-align:right;">
5.018577
</td>
<td style="text-align:right;">
2.13
</td>
<td style="text-align:left;">
9.58e-06
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CTSS
</td>
<td style="text-align:right;">
5.018577
</td>
<td style="text-align:right;">
2.59
</td>
<td style="text-align:left;">
9.58e-06
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SRGN
</td>
<td style="text-align:right;">
1.018385
</td>
<td style="text-align:right;">
2.41
</td>
<td style="text-align:left;">
9.59e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RPS27
</td>
<td style="text-align:right;">
2.018182
</td>
<td style="text-align:right;">
1.29
</td>
<td style="text-align:left;">
9.59e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
APOE
</td>
<td style="text-align:right;">
1.017937
</td>
<td style="text-align:right;">
1.25
</td>
<td style="text-align:left;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ATP2C1
</td>
<td style="text-align:right;">
1.017937
</td>
<td style="text-align:right;">
-1.13
</td>
<td style="text-align:left;">
9.60e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
BCHE
</td>
<td style="text-align:right;">
1.017937
</td>
<td style="text-align:right;">
3.35
</td>
<td style="text-align:left;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CACNB2
</td>
<td style="text-align:right;">
1.017937
</td>
<td style="text-align:right;">
-1.16
</td>
<td style="text-align:left;">
9.60e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
CCL9
</td>
<td style="text-align:right;">
1.017937
</td>
<td style="text-align:right;">
3.31
</td>
<td style="text-align:left;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CD81
</td>
<td style="text-align:right;">
1.017937
</td>
<td style="text-align:right;">
1.03
</td>
<td style="text-align:left;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CFH
</td>
<td style="text-align:right;">
1.017937
</td>
<td style="text-align:right;">
2.03
</td>
<td style="text-align:left;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CTSC
</td>
<td style="text-align:right;">
1.017937
</td>
<td style="text-align:right;">
2.40
</td>
<td style="text-align:left;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
DUSP26
</td>
<td style="text-align:right;">
1.017937
</td>
<td style="text-align:right;">
2.02
</td>
<td style="text-align:left;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ERBIN
</td>
<td style="text-align:right;">
1.017937
</td>
<td style="text-align:right;">
1.72
</td>
<td style="text-align:left;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
GJC2
</td>
<td style="text-align:right;">
1.017937
</td>
<td style="text-align:right;">
2.39
</td>
<td style="text-align:left;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
LRTM2
</td>
<td style="text-align:right;">
1.017937
</td>
<td style="text-align:right;">
1.80
</td>
<td style="text-align:left;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
PLIN4
</td>
<td style="text-align:right;">
1.017937
</td>
<td style="text-align:right;">
4.05
</td>
<td style="text-align:left;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RGS10
</td>
<td style="text-align:right;">
1.017937
</td>
<td style="text-align:right;">
2.18
</td>
<td style="text-align:left;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RNASET2B
</td>
<td style="text-align:right;">
1.017937
</td>
<td style="text-align:right;">
2.00
</td>
<td style="text-align:left;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
RPSA
</td>
<td style="text-align:right;">
1.017937
</td>
<td style="text-align:right;">
1.08
</td>
<td style="text-align:left;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SLC25A10
</td>
<td style="text-align:right;">
1.017937
</td>
<td style="text-align:right;">
2.94
</td>
<td style="text-align:left;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
SPRY1
</td>
<td style="text-align:right;">
1.017937
</td>
<td style="text-align:right;">
2.18
</td>
<td style="text-align:left;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
TANGO2
</td>
<td style="text-align:right;">
1.017937
</td>
<td style="text-align:right;">
1.80
</td>
<td style="text-align:left;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
IGHV14-2
</td>
<td style="text-align:right;">
2.017752
</td>
<td style="text-align:right;">
7.81
</td>
<td style="text-align:left;">
9.60e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
UBQLN1
</td>
<td style="text-align:right;">
1.011430
</td>
<td style="text-align:right;">
-1.06
</td>
<td style="text-align:left;">
9.74e-02
</td>
<td style="text-align:left;">
HOMO
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
<tr>
<td style="text-align:left;">
GADD45B
</td>
<td style="text-align:right;">
1.007298
</td>
<td style="text-align:right;">
1.49
</td>
<td style="text-align:left;">
9.83e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
C1QA
</td>
<td style="text-align:right;">
3.006750
</td>
<td style="text-align:right;">
2.32
</td>
<td style="text-align:left;">
9.85e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
ATRX
</td>
<td style="text-align:right;">
1.003647
</td>
<td style="text-align:right;">
-1.06
</td>
<td style="text-align:left;">
9.92e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
ITPR3
</td>
<td style="text-align:right;">
1.003647
</td>
<td style="text-align:right;">
3.12
</td>
<td style="text-align:left;">
9.92e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
CD83
</td>
<td style="text-align:right;">
4.000571
</td>
<td style="text-align:right;">
2.58
</td>
<td style="text-align:left;">
9.99e-05
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
</tbody>
</table>

### Sanes and licthman

Sanes and Lichtman 1999
<a href="https://www.nature.com/articles/nn0799_597" class="uri">https://www.nature.com/articles/nn0799_597</a>
is a review paper that discusses a bunch of genes that had been
implicated in long-term potentiation (LTP). They have this one gigantic
table of protein names, organized into catagories (e.g. calcium
channels, enzymes, glutamate receptors). I obtained the gene names for
as many of these molecules as I could, and put those genes into a list
called `sanesLichtman`.

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
genes (DEGs) to my experimental results to identify the interaction
between genes that are differentially expressed following fear condition
and following a chemical manipulation. *Note: log fold change in the Cho
dataset is very, very small. Only two genes have a log fold change
greater than one, and only about 10 have a log fold change less than
one. So, I use a liberal cutoff of +/- 0.3 for fold change differences.*

The Cho data is not tidy It is very wide, with fold change and pvalue
scores for gene expression at four different timepoints. I want to
subset this one dataframe into four smaller data frames, one for each
timepoint, each with column headings: gene, lfc, and pvalue.

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
