#### Differential Gene Expression

    ## class: DESeqDataSet 
    ## dim: 21871 18 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(21871): ENSMUSG00000000001 ENSMUSG00000000028 ...
    ##   __ambiguous __alignment_not_unique
    ## rowData names(0):
    ## colnames(18): dg_d_1 dg_d_2 ... ca1_v_2 ca1_v_3
    ## colData names(3): RNAseqID region location

    ## class: DESeqDataSet 
    ## dim: 21871 18 
    ## metadata(1): version
    ## assays(3): counts mu cooks
    ## rownames(21871): ENSMUSG00000000001 ENSMUSG00000000028 ...
    ##   __ambiguous __alignment_not_unique
    ## rowData names(37): baseMean baseVar ... deviance maxCooks
    ## colnames(18): dg_d_1 dg_d_2 ... ca1_v_2 ca1_v_3
    ## colData names(4): RNAseqID region location sizeFactor

    ## 
    ## out of 21871 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)     : 1427, 6.5% 
    ## LFC < 0 (down)   : 1348, 6.2% 
    ## outliers [1]     : 270, 1.2% 
    ## low counts [2]   : 0, 0% 
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    ## [1] 2775

    ## 
    ## out of 21871 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)     : 1265, 5.8% 
    ## LFC < 0 (down)   : 1226, 5.6% 
    ## outliers [1]     : 270, 1.2% 
    ## low counts [2]   : 6316, 29% 
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    ## [1] 2491

![](../figures/cembrowski/DifferentialGeneExpressionAnalysis-1.png)![](../figures/cembrowski/DifferentialGeneExpressionAnalysis-2.png)![](../figures/cembrowski/DifferentialGeneExpressionAnalysis-3.png)

    ## [1] 1775

    ## [1] 5879

    ## [1] 2557

    ## [1] 6182

    ## null device 
    ##           1

    ## Warning in rm(ann_colors): object 'ann_colors' not found

![](../figures/cembrowski/Heatmap100DEgenes-1.png)

    ##               PC1          PC2   group region location    name
    ## dg_d_1  -51.47509   5.36254855  dg : d     dg        d  dg_d_1
    ## dg_d_2  -49.24331   0.38042492  dg : d     dg        d  dg_d_2
    ## dg_d_3  -48.86525   1.24419689  dg : d     dg        d  dg_d_3
    ## dg_v_1  -39.61510   0.37638264  dg : v     dg        v  dg_v_1
    ## dg_v_2  -40.77015   0.17435536  dg : v     dg        v  dg_v_2
    ## dg_v_3  -39.94330   0.01999851  dg : v     dg        v  dg_v_3
    ## ca3_d_1  19.57195 -27.12475285 ca3 : d    ca3        d ca3_d_1
    ## ca3_d_2  23.35361 -32.79541002 ca3 : d    ca3        d ca3_d_2
    ## ca3_d_3  25.66254 -31.05716239 ca3 : d    ca3        d ca3_d_3
    ## ca3_v_1  17.35069 -20.15443305 ca3 : v    ca3        v ca3_v_1
    ## ca3_v_2  18.07738 -21.73212751 ca3 : v    ca3        v ca3_v_2
    ## ca3_v_3  21.75525 -22.72515377 ca3 : v    ca3        v ca3_v_3
    ## ca1_d_1  31.40685  27.39723618 ca1 : d    ca1        d ca1_d_1
    ## ca1_d_2  31.60661  28.17490738 ca1 : d    ca1        d ca1_d_2
    ## ca1_d_3  31.52614  28.55147154 ca1 : d    ca1        d ca1_d_3
    ## ca1_v_1  15.99951  21.95579937 ca1 : v    ca1        v ca1_v_1
    ## ca1_v_2  17.46806  22.41082658 ca1 : v    ca1        v ca1_v_2
    ## ca1_v_3  16.13362  19.54089168 ca1 : v    ca1        v ca1_v_3

![](../figures/cembrowski/PCA-1.png)

    library(edgeR)

    ## Warning: package 'edgeR' was built under R version 3.3.2

    ## Loading required package: limma

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:DESeq2':
    ## 
    ##     plotMA

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

    counts <- countData
    dim( counts )

    ## [1] 37717    18

    colSums( counts ) / 1e06  # in millions of reads

    ##    dg_d_1    dg_d_2    dg_d_3    dg_v_1    dg_v_2    dg_v_3   ca3_d_1 
    ## 15.205116 14.467908 14.523447 34.146937 29.420109 35.694203 33.670525 
    ##   ca3_d_2   ca3_d_3   ca3_v_1   ca3_v_2   ca3_v_3   ca1_d_1   ca1_d_2 
    ## 36.191022 34.646081 39.665186 39.704729 46.500576 42.928619 48.896815 
    ##   ca1_d_3   ca1_v_1   ca1_v_2   ca1_v_3 
    ## 38.984435 13.574680  8.312586 12.617995

    table( rowSums( counts ) )[ 1:30 ] # Number of genes with low counts

    ## 
    ##     0     2     3     4     5     6     7     8     9    10    11    12 
    ## 15081   765   357   370   322   249   259   223   204   167   180   181 
    ##    13    14    15    16    17    18    19    20    21    22    23    24 
    ##   153   126   131   118   121   125   104   102    93   117    84    80 
    ##    25    26    27    28    29    30 
    ##    76    84    78    61    82    87
