This R Markdown document will walk through the analysis of hippocampal tissue prepared with two different methods. The "homogenized" samples were collected by punch then homogenized in homogenization buffer from the Promega Maxwell kit. The "dissociated samples" were also collected similarily but the cells was dissociated after being punch and before being homogenized.

#### Differential Gene Expression Plots

    ## class: DESeqDataSet 
    ## dim: 16919 14 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(16919): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(14): 100-CA1-1 100-CA1-2 ... 101-DG-3 101-DG-4
    ## colData names(11): RNAseqID Method ... Punch.Collector jobnumber

    ## class: DESeqDataSet 
    ## dim: 16919 14 
    ## metadata(1): version
    ## assays(3): counts mu cooks
    ## rownames(16919): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(37): baseMean baseVar ... deviance maxCooks
    ## colnames(14): 100-CA1-1 100-CA1-2 ... 101-DG-3 101-DG-4
    ## colData names(12): RNAseqID Method ... jobnumber sizeFactor

    ## 
    ## out of 16919 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)     : 72, 0.43% 
    ## LFC < 0 (down)   : 12, 0.071% 
    ## outliers [1]     : 53, 0.31% 
    ## low counts [2]   : 0, 0% 
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    ## [1] 84

    ## 
    ## out of 16919 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)     : 63, 0.37% 
    ## LFC < 0 (down)   : 11, 0.065% 
    ## outliers [1]     : 53, 0.31% 
    ## low counts [2]   : 7539, 45% 
    ## (mean count < 15)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    ## [1] 74

![](../figures/allregions_onlyhomodiss/DifferentialGeneExpressionAnalysis-1.png)

    ## NULL

![](../figures/allregions_onlyhomodiss/DifferentialGeneExpressionAnalysis-2.png)![](../figures/allregions_onlyhomodiss/DifferentialGeneExpressionAnalysis-3.png)

    ## [1] 213

    ## [1] 825

    ## [1] 4

    ## [1] 603

    ## null device 
    ##           1

![](../figures/allregions_onlyhomodiss/Heatmap100DEgenes-1.png)![](../figures/allregions_onlyhomodiss/Heatmap100DEgenes-2.png)

    ##                  PC1        PC2             group      Method Punch
    ## 100-CA1-1 -18.375906  10.383356 Homogenized : CA1 Homogenized   CA1
    ## 100-CA1-2 -19.350417  12.248572 Homogenized : CA1 Homogenized   CA1
    ## 100-CA1-3 -17.846489  12.875882 Homogenized : CA1 Homogenized   CA1
    ## 100-CA3-1 -12.838096   2.228949 Homogenized : CA3 Homogenized   CA3
    ## 100-CA3-4  -9.622264  -3.913478 Homogenized : CA3 Homogenized   CA3
    ## 100-DG-2   12.195637  -2.324017  Homogenized : DG Homogenized    DG
    ## 100-DG-3   32.127064  19.097126  Homogenized : DG Homogenized    DG
    ## 101-CA1-1 -16.020723   5.436979 Dissociated : CA1 Dissociated   CA1
    ## 101-CA1-2  -5.940795  -5.422583 Dissociated : CA1 Dissociated   CA1
    ## 101-CA1-3   4.082342 -18.646593 Dissociated : CA1 Dissociated   CA1
    ## 101-CA3-1  -9.858726  -5.867701 Dissociated : CA3 Dissociated   CA3
    ## 101-CA3-4   3.617140 -21.517413 Dissociated : CA3 Dissociated   CA3
    ## 101-DG-3   19.282020 -20.287854  Dissociated : DG Dissociated    DG
    ## 101-DG-4   38.549214  15.708775  Dissociated : DG Dissociated    DG
    ##                name
    ## 100-CA1-1 100-CA1-1
    ## 100-CA1-2 100-CA1-2
    ## 100-CA1-3 100-CA1-3
    ## 100-CA3-1 100-CA3-1
    ## 100-CA3-4 100-CA3-4
    ## 100-DG-2   100-DG-2
    ## 100-DG-3   100-DG-3
    ## 101-CA1-1 101-CA1-1
    ## 101-CA1-2 101-CA1-2
    ## 101-CA1-3 101-CA1-3
    ## 101-CA3-1 101-CA3-1
    ## 101-CA3-4 101-CA3-4
    ## 101-DG-3   101-DG-3
    ## 101-DG-4   101-DG-4

![](../figures/allregions_onlyhomodiss/PCA-1.png)

``` r
library(edgeR)
```

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

``` r
counts <- countData
dim( counts )
```

    ## [1] 22485    14

``` r
colSums( counts ) / 1e06  # in millions of reads
```

    ## 100-CA1-1 100-CA1-2 100-CA1-3 100-CA3-1 100-CA3-4  100-DG-2  100-DG-3 
    ##  2.310696  6.646222  2.276635  1.974208  2.351650  1.285004  6.086292 
    ## 101-CA1-1 101-CA1-2 101-CA1-3 101-CA3-1 101-CA3-4  101-DG-3  101-DG-4 
    ##  4.782463  0.133622  0.300000  2.498531  1.192730  0.063507  0.598340

``` r
table( rowSums( counts ) )[ 1:30 ] # Number of genes with low counts
```

    ## 
    ##    0    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
    ## 5566  363  242  174  191  155  127  131  127  114  104  104   82   72   86 
    ##   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30 
    ##   67   76   72   53   68   45   56   63   54   54   56   54   57   55   39
