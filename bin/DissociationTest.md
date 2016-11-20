This Markdown document will walk through the analysis of hippocampal tissue prepared with two different methods. The "homogenized" samples were collected by punch then homogenized in homogenization buffer from the Promega Maxwell kit. The "dissociated samples" were also collected similarily but the cells was dissociated after being punch and before being homogenized.

#### The sample information file

#### Raw count summary stats for all 14 samples

``` r
summary(countbygene)
```

    ##    100-CA1-1         100-CA1-2         100-CA1-3         100-CA3-1       
    ##  Min.   :    0.0   Min.   :    0.0   Min.   :    0.0   Min.   :    0.00  
    ##  1st Qu.:    0.0   1st Qu.:    0.0   1st Qu.:    0.0   1st Qu.:    0.00  
    ##  Median :    9.0   Median :   28.0   Median :    9.0   Median :   10.00  
    ##  Mean   :  102.8   Mean   :  295.6   Mean   :  101.3   Mean   :   87.83  
    ##  3rd Qu.:   73.0   3rd Qu.:  226.0   3rd Qu.:   74.0   3rd Qu.:   71.00  
    ##  Max.   :31720.0   Max.   :95996.0   Max.   :24445.0   Max.   :24878.00  
    ##    100-CA3-4          100-DG-2           100-DG-3       
    ##  Min.   :    0.0   Min.   :    0.00   Min.   :     0.0  
    ##  1st Qu.:    0.0   1st Qu.:    0.00   1st Qu.:     0.0  
    ##  Median :   12.0   Median :    7.00   Median :    32.0  
    ##  Mean   :  104.6   Mean   :   57.18   Mean   :   270.7  
    ##  3rd Qu.:   83.0   3rd Qu.:   48.00   3rd Qu.:   231.0  
    ##  Max.   :42838.0   Max.   :22711.00   Max.   :100671.0  
    ##    101-CA1-1          101-CA1-2          101-CA1-3         101-CA3-1      
    ##  Min.   :     0.0   Min.   :   0.000   Min.   :   0.00   Min.   :    0.0  
    ##  1st Qu.:     0.0   1st Qu.:   0.000   1st Qu.:   0.00   1st Qu.:    0.0  
    ##  Median :    21.0   Median :   0.000   Median :   0.00   Median :   12.0  
    ##  Mean   :   212.7   Mean   :   6.007   Mean   :  13.38   Mean   :  111.1  
    ##  3rd Qu.:   145.0   3rd Qu.:   5.000   3rd Qu.:  10.00   3rd Qu.:   85.0  
    ##  Max.   :183815.0   Max.   :3478.000   Max.   :6174.00   Max.   :86004.0  
    ##    101-CA3-4           101-DG-3          101-DG-4      
    ##  Min.   :    0.00   Min.   :   0.00   Min.   :   0.00  
    ##  1st Qu.:    0.00   1st Qu.:   0.00   1st Qu.:   0.00  
    ##  Median :    5.00   Median :   0.00   Median :   0.00  
    ##  Mean   :   53.06   Mean   :   2.93   Mean   :  26.63  
    ##  3rd Qu.:   37.00   3rd Qu.:   2.00   3rd Qu.:  20.00  
    ##  Max.   :37665.00   Max.   :4563.00   Max.   :9988.00

#### Differential Gene Expression Plots

    ## class: DESeqDataSet 
    ## dim: 17013 14 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(17013): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(14): 100-CA1-1 100-CA1-2 ... 101-DG-3 101-DG-4
    ## colData names(11): RNAseqID Method ... Date jobnumber

![](../figures/Differential%20Gene%20Expression%20Analysis-1.png)

    ## NULL

![](../figures/Differential%20Gene%20Expression%20Analysis-2.png)

![](../figures/VennDiagram-1.png)

    ## 'data.frame':    14 obs. of  2 variables:
    ##  $ Method: Factor w/ 2 levels "dissociated",..: 2 2 2 2 2 2 2 1 1 1 ...
    ##  $ Punch : Factor w/ 3 levels "CA1","CA3","DG": 1 1 1 2 2 3 3 1 1 1 ...

![](../figures/Heatmap100DEgenes-1.png)

    ##                  PC1         PC2             group      Method Punch
    ## 100-CA1-1 -17.622093   9.1439946 homogenized : CA1 homogenized   CA1
    ## 100-CA1-2 -19.159655  11.3425943 homogenized : CA1 homogenized   CA1
    ## 100-CA1-3 -17.382752  11.6431165 homogenized : CA1 homogenized   CA1
    ## 100-CA3-1 -12.090892   1.3650739 homogenized : CA3 homogenized   CA3
    ## 100-CA3-4  -8.341418  -3.5590556 homogenized : CA3 homogenized   CA3
    ## 100-DG-2   11.632917  -0.5957455  homogenized : DG homogenized    DG
    ## 100-DG-3   29.091124  21.2715840  homogenized : DG homogenized    DG
    ## 101-CA1-1 -14.885865   3.9166594 dissociated : CA1 dissociated   CA1
    ## 101-CA1-2  -4.048188  -6.9240528 dissociated : CA1 dissociated   CA1
    ## 101-CA1-3   5.474944 -15.9714610 dissociated : CA1 dissociated   CA1
    ## 101-CA3-1  -8.468865  -6.6823684 dissociated : CA3 dissociated   CA3
    ## 101-CA3-4   4.597856 -18.3067083 dissociated : CA3 dissociated   CA3
    ## 101-DG-3   18.039195 -21.2504238  dissociated : DG dissociated    DG
    ## 101-DG-4   33.163694  14.6067925  dissociated : DG dissociated    DG
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

![](../figures/PCA-1.png)
