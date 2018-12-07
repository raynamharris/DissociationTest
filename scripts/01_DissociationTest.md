### Identifying the effects of cellular dissociation on hippocampal transcriptomes

The sample and count information for this part is found in
`../data/GSE99765_DissociationColData.csv` and
`../data/GSE99765_DissociationCountData.csv`. You can also download
these two files from [GEO
GSE99765](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99765).

Do a little data cleaning and calculate sample size and number of genes
measured.

    ##       
    ##        CA1 CA3 DG
    ##   HOMO   3   2  2
    ##   DISS   3   2  2

    ## [1] 22485    14

I used DESeq2 (Love et al., 2014) for gene expression normalization and
quantification using the following experimental design:
`Treatment + Subfield + Treatment * Subfield`. Genes with less than 2
counts across all samples were filtered, leaving us with `dim(rld)`
number of genes for analysis of differntial expression.

    ## [1] 16709    14

    ##               100-CA1-1 100-CA1-2 100-CA1-3 100-CA3-1 100-CA3-4 100-DG-2
    ## 0610007P14Rik  4.588504  4.776456  4.853058  5.079031  5.171332 5.186176
    ## 0610009B22Rik  3.186433  3.699918  3.130234  3.603000  3.097205 3.642543
    ## 0610009L18Rik  1.776915  2.360122  1.673934  2.380225  2.128375 1.829245
    ##               100-DG-3 101-CA1-1 101-CA1-2 101-CA1-3 101-CA3-1 101-CA3-4
    ## 0610007P14Rik 5.030620  4.955177  4.217255  3.466887  5.086616  5.019215
    ## 0610009B22Rik 3.482690  3.304830  4.260588  2.711418  3.820074  3.477413
    ## 0610009L18Rik 2.487691  2.038995  1.952132  1.795738  2.240395  1.916895
    ##               101-DG-3 101-DG-4
    ## 0610007P14Rik 5.439697 4.155559
    ## 0610009B22Rik 4.200700 2.744482
    ## 0610009L18Rik 3.161307 3.000024

    ##               100-CA1-1 100-CA1-2 100-CA1-3 100-CA3-1 100-CA3-4 100-DG-2
    ## 0610007P14Rik  6.024282  6.147172  6.198913  6.358510  6.425746 6.439969
    ## 0610009B22Rik  5.500732  5.801174  5.469570  5.744608  5.451894 5.771328
    ## 0610009L18Rik  5.099770  5.428706  5.033775  5.444353  5.302348 5.116703
    ##               100-DG-3 101-CA1-1 101-CA1-2 101-CA1-3 101-CA3-1 101-CA3-4
    ## 0610007P14Rik 6.321325  6.268696  5.582529  4.748434  6.363481  6.318356
    ## 0610009B22Rik 5.671124  5.569662  6.355559  4.748434  5.880938  5.669040
    ## 0610009L18Rik 5.499883  5.255534  4.748434  4.748434  5.364748  5.165706
    ##               101-DG-3 101-DG-4
    ## 0610007P14Rik 6.872924 5.722956
    ## 0610009B22Rik 6.528940 5.152565
    ## 0610009L18Rik 6.528940 5.867002

We identified 344 genes that were differentially expressed between the
homogenized and dissociated samples at FDR p-value &lt; 0.1 (Fig 1B).

    ## [1] 484

    ## [1] 98

    ## [1] 18

    ## [1] 344

A hierarchical clustering analysis of all differentially expressed genes
does not give rise to distinct clusters that are separated by subfield
or method; however, when examining the control, homogenized samples
alone (identified with light grey boxes), the three subfields form
distinct clusters, while the dissociated samples do not cluster by
subfield (Fig. 1C).

    ## [1] 67

Volcano Plots
-------------

Craete new data frames that include fold change, pvalue, and a column
describing the direction for differential gene expression. This
“direction” will be used to color code the dots on the volcano plot.
Will also save a list of DEGs at the end.

    ## 
    ## out of 16709 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 288, 1.7%
    ## LFC < 0 (down)     : 56, 0.34%
    ## outliers [1]       : 18, 0.11%
    ## low counts [2]     : 4534, 27%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    ## [1] 344

    ## [1] 2.058771

    ##             gene           pvalue              lfc          
    ##  0610007P14Rik:    1   Min.   :0.000003   Min.   :-5.06587  
    ##  0610009B22Rik:    1   1st Qu.:0.026323   1st Qu.:-0.40829  
    ##  0610009L18Rik:    1   Median :0.073080   Median : 0.04498  
    ##  0610009O20Rik:    1   Mean   :0.185088   Mean   : 0.15206  
    ##  0610010F05Rik:    1   3rd Qu.:0.204839   3rd Qu.: 0.56869  
    ##  0610010K14Rik:    1   Max.   :6.274558   Max.   : 9.47422  
    ##  (Other)      :12151                                        
    ##       padj           direction   
    ##  Min.   :0.0000005   DISS:  138  
    ##  1st Qu.:0.6239665   HOMO:   11  
    ##  Median :0.8451229   none:12008  
    ##  Mean   :0.7499809               
    ##  3rd Qu.:0.9411891               
    ##  Max.   :0.9999928               
    ## 

    ##       gene   pvalue      lfc       padj direction
    ## 5    Asap3 1.304370 3.892399 0.04961699      DISS
    ## 66   Itgam 1.304370 1.746838 0.04961699      DISS
    ## 73    Lcp1 1.304370 2.733823 0.04961699      DISS
    ## 37 Cyp27a1 1.308201 3.582689 0.04918113      DISS
    ## 3    Arl4c 1.310489 1.794929 0.04892278      DISS
    ## 39   Dhrs3 1.310489 2.947787 0.04892278      DISS

    ## 
    ## out of 16709 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 262, 1.6%
    ## LFC < 0 (down)     : 222, 1.3%
    ## outliers [1]       : 18, 0.11%
    ## low counts [2]     : 4210, 25%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    ## log2 fold change (MLE): Subfield CA1 vs DG 
    ## Wald test p-value: Subfield CA1 vs DG 
    ## DataFrame with 3 rows and 6 columns
    ##                baseMean    log2FoldChange             lfcSE
    ##               <numeric>         <numeric>         <numeric>
    ## C1ql2  130.705376021756 -7.65382313740708 0.786125207627745
    ## Stxbp6 143.438808384453 -5.19396604691866 0.558924222312385
    ## Crlf1  40.1931561676361 -7.69969383758498 0.843507316553684
    ##                     stat               pvalue                 padj
    ##                <numeric>            <numeric>            <numeric>
    ## C1ql2  -9.73613753018261 2.11438624194895e-22 2.63896546857649e-18
    ## Stxbp6 -9.29279111474208 1.50294385386716e-20 9.37912112005798e-17
    ## Crlf1  -9.12818856040704 6.96549786332862e-20 2.89787929440682e-16

    ##             gene           pvalue               lfc         
    ##  0610007P14Rik:    1   Min.   : 0.000003   Min.   :-9.3376  
    ##  0610009B22Rik:    1   1st Qu.: 0.003253   1st Qu.:-0.5547  
    ##  0610009L18Rik:    1   Median : 0.007490   Median :-0.1297  
    ##  0610009O20Rik:    1   Mean   : 0.155269   Mean   :-0.1413  
    ##  0610010F05Rik:    1   3rd Qu.: 0.052268   3rd Qu.: 0.2950  
    ##  0610010K14Rik:    1   Max.   :17.578566   Max.   : 8.4434  
    ##  (Other)      :12475                                        
    ##       padj        direction   
    ##  Min.   :0.0000   CA1 :  138  
    ##  1st Qu.:0.8866   DG  :  138  
    ##  Median :0.9829   none:12205  
    ##  Mean   :0.8577               
    ##  3rd Qu.:0.9925               
    ##  Max.   :1.0000               
    ## 

    ## 
    ## out of 16709 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 53, 0.32%
    ## LFC < 0 (down)     : 45, 0.27%
    ## outliers [1]       : 18, 0.11%
    ## low counts [2]     : 5178, 31%
    ## (mean count < 6)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    ## 
    ## out of 16709 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1, 0.006%
    ## LFC < 0 (down)     : 17, 0.1%
    ## outliers [1]       : 18, 0.11%
    ## low counts [2]     : 8415, 50%
    ## (mean count < 21)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

Table 1: % of DEGs
------------------

    ## [1] 2.896643

    ## [1] 0.5865103

    ## [1] 0.1077264

    ## [1] 2.058771

PCA
---

This PCA gives an overview of the variability between samples using the
a large matrix of log transformed gene expression data. You can see that
the biggest difference is between DG punches and the CA1 and CA3
punches. CA1 and CA3 samples have similar transcriptomes. The control
CA1 samples have the most similar transcriptonal profiles as evidenced
by their tight clustering.

    ## [1] 40 22 14  5  4  3

![](../figures/01_dissociationtest/PCA-1.png)

    ## quartz_off_screen 
    ##                 2

PCA statistics

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Subfield     2 2812.7  1406.4   17.69 0.000365 ***
    ## Residuals   11  874.3    79.5                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ Subfield, data = pcadata)
    ## 
    ## $Subfield
    ##              diff       lwr      upr     p adj
    ## CA3-CA1  5.223963 -10.31904 20.76696 0.6467960
    ## DG-CA1  33.098277  17.55528 48.64128 0.0003454
    ## DG-CA3  27.874315  10.84781 44.90082 0.0027012

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Subfield     2  243.8   121.9   0.744  0.498
    ## Residuals   11 1801.4   163.8

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ Subfield, data = pcadata)
    ## 
    ## $Subfield
    ##              diff       lwr      upr     p adj
    ## CA3-CA1 -8.297758 -30.60826 14.01275 0.5893115
    ## DG-CA1   1.924170 -20.38633 24.23468 0.9706111
    ## DG-CA3  10.221928 -14.21801 34.66186 0.5166947

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment    1    335   335.2     1.2  0.295
    ## Residuals   12   3352   279.3

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ Treatment, data = pcadata)
    ## 
    ## $Treatment
    ##              diff       lwr     upr     p adj
    ## DISS-HOMO 9.78567 -9.678756 29.2501 0.2948438

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## Treatment    1  691.2   691.2   6.125 0.0292 *
    ## Residuals   12 1354.1   112.8                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ Treatment, data = pcadata)
    ## 
    ## $Treatment
    ##                diff       lwr       upr     p adj
    ## DISS-HOMO -14.05249 -26.42385 -1.681127 0.0292306

Next, save files for dowstream GO analysis.

    ## 
    ## FALSE  TRUE 
    ## 11813   344

    ## sign
    ##   -1    1 
    ## 6989 9720

To view a histogram of the p-value distibution for each constrast,
change the Rmd file to `include=TRUE` for this chunck.

Here is the corresponding Adobe Illustrator file that combines many of
the above plots.

<img src="../figures/figure1.png" width="1370" />
