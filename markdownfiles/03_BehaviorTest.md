All together now
----------------

Combining the two previous analyses

PCA

![](../figures/03_behaviortest/PCA-1.png)![](../figures/03_behaviortest/PCA-2.png)![](../figures/03_behaviortest/PCA-3.png)

    ## DEG by contrasts
    source("resvalsfunction.R")
    contrast1 <- resvals(contrastvector = c('Region', 'CA1', 'DG'), mypval = 0.1)

    ## [1] 4657
    ## [1] 2898

    contrast2 <- resvals(contrastvector = c('Region', 'CA3', 'DG'), mypval = 0.1)

    ## [1] 5257
    ## [1] 3224

    contrast3 <- resvals(contrastvector = c('Region', 'CA1', 'CA3'), mypval = 0.1)

    ## [1] 3982
    ## [1] 1618

    contrast4 <- resvals(contrastvector = c('Group', 'trained', 'control'), mypval = 0.1)

    ## [1] 3041
    ## [1] 285

Now, we can view a histogram of the distribution

![](../figures/03_behaviortest/histogram-1.png)

    ## [1] 1

![](../figures/03_behaviortest/histogram-2.png)

    ## [1] 1

![](../figures/03_behaviortest/histogram-3.png)

    ## [1] 1

![](../figures/03_behaviortest/histogram-4.png)

    ## [1] 1

This Venn Diagram sthe overlap of differentailly expression genes by
Region and method. This shows all genes with *uncorrected* pvalue
&lt;0.1.

![](../figures/03_behaviortest/VennDiagramPVal-1.png)

This Venn Diagram sthe overlap of differentailly expression genes by
Region and method. This shows all genes with *adjusted* pvalue &lt;0.1.

![](../figures/03_behaviortest/VennDiagramPadj-1.png)

Heatmaps

![](../figures/03_behaviortest/HeatmapPadj-1.png)

![](../figures/03_behaviortest/HeatmapPvalue-1.png)

    FALSE 
    FALSE FALSE  TRUE 
    FALSE  8889   504

    FALSE 
    FALSE FALSE  TRUE 
    FALSE 14951  2000

<<<<<<< HEAD
    FALSE log2 fold change (MLE): Group trained vs control 
    FALSE Wald test p-value: Group trained vs control 
    FALSE DataFrame with 6 rows and 6 columns
    FALSE                baseMean log2FoldChange     lfcSE       stat    pvalue
    FALSE               <numeric>      <numeric> <numeric>  <numeric> <numeric>
    FALSE 0610007P14Rik 28.312171     -0.4917445 0.5416913 -0.9077947 0.3639867
    FALSE 0610009B22Rik  7.764007      1.1336107 1.1791409  0.9613870 0.3363576
    FALSE 0610009L18Rik  3.652118      2.3345425 1.9671167  1.1867840 0.2353128
    FALSE 0610009O20Rik 58.892418      0.6145420 0.3810067  1.6129426 0.1067570
    FALSE 0610010F05Rik 10.601326      0.1150367 0.7181408  0.1601868 0.8727339
    FALSE 0610010K14Rik  1.780900     -1.3461345 1.1956164 -1.1258917 0.2602114
=======
    FALSE log2 fold change (MLE): Group control vs trained 
    FALSE Wald test p-value: Group control vs trained 
    FALSE DataFrame with 6 rows and 6 columns
    FALSE                baseMean log2FoldChange     lfcSE       stat    pvalue
    FALSE               <numeric>      <numeric> <numeric>  <numeric> <numeric>
    FALSE 0610007P14Rik 28.312171      0.4917445 0.5416913  0.9077947 0.3639867
    FALSE 0610009B22Rik  7.764007     -1.1336107 1.1791409 -0.9613870 0.3363576
    FALSE 0610009L18Rik  3.652118     -2.3345425 1.9671167 -1.1867840 0.2353128
    FALSE 0610009O20Rik 58.892418     -0.6145420 0.3810067 -1.6129426 0.1067570
    FALSE 0610010F05Rik 10.601326     -0.1150367 0.7181408 -0.1601868 0.8727339
    FALSE 0610010K14Rik  1.780900      1.3461345 1.1956164  1.1258917 0.2602114
>>>>>>> 06a9ef068c8c11ec58b2295f4e31fb03700f4607
    FALSE                    padj
    FALSE               <numeric>
    FALSE 0610007P14Rik 0.6743446
    FALSE 0610009B22Rik        NA
    FALSE 0610009L18Rik        NA
    FALSE 0610009O20Rik 0.3933641
    FALSE 0610010F05Rik 0.9513276
    FALSE 0610010K14Rik        NA

    FALSE sign
    FALSE   -1    1 
<<<<<<< HEAD
    FALSE 7986 8984
=======
    FALSE 8984 7986
>>>>>>> 06a9ef068c8c11ec58b2295f4e31fb03700f4607
