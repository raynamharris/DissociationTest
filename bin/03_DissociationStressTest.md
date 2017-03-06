Methods for Dorsal Hippocampal Gene Expression Profiling
--------------------------------------------------------

#### Part 1: Examining the influence of dissasociation on gene expression in the CA1, CA3, and DG

![](../figures/03_dissociationstresstest/DifferentialGeneExpressionAnalysis-1.png)

This PCA gives an overview of the variability between samples using the
a large matrix of log transformed gene expression data. You can see that
the bigges difference is between DG punches and the CA1 and CA3 punches.
CA1 and CA3 samples have similar transcriptomes. The homogenized CA1
samples have the most similar transcriptonal profiles as evidenced by
their tight clustering.

![](../figures/03_dissociationstresstest/PCA-1.png)![](../figures/03_dissociationstresstest/PCA-2.png)

This Venn Diagram shows the number of differentially expressed by
contrast described above each oval. The most number of genes are
differntially expressed between DG and the CAs (nearly 1000) wheras only
about 200 were differntailly regulated as a result of of technical
maniplulation comparing homogenized and dissociated samples.

The first is with padj values. The second with p values

![](../figures/03_dissociationstresstest/VennDiagram1-1.png)

![](../figures/03_dissociationstresstest/VennDiagram2-1.png)

Here, the goal is the analyze the distribution of pvalues to see if they
are randomly distributed or if that is a tendency towards and increase
or decrease of low pvalues. There, I'm showing the pval and adjusted
pvale (padj) for all for two-way comparision.

    head(rldpvals)

    ##               pvalPunchCA1DG padjPunchCA1DG pvalPunchCA3DG padjPunchCA3DG
    ## 0610007P14Rik      0.4379056      1.0000000     0.30445691      0.9181681
    ## 0610009B22Rik      0.5181528      1.0000000     0.08780957      0.5721286
    ## 0610009L18Rik      0.9850306      1.0000000     0.69597687      1.0000000
    ## 0610009O20Rik      0.8413827      1.0000000     0.46609426      1.0000000
    ## 0610010F05Rik      0.3899734      0.9958327     0.88163051      1.0000000
    ## 0610010K14Rik      0.5136582      1.0000000     0.46035764      1.0000000
    ##               pvalPunchCA1CA3 padjPunchCA1CA3 pvalGrouphomecagestressed
    ## 0610007P14Rik       0.7271586       1.0000000                 0.6021758
    ## 0610009B22Rik       0.2157030       0.9859284                 0.8163708
    ## 0610009L18Rik       0.6627063       1.0000000                 0.5685804
    ## 0610009O20Rik       0.5548065       1.0000000                 0.8392834
    ## 0610010F05Rik       0.2967238       1.0000000                 0.4453358
    ## 0610010K14Rik       0.8768591       1.0000000                 0.9522420
    ##               padjGrouphomecagestressed
    ## 0610007P14Rik                 0.9999077
    ## 0610009B22Rik                 0.9999077
    ## 0610009L18Rik                 0.9999077
    ## 0610009O20Rik                 0.9999077
    ## 0610010F05Rik                 0.9999077
    ## 0610010K14Rik                 0.9999077

    rldpvalslong <- rldpvals
    rldpvalslong$gene <- row.names(rldpvalslong) 
    rldpvalslong <- melt(rldpvalslong, id=c("gene"))
    head(rldpvalslong)

    ##            gene       variable     value
    ## 1 0610007P14Rik pvalPunchCA1DG 0.4379056
    ## 2 0610009B22Rik pvalPunchCA1DG 0.5181528
    ## 3 0610009L18Rik pvalPunchCA1DG 0.9850306
    ## 4 0610009O20Rik pvalPunchCA1DG 0.8413827
    ## 5 0610010F05Rik pvalPunchCA1DG 0.3899734
    ## 6 0610010K14Rik pvalPunchCA1DG 0.5136582

    qplot(value, data=rldpvalslong, geom="histogram") + 
      facet_grid( ~ variable) +
      scale_y_log10()

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 128 rows containing non-finite values (stat_bin).

    ## Warning: Stacking not well defined when ymin != 0

![](../figures/03_dissociationstresstest/pvaluedistribution-1.png)

I'm not really happy with these two heat maps. Here's how I created
them. Top heatmap: subset the data to give only the gene with an
adjusted p value &lt; 0.05 for the homogenized vs dissociated
comparisonany two-way comparsion. Bottom heatmap: subset the data to
give only the gene with an adjusted p value &lt; 0.05 for two way brain
region comparision (CA1 vs DG, CA3, vs DG, or CA1 vs DG)

Here, you can see that the differences between samples is not as clear
cut for all comparisions. What other mechanisms would be useful for
subseting the data to identify genes of interest?

![](../figures/03_dissociationstresstest/Heatmap100DEgenes-1.png)

This is a data validation check plot. Here, I'm showing how many
millions of reads were present in each sample. On average, each sample
had 5 million reads, but the range was from 0.8 to 10 millino reads.

    FALSE [1] 22485    32

    FALSE  100-CA1-1  100-CA1-2  100-CA1-3  100-CA3-1  100-CA3-4   100-DG-2 
    FALSE   1.136597   3.311998   1.114747   0.966391   1.205348   0.658410 
    FALSE   100-DG-3  101-CA1-1  101-CA1-2  101-CA1-3  101-CA3-1  101-CA3-4 
    FALSE   3.055740   2.668415   0.072040   0.154427   1.361076   0.639942 
    FALSE   101-DG-3   101-DG-4 143B-CA1-1  143B-DG-1 144B-CA1-1 144B-CA3-1 
    FALSE   0.036498   0.300618   0.874614   1.019113   1.275137   0.506698 
    FALSE 145B-CA1-1  145B-DG-1 146B-CA1-2 146B-CA3-2  146B-DG-2  147-CA1-4 
    FALSE   1.034066   0.720798   0.506014   1.056001   0.055549   0.080721 
    FALSE  147-CA3-4   147-DG-4  148-CA1-2  148-CA3-2   148-DG-2 148B-CA1-4 
    FALSE   0.344588   0.069648   0.938866   1.148136   1.067185   0.185637 
    FALSE 148B-CA3-4  148B-DG-4 
    FALSE   1.724144   0.398258

    FALSE 
    FALSE    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
    FALSE 4907  391  378  280  231  208  176  153  117  127  123   99  116   95  102 
    FALSE   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29 
    FALSE   83   73   90   82  102   57   76   68   78   60   68   70   50   47   49

These next graphs show the correlation between samples of CA1, CA3, and
DG.

![](../figures/03_dissociationstresstest/scatterplots-1.png)![](../figures/03_dissociationstresstest/scatterplots-2.png)

This next plot shows the stregnth of the correlation between samples.

![](../figures/03_dissociationstresstest/correlationmatrix-1.png)

Save files for GO analysis. A total of 217 DEGs with unadjusted p-value
&lt; 0.1 were input into the GO anlaysis.
