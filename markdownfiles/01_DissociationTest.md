Methods for Dorsal Hippocampal Gene Expression Profiling
--------------------------------------------------------

#### Part 1: Examining the influence of dissasociation on gene expression in the CA1, CA3, and DG

Subset to just look homogenized and dissociated samples
-------------------------------------------------------

    colData <- colData %>%
      filter(Mouse %in% c("15-100")) %>% droplevels()
    savecols <- as.character(colData$RNAseqID) #selects all good samples
    savecols <- as.vector(savecols) # make it a vector
    countData <- countData %>% select(one_of(savecols)) # keep good samples

![](../figures/01_dissociationtest/DEG-1.png)

![](../figures/01_dissociationtest/MA-1.png)![](../figures/01_dissociationtest/MA-2.png)

This PCA gives an overview of the variability between samples using the
a large matrix of log transformed gene expression data. You can see that
the bigges difference is between DG punches and the CA1 and CA3 punches.
CA1 and CA3 samples have similar transcriptomes. The homogenized CA1
samples have the most similar transcriptonal profiles as evidenced by
their tight clustering.

![](../figures/01_dissociationtest/PCA-1.png)![](../figures/01_dissociationtest/PCA-2.png)

This Venn Diagram shows the number of differentially expressed by
contrast described above each oval. The most number of genes are
differntially expressed between DG and the CAs (nearly 1000) wheras only
about 200 were differntailly regulated as a result of of technical
maniplulation comparing homogenized and dissociated samples.

The first is with padj values. The second with p values

![](../figures/01_dissociationtest/VennDiagram1-1.png)

![](../figures/01_dissociationtest/VennDiagram2-1.png)

Here, the goal is the analyze the distribution of pvalues to see if they
are randomly distributed or if that is a tendency towards and increase
or decrease of low pvalues. There, I'm showing the pval and adjusted
pvale (padj) for all for two-way comparision.

    head(rldpvals)

    ##               pvalPunchCA1DG padjPunchCA1DG pvalPunchCA3DG padjPunchCA3DG
    ## 0610007P14Rik     0.56823853      1.0000000     0.74481477      0.9823162
    ## 0610009B22Rik     0.64786125      1.0000000     0.45021525      0.9288071
    ## 0610009L18Rik     0.01596479      0.3019016     0.04667271      0.5130170
    ## 0610009O20Rik     0.97898532      1.0000000     0.74979209      0.9836092
    ## 0610010F05Rik     0.61887935      1.0000000     0.56833596      0.9294782
    ## 0610010K14Rik     0.50933468      1.0000000     0.05782676      0.5562383
    ##               pvalPunchCA1CA3 padjPunchCA1CA3
    ## 0610007P14Rik       0.2833914               1
    ## 0610009B22Rik       0.7017419               1
    ## 0610009L18Rik       0.6256568               1
    ## 0610009O20Rik       0.7142028               1
    ## 0610010F05Rik       0.9220665               1
    ## 0610010K14Rik       0.1486267               1
    ##               pvalmethodhomogenizeddissociated
    ## 0610007P14Rik                        0.5375223
    ## 0610009B22Rik                        0.8919044
    ## 0610009L18Rik                        0.7574899
    ## 0610009O20Rik                        0.5562629
    ## 0610010F05Rik                        0.6300083
    ## 0610010K14Rik                        0.8105541
    ##               padjmethodhomogenizeddissociated
    ## 0610007P14Rik                         0.999948
    ## 0610009B22Rik                         0.999948
    ## 0610009L18Rik                         0.999948
    ## 0610009O20Rik                         0.999948
    ## 0610010F05Rik                         0.999948
    ## 0610010K14Rik                         0.999948

    rldpvalslong <- rldpvals
    rldpvalslong$gene <- row.names(rldpvalslong) 
    rldpvalslong <- melt(rldpvalslong, id=c("gene"))
    head(rldpvalslong)

    ##            gene       variable      value
    ## 1 0610007P14Rik pvalPunchCA1DG 0.56823853
    ## 2 0610009B22Rik pvalPunchCA1DG 0.64786125
    ## 3 0610009L18Rik pvalPunchCA1DG 0.01596479
    ## 4 0610009O20Rik pvalPunchCA1DG 0.97898532
    ## 5 0610010F05Rik pvalPunchCA1DG 0.61887935
    ## 6 0610010K14Rik pvalPunchCA1DG 0.50933468

    qplot(value, data=rldpvalslong, geom="histogram") + 
      facet_grid( ~ variable) +
      scale_y_log10()

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 72 rows containing non-finite values (stat_bin).

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 15 rows containing missing values (geom_bar).

![](../figures/01_dissociationtest/pvaluedistribution-1.png)

I'm not really happy with these two heat maps. Here's how I created
them. Top heatmap: subset the data to give only the gene with an
adjusted p value &lt; 0.05 for the homogenized vs dissociated
comparisonany two-way comparsion. Bottom heatmap: subset the data to
give only the gene with an adjusted p value &lt; 0.05 for two way brain
region comparision (CA1 vs DG, CA3, vs DG, or CA1 vs DG)

Here, you can see that the differences between samples is not as clear
cut for all comparisions. What other mechanisms would be useful for
subseting the data to identify genes of interest?

![](../figures/01_dissociationtest/Heatmap100DEgenes-1.png)

This is a data validation check plot. Here, I'm showing how many
millions of reads were present in each sample. On average, each sample
had 5 million reads, but the range was from 0.8 to 10 millino reads.

    FALSE [1] 22485    14

    FALSE 100-CA1-1 100-CA1-2 100-CA1-3 100-CA3-1 100-CA3-4  100-DG-2  100-DG-3 
    FALSE  1.136597  3.311998  1.114747  0.966391  1.205348  0.658410  3.055740 
    FALSE 101-CA1-1 101-CA1-2 101-CA1-3 101-CA3-1 101-CA3-4  101-DG-3  101-DG-4 
    FALSE  2.668415  0.072040  0.154427  1.361076  0.639942  0.036498  0.300618

    FALSE 
    FALSE    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
    FALSE 5353  455  405  314  258  226  229  158  160  147  149  131  110  121  116 
    FALSE   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29 
    FALSE   96   91   81   84   80   82   82   73   67   68   68   62   56   60   68

These next graphs show the correlation between samples of CA1, CA3, and
DG.

![](../figures/01_dissociationtest/scatterplots-1.png)![](../figures/01_dissociationtest/scatterplots-2.png)

This next plot shows the stregnth of the correlation between samples.

![](../figures/01_dissociationtest/correlationmatrix-1.png)

Save files for GO analysis. A total of 217 DEGs with unadjusted p-value
&lt; 0.1 were input into the GO anlaysis.

    rld <- rlogTransformation(dds)
    head(assay(rld))

    ##               100-CA1-1 100-CA1-2 100-CA1-3 100-CA3-1 100-CA3-4 100-DG-2
    ## 0610007P14Rik 3.5995351  3.790404 3.8704858 4.1208806 4.1989841 4.255601
    ## 0610009B22Rik 2.1590559  2.707684 2.0431224 2.6448498 2.0159499 2.651163
    ## 0610009L18Rik 1.4931716  2.425456 1.2970543 2.4611627 2.0847477 1.594549
    ## 0610009O20Rik 4.7233740  5.314397 5.1705766 5.6247761 5.4613965 5.641163
    ## 0610010F05Rik 2.3800261  2.479708 2.4360218 2.2917875 2.4704070 2.488071
    ## 0610010K14Rik 0.3022043  0.332851 0.2985802 0.1973599 0.3903414 0.269161
    ##                100-DG-3 101-CA1-1 101-CA1-2 101-CA1-3  101-CA3-1
    ## 0610007P14Rik 4.0935099 3.9721358 2.9341147 2.5612903 4.13618006
    ## 0610009B22Rik 2.4832752 2.2845668 3.1979993 1.7277362 2.85895766
    ## 0610009L18Rik 2.6263074 1.9406482 1.7805353 1.5203779 2.25283859
    ## 0610009O20Rik 5.2867893 5.5066598 4.9707980 5.1562400 5.60263405
    ## 0610010F05Rik 2.3849418 2.1220065 2.6067576 2.2536438 2.01566965
    ## 0610010K14Rik 0.4326066 0.4734489 0.3929896 0.2469774 0.04369077
    ##                101-CA3-4  101-DG-3  101-DG-4
    ## 0610007P14Rik 4.07201279 4.6041895 3.1841046
    ## 0610009B22Rik 2.46109194 3.0350883 1.4928936
    ## 0610009L18Rik 1.73221613 3.5545058 3.2983879
    ## 0610009O20Rik 5.28635863 6.0349705 4.9438768
    ## 0610010F05Rik 2.54003767 2.2344069 1.9424516
    ## 0610010K14Rik 0.07695913 0.5530437 0.7417628

    colnames(rld) <- paste(colData$Punch, colData$method, colData$RNAseqID, sep = "")
    head(assay(rld))

    ##               CA1homogenized100-CA1-1 CA1homogenized100-CA1-2
    ## 0610007P14Rik               3.5995351                3.790404
    ## 0610009B22Rik               2.1590559                2.707684
    ## 0610009L18Rik               1.4931716                2.425456
    ## 0610009O20Rik               4.7233740                5.314397
    ## 0610010F05Rik               2.3800261                2.479708
    ## 0610010K14Rik               0.3022043                0.332851
    ##               CA1homogenized100-CA1-3 CA3homogenized100-CA3-1
    ## 0610007P14Rik               3.8704858               4.1208806
    ## 0610009B22Rik               2.0431224               2.6448498
    ## 0610009L18Rik               1.2970543               2.4611627
    ## 0610009O20Rik               5.1705766               5.6247761
    ## 0610010F05Rik               2.4360218               2.2917875
    ## 0610010K14Rik               0.2985802               0.1973599
    ##               CA3homogenized100-CA3-4 DGhomogenized100-DG-2
    ## 0610007P14Rik               4.1989841              4.255601
    ## 0610009B22Rik               2.0159499              2.651163
    ## 0610009L18Rik               2.0847477              1.594549
    ## 0610009O20Rik               5.4613965              5.641163
    ## 0610010F05Rik               2.4704070              2.488071
    ## 0610010K14Rik               0.3903414              0.269161
    ##               DGhomogenized100-DG-3 CA1dissociated101-CA1-1
    ## 0610007P14Rik             4.0935099               3.9721358
    ## 0610009B22Rik             2.4832752               2.2845668
    ## 0610009L18Rik             2.6263074               1.9406482
    ## 0610009O20Rik             5.2867893               5.5066598
    ## 0610010F05Rik             2.3849418               2.1220065
    ## 0610010K14Rik             0.4326066               0.4734489
    ##               CA1dissociated101-CA1-2 CA1dissociated101-CA1-3
    ## 0610007P14Rik               2.9341147               2.5612903
    ## 0610009B22Rik               3.1979993               1.7277362
    ## 0610009L18Rik               1.7805353               1.5203779
    ## 0610009O20Rik               4.9707980               5.1562400
    ## 0610010F05Rik               2.6067576               2.2536438
    ## 0610010K14Rik               0.3929896               0.2469774
    ##               CA3dissociated101-CA3-1 CA3dissociated101-CA3-4
    ## 0610007P14Rik              4.13618006              4.07201279
    ## 0610009B22Rik              2.85895766              2.46109194
    ## 0610009L18Rik              2.25283859              1.73221613
    ## 0610009O20Rik              5.60263405              5.28635863
    ## 0610010F05Rik              2.01566965              2.54003767
    ## 0610010K14Rik              0.04369077              0.07695913
    ##               DGdissociated101-DG-3 DGdissociated101-DG-4
    ## 0610007P14Rik             4.6041895             3.1841046
    ## 0610009B22Rik             3.0350883             1.4928936
    ## 0610009L18Rik             3.5545058             3.2983879
    ## 0610009O20Rik             6.0349705             4.9438768
    ## 0610010F05Rik             2.2344069             1.9424516
    ## 0610010K14Rik             0.5530437             0.7417628

    pheatmap(cor(assay(rld)),border_color=NA, main="SampleHeatmap")

![](../figures/01_dissociationtest/sampleheatmap-1.png)
