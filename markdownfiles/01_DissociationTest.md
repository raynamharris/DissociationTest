Examining Factors that influence dorsal hippocampal gene expression profiling
-----------------------------------------------------------------------------

In this analysis, I examine the effect that cells dissasociation has on
CA1, CA3, and DG gene expression relative to homogenized tissue samples.

Here is a brief overview of the samples being compared.

    str(colData)

    ## 'data.frame':    14 obs. of  14 variables:
    ##  $ RNAseqID : Factor w/ 14 levels "100-CA1-1","100-CA1-2",..: 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ Mouse    : Factor w/ 1 level "15-100": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ year     : int  2015 2015 2015 2015 2015 2015 2015 2015 2015 2015 ...
    ##  $ Genotype : Factor w/ 1 level "WT": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ jobnumber: Factor w/ 1 level "JA16444": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Region   : Factor w/ 3 levels "CA1","CA3","DG": 1 1 1 2 2 3 3 1 1 1 ...
    ##  $ Group    : Factor w/ 1 level "homecage": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Conflict : Factor w/ 0 levels: NA NA NA NA NA NA NA NA NA NA ...
    ##  $ APA      : Factor w/ 0 levels: NA NA NA NA NA NA NA NA NA NA ...
    ##  $ Method   : Factor w/ 2 levels "control","dissociated": 1 1 1 1 1 1 1 2 2 2 ...
    ##  $ dodgy    : Factor w/ 1 level "allgood": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ daytime  : Factor w/ 1 level "norecord": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Slice    : int  1 2 3 1 4 2 3 1 2 3 ...
    ##  $ Date     : Factor w/ 1 level "9/28/15": 1 1 1 1 1 1 1 1 1 1 ...

    summary(colData)

    ##       RNAseqID    Mouse         year      Genotype   jobnumber  Region 
    ##  100-CA1-1:1   15-100:14   Min.   :2015   WT:14    JA16444:14   CA1:6  
    ##  100-CA1-2:1               1st Qu.:2015                         CA3:4  
    ##  100-CA1-3:1               Median :2015                         DG :4  
    ##  100-CA3-1:1               Mean   :2015                                
    ##  100-CA3-4:1               3rd Qu.:2015                                
    ##  100-DG-2 :1               Max.   :2015                                
    ##  (Other)  :8                                                           
    ##       Group    Conflict    APA             Method      dodgy   
    ##  homecage:14   NA's:14   NA's:14   control    :7   allgood:14  
    ##                                    dissociated:7               
    ##                                                                
    ##                                                                
    ##                                                                
    ##                                                                
    ##                                                                
    ##      daytime       Slice            Date   
    ##  norecord:14   Min.   :1.000   9/28/15:14  
    ##                1st Qu.:1.250               
    ##                Median :2.500               
    ##                Mean   :2.429               
    ##                3rd Qu.:3.000               
    ##                Max.   :4.000               
    ## 

This PCA gives an overview of the variability between samples using the
a large matrix of log transformed gene expression data. You can see that
the bigges difference is between DG punches and the CA1 and CA3 punches.
CA1 and CA3 samples have similar transcriptomes. The homogenized CA1
samples have the most similar transcriptonal profiles as evidenced by
their tight clustering.

![](../figures/01_dissociationtest/PCA-1.png)

Now, we can calulate the number of significant genes by contrast by
contrast. The first number displayed is not corrected for mutiple
hypothesis testing but the second one is.

    ## DEG by contrasts
    source("resvalsfunction.R")
    contrast1 <- resvals(contrastvector = c('Region', 'CA1', 'DG'), mypval = 0.1)

    ## [1] 1723
    ## [1] 281

    contrast2 <- resvals(contrastvector = c('Region', 'CA3', 'DG'), mypval = 0.1)

    ## [1] 903
    ## [1] 46

    contrast3 <- resvals(contrastvector = c('Region', 'CA1', 'CA3'), mypval = 0.1)

    ## [1] 551
    ## [1] 5

    contrast4 <- resvals(contrastvector = c('Method', 'control', 'dissociated'), mypval = 0.1)

    ## [1] 1662
    ## [1] 129

Now, we can view a histogram of the distribution

![](../figures/01_dissociationtest/histogram-1.png)

    ## [1] 1

![](../figures/01_dissociationtest/histogram-2.png)

    ## [1] 1

![](../figures/01_dissociationtest/histogram-3.png)

    ## [1] 1

![](../figures/01_dissociationtest/histogram-4.png)

    ## [1] 1

This Venn Diagram sthe overlap of differentailly expression genes by
Region and method. This shows all genes with *uncorrected* pvalue
&lt;0.1.

This Venn Diagram sthe overlap of differentailly expression genes by
Region and method. This shows all genes with *adjusted* pvalue &lt;0.1.

    ## null device 
    ##           1

I'm not really happy with these two heat maps. Here's how I created
them. Top heatmap: subset the data to give only the gene with an
adjusted p value &lt; 0.05 for the homogenized vs dissociated
comparisonany two-way comparsion. Bottom heatmap: subset the data to
give only the gene with an adjusted p value &lt; 0.05 for two way brain
region comparision (CA1 vs DG, CA3, vs DG, or CA1 vs DG)

Here, you can see that the differences between samples is not as clear
cut for all comparisions. What other mechanisms would be useful for
subseting the data to identify genes of interest?

![](../figures/01_dissociationtest/HeatmapPadj-1.png)

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

    colnames(rld) <- paste(colData$Region, colData$Method, colData$RNAseqID, sep = "")
    head(assay(rld))

    ##               CA1control100-CA1-1 CA1control100-CA1-2 CA1control100-CA1-3
    ## 0610007P14Rik           3.5995351            3.790404           3.8704858
    ## 0610009B22Rik           2.1590559            2.707684           2.0431224
    ## 0610009L18Rik           1.4931716            2.425456           1.2970543
    ## 0610009O20Rik           4.7233740            5.314397           5.1705766
    ## 0610010F05Rik           2.3800261            2.479708           2.4360218
    ## 0610010K14Rik           0.3022043            0.332851           0.2985802
    ##               CA3control100-CA3-1 CA3control100-CA3-4 DGcontrol100-DG-2
    ## 0610007P14Rik           4.1208806           4.1989841          4.255601
    ## 0610009B22Rik           2.6448498           2.0159499          2.651163
    ## 0610009L18Rik           2.4611627           2.0847477          1.594549
    ## 0610009O20Rik           5.6247761           5.4613965          5.641163
    ## 0610010F05Rik           2.2917875           2.4704070          2.488071
    ## 0610010K14Rik           0.1973599           0.3903414          0.269161
    ##               DGcontrol100-DG-3 CA1dissociated101-CA1-1
    ## 0610007P14Rik         4.0935099               3.9721358
    ## 0610009B22Rik         2.4832752               2.2845668
    ## 0610009L18Rik         2.6263074               1.9406482
    ## 0610009O20Rik         5.2867893               5.5066598
    ## 0610010F05Rik         2.3849418               2.1220065
    ## 0610010K14Rik         0.4326066               0.4734489
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
