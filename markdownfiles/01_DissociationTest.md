Cellular Stress
---------------

In this analysis, I examine the effect that cell dissasociation has on
CA1, CA3, and DG gene expression relative to control tissue samples.

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
    ##  $ Treatment: Factor w/ 2 levels "control","dissociated": 1 1 1 1 1 1 1 2 2 2 ...
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
    ##       Group    Conflict    APA           Treatment     dodgy   
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
CA1 and CA3 samples have similar transcriptomes. The control CA1 samples
have the most similar transcriptonal profiles as evidenced by their
tight clustering.

    source("DESeqPCAfunction.R")
    source("figureoptions.R")

    # create the dataframe using my function pcadataframe
    pcadata <- pcadataframe(rld, intgroup=c("Region", "Treatment"), returnData=TRUE)
    percentVar <- round(100 * attr(pcadata, "percentVar"))

    ## for markdown
    plotPC2PC1(aescolor = pcadata$Region, colorname = "Region", colorvalues = colorvalRegion, aesshape = pcadata$Treatment, shapename = "Treatment")

![](../figures/01_dissociationtest/PCA-1.png)

    # for adobe
    myplot <- plotPC2PC1(aescolor = pcadata$Region, colorname = "Region", aesshape = pcadata$Treatment, shapename = "Treatment", colorvalues = colorvalRegion)
    pdf(file="../figures/01_dissociationtest/PCA-1.pdf", width=4.5, height=3)
    plot(myplot)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    ## statistics
    aov1 <- aov(PC1 ~ Region, data=pcadata)
    summary(aov1) 

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Region       2 2812.7  1406.4   17.69 0.000365 ***
    ## Residuals   11  874.3    79.5                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov1, which = "Region") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ Region, data = pcadata)
    ## 
    ## $Region
    ##              diff       lwr      upr     p adj
    ## CA3-CA1  5.223942 -10.31899 20.76687 0.6467956
    ## DG-CA1  33.098083  17.55515 48.64101 0.0003454
    ## DG-CA3  27.874142  10.84772 44.90057 0.0027013

    aov2 <- aov(PC2 ~ Region, data=pcadata)
    summary(aov2) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Region       2  243.8   121.9   0.744  0.498
    ## Residuals   11 1801.4   163.8

    TukeyHSD(aov2, which = "Region") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ Region, data = pcadata)
    ## 
    ## $Region
    ##              diff       lwr      upr     p adj
    ## CA3-CA1 -8.297717 -30.60810 14.01267 0.5893113
    ## DG-CA1   1.924204 -20.38618 24.23459 0.9706097
    ## DG-CA3  10.221920 -14.21788 34.66172 0.5166917

    aov3 <- aov(PC1 ~ Treatment, data=pcadata)
    summary(aov3) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment    1    335   335.2     1.2  0.295
    ## Residuals   12   3352   279.3

    TukeyHSD(aov3, which = "Treatment")

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ Treatment, data = pcadata)
    ## 
    ## $Treatment
    ##                         diff       lwr      upr     p adj
    ## dissociated-control 9.785654 -9.678654 29.24996 0.2948417

    aov4 <- aov(PC2 ~ Treatment, data=pcadata)
    summary(aov4) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## Treatment    1  691.1   691.1   6.125 0.0292 *
    ## Residuals   12 1354.1   112.8                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov4, which = "Treatment") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ Treatment, data = pcadata)
    ## 
    ## $Treatment
    ##                          diff       lwr       upr     p adj
    ## dissociated-control -14.05242 -26.42372 -1.681116 0.0292306

    lm1 <- lm(PC1~Region*Treatment, data=pcadata)
    summary(lm1)

    ## 
    ## Call:
    ## lm(formula = PC1 ~ Region * Treatment, data = pcadata)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -9.9926 -5.2517  0.4166  5.2517  9.6799 
    ## 
    ## Coefficients:
    ##                                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                     -17.672      4.566  -3.870 0.004740 ** 
    ## RegionCA3                         7.758      7.220   1.075 0.313934    
    ## RegionDG                         36.970      7.220   5.121 0.000907 ***
    ## Treatmentdissociated             13.446      6.458   2.082 0.070874 .  
    ## RegionCA3:Treatmentdissociated   -5.068     10.210  -0.496 0.633022    
    ## RegionDG:Treatmentdissociated    -7.743     10.210  -0.758 0.469970    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 7.909 on 8 degrees of freedom
    ## Multiple R-squared:  0.8643, Adjusted R-squared:  0.7795 
    ## F-statistic: 10.19 on 5 and 8 DF,  p-value: 0.002577

    anova(lm1) 

    ## Analysis of Variance Table
    ## 
    ## Response: PC1
    ##                  Df  Sum Sq Mean Sq F value    Pr(>F)    
    ## Region            2 2812.72 1406.36 22.4834 0.0005204 ***
    ## Treatment         1  335.16  335.16  5.3581 0.0493207 *  
    ## Region:Treatment  2   38.75   19.37  0.3097 0.7420561    
    ## Residuals         8  500.41   62.55                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    lm2 <- lm(PC2~Region*Treatment, data=pcadata)
    summary(lm2)

    ## 
    ## Call:
    ## lm(formula = PC2 ~ Region * Treatment, data = pcadata)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -17.7457  -4.4746   0.0085   4.4746  17.7457 
    ## 
    ## Coefficients:
    ##                                Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                      9.9925     6.7495   1.480    0.177
    ## RegionCA3                      -10.9299    10.6719  -1.024    0.336
    ## RegionDG                         0.5478    10.6719   0.051    0.960
    ## Treatmentdissociated           -16.3430     9.5453  -1.712    0.125
    ## RegionCA3:Treatmentdissociated   5.2643    15.0924   0.349    0.736
    ## RegionDG:Treatmentdissociated    2.7528    15.0924   0.182    0.860
    ## 
    ## Residual standard error: 11.69 on 8 degrees of freedom
    ## Multiple R-squared:  0.4654, Adjusted R-squared:  0.1313 
    ## F-statistic: 1.393 on 5 and 8 DF,  p-value: 0.3218

    anova(lm2) 

    ## Analysis of Variance Table
    ## 
    ## Response: PC2
    ##                  Df  Sum Sq Mean Sq F value  Pr(>F)  
    ## Region            2  243.79  121.90  0.8919 0.44702  
    ## Treatment         1  691.15  691.15  5.0571 0.05467 .
    ## Region:Treatment  2   16.93    8.46  0.0619 0.94040  
    ## Residuals         8 1093.35  136.67                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Now, we can calulate the number of significant genes by contrast by
contrast. The first number displayed is not corrected for mutiple
hypothesis testing but the second one is.

    ## DEG by contrasts
    source("resvalsfunction.R")
    contrast1 <- resvals(contrastvector = c('Region', 'CA1', 'DG'), mypval = 0.1)

    ## [1] 1959
    ## [1] 404

    contrast2 <- resvals(contrastvector = c('Region', 'CA3', 'DG'), mypval = 0.1)

    ## [1] 1042
    ## [1] 82

    contrast3 <- resvals(contrastvector = c('Region', 'CA1', 'CA3'), mypval = 0.1)

    ## [1] 693
    ## [1] 9

    contrast4 <- resvals(contrastvector = c('Treatment', 'dissociated', 'control'), mypval = 0.1)

    ## [1] 2349
    ## [1] 260

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
Region and Treatment. This shows all genes with *adjusted* pvalue
&lt;0.1.

![](../figures/01_dissociationtest/VennDiagramPadj-1.png)

heatmaps
========

![](../figures/01_dissociationtest/HeatmapPadj-1.png)

This is a data validation check plot. Here, I'm showing how many
millions of reads were present in each sample. On average, each sample
had 5 million reads, but the range was from 0.8 to 10 millino reads.

    FALSE [1] 22485    14

    FALSE 100-CA1-1 100-CA1-2 100-CA1-3 100-CA3-1 100-CA3-4  100-DG-2  100-DG-3 
    FALSE  2.311086  6.646655  2.277596  1.974845  2.352153  1.285654  6.086605 
    FALSE 101-CA1-1 101-CA1-2 101-CA1-3 101-CA3-1 101-CA3-4  101-DG-3  101-DG-4 
    FALSE  4.782767  0.135065  0.300812  2.498914  1.193153  0.065887  0.598775

    FALSE 
    FALSE    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
    FALSE 5099  373  304  256  193  189  161  132  137  113  131  110  107   85   75 
    FALSE   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29 
    FALSE   83   82   57   77   67   61   69   48   52   63   60   58   61   50   61

Save files for GO analysis. A total of 217 DEGs with unadjusted p-value
&lt; 0.1 were input into the GO anlaysis.

    FALSE 
    FALSE FALSE  TRUE 
    FALSE 11813   344

    FALSE 
    FALSE FALSE  TRUE 
    FALSE 15163  1528

    FALSE log2 fold change (MLE): Treatment dissociated vs control 
    FALSE Wald test p-value: Treatment dissociated vs control 
    FALSE DataFrame with 6 rows and 6 columns
    FALSE                baseMean log2FoldChange     lfcSE       stat    pvalue
    FALSE               <numeric>      <numeric> <numeric>  <numeric> <numeric>
    FALSE 0610007P14Rik 30.740657     -0.5281136 0.7400694 -0.7136001 0.4754745
    FALSE 0610009B22Rik 14.127818      0.2891053 0.9835019  0.2939550 0.7687923
    FALSE 0610009L18Rik  7.310401     -0.3916629 1.3093270 -0.2991330 0.7648386
    FALSE 0610009O20Rik 42.651762      0.2571289 0.5017330  0.5124815 0.6083141
    FALSE 0610010F05Rik 41.148916     -0.2712032 0.5042840 -0.5377985 0.5907161
    FALSE 0610010K14Rik 13.717020      0.2983853 0.7673449  0.3888542 0.6973840
    FALSE                    padj
    FALSE               <numeric>
    FALSE 0610007P14Rik 0.8680052
    FALSE 0610009B22Rik 0.9591621
    FALSE 0610009L18Rik 0.9574803
    FALSE 0610009O20Rik 0.9172761
    FALSE 0610010F05Rik 0.9122632
    FALSE 0610010K14Rik 0.9401399

    FALSE sign
    FALSE   -1    1 
    FALSE 6989 9720

    FALSE Bootstrap (r = 0.5)... Done.
    FALSE Bootstrap (r = 0.6)... Done.
    FALSE Bootstrap (r = 0.7)... Done.
    FALSE Bootstrap (r = 0.8)... Done.
    FALSE Bootstrap (r = 0.9)... Done.
    FALSE Bootstrap (r = 1.0)... Done.
    FALSE Bootstrap (r = 1.1)... Done.
    FALSE Bootstrap (r = 1.2)... Done.
    FALSE Bootstrap (r = 1.3)... Done.
    FALSE Bootstrap (r = 1.4)... Done.

![](../figures/01_dissociationtest/pvclust-1.png)
