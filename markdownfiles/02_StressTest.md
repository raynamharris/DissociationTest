Behavioral Stress
-----------------

In this analysis, I examine the effect that behavioral stress has on
CA1, CA3, and DG gene expression relative to homogenized tissue samples.

    #source("http://www.bioconductor.org/biocLite.R")
    #biocLite("DESeq2")
    library(DESeq2)
    library(magrittr)
    library(tidyverse)
    library(reshape2)
    library(VennDiagram)
    library(genefilter)
    library(pheatmap)
    library(cowplot)
    library(RColorBrewer)
    library(ggcorrplot)
    library(dplyr)
    library(plyr)
    library(ggplot2)
    library(colorRamps)

    # set output file for figures 
    knitr::opts_chunk$set(fig.path = '../figures/02_stresstest/')

    # this starts with data genearated from code described in KallistoGather.Rmd
    colData <- read.csv('../data/DissociationColData.csv')
    rownames(colData) <- colData$RNAseqID
    countData <-  read.csv('../data/DissociationCountData.csv', check.names = F, row.names = 1)

Subset to just look homogenized and dissociated samples
-------------------------------------------------------

    colData <- colData %>%
      filter(Mouse != "15-100") %>% droplevels()
    savecols <- as.character(colData$RNAseqID) #selects all good samples
    savecols <- as.vector(savecols) # make it a vector
    countData <- countData %>% select(one_of(savecols)) # keep good samples

    ## rename and relevel things
    colData <- rename(colData, c("Group"="Treatment"))
    colData$Treatment <- plyr::revalue(colData$Treatment, c("control"="shocked"))
    colData$Treatment <- factor(colData$Treatment, levels = c("homecage", "shocked"))

    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ Treatment + Region + Treatment * Region )
    dds <- dds[ rowSums(counts(dds)) > 2, ] ## filter genes with 0 counts
    dds <- DESeq(dds) # Differential expression analysis
    dds

    ## class: DESeqDataSet 
    ## dim: 15722 18 
    ## metadata(1): version
    ## assays(3): counts mu cooks
    ## rownames(15722): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(37): baseMean baseVar ... deviance maxCooks
    ## colnames(18): 143B-CA1-1 143B-DG-1 ... 148B-CA3-4 148B-DG-4
    ## colData names(15): RNAseqID Mouse ... Date sizeFactor

    ## for variance stablized gene expression and log transformed data
    rld <- rlog(dds, blind=FALSE)
    vsd <- getVarianceStabilizedData(dds)

This PCA gives an overview of the variability between samples using the
a large matrix of log transformed gene expression data. You can see that
the bigges difference is between DG punches and the CA1 and CA3 punches.
CA1 and CA3 samples have similar transcriptomes. The homogenized CA1
samples have the most similar transcriptonal profiles as evidenced by
their tight clustering.

    source("DESeqPCAfunction.R")
    source("figureoptions.R")

    # create the dataframe using my function pcadataframe
    pcadata <- pcadataframe(rld, intgroup=c("Treatment", "Region"), returnData=TRUE)
    percentVar <- round(100 * attr(pcadata, "percentVar"))

    pcadata$Treatment <- factor(pcadata$Treatment, levels = c("homecage", "shocked"))

    ## plot a bunch of pca plots using my ggplot functions DESeqPCAfunction.R, with the color defined infigureoptions.R
    plotPC1PC2(aescolor = pcadata$Region, colorname = "Region", aesshape = pcadata$Treatment, shapename = "Treatment", colorvalues = colorvalRegion)

![](../figures/02_stresstest/PCA-1.png)

    plotPC2PC3(aescolor = pcadata$Region, colorname = "Region", aesshape = pcadata$Treatment, shapename = "Treatment", colorvalues = colorvalRegion)

![](../figures/02_stresstest/PCA-2.png)

    # for adobe
    myplot <- plotPC1PC2(aescolor = pcadata$Region, colorname = "Region", aesshape = pcadata$Treatment, shapename = "Treatment", colorvalues = colorvalRegion)
    pdf(file="../figures/../figures/02_stresstest/PCA-1.pdf", width=4.5, height=3)
    plot(myplot)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    ## statistics
    aov1 <- aov(PC1 ~ Region, data=pcadata)
    summary(aov1) 

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Region       2   6559    3279   93.99 3.27e-09 ***
    ## Residuals   15    523      35                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov1, which = "Region") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ Region, data = pcadata)
    ## 
    ## $Region
    ##              diff      lwr      upr     p adj
    ## CA3-CA1  4.481862 -4.50180 13.46552 0.4189128
    ## DG-CA1  42.179395 33.64360 50.71519 0.0000000
    ## DG-CA3  37.697533 28.40717 46.98789 0.0000001

    aov2 <- aov(PC2 ~ Region, data=pcadata)
    summary(aov2) 

    ##             Df Sum Sq Mean Sq F value  Pr(>F)    
    ## Region       2   3761  1880.3   48.59 2.8e-07 ***
    ## Residuals   15    581    38.7                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov2, which = "Region") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ Region, data = pcadata)
    ## 
    ## $Region
    ##              diff       lwr       upr     p adj
    ## CA3-CA1  35.82874  26.36699  45.29049 0.0000002
    ## DG-CA1   12.89803   3.90798  21.88807 0.0054054
    ## DG-CA3  -22.93071 -32.71548 -13.14595 0.0000584

    aov5 <- aov(PC3 ~ Region, data=pcadata)
    summary(aov5) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Region       2  364.3   182.2   0.993  0.393
    ## Residuals   15 2751.2   183.4

    TukeyHSD(aov5, which = "Region") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC3 ~ Region, data = pcadata)
    ## 
    ## $Region
    ##               diff       lwr      upr     p adj
    ## CA3-CA1 -10.973777 -31.57165  9.62410 0.3737501
    ## DG-CA1   -2.762631 -22.33363 16.80836 0.9289204
    ## DG-CA3    8.211147 -13.08993 29.51222 0.5873230

    aov3 <- aov(PC1 ~ Treatment, data=pcadata)
    summary(aov3) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment    1     50    50.4   0.115  0.739
    ## Residuals   16   7032   439.5

    TukeyHSD(aov3, which = "Treatment")

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ Treatment, data = pcadata)
    ## 
    ## $Treatment
    ##                       diff       lwr      upr     p adj
    ## shocked-homecage -3.549567 -25.77028 18.67114 0.7392868

    aov4 <- aov(PC2 ~ Treatment, data=pcadata)
    summary(aov4) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment    1      5    5.19   0.019  0.892
    ## Residuals   16   4336  271.00

    TukeyHSD(aov4, which = "Treatment") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ Treatment, data = pcadata)
    ## 
    ## $Treatment
    ##                      diff       lwr      upr     p adj
    ## shocked-homecage 1.139539 -16.30937 18.58845 0.8916162

    aov6 <- aov(PC3 ~ Treatment, data=pcadata)
    summary(aov6) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment    1  451.5   451.5   2.712  0.119
    ## Residuals   16 2664.0   166.5

    TukeyHSD(aov6, which = "Treatment") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC3 ~ Treatment, data = pcadata)
    ## 
    ## $Treatment
    ##                      diff       lwr      upr     p adj
    ## shocked-homecage 10.62449 -3.052643 24.30161 0.1191049

    aov7 <- aov(PC4 ~ Treatment, data=pcadata)
    summary(aov7) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment    1   71.4   71.44   0.641  0.435
    ## Residuals   16 1782.4  111.40

    TukeyHSD(aov7, which = "Treatment") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC4 ~ Treatment, data = pcadata)
    ## 
    ## $Treatment
    ##                      diff       lwr      upr     p adj
    ## shocked-homecage 4.225979 -6.961505 15.41346 0.4349914

    lm1 <- lm(PC1~Region*Treatment, data=pcadata)
    summary(lm1)

    ## 
    ## Call:
    ## lm(formula = PC1 ~ Region * Treatment, data = pcadata)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -9.6733 -2.3861 -0.5788  3.0613  8.3686 
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                -14.6987     3.8260  -3.842  0.00234 ** 
    ## RegionCA3                    2.2656     5.4108   0.419  0.68282    
    ## RegionDG                    48.9297     5.4108   9.043 1.05e-06 ***
    ## Treatmentshocked            -0.8485     4.5270  -0.187  0.85446    
    ## RegionCA3:Treatmentshocked   3.5322     6.7001   0.527  0.60767    
    ## RegionDG:Treatmentshocked  -10.1860     6.5155  -1.563  0.14394    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 5.411 on 12 degrees of freedom
    ## Multiple R-squared:  0.9504, Adjusted R-squared:  0.9297 
    ## F-statistic: 45.98 on 5 and 12 DF,  p-value: 2.047e-07

    anova(lm1) 

    ## Analysis of Variance Table
    ## 
    ## Response: PC1
    ##                  Df Sum Sq Mean Sq  F value    Pr(>F)    
    ## Region            2 6558.8  3279.4 112.0137 1.727e-08 ***
    ## Treatment         1   40.7    40.7   1.3915    0.2610    
    ## Region:Treatment  2  131.3    65.6   2.2421    0.1488    
    ## Residuals        12  351.3    29.3                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    lm2 <- lm(PC2~Region*Treatment, data=pcadata)
    summary(lm2)

    ## 
    ## Call:
    ## lm(formula = PC2 ~ Region * Treatment, data = pcadata)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -12.732  -1.203   0.802   1.526  12.732 
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                 -16.598      4.150  -4.000  0.00176 ** 
    ## RegionCA3                    39.288      5.868   6.695 2.21e-05 ***
    ## RegionDG                      8.227      5.868   1.402  0.18629    
    ## Treatmentshocked              3.285      4.910   0.669  0.51617    
    ## RegionCA3:Treatmentshocked   -5.140      7.267  -0.707  0.49287    
    ## RegionDG:Treatmentshocked     7.242      7.066   1.025  0.32566    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 5.868 on 12 degrees of freedom
    ## Multiple R-squared:  0.9048, Adjusted R-squared:  0.8651 
    ## F-statistic: 22.81 on 5 and 12 DF,  p-value: 9.603e-06

    anova(lm2)

    ## Analysis of Variance Table
    ## 
    ## Response: PC2
    ##                  Df Sum Sq Mean Sq F value    Pr(>F)    
    ## Region            2 3760.6 1880.31 54.6018 9.419e-07 ***
    ## Treatment         1   68.7   68.73  1.9957    0.1832    
    ## Region:Treatment  2   98.6   49.28  1.4310    0.2771    
    ## Residuals        12  413.2   34.44                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## DEG by contrasts
    source("resvalsfunction.R")
    contrast1 <- resvals(contrastvector = c('Region', 'CA1', 'DG'), mypval = 0.1)

    ## [1] 2674
    ## [1] 674

    contrast2 <- resvals(contrastvector = c('Region', 'CA3', 'DG'), mypval = 0.1)

    ## [1] 2849
    ## [1] 871

    contrast3 <- resvals(contrastvector = c('Region', 'CA1', 'CA3'), mypval = 0.1)

    ## [1] 2510
    ## [1] 384

    contrast4 <- resvals(contrastvector = c('Treatment', 'shocked', 'homecage'), mypval = 0.1)

    ## [1] 1550
    ## [1] 6

Now, we can view a histogram of the distribution

    source("resvalsfunction.R")
    myhistogram(contrastvector = c('Region', 'CA1', 'DG'), mypval = 0.1)

![](../figures/02_stresstest/histogram-1.png)

    ## [1] 1

    myhistogram(contrastvector = c('Region', 'CA3', 'DG'), mypval = 0.1)

![](../figures/02_stresstest/histogram-2.png)

    ## [1] 1

    myhistogram(contrastvector = c('Region', 'CA1', 'CA3'), mypval = 0.1)

![](../figures/02_stresstest/histogram-3.png)

    ## [1] 1

    myhistogram(contrastvector = c('Treatment', 'shocked', 'homecage'), mypval = 0.1)

![](../figures/02_stresstest/histogram-4.png)

    ## [1] 1

This Venn Diagram sthe overlap of differentailly expression genes by
Region and method. This shows all genes with *adjusted* pvalue &lt;0.1.

    #create a new DF with the gene counts
    rldpvals <- assay(rld)
    rldpvals <- cbind(rldpvals, contrast1, contrast2, contrast3, contrast4)
    rldpvals <- as.data.frame(rldpvals)
    rldpvals <- rldpvals[ , grepl( "padj|pval" , names( rldpvals ) ) ]


    # venn with padj values
    venn1 <- row.names(rldpvals[rldpvals[2] <0.1 & !is.na(rldpvals[2]),])
    venn2 <- row.names(rldpvals[rldpvals[4] <0.1 & !is.na(rldpvals[4]),])
    venn3 <- row.names(rldpvals[rldpvals[6] <0.1 & !is.na(rldpvals[6]),])
    venn4 <- row.names(rldpvals[rldpvals[8] <0.1 & !is.na(rldpvals[8]),])
    venn12 <- union(venn1,venn2)
    venn123 <- union(venn12,venn3)

    # save files for big venn diagram
    write(venn123, "../results/02_stress_venn123.txt")
    write(venn4, "../results/02_stress_venn4.txt")


    ## check order for correctness
    candidates <- list("Region" = venn123, "Method" = venn4)

    prettyvenn <- venn.diagram(
      scaled=T,
      x = candidates, filename=NULL, 
      col = "black",
      fill = c( "white", "white"),
      alpha = 0.5,
      cex = 1, fontfamily = "sans", #fontface = "bold",
      cat.default.pos = "text",
      cat.dist = c(0.07, 0.07), cat.pos = 1,
      cat.cex = 1, cat.fontfamily = "sans")
    #dev.off()
    grid.draw(prettyvenn)

![](../figures/02_stresstest/VennDiagramPadj-1.png)

    ## Any padj <0.1
    DEGes <- assay(rld)
    DEGes <- cbind(DEGes, contrast1, contrast2, contrast3, contrast4)
    DEGes <- as.data.frame(DEGes) # convert matrix to dataframe
    DEGes$rownames <- rownames(DEGes)  # add the rownames to the dataframe

    DEGes$padjmin <- with(DEGes, pmin(padjTreatmentshockedhomecage, padjRegionCA1DG ,padjRegionCA3DG, padjRegionCA1CA3 )) # put the min pvalue in a new column
    DEGes <- DEGes %>% filter(padjmin < 0.1)

    rownames(DEGes) <- DEGes$rownames
    drop.cols <-colnames(DEGes[,grep("padj|pval|rownames", colnames(DEGes))])
    DEGes <- DEGes %>% select(-one_of(drop.cols))
    DEGes <- as.matrix(DEGes)
    DEGes <- DEGes - rowMeans(DEGes)


    # setting color options
    source("figureoptions.R")
    ann_colors <- ann_colorsstress
    colorpalette <- cembrowskicolors
    df <- as.data.frame(colData(dds)[,c("Treatment", "Region")])
    df$Treatment <- factor(df$Treatment, levels = c("homecage", "shocked"))


    paletteLength <- 30
    myBreaks <- c(seq(min(DEGes), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(DEGes)/paletteLength, max(DEGes), length.out=floor(paletteLength/2)))


    pheatmap(DEGes, show_colnames=T, show_rownames = F,
             annotation_col=df, annotation_colors = ann_colors,
             treeheight_row = 0, treeheight_col = 25,
             fontsize = 8, 
             #width=4.5, height=3,
             border_color = "grey60" ,
             color = colorpalette,
             #cellwidth = 12, 
             #main = "Any Padj < 0.1",
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation" 
             )

![](../figures/02_stresstest/HeatmapPadj-1.png)

    # for adobe
    pheatmap(DEGes, show_colnames=F, show_rownames = F,
             annotation_col=df, annotation_colors = ann_colors,
             treeheight_row = 0, treeheight_col = 25,
             fontsize = 8, 
             width=4.5, height=3,
             border_color = "grey60" ,
             color = colorpalette,
             cellwidth = 7, 
             filename = "../figures/02_stresstest/HeatmapPadj-1.pdf",
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation" 
             )

    library(edgeR)
    counts <- countData
    dim( counts )

    ## [1] 22485    18

    colSums( counts ) / 1e06  # in millions of gene counts

    ## 143B-CA1-1  143B-DG-1 144B-CA1-1 144B-CA3-1 145B-CA1-1  145B-DG-1 
    ##   0.874614   1.019113   1.275137   0.506698   1.034066   0.720798 
    ## 146B-CA1-2 146B-CA3-2  146B-DG-2  147-CA1-4  147-CA3-4   147-DG-4 
    ##   0.506014   1.056001   0.055549   0.080721   0.344588   0.069648 
    ##  148-CA1-2  148-CA3-2   148-DG-2 148B-CA1-4 148B-CA3-4  148B-DG-4 
    ##   0.938866   1.148136   1.067185   0.185637   1.724144   0.398258

    table( rowSums( counts ) )[ 1:30 ] # Number of genes with low counts

    ## 
    ##    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
    ## 5881  459  423  280  261  209  190  159  178  156  136  141  139  107  104 
    ##   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29 
    ##   95   91   94   81   79   84   69   83   81   68   80   77   58   49   46

    # from https://github.com/rachelwright8/Ahya-White-Syndromes/blob/master/deseq2_Ahya.R

    res <- results(dds, contrast=c('Treatment', 'shocked', 'homecage'))
    table(res$padj<0.1)

    ## 
    ## FALSE  TRUE 
    ##  9295    18

    table(res$pvalue<0.05)

    ## 
    ## FALSE  TRUE 
    ## 14851   862

    head(res)

    ## log2 fold change (MLE): Treatment shocked vs homecage 
    ## Wald test p-value: Treatment shocked vs homecage 
    ## DataFrame with 6 rows and 6 columns
    ##                baseMean log2FoldChange     lfcSE       stat    pvalue
    ##               <numeric>      <numeric> <numeric>  <numeric> <numeric>
    ## 0610007P14Rik 12.542367      0.1356020  0.897367  0.1511110 0.8798882
    ## 0610009B22Rik  6.112886      0.4614292  1.226013  0.3763657 0.7066451
    ## 0610009L18Rik  1.423716     -1.2901378  2.187262 -0.5898416 0.5552969
    ## 0610009O20Rik 24.090980      1.3885828  0.962880  1.4421141 0.1492702
    ## 0610010F05Rik  5.411247     -1.0611824  1.108837 -0.9570231 0.3385556
    ## 0610010K14Rik  1.218209      1.4258700  1.908507  0.7471129 0.4549954
    ##                    padj
    ##               <numeric>
    ## 0610007P14Rik 0.9979476
    ## 0610009B22Rik 0.9905416
    ## 0610009L18Rik        NA
    ## 0610009O20Rik 0.7973706
    ## 0610010F05Rik 0.9317371
    ## 0610010K14Rik        NA

    logs <- data.frame(cbind("gene"=row.names(res),"logP"=round(-log(res$pvalue+1e-10,10),1)))
    logs$logP <- as.numeric(as.character(logs$logP))
    sign <- rep(1,nrow(logs))
    sign[res$log2FoldChange<0]=-1  ##change to correct model
    table(sign)

    ## sign
    ##   -1    1 
    ## 7250 8472

    logs$logP <- logs$logP*sign

    write.csv(logs, file = "./06_GO_MWU/02_stress_GOpvals.csv", row.names = F)

summary stats
=============

    ## now, to show how many shocks each animal received.
    unavoidableshock <- read.csv("../data/unavoidableshock.csv", header = T)

    unavoidableshock %>% 
      ggplot(aes(x=TrainSessionCombo.x, y=NumShock.y)) + 
      geom_point(position = position_jitter(width = 0.2)) + 
      #geom_point()+ 
      theme_bw() + 
      scale_x_discrete(labels=c("T1" = "1", "T2" = "2",
                                "T3" = "3", "T4_C1" = "4",
                                "T5_C2" = "5", "T6_C3" = "6")) + 
      labs(list(x = "Day", y = "Number of Shocks")) +
      theme(axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=12),
            axis.text.x = element_text(size = 11),
            axis.text.y = element_text(size = 11))

![](../figures/02_stresstest/numshocks-1.png)

    library(pvclust)
    # clustering just the degs
    result <- pvclust(DEGes, method.dist="cor", method.hclust="average", nboot=1000)

    ## Bootstrap (r = 0.5)... Done.
    ## Bootstrap (r = 0.6)... Done.
    ## Bootstrap (r = 0.7)... Done.
    ## Bootstrap (r = 0.8)... Done.
    ## Bootstrap (r = 0.9)... Done.
    ## Bootstrap (r = 1.0)... Done.
    ## Bootstrap (r = 1.1)... Done.
    ## Bootstrap (r = 1.2)... Done.
    ## Bootstrap (r = 1.3)... Done.
    ## Bootstrap (r = 1.4)... Done.

    plot(result)

![](../figures/02_stresstest/pvclust-1.png)
