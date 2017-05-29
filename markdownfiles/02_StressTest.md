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

Here is a brief overview of the samples being compared.

    ##     Treatment  Region 
    ##  homecage: 6   CA1:7  
    ##  shocked :12   CA3:5  
    ##                DG :6

18 Samples, 22485 genes.

    ## [1] 22485    18

    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ Treatment + Region + Treatment * Region )
    dds <- dds[ rowSums(counts(dds)) > 2, ] ## filter genes with 0 counts
    dds <- DESeq(dds) # Differential expression analysis
    rld <- rlog(dds, blind=FALSE) # log transformed
    dim(rld) #print total genes analyzed

    ## [1] 16229    18

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
    plotPC2PC1(aescolor = pcadata$Region, colorname = "Region", aesshape = pcadata$Treatment, shapename = "Treatment", colorvalues = colorvalRegion)

![](../figures/02_stresstest/PCA-1.png)

    plotPC2PC3(aescolor = pcadata$Region, colorname = "Region", aesshape = pcadata$Treatment, shapename = "Treatment", colorvalues = colorvalRegion)

![](../figures/02_stresstest/PCA-2.png)

    # for adobe
    myplot <- plotPC2PC1(aescolor = pcadata$Region, colorname = "Region", aesshape = pcadata$Treatment, shapename = "Treatment", colorvalues = colorvalRegion)
    pdf(file="../figures/../figures/02_stresstest/PCA-1.pdf", width=4.5, height=3)
    plot(myplot)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    ## statistics
    aov1 <- aov(PC1 ~ Region, data=pcadata)
    summary(aov1) 

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Region       2   7334    3667   42.89 6.24e-07 ***
    ## Residuals   15   1282      85                     
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
    ## CA3-CA1  3.127434 -10.93525 17.19012 0.8339245
    ## DG-CA1  44.037806  30.67620 57.39941 0.0000011
    ## DG-CA3  40.910372  26.36759 55.45315 0.0000073

    aov2 <- aov(PC2 ~ Region, data=pcadata)
    summary(aov2) 

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Region       2   2976  1487.9   11.41 0.000971 ***
    ## Residuals   15   1955   130.4                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov2, which = "Region") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ Region, data = pcadata)
    ## 
    ## $Region
    ##              diff        lwr       upr     p adj
    ## CA3-CA1  31.34187  13.976505 48.707233 0.0008014
    ## DG-CA1   18.32018   1.820548 34.819817 0.0288392
    ## DG-CA3  -13.02169 -30.979895  4.936522 0.1776538

    aov5 <- aov(PC3 ~ Region, data=pcadata)
    summary(aov5) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## Region       2   1928   963.8   6.315 0.0102 *
    ## Residuals   15   2289   152.6                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov5, which = "Region") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC3 ~ Region, data = pcadata)
    ## 
    ## $Region
    ##               diff       lwr       upr     p adj
    ## CA3-CA1 -23.912642 -42.70222 -5.123068 0.0125516
    ## DG-CA1   -1.905862 -19.75870 15.946981 0.9586161
    ## DG-CA3   22.006780   2.57574 41.437821 0.0257680

    aov3 <- aov(PC1 ~ Treatment, data=pcadata)
    summary(aov3) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment    1    142   142.1   0.268  0.612
    ## Residuals   16   8474   529.6

    aov4 <- aov(PC2 ~ Treatment, data=pcadata)
    summary(aov4) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment    1     53   52.73   0.173  0.683
    ## Residuals   16   4879  304.91

    aov6 <- aov(PC3 ~ Treatment, data=pcadata)
    summary(aov6) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment    1    339   339.4     1.4  0.254
    ## Residuals   16   3877   242.3

    aov7 <- aov(PC6 ~ Treatment, data=pcadata)
    summary(aov7) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## Treatment    1  319.1   319.1   4.774 0.0441 *
    ## Residuals   16 1069.5    66.8                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    lm1 <- lm(PC1~Region*Treatment, data=pcadata)
    summary(lm1)

    ## 
    ## Call:
    ## lm(formula = PC1 ~ Region * Treatment, data = pcadata)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -14.6047  -2.6399  -0.3964   3.2833  14.6047 
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                -15.1706     5.4511  -2.783   0.0166 *  
    ## RegionCA3                    0.2547     7.7090   0.033   0.9742    
    ## RegionDG                    57.1767     7.7090   7.417 8.09e-06 ***
    ## Treatmentshocked            -0.5283     6.4498  -0.082   0.9361    
    ## RegionCA3:Treatmentshocked   4.6873     9.5459   0.491   0.6323    
    ## RegionDG:Treatmentshocked  -19.7460     9.2829  -2.127   0.0548 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 7.709 on 12 degrees of freedom
    ## Multiple R-squared:  0.9172, Adjusted R-squared:  0.8827 
    ## F-statistic:  26.6 on 5 and 12 DF,  p-value: 4.223e-06

    anova(lm1) 

    ## Analysis of Variance Table
    ## 
    ## Response: PC1
    ##                  Df Sum Sq Mean Sq F value    Pr(>F)    
    ## Region            2 7333.5  3666.8 61.7000 4.846e-07 ***
    ## Treatment         1  131.2   131.2  2.2072   0.16316    
    ## Region:Treatment  2  438.1   219.0  3.6855   0.05652 .  
    ## Residuals        12  713.1    59.4                      
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
    ## -26.3267  -1.3168   0.2702   1.4777  26.3267 
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)                 -16.942      8.011  -2.115  0.05606 . 
    ## RegionCA3                    34.630     11.330   3.057  0.00996 **
    ## RegionDG                      8.933     11.330   0.788  0.44573   
    ## Treatmentshocked              2.980      9.479   0.314  0.75860   
    ## RegionCA3:Treatmentshocked   -4.913     14.030  -0.350  0.73227   
    ## RegionDG:Treatmentshocked    14.293     13.643   1.048  0.31543   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 11.33 on 12 degrees of freedom
    ## Multiple R-squared:  0.6876, Adjusted R-squared:  0.5575 
    ## F-statistic: 5.283 on 5 and 12 DF,  p-value: 0.008546

    anova(lm2)

    ## Analysis of Variance Table
    ## 
    ## Response: PC2
    ##                  Df  Sum Sq Mean Sq F value   Pr(>F)   
    ## Region            2 2975.79 1487.90 11.5909 0.001575 **
    ## Treatment         1  157.38  157.38  1.2260 0.289883   
    ## Region:Treatment  2  257.64  128.82  1.0035 0.395373   
    ## Residuals        12 1540.42  128.37                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## DEG by contrasts
    source("resvalsfunction.R")
    contrast1 <- resvals(contrastvector = c('Region', 'CA1', 'DG'), mypval = 0.05)

    ## [1] 918

    contrast2 <- resvals(contrastvector = c('Region', 'CA3', 'DG'), mypval = 0.05)

    ## [1] 1024

    contrast3 <- resvals(contrastvector = c('Region', 'CA1', 'CA3'), mypval = 0.05)

    ## [1] 357

    contrast4 <- resvals(contrastvector = c('Treatment', 'shocked', 'homecage'), mypval = 0.05)

    ## [1] 0

Now, we can view a histogram of the distribution

    source("resvalsfunction.R")
    myhistogram(contrastvector = c('Region', 'CA1', 'DG'), mypval = 0.05)

![](../figures/02_stresstest/histogram-1.png)

    ## [1] 1

    myhistogram(contrastvector = c('Region', 'CA3', 'DG'), mypval = 0.05)

![](../figures/02_stresstest/histogram-2.png)

    ## [1] 1

    myhistogram(contrastvector = c('Region', 'CA1', 'CA3'), mypval = 0.05)

![](../figures/02_stresstest/histogram-3.png)

    ## [1] 1

    myhistogram(contrastvector = c('Treatment', 'shocked', 'homecage'), mypval = 0.05)

![](../figures/02_stresstest/histogram-4.png)

    ## [1] 1

This Venn Diagram sthe overlap of differentailly expression genes by
Region and method. This shows all genes with *adjusted* pvalue &lt;0.05.

    #create a new DF with the gene counts
    rldpvals <- assay(rld)
    rldpvals <- cbind(rldpvals, contrast1, contrast2, contrast3, contrast4)
    rldpvals <- as.data.frame(rldpvals)
    rldpvals <- rldpvals[ , grepl( "padj|pval" , names( rldpvals ) ) ]


    # venn with padj values
    venn1 <- row.names(rldpvals[rldpvals[2] <0.05 & !is.na(rldpvals[2]),])
    venn2 <- row.names(rldpvals[rldpvals[4] <0.05 & !is.na(rldpvals[4]),])
    venn3 <- row.names(rldpvals[rldpvals[6] <0.05 & !is.na(rldpvals[6]),])
    venn4 <- row.names(rldpvals[rldpvals[8] <0.05 & !is.na(rldpvals[8]),])
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

    ## Any padj <0.05
    DEGes <- assay(rld)
    DEGes <- cbind(DEGes, contrast1, contrast2, contrast3, contrast4)
    DEGes <- as.data.frame(DEGes) # convert matrix to dataframe
    DEGes$rownames <- rownames(DEGes)  # add the rownames to the dataframe

    DEGes$padjmin <- with(DEGes, pmin(padjTreatmentshockedhomecage, padjRegionCA1DG ,padjRegionCA3DG, padjRegionCA1CA3 )) # put the min pvalue in a new column
    DEGes <- DEGes %>% filter(padjmin < 0.05)

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
             #main = "Any Padj < 0.05",
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

    # from https://github.com/rachelwright8/Ahya-White-Syndromes/blob/master/deseq2_Ahya.R

    res <- results(dds, contrast=c('Treatment', 'shocked', 'homecage'), independentFiltering = F)
    table(res$padj<0.05)

    ## 
    ## FALSE 
    ## 16220

    logs <- data.frame(cbind("gene"=row.names(res),"logP"=round(-log(res$pvalue+1e-10,10),1)))
    logs$logP <- as.numeric(as.character(logs$logP))
    sign <- rep(1,nrow(logs))
    sign[res$log2FoldChange<0]=-1  ##change to correct model
    table(sign)

    ## sign
    ##   -1    1 
    ## 7045 9184

    logs$logP <- logs$logP*sign

    write.csv(logs, file = "./06_GO_MWU/02_stress_GOpvals.csv", row.names = F)

res &lt;- results(dds, contrast=c('Treatment', 'shocked', 'homecage'),
independentFiltering = F) table(res$padj&lt;0.05)
