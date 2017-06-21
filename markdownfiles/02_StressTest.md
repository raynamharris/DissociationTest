Identifying the effects of habituation to a stressor on hippocampal transcriptomes
----------------------------------------------------------------------------------

We examined the expression patterns of 16,229 genes. We identified 0
genes that were significantly expressed between homecage and shocked
samples; 1669 genes that were were differentially expressed between any
of the three brain regions at PDF p-value &lt; 0.05 (Fig. 3B).

    #source("http://www.bioconductor.org/biocLite.R")
    library(DESeq2)
    library(VennDiagram)
    library(genefilter)
    library(pheatmap)
    library(cowplot)
    library(RColorBrewer)
    library(dplyr)
    library(plyr)
    library(ggplot2)
    library(edgeR)
    library(pvclust)
    library(knitr) 

    # set output file for figures 
    knitr::opts_chunk$set(fig.path = '../figures/02_stresstest/')

    # user defined funcitons and options
    source("resvalsfunction.R") 
    source("DESeqPCAfunction.R")
    source("figureoptions.R")

The sample and count information for this part is found in
`../data/GSE100225_IntegrativeWT2015ColData.csv` and
`../data/GSE100225_IntegrativeWT2015CountData.csv`. You can also
download these two files (with a different name but same content) from
[GEO
GSE100225](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100225).

    colData <- read.csv('../data/GSE100225_IntegrativeWT2015ColData.csv')
    rownames(colData) <- colData$RNAseqID
    countData <-  read.csv('../data/GSE100225_IntegrativeWT2015CountData.csv', check.names = F)

    colData <- colData %>%
      filter(Treatment %in% c("homecage", "shocked")) %>% droplevels()

    savecols <- as.character(colData$RNAseqID) #selects all good samples
    savecols <- as.vector(savecols) # make it a vector
    countData <- countData %>% select(one_of(savecols)) # keep good samples

    ## rename and relevel things
    colData$Treatment <- factor(colData$Treatment, levels = c("homecage", "shocked"))

Here is a brief overview of the samples being compared.

    ##     Treatment  Region 
    ##  homecage: 6   CA1:7  
    ##  shocked :12   CA3:5  
    ##                DG :6

18 Samples, 22485 genes.

    dim(countData)

    ## [1] 22485    18

    counts <- countData
    dim( counts )

    ## [1] 22485    18

    colSums( counts ) / 1e06  # in millions of gene counts

    ## 143B-CA1-1  143B-DG-1 144B-CA1-1 144B-CA3-1 145B-CA1-1  145B-DG-1 
    ##   1.719498   2.085031   2.555909   1.027388   2.020114   1.509310 
    ## 146B-CA1-2 146B-CA3-2  146B-DG-2  147-CA1-4  147-CA3-4   147-DG-4 
    ##   1.063417   2.144771   0.116106   0.159069   0.689232   0.139276 
    ##  148-CA1-2  148-CA3-2   148-DG-2 148B-CA1-4 148B-CA3-4  148B-DG-4 
    ##   1.901256   2.343035   2.231849   0.337174   3.486840   0.798668

    table( rowSums( counts ) )[ 1:30 ] # Number of genes with low counts

    ## 
    ##    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
    ## 5664  331  261  236  209  167  130  147  132  123  115  125   85   69   92 
    ##   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29 
    ##   95   79   65   67   65   64   71   77   60   54   57   52   48   40   46

    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ Treatment + Region + Treatment * Region )
    dds <- dds[ rowSums(counts(dds)) > 2, ] ## filter genes with 0 counts
    dds <- DESeq(dds) # Differential expression analysis

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    rld <- rlog(dds, blind=FALSE) # log transformed
    dim(rld) #print total genes analyzed

    ## [1] 16229    18

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

We examined the expression patterns of 16,229 genes. We identified 0
genes that were significantly expressed between homecage and shocked
samples; 1669 genes that were were differentially expressed between any
of the three brain regions at PDF p-value &lt; 0.05 (Fig. 3B).

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

Supplementary histogram of p-value distributions

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

Hierarchical clustering of the differentially expressed genes gives rise
to three distinct clusters corresponding to the three subfields, with
CA1 (purple) and CA3 (green) being more similar to one another than to
DG (orange) (Fig. 3C).

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

Next, we conducted a principal component analysis of all genes measured.
PC1 accounts for 31% of the variation and visually separates the DG
samples from the CA1 and CA3 samples ANOVA (PC1 ~ Region, F2,15= 42.89;
p &lt; 0.001) (Fig. 3D). A post hoc Tukey test showed that DG samples
are significantly different from both CA1 and CA3 samples (CA1-DG, p
&lt; 0.001; CA3-DG, p &lt; 0.001; CA1-CA3, p = 0.83). PC2 accounts for
18% of the variation and varies significantly between CA1 and CA3 and
CA1 and DG (PC2 ~ Region, ANOVA, F2, 15= 11.41; p &lt; 0.001; Tukey
test, CA1-DG, p = 0.03; CA3-DG, p = 0.18; CA1-CA3, p &lt; 0.001). PC2
accounts for 15% of the variation and also explains some brain region
specific differences (PC3 ~ Region, ANOVA, F2, 15= 6.315; p &lt; 0.01;
Tukey test, CA1-DG, p = 0.95; CA3-DG, p = 0.03; CA1-CA3, p = 0.01). PC7
is the lowest PC to explain any variance associated with treatment (PC6
~ Region, ANOVA, F1, 16= 4.774; p = 0.04

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
    myplot <- plotPC1PC2(aescolor = pcadata$Region, colorname = "Region", aesshape = pcadata$Treatment, shapename = "Treatment", colorvalues = colorvalRegion)
    pdf(file="../figures/../figures/02_stresstest/PCA-1.pdf", width=4.5, height=3)
    plot(myplot)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    ## statistics
    aov1R <- aov(PC1 ~ Region, data=pcadata)
    summary(aov1R) 

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Region       2   7334    3667   42.89 6.24e-07 ***
    ## Residuals   15   1282      85                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov1R, which = "Region") 

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

    aov2R <- aov(PC2 ~ Region, data=pcadata)
    summary(aov2R) 

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Region       2   2976  1487.9   11.41 0.000971 ***
    ## Residuals   15   1955   130.4                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov2R, which = "Region") 

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

    aov3R <- aov(PC3 ~ Region, data=pcadata)
    summary(aov3R) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## Region       2   1928   963.8   6.315 0.0102 *
    ## Residuals   15   2289   152.6                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov3R, which = "Region") 

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

    aov4R <- aov(PC3 ~ Region, data=pcadata)
    summary(aov4R) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## Region       2   1928   963.8   6.315 0.0102 *
    ## Residuals   15   2289   152.6                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov4R, which = "Region") 

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

    aov5R <- aov(PC5 ~ Region, data=pcadata)
    summary(aov5R) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Region       2   17.7    8.86   0.092  0.912
    ## Residuals   15 1439.5   95.97

    TukeyHSD(aov5R, which = "Region") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC5 ~ Region, data = pcadata)
    ## 
    ## $Region
    ##               diff       lwr      upr     p adj
    ## CA3-CA1 -0.2227589 -15.12229 14.67677 0.9991689
    ## DG-CA1  -2.1893719 -16.34610 11.96736 0.9153862
    ## DG-CA3  -1.9666130 -17.37480 13.44158 0.9414470

    aov6R <- aov(PC6 ~ Region, data=pcadata)
    summary(aov6R) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Region       2    1.3    0.66   0.007  0.993
    ## Residuals   15 1387.3   92.49

    TukeyHSD(aov6R, which = "Region") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC6 ~ Region, data = pcadata)
    ## 
    ## $Region
    ##              diff       lwr      upr     p adj
    ## CA3-CA1 0.4724674 -14.15423 15.09916 0.9961274
    ## DG-CA1  0.6061285 -13.29137 14.50363 0.9929528
    ## DG-CA3  0.1336611 -14.99238 15.25970 0.9997096

    aov1T <- aov(PC1 ~ Treatment, data=pcadata)
    summary(aov1T) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment    1    142   142.1   0.268  0.612
    ## Residuals   16   8474   529.6

    aov2T <- aov(PC2 ~ Treatment, data=pcadata)
    summary(aov2T) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment    1     53   52.73   0.173  0.683
    ## Residuals   16   4879  304.91

    aov3T <- aov(PC3 ~ Treatment, data=pcadata)
    summary(aov3T) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment    1    339   339.4     1.4  0.254
    ## Residuals   16   3877   242.3

    aov4T <- aov(PC4 ~ Treatment, data=pcadata)
    summary(aov4T) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment    1  172.9   172.9   1.435  0.248
    ## Residuals   16 1927.7   120.5

    aov5T <- aov(PC5 ~ Treatment, data=pcadata)
    summary(aov5T) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment    1  122.8  122.75   1.472  0.243
    ## Residuals   16 1334.5   83.41

    aov6T <- aov(PC6 ~ Treatment, data=pcadata)
    summary(aov6T) 

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

Next, save files for dowstream GO analysis.

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

Supplementary behavior file about timesheries of shocks.

Here is the corresponding Adobe Illustrator file that combines many of
the above plots.

<img src="../figures/02_stresstest/02_stress-01.png" width="1370" />
