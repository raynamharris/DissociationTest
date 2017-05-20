Examining of cognitive training on hippocampal transcriptomes
-------------------------------------------------------------

### Experimental Design: Treatment + Region

Next, we examined the effects of cognitiving training on hippocampal
gene expression. Mice trained in the active place avoidance task alter
their behavior to avoid footshocks. The sequence of unavoidable shock
delivered to the yoked mice mimicked the time series of shocks received
by the trained mice that were being conditioned to avoid localized, mild
shocks. While the trained and yoked animals received the same number of
shocks, only the trained animals exhibitied an avoidance response
(Supplementary figures showing the number of shocks and the avoidance
behaviors can be viewed by using 'include=TRUE' in the corresponding Rmd
file).

The animals used in experiments 2 and 3 were 3-4–month-old male C57BL/6J
mice (Jackson Laboratory). Following killed trained mice (N=4) were
yoked control mice (N=4). Mice were killed and transverse brain slices
were prepared. The DG, CA3, CA1 regions were microdissected using a 0.25
mm punch (Electron Microscopy Systems) and a dissecting scope (Zeiss).
RNA was isolated using the Maxwell 16 LEV RNA Isolation Kit (Promega).
RNA libraries were prepared by the Genomic Sequencing and Analysis
Facility at the University of Texas at Austin using the Illumina HiSeq
platform.

<img src="../figures/03_cognitiontest/03_biologicalsamples-01.png" width="297" />

    ##    Treatment  Region 
    ##  yoked  : 9   CA1:8  
    ##  trained:13   CA3:5  
    ##               DG :9

### Gene counts

Raw reads were downloaded from the Amazon cloud server to the Stampede
Cluster at the Texas Advanced Computing Facility for processing and
analysis. RNA quality was checked using the bioinformatic program FASTQC
(citation). Low quality reads and adapter sequences were removed using
the program Cutadapt (Martin, 2011). Kallisto was use for fast read
mapping and counting (Bray et al., 2016). Transcript from a single gene
were combined into a count total for each gene. In the end, we meausred
the expression of 22,485 genes in 22 samples.

    dim(countData)

    ## [1] 22485    22

### Differential gene expresssion analysis

We used DESeq2 for gene expression normalization and quantification
(Love et al., 2014). We compared the effects of treatment, region, and
the interaction with the formal
`design = ~ Treatment + Region + Treatment * Region`. After removing
genes with less than 2 counts across all samples with
`dds <- dds[ rowSums(counts(dds)) > 2, ]`, we were left with 16,970
genes.

    dds

    FALSE class: DESeqDataSet 
    FALSE dim: 16970 22 
    FALSE metadata(1): version
    FALSE assays(3): counts mu cooks
    FALSE rownames(16970): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    FALSE rowData names(37): baseMean baseVar ... deviance maxCooks
    FALSE colnames(22): 142C_CA1 142C_DG ... 147D-CA3-1 147D-DG-1
    FALSE colData names(15): RNAseqID Mouse ... Date sizeFactor

This experiment produced a larger effect on gene expression, 285 genes
were differentially expressed between the two groups, and a large
portion of those genes were also being differentially expressed between
the hippocampal subfields. Hierarchical clustering of the differentially
expressed genes separates samples by both subfield and treatment (Fig.
3C).

This Venn Diagram sthe overlap of differentailly expression genes by
Region and method. This shows all genes with adjusted pvalue according
to the set pvalue.

![](../figures/03_cognitiontest/VennDiagramPadj-1.png)

Then, we visuazlied the data as a heatmap showing the relative log fold
change of gene expression across samples. Genes were filtered for a
minimimum adjust p value &lt; 0.1 in any two-way contrast. The row mean
for each gene was subtracted for the raw value to allow for analysis of
fold change rather than raw magnitudes. The samples cluster primarily by
brain region with some small treatment-driven.

![](../figures/03_cognitiontest/HeatmapPadj-1.png)

Then, we used pvclust to obtain bootstrap values for the heatmap sample
dendrogram.

`{r pvclust, echo=FALSE, message=FALSE, comment=FALSE, warning=FALSE} library(pvclust) #clustering just the degs result <- pvclust(DEGes, method.dist="cor", method.hclust="average", nboot=1000) plot(result) #`
=================================================================================================================================================================================================================

(Supplementary figures showing the distibution of pvalues can be viewed
by using 'include=TRUE' in the corresponding Rmd file).

### Analysis of variance

A principal component analysis of all gene expression data revealed that
PC1 explains 50% of the variance in gene expression and distinguishes
between the DG and CA regions. PC2 accounts for 18% of the variance and
distinguishes the three subfields.

    source("DESeqPCAfunction.R")
    source("figureoptions.R")

    # create the dataframe using my function pcadataframe
    pcadata <- pcadataframe(rld, intgroup=c("Treatment", "Region"), returnData=TRUE)
    percentVar <- round(100 * attr(pcadata, "percentVar"))

    pcadata$Treatment <- factor(pcadata$Treatment, levels = c("yoked", "trained"))

    ## PC2 vs PC1
    plotPC2PC1(aescolor = pcadata$Region, colorname = "Region", aesshape = pcadata$Treatment, shapename = "Treatment", colorvalues = colorvalRegion)

![](../figures/03_cognitiontest/unnamed-chunk-7-1.png)

    # PC1 vs PC2 for adobe
    myplot <- plotPC2PC1(aescolor = pcadata$Region, colorname = "Region", aesshape = pcadata$Treatment, shapename = "Treatment", colorvalues = colorvalRegion)
    pdf(file="../figures/03_cognitiontest/PCA-1.pdf", width=4.5, height=3)
    plot(myplot)
    dev.off()

    FALSE quartz_off_screen 
    FALSE                 2

To confirm statistical significance of this visual pattern, we conducted
a two-way ANOVA for PC1 ~ Region: F2,19= 199.3; p = 1.78e-13). A post
hoc Tukey test showed that DG samples are significantly different from
both CA1 and CA3 samples (CA1-DG, p = 1.0e-07; CA3-DG, p = 1.0e-07;
CA1-CA3, p = 0.7002).

    aov1 <- aov(PC1 ~ Region, data=pcadata)
    summary(aov1) 

    FALSE             Df Sum Sq Mean Sq F value   Pr(>F)    
    FALSE Region       2  10887    5443   199.3 1.78e-13 ***
    FALSE Residuals   19    519      27                     
    FALSE ---
    FALSE Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov1, which = "Region")

    FALSE   Tukey multiple comparisons of means
    FALSE     95% family-wise confidence level
    FALSE 
    FALSE Fit: aov(formula = PC1 ~ Region, data = pcadata)
    FALSE 
    FALSE $Region
    FALSE              diff       lwr       upr     p adj
    FALSE CA3-CA1  2.419769 -5.149283  9.988821 0.7002216
    FALSE DG-CA1  46.138194 39.686734 52.589654 0.0000000
    FALSE DG-CA3  43.718425 36.312871 51.123979 0.0000000

The strongest contributor to PC2 is brain regions (PC2 ~ Region ANOVA:
F2,19= 220.4; p = 7.15e-14; Tukey test, p&lt;&lt;&lt;0.001 for all three
comparisons), with some influence of treatment on PC2 (PC2 ~ Treatment
ANOVA: F1,20=3.389; p = 0.0805).

    aov2 <- aov(PC2 ~ Region, data=pcadata)
    summary(aov2) 

    FALSE             Df Sum Sq Mean Sq F value   Pr(>F)    
    FALSE Region       2   3934  1966.9   220.4 7.15e-14 ***
    FALSE Residuals   19    170     8.9                     
    FALSE ---
    FALSE Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov2, which = "Region") 

    FALSE   Tukey multiple comparisons of means
    FALSE     95% family-wise confidence level
    FALSE 
    FALSE Fit: aov(formula = PC2 ~ Region, data = pcadata)
    FALSE 
    FALSE $Region
    FALSE              diff        lwr       upr p adj
    FALSE CA3-CA1  35.74129  31.414508  40.06807 0e+00
    FALSE DG-CA1   12.96556   9.277639  16.65348 1e-07
    FALSE DG-CA3  -22.77573 -27.009049 -18.54241 0e+00

PC3 and PC4 account for 7% and 4.5 of the variation in gene expression
respectively.PC3 are PC4 are influenced by variation due to treatment
(PC3 ~ Treatment ANOVA: F1,18=5.622; p = 0.0291, PC4 ~ Treatment ANOVA:
F1,18=12.01; p = 0.00276).

    source("DESeqPCAfunction.R")
    source("figureoptions.R")
    ## PC2 vs PC3
    A <- plotPC2PC3(aescolor = pcadata$Region, colorname = "Region", aesshape = pcadata$Treatment, shapename = "Treatment", colorvalues = colorvalRegion)

    B <- plotPC2PC4(aescolor = pcadata$Region, colorname = "Region", aesshape = pcadata$Treatment, shapename = "Treatment", colorvalues = colorvalRegion)

    plot_grid(A, B, rel_widths = c(1,1))

![](../figures/03_cognitiontest/PCA234figures-1.png)

    aov3 <- aov(PC3 ~ Treatment, data=pcadata)
    summary(aov3) 

    FALSE             Df Sum Sq Mean Sq F value Pr(>F)  
    FALSE Treatment    1  417.6   417.6   7.629  0.012 *
    FALSE Residuals   20 1094.8    54.7                 
    FALSE ---
    FALSE Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    aov4 <- aov(PC4 ~ Treatment, data=pcadata)
    summary(aov4) 

    FALSE             Df Sum Sq Mean Sq F value Pr(>F)   
    FALSE Treatment    1  401.3   401.3   10.43 0.0042 **
    FALSE Residuals   20  769.5    38.5                  
    FALSE ---
    FALSE Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

The gene expression data were exported to csv files for importing into
the GOMMU analysis package for subsequent analysis.
