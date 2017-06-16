Examining of cognitive training on hippocampal
----------------------------------------------

The goals of the subsequent analysis are 1) to determine the effects of
cognitiving training on hippocampal gene expression and 2) related any
detectable changes to variation cause by other technical and biological
treatements.

### Experimental Design

We use 3-4–month-old male C57BL/6J mice fro the Jackson Laboratory and
housed at the Marine Biological Laboratory. Mice (N=4) trained in the
active place avoidance task are conditioned to avoid mild shocks that
can be localized by visual cues in the enviornment. Yoked control mice
(N=4) are delivered sequence of unavoidable shock that mimickes the time
series of shocks received by the trained mice. While the trained and
yoked animals received the same number of shocks, only the trained
animals exhibitied an avoidance response. (Supplementary figures showing
the number of shocks and the avoidance behaviors can be viewed by using
'include=TRUE' in the corresponding Rmd file).

Thirty minutes after the last cognitive training session, mice were
killed and transverse brain slices were prepared. The DG, CA3, CA1
subregions were microdissected using a 0.25 mm punch (Electron
Microscopy Systems) and a dissecting scope (Zeiss). RNA was isolated
using the Maxwell 16 LEV RNA Isolation Kit (Promega). RNA libraries were
prepared by the Genomic Sequencing and Analysis Facility at the
University of Texas at Austin using the Illumina HiSeq platform.

<img src="../figures/03_cognitiontest/03_biologicalsamples-01.png" width="297" />

The orginal design was 4 animals per treament and 3 hippocampal sub
regions per animals, which would give 24 samples. After excluding
compromized samples, the final sample sizes are:

    ##    Treatment  Region 
    ##  yoked  : 9   CA1:8  
    ##  trained:13   CA3:5  
    ##               DG :9

### Differential gene expresssion analysis

Raw reads were downloaded from the Amazon cloud server to the Stampede
Cluster at the Texas Advanced Computing Facility for processing and
analysis. RNA quality was checked using the bioinformatic program FASTQC
(citation). Low quality reads and adapter sequences were removed using
the program Cutadapt (Martin, 2011). Kallisto was use for fast read
mapping and counting (Bray et al., 2016). Transcript from a single gene
were combined into a count total for each gene. In the end, we meausred
the expression of 22,485 genes in 22 samples.

    ## [1] 22485    22

We used DESeq2 (Love et al., 2014) for gene expression normalization and
quantification using the following experimental design:
`Treatment + Region + Treatment * Region`. Genes with less than 2 counts
across all samples were filtered, leaving us with `dim(rld)` genes for
analysis of differntial expression.

    dim(rld)

    FALSE [1] 17320    22

We identified 423 genes were differentially expressed between the yoked
control and cognitively trained animals, 3485 genes that were
differentially expressed across subfields, and 324 showed an interaction
at FDR p &lt; 0.05 (Fig. 4B). We see a large effect of brain region on
gene expression, with 20% of detectable genes begin differentially
expressed between one or more brain-region comparisons (3485
differentially expressed genes /17320 measured genes). This is an order
of magnitude greater than the 2% of the transcriptome that changed in
response to learning (423 DEGs /17320 genes measured).

![](../figures/03_cognitiontest/VennDiagramPadj-1.png)

Hierarchical clustering of the differentially expressed genes separates
samples by both subfield and treatment (Fig. 4C).

Then, we visuazlied the data as a heatmap showing the relative log fold
change of gene expression across samples. Genes were filtered for a
minimimum adjust p value &lt; 0.05 in any two-way contrast. The row mean
for each gene was subtracted for the raw value to allow for analysis of
fold change rather than raw magnitudes. The samples cluster primarily by
brain region with some small treatment-driven.

![](../figures/03_cognitiontest/HeatmapPadj-1.png)

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

![](../figures/03_cognitiontest/pvclust-1.png)

A principal component analysis of all gene expression data revealed that
brain region explains 75% of the variance in the data (Fig. 4D). PC1
accounts for 56% of the variance and distinguishes DG from not-DG
samples (ANOVA for PC1 ~ Region: F2,19= 226.1; p &lt; 0.001). A post hoc
Tukey test showed that DG samples are significantly different from both
CA1 and CA3 samples (CA1-DG, p &lt; 0.001; CA3-DG, p &lt; 0.001;
CA1-CA3, p = 0.23). PC2 accounts for 19% of the variance and
distinguishes the three subfields (PC2 ~ Region ANOVA: F2,19= 255.3; p
&lt; 0.001; Tukey test, p &lt; 0.001for all three comparisons). PC3 are
PC4 are influenced by variation due to cognitive training (PC3 ~
Treatment ANOVA: F1,20=7.451; p = 0.012, PC4 ~ Treatment ANOVA:
F1,20=10.11; p = 0.0047). An analysis of Gene Ontology identified 91 GO
terms at 10% FDR significant. Among the top 10 are glutamate signaling
and membrane transport systems and a downregulation of oxidoreductase
and ribosomal activity (Fig. 2C).

![](../figures/03_cognitiontest/PCA21-1.png)

    aov1 <- aov(PC1 ~ Region, data=pcadata)
    summary(aov1) 

    FALSE             Df Sum Sq Mean Sq F value   Pr(>F)    
    FALSE Region       2  12615    6307   226.1 5.65e-14 ***
    FALSE Residuals   19    530      28                     
    FALSE ---
    FALSE Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov1, which = "Region")

    FALSE   Tukey multiple comparisons of means
    FALSE     95% family-wise confidence level
    FALSE 
    FALSE Fit: aov(formula = PC1 ~ Region, data = pcadata)
    FALSE 
    FALSE $Region
    FALSE              diff       lwr      upr    p adj
    FALSE CA3-CA1  5.100214 -2.548692 12.74912 0.233216
    FALSE DG-CA1  50.510511 43.990986 57.03003 0.000000
    FALSE DG-CA3  45.410296 37.926612 52.89398 0.000000

    aov2 <- aov(PC2 ~ Region, data=pcadata)
    summary(aov2) 

    FALSE             Df Sum Sq Mean Sq F value   Pr(>F)    
    FALSE Region       2   4255  2127.4   255.3 1.86e-14 ***
    FALSE Residuals   19    158     8.3                     
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
    FALSE CA3-CA1  37.09928  32.918856  41.27970 0e+00
    FALSE DG-CA1   12.33263   8.769462  15.89581 1e-07
    FALSE DG-CA3  -24.76665 -28.856769 -20.67652 0e+00

    aov3 <- aov(PC3 ~ Treatment, data=pcadata)
    summary(aov3) 

    FALSE             Df Sum Sq Mean Sq F value Pr(>F)  
    FALSE Treatment    1  404.2   404.2   7.451 0.0129 *
    FALSE Residuals   20 1085.0    54.2                 
    FALSE ---
    FALSE Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    aov4 <- aov(PC4 ~ Treatment, data=pcadata)
    summary(aov4) 

    FALSE             Df Sum Sq Mean Sq F value  Pr(>F)   
    FALSE Treatment    1  324.1   324.1   10.11 0.00472 **
    FALSE Residuals   20  641.5    32.1                   
    FALSE ---
    FALSE Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

The gene expression data were exported to csv files for importing into
the GOMMU analysis package for subsequent analysis.

    # from https://github.com/rachelwright8/Ahya-White-Syndromes/blob/master/deseq2_Ahya.R
    res <- results(dds, contrast=c('Treatment', 'trained', 'yoked'), independentFiltering = FALSE)
    table(res$padj<0.05)

    ## 
    ## FALSE  TRUE 
    ## 16870   423

    logs <- data.frame(cbind("gene"=row.names(res),"logP"=round(-log(res$pvalue+1e-10,10),1)))
    logs$logP <- as.numeric(as.character(logs$logP))
    sign <- rep(1,nrow(logs))
    sign[res$log2FoldChange<0]=-1  ##change to correct model
    table(sign)

    ## sign
    ##   -1    1 
    ## 7759 9561

    logs$logP <- logs$logP*sign

    write.csv(logs, file = "./06_GO_MWU/03_behavior_GOpvals.csv", row.names = F)

Supplementary figures showing the distibution of pvalues.

    myhistogram(contrastvector = c('Region', 'CA1', 'DG'), mypval = 0.05)

![](../figures/03_cognitiontest/histogram-1.png)

    ## [1] 1

    myhistogram(contrastvector = c('Region', 'CA3', 'DG'), mypval = 0.05)

![](../figures/03_cognitiontest/histogram-2.png)

    ## [1] 1

    myhistogram(contrastvector = c('Region', 'CA1', 'CA3'), mypval = 0.05)

![](../figures/03_cognitiontest/histogram-3.png)

    ## [1] 1

    myhistogram(contrastvector = c('Treatment', 'trained', 'yoked'), mypval = 0.05)

![](../figures/03_cognitiontest/histogram-4.png)

    ## [1] 1

Supplementary figures of showing PC3 and PC4 contrasted against PC2.)

    ## PC2 vs PC3
    A <- plotPC2PC3(aescolor = pcadata$Region, colorname = "Region", aesshape = pcadata$Treatment, shapename = "Treatment", colorvalues = colorvalRegion)

    B <- plotPC2PC4(aescolor = pcadata$Region, colorname = "Region", aesshape = pcadata$Treatment, shapename = "Treatment", colorvalues = colorvalRegion)

    plot_grid(A, B, rel_widths = c(2.25,3))

![](../figures/03_cognitiontest/PCA34-1.png)
