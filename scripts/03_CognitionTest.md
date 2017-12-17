Examining of cognitive training on hippocampal
----------------------------------------------

The goals of the subsequent analysis are 1) to determine the effects of
cognitiving training on hippocampal gene expression and 2) related any
detectable changes to variation cause by other technical and biological
treatements.

The sample and count information for this part is found in
`../data/GSE100225_IntegrativeWT2015ColData.csv` and
`../data/GSE100225_IntegrativeWT2015CountData.csv`. You can also
download these two files (with a different name but same content) from
[GEO
GSE100225](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100225).

### Experimental Design

We use 3-4–month-old male C57BL/6J mice fro the Jackson Laboratory and
housed at the Marine Biological Laboratory. Mice (N=4) trained in the
active place avoidance task are conditioned to avoid mild shocks that
can be localized by visual cues in the enviornment. Yoked control mice
(N=4) are delivered sequence of unavoidable shock that mimickes the time
series of shocks received by the trained mice. While the trained and
yoked animals received the same number of shocks, only the trained
animals exhibitied an avoidance response.

Thirty minutes after the last cognitive training session, mice were
killed and transverse brain slices were prepared. The DG, CA3, CA1
subregions were microdissected using a 0.25 mm punch (Electron
Microscopy Systems) and a dissecting scope (Zeiss). RNA was isolated
using the Maxwell 16 LEV RNA Isolation Kit (Promega). RNA libraries were
prepared by the Genomic Sequencing and Analysis Facility at the
University of Texas at Austin using the Illumina HiSeq platform.

The orginal design was 4 animals per treament and 3 hippocampal sub
regions per animals, which would give 24 samples. After excluding
compromized samples, the final sample sizes are:

    ##    Treatment  Subfield
    ##  yoked  : 9   CA1:8   
    ##  trained:13   CA3:5   
    ##               DG :9

    ##          
    ##           CA1 CA3 DG
    ##   yoked     2   3  4
    ##   trained   6   2  5

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
`Treatment + Subfield + Treatment * Subfield`. Genes with less than 2
counts across all samples were filtered, leaving us with `dim(rld)`
genes for analysis of differntial expression.

We identified 423 genes were differentially expressed between the yoked
control and cognitively trained animals, 3485 genes that were
differentially expressed across subfields, and 324 showed an interaction
at FDR p &lt; 0.05 (Fig. 4B). We see a large effect of brain region on
gene expression, with 20% of detectable genes begin differentially
expressed between one or more brain-region comparisons (3485
differentially expressed genes /17320 measured genes). This is an order
of magnitude greater than the 2% of the transcriptome that changed in
response to learning (423 DEGs /17320 genes measured).

    res <- results(dds, contrast =c('Treatment', 'trained', 'yoked'), independentFiltering = T, alpha = 0.1)
    summary(res)

    ## 
    ## out of 17320 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)     : 838, 4.8% 
    ## LFC < 0 (down)   : 445, 2.6% 
    ## outliers [1]     : 27, 0.16% 
    ## low counts [2]   : 4689, 27% 
    ## (mean count < 5)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    table(res$padj<0.1)

    ## 
    ## FALSE  TRUE 
    ## 11321  1283

    head((res[order(res$padj),]), 10)

    ## log2 fold change (MLE): Treatment trained vs yoked 
    ## Wald test p-value: Treatment trained vs yoked 
    ## DataFrame with 10 rows and 6 columns
    ##         baseMean log2FoldChange     lfcSE      stat       pvalue
    ##        <numeric>      <numeric> <numeric> <numeric>    <numeric>
    ## Agap1  302.59405       2.891335 0.4286547  6.745137 1.528823e-11
    ## Mga    179.38078       2.979395 0.4665414  6.386133 1.701330e-10
    ## Sdhaf2  82.53733      -2.036752 0.3469272 -5.870835 4.336062e-09
    ## Ncoa4  144.40068       2.432693 0.4637984  5.245153 1.561527e-07
    ## Plxna4 965.33767       1.546712 0.2963598  5.219033 1.798594e-07
    ## Lhfpl4 283.85299       1.755736 0.3508599  5.004095 5.612521e-07
    ## Sdc3   289.35711       1.644503 0.3274979  5.021415 5.129225e-07
    ## Ttll7  418.67334       1.573042 0.3181860  4.943781 7.662192e-07
    ## Pcdh15  34.94756      -4.111590 0.8375541 -4.909044 9.152135e-07
    ## Pdxk   263.55895       1.619988 0.3377568  4.796316 1.616104e-06
    ##                padj
    ##           <numeric>
    ## Agap1  1.926929e-07
    ## Mga    1.072178e-06
    ## Sdhaf2 1.821724e-05
    ## Ncoa4  4.533897e-04
    ## Plxna4 4.533897e-04
    ## Lhfpl4 1.010575e-03
    ## Sdc3   1.010575e-03
    ## Ttll7  1.207178e-03
    ## Pcdh15 1.281706e-03
    ## Pdxk   2.036938e-03

    results <- data.frame(cbind("gene"=row.names(res), 
                             "baseMean" = res$baseMean,
                             "log2FoldChange" = res$log2FoldChange,
                             "lfcSE" = res$lfcSE,
                             "pvalue" = res$pvalue, "padj" = res$padj,
                             "logP"=round(-log(res$pvalue+1e-10,10),1)))

    write.csv(results, file = "../results/03_cognition_results.csv", row.names = F)

<table>
<thead>
<tr class="header">
<th>Contrast</th>
<th>Number of DEGs</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>CA1 vs. DG</td>
<td>2820</td>
</tr>
<tr class="even">
<td>CA3 vs. DG</td>
<td>3014</td>
</tr>
<tr class="odd">
<td>CA1 vs. CA3</td>
<td>2293</td>
</tr>
<tr class="even">
<td>CA1 vs. DG</td>
<td>1283</td>
</tr>
</tbody>
</table>

![](../figures/03_cognitiontest/VennDiagramPadj-1.png)

    candidates <- list("CA1 vs. DG" = venn1, "CA3 vs. DG" = venn2, "CA1 vs. CA3" = venn3)

    prettyvenn <- venn.diagram(
      scaled=T,
      x = candidates, filename=NULL, 
      col = "black",
      fill = c( "white", "white", "white"),
      alpha = 0.5,
      cex = 1, fontfamily = "sans", #fontface = "bold",
      cat.default.pos = "text",
      cat.dist = c(0.07, 0.07, 0.07), cat.pos = 1,
      cat.cex = 1, cat.fontfamily = "sans")
    #dev.off()
    grid.draw(prettyvenn)

![](../figures/03_cognitiontest/VennDiagramPadj2-1.png)

Hierarchical clustering of the differentially expressed genes separates
samples by both subfield and treatment.

Then, I visuazlied the data as a heatmap showing the relative log fold
change of gene expression across samples. Genes were filtered for a
minimimum adjust p value &lt; 0.05 in any two-way contrast. The row mean
for each gene was subtracted for the raw value to allow for analysis of
fold change rather than raw magnitudes. The samples cluster primarily by
brain region with some small treatment-driven.

![](../figures/03_cognitiontest/HeatmapPadj-1.png)

volcano plots yea!
==================

    res <- results(dds, contrast =c('Treatment', 'trained', 'yoked'), independentFiltering = T, alpha = 0.1)
    summary(res)

    ## 
    ## out of 17320 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)     : 838, 4.8% 
    ## LFC < 0 (down)   : 445, 2.6% 
    ## outliers [1]     : 27, 0.16% 
    ## low counts [2]   : 4689, 27% 
    ## (mean count < 5)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    resOrdered <- res[order(res$padj),]
    head(resOrdered, 3)

    ## log2 fold change (MLE): Treatment trained vs yoked 
    ## Wald test p-value: Treatment trained vs yoked 
    ## DataFrame with 3 rows and 6 columns
    ##         baseMean log2FoldChange     lfcSE      stat       pvalue
    ##        <numeric>      <numeric> <numeric> <numeric>    <numeric>
    ## Agap1  302.59405       2.891335 0.4286547  6.745137 1.528823e-11
    ## Mga    179.38078       2.979395 0.4665414  6.386133 1.701330e-10
    ## Sdhaf2  82.53733      -2.036752 0.3469272 -5.870835 4.336062e-09
    ##                padj
    ##           <numeric>
    ## Agap1  1.926929e-07
    ## Mga    1.072178e-06
    ## Sdhaf2 1.821724e-05

    data <- data.frame(gene = row.names(res),
                       pvalue = -log10(res$padj), 
                       lfc = res$log2FoldChange)
    data <- na.omit(data)

    data <- data %>%
      mutate(color = ifelse(data$lfc > 0 & data$pvalue > 1, 
                            yes = "trained", 
                            no = ifelse(data$lfc < 0 & data$pvalue > 1, 
                                        yes = "yoked", 
                                        no = "none")))
    data$color <- as.factor(data$color)
    summary(data)

    ##             gene           pvalue              lfc         
    ##  0610007P14Rik:    1   Min.   :0.000137   Min.   :-9.0497  
    ##  0610009B22Rik:    1   1st Qu.:0.098852   1st Qu.:-0.3979  
    ##  0610009O20Rik:    1   Median :0.287291   Median : 0.1503  
    ##  0610010F05Rik:    1   Mean   :0.423241   Mean   : 0.3640  
    ##  0610010K14Rik:    1   3rd Qu.:0.626286   3rd Qu.: 0.7999  
    ##  0610012G03Rik:    1   Max.   :6.715134   Max.   : 6.7936  
    ##  (Other)      :12598                                       
    ##      color      
    ##  none   :11321  
    ##  trained:  838  
    ##  yoked  :  445  
    ##                 
    ##                 
    ##                 
    ## 

    write.csv(data, "../results/03_cognition_volcanoTreatment.csv")



    res <- results(dds, contrast =c("Subfield", "CA1", "DG"), independentFiltering = T, alpha = 0.05)
    resOrdered <- res[order(res$padj),]
    head(resOrdered, 3)

    ## log2 fold change (MLE): Subfield CA1 vs DG 
    ## Wald test p-value: Subfield CA1 vs DG 
    ## DataFrame with 3 rows and 6 columns
    ##         baseMean log2FoldChange     lfcSE      stat       pvalue
    ##        <numeric>      <numeric> <numeric> <numeric>    <numeric>
    ## Pou3f1  295.1569       5.938131 0.3913099  15.17501 5.177060e-52
    ## Prkcg  1817.7550       3.023026 0.2082006  14.51977 9.081344e-48
    ## Wfs1    775.6657       6.372533 0.4522933  14.08938 4.414100e-45
    ##                padj
    ##           <numeric>
    ## Pou3f1 6.525166e-48
    ## Prkcg  5.723063e-44
    ## Wfs1   1.854511e-41

    data <- data.frame(gene = row.names(res), pvalue = -log10(res$padj), lfc = res$log2FoldChange)
    data <- na.omit(data)
    data <- data %>%
      mutate(color = ifelse(data$lfc > 0 & data$pvalue > 1, 
                            yes = "CA1", 
                            no = ifelse(data$lfc < 0 & data$pvalue > 1, 
                                        yes = "DG", 
                                        no = "none")))
    data$color <- as.factor(data$color)
    summary(data)

    ##             gene           pvalue              lfc            color     
    ##  0610007P14Rik:    1   Min.   : 0.00002   Min.   :-8.84062   CA1 :1289  
    ##  0610009B22Rik:    1   1st Qu.: 0.13324   1st Qu.:-0.92526   DG  :1546  
    ##  0610009O20Rik:    1   Median : 0.38749   Median :-0.08317   none:9769  
    ##  0610010F05Rik:    1   Mean   : 0.84454   Mean   :-0.29729              
    ##  0610010K14Rik:    1   3rd Qu.: 0.91357   3rd Qu.: 0.59402              
    ##  0610012G03Rik:    1   Max.   :47.18541   Max.   : 8.20228              
    ##  (Other)      :12598

    write.csv(data, "../results/03_cognition_volcanoCA1DG.csv")

A principal component analysis of all gene expression data revealed that
brain region explains 75% of the variance in the data. PC1 accounts for
56% of the variance and distinguishes DG from not-DG samples (ANOVA for
PC1 ~ Subfield: F2,19= 226.1; p &lt; 0.001). A post hoc Tukey test
showed that DG samples are significantly different from both CA1 and CA3
samples (CA1-DG, p &lt; 0.001; CA3-DG, p &lt; 0.001; CA1-CA3, p = 0.23).
PC2 accounts for 19% of the variance and distinguishes the three
subfields (PC2 ~ Subfield ANOVA: F2,19= 255.3; p &lt; 0.001; Tukey test,
p &lt; 0.001for all three comparisons). PC3 are PC4 are influenced by
variation due to cognitive training (PC3 ~ Treatment ANOVA: F1,20=7.451;
p = 0.012, PC4 ~ Treatment ANOVA: F1,20=10.11; p = 0.0047). An analysis
of Gene Ontology identified 91 GO terms at 10% FDR significant. Among
the top 10 are glutamate signaling and membrane transport systems and a
downregulation of oxidoreductase and ribosomal activity (Fig. 2C).

    ## statistics
    aov1R <- aov(PC1 ~ Subfield, data=pcadata)
    summary(aov1R) 

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Subfield     2  12615    6307   226.1 5.65e-14 ***
    ## Residuals   19    530      28                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov1R, which = "Subfield") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ Subfield, data = pcadata)
    ## 
    ## $Subfield
    ##              diff       lwr      upr    p adj
    ## CA3-CA1  5.100214 -2.548692 12.74912 0.233216
    ## DG-CA1  50.510511 43.990986 57.03003 0.000000
    ## DG-CA3  45.410296 37.926612 52.89398 0.000000

    aov2R <- aov(PC2 ~ Subfield, data=pcadata)
    summary(aov2R) 

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Subfield     2   4255  2127.4   255.3 1.86e-14 ***
    ## Residuals   19    158     8.3                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov2R, which = "Subfield") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ Subfield, data = pcadata)
    ## 
    ## $Subfield
    ##              diff        lwr       upr p adj
    ## CA3-CA1  37.09928  32.918856  41.27970 0e+00
    ## DG-CA1   12.33263   8.769462  15.89581 1e-07
    ## DG-CA3  -24.76665 -28.856769 -20.67652 0e+00

    aov3R <- aov(PC3 ~ Subfield, data=pcadata)
    summary(aov3R) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Subfield     2   34.4   17.21   0.225  0.801
    ## Residuals   19 1454.8   76.57

    TukeyHSD(aov3R, which = "Subfield") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC3 ~ Subfield, data = pcadata)
    ## 
    ## $Subfield
    ##               diff       lwr       upr     p adj
    ## CA3-CA1  2.4937376 -10.17909 15.166567 0.8722208
    ## DG-CA1  -0.7357084 -11.53736 10.065942 0.9836439
    ## DG-CA3  -3.2294460 -15.62853  9.169639 0.7880961

    aov4R <- aov(PC3 ~ Subfield, data=pcadata)
    summary(aov4R) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Subfield     2   34.4   17.21   0.225  0.801
    ## Residuals   19 1454.8   76.57

    TukeyHSD(aov4R, which = "Subfield") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC3 ~ Subfield, data = pcadata)
    ## 
    ## $Subfield
    ##               diff       lwr       upr     p adj
    ## CA3-CA1  2.4937376 -10.17909 15.166567 0.8722208
    ## DG-CA1  -0.7357084 -11.53736 10.065942 0.9836439
    ## DG-CA3  -3.2294460 -15.62853  9.169639 0.7880961

    aov5R <- aov(PC5 ~ Subfield, data=pcadata)
    summary(aov5R) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Subfield     2    8.0   4.014   0.154  0.858
    ## Residuals   19  495.5  26.082

    TukeyHSD(aov5R, which = "Subfield") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC5 ~ Subfield, data = pcadata)
    ## 
    ## $Subfield
    ##               diff       lwr      upr     p adj
    ## CA3-CA1 -0.2658585 -7.662234 7.130517 0.9954145
    ## DG-CA1  -1.3140576 -7.618338 4.990222 0.8579143
    ## DG-CA3  -1.0481991 -8.284807 6.188409 0.9283575

    aov6R <- aov(PC6 ~ Subfield, data=pcadata)
    summary(aov6R) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Subfield     2    1.1   0.549   0.027  0.973
    ## Residuals   19  379.9  19.993

    TukeyHSD(aov6R, which = "Subfield") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC6 ~ Subfield, data = pcadata)
    ## 
    ## $Subfield
    ##                diff       lwr      upr     p adj
    ## CA3-CA1 -0.51051373 -6.986232 5.965204 0.9781557
    ## DG-CA1   0.03953805 -5.480022 5.559098 0.9998174
    ## DG-CA3   0.55005178 -5.785786 6.885889 0.9735779

    aov1T <- aov(PC1 ~ Treatment, data=pcadata)
    summary(aov1T) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment    1     28    27.8   0.042  0.839
    ## Residuals   20  13117   655.9

    aov2T <- aov(PC2 ~ Treatment, data=pcadata)
    summary(aov2T) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## Treatment    1    599   599.2   3.142 0.0915 .
    ## Residuals   20   3814   190.7                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    aov3T <- aov(PC3 ~ Treatment, data=pcadata)
    summary(aov3T) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## Treatment    1  404.2   404.2   7.451 0.0129 *
    ## Residuals   20 1085.0    54.2                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    aov4T <- aov(PC4 ~ Treatment, data=pcadata)
    summary(aov4T) 

    ##             Df Sum Sq Mean Sq F value  Pr(>F)   
    ## Treatment    1  324.1   324.1   10.11 0.00472 **
    ## Residuals   20  641.5    32.1                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    aov5T <- aov(PC5 ~ Treatment, data=pcadata)
    summary(aov5T) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment    1    2.7   2.743    0.11  0.744
    ## Residuals   20  500.8  25.042

    aov6T <- aov(PC6 ~ Treatment, data=pcadata)
    summary(aov6T) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment    1   25.7   25.66   1.444  0.243
    ## Residuals   20  355.3   17.77

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

    myhistogram(contrastvector = c('Subfield', 'CA1', 'DG'), mypval = 0.1)

![](../figures/03_cognitiontest/histogram-1.png)

    ## [1] 1

    myhistogram(contrastvector = c('Subfield', 'CA3', 'DG'), mypval = 0.1)

![](../figures/03_cognitiontest/histogram-2.png)

    ## [1] 1

    myhistogram(contrastvector = c('Subfield', 'CA1', 'CA3'), mypval = 0.1)

![](../figures/03_cognitiontest/histogram-3.png)

    ## [1] 1

    myhistogram(contrastvector = c('Treatment', 'trained', 'yoked'), mypval = 0.1)

![](../figures/03_cognitiontest/histogram-4.png)

    ## [1] 1

Here is the corresponding Adobe Illustrator file that combines many of
the above plots.

<img src="../figures/fig_03-memory.png" width="1370" />
