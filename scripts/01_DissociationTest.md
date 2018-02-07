### Identifying the effects of cellular dissociation on hippocampal transcriptomes

The sample and count information for this part is found in
`../data/GSE99765_DissociationColData.csv` and
`../data/GSE99765_DissociationCountData.csv`. You can also download
these two files (with a different name but same content) from [GEO
GSE99765](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99765).

    colData <- read.csv('../data/GSE99765_DissociationColData.csv')
    rownames(colData) <- colData$RNAseqID
    countData <-  read.csv('../data/GSE99765_DissociationCountData.csv', check.names = F, row.names = 1)

Sample sizes

    colData <- rename(colData, c("Region"="Subfield"))
    table(colData$Treatment,colData$Subfield) 

    ##              
    ##               CA1 CA3 DG
    ##   control       3   2  2
    ##   dissociated   3   2  2

    dim(countData)

    ## [1] 22485    14

    write.csv(colData, "../results/01_dissociation_colData.csv", row.names = F)
    write.csv(countData, "../results/01_dissociation_countData.csv", row.names = T)

I used DESeq2 (Love et al., 2014) for gene expression normalization and
quantification using the following experimental design:
`Treatment + Subfield + Treatment * Subfield`. Genes with less than 2
counts across all samples were filtered, leaving us with `dim(rld)`
number of genes for analysis of differntial expression.

    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ Treatment + Subfield + Treatment * Subfield )
    dds <- dds[ rowSums(counts(dds)) > 2, ] ## filter genes 
    dds <- DESeq(dds) # Differential expression analysis

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    rld <- rlog(dds, blind=FALSE) ## log transformed data
    dim(rld) #print total genes analyzed

    ## [1] 16709    14

    vsd <- vst(dds, blind=FALSE) # variance stabilized
    head(assay(rld), 3)

    ##               100-CA1-1 100-CA1-2 100-CA1-3 100-CA3-1 100-CA3-4 100-DG-2
    ## 0610007P14Rik  4.588506  4.776457  4.853057  5.079028  5.171328 5.186172
    ## 0610009B22Rik  3.186440  3.699914  3.130243  3.602998  3.097215 3.642540
    ## 0610009L18Rik  1.776934  2.360120  1.673957  2.380223  2.128380 1.829261
    ##               100-DG-3 101-CA1-1 101-CA1-2 101-CA1-3 101-CA3-1 101-CA3-4
    ## 0610007P14Rik 5.030617  4.955175  4.217261  3.466905  5.086613  5.019213
    ## 0610009B22Rik 3.482691  3.304835  4.260577  2.711435  3.820068  3.477413
    ## 0610009L18Rik 2.487685  2.039003  1.952140  1.795752  2.240397  1.916907
    ##               101-DG-3 101-DG-4
    ## 0610007P14Rik 5.439691 4.155567
    ## 0610009B22Rik 4.200692 2.744500
    ## 0610009L18Rik 3.161293 3.000009

    head(assay(vsd), 3)

    ##               100-CA1-1 100-CA1-2 100-CA1-3 100-CA3-1 100-CA3-4 100-DG-2
    ## 0610007P14Rik  6.024282  6.147172  6.198913  6.358510  6.425746 6.439969
    ## 0610009B22Rik  5.500732  5.801174  5.469570  5.744608  5.451894 5.771328
    ## 0610009L18Rik  5.099770  5.428706  5.033775  5.444353  5.302348 5.116703
    ##               100-DG-3 101-CA1-1 101-CA1-2 101-CA1-3 101-CA3-1 101-CA3-4
    ## 0610007P14Rik 6.321325  6.268696  5.582529  4.748434  6.363481  6.318356
    ## 0610009B22Rik 5.671124  5.569662  6.355559  4.748434  5.880938  5.669040
    ## 0610009L18Rik 5.499883  5.255534  4.748434  4.748434  5.364748  5.165706
    ##               101-DG-3 101-DG-4
    ## 0610007P14Rik 6.872924 5.722956
    ## 0610009B22Rik 6.528940 5.152565
    ## 0610009L18Rik 6.528940 5.867002

    write.csv(assay(vsd), "../results/01_dissociation_vsd.csv")
    write.csv(assay(rld), "../results/01_dissociation_rld.csv")

We identified 162 genes that were differentially expressed between the
control and dissociated samples, 331 genes that were differentially
expressed genes (DEGs) between any of the three hippocampus subfields,
and 30 genes were shared between both sets of differentially expressed
genes at FDR p-value &lt; 0.05 (Fig 1B).

    ## DEG by contrasts
    source("resvalsfunction.R")
    contrast1 <- resvals(contrastvector = c('Subfield', 'CA1', 'DG'), mypval = 0.1)

    ## [1] 484

    contrast2 <- resvals(contrastvector = c('Subfield', 'CA3', 'DG'), mypval = 0.1)

    ## [1] 98

    contrast3 <- resvals(contrastvector = c('Subfield', 'CA1', 'CA3'), mypval = 0.1)

    ## [1] 18

    contrast4 <- resvals(contrastvector = c('Treatment', 'dissociated', 'control'), mypval = 0.1)

    ## [1] 344

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

    ## check order for correctness
    candidates <- list("Subfield" = venn123, "Treatment" = venn4)

    prettyvenn <- venn.diagram(
      scaled=T,
      x = candidates, filename=NULL, 
      col = "black",
      fill = c( "white", "white"),
      alpha = 0.5,
      cex = 1, fontfamily = "sans", #fontface = "bold",
      cat.default.pos = "text",
      cat.dist = c(0.08, 0.08), cat.pos = 1,
      cat.cex = 1, cat.fontfamily = "sans")
    #dev.off()
    grid.draw(prettyvenn)

![](../figures/01_dissociationtest/VennDiagramPadj-1.png)

    # save files for metanalysis
    write(venn123, "../results/01_dissociation_venn123.txt")
    write(venn4, "../results/01_dissociation_venn4.txt")
    write(venn1, "../results/01_dissociation_venn1.txt")

A hierarchical clustering analysis of all differentially expressed genes
does not give rise to distinct clusters that are separated by subfield
or method; however, when examining the control, homogenized samples
alone (identified with light grey boxes), the three subfields form
distinct clusters, while the dissociated samples do not cluster by
subfield (Fig. 1C).

    contrast4 <- resvals(contrastvector = c('Treatment', 'dissociated', 'control'), mypval = 0.001)

    ## [1] 32

    DEGes <- assay(rld)
    DEGes <- cbind(DEGes, contrast4)
    DEGes <- as.data.frame(DEGes) # convert matrix to dataframe
    DEGes$rownames <- rownames(DEGes)  # add the rownames to the dataframe

    DEGes$padjmin <- with(DEGes, pmin(padjTreatmentdissociatedcontrol)) # put the min pvalue in a new column

    write.csv(as.data.frame(DEGes), "../results/01_dissociation_DEGes.csv", row.names = F)

volcano plots yea!
==================

    res <- results(dds, contrast =c('Treatment', 'dissociated', 'control'), independentFiltering = T, alpha = 0.1)
    summary(res)

    ## 
    ## out of 16709 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)     : 288, 1.7% 
    ## LFC < 0 (down)   : 56, 0.34% 
    ## outliers [1]     : 18, 0.11% 
    ## low counts [2]   : 4534, 27% 
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    resOrdered <- res[order(res$padj),]
    head(resOrdered, 3)

    ## log2 fold change (MLE): Treatment dissociated vs control 
    ## Wald test p-value: Treatment dissociated vs control 
    ## DataFrame with 3 rows and 6 columns
    ##         baseMean log2FoldChange     lfcSE      stat       pvalue
    ##        <numeric>      <numeric> <numeric> <numeric>    <numeric>
    ## Trf    434.21998       2.724762 0.4133952  6.591178 4.363497e-11
    ## Hexb   218.64390       2.348231 0.3655736  6.423415 1.332509e-10
    ## Selplg  69.25758       2.969443 0.4682601  6.341439 2.276293e-10
    ##                padj
    ##           <numeric>
    ## Trf    5.304704e-07
    ## Hexb   8.099655e-07
    ## Selplg 9.224297e-07

    data <- data.frame(gene = row.names(res),
                       pvalue = -log10(res$padj), 
                       lfc = res$log2FoldChange)
    data <- na.omit(data)
    data <- data %>%
      mutate(color = ifelse(data$lfc > 0 & data$pvalue > 1, 
                            yes = "dissociated", 
                            no = ifelse(data$lfc < 0 & data$pvalue > 1, 
                                        yes = "control", 
                                        no = "none")))
    data$color <- as.factor(data$color)
    summary(data)

    ##             gene           pvalue              lfc          
    ##  0610007P14Rik:    1   Min.   :0.000003   Min.   :-5.06587  
    ##  0610009B22Rik:    1   1st Qu.:0.026322   1st Qu.:-0.40829  
    ##  0610009L18Rik:    1   Median :0.073082   Median : 0.04495  
    ##  0610009O20Rik:    1   Mean   :0.185069   Mean   : 0.15206  
    ##  0610010F05Rik:    1   3rd Qu.:0.204863   3rd Qu.: 0.56869  
    ##  0610010K14Rik:    1   Max.   :6.275339   Max.   : 9.47422  
    ##  (Other)      :12151                                        
    ##          color      
    ##  control    :   56  
    ##  dissociated:  288  
    ##  none       :11813  
    ##                     
    ##                     
    ##                     
    ## 

    write.csv(data, "../results/01_dissociation_volcanoTreatment.csv")

    res <- results(dds, contrast =c("Subfield", "CA1", "DG"), independentFiltering = T, alpha = 0.05)
    summary(res)

    ## 
    ## out of 16709 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)     : 194, 1.2% 
    ## LFC < 0 (down)   : 155, 0.93% 
    ## outliers [1]     : 18, 0.11% 
    ## low counts [2]   : 4534, 27% 
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    resOrdered <- res[order(res$padj),]
    head(resOrdered, 3)

    ## log2 fold change (MLE): Subfield CA1 vs DG 
    ## Wald test p-value: Subfield CA1 vs DG 
    ## DataFrame with 3 rows and 6 columns
    ##         baseMean log2FoldChange     lfcSE      stat       pvalue
    ##        <numeric>      <numeric> <numeric> <numeric>    <numeric>
    ## C1ql2  130.70538      -7.653838 0.7860552 -9.737024 2.096033e-22
    ## Stxbp6 143.43881      -5.193969 0.5589022 -9.293163 1.497702e-20
    ## Crlf1   40.19316      -7.699695 0.8435166 -9.128090 6.971829e-20
    ##                padj
    ##           <numeric>
    ## C1ql2  2.548147e-18
    ## Stxbp6 9.103783e-17
    ## Crlf1  2.825218e-16

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

    ##             gene           pvalue               lfc           color      
    ##  0610007P14Rik:    1   Min.   : 0.000002   Min.   :-9.3376   CA1 :  264  
    ##  0610009B22Rik:    1   1st Qu.: 0.004630   1st Qu.:-0.5435   DG  :  222  
    ##  0610009L18Rik:    1   Median : 0.008930   Median :-0.1267   none:11671  
    ##  0610009O20Rik:    1   Mean   : 0.161534   Mean   :-0.1334               
    ##  0610010F05Rik:    1   3rd Qu.: 0.058317   3rd Qu.: 0.2929               
    ##  0610010K14Rik:    1   Max.   :17.593776   Max.   : 8.4434               
    ##  (Other)      :12151

    write.csv(data, "../results/01_dissociation_volcanoCA1DG.csv")

This PCA gives an overview of the variability between samples using the
a large matrix of log transformed gene expression data. You can see that
the biggest difference is between DG punches and the CA1 and CA3
punches. CA1 and CA3 samples have similar transcriptomes. The control
CA1 samples have the most similar transcriptonal profiles as evidenced
by their tight clustering.

    colorvalSubfield <- c("#7570b3", "#1b9e77", "#d95f02")
    colorvalTreatment <- c("#ffffff", "#525252")


    # create the dataframe using my function pcadataframe
    pcadata <- pcadataframe(rld, intgroup=c("Subfield", "Treatment"), returnData=TRUE)
    percentVar <- round(100 * attr(pcadata, "percentVar"))

    PCA12 <- ggplot(pcadata, aes(PC1, PC2, shape = Treatment, color = Subfield)) + 
      geom_point(size = 3, alpha = 1) +
        xlab(paste0("PC1: ", percentVar[1],"% variance")) +
        ylab(paste0("PC2: ", percentVar[2],"% variance")) +
        scale_color_manual(values = colorvalSubfield) +
        theme_cowplot(font_size = 8, line_size = 0.25)  +
        theme(legend.position="none") +
        scale_shape_manual(values=c(16, 1)) 
    PCA12

![](../figures/01_dissociationtest/PCA-1.png)

    pdf(file="../figures/01_dissociationtest/PCA-1.pdf", width=1.75, height=2)
    plot(PCA12)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    ## statistics
    aov1 <- aov(PC1 ~ Subfield, data=pcadata)
    summary(aov1) 

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Subfield     2 2812.7  1406.4   17.69 0.000365 ***
    ## Residuals   11  874.3    79.5                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov1, which = "Subfield") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ Subfield, data = pcadata)
    ## 
    ## $Subfield
    ##              diff       lwr      upr     p adj
    ## CA3-CA1  5.223942 -10.31899 20.76687 0.6467956
    ## DG-CA1  33.098083  17.55515 48.64101 0.0003454
    ## DG-CA3  27.874142  10.84772 44.90057 0.0027013

    aov2 <- aov(PC2 ~ Subfield, data=pcadata)
    summary(aov2) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Subfield     2  243.8   121.9   0.744  0.498
    ## Residuals   11 1801.4   163.8

    TukeyHSD(aov2, which = "Subfield") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ Subfield, data = pcadata)
    ## 
    ## $Subfield
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

Next, save files for dowstream GO analysis.

    # from https://github.com/rachelwright8/Ahya-White-Syndromes/blob/master/deseq2_Ahya.R

    resCD=results(dds, contrast=c('Treatment', 'dissociated', 'control'), independentFiltering = F)
    table(resCD$padj<0.05)

    ## 
    ## FALSE  TRUE 
    ## 16529   162

    logs <- data.frame(cbind("gene"=row.names(resCD),"logP"=round(-log(resCD$pvalue+1e-10,10),1)))
    logs$logP=as.numeric(as.character(logs$logP))
    sign <- rep(1,nrow(logs))
    sign[resCD$log2FoldChange<0]=-1  ##change to correct model
    table(sign)

    ## sign
    ##   -1    1 
    ## 6989 9720

    logs$logP <- logs$logP*sign

    write.csv(logs, file = "./06_GO_MWU/01_dissociation_GOpvals.csv", row.names = F)

To view a histogram of the p-value distibution for each constrast,
change the Rmd file to `include=TRUE` for this chunck.

Here is the corresponding Adobe Illustrator file that combines many of
the above plots.

<img src="../figures/fig_01-dissociation.png" width="1370" />
