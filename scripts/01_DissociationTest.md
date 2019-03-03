### Identifying the effects of cellular dissociation on hippocampal transcriptomes

The sample and count information for this part is found in
`../data/GSE99765_DissociationColData.csv` and
`../data/GSE99765_DissociationCountData.csv`. You can also download
these two files from [GEO
GSE99765](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99765).

    colData <- read.csv('../data/GSE99765_DissociationColData.csv')
    rownames(colData) <- colData$RNAseqID
    countData <-  read.csv('../data/GSE99765_DissociationCountData.csv', check.names = F, row.names = 1)

Do a little data cleaning and calculate sample size and number of genes
measured.

    # rename column and samples. 
    colData <- rename(colData, c("Region"="Subfield"))
    colData$Treatment <- revalue(colData$Treatment, c("control"="HOMO", "dissociated"="DISS"))

    # calculate samples size and number of genes for which we have expression data
    table(colData$Treatment,colData$Subfield) 

    ##       
    ##        CA1 CA3 DG
    ##   HOMO   3   2  2
    ##   DISS   3   2  2

    dim(countData)

    ## [1] 22485    14

I used DESeq2 (Love et al., 2014) for gene expression normalization and
quantification using the following experimental design:
`Treatment + Subfield + Treatment * Subfield`. Genes with less than 2
counts across all samples were filtered, leaving us with `dim(rld)`
number of genes for analysis of differntial expression.

    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ Treatment + Subfield + Treatment * Subfield )
    dds <- dds[ rowSums(counts(dds)) > 2, ] ## pre-filter genes 
    dds <- DESeq(dds) # Differential expression analysis
    rld <- rlog(dds, blind=FALSE) ## log transformed data
    dim(rld) #print total genes analyzed

    ## [1] 16709    14

    vsd <- vst(dds, blind=FALSE) # variance stabilized
    head(assay(rld), 3)

    ##               100-CA1-1 100-CA1-2 100-CA1-3 100-CA3-1 100-CA3-4 100-DG-2
    ## 0610007P14Rik  4.588504  4.776456  4.853058  5.079031  5.171332 5.186176
    ## 0610009B22Rik  3.186433  3.699918  3.130234  3.603000  3.097205 3.642543
    ## 0610009L18Rik  1.776915  2.360122  1.673934  2.380225  2.128375 1.829245
    ##               100-DG-3 101-CA1-1 101-CA1-2 101-CA1-3 101-CA3-1 101-CA3-4
    ## 0610007P14Rik 5.030620  4.955177  4.217255  3.466887  5.086616  5.019215
    ## 0610009B22Rik 3.482690  3.304830  4.260588  2.711418  3.820074  3.477413
    ## 0610009L18Rik 2.487691  2.038995  1.952132  1.795738  2.240395  1.916895
    ##               101-DG-3 101-DG-4
    ## 0610007P14Rik 5.439697 4.155559
    ## 0610009B22Rik 4.200700 2.744482
    ## 0610009L18Rik 3.161307 3.000024

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

    # save results 
    write.csv(assay(vsd), "../results/vsd.csv")
    write.csv(assay(rld), "../results/rld.csv")

We identified 344 genes that were differentially expressed between the
homogenized and dissociated samples at FDR p-value &lt; 0.1.

    ## DEG by contrasts at 0.1 pvalue
    contrast1 <- resvals(contrastvector = c('Subfield', 'CA1', 'DG'), mypval = 0.1) #484

    ## [1] 484

    contrast2 <- resvals(contrastvector = c('Subfield', 'CA3', 'DG'), mypval = 0.1) #98

    ## [1] 98

    contrast3 <- resvals(contrastvector = c('Subfield', 'CA1', 'CA3'), mypval = 0.1) #18

    ## [1] 18

    contrast4 <- resvals(contrastvector = c('Treatment', 'DISS', 'HOMO'), mypval = 0.1) #344

    ## [1] 344

    # % transcrptiome altered by treatment
    344/16709 * 100

    ## [1] 2.058771

A hierarchical clustering analysis of all differentially expressed genes
does not give rise to distinct clusters that are separated by subfield
or method; however, when examining the control, homogenized samples
alone (identified with light grey boxes), the three subfields form
distinct clusters, while the dissociated samples do not cluster by
subfield.

    contrast4 <- resvals(contrastvector = c('Treatment', 'DISS', 'HOMO'), mypval = 0.01)

    ## [1] 67

    DEGs <- assay(rld)
    DEGs <- cbind(DEGs, contrast4)
    DEGs <- as.data.frame(DEGs) # convert matrix to dataframe
    DEGs$rownames <- rownames(DEGs)  # add the rownames to the dataframe

    DEGs$padjmin <- with(DEGs, pmin(padjTreatmentDISSHOMO)) # put the min pvalue in a new column

    write.csv(as.data.frame(DEGs), "../results/heatmap_DEGs.csv", row.names = F)
    write.csv(colData, "../results/heatmap_colData.csv", row.names = F) 

Volcano Plots
-------------

Craete new data frames that include fold change, pvalue, and a column
describing the direction for differential gene expression. This
“direction” will be used to color code the dots on the volcano plot.
Will also save a list of DEGs at the end.

    res <- results(dds, contrast =c('Treatment', 'DISS', 'HOMO'), independentFiltering = T, alpha = 0.01)
    summary(res)

    ## 
    ## out of 16709 with nonzero total read count
    ## adjusted p-value < 0.01
    ## LFC > 0 (up)       : 65, 0.39%
    ## LFC < 0 (down)     : 4, 0.024%
    ## outliers [1]       : 18, 0.11%
    ## low counts [2]     : 6474, 39%
    ## (mean count < 11)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    res <- results(dds, contrast =c('Treatment', 'DISS', 'HOMO'), independentFiltering = T, alpha = 0.1)
    summary(res)

    ## 
    ## out of 16709 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 288, 1.7%
    ## LFC < 0 (down)     : 56, 0.34%
    ## outliers [1]       : 18, 0.11%
    ## low counts [2]     : 4534, 27%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    288+56 # tolal number of DEGs = 344

    ## [1] 344

    (344/16709)*100 # percent of DEGs out of total measured

    ## [1] 2.058771

    data <- data.frame(gene = row.names(res),
                       pvalue = -log10(res$padj), 
                       lfc = res$log2FoldChange,
                       padj = res$padj)
    data <- na.omit(data)
    data <- data %>%
        mutate(direction = ifelse(data$lfc > 1 & data$pvalue > 1, 
                            yes = "DISS", 
                            no = ifelse(data$lfc < -1 & data$pvalue > 1, 
                                        yes = "HOMO", 
                                        no = "neither")))
    data$direction <- as.factor(data$direction)
    summary(data)

    ##             gene           pvalue              lfc          
    ##  0610007P14Rik:    1   Min.   :0.000003   Min.   :-5.06587  
    ##  0610009B22Rik:    1   1st Qu.:0.026323   1st Qu.:-0.40829  
    ##  0610009L18Rik:    1   Median :0.073080   Median : 0.04498  
    ##  0610009O20Rik:    1   Mean   :0.185088   Mean   : 0.15206  
    ##  0610010F05Rik:    1   3rd Qu.:0.204839   3rd Qu.: 0.56869  
    ##  0610010K14Rik:    1   Max.   :6.274558   Max.   : 9.47422  
    ##  (Other)      :12151                                        
    ##       padj             direction    
    ##  Min.   :0.0000005   DISS   :  274  
    ##  1st Qu.:0.6239665   HOMO   :   47  
    ##  Median :0.8451229   neither:11836  
    ##  Mean   :0.7499809                  
    ##  3rd Qu.:0.9411891                  
    ##  Max.   :0.9999928                  
    ## 

    # note, there are fewer DEGs now because they have been filtered for higher than lfc 1.5

    write.csv(data, "../results/volcanoTreatment.csv")

    # save the list of just DEGs with their pvalu and lfc and direction
    dissocDEGs <- data %>%
      filter(direction != "neither")
    dissocDEGs <- dissocDEGs[order(dissocDEGs$padj),]
    dissocDEGs$pvalue <- NULL

    dissocDEGs$padj <-formatC(dissocDEGs$padj, format = "e", digits = 2) # round to scientific
    dissocDEGs <- dissocDEGs %>% mutate_if(is.numeric, ~round(., 2)) # round to 2 dec

    head(dissocDEGs)

    ##     gene  lfc     padj direction
    ## 1    Trf 2.72 5.31e-07      DISS
    ## 2   Hexb 2.35 8.10e-07      DISS
    ## 3 Selplg 2.97 9.22e-07      DISS
    ## 4   C1qb 2.28 7.07e-06      DISS
    ## 5  Csf1r 2.13 9.58e-06      DISS
    ## 6   Ctss 2.59 9.58e-06      DISS

    write.csv(dissocDEGs, "../results/dissociationDEGs.csv", row.names = F)

    res <- results(dds, contrast =c("Subfield", "CA1", "DG"), independentFiltering = T, alpha = 0.1)
    summary(res)

    ## 
    ## out of 16709 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 262, 1.6%
    ## LFC < 0 (down)     : 222, 1.3%
    ## outliers [1]       : 18, 0.11%
    ## low counts [2]     : 4210, 25%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    resOrdered <- res[order(res$padj),]
    head(resOrdered, 3)

    ## log2 fold change (MLE): Subfield CA1 vs DG 
    ## Wald test p-value: Subfield CA1 vs DG 
    ## DataFrame with 3 rows and 6 columns
    ##                baseMean    log2FoldChange             lfcSE
    ##               <numeric>         <numeric>         <numeric>
    ## C1ql2  130.705376021756 -7.65382313740708 0.786125207627745
    ## Stxbp6 143.438808384453 -5.19396604691866 0.558924222312385
    ## Crlf1  40.1931561676361 -7.69969383758498 0.843507316553684
    ##                     stat               pvalue                 padj
    ##                <numeric>            <numeric>            <numeric>
    ## C1ql2  -9.73613753018261 2.11438624194895e-22 2.63896546857649e-18
    ## Stxbp6 -9.29279111474208 1.50294385386716e-20 9.37912112005798e-17
    ## Crlf1  -9.12818856040704 6.96549786332862e-20 2.89787929440682e-16

    data <- data.frame(gene = row.names(res), 
                       pvalue = -log10(res$padj), 
                       lfc = res$log2FoldChange, 
                       padj = res$padj )
    data <- na.omit(data)
    data <- data %>%
      mutate(direction = ifelse(data$lfc > 1 & data$pvalue > 1, 
                            yes = "CA1", 
                            no = ifelse(data$lfc < -1 & data$pvalue > 1, 
                                        yes = "DG", 
                                        no = "neither")))

    data$direction <- as.factor(data$direction)
    summary(data)

    ##             gene           pvalue               lfc         
    ##  0610007P14Rik:    1   Min.   : 0.000003   Min.   :-9.3376  
    ##  0610009B22Rik:    1   1st Qu.: 0.003253   1st Qu.:-0.5547  
    ##  0610009L18Rik:    1   Median : 0.007490   Median :-0.1297  
    ##  0610009O20Rik:    1   Mean   : 0.155269   Mean   :-0.1413  
    ##  0610010F05Rik:    1   3rd Qu.: 0.052268   3rd Qu.: 0.2950  
    ##  0610010K14Rik:    1   Max.   :17.578566   Max.   : 8.4434  
    ##  (Other)      :12475                                        
    ##       padj          direction    
    ##  Min.   :0.0000   CA1    :  239  
    ##  1st Qu.:0.8866   DG     :  208  
    ##  Median :0.9829   neither:12034  
    ##  Mean   :0.8577                  
    ##  3rd Qu.:0.9925                  
    ##  Max.   :1.0000                  
    ## 

    write.csv(data, "../results/volcanoCA1DG.csv")

    # save the list of just DEGs with their pvalu and lfc and direction
    CA1DG_DEGs <- data %>%
      filter(direction != "neither")
    CA1DG_DEGs <- CA1DG_DEGs[order(CA1DG_DEGs$pvalue),]

    write.csv(CA1DG_DEGs, "../results/CA1DG_DEGs.csv")

    res <- results(dds, contrast =c("Subfield", "CA3", "DG"), independentFiltering = T, alpha = 0.1)
    summary(res)

    ## 
    ## out of 16709 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 53, 0.32%
    ## LFC < 0 (down)     : 45, 0.27%
    ## outliers [1]       : 18, 0.11%
    ## low counts [2]     : 5178, 31%
    ## (mean count < 6)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    res <- results(dds, contrast =c("Subfield", "CA1", "CA3"), independentFiltering = T, alpha = 0.1)
    summary(res)

    ## 
    ## out of 16709 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1, 0.006%
    ## LFC < 0 (down)     : 17, 0.1%
    ## outliers [1]       : 18, 0.11%
    ## low counts [2]     : 8415, 50%
    ## (mean count < 21)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

Table 1: % of DEGs
------------------

    (222+262)/16709*100

    ## [1] 2.896643

    (45+53)/16709*100

    ## [1] 0.5865103

    (17+1)/16709*100

    ## [1] 0.1077264

    (56+288)/16709*100

    ## [1] 2.058771

PCA
---

This PCA gives an overview of the variability between samples using the
a large matrix of log transformed gene expression data. You can see that
the biggest difference is between DG punches and the CA1 and CA3
punches. CA1 and CA3 samples have similar transcriptomes. The control
CA1 samples have the most similar transcriptonal profiles as evidenced
by their tight clustering.

    colorvalSubfield <- c("#7570b3", "#1b9e77", "#d95f02")
    colorvalTreatment <- c("#ffffff", "#525252")

    #rowVars(assay(rld))

    # create the dataframe using my function pcadataframe
    pcadata <- pcadataframe(rld, intgroup=c("Subfield", "Treatment"), returnData=TRUE)
    percentVar <- round(100 * attr(pcadata, "percentVar"))
    #percentVar
    write.csv(pcadata, "../results/pcadata.csv")


    PCA12 <- ggplot(pcadata, aes(PC1, PC2, shape = pcadata$Treatment)) + 
      geom_point(size = 5, alpha = 1, aes(color = Subfield)) +
      stat_ellipse(type = "t", aes(lty=pcadata$Treatment)) +
      scale_linetype_manual(values=c(3,1)) +
        xlab(paste0("PC1: ", percentVar[1],"% variance")) +
        ylab(paste0("PC2: ", percentVar[2],"% variance")) +
        scale_color_manual(values = colorvalSubfield) +
       theme_cowplot(font_size = 12, line_size = 0.25)  +
       theme_minimal() +
      xlim(-56,36) +
      ylim(-56,36) +
        scale_shape_manual(values=c(1, 16))  +
        theme(legend.position = "right",
              legend.key.width = unit(0.1,"mm"),
              legend.key.height = unit(0.1,"cm"),
              legend.title = element_blank(),
              legend.text = element_text(size = 8),
              panel.grid.minor=element_blank())

    # thanks to https://stackoverflow.com/questions/31295382/how-to-change-the-linetype-for-ellipses-in-ggplot2-with-stat-ellipse for help with elipse

    a <- ggdraw() + draw_image("../figures/00_methodsoverview/expdesign.png", scale = 1)

    figure1 <- plot_grid(a, PCA12,  
                         nrow = 1, labels = c('A', 'B'), 
                         #align = 'h',
                         rel_widths = c(2,4))

    figure1

![](../figures/01_dissociationtest/PCA-1.png)

    pdf("../figures/figure1.pdf", width=7, height=4)
    print(figure1)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    ggsave(
      "../figures/figure1.png",
      figure1,
      width = 6,
      height = 3,
      dpi = 1200
    )

PCA statistics

    aov1 <- aov(PC1 ~ Subfield * Treatment, data=pcadata)
    summary(aov1) 

    ##                    Df Sum Sq Mean Sq F value  Pr(>F)    
    ## Subfield            2 2812.7  1406.4  22.483 0.00052 ***
    ## Treatment           1  335.2   335.2   5.358 0.04932 *  
    ## Subfield:Treatment  2   38.7    19.4   0.310 0.74205    
    ## Residuals           8  500.4    62.6                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov1, which = "Subfield") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ Subfield * Treatment, data = pcadata)
    ## 
    ## $Subfield
    ##              diff       lwr      upr     p adj
    ## CA3-CA1  5.223963 -9.363897 19.81182 0.5838919
    ## DG-CA1  33.098277 18.510417 47.68614 0.0004937
    ## DG-CA3  27.874315 11.894115 43.85451 0.0027216

    aov2 <- aov(PC2 ~ Subfield * Treatment, data=pcadata)
    summary(aov2) 

    ##                    Df Sum Sq Mean Sq F value Pr(>F)  
    ## Subfield            2  243.8   121.9   0.892 0.4470  
    ## Treatment           1  691.2   691.2   5.057 0.0547 .
    ## Subfield:Treatment  2   16.9     8.5   0.062 0.9404  
    ## Residuals           8 1093.4   136.7                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov2, which = "Subfield") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ Subfield * Treatment, data = pcadata)
    ## 
    ## $Subfield
    ##              diff       lwr      upr     p adj
    ## CA3-CA1 -8.297758 -29.86072 13.26520 0.5406312
    ## DG-CA1   1.924170 -19.63879 23.48713 0.9649491
    ## DG-CA3  10.221928 -13.39911 33.84297 0.4664055

Next, save files for dowstream GO analysis.

    # from https://github.com/rachelwright8/Ahya-White-Syndromes/blob/master/deseq2_Ahya.R

    resCD=results(dds, contrast=c('Treatment', 'DISS', 'HOMO'), independentFiltering = T)
    table(resCD$padj<0.1)

    ## 
    ## FALSE  TRUE 
    ## 11813   344

    logs <- data.frame(cbind("gene"=row.names(resCD),"logP"=round(-log(resCD$pvalue+1e-10,10),1)))
    logs$logP=as.numeric(as.character(logs$logP))
    sign <- rep(1,nrow(logs))
    sign[resCD$log2FoldChange<0]=-1  ##change to correct model
    table(sign)

    ## sign
    ##   -1    1 
    ## 6989 9720

    logs$logP <- logs$logP*sign

    write.csv(logs, file = "./05_GO_MWU/GOpvals.csv", row.names = F)

To view a histogram of the p-value distibution for each constrast,
change the Rmd file to `include=TRUE` for this chunck.
