All together now
----------------

Combining the two previous analyses

    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ Treatment + Region + Treatment * Region )
    dds <- dds[ rowSums(counts(dds)) > 2, ] ## filter genes with 0 counts
    dds <- DESeq(dds) # Differential expression analysis

    FALSE estimating size factors

    FALSE estimating dispersions

    FALSE gene-wise dispersion estimates

    FALSE mean-dispersion relationship

    FALSE final dispersion estimates

    FALSE fitting model and testing

    dds

    FALSE class: DESeqDataSet 
    FALSE dim: 16970 22 
    FALSE metadata(1): version
    FALSE assays(3): counts mu cooks
    FALSE rownames(16970): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    FALSE rowData names(37): baseMean baseVar ... deviance maxCooks
    FALSE colnames(22): 142C_CA1 142C_DG ... 147D-CA3-1 147D-DG-1
    FALSE colData names(15): RNAseqID Mouse ... Date sizeFactor

    ## log transformed data
    rld <- rlog(dds, blind=FALSE)

PCA

    source("DESeqPCAfunction.R")
    source("figureoptions.R")

    # create the dataframe using my function pcadataframe
    pcadata <- pcadataframe(rld, intgroup=c("Treatment", "Region"), returnData=TRUE)
    percentVar <- round(100 * attr(pcadata, "percentVar"))

    pcadata$Treatment <- factor(pcadata$Treatment, levels = c("unavoidable", "avoidable"))

    ## PC1 vs PC2
    plotPC1PC2(aescolor = pcadata$Region, colorname = "Region", aesshape = pcadata$Treatment, shapename = "Treatment", colorvalues = colorvalRegion)

![](../figures/03_cognitiontest/PCA-1.png)

    # PC1 vs PC2 for adobe
    myplot <- plotPC1PC2(aescolor = pcadata$Region, colorname = "Region", aesshape = pcadata$Treatment, shapename = "Treatment", colorvalues = colorvalRegion)
    pdf(file="../figures/03_cognitiontest/PCA-1.pdf", width=4.5, height=3)
    plot(myplot)
    dev.off()

    FALSE quartz_off_screen 
    FALSE                 2

    ## PC2 vs PC3
    plotPC2PC3(aescolor = pcadata$Region, colorname = "Region", aesshape = pcadata$Treatment, shapename = "Treatment", colorvalues = colorvalRegion)

![](../figures/03_cognitiontest/PCA-2.png)

    plotPC3PC4(aescolor = pcadata$Region, colorname = "Region", aesshape = pcadata$Treatment, shapename = "Treatment", colorvalues = colorvalRegion)

![](../figures/03_cognitiontest/PCA-3.png)

    ## DEG by contrasts
    source("resvalsfunction.R")
    contrast1 <- resvals(contrastvector = c('Region', 'CA1', 'DG'), mypval = 0.1)

    ## [1] 3919
    ## [1] 1546

    contrast2 <- resvals(contrastvector = c('Region', 'CA3', 'DG'), mypval = 0.1)

    ## [1] 4690
    ## [1] 2664

    contrast3 <- resvals(contrastvector = c('Region', 'CA1', 'CA3'), mypval = 0.1)

    ## [1] 3745
    ## [1] 1194

    contrast4 <- resvals(contrastvector = c('Treatment', 'avoidable', 'unavoidable'), mypval = 0.1)

    ## [1] 3041
    ## [1] 285

Now, we can view a histogram of the distribution

![](../figures/03_cognitiontest/histogram-1.png)

    ## [1] 1

![](../figures/03_cognitiontest/histogram-2.png)

    ## [1] 1

![](../figures/03_cognitiontest/histogram-3.png)

    ## [1] 1

![](../figures/03_cognitiontest/histogram-4.png)

    ## [1] 1

This Venn Diagram sthe overlap of differentailly expression genes by
Region and method. This shows all genes with *adjusted* pvalue &lt;0.1.

![](../figures/03_cognitiontest/VennDiagramPadj-1.png)

Heatmaps

    ## Any padj <0.1
    DEGes <- assay(rld)
    DEGes <- cbind(DEGes, contrast1, contrast2, contrast3, contrast4)
    DEGes <- as.data.frame(DEGes) # convert matrix to dataframe
    DEGes$rownames <- rownames(DEGes)  # add the rownames to the dataframe

    DEGes$padjmin <- with(DEGes, pmin(padjTreatmentavoidableunavoidable, padjRegionCA1DG ,padjRegionCA3DG, padjRegionCA1CA3 )) # put the min pvalue in a new column
    DEGes <- DEGes %>% filter(padjmin < 0.1)

    rownames(DEGes) <- DEGes$rownames
    drop.cols <-colnames(DEGes[,grep("padj|pval|rownames", colnames(DEGes))])
    DEGes <- DEGes %>% select(-one_of(drop.cols))
    DEGes <- as.matrix(DEGes)
    DEGes <- DEGes - rowMeans(DEGes)


    # setting color options
    source("figureoptions.R")
    ann_colors <- ann_colorsbehavior
    colorpalette <- cembrowskicolors
    df <- as.data.frame(colData(dds)[,c("Treatment", "Region")])
    paletteLength <- 30
    myBreaks <- c(seq(min(DEGes), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(DEGes)/paletteLength, max(DEGes), length.out=floor(paletteLength/2)))


    pheatmap(DEGes, show_colnames=F, show_rownames = F,
             annotation_col=df, annotation_colors = ann_colors,
             fontsize = 12, fontsize_row = 7, 
             cellwidth=10, 
             border_color = "grey60" ,
             color = colorpalette,
             clustering_distance_cols="correlation" ,
             breaks=myBreaks
             )

![](../figures/03_cognitiontest/HeatmapPadj-1.png)

    # for adobe
    pheatmap(DEGes, show_colnames=F, show_rownames = F,
             annotation_col=df, annotation_colors = ann_colors,
             treeheight_row = 0, treeheight_col = 25,
             fontsize = 11, 
             width=4.5, height=3,
             border_color = "grey60" ,
             color = colorpalette,
             cellwidth = 9, 
             filename = "../figures/03_cognitiontest/HeatmapPadj-1.pdf",
             clustering_distance_cols="correlation" ,
             breaks=myBreaks
             )

![](../figures/03_cognitiontest/HeatmapPvalue-1.png)

    FALSE 
    FALSE FALSE  TRUE 
    FALSE  8889   504

    FALSE 
    FALSE FALSE  TRUE 
    FALSE 14951  2000

    FALSE log2 fold change (MLE): Treatment avoidable vs unavoidable 
    FALSE Wald test p-value: Treatment avoidable vs unavoidable 
    FALSE DataFrame with 6 rows and 6 columns
    FALSE                baseMean log2FoldChange     lfcSE       stat    pvalue
    FALSE               <numeric>      <numeric> <numeric>  <numeric> <numeric>
    FALSE 0610007P14Rik 28.312171     -0.4917430 0.5416911 -0.9077924 0.3639879
    FALSE 0610009B22Rik  7.764007      1.1336111 1.1791370  0.9613905 0.3363559
    FALSE 0610009L18Rik  3.652118      2.3345187 1.9670998  1.1867820 0.2353136
    FALSE 0610009O20Rik 58.892418      0.6145425 0.3810066  1.6129444 0.1067566
    FALSE 0610010F05Rik 10.601326      0.1150379 0.7181401  0.1601886 0.8727325
    FALSE 0610010K14Rik  1.780900     -1.3461274 1.1956145 -1.1258875 0.2602132
    FALSE                    padj
    FALSE               <numeric>
    FALSE 0610007P14Rik 0.6743469
    FALSE 0610009B22Rik        NA
    FALSE 0610009L18Rik        NA
    FALSE 0610009O20Rik 0.3933625
    FALSE 0610010F05Rik 0.9513260
    FALSE 0610010K14Rik        NA

    FALSE sign
    FALSE   -1    1 
    FALSE 7986 8984

![](../figures/03_cognitiontest/numshocks-1.png)

    FALSE quartz_off_screen 
    FALSE                 2

![](../figures/03_cognitiontest/ephys-1.png)![](../figures/03_cognitiontest/ephys-2.png)

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

![](../figures/03_cognitiontest/pvclust-1.png)
