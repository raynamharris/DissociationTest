This R Markdown document will walk through the analysis of hippocampal tissue prepared with two different methods. The "homogenized" samples were collected by punch then homogenized in homogenization buffer from the Promega Maxwell kit. The "dissociated samples" were also collected similarily but the cells was dissociated after being punch and before being homogenized.

#### Differential Gene Expression Plots

    ## class: DESeqDataSet 
    ## dim: 17695 26 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(17695): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(26): 100-CA1-1 100-CA1-2 ... 147D-CA3-1 147D-DG-1
    ## colData names(11): RNAseqID Method ... Punch.Collector jobnumber

    ## class: DESeqDataSet 
    ## dim: 17695 26 
    ## metadata(1): version
    ## assays(5): counts mu cooks replaceCounts replaceCooks
    ## rownames(17695): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(38): baseMean baseVar ... maxCooks replace
    ## colnames(26): 100-CA1-1 100-CA1-2 ... 147D-CA3-1 147D-DG-1
    ## colData names(13): RNAseqID Method ... sizeFactor replaceable

    ## 
    ## out of 17690 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)     : 25, 0.14% 
    ## LFC < 0 (down)   : 143, 0.81% 
    ## outliers [1]     : 327, 1.8% 
    ## low counts [2]   : 0, 0% 
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    ## 
    ## out of 17690 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)     : 16, 0.09% 
    ## LFC < 0 (down)   : 120, 0.68% 
    ## outliers [1]     : 327, 1.8% 
    ## low counts [2]   : 6588, 37% 
    ## (mean count < 7)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

![](../figures/allregions_allgroups/DifferentialGeneExpressionAnalysis-1.png)

    ## NULL

![](../figures/allregions_allgroups/DifferentialGeneExpressionAnalysis-2.png)![](../figures/allregions_allgroups/DifferentialGeneExpressionAnalysis-3.png)

    ## [1] 32

    ## [1] 7

resPunchCA1DG \<- results(dds, contrast = c("Punch", "CA1", "DG"), independentFiltering = F) \#sum(resPunchCA1DG\(padj < 0.1, na.rm = TRUE) # 4170 #1127 valsPunchCA1DG <- cbind(resPunchCA1DG\)pvalue, resPunchCA1DG$padj) colnames(valsPunchCA1DG)=c("pval.CA1DG", "padj.CA1DG")

resPunchCA1CA3 \<- results(dds, contrast = c("Punch", "CA1", "CA3"), independentFiltering = F) \#sum(resPunchCA1CA3\(padj < 0.1, na.rm = TRUE) #2240 # 70 valsPunchCA1CA3 <- cbind(resPunchCA1CA3\)pvalue, resPunchCA1CA3$padj) colnames(valsPunchCA1CA3)=c("pval.CA1CA3", "padj.CA1CA3")

resPunchCA3DG \<- results(dds, contrast = c("Punch", "CA3", "DG"), independentFiltering = F) \#sum(resPunchCA3DG\(padj < 0.1, na.rm = TRUE) #4785 #591 valsPunchCA3DG <- cbind(resPunchCA3DG\)pvalue, resPunchCA3DG$padj) colnames(valsPunchCA3DG)=c("pval.CA3DG", "padj.CA3DG")

\`\`\`{r VennDiagram, echo=FALSE, message=FALSE}
================================================

rldpvals \<- as.data.frame(rldpvals)

MethodHomogDiss \<- row.names(rldpvals[rldpvals\(padj.MethodHomogDiss<0.1 & !is.na(rldpvals\)padj.MethodHomogDiss),]) \#MethodYokedTrained \<- row.names(rldpvals[rldpvals\(padj.valsMethodYokedTrained<0.1 & !is.na(rldpvals\)padj.valsMethodYokedTrained),]) PunchCA1DG \<- row.names(rldpvals[rldpvals\(padj.CA1DG<0.1 & !is.na(rldpvals\)padj.CA1DG),]) PunchCA1CA3 \<- row.names(rldpvals[rldpvals\(padj.CA1CA3<0.1 & !is.na(rldpvals\)padj.CA1CA3),]) PunchCA3DG \<- row.names(rldpvals[rldpvals\(padj.CA3DG<0.1 & !is.na(rldpvals\)padj.CA3DG),])

four way grid
-------------

candidates \<- list("CA1 v. DG" = PunchCA1DG, "CA1 v. CA3" = PunchCA1CA3, "CA3 v. DG" = PunchCA3DG, "Homogenized v. Dissociated" = MethodHomogDiss ) dev.off() prettyvenn \<- venn.diagram( x = candidates, filename=NULL, lwd=4, col = "transparent", fill = (values=c("\#00441b", "\#00441b","\#238b45", "\#238b45")), alpha = 0.5, cex = 1, fontfamily = "sans", \#fontface = "bold", cat.default.pos = "text", \#cat.col = c("darkred", "darkgreen", "blue4", "orange"), \#cat.dist = c(0.08, 0.08, 0.08, 0.08), cat.pos = 1, cat.cex = 1, cat.fontfamily = "sans") grid.draw(prettyvenn)

ca1 ca3 homo diss
-----------------

candidates \<- list("CA1 v. CA3" = PunchCA1CA3, "Homogenized v. Dissociated" = MethodHomogDiss ) dev.off() prettyvenn \<- venn.diagram( x = candidates, filename=NULL, lwd=2, col = "transparent", fill = (values=c("\#00441b", "\#00441b")), alpha = 0.5, cex = 1, fontfamily = "sans", \#fontface = "bold", cat.default.pos = "text", \#cat.col = c("darkred", "darkgreen", "blue4", "orange"), \#cat.dist = c(0.08, 0.08, 0.08, 0.08), cat.pos = 1, cat.cex = 1, cat.fontfamily = "sans") grid.draw(prettyvenn)

ca1 dg homo diss
----------------

candidates \<- list("CA1 v. DG" = PunchCA1DG, "Homogenized v. Dissociated" = MethodHomogDiss ) dev.off() prettyvenn \<- venn.diagram( x = candidates, filename=NULL, lwd=4, col = "transparent", fill = (values=c("\#00441b", "\#00441b")), alpha = 0.5, cex = 1, fontfamily = "sans", \#fontface = "bold", cat.default.pos = "text", \#cat.col = c("darkred", "darkgreen", "blue4", "orange"), \#cat.dist = c(0.08, 0.08, 0.08, 0.08), cat.pos = 1, cat.cex = 1, cat.fontfamily = "sans") grid.draw(prettyvenn)

ca3 dg homo diss
----------------

candidates \<- list("CA3 v. DG" = PunchCA3DG, "Homogenized v. Dissociated" = MethodHomogDiss ) dev.off() prettyvenn \<- venn.diagram( x = candidates, filename=NULL, lwd=2, col = "transparent", fill = (values=c("\#00441b", "\#00441b")), alpha = 0.5, cex = 1, fontfamily = "sans", \#fontface = "bold", cat.default.pos = "text", \#cat.col = c("darkred", "darkgreen", "blue4", "orange"), \#cat.dist = c(0.08, 0.08, 0.08, 0.08), cat.pos = 1, cat.cex = 1, cat.fontfamily = "sans") grid.draw(prettyvenn)

\`\`\`
======

![](../figures/allregions_allgroups/Heatmap100DEgenes-1.png) DEGes \<- as.data.frame(rldpvals) \# convert matrix to dataframe DEGes\(rownames <- rownames(DEGes) # add the rownames to the dataframe DEGes\)padjmin \<- with(DEGes, pmin(padj.MethodHomogDiss, padj.MethodYokedTrained)) \# put the min pvalue in a new column DEGes \<- DEGes %\>% filter(padj.MethodHomogDiss \< 0.05) rownames(DEGes) \<- DEGes$rownames drop.cols \<- c("padj.MethodHomogDiss", "pval.MethodHomogDiss","padj.MethodYokedTrained", "pval.MethodYokedTrained", "rownames", "padjmin") DEGes \<- DEGes %\>% select(-one\_of(drop.cols)) DEGes \<- as.matrix(DEGes) DEGes \<- DEGes - rowMeans(DEGes)

pheatmap(DEGes, show\_colnames=F, show\_rownames = TRUE, annotation\_col=df, annotation\_colors = ann\_colors, fontsize = 12, fontsize\_row = 10, \#cellwidth=10, cellheight=10, width = 10, border\_color = "grey60" , color = matlabcolors, main = "Homogenized vs Dissociated padj \< 0.05" )

DEGes \<- as.data.frame(rldpvals) \# convert matrix to dataframe DEGes\(rownames <- rownames(DEGes) # add the rownames to the dataframe DEGes\)padjmin \<- with(DEGes, pmin(padj.MethodHomogDiss, padj.MethodYokedTrained)) \# put the min pvalue in a new column DEGes \<- DEGes %\>% filter(padjmin \< 0.05) rownames(DEGes) \<- DEGes$rownames drop.cols \<- c("padj.MethodHomogDiss", "pval.MethodHomogDiss","padj.MethodYokedTrained", "pval.MethodYokedTrained", "rownames", "padjmin") DEGes \<- DEGes %\>% select(-one\_of(drop.cols)) DEGes \<- as.matrix(DEGes) DEGes \<- DEGes - rowMeans(DEGes)

pheatmap(DEGes, show\_colnames=F, show\_rownames = TRUE, annotation\_col=df, annotation\_colors = ann\_colors, fontsize = 12, fontsize\_row = 10, \#cellwidth=10, cellheight=10, width = 10, border\_color = "grey60" , color = matlabcolors, main = "Any padj \< 0.05" )

    ##                   PC1        PC2             group      Method Punch
    ## 100-CA1-1   21.844526 -12.956295 Homogenized : CA1 Homogenized   CA1
    ## 100-CA1-2   21.750539 -18.816274 Homogenized : CA1 Homogenized   CA1
    ## 100-CA1-3   19.586344 -15.902778 Homogenized : CA1 Homogenized   CA1
    ## 100-CA3-1   16.319474  -9.215784 Homogenized : CA3 Homogenized   CA3
    ## 100-CA3-4   17.557805  -8.314872 Homogenized : CA3 Homogenized   CA3
    ## 100-DG-2    -8.590330  -2.000619  Homogenized : DG Homogenized    DG
    ## 100-DG-3   -37.286963 -17.483685  Homogenized : DG Homogenized    DG
    ## 101-CA1-1   21.107411   3.057414 Dissociated : CA1 Dissociated   CA1
    ## 101-CA1-2   11.525709   8.553547 Dissociated : CA1 Dissociated   CA1
    ## 101-CA1-3    7.427149  11.091593 Dissociated : CA1 Dissociated   CA1
    ## 101-CA3-1   18.480723  -9.032546 Dissociated : CA3 Dissociated   CA3
    ## 101-CA3-4   10.678874   3.242538 Dissociated : CA3 Dissociated   CA3
    ## 101-DG-3   -11.588005  40.257655  Dissociated : DG Dissociated    DG
    ## 101-DG-4   -40.511925  -2.210746  Dissociated : DG Dissociated    DG
    ## 146C-CA1-4  16.549671  -5.425457     Trained : CA1     Trained   CA1
    ## 146C-CA3-4   5.814773  50.753192     Trained : CA3     Trained   CA3
    ## 146C-DG-4  -39.093235  -4.707806      Trained : DG     Trained    DG
    ## 146D-CA1-3  10.186217   4.951797       Yoked : CA1       Yoked   CA1
    ## 146D-CA3-3  10.891473  -1.448810       Yoked : CA3       Yoked   CA3
    ## 146D-DG-3  -28.433705  22.832486        Yoked : DG       Yoked    DG
    ## 147C-CA1-3  19.241844 -13.299751     Trained : CA1     Trained   CA1
    ## 147C-CA3-3  13.445384  -1.827744     Trained : CA3     Trained   CA3
    ## 147C-DG-3  -37.389087 -16.480214      Trained : DG     Trained    DG
    ## 147D-CA1-1  -7.257385  22.879376       Yoked : CA1       Yoked   CA1
    ## 147D-CA3-1  11.967457  -5.515119       Yoked : CA3       Yoked   CA3
    ## 147D-DG-1  -44.224736 -22.981098        Yoked : DG       Yoked    DG
    ##                  name
    ## 100-CA1-1   100-CA1-1
    ## 100-CA1-2   100-CA1-2
    ## 100-CA1-3   100-CA1-3
    ## 100-CA3-1   100-CA3-1
    ## 100-CA3-4   100-CA3-4
    ## 100-DG-2     100-DG-2
    ## 100-DG-3     100-DG-3
    ## 101-CA1-1   101-CA1-1
    ## 101-CA1-2   101-CA1-2
    ## 101-CA1-3   101-CA1-3
    ## 101-CA3-1   101-CA3-1
    ## 101-CA3-4   101-CA3-4
    ## 101-DG-3     101-DG-3
    ## 101-DG-4     101-DG-4
    ## 146C-CA1-4 146C-CA1-4
    ## 146C-CA3-4 146C-CA3-4
    ## 146C-DG-4   146C-DG-4
    ## 146D-CA1-3 146D-CA1-3
    ## 146D-CA3-3 146D-CA3-3
    ## 146D-DG-3   146D-DG-3
    ## 147C-CA1-3 147C-CA1-3
    ## 147C-CA3-3 147C-CA3-3
    ## 147C-DG-3   147C-DG-3
    ## 147D-CA1-1 147D-CA1-1
    ## 147D-CA3-1 147D-CA3-1
    ## 147D-DG-1   147D-DG-1

![](../figures/allregions_allgroups/PCA-1.png)

``` r
library(edgeR)
```

    ## Warning: package 'edgeR' was built under R version 3.3.2

    ## Loading required package: limma

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:DESeq2':
    ## 
    ##     plotMA

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

``` r
counts <- countData
dim( counts )
```

    ## [1] 22485    26

``` r
colSums( counts )
```

    ##  100-CA1-1  100-CA1-2  100-CA1-3  100-CA3-1  100-CA3-4   100-DG-2 
    ##    2311086    6646655    2277596    1974845    2352153    1285654 
    ##   100-DG-3  101-CA1-1  101-CA1-2  101-CA1-3  101-CA3-1  101-CA3-4 
    ##    6086605    4782767     135065     300812    2498914    1193153 
    ##   101-DG-3   101-DG-4 146C-CA1-4 146C-CA3-4  146C-DG-4 146D-CA1-3 
    ##      65887     598775    1360004     257822     492145     391369 
    ## 146D-CA3-3  146D-DG-3 147C-CA1-3 147C-CA3-3  147C-DG-3 147D-CA1-1 
    ##    2994536      90417    3072308    5754581    4350647        213 
    ## 147D-CA3-1  147D-DG-1 
    ##    4624995   11700703

``` r
colSums( counts ) / 1e06  # in millions of reads
```

    ##  100-CA1-1  100-CA1-2  100-CA1-3  100-CA3-1  100-CA3-4   100-DG-2 
    ##   2.311086   6.646655   2.277596   1.974845   2.352153   1.285654 
    ##   100-DG-3  101-CA1-1  101-CA1-2  101-CA1-3  101-CA3-1  101-CA3-4 
    ##   6.086605   4.782767   0.135065   0.300812   2.498914   1.193153 
    ##   101-DG-3   101-DG-4 146C-CA1-4 146C-CA3-4  146C-DG-4 146D-CA1-3 
    ##   0.065887   0.598775   1.360004   0.257822   0.492145   0.391369 
    ## 146D-CA3-3  146D-DG-3 147C-CA1-3 147C-CA3-3  147C-DG-3 147D-CA1-1 
    ##   2.994536   0.090417   3.072308   5.754581   4.350647   0.000213 
    ## 147D-CA3-1  147D-DG-1 
    ##   4.624995  11.700703

``` r
table( rowSums( counts ) )[ 1:30 ] # Number of genes with low counts
```

    ## 
    ##    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
    ## 4461  329  283  233  179  152  150  117  117   97   99   80   75   77   67 
    ##   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29 
    ##   76   65   62   74   63   49   64   40   52   50   51   53   52   31   51
