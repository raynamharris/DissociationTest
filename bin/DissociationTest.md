This R Markdown document will walk through the analysis of hippocampal tissue prepared with two different methods. The "homogenized" samples were collected by punch then homogenized in homogenization buffer from the Promega Maxwell kit. The "dissociated samples" were also collected similarily but the cells was dissociated after being punch and before being homogenized.

Save Intermediate Data file types, so they could be loaded from here with StoreReadData write.table(countData, '../intermediatefiles/countData.csv', row.names = TRUE, sep=",", col.names = T) write.table(colData, '../intermediatefiles/colData.csv', row.names = TRUE, sep=",", col.names = T) write.table(geneids, '../intermediatefiles/geneids.csv', row.names = TRUE, sep=",", col.names = T) read.csv('../intermediatefiles/colData.csv') read.csv('../intermediatefiles/countData.csv')

write.table(countData, '../intermediatefiles/countData.csv', row.names = TRUE, sep=",", col.names = T) write.table(colData, '../intermediatefiles/colData.csv', row.names = TRUE, sep=",", col.names = T) write.table(geneids, '../intermediatefiles/geneids.csv', row.names = TRUE, sep=",", col.names = T) colData \<- read.csv('../intermediatefiles/colData.csv') countData \<- read.csv('../intermediatefiles/countData.csv')

#### Differential Gene Expression Plots

![](../figures/panel1DifferentialGeneExpressionAnalysis-1.png)

    ## NULL

![](../figures/panel1DifferentialGeneExpressionAnalysis-2.png)![](../figures/panel1DifferentialGeneExpressionAnalysis-3.png)

![](../figures/panel1VennDiagram-1.png)

![](../figures/panel1Heatmap100DEgenes-1.png)

    ##                  PC1         PC2    group Punch.Collector Punch      name
    ## 100-CA1-1 -17.116137   8.0509162 MK : CA1              MK   CA1 100-CA1-1
    ## 100-CA1-2 -18.755714  10.6164908 MK : CA1              MK   CA1 100-CA1-2
    ## 100-CA1-3 -17.143935  11.1990267 MK : CA1              MK   CA1 100-CA1-3
    ## 100-CA3-1 -11.752472   1.1049952 MK : CA3              MK   CA3 100-CA3-1
    ## 100-CA3-4  -8.073388  -2.9659732 MK : CA3              MK   CA3 100-CA3-4
    ## 100-DG-2   11.180665   0.2617668  MK : DG              MK    DG  100-DG-2
    ## 100-DG-3   27.330910  21.0396030  MK : DG              MK    DG  100-DG-3
    ## 101-CA1-1 -14.236503   3.2659422 MK : CA1              MK   CA1 101-CA1-1
    ## 101-CA1-2  -3.873025  -7.0944387 MK : CA1              MK   CA1 101-CA1-2
    ## 101-CA1-3   5.481191 -15.2960486 MK : CA1              MK   CA1 101-CA1-3
    ## 101-CA3-1  -7.898660  -6.7391440 MK : CA3              MK   CA3 101-CA3-1
    ## 101-CA3-4   4.902719 -17.2677898 MK : CA3              MK   CA3 101-CA3-4
    ## 101-DG-3   18.535221 -20.8252619  MK : DG              MK    DG  101-DG-3
    ## 101-DG-4   31.419128  14.6499151  MK : DG              MK    DG  101-DG-4

![](../figures/panel1PCA-1.png)

    ##                  PC1         PC2             group      Method Punch
    ## 100-CA1-1 -17.116137   8.0509162 Homogenized : CA1 Homogenized   CA1
    ## 100-CA1-2 -18.755714  10.6164908 Homogenized : CA1 Homogenized   CA1
    ## 100-CA1-3 -17.143935  11.1990267 Homogenized : CA1 Homogenized   CA1
    ## 100-CA3-1 -11.752472   1.1049952 Homogenized : CA3 Homogenized   CA3
    ## 100-CA3-4  -8.073388  -2.9659732 Homogenized : CA3 Homogenized   CA3
    ## 100-DG-2   11.180665   0.2617668  Homogenized : DG Homogenized    DG
    ## 100-DG-3   27.330910  21.0396030  Homogenized : DG Homogenized    DG
    ## 101-CA1-1 -14.236503   3.2659422 Dissociated : CA1 Dissociated   CA1
    ## 101-CA1-2  -3.873025  -7.0944387 Dissociated : CA1 Dissociated   CA1
    ## 101-CA1-3   5.481191 -15.2960486 Dissociated : CA1 Dissociated   CA1
    ## 101-CA3-1  -7.898660  -6.7391440 Dissociated : CA3 Dissociated   CA3
    ## 101-CA3-4   4.902719 -17.2677898 Dissociated : CA3 Dissociated   CA3
    ## 101-DG-3   18.535221 -20.8252619  Dissociated : DG Dissociated    DG
    ## 101-DG-4   31.419128  14.6499151  Dissociated : DG Dissociated    DG
    ##                name
    ## 100-CA1-1 100-CA1-1
    ## 100-CA1-2 100-CA1-2
    ## 100-CA1-3 100-CA1-3
    ## 100-CA3-1 100-CA3-1
    ## 100-CA3-4 100-CA3-4
    ## 100-DG-2   100-DG-2
    ## 100-DG-3   100-DG-3
    ## 101-CA1-1 101-CA1-1
    ## 101-CA1-2 101-CA1-2
    ## 101-CA1-3 101-CA1-3
    ## 101-CA3-1 101-CA3-1
    ## 101-CA3-4 101-CA3-4
    ## 101-DG-3   101-DG-3
    ## 101-DG-4   101-DG-4

![](../figures/panel1PCA-2.png)
