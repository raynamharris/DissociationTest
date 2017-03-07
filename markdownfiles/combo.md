Methods for Dorsal Hippocampal Gene Expression Profiling
--------------------------------------------------------

#### Part 4: Comparing Rayna's and the Cembrowski data.

This histogram shows that my samples are in the bottom bin with less
than 100 million gene counts per sample while the cembrowksi datat has
up to 900 milion reads per sample.

![](../figures/combo/edgeR-1.png)

However, this first plot shows that normalization with DESeq2 can put my
samples in a comparble range as the cembrowski samples on a gene by gene
basis.

![](../figures/combo/DifferentialGeneExpressionAnalysis-1.png)![](../figures/combo/DifferentialGeneExpressionAnalysis-2.png)

Here are the number of differentially expressed genes. ~<sub>~</sub>
sum(rescembrowskiharris$padj \< 0.1, na.rm = TRUE) \#7487 sum(resCA1DG$padj
\< 0.1, na.rm = TRUE) \#4492
sum(resCA3DG$padj \< 0.1, na.rm = TRUE) \#4659 sum(resCA1CA3$padj \<
0.1, na.rm = TRUE) \#1732 ~<sub>~</sub>

![](../figures/combo/VennDiagram-1.png)

This plot shows the clear separation between the cembrowski and the
harris data. Next, the regions and the the dorsal ventral patterning
seaparates in the cembrowski dataset. Next, DG in my data appeas
distinct but CA1 and CA3 are mixed. In the bottom heatmap, you can see
that including more genes helps separate CA1 and CA3 in my data.

![](../figures/combo/Heatmap100DEgenes-1.png)![](../figures/combo/Heatmap100DEgenes-2.png)

I really love this PCA analysis. PC1 separates the two experiments. PC2
separates DG from the CAs.

![](../figures/combo/PCA-1.png)
