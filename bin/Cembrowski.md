Methods for Dorsal Hippocampal Gene Expression Profiling
--------------------------------------------------------

#### Part 3: Examining fac-sorted pyramindal neurons from CA1 and CA3 as well as granular cells from DG from Cembrowski et al 2016.

This ([2016 Cembrowski
paper](https://elifesciences.org/content/5/e14997#fig1s30)) is very
similar to my experiment, so I want to compare the two. Like mine, they
compare hippocampal gene expression from dorsal CA1, CA3, and DG sub
regions. These cells were identifed through fac sorting to isolate
genetically labeled CA1 and CA3 pyramical neurons and DG granular cells.

Here are three figures from their analysis.
![](../figures/cembrowski/elife-14997-fig1-v1-download.jpg)
![](../figures/cembrowski/elife-14997-fig2-v1.jpg)
![](../figures/cembrowski/elife-14997-fig6-v1-download.jpg)

The first thing I notice is that they have waay more reads per sample
and thus gene counts per sample than I do. They have a mean gene counts
per sample around 400 million counts per gene. My data had 5 million
counts per gene.

![](../figures/cembrowski/edgeR-1.png)

#### Differential Gene Expression

![](../figures/cembrowski/DifferentialGeneExpressionAnalysis-1.png)![](../figures/cembrowski/DifferentialGeneExpressionAnalysis-2.png)![](../figures/cembrowski/DifferentialGeneExpressionAnalysis-3.png)

This Venn Diagram shows the number of differentially expressed by
contrast described above each oval. The most number of genes are
differntially expressed between DG and the CAs (nearly 1000) wheras only
about 200 were differntailly regulated as a result of of technical
maniplulation comparing homogenized and dissociated samples.

![](../figures/cembrowski/VennDiagram.png)

I'm not really happy with these two heat maps. Here's how I created
them. Top heatmap: subset the data to give only the gene with an
adjusted p value \< 0.05 for the homogenized vs dissociated
comparisonany two-way comparsion. Bottom heatmap: subset the data to
give only the gene with an adjusted p value \< 0.05 for two way brain
region comparision (CA1 vs DG, CA3, vs DG, or CA1 vs DG)

Here, you can see that the differences between samples is not as clear
cut for all comparisions. What other mechanisms would be useful for
subseting the data to identify genes of interest?

![](../figures/cembrowski/Heatmap100DEgenes-1.png)![](../figures/cembrowski/Heatmap100DEgenes-2.png)

![](../figures/cembrowski/PCA-1.png)
