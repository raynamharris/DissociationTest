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

Rayna's Heat maps.

The top heatmap shows clean separation of each of the 6 groups. Dorsal
CA1 are most different from the rest. Ventral CA1 and CA3 are similar to
one another and to ventral CA3. DGs cluster well.

The bottom heat map is a much less stringent cutoff and this one cleanly
separates first by brain region and then by dorsal ventral location.

![](../figures/cembrowski/Heatmap100DEgenes-1.png)![](../figures/cembrowski/Heatmap100DEgenes-2.png)

![](../figures/cembrowski/PCA-1.png)
