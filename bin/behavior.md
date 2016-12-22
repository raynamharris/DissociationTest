Methods for Dorsal Hippocampal Gene Expression Profiling
--------------------------------------------------------

#### Part 2: Comparing the technial manipulation described in Part 1 with samples collected after a behavioral manipulation

For this section of the results, I added 2 or 3 samples from CA1, CA3, and DG collected from male mice that had the exprience of being in the active place avoidance arena, including spatially-trained and yoked animals. 

In this principle component anlaysis, PC1 separates the DG samples from the CA1 samples (which seems to happen in every analysis ever). The DG samples don't overlap with an other groups, but the CAs are not completeley separated by PC2. The CA1 homogenized and and trained samples have the least amount of variance of all CA1 samples. 

![](../figures/allregions_allgroups/PCA-2.png)

If you select genes that are differentially expressed in one contrast or another, you don't see clear patterns of up and down regulation. In this heatmap, I took all genes that were significantly different based on any two-way comparison. 

![](../figures/allregions_allgroups/Heatmap100DEgenes-3.png)
![](../figures/allregions_allgroups/Heatmap100DEgenes-1.png)
![](../figures/allregions_allgroups/Heatmap100DEgenes-2.png)


This is a data validation check plot. Here, I'm showing how many millions of reads were present in each sample. On average, each sample had xxmillion reads, but the range was from x to x millino reads.

![](../figures/allregions_allgroups/readcounts-1.png)

This graph examines the magnitude of gene expression differences (shown as log fold change on the y axis) as a function of read abundance (shown as mean normalized counts on the x axis.
![](../figures/allregions_allgroups/MAplot-1.png)

This is the gene with the most significant p value 

![](../figures/allregions_allgroups/Heatmap100DEgenes-1.png)
![](../figures/allregions_allgroups/Heatmap100DEgenes-2.png)
![](../figures/allregions_allgroups/Heatmap100DEgenes-3.png)

