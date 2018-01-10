# Analysis of hippocampal transcriptomic responses to technical and biological perturbations
Rayna M. Harris, Hsin-Yi Kao, Juan Marcos Alarcon, Hans A. Hofmann, and André A. Fenton

## Overview
Cost-effective next-generation sequencing has made unbiased gene expression investigations possible. Gene expression studies at the level of single neurons may be especially important for understanding nervous system structure and function because of neuron-specific functionality and plasticity. While cellular dissociation is a prerequisite technical manipulation for such single-cell studies, the extent to which the process of dissociating cells affects neural gene expression has not been determined. Here, we examine the effect of cellular dissociation on gene expression in the mouse hippocampus. We also determine to which extent such changes might confound studies on the behavioral and physiological functions of hippocampus. 

## Video Abstract

![screenshot](./figures/screenshot.png)
https://www.youtube.com/watch?v=taeAqimxXWo

## Approach

We first processed tissue punch samples from the dentate gyrus, CA3, and CA1 hippocampus subfields using either a tissue homogenization protocol or a cellular dissociation protocol. Then we used RNA sequencing to quantify sub-field specific gene expression differences between the preparation methods, while holding the tissue content constant within a region. To evaluate the impact of the tissue preparation method on the results of biological manipulations, we compared the gene expression differences to differences in homogenate tissue samples after mouse subjects were either exposed to a stressful experience or cognitive training, which are both common in behavioral research. The bioinformatic workflow can be broken down into multiple steps: High performance computing (HPC) for RNA-sequencing (RNA-seq) analysis; reproducible research workflows in R; statistical anlysis, data visualization, and data sharing. We assessed the extent to which the subfield-specific gene expression patterns are consistent with those identified in a recently published hippocampus-subfield specific gene expression database. 

<img src="https://github.com/raynamharris/DissociationTest/blob/master/figures/00_methodsoverview/00_dissociationmethods-02.png" width="400px" align="middle"/>

## Repo Contents
- [**data**](./data/): contains most of the input processed data files. Large data fiels are stored in the Gene Expression Omnibus at [GSE99765](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99765) and [GSE100225](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100225). Raw kallisto abundance files are also stored in my other GitHub repo called [MouseHippocampusRNAseqData](https://github.com/raynamharris/MouseHippocampusRNAseqData).
- [**UNIXworkflow**](./UNIXworkflow/): This descirbes the process I used to process my files using the Stampede Cluster at the Texas Advanced computing facility
- [**scripts**](./scripts/): this contains all the .R and .Rmd scripts as well as the .md output files. They have prefixes to hint at the order of operation. The workflow is described in more detail below. Note: The files are numbered consequtively starting from 00 to indidate the order of operations. To reproduce the 00_KallistoGeneCounts.Rmd file you will first have to download some pretty large files from the internet. However, if you start from 01_DissociationTest, all the data and compute power you need is there to explore at your leisure.
- [**figures**](./figures/): Contains all output for all files from the Rmarkdown scripts and my adobe-created images. 

## Results

### Figure 1: The effects of cellular dissociation on hippocampal transcriptomes

[Here are the analyses](./scripts/01_DissociationTest.Rmd) and [the results.](./scripts/01_DissociationTest.md) 

<img src="https://github.com/raynamharris/DissociationTest/blob/master/figures/fig_01-dissociation.png" width="600px" align="middle"/>

**Figure 1. The effect of cellular dissociation on hippocampal transcriptomes.** A) From a single female mouse, we collected 2 CA1, CA3, and DG hippocampal tissue samples. One sample was subjected to a cellular dissociation treatment (dissociated) whereas the control samples (control) were standardly homogenized. B) We identified 162 dissociation-induced changes in gene expression, 331 genes with region-specific expression patterns, and 30 genes differentially expressed by both region and treatment (FDR p-value < 0.05). C) Hierarchical clustering separates the hippocampal sub-fields of the homogenized samples (light gray) but not the dissociated samples (dark gray). D) PC1 accounts for 40% of all gene expression variation and by inspection, separates the DG samples from the CA1 and CA3 samples. PC2 accounts for 22% of the variation in gene expression and varies significantly with treatment. Ellipses are hand-drawn.


### Figure 2: The effects of stressful experience on hippocampal transcriptomes
[Here are the analyses](./scripts/02_StressTest.Rmd) and [the results.](./scripts/02_StressTest.md)

<img src="https://github.com/raynamharris/DissociationTest/blob/master/figures/fig_02-stress.png" width="600px" align="middle"/>

**Figure 2. The effects of a stressful experience on hippocampal transcriptomes.** A) We compared CA1, CA3, and DG tissue samples from control mice taken directly from their homecage to mice that were subjected to a mild footshock. B) We identified 0 genes that responded to treatment, and 1669 genes that were differentially regulated across regions of the hippocampus (FDR p-value < 0.05). C) Hierarchical clusters groups samples by brain region but distinct treatments clusters are not present. D) PC1 accounts for 31% of the variation and visually separates the DG samples from the CA1 and CA3 samples. PC2 accounts for 18% of the variation and distinguish the three subfields. Ellipses were hand-drawn.


### Figure 3: The effects of cognitive training on hippocampal transcriptomes.
[Here are the analyses](./scripts/03_CognitionTest.Rmd) and [the results.](./scripts/03_CognitionTest.md)

<img src="https://github.com/raynamharris/DissociationTest/blob/master/figures/fig_03-memory.png" width="600px" align="middle"/>

**Figure 3. Effects of a learned avoidance behavior on hippocampal transcriptomes.** A) Mice used in this study were either subjected to random but mild foot shocks (control) or subjected to mild foot shocks conditioned with spatial cues (trained). Tissue samples were collected from CA1, CA3, and DG. B) We identified only 423 genes that were significantly expressed according to cognitive training and identified 3485 genes that were differentially expressed between any of the three brain regions (FDR p-value <0.05). C) Hierarchical clustering of the differentially expressed genes gives rise to three distinct clusters corresponding to the three brain regions, with CA1 and CA3 being more similar to one another than to DG. D) A principal component analysis of all genes in analysis (regardless of level of significance) shows that PC1 accounts for 56% of the variation and distinguishes the DG samples and the CA1 and CA3. PC2 accounts for 19% of the variation and distinguishes all three subfields. Ellipses were hand-drawn.


### Figure 4 and 5: Meta analyses
Here are [the analyses](./scripts/05_metaanlyses.Rmd) conducted to calculate the overlapp in significantly expressed genes and [the results](./scripts/05_metaanlyses.md).
 
Here are [the analyses of Gene Ontology](./scripts/06_GO_MWU/06_GO_MWU.Rmd) and [the results.](./markdownfiles/06_GO_MWU/06_GO_MWU.md)

<img src="https://github.com/raynamharris/DissociationTest/blob/master/figures/05_metaanalyses/05_meta123-01.png" width="600px" align="middle"/>

**Figure 4. Unique and shared responses to technical treatments and biological perturbations.** A) The number of genes that responded to chemical dissociation (163 genes), a stressful experience (0 genes), and cognitive training (423 genes). The three genes that respond to both technical and biological perturbation are Epha6, Grin2a, and Itbp3. B, C) The molecular function of gene ontology (GO) categories that are significantly enriches with either up- or down-regulated genes in response to cellular dissociation (B) or cognitive training (C). The top 10 most significant GO terms are visualized, each with a p-value < 0.001. The fraction next to GO category name indicates the fraction of genes in that category that survived a 10% FDR threshold for significance. Zero terms survived a 10% FDR threshold in response to a stressful experience.


<img src="https://github.com/raynamharris/DissociationTest/blob/master/figures/05_metaanalyses/05_meta1234-01.png" width="600px" align="middle"/>

**Figure 5. Meta analysis to incorporate public data.** A) This Venn diagram shows the overlap in brain-region specific gene expression across all four experiments (cellular dissociation, stressor habitation, cognitive training, and a public dataset examining subfield comparisons). Using this approach, we identified 146 genes that were differentially expressed between any two subfields of the hippocampus in all four experiments. B) Those 146 provide robust brain-region specific markers of gene expression belong to molecular function and cellular compartment GO. The top 10 most significant GO terms are visualized, each with a p-value < 0.05. The fraction next to GO category name indicates the fraction of genes in that category that survived a 10% FDR threshold for significance. 


### Supplementary analyses: Reproducing the Cembrowski results.  
I reanalyze the [Cembrowski et al 2016](https://elifesciences.org/content/5/e14997#fig1s3) data, which has been used to create a database of sub-region specific hippocampalgene expression. There data and mine share some commonalities, so I wanted to compare the two. Like mine, they compare hippocampal gene expression from dorsal CA1, CA3, and DG sub regions. These cells were identifed through fac sorting to isolate genetically labeled CA1 and CA3 pyramical neurons and DG granular cells. [Here are the analyses.](./scripts/04_Cembrowski.Rmd) [Here are the results](./scripts/04_Cembrowski.md).

(Note: These figures are archived on [FigShare](https://figshare.com/articles/Integrative_analysis_of_hippocampal_genomic_plasticity/5116192).)

## Conclusions
Many factors contribute to variation in gene expression. We set out to identify the extent to which the process of cellular dissociation – which allows for single cell analysis of neurons – has an appreciable effect on our ability to detect biologically meaningful variation in hippocampal gene expression. We conclude that there are specific dissociation-induced and cognitive training-induced changes in gene expression that are largely non-overlapping. It is encouraging that the overlap between cellular dissociation and cognitive training is small, indicating that these technical and biological processes affect different transcriptional processes. It is also encouraging to know that the stressful experience had no substantial effect on hippocampal gene expression, which if generalizable to other tasks will allow for using behavioral control groups and behavioral manipulations that also induce modest, potentially confounding stress. These findings provide insight into how cellular and biological manipulations influence gene expression. Through meta analysis and comparison to the published literature, we have identified a subset of robust sub-region specific markers of gene expression profiling. Further research is clearly needed to uncover the influence of other variables on variation in hippocampal gene expression. 
