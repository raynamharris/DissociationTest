# Analysis of hippocampal transcriptomic responses to technical and biological perturbations
Rayna M. Harris, Hsin-Yi Kao, Juan Marcos Alarcon, Hans A. Hofmann, and André A. Fenton

## Overview
Cost-effective next-generation sequencing has made unbiased gene expression investigations possible. Gene expression studies at the level of single neurons may be especially important for understanding nervous system structure and function because of neuron-specific functionality and plasticity. While cellular dissociation is a prerequisite technical manipulation for such single-cell studies, the extent to which the process of dissociating cells affects neural gene expression has not been determined. Here, we examine the effect of cellular dissociation on gene expression in the mouse hippocampus. We also determine to which extent such changes might confound studies on the behavioral and physiological functions of hippocampus. 

## Video Abstract

I made this short video explaining how to use this GitHub repo when I submitted the first draft to the journal Hippocampus and posted a pre-print on BioRxiv.

[![screenshot](./figures/screenshot.png)](https://www.youtube.com/watch?v=taeAqimxXWo)

https://www.youtube.com/watch?v=taeAqimxXWo

## Approach

We first processed tissue punch samples from the dentate gyrus, CA3, and CA1 hippocampus subfields using either a tissue homogenization protocol or a cellular dissociation protocol. Then we used RNA sequencing to quantify sub-field specific gene expression differences between the preparation methods, while holding the tissue content constant within a region. To evaluate the impact of the tissue preparation method on the results of biological manipulations, we compared the gene expression differences to differences in homogenate tissue samples after mouse subjects were either exposed to a stressful experience or cognitive training, which are both common in behavioral research. The bioinformatic workflow can be broken down into multiple steps: High performance computing (HPC) for RNA-sequencing (RNA-seq) analysis; reproducible research workflows in R; statistical anlysis, data visualization, and data sharing. We assessed the extent to which the subfield-specific gene expression patterns are consistent with those identified in a recently published hippocampus-subfield specific gene expression database. 

<img src="https://github.com/raynamharris/DissociationTest/blob/master/figures/00_methodsoverview/00_dissociationmethods-02.png" width="400px" align="middle"/>

## Repo Contents
- [**data**](./data/): contains most of the input processed data files. Large data fiels are stored in the Gene Expression Omnibus at [GSE99765](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99765) and [GSE100225](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100225). Raw kallisto abundance files are also stored in my other GitHub repo called [MouseHippocampusRNAseqData](https://github.com/raynamharris/MouseHippocampusRNAseqData).
- [**UNIXworkflow**](./UNIXworkflow/): This descirbes the process I used to process my files using the Stampede Cluster at the Texas Advanced computing facility
- [**scripts**](./scripts/): this contains all the .R and .Rmd scripts as well as the .md output files. They have prefixes to hint at the order of operation. The workflow is described in more detail below. Note: The files are numbered consequtively starting from 00 to indidate the order of operations. To reproduce the 00_KallistoGeneCounts.Rmd file you will first have to download some pretty large files from the internet. However, if you start from 01_DissociationTest, all the data and compute power you need is there to explore at your leisure.
- [**figures**](./figures/): Contains all output for all files from the Rmarkdown scripts and my adobe-created images. 

