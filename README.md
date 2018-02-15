# Analysis of hippocampal transcriptomic responses to technical and biological perturbations
Rayna M. Harris, Hsin-Yi Kao, Juan Marcos Alarcon, Hans A. Hofmann, and André A. Fenton

## Overview
Cost-effective next-generation sequencing has made unbiased gene expression investigations possible. Gene expression studies at the level of single neurons may be especially important for understanding nervous system structure and function because of neuron-specific functionality and plasticity. While cellular dissociation is a prerequisite technical manipulation for such single-cell studies, the extent to which the process of dissociating cells affects neural gene expression has not been determined. Here, we examine the effect of cellular dissociation on gene expression in the mouse hippocampus. We also determine to which extent such changes might confound studies on the behavioral and physiological functions of hippocampus. 

## Video Abstract

I made this [short video](https://www.youtube.com/watch?v=taeAqimxXWo) explaining how to use this GitHub repo when I submitted the first draft to the journal Hippocampus and posted a pre-print on BioRxiv. 

[![screenshot](./figures/screenshot.png)](https://www.youtube.com/watch?v=taeAqimxXWo)

## Repo Contents
- [**data**](./data/): contains most of the input processed data files. Large data fiels are stored in the Gene Expression Omnibus at [GSE99765](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99765) and [GSE100225](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100225). Raw kallisto abundance files are also stored in my other GitHub repo called [MouseHippocampusRNAseqData](https://github.com/raynamharris/MouseHippocampusRNAseqData).
- [**UNIXworkflow**](./UNIXworkflow/): This descirbes the process I used to process my files using the Stampede Cluster at the Texas Advanced computing facility
- [**scripts**](./scripts/): this contains all the .R and .Rmd scripts as well as the .md output files. 
  - They have prefixes to hint at the order of operation.
  - The order was dramatically differnt when this work was first submitted for publication (version 1).
  - The current order is broken down by each figure.
- [**figures**](./figures/): Contains all output for all files from the Rmarkdown scripts and my adobe-created images. 

