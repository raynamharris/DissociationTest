# DissociationTest

In this repo, 

## Dissociation vs Homogenized
This data was collected by me in collaboration with Andre Fenton and Maddy Kao. I compare hippocampal gene expression from CA1, CA3, and DG tissue samples. These samples were collected with a circular punch centered on the pyramidal and granular cells to enrich for that cell type. Biological replicates were tried in two ways before RNA extraction: one was dissociated before cell lysis while the other was homogenized before cell lysis. We wanted to identify the effects of treatment on gene expression. 

For a summary of the analysis, click [here](bin/DissociationTest.md)

## Dissociation vs Homogenized vs Yoked versus Trained 

For this part, I pulled in a few samples from another dataset to compare the difference between biological samples with a behavioral manipulation. 

For a summary of the analysis, click [here](bin/behavior.md)

## Cembrowski


This paper ([Cembrowski et al 2016](https://elifesciences.org/content/5/e14997#fig1s3)) is very similar to my experiment, so I want to compare the two. Like mine, they compare hippocampal gene expression from dorsal CA1, CA3, and DG sub regions. These cells were identifed through fac sorting to isolate genetically labeled CA1 and CA3 pyramical neurons and DG granular cells. 

Goals:
- [ ] I'd like to recreate the Cembrowski plots with the Cembrowski data and with my data.
- [x] Plot the Cembrowski data with my R scipts.
- [ ] I'd like to see how my data compares with the Cembrowski data. 


This data was made available here [open source data](https://www.janelia.org/lab/spruston-lab/resources/source-data-simulation-code-other-resources), but I downloaded it from the GenBank archive using the following commands: 

~~~~
wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74985/suppl/GSE74985_gene_exp.diff.gz'
wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74985/suppl/GSE74985_genes.fpkm_tracking.gz'
wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74985/suppl/GSE74985_genes.read_group_tracking.txt.gz'
wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74985/suppl/GSE74985_mergedCount.txt.gz'
gunzip *.gz
gzip GSE74985_genes.fpkm_tracking
~~~~

The GSE74985_genes.fpkm_tracking file is used to extract the gene names with the corresponding ensembl gene id. The file must be unzipped for use in R, but it must be zipped in order to store it on GitHub. The 4985_mergedCount.txt file is used for gene expression analyis in R.

For a summary of my current processing and analysis of the Cembrowski, click [here](bin/Cembrowski.md).

## Combo

I've made a stab at comparing their data to min and the results are [here](bin/combo.md).
