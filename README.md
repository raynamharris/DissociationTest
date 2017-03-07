# DissociationTest

The goal of this project is to determine the influence tissue processing techniques on hippocampal gene expression. This data was collected in collaboration with Andre Fenton and Maddy Kao. I compare hippocampal gene expression from CA1, CA3, and DG tissue samples. These samples were collected with a circular punch centered on the pyramidal and granular cells to enrich for that cell type. Biological replicates were tried in two ways before RNA extraction: one was dissociated before cell lysis while the other was homogenized before cell lysis. We wanted to identify the effects of treatment on gene expression. 

This ([Cembrowski et al 2016](https://elifesciences.org/content/5/e14997#fig1s3)) paper is very similar to my experiment, so I want to compare the two. Like mine, they compare hippocampal gene expression from dorsal CA1, CA3, and DG sub regions. These cells were identifed through fac sorting to isolate genetically labeled CA1 and CA3 pyramical neurons and DG granular cells. 

## Repo Contents
- [**data**]((./data/)): contains all the raw and processed data files. They are broken up into sub folders. Some of the data was cleaned in using the process described in my other repo called [BehavEphysRNAseq](https://github.com/raynamharris/BehavEphyRNAseq)
- [**markdownfiles**](./markdownfiles/): this contains all the .R and .Rmd scripts as well as the .md output files. They have prefixes to hint at the order of operation. The workflow is described in more detail below
- [**GO_MWU**](./GO_MWU/): This is a work in progress. I'm still trying to figure out how to perfect the go analysis.
- [**figures**](./figures/): The output for all files from the Rmarkdown scripts

## Workflow
All the data analyses for this project were conducted in R. Here is a brief overview of each of the .R or .Rmd files

[**01_DissociationTest**](./markdownfiles/01_DissociationTest.md)
First, I compare CA1 samples a single individual that were prepared either by homogenization or dissociation. [Here are the results](./markdownfiles/01_DissociationTest.md)

[**02_StressTest**](./markdownfiles/02_StressTest.md)
Next, I look to see how gene expression varies across individuals that were stress (as opposites to the tissue level stress examined above. [Here are the results](./markdownfiles/02_StressTest.md)

[**03_DissociationStressTest**](./markdownfiles/03_DissociationStressTest.md)
Then, I'm wondering what patterns hold up when I look at all these samples compbined.. [Here are the results](./markdownfiles/03_DissociationStressTest.md)

[**04_Cembrowski**](./markdownfiles/04_Cembrowski.Rmd)
Then I do some analyses of the cembrowski data. For a summary of my current processing and analysis of the Cembrowski, click [here](./markdownfiles/04_Cembrowski.Rmd).

[**05_Combo**](./markdownfiles/05_combo.Rmd)
I've made a stab at comparing their data to min and the results are [here](./markdownfiles/05_combo.Rmd).
