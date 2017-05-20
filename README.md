# DissociationTest
The objective of this study was to determine the effect of cellular dissociation on gene expression in the mouse hippocampus. To address this question we prepared hippocampal samples using a tissue homogenization protocol and a cell disassociation protocol. Then we used RNA sequencing to quantify gene  expression differences. To put the detected gene expression differences in context, we compared the gene expression differences in relation to two biological manipulations: one behavioral stressor and one cognitive manipulation. Finally, we asked the extent to which brain-region specific gene expression patterns are consistent with those identified in a recently published hippocampal gene expression database.  Our results demonstrate that 1\% of the genes are differentially expressed following chemical dissociation; however, the number of genes that distinguish hippocampal subregions is diminished. We show that the process of dissociation has broad effects on many diverse cellular and molecular processes while cognitive training influences very specific molecular and cellular processes. Finally, we show that genes related to calcium signaling are differentially expressed between subregions of the hippocampus in all four gene expression studies. Although numerous studies have examined gene expression in the hippocampus, little analytic attention has been paid to the systematic analysis of how the very process to the concordant and discordant effects that technical and behavioral manipulations have on transcriptional patterns. This research is an important step toward understanding molecular function in the hippocampus.

## Repo Contents
- [**data**]((./data/)): contains all the raw and processed data files. They are broken up into sub folders. Some of the data was cleaned in using the process described in my other repo called [BehavEphysRNAseq](https://github.com/raynamharris/BehavEphyRNAseq)
- [**markdownfiles**](./markdownfiles/): this contains all the .R and .Rmd scripts as well as the .md output files. They have prefixes to hint at the order of operation. The workflow is described in more detail below
- [**figures**](./figures/): The output for all files from the Rmarkdown scripts

## The experiments, analyses, and results
All the data analyses for this project were conducted in R. Here is a brief overview of each of the project. Below are links to the .Rmd and .md files.

<img src="https://github.com/raynamharris/DissociationTest/blob/master/figures/00_methodsoverview/00_dissociationmethods-01.png" width="400px" align="middle"/>

**Experiment 1: Effect of cellular dissociation.**  First, I compare CA1, CA3, and DG hippocampal samples from single individual that were prepared either by homogenization or dissociation. [Here are the analyses.](./markdownfiles/01_DissociationTest.Rmd) [Here are the results](./markdownfiles/01_DissociationTest.md)

**Experiment 2: Effect of an organismal stressor.**  Next, I look to see how gene expression varies across individuals that exposed to a mild shock (a potential behavioral stressor). [Here are the analyses.](./markdownfiles/02_StressTest.Rmd) [Here are the results](./markdownfiles/02_StressTest.md)

**Experiment 3: Effect of a cognitive task.** Then, I'm wondering what patterns hold up when I look at all these samples compbined. [Here are the analyses.](./markdownfiles/03_CognitionTest.Rmd) [Here are the results](./markdownfiles/03_CognitionTest.md)

**Experiment 4: Reproducing the Cembrowski results.** Then I reanalyze the [Cembrowski et al 2016](https://elifesciences.org/content/5/e14997#fig1s3) data, which has been used to create a database of sub-region specific hippocampalgene expression. There data and mine share some commonalities, so I wanted to compare the two. Like mine, they compare hippocampal gene expression from dorsal CA1, CA3, and DG sub regions. These cells were identifed through fac sorting to isolate genetically labeled CA1 and CA3 pyramical neurons and DG granular cells. [Here are the analyses.](./markdownfiles/04_Cembrowski.Rmd) [Here are the results](./markdownfiles/04_Cembrowski.md).

**Meta analyses**. Then, I compare across experiments to look for shared and unique patterns of differential gene expression.
[Here are the analyses.](./markdownfiles/05_metaanlyses.Rmd) [Here are the results](./markdownfiles/05_metaanlyses.md).

**Gene Ontology (GO) analyses**. After identifying shared and unique patterns of differential gene expression, I conduct a GO analyses on these lists of genes
[Here are the analyses.](./markdownfiles/06_GO_MWU/06_GO_MWU.Rmd) [Here are the results](./markdownfiles/06_GO_MWU/06_GO_MWU.md).
