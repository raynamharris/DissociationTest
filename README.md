# Integrative analysis of the effects of biological and technical perturbations on hippocampal subfields transcriptomes
Rayna M. Harris, Hsin-Yi Kao, Juan Marcos Alarcon, Hans A. Hofmann, and André A. Fenton

## Overview
Cost-effective next-generation sequencing has made unbiased gene expression investigations possible. Gene expression studies at the level of single neurons may be especially important for understanding nervous system structure and function because of neuron-specific functionality and plasticity. While cellular dissociation is a prerequisite technical manipulation for such single-cell studies, the extent to which the process of dissociating cells affects neural gene expression has not been determined. Here, we examine the effect of cellular dissociation on gene expression in the mouse hippocampus. We also determine to which extent such changes might confound studies on the behavioral and physiological functions of hippocampus. 

## Approach

We first processed tissue punch samples from the dentate gyrus, CA3, and CA1 hippocampus subfields using either a tissue homogenization protocol or a cellular dissociation protocol. Then we used RNA sequencing to quantify sub-field specific gene expression differences between the preparation methods, while holding the tissue content constant within a region. To evaluate the impact of the tissue preparation method on the results of biological manipulations, we compared the gene expression differences to differences in homogenate tissue samples after mouse subjects were either exposed to a stressful experience or cognitive training, which are both common in behavioral research. Finally, we assessed the extent to which the subfield-specific gene expression patterns are consistent with those identified in a recently published hippocampus-subfield specific gene expression database. 

<img src="https://github.com/raynamharris/DissociationTest/blob/master/figures/00_methodsoverview/00_dissociationmethods-02.png" width="400px" align="middle"/>

## Repo Contents
- [**data**]((./data/)): contains all the raw and processed data files. They are broken up into sub folders. Some of the data was cleaned in using the process described in my other repo called [BehavEphysRNAseq](https://github.com/raynamharris/BehavEphyRNAseq)
- [**markdownfiles**](./markdownfiles/): this contains all the .R and .Rmd scripts as well as the .md output files. They have prefixes to hint at the order of operation. The workflow is described in more detail below
- [**figures**](./figures/): The output for all files from the Rmarkdown scripts

## Results

**Figure 1: The effects of cellular dissociation on hippocampal transcriptomes**  First, I compare CA1, CA3, and DG hippocampal samples from single individual that were prepared either by homogenization or dissociation. [Here are the analyses.](./markdownfiles/01_DissociationTest.Rmd) [Here are the results](./markdownfiles/01_DissociationTest.md)

**Figure 2: The effects of stressful experience on hippocampal transcriptomes**  Next, I look to see how gene expression varies across individuals that exposed to mild shocks, which are often used as behavioral stressor). [Here are the analyses.](./markdownfiles/02_StressTest.Rmd) [Here are the results](./markdownfiles/02_StressTest.md)

**Figure 3: The effects of cognitive training on hippocampal transcriptomes** Then, I identified signatures of genomic plasticity in response to cognitive training. [Here are the analyses.](./markdownfiles/03_CognitionTest.Rmd) [Here are the results](./markdownfiles/03_CognitionTest.md)

**Figure 4 and 5: Integrative analyses**. Then, I compare across experiments to look for shared and unique patterns of differential gene expression.
[Here are the analyses.](./markdownfiles/05_metaanlyses.Rmd) [Here are the results](./markdownfiles/05_metaanlyses.md).

**Figure 4 and 5: Gene Ontology (GO) analyses**. After identifying shared and unique patterns of differential gene expression, I conduct a GO analyses on these lists of genes
[Here are the analyses.](./markdownfiles/06_GO_MWU/06_GO_MWU.Rmd) [Here are the results](./markdownfiles/06_GO_MWU/06_GO_MWU.md).

**Supplementary analyses: Reproducing the Cembrowski results.** I reanalyze the [Cembrowski et al 2016](https://elifesciences.org/content/5/e14997#fig1s3) data, which has been used to create a database of sub-region specific hippocampalgene expression. There data and mine share some commonalities, so I wanted to compare the two. Like mine, they compare hippocampal gene expression from dorsal CA1, CA3, and DG sub regions. These cells were identifed through fac sorting to isolate genetically labeled CA1 and CA3 pyramical neurons and DG granular cells. [Here are the analyses.](./markdownfiles/04_Cembrowski.Rmd) [Here are the results](./markdownfiles/04_Cembrowski.md).

(Note: These figures are archived on [FigShare](https://figshare.com/articles/Integrative_analysis_of_hippocampal_genomic_plasticity/5116192).)

## Conclusions

<img src="https://github.com/raynamharris/DissociationTest/blob/master/figures/00_methodsoverview/00_dissociationmethods-01.png" width="400px" align="middle"/>

We find that 1% of the hippocampal transcriptome responds to the process of cellular dissociation and this cellular dissociation-induced response is distinct from the response of the hippocampus transcriptome to stress. There was some overlap in the cellular dissociation-induced response and the response to cognitive training. These findings of the concordant and discordant effects of technical and behavioral manipulations should inform the design of future neural transcriptome studies and thus facilitate a more comprehensive understanding of hippocampal function.