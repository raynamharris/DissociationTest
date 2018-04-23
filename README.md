Hippocampal transcriptomic responses to cellular dissociation
=============================================================

Rayna M. Harris, Hsin-Yi Kao, Juan Marcos Alarcon, Hans A. Hofmann, and
André A. Fenton

Overview
--------

Cost-effective next-generation sequencing has made unbiased gene
expression investigations possible. Gene expression studies at the level
of single neurons may be especially important for understanding nervous
system structure and function because of neuron-specific functionality
and plasticity. While cellular dissociation is a prerequisite technical
manipulation for such single-cell studies, the extent to which the
process of dissociating cells affects neural gene expression has not
been determined. Here, we examine the effect of cellular dissociation on
gene expression in the mouse hippocampus. We also determine to which
extent such changes might confound studies on the behavioral and
physiological functions of hippocampus.

Repo Contents
-------------

-   [**data**](./data/): contains most of the input processed data
    files. Large data fiels are stored in the Gene Expression Omnibus at
    [GSE99765](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99765)
    and
    [GSE100225](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100225).
    Raw kallisto abundance files are also stored in my other GitHub repo
    called
    [MouseHippocampusRNAseqData](https://github.com/raynamharris/MouseHippocampusRNAseqData).
-   [**UNIXworkflow**](./UNIXworkflow/): This descirbes the process I
    used to process my files using the Stampede Cluster at the Texas
    Advanced computing facility
-   [**scripts**](./scripts/): this contains all the .R and .Rmd scripts
    as well as the .md output files.
-   They have prefixes to hint at the order of operation.
-   The order was dramatically differnt when this work was first
    submitted for publication (version 1).
-   The current order is broken down by each figure.
-   [**figures**](./figures/): Contains all output for all files from
    the Rmarkdown scripts and my adobe-created images.

Current state of the analysis
-----------------------------

The following descriptions are not ready for publication, rather they
are converational descriptions of the current state.

<img src="./figures/fig_fig1.png" style="width:75.0%" />

**Figure 1.** General expression patterns show no major pattern of gene
expression alteration. 1A. Experimental design. 1B. Volcano plot showing
an exploratory(?) analysis of regional differences, for comparison to
published literature. 1C. Volcano plot showing a weak but asymmetric
response of gene expression to dissociation. 1D. PCA showing that no
major pattern of an effect of dissociation, but also showing that
regions are as clearly separated as one would expect given the
literature.

<table>
<thead>
<tr class="header">
<th align="center">Two-way contrast</th>
<th align="center">Up-regulated</th>
<th align="center">Down-regulated</th>
<th align="center">% DEGs</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center">CA1 vs. CA1</td>
<td align="center">222</td>
<td align="center">262</td>
<td align="center">2.9%</td>
</tr>
<tr class="even">
<td align="center">CA3 vs. DG</td>
<td align="center">45</td>
<td align="center">53</td>
<td align="center">0.5%</td>
</tr>
<tr class="odd">
<td align="center">CA1 vs. CA3</td>
<td align="center">17</td>
<td align="center">1</td>
<td align="center">0.1%</td>
</tr>
<tr class="even">
<td align="center">DISS vs. HOMO</td>
<td align="center">288</td>
<td align="center">56</td>
<td align="center">2.1%</td>
</tr>
</tbody>
</table>

**Table 1.** Number of differentially expressed genes by 2-way contrast.
"Down-regulated" means that expression is higher in the term on the far
left of column 1. "Up-regulated" means that expression is higher in the
term on the right of column 1. This shows that there is more variation
due to subfield than treatment.

![](./figures/fig_heatmapGO.png)

**Figure 2.** The gene list and go terms that everyone wants to know. 2A
The top 30 most differentially expressed genes. Genes are clustered by
correlation but samples are NOT clustered! Only "Jun" jumps out as a
gene related to learning and memory. 2B List of molecular function
categories that are enriched or depleted in the dissociated tissues
relative to controls. Again, nothing jumps out as classic memory
pathways, but there are some interesting affects on DNA regulation
(methylation, chromatin, histone, RNA binding, helicase) and metabolism
(oxidoreductase, cytokin, growth factors, ligase).

**Supplemental Table 1.** Expression level and fold change of of
significant genes (p &lt; 0.1) between dissociated tissue and
homogenized tissue. This table shows the log fold change (lfc), p-value
(padj), and direction of upregulation for each gene analyzed.

For now at:
<https://github.com/raynamharris/DissociationTest/blob/master/results/SuppTable1.csv>

    suptable <- read.csv("./results/SuppTable1.csv")
    tail(suptable, 10)

    ##         gene   lfc   padj upregulated.in
    ## 335     Rpsa  1.10 0.0959           DISS
    ## 336 Slc25a10  2.90 0.0959           DISS
    ## 337    Spry1  2.20 0.0959           DISS
    ## 338   Tango2  1.80 0.0959           DISS
    ## 339   Ubqln1 -1.10 0.0974           HOMO
    ## 340  Gadd45b  1.50 0.0984           DISS
    ## 341    Gsk3b -0.90 0.0984           HOMO
    ## 342     Atrx -1.10 0.0992           HOMO
    ## 343    Itpr3  3.10 0.0992           DISS
    ## 344    Gria2 -0.84 0.0997           HOMO

**Supplemental Table 2.** Gene ontologies of enriched genes. The first
row contains the GO category (either MF or CC). The second is the GO
term. Also shown are directionally, unumber of enriched genes in that
catory out of the total (ratio), and p-value.

<https://github.com/raynamharris/DissociationTest/blob/master/results/GOsignificantcatagories.csv>

    GOtable <- read.csv("./results/GOsignificantcatagories.csv")
    head(GOtable, 5)

    ##   GO                             GOterm enriched   ratio    pval
    ## 1 CC  cytosolic large ribosomal subunit     DISS   33/53 1.0e-15
    ## 2 CC                           ribosome     DISS  60/208 1.0e-15
    ## 3 CC               extracellular region     DISS 101/985 4.1e-15
    ## 4 CC                extracellular space     DISS 122/957 3.4e-14
    ## 5 CC                   nucleoplasm part     HOMO  38/529 6.3e-14

Outdated Video Abstract
-----------------------

I made this [short video](https://www.youtube.com/watch?v=taeAqimxXWo)
explaining how to use this GitHub repo when I submitted the first draft
to the journal Hippocampus and posted a pre-print on BioRxiv.

[![screenshot](./figures/screenshot.png)](https://www.youtube.com/watch?v=taeAqimxXWo)
