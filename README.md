[![Binder](http://mybinder.org/badge.svg)](http://beta.mybinder.org/v2/gh/raynamharris/DissociationTest/master?urlpath=rstudio)
*Click the button to launch a Binder R session. Navigate to the
`scripts` directory and open any `.Rmd` file. Note: the first two are
long and slow. The rest are quick scripts that make figures.*

Hippocampal transcriptomic responses to cellular dissociation
=============================================================

About
-----

This repository contains the R scripts, data, results that for a study
about *Hippocampal transcriptomic responses to cellular dissociation*.
This research was submitted to the journal *Hippocampus* and was
accepted fro publication on March 15, 2019. A preprint is [available on
BioRxiv](https://www.biorxiv.org/content/early/2019/01/21/153585). The
authors on the manuscript are Rayna M. Harris, Hsin-Yi Kao, Juan Marcos
Alarcon, Hans A. Hofmann, and André A. Fenton

### Repo Contents

-   [**data**](./data/): contains most of the input processed data
    files. *Note: Raw reads and differential gene expression data are
    also available on Gene Expression Omnibus at
    [GSE99765](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99765).
    Raw kallisto abundance files are also stored another GitHub repo
    called
    [MouseHippocampusRNAseqData](https://github.com/raynamharris/MouseHippocampusRNAseqData).*
-   [**results**](./results): output files from the R scripts
-   [**UNIXworkflow**](./UNIXworkflow/): This descirbes the process I
    used to process my files using the Stampede Cluster at the Texas
    Advanced computing facility
-   [**scripts**](./scripts/): this contains all the .R and .Rmd scripts
    as well as the .md output files. They have prefixes to convey the
    order of operation.
-   [**figures**](./figures/): Contains all output for all files from
    the Rmarkdown scripts and my adobe-created images.

Abstract
--------

Single-neuron gene expression studies may be especially important for
understanding nervous system structure and function because of the
neuron-specific functionality and plasticity that defines functional
neural circuits. Cellular dissociation is a prerequisite technical
manipulation for single-cell and single cell-population studies, but the
extent to which the cellular dissociation process affects neural gene
expression has not been determined. This information is necessary for
interpreting the results of experimental manipulations that affect
neural function such as learning and memory. The goal of this research
was to determine the impact of chemical cell dissociation on brain
transcriptomes. We compared gene expression of microdissected samples
from the dentate gyrus (DG), CA3, and CA1 subfields of the mouse
hippocampus either prepared by a standard tissue homogenization protocol
or subjected to a chemical cellular dissociation procedure. We report
that compared to homogenization, chemical cellular dissociation alters
about 350 genes or 2% of the hippocampal transcriptome. While only a few
genes canonically implicated in long-term potentiation (LTP) and fear
memory change expression levels in response to the dissociation
procedure, these data indicate that sample preparation can affect gene
expression profiles, which might confound interpretation of results
depending on the research question. This study is important for the
investigation of any complex tissues as research effort moves from
subfield level analysis to single cell analysis of gene expression.

### General workflow and approach

1.  Experimental design (treatment \* hippocampal subfield)
2.  RNA-seq (Illumina, GSAF)
3.  Bioinformatics (TACC, FASTQC, cutadapt, kallisto)
4.  Data viz and stats (DESeq2, GOMWU, R)
5.  Version control and sharing (Git, GitHub, NCBI)

Figures and Tables
------------------

![](./figures/figure1.png)

**Figure 1. Experimental design and global expression gene expression
patterns. A)** Experimental design. Two tissue samples were taken from
three hippocampal subfields (CA1, CA3, and DG) from 300 um brain slices.
Two adjacent samples were processed using a homogenization (HOMO)
protocol or dissociated (DISS) before processing for tissue level gene
expression profiling. **B)** Dissociation does not yield
subfield-specific changes in gene expression between homogenized (HOMO,
open circles, dotted ellipse) and dissociated tissues (DISS, filled
circles, solid ellipse). PC1 accounts for 40% of all gene expression
variation and by inspection, separates the DG samples (orange circles)
from the CA1 (purple circles) and CA3 samples (green circles). PC2
accounts for 22% of the variation in gene expression and varies
significantly with treatment. The ellipses estimate the 95% confidence
interval for a multivariate t-distribution for homogenized (dashed line)
and dissociated (solid line) samples.

![](./figures/figure2.png)

**Figure 2. Enzymatic dissociation has a moderate effect on hippocampal
gene expression patterns compared to homogenized tissue. A)** Volcano
plot showing gene expression fold-difference and significance between
treatment groups. We found that 56 genes are up-regulated in the
homogenization control group (open circles) while 288 genes are
up-regulated in the dissociated treatment group (filled dark grey
circles). Genes below the p-value &lt; 0.1 (or –log p-value &lt; 1) are
shown in light grey. **B)** Heatmap showing the top 30 differentially
expressed genes between dissociated and homogenized tissue. Square boxes
at the top are color coded by sample (white: homogenized, grey:
dissociated, purple: CA1, green: CA3, orange: DG. Within the heatmap,
log fold difference levels of expression are indicated by the
blue-green-yellow gradient with lighter colors indicating increased
expression.

<table>
<thead>
<tr class="header">
<th style="text-align: center;">Two-way contrast</th>
<th style="text-align: center;">Up-regulated</th>
<th style="text-align: center;">Down-regulated</th>
<th style="text-align: center;">% DEGs</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: center;">CA1 vs. CA1</td>
<td style="text-align: center;">222</td>
<td style="text-align: center;">262</td>
<td style="text-align: center;">2.9%</td>
</tr>
<tr class="even">
<td style="text-align: center;">CA3 vs. DG</td>
<td style="text-align: center;">45</td>
<td style="text-align: center;">53</td>
<td style="text-align: center;">0.5%</td>
</tr>
<tr class="odd">
<td style="text-align: center;">CA1 vs. CA3</td>
<td style="text-align: center;">17</td>
<td style="text-align: center;">1</td>
<td style="text-align: center;">0.1%</td>
</tr>
<tr class="even">
<td style="text-align: center;">DISS vs. HOMO</td>
<td style="text-align: center;">288</td>
<td style="text-align: center;">56</td>
<td style="text-align: center;">2.1%</td>
</tr>
</tbody>
</table>

**Table 1. Differentially expressed genes by subfield and treatment.**
The total number and percent of differentially expressed genes (DEGs)
for four two-way contrasts were calculated using DESeq2. Increased
expression cutoffs are defined as log fold-change &gt; 0; p &lt; 0.1
while decreased expression is defined as log fold-change &lt; 0; p &lt;
0.1. % DEGs/Total: The sum of up and down regulated genes divided by the
total number of genes analyzed (16,709) multiplied by 100%. This table
shows that differences between dissociated (DISS) tissue and homogenized
(HOMO) tissues are on the same scale as those between the CA1 and DG
subfields of the hippocampus.

**Preview of Supplemental Table 1. Expression level and fold change of
significant genes (p &lt; 0.1) between dissociated tissue and
homogenized tissue.** This table shows the log fold change (lfc),
p-value (padj), and direction of upregulation for each gene analyzed.
*Full table available at
<a href="https://github.com/raynamharris/DissociationTest/blob/master/results/dissociationDEGs.csv" class="uri">https://github.com/raynamharris/DissociationTest/blob/master/results/dissociationDEGs.csv</a>.*

<table>
<thead>
<tr>
<th style="text-align:left;">
gene
</th>
<th style="text-align:right;">
lfc
</th>
<th style="text-align:right;">
padj
</th>
<th style="text-align:left;">
direction
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Trf
</td>
<td style="text-align:right;">
2.72
</td>
<td style="text-align:right;">
5.00e-07
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Hexb
</td>
<td style="text-align:right;">
2.35
</td>
<td style="text-align:right;">
8.00e-07
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Selplg
</td>
<td style="text-align:right;">
2.97
</td>
<td style="text-align:right;">
9.00e-07
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
C1qb
</td>
<td style="text-align:right;">
2.28
</td>
<td style="text-align:right;">
7.10e-06
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Csf1r
</td>
<td style="text-align:right;">
2.13
</td>
<td style="text-align:right;">
9.60e-06
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ctss
</td>
<td style="text-align:right;">
2.59
</td>
<td style="text-align:right;">
9.60e-06
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cnp
</td>
<td style="text-align:right;">
2.45
</td>
<td style="text-align:right;">
4.48e-05
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Il1a
</td>
<td style="text-align:right;">
3.06
</td>
<td style="text-align:right;">
4.48e-05
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Mag
</td>
<td style="text-align:right;">
3.31
</td>
<td style="text-align:right;">
4.48e-05
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cd14
</td>
<td style="text-align:right;">
3.38
</td>
<td style="text-align:right;">
4.88e-05
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Mpeg1
</td>
<td style="text-align:right;">
2.42
</td>
<td style="text-align:right;">
6.47e-05
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Tmem88b
</td>
<td style="text-align:right;">
3.14
</td>
<td style="text-align:right;">
6.78e-05
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Nfkbia
</td>
<td style="text-align:right;">
2.10
</td>
<td style="text-align:right;">
6.84e-05
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Slc15a3
</td>
<td style="text-align:right;">
3.45
</td>
<td style="text-align:right;">
9.32e-05
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cd83
</td>
<td style="text-align:right;">
2.58
</td>
<td style="text-align:right;">
9.99e-05
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cldn11
</td>
<td style="text-align:right;">
3.14
</td>
<td style="text-align:right;">
1.05e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Laptm5
</td>
<td style="text-align:right;">
2.31
</td>
<td style="text-align:right;">
2.24e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Plau
</td>
<td style="text-align:right;">
3.90
</td>
<td style="text-align:right;">
2.24e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Plp1
</td>
<td style="text-align:right;">
2.71
</td>
<td style="text-align:right;">
2.24e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Mal
</td>
<td style="text-align:right;">
3.20
</td>
<td style="text-align:right;">
2.32e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Plek
</td>
<td style="text-align:right;">
2.50
</td>
<td style="text-align:right;">
2.32e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Jun
</td>
<td style="text-align:right;">
1.60
</td>
<td style="text-align:right;">
3.20e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Mobp
</td>
<td style="text-align:right;">
2.60
</td>
<td style="text-align:right;">
4.41e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Pcsk1n
</td>
<td style="text-align:right;">
1.30
</td>
<td style="text-align:right;">
4.80e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Socs3
</td>
<td style="text-align:right;">
2.26
</td>
<td style="text-align:right;">
4.80e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Pllp
</td>
<td style="text-align:right;">
2.76
</td>
<td style="text-align:right;">
5.46e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ndrg1
</td>
<td style="text-align:right;">
2.87
</td>
<td style="text-align:right;">
5.80e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Phldb1
</td>
<td style="text-align:right;">
2.27
</td>
<td style="text-align:right;">
6.85e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Tyrobp
</td>
<td style="text-align:right;">
3.34
</td>
<td style="text-align:right;">
8.22e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cx3cr1
</td>
<td style="text-align:right;">
2.10
</td>
<td style="text-align:right;">
9.50e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Tmem119
</td>
<td style="text-align:right;">
2.54
</td>
<td style="text-align:right;">
9.50e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
C1qa
</td>
<td style="text-align:right;">
2.32
</td>
<td style="text-align:right;">
9.85e-04
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Slc2a5
</td>
<td style="text-align:right;">
4.08
</td>
<td style="text-align:right;">
1.23e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Fcrls
</td>
<td style="text-align:right;">
2.67
</td>
<td style="text-align:right;">
1.40e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Lgmn
</td>
<td style="text-align:right;">
1.84
</td>
<td style="text-align:right;">
1.57e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Clic4
</td>
<td style="text-align:right;">
2.01
</td>
<td style="text-align:right;">
1.81e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Sema5a
</td>
<td style="text-align:right;">
3.25
</td>
<td style="text-align:right;">
2.98e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
C5ar1
</td>
<td style="text-align:right;">
3.66
</td>
<td style="text-align:right;">
3.01e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cmtm5
</td>
<td style="text-align:right;">
2.66
</td>
<td style="text-align:right;">
3.01e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ctcfl
</td>
<td style="text-align:right;">
3.08
</td>
<td style="text-align:right;">
3.05e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cryab
</td>
<td style="text-align:right;">
2.81
</td>
<td style="text-align:right;">
3.07e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rpl26
</td>
<td style="text-align:right;">
1.27
</td>
<td style="text-align:right;">
3.22e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Gab1
</td>
<td style="text-align:right;">
2.90
</td>
<td style="text-align:right;">
3.44e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
C1qc
</td>
<td style="text-align:right;">
2.60
</td>
<td style="text-align:right;">
4.40e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Fa2h
</td>
<td style="text-align:right;">
3.30
</td>
<td style="text-align:right;">
4.61e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Clec2l
</td>
<td style="text-align:right;">
1.94
</td>
<td style="text-align:right;">
5.01e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Spry2
</td>
<td style="text-align:right;">
1.79
</td>
<td style="text-align:right;">
5.08e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Mcam
</td>
<td style="text-align:right;">
2.86
</td>
<td style="text-align:right;">
5.12e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Olfml3
</td>
<td style="text-align:right;">
2.09
</td>
<td style="text-align:right;">
5.18e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Sepp1
</td>
<td style="text-align:right;">
1.54
</td>
<td style="text-align:right;">
5.78e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ptma
</td>
<td style="text-align:right;">
1.19
</td>
<td style="text-align:right;">
6.21e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Dusp1
</td>
<td style="text-align:right;">
1.35
</td>
<td style="text-align:right;">
8.03e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Lag3
</td>
<td style="text-align:right;">
3.32
</td>
<td style="text-align:right;">
8.03e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Mbp
</td>
<td style="text-align:right;">
1.95
</td>
<td style="text-align:right;">
8.03e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Atf3
</td>
<td style="text-align:right;">
2.17
</td>
<td style="text-align:right;">
8.09e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
C1ql1
</td>
<td style="text-align:right;">
3.10
</td>
<td style="text-align:right;">
8.09e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Adgrg1
</td>
<td style="text-align:right;">
1.86
</td>
<td style="text-align:right;">
8.62e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Csrp1
</td>
<td style="text-align:right;">
1.70
</td>
<td style="text-align:right;">
9.07e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ltbp3
</td>
<td style="text-align:right;">
2.29
</td>
<td style="text-align:right;">
9.07e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Mcl1
</td>
<td style="text-align:right;">
1.41
</td>
<td style="text-align:right;">
9.07e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rpl29
</td>
<td style="text-align:right;">
1.16
</td>
<td style="text-align:right;">
9.07e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rpl36a
</td>
<td style="text-align:right;">
1.57
</td>
<td style="text-align:right;">
9.07e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rps27
</td>
<td style="text-align:right;">
1.29
</td>
<td style="text-align:right;">
9.59e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ighv14-2
</td>
<td style="text-align:right;">
7.81
</td>
<td style="text-align:right;">
9.60e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Mapk3
</td>
<td style="text-align:right;">
1.61
</td>
<td style="text-align:right;">
9.80e-03
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Trem2
</td>
<td style="text-align:right;">
3.25
</td>
<td style="text-align:right;">
1.14e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Syngr2
</td>
<td style="text-align:right;">
3.35
</td>
<td style="text-align:right;">
1.20e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Irf8
</td>
<td style="text-align:right;">
2.77
</td>
<td style="text-align:right;">
1.26e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Jund
</td>
<td style="text-align:right;">
1.06
</td>
<td style="text-align:right;">
1.26e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Plekhb1
</td>
<td style="text-align:right;">
1.49
</td>
<td style="text-align:right;">
1.26e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cst3
</td>
<td style="text-align:right;">
1.53
</td>
<td style="text-align:right;">
1.45e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rpl28
</td>
<td style="text-align:right;">
1.17
</td>
<td style="text-align:right;">
1.55e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Unc93b1
</td>
<td style="text-align:right;">
2.28
</td>
<td style="text-align:right;">
1.55e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
D6Wsu163e
</td>
<td style="text-align:right;">
2.32
</td>
<td style="text-align:right;">
1.59e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Sstr1
</td>
<td style="text-align:right;">
4.26
</td>
<td style="text-align:right;">
1.59e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Gm10401
</td>
<td style="text-align:right;">
7.62
</td>
<td style="text-align:right;">
1.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Itga5
</td>
<td style="text-align:right;">
3.05
</td>
<td style="text-align:right;">
1.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ly86
</td>
<td style="text-align:right;">
2.81
</td>
<td style="text-align:right;">
1.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Olig1
</td>
<td style="text-align:right;">
1.67
</td>
<td style="text-align:right;">
1.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Zfp36
</td>
<td style="text-align:right;">
1.79
</td>
<td style="text-align:right;">
1.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rps8
</td>
<td style="text-align:right;">
1.04
</td>
<td style="text-align:right;">
1.70e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cd9
</td>
<td style="text-align:right;">
2.36
</td>
<td style="text-align:right;">
1.72e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rph3a
</td>
<td style="text-align:right;">
2.26
</td>
<td style="text-align:right;">
1.72e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Slc2a1
</td>
<td style="text-align:right;">
1.70
</td>
<td style="text-align:right;">
1.77e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Bin2
</td>
<td style="text-align:right;">
2.49
</td>
<td style="text-align:right;">
1.83e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Gpr84
</td>
<td style="text-align:right;">
2.86
</td>
<td style="text-align:right;">
1.87e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
H2-K1
</td>
<td style="text-align:right;">
2.09
</td>
<td style="text-align:right;">
1.87e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Slco2b1
</td>
<td style="text-align:right;">
2.17
</td>
<td style="text-align:right;">
1.87e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Tgfbr2
</td>
<td style="text-align:right;">
2.44
</td>
<td style="text-align:right;">
1.87e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cdh9
</td>
<td style="text-align:right;">
3.65
</td>
<td style="text-align:right;">
1.99e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ctsd
</td>
<td style="text-align:right;">
1.26
</td>
<td style="text-align:right;">
1.99e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ptgds
</td>
<td style="text-align:right;">
2.53
</td>
<td style="text-align:right;">
1.99e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Mt1
</td>
<td style="text-align:right;">
1.42
</td>
<td style="text-align:right;">
2.01e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Itgb5
</td>
<td style="text-align:right;">
1.98
</td>
<td style="text-align:right;">
2.01e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Inpp5d
</td>
<td style="text-align:right;">
2.84
</td>
<td style="text-align:right;">
2.04e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Icam1
</td>
<td style="text-align:right;">
2.87
</td>
<td style="text-align:right;">
2.06e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Myrf
</td>
<td style="text-align:right;">
2.21
</td>
<td style="text-align:right;">
2.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Sparc
</td>
<td style="text-align:right;">
1.79
</td>
<td style="text-align:right;">
2.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Tnf
</td>
<td style="text-align:right;">
2.40
</td>
<td style="text-align:right;">
2.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Siglech
</td>
<td style="text-align:right;">
2.11
</td>
<td style="text-align:right;">
2.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Serpinb1a
</td>
<td style="text-align:right;">
4.36
</td>
<td style="text-align:right;">
2.26e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
mt-Co3
</td>
<td style="text-align:right;">
1.40
</td>
<td style="text-align:right;">
2.26e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Mog
</td>
<td style="text-align:right;">
2.48
</td>
<td style="text-align:right;">
2.27e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Stat3
</td>
<td style="text-align:right;">
1.56
</td>
<td style="text-align:right;">
2.41e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Tlr7
</td>
<td style="text-align:right;">
3.89
</td>
<td style="text-align:right;">
2.42e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Lgals9
</td>
<td style="text-align:right;">
2.74
</td>
<td style="text-align:right;">
2.62e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Tln1
</td>
<td style="text-align:right;">
1.34
</td>
<td style="text-align:right;">
2.62e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cyba
</td>
<td style="text-align:right;">
3.80
</td>
<td style="text-align:right;">
2.63e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ermn
</td>
<td style="text-align:right;">
2.53
</td>
<td style="text-align:right;">
2.69e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Fn1
</td>
<td style="text-align:right;">
1.77
</td>
<td style="text-align:right;">
2.69e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
mt-Co1
</td>
<td style="text-align:right;">
1.24
</td>
<td style="text-align:right;">
2.69e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Sipa1
</td>
<td style="text-align:right;">
2.71
</td>
<td style="text-align:right;">
2.69e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Vasp
</td>
<td style="text-align:right;">
2.31
</td>
<td style="text-align:right;">
2.73e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Shroom3
</td>
<td style="text-align:right;">
6.23
</td>
<td style="text-align:right;">
2.80e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
mt-Cytb
</td>
<td style="text-align:right;">
1.27
</td>
<td style="text-align:right;">
2.86e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Qdpr
</td>
<td style="text-align:right;">
1.55
</td>
<td style="text-align:right;">
2.86e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Slc12a2
</td>
<td style="text-align:right;">
2.07
</td>
<td style="text-align:right;">
2.86e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Tmem170
</td>
<td style="text-align:right;">
2.53
</td>
<td style="text-align:right;">
2.93e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Apod
</td>
<td style="text-align:right;">
2.05
</td>
<td style="text-align:right;">
3.00e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ccl4
</td>
<td style="text-align:right;">
1.94
</td>
<td style="text-align:right;">
3.00e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Erdr1
</td>
<td style="text-align:right;">
1.33
</td>
<td style="text-align:right;">
3.00e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Il1b
</td>
<td style="text-align:right;">
2.41
</td>
<td style="text-align:right;">
3.00e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Plekhg3
</td>
<td style="text-align:right;">
2.65
</td>
<td style="text-align:right;">
3.00e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Sept4
</td>
<td style="text-align:right;">
1.84
</td>
<td style="text-align:right;">
3.00e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cd82
</td>
<td style="text-align:right;">
2.18
</td>
<td style="text-align:right;">
3.20e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cspg5
</td>
<td style="text-align:right;">
1.28
</td>
<td style="text-align:right;">
3.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ucp2
</td>
<td style="text-align:right;">
2.49
</td>
<td style="text-align:right;">
3.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Smoc2
</td>
<td style="text-align:right;">
3.81
</td>
<td style="text-align:right;">
3.27e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
mt-Nd3
</td>
<td style="text-align:right;">
1.70
</td>
<td style="text-align:right;">
3.28e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Fosb
</td>
<td style="text-align:right;">
1.59
</td>
<td style="text-align:right;">
3.32e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Mfng
</td>
<td style="text-align:right;">
4.05
</td>
<td style="text-align:right;">
3.32e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Pcdhga8
</td>
<td style="text-align:right;">
1.42
</td>
<td style="text-align:right;">
3.32e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Emp2
</td>
<td style="text-align:right;">
2.35
</td>
<td style="text-align:right;">
3.33e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cd37
</td>
<td style="text-align:right;">
4.14
</td>
<td style="text-align:right;">
3.36e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Scrg1
</td>
<td style="text-align:right;">
2.63
</td>
<td style="text-align:right;">
3.36e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Prr18
</td>
<td style="text-align:right;">
2.29
</td>
<td style="text-align:right;">
3.57e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
mt-Nd1
</td>
<td style="text-align:right;">
1.21
</td>
<td style="text-align:right;">
3.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cd274
</td>
<td style="text-align:right;">
4.24
</td>
<td style="text-align:right;">
4.12e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rpl32
</td>
<td style="text-align:right;">
1.33
</td>
<td style="text-align:right;">
4.12e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Synpr
</td>
<td style="text-align:right;">
2.31
</td>
<td style="text-align:right;">
4.12e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Pold1
</td>
<td style="text-align:right;">
2.23
</td>
<td style="text-align:right;">
4.17e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Btg2
</td>
<td style="text-align:right;">
1.39
</td>
<td style="text-align:right;">
4.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Crtc2
</td>
<td style="text-align:right;">
1.78
</td>
<td style="text-align:right;">
4.37e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Msn
</td>
<td style="text-align:right;">
1.75
</td>
<td style="text-align:right;">
4.55e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
mt-Nd4
</td>
<td style="text-align:right;">
1.15
</td>
<td style="text-align:right;">
4.55e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rtn2
</td>
<td style="text-align:right;">
1.47
</td>
<td style="text-align:right;">
4.55e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Myc
</td>
<td style="text-align:right;">
2.35
</td>
<td style="text-align:right;">
4.56e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Bcas1
</td>
<td style="text-align:right;">
1.70
</td>
<td style="text-align:right;">
4.70e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Gatm
</td>
<td style="text-align:right;">
2.12
</td>
<td style="text-align:right;">
4.74e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Gjb1
</td>
<td style="text-align:right;">
4.76
</td>
<td style="text-align:right;">
4.74e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Pros1
</td>
<td style="text-align:right;">
2.29
</td>
<td style="text-align:right;">
4.74e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Psat1
</td>
<td style="text-align:right;">
1.31
</td>
<td style="text-align:right;">
4.74e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rab3il1
</td>
<td style="text-align:right;">
2.76
</td>
<td style="text-align:right;">
4.78e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Arl4c
</td>
<td style="text-align:right;">
1.79
</td>
<td style="text-align:right;">
4.89e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Dhrs3
</td>
<td style="text-align:right;">
2.95
</td>
<td style="text-align:right;">
4.89e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Enpp2
</td>
<td style="text-align:right;">
1.77
</td>
<td style="text-align:right;">
4.89e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Gpr34
</td>
<td style="text-align:right;">
2.08
</td>
<td style="text-align:right;">
4.89e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ntng1
</td>
<td style="text-align:right;">
3.08
</td>
<td style="text-align:right;">
4.89e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
P2ry12
</td>
<td style="text-align:right;">
1.74
</td>
<td style="text-align:right;">
4.89e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Tlr2
</td>
<td style="text-align:right;">
2.65
</td>
<td style="text-align:right;">
4.89e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cyp27a1
</td>
<td style="text-align:right;">
3.58
</td>
<td style="text-align:right;">
4.92e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Gna12
</td>
<td style="text-align:right;">
1.30
</td>
<td style="text-align:right;">
4.92e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Asap3
</td>
<td style="text-align:right;">
3.89
</td>
<td style="text-align:right;">
4.96e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Itgam
</td>
<td style="text-align:right;">
1.75
</td>
<td style="text-align:right;">
4.96e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Lcp1
</td>
<td style="text-align:right;">
2.73
</td>
<td style="text-align:right;">
4.96e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Frmd8
</td>
<td style="text-align:right;">
2.79
</td>
<td style="text-align:right;">
5.03e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ncf1
</td>
<td style="text-align:right;">
1.80
</td>
<td style="text-align:right;">
5.07e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Dnajb2
</td>
<td style="text-align:right;">
1.52
</td>
<td style="text-align:right;">
5.14e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
mt-Atp6
</td>
<td style="text-align:right;">
1.25
</td>
<td style="text-align:right;">
5.15e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
mt-Co2
</td>
<td style="text-align:right;">
1.36
</td>
<td style="text-align:right;">
5.15e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Wasf2
</td>
<td style="text-align:right;">
1.67
</td>
<td style="text-align:right;">
5.20e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ccdc3
</td>
<td style="text-align:right;">
5.38
</td>
<td style="text-align:right;">
5.32e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Smox
</td>
<td style="text-align:right;">
1.58
</td>
<td style="text-align:right;">
5.40e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Tmem63a
</td>
<td style="text-align:right;">
2.46
</td>
<td style="text-align:right;">
5.40e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Vsir
</td>
<td style="text-align:right;">
2.06
</td>
<td style="text-align:right;">
5.40e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Gpr17
</td>
<td style="text-align:right;">
2.60
</td>
<td style="text-align:right;">
5.46e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cmtm6
</td>
<td style="text-align:right;">
1.70
</td>
<td style="text-align:right;">
5.75e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Dock8
</td>
<td style="text-align:right;">
3.07
</td>
<td style="text-align:right;">
5.75e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Gp1bb
</td>
<td style="text-align:right;">
3.10
</td>
<td style="text-align:right;">
5.75e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ier2
</td>
<td style="text-align:right;">
1.33
</td>
<td style="text-align:right;">
5.75e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Il16
</td>
<td style="text-align:right;">
3.32
</td>
<td style="text-align:right;">
5.75e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rftn1
</td>
<td style="text-align:right;">
3.55
</td>
<td style="text-align:right;">
5.75e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Sox10
</td>
<td style="text-align:right;">
2.16
</td>
<td style="text-align:right;">
5.75e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Junb
</td>
<td style="text-align:right;">
1.03
</td>
<td style="text-align:right;">
5.81e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
D630003M21Rik
</td>
<td style="text-align:right;">
3.12
</td>
<td style="text-align:right;">
5.89e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Kcnk6
</td>
<td style="text-align:right;">
2.84
</td>
<td style="text-align:right;">
6.10e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rela
</td>
<td style="text-align:right;">
1.52
</td>
<td style="text-align:right;">
6.18e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
mt-Nd4l
</td>
<td style="text-align:right;">
1.06
</td>
<td style="text-align:right;">
6.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Itgb4
</td>
<td style="text-align:right;">
2.93
</td>
<td style="text-align:right;">
6.30e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Zfp521
</td>
<td style="text-align:right;">
2.72
</td>
<td style="text-align:right;">
6.30e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ilkap
</td>
<td style="text-align:right;">
1.61
</td>
<td style="text-align:right;">
6.31e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Nckap1l
</td>
<td style="text-align:right;">
2.11
</td>
<td style="text-align:right;">
6.31e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ccrl2
</td>
<td style="text-align:right;">
3.01
</td>
<td style="text-align:right;">
6.48e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Gjc3
</td>
<td style="text-align:right;">
1.85
</td>
<td style="text-align:right;">
6.48e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Csrnp1
</td>
<td style="text-align:right;">
1.52
</td>
<td style="text-align:right;">
6.62e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Klf4
</td>
<td style="text-align:right;">
1.74
</td>
<td style="text-align:right;">
6.62e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Tbc1d16
</td>
<td style="text-align:right;">
1.48
</td>
<td style="text-align:right;">
6.65e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Wscd1
</td>
<td style="text-align:right;">
1.50
</td>
<td style="text-align:right;">
6.69e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Atf5
</td>
<td style="text-align:right;">
1.58
</td>
<td style="text-align:right;">
6.85e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Dusp10
</td>
<td style="text-align:right;">
2.49
</td>
<td style="text-align:right;">
6.85e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Prdm5
</td>
<td style="text-align:right;">
3.59
</td>
<td style="text-align:right;">
6.85e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Arpc1b
</td>
<td style="text-align:right;">
2.05
</td>
<td style="text-align:right;">
6.89e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
mt-Nd2
</td>
<td style="text-align:right;">
1.19
</td>
<td style="text-align:right;">
6.95e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Grap
</td>
<td style="text-align:right;">
2.96
</td>
<td style="text-align:right;">
7.08e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rpl18
</td>
<td style="text-align:right;">
1.04
</td>
<td style="text-align:right;">
7.08e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Suv420h2
</td>
<td style="text-align:right;">
2.25
</td>
<td style="text-align:right;">
7.08e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Il6ra
</td>
<td style="text-align:right;">
2.21
</td>
<td style="text-align:right;">
7.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Npr1
</td>
<td style="text-align:right;">
3.34
</td>
<td style="text-align:right;">
7.24e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Smtnl2
</td>
<td style="text-align:right;">
5.83
</td>
<td style="text-align:right;">
7.25e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
mt-Atp8
</td>
<td style="text-align:right;">
1.11
</td>
<td style="text-align:right;">
7.37e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Piezo1
</td>
<td style="text-align:right;">
3.06
</td>
<td style="text-align:right;">
7.37e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rpl23a
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
7.37e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rpl31
</td>
<td style="text-align:right;">
1.09
</td>
<td style="text-align:right;">
7.37e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Pim1
</td>
<td style="text-align:right;">
2.41
</td>
<td style="text-align:right;">
7.37e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cd63
</td>
<td style="text-align:right;">
1.58
</td>
<td style="text-align:right;">
7.40e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Csf1
</td>
<td style="text-align:right;">
1.88
</td>
<td style="text-align:right;">
7.40e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rplp1
</td>
<td style="text-align:right;">
1.13
</td>
<td style="text-align:right;">
7.40e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Myo1f
</td>
<td style="text-align:right;">
4.78
</td>
<td style="text-align:right;">
7.47e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Insig1
</td>
<td style="text-align:right;">
1.61
</td>
<td style="text-align:right;">
7.50e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Itpk1
</td>
<td style="text-align:right;">
1.53
</td>
<td style="text-align:right;">
7.58e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Lhfpl2
</td>
<td style="text-align:right;">
2.15
</td>
<td style="text-align:right;">
7.69e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Chadl
</td>
<td style="text-align:right;">
1.82
</td>
<td style="text-align:right;">
7.69e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Nlrp3
</td>
<td style="text-align:right;">
2.20
</td>
<td style="text-align:right;">
7.69e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Sh3gl1
</td>
<td style="text-align:right;">
1.56
</td>
<td style="text-align:right;">
7.69e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cidea
</td>
<td style="text-align:right;">
5.46
</td>
<td style="text-align:right;">
7.85e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rack1
</td>
<td style="text-align:right;">
1.03
</td>
<td style="text-align:right;">
7.85e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Sema4g
</td>
<td style="text-align:right;">
2.30
</td>
<td style="text-align:right;">
7.95e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rpl23
</td>
<td style="text-align:right;">
1.23
</td>
<td style="text-align:right;">
8.02e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Adam17
</td>
<td style="text-align:right;">
2.27
</td>
<td style="text-align:right;">
8.05e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cpm
</td>
<td style="text-align:right;">
2.85
</td>
<td style="text-align:right;">
8.15e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cox4i1
</td>
<td style="text-align:right;">
1.04
</td>
<td style="text-align:right;">
8.27e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Nectin1
</td>
<td style="text-align:right;">
1.84
</td>
<td style="text-align:right;">
8.27e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ostf1
</td>
<td style="text-align:right;">
2.36
</td>
<td style="text-align:right;">
8.27e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Maml2
</td>
<td style="text-align:right;">
2.19
</td>
<td style="text-align:right;">
8.28e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ap1s1
</td>
<td style="text-align:right;">
1.07
</td>
<td style="text-align:right;">
8.38e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Fam124a
</td>
<td style="text-align:right;">
2.36
</td>
<td style="text-align:right;">
8.38e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Hr
</td>
<td style="text-align:right;">
2.45
</td>
<td style="text-align:right;">
8.38e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Nfatc1
</td>
<td style="text-align:right;">
3.56
</td>
<td style="text-align:right;">
8.38e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rps6kb2
</td>
<td style="text-align:right;">
1.65
</td>
<td style="text-align:right;">
8.38e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Stk32b
</td>
<td style="text-align:right;">
5.43
</td>
<td style="text-align:right;">
8.39e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Unc5b
</td>
<td style="text-align:right;">
1.96
</td>
<td style="text-align:right;">
8.52e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cdkn1a
</td>
<td style="text-align:right;">
2.51
</td>
<td style="text-align:right;">
8.54e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Fam102b
</td>
<td style="text-align:right;">
1.53
</td>
<td style="text-align:right;">
8.75e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Sh3d19
</td>
<td style="text-align:right;">
2.43
</td>
<td style="text-align:right;">
9.07e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Tspan2
</td>
<td style="text-align:right;">
1.93
</td>
<td style="text-align:right;">
9.07e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Tnfaip6
</td>
<td style="text-align:right;">
2.76
</td>
<td style="text-align:right;">
9.09e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cd68
</td>
<td style="text-align:right;">
2.35
</td>
<td style="text-align:right;">
9.11e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Mmp15
</td>
<td style="text-align:right;">
1.94
</td>
<td style="text-align:right;">
9.19e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Lamb2
</td>
<td style="text-align:right;">
2.70
</td>
<td style="text-align:right;">
9.21e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Itpkb
</td>
<td style="text-align:right;">
1.53
</td>
<td style="text-align:right;">
9.25e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Il10ra
</td>
<td style="text-align:right;">
2.49
</td>
<td style="text-align:right;">
9.29e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rps27a
</td>
<td style="text-align:right;">
1.13
</td>
<td style="text-align:right;">
9.33e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rpl35a
</td>
<td style="text-align:right;">
1.47
</td>
<td style="text-align:right;">
9.45e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Tbx1
</td>
<td style="text-align:right;">
5.68
</td>
<td style="text-align:right;">
9.56e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Srgn
</td>
<td style="text-align:right;">
2.41
</td>
<td style="text-align:right;">
9.59e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Apoe
</td>
<td style="text-align:right;">
1.25
</td>
<td style="text-align:right;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Bche
</td>
<td style="text-align:right;">
3.35
</td>
<td style="text-align:right;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ccl9
</td>
<td style="text-align:right;">
3.31
</td>
<td style="text-align:right;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cd81
</td>
<td style="text-align:right;">
1.03
</td>
<td style="text-align:right;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Cfh
</td>
<td style="text-align:right;">
2.03
</td>
<td style="text-align:right;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Ctsc
</td>
<td style="text-align:right;">
2.40
</td>
<td style="text-align:right;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Dusp26
</td>
<td style="text-align:right;">
2.02
</td>
<td style="text-align:right;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Erbin
</td>
<td style="text-align:right;">
1.72
</td>
<td style="text-align:right;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Gjc2
</td>
<td style="text-align:right;">
2.39
</td>
<td style="text-align:right;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Lrtm2
</td>
<td style="text-align:right;">
1.80
</td>
<td style="text-align:right;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Plin4
</td>
<td style="text-align:right;">
4.05
</td>
<td style="text-align:right;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rgs10
</td>
<td style="text-align:right;">
2.18
</td>
<td style="text-align:right;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rnaset2b
</td>
<td style="text-align:right;">
2.00
</td>
<td style="text-align:right;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Rpsa
</td>
<td style="text-align:right;">
1.08
</td>
<td style="text-align:right;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Slc25a10
</td>
<td style="text-align:right;">
2.94
</td>
<td style="text-align:right;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Spry1
</td>
<td style="text-align:right;">
2.18
</td>
<td style="text-align:right;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Tango2
</td>
<td style="text-align:right;">
1.80
</td>
<td style="text-align:right;">
9.60e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Gadd45b
</td>
<td style="text-align:right;">
1.49
</td>
<td style="text-align:right;">
9.83e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Itpr3
</td>
<td style="text-align:right;">
3.12
</td>
<td style="text-align:right;">
9.92e-02
</td>
<td style="text-align:left;">
DISS
</td>
</tr>
<tr>
<td style="text-align:left;">
Csrnp3
</td>
<td style="text-align:right;">
-1.46
</td>
<td style="text-align:right;">
6.55e-03
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Stxbp5l
</td>
<td style="text-align:right;">
-1.94
</td>
<td style="text-align:right;">
9.07e-03
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Gabrb1
</td>
<td style="text-align:right;">
-1.07
</td>
<td style="text-align:right;">
1.03e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Rgs7bp
</td>
<td style="text-align:right;">
-1.36
</td>
<td style="text-align:right;">
1.14e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Grin2b
</td>
<td style="text-align:right;">
-1.66
</td>
<td style="text-align:right;">
1.32e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Sorl1
</td>
<td style="text-align:right;">
-1.41
</td>
<td style="text-align:right;">
1.35e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Birc6
</td>
<td style="text-align:right;">
-1.32
</td>
<td style="text-align:right;">
1.60e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Cdkl5
</td>
<td style="text-align:right;">
-1.67
</td>
<td style="text-align:right;">
1.60e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Ralgapa1
</td>
<td style="text-align:right;">
-1.19
</td>
<td style="text-align:right;">
1.63e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Rc3h2
</td>
<td style="text-align:right;">
-1.97
</td>
<td style="text-align:right;">
1.87e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Epha6
</td>
<td style="text-align:right;">
-1.70
</td>
<td style="text-align:right;">
2.01e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Kcnma1
</td>
<td style="text-align:right;">
-1.92
</td>
<td style="text-align:right;">
2.06e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Nos1ap
</td>
<td style="text-align:right;">
-1.23
</td>
<td style="text-align:right;">
2.21e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Plppr4
</td>
<td style="text-align:right;">
-1.25
</td>
<td style="text-align:right;">
2.49e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Celsr2
</td>
<td style="text-align:right;">
-1.13
</td>
<td style="text-align:right;">
2.54e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Grin2a
</td>
<td style="text-align:right;">
-1.66
</td>
<td style="text-align:right;">
2.69e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Lrrc7
</td>
<td style="text-align:right;">
-1.78
</td>
<td style="text-align:right;">
2.73e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Cadm2
</td>
<td style="text-align:right;">
-1.21
</td>
<td style="text-align:right;">
3.00e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Magi2
</td>
<td style="text-align:right;">
-1.27
</td>
<td style="text-align:right;">
3.21e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Rasgrp1
</td>
<td style="text-align:right;">
-1.08
</td>
<td style="text-align:right;">
3.78e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Rock2
</td>
<td style="text-align:right;">
-1.05
</td>
<td style="text-align:right;">
3.88e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Asap1
</td>
<td style="text-align:right;">
-1.88
</td>
<td style="text-align:right;">
4.02e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Cacna2d1
</td>
<td style="text-align:right;">
-1.22
</td>
<td style="text-align:right;">
4.17e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Gm42878
</td>
<td style="text-align:right;">
-3.80
</td>
<td style="text-align:right;">
4.37e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Kcnj6
</td>
<td style="text-align:right;">
-2.39
</td>
<td style="text-align:right;">
4.37e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Adcy9
</td>
<td style="text-align:right;">
-1.42
</td>
<td style="text-align:right;">
4.74e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Fat3
</td>
<td style="text-align:right;">
-2.18
</td>
<td style="text-align:right;">
5.40e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Taok1
</td>
<td style="text-align:right;">
-1.31
</td>
<td style="text-align:right;">
5.40e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Adgrl3
</td>
<td style="text-align:right;">
-1.35
</td>
<td style="text-align:right;">
5.70e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Slc4a8
</td>
<td style="text-align:right;">
-1.04
</td>
<td style="text-align:right;">
6.10e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Atxn1
</td>
<td style="text-align:right;">
-1.52
</td>
<td style="text-align:right;">
6.74e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Arfgef3
</td>
<td style="text-align:right;">
-1.52
</td>
<td style="text-align:right;">
6.84e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Arhgap32
</td>
<td style="text-align:right;">
-1.28
</td>
<td style="text-align:right;">
6.85e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Pde4d
</td>
<td style="text-align:right;">
-1.75
</td>
<td style="text-align:right;">
6.95e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Taf1
</td>
<td style="text-align:right;">
-1.21
</td>
<td style="text-align:right;">
7.15e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Clock
</td>
<td style="text-align:right;">
-1.24
</td>
<td style="text-align:right;">
7.25e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Atp2b1
</td>
<td style="text-align:right;">
-1.21
</td>
<td style="text-align:right;">
7.40e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Prr36
</td>
<td style="text-align:right;">
-1.76
</td>
<td style="text-align:right;">
7.40e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Cacna1e
</td>
<td style="text-align:right;">
-1.28
</td>
<td style="text-align:right;">
8.05e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Gnaq
</td>
<td style="text-align:right;">
-1.18
</td>
<td style="text-align:right;">
8.28e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Syt14
</td>
<td style="text-align:right;">
-1.50
</td>
<td style="text-align:right;">
8.35e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Garem
</td>
<td style="text-align:right;">
-1.39
</td>
<td style="text-align:right;">
8.88e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Pde4c
</td>
<td style="text-align:right;">
-4.00
</td>
<td style="text-align:right;">
9.01e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Atp2c1
</td>
<td style="text-align:right;">
-1.13
</td>
<td style="text-align:right;">
9.60e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Cacnb2
</td>
<td style="text-align:right;">
-1.16
</td>
<td style="text-align:right;">
9.60e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Ubqln1
</td>
<td style="text-align:right;">
-1.06
</td>
<td style="text-align:right;">
9.74e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
<tr>
<td style="text-align:left;">
Atrx
</td>
<td style="text-align:right;">
-1.06
</td>
<td style="text-align:right;">
9.92e-02
</td>
<td style="text-align:left;">
HOMO
</td>
</tr>
</tbody>
</table>

**Preview of Supplemental Table 2. Molecules implicated in hippocampal
LTP from Sanes and Lichtman 1999.** This table list the molecules review
by Sanes and Lichtman in their 1999 review article and the related
transcripts that were investigated in this study. *Full table available
at
<a href="https://github.com/raynamharris/DissociationTest/blob/master/data/SanesLichtman.csv" class="uri">https://github.com/raynamharris/DissociationTest/blob/master/data/SanesLichtman.csv</a>.*

<table>
<thead>
<tr>
<th style="text-align:left;">
Sanes…Lichtman.Molecules
</th>
<th style="text-align:left;">
Related.Transcripts
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
GLUTAMATE RECEPTORS
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
GluR1; GluR2
</td>
<td style="text-align:left;">
Gria1; Gria2
</td>
</tr>
<tr>
<td style="text-align:left;">
mGluR1; mGluR4; mGluR5; mGluR7
</td>
<td style="text-align:left;">
Grm1; Grm4; Grm5; Grm7
</td>
</tr>
<tr>
<td style="text-align:left;">
NMDA NR2A; NMDA NR2D; NMDA NR1
</td>
<td style="text-align:left;">
Grin1; Grin2a; Grin2d
</td>
</tr>
<tr>
<td style="text-align:left;">
OTHER NEUROTRANSMITTERS
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
norepinephrine and b-adrenergic receptors
</td>
<td style="text-align:left;">
Adrb1; Adrb2; Adrb3
</td>
</tr>
</tbody>
</table>

**Preview of Supplemental Table 3. Marker genes for astrocytes,
oligodendrocytes, microglia, and neurons.** This table lists the genes
from Cahoy et al., 2008 that we investigated to estimate the relative
abundance of cell types in the examined tissue. LFC: Limit fold change.
*Full table available at
<a href="https://github.com/raynamharris/DissociationTest/blob/master/results/markergenes.csv" class="uri">https://github.com/raynamharris/DissociationTest/blob/master/results/markergenes.csv</a>.*

<table>
<thead>
<tr>
<th style="text-align:left;">
marker
</th>
<th style="text-align:left;">
gene
</th>
<th style="text-align:right;">
lfc
</th>
<th style="text-align:right;">
padj
</th>
<th style="text-align:left;">
direction
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
astrocyte
</td>
<td style="text-align:left;">
ALDOC
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.348
</td>
<td style="text-align:left;">
neither
</td>
</tr>
<tr>
<td style="text-align:left;">
astrocyte
</td>
<td style="text-align:left;">
AQP4
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.964
</td>
<td style="text-align:left;">
neither
</td>
</tr>
<tr>
<td style="text-align:left;">
astrocyte
</td>
<td style="text-align:left;">
FGFR3
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.951
</td>
<td style="text-align:left;">
neither
</td>
</tr>
<tr>
<td style="text-align:left;">
astrocyte
</td>
<td style="text-align:left;">
GFAP
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.615
</td>
<td style="text-align:left;">
neither
</td>
</tr>
<tr>
<td style="text-align:left;">
astrocyte
</td>
<td style="text-align:left;">
GJB6
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.799
</td>
<td style="text-align:left;">
neither
</td>
</tr>
<tr>
<td style="text-align:left;">
astrocyte
</td>
<td style="text-align:left;">
SLC1A2
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.987
</td>
<td style="text-align:left;">
neither
</td>
</tr>
</tbody>
</table>

Outdated Video Abstract
-----------------------

I made this [short video](https://www.youtube.com/watch?v=taeAqimxXWo)
explaining how to use this GitHub repo when I submitted the first draft
to the journal Hippocampus and posted a pre-print on BioRxiv.

[![screenshot](./figures/screenshot.png)](https://www.youtube.com/watch?v=taeAqimxXWo)
