[![Binder](http://mybinder.org/badge.svg)](http://beta.mybinder.org/v2/gh/raynamharris/DissociationTest/master?urlpath=rstudio)
*Click the button to launch a Binder R session. Navigate to the
`scripts` directory and open any `.Rmd` file. Note: the first two are
long and slow. The rest are quick scripts that make figures.*

Hippocampal transcriptomic responses to cellular dissociation
=============================================================

Rayna M. Harris, Hsin-Yi Kao, Juan Marcos Alarcon, Hans A. Hofmann, and
André A. Fenton

Overview
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
about 350 genes or 2% of the hippocampal transcriptome. Few genes
canonically implicated in long-term potentiation (LTP) and fear memory
change expression levels in response to the dissociation procedure.
Nevertheless, sample preparation did affect gene expression profiles,
which might confound interpretation of results depending on the research
question. This study is important for the investigation of any complex
tissues as research effort moves from subfield level analysis to single
cell analysis of gene expression.

Repo Contents
-------------

-   [**data**](./data/): contains most of the input processed data
    files. Large data fiels are stored in the Gene Expression Omnibus at
    [GSE99765](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99765)
    . Raw kallisto abundance files are also stored in my other GitHub
    repo called
    [MouseHippocampusRNAseqData](https://github.com/raynamharris/MouseHippocampusRNAseqData).
-   [**UNIXworkflow**](./UNIXworkflow/): This descirbes the process I
    used to process my files using the Stampede Cluster at the Texas
    Advanced computing facility
-   [**scripts**](./scripts/): this contains all the .R and .Rmd scripts
    as well as the .md output files. They have prefixes to convey the
    order of operation.
-   [**figures**](./figures/): Contains all output for all files from
    the Rmarkdown scripts and my adobe-created images.

General workflow and approach
-----------------------------

1.  experimental design (treatment \* hippocampal subfield)
2.  RNA-seq (Illumina, GSAF)
3.  bioinformatics (TACC, FASTQC, cutadapt, kallisto)
4.  data viz and stats (DESeq2, GOMWU, R)
5.  version control and sharing (Git, GitHub, NCBI)

Figures and Tables
------------------

![](./figures/figure1.png)

**Figure 1. Experimental design and global expression gene expression
patterns.** A) Experimental design. Two tissue samples were taken from
three hippocampal subfields (CA1, CA3, and DG) from 300 um brain slices.
Two adjacent samples were processed using a homogenization (HOMO)
protocol or dissociated (DISS) before processing for tissue level gene
expression profiling. B) Dissociation does not yield subfield-specific
changes in gene expression between homogenized (HOMO, open circles,
dashed ellipse) and dissociated tissues (DISS, filled circles, solid
ellipse). PC1 accounts for 40% of all gene expression variation and by
inspection, separates the DG samples (orange circles) from the CA1
(purple circles) and CA3 samples (green circles). PC2 accounts for 22%
of the variation in gene expression and varies significantly with
treatment. The ellipses estimate the 95% confidence interval for a
multivariate t-distribution for homogenized (dashed line) and
dissociated (solid line) samples.

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

**Table 1.** Number of differentially expressed genes by 2-way contrast.
“Down-regulated” means that expression is higher in the term on the far
left of column 1. “Up-regulated” means that expression is higher in the
term on the right of column 1. This shows that there is more variation
due to subfield than treatment.

![](./figures/figure2.png)

**Figure 2. Enzymatic dissociation has a moderate effect on hippocampal
gene expression patterns compared to homogenized tissue.** A) Volcano
plot showing gene expression fold-difference and significance between
treatment groups. We found that 56 genes are up-regulated in the
homogenization control group (open circles) while 288 genes are
up-regulated in the dissociated treatment group (filled dark grey
circles). Genes below the p-value &lt; 0.1 (or –log p-value &lt; 1) are
shown in light grey. B) Heatmap showing the top 30 differentially
expressed genes between dissociated and homogenized tissue. Square boxes
at the top are color coded by sample (white: homogenized, grey:
dissociated, purple: CA1, green: CA3, orange: DG. Within the heatmap,
log fold difference levels of expression are indicated by the
blue-green-yellow gradient with lighter colors indicating increased
expression.

![](./figures/figure3.png)

**Preview of Supplemental Table 1.** Expression level and fold change of
significant genes (p &lt; 0.1) between dissociated tissue and
homogenized tissue. This table shows the log fold change (lfc), p-value
(padj), and direction of upregulation for each gene analyzed.

Full table available at
<https://github.com/raynamharris/DissociationTest/blob/master/results/dissociationDEGs.csv>.

    ##     gene      lfc         padj direction
    ## 1    Trf 2.724763 5.314245e-07      DISS
    ## 2   Hexb 2.348231 8.098639e-07      DISS
    ## 3 Selplg 2.969442 9.216085e-07      DISS
    ## 4   C1qb 2.276248 7.066907e-06      DISS
    ## 5  Csf1r 2.133675 9.581268e-06      DISS
    ## 6   Ctss 2.587431 9.581268e-06      DISS

    ##      gene       lfc        padj direction
    ## 1  Csrnp3 -1.457990 0.006552580      HOMO
    ## 2 Stxbp5l -1.940397 0.009065359      HOMO
    ## 3  Gabrb1 -1.074515 0.010299888      HOMO
    ## 4  Rgs7bp -1.356723 0.011416902      HOMO
    ## 5  Grin2b -1.656666 0.013241261      HOMO
    ## 6   Sorl1 -1.408413 0.013455971      HOMO

**Preview of Supplemental Table 2. Molecules implicated in hippocampal
LTP from Sanes and Lichtman 1999.**

Full table available at
<https://github.com/raynamharris/DissociationTest/blob/master/data/SanesLichtman.csv>.

    ##                  Sanes...Lichtman.molecules      related.transcripts
    ## 1                       GLUTAMATE RECEPTORS                         
    ## 2                              GluR1; GluR2            Gria1; Gria2 
    ## 3            mGluR1; mGluR4; mGluR5; mGluR7  Grm1; Grm4; Grm5; Grm7 
    ## 4            NMDA NR2A; NMDA NR2D; NMDA NR1   Grin1; Grin2a; Grin2d 
    ## 5                   OTHER NEUROTRANSMITTERS                         
    ## 6 norepinephrine and b-adrenergic receptors     Adrb1; Adrb2; Adrb3

Outdated Video Abstract
-----------------------

I made this [short video](https://www.youtube.com/watch?v=taeAqimxXWo)
explaining how to use this GitHub repo when I submitted the first draft
to the journal Hippocampus and posted a pre-print on BioRxiv.

[![screenshot](./figures/screenshot.png)](https://www.youtube.com/watch?v=taeAqimxXWo)
