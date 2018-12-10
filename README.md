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


**Background**: Single-neuron gene expression studies may be especially
important for understanding nervous system structure and function
because of the neuron-specific functionality and plasticity that defines
functional neural circuits. Cellular dissociation is a prerequisite
technical manipulation for single-cell and single cell-population
studies, but the extent to which the cellular dissociation process cells
affects neural gene expression has not been determined. This information
is necessary for interpreting the results of experimental manipulations
that affect neural function such as learning and memory. The goal of
this research was to determine the impact of chemical cell dissociation
on transcriptome. **Results**: This report shows that the process of
chemical cellular dissociation compared to homogenization alters about
2% of the tissue-level transcriptome of hippocampal subfields. Genes
related to cellular stress response pathways are activated by
dissociation compared to homogenization. Few genes canonically
implicated in learning and memory are found to change expression levels
in response to the dissociation procedure. **Discussion**: This study
suggests that chemical cellular dissociation has minimal but specific
affect genes or molecular functions canonically related to learning and
memory. However, sample preparation can affect gene expression profiles,
which might confound interpretation of results depending on the research
question. This study is important for the investigation of any complex
tissues as research effort moves from subfield level analysis to single
cell analysis of gene expression. **Methods**: We compared tissue level
expression of microdissected samples from the dentate gyrus (DG), CA3,
and CA1 subfields of the mouse hippocampus either prepared by a standard
tissue homogenization protocol or subjected to a cellular dissociation
procedure. We used the Illumina HiSeq platform for sequencing, Kallisto
for transcript abundance estimation, DESeq2 for differential gene
expression profiling, and GO\_MWU for analysis of gene ontology. Raw
reads, results, and code are available at GEO and on GitHub.

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

Results
-------

The following descriptions are not ready for publication, rather they
are converational descriptions of the current state.

![](./figures/figure1.png)

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

**Figure 2.** The

![](./figures/figure3.png)

**Figure 3.** The gene list and go terms that everyone wants to know. 2A
The top 30 most differentially expressed genes. Genes are clustered by
correlation but samples are NOT clustered! Only “Jun” jumps out as a
gene related to learning and memory. 2B List of molecular function
categories that are enriched or depleted in the dissociated tissues
relative to controls. Again, nothing jumps out as classic memory
pathways, but there are some interesting affects on DNA regulation
(methylation, chromatin, histone, RNA binding, helicase) and metabolism
(oxidoreductase, cytokin, growth factors, ligase).

**Supplemental Table 1.** Expression level and fold change of of
significant genes (p &lt; 0.05) between dissociated tissue and
homogenized tissue. This table shows the log fold change (lfc), p-value
(padj), and direction of upregulation for each gene analyzed.

For now at:
<https://github.com/raynamharris/DissociationTest/blob/master/results/dissociationDEGs.csv>

    suptable <- read.csv("~/Github/DissociationTest/results/dissociationDEGs.csv", header = T)
    suptable %>% dplyr::filter(direction == "DISS")

    ##          gene      lfc         padj direction
    ## 1         Trf 2.724763 5.314245e-07      DISS
    ## 2        Hexb 2.348231 8.098639e-07      DISS
    ## 3      Selplg 2.969442 9.216085e-07      DISS
    ## 4        C1qb 2.276248 7.066907e-06      DISS
    ## 5       Csf1r 2.133675 9.581268e-06      DISS
    ## 6        Ctss 2.587431 9.581268e-06      DISS
    ## 7         Cnp 2.452742 4.476356e-05      DISS
    ## 8        Il1a 3.060196 4.476356e-05      DISS
    ## 9         Mag 3.305352 4.476356e-05      DISS
    ## 10       Cd14 3.379471 4.878753e-05      DISS
    ## 11      Mpeg1 2.416742 6.466041e-05      DISS
    ## 12    Tmem88b 3.141197 6.775210e-05      DISS
    ## 13     Nfkbia 2.104497 6.838759e-05      DISS
    ## 14    Slc15a3 3.447553 9.319157e-05      DISS
    ## 15       Cd83 2.583078 9.986852e-05      DISS
    ## 16     Cldn11 3.139729 1.048558e-04      DISS
    ## 17     Laptm5 2.307328 2.239322e-04      DISS
    ## 18       Plau 3.896819 2.239322e-04      DISS
    ## 19       Plp1 2.708819 2.239322e-04      DISS
    ## 20        Mal 3.197771 2.318480e-04      DISS
    ## 21       Plek 2.503064 2.318480e-04      DISS
    ## 22        Jun 1.596999 3.198344e-04      DISS
    ## 23       Mobp 2.595810 4.414573e-04      DISS
    ## 24      Socs3 2.263144 4.796366e-04      DISS
    ## 25       Pllp 2.762201 5.457963e-04      DISS
    ## 26      Ndrg1 2.874903 5.795321e-04      DISS
    ## 27     Phldb1 2.268954 6.850375e-04      DISS
    ## 28     Tyrobp 3.339132 8.216980e-04      DISS
    ## 29     Cx3cr1 2.095318 9.504111e-04      DISS
    ## 30    Tmem119 2.543998 9.504111e-04      DISS
    ## 31       C1qa 2.324221 9.845767e-04      DISS
    ## 32     Slc2a5 4.077030 1.225542e-03      DISS
    ## 33      Fcrls 2.670275 1.397224e-03      DISS
    ## 34       Lgmn 1.839669 1.570913e-03      DISS
    ## 35      Clic4 2.008972 1.808469e-03      DISS
    ## 36     Sema5a 3.251563 2.975420e-03      DISS
    ## 37      C5ar1 3.656770 3.007771e-03      DISS
    ## 38      Cmtm5 2.664280 3.007771e-03      DISS
    ## 39      Ctcfl 3.079221 3.046426e-03      DISS
    ## 40      Cryab 2.805439 3.068789e-03      DISS
    ## 41       Gab1 2.897868 3.443607e-03      DISS
    ## 42       C1qc 2.599005 4.398792e-03      DISS
    ## 43       Fa2h 3.301692 4.606123e-03      DISS
    ## 44     Clec2l 1.940032 5.012805e-03      DISS
    ## 45      Spry2 1.794098 5.077573e-03      DISS
    ## 46       Mcam 2.858219 5.123514e-03      DISS
    ## 47     Olfml3 2.087093 5.184144e-03      DISS
    ## 48      Sepp1 1.537069 5.780364e-03      DISS
    ## 49       Lag3 3.323814 8.025488e-03      DISS
    ## 50        Mbp 1.950857 8.025488e-03      DISS
    ## 51       Atf3 2.173984 8.091342e-03      DISS
    ## 52      C1ql1 3.097264 8.091342e-03      DISS
    ## 53     Adgrg1 1.860419 8.617822e-03      DISS
    ## 54      Csrp1 1.697590 9.065359e-03      DISS
    ## 55      Ltbp3 2.291751 9.065359e-03      DISS
    ## 56     Rpl36a 1.569502 9.065359e-03      DISS
    ## 57   Ighv14-2 7.813700 9.599492e-03      DISS
    ## 58      Mapk3 1.606075 9.799875e-03      DISS
    ## 59      Trem2 3.251647 1.141690e-02      DISS
    ## 60     Syngr2 3.353063 1.196476e-02      DISS
    ## 61       Irf8 2.772496 1.264753e-02      DISS
    ## 62       Cst3 1.529705 1.452579e-02      DISS
    ## 63    Unc93b1 2.280418 1.551209e-02      DISS
    ## 64  D6Wsu163e 2.320522 1.591633e-02      DISS
    ## 65      Sstr1 4.261714 1.591633e-02      DISS
    ## 66    Gm10401 7.622706 1.602914e-02      DISS
    ## 67      Itga5 3.054466 1.602914e-02      DISS
    ## 68       Ly86 2.805981 1.602914e-02      DISS
    ## 69      Olig1 1.672758 1.602914e-02      DISS
    ## 70      Zfp36 1.788852 1.602914e-02      DISS
    ## 71        Cd9 2.358946 1.715204e-02      DISS
    ## 72      Rph3a 2.264486 1.715204e-02      DISS
    ## 73     Slc2a1 1.698929 1.769124e-02      DISS
    ## 74       Bin2 2.492379 1.827866e-02      DISS
    ## 75      Gpr84 2.858254 1.868563e-02      DISS
    ## 76      H2-K1 2.091260 1.868563e-02      DISS
    ## 77    Slco2b1 2.171535 1.868563e-02      DISS
    ## 78     Tgfbr2 2.440158 1.868563e-02      DISS
    ## 79       Cdh9 3.653615 1.992962e-02      DISS
    ## 80      Ptgds 2.528494 1.992962e-02      DISS
    ## 81      Itgb5 1.978733 2.014215e-02      DISS
    ## 82     Inpp5d 2.844983 2.038075e-02      DISS
    ## 83      Icam1 2.866966 2.064192e-02      DISS
    ## 84       Myrf 2.208911 2.210619e-02      DISS
    ## 85      Sparc 1.789037 2.210619e-02      DISS
    ## 86        Tnf 2.401661 2.210619e-02      DISS
    ## 87    Siglech 2.112824 2.214786e-02      DISS
    ## 88  Serpinb1a 4.357837 2.259866e-02      DISS
    ## 89        Mog 2.483406 2.271517e-02      DISS
    ## 90      Stat3 1.558550 2.408164e-02      DISS
    ## 91       Tlr7 3.888201 2.423416e-02      DISS
    ## 92     Lgals9 2.744048 2.622194e-02      DISS
    ## 93       Cyba 3.804768 2.632412e-02      DISS
    ## 94       Ermn 2.529048 2.690668e-02      DISS
    ## 95        Fn1 1.773151 2.690668e-02      DISS
    ## 96      Sipa1 2.713651 2.690668e-02      DISS
    ## 97       Vasp 2.310592 2.731750e-02      DISS
    ## 98    Shroom3 6.227512 2.796124e-02      DISS
    ## 99       Qdpr 1.549092 2.856534e-02      DISS
    ## 100   Slc12a2 2.071560 2.856534e-02      DISS
    ## 101   Tmem170 2.526577 2.934837e-02      DISS
    ## 102      Apod 2.047676 2.996745e-02      DISS
    ## 103      Ccl4 1.940311 2.996745e-02      DISS
    ## 104      Il1b 2.405914 2.996745e-02      DISS
    ## 105   Plekhg3 2.646614 2.996745e-02      DISS
    ## 106     Sept4 1.841197 2.996745e-02      DISS
    ## 107      Cd82 2.176665 3.203070e-02      DISS
    ## 108      Ucp2 2.489523 3.211818e-02      DISS
    ## 109     Smoc2 3.813698 3.266403e-02      DISS
    ## 110    mt-Nd3 1.696531 3.284794e-02      DISS
    ## 111      Fosb 1.585131 3.322296e-02      DISS
    ## 112      Mfng 4.045092 3.322296e-02      DISS
    ## 113      Emp2 2.347526 3.326797e-02      DISS
    ## 114      Cd37 4.135444 3.356350e-02      DISS
    ## 115     Scrg1 2.632423 3.362097e-02      DISS
    ## 116     Prr18 2.292415 3.567104e-02      DISS
    ## 117     Cd274 4.237913 4.119354e-02      DISS
    ## 118     Synpr 2.305883 4.119354e-02      DISS
    ## 119     Pold1 2.225724 4.165158e-02      DISS
    ## 120     Crtc2 1.779653 4.371455e-02      DISS
    ## 121       Msn 1.747697 4.545368e-02      DISS
    ## 122       Myc 2.346465 4.563485e-02      DISS
    ## 123     Bcas1 1.697616 4.698784e-02      DISS
    ## 124      Gatm 2.121835 4.742695e-02      DISS
    ## 125      Gjb1 4.758855 4.742695e-02      DISS
    ## 126     Pros1 2.289380 4.742695e-02      DISS
    ## 127   Rab3il1 2.755303 4.777123e-02      DISS
    ## 128     Arl4c 1.794929 4.892278e-02      DISS
    ## 129     Dhrs3 2.947787 4.892278e-02      DISS
    ## 130     Enpp2 1.768966 4.892278e-02      DISS
    ## 131     Gpr34 2.079302 4.892278e-02      DISS
    ## 132     Ntng1 3.081533 4.892278e-02      DISS
    ## 133    P2ry12 1.741791 4.892278e-02      DISS
    ## 134      Tlr2 2.646241 4.892278e-02      DISS
    ## 135   Cyp27a1 3.582689 4.918113e-02      DISS
    ## 136     Asap3 3.892399 4.961699e-02      DISS
    ## 137     Itgam 1.746838 4.961699e-02      DISS
    ## 138      Lcp1 2.733823 4.961699e-02      DISS

    suptable %>% dplyr::filter(direction == "HOMO")

    ##       gene       lfc        padj direction
    ## 1  Stxbp5l -1.940397 0.009065359      HOMO
    ## 2   Grin2b -1.656666 0.013241261      HOMO
    ## 3    Cdkl5 -1.665223 0.016029141      HOMO
    ## 4    Rc3h2 -1.965015 0.018685627      HOMO
    ## 5    Epha6 -1.702348 0.020058278      HOMO
    ## 6   Kcnma1 -1.918591 0.020641917      HOMO
    ## 7   Grin2a -1.659562 0.026906675      HOMO
    ## 8    Lrrc7 -1.780252 0.027317499      HOMO
    ## 9    Asap1 -1.880316 0.040158237      HOMO
    ## 10 Gm42878 -3.795602 0.043652913      HOMO
    ## 11   Kcnj6 -2.392068 0.043652913      HOMO

**Supplemental Table 2.** Gene ontologies of enriched genes. The first
row contains the GO category (either MF or CC). The second is the GO
term. Also shown are directionally, unumber of enriched genes in that
catory out of the total (ratio), and p-value.

<https://github.com/raynamharris/DissociationTest/blob/master/results/GOsignificantcatagories.csv>

    GOtable <- read.csv("~/Github/DissociationTest/results/GOsignificantcatagories.csv", header = T)
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
