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

**Supplemental Table 1.** Expression level and fold change of
significant genes (p &lt; 0.1) between dissociated tissue and
homogenized tissue. This table shows the log fold change (lfc), p-value
(padj), and direction of upregulation for each gene analyzed.

For now at:
<https://github.com/raynamharris/DissociationTest/blob/master/results/dissociationDEGs.csv>

    suptable <- read.csv("~/Github/DissociationTest/results/dissociationDEGs.csv", header = T)
    suptable %>% dplyr::filter(direction == "DISS")

    ##              gene      lfc         padj direction
    ## 1             Trf 2.724763 5.314245e-07      DISS
    ## 2            Hexb 2.348231 8.098639e-07      DISS
    ## 3          Selplg 2.969442 9.216085e-07      DISS
    ## 4            C1qb 2.276248 7.066907e-06      DISS
    ## 5           Csf1r 2.133675 9.581268e-06      DISS
    ## 6            Ctss 2.587431 9.581268e-06      DISS
    ## 7             Cnp 2.452742 4.476356e-05      DISS
    ## 8            Il1a 3.060196 4.476356e-05      DISS
    ## 9             Mag 3.305352 4.476356e-05      DISS
    ## 10           Cd14 3.379471 4.878753e-05      DISS
    ## 11          Mpeg1 2.416742 6.466041e-05      DISS
    ## 12        Tmem88b 3.141197 6.775210e-05      DISS
    ## 13         Nfkbia 2.104497 6.838759e-05      DISS
    ## 14        Slc15a3 3.447553 9.319157e-05      DISS
    ## 15           Cd83 2.583078 9.986852e-05      DISS
    ## 16         Cldn11 3.139729 1.048558e-04      DISS
    ## 17         Laptm5 2.307328 2.239322e-04      DISS
    ## 18           Plau 3.896819 2.239322e-04      DISS
    ## 19           Plp1 2.708819 2.239322e-04      DISS
    ## 20            Mal 3.197771 2.318480e-04      DISS
    ## 21           Plek 2.503064 2.318480e-04      DISS
    ## 22            Jun 1.596999 3.198344e-04      DISS
    ## 23           Mobp 2.595810 4.414573e-04      DISS
    ## 24         Pcsk1n 1.298848 4.796366e-04      DISS
    ## 25          Socs3 2.263144 4.796366e-04      DISS
    ## 26           Pllp 2.762201 5.457963e-04      DISS
    ## 27          Ndrg1 2.874903 5.795321e-04      DISS
    ## 28         Phldb1 2.268954 6.850375e-04      DISS
    ## 29         Tyrobp 3.339132 8.216980e-04      DISS
    ## 30         Cx3cr1 2.095318 9.504111e-04      DISS
    ## 31        Tmem119 2.543998 9.504111e-04      DISS
    ## 32           C1qa 2.324221 9.845767e-04      DISS
    ## 33         Slc2a5 4.077030 1.225542e-03      DISS
    ## 34          Fcrls 2.670275 1.397224e-03      DISS
    ## 35           Lgmn 1.839669 1.570913e-03      DISS
    ## 36          Clic4 2.008972 1.808469e-03      DISS
    ## 37         Sema5a 3.251563 2.975420e-03      DISS
    ## 38          C5ar1 3.656770 3.007771e-03      DISS
    ## 39          Cmtm5 2.664280 3.007771e-03      DISS
    ## 40          Ctcfl 3.079221 3.046426e-03      DISS
    ## 41          Cryab 2.805439 3.068789e-03      DISS
    ## 42          Rpl26 1.270111 3.217652e-03      DISS
    ## 43           Gab1 2.897868 3.443607e-03      DISS
    ## 44           C1qc 2.599005 4.398792e-03      DISS
    ## 45           Fa2h 3.301692 4.606123e-03      DISS
    ## 46         Clec2l 1.940032 5.012805e-03      DISS
    ## 47          Spry2 1.794098 5.077573e-03      DISS
    ## 48           Mcam 2.858219 5.123514e-03      DISS
    ## 49         Olfml3 2.087093 5.184144e-03      DISS
    ## 50          Sepp1 1.537069 5.780364e-03      DISS
    ## 51           Ptma 1.189007 6.207166e-03      DISS
    ## 52          Dusp1 1.345553 8.025488e-03      DISS
    ## 53           Lag3 3.323814 8.025488e-03      DISS
    ## 54            Mbp 1.950857 8.025488e-03      DISS
    ## 55           Atf3 2.173984 8.091342e-03      DISS
    ## 56          C1ql1 3.097264 8.091342e-03      DISS
    ## 57         Adgrg1 1.860419 8.617822e-03      DISS
    ## 58          Csrp1 1.697590 9.065359e-03      DISS
    ## 59          Ltbp3 2.291751 9.065359e-03      DISS
    ## 60           Mcl1 1.406205 9.065359e-03      DISS
    ## 61          Rpl29 1.162075 9.065359e-03      DISS
    ## 62         Rpl36a 1.569502 9.065359e-03      DISS
    ## 63          Rps27 1.293098 9.589995e-03      DISS
    ## 64       Ighv14-2 7.813700 9.599492e-03      DISS
    ## 65          Mapk3 1.606075 9.799875e-03      DISS
    ## 66          Trem2 3.251647 1.141690e-02      DISS
    ## 67         Syngr2 3.353063 1.196476e-02      DISS
    ## 68           Irf8 2.772496 1.264753e-02      DISS
    ## 69           Jund 1.063997 1.264753e-02      DISS
    ## 70        Plekhb1 1.488307 1.264753e-02      DISS
    ## 71           Cst3 1.529705 1.452579e-02      DISS
    ## 72          Rpl28 1.169842 1.551209e-02      DISS
    ## 73        Unc93b1 2.280418 1.551209e-02      DISS
    ## 74      D6Wsu163e 2.320522 1.591633e-02      DISS
    ## 75          Sstr1 4.261714 1.591633e-02      DISS
    ## 76        Gm10401 7.622706 1.602914e-02      DISS
    ## 77          Itga5 3.054466 1.602914e-02      DISS
    ## 78           Ly86 2.805981 1.602914e-02      DISS
    ## 79          Olig1 1.672758 1.602914e-02      DISS
    ## 80          Zfp36 1.788852 1.602914e-02      DISS
    ## 81           Rps8 1.038424 1.697244e-02      DISS
    ## 82            Cd9 2.358946 1.715204e-02      DISS
    ## 83          Rph3a 2.264486 1.715204e-02      DISS
    ## 84         Slc2a1 1.698929 1.769124e-02      DISS
    ## 85           Bin2 2.492379 1.827866e-02      DISS
    ## 86          Gpr84 2.858254 1.868563e-02      DISS
    ## 87          H2-K1 2.091260 1.868563e-02      DISS
    ## 88        Slco2b1 2.171535 1.868563e-02      DISS
    ## 89         Tgfbr2 2.440158 1.868563e-02      DISS
    ## 90           Cdh9 3.653615 1.992962e-02      DISS
    ## 91           Ctsd 1.257429 1.992962e-02      DISS
    ## 92          Ptgds 2.528494 1.992962e-02      DISS
    ## 93            Mt1 1.417219 2.005828e-02      DISS
    ## 94          Itgb5 1.978733 2.014215e-02      DISS
    ## 95         Inpp5d 2.844983 2.038075e-02      DISS
    ## 96          Icam1 2.866966 2.064192e-02      DISS
    ## 97           Myrf 2.208911 2.210619e-02      DISS
    ## 98          Sparc 1.789037 2.210619e-02      DISS
    ## 99            Tnf 2.401661 2.210619e-02      DISS
    ## 100       Siglech 2.112824 2.214786e-02      DISS
    ## 101     Serpinb1a 4.357837 2.259866e-02      DISS
    ## 102        mt-Co3 1.401657 2.262518e-02      DISS
    ## 103           Mog 2.483406 2.271517e-02      DISS
    ## 104         Stat3 1.558550 2.408164e-02      DISS
    ## 105          Tlr7 3.888201 2.423416e-02      DISS
    ## 106        Lgals9 2.744048 2.622194e-02      DISS
    ## 107          Tln1 1.336453 2.622194e-02      DISS
    ## 108          Cyba 3.804768 2.632412e-02      DISS
    ## 109          Ermn 2.529048 2.690668e-02      DISS
    ## 110           Fn1 1.773151 2.690668e-02      DISS
    ## 111        mt-Co1 1.240495 2.690668e-02      DISS
    ## 112         Sipa1 2.713651 2.690668e-02      DISS
    ## 113          Vasp 2.310592 2.731750e-02      DISS
    ## 114       Shroom3 6.227512 2.796124e-02      DISS
    ## 115       mt-Cytb 1.270366 2.856534e-02      DISS
    ## 116          Qdpr 1.549092 2.856534e-02      DISS
    ## 117       Slc12a2 2.071560 2.856534e-02      DISS
    ## 118       Tmem170 2.526577 2.934837e-02      DISS
    ## 119          Apod 2.047676 2.996745e-02      DISS
    ## 120          Ccl4 1.940311 2.996745e-02      DISS
    ## 121         Erdr1 1.327568 2.996745e-02      DISS
    ## 122          Il1b 2.405914 2.996745e-02      DISS
    ## 123       Plekhg3 2.646614 2.996745e-02      DISS
    ## 124         Sept4 1.841197 2.996745e-02      DISS
    ## 125          Cd82 2.176665 3.203070e-02      DISS
    ## 126         Cspg5 1.276950 3.211818e-02      DISS
    ## 127          Ucp2 2.489523 3.211818e-02      DISS
    ## 128         Smoc2 3.813698 3.266403e-02      DISS
    ## 129        mt-Nd3 1.696531 3.284794e-02      DISS
    ## 130          Fosb 1.585131 3.322296e-02      DISS
    ## 131          Mfng 4.045092 3.322296e-02      DISS
    ## 132       Pcdhga8 1.424150 3.322296e-02      DISS
    ## 133          Emp2 2.347526 3.326797e-02      DISS
    ## 134          Cd37 4.135444 3.356350e-02      DISS
    ## 135         Scrg1 2.632423 3.362097e-02      DISS
    ## 136         Prr18 2.292415 3.567104e-02      DISS
    ## 137        mt-Nd1 1.208515 3.600788e-02      DISS
    ## 138         Cd274 4.237913 4.119354e-02      DISS
    ## 139         Rpl32 1.334485 4.119354e-02      DISS
    ## 140         Synpr 2.305883 4.119354e-02      DISS
    ## 141         Pold1 2.225724 4.165158e-02      DISS
    ## 142          Btg2 1.393560 4.214353e-02      DISS
    ## 143         Crtc2 1.779653 4.371455e-02      DISS
    ## 144           Msn 1.747697 4.545368e-02      DISS
    ## 145        mt-Nd4 1.150734 4.545368e-02      DISS
    ## 146          Rtn2 1.473224 4.545368e-02      DISS
    ## 147           Myc 2.346465 4.563485e-02      DISS
    ## 148         Bcas1 1.697616 4.698784e-02      DISS
    ## 149          Gatm 2.121835 4.742695e-02      DISS
    ## 150          Gjb1 4.758855 4.742695e-02      DISS
    ## 151         Pros1 2.289380 4.742695e-02      DISS
    ## 152         Psat1 1.312726 4.742695e-02      DISS
    ## 153       Rab3il1 2.755303 4.777123e-02      DISS
    ## 154         Arl4c 1.794929 4.892278e-02      DISS
    ## 155         Dhrs3 2.947787 4.892278e-02      DISS
    ## 156         Enpp2 1.768966 4.892278e-02      DISS
    ## 157         Gpr34 2.079302 4.892278e-02      DISS
    ## 158         Ntng1 3.081533 4.892278e-02      DISS
    ## 159        P2ry12 1.741791 4.892278e-02      DISS
    ## 160          Tlr2 2.646241 4.892278e-02      DISS
    ## 161       Cyp27a1 3.582689 4.918113e-02      DISS
    ## 162         Gna12 1.298098 4.918946e-02      DISS
    ## 163         Asap3 3.892399 4.961699e-02      DISS
    ## 164         Itgam 1.746838 4.961699e-02      DISS
    ## 165          Lcp1 2.733823 4.961699e-02      DISS
    ## 166         Frmd8 2.793162 5.033222e-02      DISS
    ## 167          Ncf1 1.797278 5.066328e-02      DISS
    ## 168        Dnajb2 1.520992 5.142391e-02      DISS
    ## 169       mt-Atp6 1.246339 5.145182e-02      DISS
    ## 170        mt-Co2 1.363810 5.152893e-02      DISS
    ## 171         Wasf2 1.670352 5.198943e-02      DISS
    ## 172         Ccdc3 5.378969 5.323581e-02      DISS
    ## 173          Smox 1.583444 5.400312e-02      DISS
    ## 174       Tmem63a 2.456278 5.400312e-02      DISS
    ## 175          Vsir 2.055020 5.400312e-02      DISS
    ## 176         Gpr17 2.599991 5.455826e-02      DISS
    ## 177         Cmtm6 1.701083 5.753483e-02      DISS
    ## 178         Dock8 3.069416 5.753483e-02      DISS
    ## 179         Gp1bb 3.103505 5.753483e-02      DISS
    ## 180          Ier2 1.325464 5.753483e-02      DISS
    ## 181          Il16 3.318057 5.753483e-02      DISS
    ## 182         Rftn1 3.552850 5.753483e-02      DISS
    ## 183         Sox10 2.157168 5.753483e-02      DISS
    ## 184          Junb 1.031711 5.812471e-02      DISS
    ## 185 D630003M21Rik 3.123476 5.892826e-02      DISS
    ## 186         Kcnk6 2.839044 6.099169e-02      DISS
    ## 187          Rela 1.518490 6.180384e-02      DISS
    ## 188       mt-Nd4l 1.063672 6.210145e-02      DISS
    ## 189         Itgb4 2.929515 6.295265e-02      DISS
    ## 190        Zfp521 2.722945 6.295265e-02      DISS
    ## 191         Ilkap 1.612445 6.310124e-02      DISS
    ## 192       Nckap1l 2.107164 6.310124e-02      DISS
    ## 193         Ccrl2 3.008746 6.484511e-02      DISS
    ## 194          Gjc3 1.853694 6.484511e-02      DISS
    ## 195        Csrnp1 1.522313 6.624887e-02      DISS
    ## 196          Klf4 1.744125 6.624887e-02      DISS
    ## 197       Tbc1d16 1.476871 6.654471e-02      DISS
    ## 198         Wscd1 1.502189 6.686567e-02      DISS
    ## 199          Atf5 1.576562 6.854666e-02      DISS
    ## 200        Dusp10 2.488769 6.854666e-02      DISS
    ## 201         Prdm5 3.586916 6.854666e-02      DISS
    ## 202        Arpc1b 2.049508 6.888467e-02      DISS
    ## 203        mt-Nd2 1.188262 6.945136e-02      DISS
    ## 204          Grap 2.957468 7.075224e-02      DISS
    ## 205         Rpl18 1.035095 7.075224e-02      DISS
    ## 206      Suv420h2 2.246092 7.077014e-02      DISS
    ## 207         Il6ra 2.213044 7.210885e-02      DISS
    ## 208          Npr1 3.343268 7.244762e-02      DISS
    ## 209        Smtnl2 5.831492 7.253338e-02      DISS
    ## 210       mt-Atp8 1.105554 7.367276e-02      DISS
    ## 211        Piezo1 3.057865 7.367276e-02      DISS
    ## 212        Rpl23a 1.002800 7.367276e-02      DISS
    ## 213         Rpl31 1.085272 7.367276e-02      DISS
    ## 214          Pim1 2.413148 7.373646e-02      DISS
    ## 215          Cd63 1.580238 7.396167e-02      DISS
    ## 216          Csf1 1.881976 7.396167e-02      DISS
    ## 217         Rplp1 1.127254 7.396167e-02      DISS
    ## 218         Myo1f 4.776742 7.470057e-02      DISS
    ## 219        Insig1 1.611643 7.501089e-02      DISS
    ## 220         Itpk1 1.531741 7.577718e-02      DISS
    ## 221        Lhfpl2 2.151208 7.688014e-02      DISS
    ## 222         Chadl 1.823962 7.692316e-02      DISS
    ## 223         Nlrp3 2.196995 7.692316e-02      DISS
    ## 224        Sh3gl1 1.563999 7.692316e-02      DISS
    ## 225         Cidea 5.455145 7.846243e-02      DISS
    ## 226         Rack1 1.028424 7.846243e-02      DISS
    ## 227        Sema4g 2.300745 7.952104e-02      DISS
    ## 228         Rpl23 1.225101 8.018705e-02      DISS
    ## 229        Adam17 2.269214 8.050104e-02      DISS
    ## 230           Cpm 2.854974 8.150476e-02      DISS
    ## 231        Cox4i1 1.043149 8.273573e-02      DISS
    ## 232       Nectin1 1.841711 8.273573e-02      DISS
    ## 233         Ostf1 2.363306 8.273573e-02      DISS
    ## 234         Maml2 2.194435 8.281319e-02      DISS
    ## 235         Ap1s1 1.065627 8.377556e-02      DISS
    ## 236       Fam124a 2.361485 8.377556e-02      DISS
    ## 237            Hr 2.451564 8.377556e-02      DISS
    ## 238        Nfatc1 3.556189 8.377556e-02      DISS
    ## 239       Rps6kb2 1.652014 8.377556e-02      DISS
    ## 240        Stk32b 5.430373 8.388954e-02      DISS
    ## 241         Unc5b 1.955282 8.518824e-02      DISS
    ## 242        Cdkn1a 2.505978 8.541630e-02      DISS
    ## 243       Fam102b 1.525074 8.753250e-02      DISS
    ## 244        Sh3d19 2.429174 9.072159e-02      DISS
    ## 245        Tspan2 1.927555 9.072159e-02      DISS
    ## 246       Tnfaip6 2.763963 9.087126e-02      DISS
    ## 247          Cd68 2.346262 9.111811e-02      DISS
    ## 248         Mmp15 1.939710 9.192635e-02      DISS
    ## 249         Lamb2 2.699247 9.211781e-02      DISS
    ## 250         Itpkb 1.529801 9.248846e-02      DISS
    ## 251        Il10ra 2.490959 9.287823e-02      DISS
    ## 252        Rps27a 1.128669 9.333105e-02      DISS
    ## 253        Rpl35a 1.468246 9.445987e-02      DISS
    ## 254          Tbx1 5.677971 9.557721e-02      DISS
    ## 255          Srgn 2.405334 9.585503e-02      DISS
    ## 256          Apoe 1.252961 9.595410e-02      DISS
    ## 257          Bche 3.349739 9.595410e-02      DISS
    ## 258          Ccl9 3.308189 9.595410e-02      DISS
    ## 259          Cd81 1.027049 9.595410e-02      DISS
    ## 260           Cfh 2.026031 9.595410e-02      DISS
    ## 261          Ctsc 2.399444 9.595410e-02      DISS
    ## 262        Dusp26 2.023939 9.595410e-02      DISS
    ## 263         Erbin 1.719123 9.595410e-02      DISS
    ## 264          Gjc2 2.388262 9.595410e-02      DISS
    ## 265         Lrtm2 1.800672 9.595410e-02      DISS
    ## 266         Plin4 4.051876 9.595410e-02      DISS
    ## 267         Rgs10 2.175526 9.595410e-02      DISS
    ## 268      Rnaset2b 1.999401 9.595410e-02      DISS
    ## 269          Rpsa 1.076681 9.595410e-02      DISS
    ## 270      Slc25a10 2.942763 9.595410e-02      DISS
    ## 271         Spry1 2.184685 9.595410e-02      DISS
    ## 272        Tango2 1.800365 9.595410e-02      DISS
    ## 273       Gadd45b 1.485754 9.833361e-02      DISS
    ## 274         Itpr3 3.118733 9.916376e-02      DISS

    suptable %>% dplyr::filter(direction == "HOMO")

    ##        gene       lfc        padj direction
    ## 1    Csrnp3 -1.457990 0.006552580      HOMO
    ## 2   Stxbp5l -1.940397 0.009065359      HOMO
    ## 3    Gabrb1 -1.074515 0.010299888      HOMO
    ## 4    Rgs7bp -1.356723 0.011416902      HOMO
    ## 5    Grin2b -1.656666 0.013241261      HOMO
    ## 6     Sorl1 -1.408413 0.013455971      HOMO
    ## 7     Birc6 -1.324620 0.016029141      HOMO
    ## 8     Cdkl5 -1.665223 0.016029141      HOMO
    ## 9  Ralgapa1 -1.186125 0.016300191      HOMO
    ## 10    Rc3h2 -1.965015 0.018685627      HOMO
    ## 11    Epha6 -1.702348 0.020058278      HOMO
    ## 12   Kcnma1 -1.918591 0.020641917      HOMO
    ## 13   Nos1ap -1.233896 0.022106195      HOMO
    ## 14   Plppr4 -1.247816 0.024923161      HOMO
    ## 15   Celsr2 -1.134917 0.025399573      HOMO
    ## 16   Grin2a -1.659562 0.026906675      HOMO
    ## 17    Lrrc7 -1.780252 0.027317499      HOMO
    ## 18    Cadm2 -1.209343 0.029967447      HOMO
    ## 19    Magi2 -1.273744 0.032118182      HOMO
    ## 20  Rasgrp1 -1.083835 0.037815551      HOMO
    ## 21    Rock2 -1.049676 0.038789415      HOMO
    ## 22    Asap1 -1.880316 0.040158237      HOMO
    ## 23 Cacna2d1 -1.223018 0.041651581      HOMO
    ## 24  Gm42878 -3.795602 0.043652913      HOMO
    ## 25    Kcnj6 -2.392068 0.043652913      HOMO
    ## 26    Adcy9 -1.420367 0.047426950      HOMO
    ## 27     Fat3 -2.180080 0.054003118      HOMO
    ## 28    Taok1 -1.308505 0.054003118      HOMO
    ## 29   Adgrl3 -1.354263 0.056973405      HOMO
    ## 30   Slc4a8 -1.042741 0.060991693      HOMO
    ## 31    Atxn1 -1.519674 0.067359145      HOMO
    ## 32  Arfgef3 -1.517249 0.068353199      HOMO
    ## 33 Arhgap32 -1.277570 0.068546662      HOMO
    ## 34    Pde4d -1.752131 0.069451358      HOMO
    ## 35     Taf1 -1.210848 0.071459179      HOMO
    ## 36    Clock -1.240066 0.072533381      HOMO
    ## 37   Atp2b1 -1.206102 0.073961672      HOMO
    ## 38    Prr36 -1.756586 0.073961672      HOMO
    ## 39  Cacna1e -1.277796 0.080501037      HOMO
    ## 40     Gnaq -1.176176 0.082813194      HOMO
    ## 41    Syt14 -1.498771 0.083519563      HOMO
    ## 42    Garem -1.385046 0.088762037      HOMO
    ## 43    Pde4c -3.998829 0.090121734      HOMO
    ## 44   Atp2c1 -1.129268 0.095954102      HOMO
    ## 45   Cacnb2 -1.161053 0.095954102      HOMO
    ## 46   Ubqln1 -1.057072 0.097402529      HOMO
    ## 47     Atrx -1.058999 0.099163760      HOMO

**Supplemental Table 2. Molecules implicated in hippocampal LTP from
Sanes and Lichtman 1999.**

<https://github.com/raynamharris/DissociationTest/blob/master/data/SanesLichtman.csv>

    SanesTable <- read.csv("~/Github/DissociationTest/data/SanesLichtman.csv", header = T)
    SanesTable

    ##                                                Sanes...Lichtman.molecules
    ## 1                                                     GLUTAMATE RECEPTORS
    ## 2                                                            GluR1; GluR2
    ## 3                                          mGluR1; mGluR4; mGluR5; mGluR7
    ## 4                                          NMDA NR2A; NMDA NR2D; NMDA NR1
    ## 5            OTHER NEUROTRANSMITTERS; NEUROMODULATORS AND THEIR RECEPTORS
    ## 6                               norepinephrine and b-adrenergic receptors
    ## 7                                    adenosine and adenosine 2A receptors
    ## 8                                      dopamine and D1 dopamine receptors
    ## 9                                           mu and delta opioid receptors
    ## 10                                                acetylcholine receptors
    ## 11                                                   muscarinic receptors
    ## 12                                                         GABA receptors
    ## 13                                                       GABA-B receptors
    ## 14                                                   cannabinoid receptor
    ## 15                                  orphanin NQ and nocioceptin receptors
    ## 16                                                    serotonin receptors
    ## 17                                                           endothelin-1
    ## 18                                gamma-aminobutyric acid (GHB) receptors
    ## 19  INTERCELLULAR MESSENGERS; THEIR SYNTHETIC ENZYMES AND THEIR RECEPTORS
    ## 20                                                                     CO
    ## 21                                                                     NO
    ## 22                                                                    EGF
    ## 23                                                              basic FGF
    ## 24                                                             superoxide
    ## 25                                                             neuregulin
    ## 26                                                                  erbB4
    ## 27                                                                    NGF
    ## 28                                                                   BDNF
    ## 29                                                                   TrkB
    ## 30                                                             nNOS; eNOS
    ## 31                                                       arachidonic acid
    ## 32                                             platelet activating factor
    ## 33                                                     interleukin 1 beta
    ## 34                                                                    H2S
    ## 35                                                           beta activin
    ## 36                                    CALCIUM/CALMODULIN BINDING PROTEINS
    ## 37                                                             calmodulin
    ## 38                                                        RC3/neurogranin
    ## 39                                                             calretinin
    ## 40                                                 GAP43/B50/neuromodulin
    ## 41                                                                   S100
    ## 42                                                           ION CHANNELS
    ## 43                                                L-type calcium channels
    ## 44                              olfactory cyclic nucleotide-gated channel
    ## 45                               VESICLE- AND SYNAPSE-ASSOCIATED PROTEINS
    ## 46                                                          synaptophysin
    ## 47                                                                 a-SNAP
    ## 48                                                                   VAMP
    ## 49                                                                  rab3a
    ## 50                                                            syntaxin 1B
    ## 51                                                             Synapsin I
    ## 52                                                                SNAP 25
    ## 53                                                                 PSD-95
    ## 54                                                  TRANSCRIPTION FACTORS
    ## 55                                            Retinoic acid receptor beta
    ## 56                                                                   CREB
    ## 57                                                       Krox 20; Krox 24
    ## 58                                                     ADHESION MOLECULES
    ## 59                                                                  ephA5
    ## 60                                                               ephrinA5
    ## 61                                                                   NCAM
    ## 62                                                 E-cadherin; N-cadherin
    ## 63                                                                  thy-1
    ## 64                                                          telencephalin
    ## 65                                                               L1/NgCAM
    ## 66                                                     HB-GAM/pleitrophin
    ## 67                                                              integrins
    ## 68                                            integrin-associated protein
    ## 69                                                             tenascin-C
    ## 70                                                          CARBOHYDRATES
    ## 71                                                        Polysialic acid
    ## 72                                                        Ganglioside GM1
    ## 73                                                       Ganglioside GQ1B
    ## 74                                                                KINASES
    ## 75                                         inositol-triphosphate-3-kinase
    ## 76                                                                   MAPK
    ## 77                                                                    src
    ## 78                                                                    fyn
    ## 79                                      protein kinase A C beta 1 subunit
    ## 80                                       protein kinase A RI beta subunit
    ## 81                                                 protein kinase C-gamma
    ## 82                                                       protein kinase G
    ## 83                                                  protein kinase M-zeta
    ## 84                                                   CaM kinase I; II; IV
    ## 85                                                    ecto-protein kinase
    ## 86                                         PROTEASES AND THEIR INHIBITORS
    ## 87                                                                calpain
    ## 88                                                            calpastatin
    ## 89                                                       protease nexin 1
    ## 90                                           tissue plasminogen activator
    ## 91                                                                plasmin
    ## 92                                                 E6-AP ubiquitin ligase
    ## 93                                                          OTHER ENZYMES
    ## 94                                                       phospholipase A2
    ## 95                                                   phospholipase C beta
    ## 96                                                  phospholipase C gamma
    ## 97                                                ADP ribosyl transferase
    ## 98                                                            calcineurin
    ## 99                                                  protein phosphatase I
    ## 100                                                  acetylcholinesterase
    ## 101                                                     adenylate cyclase
    ## 102                                                     guanylate cyclase
    ## 103                                                         MISCELLANEOUS
    ## 104                                                       Spectrin/fodrin
    ## 105                                                                  GFAP
    ## 106                                                      Stathmin RB3/XB3
    ## 107                                      EBI-1 G protein-coupled receptor
    ## 108                                        Mas G-protein coupled receptor
    ## 109                                                                  Vesl
    ##                                                                                                                                                                                                                              related.transcripts
    ## 1                                                                                                                                                                                                                                               
    ## 2                                                                                                                                                                                                                                  Gria1; Gria2 
    ## 3                                                                                                                                                                                                                        Grm1; Grm4; Grm5; Grm7 
    ## 4                                                                                                                                                                                                                         Grin1; Grin2a; Grin2d 
    ## 5                                                                                                                                                                                                                                               
    ## 6                                                                                                                                                                                                                           Adrb1; Adrb2; Adrb3 
    ## 7                                                                                                                                                                                                 Adra1a; Adra1b; Adra1d; Adra2a; Adra2b; Adra2c
    ## 8                                                                                                                                                                                                                                      Th; Drd1 
    ## 9                                                                                                                                                                                                                                   Oprm1; Oprd1
    ## 10                                                                                                                                                                                                Chrna1; Chrna7; Chrna3; Chrnb1; Chrnb2; Chrnb3
    ## 11                                                                                                                                                                                                             Chrm1; Chrm2; Chrm3; Chrm4; Chrm5
    ## 12                                                                                                                                                                                                       Gabra1; Gabra2;  Gabra3; Gabra5; Gabra6
    ## 13                                                                                                                                                                                                                        Gabrb1; Gabrb2; Gabrb3
    ## 14                                                                                                                                                                                                                                    Cnr1; Cnr2
    ## 15                                                                                                                                                                                                                                  Pnoc; Oprl1;
    ## 16                                                                                                                                                      Htr1a; Htr1b; Htr1f; Htr2a; Htr2c; Htr2b; Htr3a; Htr3b; Htr5a; Htr5b; Htr7; Htr6; Htr4; 
    ## 17                                                                                                                                                                                                                                          Edn1
    ## 18                                                                                                                                                                                                                               Gabrr1;  Gabbr1
    ## 19                                                                                                                                                                                                                                              
    ## 20                                                                                                                                                                                                                                            NA
    ## 21                                                                                                                                                                                                                                            NA
    ## 22                                                                                                                                                                                                                                           Egf
    ## 23                                                                                                                                                                                                                                          Fgf2
    ## 24                                                                                                                                                                                                                           NA                 
    ## 25                                                                                                                                                                                                                              Nrg1; Nrg2; Nrg3
    ## 26                                                                                                                                                                                                                                         Erbb4
    ## 27                                                                                                                                                                                                                                           Ngf
    ## 28                                                                                                                                                                                                                                          Bdnf
    ## 29                                                                                                                                                                                                                                        Ntrk2 
    ## 30                                                                                                                                                                                                                                    Nos1; Nos3
    ## 31                                                                                                                                                                                                                                            NA
    ## 32                                                                                                                                                                                                                                            NA
    ## 33                                                                                                                                                                                                                                         Il1b 
    ## 34                                                                                                                                                                                                                                            NA
    ## 35                                                                                                                                                                                                                                         Inhba
    ## 36                                                                                                                                                                                                                                              
    ## 37                                                                                                                                                                                                                           Calm1; Calm2; Calm3
    ## 38                                                                                                                                                                                                                                          Nrgn
    ## 39                                                                                                                                                                                                                                 Calb1; Calb2 
    ## 40                                                                                                                                                                                                                                        Gap43 
    ## 41                                                                                                                                                                                                                                         S100b
    ## 42                                                                                                                                                                                                                                              
    ## 43                                                                                                                                                                                 Cacna1c; Cacna1d; Cacna1s; Cacna1f; Cacna1b; Cacna1a; Cacna1e
    ## 44                                                                                                                                                                                                                                        Cnga2 
    ## 45                                                                                                                                                                                                                                              
    ## 46                                                                                                                                                                                                                                          Syp 
    ## 47                                                                                                                                                                                                                                          Napa
    ## 48                                                                                                                                                                                                      Vamp1; Vamp2; Vamp3; Vamp4; Vamp5; Vamp8
    ## 49                                                                                                                                                                                                                                         Rab3a
    ## 50                                                                                                                                                                                                                                        Stx1b;
    ## 51                                                                                                                                                                                                                                         Syn1 
    ## 52                                                                                                                                                                                                                                        Snap25
    ## 53                                                                                                                                                                                                                                          Dlg4
    ## 54                                                                                                                                                                                                                                              
    ## 55                                                                                                                                                                                                                                         Rarb 
    ## 56                                                                                                                                                                                                                                         Creb1
    ## 57                                                                                                                                                                                                                                   Egr1; Egr2 
    ## 58                                                                                                                                                                                                                                              
    ## 59                                                                                                                                                                                                                                        Epha5 
    ## 60                                                                                                                                                                                                                                        Efna5 
    ## 61                                                                                                                                                                                                                                  Ncam1; Ncam2
    ## 62                                                                                                                                                                                                                                   Cdh1; Cdh2 
    ## 63                                                                                                                                                                                                                                         Thy1 
    ## 64                                                                                                                                                                                                                                         Icam5
    ## 65                                                                                                                                                                                                                                        L1cam 
    ## 66                                                                                                                                                                                                                                          Ptn 
    ## 67   Itga1; Itga10;Itga11; Itga2;  Itga2b;Itga3;  Itga4; Itga5;  Itga6;  Itga7;  Itga8;  Itga9;  Itgad;  Itgae; Itgal;  Itgam;  Itgav;  Itgax;  Itgb1;  Itgb1bp1; Itgb2;  Itgb2l; Itgb3;  Itgb3bp;  Itgb4; Itgb5;  Itgb6; Itgb7;  Itgb8;  Itgbl1
    ## 68                                                                                                                                                                                                                                          Cd47
    ## 69                                                                                                                                                                                                                                           Tnc
    ## 70                                                                                                                                                                                                                                              
    ## 71                                                                                                                                                                                                                                            NA
    ## 72                                                                                                                                                                                                                                            NA
    ## 73                                                                                                                                                                                                                                            NA
    ## 74                                                                                                                                                                                                                                              
    ## 75                                                                                                                                                                                                                           Itpka; Itpkb; Itpkc
    ## 76                                                                                                                                                              Mapk1; Mapk10; Mapk11; Mapk12; Mapk14;  Mapk3; Mapk4; Mapk6; Mapk7; Mapk8; Mapk9
    ## 77                                                                                                                                                                                                                                          Src 
    ## 78                                                                                                                                                                                                                                           Fyn
    ## 79                                                                                                                                                                                                                                        Prkacb
    ## 80                                                                                                                                                                                                                                       Prkar1b
    ## 81                                                                                                                                                                                                                                        Prkcg 
    ## 82                                                                                                                                                                                                                                        Prkg1 
    ## 83                                                                                                                                                                                                                                        Prkcz 
    ## 84                                                                                                                                                                                                                           Camk1; Camk2; Camk4
    ## 85                                                                                                                                                                                                                                            NA
    ## 86                                                                                                                                                                                                                                              
    ## 87                                                                                                                                                Capn1; Capn10; Capn11; Capn12; Capn13; Capn15; Capn2; Capn3; Capn5; Capn6; Capn7; Capn8; Capn9
    ## 88                                                                                                                                                                                                                                         Cast 
    ## 89                                                                                                                                                                                                                                     Serpine2 
    ## 90                                                                                                                                                                                                                                         Plat 
    ## 91                                                                                                                                                                                                                                           Plg
    ## 92                                                                                                                                                                                                                                         Ube3a
    ## 93                                                                                                                                                                                                                                              
    ## 94                                                                      Pla2g10; Pla2g12a; Pla2g12b; Pla2g15; Pla2g16; Pla2g1b; Pla2g2a; Pla2g2c; Pla2g2d; Pla2g2e; Pla2g2f; Pla2g3;  Pla2g4a; Pla2g4b; Pla2g4e; Pla2g4f; Pla2g5; Pla2g6; Pla2g7
    ## 95                                                                                                                                                                                                                    Plcb1; Plcb2; Plcb3; Plcb4
    ## 96                                                                                                                                                                                                                                  Plcg1; Plcg2
    ## 97                                                                                                                                                                                                                                        Parp1 
    ## 98                                                                                                                                                                                                                        Ppp3ca; Ppp3cb; Ppp3cc
    ## 99                                                                                                                                                                                                                                        Phpt1 
    ## 100                                                                                                                                                                                                                                         Ache
    ## 101                                                                                                                                                                        Adcy1; Adcy10; Adcy2; Adcy3; Adcy4; Adcy5; Adcy6; Adcy7; Adcy8; Adcy9
    ## 102                                                                                                                                                                           Gucy1a2; Gucy1a3; Gucy1b2; Gucy1b3; Gucy2c; Gucy2d; Gucy2e; Gucy2g
    ## 103                                                                                                                                                                                                                                             
    ## 104                                                                                                                                                                                                                              Sptan1; Sptbn1 
    ## 105                                                                                                                                                                                                                                         Gfap
    ## 106                                                                                                                                                                                                                                        Stmn4
    ## 107                                                                                                                                                                                                                                         Ccr7
    ## 108                                                                                                                                                                                                                                         Mas1
    ## 109                                                                                                                                                                                                                      Homer1; Homer2; Homer3

Outdated Video Abstract
-----------------------

I made this [short video](https://www.youtube.com/watch?v=taeAqimxXWo)
explaining how to use this GitHub repo when I submitted the first draft
to the journal Hippocampus and posted a pre-print on BioRxiv.

[![screenshot](./figures/screenshot.png)](https://www.youtube.com/watch?v=taeAqimxXWo)
