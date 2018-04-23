# Hippocampal transcriptomic responses to cellular dissociation
 
## Title Page 
Rayna M. Harris1,2, Hsin-Yi Kao2,3, Juan Marcos Alarcon2,4,5, Hans A. Hofmann1,2, and André A. Fenton2,3,5,6
 
1 Department of Integrative Biology, Center for Computational Biology and Bioinformatics, Institute for Cellular and Molecular Biology, The University of Texas at Austin, Austin, TX, USA
2 Neural Systems & Behavior Course, Marine Biological Laboratory, Woods Hole, MA, USA
3 Center for Neural Science, New York University, New York, NY, USA
4 Department of Pathology, State University of New York, Downstate Medical Center, Brooklyn, NY, USA New York, NY, USA
5 The Robert F. Furchgott Center for Neural and Behavioral Science, State University of New York, Downstate Medical Center, Brooklyn, NY, USA New York, NY, USA
6 Neuroscience Institute at the New York University Langone Medical Center, New York University, New York, NY, USA

Number of Pages: 10
Number of Figures: 2
Number of Tables: 2
Number of Words: Abstract: ; Main Text: 3123
Corresponding Author
André A. Fenton: afenton@nyu.edu
Orcid ID: 0000-0002-5063-1156
Grant sponsor: NINDS; Grant number: NS091830 to JMA.
Grant sponsor: NSF; Grant number: IOS-1501704 to HAH
Grant sponsor: NIMH; Grant number: 5R25MH059472-18
Grant sponsor: Helmsley Foundation Advanced Training at the Interface of Biology and Computational Science to MBL - Helmsley Innovation Award to AAF and HAH
Grant sponsor: Grass Foundation to MBL
Grant sponsor: Michael Vasinkevich to AAF
Competing interests: The authors declare no competing interests.
Keywords: hippocampus, transcriptomics, genomics, reproducible research, learning, memory, gene, mouse, cellular dissociation

 
## ABSTRACT

Background: Single-neuron gene expression studies may be especially important for understanding nervous system structure and function because of the neuron-specific functionality and plasticity that defines functional neural circuits. Cellular dissociation is a prerequisite technical manipulation for single-cell and single cell-population studies, but the extent to which the cellular dissociation process cells affects neural gene expression has not been determined. This information is necessary for interpreting the results of experimental manipulations that affect neural function such as learning and memory.
 
Methods: The goal of this research was to determine if dissociation itself modifies the transcriptome. We compared tissue level expression of microdissected samples from the dentate gyrus (DG), CA3, and CA1 subfields of the mouse hippocampus either prepared by a standard tissue homogenization protocol or subjected to a cellular dissociation procedure. We used the Illumina HiSeq platform for sequencing, Kallisto for transcript abundance estimation, DESeq2 for differential gene expression profiling, and GO_MWU for analysis of gene ontology. Raw reads, results, and code are available at GEO and on GitHub.
 
Results: This report shows that the process of cellular dissociation compared to homogenization alters only 2% of the of tissue-level transcriptome of hippocampal subfields. Genes related to cellular stress are affected by dissociation relative to homogenization. Few canonical genes underlying learning and memory are not found to be in the list differentially expressed genes.
 
Discussion: This study suggests that cellular dissociation does not affect genes or molecular functions canonically related to learning and memory. However, sample preparation can affect gene expression profiles, which might confound interpretation of results depending on the research question. This study is important for the investigation of any complex tissues as research effort moves from subfield level analysis to single cell analysis of gene expression.

 
## BACKGROUND 
Nervous systems are comprised of diverse cell types that include neurons, glia, and vascular cells, each serving distinct functions and thus expressing different genes. Even within anatomically-defined subfields of the brain, there are identifiable subclasses of neurons that belong to distinct functional circuits (Mizuseki et al. 2011; Danielson et al. 2016; Namburi et al. 2015). This diversity is evident and documented in the Allen Brain Atlas, and is even greater when we consider that specific cells within a functional class can be selectively altered by neural activity in the recent or distant past. All this diversity implies distinctive gene expression, at the level of single neurons, and such considerations may strongly curtail interpretations of gene expression studies that use mixtures of cells or microdissected tissue samples.
 
Recent advances in tissue harvesting and processing, as well as in sequencing technologies have allowed detailed analyses of genome-scale gene expression profiles at the level of single cell populations, in the context of brain and behavior studies (Chalancon et al., 2012; Harris and Hofmann, 2014; Mo et al., 2015; Lacar et al., 2016). These approaches have led to systems-level insights into the molecular substrates of neural function, along with the discovery and validation of candidate pathways regulating physiology and behavior (Cembrowski et al., 2016). The complexity of some tissues can confound the interpretation of transcriptome data collected from bulk samples containing hundreds to tens of thousands of cells that represent numerous cellular subclasses at different levels of diversity. These difficulties can be minimized by careful experimental design governing both data collection and data analysis. To complement this effort, and optimize experimental designs, it is necessary to understand the extent to which the treatment of tissue samples prior to transcriptome analysis might confound interpretation of the results.
 

The goal of this research was to determine if dissociation itself alters the transcriptome. We did not compare single-cell RNA-seq data to bulk tissue RNA-seq data because that is orthogonal to the present research question. Instead, we compared tissue level expression of microdissected samples from the dentate gyrus (DG), CA3, and CA1 hippocampal subfields prepared by a standard homogenization protocol to corresponding samples that were dissociated as if they were being prepared for single-cell sequencing (Fig 1A). We used the Illumina HiSeq platform for sequencing, Kallisto for transcript abundance estimation (Bray et al., 2016), DESeq2 for differential gene expression profiling (Love et al., 2014), and GO_MWU for analysis of gene ontology (Wright et al., 2015). Data and code are available at NCBI’s Gene Expression Omnibus Database (accession number GSE99765). The data and code are available on GitHub (https://github.com/raynamharris/DissociationTest). A more detailed description of the methods are provided in the supplementary “detailed methods” section.
 
## RESULTS AND DISCUSSION
This dataset contains of subfield-specific transcriptome data from three subfields of the hippocampus (CA1, CA3, DG) subjected to one of two treatments (homogenize (HOMO) or dissociated (DISS).  The hypothesis that treatment effects will affect all brain cell types similarly, predicts that gene expression differences between dissociated and homogenized samples would not be specific to the hippocampal subfields, however it is known for, example that the CA1 region is more vulneralble to anoxia than other hippocampus cell regions (REFs) , in which case region-specific differences in the influence of treatment type might be expected. If there are preparation induced sub-field specific differences the transcriptomes would separate into six distinct clusters, one for each region and treatment. However, principal component analysis does not suggest region-specific responses to tissue preparation (Fig. 1B). In this analysis PC1 clearly separates the DG samples from the CA1 and CA3 samples. A two-way treatment-by-region ANOVA confirmed a significant effect of region (F2,11= 17.69; p = 0.0004). Post hoc Tukey tests confirmed CA1 = CA3 < DG. PC2 does vary significantly with treatment (F1,12=6.125; p = 0.03) and accounts for 22% of the variation in gene expression; however, there is no treatment x region interaction. None of the other PCs showed significant variation according to either region or treatment.
 
Differences in hippocampal subfield gene expression are well known. Lein et al., 2004 provide a table with top 100 differnetially expressed genes by region as iddentified through in situ hybridization. Hawrylycz et al., 2012 used hierarchical clustering to visuzalize the top 5000 diffentially expressed genes (P < 0.01) across hippocampal subfields. Cembrowksi show that spatially explicity RNA-seq experiments gave good agreement with histological data, correctly predicting the enriched populations in ~81% of cases (124/153 genes where coronal histological images were available).
 
Thus we expected dozens to thousands of differentially expressed genes according to subfield. We only observed hundreds of differentially expressed genes (Table 1). The main reason for detecting fewer DEGs may be a lack of power (small sample size). The homogenized samples appear to have less variance than the dissociated samples (Fig. 1B). This increased variance among dissociated samples may have contributed to the reduced sub-field specific differences.
 
We found that 2.9% of the transcriptome is differentially expressed between CA1 and DG, with a roughly symmetric distribution of differential gene expression (Fig. 1C). In comparison, we found that 2.1% of the genes measured changed in response to cellular dissociation (Table 1). We found that most differentially expressed genes were upregulated (288 genes) rather than down-regulated (56 genes) in response to dissociation (Fig. 1C). The asymmetry could be due to activation of cellular stress response.  The directionality suggests that dissociated cells are not simply dying at larger numbers, since otherwise one might expect a depletion of gene expression.

Of the 344 differentially expressed genes, we identified only 3 genes (Gria2, Grin2a, Grin2b) of of a list of 49 variants of these well known learning and memory genes (Ncs, Nsf, Gria, Grin, Grim, Dlg, Prkc, Camk2, Fmr1, Creb). Of the top 30 most differentially expressed genes (Fig. 2A), none appear to be involved in the glutamatergic or calcium signally pathways that are well known to be involved in learning and memory (Huber et al., 2000). 
 
An analysis of functional enrichment allows a general examination of the molecular function and cellular components that are affected by the techniques used to process tissue samples (Fig. 2B. We found that the process of chemical dissociation relative to tissue homogenization results in a significant upregulation of molecular processes related to ribosomal activity and oxidoreductase activity and a down regulation of ligase activity). Upregulation of ribosomal activity can be used as a proxy for translation and oxidoreductase activity as a proxy for oxidative stress response (Ezraty 2017). Thus, it appears that dissociation results in an upregulation of ribosomal activity (a proxy for translation) and oxidoreductase activity (oxidative stress response). The significance levels of all GO terms is provided in Supplementary Table (2) so that researchers can query the data frame for molecular functions of interest. 
 
Importantly, we did not conduct single-neuron sequencing in this study. Even though single-cell RNA-sequencing studies are increasingly routine in basic research, a single-cell study would not have made it possible to address our present research question of how the process of cellular dissociation affects gene expression relative to tissue homogenization. 

In conclusion, we set out to identify the extent to which the process of cellular dissociation – which necessarily precedes single cell analysis of complex tissues – has an appreciable effect on our ability to detect biologically meaningful variation in hippocampal gene expression. The purpose of this study was to test whether analysis of gene expression in hippocampus subfields is changed by tissue preparation procedures (cellular dissociation versus homogenization). This is potentially important because it is increasingly necessary to dissociate cells in tissue samples for single cell or single cell-type studies. These findings provide insight into how cellular manipulations influence gene expression. It may be useful to first prepare tissues with transcription and translation blockers like puromycin and actinomycin to arrest gene expression activity before cellular dissociation (Flexner et al., 1963; Solntseva and Nikitin,  2012). Further research is clearly needed to uncover the influence of other variables and their interaction with methodological factors on gene expression variation in the hippocampus.

## ACKNOWLEDGMENTS
We thank members of the Hofmann and Fenton Labs, Boris Zemelman, Laura Colgin, and Misha Matz for helpful discussions. We thank Dennis Wylie for insightful comments on earlier versions of this manuscript. We thank Promega Corporation for generously donating molecular reagents for RNA isolation. We thank the GSAF for library preparation and sequencing. The bioinformatic workflow was inspired heavily by Center for Computational Biology’s Bioinformatics Curriculum and Software Carpentry Curriculum on the Unix Shell, Git for Version Control, and R for Reproducible Research. This work is supported by a Society for Integrative Biology (SICB) Grant in Aid of Research (GIAR) grant and a UT Austin Graduate School Continuing Fellowship to RMH; a generous gift from Michael Vasinkevich to AAF; NIH-NS091830 to JMA, IOS-1501704 to HAH; NIMH-5R25MH059472-18, the Grass Foundation, and the Helmsley Charitable Trust. The authors declare no competing interests.
 
## DETAILED METHODS
All animal care and use complies with the Public Health Service Policy on Humane Care and Use of Laboratory Animals and were approved by the New York University Animal Welfare Committee and the Marine Biological Laboratory Institutional Animal Care and Use Committee.
 
A 1-year-old female C57BL/6J mouse was taken from its cage, anesthetized with 2% (vol/vol) isoflurane for 2 minutes and decapitated. Transverse 300 μm brain slices were cut using a vibratome (model VT1000 S, Leica Biosystems, Buffalo Grove, IL) and incubated at 36°C for 30 min and then at room temperature for 90 min in oxygenated artificial cerebrospinal fluid (aCSF in mM: 125 NaCl, 2.5 KCl, 1 MgSO4, 2 CaCl2, 25 NaHCO3, 1.25 NaH2PO4 and 25 Glucose) as in (Pavlowsky and Alarcon, 2012). Tissue adjacent samples were collected from CA1, CA3, and DG, respectively in the dorsal hippocampus by punch (0.25 mm, P/N: 57391; Electron Microscopy Sciences, Hatfield, PA) (Fig 1A).
 
The homogenized (HOMO) samples were processed using the manufacturer instructors for the Maxwell 16 LEV RNA Isolation Kit (Promega, Madison, WI). The dissociated (DISS) samples were incubated for 75 minutes in aCSF containing 1 mg/ml pronase at room temperature, then vortexed and centrifuged. The incubation was terminated by replacing aCSF containing pronase with aCSF. The sample was then vortexed, centrifuged, and gently triturated by 200-μl pipette tip twenty times in aCSF containing 1% FBS. The sample was centrifuged and used as input for RNA isolation using the Maxwell 16 LEV RNA Isolation Kit (Promega, Madison, WI).
 
RNA libraries were prepared by the Genomic Sequencing and Analysis Facility at the University of Texas at Austin using the Illumina HiSeq platform. Raw reads were processed and analyzed on the Stampede Cluster at the Texas Advanced Computing Facility (TACC). Samples yielded an average of 4.9 +/- 2.6 million reads. Quality of the data was checked using the program FASTQC. Low quality reads and adapter sequences were removed using the program Cutadapt (Martin, 2011). We used Kallisto for read pseudoalignment to the Gencode MV11 mouse transcriptome and for transcript counting (Mudge and Harrow, 2015; Bray et al., 2016). (Wickham, 2016). On average, 61.2% +/- 20.8% of the trimmed reads were pseudoaligned to the mouse transcriptome.
 
Kallisto transcript counts were imported into R (R Development Core Team, 2013) and aggregated to yield gene counts using the ‘gene’ identifier from the Gencode reference transcriptome. We used DESeq2 for gene expression normalization and quantification of gene level counts (Love et al., 2014). We used a threshold of a false discovery corrected (FDR) p-value < 0.1. Statistics on the principal component analysis (PCA) were conducted in R (Wickham, 2009; Wickham and Francois, 2015). The hierarchical clustering analysis was conducted and visualized using the R package pheatmap (Kolde, 2015) with the RColorBrewer R packages for color modifications (Neuwirth, 2014). PCA was conducted in R using the DESeq2 and genefilter R packages (Love et al., 2014; Gentleman R, Carey V et al., 2017) and visualized using the ggplot2, cowplot, and RColorBrewer R packages (Wickham, 2009; Neuwirth, 2014; Wilke, 2016). We used GO_MWU for analysis of GO ontology using –log(p-value) as a continuous measure of significance to identify GO categories that are significantly enriched with either up- or down-regulated genes (Wright et al., 2015). No significance cutoff is required for the analysis, but an arbitrary p-value was set to visualize the top 10 most significant GO terms.
 
The raw sequence data and intermediate data files are archived in NCBI’s Gene Expression Omnibus Database (accession numbers GSE99765). The data and code are available on GitHub (https://github.com/raynamharris/DissociationTest), with an archived version at the time of publication available at Zenodo (Harris et al., 2017).


## REFERENCES 
Bray NL, Pimentel H, Melsted P, Pachter L. 2016. Near-optimal probabilistic RNA-seq quantification. Nat Biotechnol [Internet] 34:525–527. Available from: http://www.nature.com/doifinder/10.1038/nbt.3519
Cembrowski MS, Bachman JL, Wang L, Sugino K, Shields BC, Spruston N. 2016a. Spatial Gene-Expression Gradients Underlie Prominent Heterogeneity of CA1 Pyramidal Neurons. Neuron [Internet] 89:351–368. Available from: http://www.ncbi.nlm.nih.gov/pubmed/26777276
Cembrowski MS, Wang L, Sugino K, Shields BC, Spruston N. 2016b. Hipposeq: A comprehensive RNA-seq database of gene expression in hippocampal principal neurons. Elife [Internet] 5:e14997. Available from: http://elifesciences.org/lookup/doi/10.7554/eLife.14997
Chalancon G, Ravarani CNJ, Balaji S, Martinez-Arias A, Aravind L, Jothi R, Babu MM. 2012. Interplay between gene expression noise and regulatory network architecture. Trends Genet [Internet] 28:221–232. Available from: http://www.ncbi.nlm.nih.gov/pubmed/22365642
Danielson NB, Zaremba JD, Kaifosh P, Bowler J, Ladow M, Losonczy A. 2016. Sublayer-Specific Coding Dynamics during Spatial Navigation and Learning in Hippocampal Area CA1. Neuron [Internet] 91:652–665. Available from: http://linkinghub.elsevier.com/retrieve/pii/S0896627316302987
Flexner JB, Flexner LB, Stellar E. 1963. Memory in Mice as Affected by Intracerebral Puromycin. Science (80- ) [Internet] 141:57–59. Available from: http://www.sciencemag.org/cgi/doi/10.1126/science.141.3575.57
Gentleman R, Carey V HW and HF, Gentleman R, Carey V, Huber W, Hahne F. 2017. genefilter: genefilter: methods for filtering genes from high-throughput experiments.
Harris RM, Hofmann HA. 2014. Neurogenomics of Behavioral Plasticity. Adv Exp Med Biol 781:149–168.
Harris RM, Kao H-Y, Alarcon JM, Hofmann HA, Fenton AA, Rayna M Harris Hsin-Yi Kao JMAHAHAAF. 2017. GitHub repository for analyses of hippocampal transcriptomic responses to technical and biological perturbations. Available from: https://doi.org/10.5281/zenodo.815081
Hawrylycz MJ, Lein ES, Guillozet-Bongaarts AL, Shen EH, Ng L, Miller JA, van de Lagemaat LN, Smith KA, Ebbert A, Riley ZL, Abajian C, Beckmann CF, Bernard A, Bertagnolli D, Boe AF, Cartagena PM, Chakravarty MM, Chapin M, Chong J, Dalley RA, Daly BD, Dang C, Datta S, Dee N, Dolbeare TA, Faber V, Feng D, Fowler DR, Goldy J, Gregor BW, Haradon Z, Haynor DR, Hohmann JG, Horvath S, Howard RE, Jeromin A, Jochim JM, Kinnunen M, Lau C, Lazarz ET, Lee C, Lemon TA, Li L, Li Y, Morris JA, Overly CC, Parker PD, Parry SE, Reding M, Royall JJ, Schulkin J, Sequeira PA, Slaughterbeck CR, Smith SMSCSM, Sodt AJ, Sunkin SM, Swanson BE, Vawter MP, Williams D, Wohnoutka P, Zielke HR, Geschwind DH, Hof PR, Smith SMSCSM, Koch C, Grant SGN, Jones AR. 2012. An anatomically comprehensive atlas of the adult human brain transcriptome. Nature [Internet] 489:391–399. Available from: http://www.nature.com/nature/journal/v489/n7416/full/nature11405.html?WT.ec_id=NATURE-20120920#global-mapping-of-transcript-distributions
Huber KM, Kayser MS, Bear MF. 2000. Role for rapid dendritic protein synthesis in hippocampal mGluR-dependent long-term depression. Science 288:1254–1256.
Kolde R. 2015. pheatmap: Pretty Heatmaps.
Lacar B, Linker SB, Jaeger BN, Krishnaswami S, Barron J, Kelder M, Parylak S, Paquola A, Venepally P, Novotny M, O’Connor C, Fitzpatrick C, Erwin J, Hsu JY, Husband D, McConnell MJ, Lasken R, Gage FH. 2016. Nuclear RNA-seq of single neurons reveals molecular signatures of activation. Nat Commun [Internet] 7:11022. Available from: http://www.ncbi.nlm.nih.gov/pubmed/27090946
Lein ES, Zhao X, Gage FH. 2004. Defining a molecular atlas of the hippocampus using DNA microarrays and high-throughput in situ hybridization. J Neurosci [Internet] 24:3879–89. Available from: http://www.jneurosci.org/content/24/15/3879
Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol [Internet] 15:550. Available from: http://genomebiology.com/2014/15/12/550
Martin M. 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal 17:10–12.
Mizuseki K, Diba K, Pastalkova E, Buzsáki G. 2011. Hippocampal CA1 pyramidal cells form functionally distinct sublayers. Nat Neurosci [Internet] 14:1174–1181. Available from: http://www.ncbi.nlm.nih.gov/pubmed/21822270
Mo A, Mukamel EA, Davis FP, Luo C, Henry GL, Picard S, Urich MA, Nery JR, Sejnowski TJ, Lister R, Eddy SR, Ecker JR, Nathans J. 2015. Epigenomic Signatures of Neuronal Diversity in the Mammalian Brain. Neuron [Internet] 86:1369–1384. Available from: http://www.ncbi.nlm.nih.gov/pubmed/26087164
Mudge JM, Harrow J. 2015. Creating reference gene annotation for the mouse C57BL6/J genome assembly. Mamm Genome [Internet] 26:366–378. Available from: http://link.springer.com/10.1007/s00335-015-9583-x
Neuwirth E. 2014. RColorBrewer: ColorBrewer Palettes.
Pavlowsky A, Alarcon JM. 2012. Interaction between Long-Term Potentiation and Depression in CA1 Synapses: Temporal Constrains, Functional Compartmentalization and Protein Synthesis. PLoS One [Internet] 7:e29865. Available from: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0029865
Pavlowsky A, Wallace E, Fenton AA, Alarcon JM. 2017. Persistent modifications of hippocampal synaptic function during remote spatial memory. Neurobiol Learn Mem [Internet] 138:182–197. Available from: http://linkinghub.elsevier.com/retrieve/pii/S1074742716301587
R Development Core Team. 2013. R: a language and environment for statistical computing | GBIF.ORG. Available from: http://www.r-project.org/
Solntseva S, Nikitin V. 2012. Conditioned food aversion reconsolidation in snails is impaired by translation inhibitors but not by transcription inhibitors. Brain Res 1467:42–47.
Weiler IJ, Spangler CC, Klintsova AY, Grossman AW, Kim SH, Bertaina-Anglade V, Khaliq H, De Vries FE, Lambers FAE, Hatia F, Base CK, Greenough WT. 2004. Fragile X mental retardation protein is necessary for neurotransmitter-activated protein translation at synapses. Proc Natl Acad Sci U S A [Internet] 101:17504–17509. Available from: http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=536018&tool=pmcentrez&rendertype=abstract
Wickham H. 2009. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. Available from: http://ggplot2.org
Wickham H. 2016. Flexibly Reshape Data: A Reboot of the Reshape Package.
Wickham H, Francois R. 2015. dplyr: A Grammar of Data Manipulation. R Packag version 050 [Internet]:3. Available from: https://cran.r-project.org/package=dplyr
Wilke CO. 2016. cowplot: Streamlined Plot Theme and Plot Annotations for “ggplot2.”
Wright RM, Aglyamova G V, Meyer E, Matz M V. 2015. Gene expression associated with white syndromes in a reef building coral, Acropora hyacinthus. BMC Genomics [Internet] 16:371. Available from: http://www.ncbi.nlm.nih.gov/pubmed/25956907
Zhong J, Chuang S-C, Bianchi R, Zhao W, Lee H, Fenton AA, Wong RKS, Tiedge H. 2009. BC1 regulation of metabotropic glutamate receptor-mediated neuronal excitability. J Neurosci [Internet] 29:9977–9986. Available from: http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2866649&tool=pmcentrez&rendertype=abstract

 
 
## FIGURES AND TABLES

  
![](../figures/fig_fig1.png)  

#### Figure 1. Cellular dissociation has minors effect on hippocampal gene expression patterns.
A) Experimental design. Two tissue samples were taken from three subfields (CA1, CA3, and DG) from 300 uM brain slices. Two adjacent sample were processed using a homogenization (HOMO) protocol or dissociated (DISS) before processing for tissue level gene expression profiling. B) Dissociation does not yield subfield-specific changes in gene expression. PC1 accounts for 40% of all gene expression variation and by inspection, separates the DG samples from the CA1 and CA3 samples. PC2 accounts for 22% of the variation in gene expression and varies significantly with treatment. C) Fold change and significance values across subfiacross subfieldelds. Genes below the p-value < 0.1 (or –log p-value < 1) are shown in light grey. We found that 222 genes are up-regulated in CA1 (purple circles) while 262 are upregulated in DG (orange circles). D) Fold change and significance treatments. We found that 288 genes are up-regulated in the dissociated treatment group (filled dark grey circles) while 56 are up-regulated in the homogenization control group (open circles). Three genes (Grin2a, Grin2b, and Gria2) are highlighted for their role in learning and memory.
 
 
![](../figures/fig_heatmapGO.png)
 
#### Figure 2. Cellular dissociation affects genes and gene ontologies related to cellular stress.
A) Top 30 differentially expressed genes between dissociated and homogenized tissue. Square boxes at the top color coded by sample (white: homogenized, grey: dissociated, purple: CA1, green: CA, orange: DG. Within the heatmap, level of expression is indicated by the blue-green-yellow gradient with lighter colors indicating increased expression. B) Functional analysis with Gene Ontology. Different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated genes in the dissociated samples (DISS) relative to the homogenized control samples. The tree on the plot is hierarchical clustering of GO categories based on shared genes. Categories with no branch length between them are subsets of each other. The fraction next to GO category name indicates the fraction of genes exceeding the arbitrary absolute value cutoff of 0.05.

#### Table 1. Differentially expressed genes by subfield and treatment.
The total number and percent of differentially expressed genes (DEGs) for four two-way contrasts were calculated using DESeq2. Up-regulated: gene expression is significantly higher (log fold-change > 0; p,0.1) in the first term listed in the contrast. Down-regulated: gene expression is significantly lower (log fold-change < 0; p,0.1)) in the first term listed in the contrast. % DEGs/Total: The sum of up and down regulated genes divided by the total number of genes analyzed (16,709) multiplied by 100%. This table shows that differences between dissociated (DISS) tissue and homogenized (HOMO) tissues are on the same scale as those between the CA1 and DG subfields of the hippocampus.
Two-way contrast	Up-regulated	Down-regulated	% DEGs/Total
CA1 vs DG	222	262	2.9%
CA3 vs DG	45	53	0.5%
CA1 v. CA3	17	1	0.1%
DISS vs HOMO	288	56	2.1%
 
#### Supplemental Table 1. Expression level and fold change of of significant genes (p < 0.1) between dissociated tissue and homogenized tissue. This table shows the log fold change (lfc), p-value (padj), and direction of upregulation for each gene analyzed.

For now at: https://github.com/raynamharris/DissociationTest/blob/master/results/SuppTable1.csv

#### Supplemental Table 2. Gene ontologies of enriched genes. The first row contains the GO category (either MF or CC). The second is the GO term. Also shown are directionally, unumber of enriched genes in that catory out of the total (ratio), and p-value. 

https://github.com/raynamharris/DissociationTest/blob/master/results/GOsignificantcatagories.csv

