Reproducing the Cembrowski et al. study (Part 1)
------------------------------------------------

### Get data

This ([2016 Cembrowski
paper](https://elifesciences.org/content/5/e14997#fig1s30)) is very
similar to my experiment, so I want to compare the two. Like mine, they
compare hippocampal gene expression from dorsal CA1, CA3, and DG sub
Subfields. These cells were identifed through fac sorting to isolate
genetically labeled CA1 and CA3 pyramical neurons and DG granular cells.

Before beginning, I used the following UNIX commands to get their data.

This data was made available here [open source
data](https://www.janelia.org/lab/spruston-lab/resources/source-data-simulation-code-other-resources),
but I downloaded it from the GenBank archive using the following
commands:

    mkdir ../data/Cembrowski
    cd ../data/Cembrowski
    wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74985/suppl/GSE74985_gene_exp.diff.gz'
    wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74985/suppl/GSE74985_genes.fpkm_tracking.gz'
    wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74985/suppl/GSE74985_genes.read_group_tracking.txt.gz'
    wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74985/suppl/GSE74985_mergedCount.txt.gz'
    gunzip *.gz
    gzip GSE74985_genes.fpkm_tracking
    cd ../../bin

The GSE74985\_genes.fpkm\_tracking file is used to extract the gene
names with the corresponding ensembl gene id. The file must be unzipped
for use in R, but it must be zipped in order to store it on GitHub. The
4985\_mergedCount.txt file is used for gene expression analyis in R.

Refine data
-----------

I got this word from "Open Refine", which is a software tool to clean
data before analysis. In order to compare the data to mine directly, I
want to reanalyze the cembrowksi data using my workflow. So, I need to
remove some samples that don't address my question, and re-arrange the
data structure so that it is compatable with my workflow.

Summary of the Cembrowski Data
==============================

Here is an overview of the gene expression counts file from the
Cembrowksi study. Here, the column names contains Ensemble identifies
for genes while the row names contain the sample names.

In these next steps, I use a "geneid" file that I created a while back
(see ) to average the transcript counts for each gene and to have a
column with gene names instead of Ensembl IDs for reference. I also drop
some "non-count" columns and ERCC gene counts.

    # contains a file with the gene name and transcript id
    geneids <- read.table("../data/Cembrowski/geneids.tsv", header=T)
    head(geneids)

    ##              gene_id  gene
    ## 1 ENSMUSG00000000001 Gnai3
    ## 2 ENSMUSG00000000003  Pbsn
    ## 3 ENSMUSG00000000028 Cdc45
    ## 4 ENSMUSG00000000031   H19
    ## 5 ENSMUSG00000000037 Scml2
    ## 6 ENSMUSG00000000049  Apoh

    ## join with geneids so we can look at gene level stuff
    counts$gene_id <- row.names(counts)
    countbygene <- full_join(geneids, counts)

    ## Joining, by = "gene_id"

    ## Warning in full_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
    ## character vector and factor, coercing into character vector

    countbygene <- countbygene %>% 
      filter(gene != "-")
    countbygene <- countbygene[-c(1)] ## keep gene name and counts for samples)

    ## lengthen the dataframe, then wide with gene level sums, then make gene the row name, then round the value to nearest integer
    countbygene <- melt(countbygene, id=c("gene")) 
    countbygene  <- dcast(countbygene, gene ~ variable, value.var= "value", fun.aggregate=sum)
    row.names(countbygene) <- countbygene$gene
    countbygene[1] <- NULL
    countbygene <- round(countbygene)

    # getting ready for DESeq2
    countData <- countbygene 
    head(countData, 9305)  %>% tail()

    ##       dg_d_1 dg_d_2 dg_d_3 dg_v_1 dg_v_2 dg_v_3 ca4_1 ca4_2 ca4_3 ca3_d_1
    ## Fos    22008  21392   7392  16464  17164   6888  4760  4004  4816    1204
    ## Fosb   25116  27104  12040  32200  17612   9828  2016  2800  2912    1708
    ## Fosl1      0      0      0      0      0      0   112     0     0       0
    ## Fosl2  35868   9688  14980  30436  25284  19628 28896 31808 26712   10976
    ## Foxa1      0      0      0      0      0      0     0     0     0       0
    ## Foxa2      0      0      0      0      0      0     0     0     0       0
    ##       ca3_d_2 ca3_d_3 ca3_v_1 ca3_v_2 ca3_v_3 ca2_1 ca2_2 ca2_3 ca1_d_1
    ## Fos      4368    4172    3024    1288    1652  2576  2548   224    1652
    ## Fosb     2632     924    4900    2324    2660  1120  1456   952    4760
    ## Fosl1       0       0       0     112      28     0     0     0       0
    ## Fosl2   11144    8484  100380   91280   47376  9716 15176 12488    6832
    ## Foxa1       0       0       0       0       0     0     0     0       0
    ## Foxa2       0       0       0       0       0     0     0     0       0
    ##       ca1_d_2 ca1_d_3 ca1_v_1 ca1_v_2 ca1_v_3
    ## Fos      2604    2156    1036     308     812
    ## Fosb     3192   14196    1904     924    1988
    ## Fosl1      84      28       0       0       0
    ## Fosl2    8988   13608   21980   18452   44212
    ## Foxa1       0       0       0       0       0
    ## Foxa2       0       0       0       0       0

Now, I use the column names to create a dataframe with columns for
sample name, hippocampal subfield, and dorsal-ventral gradient location.
I however, rename this last column as "Treatment" for ease of use later.
Then, I remove the CA2 and CA4 samples.

    # extract the sample information
    colData <- as.data.frame(colnames(countData))
    names(colData)[1] <- "RNAseqID"
    colData$region <- sapply(strsplit(as.character(colData$RNAseqID),'\\_'), "[", 1)
    colData$location <- sapply(strsplit(as.character(colData$RNAseqID),'\\_'), "[", 2)

    # rename variables and columns
    colData <- rename(colData, c("region"="Subfield"))
    colData <- rename(colData, c("location"="Treatment"))
    colData$Subfield <- plyr::revalue(colData$Subfield, c("ca1"="CA1"))
    colData$Subfield <- plyr::revalue(colData$Subfield, c("ca3"="CA3"))
    colData$Subfield <- plyr::revalue(colData$Subfield, c("dg"="DG"))
    colData$Treatment <- plyr::revalue(colData$Treatment, c("d"="dorsal"))
    colData$Treatment <- plyr::revalue(colData$Treatment, c("v"="ventral"))

    # drop CA2 and CA4 columns from both colData and countData
    colData <- colData %>% 
      filter(grepl("DG|CA1|CA3", Subfield))  %>% 
      droplevels() ## subsets data
    savecols <- as.character(colData$RNAseqID) 
    savecols <- as.vector(savecols) # make it a vector
    countData <- countData %>% select(one_of(savecols))

    # filter lowcounts and NAs
    countData[countData < 2] <- 0
    countData[is.na(countData)] <- 0

    colData$Treatment <- as.factor(colData$Treatment)
    colData$Subfield <- as.factor(colData$Subfield)
    summary(colData)

    ##     RNAseqID  Subfield   Treatment
    ##  ca1_d_1: 1   CA1:6    dorsal :9  
    ##  ca1_d_2: 1   CA3:6    ventral:9  
    ##  ca1_d_3: 1   DG :6               
    ##  ca1_v_1: 1                       
    ##  ca1_v_2: 1                       
    ##  ca1_v_3: 1                       
    ##  (Other):12

    table(colData$Treatment, colData$Subfield)

    ##          
    ##           CA1 CA3 DG
    ##   dorsal    3   3  3
    ##   ventral   3   3  3

    dim(countData)

    ## [1] 34262    18

    write.csv(colData, "../results/04_cembrowksi_colData.csv", row.names = F)
    write.csv(countData, "../results/04_cembrowksi_countData.csv", row.names = T)
