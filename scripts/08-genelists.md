This script identifies to differntially expreseed genes from Cho et al.
(`../data/aac7368-Cho.SM.Table.S2.xls` from Andre Fenton) and asks if
there is any overlap in the the differentially expressed genes in their
study compared to mine.

### functions for data viz

    plotvolcano <- function(df, lfc, log10p, plottitle){
      ggplot(df, aes(x = lfc, y = log10p)) + 
      geom_point(aes(color = factor(direction), 
                     shape = factor(direction)), 
                 size = 1, alpha = 0.8, na.rm = T) + 
      scale_x_continuous(name="log fold change") +
      scale_y_continuous(name="-log10 p-value") +
      geom_hline(yintercept = 1,  size = 0.25, linetype = 2 ) +
      theme( legend.title = element_blank(),
             legend.position = "bottom")  +
      labs(title = plottitle)
    }


    plotboxplot <- function(df, log10p, direction, plottitle){
      ggplot(df, aes(x = direction, y = log10p, colour = direction)) + 
      geom_boxplot() +
      theme( legend.title = element_blank(),
             legend.position = "bottom")  +
      scale_x_discrete(name="Upregulated in") +
      scale_y_continuous(name="-log10 p-value") +
      labs(title = plottitle) +
      geom_hline(yintercept = 1,  size = 0.25, linetype = 2 )
    } 
      
    # format Cho data for plotting and data frame joining
    wrangleCho <- function(df){
      data <- df  
      data$log10p <- -log10(data$pvalue) 
      data <- data %>% select(gene, lfc, log10p, pvalue, Description) %>%
        mutate(direction = ifelse(data$lfc > 0 & data$log10p > 1, 
                            yes = "fear-conditioned", 
                            no = ifelse(data$lfc < 0 & data$log10p > 1, 
                                        yes = "control", 
                                        no = "none")))
      data$direction <- as.factor(data$direction)
      data$direction <- factor(data$direction, c("control", "none", "fear-conditioned"))
      return(data)
    }


    top100overlap <- function(Chodataframe){
      df <- Chodataframe %>%
      filter(direction != "none")  %>%
      arrange(lfc) %>%
      tail(n=100)
      
      ij <- inner_join(df, dissociation, by = "gene")
      names(ij)[6] <- "Cho"
      names(ij)[10] <- "Harris"
      ij <-  ij %>% select(gene, Cho, Harris, log10p, Description) 
      
      print(ij)
    }

### Harris et al. data

    dissociation <- read.csv("../results/01_dissociation_volcanoTreatment.csv", 
                             header = T, row.names = 1)
    names(dissociation)[5] <- "direction"
    dissociation$direction <- factor(dissociation$direction, c("HOMO", "none", "DISS"))
    summary(dissociation$direction)

    ##  HOMO  none  DISS 
    ##    56 11813   288

    myboxplot <- plotboxplot(dissociation, dissociation$pvalue, 
                             plottitle = "Harris et. al DEGs")

    volcanoplot <- plotvolcano(dissociation, dissociation$lfc, 
                               dissociation$pvalue, plottitle = "Harris et. al DEGs")

    plot_grid(volcanoplot, myboxplot)

![](../figures/08-genelists/dissociation-1.png)

    # differential gene expression in the harris data set
    dissociation <- dissociation %>% filter(direction != "none") 

### Cho et al. data at 4 hours

    S2 <- as.data.frame(readxl::read_excel("../data/aac7368-Cho.SM.Table.S2.xls", skip = 1 ))

    # RNA at 4 hours
    fourhoursRNA <- rename(S2, c(`RNA fold change (4 h/control), log2` ="lfc", 
                       `p-value (4 h)` = "pvalue",
                       `Gene Symbol` = "gene"))
    fourhoursRNA <- wrangleCho(fourhoursRNA)
    summary(fourhoursRNA$direction)

    ##          control             none fear-conditioned 
    ##              435            10503              593

    volcanoplot1 <- plotvolcano(fourhoursRNA, fourhoursRNA$lfc, fourhoursRNA$log10p, plottitle = "RNA-seq at 4 hours")


    # RPF at 4 hours
    fourhoursRPF <- rename(S2, c(`RPF fold change (4 h/control), log2` ="lfc", 
                       `p-value (4 h)` = "pvalue",
                       `Gene Symbol` = "gene"))
    fourhoursRPF <- wrangleCho(fourhoursRPF)
    summary(fourhoursRPF$direction)

    ##          control             none fear-conditioned 
    ##              517            10503              511

    volcanoplot2 <- plotvolcano(fourhoursRPF, fourhoursRPF$lfc, fourhoursRPF$log10p, plottitle = "RFP at 4 hours")


    plot_grid(volcanoplot1, volcanoplot2)

![](../figures/08-genelists/fourhours-1.png)

### Overlaping DEGs between Harris et al. and Cho et al at 4 hours post treatment

Top 100 DEGs from Cho et al. 4 hours post contextual fear conditioning
versus control. All DEGs from Harris et al. dissocated versus
homogenized control

    top100overlap(fourhoursRNA)

    ##     gene              Cho Harris   log10p
    ## 1 Selplg fear-conditioned   DISS 1.556865
    ## 2 Slc2a5 fear-conditioned   DISS 1.212306
    ## 3   Plau fear-conditioned   DISS 1.191287
    ## 4    Fn1 fear-conditioned   DISS 3.416128
    ##                                                         Description
    ## 1                        P-selectin glycoprotein ligand 1 precursor
    ## 2 solute carrier family 2, facilitated glucose transporter member 5
    ## 3                    urokinase-type plasminogen activator precursor
    ## 4                                   fibronectin isoform a precursor

    top100overlap(fourhoursRPF)

    ##     gene              Cho Harris   log10p
    ## 1    Fn1 fear-conditioned   DISS 3.416128
    ## 2 Slc2a5 fear-conditioned   DISS 1.212306
    ## 3   Plau fear-conditioned   DISS 1.191287
    ##                                                         Description
    ## 1                                   fibronectin isoform a precursor
    ## 2 solute carrier family 2, facilitated glucose transporter member 5
    ## 3                    urokinase-type plasminogen activator precursor

### Cho et al. data at 30 min

    # 30 min 
    S2 <- as.data.frame(readxl::read_excel("../data/aac7368-Cho.SM.Table.S2.xls", skip = 1 ))

    ## RNAseq at 30 min
    thirtyminRNA <- rename(S2, c(`RNA fold change (30 min/control), log2` ="lfc", 
                       `p-value (30 min)` = "pvalue",
                       `Gene Symbol` = "gene"))

    thirtyminRNA <- wrangleCho(thirtyminRNA)
    summary(thirtyminRNA$direction)

    ##          control             none fear-conditioned 
    ##              338            10932              261

    volcanoplot3 <- plotvolcano(thirtyminRNA, thirtyminRNA$lfc, thirtyminRNA$log10p, plottitle = "RNA 30 min")

    ## RPF at 30 min
    thirtyminRPF <- rename(S2, c(`RPF fold change (30 min/control), log2` ="lfc", 
                       `p-value (30 min)` = "pvalue",
                       `Gene Symbol` = "gene"))

    thirtyminRPF <- wrangleCho(thirtyminRPF)
    summary(thirtyminRPF$direction)

    ##          control             none fear-conditioned 
    ##              410            10932              189

    volcanoplot4 <- plotvolcano(thirtyminRPF, thirtyminRPF$lfc, thirtyminRPF$log10p, plottitle = "RFP 30 min")


    plot_grid(volcanoplot3, volcanoplot4)

![](../figures/08-genelists/thirtymin-1.png)

### Overlaping DEGs between Harris et al. and Cho et al at 30 min post treatment

Top 100 DEGs from Cho et al. 30 min post contextual fear conditioning
versus control. All DEGs from Harris et al. dissocated versus
homogenized control

    top100overlap(thirtyminRNA)

    ##      gene              Cho Harris    log10p
    ## 1  Csrnp1 fear-conditioned   DISS  1.518961
    ## 2    Ctss fear-conditioned   DISS  1.485500
    ## 3    Cdh9 fear-conditioned   DISS  1.072904
    ## 4   Ostf1 fear-conditioned   DISS  1.315947
    ## 5     Fn1 fear-conditioned   DISS  1.180076
    ## 6  Nfkbia fear-conditioned   DISS  3.666476
    ## 7   Dusp1 fear-conditioned   DISS  4.671673
    ## 8    Btg2 fear-conditioned   DISS  2.254590
    ## 9    Junb fear-conditioned   DISS  9.499079
    ## 10   Fosb fear-conditioned   DISS 12.831460
    ##                               Description
    ## 1  cysteine/serine-rich nuclear protein 1
    ## 2     cathepsin S isoform 1 preproprotein
    ## 3                    cadherin-9 precursor
    ## 4         osteoclast-stimulating factor 1
    ## 5         fibronectin isoform a precursor
    ## 6              NF-kappa-B inhibitor alpha
    ## 7  dual specificity protein phosphatase 1
    ## 8                            protein BTG2
    ## 9              transcription factor jun-B
    ## 10                           protein fosB

    top100overlap(thirtyminRPF)

    ##     gene              Cho Harris    log10p
    ## 1 Rpl36a fear-conditioned   DISS  1.836303
    ## 2 Csrnp1 fear-conditioned   DISS  1.518961
    ## 3 Nfkbia fear-conditioned   DISS  3.666476
    ## 4   Btg2 fear-conditioned   DISS  2.254590
    ## 5  Dusp1 fear-conditioned   DISS  4.671673
    ## 6   Junb fear-conditioned   DISS  9.499079
    ## 7   Fosb fear-conditioned   DISS 12.831460
    ##                              Description
    ## 1             60S ribosomal protein L36a
    ## 2 cysteine/serine-rich nuclear protein 1
    ## 3             NF-kappa-B inhibitor alpha
    ## 4                           protein BTG2
    ## 5 dual specificity protein phosphatase 1
    ## 6             transcription factor jun-B
    ## 7                           protein fosB

### Details

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72064>

Genome-wide profilings of transcriptome and translatome in mouse
hippocampi after contextual fear conditioning

Memory stabilization after learning requires transcriptional and
translational regulations in the brain, yet the temporal molecular
changes following learning have not been explored at the genomic scale.
We here employed ribosome profiling and RNA sequencing to quantify the
translational status and transcript levels in mouse hippocampus
following contextual fear conditioning. We identified 104 genes that are
dynamically regulated. Intriguingly, our analysis revealed novel
repressive regulations in the hippocampus: translational suppression of
ribosomal protein-coding genes at basal state; learning-induced early
translational repression of specific genes; and late persistent
suppression of a subset of genes via inhibition of ESR1/ERα signaling.
Further behavioral analyses revealed that Nrsn1, one of the newly
identified genes undergoing rapid translational repression, can act as a
memory suppressor gene. This study unveils the yet unappreciated
importance of gene repression mechanisms in memory formation.
