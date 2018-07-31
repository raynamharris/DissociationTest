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

    fourhours <- rename(S2, c(`RPF fold change (4 h/control), log2` ="lfc", 
                       `p-value (4 h)` = "pvalue",
                       `Gene Symbol` = "gene"))

    fourhours <- wrangleCho(fourhours)
    head(fourhours)

    ##      gene         lfc     log10p       pvalue
    ## 1   Npas4  0.38497412 1.84887583 1.416199e-02
    ## 2     Fos -0.74509404 6.15789760 6.951882e-07
    ## 3    Junb -0.18211979 1.44371594 3.599847e-02
    ## 4     Arc -0.13624835 1.12587393 7.483867e-02
    ## 5 Dync1h1  0.02104390 0.10257179 7.896383e-01
    ## 6    Btg2 -0.02308444 0.07133081 8.485339e-01
    ##                                          Description        direction
    ## 1           neuronal PAS domain-containing protein 4 fear-conditioned
    ## 2                               proto-oncogene c-Fos          control
    ## 3                         transcription factor jun-B          control
    ## 4 activity-regulated cytoskeleton-associated protein          control
    ## 5                 cytoplasmic dynein 1 heavy chain 1             none
    ## 6                                       protein BTG2             none

    data <- fourhours

    summary(data$direction)

    ##          control             none fear-conditioned 
    ##              517            10503              511

    volcanoplot <- plotvolcano(data, data$lfc, data$log10p, plottitle = "Cho et. al DEGs at 4 hours")

    myboxplot <- plotboxplot(df = data, data$log10p, data$direction, "Cho et. al DEGs at 4 hours")

    plot_grid(volcanoplot, myboxplot)

![](../figures/08-genelists/fourhours-1.png)

### Overlaping DEGs between Harris et al. and Cho et al at 4 hours post treatment

Top 100 DEGs from Cho et al. 4 hours post contextual fear conditioning
versus control. All DEGs from Harris et al. dissocated versus
homogenized control

    top100overlap(fourhours)

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

    # select and rename the columns of interest
    thirtymin <- rename(S2, c(`RNA fold change (30 min/control), log2` ="lfc", 
                       `p-value (30 min)` = "pvalue",
                       `Gene Symbol` = "gene"))

    thirtymin <- wrangleCho(thirtymin)

    data <- thirtymin

    head(data)

    ##      gene         lfc    log10p       pvalue
    ## 1   Npas4  0.67685150 14.827528 1.487551e-15
    ## 2     Fos  0.83418904 26.810336 1.547619e-27
    ## 3    Junb  0.44526307  9.499079 3.168989e-10
    ## 4     Arc  0.44727520 18.742271 1.810211e-19
    ## 5 Dync1h1 -0.08128958  1.290367 5.124281e-02
    ## 6    Btg2  0.41551260  2.254590 5.564291e-03
    ##                                          Description        direction
    ## 1           neuronal PAS domain-containing protein 4 fear-conditioned
    ## 2                               proto-oncogene c-Fos fear-conditioned
    ## 3                         transcription factor jun-B fear-conditioned
    ## 4 activity-regulated cytoskeleton-associated protein fear-conditioned
    ## 5                 cytoplasmic dynein 1 heavy chain 1          control
    ## 6                                       protein BTG2 fear-conditioned

    summary(data$direction)

    ##          control             none fear-conditioned 
    ##              338            10932              261

    volcanoplot <- plotvolcano(data, data$lfc, data$log10p, plottitle = "Cho et. al DEGs 30 min")

    myboxplot <- plotboxplot(df = data, data$log10p, data$direction, "Cho et. al DEGs at 30 min")

    plot_grid(volcanoplot, myboxplot)

![](../figures/08-genelists/thirtymin-1.png)

### Overlaping DEGs between Harris et al. and Cho et al at 30 min post treatment

Top 100 DEGs from Cho et al. 30 min post contextual fear conditioning
versus control. All DEGs from Harris et al. dissocated versus
homogenized control

    top100overlap(thirtymin)

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
