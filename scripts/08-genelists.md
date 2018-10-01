Cho et al 2015 used RNA sequencing to quantify transcript levels in the
mouse hippocampus after contextual fear conditioning. The Cho dataset
provides a snapshot of gene expression changes associated with
hippocampal learning and memory 30 min and 4 h after an experiment. The
Cho data are available at
<http://science.sciencemag.org/content/suppl/2015/09/30/350.6256.82.DC1>
and <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72064>. The
file `../data/aac7368-Cho.SM.Table.S2.xls` of differentially expressed
genes was used as a representative dataset for learning and memory
associated gene expression patterns.

In this analysis, I compared the Cho et al differentially expressed
genes (DEGs) to my experimental results (referred to as the Harris data)
to identify the interaction between genes that are differentially
expressed following fear condition and following a chemical
manipulation.

This analysis prints the list of genes that are differentially expressed
in both experiments. The images show that only a few (red dots) of genes
that respond to chemical dissociation are up-regulated or down-regulated
following fear-conditioning.

[Click here to view the source code.](./08-genelists.Rmd)

    library(readxl)
    library(ggplot2)
    library(cowplot)
    library(dplyr)
    library(ggrepel)
    library(plyr)


    # set output file for figures 
    knitr::opts_chunk$set(fig.path = '../figures/08-genelists/', echo = T, message = F, warning=F)

    ### functions for data viz and wrangling

    plotvolcano <- function(df, lfc, log10p, plottitle){
      ggplot(df, aes(x = lfc, y = log10p)) + 
      geom_point(aes(color = factor(direction)), size = 1, alpha = 0.8, na.rm = T) + 
      theme_cowplot(font_size = 12) +
      scale_x_continuous(name="log fold change") +
      scale_y_continuous(name="-log10 p-value",
                         breaks = c(0,5,10,15),
                         limits= c(-1,16)) +
      geom_hline(yintercept = 1,  size = 0.25, linetype = 2 ) +
      theme( legend.title = element_blank(),
             legend.position = "top",
             legend.text=element_text(size=8))   + 
      labs(title = plottitle) 
    }

      
    # format Cho data for plotting and data frame joining
    wrangleCho <- function(df){
      data <- df  
      data$log10p <- -log10(data$pvalue) 
      data <- data %>% select(gene, lfc, log10p, pvalue) %>%
        mutate(direction = ifelse(data$lfc > 1 & data$pvalue < 0.05, 
                            yes = "fear-conditioned", 
                            no = ifelse(data$lfc < -1 & data$pvalue < 0.05, 
                                        yes = "control", 
                                        no = "none")))
      data$direction <- as.factor(data$direction)
      data$direction <- factor(data$direction, c("control", "none", "fear-conditioned"))
      return(data)
    }

    wranglecembroski <- function(df){
      data <- df  
      data$log10p <- -log10(data$pvalue) 
      data <- data %>% select(gene, lfc, log10p, pvalue) %>%
        mutate(direction = ifelse(data$lfc > 1 & data$log10p > 1.3, 
                            yes = "DG", 
                            no = ifelse(data$lfc < -1 & data$log10p > 1.3, 
                                        yes = "CA1", 
                                        no = "none")))
      data$direction <- as.factor(data$direction)
      data$direction <- factor(data$direction, c("CA1", "none", "DG"))
      return(data)
    }


    # overlap with significance for each experiment
    overlap <- function(Chodataframe){
      df <- Chodataframe %>%
      arrange(lfc) 
      
      ij <- full_join(df, dissociation, by = "gene")

      names(ij)[5] <- "Cho"
      names(ij)[9] <- "Harris"
      ij$Harris <- as.character(ij$Harris)
      ij$Harris[is.na(ij$Harris)] <- "absent"
      ij$Cho <- as.character(ij$Cho)
      ij$Cho[is.na(ij$Cho)] <- "absent"

      ij$color <- ifelse(grepl("absent", ij$Harris), ij$Cho, 
                         ifelse(grepl("none", ij$Harris), ij$Cho,
                             ifelse(grepl("HOMO", ij$Harris), "HOMO",
                                    ifelse(grepl("DISS", ij$Harris), "DISS",
                                           NA))))
      ij$color <- as.factor(ij$color) 
      ij$color <- factor(ij$color, levels = c( "absent" , "none", 
                                             "control","fear-conditioned",
                                             "HOMO" , "DISS" ))
     #names(useful)
      
      return(ij)
    }


    ## suzy volcano
    suzyvolcano <- function(df, lfc, log10p, plottitle){
      ggplot(df, aes(x = lfc, y = log10p)) + 
      geom_point(aes(color = factor(color), size = factor(color)), 
                 alpha = 0.8, na.rm = T) +
      theme_cowplot(font_size = 12) +
      scale_size_manual(values=c(1, 1, 1, 1, 2, 2)) +
      scale_x_continuous(name="log fold change") +
      scale_y_continuous(name="-log10 p-value",
                         limits = c(-1, 18),
                         breaks = c(1,6,12,18)) +
      geom_hline(yintercept = 1,  size = 0.25, linetype = 2 ) +
      theme( legend.title = element_blank(),
             legend.position = "bottom",
             legend.text=element_text(size=8))  +
      scale_color_manual(values = c("none" = "grey",
                                    "absent" = "grey",
                                    "HOMO" = "dodgerblue2",
                                    "DISS" = "firebrick1",
                                    "fear-conditioned" = "#018571",
                                    "control" = "#a6611a")) + 
      
      labs(title = plottitle)  
    }

The number of differential gene expression in the harris data set and
then the Cho data sets at 4 h and 30 min.

    ### Harris et al. data wranling
    dissociation <- read.csv("../results/01_dissociation_volcanoTreatment.csv", 
                             header = T, row.names = 1)
    names(dissociation)[5] <- "direction"
    dissociation$direction <- factor(dissociation$direction, c("HOMO", "none", "DISS"))
    summary(dissociation$direction)

    ##  HOMO  none  DISS 
    ##    11 12008   138

    # differential gene expression in the harris data set
    dissociationDEGs <- dissociation %>% filter(direction != "none") 

    # table downloaded from biojupies http://amp.pharm.mssm.edu/biojupies/notebook/zGZJvFUQF

    cembrowski <- read.table("~/Downloads/signature.txt", header = 1)
    head(cembrowski)

    ##   gene_symbol     logFC    AveExpr         t      P.Value    adj.P.Val
    ## 1       PROX1 13.876409  3.6431261 10.811670 2.600771e-09 7.296508e-07
    ## 2      IGFBP5 10.655410  4.0444132  6.024944 1.062618e-05 4.829865e-04
    ## 3       C1QL2 10.009691  1.0010774 14.536638 2.117321e-11 2.088064e-08
    ## 4         DSP  9.995924  0.8950099 10.403023 4.755595e-09 1.199737e-06
    ## 5      STXBP6  9.664500  3.8723598 11.969999 5.138066e-10 2.039186e-07
    ## 6      MFSD2A  9.315131 -0.6966044 13.533442 6.935708e-11 4.614397e-08
    ##           B
    ## 1 11.551588
    ## 2  3.492611
    ## 3 14.128875
    ## 4 10.111786
    ## 5 13.103359
    ## 6 12.302552

    cembrowski <- rename(cembrowski, c(`logFC` ="lfc", 
                       `P.Value` = "pvalue",
                       `gene_symbol` = "gene"))
    cembrowski <- wranglecembroski(cembrowski)
    head(cembrowski)

    ##     gene       lfc    log10p       pvalue direction
    ## 1  PROX1 13.876409  8.584898 2.600771e-09        DG
    ## 2 IGFBP5 10.655410  4.973623 1.062618e-05        DG
    ## 3  C1QL2 10.009691 10.674213 2.117321e-11        DG
    ## 4    DSP  9.995924  8.322795 4.755595e-09        DG
    ## 5 STXBP6  9.664500  9.289200 5.138066e-10        DG
    ## 6 MFSD2A  9.315131 10.158909 6.935708e-11        DG

    volcanoplot1b <- plotvolcano(cembrowski, fourhoursRNA$lfc, fourhoursRNA$log10p, 
                                plottitle = "Cembrowski") +
      scale_color_manual(values = c("none" = "grey",
                                    "DG" = "#018571",
                                    "CA1" = "#a6611a"))
    volcanoplot1b

![](../figures/08-genelists/cembrowksi-1.png)

    ### Cho et al. data at 4 hours 

    S2 <- as.data.frame(readxl::read_excel("../data/aac7368-Cho.SM.Table.S2.xls", skip = 1 ))

    fourhoursRNA <- rename(S2, c(`RNA fold change (4 h/control), log2` ="lfc", 
                       `p-value (4 h)` = "pvalue",
                       `Gene Symbol` = "gene"))

    fourhoursRNA <- wrangleCho(fourhoursRNA)

    volcanoplot1 <- plotvolcano(fourhoursRNA, fourhoursRNA$lfc, fourhoursRNA$log10p, 
                                plottitle = "Cho DEGs - 4 h") +
      scale_color_manual(values = c("none" = "grey",
                                    "fear-conditioned" = "#018571",
                                    "control" = "#a6611a"))
    volcanoplot1

![](../figures/08-genelists/fourhours-1.png)

    ### Cho et al. data at 30 min 

    S2 <- as.data.frame(readxl::read_excel("../data/aac7368-Cho.SM.Table.S2.xls", skip = 1 ))

    thirtyminRNA <- rename(S2, c(`RNA fold change (30 min/control), log2` ="lfc", 
                       `p-value (30 min)` = "pvalue",
                       `Gene Symbol` = "gene"))

    thirtyminRNA <- wrangleCho(thirtyminRNA)

    volcanoplot3 <- plotvolcano(thirtyminRNA, thirtyminRNA$lfc, thirtyminRNA$log10p, 
                                plottitle = "Cho DEGs - 30 min")  +
      scale_color_manual(values = c("none" = "grey",
                                    "fear-conditioned" = "#018571",
                                    "control" = "#a6611a"))
    volcanoplot3

![](../figures/08-genelists/thirtymin-1.png)

    ### Plotting their data with my differential exprssion

    overlap30min <- overlap(thirtyminRNA) # no DEGs
    overlap4h <- overlap(fourhoursRNA)
    summary(overlap4h$color)

    ##           absent             none          control fear-conditioned 
    ##             1566            11403                9                0 
    ##             HOMO             DISS 
    ##               11              138

    fivemin <- rename(S2, c(`RNA fold change (5 min/control), log2` ="lfc", 
                       `p-value (5 min)` = "pvalue",
                       `Gene Symbol` = "gene"))
    fivemin <- wrangleCho(fivemin)

    tenmin <- rename(S2, c(`RNA fold change (10 min/control), log2` ="lfc", 
                       `p-value (10 min)` = "pvalue",
                       `Gene Symbol` = "gene"))
    tenmin <- wrangleCho(tenmin)

    overlaptenmin <- overlap(tenmin) # no DEGs
    overlapfivemin <- overlap(fivemin)
    summary(tenmin$color)

    ## Length  Class   Mode 
    ##      0   NULL   NULL

    summary(fivemin$color)

    ## Length  Class   Mode 
    ##      0   NULL   NULL

    overlap4h %>%
      filter(Cho == "control",
             Harris =="HOMO")  %>%
      select(gene, lfc.x, log10p) %>%
      arrange(gene)

    ## [1] gene   lfc.x  log10p
    ## <0 rows> (or 0-length row.names)

    overlap4h %>%
      filter(Cho == "control",
             Harris =="DISS")  %>%
      select(gene, lfc.x, log10p) %>%
      arrange(gene)

    ## [1] gene   lfc.x  log10p
    ## <0 rows> (or 0-length row.names)

    overlap4h %>%
      filter(Cho == "control")  %>%
      select(gene, lfc.x, log10p) %>%
      arrange(gene)

    ##     gene     lfc.x    log10p
    ## 1   Aqp1 -1.382172 12.836310
    ## 2 Calml4 -1.087474  5.187103
    ## 3  Cldn2 -1.003313  8.527052
    ## 4  Clic6 -1.017852  7.843026
    ## 5  Folr1 -1.102123 12.302898
    ## 6  Kcne2 -1.169621 12.652911
    ## 7 Slc4a5 -1.235330  6.039947
    ## 8 Tmem72 -1.654890  7.715675
    ## 9    Ttr -1.356586  5.464636

    suzyvolcano1 <- suzyvolcano(overlap30min, overlap30min$lfc.x,overlap30min$log10p, 
                                plottitle = NULL)
    suzyvolcano1

![](../figures/08-genelists/overlap-1.png)

    suzyvolcano1 <- suzyvolcano1 + 
        annotate("text", label = "338", x = -0.9, y = 2, size = 3, color = "black") + 
        annotate("text", label = "261", x = 0.9, y = 2, size = 3, color = "black")


    suzyvolcano2 <- suzyvolcano(overlap4h, overlap4h$lfc.x, overlap4h$log10p, 
                                plottitle = NULL)

    suzyvolcano2

![](../figures/08-genelists/overlap-2.png)

    suzyvolcano2 <- suzyvolcano2 + 
        annotate("text", label = "435", x = -0.9, y = 2, size = 3, color = "black") + 
        annotate("text", label = "593", x = 0.9, y = 2, size = 3, color = "black")


    # side by side plots 
    figure4 <- plot_grid(suzyvolcano1, suzyvolcano2, nrow = 1, labels = c('A', 'B'))

    figure4

![](../figures/08-genelists/overlap-3.png)

    ggsave(
      "../figures/figure4.png",
      figure4,
      width = 6,
      height = 4,
      dpi = 1200
    )
