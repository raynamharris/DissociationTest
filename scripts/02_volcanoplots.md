These are the packages I need for my volcano plots.

Here I load the relevant dataframes and set the color palette.

    # dissocation DEGes
    dissociation <- read.csv("../results/volcanoTreatment.csv", header = T, row.names = 1)

    #set levels
    dissociation$direction <- factor(dissociation$direction,
                                 levels = c("HOMO", "DISS", "neither"))

    #set colors
    dissociationcolor <-  c("HOMO" = "black", "DISS" = "black", "neither" = "#d9d9d9")

    # subfield specific degs
    subfield <- read.csv("../results/volcanoCA1DG.csv", header = T, row.names = 1)
    head(subfield)

    ##            gene       pvalue        lfc      padj direction
    ## 1 0610007P14Rik 0.0075257104 -0.5318386 0.9828207   neither
    ## 2 0610009B22Rik 0.0039837752 -0.3593416 0.9908690   neither
    ## 3 0610009L18Rik 0.0067001666 -0.6633154 0.9846907   neither
    ## 4 0610009O20Rik 0.0107561961 -0.4239904 0.9755371   neither
    ## 5 0610010F05Rik 0.0005954242  0.0301684 0.9986299   neither
    ## 6 0610010K14Rik 0.0071664863 -0.4609783 0.9836340   neither

Treatement volcano

    volcanoplot <- ggplot(dissociation, aes(x = lfc, y = pvalue)) + 
      geom_point(aes(color = factor(direction), shape = factor(direction)), 
                 size = 1.5, alpha = 0.8, na.rm = T) + 
      theme_cowplot(font_size = 12, line_size = 0.25) +
      theme(legend.title = element_blank(),
            legend.position=c(.7,.75),
            panel.grid.minor=element_blank(),
            panel.grid.major=element_blank()) +
      scale_color_manual(values = dissociationcolor) +
      scale_x_continuous(name="log fold change",
                          limits = c(-11, 11)) +
      scale_y_continuous(name="-log10 p-value",
                         limits = c(0, 18),
                         breaks = c(1,6,12,18)) +
      geom_hline(yintercept = 1,  size = 0.25, linetype = 2 ) + 
      scale_shape_manual(values = c(1,16,16)) +
        annotate("text", label = "56", x = -10, y = 2, size = 3, color = "black") + 
        annotate("text", label = "288", x = 10, y = 2, size = 3, color = "black")

    volcanoplot  

![](../figures/02_volcanoplots/plot-1.png)

    pdf(file = "../figures/02_volcanoplots/Treatment_volcano_old.pdf", width=2.5, height=2.25)
    plot(volcanoplot)
    dev.off()

    ## quartz_off_screen 
    ##                 2

Caption: Differntial gene expression according to treatment is
asymetric, with more genes enrighted in DISS. only 3 canddiate learning
and memory genes identified.

Plotting CA1 vs. DG volcano plots. The color here is set inside.

    volcanoplot2 <- ggplot(subfield, aes(x = lfc, y = pvalue)) + 
      geom_point(aes(color = factor(direction), shape = factor(direction)), 
                 size = 1.5, alpha = 0.8, na.rm = T) + 
      theme_cowplot(font_size = 12, line_size = 0.25) +
      theme(legend.title=element_blank(),
            legend.position=c(.7,.75),
            panel.grid.minor=element_blank(),
            panel.grid.major=element_blank()) + 
      scale_color_manual(values = c("CA1" = "#7570b3",
                                    "DG" = "#d95f02", 
                                    "neither" = "#d9d9d9")) +   
      scale_x_continuous(name="log fold change",
                         limits = c(-11, 11)) +
      scale_y_continuous(name="-log10 p-value",
                         limits = c(0, 18),
                         breaks = c(1,6,12,18)) +
      geom_hline(yintercept = 1,  size = 0.25, linetype = 2 ) + 
      scale_shape_manual(values = c(16,16,16)) +
      
          annotate("text", label = "222", x = -10, y = 2, size = 3, color = "black") + 
        annotate("text", label = "262", x = 10, y = 2, size = 3, color = "black")
    volcanoplot2 

![](../figures/02_volcanoplots/subfield-1.png)

    pdf(file = "../figures/02_volcanoplots/CA1DG_volcano.pdf", width=2.5, height=2.25)
    plot(volcanoplot2)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    figure2 <- plot_grid(volcanoplot, volcanoplot2,  
                         nrow = 1, labels = c('A', 'B'), 
                         #align = 'h',
                         rel_widths = c(1,1))

    figure2

![](../figures/02_volcanoplots/figure2old-1.png)

    pdf("../figures/figure2old.pdf", width=5.5, height=3)
    print(figure2)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    ggsave(
      "../figures/figure2.png",
      figure2,
      width = 6,
      height = 3,
      dpi = 1200
    )

Treatement volcano for figure 2 with heatmap

    volcanoplot <- ggplot(dissociation, aes(x = lfc, y = pvalue)) + 
      geom_point(aes(color = factor(direction), shape = factor(direction)), 
                 size = 1.5, alpha = 0.8, na.rm = T) + 
      theme_cowplot(font_size = 12, line_size = 0.25) +
      theme(legend.title = element_blank(),
            legend.position=c(.7,.75),
            panel.grid.minor=element_blank(),
            panel.grid.major=element_blank()) +
      scale_color_manual(values = dissociationcolor) +
      scale_x_continuous(name="log fold change") +
      scale_y_continuous(name="-log10 p-value") +
      geom_hline(yintercept = 1,  size = 0.25, linetype = 2 ) + 
      geom_hline(yintercept = 2,  size = 0.25, linetype = 2 ) + 
      scale_shape_manual(values = c(1,16,16)) +
        annotate("text", label = "4", x = -10, y = 2.2, size = 3, color = "black") + 
        annotate("text", label = "65", x = 10, y = 2.2, size = 3, color = "black") +
        annotate("text", label = "56", x = -10, y = 1.2, size = 3, color = "black") + 
        annotate("text", label = "288", x = 10, y = 1.2, size = 3, color = "black")

    volcanoplot  

![](../figures/02_volcanoplots/plot3-1.png)

    pdf(file = "../figures/02_volcanoplots/Treatment_volcano.pdf", width=3.5, height=3.25)
    plot(volcanoplot)
    dev.off()

    ## quartz_off_screen 
    ##                 2

Useful R tutorials
------------------

-   [ggplot axis
    help](http://ggplot2.tidyverse.org/reference/scale_continuous.html)
-   [grepply
    help](http://www.gettinggeneticsdone.com/2016/01/repel-overlapping-text-labels-in-ggplot2.html)
