
## for the pca plot
colorvalRegion <- c("#7570b3", "#1b9e77", "#d95f02")
colorvalTreatment <- c("#ffffff", "#525252")

## pheatmapt useage: set 'ann_colors' to be one of these
ann_colorsdissociation = list(Treatment = c(control = (values=c("#ffffff")), 
                              dissociated = (values=c("#525252"))),
                  Region = c(CA1 = (values=c("#7570b3")),
                            CA3 = (values=c("#1b9e77")), 
                            DG = (values=c("#d95f02"))))

ann_colorsstress = list(Treatment = c(homecage = (values=c("#ffffff")),
                                  shocked = (values=c("#525252"))),
                  Region = c(CA1 = (values=c("#7570b3")),
                            CA3 = (values=c("#1b9e77")), 
                            DG = (values=c("#d95f02"))))

ann_colorsbehavior = list(Treatment = c(yoked = (values=c("#ffffff")),
                                  trained = (values=c("#525252"))),
                        Region = c(CA1 = (values=c("#7570b3")),
                                   CA3 = (values=c("#1b9e77")), 
                                   DG = (values=c("#d95f02"))))

ann_colorscembrowksi = list(Location = c(dorsal = (values=c("#ffffff")),
                                         ventral = (values=c("#525252"))),
                  Region =  c(CA1 = (values=c("#7570b3")),  
                              CA3 = (values=c("#1b9e77")), 
                              DG = (values=c("#d95f02"))))


# for the pheatmap palette. call with 'colorpalette'
cembrowskicolors <-  colorRampPalette(c("Deep Sky Blue 3", "white", "red"))( 30 )


volcano1 <-  c("dissociated" = "#525252",
               "control" = "#525252",
               "none" = "#f0f0f0")

