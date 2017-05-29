## this R script loads some themes for ggplots

# color blind safe red blue green -ish colors
# orange #d95f02 - DG
# purple #7570b3 - CA1
# green #1b9e77  - CA3

## for the pca plot
colorvalRegion <- c("#7570b3", "#1b9e77", "#d95f02")
colorvalTreatment <- c("#969696", "#525252")

## pheatmapt useage: set 'ann_colors' to be one of these
ann_colorsdissociation = list(Treatment = c(control = (values=c("#969696")), 
                              dissociated = (values=c("#525252"))),
                  Region = c(CA1 = (values=c("#7570b3")),
                            CA3 = (values=c("#1b9e77")), 
                            DG = (values=c("#d95f02"))))

ann_colorsstress = list(Treatment = c(homecage = (values=c("#969696")),
                                  shocked = (values=c("#525252"))),
                  Region = c(CA1 = (values=c("#7570b3")),
                            CA3 = (values=c("#1b9e77")), 
                            DG = (values=c("#d95f02"))))

ann_colorsbehavior = list(Treatment = c(yoked = (values=c("#969696")),
                                  trained = (values=c("#525252"))),
                        Region = c(CA1 = (values=c("#7570b3")),
                                   CA3 = (values=c("#1b9e77")), 
                                   DG = (values=c("#d95f02"))))

ann_colorscembrowksi = list(Location = c(dorsal = (values=c("#525252")),
                                         ventral = (values=c("#969696"))),
                  Region =  c(CA1 = (values=c("#7570b3")),  
                              CA3 = (values=c("#1b9e77")), 
                              DG = (values=c("#d95f02"))))


# for the pheatmap palette. call with 'colorpalette'
cembrowskicolors <-  colorRampPalette(c("Deep Sky Blue 3", "white", "red"))( 30 )