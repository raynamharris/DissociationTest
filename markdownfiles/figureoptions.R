## this R script loads some themes for ggplots

## should be set as 'ann_colors' to be called in all pheatmapts
ann_colorsdissociation = list(Method = c(control = (values=c("#969696")), 
                              dissociated = (values=c("#525252"))),
                  Region = c(CA1 = (values=c("#7570b3")),
                            CA3 = (values=c("#1b9e77")), 
                            DG = (values=c("#d95f02"))))

ann_colorsstress = list(Group = c(stressed = (values=c("#969696")), 
                                    homecage = (values=c("#525252"))),
                  Region = c(CA1 = (values=c("#7570b3")),
                            CA3 = (values=c("#1b9e77")), 
                            DG = (values=c("#d95f02"))))

ann_colorscembrowksi = list(location = c(v = (values=c("#969696")), 
                               d = (values=c("#525252"))),
                  region =  c(ca1 = (values=c("#08306b")),  
                              ca3 = (values=c("#238443")), 
                              dg = (values=c("#a50f15"))))

ann_colorscombo = list(exp = c(cembrowski = (values=c("#969696")), 
                               harris = (values=c("#525252"))),
                  region =  c(dg = (values=c("#d95f02")),  
                              ca3 = (values=c("#1b9e77")), 
                              ca1 = (values=c("#7570b3"))))


# for the heatmap palette. call with 'colorpalette'
matlabcolors <-  matlab.like2(100)  #color schelem
cembrowskicolors <-  colorRampPalette(c("Deep Sky Blue 3", "white", "red"))( 30 )



## for the pca plot
colorvalRegion <- c("#7570b3", "#1b9e77", "#d95f02")
colorvalMethod <- c("#969696", "#525252")
colorvalGroup <- c("#969696", "#525252")

# color blind safe red blue green -ish colors
# orange #d95f02
# purple #7570b3
# green #1b9e77
