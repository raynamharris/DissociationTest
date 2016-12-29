## this R script loads some themes for R

## should be set as 'ann_colors' to be called in all pheatmapts
ann_colorsdissociation = list(Method = c(Homogenized = (values=c("#969696")), 
                              Dissociated = (values=c("#525252"))),
                  Punch =  c(DG = (values=c("#a50f15")),  
                             CA3 = (values=c("#238443")),
                             CA1 = (values=c("#08306b"))))

ann_colorsbehavior = list(Method = c(Yoked = (values=c("#f1a340")), 
                                    Trained = (values=c("#9970ab"))),
                  Punch = c(CA1 = (values=c("#08306b")),
                            CA3 = (values=c("#238443")), 
                            DG = (values=c("#a50f15"))))

ann_colorscembrowksi = list(location = c(v = (values=c("#000000")), 
                               d = (values=c("#525252"))),
                  region =  c(ca1 = (values=c("#08306b")),  
                              ca3 = (values=c("#238443")), 
                              dg = (values=c("#a50f15"))))

ann_colorscombo = list(exp = c(cembrowski = (values=c("#969696")), 
                               harris = (values=c("#525252"))),
                  region =  c(dg = (values=c("#980043")),  
                              ca3 = (values=c("#c994c7")), 
                              ca1 = (values=c("#dd1c77"))))



## for the heatmap pallet. calle with 'colorpalette'
matlabcolors <-  matlab.like2(100)  #color schelem
cembrowskicolors <-  colorRampPalette(c("Deep Sky Blue 3", "white", "red"))( 30 )
