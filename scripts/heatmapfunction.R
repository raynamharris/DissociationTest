# two heatmap functions to create custom pdf and a png
# uses the function pheatmap https://github.com/raivokolde/pheatmap

# usage: name of matrix, color palette, lengend, title
heatmap_png <- function(DEGs, ann_colors, df, title, clustercolsmethod){

  DEGs <- DEGs[order(DEGs$padjmin),]
  DEGs <- head(DEGs, 30)
  print(head(DEGs, 5))
  
  rownames(DEGs) <- DEGs$rownames
  drop.cols <-colnames(DEGs[,grep("padj|pval|rownames", colnames(DEGs))])
  DEGs <- DEGs %>% select(-one_of(drop.cols))
  DEGs <- as.matrix(DEGs)
  DEGs <- DEGs - rowMeans(DEGs)
  
  paletteLength <- 30
  myBreaks <- c(seq(min(DEGs), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(DEGs)/paletteLength, max(DEGs), length.out=floor(paletteLength/2)))
  
  pheatmap(DEGs, show_colnames=F, show_rownames = T,
           annotation_col=df, annotation_colors = ann_colors, 
           annotation_row = NA, annotation_legend = TRUE,
           annotation_names_row = FALSE, annotation_names_col = TRUE,
           treeheight_row = 0, treeheight_col = 15,
           fontsize = 7,
           border_color = "grey60" ,
           color = viridis(30),
           clustering_method="median",
           breaks=myBreaks,
           clustering_distance_cols=clustercolsmethod, 
           cluster_cols = T,
           main = title)  
}


heatmap_pdf <- function(DEGs, ann_colors, df, title, clustercolsmethod){
  
  myfile <-  paste("../figures/03_heatmaps/", substitute(DEGs), substitute(clustercolsmethod), ".pdf", sep="")

  DEGs <- DEGs[order(DEGs$padjmin),]
  DEGs <- head(DEGs, 30)

  rownames(DEGs) <- DEGs$rownames
  drop.cols <-colnames(DEGs[,grep("padj|pval|rownames", colnames(DEGs))])
  DEGs <- DEGs %>% select(-one_of(drop.cols))
  DEGs <- as.matrix(DEGs)
  DEGs <- DEGs - rowMeans(DEGs)
  
  paletteLength <- 30
  myBreaks <- c(seq(min(DEGs), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(DEGs)/paletteLength, max(DEGs), length.out=floor(paletteLength/2)))
  
  pheatmap(DEGs, show_colnames=F, show_rownames = T,
           annotation_col=df, annotation_colors = ann_colors, 
           annotation_row = NA, annotation_legend = TRUE,
           annotation_names_row = FALSE, annotation_names_col = TRUE,
           treeheight_row = 0, treeheight_col = 5,
           fontsize = 7, 
           border_color = "grey60" ,
           color = viridis(30),
           width=3.5, height=3.5,
           clustering_method="median",
           breaks=myBreaks,
           clustering_distance_cols=clustercolsmethod, 
           cluster_cols = T,
           main = title,
           filename =  myfile)
}