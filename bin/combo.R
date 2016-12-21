## first read in cembrowski from Cembrowski.Rmd
cembrowskicountdata <- countData
cembrowskicountdata$rownames <- row.names(cembrowskicountdata)
cembrowskicoldata <- colData
cembrowskicoldata$exp <- rep("cembrowski",nrow(cembrowskicoldata))
head(cembrowskicountdata)

cembrowskicountdata$rownames <- NULL
normlizedcembrowskicountdata <- cembrowskicountdata/28
normlizedcembrowskicountdata$rownames <- row.names(normlizedcembrowskicountdata)

## then read in data from DissociationTest.Rmd
library(plyr)
raynacountdata <- countData
raynacountdata$rownames <- row.names(raynacountdata)
raynacoldata <- colData
raynacoldata$exp <- rep("harris",nrow(raynacoldata))
raynacoldata$location <- rep("d",nrow(raynacoldata))
raynacoldata <- raynacoldata %>% select(RNAseqID, Punch, location, exp)
names(raynacoldata)[2] <- "region"
raynacoldata$region <- revalue(raynacoldata$region, c("CA1" = "ca1")) 
raynacoldata$region <- revalue(raynacoldata$region, c("CA3" = "ca3")) 
raynacoldata$region <- revalue(raynacoldata$region, c("DG" = "dg")) 

raynacountdata$rownames <- NULL
normlizedraynacountdata <- raynacountdata*10
normlizedraynacountdata$rownames <- row.names(normlizedraynacountdata)


allcolData <- rbind(cembrowskicoldata,raynacoldata)
allcountData <- inner_join(normlizedcembrowskicountdata,normlizedraynacountdata)
row.names(allcountData) <- allcountData$rownames
allcountData$rownames <- NULL


dds <- DESeqDataSetFromMatrix(countData = allcountData,
                              colData = allcolData,
                              design = ~ region + exp + region * exp )
## filter genes with 0 counts
dds <- dds[ rowSums(counts(dds)) > 2, ]
dds

# Differential expression analysis
dds <- DESeq(dds)
dds

# general deseq
res <- results(dds, independentFiltering = F)
#res
summary(res)
resOrdered <- res[order(res$padj),]
sum(res$padj < 0.1, na.rm = TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

## make a plot of fold change as function of expression

plot <- plotMA(res05, main="MA Plot")

plotCounts(dds, gene=which.min(res$padj), intgroup="region")
plotCounts(dds, gene=which.min(res$padj), intgroup="exp")

rld <- rlog(dds, blind=FALSE)

resdorsalventral <- results(dds, contrast = c("exp", "cembrowski", "harris"), independentFiltering = F)
sum(resdorsalventral$padj < 0.1, na.rm = TRUE) #7487
valsdorsalventral <- cbind(resdorsalventral$pvalue, resdorsalventral$padj)
colnames(valsdorsalventral)=c("pval.dorsalventral", "padj.dorsalventral")

resCA1DG <- results(dds, contrast = c("region", "ca1", "dg"), independentFiltering = F)
sum(resCA1DG$padj < 0.1, na.rm = TRUE) #4492
valsCA1DG <- cbind(resCA1DG$pvalue, resCA1DG$padj) 
colnames(valsCA1DG)=c("pval.CA1DG", "padj.CA1DG")

resCA1CA3 <- results(dds, contrast = c("region", "ca1", "ca3"), independentFiltering = F)
sum(resCA1CA3$padj < 0.1, na.rm = TRUE) #1732
valsCA1CA3 <- cbind(resCA1CA3$pvalue, resCA1CA3$padj) 
colnames(valsCA1CA3)=c("pval.CA1CA3", "padj.CA1CA3")

resCA3DG <- results(dds, contrast = c("region", "ca3", "dg"), independentFiltering = F)
sum(resCA3DG$padj < 0.1, na.rm = TRUE) #4659
valsCA3DG <- cbind(resCA3DG$pvalue, resCA3DG$padj) 
colnames(valsCA3DG)=c("pval.CA3DG", "padj.CA3DG")

rldd <- assay(rld)
rldpvals <- cbind(rldd, valsdorsalventral,valsCA1DG, valsCA1CA3, valsCA3DG)

rldpvals <- as.data.frame(rldpvals)

dorsalventral <- row.names(rldpvals[rldpvals$padj.dorsalventral<0.1 & !is.na(rldpvals$padj.dorsalventral),])
CA1DG <- row.names(rldpvals[rldpvals$padj.CA1DG<0.1 & !is.na(rldpvals$padj.CA1DG),])
CA1CA3 <- row.names(rldpvals[rldpvals$padj.CA1CA3<0.1 & !is.na(rldpvals$padj.CA1CA3),])
CA3DG <- row.names(rldpvals[rldpvals$padj.CA3DG<0.1 & !is.na(rldpvals$padj.CA3DG),])

## four way grid
candidates <- list("CA1 v. DG" = CA1DG, "CA1 v. CA3" = CA1CA3, "CA3 v. DG" = CA3DG,  "Cembrowski v. Harris" = dorsalventral )
dev.off()
prettyvenn <- venn.diagram(
  x = candidates, filename=NULL, lwd=4,
  col = "transparent",
  fill = (values=c("#00441b", "#00441b","#238b45", "#238b45")),
  alpha = 0.5,
  cex = 1, fontfamily = "sans", #fontface = "bold",
  cat.default.pos = "text",
  #cat.col = c("darkred", "darkgreen", "blue4", "orange"),
  cat.dist = c(0.08, 0.08, 0.08, 0.08), cat.pos = 1,
  cat.cex = 1, cat.fontfamily = "sans")
grid.draw(prettyvenn)




nt <- normTransform(dds) # defaults to log2(x+1) 
df <- as.data.frame(colData(dds)[,c("region", "exp")])

ann_colors = list(exp = c(cembrowski = (values=c("#969696")), harris = (values=c("#525252"))),
                  region =  c(dg = (values=c("#980043")),  ca3 = (values=c("#c994c7")), 
                              ca1 = (values=c("#dd1c77"))))
matlabcolors <-  matlab.like2(100)  #color scheme

DEGes <- as.data.frame(rldpvals) # convert matrix to dataframe
DEGes$rownames <- rownames(DEGes)  # add the rownames to the dataframe
DEGes$padjmin <- with(DEGes, pmin(padj.dorsalventral)) # put the min pvalue in a new column
DEGes <- DEGes %>% filter(padjmin < 0.00000000000000001)
rownames(DEGes) <- DEGes$rownames
drop.cols <- c("padj.dorsalventral", "pval.dorsalventral", "padj.CA1DG" , "padj.CA1CA3" , "padj.CA3DG" , "pval.CA1DG" , "pval.CA1CA3" , "pval.CA3DG" , "rownames", "padjmin")
DEGes <- DEGes %>% select(-one_of(drop.cols))
DEGes <- as.matrix(DEGes)
DEGes <- DEGes - rowMeans(DEGes)

pheatmap(DEGes, show_colnames=T, show_rownames = F,
         annotation_col=df, annotation_colors = ann_colors,
         fontsize = 12, fontsize_row = 7, 
         #cellwidth=10, cellheight=10, width = 10,
         border_color = "grey60" ,
         color = matlabcolors,
         main = "top DE genes p <<<<< 0.01"
)

DEGes <- as.data.frame(rldpvals) # convert matrix to dataframe
DEGes$rownames <- rownames(DEGes)  # add the rownames to the dataframe
DEGes$padjmin <- with(DEGes, pmin(padj.CA1DG)) # put the min pvalue in a new column
DEGes <- DEGes %>% filter(padjmin < 0.000001)
rownames(DEGes) <- DEGes$rownames
drop.cols <- c("padj.dorsalventral", "pval.dorsalventral", "padj.CA1DG" , "padj.CA1CA3" , "padj.CA3DG" , "pval.CA1DG" , "pval.CA1CA3" , "pval.CA3DG" , "rownames", "padjmin")
DEGes <- DEGes %>% select(-one_of(drop.cols))
DEGes <- as.matrix(DEGes)
DEGes <- DEGes - rowMeans(DEGes)

pheatmap(DEGes, show_colnames=T, show_rownames = F,
         annotation_col=df, annotation_colors = ann_colors,
         fontsize = 12, fontsize_row = 7, 
         #cellwidth=10, cellheight=10, width = 10,
         border_color = "grey60" ,
         color = matlabcolors,
         main = "top DE genes CA1DG p < 0.000001"
)


plotPCA(rld, intgroup=c("region", "exp"), returnData=TRUE)
pcadata <- plotPCA(rld, intgroup=c("region", "exp"), returnData=TRUE)
percentVar <- round(100 * attr(pcadata, "percentVar"))

ggplot(pcadata, aes(PC1, PC2, color=region, shape=exp, label=name)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_classic() +
  geom_text(aes(label=name),vjust=2)


counts <- allcountData
dim( counts )
colSums( counts ) / 1e06  # in millions of reads
table( rowSums( counts ) )[ 1:30 ] # Number of genes with low counts

rowsum <- as.data.frame(colSums( counts ) / 1e06 )
names(rowsum)[1] <- "millioncounts"
rowsum$sample <- row.names(rowsum)

ggplot(rowsum, aes(x=millioncounts)) + 
  geom_histogram(bins = 20, colour = "black", fill = "darkgrey") +
  theme_classic() +
  scale_x_continuous(name = "Millions of Gene Counts per Sample",
                     breaks = seq(0, 30, 5),
                     limits=c(0, 30)) +
  scale_y_continuous(name = "Number of Samples")
