# prep contrasts
prepcontrasttable <- function(factor, level1, level2, resname, valsname){
  results(dds, contrast = c(factor, level1, level2), independentFiltering = F)
  sum(resname$padj < 0.1, na.rm = TRUE)
  valsname <- cbind(resname$pvalue, resname$padj)
  colnames(valsname)=c("pval.factor", "padj.factor")
  return(valsname)
}

