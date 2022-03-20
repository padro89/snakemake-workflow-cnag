# Guardar el workspace per debugging amb snakemake
# Netejar el workspace abans
#eliminats <- ls()
#eliminats <- eliminats[eliminats!="snakemake" & eliminats!="Snakemake"]
#rm(list=eliminats)
#rm(eliminats)
# Guardar-lo amb només l'objecte snakemake actualitzat
save.image(file="workspace")

# Load required libraries
library("DESeq2")
library("BiocParallel")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")
library("gridExtra")
require("grid")
library("methods")

## Import parameters from the config file ##
# S'utilitzen
plot <- snakemake@config$plot
shrinkage_method <-snakemake@config$shrinkage_method
heatmap_atr <- snakemake@config$plot_atr$heatmap_ann
de_genes_n <- snakemake@config$plot_atr$de_genes_n

## Load the data

load(snakemake@input$dds_design)
rlogMat <- as.matrix(read.table(snakemake@input$rlogmat,
                                header = T, row.names = 1))
contrast <- c("group",snakemake@params$contrast)


## EXTRACT RESULTS ##
resAll <- results(dds, cooksCutoff=TRUE, contrast=contrast, parallel=TRUE)
res2 <- lfcShrink(dds, contrast=contrast, res=resAll, type=shrinkage_method)
res <- subset(res2, abs(log2FoldChange) > log2(1.5))
resOrdered<- res[order(res$padj),]
resAllOrdered <- resAll[order(resAll$padj),]
  
## EXTRACT DESCRIPTION AND SUMMARY ##
description_stats<-mcols(res)$description
summary_stats<-summary(res)
stats<-rbind(c(description_stats, summary_stats))
sink(snakemake@output$stats)
mcols(res)$description
summary(res)
sink()
  
## WRITE TABLE RESULTS ##
pass_filter <- as.numeric(as.numeric(rownames(resAllOrdered) %in% rownames(resOrdered)) & (resAllOrdered$padj<0.05) )
df_all <- as.data.frame(resAllOrdered)
df_all["filter"]<- pass_filter
df_all["shrunkenlfc"] = res2[rownames(df_all),"log2FoldChange"]
df_all <- df_all[,c("baseMean","log2FoldChange","shrunkenlfc",
                    "lfcSE","stat", "filter", "pvalue", "padj")]
write.table(df_all,file = snakemake@output$deg_results, quote=FALSE)
print("DE analysis finished")

if ((plot == T) & (length(rownames(resOrdered))>=de_genes_n)){  
  select=rlogMat[rownames(resOrdered),][1:de_genes_n,]
  data_subset=subset(colData(dds), dds[[contrast[1]]] %in% c(contrast[3], contrast[2]))
  select2=select[,row.names(data_subset)]
  
  df=as.data.frame(data_subset[,heatmap_atr], row.names=row.names(data_subset))
  colnames(df)=c(heatmap_atr)
    
  pdf(file = snakemake@output$topDEgenes_heatmap,onefile=FALSE)
    
  pheatmap(select2, fontsize_row=6,show_rownames=TRUE, cluster_rows=TRUE, cluster_cols=TRUE, annotation_col=df,scale="row", fontsize=4, show_colnames = TRUE)
  dev.off()
  print(paste("top", de_genes_n, "DE genes heatmap done"))

  tiff(snakemake@output$heatmap_custom, width=2000, height=2000, res=300)
  pheatmap(select2, fontsize_row=6,show_rownames=TRUE, cluster_rows=TRUE, cluster_cols=FALSE, annotation_col=df,scale="row", fontsize=4, show_colnames = TRUE)
  dev.off()

} else {
  if ((plot == T) & (length(rownames(resOrdered))>= 25)){
    select=rlogMat[rownames(resOrdered),][1:length(rownames(resOrdered)),]
    data_subset=subset(colData(dds), dds[[contrast[1]]] %in% c(contrast[3], contrast[2]))
    select2=select[,row.names(data_subset)]
    
    df=as.data.frame(data_subset[,heatmap_atr], row.names=row.names(data_subset))
    colnames(df)=c(heatmap_atr)
      
    pdf(snakemake@output$topDEgenes_heatmap,onefile=FALSE)
    
    pheatmap(select2, fontsize_row=6,show_rownames=TRUE, cluster_rows=TRUE, cluster_cols=TRUE, annotation_col=df,scale="row", fontsize=4, show_colnames = TRUE)
    dev.off()
    print(paste("top", de_genes_n, "DE genes heatmap done"))
      
    tiff(snakemake@output$heatmap_custom, width=2000, height=2000, res=300)
    pheatmap(select2, fontsize_row=6,show_rownames=TRUE, cluster_rows=TRUE, cluster_cols=FALSE, annotation_col=df,scale="row", fontsize=4, show_colnames = TRUE)
    dev.off()
  }
}




