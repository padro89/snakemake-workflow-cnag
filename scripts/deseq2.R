# Guardar el workspace per debugging amb snakemake
# Netejar el workspace abans
#eliminats <- ls()
#eliminats <- eliminats[eliminats!="snakemake" & eliminats!="Snakemake"]
#rm(list=eliminats)
#rm(eliminats)
# Guardar-lo amb només l'objecte snakemake actualitzat
# save.image(file="workspace",)


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
# S'han de crear wildcards
group <- snakemake@config$group
control <- snakemake@config$control
project <- snakemake@config$project
plot <- snakemake@config$plot
onlypca <- snakemake@config$onlypca
factors <- snakemake@config$factors
continuous <- snakemake@config$continuous
formula <- snakemake@config$formula
contrastos <- snakemake@config$contrast
shrinkage_method <-snakemake@config$shrinkage_method
pca_atr <- snakemake@config$plot_atr$pca
heatmap_atr <- snakemake@config$plot_atr$heatmap_ann
de_genes_n <- snakemake@config$plot_atr$de_genes_n
directory <- snakemake@config$directory


################## CONTRASTS FROM NOW ###############################################

rlogMat <- as.matrix(read.table(snakemake@input$rlogmat,
                                header = T, row.names = 1))
load(snakemake@input$dds)
dds@design <- formula(snakemake@config$formula)
dds <- DESeq(dds, parallel=TRUE)

## EXTRACT RESULTS ##
process_contrast <- function(title, contrastos){
  resAll <- results(dds, cooksCutoff=TRUE, contrast=contrastos, parallel=TRUE)
  res2 <- lfcShrink(dds, contrast=contrastos , res=resAll, type=shrinkage_method)
  res <- subset(res2, abs(log2FoldChange) > log2(1.5))
  resOrdered<- res[order(res$padj),]
  resAllOrdered <- resAll[order(resAll$padj),]
  ## EXTRACT COUNTS NORMALIZED ##
  c=counts(dds,normalized=TRUE)
  c_ordered=c[rownames(resAllOrdered),]
  colnames(c_ordered)=paste(dds[[contrastos[1]]],",",colnames(c_ordered),sep="")
  counts_dds=(c_ordered[,order(colnames(c_ordered))])
  cc=round(counts_dds,digits=2)
  
  ## EXTRACT DESCRIPTION AND SUMMARY ##
  description_stats<-mcols(res)$description
  summary_stats<-summary(res)
  stats<-rbind(c(description_stats, summary_stats))
  sink(snakemake@output$stats)
  mcols(res)$description
  summary(res)
  sink()
  
  ## WRITE TABLE RESULTS ##
  write.table(cc, paste(file.path(directory,"Results/"),project,"_",title,"_norm_counts.txt",sep=""),quote=FALSE)
  pass_filter <- as.numeric(as.numeric(rownames(resAllOrdered) %in% rownames(resOrdered)) & (resAllOrdered$padj<0.05) )
  df_all <- as.data.frame(resAllOrdered)
  df_all["filter"]<- pass_filter
  df_all["shrunkenlfc"] = res2[rownames(df_all),"log2FoldChange"]
  df_all <- df_all[,c("baseMean","log2FoldChange","shrunkenlfc",
                      "lfcSE","stat", "filter", "pvalue", "padj")]
  write.table(df_all,paste(file.path(directory,"Results/"),project,"_", title,"_results.txt",sep=""),quote=FALSE)
  print("DE analysis finished")

  if ((plot == T) & (length(rownames(resOrdered))>=de_genes_n)){  
    select=rlogMat[rownames(resOrdered),][1:de_genes_n,]
    data_subset=subset(colData(dds), dds[[contrastos[1]]] %in% c(contrastos[3], contrastos[2]))
    select2=select[,row.names(data_subset)]
    
    df=as.data.frame(data_subset[,heatmap_atr], row.names=row.names(data_subset))
    colnames(df)=c(heatmap_atr)
    
    pdf(paste(file.path(directory,"Results/"),project,"_",title,"_top50DEgenes_heatmap.pdf",sep=""),onefile=FALSE)
    
    pheatmap(select2, fontsize_row=6,show_rownames=TRUE, cluster_rows=TRUE, cluster_cols=TRUE, annotation_col=df,scale="row", fontsize=4, show_colnames = FALSE)
    dev.off()
    print("top 50 DE genes heatmap done")

    tiff(paste(file.path(directory,"Results/"),project,"_heatmap_custom.tiff", sep=""), width=2000, height=2000, res=300)
    pheatmap(select2, fontsize_row=6,show_rownames=TRUE, cluster_rows=TRUE, cluster_cols=FALSE, annotation_col=df,scale="row", fontsize=4, show_colnames = TRUE)
    dev.off()

  }
}

sapply(names(contrastos), function(x) process_contrast(x, contrastos[[x]]))  




