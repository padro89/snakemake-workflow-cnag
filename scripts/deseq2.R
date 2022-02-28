# Guardar el workspace per debugging amb snakemake
# Netejar el workspace abans
eliminats <- ls()
eliminats <- eliminats[eliminats!="snakemake" & eliminats!="Snakemake"]
rm(list=eliminats)
rm(eliminats)
# Guardar-lo amb només l'objecte snakemake actualitzat
save.image(file="workspace",)


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
group <- "GROUP"
control <- "Unaffected"
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

# Create a results directory
if(!file.exists(file.path(directory,"Results"))){
  dir.create(file.path(directory,"Results"))
}

## Import counts matrix and sample info ##
counts_raw <- read.table(snakemake@input[[1]],
                         header = T, row.names = 1,
                         check.names = F)
info_raw <- read.table(snakemake@input[[2]],
                       header = T, row.names = 1)
coldata=subset(info_raw, row.names(info_raw) %in% colnames(counts_raw))
countdata=counts_raw[,colnames(counts_raw) %in% rownames(coldata)]
coldata=subset(coldata, row.names(coldata) %in% colnames(countdata))
names(coldata)[which(names(coldata)==group)] = "group"

## RELEVEL CONTROL GROUP ##
coldata$group <- relevel(factor(coldata$group), ref=control)
#dds[["group"]] <- relevel(dds[["group"]], control)

##READ FILES ##

# Adapt info
if (!is.null(factors)){
  for (i in seq(length(factors))){
    coldata[,factors[i]]=as.factor(coldata[,factors[i]])
  }
}

## Funció prèvia. Probablement s'hagi d'adaptar l'altra.
#if (!is.null(continuous)){
#  for (i in seq(length(names(continuous)))){
#    levels_tmp <- paste(names(continuous)[i],seq(as.numeric(continuous[i][2])), sep="")
#    coldata[,continuous[i][1]] = cut(coldata[,continuous[i][1]],
#                                                as.numeric(continuous[i][2]),
#                                                labels=levels_tmp)
#  }
#}

#if (!is.null(continuous)){
#  for (i in seq(length(continuous))){
#    coldata[,continuous[i]] <- as.numeric(coldata[,continuous[i]])
#  }
#  
#}

#Sort names countdata as rownames in coldata
countdata <- countdata[rownames(coldata)]
## BATCH EFFECT ##
dds <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = formula(eval(parse(text=formula))))

dds <- estimateSizeFactors(dds)
bc_per_group <- lapply(unique(coldata$group), function (x) rownames(coldata)[coldata$group==x])
all_members <- function(x){
  if (length(x)>1){
    return(rowSums(counts(dds,normalized=TRUE)[,unlist(x)] >= 0)==length(unlist(x)))
  }else{
    return(counts(dds,normalized=TRUE)[,unlist(x)] >= 0)
  }
}
keep <- Reduce("|", lapply(bc_per_group, all_members))
#keep <- rowSums(counts(dds,normalized=TRUE) >= 10) >= min(rle(as.vector(coldata$group))$lengths)
dds <- dds[keep,] 

## DIFFERENTIAL ANALYSIS ##
if (onlypca == F){
  dds <- DESeq(dds, parallel=TRUE)  
}



############## POLTS SENSE CONTRASTOS ############################################################
## rlog transformation and variance stabilization ##
#other normalization methods for plotting#

rld <- rlog(dds)
#vsd <-varianceStabilizingTransformation(dds)
rlogMat<-assay(rld)
write.table(rlogMat, file = paste(file.path(directory,"Results/"),project,"_rlogMat.txt",sep=""), quote = F)
#vstMat<-assay(vsd)

#heatmap samples"
if (plot == T){
  sampleDists <- dist(t(rlogMat))
  sampleDistMatrix <- as.matrix(sampleDists)
  #colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pdf(paste(file.path(directory,"Results/"),project,"_sampletosample_heatmap.pdf",sep=""),onefile=FALSE)
 #annotation_heatmap = as.data.frame(sapply(plot_atr$heatmap_ann,function(x) eval(parse(text=x))), row.names = colnames(coldata))
  annotation_heatmap=as.data.frame(colData(dds)[,heatmap_atr], row.names=row.names(colData(dds)))
  names(annotation_heatmap) = heatmap_atr
## SAVED: annotation_col = as.data.frame(coldata[,plot_atr$heatmap_ann]
  pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, 
           fontsize = 4, border_color = NA, annotation_col=annotation_heatmap, row.names = row.names(coldata))
  dev.off()
  print("sample correlation heatmap done")
}

## PCA!!!
#Select

rv <- rowVars(assay(rld))
if (plot == T){
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
pca <- prcomp(t(assay(rld)[select,]))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
dims <- combn(c(1:4),2)
tmp_color = coldata[,pca_atr]
if (length(pca_atr)>1){
  color_factor = factor(apply(tmp_color, 1 ,paste, collapse=" : "))  
}else{
  color_factor = tmp_color
}

#  tiff(paste(args$project,"_pca_custom.tiff", sep=""), width=2000, height=2000, res=300)
#  print(ggplot()+geom_text(aes_string(pca$x[,dims[1,1]], pca$x[,dims[2,1]], color=color_factor, label = as.factor(coldata$CODE)))+
#        xlab(paste(colnames(pca$x)[1]," (",round(percentVar[1]*100,digits = 2),"%)",sep=''))+
#        ylab(paste(colnames(pca$x)[2]," (",round(percentVar[2]*100,digits=2),"%)",sep=''))+
#        theme(legend.title = element_blank()))
#  dev.off()

  pdf(paste(file.path(directory,"Results/"),project,"_pca.pdf", sep=""))
  print(ggplot()+geom_text(aes_string(x=pca$x[,dims[1,1]], y=pca$x[,dims[2,1]], color = color_factor, label = as.factor(rownames(pca$x))))+
    xlab(paste(colnames(pca$x)[1]," (",round(percentVar[1]*100,digits = 2),"%)",sep=''))+
    ylab(paste(colnames(pca$x)[2]," (",round(percentVar[2]*100,digits=2),"%)",sep=''))+
    theme(legend.title = element_blank()))


  grid_arrange_shared_legend <- function(plots) {
  #  plots <- list(...)
    g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    grid.arrange(
      do.call(arrangeGrob, lapply(plots, function(x)
        x + theme(legend.position="none"))),
      legend,
      ncol = 1,
      heights = unit.c(unit(1, "npc") - lheight, lheight))
  }
  
  toplot <- lapply(c(2:dim(dims)[2]), function(i) ggplot()+geom_point(aes_string(x=pca$x[,dims[1,i]], y=pca$x[,dims[2,i]], color = color_factor))+
                     xlab(paste(colnames(pca$x)[dims[1,i]]," (",round(percentVar[dims[1,i]]*100,digits = 2),"%)",sep=''))+
                     ylab(paste(colnames(pca$x)[dims[2,i]]," (",round(percentVar[dims[2,i]]*100,digits=2),"%)",sep=''))+
                     theme(legend.title = element_blank(), legend.position = 'none'))
  
  #grid.arrange(grobs = toplot, nrow=3)
  grid_arrange_shared_legend(toplot)
  
  dims <- combn(c(5:ifelse(dim(pca$x)[2]<8, dim(pca$x)[2],8)),2)
  toplot <- lapply(c(1:dim(dims)[2]), function(i) ggplot()+geom_point(aes_string(x=pca$x[,dims[1,i]], y=pca$x[,dims[2,i]], color = color_factor))+
                     xlab(paste(colnames(pca$x)[dims[1,i]]," (",round(percentVar[dims[1,i]]*100,digits = 2),"%)",sep=''))+
                     ylab(paste(colnames(pca$x)[dims[2,i]]," (",round(percentVar[dims[2,i]]*100,digits=2),"%)",sep=''))+
                     theme(legend.title = element_blank(),  legend.position = 'none'))
  
  
  #grid.arrange(grobs = toplot, nrow=3)
  grid_arrange_shared_legend(toplot)
  
  dev.off()

write.table( sweep(abs(pca$rotation), 2, colSums(abs(pca$rotation)), "/"), paste(file.path(directory,"Results/"),project,'_pc_contribution.txt', sep=""), 
             sep = "\t", quote = F,row.names = T, col.names = T)

}


################## CONTRASTS FROM NOW ###############################################


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
  sink(paste(file.path(directory,"Results/"),project,'_',title,"_stats.txt",sep=""))
  mcols(res)$description
  summary(res)
  sink()
  
  ## WRITE TABLE RESULTS ##
  write.table(cc, paste(file.path(directory,"Results/"),project,"_",title,"_norm_counts.txt",sep=""),quote=FALSE)
  pass_filter <- as.numeric(as.numeric(rownames(resAllOrdered) %in% rownames(resOrdered)) & (resAllOrdered$padj<0.05) )
  df_all <- as.data.frame(resAllOrdered)
  df_all["filter"]<- pass_filter
  df_all["shrunkenlfc"] = res2[rownames(df_all),"log2FoldChange"]
  df_all <- df_all[,c("baseMean","log2FoldChange",#"shrunkenlfc",
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

if (onlypca == F){
  sapply(names(contrastos), function(x) process_contrast(x, contrastos[[x]]))  
}
 




