#####ARGPARSER

args <- commandArgs(TRUE)
# args <- c('--counts=tmp_dir/COUNTS_genes_IBD_02','--info=tmp_dir/info_IBD_02_merged',
#          '--group=dx_UC_extension_CD_Location_location_inflammation',
#          '--control=Healthy_control_na_na_l4_c',
#          '--config=tmp_dir/deseq2.example.config.R',
#          '--project=IDB_02')
## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      deseq2.R
      
      Arguments:
      --counts=table   
      --info=table
      --group=column_name    -Column where groups are
      --contol=control_group
      --config=config.R      -List in R syntaxis
      --project=project_name
      --plot=boolean         -Set FALSE if you don't want plots (default: TRUE)
      --onlypca=boolean
      --help        
      
      Example:
      deseq2.R --counts=\"table.txt\" --info=\"output.txt\" --group=column_name --control=control_group --config=config.R --onlypca=FALSE --plot=TRUE \n\n")
  
  q(save="no")
}

parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
args <- as.list(as.character(argsDF$V2))
names(args) <- argsDF$V1

if(is.null(args$counts)) {
  cat('--counts is mandatory')
  q(save='no')
}

if(is.null(args$info)) {
  cat('--info is mandatory')
  q(save='no')
}

if(is.null(args$group)) {
  cat('--group is mandatory')
  q(save='no')
}


if(is.null(args$control)) {
  cat('--control is mandatory')
  q(save='no')
}

if(is.null(args$config)) {
  cat('--config is mandatory')
  q(save='no')
}

if(is.null(args$project)) {
  cat('--project is mandatory')
  q(save='no')
}

if(is.null(args$plot)) {
  args$plot=TRUE
}else{
  if (args$plot == "TRUE"){
    args$plot=TRUE  
  }else{
    args$plot=FALSE  
  }
}

if(is.null(args$onlypca)) {
  args$onlypca=FALSE
}else{
  if (args$onlypca == "TRUE"){
    args$onlypca=TRUE  
  }else{
    args$onlypca=FALSE  
  }
}

print(args)



##LIBRARIES ##
library("DESeq2")
library("BiocParallel")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")
library("gridExtra")
require("grid")
library("methods")

##READ FILES ##
group = args$group
control = args$control
counts_raw=read.table(args$counts,h=T,row.names=1,check.names=FALSE)
info_raw=read.table(args$info,h=T,row.names=1)
coldata=subset(info_raw, row.names(info_raw) %in% colnames(counts_raw))
countdata=counts_raw[,colnames(counts_raw) %in% rownames(coldata) ]
coldata=subset(coldata, row.names(coldata) %in% colnames(countdata))
names(coldata)[which(names(coldata)==group)] = "group"
# Adapt info
source(args$config)
if (!is.null(covariants$factor)){
  for (i in seq(length(names(covariants$factor)))){
    coldata[,covariants$factor[[i]]]=as.factor(coldata[,covariants$factor[[i]]])
  }
}

if (!is.null(covariants$continuous)){
  for (i in seq(length(names(covariants$continuous)))){
    levels_tmp <- paste(names(covariants$continuous)[i],seq(as.numeric(covariants$continuous[[i]][2])), sep="")
    coldata[,covariants$continuous[[i]][1]] = cut(coldata[,covariants$continuous[[i]][1]],
                                                as.numeric(covariants$continuous[[i]][2]),
                                                labels=levels_tmp)
  }
}

#Sort names countdata as rownames in coldata
countdata <- countdata[rownames(coldata)]
## BATCH EFFECT ##
dds <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = formula(eval(parse(text=mod))))

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

## RELEVEL CONTROL GROUP ##
dds[[group]] <- relevel(dds[["group"]], control)

## DIFFERENTIAL ANALYSIS ##
if (!args$onlypca){
  dds <- DESeq(dds, parallel=TRUE)  
}



############## POLTS SENSE CONTRASTOS ############################################################
## rlog transformation and variance stabilization ##
#other normalization methods for plotting#

rld <- rlog(dds)
#vsd <-varianceStabilizingTransformation(dds)
rlogMat<-assay(rld)
write.table(rlogMat, file = paste(args$project,"_rlogMat.txt",sep=""), quote = F)
#vstMat<-assay(vsd)

#heatmap samples"
if ((args$plot)){
  sampleDists <- dist(t(rlogMat))
  sampleDistMatrix <- as.matrix(sampleDists)
  #colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pdf(paste(args$project,"_sampletosample_heatmap.pdf",sep=""),onefile=FALSE)
 #annotation_heatmap = as.data.frame(sapply(plot_atr$heatmap_ann,function(x) eval(parse(text=x))), row.names = colnames(coldata))
  annotation_heatmap=as.data.frame(colData(dds)[,plot_atr$heatmap_ann], row.names=row.names(colData(dds)))
  names(annotation_heatmap) = plot_atr$heatmap_ann
## SAVED: annotation_col = as.data.frame(coldata[,plot_atr$heatmap_ann]
  pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, 
           fontsize = 4, border_color = NA, annotation_col=annotation_heatmap, row.names = row.names(coldata))
  dev.off()
  print("sample correlation heatmap done")
}

## PCA!!!
#Select

rv <- rowVars(assay(rld))
if (args$plot){
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
pca <- prcomp(t(assay(rld)[select,]))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
dims <- combn(c(1:4),2)
tmp_color = coldata[,plot_atr$pca]
if (length(plot_atr$pca)>1){
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

  pdf(paste(args$project,"_pca.pdf", sep=""))
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

write.table( sweep(abs(pca$rotation), 2, colSums(abs(pca$rotation)), "/"), paste(args$project,'_pc_contribution.txt', sep=""), 
             sep = "\t", quote = F,row.names = T, col.names = T)

}


################## CONTRASTS FROM NOW ###############################################


## EXTRACT RESULTS ##
process_contrast <- function(title, contrastos){
  resAll <- results(dds, cooksCutoff=TRUE,contrast=contrastos, parallel=TRUE)
  res2 <- lfcShrink(dds, contrast = contrastos, res=resAll)
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
  sink(paste(args$project,'_',title,"_stats.txt",sep=""))
  mcols(res)$description
  summary(res)
  sink()
  
  ## WRITE TABLE RESULTS ##
  write.table(cc, paste(args$project,"_",title,"_norm_counts.txt",sep=""),quote=FALSE)
  pass_filter <- as.numeric(as.numeric(rownames(resAllOrdered) %in% rownames(resOrdered)) & (resAllOrdered$padj<0.05) )
  df_all <- as.data.frame(resAllOrdered)
  df_all["filter"]<- pass_filter
  df_all["shrunkenlfc"] = res2[rownames(df_all),"log2FoldChange"]
  df_all <- df_all[,c("baseMean","log2FoldChange","shrunkenlfc","lfcSE","stat", "filter", "pvalue", "padj")]
  write.table(df_all,paste(args$project,"_", title,"_results.txt",sep=""),quote=FALSE)
  print("DE analysis finished")

  if ((args$plot) & (length(rownames(resOrdered))>=plot_atr$de_genes_n)){  
    select=rlogMat[rownames(resOrdered),][1:plot_atr$de_genes_n,]
    data_subset=subset(colData(dds), dds[[contrastos[1]]] %in% c(contrastos[3], contrastos[2]))
    select2=select[,row.names(data_subset)]
    
    df=as.data.frame(data_subset[,plot_atr$heatmap_ann], row.names=row.names(data_subset))
    colnames(df)=c(plot_atr$heatmap_ann)
    
    pdf(paste(args$project,"_",title,"_top50DEgenes_heatmap.pdf",sep=""),onefile=FALSE)
    
    pheatmap(select2, fontsize_row=6,show_rownames=TRUE, cluster_rows=TRUE, cluster_cols=TRUE, annotation_col=df,scale="row", fontsize=4, show_colnames = FALSE)
    dev.off()
    print("top 50 DE genes heatmap done")

    tiff(paste(args$project,"_heatmap_custom.tiff", sep=""), width=2000, height=2000, res=300)
    pheatmap(select2, fontsize_row=6,show_rownames=TRUE, cluster_rows=TRUE, cluster_cols=FALSE, annotation_col=df,scale="row", fontsize=4, show_colnames = TRUE)
    dev.off()

  }
}

if (!args$onlypca){
  sapply(names(contrast), function(x) process_contrast(x, contrast[[x]]))  
}
# 




