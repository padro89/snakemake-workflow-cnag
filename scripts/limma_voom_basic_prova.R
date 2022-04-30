save.image(file="workspace")

# Load necessary libraries
library(limma)
library(gridExtra)
library(grid)
library(edgeR)
library(sva)
library(pheatmap)
library(BiocParallel)
library(genefilter)
library(RColorBrewer)
library(ggplot2)

# Load config
group <- snakemake@config$group
plot <- snakemake@config$plot
if(snakemake@config$plot_atr$pca == snakemake@config$group |
   is.null(snakemake@config$plot_atr$pca)){
  pca_atr <- "group"
}else{
  pca_atr <- snakemake@config$plot_atr$pca
}
heatmap_atr <- snakemake@config$plot_atr$heatmap_an

# Load count matrix and info
counts <- read.table(snakemake@input$counts,
                     header = T, row.names = 1)
info <- read.table(snakemake@input$info,
                   header = T, row.names = 1)

# Adapt counts and info
counts <- counts[,names(counts) %in% row.names(info)]
info <- info[rownames(info) %in% names(counts),]
info <- info[names(counts),]
names(info)[which(names(info)==group)] = "group"
counts <- counts[,rownames(info)]
group <- factor(info[,"group"])


# Non-specific filtering
y <- DGEList(counts = counts)
isexpr <- rowSums(cpm(y) > 1) >= 3
y <- y[isexpr,keep.lib.size=FALSE]
y=calcNormFactors(y)
dim(y)

#heatmap samples"
if (plot == T){
  sampleDists <- dist(t(y$counts))
  sampleDistMatrix <- as.matrix(sampleDists)
  #colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pdf(snakemake@output$sampletosample,
      # paste(file.path(directory,"Results/"),project,"_sampletosample_heatmap.pdf",sep=""),
      onefile=FALSE)
  #annotation_heatmap = as.data.frame(sapply(plot_atr$heatmap_ann,function(x) eval(parse(text=x))), row.names = colnames(coldata))
  annotation_heatmap <- as.data.frame(info[,heatmap_atr], row.names=row.names(info))
  names(annotation_heatmap) <- heatmap_atr
  ## SAVED: annotation_col = as.data.frame(coldata[,plot_atr$heatmap_ann]
  pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, 
           fontsize = 4, border_color = NA, annotation_col=annotation_heatmap, row.names = row.names(info))
  dev.off()
  print("sample correlation heatmap done")
}


rv <- rowVars(y$counts)
if (plot == T){
  select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
  pca <- prcomp(t(y$counts[select,]))
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  dims <- combn(c(1:4),2)
  tmp_color <- info[,pca_atr]
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
  
  pdf(snakemake@output$pcas)
  print(ggplot()+geom_text(aes(x=pca$x[,dims[1,1]], y=pca$x[,dims[2,1]], color = color_factor, label = as.factor(rownames(pca$x))))+
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
  
  toplot <- lapply(c(2:dim(dims)[2]), function(i) ggplot()+geom_point(aes(x=pca$x[,dims[1,i]], y=pca$x[,dims[2,i]], color = color_factor))+
                     xlab(paste(colnames(pca$x)[dims[1,i]]," (",round(percentVar[dims[1,i]]*100,digits = 2),"%)",sep=''))+
                     ylab(paste(colnames(pca$x)[dims[2,i]]," (",round(percentVar[dims[2,i]]*100,digits=2),"%)",sep=''))+
                     theme(legend.title = element_blank(), legend.position = 'none'))
  
  #grid.arrange(grobs = toplot, nrow=3)
  grid_arrange_shared_legend(toplot)
  
  dims <- combn(c(5:ifelse(dim(pca$x)[2]<8, dim(pca$x)[2],8)),2)
  toplot <- lapply(c(1:dim(dims)[2]), function(i) ggplot()+geom_point(aes(x=pca$x[,dims[1,i]], y=pca$x[,dims[2,i]], color = color_factor))+
                     xlab(paste(colnames(pca$x)[dims[1,i]]," (",round(percentVar[dims[1,i]]*100,digits = 2),"%)",sep=''))+
                     ylab(paste(colnames(pca$x)[dims[2,i]]," (",round(percentVar[dims[2,i]]*100,digits=2),"%)",sep=''))+
                     theme(legend.title = element_blank(),  legend.position = 'none'))
  
  
  #grid.arrange(grobs = toplot, nrow=3)
  grid_arrange_shared_legend(toplot)
  
  dev.off()
  
  write.table(sweep(abs(pca$rotation), 2, colSums(abs(pca$rotation)), "/"),
              file = snakemake@output$pc_contribution,
              # paste(file.path(directory,"Results/"),project,'_pc_contribution.txt', sep=""), 
              sep = "\t", quote = F,row.names = T, col.names = T)
  
}
# 
# # pdf(snakemake@output$mds)
# # plotMDS(y,col=as.numeric(group),dim.plot = c(1,2), labels=group)
# # dev.off()
#  
# # Model matrix
# mod <- model.matrix(formula(snakemake@config$formula))
# colnames(mod)=gsub("group","",colnames(mod))
# 
# # Contrasts
# contrasts_string <- vector(length=length(snakemake@config$contrasts))
# for(i in 1:length(snakemake@config$contrasts)){
#   contrasts_string[i] <- paste(names(snakemake@config$contrasts)[i],
#                                "=",
#                                snakemake@config$contrasts[[i]])
# }
# # 
# # for(i in 1:length(snakemake@config$contrasts)){
# #   contrasts_string[i] <- paste(names(snakemake@config$contrasts)[i],
# #                                "=",
# #                                snakemake@config$contrasts[[i]][1],
# #                                "-",
# #                                snakemake@config$contrasts[[i]][2])
# # }
# # 
# # 
# # if(!is.null(snakemake@config$interaction)){
# #   interactions_string <- vector(length=length(snakemake@config$interaction))
# #   for(i in 1:length(snakemake@config$interaction)){
# #     interactions_string[i] <- paste(names(snakemake@config$interaction)[i],
# #                                     "= (",
# #                                     snakemake@config$interaction[[i]][1],
# #                                     "-",
# #                                     snakemake@config$interaction[[i]][2],
# #                                     ") - (",
# #                                     snakemake@config$interaction[[i]][3],
# #                                     "-",
# #                                     snakemake@config$interaction[[i]][4],
# #                                     ")")
# #   }
# # }
# # 
# # contrasts_string <- append(contrasts_string, interactions_string)
# contr.matrix <- makeContrasts(contrasts = contrasts_string,
#                               levels = colnames(mod))
# colnames(contr.matrix) <- names(snakemake@config$contrasts)
# 
# # Voom transformation and model fitting with Bayesian error estimation
# v=voom(y,mod)
# fit=lmFit(v,mod)
# fit=contrasts.fit(fit, contrasts=contr.matrix)
# fit2=eBayes(fit)
# summary(decideTests(fit2))
# 
# for (i in colnames(fit2$coefficients)){
#   #snakemake@output[i] <- write.table(
#   top <- topTable(fit2,coef=i,sort="p",n=Inf)
#   dir.create(paste0(snakemake@config$path$dge,"/", i))
#   write.table(top, 
#               file = paste0(snakemake@config$path$dge,
#                             "/", i, "/", i, "_deg_results.txt"),
#               quote=F)
# }
# 
# # for (i in  colnames(fit2$coefficients)){
# # top=topTable(fit2,coef=i,sort="p", n=Inf)
# # genes=rownames(top[which(top$adj.P.Val<0.05),])
# # write.table(top,paste(i,"_limma_voom_NADIFNAE_02_results",sep=""),quote=F)}
# 
# for (i in  colnames(fit2$coefficients)[1:4]){
# top=topTable(fit2,coef=i,sort="p", n=Inf)
# genes=rownames(top[which(top$adj.P.Val<0.05),])
# term1=strsplit(i,split="_vs_")[[1]][1]
# term2=strsplit(i,split="_vs_")[[1]][2]
# samples=rownames(subset(info,GROUP==term1 | GROUP==term2))
# expr=v$E[genes,samples]
# rownames(expr)=do.call(rbind, strsplit(genes, ','))[,2]
# if (length(genes) >1) {
# pdf(paste("pheatmap_top_DE_genes_NADIFANE_02",i,".pdf",sep=""))
# pheatmap(expr,scale="row",annotation_col=info[,c("GROUP","SAMPLE")], border_color = "NA",show_rownames = F)
# dev.off()
# }}
# 
# write.table(v$E,"logcpm_NADIFNAE_02.txt",quote=F)
