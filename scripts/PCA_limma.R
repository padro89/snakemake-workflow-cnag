# save.image(file="workspace")

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

## Load config
group <- snakemake@config$group
plot <- snakemake@config$plot
if(snakemake@config$plot_atr$pca == snakemake@config$group |
   is.null(snakemake@config$plot_atr$pca)){
  pca_atr <- "group"
}else{
  pca_atr <- snakemake@config$plot_atr$pca
}
heatmap_atr <- snakemake@config$plot_atr$heatmap_an
factors <- snakemake@config$factors
continuous <- snakemake@config$continuous
n_continuous_splits <- snakemake@config$n_continuous_splits

# Load count matrix and info
counts <- read.table(snakemake@input$counts,
                     header = T, row.names = 1)
info <- read.table(snakemake@input$info,
                   header = T, row.names = 1)

# Specify factors
if (!is.null(factors)){
  for (i in seq(length(factors))){
    info[,factors[i]]=as.factor(info[,factors[i]])
  }
}

# Convert continuous to factors
if (!is.null(continuous)){
  for(i in seq(length(continuous))){
    levels_temp <- paste(continuous[i],seq(n_continuous_splits[i]),sep="_")
    info[,continuous[i]] <- cut(info[,continuous[i]],
                                breaks = n_continuous_splits[i],
                                labels=levels_temp)
  }
}

# Adapt counts and info
counts <- counts[,names(counts) %in% row.names(info)]
info <- info[rownames(info) %in% names(counts),]
info <- info[names(counts),]
names(info)[which(names(info)==group)] = "group"
counts <- counts[,rownames(info)]
group <- factor(info[,"group"])



# Non-specific filtering
y <- DGEList(counts = counts,
             samples = info)
isexpr <- rowSums(cpm(y) > 1) >= 3
y <- y[isexpr,keep.lib.size=FALSE]
y=calcNormFactors(y)
dim(y)

#heatmap samples"
if (plot == T){
  sampleDists <- dist(t(y$counts))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pdf(file = snakemake@output$sampletosample,
      onefile=FALSE)
  annotation_heatmap <- as.data.frame(info[,heatmap_atr], row.names=row.names(info))
  names(annotation_heatmap) <- heatmap_atr
  pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, 
           clustering_distance_cols=sampleDists, col=colors, 
           fontsize = 4, border_color = NA, 
           annotation_col=annotation_heatmap, row.names = row.names(info))
  dev.off()
  print("sample correlation heatmap done")
}

# PCA
rv <- rowVars(y$counts)
if (plot == T){
  select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
  pca <- prcomp(t(y$counts[select,]))
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  dims <- combn(c(1:4),2)
  tmp_color <- info[,pca_atr]
  if (length(pca_atr)>1){
    ## Aquí es podria posar que guardés un color_factor numerat per cadascun
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
  
  pdf(file = snakemake@output$pcas)
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
  
  write.table(x = sweep(abs(pca$rotation), 2, colSums(abs(pca$rotation)), "/"),
             file = snakemake@output$pc_contribution,
             sep = "\t", quote = F,row.names = T, col.names = T)

  
}

save(y, file = snakemake@output$dge)

