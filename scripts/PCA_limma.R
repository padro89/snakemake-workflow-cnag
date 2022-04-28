save.image(file="workspace")

# Load necessary libraries
library(limma)
library(edgeR)
library(sva)
library(pheatmap)
library(BiocParallel)

# Load config
group <- snakemake@config$group

# Load count matrix and info
counts <- read.table(snakemake@input$counts,
                     header = T, row.names = 1)
info <- read.table(snakemake@input$info,
                   header = T, row.names = 1)

# Adapt counts and info
counts <- counts[,names(counts) %in% row.names(info)]
info <- info[rownames(info) %in% names(counts),]
info <- info[names(counts),]
counts <- counts[,rownames(info)]

group <- factor(info[,group])

y <- DGEList(counts = counts)
isexpr <- rowSums(cpm(y) > 1) >= 3
y <- y[isexpr,keep.lib.size=FALSE]
y=calcNormFactors(y)
dim(y)

#heatmap samples"
if (plot == T){
  sampleDists <- dist(t(rlogMat))
  sampleDistMatrix <- as.matrix(sampleDists)
  #colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pdf(snakemake@output$sampletosample,
      # paste(file.path(directory,"Results/"),project,"_sampletosample_heatmap.pdf",sep=""),
      onefile=FALSE)
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
  
  pdf(snakemake@output$pcas)    #(paste(file.path(directory,"Results/"),project,"_pca.pdf", sep=""))
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
  
  write.table(sweep(abs(pca$rotation), 2, colSums(abs(pca$rotation)), "/"),
              file = snakemake@output$pc_contribution,
              # paste(file.path(directory,"Results/"),project,'_pc_contribution.txt', sep=""), 
              sep = "\t", quote = F,row.names = T, col.names = T)
  
}

pdf(snakemake@output$mds)
plotMDS(y,col=as.numeric(group),dim.plot = c(1,2), labels=group)
dev.off()