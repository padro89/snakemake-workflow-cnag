# Guardar el workspace per debugging amb snakemake
# Netejar el workspace abans
#eliminats <- ls()
#eliminats <- eliminats[eliminats!="snakemake" & eliminats!="Snakemake"]
#rm(list=eliminats)
#rm(eliminats)
# Guardar-lo amb només l'objecte snakemake actualitzat
#save.image(file="workspace")

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
group <- snakemake@config$group
control <- snakemake@config$control
project <- snakemake@config$project
plot <- snakemake@config$plot
factors <- snakemake@config$factors
continuous <- snakemake@config$continuous
pca_atr <- snakemake@config$plot_atr$pca
directory <- snakemake@config$path$pca
heatmap_atr <- snakemake@config$plot_atr$heatmap_ann

# No s'utilitzen
# contrastos <- snakemake@config$contrast
# shrinkage_method <-snakemake@config$shrinkage_method
# formula <- snakemake@config$formula

## Import counts matrix and sample info ##
counts_raw <- read.table(snakemake@input[["counts"]],
                         header = T, row.names = 1,
                         check.names = F)
info_raw <- read.table(snakemake@input[["info"]],
                       header = T, row.names = 1)
coldata=subset(info_raw, row.names(info_raw) %in% colnames(counts_raw))
countdata=counts_raw[,colnames(counts_raw) %in% rownames(coldata)]
coldata=subset(coldata, row.names(coldata) %in% colnames(countdata))

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

names(coldata)[which(names(coldata)==group)] = "group"

## RELEVEL CONTROL GROUP ##
coldata$group <- relevel(factor(coldata$group), ref=control)


#Sort names countdata as rownames in coldata
countdata <- countdata[rownames(coldata)]

## BATCH EFFECT ##
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = formula("~1"))

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

############## POLTS SENSE CONTRASTOS ############################################################
## rlog transformation and variance stabilization ##
#other normalization methods for plotting#

rld <- rlog(dds)
#vsd <-varianceStabilizingTransformation(dds)
rlogMat<-assay(rld)
write.table(rlogMat, file = snakemake@output$rlog # paste(file.path(directory,"Results/"),project,"_rlogMat.txt",sep="")
            , quote = F)
#vstMat<-assay(vsd)

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
save(dds, file = snakemake@output$dds)

