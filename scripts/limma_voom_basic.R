save.image(file="workspace")

# Load necessary libraries
library(limma)
require(statmod)
library(BiocParallel)
library(pheatmap)
# library(gridExtra)
# library(grid)
library(edgeR)
# library(sva)
# 
# 
# library(genefilter)
# library(RColorBrewer)
# library(ggplot2)

load(snakemake@input$dge)

interaction_id <- snakemake@config$interaction_id
info <- y$samples
if(is.null(snakemake@config$samples)){
  info$samples <- rownames(y$samples)
  samples_var <- "samples"
}else{
  samples_var <- snakemake@config$samples
}

# Model matrix
mod <- model.matrix(formula(snakemake@config$formula), data = y$samples)
colnames(mod)=gsub("group","",colnames(mod))

# Contrasts from config
contrasts_string <- vector(length=length(snakemake@config$contrasts))
for(i in 1:length(snakemake@config$contrasts)){
  contrasts_string[i] <- paste(names(snakemake@config$contrasts)[i],
                               "=",
                               snakemake@config$contrasts[[i]])
}

# Make contrasts
contr.matrix <- makeContrasts(contrasts = contrasts_string,
                              levels = colnames(mod))
colnames(contr.matrix) <- names(snakemake@config$contrasts)


# Voom transformation and model fitting with Bayesian error estimation
v=voom(y,mod)

# Blocking
if(!is.null(snakemake@config$blocking)){
  blocking <- as.factor(info[,snakemake@config$blocking])
  corfit <- duplicateCorrelation(v,mod, block=blocking)
  v <- voom(y, mod, block = blocking, correlation = corfit$consensus)
  fit <- lmFit(v, mod, block = blocking, correlation = corfit$consensus)
}else{
  fit <- lmFit(v,mod)
}

# Contrasts and bayesian correction of the errors.
fit=contrasts.fit(fit, contrasts=contr.matrix)
fit2=eBayes(fit)
summary(decideTests(fit2))

# Create and save topTables
for (i in colnames(fit2$coefficients)){
  top <- topTable(fit2,coef=i,sort="p",n=Inf)
  dir.create(paste0(snakemake@config$path$dge,"/", i),
             showWarnings = F)
  write.table(top,
              file = paste0(snakemake@config$path$dge,
                            "/", i, "/", i, "_deg_results.txt"),
              quote=F)
}

# Create heatmaps for each comparation
interactions <- grep(interaction_id,
                     colnames(fit2$coefficients),
                     ignore.case = T)
if(length(interactions) < 1){
  for (i in colnames(fit2$coefficients)){
    top <- topTable(fit2,coef=i,sort="p", n=Inf)
    genes <- rownames(top[which(top$adj.P.Val<0.05),])
    term1 <- strsplit(i,split="_vs_")[[1]][1]
    term2 <- strsplit(i,split="_vs_")[[1]][2]
    samples <- rownames(subset(y$samples,group==term1 | group==term2))
    expr <- v$E[genes,samples]
    rownames(expr) <- do.call(rbind, strsplit(genes, ','))[,2]
    if (length(genes) >1) {
      pdf(file = paste0(snakemake@config$path$dge,
                        "/", i, "/", i, "_heatmap.pdf"))
      pheatmap(expr,scale="row",annotation_col=info[,c("group",samples_var)], 
               border_color = "NA",show_rownames = F)
      dev.off()
    } else {
      # Funky function to make an empty PDF if not enough genes
      pdf(file = paste0(snakemake@config$path$dge,
                        "/", i, "/", i, "_heatmap.pdf"))
      plot(NA, xlim= c(0,5), ylim=c(0,5), bty = "n",
           xaxt = "n", yaxt = "n", xlab ="", ylab = "")
      text(2,5,print("There were no differentially expressed genes for this comparation"))
      dev.off()
    }
  }
} else {
  for (i in colnames(fit2$coefficients)[-interactions]){
    top <- topTable(fit2,coef=i,sort="p", n=Inf)
    genes <- rownames(top[which(top$adj.P.Val<0.05),])
    term1 <- strsplit(i,split="_vs_")[[1]][1]
    term2 <- strsplit(i,split="_vs_")[[1]][2]
    samples <- rownames(subset(y$samples,group==term1 | group==term2))
    expr <- v$E[genes,samples]
    rownames(expr) <- do.call(rbind, strsplit(genes, ','))[,2]
    if (length(genes) >1) {
      pdf(file = paste0(snakemake@config$path$dge,
                        "/", i, "/", i, "_heatmap.pdf"))
      pheatmap(expr,scale="row",annotation_col=info[,c("group",samples_var)], 
               border_color = "NA",show_rownames = F)
      dev.off()
    } else {
      # Funky function to make an empty PDF if not enough genes
      pdf(file = paste0(snakemake@config$path$dge,
                        "/", i, "/", i, "_heatmap.pdf"))
      plot(NA, xlim= c(0,5), ylim=c(0,5), bty = "n",
           xaxt = "n", yaxt = "n", xlab ="", ylab = "")
      text(2,5,print("There were no differentially expressed genes for this comparation"))
      dev.off()
    }
    for (i in colnames(fit2$coefficients)[interactions]){
      pdf(file = paste0(snakemake@config$path$dge,
                        "/", i, "/", i, "_heatmap.pdf"))
      plot(NA, xlim= c(0,5), ylim=c(0,5), bty = "n",
           xaxt = "n", yaxt = "n", xlab ="", ylab = "")
      text(2,5,print("No Heatmap was created for the interaction"))
      dev.off()
    }
  }
}



# Create empty PDF for interactions
# 
# if(length(interactions) > 1){
#     
# }

# Save logCPM file
write.table(v$E, snakemake@output$logcpm,quote=F)
