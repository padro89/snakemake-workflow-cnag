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

# Non-specific filtering
y <- DGEList(counts = counts)
isexpr <- rowSums(cpm(y) > 1) >= 3
y <- y[isexpr,keep.lib.size=FALSE]
y=calcNormFactors(y)
dim(y)

pdf(snakemake@output$mds)
plotMDS(y,col=as.numeric(group),dim.plot = c(1,2), labels=group)
dev.off()
 
# Model matrix
mod <- model.matrix(formula(snakemake@config$formula))
colnames(mod)=gsub("group","",colnames(mod))

# Contrasts
contrasts_string <- vector(length=length(snakemake@config$contrasts))
for(i in 1:length(snakemake@config$contrasts)){
  contrasts_string[i] <- paste(names(snakemake@config$contrasts)[i],
                               "=",
                               snakemake@config$contrasts[[i]])
}
# 
# for(i in 1:length(snakemake@config$contrasts)){
#   contrasts_string[i] <- paste(names(snakemake@config$contrasts)[i],
#                                "=",
#                                snakemake@config$contrasts[[i]][1],
#                                "-",
#                                snakemake@config$contrasts[[i]][2])
# }
# 
# 
# if(!is.null(snakemake@config$interaction)){
#   interactions_string <- vector(length=length(snakemake@config$interaction))
#   for(i in 1:length(snakemake@config$interaction)){
#     interactions_string[i] <- paste(names(snakemake@config$interaction)[i],
#                                     "= (",
#                                     snakemake@config$interaction[[i]][1],
#                                     "-",
#                                     snakemake@config$interaction[[i]][2],
#                                     ") - (",
#                                     snakemake@config$interaction[[i]][3],
#                                     "-",
#                                     snakemake@config$interaction[[i]][4],
#                                     ")")
#   }
# }
# 
# contrasts_string <- append(contrasts_string, interactions_string)
contr.matrix <- makeContrasts(contrasts = contrasts_string,
                              levels = colnames(mod))
colnames(contr.matrix) <- names(snakemake@config$contrasts)

# Voom transformation and model fitting with Bayesian error estimation
v=voom(y,mod)
fit=lmFit(v,mod)
fit=contrasts.fit(fit, contrasts=contr.matrix)
fit2=eBayes(fit)
summary(decideTests(fit2))

for (i in colnames(fit2$coefficients)){
  #snakemake@output[i] <- write.table(
  top <- topTable(fit2,coef=i,sort="p",n=Inf)
  dir.create(paste0(snakemake@config$path$dge,"/", i))
  write.table(top, 
              file = paste0(snakemake@config$path$dge,
                            "/", i, "/", i, "_deg_results.txt"),
              quote=F)
}

# for (i in  colnames(fit2$coefficients)){
# top=topTable(fit2,coef=i,sort="p", n=Inf)
# genes=rownames(top[which(top$adj.P.Val<0.05),])
# write.table(top,paste(i,"_limma_voom_NADIFNAE_02_results",sep=""),quote=F)}

for (i in  colnames(fit2$coefficients)[1:4]){
top=topTable(fit2,coef=i,sort="p", n=Inf)
genes=rownames(top[which(top$adj.P.Val<0.05),])
term1=strsplit(i,split="_vs_")[[1]][1]
term2=strsplit(i,split="_vs_")[[1]][2]
samples=rownames(subset(info,GROUP==term1 | GROUP==term2))
expr=v$E[genes,samples]
rownames(expr)=do.call(rbind, strsplit(genes, ','))[,2]
if (length(genes) >1) {
pdf(paste("pheatmap_top_DE_genes_NADIFANE_02",i,".pdf",sep=""))
pheatmap(expr,scale="row",annotation_col=info[,c("GROUP","SAMPLE")], border_color = "NA",show_rownames = F)
dev.off()
}}

write.table(v$E,"logcpm_NADIFNAE_02.txt",quote=F)
