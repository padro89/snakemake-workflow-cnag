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

pdf(snakemake@output$mds)
plotMDS(y,col=as.numeric(group),dim.plot = c(1,2), labels=group)
dev.off()
 
mod <- model.matrix(~0+group)
colnames(mod)=gsub("group","",colnames(mod))

contr.matrix <- makeContrasts(
WTC_TTX_vs_WTC_NT=WTC_TTX-WTC_NT,
X409B2_TTX_vs_X409B2_NT=X409B2_TTX-X409B2_NT,
X409B2_TTX_vs_WTC_TTX=X409B2_TTX-WTC_TTX,
X409B2_NT_vs_WTC_NT=X409B2_NT-WTC_NT,
diff_of_diff=(X409B2_TTX-X409B2_NT)-(WTC_TTX-WTC_NT),
levels=colnames(mod))

v=voom(y,mod)
fit=lmFit(v,mod)
fit=contrasts.fit(fit, contrasts=contr.matrix)
fit2=eBayes(fit)
summary(decideTests(fit2))

for (i in  colnames(fit2$coefficients)[1:5]){
top=topTable(fit2,coef=i,sort="p", n=Inf)
genes=rownames(top[which(top$adj.P.Val<0.05),])
write.table(top,paste(i,"_limma_voom_NADIFNAE_02_results",sep=""),quote=F)}

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