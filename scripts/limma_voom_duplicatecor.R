save.image(file="workspace")

library(edgeR)
library(limma)
library(ggrepel)
library(ggplot2)
library(pheatmap)

counts <- read.table(snakemake@input$counts,
                     header = T, row.names = 1,
                     check.names=FALSE)
info <- read.table(snakemake@input$info,
                   header = T, row.names = 1)

counts=counts[,names(counts) %in% row.names(info)]
info=info[rownames(info) %in% names(counts),]
info=info[names(counts),]
counts=counts[,rownames(info)]

y=DGEList(counts=counts)
A<-rowSums(y$counts)
isexpr<-A>500
y=y[isexpr,keep.lib.size=FALSE]
dim(y)

y=calcNormFactors(y)

group=as.factor(info$GROUP)
individual=as.factor(info$PACIENT)
	  

### model1 (adjusting for patient)
mod=model.matrix(~0+group)
colnames(mod)=gsub("group","",colnames(mod))
contr.matrix <- makeContrasts(iDC_C4BPb_vs_iDC_V4=iDC_C4BPb-iDC_V4,
                              iDC_V4_vs_iDC=iDC_V4-iDC,
							  iDC_C4BPb_vs_iDC=iDC_C4BPb-iDC,
							  M0_C4BPb_vs_M0_V4=M0_C4BPb-M0_V4,
                              M0_V4_vs_M0=M0_V4-M0,
							  M0_C4BPb_vs_M0=M0_C4BPb-M0,
							  OC_C4BPb_vs_OC_V4=OC_C4BPb-OC_V4,
                              OC_V4_vs_OC=OC_V4-OC,
							  OC_C4BPb_vs_OC=OC_C4BPb-OC,
							  levels=mod) 
v=voom(y,mod)
corfit <- duplicateCorrelation(v,mod, block=individual)
v <- voom(y, mod, block = individual, correlation = corfit$consensus)
#fit <- lmFit(v, mod)
fit <- lmFit(v, mod, block = individual, correlation = corfit$consensus)
fit=contrasts.fit(fit, contrasts=contr.matrix)
fit2=eBayes(fit)
summary(decideTests(fit2))

write.table(v$E,"logcpm_aranjos_01.txt",quote=F)

for (i in  colnames(fit2$coefficients)){
top=topTable(fit2,coef=i,sort="p", n=Inf)
write.table(top,paste(i,"_limma_voom_dup_aranjos_01",sep=""),quote=F)
genes=rownames(top[which(top$adj.P.Val<0.05),])

 term1=strsplit(i,split="_vs_")[[1]][1]
 term2=strsplit(i,split="_vs_")[[1]][2]
 samples=rownames(subset(info,GROUP==term1 | GROUP==term2))
 expr=v$E[genes,samples]
 rownames(expr)=do.call(rbind, strsplit(genes, ','))[,2]
 if (length(genes) >1) {
 pdf(paste("pheatmap_top50_DE_genes_aranjos_01_",i,".pdf",sep=""))
 pheatmap(expr,scale="row",annotation_col=info[,c("GROUP","PACIENT")], border_color = "NA",show_rownames = F)
 dev.off()
 }}