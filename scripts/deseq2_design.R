# Load required libraries
library("DESeq2")
library("BiocParallel")

# Import the design
formula <- snakemake@config$formula

# DESeq with the design
load(snakemake@input$dds)
dds@design <- formula(snakemake@config$formula)
dds <- DESeq(dds, parallel=TRUE)
save(dds, file = snakemake@output$dds_design)


# Extract normalized counts

norm_counts <- counts(dds,normalized=TRUE)
row_sums <- rowSums(norm_counts)
norm_counts_ordered <- norm_counts[order(row_sums, decreasing = T),]
colnames(norm_counts_ordered) <- paste0(dds$group,",",
                                        colnames(norm_counts_ordered))
counts_dds=(norm_counts_ordered[,order(colnames(norm_counts_ordered))])
rounded_counts=round(counts_dds,digits=2)
write.table(rounded_counts, file = snakemake@output$norm_counts,quote=FALSE)