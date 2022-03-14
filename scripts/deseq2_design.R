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
