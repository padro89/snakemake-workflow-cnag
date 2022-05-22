#save.image(file="workspace")

# Get species
species <- snakemake@config$species
species <- unlist(strsplit(species, " "))
species <- paste0(tolower(substr(species[1],1,1)),
                  species[2])

# Get GMT
if(is.null(snakemake@config$gmt)){
  gmt <- snakemake@input$gmt
}else{
  gmt <- snakemake@config$gmt
}

# Get rank column
if(is.null(snakemake@config$rank)){
  if(snakemake@config$dge_method == "deseq2"){
    rank <- "stat"
  }else{
    rank <- "t"
  }
}else{
  rank <- snakemake@config$rank
}

# Print ranking criteria
if(is.character(rank)){
  print(paste("Ranking by:", rank))
}else{
  rankName <- colnames(raw_table)[rank]
  print(paste("Ranking by:", rankName))
}

# Comparation name
comp_name <- unlist(strsplit(snakemake@input$deg_results, "/"))
comp_name <- strsplit(comp_name[length(comp_name)], "\\.")[[1]][1]
comp_name <- sub("_deg_results", "", comp_name)
print(paste("Running FGSEA for:", comp_name))

# Setting ENSEMBL or SYMBOL
if(snakemake@config$id_type == "ENSEMBL" || is.null(snakemake@config$id_type)){
  print("Using ENSEMBLID")
  colid <- 1
}else{
  if(snakemake@config$id_type == "SYMBOL"){
    print("Using SYMBOL")
    colid <- 2
  }else{
    print("Unknown gene identifier. Defaulting to ENSEMBLID")
    colid <- 1
  }
}

### LIBRARIES
library('fgsea')
#library("repr")
library('data.table')
library("BiocParallel")
library("gprofiler2")

# Check need of orthologous
gmt_species <- c("Bos taurus", "Caenorhabditis elegans", "Canis familiaris",
                 "Danio rerio", "Dictyostelium discoideum", "Drosophila melanogaster",
                 "Gallus gallus","Homo sapiens", "Mus musculus", "Plasmodium falciparum", 
                 "Rattus norvegicus", "Saccharomyces cerevisiae",
                 "Schizosaccharomyces pombe", "Sus scrofa", "Xenopus tropicalis")

if(snakemake@config$species %in% gmt_species && snakemake@config$orth_species == F){
  orthologous <- F
}else{
  orthologous <- T
}
  
### FGSEA
## Get ranked list
raw_table <- read.table(snakemake@input$deg_results, sep=" ", header=T)
# Get the desired column
ranks <- data.frame(list(ID=sapply(strsplit(rownames(raw_table),
                                            split = "\\,"),
                                   function(x) unlist(x)[colid]),
                         Rank=raw_table[,rank]))

# Get everything before the point in the identifier, if it exists
ranks$ID <- sub("\\.[^.]*$", "", ranks$ID)
ranks <- unique(ranks, by = "ID")

# Use orthologous genes if needed (defaults to human)
if(orthologous == T){
  print("Using orthologous genes.")
  # Use alternative orthologous species if one is provided
  if(snakemake@config$orth_species != F){
    target_organism <- snakemake@config$orth_species
  }else{
    target_organism <- "hsapiens"
  }
  genes <- gorth(ranks$ID, source_organism=species, target_organism=target_organism)
  # Merge the result in a dataframe
  merged <- merge(ranks, genes, by.x = "ID", by.y = "input_ensg")
  # Create the ranked list
  ranks <- data.frame(ID = merged$ortholog_ensg,
                      Rank = merged$Rank)
}

# Make the ranked list
ranks <- setNames(ranks$Rank, ranks$ID)
ranks <- ranks[!is.na(ranks)]

## Get GMT pathways
allLevels <- gmtPathways(gmt)

## Run FGSEA
fgseaRes <- fgsea(allLevels, ranks, eps = 0.0,
                  minSize = 15, maxSize=500)


## Plot results
pdf(snakemake@output$fgsea_pdf, width = 15, height = 10)
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=20), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=20), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(allLevels[topPathways], ranks, fgseaRes,
              gseaParam = 0.5, colwidths = c(5, 3,0.8, 1.2, 1.2))
dev.off()

# Collapsed pathways
pdf(snakemake@output$fgsea_main_pathways_pdf, width = 15, height = 10)
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01], allLevels, ranks)
mainPathways <-fgseaRes[pathway %in% collapsedPathways$mainPathways][order(-NES),pathway]
plotGseaTable(allLevels[mainPathways], ranks, fgseaRes, gseaParam = 0.5, 
              colwidths = c(5, 3,0.8, 1.2, 1.2))
dev.off()

## Save table of results
fwrite(fgseaRes[padj<0.05][order(pval,NES)], 
          file = snakemake@output$fgsea,
          sep = "\t",
          sep2 = c("", " ", ""))

#sink(format(Sys.time(), "sessionInfo-fgsea_%Y_%b_%d_%H%M%S.txt"))
#sessionInfo()
#sink()