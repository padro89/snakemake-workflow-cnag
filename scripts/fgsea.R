save.image(file="workspace")

# Get config
# species <- snakemake@config$species
gmt <- snakemake@config$gmt
# 
# # Adapt species name
# if(is.null(species)){
#   species<-"hsapiens"
#   print("No species specified. Defaulting to H. sapiens")
# }else{
#   species <- unlist(strsplit(species, " "))
#   species <- paste0(tolower(substr(species[1],1,1)),
#                     species[2])
# }


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

# 
# # Getting the split character:
# # Defaults to coma
# split <- "\\,"
# # If ENSEMBLEID is used
# if(colid == 1){
#   # If human or mouse
#   if(snakemake@config$species == "hsapiens" |
#      snakemake@config$species == "mmusculus"){
#     # Sets it to point
#     split <-"\\."
#   }
# }




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

# Remove duplicates
ranks <- unique(ranks, by = "ID")
ranks <- setNames(ranks$Rank, ranks$ID)
ranks <- ranks[!is.na(ranks)]

## Get GMT pathways
allLevels <- gmtPathways(snakemake@input$gmt)

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

## Save table of results
fwrite(fgseaRes[padj<0.05][order(pval,NES)], 
          file = snakemake@output$fgsea,
          sep = "\t",
          sep2 = c("", " ", ""))

#sink(format(Sys.time(), "sessionInfo-fgsea_%Y_%b_%d_%H%M%S.txt"))
#sessionInfo()
#sink()