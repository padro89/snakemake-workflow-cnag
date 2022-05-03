#save.image(file="workspace")

# Load library
library('gprofiler2')

# Get config
species <- snakemake@config$species
baseurl <- snakemake@config$baseurl
orth_species <- snakemake@config$orth_species

# Adapt species name
if(is.null(species)){
  species<-"hsapiens"
  print("No species specified. Defaulting to H. sapiens")
}else{
  species <- unlist(strsplit(species, " "))
  species <- paste0(tolower(substr(species[1],1,1)),
                    species[2])
}


# Check some defaults
if(is.null(baseurl)){
  baseurl<-"https://biit.cs.ut.ee/gprofiler"
  print("No baseurl specified. Defaulting to https://biit.cs.ut.ee/gprofiler")
}

if(is.null(orth_species)){
  orth_species<-FALSE
}

# Setting base URL
set_base_url(baseurl)

# Select columns to export
colsExport <- c(
  "query",
  "significant",
  "p_value",
  "term_size",
  "query_size",
  "intersection_size",
  "precision",
  "recall",
  "term_id",
  "source",
  "term_name",
  #  "parents",
  "intersection"
)

# Get DEGs
genes <- try(read.table(snakemake@input$deg_list, header = F)[,1],
             silent = F)

if(class(genes) != "try-error"){
  if(orth_species != FALSE){
    genes <- gorth(gene, source_organism=species, target_organism=orth_species)$ortholog_ensg
    species <- orth_species
  }

  result <- gost(
    genes,
    organism = species,
    evcodes = TRUE)

  write.table(result$result[,colsExport],snakemake@output$ora,
              quote= FALSE, sep="\t", row.names = F)
} else {
  write.table("Not enough differentialy expressed genes to perform the analysis",
              snakemake@output$ora,
              quote= FALSE, row.names = F)
}
# Session info. Should I use it somewhere else?
#sink(format(Sys.time(), "sessionInfo-GetGO_gProfiler2_%Y_%b_%d_%H%M%S.txt"))
#sessionInfo()
#sink()