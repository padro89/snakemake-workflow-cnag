#save.image(file="workspace")

# Load library
library('gprofiler2')

# Get config
species <- snakemake@config$species
baseurl <- snakemake@config$baseurl
orth_species <- snakemake@config$orth_species

# Check some defaults
if(is.null(species)){
  species<-"hsapiens"
  print("No species specified. Defaulting to H. sapiens")
}

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

# This is not necessary, as snakemake gives one file at a time:
#files=Sys.glob(opt$input)
#for (i in files){
#  if (as.numeric(unlist(strsplit(system(paste("wc -l ", i, sep=""), intern=TRUE)," "))[1] )>0){
#    name_f=basename(i)
#   cat(paste(name_f,"\n",sep=""))
#  gene=read.table(i, header = F)$V1
# 
#  if (opt$orth_species != FALSE){
#    gene = gorth(gene, source_organism=opt$species, target_organism=opt$orth_species)$ortholog_ensg
#    opt$species <- opt$orth_species
#  }
#    
#    result=gost(
#      gene,
#      organism = opt$species,
#      evcodes = TRUE
#    )
#    #result$result$parents <- lapply(result$result$parents, function(x) paste(x, collapse=","))
#    write.table(result$result[,colsExport],paste(opt$outdir,"/",name_f,'.table',sep=""), quote= FALSE, sep="\t", row.names = F)
#  }
#}

# Get DEGs
genes <- read.table(snakemake@input$deg_list, header = F)[,1]

if(length(genes > 0)){
  # Use specified organism if Orth_species == T
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
}
# Session info. Should I use it somewhere else?
#sink(format(Sys.time(), "sessionInfo-GetGO_gProfiler2_%Y_%b_%d_%H%M%S.txt"))
#sessionInfo()
#sink()