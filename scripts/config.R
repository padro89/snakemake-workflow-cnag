covariants  <- list(
  factors = list(
    sex='SEX'
    )

  )


#mod <- "~group"
mod <- "~SEX + group"

contrast <- list(affected_vs_unaffected = c("group", "Affected", "Unaffected")		  
                )



plot_atr <- list(
  pca =  c('group'),
  heatmap_ann = c('group'),
  de_genes_n = 50
)

