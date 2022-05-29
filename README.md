# Snakemake whole-RNA-seq pipeline

Repository used to create a Snakemake transcriptomics pipeline for CNAG.

To use the pipeline, simply copy the folder to the desired computer, install Snakemake and Conda (preferably with Mamba) package managers and instal the complete environment from the mamba_envir.yaml file. 

Every aspect of the pipeline should be managed through the edition of the config.yml, that contains comments and examples. 

The Pipeline starts from the tab separated tables containing the counts matrix and the sample information. If the option `onlypca` is set to true, it outputs a graphical representation of the PCA and a sample-to-sample heatmap. After proper analysis of the results, the model formula can be specified as a string in the config.yaml, as well as the package to perform the differential expression analysis (DESeq2 or limma-voom) and whether to perform over-repesentation analysis with g:Profiler2, gene set enrichment analysis with fgsea or both. The species of the samples should be specified in the config.yaml to perform the enrichment analysis, as well as a GMT file for the fgsea package and an orthologous species if the target species is not sufficiently documented. If no GMT file is provided, the pipeline will automatically make one from the latest Reactome pathways, filtering by the species or using orthologous human genes if the species is not in the Reactome database and no orthologous species is provided.

To run Snakemake from the folder where the Snakefile is saved, simply run:

```
snakemake --cores <desired number of cores>
```

And the desired pipeline will take place according to the settings in the config.yaml.

# What's new
- Corrected some bugs (see below).
- Added plot of the main pathways from the FGSEA.

# Bugs
- If plot is false there will be errors for sure. I should make conditional output in each rule according to this parameter.
- I should check the limma PCA and sample to sample heatmap. I don't think it works properly. *Updated:* I should run voom transformation before plotting the PCA.
- The 'tiff' argument doesn't do anything. I should create conditional output in each rule according to this parameter.
- Orthologous species only work with ENSEMBL ID, not SYMBOL, in FGSEA.

## Fixed bugs
- There was a bug that caused that if one wrote "false" in the orth species field instead of leaving it blank, it tried to find orth species. I changed the line in the fgsea script and now it works.
- Fixed a bug that caused gProfiler to only get ENSEMBL ids (it was a problem with the awk command). Now it gets both ENSEMBL and SYMBOL and uses the one specified. 

# Planned activites
- Fix the bugs.
- Add conda dependencies for each rule in case no prefix is specified.
- Integrate an alignment and quantification module.
- Add the project name to each file.
- Add volcano plots of the DGE results.
- Add ranked graphs to the FGSEA.
- Make different PCAs for different variables (a for look inside the PDF, everything in the same PDF to avoid errors).