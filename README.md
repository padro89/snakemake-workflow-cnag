# Transcriptomics workflow
Repository used to create a Snakemake transcriptomics workflow for CNAG


# Some stuff I'm using
- Conda env to download specific R and Python packages versions and test the scripts.
- Version control by use of github.

# What I'd like to do
## Implement the DESEQ2 1.26.0 script by modularizing it. 
- Arguments should be removed in favor of config files.
- Ideally, one script per step.
- Snakemake S4 objects should be used inside the scripts to get results from scripts run in other rules.
- PCA should be shown first, and only afterwards should one run the complete Snakemake.

## Implement the option to use Limma instead of DESEQ2

## Some kind of comparation or Benchmark
