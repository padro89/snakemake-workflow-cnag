# Transcriptomics workflow
Repository used to create a Snakemake transcriptomics workflow for CNAG


# Some stuff I'm using
## Conda env to download specific R and Python packages versions and test the scripts.

- I created a .yaml file and and environment adapted to de DESeq2 version (1.26). Conda (or Mamba) automatically resolved other packages versions.
- The .yaml file had several fields. `name`, `channels` (where I specified conda-forge and bioconda) and `dependencies`, where I specified the packages to install into the new environment.
- To create the environment, I ran `mamba env create -f <yaml_file>`.
- I'm going to try to run the script.
- Any time I manually add a package, I'll add it to my conda environment yaml file.

## Version control by use of github

- I created a token from my github account to access from git.
- I created an empty repository called snakemake-workflow-cnag
- I cloned my empty github repository to a local folder.
- I created the project's first branch, setup, with `git branch setup` where I'm trying to implement everything I'm doing.
- To swap between branches, I used `git checkout <branch>` altough I think I could just `git push --set-upstream origin <branch>`.
- Each time I finish working, I `git add .` (or at least whatever I want to upload) and then `git commit -m "Description"` to commit the changes. Before comitting, I can check which files will be updated by `git status`. Finally, I `git push` to commit the changes to the repository. 
- Each time I resume working in my PC, I `git pull` and all the folders in my PC are automatically updated from the repository.
- Eventually, when setup is finished, I may need to merge it to master.
- Then, I'll just create another branch for anything I'm trying to do, develop it, test it and merge it.

# What I'd like to do
## Implement the DESEQ2 1.26.0 script by modularizing it

- Arguments should be removed in favor of config files.
- Ideally, one script per step.
- Snakemake S4 objects should be used inside the scripts to get results from scripts run in other rules.
- PCA should be shown first, and only afterwards should one run the complete Snakemake.

## Implement the pathway analysis using Snakemake

## Update the  script to include DESeq latest version

## Implement the option to use Limma instead of DESEQ2

## Some kind of comparation or Benchmark
