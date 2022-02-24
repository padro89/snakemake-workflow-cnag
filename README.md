# Transcriptomics workflow
Repository used to create a Snakemake transcriptomics workflow for CNAG


# Some stuff I'm using
## Conda env to download specific R and Python packages versions and test the scripts.

- I created a .yaml file and and environment adapted to de DESeq2 version (1.26). Conda (or Mamba) automatically resolved other packages versions.
- The .yaml file had several fields. `name`, `channels` (where I specified conda-forge and bioconda) and `dependencies`, where I specified the packages to install into the new environment.
- To create the environment, I ran `mamba env create -f <yaml_file>`.
- I ran the deseq2.R script using Rscript and the arguments, and it worked flawlesly without any need of adding any packages.
- If I need to  manually add a package, I'll add it to my conda environment yaml file.

## Version control by use of github

For every step of this project I'm going to use github for version control. Below I specify what I learned so far.

- I created a token from my github account to access from git.
- I created an empty repository called snakemake-workflow-cnag
- I cloned my empty github repository to a local folder.
- I created the project's first branch, setup, with `git branch setup` where I'm trying to implement everything I'm doing.
- To swap between branches, I used `git checkout <branch>` altough I think I could just `git push --set-upstream origin <branch>`.
- Each time I finish working, I `git add .` (or at least whatever I want to upload) and then `git commit -m "Description"` to commit the changes. Before comitting, I can check which files will be updated by `git status`. Finally, I `git push` to commit the changes to the repository. 
- Each time I resume working in my PC, I `git pull` and all the folders in my PC are automatically updated from the repository.
- I can delete the folder where I've cloned the repository, and then just clone it again.
- If I finish working on a branch, I can `git merge <branca_a_afegir_a_main>` to merge the branch into main. This must be done from master branch -otherwise it would just copy master into the other branch. It is important to update the branch before merging.
- Then, I'll just create another branch for anything I'm trying to do, develop it, test it and merge it.

### Branches I create

- `setup`: A branch where I created the conda environment and tested the script.
- `deseq2`: A branch where I try to implement the deseq2.R script in Snakemake.

# What I'd like to do
## Implement the DESEQ2 1.26.0 script by modularizing it

- Arguments should be removed in favor of config files.
- Ideally, one script per step.
- Snakemake S4 objects should be used inside the scripts to get results from scripts run in other rules.
- PCA should be shown first, and only afterwards should one run the complete Snakemake.

## Implement the pathway analysis using Snakemake

## Update the  script to include DESeq latest version

- There is some problem with the `lfcShrinkage()` part. Second argument should be `coefs`, not `contrasts`. Maybe if I fix this part it will work with the latest version.

## Implement the option to use Limma instead of DESEQ2

## Some kind of comparation or Benchmark
