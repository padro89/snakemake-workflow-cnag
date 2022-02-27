# Transcriptomics workflow

Repository used to create a Snakemake transcriptomics workflow for CNAG

In this file I will reflect any ideas I have concerning this project, and what I am using. Eventually it will become a description of the process and it will explain how to use, functioning as the pipeline ducumentation.

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
- `deseq2-script`: A branch where I try to implement the deseq2.R script in Snakemake.
- `update-deseq2-version`: A branch where I will update the DESeq2 version in the script.

## What I'd like to do

### Install conda/mamba and create an environment to run the deseq2 script

- I created a .yaml file `mamba-envir.yaml` and environment adapted to de DESeq2 version (1.26). Conda (or Mamba) automatically resolved other packages versions.
- The .yaml file had several fields. `name`, `channels` (where I specified conda-forge and bioconda) and `dependencies`, where I specified the packages to install into the new environment.
- To create the environment, I ran `mamba env create -f <yaml_file>`.
- I ran the deseq2.R script using Rscript and the arguments, and it worked flawlesly without any need of adding any packages.
- If I need to  manually add a package, I'll add it to my conda environment yaml file.

### Implement the deseq2 script

- I used the parameters defined in the config.R file in the config.yaml Snakemake file.
- All the arguments were included in the config.ymal.
- I implemented the R script using S4 objects in the script to make calls to Snakemake config.yaml.
- Should I change anything about the output?

### Modularize the deseq2 script inside snakemake

- Ideally, one script per step.
- Snakemake S4 objects should be used inside the scripts to get results from scripts run in other rules.
- PCA should be shown first, and only afterwards should one run the complete Snakemake.
- Output should be generated in Snakemake S4 objects so they can be properly organized.

### Implement the pathway analysis using Snakemake

### Update the  script to include DESeq latest version

- There is some problem with the `lfcShrinkage()` part. Second argument should be `coefs`, not `contrasts`. Maybe if I fix this part it will work with the latest version.
- I don't need a specific environment to test this, I can just run and correct the script in my system deactivating conda/mamba.

### Change plots in the output

- Add volcano plots.
- Modify the PCA using my own function.

### Implement the option to use limma instead of DESeq2

- It should be specified in the config.yaml whether DESeq2 or limma-voom should be used.

### Implement the option to create a "Materials and Methods" file

- It could export software versions.
- It could generate the relevant bibliography in a .bib file.

### Some kind of comparation or Benchmark

### Implement a different README.md file explaining how everything works.
