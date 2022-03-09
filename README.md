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
- `deseq2-version`: A branch where I will update the DESeq2 version in the script.

## What I'd like to do

### Install conda/mamba and create an environment to run the deseq2 script

- I created a .yaml file `mamba-envir.yaml` and environment adapted to de DESeq2 version (1.26). Conda (or Mamba) automatically resolved other packages versions.
- The .yaml file had several fields. `name`, `channels` (where I specified conda-forge and bioconda) and `dependencies`, where I specified the packages to install into the new environment.
- To create the environment, I ran `mamba env create -f <yaml_file>`.
- I ran the deseq2.R script using Rscript and the arguments, and it worked flawlesly without any need of adding any packages.
- If I need to  manually add a package, I'll add it to my conda environment yaml file.
- After updating the deseq2 script I created a new environment with the up-to-date packages.

### Implement the deseq2 script

- I used the parameters defined in the config.R file in the config.yaml Snakemake file.
- All the arguments were included in the config.ymal.
- I implemented the R script using S4 objects in the script to make calls to Snakemake config.yaml and to import the data from the input in the rule.
- Should I change anything about the output?
- I MUST ADAPT THE PROCESSING OF THE CONTINUOUS VARIABLES, BUT I HAVE NONE. I SHOULD ASK FOR ONE EXAMPLE SO I GET WHAT THEY DO. They are probably adapted for graphics, so I could make a random continuous variable in the coldata and run the script to see where it fails and fix it.

### Modularize the deseq2 script inside snakemake

- I can save anything, for example, a PDF, as `snakemake@output$name`, and then call it from the output rule:

Script:
```
pdf(snakemake@output$result)
do_something()
dev.off()
```

Snakemake rule:
```
    output:
        result = "some/location/file.pdf"
```
And it works. Now I just need to implement every output like this. 

- ~If I run this with the output created in the snakemake folder, everything works fine. If I don't, somehow it creates the files but then says the files are missing.~ It turns out that the only problem was writing "~/". Full paths must be used.
- ~For some reason, this didn't work with the function `write.table()`~. It does, I just have to specify the argument `file=`.
- I should try to use the same method to save objects in temporary locations like the dds object.
- If this is possible, I would like to create a first script that creates a dds file, and the rest of scripts just import it (or update the formula). **A disadvantage of this is that I may have to load DESeq2 library many times**.
- PCA should be shown first, and only afterwards should one run the complete Snakemake.
- When running the PCA, if `pca_atr` is blank it should run groups by default.
```
if(pca_atr==NULL){
    pca_atr=="group"
}
```
- I have to update the libraries needed in the PCA script so it takes less time to charge.

### Implement the pathway analysis using Snakemake

- The pathway analysis is performed with g:Profiler.
- Both ORA and GSEA are carried out, I have a script for each of them.
- ORA: 
  - For the DGE list, I take the genes that are differentialy expressed (they may have a column "Filter" if the analysis is run with DESeq2).
  - I then get the ENSEMBL ID. It may have a last number after a point I should remove, this happens with *Mus musculus* and *Homo sapiens* samples. I can use awk for this.
  - With the ENSEMBL ID I can use a GMT file that contains all the pathways for the organism.
  - There is no need to select a Universe, g:Profiler uses one automatically.
- GSEA:
  - I need to indicate the column that contains the statistic in the DGE list in order to rank the genes. 
  - Maybe it would be interesting to create some graphs, like in the GSEA, like those in the REVIGO.

### Integrate Marc's alignment workflow

- He uses a FLI file with the samples. I'm not sure if I can include it in the config file, but it would probably be more useful as an input file.

### Integrate everything toghether
- I should probably make subworkflows for different options in the config.
- Should I use the config for everything, or should I use params?

### Add conda dependencies to the each step in the Snakemake workflow

- There's a specific option in each rule.

### Update the  script to include DESeq latest version

- There is some problem with the `lfcShrinkage()` part. Second argument should be `coefs`, not `contrasts` if "apeglem" is to be used. "normal" should work with `contrasts`, as well as "ashr", and should be specified in the argument `type` of the funciton.
- When I tried to use the `coefs`argument it failed.
- This was because there was an error when releveling the group in the DDS object. "" were missing, so instead of releveling the "group" variable, it created a new one using the group variable, which then releveled, leaving the "group" variable to be converted to a factor automatically by DESeq2 when the function was run. This caused that the names of the coefs in `resultsNames(dds)` did not match the name of the coef to use. I fixed this by adding "group".
- Using "apeglm" with `coefs` did not allow for different comparisions, only for the ones in `resultsNames(dds)`, otherwise the dds object has to be created with the new reference level of interest (which is not very practical).
- I have found a post discussing this (https://support.bioconductor.org/p/123247/) and I can use something like this to avoid running DSeq2 again:
```
dds <- makeExampleDESeqDataSet()
dds$condition <- factor(rep(1:3,each=4)) # suppose 3 groups
dds <- DESeq(dds)
resultsNames(dds)
[1] "Intercept"        "condition_2_vs_1" "condition_3_vs_1"
dds$condition <- relevel(dds$condition, "2")
dds <- nbinomWaldTest(dds)
resultsNames(dds)
[1] "Intercept"        "condition_1_vs_2" "condition_3_vs_2"
lfc <- lfcShrink(dds, coef="condition_3_vs_2", type="apeglm")
```
- "ashr" should be just fine, so I added an option to use either "normal" or "ashr" in the `config.yaml`. If I can implement "apeglm" I should make an additional option.
- It now works with the latest version.

### Change plots in the output

- Add volcano plots.
- Modify the PCA using my own function.
- Add the option to choose any variables for coloring (including continuous variables, that should be factorized for this).

### Implement the option to use limma instead of DESeq2

- It should be specified in the config.yaml whether DESeq2 or limma-voom should be used.

### Implement the option to create a "Materials and Methods" file

- It could export software versions. `sink(sessionInfo)`
- It could generate the relevant bibliography in a .bib file.

### Some kind of comparation or Benchmark

- Maybe I can compare DESeq2 and limma, or GSEA and ORA.
- Maybe I can compare different programs.

### Implement a different README.md file explaining how everything works.

- This file, but for the moment it only contains everything I do and learn.