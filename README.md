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
- `deseq2-modular`: A branch where I try to modularize the DESeq2 script.
- `pathway_analysis`: A branch where I will try to create a pathway analysis.


## Sources:

I found this link where a similar analysis is run, and where I can get ideas:

https://github.com/snakemake-workflows/rna-seq-star-deseq2

## What I'd like to do

### Install conda/mamba and create an environment to run the deseq2 script

- I created a .yaml file `mamba-envir.yaml` and environment adapted to de DESeq2 version (1.26). Conda (or Mamba) automatically resolved other packages versions.
- The .yaml file had several fields. `name`, `channels` (where I specified conda-forge and bioconda) and `dependencies`, where I specified the packages to install into the new environment.
- To create the environment, I ran `mamba env create -f <yaml_file>`.
- I ran the deseq2.R script using Rscript and the arguments, and it worked flawlesly without any need of adding any packages.
- If I need to  manually add a package, I'll add it to my conda environment yaml file.
- After updating the deseq2 script I created a new environment with the up-to-date packages.

### Implement the deseq2 script

- I ported the parameters defined in the config.R file in the config.yaml Snakemake file.
- All the arguments were included in the config.ymal.
- I implemented the R script using S4 objects in the script to make calls to Snakemake config.yaml and to import the data from the input in the rule.
- Should I change anything about the output?
- Continuous variables are converted to factors. They must be specified in the config file in the `continuous` argument, along with the number of groups to do in the `n_continuous splits` package. I made a little loop in the PCA.R script that then converts them to factors:

```
if (!is.null(continuous)){
  for(i in seq(length(continuous))){
    levels_temp <- paste(continuous[i],seq(n_continuous_splits[i]),sep="_")
    coldata[,continuous[i]] <- cut(coldata[,continuous[i]],
                                  breaks = n_continuous_splits[i],
                                  labels=levels_temp)
 }
}
```

### Modularize the deseq2 script inside snakemake

#### Saving the output

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

#### Saving temporary files

- I should try to use the same method to save objects in temporary locations like the dds object. Teoretically, I should just:

In R
```
save(dds, file = snakemake@output$dds)
```

In Snakemake
```
rule A:
  input: dds = "/some/location/dds"
  script: someRscript.R

rule B:
  input: "/some/input"
  output:
    temp(dds = "/some/location/dds")

## This doesn't work either
rule B:
  input: "/some/input"
  output:
    temp(dds = "/some/location/dds")
```

~For some reason this doesn't work. The `temp()` function interprets de `dds =` as an argument. I don't know how to use the output inside this function.~ The correct syntax in the output is `dds = temp("/file/location")`.

#### Starting modularization

##### PCA

- PCA should be shown first, and only afterwards should one run the complete Snakemake. In order to do this, I created a different script, called `PCA.R`.
- This script creates the dds file with a null formula to create the PCA from it, and the rest of scripts just import it (or update the formula). **A disadvantage of this is that I may have to load DESeq2 library many times**.
- When running the PCA, if `pca_atr` is blank it should run groups by default.
```
if(is.null(pca_atr) | pca_atr ! %in% colnames(dds@coldata)){
    pca_atr=="group"
}
```
- I should be able to specify several variables for the PCA so it runs the PCA coloring for each variable FROM THE START. Currently, if several variables are specified, it creates subgroups.
- I have to update the libraries needed in the PCA script so it takes less time to charge.
- It works like this: First, a PCA is run, generating the graphs and the rlogMat.txt file, as well as other relevant files and the dds.
- PCA .tiff graph is not created. There is a piece of code to create it I just have to uncomment and add to the output.
- In the script, the interest variable name is changed to "group". This causes that if one specifies the plot attribute to be the name of the interest variable instead of "group", it cannot do the PCA. To solve this, I put a conditional statement when loading the snakemake object:
  
```
if(snakemake@config$plot_atr$pca == snakemake@config$group |
   is.null(snakemake@config$plot_atr$pca)){
  pca_atr <- "group"
}else{
  pca_atr <- snakemake@config$plot_atr$pca
}
```

##### Updating the design formula.

- A script called `deseq2_design.R` is used for this.
- The dds file is passed as input and it runs `DESeq()` function using the specified formula. 
- It creates a new temporary file called `dds_design`.

##### Running the DGEA

- A third script, simply called `deseq2.R` uses the `dds_design` object and the `rlogMat.txt` file.
- In order to run multiple comparisions, the contrasts must be specified in the `config.yaml` like this:

```
contrasts:
  contrast1:
    - Group 1
    - Group 2
```
- Then I can use each contrast in the script by using params:

```
rule dge:
  input: 
    "/some/input",
  output:
    "/some/{contrast}_output",
  params:
    contrast = get_contrast
  script:
    "some_script.R"
  ```

- To use the contrasts in the params, I created a specific function that gets each contrast:

```
def get_contrast(wildcards):
  return config["contrasts"][wildcards.contrast]
```
- The problem with wildcards is, if I want to check if it works, I must run the all rule with the wildcards in it, or otherwise it won't work:

```
rule all:
  input:
    for_every_file = expand("folder/{contrast}_file",
                            contrast = config["contrasts"])
```
- This could be a problem to run only the PCA, but for the moment it is not because it has no wildcards. If I put wildcards in it later, I could do something like:

```
if config["only_pca"] == True:
  rule all:
    only the PCA files
```

- Everytime a comparation is performed, a new normalized counts file is created and saved. I cannot output some files with wildcards and others without, so I should create a new script that only saves the normalized counts (sorted).
- I created an if statement inside the deseq2.R script that should create a heatmap even if there are not enough differentially expressed genes according to de_genes_n, but at least 25 DEG are present. I should test if it works.
- I just realized that in the forumla the group variable must always be "group", as the name is changed in the first script. I should really try to simplify the script.

### Implement the pathway analysis using Snakemake

#### ORA

- The pathway analysis is performed with g:Profiler.
- For the DGE list, I take the genes that are differentialy expressed (they may have a column "Filter" if the analysis is run with DESeq2).
- I then get the ENSEMBL ID. It may have a last number after a point I should remove, this happens with *Mus musculus* and *Homo sapiens* samples. I can use awk for this.
- To do this, I use a new rule that, using awk, selects the genes with `Filter==1` and then removes everything before the comma with `cut -d "," -f 1` (-d specifies the delimiter, -f the column) and everything before the point (in case its from mouse or human).
- To avoid using the filter column (which is not present always), I can create the same awk command using the last column (adj p. value) and the shrunken logFC. To do this, I use an if statement comparing the adjusted p-value and an absolute shrunken-logFC (which I get with the square root of the square of the column) greater than log2(1.5). awk only provides natural logarithm, but I have to consider that:

log2(x) = ln(x)/ln(2) (base change rule of logarithms).

The resulting command is: 
```
awk '{{if($NF<0.05 & sqrt($4^2) > log(1.5)/log(2)) print $1}}'
```
I have tested this result many times, and it gives the same as the filter column.

Note: I MUST remember that when using awk inside snakemake, I must use double brackets{{}}.

- Maybe it would be more flexible to just transform the results inside an R script, where I can use conditional statements more easily.
- There is no need to select a Universe, g:Profiler uses one automatically.

#### GSEAs
- I need to indicate the column that contains the statistic in the DGE list in order to rank the genes. They are using the shrunkenlFC. Should I use the statistic?
- If instead of ENSEMBL, GENEID is used, the split should be done always before the comma. If human or mouse is used with ENSEMBL ID, the split should be before the point. I solved this with this statement inside the fgsea script:

```
# Getting the split character:
# Defaults to comma
split <- "\\,"
# If ENSEMBLEID is used
if(colid == 1){
  # If human or mouse
  if(snakemake@config$species == "hsapiens" |
     snakemake@config$species == "mmusculus"){
    # Sets it to point
    split <-"\\."
  }
}
```
- I can make changes to avoid the COLID variable, which states which of the names should be used (ENSEMBL, GENEID), to something more intuitive in the config file like identifier.
- With the ENSEMBL ID I can use a GMT file that contains all the pathways for the organism.
- Maybe it would be interesting to create some graphs, like in the GSEA, like those in the REVIGO.


Note: IF THERE ARE NO DEGs, I SHOULD CREATE AN EMPTY FILE TO AVOID SNAKEMAKE FAILURE LACKING OUTPUT.

### Integrate Marc's alignment workflow

- He uses a FLI file with the samples. I'm not sure if I can include it in the config file, but it would probably be more useful as an input file.
- He told me it does not work. I should revise it. I can connect to the CNAG cluster and try it in interactive mode. I created aliases. Once in the cluster I can enter interactive mode with `mnsh -C 8 -x`.

### Integrate everything toghether
- I should probably make subworkflows for different options in the config.
- I can get the output for rule all with different methods. One is to define a function that gets it. I can use something like this:

```
def get_final_output(wildcards):
    deseq2_output = expand(rules.<rule>.output, contrast = config["contrasts"]
    pca_output = rules.PCA.output

# Then try something like
    final_output = []
    if config["some_option"] == value:
        final_output.append(deseq2_output)

    return final_output

# And call the function in the final rule.
rule all:
    input:
        get_final_output
```
**Note:** This does not work. I have changed it in the code so it works. 
- I can use Python `if()` inside any rule to change the workflow according to specified parameters.

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
- I can include a conditional statement in the R scripts depending on the config file deciding wether to make tiff output. This would give an error whenever the output is missing. To avoid this, I can make a python function to get the snakemake output including the same conditional statements.

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