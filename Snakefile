# Mapping workflow
"""# Including rules
include: "rules/common.smk"
include: "rules/trim.smk"
include: "rules/map.smk"
include: "rules/quant.smk"

# Setting targets 

rule all:
    input:
        expand(config['path']['map']+"/{FLI}.Aligned.sortedByCoord.out.bam",FLI=FLI),
        expand(config['path']['quant']+"/{FLI}.genes.results",FLI=FLI),
        expand(config['path']['quant']+"/{FLI}.isoforms.results",FLI=FLI),
        config['project']+".isoforms.tsv",
        config['project']+".genes.tsv"""

# Loading config file

configfile: "config.yaml"

rule dge:
    input:
        "../testdata/counts_3",
        "../testdata/info_3.txt"
    script:
        "scripts/deseq2.R"
    
# Running the PCA
rule PCA:
    input:
        "../testdata/counts_3",
        "../testdata/info_3.txt"
    output:
        countdata = output.countdata,
        coldata = output.coldata
# Executing the Deseq2 script
    script:
        "scripts/PCA.R"