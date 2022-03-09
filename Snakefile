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

rule all:
    input:
        "../testdata/counts_3",
        "../testdata/info_3.txt"
    
# Running the PCA
rule PCA:
    input:
        counts = "../testdata/counts_3",
        info = "../testdata/info_3.txt"
    # These directories should make use of wildcards from the config file.
    output:
        pcas = "/home/joan/Test/PCA.pdf",
        sampletosample = "/home/joan/Test/Results/sampletosample_heatmap.pdf",
        rlog = "/home/joan/Test/Results/rlogMat.txt",
        pc_contribution = "/home/joan/Test/Results/pc_contribution.txt",
        temp("/home/joan/Test/Results/dds.Rds")
    script:
        "scripts/PCA.R"
rule dge: 
    input:
        dds = "/home/joan/Test/Results/dds.Rds"
    script:
        "scripts/deseq2.R"