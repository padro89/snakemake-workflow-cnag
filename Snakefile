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

# Running de DGE analysis with DESeq2
rule dge: 
    input:
        dds = "/home/joan/Test/Results/dds",
        rlogmat = "/home/joan/Test/Results/rlogMat.txt"
    output:
        stats = "/home/joan/Test/Results/stats.txt"
    script:
        "scripts/deseq2.R"
    
# Running the PCA 
rule PCA:
    input:
        counts = "../testdata/counts_3",
        info = "../testdata/info_3.txt"
    # These directories should make use of wildcards from the config file.
    output:
        pcas = "/home/joan/Test/Results/PCA.pdf",
        sampletosample = "/home/joan/Test/Results/sampletosample_heatmap.pdf",
        rlog = "/home/joan/Test/Results/rlogMat.txt",
        pc_contribution = "/home/joan/Test/Results/pc_contribution.txt",
        dds = "/home/joan/Test/Results/dds"
    script:
        "scripts/PCA.R"