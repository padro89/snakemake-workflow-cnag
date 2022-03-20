def get_contrast(wildcards):
    return config["contrasts"][wildcards.contrast]

def get_final_output(wildcards):
    pass

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
        stats = expand("/home/joan/Test/Results/{contrast}_stats.txt",
            contrast=config["contrasts"]),
        norm_counts = expand("/home/joan/Test/Results/{contrast}_norm_counts.txt",
            contrast=config["contrasts"]),
        deg_results = expand("/home/joan/Test/Results/{contrast}_deg_results.txt",
            contrast=config["contrasts"]),
        topDEgenes_heatmap = expand("/home/joan/Test/Results/{contrast}_top50DEgenes_heatmap.pdf",
            contrast=config["contrasts"]),
        heatmap_custom = expand("/home/joan/Test/Results/{contrast}_topgenes_heatmap.tiff",
            contrast=config["contrasts"])

# Running de DGE analysis with DESeq2

rule dge: 
    input:
        dds_design = "/home/joan/Test/Results/dds_design",
        rlogmat = "/home/joan/Test/Results/rlogMat.txt",
    output:
        stats = "/home/joan/Test/Results/{contrast}_stats.txt",
        norm_counts = "/home/joan/Test/Results/{contrast}_norm_counts.txt",
        deg_results = "/home/joan/Test/Results/{contrast}_deg_results.txt",
        topDEgenes_heatmap = "/home/joan/Test/Results/{contrast}_top50DEgenes_heatmap.pdf",
        heatmap_custom = "/home/joan/Test/Results/{contrast}_topgenes_heatmap.tiff",
    params: 
        contrast = get_contrast,
    script:
        "scripts/deseq2.R"

# Calculating dispersion with the design

rule design:
    input:
        dds = "/home/joan/Test/Results/dds"
    output:
        dds_design = "/home/joan/Test/Results/dds_design"
    script:
        "scripts/deseq2_design.R"

# Running the PCA 
rule PCA:
    input:
        counts = "../testdata/counts_3",
        info = "../testdata/info_3_several.txt"
    # These directories should make use of wildcards from the config file.
    output:
        pcas = "/home/joan/Test/Results/PCA.pdf",
        sampletosample = "/home/joan/Test/Results/sampletosample_heatmap.pdf",
        rlog = "/home/joan/Test/Results/rlogMat.txt",
        pc_contribution = "/home/joan/Test/Results/pc_contribution.txt",
        dds = temp("/home/joan/Test/Results/dds")
    script:
        "scripts/PCA.R"