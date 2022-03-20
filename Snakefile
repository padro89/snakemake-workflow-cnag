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
        stats = expand(config["path"]["dge"]+"/{contrast}/{contrast}_stats.txt",
            contrast=config["contrasts"]),
        norm_counts = config["path"]["dge"]+"/norm_counts.txt",
        deg_results = expand(config["path"]["dge"]+"/{contrast}/{contrast}_deg_results.txt",
            contrast=config["contrasts"]),
        topDEgenes_heatmap = expand(config["path"]["dge"]+"/{contrast}/{contrast}_top50DEgenes_heatmap.pdf",
            contrast=config["contrasts"]),
        heatmap_custom = expand(config["path"]["dge"]+"/{contrast}/{contrast}_topgenes_heatmap.tiff",
            contrast=config["contrasts"])

# Running de DGE analysis with DESeq2

rule dge: 
    input:
        dds_design = config["path"]["dge"]+"/dds_design",
        rlogmat = config["path"]["dge"]+"/rlogMat.txt",
    output:
        stats = config["path"]["dge"]+"/{contrast}/{contrast}_stats.txt",
        deg_results = config["path"]["dge"]+"/{contrast}/{contrast}_deg_results.txt",
        topDEgenes_heatmap = config["path"]["dge"]+"/{contrast}/{contrast}_top50DEgenes_heatmap.pdf",
        heatmap_custom = config["path"]["dge"]+"/{contrast}/{contrast}_topgenes_heatmap.tiff",
    params: 
        contrast = get_contrast,
    script:
        "scripts/deseq2.R"

# Calculating dispersion with the design

rule design:
    input:
        dds = config["path"]["dge"]+"/dds"
    output:
        dds_design = temp(config["path"]["dge"]+"/dds_design"),
        norm_counts = config["path"]["dge"]+"/norm_counts.txt"
    script:
        "scripts/deseq2_design.R"

# Running the PCA 
rule PCA:
    input:
        counts = "../testdata/counts_3",
        info = "../testdata/info_3_several.txt"
    # These directories should make use of wildcards from the config file.
    output:
        pcas = config["path"]["dge"]+"/PCA.pdf",
        sampletosample = config["path"]["dge"]+"/sampletosample_heatmap.pdf",
        rlog = config["path"]["dge"]+"/rlogMat.txt",
        pc_contribution = config["path"]["dge"]+"/pc_contribution.txt",
        dds = temp(config["path"]["dge"]+"/dds")
    script:
        "scripts/PCA.R"