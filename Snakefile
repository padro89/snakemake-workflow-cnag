def get_contrast(wildcards):
    return config["contrasts"][wildcards.contrast]

def get_pca_output(wildcards):
    return rules.PCA.output

def get_deseq_output(wildcards):
    # mapping_output = []
    # dge_limma_output = [expand(config["path"]["dge"]+"/{contrast}/{contrast}_stats.txt",
    #                         contrast=config["contrasts"]),
    #                     config["path"]["dge"]+"/norm_counts.txt",
    #                     expand(config["path"]["dge"]+"/{contrast}/{contrast}_deg_results.txt",
    #                         contrast=config["contrasts"]),
    #                     expand(config["path"]["dge"]+"/{contrast}/{contrast}_top50DEgenes_heatmap.pdf",
    #                         contrast=config["contrasts"]),
    #                     expand(config["path"]["dge"]+"/{contrast}/{contrast}_topgenes_heatmap.tiff",
    #                         contrast=config["contrasts"]),
    #                     expand(config["path"]["dge"]+"/{contrast}/{contrast}_deg_list.txt",
    #                         contrast=config["contrasts"])]
    # ora_output = [expand(config["path"]["dge"]+"/{contrast}/{contrast}_ora.table",
    #                 contrast=config["contrasts"])]
    # final_output.append(mapping_output)
    # final_output.append(dge_limma_output)
    # final_output = rules.PCA.output
    if config["onlypca"] == False:
        deseq_output = expand(rules.ora.output, contrast = config["contrasts"])
    else:
        deseq_output = []
    #final_output.append(rules.select_deg.output)
    return deseq_output

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
        get_pca_output,
        get_deseq_output,
        #get_limma_output
        #expand(config["path"]["dge"]+"/{contrast}/{contrast}_ora.table",
        #   contrast=config["contrasts"])

# Running the GSEA



# Running the ORA

rule ora:
    input:
        deg_list = config["path"]["dge"]+"/{contrast}/{contrast}_deg_list.txt"
    output:
        ora = config["path"]["dge"]+"/{contrast}/{contrast}_ora.table"
    script:
        "scripts/gprofiler.R"

# Selecting the DEGs.
rule select_deg:
    input:
        deg_results = config["path"]["dge"]+"/{contrast}/{contrast}_deg_results.txt"
    output:
        deg_list = temp(config["path"]["dge"]+"/{contrast}/{contrast}_deg_list.txt"),
    shell:
        "cat {input} | awk '{{if($NF < 0.05 && sqrt($4^2) > log(1.5)/log(2))"
        " print $1}}' | cut -d ',' -f 1 | cut -d '.' -f 1 > {output}"

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
        info = "../testdata/info_3_continua.txt"
    # These directories should make use of wildcards from the config file.
    output:
        pcas = config["path"]["dge"]+"/PCA.pdf",
        sampletosample = config["path"]["dge"]+"/sampletosample_heatmap.pdf",
        rlog = config["path"]["dge"]+"/rlogMat.txt",
        pc_contribution = config["path"]["dge"]+"/pc_contribution.txt",
        dds = temp(config["path"]["dge"]+"/dds")
    script:
        "scripts/PCA.R"