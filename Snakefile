# Loading config file

configfile: "config.yaml"

# Functions to get output and some parameters
def get_contrast(wildcards):
    return config["contrasts"][wildcards.contrast]

def get_pca_output(wildcards):
    pca_output = []
    if config["dge_method"] == "deseq2":
        pca_output = rules.pca_deseq2.output
    if config["dge_method"] == "limma":
        pca_output = rules.pca_limma.output
    return pca_output

def get_deseq2_output(wildcards):
    deseq2_output = []
    if config["onlypca"] == False:
        if config["dge_method"] == "deseq2":
            deseq2_output = expand(rules.dge_deseq2.output, contrast = config["contrasts"])
    return deseq2_output

def get_limma_output(wildcards):
    limma_output = []
    if config["onlypca"] == False:
        if config["dge_method"] == "limma":
            limma_output = expand(rules.dge_limma.output, contrast = config["contrasts"])
    return limma_output

def get_ora_output(wildcards):
    ora_output = []
    if config["onlypca"] == False:
        if config["pathway_method"] == "both" or config["pathway_method"] == "gprofiler":
            ora_output = expand(rules.ora.output, contrast = config["contrasts"])
    return ora_output

def get_fgsea_output(wildcards):
    fgsea_output = []
    if config["onlypca"] == False:
        if config["pathway_method"] == "both" or config["pathway_method"] == "fgsea":
            fgsea_output = expand(rules.fgsea.output, contrast = config["contrasts"])
    return fgsea_output    

# Rules to include
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


# Target rule
rule all:
    input:
        get_pca_output,
        get_deseq2_output,
        get_limma_output,
        get_ora_output,
        get_fgsea_output

# Other rules (to be put in other files)
if config["gmt"] is not None:
    rule fgsea:
        input:
            deg_results = config["path"]["dge"]+"/{contrast}/{contrast}_deg_results.txt",
            gmt = config["gmt"],
        output:
            fgsea_pdf = config["path"]["dge"]+"/{contrast}/{contrast}_fgsea.pdf",
            fgsea_main_pathways_pdf = config["path"]["dge"]+"/{contrast}/{contrast}_fgsea_main_pathways.pdf",
            fgsea = config["path"]["dge"]+"/{contrast}/{contrast}_fgsea.tsv"
        script:
            "scripts/fgsea.R"

if config["gmt"] is None:
    from datetime import datetime
    current_date = datetime.now().strftime("%Y_%m_%d")

    rule get_reactome_pathways:
        output:
            reactome = temp(config["path"]["dge"]+"/ensemble2reactome.txt")
        shell:
            "echo 'Downloading patwhays from reactome, it may take a while' &&"
            "GMT_URL=https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt &&"
            "wget $GMT_URL -O {output} && unset $GMT_URL"

    rule make_gmt:
        input:
            reactome = config["path"]["dge"]+"/ensemble2reactome.txt"
        output:
            gmt = config["path"]["dge"]+f"/ensemble2reactome_{current_date}.gmt"
        script:
            "scripts/create_gmt.py"

    rule fgsea:
        input:
            deg_results = config["path"]["dge"]+"/{contrast}/{contrast}_deg_results.txt",
            gmt = config["path"]["dge"]+f"/ensemble2reactome_{current_date}.gmt"
        output:
            fgsea_pdf = config["path"]["dge"]+"/{contrast}/{contrast}_fgsea.pdf",
            fgsea = config["path"]["dge"]+"/{contrast}/{contrast}_fgsea.tsv"
        script:
            "scripts/fgsea.R"    

rule ora:
    input:
        deg_list = config["path"]["dge"]+"/{contrast}/{contrast}_deg_list.txt"
    output:
        ora = config["path"]["dge"]+"/{contrast}/{contrast}_ora.table"
    script:
        "scripts/gprofiler.R"

# Selecting the DEGs.
if config["dge_method"] == "deseq2":
    rule select_deg_deseq2:
        input:
            deg_results = config["path"]["dge"]+"/{contrast}/{contrast}_deg_results.txt"
        output:
            deg_list = temp(config["path"]["dge"]+"/{contrast}/{contrast}_deg_list.txt"),
        shell:
            "cat {input} | awk '{{if($NF < 0.05 && sqrt($4^2) > log(1.5)/log(2))"
            " print $1}}' > {output}"

if config["dge_method"] == "limma":
    rule select_deg_limma:
        input:
            deg_results = config["path"]["dge"]+"/{contrast}/{contrast}_deg_results.txt"
        output:
            deg_list = temp(config["path"]["dge"]+"/{contrast}/{contrast}_deg_list.txt"),
        shell:
            "cat {input} | awk '{{if($6 < 0.05 && sqrt($2^2) > log(1.5)/log(2))"
            " print $1}}' > {output}"


# Running the PCA 

if config["dge_method"] == "deseq2":
    rule pca_deseq2:
        input:
            counts = config["counts"],
            info = config["info"]
        output:
            pcas = config["path"]["dge"]+"/PCA.pdf",
            sampletosample = config["path"]["dge"]+"/sampletosample_heatmap.pdf",
            rlog = config["path"]["dge"]+"/rlogMat.txt",
            pc_contribution = config["path"]["dge"]+"/pc_contribution.txt",
            dds = temp(config["path"]["dge"]+"/dds"),
            norm_counts = config["path"]["dge"]+"/norm_counts.txt"
        script:
            "scripts/PCA_deseq2.R"

    rule dge_deseq2: 
        input:
            dds = config["path"]["dge"]+"/dds",
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

if config["dge_method"] == "limma":
    rule pca_limma:
        input:
            counts = config["counts"],
            info = config["info"]
        output:
            sampletosample = config["path"]["dge"]+"/sampletosample_heatmap.pdf",
            pcas = config["path"]["dge"]+"/PCA.pdf",
            pc_contribution = config["path"]["dge"]+"/pc_contribution.txt",
            dge = temp(config["path"]["dge"]+"/dge")
        script:
            "scripts/PCA_limma.R"

    rule dge_limma:
        input:
            dge = config["path"]["dge"]+"/dge"
        output:
            expand(config["path"]["dge"]+"/{contrast}/{contrast}_deg_results.txt",
            contrast = config["contrasts"]),
            expand(config["path"]["dge"]+"/{contrast}/{contrast}_heatmap.pdf",
            contrast = config["contrasts"]),
            logcpm = config["path"]["dge"]+"/logCPM.txt"
        script:
            "scripts/limma_voom.R"
