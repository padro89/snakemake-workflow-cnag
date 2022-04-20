def get_contrast(wildcards):
    return config["contrasts"][wildcards.contrast]

def get_pca_output(wildcards):
    return rules.PCA.output

def get_deseq_output(wildcards):
    if config["onlypca"] == False:
        deseq_output = expand(rules.ora.output, contrast = config["contrasts"])
    else:
        deseq_output = []
    #final_output.append(rules.select_deg.output)
    return deseq_output

def get_limma_output(wildcards):
    return rules.limma.output

def get_fgsea_output(wildcards):
    if config["onlypca"] == False:
        fgsea_output = expand(rules.fgsea.output, contrast = config["contrasts"])
    else:
        fgsea_output = []
    return fgsea_output

# Mapping workflows
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
        #get_deseq_output,
        #get_fgsea_output,
        get_limma_output

# Running the GSEA

rule fgsea:
    input:
        deg_results = config["path"]["dge"]+"/{contrast}/{contrast}_deg_results.txt",
        gmt = config["gmt"],
    output:
        fgsea_pdf = config["path"]["dge"]+"/{contrast}/{contrast}_fgsea.pdf",
        fgsea = config["path"]["dge"]+"/{contrast}/{contrast}_fgsea.tsv"
    script:
        "scripts/fgsea.R"


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

rule limma:
    input:
        counts = "../testdata/counts_3",
        info = "../testdata/info_3_continua.txt"
    output:
        mds = config["path"]["dge"]+"/MDS.pdf",
    #params:
    #    contrast = get_contrast,
    script:
        "scripts/limma_voom_basic.R"
