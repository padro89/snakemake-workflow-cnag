
def get_strand_param(wildcards):

    if config['strand'] == 'forw':
        return '--forward-prob 1'
    elif config['strand'] == 'rev':
        return '--forward-prob 0'
    else:
        return '--forward-prob 0.5'

def get_paired_param(wildcards):
    param = ''
    if config['paired']:
        param = '--paired-end'
    return param

def get_prefix_quant(wildcards):
    return config['path']['quant']+"/"+wildcards.FLI

rule quantification:
    input:
        config['path']['map']+"/{FLI}.Aligned.Transcriptome.out.bam"
    
    conda:
        "env/rsem."+config["versions"]['rsem_version']+".yaml"
    
    params:
        strand = get_strand_param,
        paired = get_paired_param,
        basename = get_prefix_quant,
        ref = config["index"]["rsem"]   
    threads: 8
    output:
        config['path']['quant']+"/{FLI}.genes.results",
        config['path']['quant']+"/{FLI}.isoforms.results"

    shell:
        "rsem-calculate-expression -p {threads} "
        "--no-bam-output --bam "
        "{params.strand} "
        "--temporary-folder $TMPDIR/rsem  "
        "{params.paired} "
        "{input} "
        "{params.ref} "
        "{params.basename}"

rule generate_table_genes:
    input:
        expand(config['path']['quant']+"/{FLI}.genes.results", FLI=FLI)

    output:
        config['project']+".genes.tsv"

    params:
        gtf = config["index"]['gtf']

    conda:
        "env/python.yaml"
    group: "table_counts"    
    shell:
        "scripts/gtf2genes.py -t counts "
        "--field gene "
        "-g {params.gtf} "
        "--feature gene_name gene_type "
        "-f {input} > {output}"
    

rule generate_table_isoforms:
    input:
        expand(config['path']['quant']+"/{FLI}.isoforms.results", FLI=FLI)

    output:
        config['project']+".isoforms.tsv"

    params:
        gtf = config["index"]['gtf']

    conda:
        "env/python.yaml"
    group: "table_counts"    
    
    shell:
        "scripts/gtf2genes.py -t counts "
        "--field transcript "
        "-g {params.gtf} "
        "--feature gene_name transcript_type "
        "-f {input} > {output}"