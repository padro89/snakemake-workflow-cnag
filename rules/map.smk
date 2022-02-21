#Functions
def get_map_params(wildcards):
    return ' '.join(map_config[config['map_mode']])

def set_prefix_map(wildcards):
    return config['path']['map']+"/"+wildcards.FLI+"."

def get_fq_map(wildcards):
    fq_a = [meta_seq.F1.loc[wildcards.FLI]]
    if config['paired']:
        fq_a.append(meta_seq.F2.loc[wildcards.FLI])
    if config['trimming']:
        fq_a = [config['path']['trim']+"/{FLI}_trimmed.fq.gz"]
        if config['paired']:
            fq_a = [
                config['path']['trim']+"/{FLI}_val_1.fq.gz",
                config['path']['trim']+"/{FLI}_val_2.fq.gz"
            ]
    return fq_a
    
# Rules
rule mappings:
    input:
        get_fq_map

    output: 
        gbam = config['path']['map']+"/{FLI}.Aligned.sortedByCoord.out.bam",
        tbam = config['path']['map']+"/{FLI}.Aligned.Transcriptome.out.bam"

    threads: 8
    
    params: 
        time = "02:00:00",
        map_params = get_map_params,
        ref = config["index"]['star'],
        FLI = get_fli,
        Barcode = get_barcode,
        prefix = set_prefix_map

    conda:
        "env/star."+config["versions"]['star_version']+".yml"
    shell:
        "STAR --runThreadN {threads} "
        "--outSAMunmapped Within "
        "--genomeDir {params.ref} "
        "--readFilesIn {input} "
        "--outFileNamePrefix {params.prefix} "
        "--readFilesCommand zcat "
        "--quantMode TranscriptomeSAM --outFilterType BySJout "
        "{params.map_params} "
        "--outTmpDir $TMPDIR/STAR_ "
        "--outSAMType BAM Sorted "
        "--outSAMattrIHstart 0 "
        "--outSAMattributes NH HI NM MD AS nM "
        "--outSAMattrRGline ID:{params.FLI} SM:{params.Barcode}"
