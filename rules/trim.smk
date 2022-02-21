
def set_length(wildcards):
    if config['mode'] == 'small':
        return 16
    else:
        return 18

def get_prefix_trim(wildcards):
    return wildcards.FLI

def get_fastq(wildcards):
    fq_a = [meta_seq.F1.loc[wildcards.FLI]]
    if config['paired']:
        fq_a.append(meta_seq.F2.loc[wildcards.FLI])
    return fq_a

def set_params(wildcards):
    params = ""
    if config['paired']:
        params += "--paired"
    return params

        
if config['paired']:
    rule trim_pe:
        input:
            get_fastq
        params:
            length = set_length,
            basename = get_prefix_trim,
            params = set_params,
            outdir = config['path']["trim"]

        output: 
            config['path']['trim']+"/{FLI}_val_1.fq.gz",
            config['path']['trim']+"/{FLI}_val_2.fq.gz"

        threads: 8
        
        shell:
            "trim-galore "
            "--fastqc "
            "--length {params.length} "
            "{params.params} "
            "-j {threads} "
            "-o {params.outdir} "
            "--basename {params.basename} "
            "{input}"
else:
    rule trim_se:
        input:
            get_fastq
        params:
            length = set_length,
            basename = get_prefix_trim,
            params = set_params,
            outdir = config['path']["trim"]

        output: 
            config['path']['trim']+"/{FLI}_trimmed.fq.gz"
        
        threads: 8

        shell:
            "trim-galore "
            "--fastqc "
            "--length {params.length} "
            "{params.params} "
            "-j {threads} "
            "-o {params.outdir} "
            "--basename {params.basename} "
            "{input}"