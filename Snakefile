
# Including rules
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
        config['project']+".genes.tsv"