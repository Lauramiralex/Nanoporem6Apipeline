
rule align: 
    input:
        genome = config["genome"],
        bam ="{output_dir}/{name}{replicate}/basecalled.bam"
    output:
        "{output_dir}/{name}{replicate}/aligned/basecalled.bam"
    conda:
        "../envs/environment_base.yml"
    log:
        "{output_dir}/{name}{replicate}/log/aligned.log"
    shell:
        "dorado aligner {input.genome} {input.bam} --output-dir {wildcards.output_dir}/{wildcards.name}{wildcards.replicate}/aligned  --emit-summary 2> {log}"
        

rule summarize_aligner:
    input:
        "{output_dir}/{name}{replicate}/aligned/basecalled.bam"
    output:
        "{output_dir}/{name}{replicate}/summary/aligner_summary.tsv"
    conda:
        "../envs/environment_base.yml"
    log:
        "{output_dir}/{name}{replicate}/log/aligner_summary.log"
    shell:
        "dorado summary {input} > {output}"
