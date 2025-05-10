rule call_to_base: 
    input:
        get_input_files
    output:
        txt="{output_dir}/{name}{replicate}/basecalled.bam"
    log:
        "{output_dir}/{name}{replicate}/log/basecalled.log"
    resources:
        reparation_instances = 1
    conda:
        "../envs/environment_base.yml"
    shell:
        "dorado basecaller {basecall_model} {input} --modified-bases {base} --batchsize {batchsize} > {output.txt}  2> {log}"

rule summarize_basecaller:
    input:
        "{output_dir}/{name}{replicate}/basecalled.bam"
    output:
        "{output_dir}/{name}{replicate}/summary/basecaller_summary.tsv"
    conda:
        "../envs/environment_base.yml"
    log:
        "{output_dir}/{name}{replicate}/log/basecaller_summary.log"
    shell:
        "dorado summary {input} > {output}"
