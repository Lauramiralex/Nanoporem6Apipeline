rule sort:
    input:
        "{output_dir}/{name}{replicate}/aligned/basecalled.bam"
    output:
        "{output_dir}/{name}{replicate}/aligned/sorted.bam"
    conda:
        "../envs/environment_adv.yml"
    log:
        "{output_dir}/{name}{replicate}/log/sorted.log"
    shell:
        "samtools sort {input} > {output} 2> {log}"

rule summarize:
    input:
        "{output_dir}/{name}{replicate}/aligned/sorted.bam"
    output:
        "{output_dir}/{name}{replicate}/summary/summary_reads.txt"
    conda:
        "../envs/environment_adv.yml"
    shell:
        "modkit summary {input} > {output}"


rule index:
    input:
        "{output_dir}/{name}{replicate}/aligned/sorted.bam"
    output:
        "{output_dir}/{name}{replicate}/aligned/indexed.bam.bai"
    conda:
        "../envs/environment_adv.yml"
    log:
        "{output_dir}/{name}{replicate}/log/indexed.log"
    shell:
        "samtools index {input} > {output} 2> {log}"
