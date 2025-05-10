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

rule base_count_sum:
    input:
        bam="{output_dir}/{name}{replicate}/aligned/sorted.bam",
        bam_index="{output_dir}/{name}{replicate}/aligned/indexed.bam.bai"
    output:
        "{output_dir}/{name}{replicate}/pileup/pileup.bed"
    conda:
        "../envs/environment_adv.yml"
    log:
        "{output_dir}/{name}{replicate}/log/pileup.log"
    shell:
        "modkit pileup {input.bam} {output} --log-filepath pileup.log 2> {log}"
