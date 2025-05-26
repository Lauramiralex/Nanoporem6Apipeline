rule extract:
    input:
        inp= "{output_dir}/{name}{replicate}/aligned/sorted.bam",
        ref= config["genome"]
    output:
        "{output_dir}/{name}{replicate}/pileup/extracted.tsv"
    conda:
        "../envs/environment_adv.yml"
    shell:
        "modkit extract full {input.inp} {output}"

