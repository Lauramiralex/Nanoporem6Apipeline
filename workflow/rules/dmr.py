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

rule zipup:
    input:
        "{output_dir}/{name}{replicate}/pileup/pileup.bed"
    output:
        "{output_dir}/{name}{replicate}/pileup/pileup.bed.gz"
    conda:
        "../envs/environment_adv.yml"
    log:
        "{output_dir}/{name}{replicate}/log/zip.log"
    shell:
        "bgzip -f -k {input}"

rule index_zipped:
    input:
        "{output_dir}/{name}{replicate}/pileup/pileup.bed.gz"
    output:
        "{output_dir}/{name}{replicate}/pileup/pileup.bed.gz.tbi"
    conda:
        "../envs/environment_adv.yml"
    log:
        "{output_dir}/{name}{replicate}/log/index.log"
    shell:
        "tabix -f -p bed {input}"



rule dmr:
    input:
        exp =  lambda wildcards:get_path_pairs_exp(samples),
        contr  = lambda wildcards: get_path_pairs_contr(samples),
        ref= config["genome"],
        reg= config["genome_fai"]
    output:
        "{output_dir}/{name}{replicate}/dmr/{name}{replicate}_dmr.txt"
    conda:
        "../envs/environment_adv.yml"
    log:
        "{output_dir}/{name}{replicate}/log/{name}{replicate}_dmr.log"
    params:
        control_args= lambda wildcards, input: " ".join(f"-a {c}" for c in input.contr),
        experiment_args= lambda wildcards, input: " ".join(f"-b {e}" for e in input.exp)
    shell:
        "modkit dmr pair {params.experiment_args} {params.control_args} -o {output} --ref {input.ref} --base A --threads {threads}"
