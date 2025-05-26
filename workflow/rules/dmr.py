rule dmr:
    input:
        exp =  lambda wildcards:get_path_pairs_exp(samples),
        contr  = lambda wildcards:get_path_pairs_contr(samples),
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
