rule dmr:
    input:
        exp =  get_path_pairs_exp("config/samples.tsv"),
        contr  = get_path_pairs_contr("config/samples.tsv"),
        #dest = get_exp_names("config/samples.tsv"),
        ref= config["genome"],
        reg= config["genome_fai"]
    output:
        "{output_dir}/results/dmr_out.txt"
    conda:
        "../envs/environment_adv.yml"
    log:
        "{output_dir}/results/log/dmr_out.log"
    params:
        control_args= lambda wildcards, input: " ".join(f"-a {c}" for c in input.contr),
        experiment_args= lambda wildcards, input: " ".join(f"-b {e}" for e in input.exp)
    shell:
        "modkit dmr pair {params.experiment_args} {params.control_args} -o {output} --ref {input.ref} --base A --threads {threads}"
