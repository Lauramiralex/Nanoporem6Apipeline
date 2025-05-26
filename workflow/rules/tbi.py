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

