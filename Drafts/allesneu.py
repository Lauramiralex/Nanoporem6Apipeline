import pandas as pd
configfile: "config/config.yaml"

experiments =pd.read_csv(config["experiments"], dtype=str, sep="\t")
controls =pd.read_csv(config["controlgroup"], dtype=str, sep="\t")
allfiles = pd.read_csv(config["allfiles"], dtype= str, sep="\t")
basecall_model = config["basecall_model"]
threads = config["threads"]
region = config["extract_region"]

experiments.sort_values(by=["experiment"], key=lambda col: col.str.lower() if col.dtype==str else col, inplace=True)
controls.sort_values(by=["control"], key=lambda col: col.str.lower() if col.dtype==str else col, inplace=True)
allfiles.sort_values(by=["allfiles"], key=lambda col: col.str.lower() if col.dtype==str else col, inplace=True)



rule all:
    input:
        expand("basecalled/{allfiles}.bam",zip, allfiles= allfiles["allfiles"]),
        expand("aligned/{allfiles}/pileup.bed",zip, allfiles= allfiles["allfiles"]),
        expand("aligned/{allfiles}/indexed.bam.bai",zip, allfiles= allfiles["allfiles"]),
        expand("aligned/{allfiles}/sorted.bam",zip,allfiles= allfiles["allfiles"]),
        expand("aligned/{allfiles}/pileup.bed.gz",zip,allfiles= allfiles["allfiles"]),
        expand("aligned/{allfiles}/pileup.bed.gz.tbi",zip,allfiles= allfiles["allfiles"]),
        expand("align2_out/{allfiles}.bam", zip,allfiles= allfiles["allfiles"]),
        expand("dmr/{experiment}-{control}.txt", zip, experiment = experiments["experiment"], control = controls["control"]),
        expand("aligned/{control}/summary_reads.txt", control= controls["control"]),
        expand("aligned/{experiment}/extracted.tsv" ,zip, experiment = experiments["experiment"])


        
rule call_to_base: 
    input:
        "pod5/{allfiles}"
    output:
        "basecalled/{allfiles}.bam"
    log:
        "log/{allfiles}_basecalled.log"
    resources:
        reparation_instances = 1
    shell:
        "dorado basecaller {basecall_model} {input} --modified-bases m6A_DRACH --batchsize 4000 > {output}  2> {log}"

rule align: 
    input:
        "mm39_fasta/Mus_musculus.GRCm39.cdna.all.fa",
        "basecalled/{allfiles}.bam"
    output:
        "align2_out/{allfiles}.bam"
    log:
        "log/{allfiles}_aligned.log"
    shell:
        "dorado aligner {input} --output-dir align2_out --emit-summary 2> {log}"

rule sort:
    input:
        "align2_out/{allfiles}.bam"
    output:
        "aligned/{allfiles}/sorted.bam"
    log:
        "log/{allfiles}_sorted.log"
    shell:
        "samtools sort {input} > {output} 2> {log}"        

rule index:
    input:
        "aligned/{allfiles}/sorted.bam"
    output:
        "aligned/{allfiles}/indexed.bam.bai"
    log:
        "log/{allfiles}_indexed.log"
    shell:
        "samtools index {input} > {output} 2> {log}"


rule base_count_sum:
    input:
        bam="aligned/{allfiles}/sorted.bam",
        bam_index="aligned/{allfiles}/indexed.bam.bai"
    output:
        "aligned/{allfiles}/pileup.bed"
    log:
        "log/{allfiles}_pileup.log"
    shell:
        "modkit pileup {input.bam} {output} --log-filepath pileup.log 2> {log}"


rule summarize:
    input:
        "aligned/{allfiles}/sorted.bam"
    output:
        "aligned/{allfiles}/summary_reads.txt"
    shell:
        "modkit summary {input} > {output}"

rule extract:
    input:
        inp= "aligned/{experiment}/sorted.bam",
        ref= "mm39_fasta/Mus_musculus.GRCm39.cdna.all.fa"
    output:
        "aligned/{experiment}/extracted.tsv"
    shell:
        "modkit extract full {input.inp} {output}"
        #"modkit extract full {input.inp} {output} --region {region} --ref {input.ref}"
        
rule zipup:
    input:
        "aligned/{allfiles}/pileup.bed"
    output:
        "aligned/{allfiles}/pileup.bed.gz"
    log:
        "log/{allfiles}_zip.log"
    shell:
        "bgzip -f -k {input}"


rule index_zipped:
    input:
        "aligned/{allfiles}/pileup.bed.gz"
    output:
        "aligned/{allfiles}/pileup.bed.gz.tbi"
    log:
        "log/{allfiles}_index.log"
    shell:
        "tabix -f -p bed {input}"


rule dmr:
    input:
        unused= "aligned/{experiment}/pileup.bed.gz.tbi",
        unusedd= "aligned/{control}/pileup.bed.gz.tbi",
        unuseddd= "aligned/{experiment}/summary_reads.txt",
        unusedddd= "aligned/{control}/summary_reads.txt",
        exp= lambda wildcards: expand("aligned/{experiment}/pileup.bed.gz", experiment= experiments["experiment"]),
        contr= lambda wildcards: expand("aligned/{control}/pileup.bed.gz", control = controls["control"]),
        ref= "mm39_fasta/Mus_musculus.GRCm39.cdna.all.fa",
        reg= "mm39_fasta/Mus_musculus.GRCm39.cdna.all.fa.fai"
    output:
        "dmr/{experiment}-{control}.txt"
    log:
        "log/{experiment}-{control}_dmr.log"
    params:
        control_args= lambda wildcards, input: " ".join(f"-a {c}" for c in input.contr),
        experiment_args= lambda wildcards, input: " ".join(f"-b {e}" for e in input.exp)
    shell:
        "modkit dmr pair {params.experiment_args} {params.control_args} -o {output} --ref {input.ref} --base A --threads {threads}"







