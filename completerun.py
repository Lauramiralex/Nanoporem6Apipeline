import pandas as pd
import datetime
configfile: "config/config.yaml"
os.chdir(config["workdirection"])

experiments =pd.read_csv(config["experiments"], dtype=str, sep="\t")
controls =pd.read_csv(config["controlgroup"], dtype=str, sep="\t")
allfiles = pd.read_csv(config["allfiles"], dtype= str, sep="\t")
path = pd.read_csv(config["path"], dtype= str, sep="\t")
basecall_model = config["basecall_model"]
threads = config["threads"]
region = config["extract_region"]
base = config["modbase"]
batchsize = config["batchsize"]
genome = config["genome"]
genome_fai= config["genome_fai"]
timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")




experiments.sort_values(by=["experiment"], key=lambda col: col.str.lower() if col.dtype==str else col, inplace=True)
controls.sort_values(by=["control"], key=lambda col: col.str.lower() if col.dtype==str else col, inplace=True)
allfiles.sort_values(by=["allfiles"], key=lambda col: col.str.lower() if col.dtype==str else col, inplace=True)
path.sort_values(by=["path"], key=lambda col: col.str.lower() if col.dtype==str else col, inplace=True)



rule all:
    input:
        expand("{path}/{allfiles}/results/basecalled/{allfiles}.bam",zip, path = path["path"],allfiles= allfiles["allfiles"]),
        expand("{path}/{allfiles}/results/basecalled/basecaller_summary.tsv", zip,path = path["path"], allfiles= allfiles["allfiles"]),
        expand("{path}/{allfiles}/results/aligned/aligner_summary.tsv",zip,path = path["path"], allfiles= allfiles["allfiles"]),
        expand("{path}/{allfiles}/results/aligned/{allfiles}.bam", zip,path = path["path"],allfiles= allfiles["allfiles"]),
        expand("{path}/{allfiles}/results/aligned/sorted.bam",zip,path = path["path"], allfiles= allfiles["allfiles"]),
        expand("{path}/{allfiles}/results/aligned/summary_reads.txt", path = path["path"], allfiles= allfiles["allfiles"]),
        expand("{path}/{allfiles}/results/aligned/indexed.bam.bai",zip, path = path["path"], allfiles= allfiles["allfiles"]),
        expand("{path}/{allfiles}/results/pileup/pileup.bed",zip, path = path["path"], allfiles= allfiles["allfiles"]),
        expand("{path}/{experiment}/results/pileup/extracted.tsv" ,zip,path = path["path"], experiment = experiments["experiment"]),
        expand("{path}/{allfiles}/results/pileup/pileup.bed.gz",zip,path = path["path"],allfiles= allfiles["allfiles"]),
        expand("{path}/{allfiles}/results/pileup/pileup.bed.gz.tbi",zip,path = path["path"], allfiles= allfiles["allfiles"]),
        expand("{path}/{experiment}/results/dmr/{experiment}-{control}.txt", zip, path= path["path"], experiment = experiments["experiment"], control = controls["control"])
    output:
        f"{timestamp}_final_marker_completerun.done"
    shell:
        """
        python config/move.py
        touch {output}
        """ 
        
rule call_to_base: 
    input:
        "{path}/{allfiles}"
    output:
        "{path}/{allfiles}/results/basecalled/{allfiles}.bam"
    log:
        "{path}/{allfiles}/results/log/basecalled.log"
    resources:
        reparation_instances = 1
    shell:
        "dorado basecaller {basecall_model} {input} --modified-bases {base} --batchsize {batchsize} > {output}  2> {log}"
       

rule summarize_basecaller:
    input:
        "{path}/{allfiles}/results/basecalled/{allfiles}.bam"
    output:
        "{path}/{allfiles}/results/basecalled/basecaller_summary.tsv"
    log:
        "{path}/{allfiles}/results/log/basecaller_summary.log"
    shell:
        "dorado summary {input} > {output}"
        
        
rule align: 
    input:
        genome = config["genome"],
        bam ="{path}/{allfiles}/results/basecalled/{allfiles}.bam"
    output:
        "{path}/{allfiles}/results/aligned/{allfiles}.bam"
    log:
        "{path}/{allfiles}/results/log/aligned.log"
    shell:
        "dorado aligner {input.genome} {input.bam} --output-dir {wildcards.path}/{wildcards.allfiles}/results/aligned  --emit-summary 2> {log}"
        

rule summarize_aligner:
    input:
        "{path}/{allfiles}/results/aligned/{allfiles}.bam"
    output:
        "{path}/{allfiles}/results/aligned/aligner_summary.tsv"
    log:
        "{path}/{allfiles}/results/log/aligner_summary.log"
    shell:
        "dorado summary {input} > {output}"

rule sort:
    input:
        "{path}/{allfiles}/results/aligned/{allfiles}.bam"
    output:
        "{path}/{allfiles}/results/aligned/sorted.bam"
    log:
        "{path}/{allfiles}/results/log/sorted.log"
    shell:
        "samtools sort {input} > {output} 2> {log}"

rule summarize:
    input:
        "{path}/{allfiles}/results/aligned/sorted.bam"
    output:
        "{path}/{allfiles}/results/aligned/summary_reads.txt"
    shell:
        "modkit summary {input} > {output}"


rule index:
    input:
        "{path}/{allfiles}/results/aligned/sorted.bam"
    output:
        "{path}/{allfiles}/results/aligned/indexed.bam.bai"
    log:
        "{path}/{allfiles}/results/log/indexed.log"
    shell:
        "samtools index {input} > {output} 2> {log}"

rule base_count_sum:
    input:
        bam="{path}/{allfiles}/results/aligned/sorted.bam",
        bam_index="{path}/{allfiles}/results/aligned/indexed.bam.bai"
    output:
        "{path}/{allfiles}/results/pileup/pileup.bed"
    log:
        "{path}/{allfiles}/results/log/pileup.log"
    shell:
        "modkit pileup {input.bam} {output} --log-filepath pileup.log 2> {log}"

rule extract:
    input:
        inp= "{path}/{experiment}/results/aligned/sorted.bam",
        ref= config["genome"]
    output:
        "{path}/{experiment}/results/pileup/extracted.tsv"
    shell:
        "modkit extract full {input.inp} {output}"

rule zipup:
    input:
        "{path}/{allfiles}/results/pileup/pileup.bed"
    output:
        "{path}/{allfiles}/results/pileup/pileup.bed.gz"
    log:
        "{path}/{allfiles}/results/log/zip.log"
    shell:
        "bgzip -f -k {input}"

rule index_zipped:
    input:
        "{path}/{allfiles}/results/pileup/pileup.bed.gz"
    output:
        "{path}/{allfiles}/results/pileup/pileup.bed.gz.tbi"
    log:
        "{path}/{allfiles}/results/log/index.log"
    shell:
        "tabix -f -p bed {input}"



rule dmr:
    input:
        unused= lambda wildcards: expand("{path}/{experiment}/results/pileup/pileup.bed.gz.tbi", path = path["path"], experiment= experiments["experiment"]),
        unusedd= lambda wildcards: expand("{path}/{control}/results/pileup/pileup.bed.gz.tbi",path=path["path"], control = controls["control"]),
        unuseddd= "{path}/{experiment}/results/aligned/summary_reads.txt",
        unusedddd= "{path}/{control}/results/aligned/summary_reads.txt",
        exp= lambda wildcards: expand("{path}/{experiment}/results/pileup/pileup.bed.gz", path=path["path"], experiment= experiments["experiment"]),
        contr= lambda wildcards: expand("{path}/{control}/results/pileup/pileup.bed.gz",path=path["path"], control = controls["control"]),
        ref= config["genome"],
        reg= config["genome_fai"]
    output:
        "{path}/{experiment}/results/dmr/{experiment}-{control}.txt"
    log:
        "{path}/{experiment}/results/log/{experiment}-{control}_dmr.log"
    params:
        control_args= lambda wildcards, input: " ".join(f"-a {c}" for c in input.contr),
        experiment_args= lambda wildcards, input: " ".join(f"-b {e}" for e in input.exp)
    shell:
        "modkit dmr pair {params.experiment_args} {params.control_args} -o {output} --ref {input.ref} --base A --threads {threads}"

