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
timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")




experiments.sort_values(by=["experiment"], key=lambda col: col.str.lower() if col.dtype==str else col, inplace=True)
controls.sort_values(by=["control"], key=lambda col: col.str.lower() if col.dtype==str else col, inplace=True)
allfiles.sort_values(by=["allfiles"], key=lambda col: col.str.lower() if col.dtype==str else col, inplace=True)
path.sort_values(by=["path"], key=lambda col: col.str.lower() if col.dtype==str else col, inplace=True)



rule all:
    input:
        expand("{path}/{allfiles}/results/basecalled/{allfiles}.bam",zip, path = path["path"],allfiles= allfiles["allfiles"]),
        expand("{path}/{allfiles}/results/basecalled/basecaller_summary.tsv", zip,path = path["path"], allfiles= allfiles["allfiles"])
    output:
        f"{timestamp}_final_marker_basecall.done"
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


        
