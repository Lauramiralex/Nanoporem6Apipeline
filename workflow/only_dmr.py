import os
import sys
import pandas as pd
from pathlib import Path
import datetime

configfile: "config/config.yaml"

include: "rules/common.py"
include: "rules/dmr.py"

samples = pd.read_csv(config["samples"], sep="\t", header=0)
exp_names = get_exp_names("config/samples.tsv")
contr_names = get_contr_names("config/samples.tsv")

samples.sort_values(
    by=["name", "replicate", "filepath"],
    key=lambda col: col.astype(str).str.lower() if col.dtype == "object" else col,
    inplace=True
)

Path(config["output_dir"]).mkdir(parents=True, exist_ok=True)

threads = config["threads"]
genome = config["genome"]
base = config["modbase"]
timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

rule all:
    input:
        expand("{output_dir}/results/{dest}{contr}_dmr.txt",zip, dest = exp_names, contr = contr_names, output_dir=[config["output_dir"]] * len(exp_names))
    output:
        f"{timestamp}_final_marker_align.done"  
    shell:
        """
        touch {output}
        """  

