import os
import sys
import pandas as pd
from pathlib import Path
import datetime

configfile: "config/config.yaml"

include: "rules/common.py"
include: "rules/dmr.py"

samples = pd.read_csv(config["samples"], sep="\t", header=0)

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
        f"{config["output_dir"]}/results/dmr_out.txt"
    output:
        f"{timestamp}_final_marker_align.done"  
    shell:
        """
        touch {output}
        """  

