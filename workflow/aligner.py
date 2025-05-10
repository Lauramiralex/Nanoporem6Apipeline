import os
import sys
import pandas as pd
from pathlib import Path
import datetime

configfile: "config/config.yaml"

samples = pd.read_csv(config["samples"], sep="\t", header=0)

samples.sort_values(
    by=["name", "replicate", "filepath"],
    key=lambda col: col.astype(str).str.lower() if col.dtype == "object" else col,
    inplace=True
)

Path(config["output_dir"]).mkdir(parents=True, exist_ok=True)

genome = config["genome"]
timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

include: "rules/common.py"
include: "rules/aligner.py"

rule all:
    input:
        expand("{output_dir}/{name}{replicate}/aligned/basecalled.bam",zip, name=samples["name"],replicate= samples["replicate"],output_dir=[config["output_dir"]] * len(samples)),
        expand("{output_dir}/{name}{replicate}/summary/aligner_summary.tsv",zip, name=samples["name"],replicate= samples["replicate"],output_dir=[config["output_dir"]] * len(samples))
    output:
        f"{timestamp}_final_marker_align.done"  
    shell:
        """
        touch {output}
        """  

