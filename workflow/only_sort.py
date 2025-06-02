import os
import sys
import pandas as pd
from pathlib import Path
import datetime

configfile: "config/config.yaml"

include: "rules/common.py"
include: "rules/sort.py"

samples = pd.read_csv(config["samples"], sep="\t", header=0)

samples.sort_values(
    by=["name", "replicate", "filepath"],
    key=lambda col: col.astype(str).str.lower() if col.dtype == "object" else col,
    inplace=True
)

Path(config["output_dir"]).mkdir(parents=True, exist_ok=True)

timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")



rule all:
    input:
        expand("{output_dir}/{name}{replicate}/aligned/sorted.bam",zip, name=samples["name"],replicate= samples["replicate"],output_dir=[config["output_dir"]] * len(samples)),
        expand("{output_dir}/{name}{replicate}/summary/summary_reads.txt",zip, name=samples["name"],replicate= samples["replicate"],output_dir=[config["output_dir"]] * len(samples)),
        expand("{output_dir}/{name}{replicate}/aligned/indexed.bam.bai",zip, name=samples["name"],replicate= samples["replicate"],output_dir=[config["output_dir"]] * len(samples))  
    output:
        f"{timestamp}_final_marker_align.done"  
    shell:
        """
        touch {output}
        """  

