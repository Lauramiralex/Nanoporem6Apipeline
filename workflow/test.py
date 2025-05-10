import os
import sys
import pandas as pd
from pathlib import Path
import datetime

configfile: "config/config.yaml"
#os.chdir(config["workdirection"])
samples = pd.read_csv("samples.tsv", sep="\t")
filename = "testlauf"


def get_path_pairs(samples):
    pairs = []
    grouped = samples.groupby("replicate")

    for rep, group in grouped:
        name = group["name"].tolist()
        paths = group["filepath"].tolist()
        
        if len(name) != 2:
            raise ValueError(f"Replicate {rep} has {len(name)} conditions (expected 2).")
        
        # Erstelle das Paar der Paths für die beiden Conditions
        path_pair = (paths[0], paths[1])
        pairs.append(path_pair)

    return pairs

def save_path_pairs_to_txt(samples, filename):
    pairs = get_path_pairs(samples)
    
    with open(filename, 'w') as f:
        for pair in pairs:
            f.write(f"{pair[0]} {pair[1]}\n")

#get_path_pairs(samples)
save_path_pairs_to_txt(samples, filename)
