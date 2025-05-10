#!/usr/bin/env python3

import os
import shutil
import pandas as pd
import yaml

#Load config.yaml
with open("config/config.yaml", "r") as f:
    config = yaml.safe_load(f)

os.chdir(config["workdirection"])

print("Loaded config:", config)

destination = config["resultdirectory"]
os.makedirs(destination, exist_ok=True)

# Load TSV (containing sample names)
# Update this path if needed
df = pd.read_csv("config/allfiles.tsv", sep="\t", dtype=str)
os.makedirs(destination, exist_ok=True)



#Loop through sample names and copy folders
#sample_names = [os.path.join(row['path'], row['allfiles']) for _, row in df.iterrows()]


#for sample in sample_names:
#    src_folder = os.path.join(sample, "results")
 #   dest_folder = os.path.join(destination, row['allfiles'])
for _, row in df.iterrows():
    sample = os.path.join(row['path'], row['allfiles'])
    src_folder = os.path.join(sample, "results")
    dest_folder = os.path.join(destination, row['allfiles'])
    
    if os.path.isdir(src_folder):
        shutil.copytree(src_folder, dest_folder, dirs_exist_ok=True)
        print(f"✅ Copied {src_folder} → {dest_folder}")
    else:
        print(f"⚠️  Source folder not found: {src_folder}")

