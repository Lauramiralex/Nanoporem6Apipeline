def get_input_files(wildcards):
        row = samples[(samples["name"] == wildcards.name) & (samples["replicate"]== int(wildcards.replicate))]
        return row["filepath"].iloc[0]

def get_path_pairs_exp(samples):
    pairs = []
    first_paths= []
    grouped = samples.groupby("replicate")

    for rep, group in grouped:
        name = group["name"].tolist()
        paths = group["filepath"].tolist()
        
        if len(name) != 2:
            raise ValueError(f"Replicate {rep} has {len(name)} conditions (expected 2).")

        first_paths.append("{output_dir}/"+ name[0] +"{replicate}/pileup/pileup.bed.gz")

    return first_paths

def get_path_pairs_contr(samples):
    pairs = []
    second_paths= []
    grouped = samples.groupby("replicate")

    for rep, group in grouped:
        name = group["name"].tolist()
        paths = group["filepath"].tolist()
        
        if len(name) != 2:
            raise ValueError(f"Replicate {rep} has {len(name)} conditions (expected 2).")

        second_paths.append("{output_dir}/"+ name[0] +"{replicate}/pileup/pileup.bed.gz")

    return second_paths
