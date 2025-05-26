def get_input_files(wildcards):
        row = samples[(samples["name"] == wildcards.name) & (samples["replicate"]== int(wildcards.replicate))]
        return row["filepath"].iloc[0]

def get_path_pairs_exp(samples):
    treated_paths= []
    grouped = samples.groupby("replicate")

    for rep, group in grouped:
        treated_rows = group[group["group"] == "treated"]
        paths = treated_rows["filepath"].tolist()
        for _,row in treated_rows["name"]:
                name = row["name"]
                treated_paths.append("{output_dir}/"+ name +"{replicate}/pileup/pileup.bed.gz")
    print(treated_paths)
    return treated_paths

def get_path_pairs_contr(samples):
    control_paths= []
    grouped = samples.groupby("replicate")

    for rep, group in grouped:
        control_rows = group[group["group"] == "control"]
        paths = control_rows["filepath"].tolist()
        for _,row in control_paths["name"]:
                name = row["name"]
                control_paths.append("{output_dir}/"+ name +"{replicate}/pileup/pileup.bed.gz")
                
        #name = control_rows["name"].values[0]
        #control_paths.append("{output_dir}/"+ name +"{replicate}/pileup/pileup.bed.gz")
    print(control_paths)
    return control_paths

