def get_input_files(wildcards):
        row = samples[(samples["name"] == wildcards.name) & (samples["replicate"]== int(wildcards.replicate))]
        return row["filepath"].iloc[0]


def get_path_pairs_exp(samples):
    treated_pairs = []
    treated_df = samples[samples["group"] == "treated"]
    for _, row in treated_df.iterrows():
        replicate = row["replicate"]
        path = "{output_dir}/" + row["name"] + str(row["replicate"]) + "/pileup/pileup.bed.gz"
        treated_pairs.append(path)

    print(treated_pairs)
    return treated_pairs

def get_path_pairs_contr(samples):
    untreated_pairs = []

    untreated_df = samples[samples["group"] == "control"]
    for _, row in untreated_df.iterrows():
        replicate = row["replicate"]
        path = "{output_dir}/" + row["name"] + str(row["replicate"]) + "/pileup/pileup.bed.gz"
        untreated_pairs.append(path)
    print(untreated_pairs)
    return untreated_pairs


