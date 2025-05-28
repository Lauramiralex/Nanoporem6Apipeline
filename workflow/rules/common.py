import pandas

def get_input_files(wildcards):
        row = samples[(samples["name"] == wildcards.name) & (samples["replicate"]== int(wildcards.replicate))]
        return row["filepath"].iloc[0]

def get_path_pairs_exp(tsv_path):
    treated_pairs = []
    samples = pandas.read_csv(tsv_path, sep="\t")
    treated_df = samples[samples["group"] == "treated"]

    for _, row in treated_df.iterrows():
        path = "{output_dir}/" + row["name"] + str(row["replicate"]) + "/pileup/pileup.bed.gz"
        treated_pairs.append(path)
    return treated_pairs


def get_path_pairs_contr(tsv_path):
    untreated_pairs = []
    samples = pandas.read_csv(tsv_path, sep="\t")
    untreated_df = samples[samples["group"] == "control"]

    for _, row in untreated_df.iterrows():
        path = "{output_dir}/" + row["name"] + str(row["replicate"]) + "/pileup/pileup.bed.gz"
        untreated_pairs.append(path)
    return untreated_pairs

