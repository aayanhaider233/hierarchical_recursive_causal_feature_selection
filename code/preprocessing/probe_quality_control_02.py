import pandas as pd


def remove_non_cpg_probes(methylation_matrix_path, batch):

    methylation_df = pd.read_csv(methylation_matrix_path, index_col=0)
    transposed_methylation_df = methylation_df.T

    cpg_only_methylation_df = transposed_methylation_df.loc[:, methylation_df.columns.str.startswith("cg")]

    cpg_only_methylation_df.insert(0, "sample_id", batch + "_" + methylation_df.index.astype(str))
    final_methylation_df = cpg_only_methylation_df.reset_index(drop=True)

    return final_methylation_df


def filter_by_metadata(methylation_df, metadata_df):

    valid_samples = set(methylation_df["sample_id"])

    filtered_metadata_df = metadata_df[metadata_df["sample_id"].isin(valid_samples)]
    filtered_methylation_df = methylation_df[methylation_df["sample_id"].isin(valid_samples)]

    return filtered_methylation_df, filtered_metadata_df


def filter_by_detection_pvalue(methylation_df, pval_threshold=0.05, label="_Detection_Pval"):

    pval_cols = [c for c in methylation_df.columns if c.endswith(label)]

    if len(pval_cols) == 0:
        return methylation_df

    probe_names = [c.replace(label, "") for c in pval_cols]

    good_mask = (methylation_df[pval_cols] <= pval_threshold).all(axis=0).values
    good_probes = [probe for probe, keep in zip(probe_names, good_mask) if keep]

    probe_filtered_methylation_df = methylation_df[["sample_id"] + good_probes]

    return probe_filtered_methylation_df


def clean_cpg_probes(methylation_df):

    sample_ids = methylation_df[["sample_id"]]
    values = methylation_df.drop(columns=["sample_id"])

    values = values.clip(0, 1)
    values = values.dropna(axis=1, thresh=int(0.9 * len(values)))

    return pd.concat([sample_ids, values], axis=1)