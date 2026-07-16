import pandas as pd

def concatenate_metadata(*metadata_dfs):
    return pd.concat(metadata_dfs, axis=0, ignore_index=True)

def generate_batch_info(metadata_df):
    batch_info = pd.DataFrame({
        "sample": metadata_df["sample_id"],
        "batch": metadata_df["sample_id"].str.split("_", n=1).str[0]
    })
    return batch_info

def concatenate_methylation_matrices(*methylation_dfs):
    return pd.concat(methylation_dfs, axis=0, ignore_index=True)