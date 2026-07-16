import pandas as pd
from neuroCombat import neuroCombat

def batch_effect_correction(methylation_df, metadata_df, batch_info, chunk_size=10000):
    metadata_df = metadata_df.set_index("sample_id")
    batch_info = batch_info.set_index("sample")

    sample_ids = methylation_df["sample_id"]

    covariates = pd.DataFrame({
        "batch": batch_info.loc[sample_ids, "batch"].values,
        "disease": metadata_df.loc[sample_ids, "disease"].values
    }, index=sample_ids)

    values = methylation_df.drop(columns=["sample_id"]).T

    n_cpgs = values.shape[0]
    corrected_chunks = []

    for start in range(0, n_cpgs, chunk_size):
        end = min(start + chunk_size, n_cpgs)
        print(f"Processing CpGs {start}-{end}")

        chunk = values.iloc[start:end, :]

        combat_out = neuroCombat(
            dat=chunk,
            covars=covariates,
            batch_col="batch"
        )["data"]

        corrected_chunk = pd.DataFrame(
            combat_out,
            index=chunk.index,
            columns=chunk.columns
        )

        corrected_chunks.append(corrected_chunk)

    corrected_values = pd.concat(corrected_chunks, axis=0)

    corrected_values = corrected_values.T

    corrected_methylation_df = pd.concat([
        sample_ids.reset_index(drop=True),
        corrected_values.reset_index(drop=True)
    ], axis=1)

    return corrected_methylation_df