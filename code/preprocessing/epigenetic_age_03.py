import pyaging as pya

def generate_epigenetic_age(methylation_df):
    methylation_df = methylation_df.T
    if methylation_df.columns.duplicated().any():
        methylation_df = methylation_df.groupby(methylation_df.columns, axis=1).mean()

    adata = pya.pp.df_to_adata(methylation_df, imputer_strategy = 'knn')
    pya.pred.predict_age(adata, ['Horvath2013', 'Hannum', 'GrimAge2'])

    epi_df = adata.obs.reset_index().rename(columns={"index": "sample_id"})
    epi_df = epi_df[["sample_id", "horvath2013", "hannum", "grimage2"]]

    return epi_df

def add_epiage_to_metadata(metadata_df, epi_df):
    merged_df = metadata_df.merge(epi_df, on="sample_id", how="inner")

    column_order = [
        "sample_id",
        "age",
        "sex",
        "horvath2013",
        "hannum",
        "grimage2",
        "disease"
    ]

    merged_df = merged_df[column_order]

    return merged_df