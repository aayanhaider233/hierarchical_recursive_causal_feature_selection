from sklearn.model_selection import train_test_split

def set_partition(methylation_df, metadata_df, batch_info, test_size=0.3, stratify_class="disease"):
    common_samples = sorted(set(methylation_df["sample_id"]) & set(metadata_df["sample_id"]))

    methylation_df = methylation_df.set_index("sample_id").loc[common_samples]
    metadata_df = metadata_df.set_index("sample_id").loc[common_samples]

    labels = metadata_df[stratify_class]

    train_samples, test_samples = train_test_split(
        common_samples,
        test_size=test_size,
        random_state=0,
        stratify=labels
    )

    train_methylation = methylation_df.loc[train_samples]
    test_methylation  = methylation_df.loc[test_samples]
    train_metadata = metadata_df.loc[train_samples]
    test_metadata  = metadata_df.loc[test_samples]
    train_batch = batch_info[train_samples]
    test_batch = batch_info[test_samples]

    return train_methylation, test_methylation, train_metadata, test_metadata, train_batch, test_batch