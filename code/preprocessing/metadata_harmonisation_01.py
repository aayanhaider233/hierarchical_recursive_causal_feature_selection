import pandas as pd
import re
import gzip
from difflib import get_close_matches

COLUMN_CANONICAL = {
    'disease': [
        'disease_status', 'diagnosis', 'condition', 'group', 'status',
        'disease state', 'disease type', 'diseasegroup', 'diseasecategory'
    ],
    'age': [
        'age_years', 'chronological age', 'age (years)', 'subject age',
        'age_year', 'age_at_sampling'
    ],
    'sex': [
        'gender', 'biological sex', 'sex of individual', 'sex/gender'
    ],
    'source_tissue': [
        'source tissue', 'tissue', 'source (tissue)'
    ],
}

def normalise_columns(df, similarity_threshold = 0.8):
    df = df.copy()
    cleaned_cols = (
        df.columns
        .str.strip()
        .str.lower()
        .str.replace(r'[\s\-]+', '_', regex=True)
        .str.replace(r'[^0-9a-zA-Z_]', '', regex=True)
    )

    normalised_cols = []
    
    for col in cleaned_cols:
        mapped_name = None
        for canonical, variants in COLUMN_CANONICAL.items():
            if col in variants or col == canonical:
                mapped_name = canonical
                break
        if not mapped_name:
            all_variants = {v: k for k, vals in COLUMN_CANONICAL.items() for v in vals}
            close = get_close_matches(col, all_variants.keys(), n=1, cutoff=similarity_threshold)
            if close:
                mapped_name = all_variants[close[0]]
        normalised_cols.append(mapped_name if mapped_name else col)

    df.columns = normalised_cols
    return df

def parse_metadata(metadata_path, batch):
    metadata = {}
    gsm_list = []

    open_func = gzip.open if str(metadata_path).endswith(".gz") else open

    with open_func(metadata_path, "rt", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("!Sample_geo_accession"):
                gsm_list = re.findall(r"GSM\d+", line)
                for gsm in gsm_list:
                    if gsm not in metadata:
                        metadata[gsm] = {}
                continue
            if line.startswith("!Sample_description"):
                values = re.findall(r'"([^"]*)"', line)
                for i, val in enumerate(values):
                    if i < len(gsm_list):
                        metadata[gsm_list[i]]['barcode'] = val.strip()
            if line.startswith("!Sample_characteristics_ch1"):
                values = re.findall(r'"([^"]+)"', line)
                key_values = [v.split(': ', 1) for v in values]
                for gsm, kv in zip(gsm_list, key_values):
                    if len(kv) == 2:
                        key, value = kv
                        metadata[gsm][key.strip()] = value.strip()

    meta_df = pd.DataFrame.from_dict(metadata, orient='index')
    meta_df = meta_df.reset_index(drop=True)
    
    meta_df['sample_id'] = batch + '_' + meta_df['barcode']
    
    meta_df = normalise_columns(meta_df)

    return meta_df 

def clean_metadata(meta_df, control_label = "0", case_label = "1", male_label = "M", female_label = "F", drop_columns = "source_tissue"):
    meta_df = meta_df[~meta_df['age'].isin(['NA', 'NaN', '', None])]
    meta_df = meta_df.dropna(subset = ['age'])
    meta_df['age'] = meta_df['age'].astype(float)
    
    meta_df['sex'] = meta_df['sex'].map({male_label: 0, female_label: 1})
    meta_df['disease'] = meta_df['disease'].map({control_label: 0, case_label: 1})
    
    missing_mask = meta_df.isna().any(axis=1)
    meta_df = meta_df[~missing_mask]
    
    if 'age' in meta_df.columns:
        invalid_age_mask = ~meta_df['age'].between(0, 120)
        meta_df = meta_df[~invalid_age_mask]

    if drop_columns is not None and drop_columns in meta_df.columns:
        meta_df = meta_df.drop(columns=drop_columns)
    return meta_df