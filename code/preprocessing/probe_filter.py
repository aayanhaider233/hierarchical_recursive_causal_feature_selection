import pandas as pd

def filter_probes(methylation_df, manifest_path, cross_reactive_probes_path, multimap_probes_path):
    
    cross = pd.read_csv(cross_reactive_probes_path, usecols=['TargetID'])['TargetID'].tolist()
    
    multi = pd.read_csv(multimap_probes_path, header=None)[0].tolist()
    
    manifest_df = pd.read_csv(manifest_path, skiprows=7)
    manifest_df.columns = manifest_df.columns.str.strip()
    manifest_df['CHR'] = manifest_df['CHR'].astype(str)
    sex_probes = manifest_df.loc[manifest_df['CHR'].isin(['X','Y']), 'IlmnID'].tolist()
    
    to_remove = set(cross) | set(multi) | set(sex_probes)

    filtered_methylation_df = methylation_df.loc[~methylation_df.index.isin(to_remove)]

    return filtered_methylation_df