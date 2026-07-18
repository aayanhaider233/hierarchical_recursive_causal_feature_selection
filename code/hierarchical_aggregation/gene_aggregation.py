import pandas as pd

def load_dmr_matrix(dmr_matrix_path):
    return pd.read_csv(dmr_matrix_path)

def build_gene_dmr_lookup(dmr_gene_df):
    gene_dmr_lookup = (
        dmr_gene_df.groupby("Gene")["DMR"]
        .apply(list)
        .to_dict()
    )

    return gene_dmr_lookup

def aggregate_gene_methylation(dmr_matrix_df, gene_dmr_lookup):
    dmr_matrix_df = dmr_matrix_df.copy()

    sample_col = dmr_matrix_df.columns[0]

    gene_df = pd.DataFrame()
    gene_df[sample_col] = dmr_matrix_df[sample_col]

    available_dmrs = set(dmr_matrix_df.columns)

    for gene, dmr_list in gene_dmr_lookup.items():
        valid_dmrs = [dmr for dmr in dmr_list if dmr in available_dmrs]

        if len(valid_dmrs) == 0:
            continue

        gene_df[gene] = dmr_matrix_df[valid_dmrs].mean(axis=1)

    return gene_df