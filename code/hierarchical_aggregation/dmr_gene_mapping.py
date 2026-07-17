import pandas as pd

PROMOTER_UPSTREAM = 2000  
PROMOTER_DOWNSTREAM = 500 

def load_manifest(manifest_path):
    manifest = pd.read_csv(manifest_path, skiprows=7, engine='python', encoding='utf-8-sig')
    manifest.columns = manifest.columns.str.strip()
    return manifest 

def build_cpg_gene_lookup(manifest_df):
    manifest_genic = manifest_df[
        manifest_df['UCSC_RefGene_Name'].notna() 
        & (manifest_df['UCSC_RefGene_Name'] != '')
    ]
    return {
        row['IlmnID'] : row['UCSC_RefGene_Name'].split(';') 
        for _, row in manifest_genic.iterrows()
    } 

def load_map(map_path):
    map_df = pd.read_csv(map_path)
    map_df['cpg_list'] = map_df['CpGs'].str.split(';')
    return map_df

def constituent_cpg_map(cpg_list, lookup):
    genes = set()
    for cpg in cpg_list:
        if cpg in lookup:
            genes.update(lookup[cpg])
    return list(genes)

def assign_constituent_cpg_genes(map_df, lookup):
    map_df = map_df.copy()
    map_df["Genes"] = map_df["cpg_list"].apply(
        lambda x: constituent_cpg_map(x, lookup)
    )
    # map_df["Assignment_Type"] = map_df["Genes"].apply(lambda x: "CpG" if len(x) > 0 else "")
    map_df["Assignment_Type"] = map_df["Genes"].str.len().gt(0).map({True: "CpG", False: ""})

    return map_df

def build_cpg_coordinates(manifest_df):
    return manifest_df[['IlmnID','CHR','MAPINFO']].set_index('IlmnID')

def get_dmr_coords(cpg_list, cpg_coords):
    coords = cpg_coords.loc[cpg_list]
    chr_ = coords['CHR'].unique()
    if len(chr_) != 1:
        raise ValueError(f"DMR contains CpGs on multiple chromosomes: {chr_}")
    return pd.Series({
        'Chr': chr_[0],
        'Start': coords['MAPINFO'].min(),
        'End': coords['MAPINFO'].max()
    })

def add_dmr_coordinates(dmr_df, cpg_coords):
    dmr_coords = dmr_df['cpg_list'].apply(
        lambda x: get_dmr_coords(
            x, 
            cpg_coords
        )
    )
    return pd.concat([dmr_df, dmr_coords], axis=1)

def build_tss_df(manifest_df):
    tss_df = manifest_df[["CHR", "MAPINFO", "UCSC_RefGene_Name", "Strand"]].drop_duplicates()
    tss_df = tss_df.rename(columns={'CHR':'Chr','MAPINFO':'TSS','UCSC_RefGene_Name':'Gene'})
    tss_df = tss_df[tss_df["Gene"].notna() & (tss_df["Gene"] != "")]
    tss_df['promoter_start'] = tss_df['TSS'] - PROMOTER_UPSTREAM
    tss_df['promoter_end'] = tss_df['TSS'] + PROMOTER_DOWNSTREAM
    return tss_df

def assign_promoter_overlap_genes(dmr_df, tss_df):
    dmr_unmapped = dmr_df[dmr_df['Assignment_Type'] == ""].copy()
    promoter_assignments = []
    tss_by_chr = {chr_: df for chr_, df in tss_df.groupby("Chr")}
    for _, dmr in dmr_unmapped.iterrows():
        chr_tss = tss_by_chr.get(dmr["Chr"])
        if chr_tss is None:
            continue
        
        hits = chr_tss[ (chr_tss["promoter_start"] <= dmr["End"]) & (chr_tss["promoter_end"] >= dmr["Start"]) ]

        for gene in hits["Gene"]:
            promoter_assignments.append({
                "DMR": dmr["DMR"],
                "Gene": gene,
                "Assignment_Type": "Promoter",
                "Distance_to_TSS": 0
            })

    return dmr_unmapped, pd.DataFrame(promoter_assignments)

def assign_nearest_tss_genes(promoter_df, dmr_unmapped, tss_df):
    mapped_dmrs = set(promoter_df['DMR']) if not promoter_df.empty else set()
    nearest_dmrs = dmr_unmapped[~dmr_unmapped['DMR'].isin(mapped_dmrs)]

    nearest_assignments = []

    tss_by_chr = {chr_: df for chr_, df in tss_df.groupby("Chr")}

    for _, dmr in nearest_dmrs.iterrows():
        chr_tss = tss_by_chr.get(dmr["Chr"])
        if chr_tss is None:
            continue
        
        distances = chr_tss["TSS"].apply( lambda x: min(abs(dmr["Start"] - x), abs(dmr["End"] - x)) )

        min_idx = distances.idxmin()

        nearest_assignments.append({
            "DMR": dmr["DMR"],
            "Gene": chr_tss.loc[min_idx, "Gene"],
            "Assignment_Type": "Nearest",
            "Distance_to_TSS": distances[min_idx]
        })

    return pd.DataFrame(nearest_assignments)

def combine_gene_assignments(dmr_df, promoter_df, nearest_df):
    dmr_cpg_df = dmr_df[dmr_df["Assignment_Type"] == "CpG"].explode("Genes")

    dmr_cpg_df = dmr_cpg_df[["DMR", "Genes"]].rename(columns={"Genes": "Gene"})

    dmr_cpg_df["Assignment_Type"] = "CpG"
    dmr_cpg_df["Distance_to_TSS"] = 0

    final_df = pd.concat([dmr_cpg_df, promoter_df, nearest_df], ignore_index=True)

    final_df["Gene"] = final_df["Gene"].str.split(";")
    final_df = final_df.explode("Gene")

    final_df = final_df.drop_duplicates(subset=["DMR", "Gene"]).reset_index(drop=True)

    return final_df