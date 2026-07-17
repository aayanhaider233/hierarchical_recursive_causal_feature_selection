# --------------------------------------------------
# Imports
# --------------------------------------------------

import subprocess
from pathlib import Path
import dmr_gene_mapping as dgm


# --------------------------------------------------
# Dataset directories
# --------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

METHYLATION_PATH = DATA / "methylation"
TRAIN_CPG = METHYLATION_PATH / "train_cpg_mval_matrix_corrected.csv.csv"
TEST_CPG  = METHYLATION_PATH / "test_cpg_mval_matrix_corrected.csv.csv"

METADATA_PATH = DATA / "metadata"
TRAIN_META = METADATA_PATH / "train_metadata.csv"

MANIFEST_PATH = ROOT / "annotations" / "methylation" / "humanmethylation450_15017482_v1-2.csv"


# --------------------------------------------------
# Save directories
# --------------------------------------------------

INTERMEDIATE_PATH = DATA / "intermediate"
INTERMEDIATE_PATH.mkdir(parents=True, exist_ok=True)
DMR_RDS = INTERMEDIATE_PATH / "dmr_groups.rds"
DMR_CPG_MAP = INTERMEDIATE_PATH / "dmr_cpg_map.csv"
TRAIN_DMR = METHYLATION_PATH / "train_dmr_matrix.csv"
TEST_DMR  = METHYLATION_PATH / "test_dmr_matrix.csv"


# --------------------------------------------------
# Present working directory
# --------------------------------------------------

AGGREGATION_PATH = ROOT / "code" / "hierarchical_aggregation"


# --------------------------------------------------
# 1. DMR detection & aggregation (CpG to DMR)
# --------------------------------------------------

subprocess.run([
    "Rscript",
    str(AGGREGATION_PATH / "dmr_detection.R"),
    str(TRAIN_CPG),
    str(TRAIN_META),
    str(DMR_RDS),
    str(DMR_CPG_MAP)
], check=True)

subprocess.run([
    "Rscript",
    str(AGGREGATION_PATH / "dmr_aggregation.R"),
    str(TRAIN_CPG),
    str(DMR_RDS),
    str(TRAIN_DMR)
], check=True)

subprocess.run([
    "Rscript",
    str(AGGREGATION_PATH / "dmr_aggregation.R"),
    str(TEST_CPG),
    str(DMR_RDS),
    str(TEST_DMR)
], check=True)


# --------------------------------------------------
# 2. DMR to gene mapping
# --------------------------------------------------

manifest_df = dgm.load_manifest(MANIFEST_PATH)
lookup = dgm.build_cpg_gene_lookup(manifest_df)
cpg_coords = dgm.build_cpg_coordinates(manifest_df)
tss_df = dgm.build_tss_df(manifest_df)
dmr_df = dgm.load_map(DMR_CPG_MAP)
dmr_df = dgm.assign_constituent_cpg_genes(dmr_df, lookup)
dmr_df = dgm.add_dmr_coordinates(dmr_df, cpg_coords)
dmr_unmapped, promoter_df = dgm.assign_promoter_overlap_genes(dmr_df, tss_df)
nearest_df = dgm.assign_nearest_tss_genes(promoter_df, dmr_unmapped, tss_df)
final_df = dgm.combine_gene_assignments(dmr_df, promoter_df, nearest_df)
final_df.to_csv(INTERMEDIATE_PATH / "dmr_to_gene.csv", index=False)


# --------------------------------------------------
# 3. WGCNA (Genes to modules represented by MEs)
# --------------------------------------------------

