# --------------------------------------------------
# Imports
# --------------------------------------------------

import subprocess
from pathlib import Path
import dmr_to_gene_mapping as dgm
import gene_aggregation as ga


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
CAUSAL_DISCOVERY_INPUTS = DATA / "causal_discovery"
INTERMEDIATE_PATH.mkdir(parents=True, exist_ok=True)
DMR_RDS = INTERMEDIATE_PATH / "dmr_groups.rds"
DMR_CPG_MAP = INTERMEDIATE_PATH / "dmr_cpg_map.csv"
TRAIN_DMR = METHYLATION_PATH / "train_dmr_matrix.csv"
TEST_DMR  = METHYLATION_PATH / "test_dmr_matrix.csv"
DMR_GENE_MAP = INTERMEDIATE_PATH / "dmr_gene_map.csv"
GENE_MATRIX = INTERMEDIATE_PATH / "gene_methylation_matrix.csv"
GENE_MODULE_ASSIGNMENT = INTERMEDIATE_PATH / "gene_module_assignments.csv" 
ME_MATRIX = CAUSAL_DISCOVERY_INPUTS / "me_matrix.csv"
GENE_MODULE_MEMBERSHIP = INTERMEDIATE_PATH / "gene_module_kme_matrix.csv"

# --------------------------------------------------
# R code files
# --------------------------------------------------

AGGREGATION_PATH = ROOT / "code" / "hierarchical_aggregation"
DMR_DETECTION = AGGREGATION_PATH / "dmr_detection.R"
DMR_AGGREGATION = AGGREGATION_PATH / "dmr_aggregation.R"
WGCNA = AGGREGATION_PATH / "wgcna.R"


def run_aggregation():

    # --------------------------------------------------
    # 1. DMR detection & aggregation (CpG to DMR)
    # --------------------------------------------------

    subprocess.run([
        "Rscript",
        str(DMR_DETECTION),
        str(TRAIN_CPG),
        str(TRAIN_META),
        str(DMR_RDS),
        str(DMR_CPG_MAP)
    ], check=True)

    subprocess.run([
        "Rscript",
        str(DMR_AGGREGATION),
        str(TRAIN_CPG),
        str(DMR_RDS),
        str(TRAIN_DMR)
    ], check=True)

    subprocess.run([
        "Rscript",
        str(DMR_AGGREGATION),
        str(TEST_CPG),
        str(DMR_RDS),
        str(TEST_DMR)
    ], check=True)


    # --------------------------------------------------
    # 2. DMR to gene mapping
    # --------------------------------------------------

    manifest_df = dgm.load_manifest(
        manifest_path=MANIFEST_PATH
    )

    cpg_gene_lookup = dgm.build_cpg_gene_lookup(
        manifest_df=manifest_df
    )

    cpg_coords = dgm.build_cpg_coordinates(
        manifest_df=manifest_df
    )

    tss_df = dgm.build_tss_df(
        manifest_df=manifest_df
    )

    dmr_df = dgm.load_map(
        map_path=DMR_CPG_MAP
    )

    dmr_df = dgm.assign_constituent_cpg_genes(
        map_df=dmr_df, 
        lookup=cpg_gene_lookup
    )

    dmr_df = dgm.add_dmr_coordinates(
        dmr_df=dmr_df, 
        cpg_coords=cpg_coords
    )

    dmr_unmapped, promoter_df = dgm.assign_promoter_overlap_genes(
        dmr_df=dmr_df, 
        tss_df=tss_df
    )

    nearest_df = dgm.assign_nearest_tss_genes(
        promoter_df=promoter_df, 
        dmr_unmapped=dmr_unmapped, 
        tss_df=tss_df
    )

    dmr_gene_df = dgm.combine_gene_assignments(
        dmr_df=dmr_df, 
        promoter_df=promoter_df, 
        nearest_df=nearest_df
    )

    dmr_matrix = ga.load_dmr_matrix(
        dmr_matrix_path=TRAIN_DMR
    )

    gene_dmr_lookup = ga.build_gene_dmr_lookup(
        dmr_gene_df=dmr_gene_df
    )

    gene_matrix = ga.aggregate_gene_methylation(
        dmr_matrix_df=dmr_matrix,
        gene_dmr_lookup=gene_dmr_lookup
    )

    dmr_gene_df.to_csv(
        DMR_GENE_MAP, 
        index=False
    )

    gene_matrix.to_csv(
        GENE_MATRIX,
        index=False
    )


    # --------------------------------------------------
    # 3. WGCNA (Genes to modules represented by MEs)
    # --------------------------------------------------

    subprocess.run([
        "Rscript",
        str(WGCNA),
        str(GENE_MATRIX),
        str(GENE_MODULE_ASSIGNMENT),
        str(ME_MATRIX),
        str(GENE_MODULE_MEMBERSHIP)
    ], check=True)


if __name__ == "__main__":
    run_aggregation()