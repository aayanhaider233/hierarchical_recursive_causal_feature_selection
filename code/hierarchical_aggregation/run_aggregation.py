import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

METHYLATION_PATH = DATA / "methylation"
METADATA_PATH = DATA / "metadata"

INTERMEDIATE_PATH = DATA / "intermediate"
INTERMEDIATE_PATH.mkdir(parents=True, exist_ok=True)

AGGREGATION_PATH = ROOT / "code" / "hierarchical_aggregation"

train_cpg = METHYLATION_PATH / "train_cpg_mval_matrix_corrected.csv.csv"
test_cpg  = METHYLATION_PATH / "test_cpg_mval_matrix_corrected.csv.csv"
train_meta = METADATA_PATH / "train_metadata.csv"

dmr_rds = INTERMEDIATE_PATH / "dmr_groups.rds"
dmr_map = INTERMEDIATE_PATH / "dmr_cpg_map.csv"

train_dmr = METHYLATION_PATH / "train_dmr_matrix.csv"
test_dmr  = METHYLATION_PATH / "test_dmr_matrix.csv"

subprocess.run([
    "Rscript",
    str(AGGREGATION_PATH / "dmr_detection.R"),
    str(train_cpg),
    str(train_meta),
    str(dmr_rds),
    str(dmr_map)
], check=True)

subprocess.run([
    "Rscript",
    str(AGGREGATION_PATH / "dmr_aggregation.R"),
    str(train_cpg),
    str(dmr_rds),
    str(train_dmr)
], check=True)

subprocess.run([
    "Rscript",
    str(AGGREGATION_PATH / "dmr_aggregation.R"),
    str(test_cpg),
    str(dmr_rds),
    str(test_dmr)
], check=True)
