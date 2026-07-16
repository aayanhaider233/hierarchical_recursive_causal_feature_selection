import metadata_harmonisation_01 as meta
import probe_quality_control_02 as pqc
import epigenetic_age_03 as epiage
import probe_filter_04 as pfilter
import concatenate_datasets_05 as cdata
import methylation_beta_to_M_06 as btm
import partition_data_07 as partition
import batch_effect_correction_08 as bcorr
import epigenetic_age_acceleration_09 as eaa

from pathlib import Path 

ROOT = Path(__file__).resolve().parent.parent 
DATA = ROOT / "data"
ANNOTATIONS = ROOT / "annotations"

RAW = DATA / "raw"

BATCHES = [
    "GSE80417", 
    "GSE84727", 
    "GSE147221", 
    "GSE152027"
]

METADATA_PATHS = {batch : RAW / f"{batch}_series_matrix.txt.gz" for batch in BATCHES}
METHYLATION_MATRIX_PATHS = {
    "GSE80417" : RAW / "GSE80417_normalizedBetas.csv.gz",
    "GSE84727" : RAW / "GSE84727_normalisedBetas.csv.gz",
    "GSE147221" : RAW / "GSE147221_Dublin_blood_processed_signals.csv.gz",
    "GSE152027" : RAW / "GSE152027_IOP_processed_signals.csv.gz" 
}


metadata_dfs = {}

for batch, path in METADATA_PATHS.items():
    metadata_dfs[batch] = meta.parse_metadata(path, batch)

cleaning_params = {
    "GSE80417": {
        "control_label": "1",
        "case_label": "2"
    },
    "GSE84727": {
        "control_label": "1",
        "case_label": "2"
    },
    "GSE147221": {},
    "GSE152027": {
        "control_label": "CON",
        "case_label": "SCZ"
    }
}

for batch, params in cleaning_params.items():
    metadata_dfs[batch] = meta.clean_metadata(
        metadata_dfs[batch],
        **params
    )


methylation_dfs = {}

for batch, path in METHYLATION_MATRIX_PATHS.items():
    
    methylation_dfs[batch] = pqc.remove_non_cpg_probes(path, batch)
    
    methylation_dfs[batch], metadata_dfs[batch] = pqc.filter_by_metadata(
        methylation_dfs[batch],
        metadata_dfs[batch]
    )
    
    methylation_dfs[batch] = pqc.filter_by_detection_pvalue(methylation_dfs[batch])
    
    methylation_dfs[batch] = pqc.clean_cpg_probes(methylation_dfs[batch])
    
    assert set(methylation_dfs[batch]["sample_id"]) == set(metadata_dfs[batch]["sample_id"])