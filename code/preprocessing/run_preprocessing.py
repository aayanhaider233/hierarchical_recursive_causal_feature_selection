# --------------------------------------------------
# Imports
# --------------------------------------------------

import metadata_harmonisation as meta
import probe_quality_control as pqc
import epigenetic_age as epiage
import probe_filter as pfilter
import concatenate_datasets as cdata
import methylation_beta_to_M as btm
import batch_effect_correction as bcorr
import partition_data as partition
import epigenetic_age_acceleration as eaa

from pathlib import Path 


# --------------------------------------------------
# Dataset directories
# --------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent 
DATA = ROOT / "data"
ANNOTATIONS = ROOT / "annotations" / "methylation"
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
ANNOTATION_PATHS = {
    "manifest_path" : ANNOTATIONS / "humanmethylation450_15017482_v1-2.csv",
    "cross_reactive_probes_path" : ANNOTATIONS / "48639-non-specific-probes-Illumina450k.csv",
    "multimap_probes_path" : ANNOTATIONS / "HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.csv" 
}


# --------------------------------------------------
# Save directories
# --------------------------------------------------

METADATA_SAVE = DATA / "metadata"
METADATA_SAVE.mkdir(parents=True, exist_ok=True)
METHYLATION_SAVE = DATA / "methylation"
METHYLATION_SAVE.mkdir(parents=True, exist_ok=True)
INTERMEDIATE_SAVE = DATA / "intermediate"
INTERMEDIATE_SAVE.mkdir(parents=True, exist_ok=True)


def run_preprocessing():

    print("Preprocessing started.")
    # --------------------------------------------------
    # 1. Metadata harmonisation
    # --------------------------------------------------

    metadata_dfs = {}

    for batch, path in METADATA_PATHS.items():
        metadata_dfs[batch] = meta.parse_metadata(
            metadata_path=path, 
            batch=batch
        )

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
            meta_df=metadata_dfs[batch],
            **params
        )

    # --------------------------------------------------
    # 2. Non-CpG probe filter & probe QC
    # --------------------------------------------------

    methylation_dfs = {}

    for batch, path in METHYLATION_MATRIX_PATHS.items():
        
        methylation_dfs[batch] = pqc.remove_non_cpg_probes(
            methylation_matrix_path=path, 
            batch=batch
        )
        
        methylation_dfs[batch], metadata_dfs[batch] = pqc.filter_by_metadata(
            methylation_df=methylation_dfs[batch],
            metadata_df=metadata_dfs[batch]
        )
        
        methylation_dfs[batch] = pqc.filter_by_detection_pvalue(
            methylation_df=methylation_dfs[batch]
        )
        
        methylation_dfs[batch] = pqc.clean_cpg_probes(
            methylation_df=methylation_dfs[batch]
        )
        
        assert set(methylation_dfs[batch]["sample_id"]) == set(metadata_dfs[batch]["sample_id"])

    # --------------------------------------------------
    # 3. Epigenetic age calculation
    # --------------------------------------------------

    epi_dfs = {}

    for batch, df in methylation_dfs.items():
        epi_dfs[batch] = epiage.generate_epigenetic_age(
            methylation_df=df
        )

    for batch, df in metadata_dfs.items():
        metadata_dfs[batch] = epiage.add_epiage_to_metadata(
            metadata_df=metadata_dfs[batch], 
            epi_df=epi_dfs[batch]
        )


    # --------------------------------------------------
    # 4. Probe filter
    # --------------------------------------------------

    for batch, df in methylation_dfs.items():
        methylation_dfs[batch] = pfilter.filter_probes(
            methylation_df=methylation_dfs[batch], 
            **ANNOTATION_PATHS
        )


    # --------------------------------------------------
    # 5. Dataset concatenation & batch info generation
    # --------------------------------------------------

    metadata_full = cdata.concatenate_metadata(
        *metadata_dfs.values()
    )

    batch_info = cdata.generate_batch_info(
        metadata_df=metadata_full
    )

    methylation_matrix_full = cdata.concatenate_methylation_matrices(
        *methylation_dfs.values()
    )

    assert set(methylation_matrix_full["sample_id"]) == set(metadata_full["sample_id"])


    # --------------------------------------------------
    # 6. Beta- to M-value conversion
    # --------------------------------------------------

    methylation_mval_matrix_full = btm.beta_to_m(
        methylation_df=methylation_matrix_full
    )


    # --------------------------------------------------
    # 7. Batch effect correction
    # --------------------------------------------------

    methylation_mval_matrix_corrected = bcorr.batch_effect_correction(
        methylation_df=methylation_mval_matrix_full,
        metadata_df=metadata_full,
        batch_info=batch_info
    )


    # --------------------------------------------------
    # 8. Train-test split
    # --------------------------------------------------

    train_methylation, test_methylation, train_metadata, test_metadata, train_batch, test_batch = partition.partition_data(
        methylation_df=methylation_mval_matrix_corrected, 
        metadata_df=metadata_full, 
        batch_info=batch_info
    )


    # --------------------------------------------------
    # 9. Epigenetic age acceleration computation
    # --------------------------------------------------

    train_metadata, test_metadata, train_eaa, test_eaa = eaa.compute_eaa(
        train_metadata=train_metadata, 
        test_metadata=test_metadata
    )


    # --------------------------------------------------
    # Save
    # --------------------------------------------------

    batch_info.to_csv(
        INTERMEDIATE_SAVE / "batch_info.csv",
        index=False
    )

    methylation_mval_matrix_corrected.to_csv(
        METHYLATION_SAVE / "cpg_mval_matrix_full_corrected.csv",
        index=False
    )

    train_methylation.to_csv(
        METHYLATION_SAVE / "train_cpg_mval_matrix_corrected.csv", 
        index=False
    )
    test_methylation.to_csv(
        METHYLATION_SAVE / "test_cpg_mval_matrix_corrected.csv", 
        index=False
    )

    train_metadata.to_csv(
        METADATA_SAVE / "train_metadata.csv", 
        index=False
    )
    test_metadata.to_csv(
        METADATA_SAVE / "test_metadata.csv", 
        index=False
    )

    train_eaa.to_csv(
        METADATA_SAVE / "train_eaa.csv", 
        index=False
    )
    test_eaa.to_csv(
        METADATA_SAVE / "test_eaa.csv", 
        index=False
    )

    print("Preprocessing complete.")


if __name__ == "__main__":
    run_preprocessing()