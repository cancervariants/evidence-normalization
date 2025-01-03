"""Module for testing cbioportal"""

import pytest

from evidence.data_sources import CBioPortal


@pytest.fixture(scope="module")
def cbioportal(evidence_data_dir):
    """Create test fixture for cbioportal class"""
    cbioportal_data_dir = evidence_data_dir / "cbioportal"
    return CBioPortal(
        transformed_mutations_data_path=cbioportal_data_dir
        / "msk_impact_2017_mutations.csv",
        transformed_case_lists_data_path=cbioportal_data_dir
        / "msk_impact_2017_case_lists.csv",
    )


@pytest.fixture(scope="module")
def braf():
    """Create test fixture for braf"""
    return {
        "Colorectal Cancer": {
            "count": 119,
            "total": 1006,
            "percent_altered": 119 / 1006 * 100,
        },
        "Anal Cancer": {"count": 1, "total": 32, "percent_altered": 1 / 32 * 100},
        "Skin Cancer, Non-Melanoma": {
            "count": 4,
            "total": 148,
            "percent_altered": 4 / 148 * 100,
        },
        "Miscellaneous Brain Tumor": {"count": 0, "total": 6, "percent_altered": 0.0},
        "Bladder Cancer": {
            "count": 17,
            "total": 423,
            "percent_altered": 17 / 423 * 100,
        },
        "Small Bowel Cancer": {
            "count": 6,
            "total": 35,
            "percent_altered": 6 / 35 * 100,
        },
        "Renal Cell Carcinoma": {"count": 0, "total": 360, "percent_altered": 0.0},
        "Uterine Sarcoma": {"count": 2, "total": 93, "percent_altered": 2 / 93 * 100},
        "Embryonal Tumor": {  # Peripheral Nervous System on cBioPortal UI
            "count": 2,
            "total": 86,
            "percent_altered": 2 / 86 * 100,
        },
        "Retinoblastoma": {"count": 0, "total": 4, "percent_altered": 0.0},
        "Endometrial Cancer": {
            "count": 9,
            "total": 218,
            "percent_altered": 9 / 218 * 100,
        },
        "Glioma": {"count": 17, "total": 553, "percent_altered": 17 / 553 * 100},
        "Melanoma": {"count": 114, "total": 350, "percent_altered": 114 / 350 * 100},
        "Cervical Cancer": {"count": 1, "total": 50, "percent_altered": 1 / 50 * 100},
        "Soft Tissue Sarcoma": {
            "count": 6,
            "total": 435,
            "percent_altered": 6 / 435 * 100,
        },
        "NA": {"count": 1, "total": 57, "percent_altered": 1 / 57 * 100},
        "Multiple Myeloma": {"count": 0, "total": 1, "percent_altered": 0},
        "Small Cell Lung Cancer": {
            "count": 3,
            "total": 82,
            "percent_altered": 3 / 82 * 100,
        },
        "Miscellaneous Neuroepithelial Tumor": {
            "count": 0,
            "total": 9,
            "percent_altered": 0.0,
        },
        "Pancreatic Cancer": {
            "count": 13,
            "total": 502,
            "percent_altered": 13 / 502 * 100,
        },
        "Salivary Gland Cancer": {
            "count": 2,
            "total": 114,
            "percent_altered": 2 / 114 * 100,
        },
        "Penile Cancer": {"count": 1, "total": 7, "percent_altered": 1 / 7 * 100},
        "Hodgkin Lymphoma": {"count": 0, "total": 5, "percent_altered": 0.0},
        "Pheochromocytoma": {"count": 0, "total": 4, "percent_altered": 0.0},
        "Breast Cancer": {"count": 4, "total": 1324, "percent_altered": 4 / 1324 * 100},
        "Breast Sarcoma": {"count": 0, "total": 13, "percent_altered": 0.0},
        "Adrenocortical Carcinoma": {
            "count": 0,
            "total": 26,  # UI says 25, but API combines Adrenocortical Adenoma (1) and  Adrenocortical Carcinoma (25), making it 26
            "percent_altered": 0.0,
        },
        "Germ Cell Tumor": {"count": 2, "total": 288, "percent_altered": 2 / 288 * 100},
        "Head and Neck Cancer": {
            "count": 1,
            "total": 173,
            "percent_altered": 1 / 173 * 100,
        },
        "Hepatobiliary Cancer": {
            "count": 11,
            "total": 338,
            "percent_altered": 11 / 338 * 100,
        },
        "Esophagogastric Cancer": {
            "count": 5,
            "total": 341,
            "percent_altered": 5 / 341 * 100,
        },
        "Prostate Cancer": {
            "count": 12,
            "total": 717,
            "percent_altered": 12 / 717 * 100,
        },
        "Ampullary Carcinoma": {
            "count": 0,
            "total": 25,  # UI says 9, but API combines Ampullary Cancer (16) and Ampullary Carcinoma (9), making it 25
            "percent_altered": 0.0,
        },
        "Histiocytosis": {"count": 6, "total": 22, "percent_altered": 6 / 22 * 100},
        "Gastrointestinal Neuroendocrine Tumor": {
            "count": 4,
            "total": 46,
            "percent_altered": 4 / 46 * 100,
        },
        "Non-Small Cell Lung Cancer": {
            "count": 87,
            "total": 1668,
            "percent_altered": 87 / 1668 * 100,
        },
        "Non-Hodgkin Lymphoma": {
            "count": 5,
            "total": 165,
            "percent_altered": 5 / 165 * 100,
        },
        "Cancer of Unknown Primary": {
            "count": 10,
            "total": 184,
            "percent_altered": 10 / 184 * 100,
        },
        "Mesothelioma": {"count": 0, "total": 107, "percent_altered": 0.0},
        "Thyroid Cancer": {
            "count": 88,
            "total": 231,
            "percent_altered": 88 / 231 * 100,
        },
        "Pineal Tumor": {"count": 0, "total": 3, "percent_altered": 0.0},
        "Wilms Tumor": {"count": 0, "total": 2, "percent_altered": 0.0},
        "Thymic Tumor": {"count": 0, "total": 18, "percent_altered": 0.0},
        "Appendiceal Cancer": {
            "count": 1,
            "total": 79,
            "percent_altered": 1 / 79 * 100,
        },
        "Gestational Trophoblastic Disease": {
            "count": 0,
            "total": 3,
            "percent_altered": 0.0,
        },
        "Leukemia": {"count": 1, "total": 4, "percent_altered": 1 / 4 * 100},
        "Mastocytosis": {"count": 0, "total": 1, "percent_altered": 0.0},
        "Vaginal Cancer": {"count": 0, "total": 4, "percent_altered": 0.0},
        "Sex Cord Stromal Tumor": {"count": 0, "total": 19, "percent_altered": 0.0},
        "Nerve Sheath Tumor": {"count": 0, "total": 16, "percent_altered": 0.0},
        "Sellar Tumor": {"count": 0, "total": 5, "percent_altered": 0.0},
        "Bone Cancer": {"count": 1, "total": 134, "percent_altered": 1 / 134 * 100},
        "Gastrointestinal Stromal Tumor": {
            "count": 2,
            "total": 137,
            "percent_altered": 2 / 137 * 100,
        },
        "CNS Cancer": {"count": 1, "total": 48, "percent_altered": 1 / 48 * 100},
        "Ovarian Cancer": {"count": 3, "total": 224, "percent_altered": 3 / 224 * 100},
    }


def test_get_mutation_data(cbioportal, braf):
    """Test that get_mutation_data method works correctly."""
    resp = cbioportal.cancer_types_summary("braf").model_dump(by_alias=True)
    data = resp["data"]
    assert data.keys() == braf.keys()

    total_cancer_type_count = sum(v["count"] for v in data.values())
    assert total_cancer_type_count == 562

    assert data == braf
    assert resp["_id"] == "normalize.evidence:0d76601b070d3c9f36e9c483454adf48"
    assert resp["source_meta_"]["label"] == "cBioPortal"
    assert resp["source_meta_"]["version"] == "msk_impact_2017"

    resp = cbioportal.cancer_types_summary("dummy").model_dump()
    assert resp["data"] == {}
    assert resp["source_meta_"]["label"] == "cBioPortal"
    assert resp["source_meta_"]["version"] == "msk_impact_2017"
