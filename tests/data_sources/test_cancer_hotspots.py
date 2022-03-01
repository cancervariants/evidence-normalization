"""Module for testing cancer hotspots"""
import pytest

from evidence.data_sources import CancerHotspots


@pytest.fixture(scope="module")
def cancer_hotspots():
    """Create test fixture for cancer hotspots class"""
    return CancerHotspots()


@pytest.fixture(scope="module")
def braf_v600e():
    """Create test fixture for braf v600e"""
    return {
        "codon": "V600",
        "mutation": "V600E",
        "q_value": 0.0,
        "observations": 833,
        "total_observations": 897
    }


def check_source_meta(response):
    """Check that source meta is correct"""
    assert response["source_meta_"]["label"] == "Cancer Hotspots"
    assert response["source_meta_"]["version"] == "2"


def test_query_snv_hotspots(cancer_hotspots, braf_v600e):
    """Test that query_snv_hotspots method works correctly."""
    resp = cancer_hotspots.query_snv_hotspots(
        vrs_variation_id="ga4gh:VA.8JkgnqIgYqufNl-OV_hpRG_aWF9UFQCE"
    )
    assert resp == braf_v600e

    resp = cancer_hotspots.query_snv_hotspots(
        vrs_variation_id="ga4gh:VA.8JkgnqIgYqufNl-OV_hpRG_aWF9UFQCEssdfs"
    )
    assert resp is None


def test_query_indel_hotspots(cancer_hotspots):
    """Test that query_indel_hotspots method works correctly."""
    # TODO: BRAF T599dup not supported yet in variation-normalizer
    # resp = cancer_hotspots.query_indel_hotspots()
    # assert resp == {
    #     "codon": "592-604",
    #     "mutation": "T599dup",
    #     "q_value": 6.86065772143357e-11,
    #     "observations": 3,
    #     "total_observations": 8
    # }

    # BRAF N486_A489delinsK
    resp = cancer_hotspots.query_indel_hotspots(
        "ga4gh:VA.elb_blzX0l5YHbSEHHdFYIoJmAEbJ2yv")
    assert resp == {
        "codon": "486-494",
        "mutation": "N486_A489delinsK",
        "q_value": 3.950065372677611e-09,
        "observations": 1,
        "total_observations": 7
    }

    # BRAF N486_P490del
    resp = cancer_hotspots.query_indel_hotspots(
        "ga4gh:VA.VNzgFrkk2GFvUFAm7gtB6BxdN5H8uuKC")
    assert resp == {
        "codon": "486-494",
        "mutation": "N486_P490del",
        "q_value": 3.950065372677611e-09,
        "observations": 3,
        "total_observations": 7
    }

    # TP53 I255del
    resp = cancer_hotspots.query_indel_hotspots(
        "ga4gh:VA._vUTgjJ8TUV6lh0ol0fHkf3ptwOWYU_W")
    assert resp == {
        "codon": "229-292",
        "mutation": "I255del",
        "q_value": 1.46521691051924e-64,
        "observations": 9,
        "total_observations": 76
    }

    # CEBPA T310_Q311insKQNP
    resp = cancer_hotspots.query_indel_hotspots(
        "ga4gh:VA.pDj9pPKvaZ-BjhdtZ2wi717V4ZAuf-OD")
    assert resp == {
        "codon": "307-311",
        "mutation": "T310_Q311insKQNP",
        "q_value": 0.0050613902181432,
        "observations": 1,
        "total_observations": 4
    }


def test_hotspot_data(cancer_hotspots, braf_v600e):
    """Test that hotspot_data method works correctly."""
    resp = cancer_hotspots.mutation_hotspots(
        so_id="SO:0001606", vrs_variation_id="ga4gh:VA.8JkgnqIgYqufNl-OV_hpRG_aWF9UFQCE"
    ).dict(by_alias=True)
    assert resp["_id"] == "normalize.evidence:f14dcca46895ceda70ea901452dfe1d4"
    assert resp["data"] == braf_v600e
    check_source_meta(resp)

    # invalid vrs_variation_id
    resp = cancer_hotspots.mutation_hotspots(
        so_id="SO:0001606", vrs_variation_id="ga4ghVA8JkgnqIgYqufNl-OV_hpRG_aWF9UFQCE"
    ).dict(by_alias=True)
    assert resp["_id"] is None
    assert resp["data"] == dict()
    check_source_meta(resp)

    # invalid so_id
    resp = cancer_hotspots.mutation_hotspots(
        so_id="SO0001606", vrs_variation_id="ga4gh:VA.8JkgnqIgYqufNl-OV_hpRG_aWF9UFQCE"
    ).dict(by_alias=True)
    assert resp["_id"] is None
    assert resp["data"] == dict()
    check_source_meta(resp)
