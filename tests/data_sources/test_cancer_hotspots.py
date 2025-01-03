"""Module for testing cancer hotspots"""

import pytest

from evidence.data_sources import CancerHotspots


@pytest.fixture(scope="module")
def cancer_hotspots(evidence_data_dir):
    """Create test fixture for cancer hotspots class"""
    globbed = (evidence_data_dir / "cancer_hotspots").glob("cancer_hotspots_*.json")
    return CancerHotspots(transformed_data_path=sorted(globbed)[-1])


@pytest.fixture(scope="module")
def braf_v600e():
    """Create test fixture for braf v600e"""
    return {
        "variation": "BRAF V600E",
        "codon": "V600",
        "mutation": "V600E",
        "q_value": 0.0,
        "observations": 833,
        "total_observations": 897,
    }


def check_source_meta(response):
    """Check that source meta is correct"""
    assert response["source_meta_"]["label"] == "Cancer Hotspots"
    assert response["source_meta_"]["version"] == "2"


def test_mutation_hotspots(cancer_hotspots, braf_v600e):
    """Test that mutation_hotspots method works correctly."""
    resp = cancer_hotspots.mutation_hotspots(
        "ga4gh:VA.urVNVupVvzqE57gZFBu6vsAFScIBNojG"
    )
    assert resp.data == {
        "variation": "BRAF N486_A489delinsK",
        "codon": "486-494",
        "mutation": "N486_A489delinsK",
        "q_value": 3.9500653726776106e-09,
        "observations": 1,
        "total_observations": 7,
    }

    resp = cancer_hotspots.mutation_hotspots(
        "ga4gh:VA.iKBIyhCE-3Cnh5XJpVNZ4t0Ju8Ov8Ine"
    )
    assert resp.data == {
        "variation": "BRAF N486_P490del",
        "codon": "486-494",
        "mutation": "N486_P490del",
        "q_value": 3.9500653726776106e-09,
        "observations": 3,
        "total_observations": 7,
    }

    resp = cancer_hotspots.mutation_hotspots(
        "ga4gh:VA.uIzxca0MnfZ1ni3NYNkvNh6Isw-5QXG9"
    )
    assert resp.data == {
        "variation": "TP53 I255del",
        "codon": "229-292",
        "mutation": "I255del",
        "q_value": 1.46521691051924e-64,
        "observations": 9,
        "total_observations": 76,
    }

    resp = cancer_hotspots.mutation_hotspots(
        "ga4gh:VA.AM5oy_X7xBxJmWyUPJ0ZSLpRaTgu2kHt"
    )
    assert resp.data == {
        "variation": "CEBPA T310_Q311insKQNP",
        "codon": "307-311",
        "mutation": "T310_Q311insKQNP",
        "q_value": 0.0050613902181432,
        "observations": 1,
        "total_observations": 4,
    }

    resp = cancer_hotspots.mutation_hotspots(
        vrs_variation_id="ga4gh:VA.j4XnsLZcdzDIYa5pvvXM7t1wn9OITr0L"
    ).model_dump(by_alias=True)
    assert resp["_id"] == "normalize.evidence:92f3db383a79d855323a71d65d860ec3"
    assert resp["data"] == braf_v600e
    check_source_meta(resp)

    # invalid vrs_variation_id
    resp = cancer_hotspots.mutation_hotspots(
        vrs_variation_id="ga4ghVA8JkgnqIgYqufNl-OV_hpRG_aWF9UFQCE"
    ).model_dump(by_alias=True)
    assert resp["_id"] is None
    assert resp["data"] == {}
    check_source_meta(resp)
