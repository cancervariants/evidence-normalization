"""Module for testing gnomad"""
import time

import pytest

from evidence.data_sources import gnomAD


@pytest.fixture(scope="module")
def gnomad():
    """Create test fixture for gnomad class"""
    return gnomAD()


@pytest.fixture(scope="module")
def braf():
    """Create test fixture for BRAF"""
    return {
        "variant": "7-140453136-A-T",
        "assembly": "GRCh37",
        "total_observations": {
            "allele_count": 1,
            "allele_number": 251260,
            "decimal": "0.000003980"
        },
        "max_pop_freq": {
            "population": "South Asian",
            "allele_count": 1,
            "allele_number": 30612,
            "decimal": "0.000032667"
        },
        "gnomad_url": "https://gnomad.broadinstitute.org/variant/7-140453136-A-T?dataset=gnomad_r2_1"  # noqa: E501
    }


@pytest.fixture(scope="module")
def tpm3():
    """Create test fixture for TPM3"""
    return {
        "variant": "1-154157570-C-T",
        "assembly": "GRCh38",
        "total_observations": {
            "allele_count": 1,
            "allele_number": 152152,
            "decimal": "0.000006572"
        },
        "max_pop_freq": {
            "population": "African/African American",
            "allele_count": 1,
            "allele_number": 41434,
            "decimal": "0.000024135"
        },
        "gnomad_url": "https://gnomad.broadinstitute.org/variant/1-154157570-C-T?dataset=gnomad_r3"  # noqa: E501
    }


@pytest.fixture(scope="module")
def egfr_37():
    """Create test fixture for EGFR on 37"""
    return {
        "variant": "7-55209955-G-A",
        "assembly": "GRCh37",
        "total_observations": {
            "allele_count": 16,
            "allele_number": 282848,
            "decimal": "0.000056567"
        },
        "max_pop_freq": {
            "population": "Latino/Admixed American",
            "allele_count": 14,
            "allele_number": 35436,
            "decimal": "0.000395078"
        },
        "gnomad_url": "https://gnomad.broadinstitute.org/variant/7-55209955-G-A?dataset=gnomad_r2_1"  # noqa: E501
    }


@pytest.fixture(scope="module")
def egfr_38():
    """Create test fixture for EGFR on 38"""
    return {
        "variant": "7-55142262-G-A",
        "assembly": "GRCh38",
        "total_observations": {
            "allele_count": 7,
            "allele_number": 152192,
            "decimal": "0.000045995"
        },
        "max_pop_freq": {
            "population": "Latino/Admixed American",
            "allele_count": 3,
            "allele_number": 15278,
            "decimal": "0.000196361"
        },
        "gnomad_url": "https://gnomad.broadinstitute.org/variant/7-55142262-G-A?dataset=gnomad_r3"  # noqa: E501
    }


@pytest.fixture(scope="module")
def variant_id():
    """Create test fixture for variant_id"""
    return {
        "data": {
            "variant_search": [
                {
                    "variant_id": "7-140453136-A-T"
                }
            ],
            "dataset": "gnomad_r2_1"
        }
    }


def test_liftover_37_to_38(gnomad):
    """Test that liftover_37_to_38 method works correctly"""
    resp = gnomad.liftover_37_to_38("7-140453136-A-T").dict()
    assert resp["data"] == {
        "gnomad_variant_id": "7-140753336-A-T",
        "reference_genome": "GRCh38",
        "dataset": "gnomad_r3"
    }

    resp = gnomad.liftover_37_to_38("7-140753336-A-T").dict()
    assert resp["data"] == dict()
    time.sleep(5)


def test_liftover_38_to_37(gnomad):
    """Test that liftover_38_to_37 method works correctly"""
    resp = gnomad.liftover_38_to_37("7-140753336-A-T").dict()
    assert resp["data"] == {
        "gnomad_variant_id": "7-140453136-A-T",
        "reference_genome": "GRCh37",
        "dataset": "gnomad_r2_1"
    }

    resp = gnomad.liftover_38_to_37("7-140453136-A-T").dict()
    assert resp["data"] == dict()
    time.sleep(5)


def test_clinvar_variation_id(gnomad):
    """Test that clinvar_variation_id method works correctly."""
    resp = gnomad.clinvar_variation_id("7-140453136-A-T", "GRCh37").dict()
    assert resp["data"]["clinvar_variation_id"] == "13961"

    resp = gnomad.clinvar_variation_id("7-140453136-A-T").dict()
    assert resp["data"]["clinvar_variation_id"] == "13961"

    resp = gnomad.clinvar_variation_id("7-140753336-A-T", "GRCh38").dict()
    assert resp["data"]["clinvar_variation_id"] == "13961"

    resp = gnomad.clinvar_variation_id("7-140453136-A-T", "GRCh38").dict()
    assert resp["data"] == dict()
    time.sleep(5)


def test_variant_id_to_gnomad_id(gnomad, variant_id):
    """Test that variant_id_to_gnomad_id method works correctly."""
    resp = gnomad.variant_id_to_gnomad_id("13961")
    assert resp == variant_id

    resp = gnomad.variant_id_to_gnomad_id("7-140453136-A-T")
    assert resp == variant_id

    resp = gnomad.variant_id_to_gnomad_id("7-140453136-A-T", "GRCh37")
    assert resp == variant_id

    resp = gnomad.variant_id_to_gnomad_id("rs113488022")
    assert resp == variant_id

    resp = gnomad.variant_id_to_gnomad_id("CA123643")
    assert resp == variant_id

    # clinvar id not a string
    resp = gnomad.variant_id_to_gnomad_id(13961)
    assert resp["errors"]

    # Invalid query
    resp = gnomad.variant_id_to_gnomad_id("fake")
    assert resp["errors"]
    time.sleep(5)


def test_variant(gnomad):
    """Test that variant method works correctly."""
    resp = gnomad.variant("7-140453136-A-T", "GRCh37")
    assert resp["data"]["variant"]["exome"]
    assert resp["data"]["variant"]["genome"] is None

    resp = gnomad.variant("7-140426282-C-A", "GRCh37")
    assert resp["data"]["variant"]["exome"]
    assert resp["data"]["variant"]["genome"]

    # wrong dataset
    resp = gnomad.variant("7-140453136-A-T", "GRCh38")
    assert resp["errors"]

    # invalid variant_id
    resp = gnomad.variant("dummy", "GRCh38")
    assert resp["errors"]
    time.sleep(10)


def test_frequency_data(gnomad, braf, tpm3, egfr_37, egfr_38):
    """Test that frequency_data method works correctly."""
    # Asssembly 37
    resp = gnomad.frequency_data("13961").dict()
    assert resp["data"] == braf
    assert resp["source_meta_"]["label"] == "gnomAD"
    assert resp["source_meta_"]["version"] == "gnomad_r2_1"
    time.sleep(10)

    resp = gnomad.frequency_data("7-140453136-A-T").dict()
    assert resp["data"] == braf
    assert resp["source_meta_"]["label"] == "gnomAD"
    assert resp["source_meta_"]["version"] == "gnomad_r2_1"
    time.sleep(10)

    resp = gnomad.frequency_data("rs113488022").dict()
    assert resp["data"] == braf
    assert resp["source_meta_"]["label"] == "gnomAD"
    assert resp["source_meta_"]["version"] == "gnomad_r2_1"
    time.sleep(10)

    resp = gnomad.frequency_data("CA123643").dict()
    assert resp["data"] == braf
    assert resp["source_meta_"]["label"] == "gnomAD"
    assert resp["source_meta_"]["version"] == "gnomad_r2_1"
    time.sleep(10)

    # Assembly 38
    resp = gnomad.frequency_data("1-154157570-C-T").dict()
    assert resp["data"] == tpm3
    assert resp["source_meta_"]["label"] == "gnomAD"
    assert resp["source_meta_"]["version"] == "gnomad_r3"
    time.sleep(10)

    resp = gnomad.frequency_data("CA889563985").dict()
    assert resp["data"] == tpm3
    assert resp["source_meta_"]["label"] == "gnomAD"
    assert resp["source_meta_"]["version"] == "gnomad_r3"
    time.sleep(10)

    resp = gnomad.frequency_data("rs1443960419").dict()
    assert resp["data"] == tpm3
    assert resp["source_meta_"]["label"] == "gnomAD"
    assert resp["source_meta_"]["version"] == "gnomad_r3"
    time.sleep(10)

    # both exome and genome
    resp = gnomad.frequency_data("7-55209955-G-A").dict()
    assert resp["data"] == egfr_37
    assert resp["source_meta_"]["label"] == "gnomAD"
    assert resp["source_meta_"]["version"] == "gnomad_r2_1"
    time.sleep(10)

    resp = gnomad.frequency_data("7-55142262-G-A").dict()
    assert resp["data"] == egfr_38
    assert resp["source_meta_"]["label"] == "gnomAD"
    assert resp["source_meta_"]["version"] == "gnomad_r3"
    time.sleep(10)

    resp = gnomad.frequency_data("rs368926971").dict()
    assert resp["data"] == egfr_38
    assert resp["source_meta_"]["label"] == "gnomAD"
    assert resp["source_meta_"]["version"] == "gnomad_r3"
    time.sleep(10)

    resp = gnomad.frequency_data("CA4265128").dict()
    assert resp["data"] == egfr_38
    assert resp["source_meta_"]["label"] == "gnomAD"
    assert resp["source_meta_"]["version"] == "gnomad_r3"
    time.sleep(10)

    resp = gnomad.frequency_data("fake").dict()
    assert resp["data"] == dict()
    assert resp["source_meta_"] == {"label": "gnomAD", "version": None}
