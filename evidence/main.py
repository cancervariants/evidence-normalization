"""Main application for FastAPI."""
from typing import Dict, Optional
from urllib.parse import unquote

from fastapi import FastAPI, Query
from fastapi.openapi.utils import get_openapi

from evidence.data_sources import CancerHotspots, CBioPortal, GnomAD
from evidence.version import __version__
from evidence.schemas import ReferenceGenome, Response

gnomad = GnomAD()
cbioportal = CBioPortal()
cancer_hotspots = CancerHotspots()
app = FastAPI(docs_url="/evidence", openapi_url="/evidence/openapi.json",
              swagger_ui_parameters={"tryItOutEnabled": True})


def custom_openapi() -> Dict:
    """Generate custom fields for OpenAPI response"""
    if app.openapi_schema:
        return app.openapi_schema
    openapi_schema = get_openapi(
        title="The VICC Evidence Normalizer",
        version=__version__,
        description="Services for normalizing evidence.",
        routes=app.routes
    )
    openapi_schema["info"]["contact"] = {
        "name": "Alex H. Wagner",
        "email": "Alex.Wagner@nationwidechildrens.org",
        "url": "https://www.nationwidechildrens.org/specialties/institute-for-genomic-medicine/research-labs/wagner-lab"  # noqa: E501
    }
    app.openapi_schema = openapi_schema
    return app.openapi_schema


app.openapi = custom_openapi
RESPONSE_DESCRIPTION = "A response to a validly-formed query."


@app.get("/evidence/cancer_hotspots/mutation_hotspots",
         summary="Given variant data, return cancer hotspots mutation data.",
         response_description=RESPONSE_DESCRIPTION,
         description="Return mutation hotspots data for variant.",
         response_model=Response)
def get_cancer_hotspots(
    so_id: str = Query("SO:0001606", enum=["SO:0001606", "SO:0001017"],
                       description="The structural type of the variant"),
    vrs_variation_id: str = Query(..., description="The VRS digest for the variation")
) -> Response:
    """Return cancer types with gene mutations.

    :param str so_id: The structural type of the variant
    :param str vrs_variation_id: The VRS digest for the variation

    :return: Cancer Hotspot data
    """
    return cancer_hotspots.mutation_hotspots(unquote(so_id), unquote(vrs_variation_id))


@app.get("/evidence/cbioportal/cancer_types_summary",
         summary="Given HGNC gene symbol, return cancer types summary.",
         response_description=RESPONSE_DESCRIPTION,
         description="Return cancer types with gene mutations.",
         response_model=Response)
def get_cancer_types_summary(
    hgnc_symbol: str = Query(..., description="HGNC gene symbol.")
) -> Response:
    """Return cancer types with gene mutations.

    :param str hgnc_symbol: HGNC gene symbol
    :return: Return cancer types with `hgnc_symbol` mutations
    """
    return cbioportal.cancer_types_summary(hgnc_symbol)


@app.get("/evidence/gnomad/liftover/38_to_37",
         summary="Liftover gnomad variant id from GRCh38 to GRCh37",
         response_description=RESPONSE_DESCRIPTION,
         description="Return GRCh37 gnomad variant id",
         response_model=Response)
def gnomad_38_to_37(
    gnomad_variant_id: str = Query(..., decsription="gnomAD variant ID on 38 assembly")
) -> Response:
    """Given 38 gnomad variant id, liftover to 37

    :param str gnomad_variant_id: gnomAD variant ID on 38 assembly
    :return: gnomad 37 assembly data
    """
    return gnomad.liftover_38_to_37(gnomad_variant_id)


@app.get("/evidence/gnomad/liftover/37_to_38",
         summary="Liftover gnomad variant id from GRCh37 to GRCh38",
         response_description=RESPONSE_DESCRIPTION,
         description="Return GRCh38 gnomad variant id",
         response_model=Response)
def gnomad_37_to_38(
    gnomad_variant_id: str = Query(..., decsription="gnomAD variant ID on 37 assembly")
) -> Response:
    """Given 37 gnomad variant id, liftover to 38

    :param str gnomad_variant_id: gnomAD variant ID on 37 assembly
    :return: gnomad 38 assembly data
    """
    return gnomad.liftover_37_to_38(gnomad_variant_id)


@app.get("/evidence/gnomad/clinvar_variation_id",
         summary="Given gnomad variant id, return clinvar variation id",
         response_description=RESPONSE_DESCRIPTION,
         description="Return clinvar variation id",
         response_model=Response)
def get_clinvar_variation_id(
    gnomad_variant_id: str = Query(..., description="gnomAD variant ID"),
    reference_genome: Optional[str] = Query(None, description="GRCh37 or GRCh38")
) -> Response:
    """Given gnomad variant id, return clinvar variant id

    :param str gnomad_variant_id: gnomad variant ID
    :param Optional[str] reference_genome: GRCh37 or GRCh38
    :return: Clinvar variant ID
    """
    return gnomad.clinvar_variation_id(unquote(gnomad_variant_id), reference_genome)


@app.get("/evidence/gnomad/frequency_data",
         summary="Given variant id, return gnomAD Frequency.",
         response_description=RESPONSE_DESCRIPTION,
         description="Return gnomAD population frequency data for variant.",
         response_model=Response)
def get_gnomad_frequency(
    variant_id: str = Query(
        None,
        description="gnomAD variant ID, rsID, Clin Gen Allele Registry ID, or ClinVar variation ID"),  # noqa: E501
    reference_genome: Optional[ReferenceGenome] = Query(
        None, description="Reference genome for `variant_id`. Must be either `GRCh38` or `GRCh37`")  # noqa: E501
) -> Response:
    """Return gnomAD population frequency data for variant.

    :param str variant_id: variation id
    :return: Return gnomAD population frequency data for variant.
    """
    return gnomad.frequency_data(unquote(variant_id), reference_genome)
