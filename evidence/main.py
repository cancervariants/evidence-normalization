"""Main application for FastAPI."""
import html
from typing import Dict, Optional

from fastapi import FastAPI, HTTPException, Query
from fastapi.openapi.utils import get_openapi

from evidence.data_sources import CancerHotspots, cBioPortal, gnomAD
from evidence.version import __version__
from evidence.schemas import Response


class InvalidParameterException(Exception):
    """Exception for invalid parameter args provided by the user."""

    def __init__(self, message: str) -> None:
        """Create new instance
        :param str message: string describing the nature of the error
        """
        super().__init__(message)


gnomad = gnomAD()
cbioportal = cBioPortal()
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


@app.get("/evidence/cancer_hotspots/mutation_by_disease",
         summary="Given variant data, return cancer hotspots data.",
         response_description="A response to a validly-formed query.",
         description="Return cancer hotspots data for variant..")
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
    try:
        resp = cancer_hotspots.hotspot_data(html.unescape(so_id), vrs_variation_id)
    except InvalidParameterException as e:
        raise HTTPException(status_code=422, detail=str(e))
    return resp


@app.get("/evidence/cbioportal/cancer_types_summary",
         summary="Given entrez ID for a gene, return cancer types summary.",
         response_description="A response to a validly-formed query.",
         description="Return cancer types with gene mutations.")
def get_cancer_types_summary(
    entrez_gene_id: str = Query(..., description="Entrez ID for gene.")
) -> Response:
    """Return cancer types with gene mutations.

    :param str entrez_gene_id: Enetrez ID for gene
    :return: Return cancer types with `entrez_gene_id` mutations
    """
    try:
        resp = cbioportal.get_mutation_data(entrez_gene_id)
    except InvalidParameterException as e:
        raise HTTPException(status_code=422, detail=str(e))
    return resp


@app.get("/evidence/gnomad/liftover/38_to_37")
def gnomad_38_to_37(
    gnomad_variant_id: str = Query(..., decsription="gnomAD variant ID on 38 assembly")
) -> Response:
    """Given 38 gnomad variant id, liftover to 37

    :param str gnomad_variant_id: gnomAD variant ID on 38 assembly
    :return: gnomad 37 assembly data
    """
    return gnomad.liftover_38_to_37(gnomad_variant_id)


@app.get("/evidence/gnomad/liftover/37_to_38")
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
         response_description="A response to a validly-formed query.",
         description="Return clinvar variation id")
def get_clinvar_variation_id(
    gnomad_variant_id: str = Query(..., description="gnomAD variant ID"),
    reference_genome: Optional[str] = Query(None, description="GRCh37 or GRCh38")
) -> Response:
    """Given gnomad variant id, return clinvar variant id

    :param str gnomad_variant_id: gnomad variant ID
    :param Optional[str] reference_genome: GRCh37 or GRCh38
    :return: Clinvar variant ID
    """
    return gnomad.clinvar_variation_id(html.unescape(gnomad_variant_id),
                                       reference_genome)


@app.get("/evidence/gnomad/frequency_data",
         summary="Given variant id, return gnomAD Frequency.",
         response_description="A response to a validly-formed query.",
         description="Return gnomAD population frequency data for variant.")
def get_gnomad_frequency(
    variant_id: str = Query(
        None,
        description="gnomAD variant ID, rsID, Clin Gen Allele Registry ID, or ClinVar variation ID"),  # noqa: E501
) -> Response:
    """Return gnomAD population frequency data for variant.

    :param str variant_id: variation id
    :return: Return gnomAD population frequency data for variant.
    """
    if not variant_id:
        raise HTTPException(status_code=422, detail="Must provide `variant_id`")
    try:
        resp = gnomad.frequency_data(html.unescape(variant_id))
    except InvalidParameterException as e:
        raise HTTPException(status_code=422, detail=str(e))
    return resp