"""Main application for FastAPI."""
from typing import Dict
from urllib.parse import unquote

from fastapi import FastAPI, Query
from fastapi.openapi.utils import get_openapi

from evidence.data_sources import CancerHotspots, CBioPortal
from evidence.version import __version__
from evidence.schemas import Response

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
