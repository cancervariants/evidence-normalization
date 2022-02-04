"""Module containing schemas"""
from enum import Enum
from typing import Dict, Optional, Union, List

from pydantic import BaseModel
from pydantic.types import StrictStr


class Sources(str, Enum):
    """Define data sources"""

    GNOMAD = "gnomAD"
    CANCER_HOTSPOTS = "Cancer Hotspots"
    CBIOPORTAL = "cBioPortal"


class SourceMeta(BaseModel):
    """Metadata for sources"""

    label: StrictStr
    version: StrictStr


class Response(BaseModel):
    """Response model"""

    data: Optional[Union[Dict, List[Dict]]]
    source_meta_: Optional[SourceMeta]


class GnomadDataset(str, Enum):
    """Define datasets used in gnomad"""

    GNOMAD_R3 = "gnomad_r3"
    GNOMAD_R2_1 = "gnomad_r2_1"


class ReferenceGenome(str, Enum):
    """Define reference genome assemblies"""

    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"
