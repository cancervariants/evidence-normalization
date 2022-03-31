"""Module containing schemas"""
from enum import Enum
from typing import Dict, Optional

from pydantic import BaseModel, Field, Extra
from pydantic.types import StrictStr


class Base(BaseModel):
    """Base class for pydantic models"""

    class Config:
        """Class configs"""

        extra = Extra.forbid


class Sources(str, Enum):
    """Define data sources"""

    GNOMAD = "gnomAD"
    CANCER_HOTSPOTS = "Cancer Hotspots"
    CBIOPORTAL = "cBioPortal"


class SourceMeta(Base):
    """Metadata for sources"""

    label: Sources
    version: Optional[StrictStr] = None


class Response(Base):
    """Response model"""

    class Config(Base.Config):
        """Class configs"""

        allow_population_by_field_name = True

    id: Optional[str] = Field(alias="_id")
    data: Dict = dict()
    source_meta_: SourceMeta


class GnomadDataset(str, Enum):
    """Define datasets used in gnomad"""

    # MUST put in descending dataset
    GNOMAD_R3 = "gnomad_r3"
    GNOMAD_R2_1 = "gnomad_r2_1"


class ReferenceGenome(str, Enum):
    """Define reference genome assemblies"""

    # MUST put in descending assembly
    GRCH38 = "GRCh38"
    GRCH37 = "GRCh37"


class SourceDataType(str, Enum):
    """Define constraints for source data type for data file"""

    CANCER_HOTSPOTS_SNV = "snv"
    CANCER_HOTSPOTS_INDEL = "indel"
    CBIOPORTAL_CASE_LISTS = "case_lists"
    CBIOPORTAL_MUTATIONS = "mutations"
