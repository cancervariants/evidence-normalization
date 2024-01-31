"""Module containing schemas"""
from enum import Enum
from typing import Dict, Optional

from pydantic import BaseModel, Field, StrictStr, ConfigDict


class Base(BaseModel, extra="forbid"):
    """Base class for pydantic models"""


class Sources(str, Enum):
    """Define data sources"""

    CANCER_HOTSPOTS = "Cancer Hotspots"
    CBIOPORTAL = "cBioPortal"


class SourceMeta(Base):
    """Metadata for sources"""

    label: Sources
    version: Optional[StrictStr] = None


class Response(Base):
    """Response model"""

    model_config = ConfigDict(populate_by_name=True)

    id: Optional[str] = Field(default=None, alias="_id")
    data: Dict = dict()
    source_meta_: SourceMeta


class SourceDataType(str, Enum):
    """Define constraints for source data type for data file"""

    CANCER_HOTSPOTS = "mutation_hotspots"
    CBIOPORTAL_CASE_LISTS = "case_lists"
    CBIOPORTAL_MUTATIONS = "mutations"
