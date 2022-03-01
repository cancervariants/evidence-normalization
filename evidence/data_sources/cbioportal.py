"""Module for accessing python client for cBioPortal."""
from typing import Dict, Optional, Union

import requests
from bravado.client import SwaggerClient
from bravado.exception import HTTPNotFound

from evidence.data_sources.base import DataSource
from evidence.schemas import SourceMeta, Response, Sources


class CBioPortal(DataSource):
    """cBioPortal class."""

    def __init__(self, study_id: str = "msk_impact_2017",
                 api_docs_url: str = "https://www.cbioportal.org/api/api-docs") -> None:
        """Initialize cBioPortal class

        :param str study_id: The id for the study to retrieve mutation data from
        :param str api_docs_url: The url to api docs
        """
        self.api_docs_url = api_docs_url
        self.cbioportal = SwaggerClient.from_url(
            self.api_docs_url,
            config={
                "validate_requests": False,
                "validate_responses": False,
                "validate_swagger_spec": False
            }
        )
        self.cbioportal_dir = dir(self.cbioportal)
        for a in self.cbioportal_dir:
            self.cbioportal.__setattr__(a.replace(" ", "_").lower(),
                                        self.cbioportal.__getattr__(a))
        self.study_id = study_id
        self.source_meta = self.source_meta()

    def source_meta(self) -> Optional[SourceMeta]:
        """Return source meta for cBioPortal"""
        r = requests.get(self.api_docs_url)
        if r.status_code == 200:
            resp = r.json()
            version = resp["info"]["version"]
            return SourceMeta(
                label=Sources.CBIOPORTAL,
                version=version[:version.index(".", 2)]
            )

    def cancer_types_summary(self, gene_id: int) -> Response:
        """Get cancer types with gene mutations data

        :param int gene_id: Entrez ID for gene
        :return: Cancer types summary for gene
        """
        try:
            self.cbioportal.genes.getGeneUsingGET(geneId=gene_id).result()
        except HTTPNotFound:
            return self.format_response(
                Response(data=dict(), source_meta_=self.source_meta)
            )

        mutations = self.cbioportal.mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(  # noqa: E501
            molecularProfileId=f"{self.study_id}_mutations",
            sampleListId=f"{self.study_id}_all",
            entrezGeneId=gene_id,
            projection="DETAILED"
        ).result()
        mutations_sample_ids = {m.sampleId for m in mutations}

        samples_lists = self.cbioportal.Sample_Lists.getAllSampleListsInStudyUsingGET(
            studyId=self.study_id,
            projection="DETAILED"
        ).result()

        tumor_type_totals: Dict[str, Dict[str, Union[int, float]]] = dict()
        for sample_list in samples_lists:
            if ":" in sample_list.name:
                tumor_type = sample_list.name.split(": ")[-1]
                tumor_type_totals[tumor_type] = {
                    "count": 0,
                    "total": len(sample_list.sampleIds)
                }
                for sample_id in sample_list.sampleIds:
                    if sample_id in mutations_sample_ids:
                        tumor_type_totals[tumor_type]["count"] += 1
                tumor_type_totals[tumor_type]["percent_altered"] = (tumor_type_totals[tumor_type]["count"] / tumor_type_totals[tumor_type]["total"]) * 100  # noqa: E501
        return self.format_response(
            Response(data=tumor_type_totals, source_meta_=self.source_meta)
        )
