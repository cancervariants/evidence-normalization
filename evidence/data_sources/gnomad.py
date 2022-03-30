"""Module for retrieving gnomAD Frequency data."""
import json
from typing import Dict, Optional

import requests

from evidence.data_sources.base import DataSource
from evidence.schemas import SourceMeta, Response, GnomadDataset, \
    ReferenceGenome, Sources


class GnomAD(DataSource):
    """gnomAD class."""

    def __init__(self) -> None:
        """Initialize gnomAD class"""
        self.api = "https://gnomad.broadinstitute.org/api"
        self.datasets = [v.value for v in GnomadDataset.__members__.values()]
        self.reference_genomes = [v.value for v in ReferenceGenome.__members__.values()]
        self.population_ids = {
            "afr": "African/African American",
            "ami": "Amish",
            "amr": "Latino/Admixed American",
            "asj": "Ashkenazi Jewish",
            "eas": "East Asian",
            "jpn": "Japanese",
            "kor": "Korean",
            "oea": "Other East Asian",
            "fin": "European (Finnish)",
            "mid": "Middle Eastern",
            "nfe": "European (non-Finnish)",
            "bgr": "Bulgarian",
            "est": "Estonian",
            "nwe": "North-western European",
            "onf": "Other non-Finnish European",
            "seu": "Southern European",
            "swe": "Swedish",
            "oth": "Other",
            "sas": "South Asian"
        }

    @staticmethod
    def reference_genome_to_dataset(ref_genome: ReferenceGenome) -> Optional[str]:
        """Return dataset associated to reference genome

        :param ReferenceGenome ref_genome: Reference genome
        :return: gnomad dataset
        """
        dataset = None
        if ref_genome == ReferenceGenome.GRCH38:
            dataset = GnomadDataset.GNOMAD_R3.value
        elif ref_genome == ReferenceGenome.GRCH37:
            dataset = GnomadDataset.GNOMAD_R2_1.value
        return dataset

    @staticmethod
    def dataset_to_reference_genome(dataset: GnomadDataset) -> Optional[str]:
        """Return reference genome associated to dataset

        :param GnomadDataset dataset: gnomad dataset
        :return: Reference genome
        """
        ref_genome = None
        if dataset == GnomadDataset.GNOMAD_R3:
            ref_genome = ReferenceGenome.GRCH38
        elif dataset == GnomadDataset.GNOMAD_R2_1:
            ref_genome = ReferenceGenome.GRCH37
        return ref_genome

    def query(self, query: str, variables: Dict) -> Dict:
        """Return response for API query.

        :param str query: Query
        :param Dict variables: Variables in query
        :return: Response for API query
        """
        return requests.post(
            self.api,
            data=json.dumps({
                "query": query,
                "variables": variables,
            }),
            headers={
                "Content-Type": "application/json",
            },
        ).json()

    def liftover_37_to_38(self, gnomad_variant_id: str) -> Response:
        """Liftover 37 to 38

        :param str gnomad_variant_id: gnomad variant id to liftover
        :return: gnomad 38 data
        """
        query = """
        query liftover37GnomadVariant(
            $source_variant_id: String!, $reference_genome: ReferenceGenomeId!) {
            liftover(
                source_variant_id: $source_variant_id,
                reference_genome: $reference_genome
            ) {
                liftover {variant_id, reference_genome},
                datasets
            }
        }
        """
        variables = {
            "source_variant_id": gnomad_variant_id,
            "reference_genome": self.reference_genomes[1]
        }
        resp = self.query(query, variables)
        if "errors" in resp or len(resp["data"]["liftover"]) == 0:
            return self.format_response(
                Response(source_meta_=SourceMeta(label=Sources.GNOMAD))
            )
        else:
            liftover = resp["data"]["liftover"][0]["liftover"]
            data = {
                "gnomad_variant_id": liftover["variant_id"],
                "reference_genome": liftover["reference_genome"],
                "dataset": self.datasets[0]
            }
            return self.format_response(
                Response(data=data, source_meta_=SourceMeta(label=Sources.GNOMAD))
            )

    def liftover_38_to_37(self, gnomad_variant_id: str) -> Response:
        """Liftover 38 to 37

        :param str gnomad_variant_id: gnomad variant id to liftover
        :return: gnomad 37 data
        """
        query = """
        query liftover38GnomadVariant(
            $liftover_variant_id: String!, $reference_genome: ReferenceGenomeId!) {
                liftover(
                    liftover_variant_id: $liftover_variant_id,
                    reference_genome: $reference_genome
                ) {
                    source {variant_id, reference_genome},
                    datasets
                }
        }
        """
        variables = {
            "liftover_variant_id": gnomad_variant_id,
            "reference_genome": self.reference_genomes[0]
        }
        resp = self.query(query, variables)
        if "errors" in resp or len(resp["data"]["liftover"]) == 0:
            return self.format_response(
                Response(source_meta_=SourceMeta(label=Sources.GNOMAD))
            )
        else:
            liftover = resp["data"]["liftover"][0]["source"]
            data = {
                "gnomad_variant_id": liftover["variant_id"],
                "reference_genome": liftover["reference_genome"],
                "dataset": self.datasets[1]
            }
            return self.format_response(
                Response(data=data, source_meta_=SourceMeta(label=Sources.GNOMAD))
            )

    def clinvar_variation_id(self, gnomad_variant_id: str,
                             reference_genome: ReferenceGenome = None) -> Response:
        """Return clinvar variant id given gnomad variant id

        :param str gnomad_variant_id: gnomad variant id
        :param Optional[ReferenceGenome] reference_genome: reference genome assembly
        :return: Clinvar variant id if it exists
        """
        query = """
        query getClinVarVariant($variant_id: String!,
            $reference_genome: ReferenceGenomeId!) {
            clinvar_variant(variant_id: $variant_id,
            reference_genome: $reference_genome) {
                clinvar_variation_id
            }
        }
        """
        resp = dict()
        if reference_genome:
            variables = {
                "variant_id": gnomad_variant_id,
                "reference_genome": reference_genome
            }
            resp = self.query(query, variables)
        else:
            for reference_genome in self.reference_genomes:
                variables = {
                    "variant_id": gnomad_variant_id,
                    "reference_genome": reference_genome
                }
                resp = self.query(query, variables)
                if "errors" not in resp:
                    break
        if resp and "errors" not in resp:
            return self.format_response(
                Response(data=resp["data"]["clinvar_variant"],
                         source_meta_=SourceMeta(label=Sources.GNOMAD))
            )
        else:
            return self.format_response(
                Response(source_meta_=SourceMeta(label=Sources.GNOMAD))
            )

    def variant_id_to_gnomad_id(
            self, variant_id: str,
            reference_genome: Optional[ReferenceGenome] = None) -> Dict:
        """Convert variant_id to gnomAD variation id.

        :param str variant_id: ID for variant to search
        :param Optional[ReferenceGenome] reference_genome: Reference genome
        :return: Dictionary containing variation ID and dataset used
        """
        query = """
        query getVariantSearch($query: String!, $dataset: DatasetId!) {
            variant_search(query: $query, dataset: $dataset) {
                variant_id
            }
        }
        """
        response = dict()
        if reference_genome:
            dataset = self.reference_genome_to_dataset(reference_genome)
            variables = {
                "query": variant_id,
                "dataset": dataset
            }
            response = self.query(query, variables)
            if "errors" not in response:
                if response.get("data").get("variant_search"):
                    response["data"]["dataset"] = dataset
                    return response
        else:
            for dataset in self.datasets:
                variables = {
                    "query": variant_id,
                    "dataset": dataset
                }
                response = self.query(query, variables)
                if "errors" not in response:
                    if response.get("data").get("variant_search"):
                        response["data"]["dataset"] = dataset
                        return response
        return response

    def variant(self, gnomad_variant_id: str,
                reference_genome: Optional[ReferenceGenome] = None) -> Dict:
        """Get gnomAD data for variant.

        :param str gnomad_variant_id: gnomAD variant id
        :param Optional[ReferenceGenome] reference_genome: Reference genome
        :return: gnomAD data for variant
        """
        query = """
        query getVariant($variantId: String!, $dataset: DatasetId!) {
            variant(variantId: $variantId, dataset: $dataset) {
                variantId
                reference_genome
                genome {
                    ac
                    an
                    populations {
                        id
                        ac
                        an
                    }
                }
                exome {
                    ac
                    an
                    populations {
                        id
                        ac
                        an
                    }
                }
            }
        }
        """
        if reference_genome:
            datasets = [self.reference_genome_to_dataset(reference_genome)]
        else:
            datasets = self.datasets
        for dataset in datasets:
            variables = {
                "variantId": gnomad_variant_id,
                "dataset": dataset
            }
            response = self.query(query, variables)
            if "errors" not in response:
                return response
        return {"errors": {"message": f"Variant not found: {gnomad_variant_id}"}}

    def population_data(self, data: Dict, result: Dict = None) -> Dict:
        """Merge population data into result.

        :param Dict data: Population data for exome or genome
        :param Optional[Dict] result: Merged result object
        :return: Merged result object containing population data
        """
        if not result:
            result = {
                "ac": 0,
                "an": 0,
                "populations": list()
            }
            for k in self.population_ids:
                result["populations"].append({
                    "id": k,
                    "ac": 0,
                    "an": 0
                })

        if not data:
            return result

        for k, v in data.items():
            if k in ["ac", "an"]:
                result[k] += v
            elif k == "populations":
                for pop_data in v:
                    for result_data in result["populations"]:
                        if pop_data["id"] == result_data["id"]:
                            for field in ["ac", "an"]:
                                result_data[field] += pop_data[field]
        return result

    def frequency_data(self, variant_id: str,
                       reference_genome: Optional[ReferenceGenome] = None) -> Response:
        """Get gnomad frequency data for variant

        :param str variant_id: gnomAD variant id, clingen allele registry id,
            rsid, or gnomad variant id
        :param Optional[ReferenceGenome] reference_genome: Reference genome
        :return: Population data for variant
        """
        data = {
            "variant": None,
            "assembly": None,
            "total_observations": {
                "allele_count": None,
                "allele_number": None,
                "decimal": None
            },
            "max_pop_freq": {
                "population": None,
                "allele_count": None,
                "allele_number": None,
                "decimal": None
            },
            "gnomad_url": None
        }
        variant_id_resp = self.variant_id_to_gnomad_id(variant_id, reference_genome)
        if "errors" in variant_id_resp or len(variant_id_resp["data"]["variant_search"]) == 0:  # noqa: E501
            return self.format_response(
                Response(source_meta_=SourceMeta(label=Sources.GNOMAD))
            )

        variant_id = variant_id_resp["data"]["variant_search"][0]["variant_id"]
        if not reference_genome:
            reference_genome = self.dataset_to_reference_genome(
                variant_id_resp["data"]["dataset"])

        variant_resp = self.variant(variant_id, reference_genome)
        if "errors" in variant_resp:
            return self.format_response(
                Response(source_meta_=SourceMeta(label=Sources.GNOMAD))
            )
        gnomad_frequency_resp = variant_resp["data"]["variant"]
        data["variant"] = gnomad_frequency_resp["variantId"]
        data["assembly"] = gnomad_frequency_resp["reference_genome"]
        dataset = self.reference_genome_to_dataset(reference_genome)

        data["gnomad_url"] = f"https://gnomad.broadinstitute.org/variant/{data['variant']}?dataset={dataset}"  # noqa: E501

        result = dict()
        for k in ["genome", "exome"]:
            result = self.population_data(gnomad_frequency_resp[k], result)
        data["total_observations"]["allele_count"] = result["ac"]
        data["total_observations"]["allele_number"] = result["an"]
        data["total_observations"]["decimal"] = f"{(result['ac']/result['an']):.9f}"

        max_pop_freq = {"id": None, "ac": 0, "an": 0, "decimal": 0}
        for pop_data in result["populations"]:
            try:
                pop_data["decimal"] = pop_data["ac"] / pop_data["an"]
            except ZeroDivisionError:
                continue
            else:
                if pop_data["decimal"] > max_pop_freq["decimal"]:
                    max_pop_freq = pop_data

        for response_label, max_pop_freq_label in [("population", "id"),
                                                   ("allele_count", "ac"),
                                                   ("allele_number", "an"),
                                                   ("decimal", "decimal")]:
            if response_label == "population":
                max_pop_freq[max_pop_freq_label] = self.population_ids[max_pop_freq[max_pop_freq_label]]  # noqa: E501
            elif response_label == "decimal":
                max_pop_freq[response_label] = f"{(max_pop_freq[response_label]):.9f}"
            data["max_pop_freq"][response_label] = max_pop_freq[max_pop_freq_label]

        return self.format_response(
            Response(data=data,
                     source_meta_=SourceMeta(label=Sources.GNOMAD, version=dataset))
        )
