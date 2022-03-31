"""The VICC library for normalizing evidence"""
from pathlib import Path
from os import environ
import logging

APP_ROOT = Path(__file__).resolve().parents[0]
DATA_DIR_PATH = environ.get("DATA_DIR_PATH", APP_ROOT / "data")
SEQREPO_DATA_PATH = environ.get("SEQREPO_DATA_PATH", "/usr/local/share/seqrepo/latest")

if environ.get("EVIDENCE_PROD") == "True":
    ENV_NAME = "PROD"
    environ["VARIATION_NORM_EB_PROD"] = "True"
    environ["GENE_NORM_EB_PROD"] = "True"
else:
    ENV_NAME = "DEV"

logging.basicConfig(
    filename="evidence-normalizer.log",
    format="[%(asctime)s] - %(name)s - %(levelname)s : %(message)s"
)
logger = logging.getLogger("evidence")
logger.setLevel(logging.DEBUG)

logging.getLogger("boto3").setLevel(logging.INFO)
logging.getLogger("botocore").setLevel(logging.INFO)
logging.getLogger("python_jsonschema_objects").setLevel(logging.INFO)
logging.getLogger("swagger_spec_validator").setLevel(logging.INFO)
logging.getLogger("urllib3").setLevel(logging.INFO)
logging.getLogger("s3transfer.utils").setLevel(logging.INFO)
logging.getLogger("s3transfer.tasks").setLevel(logging.INFO)
logging.getLogger("s3transfer.futures").setLevel(logging.INFO)
logging.getLogger("hgvs.parser").setLevel(logging.INFO)
logging.getLogger("biocommons.seqrepo.seqaliasdb.seqaliasdb").setLevel(logging.INFO)
logging.getLogger("biocommons.seqrepo.fastadir.fastadir").setLevel(logging.INFO)
