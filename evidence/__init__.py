"""The VICC library for normalizing evidence"""
from pathlib import Path
from os import environ
import logging

APP_ROOT = Path(__file__).resolve().parents[0]
DATA_ROOT = APP_ROOT / "data"
SEQREPO_DATA_PATH = environ.get("SEQREPO_DATA_PATH", "/usr/local/share/seqrepo/latest")

logging.basicConfig(
    filename="evidence-normalizer.log",
    format="[%(asctime)s] - %(name)s - %(levelname)s : %(message)s"
)
logger = logging.getLogger("evidence")
logger.setLevel(logging.DEBUG)

logging.getLogger("bravado").setLevel(logging.INFO)
logging.getLogger("bravado_core").setLevel(logging.INFO)
logging.getLogger("boto3").setLevel(logging.INFO)
logging.getLogger("botocore").setLevel(logging.INFO)
logging.getLogger("python_jsonschema_objects").setLevel(logging.INFO)
logging.getLogger("swagger_spec_validator").setLevel(logging.INFO)
logging.getLogger("urllib3").setLevel(logging.INFO)
