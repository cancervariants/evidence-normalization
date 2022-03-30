"""Module for the base data source class"""
import json
import hashlib

from evidence.schemas import Response


class DataSource:
    """A base class for data sources"""

    @staticmethod
    def format_response(resp: Response) -> Response:
        """Add `id` to resp object if data exists

        :param Response resp: Response object
        :return: Response object with `id` field added if data exists
        """
        if resp.data:
            blob = json.dumps(
                resp.dict(), sort_keys=True, separators=(",", ":"), indent=None
            ).encode("utf-8")
            digest = hashlib.md5(blob)
            resp.id = f"normalize.evidence:{digest.hexdigest()}"
        return resp
