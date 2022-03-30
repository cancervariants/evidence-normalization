"""Dev CLI"""
import click

from evidence.dev.etl.cancer_hotspots import CancerHotspotsETL, \
    CancerHotspotsETLException
from evidence.dev.etl.cbioportal import CBioPortalETL, CBioPortalETLException


@click.command()
@click.option(
    "--normalize_cancer_hotspots",
    is_flag=True,
    default=False,
    help="Normalize Cancer Hotspots data"
)
@click.option(
    "--transform_cbioportal",
    is_flag=True,
    default=False,
    help="Transform cBioPortal data"
)
def cli(normalize_cancer_hotspots: bool, transform_cbioportal: bool) -> None:
    """Execute CLI methods

    :param bool normalize_cancer_hotspots: Determines whether or not to normalize
        Cancer Hotspots data
    :param bool transform_cbioportal: Determines whether or not to transform cBioPortal
        data
    """
    if transform_cbioportal:
        c = CBioPortalETL()
        try:
            c.transform_data()
        except CBioPortalETLException as e:
            click.echo(e)
    if normalize_cancer_hotspots:
        c = CancerHotspotsETL()
        try:
            c.add_vrs_identifier_to_data()
        except CancerHotspotsETLException as e:
            click.echo(e)


if __name__ == "__main__":
    cli()
