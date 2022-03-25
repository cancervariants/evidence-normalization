"""Dev CLI"""
import click

from evidence.dev.etl.cancer_hotspots import CancerHotspotsETL, \
    CancerHotspotsETLException


@click.command()
@click.option(
    "--normalize_cancer_hotspots",
    is_flag=True,
    default=False,
    help="Normalize Cancer Hotspots data"
)
def cli(normalize_cancer_hotspots: bool) -> None:
    """Execute CLI methods

    :param bool normalize_cancer_hotspots: Determines whether or not to normalize
        Cancer Hotspots data
    """
    if normalize_cancer_hotspots:
        c = CancerHotspotsETL()
        try:
            c.add_vrs_identifier_to_data()
        except CancerHotspotsETLException as e:
            click.echo(e)


if __name__ == "__main__":
    cli()
