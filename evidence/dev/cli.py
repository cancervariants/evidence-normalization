"""Dev CLI"""
import asyncclick as click

from evidence.dev.etl.cancer_hotspots import CancerHotspotsETL, \
    CancerHotspotsETLException
from evidence.dev.etl.cbioportal import CBioPortalETL, CBioPortalETLException


@click.command()
@click.option(
    "--transform_cancer_hotspots",
    is_flag=True,
    default=False,
    help="Transform Cancer Hotspots data"
)
@click.option(
    "--transform_cbioportal",
    is_flag=True,
    default=False,
    help="Transform cBioPortal data"
)
@click.option(
    "--transform_all",
    is_flag=True,
    default=False,
    help="Transforms all source data, currently cBioPortal and Cancer Hotspots"
)
async def cli(transform_cancer_hotspots: bool, transform_cbioportal: bool,
              transform_all: bool) -> None:
    """Execute CLI methods

    :param bool transform_cancer_hotspots: Determines whether or not to transform
        Cancer Hotspots data
    :param bool transform_cbioportal: Determines whether or not to transform cBioPortal
        data
    :param bool transform_all: Transforms all source data
    """
    if transform_all:
        transform_cbioportal_data()
        await transform_cancer_hotspots_data()
    else:
        if transform_cbioportal:
            transform_cbioportal_data()
        if transform_cancer_hotspots:
            await transform_cancer_hotspots_data()


def transform_cbioportal_data() -> None:
    """Transform cBioPortal data"""
    c = CBioPortalETL()
    try:
        c.transform_data()
    except CBioPortalETLException as e:
        click.echo(e)


async def transform_cancer_hotspots_data() -> None:
    """Transform Cancer Hotspots data"""
    c = CancerHotspotsETL()
    try:
        await c.add_vrs_identifier_to_data()
    except CancerHotspotsETLException as e:
        click.echo(e)


if __name__ == "__main__":
    cli(_anyio_backend="asyncio")
