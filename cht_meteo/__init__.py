"""cht_meteo: Meteorological dataset management and download utilities.

Provides the MeteoDatabase and MeteoDataset classes along with a convenience
``dataset`` factory function for constructing datasets from named sources.
"""

__version__ = "0.3.1"


from cht_meteo.coamps_tc_forecast_s3 import MeteoDatasetCOAMPSTCForecastS3
from cht_meteo.database import MeteoDatabase
from cht_meteo.dataset import MeteoDataset
from cht_meteo.gfs_anl_0p50 import MeteoDatasetGFSAnalysis0p50
from cht_meteo.gfs_forecast_0p25 import MeteoDatasetGFSForecast0p25

__all__ = ["MeteoDatabase"]


def dataset(
    name: str | None = None,
    source: str | None = None,
    path: str | None = None,
    **kwargs,
) -> MeteoDataset:
    """Create a MeteoDataset for the given source.

    Parameters
    ----------
    name : str, optional
        Name to assign to the dataset.
    source : str, optional
        Source identifier, e.g. ``"gfs_forecast_0p25"``, ``"gfs_analysis_0p50"``,
        or ``"coamps_tc_forecast_s3"``.  If ``None`` a generic
        :class:`MeteoDataset` (no download support) is returned.
    path : str, optional
        Local path used to store downloaded files.
    **kwargs
        Additional keyword arguments forwarded to the dataset constructor.

    Returns
    -------
    MeteoDataset
        A concrete dataset instance for the requested source.
    """
    if source is not None:
        if source == "coamps_tc_forecast_s3":
            md = MeteoDatasetCOAMPSTCForecastS3(name=name, path=path, **kwargs)
        elif source == "gfs_forecast_0p25":
            md = MeteoDatasetGFSForecast0p25(name=name, path=path, **kwargs)
        elif source == "gfs_analysis_0p50":
            md = MeteoDatasetGFSAnalysis0p50(name=name, path=path, **kwargs)
        else:
            md = MeteoDataset(name=name)
            print(
                f"Error while reading meteo database : source {source} not recognized"
            )

    else:
        # Use generic meteo dataset (this does not have download functionality)
        md = MeteoDataset(name=name)

    return md
