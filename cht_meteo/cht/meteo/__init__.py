# -*- coding: utf-8 -*-
"""
This is the old meteo toolkit that was originally in CoastalHazardsToolkit
Created on Sun Apr 25 10:58:08 2021

@author: ormondt
"""

from cht_meteo.cht.meteo.coamps_tc_forecast_s3 import MeteoDatasetCOAMPSTCForecastS3
from cht_meteo.cht.meteo.database import MeteoDatabase
from cht_meteo.cht.meteo.dataset import MeteoDataset
from cht_meteo.cht.meteo.gfs_anl_0p50 import MeteoDatasetGFSAnalysis0p50
from cht_meteo.cht.meteo.gfs_forecast_0p25 import MeteoDatasetGFSForecast0p25

__all__ = ["MeteoDatabase"]


def dataset(name=None, source=None, path=None, **kwargs):
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
