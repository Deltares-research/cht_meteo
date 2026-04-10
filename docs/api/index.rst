API reference
=============

MeteoDatabase
-------------

.. autoclass:: cht_meteo.database.MeteoDatabase
   :members:
   :undoc-members:
   :show-inheritance:

MeteoDataset
------------

.. autoclass:: cht_meteo.dataset.MeteoDataset
   :members:
   :undoc-members:
   :show-inheritance:

Export functions
----------------

.. autofunction:: cht_meteo.dataset_to_delft3d.write_to_delft3d_ascii

.. autofunction:: cht_meteo.dataset_to_delft3d.write_to_delft3d_netcdf

.. autofunction:: cht_meteo.dataset_to_json_wind.write_wind_to_json

Data source classes
-------------------

.. autoclass:: cht_meteo.gfs_forecast_0p25.MeteoDatasetGFSForecast0p25
   :members:
   :show-inheritance:

.. autoclass:: cht_meteo.gfs_anl_0p50.MeteoDatasetGFSAnalysis0p50
   :members:
   :show-inheritance:

.. autoclass:: cht_meteo.ecmwf_forecast_0p25.MeteoDatasetECMWFForecast0p25
   :members:
   :show-inheritance:

.. autoclass:: cht_meteo.era5_reanalysis_0p25.MeteoDatasetERA5Reanalysis0p25
   :members:
   :show-inheritance:

.. autoclass:: cht_meteo.coamps_tc_forecast_s3.MeteoDatasetCOAMPSTCForecastS3
   :members:
   :show-inheritance:

.. autoclass:: cht_meteo.matroos_forecast.MeteoDatasetMatroos
   :members:
   :show-inheritance:
