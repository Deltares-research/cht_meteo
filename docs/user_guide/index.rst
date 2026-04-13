User guide
==========

The database / dataset pattern
-------------------------------

``cht_meteo`` is built around two classes:

- :class:`~cht_meteo.database.MeteoDatabase` -- a registry holding one or
  more named datasets. Reads/writes a TOML configuration file.
- :class:`~cht_meteo.dataset.MeteoDataset` -- a container for a single gridded
  dataset. Concrete subclasses add source-specific download logic.

Typical workflow:

1. Create a ``MeteoDatabase`` and register datasets (or load from TOML).
2. Call ``dataset.download(time_range)`` to fetch and cache raw data.
3. Call ``dataset.collect(time_range)`` to assemble cached files into
   ``dataset.ds`` (an ``xr.Dataset``).
4. Pass ``dataset`` to an export function.

Supported data sources
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Source identifier
     - Description
   * - ``gfs_forecast_0p25``
     - NCEP GFS 0.25 degree global forecast via UCAR THREDDS
   * - ``gfs_analysis_0p50``
     - NCEP GFS 0.5 degree analysis via UCAR THREDDS
   * - ``gfs_forecast_0p25_ncar_archive``
     - GFS 0.25 degree historical archive via NCAR RDA
   * - ``ecmwf_forecast_0p25``
     - ECMWF open-data forecast via ``ecmwf-opendata``
   * - ``era5_reanalysis_0p25``
     - ERA5 hourly reanalysis via Copernicus CDS (requires ``cdsapi``)
   * - ``coamps_tc_forecast``
     - COAMPS-TC tropical-cyclone forecast from AWS S3
   * - ``matroos``
     - Matroos (Rijkswaterstaat) operational forecast

Variables
---------

Each dataset uses four standard variable names:

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Variable
     - Units
     - Description
   * - ``wind_u``
     - m/s
     - Eastward wind at 10 m
   * - ``wind_v``
     - m/s
     - Northward wind at 10 m
   * - ``barometric_pressure``
     - Pa
     - Mean sea-level pressure
   * - ``precipitation``
     - mm/h
     - Precipitation rate

Exporting to model input
-------------------------

**Delft3D ASCII format** (produces ``.amu``, ``.amv``, ``.amp``, ``.ampr``)::

   from cht_meteo.dataset_to_delft3d import write_to_delft3d_ascii

   write_to_delft3d_ascii(
       dataset, file_name="wind", path="/model/meteo",
       parameters=["wind", "barometric_pressure", "precipitation"],
   )

**Delft3D NetCDF format**::

   from cht_meteo.dataset_to_delft3d import write_to_delft3d_netcdf

   write_to_delft3d_netcdf(dataset, file_name="meteo.nc", path="/model/meteo")

**JSON wind format** (for web visualisation)::

   from cht_meteo.dataset_to_json_wind import write_wind_to_json

   write_wind_to_json(dataset, file_name="wind.json")

Using a TOML database file
---------------------------

For operational workflows, store dataset definitions in TOML:

.. code-block:: toml

   # meteo_database.toml
   [[meteo_dataset]]
   name   = "gfs_forecast"
   source = "gfs_forecast_0p25"
   x_range = [-80.0, -60.0]
   y_range = [20.0, 35.0]

Load with:

.. code-block:: python

   from cht_meteo import MeteoDatabase

   db = MeteoDatabase()
   db.read_datasets("meteo_database.toml")
