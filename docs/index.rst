cht_meteo
#########

``cht_meteo`` downloads, caches, and serves gridded meteorological data
(wind, pressure, precipitation) from operational and reanalysis sources
for use with coastal and hydrodynamic models (SFINCS, Delft3D, HurryWave).

Key features:

* Download from GFS, ECMWF, ERA5, COAMPS-TC, GFS-NCAR archive, and Matroos
* Automatic caching as per-timestep NetCDF files
* Export to Delft3D ASCII and NetCDF meteo format
* Export wind fields to JSON for web visualisation
* Registry-based :class:`~cht_meteo.database.MeteoDatabase` for managing
  multiple named datasets

.. toctree::
   :maxdepth: 2
   :caption: Contents

   getting_started
   user_guide/index
   api/index
   changelog
