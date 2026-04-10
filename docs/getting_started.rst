Getting started
===============

Installation
------------

Install from GitHub:

.. code-block:: bash

   pip install git+https://github.com/Deltares-research/cht_meteo.git

For ERA5 downloads, also install ``cdsapi`` and configure your CDS API key
in ``~/.cdsapirc``.

For development:

.. code-block:: bash

   git clone https://github.com/Deltares-research/cht_meteo.git
   cd cht_meteo
   pip install -e ".[tests]"

Quick example -- download GFS forecast
----------------------------------------

.. code-block:: python

   import datetime
   from cht_meteo import MeteoDatabase

   # Create a database pointing at a local storage folder
   db = MeteoDatabase(path="/data/meteo")

   # Register a GFS 0.25-degree forecast dataset
   db.add_dataset(
       dataset_name="gfs_forecast",
       source_name="gfs_forecast_0p25",
       lon_range=[-80.0, -60.0],
       lat_range=[20.0, 35.0],
   )

   ds = db.dataset["gfs_forecast"]

   # Download data for the next 48 hours
   now = datetime.datetime.now(datetime.timezone.utc)
   time_range = [now, now + datetime.timedelta(hours=48)]
   ds.download(time_range)

   # Merge downloaded files into memory
   ds.collect(time_range)

Quick example -- export to Delft3D
-----------------------------------

.. code-block:: python

   from cht_meteo.dataset_to_delft3d import write_to_delft3d_ascii

   write_to_delft3d_ascii(
       ds,
       file_name="meteo",
       path="/data/model_input",
       parameters=["wind", "barometric_pressure"],
   )

This produces ``meteo.amu`` (u-wind), ``meteo.amv`` (v-wind), and
``meteo.amp`` (pressure).
