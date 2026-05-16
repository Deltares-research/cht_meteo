"""Example 'custom' meteo dataset for cht_meteo.

A *custom* dataset is used when meteorological data has already been
downloaded (or produced) by some external process and stored somewhere on
disk. The dataset does not download anything itself - it only converts the
already-available raw data into the cosmos netCDF format: one netCDF file per
time step, containing the variables ``wind_u``, ``wind_v``,
``barometric_pressure`` and ``precipitation`` on a regular ``lat`` / ``lon``
grid.

This file is meant as a template. Copy it, give the class a descriptive name,
and replace the body of :meth:`MeteoDatasetCustom.read_raw_cycle_data` with
code that reads your own data format. The rest (cycle-folder handling and
writing the per-time-step netCDF files) can usually stay as-is.
"""

import datetime
import os

import pandas as pd
import xarray as xr

from .dataset import MeteoDataset


class MeteoDatasetCustom(MeteoDataset):
    """Meteo dataset that converts already-downloaded data to cosmos format.

    Parameters
    ----------
    **kwargs
        Forwarded to :class:`MeteoDataset`. The location of the raw,
        already-downloaded data is expected in ``source_path``.
    """

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)

        # Source information
        self.source_name = "custom"
        self.source_type = "forecast"
        self.source_delay = 0  # hours after cycle time before data is available
        self.source_cycle_interval = 6  # hours between forecast cycles
        self.source_time_interval = 3  # hours between time steps within a cycle

        # Path to the raw, already-downloaded data. Set this via kwargs (or
        # after construction) to point at wherever your external process
        # stores its output.
        if not hasattr(self, "source_path"):
            self.source_path = None

    def download_forecast_cycle(self, **kwargs) -> None:
        """Convert one already-downloaded forecast cycle to cosmos netCDF format.

        Despite the name, nothing is downloaded here. The raw data for the
        requested cycle is read from disk and written out as one netCDF file
        per time step in ``<self.path>/<cycle_name>/``.

        Parameters
        ----------
        **kwargs
            Required keys:

            ``cycle_time`` : datetime
                Forecast initialisation time.

            Optional keys:

            ``time_range`` : list of datetime
                ``[start, end]`` window for which data is needed.  Defaults to
                the full forecast duration defined on the dataset.
        """

        if "cycle_time" in kwargs:
            cycle_time = kwargs["cycle_time"]
        else:
            # Throw error if cycle_time is not provided
            print("Error: cycle_time not provided")
            return

        if "time_range" in kwargs:
            time_range = kwargs["time_range"]
        else:
            # Get all data from this cycle
            time_range = [
                cycle_time,
                cycle_time + datetime.timedelta(hours=self.source_forecast_duration),
            ]

        cycle_name = cycle_time.strftime("%Y%m%d_%Hz")

        # Folder where the converted cosmos netcdf files are stored
        forecast_path = os.path.join(self.path, cycle_name)
        os.makedirs(forecast_path, exist_ok=True)

        # Read the raw, already-downloaded data and convert it to a single
        # xarray Dataset in cosmos format (see read_raw_cycle_data).
        ds = self.read_raw_cycle_data(cycle_time, time_range)
        if ds is None:
            print(f"No data found for cycle {cycle_name}")
            return

        # Write one netcdf file per time step
        write2nc(ds, self.name, forecast_path)

        ds.close()

        self.ds = ds

    def read_raw_cycle_data(self, cycle_time, time_range) -> xr.Dataset | None:
        """Read the already-downloaded raw data for one forecast cycle.

        **This is the method to customise.** Replace the body with code that
        reads your own data format from ``self.source_path`` (or wherever your
        external process stores its output).

        It must return an :class:`xarray.Dataset` with:

        * dimensions ``("time", "lat", "lon")``
        * 1-D coordinates ``time``, ``lat`` (south-to-north) and ``lon``
        * data variables ``wind_u`` and ``wind_v`` (m/s),
          ``barometric_pressure`` (Pa) and ``precipitation`` (mm/h)

        Return ``None`` if no data is available for the requested cycle.

        Parameters
        ----------
        cycle_time : datetime
            Forecast initialisation time.
        time_range : list of datetime
            ``[start, end]`` window for which data is needed.
        """
        raise NotImplementedError(
            "MeteoDatasetCustom.read_raw_cycle_data must be implemented to read "
            "your own data format. See the docstring for the expected return "
            "format."
        )


def write2nc(ds: xr.Dataset, meteo_name: str, meteo_path: str) -> None:
    """Write one netCDF file per time step in the dataset.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing a ``time`` dimension.
    meteo_name : str
        Dataset name used as the file-name prefix.
    meteo_path : str
        Output directory.
    """
    # Loop though times in ds
    times = ds["time"].to_numpy()
    for it, t in enumerate(times):
        time_string = pd.to_datetime(t).strftime("%Y%m%d_%H%M")
        file_name = f"{meteo_name}.{time_string}.nc"
        full_file_name = os.path.join(meteo_path, file_name)
        ds_time = ds.isel(time=it)
        # Remove time and reftime
        ds_time = ds_time.drop_vars(["time", "reftime"], errors="ignore")
        ds_time.to_netcdf(path=full_file_name)
        ds_time.close()
