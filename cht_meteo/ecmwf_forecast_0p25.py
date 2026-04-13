"""ECMWF open-data IFS forecast dataset at 0.25-degree resolution."""

import os
import shutil

import pandas as pd
import xarray as xr
from ecmwf.opendata import Client

from .dataset import MeteoDataset


class MeteoDatasetECMWFForecast0p25(MeteoDataset):
    """ECMWF IFS open-data forecast at 0.25-degree resolution.

    Downloads wind, mean-sea-level pressure and precipitation from the ECMWF
    open-data service using the ``ecmwf-opendata`` Python client.

    Parameters
    ----------
    **kwargs
        Forwarded to :class:`MeteoDataset`.
    """

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)

        # Set some source information
        self.source_name = "ecmwf_forecast_0p25"
        self.source_type = "forecast"
        self.source_delay = 8  # hours after cycle time when data becomes available
        self.source_cycle_interval = 6
        self.source_time_interval = 3

    def download_forecast_cycle(self, **kwargs) -> None:
        """Download a single ECMWF IFS forecast cycle.

        Parameters
        ----------
        **kwargs
            Required keys:

            ``cycle_time`` : datetime
                Forecast initialisation time.
        """

        if "cycle_time" in kwargs:
            cycle_time = kwargs["cycle_time"]
        else:
            # Throw error if cycle_time is not provided
            print("Error: cycle_time not provided")
            return

        cycle_name = cycle_time.strftime("%Y%m%d_%Hz")

        # Make folder for the forecast
        forecast_path = os.path.join(
            self.path, cycle_name
        )  # Folder to store the forecast netcdf files

        # Check if path already exists and holds all the required files
        if os.path.exists(forecast_path):
            existing_files = [f for f in os.listdir(forecast_path) if f.endswith(".nc")]
            expected_file_count = (72 // self.source_time_interval) + 1  # +1 for 0h
            if len(existing_files) >= expected_file_count:
                print(
                    f"Forecast cycle {cycle_name} already exists with {len(existing_files)} files. Skipping download."
                )
                return
            else:
                print(
                    f"Forecast cycle {cycle_name} exists but incomplete ({len(existing_files)}/{expected_file_count} files). Re-downloading."
                )

        # Make folder for the forecast
        os.makedirs(forecast_path, exist_ok=True)

        # =============================================================
        # SETTINGS
        # =============================================================

        # Forecast hours (0–72 every 3h)
        forecast_hours = list(range(0, 73, 3))

        # ECMWF short names for variables
        variables = ["10u", "10v", "msl", "tp"]  # include total precipitation

        # Output paths
        output_dir = os.path.join(forecast_path, "_TMP_ecmwf_data")

        lon_lim = self.lon_range  # degrees east
        lat_lim = self.lat_range  # degrees north

        os.makedirs(output_dir, exist_ok=True)

        # =============================================================
        # DOWNLOAD FORECAST FILES
        # =============================================================

        client = Client(source="ecmwf", model="ifs")

        date = cycle_time.strftime("%Y-%m-%d")
        time = cycle_time.strftime("%H:%M")

        print(f"Downloading ECMWF open forecast: {date} {time} UTC")

        for var in variables:
            print(f"Downloading {var}")
            target = os.path.join(output_dir, f"{var}.grib2")
            if not os.path.exists(target):
                try:
                    client.retrieve(
                        date=date,
                        time=cycle_time.hour,
                        step=forecast_hours,
                        param=var,
                        stream="oper",
                        type="fc",
                        target=target,
                    )
                except Exception as e:
                    print(f"Failed to get {var}: {e}")

        print("All downloads completed.")

        # =============================================================
        # LOAD, CLIP, AND MERGE DATASETS
        # =============================================================

        datasets = []

        for var in variables:
            var_datasets = []
            file = os.path.join(output_dir, f"{var}.grib2")
            if os.path.exists(file):
                ds = xr.open_dataset(file, engine="cfgrib")
                # Clip to requested region (after load)
                if lon_lim is not None and lat_lim is not None:
                    ds = ds.sel(
                        latitude=slice(lat_lim[1], lat_lim[0]),  # north to south
                        longitude=slice(lon_lim[0], lon_lim[1]),
                    )
                var_datasets.append(ds)

            if var_datasets:
                datasets.append(xr.concat(var_datasets, dim="step"))

        # Merge variables
        ds_merged = xr.merge(datasets)

        # =============================================================
        # COMPUTE RAINFALL RATE (mm/hour)
        # =============================================================

        # ECMWF total precipitation (tp) is in meters since forecast start
        # We compute 3-hour increments and convert to mm/hour
        if "tp" in ds_merged:
            tp = ds_merged["tp"]

            # Compute increments over 3-hour intervals
            tp_diff = tp.diff("step")  # tp[1]-tp[0], tp[2]-tp[1], etc.

            # Convert to mm/h (ECMWF tp is in meters since forecast start)
            tp_rate = tp_diff * 1000 / 3.0

            # Shift step coordinate backward by one interval, so 0–3h rain stored at 0h
            tp_rate = tp_rate.assign_coords(step=tp.step[:-1])
            tp_rate.name = "rain_rate"
            tp_rate.attrs["units"] = "mm/h"
            tp_rate.attrs["description"] = (
                "3-hourly mean precipitation rate (interval start time)"
            )

            # Merge back (drop cumulative field)
            ds_merged = ds_merged.drop_vars("tp")

            ds_merged = xr.merge([ds_merged, tp_rate])

            # Replace any nans in ds_merged["rain_rate"] with zeros (no rain)
            ds_merged["rain_rate"] = ds_merged["rain_rate"].fillna(0.0)

        # Change variable names to more standard ones
        # Rename latitude and longitude to lat and lon
        ds_merged = ds_merged.rename({"latitude": "lat", "longitude": "lon"})
        rename_dict = {
            "u10": "wind_u",
            "v10": "wind_v",
            "msl": "barometric_pressure",
            "rain_rate": "precipitation",
        }
        ds_merged = ds_merged.rename(rename_dict)

        # Remove "time" (this is the cycle time, not the forecast time) and rename "valid_time" to "time"
        # Also remove step, as it has an attribute "dtype"
        ds_merged = ds_merged.drop_vars(["time", "step"], errors="ignore")
        ds_merged = ds_merged.rename({"valid_time": "time"})

        # Now set the first dimension to be "time" (forecast valid time)
        ds_merged = ds_merged.swap_dims({"step": "time"})

        # And finalize flip of latitude to be south to north
        ds_merged = ds_merged.sortby("lat")

        # =============================================================
        # SAVE TO NETCDF
        # =============================================================

        write2nc(ds_merged, self.name, os.path.join(self.path, cycle_name))

        ds_merged.close()

        self.ds = ds_merged

        # Delete temporary files (grib files in output_dir)
        shutil.rmtree(output_dir, ignore_errors=True)


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
    # Loop though times in ds and save each to a separate netcdf file
    times = ds["time"].to_numpy()
    for it, t in enumerate(times):
        time_string = pd.to_datetime(t).strftime("%Y%m%d_%H%M")
        file_name = f"{meteo_name}.{time_string}.nc"
        full_file_name = os.path.join(meteo_path, file_name)
        ds_time = ds.isel(time=it)
        # Remove time and reftime
        ds_time = ds_time.drop_vars(["time"])
        ds_time.to_netcdf(path=full_file_name)
        ds_time.close()
