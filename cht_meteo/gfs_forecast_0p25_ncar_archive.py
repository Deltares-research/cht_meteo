# -*- coding: utf-8 -*-
"""
Created on Thu May 20 10:32:33 2021

@author: ormondt
"""

import datetime
import os

import numpy as np
import pandas as pd
import xarray as xr
from siphon.catalog import TDSCatalog
from xarray.backends import NetCDF4DataStore

from .dataset import MeteoDataset


class MeteoDatasetGFSForecast0p25NCARArchive(MeteoDataset):
    # Inherit from MeteoDomain
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Set some source information
        self.source_name = "gfs_forecast_0p25_ncar_archive"
        self.source_type = "forecast"
        self.source_delay = 6
        self.source_cycle_interval = 6
        self.source_time_interval = 3

    def download_forecast_cycle(self, **kwargs):
        """Downloads COAMPS-TC forecast cycle for a given storm number and cycle time"""

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

        lon0 = self.lon_range[0]
        lon1 = self.lon_range[1]
        if lon0 < 0.0 and lon1 < 0.0:
            lon0 = lon0 + 360.0
            lon1 = lon1 + 360.0
        lat0 = self.lat_range[0]
        lat1 = self.lat_range[1]

        cycle_year_string = cycle_time.strftime("%Y")
        cycle_day_string = cycle_time.strftime("%Y%m%d")
        cycle_string = cycle_time.strftime("%Y%m%d%H")
        cycle_name = cycle_time.strftime("%Y%m%d_%Hz")

        # Make folder for the forecast
        forecast_path = os.path.join(
            self.path, cycle_name
        )  # Folder to store the forecast netcdf files

        # Make folder for the forecast
        os.makedirs(forecast_path, exist_ok=True)

        base_url = "https://thredds.rda.ucar.edu/thredds/dodsC/files/g/d084001"
        base_url += f"/{cycle_year_string}/{cycle_day_string}"

        url = "https://thredds.rda.ucar.edu/thredds/dodsC/files/g/d084001/2024/20240927/gfs.0p25.2024092700.f000.grib2"

        # Loop through forecast times
        # Get time at 0, 3, 6, 9, 12, 15, 18, 21, 24 hours
        forecast_hours = np.arange(0, 25, 3).tolist()

        # Loop through forecast hours

        for forecast_hour in forecast_hours:

            # print(f"Downloading forecast hour {forecast_hour} of cycle {cycle_name}")

            # Create the URL for the forecast time
            url = base_url + f"/gfs.0p25.{cycle_string}.f{forecast_hour:03d}.grib2"

            ds0 = xr.open_dataset(url)

            param_list = ["wind", "barometric_pressure", "precipitation"]

            # Create a new dataset
            ds = xr.Dataset()

            # Loop through requested parameters
            for param in param_list:
                if param == "wind":

                    dau = ds0["u-component_of_wind_height_above_ground"]
                    # Clip data to lon, lat range
                    dau = dau.sel(
                        lat=slice(lat1, lat0),
                        lon=slice(lon0, lon1),
                        height_above_ground2=10.0,
                    )

                    dav = ds0["v-component_of_wind_height_above_ground"]
                    # Clip data to lon, lat range
                    dav = dav.sel(
                        lat=slice(lat1, lat0),
                        lon=slice(lon0, lon1),
                        height_above_ground2=10.0,
                    )

                    # Drop the height_above_ground2 dimension
                    dau = dau.drop("height_above_ground2")
                    dav = dav.drop("height_above_ground2")

                    lat = np.flip(np.array(dau["lat"]))
                    ds["lon"] = np.array(dau["lon"])
                    ds["lat"] = lat

                    if "time" in dau.dims:
                        tdim = "time"
                    elif "time1" in dau.dims:
                        tdim = "time1"
                    elif "time2" in dau.dims:
                        tdim = "time2"
                    elif "time3" in dau.dims:
                        tdim = "time3"

                    # Make new time dimension in ds
                    ds["time"] = xr.DataArray(
                        # flip the lat dimension
                        np.flip(np.array(dau[tdim]), axis=0),
                        dims=("time"),
                    )

                    ds["wind_u"] = xr.DataArray(
                        # flip the lat dimension
                        np.flip(dau.to_numpy(), axis=1),
                        dims=("time", "lat", "lon"),
                    )

                    ds["wind_v"] = xr.DataArray(
                        # flip the lat dimension
                        np.flip(dav.to_numpy(), axis=1),
                        dims=("time", "lat", "lon"),
                    )

                else:
                    # Other scalar variables
                    fac = 1.0

                    if param == "barometric_pressure":
                        var_name = "Pressure_reduced_to_MSL_msl"
                    elif param == "precipitation":
                        var_name = "Precipitation_rate_surface"
                        fac = 3600.0

                    da = ds0[var_name]
                    # Clip data to lon, lat range
                    da = da.sel(
                        lat=slice(lat1, lat0),
                        lon=slice(lon0, lon1),
                    )

                    ds[param] = xr.DataArray(
                        # flip the lat dimension
                        np.flip(fac * da.to_numpy(), axis=1),
                        dims=("time", "lat", "lon"),
                    )

            ds0.close()

            # Rename time dimension if the name is "time1"
            if "time1" in ds.dims:
                ds = ds.rename({"time1": "time"})

            write2nc(ds, self.name, os.path.join(self.path, cycle_name))

            ds.close()

            self.ds = ds


def write2nc(ds, meteo_name, meteo_path):
    # Loop though times in ds
    times = ds["time"].to_numpy()
    for it, t in enumerate(times):
        time_string = pd.to_datetime(t).strftime("%Y%m%d_%H%M")
        file_name = meteo_name + "." + time_string + ".nc"
        full_file_name = os.path.join(meteo_path, file_name)
        ds_time = ds.isel(time=it)
        ds_time.to_netcdf(path=full_file_name)
        ds_time.close()
