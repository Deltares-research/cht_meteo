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


class MeteoDatasetGFSForecast0p25(MeteoDataset):
    # Inherit from MeteoDomain
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Set some source information
        self.source_name = "gfs_forecast_0p25"
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

        cycle_string = cycle_time.strftime("%Y%m%d_%H%M")
        cycle_name = cycle_time.strftime("%Y%m%d_%Hz")

        # Make folder for the forecast
        forecast_path = os.path.join(
            self.path, cycle_name
        )  # Folder to store the forecast netcdf files

        # Make folder for the forecast
        os.makedirs(forecast_path, exist_ok=True)

        base_url = (
            "https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/"
        )
        url = base_url + "GFS_Global_0p25deg_" + cycle_string + ".grib2/catalog.xml"
        url = "https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml?dataset=grib/NCEP/GFS/Global_0p25deg/Best"
        # TODO right now the url is replaced to use the best (where we don't know the cycle used)
        # TODO what about using S3 bucket (e.g. https://noaa-gfs-bdp-pds.s3.amazonaws.com/index.html#gfs.20240627/)

        # We assume that the best uses the latest
        latest_xml = "https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/latest.xml"
        latest_file = TDSCatalog(latest_xml).catalog_name
        latest_date = "".join(latest_file.split(".")[0].split("_")[-2:])
        latest_cycle_time = datetime.datetime.strptime(
            latest_date, "%Y%m%d%H%M"
        ).replace(tzinfo=datetime.timezone.utc)
        # If we request times later than the latest available cycle the latest has been used
        if time_range[0] > latest_cycle_time:
            cycle_time = latest_cycle_time

        gfs = TDSCatalog(url)
        ncss = list(gfs.datasets.values())[0].subset()

        param_list = ["wind", "barometric_pressure", "precipitation"]

        def query_ncss(var_names, vertical_level=None):
            query = ncss.query()
            query.lonlat_box(
                north=self.lat_range[1],
                south=self.lat_range[0],
                east=360,  # full longitude
                west=0
            ).time_range(time_range[0], time_range[1])
            if vertical_level is not None:
                query.vertical_level(vertical_level)
            query.variables(*var_names)
            return ncss.get_data(query)
        
        ds_full  = xr.Dataset()

        # Loop through requested parameters
        for param in param_list:
            if param == "wind":
                # dataset.quantity = param
                ncss_data = query_ncss(
                    ["u-component_of_wind_height_above_ground",
                    "v-component_of_wind_height_above_ground"],
                    vertical_level=10.0
                )

                with xr.open_dataset(NetCDF4DataStore(ncss_data)) as ds0:
                    lat = np.array(ds0["latitude"])[::-1]
                    lon = np.array(ds0["longitude"])
                    ds_full["lat"] = xr.DataArray(lat, dims="lat")
                    ds_full["lon"] = xr.DataArray(lon, dims="lon")

                    if "time1" in ds0:
                        ds0 = ds0.rename({"time1": "time"})
                    elif "time2" in ds0:
                        ds0 = ds0.rename({"time2": "time"})
                    ds_full["time"] =  xr.DataArray(np.array(ds0["time"]), dims="time")

                    ds_full["wind_u"] = xr.DataArray(
                        np.flip(
                            np.squeeze(
                                ds0[
                                    "u-component_of_wind_height_above_ground"
                                ].to_numpy()
                            ),
                        axis=1,),
                        dims=("time", "lat", "lon"),
                    )
                    ds_full["wind_v"] = xr.DataArray(
                        np.flip(
                            np.squeeze(
                                ds0[
                                    "v-component_of_wind_height_above_ground"
                                ].to_numpy()
                            ),
                        axis=1,),
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

                ncss_data = query_ncss([var_name])

                with xr.open_dataset(NetCDF4DataStore(ncss_data)) as ds0:
                    # Check if lon, lat and time are already in the dataset
                    if "lon" not in ds_full or "lat" not in ds_full or "time" not in ds_full:
                        ds_full["lon"] = xr.DataArray(np.array(ds0["longitude"]), dims="lon")
                        ds_full["lat"] = xr.DataArray(np.array(ds0["latitude"])[::-1], dims="lat")

                        if "time1" in ds0:
                            ds0 = ds0.rename({"time1": "time"})
                        elif "time2" in ds0:
                            ds0 = ds0.rename({"time2": "time"})
                        ds_full["time"] = xr.DataArray(np.array(ds0["time"]), dims="time")

                    v = np.flip(np.squeeze(ds0[var_name].to_numpy()), axis=1)

                    if param == "precipitation":
                        v = v * fac

                    ds_full[param] = xr.DataArray(v, dims=("time", "lat", "lon"))

        # --- Subset and reorder longitude based on self.lon_range ---
        lon_min, lon_max = self.lon_range
        # ds_full = ds_full.assign_coords(lon=((ds_full.lon + 360) % 360))  # ensure 0-360

        if lon_min <= lon_max:
            ds_sel = ds_full.sel(lon=slice(lon_min, lon_max))
        else:
            ds_part1 = ds_full.sel(lon=slice(lon_min, 360))
            ds_part2 = ds_full.sel(lon=slice(0, lon_max))
            ds_sel = xr.concat([ds_part1, ds_part2], dim="lon")

        write2nc(ds_sel, self.name, os.path.join(self.path, cycle_name))

        ds_sel.close()

        self.ds = ds_sel


def write2nc(ds, meteo_name, meteo_path):
    # Loop though times in ds
    times = ds["time"].to_numpy()
    for it, t in enumerate(times):
        time_string = pd.to_datetime(t).strftime("%Y%m%d_%H%M")
        file_name = meteo_name + "." + time_string + ".nc"
        full_file_name = os.path.join(meteo_path, file_name)
        ds_time = ds.isel(time=it)
        # Remove time and reftime
        ds_time = ds_time.drop_vars(["time"])
        ds_time.to_netcdf(path=full_file_name)
        ds_time.close()
