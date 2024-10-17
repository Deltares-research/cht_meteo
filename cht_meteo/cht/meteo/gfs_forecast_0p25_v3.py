# -*- coding: utf-8 -*-
"""
Created on Thu May 20 10:32:33 2021

@author: ormondt
"""
import os
import datetime
from siphon.catalog import TDSCatalog
from xarray.backends import NetCDF4DataStore
import xarray as xr
import pandas as pd
import numpy as np

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
            time_range = [cycle_time, cycle_time + datetime.timedelta(hours=self.source_forecast_duration)]

        cycle_string = cycle_time.strftime("%Y%m%d_%H%M")
        cycle_name   = cycle_time.strftime("%Y%m%d_%Hz")

        # Make folder for the forecast
        forecast_path = os.path.join(self.path, cycle_name) # Folder to store the forecast netcdf files

        # Make folder for the forecast
        os.makedirs(forecast_path, exist_ok=True) 

        base_url = "https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/"
        url = base_url + "GFS_Global_0p25deg_" + cycle_string + ".grib2/catalog.xml"
        url = "https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml?dataset=grib/NCEP/GFS/Global_0p25deg/Best"
        # TODO right now the url is replaced to use the best (where we don't know the cycle used)
        # TODO what about using S3 bucket (e.g. https://noaa-gfs-bdp-pds.s3.amazonaws.com/index.html#gfs.20240627/)
        
        # We assume that the best uses the latest
        latest_xml  = "https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/latest.xml"
        latest_file = TDSCatalog(latest_xml).catalog_name
        latest_date = "".join(latest_file.split(".")[0].split("_")[-2:])
        latest_cycle_time = datetime.datetime.strptime(latest_date, "%Y%m%d%H%M").replace(tzinfo=datetime.timezone.utc)
        # If we request times later than the latest available cycle the latest has been used
        if time_range[0] > latest_cycle_time:
            cycle_time = latest_cycle_time
        
        gfs   = TDSCatalog(url)
        ncss    = list(gfs.datasets.values())[0].subset()

        param_list = ["wind", "barometric_pressure", "precipitation"]

        ds = xr.Dataset()

        # Loop through requested parameters     
        for param in param_list:
    
            if param == "wind":
                
                # dataset.quantity = param
                
                query = ncss.query()
                query.lonlat_box(north=self.lat_range[1], \
                                south=self.lat_range[0], \
                                east=self.lon_range[1],  \
                                west=self.lon_range[0]).time_range(time_range[0], time_range[1]).vertical_level(10.0)
                query.variables("u-component_of_wind_height_above_ground",
                                "v-component_of_wind_height_above_ground")
                ncss_data = ncss.get_data(query)
                ds0 = xr.open_dataset(NetCDF4DataStore(ncss_data))

                lat = np.array(ds0["latitude"])                 
                lat = np.flip(lat)

                ds["lon"] = np.array(ds0["longitude"]) - 360.0
                ds["lat"] = lat
                ds["time"] = ds0["time"]

                ds["wind_u"] = xr.DataArray(np.flip(np.squeeze(ds0['u-component_of_wind_height_above_ground'].values)), dims=("time", "lat", "lon"))
                ds["wind_v"] = xr.DataArray(np.flip(np.squeeze(ds0['v-component_of_wind_height_above_ground'].values)), dims=("time", "lat", "lon"))

                ds0.close()
            
            else:
                
                # Other scalar variables                
                fac = 1.0

                if param == "barometric_pressure":
                    var_name = "Pressure_reduced_to_MSL_msl"
                elif param == "precipitation":    
                    var_name = "Precipitation_rate_surface"
                    fac = 3600.0
            
                query = ncss.query()            
                query.lonlat_box(north=self.lat_range[1], \
                                south=self.lat_range[0], \
                                east=self.lon_range[1],  \
                                west=self.lon_range[0]).time_range(time_range[0], time_range[1])
                query.variables(var_name)
                ncss_data = ncss.get_data(query)
                ds0 = xr.open_dataset(NetCDF4DataStore(ncss_data))

                # Check if lon, lat and time are already in the dataset
                if "lon" not in ds or "lat" not in ds or "time" not in ds:
                    ds["lon"] = np.array(ds0["longitude"]) - 360.0
                    lat = np.array(ds0["latitude"])                 
                    lat = np.flip(lat)
                    ds["lat"] = lat
                    ds["time"] = ds0["time"]

                v = np.flip(np.squeeze(ds0[var_name].values))

                if param == "precipitation":
                    v = v * fac

                ds[param] = xr.DataArray(v, dims=("time", "lat", "lon"))

                ds0.close()

        write2nc(ds, self.name, os.path.join(self.path, cycle_name))

        ds.close()

        self.ds = ds

def write2nc(ds, meteo_name, meteo_path):

    # Loop though times in ds
    times = ds["time"].values
    for it, t in enumerate(times):        
        time_string = pd.to_datetime(t).strftime("%Y%m%d_%H%M")    
        file_name = meteo_name + "." + time_string + ".nc"
        full_file_name = os.path.join(meteo_path, file_name)        
        ds_time = ds.isel(time=it)
        # Remove time and reftime
        ds_time = ds_time.drop_vars(["time", "reftime"])
        ds_time.to_netcdf(path=full_file_name)
        ds_time.close()
