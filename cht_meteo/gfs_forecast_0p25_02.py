# -*- coding: utf-8 -*-
"""
Created on Thu May 20 10:32:33 2021

@author: ormondt
"""

import os
from siphon.catalog import TDSCatalog
from xarray.backends import NetCDF4DataStore
import xarray as xr
import pandas as pd
import numpy as np
from pyproj import CRS
from metpy.units import units

class Dataset():

    def __init__(self):

        self.quantity = None
        self.unit     = None
        self.time     = None
        self.x        = None
        self.y        = None
        self.crs      = None
        self.val      = None

def download(param_list, lon_range, lat_range, time_range, cycle_time,
             path, prefix):

    cycle_string = cycle_time.strftime("%Y%m%d_%H%M")

    base_url = "https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/"
    url = base_url + "GFS_Global_0p25deg_" + cycle_string + ".grib2/catalog.xml"

    gfs   = TDSCatalog(url)
    ds    = list(gfs.datasets.values())[0]
    ncss  = ds.subset()
        
    datasets = {}

    # Loop through requested parameters     
    for param in param_list:
        
        dataset = Dataset()
        dataset.crs = CRS.from_epsg(4326)
 
        fac = 1.0

        if param == "wind_u":
            var_name = "u-component_of_wind_height_above_ground"
        elif param == "wind_v":    
            var_name = "v-component_of_wind_height_above_ground"
        elif param == "barometric_pressure":
            var_name = "Pressure_reduced_to_MSL_msl"
        elif param == "precipitation":    
            var_name = "Precipitation_rate_surface"
            fac = 3600.0
        
        dataset.quantity = param            
        query = ncss.query()            
                    
        if param == "wind_u" or param == "wind_v":
            query.lonlat_box(north=lat_range[1], \
                             south=lat_range[0], \
                             east=lon_range[1],  \
                             west=lon_range[0]).time_range(time_range[0], time_range[1]).vertical_level(10.0)
        else:
            query.lonlat_box(north=lat_range[1], \
                             south=lat_range[0], \
                             east=lon_range[1],  \
                             west=lon_range[0]).time_range(time_range[0], time_range[1])
            
        query.variables(var_name)
        data = ncss.get_data(query)
        data = xr.open_dataset(NetCDF4DataStore(data))

        dataset.x    = np.array(data['lon'])
        dataset.y    = np.array(data['lat'])

        val = data[var_name]
        time   = find_time_var(val)
        dataset.time = pd.to_datetime(time.data).to_pydatetime()
        dataset.unit = val.units
        val = val.metpy.unit_array.squeeze()
        dataset.val = np.array(val)*fac        
        
        datasets[param] = dataset

    # Write data to netcdf
    
    for it in range(len(datasets[param_list[0]].time)):

        time = datasets[param_list[0]].time[it]  
        time_string    = time.strftime("%Y%m%d_%H%M")
        file_name      = prefix + "." + time_string + ".nc"
        full_file_name = os.path.join(path, file_name)
        ds = xr.Dataset()
         
        okay = False
        
        for param in param_list:        
            
            val = datasets[param].val[it,:,:]    

            if not np.any(np.isnan(val)):                    
                okay = True                    
                da = xr.DataArray(val,
                                  coords=[("lat", datasets[param].y),
                                          ("lon", datasets[param].x)])
                ds[param] = da

        if okay:
            # Only write to file if there is any data
            ds.to_netcdf(path=full_file_name)
    
    return datasets

# Helper function for finding proper time variable
def find_time_var(var, time_basename='time'):
    for coord_name in var.coords:
        if coord_name.startswith(time_basename):
            return var.coords[coord_name]
    raise ValueError('No time variable found for ' + var.name)
