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
        self.time     = []
        self.x        = None
        self.y        = None
        self.crs      = None
        self.val      = None
#        self.u        = None
#        self.v        = None

def download(param_list, lon_range, lat_range, path, prefix,
             time_range=None, times=None):

    base_url = "https://www.ncei.noaa.gov/thredds/catalog/model-gfs-g4-anl-files-old/"

    if times is not None:
        requested_times = times
        time_range = [times[0], times[-1]]
    else:    
        requested_times = pd.date_range(start=time_range[0],
                          end=time_range[1],
                          freq='3H').to_pydatetime().tolist()

    ntime = len(requested_times)

    datasets = {}

    for param in param_list:        
        dataset = Dataset()
        dataset.crs = CRS.from_epsg(4326)
        dataset.quantity = param
        datasets[param] = dataset

    icont = False
    # Get lat,lon
    for it, time in enumerate(requested_times):
        h            = requested_times[it].hour
        month_string = requested_times[it].strftime("%Y%m")
        date_string  = requested_times[it].strftime("%Y%m%d")
        url = base_url + month_string + "/" + date_string + "/catalog.xml"
        try:
            gfs   = TDSCatalog(url)
        except:
            gfs   = []
        cstr  = "0000_000"
        name = "gfsanl_4_" + date_string + "_" + cstr + ".grb2"
        okay = False               
        if gfs:
            for j, ds in enumerate(gfs.datasets):
                if ds==name:
                    okay = True
                    break
        if not okay:
            # Try the next time
            continue
        ncss  = gfs.datasets[j].subset()                
        query = ncss.query()            
        query.lonlat_box(north=lat_range[1], \
                         south=lat_range[0], \
                          east=lon_range[1],  \
                          west=lon_range[0]).vertical_level(10.0)
        query.variables("u-component_of_wind_height_above_ground")
        data = ncss.get_data(query)
        data = xr.open_dataset(NetCDF4DataStore(data))        
        lon = np.array(data['lon'])
        lat = np.array(data['lat'])  
        nrows = len(lat)
        ncols = len(lon)
        # Latitude and longitude found, so we can stop now
        icont = True
        break
    
    if not icont:
        # Could not find any data
        print("Could not find any data in requested range !")
        datasets = {}
        return datasets
    
    # initialize matrices
    for param in param_list:        
        datasets[param].x = lon
        datasets[param].y = lat
        datasets[param].val    = np.empty((ntime, nrows, ncols))
        datasets[param].val[:] = np.NaN            
    
    for it, time in enumerate(requested_times):
        
        h            = time.hour
        month_string = time.strftime("%Y%m")
        date_string  = time.strftime("%Y%m%d")
        url = base_url + month_string + "/" + date_string + "/catalog.xml"
        
        try:
            gfs   = TDSCatalog(url)
        except:
            print("Could not fetch catalogue")
            continue
        
        if h==0:
            cstr  = "0000_000"
            crstr = "0000_003"
            var_prcp = "Total_precipitation_surface_3_Hour_Accumulation"
        elif h==3:    
            cstr  = "0000_003"
            crstr = "0000_006"
            var_prcp = "Total_precipitation_surface_6_Hour_Accumulation"
        elif h==6:    
            cstr  = "0600_000"
            crstr = "0600_003"
            var_prcp = "Total_precipitation_surface_3_Hour_Accumulation"
        elif h==9:    
            cstr  = "0600_003"
            crstr = "0600_006"
            var_prcp = "Total_precipitation_surface_6_Hour_Accumulation"
        elif h==12:    
            cstr  = "1200_000"
            crstr = "1200_003"
            var_prcp = "Total_precipitation_surface_3_Hour_Accumulation"
        elif h==15:    
            cstr  = "1200_003"
            crstr = "1200_006"
            var_prcp = "Total_precipitation_surface_6_Hour_Accumulation"
        elif h==18:    
            cstr  = "1800_000"
            crstr = "1800_003"
            var_prcp = "Total_precipitation_surface_3_Hour_Accumulation"
        elif h==21:    
            cstr  = "1800_003"
            crstr = "1800_006"
            var_prcp = "Total_precipitation_surface_6_Hour_Accumulation"

        # Loop through requested parameters
        for ind, param in enumerate(param_list):
            
            datasets[param].time.append(time)

            if param == "precipitation":
                name = "gfsanl_4_" + date_string + "_" + crstr + ".grb2" 
            else:    
                name = "gfsanl_4_" + date_string + "_" + cstr + ".grb2" 

            try:
                            
                okay = False
                for j, ds in enumerate(gfs.datasets):
                    if ds==name:
                        okay = True
                        break
                
                if not okay:                    
                    # File not found, on to the next parameter
                    print(name + " was not found on server ...")
                    continue
    
                print(name + " : " + param)    
    
                ncss  = gfs.datasets[j].subset()
                
                if param == "wind_u":
                    var_name = "u-component_of_wind_height_above_ground"
                elif param == "wind_v":    
                    var_name = "v-component_of_wind_height_above_ground"
                elif param == "barometric_pressure":
                    var_name = "Pressure_reduced_to_MSL_msl"
                elif param == "precipitation":    
                    var_name = var_prcp
                            
                query = ncss.query()

                if param == "wind_u" or param == "wind_v":
                    query.lonlat_box(north=lat_range[1], \
                                     south=lat_range[0], \
                                      east=lon_range[1],  \
                                      west=lon_range[0]).vertical_level(10.0)
                else:
                    query.lonlat_box(north=lat_range[1], \
                                     south=lat_range[0], \
                                     east=lon_range[1],  \
                                     west=lon_range[0])

                query.variables(var_name)                    
                data = ncss.get_data(query)
                data = xr.open_dataset(NetCDF4DataStore(data))
                val  = data[var_name]
                datasets[param].unit = val.units
                val  = np.array(val.metpy.unit_array.squeeze())

                if param == "precipitation":
                    # Data is stored either as 3-hourly (at 03h) or 6-hourly (at 06h) accumulated rainfall
                    # For the first, just divide by 3 to get hourly precip
                    # For the second, first subtract volume that fell in the first 3 hours
                    if h==0 or h==6 or h==12 or h==18:
                        val = val/3 # Convert to mm/h
                    else:
                        val = (val - 3*np.squeeze(dataset.val[it-1,:,:]))/3

                datasets[param].val[it,:,:]  = val

            except:
                
                print("Could not download data")

        # Write data to netcdf
        time_string = time.strftime("%Y%m%d_%H%M")
        file_name = prefix + "." + time_string + ".nc"
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
