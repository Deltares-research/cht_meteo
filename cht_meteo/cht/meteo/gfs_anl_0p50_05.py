# -*- coding: utf-8 -*-
"""
Created on Thu May 20 10:32:33 2021

@author: ormondt
"""

import os
# from siphon.catalog import TDSCatalog
# from xarray.backends import NetCDF4DataStore
import xarray as xr
import numpy as np
# from pyproj import CRS
# from metpy.units import units
from datetime import datetime

from .dataset import MeteoDataset

class MeteoDatasetGFSAnalysis0p50(MeteoDataset):

    # Inherit from MeteoDomain
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Set some source information
        self.source_name = "gfs_analysis_0p50"
        self.source_type = "analysis"
        self.source_delay = 6
        self.source_cycle_interval = 6
        self.source_time_interval = 3
    
    def download_analysis_times(self, requested_times, **kwargs):
        """Downloads COAMPS-TC forecast cycle for a given storm number and cycle time"""

        # ntime = len(requested_times)
        toldnew = datetime(2020, 5, 15, 0, 0, 0, 0)

        lon_range = self.lon_range
        lat_range = self.lat_range
        # GFS longitude is defined in degrees east (0 - 360)
        # If the lon_range is defined in degrees west (-180 - 180), we need to convert it
        londeg = "east"
        if lon_range[0] < 0:
            londeg = "west"
            lon_range[0] = 360.0 + lon_range[0]
            lon_range[1] = 360.0 + lon_range[1]
        # lon_range[0] = np.mod(lon_range[0], 360.0)
        # lon_range[1] = np.mod(lon_range[1], 360.0)

        icont = False

        # Get lat,lon
        for it, time in enumerate(requested_times):
            
            try:
                
                h            = requested_times[it].hour
                month_string = requested_times[it].strftime("%Y%m")
                date_string  = requested_times[it].strftime("%Y%m%d")
                
                # cstr  = "0000_000"
            
                if h==0:
                    cstr  = "0000_000"
                elif h==3:    
                    cstr  = "0000_003"
                elif h==6:    
                    cstr  = "0600_000"
                elif h==9:    
                    cstr  = "0600_003"
                elif h==12:    
                    cstr  = "1200_000"
                elif h==15:    
                    cstr  = "1200_003"
                elif h==18:    
                    cstr  = "1800_000"
                elif h==21:    
                    cstr  = "1800_003"
                
                if time<toldnew:
                    # Analysis data before May 15th, 2020 stored in different url
                    base_url = "https://www.ncei.noaa.gov/thredds/dodsC/model-gfs-g4-anl-files-old/"
                    url = base_url + month_string + "/" + date_string + "/"
                    name = "gfsanl_4_" + date_string + "_" + cstr + ".grb2"
                else:     
                    base_url = "https://www.ncei.noaa.gov/thredds/dodsC/model-gfs-g4-anl-files/"
                    url = base_url + month_string + "/" + date_string + "/"
                    name = "gfs_4_" + date_string + "_" + cstr + ".grb2"

                full_url = url + name    
                ds0 = xr.open_dataset(full_url)

                # Latitude and longitude

                # Data will be stored in ascending lat order !
                lon = ds0.lon.values[:]
                lat = ds0.lat.values[:]

                # Get lat, lon indices        

                j1 = np.where(lon<lon_range[0])[0]
                if len(j1) > 0:
                    j1 = j1[-1]
                else:
                    j1 = 0

                j2 = np.where(lon>lon_range[1])[0]
                if len(j2) > 0:
                    j2 = j2[0]
                else:
                    j2 = len(lon)

                i1 = np.where(lat>lat_range[1])[0]
                if len(i1) > 0:
                    i1 = i1[-1]
                else:
                    i1 = 0

                i2 = np.where(lat<lat_range[0])[0]
                if len(i2) > 0:
                    i2 = i2[0]
                else:
                    i2 = len(lat)

                if i2 <= i1 or j2 <= j1:
                    print("Error: cut-out is empty with given x_range and y_range")
                    return    

                # Latitude and longitude
                lat = lat[i1:i2]
                lon = lon[j1:j2]

                lat = np.flip(lat)

                if londeg == "west":
                    lon = lon - 360.0
        
                # Latitude and longitude found, so we can stop now
                icont = True
                ds0.close()

                break

            except Exception as e:
                print(e)                
                # Try another time
                print("Could not read " + full_url + " !")
        
        if not icont:
            # Could not find any data
            print("Could not find any data in requested range !")
            return
        
        r0 = np.zeros((np.size(lat), np.size(lon)))
        
        # Let's try to parallelize this
        for it, time in enumerate(requested_times):

            # Make new xarray dataset
            ds = xr.Dataset()
            ds["lon"] = xr.DataArray(lon, dims=("lon"))
            ds["lat"] = xr.DataArray(lat, dims=("lat"))

            h            = time.hour
            month_string = time.strftime("%Y%m")
            date_string  = time.strftime("%Y%m%d")
            url = base_url + month_string + "/" + date_string + "/"
            
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
            # else:
            #     print('ERROR: hours need to be round of to 3 hour values!')

            # Loop through requested parameters
            param_list = ["wind", "barometric_pressure", "precipitation"]

            for ind, param in enumerate(param_list):
                
                if param == "precipitation":
                    if time < toldnew:                
                        name = "gfsanl_4_" + date_string + "_" + crstr + ".grb2" 
                    else:
                        name = "gfs_4_" + date_string + "_" + crstr + ".grb2" 
                        
                else:    
                    if time < toldnew:                
                        name = "gfsanl_4_" + date_string + "_" + cstr + ".grb2" 
                    else:
                        name = "gfs_4_" + date_string + "_" + cstr + ".grb2" 

                try:
                                
                    okay = False
        
                    print(name + " : " + param)    
                    
                    for iattempt in range(10):    
                        try:
                            ds0 = xr.open_dataset(url + name)
                            if iattempt>0:
                                print("Success at attempt no " + int(iattempt + 1))
                            okay = True
                            break
                        except:
                            # Try again
                            pass

                    if okay:

                        if param == "wind":

                            u   = ds0["u-component_of_wind_height_above_ground"][0, 0, i1:i2, j1:j2].values
                            v   = ds0["v-component_of_wind_height_above_ground"][0, 0, i1:i2, j1:j2].values

                            u = np.flip(u, axis=0)
                            v = np.flip(v, axis=0)

                            ds["wind_u"] = xr.DataArray(u, dims=("lat", "lon"))
                            ds["wind_v"] = xr.DataArray(v, dims=("lat", "lon"))

                        else:

                            # Other scalar variables
                            if param == "barometric_pressure":
                                var_name = "Pressure_reduced_to_MSL_msl"
                            elif param == "precipitation":
                                var_name = var_prcp

                            val = ds0[var_name][0, i1:i2, j1:j2].values
                            val = np.flip(val, axis=0)

                            if param == "precipitation":
                                # Data is stored either as 3-hourly (at 03h) or 6-hourly (at 06h) accumulated rainfall
                                # For the first, just divide by 3 to get hourly precip
                                # For the second, first subtract volume that fell in the first 3 hours
                                if h==0 or h==6 or h==12 or h==18:
                                    val = val / 3 # Convert to mm/h
                                else:
                                    val = (val - 3 * np.squeeze(r0)) / 3

                                # Update r0
                                r0  = val    

                            ds[param] = xr.DataArray(val, dims=("lat", "lon"))

                    else:        
                        print("Could not get data ...")


                except:
                    
                    print("Could not download data")

            # Write to netcdf file
            ds.to_netcdf(path=os.path.join(self.path, self.name + "." + time.strftime("%Y%m%d_%H%M") + ".nc"))
            ds.close()



    # def download_time(time, lon, lat, base_url, toldnew):

    #     # Make new xarray dataset
    #     ds = xr.Dataset()
    #     ds["lon"] = xr.DataArray(lon, dims=("lon"))
    #     ds["lat"] = xr.DataArray(lat, dims=("lat"))

    #     h            = time.hour
    #     month_string = time.strftime("%Y%m")
    #     date_string  = time.strftime("%Y%m%d")
    #     url = base_url + month_string + "/" + date_string + "/"
        
    #     if h==0:
    #         cstr  = "0000_000"
    #         crstr = "0000_003"
    #         var_prcp = "Total_precipitation_surface_3_Hour_Accumulation"
    #     elif h==3:    
    #         cstr  = "0000_003"
    #         crstr = "0000_006"
    #         var_prcp = "Total_precipitation_surface_6_Hour_Accumulation"
    #     elif h==6:    
    #         cstr  = "0600_000"
    #         crstr = "0600_003"
    #         var_prcp = "Total_precipitation_surface_3_Hour_Accumulation"
    #     elif h==9:    
    #         cstr  = "0600_003"
    #         crstr = "0600_006"
    #         var_prcp = "Total_precipitation_surface_6_Hour_Accumulation"
    #     elif h==12:    
    #         cstr  = "1200_000"
    #         crstr = "1200_003"
    #         var_prcp = "Total_precipitation_surface_3_Hour_Accumulation"
    #     elif h==15:    
    #         cstr  = "1200_003"
    #         crstr = "1200_006"
    #         var_prcp = "Total_precipitation_surface_6_Hour_Accumulation"
    #     elif h==18:    
    #         cstr  = "1800_000"
    #         crstr = "1800_003"
    #         var_prcp = "Total_precipitation_surface_3_Hour_Accumulation"
    #     elif h==21:    
    #         cstr  = "1800_003"
    #         crstr = "1800_006"
    #         var_prcp = "Total_precipitation_surface_6_Hour_Accumulation"
    #     # else:
    #     #     print('ERROR: hours need to be round of to 3 hour values!')

    #     # Loop through requested parameters
    #     param_list = ["wind", "barometric_pressure", "precipitation"]

    #     for ind, param in enumerate(param_list):
            
    #         if param == "precipitation":
    #             if time < toldnew:                
    #                 name = "gfsanl_4_" + date_string + "_" + crstr + ".grb2" 
    #             else:
    #                 name = "gfs_4_" + date_string + "_" + crstr + ".grb2" 
                    
    #         else:    
    #             if time < toldnew:                
    #                 name = "gfsanl_4_" + date_string + "_" + cstr + ".grb2" 
    #             else:
    #                 name = "gfs_4_" + date_string + "_" + cstr + ".grb2" 

    #         try:
                            
    #             okay = False
    
    #             print(name + " : " + param)    
                
    #             for iattempt in range(10):    
    #                 try:
    #                     ds0 = xr.open_dataset(url + name)
    #                     if iattempt>0:
    #                         print("Success at attempt no " + int(iattempt + 1))
    #                     okay = True
    #                     break
    #                 except:
    #                     # Try again
    #                     pass

    #             if okay:

    #                 if param == "wind":

    #                     u   = ds0["u-component_of_wind_height_above_ground"][0, 0, i1:i2, j1:j2].values
    #                     v   = ds0["v-component_of_wind_height_above_ground"][0, 0, i1:i2, j1:j2].values

    #                     u = np.flip(u, axis=0)
    #                     v = np.flip(v, axis=0)

    #                     ds["wind_u"] = xr.DataArray(u, dims=("lat", "lon"))
    #                     ds["wind_v"] = xr.DataArray(v, dims=("lat", "lon"))

    #                 else:

    #                     # Other scalar variables
    #                     if param == "barometric_pressure":
    #                         var_name = "Pressure_reduced_to_MSL_msl"
    #                     elif param == "precipitation":
    #                         var_name = var_prcp

    #                     val = ds0[var_name][0, i1:i2, j1:j2].values
    #                     val = np.flip(val, axis=0)

    #                     if param == "precipitation":
    #                         # Data is stored either as 3-hourly (at 03h) or 6-hourly (at 06h) accumulated rainfall
    #                         # For the first, just divide by 3 to get hourly precip
    #                         # For the second, first subtract volume that fell in the first 3 hours
    #                         if h==0 or h==6 or h==12 or h==18:
    #                             val = val / 3 # Convert to mm/h
    #                         else:
    #                             val = (val - 3 * np.squeeze(r0)) / 3

    #                         # Update r0
    #                         r0  = val    

    #                     ds[param] = xr.DataArray(val, dims=("lat", "lon"))

    #             else:        
    #                 print("Could not get data ...")


    #         except:
                
    #             print("Could not download data")




# # Helper function for finding proper time variable
# def find_time_var(var, time_basename='time'):
#     for coord_name in var.coords:
#         if coord_name.startswith(time_basename):
#             return var.coords[coord_name]
#     raise ValueError('No time variable found for ' + var.name)

# # Helper function for finding proper time variable
# def find_height_var(var, time_basename='height'):
#     for coord_name in var.coords:
#         if coord_name.startswith(time_basename):
#             return var.coords[coord_name]
#     raise ValueError('No height variable found for ' + var.name)
