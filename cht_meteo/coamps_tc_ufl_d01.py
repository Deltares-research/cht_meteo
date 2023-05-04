# -*- coding: utf-8 -*-
"""
Created on Thu May 20 10:32:33 2021

@author: ormondt
"""

from siphon.catalog import TDSCatalog
from xarray.backends import NetCDF4DataStore
import xarray as xr
import pandas as pd
import numpy as np
from pyproj import CRS
from metpy.units import units
import netCDF4

class Dataset():

    def __init__(self):

        self.quantity = None
        self.unit     = None
        self.time     = None
        self.x        = None
        self.y        = None
        self.crs      = None
        self.val      = None

def download(param_list, lon_range, lat_range, time_range, cycle_time):

    cycle_string = cycle_time.strftime("%Y%m%d%H")

    base_url = "https://icoast.rc.ufl.edu/thredds/dodsC/coamps/"

    times = pd.date_range(start=time_range[0],
                          end=time_range[1],
                          freq='1H').to_pydatetime().tolist()

    datasets = []

    # Get lon and lat from first time
    nc_file = "coamps-tc_d01_" + cycle_string + "_tau000.nc"
    try:        
        nc4 = netCDF4.Dataset(base_url + nc_file)
    except:
        print("Not available: " + base_url + nc_file)
        return datasets
        
    lon0 = np.squeeze(nc4["lon"][0,:].astype('float64')) - 360.0
    lat0 = np.squeeze(nc4["lat"][:,0].astype('float64'))

    ilon0 = np.where(lon0<lon_range[0])[0][-1]
    ilon1 = np.where(lon0>lon_range[1])[0][0] + 1
    ilat0 = np.where(lat0<lat_range[0])[0][-1]
    ilat1 = np.where(lat0>lat_range[1])[0][0] + 1
    lon = lon0[ilon0:ilon1]
    lat = np.flipud(lat0[ilat0:ilat1])    
#    lat = lat0[ilat0:ilat1]    

    # Loop through requested parameters     
    for param in param_list:
        
        dataset = Dataset()
        dataset.crs = CRS.from_epsg(4326)
 
        if param == "wind":
            
            dataset.time     = times
            dataset.x        = lon
            dataset.y        = lat
            dataset.quantity = param
            dataset.unit     = "m/s"
            dataset.u        = np.zeros((np.size(times), np.size(lat), np.size(lon)))
            dataset.v        = np.zeros((np.size(times), np.size(lat), np.size(lon)))

            iok = -1            
            for it, time in enumerate(times):
                hrs = int((time.replace(tzinfo=None) - cycle_time.replace(tzinfo=None)).total_seconds()/3600)
                nc_file = "coamps-tc_d01_" + cycle_string + "_tau" + str(hrs).zfill(3) + ".nc"
                try:
                    nc4 = netCDF4.Dataset(base_url + nc_file)
                    dataset.u[it,:,:] = np.flipud(nc4["uuwind"][ilat0:ilat1, ilon0:ilon1])
                    dataset.v[it,:,:] = np.flipud(nc4["vvwind"][ilat0:ilat1, ilon0:ilon1])
                    iok = it

                except:
                    pass
            if iok<np.size(times) - 1:
                dataset.time = dataset.time[0:iok + 1]
                dataset.u    = dataset.u[0:iok + 1, :, :]
                dataset.v    = dataset.v[0:iok + 1, :, :]
        
        else:
            
            # Other scalar variables
            
            fac = 1.0
            
            dataset.quantity = param

            if param == "barometric_pressure":
                var_name = "slpres"
                dataset.unit     = "Pa"
                fac = 100.0
            elif param == "precipitation":    
                var_name = None
                dataset.unit     = "mm/h"
                fac = 3600.0

            dataset.time     = times
            dataset.x        = lon
            dataset.y        = lat
            dataset.quantity = param
            dataset.val      = np.zeros((np.size(times), np.size(lat), np.size(lon)))

            if var_name:                
                iok = -1            
                for it, time in enumerate(times):    
                    hrs = int((time.replace(tzinfo=None) - cycle_time.replace(tzinfo=None)).total_seconds()/3600)
                    nc_file = "coamps-tc_d01_" + cycle_string + "_tau" + str(hrs).zfill(3) + ".nc"
                    try:
                        nc4 = netCDF4.Dataset(base_url + nc_file)
                        dataset.val[it,:,:] = fac*np.flipud(nc4[var_name][ilat0:ilat1, ilon0:ilon1])
                        iok = it
                    except:
                        pass
                if iok<np.size(times) - 1:
                    dataset.time = dataset.time[0:iok + 1]
                    dataset.val  = dataset.val[0:iok + 1, :, :]
        if iok>=0:
            datasets.append(dataset)    
    
    return datasets
