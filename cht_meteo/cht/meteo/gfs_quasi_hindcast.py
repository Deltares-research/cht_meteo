# -*- coding: utf-8 -*-
import datetime
import numpy as np
import pandas as pd
from pyproj import CRS
from metpy import xarray
from .metget_utils import *


class Dataset():

    def __init__(self):

        self.quantity = None
        self.unit     = None
        self.time     = []
        self.x        = None
        self.y        = None
        self.crs      = None
        self.val      = None
        self.source   = []
        self.metget   = None


def download(param_list, lon_range, lat_range, time_range, path, prefix=None, resolution=0.25, dt=3, config_path=None):
    """Function to download coamps-tc re-forecasts using the metget tool api"""

    # Create an instance with the MetGet api
    metget = MetGet(config_path=config_path)

    # Get info
    df = metget.get_coamps_quasi_hindcast(time_range, dt=dt)
    requested_times = df.index
    ntime = len(requested_times)
    
    dss = {}  # prepare dictionary to save the datasets
    metadata = {}
    
    # Get variable names in metget
    variables = list(np.unique([var[1].metget_name for var in metget.meteo_variables.items()])) # ['wind_pressure', 'rain']
    
    # Prepare metget domain
    domain = ["gfs"] + [resolution] + [lon_range[0], lat_range[0], lon_range[1], lat_range[1]]
    for var in variables:
        ds, meta  = metget.get_meteo(time_range=time_range, dt=dt, domain=domain, var=var, multiple_forecasts=True)
        dss[var] = ds
        metadata[var] = meta

    # Do checks
    for key in dss.keys():
        if dss[key] is not None:
            data0 = dss[key]
            lon = np.array(data0['lon'])
            lat = np.array(data0['lat'])
            if lat[1] - lat[0] > 0:  # lat should be in descending order
                lat = lat[::-1]
                reverse = True
            else:
                reverse = False
            lon = (lon + 360) % 360  # added this to be consistent with gfs
            nrows = len(lat)
            ncols = len(lon)

    # Prepare datasets to save
    datasets = []
    for param in param_list:
        dataset = Dataset()
        dataset.crs = CRS.from_epsg(4326)
        dataset.quantity = param
        dataset.x = lon
        dataset.y = lat
        dataset.time = pd.to_datetime([tt.replace(tzinfo=None) for tt in requested_times]).to_pydatetime()
        if dataset.quantity == "wind":
            dataset.u    = np.empty((ntime, nrows, ncols))
            dataset.u[:] = np.NaN
            dataset.v    = np.empty((ntime, nrows, ncols))
            dataset.v[:] = np.NaN
        else:
            dataset.val    = np.empty((ntime, nrows, ncols))
            dataset.val[:] = np.NaN
        datasets.append(dataset)

    # Loop through requested parameters
    for ind, param in enumerate(param_list):
        dataset = datasets[ind]
        #  First check if there are data available for this parameter
        if dss[metget.meteo_variables[param].metget_name] is None:
            raise KeyError
        
        for it, time_i in enumerate(requested_times):
            source = metadata[metget.meteo_variables[param].metget_name]["input_files"]["gfs"][it]
            # Check naming of files to get information on cycle used
            tau = source.split(".")[-1].split("f")[-1]
            # get cycle info
            cycle_used = time_i - datetime.timedelta(hours=int(tau))
            dataset.source.append(f'gfs_{cycle_used.strftime("%Y%m%d%H")}z_{tau}')    

        # Get metadata of input used
        dataset.config = metadata[metget.meteo_variables[param].metget_name]["input"]
        
        if param == "wind":
            data = dss[metget.meteo_variables[param].metget_name]
            u = data['wind_u']
            v = data['wind_v']
            dataset.unit = u.units
            
            if np.any(np.isnan(u)) or np.any(np.isnan(v)):  # check if there are nans and fill them up
                u = u.fillna(metget.meteo_variables[param].fill_value)
                v = v.fillna(metget.meteo_variables[param].fill_value)
            
            u = u.metpy.unit_array.squeeze()
            v = v.metpy.unit_array.squeeze()

            if reverse:
                dataset.u = np.array(u[:, ::-1, :])
                dataset.v = np.array(v[:, ::-1, :])
            else:
                dataset.u = np.array(u)
                dataset.v = np.array(v)
        else:
            if param == "barometric_pressure":
                var_name = "mslp"
            elif param == "precipitation":
                var_name = 'precipitation'
            data = dss[metget.meteo_variables[param].metget_name]
            val          = data[var_name]

            # Added this check to ensure that pressure is in Pa
            if param == "barometric_pressure":
                if val.units == 'mb':
                    val = val * 100
                    val.attrs['units'] = 'Pa'
            dataset.unit = val.units

            if param == "precipitation":  # Added this check to ensure that nan in precipitation are correctly interpreted
                val = val.where(val>=0, np.nan)

            if np.any(np.isnan(val)):  # check if there are nans and fill them up
                val = val.fillna(metget.meteo_variables[param].fill_value)
                
            val = np.array(val.metpy.unit_array.squeeze())

            if reverse:
                dataset.val = np.array(val[:, ::-1, :])
            else:
                dataset.val = np.array(val)

                        
    save_to_nc(prefix, path, datasets)


    return datasets


