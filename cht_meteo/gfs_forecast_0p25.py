# -*- coding: utf-8 -*-
"""
Created on Thu May 20 10:32:33 2021

@author: ormondt
"""

import numpy as np
import pandas as pd
import xarray as xr
from pyproj import CRS
from siphon.catalog import TDSCatalog
from xarray.backends import NetCDF4DataStore


class Dataset:
    def __init__(self):
        self.quantity = None
        self.unit = None
        self.time = None
        self.x = None
        self.y = None
        self.crs = None
        self.val = None
        self.src = "gfs_forecast_0p25"


def download(param_list, lon_range, lat_range, time_range, cycle_time):
    cycle_string = cycle_time.strftime("%Y%m%d_%H%M")

    base_url = "https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/"
    url = base_url + "GFS_Global_0p25deg_" + cycle_string + ".grib2/catalog.xml"
    url = "https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml?dataset=grib/NCEP/GFS/Global_0p25deg/Best"

    gfs = TDSCatalog(url)
    ds = list(gfs.datasets.values())[0]
    ncss = ds.subset()

    datasets = []

    # Loop through requested parameters
    for param in param_list:
        dataset = Dataset()
        dataset.crs = CRS.from_epsg(4326)

        if param == "wind":
            dataset.quantity = param

            query = ncss.query()
            query.lonlat_box(
                north=lat_range[1],
                south=lat_range[0],
                east=lon_range[1],
                west=lon_range[0],
            ).time_range(time_range[0], time_range[1]).vertical_level(10.0)
            query.variables(
                "u-component_of_wind_height_above_ground",
                "v-component_of_wind_height_above_ground",
            )
            data = ncss.get_data(query)
            data = xr.open_dataset(NetCDF4DataStore(data))
            try:
                dataset.x = np.array(data["lon"])
                dataset.y = np.array(data["lat"])
            except Exception:
                dataset.x = np.array(data["longitude"])
                dataset.y = np.array(data["latitude"])
            #           time   = data['time']

            u = data["u-component_of_wind_height_above_ground"]
            time = find_time_var(u)
            dataset.time = pd.to_datetime(time.data).to_pydatetime()
            dataset.unit = u.units
            u = u.metpy.unit_array.squeeze()
            dataset.u = np.array(u)
            v = data["v-component_of_wind_height_above_ground"]
            v = v.metpy.unit_array.squeeze()
            dataset.v = np.array(v)

        else:
            # Other scalar variables

            fac = 1.0

            dataset.quantity = param

            if param == "barometric_pressure":
                var_name = "Pressure_reduced_to_MSL_msl"
            elif param == "precipitation":
                var_name = "Precipitation_rate_surface"
                fac = 3600.0

            query = ncss.query()
            query.lonlat_box(
                north=lat_range[1],
                south=lat_range[0],
                east=lon_range[1],
                west=lon_range[0],
            ).time_range(time_range[0], time_range[1])
            query.variables(var_name)
            data = ncss.get_data(query)
            data = xr.open_dataset(NetCDF4DataStore(data))

            try:
                dataset.x = np.array(data["lon"])
                dataset.y = np.array(data["lat"])
            except Exception:
                dataset.x = np.array(data["longitude"])
                dataset.y = np.array(data["latitude"])
            #            time   = data['time']
            #            d.time = pd.to_datetime(time.data).to_pydatetime()

            val = data[var_name]
            time = find_time_var(val)
            dataset.time = pd.to_datetime(time.data).to_pydatetime()
            dataset.unit = val.units
            val = val.metpy.unit_array.squeeze()
            dataset.val = np.array(val) * fac

        datasets.append(dataset)

    return datasets


# Helper function for finding proper time variable
def find_time_var(var, time_basename="time"):
    for coord_name in var.coords:
        if coord_name.startswith(time_basename):
            return var.coords[coord_name]
    raise ValueError("No time variable found for " + var.name)


# def write_wind_to_nc(lon, lat, time, u, v, meteo_name, meteo_path):

#     for it, t in enumerate(time):

#         time_string = t.strftime("%Y%m%d_%H%M")
#         file_name = meteo_name + ".wind." + time_string + ".f.nc"
#         full_file_name = os.path.join(meteo_path, file_name)
#         uu = u[it,:,:]
#         vv = v[it,:,:]
#         ds = xr.Dataset({"u": (("lat", "lon"), uu),
#                          "v": (("lat", "lon"), vv)},
#                         coords={
#                         "lon": lon,
#                         "lat": lat})
#         ds.to_netcdf(path=full_file_name)

# def write_p_to_nc(lon, lat, time, p, meteo_name, meteo_path):

#     for it, t in enumerate(time):

#         time_string = t.strftime("%Y%m%d_%H%M")
#         file_name = meteo_name + ".p." + time_string + ".f.nc"
#         full_file_name = os.path.join(meteo_path, file_name)
#         pp = p[it,:,:]
#         ds = xr.Dataset({"p": (("lat", "lon"), pp)},
#                         coords={
#                         "lon": lon,
#                         "lat": lat})
#         ds.to_netcdf(path=full_file_name)

# def write_pr_to_nc(lon, lat, time, pr, meteo_name, meteo_path):

#     for it, t in enumerate(time):

#         time_string = t.strftime("%Y%m%d_%H%M")
#         file_name = meteo_name + ".precip." + time_string + ".f.nc"
#         full_file_name = os.path.join(meteo_path, file_name)
#         pp = pr[it,:,:]
#         ds = xr.Dataset({"precip": (("lat", "lon"), pp)},
#                         coords={
#                         "lon": lon,
#                         "lat": lat})
#         ds.to_netcdf(path=full_file_name)
