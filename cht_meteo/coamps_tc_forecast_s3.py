# -*- coding: utf-8 -*-
"""
@author:    Panos Athanasiou
Date:       04/07/2022
"""

import datetime
import json
import os
import time

import numpy as np
import pandas as pd
import requests
from pyproj import CRS

from urllib.request import urlretrieve

# from cht_meteo import gfs_forecast_0p25

class Dataset:
    def __init__(self):
        self.quantity = None
        self.unit = None
        self.time = []
        self.x = None
        self.y = None
        self.crs = None
        self.val = None
        self.src = []


def download(param_list, lon_range, lat_range, time_range, cycle_time):
    """Function to download coamps-tc forecasts using the scrubber tool functions
    Right now resolution is hardcoded to 0.25 degrees but this can change"""

    base_url = "https://coamps-tc-data.s3.us-east-2.amazonaws.com/index.html#deterministic/realtime"
    base_url = base_url + "#deterministic/realtime"

    # Get year from cycle_time
    year = cycle_time.year
    cycle_time = cycle_time.strftime("%Y%m%d%H")
    storm_number = "02L"

    tar_file = storm_number + "_" + cycle_time + "_netcdf.tar"

    url = base_url + "/" + str(year) + "/" + storm_number + "/" + cycle_time + "/" + tar_file

    # Download the tar file from url to a temporary directory
    urlretrieve(url, tar_file)


    





#     fill_values = dict(
#         wind=0.0, precipitation=0.0, barometric_pressure=102000.0
#     )  # fill values to use in case of nan
#     units = dict(
#         wind="m/s", precipitation="kg.m-2.hour-1", barometric_pressure="Pa"
#     )  # units

#     apikey = "GrPBTsclg71qgDvw4ojGf4mTwZmkbqNu6jcq7T9V"  # Deltares apikey hardcoded for the time being
#     endpoint = "https://api.metget.zachcobell.com"

#     # Connect parameters names with coamps-tc names
#     param_names = {
#         "wind": "wind_pressure",
#         "barometric_pressure": "wind_pressure",
#         "precipitation": "rain",
#     }
#     variables = list(
#         np.unique([param_names[name] for name in param_list])
#     )  # ['wind_pressure', 'rain']

#     requested_times = (
#         pd.date_range(start=time_range[0], end=time_range[1], freq="3H")
#         .to_pydatetime()
#         .tolist()
#     )

#     timestep = int(
#         (requested_times[1] - requested_times[0]).total_seconds()
#     )  # get time-step of requested time assuming that the time-step stays the same
#     ntime = len(requested_times)

#     # Get storms from endpoint metadata
#     storms = check_coamps(apikey, endpoint)

#     # Check if forecast exists for the selected dates and make a dataframe
#     fr = pd.DataFrame(index=requested_times)
#     fr["storm_id"] = None

#     for i, ti in enumerate(requested_times):
#         fr.loc[ti, "storm_id"] = []
#         for j, storm_id in enumerate(storms):
#             # The MetGet tool seems to have changed ...
#             storm = storms[storm_id]
#             ####### CHECK WHICH times to use for range ####################
#             t1 = datetime.datetime.strptime(
#                 storm["min_forecast_date"], "%Y-%m-%d %H:%M:%S"
#             )
#             t2 = datetime.datetime.strptime(
#                 storm["max_forecast_date"], "%Y-%m-%d %H:%M:%S"
#             )
#             # t1 = datetime.datetime.strptime(storm['latest_complete_forecast_start'], "%Y-%m-%d %H:%M:%S")
#             # t2 = datetime.datetime.strptime(storm['latest_complete_forecast_end'], "%Y-%m-%d %H:%M:%S")
#             belongs = (
#                 t1.replace(tzinfo=datetime.timezone.utc)
#                 <= ti.replace(tzinfo=datetime.timezone.utc)
#                 <= t2.replace(tzinfo=datetime.timezone.utc)
#             )
#             if belongs:
#                 fr.loc[ti, "storm_id"].append(storm_id)

#     # Check if there is a txt file with the priority storm
#     file_storm = r"p:\11206085-onr-fhics\03_cosmos\Priority_Storm.txt"  # !!! Hard-coded
#     if os.path.exists(file_storm):
#         with open(file_storm, "r") as f:
#             lines = f.readlines()
#         f.close()
#         pr_st_txt = lines[0]
#     else:
#         pr_st_txt = ""
#     # Get storm priority from noaa
#     pr_st_noaa = get_priority_storm()
#     method = "first"
#     for i, ti in enumerate(requested_times):
#         if not fr.loc[
#             ti, "storm_id"
#         ]:  # if there is not coamps storm make None so gfs is used
#             fr.loc[ti, "storm_id"] = None
#         else:
#             if (
#                 pr_st_txt in fr.loc[ti, "storm_id"]
#             ):  # if the priority storm is available that timestep use it
#                 fr.loc[ti, "storm_id"] = pr_st_txt
#                 continue
#             if (
#                 pr_st_noaa in fr.loc[ti, "storm_id"]
#             ):  # if the priority storm is available that timestep use it
#                 fr.loc[ti, "storm_id"] = pr_st_noaa
#                 continue
#             # else used method defined
#             if method == "first":
#                 fr.loc[ti, "storm_id"] = fr.loc[ti, "storm_id"][0]
#             elif method == "last":
#                 fr.loc[ti, "storm_id"] = fr.loc[ti, "storm_id"][-1]

#     storm_log = np.array(
#         [v is not None for v in fr["storm_id"].to_numpy()]
#     )  # check for which times there are forecasts from coamps-tc
#     if all(~storm_log):  # If no coamps-tc data are available, gfs data are downloaded
#         datasets = gfs_forecast_0p25.download(
#             param_list, lon_range, lat_range, time_range, cycle_time
#         )
#         print(
#             "Could not find any coamps-tc data in requested time range! Using gfs instead"
#         )
#         return datasets

#     # Else if there are coamps-tc data see which storms (models) to use
#     storms_id = np.unique(fr["storm_id"][storm_log])
#     dss = {}  # prepare dictionary to save the datasets

#     # Loop for the case there is more than one storm in requested times and get storm forecast

#     for i, st in enumerate(storms_id):
#         dss[st] = {}
#         #        storm = storms[[storm_i for storm_i in range(len(storms)) if storms[storm_i]['storm'] == st][0]]
#         storm = storms[st]
#         t1 = datetime.datetime.strptime(
#             storm["min_forecast_date"], "%Y-%m-%d %H:%M:%S"
#         ).replace(tzinfo=datetime.timezone.utc)
#         t2 = datetime.datetime.strptime(
#             storm["max_forecast_date"], "%Y-%m-%d %H:%M:%S"
#         ).replace(tzinfo=datetime.timezone.utc)
#         # This parts ensures that requested data are within the model range
#         t1 = max(requested_times[0], t1)
#         t2 = min(requested_times[-1], t2)
#         if t1 == t2:  # scrubber tool cannot handle single time
#             t1 = t1 - datetime.timedelta(hours=3)
#         for ii, var in enumerate(variables):
#             ds = met_get(
#                 domain=[
#                     [
#                         "coamps-{}".format(st),
#                         resolution,
#                         lon_range[0],
#                         lat_range[0],
#                         lon_range[1],
#                         lat_range[1],
#                     ]
#                 ],  # list is needed here!
#                 start=t1.strftime("%Y-%m-%d %H:%M:%S"),
#                 end=t2.strftime("%Y-%m-%d %H:%M:%S"),
#                 timestep=timestep,
#                 variable=var,
#                 apikey=apikey,
#                 endpoint=endpoint,
#                 output="coamps_{}_{}".format(
#                     st, var
#                 ),  # this is used for the online file creation
#             )
#             if ds is not None:
#                 dss[st][var] = ds[0]
#             else:
#                 dss[st][var] = None

#     # if no data could be downloaded (CHECK AGAIN WHAT SHOULD HAPPEN IN THIS CASE)
#     if all(
#         [dss[k][kk] is None for k in list(dss.keys()) for kk in list(dss[k].keys())]
#     ):
#         datasets = gfs_forecast_0p25.download(
#             param_list, lon_range, lat_range, time_range, cycle_time
#         )
#         print(
#             "Could not find any coamps-tc data in requested time range! Using gfs instead"
#         )
#         return datasets
#     else:  # Else if there are available data get lon lat info
#         for key in dss[storms_id[0]].keys():
#             if dss[storms_id[0]][key] is not None:
#                 data0 = dss[storms_id[0]][key]
#                 lon = np.array(data0["lon"])
#                 lat = np.array(data0["lat"])
#                 if lat[1] - lat[0] > 0:  # lat should be in descending order
#                     lat = lat[::-1]
#                     reverse = True
#                 else:
#                     reverse = False
#                 lon = (lon + 360) % 360  # added this to be consistent with gfs
#                 nrows = len(lat)
#                 ncols = len(lon)

#     # Prepare datasets
#     datasets = []
#     for param in param_list:
#         dataset = Dataset()
#         dataset.crs = CRS.from_epsg(4326)
#         dataset.quantity = param
#         dataset.x = lon
#         dataset.y = lat
#         dataset.time = pd.to_datetime(
#             [tt.replace(tzinfo=None) for tt in requested_times]
#         ).to_pydatetime()
#         if dataset.quantity == "wind":
#             dataset.u = np.empty((ntime, nrows, ncols))
#             dataset.u[:] = np.NaN
#             dataset.v = np.empty((ntime, nrows, ncols))
#             dataset.v[:] = np.NaN
#         else:
#             dataset.val = np.empty((ntime, nrows, ncols))
#             dataset.val[:] = np.NaN
#         datasets.append(dataset)

#     for it, time_i in enumerate(requested_times):
#         storm_id = fr["storm_id"][it]

#         if storm_id is None:
#             try:
#                 #  Currently data are used directly from gfs since grid is the same, but  resampling should be done in a new version
#                 datasets_gfs = gfs_forecast_0p25.download(
#                     param_list, lon_range, lat_range, [time_i, time_i], cycle_time
#                 )
#                 for ind, param in enumerate(param_list):
#                     if param == "wind":
#                         datasets[ind].u[it, :, :] = datasets_gfs[ind].u
#                         datasets[ind].v[it, :, :] = datasets_gfs[ind].v
#                     else:
#                         datasets[ind].val[it, :, :] = datasets_gfs[ind].val
#                     datasets[ind].unit = datasets_gfs[ind].unit
#                     datasets[ind].src.append("gfs")
#                 print(
#                     "No coamps-tc forecast available for {} --> using gfs forecast instead !!!".format(
#                         time_i
#                     )
#                 )

#             except Exception:
#                 print("gfs data not available")
#             continue

#         # Loop through requested parameters
#         for ind, param in enumerate(param_list):
#             dataset = datasets[ind]

#             try:
#                 okay = False
#                 makezeros = False
#                 #  First check if there are coamps-tc data available for this parameter
#                 if dss[storm_id][param_names[param]] is not None:
#                     okay = True
#                     model_t = pd.to_datetime(
#                         dss[storm_id][param_names[param]].time.values
#                     ).to_pydatetime()
#                     model_t_ind = np.where(model_t == time_i.replace(tzinfo=None))[0][0]
#                 else:
#                     makezeros = True

#                 if okay:
#                     dataset.coampsflag = True
#                     dataset.src.append("coamps_{}".format(storm_id))
#                     if param == "wind":
#                         data = dss[storm_id][param_names[param]]
#                         u = data["wind_u"]
#                         v = data["wind_v"]
#                         dataset.unit = u.units

#                         if np.any(np.isnan(u)) or np.any(
#                             np.isnan(v)
#                         ):  # check if there are nans and fill them up
#                             u = u.fillna(fill_values[param])
#                             v = v.fillna(fill_values[param])

#                         u = u.metpy.unit_array.squeeze()
#                         v = v.metpy.unit_array.squeeze()

#                         if reverse:
#                             dataset.u[it, :, :] = np.array(u[model_t_ind, ::-1, :])
#                             dataset.v[it, :, :] = np.array(v[model_t_ind, ::-1, :])
#                         else:
#                             dataset.u[it, :, :] = np.array(u[model_t_ind, :, :])
#                             dataset.v[it, :, :] = np.array(v[model_t_ind, :, :])
#                     else:
#                         if param == "barometric_pressure":
#                             var_name = "mslp"
#                         elif param == "precipitation":
#                             var_name = "rain"  # Check when available in scrubber tool!
#                         data = dss[storm_id][param_names[param]]
#                         val = data[var_name]

#                         # Added this check to ensure that pressure is in Pa
#                         if param == "barometric_pressure":
#                             if val.units == "mb":
#                                 val = val * 100
#                                 val.attrs["units"] = "Pa"
#                         dataset.unit = val.units

#                         if (
#                             param == "precipitation"
#                         ):  # Added this check to ensure that nan in precipitation are correctly interpreted
#                             val = val.where(val >= 0, np.nan)

#                         if np.any(
#                             np.isnan(val)
#                         ):  # check if there are nans and fill them up
#                             val = val.fillna(fill_values[param])

#                         val = np.array(val.metpy.unit_array.squeeze())

#                         if reverse:
#                             dataset.val[it, :, :] = np.array(val[model_t_ind, ::-1, :])
#                         else:
#                             dataset.val[it, :, :] = np.array(val[model_t_ind, :, :])

#                 elif makezeros:  # add zeros
#                     try:
#                         #  Currently data are used directly from gfs since grid is the same, but  resampling should be done in a new version
#                         datasets_gfs = gfs_forecast_0p25.download(
#                             param_list,
#                             lon_range,
#                             lat_range,
#                             [time_i, time_i],
#                             cycle_time,
#                         )
#                         if param == "wind":
#                             dataset.u[it, :, :] = datasets_gfs[ind].u
#                             dataset.v[it, :, :] = datasets_gfs[ind].v
#                         else:
#                             dataset.val[it, :, :] = datasets_gfs[ind].val
#                         dataset.unit = datasets_gfs[ind].unit
#                         dataset.src.append("gfs")
#                         print(
#                             param
#                             + " was not found on server for {} --> using gfs forecast data instead".format(
#                                 time_i
#                             )
#                         )
#                     except Exception:
#                         if param == "wind":
#                             dataset.u[:] = fill_values[param]
#                             dataset.v[:] = fill_values[param]
#                             dataset.unit = units[param]
#                             dataset.src.append("default_value")
#                             print(
#                                 param
#                                 + " was not found on server for {} --> using {} {} instead !!!".format(
#                                     time_i, fill_values[param], units[param]
#                                 )
#                             )

#                         else:
#                             dataset.val[:] = fill_values[param]
#                             dataset.unit = units[param]
#                             dataset.src.append("default_value")
#                             print(
#                                 param
#                                 + " was not found on server ... --> using {} {} instead !!!".format(
#                                     fill_values[param], units[param]
#                                 )
#                             )

#             except Exception:
#                 print("Could not download data")

#     return datasets


# def get_priority_storm():
#     """Find the storm with the highest priority from NOAA"""
#     tcvitals = requests.get("https://ftp.nhc.noaa.gov/atcf/com/tcvitals").text
#     splits = [line.split() for line in tcvitals.split("\n")[:-1]]  # last line is empty
#     priority = pd.DataFrame(splits)
#     priority = priority.iloc[
#         [i for i in range(len(priority)) if "L" in priority.iloc[i, 1]], :
#     ]  # check only in atlantic basin storms
#     try:
#         pr_st_noaa = priority.iloc[0, 1]  # define the top L storm as the priority one
#         return pr_st_noaa
#     except Exception:
#         return None


# def check_coamps(apikey, endpoint):
#     """Read metadata of available forecasts from endpoint and returns a dictionary of the available storms"""
#     headers = {"x-api-key": apikey}
#     response_API = requests.get("{}/status".format(endpoint), headers=headers)
#     data = json.loads(response_API.text)
#     storms = data["body"]["data"]["metget"]["coamps-tc"]

#     return storms


# def parse_domain_data(domain_list: list, level) -> dict:
#     import warnings

#     AVAILABLE_MODELS = {
#         "gfs": "gfs-ncep",
#         "nam": "nam-ncep",
#         "hwrf": "hwrf",
#         "coamps": "coamps-tc",
#     }

#     model = domain_list[0]
#     if "hwrf" in model:
#         storm = model.split("-")[1]
#         model = "hwrf"
#         warnings.warn(
#             "HWRF not fully supported yet. Use at your own risk.", RuntimeWarning
#         )
#     elif "coamps" in model:
#         storm = model.split("-")[1]
#         model = "coamps"
#         # warnings.warn(
#         #     "COAMPS not fully supported yet. Use at your own risk.", RuntimeWarning
#         # )

#     res = float(domain_list[1])
#     x0 = float(domain_list[2])
#     y0 = float(domain_list[3])
#     x1 = float(domain_list[4])
#     y1 = float(domain_list[5])

#     if model not in AVAILABLE_MODELS.keys():
#         raise RuntimeError("Specified model '" + model + "' is not available")

#     xmax = max(x0, x1)
#     xmin = min(x0, x1)
#     ymax = max(y0, y1)
#     ymin = min(y0, y1)
#     res = abs(res)
#     if res <= 0:
#         raise RuntimeError("Specified model resolution is invalid")

#     if model == "hwrf" or model == "coamps":
#         return {
#             "name": AVAILABLE_MODELS[model] + "-" + storm,
#             "service": AVAILABLE_MODELS[model],
#             "storm": storm,
#             "x_init": xmin,
#             "y_init": ymin,
#             "x_end": xmax,
#             "y_end": ymax,
#             "di": res,
#             "dj": res,
#             "level": level,
#         }
#     else:
#         return {
#             "name": model,
#             "service": AVAILABLE_MODELS[model],
#             "x_init": xmin,
#             "y_init": ymin,
#             "x_end": xmax,
#             "y_end": ymax,
#             "di": res,
#             "dj": res,
#             "level": level,
#         }


# def make_metget_request(endpoint, apikey, request_json):
#     headers = {"x-api-key": apikey}
#     r = requests.post(endpoint + "/build", headers=headers, json=request_json)
#     if r.status_code != 200:
#         raise RuntimeError(
#             "Request to MetGet was returned status code = " + str(r.status_code)
#         )
#     return_data = json.loads(r.text)
#     data_id = return_data["body"]["request_id"]
#     status_code = return_data["statusCode"]
#     if status_code != 200:
#         with open("metget.debug", "a") as f:
#             f.write("[WARNING]: MetGet returned status code " + str(status_code) + "\n")
#             f.write(str(return_data["body"]["error_text"]))
#     return data_id, status_code


# def download_metget_data(data_id, endpoint, apikey, sleeptime, max_wait):
#     from datetime import datetime, timedelta

#     # from netCDF4 import Dataset
#     import xarray as xr

#     # ...Status check
#     headers = {"x-api-key": apikey}
#     request_json = {"request": data_id}

#     # ...Wait time
#     end_time = datetime.utcnow() + timedelta(hours=max_wait)

#     # ...Wait for request data to appear
#     tries = 0
#     data_ready = False
#     status = None
#     # print("Waiting for request id: ", data_id, flush=True)
#     while datetime.utcnow() <= end_time:
#         tries += 1
#         try:
#             # print(
#             #     "["
#             #     + datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")
#             #     + "]: Checking request status...(n="
#             #     + str(tries)
#             #     + "): ",
#             #     flush=True,
#             #     end="",
#             # )
#             response = requests.post(
#                 endpoint + "/check", headers=headers, json=request_json
#             )
#             json_response = json.loads(response.text)
#             status = json_response["body"]["status"]
#             data_url = json_response["body"]["destination"]
#             # print(status, flush=True)
#             if status == "completed":
#                 # ...Parse the return to get data
#                 data_ready = True
#                 flist_url = data_url + "/filelist.json"
#                 u = requests.get(flist_url)
#                 if u.status_code == 200:
#                     return_data = json.loads(u.text)
#                     # Turned this off to avoid saving this file
#                     # with open("filelist.json", "w") as jsonfile:
#                     #     jsonfile.write(
#                     #         json.dumps(return_data, indent=2, sort_keys=True)
#                     #     )
#                     break
#             elif status == "error":
#                 # print("Request could not be completed")
#                 return
#             else:
#                 time.sleep(sleeptime)
#                 continue
#         except KeyboardInterrupt:
#             # print("[ERROR]: Process was ended by the user")
#             raise

#     # ...Download files
#     if data_ready:
#         file_list = return_data["output_files"]
#         ds_list = []
#         for f in file_list:
#             # print("[{:s}]: Getting file: {:s}".format(time_stamp, f), flush=True)
#             url = data_url + "/" + f
#             # ds = Dataset('name', memory=requests.get(url).content)
#             try:
#                 ds = xr.open_dataset(
#                     url + "#mode=bytes"
#                 )  # Added this last part to allow opening with xarray
#             except Exception:
#                 from netCDF4 import Dataset

#                 data = requests.get(url).content
#                 ds0 = Dataset("temp", memory=data)
#                 ds = xr.open_dataset(xr.backends.NetCDF4DataStore(ds0))
#             ds_list.append(ds)

#         return ds_list
#     else:
#         if status == "restore":
#             print(
#                 "[WARNING]: Data for request "
#                 + data_id
#                 + " did not become ready before the max-wait time expired. You can rerun and ask for this request by id"
#             )
#         elif status == "running":
#             print(
#                 "[WARNING]: Data for request "
#                 + data_id
#                 + " is still being constructed when the max-wait time expired. Please check on it later"
#             )
#         elif status == "queued":
#             print(
#                 "[WARNING]: Data for request "
#                 + data_id
#                 + " is still queued. If this does not change soon, please contact an administrator"
#             )
#         else:
#             print("[ERROR]: Data has not become available due to an unknown error")
#         return


# def met_get(**kwargs):
#     """ "based on the get_metget_data.py script but instead of downloading a ntecdf files it returns an xr.dataset
#     in memory"""
#     import getpass
#     import socket
#     import types

#     AVAILABLE_VARIABLES = {"wind_pressure", "rain", "temperature", "humidity", "ice"}
#     AVAILABLE_FORMATS = {"ascii", "owi-ascii", "adcirc-netcdf", "hec-netcdf", "delft3d"}

#     args = types.SimpleNamespace()

#     for name in kwargs.keys():
#         exec('args.{} = kwargs["{}"]'.format(name, name))

#     optional = {
#         "analysis": False,
#         "multiple_forecasts": True,
#         "format": "hec-netcdf",
#         "variable": "wind_pressure",
#         "check_interval": 30,
#         "max_wait": 24,
#         "strict": True,
#         "backfile": True,
#         "epsg": 4326,
#         "dryrun": False,
#         "request": None,
#         "compression": False,
#         "endpoint": None,
#         "apikey": None,
#     }

#     for key in optional.keys():
#         if key not in kwargs:
#             exec('args.{} = optional["{}"]'.format(key, key))

#     if not args.endpoint:
#         if "METGET_ENDPOINT" not in os.environ:
#             raise RuntimeError("No endpoint found.")
#         else:
#             endpoint = os.environ["METGET_ENDPOINT"]
#     else:
#         endpoint = args.endpoint

#     if not args.apikey:
#         if "METGET_API_KEY" not in os.environ:
#             raise RuntimeError("No API key was found.")
#         else:
#             apikey = os.environ["METGET_API_KEY"]
#     else:
#         apikey = args.apikey

#     # ...Check for required arguments
#     if not args.request:
#         if not args.start:
#             print("[ERROR]: Must provide '--start'")
#             exit(1)
#         if not args.end:
#             print("[ERROR]: Must provide '--end'")
#             exit(1)
#         if not args.timestep:
#             print("[ERROR]: Must provide '--timestep'")
#             exit(1)

#         # ...Building the request
#         domains = []
#         idx = 0
#         for d in args.domain:
#             j = parse_domain_data(d, idx)
#             domains.append(j)
#             idx += 1

#         if args.format not in AVAILABLE_FORMATS:
#             print("ERROR: Invalid output format selected")
#             exit(1)

#         if args.format == "delft3d" or args.format == "hec-netcdf":
#             if len(domains) > 1:
#                 print(
#                     "[ERROR]: "
#                     + args.format
#                     + " does not support more than one domain."
#                 )
#                 exit(1)

#         if args.variable not in AVAILABLE_VARIABLES:
#             print("ERROR: Invalid variable selected")
#             exit(1)

#         request_from = getpass.getuser() + "." + socket.gethostname()
#         request_data = {
#             "version": "0.0.1",
#             "creator": request_from,
#             "background_pressure": 1013.0,
#             "backfill": True,
#             "nowcast": args.analysis,
#             "multiple_forecasts": args.multiple_forecasts,
#             "start_date": str(args.start),
#             "end_date": str(args.end),
#             "format": args.format,
#             "data_type": args.variable,
#             "time_step": args.timestep,
#             "domains": domains,
#             "compression": args.compression,
#             "epsg": args.epsg,
#             "filename": args.output,
#         }
#         if args.strict:
#             request_data["strict"] = True
#         if args.dryrun:
#             request_data["dry_run"] = True

#         data_id, status_code = make_metget_request(endpoint, apikey, request_data)
#         if not args.dryrun and status_code == 200:
#             ds_list = download_metget_data(
#                 data_id, endpoint, apikey, args.check_interval, args.max_wait
#             )
#         else:
#             print(status_code)
#             return None

#     else:
#         ds_list = download_metget_data(
#             args.request, endpoint, apikey, args.check_interval, args.max_wait
#         )

#     return ds_list
