# -*- coding: utf-8 -*-
from metget.metget_build import MetGetBuildRest
import requests
import tomli
import pandas as pd
import numpy as np
from pyproj import CRS
import datetime
from cht_meteo import gfs_forecast_0p25
from .coamps_utils import *

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


def download(param_list, lon_range, lat_range, time_range, cycle_time, resolution=0.25, config_path=None):
    """Function to download coamps-tc forecasts using the metget tool api
    Right now resolution is hardcoded to 0.25 degrees but this can change"""

    fill_values = dict(wind=0.0, precipitation=0.0, barometric_pressure=102000.0)  # fill values to use in case of nan
    units = dict(wind='m/s', precipitation='kg.m-2.hour-1', barometric_pressure='Pa')  # units

    if not config_path:
        raise ValueError("Configuration data need to be provided")
    else:
        with open(config_path, mode="rb") as fp:
            config = tomli.load(fp)
            
    # Read data from the config file
    apikey = config["apikey"]
    endpoint = config["endpoint"]
    api_version = config["api_version"]
    if "tau" in config:
        tau = config["tau"]
    else:
        tau = 0
    if "override" in config:
        override = True
    else:
        override = False

    # Connect parameters names with coamps-tc names
    param_names = {'wind': 'wind_pressure', 'barometric_pressure': 'wind_pressure', 'precipitation': 'rain'}
    variables = list(np.unique([param_names[name] for name in param_list])) # ['wind_pressure', 'rain']

    # Get individual times requested
    # TODO how is time-step handled in cosmos?
    requested_times = pd.date_range(start=time_range[0],
                      end=time_range[1],
                      freq='3H').to_pydatetime().tolist()
    requested_times = [ti.replace(tzinfo=datetime.timezone.utc) for ti in requested_times]

    timestep = int((requested_times[1] - requested_times[0]).total_seconds())  # get time-step of requested time assuming that the time-step stays the same
    ntime = len(requested_times)

    # Get the available storms in the Coamps-TC forecast
    try:
        coamps_storms = check_coamps(apikey, endpoint, time_range[0])
    except:
        coamps_storms = {}
    
    storms = [name for name in coamps_storms.keys() if int(name.split("L")[0]) < 90]
    storms.sort()
    
    # Check if forecast exists for the selected dates and make a dataframe
    fr = pd.DataFrame(index=requested_times)
    fr["storms"] = None

    for i, ti in enumerate(requested_times):
        fr.loc[ti, "storms"] = []
        for j, storm_id in enumerate(storms):
            storm = coamps_storms[storm_id]                 
            ####### CHECK WHICH times to use for range ####################
            t1 = date_transform(storm['min_forecast_date'])
            t2 = date_transform(storm['max_forecast_date'])
            # t1 = date_transform(storm['latest_complete_forecast_start'])
            # t2 = date_transform(storm['latest_complete_forecast_end'])
            belongs = t1 <= ti <= t2
            if belongs:
                fr.loc[ti, "storms"].append(storm_id)

    # Check if there is a priority storm given
    if "priority_storm" in config:
        if config["priority_storm"] == "tc_vitals":
            # Get storm priority from noaa
            priority_storm = tc_vitals_storm()
        else:
            priority_storm = config["priority_storm"]
    else:
        priority_storm=None

    # Method to use when multiple storm are available without priority storm
    if "multi_storms_method" in config:
        method = config["multi_storms_method"]
    else:
        method = "last"

    fr["storm_id"] = None
    # Loop through time steps and choose a single storm to use 
    for i, ti in enumerate(requested_times):
        # # Check if for the chosen storm the time is during the skip time (tau) of the first forecast
        # for sn in fr.loc[ti, "storms"]:
        #     t0 = date_transform(coamps_storms[sn]["first_available_cycle"])
        #     t0_tau = t0 + datetime.timedelta(hours=tau)
        #     if ti <= t0_tau:
        #         fr.loc[ti, "storms"].remove(sn)
        storm_name = None # if there is no coamps storm make None so gfs is used
        if priority_storm and override:
            storm_name = priority_storm
        elif priority_storm in fr.loc[ti, "storms"]:  # if the priority storm is available that timestep use it
            storm_name = priority_storm
        else:
            # else used method defined
            if method == 'first':
                storm_name = fr.loc[ti, "storms"][0]
            elif method == 'last':
                storm_name = fr.loc[ti, "storms"][-1]

        # Assign storm
        fr.loc[ti, "storm_id"] = storm_name
        
    storm_log = np.array([v is not None for v in fr['storm_id'].values]) # check for which times there are forecasts from coamps-tc
    if all(~storm_log):  # If no coamps-tc data are available, gfs data are downloaded
        datasets = gfs_forecast_0p25.download(param_list, lon_range, lat_range, time_range, cycle_time)
        print("Could not find any coamps-tc data in requested time range! Using gfs instead")
        return datasets
    
    # Else if there are coamps-tc data see which storms (models) to use
    storms_id = np.unique(fr['storm_id'][storm_log])
    dss = {}  # prepare dictionary to save the datasets
    metadata = {}
    # Loop for the case there is more than one storm in requested times and get storm forecast
    for i, st in enumerate(storms_id):
        dss[st] = {}
        metadata[st] = {}
        if not override:
            storm = coamps_storms[st]
            t1 = datetime.datetime.strptime(storm['min_forecast_date'], "%Y-%m-%d %H:%M:%S").replace(tzinfo=datetime.timezone.utc)
            t2 = datetime.datetime.strptime(storm['max_forecast_date'], "%Y-%m-%d %H:%M:%S").replace(tzinfo=datetime.timezone.utc)
            # This parts ensures that requested data are within the model range
            t1 = max(requested_times[0], t1)
            t2 = min(requested_times[-1], t2)
        else:
            t1 = requested_times[0]
            t2 = requested_times[-1]
        if t1 == t2:  # scrubber tool cannot handle single time
            t1 = t1 - datetime.timedelta(hours=3)
        # Prepare metget domain
        domain = [f"coamps-{st}"] + [resolution] + [lon_range[0], lat_range[0], lon_range[1], lat_range[1]]

        for ii, var in enumerate(variables):
            # ...Building the request
            request_data = MetGetBuildRest.generate_request_json(
                start_date=t1.strftime("%Y%m%d %H%M%S"),
                end_date=t2.strftime("%Y%m%d %H%M%S"),
                format='hec-netcdf',
                data_type=var,
                time_step=timestep,
                domains=MetGetBuildRest.parse_command_line_domains([domain], tau),
                epsg=4326,
                filename=f"{st}_{var}",
                backfill=True,
                # nowcast=True,
                multiple_forecasts=True,
                # compression=True,
                save_json_request=False,
                # dry_run=True,
                strict=True,
            )
            client = MetGetBuildRest(endpoint, apikey, api_version)
            data_id, status_code = client.make_metget_request(request_data)
            urls, meta = client.download_metget_data(
                data_id,
                30, # sleep_time in seconds
                3, # max_wait in hours
                output_directory=None,
                return_only_url=True,
            )
            ds_list = []
            for url in urls:
                ds = get_da_from_url(url)
                ds_list.append(ds)
            if ds_list is not None:
                dss[st][var] = ds_list[0]
                metadata[st][var] = meta
            else:
                dss[st][var] = None

    # if no data could be downloaded (CHECK AGAIN WHAT SHOULD HAPPEN IN THIS CASE)
    if all([dss[k][kk] is None for k in list(dss.keys()) for kk in list(dss[k].keys())]):
        datasets = gfs_forecast_0p25.download(param_list, lon_range, lat_range, time_range, cycle_time)
        print("Could not find any coamps-tc data in requested time range! Using gfs instead")
        return datasets
    else: # Else if there are available data get lon lat info
        for key in dss[storms_id[0]].keys():
            if dss[storms_id[0]][key] is not None:
                data0 = dss[storms_id[0]][key]
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

    # Prepare datasets
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


    for it, time_i in enumerate(requested_times):
        storm_id = fr['storm_id'][it]

        if storm_id is None:
            try:
                #  Currently data are used directly from gfs since grid is the same, but  resampling should be done in a new version
                datasets_gfs = gfs_forecast_0p25.download(param_list, lon_range, lat_range,
                                                          [time_i, time_i], cycle_time)
                for ind, param in enumerate(param_list):
                    if param == "wind":
                        datasets[ind].u[it, :, :] = datasets_gfs[ind].u
                        datasets[ind].v[it, :, :] = datasets_gfs[ind].v
                    else:
                        datasets[ind].val[it, :, :] = datasets_gfs[ind].val
                    datasets[ind].unit = datasets_gfs[ind].unit
                    datasets[ind].src.append('gfs')
                print("No coamps-tc forecast available for {} --> using gfs forecast instead !!!".format(time_i))

            except:
                print('gfs data not available')
            continue

        # Loop through requested parameters
        for ind, param in enumerate(param_list):
            dataset = datasets[ind]

            try:
                okay = False
                makezeros = False
                #  First check if there are coamps-tc data available for this parameter
                if dss[storm_id][param_names[param]] is not None:
                    okay = True
                    model_t = pd.to_datetime(dss[storm_id][param_names[param]].time.values).to_pydatetime()
                    model_t_ind = np.where(model_t == time_i.replace(tzinfo=None))[0][0]
                else:
                    makezeros = True

                if okay:
                    # Get metadata of input used
                    dataset.config = metadata[storm_id][param_names[param]]["input"]
                    inputs = [name for name in metadata[storm_id][param_names[param]]["input_files"][f"coamps-tc-{storm_id}"] if time_i.strftime("%Y%m%d%H%M") in name]
                    input_name = inputs[0].split(".")[2]
                    cycle_used = input_name.split("_")[2]
                    cycle_hour = input_name.split("_")[3]
                    dataset.source.append(f'coamps_tc_{storm_id}_{cycle_used}z_{cycle_hour}')
                    if param == "wind":
                        data = dss[storm_id][param_names[param]]
                        u = data['wind_u']
                        v = data['wind_v']
                        dataset.unit = u.units
                        
                        if np.any(np.isnan(u)) or np.any(np.isnan(v)):  # check if there are nans and fill them up
                            u = u.fillna(fill_values[param])
                            v = v.fillna(fill_values[param])
                        
                        u = u.metpy.unit_array.squeeze()
                        v = v.metpy.unit_array.squeeze()

                        if reverse:
                            dataset.u[it,:,:] = np.array(u[model_t_ind, ::-1, :])
                            dataset.v[it,:,:] = np.array(v[model_t_ind, ::-1, :])
                        else:
                            dataset.u[it,:,:] = np.array(u[model_t_ind, :, :])
                            dataset.v[it,:,:] = np.array(v[model_t_ind, :, :])
                    else:
                        if param == "barometric_pressure":
                            var_name = "mslp"
                        elif param == "precipitation":
                            var_name = 'rain' # Check when available in scrubber tool!
                        data = dss[storm_id][param_names[param]]
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
                            val = val.fillna(fill_values[param])
                            
                        val = np.array(val.metpy.unit_array.squeeze())

                        if reverse:
                            dataset.val[it, :, :] = np.array(val[model_t_ind, ::-1, :])
                        else:
                            dataset.val[it, :, :] = np.array(val[model_t_ind, :, :])

                elif makezeros: # add zeros
                    try:
                        #  Currently data are used directly from gfs since grid is the same, but  resampling should be done in a new version
                        datasets_gfs = gfs_forecast_0p25.download(param_list, lon_range, lat_range,
                                                               [time_i, time_i], cycle_time)
                        if param == "wind":
                            dataset.u[it, :, :] = datasets_gfs[ind].u
                            dataset.v[it, :, :] = datasets_gfs[ind].v
                        else:
                            dataset.val[it, :, :]  = datasets_gfs[ind].val
                        dataset.unit = datasets_gfs[ind].unit
                        dataset.src.append('gfs')
                        print(param + " was not found on server for {} --> using gfs forecast data instead".format(time_i))
                    except:
                        if param == "wind":
                            dataset.u[:] = fill_values[param]
                            dataset.v[:] = fill_values[param]
                            dataset.unit = units[param]
                            dataset.src.append('default_value')
                            print(param + " was not found on server for {} --> using {} {} instead !!!".format(time_i, fill_values[param], units[param]))

                        else:
                            dataset.val[:] = fill_values[param]
                            dataset.unit = units[param]
                            dataset.src.append('default_value')
                            print(param + " was not found on server ... --> using {} {} instead !!!".format(fill_values[param], units[param]))

            except:
                print("Could not download data")

    return datasets
