# -*- coding: utf-8 -*-
import datetime

import numpy as np
import pandas as pd
import tomli
from metget.metget_build import MetGetBuildRest
from pyproj import CRS

from cht_meteo import gfs_forecast_0p25

from .metget_utils import *


class Dataset:
    def __init__(self):
        self.quantity = None
        self.unit = None
        self.time = []
        self.x = None
        self.y = None
        self.crs = None
        self.val = None
        self.source = []


def download(
    param_list,
    lon_range,
    lat_range,
    time_range,
    cycle_time,
    resolution=0.25,
    dt=3,
    config_path=None,
):
    """Function to download coamps-tc forecasts using the metget tool api"""

    metget = MetGet(config_path=config_path)

    # Get info
    df = metget.get_coamps_quasi_hindcast(time_range, dt=dt)
    requested_times = df.index
    ntime = len(requested_times)

    coamps_storms = metget.get_forecasts(
        model="coamps", start=time_range[0], end=time_range[1]
    )
    # Get the storm names that have a forecast available for the requested cycle
    storms = [
        name
        for name in coamps_storms.keys()
        if cycle_time.strftime("%Y-%m-%d %H:%M:%S")
        in coamps_storms[name]["cycles_complete"]
    ]

    # If there are no storms available raise Error
    if len(storms) == 0:
        raise ImportError(
            f"There are no COAMPS-TC forecasts available for cycle: {cycle_time.strftime('%Y-%m-%d %H:%M:%S')}!"
        )

    # Get priority storm
    priority_storm = metget.priority_storm
    if priority_storm not in storms:
        raise ValueError(
            f"There is no forecast available for priority storm {priority_storm} for cycle: {cycle_time.strftime('%Y-%m-%d %H:%M:%S')}!"
        )

    # Get variable names in metget
    variables = list(
        np.unique([var[1].metget_name for var in metget.meteo_variables.items()])
    )  # ['wind_pressure', 'rain']

    dss = {}  # prepare dictionary to save the datasets
    metadata = {}

    # Prepare metget domain
    domain = (
        [f"coamps-{priority_storm}"]
        + [resolution]
        + [lon_range[0], lat_range[0], lon_range[1], lat_range[1]]
    )
    for var in variables:
        ds, meta = metget.get_meteo(
            time_range=time_range,
            dt=dt,
            domain=domain,
            var=var,
            multiple_forecasts=False,
        )
        dss[var] = ds
        metadata[var] = meta

    # if no data could be downloaded (CHECK AGAIN WHAT SHOULD HAPPEN IN THIS CASE)
    if all([dss[k] is None for k in list(dss.keys())]):
        raise ValueError("Could not find any coamps-tc data in requested time range!")

    # Do checks
    for key in dss.keys():
        if dss[key] is not None:
            data0 = dss[key]
            lon = np.array(data0["lon"])
            lat = np.array(data0["lat"])
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
        dataset.time = pd.to_datetime(
            [tt.replace(tzinfo=None) for tt in requested_times]
        ).to_pydatetime()
        if dataset.quantity == "wind":
            dataset.u = np.empty((ntime, nrows, ncols))
            dataset.u[:] = np.nan
            dataset.v = np.empty((ntime, nrows, ncols))
            dataset.v[:] = np.nan
        else:
            dataset.val = np.empty((ntime, nrows, ncols))
            dataset.val[:] = np.nan
        datasets.append(dataset)

    # Loop through requested parameters
    for ind, param in enumerate(param_list):
        dataset = datasets[ind]
        #  First check if there are data available for this parameter
        if dss[metget.meteo_variables[param].metget_name] is None:
            raise KeyError

        # Check naming of files to get information on cycle used
        naming_format = metadata[metget.meteo_variables[param].metget_name][
            "input_files"
        ][f"coamps-tc-{priority_storm}"][0]
        if naming_format.split("_")[0] == "coamps-tc":
            cycles_used = [
                name.split("_")[2]
                for name in metadata[metget.meteo_variables[param].metget_name][
                    "input_files"
                ][f"coamps-tc-{priority_storm}"]
            ]
            taus = [
                name.split("_")[3].split(".")[0].split("tau")[1]
                for name in metadata[metget.meteo_variables[param].metget_name][
                    "input_files"
                ][f"coamps-tc-{priority_storm}"]
            ]
            times = [
                datetime.datetime.strptime(cycle, "%Y%m%d%H")
                + datetime.timedelta(hours=int(tau))
                for (cycle, tau) in zip(cycles_used, taus)
            ]
            times = [t.strftime("%Y%m%d%H%M") for t in times]
        else:
            raise ValueError

        for it, time_i in enumerate(requested_times):
            inputs = [
                name
                for i, name in enumerate(
                    metadata[metget.meteo_variables[param].metget_name]["input_files"][
                        f"coamps-tc-{priority_storm}"
                    ]
                )
                if time_i.strftime("%Y%m%d%H%M") in times[i]
            ]
            # get cycle info
            cycle_used = inputs[0].split("_")[2]
            cycle_hour = inputs[0].split("_")[3].split(".")[0].split("tau")[1]
            dataset.source.append(
                f"coamps_tc_{priority_storm}_{cycle_used}z_{cycle_hour}"
            )

        # Get metadata of input used
        dataset.config = metadata[metget.meteo_variables[param].metget_name]["input"]

        if param == "wind":
            data = dss[metget.meteo_variables[param].metget_name]
            u = data["wind_u"]
            v = data["wind_v"]
            dataset.unit = u.units

            if np.any(np.isnan(u)) or np.any(
                np.isnan(v)
            ):  # check if there are nans and fill them up
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
                var_name = "precipitation"
            data = dss[metget.meteo_variables[param].metget_name]
            val = data[var_name]

            # Added this check to ensure that pressure is in Pa
            if param == "barometric_pressure":
                if val.units == "mb":
                    val = val * 100
                    val.attrs["units"] = "Pa"
            dataset.unit = val.units

            if (
                param == "precipitation"
            ):  # Added this check to ensure that nan in precipitation are correctly interpreted
                val = val.where(val >= 0, np.nan)

            if np.any(np.isnan(val)):  # check if there are nans and fill them up
                val = val.fillna(metget.meteo_variables[param].fill_value)

            val = np.array(val.metpy.unit_array.squeeze())

            if reverse:
                dataset.val = np.array(val[:, ::-1, :])
            else:
                dataset.val = np.array(val)

    return datasets
