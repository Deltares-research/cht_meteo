import os
import pandas as pd
import numpy as np

def write_to_delft3d(
    dataset,
    file_name,
    version="1.03",
    path=None,
    header_comments=False,
    refdate=None,
    parameters=None,
    time_range=None,
):
    
    
    # Convert to datetime
    time = pd.to_datetime(dataset.ds.time.values)

    if not refdate:
        refdate = time[0]

    if not time_range:
        time_range = [time[0], time[-1]]

    if not parameters:
        parameters = []
        for var_name in dataset.ds.data_vars:
            if var_name not in ["lon", "lat", "time"]:
                parameters.append(var_name)

    if dataset.crs.is_geographic:
        grid_unit = "degrees"
    else:
        grid_unit = "m"

    files = []
    for param in parameters:
        # Look up index of this parameter
        for ind, var_name in dataset.ds.data_vars:
            if param == var_name:
                q = var_name
                break
        if param == "wind":
            file = {}
            file["data"] = q.u
            file["ext"] = "amu"
            file["quantity"] = "x_wind"
            file["unit"] = "m s-1"
            file["fmt"] = "%6.1f"
            files.append(file)
            file = {}
            file["data"] = q.v
            file["ext"] = "amv"
            file["quantity"] = "y_wind"
            file["unit"] = "m s-1"
            file["fmt"] = "%6.1f"
            files.append(file)
        elif param == "barometric_pressure":
            file = {}
            file["data"] = q.val
            file["ext"] = "amp"
            file["quantity"] = "air_pressure"
            file["unit"] = "Pa"
            file["fmt"] = "%7.0f"
            files.append(file)
        elif param == "precipitation":
            file = {}
            file["data"] = q.val
            file["ext"] = "ampr"
            file["quantity"] = "precipitation"
            file["unit"] = "mm h-1"
            file["fmt"] = "%7.1f"
            files.append(file)

    # if dataset.quantity == "x_wind":
    #     unit = "m s-1"
    #     ext  = "amu"
    #     fmt  = "%6.1f"
    # elif dataset.quantity == "y_wind":
    #     unit = "m s-1"
    #     ext  = "amv"
    #     fmt  = "%6.1f"
    # elif dataset.quantity == "air_pressure":
    #     unit = "Pa"
    #     ext  = "amp"
    #     fmt  = "%7.0f"
    # elif dataset.quantity == "air_temperature":
    #     unit = "Celsius"
    #     ext  = "amt"
    #     fmt  = "%7.1f"
    # elif dataset.quantity == "relative_humidity":
    #     unit = "%"
    #     ext  = "amr"
    #     fmt  = "%7.1f"
    # elif dataset.quantity == "cloudiness":
    #     unit = "%"
    #     ext  = "amc"
    #     fmt  = "%7.1f"
    # elif dataset.quantity == "sw_radiation_flux":
    #     unit = "W/m2"
    #     ext  = "ams"
    #     fmt  = "%7.1f"
    # elif dataset.quantity == "precipitation":
    #     unit = "mm/h"
    #     ext  = "ampr"
    #     fmt  = "%7.1f"

    for file in files:
        ncols = len(dataset.x)
        nrows = len(dataset.y)

        dx = (dataset.x[-1] - dataset.x[0]) / (len(dataset.x) - 1)
        dy = (dataset.y[-1] - dataset.y[0]) / (len(dataset.y) - 1)

        if path:
            full_file_name = os.path.join(path, file_name + "." + file["ext"])
        else:
            full_file_name = file_name + "." + file["ext"]

        fid = open(full_file_name, "w")

        if header_comments:
            fid.write("### START OF HEADER\n")
            fid.write(
                "### All text on a line behind the first # is parsed as commentary\n"
            )
            fid.write("### Additional comments\n")

        fid.write(
            "FileVersion      =   "
            + version
            + "                                               # Version of meteo input file, to check if the newest file format is used\n"
        )
        fid.write(
            "filetype         =   meteo_on_equidistant_grid                          # Type of meteo input file: meteo_on_flow_grid, meteo_on_equidistant_grid, meteo_on_curvilinear_grid or meteo_on_spiderweb_grid\n"
        )
        fid.write(
            "NODATA_value     =   -999                                               # Value used for undefined or missing data\n"
        )
        fid.write("n_cols           =   " + str(ncols) + "\n")
        fid.write("n_rows           =   " + str(nrows) + "\n")
        fid.write("grid_unit        =   " + grid_unit + "\n")
        #            fid.write("x_llcorner       =   " + str(min(dataset.x)) + "\n")
        #            fid.write("y_llcorner       =   " + str(min(dataset.y)) + "\n")
        fid.write("x_llcorner       =   " + str(min(dataset.x) - 0.5 * dx) + "\n")
        fid.write("y_llcorner       =   " + str(min(dataset.y) - 0.5 * dy) + "\n")
        if version == "1.02":
            fid.write("value_pos       =    corner\n")
        fid.write("dx               =   " + str(dx) + "\n")
        fid.write("dy               =   " + str(dy) + "\n")
        fid.write(
            "n_quantity       =   1                                                  # Number of quantities prescribed in the file\n"
        )
        fid.write("quantity1        =   " + file["quantity"] + "\n")
        fid.write("unit1            =   " + file["unit"] + "\n")
        if header_comments:
            fid.write("### END OF HEADER\n")

        # Add extra blocks if data does not cover time range
        if time[0] > time_range[0]:
            dt = time_range[0] - refdate
            tim = dt.total_seconds() / 60
            val = np.flipud(file["data"][0, :, :])
            # Skip blocks with only nans
            if not np.all(np.isnan(val)):
                val[val == np.nan] = -999.0
                fid.write(
                    "TIME = "
                    + str(tim)
                    + " minutes since "
                    + refdate.strftime("%Y-%m-%d %H:%M:%S")
                    + " +00:00\n"
                )
                np.savetxt(fid, val, fmt=file["fmt"])

        for it, time in enumerate(time):
            dt = time - refdate
            tim = dt.total_seconds() / 60
            val = np.flipud(file["data"][it, :, :])

            if param == "wind":
                if np.max(val) > 1000.0:
                    val = np.zeros_like(
                        val
                    )  # Weird array, don't trust. Set everything to zeros.
                    val[np.where(val == 0.0)] = np.nan
                    print(
                        "Warning! Wind speed > 1000 m/s at "
                        + time.strftime("%Y-%m-%d %H:%M:%S")
                        + " !"
                    )
                if np.min(val) < -1000.0:
                    val = np.zeros_like(
                        val
                    )  # Weird array, don't trust. Set everything to zeros.
                    print(
                        "Warning! Wind speed > 1000 m/s at "
                        + time.strftime("%Y-%m-%d %H:%M:%S")
                        + " !"
                    )
                    val[np.where(val == 0.0)] = np.nan
            if param == "barometric_pressure":
                if np.max(val) > 200000.0:
                    val = np.zeros_like(
                        val
                    )  # Weird array, don't trust. Set everything to zeros.
                    val[np.where(val == 0.0)] = np.nan
                if np.min(val) < 10000.0:
                    val = np.zeros_like(
                        val
                    )  # Weird array, don't trust. Set everything to zeros.
                    val[np.where(val == 0.0)] = np.nan
            if param == "precipitation":
                if np.nanmax(val) > 1000.0:
                    val = np.zeros_like(
                        val
                    )  # Weird array, don't trust. Set everything to zeros.
                    print(
                        "Warning! Precipitation exceeds 1000 mm/h at "
                        + time.strftime("%Y-%m-%d %H:%M:%S")
                        + " !"
                    )
                    val[np.where(val == 0.0)] = np.nan
                if np.nanmin(val) < 0.0:
                    val[np.where(val < 0.0)] = 0.0

            if np.all(np.isnan(val)):
                if it > 0:
                    print(
                        "Warning! Only NaNs found for "
                        + param
                        + " at "
                        + time.strftime("%Y-%m-%d %H:%M:%S")
                        + " ! Using data from previous time."
                    )

                else:
                    if param == "wind" or param == "precipitation":
                        print(
                            "Warning! Only NaNs found for "
                            + param
                            + " at "
                            + time.strftime("%Y-%m-%d %H:%M:%S")
                            + " ! Setting values to 0.0 !"
                        )
                        val = np.zeros_like(val)
                    elif param == "barometric_pressure":
                        print(
                            "Warning! Only NaNs found for "
                            + param
                            + " at "
                            + time.strftime("%Y-%m-%d %H:%M:%S")
                            + " ! Setting values to 101300.0 !"
                        )
                        val = np.zeros_like(val) + 101300.0

            fid.write(
                "TIME = "
                + str(tim)
                + " minutes since "
                + refdate.strftime("%Y-%m-%d %H:%M:%S")
                + " +00:00\n"
            )
            np.savetxt(fid, val, fmt=file["fmt"])

        # Add extra blocks if data does not cover time range
        if time[-1] < time_range[1]:
            dt = time_range[1] - refdate
            tim = dt.total_seconds() / 60
            val = np.flipud(file["data"][-1, :, :])
            # Skip blocks with only nans
            if not np.all(np.isnan(val)):
                val[val == np.nan] = -999.0
                fid.write(
                    "TIME = "
                    + str(tim)
                    + " minutes since "
                    + refdate.strftime("%Y-%m-%d %H:%M:%S")
                    + " +00:00\n"
                )
                np.savetxt(fid, val, fmt=file["fmt"])

        fid.close()

def write_wind_to_json(self, file_name, time_range=None, iref=1, js=False):
    import json

    if not time_range:
        time_range = []
        time_range.append(time[0])
        time_range.append(time[-1])

    data = []

    header = {
        "discipline": 0,
        "disciplineName": "Meteorological products",
        "gribEdition": 2,
        "gribLength": 76420,
        "center": 7,
        "centerName": "US National Weather Service - NCEP(WMC)",
        "subcenter": 0,
        "refTime": "2016-04-30T06:00:00.000Z",
        "significanceOfRT": 1,
        "significanceOfRTName": "Start of forecast",
        "productStatus": 0,
        "productStatusName": "Operational products",
        "productType": 1,
        "productTypeName": "Forecast products",
        "productDefinitionTemplate": 0,
        "productDefinitionTemplateName": "Analysis/forecast at horizontal level/layer at a point in time",
        "parameterCategory": 2,
        "parameterCategoryName": "Momentum",
        "parameterNumber": 2,
        "parameterNumberName": "U-component_of_wind",
        "parameterUnit": "m.s-1",
        "genProcessType": 2,
        "genProcessTypeName": "Forecast",
        "forecastTime": 0,
        "surface1Type": 103,
        "surface1TypeName": "Specified height level above ground",
        "surface1Value": 10.0,
        "surface2Type": 255,
        "surface2TypeName": "Missing",
        "surface2Value": 0.0,
        "gridDefinitionTemplate": 0,
        "gridDefinitionTemplateName": "Latitude_Longitude",
        "numberPoints": 65160,
        "shape": 6,
        "shapeName": "Earth spherical with radius of 6,371,229.0 m",
        "gridUnits": "degrees",
        "resolution": 48,
        "winds": "true",
        "scanMode": 0,
        "nx": 360,
        "ny": 181,
        "basicAngle": 0,
        "subDivisions": 0,
        "lo1": 0.0,
        "la1": 90.0,
        "lo2": 359.0,
        "la2": -90.0,
        "dx": 1.0,
        "dy": 1.0,
    }

    header["lo1"] = float(min(dataset.x) + 360.0)
    header["lo2"] = float(max(dataset.x) + 360.0)
    header["la1"] = float(max(dataset.y))
    header["la2"] = float(min(dataset.y))
    header["dx"] = float(dataset.x[1] - dataset.x[0])
    header["dy"] = float(dataset.y[1] - dataset.y[0])
    header["nx"] = len(dataset.x)
    header["ny"] = len(dataset.y)
    header["numberPoints"] = len(dataset.x) * len(dataset.y)

    header_u = header.copy()
    header_v = header.copy()

    header_u["parameterNumberName"] = "U-component_of_wind"
    header_u["parameterNumber"] = 2
    header_v["parameterNumberName"] = "V-component_of_wind"
    header_v["parameterNumber"] = 3

    for it, t in enumerate(time):
        if t >= time_range[0] and t <= time_range[1]:
            dd = []

            tstr = t.strftime("%Y-%m-%dT%H:%M:%SZ")

            u_list = (
                np.flipud(np.around(dataset.quantity[0].u[it, :, :], decimals=1))
                .flatten()
                .tolist()
            )
            data0 = {"header": header_u.copy(), "data": u_list}
            data0["header"]["refTime"] = tstr
            dd.append(data0)

            v_list = (
                np.flipud(np.around(dataset.quantity[0].v[it, :, :], decimals=1))
                .flatten()
                .tolist()
            )
            data0 = {"header": header_v.copy(), "data": v_list}
            data0["header"]["refTime"] = tstr
            dd.append(data0)

            data.append(dd)

    json_string = json.dumps(data, separators=(",", ":"))
    fid = open(file_name, "w")
    if js:
        fid.write("wind = ")
    fid.write(json_string)
    fid.close()
