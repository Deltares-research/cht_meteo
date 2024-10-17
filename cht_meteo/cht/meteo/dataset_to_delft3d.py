import os
import numpy as np

def write_to_delft3d_ascii(
    dataset,
    file_name,
    version="1.03",
    path=None,
    header_comments=False,
    refdate=None,
    parameters=None,
    time_range=None,
):
        
    # Convert numpy.datetime64 to datetime
    time = dataset.ds.time.values.astype("M8[s]").astype("O")

    if not refdate:
        refdate = time[0]

    if not time_range:
        time_range = [time[0], time[-1]]

    if not parameters:
        parameters = ["wind", "barometric_pressure", "precipitation"]

    if dataset.crs.is_geographic:
        grid_unit = "degrees"
    else:
        grid_unit = "m"

    files = []
    for param in parameters:
        if param == "wind":
            file = {}
            file["data"] = dataset.ds["wind_u"].values[:]
            file["ext"] = "amu"
            file["quantity"] = "x_wind"
            file["unit"] = "m s-1"
            file["fmt"] = "%6.1f"
            files.append(file)
            file = {}
            file["data"] = dataset.ds["wind_v"].values[:]
            file["ext"] = "amv"
            file["quantity"] = "y_wind"
            file["unit"] = "m s-1"
            file["fmt"] = "%6.1f"
            files.append(file)
        elif param == "barometric_pressure":
            file = {}
            file["data"] = dataset.ds["barometric_pressure"].values[:]
            file["ext"] = "amp"
            file["quantity"] = "air_pressure"
            file["unit"] = "Pa"
            file["fmt"] = "%7.0f"
            files.append(file)
        elif param == "precipitation":
            file = {}
            file["data"] = dataset.ds["precipitation"].values[:]
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

        if "lon" in dataset.ds:
            ncols = len(dataset.ds.lon)
            x = dataset.ds.lon.values[:]
        else:
            ncols = len(dataset.ds.x)
            x = dataset.ds.x.values[:]
        if "lat" in dataset.ds:
            nrows = len(dataset.ds.lat)
            y = dataset.ds.lat.values[:]
        else:
            nrows = len(dataset.ds.y)
            y = dataset.ds.y.values[:]

        dx = (x[-1] - x[0]) / (ncols - 1)
        dy = (y[-1] - y[0]) / (nrows - 1)

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
        fid.write("x_llcorner       =   " + str(min(x) - 0.5 * dx) + "\n")
        fid.write("y_llcorner       =   " + str(min(y) - 0.5 * dy) + "\n")
        if version == "1.02":
            fid.write("value_pos       =    corner\n")
        # Use max 4 decimals for dx and dy
        if dataset.crs.is_geographic:
            fid.write("dx               =   " + "{:.5f}".format(dx) + "\n")
            fid.write("dy               =   " + "{:.5f}".format(dy) + "\n")
        else:
            fid.write("dx               =   " + "{:.1f}".format(dx) + "\n")
            fid.write("dy               =   " + "{:.1f}".format(dy) + "\n")
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

        for it, t in enumerate(time):

            dt = t - refdate
            tim = dt.total_seconds() / 60
            val = np.flipud(file["data"][it, :, :])

            if file["quantity"] == "x_wind" or file["quantity"] == "y_wind":
                if np.max(val) > 1000.0:
                    val = np.zeros_like(
                        val
                    )  # Weird array, don't trust. Set everything to zeros.
                    val[np.where(val == 0.0)] = np.nan
                    print(
                        "Warning! Wind speed > 1000 m/s at "
                        + t.strftime("%Y-%m-%d %H:%M:%S")
                        + " !"
                    )
                if np.min(val) < -1000.0:
                    val = np.zeros_like(
                        val
                    )  # Weird array, don't trust. Set everything to zeros.
                    print(
                        "Warning! Wind speed > 1000 m/s at "
                        + t.strftime("%Y-%m-%d %H:%M:%S")
                        + " !"
                    )
                    val[np.where(val == 0.0)] = np.nan
            if file["quantity"] == "air_pressure":        
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
            if file["quantity"] == "precipitation":
                if np.nanmax(val) > 1000.0:
                    val = np.zeros_like(
                        val
                    )  # Weird array, don't trust. Set everything to zeros.
                    print(
                        "Warning! Precipitation exceeds 1000 mm/h at "
                        + t.strftime("%Y-%m-%d %H:%M:%S")
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
                        + t.strftime("%Y-%m-%d %H:%M:%S")
                        + " ! Using data from previous time."
                    )
                    val = val0

                else:
                    if file["quantity"] == "x_wind" or file["quantity"] == "y_wind" or file["quantity"] == "precipitation":
                        print(
                            "Warning! Only NaNs found for "
                            + param
                            + " at "
                            + t.strftime("%Y-%m-%d %H:%M:%S")
                            + " ! Setting values to 0.0 !"
                        )
                        val = np.zeros_like(val)
                    elif file["quantity"] == "air_pressure":
                        print(
                            "Warning! Only NaNs found for "
                            + param
                            + " at "
                            + t.strftime("%Y-%m-%d %H:%M:%S")
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
            val0 = val # Saved in case next time only has NaNs

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

