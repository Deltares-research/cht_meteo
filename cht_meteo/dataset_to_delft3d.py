"""Export routines for writing MeteoDataset objects to Delft3D meteo format.

Supports both legacy ASCII (``.amu``/``.amv``/``.amp``/``.ampr``) and
NetCDF output.
"""

import os

import numpy as np
import pandas as pd
import xarray as xr


def write_to_delft3d_ascii(
    dataset,
    file_name: str,
    version: str = "1.03",
    path: str | None = None,
    header_comments: bool = False,
    refdate=None,
    parameters: list | None = None,
    time_range: list | None = None,
) -> None:
    """Write a MeteoDataset to Delft3D ASCII meteo files.

    One file is written per meteorological parameter using the standard
    Delft3D ``meteo_on_equidistant_grid`` format.

    Parameters
    ----------
    dataset : MeteoDataset
        Source dataset to export.
    file_name : str
        Base file name (without extension).
    version : str, optional
        Delft3D meteo file version string.  Default ``"1.03"``.
    path : str, optional
        Output directory.  If ``None`` the current directory is used.
    header_comments : bool, optional
        When ``True``, write verbose ``### ...`` comment lines in the header.
    refdate : datetime, optional
        Reference date for the time axis.  Defaults to the first time step.
    parameters : list of str, optional
        Parameters to export.  Defaults to
        ``["wind", "barometric_pressure", "precipitation"]``.
    time_range : list of datetime, optional
        ``[start, end]`` time window.  Defaults to the full dataset extent.
    """
    # Convert numpy.datetime64 to datetime
    time = dataset.ds.time.to_numpy().astype("M8[s]").astype("O")

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
            file["data"] = dataset.ds["wind_u"].to_numpy()[:]
            file["ext"] = "amu"
            file["quantity"] = "x_wind"
            file["unit"] = "m s-1"
            file["fmt"] = "%6.1f"
            files.append(file)
            file = {}
            file["data"] = dataset.ds["wind_v"].to_numpy()[:]
            file["ext"] = "amv"
            file["quantity"] = "y_wind"
            file["unit"] = "m s-1"
            file["fmt"] = "%6.1f"
            files.append(file)
        elif param == "barometric_pressure":
            file = {}
            file["data"] = dataset.ds["barometric_pressure"].to_numpy()[:]
            file["ext"] = "amp"
            file["quantity"] = "air_pressure"
            file["unit"] = "Pa"
            file["fmt"] = "%7.0f"
            files.append(file)
        elif param == "precipitation":
            file = {}
            file["data"] = dataset.ds["precipitation"].to_numpy()[:]
            file["ext"] = "ampr"
            file["quantity"] = "precipitation"
            file["unit"] = "mm h-1"
            file["fmt"] = "%7.1f"
            files.append(file)

    for file in files:
        if "lon" in dataset.ds:
            ncols = len(dataset.ds.lon)
            x = dataset.ds.lon.to_numpy()[:]
        else:
            ncols = len(dataset.ds.x)
            x = dataset.ds.x.to_numpy()[:]
        if "lat" in dataset.ds:
            nrows = len(dataset.ds.lat)
            y = dataset.ds.lat.to_numpy()[:]
        else:
            nrows = len(dataset.ds.y)
            y = dataset.ds.y.to_numpy()[:]

        dx = (x[-1] - x[0]) / (ncols - 1)
        dy = (y[-1] - y[0]) / (nrows - 1)

        if path:
            full_file_name = os.path.join(path, f"{file_name}.{file['ext']}")
        else:
            full_file_name = f"{file_name}.{file['ext']}"

        with open(full_file_name, "w") as fid:
            if header_comments:
                fid.write("### START OF HEADER\n")
                fid.write(
                    "### All text on a line behind the first # is parsed as commentary\n"
                )
                fid.write("### Additional comments\n")

            fid.write(
                f"FileVersion      =   {version}"
                "                                               # Version of meteo input file, to check if the newest file format is used\n"
            )
            fid.write(
                "filetype         =   meteo_on_equidistant_grid                          # Type of meteo input file: meteo_on_flow_grid, meteo_on_equidistant_grid, meteo_on_curvilinear_grid or meteo_on_spiderweb_grid\n"
            )
            fid.write(
                "NODATA_value     =   -999                                               # Value used for undefined or missing data\n"
            )
            fid.write(f"n_cols           =   {ncols}\n")
            fid.write(f"n_rows           =   {nrows}\n")
            fid.write(f"grid_unit        =   {grid_unit}\n")
            fid.write(f"x_llcorner       =   {min(x) - 0.5 * dx}\n")
            fid.write(f"y_llcorner       =   {min(y) - 0.5 * dy}\n")
            if version == "1.02":
                fid.write("value_pos       =    corner\n")
            # Use max 4 decimals for dx and dy
            if dataset.crs.is_geographic:
                fid.write(f"dx               =   {dx:.5f}\n")
                fid.write(f"dy               =   {dy:.5f}\n")
            else:
                fid.write(f"dx               =   {dx:.1f}\n")
                fid.write(f"dy               =   {dy:.1f}\n")
            fid.write(
                "n_quantity       =   1                                                  # Number of quantities prescribed in the file\n"
            )
            fid.write(f"quantity1        =   {file['quantity']}\n")
            fid.write(f"unit1            =   {file['unit']}\n")
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
                        f"TIME = {tim} minutes since "
                        f"{refdate.strftime('%Y-%m-%d %H:%M:%S')} +00:00\n"
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
                            f"Warning! Wind speed > 1000 m/s at {t.strftime('%Y-%m-%d %H:%M:%S')} !"
                        )
                    if np.min(val) < -1000.0:
                        val = np.zeros_like(
                            val
                        )  # Weird array, don't trust. Set everything to zeros.
                        print(
                            f"Warning! Wind speed > 1000 m/s at {t.strftime('%Y-%m-%d %H:%M:%S')} !"
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
                            f"Warning! Precipitation exceeds 1000 mm/h at {t.strftime('%Y-%m-%d %H:%M:%S')} !"
                        )
                        val[np.where(val == 0.0)] = np.nan
                    if np.nanmin(val) < 0.0:
                        val[np.where(val < 0.0)] = 0.0

                if np.all(np.isnan(val)):
                    if it > 0:
                        print(
                            f"Warning! Only NaNs found for {param} at "
                            f"{t.strftime('%Y-%m-%d %H:%M:%S')} ! Using data from previous time."
                        )
                        val = val0  # noqa: F821

                    else:
                        if (
                            file["quantity"] == "x_wind"
                            or file["quantity"] == "y_wind"
                            or file["quantity"] == "precipitation"
                        ):
                            print(
                                f"Warning! Only NaNs found for {param} at "
                                f"{t.strftime('%Y-%m-%d %H:%M:%S')} ! Setting values to 0.0 !"
                            )
                            val = np.zeros_like(val)
                        elif file["quantity"] == "air_pressure":
                            print(
                                f"Warning! Only NaNs found for {param} at "
                                f"{t.strftime('%Y-%m-%d %H:%M:%S')} ! Setting values to 101300.0 !"
                            )
                            val = np.zeros_like(val) + 101300.0

                fid.write(
                    f"TIME = {tim} minutes since "
                    f"{refdate.strftime('%Y-%m-%d %H:%M:%S')} +00:00\n"
                )
                np.savetxt(fid, val, fmt=file["fmt"])
                val0 = val  # Saved in case next time only has NaNs  # noqa: F841

            # Add extra blocks if data does not cover time range
            if time[-1] < time_range[1]:
                dt = time_range[1] - refdate
                tim = dt.total_seconds() / 60
                val = np.flipud(file["data"][-1, :, :])
                # Skip blocks with only nans
                if not np.all(np.isnan(val)):
                    val[val == np.nan] = -999.0
                    fid.write(
                        f"TIME = {tim} minutes since "
                        f"{refdate.strftime('%Y-%m-%d %H:%M:%S')} +00:00\n"
                    )
                    np.savetxt(fid, val, fmt=file["fmt"])


def write_to_delft3d_netcdf(
    dataset,
    file_name: str,
    path: str | None = None,
    refdate=None,
    parameters: list | None = None,
) -> None:
    """Write a MeteoDataset to Delft3D NetCDF meteo files.

    One NetCDF file is written per parameter group.

    Parameters
    ----------
    dataset : MeteoDataset
        Source dataset to export.
    file_name : str
        Base file name (without extension and without the parameter suffix).
    path : str, optional
        Output directory.  If ``None`` the current directory is used.
    refdate : datetime, optional
        Reference date for time encoding.  Defaults to the first time step.
    parameters : list of str, optional
        Parameters to export.  Defaults to
        ``["wind", "barometric_pressure", "precipitation"]``.
    """
    # Convert numpy.datetime64 to datetime
    time = dataset.ds.time.to_numpy().astype("M8[s]").astype("O")

    if not refdate:
        refdate = time[0]

    if not parameters:
        parameters = ["wind", "barometric_pressure", "precipitation"]

    files = []
    for param in parameters:
        if param == "wind":
            file = {}
            file["davars"] = ["wind_u", "wind_v"]
            file["ncvars"] = ["eastward_wind", "northward_wind"]
            file["ext"] = "_wind"
            file["unit"] = "m s-1"
            files.append(file)
        elif param == "barometric_pressure":
            file = {}
            file["davars"] = ["barometric_pressure"]
            file["ncvars"] = ["barometric_pressure"]
            file["ext"] = "_barometric_pressure"
            file["unit"] = "Pa"
            files.append(file)
        elif param == "precipitation":
            file = {}
            file["davars"] = ["precipitation"]
            file["ncvars"] = ["precipitation"]
            file["ext"] = "_precipitation"
            file["unit"] = "mm h-1"
            files.append(file)

    # Convert times to float minutes since 1970-01-01
    # subtract np.datetime64 from datetime.datetime object
    time = pd.to_datetime(time)
    float_minutes = (
        (time - np.datetime64("1970-01-01T00:00:00")) / pd.Timedelta(seconds=60)
    ).to_numpy()

    for file in files:
        if "lon" in dataset.ds:
            x = dataset.ds.lon.to_numpy()[:]
        else:
            x = dataset.ds.x.to_numpy()[:]
        if "lat" in dataset.ds:
            y = dataset.ds.lat.to_numpy()[:]
        else:
            y = dataset.ds.y.to_numpy()[:]

        if path:
            full_file_name = os.path.join(path, f"{file_name}{file['ext']}.nc")
        else:
            full_file_name = f"{file_name}{file['ext']}.nc"

        ds = xr.Dataset()

        ds["time"] = xr.DataArray(float_minutes, dims=["time"])
        ds["time"].attrs["units"] = "minutes since 1970-01-01 00:00:00.0 +0000"
        ds["x"] = xr.DataArray(x, dims=["x"])
        ds["y"] = xr.DataArray(y, dims=["y"])
        for davar, ncvar in zip(file["davars"], file["ncvars"]):
            if davar in dataset.ds:
                ds[ncvar] = xr.DataArray(
                    dataset.ds[davar].to_numpy()[:], dims=["time", "y", "x"]
                )
                ds[ncvar].attrs["units"] = file["unit"]
                ds[ncvar].attrs["long_name"] = davar.replace("_", " ").capitalize()
            else:
                print(f"Warning: {davar} not found in dataset. Skipping {ncvar}.")

        # Add attributes
        ds.attrs["description"] = "SFINCS meteo forcing"

        # Write netcdf file
        ds.to_netcdf(full_file_name, mode="w", format="NETCDF4", engine="netcdf4")
