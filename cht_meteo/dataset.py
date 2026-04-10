"""Core MeteoDataset class and supporting helpers for meteorological data handling.

Provides storage, download orchestration, spatial cut-out / interpolation,
and export routines (Delft3D, JSON wind) for gridded meteorological datasets.
"""

import copy
import datetime
import os
from typing import Optional

import cht_utils.fileops as fo
import numpy as np
import pandas as pd
import xarray as xr
from pyproj import CRS, Transformer
from scipy.interpolate import RegularGridInterpolator

from cht_meteo.dataset_to_delft3d import write_to_delft3d_ascii, write_to_delft3d_netcdf
from cht_meteo.dataset_to_json_wind import write_wind_to_json


class MeteoDataset:
    """Container for a gridded meteorological dataset.

    Parameters
    ----------
    **kwargs
        Arbitrary keyword arguments used to set instance attributes.
        Special keys ``time``, ``x``/``lon``, and ``y``/``lat`` initialise
        the underlying xarray Dataset via :meth:`init_ds`.
    """

    def __init__(self, **kwargs) -> None:
        # Default values set to None
        self.name = None  # Name of the dataset
        self.path = None  # Path to the dataset
        self.long_name = None  # Long name of the dataset (not really needed)
        self.source = None  # Source of the dataset (currently one of "gfs_forecast_0p25", "gfs_analysis_0p50", "coamps_tc_forecast" or None)
        self.parameters = (
            None  # Parameters in the dataset (a list with strings). Not currently used.
        )
        self.lon_range = None  # Longitude range of the dataset
        self.lat_range = None  # Latitude range of the dataset
        self.var_names = [
            "wind_u",
            "wind_v",
            "barometric_pressure",
            "precipitation",
        ]  # Variable names in the dataset
        self.crs = CRS(4326)  # Coordinate reference system of the dataset
        self.tau = 0  # Time interval in hours between cycle and data
        self.last_analysis_time = None  # Time of last analysis in the dataset
        self.last_forecast_cycle_time = (
            None  # Time of last forecast cycle in the dataset
        )

        # Loop through kwargs to set attributes
        reserved_keys = ["x", "y", "lon", "lat", "time"]
        for key, value in kwargs.items():
            if key not in reserved_keys:
                setattr(self, key, value)
        self.resolution = 9999.0
        self.lazy = False

        # Set some source information
        self.source_name = ""  # Name of the source
        self.source_type = "forecast"  # Can be "forecast" or "analysis"
        self.source_delay = 0  # Delay in hours before data is available
        self.source_cycle_interval = 6  # Interval in hours between cycles
        self.source_time_interval = 3  # Output interval in hours
        self.source_forecast_duration = (
            240  # The length of the forecast in hours that will be downloaded
        )

        # Loop through kwargs to set attributes
        reserved_keys = ["x", "y", "lon", "lat", "time"]
        for key, value in kwargs.items():
            if key not in reserved_keys:
                setattr(self, key, value)

        # Empty list for subsets (only for COAMPS-TC for now)
        self.subset = []

        # Create empty xarray dataset
        self.ds = xr.Dataset()

        time = None
        x = None
        y = None

        if "time" in kwargs:
            time = kwargs["time"]
        if "x" in kwargs:
            x = kwargs["x"]
        if "lon" in kwargs:
            x = kwargs["lon"]
        if "y" in kwargs:
            y = kwargs["y"]
        if "lat" in kwargs:
            y = kwargs["lat"]

        if time is not None and x is not None and y is not None:
            self.init_ds(time, x, y)

    def init_ds(self, time: np.ndarray, x: np.ndarray, y: np.ndarray) -> None:
        """Initialise the xarray Dataset with NaN-filled arrays.

        Parameters
        ----------
        time : array-like
            Time coordinate values.
        x : array-like
            Longitude (or x) coordinate values.
        y : array-like
            Latitude (or y) coordinate values.
        """

        x = np.array(x)
        y = np.array(y)
        time = np.array(time)

        # Create empty dataset
        self.ds["time"] = xr.DataArray(time, dims=("time"))
        if self.crs.is_geographic:
            self.ds["lon"] = xr.DataArray(x, dims=("lon"))
            self.ds["lat"] = xr.DataArray(y, dims=("lat"))
        else:
            self.ds["x"] = xr.DataArray(x, dims=("x"))
            self.ds["y"] = xr.DataArray(y, dims=("y"))
        # empty_data = np.empty((len(time), len(y), len(x))) + np.nan
        for var_name in self.var_names:
            if self.crs.is_geographic:
                self.ds[var_name] = xr.DataArray(
                    np.empty((len(time), len(y), len(x))) + np.nan,
                    dims=("time", "lat", "lon"),
                )
            else:
                self.ds[var_name] = xr.DataArray(
                    np.empty((len(time), len(y), len(x))) + np.nan,
                    dims=("time", "y", "x"),
                )

    def fill(self, u: float = 0.0, v: float = 0.0, p: float = 101300.0, pr: float = 0.0) -> None:
        """Fill the dataset with constant values.

        Parameters
        ----------
        u : float, optional
            East-ward wind component (m/s).  Default ``0.0``.
        v : float, optional
            North-ward wind component (m/s).  Default ``0.0``.
        p : float, optional
            Barometric pressure (Pa).  Default ``101300.0``.
        pr : float, optional
            Precipitation rate (mm/h).  Default ``0.0``.
        """

        if "wind_u" in self.ds:
            self.ds["wind_u"][:, :, :] = u
        if "wind_v" in self.ds:
            self.ds["wind_v"][:, :, :] = v
        if "barometric_pressure" in self.ds:
            self.ds["barometric_pressure"][:, :, :] = p
        if "precipitation" in self.ds:
            self.ds["precipitation"][:, :, :] = pr

    def download(self, time_range: list, **kwargs) -> None:
        """Download data for the given time range.

        Parameters
        ----------
        time_range : list of datetime
            ``[start, end]`` time interval to download.
        **kwargs
            Forwarded to :meth:`download_forecast` or
            :meth:`download_analysis`.
        """
        # Make path
        os.makedirs(self.path, exist_ok=True)

        if self.source_type == "forecast":
            # Download forecast from cycles
            self.download_forecast(time_range, **kwargs)
        else:
            # Download analysis
            self.download_analysis(time_range, **kwargs)

    def download_forecast(self, time_range: list, **kwargs) -> None:
        """Download all forecast cycles that cover the requested time range.

        Parameters
        ----------
        time_range : list of datetime
            ``[start, end]`` time interval.
        **kwargs
            Optional key ``last_cycle`` (datetime) caps the latest cycle
            that will be downloaded.
        """

        # Check if this dataset has a download_forecast_cycle method
        if not hasattr(self, "download_forecast_cycle"):
            print(
                f"Error: download_forecast_cycle method not implemented for dataset {self.name}"
            )
            return

        last_cycle_time = None
        if "last_cycle" in kwargs:
            if kwargs["last_cycle"] is not None:
                last_cycle_time = kwargs["last_cycle"]

        # Round first time in range down to hour
        h0 = time_range[0].hour
        # Round down to cycle interval
        h0 = h0 - np.mod(h0, self.source_cycle_interval)
        t0 = time_range[0].replace(
            microsecond=0, second=0, minute=0, hour=h0, tzinfo=datetime.timezone.utc
        )

        # Round last time in range down to hour
        h1 = time_range[1].hour
        # Round down to cycle interval
        h1 = h1 - np.mod(h1, self.source_cycle_interval)
        t1 = time_range[1].replace(
            microsecond=0, second=0, minute=0, hour=h1, tzinfo=datetime.timezone.utc
        )

        print(
            f"Now : {datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d_%H%M%S')}"
        )

        # Determine latest available cycle
        t_latest = datetime.datetime.now(datetime.timezone.utc) - datetime.timedelta(
            hours=self.source_delay
        )
        h0 = t_latest.hour
        h0 = h0 - np.mod(h0, self.source_cycle_interval)
        t_latest = t_latest.replace(microsecond=0, second=0, minute=0, hour=h0)
        print(f"t_latest : {t_latest.strftime('%Y%m%d_%H%M%S')}")

        # Make sure t0 and t1 are not later than t_latest
        t0 = min(t0, t_latest)
        t1 = min(t1, t_latest)

        print(
            f"Downloading from {self.name} - cycles : "
            f"{t0.strftime('%Y%m%d_%Hz')} to {t1.strftime('%Y%m%d_%Hz')}"
        )

        # Make list with cycle times
        cycle_times = (
            pd.date_range(start=t0, end=t1, freq=f"{self.source_cycle_interval}h")
            .to_pydatetime()
            .tolist()
        )

        # If last_cycle has been provided, do not get data from later cycles
        if last_cycle_time is not None:
            print(
                f"Downloading up to last_cycle : {last_cycle_time.strftime('%Y%m%d_%Hz')}"
            )
            cycle_times = [t for t in cycle_times if t <= last_cycle_time]

        # Loop through all cycle times
        for it, t in enumerate(cycle_times):
            print(
                f"Downloading {it + 1} of {len(cycle_times)} - cycle : "
                f"{t.strftime('%Y%m%d_%Hz')}"
            )
            try:
                if it < len(cycle_times) - 1:
                    # Earlier cycle, so only get data up to next cycle
                    self.download_forecast_cycle(
                        cycle_time=t,
                        time_range=[cycle_times[it], cycle_times[it + 1]],
                        **kwargs,
                    )
                else:
                    self.download_forecast_cycle(cycle_time=t, **kwargs)
            except Exception as e:
                print(
                    f"Error downloading data from dataset {self.name} - cycle : {t.strftime('%Y%m%d_%Hz')}"
                )
                print(e)

    def download_analysis(self, time_range: list, **kwargs) -> None:
        """Download analysis or reanalysis data for the requested time range.

        Parameters
        ----------
        time_range : list of datetime
            ``[start, end]`` time interval.
        **kwargs
            Forwarded to the concrete download method.
        """

        # Check if this dataset has a download_analysis_times method (e.g. gfs_analysis_0p50)

        if hasattr(self, "download_analysis_times"):
            # Make list of requested times
            freqstr = f"{self.source_time_interval}h"
            requested_times = (
                pd.date_range(start=time_range[0], end=time_range[1], freq=freqstr)
                .to_pydatetime()
                .tolist()
            )

            # Loop through all requested times and see which data is already available
            rtimes = []
            # Check which files do not yet exist
            for t in requested_times:
                file_name = os.path.join(
                    self.path, f"{self.name}.{t.strftime('%Y%m%d_%H%M')}.nc"
                )
                if not os.path.exists(file_name):
                    rtimes.append(t)

            if rtimes:
                try:
                    self.download_analysis_times(rtimes)
                except Exception as e:
                    print(f"Error downloading data from dataset {self.name}")
                    print(e)
            else:
                print("Requested meteo data already available")

        elif hasattr(self, "download_reanalysis"):  # e.g. era5_reanalysis_0p25
            try:
                self.download_reanalysis(time_range, **kwargs)
            except Exception as e:
                print(f"Error downloading data from dataset {self.name}")
                print(e)

    def collect(self, time_range: list, **kwargs) -> xr.Dataset:
        """Merge per-cycle netCDF files into a single in-memory xarray Dataset.

        The actual data is stored in ``self.ds`` (or in each element of
        ``self.subset`` for multi-resolution datasets such as COAMPS-TC).

        Parameters
        ----------
        time_range : list of datetime
            ``[start, end]`` time interval to collect.
        **kwargs
            Optional keys: ``tau`` (int, hours), ``last_cycle`` (datetime).

        Returns
        -------
        xr.Dataset
            The assembled dataset (also stored on ``self.ds``).
        """

        if "tau" in kwargs:
            tau = kwargs["tau"]
        else:
            tau = self.tau

        last_cycle_time = None
        if "last_cycle" in kwargs:
            if kwargs["last_cycle"] is not None:
                # last_cycle_time = datetime.datetime.strptime(kwargs["last_cycle"], "%Y%m%d_%H")
                last_cycle_time = kwargs["last_cycle"]
                # FIXME make timezone naive
                last_cycle_time = last_cycle_time.replace(tzinfo=None)

        # Subsets are only used when there are subsets with different resolutions (e.g. as in COAMPS-TC)
        if len(self.subset) > 0:
            # There are subsets in this dataset. For now we automatically collect data of each subset.
            subsets = True
            subsets_to_get = []
            for subset in self.subset:
                subsets_to_get.append(subset.name)

        else:
            subsets = False
            subsets_to_get = ["dummy"]

        # Loop through all subsets (i.e. nesting levels)
        for isub, subset_name in enumerate(subsets_to_get):
            # Make empty lists for time and files
            time_list = []
            file_list = []

            if subsets:
                subsetstr = "." + self.subset[isub].name + "."
                moving = self.subset[isub].moving
            else:
                subsetstr = ""
                moving = False

            if self.source_type == "forecast":
                # Make list of all cycles
                all_cycle_paths = fo.list_folders(os.path.join(self.path, "*"))

                icycle = -1

                # Loop through all cycle paths
                for cycle_path in all_cycle_paths:
                    # Get time from cycle path
                    t_cycle = datetime.datetime.strptime(
                        cycle_path[-12:-1], "%Y%m%d_%H"
                    )

                    # If last_cycle has been provided, do not get data from later cycles
                    if last_cycle_time is not None:
                        if t_cycle > last_cycle_time:
                            # We should stop now
                            break

                    # Check if time of this cycle falls within requested range
                    if t_cycle < time_range[0] or t_cycle > time_range[1]:
                        # We can skip this cycle
                        continue

                    # Set this as the last cycle that is collected (used e.g. to display in CoSMoS webviewer)
                    self.last_forecast_cycle_time = t_cycle
                    self.last_analysis_time = t_cycle

                    # Find all times available in this cycle as it may contain our data
                    files_in_cycle = fo.list_files(
                        os.path.join(cycle_path, f"*{subsetstr}*.nc")
                    )

                    icycle += 1

                    # Loop through all files in this cycle
                    for ifile, file in enumerate(files_in_cycle):
                        t_file = datetime.datetime.strptime(file[-16:-3], "%Y%m%d_%H%M")

                        if ifile == 0:
                            self.last_analysis_time = t_file

                        if tau > 0 and icycle > 0:
                            # Compute time interval between cycle and file
                            th = int((t_file - t_cycle).total_seconds() / 3600)
                            if th < tau:
                                # We can skip this file
                                continue

                        if t_file < time_range[0] or t_file > time_range[1]:
                            # We can skip this file
                            continue

                        # We want the data in this file !

                        # Check if time is already available in time_list. If so, insert it in the correct place. If not, append it.
                        if t_file in time_list:
                            # Find index
                            ind = time_list.index(t_file)
                            time_list[ind] = t_file
                            file_list[ind] = file
                        else:
                            time_list.append(t_file)
                            file_list.append(file)

            else:
                # A lot easier
                files_in_cycle = fo.list_files(
                    os.path.join(self.path, f"*{subsetstr}*.nc")
                )
                for file in files_in_cycle:
                    try:
                        t_file = datetime.datetime.strptime(file[-15:-3], "%Y%m%d%H%M")
                    except Exception:
                        t_file = datetime.datetime.strptime(file[-16:-3], "%Y%m%d_%H%M")
                    if t_file >= time_range[0] and t_file <= time_range[1]:
                        file_list.append(os.path.join(self.path, file))
                        time_list.append(t_file)

            # Now we loop through the files, read them and store them in large lazy array
            if not time_list:
                print("No meteo data files found within requested time range")
                return

            # Create datasets with injected time coordinate
            datasets = [
                add_time_coord(
                    xr.open_dataset(f, chunks={"lat": 512, "lon": 512}), t, moving
                )
                for f, t in zip(file_list, time_list)
            ]

            ds = xr.concat(datasets, dim="time")

            # add physically logical fill values to the variables
            # NOTE is this needed here, also happens when writing it seems?
            fill_values = {
                "wind_u": 0.0,
                "wind_v": 0.0,
                "precipitation": 0.0,
                "barometric_pressure": 101300.0,
            }

            for var in self.var_names:
                fill_value = fill_values.get(var, 0.0)  # Default fallback
                if var in ds:
                    ds[var] = ds[var].where(~np.isnan(ds[var]), other=fill_value)

            # Make sure files are oriented S-N
            lat = ds["lat"]
            if lat.ndim == 1:
                if lat[0] > lat[-1]:
                    ds = ds.isel(lat=slice(None, None, -1))
            elif lat.ndim == 2:
                if np.nanmean(lat[0]) > np.nanmean(lat[-1]):
                    ds = ds.isel(lat=slice(None, None, -1))

            # Only load the dataset if lazy is False (in-memory is typically faster)
            if not self.lazy:
                ds = ds.load()

            if subsets:
                # Store the data in the subset
                self.subset[isub].ds = ds
            else:
                self.ds = ds

        return ds

    def cut_out(
        self,
        name: str | None = None,
        lon_range: list | None = None,
        lat_range: list | None = None,
        time_range: list | None = None,
        x_range: list | None = None,
        y_range: list | None = None,
        lon: np.ndarray | None = None,
        lat: np.ndarray | None = None,
        x: np.ndarray | None = None,
        y: np.ndarray | None = None,
        dx: float | None = None,
        dy: float | None = None,
        copy_time: bool = False,
        crs: CRS | None = None,
    ) -> "MeteoDataset":
        """Return a new dataset cut out and optionally re-gridded from this one.

        If the original dataset has subsets, they are first interpolated onto
        the coarsest subset before the cut-out is applied.

        Parameters
        ----------
        name : str, optional
            Name for the new dataset.
        lon_range : list, optional
            ``[west, east]`` longitude bounds.
        lat_range : list, optional
            ``[south, north]`` latitude bounds.
        time_range : list of datetime, optional
            ``[start, end]`` time bounds.
        x_range : list, optional
            Alternative to ``lon_range`` for projected CRS.
        y_range : list, optional
            Alternative to ``lat_range`` for projected CRS.
        lon : array-like, optional
            Explicit longitude array for the output grid.
        lat : array-like, optional
            Explicit latitude array for the output grid.
        x : array-like, optional
            Explicit x-coordinate array for the output grid.
        y : array-like, optional
            Explicit y-coordinate array for the output grid.
        dx : float, optional
            Grid spacing in x; triggers interpolation when combined with
            ``x_range`` / ``y_range``.
        dy : float, optional
            Grid spacing in y.
        copy_time : bool, optional
            Copy the time coordinate from the source dataset.
        crs : CRS, optional
            Target coordinate reference system.  Defaults to ``self.crs``.

        Returns
        -------
        MeteoDataset
            New dataset containing the cut-out data.
        """

        # Start with making a copy of self

        if crs is None:
            crs = self.crs

        # Check the options (call everything x and y for now)
        if lon_range is not None:
            x_range = lon_range
        if lat_range is not None:
            y_range = lat_range
        if lon is not None:
            x = lon
        if lat is not None:
            y = lat

        # There are now 3 options for the cut-out:
        # 1) none are given (no cut-out, but still interpolation), so set x_range and y_range and continue to option 2
        # 2) x_range and y_range are given, but not dx and dy (just cut-out, no interpolation) - this will only work if the CRS is the same for both datasets
        # 3) x_range and y_range are given, and dx and dy are given (cut-out and interpolate)
        # 4) x and y are given, but not x_range and y_range (cut-out and interpolate)

        # Option 1
        if (
            x_range is None
            and y_range is None
            and dx is None
            and dy is None
            and x is None
            and y is None
        ):
            # No cut-out, just interpolation
            x_range = [-1.0e9, 1.0e9]
            y_range = [-1.0e9, 1.0e9]

        # Option 2
        if x_range is not None and y_range is not None and dx is None and dy is None:
            # Cut-out only (but possible interpolation)
            if crs != self.crs:
                print(
                    "Error: cut-out with different CRS requires arrays x and y, or dx and dy in combination with x_range and y_range"
                )
                return
            # Get coordinates (in case of subsets, we take the largest (coarsest) subset)
            x, y = self.get_coordinates(0)
            # Let's see if x and y are regularly spaced
            if np.min(np.diff(x)) == np.max(np.diff(x)) and np.min(
                np.diff(y)
            ) == np.max(np.diff(y)):
                # Regular spacing, all good
                pass
            else:
                # Irregular spacing
                # Create new regular x and y arrays with similar spacing
                dx = np.mean(np.diff(x))
                dy = np.mean(np.diff(y))
                x = np.arange(x[0], x[-1], dx)
                y = np.arange(y[0], y[-1], dy)
                # Reset to None
                dx, dy = None, None
            # Limit x, y
            ix0 = np.where(x <= x_range[0])[0]
            if len(ix0) > 0:
                ix0 = ix0[-1]
            else:
                ix0 = 0

            ix1 = np.where(x >= x_range[1])[0]
            if len(ix1) > 0:
                ix1 = ix1[0]
            else:
                ix1 = len(x)

            it0 = np.where(y <= y_range[0])[0]
            if len(it0) > 0:
                it0 = it0[-1]
            else:
                it0 = 0

            it1 = np.where(y >= y_range[1])[0]
            if len(it1) > 0:
                it1 = it1[0]
            else:
                it1 = len(y)

            if ix1 <= ix0 or it1 <= it0:
                print("Error: cut-out is empty with given x_range and y_range")
                return

            x = x[ix0:ix1]
            y = y[it0:it1]

        # Option 3
        elif (
            x_range is not None
            and y_range is not None
            and dx is not None
            and dy is not None
        ):
            # Create new x and y arrays
            x = np.arange(x_range[0] - dx, x_range[1] + dx, dx)
            y = np.arange(y_range[0] - dy, y_range[1] + dy, dy)

        # Option 4
        elif x is not None and y is not None:
            # No need to do anything
            pass

        if len(self.subset) > 0:
            t = self.subset[0].ds["time"].to_numpy()
        else:
            t = self.ds["time"].to_numpy()

        if time_range is not None:
            # convert values in time_range to np.datetime64
            time_range = [np.datetime64(t) for t in time_range]
            t = t[(t >= time_range[0]) & (t <= time_range[1])]

        # Create new dataset
        if len(self.subset) > 0 or crs != self.crs:
            dataset = MeteoDataset(name=name, x=x, y=y, time=t, crs=crs)
            # make deepcopy of self to avoid modifying the original dataset and clip to new extent
            self_clipped = copy.deepcopy(self)
            # Transformer from new dataset to original dataset CRS
            transformer = Transformer.from_crs(dataset.crs, self.crs, always_xy=True)
            xg, yg = np.meshgrid(x, y)
            xg_t, yg_t = transformer.transform(xg, yg)
            if len(self.subset) > 0:
                # Clip all datasets in the subsets
                for isub in range(len(self.subset)):
                    lon_slice = get_buffered_slice(
                        self.subset[isub].ds.lon.values, xg_t.min(), xg_t.max()
                    )
                    lat_slice = get_buffered_slice(
                        self.subset[isub].ds.lat.values, yg_t.min(), yg_t.max()
                    )
                    self_clipped.subset[isub].ds = (
                        self_clipped.subset[isub]
                        .ds.sel(
                            lon=lon_slice,
                            lat=lat_slice,
                            time=slice(t[0], t[-1]),
                        )
                        .load()
                    )
            else:
                # Clip the main dataset
                lon_slice = get_buffered_slice(
                    self_clipped.ds.lon.values, xg_t.min(), xg_t.max()
                )
                lat_slice = get_buffered_slice(
                    self_clipped.ds.lat.values, yg_t.min(), yg_t.max()
                )
                self_clipped.ds = self_clipped.ds.sel(
                    lon=lon_slice,
                    lat=lat_slice,
                    time=slice(t[0], t[-1]),
                ).load()
            dataset.interpolate_dataset(self_clipped)
        else:
            dataset = copy.deepcopy(self)
            lon_slice = get_buffered_slice(dataset.ds.lon.values, x[0], x[-1])
            lat_slice = get_buffered_slice(dataset.ds.lat.values, y[0], y[-1])
            dataset.ds = dataset.ds.sel(
                lon=lon_slice, lat=lat_slice, time=slice(t[0], t[-1])
            )
            # Interpolate to new grid when dx and dy were initially provided
            if dx is not None and dy is not None:
                dataset.ds = dataset.ds.interp(lat=y, method="linear")
                dataset.ds = dataset.ds.interp(lon=x, method="linear")

            # load the dataset for faster writing
            dataset.ds.load()

            # Fill nans
            fill_values = {
                "wind_u": 0.0,
                "wind_v": 0.0,
                "precipitation": 0.0,
                "barometric_pressure": 101300.0,
            }

            # Loop over the dictionary and apply fill values where needed
            for var, fill_value in fill_values.items():
                if var in dataset.ds:
                    dataset.ds[var] = dataset.ds[var].where(
                        ~np.isnan(dataset.ds[var]), other=fill_value
                    )

        return dataset

    def interpolate_dataset(self, dataset: "MeteoDataset", copy_time: bool = False) -> None:
        """Interpolate data from another dataset onto this dataset's grid.

        Parameters
        ----------
        dataset : MeteoDataset
            Source dataset to interpolate from.
        copy_time : bool, optional
            If ``True``, also copy the time coordinate from *dataset*.
        """
        # Dimensions have already been set (also time?)
        if copy_time:
            self.ds["time"] = dataset.ds["time"]
        # Get horizontal coordinates
        if "x" in self.ds:
            x = self.ds["x"].to_numpy()
            y = self.ds["y"].to_numpy()
        else:
            x = self.ds["lon"].to_numpy()
            y = self.ds["lat"].to_numpy()
        xg, yg = np.meshgrid(x, y)
        # Loop through variables
        for var_name in self.var_names:
            da = self.ds[var_name].copy()
            if var_name in dataset.ds:
                for it, t in enumerate(self.ds["time"].to_numpy()):
                    # Get data
                    v = dataset.interpolate_variable(var_name, t, xg, yg, crs=self.crs)
                    # Set points in v equal to points in original data vori where vori already has a value
                    vori = da.loc[dict(time=t)].to_numpy()[:]
                    not_nan = np.where(~np.isnan(vori))
                    v[not_nan] = vori[not_nan]
                    da.loc[dict(time=t)] = v

                self.ds[var_name] = da

    def merge_datasets(self, datasets: list, **kwargs) -> None:
        """Merge a list of datasets into this one by interpolation.

        Parameters
        ----------
        datasets : list of MeteoDataset
            Datasets to merge, typically at different resolutions.
        **kwargs
            Forwarded to :meth:`interpolate_dataset`.
        """
        for dataset in datasets:
            self.interpolate_dataset(dataset)

    def interpolate_variable(
        self,
        var_name: str,
        t: datetime.datetime,
        x: np.ndarray,
        y: np.ndarray,
        crs: CRS | None = None,
        fill_missing_data: bool = False,
    ) -> np.ndarray:
        """Interpolate a single variable to arbitrary coordinates at a given time.

        Parameters
        ----------
        var_name : str
            Variable name present in ``self.ds``.
        t : datetime or np.datetime64
            Target time.
        x : np.ndarray
            Target x (or longitude) coordinates; 1-D or 2-D.
        y : np.ndarray
            Target y (or latitude) coordinates; 1-D or 2-D.
        crs : CRS, optional
            CRS of the *x*/*y* coordinates.  If different from ``self.crs``
            the points are re-projected before interpolation.
        fill_missing_data : bool, optional
            When ``True``, fill out-of-domain points with a physical default
            instead of NaN.

        Returns
        -------
        np.ndarray
            Interpolated values with the same shape as *x*/*y*.
        """

        # Check shape of x and y. They can be either:
        # 1) 1D arrays with regular spacing (turn into grid)
        # 2) 1D arrays with irregular spacing (treat them as point cloud)
        # 3) 2D arrays with regular or irregular spacing
        if len(np.shape(x)) == 1 and len(np.shape(y)) == 1:
            # 1D arrays
            # Check if x and y are regularly spaced
            if np.min(np.diff(x)) == np.max(np.diff(x)) and np.min(
                np.diff(y)
            ) == np.max(np.diff(y)):
                # 1D arrays with regular spacing
                xg, yg = np.meshgrid(x, y)
            else:
                # 1D arrays with irregular spacing
                xg = x
                yg = y
        else:
            # 2D arrays
            xg = x
            yg = y

        if crs is not None:
            # Transformer for x, y to lon, lat
            transformer = Transformer.from_crs(crs, self.crs, always_xy=True)
            xg, yg = transformer.transform(xg, yg)

        # Check if there are subsets
        if len(self.subset) > 0:
            # There are subsets in this dataset. For now we automatically collect data of each subset.
            subsets = True
            nsub = len(self.subset)
        else:
            subsets = False
            nsub = 1

        # Make array of nans
        v = np.empty(np.shape(xg))
        v[:] = np.nan

        # Loop in reverse order to get the highest resolution data first
        for isub in range(nsub - 1, -1, -1):
            if subsets:
                ds = self.subset[isub].ds
            else:
                ds = self.ds

            # Check if t is a numpy.datetime64
            if not isinstance(t, np.datetime64):
                t = np.datetime64(t)

            # Get data
            if t in ds["time"].to_numpy()[:]:
                # Get data at time t
                da = ds[var_name].sel(time=t)
            else:
                # Interpolate data at time t
                da = ds[var_name].interp(time=t)

            # Get horizontal coordinates
            if self.crs.is_geographic:
                x = da["lon"].to_numpy()
                y = da["lat"].to_numpy()
            else:
                x = da["x"].to_numpy()
                y = da["y"].to_numpy()

            # Make interpolator
            interp = RegularGridInterpolator((y, x), da.to_numpy()[:])

            # Find points outside of grid
            iout = np.where(
                (xg < np.min(x))
                | (xg > np.max(x))
                | (yg < np.min(y))
                | (yg > np.max(y))
            )
            x1 = np.maximum(xg, np.min(x))
            x1 = np.minimum(x1, np.max(x))
            y1 = np.maximum(yg, np.min(y))
            y1 = np.minimum(y1, np.max(y))

            # Find points outside of grid
            v1 = interp((y1, x1))

            # Set values outside of grid to no_data_value
            if fill_missing_data:
                if var_name == "barometric_pressure":
                    no_data_value = 101300.0
                else:
                    no_data_value = 0.0
                v1[iout] = no_data_value
            else:
                v1[iout] = np.nan

            # Where v is nan, replace with v1
            v[np.isnan(v)] = v1[np.isnan(v)]

        return v

    def to_netcdf(self, file_name: Optional[os.PathLike] = None, **kwargs) -> None:
        """Write the dataset to netCDF file(s).

        Parameters
        ----------
        file_name : path-like, optional
            If given, write the entire dataset to a single file; otherwise
            write one file per time step to ``self.path``.
        **kwargs
            Currently unused; reserved for future use.
        """
        if file_name:
            # Write to single file
            self.ds.to_netcdf(path=file_name)
        else:
            # Write to database
            os.makedirs(self.path, exist_ok=True)
            # Loop through times in ds
            time = self.ds["time"].to_numpy()
            for it, t in enumerate(time):
                time_string = pd.to_datetime(t).strftime("%Y%m%d_%H%M")
                file_name = f"{self.name}.{time_string}.nc"
                full_file_name = os.path.join(self.path, file_name)
                ds_time = self.ds.isel(time=it)
                # Remove time and reftime
                ds_time = ds_time.drop_vars(["time", "reftime"])
                ds_time.to_netcdf(path=full_file_name)
                ds_time.close()

    def to_delft3d(
        self,
        file_name: Optional[os.PathLike] = None,
        version: str = "1.03",
        path: str | None = None,
        header_comments: bool = False,
        refdate: datetime.datetime | None = None,
        parameters: list | None = None,
        time_range: list | None = None,
        format: str = "ascii",
    ) -> None:
        """Export the dataset in Delft3D meteo format.

        Parameters
        ----------
        file_name : path-like, optional
            Base file name (without extension).
        version : str, optional
            Delft3D meteo file version string.  Default ``"1.03"``.
        path : str, optional
            Output directory.
        header_comments : bool, optional
            Write verbose header comment lines.
        refdate : datetime, optional
            Reference date for the time axis.
        parameters : list of str, optional
            Parameters to export (e.g. ``["wind", "barometric_pressure"]``).
        time_range : list of datetime, optional
            ``[start, end]`` subset for output.
        format : str, optional
            ``"ascii"`` (default) or ``"netcdf"``.
        """
        if len(self.subset) > 0:
            # There are subsets in this dataset. For now we automatically collect data of each subset.
            print(
                "Warning! Cannot write meteo dataset with subsets to delft3d format ! Please consider merging the subsets first."
            )
            return

        if format == "ascii":
            write_to_delft3d_ascii(
                self,
                str(file_name),
                version,
                path,
                header_comments,
                refdate,
                parameters,
                time_range,
            )
        elif format == "netcdf":
            write_to_delft3d_netcdf(
                self, file_name, path=path, refdate=refdate, parameters=parameters
            )

    def wind_to_json(self, file_name: str, time_range: list | None = None, js: bool = True, iref: int = 1) -> None:
        """Write wind data to a JSON file suitable for web-based wind visualisation.

        Parameters
        ----------
        file_name : str
            Output file path.
        time_range : list of datetime, optional
            ``[start, end]`` subset to write.
        js : bool, optional
            Prepend ``wind = `` so the file can be loaded as a JS module.
        iref : int, optional
            Reference index (passed through to the writer).
        """
        write_wind_to_json(self, file_name, time_range=time_range, iref=iref, js=js)

    def get_coordinates(self, *args) -> tuple[np.ndarray, np.ndarray]:
        """Return the horizontal coordinate arrays of the dataset.

        Parameters
        ----------
        *args
            Optional subset index (int).  Defaults to ``0``.

        Returns
        -------
        x : np.ndarray
            Longitude or x coordinates.
        y : np.ndarray
            Latitude or y coordinates.
        """
        if len(self.subset) > 0:
            if len(args) > 0:
                isub = args[0]
            else:
                isub = 0
            ds = self.subset[isub].ds
        else:
            ds = self.ds
        if self.crs.is_geographic:
            x = ds["lon"].to_numpy()
            y = ds["lat"].to_numpy()
        else:
            x = ds["x"].to_numpy()
            y = ds["y"].to_numpy()
        return x, y

    def get_times(self, *args) -> np.ndarray:
        """Return the time coordinate array of the dataset.

        Parameters
        ----------
        *args
            Optional subset index (int).  Defaults to ``0``.

        Returns
        -------
        np.ndarray
            Array of ``np.datetime64`` values.
        """
        if len(self.subset) > 0:
            if len(args) > 0:
                isub = args[0]
            else:
                isub = 0
            ds = self.subset[isub].ds
        else:
            ds = self.ds
        t = ds["time"].to_numpy()
        return t


def add_time_coord(ds: xr.Dataset, date: datetime.datetime, moving: bool) -> xr.Dataset:
    """Add a time coordinate to an xarray Dataset.

    Parameters
    ----------
    ds : xr.Dataset
        Input dataset (typically a single time-slice read from file).
    date : datetime
        The timestamp to assign as the time coordinate.
    moving : bool
        If ``True``, promote ``lat``/``lon`` from coordinates to data
        variables so they can vary with time (required for moving grids).

    Returns
    -------
    xr.Dataset
        Dataset with a new ``time`` dimension of length 1.
    """
    # Inject a new time coordinate
    ds = ds.expand_dims(time=[date])  # Make sure time has length 1

    # Promote lat/lon from coordinates to variables if needed for moving grids
    if moving:
        for coord in ["lat", "lon"]:
            # Expand along time and assign back — even if it's a coord
            ds[coord] = ds[coord].expand_dims(time=[date])
        # Move from coordinates to variables
        ds = ds.reset_coords(["lat", "lon"])
    return ds


def get_buffered_slice(coord_array: np.ndarray, min_val: float, max_val: float, buffer: int = 2) -> slice:
    """Return a coordinate slice with a cell buffer around the given value range.

    Parameters
    ----------
    coord_array : np.ndarray
        Sorted 1-D coordinate array.
    min_val : float
        Lower bound of the desired range.
    max_val : float
        Upper bound of the desired range.
    buffer : int, optional
        Number of extra cells to include on each side.  Default ``2``.

    Returns
    -------
    slice
        A ``slice`` object using coordinate values (not integer indices).
    """
    idx_min = max(np.searchsorted(coord_array, min_val) - buffer, 0)
    idx_max = min(
        np.searchsorted(coord_array, max_val, side="right") + buffer,
        len(coord_array) - 1,
    )
    return slice(coord_array[idx_min], coord_array[idx_max])
