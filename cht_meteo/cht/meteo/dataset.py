import os
import datetime
import xarray as xr
import numpy as np
import pandas as pd
from pyproj import CRS
from pyproj import Transformer
from scipy.interpolate import RegularGridInterpolator
import time

from .dataset_to_delft3d import write_to_delft3d_ascii
from .dataset_to_json_wind import write_wind_to_json

import cht_utils.fileops as fo

class MeteoDataset():
    def __init__(self, **kwargs):

        # Default values set to None 
        self.name = None               # Name of the dataset  
        self.path = None               # Path to the dataset
        self.long_name = None          # Long name of the dataset (not really needed)
        self.source = None             # Source of the dataset
        self.parameters = None         # Parameters in the dataset (a list with strings). Not currently used.
        self.lon_range = None          # Longitude range of the dataset
        self.lat_range = None          # Latitude range of the dataset
        self.var_names = ["wind_u", "wind_v", "barometric_pressure", "precipitation"] # Variable names in the dataset
        self.crs = CRS(4326)           # Coordinate reference system of the dataset
        self.tau = 0                   # Time interval in hours between cycle and data
        self.last_analysis_time = None # Time of last analysis in the dataset

        # Loop through kwargs to set attributes
        reserved_keys = ["x", "y", "lon", "lat", "time"]    
        for key, value in kwargs.items():            
            if key not in reserved_keys:
                setattr(self, key, value)  

        # Set some source information
        self.source_name = ""               # Name of the source
        self.source_type = "forecast"       # Can be "forecast" or "analysis" 
        self.source_delay = 0               # Delay in hours before data is available
        self.source_cycle_interval = 6      # Interval in hours between cycles
        self.source_time_interval = 3       # Output interval in hours
        self.source_forecast_duration = 240 # The length of the forecast in hours that will be downloaded

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

    def init_ds(self, time, x, y):
        """Fill the dataset with nodata values."""

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
        empty_data = np.empty((len(time), len(y), len(x))) + np.nan  
        for var_name in self.var_names:
            if self.crs.is_geographic:
                self.ds[var_name] = xr.DataArray(empty_data, dims=("time", "lat", "lon"))    
            else:
                self.ds[var_name] = xr.DataArray(empty_data, dims=("time", "y", "x"))

    def download(self, time_range, **kwargs):

        # Make path
        os.makedirs(self.path, exist_ok=True)

        if self.source_type == "forecast":
            # Download forecast from cycles
            self.download_forecast(time_range, **kwargs)
        else:
            # Download analysis
            self.download_analysis(time_range, **kwargs)

    def download_forecast(self, time_range, **kwargs):            # Need to check on previous cycles
        """"Download all (!) data from a forecast dataset. This is done by downloading all cycles within the time range."""

        # Check if this dataset has a download_forecast_cycle method
        if not hasattr(self, "download_forecast_cycle"):
            print(f"Error: download_forecast_cycle method not implemented for dataset {self.name}")
            return

        # Round first time in range down to hour        
        h0 = time_range[0].hour
        # Round down to cycle interval
        h0 = h0 - np.mod(h0, self.source_cycle_interval)
        t0 = time_range[0].replace(microsecond=0, second=0, minute=0, hour=h0,
                                    tzinfo=datetime.timezone.utc)

        # Round last time in range down to hour        
        h1 = time_range[1].hour
        # Round down to cycle interval
        h1 = h1 - np.mod(h1, self.source_cycle_interval)
        t1 = time_range[1].replace(microsecond=0, second=0, minute=0, hour=h1,
                                    tzinfo=datetime.timezone.utc)
        
        print("Now : " + datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%d_%H%M%S"))

        # Determine latest available cycle  
        t_latest = datetime.datetime.now(datetime.timezone.utc) - \
            datetime.timedelta(hours=self.source_delay)
        h0 = t_latest.hour
        h0 = h0 - np.mod(h0, self.source_cycle_interval)
        t_latest = t_latest.replace(microsecond=0, second=0, minute=0, hour=h0)
        print("t_latest : " + t_latest.strftime("%Y%m%d_%H%M%S"))

        # Make sure t0 and t1 are not later than t_latest        
        t0 = min(t0, t_latest)
        t1 = min(t1, t_latest)

        print(f"Downloading from {self.name} - cycles : " + t0.strftime("%Y%m%d_%Hz") + " to " + t1.strftime("%Y%m%d_%Hz"))

        # Make list with cycle times
        cycle_times = pd.date_range(start=t0,
                                    end=t1,
                                    freq=str(self.source_cycle_interval) + 'h').to_pydatetime().tolist()
        
        # Loop through all cycle times
        for it, t in enumerate(cycle_times):
            print(f"Downloading {it + 1} of {len(cycle_times)} - cycle : " + t.strftime("%Y%m%d_%Hz"))
            try:
                self.download_forecast_cycle(cycle_time=t,
                                             **kwargs)
            except Exception as e:
                print(f"Error downloading data from dataset {self.name} - cycle : {t.strftime('%Y%m%d_%Hz')}")
                print(e)

    def download_analysis(self, time_range, **kwargs):

        # Check if this dataset has a download_forecast_cycle method
        if not hasattr(self, "download_analysis_times"):
            print(f"Error: download_analysis_times method not implemented for dataset {self.name}")
            return

        # Make list of requested times
        freqstr = str(self.source_time_interval) + 'h'
        requested_times = pd.date_range(start=time_range[0],
                            end=time_range[1],
                            freq=freqstr).to_pydatetime().tolist()
        
        # Loop through all requested times and see which data is already available
        rtimes = []
        # Check which files do not yet exist
        for t in requested_times:
            file_name = os.path.join( self.path, f"{ self.name }.{ t.strftime('%Y%m%d_%H%M')}.nc" )
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

    def collect(self, time_range, **kwargs):
        """Merge data from separate netcdf files. The actual data is stored in the xarray dataset self.ds. This is a list, so we can optionally store multiple subsets of the data."""

        if "tau" in kwargs:
            tau = kwargs["tau"]
        else:
            tau = self.tau

        last_cycle_time = None
        if "last_cycle" in kwargs:
            if kwargs["last_cycle"] is not None:
                # last_cycle_time = datetime.datetime.strptime(kwargs["last_cycle"], "%Y%m%d_%H")
                last_cycle_time = kwargs["last_cycle"]

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
                    t_cycle = datetime.datetime.strptime(cycle_path[-12:-1], "%Y%m%d_%H")

                    # If last_cycle has been provided, do not get data from later cycles
                    if last_cycle_time is not None:
                        if t_cycle > last_cycle_time:
                            # We should stop now
                            break

                    # Check if time of this cycle falls within requested range
                    if t_cycle < time_range[0] or t_cycle > time_range[1]:
                        # We can skip this cycle
                        continue

                    # Find all times available in this cycle as it may contain our data
                    files_in_cycle = fo.list_files(os.path.join(cycle_path, "*" + subsetstr + "*.nc"))

                    icycle += 1

                    # Loop through all files in this cycle
                    for ifile, file in enumerate(files_in_cycle):

                        t_file = datetime.datetime.strptime(file[-16:-3],
                                                            "%Y%m%d_%H%M")
                        
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

                # fill_missing_at_start = True
                # if fill_missing_at_start:
                #     # Use earliest forecast file available for spinup (force earlier timesteps with this file)
                #     id1 = next(i for i,v in enumerate(requested_files) if v is not None)
                #     for ind in range(id1):
                #         requested_files[ind] = requested_files[id1]

                # # Get rid of None values
                # times_to_remove = []
                # for ind, file in enumerate(requested_files):
                #     if not file:
                #         times_to_remove.append(requested_times[ind])
                # if times_to_remove:  
                #     for tr in times_to_remove:
                #         requested_times.remove(tr)
                # requested_files = [value for value in requested_files if value != None]

                # # Turn time array into nump array
                # requested_times = np.array(requested_times)

            else:

                # A lot easier
                files_in_cycle = fo.list_files(os.path.join(self.path, "*" + subsetstr + "*.nc"))
                for file in files_in_cycle:
                    t_file = datetime.datetime.strptime(file[-16:-3],
                                                        "%Y%m%d_%H%M")
                    if t_file>=time_range[0] and t_file<=time_range[1]:
                        file_list.append(os.path.join(self.path, file))
                        time_list.append(t_file)
                        
            # Now we loop through the files, read them and store them in large array
            time = np.array(time_list)
            
            if not time_list:
                print("No meteo data files found within requested time range")
                return

            ntime = len(time_list)

            # Read in first file to get lons and lats
            ds0 = xr.open_dataset(file_list[0])
            lon = ds0["lon"].values[:]
            if lon[0]>180.0:
                lon = lon - 360.0
            lat = ds0["lat"].values[:]

            # Create new dataset
            ds = xr.Dataset()
            # Add time dimension
            ds["time"] = xr.DataArray(time, dims=("time"))
            # Add the lon and lat
            if moving:
                lons = np.empty((ntime, len(lon)))
                lats = np.empty((ntime, len(lat)))
                ds["lon"] = xr.DataArray(lons, dims=("time", "lon"))
                ds["lat"] = xr.DataArray(lats, dims=("time", "lat"))
            else:
                ds["lon"] = xr.DataArray(lon, dims=("lon"))
                ds["lat"] = xr.DataArray(lat, dims=("lat"))

            # First we create empty data arrays
            for var in self.var_names:
                if var == "wind_u" or var == "wind_v" or var == "precipitation":
                    vdefault = 0.0
                elif var == "barometric_pressure":
                    vdefault = 101300.0
                v = np.zeros((ntime, len(lat), len(lon))) + vdefault    
                ds[var] = xr.DataArray(v, dims=("time", "lat", "lon"))

            # Now loop through times
            for it, file in enumerate(file_list):
                # Read in file
                dsin = xr.open_dataset(file)
                if moving:
                    lon = dsin["lon"].values
                    lat = dsin["lat"].values
                    ds["lon"][it,:] = lon
                    ds["lat"][it,:] = lat
                for var in self.var_names:
                    if var in dsin:
                        # Get the data
                        ds[var][it,:,:] = dsin[var].values

            if subsets:
                # Store the data in the subset
                self.subset[isub].ds = ds
            else:
                self.ds = ds

        return ds

    def cut_out(self,
                name=None,
                lon_range=None,
                lat_range=None,
                time_range=None,
                x_range=None,
                y_range=None,
                lon=None,
                lat=None,
                x=None,
                y=None,
                dx=None,
                dy=None,
                copy_time=False,
                crs=None):
        """Returns a new meteo dataset with a cut-out of this dataset. If the original dataset has subsets, these will first be interpolated onto to largest (coarsest) subset."""

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
        if x_range is None and y_range is None and dx is None and dy is None and x is None and y is None:
            # No cut-out, just interpolation
            x_range = [-1.0e9, 1.0e9]
            y_range = [-1.0e9, 1.0e9]

        # Option 2
        if x_range is not None and y_range is not None and dx is None and dy is None:
            # Cut-out only (but possible interpolation)
            if crs != self.crs:
                print("Error: cut-out with different CRS requires arrays x and y, or dx and dy in combination with x_range and y_range")
                return
            # Get coordinates (in case of subsets, we take the largest (coarsest) subset)
            x, y = self.get_coordinates(0)
            # Let's see if x and y are regularly spaced
            if np.min(np.diff(x)) == np.max(np.diff(x)) and np.min(np.diff(y)) == np.max(np.diff(y)):
                # Regular spacing, all good
                pass
            else:
                # Irregular spacing
                # Create new regular x and y arrays with similar spacing
                dx = np.mean(np.diff(x))
                dy = np.mean(np.diff(y))
                x = np.arange(x[0], x[-1], dx)
                y = np.arange(y[0], y[-1], dy)

            # Limit x, y
            ix0 = np.where(x<=x_range[0])[0]
            if len(ix0) > 0:
                ix0 = ix0[-1]
            else:
                ix0 = 0

            ix1 = np.where(x>=x_range[1])[0]
            if len(ix1) > 0:
                ix1 = ix1[0]
            else:
                ix1 = len(x)

            iy0 = np.where(y<=y_range[0])[0]
            if len(iy0) > 0:
                iy0 = iy0[-1]
            else:
                iy0 = 0

            iy1 = np.where(y>=y_range[1])[0]
            if len(iy1) > 0:
                iy1 = iy1[0]
            else:
                iy1 = len(y)

            if ix1<=ix0 or iy1<=iy0:
                print("Error: cut-out is empty with given x_range and y_range")
                return    

            x = x[ix0:ix1]
            y = y[iy0:iy1]

        # Option 3
        elif x_range is not None and y_range is not None and dx is not None and dy is not None:
            # Create new x and y arrays
            x = np.arange(x_range[0], x_range[1], dx)
            y = np.arange(y_range[0], y_range[1], dy)
        
        # Option 4
        elif x is not None and y is not None:
            # No need to do anything
            pass

        if len(self.subset) > 0:
            t = self.subset[0].ds["time"].values
        else:
            t = self.ds["time"].values

        if time_range is not None:
            # convert values in time_range to np.datetime64
            time_range = [np.datetime64(t) for t in time_range]
            t = t[(t>=time_range[0]) & (t<=time_range[1])]

        # Create new dataset
        dataset = MeteoDataset(name=name, x=x, y=y, time=t, crs=crs)
        dataset.interpolate_dataset(self)

        return dataset

    def interpolate_dataset(self, dataset, copy_time=False):
        """Interpolate data from another dataset to this dataset."""
        # Dimensions have already been set (also time?)
        if copy_time:
            self.ds["time"] = dataset.ds["time"]
        # Get horizontal coordinates
        if "x" in self.ds:
            x = self.ds["x"].values
            y = self.ds["y"].values
        else: 
            x = self.ds["lon"].values
            y = self.ds["lat"].values
        xg, yg = np.meshgrid(x, y)
        # Loop through variables
        for var_name in self.var_names:
            da = self.ds[var_name].copy()
            for it, t in enumerate(self.ds["time"].values):
                # Get data
                v = dataset.interpolate_variable(var_name, t, xg, yg, crs=self.crs)
                # Set points in v equal to points in original data vori where vori already has a value
                vori = da.loc[dict(time=t)].values[:]
                not_nan = np.where(~np.isnan(vori))
                v[not_nan] = vori[not_nan]
                da.loc[dict(time=t)] = v

            self.ds[var_name] = da

    def merge_datasets(self, datasets, **kwargs):
        """Merge datasets. This is useful when we have multiple datasets with different resolutions."""
        # For now we just return the first dataset
        for dataset in datasets:
            self.interpolate_dataset(dataset)
    
    def interpolate_variable(self, var_name: str, t: datetime.datetime, x: np.array, y: np.array, crs=None, fill_missing_data=False):
        """Returns numpy array with interpolated values of quantity at requested_time and lon, lat. If quantity is a vector, the function returns the x-component of the vector. If the quantity is not found, the function returns an array with zeros."""

        # Check shape of x and y. They can be either:
        # 1) 1D arrays with regular spacing (turn into grid)
        # 2) 1D arrays with irregular spacing (treat them as point cloud) 
        # 3) 2D arrays with regular or irregular spacing
        if len(np.shape(x)) == 1 and len(np.shape(y)) == 1:
            # 1D arrays
            # Check if x and y are regularly spaced
            if np.min(np.diff(x)) == np.max(np.diff(x)) and np.min(np.diff(y)) == np.max(np.diff(y)):
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
            transformer = Transformer.from_crs(crs,
                                               self.crs,
                                               always_xy=True)
            xg, yg = transformer.transform(xg, yg)

        # Check if there are subsets 
        if len(self.subset) > 0:
            # There are subsets in this dataset. For now we automatically collect data of each subset.
            subsets = True
            nsub    = len(self.subset)
        else:
            subsets = False
            nsub    = 1

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
            if t in ds["time"].values[:]:
                # Get data at time t
                da = ds[var_name].sel(time=t)
            else:
                # Interpolate data at time t                
                da = ds[var_name].interp(time=t)

            # Get horizontal coordinates
            if self.crs.is_geographic:
                x = da["lon"].values
                y = da["lat"].values
            else:
                x = da["x"].values
                y = da["y"].values

            # Make interpolator
            interp = RegularGridInterpolator((y, x), da.values[:])

            # Find points outside of grid
            iout = np.where((xg<np.min(x)) | (xg>np.max(x)) | (yg<np.min(y)) | (yg>np.max(y)))
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

    def to_netcdf(self, **kwargs):
        """Write to netcdf files. This is not implemented yet."""
        if "filename" in kwargs:
            # Write to single file
            self.ds.to_netcdf(path=kwargs["filename"])        
        else:    
            # Write to database
            os.makedirs(self.path, exist_ok=True)
            # Loop through times in ds
            time = self.ds["time"].values
            for it, t in enumerate(time):
                time_string = pd.to_datetime(t).strftime("%Y%m%d_%H%M")
                file_name = self.name + "." + time_string + ".nc"
                full_file_name = os.path.join(self.path, file_name)
                ds_time = self.ds.isel(time=it)
                # Remove time and reftime
                ds_time = ds_time.drop_vars(["time", "reftime"])
                ds_time.to_netcdf(path=full_file_name)
                ds_time.close()

    def to_delft3d(self, file_name, version="1.03", path=None, header_comments=False, refdate=None, parameters=None, time_range=None, format="ascii"):

        if len(self.subset) > 0:
            # There are subsets in this dataset. For now we automatically collect data of each subset.
            print("Warning! Cannot write meteo dataset with subsets to delft3d format ! Please consider merging the subsets first.")
            return

        if format == "ascii":
            write_to_delft3d_ascii(self, file_name, version, path, header_comments, refdate, parameters, time_range)
        # else:
        #     write_to_delft3d_netcdf(self, file_name, version, path, header_comments, refdate, parameters, time_range)

    def wind_to_json(self, file_name, time_range=None, js=True, iref=1):
        write_wind_to_json(self, file_name, time_range=time_range, iref=iref, js=js)

    def get_coordinates(self, *args):
        """Returns the horizontal coordinates of the dataset."""
        if len(self.subset) > 0:
            if len(args) > 0:
                isub = args[0]
            else:
                isub = 0
            ds = self.subset[isub].ds
        else:
            ds = self.ds
        if self.crs.is_geographic:
            x = ds["lon"].values
            y = ds["lat"].values
        else:
            x = ds["x"].values
            y = ds["y"].values
        return x, y 

    def get_times(self, *args):
        """Returns the times of the dataset. In np64."""
        if len(self.subset) > 0:
            if len(args) > 0:
                isub = args[0]
            else:
                isub = 0
            ds = self.subset[isub].ds
        else:
            ds = self.ds
        t = ds["time"].values
        return t

