# -*- coding: utf-8 -*-
import os
import datetime
import urllib
import tarfile as trf
import xarray as xr
import requests

from cht_utils import fileops as fo
from .dataset import MeteoDataset

class MeteoSubset:
    def __init__(self, name, moving):
        self.name = name
        self.moving = moving
        self.ds = xr.Dataset()

class MeteoDatasetCOAMPSTCForecastS3(MeteoDataset):

    # Inherit from MeteoDomain
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Set some source information
        self.source_name = "coamps_tc_forecast"
        self.source_type = "forecast"
        self.source_delay = 6
        self.source_cycle_interval = 6
        self.source_time_interval = 1

        # List subsets
        self.subset = []
        self.subset.append(MeteoSubset("d01", False))
        self.subset.append(MeteoSubset("d02", True))
        self.subset.append(MeteoSubset("d03", True))
    

    def download_forecast_cycle(self, **kwargs):
        """Downloads COAMPS-TC forecast cycle from S3 for a given storm number and cycle time"""

        if "cycle_time" in kwargs:
            cycle_time = kwargs["cycle_time"]
        else:
            # Throw error if cycle_time is not provided            
            print("Error: cycle_time not provided")
            return
        
        if "storm_number" in kwargs:
            storm_number = kwargs["storm_number"]
        else:
            # Throw error if storm_number is not provided            
            print("Error: storm_number not provided")
            return
        
        if "only_track" in kwargs:
            only_track = kwargs["only_track"]
        else:
            only_track = False

        # For COAMPS-TC we always get the entire dataset (lon_range, lat_range time_range are not used)

        base_url = "https://coamps-tc-data.s3.us-east-2.amazonaws.com"
        base_url = base_url + "/deterministic/realtime"

        # Get year from cycle_time
        year = cycle_time.year
        cycle_time_coamps = cycle_time.strftime("%Y%m%d%H")
        cycle_time_meteo  = cycle_time.strftime("%Y%m%d_%HZ")

        # Make folder for the forecast
        forecast_path = os.path.join(self.path, cycle_time_meteo) # Folder to store the forecast netcdf files
        os.makedirs(forecast_path, exist_ok=True) 

        # Start with the track file
        get_storm_track(forecast_path, year, storm_number, cycle_time_coamps)

        if only_track:
            # No need to download the gridded forecast, apparently
            return

        # Make folder for the temporary files
        tar_file_path = os.path.join(forecast_path, "_TMP")  # Temporary folder to store the tar file
        os.makedirs(tar_file_path, exist_ok=True)

        tar_file_name = storm_number + "_" + cycle_time_coamps + "_netcdf.tar"

        tar_file_full_path = os.path.join(tar_file_path, tar_file_name)

        url = base_url + "/" + str(year) + "/" + storm_number + "/" + cycle_time_coamps + "/" + tar_file_name

        # Download the tar file
        print(f"Downloading {url}")
        try:    
            urllib.request.urlretrieve(url, tar_file_full_path)
        except Exception as e:
            print(f"Error downloading {url}")
            print(e)
            # Remove the forecast folder
            fo.delete_folder(forecast_path)
            return

        # Extract the tar file
        with trf.open(tar_file_full_path, "r") as tar:
            for member in tar.getmembers():
                tar.extract(member, path=tar_file_path)

        # Convert all three resolutions to meteo format
        print("Converting COAMPS-TC netcdf files to meteo netcdf files")
        for subset in self.subset:
            res = subset.name
            # Get list of all the files in the _TMP/netcdf folder
            file_list = fo.list_files(os.path.join(tar_file_path, "netcdf", "*_" + res + "_*"))
            for file_name in file_list:
                # Read tau from file name
                tau = file_name.split("_")[-1].split(".")[0][3:]
                t = cycle_time + datetime.timedelta(hours=int(tau))
                tstr = t.strftime("%Y%m%d_%H%M")
                output_file = os.path.join(forecast_path, self.name + "." + res + "." + tstr + ".nc")
                convert_coamps_nc_to_meteo_nc(file_name, output_file)

        # Remove the temporary folder
        fo.delete_folder(tar_file_path)

        # Also try to download the track file
        track_file_name = storm_number + "_" + cycle_time_coamps + "_track.txt"

def convert_coamps_nc_to_meteo_nc(inpfile, outfile):
 
    # Open the COAMPS-TC netcdf file
    dsin = xr.open_dataset(inpfile)

    # Get the lon and lat
    lon = dsin["lon"].values[0, :] - 360.0
    lat = dsin["lat"].values[:, 0]

    # Create new dataset
    ds = xr.Dataset()
    # Add the lon and lat
    ds["lon"] = xr.DataArray(lon, dims=("lon"))
    ds["lat"] = xr.DataArray(lat, dims=("lat"))
    # Add the variables
    variables_coamps = ["uuwind", "vvwind", "slpres", "precip"]
    variables = ["wind_u", "wind_v", "barometric_pressure", "precipitation"]
    # units = ["m/s", "m/s", "hPa", "mm"]
    units = ["m/s", "m/s", "Pa", "mm/h"]
    for ivar, var in enumerate(variables_coamps):
        # Conversion factors
        fconv = 1.0
        if var == "slpres":
            fconv = 100.0
        if var in dsin:
            ds[variables[ivar]] = dsin[var] * fconv
            ds[variables[ivar]].attrs["units"] = units[ivar]

    # Write output file
    ds.to_netcdf(outfile)

    # Close the raw netcdf file
    dsin.close()

def get_storm_track(track_path:str, year:int, storm:str, cycle:str):
    """
    Retrieves the storm track data for a given year, storm, and cycle.

    Parameters:
    year (int): The year of the storm track data.
    storm (str): The name of the storm.
    cycle (str): The cycle of the storm track data.

    Returns:
    bytes: The content of the storm track data.

    """
    try:
        filename = os.path.join(track_path, f"TRK_COAMPS_CTCX_3_{cycle}_{storm}.trk")
        url = f"https://coamps-tc-data.s3.us-east-2.amazonaws.com/deterministic/realtime/{year}/{storm}/{cycle}/TRK_COAMPS_CTCX_3_{cycle}_{storm}"
        urllib.request.urlretrieve(url, filename)
    except Exception as e:
        print(f"Error downloading {url}")
        print(e)
