import os
import shutil
import zipfile

import pandas as pd
import xarray as xr
import cdsapi
# from ecmwf.opendata import Client
# from tqdm import tqdm
import cfgrib 

from .dataset import MeteoDataset


class MeteoDatasetERA5Reanalysis0p25(MeteoDataset):
    # Inherit from MeteoDomain
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Set some source information
        self.source_name = "era5_reanalysis_0p25"
        self.source_type = "reanalysis"
        self.source_time_interval = 1


    def download_reanalysis(self, time_range, chunk="day"):
        """
        Downloads ERA5 reanalysis data for the given time range.

        Parameters
        ----------
        time_range : list of datetimes
            [start_time, end_time]
        chunk : str
            "day" (default) or "month" — determines download chunk size
        """

        tmp_path = os.path.join(self.path, "_TMP_era5_data")
        os.makedirs(tmp_path, exist_ok=True)

        lon_lim = self.lon_range
        lat_lim = self.lat_range

        # Define list of time chunks
        if chunk == "month":
            periods = pd.date_range(start=time_range[0], end=time_range[1], freq="MS").to_pydatetime()
            get_chunk = lambda dt: (dt.year, dt.month, None)
        else:
            periods = pd.date_range(start=time_range[0], end=time_range[1], freq="D").to_pydatetime()
            get_chunk = lambda dt: (dt.year, dt.month, dt.day)

        for dt in periods:
            year, month, day = get_chunk(dt)

            if chunk == "month":
                nc_zip_file = os.path.join(tmp_path, f"_TMP_era5_{year}{month:02d}.zip")
            else:
                nc_zip_file = os.path.join(tmp_path, f"_TMP_era5_{year}{month:02d}{day:02d}.zip")

            # Check if already exists
            if os.path.exists(nc_zip_file.replace(".zip", ".nc")):
                print(f"ERA5 reanalysis for {year}-{month:02d}" + (f"-{day:02d}" if day else "") + " already exists. Skipping.")
                continue

            print(f"📥 Downloading ERA5 {chunk}: {year}-{month:02d}" + (f"-{day:02d}" if day else ""))

            # Call appropriate downloader
            if chunk == "month":
                download_era5_month(year, month, lon_lim, lat_lim, nc_zip_file)
            else:
                download_era5_day(year, month, day, lon_lim, lat_lim, nc_zip_file)

            # =============================================================
            # PROCESS DOWNLOADED DATA
            # =============================================================
            with zipfile.ZipFile(nc_zip_file, 'r') as zip_ref:
                zip_ref.extractall(tmp_path)

            # Load instant fields
            nc_file = os.path.join(tmp_path, "data_stream-oper_stepType-instant.nc")
            ds = xr.load_dataset(nc_file)
            ds = ds.rename({"latitude": "lat", "longitude": "lon", "valid_time": "time"})
            ds = ds.rename({"u10": "wind_u", "v10": "wind_v", "msl": "barometric_pressure"})

            # Load precipitation
            nc_file = os.path.join(tmp_path, "data_stream-oper_stepType-accum.nc")
            ds_precip = xr.load_dataset(nc_file)
            ds_precip = ds_precip.rename({"latitude": "lat", "longitude": "lon", "valid_time": "time"})
            ds["total_precipitation"] = ds_precip["tp"]

            # Sort lat ascending (south → north)
            ds = ds.sortby("lat")

            # Save
            write2nc(ds, self.name, os.path.join(self.path))

            # Cleanup
            for f in ["data_stream-oper_stepType-instant.nc", "data_stream-oper_stepType-accum.nc"]:
                try:
                    os.remove(os.path.join(tmp_path, f))
                except FileNotFoundError:
                    pass

    def download_reanalysis_ori(self, time_range):

        """Downloads ERA5 reanalysis data for a given storm number and cycle time"""

        tmp_path = os.path.join(self.path, "_TMP_era5_data")

        # Make folder for the forecast
        os.makedirs(tmp_path, exist_ok=True)

        lon_lim = self.lon_range  # degrees east
        lat_lim = self.lat_range  # degrees north

        # Loop through days
        days = pd.date_range(
            start=time_range[0], end=time_range[1], freq="D"
        ).to_pydatetime()

        download = False

        for dy in days:

            # Check if data for this day already exists
            ok = True
            for it in range(24):
                time_string = dy.strftime("%Y%m%d") + f"_{it:02d}00"
                file_name = os.path.join(self.path,
                                         self.name + "." + time_string + ".nc")
                if not os.path.exists(file_name):
                    ok = False
                    break
            if ok:
                # All files for this day already exist. Skip to next day.
                print(f"ERA5 reanalysis data for {dy.strftime('%Y-%m-%d')} already exists. Skipping download.")
                continue    

            nc_zip_file = os.path.join(tmp_path, f"_TMP_era5_{dy.year}{dy.month:02d}{dy.day:02d}.zip")

            download = True

            if download:
                download_era5_day(
                    dy.year,
                    dy.month,
                    dy.day,
                    lon_lim,
                    lat_lim,
                    nc_zip_file,
                )

            # =============================================================
            # LOAD, CLIP, AND MERGE DATASETS
            # =============================================================
            
            # Unzip files
            with zipfile.ZipFile(nc_zip_file, 'r') as zip_ref:
                zip_ref.extractall(tmp_path)

            # If everything went well, we should have two netcdf files now
            # 1) data_stream-oper_stepType-instant.nc (wind + pressure)
            # 2) data_stream-oper_stepType-accum.nc (total precipitation)

            # Start with first one
            nc_file = os.path.join(tmp_path, "data_stream-oper_stepType-instant.nc")

            ds = xr.load_dataset(nc_file)
            ds = ds.rename({"latitude": "lat", "longitude": "lon"})
            rename_dict = {
                "u10": "wind_u",
                "v10": "wind_v",
                "msl": "barometric_pressure",
            }
            ds = ds.rename(rename_dict)
            ds = ds.rename({"valid_time": "time"})

            # =============================================================
            # COMPUTE RAINFALL RATE (mm/hour)
            # =============================================================

            # # ECMWF total precipitation (tp) is in meters since forecast start
            # # We compute 3-hour increments and convert to mm/hour
            # if "tp" in ds_merged:
            #     tp = ds_merged["tp"]

            #     # Compute increments over 3-hour intervals
            #     tp_diff = tp.diff("step")  # tp[1]-tp[0], tp[2]-tp[1], etc.

            #     # Convert to mm/h (ECMWF tp is in meters since forecast start)
            #     tp_rate = tp_diff * 1000 / 3.0

            #     # Shift step coordinate backward by one interval, so 0–3h rain stored at 0h
            #     tp_rate = tp_rate.assign_coords(step=tp.step[:-1])
            #     tp_rate.name = "rain_rate"
            #     tp_rate.attrs["units"] = "mm/h"
            #     tp_rate.attrs["description"] = (
            #         "3-hourly mean precipitation rate (interval start time)"
            #     )

            #     # Merge back (drop cumulative field)
            #     ds_merged = ds_merged.drop_vars("tp")

            #     ds_merged = xr.merge([ds_merged, tp_rate])

            #     # Replace any nans in ds_merged["rain_rate"] with zeros (no rain)
            #     ds_merged["rain_rate"] = ds_merged["rain_rate"].fillna(0.0)

            # Change variable names to more standard ones
            # Rename latitude and longitude to lat and lon
            # ds = ds.rename({"latitude": "lat", "longitude": "lon"})
            # rename_dict = {
            #     "u10": "wind_u",
            #     "v10": "wind_v",
            #     "msl": "barometric_pressure",
            #     "rain_rate": "precipitation",
            # }
            # ds = ds.rename(rename_dict)

            # Now the precip file
            nc_file = os.path.join(tmp_path, "data_stream-oper_stepType-accum.nc")
            ds_precip = xr.load_dataset(nc_file)
            ds_precip = ds_precip.rename({"latitude": "lat", "longitude": "lon"})
            ds_precip = ds_precip.rename({"valid_time": "time"})
            # Now copy ds_precip's total precipitation into ds
            ds["total_precipitation"] = ds_precip["tp"]

            # And finalize flip of latitude to be south to north
            ds = ds.sortby("lat")

            # =============================================================
            # SAVE TO NETCDF
            # =============================================================

            write2nc(ds, self.name, os.path.join(self.path))

            # ds.close()
            # ds_precip.close()

            # Delete temporary files (nc files in tmp_path)
            try:
                os.remove(os.path.join(tmp_path, "data_stream-oper_stepType-instant.nc"))
                os.remove(os.path.join(tmp_path, "data_stream-oper_stepType-accum.nc"))
            except Exception as e:
                print(f"Error deleting temporary files: {e}")
            # os.remove(nc_zip_file)
            
            # shutil.rmtree(tmp_path, ignore_errors=True)


def download_era5_day(year, month, day, lon_range, lat_range, pth):

    dataset = "reanalysis-era5-single-levels"
    request = {
        "product_type": ["reanalysis"],
        "variable": [
            "10m_u_component_of_wind",
            "10m_v_component_of_wind",
            "mean_sea_level_pressure",
            "total_precipitation"
        ],
        "year": [str(year)],
        "month": [f"{month:02d}"],
        "day": [f"{day:02d}"],
        "time": [
            "00:00", "01:00", "02:00",
            "03:00", "04:00", "05:00",
            "06:00", "07:00", "08:00",
            "09:00", "10:00", "11:00",
            "12:00", "13:00", "14:00",
            "15:00", "16:00", "17:00",
            "18:00", "19:00", "20:00",
            "21:00", "22:00", "23:00"
        ],
        #"format": "grib",
        "data_format": "netcdf",
        # "download_format": "unarchived",
        "area": [lat_range[1], lon_range[0], lat_range[0], lon_range[1]]
    }
    
    client = cdsapi.Client()
    client.retrieve(dataset, request, target=pth)

def download_era5_month(year, month, lon_range, lat_range, pth):
    dataset = "reanalysis-era5-single-levels"
    # Get the available days in the month (january always has 31 days, february 28 or 29, etc.)
    if month == 2:
        if (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0):
            days_in_month = 29
        else:
            days_in_month = 28
    elif month in [4, 6, 9, 11]:
        days_in_month = 30
    else:
        days_in_month = 31
    # Make day list
    day_list = [f"{day:02d}" for day in range(1, days_in_month + 1)]    

    request = {
        "product_type": ["reanalysis"],
        "variable": [
            "10m_u_component_of_wind",
            "10m_v_component_of_wind",
            "mean_sea_level_pressure",
            "total_precipitation"
        ],
        "year": [str(year)],
        "month": [f"{month:02d}"],
        "day": day_list,
        "time": [
            "00:00", "01:00", "02:00",
            "03:00", "04:00", "05:00",
            "06:00", "07:00", "08:00",
            "09:00", "10:00", "11:00",
            "12:00", "13:00", "14:00",
            "15:00", "16:00", "17:00",
            "18:00", "19:00", "20:00",
            "21:00", "22:00", "23:00"
        ],
        #"format": "grib",
        "data_format": "netcdf",
        # "download_format": "unarchived",
        "area": [lat_range[1], lon_range[0], lat_range[0], lon_range[1]]
    }
    
    client = cdsapi.Client()
    client.retrieve(dataset, request, target=pth)

def write2nc(ds, meteo_name, meteo_path):
    # Loop though times in ds and save each to a separate netcdf file
    times = ds["time"].to_numpy()
    for it, t in enumerate(times):
        time_string = pd.to_datetime(t).strftime("%Y%m%d_%H%M")
        file_name = meteo_name + "." + time_string + ".nc"
        full_file_name = os.path.join(meteo_path, file_name)
        ds_time = ds.isel(time=it)
        # Remove time and reftime
        ds_time = ds_time.drop_vars(["time"])
        ds_time.to_netcdf(path=full_file_name)
        ds_time.close()
