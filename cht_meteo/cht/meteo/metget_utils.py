import datetime
import pandas as pd
import requests
import tomli
import xarray as xr
from typing import Tuple
from metget.metget_build import MetGetBuildRest
from pathlib import Path

def date_transform(date):
    """
    Transforms a string representation of a date into a datetime object with UTC timezone.

    Args:
        date (str): A string representing a date in the format "%Y-%m-%d %H:%M:%S".

    Returns:
        datetime.datetime: A datetime object with the transformed date and UTC timezone.
    """
    return datetime.datetime.strptime(date, "%Y-%m-%d %H:%M:%S").replace(tzinfo=datetime.timezone.utc)

def get_da_from_url(url):
    """
    Retrieves a dataset from a given URL and returns it as an xarray DataArray.

    Parameters:
    url (str): The URL of the dataset to retrieve.

    Returns:
    xr.DataArray: The dataset as an xarray DataArray.
    """
    from netCDF4 import Dataset as nc_Dataset
    data = requests.get(url).content
    ds0 = nc_Dataset('temp', memory=data)
    ds = xr.open_dataset(xr.backends.NetCDF4DataStore(ds0))
    # Alternative
    # ds = xr.open_dataset(url + '#mode=bytes') # Added this last part to allow opening with xarray
    return ds

def tc_vitals_storm():
    """Find the storm with the highest priority from NOAA"""
    try:
        tcvitals = requests.get('https://ftp.nhc.noaa.gov/atcf/com/tcvitals').text
        splits = [line.split() for line in tcvitals.split('\n')[:-1]] # last line is empty
        priority = pd.DataFrame(splits)
        priority = priority.iloc[[i for i in range(len(priority)) if 'L' in priority.iloc[i, 1]], :] # check only in atlantic basin storms
        pr_st_noaa = priority.iloc[0, 1] # define the top L storm as the priority one
        return pr_st_noaa
    except:
        return None

class MeteoVariable:
    def __init__(self, unit, fill_value, metget_name):
        self.unit = unit
        self.fill_value = fill_value
        self.metget_name = metget_name

class MetGet:
    
    def __init__(self, apikey:str = None, endpoint:str = None, api_version:int = 2, config_path = None):
        if config_path:
            with open(config_path, mode="rb") as fp:
                config = tomli.load(fp)
            self.apikey = config["apikey"]
            self.endpoint = config["endpoint"]
            self.api_version = config["api_version"]
            # Check if there is a priority storm given
            if "priority_storm" in config:
                if config["priority_storm"] == "tc_vitals":
                    # Get storm priority from noaa
                    self.priority_storm = tc_vitals_storm()
                    if not self.priority_storm:
                        print("Priority storm could not be loaded from https://ftp.nhc.noaa.gov/atcf/com/tcvitals")
                    else:
                        print(f"Priority storm {self.priority_storm} found from https://ftp.nhc.noaa.gov/atcf/com/tcvitals")
                else:
                    self.priority_storm = config["priority_storm"]
                    print(f"Priority storm {self.priority_storm} provided in config file.")  
            else:
                self.priority_storm=None
            if "tau" in config:
                self.tau = config["tau"]
            else:
                self.tau = 0
        elif apikey and endpoint:
            self.apikey = apikey
            self.endpoint = endpoint
            self.api_version = api_version
        else:
            raise ValueError("Either 'config_path' or both 'apikey' and 'endpoint' must be provided")
        
        self.meteo_variables = {}
        self.meteo_variables["wind"] = MeteoVariable(unit='m/s', fill_value=0.0, metget_name="wind_pressure")
        self.meteo_variables["barometric_pressure"] = MeteoVariable(unit='Pa', fill_value=102000.0, metget_name="wind_pressure")
        self.meteo_variables["precipitation"] = MeteoVariable(unit='kg.m-2.hour-1', fill_value=0.0, metget_name="rain")

    def get_forecasts(self, model:str = "coamps", start:datetime = None, end:datetime = None) -> dict:
        """
        Read metadata of available forecasts from endpoint and returns a dictionary of the available storms.

        Parameters:
        - start (datetime): The start date for the forecast data.
        - end (datetime): The end date for the forecast data.

        Returns:
        - forecasts (dict): A dictionary containing the available storms and their metadata.
        """
        url = "{:s}/status?model={:s}".format(self.endpoint, model)
        if start:
            # Ensure that we are getting the previous available forecasts as well
            start = start - datetime.timedelta(hours=6)
            url += "&start={:s}".format(start.strftime("%Y-%m-%d"))
        if end:
            url += "&end={:s}".format(end.strftime("%Y-%m-%d"))
        # ...Get the json from the endpoint
        response = requests.get(
            url, headers={"x-api-key": self.apikey}
        )
        data = response.json()["body"]
        
        if not data:
            return {}
        # TODO check if this is consistent with the latest endpoint
        if "data" in data:
            forecasts = data["data"]["metget"]["coamps-tc"]
        else:
            last_year = list(data.keys())[0]
            forecasts = data[last_year]
        
        return forecasts

    def get_coamps_available_storms(self, time_range:Tuple[datetime.datetime, datetime.datetime], dt:int = 1):
        # Get available forecast for the time_range
        coamps_forecasts = self.get_forecasts(model="coamps", start=time_range[0], end=time_range[1])
        
        # Create time-series based on time_range and dt
        requested_times = pd.date_range(start=time_range[0],
                      end=time_range[1],
                      freq=f'{dt}H').to_pydatetime().tolist()
        requested_times = [ti.replace(tzinfo=datetime.timezone.utc) for ti in requested_times]
  
        # Check which of the available storm forecast all named storms 
        storms = [name for name in coamps_forecasts.keys()]
        storms.sort()
        
        # Check if forecast exists for the selected dates and make a dataframe
        df = pd.DataFrame(index=requested_times)
        df["storms"] = None

        for ti in requested_times:
            df.loc[ti, "storms"] = []
            for storm_id in storms:
                storm = coamps_forecasts[storm_id]                 
                t1 = date_transform(storm['min_forecast_date'])
                t2 = date_transform(storm['max_forecast_date'])
                belongs = t1 <= ti <= t2
                if belongs:
                    df.loc[ti, "storms"].append(storm_id)
                    
        return df

    def get_coamps_quasi_hindcast(self, time_range:Tuple[datetime.datetime, datetime.datetime], 
                           dt:int = 1, 
                           priority_storm:str = None,
                           multi_storms_method:str = "last",
                           tau:int = None,
                           l_storms_only = True
                           ):
        if priority_storm is None:
            priority_storm = self.priority_storm
        if tau is None:
            tau = self.tau
        coamps_forecasts = self.get_forecasts(model="coamps", start=time_range[0], end=time_range[1])
        df = self.get_coamps_available_storms(time_range=time_range, dt=dt)
        df["storm_id"] = None
        df["cycle"] = None
        df["hour"] = None
        # Loop through time steps and choose a single storm to use 
        for ti in df.index:
            # Check if for the chosen storm the time is during the skip time (tau) of the first forecast
            for sn in df.loc[ti, "storms"]:
                t0 = date_transform(coamps_forecasts[sn]["first_available_cycle"])
                t0_tau = t0 + datetime.timedelta(hours=tau)
                if ti < t0_tau:
                    df.loc[ti, "storms"].remove(sn)
            storm_name = None # if there is no coamps storm make None so gfs is used
            if l_storms_only:
                storms = [name for name in df.loc[ti, "storms"] if int(name.split("L")[0]) < 90]
            else:
                storms = df.loc[ti, "storms"]
            if len(storms)>0:
                if priority_storm in storms:  # if the priority storm is available that timestep use it
                    storm_name = priority_storm
                else:
                    # else used method defined
                    if multi_storms_method == 'first':
                        storm_name = storms[0]
                    elif multi_storms_method == 'last':
                        storm_name = storms[-1]
            # Assign storm
            df.loc[ti, "storm_id"] = storm_name
        
        # Get available cycles per storm forecast
        cycles = {}
        for storm in df["storm_id"].dropna().unique():
            cycles[storm] = pd.date_range(
                start=date_transform(coamps_forecasts[storm]["first_available_cycle"]),
                end=date_transform(coamps_forecasts[storm]["latest_available_cycle"]),
                freq='6H').to_pydatetime()
            
        # Get forecast used for each time-step
        for ti in df.index:
            if df.loc[ti, "storm_id"] is not None:
                cycle = cycles[df.loc[ti, "storm_id"]][cycles[df.loc[ti, "storm_id"]] <= ti][-1]
                hour = int((ti - cycle).total_seconds() / 3600)
                if hour < tau:
                    cycle = cycles[df.loc[ti, "storm_id"]][cycles[df.loc[ti, "storm_id"]] <= ti][-2]
                    hour = int((ti - cycle).total_seconds() / 3600)
                df.loc[ti, "cycle"] = cycle.strftime("%Y%m%d%Hz")
                df.loc[ti, "hour"] = hour
                
        return df
      
    def get_meteo(self, time_range:Tuple[datetime.datetime, datetime.datetime], var:str, dt:float, domain:list, tau=None, multiple_forecasts=False):
        if tau is None:
            tau = self.tau
        t1, t2 = time_range[0], time_range[1]
        request_data = MetGetBuildRest.generate_request_json(
            start_date=t1.strftime("%Y%m%d %H%M%S"),
            end_date=t2.strftime("%Y%m%d %H%M%S"),
            format='hec-netcdf',
            data_type=var,
            time_step=dt*3600,
            domains=MetGetBuildRest.parse_command_line_domains([domain], tau),
            epsg=4326,
            filename=f"{domain[0]}_{var}",
            backfill=True,
            multiple_forecasts=multiple_forecasts,
            # compression=True,
            save_json_request=False,
            # dry_run=True,
            strict=True,
        )
        client = MetGetBuildRest(self.endpoint, self.apikey, self.api_version)
        data_id, status_code = client.make_metget_request(request_data) # Status code should be 200 if everything is ok
        if status_code != 200:
            raise RuntimeError(f"Unexpected status code: {status_code}, when downloading data.")
        urls, meta = client.download_metget_data(
            data_id,
            30, # sleep_time in seconds
            3, # max_wait in hours
            output_directory=None,
            return_only_url=True,
        )
        
        ds = get_da_from_url(urls[0])
        
        return ds, meta     
            
def get_storm_track(year:int, storm:str, cycle:str):
    """
    Retrieves the storm track data for a given year, storm, and cycle.

    Parameters:
    year (int): The year of the storm track data.
    storm (str): The name of the storm.
    cycle (str): The cycle of the storm track data.

    Returns:
    bytes: The content of the storm track data.

    """
    url = f"https://coamps-tc-data.s3.us-east-2.amazonaws.com/deterministic/realtime/{year}/{storm}/{cycle}/TRK_COAMPS_CTCX_3_{cycle}_{storm}"
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.content
    else:
        data = None
    return data

def save_to_nc(name, path, data):
    path = Path(path)
    for it, t in enumerate(data[0].time):
        time_string = t.strftime("%Y%m%d_%H%M")
        file_name = name + "." + time_string + ".nc"
        full_file_name = path.joinpath(file_name)
        ds = xr.Dataset()

        for dd in data:
            if dd.quantity == "wind":
                uu = dd.u[it, :, :]
                da = xr.DataArray(uu,
                                  coords=[("lat", dd.y),
                                          ("lon", dd.x)])
                if dd.source[it]:
                    da.attrs['source'] = dd.source[it]
                ds["wind_u"] = da
                vv = dd.v[it, :, :]
                da = xr.DataArray(vv,
                                  coords=[("lat", dd.y),
                                          ("lon", dd.x)])
                if dd.source[it]:
                    da.attrs['source'] = dd.source[it]
                ds["wind_v"] = da
            else:
                try:
                    val = dd.val[it, :, :]
                    da = xr.DataArray(val,
                                      coords=[("lat", dd.y),
                                              ("lon", dd.x)])
                    if dd.source[it]:
                        da.attrs['source'] = dd.source[it]
                    ds[dd.quantity] = da
                except:
                    print("Could not write " + dd.quantity + " to file ...")

        ds.to_netcdf(path=full_file_name)

        ds.close()