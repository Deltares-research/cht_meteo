import datetime
import pandas as pd
import requests
import xarray as xr

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


def check_coamps(apikey:str, endpoint:str, start:datetime=None, end:datetime=None) -> dict:
    """
    Read metadata of available forecasts from endpoint and returns a dictionary of the available storms.

    Parameters:
    - apikey (str): The API key for accessing the endpoint.
    - endpoint (str): The URL endpoint for retrieving the forecast data.
    - start (datetime): The start date for the forecast data.
    - end (datetime): The end date for the forecast data.

    Returns:
    - storms (dict): A dictionary containing the available storms and their metadata.
    """
    url = "{:s}/status?model={:s}".format(endpoint, "coamps")
    if start:
        url += "&start={:s}".format(start.strftime("%Y-%m-%d"))
    if end:
        url += "&end={:s}".format(end.strftime("%Y-%m-%d"))
    # ...Get the json from the endpoint
    response = requests.get(
        url, headers={"x-api-key": apikey}
    )
    data = response.json()["body"]
    
    if not data:
        return {}
    # TODO check if this is consistent with the latest endpoint
    if "data" in data:
        storms = data["data"]["metget"]["coamps-tc"]
    else:
        last_year = list(data.keys())[0]
        storms = data[last_year]
    return storms

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