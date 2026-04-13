"""KNMI Harmonie forecast dataset from the Matroos THREDDS service."""

import datetime
import os

import pandas as pd
import xarray as xr
from siphon.catalog import TDSCatalog

from .dataset import MeteoDataset


class MeteoDatasetMatroos(MeteoDataset):
    """KNMI Harmonie forecast dataset served by Deltares Matroos.

    Downloads wind, pressure and precipitation via the Matroos OPeNDAP
    THREDDS catalogue.

    Parameters
    ----------
    **kwargs
        Forwarded to :class:`MeteoDataset`.
    """

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)

        # Set some source information
        self.source_name = "matroos"
        self.source_type = "forecast"
        self.source_delay = 4
        self.source_cycle_interval = (
            1  # real value 1, but maybe increase this to limit downloads ...
        )
        self.source_time_interval = 1
        self.source_forecast_duration = 60

    def download_forecast_cycle(self, **kwargs) -> None:
        """Download a single Matroos forecast cycle.

        Parameters
        ----------
        **kwargs
            Required keys:

            ``cycle_time`` : datetime
                Forecast initialisation time.

            Optional keys:

            ``time_range`` : list of datetime
                ``[start, end]`` window for which data is needed.  Defaults to
                the full forecast duration defined on the dataset.

        Raises
        ------
        ValueError
            If the requested cycle time is earlier than the earliest available
            cycle in the Matroos catalogue.
        """

        if "cycle_time" in kwargs:
            cycle_time = kwargs["cycle_time"]
        else:
            # Throw error if cycle_time is not provided
            print("Error: cycle_time not provided")
            return

        if "time_range" in kwargs:
            time_range = kwargs["time_range"]
        else:
            # Get all data from this cycle
            time_range = [
                cycle_time,
                cycle_time + datetime.timedelta(hours=self.source_forecast_duration),
            ]

        cycle_name = cycle_time.strftime("%Y%m%d_%Hz")

        # Make folder for the forecast
        forecast_path = os.path.join(
            self.path, cycle_name
        )  # Folder to store the forecast netcdf files

        # Make folder for the forecast
        os.makedirs(forecast_path, exist_ok=True)

        # We assume that the best uses the latest
        url = "https://rwsos-dataservices-prod.avi.deltares.nl//thredds//catalog//data_matroos//maps//normal//knmi_harmonie43"
        catalog_xml = f"{url}//catalog.xml"
        catalog = TDSCatalog(catalog_xml)

        # If we request times earlier than the earliest available cycle the earliest has been used
        first_file = catalog.datasets[0].name
        first_date = "".join(first_file.split(".")[0].split("_")[-2:])
        first_cycle_time = datetime.datetime.strptime(
            first_date, "%Y%m%d%H%M%S"
        ).replace(tzinfo=datetime.timezone.utc)
        if time_range[0] < first_cycle_time:
            raise ValueError(
                f"Matroos data not available for times earlier than {first_cycle_time.strftime('%Y%m%d_%H%M')}"
            )

        # If we request times later than the latest available cycle the latest has been used
        latest_file = catalog.datasets[-1].name
        latest_date = "".join(latest_file.split(".")[0].split("_")[-2:])
        latest_cycle_time = datetime.datetime.strptime(
            latest_date, "%Y%m%d%H%M%S"
        ).replace(tzinfo=datetime.timezone.utc)
        if time_range[0] > latest_cycle_time:
            cycle_time = latest_cycle_time
            print(
                f"Warning: Matroos data not available for times later than {latest_cycle_time.strftime('%Y%m%d_%H%M')}, using latest available cycle"
            )

        # now get the actual cycle
        matroos_cycle_string = cycle_time.strftime("%Y%m%d%H%M")
        access_url = catalog.datasets[f"{matroos_cycle_string}.nc"].access_urls[
            "OPENDAP"
        ]
        ds = xr.open_dataset(access_url)

        # subset variables
        variables = [
            "eastward_wind_fixed_height",
            "northward_wind_fixed_height",
            "air_pressure_fixed_height",
            "precipitation_amount",
        ]
        ds = ds[variables]

        # Rename
        ds = ds.rename(
            {
                "eastward_wind_fixed_height": "wind_u",
                "northward_wind_fixed_height": "wind_v",
                "air_pressure_fixed_height": "barometric_pressure",
                "precipitation_amount": "precipitation",
                "x": "lon",
                "y": "lat",
            }
        )

        # Subset in space
        ds = ds.sel(
            lon=slice(self.lon_range[0], self.lon_range[1]),
            lat=slice(self.lat_range[0], self.lat_range[1]),
        )

        # Subset in time
        # NOTE: times in matroos are in UTC without timezone, so drop the timezone info for the comparison
        time_range = [
            t.astimezone(datetime.timezone.utc).replace(tzinfo=None) for t in time_range
        ]
        ds = ds.sel(time=slice(time_range[0], time_range[1]))

        # Write to netcdf files
        write2nc(ds, self.name, os.path.join(self.path, cycle_name))

        ds.close()


def write2nc(ds: xr.Dataset, meteo_name: str, meteo_path: str) -> None:
    """Write one netCDF file per time step in the dataset.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing a ``time`` dimension.
    meteo_name : str
        Dataset name used as the file-name prefix.
    meteo_path : str
        Output directory.
    """
    # Loop though times in ds
    times = ds["time"].to_numpy()
    for it, t in enumerate(times):
        time_string = pd.to_datetime(t).strftime("%Y%m%d_%H%M")
        file_name = f"{meteo_name}.{time_string}.nc"
        full_file_name = os.path.join(meteo_path, file_name)
        ds_time = ds.isel(time=it)
        # Remove time and reftime
        ds_time = ds_time.drop_vars(["time"])
        ds_time.to_netcdf(path=full_file_name)
        ds_time.close()
