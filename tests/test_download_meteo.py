# import os
from datetime import datetime
from pathlib import Path

from pyproj import CRS

from cht_meteo.meteo import MeteoGrid, MeteoSource

path = Path().absolute() / "tests"


def test_download_meteo():
    params = ["wind", "barometric_pressure", "precipitation"]
    lat = 32.77
    lon = -79.95

    # Download the actual datasets
    gfs_source = MeteoSource("gfs_anl_0p50", "gfs_anl_0p50_04", "hindcast", delay=None)

    # Create subset
    name = "gfs_anl_0p50_us_southeast"
    gfs_conus = MeteoGrid(
        name=name,
        source=gfs_source,
        parameters=params,
        path=path,
        x_range=[lon - 1, lon + 1],
        y_range=[lat - 1, lat + 1],
        crs=CRS.from_epsg(4326),
    )

    # Download and collect data
    t0 = datetime.strptime("20230101 000000", "%Y%m%d %H%M%S")
    t1 = datetime.strptime("20230101 020000", "%Y%m%d %H%M%S")
    time_range = [t0, t1]

    gfs_conus.download(time_range)
    assert path.joinpath("gfs_anl_0p50_us_southeast.20230101_0000.nc").is_file()

    gfs_conus.collect(time_range)

    assert gfs_conus.quantity[1].name == "barometric_pressure"
    assert gfs_conus.quantity[0].u.dtype == "float64"

    # Remove downloaded file
    # os.remove(path.joinpath("gfs_anl_0p50_us_southeast.20230101_0000.nc"))
