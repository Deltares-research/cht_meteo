import shutil
import tempfile
import os
from datetime import datetime
from pathlib import Path

# import pytest
from pyproj import CRS

import cht_meteo.cht.meteo as meteo

meteo_path = "c:/work/temp/test_meteo_database"

time_range = [datetime(2024, 9, 26, 0, 0, 0), datetime(2024, 9, 27, 0, 0, 0)]
lon_range = [-80.0, -70.0]
lat_range = [30.0, 40.0]
storm_number = "10L"

# Create GFS dataset
gfs_anl = meteo.dataset(name="gfs_anl_0p50_atlantic",
                        source="gfs_analysis_0p50",                      
                        path=os.path.join(meteo_path, "gfs_anl_0p50_atlantic"),
                        lon_range=lon_range,
                        lat_range=lat_range,)
    
# Alternatively, we could use the following:
# from cht_meteo.cht.meteo.gfs_anl_0p50 import MeteoDatasetGFSAnalysis0p50
# d = meteo.MeteoDatasetGFSAnalysis0p50(name="gfs_anl_0p50_atlantic",
#                                       path=os.path.join(meteo_path, "gfs_anl_0p50_atlantic"),
#                                       lon_range=[-80.0, -70.0],
#                                       lat_range=[30.0, 40.0],)

# Download the data for the specified time range    
gfs_anl.download(time_range)
gfs_anl.collect(time_range)
    
# And now for download COAMPS-TC forecast data
ctc = meteo.dataset(name="coamps_tc_forecast_atlantic",
                    source="coamps_tc_forecast_s3",
                    path=os.path.join(meteo_path, "coamps_tc_forecast_atlantic"),
                    lon_range=lon_range,
                    lat_range=lat_range,)

# Download the data for the specified time range
# Since COAMPS-TC is a forecasts product, this should download several cycles

# In addition to the time range, we also need to specify the storm number for COAMPS-TC forecasts
ctc.download(time_range, storm_number=storm_number)

# Collect the downloaded data
ctc.collect(time_range)

# Since ctc has 3 subsets, we need to interpolate onto a fixed grid
# Create a new meteo dataset
ctcf = ctc.cut_out(x_range=lon_range,
                   y_range=lat_range,
                   dx=0.1, dy=0.1,
                  )

# Fill gaps with GFS data
ctcf.merge_datasets([gfs_anl])

# Write to netcdf and delft3d
ctcf.to_netcdf(filename=os.path.join(meteo_path, "coamps_tc_forecast_atlantic", "coamps_tc_forecast_atlantic.nc"))
ctcf.to_delft3d(os.path.join(meteo_path, "coamps_tc_forecast_atlantic", "coamps_tc_forecast_atlantic"))

#     gfs_conus.download(time_range)
#     assert (
#         setup_temp_test_dir / "gfs_anl_0p50_us_southeast.20230101_0000.nc"
#     ).is_file()

#     gfs_conus.collect(time_range)

#     assert gfs_conus.quantity[1].name == "barometric_pressure"
#     assert gfs_conus.quantity[0].u.dtype == "float64"

#     del gfs_conus, gfs_source


# @pytest.fixture()
# def setup_temp_test_dir():
#     test_path = Path(tempfile.gettempdir()) / "test_download_meteo"
#     if test_path.exists():
#         shutil.rmtree(test_path)
#     test_path.mkdir(parents=True)

#     yield test_path

#     shutil.rmtree(test_path)
