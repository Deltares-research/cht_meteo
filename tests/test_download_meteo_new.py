import shutil
import tempfile
from datetime import datetime, timedelta
from pathlib import Path

import pytest

import cht_meteo.cht.meteo as meteo

lon_range = [-80.0, -70.0]
lat_range = [30.0, 40.0]
time_range = [datetime(2024, 9, 26, 0, 0, 0), datetime(2024, 9, 26, 3, 0, 0)]


def test_gfs_anl_0p5(setup_temp_test_dir):
    meteo_name = "gfs_anl_0p50"

    # Create GFS dataset
    gfs_anl = meteo.dataset(
        name=meteo_name,
        source="gfs_analysis_0p50",
        path=setup_temp_test_dir,
        lon_range=lon_range,
        lat_range=lat_range,
    )

    gfs_anl.download(time_range)

    assert (
        setup_temp_test_dir / f"{meteo_name}.{time_range[0].strftime('%Y%m%d_%H%M')}.nc"
    ).is_file()

    gfs_anl.collect(time_range)

    assert list(gfs_anl.ds.data_vars.keys())[2] == "barometric_pressure"
    assert gfs_anl.ds["barometric_pressure"].dtype == "float64"

    del gfs_anl


def test_gfs_forecast_0p25(setup_temp_test_dir):
    meteo_name = "gfs_forecast_0p25"

    # Create GFS dataset
    gfs_for = meteo.dataset(
        name=meteo_name,
        source="gfs_forecast_0p25",
        path=setup_temp_test_dir,
        lon_range=lon_range,
        lat_range=lat_range,
    )

    # Get current time
    t_now = datetime.now()
    hours_to_add = (6 - t_now.hour % 6) % 6
    t0 = (t_now + timedelta(hours=hours_to_add)).replace(
        minute=0, second=0, microsecond=0
    )

    time_range = [t0, t0 + timedelta(hours=3)]

    gfs_for.download(time_range)

    gfs_for.collect(time_range)

    assert list(gfs_for.ds.data_vars.keys())[2] == "barometric_pressure"
    assert gfs_for.ds["barometric_pressure"].dtype == "float64"

    del gfs_for


def test_coamps_tc_forecast_s3(setup_temp_test_dir):
    meteo_name = "coamps_tc_forecast_s3"
    storm_number = "10L"

    # And now for download COAMPS-TC forecast data
    ctc = meteo.dataset(
        name=meteo_name,
        source="coamps_tc_forecast_s3",
        path=setup_temp_test_dir,
        lon_range=lon_range,
        lat_range=lat_range,
    )

    # Download the data for the specified time range
    # Since COAMPS-TC is a forecasts product, this should download several cycles

    # In addition to the time range, we also need to specify the storm number for COAMPS-TC forecasts
    ctc.download(time_range, storm_number=storm_number)

    # Collect the downloaded data
    ctc.collect(time_range)

    # Since ctc has 3 subsets, we need to interpolate onto a fixed grid
    # Create a new meteo dataset
    ctcf = ctc.cut_out(
        x_range=lon_range,
        y_range=lat_range,
        dx=0.1,
        dy=0.1,
    )
    # Write to netcdf
    filename = setup_temp_test_dir / f"{meteo_name}.nc"
    ctcf.to_netcdf(filename=filename)

    assert filename.is_file()
    assert list(ctcf.ds.data_vars.keys())[2] == "barometric_pressure"
    assert ctcf.ds["barometric_pressure"].dtype == "float64"

    # Write to delft3d
    filename = setup_temp_test_dir / f"{meteo_name}"
    ctcf.to_delft3d(file_name=str(filename))

    assert all(
        [Path(f"{filename}.{ext}").is_file() for ext in ["amp", "ampr", "amu", "amv"]]
    )

    del ctcf


@pytest.fixture()
def setup_temp_test_dir(request):
    test_name = request.node.name
    test_path = Path(tempfile.gettempdir()) / f"test_download_meteo_{test_name}"
    if test_path.exists():
        shutil.rmtree(test_path)
    test_path.mkdir(parents=True)

    yield test_path

    shutil.rmtree(test_path)
