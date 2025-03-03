import shutil
import tempfile
from datetime import datetime, timedelta
from pathlib import Path

import pytest

import cht_meteo


@pytest.fixture(scope="session")
def lon_range():
    return (-80.0, -70.0)


@pytest.fixture(scope="session")
def lat_range():
    return (30.0, 40.0)


@pytest.fixture(scope="session")
def time_range():
    return (datetime(2024, 9, 26, 0, 0, 0), datetime(2024, 9, 26, 3, 0, 0))


@pytest.fixture(scope="session")
def time_range_now():
    # Get current time
    t_now = datetime.now()
    hours_to_add = (6 - t_now.hour % 6) % 6
    t0 = (t_now + timedelta(hours=hours_to_add)).replace(
        minute=0, second=0, microsecond=0
    )

    time_range = (t0, t0 + timedelta(hours=3))
    return time_range


@pytest.fixture(scope="session")
def gfs_anl_dataset(request, lon_range, lat_range, time_range):
    temp_dir = Path(tempfile.gettempdir()) / "gfs_anl_data"
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    temp_dir.mkdir(parents=True)

    meteo_name = "gfs_anl_0p50"
    gfs_anl = cht_meteo.dataset(
        name=meteo_name,
        source="gfs_analysis_0p50",
        path=temp_dir,
        lon_range=lon_range,
        lat_range=lat_range,
    )
    gfs_anl.download(time_range)

    yield gfs_anl

    shutil.rmtree(temp_dir)


@pytest.fixture(scope="session")
def gfs_fc_dataset(request, lon_range, lat_range, time_range_now):
    temp_dir = Path(tempfile.gettempdir()) / "gfs_fc_data"
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    temp_dir.mkdir(parents=True)

    meteo_name = "gfs_forecast_0p25"
    # Create GFS dataset
    gfs_for = cht_meteo.dataset(
        name=meteo_name,
        source="gfs_forecast_0p25",
        path=temp_dir,
        lon_range=lon_range,
        lat_range=lat_range,
    )
    gfs_for.download(time_range_now)

    yield gfs_for

    shutil.rmtree(temp_dir)


@pytest.fixture(scope="session")
def coamps_tc_dataset(request, lon_range, lat_range, time_range):
    temp_dir = Path(tempfile.gettempdir()) / "coamps_tc_data"
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    temp_dir.mkdir(parents=True)

    meteo_name = "coamps_tc_forecast_s3"
    storm_number = "10L"
    ctc = cht_meteo.dataset(
        name=meteo_name,
        source="coamps_tc_forecast_s3",
        path=temp_dir,
        lon_range=lon_range,
        lat_range=lat_range,
    )
    ctc.download(time_range, storm_number=storm_number)

    yield ctc

    shutil.rmtree(temp_dir)
