import shutil
import tempfile
from pathlib import Path

import pytest
import xarray as xr


@pytest.fixture
def coamps_tc_dataset_cutout(coamps_tc_dataset, time_range, lon_range, lat_range):
    dx = 0.1
    dy = 0.1
    coamps_tc_dataset.collect(time_range)
    dataset = coamps_tc_dataset.cut_out(
        dx=dx, dy=dy, lon_range=lon_range, lat_range=lat_range
    )
    return dataset


@pytest.mark.parametrize(
    "dataset_fixture", ["coamps_tc_dataset_cutout", "gfs_anl_dataset", "gfs_fc_dataset"]
)
def test_save_netcdf(request, dataset_fixture, time_range, time_range_now):
    dataset = request.getfixturevalue(dataset_fixture)
    if dataset_fixture == "gfs_anl_dataset":
        dataset.collect(time_range)
    elif dataset_fixture == "gfs_fc_dataset":
        dataset.collect(time_range_now)
    temp_dir = Path(tempfile.gettempdir()) / "save_netcdf"
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    temp_dir.mkdir(parents=True)
    file = temp_dir / "dataset.nc"
    dataset.to_netcdf(file_name=file)
    assert file.is_file()
    assert xr.open_dataset(file).equals(dataset.ds)
    shutil.rmtree(temp_dir)


@pytest.mark.parametrize(
    "dataset_fixture", ["coamps_tc_dataset_cutout", "gfs_anl_dataset", "gfs_fc_dataset"]
)
def test_save_delft3d(request, dataset_fixture, time_range, time_range_now):
    dataset = request.getfixturevalue(dataset_fixture)
    if dataset_fixture == "gfs_anl_dataset":
        dataset.collect(time_range)
    elif dataset_fixture == "gfs_fc_dataset":
        dataset.collect(time_range_now)
    temp_dir = Path(tempfile.gettempdir()) / "save_test_coamps_tc_data_delft3d"
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    temp_dir.mkdir(parents=True)
    file = temp_dir / "dataset"
    dataset.to_delft3d(file_name=str(file))

    assert all(
        [Path(f"{file}.{ext}").is_file() for ext in ["amp", "ampr", "amu", "amv"]]
    )
    shutil.rmtree(temp_dir)
