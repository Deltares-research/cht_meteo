import numpy as np


def test_coamps_tc_forecast_cutout(coamps_tc_dataset, time_range, lon_range, lat_range):
    # Collect the data
    coamps_tc_dataset.collect(time_range)
    # Since ctc has 3 subsets, we need to interpolate onto a fixed grid
    # Create a new meteo dataset
    dx = 0.1
    dy = 0.1
    dataset = coamps_tc_dataset.cut_out(
        dx=dx, dy=dy, lon_range=lon_range, lat_range=lat_range
    )

    # Check if the data was collected correctly
    assert list(dataset.ds.data_vars.keys())[2] == "barometric_pressure"
    assert dataset.ds["barometric_pressure"].dtype == "float64"
    assert np.isclose(
        dataset.ds["lon"][1].to_numpy() - dataset.ds["lon"][0].to_numpy(), dx, atol=1e-6
    )
    assert np.isclose(
        dataset.ds["lat"][1].to_numpy() - dataset.ds["lat"][0].to_numpy(), dy, atol=1e-6
    )


def test_forecast_merge(
    coamps_tc_dataset, gfs_anl_dataset, time_range, lon_range, lat_range
):
    # Collect the data
    coamps_tc_dataset.collect(time_range)
    gfs_anl_dataset.collect(time_range)
    # Since ctc has 3 subsets, we need to interpolate onto a fixed grid
    # Create a new meteo dataset
    dx = 0.1
    dy = 0.1
    dataset = coamps_tc_dataset.cut_out(
        dx=dx, dy=dy, lon_range=lon_range, lat_range=lat_range
    )

    # Force some values to NaN
    dataset.ds = dataset.ds.where(dataset.ds["lat"] < 36, np.nan)
    # Now use GFS to fill the nan values
    dataset.merge_datasets([gfs_anl_dataset])
    # Check if the data was filled correctly
    for var in dataset.ds.data_vars:
        assert not np.isnan(
            dataset.ds[var].to_numpy()
        ).any(), f"NaN values found in variable {var}"
