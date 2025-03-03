def test_gfs_anl_0p50_collect(gfs_anl_dataset, time_range):
    dataset = gfs_anl_dataset
    # Collect the data
    dataset.collect(time_range)

    # Check if the data was collected correctly
    assert list(dataset.ds.data_vars.keys())[2] == "barometric_pressure"
    assert dataset.ds["barometric_pressure"].dtype == "float64"


# @pytest.mark.skip(reason="issue with time variable when downloading forecast.")
def test_gfs_fc_0p25_collect(gfs_fc_dataset, time_range_now):
    dataset = gfs_fc_dataset
    # Collect the data
    dataset.collect(time_range_now)

    # Check if the data was collected correctly
    assert list(dataset.ds.data_vars.keys())[2] == "barometric_pressure"
    assert dataset.ds["barometric_pressure"].dtype == "float32"


def test_coamps_tc_forecast_collect(coamps_tc_dataset, time_range):
    dataset = coamps_tc_dataset
    # Collect the data
    dataset.collect(time_range)

    # Check if the data was collected correctly
    assert len(dataset.subset) == 3
    assert list(dataset.subset[0].ds.data_vars.keys())[2] == "barometric_pressure"
    assert dataset.subset[0].ds["barometric_pressure"].dtype == "float64"
