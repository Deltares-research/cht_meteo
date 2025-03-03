def test_gfs_anl_0p50_download(gfs_anl_dataset, time_range):
    # Check if the GFS analysis data file exists
    meteo_name = "gfs_anl_0p50"
    expected_file = (
        gfs_anl_dataset.path
        / f"{meteo_name}.{time_range[0].strftime('%Y%m%d_%H%M')}.nc"
    )
    assert expected_file.is_file()


# @pytest.mark.skip(reason="issue with time variable when downloading forecast.")
def test_gfs_fc_0p25_download(gfs_fc_dataset, time_range_now):
    # Check if the GFS analysis data file exists
    meteo_name = "gfs_forecast_0p25"
    cycle_dir = list(gfs_fc_dataset.path.glob("*z"))[0]
    expected_file = (
        cycle_dir / f"{meteo_name}.{time_range_now[0].strftime('%Y%m%d_%H%M')}.nc"
    )
    assert expected_file.is_file()


def test_coamps_tc_forecast_download(coamps_tc_dataset, time_range):
    # Check if the COAMPS-TC forecast data file exists
    meteo_name = "coamps_tc_forecast_s3"
    # Get download cycle name and check if folder exists
    cycle_name = time_range[0].strftime("%Y%m%d_%Hz")
    expected_cycle_folder = coamps_tc_dataset.path / cycle_name
    assert expected_cycle_folder.exists()
    # Check if the different levels of forecasts exist for the first time step
    expected_files = [
        expected_cycle_folder
        / f"{meteo_name}.{level}.{time_range[0].strftime('%Y%m%d_%H%M')}.nc"
        for level in ["d01", "d02", "d03"]
    ]
    assert all([file.is_file() for file in expected_files])
    # Check if the track file was downloaded
    track_file = list(expected_cycle_folder.glob("*.trk"))
    assert len(track_file) == 1
