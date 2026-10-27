"""Tests for longitude handling: ranges west of Greenwich and across the dateline."""

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from cht_meteo.dataset import (
    MeteoDataset,
    is_global_lon,
    lon_extent,
    lon_window,
    normalize_lon_range,
)

TIMES = pd.date_range("2026-01-01", periods=2, freq="3h")


def pressure(lon, lat):
    """Smooth field that is continuous across the dateline."""
    lon = np.asarray(lon)
    lat = np.asarray(lat)
    return 101000.0 + 500.0 * np.cos(np.deg2rad(lon)) + 100.0 * np.sin(np.deg2rad(lat))


def make_dataset(lon):
    """MeteoDataset with a pressure field on the given longitudes."""
    lat = np.arange(-30.0, 30.01, 0.25)
    lon2, lat2 = np.meshgrid(lon, lat)
    values = np.repeat(pressure(lon2, lat2)[None, :, :], len(TIMES), axis=0)
    dataset = MeteoDataset()
    dataset.ds = xr.Dataset(
        {"barometric_pressure": (("time", "lat", "lon"), values)},
        coords={"time": TIMES, "lat": lat, "lon": lon},
    )
    return dataset


GLOBAL_0_360 = np.arange(0.0, 360.0, 0.25)
GLOBAL_180 = np.arange(-180.0, 180.0, 0.25)


def test_normalize_lon_range():
    assert normalize_lon_range([170, -170]) == [170.0, 190.0]
    assert normalize_lon_range([-10, 10]) == [-10.0, 10.0]
    assert normalize_lon_range([350, 10]) == [350.0, 370.0]
    assert normalize_lon_range([170, 190]) == [170.0, 190.0]


def test_lon_extent():
    assert lon_extent(np.array([-179.5, 179.5])) == [179.5, 180.5]
    assert lon_extent(np.array([-10.0, 10.0])) == [-10.0, 10.0]
    assert lon_extent(np.array([100.0, 120.0])) == [100.0, 120.0]
    assert lon_extent(np.array([350.0, 5.0])) == [-10.0, 5.0]


def test_is_global_lon():
    assert is_global_lon(GLOBAL_0_360)
    assert is_global_lon(GLOBAL_180)
    assert not is_global_lon(np.arange(100.0, 120.0, 0.25))


@pytest.mark.parametrize("source_lon", [GLOBAL_0_360, GLOBAL_180], ids=["0-360", "180"])
@pytest.mark.parametrize(
    "lon_range, west, east",
    [
        ([-10.0, 10.0], -10.0, 10.0),  # west of Greenwich
        ([100.0, 120.0], 100.0, 120.0),  # plain range
        ([170.0, -170.0], 170.0, 190.0),  # across the dateline
        ([170.0, 190.0], 170.0, 190.0),  # across the dateline, continuous
        ([-190.0, -170.0], -190.0, -170.0),  # across the dateline, negative
    ],
)
def test_cut_out_lon_range(source_lon, lon_range, west, east):
    out = make_dataset(source_lon).cut_out(lon_range=lon_range, lat_range=[-5.0, 5.0])
    assert out is not None
    lon = out.ds["lon"].to_numpy()
    lat = out.ds["lat"].to_numpy()
    # Continuous, ascending longitudes covering the requested range
    assert np.all(np.diff(lon) > 0)
    assert lon[0] <= west and lon[-1] >= east
    assert lon[-1] - lon[0] < 30.0
    # Values match the field at the (shifted) longitudes
    lon2, lat2 = np.meshgrid(lon, lat)
    values = out.ds["barometric_pressure"].isel(time=0).to_numpy()
    np.testing.assert_allclose(values, pressure(lon2, lat2), atol=1e-6)


def test_cut_out_dx_dy_across_dateline():
    out = make_dataset(GLOBAL_0_360).cut_out(
        lon_range=[170.0, -170.0], lat_range=[-5.0, 5.0], dx=0.5, dy=0.5
    )
    lon = out.ds["lon"].to_numpy()
    lat = out.ds["lat"].to_numpy()
    assert np.all(np.diff(lon) > 0)
    assert lon[0] <= 170.0 and lon[-1] >= 190.0
    lon2, lat2 = np.meshgrid(lon, lat)
    values = out.ds["barometric_pressure"].isel(time=0).to_numpy()
    np.testing.assert_allclose(values, pressure(lon2, lat2), atol=0.5)


def test_cut_out_regional_source_shifted():
    # Regional source stored 0-360, requested in negative longitudes
    source = make_dataset(np.arange(160.0, 200.0, 0.25))
    out = source.cut_out(lon_range=[-190.0, -170.0], lat_range=[-5.0, 5.0])
    lon = out.ds["lon"].to_numpy()
    assert lon[0] <= -190.0 and lon[-1] >= -170.0
    values = out.ds["barometric_pressure"].isel(time=0, lat=0).to_numpy()
    lat0 = out.ds["lat"].to_numpy()[0]
    np.testing.assert_allclose(values, pressure(lon, lat0), atol=1e-6)


@pytest.mark.parametrize("source_lon", [GLOBAL_0_360, GLOBAL_180], ids=["0-360", "180"])
def test_interpolate_variable_any_longitude(source_lon):
    dataset = make_dataset(source_lon)
    lon = np.array([-179.5, 179.5, 180.5, 359.9, -0.1, 0.1, 10.0])
    lat = np.full_like(lon, 1.0)
    values = dataset.interpolate_variable(
        "barometric_pressure", TIMES[0], lon, lat, crs=dataset.crs
    )
    assert not np.any(np.isnan(values))
    np.testing.assert_allclose(values, pressure(lon, lat), atol=0.5)


def test_lon_window_drops_duplicate_seam():
    # Grid that includes both 0 and 360
    ds = make_dataset(np.arange(0.0, 360.01, 0.25)).ds
    out = lon_window(ds, -5.0, 5.0)
    assert np.all(np.diff(out["lon"].to_numpy()) > 0)


@pytest.mark.parametrize(
    "server_lon, lon_range, expected",
    [
        # The GFS THREDDS server returns 359..361 for a request of -1..1
        (np.arange(359.0, 361.01, 0.25), [-1.0, 1.0], (-1.0, 1.0)),
        (np.arange(179.0, 181.01, 0.25), [179.0, -179.0], (179.0, 181.0)),
        (np.arange(170.0, 190.01, 0.25), [-190.0, -170.0], (-190.0, -170.0)),
    ],
)
def test_clip_lon_regional_server_response(server_lon, lon_range, expected):
    reader = MeteoDataset(lon_range=lon_range)
    out = reader.clip_lon(make_dataset(server_lon).ds)
    lon = out["lon"].to_numpy()
    assert lon[0] == pytest.approx(expected[0])
    assert lon[-1] == pytest.approx(expected[1])


@pytest.mark.parametrize("lon_range", [[170.0, -170.0], [-10.0, 10.0], [350.0, 370.0]])
def test_clip_lon_global_dataarray(lon_range):
    # As used by the readers on (lazily loaded) global 0-360 GFS fields
    reader = MeteoDataset(lon_range=lon_range)
    da = make_dataset(GLOBAL_0_360).ds["barometric_pressure"].isel(time=0)
    out = reader.clip_lon(da)
    lon = out["lon"].to_numpy()
    west, east = normalize_lon_range(lon_range)
    assert np.all(np.diff(lon) > 0)
    assert lon[0] <= west and lon[-1] >= east and lon[-1] - lon[0] < 30.0
    lon2, lat2 = np.meshgrid(lon, out["lat"].to_numpy())
    np.testing.assert_allclose(out.to_numpy(), pressure(lon2, lat2), atol=1e-6)


def test_clip_lon_with_other_coordinate_name():
    reader = MeteoDataset(lon_range=[170.0, -170.0])
    ds = make_dataset(GLOBAL_180).ds.rename({"lon": "longitude"})
    out = reader.clip_lon(ds, lon_name="longitude")
    lon = out["longitude"].to_numpy()
    assert lon[0] <= 170.0 and lon[-1] >= 190.0


def test_lon_window_selects_only_needed_columns():
    # Global data is cut before joining, so the result is small (not 3 x global)
    ds = make_dataset(GLOBAL_0_360).ds
    out = lon_window(ds, 170.0, 190.0)
    assert out.sizes["lon"] <= (190.0 - 170.0) / 0.25 + 1 + 2 * 2
