import numpy as np
import xarray as xr

from itertools import combinations
import dask.array as dsa

from ..core import calc_cape
from .fixtures import empty_dask_array, dataset_soundings

import pytest

@pytest.fixture(scope='module')
def p_t_td_1d(nlevs=20):
    p = np.random.rand(nlevs)
    t = np.random.rand(nlevs)
    td = np.random.rand(nlevs)
    return p, t, td

@pytest.fixture(scope='module')
def p_t_td_3d(nlevs=20, nx=10, ny=5):
    p = np.random.rand(nlevs, ny, nx)
    t = np.random.rand(nlevs, ny, nx)
    td = np.random.rand(nlevs, ny, nx)
    return p, t, td

@pytest.fixture(scope='module')
def p_t_td_surface(nx=10, ny=5):
    ps = np.random.rand(ny, nx)
    ts = np.random.rand(ny, nx)
    tds = np.random.rand(ny, nx)
    return ps, ts, tds

# surface mode returns cape, cin
# most-unstable mode returns cape, cin, mulev, zmulev
@pytest.mark.parametrize('sourcein,n_returns',
                         [('surface', 2), ('most-unstable', 4)])
def test_calc_cape_shape_3d(p_t_td_3d, p_t_td_surface, sourcein, n_returns):
    p, t, td = p_t_td_3d
    ps, ts, tds = p_t_td_surface
    result = calc_cape(p, t, td, ps, ts, tds, source=sourcein, method='dummy')
    assert len(result) == n_returns

    for data in result:
        assert data.shape == (1, p.shape[1], p.shape[2])

# tolerance for tests
decimal = 3

def test_calc_surface_cape_model_lev(dataset_soundings):
    """Test Surface Cape based on previously calculated using George Bryans code"""
    ds = dataset_soundings

    cape, cin = calc_cape(ds.pressure.values[1:],
                          ds.temperature.values[1:],
                          ds.dewpoint.values[1:],
                          ds.pressure.values[0],
                          ds.temperature.values[0],
                          ds.dewpoint.values[0],
                          source='surface', ml_depth=500., adiabat='pseudo-liquid',
                          pinc=100.,
                          method='fortran', vertical_lev='sigma', pres_lev_pos=1)

    np.testing.assert_almost_equal(cape[0], ds.SB_CAPE_pinc100.values, decimal)
    np.testing.assert_almost_equal(cin[0], ds.SB_CIN_pinc100.values, decimal)

def test_calc_most_unstable_cape_model_lev(dataset_soundings):
    """Test Surface Cape based on previously calculated using George Bryans code"""
    ds = dataset_soundings

    cape, cin, mulv, zmulv = calc_cape(ds.pressure.values[1:],
                          ds.temperature.values[1:],
                          ds.dewpoint.values[1:],
                          ds.pressure.values[0],
                          ds.temperature.values[0],
                          ds.dewpoint.values[0],
                          source='most-unstable', ml_depth=500., adiabat='pseudo-liquid',
                          pinc=100.,
                          method='fortran', vertical_lev='sigma', pres_lev_pos=1)

    np.testing.assert_almost_equal(cape[0], ds.MU_CAPE_pinc100.values, decimal)
    np.testing.assert_almost_equal(cin[0], ds.MU_CIN_pinc100.values, decimal)
    np.testing.assert_almost_equal(mulv[0], ds.MU_lv_pinc100.values.astype('int32'), decimal)
    np.testing.assert_almost_equal(zmulv[0], ds.MU_z_pinc100.values, decimal)

def test_calc_mixed_layer_cape_model_lev(dataset_soundings):
    """Test Surface Cape based on previously calculated using George Bryans code"""
    ds = dataset_soundings

    cape, cin = calc_cape(ds.pressure.values[1:],
                          ds.temperature.values[1:],
                          ds.dewpoint.values[1:],
                          ds.pressure.values[0],
                          ds.temperature.values[0],
                          ds.dewpoint.values[0],
                          source='mixed-layer', ml_depth=500., adiabat='pseudo-liquid',
                          pinc=1000.,
                          method='fortran', vertical_lev='sigma', pres_lev_pos=1)

    np.testing.assert_almost_equal(cape[0], ds.ML_CAPE_pinc1000_mldepth500.values, decimal)
    np.testing.assert_almost_equal(cin[0], ds.ML_CIN_pinc1000_mldepth500.values, decimal)
