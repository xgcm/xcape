import numpy as np
import xarray as xr

from itertools import combinations
import dask.array as dsa
import dask

from ..core import calc_cape
from ..core import calc_srh
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
    p = np.random.rand(ny, nx, nlevs)
    t = np.random.rand(ny, nx, nlevs)
    td = np.random.rand(ny, nx, nlevs)
    return p, t, td

@pytest.fixture(scope='module')
def p_t_td_surface(nx=10, ny=5):
    ps = np.random.rand(ny, nx)
    ts = np.random.rand(ny, nx)
    tds = np.random.rand(ny, nx)
    return ps, ts, tds

# surface mode returns cape, cin
# most-unstable mode returns cape, cin, mulev, zmulev
@pytest.mark.parametrize('use_dask', [False, True])
@pytest.mark.parametrize('sourcein,n_returns',
                         [('surface', 2), ('most-unstable', 4)])
def test_calc_cape_shape_3d(p_t_td_3d, p_t_td_surface, sourcein, n_returns, use_dask):
    p, t, td = p_t_td_3d
    ps, ts, tds = p_t_td_surface
    args = (p, t, td, ps, ts, tds)
    if use_dask:
        args = [dsa.from_array(a) for a in args]
    result = calc_cape(*args, source=sourcein, method='dummy')
    assert len(result) == n_returns
    for data in result:
        assert data.shape == p.shape[:-1]
        if use_dask:
            assert isinstance(data, dsa.Array)
            data.compute()
    

# tolerance for tests
decimal_cape = 0
decimal_cin = 0
decimal_mulv = 0
decimal_zmulv = 0

@pytest.mark.parametrize('use_dask', [False, True])
def test_calc_surface_cape_model_lev(dataset_soundings, use_dask):
    """Test Surface Cape based on previously calculated using George Bryans code"""
    ds = dataset_soundings
    if use_dask:
        ds = ds.chunk()

    cape, cin = calc_cape(ds.pressure.data[:, 1:],
                          ds.temperature.data[:, 1:],
                          ds.dewpoint.data[:, 1:],
                          ds.pressure.data[:, 0],
                          ds.temperature.data[:, 0],
                          ds.dewpoint.data[:, 0],
                          source='surface', ml_depth=500., adiabat='pseudo-liquid',
                          pinc=100.,
                          method='fortran', vertical_lev='sigma', pres_lev_pos=1)
    if use_dask:
        assert isinstance(cape, dsa.Array)
        assert isinstance(cin, dsa.Array)
        cape, cin = dask.compute(cape, cin)

    np.testing.assert_almost_equal(cape, ds.SB_CAPE_pinc100.values, decimal_cape)
    np.testing.assert_almost_equal(cin, ds.SB_CIN_pinc100.values, decimal_cin)

def test_calc_most_unstable_cape_model_lev(dataset_soundings):
    """Test Surface Cape based on previously calculated using George Bryans code"""
    ds = dataset_soundings

    # in real data, the surface values will come in separate variables
    cape, cin, mulv, zmulv = calc_cape(ds.pressure.values[:, 1:],
                          ds.temperature.values[:, 1:],
                          ds.dewpoint.values[:, 1:],
                          ds.pressure.values[:, 0],
                          ds.temperature.values[:, 0],
                          ds.dewpoint.values[:, 0],
                          source='most-unstable', ml_depth=500., adiabat='pseudo-liquid',
                          pinc=100.,
                          method='fortran', vertical_lev='sigma', pres_lev_pos=1)

    np.testing.assert_almost_equal(cape, ds.MU_CAPE_pinc100.values, decimal_cape)
    np.testing.assert_almost_equal(cin, ds.MU_CIN_pinc100.values, decimal_cin)
    np.testing.assert_almost_equal(mulv, ds.MU_lv_pinc100.values.astype('int32'), decimal_mulv)
    np.testing.assert_almost_equal(zmulv, ds.MU_z_pinc100.values, decimal_zmulv)

def test_calc_mixed_layer_cape_model_lev(dataset_soundings):
    """Test Surface Cape based on previously calculated using George Bryans code"""
    ds = dataset_soundings

    cape, cin = calc_cape(ds.pressure.values[:, 1:],
                          ds.temperature.values[:, 1:],
                          ds.dewpoint.values[:, 1:],
                          ds.pressure.values[:, 0],
                          ds.temperature.values[:, 0],
                          ds.dewpoint.values[:, 0],
                          source='mixed-layer', ml_depth=500., adiabat='pseudo-liquid',
                          pinc=1000.,
                          method='fortran', vertical_lev='sigma', pres_lev_pos=1)

    np.testing.assert_almost_equal(cape, ds.ML_CAPE_pinc1000_mldepth500.values, decimal_cape)
    np.testing.assert_almost_equal(cin, ds.ML_CIN_pinc1000_mldepth500.values, decimal_cin)

def test_calc_surface_cape_pressure_lev(dataset_soundings):
    """Test Surface Cape based on previously calculated using George Bryans code"""
    ds = dataset_soundings
    
    cape, cin = calc_cape(ds.pressure.values[:, 1:],
                          ds.temperature.values[:, 1:],
                          ds.dewpoint.values[:, 1:],
                          ds.pressure.values[:, 0],
                          ds.temperature.values[:, 0],
                          ds.dewpoint.values[:, 0],
                          source='surface', ml_depth=500., adiabat='pseudo-liquid',
                          pinc=100.,
                          method='fortran', vertical_lev='pressure', 
                          pres_lev_pos=ds.pressure.values[:,0]*0+1)

    np.testing.assert_almost_equal(cape, ds.SB_CAPE_pinc100.values, decimal_cape)
    np.testing.assert_almost_equal(cin, ds.SB_CIN_pinc100.values, decimal_cin)

def test_calc_most_unstable_cape_pressure_lev(dataset_soundings):
    """Test Surface Cape based on previously calculated using George Bryans code"""
    ds = dataset_soundings

    # in real data, the surface values will come in separate variables
    cape, cin, mulv, zmulv = calc_cape(ds.pressure.values[:, 1:],
                          ds.temperature.values[:, 1:],
                          ds.dewpoint.values[:, 1:],
                          ds.pressure.values[:, 0],
                          ds.temperature.values[:, 0],
                          ds.dewpoint.values[:, 0],
                          source='most-unstable', ml_depth=500., adiabat='pseudo-liquid',
                          pinc=100.,
                          method='fortran', vertical_lev='pressure', 
                          pres_lev_pos=ds.pressure.values[:,0]*0+1)

    np.testing.assert_almost_equal(cape, ds.MU_CAPE_pinc100.values, decimal_cape)
    np.testing.assert_almost_equal(cin, ds.MU_CIN_pinc100.values, decimal_cin)
    np.testing.assert_almost_equal(mulv, ds.MU_lv_pinc100.values.astype('int32'), decimal_mulv)
    np.testing.assert_almost_equal(zmulv, ds.MU_z_pinc100.values, decimal_zmulv)

def test_calc_mixed_layer_cape_pressure_lev(dataset_soundings):
    """Test Surface Cape based on previously calculated using George Bryans code"""
    ds = dataset_soundings

    cape, cin = calc_cape(ds.pressure.values[:, 1:],
                          ds.temperature.values[:, 1:],
                          ds.dewpoint.values[:, 1:],
                          ds.pressure.values[:, 0],
                          ds.temperature.values[:, 0],
                          ds.dewpoint.values[:, 0],
                          source='mixed-layer', ml_depth=500., adiabat='pseudo-liquid',
                          pinc=1000.,
                          method='fortran', vertical_lev='pressure', 
                          pres_lev_pos=ds.pressure.values[:, 0]*0+1)

    np.testing.assert_almost_equal(cape, ds.ML_CAPE_pinc1000_mldepth500.values, decimal_cape)
    np.testing.assert_almost_equal(cin, ds.ML_CIN_pinc1000_mldepth500.values, decimal_cin)

    
def test_calc_srh_model_lev(dataset_soundings):
    """Test SRH code"""
    ds = dataset_soundings

    srh, rm, lm, mean_6km = calc_srh(ds.pressure.values[:,1:],
                          ds.temperature.values[:,1:],
                          ds.dewpoint.values[:,1:],
                          ds.u_wind_ms.values[:,1:],
                          ds.v_wind_ms.values[:,1:],                         
                          ds.pressure.values[:,0],
                          ds.temperature.values[:,0],
                          ds.dewpoint.values[:,0],
                          ds.u_wind_ms.values[:,0],
                          ds.v_wind_ms.values[:,0],
                          depth = 3000,
                          vertical_lev='sigma', pres_lev_pos=1,
                          output_var='all')
    srh2 = calc_srh(ds.pressure.values[:,1:],
                          ds.temperature.values[:,1:],
                          ds.dewpoint.values[:,1:],
                          ds.u_wind_ms.values[:,1:],
                          ds.v_wind_ms.values[:,1:],                         
                          ds.pressure.values[:,0],
                          ds.temperature.values[:,0],
                          ds.dewpoint.values[:,0],
                          ds.u_wind_ms.values[:,0],
                          ds.v_wind_ms.values[:,0],
                          depth = 3000,
                          vertical_lev='sigma', pres_lev_pos=1,
                          output_var='srh')
    np.testing.assert_almost_equal(srh, ds.SRH03_model_lev.values, 5)
    np.testing.assert_almost_equal(srh2, ds.SRH03_model_lev.values, 5)
    
def test_calc_srh_pressure_lev(dataset_soundings):
    """Test SRH code"""
    ds = dataset_soundings

    srh, rm, lm, mean_6km = calc_srh(ds.pressure.values[:,1:],
                          ds.temperature.values[:,1:],
                          ds.dewpoint.values[:,1:],
                          ds.u_wind_ms.values[:,1:],
                          ds.v_wind_ms.values[:,1:],                         
                          ds.pressure.values[:,0],
                          ds.temperature.values[:,0],
                          ds.dewpoint.values[:,0],
                          ds.u_wind_ms.values[:,0],
                          ds.v_wind_ms.values[:,0],
                          depth = 3000,
                          vertical_lev='pressure', 
                          pres_lev_pos=ds.pressure.values[:,0]*0+1,
                          output_var='all')
    srh2 = calc_srh(ds.pressure.values[:,1:],
                          ds.temperature.values[:,1:],
                          ds.dewpoint.values[:,1:],
                          ds.u_wind_ms.values[:,1:],
                          ds.v_wind_ms.values[:,1:],                         
                          ds.pressure.values[:,0],
                          ds.temperature.values[:,0],
                          ds.dewpoint.values[:,0],
                          ds.u_wind_ms.values[:,0],
                          ds.v_wind_ms.values[:,0],
                          depth = 3000,
                          vertical_lev='pressure',
                          pres_lev_pos=ds.pressure.values[:,0]*0+1,
                          output_var='srh')
    np.testing.assert_almost_equal(srh, ds.SRH03_pressure_lev.values, 5)
    np.testing.assert_almost_equal(srh2, ds.SRH03_pressure_lev.values, 5)
