import numpy as np
import xarray as xr

from itertools import combinations
import dask.array as dsa
import dask

from ..core import calc_cape
from ..core import calc_srh
from .fixtures import empty_dask_array, dataset_soundings, dataset_ERA5pressurelevel

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
@pytest.mark.parametrize('vertical_levin', ['sigma', 'pressure'])
def test_calc_cape_shape_3d(p_t_td_3d, p_t_td_surface, sourcein, n_returns, use_dask,vertical_levin):
    p, t, td = p_t_td_3d
    ps, ts, tds = p_t_td_surface
    if vertical_levin=='sigma':
        args = (p, t, td, ps, ts, tds)
    elif vertical_levin =='pressure':        
        args = (p, t, td, ps, ts, tds,ps*0+1)
    if use_dask:
        args = [dsa.from_array(a) for a in args]
    result = calc_cape(*args, source=sourcein,vertical_lev=vertical_levin, method='dummy')
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
@pytest.mark.parametrize('sourcein, pinc_used',
                         [('surface', 100),('mixed-layer',1000),('most-unstable', 100)])
@pytest.mark.parametrize('vertical_levin', ['sigma', 'pressure'])
def test_calc_cape(dataset_soundings, sourcein, pinc_used, use_dask,vertical_levin):
    """Test Surface Cape based on previously calculated using George Bryans code"""
    ds = dataset_soundings
    if use_dask:
        ds = ds.chunk()
    if vertical_levin=='sigma':
        returns = calc_cape(ds.pressure.data[:, 1:],
                              ds.temperature.data[:, 1:],
                              ds.dewpoint.data[:, 1:],
                              ds.pressure.data[:, 0],
                              ds.temperature.data[:, 0],
                              ds.dewpoint.data[:, 0],
                              source=sourcein, ml_depth=500., adiabat='pseudo-liquid',
                              pinc=pinc_used,
                              method='fortran', vertical_lev=vertical_levin)
    elif  vertical_levin=='pressure':
        returns = calc_cape(ds.pressure.data[:, 1:],
                              ds.temperature.data[:, 1:],
                              ds.dewpoint.data[:, 1:],
                              ds.pressure.data[:, 0],
                              ds.temperature.data[:, 0],
                              ds.dewpoint.data[:, 0],
                              ds.pressure.data[:,0]*0, #pres_lev_pos
                              source=sourcein, ml_depth=500., adiabat='pseudo-liquid',
                              pinc=pinc_used,
                              method='fortran', vertical_lev=vertical_levin)
    
    if sourcein=='most-unstable':
        cape = returns[0]
        cin = returns[1]
        mulv = returns[2] 
        zmulv = returns[3] 
        if use_dask:
            assert isinstance(cape, dsa.Array)
            assert isinstance(cin, dsa.Array)
            assert isinstance(mulv, dsa.Array)
            assert isinstance(zmulv, dsa.Array)
            cape, cin, mulv, zmulv = dask.compute(cape, cin, mulv, zmulv)

        np.testing.assert_almost_equal(cape, ds.MU_CAPE_pinc100.values, decimal_cape)
        np.testing.assert_almost_equal(cin, ds.MU_CIN_pinc100.values, decimal_cin)
        np.testing.assert_almost_equal(mulv, ds.MU_lv_pinc100.values.astype('int32'), decimal_mulv)
        np.testing.assert_almost_equal(zmulv, ds.MU_z_pinc100.values, decimal_zmulv)
    else:
        cape = returns[0]
        cin = returns[1]
        if use_dask:
            assert isinstance(cape, dsa.Array)
            assert isinstance(cin, dsa.Array)
            cape, cin = dask.compute(cape, cin)
            print(type(cape))
        if sourcein=='surface':
            np.testing.assert_almost_equal(cape, ds.SB_CAPE_pinc100.values, decimal_cape)
            np.testing.assert_almost_equal(cin, ds.SB_CIN_pinc100.values, decimal_cin)
        elif sourcein=='mixed-layer':
            np.testing.assert_almost_equal(cape, ds.ML_CAPE_pinc1000_mldepth500.values, decimal_cape)
            np.testing.assert_almost_equal(cin, ds.ML_CIN_pinc1000_mldepth500.values, decimal_cin)

@pytest.mark.parametrize('use_dask', [False, True])
@pytest.mark.parametrize('sourcein', ['surface', 'mixed-layer', 'most-unstable'])
@pytest.mark.parametrize('vertical_levin', ['pressure'])
@pytest.mark.parametrize('pinc_used', [500])
@pytest.mark.parametrize('ml_depthin', [300])
def test_calc_capepressure1d(dataset_ERA5pressurelevel, sourcein, pinc_used, ml_depthin, use_dask,vertical_levin):
    """Test Surface Cape based on previously calculated using George Bryans code"""
    ds3d, dssurf = dataset_ERA5pressurelevel
    if use_dask:
        ds3d = ds3d.chunk()
        dssurf = dssurf.chunk()
    returns = calc_cape(ds3d.level.data,
                          ds3d.t.data,
                          ds3d.td.data,
                          dssurf.p.data,
                          dssurf.t.data,
                          dssurf.td.data,
                          dssurf.start3d.data, #pres_lev_pos
                          source=sourcein, ml_depth=ml_depthin, adiabat='pseudo-liquid',
                          pinc=pinc_used,
                          method='fortran', 
                        vertical_lev=vertical_levin)
    
    if sourcein=='most-unstable':
        cape = returns[0]
        cin = returns[1]
        mulv = returns[2] 
        zmulv = returns[3] 
        if use_dask:
            assert isinstance(cape, dsa.Array)
            assert isinstance(cin, dsa.Array)
            assert isinstance(mulv, dsa.Array)
            assert isinstance(zmulv, dsa.Array)
            cape, cin, mulv, zmulv = dask.compute(cape, cin, mulv, zmulv)

        np.testing.assert_almost_equal(cape, dssurf.capemup500.values, decimal_cape)
        np.testing.assert_almost_equal(cin, dssurf.cinmup500.values, decimal_cin)
    else:
        cape = returns[0]
        cin = returns[1]
        if use_dask:
            assert isinstance(cape, dsa.Array)
            assert isinstance(cin, dsa.Array)
            cape, cin = dask.compute(cape, cin)
            print(type(cape))
        if sourcein=='surface':
            np.testing.assert_almost_equal(cape, dssurf.capesp500.values, decimal_cape)
            np.testing.assert_almost_equal(cin, dssurf.cinsp500.values, decimal_cin)
        elif sourcein=='mixed-layer':
            np.testing.assert_almost_equal(cape, dssurf.capeml300p500.values, decimal_cape)
            np.testing.assert_almost_equal(cin, dssurf.cinml300p500.values, decimal_cin)
        

@pytest.mark.parametrize('use_dask', [False, True])
@pytest.mark.parametrize('output_var_in, n_returns',
                         [('srh',1)])
#                          [('all', 4),('srh',1)])
@pytest.mark.parametrize('vertical_levin', ['sigma', 'pressure'])
def test_calc_srh(dataset_soundings, output_var_in, n_returns, use_dask, vertical_levin):
    """Test SRH code"""
    ds = dataset_soundings
    if use_dask:
        ds = ds.chunk()
    if vertical_levin=='sigma':
        returns = calc_srh(ds.pressure.data[:,1:],
                              ds.temperature.data[:,1:],
                              ds.dewpoint.data[:,1:],
                              ds.u_wind_ms.data[:,1:],
                              ds.v_wind_ms.data[:,1:],                         
                              ds.pressure.data[:,0],
                              ds.temperature.data[:,0],
                              ds.dewpoint.data[:,0],
                              ds.u_wind_ms.data[:,0],
                              ds.v_wind_ms.data[:,0],
                              depth = 3000,
                              vertical_lev=vertical_levin, 
                              output_var=output_var_in)
    elif  vertical_levin=='pressure':
        returns = calc_srh(ds.pressure.data[:,1:],
                              ds.temperature.data[:,1:],
                              ds.dewpoint.data[:,1:],
                              ds.u_wind_ms.data[:,1:],
                              ds.v_wind_ms.data[:,1:],                         
                              ds.pressure.data[:,0],
                              ds.temperature.data[:,0],
                              ds.dewpoint.data[:,0],
                              ds.u_wind_ms.data[:,0],
                              ds.v_wind_ms.data[:,0],
                              ds.pressure.data[:,0]*0, #pres_lev_pos
                              depth = 3000,
                              vertical_lev=vertical_levin, 
                              output_var=output_var_in)
        

#     print(returns)
#     print(type(returns))
    if output_var_in=='all':
        srh = returns[0]
        rm = returns[1]
        lm = returns[2] 
        mean_6km = returns[3] 
        if use_dask:
            assert isinstance(srh, dsa.Array)
            assert isinstance(rm, dsa.Array)
            assert isinstance(lm, dsa.Array)
            assert isinstance(mean_6km, dsa.Array)
            srh, rm, lm, mean_6km = dask.compute(srh, rm, lm, mean_6km)
            np.testing.assert_almost_equal(srh, ds.SRH03_model_lev.values, 5)
    else:
        srh = returns[0]
        print(type(srh))
        if use_dask:
            srh=returns
            assert isinstance(srh, dsa.Array)
            srh = dask.compute(srh)
            print(type(srh))
            srh = srh[0]
            np.testing.assert_almost_equal(srh, ds.SRH03_model_lev.values, 5)

@pytest.mark.parametrize('use_dask', [False, True])
@pytest.mark.parametrize('output_var_in, n_returns',
                         [('srh',1)])
#                          [('all', 4),('srh',1)])
@pytest.mark.parametrize('vertical_levin', ['pressure'])
def test_calc_srh_pressure1d(dataset_ERA5pressurelevel, output_var_in, n_returns, use_dask, vertical_levin):
    """Test SRH code"""
    ds3d, dssurf = dataset_ERA5pressurelevel
    if use_dask:
        ds3d = ds3d.chunk()
        dssurf = dssurf.chunk()
    returns = calc_srh(ds3d.level.data,
                          ds3d.t.data,
                          ds3d.td.data,
                          ds3d.u.data,
                          ds3d.v.data,
                          dssurf.p.data,
                          dssurf.t.data,
                          dssurf.td.data,
                          dssurf.u.data,
                          dssurf.v.data,
                          dssurf.start3d.data, #pres_lev_pos
                          depth = 3000,
                          vertical_lev=vertical_levin,
                          output_var=output_var_in)
    if output_var_in=='all':
        srh = returns[0]
        rm = returns[1]
        lm = returns[2] 
        mean_6km = returns[3] 
        if use_dask:
            assert isinstance(srh, dsa.Array)
            assert isinstance(rm, dsa.Array)
            assert isinstance(lm, dsa.Array)
            assert isinstance(mean_6km, dsa.Array)
            srh, rm, lm, mean_6km = dask.compute(srh, rm, lm, mean_6km)
            np.testing.assert_almost_equal(srh, dssurf.srh.values, 5)
    else:
        srh = returns[0]
        if use_dask:
            srh=returns
            assert isinstance(srh, dsa.Array)
            srh = dask.compute(srh)
            print(type(srh))
            srh = srh[0]
            np.testing.assert_almost_equal(srh, dssurf.srh.values, 5)