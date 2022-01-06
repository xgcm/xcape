import numpy as np
import xarray as xr

from itertools import combinations
import dask.array as dsa
import dask

from xcape.core import calc_cape
from xcape.core import calc_srh
from .fixtures import empty_dask_array, dataset_soundings, dataset_ERA5pressurelevel

import pytest


def _random_data(*shape):
    # p, t, td
    return [np.random.rand(*shape),
            np.random.rand(*shape),
            np.random.rand(*shape)]

@pytest.fixture(scope='module')
def p_t_td_1d(nlevs=20):
    p, t, td = _random_data(nlevs)
    return p, t, td

@pytest.fixture(scope='module')
def p_t_td_3d(nlevs=20, nx=10, ny=5):
    p, t, td = _random_data(nx, ny, nlevs)
    return p, t, td

@pytest.fixture(scope='module')
def p_t_td_3d_xr(nlevs=20, nx=10, ny=5):
    p, t, td = _random_data(nx, ny, nlevs)
    dims = ['lat', 'lon', 'lev']
    ds = xr.Dataset({'p': (dims, p),
                     't': (dims, t),
                     'td': (dims, td)})
    return ds

@pytest.fixture(scope='module')
def p_t_td_surface(nx=10, ny=5):
    ps, ts, tds = _random_data(nx, ny)
    return ps, ts, tds


@pytest.fixture(scope='module')
def p_t_td_surface_xr(nx=10, ny=5):
    ps, ts, tds = _random_data(nx, ny)
    dims = ['lat', 'lon']
    ds = xr.Dataset({'ps': (dims, ps),
                     'ts': (dims, ts),
                     'tds': (dims, tds)})
    return ds

# surface mode returns cape, cin
# most-unstable mode returns cape, cin, mulev, zmulev

# TEST CHANGE IN SHAPE OF VARIABLES IN
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
        args = (np.ones(t.shape[-1]), t, td, ps, ts, tds)
    if use_dask:
        args = [dsa.from_array(a) for a in args]
    result = calc_cape(*args, source=sourcein,vertical_lev=vertical_levin)
    assert len(result) == n_returns
    for data in result:
        assert data.shape == p.shape[:-1]
        if use_dask:
            assert isinstance(data, dsa.Array)
            data.compute()


@pytest.mark.parametrize('use_dask', [False, True])
@pytest.mark.parametrize('sourcein, n_returns',
                         [('surface', 2), ('most-unstable', 4)])
@pytest.mark.parametrize('transpose', [True, False])
@pytest.mark.parametrize('vertical_levin', ['sigma', 'pressure'])
def test_calc_cape_shape_3d_xr(p_t_td_3d_xr, p_t_td_surface_xr, sourcein, n_returns,
                               use_dask, vertical_levin, transpose):
    ds = p_t_td_3d_xr
    ds_surf = p_t_td_surface_xr
    if transpose:
        # make sure we can put lev on any axis
        ds = ds.transpose('lev', 'lat', 'lon')
    if use_dask:
        ds = ds.chunk({'lat': 1})
        ds_surf = ds_surf.chunk({'lat': 1})
    if vertical_levin=='sigma':
        args = (ds.p, ds.t, ds.td, ds_surf.ps, ds_surf.ts, ds_surf.tds)
    elif vertical_levin =='pressure':
        args = (ds.lev, ds.t, ds.td, ds_surf.ps, ds_surf.ts, ds_surf.tds)

    result = calc_cape(*args, dim='lev', source=sourcein, vertical_lev=vertical_levin)

    if n_returns == 2:
        assert list(result) == ['CAPE', 'CIN']
    else:
        assert list(result) == ['CAPE', 'CIN', 'mulev', 'z_mulev']

    expected_dims = list(ds.t.dims).copy()
    expected_dims.remove('lev')
    for vname in result:
        da = result[vname]
        assert da.dims == tuple(expected_dims)
        if use_dask:
            assert isinstance(da.data, dsa.Array)
            da.compute()


# tolerance for tests
decimal_cape = 0
decimal_cin = 0
decimal_mulv = 0
decimal_zmulv = 0
# TEST FOR CAPE CALCULATION FOR MODEL LEVEL USING THE DATASET_SOUNDING DATA
# THAT MIMIC THE MODEL LEVELS
@pytest.mark.parametrize('use_dask', [False, True])
@pytest.mark.parametrize('sourcein, pinc_used',
                         [('surface', 100),('mixed-layer',1000),('most-unstable', 100)])
@pytest.mark.parametrize('vertical_levin', ['sigma'])
def test_calc_cape_sigma(dataset_soundings, sourcein, pinc_used, use_dask,vertical_levin):
    """Test Surface Cape based on previously calculated using George Bryans code"""
    ds = dataset_soundings
    if use_dask:
        ds = ds.chunk()
    returns = calc_cape(ds.pressure.data[:, 1:],
                          ds.temperature.data[:, 1:],
                          ds.dewpoint.data[:, 1:],
                          ds.pressure.data[:, 0],
                          ds.temperature.data[:, 0],
                          ds.dewpoint.data[:, 0],
                          source=sourcein, ml_depth=500., adiabat='pseudo-liquid',
                          pinc=pinc_used,
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

# TEST FOR CAPE CALCULATION FOR PRESSURE LEVEL USING THE DATASET_ERA5pressurelevel DATA
# THAT ARE PRESSURE LEVEL DATA
@pytest.mark.parametrize('use_dask', [False, True])
@pytest.mark.parametrize('sourcein', ['surface', 'mixed-layer', 'most-unstable'])
@pytest.mark.parametrize('vertical_levin', ['pressure'])
@pytest.mark.parametrize('pinc_used', [500])
@pytest.mark.parametrize('ml_depthin', [300])
def test_calc_cape_pressure(dataset_ERA5pressurelevel, sourcein, pinc_used, ml_depthin, use_dask,vertical_levin):
    """Test Surface Cape based on previously calculated using George Bryans code"""
    ds3d, dssurf = dataset_ERA5pressurelevel
    if use_dask:
        ds3d = ds3d.chunk()
        dssurf = dssurf.chunk()
    # here i use ds3d.level.data - which is 1D as it should be, I don't use the ds3d.p variable in fixture
    returns = calc_cape(ds3d.level.data,
                          ds3d.t.data,
                          ds3d.td.data,
                          dssurf.p.data,
                          dssurf.t.data,
                          dssurf.td.data,
                          source=sourcein, ml_depth=ml_depthin, adiabat='pseudo-liquid',
                          pinc=pinc_used,
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
#                          [('srh',1)])
                         [('all', 8),('srh',1)])
@pytest.mark.parametrize('vertical_levin', ['sigma'])
def test_calc_srh_sigma(dataset_soundings, output_var_in, n_returns, use_dask, vertical_levin):
    """Test SRH code sigma/model level"""
    ds = dataset_soundings
    if use_dask:
        ds = ds.chunk()
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


#     print(returns)
#     print(type(returns))
    if output_var_in=='all':
        srh_rm = returns[0]
        srh_lm = returns[1]
        rm0 = returns[2]
        rm1 = returns[3]
        lm0 = returns[4]
        lm1 = returns[5]
        mean_6km0 = returns[6]
        mean_6km1 = returns[7]
        if use_dask:
            assert isinstance(srh_rm, dsa.Array)
            assert isinstance(srh_lm, dsa.Array)
            assert isinstance(rm0, dsa.Array)
            assert isinstance(rm1, dsa.Array)
            assert isinstance(lm0, dsa.Array)
            assert isinstance(lm1, dsa.Array)
            assert isinstance(mean_6km0, dsa.Array)
            assert isinstance(mean_6km1, dsa.Array)
            srh_rm, srh_lm, rm0, rm1, lm0, lm1, mean_6km0,mean_6km1 = dask.compute(srh_rm, srh_lm, rm0, rm1, lm0, lm1, mean_6km0,mean_6km1)
        np.testing.assert_almost_equal(srh_rm, ds.SRH03_model_lev_rm.values, 5)
        np.testing.assert_almost_equal(srh_lm, ds.SRH03_model_lev_lm.values, 5)
    else:
        srh_rm = returns[0]
        srh_lm = returns[1]
        if use_dask:
            assert isinstance(srh_rm, dsa.Array)
            assert isinstance(srh_lm, dsa.Array)
            srh_rm, srh_lm = dask.compute(srh_rm, srh_lm)
        np.testing.assert_almost_equal(srh_rm, ds.SRH03_model_lev_rm.values, 5)
        np.testing.assert_almost_equal(srh_lm, ds.SRH03_model_lev_lm.values, 5)

@pytest.mark.parametrize('use_dask', [False, True])
@pytest.mark.parametrize('output_var_in, n_returns',
#                          [('srh',1)])
                         [('all', 8),('srh',1)])
@pytest.mark.parametrize('vertical_levin', ['pressure'])
def test_calc_srh_pressure(dataset_ERA5pressurelevel, output_var_in, n_returns, use_dask, vertical_levin):
    """Test SRH code pressure level"""
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
                          depth = 3000,
                          vertical_lev=vertical_levin,
                          output_var=output_var_in)
    if output_var_in=='all':
        srh_rm = returns[0]
        srh_lm = returns[1]
        rm0 = returns[2]
        rm1 = returns[3]
        lm0 = returns[4]
        lm1 = returns[5]
        mean_6km0 = returns[6]
        mean_6km1 = returns[7]
        if use_dask:
            assert isinstance(srh_rm, dsa.Array)
            assert isinstance(srh_lm, dsa.Array)
            assert isinstance(rm0, dsa.Array)
            assert isinstance(rm1, dsa.Array)
            assert isinstance(lm0, dsa.Array)
            assert isinstance(lm1, dsa.Array)
            assert isinstance(mean_6km0, dsa.Array)
            assert isinstance(mean_6km1, dsa.Array)
            srh_rm, srh_lm, rm0, rm1, lm0, lm1, mean_6km0,mean_6km1 = dask.compute(srh_rm, srh_lm, rm0, rm1, lm0, lm1, mean_6km0,mean_6km1)

        np.testing.assert_almost_equal(srh_rm, dssurf.srh_rm.values, 5)
        np.testing.assert_almost_equal(srh_lm, dssurf.srh_lm.values, 5)
    else:
        srh_rm = returns[0]
        srh_lm = returns[1]
        if use_dask:
            assert isinstance(srh_rm, dsa.Array)
            assert isinstance(srh_lm, dsa.Array)
            srh_rm, srh_lm = dask.compute(srh_rm, srh_lm)

        np.testing.assert_almost_equal(srh_rm, dssurf.srh_rm.values, 5)
        np.testing.assert_almost_equal(srh_lm, dssurf.srh_lm.values, 5)
