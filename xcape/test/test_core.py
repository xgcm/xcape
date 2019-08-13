import numpy as np

from itertools import combinations
import dask.array as dsa

from ..core import calc_cape
from .fixtures import empty_dask_array

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

# just see that the function can be called without an error and
# returns the correct shaped object
def test_calc_cape_shape_1d(p_t_td_1d):
    p, t, td = p_t_td_1d
    cape, cin = calc_cape(p, t, td, method='dummy')
    assert cape.shape == (1,)
    assert cin.shape == (1,)

def test_calc_cape_shape_3d(p_t_td_3d):
    p, t, td = p_t_td_3d
    cape, cin = calc_cape(p, t, td, method='dummy')
    assert cape.shape == (p.shape[0], p.shape[1], 1)
    assert cin.shape == (p.shape[0], p.shape[1], 1)
