"""
Numpy API for xcape.
"""

from functools import reduce
import numpy as np
from .duck_array_ops import (reshape, ravel_multi_index, concatenate,
                             broadcast_arrays)

from .cape_fortran import cape as _cape_fortran
from .cape_numba import cape as _cape_numba

def _prod(v):
    return reduce(lambda x, y: x*y, v)

def _reshape_inputs(*args):
    a0 = args[0]
    shape = a0.shape
    for a in args:
        if a.shape != shape:
            raise ValueError('Input arrays must have the same shape.')

    args_al2d = [np.atleast_2d(a) for a in args]
    original_shape = args_al2d[0].shape
    # for now CAPE works with shape [lev , prod(lat*lon*time)]
    # lev is the first dimension
    #new_shape = (_prod(original_shape[:-1]),) + (original_shape[-1],)
    new_shape = (original_shape[0],) + (_prod(original_shape[1:]),) 
    print(new_shape)
    args_2d = [np.reshape(a, new_shape) for a in args_al2d]
    return args_2d

def _reshape_outputs(*args, shape=None):
    # for now CAPE works with shape [lev , prod(lat*lon*time)]
    # lev is the first dimension
    
    if len(shape)==1:
        target_shape = (1,)
    else:
#         target_shape = shape[:-1] + (1,)
        target_shape = (1,) + shape[1:]
    return [np.reshape(a, target_shape) for a in args]

def _cape_dummy(p, t, td, **kwargs):
    # cape is a reduction along the first axis.
    assert p.ndim == 2
    shape = p.shape
    cape = np.ones((shape[0], 1))
    cin = np.ones((shape[0], 1))
    return cape, cin

def calc_cape(p, t, td, source='surface', ml_depth=500., adiabat='pseudo-liquid',
         pinc=1000., method='fortran', vertical_lev='model'):
    """
    Calculate cape for a set of profiles over the final axis of the arrays.

    Parameters
    ----------
    p : array-like
        Pressure in mb.
    t : array-like
        Temperature in Celsius
    td : array-like
        Dew point temperature in Celsius
    source : {'surface', 'most-unstable', 'mixed-layer'}
    ml_depth : float, optional
        Depth (m) of mixed layer. Only applies when source='mixed-layer'
    adiabat : {'pseudo-liquid', 'reversible-liquid','pseudo-ice', 'reversible-ice'}
    pinc : float, optional
        Pressure increment (Pa) - Recommended between 1000. (faster) and 100 (slower)
    method : {'fortran', 'numba'}
        Which numerical routine to use
    vertical_lev : {'model', 'pressure'}
        Which vertical grid is used

    Returns
    -------
    cape : array-like
        Convective available potential energy (J/Kg)
    cin : array-like
        Convective inhibition (J/Kg)
    """

    original_shape = p.shape
    print(original_shape)

    # after this, all arrays are 2d shape (npoints, nlevs)
    p_2d, t_2d, td_2d = _reshape_inputs(p, t, td)
    _source_options_ ={'surface':1, 'most-unstable':2, 'mixed-layer':3}
    _adiabat_options_ ={'pseudo-liquid':1, 'reversible-liquid':2,
                        'pseudo-ice':3, 'reversible-ice':4}
            
    kwargs = dict(source=_source_options_[source], 
                  ml_depth=ml_depth, 
                  adiabat=_adiabat_options_[adiabat], 
                  pinc=pinc)

    if method == 'fortran':
        cape_2d, cin_2d, mulev, zmulev = _cape_fortran(p_2d, t_2d, td_2d, **kwargs)
    elif method == 'numba':
        cape_2d, cin_2d = _cape_numba(p_2d, t_2d, td_2d, **kwargs)
    elif method == 'dummy':
        cape_2d, cin_2d = _cape_dummy(p_2d, t_2d, td_2d, **kwargs)
    else:
        raise ValueError('invalid method')
    
    cape, cin = _reshape_outputs(cape_2d, cin_2d, shape=original_shape)
    return cape, cin, mulev, zmulev
