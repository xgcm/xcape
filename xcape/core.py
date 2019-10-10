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
    # calc_cape needs input in shape (nlevs, npoints)
    new_shape = (original_shape[0],) + (_prod(original_shape[1:]),)
    args_2d = [np.reshape(a, new_shape) for a in args_al2d]
    return args_2d

def _reshape_outputs(*args, shape=None):
    # calc_cape needs input in shape (nlevs, npoints)

    if len(shape)==1:
        target_shape = (1,)
    else:
        target_shape = (1,) + shape[1:]
    return [np.reshape(a, target_shape) for a in args]

def _reshape_surface_inputs(*args):
    a0 = args[0]
    shape = a0.shape
    for a in args:
        if a.shape != shape:
            raise ValueError('Input arrays must have the same shape.')

    args_al1d = [np.atleast_1d(a) for a in args]
    original_shape = args_al1d[0].shape

    new_shape = (_prod(original_shape),)

    args_1d = [np.reshape(a, new_shape) for a in args_al1d]
    return args_1d

def _reshape_surface_outputs(*args, shape=None):
    # calc_cape needs input in shape (nlevs, npoints)

    if len(shape)==1:
        target_shape = (1,)
    else:
        target_shape = (1,) + shape[1:]
    return [np.reshape(a, target_shape) for a in args]


def _cape_dummy(p, t, td, ps, ts, tds, pres_lev_pos,
                source, ml_depth, adiabat, pinc, type_grid):
    # cape is a reduction along the second axis.
    # this tests that reshaping works.
    # calc_cape needs input in shape (nlevs, npoints)
    assert p.ndim == 2
    shape = p.shape
    cape = np.ones((1, shape[1]))
    cin = np.ones((1, shape[1]))
    mulev = np.ones((1, shape[1]))
    zmulev = np.ones((1, shape[1]))
    return cape, cin, mulev, zmulev

def calc_cape(p, t, td, ps, ts, tds, source='surface', ml_depth=500., adiabat='pseudo-liquid',
         pinc=1000., method='fortran', vertical_lev='sigma', pres_lev_pos=1):
    """
    Calculate cape for a set of profiles over the first axis of the arrays.

    Parameters
    ----------
    p : array-like
        Pressure in mb.
        When vertical_lev='model', p.shape = t.shape = (nlev, x, y, ...)
        When vertical_lev='pressure', p.shape = t.shape[0] = (nlev)
    t : array-like
        Temperature in Celsius
    td : array-like
        Dew point temperature in Celsius
    ps : array-like
        Surface Pressure in mb.
    ts : array-like
        Surface Temperature in Celsius
    tds : array-like
        Surface Dew point temperature in Celsius
    source : {'surface', 'most-unstable', 'mixed-layer'}
    ml_depth : float, optional
        Depth (m) of mixed layer. Only applies when source='mixed-layer'
    adiabat : {'pseudo-liquid', 'reversible-liquid','pseudo-ice', 'reversible-ice'}
    pinc : float, optional
        Pressure increment (Pa) - Recommended between 1000. (faster) and 100 (slower)
    method : {'fortran', 'numba'}
        Which numerical routine to use
    vertical_lev : {'sigma', 'pressure'}
        Which vertical grid is used
    pres_lev_pos :  array-like,
        location in fortran values (1: nlev) of where p <= ps.
        When vertical_lev='model', pres_lev_pos = 1
        When vertical_lev='pressure', pres_lev_pos.shape = ps.shape
    Returns
    -------
    cape : array-like
        Convective available potential energy (J/Kg)
    cin : array-like
        Convective inhibition (J/Kg)
    MUlev : array-like
        Most Unstable level location index (-)
    zMUlev : array-like
        height of MUlev (m)
     """

    original_shape = p.shape
    original_surface_shape = ps.shape

    # after this, all arrays are 2d shape (nlevs, npoints)
    p_2d, t_2d, td_2d = _reshape_inputs(p, t, td)
    # after this, all surface arrays are 1d shape (npoints)
    p_s1d, t_s1d, td_s1d = _reshape_surface_inputs(ps, ts, tds)

    _source_options_ ={'surface':1, 'most-unstable':2, 'mixed-layer':3}
    _adiabat_options_ ={'pseudo-liquid':1, 'reversible-liquid':2,
                        'pseudo-ice':3, 'reversible-ice':4}
    _vertical_lev_options_ ={'sigma':1, 'pressure':2}

    kwargs = dict(source=_source_options_[source],
                  ml_depth=ml_depth,
                  adiabat=_adiabat_options_[adiabat],
                  pinc=pinc,
                  type_grid=_vertical_lev_options_[vertical_lev])

    if method == 'fortran':
        cape_2d, cin_2d, mulev, zmulev = _cape_fortran(p_2d, t_2d, td_2d,
                                                       p_s1d, t_s1d, td_s1d,
                                                       pres_lev_pos,
                                                       **kwargs)
    elif method == 'numba':
        cape_2d, cin_2d = _cape_numba(p_2d, t_2d, td_2d, **kwargs)
    elif method == 'dummy':
        cape_2d, cin_2d, mulev, zmulev = _cape_dummy(p_2d, t_2d, td_2d,
                                                     p_s1d, t_s1d, td_s1d,
                                                     pres_lev_pos,
                                                     **kwargs)
    else:
        raise ValueError('invalid method')


    if _source_options_[source]==2:
        cape, cin, mulev, zmulev = _reshape_outputs(cape_2d, cin_2d, mulev, zmulev, shape=original_shape)
        return cape, cin, mulev, zmulev
    else:
        cape, cin = _reshape_outputs(cape_2d, cin_2d, shape=original_shape)
        return cape, cin
    

def calc_srh(p, t, td, u, v,  ps, ts, tds, us, vs, depth = 3000, 
             vertical_lev='sigma', pres_lev_pos=1):
    """
    Calculate cape for a set of profiles over the first axis of the arrays.

    Parameters
    ----------
    p : array-like
        Pressure in mb.
        When vertical_lev='model', p.shape = t.shape = (nlev, x, y, ...)
        When vertical_lev='pressure', p.shape = t.shape[0] = (nlev)
    t : array-like
        Temperature in Celsius
    td : array-like
        Dew point temperature in Celsius
        
    ps : array-like
        Surface Pressure in mb.
    ts : array-like
        Surface Temperature in Celsius
    tds : array-like
        Surface Dew point temperature in Celsius
        
    depth : float, optional
        Depth (m) of SRH layer.

    vertical_lev : {'sigma', 'pressure'}
        Which vertical grid is used
    pres_lev_pos :  array-like,
        location in fortran values (1: nlev) of where p <= ps. 
        When vertical_lev='model', pres_lev_pos = 1
        When vertical_lev='pressure', pres_lev_pos.shape = ps.shape
    Returns
    -------
    srh : array-like
    """

    original_shape = p.shape
    print(original_shape)
    original_surface_shape = ps.shape
    print(original_surface_shape)

    # after this, all arrays are 2d shape (nlevs, npoints)
    p_2d, t_2d, td_2d, u2d, v2d = _reshape_inputs(p, t, td, u, v)
    # after this, all surface arrays are 1d shape (npoints)    
    p_s1d, t_s1d, td_s1d, u_s1d, v_s1d = _reshape_surface_inputs(ps, ts, tds, us, vs)
    print(p_2d.shape, p_s1d.shape)
    
    _vertical_lev_options_ ={'sigma':1, 'pressure':2}
            
    kwargs = dict(depth=ml_depth, 
                  type_grid=_vertical_lev_options_[vertical_lev])
    print(kwargs)
    aglh_2d, aglh_s1d = _stdheight(p_2d, t_2d, td_2d,
                                   p_s1d, t_s1d, td_s1d,
                                   pres_lev_pos, aglh0 = 2.,
                                   **kwargs)
    
    srh_2d = _srh(u_2d, v_2d, aglh_2d, 
                  u_s1d, v_s1d, aglh_s1d,
                  pres_lev_pos,
                  **kwargs)
    
    
    srh = _reshape_outputs(srh_2d, shape=original_shape)
    return srh
    
    
    from conv_parameters import SREH_model_lev
from conv_parameters import Bunkers_model_lev
from conv_parameters import stdheight_2D_model_lev
start = timer()

H2D_split, H2D_surface  = stdheight_2D_model_lev.loop_stdheight(Pst[1:],Tst[1:],Tdst[1:],
                                                                Pst[0],Tst[0],Tdst[0],
                                                                Hinst*0+2.,
                                                                Pst.shape[0]-1,Pst.shape[1])
end = timer()
print('stdheight fortran 2D split',end-start)

rm_sup3,lm_sup3,meanwind_6km3 = Bunkers_model_lev.loop(U.reshape(P.shape[0],241*480)[1:],
                                    V.reshape(P.shape[0],241*480)[1:],
                                    H2D_split,
                                    U.reshape(P.shape[0],241*480)[0],
                                    V.reshape(P.shape[0],241*480)[0],
                                    H2D_surface)


srh3kmb4 = SREH_model_lev.loop_sreh(U.reshape(P.shape[0],241*480)[1:],
                                    V.reshape(P.shape[0],241*480)[1:],
                                    AGLH.reshape(P.shape[0],241*480)[1:],
                                    U.reshape(P.shape[0],241*480)[0],
                                    V.reshape(P.shape[0],241*480)[0],
                                    AGLH.reshape(P.shape[0],241*480)[0],
                          rm_sup2[0,:],
                          rm_sup2[1,:],3000.)
end = timer()
print('loop fortran flat4 ',end-start)
