"""
Numpy API for xcape.
"""

from functools import reduce
import numpy as np
import dask.array as da
from .duck_array_ops import (reshape, ravel_multi_index, concatenate,
                             broadcast_arrays)

from .cape_fortran import cape as _cape_fortran
from .cape_numba import cape as _cape_numba
from .srh import srh as _srh
from .stdheight import stdheight as _stdheight
from .Indices_calc import Indices_calc as _Indices_calc

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
    new_shape = (_prod(original_shape[:-1]),) + (original_shape[-1],)
    # transpose to move nlevs to first axis
    args_2d = [np.reshape(a, new_shape).transpose() for a in args_al2d]
    return args_2d

def _reshape_outputs(*args, shape=None):
    if len(shape)==1:
        target_shape = (1,)
    else:
        target_shape = shape[:-1]
    return [np.reshape(a.transpose(), target_shape) for a in args]

def _reshape_outputs_uv_components(*args, shape=None):
    if len(shape)==1:
        target_shape = (2,)
    else:
        # 1 is in place of the level dimension
        # shape[1:] is the remaining shape
        target_shape = (2,) + shape[:-1]
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

# def _reshape_outputs_stdheight(*args, shape=None):
#     # calc_cape needs input in shape (nlevs, npoints)
#     target_shape = shape
#     return [np.reshape(a, target_shape) for a in args]

def _cape_dummy(*args,**kwargs):
    p, t, td, ps, ts, tds = args
    # cape is a reduction along the second axis.
    # this tests that reshaping works.
    # calc_cape needs input in shape (nlevs, npoints)
    assert p.ndim == 2
    shape = t.shape
    cape = np.ones((1, shape[1]))
    cin = np.ones((1, shape[1]))
    mulev = np.ones((1, shape[1]))
    zmulev = np.ones((1, shape[1]))
    return cape, cin, mulev, zmulev


def _any_dask_array(*args):
    return any([isinstance(a, da.Array) for a in args])

def calc_cape(*args, **kwargs):
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
    if len(args)<6:
        raise ValueError("Too little arguments.")     
        
    if len(args)>6:
        raise ValueError("Too many arguments.")     
        
    allowed_vertical_levs = ['sigma', 'pressure']
    if kwargs['vertical_lev'] not in allowed_vertical_levs:
        raise ValueError(f"`vertrical_lev` must be one of: {allowed_vertical_levs}")
        
    if _any_dask_array(*args):
        return _calc_cape_gufunc(*args, **kwargs)
    else:
        return _calc_cape_numpy(*args, **kwargs)
    

def _calc_cape_gufunc(*args, **kwargs):
    
    if (kwargs['vertical_lev']=='sigma'):
        signature = "(i),(i),(i),(),(),()->(),()"
        output_dtypes = ('f4','f4')
    elif (kwargs['vertical_lev']=='pressure'):
        signature = "(i),(i),(i),(),(),(),()->(),()"
        output_dtypes = ('f4','f4')
        
    if kwargs['source']=='most-unstable':
        signature += ",(),()"
        output_dtypes = output_dtypes + ('i4','f4')

    return da.apply_gufunc(_calc_cape_numpy, signature,
                               *args, 
                               output_dtypes=output_dtypes,
                               axis=-1,
                               vectorize=False,
                               **kwargs)
    
# the numpy version of the algorithm
def _calc_cape_numpy(*args, 
                     source='surface', ml_depth=500.,
                     adiabat='pseudo-liquid',pinc=500., method='fortran',
                     vertical_lev='sigma'):

    p, t, td, ps, ts, tds = args
    
    original_shape = t.shape #shape of 3D variable, i.e. t (p could be 1d)
    original_surface_shape = ts.shape #shape of surface variable, i.e. ps

    # after this, all arrays are 2d shape (nlevs, npoints)
    p_s1d, t_s1d, td_s1d = _reshape_surface_inputs(ps, ts, tds) 
    if len(p.shape) == 1:
        t_2d, td_2d = _reshape_inputs(t, td)
        p_2d = _reshape_inputs(p)[0]
        flag_1d = 1
        # calculate pres_lev_pos
        temp_index = (p_s1d-p_2d)
        pres_lev_pos = np.ma.masked_less(temp_index,0).argmin(axis=0)
        # fortran convention
        pres_lev_pos = pres_lev_pos+1
        
    elif (p.shape == t.shape)&(vertical_lev=='sigma'):
        p_2d, t_2d, td_2d = _reshape_inputs(p, t, td)
        flag_1d = 0
        pres_lev_pos = 1

    elif (p.shape == t.shape)&(vertical_lev=='pressure'):
        raise ValueError("P should be 1d")     

    
        
    
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
                                                       flag_1d,
                                                       pres_lev_pos,
                                                       **kwargs)
    elif method == 'dummy':
        cape_2d, cin_2d, mulev, zmulev = _cape_dummy(p_2d, t_2d, td_2d,
                                                     p_s1d, t_s1d, td_s1d,
                                                     **kwargs)
    else:
        raise ValueError('invalid method')

    
    if _source_options_[source]==2:
        cape, cin, mulev, zmulev = _reshape_outputs(cape_2d, cin_2d, mulev, zmulev, shape=original_shape)
        return cape, cin, mulev, zmulev
    else:
        cape, cin = _reshape_outputs(cape_2d, cin_2d, shape=original_shape)
        return cape, cin
    

def calc_srh(*args, **kwargs):
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
    output_var : {'srh', 'all'}
        'srh' = for only srh
        'all' = for srh, Bunkers' right-moving and left-moving storm component, 
                mean not pressure averaged 6km wind
    Returns
    -------
    srh : array-like
    """
    if _any_dask_array(*args):
        return _calc_srh_gufunc(*args, **kwargs)
    else:
        return _calc_srh_numpy(*args, **kwargs)

def _calc_srh_gufunc(*args, **kwargs):

    
    if (kwargs['vertical_lev']=='sigma'):
        signature = "(i),(i),(i),(i),(i),(),(),(),(),()->(),()"
        output_dtypes = ('f4','f4')
    elif (kwargs['vertical_lev']=='pressure'):
        signature = "(i),(i),(i),(i),(i),(),(),(),(),(),()->(),()"
        output_dtypes = ('f4','f4')
    if kwargs['output_var']=='all':
        signature += ",(2),(2),(2)"
        output_dtypes = output_dtypes + ('f4','f4','f4')
        
    return da.apply_gufunc(_calc_srh_numpy, signature,
                               *args,
                               output_dtypes=output_dtypes,
                               axis=-1,
                               vectorize=False,
                               **kwargs)
    

# the numpy version of the algorithm
def _calc_srh_numpy(*args,
                    depth = 3000, vertical_lev='sigma',output_var='srh'):    
    
    p, t, td, u, v,  ps, ts, tds, us, vs = args
    original_shape = t.shape #shape of 3D variable, i.e. p
    original_surface_shape = ts.shape #shape of surface variable, i.e. ps
        
    # after this, all arrays are 2d shape (nlevs, npoints)
    p_s1d, t_s1d, td_s1d, u_s1d, v_s1d = _reshape_surface_inputs(ps, ts,tds, us, vs) 
    if len(p.shape) == 1:
        t_2d, td_2d, u_2d, v_2d = _reshape_inputs(t, td, u, v)
        p_2d = _reshape_inputs(p)[0]
        flag_1d = 1
        # calculate pres_lev_pos
        temp_index = (p_s1d-p_2d)
        pres_lev_pos = np.ma.masked_less(temp_index,0).argmin(axis=0)
        # fortran convention
        pres_lev_pos = pres_lev_pos+1

    elif (p.shape == t.shape)&(vertical_lev=='sigma'):
        p_2d, t_2d, td_2d, u_2d, v_2d = _reshape_inputs(p, t, td, u, v)
        flag_1d = 0
        pres_lev_pos = 1
    elif (p.shape == t.shape)&(vertical_lev=='pressure'):
        raise ValueError("P should be 1d")     
    
    _vertical_lev_options_ ={'sigma':1, 'pressure':2}
    _output_var_options = {'srh':1, 'all':2}        

    kwargs_stdh = dict(type_grid=_vertical_lev_options_[vertical_lev])

    aglh_2d, aglh_s1d = _stdheight(p_2d, t_2d, td_2d,
                                   p_s1d, t_s1d, td_s1d,
                                   flag_1d,
                                   pres_lev_pos, aglh0 = 2.,
                                   **kwargs_stdh)
    
    kwargs = dict(type_grid=_vertical_lev_options_[vertical_lev],
                  output = _output_var_options[output_var])
    
    
    if _output_var_options[output_var] == 1:
        srh_2d_rm, srh_2d_lm = _srh(u_2d, v_2d, aglh_2d, 
                      u_s1d, v_s1d, aglh_s1d,
                      pres_lev_pos, depth, 
                      **kwargs)
        #_reshape_outputs returns a list
        srh_rm, srh_lm = _reshape_outputs(srh_2d_rm, srh_2d_lm, shape=original_shape)
        return srh_rm, srh_lm
    else:
        srh_2d_rm, srh_2d_lm, rm_2d, lm_2d, mean_6km_2d = _srh(u_2d, v_2d, aglh_2d, 
                                                          u_s1d, v_s1d, aglh_s1d,
                                                          pres_lev_pos, depth, 
                                                          **kwargs)

        srh_rm, srh_lm = _reshape_outputs(srh_2d_rm, srh_2d_lm, shape=original_shape)
        rm, lm, mean_6km = _reshape_outputs_uv_components(rm_2d, lm_2d, mean_6km_2d, shape=original_shape)
        return srh_rm, srh_lm, rm, lm, mean_6km
    
def calc_indices(*args, **kwargs):
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
    u
    v
        
    ps : array-like
        Surface Pressure in mb.
    ts : array-like
        Surface Temperature in Celsius
    tds : array-like
        Surface Dew point temperature in Celsius
    us,
    vs
    
    vertical_lev : {'sigma', 'pressure'}
        Which vertical grid is used
    Returns
    -------
    LAPSE24, LAPSE3,  LAPSE700_500, THGZ, S06, SBLCL, T500, FZL : array-like
    """
    if _any_dask_array(*args):
        return _calc_srh_gufunc(*args, **kwargs)
    else:
        return _calc_srh_numpy(*args, **kwargs)

def _calc_indices_gufunc(*args, **kwargs):

    
    signatureA = "(i),(i),(i),(i),(i),(),(),(),(),()"
    signatureB = "->(),()"
    output_dtypes = ('f4','f4','f4','f4','f4','f4','f4','f4')
    if (kwargs['vertical_lev']=='pressure'):
        #ADD 1 INPUT FOR PRES_LEV_POS
        signatureA = signatureA + ",()"
    signature = signatureA + signatureB
    return da.apply_gufunc(_calc_indices_numpy, signature,
                               *args,
                               output_dtypes=output_dtypes,
                               axis=-1,
                               vectorize=False,
                               **kwargs)
    

# the numpy version of the algorithm
def _calc_indices_numpy(*args,
                        vertical_lev='sigma'):    
    
    #21            
    p, t, td, u, v,  ps, ts, tds, us, vs = args
    
    original_shape = t.shape #shape of 3D variable, i.e. p
    original_surface_shape = ts.shape #shape of surface variable, i.e. ps
        
    # after this, all arrays are 2d shape (nlevs, npoints)
    p_s1d, t_s1d, td_s1d, u_s1d, v_s1d = _reshape_surface_inputs(ps, ts, tds, us, vs) 
    
    if len(p.shape) == 1:
        t_2d, td_2d, u_2d, v_2d = _reshape_inputs(t, td, u, v)
        p_2d = _reshape_inputs(p)[0]
        flag_1d = 1
        # calculate pres_lev_pos
        temp_index = (p_s1d-p_2d)
        pres_lev_pos = np.ma.masked_less(temp_index,0).argmin(axis=0)
        # fortran convention
        pres_lev_pos = pres_lev_pos+1

    elif (p.shape == t.shape)&(vertical_lev=='sigma'):
        p_2d, t_2d, td_2d, u_2d, v_2d = _reshape_inputs(p, t, td, u, v)
        flag_1d = 0
        pres_lev_pos = 1
    elif (p.shape == t.shape)&(vertical_lev=='pressure'):
        raise ValueError("P should be 1d")     
    
    
    kwargs = dict(type_grid=_vertical_lev_options_[vertical_lev])
    
    
    LAPSE24_1d, LAPSE3_1d,  LAPSE700_500_1d, THGZ_1d, S06_1d, SBLCL_1d, T500_1d, FZL_1d = _Indices_calc(
                                                                                 p_2d, t_2d, td_2d, u_2d, v_2d, 
                                                                                 p_s1d, t_s1d, td_s1d, u_s1d, v_s1d, 
                                                                                 pres_lev_pos, 
                                                                                  **kwargs)
        #_reshape_outputs returns a list
    LAPSE24, LAPSE3,  LAPSE700_500, THGZ, S06, SBLCL, T500, FZL = _reshape_outputs(LAPSE24_1d, LAPSE3_1d,  
                                                                 LAPSE700_500_1d, THGZ_1d, 
                                                                 S06_1d, SBLCL_1d, T500_1d, FZL_1d, 
                                                                 shape=original_shape)
    
    return LAPSE24, LAPSE3,  LAPSE700_500, THGZ, S06, SBLCL, T500, FZL

