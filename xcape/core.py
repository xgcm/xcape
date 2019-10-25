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

def _cape_dummy(p, t, td, ps, ts, tds,
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
    pres_lev_pos :  array-like, optional
        location in python values (0: nlev-1) of where p <= ps.
        When vertical_lev='model', pres_lev_pos is set to a flag = 1
        When vertical_lev='pressure', pres_lev_pos.shape = ps.shape
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
        
    if len(args)>7:
        raise ValueError("Too many arguments.")     
        
    allowed_vertical_levs = ['sigma', 'pressure']
    if kwargs['vertical_lev'] not in allowed_vertical_levs:
        raise ValueError(f"`vertrical_lev` must be one of: {allowed_vertical_levs}")
        
    if (len(args)==6)&(kwargs['vertical_lev'] == 'pressure'):
        raise ValueError(f"`vertical grid` is set to `pressure`,\n",
                "but location in fortran values (1: nlev) of where p <= ps is missing")
        
    if (len(args)==7):
        if (args[6].shape!=args[3].shape)&(kwargs['vertical_lev'] == 'pressure'):
            raise ValueError(f'`vertical grid` is set to `pressure`,\n',
                    'but location in python values (0: nlev-1) of where p <= ps is not the correct shape\n',
                    '(pres_lev_pos.shape should be equal to ps.shape)')

        
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

    p, t, td, ps, ts, tds, *pres_lev_pos_in = args
    
    original_shape = t.shape #shape of 3D variable, i.e. t (p could be 1d)
    original_surface_shape = ts.shape #shape of surface variable, i.e. ps

    # after this, all arrays are 2d shape (nlevs, npoints)

    if len(p.shape) == 1:
        t_2d, td_2d = _reshape_inputs(t, td)
        p_2d = _reshape_inputs(p)[0]
        flag_1d = 1
    elif p.shape[-1] == t.shape[-1]:
        p_2d, t_2d, td_2d = _reshape_inputs(p, t, td)
        flag_1d = 0
        
   
    p_s1d, t_s1d, td_s1d, *pres_lev_pos = _reshape_surface_inputs(ps, ts, tds, *pres_lev_pos_in) 
    
    if len(args)==7:
        pres_lev_pos = pres_lev_pos[0]+1 # to fortran convention
    elif len(args)==6:
        pres_lev_pos = 1
        
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
#     elif method == 'numba':
#         cape_2d, cin_2d = _cape_numba(p_2d, t_2d, td_2d, **kwargs)
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
    pres_lev_pos :  array-like, optional
        location in python values (0: nlev-1) of where p <= ps.
        When vertical_lev='model', pres_lev_pos is set to a flag = 1
        When vertical_lev='pressure', pres_lev_pos.shape = ps.shape        
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
        signature = "(i),(i),(i),(i),(i),(),(),(),(),()->()"
        output_dtypes = ('f4',)
    elif (kwargs['vertical_lev']=='pressure'):
        signature = "(i),(i),(i),(i),(i),(),(),(),(),(),()->()"
        output_dtypes = ('f4',)
    if kwargs['output_var']=='all':
        signature += ",(),(),()"
        output_dtypes = output_dtypes + ('f4','f4','f4')
        
    return da.apply_gufunc(_calc_srh_numpy, signature,
                               *args,
                               output_dtypes=output_dtypes,
                               axis=-1,
                               vectorize=False,
                               **kwargs)
    

# the numpy version of the algorithm
def _calc_srh_numpy( *args,
                    depth = 3000, vertical_lev='sigma',output_var='all'):    
    
    p, t, td, u, v,  ps, ts, tds, us, vs, *pres_lev_pos_in = args
    original_shape = p.shape
    original_surface_shape = ps.shape

    # after this, all arrays are 2d shape (nlevs, npoints)
    p_2d, t_2d, td_2d, u_2d, v_2d = _reshape_inputs(p, t, td, u, v)
    # after this, all surface arrays are 1d shape (npoints)    
    p_s1d, t_s1d, td_s1d, u_s1d, v_s1d, *pres_lev_pos = _reshape_surface_inputs(ps, ts, 
                                                                        tds, us, vs, 
                                                                        *pres_lev_pos_in) # to fortran convention
    if len(args)==11:
        pres_lev_pos = pres_lev_pos[0]+1 # to fortran convention
    elif len(args)==10:
        pres_lev_pos = 1
    
    _vertical_lev_options_ ={'sigma':1, 'pressure':2}
    _output_var_options = {'srh':1, 'all':2}        

    kwargs_stdh = dict(type_grid=_vertical_lev_options_[vertical_lev])
    
    aglh_2d, aglh_s1d = _stdheight(p_2d, t_2d, td_2d,
                                   p_s1d, t_s1d, td_s1d,
                                   pres_lev_pos, aglh0 = 2.,
                                   **kwargs_stdh)
    
    kwargs = dict(type_grid=_vertical_lev_options_[vertical_lev],
                  output = _output_var_options[output_var])
    
    
    if _output_var_options[output_var] == 1:
        srh_2d = _srh(u_2d, v_2d, aglh_2d, 
                      u_s1d, v_s1d, aglh_s1d,
                      pres_lev_pos, depth, 
                      **kwargs)
        #_reshape_outputs returns a list
        srh = _reshape_outputs(srh_2d, shape=original_shape)[0]
        print('h', type(srh))
        return srh
    else:
        srh_2d, rm_2d, lm_2d, mean_6km_2d = _srh(u_2d, v_2d, aglh_2d, 
                      u_s1d, v_s1d, aglh_s1d,
                      pres_lev_pos, depth, 
                      **kwargs)

        srh = _reshape_outputs(srh_2d, shape=original_shape)[0]
        rm, lm, mean_6km = _reshape_outputs_uv_components(rm_2d, lm_2d, mean_6km_2d, shape=original_shape)
        print('j', type(srh))
        return srh, rm, lm, mean_6km
