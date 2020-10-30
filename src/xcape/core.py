#Copyright (c) 2020 xcape Developers.
"""
Numpy API for xcape, for calculation of Convective Available Potential Energy (CAPE)
and Storm Relative Helicity (SRH).
"""

from functools import reduce
import numpy as np
import dask.array as da
import xarray as xr

# put this statement within try so that i don't import need to import
# and compile the modules when I build the docs
try:
    from .cape import cape as _cape
except ImportError:
    # allows us to import on readthedocs
    cape = None

from .srh import srh as _srh
from .stdheight import stdheight as _stdheight

def _any_dask_array(*args):
    return any([isinstance(a, da.Array) for a in args])

def _any_xarray_da(*args):
    return any([isinstance(a, xr.DataArray) for a in args])


def calc_cape(*args, **kwargs):
    """
    Calculate Convective Available Potential Energy (CAPE) and Convective Inhibition (CIN).


    Calculate the CAPE and CIN over a a 3D gridded field, by iterating over a point profile
    integration of the area between an environmental vertical profile and a specified parcel profile.
    The integration is performed by a trapezoidal approach iteratively over a specified pressure
    increment. Parcel properties are user selected between a surface-based, a specified
    mixed-layer depth and most unstable. Adiabat for parcel trajectory is specifiable as one of
    pseudo-liquid, reversible-liquid, pseudo-ice and reversible-ice depending on the end application.
    Vertical level options should be specified based on the input model data, whether defined on
    pressure levels or model levels.


    Parameters
    ----------
    p : array-like
        Atmospheric pressure at each vertical level in hPa.
        When vertical_lev='model', p.shape = t.shape = (nlev, x, y, ...)
        When vertical_lev='pressure', p.shape = t.shape[0] = (nlev)
    t : array-like
        Atmospheric temperature in Celsius. Vertical shape should be identical to pressure.
    td : array-like
        Atmospheric dew point temperature in Celsius. Vertical shape should be identical to pressure.
    ps : array-like
        Surface Pressure in hPa.
    ts : array-like
        Surface Temperature in Celsius.
    tds : array-like
        Surface dew point temperature in Celsius.
    source : str, optional, default is 'surface'
        Select parcel based on desired assumptions under parcel theory
            - 'surface' = Surface-based parcels are subject to substantial errors depending on surface heating and source data, and can be influenced by moisture depth.
            - 'mixed-layer' = Mixed-layer parcels are generally a good assumption for profiles at peak heating when the boundary layer is deeply mixed to approximately the boundary layer depth.
            - 'most-unstable' = Most-unstable is defined by the layer below 500hPa with the highest equivalent potential temperature.
    ml_depth : float, optional (default is 500)
        Depth (m) of mixed layer. Only applies when the source='mixed-layer'
    adiabat : str, optional, default is 'pseudo-liquid'
        options include: 'pseudo-liquid', 'reversible-liquid','pseudo-ice', 'reversible-ice'
    pinc : float, optional, default is 500
        Pressure increment for integration (Pa) - Recommended usage (between 1000 and 100) is
        based on desired speed, with accuracy of the calculation increasing with smaller
        integration increments.
    vertical_lev : str, optional, default is 'sigma'
        Option to select vertical grid, between model coordinates, 'sigma', and pressure levels, 'pressure'.

    Returns
    -------
    cape : array-like
        Convective available potential energy (J/Kg)
    cin : array-like
        Convective inhibition (J/Kg)
    MUlev : array-like
        Most Unstable level location index (only returned for source: {'most-unstable'})
    zMUlev : array-like
        height of MUlev (m) (only returned for source: {'most-unstable'}

    Notes
    -----
    CAPE is calculated on a user specified set of parcel options based on the integration  [1]_:


    .. math:: CAPE = g \\int_{LFC}^{EL} (T_{v,parcel} - T_{v,env}/(T_{v,env}) \\text{dz}

    .. math:: CIN = g \\int_{SFC}^{LFC} (T_{v,parcel} - T_{v,env})/(T_{v,env}) \\text{dz}

    * :math:`CAPE` = Convective available potential energy
    * :math:`CIN` = Convective inhibition
    * :math:`LFC` = Level of free convection
    * :math:`EL` = Equilibrium level
    * :math:`g` = Gravitational acceleration
    * :math:`T_{v,parcel}` = Virtual potential temperature of the parcel
    * :math:`T_{v,env}` = Virtual potential temperature of the environment
    * :math:`z` = height above ground

    Examples
    -------
    Example of usage:

    >>> cape, cin = core.calc_cape(p, t, td, ps, ts, tds, source ='mixed-layer',
                    mldepth=500., adiabat='pseudo-liquid', pinc = 500.,
                    vertical_lev='sigma')

    References
    ----------
    .. [1] Emanuel, K. A. (1994). Atmospheric convection. Oxford University Press on Demand.

    """

    if len(args)<6:
        raise ValueError("Too few arguments.")

    if len(args)>6:
        raise ValueError("Too many arguments.")

    allowed_vertical_levs = ['sigma', 'pressure']
    if kwargs['vertical_lev'] not in allowed_vertical_levs:
        raise ValueError(f"`vertical_lev` must be one of: {allowed_vertical_levs}")

    if _any_xarray_da(*args):
        return xr.apply_ufunc(_calc_cape_gufunc,
                              *args)
    elif _any_dask_array(*args):
        return _calc_cape_gufunc(*args, **kwargs)
    else:
        return _calc_cape_numpy(*args, **kwargs)


def _calc_cape_gufunc(*args, **kwargs):
    ''' Wrapped function for cape calculation for dask arrays to leverage parallelized calculation
        over the grid.
    '''

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


def _prod(v):
    '''
    Lambda function for products between two arrays.
    '''
    return reduce(lambda x, y: x*y, v)

def _reshape_inputs(*args):
    '''
    Function to reshape multiple arrays of input data from, i.e., native 4D (X,Y,T,Z) format
    to a 2D vertical array (XYT,Z), or any (X1,X2,X3,X4...,Z) to (X1*X2*X3*X4...,Z) to allow
    for parallelization of calculation for desired parameter.
    Called by functions _calc_*_numpy.
    '''
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
    '''
    Function to reshape arrays of calculated output data (X1*X2*X3*X4...Z) to the original
    input shape minus the Z dimension (X1,X2,X3,X4,...) to allow
    for parallelization of calculation for desired parameter.
    Called by _calc_*_numpy.
    '''
    if len(shape)==1:
        target_shape = (1,)
    else:
        target_shape = shape[:-1]
    return [np.reshape(a.transpose(), target_shape) for a in args]

def _reshape_outputs_uv_components(*args, shape=None):
    '''
    Function to reshape arrays of calculated output data with 2 components (X1*X2*X3*X4...Z,2)
    to the original input shape minus the Z dimension (X1,X2,X3,X4,...,2) to allow
    for parallelization of calculation for desired parameter. This is most commonly
    applicable to wind or storm motion components.
    Called by _calc_srh_numpy.
    '''
    if len(shape)==1:
        target_shape = (2,)
    else:
        # 1 is in place of the level dimension
        # shape[1:] is the remaining shape
        target_shape = (2,) + shape[:-1]
    return [np.reshape(a, target_shape) for a in args]

def _reshape_surface_inputs(*args):
    '''
    Function to reshape multiple arrays of input data from, i.e., native 3D (X,Y,T) format
    to a 2D vertical array (XYT), or any (X1,X2,X3,X4...) to (X1*X2*X3*X4...) to allow
    for parallelization of calculation for desired parameter.
    Called by functions _calc_*_numpy.
    '''
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

# the numpy version of the algorithm
def _calc_cape_numpy(*args,
                     source='surface', ml_depth=500.,
                     adiabat='pseudo-liquid',pinc=500.,
                     vertical_lev='sigma'):
    '''
    Wrapper function for cape calculation to setup optional kwargs and ensure data
    is provided in a format suitable output to call CAPE and CIN calculation.
    '''


    p, t, td, ps, ts, tds = args

    original_shape = t.shape #shape of 3D variable, i.e. t (p could be 1d)
    original_surface_shape = ts.shape #shape of surface variable, i.e. ps

    # after this, all arrays are 2d shape (nlevs, npoints)
    p_s1d, t_s1d, td_s1d = _reshape_surface_inputs(ps, ts, tds)
    if len(p.shape) == 1:
        t_2d, td_2d = _reshape_inputs(t, td)
        p_2d = _reshape_inputs(p)[0]
        flag_1d = 1
        # calculate pres_lev_pos to use only levels where p_2d <= p_s1d
        # in pressure level data, data are regridded and back-propagated
        # for the whole vertical grid, even when surface pressure is lower than
        # the first level of the grid.
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

    cape_2d, cin_2d, mulev, zmulev = _cape(p_2d, t_2d, td_2d,
                                                   p_s1d, t_s1d, td_s1d,
                                                   flag_1d,
                                                   pres_lev_pos,
                                                   **kwargs)

    # reshaping outputs
    if _source_options_[source]==2:
        cape, cin, mulev, zmulev = _reshape_outputs(cape_2d, cin_2d, mulev, zmulev, shape=original_shape)
        return cape, cin, mulev, zmulev
    else:
        cape, cin = _reshape_outputs(cape_2d, cin_2d, shape=original_shape)
        return cape, cin


def calc_srh(*args, **kwargs):
    """
    Calculate storm relative helicity for a vectorized set of profiles over the first axis of the arrays.

    Calculate SRH over a predefined depth for a N-dimensional gridded field.
    Performs a point profile integration of the area between the hodograph and
    the vectors between the estimated storm motion vector at the bottom and top
    layer. The calculation assumes that the storm motion is described by the
    un-weighted mean 0-6 km wind vector subsampled at 500m intervals, with a
    left or right moving storm deviating by 7.5 m/s perpendicular to this
    mean wind, following the empirically estimated approach defined by
    Bunkers et al. (2000). Output can be specified to include the assumed
    motion, or to return only the helicity as calculated for the user-defined
    layer. Vertical level option should be specified based on the input model
    data, whether defined on pressure or model levels.

    Parameters
    ----------
    p : array-like
        Atmospheric pressure in hPa.
        - When vertical_lev='model', p.shape = t.shape = (nlev, x, y, ...)
        - When vertical_lev='pressure', p.shape = t.shape[0] = (nlev)
    t : array-like
        Atmospheric temperature in Celsius. Vertical shape should be identical to pressure.
    td : array-like
        Atmospheric dew point temperature in Celsius. Vertical shape should be identical to pressure.
    u  : array-like
        Zonal wind component in meters per second. Vertical shape should be identical to pressure.
    v  : array-like
        Meridional wind component in meters per second. Vertical shape should be identical to pressure.
    ps : array-like
        Surface Pressure in hPa.
    ts : array-like
        Surface Temperature in Celsius.
    tds : array-like
        Surface dew point temperature in Celsius.
    us : array-like
        10m zonal wind field in meters per second.
    us : array-like
        10m meridional wind field in meters per second.
    depth : float, optional
        Depth in meters (m) of the layer to calculate SRH.
    vertical_lev : str, optional, default is 'sigma'
        Option to select vertical grid, between model coordinates, 'sigma', and pressure levels, 'pressure'.
    output_var : str, optional, default is 'srh'
        - 'srh' = for only srh
        - 'all' = for srh, Bunkers' right-moving and left-moving storm component, mean not pressure averaged 6km wind
        Option to return either calculated SRH only, or also the Bunker's storm
        motion for the right-moving and left-moving storms, along with the not
        pressure-weighted mean 0-6 km wind.

    Returns
    -------
    srh_rm : array-like
          Storm relative helicity for the right-moving storm.
    srh_lm : array-like
          Storm relative helicity for the left-moving storm.
    bunkers_rm_u: array-like
          Estimated storm motion for the right-moving storm, u component.
    bunkers_rm_v: array-like
          Estimated storm motion for the right-moving storm, v component.
    bunkers_lm_u: array-like
          Estimated storm motion for the left-moving storm, u component.
    bunkers_lm_v: array-like
          Estimated storm motion for the left-moving storm, v component.
    mean_6km_u: array-like
          Mean wind between 0 and 6km without pressure weighting, u component.
    mean_6km_v: array-like
          Mean wind between 0 and 6km without pressure weighting, v component.

    Notes
    -----

    Calculates SRH for a user specified depth based on the numeric integration:

    .. math:: SRH = \int_{0}^{h} (\\bar v - c) \cdot \\bar \omega_{h} \,dz

    * :math:`SRH` Storm Relative Helicity
    * :math:`\\bar v` Environmental Wind Vector
    * :math:`\\bar c` Storm Motion Vector (left and right components) following Bunkers et al. (2000) [1]_
    * :math:`\\bar\omega_{h}` Horizontal Vorticity Vector
    * :math:`z` height above ground

    For further details see Markowski and Richardson [2010], pages 230-231 [2]_


    Examples
    -------
    Example of usage:

    >>> srh_rm, srh_lm = core.calc_srh(p, t, td, u, v, ps, ts, tds, us, vs,
                                       depth=1000,
                                       output_var='srh',
                                       vertical_lev='sigma')

    >>> srh_rm, srh_lm, rm_u, rm_v, lm_u, lm_v, mean_6km_u, mean_6km_v = core.calc_srh(
                                       p, t, td, u, v, ps, ts, tds, us, vs,
                                       depth=1000,
                                       output_var='srh',
                                       vertical_lev='sigma')

    References
    ----------
    .. [1] Bunkers, M. J., B. A. Klimowski, J. W. Zeitler, R. L. Thompson, and M. L.
       Weisman, 2000: Predicting supercell motion using a new hodograph technique.
       Wea. Forecasting, 15, 61-79.
    .. [2] Markowski, P., & Richardson, Y. (2011). Mesoscale meteorology in midlatitudes (Vol. 2). John Wiley & Sons.

    """
    if _any_dask_array(*args):
        return _calc_srh_gufunc(*args, **kwargs)
    else:
        return _calc_srh_numpy(*args, **kwargs)

def _calc_srh_gufunc(*args, **kwargs):
    ''' Wrapped function for srh calculation for dask arrays to leverage parallelized calculation
        over the grid.
    '''

    if (kwargs['vertical_lev']=='sigma'):
        signature = "(i),(i),(i),(i),(i),(),(),(),(),()->(),()"
        output_dtypes = ('f4','f4')
    elif (kwargs['vertical_lev']=='pressure'):
        signature = "(i),(i),(i),(i),(i),(),(),(),(),(),()->(),()"
        output_dtypes = ('f4','f4')
    if kwargs['output_var']=='all':
        signature += ",(),(),(),(),(),()" #",(2),(2),(2)"
        output_dtypes = output_dtypes +  ('f4','f4','f4','f4','f4','f4') #('f4','f4','f4')

    return da.apply_gufunc(_calc_srh_numpy, signature,
                               *args,
                               output_dtypes=output_dtypes,
                               axis=-1,
                               vectorize=False,
                               **kwargs)


# the numpy version of the algorithm
def _calc_srh_numpy(*args,
                    depth = 3000, vertical_lev='sigma',output_var='srh'):
    '''
    Wrapper function for SRH calculation to setup optional kwargs and
    ensure data is provided in a format suitable output to call the
    Fortran implementation of SRH calculation. Also calculates standard
    height using the hypsometric equation, which necessitates the passing
    of pressure, temperature and dewpoint temperature. Called by
    _calc_srh_gufunc.
    '''

    p, t, td, u, v,  ps, ts, tds, us, vs = args
    original_shape = t.shape #shape of 3D variable, i.e. p
    original_surface_shape = ts.shape #shape of surface variable, i.e. ps

    # after this, all arrays are 2d shape (nlevs, npoints)
    p_s1d, t_s1d, td_s1d, u_s1d, v_s1d = _reshape_surface_inputs(ps, ts,tds, us, vs)
    if len(p.shape) == 1:
        t_2d, td_2d, u_2d, v_2d = _reshape_inputs(t, td, u, v)
        p_2d = _reshape_inputs(p)[0]
        flag_1d = 1
        # calculate pres_lev_pos to use only levels where p_2d <= p_s1d
        # in pressure level data, data are regridded and back-propagated
        # for the whole vertical grid, even when surface pressure is lower than
        # the first level of the grid.
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

    # calculates standard height using the hypsometric equation
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
        return  srh_rm, srh_lm, rm[0],rm[1], lm[0],lm[1], mean_6km[0], mean_6km[1] #srh_rm, srh_lm, rm, lm, mean_6km
