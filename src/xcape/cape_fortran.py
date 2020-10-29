from .fortran import CAPE_CODE_model_lev, CAPE_CODE_pressure_lev

def cape(p_2d, t_2d, td_2d, p_s, t_s, td_s, flag_1d, pres_lev_pos, source, ml_depth, adiabat, pinc, type_grid):

    """
     Description:
     ------------
     Function called by calc_cape (core.py) to perform the integration
     of CAPE for the specified parcel using the wrapped compiled fortran
     code for either pressure level or model level data on 2-dimensions
     (collapsed grid and vertical) which are specified in the call to
     this function. For full description, see function calc_cape in core.py,
     or CAPE_CODE_model_lev.f90, CAPE_CODE_press_lev.f90.

     Parameters:
     ------------
     p_2d:          'array-like'
                    Gridded component pressure in hPa (nk, ngrid), or (nk).
     t_2d:          'array-like'
                    Gridded temperature in K (nk, ngrid).
     td_2d:         'array-like'
                    Gridded dewpoint temperature in K (nk, ngrid).
     p_s:           'array-like'
                    Gridded surface pressure in hPa (ngrid).
     t_s:           'array-like'
                    Gridded temperature at the surface (2m) in K (ngrid).
     td_s:          'array-like'
                    Gridded dewpoint temperature at the surface (2m) in K (ngrid).
     flag_1d:       'integer'
                    Flag set to 1 if pressure is a 1D vector, 0 if 2D or greater.
     press_lev_pos: 'integer'
                    Specifies whether pressure level data is on the edge or center of a vertical grid.
     source:        'string'
                    Optional, default is 'surface', select parcel based on desired assumptions under parcel theory.
     ml_depth:      'float'
                    Depth (m) of mixed layer. Only applies when the source='mixed-layer'.
     adiabat:       'string'
                    Optional, default is 'pseudo-liquid', select parcel adiabat.
     pinc : 	    'float'
                    Optional, default is 500. Pressure increment for integration (Pa).
     type_grid:     'Integer'
                    Whether the input data is 'sigma':1 or 'pressure':2 level.
    """

    # nlev has to be the first dimension
    # nlev here is the number of levels in 3d variables (without surface level)
    nlev, ngrid = t_2d.shape
    # type_grid  type of vertical grid: 1 for model levels, 2 for pressure levels:
    if type_grid == 1:
        CAPE, CIN, MUlev, zMUlev = CAPE_CODE_model_lev.loopcape_ml(p_2d, t_2d, td_2d,
                                                      p_s, t_s, td_s,
                                                      pinc, source, ml_depth, adiabat,
                                                      nlev, ngrid)
    elif type_grid == 2:
        # this option should not be triggered because of what is in core.py
#         if flag_1d == 0:
#             CAPE, CIN, MUlev, zMUlev = CAPE_CODE_pressure_lev.loopcape_pl(p_2d, t_2d, td_2d,
#                                                               p_s, t_s, td_s,
#                                                               pinc, source, ml_depth, adiabat,
#                                                               pres_lev_pos,
#                                                               nlev, ngrid)
        if flag_1d == 1:
            CAPE, CIN, MUlev, zMUlev = CAPE_CODE_pressure_lev.loopcape_pl1d(t_2d, td_2d, p_2d,
                                                              p_s, t_s, td_s,
                                                              pinc, source, ml_depth, adiabat,
                                                              pres_lev_pos,
                                                              nlev, ngrid)

    return CAPE, CIN, MUlev, zMUlev
    #pass
