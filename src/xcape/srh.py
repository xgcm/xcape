from .fortran import Bunkers_model_lev, Bunkers_pressure_lev
from .fortran import SREH_model_lev, SREH_pressure_lev

def srh(u_2d, v_2d, aglh_2d, u_s, v_s, aglh_s, pres_lev_pos, depth, type_grid, output):
    """
     Description:
     ------------
     Function called by calc_srh (core.py) to perform Bunkers storm motion
     and derived SRH calculations for specified depth using the wrapped
     compiled fortran code for either pressure level or model level data
     on 2-dimensions (collapsed grid and vertical) which are specified
     in the call to this function. For full description, see function
     calc_srh in core.py, or in Bunkers_model_lev.f90.

     Parameters:
     ------------
     u_2d,v_2d :    'array-like'
                    Gridded component winds in m/s (nk, ngrid).
     aglh_2d:       'array-like'
                    Gridded height above ground level in m (nk, ngrid).
     u_s,v_s:       'array-like'
                    Component winds at the lowest model level or surface in m/s.
                    (ngrid).
     aglh_s:        'array-like'
                    Height of the surface level above ground level in m (ngrid)
     press_lev_pos: 'Integer'
                    For type_grid='pressure' indicates the levels where p_2d <= p_s1d
     depth:         'Integer'
                    Depth for SRH integration to be performed in m.
     type_grid:     'Integer'
                    Whether the input data is 'sigma':1 or 'pressure':2 level.
     output:        'Integer'
                    If set to 1, only return SRH, otherwise return SRH and
                    motions.
    """

    # nlev has to be the first dimension
    # nlev here is the number of levels in 3d variables (without surface level)
    nlev, ngrid = u_2d.shape
    # type_grid  type of vertical grid: 1 for model levels, 2 for pressure levels:
    if type_grid == 1:
        rm_sup,lm_sup,mean_6km = Bunkers_model_lev.bunkers_loop_ml(u_2d, v_2d, aglh_2d,
                                                                   u_s, v_s, aglh_s)
        srh_rm, srh_lm = SREH_model_lev.loop_sreh_ml(u_2d, v_2d, aglh_2d,
                                                    u_s, v_s, aglh_s,
                                                    rm_sup[0,:],
                                                    rm_sup[1,:],
                                                    lm_sup[0,:],
                                                    lm_sup[1,:], depth)

    elif type_grid == 2:
        rm_sup,lm_sup,mean_6km = Bunkers_pressure_lev.bunkers_loop_pl(u_2d, v_2d, aglh_2d,
                                                                      u_s, v_s, aglh_s,
                                                                      pres_lev_pos)
        srh_rm, srh_lm = SREH_pressure_lev.loop_sreh_pl(u_2d, v_2d, aglh_2d,
                                                        u_s, v_s, aglh_s,
                                                        rm_sup[0,:],
                                                        rm_sup[1,:],
                                                        lm_sup[0,:],
                                                        lm_sup[1,:],depth,
                                                        pres_lev_pos)

    if output ==1:
        return srh_rm, srh_lm
    else:
        return srh_rm, srh_lm, rm_sup, lm_sup, mean_6km
    #pass
