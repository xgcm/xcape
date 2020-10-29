import numpy as np
from .fortran import stdheight_2D_model_lev, stdheight_2D_pressure_lev


def stdheight(p_2d, t_2d, td_2d, p_s, t_s, td_s, flag_1d, pres_lev_pos, aglh0, type_grid):

    # nlev has to be the first dimension
    # nlev here is the number of levels in 3d variables (without surface level)
    nlev, ngrid = t_2d.shape

    # if above ground level height is 1 value
    # (i.e. 2m or 10m) generate ngrid-long array with aglh0
    if np.isscalar(aglh0):
        aglh_in = np.ones(ngrid)*aglh0
    else:
        # add check on shape
        aglh_in = aglh0

    # type_grid  type of vertical grid: 1 for model levels, 2 for pressure levels:
    if type_grid == 1:
        H2D, H_s  = stdheight_2D_model_lev.loop_stdheight_ml(p_2d, t_2d, td_2d,
                                                      p_s, t_s, td_s,
                                                      aglh_in,
                                                      nlev, ngrid)
    elif type_grid == 2:
        # this option should not be triggered because of what is in core.py
#         if flag_1d ==0:
#             H2D, H_s  = stdheight_2D_pressure_lev.loop_stdheight_pl(p_2d, t_2d, td_2d,
#                                                       p_s, t_s, td_s,
#                                                       aglh_in,pres_lev_pos,
#                                                       nlev, ngrid)
        if flag_1d ==1:
            H2D, H_s  = stdheight_2D_pressure_lev.loop_stdheight_pl1d(t_2d, td_2d, p_2d,
                                                      p_s, t_s, td_s,
                                                      aglh_in,pres_lev_pos,
                                                      nlev, ngrid)

    return H2D, H_s
