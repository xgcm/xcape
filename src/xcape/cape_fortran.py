from . import CAPE_CODE_model_lev, CAPE_CODE_pressure_lev

def cape(p_2d, t_2d, td_2d, p_s, t_s, td_s, flag_1d, pres_lev_pos, source, ml_depth, adiabat, pinc, type_grid):


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
        elif flag_1d == 1:
            CAPE, CIN, MUlev, zMUlev = CAPE_CODE_pressure_lev.loopcape_pl1d(t_2d, td_2d, p_2d,
                                                              p_s, t_s, td_s,
                                                              pinc, source, ml_depth, adiabat,
                                                              pres_lev_pos,
                                                              nlev, ngrid)

    return CAPE, CIN, MUlev, zMUlev
    #pass
