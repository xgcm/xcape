
def cape(p_2d, t_2d, td_2d, p_s, t_s, td_s, source, ml_depth, adiabat, pinc, type_grid):
    from xcape import CAPE_CODE
    # nlev has to be the first dimension

    nlev, ngrid = p_2d.shape
    
    CAPE, CIN, MUlev, zMUlev = CAPE_CODE.loopcape(p_2d, t_2d, td_2d,
                                                  p_s, t_s, td_s,
                                                    pinc, source, ml_depth, adiabat, type_grid,
                                                    nlev, ngrid)
#     print(MUlev)
    return CAPE, CIN, MUlev, zMUlev
    #pass
