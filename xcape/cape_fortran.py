
def cape(p_2d, t_2d, td_2d, source, ml_depth, adiabat, pinc):
    from xcape import ALLCAPELOOP
    # nlev has to be the first dimension

    nlev, ngrid = p_2d.shape
    
    CAPE, CIN, MUlev, zMUlev = ALLCAPELOOP.loopcape(p_2d, t_2d, td_2d,
                                                    pinc, source, ml_depth, adiabat,
                                                    nlev, ngrid)
    print(MUlev)
    return CAPE, CIN, MUlev, zMUlev
    #pass
