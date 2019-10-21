
def srh(u_2d, v_2d, aglh_2d, u_s, v_s, aglh_s, pres_lev_pos, depth, type_grid, output):
    
                
    # nlev has to be the first dimension
    # nlev here is the number of levels in 3d variables (without surface level)
    nlev, ngrid = u_2d.shape
    # type_grid  type of vertical grid: 1 for model levels, 2 for pressure levels:
    if type_grid == 1:
        from xcape import Bunkers_model_lev
        from xcape import SREH_model_lev
    
        rm_sup,lm_sup,mean_6km = Bunkers_model_lev.bunkers_loop_ml(u_2d, v_2d, aglh_2d, 
                                                 u_s, v_s, aglh_s)
        srh = SREH_model_lev.loop_sreh_ml(u_2d, v_2d, aglh_2d, 
                                       u_s, v_s, aglh_s,
                                  rm_sup[0,:],
                                  rm_sup[1,:],depth)

    elif type_grid == 2:
        from xcape import Bunkers_pressure_lev
        from xcape import SREH_pressure_lev
    
        rm_sup,lm_sup,mean_6km = Bunkers_pressure_lev.bunkers_loop_pl(u_2d, v_2d, aglh_2d, 
                                                 u_s, v_s, aglh_s,
                                                          pres_lev_pos)
        
        srh = SREH_pressure_lev.loop_sreh_pl(u_2d, v_2d, aglh_2d, 
                                       u_s, v_s, aglh_s,
                                  rm_sup[0,:],
                                  rm_sup[1,:],depth,
                                         pres_lev_pos)
        
                
    if output ==1:
        return srh
    else:
        return srh, rm_sup, lm_sup, mean_6km
    #pass

