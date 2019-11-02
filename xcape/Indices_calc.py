
def Indices_calc(p_2d, t_2d, td_2d, u_2d, v_2d,  
                 p_s, t_s, td_s,  u_s, v_s, 
                 pres_lev_pos, aglh0, type_grid):
    
#     def Indices_calc(p_2d, t_2d, td_2d, u_2d, v_2d, q_2d, 
#                  p_s, t_s, td_s,  u_s, v_s, q_s,
#                  MUlevs , camu_s, cas_s, cin_s, caml_s, srh_rm1_s, srh_lm1_s, srh_rm3_s, srh_lm3_s,
#                  pres_lev_pos, aglh0, type_grid):
    
    import numpy as np
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
        from xcape import stdheight_2D_model_lev
        H2D, H_s  = stdheight_2D_model_lev.loop_stdheight_ml(p_2d, t_2d, td_2d,
                                                      p_s, t_s, td_s,
                                                      aglh_in,
                                                      nlev, ngrid)
        from xcape import Interp_model_lev
        T2km = Interp_model_lev.interp_loop_ml(t_2d, H2D, t_s, H_s, 2000)
        T3km = Interp_model_lev.interp_loop_ml(t_2d, H2D, t_s, H_s, 3000)
        T4km = Interp_model_lev.interp_loop_ml(t_2d, H2D, t_s, H_s, 4000)

        T700 = Interp_model_lev.interp_loop_ml(t_2d, p_2d, t_s, p_s, 700)
        T500 = Interp_model_lev.interp_loop_ml(t_2d, p_2d, t_s, p_s, 500)
        z700 = Interp_model_lev.interp_loop_ml(H2D, p_2d, H_s, p_s, 700)
        z500 = Interp_model_lev.interp_loop_ml(H2D, p_2d, H_s, p_s, 500)


        FZL = Interp_model_lev.interp_loop_ml(H2D, t_2d, H_s, t_s, -0.001)
        HT30 = Interp_model_lev.interp_loop_ml(H2D, t_2d, H_s, t_s, -30)
        HT10 = Interp_model_lev.interp_loop_ml(H2D, t_2d, H_s, t_s, -10)

        U6km = Interp_model_lev.interp_loop_ml(u_2d, H2D, u_s, H_s, 6000)
        V6km = Interp_model_lev.interp_loop_ml(v_2d, H2D, v_s, H_s, 6000)
        U1km = Interp_model_lev.interp_loop_ml(u_2d, H2D, u_s, H_s, 1000)
        V1km = Interp_model_lev.interp_loop_ml(v_2d, H2D, v_s, H_s, 1000)

        
        
    elif type_grid == 2:
        from xcape import stdheight_2D_pressure_lev
#         if flag_1d ==0:
#             H2D, H_s  = stdheight_2D_pressure_lev.loop_stdheight_pl(p_2d, t_2d, td_2d,
#                                                       p_s, t_s, td_s,
#                                                       aglh_in,pres_lev_pos,
#                                                       nlev, ngrid)
#         if flag_1d ==1:
        H2D, H_s  = stdheight_2D_pressure_lev.loop_stdheight_pl1d(t_2d, td_2d, p_2d, 
                                                      p_s, t_s, td_s,
                                                      aglh_in,pres_lev_pos,
                                                      nlev, ngrid)        
        from xcape import Interp_pressure_lev
        T2km = Interp_pressure_lev.interp_loop_pl(t_2d, H2D, t_s, H_s, 2000,pres_lev_pos)
        T3km = Interp_pressure_lev.interp_loop_pl(t_2d, H2D, t_s, H_s, 3000,pres_lev_pos)
        T4km = Interp_pressure_lev.interp_loop_pl(t_2d, H2D, t_s, H_s, 4000,pres_lev_pos)

        T700 = Interp_pressure_lev.interp_loop_pl_y1d(t_2d, p_2d, t_s, p_s, 700,pres_lev_pos)
        T500 = Interp_pressure_lev.interp_loop_pl_y1d(t_2d, p_2d, t_s, p_s, 500,pres_lev_pos)
        z700 = Interp_pressure_lev.interp_loop_pl_y1d(H2D, p_2d, H_s, p_s, 700,pres_lev_pos)
        z500 = Interp_pressure_lev.interp_loop_pl_y1d(H2D, p_2d, H_s, p_s, 500,pres_lev_pos)


        FZL = Interp_pressure_lev.interp_loop_pl(H2D, t_2d, H_s, t_s, -0.00001,pres_lev_pos)
        HT30 = Interp_pressure_lev.interp_loop_pl(H2D, t_2d, H_s, t_s, -30,pres_lev_pos)
        HT10 = Interp_pressure_lev.interp_loop_pl(H2D, t_2d, H_s, t_s, -10,pres_lev_pos)

        U6km = Interp_pressure_lev.interp_loop_pl(u_2d, H2D, u_s, H_s, 6000,pres_lev_pos)
        V6km = Interp_pressure_lev.interp_loop_pl(v_2d, H2D, v_s, H_s, 6000,pres_lev_pos)
        U1km = Interp_pressure_lev.interp_loop_pl(u_2d, H2D, u_s, H_s, 1000,pres_lev_pos)
        V1km = Interp_pressure_lev.interp_loop_pl(v_2d, H2D, v_s, H_s, 1000,pres_lev_pos)

    # Calculate LAPSE RATES  - (Tsuperior-Tinferior)/(zsuperior-zinferior)
    LAPSE24 = - (T4km-T2km)/2
    LAPSE3 = - (T3km-t_s)/3
    LAPSE700_500 = - (T500-T700)/((z500-z700)/1000)
    THGZ = HT30-HT10
    # CALCULATE SBLCL
    SBLCL = 125.* ( (t_s-t_2d[0])/2. -td_s)
    # S06
    S06 = ( (U6km-u_s)**2 + (V6km-v_s)**2 )**0.5
    # S01
    S01 = ( (U1km-u_s)**2 + (V1km-v_s)**2 )**0.5
    
    return LAPSE24, LAPSE3,  LAPSE700_500, THGZ, S06, S01, SBLCL, T500, FZL

#     ################
#     ##### SHIP #####    
#     q_mulev = np.where(MUlevs==1, q_s, q_2d[MUlevs-1,np.arange(q_2d.shape[1])])
    
#     #SHIP
#     SHIP = (camu_s* q_mulev* LAPSE700_500 * (-T500) * S06)/(44_000_000)
#     SHIP = np.where(camu_s<1300, SHIP*(camu_s/1300), SHIP)
#     SHIP = np.where(LAPSE700_500<5.8, SHIP*(LAPSE700_500/5.8), SHIP)
#     SHIP = np.where(FZL<2400, SHIP*(FZL/2400), SHIP)
    
#     ################
#     ##### STP
#     Aterm = cas_s/1500.
#     Bterm = (2000.-SBLCL)/1000.
#     Bterm[SBLCL<1000] = 1
#     Bterm[SBLCL>2000] = 0
#     Cterm_rm = srh_rm1_s/150.
#     Cterm_lm = srh_lm1_s/150.
#     Dterm = S06/20.
#     Dterm[S06>30] =  1.5
#     Dterm[S06 < 12.5] = 0.
#     Eterm = np.fabs(cin_s)
#     Eterm[np.fabs(cin_s)>125]=0.
#     Eterm[np.fabs(cin_s)<=125]=1.
#     STP_rm = Aterm * Bterm * Cterm_rm * Dterm * Eterm
#     STP_lm = Aterm * Bterm * Cterm_lm * Dterm * Eterm
#     print(STP.shape)

#     ################
#     ##### SCP#
#     scp1 = cas_s/1000.
#     scp2_rm = srh_rm3_s/100.
#     scp2_lm = srh_lm3_s/100.
#     scp3 = S06/20.
#     SCP_rm = scp1*scp2_rm*scp3*Eterm
#     SCP_lm = scp1*scp2_lm*scp3*Eterm

#     ################
#     ##### EHI
#     EHI3_rm = ((caml_s)*(srh_rm3_s))/(1.6*10**5)
#     EHI1_rm = ((caml_s)*(srh_rm1_s))/(1.6*10**5)
#     EHI3_lm = ((caml_s)*(srh_lm3_s))/(1.6*10**5)
#     EHI1_lm = ((caml_s)*(srh_lm1_s))/(1.6*10**5)


#     return LAPSE24, LAPSE3,  LAPSE700_500, THGZ, S06, SHIP, STP_rm, STP_lm, SCP_rm, SCP_lm, EHI1_rm, EHI1_lm, EHI3_rm, EHI3_lm

