!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
    subroutine bunkers_loop_pl(U3d, V3d, AGLH3d, Us, Vs, AGLHs, start_3d, nk, n2, RM, LM, Mean6kmwind)
!-----------------------------------------------------------------------
!  bunkers_loop - loop along n2 dimension to calculate RM,LM,Mean6kmwind.
!
!  Author:  Chiara Lepore @chiaral - John Allen @xebadir
!
!  Disclaimer:  This code is made available WITHOUT WARRANTY.
!-----------------------------------------------------------------------
!  U3d,V3d = zonal and meridional wind
!  AGLH3d = above ground level height
!  nk = number of pressure levels
!  n2 = number of grid point for which calculate RM,LM,Mean6kmwind
!-----------------------------------------------------------------------
    implicit none

    !f2py threadsafe
    !f2py intent(out) :: RM3d, LM3d, Mean6kmwind3d

    integer, intent(in) :: nk,n2
    integer :: i, nk_all, nk_pl_in, nk_start
    real, dimension(nk,n2), intent(in) :: U3d, V3d, AGLH3d
    real, dimension(n2), intent(in) :: Us, Vs, AGLHs, start_3d
    real, dimension(2,n2), intent(out) :: RM, LM, Mean6kmwind
    ! at most, the number of levels will be 3d +  surface
    real, dimension(nk+1) :: U_all, V_all, AGLH_all

    do i = 1, n2
        nk_start = start_3d(i)
        ! compute number of used levels in 3d
        nk_pl_in = nk - nk_start + 1
        ! compute number of used levels in 3d + surface
        nk_all = nk_pl_in+1

        U_all(1) = Us(i)
        v_all(1) = Vs(i)
        AGLH_all(1) = AGLHs(i)
        U_all(2:nk_all) = U3d(nk_start:nk,i)
        V_all(2:nk_all) = V3d(nk_start:nk,i)
        AGLH_all(2:nk_all) = AGLH3d(nk_start:nk,i)

        call bunkers_calc_pl(U_all(1:nk_all), V_all(1:nk_all), AGLH_all(1:nk_all), nk_all, &
                          &RM(:,i), LM(:,i), Mean6kmwind(:,i))
    enddo
    return
    end subroutine bunkers_loop_pl


!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
    subroutine bunkers_calc_pl(U, V, AGLH, nk, RM, LM, Mean6kmwind)
!-----------------------------------------------------------------------
!  bunkers_calc - loop along n2 dimension to calculate RM,LM,Mean6kmwind.
!
!  Author:  Chiara Lepore @chiaral - John Allen @xebadir
!
!  Disclaimer:  This code is made available WITHOUT WARRANTY.
!-----------------------------------------------------------------------
!  U,V = zonal and meridional wind
!  AGLH = above ground level height
!  nk = number of pressure levels
!-----------------------------------------------------------------------
    implicit none

    !f2py threadsafe
    !f2py intent(out) :: RM, LM, Mean6kmwind

    integer, intent(in) :: nk
    integer :: i, k1, k2
    real, dimension(nk), intent(in) :: U, V, AGLH
    real, dimension(2), intent(out) :: RM, LM, Mean6kmwind
    real :: ushr, vshr, ushr_u, vshr_u, ushr_d, vshr_d
    real, dimension(12) :: inter_lvls
    real, dimension(13) :: u_surf6km, v_surf6km

    inter_lvls = (/500.,1000.,1500.,2000.,2500.,3000., 3500.,4000.,4500.,5000.,5500.,6000. /)

    !#Perform interpolation for layers 500-6000m
    do i=1,13
      IF (i .EQ. 1) THEN
        u_surf6km(i)=U(i)
        v_surf6km(i)=V(i)
      ELSE
        call DINTERP2DZ(U,AGLH,inter_lvls(i-1),u_surf6km(i),1,nk)
        call DINTERP2DZ(V,AGLH,inter_lvls(i-1),v_surf6km(i),1,nk)
      ENDIF

    enddo

    !#Take the mean wind 0-6km - non-pressure weighted.
    k1=0
    k2=0
    do i=1,13
      !IF (u_surf6km(i).EQ.u_surf6km(i)) THEN !check for NaN if this fails it's NaN
        Mean6kmwind(1) = Mean6kmwind(1)+ u_surf6km(i)
        k1 = k1+1
      !ENDIF
      !IF (v_surf6km(i).EQ.v_surf6km(i)) THEN  !check for NaN if this fails it's NaN
        Mean6kmwind(2) = Mean6kmwind(2)+ v_surf6km(i)
        k2 = k2+1
      !ENDIF
    enddo
    Mean6kmwind(1) = Mean6kmwind(1)/k1
    Mean6kmwind(2) = Mean6kmwind(2)/k2

    !#Calculate shear between mean between 0-500m and 5500-6000m.
    k1=0
    k2=0
    do i=12,13
      !IF (u_surf6km(i).EQ.u_surf6km(i)) THEN  !check for NaN if this fails it's NaN
        ushr_u = ushr_u+ u_surf6km(i)
        k1 = k1+1
      !ENDIF
      !IF (v_surf6km(i).EQ.v_surf6km(i)) THEN  !check for NaN if this fails it's NaN
        vshr_u = vshr_u+ v_surf6km(i)
        k2 = k2+1
      !ENDIF
    enddo
    ushr_u = ushr_u/k1
    vshr_u = vshr_u/k2

    k1=0
    k2=0
    do i=1,2
      !IF (u_surf6km(i).EQ.u_surf6km(i)) THEN  !check for NaN if this fails it's NaN
        ushr_d = ushr_d+ u_surf6km(i)
        k1 = k1+1
      !ENDIF
      !IF (v_surf6km(i).EQ.v_surf6km(i)) THEN  !check for NaN if this fails it's NaN
        vshr_d = vshr_d+ v_surf6km(i)
        k2 = k2+1
      !ENDIF
    enddo
    ushr_d = ushr_d/k1
    vshr_d = vshr_d/k2

    ushr = ushr_u-ushr_d
    vshr = vshr_u-vshr_d

    !#Calculate Left and Right Mover U and V components, Bunkers et al. [2000] method.
    RM(1) = Mean6kmwind(1) + 7.5*vshr/(ushr**2+vshr**2)**0.5
    RM(2) = Mean6kmwind(2) - 7.5*ushr/(ushr**2+vshr**2)**0.5
    LM(1) = Mean6kmwind(1) - 7.5*vshr/(ushr**2+vshr**2)**0.5
    LM(2) = Mean6kmwind(2) + 7.5*ushr/(ushr**2+vshr**2)**0.5

    return
    end subroutine bunkers_calc_pl

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    SUBROUTINE DINTERP2DZ(V3D,Z,LOC,V2D,N2,NZ)
    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(out) :: V2D

    INTEGER N2,NZ
    REAL, INTENT(IN)  :: V3D(NZ,N2)
    REAL, INTENT(OUT)  ::  V2D(N2)
    REAL, INTENT(IN)  ::  Z(NZ,N2)
    REAL, INTENT(IN)  ::  LOC

    REAL VMSG
    INTEGER J,KP,IP,IM
    LOGICAL INTERP
    REAL HEIGHT,W1,W2

    HEIGHT = LOC
    VMSG = -999999.
    ! does vertical coordinate increase or decrease with increasing k?
    ! set offset appropriately

    IP = 0
    IM = 1
    IF (Z(1,1).GT.Z(NZ,1)) THEN
        IP = 1
        IM = 0
    END IF
    !print *,'  V3D,Z,LOC,V2D,N2,NZ in interp = ',V3D,Z,LOC,V2D,N2,NZ

    DO J = 1,N2
    ! Initialize to missing.  Was initially hard-coded to -999999.
            V2D(J) = VMSG
            INTERP = .false.
            KP = NZ
            DO WHILE ((.NOT.INTERP) .AND. (KP.GE.2))
              IF (((Z(KP-IM,J).LE.HEIGHT).AND. (Z(KP-IP,J).GT.HEIGHT))) THEN
                    W2 = (HEIGHT-Z(KP-IM,J))/(Z(KP-IP,J)-Z(KP-IM,J))
                    W1 = 1.D0 - W2
                    V2D(J) = W1*V3D(KP-IM,J) + W2*V3D(KP-IP,J)
                    INTERP = .true.
              ENDIF
              KP = KP - 1
            ENDDO
    ENDDO
    !print *,'  V2D in interp = ',V2D

    RETURN
    END SUBROUTINE DINTERP2DZ
