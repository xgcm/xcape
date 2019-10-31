!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
    subroutine interp_loop_pl(INPUT, Y_TO_USE, INPUTs, Y_TO_USEs, LOC_in, start_3d, nk, n2, OUTPUT)
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
    !f2py intent(out) :: OUTPUT

    integer, intent(in) :: nk,n2
    integer :: i, nk_all, nk_pl_in, nk_start

    real, dimension(nk,n2), intent(in) :: INPUT, Y_TO_USE
    real, dimension(n2), intent(in) :: INPUTs, Y_TO_USEs
    real, dimension(1), intent(in) :: LOC_in
    
    real, dimension(n2), intent(out) :: OUTPUT
    real, dimension(nk+1) :: INPUT_all, Y_TO_USE_all

    do i = 1, n2
        nk_start = start_3d(i)
        ! compute number of used levels in 3d
        nk_pl_in = nk - nk_start + 1
        ! compute number of used levels in 3d + surface
        nk_all = nk_pl_in+1

        INPUT_all(1) = INPUTs(i)
        Y_TO_USE_all(1) = Y_TO_USEs(i)
        INPUT_all(2:nk_all) = INPUT(nk_start:nk,i)
        Y_TO_USE_all(2:nk_all) = Y_TO_USE(nk_start:nk,i)

        call DINTERP_DZ(INPUT_all,Y_TO_USE_all,LOC_in,OUTPUT(i),1,nk_all)
    enddo
    return
    end subroutine Interp_loop_pl

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
    subroutine interp_loop_pl_Y1d(INPUT, Y_TO_USE, INPUTs, Y_TO_USEs, LOC_in, start_3d, nk, n2, OUTPUT)
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
    !f2py intent(out) :: OUTPUT

    integer, intent(in) :: nk,n2
    integer :: i, nk_all, nk_pl_in, nk_start

    real, dimension(nk,n2), intent(in) :: INPUT, 
    real, dimension(nk_in,1), intent(in) :: Y_TO_USE
    real, dimension(n2), intent(in) :: INPUTs, Y_TO_USEs
    real, dimension(1), intent(in) :: LOC_in
    
    real, dimension(n2), intent(out) :: OUTPUT
    real, dimension(nk+1) :: INPUT_all, Y_TO_USE_all

    do i = 1, n2
        nk_start = start_3d(i)
        ! compute number of used levels in 3d
        nk_pl_in = nk - nk_start + 1
        ! compute number of used levels in 3d + surface
        nk_all = nk_pl_in+1

        INPUT_all(1) = INPUTs(i)
        Y_TO_USE_all(1) = Y_TO_USEs(i)
        INPUT_all(2:nk_all) = INPUT(nk_start:nk,i)
        Y_TO_USE_all(2:nk_all) = Y_TO_USE(nk_start:nk,1)

        call DINTERP_DZ(INPUT_all,Y_TO_USE_all,LOC_in,OUTPUT(i),1,nk_all)
    enddo
    return
    end subroutine Interp_loop_pl

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    SUBROUTINE DINTERP_DZ(V3D,Z,LOC,V2D,N2,NZ)
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
    END SUBROUTINE DINTERP_DZ
