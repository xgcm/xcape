!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
      SUBROUTINE loop_sreh_ml(u, v, aglh, us, vs, aglhs, &
                          &cu_rm, cv_rm, cu_lm, cv_lm, top, nk, n2, &
                          &sreh_rm, sreh_lm)
!-----------------------------------------------------------------------
!  loop_sreh_ml - loop along n2 dimension to calculate Storm Relative
!            Helicity (SRH) for model/sigma level data
!
!  Author:  Chiara Lepore @chiaral
!
!  Disclaimer:  This code is made available WITHOUT WARRANTY.
!-----------------------------------------------------------------------
!  Input:
!  nk = number of pressure levels
!  n2 = number of grid point for which calculate CAPE
!  u, v, us, vs = zonal and meridional winds (vertical and surface)
!  aglh, aglhs = standard above ground height for vertical and surface levels
!  cu_rm, cv_rm, cu_lm, cv_lm = zonal and meridional, right moving and 
!                               left moving storm flow
!  top = height for which srh is calculated
!  Output:
!  sreh_rm, sreh_lm = storm relative helicity right and left moving
!-----------------------------------------------------------------------                          
                          
      IMPLICIT NONE

      !f2py threadsafe
      !f2py intent(out) :: sreh

      INTEGER, INTENT(IN) :: nk, n2
      REAL(KIND=8), DIMENSION(nk ,n2), INTENT(IN) :: u, v,aglh
      REAL(KIND=8), INTENT(IN) :: top
      REAL(KIND=8), DIMENSION(n2), INTENT(IN) :: us, vs, aglhs
      REAL(KIND=8), DIMENSION(n2), INTENT(IN) :: cu_rm, cv_rm, cu_lm, cv_lm 
      REAL(KIND=8), DIMENSION(n2), INTENT(OUT) :: sreh_rm, sreh_lm
      INTEGER :: i, nk_all
      REAL(KIND=8), DIMENSION(nk+1) :: u_all, v_all, aglh_all

      nk_all = nk+1
      do i = 1, n2
          u_all(1) = us(i)
          v_all(1) = vs(i)
          aglh_all(1) = aglhs(i)
          u_all(2:nk_all) = u(:,i)
          v_all(2:nk_all) = v(:,i)
          aglh_all(2:nk_all) = aglh(:,i)

          call DCALRELHL_ml(u_all, v_all, aglh_all, &
          &cu_rm(i), cv_rm(i), cu_lm(i), cv_lm(i), top, &
          &nk_all, sreh_rm(i), sreh_lm(i) )
      enddo
      return
      end subroutine loop_sreh_ml

      !-----------------------------------------------------------------------
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !-----------------------------------------------------------------------

      SUBROUTINE DCALRELHL_ml(u, v, aglh, cu_rm, cv_rm, cu_lm, cv_lm, top, nk, &
                              &sreh_rm, sreh_lm )

      IMPLICIT NONE

      !f2py threadsafe
      !f2py intent(out) :: sreh

      INTEGER, INTENT(IN) :: nk
      REAL(KIND=8), DIMENSION(nk), INTENT(IN) :: u, v,aglh
      REAL(KIND=8), INTENT(IN) :: top
      REAL(KIND=8), DIMENSION(1), INTENT(IN) :: cu_rm, cv_rm, cu_lm, cv_lm 
      REAL(KIND=8), DIMENSION(1), INTENT(OUT) :: sreh_rm, sreh_lm
      REAL(KIND=8) :: INTERP1



      REAL(KIND=8), DIMENSION(1) :: x, sum
      REAL(KIND=8), DIMENSION(nk) :: utop,vtop
      INTEGER :: k, ktop

      utop = u
      vtop = v
      ktop = 0
      DO k = 2 , nk
          IF (((aglh(k)) .GT. top) .AND. &
             (ktop .EQ.0)) THEN
              ktop = k
              ! Remember the arrays are inverted hence ktop+1
              utop(ktop) = INTERP1(u(ktop),u(ktop-1), &
                     aglh(ktop),top,aglh(ktop-1))
              vtop(ktop) = INTERP1(v(ktop),v(ktop-1), &
                     aglh(ktop),top,aglh(ktop-1))

          ENDIF
      END DO
      !u(ktop,j)=utop(ktop,j)
      !v(ktop,j)=vtop(ktop,j)

      sum = 0.D0
      DO k = 2, ktop
          x = ((utop(k) - cu_rm)*(vtop(k) - vtop(k-1))) - &
              ((vtop(k) - cv_rm)*(utop(k) - utop(k-1)))
          sum = sum + x
      END DO
      sreh_rm = -sum

      sum = 0.D0
      DO k = 2, ktop
          x = ((utop(k) - cu_lm)*(vtop(k) - vtop(k-1))) - &
              ((vtop(k) - cv_lm)*(utop(k) - utop(k-1)))
          sum = sum + x
      END DO
      sreh_lm = -sum

      RETURN

      END SUBROUTINE DCALRELHL_ml

!-------Functions - Interpolation
!-------Given Y1 at X1, Y3 at X3, and X2, Calculate Y2(Interp) at X2

      REAL(KIND=8) FUNCTION INTERP1(Y1,Y3,X1,X2,X3)
      IMPLICIT NONE
      REAL(KIND=8) :: Y1,Y3,X1,X2,X3
      IF (X3.EQ.X1) X1=X1-0.01
      INTERP1 = Y1+((Y3-Y1)*((X2-X1)/(X3-X1)))
      RETURN
      END FUNCTION INTERP1
