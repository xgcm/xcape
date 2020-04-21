!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
    subroutine loopcape_pl(p3d_in , t3d_in , td3d_in, ps1d_in , ts1d_in , tds1d_in, &
                        &pinc, source, ml_depth, adiabat, start_3d, &
                        &nk_in, n2, cape3d, cin3d, MUlvl3d, z_out3d)
!-----------------------------------------------------------------------
!  loopcape - loop along n2 dimension to calculate Convective Available
!            Potential Energy (CAPE) from a sounding.
!
!  Author:  Chiara Lepore @chiaral
!
!  Disclaimer:  This code is made available WITHOUT WARRANTY.
!-----------------------------------------------------------------------
!  nk = number of pressure levels
!  n2 = number of grid point for which calculate CAPE
!  User options:

!  pinc   ! Pressure increment (Pa)
          ! (smaller number yields more accurate
          !  results,larger number makes code
          !  go faster) good value 1000

!  source     ! Source parcel:
              ! 1 = surface
              ! 2 = most unstable (max theta-e)
              ! 3 = mixed-layer (specify ml_depth)

! ml_depth  ! depth (m) of mixed layer for source=3 - good value 500m

!  adiabat    ! Formulation of moist adiabat:
              ! 1 = pseudoadiabatic, liquid only
              ! 2 = reversible, liquid only
              ! 3 = pseudoadiabatic, with ice
              ! 4 = reversible, with ice

! type_grid   ! type of vertical grid
              ! 1 = model level grid - uses all values in 3d files
              ! 2 = fixed pressure level grid
              !                 - p3d_in is 1d
              !                 - values below ps height needs to be deleted
!-----------------------------------------------------------------------

    implicit none

    !f2py threadsafe
    !f2py intent(out) :: cape3d,cin3d,z_out3d,MUlvl3d

    integer, intent(in) :: nk_in,n2
    integer :: i,j
    real, dimension(nk_in,n2), intent(in) :: p3d_in , t3d_in , td3d_in
    real, dimension(n2), intent(in) :: ps1d_in , ts1d_in , tds1d_in
    real, dimension(n2), intent(out) :: cape3d,cin3d,z_out3d
    integer, dimension(n2), intent(out) :: MUlvl3d
    real, intent(in) :: pinc ! Pressure increment (Pa) for interpolation
    integer, intent(in) :: source ! Source parcel
    real, intent(in) :: ml_depth ! depth (m) of mixed layer,for source=3
    integer, intent(in) :: adiabat   ! Formulation of moist adiabat:
    integer, dimension(n2), intent(in) ::  start_3d  ! for type_grid = 1 : equal to 1
                                                               ! for type_grid = 2 " array of position of
                                                               ! valid values in pressure levels
                                                               ! (where p3d_in <= ps_in )
    integer :: nk_start, nk_pl_in


    do i = 1, n2
      IF(ts1d_in(i).gt.0.0)THEN
        ! use only levels where p3d_in <= ps_in
        nk_start = start_3d(i)
        ! compute number of used levels
        nk_pl_in = nk_in - nk_start + 1
        call getcape_pl(p3d_in(nk_start:nk_in,i) , t3d_in(nk_start:nk_in,i) , td3d_in(nk_start:nk_in,i), &
        &ps1d_in(i), ts1d_in(i), tds1d_in(i), &
        &pinc, source, ml_depth, adiabat, nk_pl_in, &
        &cape3d(i) , cin3d(i), MUlvl3d(i), z_out3d(i))
      ELSE
        cape3d(i) = 0.0
        cin3d(i)  = 0.0
        z_out3d(i)= 0.0
        MUlvl3d(i)= 0
      ENDIF
    enddo
    return
    end subroutine loopcape_pl
!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
    subroutine loopcape_pl1d(t3d_in , td3d_in, p3d_in , ps1d_in , ts1d_in , tds1d_in, &
                        &pinc, source, ml_depth, adiabat, start_3d, &
                        &nk_in, n2, cape3d, cin3d, MUlvl3d, z_out3d)
!-----------------------------------------------------------------------
!  loopcape - loop along n2 dimension to calculate Convective Available
!            Potential Energy (CAPE) from a sounding.
!
!  Author:  Chiara Lepore @chiaral
!
!  Disclaimer:  This code is made available WITHOUT WARRANTY.
!-----------------------------------------------------------------------
!  nk = number of pressure levels
!  n2 = number of grid point for which calculate CAPE
!  User options:

!  pinc   ! Pressure increment (Pa)
          ! (smaller number yields more accurate
          !  results,larger number makes code
          !  go faster) good value 1000

!  source     ! Source parcel:
              ! 1 = surface
              ! 2 = most unstable (max theta-e)
              ! 3 = mixed-layer (specify ml_depth)

! ml_depth  ! depth (m) of mixed layer for source=3 - good value 500m

!  adiabat    ! Formulation of moist adiabat:
              ! 1 = pseudoadiabatic, liquid only
              ! 2 = reversible, liquid only
              ! 3 = pseudoadiabatic, with ice
              ! 4 = reversible, with ice

! type_grid   ! type of vertical grid
              ! 1 = model level grid - uses all values in 3d files
              ! 2 = fixed pressure level grid
              !                 - p3d_in is 1d
              !                 - values below ps height needs to be deleted
!-----------------------------------------------------------------------

    implicit none

    !f2py threadsafe
    !f2py intent(out) :: cape3d,cin3d,z_out3d,MUlvl3d

    integer, intent(in) :: nk_in,n2
    integer :: i,j
    real, dimension(nk_in,n2), intent(in) ::t3d_in , td3d_in
    real, dimension(nk_in,1), intent(in) :: p3d_in
    real, dimension(n2), intent(in) :: ps1d_in , ts1d_in , tds1d_in
    real, dimension(n2), intent(out) :: cape3d,cin3d,z_out3d
    integer, dimension(n2), intent(out) :: MUlvl3d
    real, intent(in) :: pinc ! Pressure increment (Pa) for interpolation
    integer, intent(in) :: source ! Source parcel
    real, intent(in) :: ml_depth ! depth (m) of mixed layer,for source=3
    integer, intent(in) :: adiabat   ! Formulation of moist adiabat:
    integer, dimension(n2), intent(in) ::  start_3d  ! for type_grid = 1 : equal to 1
                                                               ! for type_grid = 2 " array of position of
                                                               ! valid values in pressure levels
                                                               ! (where p3d_in <= ps_in )
    integer :: nk_start, nk_pl_in


    do i = 1, n2
      IF(ts1d_in(i).gt.0.0)THEN
        ! use only levels where p3d_in <= ps_in
        nk_start = start_3d(i)
        ! compute number of used levels
        nk_pl_in = nk_in - nk_start + 1
        call getcape_pl(p3d_in(nk_start:nk_in,1) , t3d_in(nk_start:nk_in,i) , td3d_in(nk_start:nk_in,i), &
        &ps1d_in(i), ts1d_in(i), tds1d_in(i), &
        &pinc, source, ml_depth, adiabat, nk_pl_in, &
        &cape3d(i) , cin3d(i), MUlvl3d(i), z_out3d(i))
      ELSE
        cape3d(i) = 0.0
        cin3d(i)  = 0.0
        z_out3d(i)= 0.0
        MUlvl3d(i)= 0
      ENDIF
    enddo
    return
    end subroutine loopcape_pl1d
!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    subroutine getcape_pl(p_in_A , t_in_A , td_in_A, ps_in_B , ts_in_B , tds_in_B, &
      &pinc, source, ml_depth, adiabat, nk_in, cape , cin, MUlvl, z_out)
    implicit none

    integer, intent(in) :: nk_in ! numbers of levels of 3d array
    real, dimension(nk_in), intent(in) :: p_in_A,t_in_A,td_in_A
    real, dimension(1), intent(in) :: ps_in_B,ts_in_B,tds_in_B
    real, intent(out) :: cape,cin,z_out
    integer, intent(out) :: MUlvl
    real, intent(in) :: pinc ! Pressure increment (Pa) for interpolation
    integer, intent(in) :: source ! Source parcel
    real, intent(in) :: ml_depth ! depth (m) of mixed layer,for source=3
    integer, intent(in) :: adiabat   ! Formulation of moist adiabat:
!-----------------------------------------------------------------------
!
!  getcape - a fortran90 subroutine to calculate Convective Available
!            Potential Energy (CAPE) from a sounding.
!
!  Version 1.02                           Last modified:  10 October 2008
!
!  Author:  George H. Bryan
!           Mesoscale and Microscale Meteorology Division
!           National Center for Atmospheric Research
!           Boulder, Colorado, USA
!           gbryan@ucar.edu
!
!  Disclaimer:  This code is made available WITHOUT WARRANTY.
!
!  References:  Bolton (1980, MWR, p. 1046) (constants and definitions)
!               Bryan and Fritsch (2004, MWR, p. 2421) (ice processes)
!
!-----------------------------------------------------------------------
!
!  Input:     nk - number of levels in the sounding (integer)
!
!           p_in - one-dimensional array of pressure (mb) (real)
!
!           t_in - one-dimensional array of temperature (C) (real)
!
!          td_in - one-dimensional array of dewpoint temperature (C) (real)
!
!           z - one dimensional array of heights above ground level (m) (real)
!  Output:  cape - Convective Available Potential Energy (J/kg) (real)
!
!            cin - Convective Inhibition (J/kg) (real)
!
!-----------------------------------------------------------------------

! FROM ORIGINAL getcape.F available at 
! https://www2.mmm.ucar.edu/people/bryan/Code/getcape.F
! Â©2019 - University Corporation for Atmospheric Research
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is furnished
! to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
!
!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!            No need to modify anything below here:
!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    logical :: doit,ice,cloud,not_converged
    integer :: k,kmax,n,nloop,i,orec

    real, dimension(nk_in+1) :: p,t,td,pi,q,z,th,thv,pt,pb,pc,pn,ptv ! use nk_in +1 in place of nk

    real :: the,maxthe,parea,narea,lfc
    real :: th1,p1,t1,qv1,ql1,qi1,b1,pi1,thv1,qt,dp,dz,ps,frac
    real :: th2,p2,t2,qv2,ql2,qi2,b2,pi2,thv2
    real :: thlast,fliq,fice,tbar,qvbar,qlbar,qibar,lhv,lhs,lhf,rm,cpm
    real*8 :: avgth,avgqv
    real :: getqvs,getqvi,getthe

!-----------------------------------------------------------------------

    real, parameter :: g     = 9.81
    real, parameter :: p00   = 100000.0
    real, parameter :: cp    = 1005.7
    real, parameter :: rd    = 287.04
    real, parameter :: rv    = 461.5
    real, parameter :: xlv   = 2501000.0
    real, parameter :: xls   = 2836017.0
    real, parameter :: t0    = 273.15
    real, parameter :: cpv   = 1875.0
    real, parameter :: cpl   = 4190.0
    real, parameter :: cpi   = 2118.636
    real, parameter :: lv1   = xlv+(cpl-cpv)*t0
    real, parameter :: lv2   = cpl-cpv
    real, parameter :: ls1   = xls+(cpi-cpv)*t0
    real, parameter :: ls2   = cpi-cpv

    real, parameter :: rp00  = 1.0/p00
    real, parameter :: eps   = rd/rv
    real, parameter :: reps  = rv/rd
    real, parameter :: rddcp = rd/cp
    real, parameter :: cpdrd = cp/rd
    real, parameter :: cpdg  = cp/g

    real, parameter :: converge = 0.0002

    integer, parameter :: debug_level =   0

!-----------------------------------------------------------------------
    real, dimension(nk_in+1) :: p_in,t_in,td_in
    integer :: nk  ! number of levels of 3d array PLUS surface level

    nk = nk_in + 1 ! number of levels of 3d array PLUS surface level
    p_in(1) = ps_in_B(1)
    t_in(1) = ts_in_B(1)
    td_in(1) = tds_in_B(1)
    p_in(2:nk) = p_in_A
    t_in(2:nk) = t_in_A
    td_in(2:nk) = td_in_A
    !print *,'  ts_in_B = ',ts_in_B
    !print *,'  nk_in, t_in_A = ', nk_in, t_in_A
    !print *,'  nk, t_in = ',nk, t_in

!---- convert p,t,td to mks units; get pi,q,th,thv ----!

    do k=1,nk
        p(k) = 100.0*p_in(k)
        t(k) = 273.15+t_in(k)
       td(k) = 273.15+td_in(k)
       pi(k) = (p(k)*rp00)**rddcp
        q(k) = getqvs(p(k),td(k))
       th(k) = t(k)/pi(k)
      thv(k) = th(k)*(1.0+reps*q(k))/(1.0+q(k))
    enddo

! ---- get height using the hydrostatic equation ----!

      z(1) = 0.0
      do k=2,nk
        dz = -cpdg*0.5*(thv(k)+thv(k-1))*(pi(k)-pi(k-1))
        z(k) = z(k-1) + dz
      enddo

! ---- Set MUlvl and z_out to -999999 for source = 1,3  ----!

  MUlvl = -999999
  z_out = -999999.

!---- find source parcel ----!

  IF(source.eq.1)THEN
    ! use surface parcel
    kmax = 1

  ELSEIF(source.eq.2)THEN
    ! use most unstable parcel (max theta-e)

    IF(p(1).lt.50000.0)THEN
      ! first report is above 500 mb ... just use the first level reported
      kmax = 1
      maxthe = getthe(p(1),t(1),td(1),q(1))
    ELSE
      ! find max thetae below 500 mb
      maxthe = 0.0
      do k=1,nk
        if(p(k).ge.50000.0)then
          the = getthe(p(k),t(k),td(k),q(k))
          if( the.gt.maxthe )then
            MUlvl = k
            maxthe = the
            kmax = k
          endif
        endif
      enddo
    ENDIF
    if(debug_level.ge.100) print *,'  kmax,maxthe = ',kmax,maxthe

  ELSEIF(source.eq.3)THEN
    ! use mixed layer

    IF( (z(2)-z(1)).gt.ml_depth )THEN
      ! the second level is above the mixed-layer depth:  just use the
      ! lowest level

      avgth = th(1)
      avgqv = q(1)
      kmax = 1

    ELSEIF( z(nk).lt.ml_depth )THEN
      ! the top-most level is within the mixed layer:  just use the
      ! upper-most level

      avgth = th(nk)
      avgqv = q(nk)
      kmax = nk

    ELSE
      ! calculate the mixed-layer properties:

      avgth = 0.0
      avgqv = 0.0
      k = 2
      if(debug_level.ge.100) print *,'  ml_depth = ',ml_depth
      if(debug_level.ge.100) print *,'  k,z,th,q:'
      if(debug_level.ge.100) print *,1,z(1),th(1),q(1)

      do while( (z(k).le.ml_depth) .and. (k.le.nk) )

        if(debug_level.ge.100) print *,k,z(k),th(k),q(k)

        avgth = avgth + 0.5*(z(k)-z(k-1))*(th(k)+th(k-1))
        avgqv = avgqv + 0.5*(z(k)-z(k-1))*(q(k)+q(k-1))

        k = k + 1

      enddo

      th2 = th(k-1)+(th(k)-th(k-1))*(ml_depth-z(k-1))/(z(k)-z(k-1))
      qv2 =  q(k-1)+( q(k)- q(k-1))*(ml_depth-z(k-1))/(z(k)-z(k-1))

      if(debug_level.ge.100) print *,999,ml_depth,th2,qv2

      avgth = avgth + 0.5*(ml_depth-z(k-1))*(th2+th(k-1))
      avgqv = avgqv + 0.5*(ml_depth-z(k-1))*(qv2+q(k-1))

      if(debug_level.ge.100) print *,k,z(k),th(k),q(k)

      avgth = avgth/ml_depth
      avgqv = avgqv/ml_depth

      kmax = 1

    ENDIF

    if(debug_level.ge.100) print *,avgth,avgqv

  ELSE

    print *
    print *,'  Unknown value for source'
    print *
    print *,'  source = ',source
    print *
    stop

  ENDIF

!---- define parcel properties at initial location ----!
    narea = 0.0

  if( (source.eq.1).or.(source.eq.2) )then
    k    = kmax
    th2  = th(kmax)
    pi2  = pi(kmax)
    p2   = p(kmax)
    t2   = t(kmax)
    thv2 = thv(kmax)
    qv2  = q(kmax)
    b2   = 0.0
  elseif( source.eq.3 )then
    k    = kmax
    th2  = avgth
    qv2  = avgqv
    thv2 = th2*(1.0+reps*qv2)/(1.0+qv2)
    pi2  = pi(kmax)
    p2   = p(kmax)
    t2   = th2*pi2
    b2   = g*( thv2-thv(kmax) )/thv(kmax)
  endif

    ql2 = 0.0
    qi2 = 0.0
    qt  = qv2

    cape = 0.0
    cin  = 0.0
    lfc  = 0.0

    doit = .true.
    cloud = .false.
    if(adiabat.eq.1.or.adiabat.eq.2)then
      ice = .false.
    else
      ice = .true.
    endif

      the = getthe(p2,t2,t2,qv2)
      if(debug_level.ge.100) print *,'  the = ',the

!---- begin ascent of parcel ----!

      if(debug_level.ge.100)then
        print *,'  Start loop:'
        print *,'  p2,th2,qv2 = ',p2,th2,qv2
      endif

    do while( doit .and. (k.lt.nk) )

        k = k+1
       b1 =  b2

       dp = p(k-1)-p(k)

      if( dp.lt.pinc )then
        nloop = 1
      else
        nloop = 1 + int( dp/pinc )
        dp = dp/float(nloop)
      endif

      do n=1,nloop

         p1 =  p2
         t1 =  t2
        pi1 = pi2
        th1 = th2
        qv1 = qv2
        ql1 = ql2
        qi1 = qi2
        thv1 = thv2

        p2 = p2 - dp
        pi2 = (p2*rp00)**rddcp

        thlast = th1
        i = 0
        not_converged = .true.

        do while( not_converged )
          i = i + 1
          t2 = thlast*pi2
          if(ice)then
            fliq = max(min((t2-233.15)/(273.15-233.15),1.0),0.0)
            fice = 1.0-fliq
          else
            fliq = 1.0
            fice = 0.0
          endif
          qv2 = min( qt , fliq*getqvs(p2,t2) + fice*getqvi(p2,t2) )
          qi2 = max( fice*(qt-qv2) , 0.0 )
          ql2 = max( qt-qv2-qi2 , 0.0 )

          tbar  = 0.5*(t1+t2)
          qvbar = 0.5*(qv1+qv2)
          qlbar = 0.5*(ql1+ql2)
          qibar = 0.5*(qi1+qi2)

          lhv = lv1-lv2*tbar
          lhs = ls1-ls2*tbar
          lhf = lhs-lhv

          rm=rd+rv*qvbar
          cpm=cp+cpv*qvbar+cpl*qlbar+cpi*qibar
          th2=th1*exp(  lhv*(ql2-ql1)/(cpm*tbar)     &
                       +lhs*(qi2-qi1)/(cpm*tbar)     &
                       +(rm/cpm-rd/cp)*alog(p2/p1) )

          if(i.gt.90) print *,i,th2,thlast,th2-thlast
          if(i.gt.100)then
            print *
            print *,'  Error:  lack of convergence'
            print *
            print *,'  ... stopping iteration '
            print *
            cape = 0.0
            cin = 0.0
            return
          endif
          if( abs(th2-thlast).gt.converge )then
            thlast=thlast+0.3*(th2-thlast)
          else
            not_converged = .false.
          endif
        enddo

        ! Latest pressure increment is complete.  Calculate some
        ! important stuff:

        if( ql2.ge.1.0e-10 ) cloud = .true.

        IF(adiabat.eq.1.or.adiabat.eq.3)THEN
          ! pseudoadiabat
          qt  = qv2
          ql2 = 0.0
          qi2 = 0.0
        ELSEIF(adiabat.le.0.or.adiabat.ge.5)THEN
          print *
          print *,'  Undefined adiabat'
          print *
          stop 10000
        ENDIF

      enddo

      thv2 = th2*(1.0+reps*qv2)/(1.0+qv2+ql2+qi2)
        b2 = g*( thv2-thv(k) )/thv(k)
        dz = -cpdg*0.5*(thv(k)+thv(k-1))*(pi(k)-pi(k-1))

      the = getthe(p2,t2,t2,qv2)

      ! Get contributions to CAPE and CIN:

      if( (b2.ge.0.0) .and. (b1.lt.0.0) )then
        ! first trip into positive area
        ps = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1)
        frac = b2/(b2-b1)
        parea =  0.5*b2*dz*frac
        narea = narea-0.5*b1*dz*(1.0-frac)
        if(debug_level.ge.200)then
          print *,'      b1,b2 = ',b1,b2
          print *,'      p1,ps,p2 = ',p(k-1),ps,p(k)
          print *,'      frac = ',frac
          print *,'      parea = ',parea
          print *,'      narea = ',narea
        endif
        cin  = cin  + narea
        narea = 0.0
      elseif( (b2.lt.0.0) .and. (b1.gt.0.0) )then
        ! first trip into neg area
        ps = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1)
        frac = b1/(b1-b2)
        parea =  0.5*b1*dz*frac
        narea = -0.5*b2*dz*(1.0-frac)
        if(debug_level.ge.200)then
          print *,'      b1,b2 = ',b1,b2
          print *,'      p1,ps,p2 = ',p(k-1),ps,p(k)
          print *,'      frac = ',frac
          print *,'      parea = ',parea
          print *,'      narea = ',narea
        endif
      elseif( b2.lt.0.0 )then
        ! still collecting negative buoyancy
        parea =  0.0
        narea = narea-0.5*dz*(b1+b2)
      else
        ! still collecting positive buoyancy
        parea =  0.5*dz*(b1+b2)
        narea =  0.0
      endif

      cape = cape + max(0.0,parea)

      if(debug_level.ge.200)then
        write(6,102) p2,b1,b2,cape,cin,cloud
102     format(5(f13.4),2x,l1)
      endif

      if( (p(k).le.10000.0).and.(b2.lt.0.0) )then
        ! stop if b < 0 and p < 100 mb
        doit = .false.
      endif
      z_out=z(k)
    enddo

!---- All done ----!

    return
    end subroutine getcape_pl

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    real function getqvs(p,t)
    implicit none

    real :: p,t,es

    real, parameter :: eps = 287.04/461.5

    es = 611.2*exp(17.67*(t-273.15)/(t-29.65))
    getqvs = eps*es/(p-es)

    return
    end function getqvs

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    real function getqvi(p,t)
    implicit none

    real :: p,t,es

    real, parameter :: eps = 287.04/461.5

    es = 611.2*exp(21.8745584*(t-273.15)/(t-7.66))
    getqvi = eps*es/(p-es)

    return
    end function getqvi

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    real function getthe(p,t,td,q)
    implicit none

    real :: p,t,td,q
    real :: tlcl

    if( (td-t).ge.-0.1 )then
      tlcl = t
    else
      tlcl = 56.0 + ( (td-56.0)**(-1) + 0.00125*alog(t/td) )**(-1)
    endif

    getthe=t*( (100000.0/p)**(0.2854*(1.0-0.28*q)) )   &
            *exp( ((3376.0/tlcl)-2.54)*q*(1.0+0.81*q) )

    return
    end function getthe

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
