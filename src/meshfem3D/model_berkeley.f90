!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!--------------------------------------------------------------------------------------------------
! SEMUCB - WM1
!
! A global radially anisotropic shear-wave speed model developed by French & Romanowicz [2014].
!
! reference:
!   Whole-mantle radially anisotropic shear velocity structure from spectral-element waveform tomography
!   S. W. French,   B. A. Romanowicz
!   Geophysical Journal International, Volume 199, Issue 3, December 2014,
!   Pages 1303--1327, https://doi.org/10.1093/gji/ggu334
!
! note: uses 1D Berkeley model as reference model together with a smooth crustal model
!--------------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------------
! Berkeley 3D model data
!--------------------------------------------------------------------------------------------------

module model_berkeley_par

  implicit none

  ! spline arrays
  double precision, dimension(:), allocatable :: aknot,oknot,aknot2,oknot2
  double precision, dimension(:), allocatable :: mdl
  integer, dimension(:), allocatable :: level,level2

  ! spline node radii
  real, dimension(:), allocatable :: kntrad,kntrad_hh     ! real arrays - to avoid conversion in fspl() to float

  ! parameterization
  integer, parameter :: MAXPAR = 4
  integer, dimension(MAXPAR) :: nknotA1 = 0,nknotA2 = 0

  character(1), dimension(:), allocatable :: parblock
  integer :: npar,ndisc,surface,NBPARAM

  logical :: hknots2_exist = .false.,unconformal = .false.

  ! radius of moho discontinuity from reference 1D model (in km)
  double precision :: moho1D_radius = -1.0

  ! work arrays
  double precision, dimension(:), allocatable :: work_dh
  integer, dimension(:), allocatable :: work_kindex

end module model_berkeley_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_berkeley_broadcast()

! standard routine to setup model

  use model_berkeley_par
  use constants, only: A3d_folder,myrank,IMAIN,EARTH_R_KM

  implicit none

  integer,parameter :: unit1 = 51,unit2 = 52
  integer :: dum2,i,j,k,n,mdim,ier
  integer :: size_work
  character :: trash
  character(len=100), parameter :: A3d_dat            = trim(A3d_folder) // 'A3d.dat'
  character(len=100), parameter :: hknots_dat         = trim(A3d_folder) // 'hknots.dat'
  character(len=100), parameter :: hknots2_dat        = trim(A3d_folder) // 'hknots2.dat'

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) 'broadcast model: SEMUCB Berkeley'
    call flush_IMAIN()
  endif

  ! determine moho radius from 1D reference model
  if (myrank == 0) then
    ! gets exact 1D moho radius (in km)
    call determine_1dberkeley_moho_radius(moho1D_radius)

    ! adjust radius slightly to be in mantle
    moho1D_radius = moho1D_radius - 0.1d0

    ! use non-dimensionalized moho radius for comparisons
    moho1D_radius = moho1D_radius / EARTH_R_KM

    !debug
    !print *,'debug: [model_berkeley_broadcast] moho radius = ',moho1D_radius
  endif
  ! broadcasts to all other processes
  call bcast_all_singledp(moho1D_radius)

  !
  ! reads 3D model
  !
  if (myrank == 0) then
    ! spline knots
    open(unit1,file=trim(hknots_dat),status='old',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(hknots_dat)
      stop 'Error opening file hknots.dat'
    endif

    read(unit1,*) nknotA2(1)

    if (allocated(oknot) .or. allocated(aknot) .or. allocated(level)) then
      print *,'[model_berkeley_broadcast] A3d already initiated'
      return
    endif

    allocate(oknot(nknotA2(1)),aknot(nknotA2(1)),level(nknotA2(1)),stat=ier)
    if (ier /= 0) stop 'Error allocating oknot,.. arrays'
    oknot(:) = 0.d0
    aknot(:) = 0.d0
    level(:) = 0

    do i = 1,nknotA2(1)
      read(unit1,*) oknot(i),aknot(i),level(i)
    enddo
    close(unit1)

    inquire(file=trim(hknots2_dat),exist=hknots2_exist)

    if (hknots2_exist) then
      open(unit1,file=trim(hknots2_dat),status='old',iostat=ier)
      if (ier /= 0) then
        print *,'Error opening file: ',trim(hknots2_dat)
        stop 'Error opening file hknots2.dat'
      endif

      read(unit1,*) nknotA2(2)

      if (allocated(oknot2) .or. allocated(aknot2) .or. allocated(level2)) then
        print *,'[model_berkeley_broadcast] A3d hnots2 already initiated'
        return
      endif

      allocate(oknot2(nknotA2(2)),aknot2(nknotA2(2)),level2(nknotA2(2)),stat=ier)
      if (ier /= 0) stop 'Error allocating oknot2,.. arrays'
      oknot2(:) = 0.d0
      aknot2(:) = 0.d0
      level2(:) = 0

      do i = 1,nknotA2(2)
        read(unit1,*) oknot2(i),aknot2(i),level2(i)
      enddo
      close(unit1)
    else
      nknotA2(2) = nknotA2(1)
    endif

    !open model file and read in model
    open(unit2,file=trim(A3d_dat),status='old',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(A3d_dat)
      stop 'Error opening file A3d.dat'
    endif

    read(unit2,*) npar

    NBPARAM = npar
    if (npar > MAXPAR)   stop 'npar greater than MAXPAR'

    allocate(parblock(npar),stat=ier)
    if (ier /= 0) stop 'Error allocating parblock array'
    parblock(:) = ''

    do i = 1,npar
      read(unit2,*) dum2,nknotA1(i),parblock(i)

      ! same number of splines for all parameters
      if (i > 1 .and. nknotA1(i) /= nknotA1(1)) then
        stop 'Inconsistent A1 splines between parameters'
      endif

      if (i > 2) then
        nknotA2(i) = dum2
        if (nknotA2(i) /= nknotA2(2)) stop 'Param 3 and 4 need the same A2 splines than param 2'
      else if (dum2 /= nknotA2(i)) then
        stop 'Inconsistent hknots.dat and A3d.dat'
      endif
      if (i == 2 .and. nknotA2(i) /= nknotA2(1)) then
        unconformal = .true.
        if (.not. hknots2_exist) stop 'unconformal grid requires hknots2.dat'
      endif
    enddo

    read(unit2,*) ndisc

    surface = 0
    if (ndisc > 0) then
      surface = ndisc
      if (unconformal) print *,'discontinuities assumed same grid as first par'
      do i = 1,surface
        read(unit2,*) dum2, trash
      enddo
    endif

    allocate(kntrad(nknotA1(1)), &
             kntrad_hh(nknotA1(1)-1),stat=ier)
    if (ier /= 0) stop 'Error allocating kntrad,.. arrays'
    kntrad(:) = 0.e0
    kntrad_hh(:) = 0.e0

    read(unit2,*) (kntrad(i),i=1,nknotA1(1))

    ! takes spacings between spline radii
    ! (used for spline evaluations fspl(..) in spl_A3d.c, see routine fill_hh_A3d())
    do i = 1,nknotA1(1)-1
      kntrad_hh(i) = kntrad(i+1) - kntrad(i)
    enddo

    mdim = 0
    do i = 1,npar
      mdim = mdim + nknotA1(i) * nknotA2(i)
    enddo
    mdim = mdim + ndisc * nknotA2(1)

    allocate(mdl(mdim),stat=ier)
    if (ier /= 0) stop 'Error allocating mdl array'
    mdl(:) = 0.d0

    n = 0
    do i = 1,npar
      do j = 1,nknotA1(i)
        read(unit2,*) (mdl(k+n),k=1,nknotA2(i))
        n = n + nknotA2(i)
      enddo
    enddo
    do i = 1,ndisc
      read(unit2,*) (mdl(k+n),k=1,nknotA2(1))
      n = n + nknotA2(1)
    enddo

    if (n /= mdim) stop 'init_A3d dimension error'
    close(unit2)

  endif

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) '  number of parameters: NBPARAM = ',NBPARAM
    call flush_IMAIN()
  endif

  !
  ! broadcast
  !
  !
  call bcast_all_singlei(nknotA2(1))

  if (.not. allocated(oknot)) allocate(oknot(nknotA2(1)))
  if (.not. allocated(aknot)) allocate(aknot(nknotA2(1)))
  if (.not. allocated(level)) allocate(level(nknotA2(1)))

  call bcast_all_i(level,nknotA2(1))
  call bcast_all_dp(oknot,nknotA2(1))
  call bcast_all_dp(aknot,nknotA2(1))

  call bcast_all_singlel(hknots2_exist)

  if (hknots2_exist) then
    call bcast_all_singlei(nknotA2(2))

    if (.not. allocated(oknot2)) allocate(oknot2(nknotA2(2)))
    if (.not. allocated(aknot2)) allocate(aknot2(nknotA2(2)))
    if (.not. allocated(level2)) allocate(level2(nknotA2(2)))

    call bcast_all_i(level2,nknotA2(2))
    call bcast_all_dp(oknot2,nknotA2(2))
    call bcast_all_dp(aknot2,nknotA2(2))
  else
    nknotA2(2) = nknotA2(1)
  endif

  call bcast_all_singlei(npar)
  NBPARAM = npar

  if (.not. allocated(parblock)) allocate(parblock(npar))

  call bcast_all_ch_array(parblock,npar,1)
  call bcast_all_i(nknotA1,MAXPAR)
  call bcast_all_i(nknotA2,MAXPAR)
  call bcast_all_singlel(unconformal)
  call bcast_all_singlei(ndisc)

  ! spline radii
  if (.not. allocated(kntrad)) allocate(kntrad(nknotA1(1)))
  if (.not. allocated(kntrad_hh)) allocate(kntrad_hh(nknotA1(1)-1))

  call bcast_all_r(kntrad,nknotA1(1))
  call bcast_all_r(kntrad_hh,nknotA1(1)-1)

  call bcast_all_singlei(mdim)

  if (.not. allocated(mdl)) allocate(mdl(mdim))

  call bcast_all_dp(mdl,mdim)

  ! allocates temporary work arrays for model_berkeley_shsv() routine
  ! to avoid re-allocations of the same arrays for each call
  size_work = maxval(nknotA2(:))
  allocate(work_dh(size_work),work_kindex(size_work),stat=ier)
  if (ier /= 0) stop 'Error allocating work arrays for Berkeley model'
  work_dh(:) = 0.d0
  work_kindex(:) = 0

  end subroutine model_berkeley_broadcast


!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_berkeley_shsv(r,theta,phi,dvsh,dvsv,dvph,dvpv,drho,eta_aniso,iregion_code,CRUSTAL)

! returns isotropic vs, vp, and rho assuming scaling dlnVs/dlnVp=2 dlnVs/dlnrho=3
! also returns anisotropic parameters xi,fi,eta,Gc,Gs,Hc,Hs,Bc,Bs if ifanis=1

  use model_berkeley_par
  use constants

  implicit none

  double precision, intent(in) :: r,theta,phi
  double precision, intent(out) :: dvsv,dvsh,dvpv,dvph,drho,eta_aniso

  integer, intent(in) :: iregion_code
  logical, intent(in) :: CRUSTAL

  ! local parameters
  double precision :: x,rho1d,vpv1d,vph1d,vsv1d,vsh1d,eta1d,Qmu1d,Qkappa1d

  double precision :: vs,vp,rho,Qs
  double precision :: xi,fi,eta,Gc,Gs
  double precision :: fi_inv

  integer :: jump,effnknot,i,j,k
  !integer, dimension(:), allocatable :: kindex

  double precision :: lat,lon
  double precision :: del,dr,dv
  double precision :: AA,CC,FF,LL,NN
  double precision :: eta1,adel1,r_

  ! Include level 8 (FM, April 2021 - Mod. courtesy of Dan Frost)
  double precision, dimension(8), parameter :: adel = (/63.4, 31.7, 15.8, 7.9, 3.95, 1.98, 0.99, 0.495/)
  !adel=(/63.4, 31.7, 15.8, 7.9, 3.95, 1.98, 0.99/)

  !double precision, dimension(:), allocatable :: dh

  double precision,external :: getdel
  double precision,external :: spbsp

  double precision :: vsv,vsh,vpv,vph,scaleval
  double precision :: aa1, bb1

  double precision, parameter :: rad2deg = 180.d0/PI

  ! initializes model perturbations
  dvsv = 0.d0
  dvsh = 0.d0
  dvpv = 0.d0
  dvph = 0.d0
  drho = 0.d0
  eta_aniso = 1.d0

  xi = 1.d0
  fi = 1.d0

  !if (.not.model1D_exist) &
  !    stop 'no 1D model file'

  !if (ifanis_berk==1) then
  !    if (.not.(present(xi).and.present(fi).and.present(eta))) &
  !        stop 'A3d_full: ifanis_berk inconsistent'
  !endif

  ! limits radius to stay below 1D moho depth
  ! note: r is non-dimensionalized/normalized input radius between [0,1]
  if (r > moho1D_radius) then
    r_ = moho1D_radius  ! * EARTH_R_KM
  else
    r_ = r              ! * EARTH_R_KM
  endif

  ! non-dimensionalized radius
  x = r_                ! / EARTH_R_KM

  call model_1dberkeley(x,rho1d,vpv1d,vph1d,vsv1d,vsh1d,eta1d,Qkappa1d,Qmu1d,iregion_code,CRUSTAL)

  if (rho1d < 1200.d0) then
    ! No water in RegSEM please
    x = r_ * EARTH_R_KM / 6367.999d0
    call model_1dberkeley(x,rho1d,vpv1d,vph1d,vsv1d,vsh1d,eta1d,Qkappa1d,Qmu1d,iregion_code,CRUSTAL)
  endif

  ! below use r_ as radius in km
  r_ = r_ * EARTH_R_KM

  ! re-dimensionalizes 1d model values
  scaleval = EARTH_R * dsqrt(PI*GRAV*EARTH_RHOAV)

  rho1d = rho1d * EARTH_RHOAV
  vpv1d = vpv1d * scaleval
  vph1d = vph1d * scaleval
  vsv1d = vsv1d * scaleval
  vsh1d = vsh1d * scaleval

  rho = rho1d
  AA  = vph1d*vph1d * rho
  CC  = vpv1d*vpv1d * rho
  LL  = vsv1d*vsv1d * rho
  NN  = vsh1d*vsh1d * rho
  FF  = eta1d * (AA - 2.d0 * LL)
  qs  = qmu1d

  !call get_1Dmodel_TI(AA,CC,FF,LL,NN,rho,Qs,r_*1000.d0)
  !if (rho < 1200.d0)   call get_1Dmodel_TI(AA,CC,FF,LL,NN,rho,Qs,6367999.d0)   ! No water in RegSEM please

  eta1 = FF / (AA - 2.d0 * LL)

  ! Voigt average
  vp = sqrt((3.d0 * CC + (8.d0 + 4.d0 * eta1) * AA + 8.d0 * (1.d0 - eta1) * LL) / (15.d0 * rho))
  vs = sqrt((CC + (1.d0 - 2.d0 * eta1) * AA + (6.d0 + 4.d0 * eta1) * LL + 5.d0 * NN) / (15.d0 * rho))

  eta = eta1
  xi = NN / LL
  fi = CC / AA

  ! point lat/lon in degrees
  lat = 90.d0 - rad2deg * theta
  lon = rad2deg * phi

  ! Vs perturbation
  if (r_ > kntrad(nknotA1(1)) .or. r_ < kntrad(1)) then
    dv = 0.d0
  else
    jump = 0
    !allocate(dh(nknotA2(1)),kindex(nknotA2(1)))
    work_dh(:) = 0.d0
    work_kindex(:) = 0

    effnknot = 0
    do i = 1,nknotA2(1)
      del = getdel(lat,lon,aknot(i),oknot(i))

      if (del <= adel(level(i)) * 2.d0) then
        effnknot = effnknot+1
        work_kindex(effnknot) = i
        work_dh(effnknot) = spbsp(del,adel(level(i)))
      endif
    enddo

    dv = 0.d0
    do i = 1,nknotA1(1)
      ! spline value
      call fspl(i,nknotA1(1),kntrad,kntrad_hh,r_,dr)

      do j = 1,effnknot
        dv = dv + dr * work_dh(j) * mdl(jump + work_kindex(j) + nknotA2(1) * (i-1))
      enddo
    enddo
    !deallocate(dh,kindex)
  endif

  ! Scaling
  vs = vs + dv * vs
  vp = vp + 0.5d0 * dv * vp
  rho = rho + 0.33333d0 * dv * rho

  ! Perturbation of other parameters
  ! note: default in A3d.dat is npar == 2
  if (npar > 1) then
    do k = 2,npar
      if (r_ > kntrad(nknotA1(k)) .or. r_ < kntrad(1)) then
        dv = 0.d0
      else
        jump = jump + nknotA1(k-1)*nknotA2(k-1)

        !allocate(dh(nknotA2(k)),kindex(nknotA2(k)))
        work_dh(:) = 0.d0
        work_kindex(:) = 0

        effnknot = 0
        do i = 1,nknotA2(k)
          if (unconformal) then
            del = getdel(lat,lon,aknot2(i),oknot2(i))
            adel1 = adel(level2(i))
          else
            del = getdel(lat,lon,aknot(i),oknot(i))
            adel1 = adel(level(i))
          endif
          if (del <= adel1 * 2.d0) then
            effnknot = effnknot+1
            work_kindex(effnknot) = i
            work_dh(effnknot) = spbsp(del,adel1)
          endif
        enddo

        dv = 0.d0
        do i = 1,nknotA1(k)
          ! spline value
          call fspl(i,nknotA1(k),kntrad,kntrad_hh,r_,dr)

          do j = 1,effnknot
            dv = dv + dr * work_dh(j) * mdl(jump + work_kindex(j) + nknotA2(k) * (i-1))
          enddo
        enddo
        !deallocate(dh,kindex)
      endif

      if (k == 2) then
        xi = xi + dv * xi
        fi = fi - 1.5d0 * dv * fi
        eta = eta - 2.5d0 * dv * eta
      else if (k == 3) then
        Gc = dv
      else if (k == 4) then
        Gs = dv
      endif
      ! Here we can add a scaling to get Hc, Hs, Bc and Bs
    enddo
  else
    ! no other parameters
    dv = 0.d0
  endif

  ! ============= commented by < FM> on Feb 3, 2020 =====
  !vsv = sqrt(3.d0/(xi+2.d0))*vs
  !vsh = sqrt(xi)*vsv
  !vph = sqrt(5.d0/(fi+4.d0))*vp
  !vpv = sqrt(fi)*vph
  !rho = rho

  ! ====================================================
  ! New conversion relationships < FM> - Feb 3, 2020
  ! Auxiliar values
  fi_inv = 1.d0 / fi

  aa1 = 3.d0 + ( 8.d0 + 4.d0 * eta ) * fi_inv
  bb1 = 1.d0 + ( 1.d0 - 2.d0 * eta ) * fi_inv
  vsv = sqrt( 15.d0 * (vp*vp * bb1 - vs*vs * aa1) / &
              (8.d0 * (1.d0-eta) * bb1 - (6.d0 + 4.d0 * eta + 5.d0 * xi) * aa1 ) )
  vsh = sqrt( xi ) * vsv
  vpv = sqrt( (15.d0 * vp*vp - 8.d0 * (1.d0-eta) * vsv*vsv ) / &
              (3.d0 + ( 8.d0 + 4.d0 * eta ) * fi_inv ) )
  vph = vpv * sqrt( fi_inv )
  rho = rho
  ! ===================================================

  ! perturbations
  if (vsv1d == 0.d0) then
    dvsv = 0.d0
  else
    dvsv = vsv / vsv1d - 1.d0
  endif
  if (vsh1d == 0.d0) then
    dvsh = 0.d0
  else
    dvsh = vsh / vsh1d - 1.d0
  endif

  dvpv = vpv / vpv1d - 1.d0
  dvph = vph / vph1d - 1.d0

  drho  = rho / rho1d - 1.d0
  eta_aniso = eta

  end subroutine model_berkeley_shsv

!
!--------------------------------------------------------------------------------------------------
!

  double precision function getdel(a0,o0,a,o)

  implicit none
  double precision,intent(in) :: a0,o0,a,o
  ! local parameters
  double precision :: q0,sq0,cq0,q,sq,cq,ff,cff,arg  ! sff

  double precision, parameter :: deg2rad = (4.d0 * datan(1.d0)) / 180.d0
  double precision, parameter :: rad2deg = 180.d0 / (4.d0 * datan(1.d0))

  q0 = (90.d0 - a0) * deg2rad
  sq0 = sin(q0)
  cq0 = cos(q0)

  q = (90.d0 - a) * deg2rad
  sq = sin(q)
  cq = cos(q)

  ff = (o - o0) * deg2rad
  !sff = sin(ff)  ! not used
  cff = cos(ff)

  arg = cq * cq0 + sq * sq0 * cff

  if (arg > 1.d0) arg = 1.d0
  if (arg < -1.d0) arg = -1.d0

  getdel = rad2deg * acos(arg)

  end function getdel

!
!--------------------------------------------------------------------------------------------------
!

  double precision function spbsp(hdel,ahdel)

  implicit none
  double precision,intent(in) :: hdel,ahdel
  ! local parameters
  double precision :: ratio,two_minus_ratio

  ! factor
  ratio = hdel / ahdel

  if (hdel < ahdel) then
    spbsp = (0.75d0 * ratio - 1.5d0) * ratio * ratio + 1.d0
  else if (hdel <= ahdel * 2.d0) then
    two_minus_ratio = 2.d0 - ratio
    if (two_minus_ratio < 0.d0) then
      spbsp = 0.d0
    else
      spbsp = 0.25d0 * two_minus_ratio * two_minus_ratio * two_minus_ratio
    endif
    !if (spbsp < 0.d0) spbsp = 0.d0    ! spbsp is negative if two_minus_ratio is negative
  else
    spbsp = 0.d0
  endif

  end function spbsp

