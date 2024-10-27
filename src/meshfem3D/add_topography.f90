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

  subroutine add_topography(xelm,yelm,zelm,ibathy_topo)

  use constants, only: myrank,NGNOD,R_UNIT_SPHERE,ONE
  use shared_parameters, only: REGIONAL_MESH_CUTOFF,REGIONAL_MESH_CUTOFF_DEPTH,USE_LOCAL_MESH,ELLIPTICITY
  use meshfem_par, only: R220,NX_BATHY,NY_BATHY,R_PLANET

  ! for old version Berkeley compatibility
  use constants, only: USE_OLD_VERSION_FORMAT,ICRUST_BERKELEY,THREE_D_MODEL_BERKELEY
  use meshfem_models_par, only: THREE_D_MODEL,REFERENCE_CRUSTAL_MODEL

  implicit none

  double precision,dimension(NGNOD), intent(inout) :: xelm,yelm,zelm

  ! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY), intent(in) :: ibathy_topo

  ! local parameters
  double precision :: r,lat,lon,elevation,rbottom
  double precision :: x,y,z
  double precision :: gamma

  integer :: ia

  ! for compatibility
  double precision :: theta,phi
  double precision :: vpvc,vphc,vsvc,vshc,etac,rhoc
  double precision :: moho
  double precision :: rmoho
  logical :: found_crust,elem_in_crust,moho_only

  ! we loop on all the points of the element
  do ia = 1,NGNOD

    x = xelm(ia)
    y = yelm(ia)
    z = zelm(ia)

    ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
    !
    ! note: at this point, the mesh is still spherical. converting from x/y/z to r/lat/lon would give a geocentric position.
    !       however, we want to add the topography of the geographic position if the mesh should become elliptical.
    !       in SPECFEM, we add an elliptical stretch factor to the radius after adding topography, which only corrects
    !       the radial direction. we thus have to correct the colatitude here for topography stretching,
    !       depending on the user defined parameter ELLIPTICITY, read in from the Par_file, to end up with an
    !       elliptical and topographic Earth.
    call xyz_2_rlatlon_dble(x,y,z,r,lat,lon,ELLIPTICITY)

    ! compute elevation at current point
    call get_topo_bathy(lat,lon,elevation,ibathy_topo)

    ! non-dimensionalize the elevation, which is in meters
    elevation = elevation / R_PLANET

    ! stretching topography between d220 and the surface
    if (REGIONAL_MESH_CUTOFF .and. USE_LOCAL_MESH) then
      rbottom = (R_PLANET - REGIONAL_MESH_CUTOFF_DEPTH*1000.d0) / R_PLANET
    else
      rbottom = R220 / R_PLANET
    endif
    gamma = (r - rbottom) / (R_UNIT_SPHERE - rbottom)

    ! old version compatility
    if (USE_OLD_VERSION_FORMAT) then
      ! Berkeley model
      if (THREE_D_MODEL == THREE_D_MODEL_BERKELEY .and. &
          REFERENCE_CRUSTAL_MODEL == ICRUST_BERKELEY) then
        ! convert lat/lon to theta/phi
        call latlon_2_thetaphi_dble(lat,lon,theta,phi)

        ! gets smoothed moho depth
        elem_in_crust = .true.
        moho_only = .true.
        call model_berkeley_crust_aniso(r,theta,phi,vpvc,vphc,vsvc,vshc,etac,rhoc,moho,found_crust,elem_in_crust,moho_only)

        rmoho = R_UNIT_SPHERE - moho

        ! if point above moho then move points, otherwise skip
        if (r <= rmoho) cycle

        ! adjust gamma stretching to moho boundary distance
        gamma = (r - rmoho) / (R_UNIT_SPHERE - rmoho)
      endif
    endif

    ! add elevation to all the points of that element
    ! also make sure gamma makes sense
    if (gamma < -0.02 .or. gamma > 1.02) then
      print '(a, 2F9.2,2F9.4)','DEBUG lat/lon/r/gamma:',lat,lon,r,gamma
      call exit_MPI(myrank,'incorrect value of gamma for topography')
    endif

    xelm(ia) = x*(ONE + gamma * elevation / r)
    yelm(ia) = y*(ONE + gamma * elevation / r)
    zelm(ia) = z*(ONE + gamma * elevation / r)

  enddo

  end subroutine add_topography

!
!-------------------------------------------------------------------------------------------------
!

  !> Hejun
  ! This subroutine uses GLL points to capture topography variation rather
  ! than using control nodes
  ! Hejun Zhu, OCT16, 2009

  subroutine add_topography_gll(xstore,ystore,zstore,ispec,nspec,ibathy_topo)

  use constants
  use shared_parameters, only: R_PLANET
  use shared_parameters, only: REGIONAL_MESH_CUTOFF,REGIONAL_MESH_CUTOFF_DEPTH,USE_LOCAL_MESH,ELLIPTICITY
  use meshfem_par, only: R220,NX_BATHY,NY_BATHY

  implicit none

  ! input parameters
  integer, intent(in) :: ispec,nspec

  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(inout):: xstore,ystore,zstore

  integer, dimension(NX_BATHY,NY_BATHY),intent(in) :: ibathy_topo

  ! local parameters used in this subroutine
  integer :: i,j,k
  double precision :: r,lat,lon,elevation,gamma,rbottom
  double precision :: x,y,z

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        x = xstore(i,j,k,ispec)
        y = ystore(i,j,k,ispec)
        z = zstore(i,j,k,ispec)

        ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
        ! note: at this point, the mesh is still spherical (no need to correct latitude for ellipticity)
        call xyz_2_rlatlon_dble(x,y,z,r,lat,lon,ELLIPTICITY)

        ! compute elevation at current point
        call get_topo_bathy(lat,lon,elevation,ibathy_topo)

        ! non-dimensionalize the elevation, which is in meters
        elevation = elevation / R_PLANET

        ! stretching topography between d220 and the surface
        if (REGIONAL_MESH_CUTOFF .and. USE_LOCAL_MESH) then
          rbottom = (R_PLANET - REGIONAL_MESH_CUTOFF_DEPTH*1000.d0) / R_PLANET
        else
          rbottom = R220 / R_PLANET
        endif
        gamma = (r - rbottom) / (R_UNIT_SPHERE - rbottom)

        ! add elevation to all the points of that element
        ! also make sure factor makes sense
        if (gamma < -0.02 .or. gamma > 1.02) then
          call exit_MPI(myrank,'incorrect value of factor for topography GLL points')
        endif

        ! since not all GLL points are exactly at R220, use a small
        ! tolerance for R220 detection
        if (abs(gamma) < SMALLVAL) then
          gamma = 0.d0
        endif

        xstore(i,j,k,ispec) = x*(ONE + gamma * elevation / r)
        ystore(i,j,k,ispec) = y*(ONE + gamma * elevation / r)
        zstore(i,j,k,ispec) = z*(ONE + gamma * elevation / r)
      enddo
    enddo
  enddo

  end subroutine add_topography_gll
