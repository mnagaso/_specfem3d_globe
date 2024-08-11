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

  subroutine write_movie_volume_mesh_hdf5(nu_3dmovie,num_ibool_3dmovie,mask_3dmovie,mask_ibool_3dmovie, &
                                          muvstore_crust_mantle_3dmovie,npoints_3dmovie)

  use specfem_par
  use specfem_par_crustmantle, only: ibool_crust_mantle,rstore_crust_mantle

#ifdef USE_HDF5
  use specfem_par_movie_hdf5
#endif

  implicit none

  integer,intent(in) :: npoints_3dmovie
  integer, dimension(NGLOB_CRUST_MANTLE_3DMOVIE),intent(in) :: num_ibool_3dmovie

  real(kind=CUSTOM_REAL), dimension(3,3,npoints_3dmovie),intent(inout) :: nu_3dmovie

  logical, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE),intent(in) :: mask_3dmovie
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE),intent(in) :: muvstore_crust_mantle_3dmovie
  logical, dimension(NGLOB_CRUST_MANTLE_3DMOVIE),intent(in) :: mask_ibool_3dmovie

#ifdef USE_HDF5
  ! local parameters
  integer :: ipoints_3dmovie,ispec,i,j,k,iNIT
  integer :: iglob,iglob_center
  real(kind=CUSTOM_REAL) :: rval,thetaval,phival,xval,yval,zval,st,ct,sp,cp
  real(kind=CUSTOM_REAL), dimension(npoints_3dmovie) :: store_val3D_x,store_val3D_y, store_val3D_z
  real(kind=CUSTOM_REAL), dimension(npoints_3dmovie) :: store_val3D_mu
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_poin_vol
  integer :: npoints_vol_mov_all_proc
  integer, dimension(9,npoints_3dmovie) :: elm_conn

  ! safety check
  if (NDIM /= 3) stop 'movie volume output requires NDIM = 3'

  ! prepare offset array
  call gather_all_all_singlei(npoints_3dmovie, offset_poin_vol, NPROCTOT_VAL)
  npoints_vol_mov_all_proc = sum(offset_poin_vol)

  ! stepping
  if (MOVIE_COARSE) then
    iNIT = NGLLX-1
  else
    iNIT = 1
  endif

  ! loops over all elements
  ipoints_3dmovie = 0
  do ispec = 1,NSPEC_CRUST_MANTLE

    ! checks center of element for movie flag
    iglob_center = ibool_crust_mantle((NGLLX+1)/2,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)

    ! checks if movie element
    if (mask_ibool_3dmovie(iglob_center)) then

      ! stores element coordinates
      do k = 1,NGLLZ,iNIT
        do j = 1,NGLLY,iNIT
          do i = 1,NGLLX,iNIT
            ! only store points once
            if (mask_3dmovie(i,j,k,ispec)) then
              ! point increment
              ipoints_3dmovie = ipoints_3dmovie + 1

              ! gets point position
              iglob = ibool_crust_mantle(i,j,k,ispec)

              rval = rstore_crust_mantle(1,iglob)
              thetaval = rstore_crust_mantle(2,iglob)
              phival = rstore_crust_mantle(3,iglob)

              !x,y,z store have been converted to r theta phi already, need to revert back for xyz output
              call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)

              store_val3D_x(ipoints_3dmovie) = xval
              store_val3D_y(ipoints_3dmovie) = yval
              store_val3D_z(ipoints_3dmovie) = zval
              store_val3D_mu(ipoints_3dmovie) = muvstore_crust_mantle_3dmovie(i,j,k,ispec)

              st = sin(thetaval)
              ct = cos(thetaval)
              sp = sin(phival)
              cp = cos(phival)

              nu_3dmovie(1,1,ipoints_3dmovie) = -ct*cp
              nu_3dmovie(1,2,ipoints_3dmovie) = -ct*sp
              nu_3dmovie(1,3,ipoints_3dmovie) = st
              nu_3dmovie(2,1,ipoints_3dmovie) = -sp
              nu_3dmovie(2,2,ipoints_3dmovie) = cp
              nu_3dmovie(2,3,ipoints_3dmovie) = 0.d0
              nu_3dmovie(3,1,ipoints_3dmovie) = st*cp
              nu_3dmovie(3,2,ipoints_3dmovie) = st*sp
              nu_3dmovie(3,3,ipoints_3dmovie) = ct
            endif !mask_3dmovie
          enddo  !i
        enddo  !j
      enddo  !k
    endif

  enddo !ispec

  ! create elm_conn for movie
  call get_conn_for_movie(elm_conn, sum(offset_poin_vol(0:myrank-1)), iNIT, npoints_3dmovie, num_ibool_3dmovie, mask_ibool_3dmovie)

  ! TODO ADD IOSERVER

  ! initialize h5 file for volume movie
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  ! file name and group name
  file_name = 'movie_volume.h5'
  group_name = 'mesh'

  ! create the file, group and dataset
  if (myrank == 0) then
    call h5_create_file(file_name)
    call h5_open_or_create_group(group_name)

    ! create the dataset
    call h5_create_dataset_gen_in_group('elm_conn', (/9, npoints_vol_mov_all_proc/), 2, 1)
    call h5_create_dataset_gen_in_group('x', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group('y', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group('z', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)

    ! close the group
    call h5_close_group()
    ! close the file
    call h5_close_file()

  endif

  ! write the data
  call h5_open_file_p_collect(file_name)
  call h5_open_group(group_name)

  ! write the data
  call h5_write_dataset_collect_hyperslab_in_group('elm_conn', elm_conn, (/0, sum(offset_poin_vol(0:myrank-1))/), if_col)
  call h5_write_dataset_collect_hyperslab_in_group('x', store_val3D_x, (/sum(offset_poin_vol(0:myrank-1))/), if_col)
  call h5_write_dataset_collect_hyperslab_in_group('y', store_val3D_y, (/sum(offset_poin_vol(0:myrank-1))/), if_col)
  call h5_write_dataset_collect_hyperslab_in_group('z', store_val3D_z, (/sum(offset_poin_vol(0:myrank-1))/), if_col)

  ! close the group
  call h5_close_group()
  ! close the file
  call h5_close_file_p()

  ! write xdmf for all timesteps
  call write_xdmf_vol_hdf5(.true., .false., .false., .false., .false., .false.)

#else

    print*, 'Error: HDF5 is not enabled in this version of the code.'
    print*, 'Please recompile with the HDF5 option enabled with --with-hdf5'
    stop

#endif

  end subroutine write_movie_volume_mesh_hdf5


  subroutine write_movie_volume_strains_hdf5(vnspec_eps_cm, &
                                        eps_trace_over_3_crust_mantle, &
                                        vnspec_cm, &
                                        epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                                        epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle)

  use constants_solver
  use shared_parameters, only: OUTPUT_FILES,MOVIE_VOLUME_TYPE,MOVIE_COARSE
  use specfem_par, only: it
  use specfem_par_movie, only: npoints_3dmovie,muvstore_crust_mantle_3dmovie,mask_3dmovie,nu_3dmovie

#ifdef USE_HDF5
  use specfem_par_movie_hdf5
#endif

  implicit none

  ! input
  integer,intent(in) :: vnspec_eps_cm
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec_eps_cm),intent(in) :: eps_trace_over_3_crust_mantle

  integer,intent(in) :: vnspec_cm
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec_cm),intent(in) :: &
    epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
    epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle

#ifdef USE_HDF5
  ! variables
  real(kind=CUSTOM_REAL) :: muv_3dmovie
  real(kind=CUSTOM_REAL),dimension(3,3) :: eps_loc,eps_loc_new
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: store_val3d_NN,store_val3d_EE,store_val3d_ZZ, &
                                                     store_val3d_NE,store_val3d_NZ,store_val3d_EZ
  integer :: ipoints_3dmovie,i,j,k,ispec,iNIT,ier
  character(len=1) :: movie_prefix
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_poin_vol
  integer :: npoints_vol_mov_all_proc

  ! prepare offset array
  call gather_all_all_singlei(npoints_3dmovie, offset_poin_vol, NPROCTOT_VAL)
  npoints_vol_mov_all_proc = sum(offset_poin_vol)

  ! check
  if (NDIM /= 3) call exit_MPI(myrank, 'write_movie_volume_strains() requires NDIM = 3')
  if (vnspec_cm /= NSPEC_CRUST_MANTLE) call exit_MPI(myrank,'Invalid vnspec_cm value for write_movie_volume_strains() routine')

  ! allocates arrays
  allocate(store_val3d_NN(npoints_3dmovie), &
           store_val3d_EE(npoints_3dmovie), &
           store_val3d_ZZ(npoints_3dmovie), &
           store_val3d_NE(npoints_3dmovie), &
           store_val3d_NZ(npoints_3dmovie), &
           store_val3d_EZ(npoints_3dmovie), &
           stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating store_val3d_ .. arrays')

  if (MOVIE_VOLUME_TYPE == 1) then
    movie_prefix='E' ! strain
  else if (MOVIE_VOLUME_TYPE == 2) then
    movie_prefix='S' ! time integral of strain
  else if (MOVIE_VOLUME_TYPE == 3) then
    movie_prefix='P' ! potency, or integral of strain x \mu
  endif

  ! stepping
  if (MOVIE_COARSE) then
   iNIT = NGLLX-1
  else
   iNIT = 1
  endif

  ipoints_3dmovie = 0
  do ispec = 1,NSPEC_CRUST_MANTLE
    do k = 1,NGLLZ,iNIT
      do j = 1,NGLLY,iNIT
        do i = 1,NGLLX,iNIT
          if (mask_3dmovie(i,j,k,ispec)) then
            ipoints_3dmovie = ipoints_3dmovie + 1
            muv_3dmovie = muvstore_crust_mantle_3dmovie(i,j,k,ispec)

            eps_loc(1,1) = eps_trace_over_3_crust_mantle(i,j,k,ispec) + epsilondev_xx_crust_mantle(i,j,k,ispec)
            eps_loc(2,2) = eps_trace_over_3_crust_mantle(i,j,k,ispec) + epsilondev_yy_crust_mantle(i,j,k,ispec)
            eps_loc(3,3) = eps_trace_over_3_crust_mantle(i,j,k,ispec) &
                           - epsilondev_xx_crust_mantle(i,j,k,ispec) &
                           - epsilondev_yy_crust_mantle(i,j,k,ispec)

            eps_loc(1,2) = epsilondev_xy_crust_mantle(i,j,k,ispec)
            eps_loc(1,3) = epsilondev_xz_crust_mantle(i,j,k,ispec)
            eps_loc(2,3) = epsilondev_yz_crust_mantle(i,j,k,ispec)

            eps_loc(2,1) = eps_loc(1,2)
            eps_loc(3,1) = eps_loc(1,3)
            eps_loc(3,2) = eps_loc(2,3)

            ! rotate eps_loc to spherical coordinates
            eps_loc_new(:,:) = matmul(matmul(nu_3dmovie(:,:,ipoints_3dmovie),eps_loc(:,:)), &
                                      transpose(nu_3dmovie(:,:,ipoints_3dmovie)))
            if (MOVIE_VOLUME_TYPE == 3) eps_loc_new(:,:) = eps_loc(:,:)*muv_3dmovie

            store_val3d_NN(ipoints_3dmovie) = eps_loc_new(1,1)
            store_val3d_EE(ipoints_3dmovie) = eps_loc_new(2,2)
            store_val3d_ZZ(ipoints_3dmovie) = eps_loc_new(3,3)
            store_val3d_NE(ipoints_3dmovie) = eps_loc_new(1,2)
            store_val3d_NZ(ipoints_3dmovie) = eps_loc_new(1,3)
            store_val3d_EZ(ipoints_3dmovie) = eps_loc_new(2,3)
          endif
        enddo
      enddo
    enddo
  enddo
  if (ipoints_3dmovie /= npoints_3dmovie) stop 'did not find the right number of points for 3D movie'

  ! initialize h5 file for volume movie
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  ! create group and datasets
  file_name = trim(OUTPUT_FILES) // '/movie_volume.h5'
  group_name = 'it_' // trim(i2c(it))

  if (myrank == 0) then

    ! create the file, group and dataset
    call h5_open_file(file_name)
    call h5_open_or_create_group(group_name)

    ! create the dataset
    call h5_create_dataset_gen_in_group(trim(movie_prefix)//'NN', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group(trim(movie_prefix)//'EE', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group(trim(movie_prefix)//'ZZ', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group(trim(movie_prefix)//'NE', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group(trim(movie_prefix)//'NZ', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group(trim(movie_prefix)//'EZ', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)

    ! close the group
    call h5_close_group()
    ! close the file
    call h5_close_file()

  endif

  ! write the data
  call h5_open_file_p_collect(file_name)
  call h5_open_group(group_name)

  call h5_write_dataset_collect_hyperslab_in_group(trim(movie_prefix)//'NN', store_val3d_NN, (/sum(offset_poin_vol(0:myrank-1))/), if_col)
  call h5_write_dataset_collect_hyperslab_in_group(trim(movie_prefix)//'EE', store_val3d_EE, (/sum(offset_poin_vol(0:myrank-1))/), if_col)
  call h5_write_dataset_collect_hyperslab_in_group(trim(movie_prefix)//'ZZ', store_val3d_ZZ, (/sum(offset_poin_vol(0:myrank-1))/), if_col)
  call h5_write_dataset_collect_hyperslab_in_group(trim(movie_prefix)//'NE', store_val3d_NE, (/sum(offset_poin_vol(0:myrank-1))/), if_col)
  call h5_write_dataset_collect_hyperslab_in_group(trim(movie_prefix)//'NZ', store_val3d_NZ, (/sum(offset_poin_vol(0:myrank-1))/), if_col)
  call h5_write_dataset_collect_hyperslab_in_group(trim(movie_prefix)//'EZ', store_val3d_EZ, (/sum(offset_poin_vol(0:myrank-1))/), if_col)

  call h5_close_group()
  call h5_close_file_p()

  deallocate(store_val3d_NN,store_val3d_EE,store_val3d_ZZ, &
             store_val3d_NE,store_val3d_NZ,store_val3d_EZ)


#else

    print*, 'Error: HDF5 is not enabled in this version of the code.'
    print*, 'Please recompile with the HDF5 option enabled with --with-hdf5'
    stop

#endif

  end subroutine write_movie_volume_strains_hdf5




subroutine write_movie_volume_divcurl_hdf5(vnspec_eps_cm,eps_trace_over_3_crust_mantle, &
                                       div_displ_outer_core, &
                                       accel_outer_core,kappavstore_outer_core,rhostore_outer_core,ibool_outer_core, &
                                       vnspec_eps_ic,eps_trace_over_3_inner_core, &
                                       vnspec_cm,epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                                       epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
                                       vnspec_ic,epsilondev_xx_inner_core,epsilondev_yy_inner_core,epsilondev_xy_inner_core, &
                                       epsilondev_xz_inner_core,epsilondev_yz_inner_core)

! outputs divergence and curl: MOVIE_VOLUME_TYPE == 4

  use constants_solver
  use shared_parameters, only: OUTPUT_FILES
  use specfem_par, only: it

#ifdef USE_HDF5
  use specfem_par_movie_hdf5
#endif

  implicit none

  integer,intent(in) :: vnspec_eps_cm
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec_eps_cm) :: eps_trace_over_3_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_3DMOVIE) :: div_displ_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: accel_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: rhostore_outer_core,kappavstore_outer_core
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

  integer,intent(in) :: vnspec_eps_ic
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec_eps_ic) :: eps_trace_over_3_inner_core

  integer,intent(in) :: vnspec_cm
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec_cm) :: &
    epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
    epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle

  integer,intent(in) :: vnspec_ic
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec_ic) :: &
    epsilondev_xx_inner_core,epsilondev_yy_inner_core,epsilondev_xy_inner_core, &
    epsilondev_xz_inner_core,epsilondev_yz_inner_core

#ifdef USE_HDF5

  ! local parameters
  real(kind=CUSTOM_REAL) :: rhol,kappal
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: div_s_outer_core
  integer :: ispec,iglob,i,j,k,ier
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: tmp_data

  ! output parameters
  logical,parameter :: MOVIE_OUTPUT_DIV = .true.          ! divergence
  logical,parameter :: MOVIE_OUTPUT_CURL = .true.         ! curl
  logical,parameter :: MOVIE_OUTPUT_CURLNORM = .true.     ! Frobenius norm of curl

  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_vol_eps_cm
  integer :: nspec_vol_mov_all_proc_eps_cm
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_vol_eps_ic
  integer :: nspec_vol_mov_all_proc_eps_ic
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_vol_cm
  integer :: nspec_vol_mov_all_proc_cm
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_vol_oc
  integer :: nspec_vol_mov_all_proc_oc
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_vol_ic
  integer :: nspec_vol_mov_all_proc_ic
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_vol_cm_norm
  integer :: nspec_vol_mov_all_proc_cm_norm
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_vol_ic_norm
  integer :: nspec_vol_mov_all_proc_ic_norm

  ! prepare offset array
  call gather_all_all_singlei(vnspec_eps_cm, offset_nspec_vol_eps_cm, NPROCTOT_VAL)
  nspec_vol_mov_all_proc_eps_cm = sum(offset_nspec_vol_cm)
  call gather_all_all_singlei(vnspec_eps_ic, offset_nspec_vol_eps_ic, NPROCTOT_VAL)
  nspec_vol_mov_all_proc_eps_ic = sum(offset_nspec_vol_ic)
  call gather_all_all_singlei(vnspec_cm, offset_nspec_vol_cm, NPROCTOT_VAL)
  nspec_vol_mov_all_proc_cm = sum(offset_nspec_vol_cm)
  call gather_all_all_singlei(vnspec_ic, offset_nspec_vol_ic, NPROCTOT_VAL)
  nspec_vol_mov_all_proc_ic = sum(offset_nspec_vol_ic)
  if (NSPEC_OUTER_CORE_3DMOVIE > 1) then
    call gather_all_all_singlei(NSPEC_OUTER_CORE_3DMOVIE, offset_nspec_vol_oc, NPROCTOT_VAL)
  else
    call gather_all_all_singlei(NSPEC_OUTER_CORE, offset_nspec_vol_oc, NPROCTOT_VAL)
  endif
  nspec_vol_mov_all_proc_oc = sum(offset_nspec_vol_oc)
  if (MOVIE_OUTPUT_CURLNORM) then
    call gather_all_all_singlei(NSPEC_CRUST_MANTLE, offset_nspec_vol_cm_norm, NPROCTOT_VAL)
    nspec_vol_mov_all_proc_cm_norm = sum(offset_nspec_vol_cm_norm)
    call gather_all_all_singlei(NSPEC_INNER_CORE, offset_nspec_vol_ic_norm, NPROCTOT_VAL)
    nspec_vol_mov_all_proc_ic_norm = sum(offset_nspec_vol_ic_norm)
  endif

  ! checks
  if (vnspec_cm /= NSPEC_CRUST_MANTLE) call exit_MPI(myrank,'Invalid vnspec_cm value for write_movie_volume_divcurl() routine')

  ! create group and datasets
  file_name = trim(OUTPUT_FILES) // '/movie_volume.h5'
  group_name = 'it_' // trim(i2c(it))

  if (myrank == 0) then
    call h5_open_file(file_name)
    call h5_open_or_create_group(group_name)

    if (MOVIE_OUTPUT_DIV) then
      call h5_create_dataset_gen_in_group('reg1_div_displ', (/NGLLX,NGLLY,NGLLZ,nspec_vol_mov_all_proc_eps_cm/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('reg2_div_displ', (/NGLLX,NGLLY,NGLLZ,nspec_vol_mov_all_proc_oc/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('reg3_div_displ', (/NGLLX,NGLLY,NGLLZ,nspec_vol_mov_all_proc_eps_ic/), 4, CUSTOM_REAL)
    endif

    if (MOVIE_OUTPUT_CURL) then
      call h5_create_dataset_gen_in_group('crust_mantle_epsdev_displ_xx', (/NGLLX,NGLLY,NGLLZ,nspec_vol_mov_all_proc_cm/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('crust_mantle_epsdev_displ_yy', (/NGLLX,NGLLY,NGLLZ,nspec_vol_mov_all_proc_cm/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('crust_mantle_epsdev_displ_xy', (/NGLLX,NGLLY,NGLLZ,nspec_vol_mov_all_proc_cm/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('crust_mantle_epsdev_displ_xz', (/NGLLX,NGLLY,NGLLZ,nspec_vol_mov_all_proc_cm/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('crust_mantle_epsdev_displ_yz', (/NGLLX,NGLLY,NGLLZ,nspec_vol_mov_all_proc_cm/), 4, CUSTOM_REAL)

      call h5_create_dataset_gen_in_group('inner_core_epsdev_displ_xx', (/NGLLX,NGLLY,NGLLZ,nspec_vol_mov_all_proc_ic/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('inner_core_epsdev_displ_yy', (/NGLLX,NGLLY,NGLLZ,nspec_vol_mov_all_proc_ic/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('inner_core_epsdev_displ_xy', (/NGLLX,NGLLY,NGLLZ,nspec_vol_mov_all_proc_ic/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('inner_core_epsdev_displ_xz', (/NGLLX,NGLLY,NGLLZ,nspec_vol_mov_all_proc_ic/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('inner_core_epsdev_displ_yz', (/NGLLX,NGLLY,NGLLZ,nspec_vol_mov_all_proc_ic/), 4, CUSTOM_REAL)
    endif

    if (MOVIE_OUTPUT_CURLNORM) then
      call h5_create_dataset_gen_in_group('reg1_epsdev_displ_norm', (/NGLLX,NGLLY,NGLLZ,nspec_vol_mov_all_proc_cm_norm/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group('reg3_epsdev_displ_norm', (/NGLLX,NGLLY,NGLLZ,nspec_vol_mov_all_proc_ic_norm/), 4, CUSTOM_REAL)
    endif

    call h5_close_group()
    call h5_close_file()
  endif

  call synchronize_all()

  ! write the data
  call h5_open_file_p_collect(file_name)
  call h5_open_group(group_name)

  if (MOVIE_OUTPUT_DIV) then
    call h5_write_dataset_collect_hyperslab_in_group('reg1_div_displ', eps_trace_over_3_crust_mantle, (/0,0,0,offset_nspec_vol_eps_cm(0:myrank-1)/), if_col)

    if (NSPEC_OUTER_CORE_3DMOVIE > 1) then
      call h5_write_dataset_collect_hyperslab_in_group('reg2_div_displ', ONE_THIRD*div_displ_outer_core, (/0,0,0,offset_nspec_vol_oc(0:myrank-1)/), if_col)
    else
      allocate(div_s_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array div_s_outer_core')
      do ispec = 1, NSPEC_OUTER_CORE
        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i = 1, NGLLX
              iglob = ibool_outer_core(i,j,k,ispec)
              rhol = rhostore_outer_core(i,j,k,ispec)
              kappal = kappavstore_outer_core(i,j,k,ispec)
              div_s_outer_core(i,j,k,ispec) = rhol * accel_outer_core(iglob) / kappal
            enddo
          enddo
        enddo
      enddo
      call h5_write_dataset_collect_hyperslab_in_group('reg2_div_displ', div_s_outer_core, (/0,0,0,offset_nspec_vol_oc(0:myrank-1)/), if_col)
      deallocate(div_s_outer_core)
    endif

    call h5_write_dataset_collect_hyperslab_in_group('reg3_div_displ', eps_trace_over_3_inner_core, (/0,0,0,offset_nspec_vol_eps_ic(0:myrank-1)/), if_col)

  endif

  if (MOVIE_OUTPUT_CURL) then
    call h5_write_dataset_collect_hyperslab_in_group('crust_mantle_epsdev_displ_xx', epsilondev_xx_crust_mantle, (/0,0,0,offset_nspec_vol_cm(0:myrank-1)/), if_col)
    call h5_write_dataset_collect_hyperslab_in_group('crust_mantle_epsdev_displ_yy', epsilondev_yy_crust_mantle, (/0,0,0,offset_nspec_vol_cm(0:myrank-1)/), if_col)
    call h5_write_dataset_collect_hyperslab_in_group('crust_mantle_epsdev_displ_xy', epsilondev_xy_crust_mantle, (/0,0,0,offset_nspec_vol_cm(0:myrank-1)/), if_col)
    call h5_write_dataset_collect_hyperslab_in_group('crust_mantle_epsdev_displ_xz', epsilondev_xz_crust_mantle, (/0,0,0,offset_nspec_vol_cm(0:myrank-1)/), if_col)
    call h5_write_dataset_collect_hyperslab_in_group('crust_mantle_epsdev_displ_yz', epsilondev_yz_crust_mantle, (/0,0,0,offset_nspec_vol_cm(0:myrank-1)/), if_col)

    call h5_write_dataset_collect_hyperslab_in_group('inner_core_epsdev_displ_xx', epsilondev_xx_inner_core, (/0,0,0,offset_nspec_vol_ic(0:myrank-1)/), if_col)
    call h5_write_dataset_collect_hyperslab_in_group('inner_core_epsdev_displ_yy', epsilondev_yy_inner_core, (/0,0,0,offset_nspec_vol_ic(0:myrank-1)/), if_col)
    call h5_write_dataset_collect_hyperslab_in_group('inner_core_epsdev_displ_xy', epsilondev_xy_inner_core, (/0,0,0,offset_nspec_vol_ic(0:myrank-1)/), if_col)
    call h5_write_dataset_collect_hyperslab_in_group('inner_core_epsdev_displ_xz', epsilondev_xz_inner_core, (/0,0,0,offset_nspec_vol_ic(0:myrank-1)/), if_col)
    call h5_write_dataset_collect_hyperslab_in_group('inner_core_epsdev_displ_yz', epsilondev_yz_inner_core, (/0,0,0,offset_nspec_vol_ic(0:myrank-1)/), if_col)
  endif

  if (MOVIE_OUTPUT_CURLNORM) then
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')
    ! Frobenius norm
    do ispec = 1, NSPEC_CRUST_MANTLE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            tmp_data(i,j,k,ispec) = sqrt( epsilondev_xx_crust_mantle(i,j,k,ispec)**2 &
                                          + epsilondev_yy_crust_mantle(i,j,k,ispec)**2 &
                                          + epsilondev_xy_crust_mantle(i,j,k,ispec)**2 &
                                          + epsilondev_xz_crust_mantle(i,j,k,ispec)**2 &
                                          + epsilondev_yz_crust_mantle(i,j,k,ispec)**2)
          enddo
        enddo
      enddo
    enddo
    call h5_write_dataset_collect_hyperslab_in_group('reg1_epsdev_displ_norm', tmp_data, (/0,0,0,offset_nspec_vol_cm_norm(0:myrank-1)/), if_col)
    deallocate(tmp_data)

    ! Frobenius norm
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')
    do ispec = 1, NSPEC_INNER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            tmp_data(i,j,k,ispec) = sqrt( epsilondev_xx_inner_core(i,j,k,ispec)**2 &
                                          + epsilondev_yy_inner_core(i,j,k,ispec)**2 &
                                          + epsilondev_xy_inner_core(i,j,k,ispec)**2 &
                                          + epsilondev_xz_inner_core(i,j,k,ispec)**2 &
                                          + epsilondev_yz_inner_core(i,j,k,ispec)**2)
          enddo
        enddo
      enddo
    enddo
    call h5_write_dataset_collect_hyperslab_in_group('reg3_epsdev_displ_norm', tmp_data, (/0,0,0,offset_nspec_vol_ic_norm(0:myrank-1)/), if_col)
    deallocate(tmp_data)
  endif

  call h5_close_group()
  call h5_close_file()

  !if (myrank == 0) call write_xdmf...

#else

    print*, 'Error: HDF5 is not enabled in this version of the code.'
    print*, 'Please recompile with the HDF5 option enabled with --with-hdf5'
    stop

#endif

  end subroutine write_movie_volume_divcurl_hdf5



  subroutine write_movie_volume_vector_hdf5(npoints_3dmovie, &
                                       ibool_crust_mantle,vector_crust_mantle, &
                                       scalingval,mask_3dmovie,nu_3dmovie)

! outputs displacement/velocity: MOVIE_VOLUME_TYPE == 5 / 6

  use constants_solver
  use shared_parameters, only: OUTPUT_FILES,MOVIE_VOLUME_TYPE,MOVIE_COARSE
  use specfem_par, only: it

#ifdef USE_HDF5
  use specfem_par_movie_hdf5
#endif

  implicit none

  ! input
  integer :: npoints_3dmovie
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  ! displacement or velocity array
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: vector_crust_mantle
  logical, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE) :: mask_3dmovie

  double precision :: scalingval
  real(kind=CUSTOM_REAL), dimension(NDIM,NDIM,npoints_3dmovie) :: nu_3dmovie

#ifdef USE_HDF5

  ! local variables
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val3d_N,store_val3d_E,store_val3d_Z
  real(kind=CUSTOM_REAL), dimension(NDIM) :: vector_local,vector_local_new

  integer :: ipoints_3dmovie,i,j,k,ispec,iNIT,iglob,ier
  character(len=2) :: movie_prefix

  integer, dimension(0:NPROCTOT_VAL-1) :: offset_npoints_3dmovie
  integer :: npoints_3dmovie_all_proc


  ! check
  if (NDIM /= 3) call exit_MPI(myrank,'write_movie_volume requires NDIM = 3')

  ! gather number of points
  call gather_all_all_singlei(npoints_3dmovie, offset_npoints_3dmovie, NPROCTOT_VAL)
  npoints_3dmovie_all_proc = sum(offset_npoints_3dmovie)

  ! initialize h5 file for volume movie
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = trim(OUTPUT_FILES) // '/movie_volume.h5'
  group_name = 'it_' // trim(i2c(it))

  ! allocates arrays
  allocate(store_val3d_N(npoints_3dmovie), &
           store_val3d_E(npoints_3dmovie), &
           store_val3d_Z(npoints_3dmovie), &
           stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating store_val3d_N,.. movie arrays')

  if (MOVIE_VOLUME_TYPE == 5) then
    movie_prefix='DI' ! displacement
  else if (MOVIE_VOLUME_TYPE == 6) then
    movie_prefix='VE' ! velocity
  endif

  if (MOVIE_COARSE) then
   iNIT = NGLLX-1
  else
   iNIT = 1
  endif

  ipoints_3dmovie = 0

  ! stores field in crust/mantle region
  do ispec = 1,NSPEC_CRUST_MANTLE
    do k = 1,NGLLZ,iNIT
      do j = 1,NGLLY,iNIT
        do i = 1,NGLLX,iNIT
          if (mask_3dmovie(i,j,k,ispec)) then
            ipoints_3dmovie = ipoints_3dmovie + 1
            iglob = ibool_crust_mantle(i,j,k,ispec)

            ! dimensionalizes field by scaling
            vector_local(:) = vector_crust_mantle(:,iglob)*real(scalingval,kind=CUSTOM_REAL)

            ! rotate eps_loc to spherical coordinates
            vector_local_new(:) = matmul(nu_3dmovie(:,:,ipoints_3dmovie), vector_local(:))

            ! stores field
            store_val3d_N(ipoints_3dmovie) = vector_local_new(1)
            store_val3d_E(ipoints_3dmovie) = vector_local_new(2)
            store_val3d_Z(ipoints_3dmovie) = vector_local_new(3)
          endif
        enddo
      enddo
   enddo
  enddo
  close(IOUT)

  ! checks number of processed points
  if (ipoints_3dmovie /= npoints_3dmovie) stop 'did not find the right number of points for 3D movie'

  ! create group and datasets
  if (myrank == 0) then
    call h5_open_file(file_name)
    call h5_open_or_create_group(group_name)

    call h5_create_dataset_gen_in_group(trim(movie_prefix)//'N', (/npoints_3dmovie_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group(trim(movie_prefix)//'E', (/npoints_3dmovie_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group(trim(movie_prefix)//'Z', (/npoints_3dmovie_all_proc/), 1, CUSTOM_REAL)

    call h5_close_group()
    call h5_close_file()
  endif

  call synchronize_all()

  ! write the data
  call h5_open_file_p_collect(file_name)
  call h5_open_group(group_name)

  call h5_write_dataset_collect_hyperslab_in_group(trim(movie_prefix)//'N', store_val3d_N(1:npoints_3dmovie), (/offset_npoints_3dmovie(0:myrank-1)/), if_col)
  call h5_write_dataset_collect_hyperslab_in_group(trim(movie_prefix)//'E', store_val3d_E(1:npoints_3dmovie), (/offset_npoints_3dmovie(0:myrank-1)/), if_col)
  call h5_write_dataset_collect_hyperslab_in_group(trim(movie_prefix)//'Z', store_val3d_Z(1:npoints_3dmovie), (/offset_npoints_3dmovie(0:myrank-1)/), if_col)

  call h5_close_group()
  call h5_close_file()

  deallocate(store_val3d_N,store_val3d_E,store_val3d_Z)

#else

    print*, 'Error: HDF5 is not enabled in this version of the code.'
    print*, 'Please recompile with the HDF5 option enabled with --with-hdf5'
    stop

#endif


  end subroutine write_movie_volume_vector_hdf5


  subroutine write_movie_volume_displnorm_hdf5(displ_crust_mantle,displ_inner_core,displ_outer_core, &
                                         ibool_crust_mantle,ibool_inner_core,ibool_outer_core)

! outputs norm of displacement: MOVIE_VOLUME_TYPE == 7

  use constants_solver
  use shared_parameters, only: OUTPUT_FILES
  use specfem_par, only: it, scale_displ

#ifdef USE_HDF5
  use specfem_par_movie_hdf5
#endif

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: displ_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: displ_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: displ_outer_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

#ifdef USE_HDF5

  ! local parameters
  integer :: ispec,iglob,i,j,k,ier
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: tmp_data

  logical,parameter :: OUTPUT_CRUST_MANTLE = .true.
  logical,parameter :: OUTPUT_OUTER_CORE = .true.
  logical,parameter :: OUTPUT_INNER_CORE = .true.

  ! offset arrays
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_cm
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_oc
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_ic
  integer :: nspec_cm_all_proc,nspec_oc_all_proc,nspec_ic_all_proc

  ! gather number of points
  call gather_all_all_singlei(NSPEC_CRUST_MANTLE, offset_nspec_cm, NPROCTOT_VAL)
  call gather_all_all_singlei(NSPEC_OUTER_CORE, offset_nspec_oc, NPROCTOT_VAL)
  call gather_all_all_singlei(NSPEC_INNER_CORE, offset_nspec_ic, NPROCTOT_VAL)
  nspec_cm_all_proc = sum(offset_nspec_cm)
  nspec_oc_all_proc = sum(offset_nspec_oc)
  nspec_ic_all_proc = sum(offset_nspec_ic)

  ! initialize h5 file for volume movie
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = trim(OUTPUT_FILES) // '/movie_volume.h5'
  group_name = 'it_' // trim(i2c(it))

  ! create group and datasets
  if (myrank == 0) then
    call h5_open_file(file_name)
    call h5_open_or_create_group(group_name)

    if (OUTPUT_CRUST_MANTLE) then
      call h5_create_dataset_gen_in_group('reg1_displ', (/NGLLX,NGLLY,NGLLZ,nspec_cm_all_proc/), 4, CUSTOM_REAL)
    endif
    if (OUTPUT_OUTER_CORE) then
      call h5_create_dataset_gen_in_group('reg2_displ', (/NGLLX,NGLLY,NGLLZ,nspec_oc_all_proc/), 4, CUSTOM_REAL)
    endif
    if (OUTPUT_INNER_CORE) then
      call h5_create_dataset_gen_in_group('reg3_displ', (/NGLLX,NGLLY,NGLLZ,nspec_ic_all_proc/), 4, CUSTOM_REAL)
    endif

    call h5_close_group()
    call h5_close_file()
  endif

  call synchronize_all()

  ! write the data
  call h5_open_file_p_collect(file_name)
  call h5_open_group(group_name)

  ! outputs norm of displacement
  if (OUTPUT_CRUST_MANTLE) then
    ! crust mantle
    ! these binary arrays can be converted into mesh format using the utility ./bin/xcombine_vol_data
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_CRUST_MANTLE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_crust_mantle(i,j,k,ispec)
            ! norm
            tmp_data(i,j,k,ispec) = real(scale_displ,kind=CUSTOM_REAL) * sqrt( displ_crust_mantle(1,iglob)**2 &
                                          + displ_crust_mantle(2,iglob)**2 &
                                          + displ_crust_mantle(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo

    call h5_write_dataset_collect_hyperslab_in_group('reg1_displ', tmp_data, (/0,0,0,offset_nspec_cm(0:myrank-1)/), if_col)

    deallocate(tmp_data)
  endif

  if (OUTPUT_OUTER_CORE) then
    ! outer core
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_OUTER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_outer_core(i,j,k,ispec)
            ! norm
            ! note: disp_outer_core is potential, this just outputs the potential,
            !          not the actual displacement u = grad(rho * Chi) / rho
            tmp_data(i,j,k,ispec) = abs(displ_outer_core(iglob))
          enddo
        enddo
      enddo
    enddo

    call h5_write_dataset_collect_hyperslab_in_group('reg2_displ', tmp_data, (/0,0,0,offset_nspec_oc(0:myrank-1)/), if_col)

    deallocate(tmp_data)
  endif

  if (OUTPUT_INNER_CORE) then
    ! inner core
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_INNER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_inner_core(i,j,k,ispec)
            ! norm
            tmp_data(i,j,k,ispec) = real(scale_displ,kind=CUSTOM_REAL) * sqrt( displ_inner_core(1,iglob)**2 &
                                          + displ_inner_core(2,iglob)**2 &
                                          + displ_inner_core(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo

    call h5_write_dataset_collect_hyperslab_in_group('reg3_displ', tmp_data, (/0,0,0,offset_nspec_ic(0:myrank-1)/), if_col)

    deallocate(tmp_data)
  endif

  call h5_close_group()
  call h5_close_file_p()


#else

    print*, 'Error: HDF5 is not enabled in this version of the code.'
    print*, 'Please recompile with the HDF5 option enabled with --with-hdf5'
    stop

#endif

  end subroutine write_movie_volume_displnorm_hdf5


subroutine write_movie_volume_velnorm_hdf5(veloc_crust_mantle,veloc_inner_core,veloc_outer_core, &
                                       ibool_crust_mantle,ibool_inner_core,ibool_outer_core)

! outputs norm of velocity: MOVIE_VOLUME_TYPE == 8

  use constants_solver
  use shared_parameters, only: OUTPUT_FILES
  use specfem_par, only: it, scale_veloc

#ifdef USE_HDF5
  use specfem_par_movie_hdf5
#endif

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: veloc_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: veloc_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: veloc_inner_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

#ifdef USE_HDF5

  ! local parameters
  integer :: ispec,iglob,i,j,k,ier
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: tmp_data

  logical,parameter :: OUTPUT_CRUST_MANTLE = .true.
  logical,parameter :: OUTPUT_OUTER_CORE = .true.
  logical,parameter :: OUTPUT_INNER_CORE = .true.

  ! offset arrays
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_cm
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_oc
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_ic
  integer :: nspec_cm_all_proc,nspec_oc_all_proc,nspec_ic_all_proc

  ! gather number of points
  call gather_all_all_singlei(NSPEC_CRUST_MANTLE, offset_nspec_cm, NPROCTOT_VAL)
  call gather_all_all_singlei(NSPEC_OUTER_CORE, offset_nspec_oc, NPROCTOT_VAL)
  call gather_all_all_singlei(NSPEC_INNER_CORE, offset_nspec_ic, NPROCTOT_VAL)
  nspec_cm_all_proc = sum(offset_nspec_cm)
  nspec_oc_all_proc = sum(offset_nspec_oc)
  nspec_ic_all_proc = sum(offset_nspec_ic)

  ! initialize h5 file for volume movie
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = trim(OUTPUT_FILES) // '/movie_volume.h5'
  group_name = 'it_' // trim(i2c(it))

  ! create group and datasets
  if (myrank == 0) then
    call h5_open_file(file_name)
    call h5_open_or_create_group(group_name)

    if (OUTPUT_CRUST_MANTLE) then
      call h5_create_dataset_gen_in_group('reg1_veloc', (/NGLLX,NGLLY,NGLLZ,nspec_cm_all_proc/), 4, CUSTOM_REAL)
    endif
    if (OUTPUT_OUTER_CORE) then
      call h5_create_dataset_gen_in_group('reg2_veloc', (/NGLLX,NGLLY,NGLLZ,nspec_oc_all_proc/), 4, CUSTOM_REAL)
    endif
    if (OUTPUT_INNER_CORE) then
      call h5_create_dataset_gen_in_group('reg3_veloc', (/NGLLX,NGLLY,NGLLZ,nspec_ic_all_proc/), 4, CUSTOM_REAL)
    endif

    call h5_close_group()
    call h5_close_file()
  endif

  call synchronize_all()

  ! write the data
  call h5_open_file_p_collect(file_name)
  call h5_open_group(group_name)

  ! outputs norm of velocity
  if (OUTPUT_CRUST_MANTLE) then
    ! crust mantle
    ! these binary arrays can be converted into mesh format using the utility ./bin/xcombine_vol_data
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_CRUST_MANTLE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_crust_mantle(i,j,k,ispec)
            ! norm of velocity
            tmp_data(i,j,k,ispec) = real(scale_veloc,kind=CUSTOM_REAL) * sqrt( veloc_crust_mantle(1,iglob)**2 &
                                          + veloc_crust_mantle(2,iglob)**2 &
                                          + veloc_crust_mantle(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo

    call h5_write_dataset_collect_hyperslab_in_group('reg1_veloc', tmp_data, (/0,0,0,offset_nspec_cm(0:myrank-1)/), if_col)

    deallocate(tmp_data)
  endif

  if (OUTPUT_OUTER_CORE) then
    ! outer core
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_OUTER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_outer_core(i,j,k,ispec)
            ! norm of velocity
            ! note: this outputs only the first time derivative of the potential,
            !          not the actual velocity v = grad(Chi_dot)
            tmp_data(i,j,k,ispec) = abs(veloc_outer_core(iglob))
          enddo
        enddo
      enddo
    enddo

    call h5_write_dataset_collect_hyperslab_in_group('reg2_veloc', tmp_data, (/0,0,0,offset_nspec_oc(0:myrank-1)/), if_col)

    deallocate(tmp_data)
  endif

  if (OUTPUT_INNER_CORE) then
    ! inner core
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_INNER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_inner_core(i,j,k,ispec)
            ! norm of velocity
            tmp_data(i,j,k,ispec) = real(scale_veloc,kind=CUSTOM_REAL) * sqrt( veloc_inner_core(1,iglob)**2 &
                                          + veloc_inner_core(2,iglob)**2 &
                                          + veloc_inner_core(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo

    call h5_write_dataset_collect_hyperslab_in_group('reg3_veloc', tmp_data, (/0,0,0,offset_nspec_ic(0:myrank-1)/), if_col)

    deallocate(tmp_data)
  endif

  call h5_close_group()
  call h5_close_file_p()

#else

    print*, 'Error: HDF5 is not enabled in this version of the code.'
    print*, 'Please recompile with the HDF5 option enabled with --with-hdf5'
    stop

#endif

  end subroutine write_movie_volume_velnorm_hdf5



 subroutine write_movie_volume_accelnorm_hdf5(accel_crust_mantle,accel_inner_core,accel_outer_core, &
                                         ibool_crust_mantle,ibool_inner_core,ibool_outer_core)

! outputs norm of acceleration: MOVIE_VOLUME_TYPE == 1

  use constants_solver
  use shared_parameters, only: OUTPUT_FILES
  use specfem_par, only: it, scale_t_inv,scale_veloc

#ifdef USE_HDF5
  use specfem_par_movie_hdf5
#endif

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: accel_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: accel_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: accel_outer_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

#ifdef USE_HDF5
  ! local parameters
  integer :: ispec,iglob,i,j,k,ier
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: tmp_data
  real(kind=CUSTOM_REAL) :: scale_accel

  logical,parameter :: OUTPUT_CRUST_MANTLE = .true.
  logical,parameter :: OUTPUT_OUTER_CORE = .true.
  logical,parameter :: OUTPUT_INNER_CORE = .true.

  ! offset arrays
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_cm
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_oc
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_ic

  integer :: nspec_cm_all_proc,nspec_oc_all_proc,nspec_ic_all_proc

  ! gather number of points
  call gather_all_all_singlei(NSPEC_CRUST_MANTLE, offset_nspec_cm, NPROCTOT_VAL)
  call gather_all_all_singlei(NSPEC_OUTER_CORE, offset_nspec_oc, NPROCTOT_VAL)
  call gather_all_all_singlei(NSPEC_INNER_CORE, offset_nspec_ic, NPROCTOT_VAL)
  nspec_cm_all_proc = sum(offset_nspec_cm)
  nspec_oc_all_proc = sum(offset_nspec_oc)
  nspec_ic_all_proc = sum(offset_nspec_ic)

  ! dimensionalized scaling
  scale_accel = real(scale_veloc * scale_t_inv,kind=CUSTOM_REAL)

  ! initialize h5 file for volume movie
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = trim(OUTPUT_FILES) // '/movie_volume.h5'
  group_name = 'it_' // trim(i2c(it))

  ! create group and datasets
  if (myrank == 0) then
    call h5_open_file(file_name)
    call h5_open_or_create_group(group_name)

    if (OUTPUT_CRUST_MANTLE) then
      call h5_create_dataset_gen_in_group('reg1_accel', (/NGLLX,NGLLY,NGLLZ,nspec_cm_all_proc/), 4, CUSTOM_REAL)
    endif
    if (OUTPUT_OUTER_CORE) then
      call h5_create_dataset_gen_in_group('reg2_accel', (/NGLLX,NGLLY,NGLLZ,nspec_oc_all_proc/), 4, CUSTOM_REAL)
    endif
    if (OUTPUT_INNER_CORE) then
      call h5_create_dataset_gen_in_group('reg3_accel', (/NGLLX,NGLLY,NGLLZ,nspec_ic_all_proc/), 4, CUSTOM_REAL)
    endif

    call h5_close_group()
    call h5_close_file()
  endif

  call synchronize_all()

  ! write the data
  call h5_open_file_p_collect(file_name)
  call h5_open_group(group_name)

  ! outputs norm of acceleration
  if (OUTPUT_CRUST_MANTLE) then
    ! acceleration
    ! these binary arrays can be converted into mesh format using the utility ./bin/xcombine_vol_data
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_CRUST_MANTLE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_crust_mantle(i,j,k,ispec)
            ! norm
            tmp_data(i,j,k,ispec) = scale_accel * sqrt( accel_crust_mantle(1,iglob)**2 &
                                          + accel_crust_mantle(2,iglob)**2 &
                                          + accel_crust_mantle(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo

    call h5_write_dataset_collect_hyperslab_in_group('reg1_accel', tmp_data, (/0,0,0,offset_nspec_cm(0:myrank-1)/), if_col)

    deallocate(tmp_data)
  endif

  if (OUTPUT_OUTER_CORE) then
    ! outer core acceleration
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_OUTER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_outer_core(i,j,k,ispec)
            ! norm
            ! note: this outputs only the second time derivative of the potential,
            !          not the actual acceleration or pressure p = - rho * Chi_dot_dot
            tmp_data(i,j,k,ispec) = abs(accel_outer_core(iglob))
          enddo
        enddo
      enddo
    enddo

    call h5_write_dataset_collect_hyperslab_in_group('reg2_accel', tmp_data, (/0,0,0,offset_nspec_oc(0:myrank-1)/), if_col)

    deallocate(tmp_data)
  endif

  if (OUTPUT_INNER_CORE) then
    ! inner core
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_INNER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_inner_core(i,j,k,ispec)
            ! norm of acceleration
            tmp_data(i,j,k,ispec) = scale_accel * sqrt( accel_inner_core(1,iglob)**2 &
                                          + accel_inner_core(2,iglob)**2 &
                                          + accel_inner_core(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo

    call h5_write_dataset_collect_hyperslab_in_group('reg3_accel', tmp_data, (/0,0,0,offset_nspec_ic(0:myrank-1)/), if_col)

    deallocate(tmp_data)
  endif

#else

    print*, 'Error: HDF5 is not enabled in this version of the code.'
    print*, 'Please recompile with the HDF5 option enabled with --with-hdf5'
    stop

#endif

  end subroutine write_movie_volume_accelnorm_hdf5

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_HDF5
  subroutine get_conn_for_movie(elm_conn,offset,iNIT,npoints_3dmovie,num_ibool_3dmovie,mask_ibool_3dmovie)

  use specfem_par
  use specfem_par_crustmantle, only: ibool_crust_mantle
  use constants_solver


  implicit none

  integer, dimension(9,npoints_3dmovie), intent(out) :: elm_conn
  integer, intent(in) :: offset ! node id offset (starting global element id of each proc)
  integer, intent(in) :: iNIT
  integer, intent(in) :: npoints_3dmovie
  integer, dimension(NGLOB_CRUST_MANTLE_3DMOVIE),intent(in) :: num_ibool_3dmovie
  logical, dimension(NGLOB_CRUST_MANTLE_3DMOVIE),intent(in) :: mask_ibool_3dmovie


  ! local parameters
  integer :: ispec,ii
  integer,parameter :: cell_type = 9
  integer :: iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
  integer :: n1,n2,n3,n4,n5,n6,n7,n8
  integer :: i,j,k
  integer :: ispecele, iglob_center

  ispecele = 0
  do ispec = 1,NSPEC_CRUST_MANTLE

    ! checks center of element for movie flag
    iglob_center = ibool_crust_mantle((NGLLX+1)/2,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)

    ! checks if movie element
    if (mask_ibool_3dmovie(iglob_center)) then

      ! this element is in the movie region
      ispecele = ispecele+1

      do k = 1,NGLLZ-1,iNIT
        do j = 1,NGLLY-1,iNIT
          do i = 1,NGLLX-1,iNIT
            ! defines corners of a vtk element
            iglob1 = ibool_crust_mantle(i,j,k,ispec)
            iglob2 = ibool_crust_mantle(i+iNIT,j,k,ispec)
            iglob3 = ibool_crust_mantle(i+iNIT,j+iNIT,k,ispec)
            iglob4 = ibool_crust_mantle(i,j+iNIT,k,ispec)
            iglob5 = ibool_crust_mantle(i,j,k+iNIT,ispec)
            iglob6 = ibool_crust_mantle(i+iNIT,j,k+iNIT,ispec)
            iglob7 = ibool_crust_mantle(i+iNIT,j+iNIT,k+iNIT,ispec)
            iglob8 = ibool_crust_mantle(i,j+iNIT,k+iNIT,ispec)

            ! vtk indexing starts at 0 -> adds minus 1
            n1 = num_ibool_3dmovie(iglob1)-1
            n2 = num_ibool_3dmovie(iglob2)-1
            n3 = num_ibool_3dmovie(iglob3)-1
            n4 = num_ibool_3dmovie(iglob4)-1
            n5 = num_ibool_3dmovie(iglob5)-1
            n6 = num_ibool_3dmovie(iglob6)-1
            n7 = num_ibool_3dmovie(iglob7)-1
            n8 = num_ibool_3dmovie(iglob8)-1

            elm_conn(1, ispecele)  = cell_type
            elm_conn(2, ispecele)  = n1 + offset   ! node id starts 0 in xdmf rule
            elm_conn(3, ispecele)  = n2 + offset
            elm_conn(4, ispecele)  = n3 + offset
            elm_conn(5, ispecele)  = n4 + offset
            elm_conn(6, ispecele)  = n5 + offset
            elm_conn(7, ispecele)  = n6 + offset
            elm_conn(8, ispecele)  = n7 + offset
            elm_conn(9, ispecele)  = n8 + offset

            ! checks indices
            if (n1 < 0 .or. n2 < 0 .or. n3 < 0 .or. n4 < 0 .or. n5 < 0 .or. n6 < 0 .or. n7 < 0 .or. n8 < 0) then
              print *,'Error: movie element ',ispec,ispecele,'has invalid node index:',n1,n2,n3,n4,n5,n6,n7,n8
              call exit_mpi(myrank,'Error invalid movie element node index')
            endif

          enddo !i
        enddo !j
      enddo !k
    endif

  enddo !ispec

  end subroutine get_conn_for_movie

#endif



!
!-------------------------------------------------------------------------------------------------
!
#ifdef USE_HDF5

  subroutine write_xdmf_vol_hdf5(pressure_io, divglob_io, div_io, veloc_io, curl_io, stress_io)

  use specfem_par
  use specfem_par_movie_hdf5

  implicit none

  logical, intent(in) :: pressure_io, divglob_io, div_io, veloc_io, curl_io, stress_io

  ! local parameters
  integer :: i, ii
  character(len=20) :: it_str,nelm_str,nglo_str
  character(len=MAX_STRING_LEN) :: fname_xdmf_vol
  character(len=MAX_STRING_LEN) :: fname_h5_data_vol_xdmf

  ! checks if anything do, only main process writes out xdmf file
  if (myrank /= 0) return

  ! writeout xdmf file for volume movie
  fname_xdmf_vol = trim(OUTPUT_FILES) // "/movie_volume.xmf"
  fname_h5_data_vol_xdmf = "./movie_volume.h5"  ! relative to movie_volume.xmf file
  ! this seems to point to a wrong director:
  !   fname_h5_data_vol_xdmf = trim(OUTPUT_FILES) // "/movie_volume.h5"

  open(unit=xdmf_vol, file=trim(fname_xdmf_vol), recl=256)

  ! definition of topology and geometry
  ! refer only control nodes (8 or 27) as a coarse output
  ! data array need to be extracted from full data array on GLL points
  write(xdmf_vol,'(a)') '<?xml version="1.0" ?>'
  write(xdmf_vol,*) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(xdmf_vol,*) '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="3.0">'
  write(xdmf_vol,*) '<Domain name="mesh">'

  ! loop for writing information of mesh partitions
!  nelm_str = i2c(sum(nelm_par_proc_nio(:))*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1))
!  nglo_str = i2c(sum(nglob_par_proc_nio(:)))
!
!  write(xdmf_vol,*) '<Topology TopologyType="Mixed" NumberOfElements="'//trim(nelm_str)//'">'
!  write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Int" Precision="4" Dimensions="'&
!                         //trim(nelm_str)//' 9">'
!  write(xdmf_vol,*) '       '//trim(fname_h5_data_vol_xdmf)//':/mesh/elm_conn'
!  write(xdmf_vol,*) '</DataItem>'
!  write(xdmf_vol,*) '</Topology>'
!  write(xdmf_vol,*) '<Geometry GeometryType="X_Y_Z">'
!  write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
!                      //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
!  write(xdmf_vol,*) '       '//trim(fname_h5_data_vol_xdmf)//':/mesh/x'
!  write(xdmf_vol,*) '</DataItem>'
!  write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
!                      //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
!  write(xdmf_vol,*) '       '//trim(fname_h5_data_vol_xdmf)//':/mesh/y'
!  write(xdmf_vol,*) '</DataItem>'
!  write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
!                      //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
!  write(xdmf_vol,*) '       '//trim(fname_h5_data_vol_xdmf)//':/mesh/z'
!  write(xdmf_vol,*) '</DataItem>'
!  write(xdmf_vol,*) '</Geometry>'
!  write(xdmf_vol,*) '<Grid Name="time_col" GridType="Collection" CollectionType="Temporal">'
!
!  do i = 1, int(NSTEP/NTSTEP_BETWEEN_FRAMES)
!
!    ii = i*NTSTEP_BETWEEN_FRAMES
!    !write(it_str, "(i6.6)") ii
!    it_str = i2c(ii)
!
!    write(xdmf_vol,*) '<Grid Name="vol_mov" GridType="Uniform">'
!    write(xdmf_vol,*) '<Time Value="'//trim(r2c(sngl((ii-1)*DT-t0)))//'" />'
!    write(xdmf_vol,*) '<Topology Reference="/Xdmf/Domain/Topology" />'
!    write(xdmf_vol,*) '<Geometry Reference="/Xdmf/Domain/Geometry" />'
!
!    if (pressure_io) then
!      write(xdmf_vol,*) '<Attribute Name="pressure" AttributeType="Scalar" Center="Node">'
!      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
!                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
!      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/pressure'
!      write(xdmf_vol,*) '</DataItem>'
!      write(xdmf_vol,*) '</Attribute>'
!    endif
!
!    if (divglob_io) then
!      write(xdmf_vol,*) '<Attribute Name="div_glob" AttributeType="Scalar" Center="Node">'
!      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
!                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
!      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/div_glob'
!      write(xdmf_vol,*) '</DataItem>'
!      write(xdmf_vol,*) '</Attribute>'
!    endif
!
!    if (div_io) then
!      write(xdmf_vol,*) '<Attribute Name="div" AttributeType="Scalar" Center="Node">'
!      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
!                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
!      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/div'
!      write(xdmf_vol,*) '</DataItem>'
!      write(xdmf_vol,*) '</Attribute>'
!    endif
!
!    if (veloc_io) then
!      write(xdmf_vol,*) '<Attribute Name="veloc_x" AttributeType="Scalar" Center="Node">'
!      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
!                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
!      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/veloc_x'
!      write(xdmf_vol,*) '</DataItem>'
!      write(xdmf_vol,*) '</Attribute>'
!      write(xdmf_vol,*) '<Attribute Name="veloc_y" AttributeType="Scalar" Center="Node">'
!      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
!                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
!      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/veloc_y'
!      write(xdmf_vol,*) '</DataItem>'
!      write(xdmf_vol,*) '</Attribute>'
!      write(xdmf_vol,*) '<Attribute Name="veloc_z" AttributeType="Scalar" Center="Node">'
!      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
!                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
!      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/veloc_z'
!      write(xdmf_vol,*) '</DataItem>'
!      write(xdmf_vol,*) '</Attribute>'
!    endif
!
!    if (curl_io) then
!      write(xdmf_vol,*) '<Attribute Name="curl_x" AttributeType="Scalar" Center="Node">'
!      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
!                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
!      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/curl_x'
!      write(xdmf_vol,*) '</DataItem>'
!      write(xdmf_vol,*) '</Attribute>'
!      write(xdmf_vol,*) '<Attribute Name="curl_y" AttributeType="Scalar" Center="Node">'
!      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
!                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
!      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/curl_y'
!      write(xdmf_vol,*) '</DataItem>'
!      write(xdmf_vol,*) '</Attribute>'
!      write(xdmf_vol,*) '<Attribute Name="curl_z" AttributeType="Scalar" Center="Node">'
!      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
!                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
!      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/curl_z'
!      write(xdmf_vol,*) '</DataItem>'
!      write(xdmf_vol,*) '</Attribute>'
!    endif
!
!    if (stress_io) then
!      write(xdmf_vol,*) '<Attribute Name="stress_xx" AttributeType="Scalar" Center="Node">'
!      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
!                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
!      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/stress_xx'
!      write(xdmf_vol,*) '</DataItem>'
!      write(xdmf_vol,*) '</Attribute>'
!      write(xdmf_vol,*) '<Attribute Name="stress_yy" AttributeType="Scalar" Center="Node">'
!      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
!                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
!      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/stress_yy'
!      write(xdmf_vol,*) '</DataItem>'
!      write(xdmf_vol,*) '</Attribute>'
!      write(xdmf_vol,*) '<Attribute Name="stress_zz" AttributeType="Scalar" Center="Node">'
!      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
!                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
!      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/stress_zz'
!      write(xdmf_vol,*) '</DataItem>'
!      write(xdmf_vol,*) '</Attribute>'
!      write(xdmf_vol,*) '<Attribute Name="stress_xy" AttributeType="Scalar" Center="Node">'
!      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
!                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
!      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/stress_xy'
!      write(xdmf_vol,*) '</DataItem>'
!      write(xdmf_vol,*) '</Attribute>'
!
!      write(xdmf_vol,*) '<Attribute Name="stress_xz" AttributeType="Scalar" Center="Node">'
!      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
!                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
!      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/stress_xz'
!      write(xdmf_vol,*) '</DataItem>'
!      write(xdmf_vol,*) '</Attribute>'
!
!      write(xdmf_vol,*) '<Attribute Name="stress_yz" AttributeType="Scalar" Center="Node">'
!      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
!                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
!      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/stress_yz'
!      write(xdmf_vol,*) '</DataItem>'
!      write(xdmf_vol,*) '</Attribute>'
!    endif
!
!    write(xdmf_vol,*) '</Grid>'
!  enddo

  ! file finish
  write(xdmf_vol,*) '</Grid>'
  write(xdmf_vol,*) '</Domain>'
  write(xdmf_vol,*) '</Xdmf>'

  close(xdmf_vol)

  end subroutine write_xdmf_vol_hdf5

#endif
