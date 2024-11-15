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

module io_server_hdf5

  use constants, only: CUSTOM_REAL
  use shared_parameters, only: HDF5_IO_NODES, &
    IO_storage_task, IO_compute_task

  implicit none

#ifdef USE_HDF5

  ! public routines
  public :: initialize_io_server
  public :: finalize_io_server
  public :: do_io_start_idle
  public :: pass_info_to_io
  public :: get_info_from_comp
  public :: wait_all_send

  ! functions to prepare the offset arrays for the HDF5 output
  public :: initialize_hdf5_solver
  public :: finalize_hdf5_solver
  public :: allocate_offset_arrays
  public :: deallocate_offset_arrays
  public :: create_hdf5_files_and_datasets

  ! public parameters
  public :: dest_ionod

  ! MPI message tags
  ! undo attenuation
  public :: io_tag_ford_undo_d_cm, io_tag_ford_undo_v_cm, io_tag_ford_undo_a_cm
  public :: io_tag_ford_undo_d_oc, io_tag_ford_undo_v_oc, io_tag_ford_undo_a_oc
  public :: io_tag_ford_undo_d_ic, io_tag_ford_undo_v_ic, io_tag_ford_undo_a_ic
  public :: io_tag_ford_undo_eps_xx_cm, io_tag_ford_undo_eps_yy_cm, io_tag_ford_undo_eps_xy_cm
  public :: io_tag_ford_undo_eps_xz_cm, io_tag_ford_undo_eps_yz_cm
  public :: io_tag_ford_undo_eps_xx_ic, io_tag_ford_undo_eps_yy_ic, io_tag_ford_undo_eps_xy_ic
  public :: io_tag_ford_undo_eps_xz_ic, io_tag_ford_undo_eps_yz_ic
  public :: io_tag_ford_undo_A_rot, io_tag_ford_undo_B_rot
  public :: io_tag_ford_undo_R_xx_cm, io_tag_ford_undo_R_yy_cm, io_tag_ford_undo_R_xy_cm
  public :: io_tag_ford_undo_R_xz_cm, io_tag_ford_undo_R_yz_cm
  public :: io_tag_ford_undo_R_xx_ic, io_tag_ford_undo_R_yy_ic, io_tag_ford_undo_R_xy_ic
  public :: io_tag_ford_undo_R_xz_ic, io_tag_ford_undo_R_yz_ic
  public :: io_tag_ford_undo_neq, io_tag_ford_undo_neq1, io_tag_ford_undo_pgrav1
  public :: io_tag_ford_undo_nmsg

  ! MPI requests
  public :: n_req_ford_undo
  public :: req_dump_ford_undo
  public :: n_msg_ford_undo

  public :: nproc_io

  ! verbosity
  public :: VERBOSE



  private

  ! MPI tags for io server implementation
  ! for forward undo att arrays
  integer :: io_tag_ford_undo_d_cm      = 100001
  integer :: io_tag_ford_undo_v_cm      = 100002
  integer :: io_tag_ford_undo_a_cm      = 100003
  integer :: io_tag_ford_undo_d_oc      = 100004
  integer :: io_tag_ford_undo_v_oc      = 100005
  integer :: io_tag_ford_undo_a_oc      = 100006
  integer :: io_tag_ford_undo_d_ic      = 100007
  integer :: io_tag_ford_undo_v_ic      = 100008
  integer :: io_tag_ford_undo_a_ic      = 100009
  integer :: io_tag_ford_undo_eps_xx_cm = 100010
  integer :: io_tag_ford_undo_eps_yy_cm = 100011
  integer :: io_tag_ford_undo_eps_xy_cm = 100012
  integer :: io_tag_ford_undo_eps_xz_cm = 100013
  integer :: io_tag_ford_undo_eps_yz_cm = 100014
  integer :: io_tag_ford_undo_eps_xx_ic = 100015
  integer :: io_tag_ford_undo_eps_yy_ic = 100016
  integer :: io_tag_ford_undo_eps_xy_ic = 100017
  integer :: io_tag_ford_undo_eps_xz_ic = 100018
  integer :: io_tag_ford_undo_eps_yz_ic = 100019
  integer :: io_tag_ford_undo_A_rot     = 100020
  integer :: io_tag_ford_undo_B_rot     = 100021
  integer :: io_tag_ford_undo_R_xx_cm   = 100022
  integer :: io_tag_ford_undo_R_yy_cm   = 100023
  integer :: io_tag_ford_undo_R_xy_cm   = 100024
  integer :: io_tag_ford_undo_R_xz_cm   = 100025
  integer :: io_tag_ford_undo_R_yz_cm   = 100026
  integer :: io_tag_ford_undo_R_xx_ic   = 100027
  integer :: io_tag_ford_undo_R_yy_ic   = 100028
  integer :: io_tag_ford_undo_R_xy_ic   = 100029
  integer :: io_tag_ford_undo_R_xz_ic   = 100030
  integer :: io_tag_ford_undo_R_yz_ic   = 100031
  integer :: io_tag_ford_undo_neq       = 100032
  integer :: io_tag_ford_undo_neq1      = 100033
  integer :: io_tag_ford_undo_pgrav1    = 100034
  integer :: io_tag_nsubset_iterations  = 100035
  integer :: io_tag_ford_undo_nmsg      = 100036

  ! mpi_req dump (used in wait_all_send)
  integer :: n_req_ford_undo = 0
  integer, dimension(34) :: req_dump_ford_undo
  integer :: n_msg_ford_undo = 0

  ! responsible id of io node
  integer :: dest_ionod = 0
  ! number of computer nodes sending info to each io node
  integer :: nproc_io
  ! id for io node
  integer :: my_io_id

  ! verbose output (for debugging)
  logical, parameter :: VERBOSE = .false.


! USE_HDF5
#endif

contains

  subroutine initialize_io_server()

  ! initialization of IO server splits compute and io nodes
#ifdef USE_HDF5

  use constants, only: IMAIN,myrank,MAX_STRING_LEN
  use shared_parameters, only: NPROCTOT !,NUMBER_OF_SIMULTANEOUS_RUNS

  implicit none

  integer :: sizeval,irank
  integer :: mpi_comm,split_comm,inter_comm
  integer :: key,io_start,comp_start
  ! test node name
  character(len=MAX_STRING_LEN), dimension(:), allocatable :: node_names
  character(len=MAX_STRING_LEN) :: this_node_name
  integer :: node_len

  ! split comm into computation nodes and io node
  ! here we use the last HDF5_IO_NODES ranks (intra comm) as the io node
  ! for using one additional node, xspecfem3D need to be run with + HDF5_IO_NODES node
  ! thus for running mpirun, it should be like e.g.
  ! mpirun -n $((NPROC+HDF5_IO_NODES)) ./bin/xspecfem3D

  ! safety check
  if (HDF5_IO_NODES < 0) stop 'Invalid HDF5_IO_NODES, must be zero or positive'

  ! initializes
  call world_get_comm(mpi_comm)         ! mpi_comm == my_local_mpi_comm_world

  ! default MPI group, no i/o server
  call world_set_comm_inter(mpi_comm)   ! my_local_mpi_comm_inter = mpi_comm

  ! checks if anything to do
  if (HDF5_IO_NODES == 0) return

  ! select the task of this proc
  ! get the local mpi_size and rank
  call world_size(sizeval)
  call world_rank(myrank)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'HDF5 I/O server run:'
    write(IMAIN,*) '  number of dedicated HDF5 IO nodes = ',HDF5_IO_NODES
    write(IMAIN,*) '  total number of MPI processes     = ',sizeval
    write(IMAIN,*)
    write(IMAIN,*) '  separating subgroups for compute tasks with ',sizeval - HDF5_IO_NODES,'MPI processes'
    write(IMAIN,*) '                       for io tasks      with ',HDF5_IO_NODES,'MPI processes'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! checks if solver run with a number of MPI processes equal to (NPROC + HDF5_IO_NODES)
  if (sizeval /= NPROCTOT + HDF5_IO_NODES) then
    if (myrank == 0) then
      print *,"Error: HDF5 IO server needs ",HDF5_IO_NODES," additional process together with NPROCTOT",NPROCTOT,"processes"
      print *
      print *,"However, this simulation runs with ",sizeval,"MPI processes,"
      print *,"where instead it should run with ",NPROCTOT + HDF5_IO_NODES,"MPI processes"
      print *
      print *,"Please check your call to 'mpirun -np .. xspecfem3D' and rerun with the correct number of processes"
      print *
    endif
    stop 'Invalid number of MPI processes for HDF5_IO_NODES > 0'
  endif

  ! to select io node and compute nodes on the same cluster node.
  allocate(node_names(0:sizeval-1))

  call world_get_processor_name(this_node_name,node_len)

  ! share the node name between procs
  call gather_all_all_single_ch(this_node_name,node_names,sizeval,MAX_STRING_LEN)

  call synchronize_all()

  ! select the task of this proc
  call select_io_node(node_names,myrank,sizeval,key,io_start)

  deallocate(node_names)

  ! split communicator into compute_comm and io_comm
  call world_comm_split(mpi_comm, key, myrank, split_comm)

  ! create inter communicator and set as my_local_mpi_comm_inter
  if (IO_storage_task) then
    comp_start = 0
    call world_create_intercomm(split_comm, 0, mpi_comm, comp_start, 11111, inter_comm)
  else
    !io_start = sizeval - HDF5_IO_NODES !+dest_ionod
    call world_create_intercomm(split_comm, 0, mpi_comm, io_start, 11111, inter_comm)
  endif

  ! sets new comm_world subgroup
  call world_set_comm(split_comm)            ! such that: my_local_mpi_comm_world == split_comm

  ! use inter_comm as my_local_mpi_comm_world for all send/recv
  call world_set_comm_inter(inter_comm)      ! such that: my_local_mpi_comm_inter == inter_comm

  ! exclude io nodes from the other compute nodes
  !if (NUMBER_OF_SIMULTANEOUS_RUNS > 1) NPROCTOT = NPROCTOT - HDF5_IO_NODES

  ! re-define myrank (within new comm_world subgroup)
  call world_rank(myrank)
  call world_size(sizeval)

  ! debug print inter_comm and if IO_compute_task
  if (VERBOSE) then
    do irank = 0, sizeval-1
      if (myrank == irank) then
        if (IO_compute_task) then
          print *, "io_server: rank ", myrank, " compute task, my_local_mpi_comm_inter = ", &
            inter_comm, " my_local_mpi_comm_world = ", split_comm
        else
          print *, "io_server: rank ", myrank, " io task, my_local_mpi_comm_inter = ", &
            inter_comm, " my_local_mpi_comm_world = ", split_comm
        endif
      endif
      call flush_stdout()
      call synchronize_all()
    enddo
  endif

  call synchronize_inter()

  ! note: IO compute tasks cover the lower part of MPI processes, i.e., with ranks in
  !       the initial 0 to (sizeval-HDF5_IO_NODES)-1 range.
  !       the io nodes with IO_storage_task set to .true. are added at the end, with ranks in
  !       the upper range (sizeval-HDF5_IO_NODES) to (sizeval-1).
  !
  !       for user output, we only have the initial rank 0 as the main process opening the output_solver.txt file.
  !       after separating the tasks and creating new subgroup, we have a rank 0 process in the compute subgroup
  !       as well as a rank 0 process in the storage subgroup.
  !
  !       from here on, we only want the process with rank 0 in the compute group write output to IMAIN.
  !
  ! user output
  if (IO_compute_task) then
    if (myrank == 0) then
      write(IMAIN,*) '  new number of processes in compute group = ',sizeval
      write(IMAIN,*) '  new addtional processes in io group      = ',HDF5_IO_NODES
      write(IMAIN,*) '  MPI i/o server setup done'
      call flush_IMAIN()
    endif
  endif

  ! allocate the offset arrays
  !call initialize_hdf5_solver()

#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine initialize_io_server() called without HDF5 Support."
  print *, "       HDF5_IO_NODES > 0 requires HDF5 support and HDF5_ENABLED must be set to .true."
  print *
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'initialize_io_server() called without HDF5 compilation support'

#endif

  end subroutine initialize_io_server

  !
  !-------------------------------------------------------------------------------------------------
  !

  subroutine finalize_io_server()

#ifdef USE_HDF5

  implicit none

  ! checks if anything to do
  if (HDF5_IO_NODES == 0) return

  ! finish MPI subgroup
  if (HDF5_IO_NODES > 0) then
    ! deallocate offset arrays
    !call deallocate_offset_arrays()
    ! wait for all to finish
    call synchronize_inter()
    ! free subgroup
    call world_comm_free_inter()
  endif
#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine finalize_io_server() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'finalize_io_server() called without HDF5 compilation support'

#endif

  end subroutine finalize_io_server

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_HDF5

  subroutine select_io_node(node_names, myrank, sizeval, key, io_start)

  use constants, only: MAX_STRING_LEN

  implicit none

  integer, intent(in)                         :: myrank, sizeval
  integer, intent(out)                        :: key, io_start
  character(len=MAX_STRING_LEN), dimension(0:sizeval-1), intent(in) :: node_names
  ! local parameters
  integer, dimension(sizeval) :: n_procs_on_node ! number of procs on each cluster node
  integer, dimension(:), allocatable :: n_ionode_on_cluster ! number of ionode on the cluster nodes
  integer :: i,j,c,n_cluster_node,my_cluster_id,n_rest_io,n_ionode,n_comp_node
  real(kind=CUSTOM_REAL) :: io_ratio ! dum
  character(len=MAX_STRING_LEN), dimension(sizeval) :: dump_node_names ! names of cluster nodes

  ! initialize
  my_cluster_id = -1
  n_cluster_node = 0
  n_procs_on_node(:) = 0
  dump_node_names(:) = "nan"

  ! BUG: nprocs goes wrong when 2 cluster nodes and 64 compute 8 io procs
  ! (only 62 compute nodes are assigned)

  ! cluster_node_nums = [n_1,...,n_i,...n_cn,-1,-1,...] ! n_i is the number of procs on
  ! each cluster node
  do i = 1, sizeval
    ! search the node name already found and registered in the dump_node_names array
    c = 0
    do j = 1, sizeval
      if (node_names(i-1) == dump_node_names(j)) c = j
    enddo

    ! if node_name[i-1] is not registered yet
    if (c == 0) then
      ! count the number of cluster node
      n_cluster_node = n_cluster_node + 1

      dump_node_names(n_cluster_node) = node_names(i-1) ! register the name
      n_procs_on_node(n_cluster_node) = 1 ! count up the number of procs on this cluster node

      c = n_cluster_node
    else ! if node_name[i] has already be found, count up the number of procs
      n_procs_on_node(c) = n_procs_on_node(c) + 1
    endif

    ! the id of cluster which this process belongs to
    if (i-1 == myrank) then
      my_cluster_id = c
    endif
  enddo

  ! warning when a cluster node has only one single proc.
  do i = 1, n_cluster_node
    if (n_procs_on_node(i) == 1) then
      print *
      print *, "***************************************************************************"
      print *, "IO server Warning:"
      print *, "  node name: " // trim(dump_node_names(i)) // " has only one procs."
      print *, "  this may lead a io performance issue by inter-clusternode communication."
      print *, "***************************************************************************"
      print *
    endif
  enddo

  ! select HDF5_IO_NODES of io nodes
  allocate(n_ionode_on_cluster(n_cluster_node))
  !! decide the number of io node on each cluster node
  ! at least one io node on each cluster node
  n_ionode_on_cluster(:) = 1

  ! check if the total number of io node > HDF5_IO_NODES
  if (sum(n_ionode_on_cluster) > HDF5_IO_NODES) then
    print *, "Error: HDF5_IO_NODES in Par_file is too small,"
    print *, "       at least one io node for each cluster node is necessary"
    stop 'Invalid HDF5_IO_NODES value too small'
  endif

  ! share the rest of ionodes based on the ratio of compute nodes
  n_rest_io = HDF5_IO_NODES - sum(n_ionode_on_cluster)
  if (n_rest_io /= 0) then
    ! add io nodes one by one to the cluster node where
    ! the io_node/total_node ratio is rowest
    do i = 1, n_rest_io
      io_ratio = 9999.0 !initialize at each i
      do j = 1, n_cluster_node
        if (io_ratio > real(n_ionode_on_cluster(j))/real(n_procs_on_node(j))) then
          ! dump if largest
          io_ratio = real(n_ionode_on_cluster(j))/real(n_procs_on_node(j))
          ! destination of additional io node
          c = j
        endif
      enddo
      ! add one additional io node
      n_ionode_on_cluster(c) = n_ionode_on_cluster(c) + 1
    enddo
  endif

  !! choose the io node from the last rank of each cluster node
  n_ionode = 0
  do i = 1, n_cluster_node
    c = 0
    ! number of compute node on this cluster node
    n_comp_node = n_procs_on_node(i) - n_ionode_on_cluster(i)
    do j = 1, sizeval
      ! if this rank is on i cluster node
      if (dump_node_names(i) == node_names(j-1)) then
        c = c + 1
        if (n_comp_node < c) then
          ! j is io node
          if (j-1 == myrank) then
            ! check  if j is myrank
            IO_storage_task = .true. ! set io node flag
            IO_compute_task = .false.
            key          = 0
            my_io_id     = n_ionode ! id of io_node
            nproc_io     = n_comp_node/n_ionode_on_cluster(i) ! number of compute nodes which use this io node
            dest_ionod   = -1
            if (c-n_comp_node-1 < mod(n_comp_node,n_ionode_on_cluster(i))) nproc_io = nproc_io + 1
          endif

          ! rank of io_start
          if (n_ionode == 0) io_start = j-1

          n_ionode = n_ionode+1

        else
          ! j is compute node
          if (j-1 == myrank) then
            IO_storage_task = .false.
            IO_compute_task = .true.
            key          = 1
            ! set the destination of MPI communication
            dest_ionod   = mod(c-1,n_ionode_on_cluster(i)) + n_ionode
          endif
        endif
      endif
    enddo
  enddo

  ! debug
  if (VERBOSE) then
    if (myrank == 0) then
      print *, "io_server: n_procs_on_node", n_procs_on_node(:)
      print *, "io_server: n_ionode_on_cluster", n_ionode_on_cluster
      print *
    endif
    call flush_stdout()
    call synchronize_all()

    if (IO_storage_task) then
      print *, "io_server: io_task rank ", myrank, " nprocio ",nproc_io
      print *, "io_server: io_task rank ", myrank, " node ", trim(node_names(myrank)), " my_io_id", my_io_id
      print *
    endif
    call flush_stdout
    call synchronize_all()

    do i = 0, n_ionode-1
      if (.not. IO_storage_task) then
        if (dest_ionod == i) then
          print *, "io_server: rank ", myrank, " node ", trim(node_names(myrank)), " dest ", dest_ionod
        endif
      endif
      call flush_stdout()
      call synchronize_all()
    enddo
    call flush_stdout()
    call synchronize_all()

  endif

  deallocate(n_ionode_on_cluster)

  end subroutine select_io_node

#endif

!
!-------------------------------------------------------------------------------------------------
!

  subroutine do_io_start_idle()

#ifdef USE_HDF5

  use specfem_par
  use specfem_par_movie_hdf5
  use constants, only: myrank, my_status_size, my_status_source, my_status_tag

  use io_bandwidth

  implicit none

  integer :: status(my_status_size)
  integer :: tag, tag_src

  ! vars forward undo att arrays
  ! undo attenuation
  integer :: n_recv_msg_ford_undo, &      ! number of messages received for undo attenuation of one iteration
             max_ford_undo_out, &         ! number of iterations when IO happens for undo attenuation
             ford_undo_out_count          ! count the completed iterations for undo attenuation

  ! array for dumping the data array
  real(kind=CUSTOM_REAL), dimension(:),         allocatable :: dump_ford_undo_1d_glob
  real(kind=CUSTOM_REAL), dimension(:,:),       allocatable :: dump_ford_undo_2d_glob
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),   allocatable :: dump_ford_undo_4d
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: dump_ford_undo_5d

  ! maximum nglob and nspec in offset arrays
  integer :: max_nglob
  integer :: max_nspec

  integer :: i_out

  ! initialize all counters
  ! undo attenuation
  n_recv_msg_ford_undo = 0 ! number of messages received for undo attenuation of one iteration
  max_ford_undo_out    = 0 ! number of iterations when IO happens for undo attenuation
  ford_undo_out_count  = 0 ! count the completed iterations for undo attenuation

  ! create all the HDF5 files and datasets
  do i_out = 1, NSUBSET_ITERATIONS
    call create_hdf5_files_and_datasets(i_out)
  enddo

  !
  ! initialize
  !

  ! undo attenuation
  if (UNDO_ATTENUATION) then

    n_msg_ford_undo   = n_msg_ford_undo*nproc_io ! multiply by the number of compute nodes for each io node
    max_ford_undo_out = NSUBSET_ITERATIONS ! multiply by the number of snapshots

    ! find the maximum number of nglob in offset_nglob_cm, _oc, _ic
    max_nglob = maxval([maxval(offset_nglob_cm), maxval(offset_nglob_oc), maxval(offset_nglob_ic)])

    ! find the maximum number of nspec in offset_nspec_cm_soa, _ic_soa, _oc_rot, _cm_att, _ic_att
    max_nspec = maxval([maxval(offset_nspec_cm_soa), maxval(offset_nspec_ic_soa), &
                        merge(maxval(offset_nspec_oc_rot), -huge(1), ROTATION_VAL), &
                        merge(maxval(offset_nspec_cm_att), -huge(1), ATTENUATION_VAL), &
                        merge(maxval(offset_nspec_ic_att), -huge(1), ATTENUATION_VAL)])

    ! allocate the dump arrays
    allocate(dump_ford_undo_1d_glob(max_nglob), &
             dump_ford_undo_2d_glob(NDIM, max_nglob), &
             dump_ford_undo_4d(NGLLX, NGLLY, NGLLZ, max_nspec), &
             dump_ford_undo_5d(NGLLX, NGLLY, NGLLZ, N_SLS, max_nspec))

  endif ! UNDO_ATTENUATION


  ! initialize timer
  call initialize_bytes_written()

  !
  ! idling loop
  !
  do while (ford_undo_out_count < max_ford_undo_out)

    ! check the iteration counter


    ! waiting for a MPI message
    call idle_mpi_io(status)

    tag = status(my_status_tag)
    tag_src = status(my_status_source)

    ! debug output on the received message
    if (VERBOSE) then
      print *, "io_server: rank ", myrank, " received message with tag ", tag, " from rank ", tag_src, &
               " counters, ford undo: ", n_recv_msg_ford_undo, " / ", n_msg_ford_undo
    endif

    ! undo attenuation
    if (UNDO_ATTENUATION) then
      ! receive the data
      call recv_and_write_ford_undo(tag, tag_src, status, &
                                    dump_ford_undo_1d_glob, &
                                    dump_ford_undo_2d_glob, &
                                    dump_ford_undo_4d, &
                                    dump_ford_undo_5d, &
                                    ford_undo_out_count) ! use for the filename

      ! count 1 message received
      n_recv_msg_ford_undo = n_recv_msg_ford_undo + 1

    endif ! UNDO_ATTENUATION

    ! check receive counters
    if (n_recv_msg_ford_undo >= n_msg_ford_undo) then
      ! increment the iteration counter
      ford_undo_out_count = ford_undo_out_count + 1
      ! reset the receive counter
      n_recv_msg_ford_undo = 0

      ! write out the times
      call calculate_bandwidth_all_procs()

      ! re-initialize timer
      call initialize_bytes_written()

    endif


  enddo
  !
  ! end of idling loop
  !

  ! deallocate temporary arrays
  ! undo attenuation
  if (UNDO_ATTENUATION) then
    deallocate(dump_ford_undo_1d_glob, &
               dump_ford_undo_2d_glob, &
               dump_ford_undo_4d, &
               dump_ford_undo_5d)
  endif


#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine do_io_start_idle() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'pass_info_to_io() called without HDF5 compilation support'

#endif

  end subroutine do_io_start_idle

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_info_from_comp()

#ifdef USE_HDF5

  use specfem_par
  use specfem_par_movie_hdf5

  implicit none

  integer, dimension(1) :: tmp_arr

  ! pass necessary information to io node

  ! forward undo att
  if (UNDO_ATTENUATION) then
    ! receive NSUBSET_ITERATIONS
    call recv_i_inter(tmp_arr, 1, 0, io_tag_nsubset_iterations)
    NSUBSET_ITERATIONS = tmp_arr(1)

    ! receive offset arrays
    call recv_i_inter(offset_nglob_cm, NPROCTOT_VAL, 0, io_tag_ford_undo_d_cm) !!!!!!!!!!!!!!
    call recv_i_inter(offset_nglob_oc, NPROCTOT_VAL, 0, io_tag_ford_undo_d_oc)
    call recv_i_inter(offset_nglob_ic, NPROCTOT_VAL, 0, io_tag_ford_undo_d_ic)
    call recv_i_inter(offset_nspec_cm_soa, NPROCTOT_VAL, 0, io_tag_ford_undo_eps_xx_cm)
    call recv_i_inter(offset_nspec_ic_soa, NPROCTOT_VAL, 0, io_tag_ford_undo_eps_xx_ic)
    if (ROTATION_VAL) then
      call recv_i_inter(offset_nspec_oc_rot, NPROCTOT_VAL, 0, io_tag_ford_undo_A_rot)
    endif
    if (ATTENUATION_VAL) then
      call recv_i_inter(offset_nspec_cm_att, NPROCTOT_VAL, 0, io_tag_ford_undo_R_xx_cm)
      call recv_i_inter(offset_nspec_ic_att, NPROCTOT_VAL, 0, io_tag_ford_undo_R_xx_ic)
    endif
    if (FULL_GRAVITY_VAL) then
      call recv_i_inter(offset_pgrav1, NPROCTOT_VAL, 0, io_tag_ford_undo_pgrav1)
    endif

    ! receive the number of messages
    call recv_i_inter(tmp_arr, 1, 0, io_tag_ford_undo_nmsg)
    n_msg_ford_undo = tmp_arr(1)

  endif ! UNDO_ATTENUATION

#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine pass_info_to_io() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'pass_info_to_io() called without HDF5 compilation support'

#endif

  end subroutine get_info_from_comp
!
!-------------------------------------------------------------------------------------------------
!

  subroutine pass_info_to_io()

#ifdef USE_HDF5

    use specfem_par
    use specfem_par_movie_hdf5

    implicit none

    integer :: i_ionod
    integer :: tmp_int

    ! pass necessary information to io node

    ! forward undo att
    if (UNDO_ATTENUATION) then
      if (myrank == 0) then

        ! count up the number of messages
        n_msg_ford_undo = 19
        if (ROTATION_VAL) n_msg_ford_undo = n_msg_ford_undo + 2
        if (ATTENUATION_VAL) n_msg_ford_undo = n_msg_ford_undo + 10
        if (FULL_GRAVITY_VAL) n_msg_ford_undo = n_msg_ford_undo + 3

        do i_ionod = 0, HDF5_IO_NODES-1

          ! send NSUBSET_ITERATIONS
          call send_i_inter((/NSUBSET_ITERATIONS/), 1, i_ionod, io_tag_nsubset_iterations)

          ! send offset arrays
          call send_i_inter(offset_nglob_cm, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_d_cm) !!!!!!!!!!!!!!
          call send_i_inter(offset_nglob_oc, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_d_oc)
          call send_i_inter(offset_nglob_ic, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_d_ic)
          call send_i_inter(offset_nspec_cm_soa, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_eps_xx_cm)
          call send_i_inter(offset_nspec_ic_soa, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_eps_xx_ic)
          if (ROTATION_VAL) then
            call send_i_inter(offset_nspec_oc_rot, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_A_rot)
          endif
          if (ATTENUATION_VAL) then
            call send_i_inter(offset_nspec_cm_att, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_R_xx_cm)
            call send_i_inter(offset_nspec_ic_att, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_R_xx_ic)
          endif
          if (FULL_GRAVITY_VAL) then
            call send_i_inter(offset_pgrav1, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_pgrav1)
          endif

          ! send the number of messages
          tmp_int = n_msg_ford_undo
          call send_i_inter((/tmp_int/), 1, i_ionod, io_tag_ford_undo_nmsg)

        enddo ! i_ionod
      endif ! myrank == 0

    endif ! UNDO_ATTENUATION

#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine pass_info_to_io() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'pass_info_to_io() called without HDF5 compilation support'

#endif

  end subroutine pass_info_to_io

!-------------------------------------------------------------------------------------------------
!
! MPI communications
!
!-------------------------------------------------------------------------------------------------

  subroutine wait_all_send()

#ifdef USE_HDF5

  implicit none

  integer :: ireq

  ! forward undo arrays
  if (n_req_ford_undo /= 0) then
    ! wait till all mpi_isends are finished
    do ireq = 1,n_req_ford_undo
      call wait_req(req_dump_ford_undo(ireq))
    enddo
  endif
  n_req_ford_undo = 0

  ! surface movie
!  if (n_req_surf /= 0) then
!    ! wait till all mpi_isends are finished
!    do ireq = 1,n_req_surf
!      call wait_req(req_dump_surf(ireq))
!    enddo
!  endif
!
!  ! volume movie
!  if (n_req_vol /= 0) then
!    ! wait till all mpi_isends are finished
!    do ireq = 1,n_req_vol
!      call wait_req(req_dump_vol(ireq))
!    enddo
!  endif
!
!  n_req_surf = 0; n_req_vol = 0

  call synchronize_all()

#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine wait_all_send() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'wait_all_send() called without HDF5 compilation support'

#endif

  end subroutine wait_all_send

!
!-------------------------------------------------------------------------------------------------
!
  subroutine allocate_offset_arrays()

#ifdef USE_HDF5
    use specfem_par
    use specfem_par_movie_hdf5

    implicit none

    integer :: ier

    ! nglobs
    allocate(offset_nglob_cm(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nglob_cm')
    allocate(offset_nglob_oc(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nglob_oc')
    allocate(offset_nglob_ic(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nglob_ic')

    ! nelems
    allocate(offset_nspec_cm(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_cm')
    allocate(offset_nspec_oc(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_oc')
    allocate(offset_nspec_ic(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_ic')
    allocate(offset_nspec_cm_soa(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_cm_soa')
    allocate(offset_nspec_ic_soa(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_ic_soa')
    if (ROTATION_VAL) then
      allocate(offset_nspec_oc_rot(0:NPROCTOT_VAL-1),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_oc_rot')
    endif
    if (ATTENUATION_VAL) then
      allocate(offset_nspec_cm_att(0:NPROCTOT_VAL-1),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_cm_att')
      allocate(offset_nspec_ic_att(0:NPROCTOT_VAL-1),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_ic_att')
    endif

    if (FULL_GRAVITY_VAL) then
      allocate(offset_pgrav1(0:NPROCTOT_VAL-1),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_pgrav1')
    endif

    npoints_vol_mov_all_proc_cm     = 0
    npoints_vol_mov_all_proc_oc     = 0
    npoints_vol_mov_all_proc_ic     = 0

    nspec_vol_mov_all_proc_cm     = 0
    nspec_vol_mov_all_proc_oc     = 0
    nspec_vol_mov_all_proc_ic     = 0
    nspec_vol_mov_all_proc_cm_soa = 0
    nspec_vol_mov_all_proc_ic_soa = 0
    if (ROTATION_VAL) then
      nspec_vol_mov_all_proc_oc_rot = 0
    endif
    if (ATTENUATION_VAL) then
      nspec_vol_mov_all_proc_cm_att = 0
      nspec_vol_mov_all_proc_ic_att = 0
    endif

#else

  print *
  print *, "Error: HDF5 I/O server routine allocate_offset_arrays() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'allocate_offset_arrays() called without HDF5 compilation support'

#endif

  end subroutine allocate_offset_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine deallocate_offset_arrays()

#ifdef USE_HDF5

    use specfem_par
    use specfem_par_movie_hdf5

    implicit none

    deallocate(offset_nglob_cm)
    deallocate(offset_nglob_oc)
    deallocate(offset_nglob_ic)

    deallocate(offset_nspec_cm)
    deallocate(offset_nspec_oc)
    deallocate(offset_nspec_ic)
    deallocate(offset_nspec_cm_soa)
    deallocate(offset_nspec_ic_soa)
    if (ROTATION_VAL) then
      deallocate(offset_nspec_oc_rot)
    endif
    if (ATTENUATION_VAL) then
      deallocate(offset_nspec_cm_att)
      deallocate(offset_nspec_ic_att)
    endif
    if (FULL_GRAVITY_VAL) then
      deallocate(offset_pgrav1)
    endif

#else

  print *
  print *, "Error: HDF5 I/O server routine deallocate_offset_arrays() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'deallocate_offset_arrays() called without HDF5 compilation support'

#endif

  end subroutine deallocate_offset_arrays

!
!-------------------------------------------------------------------------------------------------
!



  subroutine initialize_hdf5_solver()

#ifdef USE_HDF5

    use specfem_par
    use specfem_par_movie_hdf5
    use specfem_par_full_gravity, only : neq1

    implicit none

    ! allocate offset arrays
    call allocate_offset_arrays()

    call gather_all_all_singlei(NGLOB_CRUST_MANTLE, offset_nglob_cm, NPROCTOT_VAL)
    call gather_all_all_singlei(NGLOB_OUTER_CORE,   offset_nglob_oc, NPROCTOT_VAL)
    call gather_all_all_singlei(NGLOB_INNER_CORE,   offset_nglob_ic, NPROCTOT_VAL)

    call gather_all_all_singlei(NSPEC_CRUST_MANTLE, offset_nspec_cm, NPROCTOT_VAL)
    call gather_all_all_singlei(NSPEC_OUTER_CORE,   offset_nspec_oc, NPROCTOT_VAL)
    call gather_all_all_singlei(NSPEC_INNER_CORE,   offset_nspec_ic, NPROCTOT_VAL)
    call gather_all_all_singlei(NSPEC_CRUST_MANTLE_STR_OR_ATT, offset_nspec_cm_soa, NPROCTOT_VAL)
    call gather_all_all_singlei(NSPEC_INNER_CORE_STR_OR_ATT,   offset_nspec_ic_soa, NPROCTOT_VAL)
    if (ROTATION_VAL) then
      call gather_all_all_singlei(NSPEC_OUTER_CORE_ROTATION, offset_nspec_oc_rot, NPROCTOT_VAL)
    endif
    if (ATTENUATION_VAL) then
      call gather_all_all_singlei(NSPEC_CRUST_MANTLE_ATTENUATION, offset_nspec_cm_att, NPROCTOT_VAL)
      call gather_all_all_singlei(NSPEC_INNER_CORE_ATTENUATION,   offset_nspec_ic_att, NPROCTOT_VAL)
    endif
    if (FULL_GRAVITY_VAL) then
      call gather_all_all_singlei(neq1, offset_pgrav1, NPROCTOT_VAL)
    endif

    npoints_vol_mov_all_proc_cm = sum(offset_nglob_cm)
    npoints_vol_mov_all_proc_oc = sum(offset_nglob_oc)
    npoints_vol_mov_all_proc_ic = sum(offset_nglob_ic)

    nspec_vol_mov_all_proc_cm = sum(offset_nspec_cm)
    nspec_vol_mov_all_proc_oc = sum(offset_nspec_oc)
    nspec_vol_mov_all_proc_ic = sum(offset_nspec_ic)
    nspec_vol_mov_all_proc_cm_soa = sum(offset_nspec_cm_soa)
    nspec_vol_mov_all_proc_ic_soa = sum(offset_nspec_ic_soa)
    if (ROTATION_VAL) then
      nspec_vol_mov_all_proc_oc_rot = sum(offset_nspec_oc_rot)
    endif
    if (ATTENUATION_VAL) then
      nspec_vol_mov_all_proc_cm_att = sum(offset_nspec_cm_att)
      nspec_vol_mov_all_proc_ic_att = sum(offset_nspec_ic_att)
    endif

#else

  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine initialize_hdf5_solver() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'initialize_hdf5_solver() called without HDF5 compilation support'

#endif

  end subroutine initialize_hdf5_solver

!
!-------------------------------------------------------------------------------------------------
!

  subroutine finalize_hdf5_solver()

#ifdef USE_HDF5

    use specfem_par
    use specfem_par_movie_hdf5

    implicit none

    ! deallocate offset arrays
    call deallocate_offset_arrays()

#else

  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine finalize_hdf5_solver() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'finalize_hdf5_solver() called without HDF5 compilation support'

#endif

  end subroutine finalize_hdf5_solver

!
!-------------------------------------------------------------------------------------------------
!

  subroutine create_hdf5_files_and_datasets(i_snapshot)


#ifdef USE_HDF5
    use specfem_par
    use specfem_par_movie_hdf5
#endif

    implicit none

    integer, intent(in) :: i_snapshot

#ifdef USE_HDF5
    ! forward undo arrays
    if (UNDO_ATTENUATION) then

      ! create output files and datasets
      write(file_name, '(a,i6.6,a)') 'save_frame_at',i_snapshot,'.h5'
      file_name = trim(LOCAL_PATH)//'/'//trim(file_name)

      ! get MPI parameters
      call world_get_comm(comm)
      call world_get_info_null(info)

      ! initialize HDF5
      call h5_initialize() ! called in initialize_mesher()
      ! set MPI
      call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

      ! create file and datasets by myrank==0
      if (myrank == 0) then
        call h5_create_file(file_name)

        ! create datasets
        call h5_create_dataset_gen('displ_crust_mantle', (/NDIM, sum(offset_nglob_cm)/), 2, CUSTOM_REAL)
        call h5_create_dataset_gen('veloc_crust_mantle', (/NDIM, sum(offset_nglob_cm)/), 2, CUSTOM_REAL)
        call h5_create_dataset_gen('accel_crust_mantle', (/NDIM, sum(offset_nglob_cm)/), 2, CUSTOM_REAL)
        call h5_create_dataset_gen('displ_outer_core', (/sum(offset_nglob_oc)/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('veloc_outer_core', (/sum(offset_nglob_oc)/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('accel_outer_core', (/sum(offset_nglob_oc)/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('displ_inner_core', (/NDIM, sum(offset_nglob_ic)/), 2, CUSTOM_REAL)
        call h5_create_dataset_gen('veloc_inner_core', (/NDIM, sum(offset_nglob_ic)/), 2, CUSTOM_REAL)
        call h5_create_dataset_gen('accel_inner_core', (/NDIM, sum(offset_nglob_ic)/), 2, CUSTOM_REAL)
        call h5_create_dataset_gen('epsilondev_xx_crust_mantle', &
                                   (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm_soa)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('epsilondev_yy_crust_mantle', &
                                   (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm_soa)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('epsilondev_xy_crust_mantle', &
                                   (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm_soa)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('epsilondev_xz_crust_mantle', &
                                   (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm_soa)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('epsilondev_yz_crust_mantle', &
                                   (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm_soa)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('epsilondev_xx_inner_core', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('epsilondev_yy_inner_core', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('epsilondev_xy_inner_core', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('epsilondev_xz_inner_core', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('epsilondev_yz_inner_core', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)

        if (ROTATION_VAL) then
          call h5_create_dataset_gen('A_array_rotation', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_oc_rot)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen('B_array_rotation', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_oc_rot)/), 4, CUSTOM_REAL)
        endif

        if (ATTENUATION_VAL) then
          call h5_create_dataset_gen('R_xx_crust_mantle', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_cm_att)/), 5, CUSTOM_REAL)
          call h5_create_dataset_gen('R_yy_crust_mantle', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_cm_att)/), 5, CUSTOM_REAL)
          call h5_create_dataset_gen('R_xy_crust_mantle', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_cm_att)/), 5, CUSTOM_REAL)
          call h5_create_dataset_gen('R_xz_crust_mantle', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_cm_att)/), 5, CUSTOM_REAL)
          call h5_create_dataset_gen('R_yz_crust_mantle', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_cm_att)/), 5, CUSTOM_REAL)

          call h5_create_dataset_gen('R_xx_inner_core', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_ic_att)/), 5, CUSTOM_REAL)
          call h5_create_dataset_gen('R_yy_inner_core', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_ic_att)/), 5, CUSTOM_REAL)
          call h5_create_dataset_gen('R_xy_inner_core', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_ic_att)/), 5, CUSTOM_REAL)
          call h5_create_dataset_gen('R_xz_inner_core', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_ic_att)/), 5, CUSTOM_REAL)
          call h5_create_dataset_gen('R_yz_inner_core', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_ic_att)/), 5, CUSTOM_REAL)
        endif ! ATTENUATION_VAL

        if (FULL_GRAVITY_VAL) then
          call h5_create_dataset_gen('neq', (/NPROCTOT_VAL/), 1, 1)
          call h5_create_dataset_gen('neq1', (/NPROCTOT_VAL/), 1, 1)
          call h5_create_dataset_gen('pgrav1', (/sum(offset_pgrav1)/), 1, CUSTOM_REAL)
        endif ! FULL_GRAVITY_VAL

        ! close file
        call h5_close_file()
      endif ! myrank == 0

    endif ! UNDO_ATTENUATION

    ! wait for rank 0 to finish
    call synchronize_all()

#else

  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine create_files_and_datasets_ford_undo() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'create_files_and_datasets_ford_undo() called without HDF5 compilation support'

#endif

  end subroutine create_hdf5_files_and_datasets


!-------------------------------------------------------------------------------
!
! HDF5 io routines (only available with HDF5 compilation support)
!
!-------------------------------------------------------------------------------

#if defined(USE_HDF5)
  ! only available with HDF5 compilation support

  subroutine idle_mpi_io(status)
  ! wait for an arrival of any MPI message

    use constants, only: my_status_size

    implicit none

    integer, intent(inout) :: status(my_status_size)

    call world_probe_any_inter(status)

  end subroutine idle_mpi_io

#endif
  !
  !-------------------------------------------------------------------------------------------------
  !

#if defined(USE_HDF5)

  subroutine recv_and_write_ford_undo(tag, tag_src, status, &
                                     dump_ford_undo_1d_glob, &
                                     dump_ford_undo_2d_glob, &
                                     dump_ford_undo_4d, &
                                     dump_ford_undo_5d, &
                                     i_snapshot)

    use specfem_par
    use specfem_par_movie_hdf5
    use manager_hdf5
    use constants, only: CUSTOM_REAL, my_status_size, my_status_source, my_status_tag
    use io_bandwidth

    implicit none

    integer, intent(in) :: tag, tag_src
    integer, intent(in) :: status(my_status_size)
    real(kind=CUSTOM_REAL), dimension(:),         intent(inout) :: dump_ford_undo_1d_glob
    real(kind=CUSTOM_REAL), dimension(:,:),       intent(inout) :: dump_ford_undo_2d_glob
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),   intent(inout) :: dump_ford_undo_4d
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(inout) :: dump_ford_undo_5d
    integer, intent(in) :: i_snapshot

    integer :: msg_size, ista, iend, req_dummy, data_len
    integer :: neq, neq1

    integer :: data_size ! size of one data element in bytes

    ! default datasize is for CUSTOM_REAL
    data_size = CUSTOM_REAL ! 8 for double precision, 4 for single precision

    ! get message size
    call world_get_size_msg(status, msg_size)

    ! file name to write (save_frame_atXXXXXX.h5) XXXXXX is the snapshot number + 1
    write(file_name, '(a,i6.6,a)') 'save_frame_at',i_snapshot+1,'.h5'
    file_name = trim(LOCAL_PATH)//'/'//trim(file_name)

    ! open file
    !! get MPI parameters
    !call world_get_comm(comm)
    !call world_get_info_null(info)

    !! initialize HDF5
    !call h5_initialize() ! called in initialize_mesher()
    !! set MPI
    !call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

    if (H5_COL) then
      ! open file
      call h5_open_file_p_collect(file_name)
    else
      ! open file
      call h5_open_file_p(file_name)
    endif

    ! receive the data
    if (tag == io_tag_ford_undo_d_cm) then
      ! displ_cust_mantle
      ista = sum(offset_nglob_cm(0:tag_src-1))
      iend = sum(offset_nglob_cm(0:tag_src))
      data_len = offset_nglob_cm(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_2d_glob(:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('displ_crust_mantle', dump_ford_undo_2d_glob(:,1:data_len), (/0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_v_cm) then
      ! veloc_crust_mantle
      ista = sum(offset_nglob_cm(0:tag_src-1))
      iend = sum(offset_nglob_cm(0:tag_src))
      data_len = offset_nglob_cm(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_2d_glob(:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('veloc_crust_mantle', dump_ford_undo_2d_glob(:,1:data_len), (/0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_a_cm) then
      ! accel_crust_mantle
      ista = sum(offset_nglob_cm(0:tag_src-1))
      iend = sum(offset_nglob_cm(0:tag_src))
      data_len = offset_nglob_cm(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_2d_glob(:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('accel_crust_mantle', dump_ford_undo_2d_glob(:,1:data_len), (/0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_d_oc) then
      ! displ_outer_core
      ista = sum(offset_nglob_oc(0:tag_src-1))
      iend = sum(offset_nglob_oc(0:tag_src))
      data_len = offset_nglob_oc(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_1d_glob(1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('displ_outer_core', dump_ford_undo_1d_glob(1:data_len), (/ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_v_oc) then
      ! veloc_outer_core
      ista = sum(offset_nglob_oc(0:tag_src-1))
      iend = sum(offset_nglob_oc(0:tag_src))
      data_len = offset_nglob_oc(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_1d_glob(1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('veloc_outer_core', dump_ford_undo_1d_glob(1:data_len), (/ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_a_oc) then
      ! accel_outer_core
      ista = sum(offset_nglob_oc(0:tag_src-1))
      iend = sum(offset_nglob_oc(0:tag_src))
      data_len = offset_nglob_oc(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_1d_glob(1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('accel_outer_core', dump_ford_undo_1d_glob(1:data_len), (/ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_d_ic) then
      ! displ_inner_core
      ista = sum(offset_nglob_ic(0:tag_src-1))
      iend = sum(offset_nglob_ic(0:tag_src))
      data_len = offset_nglob_ic(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_2d_glob(:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('displ_inner_core', dump_ford_undo_2d_glob(:,1:data_len), (/0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_v_ic) then
      ! veloc_inner_core
      ista = sum(offset_nglob_ic(0:tag_src-1))
      iend = sum(offset_nglob_ic(0:tag_src))
      data_len = offset_nglob_ic(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_2d_glob(:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('veloc_inner_core', dump_ford_undo_2d_glob(:,1:data_len), (/0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_a_ic) then
      ! accel_inner_core
      ista = sum(offset_nglob_ic(0:tag_src-1))
      iend = sum(offset_nglob_ic(0:tag_src))
      data_len = offset_nglob_ic(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_2d_glob(:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('accel_inner_core', dump_ford_undo_2d_glob(:,1:data_len), (/0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_xx_cm) then
      ! epsilondev_xx_crust_mantle
      ista = sum(offset_nspec_cm_soa(0:tag_src-1))
      iend = sum(offset_nspec_cm_soa(0:tag_src))
      data_len = offset_nspec_cm_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_xx_crust_mantle', dump_ford_undo_4d(:,:,:,1:data_len), &
                                                                                              (/0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_yy_cm) then
      ! epsilondev_yy_crust_mantle
      ista = sum(offset_nspec_cm_soa(0:tag_src-1))
      iend = sum(offset_nspec_cm_soa(0:tag_src))
      data_len = offset_nspec_cm_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_yy_crust_mantle', dump_ford_undo_4d(:,:,:,1:data_len), &
                                                                                              (/0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_xy_cm) then
      ! epsilondev_xy_crust_mantle
      ista = sum(offset_nspec_cm_soa(0:tag_src-1))
      iend = sum(offset_nspec_cm_soa(0:tag_src))
      data_len = offset_nspec_cm_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_xy_crust_mantle', dump_ford_undo_4d(:,:,:,1:data_len), &
                                                                                              (/0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_xz_cm) then
      ! epsilondev_xz_crust_mantle
      ista = sum(offset_nspec_cm_soa(0:tag_src-1))
      iend = sum(offset_nspec_cm_soa(0:tag_src))
      data_len = offset_nspec_cm_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_xz_crust_mantle', dump_ford_undo_4d(:,:,:,1:data_len), &
                                                                                              (/0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_yz_cm) then
      ! epsilondev_yz_crust_mantle
      ista = sum(offset_nspec_cm_soa(0:tag_src-1))
      iend = sum(offset_nspec_cm_soa(0:tag_src))
      data_len = offset_nspec_cm_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_yz_crust_mantle', dump_ford_undo_4d(:,:,:,1:data_len), &
                                                                                              (/0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_xx_ic) then
      ! epsilondev_xx_inner_core
      ista = sum(offset_nspec_ic_soa(0:tag_src-1))
      iend = sum(offset_nspec_ic_soa(0:tag_src))
      data_len = offset_nspec_ic_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_xx_inner_core', dump_ford_undo_4d(:,:,:,1:data_len), &
                                                                                            (/0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_yy_ic) then
      ! epsilondev_yy_inner_core
      ista = sum(offset_nspec_ic_soa(0:tag_src-1))
      iend = sum(offset_nspec_ic_soa(0:tag_src))
      data_len = offset_nspec_ic_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_yy_inner_core', dump_ford_undo_4d(:,:,:,1:data_len), &
                                                                                            (/0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_xy_ic) then
      ! epsilondev_xy_inner_core
      ista = sum(offset_nspec_ic_soa(0:tag_src-1))
      iend = sum(offset_nspec_ic_soa(0:tag_src))
      data_len = offset_nspec_ic_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_xy_inner_core', dump_ford_undo_4d(:,:,:,1:data_len), &
                                                                                            (/0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_xz_ic) then
      ! epsilondev_xz_inner_core
      ista = sum(offset_nspec_ic_soa(0:tag_src-1))
      iend = sum(offset_nspec_ic_soa(0:tag_src))
      data_len = offset_nspec_ic_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_xz_inner_core', dump_ford_undo_4d(:,:,:,1:data_len), &
                                                                                            (/0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_yz_ic) then
      ! epsilondev_yz_inner_core
      ista = sum(offset_nspec_ic_soa(0:tag_src-1))
      iend = sum(offset_nspec_ic_soa(0:tag_src))
      data_len = offset_nspec_ic_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_yz_inner_core', dump_ford_undo_4d(:,:,:,1:data_len), &
                                                                                            (/0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_A_rot) then
      ! A_array_rotation
      ista = sum(offset_nspec_oc_rot(0:tag_src-1))
      iend = sum(offset_nspec_oc_rot(0:tag_src))
      data_len = offset_nspec_oc_rot(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('A_array_rotation', dump_ford_undo_4d(:,:,:,1:data_len), &
                                                                                    (/0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_B_rot) then
      ! B_array_rotation
      ista = sum(offset_nspec_oc_rot(0:tag_src-1))
      iend = sum(offset_nspec_oc_rot(0:tag_src))
      data_len = offset_nspec_oc_rot(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('B_array_rotation', dump_ford_undo_4d(:,:,:,1:data_len), &
                                                                                    (/0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_xx_cm) then
      ! R_xx_crust_mantle
      ista = sum(offset_nspec_cm_att(0:tag_src-1))
      iend = sum(offset_nspec_cm_att(0:tag_src))
      data_len = offset_nspec_cm_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_xx_crust_mantle', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                                                                    (/0, 0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_yy_cm) then
      ! R_yy_crust_mantle
      ista = sum(offset_nspec_cm_att(0:tag_src-1))
      iend = sum(offset_nspec_cm_att(0:tag_src))
      data_len = offset_nspec_cm_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_yy_crust_mantle', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                                                                    (/0, 0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_xy_cm) then
      ! R_xy_crust_mantle
      ista = sum(offset_nspec_cm_att(0:tag_src-1))
      iend = sum(offset_nspec_cm_att(0:tag_src))
      data_len = offset_nspec_cm_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_xy_crust_mantle', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                                                                    (/0, 0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_xz_cm) then
      ! R_xz_crust_mantle
      ista = sum(offset_nspec_cm_att(0:tag_src-1))
      iend = sum(offset_nspec_cm_att(0:tag_src))
      data_len = offset_nspec_cm_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_xz_crust_mantle', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                                                                    (/0, 0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_yz_cm) then
      ! R_yz_crust_mantle
      ista = sum(offset_nspec_cm_att(0:tag_src-1))
      iend = sum(offset_nspec_cm_att(0:tag_src))
      data_len = offset_nspec_cm_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_yz_crust_mantle', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                                                                    (/0, 0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_xx_ic) then
      ! R_xx_inner_core
      ista = sum(offset_nspec_ic_att(0:tag_src-1))
      iend = sum(offset_nspec_ic_att(0:tag_src))
      data_len = offset_nspec_ic_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_xx_inner_core', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                                                                  (/0, 0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_yy_ic) then
      ! R_yy_inner_core
      ista = sum(offset_nspec_ic_att(0:tag_src-1))
      iend = sum(offset_nspec_ic_att(0:tag_src))
      data_len = offset_nspec_ic_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_yy_inner_core', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                                                                  (/0, 0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_xy_ic) then
      ! R_xy_inner_core
      ista = sum(offset_nspec_ic_att(0:tag_src-1))
      iend = sum(offset_nspec_ic_att(0:tag_src))
      data_len = offset_nspec_ic_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_xy_inner_core', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                                                                  (/0, 0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_xz_ic) then
      ! R_xz_inner_core
      ista = sum(offset_nspec_ic_att(0:tag_src-1))
      iend = sum(offset_nspec_ic_att(0:tag_src))
      data_len = offset_nspec_ic_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_xz_inner_core', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                                                                  (/0, 0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_yz_ic) then
      ! R_yz_inner_core
      ista = sum(offset_nspec_ic_att(0:tag_src-1))
      iend = sum(offset_nspec_ic_att(0:tag_src))
      data_len = offset_nspec_ic_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_yz_inner_core', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                                                                  (/0, 0, 0, 0, ista/), H5_COL)
      call stop_timer()

    else if (tag == io_tag_ford_undo_neq) then
      ! neq
      ! receive
      call irecv_i_inter((/neq/), 1, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('neq', (/neq/), (/tag_src/), H5_COL)
      call stop_timer()

      ! data size if for integer
      data_size = 4

    else if (tag == io_tag_ford_undo_neq1) then
      ! neq1
      ! receive
      call irecv_i_inter((/neq1/), 1, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('neq1', (/neq1/), (/tag_src/), H5_COL)
      call stop_timer()

      ! data size if for integer
      data_size = 4

    else if (tag == io_tag_ford_undo_pgrav1) then
      ! pgrav1
      ista = sum(offset_pgrav1(0:tag_src-1))
      iend = sum(offset_pgrav1(0:tag_src))
      data_len = offset_pgrav1(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_1d_glob(1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      call start_timer()
      call h5_write_dataset_collect_hyperslab('pgrav1', dump_ford_undo_1d_glob(1:data_len), (/ista/), H5_COL)
      call stop_timer()
    else
      ! unknown tag
      print *, 'Error: unknown tag in recv_and_write_ford_undo'
      stop 'Error: unknown tag in recv_and_write_ford_undo'

    end if


    ! count the bytes written
    call set_bytes_written_from_array(data_size*8, msg_size) ! converting to bits from bytes

    ! close file
    call h5_close_file_p()
    !call h5_close_file()

  end subroutine recv_and_write_ford_undo

!
!-------------------------------------------------------------------------------------------------
!

#endif

end module io_server_hdf5
