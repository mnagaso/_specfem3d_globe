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
  public :: wait_all_send

  ! public parameters
  public :: dest_ionod

  public :: nproc_io

  ! verbosity
  public :: VERBOSE

  private

  ! responsible id of io node
  integer :: dest_ionod = 0
  ! number of computer nodes sending info to each io node
  integer :: nproc_io
  ! id for io node
  integer :: my_io_id

  ! verbose output (for debugging)
  logical, parameter :: VERBOSE = .true.


! USE_HDF5
#endif

contains

  subroutine initialize_io_server()

  ! initialization of IO server splits compute and io nodes
#ifdef USE_HDF5

  use constants, only: IMAIN,myrank,MAX_STRING_LEN
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,NPROCTOT

  implicit none

  integer :: sizeval,irank
  integer :: mpi_comm,split_comm,inter_comm,ier
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
          print *, "io_server: rank ", myrank, " compute task, my_local_mpi_comm_inter = ", inter_comm, " my_local_mpi_comm_world = ", split_comm
        else
          print *, "io_server: rank ", myrank, " io task, my_local_mpi_comm_inter = ", inter_comm, " my_local_mpi_comm_world = ", split_comm
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

    implicit none

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

  subroutine pass_info_to_io()

#ifdef USE_HDF5

    implicit none

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
!
!-------------------------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------------------------
!
! MPI communications
!
!-------------------------------------------------------------------------------------------------


  subroutine wait_all_send()

#ifdef USE_HDF5

  implicit none

  integer :: ireq

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




end module io_server_hdf5
