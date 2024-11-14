module io_bandwidth
    use mpi, only: MPI_Wtime

    implicit none
    double precision :: start_time, end_time, time_delta
    integer :: bytes_written

  contains

    subroutine initialize_bytes_written()
      bytes_written = 0
      time_delta = 0.0
    end subroutine initialize_bytes_written

    subroutine start_timer()
        start_time = MPI_Wtime()
    end subroutine start_timer

    subroutine stop_timer()
        end_time = MPI_Wtime()
        time_delta = time_delta + (end_time - start_time)
    end subroutine stop_timer

    subroutine set_bytes_written(bytes)
      implicit none
      integer, intent(in) :: bytes
      bytes_written = bytes
    end subroutine set_bytes_written

    subroutine set_bytes_written_from_array(element_size, num_elements)
        implicit none
        integer, intent(in):: element_size, num_elements

        ! Calculate the total bytes written and add to the existing value
        bytes_written = bytes_written + element_size * num_elements

      end subroutine set_bytes_written_from_array


    subroutine calculate_bandwidth_this_proc()
      implicit none
      double precision :: elapsed_time, bandwidth

      elapsed_time = time_delta
      if (elapsed_time > 0.0) then
        bandwidth = bytes_written / (elapsed_time * 1.0e6)  ! Bandwidth in MB/s
        print *, 'I/O Bandwidth: ', bandwidth, ' MB/s'
      else
        print *, 'Elapsed time is zero, cannot calculate bandwidth'
      endif
    end subroutine calculate_bandwidth_this_proc

    ! calculate total bandwidth across all processes
    subroutine calculate_bandwidth_all_procs()
      use specfem_par
      use shared_parameters
      implicit none
      integer :: i,total_bytes
      double precision :: elapsed_time, max_elapsed_time, &
              bandwidth, total_bandwidth


      print *, '-------------------------------------------------------------------------------'

      elapsed_time = time_delta
      if (elapsed_time > 0.0) then
        bandwidth = bytes_written / (elapsed_time * 1.0e6)  ! Bandwidth in MB/s

        ! calculate total bytes written across all processes
        call sum_all_all_i(bytes_written, total_bytes)
        ! get the maximum elapsed time across all processes
        call max_all_dp(elapsed_time, max_elapsed_time)

        ! calculate total bandwidth
        total_bandwidth = total_bytes / (max_elapsed_time * 1.0e6) ! Bandwidth in MB/s

        !print*, "DEBUG: elapsed and max_elapsed_time = ", elapsed_time, max_elapsed_time

        do i = 0, NPROCTOT_VAL-1
          if (myrank == i) then
            print *, 'process ', myrank, ' bytes_written = ', bytes_written, ' elapsed_time (s) = ', &
                        elapsed_time, ' bandwidth = ', bandwidth, ' MB/s'
          end if
        end do

        ! write total bandwidth
        if (myrank == 0) then
          print *, 'total bytes_written = ', total_bytes, ' max_elapsed_time (s) = ', max_elapsed_time, &
                        ' total bandwidth = ', total_bandwidth, ' MB/s'
        end if
        call synchronize_all()


      else
        print *, 'Elapsed time is zero, cannot calculate bandwidth'
      endif
    end subroutine calculate_bandwidth_all_procs

  end module io_bandwidth