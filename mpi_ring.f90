! Copyright (C) 2024 Alexander Buccheri alexanderbuccheri@googlemail.com
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! MPI DEMO: Ring comms pattern using send/receive


! mpif90 -pedantic -fbacktrace -fbounds-check -Wl,-syslibroot,/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk mpi_ring.f90 -o mpi_ring
! mpirun -np 2 ./mpi_ring
module mpi_routines
  implicit none
  public


contains

  !> @brief Number of elements on process 'rank' for an array of n elements,
  !! distributed over np processes.
  integer function distribute_elements(np, n, rank)
    integer, intent(in) :: n               !< Number of elements to distribute
    integer, intent(in) :: np              !< Number of processes to distribute across
    integer, intent(in) :: rank            !< Process to compute local number of elements
    integer, allocatable :: indices(:, :)  !< Start and stop indices for each each process
    indices = distribute_elements_start_stop(np, n)
    distribute_elements = indices(2, rank + 1) - indices(1, rank + 1) + 1
  end function distribute_elements

  !> @brief Distribute N elements over NP processes.
  !!
  !! Example 1: n = 12 and np = 3
  !! [1, 2, 3, 4 | 5, 6, 7, 8 | 9, 10, 11, 12]
  !!
  !! Example 2: n = 14 and np = 4
  !! If n is not integer-divisible by np, assign the remainders equally between the first n_remainder processes
  !! [1, 2, 3, 4 | 5, 6, 7, 8 |9, 10, 11 |12, 13, 14]
  !!
  function distribute_elements_start_stop(np, n) result(indices)
    integer, intent(in) :: n       !< Number of elements to distribute
    integer, intent(in) :: np      !< Number of processes to distribute across
    integer :: elements_per_process, n_remainder
    integer :: n_assigned_processes, n_assigned_elements
    integer :: iprocess, istart, istop
    integer :: indices(2, np)

    ! When load-balanced
    elements_per_process = int(n / np)
    n_remainder = mod(n, np)

    ! Assign processes with extra elements first
    iprocess = 0
    do while (n_remainder > 0)
       iprocess = iprocess + 1
       istart = 1 + (iprocess - 1)*(elements_per_process + 1)
       istop = iprocess * (elements_per_process + 1)
       indices(:, iprocess) = [istart, istop]
       n_remainder = n_remainder - 1
    end do

    n_assigned_processes = iprocess
    n_assigned_elements = (elements_per_process + 1) * n_assigned_processes

    do iprocess = n_assigned_processes + 1, np
       istart = 1 + ((iprocess - n_assigned_processes - 1) * elements_per_process) + n_assigned_elements
       istop = ((iprocess - n_assigned_processes) * elements_per_process) + n_assigned_elements
       indices(:, iprocess) = [istart, istop]
    end do

  end function distribute_elements_start_stop

end module mpi_routines


program mpi_ring
  use mpi

  ! My modules
  use  mpi_routines
  implicit none

  integer :: rank, np, ierr, iprocess, send_id, rec_id, iround, send_index
  integer, allocatable :: collated_data(:), send_to(:), receive_from(:), initial_send_to(:), initial_receive_from(:), receive_send(:, :)
  integer, parameter :: tag = 101
  integer :: status(MPI_STATUS_SIZE)

  integer :: i, cnt, n_values, n_local, ij
  integer, allocatable :: local_state_indices(:), indices(:, :), pair_indices(:, :)

  call MPI_INIT(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

  ! Start send-tos and receive-froms by permuting this process's rank
  initial_send_to = modulo([(i, i = 0, np-1)] + 1, np)
  initial_receive_from = modulo([(i, i = 0, np-1)] - 1, np)

  ! Defines receive-from and send-to ranks for each round of the ring pattern
  ! Defined specifically for each process
  allocate(receive_send(2, np -1))
  
  ! Round comm pattern
  do iround = 1, np - 1
     send_to = modulo(initial_send_to + iround - 1, np)
     receive_from = modulo(initial_receive_from - (iround - 1), np)
     receive_send(:, iround) = [receive_from(rank+1), send_to(rank+1)]
  enddo

  do iround = 1, np - 1
     write(*, *) "** Round ", iround, " **"
     if(rank ==0) write(*, *) 'receive-from, rank, send-to', receive_send(1, iround), rank, receive_send(2, iround)
     if(rank ==1) write(*, *) 'receive-from, rank, send-to', receive_send(1, iround), rank, receive_send(2, iround)
     if(rank ==2) write(*, *) 'receive-from, rank, send-to', receive_send(1, iround), rank, receive_send(2, iround)
 enddo


  
  call MPI_FINALIZE(ierr)


end program mpi_ring
