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
! mpirun -np 3 ./mpi_ring
module mpi_routines
  implicit none
  public

contains

  !> @brief Produce an index array of receive-from and send-to ranks per iteration of a ring pattern.
  !!
  !! \p receive_send has ranks defined for the local process \p rank.
  !!
  !! Example. Running with three processes, one needs two iterations of a ring for each process
  !! to exchange data with all other processes:
  !! 
  !! iround   receive_from    my_rank    send_to
  !!   1            2            0          1
  !!   1            0            1          2
  !!   1            1            2          0
  !!   --------------------------------------------
  !!   2            1            0          2
  !!   2            2            1          0
  !!   2            0            2          1
  subroutine mpi_ring_pattern_receive_send_ranks(rank, np, receive_send)
    integer,              intent(in ) :: rank               !< MPI rank
    integer,              intent(in ) :: np                 !< Number of MPI processes
    integer, allocatable, intent(out) :: receive_send(:, :) !<  receive-from and send-to ranks
    
    integer, allocatable :: initial_receive_from(:), initial_send_to(:),  receive_from(:), send_to(:)
    integer              :: nring, iround, i

    !  Number of ring pattern iterations
    nring = np - 1 
    
    ! Start send-tos and receive-froms by permuting this process's rank
    initial_send_to = modulo([(i, i = 0, np-1)] + 1, np)
    initial_receive_from = modulo([(i, i = 0, np-1)] - 1, np)
    
    ! Defines receive-from and send-to ranks for each round of the ring pattern
    ! Defined specifically for this local MPI process
    allocate(receive_send(2, nring))
    
    ! Round comm pattern
    do iround = 1, nring
       send_to = modulo(initial_send_to + iround - 1, np)
       receive_from = modulo(initial_receive_from - (iround - 1), np)
       receive_send(:, iround) = [receive_from(rank+1), send_to(rank+1)]
    enddo

  end subroutine mpi_ring_pattern_receive_send_ranks
  

end module mpi_routines


program mpi_ring
  use mpi

  ! My modules
  use  mpi_routines
  implicit none

  integer :: reqs(2), stats(MPI_STATUS_SIZE, 2)
  integer :: rank, np, ierr

  integer :: iprocess, iround
  integer :: irecv, isend, rf, st, tag
  integer :: send_buffer, recv_buffer
  
  integer, allocatable :: receive_send(:, :)

  call MPI_INIT(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

  call  mpi_ring_pattern_receive_send_ranks(rank, np, receive_send)

  ! Start rounds with non-blocking communication
  do iround = 1, size(receive_send, 2)
     ! Prepare buffers
     send_buffer = rank
     recv_buffer = -1
     tag = iround

     ! Asynchronous receive
     rf = receive_send(1, iround)
     call MPI_IRecv(recv_buffer, 1, MPI_INTEGER, rf, tag, MPI_COMM_WORLD, reqs(1), ierr)

     ! Asynchronous send
     st = receive_send(2, iround)
     call MPI_ISend(send_buffer, 1, MPI_INTEGER, st, tag, MPI_COMM_WORLD, reqs(2), ierr)

     ! Wait for both send and receive to complete
     call MPI_Waitall(2, reqs, stats, ierr)

     ! Output the result
     write(*, *) "** Round ", iround, " Rank ", rank, " **"
     write(*, *) 'receive-from, rank, send-to', rf, rank, st
     write(*, *) 'Received value:', recv_buffer

  enddo
  
  call MPI_FINALIZE(ierr)

end program mpi_ring
