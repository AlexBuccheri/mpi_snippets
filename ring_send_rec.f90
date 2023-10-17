! Copyright (C) 2023 Alexander Buccheri alexanderbuccheri@googlemail.com
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
!
! Note, one could achieve the same thing using an allgather
! or allreduce, however this example allows the order of communication be
! to be varied, and allows for communication between a subset of processes.
!
! Compile: mpif90 sizes.f90 -o s.exe
! Run: mpirun -np <NP> s.exe
! Works for any number of processes
!
program sizes
  use mpi
  implicit none

  integer :: rank, np, ierr, iprocess, send_id, rec_id, iround, send_index
  integer, allocatable :: collated_data(:), rec_ids(:), initial_rec_ids(:)
  
  call MPI_INIT(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

  ! Initial on all processes, with some arb value >> rank
  allocate(collated_data(np), source=101)

  ! Assign from same process
  do iprocess = 0, np - 1
     if(rank==iprocess) collated_data(iprocess+1) = rank
  enddo
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  ! Start receivers by permuting the sender indices by 1
  initial_rec_ids = modulo([(i, i = 0, np-1)] + 1, np)

  ! Round comm pattern
  ! Permute the receiver id by 1, per communication round
  do iround = 1, np - 1
     rec_ids = modulo(initial_rec_ids + iround - 1, np)

     ! This loop defines sender ids, [0, 1, 2, 3, ..., np - 1]
     do send_id = 0, np - 1
        send_index = send_id + 1
        rec_id = rec_ids(send_id+1)
        ! Note, one could pass any contiguous block of data i.e.
        ! collated_data(i:j)
        call mpi_send_receive(MPI_COMM_WORLD, rank, send_id, rec_id, collated_data(send_index))
     enddo
     
  enddo
  
  ! All ranks have collated_data from one another
  ! Expect each rank to print  collated_data = [0, 1, 2, ..., np - 1]
  write(*, *) rank, collated_data
  
  call MPI_FINALIZE(ierr)

contains

  ! Send and receive from the same array, with different elements specifed on different processes
  !
  ! Array is declared for np processes
  ! Array only contains local information => Send that information to the same index of
  ! array, on all other processes (np - 1), in a ring pattern.
  subroutine mpi_send_receive(comm, rank, send_id, rec_id, data)
    integer, intent(in) :: comm, rank
    integer, intent(in) :: send_id, rec_id
    integer, contiguous, intent(inout) :: data(..)
    integer, parameter :: tag = 101
    
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)
    
    if (rank == send_id) then
       call MPI_Send(data, size(data), MPI_INTEGER, rec_id, tag, comm, ierr)
    elseif (rank == rec_id)then
       call MPI_Recv(data, size(data), MPI_INTEGER, send_id, tag, comm, status, ierr)
    endif
  end subroutine mpi_send_receive

end program sizes
