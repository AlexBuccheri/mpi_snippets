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

!> Test shared memory model for MPI
!! Compilation on my mac
!! mpif90 -pedantic -fbacktrace -fbounds-check -Wl,-syslibroot,/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk shared_mem_mpi.f90 -o mpi_shared
!! mpirun -np 2 ./mpi_shared
program shared_mem_mpi
  use mpi_f08
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: intra_root = 0 !< Root process for intra comm. I.e. the Origin process. This is the only one to call RMA routines
  type(MPI_Comm) :: intra_comm, inter_comm              !< Comms
  integer        :: base_rank, intra_rank, inter_rank   !< Ranks per comm
  integer        :: ierr                                !< Error integer
  integer        :: intra_np

  
  ! MPI Shared memory specifics
  integer(MPI_ADDRESS_KIND) :: window_size              !< Size of the memory window, in bytes
  type(c_ptr)               :: ptr                      !< ADD ME
  type(MPI_Win)             :: window                   !< Shared memory window handle
  integer                   :: disp_unit                !< Local unit size for displacements, in bytes

  ! Data
  integer          ::  number_of_elements, n_elements_local, offset, i, j, intra_size, inter_size
  integer, pointer ::  array(:)
  
  
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, base_rank, ierr)

  ! Set up some data
  number_of_elements = 48
  
  ! Create communicator for intranode communication i.e. into shared memory islands
  ! - This is used for MPI-3 shared memory windows. Number of processes is determined by the hardware
  call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, base_rank, MPI_INFO_NULL, intra_comm, ierr)
  call MPI_Comm_size(intra_comm, intra_size, ierr)
  call mpi_comm_rank(intra_comm, intra_rank, ierr)

  ! Create communicator for internode communication
  ! - We have one communicator for each local rank
  call MPI_Comm_split(MPI_COMM_WORLD, intra_rank, base_rank, inter_comm, ierr)

  call mpi_comm_rank(inter_comm, inter_rank, ierr)
  call MPI_Comm_size(inter_comm, inter_size, ierr)
  write(*,*) 'Process', base_rank, 'has intra_rank', intra_rank, &
             'intra_size', intra_size, 'inter_rank', inter_rank, 'inter_size', inter_size
  
  ! Allocate only on rank 0 of each intra_comm, else one will get the unnecessary copies we
  ! are trying to avoid
  if (intra_rank == intra_root) then
     window_size = number_of_elements * sizeof(int32)
  else
     window_size = 0_MPI_ADDRESS_KIND
  end if
    
  ! Defined a shared memory window on each island. This sets ptr and the window handle,
  ! given the specified window size on intra_comm.
  !
  ! This procedure is collective over the group of intra_comm.
  ! - On each MPI process, it allocates memory of at least size bytes that is shared among all MPI processes in intra_comm,
  !   and returns a pointer to the locally allocated segment in ptr that can be used for load/store accesses on
  !   the calling MPI process.
  ! Note: MPI_INFO_NULL could be replaced with a pointer to window info
  call MPI_WIN_ALLOCATE_SHARED(window_size, int(sizeof(int32)), MPI_INFO_NULL, intra_comm, ptr, window, ierr)

  ! Get the pointer on all ranks
  call MPI_Win_shared_query(window, intra_root, window_size, disp_unit, ptr, ierr)
  
  ! Assign C pointer to fortran pointer
  call c_f_pointer(ptr, array, [number_of_elements])
  
  n_elements_local = number_of_elements / intra_size
  offset = intra_rank * n_elements_local
  write(*, *) 'Intra rank', intra_rank, 'runs from ', offset + 1, 'to', offset + n_elements_local
  call MPI_Win_fence(0, window, ierr)

  ! Allow the calling process to access the memory regions of all other processes in that
  ! group without interference from other processes
  ! Allows one-sided communication - the target process doesn't need to explicitly participate
  ! where the target is the rank that calls shared_query
  ! call MPI_Win_lock_all(MPI_MODE_NOCHECK, window, ierr)
  
  do i = 1, n_elements_local
     j = i + offset
     array(j) = j
  enddo

   ! Synchronize before accessing shared memory
  call MPI_Win_fence(0, window, ierr)

  ! One can print the whole array from any process in intra_rank
  if (intra_rank == 5) then
     write(*, *) 'Access shared memory from (see specified process)'
     do j = 1, number_of_elements
        write(*, *) array(j) 
     enddo
  end if
  
  !call MPI_Win_unlock_all(window, ierr)

  ! Clean up
  call MPI_Win_fence(0, window, ierr)
  call MPI_Win_free(window, ierr)
  nullify(array)

  call mpi_finalize(ierr)
  
end program shared_mem_mpi
