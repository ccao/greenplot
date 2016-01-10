MODULE para
  !
#if defined __MPI
  use mpi
#endif
  !
  implicit none
  !
  integer inode, nnode
  integer first_ene, last_ene
  !
INTERFACE para_sync
  MODULE PROCEDURE para_sync_int0, para_sync_real0, para_sync_real1, para_sync_real2, para_sync_cmplx3
END INTERFACE
  !
CONTAINS
  !
SUBROUTINE init_para()
  !
  use constants
  !
  implicit none
  !
  integer ierr
  !
#if defined __MPI
  CALL mpi_init(ierr)
  !
  CALL mpi_comm_rank(mpi_comm_world, inode, ierr)
  CALL mpi_comm_size(mpi_comm_world, nnode, ierr)
  !
  if (inode.eq.0) write(stdout, *) "# Greenplot running on ", nnode, " nodes..."
  !
#else
  inode=0
  nnode=1
  write(stdout, *) "# Greenplot serial ..."
#endif
  !
END SUBROUTINE
  !
SUBROUTINE split_emesh(ne)
  !
  implicit none
  !
  integer ne
  integer part_n, res
  integer ii
  !
#if defined __MPI
  !
  res=mod(ne, nnode)
  part_n=(ne-res)/nnode
  !
  if (inode<res) then
    first_ene=inode*(part_n+1)+1
    last_ene=(inode+1)*(part_n+1)
  else
    first_ene=res*(part_n+1)+(inode-res)*part_n+1
    last_ene=res*(part_n+1)+(inode-res+1)*part_n
  endif
  !
#else
  !
  first_ene=1
  last_ene=ne
  !
#endif
  !
END SUBROUTINE
  !
SUBROUTINE para_sync_int0(dat)
  !
  implicit none
  !
  integer :: dat
  !
  integer ierr
  !
#if defined __MPI
  CALL mpi_bcast(dat, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_sync_real0(dat)
  !
  use constants, only: dp
  !
  implicit none
  !
  real(dp) :: dat
  !
  integer ierr
  !
#if defined __MPI
  CALL mpi_bcast(dat, 1, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_sync_real1(dat, dat_size)
  !
  use constants, only : dp
  !
  implicit none
  !
  integer :: dat_size
  real(dp) :: dat(1:dat_size)
  !
  integer ierr
  !
#if defined __MPI
  CALL mpi_bcast(dat, dat_size, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_sync_real2(dat, size1, size2)
  !
  use constants, only: dp
  !
  implicit none
  !
  integer :: size1, size2
  real(dp) :: dat(1:size1, 1:size2)
  !
  integer ierr, dat_size
  !
#if defined __MPI
  dat_size=size1*size2
  CALL mpi_bcast(dat, dat_size, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_sync_cmplx3(dat, size1, size2, size3)
  !
  use constants, only : dp
  !
  implicit none
  !
  integer :: size1, size2, size3
  complex(dp) :: dat(1:size1, 1:size2, 1:size3)
  !
  integer ierr, dat_size
  !
#if defined __MPI
  dat_size=size1*size2*size3
  CALL mpi_bcast(dat, dat_size, MPI_DOUBLE_COMPLEX, 0, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE
  !
SUBROUTINE para_merge(dat, dat_size)
  !
  use constants, only :dp
  !
  implicit none
  !
  integer :: dat_size
  real(dp) :: dat(1:dat_size)
  real(dp), allocatable :: buf(:)
  !
  integer ierr
  !
#if defined __MPI
  allocate(buf(1:dat_size))
  CALL mpi_allreduce(dat, buf, dat_size, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
  dat(:)=buf(:)
  deallocate(buf)
#endif
  !
END SUBROUTINE
  !
SUBROUTINE finalize_para
  !
  implicit none
  !
#if defined __MPI
  integer ierr
  CALL mpi_finalize(ierr)
#endif
  !
END SUBROUTINE

END MODULE
