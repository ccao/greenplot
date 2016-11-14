MODULE ucscinp
  !
  use constants
  !
  implicit none
  !
  real(dp),dimension(3,3) :: sc
  !  Definition of supercells (i.e., A1'=sc(1,1)*A1+sc(1,2)*A2+sc(1,3)*A3, etc)
  integer nuc
  !  Number of unit cells, = det(sc)
  real(dp),allocatable :: xorb(:, :)
  !  The lattice coordinate of orbitals in SC (in unit of unit cell lattices)
  !  e.g., in 1x2x1 SC, orbitals in the 2nd UC are xorb(:, iorb)=(0,1,0)
  integer :: nnx, nny, nnz
  !  The interaction size (in the SC)
  integer :: nnorb      !  number of orbitals in SC
  integer :: nnrpt      !  number of rvecs in SC
  !
END MODULE

PROGRAM uc2sc
  !
  ! This program reads wannier90 hamiltonian of a unit cell
  !   and writes out a super cell hamiltonian in the same format
  !
  !  Sample input
  !*******************
  ! wannier90_hr.dat
  ! 1 -1  0    A1'
  ! 1  1  0    A2'
  ! 0  0  1    A3'
  ! 12 12 12   new wannier90_hr.dat size
  !  1         Number of orbitals (in UC)
  ! 0.5 0.0 0.0    orb 1 pos  (in UC)
  !
  use constants
  use wanndata
  use ucscinp
  !
  implicit none
  !
  integer :: nx, ny, nz  ! in UC
  integer,allocatable :: rveclist(:, :, :)
  !  The (reverse) list of rvecs
  !
  ! Local variables
  !
  integer ix, iy, iz    ! in SC
  integer ib, jb, ii
  !
  character(len=80) seed
  !
  CALL read_input()
  !
  CALL read_ham(seed)
  !
  nx=nint(maxval(1,:))
  ny=nint(maxval(2,:))
  nz=nint(maxval(3,:))
  !
  allocate(rveclist(-nx:nx,-ny:ny,-nz:nz))
  do ii=1, nrpt
    ix=nint(rvec(1,ii))
    iy=nint(rvec(2,ii))
    iz=nint(rvec(3,ii))
    rveclist(ix,iy,iz)=ii
  enddo
  !
  do ix=-nnx,nnx
    do iy=-nny,nny
      do iz=-nnz,nnz
   
  !
END PROGRAM

SUBROUTINE read_input()
  !
  use constants
  use ucscinp
  !
  implicit none
  !
  integer norb
  real(dp),allocatable :: xx(:, :)
  ! Coordinates of orbitals (in UC)
  real(dp),allocatable :: xxx(:, :)
  ! Coordinates of orbitals (in SC)
  integer ii, jj, kk
  real(dp),dimension(3) :: tt
  !
  open(unit=fin, file='supercell.def')
  !
  read(fin, *) seed
  do ii=1, 3
    read(fin, *) sc(:, ii)
  enddo
  read(fin, *) nnx, nny, nnz
  read(fin, *) norb
  allocate(xx(3, norb))
  do ii=1, norb
    read(fin, *) xx(:, ii)
  enddo
  !
  close(unit=fin)
  !
  nuc=det(sc)
  !
  nnorb=norb*nuc
  allocate(xxx(3,nnorb))
  !
  CALL fillsc(xxx, xx, sc, norb, nuc)
  !
  do ii=1, nuc
    do jj=1, norb
      (tt(kk)=sum(xxx(:,(ii-1)*norb+jj)*sc(kk,:)),kk=1,3)
      xxorb(:,(ii-1)*norb+jj)=nint(tt(:)-xx(:,jj))
    enddo
  enddo
  !
  deallocate(xx)
  deallocate(xxx)
  !
END SUBROUTINE

SUBROUTINE fillsc(xxx, xx, sc, norb, nuc)
  !
  use constants
  !
  implicit none
  !
  real(dp),dimension(3,norb*nuc) :: xxx
  real(dp),dimension(3,norb) :: xx
  real(dp),dimension(3,3) :: sc
  integer norb, nuc
  !
  real(dp) :: uc(3,3)
  ! UC lattice vectors in SC lattices (e.g., A1=uc(1,1)A1'+uc(1,2)A2'+uc(1,3)A3',
  !   etc)
  integer,dimension(3) :: ipiv
  integer ix, iy, iz, ii, jj, nn
  real(dp),allocatable :: lat(:,:)
  real(dp),dimension(3) :: t,tt
  logical isold
  !
  uc(:,:)=sc(:,:)
  CALL getrf(uc, ipiv)
  CALL getri(uc, ipiv)
  !
  allocate(lat(3,nuc))
  nn=0
  do ix=0, ceiling(sqrt(sum(sc(:,1)**2)))-1
   do iy=0, ceiling(sqrt(sum(sc(:,2)**2)))-1
    do iz=0, ceiling(sqrt(sum(sc(:,3)**2)))-1
      if (nn==0)
        nn=nn+1
        lat(:, nn)=0.d0
      else
        t(1)=ix
        t(2)=iy
        t(3)=iz
        do ii=1,3
          tt(ii)=unitize(sum(uc(jj,:)*t(:)))
        enddo
        isold=.false.
        do ii=1,nn
          if (sum((tt(:)-lat(:,ii))**2).lt.eps4) isold=.true.
        enddo
        if (.not. isold) then
          nn=nn+1
          lat(:, nn)=tt(:)
        endif
      endif
    enddo ! iz
   enddo  ! iy
  enddo   ! ix
  !
  if (nn.lt.nuc) then
    write(*,*) '!!! FATAL! Incorrect!'
    stop
  endif
  !
  do nn=1, nuc
    do ii=1, norb
      do jj=1, 3
        xxx(jj, (nn-1)*norb+ii)=unitize(sum(uc(jj,:)*xx(:,ii))+lat(jj,nn))
      enddo
    enddo
  enddo
  !
  deallocate(lat)
  !
END SUBROUTINE

FUNCTION unitize(x)
  use constants
  implicit none
  !
  real(dp) x, unitize
  !
  unitize=x-floor(x);
END FUNCTION

SUBROUTINE cross_prod(c, a, b)
  !
  use constants
  !
  implicit none
  !
  real(dp),dimension(3) :: c, a, b
  !
  c(1)=a(2)*b(3)-a(3)*b(2)
  c(2)=a(3)*b(1)-a(1)*b(3)
  c(3)=a(1)*b(2)-a(2)*b(1)
  !
END SUBROUTINE

FUNCTION det(m)
  !
  use constants
  !
  implicit none
  !
  real(dp) det
  real(dp), dimension(3,3) :: m
  !
  det=m(1,1)*m(2,2)*m(3,3)+m(1,2)*m(2,3)*m(3,1)+m(1,3)*m(2,1)*m(3,2)- &
      m(1,3)*m(2,2)*m(3,1)-m(1,2)*m(2,1)*m(3,3)-m(1,1)*m(2,3)*m(3,2)
  !
  if (det<0) det=-det
  !
END FUNCTION
