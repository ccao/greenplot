MODULE input
  !
  USE constants
  !
  implicit none
  !
  real(dp) :: eta, ef
  real(dp), allocatable :: emesh(:), xk(:)
  real(dp), allocatable :: kvec(:, :)
  !
  integer ne, nk
  !
  logical isbulk
  !
  character(len=80) :: seed
  !
 CONTAINS
  !
 SUBROUTINE read_input()
  !
  USE constants
  USE para
  !
  implicit none
  !
  real(dp), dimension(3, 3) :: alat, bvec
  !
  integer nseg, nkseg
  !
  real(dp), dimension(3) :: t1, t2, tt
  real(dp) omega
  integer ii, jj
  !
  if (inode.eq.0) then
    !
    open(unit=fin, file='greenplot.inp')
    !
    read(fin, *) seed
    !
    read(fin, *) t2(1), t1(1:3)
    !
    write(*, '(A,1F9.4)') ' System Fermi level : Ef=', t2(1)
    !
  endif
  !
  CALL para_sync(t1, 3)
  CALL para_sync(t2, 3)
  !
  ef=t2(1)
  !
  ne=nint((t1(2)-t1(1))/t1(3))+1
  !
  allocate(emesh(ne))
  !
  do ii=1, ne
    !
    emesh(ii)=t1(1)+(ii-1)*t1(3)
    !
  enddo
  !
  eta=t1(3)*2.d0
  !
  if (inode.eq.0) then
    !
    do ii=1, 3
      read(fin, *) alat(ii, :)
    enddo
    !
  endif
  !
  CALL para_sync(alat, 3, 3)
  !
  CALL cross_prod(bvec(1, :), alat(2, :), alat(3, :))
  CALL cross_prod(bvec(2, :), alat(3, :), alat(1, :))
  CALL cross_prod(bvec(3, :), alat(1, :), alat(2, :))
  !
  omega=sum(bvec(1, :)*alat(1, :))
  bvec(:, :)=bvec(:, :)/omega
  !
  if (inode.eq.0) then
    !
    read(fin, *) ii
    read(fin, *) nseg, nkseg
    !
    t1(1)=ii
    t1(2)=nseg
    t1(3)=nkseg
    !
  endif
  !
  CALL para_sync(t1, 3)
  !
  ii=nint(t1(1))
  nseg=nint(t1(2))
  nkseg=nint(t1(3))
  !
  if (ii.eq.3) isbulk=.true.
  !
  nk=nseg*nkseg
  !
  allocate(kvec(nk, 3), xk(nk))
  !
  if (inode.eq.0) then
    !
    do ii=0, nseg-1
      !
      read(fin, *) t1(:), t2(:)
      !
      do jj=1, 3
        tt(jj)=sum((t1(:)-t2(:))*bvec(:, jj))
      enddo
      !
      omega=dsqrt(sum(tt(:)**2))
      !
      do jj=1, nkseg
        !
        kvec(ii*nkseg+jj, :)=((jj-1)*t2(:)+(nkseg-jj)*t1(:))/(nkseg-1)
        !
        if (jj.eq.1) then
          if (ii.eq.0) then
            xk(ii*nkseg+jj)=0.d0
          else
            xk(ii*nkseg+jj)=xk(ii*nkseg+jj-1)
          endif
        else
          xk(ii*nkseg+jj)=xk(ii*nkseg+jj-1)+omega/(nkseg-1)
        endif
        !
      enddo  ! jj
      !
    enddo  ! ii
    !
    close(unit=fin)
    !
  endif   ! inode.eq.0
  !
  CALL para_sync(xk, nk)
  CALL para_sync(kvec, nk, 3)
  !
 END SUBROUTINE
  !
 SUBROUTINE cross_prod(v1, v2, v3)
  !
  use constants
  !
  implicit none
  !
  real(dp), dimension(3) :: v1, v2, v3
  !
  v1(1)=v2(2)*v3(3)-v2(3)*v3(2)
  v1(2)=v2(3)*v3(1)-v2(1)*v3(3)
  v1(3)=v2(1)*v3(2)-v2(2)*v3(1)
  !
 END SUBROUTINE
  !
END MODULE
