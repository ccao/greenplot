MODULE input
  !
  ! Example:
  !wannier90    ! seed
  ! 8.26807433 -0.4  0.4  0.001   ! ef emin emax de
  !     5.1781764    -1.8482837     0.0000000  ! a1
  !     5.1781764     1.8482837     0.0000000  ! a2
  !    -4.2580695     0.0000000     7.2229161  ! a3
  ! 1  1     ! mode  dimension
  ! 3  500   ! par1 par2
  ! G 0.00  0.00  0.00 M 0.50  0.00  0.00
  ! M 0.50  0.00  0.00 X 0.50  0.50  0.00
  ! X 0.50  0.50  0.00 G 0.00  0.00  0.00
  !
  ! Calculation modes:
  !    0 : bulk calculation
  !    1 : surface band calculation
  ! dimension:
  !    0 : single k-point calculation
  !    1 : line mode
  !      par1 is number of high symmetry line segs
  !      par2 is number of k-points in each seg
  !    2 : surface iso energy state
  !      par1xpar2 is the k-mesh size
  !
  !     For surface calculations, assumes [001] direction
  !      or more accurately, the surface is defined by
  !      A1 x A2, transform Hamiltonian if necessary
  !
  USE constants
  !
  implicit none
  !
  real(dp), dimension(3,3) :: alat, bvec
  ! alat is the real-space lattice vector
  ! bvec is the reciprocal-space lattice vector
  real(dp) :: eta, ef
  ! eta is the infinitesmal imaginary part in Green's function
  ! ef is the Fermi level
  real(dp), allocatable :: emesh(:), xk(:)
  ! emesh is the energy mesh
  ! xk is for simple plot
  real(dp), allocatable :: kvec(:, :)
  ! kvec(3, nk) are the actual k-points to be calculated
  integer ne, nk, nkx, nky
  integer mode, dmsn
  ! Calculation modes and dimension
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
  integer par1, par2
  !
  real(dp), dimension(3) :: t1, t2, tt
  real(dp) omega
  integer ii, jj
  character(len=10) label1, label2
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
    read(fin, *) par1,par2   ! mode and dimension
    t1(1)=par1
    t1(2)=par2
    read(fin, *) par1,par2
    t2(1)=par1
    t2(2)=par2
    !
  endif
  !
  CALL para_sync(t1, 3)
  CALL para_sync(t2, 3)
  !
  mode=nint(t1(1))
  dmsn=nint(t1(2))
  nkx=nint(t2(1))
  nky=nint(t2(2))
  !
  nk=nkx*nky
  !
  allocate(kvec(3, nk), xk(nk))
  !
  if (inode.eq.0) then
    !
    if (dmsn.eq.0) then                             ! Single point calculation
      !
      write(stdout, '(A,1I,A)') '  Single point calculation using ', ne, ' energies.'
      read(fin, *) kvec(:,1)
      !
    elseif (dmsn.eq.1) then                         ! Line mode is default for bulk calculation
      !
      write(stdout, '(A)') '  Line mode calculation'
      open(unit=fout, file='klabels.dat')           ! Create this file for easy plot.
      !
      do ii=0, par1-1
        !
        read(fin, *) label1, t1(:), label2, t2(:)   ! The labels of this segment
        do jj=1, 3
          tt(jj)=sum((t1(:)-t2(:))*bvec(:, jj))
        enddo
        omega=dsqrt(sum(tt(:)**2))                  ! Length of this segment
        !
        do jj=1, par2                               ! Interpolates the k-points in between
          !
          kvec(:, ii*par2+jj)=((jj-1)*t2(:)+(par2-jj)*t1(:))/(par2-1)
          if (jj.eq.1) then                         !   as well as the xk position in plot
            if (ii.eq.0) then
              xk(ii*par2+jj)=0.d0
            else
              xk(ii*par2+jj)=xk(ii*par2+jj-1)
            endif
          else
            xk(ii*par2+jj)=xk(ii*par2+jj-1)+omega/(par2-1)
          endif
          !
        enddo  ! jj
        !
        write(fout, '(1F16.9,1X,A,1X,A)') xk(ii*par2+1), label1, label2
        !
      enddo  ! ii
      !
      close(unit=fout)
      close(unit=fin)
    elseif (dmsn.eq.2) then                         ! Surface mode
      !  Create the 2d-mesh
      write(stdout, '(A)') '  2D mode calculation'
      write(stdout, '(A,I5,A,I5)') '   mesh size: ', nkx, 'x', nky
      do ii=1,nkx
        do jj=1,nky
          !
          kvec(:, (ii-1)*nky+jj)=(/ 1.d0*(ii-1)/nkx, 1.d0*(jj-1)/nky, 0.d0 /)
          !
        enddo
      enddo
      !
    endif
    !
  endif   ! inode.eq.0
  !
  CALL para_sync(xk, nk)
  CALL para_sync(kvec, 3, nk)
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
