MODULE hamk
  !
  USE constants
  !
  implicit none
  !
  real(dp),allocatable :: surfwt(:)
  complex(dp),allocatable :: ham00(:, :), ham01(:, :)
  !
  integer nbnd
  !
 CONTAINS
  !
  SUBROUTINE init_hamk()
    !
    USE constants,    ONLY : dp
    USE wanndata,     ONLY : norb, rvec, nrpt
    USE input,        ONLY : mode
    !
    implicit none
    !
    integer nl
    !
    if (mode.eq.0) then
      ! Bulk calculation, no need to worry about other stuff...
      nbnd=norb
      allocate(ham00(nbnd, nbnd))
    else
      ! Surface calculation
      !   The original wannier90 Hamiltonian contains interactions in all
      !   three dimensions. This Hamiltonian is defined with H(R), where
      !   R=(nx, ny, nz) denoting interactions between orbitals in the home
      !   cell (0,0,0) and the unit cell in (nx, ny, nz). R has a cutoff as
      !   defined by nkx, nky, nkz when producing the Hamiltonian, meaning
      !   that for all Rs that | nz | > nnz the interactions will be 0. Thus,
      !   We take nnz+1 layers to construct the 2D Hamiltonian to avoid loss
      !   of interactions.
      !
      nl=nint(maxval(abs(rvec(3, :))))+1
      nbnd=norb*nl
      allocate(ham01(nbnd, nbnd), ham00(nbnd, nbnd))
      !
      ! After that, we have to recalculate the weight of real space lattice
      !   vectors, so that the Hamiltonian in 2D k-space can be successfully
      !   constructed.
      !
      allocate(surfwt(nrpt))
      CALL find_2dws()
      !
    endif ! mode.eq.0
    !
  END SUBROUTINE

  SUBROUTINE find_2dws()
    ! FIXME: did not consider the possibility that a3 is not parallel to a1xa2
    use constants,  only : dp, eps7, eps5
    use input,      only : alat
    use wanndata,   only : rvec, nrpt
    !
    implicit none
    !
    integer nx, ny, ix, iy, i1, i2
    integer ii
    integer i, j
    integer nn, ncount
    integer, dimension(1:2) :: rtmp
    real(dp), dimension(1:25) :: dd
    real(dp) dd_min, sum_tot
    !
    surfwt(:)=0.d0
    nx=nint(maxval(rvec(1, :))-minval(rvec(1, :)))
    ny=nint(maxval(rvec(2, :))-minval(rvec(2, :)))
    !
    ncount=0
    sum_tot=0.d0
    do ix=-nx, nx
      do iy=-ny, ny
        ii=0
        do i1=-2, 2
          do i2=-2, 2
            ii=ii+1
            rtmp(1)=ix-i1*nx
            rtmp(2)=iy-i2*ny
            dd(ii)=0.d0
            do i=1, 3
              dd(ii)=dd(ii)+sum(rtmp(:)*alat(1:2,i))**2
            enddo
          enddo
        enddo
        dd_min=minval(dd)
        if(abs(dd(13)-dd_min).lt.eps7) then
          !
          nn=0
          do ii=1,25
            if (abs(dd(ii)-dd_min).lt.eps7) nn=nn+1
          enddo
          !
          sum_tot=sum_tot+1.d0/nn
          !
          do ii=1, nrpt
            rtmp(1)=rvec(1, ii)-ix
            rtmp(2)=rvec(2, ii)-iy
            if (sum(abs(rtmp(:))).eq.0) then
              surfwt(ii)=nn
            endif
          enddo
          !
        endif
      enddo
    enddo
    !
    if (abs(sum_tot-nx*ny).gt.eps5) then
      write(stdout, '(A,1F16.9)') '!!! ERROR: sum of weight is ', sum_tot
    endif
    !
    do ii=1, nrpt
      if (surfwt(ii).eq.0) then
        write(stdout, '(A,3I5,A)') '!!! ERROR: weight of rvec:', nint(rvec(:, ii)), ' is 0'
      endif
    enddo
    !
  END SUBROUTINE
  !
  SUBROUTINE finalize_hamk()
    !
    implicit none
    !
    if (allocated(ham00)) deallocate(ham00)
    !
    if (allocated(ham01)) deallocate(ham01)
    !
  END SUBROUTINE
  !
  SUBROUTINE calc_hamk( kv )
    !
    USE constants
    USE wanndata
    USE input,        ONLY : mode
    !
    implicit none
    !
    real(dp), dimension(3) :: kv
    !
    complex(dp) coeff
    complex(dp), allocatable :: htmp(:, :, :)
    !
    integer ir, nlayer, nl
    integer il, jl
    !
    ham00(:, :) = 0.d0
    !
    if (mode.eq.0) then
      !
      do ir=1, nrpt
        !
        coeff=exp(-cmplx_i*twopi*sum(kv(:)*rvec(:, ir)))
        ham00(:, :)=ham00(:, :)+coeff*ham(:, :, ir)/weight(ir)
        !
      enddo
      !
    else
      !
      ham01(:, :)=0.d0
      kv(3)=0.d0
      !
      nl=nint(maxval(abs(rvec(3, :))))
      nlayer=2*nl+1
      allocate(htmp(norb, norb, nlayer))
      !
      htmp(:, :, :)=0.d0
      !
      do ir=1, nrpt
        !
        il=nint(rvec(3, ir))+nl+1
        coeff=exp(-cmplx_i*twopi*sum(kv(:)*rvec(:, ir)))
        htmp(:, :, il) = htmp(:, :, il)+coeff*ham(:, :, ir)/surfwt(ir)
        !
      enddo
      !
      do il=1, nl+1
        !
        do jl=1, nl+1
          !
          ham00((il-1)*norb+1:il*norb, (jl-1)*norb+1:jl*norb) = htmp(:, :, (il-jl)+nl+1)
          !
          if (il-jl>0) then
            !
            ham01((il-1)*norb+1:il*norb, (jl-1)*norb+1:jl*norb) = htmp(:, :, (il-jl))
            !
          endif
          !
        enddo ! jl
        !
      enddo   ! il
      !
      deallocate(htmp)
      !
    endif    ! isbulk
    !
  END SUBROUTINE
  !
END MODULE
