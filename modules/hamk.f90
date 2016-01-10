MODULE hamk
  !
  USE constants
  !
  implicit none
  !
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
    USE input,        ONLY : isbulk
    !
    implicit none
    !
    integer nl
    !
    if (isbulk) then
      nbnd=norb
      allocate(ham00(nbnd, nbnd))
    else
      nl=nint(maxval(rvec(3, :)))+1
      nbnd=norb*nl
      allocate(ham01(nbnd, nbnd), ham00(nbnd, nbnd))
    endif
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
    USE input,        ONLY : isbulk
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
    if (isbulk) then
      !
      do ir=1, nrpt
        !
        coeff=exp(-cmplx_i*twopi*sum(kv(:)*rvec(:, ir)))/weight(ir)
        ham00(:, :)=ham00(:, :)+coeff*ham(:, :, ir)
        !
      enddo
      !
    else
      !
      ham01(:, :)=0.d0
      kv(3)=0.d0
      !
      nl=nint(maxval(rvec(3, :)))
      nlayer=2*nl+1
      allocate(htmp(norb, norb, nlayer))
      !
      htmp(:, :, :)=0.d0
      !
      do ir=1, nrpt
        !
        il=nint(rvec(3, ir))+nl+1
        coeff=exp(-cmplx_i*twopi*sum(kv(:)*rvec(:, ir)))/weight(ir)
        htmp(:, :, il) = htmp(:, :, il)+coeff*ham(:, :, ir)
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
