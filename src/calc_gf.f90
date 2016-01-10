!include 'lapack.f90'

SUBROUTINE calc_gf(gf, h, ene, eta)
  !
  USE lapack95,   only : getrf, getri, hetrf, hetri
  USE constants
  USE hamk,  only : nbnd
  !
  implicit none
  !
  complex(dp), dimension(nbnd, nbnd) :: gf
  complex(dp), dimension(nbnd, nbnd) :: h
  real(dp) :: ene, eta
  integer, allocatable :: ipiv (:)
  integer ii
  !
  gf(:, :)=-h(:, :)
  !
  allocate(ipiv(nbnd))
  !
  if (eta.gt.eps12) then
    !
    do ii=1, nbnd
      gf(ii, ii)=gf(ii, ii)+ene+eta*cmplx_i
    enddo
    !
    CALL getrf(gf, ipiv)
    CALL getri(gf, ipiv)
    !
  else
    !
    do ii=1, nbnd
      gf(ii, ii)=gf(ii, ii)+ene
    enddo
    !
    CALL hetrf(gf, 'U', ipiv)
    CALL hetri(gf, ipiv)
    !
  endif
  !
  deallocate(ipiv)
  !
END SUBROUTINE
