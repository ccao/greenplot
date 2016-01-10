SUBROUTINE calc_bulkgf(gf, ene, eta)
  !
  USE constants
  USE hamk
  !
  implicit none
  !
  complex(dp), dimension(nbnd, nbnd) :: gf
  real(dp) :: ene, eta
  complex(dp), allocatable :: work (:)
  integer, allocatable :: ipiv(:)
  integer info, ii
  !
  allocate(work(nbnd), ipiv(nbnd))
  !
  gf(:, :)=-ham00(:, :)
  !
  do ii=1, nbnd
    gf(ii, ii)=gf(ii, ii)+ene+eta*cmplx_i
  enddo
  !
  CALL zgetrf(nbnd, nbnd, gf, nbnd, ipiv, info)
  !
  if (info.ne.0) then
    write(*, *) '!!! FATAL ERROR : zgetrf failed!'
    stop
  endif
  !
  CALL zgetri(nbnd, gf, nbnd, ipiv, work, nbnd, info)
  !
  if (info.ne.0) then
    write(*, *) '!!! FATAL ERROR: zgetri failed!'
    stop
  endif
  !
  deallocate(work, ipiv)
  !
END SUBROUTINE
