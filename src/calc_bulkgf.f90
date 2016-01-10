SUBROUTINE calc_bulkgf(gf, ene, eta)
  !
  USE constants
  USE hamk
  !
  implicit none
  !
  complex(dp), dimension(nbnd, nbnd) :: gf
  real(dp) :: ene, eta
  !
  CALL calc_gf(gf, ham00, ene, eta)
  !
END SUBROUTINE
