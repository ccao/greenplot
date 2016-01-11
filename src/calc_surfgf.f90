!include 'blas.f90'

SUBROUTINE calc_surfgf(gf, ene, eta)
  !
  USE blas95,   only : gemm
  USE constants
  USE hamk
  !
  implicit none
  !
  complex(dp), dimension(nbnd, nbnd) :: gf
  real(dp) :: ene, eta
  !
  complex(dp), allocatable :: a0(:, :), b0(:, :)
  complex(dp), allocatable :: a1(:, :), b1(:, :)
  complex(dp), allocatable :: h0(:, :), g0(:, :)
  complex(dp), allocatable :: hs(:, :)
  complex(dp), allocatable :: ga(:, :), gb(:, :)
  !
  integer iter, ii, jj
  integer info
  !
  allocate(a0(nbnd, nbnd), a1(nbnd, nbnd))
  allocate(b0(nbnd, nbnd), b1(nbnd, nbnd))
  allocate(h0(nbnd, nbnd), g0(nbnd, nbnd))
  allocate(hs(nbnd, nbnd))
  allocate(ga(nbnd, nbnd), gb(nbnd, nbnd))
  !
  hs(:, :)=ham00(:, :) ! iteration 0 : h00
  h0(:, :)=ham00(:, :) ! iteration 0 :
  !
  CALL calc_gf(g0, h0, ene, eps3) ! calculates the green's function for iteration 0
  !
  a0(:, :)=ham01(:, :) ! iteration 0 : a0 and b0
  !
  do ii=1, nbnd
    do jj=1, nbnd
      b0(ii, jj)=conjg(a0(jj, ii))
    enddo
  enddo
  !
  do iter=1, maxiter
    !
    ! Calculates the temporary matrices
    !  ga=g_i*a_i, gb=g_i*b_i
    !
    CALL gemm(g0, a0, ga)
    CALL gemm(g0, b0, gb)
    !
    ! Calculates the new Hamiltonian
    !   h_{i+1)=h_i+a_i*g_i*cmplx_0_i+cmplx_0_i*g_i*cmplx_1_i
    !
    CALL gemm(a0, gb, h0, 'N', 'N', cmplx_1, cmplx_1)
    CALL gemm(b0, ga, h0, 'N', 'N', cmplx_1, cmplx_1)
    !
    ! Calculates the new Surface hamiltonian
    !  hs_{i+1}=hs_i+cmplx_1_i*g_i*cmplx_0_i
    !
    CALL gemm(a0, gb, hs, 'N', 'N', cmplx_1, cmplx_1)
    !
    ! Calculates the new cmplx_1, cmplx_0
    !  cmplx_1_{i+1}=cmplx_1_i*g_i*cmplx_1_i, cmplx_0_{i+1}=cmplx_0_i*g_i*cmplx_0_i
    !
    CALL gemm(a0, ga, a1)
    CALL gemm(b0, gb, b1)
    !
    if ((maxval(abs(a1)).le.eps6).and.(maxval(abs(b1)).le.eps6)) exit ! do iter
    !
    a0(:,:)=a1(:,:)
    b0(:,:)=b1(:,:)
    !
    CALL calc_gf(g0, h0, ene, eps3)
    !
  enddo
  !
  CALL gemm(a0, gb, hs, 'N', 'N', cmplx_1, cmplx_1)
  !
  CALL calc_gf(gf, hs, ene, eta)
  !
  deallocate(a0, a1, b0, b1, h0, g0, hs, ga, gb)
  !
END SUBROUTINE
