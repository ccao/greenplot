MODULE constants
  !
  implicit none
  !
  integer, parameter :: dp=selected_real_kind(14, 200)
  real(dp), parameter :: twopi=3.141592653589793*2.d0
  real(dp), parameter :: sqrtpi=1.7724538509055158819194275566d0
  real(dp), parameter :: logpi_2=0.572364942924700087071713675675d0
  complex(dp), parameter :: cmplx_i=cmplx(0.d0, 1.d0)
  complex(dp), parameter :: cmplx_0=cmplx(0.d0, 0.d0)
  complex(dp), parameter :: cmplx_1=cmplx(1.d0, 0.d0)
  integer, parameter :: stdin=5
  integer, parameter :: stdout=6
  integer, parameter :: fin=10
  integer, parameter :: fout=11
  real(dp), parameter :: eps3=1.0d-3
  real(dp), parameter :: eps4=1.0d-4
  real(dp), parameter :: eps5=1.0d-5
  real(dp), parameter :: eps6=1.0d-6
  real(dp), parameter :: eps7=1.0d-7
  real(dp), parameter :: eps8=1.0d-8
  real(dp), parameter :: eps9=1.0d-9
  real(dp), parameter :: eps12=1.0d-12
  !
  integer, parameter :: maxiter = 10
  !
END MODULE
