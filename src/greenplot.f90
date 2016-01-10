PROGRAM greenplot
  !
  USE constants
  USE input
  USE wanndata
  USE hamk
  USE para
  !
  implicit none
  !
  integer ik, ie, ii
  !
  complex(dp), allocatable :: gf(:, :)
  real(dp), allocatable :: dos(:)
  !
  CALL init_para()
  !
  CALL read_input()
  !
  CALL split_emesh(ne)
  !
  CALL read_ham(seed)
  !
  CALL init_hamk()
  !
  if (inode.eq.0) then
    !
    open(unit=fout, file='plot.dat')
    write(fout, '(2I10)') nk, ne
    write(fout, '(10F16.9)') xk(:)
    write(fout, '(10F16.9)') emesh(:)
    !
  endif
  !
  allocate(dos(ne))
  !
  allocate(gf(nbnd, nbnd))
  !
  do ik=1, nk
    !
    CALL calc_hamk( kvec(ik, :) )
    !
    dos(:)=0.d0
    !
    do ie=first_ene, last_ene
      !
      if (isbulk) then
        !
        CALL calc_bulkgf(gf, emesh(ie)+ef, eta)
        !
      else
        !
        CALL calc_surfgf(gf, emesh(ie)+ef, eta)
        !
      endif
      !
      do ii=1, nbnd
        dos(ie)=dos(ie)+aimag(gf(ii,ii))/twopi
      enddo
      !
      dos(ie) = -dos(ie)
      !
    enddo ! ie
    !
    CALL para_merge(dos, ne)
    !
    write(fout, '(10F16.9)') dos
    !
  enddo ! ik
  !
  deallocate(gf)
  !
  deallocate(dos)
  !
  if (inode.eq.0) close(unit=fout)
  !
  CALL finalize_para()
  !
END PROGRAM
