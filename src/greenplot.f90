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
    !
    write(fout, '(2I10)') nk, ne
    write(fout, '(10F16.9)') xk(:)
    write(fout, '(10F16.9)') emesh(:)
    !
    if (dir.eq.0) then
      write(stdout, '(A)') ' Bulk calculation...'
    else
      write(stdout, '(A,3I5,A)') ' Surface calculation of [', dirvec(:), ']'
    endif
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
      if (dir.eq.0) then
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
    if (inode.eq.0) then
      !
      write(fout, '(10F16.9)') dos
      !
      write(stdout, '(A,1I5,A,1I5,A)') 'Kpt #', ik, ' out of', nk, ' Done'
      !
    endif
    !
  enddo ! ik
  !
  deallocate(gf)
  !
  deallocate(dos)
  !
  if (inode.eq.0) close(unit=fout)
  !
  CALL finalize_hamk()
  CALL finalize_para()
  !
END PROGRAM
