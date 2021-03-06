PROGRAM greenplot
  !
  USE constants,  only : dp, twopi
  USE input
  USE wanndata,   only : read_ham, rvec, nrpt
  USE hamk,       only : nbnd, init_hamk, calc_hamk, finalize_hamk
  USE para
  !
  implicit none
  !
  integer ik, ie, ii
  !
  complex(dp), allocatable :: gf(:, :)
  real(dp), allocatable :: dos(:, :), buf(:)
  !
  CALL init_para()
  !
  CALL read_input()
  !
  CALL split_ekmesh(ne, nk, dmsn)
  !
  CALL read_ham(seed)
  !
  CALL init_hamk()
  !
  if (inode.eq.0) then
    !
    if (mode.eq.0) then
      write(stdout, '(A)') ' Bulk calculation...'
    else
      write(stdout, '(A,3I5,A)') ' [001] Surface calculation '
    endif
    !
    if (dmsn.eq.2) then
      write(stdout, '(A)') '....x....1....x....2....x....3....x....4....x....5'
    endif
    !
  endif
  !
  allocate(dos(ne, nk))
  dos(:, :)=0.d0
  ! 
  allocate(gf(nbnd, nbnd))
  !
  do ik=first_k, last_k
    !
    CALL calc_hamk( kvec(:, ik) )
    !
    do ie=first_ene, last_ene
      !
      if (mode.eq.0) then
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
        dos(ie, ik)=dos(ie, ik)-aimag(gf(ii,ii))/twopi
      enddo
      !
    enddo ! ie
    !
    if (dmsn.lt.2) then
      !
      CALL para_merge(dos(:, ik), ne)
      !
    endif
    !
    if (inode.eq.0) then
      !
      write(stdout, '(A,1F5.1,A)') '# ', (ik-first_k)*100.d0/(last_k-first_k), '% Done'
      flush(stdout)
      !
    endif
    !
  enddo ! ik
  !
  if (dmsn.gt.1) then
    !
    allocate(buf(nk))
    !
    do ii=first_ene, last_ene
      !
      buf(:)=dos(ii, 1:nk)
      CALL para_merge(buf, nk)
      dos(ii,:)=buf(:)
      !
    enddo
    !
    deallocate(buf)
    !
  endif
  !
  if (inode.eq.0) then
    ! output results...
    open(unit=fout, file='plot.dat')
    !
    if (dmsn.eq.0) then
      do ie=1, ne
        write(fout, '(2F16.9)') emesh(ie), dos(ie, 1)
      enddo
    elseif (dmsn.eq.1) then
      ! line mode output
      write(fout, '(2I10)') nk, ne
      write(fout, '(10F16.9)') xk(:)
      write(fout, '(10F16.9)') emesh(:)
      !
      do ik=1, nk
        write(fout, '(10F16.9)') dos(:, ik)
      enddo
    else
      ! 2D mode output
      write(fout, '(A,2I10)') '#  ', nkx, nky
      do ii=1, nky
        do ik=1, nkx
          write(fout, '(2F9.4,4X,1F16.9)') kvec(1:2, (ik-1)*nky+ii), dos(1, (ik-1)*nky+ii)
        enddo
        write(fout, *) ''
      enddo
      !
    endif
    !
    close(unit=fout)
    !
  endif
  !
  deallocate(gf)
  deallocate(dos)
  !
  CALL finalize_hamk()
  CALL finalize_para()
  !
END PROGRAM
