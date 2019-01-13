! Determine the overlap matrix for all atoms or specified atoms
program dorho

  use mod_prec, only: rp, ip, print_kind_info
  !$ use omp_lib, only: omp_get_wtime
  use mod_parallel, only: init_parallel, ichunk, itile, nthreads, &
                          info_parallel, end_parallel
  use mod_io, only: ferror, faterr, getdate, stdargs, equal, &
                    isinteger, isreal, isword, ioinit, getline,  &
                    uin, uout, mline, string, flush_unit, &
                    nwarns, ncomms, warning, lgetword
  use mod_param, only: debug, largwr, init_param, bzip
  use mod_wfn, only: init_wfn, epsdet, epsrho, dmftxc, idmft
  implicit none
 
  character(len=:), allocatable :: optv !< command-line arguments 
  character(len=:), allocatable :: fileroot !< file prefix
  character(len=mline) :: filedat
 
  logical :: ok
  integer(kind=ip) :: lp
  real(kind=rp) :: tiempo1, tiempo2
  character(len=mline) :: wdate, word2, line2  
  character(len=:), allocatable :: line, subline, word

  character(len=6) :: wtmp

  ! Begin program
  call init_parallel ()
  call ioinit ()
  call stdargs (optv, fileroot)
  call getdate (wdate)
  write (uout,'(1x,a)') string('#  =====================================================')
  write (uout,'(1x,a)') string('# |       dorho: OVERLAP MATRICES OVER SURFACES         |')
  write (uout,'(1x,a)') string('# |    (c) E. Francisco, University of Oviedo, 2017     |')
  write (uout,'(1x,a)') string('#  =====================================================')
  write (uout,'(1x,a)') string('#')
  write (uout,'(1x,a)') string('# Calculation starts at '//wdate)
  call info_parallel (uout)
  call cpu_time (tiempo1)
  !$ tiempo1 = omp_get_wtime()
  call init_param ()
  call init_wfn ()
  if (index(optv,"d") /= 0) then
    debug = .true.
    largwr = .true.
    call ferror ('dorho', 'debug mode is enabled', warning)
    write (uout,'(1x,a)') string('# WARNING: DEBUG MODE ENABLED !!')
    write (uout,'(1x,a)') string('# Lot of info will be printed !!')
    write (uout,'(1x,a)') string('# WARNING: DEBUG MODE ENABLED !!')
    write (uout,*)
    call print_kind_info (uout)
  end if
  if (index(optv,"h") /= 0) then
    write (uout,'(1x,a)') string('# +++ Help not yet available')
    goto 1
  end if

  ! First line is data file
  read (uin,'(a)') line2
  lp = 1
  if (isword(word2,line2,lp)) then
    filedat = adjustl(word2)
    filedat = filedat(1:len_trim(filedat))
  else
    call ferror('dorho', 'data file name not read in', faterr)
  end if
  write (uout,'(1x,a)') string('#')
  write (uout,'(1x,a)') string('# +++ Begin to read data file : '//filedat)
  call rdwfn (filedat)                                  
  write (uout,'(1x,a)') string('# +++ End of read data file')
  write (uout,'(1x,a)') string('#')

  ! Init all data
  write (uout,'(1x,a)') string('#')
  write (uout,'(1x,a)') string('# +++ Begin to initialize data')
  write (uout,'(1x,a)') string('# +++ End of initialize data')
  write (uout,'(1x,a)') string('#')

  ! Read the rest of input file
  call flush_unit (uout)
  write (uout,'(1x,a)') string('#')
  write (uout,'(1x,a)') string('# +++ Begin to read input')
  do while (getline(uin,line))
    lp = 1
    word = lgetword(line,lp)
    subline = line(lp:)

    if (equal(word,'#')) then
      continue

    ! Parallel options
    else if (equal(word,'threads')) then
      ok = isinteger(nthreads, line, lp)
      ok = ok .and. nthreads.ne.0
      if (.not.ok) call ferror ('dorho', 'wrong threads line', faterr) 
      !$ nthreads = abs(nthreads)
      !$ write (uout,'(1x,a,1x,i0)') string('# *** Number of threads changed to :'), nthreads
      !$ call omp_set_num_threads(nthreads)

    else if (equal(word,'chunk')) then
      ok = isinteger(ichunk, line, lp)
      ok = ok .and. ichunk.ne.0
      if (.not.ok) call ferror ('dorho', 'wrong chunk line', faterr) 
      !$ ichunk = abs(ichunk)
      !$ write (uout,'(1x,a,1x,i0)') string('# *** Number of chunks changed to :'), ichunk

    else if (equal(word,'tile')) then
      ok = isinteger(itile, line, lp)
      ok = ok .and. itile.ne.0
      if (.not.ok) call ferror ('dorho', 'wrong tile line', faterr) 
      itile = abs(itile)
      write (uout,'(1x,a,1x,i0)') string('# *** Variable tile changed to :'), itile

    ! Param options
    else if (equal(word,'verbose')) then
      largwr = .true.
      write (uout,'(1x,a)') string('# *** Verbose mode is enabled')

    else if (equal(word,'bzip')) then
      bzip = .true.
      write (uout,'(1x,a)') string('# *** Compression mode is enabled')
     
    ! WFN options
    else if (equal(word,'epsdet')) then
      ok = isreal(epsdet, line, lp)
      ok = ok .and. epsdet.ne.0.0_rp
      if (.not.ok) call ferror ('dorho', 'wrong epsdet line', faterr) 
      epsdet = abs(epsdet)
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epsdet changed to :'), epsdet

    else if (equal(word,'epsrho')) then
      ok = isreal(epsrho, line, lp)
      ok = ok .and. epsrho.ne.0.0_rp
      if (.not.ok) call ferror ('dorho', 'wrong epsrho line', faterr) 
      epsrho = abs(epsrho)
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epsrho changed to :'), epsrho

    else if (equal(word,'dmftxc')) then
      ok = isword(wtmp, line, lp)
      if (.not.ok) call ferror ('dorho', 'wrong dmftxc line', faterr) 
      dmftxc = .true.
      if (wtmp(1:2).eq.'BB') then
        idmft = 0_ip
      else if (wtmp(1:4).eq.'BBC1') then
        idmft = 1_ip
      else if (wtmp(1:4).eq.'BBC2') then
        idmft = 2_ip
      else if (wtmp(1:2).eq.'CA') then
        idmft = 3_ip
      else if (wtmp(1:2).eq.'GU') then
        idmft = 4_ip
      else if (wtmp(1:3).eq.'ICA') then
        idmft = 5_ip
      else if (wtmp(1:3).eq.'HYB') then
        idmft = 6_ip
      else if (wtmp(1:5).eq.'PNOF4') then
        idmft = 7_ip
      else if (wtmp(1:5).eq.'PNOF5') then
        idmft = 8_ip
      else
        call ferror ('dorho', 'unknown dmftxc '//string(wtmp), faterr)
      end if

    else if (equal(word,'endorho')) then
      exit

    else if (equal(word,'end')) then
      goto 1

    else
      !nothing to easy common file
                    
    end if
  end do
  write (uout,'(1x,a)') string('# +++ End of read input')
  write (uout,'(1x,a)') string('#')
  call flush_unit (uout)

  ! Some more init and print all final options for calculation
  write (uout,'(1x,a)') string('#')
  write (uout,'(1x,a)') string('# +++ Calculation info')
  call info ()
  write (uout,'(1x,a)') string('# +++ Calculation info')
  write (uout,'(1x,a)') string('#')

  ! Do the work
  write (uout,'(1x,a)') string('#')
  write (uout,'(1x,a)') string('# +++ Begin computing RDM')
  !call binrdm (filedat)
  write (uout,'(1x,a)') string('# +++ End computing RDM')
  write (uout,'(1x,a)') string('#')
 
  ! End program
  call deallocate_space_for_wfn ()
  call deallocate_space_for_shells ()
1 call cpu_time (tiempo2)
  !$ tiempo2 = omp_get_wtime()
  call getdate (wdate)
  write (uout,'(1x,a)') string('#')
  write (uout,'(1x,a,f16.6)') string('# Total elapsed time = '), tiempo2-tiempo1
  write (uout,'(" # Check : (",A," WARNINGS, ",A," COMMENTS)")') string(nwarns), string(ncomms)
  write (uout,'(1x,a)') string('# Calculation ends at '//wdate)
  call end_parallel ()

end program
