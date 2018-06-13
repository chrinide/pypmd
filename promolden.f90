program promolden

  !$ use omp_lib
  use mod_prec, only: rp, ip
  use mod_io, only: ferror, faterr, getdate, stdargs, equal, &
                    isinteger, isreal, isword, ioinit, getline,  &
                    uin, uout, mline, string, flush_unit, &
                    nwarns, ncomms, warning, lgetword
  use mod_param, only: debug, largwr, init_param
  use mod_parallel, only: init_parallel, nthreads, info_parallel, end_parallel
  implicit none
 
  character(len=:), allocatable :: optv !< command-line arguments 
  character(len=:), allocatable :: fileroot !< file prefix
 
  logical :: ok
  integer(kind=ip) :: lp
  real(kind=rp) :: tiempo1, tiempo2
  character(len=mline) :: wdate  
  character(len=:), allocatable :: line, subline, word

  ! Begin program
  call init_parallel ()
  call ioinit ()
  call stdargs (optv, fileroot)

  call getdate (wdate)
  write (uout,'(1x,a)') string('#  =====================================================')
  write (uout,'(1x,a)') string('# |            PROMOLDEN: A QTAIM/IQA code              |')
  write (uout,'(1x,a)') string('# |  (c) E. Francisco, University of Oviedo, 2018       |')
  write (uout,'(1x,a)') string('# |  (c) A. Martin Pendas, University of Oviedo, 2018   |')
  write (uout,'(1x,a)') string('# |  (c) J. L. Casals Sainz                             |')
  write (uout,'(1x,a)') string('#  =====================================================')
  write (uout,'(1x,a)') string('#')
  write (uout,'(1x,a)') string('# Calculation starts at '//wdate)
  call info_parallel (uout)

  call cpu_time (tiempo1)
  !$ tiempo1 = omp_get_wtime()

  call init_param ()
  if (index(optv,"d") /= 0) then
    debug = .true.
    largwr = .true.
    call ferror ('promolden', 'debug mode is enabled', warning)
    write (uout,'(1x,a)') string('# WARNING: DEBUG MODE ENABLED !!')
    write (uout,'(1x,a)') string('# Lot of info will be printed !!')
    write (uout,'(1x,a)') string('# WARNING: DEBUG MODE ENABLED !!')
    write (uout,*)
  end if
  if (index(optv,"h") /= 0) then
    write (uout,'(1x,a)') string('# +++ Help not yet available')
    goto 1
  end if

  ! Read the input file
  call flush_unit (uout)
  write (uout,'(1x,a)') string('# +++ Begin to read input')
  do while (getline(uin,line))
    lp = 1
    word = lgetword(line,lp)
    subline = line(lp:)

    if (equal(word,'#')) then
      continue

    ! Parallel options
    else if (equal(word,'nthreads')) then
      ok = isinteger(nthreads, line, lp)
      ok = ok .and. nthreads.ne.0
      if (.not.ok) call ferror('promolden', 'wrong nthreads line', faterr) 
      !$ nthreads = abs(nthreads)
      !$ write (uout,'(1x,a,1x,i0)') string('# *** Number of threads changed to :'), nthreads
      !$ call omp_set_num_threads(nthreads)

    ! Param options
    else if (equal(word,'verbose')) then
      largwr = .true.
      write (uout,'(1x,a)') string('# *** Verbose mode is enabled')

    else if (equal(word,'end')) then
      exit

    else if (equal(word,'exit') .or. equal(word,'quit')) then
      goto 1

    else
      !nothing yet to easy common file

    end if
  end do
  write (uout,'(1x,a)') string('# +++ End of read input')
  write (uout,'(1x,a)') string('#')
  call flush_unit (uout)

1 call cpu_time (tiempo2)
  !$ tiempo2 = omp_get_wtime()
  call getdate (wdate)
  write (uout,'(1x,a,f16.6)') string('# Total elapsed time = '), tiempo2-tiempo1
  write (uout,'(" # Check : (",A," WARNINGS, ",A," COMMENTS)")') string(nwarns), string(ncomms)
  write (uout,'(1x,a)') string('# Calculation ends at '//wdate)
  call end_parallel ()

end program
