program promolden

  !$ use omp_lib
  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: rp, ip
  use mod_io, only: ferror, faterr, getdate, stdargs, equal, &
                    isinteger, isreal, isword, ioinit, getline,  &
                    uin, mline, string, flush_unit, &
                    nwarns, ncomms, warning, lgetword
  use mod_param, only: isdata
  implicit none
 
  character(len=:), allocatable :: optv !< command-line arguments 
  character(len=:), allocatable :: fileroot !< file prefix
 
  logical :: ok
  integer(kind=ip) :: lp
  real(kind=rp) :: tiempo1, tiempo2
  character(len=mline) :: wdate  
  character(len=:), allocatable :: line, subline, word

  integer(kind=ip) :: ival
  character(len=mline) :: vchar
  real(kind=rp) :: rval, rho, grad(3), gradmod, xpoint(3)

  interface
    subroutine optssurf(var,rval,ival)
      import rp, ip
      implicit none
      character(len=*), intent(in) :: var
      real(kind=rp), optional :: rval
      integer(kind=ip), optional :: ival
    end subroutine
  end interface

  ! Begin program
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

  call cpu_time (tiempo1)
  !$ tiempo1 = omp_get_wtime()

  call init_param ()
  call init_wfn ()
  call init_geom ()
  call init_surf ()
  
  if (index(optv,"d") /= 0) then
    call optsparam('debug',.true.)
    call optsparam('verbose',.true.)
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

    ! Param options
    else if (equal(word,'verbose')) then
      call optsparam(word,.true.)

    ! WFN 
    else if (equal(word,'cuttz')) then
      ok = isreal(rval, line, lp)
      if (.not.ok) call ferror('promolden', 'wrong cuttz line', faterr) 
      call optswfn(word,rval)

    else if (equal(word,'epsocc')) then
      ok = isreal(rval, line, lp)
      if (.not.ok) call ferror('promolden', 'wrong epsocc line', faterr) 
      call optswfn(word,rval)

    else if (equal(word,'epsortho')) then
      ok = isreal(rval, line, lp)
      if (.not.ok) call ferror('promolden', 'wrong epsortho line', faterr) 
      call optswfn(word,rval)

    else if (equal(word,'rmaxatom')) then
      ok = isreal(rval, line, lp)
      if (.not.ok) call ferror('promolden', 'wrong rmaxatom line', faterr) 
      call optswfn(word,rval)

    else if (equal(word,'load')) then
      ok = isword(vchar, line, lp)
      if (.not.ok) call ferror('promolden', 'wrong load line', faterr) 
      call loadwfn(vchar)
      isdata = .true.

    else if (equal(word,'point')) then
      ok = isreal(xpoint(1), line, lp)
      ok = ok .and. isreal(xpoint(2), line, lp)
      ok = ok .and. isreal(xpoint(3), line, lp)
      if (.not.ok) call ferror('promolden', 'wrong point line', faterr) 
      if (isdata) then
        call pointr1(xpoint,rho,grad,gradmod)
        write (*,*) xpoint, rho, grad, gradmod
      else
        call ferror('promolden', 'data not loaded', faterr) 
      end if

    ! Surf options
    else if (equal(word,'steeper')) then
      ok = isinteger(ival, line, lp)
      ok = ok .and. ival.ne.0_ip
      if (.not.ok) call ferror('promolden', 'wrong steeper line', faterr) 
      if (ival.gt.3) then
        call ferror('promolden', 'wrong steeper value', faterr)  
      end if
      ival = abs(ival)
      call optssurf(word,ival=ival)

    else if (equal(word,'epsiscp')) then
      ok = isreal(rval, line, lp)
      ok = ok .and. rval.ne.0.0_rp
      if (.not.ok) call ferror('promolden', 'wrong epsiscp line', faterr) 
      rval = abs(rval)
      call optssurf(word,rval=rval)

    else if (equal(word,'epsilon')) then
      ok = isreal(rval, line, lp)
      ok = ok .and. rval.ne.0.0_rp
      if (.not.ok) call ferror('promolden', 'wrong epsilon line', faterr) 
      rval = abs(rval)
      call optssurf(word,rval=rval)

    else if (equal(word,'rmaxsurf')) then
      ok = isreal(rval, line, lp)
      ok = ok .and. rval.ne.0.0_rp
      if (.not.ok) call ferror('promolden', 'wrong rmaxsurf line', faterr) 
      rval = abs(rval)
      call optssurf(word,rval=rval)

    else if (equal(word,'ntrial')) then
      ok = isinteger(ival, line, lp)
      ok = ok .and. ival.ne.0_ip
      if (.not.ok) call ferror('promolden', 'wrong ntrial line', faterr) 
      ival = abs(ival)
      if (mod(ival,2).eq.0) ival = ival + 1_ip
      call optssurf(word,ival=ival)

    else if (equal(word,'rprimer')) then
      ok = isreal(rval, line, lp)
      ok = ok .and. rval.ne.0.0_rp
      if (.not.ok) call ferror('promolden', 'wrong rprimer line', faterr) 
      rval = abs(rval)
      call optssurf(word,rval=rval)
    
    ! Geometry analysis
    else if (equal(word,'geometry')) then
      call loadgeom()

    ! End of input
    else if (equal(word,'end')) then
      exit

    ! Exit main driver
    else if (equal(word,'exit') .or. equal(word,'quit')) then
      goto 1

    else
      !nothing yet to easy common file

    end if
  end do
  write (uout,'(1x,a)') string('# +++ End of read input')
  write (uout,'(1x,a)') string('#')
  call flush_unit (uout)

  call end_wfn ()
  call end_surf ()

1 call cpu_time (tiempo2)
  !$ tiempo2 = omp_get_wtime()
  call getdate (wdate)
  write (uout,'(1x,a,f16.6)') string('# Total elapsed time = '), tiempo2-tiempo1
  write (uout,'(" # Check : (",A," WARNINGS, ",A," COMMENTS)")') string(nwarns), string(ncomms)
  write (uout,'(1x,a)') string('# Calculation ends at '//wdate)

end program
