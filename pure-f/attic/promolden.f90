! Copyright (c) 2018 
! Jose Luis Casals Sainz <jluiscasalssainz@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and 
! Evelio Francisco Miguelez <evelio@uniovi.es>. 
!
! promolden is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
! 
! promolden is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
program promolden

  !$ use omp_lib
  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: rp, ip
  use mod_io, only: ferror, faterr, getdate, stdargs, equal, &
                    isinteger, isreal, isword, ioinit, getline,  &
                    uin, mline, string, flush_unit, &
                    nwarns, ncomms, warning, lgetword
  use mod_param, only: init_param, optsparam
  use mod_wfn, only: init_wfn, end_wfn, optswfn, isdata, loadwfn
  use mod_surf, only: init_surf, optssurf, surface, nangleb
  use mod_fields, only: pointr1, testpointr1
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
  real(kind=rp) :: rval, xpoint(3)

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

  if (index(optv,"d") /= 0) then
    call optsparam('debug',.true.)
    call optsparam('verbose',.true.)
  end if
  if (index(optv,"h") /= 0) then
    write (uout,'(1x,a)') string('# +++ Help not yet available')
    goto 1
  end if

  call init_param ()
  call init_wfn ()
  call init_surf ()
  
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

    else if (equal(word,'unload')) then
      if (isdata) then
        call end_wfn ()
        isdata = .false.
      else        
        call ferror('promolden', 'unload not data available', faterr)  
      end if
     
    ! Field options
    else if (equal(word,'pointr1')) then
      ok = isreal(xpoint(1), line, lp)
      ok = ok .and. isreal(xpoint(2), line, lp)
      ok = ok .and. isreal(xpoint(3), line, lp)
      if (.not.ok) call ferror('promolden', 'wrong point line', faterr) 
      if (isdata) then
        call testpointr1 (xpoint)
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
    
    else if (equal(word,'agrid')) then
      nangleb(1) = 0
      ok = isinteger(ival, line, lp) .and. ival.ne.0
      ival = abs(ival)
      nangleb(4) = ival
      if (.not.ok) call ferror('promolden', 'wrong agrid iqudt', faterr) 
      ok = isinteger(ival, line, lp) .and. ival.ne.0
      ival = abs(ival)
      nangleb(2) = ival
      if (.not.ok) call ferror('promolden', 'wrong agrid ntheta', faterr) 
      ok = isinteger(ival, line, lp) .and. ival.ne.0
      ival = abs(ival)
      nangleb(3) = ival
      if (.not.ok) call ferror('promolden', 'wrong agrid nphi', faterr) 

    else if (equal(word,'lebedev')) then
      ok = isinteger(ival, line, lp)
      ok = ok .and. ival.ne.0_ip
      if (.not.ok) call ferror('promolden', 'wrong agrid lebedev', faterr) 
      ival = abs(ival)
      call good_lebedev (ival)
      nangleb(1) = 1
      nangleb(2) = ival
      nangleb(3) = 0
      nangleb(4) = 0

    else if (equal(word,'surface')) then
      ok = isinteger(ival, line, lp)
      ok = ok .and. ival.ne.0_ip
      if (.not.ok) call ferror('promolden', 'wrong surface line', faterr) 
      ival = abs(ival)
      if (isdata) then
        call surface (ival)
      else        
        call ferror('promolden', 'surface not data available', faterr)  
      end if

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

  if (isdata) call end_wfn ()

1 call cpu_time (tiempo2)
  !$ tiempo2 = omp_get_wtime()
  call getdate (wdate)
  write (uout,'(1x,a,f16.6)') string('# Total elapsed time = '), tiempo2-tiempo1
  write (uout,'(" # Check : (",A," WARNINGS, ",A," COMMENTS)")') string(nwarns), string(ncomms)
  write (uout,'(1x,a)') string('# Calculation ends at '//wdate)

end program
