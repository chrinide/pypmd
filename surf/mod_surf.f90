module mod_surf

  implicit none

contains

  subroutine read_surf(verbose, debug)
  
  ! Read the rest of input file
  call flush_unit (uout)
  write (uout,'(1x,a)') string('# +++ Begin to read surf block')
  do while (getline(uin,line))
    lp = 1
    word = lgetword(line,lp)
    subline = line(lp:)

    if (equal(word,'#')) then
      continue

    ! Surf options
    else if (equal(word,'steeper')) then
      ok = isinteger(steeper, line, lp)
      ok = ok .and. steeper.ne.0_ip
      if (.not.ok) call ferror('dosurf', 'wrong steeper line', faterr) 
      steeper = abs(steeper)
      if (steeper.gt.3) then
        call ferror('dosurf', 'wrong steeper value', faterr)  
      end if

    else if (equal(word,'epsiscp')) then
      ok = isreal(epsiscp, line, lp)
      ok = ok .and. epsiscp.ne.0.0_rp
      if (.not.ok) call ferror('dosurf', 'wrong epsiscp line', faterr) 
      epsiscp = abs(epsiscp)
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epsiscp changed to :'), epsiscp

    else if (equal(word,'epsilon')) then
      ok = isreal(epsilon, line, lp)
      ok = ok .and. epsilon.ne.0.0_rp
      if (.not.ok) call ferror('dosurf', 'wrong epsilon line', faterr) 
      epsilon = abs(epsilon)
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epsilon changed to :'), epsilon

    else if (equal(word,'agrid')) then
      ok = isinteger(iqudt, line, lp)
      ok = ok .and. isinteger(ntheta, line, lp)
      ok = ok .and. isinteger(nphi, line, lp)
      ok = ok .and. ntheta.ne.0_ip .and. nphi.ne.0_ip .and. iqudt.ne.0
      if (.not.ok) call ferror('dosurf', 'wrong agrid line', faterr) 
      iqudt = abs(iqudt)
      ntheta = abs(ntheta)
      nphi = abs(nphi)
      nangleb(:,1) = 0
      nangleb(:,2) = ntheta
      nangleb(:,3) = nphi 
      nangleb(:,4) = iqudt
      write (uout,'(1x,a,3(1x,i0))') string('# *** Surface agrid (iqudt,ntheta,nphi) :'), iqudt, ntheta, nphi

    else if (equal(word,'lebedev')) then
      ok = isinteger(npang, line, lp)
      ok = ok .and. npang.ne.0_ip
      if (.not.ok) call ferror('dosurf', 'wrong lebedev line', faterr) 
      npang = abs(npang)
      call good_lebedev (npang)
      nangleb(:,1) = 1
      nangleb(:,2) = npang
      nangleb(:,3) = 0
      nangleb(:,4) = 0
      write (uout,'(1x,a,1x,i0)') string('# *** Surface lebedev changed to :'), npang

    else if (equal(word,'rmaxatom')) then
      ok = isreal(rmax, line, lp)
      ok = ok .and. rmax.ne.0.0_rp
      if (.not.ok) call ferror('dosurf', 'wrong rmaxatom line', faterr) 
      rmaxatom = abs(rmax)
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable rmaxatom changed to :'), rmax

    else if (equal(word,'rmaxsurf')) then
      ok = isreal(rmax, line, lp)
      ok = ok .and. rmax.ne.0.0_rp
      if (.not.ok) call ferror('dosurf', 'wrong rmaxsurf line', faterr) 
      rmaxsurf = abs(rmax)
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable rmaxsurf changed to :'), rmax

    else if (equal(word,'ntrial')) then
      ok = isinteger(ntrial, line, lp)
      ok = ok .and. ntrial.ne.0_ip
      if (.not.ok) call ferror('dosurf', 'wrong ntrial line', faterr) 
      ntrial = abs(ntrial)
      if (mod(ntrial,2).eq.0) ntrial = ntrial + 1_ip
      write (uout,'(1x,a,1x,i0)') string('# *** Variable ntrial changed to :'), ntrial

    else if (equal(word,'rprimer')) then
      ok = isreal(rprimer, line, lp)
      if (.not.ok) call ferror('dosurf', 'wrong rprimer line', faterr) 
      rprimer = abs(rprimer)
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable rprimer changed to :'), rprimer

    else if (equal(word,'rsearch')) then
      ok = isinteger(nsearch,line,lp)
      ok = ok .and. isinteger(nrstart,line,lp)
      ok = ok .and. nsearch.ne.0_ip
      ok = ok .and. nrstart.ne.0_ip
      if (ok) then
        if (nrstart.gt.maxstart) then
          call ferror('dosurf', 'nrstart.gt.maxstart in rsearch order', faterr)
        end if
        nrstart = abs(nrstart)
        if (nsearch.gt.ncent) then
          call ferror('dosurf', 'nsearch.gt.ncent in rsearch order', faterr)
        end if
        if (nsearch.ne.-1_ip) then
          nsearch = abs(nsearch)
          nrsearch(nsearch) = nrstart
          lstart(nsearch) = .true.
          icon = 0_ip
          do while (icon.lt.nrstart)
            if (.not.isreal(rini,line,lp)) then
              call ferror('dosurf', 'bad format in rsearch order', faterr)
            end if
            icon = icon + 1_ip
            rstart(nsearch,icon) = rini
          end do
        else
          forall (ic=1:ncent) nrsearch(ic) = nrstart
          forall (ic=1:ncent) lstart(ic) = .true.
          icon = 0_ip
          do while (icon.lt.nrstart)
            if (.not.isreal(rini,line,lp)) then
              call ferror('dosurf', 'bad format in rsearch order', faterr)
            end if
            icon = icon + 1_ip
            forall (ic=1:ncent) rstart(ic,icon) = rini
          end do
        end if
      end if

    else if (equal(word,'endsurface')) then
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

