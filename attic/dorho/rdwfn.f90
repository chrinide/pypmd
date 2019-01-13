! Read wavefunction with slightly modified AIMPAC routine.
subroutine rdwfn (wfnfile)

  use mod_prec, only: rp, ip
  use mod_memory, only: alloc, free
  use mod_io, only: ferror, faterr, mline, udat, string, warning
  use mod_wfn, only: cciqa, maxtype, ncent, nmo, nprims, mgrp, rohf, &
                     ngtoH, occ, oexp, ityp, eorb, xyz, rint, uhf, &
                     atnam, charge, icen, coef, mpiqa, rhf, totel, &
                     wfnat, nalpha, nbeta, ndou, nsin, ialpha, ibeta
  implicit none
 
  ! Arguments
  character(len=*), intent(in) :: wfnfile
 
  ! Local vars
  integer(kind=ip) :: i, iwfn, j, k 
  real(kind=rp) :: dis, gamma, occmin, occmax
  real(kind=rp) :: tote, x1, x2, y1, y2, z1, z2
  character(len=80) :: wfnttl
  character(len=4) :: mode
  character(len=17) :: label
  character(len=8) :: check

  ! Init data
  open (udat,file=wfnfile,status='old') 
  iwfn = udat
  read (iwfn,101) wfnttl
  read (iwfn,102) mode, nmo, nprims, ncent
  call allocate_space_for_wfn (ncent,nprims,nmo,mgrp,ngtoH)
  do i = 1,ncent
    read (iwfn,103) atnam(i),j,(xyz(j,k),k=1,3),charge(j)
  end do
  occmin =  2.0_rp
  occmax = -2.0_rp
  totel = 0.0_rp

  ! Evaluate internuclear distances
  do i = 1,ncent
    x1 = xyz(i,1)
    y1 = xyz(i,2)
    z1 = xyz(i,3)
    do j = 1,i
      x2 = xyz(j,1)
      y2 = xyz(j,2)
      z2 = xyz(j,3)
      dis = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
      rint(i,j) = dis
      rint(j,i) = dis
    end do
  end do

  ! Read 
  read (iwfn,104) (icen(i),i=1,nprims)
  read (iwfn,104) (ityp(i),i=1,nprims)
  read (iwfn,105) (oexp(i),i=1,nprims)
  do i = 1,nprims
    if (ityp(i).gt.maxtype) then
      call ferror('rdwfn', 'cannot work with , i- or higher primitives', faterr)
    end if
  end do
  do i = 1,nmo
    read (iwfn,106) occ(i), eorb(i) 
    if (occ(i).gt.occmax) occmax = occ(i)
    if (occ(i).lt.occmin) occmin = occ(i)
    if (occ(i).lt.0.0_rp) then
      if (abs(abs(occ(i))-1.0_rp).lt.1d-6) then
        nbeta = nbeta + 1_ip
        ibeta(nbeta) = i
        nsin = nsin + 1_ip
      else
        wfnat = .true.
      end if
    else
      if (abs(occ(i)-1.0_rp).lt.1d-6) then
        nalpha = nalpha + 1_ip
        ialpha(nalpha) = i
        nsin = nsin + 1_ip
      else if (abs(occ(i)-2.0_rp).lt.1d-6) then
        nalpha = nalpha + 1_ip
        nbeta = nbeta + 1_ip
        ialpha(nalpha) = i
        ibeta(nbeta) = i
        ndou = ndou + 1_ip
      else
        wfnat = .true.
      end if
    end if
    if (occ(i).lt.0.0_rp) then
      call ferror('rdwfn', 'nmo with negative occupation : '//string(i), warning)
    end if  
    read (iwfn,107) (coef(i,j),j=1,nprims)
    coef(i+nmo,1:nprims) = coef(i,1:nprims)
  end do
  totel = sum(occ)
  if (abs(occmin-2.0_rp).lt.1d-6) then
    rhf = .true.
  else if (abs(abs(occmin)-1.0_rp).lt.1d-6 .and. abs(occmax-1.0_rp).lt.1d-6) then
    uhf = .true.
  elseif (abs(occmin-1.0_rp).lt.1d-6 .and. abs(occmax-2.0_rp).lt.1d-6) then
    rohf = .true.
  end if
  read (iwfn,108) check
  if (check .ne. 'END DATA') then
    call ferror('rdwfn', 'end card 1 not found', faterr)
  endif

  ! Reduce and set info
  call setupwfn ()

  ! Special cases
  read (iwfn,109) label,tote,gamma

  ! CCSD calculation 
  if (label(1:5).eq.'CCIQA') then
    cciqa = .true.
  else if (label(1:5).eq.'MPIQA') then
    mpiqa = .true.
  end if
 
  close (iwfn)

  ! formats
101 format (a80)
102 format (4x,a4,10x,3(i5,15x))
103 format (a8,11x,i3,2x,3f12.8,10x,f5.1)
104 format (20x,20i3)
105 format (10x,5e14.7)
106 format (35x,f12.8,15x,f12.8)
107 format (5e16.8)
108 format (a8)
109 format (a17,f20.12,18x,f13.8)
 
  end subroutine
