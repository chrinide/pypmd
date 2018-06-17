! MOs in the WFN file are expressed in terms of unnormalized carte-
! sian gaussians (CG). It could happens that a given CG is repeated.
! In that case, we condensate all the coefficients in a single one.
subroutine setupwfn ()

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: rp, ip
  use mod_memory, only: alloc, free
  use mod_io, only: ferror, faterr
  use mod_wfn, only: nmo, nprims, ncent, coef=>coefnat, maxgrp, mgrp, &
                     numshells, coefcan, &
                     icen, oexp, ityp, npc, nlm, nuexp, oexp, ngroup, &
                     icenat, nzexp, allocate_space_for_shells
  implicit none

  real(kind=rp) :: alph
  integer(kind=ip) :: npa, i, ic, icentro, inda, isud
  integer(kind=ip) :: isuk, itd, itip, itk, j, k, m, npcant
  real(kind=rp), allocatable, dimension(:) :: oexpa
  real(kind=rp), allocatable, dimension(:,:) :: coefa
  integer(kind=ip), allocatable, dimension(:) :: icena
  integer(kind=ip), allocatable, dimension(:) :: itypa

  ! Init
  call alloc ('setupwfn', 'oexpa', oexpa, nprims)
  call alloc ('setupwfn', 'coefa', coefa, nmo, nprims)
  call alloc ('setupwfn', 'icena', icena, nprims)
  call alloc ('setupwfn', 'itypa', itypa, nprims)
  npa = 1_ip
  icena(1) = icen(1)
  oexpa(1) = oexp(1)
  itypa(1) = ityp(1)
  coefa(1:nmo,1) = coef(1:nmo,1)
  cyclej: do j= 2,nprims
    do m = 1,npa
      if (icen(j).eq.icena(m) .and. ityp(j).eq.itypa(m) .and. abs(oexp(j)-oexpa(m)).le.1d-10) then
        coefa(1:nmo,m) = coefa(1:nmo,m)+coef(1:nmo,j)
        cycle cyclej
      end if
    end do
    npa = npa + 1_ip
    icena(npa) = icen(j)
    oexpa(npa) = oexp(j)
    itypa(npa) = ityp(j)
    coefa(1:nmo,npa) = coef(1:nmo,j)
  end do cyclej

  ! Recompute the original variables
  write (uout,2) nmo
  write (uout,1) nprims, npa
  nprims = npa
  do j = 1,nprims
    icen(j) = icena(j)
    oexp(j) = oexpa(j)
    ityp(j) = itypa(j)
    coef(1:nmo,j) = coefa(1:nmo,j)
  end do

  ! Determine primitives corresponding to each center.
  do ic = 1,ncent
    npc(ic) = 0_ip
  end do
  do j = 1,nprims
    ic = icen(j)
    npc(ic) = npc(ic) + 1_ip
    inda = npc(ic)
    icenat(inda,ic) = j
  end do

  ! Classify primitives in each center by types and similar exponents.
  do ic = 1,ncent
    scyclej: do j = 1,npc(ic)
      k = icenat(j,ic)
      itk = ityp(k)
      isuk = nlm(itk,1) + nlm(itk,2) + nlm(itk,3)
      if (j.eq.1) then
        ngroup(ic) = 1_ip
        nzexp(ic,1) = 1_ip
        nuexp(ic,1,1) = k
      else
        do m = 1,ngroup(ic)
          inda = nuexp(ic,m,1)
          itd = ityp(inda)
          isud = nlm(itd,1) + nlm(itd,2) + nlm(itd,3)
          if (abs(oexp(k)-oexp(inda)).lt.1d-8) then
            if (itk.eq.1.and.itd.eq.1) then
              call ferror ('setupwfn', 'two s primitives with equal exponents', faterr)
            else
              if (isuk.eq.isud) then
                nzexp(ic,m) = nzexp(ic,m) + 1_ip
                nuexp(ic,m,nzexp(ic,m)) = k
                cycle scyclej
              end if
            end if
          end if
        end do
        ngroup(ic) = ngroup(ic) + 1_ip
        if (ngroup(ic).gt.mgrp) then
          call ferror ('setupwfn', 'increase mgrp in mod_wfn file', faterr)
        end if
        nzexp(ic,ngroup(ic)) = 1_ip
        nuexp(ic,ngroup(ic),1) = k
      end if   
    end do scyclej
  end do

  ! Reconstruct the values of ityp(),icen(),oexp(), and coef()
  i = 0_ip
  do ic = 1,ncent
    do m = 1,ngroup(ic)
      do k = 1,nzexp(ic,m)
        j = nuexp(ic,m,k)
        alph = oexp(j)
        itip = ityp(j)
        icentro = icen(j)
        i = i + 1_ip
        itypa(i) = itip
        oexpa(i) = alph
        icena(i) = icentro
        coefa(1:nmo,i) = coef(1:nmo,j)
      enddo
    end do
  end do
  
  call free ('setupwfn', 'coefnat', coef)
  call free ('setupwfn', 'coefcan', coefcan)
  call alloc ('setupwfn', 'coefnat', coef, nmo, nprims)
  call alloc ('setupwfn', 'coefcan', coefcan, nmo, nprims)
  do i = 1,nprims
    ityp(i) = itypa(i)
    oexp(i) = oexpa(i)
    icen(i) = icena(i)
    coef(1:nmo,i) = coefa(1:nmo,i)
  end do
 
  npcant = 0_ip
  do ic = 1,ncent
    do k = 1,npc(ic)
      icenat(k,ic) = k + npcant
    end do
    npcant = npcant + npc(ic)
  end do

  ! Reconstruct the values of nuexp(). Now, they are ordered.
  ! Determine also which is the maximum value of ngroup(ic)
  ! Determine the total number of shells.
  i = 0_ip
  numshells = 0_ip
  maxgrp = 0_ip
  do ic = 1,ncent
    numshells = numshells + ngroup(ic)
    if (ngroup(ic).gt.maxgrp) maxgrp = ngroup(ic)
    do m = 1,ngroup(ic)
      do k = 1,nzexp(ic,m)
        i = i + 1_ip
        nuexp(ic,m,k) = i
      end do
    end do
  end do

  ! Allocate space for nshell(), ishell(), and atcenter() arrays.
  call allocate_space_for_shells ()

  ! Deallocate arrays.
  call free ('setupwfn', 'oexpa', oexpa)
  call free ('setupwfn', 'coefa', coefa)
  call free ('setupwfn', 'icena', icena)
  call free ('setupwfn', 'itypa', itypa)

  ! Formats
1 format (1x,'# Input number of Primitives ',i0,' reduced to ',i0)
2 format (1x,'# Number of molecular orbitals ',i0)

end subroutine
