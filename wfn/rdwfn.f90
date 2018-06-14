subroutine rdwfn (wfnfile)

  use mod_prec, only: rp, ip
  use mod_memory, only: alloc, free
  use mod_io, only: ferror, faterr, mline, udat, string, warning
  use mod_linalg, only: jacobi
  use mod_wfn, only: maxtype, ncent, nmo, nprims, &
                     occ, oexp, ityp, eorb, xyz, rint, rdm, &
                     atnam, charge, icen, coef, c1et, epsocc, &
                     noccupied, occupied, nvirtual, occv, &
                     allocate_space_for_wfn, allocate_space_for_rho, &
                     allocate_space_for_rdm, deallocate_space_for_rdm
  implicit none
 
  ! Arguments
  character(len=*), intent(in) :: wfnfile
 
  ! Local vars
  integer(kind=ip) :: icounto, icountv
  integer(kind=ip) :: i, icol, ifil, iwfn, j, k, nerr, nrot
  real(kind=rp) :: cij, cji, dis, gamma, rho1val, tmp
  real(kind=rp) :: tote, x1, x2, y1, y2, z1, z2
  character(len=mline) :: rholine
  character(len=80) :: wfnttl
  character(len=4) :: mode
  character(len=17) :: label
  character(len=8) :: check

  ! Arrays needed in case of a WFN file based on a CCSD calculation
  real(kind=rp), allocatable, dimension(:,:) :: v1mata
  real(kind=rp), allocatable, dimension(:) :: d1mata

  ! Init data
  open (udat,file=wfnfile,status='old') 
  iwfn = udat
  read (iwfn,101) wfnttl
  read (iwfn,102) mode, nmo, nprims, ncent
  call allocate_space_for_wfn ()
  do i = 1,ncent
    read (iwfn,103) atnam(i),j,(xyz(j,k),k=1,3),charge(j)
  end do

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
    if (occ(i).lt.0.0_rp) then
      call ferror('rdwfn', 'nmo with negative occupation : '//string(i), warning)
    end if  
    read (iwfn,107) (coef(i,j),j=1,nprims)
    coef(i+nmo,1:nprims) = coef(i,1:nprims)
  end do
  read (iwfn,108) check
  if (check .ne. 'END DATA') then
    call ferror('rdwfn', 'end card 1 not found', faterr)
  endif

  ! Reduce and set info
  call setupwfn ()

  ! Special cases
  read (iwfn,109) label,tote,gamma

  ! RDM 
  if (label(1:3).eq.'RDM') then
    rdm = .true.
    call allocate_space_for_rdm ()
    call alloc ('rdwfn', 'v1mata', v1mata , nmo, nmo)
    call alloc ('rdwfn', 'd1mata', d1mata , nmo)
    open (57,file=trim(wfnfile)//".1rdm",status='old',iostat=nerr)
    if (nerr.ne.0) then
      call ferror('rdwfn', 'unable to open 1rdm file', faterr)
    end if
    read (57,'(a)') rholine
    do
      read (57,*,err=123) ifil, icol, rho1val
      if (ifil.gt.nmo .or. icol.gt.nmo) then
        call ferror('rdwfn', 'row and/or column > nmo', faterr)
      else
        c1et(ifil,icol) = rho1val
      end if
    end do
123 close (57)
    ! Test the symmetric character of c1et().
    do i = 2,nmo
      do j = 1,i-1
        cij = c1et(i,j)
        cji = c1et(j,i)
        c1et(i,j) = 0.5_rp*(cij+cji)
        c1et(j,i) = c1et(i,j)
      end do
    end do
    ! Diagonalize total 1st-order matrix
    call jacobi (c1et, d1mata, v1mata, nrot)
    occ(1:nmo) = d1mata(1:nmo)
    do i = 1,nmo
      do j = 1,nprims
        tmp = 0.0_rp
        do k = 1,nmo
          tmp = tmp + v1mata(k,i)*coef(k+nmo,j)
        end do
        coef(i,j) = tmp
      end do
    end do
    call free ('rdwfn', 'v1mata', v1mata)
    call free ('rdwfn', 'd1mata', d1mata)
    call deallocate_space_for_rdm ()
  end if

  noccupied = count(abs(occ)>epsocc)
  nvirtual = nmo - noccupied
  icountv = 0
  icounto = 0
  call allocate_space_for_rho ()
  do i = 1,nmo
    if (abs(occ(i))>epsocc) then
      icounto = icounto + 1_ip
      occupied(icounto) = i
      occv(icounto) = occ(i)
    end if
  end do
 
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
