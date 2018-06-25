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
subroutine rdwfx (wfnfile)

  use mod_prec, only: rp, ip
  use mod_memory, only: free, alloc
  use mod_io, only: ferror, faterr, udat, string, warning, &
                    getline_raw, mline
  use mod_wfn, only: maxtype, ncent, nmo, nprims, &
                     occ, oexp, ityp, xyz, rint, rdm, ncore, &
                     atnam, charge, icen, coefcan, c1et, epsocc, &
                     noccupied, occupied, nvirtual, occv, coefnat, &
                     allocate_space_for_wfn, allocate_space_for_rho, &
                     allocate_space_for_rdm, deallocate_space_for_rdm
  implicit none

  character(len=*), intent(in) :: wfnfile

  logical :: keyw(7)
  integer(kind=ip) :: iwfn, i, icounto, icountv
  character(len=:), allocatable :: line, line2 
  real(kind=rp), allocatable, dimension(:,:) :: v1mata
  real(kind=rp), allocatable, dimension(:) :: d1mata
  real(kind=rp), allocatable, dimension(:,:) :: xyztmp
  integer(kind=ip) :: icol, ifil, j, k, nerr, nrot
  real(kind=rp) :: cij, cji, dis, rho1val, tmp
  real(kind=rp) :: x1, x2, y1, y2, z1, z2
  character(len=mline) :: rholine

  ! Init data
  open (udat,file=wfnfile,status='old') 
  iwfn = udat

  ! Get initial data
  ncent = 0_ip
  nmo = 0_ip
  nprims = 0_ip
  ncore = 0_ip
  do while (getline_raw(iwfn,line))
    if (len_trim(line) < 1_ip) exit
    if (line(1:1) == "<" .and. line(2:2) /= "/") then
      if (trim(line) == "<Number of Occupied Molecular Orbitals>") then
        read (iwfn,*) nmo
      else if (trim(line) == "</Number of Core Electrons>") then
        read (iwfn,*) ncore
        if (ncore > 0_ip) then
          call ferror("rdwfx", "ECP not yet available", faterr)
        end if
      else if (trim(line) == "<Number of Primitives>") then
        read (iwfn,*) nprims
      else if (trim(line) == "<Number of Nuclei>") then
        read (iwfn,*) ncent
      end if
    end if
  end do

  if (ncore > 0_ip) then
    call ferror("rdwfx", "ECP not yet supported", faterr)
  end if
  if (ncent == 0_ip) then
    call ferror("rdwfx", "number of nuclei not found", faterr)
  end if
  if (nmo == 0_ip) then
    call ferror("rdwfx", "number of molecular orbitals not found", faterr)
  end if
  if (nprims == 0_ip) then
    call ferror("rdwfx", "number of primitives not found", faterr)
  end if
  call allocate_space_for_wfn ()

  ! second pass
  rewind (iwfn)
  keyw = .false.
  do while (.true.)
    read (iwfn,'(A)',end=20) line
    line2 = adjustl(line)
    line = line2
    if (line(1:1) == "<" .and. line(2:2) /= "/") then
      if (trim(line) == "<Primitive Centers>") then
        icen = wfx_read_integers(iwfn,nprims)
        keyw(1) = .true.
      else if (trim(line) == "<Primitive Types>") then
        ityp = wfx_read_integers(iwfn,nprims)
        if (any(ityp(1:nprims) > maxtype)) then
          call ferror('rdwfx', 'cannot work with , i- or higher primitives', faterr)
        end if
        keyw(2) = .true.
      else if (trim(line) == "<Primitive Exponents>") then
        oexp = wfx_read_reals1(iwfn,nprims)
        keyw(3) = .true.
      else if (trim(line) == "<1-RDM>") then
        rdm = .true.
      else if (trim(line) == "<Molecular Orbital Occupation Numbers>") then
        occ = wfx_read_reals1(iwfn,nmo)
        do i = 1,nmo
          if (occ(i).lt.0.0_rp) then
            call ferror('rdwfx', 'negative occupations', warning)
          end if  
        end do
        keyw(4) = .true.
      else if (trim(line) == "<Molecular Orbital Primitive Coefficients>") then
        read (iwfn,*)
        do i = 1,nmo
          read (iwfn,*)
          read (iwfn,*)
          coefnat(i,:) = wfx_read_reals1(iwfn,nprims)
        end do
        keyw(5) = .true.
      else if (trim(line) == "<Nuclear Charges>") then ! Check this for ECP
        charge = wfx_read_reals1(iwfn,ncent)
        do i = 1,ncent
          !atnam(i) = nameguess(int(charge(i),4))
        end do
        keyw(6) = .true.
      else if (trim(line) == "<Nuclear Cartesian Coordinates>") then
        call alloc ('rdwfx', 'xyztmp', xyztmp, 3, ncent)
        xyztmp = reshape(wfx_read_reals1(iwfn,3*ncent),shape(xyztmp))
        xyztmp = reshape(xyztmp,(/3,ncent/))
        do i = 1,ncent 
          do j = 1,3
            xyz(i,j) = xyztmp(j,i)
          end do
        end do
        call free ('rdwfx', 'xyztmp', xyztmp)
        keyw(7) = .true.
      end if
    end if
  end do
20 continue
  if (any(.not.keyw)) call ferror("rdwfx", "missing array in wfx file", faterr)

  ! Evaluate internuclear distances
  ! TODO:mv out
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

  ! Reduce and set info
  call setupwfn ()
  coefcan = coefnat

  ! RDM 
  if (rdm) then
    call allocate_space_for_rdm ()
    call alloc ('rdwfx', 'v1mata', v1mata , nmo, nmo)
    call alloc ('rdwfx', 'd1mata', d1mata , nmo)
    open (57,file=trim(wfnfile)//".1rdm",status='old',iostat=nerr)
    if (nerr.ne.0) then
      call ferror('rdwfx', 'unable to open 1rdm file', faterr)
    end if
    read (57,'(a)') rholine
    do
      read (57,*,err=123) ifil, icol, rho1val
      if (ifil.gt.nmo .or. icol.gt.nmo) then
        call ferror('rdwfx', 'row and/or column > nmo', faterr)
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
    !TODO:change to matmult
    do i = 1,nmo
      do j = 1,nprims
        tmp = 0.0_rp
        do k = 1,nmo
          tmp = tmp + v1mata(k,i)*coefcan(k+nmo,j)
        end do
        coefnat(i,j) = tmp
      end do
    end do
    call free ('rdwfx', 'v1mata', v1mata)
    call free ('rdwfx', 'd1mata', d1mata)
    call deallocate_space_for_rdm ()
  end if
  
  ! Get only occupied orbital for density
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

contains

  !> Read a list of n integers from a logical unit
  function wfx_read_integers(lu,n) result(x)

    use mod_io, only: getline_raw, isinteger, ferror, faterr
    implicit none

    integer, intent(in) :: lu, n
    integer :: x(n)

    integer :: kk, lp, idum
    character(len=:), allocatable :: line
    logical :: ok

    kk = 0
    lp = 1
    ok = getline_raw(lu,line,.true.)
    do while (.true.)
      if (.not.isinteger(idum,line,lp)) then
        lp = 1
        ok = getline_raw(lu,line)
        if (.not.ok .or. line(1:2) == "</") exit
      else
        kk = kk + 1
        if (kk > n) call ferror("rdwfx", "read integers exceeded size of the array", faterr)
        x(kk) = idum
      end if
    end do

  end function wfx_read_integers

  !> Read a list of n reals from a logical unit
  function wfx_read_reals1(lu,n) result(x)

    use mod_io, only: getline_raw, isreal, ferror, faterr
    implicit none

    integer, intent(in) :: lu, n
    real(kind=rp) :: x(n)

    integer :: kk, lp
    real(kind=rp) :: rdum
    character(len=:), allocatable :: line
    logical :: ok

    kk = 0
    lp = 1
    ok = getline_raw(lu,line,.true.)
    do while (.true.)
      if (.not.isreal(rdum,line,lp)) then
        lp = 1
        ok = getline_raw(lu,line)
        if (.not.ok .or. line(1:1) == "<") exit
      else
        kk = kk + 1
        if (kk > n) call ferror("rdwfx", "read reals exceeded size of the array", faterr)
        x(kk) = rdum
      end if
    end do

  end function wfx_read_reals1

end subroutine
