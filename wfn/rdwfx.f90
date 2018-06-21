subroutine rdwfx (wfnfile)

  use mod_prec, only: rp, ip
  use mod_io, only: ferror, faterr, udat, string, warning, &
                    getline_raw
  use mod_wfn, only: maxtype, ncent, nmo, nprims, &
                     occ, oexp, ityp, xyz, rint, rdm, ncore, &
                     atnam, charge, icen, coefcan, c1et, epsocc, &
                     noccupied, occupied, nvirtual, occv, coefnat, &
                     allocate_space_for_wfn, allocate_space_for_rho, &
                     allocate_space_for_rdm, deallocate_space_for_rdm
  implicit none

  character(len=*), intent(in) :: wfnfile

  integer(kind=ip) :: iwfn, i, icounto, icountv
  character(len=:), allocatable :: line, line2 
  logical :: keyw(7)

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
      end if
    end if
  end do
20 continue
  if (any(.not.keyw)) call ferror("rdwfx", "missing array in wfx file", faterr)

  ! Reduce and set info
  call setupwfn ()
  coefcan = coefnat
  
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
