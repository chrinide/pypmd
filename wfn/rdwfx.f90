subroutine rdwfx (wfnfile)

  use mod_prec, only: rp, ip
  use mod_wfn, only: ncent, nmo, nprims, ncore
  use mod_io, only: ferror, faterr, udat, getline_raw
  implicit none

  character(len=*), intent(in) :: wfnfile

  integer(kind=ip) :: iwfn
  character(len=:), allocatable :: line, line2 
  logical :: keyw(7)

  ! Init data
  open (udat,file=wfnfile,status='old') 
  iwfn = udat

  ! Get initial data
  ncent = 0
  nmo = 0
  nprims = 0
  ncore = 0
  do while (getline_raw(iwfn,line))
    if (len_trim(line) < 1) exit
    if (line(1:1) == "<" .and. line(2:2) /= "/") then
      if (trim(line) == "<Number of Occupied Molecular Orbitals>") then
        read (iwfn,*) nmo
      else if (trim(line) == "</Number of Core Electrons>") then
        read (iwfn,*) ncore
        if (ncore > 0) then
          call ferror("rdwfx", "ECP not yet available", faterr)
        end if
      else if (trim(line) == "<Number of Primitives>") then
        read (iwfn,*) nprims
      else if (trim(line) == "<Number of Nuclei>") then
        read (iwfn,*) ncent
      end if
    end if
  end do
  if (ncent == 0) call ferror("rdwfx", "number of nuclei not found", faterr)
  if (nmo == 0) call ferror("rdwfx", "number of molecular orbitals not found", faterr)
  if (nprims == 0) call ferror("rdwfx", "number of primitives not found", faterr)

  ! second pass
  rewind (iwfn)
  keyw = .false.
 
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
