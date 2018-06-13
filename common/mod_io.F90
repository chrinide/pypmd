module mod_io

  implicit none
  public

  private :: string_int4
  private :: string_int8
  private :: string_real8
  private :: string_real4
  private :: string_char
  private :: string_logical

  public :: fourchar
  public :: string
  public :: equal
  public :: getline
  public :: getline_raw
  public :: upper
  public :: lower
  public :: getword
  public :: lgetword
  public :: isletter
  public :: isdigit
  public :: isinteger
  public :: isreal
  public :: isword
  public :: isexpression
  public :: isexpression_or_word
  public :: ioinit
  public :: fopen_read
  public :: fopen_write
  public :: fopen_append
  public :: fopen_scratch
  public :: fclose
  public :: falloc
  public :: fdealloc
  public :: ferror
  public :: tictac
  public :: flush_unit
  public :: getdate

  ! Justification parameters
  integer, parameter :: ioj_left = -1
  integer, parameter :: ioj_right = -2
  integer, parameter :: ioj_center = -3
  integer, parameter :: ioj_atpoint = -4

  ! General string-building routine
  interface string
     module procedure string_int4
     module procedure string_int8
     module procedure string_real8
     module procedure string_real4
     module procedure string_char
     module procedure string_logical
  end interface string

  ! character parameters
  character(len=1), parameter :: tab = char(9) !< tab char (ascii 9)

  ! string parameters
  integer, parameter, private :: maxlen = 255 !< maximum length of a number
  integer, parameter :: mline = 255 !< maximum length of string

  ! main input and output, logical unit allocation
  integer :: uin !< main input lu
  integer :: uout !< main output lu
  integer :: uerr !< main error lu
  integer :: udat !< main data lu
  integer :: ucopy !< logical unit where the copy of the input is written
  logical, private :: alloc(0:200) !< allocation flag array
  character(len=:), allocatable :: filepath !< relative path to find related files
  logical :: interactive !< is this an interactive sesion?

  ! error system
  integer, parameter :: faterr = -1 !< fatal error flag
  integer, parameter :: warning = 1 !< warning flag
  integer, parameter :: noerr = 0 !< info flag
  integer :: nwarns = 0 !< Number of warnings
  integer :: ncomms = 0 !< Number of comments
  
  ! directory separator
  character(len=1), parameter :: dirsep = "/"

  ! clock
  integer, private :: stime0
  real, private :: ctime0

contains

  !> Connect input files to units with standard defaults.
  subroutine stdargs(optv,uroot)

    use iso_fortran_env, only: input_unit, output_unit
    
    character(len=:), allocatable, intent(out) :: optv !< Dash-options passed to the program
    character(len=:), allocatable, intent(out) :: uroot !< file root for the run

    integer, parameter :: arglen = 1024

    integer :: n 
    integer :: argc, idx
    character(len=arglen) :: argv, aux

    optv = ""
    uin = input_unit
    uout = output_unit
    interactive = .true.
    ucopy = uout
    uroot = "stdin"
    filepath = "."

    argc = command_argument_count()
    if (argc > 0) then
       n = 0
       do while (n < argc) 
          n = n + 1
          call getarg(n,argv)    
          argv = adjustl(argv)
          if (argv(1:1) == "-") then
             optv = trim(optv)//argv(2:2)
          else if (uin == input_unit) then
             uroot = trim(argv)
             if (index(uroot,dirsep) > 0) then
                filepath = uroot(1:index(uroot,dirsep,.true.)-1)
                if (len_trim(filepath) == 0) filepath = "."
                ! xx(note1)xx
                ! gfortran-5 and gfortran-6 have problems with the much simpler:
                !   uroot = uroot(index(uroot,dirsep,.true.)+1:)
                ! where they bail out with internal compiler error. Strangely,
                ! gfortran-4 works perfectly. -- A.
                idx = index(uroot,dirsep,.true.)
                aux = uroot(idx+1:)
                uroot = trim(aux)
             end if
             if (index(uroot,'.') > 0) then
                ! see xx(note1)xx, mod_io.f90
                idx = index(uroot,'.',.true.)
                aux = uroot(1:idx-1)
                uroot = trim(aux)
             end if
             uroot = trim(adjustl(uroot))
             uin = fopen_read(argv,abspath0=.true.)
             interactive = .false.
          else if (uout == output_unit) then
             uout = fopen_write(argv)
          end if
       end do
    end if

  end subroutine stdargs

  !> Build a string from an integer*4. length is the length of the
  !> output string. If the length given by the user is not enough to
  !> contain the minimal representation for the integer, then the
  !> smallest possible length is used instead. Justify selects the
  !> justification of the integer in the field (left, center, and
  !> right). pad0, pad the full length with zeros. padspace, pad with
  !> this number of spaces. If padspace > 0, pad in the left, if < 0,
  !> in the right.
  function string_int4(a,length,justify,pad0,padspace) result(s)

    character(len=:), allocatable :: s
    integer(kind=4), intent(in) :: a
    integer, intent(in), optional :: length
    integer, intent(in), optional :: justify
    logical, intent(in), optional :: pad0
    integer, intent(in), optional :: padspace

    logical :: ispad
    character(len=maxlen) :: aux, fmt
    integer :: i, ialen, ipad, ialen0
    
    if (present(pad0)) then
       ispad = pad0
    else 
       ispad = .false.
    end if
    ialen = 0
    if (present(length)) then
       if (length > 0) then
          ialen = length
       end if
    end if
    write (aux,'(I0)') a
    ialen0 = len_trim(aux)
    if (ialen == 0) then
       ialen = ialen0
    else
       ialen = max(ialen,ialen0)
       if (ispad) then
          write (fmt,'("(I0.",I0,")")') ialen
          write (aux,fmt) a
       else
          write (aux,'(I0)') a
       end if
    end if
    aux = adjustl(aux)
    s = aux(1:ialen)

    if (present(justify)) then
       if (justify == ioj_right) then
          s = adjustr(s)
       else if (justify == ioj_center) then
          ialen = len_trim(s)
          if (len(s) /= ialen) then
             i = (len(s) - ialen)/2 + 1
             s(i:i+ialen) = trim(s)
             s(1:i-1) = ""
             s(i+ialen+1:) = ""
          end if
       end if
    end if

    if (present(padspace)) then
       if (padspace >= 0) then
          do ipad = padspace,0,-1
             if (s(len(s)-ipad+1:) == "") exit
          end do
          ! see xx(note1)xx, mod_io.f90
          aux = s(1:len(s)-ipad)
          s(ipad+1:) = trim(adjustl(aux))
          s(1:ipad) = ""
       else
          do ipad = -padspace,0,-1
             if (s(1:ipad) == "") exit
          end do
          ialen = len(s)-ipad
          ! see xx(note1)xx, mod_io.f90
          aux = s(ipad+1:)
          s(1:ialen) = trim(adjustl(aux))
          s(ialen+1:) = ""
       end if
    end if

  end function string_int4

  !> Build a string from an integer*8. length is the length of the
  !> output string. If the length given by the user is not enough to
  !> contain the minimal representation for the integer, then the
  !> smallest possible length is used instead. Justify selects the
  !> justification of the integer in the field (left, center, and
  !> right). pad0, pad the full length with zeros. padspace, pad with
  !> this number of spaces. If padspace > 0, pad in the left, if < 0,
  !> in the right.
  function string_int8(a,length,justify,pad0,padspace) result(s)

    character(len=:), allocatable :: s
    integer(kind=8), intent(in) :: a
    integer, intent(in), optional :: length
    integer, intent(in), optional :: justify
    logical, intent(in), optional :: pad0
    integer, intent(in), optional :: padspace

    logical :: ispad
    character(len=maxlen) :: aux, fmt
    integer :: i, ialen, ipad, ialen0
    
    if (present(pad0)) then
       ispad = pad0
    else 
       ispad = .false.
    end if
    ialen = 0
    if (present(length)) then
       if (length > 0) then
          ialen = length
       end if
    end if
    write (aux,'(I0)') a
    ialen0 = len_trim(aux)
    if (ialen == 0) then
       ialen = ialen0
    else
       ialen = max(ialen,ialen0)
       if (ispad) then
          write (fmt,'("(I0.",I0,")")') ialen
          write (aux,fmt) a
       else
          write (aux,'(I0)') a
       end if
    end if
    aux = adjustl(aux)
    s = aux(1:ialen)

    if (present(justify)) then
       if (justify == ioj_right) then
          s = adjustr(s)
       elseif (justify == ioj_center) then
          ialen = len_trim(s)
          i = (len(s) - ialen)/2 + 1
          s(i:i+ialen) = trim(s)
          s(1:i-1) = ""
          s(i+ialen+1:) = ""
       end if
    end if

    if (present(padspace)) then
       if (padspace >= 0) then
          do ipad = padspace,0,-1
             if (s(len(s)-ipad+1:) == "") exit
          end do
          ! see xx(note1)xx, mod_io.f90
          aux = s(1:len(s)-ipad)
          s(ipad+1:) = trim(adjustl(aux))
          s(1:ipad) = ""
       else
          do ipad = -padspace,0,-1
             if (s(1:ipad) == "") exit
          end do
          ialen = len(s)-ipad
          ! see xx(note1)xx, mod_io.f90
          aux = s(ipad+1:)
          s(1:ialen) = trim(adjustl(aux))
          s(ialen+1:) = ""
       end if
    end if

  end function string_int8

  !> Build a string from a real*8 number. desc is the descriptor used
  !> (f, g, e, d, etc.). length is the length of the output string. If
  !> the length given by the user is not enough to contain the minimal
  !> representation for the number, then the smallest possible length is
  !> used instead. decimal is the number of decimal places (default:
  !> 7). Justify selects the justification of the integer in the field
  !> (left, center, and right). If a postive value is given to justify,
  !> then align the point to that character. padspace, pad with this
  !> number of spaces. If padspace > 0, pad in the left, if < 0, in the
  !> right.
  function string_real8(a,desc,length,decimal,justify,padspace) result(s)

    character(len=:), allocatable :: s
    real(kind=8), intent(in) :: a
    character(len=1), intent(in) :: desc
    integer, intent(in), optional :: length
    integer, intent(in), optional :: decimal
    integer, intent(in), optional :: justify
    integer, intent(in), optional :: padspace

    integer, parameter :: decimal_default = 7

    character(len=maxlen) :: aux, fmt
    integer :: i, ialen, ipad, ipad0, idx, dec
    logical :: islen, first
    
    islen = .false.
    if (present(length)) islen = (length > 0)
    if (present(decimal)) then
       dec = decimal
    else
       dec = decimal_default
    end if

    ! We'll need two passes through here if length is present but
    ! decimal is not. The first is to measure the length with the
    ! default number of decimal places. The second sets the number of
    ! decimals to match the requested length.
    first = .true.
10  continue
    if (desc == 'e' .or. desc == 'E' .or. desc == 'd' .or. desc == 'D') then
       write (fmt,'("(1p,",A1,I0,".",I0,")")') desc, maxlen, dec
    else
       write (fmt,'("(",A1,I0,".",I0,")")') desc, maxlen, dec
    end if
    write (aux,fmt) a
    aux = adjustl(aux)
    ialen = len_trim(aux)
    if (islen .and..not.present(decimal)) then
       if (ialen > length .and. first) then
          ! Calculate the number of decimals to get that length, then
          ! loop back and write the actual string.
          dec = max(decimal_default - (ialen-length),1)
          first = .false.
          goto 10
       end if
    end if

    if (islen) ialen = max(length,ialen)
    s = aux(1:ialen)

    ipad0 = 0
    if (present(justify)) then
       if (justify == ioj_right) then
          s = adjustr(s)
       else if (justify == ioj_center) then
          ialen = len_trim(s)
          i = (len(s) - ialen)/2 + 1
          s(i:i+ialen) = trim(s)
          s(1:i-1) = ""
          s(i+ialen+1:) = ""
       else if (justify > 0) then
          idx = index(s,".")
          if (idx > 0) then
             ipad0 = max(justify - idx,0)
          end if
       end if
    end if

    if (present(padspace)) then
       if (ipad0 == 0) then
          ipad0 = padspace
       else
          ipad0 = ipad0 + padspace
       end if
    end if
    if (ipad0 > 0) then
       do ipad = ipad0,0,-1
          if (s(len(s)-ipad+1:) == "") exit
       end do
       s(ipad+1:) = s(1:len(s)-ipad)
       s(1:ipad) = ""
    else if (ipad0 < 0) then
       do ipad = -ipad0,0,-1
          if (s(1:ipad) == "") exit
       end do
       ialen = len(s)-ipad
       s(1:ialen) = s(ipad+1:)
       s(ialen+1:) = ""
    end if

  end function string_real8

  !> Use the real*8 routine for real*4.
  function string_real4(a,desc,length,decimal,justify,padspace) result(s)

    character(len=:), allocatable :: s
    real(kind=4), intent(in) :: a
    character(len=1), intent(in) :: desc
    integer, intent(in), optional :: length
    integer, intent(in), optional :: decimal
    integer, intent(in), optional :: justify
    integer, intent(in), optional :: padspace
    
    s = string_real8(real(a,8),desc,length,decimal,justify,padspace)

  end function string_real4

  !> Build a string from another string. length is the length of the
  !> output string. If the length given by the user is not enough to
  !> contain the original string, then it is ignored. Justify selects
  !> the justification of the integer in the field (left, center, and
  !> right). padspace, pad with this number of spaces. If padspace > 0,
  !> pad in the left, if < 0, in the right.
  function string_char(a,length,justify,padspace) result(s)

    character(len=:), allocatable :: s
    character(len=*), intent(in) :: a
    integer, intent(in), optional :: length
    integer, intent(in), optional :: justify
    integer, intent(in), optional :: padspace

    integer :: i, ialen, ipad, i2len
    
    i2len = len_trim(adjustl(a))
    ialen = i2len
    if (present(length)) then
       if (length > ialen) then
          ialen = length
       end if
    end if
    allocate(character(ialen)::s)
    s(1:i2len) = trim(adjustl(a))
    s(i2len+1:ialen) = ""

    if (present(justify)) then
       if (justify == ioj_right) then
          s = adjustr(s)
       else if (justify == ioj_center) then
          ialen = len_trim(s)
          i = (len(s) - ialen)/2 + 1
          s(i:min(i+ialen,len(s))) = trim(s)
          s(1:i-1) = ""
          s(i+ialen+1:) = ""
       end if
    end if
    
    if (present(padspace)) then
       if (padspace >= 0) then
          do ipad = padspace,0,-1
             if (s(len(s)-ipad+1:) == "") exit
          end do
          s(ipad+1:) = s(1:len(s)-ipad)
          s(1:ipad) = ""
       else
          do ipad = -padspace,0,-1
             if (s(1:ipad) == "") exit
          end do
          ialen = len(s)-ipad
          s(1:ialen) = s(ipad+1:)
          s(ialen+1:) = ""
       end if
    end if

  end function string_char

  !> Build a string from a logical value. Returns "Yes"
  !> if input is true or "No" if it is false.
  function string_logical(a) result(s)

    character(len=:), allocatable :: s
    logical, intent(in) :: a

    if (a) then
       s = "Yes"
    else
       s = "No"
    end if

  end function string_logical

  !> Compare two strings for equality.
  logical function equal (s,t)
    
    character(len=*), intent(in) :: s !< First string
    character(len=*), intent(in) :: t !< Second string

    equal = string(s) == string(t)

  end function equal

  !> Read a line from logical unit u, and return true if read was
  !> successful. Continuation with "\", skip blank lines and comments.
  !> If eofstop is true, raise error on EOF. If ucopy (integer) exists,
  !> write a copy of the output line to that logical unit, preceded
  !> by a prefix.
  function getline(u,oline,eofstop,ucopy)
 
    character(len=:), allocatable, intent(out) :: oline
    integer, intent(in) :: u
    logical, intent(in), optional :: eofstop
    integer, intent(in), optional :: ucopy
    logical :: getline
 
    integer :: i, lenu
    character(len=:), allocatable :: line
    logical :: notfirst, ok
    character(len=1) :: inquote
 
    character(len=2), parameter :: prfx = '%%'
 
    getline = .false.
    oline = ""
    notfirst = .false.
    do while (.true.)
       ! read the line
       ok = getline_raw(u,line)
 
       ! remove tabs
       do i = 1, len(line)
          if (line(i:i) == tab) line(i:i) = " "
       end do
 
       ! remove trailing blanks
       line = adjustl(line)
       lenu = len_trim(line)
 
       ! remove comments outside quotes
       inquote = ""
       do i = 1, lenu
          if (line(i:i) == "'".or.line(i:i) == '"') then
             if (len_trim(inquote) == 0) then
                inquote = line(i:i)
                cycle
             else if (line(i:i) == inquote) then
                inquote = ""
                cycle
             end if
          end if
          if (len_trim(inquote) > 0) cycle
          if (line(i:i) == "#") then
             line = line(1:i-1)
             exit
          end if
       end do
          
       ! exit if eof
       if (.not.ok) then
          if (present(eofstop)) then
             if (eofstop) call ferror("getline","unexpected end of file",faterr)
          end if
          return
       end if
 
       ! skip blank lines
       lenu = len_trim(line)
       if (lenu == 0) then
          if (notfirst) then
             exit
          else
             cycle
          end if
       end if
 
       ! continuation !AQUI problemas con PGI
#ifdef __PGI
       if (line(lenu:lenu) .ne. "\\") then
#else
       if (line(lenu:lenu) .ne. "\") then 
#endif
          oline = trim(oline) // " " // trim(line)
          oline = adjustl(oline)
          exit
       end if
       oline = trim(oline) // " " // line(1:lenu-1)
       notfirst = .true.
    end do
    if (present(ucopy)) then
       if (ucopy >= 0) &
          write (ucopy,'(A,1x,A/)') prfx, trim(oline)
    endif
    getline = .true.
 
  end function getline

  !> Read a line from logical unit u, and return true if read was
  !> successful. Don't process it to remove comments, etc.  If eofstop
  !> is true, raise error on EOF.
  function getline_raw(u,line,eofstop) result(ok)

    logical :: ok
    integer, intent(in) :: u
    character(len=:), allocatable, intent(out) :: line
    logical, intent(in), optional :: eofstop

    integer, parameter :: blocksize = 128

    character(len=:), allocatable :: aux
    integer :: ios, nread, i

    allocate(character(len=blocksize)::aux)
    do i = 1,blocksize
       aux(i:i) = " "
    end do
    line = ""
    ok = .true.
    do while(.true.)
       read(u,'(A)',advance="no",iostat=ios,size=nread) aux
       line = trim(line) // trim(aux)
       ok = .not.is_iostat_end(ios)
       if (is_iostat_eor(ios) .or. is_iostat_end(ios)) exit
    end do
    if (.not.ok.and.present(eofstop)) then
       if (eofstop) call ferror("getline_raw","unexpected end of file",faterr)
    end if

  end function getline_raw

  !> Convert a string to uppercase, portable.
  function upper(a)

    character(len=*), intent(in) :: a
    character(len=len(a)) :: upper

    character(len=*), parameter :: lo = 'abcdefghijklmnopqrstuvwxyz'
    character(len=*), parameter :: up = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    integer :: i, idx

    upper = a
    do i = 1,len(a)
       idx = index(lo,a(i:i))
       if (idx > 0) upper(i:i) = up(idx:idx)
    end do

  end function upper

  !> Convert a string to lowercase, portable.
  function lower(a)

    character(len=*), intent(in) :: a
    character(len=len(a)) :: lower

    character(len=*), parameter :: lo = 'abcdefghijklmnopqrstuvwxyz'
    character(len=*), parameter :: up = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    integer :: i, idx

    lower = a
    do i = 1,len(a)
       idx = index(up,a(i:i))
       if (idx > 0) lower(i:i) = lo(idx:idx)
    end do

  end function lower

  !> Get next word from line at lp. If successful, increase lp. A word
  !> is defined as any sequence on nonblanks.
  function getword(line,lp)

    character(len=:), allocatable :: getword !< Word read
    character(len=*), intent(in) :: line !< Input line
    integer, intent(inout) :: lp !< Pointer to position on input line, updated after reading.

    integer :: len0, wp, i

    len0 = len(line)
    getword = ""
    if (lp > len0) return

    i = lp
    do while (line(i:i) == " ")
       i = i+1
       if (i > len0) return
    end do

    wp = i
    do while (i<len0 .and. line(i:i)/=" ")
       i = i+1
    end do
    getword = line(wp:i)
    lp = i+1

  end function getword

  ! Same as getword, but convert to lowercase
  function lgetword (line,lp)

    character(len=*), intent(in) :: line !< Input line
    integer, intent(inout) :: lp !< Pointer to position on input line, updated after reading.
    character(len=:), allocatable :: lgetword

    lgetword = getword(line,lp)
    if (len_trim(lgetword) > 0) lgetword = lower(lgetword)

  endfunction lgetword

  !> True if the input is a letter.
  function isletter(a)

    character(len=1), intent(in) :: a
    logical :: isletter

    character(len=*), parameter :: letters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

    isletter = (index(letters,a) > 0)

  end function isletter

  !> True if the input is a digit.
  function isdigit(a)

    character(len=1), intent(in) :: a
    logical :: isdigit

    character(len=*), parameter :: digits = '0123456789'

    isdigit = (index(digits,a) > 0)

  end function isdigit

  !> Get an integer value from input text. If a valid integer is not
  !> found, then return .false.
  function isinteger (ival,line,lp)

    logical :: isinteger
    character(len=*), intent(in) :: line !< Input string
    integer, intent(inout) :: lp !< Pointer to current position on string
    integer, intent(inout) :: ival !< Integer value read, if not a valid integer is found, it is not modified

    integer :: i, len0
    logical :: readit

    isinteger = .false.
    len0 = len(line)
    if (lp > len0) return
    i = lp
    do while (line(i:i)==" " .or. line(i:i)==tab)
       i = i + 1
       if (i > len0) return
    end do
    if (line(i:i) == '+' .or. line(i:i) == '-') i = i + 1
    if (i > len0) return
    if (isdigit(line(i:i))) then 
       do while (isdigit(line(i:i)))
          if (i >= len0) exit
          i = i + 1
       end do
       readit = (i == len0)
       if (.not.readit) readit = (line(i:i) == " ")
       if (readit) then
          read (line(lp:i),*) ival
          lp = i + 1
          isinteger=.true.
       end if
    end if

  end function isinteger

  !> Get a real number from line and sets rval to it. If a valid real
  !> number is not found, isreal returns .false. and rval is not changed. 
  function isreal (rval, line, lp)

    character(len=*), intent(in) :: line !< Input string
    integer, intent(inout) :: lp !< Pointer to current position on string
    real(kind=8), intent(out) :: rval !< Real value read
    logical :: isreal

    integer :: len0
    logical :: matched

    character(len=*), parameter :: exponents = 'eEdDqQ+-'
    character(len=*), parameter :: lexponents = 'eEdDqQ'

    integer :: i

    ! initialize
    isreal = .false.
    len0 = len(line)
    if (lp > len0) return

    ! trailing blanks
    i = lp
    do while (line(i:i)==" " .or. line(i:i)==tab)
       i=i+1
       if (i > len0) return
    end do

    ! sign
    matched = .false.
    if (line(i:i) == '+' .or. line(i:i) == '-') i = i + 1 
    if (i > len0) return

    ! first sequence of digits
    if (isdigit(line(i:i))) then 
       matched = .true.
       do while (isdigit(line(i:i)))
          i=i+1
          if (i > len0) goto 10
       end do
    end if

    ! decimal point
    if (line(i:i) == '.') then 
       i=i+1
       if (i > len0) goto 10
    end if

    ! second sequence of digits
    if (isdigit(line(i:i))) then 
       matched = .true.
       do while (isdigit(line(i:i)))
          i=i+1
          if (i > len0) goto 10
       end do
    end if

    ! that did not look of a real number at all!
    if (.not.matched) return

    ! get optional exponent
    if (index(exponents,line(i:i)) > 0) then
       i = i+1
       if (i > len0) return

       if (line(i:i)=='+' .or. line(i:i)=='-') then
          if (index(lexponents,line(i-1:i-1)) > 0) then
             i = i+1 
             if (i > len0) return
          else
             return
          end if
       end if

       if (.not.isdigit(line(i:i))) return
       do while (isdigit(line(i:i)))
          i=i+1
          if (i > len0) goto 10
       end do

       if (i < len0 .and..not.line(i:i) == " ") return
    end if

10  continue
    read (line(lp:i-1),*) rval
    lp = i
    isreal = .true.

  end function isreal

  !> Get a expression from the line. Expressions start with " or '.
  !> If the read fails, then lp0 is unchanged, and isexpression is
  !> false.
  function isexpression(word, line, lp0)

    character(len=:), allocatable, intent(out) :: word !< Expression read
    character(len=*), intent(in) :: line !< Input line
    integer, intent(inout) :: lp0 !< Pointer to position on input line, updated after reading.
    logical :: isexpression !< True if the read was successful

    integer :: l, lp, lpi
    character(len=1) :: sclose

    isexpression = .false.
    lp = lp0
    l = len(line)
    if (lp > l) return
    do while (line(lp:lp)==" ")   !skip blanks
       lp = lp + 1
       if (lp > l) return
    end do

    if (line(lp:lp) /= '"' .and. line(lp:lp) /= "'") return
    sclose = line(lp:lp)
    lp = lp + 1
    lpi = lp
    if (lp > l) return
    do while (line(lp:lp) /= sclose)
       lp = lp + 1
       if (lp > l) return
    end do
    if (line(lp:lp) /= sclose) return

    isexpression = .true.
    word = line(lpi:lp-1)
    lp0 = lp + 1

  end function isexpression

  !> Get a expression from line. Expressions start with " or '. If no
  !> expression is available, then get the next word.
  function isexpression_or_word(word, line, lp0)

    character(len=:), allocatable, intent(out) :: word !< Word/expression read
    character(len=*), intent(in) :: line !< Input line
    integer, intent(inout) :: lp0 !< Pointer to position on input line, updated after reading.
    logical :: isexpression_or_word !< True if the read was successful

    character(len=:), allocatable :: word0

    ! is expression...
    isexpression_or_word = isexpression(word,line,lp0)
    if (isexpression_or_word) return
    
    ! ... or word?
    word0 = getword(line,lp0)
    if (len_trim(word0) > 0) then
       word = word0
       isexpression_or_word = .true.
    else
       isexpression_or_word = .false.
    end if

  end function isexpression_or_word

  !> Get a string from line.
  function isword(word, line, lp0)

    character(len=*), intent(out) :: word !< Word/expression read
    character(len=*), intent(in) :: line !< Input line
    integer, intent(inout) :: lp0 !< Pointer to position on input line, updated after reading.
    logical :: isword !< True if the read was successful

    character(len=:), allocatable :: word0

    ! ... or word?
    word0 = getword(line,lp0)
    if (len_trim(word0) > 0) then
       word = word0
       isword = .true.
    else
       isword = .false.
    end if

  end function isword

  !> Initialize file system and connect standard units. 
  !> Use this before falloc and fdealloc.
  subroutine ioinit ()

    use iso_fortran_env, only: error_unit, input_unit, output_unit

    uin = input_unit
    uout = output_unit
    uerr = error_unit
    udat = 101     

    alloc = .false.
    alloc(error_unit) = .true.
    alloc(input_unit) = .true.
    alloc(output_unit) = .true.
    alloc(101) = .true.

  end subroutine ioinit

  !> Open a file for reading. The argument form controls the
  !> formatting, and is passed directly to open(). If abspath is
  !> present, file in input is as an absolute path.
  function fopen_read(file,form,abspath0) result(lu)

    character(len=*), intent(in) :: file
    character(len=*), intent(in), optional :: form
    logical, intent(in), optional :: abspath0
    integer :: lu
    
    integer :: ios
    character(len=:), allocatable :: ofile
    logical :: abspath

    abspath = .false.
    ofile = trim(adjustl(filepath)) // dirsep // file
    if (file(1:1) == dirsep) abspath = .true.
    if (present(abspath0)) abspath = abspath0
    if (abspath) ofile = file

    lu = falloc()
    if (present(form)) then
       open(unit=lu,file=ofile,status='old',iostat=ios,form=form)
    else
       open(unit=lu,file=ofile,status='old',iostat=ios)
    end if
    if (ios /= 0) call ferror("fopen_read","error opening file: "//string(file),faterr)

  end function fopen_read

  !> Open a file for reading. The argument form controls the
  !> formatting, and is passed directly to open(). If abspath is
  !> present, file in input is as an absolute path.
  function fopen_write(file,form,abspath0) result(lu)

    character(len=*), intent(in) :: file
    character(len=*), intent(in), optional :: form
    logical, intent(in), optional :: abspath0
    integer :: lu
    
    integer :: ios
    character(len=:), allocatable :: ofile
    logical :: abspath

    abspath = .false.
    ofile = trim(adjustl(filepath)) // dirsep // file
    if (file(1:1) == dirsep) abspath = .true.
    if (present(abspath0)) abspath = abspath0
    if (abspath) ofile = file

    lu = falloc()
    if (present(form)) then
       open(unit=lu,file=ofile,status='unknown',iostat=ios,form=form)
    else
       open(unit=lu,file=ofile,status='unknown',iostat=ios)
    end if
    if (ios /= 0) call ferror("fopen_write","error opening file: "//string(file),faterr)

  end function fopen_write

  !> Open a file for appending
  function fopen_append(file,form,abspath0) result(lu)

    character(len=*), intent(in) :: file
    character(len=*), intent(in), optional :: form
    logical, intent(in), optional :: abspath0
    integer :: lu
    
    integer :: ios
    character(len=:), allocatable :: ofile
    logical :: abspath

    abspath = .false.
    ofile = trim(adjustl(filepath)) // dirsep // file
    if (file(1:1) == dirsep) abspath = .true.
    if (present(abspath0)) abspath = abspath0
    if (abspath) ofile = file

    lu = falloc()
    if (present(form)) then
       open(unit=lu,file=ofile,status='old',access='append',iostat=ios,form=form)
    else
       open(unit=lu,file=ofile,status='old',access='append',iostat=ios)
    end if
    if (ios /= 0) call ferror("fopen_append","error opening file: "//string(file),faterr)

  end function fopen_append

  !> Open a file for writing
  function fopen_scratch(form) result(lu)

    integer :: lu
    character(len=*), intent(in), optional :: form
    
    integer :: ios

    lu = falloc()
    if (present(form)) then
       open(unit=lu,status='scratch',form=form,iostat=ios)
    else
       open(unit=lu,status='scratch',form="unformatted",iostat=ios)
    end if
    if (ios /= 0) call ferror("fopen_scratch","error opening scratch file",faterr)

  end function fopen_scratch

  !> Close a file
  subroutine fclose(lu)

    integer, intent(in) :: lu

    close(lu)
    call fdealloc(lu)

  end subroutine fclose

  !> Allocate a logical unit
  function falloc()

    integer :: falloc

    do falloc = 1, size(alloc)
       if (.not.alloc(falloc)) exit
    end do
    if (alloc(falloc)) &
       call ferror("falloc","could not allocate logical unit",faterr)
    alloc(falloc) = .true.
    
  end function falloc

  !> Deallocate a logical unit
  subroutine fdealloc(u)

    integer, intent(in) :: u !< LU to deallocate

    if (u < 0 .or. u > size(alloc)) &
       call ferror("fdealloc","could not deallocate logical unit",faterr)
    alloc(u)=.false.

  end subroutine fdealloc

  !> Send an error message 'message' to stdout, coming from routine
  !> 'routine'. errortype is the error code (see mod_param.f90).
  subroutine ferror(routine,message,errortype,inputline,syntax)

    implicit none

    character(len=*), intent(in) :: routine !< Routine calling the error
    character(len=*), intent(in) :: message !< The message
    integer, intent(in) :: errortype !< Fatal, warning or info
    character(len=*), intent(in), optional :: inputline !< Previous-line input message
    logical, intent(in), optional :: syntax !< Is this a syntax error? (do not stop if interactive)

    character(len=20) :: chtype

    ! message styles shamelessly copied from abinit.
    if (errortype.eq.faterr) then
       chtype = '# ERROR'
    else if (errortype.eq.warning) then
       chtype = '# WARNING'
    else if (errortype.eq.noerr) then
       chtype = '# COMMENT'
    else
       chtype = '# UNKNOWN'
    endif
    !$omp critical (IO)
    if (present(inputline)) &
       write (uout,'("!error! ",A)') trim(inputline)
    write (uout,100) trim(chtype),trim(routine),trim(message)
    !$omp end critical (IO)
    if (errortype.eq.faterr) then
       !$omp critical (IO)
       write (uout,100) trim(chtype),trim(routine), trim(message)
       !$omp end critical (IO)
       if (present(syntax)) then
          if (syntax .and. interactive) return
       end if
       stop 1
    else if(errortype.eq.warning) then
       nwarns = nwarns + 1
    else if(errortype.eq.noerr) then
       ncomms = ncomms + 1
    endif

100 format (1x,A,"(",A,"): ",A)

  end subroutine ferror

  !> Interface to the timer routines of different computers. Returns
  !> machine seconds at the calling time. (private)
  subroutine tictac(mesg)

    character(len=*), intent(in) :: mesg !< Prefix message
    integer :: values(8)
    
    character(len=:), allocatable :: output

    call date_and_time(values=values)

    output = string(mesg)
    output = output // "--" // string(values(1))
    output = output // "/" // string(values(2)) 
    output = output // "/" // string(values(3))
    output = output // ", " // string(values(5),2,pad0=.true.) 
    output = output // ":" // string(values(6),2,pad0=.true.)
    output = output // ":" // string(values(7),2,pad0=.true.)  
    output = output // "." // string(values(8),3,pad0=.true.)
    write (uout,'(1x,A)') '# '//output

  end subroutine tictac

  subroutine flush_unit(unit_tobeflushed)

    implicit none    

    integer(kind=4), intent(in) :: unit_tobeflushed
    logical :: opnd
    inquire (unit=unit_tobeflushed,opened=opnd)
    if (opnd) call flush(unit_tobeflushed)

  end subroutine

  function find_free_unit()

    implicit none
  
    integer(kind=4), parameter :: max_unit = 240
    integer(kind=4) :: find_free_unit
    integer(kind=4) :: iunit
    logical :: opnd
    
    find_free_unit = -1
  
    do iunit = 100, max_unit
      inquire (unit=iunit,opened=opnd)
      if (.not.opnd) then
        find_free_unit = iunit
        return
      end if
    end do 
  
    call ferror('find_free_unit','free unit not found ?!?',faterr)
  
    return
  
  end function find_free_unit

  subroutine getdate(line)

    implicit none
    character(len=*), intent(inout) :: line
    character(len=mline) :: sdate
    integer(kind=4) :: values(8)

    call date_and_time(values=values)
    write (sdate,'(I4.4,".",I2.2,".",I2.2,", ",I2.2,":",I2.2,":",I2.2,".",I3.3)') values(1:3),values(5:8)
    line = trim(adjustl(sdate))

  end subroutine

  !> Start the clock
  subroutine start_clock()
    call system_clock(stime0)
    call cpu_time(ctime0)
  end subroutine start_clock
  
  !> Print the elapsed time to output
  subroutine print_clock()
    
    integer :: irate, imax
    real :: ctime
    integer :: stime
    character(len=:), allocatable :: sout

    call system_clock(count_rate=irate,count_max=imax)
    call cpu_time(ctime)
    call system_clock(stime)

    call sectotime(real(stime-stime0)/real(irate),s0=sout)
    write (uout,'("Elapsed wall time: ",A)') string(sout)
    call sectotime(ctime-ctime0,s0=sout)
    write (uout,'("Elapsed CPU time: ",A)') string(sout)

  end subroutine print_clock

  !> Seconds to days/hours/minutes/seconds format. x(1): seconds,
  !> x(2): minutes, x(3): hours, x(4): days.
  subroutine sectotime(s,x0,s0)
    real, intent(in) :: s
    integer, intent(out), optional :: x0(4)
    character(len=:), allocatable, optional :: s0

    integer :: i
    integer, parameter :: f(3) = (/60, 60, 24/)
    integer :: x(4)
    character(len=:), allocatable :: saux

    x = 0
    x(1) = floor(s)
    do i = 1, 3
       x(i+1) = x(i) / f(i)
       x(i) = x(i) - f(i) * x(i+1)
    end do

    if (present(x0)) x0 = x
    if (present(s0)) then
       s0 = ""
       if (x(4) > 0) then
          saux = trim(s0) // " " // string(x(4)) // "d"
          s0 = saux
       end if
       if (x(3) > 0) then
          saux = trim(s0) // " " // string(x(3)) // "h"
          s0 = saux
       end if
       if (x(2) > 0) then
          saux = trim(s0) // " " // string(x(2)) // "m"
          s0 = saux
       end if
       saux = trim(s0) // " " // string(x(1)) // "s"
       s0 = saux
    end if

  end subroutine sectotime

  character (len=4) function fourchar (i)
  
    integer :: i,i1,i2,i3,i4,i5,i6
    character(len=1) :: digs(0:9)
    data digs /'0','1','2','3','4','5','6','7','8','9'/

    i1 = i/1000
    i2 = mod(i,1000)
    i3 = i2/100
    i4 = mod(i2,100)
    i5 = i4/10
    i6 = mod(i4,10)
    fourchar = digs(i1)//digs(i3)//digs(i5)//digs(i6)

    return

  end function
      
  ! get command line arguments and argument count
  ! on return, argc will contain the number of arguments and
  ! argv will be a character array containing the arguments.
  subroutine getargs (argc, argv)

    implicit none
    integer :: argc
    character(len=*) :: argv(*)

    integer :: i
    character(len=mline) :: line         !local string for parsing
    integer :: length       !total length of command line

    argc=0
 10 continue
    call getarg (argc+1,line)
    length=len_trim(line)
    if (length.gt.0) then
      argc=argc+1
      if (argc.ge.3) return
      argv(argc)=line(1:length)
    endif
    if (length.gt.0) goto 10
    do i=1,argc
      call getarg(i,argv(i))
    end do
      
  end subroutine

end module mod_io
