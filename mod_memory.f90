module mod_memory

  use mod_prec, only: rp, i4=>ip, i8=>size_t
  use mod_io, only: faterr, ferror
  implicit none
  public

  integer(kind=i8) :: memavail = 0
  integer(kind=i8) :: mallocated = 0
  integer(kind=i8) :: max_allocated = 0
  
  interface free
    module procedure ifree
    module procedure ifree2
    module procedure ifree3
    module procedure dfree
    module procedure dfree2
    module procedure dfree3
    module procedure lfree
    module procedure lfree2
    module procedure lfree3
  end interface free
  
  interface alloc
    module procedure ialloc
    module procedure ialloc2
    module procedure ialloc3
    module procedure dalloc
    module procedure dalloc2
    module procedure dalloc3
    module procedure lalloc
    module procedure lalloc2
    module procedure lalloc3
  end interface alloc

contains

  subroutine ialloc(file__, array__, array, size_t, debug, index_)

    implicit none
    integer(kind=i4), intent(in) :: size_t
    character(len=*), intent(in) :: file__, array__
    integer(kind=i4), allocatable, dimension(:), intent(out) :: array
    logical, optional, intent(in) :: debug
    integer(kind=i4), optional, intent(in) :: index_
 
    integer(kind=i4) :: ier
 
    if (allocated(array)) then
      call ferror(trim(file__), 'array '//trim(array__)//' already in use', faterr)
    else
      if (present(index_)) then
        allocate (array(index_:size_t), stat=ier) 
      else
        allocate (array(size_t), stat=ier) 
      end if
      if (ier.ne.0) then
        call ferror(trim(file__), 'error allocating '//trim(array__), faterr)
      end if
      array = 0_i4
    end if
    mallocated = mallocated + size_t 

    if (present(debug) .and. debug) then
    end if

  end subroutine

  subroutine dalloc(file__, array__, array, size_t, debug, index_)

    implicit none
    integer(kind=i4), intent(in) :: size_t
    character(len=*), intent(in) :: file__, array__
    real(kind=rp), allocatable, dimension(:), intent(out) :: array
    logical, optional, intent(in) :: debug
    integer(kind=i4), optional, intent(in) :: index_
 
    integer(kind=i4) :: ier
 
    if (allocated(array)) then
      call ferror(trim(file__), 'array '//trim(array__)//' already in use', faterr)
    else
      if (present(index_)) then
        allocate (array(index_:size_t), stat=ier) 
      else
        allocate (array(size_t), stat=ier) 
      end if
      if (ier.ne.0) then
        call ferror(trim(file__), 'error allocating '//trim(array__), faterr)
      end if
      array = 0.0_rp
    end if
    mallocated = mallocated + size_t 

    if (present(debug) .and. debug) then
    end if

  end subroutine

  subroutine lalloc(file__, array__, array, size_t, debug)

    implicit none
    integer(kind=i4), intent(in) :: size_t
    character(len=*), intent(in) :: file__, array__
    logical, allocatable, dimension(:) :: array
    logical, optional :: debug
 
    integer(kind=i4) :: ier
 
    if (allocated(array)) then
      call ferror(trim(file__), 'array '//trim(array__)//' already in use', faterr)
    else
      allocate (array(size_t), stat=ier) 
      if (ier.ne.0) then
        call ferror(trim(file__), 'error allocating '//trim(array__), faterr)
      end if
    end if
    mallocated = mallocated + size_t 

    if (present(debug) .and. debug) then
    end if

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ialloc2(file__, array__, array, size_t1, size_t2, debug)

    implicit none
    integer(kind=i4), intent(in) :: size_t1, size_t2
    character(len=*), intent(in) :: file__, array__
    integer(kind=i4), allocatable, dimension(:,:) :: array
    logical, optional :: debug
 
    integer(kind=i4) :: ier
 
    if (allocated(array)) then
      call ferror(trim(file__), 'array '//trim(array__)//' already in use', faterr)
    else
      allocate (array(size_t1,size_t2), stat=ier) 
      if (ier.ne.0) then
        call ferror(trim(file__), 'error allocating '//trim(array__), faterr)
      end if
      array = 0_i4
    end if
    mallocated = mallocated + size_t1*size_t2

    if (present(debug) .and. debug) then
    end if

  end subroutine

  subroutine dalloc2(file__, array__, array, size_t1, size_t2, debug, index_)

    implicit none
    integer(kind=i4), intent(in) :: size_t1, size_t2
    character(len=*), intent(in) :: file__, array__
    real(kind=rp), allocatable, dimension(:,:), intent(out) :: array
    logical, optional, intent(in) :: debug
    integer(kind=i4), optional, intent(in) :: index_
 
    integer(kind=i4) :: ier
 
    if (allocated(array)) then
      call ferror(trim(file__), 'array '//trim(array__)//' already in use', faterr)
    else
      if (present(index_)) then
        allocate (array(index_:size_t1,index_:size_t2), stat=ier) 
      else
        allocate (array(size_t1,size_t2), stat=ier) 
      end if
      if (ier.ne.0) then
        call ferror(trim(file__), 'error allocating '//trim(array__), faterr)
      end if
      array = 0.0_rp
    end if
    mallocated = mallocated + size_t1*size_t2

    if (present(debug) .and. debug) then
    end if

  end subroutine

  subroutine lalloc2(file__, array__, array, size_t1, size_t2, debug)

    implicit none
    integer(kind=i4), intent(in) :: size_t1, size_t2
    character(len=*), intent(in) :: file__, array__
    logical, allocatable, dimension(:,:) :: array
    logical, optional :: debug
 
    integer(kind=i4) :: ier
 
    if (allocated(array)) then
      call ferror(trim(file__), 'array '//trim(array__)//' already in use', faterr)
    else
      allocate (array(size_t1,size_t2), stat=ier) 
      if (ier.ne.0) then
        call ferror(trim(file__), 'error allocating '//trim(array__), faterr)
      end if
    end if
    mallocated = mallocated + size_t1*size_t2

    if (present(debug) .and. debug) then
    end if

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ialloc3(file__, array__, array, size_t1, size_t2, size_t3, debug)

    implicit none
    integer(kind=i4), intent(in) :: size_t1, size_t2, size_t3
    character(len=*), intent(in) :: file__, array__
    integer(kind=i4), allocatable, dimension(:,:,:) :: array
    logical, optional :: debug
 
    integer(kind=i4) :: ier
 
    if (allocated(array)) then
      call ferror(trim(file__), 'array '//trim(array__)//' already in use', faterr)
    else
      allocate (array(size_t1,size_t2,size_t3), stat=ier) 
      if (ier.ne.0) then
        call ferror(trim(file__), 'error allocating '//trim(array__), faterr)
      end if
      array = 0_i4
    end if
    mallocated = mallocated + size_t1*size_t2*size_t3

    if (present(debug) .and. debug) then
    end if

  end subroutine

  subroutine dalloc3(file__, array__, array, size_t1, size_t2, size_t3, debug)

    implicit none
    integer(kind=i4), intent(in) :: size_t1, size_t2, size_t3
    character(len=*), intent(in) :: file__, array__
    real(kind=rp), allocatable, dimension(:,:,:) :: array
    logical, optional :: debug
 
    integer(kind=i4) :: ier
 
    if (allocated(array)) then
      call ferror(trim(file__), 'array '//trim(array__)//' already in use', faterr)
    else
      allocate (array(size_t1,size_t2,size_t3), stat=ier) 
      if (ier.ne.0) then
        call ferror(trim(file__), 'error allocating '//trim(array__), faterr)
      end if
      array = 0.0_rp
    end if
    mallocated = mallocated + size_t1*size_t2*size_t3

    if (present(debug) .and. debug) then
    end if

  end subroutine

  subroutine lalloc3(file__, array__, array, size_t1, size_t2, size_t3, debug)

    implicit none
    integer(kind=i4), intent(in) :: size_t1, size_t2, size_t3
    character(len=*), intent(in) :: file__, array__
    logical, allocatable, dimension(:,:,:) :: array
    logical, optional :: debug
 
    integer(kind=i4) :: ier
 
    if (allocated(array)) then
      call ferror(trim(file__), 'array '//trim(array__)//' already in use', faterr)
    else
      allocate (array(size_t1,size_t2,size_t3), stat=ier) 
      if (ier.ne.0) then
        call ferror(trim(file__), 'error allocating '//trim(array__), faterr)
      end if
    end if
    mallocated = mallocated + size_t1*size_t2*size_t3

    if (present(debug) .and. debug) then
    end if

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ifree(file__, array__, array)
 
    implicit none
    character(len=*), intent(in) :: file__, array__
    integer(kind=i4), allocatable, dimension(:) :: array
 
    integer(kind=i4) :: ier
 
    if (.not.allocated(array)) then
      call ferror(trim(file__), 'array '//trim(array__)//' already deallocated', faterr)
    else
      deallocate (array, stat=ier) 
      if (ier.ne.0) then
        call ferror(trim(file__), 'error deallocating '//trim(array__), faterr)
      end if
    end if
 
  end subroutine

  subroutine ifree2(file__, array__, array)
 
    implicit none
    character(len=*), intent(in) :: file__, array__
    integer(kind=i4), allocatable, dimension(:,:) :: array
 
    integer(kind=i4) :: ier
 
    if (.not.allocated(array)) then
      call ferror(trim(file__), 'array '//trim(array__)//' already deallocated', faterr)
    else
      deallocate (array, stat=ier) 
      if (ier.ne.0) then
        call ferror(trim(file__), 'error deallocating '//trim(array__), faterr)
      end if
    end if
 
  end subroutine

  subroutine ifree3(file__, array__, array)
 
    implicit none
    character(len=*), intent(in) :: file__, array__
    integer(kind=i4), allocatable, dimension(:,:,:) :: array
 
    integer(kind=i4) :: ier
 
    if (.not.allocated(array)) then
      call ferror(trim(file__), 'array '//trim(array__)//' already deallocated', faterr)
    else
      deallocate (array, stat=ier) 
      if (ier.ne.0) then
        call ferror(trim(file__), 'error deallocating '//trim(array__), faterr)
      end if
    end if
 
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dfree(file__, array__, array)
 
    implicit none
    character(len=*), intent(in) :: file__, array__
    real(kind=rp), allocatable, dimension(:) :: array
    
    integer(kind=i4) :: ier
 
    if (.not.allocated(array)) then
      call ferror(trim(file__), 'array '//trim(array__)//' already deallocated', faterr)
    else
      deallocate (array, stat=ier) 
      if (ier.ne.0) then
        call ferror(trim(file__), 'error deallocating '//trim(array__), faterr)
      end if
    end if
 
  end subroutine

  subroutine dfree2(file__, array__, array)
 
    implicit none
    character(len=*), intent(in) :: file__, array__
    real(kind=rp), allocatable, dimension(:,:) :: array
    
    integer(kind=i4) :: ier
 
    if (.not.allocated(array)) then
      call ferror(trim(file__), 'array '//trim(array__)//' already deallocated', faterr)
    else
      deallocate (array, stat=ier) 
      if (ier.ne.0) then
        call ferror(trim(file__), 'error deallocating '//trim(array__), faterr)
      end if
    end if
 
  end subroutine

  subroutine dfree3(file__, array__, array)
 
    implicit none
    character(len=*), intent(in) :: file__, array__
    real(kind=rp), allocatable, dimension(:,:,:) :: array
    
    integer(kind=i4) :: ier
 
    if (.not.allocated(array)) then
      call ferror(trim(file__), 'array '//trim(array__)//' already deallocated', faterr)
    else
      deallocate (array, stat=ier) 
      if (ier.ne.0) then
        call ferror(trim(file__), 'error deallocating '//trim(array__), faterr)
      end if
    end if
 
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lfree(file__, array__, array)
 
    implicit none
    character(len=*), intent(in) :: file__, array__
    logical, allocatable, dimension(:) :: array
    
    integer(kind=i4) :: ier
 
    if (.not.allocated(array)) then
      call ferror(trim(file__), 'array '//trim(array__)//' already deallocated', faterr)
    else
      deallocate (array, stat=ier) 
      if (ier.ne.0) then
        call ferror(trim(file__), 'error deallocating '//trim(array__), faterr)
      end if
    end if
 
  end subroutine

  subroutine lfree2(file__, array__, array)
 
    implicit none
    character(len=*), intent(in) :: file__, array__
    logical, allocatable, dimension(:,:) :: array
    
    integer(kind=i4) :: ier
 
    if (.not.allocated(array)) then
      call ferror(trim(file__), 'array '//trim(array__)//' already deallocated', faterr)
    else
      deallocate (array, stat=ier) 
      if (ier.ne.0) then
        call ferror(trim(file__), 'error deallocating '//trim(array__), faterr)
      end if
    end if
 
  end subroutine

  subroutine lfree3(file__, array__, array)
 
    implicit none
    character(len=*), intent(in) :: file__, array__
    logical, allocatable, dimension(:,:,:) :: array
    
    integer(kind=i4) :: ier
 
    if (.not.allocated(array)) then
      call ferror(trim(file__), 'array '//trim(array__)//' already deallocated', faterr)
    else
      deallocate (array, stat=ier) 
      if (ier.ne.0) then
        call ferror(trim(file__), 'error deallocating '//trim(array__), faterr)
      end if
    end if
 
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine setmemavail(bytes)
  
    implicit none
    integer, intent(in) :: bytes

    memavail = bytes

  end subroutine

  subroutine memstats()

    implicit none

    max_allocated = mallocated

  end subroutine memstats

end module mod_memory
