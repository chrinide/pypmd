! A segment of range <seg_range> will be divided into <subseg_num> 
! subsegments maximally uniformly. The length of each subsegment will 
! be returned in the array <subseg_sizes(1:subseg_num)>. Any two 
! subsegments will not differ in length by more than 1, longer subsegments 
! preceding the shorter ones.
subroutine divide_segment(seg_range,subseg_num,subseg_sizes,iloop)

  use mod_prec, only: ip
  use mod_io, only: faterr, ferror
  implicit none
  integer(kind=ip), intent(in) :: seg_range, subseg_num
  integer(kind=ip), intent(out) :: subseg_sizes(1:subseg_num)
  integer(kind=ip), intent(out) :: iloop(1:subseg_num,1:2)
  integer(kind=ip) :: i,m,n,j

  subseg_sizes = 0.0_ip
  if (seg_range.gt.0 .and. subseg_num.gt.0) then
    n = seg_range/subseg_num
    m = mod(seg_range,subseg_num)
    do i=1,m  
      subseg_sizes(i) = n + 1_ip
    end do
    do i = m+1,subseg_num
      subseg_sizes(i) = n
    end do
  else
    call ferror ('divide_segment', 'wrong input params', faterr)
  end if

  j = 0_ip
  iloop = 0_ip
  do i = 1,subseg_num
    iloop(i,1) = j + 1_ip
    j = j + subseg_sizes(i)
    iloop(i,2) = j
  end do

end subroutine divide_segment
!!!!!!$omp do schedule(static,ichunk)
!do jj = 1,npang,itile
!  do j = jj,min(jj+itile-1,npang)
!do j = 1,npang
