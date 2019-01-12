subroutine atomic(npoints,z,xyz,points,output) bind(c)

  use mod_prec, only: rp, ip
  use mod_atomic, only: init_atomic, ftot
  use mod_math, only: spline
  use mod_memory, only: alloc, free
  use mod_io, only: ferror, faterr
  implicit none

  integer(kind=ip), intent(in), value :: npoints
  real(kind=rp), intent(in), dimension(3,npoints) :: points
  real(kind=rp), intent(out), dimension(npoints) :: output
  integer(kind=ip), intent(in), value :: z
  real(kind=rp), intent(in), dimension(3) :: xyz

  integer(kind=ip) :: intq, j, ndata, kk 
  real(kind=rp) :: h, q, x(3), r, dq, arho, rdata, rmid
  real(kind=rp), allocatable, dimension(:) :: a1, b1, c1, f1

  real(kind=rp), parameter :: third = 1.0_rp/3.0_rp

  call init_atomic ()
  rmid = 1.0_rp / z**third

  if (z <= 2) then
    ndata = 200
  elseif (z <= 10) then
    ndata = 400
  elseif (z <= 18) then
    ndata = 600
  elseif (z <= 36) then
    ndata = 800
  elseif (z <= 54) then
    ndata = 1000
  elseif (z <= 86) then
    ndata = 1200
  elseif (z <= 94) then
    ndata = 1400
  else
    call ferror('atomic','atomic number out of range',faterr)
  end if

  call alloc ('atomic', 'f1', f1, ndata, index_=0)
  call alloc ('atomic', 'a1', a1, ndata, index_=0)
  call alloc ('atomic', 'b1', b1, ndata, index_=0)
  call alloc ('atomic', 'c1', c1, ndata, index_=0)

  h = 1.0_rp / (ndata+1)
  f1(0) = 0.0_rp
  do j = 1,ndata
    q = h*j
    rdata = rmid*q/(1.0_rp-q)
    f1(j) = ftot(j,z)
  end do
  call spline (h,f1,a1,b1,c1,ndata,0.0_rp)

  do kk = 1,npoints
    x = points(:,kk) - xyz(:)
    r = sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
    q = r/(r + rmid)
    intq = int((ndata+1)*q)
    dq = q - intq*h
    arho = abs((f1(intq)+dq*(a1(intq)+dq*(b1(intq)+dq*c1(intq)))))/(r*r)
    output(kk) = arho
  end do

  call free ('atomic', 'f1', f1)
  call free ('atomic', 'f1', a1)
  call free ('atomic', 'f1', b1)
  call free ('atomic', 'f1', c1)

end subroutine atomic
