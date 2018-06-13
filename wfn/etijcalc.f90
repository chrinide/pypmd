subroutine etijcalc (m,lamx,la,lb,ce,a,b,ax,bx)
 
!-----This routine computes recursively the E_t^ij coefficients
!     resulting from the expansion of a cartesian gaussian product in 
!     terms of hermite gaussians. The expressions used come from the 
!     book "Molecular Electronic-Structure Theory" by T. Helgaker, 
!     P.Jorgensen, and J. Olsen. All the coefficients arise from the 
!     first one [E_0^00], taking also into account that E_t^ij=0 for 
!     t<0 or t>i+j.
!
!-----Details
!     
!     G_i (x,a,A_x) = x_A^i exp(-a x_A^2)  ;  x_A = x -  A_x
!     G_j (x,b,B_x) = x_B^j exp(-b x_B^2)  ;  x_B = x -  B_x
!
!     where
!     x   = X coordinate of an electron in a general reference frame
!     a   = Exponente of the first  Cartesian Gaussian
!     b   = Exponente of the second Cartesian Gaussian
!     A_X = X coordinate of nucleus A in a general reference frame
!     B_X = X coordinate of nucleus B in a general reference frame
!
!     Omega_ij      = G_i (x,a,A_x) * G_j (x,b,B_x)
!
!     Omega_ij      = SUM [t=0,i+j] K_ab^x E_t^ij Lambda_t (x,p,P_x)
!
!     where
!
!     Lambda_t (x,p,P_x) = [d/d P_x]^t exp(-p x_P^2)  ;  x_P = x -  P_x
!     p                  = a+b
!     P_x                = (a A_x + b B_x)/p
!     mu                 = ab/(a+b)
!     X_AB               = A_x - B_x
!     K_ab^x             = exp(-mu X_AB^2)
!
!
!     The gives E_t^ij for  (0 <= i <= la; 0 <= j <= lb; 0 <= t <= i+j)
!
!-----------------------------------------------------------------------
 
  use mod_prec, only: rp, ip
  use mod_io, only: faterr, ferror

  implicit none
  integer(kind=ip) :: la,lb,lab,i,j,t,i1,j1,t1,m,lamx
  real(kind=rp) :: ce(-1:2*lamx,-1:lamx,-1:lamx,3)
  real(kind=rp) :: a,b,ax,bx,p,ab,pa,pb,tp

  if (la.lt.0) call ferror ('etijcalc', 'fatal error, la < 0', faterr)
  if (lb.lt.0) call ferror ('etijcalc', 'fatal error, lb < 0', faterr)
  if (la.gt.lamx) call ferror ('etijcalc', 'fatal error, la > lamx', faterr)
  if (lb.gt.lamx) call ferror ('etijcalc', 'fatal error, lb > lamx', faterr)

  lab = la + lb
  ce(-1:lab,-1:la,-1:lb,m) = 0.0_rp
  ce(0,0,0,m) = 1.0_rp
  if (lab.eq.0) return
  p  = a + b
  ab = ax - bx
  pa = -b*ab/p
  pb = +a*ab/p
  tp = 1.0_rp/(2.0_rp*p)
  do i = 0,la
    i1 = i-1
    do j = 0,lb
      j1 = j-1
      do t = 1,i+j
        t1 = t-1
        ce(t,i,j,m) = tp*(i*ce(t1,i1,j,m) + j*ce(t1,i,j1,m))/real(t,rp)
      end do
      if (i.lt.la) ce(0,i+1,j,m) = pa*ce(0,i,j,m) + ce(1,i,j,m)
      if (j.lt.lb) ce(0,i,j+1,m) = pb*ce(0,i,j,m) + ce(1,i,j,m)
    end do
  end do 
 
end subroutine
