module mod_gto
  
  use mod_prec, only: ip, rp
  implicit none
  private

  public :: gtogto

contains

  ! Overlap matrix between primitive Cartesian Gaussian Functions
  subroutine gtogto (sprim)

    use mod_param, only: pi
    use mod_mole, only: ncent_, xyz=>coords_
    use mod_basis, only: nlm, ityp_, oexp_, nprims_, ngroup_, nuexp_, nzexp_
    implicit none
    integer(kind=ip), parameter :: lamx = 12
 
    real(kind=rp), intent(out) :: sprim(nprims_,nprims_)
 
    real(kind=rp), dimension(ncent_,ncent_) :: ab2
    real(kind=rp) :: ax(1:3), bx(1:3), za, zb, p, pioverp, abaux
    real(kind=rp) :: prefactor, prod, xmu
    integer(kind=ip) :: i, j, k, l, m, ica, icb, nua, nub, n
    integer(kind=ip) :: itipa, itipb, la, lb, ka, kb, ma, mb
    ! ceabx() are the coefficients, except for the factor 
    ! EXP(-XMU*R_AB^2), where XMU=a*b/(a+b), that result from the 
    ! expansion of the product of two primitive cartesian Gaussian
    real(kind=rp) :: ceabx(-1:2*lamx,-1:lamx,-1:lamx,3)
 
    do ica = 1,ncent_
      do icb = 1,ica
        abaux = 0.0_rp
        do j = 1,3
          abaux = abaux + (xyz(j,ica)-xyz(j,icb))**2
        end do
        ab2(ica,icb) = abaux
        ab2(icb,ica) = abaux
      end do
    end do
   
    ! Compute the electronic molecular electrostatic potential.
    do ica = 1,ncent_
      ax(:) = xyz(:,ica)
      do ma = 1,ngroup_(ica)   
        nua = nuexp_(ica,1,ma)
        itipa = ityp_(nua)
        la = nlm(itipa,1) + nlm(itipa,2) + nlm(itipa,3)
        za = oexp_(nua)
        do icb = 1,ica
          bx(:) = xyz(:,icb)
          do mb = 1,ngroup_(icb)
            nub = nuexp_(icb,1,mb)
            itipb = ityp_(nub)
            lb = nlm(itipb,1) + nlm(itipb,2) + nlm(itipb,3)
            zb = oexp_(nub)
            p = za + zb
            xmu = za*zb/p
            prefactor = exp(-xmu*ab2(ica,icb))
            pioverp = pi/p
            pioverp = sqrt(pioverp*pioverp*pioverp)
            do j = 1,3
              call etijcalc (j,lamx,la,lb,ceabx,za,zb,ax(j),bx(j))
            end do
            ! Compute the target functions for all the products of
            ! of Gaussian primitives.
            do ka = 1,nzexp_(ica,ma)
              nua = nuexp_(ica,ka,ma)
              itipa = ityp_(nua)
              i = nlm(itipa,1)
              k = nlm(itipa,2)
              m = nlm(itipa,3)
              do kb = 1,nzexp_(icb,mb)
                nub = nuexp_(icb,kb,mb)
                if (nua.ge.nub) then
                  itipb = ityp_(nub)
                  j = nlm(itipb,1)
                  l = nlm(itipb,2)
                  n = nlm(itipb,3)
                  prod = ceabx(0,i,j,1)*ceabx(0,k,l,2)*ceabx(0,m,n,3)
                  sprim(nua,nub) = prod*pioverp*prefactor
                  sprim(nub,nua) = sprim(nua,nub)
                end if
              end do
            end do
          end do
        end do
      end do
    end do
   
  end subroutine

  subroutine etijcalc (m,lamx,la,lb,ce,a,b,ax,bx)
 
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

end module mod_gto
