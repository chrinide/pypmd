  subroutine gtogto (sprim)
!
! Overlap matrix between primitive Cartesian Gaussian Functions
!
  use mod_prec, only: rp, ip
  use mod_param, only: pi
  use mod_wfn, only: xyz, nuexp, ncent, nlm, nprims, ngroup, &
                     ityp, oexp, nzexp  
  implicit none
  integer(kind=ip), parameter :: lamx = 12
!
  real(kind=rp), intent(out) :: sprim(nprims,nprims)
!
  real(kind=rp), dimension(ncent,ncent) :: ab2
  real(kind=rp) :: ax(1:3), bx(1:3), za, zb, p, pioverp, abaux
  real(kind=rp) :: prefactor, prod, xmu
  integer(kind=ip) :: i, j, k, l, m, ica, icb, nua, nub, n
  integer(kind=ip) :: itipa, itipb, la, lb, ka, kb, ma, mb
!
! ceabx() are the coefficients, except for the factor 
! EXP(-XMU*R_AB^2), where XMU=a*b/(a+b), that result from the 
! expansion of the product of two primitive cartesian Gaussian
!
  real(kind=rp) :: ceabx(-1:2*lamx,-1:lamx,-1:lamx,3)
!
  do ica = 1,ncent
    do icb = 1,ica
      abaux = 0.0_rp
      do j = 1,3
        abaux = abaux + (xyz(ica,j)-xyz(icb,j))**2
      end do
      ab2(ica,icb) = abaux
      ab2(icb,ica) = abaux
    end do
  end do
!
! Compute the electronic molecular electrostatic potential.
!
  do ica = 1,ncent
    ax(1:3) = xyz(ica,1:3)
    do ma = 1,ngroup(ica)   
      nua = nuexp(ica,ma,1)
      itipa = ityp(nua)
      la = nlm(itipa,1) + nlm(itipa,2) + nlm(itipa,3)
      za = oexp(nua)
      do icb = 1,ica
        bx(1:3) = xyz(icb,1:3)
        do mb = 1,ngroup(icb)
          nub = nuexp(icb,mb,1)
          itipb = ityp(nub)
          lb = nlm(itipb,1) + nlm(itipb,2) + nlm(itipb,3)
          zb = oexp(nub)
          p = za + zb
          xmu = za*zb/p
          prefactor = exp(-xmu*ab2(ica,icb))
          pioverp = pi/p
          pioverp = sqrt(pioverp*pioverp*pioverp)
          do j = 1,3
            call etijcalc (j,lamx,la,lb,ceabx,za,zb,ax(j),bx(j))
          end do
!
!.........Compute the target functions for all the products of
!         of Gaussian primitives.
!
          do ka = 1,nzexp(ica,ma)
            nua = nuexp(ica,ma,ka)
            itipa = ityp(nua)
            i = nlm(itipa,1)
            k = nlm(itipa,2)
            m = nlm(itipa,3)
            do kb = 1,nzexp(icb,mb)
              nub = nuexp(icb,mb,kb)
              if (nua.ge.nub) then
                itipb = ityp(nub)
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
!
  end subroutine
