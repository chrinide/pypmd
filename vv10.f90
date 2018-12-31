real(kind=rp) function vv10_e(npoints, coef_c, coef_b, &
                              points, rho, weights, gnorm2) bind(c) 
                    
  use mod_param, only: pi
  use mod_prec, only: rp, ip
  implicit none

  integer(kind=ip), intent(in), value :: npoints
  real(kind=rp), intent(in), value :: coef_c, coef_b
  real(kind=rp), intent(in), dimension(3,npoints) :: points
  real(kind=rp), intent(in), dimension(npoints) :: rho
  real(kind=rp), intent(in), dimension(npoints) :: weights
  real(kind=rp), intent(in), dimension(npoints) :: gnorm2
  
  real(kind=rp), parameter :: const43 = 4.0/3.0*pi 

  integer(kind=ip) :: idx, jdx
  real(kind=rp) :: coef_beta, kappa_pref, g, gp
  real(kind=rp) :: px, py, pz, rho1, weight1, wp1, wg1, w01, kappa1, gamma1
  real(kind=rp) :: kernel, r, rho2, weight2, gamma2, wp2, wg2, w02, kappa2

  vv10_e = 0.0_rp
  coef_beta = 1.0/32.0*(3.0/(coef_B*coef_B))**(3.0/4.0)
  kappa_pref = coef_B*(1.5*pi)/ (9.0*pi)**(1.0/6.0)

  do idx = 1,npoints
    px = points(1,idx)
    py = points(2,idx)
    pz = points(3,idx)
    rho1 = rho(idx)
    weight1 = weights(idx)
    gamma1 = gnorm2(idx)
    wp1 = const43*rho1
    wg1 = coef_c*(gamma1/(rho1*rho1))**2
    w01 = sqrt(wg1 + wp1)
    kappa1 = rho1**(1.0/6.0)*kappa_pref
    kernel = 0.0_rp
    do jdx = 1,npoints
      r = (px-points(1,jdx))*(px-points(1,jdx)) 
      r = r + (py-points(2,jdx))*(py-points(2,jdx)) 
      r = r + (pz-points(3,jdx))*(pz-points(3,jdx)) 
      rho2 = rho(jdx)
      weight2 = weights(jdx)
      gamma2 = gnorm2(jdx)
      wp2 = const43*rho2 
      wg2 = coef_c*(gamma2/(rho2*rho2))**2
      w02 = sqrt(wg2 + wp2) 
      kappa2 = rho2**(1.0/6.0)*kappa_pref
      g = w01*r + kappa1 
      gp = w02*r + kappa2 
      kernel = kernel - 1.5*weight2*rho2/(g*gp*(g+gp)) 
    end do
    vv10_e = vv10_e + weight1*rho1*(coef_beta + 0.5*kernel) 
  end do

  return 

end function vv10_e
