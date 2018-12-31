#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double vv10(const int n, const double coef_C,
                         const double coef_B,
                         const double *coords, 
                         const double *rho,
                         const double *weights,
                         const double *gnorm2)
{

  const double pi = 3.14159265358979323846264338328; 
  const double const43 = 4.0/3.0*pi;

  size_t idx, jdx;
  double vv10_e = 0.0;

  double coef_beta;
  coef_beta = 1.0/32.0*pow(3.0/(coef_B*coef_B),3.0/4.0);
  double kappa_pref;
  kappa_pref = coef_B*(1.5*pi)/pow(9.0*pi,1.0/6.0);

  for (idx=0; idx<n; idx++){
    double point1x = coords[idx*3+0];
    double point1y = coords[idx*3+1];
    double point1z = coords[idx*3+2];
    double rho1 = rho[idx];
    double weigth1 = weights[idx];
    double gamma1 = gnorm2[idx];
    double Wp1 = const43*rho1;
    double Wg1 = coef_C*pow(gamma1/(rho1*rho1),2);
    double W01 = sqrt(Wg1 + Wp1);
    double kappa1 = pow(rho1,1.0/6.0)*kappa_pref;
    double kernel = 0.0;
#pragma omp parallel default(none) reduction(+:kernel) \
    shared(point1x,point1y,point1z,coords,rho,weights,kappa_pref,gnorm2,kappa1,W01)
{
#pragma omp for nowait schedule(dynamic)
    for (jdx=0; jdx<n; jdx++){
      double R = (point1x-coords[jdx*3+0])*(point1x-coords[jdx*3+0]);
      R += (point1y-coords[jdx*3+1])*(point1y-coords[jdx*3+1]);
      R += (point1z-coords[jdx*3+2])*(point1z-coords[jdx*3+2]);
      double rho2 = rho[jdx];
      double weigth2 = weights[jdx];
      double gamma2 = gnorm2[jdx];
      double Wp2 = const43*rho2;
      double Wg2 = coef_C*pow(gamma2/(rho2*rho2),2);
      double W02 = sqrt(Wg2 + Wp2);
      double kappa2 = pow(rho2,1.0/6.0)*kappa_pref;
      double g = W01*R + kappa1;
      double gp = W02*R + kappa2;
      kernel += -1.5*weigth2*rho2/(g*gp*(g+gp));
    }
} 
    vv10_e += weigth1*rho1*(coef_beta + 0.5*kernel);
  }
 
  return vv10_e;

}
