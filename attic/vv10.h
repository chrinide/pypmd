#ifndef VV10_H
#define VV10_H

void driver_vv10(const int nmo, 
                 const int nprims,
                 const int *icen, 
                 const int *ityp,
                 const double *oexp, 
                 const int *ngroup,
            		 const int *nzexp,
            		 const int *nuexp,
            		 const double *rcutte,
            		 const double *mo_coeff,
            		 const double *mo_occ, 
            		 const int natm,
            		 const double *coords,
								 const int npoints,
								 const double *points,
								 double *rho,
								 double *gnorm);

inline void point_rho_grad_shell(const double *p, double *rho, double *grad, double *gradmod);

double vv10(const int n, const double coef_C,
                         const double coef_B,
                         const double *coords, 
                         const double *rho,
                         const double *weights,
                         const double *gnorm2);

#endif
