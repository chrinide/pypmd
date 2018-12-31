#ifndef FIELDS_H
#define FIELDS_H

inline void point_rho_grad(const double *p, double *rho, double *grad, double *gradmod);
inline void point_rho_grad_shell(const double *p, double *rho, double *grad, double *gradmod);
inline void point_rho_shell(const double *p, double *rho);

#endif
