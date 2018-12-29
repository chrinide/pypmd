#ifndef CSURF_H
#define CSURF_H

int ncent_;
double *__restrict__ coords_;

enum {mgrp = 200, ngtoh = 21};
int nmo_;
int nprims_;
double *__restrict__ oexp_;
double *__restrict__ mo_coeff_;
double *__restrict__ mo_occ_;
int *__restrict__ icen_;
int *__restrict__ ityp_;

int *__restrict__ ngroup_;
int *__restrict__ nzexp_;
int *__restrict__ nuexp_;
double *__restrict__ rcutte_;

int nlm[56][3];

// Surface info
#define EPS 1e-7
#define GRADEPS 1e-10
#define RHOEPS 1e-10
#define MINSTEP 1e-6
#define MAXSTEP 0.75
#define SAFETY 0.9
#define ENLARGE 1.2
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
#define TINY 1.0e-30
int inuc_;        
double epsiscp_;
double epsroot_;
int ntrial_;      
int npang_;
double rmaxsurf_;
int backend_;
double epsilon_;
double step_;
int mstep_;
int steeper_;
double xnuc_[3];
double *__restrict__ xyzrho_;
double *__restrict__ rpru_;
double *__restrict__ ct_;
double *__restrict__ st_;
double *__restrict__ cp_;
double *__restrict__ sp_;
int *nlimsurf_;
double *rsurf_;

void cerror(const char *text);

void init_nlm(void);

void print_basis(void);
void print_mole(void);

void rho_grad(const double *p, 
              double *rho, double *grad, double *gradmod);
bool checkcp(const double *x, int *nuc);

int odeint(double *ystart, double h1, double eps);
void rkqs(double *y, double *dydx, double *x, 
          double htry, double eps,
	        double *yscal, double *hnext);
void steeper_rkck(double *y, double *dydx, double h, double *yout, double *yerr);

void csurf_driver(const int nmo, 
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
                  const int npang,
                  const int inuc,
                  const double *xyzrho,
                  const double *ct,
                  const double *st,
                  const double *cp,
                  const double *sp,
                  const int backend,
                  const int ntrial,
                  const double epsiscp,
                  const double epsroot,
                  const double rmaxsurf,
                  const double epsilon,
                  const double *rpru,
                  const double step,
                  const int mstep,
                  double *nlimsurf, double *rlimsurf);
void surface(void);

#endif
