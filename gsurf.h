#ifndef GSURF_H
#define GSURF_H

__constant__ int dncent_;
double *__restrict__ dcoords_;

enum {mgrp = 200, ngtoh = 21};
__constant__ int dnmo_;
__constant__ int dnprims_;
double *__restrict__ doexp_;
double *__restrict__ dmo_coeff_;
double *__restrict__ dmo_occ_;
int *__restrict__ dicen_;
int *__restrict__ dityp_;

int *__restrict__ dngroup_;
int *__restrict__ dnzexp_;
int *__restrict__ dnuexp_;
double *__restrict__ drcutte_;

int nlm[56][3];
__constant__ int dnlm[56][3];

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
__constant__ int dinuc_;        
__constant__ double depsiscp_;
__constant__ double depsroot_;
__constant__ int dntrial_;      
__constant__ int dnpang_;
__constant__ double drmaxsurf_;
__constant__ int dbackend_;
__constant__ double depsilon_;
__constant__ double dstep_;
__constant__ int dmstep_;
__constant__ int dsteeper_;
__constant__ double dxnuc_[3];
double *__restrict__ dxyzrho_;
double *__restrict__ drpru_;
double *__restrict__ dct_;
double *__restrict__ dst_;
double *__restrict__ dcp_;
double *__restrict__ dsp_;
int *dnlimsurf_;
double *drsurf_;

void init_nlm(void);

__device__ bool checkcp(const double *x, int *nuc,
             const int *ityp,            
             const double *oexp,         
             const int *ngroup,          
             const int *nzexp,           
             const int *nuexp,           
             const double *rcutte,       
             const double *mo_coeff,     
             const double *mo_occ,       
             const double *coords,
						 const double *xyzrho);

__device__ int odeint(double *ystart, double h1, double eps,
                      const int *ityp,            
                      const double *oexp,         
                      const int *ngroup,          
                      const int *nzexp,           
                      const int *nuexp,           
                      const double *rcutte,       
                      const double *mo_coeff,     
                      const double *mo_occ,       
                      const double *coords,
            					const double *xyzrho);

__device__ void rkqs(double *y, double *dydx, double *x, 
                     double htry, double eps,
	                   double *yscal, double *hnext,
                     const int *ityp,            
                     const double *oexp,         
                     const int *ngroup,          
                     const int *nzexp,           
                     const int *nuexp,           
                     const double *rcutte,       
                     const double *mo_coeff,     
                     const double *mo_occ,       
                     const double *coords);

__device__ void steeper_rkck(double *xpoint, double *grdt, double h0, double *xout, double *xerr,
                             const int *ityp,            
                             const double *oexp,         
                             const int *ngroup,          
                             const int *nzexp,           
                             const int *nuexp,           
                             const double *rcutte,       
                             const double *mo_coeff,     
                             const double *mo_occ,       
                             const double *coords);
__device__ void rho_grad(const int *ityp,        
                         const double *oexp,     
                         const int *ngroup,      
                         const int *nzexp,       
                         const int *nuexp,       
                         const double *rcutte,   
                         const double *mo_coeff, 
                         const double *mo_occ,   
                         const double *coords,   
                         const double *point,    
                         double *rho,            
                         double *grad,           
                         double *gradmod);       
__global__ void compute(const int *ityp,         
                        const double *oexp,      
                        const int *ngroup,       
                        const int *nzexp,        
                        const int *nuexp,        
                        const double *rcutte,    
                        const double *mo_coeff,  
                        const double *mo_occ,    
                        const double *coords,    
                        const double *xyzrho,    
                        const double *rpru,
				                const double *ct,         
				                const double *st,         
				                const double *cp,         
				                const double *sp,          
				                int *nlimsurf,         
				                double *rsurf);         

/*extern "C" void gsurf_driver(const int nmo, 
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
                  int *nlimsurf, double *rlimsurf);*/

#endif
