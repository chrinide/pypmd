#ifndef GAPI_H
#define GAPI_H

enum {mgrp = 200, ngtoh = 21};

__constant__ int dnprims;
__constant__ int dnmo;
__constant__ int dncent;
__constant__ int dnpoints;

int nlm[56][3];
__constant__ int dnlm[56][3];

#define EPS 1e-7
#define RHOEPS 1e-10

void cudasafe(int error, const char *message, const char *file, const int line);
void cerror(const char *text, const char *file, const int line);
extern "C" void gpu_info(void);
void init_nlm(void);

__device__ __forceinline__ void eval_rho_grad(const int *ityp,        
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
                        const double *points,    
                        double *output);         
#endif
