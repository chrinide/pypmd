#ifndef VV10_GPU_H
#define VV10_GPU_H

extern "C"  
void driver_vv10_gpu(const int nmo, 
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

// GPU
__constant__ int dnprims_;
__constant__ int dnmo_;
double *dbrho;
double *dbcoords;
double *doexp_;
double *dmo_occ_;
double *dcoords_;
double *dmo_coeff_;
int *dicen_;
int *dityp_;
__constant__ int dnatm_;
__constant__ int dnlm_[56][3];

#endif
