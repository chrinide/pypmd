#ifndef DRIVER_H
#define DRIVER_H

// CPU
int nprim_;
int nmo_;
double *oexp_;
double *coords_;
double *mo_coeff_;
double *mo_occ_;
int *icen_;
int *ityp_;
int natm_;
int nlm[56][3];
void init_nlm();

// Functions
extern "C" void driver(const int nmo, const int nprim, const int natm, 
						const int npoints, const int *icen, const int *ityp, 
						const double *oexp, const double *mo_coeff, 
						const double *mo_occ, const double *coords, 
            const double *bcoords, double *brho);

double point_rho(double *p);

// GPU
__constant__ int dnprim_;
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

__device__ __forceinline__ double rho(double *point, double *coeff,
						          double *occ, double *oxp, int *cen, int *typ,
                      double *xyz); 

__global__ void compute(double *out, double *points, double *coeff, 
						            double *occ, double *oxp, int *cen, int *typ,
                        double *xyz, const int npoints);

#endif
