#ifndef GAPI_H
#define GAPI_H

enum {mgrp = 200, ngtoh = 21};

int ncent_;
double *__restrict__ coords_;

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

#define EPS 1e-7
#define RHOEPS 1e-10

void cudasafe(int error, const char *message, const char *file, const int line);
void cerror(const char *text, const char *file, const int line);
extern "C" void gpu_info(void);
void init_nlm(void);

#endif
