#ifdef _OPENMP
#include <omp.h>
int get_omp_threads() {
  return omp_get_max_threads();
}
int set_omp_threads(int n) {
  omp_set_num_threads(n);
  return n;
}
#else
// mimic omp_get_max_threads omp_set_num_threads function of libgomp
int get_omp_threads() {return 1;}
int set_omp_threads(int n) {return 0;}
#endif
