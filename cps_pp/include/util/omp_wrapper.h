#ifndef OMP_WRAPPER_H
#define OMP_WRAPPER_H
#include <config.h>
#ifdef USE_OMP
#warning "omp_wrapper.h: Using OpenMP"
#include <omp.h>
CPS_START_NAMESPACE
const int MAX_THREADS=64;
CPS_END_NAMESPACE
#define Pragma(x) _Pragma(#x)
#define OMP4(A) _Pragma( #A ) 
#else
#warning "omp_wrapper.h: NOT using OpenMP"
inline void omp_set_dynamic(bool val){}
inline int omp_get_num_threads(void) {return 1;}
inline int omp_get_max_threads(void) {return 1;}
inline int omp_get_thread_num(void) {return 0;}
inline void omp_set_num_threads(int n) {}
CPS_START_NAMESPACE
const int MAX_THREADS=1;
CPS_END_NAMESPACE
#endif
#endif
