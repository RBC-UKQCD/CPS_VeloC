#if (defined USE_SSE)||(defined SSE_TO_C)
/*
 *  Wilson Dslahsh  2010/03/  Taku Izubuchi
 *
 *    Acknowledgements:
           Winners of the Intel tuning content 2006, especially
           NUKADA Akira is acknowledged.
           Sinya Aoki is acknowledged for the MPI transofrmation part.
 */

#include <config.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <util/vector.h>
#include <util/wilson.h>
#include <util/dirac_op.h>
#include <util/error.h>
#include <util/gjp.h>
//#include <comms/scu.h>
#include <util/data_types.h>
//#include <omp.h>
#include <util/omp_wrapper.h>
//inline int omp_get_num_threads(void) {return 1;}
#ifdef USE_SSE
#include <pmmintrin.h>
#endif

#include "sse-defs.h"
#include "sse-subs.h"
#ifdef SSE_TO_C
#else
#define SSE_C_FLOAT Float
#endif


//#define SSE_TO_C

//#ifdef SSE_TO_C
#if 0
//#else
inline static void TOUCH(const double*a, int n)
{
  int i;
  for(i=0;i<n;i+=64/sizeof(double)){
    __builtin_prefetch (((const char*)(a+i)), 0, (_MM_HINT_T0)) ;
  }
}
#endif

// TIZB, restore later !
#pragma warning disable 592


#if defined(_OPENMP) && defined(__linux__)
/*
 *  Set CPU affinity for 2-way dual-core system.
 *
 *  Proc. ID   Phys. Proc.   Core          Thread ID
 *     0           0           0     <---     0
 *     1           1           0     <---     2
 *     2           0           1     <---     1
 *     3           1           1     <---     3
 */
#include <sys/types.h>
#include <linux/unistd.h>
#include <errno.h>

#include <sys/types.h>
#include <unistd.h>
#include <sched.h>
#include <linux/unistd.h>
#include <errno.h>

//_syscall0(pid_t, gettid)
pid_t gettid( void )
{
	return syscall( __NR_gettid );
}
void cpubind(void)
{
#pragma omp parallel
  {
    cpu_set_t *cs;
    int p;
    cs = (cpu_set_t *)malloc(sizeof(cpu_set_t));
    CPU_ZERO(cs);
    p = omp_get_thread_num();
    ::printf("cpubind thread_num=%d\n", p);
    p = ((p & 1) << 1) | ((p & 2) >> 1);
    CPU_SET(p, cs);
    sched_setaffinity(gettid(), sizeof(cpu_set_t), cs);
  }
}
#endif




// #define AXG_WILSON

//#define DEBUG_PRINT_DAG0
//#define DEBUG_PRINT_DAG1
//#define DEBUG_NOEDGE
#define DEBUG_MFLOPS (1*1);


#undef PROFILE
//#define PROF_printf(...)  printf(__VA_ARGS__);
//#define PROF_printf(...)  if (!UniqueID()) printf(__VA_ARGS__);
#define PROF_printf(...)  {}




// USE register  t00-t11 to add into
//#define ADD2REG
// use HERN2 (doesn't use the temporal registers _a,_b,_c,_d)
//#define USE_HERN2


#ifndef AXG_WILSON
#define GAUGE_SIZE 72 /* re/im * colors*colors * dim */
#else
#define GAUGE_SIZE 54 /* re/im * colors*colors * dim */
#endif

CPS_START_NAMESPACE

#if 0
#define U(r,row,col,d)  *(u+(r+2*(row+3*(col+3*d))))
#define PSI(r,c,s)      *(psi +(r+2*(c+3*s)))
//! As above, but the vector is called chi
#define CHI(r,c,s)      *(chi +(r+2*(c+3*s)))
#define TMP(r,c,s)      *(tmp +(r+2*(c+3*s))) 
#define TMP1(r,c,s)     *(tmp1+(r+2*(c+3*s))) 
#define TMP2(r,c,s)     *(tmp2+(r+2*(c+3*s))) 
#define TMP3(r,c,s)     *(tmp3+(r+2*(c+3*s))) 
#define TMP4(r,c,s)     *(tmp4+(r+2*(c+3*s))) 
#define TMP5(r,c,s)     *(tmp5+(r+2*(c+3*s))) 
#define TMP6(r,c,s)     *(tmp6+(r+2*(c+3*s))) 
#define TMP7(r,c,s)     *(tmp7+(r+2*(c+3*s))) 
#define TMP8(r,c,s)     *(tmp8+(r+2*(c+3*s))) 
#define FBUF(r,c,s)     *(fbuf+(r+2*(c+3*s))) 
#endif


Float dclock(void);


extern "C" {
void wilson_dslash_bnd_dag0(
                   IFloat *chi_p_f, 
		   IFloat *u_p_f, 
		   IFloat *psi_p_f, 
		   int cb,
		   Wilson *wilson_p);
void wilson_dslash_bnd_dag1(
                   IFloat *chi_p_f, 
		   IFloat *u_p_f, 
		   IFloat *psi_p_f, 
		   int cb,
		   Wilson *wilson_p);
  
  
void wilson_dslash_blk_dag0(
                   IFloat *chi_p_f, 
		   IFloat *u_p_f, 
		   IFloat *psi_p_f, 
		   int cb,
		   Wilson *wilson_p);

void wilson_dslash_blk_dag1(
                   IFloat *chi_p_f, 
		   IFloat *u_p_f, 
		   IFloat *psi_p_f, 
		   int cb,
		   Wilson *wilson_p);

void wilson_dslash(
                   IFloat *chi_p_f, 
		   IFloat *u_p_f, 
		   IFloat *psi_p_f, 
		   int cb,
		   int dag,
		   Wilson *wilson_p)
{
  if(GJP.Gparity()){
    if(!UniqueID()) 
      printf("sse wilson dslash not implemented. Be my guest, I guarantee you will not have much fun doing it.\n");
    exit(-1);
  }

  const int   lx = wilson_p->ptr[0];
  const int   ly = wilson_p->ptr[1];
  const int   lz = wilson_p->ptr[2];
  const int   lt = wilson_p->ptr[3];
  const int   vol = wilson_p->vol[0];


#ifdef PROFILE
  int NITR = wilson_p->NITR;
  const Float MultFlops=wilson_p->MultFlops;
  const Float MultFlops_bnd=wilson_p->MultFlops_bnd;
  const Float MultFlops_blk=wilson_p->MultFlops_blk;
  const Float CPU_GHZ = wilson_p->CPU_GHZ;
  const int  num_threads= wilson_p->num_threads;
  Float time_bnd(0), time_blk(0), time_comm1(0), time_comm2(0);
#endif

#ifdef PROFILE
  for(int itr=0;itr<NITR;++itr){
    Float time = -dclock();
#endif  

    /* Do the boundary parts */

    if(dag == 0 ) wilson_dslash_bnd_dag0(chi_p_f, u_p_f, psi_p_f, cb, wilson_p);
    else          wilson_dslash_bnd_dag1(chi_p_f, u_p_f, psi_p_f, cb, wilson_p);


    /* wait for finishing communication */
    for(int idir=0;idir<4;++idir)
      if(GJP.Nodes(idir)!=1){
	//printf("wait for communication for %d\n",i);
	//fflush(stdin);

	QMP_status_t  status = QMP_wait(wilson_p->multiple[idir]);
	if (status != QMP_SUCCESS)
	  QMP_error("Error in GlobalDataShift::Shift:%s\n",
		    QMP_error_string(status));

        status = QMP_wait(wilson_p->multiple[idir+4]);
	if (status != QMP_SUCCESS)
	  QMP_error("Error in GlobalDataShift::Shift:%s\n",
		    QMP_error_string(status));
      }

    
#ifdef PROFILE

  time_bnd+= time+dclock();
  time = -dclock();
#endif

  /* Do the bulk part */

  if(dag == 0) {
    wilson_dslash_blk_dag0(chi_p_f, u_p_f, psi_p_f, cb, wilson_p);
  }
  else{
    wilson_dslash_blk_dag1(chi_p_f, u_p_f, psi_p_f, cb, wilson_p);
  }

#ifdef PROFILE
  time_blk+= time+dclock();
  }
#endif

#ifdef PROFILE
  Float time_tot= time_bnd+time_blk;
  const Float MultFlops_percent = MultFlops/time_tot
    *1e-9/(CPU_GHZ * 4* num_threads)*100;
  const Float MultFlops_blk_percent = MultFlops_blk/time_blk
    *1e-9/(CPU_GHZ * 4* num_threads)*100;
  const Float MultFlops_bnd_percent = MultFlops_bnd/time_bnd
    *1e-9/(CPU_GHZ * 4* num_threads)*100;
  

  PROF_printf("(cb %d dag %d) %e flops/ %e seconds = %e GFlops [ %4.2f %% peak] (%d threads)\n", 
       cb,dag,(double)MultFlops,time_tot, MultFlops/time_tot*1e-9,
       MultFlops_percent, num_threads);

  PROF_printf("  blk: %4.3e sec %4.3e GFlops (%2.1f %% tm, %2.1f %%cmp) [%4.2f %% peak]\n",
	 time_blk, MultFlops_blk/time_blk*1e-9, 
	 100*time_blk/time_tot, 100*(MultFlops_blk/MultFlops),
	 MultFlops_blk_percent);
  PROF_printf("  bnd: %4.3e sec %4.3e GFlops (%2.1f %% tm, %2.1f %%cmp) [%4.2f %% peak]\n",
	 time_bnd, MultFlops_bnd/time_bnd*1e-9, 
	 100*time_bnd/time_tot, 100*(MultFlops_bnd/MultFlops),
	 MultFlops_bnd_percent);
  
#endif
DiracOp::CGflops += 1320*vol;

}


}//extern"C"

#ifdef SSE_TO_C
#define SSE_TO_C2
#undef SSE_TO_C
#endif

#define BND_COMM


#ifndef SSE_TO_C2
#define SSE_C_FLOAT Float
#include "sse-blk-dag0.h"
#include "sse-blk-dag1.h"
#include "sse-bnd-dag0.h"
#include "sse-bnd-dag1.h"
#endif


CPS_END_NAMESPACE

#endif
