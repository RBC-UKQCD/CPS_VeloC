// matmul.cpp: Matrix Multiplication Example using OpenMP Offloading 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>

#include <config.h>
#include <precision.h>
#include <util/time_cps.h>

#include <nvToolsExt.h>
// Color definitions for nvtx calls
#define CLR_RED 0xFFFF0000
#define CLR_BLUE 0xFF0000FF
#define CLR_GREEN 0xFF008000
#define CLR_YELLOW 0xFFFFFF00
#define CLR_CYAN 0xFF00FFFF
#define CLR_MAGENTA 0xFFFF00FF
#define CLR_GRAY 0xFF808080
#define CLR_PURPLE 0xFF800080

// Macro for calling nvtxRangePushEx
#define RANGE_PUSH(range_title,range_color) { \
nvtxEventAttributes_t eventAttrib = {0}; \
eventAttrib.version = NVTX_VERSION; \
eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE; \
eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;\
eventAttrib.colorType = NVTX_COLOR_ARGB; \
eventAttrib.color = range_color; \
eventAttrib.message.ascii = range_title; \
nvtxRangePushEx(&eventAttrib); \
}
// Macro for calling nvtxRangePop
#define RANGE_POP {\
nvtxRangePop();\
}


#define THREAD_WILSON_MATRIX
#include <alg/wilson_matrix.h>

#ifdef THREAD_WILSON_MATRIX
#define OMP_WM(x) Pragma(omp x)
#define OMP_LP(x)
#define VEC_LEN 1
#define NTH 144
#else
#define OMP_WM(x)
#define OMP_LP(x) Pragma(omp x)
#define VEC_LEN 4
#define NTH 4
#endif


USING_NAMESPACE_CPS
//static  int MAX=128;
#define MAX0  128
//#define VOL 131072 does not finish
#define VOL 32768
//#define VOL 16384
//const int VOL=16*16*16*8;
typedef int BOOL;
//typedef int TYPE;


template < typename TYPE >
  BOOL check_result (TYPE * actual, TYPE * expected, unsigned n, unsigned m)
{
  for (unsigned i = 0; i < n * m; i++) {
    if (actual[i] != expected[i]) {
//            printf("Value mismatch at index = %d %d. Expected: %d"
//                   ", Actual: %d.\n", i, expected[i], actual[i]);
      std::cout << "Value mismatch at index = " << i / m << " " << i %
        m << "  Expected: " << expected[i] << " " << actual[i] << std::endl;
      return 0;
    }
  }
  return 1;
}


static const WilsonMatrix & apply_op (WilsonMatrix & out,
                                      const WilsonMatrix & in, int op_id)
{
  switch (op_id) {
  case 0:
  case 1:
  case 2:
  case 3:
    return out.glV (in, op_id);
  case 4:
    return out = in;
  default:
    ERR.General ("", "apply_op", "Invalid op_id = %d\n", op_id);
    return out;
  }
}

#if 1
typedef struct
{
  Rcomplex d[144];
} wmhack;

#if 0
inline void cmad (Rcomplex & x, const Rcomplex & y, const Rcomplex & z)
{
  x += y * z;
}

inline void cmeq (Rcomplex & x, const Rcomplex & y, const Rcomplex & z)
{
  x = y * z;
}
#endif

inline void eq_mult (Rcomplex * xoff, const Rcomplex * a, const Rcomplex * b)
{
//  const Rcomplex* a(amat.d);
//  const Rcomplex* b(bmat.d);
//  Rcomplex* xoff(xmat.d);
  Rcomplex const *point;
//#pragma omp parallel for
  for (int i1 = 0; i1 < 12; ++i1) {
    point = b;
    register const Rcomplex & aval (*a);
    cmeq (xoff[0], aval, point[0]);
    cmeq (xoff[1], aval, point[1]);
    cmeq (xoff[2], aval, point[2]);
    cmeq (xoff[3], aval, point[3]);
    cmeq (xoff[4], aval, point[4]);
    cmeq (xoff[5], aval, point[5]);
    cmeq (xoff[6], aval, point[6]);
    cmeq (xoff[7], aval, point[7]);
    cmeq (xoff[8], aval, point[8]);
    cmeq (xoff[9], aval, point[9]);
    cmeq (xoff[10], aval, point[10]);
    cmeq (xoff[11], aval, point[11]);
    a++;
    point += 12;
    for (int i3 = 1; i3 < 12; ++i3) {
      register const Rcomplex & aval (*a);
      cmad (xoff[0], aval, point[0]);
      cmad (xoff[1], aval, point[1]);
      cmad (xoff[2], aval, point[2]);
      cmad (xoff[3], aval, point[3]);
      cmad (xoff[4], aval, point[4]);
      cmad (xoff[5], aval, point[5]);
      cmad (xoff[6], aval, point[6]);
      cmad (xoff[7], aval, point[7]);
      cmad (xoff[8], aval, point[8]);
      cmad (xoff[9], aval, point[9]);
      cmad (xoff[10], aval, point[10]);
      cmad (xoff[11], aval, point[11]);
      a++;
      point += 12;
    }
    xoff += 12;
  }
//  return xmat;
}

inline void eq_mult (wmhack & xmat, const wmhack & amat, const wmhack & bmat)
{
  const Rcomplex *a (amat.d);
  const Rcomplex *b (bmat.d);
  Rcomplex *xoff (xmat.d);
  eq_mult (xoff, a, b);
}
#endif

int main ()
{
//  if (argc>1) sscanf(argv[1],"%d",&VOL);
  Float dtime = -dclock ();
  int A[MAX0][MAX0], B[MAX0][MAX0], C[MAX0][MAX0], C_SERIAL[MAX0][MAX0];
  for (int all = 0; all < 0; all++) {
    for (int i = 0; i < MAX0; i++)
      for (int j = 0; j < MAX0; j++) {
        A[i][j] = i + j - 1;
        B[i][j] = i - j + 1;
        C_SERIAL[i][j] = 0;
        C[i][j] = 0;
      }

    time_t t1 = time (NULL), t2;
    dtime = -dclock ();
    if (1) {
      for (int i = 0; i < MAX0; i++)
        for (int j = 0; j < MAX0; j++)
          for (int k = 0; k < MAX0; k++)
            C_SERIAL[i][j] += A[i][k] * B[k][j];
      dtime += dclock ();
      t2 = time (NULL);
      t1 = t2;
      printf ("serial: MAX0=%d %ld %e sec\n", MAX0, t2 - t1, dtime);
//  Compute();

//void __attribute__ ((noinline)) Compute()
      dtime = -dclock ();
      {
#pragma omp target teams distribute parallel for map(to: A, B) map(tofrom: C) \
                                                   thread_limit(128)
        {
          for (int i = 0; i < MAX0; i++)
            for (int j = 0; j < MAX0; j++)
              for (int k = 0; k < MAX0; k++)
                C[i][j] += A[i][k] * B[k][j];
        }
      }
      dtime += dclock ();
      t2 = time (NULL);
      printf ("OpenMP: VOL=%d %ld %e sec\n", MAX0, t2 - t1, dtime);
      t1 = t2;

      if (!check_result ((int *) &C[0][0], (int *) &C_SERIAL[0][0], MAX0, MAX0)) {
        printf ("FAILED\n");
//      return 1;
      }
      printf ("PASSED\n");
    }
  }


#if 1
  {
    wmhack mat[VOL], mat2[VOL], mat3[VOL], mat4[VOL], mat5[VOL];
    printf ("mat[0](%p)-mat[0](%p)=%d\n", mat[1].d, mat[0].d,
            (mat[1].d - mat[0].d));

//    { int all=0;
    for (int all = 0; all < 5; all++) {
      for (int i = 0; i < VOL; i++) {
//      printf("mat %d ?= %d/%d \n",i, (mat[i].d - mat[0].d),sizeof(wilson_matrix));
        for (int c1 = 0; c1 < 144; c1++) {
          mat[i].d[c1].real (drand48 ());
          mat[i].d[c1].imag (drand48 ());
          mat2[i].d[c1].real (drand48 ());
          mat2[i].d[c1].imag (drand48 ());
          mat5[i].d[c1].real (drand48 ());
          mat5[i].d[c1].imag (drand48 ());
          mat3[i].d[c1].real (0.);
          mat3[i].d[c1].imag (0.);
          mat4[i].d[c1].real (0.);
          mat4[i].d[c1].imag (0.);
        }
      }

      dtime = -dclock ();
#pragma omp parallel for
      for (int i = 0; i < VOL; i++) {
        eq_mult (mat3[i], mat[i], mat2[i]);
//        eq_mult (mat3[i], mat[i], mat5[i]);
      }
      dtime += dclock ();
      printf ("CPU(OpenMP): all=%d VOL=%d %e sec\n", all, VOL, dtime);
//  printf("CPU(Serial): VOL=%d %e sec\n",VOL,dtime);

      Rcomplex *mat_p = (Rcomplex *) mat[0].d;
      Rcomplex *mat2_p = (Rcomplex *) mat2[0].d;
      Rcomplex *mat3_p = (Rcomplex *) mat3[0].d;
      Rcomplex *mat4_p = (Rcomplex *) mat4[0].d;

      dtime = -dclock ();
//#pragma omp target teams distribute parallel for map(to: mat, mat2) map(from: mat4) thread_limit(128)
#pragma omp target teams distribute map(to: mat, mat2) map(from: mat4) parallel for collapse(3)
//thread_limit(144)
      for (int i = 0; i < VOL; i++) {
        for (int j = 0; j < 12; j++)
          for (int k = 0; k < 12; k++) {
            Rcomplex *mat_p = (Rcomplex *) mat[i].d;
            Rcomplex *mat2_p = (Rcomplex *) mat2[i].d;
            Rcomplex *mat5_p = (Rcomplex *) mat5[i].d;
            Rcomplex *mat4_p = (Rcomplex *) mat4[i].d;

//#pragma omp parallel for collapse(2)
//#pragma omp parallel for 
            mat4_p[j * 12 + k] = 0.;
            for (int l = 0; l < 12; l++)
              mat4_p[j * 12 + k] += mat_p[j * 12 + l] * mat2_p[l * 12 + k];
          }

//       eq_mult((Rcomplex*) mat4_p,mat_p,mat2_p);
      }
      dtime += dclock ();
      printf ("GPU(OpenMP): all=%d VOL=%d %e sec\n", all, VOL, dtime);

      Float max_diff = 0.;
      for (int i = 0; i < VOL; i++) {
        Rcomplex *mat3_p = mat3[i].d;
        Rcomplex *mat4_p = mat4[i].d;
        for (int j = 0; j < 144; j++) {
          if (norm (*(mat3_p + j) - *(mat4_p + j)) > max_diff) {
            max_diff = norm (*(mat3_p + j) - *(mat4_p + j));
            std::cout << i << " " << j << " : mat3= " << *(mat3_p +
                                                           j) << " mat4= " <<
              *(mat4_p + j) << " diff= " << max_diff << std::endl;
          }
        }
      }
    }
  }
#endif

  {
    const int total_ops = 5;
    const int t_size = 192;
    WilsonMatrix ltwst[VOL], sprop[VOL], lwsnk[VOL], p0[VOL], p1[VOL], p2[VOL],
      res1[VOL], res2[VOL];
    Rcomplex tr1[total_ops][t_size][VOL], tr2[total_ops][t_size][VOL];

    printf ("mat[0](%p)-mat[0](%p)=%d\n", ltwst[1].ptr (), ltwst[0].ptr (),
            (ltwst[1].ptr () - ltwst[0].ptr ()));

    for (int all = 0; all < 5; all++) {
    RANGE_PUSH ("Initialize Arrays (CPU)", CLR_BLUE);
      for (int i = 0; i < VOL; i++) {
        for (int c1 = 0; c1 < 3; c1++)
          for (int c2 = 0; c2 < 3; c2++)
            for (int s1 = 0; s1 < 4; s1++)
              for (int s2 = 0; s2 < 4; s2++) {
                ltwst[i] (s1, c1, s2, c2) = Rcomplex (drand48 (), drand48 ());
                sprop[i] (s1, c1, s2, c2) = Rcomplex (drand48 (), drand48 ());
                lwsnk[i] (s1, c1, s2, c2) = Rcomplex (drand48 (), drand48 ());
                p0[i] (s1, c1, s2, c2) = Rcomplex (0., 0.);
                p1[i] (s1, c1, s2, c2) = Rcomplex (0., 0.);
                p2[i] (s1, c1, s2, c2) = Rcomplex (0., 0.);
              }
      }
      RANGE_POP;

      dtime = -dclock ();
      RANGE_PUSH ("OpenMP (CPU)", CLR_GREEN);
#pragma omp parallel for
      for (int i = 0; i < VOL; i++) {
        p0[i].glV (sprop[i], -5);
        p2[i].glV (lwsnk[i], -5);
        p0[i].hconj ();
        p2[i].hconj ();
        res1[i] = ltwst[i];
        res1[i] *= p2[i] * p0[i];
#if 0
        for (int op_id = 0; op_id < total_ops; ++op_id) {
          apply_op (p0[i], res1[i], op_id);
          tr1[op_id][t_glb][i] = p0[i].Trace ();
        }
#endif

      }
      RANGE_POP;
      dtime += dclock ();
      printf ("CPU(OpenMP): VOL=%d %e sec\n", VOL, dtime);
      dtime = -dclock ();
//      const int VEC_LEN = 16;
//#pragma omp target teams distribute map(to: ltwst, sprop,lwsnk) map(from: p2) thread_limit(144)
#pragma omp target teams distribute map(to: ltwst, sprop,lwsnk) map(from: res2)
//thread_limit(NTH)
      for (int index = 0; index < VOL; index += VEC_LEN) {
        OMP_LP (parallel for)
          for (int ii = 0; ii < VEC_LEN; ii++) {
            int i = index + ii;
            Rcomplex *ltwst_p = ltwst[i].ptr ();
            Rcomplex *sprop_p = sprop[i].ptr ();
            Rcomplex *lwsnk_p = lwsnk[i].ptr ();
            Rcomplex *p0_p = p0[i].ptr ();
            Rcomplex *p1_p = p1[i].ptr ();
            Rcomplex *p2_p = p2[i].ptr ();
            Rcomplex *res2_p = res2[i].ptr ();

#if 0
            p0[i].glV (sprop[i], -5);
            p2[i].glV (lwsnk[i], -5);
#else
            p1[i].glV (sprop[i], -5);
            p0[i].hconj (p1[i]);
            p1[i].glV (lwsnk[i], -5);
            p2[i].hconj (p1[i]);
#endif

#if 0
            eq_mult (p1_p, ltwst_p, p2_p);
            eq_mult (res2_p, p1_p, p0_p);
#else
//#pragma omp parallel for collapse(2)
#if 0
            OMP_WM (parallel for collapse (2))
              for (int j = 0; j < 12; j++)
                for (int k = 0; k < 12; k++)
                  p1_p[j * 12 + k] = Rcomplex (0., 0.);
#endif
//#pragma unroll
            OMP_WM (parallel for collapse (2))
              for (int j = 0; j < 12; j++)
                for (int k = 0; k < 12; k++) {
                  p1_p[j * 12 + k] = Rcomplex (0., 0.);
                  for (int l = 0; l < 12; l++)
                    p1_p[j * 12 + k] += ltwst_p[j * 12 + l] * p2_p[l * 12 + k];
//                p1_p[j * 12 + k] += ltwst_p[j * 12 + l] * sprop_p[l * 12 + k];
                }
            OMP_WM (parallel for collapse (2))
              for (int j = 0; j < 12; j++)
                for (int k = 0; k < 12; k++) {
                  res2_p[j * 12 + k] = Rcomplex (0., 0.);
//#pragma unroll
                  for (int l = 0; l < 12; l++)
                    res2_p[j * 12 + k] += p1_p[j * 12 + l] * p0_p[l * 12 + k];
//                res2_p[j * 12 + k] += p1_p[j * 12 + l] * lwsnk_p[l * 12 + k];
                }
#endif
          }
#if 0
        for (int op_id = 0; op_id < total_ops; ++op_id) {
          apply_op (p0[i], res2[i], op_id);
          tr1[op_id][t_glb][i] = p0[i].Trace ();
        }
#endif
      }
      dtime += dclock ();
      printf ("GPU(OpenMP): VOL=%d %e sec\n", VOL, dtime);
//  Complex *mat1_p = mat1.p; 

      RANGE_PUSH ("Check (CPU)", CLR_RED);
      Float max_diff = 0.;
      for (int i = 0; i < VOL; i++) {
        Rcomplex *mat3_p = res1[i].ptr ();
        Rcomplex *mat4_p = res2[i].ptr ();
        for (int j = 0; j < 144; j++) {
          if (norm (*(mat3_p + j) - *(mat4_p + j)) > max_diff) {
            max_diff = norm (*(mat3_p + j) - *(mat4_p + j));
            std::cout << i << " " << j << " : mat3= " << *(mat3_p +
                                                           j) << " mat4= " <<
              *(mat4_p + j) << " diff= " << max_diff << std::endl;
          }
        }
      }
      RANGE_POP;
    }
  }


  return 0;
}
