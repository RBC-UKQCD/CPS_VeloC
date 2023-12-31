#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//--------------------------------------------------------------------
//------------------------------------------------------------------
// This routine computes the inverse of the 5th dimension hopping term
// for MOBIUS, plus the diagonal term for the whole 5 dim MOBIUS Dirac operator.
// It is the inverse of 1 - kappa * mobius_dslash_5_plus
// where mobius_dslash_5_plus is the derivative part of the 5th direction
// part of the fermion matrix.
//
// It is part of the 4d odd-even preconditioned MOBIUS operator
//
// Storage order for MOBIUS fermions
//------------------------------------------------------------------
//  
//  |     |
//  | |r| |
//  | |i| |
//  |     |
//  | |r| | = |spin comp|
//  | |i| |
//  |     |
//  | |r| |
//  | |i| |
//  |     |
//  
//  
//  |             |
//  | |spin comp| |
//  |             |
//  |             |
//  | |spin comp| | 
//  |             | = |spinor|
//  |             |
//  | |spin comp| |
//  |             |
//  |             |
//  | |spin comp| |
//  |             |
//  
//
//  |            |
//  |  |spinor|  |
//  |  |spinor|  |
//  |  |spinor|  |
//  |  |spinor|  |
//  |  |spinor|  |
//  |     .      | = |s-block|   The spinors are arranged in Wilson
//  |     .      |               order with odd - even 4d-checkerboard
//  |     .      |               storage.



//  for 4d preconditioned operator,
//  |                |
//  | |s-block odd|  |  For odd chckerboard (same for even)
//  | |s-block odd|  |
//  | |s-block odd|  |
//  |       .        |
//  |       .        |
//  |       .        |
//  |                |

//
//------------------------------------------------------------------


CPS_END_NAMESPACE
#include<util/zmobius.h>
#include<util/gjp.h>
#include<util/dirac_op.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/smalloc.h>
#include<comms/scu.h>

#include "blas-subs.h"

CPS_START_NAMESPACE

/*
  To Chulwoo,

    There are several equivalent version for m5^{-1}.

    The before_loop_unrolling version, which has "two s-loop",
    is corresponding to the one in the note, and would be easy to read and modify.
    
    There are also unroll loop version which is trivially equivalent to the above, but lengthy.
    This is currently used. But you may want to switch to use the versions before unrolling
    as reference implementations. 

    At the end of this file, there are version with four s-loop, which corresponds to earlier version in the note.


    By changing the Float fact into array, we could easily implement the part of Mobius.
    
 */





#if 0
//----------------------------------------
//
//  dag 0 two s-loop version, after loop unrolling
//
void zmobius_m5inv_dag0(Vector *inout, 
		       const Float mass,
		       const Zmobus *mobius_lib_arg)
{
  
  int x;
  int s;

  // Initializations
  //------------------------------------------------------------------
  const int ls = GJP.SnodeSites()*GJP.Snodes();
  // const int s_nodes = GJP.Snodes();
  // const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = mobius_lib_arg->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;
  const size_t f_size = 24 * vol_4d_cb * ls;

  // two_kappa is (- 2 \kappa) in the note
  const IFloat two_kappa = -mobius_lib_arg->mobius_kappa_b/mobius_lib_arg->mobius_kappa_c;
  const IFloat inv_two_kappa = 1.0 / two_kappa;
  const Float  inv_d_last = 1.0 / ( 1.0 + pow(two_kappa, ls)*mass); // 1.0 / d_{ls-1}

  Float fact;

  Float* const f_out = (IFloat *) inout;

  //time_elapse();

  //-----------------------------------------------------
  //  The first forward ls loop both for upper  and downer spinors
  //----------------------------------------------------

  // s = 0
  //-------
  //#pragma ivdep
  for(x=0; x<vol_4d_cb; x+=4) {
    // downer part  fout[x,ls-1] *= d_last;
    VEC_TIMESEQU_FLOAT(f_out+ 12+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(24+f_out+ 12+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(48+f_out+ 12+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(72+f_out+ 12+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
  }
  
  // FIXME  :   perhaps using array d[i] is faster ?
  fact = - two_kappa * mass * inv_d_last ; // - d0 / d_{ls-1}

  // s = 1 ... ls - 2
  //-------------------
  for(s=0; s<= ls-2 ;++s) {
    //#pragma ivdep
    for(x=0; x<vol_4d_cb; x+=4) {
      // upper part  fout[x,s+1] +=  two_kappa fout[x,s]
      AXPY(12, two_kappa, f_out + 24*(x+vol_4d_cb*s),f_out + 24*(x+vol_4d_cb*(s+1)));
      // downer part fout[x,ls-1] +=  fact * fout[x,s]
      AXPY(12, fact, f_out + 12+24*(x+vol_4d_cb*s),f_out + 12+24*(x+vol_4d_cb*(ls-1)));

      // upper part  fout[x,s+1] +=  two_kappa fout[x,s]
      AXPY(12, two_kappa, 24+f_out + 24*(x+vol_4d_cb*s),24+f_out + 24*(x+vol_4d_cb*(s+1)));
      // downer part fout[x,ls-1] +=  fact * fout[x,s]
      AXPY(12, fact, 24+f_out + 12+24*(x+vol_4d_cb*s),24+f_out + 12+24*(x+vol_4d_cb*(ls-1)));

      // upper part  fout[x,s+1] +=  two_kappa fout[x,s]
      AXPY(12, two_kappa, 48+f_out + 24*(x+vol_4d_cb*s),48+f_out + 24*(x+vol_4d_cb*(s+1)));
      // downer part fout[x,ls-1] +=  fact * fout[x,s]
      AXPY(12, fact, 48+f_out + 12+24*(x+vol_4d_cb*s),48+f_out + 12+24*(x+vol_4d_cb*(ls-1)));

      // upper part  fout[x,s+1] +=  two_kappa fout[x,s]
      AXPY(12, two_kappa, 72+f_out + 24*(x+vol_4d_cb*s),72+f_out + 24*(x+vol_4d_cb*(s+1)));
      // downer part fout[x,ls-1] +=  fact * fout[x,s]
      AXPY(12, fact, 72+f_out + 12+24*(x+vol_4d_cb*s),72+f_out + 12+24*(x+vol_4d_cb*(ls-1)));
    }
    fact *= two_kappa;
  }

  //DEBUG_MOBIUS_DSLASH("dag0 first part %e\n", time_elapse());
  
  //-------------------------------------------------------------------
  //  The second backward ls loop both for upper  and downer spinors
  //-------------------------------------------------------------------
  // FIXME: perhaps x loop should also be the reverse order ?
  
  
  fact = - pow(two_kappa,ls-1) * mass * inv_d_last ; // - d_{ls-2} / d_{ls-1}

  // s = ls-2, ... ,  0
  //----------------------
  for(s=ls-2; s >=0 ; --s){
    //#pragma ivdep
    for(x=0; x<vol_4d_cb; x+=4) {
      // upper fout[x,s] += fact* fout[x,ls-1]
      AXPY(12, fact, f_out +24*(x+vol_4d_cb*(ls-1)), f_out + 24*(x+vol_4d_cb*s));
      // downer fout[x,s] += two_kappa fout[x,s+1]
      AXPY(12, two_kappa, f_out+12+24*(x+vol_4d_cb*(s+1)),f_out + 12+24*(x+vol_4d_cb*s));
      
      // upper fout[x,s] += fact* fout[x,ls-1]
      AXPY(12, fact, 24+f_out +24*(x+vol_4d_cb*(ls-1)), 24+f_out + 24*(x+vol_4d_cb*s));
      // downer fout[x,s] += two_kappa fout[x,s+1]
      AXPY(12, two_kappa, 24+f_out+12+24*(x+vol_4d_cb*(s+1)),24+f_out + 12+24*(x+vol_4d_cb*s));
      
      // upper fout[x,s] += fact* fout[x,ls-1]
      AXPY(12, fact, 48+f_out +24*(x+vol_4d_cb*(ls-1)), 48+f_out + 24*(x+vol_4d_cb*s));
      // downer fout[x,s] += two_kappa fout[x,s+1]
      AXPY(12, two_kappa, 48+f_out+12+24*(x+vol_4d_cb*(s+1)),48+f_out + 12+24*(x+vol_4d_cb*s));

      // upper fout[x,s] += fact* fout[x,ls-1]
      AXPY(12, fact, 72+f_out +24*(x+vol_4d_cb*(ls-1)), 72+f_out + 24*(x+vol_4d_cb*s));
      // downer fout[x,s] += two_kappa fout[x,s+1]
      AXPY(12, two_kappa, 72+f_out+12+24*(x+vol_4d_cb*(s+1)),72+f_out + 12+24*(x+vol_4d_cb*s));

    }
    //FIXME : which is faster ?
    //fact *= inv_two_kappa;
    fact /= two_kappa;
  }

  // s = ls - 1
  //-------------

  for(x=0; x<vol_4d_cb; x+=4) {
    // upper   fou[x,ls-1] *= inv_d_last 
    VEC_TIMESEQU_FLOAT(f_out+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(24+f_out+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);      
    VEC_TIMESEQU_FLOAT(48+f_out+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);  
    VEC_TIMESEQU_FLOAT(72+f_out+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);  
  }

  //DEBUG_MOBIUS_DSLASH("dag0 second part %e\n", time_elapse());
  DiracOp::CGflops+=vol_4d_cb*(ls*96-48);  
  //DiracOp::CGflops+=2*2*vol_4d_cb*local_ls*12;
}
#endif
//----------------------------------------
//
//  dag 0 two s-loop version, after loop unrolling
//
void zmobius_m5inv_cmplx_dag0(Vector *inout, 
			      const Float mass,
			      const Zmobus *mobius_lib_arg,
			      Complex* K)
{
  
//  int x;
//  int s;

  // Initializations
  //------------------------------------------------------------------
  const int ls = GJP.SnodeSites()*GJP.Snodes();
  // const int s_nodes = GJP.Snodes();
  // const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = mobius_lib_arg->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;
  const size_t f_size = 24 * vol_4d_cb * ls;

  // two_kappa is (- 2 \kappa) in the note
  //const IFloat two_kappa = -mobius_lib_arg->mobius_kappa_b/mobius_lib_arg->mobius_kappa_c;
  //const IFloat inv_two_kappa = 1.0 / two_kappa;
  //const Float  inv_d_last = 1.0 / ( 1.0 + pow(two_kappa, ls)*mass); // 1.0 / d_{ls-1}

  //const Complex* K = mobius_lib_arg->zmobius_kappa_ratio;
  Complex inv_d_last;
  {
    Complex prod=1.0;
    for(int i=0;i<ls;++i) prod *= -K[i];
    prod *= mass;
    inv_d_last = 1.0/ (1.0 + prod);
  }


  
  Complex factL, factR;

  Float* const f_out = (IFloat *) inout;

  //time_elapse();

  //-----------------------------------------------------
  //  The first forward ls loop both for upper  and downer spinors
  //----------------------------------------------------

  // s = 0
  //-------
  //#pragma ivdep
#pragma omp parallel for
  for(int x=0; x<vol_4d_cb; x+=4) {
    // downer part  fout[x,ls-1] *= inv_d_last;
    VEC_TIMESEQU_COMPLEX(f_out+ 12+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_COMPLEX(24+f_out+ 12+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_COMPLEX(48+f_out+ 12+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_COMPLEX(72+f_out+ 12+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
  }
  
  // FIXME  :   perhaps using array d[i] is faster ?
  //  factL = -d_s / d_Ls-1 , note the sign
  factL =  K[ls-1] * mass * inv_d_last ; // - d0 / d_{ls-1}

  // s = 1 ... ls - 2
  //-------------------
  for(int s=0; s<= ls-2 ;++s) {
    //#pragma ivdep
#pragma omp parallel for
    for(int x=0; x<vol_4d_cb; x+=4) {
      // upper part  fout[x,s+1] +=  -K[s+1] fout[x,s]
      ZAXPY(12, -K[s+1], f_out + 24*(x+vol_4d_cb*s),f_out + 24*(x+vol_4d_cb*(s+1)));
      // downer part fout[x,ls-1] +=  factL * fout[x,s]
      ZAXPY(12, factL, f_out + 12+24*(x+vol_4d_cb*s),f_out + 12+24*(x+vol_4d_cb*(ls-1)));

      // upper part  fout[x,s+1] +=  -K[s+1] fout[x,s]
      ZAXPY(12, -K[s+1], 24+f_out + 24*(x+vol_4d_cb*s),24+f_out + 24*(x+vol_4d_cb*(s+1)));
      // downer part fout[x,ls-1] +=  factL * fout[x,s]
      ZAXPY(12, factL, 24+f_out + 12+24*(x+vol_4d_cb*s),24+f_out + 12+24*(x+vol_4d_cb*(ls-1)));

      // upper part  fout[x,s+1] +=  -K[s+1] fout[x,s]
      ZAXPY(12, -K[s+1], 48+f_out + 24*(x+vol_4d_cb*s),48+f_out + 24*(x+vol_4d_cb*(s+1)));
      // downer part fout[x,ls-1] +=  factL * fout[x,s]
      ZAXPY(12, factL, 48+f_out + 12+24*(x+vol_4d_cb*s),48+f_out + 12+24*(x+vol_4d_cb*(ls-1)));

      // upper part  fout[x,s+1] +=  -K[s+1] fout[x,s]
      ZAXPY(12, -K[s+1], 72+f_out + 24*(x+vol_4d_cb*s),72+f_out + 24*(x+vol_4d_cb*(s+1)));
      // downer part fout[x,ls-1] +=  factL * fout[x,s]
      ZAXPY(12, factL, 72+f_out + 12+24*(x+vol_4d_cb*s),72+f_out + 12+24*(x+vol_4d_cb*(ls-1)));
    }
    factL *= -K[s];
  }

  //DEBUG_MOBIUS_DSLASH("dag0 first part %e\n", time_elapse());
  
  //-------------------------------------------------------------------
  //  The second backward ls loop both for upper  and downer spinors
  //-------------------------------------------------------------------
  // FIXME: perhaps x loop should also be the reverse order ?
  
  {// factR = - d_Ls-2 / d_Ls-1, note the sign 
    Complex prod=1.0;
    for(int i=0;i<ls-1;++i) prod *= -K[i];
    factR = - prod*mass*inv_d_last;
  }

  // s = ls-2, ... ,  0
  //----------------------
  for(int s=ls-2; s >=0 ; --s){
    //#pragma ivdep
#pragma omp parallel for
    for(int x=0; x<vol_4d_cb; x+=4) {
      // upper fout[x,s] += factR* fout[x,ls-1]
      ZAXPY(12, factR, f_out +24*(x+vol_4d_cb*(ls-1)), f_out + 24*(x+vol_4d_cb*s));
      // downer fout[x,s] += -K[s] fout[x,s+1]
      ZAXPY(12, -K[s], f_out+12+24*(x+vol_4d_cb*(s+1)),f_out + 12+24*(x+vol_4d_cb*s));
      
      // upper fout[x,s] += factR* fout[x,ls-1]
      ZAXPY(12, factR, 24+f_out +24*(x+vol_4d_cb*(ls-1)), 24+f_out + 24*(x+vol_4d_cb*s));
      // downer fout[x,s] += -K[s] fout[x,s+1]
      ZAXPY(12, -K[s], 24+f_out+12+24*(x+vol_4d_cb*(s+1)),24+f_out + 12+24*(x+vol_4d_cb*s));
      
      // upper fout[x,s] += factR* fout[x,ls-1]
      ZAXPY(12, factR, 48+f_out +24*(x+vol_4d_cb*(ls-1)), 48+f_out + 24*(x+vol_4d_cb*s));
      // downer fout[x,s] += -K[s] fout[x,s+1]
      ZAXPY(12, -K[s], 48+f_out+12+24*(x+vol_4d_cb*(s+1)),48+f_out + 12+24*(x+vol_4d_cb*s));

      // upper fout[x,s] += factR* fout[x,ls-1]
      ZAXPY(12, factR, 72+f_out +24*(x+vol_4d_cb*(ls-1)), 72+f_out + 24*(x+vol_4d_cb*s));
      // downer fout[x,s] += -K[s] fout[x,s+1]
      ZAXPY(12, -K[s], 72+f_out+12+24*(x+vol_4d_cb*(s+1)),72+f_out + 12+24*(x+vol_4d_cb*s));

    }
    //FIXME : which is faster ?
    //fact *= inv_two_kappa;
    factR /= -K[s];
  }

  // s = ls - 1
  //-------------

#pragma omp parallel for
  for(int x=0; x<vol_4d_cb; x+=4) {
    // upper   fou[x,ls-1] *= inv_d_last 
    VEC_TIMESEQU_COMPLEX(f_out+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_COMPLEX(24+f_out+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);      
    VEC_TIMESEQU_COMPLEX(48+f_out+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);  
    VEC_TIMESEQU_COMPLEX(72+f_out+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);  
  }

  //DEBUG_MOBIUS_DSLASH("dag0 second part %e\n", time_elapse());
  DiracOp::CGflops+=vol_4d_cb*(ls*96-48);  
  //DiracOp::CGflops+=2*2*vol_4d_cb*local_ls*12;
}

#if 0
//----------------------------------------
//
//  dag 1 two s-loop version, after loop unrolling
//
void zmobius_m5inv_dag1(Vector *inout, 
		    const Float mass,
			const Zmobus *mobius_lib_arg)
{

  int x;
  int s;

// Initializations
//------------------------------------------------------------------
  const int ls = GJP.SnodeSites()*GJP.Snodes();
  // const int s_nodes = GJP.Snodes();
  // const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = mobius_lib_arg->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;
  const size_t f_size = 24 * vol_4d_cb * ls;

  const IFloat two_kappa = -mobius_lib_arg->mobius_kappa_b/mobius_lib_arg->mobius_kappa_c;
  const IFloat inv_two_kappa = 1.0 / two_kappa;
  const Float  inv_d_last = 1.0 / ( 1.0 + pow(two_kappa, ls)*mass); // 1.0 / d_{ls-1}

  Float fact;
  
  Float* const f_out = (IFloat *) inout;

  //time_elapse();
  //-----------------------------------------------------
  //  The first forward ls loop both for upper  and downer spinors
  //----------------------------------------------------

  // s = 0
  //-------
  for(x=0; x<vol_4d_cb; x+=4) {
    // upper part  fout[x,ls-1] *= d_last;
    VEC_TIMESEQU_FLOAT(f_out+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(24+f_out+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(48+f_out+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(72+f_out+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
  }
  
  // FIXME  :   perhaps using array d[i] is faster ?
  fact = - two_kappa * mass * inv_d_last ; // - d0 / d_{ls-1}


  // s = 1 ... ls - 2
  //-------------------
  for(s=0; s<= ls-2 ;++s) {
    for(x=0; x<vol_4d_cb; x+=4) {
      // upper part fout[x,ls-1] +=  fact * fout[x,s]
      AXPY(12, fact, f_out +24*(x+vol_4d_cb*s),f_out +24*(x+vol_4d_cb*(ls-1)));
      // downer part  fout[x,s+1] +=  two_kappa fout[x,s]
      AXPY(12, two_kappa, f_out+ 12 + 24*(x+vol_4d_cb*s),f_out+ 12 + 24*(x+vol_4d_cb*(s+1)));

      // upper part fout[x,ls-1] +=  fact * fout[x,s]
      AXPY(12, fact, 24+f_out +24*(x+vol_4d_cb*s),24+f_out +24*(x+vol_4d_cb*(ls-1)));
      // downer part  fout[x,s+1] +=  two_kappa fout[x,s]
      AXPY(12, two_kappa, 24+f_out+ 12 + 24*(x+vol_4d_cb*s),24+f_out+ 12 + 24*(x+vol_4d_cb*(s+1)));

      // upper part fout[x,ls-1] +=  fact * fout[x,s]
      AXPY(12, fact, 48+f_out +24*(x+vol_4d_cb*s),48+f_out +24*(x+vol_4d_cb*(ls-1)));
      // downer part  fout[x,s+1] +=  two_kappa fout[x,s]
      AXPY(12, two_kappa, 48+f_out+ 12 + 24*(x+vol_4d_cb*s),48+f_out+ 12 + 24*(x+vol_4d_cb*(s+1)));

      // upper part fout[x,ls-1] +=  fact * fout[x,s]
      AXPY(12, fact, 72+f_out +24*(x+vol_4d_cb*s),72+f_out +24*(x+vol_4d_cb*(ls-1)));
      // downer part  fout[x,s+1] +=  two_kappa fout[x,s]
      AXPY(12, two_kappa, 72+f_out+ 12 + 24*(x+vol_4d_cb*s),72+f_out+ 12 + 24*(x+vol_4d_cb*(s+1)));
    }
    fact *= two_kappa;
  }

  //DEBUG_MOBIUS_DSLASH("dag1 first part %e\n", time_elapse());
  
  //-------------------------------------------------------------------
  //  The second backward ls loop both for upper  and downer spinors
  //-------------------------------------------------------------------
  // FIXME: perhaps x loop should also be the reverse order ?
  
  
  fact = - pow(two_kappa,ls-1) * mass * inv_d_last ; // - d_{ls-2} / d_{ls-1}

  // s = ls-2, ... ,  0
  //----------------------
  for(s=ls-2; s >=0 ; --s){
    for(x=0; x<vol_4d_cb; x+=4) {
      // upper fout[x,s] += two_kappa fout[x,s+1]
      AXPY(12, two_kappa, f_out+24*(x+vol_4d_cb*(s+1)),f_out +24*(x+vol_4d_cb*s));
      // downer fout[x,s] += fact* fout[x,ls-1]
      AXPY(12, fact, f_out +12+24*(x+vol_4d_cb*(ls-1)), f_out + 12 + 24*(x+vol_4d_cb*s));

      // upper fout[x,s] += two_kappa fout[x,s+1]
      AXPY(12, two_kappa, 24+f_out+24*(x+vol_4d_cb*(s+1)),24+f_out +24*(x+vol_4d_cb*s));
      // downer fout[x,s] += fact* fout[x,ls-1]
      AXPY(12, fact, 24+f_out +12+24*(x+vol_4d_cb*(ls-1)), 24+f_out + 12 + 24*(x+vol_4d_cb*s));

      // upper fout[x,s] += two_kappa fout[x,s+1]
      AXPY(12, two_kappa, 48+f_out+24*(x+vol_4d_cb*(s+1)),48+f_out +24*(x+vol_4d_cb*s));
      // downer fout[x,s] += fact* fout[x,ls-1]
      AXPY(12, fact, 48+f_out +12+24*(x+vol_4d_cb*(ls-1)), 48+f_out + 12 + 24*(x+vol_4d_cb*s));

      // upper fout[x,s] += two_kappa fout[x,s+1]
      AXPY(12, two_kappa, 72+f_out+24*(x+vol_4d_cb*(s+1)),72+f_out +24*(x+vol_4d_cb*s));
      // downer fout[x,s] += fact* fout[x,ls-1]
      AXPY(12, fact, 72+f_out +12+24*(x+vol_4d_cb*(ls-1)), 72+f_out + 12 + 24*(x+vol_4d_cb*s));
}
    //FIXME : which is faster ?
    //fact *= inv_two_kappa;
    fact /= two_kappa;
  }

  // s = ls - 1
  //-------------
  for(x=0; x<vol_4d_cb; x+=4) {
    // downer  fou[x,ls-1] *= inv_d_last 
    VEC_TIMESEQU_FLOAT(f_out+12+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(24+f_out+12+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(48+f_out+12+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(72+f_out+12+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
  }

  //DEBUG_MOBIUS_DSLASH("dag1 second part %e\n", time_elapse());
  
  //DiracOp::CGflops+=2*2*vol_4d_cb*local_ls*12;
  DiracOp::CGflops+=vol_4d_cb*(ls*96-48);
}
#endif

//----------------------------------------
//
//  dag 1 two s-loop version, after loop unrolling
//
void zmobius_m5inv_cmplx_dag1(Vector *inout, 
			      const Float mass,
			      const Zmobus *mobius_lib_arg,
			      Complex* K)
{

//  int x;
//  int s;

// Initializations
//------------------------------------------------------------------
  const int ls = GJP.SnodeSites()*GJP.Snodes();
  // const int s_nodes = GJP.Snodes();
  // const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = mobius_lib_arg->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;
  const size_t f_size = 24 * vol_4d_cb * ls;

  //const IFloat two_kappa = -mobius_lib_arg->mobius_kappa_b/mobius_lib_arg->mobius_kappa_c;
  //const IFloat inv_two_kappa = 1.0 / two_kappa;
  //const Float  inv_d_last = 1.0 / ( 1.0 + pow(two_kappa, ls)*mass); // 1.0 / d_{ls-1}

  //const Complex* K = mobius_lib_arg->zmobius_kappa_ratio;
  Complex factR, factL;
  
  Complex inv_d_last;
  {
    Complex prod=1.0;
    for(int i=0;i<ls;++i) prod *= K[i];
    prod *= mass;
    inv_d_last = 1.0/(1.0+prod);
    inv_d_last = conj(inv_d_last);
  }

  Float* const f_out = (IFloat *) inout;

  //time_elapse();
  //-----------------------------------------------------
  //  The first forward ls loop both for upper  and downer spinors
  //----------------------------------------------------

  // s = 0
  //-------
#pragma omp parallel for
  for(int x=0; x<vol_4d_cb; x+=4) {
    // upper part  fout[x,ls-1] *= d_last;
    VEC_TIMESEQU_COMPLEX(f_out+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_COMPLEX(24+f_out+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_COMPLEX(48+f_out+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_COMPLEX(72+f_out+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
  }
  
  // FIXME  :   perhaps using array d[i] is faster ?
  //fact = - two_kappa * mass * inv_d_last ; // - d0 / d_{ls-1}
  factR = conj(K[0])*mass * inv_d_last;

  // s = 1 ... ls - 2
  //-------------------
  for(int s=0; s<= ls-2 ;++s) {
#pragma omp parallel for
    for(int x=0; x<vol_4d_cb; x+=4) {
      // upper part fout[x,ls-1] +=  fact * fout[x,s]
      ZAXPY(12, factR, f_out +24*(x+vol_4d_cb*s),f_out +24*(x+vol_4d_cb*(ls-1)));
      // downer part  fout[x,s+1] +=  -conj(K[s])  fout[x,s]
      ZAXPY(12, -conj(K[s]), f_out+ 12 + 24*(x+vol_4d_cb*s),f_out+ 12 + 24*(x+vol_4d_cb*(s+1)));

      // upper part fout[x,ls-1] +=  fact * fout[x,s]
      ZAXPY(12, factR, 24+f_out +24*(x+vol_4d_cb*s),24+f_out +24*(x+vol_4d_cb*(ls-1)));
      // downer part  fout[x,s+1] +=  -conj(K[s]) fout[x,s]
      ZAXPY(12, -conj(K[s]), 24+f_out+ 12 + 24*(x+vol_4d_cb*s),24+f_out+ 12 + 24*(x+vol_4d_cb*(s+1)));

      // upper part fout[x,ls-1] +=  fact * fout[x,s]
      ZAXPY(12, factR, 48+f_out +24*(x+vol_4d_cb*s),48+f_out +24*(x+vol_4d_cb*(ls-1)));
      // downer part  fout[x,s+1] +=  -conj(K[s]) fout[x,s]
      ZAXPY(12, -conj(K[s]), 48+f_out+ 12 + 24*(x+vol_4d_cb*s),48+f_out+ 12 + 24*(x+vol_4d_cb*(s+1)));

      // upper part fout[x,ls-1] +=  fact * fout[x,s]
      ZAXPY(12, factR, 72+f_out +24*(x+vol_4d_cb*s),72+f_out +24*(x+vol_4d_cb*(ls-1)));
      // downer part  fout[x,s+1] +=  -conj(K[s]) fout[x,s]
      ZAXPY(12, -conj(K[s]), 72+f_out+ 12 + 24*(x+vol_4d_cb*s),72+f_out+ 12 + 24*(x+vol_4d_cb*(s+1)));
    }
    factR *= -conj(K[s+1]);
  }

  //DEBUG_MOBIUS_DSLASH("dag1 first part %e\n", time_elapse());
  
  //-------------------------------------------------------------------
  //  The second backward ls loop both for upper  and downer spinors
  //-------------------------------------------------------------------
  // FIXME: perhaps x loop should also be the reverse order ?
  
  
  {
    Complex prod=1.0;
    for(int i=0;i<ls-2;++i) prod *= -K[i];
    prod *= -K[ls-1];
    factL = - conj(prod)*mass*inv_d_last;
  }
    

  // s = ls-2, ... ,  0
  //----------------------
  for(int s=ls-2; s >=0 ; --s){
#pragma omp parallel for
    for(int x=0; x<vol_4d_cb; x+=4) {
      // upper fout[x,s] += -conj(K[s+1]) fout[x,s+1]
      ZAXPY(12, -conj(K[s+1]), f_out+24*(x+vol_4d_cb*(s+1)),f_out +24*(x+vol_4d_cb*s));
      // downer fout[x,s] += factL* fout[x,ls-1]
      ZAXPY(12, factL, f_out +12+24*(x+vol_4d_cb*(ls-1)), f_out + 12 + 24*(x+vol_4d_cb*s));

      // upper fout[x,s] += -conj(K[s+1]) fout[x,s+1]
      ZAXPY(12, -conj(K[s+1]), 24+f_out+24*(x+vol_4d_cb*(s+1)),24+f_out +24*(x+vol_4d_cb*s));
      // downer fout[x,s] += factL* fout[x,ls-1]
      ZAXPY(12, factL, 24+f_out +12+24*(x+vol_4d_cb*(ls-1)), 24+f_out + 12 + 24*(x+vol_4d_cb*s));

      // upper fout[x,s] += -conj(K[s+1]) fout[x,s+1]
      ZAXPY(12, -conj(K[s+1]), 48+f_out+24*(x+vol_4d_cb*(s+1)),48+f_out +24*(x+vol_4d_cb*s));
      // downer fout[x,s] += factL* fout[x,ls-1]
      ZAXPY(12, factL, 48+f_out +12+24*(x+vol_4d_cb*(ls-1)), 48+f_out + 12 + 24*(x+vol_4d_cb*s));

      // upper fout[x,s] += -conj(K[s+1]) fout[x,s+1]
      ZAXPY(12, -conj(K[s+1]), 72+f_out+24*(x+vol_4d_cb*(s+1)),72+f_out +24*(x+vol_4d_cb*s));
      // downer fout[x,s] += factL* fout[x,ls-1]
      ZAXPY(12, factL, 72+f_out +12+24*(x+vol_4d_cb*(ls-1)), 72+f_out + 12 + 24*(x+vol_4d_cb*s));
    }
    //FIXME : which is faster ?
    //fact *= inv_two_kappa;
    if(s>0) factL /= -conj(K[s-1]);
  }

  // s = ls - 1
  //-------------
#pragma omp parallel for
  for(int x=0; x<vol_4d_cb; x+=4) {
    // downer  fou[x,ls-1] *= inv_d_last 
    VEC_TIMESEQU_COMPLEX(f_out+12+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_COMPLEX(24+f_out+12+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_COMPLEX(48+f_out+12+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_COMPLEX(72+f_out+12+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
  }

  //DEBUG_MOBIUS_DSLASH("dag1 second part %e\n", time_elapse());
  
  //DiracOp::CGflops+=2*2*vol_4d_cb*local_ls*12;
  DiracOp::CGflops+=vol_4d_cb*(ls*96-48);
}






//----------------------------------------
#if 0
void zmobius_m5inv(Vector *inout,
	       Float mass,
	       int dag,
	       Zmobus *mobius_lib_arg)
{
  //! fixme most likely slow
  IFloat* f_inout  = (IFloat *) inout;
  
  const int local_ls = GJP.SnodeSites();
  const int s_nodes = GJP.Snodes();
  const int global_ls = local_ls * s_nodes;
  const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = mobius_lib_arg->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;

  // prepare kappa_ratio
  Complex *kappa_ratio = new Complex[global_ls];
  Complex kapR=1.0;
  for(int s=0; s < global_ls; s++){
    kappa_ratio[s] =
      mobius_lib_arg->zmobius_kappa_b[s]/mobius_lib_arg->zmobius_kappa_c[s];
        printf("kapR %d %e %e\n", s, kappa_ratio[s].real(), kappa_ratio[s].imag());
    kapR *= kappa_ratio[s];
  }
  if( fabs(kapR.imag()/kapR.real()) > 1e-10 || kapR.real() <=0){
    ERR.General("","zmobius_m5inv(V*,F,I,Zmobus*)",
		"assumption fail kapR %e %e\n", kapR.real(), kapR.imag());

  }
  kapR = exp(log(kapR.real())/global_ls);
  
  {
    Complex factR=1.0;
    Complex factL=1.0;
    Float* fio = f_inout + ls_stride;
    for(int s=1; s < local_ls; s++){
      int glb_s = s + local_ls*s_node_coor;

      factR *= kappa_ratio[s]/kapR;
      factL *= kappa_ratio[s-1]/kapR;

      Complex fR=factR;
      Complex fL=factL;      
      if(dag){ fR=conj(fR); fL=conj(fL);}
      Complex i_fR=1.0/fR;
      Complex i_fL=1.0/fL;

      
      for(int x=0; x<vol_4d_cb; x++){

	//P_R
	if(! dag)
	  vecTimesEquComplex((Complex*)fio, i_fR, 12);
	else
	  vecTimesEquComplex((Complex*)fio, fR, 12);
	fio  +=  12;

	//P_L
	if(! dag)
	  vecTimesEquComplex((Complex*)fio, fL, 12);
	else
	  vecTimesEquComplex((Complex*)fio, i_fL, 12);	  
	fio  +=  12;      
      }
    }
  }
  Float kappa_b_save =   mobius_lib_arg->mobius_kappa_b;
  Float kappa_c_save =   mobius_lib_arg->mobius_kappa_c;
  mobius_lib_arg->mobius_kappa_b = kapR.real();
  mobius_lib_arg->mobius_kappa_c = 1.0;
  

  if(dag==0)
    zmobius_m5inv_dag0(inout, mass, mobius_lib_arg);
  else 
    zmobius_m5inv_dag1(inout, mass, mobius_lib_arg);  


  mobius_lib_arg->mobius_kappa_b = kappa_b_save;
  mobius_lib_arg->mobius_kappa_c = kappa_c_save;

  {
    Complex factR=1.0;
    Complex factL=1.0;
    Float* fio = f_inout + ls_stride;
    for(int s=1; s < local_ls; s++){
      int glb_s = s + local_ls*s_node_coor;

      factR *= kappa_ratio[s]/kapR;
      factL *= kappa_ratio[s-1]/kapR;

      Complex fR=factR;
      Complex fL=factL;      
      if(dag){ fR=conj(fR); fL=conj(fL);}
      Complex i_fR=1.0/fR;
      Complex i_fL=1.0/fL;

      for(int x=0; x<vol_4d_cb; x++){

	//P_R
	if(! dag)
	  vecTimesEquComplex((Complex*)fio, fR, 12);
	else
	  vecTimesEquComplex((Complex*)fio, i_fR, 12);
	fio  +=  12;

	//P_L
	if(! dag)
	  vecTimesEquComplex((Complex*)fio, i_fL, 12);
	else
	  vecTimesEquComplex((Complex*)fio, fL, 12);	  
	fio  +=  12;      
      }
    }
  }

  delete [] kappa_ratio;
}
#else
void zmobius_m5inv(Vector *inout,
		   Float mass,
		   int dag,
		   Zmobus *mobius_lib_arg,
		   Complex* K)
{
  if(dag==0)
    zmobius_m5inv_cmplx_dag0(inout, mass, mobius_lib_arg,K);
  else
    zmobius_m5inv_cmplx_dag1(inout, mass, mobius_lib_arg,K);

}
#endif
//----------------------------------------

void zmobius_m5inv(Vector *out, Vector *in,
		  Float mass,
		  int dag,
		   Zmobus *mobius_lib_arg,
		   Complex* K)
{
  const int vol_4d_cb = mobius_lib_arg->vol_4d / 2;
  const size_t f_size = 24 * mobius_lib_arg->ls * vol_4d_cb;
  moveFloat( (Float*)out, (Float*)in, f_size );

  //zmobius_m5inv_dag0(out, mass, mobius_lib_arg);

  //return;
  
  if(dag==0)
    zmobius_m5inv_cmplx_dag0(out, mass, mobius_lib_arg,K);
  else
    zmobius_m5inv_cmplx_dag1(out, mass, mobius_lib_arg,K);

}



CPS_END_NAMESPACE
