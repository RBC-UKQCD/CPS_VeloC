#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mobius/noarch/mobius_m.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
// mobius_m.C
//
// mobius_m is the fermion matrix.  
// The in, out fields are defined on the checkerboard lattice
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/dwf.h>
#include<util/zmobius.h>
#include<util/gjp.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/dirac_op.h>
#include<util/time_cps.h>

#include "blas-subs.h"

CPS_START_NAMESPACE


#include "zmobius_m-orig.h"
#include "zmobius_m-sym1.h"
#include "zmobius_m-sym1-MIT.h"
#include "zmobius_m-sym2.h"
#include "zmobius_m-sym2-MIT.h"
#include "zmobius_m-sym3.h"

void  zmobius_m(Vector *out, 
		Matrix *gauge_field, 
		Vector *in, 
		Float mass, 
		Zmobus *mobius_lib_arg)
{
  Float dtime;
  switch( mobius_lib_arg-> pc_type ){
  case   ZMOB_PC_ORIG:
//dtime = -dclock(true);
    zmobius_m_orig(out, gauge_field, in, mass, mobius_lib_arg);
//dtime += dclock(true);
//    print_flops("","zmobius_m-orig()",0,dtime);
    break;
  case   ZMOB_PC_SYM1:
//dtime = -dclock(true);
    zmobius_m_sym1(out, gauge_field, in, mass, mobius_lib_arg);
//dtime += dclock(true);
    break;
  case   ZMOB_PC_SYM2:
    zmobius_m_sym2(out, gauge_field, in, mass, mobius_lib_arg);
    break;
  case   ZMOB_PC_SYM1_MIT:
    zmobius_m_sym1_MIT(out, gauge_field, in, mass, mobius_lib_arg);
    break;
  case   ZMOB_PC_SYM2_MIT:
    zmobius_m_sym2_MIT(out, gauge_field, in, mass, mobius_lib_arg);
    break;
  case   ZMOB_PC_SYM3:
    zmobius_m_sym3(out, gauge_field, in, mass, mobius_lib_arg);
    break;
  default:
    ERR.NotImplemented("","zmobius_mdag(...)");
  }    
}








CPS_END_NAMESPACE
