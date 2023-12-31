#ifdef USE_SSE
#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /space/cvs/cps/cps++/src/util/dirac_op/d_op_dwf/sse/dwf_init_sse.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
// 11/26/97
//
// dwf_int:
//
// This routine performs all initializations needed before dwf
// func are called. It sets the addressing related arrays and 
// reserves memory for the needed temporary buffers. It only needs
// to be called once at the begining of the program (or after a 
// dwf_end call) before any number of calls to dwf funcs are made.
//
// WARNING:
//
// This set of routines will work only if the node sublattices have
// even number of sites in each direction.
//
//------------------------------------------------------------------
// ======================================================================
/*


  
    * * * * * * *  N O T I C E  * * * * * * * *  

    On  2012-01-01  

    This code is now temporaliry modifed for the 4dim even/odd preconditioning.

    The modified parts code is distinguished by following preprocessor flag.


    
*/

#define FOUR_DIM_EVEN_ODD_PRECONDITION

// ======================================================================

CPS_END_NAMESPACE
#include<util/dwf.h>
#include<util/gjp.h>
#include<util/smalloc.h>
#include<util/verbose.h>
#include<util/error.h>
CPS_START_NAMESPACE



void dwf_init(Dwf *dwf_p)
{
  char *cname = " ";
  char *fname = "dwf_init(Dwf*)";
  VRB.Func(cname,fname);

//------------------------------------------------------------------
// Do initializations before the wilson library can be used
// Initialization involve memory allocation.
//------------------------------------------------------------------
  static Wilson wilson_struct;
  dwf_p->wilson_p = &wilson_struct;
  wilson_init(dwf_p->wilson_p);

//------------------------------------------------------------------
// Allocate memory for two temporary fermion checkerboard fields  
//------------------------------------------------------------------
  size_t f_size = 24 * GJP.VolNodeSites() * GJP.SnodeSites() / 2; 

//  dwf_p->frm_tmp1 = (IFloat *) smalloc(f_size*sizeof(IFloat));
  dwf_p->frm_tmp1 = (IFloat *) smalloc(cname,fname, "frm_tmp1", f_size*sizeof(IFloat));

  dwf_p->frm_tmp2 = (IFloat *) smalloc(cname,fname, "frm_tmp2", f_size*sizeof(IFloat));
  dwf_p->frm_tmp3 = (IFloat *) smalloc(cname,fname, "frm_tmp3", f_size*sizeof(IFloat));
//  if(dwf_p->frm_tmp2 == 0)
//    ERR.Pointer(cname,fname, "frm_tmp2");
//  VRB.Smalloc(cname,fname,
//	      "frm_tmp2", dwf_p->frm_tmp2, f_size*sizeof(IFloat));


#if 1
  dwf_p->frm_tmp3 = (IFloat *) smalloc(f_size*sizeof(IFloat));
  if(dwf_p->frm_tmp3 == 0)
    ERR.Pointer(cname,fname, "frm_tmp3");
  VRB.Smalloc(cname,fname,
	      "frm_tmp3", dwf_p->frm_tmp3, f_size*sizeof(IFloat));
#endif
  
 //------------------------------------------------------------------
// Allocate memory for a 12 word communications buffer needed
// for the spread-out case.
//------------------------------------------------------------------
  dwf_p->comm_buf = (IFloat *) smalloc(12 * sizeof(IFloat));
  if(dwf_p->comm_buf == 0)
    ERR.Pointer(cname,fname, "comm_buf");
  VRB.Smalloc(cname,fname,
	      "comm_buf", dwf_p->comm_buf, 12*sizeof(IFloat));


//------------------------------------------------------------------
// Set the dwf coefficients
//------------------------------------------------------------------
  dwf_p->vol_4d = GJP.VolNodeSites();
  dwf_p->ls = GJP.SnodeSites();
  dwf_p->dwf_kappa = 
    1.0 / ( 2 * ( 4 + GJP.DwfA5Inv() - GJP.DwfHeight() ) );
  dwf_p->mobius_kappa_b = 
    1.0 / ( 2 * (GJP.Mobius_b()*( 4 - GJP.DwfHeight() ) + GJP.DwfA5Inv()) );
  dwf_p->mobius_kappa_c = 
    1.0 / ( 2 * (GJP.Mobius_c()*( 4 - GJP.DwfHeight() ) - GJP.DwfA5Inv()) );

}


CPS_END_NAMESPACE
#endif
