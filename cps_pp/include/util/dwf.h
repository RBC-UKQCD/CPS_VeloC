#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// dwf.h
//
// C header file for the dwf fermion library
//
//------------------------------------------------------------------

#ifndef INCLUDED_DWF_H
#define INCLUDED_DWF_H

CPS_END_NAMESPACE
#include <vector> //std::vector
#include <util/wilson.h>
#include <util/vector.h>
#include <util/rcomplex.h>
#include <comms/sysfunc_cps.h> //for UniqueID()
CPS_START_NAMESPACE

//------------------------------------------------------------------
// External declarations
//------------------------------------------------------------------

extern int dwfso_wire_map[];
// Set in dwf_int. For a given index 0-1 corresponding to
// S+, S-  it gives the corresponding wire.

inline void DEBTIZB(char* str, Vector* vp, int len)
{
  Float norm = vp->NormSqGlbSum(len);
  if(!UniqueID()) printf("DEBTIZB %s %.16e\n", str, norm);
}


//------------------------------------------------------------------
// Type definitions
//------------------------------------------------------------------

// pre/post-conditioning variations
enum ZMobiusPCType {
  ZMOB_PC_ORIG,
  ZMOB_PC_SYM1,
  ZMOB_PC_SYM2,
  ZMOB_PC_SYM1_MIT,
  ZMOB_PC_SYM2_MIT,
  ZMOB_PC_SYM3};


// The Dwf structure typedef
typedef struct{
  IFloat dwf_kappa;    // 1/[2*(5-dwf_height)]
  IFloat mobius_kappa_b;    // 1/[2*(b*(4-dwf_height)+1)]
  IFloat mobius_kappa_c;    // 1/[2*(c*(4-dwf_height)-1)]
#if 0
  Rcomplex *zmobius_kappa_b;    // 1/[2*(b*(4-dwf_height)+1)]
  Rcomplex *zmobius_kappa_c;    // 1/[2*(c*(4-dwf_height)-1)]
  Rcomplex *zmobius_kappa_ratio;//    kappa_b / kappa_c
#else
  std::vector < Rcomplex > zmobius_b;    // b[i]
  std::vector < Rcomplex > zmobius_c;    // c[i]
  std::vector < Rcomplex > zmobius_kappa_b;    // 1/[2*(b*(4-dwf_height)+1)]
  std::vector < Rcomplex > zmobius_kappa_c;    // 1/[2*(c*(4-dwf_height)-1)]
  std::vector < Rcomplex > zmobius_kappa_ratio;//    kappa_b / kappa_c
#endif
  // pre/post-conditioning for zmobius
  ZMobiusPCType pc_type;

  int vol_4d;         // The 4-dimensional volume   
  int ls;             // The extent of the 5th direction
  IFloat *frm_tmp1;    // Pointer to temorary fermion field 1
  IFloat *frm_tmp2;    // Pointer to temorary fermion field 2
  IFloat *frm_tmp3;    // Pointer to temorary fermion field 3

  // doesn't need any more
  //IFloat *frm_tmp3;    // Pointer to temorary fermion field 3
  
  IFloat *comm_buf;    // Communications buffer (12 words for dwfso)
  Wilson *wilson_p;   // Pointer to the wilson structure
} Dwf;


//------------------------------------------------------------------
// Function declarations
//------------------------------------------------------------------


//------------------------------------------------------------------
// This routine performs all initializations needed before dwf 
// library funcs are called. It sets the addressing related arrays 
// and reserves memory for the needed temporary buffers. It only 
// needs to be called only once at the begining of the program 
// (or after a dwf_end call)before any number of calls to dwf 
// functions are made.
//------------------------------------------------------------------
void dwf_init(Dwf *dwf_p);             // pointer to Dwf struct.


//------------------------------------------------------------------
// This routine frees any memory reserved by dwf_init
//------------------------------------------------------------------
void dwf_end(Dwf *dwf_p);              // pointer to Dwf struct.


//------------------------------------------------------------------
// dwf_mdagm M^dag M where M is the fermion matrix.
// The in, out fields are defined on the checkerboard lattice.
// <out, in> = <dwf_mdagm*in, in> is returned in dot_prd.
//------------------------------------------------------------------
void  dwf_mdagm(Vector *out, 
		Matrix *gauge_field, 
		Vector *in, 
		Float *dot_prd,
		Float mass,
		Dwf *dwf_lib_arg);


// spectrum shifted version : (H-mu)(H-mu)
void    dwf_mdagm_shift(Vector *out, 
			Matrix *gauge_field, 
			Vector *in, 
			Float *dot_prd,
			Float mass,
			Dwf *dwf_lib_arg,
			Float shift);

//------------------------------------------------------------------
// dwf_dslash is the derivative part of the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
// cb = 0/1 <--> even/odd checkerboard of in field.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//------------------------------------------------------------------
void dwf_dslash(Vector *out, 
		Matrix *gauge_field, 
		Vector *in,
		Float mass,
		int cb,
		int dag,
		Dwf *dwf_lib_arg);


//------------------------------------------------------------------
// dwf_m is the fermion matrix.  
// The in, out fields are defined on the checkerboard lattice.
//------------------------------------------------------------------
void  dwf_m(Vector *out, 
	    Matrix *gauge_field, 
	    Vector *in,
	    Float mass,
	    Dwf *dwf_lib_arg);
void  dwf_m_4dpc(Vector *out, 
		 Matrix *gauge_field, 
		 Vector *in,
		 Float mass,
		 Dwf *dwf_lib_arg);
void  dwf_mdag_4dpc(Vector *out, 
		    Matrix *gauge_field, 
		    Vector *in,
		    Float mass,
		    Dwf *dwf_lib_arg);


//------------------------------------------------------------------
// dwf_mdag is the dagger of the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
//------------------------------------------------------------------
void dwf_mdag(Vector *out, 
	      Matrix *gauge_field, 
	      Vector *in,
	      Float mass,
	      Dwf *dwf_lib_arg);
void  dwf_mdag_4dpc(Vector *out, 
		    Matrix *gauge_field, 
		    Vector *in,
		    Float mass,
		    Dwf *dwf_lib_arg);



//------------------------------------------------------------------
// dwf_dslash_4 is the derivative part of the space-time part of
// the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
// cb = 0/1 <--> even/odd checkerboard of in field.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//------------------------------------------------------------------
void dwf_dslash_4(Vector *out, 
		  Matrix *gauge_field, 
		  Vector *in,
		  int cb,
		  int dag,
		  Dwf *dwf_lib_arg);



//------------------------------------------------------------------
// dwf_dslash_5_plus is the derivative part of the 5th direction
// part of the fermion matrix. This routine accumulates the result
// on the out field 
// The in, out fields are defined on the checkerboard lattice.
// The action of this operator is the same for even/odd
// checkerboard fields because there is no gauge field along
// the 5th direction.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//------------------------------------------------------------------
void dwf_dslash_5_plus(Vector *out, 
		       Vector *in,
		       Float mass,
		       int dag,
		       Dwf *dwf_lib_arg);
// subroutine for possible better cache management  (CJ)
void dwf_dslash_5_plus_slice(Vector *out, 
		       Vector *in,
		       Float mass,
		       int dag,
		       Dwf *dwf_lib_arg,
     		       int s_slice);
void dwf_dslash_5_plus_start(Vector *out, 
		       Vector *in,
		       Float mass,
		       int dag,
		       Dwf *dwf_lib_arg);
void dwf_dslash_all(Vector *out, 
		Matrix *gauge_field, 
		Vector *in,
		Float mass,
		int cb,
		int dag,
		Dwf *dwf_lib_arg);



void dwf_m5inv(Vector *out,
	       Vector *in, 
	       Float mass,
	       int dag,
	       Dwf *dwf_lib_arg);

// in-place version
void dwf_m5inv(Vector *inout, 
	       Float mass,
	       int dag,
	       Dwf *dwf_lib_arg);



void ReflectAndMultGamma5( Vector *out, const Vector *in,  int nodevol, int ls);


#endif

CPS_END_NAMESPACE
