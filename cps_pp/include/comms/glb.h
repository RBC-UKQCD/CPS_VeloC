#ifndef INCLUDED_GLB_H
#define INCLUDED_GLB_H

#include<stdint.h>
#include<config.h>
#include<precision.h>
//-------------------------------------------------------------------
/*!\file
  \brief  Declarations of collective communications routines

*/
/*
 *  glb.h
 */

#include <util/vector.h>
#include <util/lattice.h>
//#include <comms/nga_reg.h>
CPS_START_NAMESPACE

/*! \defgroup collectivecomms Collective communications routines
  \ingroup comms      
  @{ */  

//--------------------------------------------------------------
//! Sums a floating point number over all nodes in a 4-dimensinal grid.
//--------------------------------------------------------------
extern void glb_sum(Float * float_p);
extern void glb_sum_gimp(Float * float_p);

//--------------------------------------------------------------
//! Sums a floating point number over all nodes in a 5-dimensinal grid.
// Sum over all nodes of the "virtual" 5-dimensional volume.
// {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes(), GJP.Snodes}
// Relevant for spread-out DWF (GJP.s_nodes not 1) only.
//--------------------------------------------------------------
extern void glb_sum_five(Float * float_p);

//--------------------------------------------------------------
//! Sums a floating point number over all nodes along a single direction.
//--------------------------------------------------------------
extern void glb_sum_dir(Float * float_p, int dir);
//! Sums a vector of floating point numbers over all nodes along a single direction.
extern void glb_sum_multi_dir(const Float * float_p, const int dir, const int len);
extern void glb_sum_multi_dir(LatData &dat, const int dir);
//! Sums a Matrix over all nodes along a single direction.
extern void glb_sum_matrix_dir(Matrix * float_p, int dir);

//--------------------------------------------------------------
//! Finds the maximum floating point number over all nodes.
//--------------------------------------------------------------
extern void glb_max(Float * float_p);

//--------------------------------------------------------------
//! Finds the minimum floating point number over all nodes.
//--------------------------------------------------------------
extern void glb_min(Float * float_p);

//--------------------------------------------------------------
//! Sums a vector of floating point numbers over all nodes in a hyperplane
//--------------------------------------------------------------
extern void 
slice_sum(Float * float_p, int blcklength, int dir);
extern int glb_sum(long * float_p, const long n_elem);
extern int glb_sum(uint32_t * float_p, const long n_elem);
extern int glb_sum(int * float_p, const long n_elem);
extern int glb_sum(double * float_p, const long n_elem);

/*! @} */


extern unsigned int local_checksum(Float * float_p, int len);
extern unsigned int global_checksum(Float * float_p, int len);
extern unsigned int test_checksum(Float * float_p, int len);

CPS_END_NAMESPACE
#endif
