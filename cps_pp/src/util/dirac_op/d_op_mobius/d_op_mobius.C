#include <config.h>
#include <stdio.h>
#include <alg/no_arg.h>
#include <alg/common_arg.h>
#include <alg/alg_plaq.h>

CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOpMobius class methods.

*/
//------------------------------------------------------------------
//
// d_op_mobius.C
//
// DiracOpMobius is derived from the DiracOp base class. 
// DiracOpMobius is the front end for a library that contains
// all Dirac operators associated with Dwf fermions.
//
//------------------------------------------------------------------
  CPS_END_NAMESPACE
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/wilson.h>
#include <util/time_cps.h>
#include <util/dwf.h>
#include <util/mobius.h>
#include <util/timer.h>
#include <comms/glb.h>
#ifdef USE_BLAS
#include "noarch/blas-subs.h"
#endif
  CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!
  Only one instance of this class is allowed to be in existence at
  any time.
  \param latt The lattice on which this Dirac operator is defined
  \param f_field_out A (pointer to) a spin-colour field (optionally). 
  \param f_field_in A (pointer to) a spin-colour field (optionally).
  \param arg Parameters for the solver.
  \param cnv_frm_flag Whether the lattice fields should be converted to
  to a new storage order appropriate for the type of fermion action.
  If this is ::CNV_FRM_NO, then just the gauge field is converted.
  If this is ::CNV_FRM_YES, then the fields \a f_field_out and \a f_field_in
  are also converted: This assumes they are initially in the same order as
  the gauge field.
 */
//------------------------------------------------------------------
#if 0
static Matrix *new_gauge_field;
static Matrix *old_gauge_field;
#endif
DiracOpMobius::DiracOpMobius (Lattice & latt, Vector * f_field_out, Vector * f_field_in, CgArg * arg, CnvFrmType cnv_frm_flg):
DiracOpWilsonTypes (latt,
                    f_field_out, f_field_in, arg, cnv_frm_flg)
{
  cname = "DiracOpMobius";
  const char *fname = "DiracOpMobius(L&,V*,V*,CgArg*,CnvFrmType)";
  VRB.Func (cname, fname);


  //----------------------------------------------------------------
  // Do the necessary conversions
  //----------------------------------------------------------------
#undef PROFILE
#ifdef PROFILE
  Float dtime = -dclock (true);
#endif
  if (cnv_frm == CNV_FRM_YES)
    lat.Convert (DWF_4D_EOPREC_EE, f_out, f_in);
  else if (cnv_frm == CNV_FRM_NO)
    lat.Convert (WILSON);
#ifdef PROFILE
  dtime += dclock ();
  print_flops ("lattice", "Convert()", 0, dtime);
#endif

  //----------------------------------------------------------------
  // Initialize parameters
  //----------------------------------------------------------------
  DiracArg (dirac_arg);

  //----------------------------------------------------------------
  // Initialize the pointer to the initialized Mobius structure
  // (the structure has been initialized by the Lattice::Fmobius
  // constructor.
  //----------------------------------------------------------------
  mobius_lib_arg = lat.FdiracOpInitPtr ();

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // call to mobius_init but is done here again in case the user
  // has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the Fmobius object.
  //----------------------------------------------------------------
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;

  int Ls = GJP.SnodeSites ();
  mobius_arg->ls = Ls;
  mobius_arg->pc_type = GJP.ZMobius_PC_Type ();
  mobius_arg->mobius_kappa_b =
    1.0 / (2 * (GJP.Mobius_b () * (4 - GJP.DwfHeight ()) + GJP.DwfA5Inv ()));
  mobius_arg->mobius_kappa_c =
    1.0 / (2 * (GJP.Mobius_c () * (4 - GJP.DwfHeight ()) - GJP.DwfA5Inv ()));
  VRB.Debug (cname, fname, "pc_type=%d mobius_kappa_(b,c)=%e %e\n",
             mobius_arg->pc_type,
             mobius_arg->mobius_kappa_b, mobius_arg->mobius_kappa_c);
// to use zmobius dslash implementation
  mobius_arg->zmobius_b.resize (Ls);
  mobius_arg->zmobius_c.resize (Ls);
  mobius_arg->zmobius_kappa_b.resize (Ls);
  mobius_arg->zmobius_kappa_c.resize (Ls);
  mobius_arg->zmobius_kappa_ratio.resize (Ls);
  for (int i = 0; i < Ls; i++) {
    mobius_arg->zmobius_b[i] = GJP.Mobius_b ();
    mobius_arg->zmobius_c[i] = GJP.Mobius_c ();
    mobius_arg->zmobius_kappa_b[i] = mobius_arg->mobius_kappa_b;
    mobius_arg->zmobius_kappa_c[i] = mobius_arg->mobius_kappa_c;
    mobius_arg->zmobius_kappa_ratio[i] =
      mobius_arg->zmobius_kappa_b[i] / mobius_arg->zmobius_kappa_c[i];
  }
  size_t f_size_cb = 12 * GJP.VolNodeSites () * GJP.SnodeSites ();
  if (!mobius_arg->frm_tmp3)
    mobius_arg->frm_tmp3 =
      (IFloat *) smalloc (cname, fname, "frm_tmp3",
                          f_size_cb * sizeof (IFloat));
  if (!mobius_arg->frm_tmp2)
    mobius_arg->frm_tmp2 =
      (IFloat *) smalloc (cname, fname, "frm_tmp2",
                          f_size_cb * sizeof (IFloat));
  if (!mobius_arg->frm_tmp1)
    mobius_arg->frm_tmp1 =
      (IFloat *) smalloc (cname, fname, "frm_tmp1",
                          f_size_cb * sizeof (IFloat));

}


//------------------------------------------------------------------
/*!
  If the storage order of any fields was changed by the constructor
  then they are changed to the canonical order by the destructor.
*/
//------------------------------------------------------------------
DiracOpMobius::~DiracOpMobius ()
{
  const char *fname = "~DiracOpMobius()";
  VRB.Func (cname, fname);

  //----------------------------------------------------------------
  // Do the necessary conversions
  //----------------------------------------------------------------
#undef PROFILE
#ifdef PROFILE
  Float dtime = -dclock (true);
#endif
  if (cnv_frm == CNV_FRM_YES)
    lat.Convert (CANONICAL, f_out, f_in);
  else if (cnv_frm == CNV_FRM_NO)
    lat.Convert (CANONICAL);
#ifdef PROFILE
  dtime += dclock ();
  print_flops ("lattice", "Convert()", 0, dtime);
#endif
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  sfree (cname, fname, "frm_tmp1", mobius_arg->frm_tmp1);
  mobius_arg->frm_tmp1 = NULL;
  sfree (cname, fname, "frm_tmp2", mobius_arg->frm_tmp2);
  mobius_arg->frm_tmp2 = NULL;
  sfree (cname, fname, "frm_tmp3", mobius_arg->frm_tmp3);
  mobius_arg->frm_tmp3 = NULL;

}


//------------------------------------------------------------------
// DiracArg(CgArg *arg):
// It sets the dirac_arg pointer to arg and initializes
// mass.
//------------------------------------------------------------------
void DiracOpMobius::DiracArg (CgArg * arg)
{
  dirac_arg = arg;
  mass = dirac_arg->mass;
}


//------------------------------------------------------------------
// MatPcDagMatPc :
// MatPcDagMatPc is the fermion matrix that appears in the HMC 
// evolution. It is a Hermitian matrix.
// The in, out fields are defined on the checkerboard lattice.
// If dot_prd is not 0 then the dot product (on node)
// <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.
//------------------------------------------------------------------
void DiracOpMobius::MatPcDagMatPc (Vector * out, Vector * in, Float * dot_prd)
{

//  Float dtime = -dclock();
  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  // disabling this as zmobius is called - CJ
  //----------------------------------------------------------------

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  mobius_mdagm (out, gauge_field, in, dot_prd, mass, (Dwf *) mobius_lib_arg);
//  dtime += dclock();
//  print_time(cname,"MatPcDag()",dtime);
}

//------------------------------------------------------------------
// MatPcDagMatPc :
// MatPcDagMatPc is the fermion matrix that appears in the HMC 
// evolution. It is a Hermitian matrix.
// The in, out fields are defined on the checkerboard lattice.
// If dot_prd is not 0 then the dot product (on node)
// <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.
//
// When dirac_arg->eigen_shift is non-zero, it shift the spectrum of matrix:
//    MatPcDagMatPc = H^2  ->  (H-shift)(H-shift)
// where H = Gamma_5 MatPc
//
// For other fermions, one could also implement similar shifts.
// For wilson, H = gamma_5 MatPC .
//------------------------------------------------------------------
void DiracOpMobius::MatPcDagMatPcShift (Vector * out,
                                        Vector * in, Float * dot_prd)
{

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------

  const Float shift = dirac_arg->eigen_shift;

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------

  // we still check if shift is really needed
//  if (shift == 0.0) {
//    mobius_mdagm (out, gauge_field, in, dot_prd, mass, (Dwf *) mobius_lib_arg);
//  } else 
  {

    //mobius_arg->eigen_shift = dirac_arg->eigen_shift;

    mobius_mdagm_shift (out,
                        gauge_field,
                        in, dot_prd, mass, (Dwf *) mobius_lib_arg, shift);
  }

}


//------------------------------------------------------------------
// Dslash(Vector *out, Vector *in, ChkbType cb, DagType dag) :
// Dslash is the derivative part of the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
// cb is not used.
//------------------------------------------------------------------
void DiracOpMobius::Dslash (Vector * out, Vector * in, ChkbType cb, DagType dag)
{

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------
  VRB.Debug (cname, "Dslash", "mass=%e\n", mass);
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  mobius_dslash (out, gauge_field, in, mass, cb, dag, (Dwf *) mobius_lib_arg);
}

void DiracOpMobius::Dslash (Vector * out, Vector * in, DagType dag)
{

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------
  VRB.Debug (cname, "Dslash", "mass=%e\n", mass);
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  mobius_unprec (out, gauge_field, in, mass,
//              cb,
                 dag, (Dwf *) mobius_lib_arg);
}

//------------------------------------------------------------------
/*!
  The preconditioned matrix connects sites of even parity.
  \param out The resulting vector.
  \param in The vector to be multiplied.
*/
//------------------------------------------------------------------
void DiracOpMobius::MatPc (Vector * out, Vector * in)
{

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  const char *fname = "MatPc(*V,*V)";
  VRB.Func (cname, fname);
  mobius_m (out, gauge_field, in, mass, (Dwf *) mobius_lib_arg);
}

//------------------------------------------------------------------
/*!
  The preconditioned matrix connects sites of even parity.
  \param out The resulting vector.
  \param in The vector to be multiplied.
*/
//------------------------------------------------------------------
void DiracOpMobius::MatPcDag (Vector * out, Vector * in)
{

  const char *fname = "MatPcDag(*V,*V)";
  VRB.Func (cname, fname);
  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------
#if 0
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  mobius_arg->ls = GJP.SnodeSites ();
  mobius_arg->mobius_kappa_b =
    1.0 / (2 * (GJP.Mobius_b () * (4 - GJP.DwfHeight ()) + GJP.DwfA5Inv ()));
  mobius_arg->mobius_kappa_c =
    1.0 / (2 * (GJP.Mobius_c () * (4 - GJP.DwfHeight ()) - GJP.DwfA5Inv ()));
#endif

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
//  IFloat *tmp = (IFloat *) in;
  mobius_mdag (out, gauge_field, in, mass, (Dwf *) mobius_lib_arg);
//  tmp = (IFloat *) out;
//  VRB.FuncEnd (cname, fname);
}


//------------------------------------------------------------------
// int MatInv(Vector *out, Vector *in, 
//            Float *true_res, PreserveType prs_in);
// The inverse of the unconditioned Dirac Operator 
// using Conjugate gradient.
// If true_res !=0 the value of the true residual is returned
// in true_res.
// *true_res = |src - MatPcDagMatPc * sol| / |src|
// prs_in is used to specify if the source
// in should be preserved or not. If not the memory usage
// is less by half the size of a fermion vector.
// The function retPRESERVE_Nsurns the total number of CG iterations.
//------------------------------------------------------------------
int DiracOpMobius::MatInv (Vector * out_p,
                           Vector * in_p, Float * true_res, PreserveType prs_in)
{
  std::vector < Vector * >out;
  out.push_back (out_p);
  std::vector < Vector * >in;
  in.push_back (in_p);
  std::vector < Float > res;
  res.push_back (0.);
  std::vector < int >iters;
  iters.push_back (0);
//used later
#define DO_PRESERVE
//  assert (prs_in == PRESERVE_NO);       // should be fixed 
  this->MatInv (out, in, res, iters);
  *true_res = res[0];
  return iters[0];
}

void DiracOpMobius::MatInv (std::vector < Vector * >out,
                            std::vector < Vector * >in,
                            std::vector < Float > &true_res,
                            std::vector < int >&iters)
{
  char *fname = "MatInv(V**,V**,F*,I*,prs)";
  VRB.Func (cname, fname);
  int len = out.size ();
  VRB.Result (cname, fname, "mass=%e len=%d\n", mass, len);
  assert (in.size () == len);
  true_res.resize (len);
  iters.resize (len);

  Zmobus *mobius_arg = (Zmobus *) mobius_lib_arg;
  // (b(4-M)+a5)
  const int local_ls = GJP.SnodeSites ();
  const int s_nodes = GJP.Snodes ();
  const int global_ls = local_ls * s_nodes;
  const int s_node_coor = GJP.SnodeCoor ();
  const int vol_4d_cb = ((Zmobus *) mobius_lib_arg)->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;
  size_t temp_size = GJP.VolNodeSites () * lat.FsiteSize () / 2;
  assert (temp_size == lat.half_size);
  Float two_kappa_b = 2. * mobius_arg->mobius_kappa_b;

  Float *res = true_res.data ();

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  Float norm;
  std::vector < Vector * >temp;
  Vector *temp2;
  std::vector < Vector * >temp3;
  Vector *save_in;
  std::vector < Vector * >odd_in, even_in, odd_out, even_out;

  temp2 =
    (Vector *) smalloc (cname, fname, "temp2", temp_size * sizeof (Float));

  // points to the even part of fermion source 
  // CJ: probably need to move odd to in front. Keeping it as is for now
  //

  for (int i = 0; i < len; i++) {
    temp.push_back ((Vector *)
                    smalloc (cname, fname, "temp", temp_size * sizeof (Float)));
    odd_in.push_back ((Vector *) ((IFloat *) in[i] + temp_size));
    even_in.push_back ((Vector *) ((IFloat *) in[i] + 0));
    odd_out.push_back ((Vector *) ((IFloat *) out[i] + temp_size));
    even_out.push_back ((Vector *) ((IFloat *) out[i] + 0));
  }

#ifdef DO_PRESERVE
  // save source
  PreserveType prs_in = PRESERVE_YES;
  for (int i = 0; i < len; i++)
    if (prs_in == PRESERVE_YES) {
      temp3.push_back ((Vector *) smalloc (cname, fname, "temp3",
                                           2 * temp_size * sizeof (Float)));
      moveMem ((IFloat *) temp3[i], (IFloat *) in[i],
               2 * temp_size * sizeof (IFloat));
    }
#endif

  //-------------------
  // Preconditioning source and guess vector
  //-------------------

  // for source
  switch (mobius_arg->pc_type) {
  case ZMOB_PC_ORIG:
  case ZMOB_PC_SYM1:
  case ZMOB_PC_SYM2:
  case ZMOB_PC_SYM3:
    for (int i = 0; i < len; i++) {
      zmobius_m5inv (temp[i], even_in[i], mass, DAG_NO, mobius_arg,
                     mobius_arg->zmobius_kappa_ratio.data ());

      zmobius_dslash_4 (temp2, gauge_field, temp[i], CHKB_EVEN, DAG_NO,
                        mobius_arg, mass);
      zmobius_zvectTimesV1PlusV2 (temp[i], mobius_arg->zmobius_kappa_b.data (),
                                  temp2, odd_in[i], local_ls, ls_stride,
                                  s_node_coor);
    }
    break;
  default:
    ERR.NotImplemented (cname, fname);
    break;
  }
  // for guess vector
  switch (mobius_arg->pc_type) {
  case ZMOB_PC_ORIG:
  case ZMOB_PC_SYM1:
    // need nothing
    break;
  case ZMOB_PC_SYM2:
  case ZMOB_PC_SYM3:
    // Apply M5 to out for sym2 preconditioning
    for (int i = 0; i < len; i++) {
      moveFloat ((IFloat *) temp2, (IFloat *) odd_out[i], temp_size);
      zmobius_kappa_dslash_5_plus_cmplx (odd_out[i], temp2, mass, DAG_NO,
                                         mobius_arg,
                                         mobius_arg->zmobius_kappa_ratio.
                                         data ());
    }
    break;
  default:
    ERR.NotImplemented (cname, fname);
    break;
  }

  switch (mobius_arg->pc_type) {
  case ZMOB_PC_SYM1:
    for (int i = 0; i < len; i++) {
      zmobius_m5inv (temp2, temp[i], mass, DAG_NO, mobius_arg,
                     mobius_arg->zmobius_kappa_ratio.data ());
      moveFloat ((IFloat *) temp[i], (IFloat *) temp2, temp_size);
    }
    break;
  case ZMOB_PC_SYM3:
  case ZMOB_PC_SYM2:
  case ZMOB_PC_ORIG:
    break;
  default:
    ERR.NotImplemented (cname, fname);
    break;
  }

//  int iter;
  switch (dirac_arg->Inverter) {
  case FAKE:
//    MatPcDag (out, temp);
    for (int i = 0; i < len; i++)
      lat.RandGaussVector (odd_out[i], 1.0, 1, FIVE_D);
    break;
  case CG:
  case CG_FIXED_ITER:
    for (int i = 0; i < len; i++)
      MatPcDag (odd_in[i], temp[i]);
#ifdef USE_QUDA
    QudaInvert (odd_out, odd_in, true_res, iters, 1);
//    for (int i = 0; i < len; i++) 
//      iters[i] = QudaInvert (odd_out[i], odd_in[i], res + i, 1);
#else
    for (int i = 0; i < len; i++) {
      iters[i] = InvCg (odd_out[i], odd_in[i], res + i);
    }
#endif
    break;
  case BICGSTAB:
    for (int i = 0; i < len; i++) {
#ifdef USE_QUDA
      iters[i] = QudaInvert (odd_out[i], temp[i], res + i, 0);
#else
      iters[i] =
        BiCGstab (odd_out[i], temp[i], 0.0, dirac_arg->bicgstab_n, res + i);
#endif
    }
    break;
  case LOWMODEAPPROX:
    for (int i = 0; i < len; i++)
      MatPcDag (odd_in[i], temp[i]);
    InvLowModeApprox (odd_out, odd_in, dirac_arg->fname_eigen, dirac_arg->neig);
    break;
  case CG_LOWMODE_DEFL:
    for (int i = 0; i < len; i++)
      MatPcDag (odd_in[i], temp[i]);
    InvLowModeApprox (odd_out, odd_in, dirac_arg->fname_eigen, dirac_arg->neig);
//void DiracOpMobius::MatInv (std::vector < Vector * >out,
//                            std::vector < Vector * >in,
//                            std::vector < Float > &true_res,
//                            std::vector < int >&iters)
#ifdef USE_QUDA
#if 1
    QudaInvert (odd_out, odd_in, true_res, iters, 1);
#else
    for (int i = 0; i < len; i++) {
      iters[i] = QudaInvert (odd_out[i], odd_in[i], res + i, 1);
    }
#endif
#else
    for (int i = 0; i < len; i++) {
      iters[i] = InvCg (odd_out[i], odd_in[i], res + i);
    }
#endif
//    exit(-43);
    break;
  default:
    ERR.General (cname, fname, "InverterType %d not implemented\n",
                 dirac_arg->Inverter);
    break;
  }


  if (VRB.IsActivated (VERBOSE_DEBUG_LEVEL)) 
  {
    // check solution
    for (int i = 0; i < len; i++) {
      norm = odd_out[i]->NormSqGlbSum (temp_size);
      if (!UniqueID ())
        printf ("Norm %d out %.14e\n", i, norm);
      norm = odd_in[i]->NormSqGlbSum (temp_size);
      if (!UniqueID ())
        printf ("Norm %d in %.14e\n", i, norm);
      MatPcDagMatPc (temp[0], odd_out[0]);
      norm = temp[i]->NormSqGlbSum (temp_size);
      if (!UniqueID ())
        printf ("Norm %d MatPcDagMatPc*out %.14e\n", i, norm);
    }
  }
#ifdef DO_PRESERVE
  // restore source
  if (prs_in == PRESERVE_YES) {
    for (int i = 0; i < len; i++)
      moveMem ((IFloat *) in[i], (IFloat *) temp3[i],
               2 * temp_size * sizeof (IFloat));
  }
#endif

  switch (mobius_arg->pc_type) {
  case ZMOB_PC_SYM3:
  case ZMOB_PC_SYM2:
    for (int i = 0; i < len; i++) {
      zmobius_m5inv (temp[i], odd_out[i], mass, DAG_NO, mobius_arg,
                     mobius_arg->zmobius_kappa_ratio.data ());
      moveFloat ((IFloat *) odd_out[i], (IFloat *) temp[i], temp_size);

      // Below is the same original postconditioning (dare to write again for clarity)
      zmobius_dslash_4 (temp[i], gauge_field, odd_out[i], CHKB_ODD, DAG_NO,
                        mobius_arg, mass);
      zmobius_zvectTimesEquComplex (temp[i],
                                    mobius_arg->zmobius_kappa_b.data (),
                                    local_ls, ls_stride, s_node_coor);
      zmobius_m5inv (even_out[i], temp[i], mass, DAG_NO, mobius_arg,
                     mobius_arg->zmobius_kappa_ratio.data ());
      zmobius_m5inv (temp[i], even_in[i], mass, DAG_NO, mobius_arg,
                     mobius_arg->zmobius_kappa_ratio.data ());
      even_out[i]->VecAddEquVec (temp[i], temp_size);
    }
    break;
  case ZMOB_PC_SYM1:
  case ZMOB_PC_ORIG:
    for (int i = 0; i < len; i++) {
      zmobius_dslash_4 (temp[i], gauge_field, odd_out[i], CHKB_ODD, DAG_NO,
                        mobius_arg, mass);
      zmobius_zvectTimesEquComplex (temp[i],
                                    mobius_arg->zmobius_kappa_b.data (),
                                    local_ls, ls_stride, s_node_coor);
      zmobius_m5inv (even_out[i], temp[i], mass, DAG_NO, mobius_arg,
                     mobius_arg->zmobius_kappa_ratio.data ());
      zmobius_m5inv (temp[i], even_in[i], mass, DAG_NO, mobius_arg,
                     mobius_arg->zmobius_kappa_ratio.data ());
      even_out[i]->VecAddEquVec (temp[i], temp_size);
    }
    break;
  default:
    ERR.NotImplemented (cname, fname);
    break;
  }

  for (int i = 0; i < len; i++) {
    Float *temp_f = (Float *) out[i];
    Float *in_f = (Float *) in[i];
    Float max_diff = 0.;
#pragma omp parallel for default(shared) reduction(max:max_diff)
    for (size_t i = 0; i < 2 * temp_size; ++i) {
      *(temp_f + i) *= two_kappa_b;
//      Float temp = *(temp_f+i)-*(in_f+i);
//        temp = temp*temp;
//        if(temp > max_diff) max_diff=temp;
    }
//  VRB.Result(cname,fname,"max_diff[%d]=%e\n",i, sqrt(max_diff));
  }

  for (int i = 0; i < len; i++) {
    sfree (cname, fname, "temp", temp.back ());
    temp.pop_back ();
#ifdef DO_PRESERVE
    if (prs_in == PRESERVE_YES) {
      sfree (cname, fname, "temp3", temp3.back ());
      temp3.pop_back ();
    }
#endif
  }
  sfree (cname, fname, "temp2", temp2);


//  return iter;
}


//------------------------------------------------------------------
// Overloaded function is same as original 
// but true_res=0.
//------------------------------------------------------------------
int DiracOpMobius::MatInv (Vector * out, Vector * in, PreserveType prs_in)
{
  return MatInv (out, in, 0, prs_in);
}


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in and out = f_out.
//------------------------------------------------------------------
int DiracOpMobius::MatInv (Float * true_res, PreserveType prs_in)
{
  return MatInv (f_out, f_in, true_res, prs_in);
}


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in, out = f_out, true_res=0.
//------------------------------------------------------------------
int DiracOpMobius::MatInv (PreserveType prs_in)
{
  return MatInv (f_out, f_in, 0, prs_in);
}


//------------------------------------------------------------------
// Mat(Vector *out, Vector *in) :
// Mat is the unpreconditioned fermion matrix.  
// Mat works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpMobius::Mat (Vector * out, Vector * in)
{
  const char *fname = "Mat(V*,V*)";
  VRB.Func (cname, fname);

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  mobius_arg->ls = GJP.SnodeSites ();
  mobius_arg->mobius_kappa_b =
    1.0 / (2 * (GJP.Mobius_b () * (4 - GJP.DwfHeight ()) + GJP.DwfA5Inv ()));
  mobius_arg->mobius_kappa_c =
    1.0 / (2 * (GJP.Mobius_c () * (4 - GJP.DwfHeight ()) - GJP.DwfA5Inv ()));
  Float kappa = mobius_arg->mobius_kappa_b;
  Float minus_kappa = -kappa;
  Float kappa_ratio = mobius_arg->mobius_kappa_b / mobius_arg->mobius_kappa_c;

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  int temp_size = GJP.VolNodeSites () * lat.FsiteSize () / 2;

  // points to the even part of fermion source 
  Vector *odd_in = (Vector *) ((IFloat *) in + temp_size);
  // points to the even part of fermion solution
  Vector *odd_out = (Vector *) ((IFloat *) out + temp_size);
  // temp
  Vector *frm_tmp2 = (Vector *) mobius_arg->frm_tmp2;

  //odd part
  //mobius_dslash_4(out, gauge_field, odd_in, CHKB_EVEN, DAG_NO, mobius_arg, mass);
  mobius_dslash_4 (out, gauge_field, odd_in, CHKB_ODD, DAG_NO, mobius_arg,
                   mass);
  out->VecTimesEquFloat (minus_kappa, temp_size);

  // intialize to zero since using the "plus-equal version"
  for (int i = 0; i < temp_size; i++) {
    *((IFloat *) frm_tmp2 + i) = 0.0;
  }
  mobius_dslash_5_plus (frm_tmp2, in, mass, 0, mobius_arg);

  fTimesV1PlusV2 ((IFloat *) frm_tmp2, kappa_ratio, (IFloat *) frm_tmp2,
                  (IFloat *) in, temp_size);
  out->VecAddEquVec (frm_tmp2, temp_size);

  //even part
  mobius_dslash_4 (odd_out, gauge_field, in, CHKB_EVEN, DAG_NO, mobius_arg,
                   mass);
  odd_out->VecTimesEquFloat (minus_kappa, temp_size);
  // intialize to zero since using the "plus-equal version"
  for (int i = 0; i < temp_size; i++) {
    *((IFloat *) frm_tmp2 + i) = 0.0;
  }
  mobius_dslash_5_plus (frm_tmp2, odd_in, mass, 0, mobius_arg);
  fTimesV1PlusV2 ((IFloat *) frm_tmp2, kappa_ratio, (IFloat *) frm_tmp2,
                  (IFloat *) odd_in, temp_size);
  odd_out->VecAddEquVec (frm_tmp2, temp_size);

}


void DiracOpMobius::Dminus (Vector * out, Vector * in)
{
  const char *fname = "Dminus(V*,V*)";
  VRB.Func (cname, fname);

  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  VRB.Debug (cname, fname, "in=%p out=%p mobius_arg=%p\n", in, out, mobius_arg);
  Float kappa_c_inv_div2 =
    0.5 * (2 * (GJP.Mobius_c () * (4 - GJP.DwfHeight ()) - GJP.DwfA5Inv ()));
  VRB.Debug (cname, fname, "kappa_c_inv_div2=%e\n", kappa_c_inv_div2);

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  size_t temp_size = GJP.VolNodeSites () * lat.FsiteSize () / 2;
  out->VecZero (temp_size * 2);

  // points to the odd part of fermion source 
  Vector *odd_in = (Vector *) ((IFloat *) in + temp_size);
  // points to the odd part of fermion solution
  Vector *odd_out = (Vector *) ((IFloat *) out + temp_size);

  VRB.Flow (cname, fname, "odd_in=%p odd_out=%p mobius_arg=%p\n", odd_in,
            odd_out, mobius_arg);
  mobius_dminus (out, gauge_field, odd_in, CHKB_ODD, DAG_NO, mobius_arg);
  VRB.Flow (cname, fname, "mobius_dminus()\n");
  mobius_dminus (odd_out, gauge_field, in, CHKB_EVEN, DAG_NO, mobius_arg);
  // out = -(c*D_W-1)*in (= 1 for DWF)
  // CJ: Was there a sign change?? 
  fTimesV1PlusV2 ((IFloat *) out, kappa_c_inv_div2, (IFloat *) in,
                  (IFloat *) out, 2 * temp_size);
  out->VecTimesEquFloat (-1.0, 2 * temp_size);
  VRB.FuncEnd (cname, fname);

}


//------------------------------------------------------------------
// MatDag(Vector *out, Vector *in) :
// MatDag is the unpreconditioned fermion matrix.  
// MatDag works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpMobius::MatDag (Vector * out, Vector * in)
{
  const char *fname = "MatDag(V*,V*)";
  VRB.Func (cname, fname);

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.MobiusA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  mobius_arg->ls = GJP.SnodeSites ();
  mobius_arg->mobius_kappa_b =
    1.0 / (2 * (GJP.Mobius_b () * (4 - GJP.DwfHeight ()) + GJP.DwfA5Inv ()));
  mobius_arg->mobius_kappa_c =
    1.0 / (2 * (GJP.Mobius_c () * (4 - GJP.DwfHeight ()) - GJP.DwfA5Inv ()));
  Float kappa = mobius_arg->mobius_kappa_b;
  Float minus_kappa = -kappa;
  Float kappa_ratio = mobius_arg->mobius_kappa_b / mobius_arg->mobius_kappa_c;

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  int temp_size = GJP.VolNodeSites () * lat.FsiteSize () / 2;

  // points to the even part of fermion source 
  Vector *odd_in = (Vector *) ((IFloat *) in + temp_size);
  // points to the even part of fermion solution
  Vector *odd_out = (Vector *) ((IFloat *) out + temp_size);
  // temp
  Vector *frm_tmp2 = (Vector *) mobius_arg->frm_tmp2;

  //odd part
  mobius_dslash_4 (out, gauge_field, odd_in, CHKB_ODD, DAG_YES, mobius_arg,
                   mass);
  //mobius_dslash_4(out, gauge_field, odd_in, CHKB_EVEN, DAG_YES, mobius_arg, mass);
  ERR.General (cname, fname,
               "TIZB I do not understand why kappa instead of minus_kappa here\n");
  out->VecTimesEquFloat (kappa, temp_size);
  // intialize to zero since using the "plus-equal version"
  for (int i = 0; i < temp_size; i++) {
    *((IFloat *) frm_tmp2 + i) = 0.0;
  }
  mobius_dslash_5_plus (frm_tmp2, in, mass, DAG_YES, mobius_arg);
  fTimesV1PlusV2 ((IFloat *) frm_tmp2, kappa_ratio, (IFloat *) frm_tmp2,
                  (IFloat *) in, temp_size);
  out->VecAddEquVec (frm_tmp2, temp_size);

  //even part
  //mobius_dslash_4(odd_out, gauge_field, in, CHKB_ODD, DAG_YES, mobius_arg, mass);
  mobius_dslash_4 (odd_out, gauge_field, in, CHKB_EVEN, DAG_YES, mobius_arg,
                   mass);
  odd_out->VecTimesEquFloat (kappa, temp_size);
  // intialize to zero since using the "plus-equal version"
  for (int i = 0; i < temp_size; i++) {
    *((IFloat *) frm_tmp2 + i) = 0.0;
  }
  mobius_dslash_5_plus (frm_tmp2, odd_in, mass, DAG_YES, mobius_arg);
  fTimesV1PlusV2 ((IFloat *) frm_tmp2, kappa_ratio, (IFloat *) frm_tmp2,
                  (IFloat *) odd_in, temp_size);
  odd_out->VecAddEquVec (frm_tmp2, temp_size);

}


//------------------------------------------------------------------
// MatHerm(Vector *out, Vector *in) :
// MatHerm is gamma5*R*Mat.
// MatHerm works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpMobius::MatHerm (Vector * out, Vector * in)
{
  const char *fname = "MatHerm(V*,V*)";
  VRB.Func (cname, fname);

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  mobius_arg->ls = GJP.SnodeSites ();
  mobius_arg->mobius_kappa_b =
    1.0 / (2 * (GJP.Mobius_b () * (4 - GJP.DwfHeight ()) + GJP.DwfA5Inv ()));
  mobius_arg->mobius_kappa_c =
    1.0 / (2 * (GJP.Mobius_c () * (4 - GJP.DwfHeight ()) - GJP.DwfA5Inv ()));

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  int temp_size = GJP.VolNodeSites () * lat.FsiteSize ();
  Vector *temp =
    (Vector *) smalloc (cname, fname, "temp", temp_size * sizeof (Float));
  if (temp == 0)
    ERR.Pointer (cname, fname, "temp");

  Mat (out, in);
  lat.Freflex (temp, out);
  MultGamma (out, temp, 15, GJP.VolNodeSites () * GJP.SnodeSites ());

  sfree (cname, fname, "temp", temp);

}


//------------------------------------------------------------------
/*!
  \pre This method is to be used when the instance of this object has been
  created with \a f_field_out and \a f_field_in pointers to spin-colour
  vectors defined over the whole lattice. 

  \param chi A spin-colour vector defined on odd parity lattice sites.

  \post The vector \a f_field_out is \f$ (1+D)\chi \f$

  and the vector \a f_field_in is \f$ (D^\dagger-\kappa^2 M)\chi \f$

  where \e M is the odd-even preconditioned fermion matrix connecting odd to
  odd parity sites and \e D is the hopping term connecting odd to
  even parity sites. Recall that \a chi is defined on odd sites only.
  The new vectors are in odd-even order.
*/
//------------------------------------------------------------------

void DiracOpMobius::CalcHmdForceVecs (Vector * chi)
{
  const char *fname = "CalcHmdForceVecs(V*)";
  static Timer timer (cname, fname);
  timer.start (false);
  VRB.Func (cname, fname);

  if (f_out == 0)
    ERR.Pointer (cname, fname, "f_out");

  if (f_in == 0)
    ERR.Pointer (cname, fname, "f_in");

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.MobiusA5Inv(), GJP.MobiusHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------

//------------------------------------------------------------------
// f_out stores (chi,rho), f_in stores (psi,sigma)
//------------------------------------------------------------------

  Vector *chi_new, *rho, *psi, *sigma;

  size_t f_size_cb = 12 * GJP.VolNodeSites () * GJP.SnodeSites ();

  chi_new = f_out;

  chi_new->CopyVec (chi, f_size_cb);

  psi = f_in;

  MatPc (psi, chi);
  CalcHmdForceVecs (f_in, f_out, chi_new, psi);
  timer.stop (false);
  return;
}

#if 0
static void print (const char *name, Vector * v, size_t f_size)
{
  Float *v_p = (Float *) v;
  if (!UniqueID ())
    printf ("%s: %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g norm=%0.12g\n", name,
            v_p[0], v_p[1], v_p[2], v_p[3], v_p[4], v_p[5],
            v->NormSqGlbSum (f_size));
}
#endif

void DiracOpMobius::CalcHmdForceVecs (Vector * v1, Vector * v2, Vector * phi1,
                                      Vector * phi2)
{
  const char *fname = "CalcHmdForceVecs(V*,V*,V*,V*)";
  VRB.Func (cname, fname);
  static Timer timer (cname, fname);
  timer.start (false);

  if (f_out == 0)
    ERR.Pointer (cname, fname, "f_out");

  if (f_in == 0)
    ERR.Pointer (cname, fname, "f_in");

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.MobiusA5Inv(), GJP.MobiusHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------

//------------------------------------------------------------------
// f_out stores (chi,rho), f_in stores (psi,sigma)
//------------------------------------------------------------------

  size_t f_size_cb = 12 * GJP.VolNodeSites () * GJP.SnodeSites ();
//  Vector *phi=phi1, *eta=phi2;
  Vector *tmp1 =
    (Vector *) smalloc (cname, fname, "tmp1", f_size_cb * sizeof (Float));
  Vector *tmp2 =
    (Vector *) smalloc (cname, fname, "tmp2", f_size_cb * sizeof (Float));
  Vector *tmp3 =
    (Vector *) smalloc (cname, fname, "tmp3", f_size_cb * sizeof (Float));
  Vector *tmp4 =
    (Vector *) smalloc (cname, fname, "tmp4", f_size_cb * sizeof (Float));

  Vector *v1_o = (Vector *) ((Float *) v1 + f_size_cb);
  Vector *v1_e = (Vector *) ((Float *) v1 + 0);
  Vector *v2_o = (Vector *) ((Float *) v2 + f_size_cb);
  Vector *v2_e = (Vector *) ((Float *) v2 + 0);

// Float dtime = -dclock(true);
  switch (mobius_arg->pc_type) {
  case ZMOB_PC_SYM1:
    mobius_m5inv (v1_o, phi1, mass, DAG_YES, (Dwf *) mobius_lib_arg);
    phi1->CopyVec (v1_o, f_size_cb);
    break;
  case ZMOB_PC_SYM2:
    mobius_m5inv (v1_o, phi2, mass, DAG_NO, (Dwf *) mobius_lib_arg);
    phi2->CopyVec (v1_o, f_size_cb);
    break;
  case ZMOB_PC_ORIG:
    break;
  default:
    ERR.NotImplemented (cname, fname);
    break;
  }
  v1_o->CopyVec (phi1, f_size_cb);
//  dtime += dclock(true);
//  print_flops(fname,"m5inv+CopyVec()",0,dtime);
//  dtime =-dclock();



//  Float kb = mobius_arg->mobius_kappa_b;
  if (0) {
//  v2_o->VecTimesEquFloat(-kb*kb,f_size_cb);
    VRB.Debug (cname, fname, "phi1=%g\n", phi1->NormSqGlbSum (f_size_cb));
    VRB.Debug (cname, fname, "v1_o=%g\n", v1_o->NormSqGlbSum (f_size_cb));
    phi1->print ("phi1", f_size_cb);
    lat.Dump ("phi1", phi1, Odd);
    v1_o->print ("v1_o", f_size_cb);
    lat.Dump ("v1_o", v1_o, Odd);


//  v2_o->CopyVec(phi2,f_size_cb);
    VRB.Debug (cname, fname, "phi2=%g\n", phi2->NormSqGlbSum (f_size_cb));
    phi2->print ("phi2", f_size_cb);
    lat.Dump ("phi2", phi2, Even);      // hack to match the test program

//  Vector *rho = (Vector *) ((Float *) v2 + 0);

//  VRB.Debug(cname,fname,"mass=%g Gparity=%d\n",mass,GJP.Gparity());
//  Dslash (v2_e, v2_o, CHKB_ODD, DAG_NO);
//  dtime += dclock(true);
//  print_flops(fname,"Dump()",0,dtime);
//  dtime =-dclock();
  }
  static Timer timer_dslash (fname, "mobius_dslash");
  timer_dslash.start (false);
  mobius_dslash_4 (tmp1, gauge_field, phi2, CHKB_ODD, DAG_NO,
                   (Dwf *) mobius_lib_arg, mass);
  timer_dslash.stop (false);
//  dtime += dclock(true);
//  print_flops(fname,"mobius_dslash",0,dtime);
//  dtime =-dclock();
// lat.Dump("Meophi2",tmp1,Even); // (-2.)* BFM
  static Timer timer_m5inv (fname, "mobius_m5inv");
  timer_m5inv.start (false);
  mobius_m5inv (tmp2, tmp1, mass, DAG_NO, (Dwf *) mobius_lib_arg);
  timer_m5inv.stop (false);
//  dtime += dclock(true);
//  print_flops(fname,"mobius_m5inv",0,dtime);
//  dtime =-dclock();
//  lat.Dump("MeeInvMeophi2",tmp2,Even); // Factor?
  static Timer timer_Booee (fname, "mobius_Booee");
  timer_Booee.start (false);
  mobius_Booee (v2_e, tmp2, DAG_NO, (Dwf *) mobius_lib_arg, mass);
//  lat.Dump("v2_e",v2_e,Even); // Factor?

  VRB.Debug (cname, fname, "v2_e=%g\n", v2_e->NormSqGlbSum (f_size_cb));

  mobius_Booee (v2_o, phi2, DAG_NO, (Dwf *) mobius_lib_arg, mass);
  timer_Booee.stop (false);
//  dtime += dclock(true);
//  print_flops(fname,"mobius_Booee",0,dtime);
//  dtime =-dclock();
//  lat.Dump("v2_o",v2_o,Odd); // Factor?

  timer_dslash.start (false);
  mobius_dslash_4 (tmp4, gauge_field, phi1, CHKB_ODD, DAG_YES,
                   (Dwf *) mobius_lib_arg, mass);
  timer_dslash.stop (false);
//  lat.Dump("MeoDagphi1",tmp4,Even); // (-2.)* BFM
  timer_m5inv.start (false);
  mobius_m5inv (v1_e, tmp4, mass, DAG_YES, (Dwf *) mobius_lib_arg);
  timer_m5inv.stop (false);
//  dtime += dclock(true);
//  print_flops(fname,"mobius_m5inv()",0,dtime);
//  dtime =-dclock();
//  lat.Dump("v1_e",v1_e,Even); // Factor?

//  VRB.Debug(cname,fname,"v1_e=%g\n",v1_e->NormSqGlbSum(f_size_cb));
//  v1_e->print("v1_e",f_size_cb);

  sfree (cname, fname, "tmp4", tmp4);
  sfree (cname, fname, "tmp3", tmp3);
  sfree (cname, fname, "tmp2", tmp2);
  sfree (cname, fname, "tmp1", tmp1);
//  dtime += dclock(true);
//  print_flops(fname,"print()",0,dtime);

  timer.stop (false);

  return;
}

//------------------------------------------------------------------
// DiracOpGlbSum(Float *): 
// The global sum used by InvCg. If s_nodes = 1
// it is the usual global sum. If s_nodes > 1 it
// is the 5-dimensional globals sum glb_sum_five.
//------------------------------------------------------------------
void DiracOpMobius::DiracOpGlbSum (Float * float_p)
{
//  if(GJP.Snodes() == 1) {
//    glb_sum(float_p);
//  }
//  else {
  glb_sum_five (float_p);
//  }
}



#ifndef USE_BLAS
#define MOVE_FLOAT( pa, pb, n )  moveFloat(pa, pb, n)
#define VEC_TIMESEQU_FLOAT(py, fact, n ) vecTimesEquFloat( py, fact, n)
#define AXPY(n, fact, px, py)  fTimesV1PlusV2(py, fact, px, py, n)
#else
#define MOVE_FLOAT( pa, pb, n )  cblas_dcopy(n, pb, 1, pa, 1)
#define VEC_TIMESEQU_FLOAT(py, fact, n ) cblas_dscal( n,  fact, py,1 )
#define AXPY(n, fact, px, py)  cblas_daxpy(n, fact, px,1,py,1)
#endif


//------------------------------------------------------------------
// RitzMat(Vector *out, Vector *in) :
// RitzMat is the base operator used in in Ritz.
// RitzMat works on the full or half lattice.
// The in, out fields are defined on the full or half lattice.
//------------------------------------------------------------------
void DiracOpMobius::RitzMat (Vector * out, Vector * in)
{
  const char *fname = "RitzMat(V*,V*)";
  VRB.Func (cname, fname);

  //printf("single ritzmat %d\n",dirac_arg->RitzMatOper);
  switch (dirac_arg->RitzMatOper) {
    // Now always call MatPcDagMatPcShift even for MATPCDAG_MATPC case
    //
    //case MATPCDAG_MATPC:
    //MatPcDagMatPc(out, in);
    //break;
  case MATPCDAG_MATPC:
  case MATPCDAG_MATPC_SHIFT:
    MatPcDagMatPcShift (out, in);
    break;
  case MATPC_HERM:
    MatPcHerm (out, in);
    break;

  case MAT_HERM:
  case MATDAG_MAT:
    MatDagMat (out, in);
    break;
  case NEG_MATPCDAG_MATPC:
    MatPcDagMatPc (out, in);
    out->VecNegative (out, RitzLatSize ());
    break;

  case NEG_MATDAG_MAT:
    MatDagMat (out, in);
    out->VecNegative (out, RitzLatSize ());
    break;

  default:
    ERR.General (cname, fname, "RitzMatOper %d not implemented\n",
                 dirac_arg->RitzMatOper);
  }

#if 0
  //debug
  const int size = RitzLatSize ();      // this is number of Float, f_size .
  Complex deb = out->CompDotProductGlbSum (in, size);
  Float d_deb = in->NormSqGlbSum (size);
  printf ("single Ritz %e %e %e\n", deb.real () / d_deb, deb.imag () / d_deb,
          d_deb);
  //debug
#endif
}



// PolynomialAccerelation
//
//  Q = [ -2 RitzMat + (alpha + beta) ] / [ alpha - beta ]
//
//  Output:  out =  T_n(Q) in
//
//  T_0 = 1,    T_1 = Q
//   T_{n+1}(Q) =  2 Q T_n(Q)  - T_{n-1}(Q)
//  
// Calling virtual RitzMat(V*,V*)
//
//   alpha = param[0]^2
//   beta = ( param[1] + fabs(eigend_shift) )^2
//
void DiracOpMobius::RitzMat (Vector * out, Vector * in,
                             MatrixPolynomialArg * cheby_arg)
{
  const char *fname = "RitzMat(V*,V*,MatrixPolyArg*)";
  VRB.Func (cname, fname);

  double time_start = dclock ();

  //debug const
//  const Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  const Float shift = dirac_arg->eigen_shift;

  const int Npol = cheby_arg->Npol;
  const int size = RitzLatSize ();      // this is number of Float, f_size .
  //
  // Q = 2 / (alpha-beta)  RitzMat  -   (alpha+beta)/(alpha-beta)
  // 2 Q =   c1  (  c0 Ddag D  -   1 )
  //  c1 = 2 (alpha+beta)/(alpha-beta),   c0 =  2 / (alpha+beta)
  //
  const Float alpha = pow (cheby_arg->params.params_val[0], 2);
  const Float beta = pow (cheby_arg->params.params_val[1] + fabs (shift), 2);
  //printf("alpha=%e beta=%e\n", alpha,beta);

  const Float c1 = 2.0 * (alpha + beta) / (alpha - beta);
  const Float c0 = 2.0 / (alpha + beta);

  Vector *tmp = (Vector *) cheby_arg->tmp1;
  Vector *tmp2 = (Vector *) cheby_arg->tmp2;


  //  tmp2 =  T_0 v = v = in
  //tmp2 -> CopyVec(in, size);
  MOVE_FLOAT ((Float *) tmp2, (Float *) in, size);
  //  tmp =  T_1 v = Q v = Q in
  //  QV = 0.5* (2Q)V = 0.5 c1 ( c0 Ddag D - 1)

  RitzMat (tmp, in);

#if 0
  tmp->VecTimesEquFloat (c0, size);
  tmp->VecMinusEquVec (in, size);
  tmp->VecTimesEquFloat (0.5 * c1, size);
  //  tmp =  0.5 c1 ( c0 DdagD in - in )
#else
  VEC_TIMESEQU_FLOAT ((Float *) tmp, 0.5 * c1 * c0, size);
  AXPY (size, -0.5 * c1, (Float *) in, (Float *) tmp);
#endif

  // debug
  //out->CopyVec(tmp,size);
  //printf("cheby %f %f\n", alpha,beta);

  // loop over
  for (int i = 2; i <= Npol; ++i) {
#if 0
    // out = 2 Q tmp
    RitzMat (out, tmp);

    out->VecTimesEquFloat (c0, size);
    out->VecMinusEquVec (tmp, size);
    out->VecTimesEquFloat (c1, size);

    // out = out - tmp2
    out->VecMinusEquVec (tmp2, size);
#else
    // out = c1 (  c0 DagD tmp - tmp) -tmp2

    RitzMat (out, tmp);

    VEC_TIMESEQU_FLOAT ((Float *) out, c1 * c0, size);
    AXPY (size, -c1, (Float *) tmp, (Float *) out);
    AXPY (size, -1.0, (Float *) tmp2, (Float *) out);
#endif


    if (i != Npol) {
#if 0
      // tmp2 = tmp
      tmp2->CopyVec (tmp, size);
      // tmp = out
      tmp->CopyVec (out, size);
#else
      // tmp2 = tmp
      Vector *swap_tmp2 = tmp2;
      tmp2 = tmp;
      tmp = swap_tmp2;
      // tmp = out
      tmp->CopyVec (out, size);
#endif
    }
  }


  if (!UniqueID ())
    printf ("mpoly tot = %e\n", dclock () - time_start);

#if 0
  //debug
  Complex deb = out->CompDotProductGlbSum (in, size);
  printf ("debug %e %e %e\n", deb.real (), deb.imag (),
          in->NormSqGlbSum (size));
  //debug
#endif

}


//!! N.B. This overwrites contents of  mobius_arg->frm_tmp2
void DiracOpMobius::MatPcHerm (Vector * out, Vector * in)
{
  const char *fname = "MatPcHerm(V*,V*)";
  VRB.Func (cname, fname);
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  Vector *vtmp = (Vector *) (mobius_arg->frm_tmp1);

  MatPc (vtmp, in);
  ReflectAndMultGamma5 (out, vtmp, mobius_arg->vol_4d / 2, mobius_arg->ls);

}


// specific to dwf 
void ReflectAndMultGamma5 (Vector * out, const Vector * in, int nodevol, int ls)
{
  const char *fname = "MultGamma5(V*,V*,i)";
  VRB.Func ("", fname);
  for (int s = 0; s < ls; ++s) {
    IFloat *p = (IFloat *) out + 24 * nodevol * s;
    IFloat *q = (IFloat *) in + 24 * nodevol * (ls - 1 - s);
    for (int n = 0; n < nodevol; ++n) {
      int i;
      for (i = 0; i < 12; ++i)
        *p++ = *q++;

      for (i = 0; i < 12; ++i)
        *p++ = -*q++;
    }
  }

}

void HermicianDWF_ee (Vector * vtmp, Vector * evec, Float mass,
                      Lattice * lattice, Vector * Apsi)
{
  CgArg cg_arg;
  cg_arg.mass = mass;
  cg_arg.RitzMatOper = MATPC_HERM;      // could be MATPCDAG_MATPC;
  DiracOpDwf dop (*lattice, 0, 0, &cg_arg, CNV_FRM_NO);

  dop.MatPc (Apsi, evec);
  ReflectAndMultGamma5 (vtmp, Apsi, GJP.VolNodeSites () / 2, GJP.SnodeSites ());
}

CPS_END_NAMESPACE
