// C
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// C++
#include <cassert>
#include <string>
#include <vector>

// QMP
#include <qmp.h>
#include <cps.h>

// Measurement package
//#include "eigcg.h"
#include "my_util.h"
#include "prop_container.h"
#include "run_2pion.h"
#include "run_bk.h"
#include "run_k2pipi.h"
#include "run_kl3.h"
#include "run_meson.h"
#include "run_mres.h"
#include "run_omega.h"
#include "run_prop.h"
#include "twisted_bc.h"

static const char *cname = "";

USING_NAMESPACE_CPS
using namespace std;

CommonArg common_arg;
DoArg do_arg;
DoArgExt do_ext;
LanczosArg lanczos_arg;
FixGaugeArg fix_gauge_arg;
MeasArg meas_arg;
QPropWArg lqpropw_arg_e;
QPropWArg lqpropw_arg;
QPropWArg sqpropw_arg_e;
QPropWArg sqpropw_arg;
QPropW4DBoxArg box_arg;

// Integer arrays, for setting time locations of exact propagators
LocArg l_ut_loc_e; // untwisted light
LocArg l_ut_loc; // untwisted light
//LocArg l_tw_loc; // twisted light 
//LocArg s_ut_loc_e; // untwisted strange
//LocArg s_ut_loc; // untwisted strange
//LocArg s_tw_loc; // twisted strange

// d quark momentum for K -> pi pi, for 48^3 this should be {1, 1, 1}.
IntArray d_mom_kpp;

// l and s quark twists
FloatArray l_twist_arg;
FloatArray s_twist_arg;
//--------------------------------------------------------------------

//--------------------------------------------------------------------
// DJM: Setup routines common to both interfaces
//

#undef encode_vml
#define encode_vml(arg_name, traj) do{                                  \
        char vml_file[256];                                             \
        sprintf(vml_file, #arg_name".%d", traj);                        \
        if( !arg_name.Encode(vml_file, #arg_name) ){                    \
            ERR.General(cname, fname, #arg_name " encoding failed.\n"); \
        }                                                               \
    }while(0)

#undef decode_vml
#define decode_vml(arg_name)  do{                        \
  if ( ! arg_name.Decode(#arg_name".vml", #arg_name) )   \
  ERR.General(cname, fname, "Bad " #arg_name ".vml.\n"); \
} while(0)  



inline int Chdir(const char* dir)
{
  const char* fname = "Chdir(char*)";

  if(chdir(dir) != 0){
    ERR.General("", fname, "Changing to directory %s failed.\n", dir);
  }

  return 0;
}



static void SetZmobiusPC(int flag){
  switch(flag){
  case 0:
    GJP.ZMobius_PC_Type (ZMOB_PC_ORIG);
    break;
  case 1:
    GJP.ZMobius_PC_Type (ZMOB_PC_SYM1);
    break;
  case 2:
    GJP.ZMobius_PC_Type (ZMOB_PC_SYM2);
    break;
  default:
    ERR.General ("", "SetZmobiusPC()", "%d Not implemented\n", flag);
  }
}






void decode_vml_all(void)
{
  const char *fname = "decode_vml_all()";

  decode_vml(do_arg);
  decode_vml(do_ext);
  decode_vml(meas_arg);
  decode_vml(lqpropw_arg_e);
  decode_vml(lqpropw_arg);
  decode_vml(sqpropw_arg_e);
  decode_vml(sqpropw_arg);
  decode_vml(box_arg);
  decode_vml(fix_gauge_arg);
  decode_vml(lanczos_arg);

  decode_vml(l_ut_loc_e);
  decode_vml(l_ut_loc);
//  decode_vml(l_tw_loc);
//  decode_vml(s_tw_loc);

  decode_vml(d_mom_kpp);

  decode_vml(l_twist_arg);
  decode_vml(s_twist_arg);
#ifdef USE_QUDA
   if ( !QudaParam.Decode("quda_arg.vml","QudaParam") )
        { printf("Bum quda_arg\n"); exit(-1);}
#endif
}

void encode_vml_all(int traj)
{
  const char *fname = "encode_vml_all()";

  encode_vml(do_arg,traj);
  encode_vml(do_ext,traj);
  encode_vml(meas_arg,traj);
  encode_vml(lqpropw_arg_e,traj);
  encode_vml(lqpropw_arg,traj);
  encode_vml(sqpropw_arg_e,traj);
  encode_vml(sqpropw_arg,traj);
  encode_vml(box_arg,traj);
  encode_vml(fix_gauge_arg,traj);
  // decode_vml(l_eigcg_arg);
  // decode_vml(s_eigcg_arg);
  // decode_vml(lanc_arg);

  encode_vml(l_ut_loc_e,traj);
  encode_vml(l_ut_loc,traj);
//  encode_vml(l_tw_loc,traj);
//  encode_vml(s_ut_loc_e,traj);
//  encode_vml(s_ut_loc,traj);
//  encode_vml(s_tw_loc,traj);

  encode_vml(d_mom_kpp,traj);

  encode_vml(l_twist_arg,traj);
  encode_vml(s_twist_arg,traj);
#ifdef USE_QUDA
//  encode_vml(quda_arg,traj);
#endif
}

void load_checkpoint(int traj)
{
  const char *fname = "load_checkpoint()";

  char lat_file[256];
  GnoneFnone lat;

  sprintf(lat_file, "%s.%d", meas_arg.GaugeStem, traj);
  QioArg rd_arg(lat_file, 0.001);
  rd_arg.ConcurIONumber = meas_arg.IOconcurrency;
  ReadLatticeParallel rl;
  rl.read(lat,rd_arg);
  if(!rl.good()) ERR.General(cname,fname,"Failed read lattice %s\n",lat_file);
}

void setup(int argc, char** argv)
{
  const char *fname = "setup()";

  Start(&argc, &argv);

  CommandLine::is(argc, argv);
  Chdir(CommandLine::arg());

#if 0
  if(!QudaParam.Decode(argv[2], "QudaParam")){ ERR.General("", fname, "Bum quda_arg\n"); }
  VRB.Result("", fname, "device %d\n", QudaParam.device);
#endif

  decode_vml_all();
  encode_vml_all(0);
  VRB.Result("", fname, "Read VML files successfully.\n");
  GJP.Initialize(do_arg);
  GJP.InitializeExt(do_ext);
  //LRG.setSerial();
  LRG.Initialize();

  SetZmobiusPC(1); // set as 1, change it if needed

  int nthreads(64);
  if(getenv("OMP_NUM_THREADS")){ 
    nthreads = atoi(getenv("OMP_NUM_THREADS")); 
  } else {
    VRB.Result("", fname, "WARNING: using default %d OMP threads.\n", nthreads);
  }
  GJP.SetNthreads(nthreads);

}

// Shift the locations of where exact propagators are calculated. This
// is done by shifting the following arrays,
//
// l_ut_loc, l_tw_loc
// s_ut_loc, s_tw_loc
// Commenting out until working out what to do
#if 0
void do_shift(int traj)
{
  const char* fname = "do_shift()";
  VRB.Result(cname, fname, "Shifting locations at which we will calculate exact propagators...\n");

  const int t_size = GJP.TnodeSites() * GJP.Tnodes();

  // We shift the exact solutions by a random amount uniformly
  // distributed between [0,T).
//  int shift = drand48() * t_size;
  int shift = 3;

  // Make sure we have the same number on all nodes.
  QMP_broadcast(&shift, sizeof(int));
  assert(shift >= 0 && shift < t_size);

  static int shift_acc = 0;

  shift_acc = (shift_acc + shift) % t_size;
  VRB.Result(cname, fname, "traj = %d, Shift on the exact propagators = %d\n", traj, shift_acc);

  for(unsigned i = 0; i < l_ut_loc_e.locs.locs_len; ++i) {
    l_ut_loc.v.v_val[i] = (l_ut_loc_e.v.v_val[i] + shift) % t_size;
  }
  for(unsigned i = 0; i < l_ut_loc.locs.locs_len; ++i) {
    l_ut_loc.v.v_val[i] = (l_ut_loc.v.v_val[i] + shift) % t_size;
  }
  for(unsigned i = 0; i < l_tw_loc.locs.locs_len; ++i) {
    l_tw_loc.v.v_val[i] = (l_tw_loc.v.v_val[i] + shift) % t_size;
  }
  for(unsigned i = 0; i < s_ut_loc_e.locs.locs_len; ++i) {
    s_ut_loc.v.v_val[i] = (s_ut_loc_e.v.v_val[i] + shift) % t_size;
  }
  for(unsigned i = 0; i < s_ut_loc.locs.locs_len; ++i) {
    s_ut_loc.v.v_val[i] = (s_ut_loc.v.v_val[i] + shift) % t_size;
  }
  for(unsigned i = 0; i < s_tw_loc.locs.locs_len; ++i) {
    s_tw_loc.v.v_val[i] = (s_tw_loc.v.v_val[i] + shift) % t_size;
  }
}
#endif

//--------------------------------------------------------------------

void run_contractions(const AllProp &sprop, const AllProp &stwst,
    const AllProp &lprop, const AllProp &ltwst,
    const string &rdir,
    int traj, PROP_TYPE ptype)
{
  const char *fname = "run_contractions()";
  VRB.Result(cname, fname, "run_contractions for trajectory %d with prop type %d\n", traj, ptype);

  const string trajs = string(".") + tostring(traj);

  //////////////////////////////////////////////////////////////////////
  // 1. meson contractions

  //Greg: Can skip contractions involving twists, so I've commented them out.

  // pion and kaon
  VRB.Result(cname, fname, "Doing untwisted point sink pion contractions\n");
  run_meson_pt(lprop, lprop, GAMMA_5, GAMMA_5, rdir + "/pion-00WP" + trajs, ptype);
  run_meson_pt(lprop, ltwst, GAMMA_5, GAMMA_5, rdir + "/pion-01WP" + trajs, ptype);

  VRB.Result(cname, fname, "Doing untwisted point sink kaon contractions\n");
  run_meson_pt(sprop, lprop, GAMMA_5, GAMMA_5, rdir + "/kaon-00WP" + trajs, ptype);
  run_meson_pt(stwst, lprop, GAMMA_5, GAMMA_5, rdir + "/kaon-10WP" + trajs, ptype);
  run_meson_pt(sprop, ltwst, GAMMA_5, GAMMA_5, rdir + "/kaon-01WP" + trajs, ptype);
  run_meson_wall(lprop, lprop, ID,      ID,      rdir + "/sigma-00WW" + trajs, ptype);

  VRB.Result(cname, fname, "Doing untwisted point sink ss contractions\n");
  run_meson_pt(sprop, sprop, GAMMA_5, GAMMA_5, rdir + "/ss-00WP" + trajs, ptype);

  VRB.Result(cname, fname, "Doing untwisted wall sink pion contractions\n");
  run_meson_wall(lprop, lprop, GAMMA_5, GAMMA_5, rdir + "/pion-00WW" + trajs, ptype);
  run_meson_wall(lprop, ltwst, GAMMA_5, GAMMA_5, rdir + "/pion-01WW" + trajs, ptype);
  VRB.Result(cname, fname, "Doing untwisted wall sink kaon contractions\n");
  run_meson_wall(sprop, lprop, GAMMA_5, GAMMA_5, rdir + "/kaon-00WW" + trajs, ptype);
  run_meson_wall(stwst, lprop, GAMMA_5, GAMMA_5, rdir + "/kaon-10WW" + trajs, ptype);
  run_meson_wall(sprop, ltwst, GAMMA_5, GAMMA_5, rdir + "/kaon-01WW" + trajs, ptype);

  // scalar meson (sigma) contractions.
  VRB.Result(cname, fname, "Doing untwisted light-light scalar contractions\n");
  run_meson_wall(lprop, lprop, ID,      ID,      rdir + "/sigma-00WW" + trajs, ptype);
  run_meson_disc(lprop, lprop, ID,      ID,      rdir + "/sigma-dis-00WW" + trajs, ptype);

  // eta eta' contractions
  //
  // We share the light-light propagator with pion contractions.
  VRB.Result(cname, fname, "Doing untwisted strange-strange pseudoscalar contractions\n");
  run_meson_wall(sprop, sprop, GAMMA_5, GAMMA_5, rdir + "/ss-00WW" + trajs, ptype);
  VRB.Result(cname, fname, "Doing untwisted light-strange pseudoscalar disconnected contractions\n");
  run_meson_disc(lprop, sprop, GAMMA_5, GAMMA_5, rdir + "/ls-dis-00WW" + trajs, ptype);

  // f_K and f_pi measurements
  VRB.Result(cname, fname, "Doing a bunch of local axial current contractions for f_K and f_pi\n");
  run_meson_pt(lprop, lprop, GAMMA_35, GAMMA_5, rdir + "/fp-00WP" + trajs, ptype);
  run_meson_pt(sprop, lprop, GAMMA_35, GAMMA_5, rdir + "/fk-00WP" + trajs, ptype);

  run_meson_pt(lprop, lprop, GAMMA_5, GAMMA_35, rdir + "/fpr-00WP" + trajs, ptype);
  run_meson_pt(sprop, lprop, GAMMA_5, GAMMA_35, rdir + "/fkr-00WP" + trajs, ptype);

  run_meson_pt(lprop, lprop, GAMMA_35, GAMMA_35, rdir + "/ap-00WP" + trajs, ptype);
  run_meson_pt(sprop, lprop, GAMMA_35, GAMMA_35, rdir + "/ak-00WP" + trajs, ptype);

  run_meson_wall(lprop, lprop, GAMMA_35, GAMMA_5, rdir + "/fp-00WW" + trajs, ptype);
  run_meson_wall(sprop, lprop, GAMMA_35, GAMMA_5, rdir + "/fk-00WW" + trajs, ptype);

  // rho meson
  VRB.Result(cname, fname, "Doing a bunch of vector meson (rho) contractions\n");
  run_meson_pt(lprop, lprop, GAMMA_0, GAMMA_0, rdir + "/rho-x-00WP" + trajs, ptype);
  run_meson_pt(lprop, lprop, GAMMA_1, GAMMA_1, rdir + "/rho-y-00WP" + trajs, ptype);
  run_meson_pt(lprop, lprop, GAMMA_2, GAMMA_2, rdir + "/rho-z-00WP" + trajs, ptype);
  run_meson_wall(lprop, lprop, GAMMA_0, GAMMA_0, rdir + "/rho-x-00WW" + trajs, ptype);
  run_meson_wall(lprop, lprop, GAMMA_1, GAMMA_1, rdir + "/rho-y-00WW" + trajs, ptype);
  run_meson_wall(lprop, lprop, GAMMA_2, GAMMA_2, rdir + "/rho-z-00WW" + trajs, ptype);

  //////////////////////////////////////////////////////////////////////
  // 2. Omega baryon
  VRB.Result(cname, fname, "Doing a bunch of omega baryon contractions\n");
  run_omega_pt(sprop, GAMMA_0, rdir + "/sss-x-00WP" + trajs, ptype);
  run_omega_pt(sprop, GAMMA_1, rdir + "/sss-y-00WP" + trajs, ptype);
  run_omega_pt(sprop, GAMMA_2, rdir + "/sss-z-00WP" + trajs, ptype);
  run_omega_pt(sprop, GAMMA_3, rdir + "/sss-t-00WP" + trajs, ptype);
  run_omega_pt(sprop, GAMMA_5, rdir + "/sss-5-00WP" + trajs, ptype);

  //////////////////////////////////////////////////////////////////////
  // 3. Kl3

#if 1
  run_kl3(sprop, lprop, lprop, rdir + "/kl3-00" + trajs, ptype);
  run_kl3(sprop, lprop, ltwst, rdir + "/kl3-01" + trajs, ptype);
//  run_kl3(stwst, lprop, lprop, rdir + "/kl3-10" + trajs, ptype);

  // The following contractions are used for Z_V.
  // Useful for the double ratio method or the UKQCD method.
  run_kl3(lprop, lprop, lprop, rdir + "/zpa-00" + trajs, ptype);
  run_kl3(sprop, lprop, sprop, rdir + "/zka-00" + trajs, ptype);
  run_kl3(lprop, sprop, lprop, rdir + "/zkb-00" + trajs, ptype);
  run_kl3(sprop, sprop, sprop, rdir + "/zss-00" + trajs, ptype);

  // Since line 1 and 3 both carry momentum, are their directions
  // consistent?
  // I think they are.
  run_kl3(ltwst, lprop, ltwst, rdir + "/zpa-11" + trajs, ptype);
  run_kl3(lprop, ltwst, lprop, rdir + "/zpb-11" + trajs, ptype);
  run_kl3(stwst, lprop, stwst, rdir + "/zka-11" + trajs, ptype);
  run_kl3(lprop, stwst, lprop, rdir + "/zkb-11" + trajs, ptype);
#endif

  //////////////////////////////////////////////////////////////////////
  // 4. Bk
  VRB.Result(cname, fname, "Doing B_k contractions\n");
  run_bk(lprop, sprop, lprop, sprop, rdir + "/bk" + trajs, ptype);
}

void run_omega_box_contractions(const AllProp &sprop_box,
                                const string &rdir,
                                int traj, PROP_TYPE ptype)
{
    const char *fname = "run_omega_box_contractions()";

    const string trajs = string(".") + tostring(traj);

    //////////////////////////////////////////////////////////////////////
    // 2. meson contractions

    // eta eta' contractions
    //
    // We share the light-light propagator with pion contractions.
    run_meson_pt  (sprop_box, sprop_box, GAMMA_5, GAMMA_5, rdir + "/box-ss-00WP" + trajs, ptype);
    run_meson_wall(sprop_box, sprop_box, GAMMA_5, GAMMA_5, rdir + "/box-ss-00WW" + trajs, ptype);

    //////////////////////////////////////////////////////////////////////
    // 2. Omega baryon
    run_omega_pt(sprop_box, GAMMA_0, rdir + "/box-sss-x-00WP" + trajs, ptype);
    run_omega_pt(sprop_box, GAMMA_1, rdir + "/box-sss-y-00WP" + trajs, ptype);
    run_omega_pt(sprop_box, GAMMA_2, rdir + "/box-sss-z-00WP" + trajs, ptype);
    run_omega_pt(sprop_box, GAMMA_3, rdir + "/box-sss-t-00WP" + trajs, ptype);
    run_omega_pt(sprop_box, GAMMA_5, rdir + "/box-sss-5-00WP" + trajs, ptype);
}


void run_k2pipi_contractions(const AllProp &sprop,
    const AllProp &uprop,
    const AllProp &dprop,
    const string &rdir,
    int traj, const int mom[3],
    PROP_TYPE ptype)
{
  const string trajs = string(".") + tostring(traj);

  // zero momentum pi pi scattering
  int zmom[3] = {0, 0, 0};
  run_2pionDC(uprop, uprop, rdir + "/2pion000" + trajs, ptype, zmom);

//   Only need the zero-momentum case to compute the S-wave pi-pi scattering length
	for(unsigned i = 0; i < 8; ++i) {
    int p[3];
    p[0] = (i & 1) ? -mom[0] : mom[0];
    p[1] = (i & 2) ? -mom[1] : mom[1];
    p[2] = (i & 4) ? -mom[2] : mom[2];

    const string fn = rdir + "/2pion"
      + tostring(p[0]) + tostring(p[1]) + tostring(p[2])
      + trajs;

    run_2pionDC(uprop, dprop, fn, ptype, p);
  }

  // K->pipi without momentum.
   run_k2pipi(sprop, uprop, uprop, rdir + "/k2pipi-0" + trajs, ptype);

  // K->pipi with momentum.
   run_k2pipi(sprop, uprop, dprop, rdir + "/k2pipi-1" + trajs, ptype);
}

// stw: twisting angle of the strange quark (connecting the operator
// and the kaon).
//
// ltw: twisting angle of the light quark (connecting the operator and
// the pion).
void run_all(Lattice &lat,
    const double stw[4], // strange quark twists, for Kl3
    const double ltw[4], //   light quark twists, for Kl3
    const int mom[3],    // momentum for the d quark, used by K2pipi
    int traj)
{
  const char *fname = "run_all()";
  VRB.Result(cname, fname, "Starting run_all for trajectory %d\n", traj);

  int if_tw=0;

  EigenCache *ecache = NULL;
  const size_t MAX_LEN = 256;

  const char *evec_name = "light_evec";
  const string evec_dir = string("evec.") + tostring(traj);
  int N_evec = lqpropw_arg.cg.neig;
  VRB.Result("",fname,"evec_dir=%s N_evec=%d\n",evec_dir.c_str(),N_evec);
  if(N_evec>0) {

    const int n_fields = GJP.SnodeSites ();
    const size_t f_size_per_site = lat.FsiteSize () / n_fields / 2;     // checkerboarding
    size_t evec_size = (size_t) (GJP.VolNodeSites () / 2) * lat.FsiteSize ();
    if (evec_size != lat.half_size)
      ERR.General ("", fname, "evec_size(%d)!=half_size(%d)\n");

    ecache = new EigenCache ("light_evec");
    int data_size = sizeof (float);
    Float dtime = -dclock();
    ecache->readCompressedBlocks (evec_dir.c_str());
    dtime +=dclock();
    print_time("","readCompressedBlocks()",dtime);
    printf("Node %d: ecache(%s)=%p\n",UniqueID(), evec_dir.c_str(),ecache);
    EigenCacheList.push_back (ecache);
  }

//  sprintf(evec_dir,"evec_tw.%d",traj);
//  evec_dir.string(""); 
//  evec_dir <<"evec_tw."<<traj;



  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // 1. props: strange, strange twisted, light, light twisted

  // exact strange
  AllProp sprop(AllProp::DOUBLE), sprop_e(AllProp::DOUBLE),stwst(AllProp::DOUBLE);
  // exact light
  AllProp lprop_e(AllProp::DOUBLE), ltwst_e(AllProp::DOUBLE);
  // sloppy light, single precision
  AllProp lprop(AllProp::SINGLE), ltwst(AllProp::SINGLE);
  // exact strange, Z3 box source
  //AllProp sprop_box(AllProp::DOUBLE);

  Float dtime0 = dclock();

  //this constructor does the gauge fixing, and the destructor
  //automatically frees the gauge fixing matrices at the end of run_all
  VRB.Result(cname, fname, "Starting gauge fixing for trajectory %d\n", traj);
  //AlgFixGauge fg(lat, &common_arg, &fix_gauge_arg);
  const char *gfix_in = "gfix_lat.in";
  const char *gfix_out = "gfix_lat.out";
  CommonArg com_fg("fg.out");
  AlgFixGauge fg(lat, &com_fg, &fix_gauge_arg);
  FILE *fp = fopen(gfix_in, "r");
  if(!fp){
    fg.run();
  }else{
    fclose(fp);
    QioArg rd_arg(gfix_in);
    fg.run(&rd_arg);
  }
  fg.Save(gfix_out, 0);

  Float dtime01 = dclock();

  VRB.Result(cname, fname, "Starting Z3 box source propagator inversions for trajectory %d\n", traj);
  //run_box_prop(&sprop_box, lat, sqpropw_arg, box_arg, traj, false);

  Float dtime02 = dclock();
#if 0

  VRB.Result(cname, fname, "Doing Z3 box contractions for trajectory %d\n", traj);
  //run_omega_box_contractions(sprop_box, "../resultsPA", traj, PROP_PA);
  //run_omega_box_contractions(sprop_box, "../resultsP", traj, PROP_P);
  //run_omega_box_contractions(sprop_box, "../resultsA", traj, PROP_A);

#endif
  Float dtime1 = dclock();
  print_time("main()","omega_box",dtime1-dtime02);
//  exit(-42);
  bool do_deflation = true;
  if(N_evec<=0) do_deflation = false;
  bool do_mres = true;
  

#if 1
  VRB.Result(cname, fname, "Starting wall source propagator inversions for trajectory %d\n", traj);

  // l untwisted
  VRB.Result(cname, fname, "Doing light untwisted wall source propagators for trajectory %d\n", traj);
  lqpropw_arg.cg.fname_eigen=(char*)"light_evec";
//  int l_max_wall= 2;
  run_wall_prop(&lprop_e, l_ut_loc_e, lat, lqpropw_arg_e, "../resultsEA", traj, do_mres, do_deflation);
  run_wall_prop(&lprop,  l_ut_loc, lat, lqpropw_arg, "../resultsA", traj, do_mres, do_deflation);

  if(if_tw && N_evec>0) {
  const string evec_tw =string("evec_tw.") + tostring(traj);
  VRB.Result("",fname,"evec_dir=%s\n",evec_tw.c_str());
    const int n_fields = GJP.SnodeSites ();
    const size_t f_size_per_site = lat.FsiteSize () / n_fields / 2;     // checkerboarding
    size_t evec_size = (size_t) (GJP.VolNodeSites () / 2) * lat.FsiteSize ();
    if (evec_size != lat.half_size)
      ERR.General ("", fname, "evec_size(%d)!=half_size(%d)\n");
    ecache = new EigenCache ("light_evec_tw");
    int data_size = sizeof (float);
    Float dtime = -dclock();
//    ecache->alloc (N_evec, evec_size, data_size);
//    ecache->read_compressed (evec_tw.c_str());
    ecache->readCompressedBlocks (evec_tw.c_str());
    dtime +=dclock();
    print_time("","readCompressedBlocks()",dtime);
    printf("Node %d: ecache(%s)=%p\n",UniqueID(), evec_tw.c_str(),ecache);
    EigenCacheList.push_back (ecache);
  }
//  cps::sync();

if(if_tw){
  // l twisted
  lqpropw_arg.cg.fname_eigen=(char*)"light_evec_tw";
  twisted_bc(lat, ltw, true);
//  run_wall_prop(&ltwst_e, &ltwst, l_tw_loc, lat, lqpropw_arg, NULL,  traj, true,do_deflation );
  twisted_bc(lat, ltw, false);
}

  EigenCacheListCleanup();

  // s untwisted
  VRB.Result(cname, fname, "Doing strange untwisted wall source propagators for trajectory %d\n", traj);
//  int s_max_wall= 2;
  run_wall_prop(&sprop_e, l_ut_loc_e, lat, sqpropw_arg_e, "../resultsEA", traj, do_mres ,false);
  run_wall_prop(&sprop, l_ut_loc, lat, sqpropw_arg, "../resultsA", traj, do_mres  ,false);

  // s twisted
//  twisted_bc(lat, stw, true);
//  run_wall_prop(NULL, &stwst, s_tw_loc, lat, sqpropw_arg, NULL, traj, false, false );
//  twisted_bc(lat, stw, false);

  Float dtime2 = dclock();

  VRB.Result(cname, fname, "Starting contractions for trajectory %d\n", traj);

  VRB.Result(cname, fname, "Doing contractions with exact light propagators\n");
  run_contractions(sprop_e, stwst, lprop_e, ltwst_e, "../resultsEA",  traj, PROP_A);

  VRB.Result(cname, fname, "Doing contractions with inexact light propagators\n");
  run_contractions(sprop, stwst, lprop, ltwst, "../resultsA",  traj, PROP_A);

  Float dtime3 = dclock();

  ////////////////////////////////////////////////////////////////////////
  // I=2 K to pi pi
  // free unwanted propagators to save some memory.
  ltwst_e.clear();
  ltwst.clear();
  stwst.clear();


  //Greg: Don't need to run K to pi pi stuff, so I've them out

#if 0
  // twisted light for K -> pi pi
  run_mom_prop(&ltwst_e, &ltwst, l_tw_loc, lat, lqpropw_arg, NULL, traj, mom,false);

  Float dtime4 = dclock();

  //run_k2pipi_contractions(sprop, lprop_e, ltwst_e, "../resultsEPA", traj, mom, PROP_PA);
  //run_k2pipi_contractions(sprop, lprop_e, ltwst_e, "../resultsEP",  traj, mom, PROP_P);
  run_k2pipi_contractions(sprop, lprop_e, ltwst_e, "../resultsEA",  traj, mom, PROP_A);

  //run_k2pipi_contractions(sprop, lprop, ltwst, "../resultsPA", traj, mom, PROP_PA);
  //run_k2pipi_contractions(sprop, lprop, ltwst, "../resultsP",  traj, mom, PROP_P);
  run_k2pipi_contractions(sprop, lprop, ltwst, "../resultsA",  traj, mom, PROP_A);

#endif
  Float dtime5 = dclock();

  VRB.Result(cname, fname, "fix gauge    = %17.10e seconds\n", dtime01 - dtime0);
  VRB.Result(cname, fname, "box prop     = %17.10e seconds\n", dtime02 - dtime01);
  VRB.Result(cname, fname, "box contract = %17.10e seconds\n", dtime1 - dtime02);
  VRB.Result(cname, fname, "kl3 prop     = %17.10e seconds\n", dtime2 - dtime1);
  VRB.Result(cname, fname, "kl3          = %17.10e seconds\n", dtime3 - dtime2);
//  VRB.Result(cname, fname, "k2pipi prop  = %17.10e seconds\n", dtime4 - dtime3);
//  VRB.Result(cname, fname, "k2pipi       = %17.10e seconds\n", dtime5 - dtime4);
  VRB.Result(cname, fname, "total        = %17.10e seconds\n", dtime3 - dtime0);

  //////////////////////////////////////////////////////////////////////
  // store propagators
  //lprop_e.store_all("lprop_raw_", lqpropw_arg.cg.mass, traj);

  Float dtime6 = dclock();

  VRB.Result(cname, fname, "store prop   = %17.10e seconds\n", dtime6 - dtime5);
#endif
  //sfree(out_d);
// if(ecache) ecache->dealloc();
}

int main(int argc,char *argv[])
{
  const char *fname = "main()";

  // Seed the random number generator this is used for shifting the source 
  // times of the exact propagators.
  srand48(time(NULL));
  // srand48(123456); //FIXME: fixed seed being used for testing
  // VRB.Result(cname, fname, "WARNING!!!! Using fixed random seed!!!!! FIX ME!!!!\n");
  
  setup(argc, argv);

  VRB.Result(cname, fname, "Starting main trajectory loop\n");

  int traj = meas_arg.TrajStart;
	int ntraj = (meas_arg.TrajLessThanLimit - traj)/meas_arg.TrajIncrement + 1;

  for(int conf = 0; conf < ntraj; ++conf) {
    VRB.Result(cname, fname, "Starting to work on trajectory %d\n", traj);

    // shift the exact propagators
  //  do_shift(traj);
    load_checkpoint(traj);

    Lattice& lat = LatticeFactory::Create(meas_arg.Fermion, meas_arg.Gluon);
    VRB.Result(cname, fname, "lattice created\n");
    // NOTE: there are 4 twists but 3 momenta. This is just
    // because of how code is written. The t twist is normally
    // zero.
    assert(l_twist_arg.Floats.Floats_len == 4);
    assert(s_twist_arg.Floats.Floats_len == 4);
    assert(d_mom_kpp.v.v_len == 3);

    VRB.Result(cname, fname, "asserted\n");
    const double *ltw = l_twist_arg.Floats.Floats_val;
    const double *stw = s_twist_arg.Floats.Floats_val;
    const int *dmom = d_mom_kpp.v.v_val;

    VRB.Result(cname, fname,
        "l quark twist (kl3) = %17.10e %17.10e %17.10e %17.10e\n",
        ltw[0], ltw[1], ltw[2], ltw[3]);
    VRB.Result(cname, fname,
        "s quark twist (kl3) = %17.10e %17.10e %17.10e %17.10e\n",
        stw[0], stw[1], stw[2], stw[3]);
    VRB.Result(cname, fname,
        "d quark mom  (k2pp) = %d %d %d\n",
        dmom[0], dmom[1], dmom[2]);

    run_all(lat, stw, ltw, dmom, traj);
    exit(-42); // cheating

    traj += meas_arg.TrajIncrement;
  }

#ifdef USE_QUDA
  EigenCacheListCleanup();
#endif

  VRB.Result(cname, fname, "Program ended normally.\n");
  End();
}
