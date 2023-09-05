// -*- mode:c++; c-basic-offset:4 -*-
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <util/lattice.h>
#include <util/lattice/fbfm.h>

#include <alg/array_arg.h>
#include <alg/alg_int.h>
#include <alg/eigcg_arg.h>
#include <alg/qpropw.h>

#include <util/gjp.h>
#include <util/time_cps.h>
#include <util/verbose.h>
#include <util/qcdio.h>
#include <util/error.h>

#include <string>
#include <vector>
#include <cassert>

#include "prop_container.h"
#include "eigcg.h"
#include "my_util.h"
#include "run_mres.h"
#include "run_prop.h"

static const char *cname = "";

USING_NAMESPACE_CPS 
using namespace std;

static void run_mres_za (const QPropW & qp, const QPropWArg & qp_arg,
                         const string & rdir, int traj)
{
  string mres_fn = rdir + "/mres_"
    + tostring (qp_arg.cg.mass) + "." + tostring (traj)+"."+tostring(qp_arg.t);

  string za_fn = rdir + "/za_"
    + tostring (qp_arg.cg.mass) + "." + tostring (traj)+"."+tostring(qp_arg.t);
  VRB.Result("","run_mres_za()","measure mres and za for t=%d mass=%e to %s %s\n",qp_arg.t,qp_arg.cg.mass,mres_fn.c_str (),za_fn.c_str ());
  run_mres (qp, qp_arg.t, mres_fn.c_str ());
  run_za (qp, qp_arg.cg.mass, qp_arg.t, za_fn.c_str ());
}

// Temporary hack, solve a 4D volume source to collect low modes,
// useful for AMA.
//
// Note: How many times we solve the volume source depends on how many
// low modes we want to solve. Lattice properties also apply.
//
// For 300 low modes, 1 propagator using mixed solver will be good
// (depends on EigCGArg).
//
// On 48^3 2 solves are needed for 600 low modes.
static void collect_lowmodes (Lattice & lat,
                              QPropWArg & qp_arg, CommonArg & com_prop)
{
  const char *fname = "collect_lowmodes()";
  Float timer0 = dclock ();

  double stop_rsd = qp_arg.cg.stop_rsd;
  double true_rsd = qp_arg.cg.true_rsd;

  qp_arg.cg.stop_rsd = 1e-10;
  qp_arg.cg.true_rsd = 1e-10;

  QPropW4DBoxArg vol_arg;
  for (int mu = 0; mu < 4; ++mu) {
    vol_arg.box_start[mu] = 0;
    vol_arg.box_size[mu] = GJP.Sites (mu);
    vol_arg.mom[mu] = 0;
  }

  // 2 solves for 600 low modes.
  for (unsigned i = 0; i < 2; ++i) {
    // QPropWVolSrc(lat, &qp_arg, &com_prop);
    QPropW4DBoxSrc qp_box (lat, &qp_arg, &vol_arg, &com_prop);
  }

  qp_arg.cg.stop_rsd = stop_rsd;
  qp_arg.cg.true_rsd = true_rsd;

  Float timer1 = dclock ();
  VRB.Result (cname, fname, "Finished collecting low modes: took %e seconds\n",
              timer1 - timer0);
}

void run_wall_prop (AllProp * prop,
                    LocArg & eloc,
                    Lattice & lat,
                    QPropWArg & qp_arg, const string & rdir,
                    int traj, bool do_mres, bool do_deflation)
{
  const char *fname = "run_wall_prop()";


  char *prop_dir = qp_arg.file;
  char prop_name[256];

  char buf[256];
  CommonArg com_prop;
  sprintf (buf, "../results/%s.%d", qp_arg.ensemble_label, traj);
  com_prop.set_filename (buf);
  CgArg cg_arg_tmp = qp_arg.cg;
  int i_prop=0;
  for (int bc = 1; bc < 2; ++bc) {
//        GJP.Tbc(bc == 0 ? BND_CND_PRD : BND_CND_APRD);

    Float timer0 = dclock ();
    int status;

    // exact propagators
    if (prop != NULL) {
      double stop_rsd = qp_arg.cg.stop_rsd;
      double true_rsd = qp_arg.cg.true_rsd;

      int n_prop=0;
      if(eloc.pattern_kind==ARRAY)
      n_prop=eloc.locs.locs_len;
      else 
      n_prop=eloc.n_loc;

      for (unsigned i = 0; i < n_prop; ++i) {
        if(eloc.pattern_kind==ARRAY)
        qp_arg.t = eloc.start+eloc.locs.locs_val[i];
        else 
        qp_arg.t = eloc.start + i* eloc.step;
        if ( (qp_arg.t <0 )  )  continue;
        if ( (qp_arg.t >= GJP.Sites(3)))  qp_arg.t = qp_arg.t% GJP.Sites(3) ;
        sprintf (prop_name, "%s/t%d", prop_dir, qp_arg.t);
        std::string str = "mkdir -p ";
        str = str +prop_name;  
        const char *command = str.c_str();
        status=system ( command);
        VRB.Result ("", "", "mkdir %s returned %d\n", prop_name, status);
        sprintf (prop_name, "%s/t%d/prop_m%0.6f.vol0000", prop_dir, qp_arg.t, qp_arg.cg.mass);
        FILE *fp=fopen(prop_name,"r");
        VRB.Result ("", "", "fopen %s returned %p\n", prop_name, fp);
        if(fp){ fclose(fp);  continue; }
        sprintf (prop_name, "%s/t%d/prop_m%0.6f", prop_dir, qp_arg.t, qp_arg.cg.mass);
        fp=fopen(prop_name,"r");
        VRB.Result ("", "", "fopen %s returned %p\n", prop_name, fp);
        if(fp){ fclose(fp);  continue; }
        if(i_prop >= eloc.max_src) continue;
//        sprintf (prop_name, "%s/t%d/prop_m%0.6f", prop_dir, qp_arg.t, qp_arg.cg.mass);
        qp_arg.file = prop_name;
        VRB.Result (cname, fname,
                    "Solving propagator at %d with stop_rsd = %e saving to %s\n",
                    qp_arg.t, qp_arg.cg.stop_rsd, qp_arg.file);
        if (do_deflation) {
          qp_arg.cg.Inverter = CG_LOWMODE_DEFL;
        }
        QPropWWallSrc qp_wall (lat, &qp_arg, &com_prop);
//        QPropWPointSrc qp_wall (lat, &qp_arg, &com_prop);
        if (do_mres) {
          run_mres_za (qp_wall, qp_arg,
                       rdir, traj);
//                       string ("../results") + (bc == 0 ? "EP" : "EA"), traj);
        }
        prop->add (qp_wall, qp_arg.t, bc == 0);
        i_prop++;
      }
      qp_arg.cg.stop_rsd = stop_rsd;
      qp_arg.cg.true_rsd = true_rsd;
    }
    qp_arg.cg.Inverter = cg_arg_tmp.Inverter;
    Float timer1 = dclock ();
    VRB.Result (cname, fname,
                "Total time for   propagators = %e seconds\n",
                timer1 - timer0);

  }
  qp_arg.file = prop_dir;

}

void run_mom_prop (AllProp * prop_e,
                   AllProp * prop,
                   IntArray & eloc,
                   Lattice & lat,
                   QPropWArg & qp_arg,
                   EigCGArg * eigcg_arg,
                   int traj, const int mom[3], bool do_deflation)
{
  const char *fname = "run_mom_prop()";

  // Ensure that all 4 directions have periodic boundary condition.
  // FIXME: This check is not perfect as we have no way detecting
  // how the actual gauge field data were manipulated.
  for (int mu = 0; mu < 4; ++mu) {
    if (GJP.Bc (mu) == BND_CND_APRD) {
      ERR.General (cname, fname, "Boundary condition does not match!\n");
    }
    if (mu < 3 && mom[mu]) {
      GJP.Bc (mu, BND_CND_APRD);
    }
  }

  char buf[256];
  CommonArg com_prop;
  sprintf (buf, "../results/%s.%d", qp_arg.ensemble_label, traj);
  com_prop.set_filename (buf);
  CgArg cg_arg_tmp = qp_arg.cg;
  // P+A
  for (int bc = 1; bc < 2; ++bc) {
//        GJP.Tbc(bc == 0 ? BND_CND_PRD : BND_CND_APRD);
#ifdef USE_BFM
    lat.BondCond ();
#endif
    EigCG *eig_cg = NULL;
    if (eigcg_arg) {
#ifdef USE_BFM
      eig_cg = new EigCG (eigcg_arg, Fbfm::use_mixed_solver);
      VRB.Result (cname, fname, "Collecting low modes...\n");
      collect_lowmodes (lat, qp_arg, com_prop);

      const string fn = string ("../results") + (bc == 0 ? "EP" : "EA")
        + "/eigH_mom_" + tostring (qp_arg.cg.mass) + "." + tostring (traj);

      eig_cg->printH (fn);
#else
      ERR.General (cname, fname, "EigCG solver not implemented.\n");
#endif
    }
    // exact propagators
    if (prop_e != NULL) {
      double stop_rsd = qp_arg.cg.stop_rsd;
      double true_rsd = qp_arg.cg.true_rsd;

      qp_arg.cg.stop_rsd = 1e-8;
      qp_arg.cg.true_rsd = 1e-8;
      for (unsigned i = 0; i < eloc.v.v_len; ++i) {
        qp_arg.t = eloc.v.v_val[i];
        VRB.Result (cname, fname, "Solving exact propagator at %d\n", qp_arg.t);
        if (do_deflation) {
          qp_arg.cg.Inverter = CG_LOWMODE_DEFL;
        }
        QPropWMomCosTwistSrc qp_mom (lat, &qp_arg, mom, &com_prop);
        prop_e->add (qp_mom, qp_arg.t, bc == 0);
      }
      qp_arg.cg.stop_rsd = stop_rsd;
      qp_arg.cg.true_rsd = true_rsd;
    }
    qp_arg.cg.Inverter = cg_arg_tmp.Inverter;
    // inexact propagators
    for (int t = 0; t < GJP.Sites (3); ++t) {
      qp_arg.t = t;
      if (do_deflation) {
        qp_arg.cg.Inverter = CG_LOWMODE_DEFL;
      }
      QPropWMomCosTwistSrc qp_mom (lat, &qp_arg, mom, &com_prop);
      prop->add (qp_mom, qp_arg.t, bc == 0);
    }
    qp_arg.cg.Inverter = cg_arg_tmp.Inverter;
    delete eig_cg;
#ifdef USE_BFM
    lat.BondCond ();
#endif
  }

  // Note: If I call lat.BondCond() even times, then there is no
  // overall effect.
  for (int mu = 0; mu < 4; ++mu) {
    GJP.Bc (mu, BND_CND_PRD);
  }
}

void run_box_prop (AllProp * prop,
                   Lattice & lat,
                   QPropWArg & qp_arg,
                   QPropW4DBoxArg & box_arg, int traj, bool do_deflation)
{
  const char *fname = "run_box_prop()";

  // Check boundary condition. We need this to ensure that we are
  // doing P + A and P - A, not A + P and A - P (I think it's OK to
  // skip this check, though).
//   if(GJP.Tbc() == BND_CND_APRD) {
//ERR.General(cname, fname, "Boundary condition does not match!\n");
//}

  char buf[256];
  CommonArg com_prop;
  sprintf (buf, "../results/%s.%d", qp_arg.ensemble_label, traj);
  com_prop.set_filename (buf);
  CgArg cg_arg_tmp = qp_arg.cg;
  // A only
  for (int bc = 1; bc < 2; ++bc) {
    //GJP.Tbc(bc == 0 ? BND_CND_PRD : BND_CND_APRD);
#ifdef USE_BFM
    lat.BondCond ();
#endif
    // inexact propagators
    for (int t = 0; t < GJP.Sites (3); ++t) {
      box_arg.box_start[3] = qp_arg.t = t;
      if (do_deflation) {
        qp_arg.cg.Inverter = CG_LOWMODE_DEFL;
      }
      QPropWZ3BWallSrc qp_z3 (lat, &qp_arg, &box_arg, &com_prop);

//            run_mres_za(qp_z3, qp_arg,
//                        string("../results") + (bc == 0 ? "P" : "A"),
//                        traj);

      prop->add (qp_z3, box_arg.box_start[3], bc == 0);
    }
    qp_arg.cg.Inverter = cg_arg_tmp.Inverter;
#ifdef USE_BFM
    lat.BondCond ();
#endif
  }

  // Note: If I call lat.BondCond() even times, then there is no
  // overall effect.
//GJP.Tbc(BND_CND_PRD);
}
