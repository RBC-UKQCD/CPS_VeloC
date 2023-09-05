// -*- mode:c++; c-basic-offset:4 -*-
#ifndef INCLUDED_RUN_PROP_H_KL3
#define INCLUDED_RUN_PROP_H_KL3

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <string>
#include <vector>
#include <cassert>

//#include <util/lattice.h>
#include <alg/qpropw_arg.h>
#include <alg/array_arg.h>
#include <alg/eigcg_arg.h>
#include <alg/loc_arg.h>

#include "prop_container.h"

namespace cps {
    class Lattice;
};

void run_wall_prop(AllProp *prop,
                   cps::LocArg &eloc,
                   cps::Lattice &lat,
                   cps::QPropWArg &qp_arg,
                   const std::string & rdir,
                   int traj,
                   bool do_mres,
                   bool do_deflation);

void run_mom_prop(AllProp *prop_e,
                  AllProp *prop,
                  cps::LocArg &eloc,
                  cps::Lattice &lat,
                  cps::QPropWArg &qp_arg,
                  cps::EigCGArg *eigcg_arg,
                  int traj,
                  const int mom[3],
		  bool do_deflation);

void run_box_prop(AllProp *prop,
                  cps::Lattice &lat,
                  cps::QPropWArg &qp_arg,
                  cps::QPropW4DBoxArg &box_arg,
                  int traj,
		  bool do_deflation);

#endif
