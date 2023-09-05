/*
 * put your initializations here
 */

#include <alg/threept_arg.h>
CPS_START_NAMESPACE

ThreePtArg::ThreePtArg()
{
  t_src=0;
  t_shift=0;
  num_hits=0;
  do_susy=1;
  do_zero_mom=1;
  do_first_mom=1;
  do_second_mom=1;
  do_third_mom=1;
  do_pipi_non_zero_tot_mom=1;
  do_p_plus_a_kaon=1;
  do_kaon_at_walls=1;
  do_kaons_tK=1;
  chkpoints=1;
};

CPS_END_NAMESPACE
