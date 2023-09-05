#ifndef ALGKLKS_H
#define ALGKLKS_H

#include <config.h>
#include <util/lattice/fbfm.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/qpropw_arg.h>
#include <alg/qpropw.h>
#include <alg/qpropw_s.h>
#include <alg/alg_fix_gauge.h>
#include <util/lattice.h>
#include <util/rcomplex.h>
#include <util/qcdio.h>
#include <util/time_cps.h>
#include <alg/lanc_arg.h>
#include <omp.h>

#include <util/timer.h>
using namespace cps;
//void eq_mult_simd_2( complex<double> * xmat,
//		    const complex<double>* amat,
//				    const complex<double>* bmat );


static int Mee = 0;
/*
inline double set_arg_alpha(int Ls, double mass, double map_mass, double alpha)
{
	  Fbfm::arg_map[map_mass].ScaledShamirCayleyTanh(mass, 1.8, Ls, alpha);
		Fbfm::arg_map[map_mass].CGdiagonalMee = Mee;
}
*/
class AlgKlKs : public Alg
{
private:
	const char *cname;	
	const int minsep; // mininum time separation betweem kaon wall source and interaction box 
	int nhits; // eaxct hits of random wall sources
	const int sep; //time difference between 2 pions
	int glb_t; // time dimension of global lattice
	const int toffset;// time offset of local node
	const int size_3d; // 3d lcoal volume
	const int size_4d; // 4d local volume
	const int size_3d_glb; // 3d global volume
	int boxsize; // tb-ta 
	bool do_exact;
	int T; // tb-ta+1
	double cmass_1;
	double cmass_2;
	double cmass_3;
	double cmass_4;
	int min_3pt_sep;
	int sep_step;
	int max_3pt_sep;

	Lattice* lat;

	//Argument for lanczos
	LancArg lanc_arg;

	//Argument for propagators
	QPropWArg *lqpropw_arg;
	QPropWArg *sqpropw_arg;
	QPropWArg *cqpropw_arg;

	//wall soruce propagtors
	QPropWS **lwall;
	QPropWS **swall;
	
	//light minus charm loops
	QPropWS **lmcloop;
	QPropWS **lloop;
	QPropWS *lmcloop_avg;
	
	//point source propagators
	QPropWS ** lpoint;
	QPropWS *l_low;
	QPropWS **lmcpoint;

	//wall src wall snk propagators
	WilsonMatrix **lwwprop;
	WilsonMatrix **swwprop;

	//wall and random wall source propagators
	void LightWallProp();
	void StrangeWallProp();
	void WallSrcWallSnkProp();

	void LightCharmRandVolProp(double cmass_1);
	void LightCharmRandVolProp_A2A(double cmass_1);
	void LightRandVolProp_A2A();
	void CharmRandVolProp(double cmass);
	void AvgSelfLoop();

	//point source propagtors
	void LightPointProp();
	void CharmPointProp(double cmass);
  
	//two point functions
	void RunPion(const char* label) const;
  void RunKaon(const char* label) const;
	void Runpipi(const char* label) const;
  
	//three point functions
	void RunThreePointCorr(const char* label) const;
	void Runk2pipi(const char* label) const;
	
	//four point functions
	void RunType12(const char* label)const;
	void RunType34(const char* label);
	
	void syncNode() const {
		double x(0.0);
		QMP_sum_double_array(&x, 1);
	}

//	WilsonMatrix& eq_mult(WilsonMatrix& xmat, WilsonMatrix& amat,
//			WilsonMatrix& bmat) const {
//		eq_mult_simd_2(xmat.ptr(), amat.ptr(), bmat.ptr());
//		return xmat;
//	}

	//half of type 3 and type 4 contractions
	WilsonMatrix **rand_uppc, **rand_uppv;
	WilsonMatrix **rand_downpc, **rand_downpv;
	Rcomplex **rand_leftpc, **rand_leftpv;
	Rcomplex **rand_rightpc, **rand_rightpv;
	
	void SaveCorr(Rcomplex* corr, char* filename, int len) const;

public:
	AlgKlKs(Lattice & lat, 
			CommonArg *common_arg, 
			QPropWArg *lqpropw_arg, 
			QPropWArg *sqpropw_arg, 
			QPropWArg *cqpropw_arg, 
			int minsep, 
			int nhits, int sep, bool do_exact);
	void Run1();
	void testProp();
	void Run2();

	~AlgKlKs();
};

#endif
