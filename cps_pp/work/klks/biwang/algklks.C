//my class
#include "algklks.h"

//cps classes
#include <config.h>
#include <alg/alg_meas.h>   
#include <alg/qpropw.h>
#include <alg/common_arg.h>
#include <alg/no_arg.h>
#include <alg/qpropw_arg.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <util/qcdio.h>
#include <alg/fix_gauge_arg.h>
#include <alg/alg_fix_gauge.h>
#include <alg/alg_rnd_gauge.h>

using namespace cps;
using namespace std;

AlgKlKs::AlgKlKs(Lattice& _lat, 
		CommonArg* _common_arg, 
		QPropWArg* _lqpropw_arg, 
		QPropWArg* _sqpropw_arg, 
		QPropWArg* _cqpropw_arg,
		int _minsep, 
		int _nhits,int _sep, bool _do_exact) : 
	Alg(_lat, _common_arg), 
	cname("AlgKlKs()"), 
	lat(&_lat),
	minsep(_minsep), 
	nhits(_nhits),
	sep(_sep),
	do_exact(_do_exact),
	glb_t(GJP.TnodeSites()*GJP.Tnodes()),
	toffset(GJP.TnodeSites()*GJP.TnodeCoor()),
	size_3d(GJP.VolNodeSites()/GJP.TnodeSites()),
	size_4d(GJP.VolNodeSites()),
	size_3d_glb(GJP.VolSites()/glb_t)
{
	const char *fname = "AlgKlKs()";

	////	initialize arguments
	lqpropw_arg = _lqpropw_arg;
	sqpropw_arg = _sqpropw_arg;
	cqpropw_arg = _cqpropw_arg;
	cmass_1=0.3;
	cmass_2=0.38;

	min_3pt_sep = 12;
	sep_step = 4;
	max_3pt_sep = 28;

	//initialize propagator
	lwall = new QPropWS*[glb_t];
	lpoint = new QPropWS*[glb_t];
	swall = new QPropWS*[glb_t];
	lmcpoint = new QPropWS*[glb_t];
	for(int i = 0; i < glb_t; i++)
	{
		lwall[i] = NULL;
		swall[i] = NULL;
		lpoint[i]=NULL;
		lmcpoint[i] = NULL;
	}
	lmcloop = new QPropWS*[nhits];
	lloop = new QPropWS*[nhits];
	for(int i = 0; i < nhits; i++){
		lmcloop[i] = NULL;
		lloop[i] = NULL;
	}
	lmcloop_avg= NULL;
	l_low = NULL;
}

AlgKlKs::~AlgKlKs()
{
	char *fname = "~AlgKlKs()";

	VRB.Result(cname,fname,"destructor called. \n");
	int glb_t = GJP.TnodeSites()*GJP.Tnodes();
	for(int i = 0; i < glb_t; i++)
	{
		delete lwall[i];
		delete swall[i];
		if(lpoint[i] != NULL)
			delete lpoint[i];
		if (lmcpoint[i] != NULL)
			delete lmcpoint[i];
		delete [] lwwprop[i];
		delete [] swwprop[i];
	}
	for(int i = 0; i < nhits; i++){
		delete lloop[i];
		delete lmcloop[i];
	}

	delete [] lwall;
	delete [] lpoint;
	delete [] swall;
	delete [] lmcpoint;
	delete [] lmcloop;
	delete [] lloop;
	delete lmcloop_avg;
	delete [] lwwprop;
	delete [] swwprop;	
	delete l_low;
}	

void AlgKlKs::Run1()
{
	const char* fname="Run()";
	bool do_type12 = do_exact;
	char label[100];
	char ama_label[100];

	if (do_exact == false) {
		lqpropw_arg->cg.stop_rsd = 1.0e-4;
		sqpropw_arg->cg.stop_rsd = 1.0e-4;
		cqpropw_arg->cg.stop_rsd = 1.0e-4;
		sprintf(ama_label, "1e-4");
	} else {
		lqpropw_arg->cg.stop_rsd = 1.0e-8;
		sqpropw_arg->cg.stop_rsd = 1.0e-8;
		cqpropw_arg->cg.stop_rsd = 1.0e-8;
		sprintf(ama_label, "1e-8");
	}

	

	int Ls =  GJP.SnodeSites();
	extern double mobius_fac;

	LightWallProp();
	LightPointProp();
//	set_arg_alpha(Ls,sqpropw_arg->cg.mass,sqpropw_arg->cg.mass,mobius_fac);
	StrangeWallProp();
 	WallSrcWallSnkProp();

	RunPion(ama_label);
	RunKaon(ama_label);
	Runpipi(ama_label);

	boxsize = 16;
	T = boxsize + 1;
	cqpropw_arg->cg.mass = cmass_1;
//	set_arg_alpha(Ls,cqpropw_arg->cg.mass,cqpropw_arg->cg.mass,mobius_fac);
	
	////LightCharmRandVolProp(cmass_1);
	LightRandVolProp_A2A();
	syncNode();

	//lat -> free_lanczos();

	CharmRandVolProp(cmass_1);
	AvgSelfLoop();
	sprintf(label, "%4.2f_%s", cmass_1,ama_label);
	RunThreePointCorr(label);
	Runk2pipi(label);
	RunType34(label);
	CharmPointProp(cmass_1);
	RunType12(label);	
//	CharmRandVolProp(cmass_2);
//	AvgSelfLoop();
//	sprintf(label, "%4.2f_%s", cmass_2,ama_label);
//	RunThreePointCorr(label);
//	Runk2pipi(label);
//	RugtnType34(label);

}

void AlgKlKs::Run2(){
	const char * fname="Run Point ()";
	char label[100];

	LightPointProp();

	CharmPointProp(cmass_1);
	sprintf(label, "%4.2f", cmass_1);
	RunType12(label);

	CharmPointProp(cmass_2);
	sprintf(label, "%4.2f", cmass_2);
	RunType12(label);

	CharmPointProp(cmass_3);
	sprintf(label, "%4.2f", cmass_3);
	RunType12(label);

	CharmPointProp(cmass_4);
	sprintf(label, "%4.2f", cmass_4);
	RunType12(label);
}



void AlgKlKs::SaveCorr(Rcomplex* corr, char* filename, int len) const
{
	const char* fname = "Save()";
	FILE *fp;
	if((fp=Fopen(filename,"w"))==NULL)
		ERR.General(cname, fname, "Cann't open %s\n", filename);
	else {
		for(int i = 0; i < len; i++)
			Fprintf(fp, "%17.10e\n", corr[i].real()); 
	}
	Fclose(fp);
}

void AlgKlKs::WallSrcWallSnkProp()
{
	const char* fname = "WallSrcWallSnkProp()";
	double cost = -dclock();

	//initialize wall src wall snk propagators
	lwwprop = new WilsonMatrix*[glb_t];
	swwprop = new WilsonMatrix*[glb_t];
	for(int i = 0; i < glb_t; i++)
	{
		lwwprop[i] = new WilsonMatrix[glb_t];
		swwprop[i] = new WilsonMatrix[glb_t];
	}

	for(int tsrc = 0; tsrc < glb_t; tsrc++)
		for(int tsnk = 0; tsnk < glb_t; tsnk++)
		{
			lwwprop[tsrc][tsnk] = lwall[tsrc]->WallSinkProp(tsnk);
			swwprop[tsrc][tsnk] = swall[tsrc]->WallSinkProp(tsnk);
		}
	
	cost += dclock();
	print_time(cname, fname, cost);
}

void AlgKlKs::AvgSelfLoop()
{

	QPropW loop(*lat, lqpropw_arg, common_arg);
	loop.Allocate(0);

#pragma omp parallel for
	for(int site = 0; site < GJP.VolNodeSites(); site++)
	{
		loop[site] = 0.0;
		for(int k = 0; k < nhits; k++)
		{
			loop[site] += (*lmcloop[k])[site] * (1.0/nhits);
		}
	}
	if (lmcloop_avg != NULL ) delete lmcloop_avg;
	lmcloop_avg = new QPropWS(loop, 0);

//	WilsonMatrix lmcloop_sum = 0;
//	for(int site = 0; site < GJP.VolNodeSites(); site++){
//		lmcloop_sum += (*lmcloop_avg)[site];
//	}
//	QMP_sum_double_array((double*)&lmcloop_sum, 144*2);
	

//	if (UniqueID() == 0) {
//		for (int j = 0; j < 144; j++) {
//			printf("%e,%e\n", ((double*)&lmcloop_sum)[2*j+0], ((double*)&lmcloop_sum)[2*j+1]);
//		}
//	}

	syncNode();
}

