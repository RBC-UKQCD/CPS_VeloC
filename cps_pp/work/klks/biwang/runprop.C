#include "algklks.h"
//#include "lanc_io.h"


//extern RandomNumberGenerator  rng;
//light quark wall source propagators
void AlgKlKs::LightWallProp()
{
	//Lanczos_5d<lanc_type>* eig = (Lanczos_5d<lanc_type>*)(lat->get_eig());
//	BFM_Krylov::Lanczos_5d<lanc_type>* eig = (BFM_Krylov::Lanczos_5d<lanc_type>* )(lat->get_eig());
//	Fbfm *lat_bfm = dynamic_cast<Fbfm *>(lat);
//	lat_bfm->set_deflation(&(eig->bq), &(eig->bl), eig->get);

	char fname [500];
	double total_cost=-dclock();

	CommonArg carg("prop", "");
  char filename[1024];
  sprintf(filename, "%s_prop.txt", common_arg->filename);
  carg.set_filename(filename);
	
	//wall source propagators
	//int glb_t = GJP.TnodeSites()*GJP.Tnodes();
	for(int i = 0; i < glb_t; i++)
	{
		double cost = -dclock();
		sprintf(fname,"LightWallProp(),t=%d",i);
		lqpropw_arg->t = i;
		QPropWWallSrc prop(*lat, lqpropw_arg, &carg);
		if(lwall[i] != NULL) delete lwall[i];
		lwall[i] = new QPropWS(prop, 0);
		cost += dclock();
		print_time(cname, fname, cost);
		syncNode();
	}
	total_cost+=dclock();
	sprintf(fname,"Light Wall Prop() Total");
	print_time(cname,fname,total_cost);
	if (UniqueID() == 0)
		printf("After Light Wall Prop\n");
//	lat_bfm->unset_deflation();
	syncNode();

}

//strange wall source propagators
void AlgKlKs::StrangeWallProp()
{
	char fname[500];
	double total_cost=-dclock();
	CommonArg carg("prop", "");
	char filename[1024];
	sprintf(filename, "%s_prop.txt", common_arg->filename);
	carg.set_filename(filename);

//	Lattice &lat = AlgLattice();

	//wall source propagators
	//int glb_t = GJP.TnodeSites()*GJP.Tnodes();
	for(int i = 0; i < glb_t; i++)
	{
	double cost = -dclock();
	sprintf(fname,"StrangeWallProp(),t=%d",i);
	sqpropw_arg->t = i;
	QPropWWallSrc prop(*lat, sqpropw_arg, &carg);
	if(swall[i] != NULL) delete swall[i];
	swall[i] = new QPropWS(prop, 0);
	cost += dclock();
	print_time(cname, fname, cost);
	}
	total_cost+=dclock();
	sprintf(fname,"Strange Wall Prop() Total");
	print_time(cname,fname,total_cost);

}

//light and charm quark random Volume source propagators
void AlgKlKs::LightCharmRandVolProp(double cmass1)
{
	char fname[500]; 
	double total_cost=-dclock();
	CommonArg carg("prop", "");
	char filename[1024];
	sprintf(filename, "%s_prop.txt", common_arg->filename);
	carg.set_filename(filename);

//	Lattice & lat = AlgLattice();

	QPropWRandArg rand_arg;
	QPropWRandVolSrc *lrqprop;
	QPropWRandVolSrc *crqprop;

	//use ZTWO random src
	rand_arg.rng = TEST;

	//high mode part
	for(int hit = 0; hit < nhits; hit++)
	{
		double cost = -dclock();
		lrqprop = new QPropWRandVolSrc(*lat, lqpropw_arg, &rand_arg, &carg);
		cost+=dclock();
		sprintf(fname,"Light Random Volume Prop(),hit=%d",hit);
		print_time(cname,fname,cost);

		QPropWRandVolSrc selfloop(*lrqprop);

		(cqpropw_arg->cg).mass=cmass1;
		crqprop = new QPropWRandVolSrc(*lat, cqpropw_arg, lrqprop, &carg);

		for(int ind = 0; ind < size_4d; ind++)
		{
			Site s(ind);
			WilsonMatrix loop((*lrqprop)[ind]);
			loop -= (*crqprop)[ind];
			loop *= conj(lrqprop->rand_src(ind));
			selfloop[ind] = loop;
		}

		delete lrqprop;
		delete crqprop;
		if(lmcloop[hit] != NULL) delete lmcloop[hit];
		lmcloop[hit] = new QPropWS(selfloop, 1);
	}

	sprintf(fname,"Light  charm Rondom Vol Prop(), total");
	total_cost += dclock();
	print_time(cname, fname, total_cost);

	return;
}




//void AlgKlKs::LightCharmRandVolProp_A2A(double cmass1)
//{
//	char fname[500]; 
//	double total_cost=-dclock();
//	CommonArg carg("prop", "");
//	char filename[1024];
//	sprintf(filename, "%s_prop.txt", common_arg->filename);
//	carg.set_filename(filename);
//	syncNode();
//	VRB.Result(cname,"LightRandVol A2A","begin low mode propagator. \n");
//
//	//Lattice & lat = AlgLattice();
//
//	QPropWRandArg rand_arg;
//	QPropWRandVolSrc *lrqprop;
//	QPropWRandVolSrc *crqprop;
//	QPropWA2ALoopLow loop_low(*lat,lqpropw_arg, &carg);
//	syncNode();
//	VRB.Result(cname,"LightRandVol A2A","After finding low mode propagator. \n");
//
//	//use ZTWO random src
//	rand_arg.rng = ZTWO;
//	InverterType tmpInverter = lqpropw_arg->cg.Inverter;
//	lqpropw_arg->cg.Inverter = DEFLATED_MIXCG_H;
//
//	//high mode part
//	for(int hit = 0; hit < nhits; hit++)
//	{
//		double cost = -dclock();
//		QPropW selfloop(*lat, &carg);
//		selfloop.Allocate(0);
//		lrqprop = new QPropWRandVolSrc(*lat, lqpropw_arg, &rand_arg, &carg);
//		cost+=dclock();
//		sprintf(fname,"Light Random Volume Prop(),hit=%d",hit);
//		print_time(cname,fname,cost);
//		(cqpropw_arg->cg).mass=cmass1;
//		crqprop = new QPropWRandVolSrc(*lat, cqpropw_arg, lrqprop, &carg);
//
//		for(int ind = 0; ind < size_4d; ind++)
//		{
//			Site s(ind);
//			WilsonMatrix loop((*lrqprop)[ind]);
//			loop -= (*crqprop)[ind];
//			loop *= conj(lrqprop->rand_src(ind));
//			loop += loop_low[ind];
//			selfloop[ind] = loop;
//		}
//
//		delete lrqprop;
//		delete crqprop;
//		if(lmcloop[hit] != NULL) delete lmcloop[hit];
//		lmcloop[hit] = new QPropWS(selfloop, 0);
//	}
//
//	sprintf(fname,"Light  charm Rondom Vol Prop A2A, total");
//	total_cost += dclock();
//	print_time(cname, fname, total_cost);
//	lqpropw_arg->cg.Inverter = tmpInverter;
//
//	return;
//}

void AlgKlKs::LightRandVolProp_A2A()
{
//	BFM_Krylov::Lanczos_5d<lanc_type>* eig = (BFM_Krylov::Lanczos_5d<lanc_type>*)(lat->get_eig());
//	Fbfm *lat_bfm = dynamic_cast<Fbfm *>(lat);
//	lat_bfm->set_deflation(&(eig->bq), &(eig->bl), eig->get);

	char fname[500]; 
	double total_cost=-dclock();
	CommonArg carg("prop", "");
	char filename[1024];
	sprintf(filename, "%s_prop.txt", common_arg->filename);
	carg.set_filename(filename);
	syncNode();
	VRB.Result(cname,"LightRandVol A2A","begin low mode propagator. \n");

	QPropWRandArg rand_arg;
	QPropWRandVolSrc *lrqprop;
	//QPropW loop_low(*lat, &carg);
	//loop_low.Allocate(0);
	/*
	QPropWA2ALoopLow loop_low(*lat,lqpropw_arg, &carg);
	l_low = new QPropWS(loop_low, 0);
	syncNode();
	VRB.Result(cname,"LightRandVol A2A","After finding low mode propagator. \n");

	if (UniqueID() == 0) {
			printf("printing loop_low\n");
		for (int site = 0; site < 10; site++)
			for (int j = 0; j < 144; j++)
			printf("%e,%e\n", ((double*)&loop_low[site])[2*j+0], ((double*)&loop_low[site])[2*j+1]);
	}
	*/	//Comment this part to nullify the low mode, By Bigeng
	//if (UniqueID() == 0) {
	//		printf("printing l_low\n");
	//	for (int site = 0; site < 10; site++){
	//		WilsonMatrix tmp = 	(*l_low)[site];
	//		for (int j = 0; j < 144; j++){
	//		printf("%e,%e\n", ((double*)&tmp)[2*j+0], ((double*)&tmp)[2*j+1]);
	//		}
	//	}
	//}

	////use ZTWO random src
	rand_arg.rng = TEST;

	
	
//	lat_bfm->set_high_mode_only(true);

	//high mode part
	for(int hit = 0; hit < nhits; hit++)
	{
		double cost = -dclock();
		lrqprop = new QPropWRandVolSrc(*lat, lqpropw_arg, &rand_arg, &carg);
		cost+=dclock();
		sprintf(fname,"Light Random Volume Prop(),hit=%d",hit);
		print_time(cname,fname,cost);

		if(lloop[hit] != NULL) delete lloop[hit];
		lloop[hit] = new QPropWS(*lrqprop, 1);
		delete lrqprop;
	}
	
	total_cost += dclock();
//	lat_bfm->set_high_mode_only(false);
	
	
	if (UniqueID() == 0)
		cout << "after unset high_mode_only" << endl;
//	lat_bfm->unset_deflation();
	syncNode();
	print_time(cname, fname, total_cost);
	return;
}

void AlgKlKs::CharmRandVolProp(double cmass1)
{
	char fname[500]; 
	double total_cost=-dclock();
	CommonArg carg("prop", "");
	char filename[1024];
	sprintf(filename, "%s_prop.txt", common_arg->filename);
	carg.set_filename(filename);

//	Lattice & lat = AlgLattice();

	QPropWRandArg rand_arg;
	QPropWRandVolSrc *crqprop;

	(cqpropw_arg->cg).mass=cmass1;

	//high mode part
	for(int hit = 0; hit < nhits; hit++)
	{
		double cost = -dclock();
		QPropW selfloop(*lat, &carg);
		selfloop.Allocate(0);
		sprintf(fname,"charm Random Volume Prop(),hit=%d",hit);
		crqprop = new QPropWRandVolSrc(*lat, cqpropw_arg, lloop[hit], &carg);
		cost+=dclock();
		print_time(cname,fname,cost);

		for(int ind = 0; ind < size_4d; ind++)
		{
			Site s(ind);
			//WilsonMatrix loop(0);
			WilsonMatrix loop((*lloop[hit])[ind]);
			loop -= (*crqprop)[ind];
			loop *= conj(crqprop->rand_src(ind));
			//loop += (*l_low)[ind];	//Comment this line to nullify the low mode, By Bigeng
			selfloop[ind] = loop;
		}

		delete crqprop;
		if(lmcloop[hit] != NULL) delete lmcloop[hit];
		lmcloop[hit] = new QPropWS(selfloop, 0);
	}
	
	sprintf(fname,"charm Rondom Vol Prop , total");
	total_cost += dclock();
	print_time(cname, fname, total_cost);
	return;
}


//light and charm quark point source propagators
void AlgKlKs::CharmPointProp(double cmass)
{
	char fname [500];
	double total_cost=-dclock();
	CommonArg carg("prop", "");
	char filename[1024];
	sprintf(filename, "%s_prop.txt", common_arg->filename);
	carg.set_filename(filename);
//	Lattice &lat = AlgLattice();

	(cqpropw_arg->cg).mass=cmass;

	QPropW loop(*lat, &carg);
	loop.Allocate(0);

	//wall source propagators
//	int glb_t = GJP.TnodeSites()*GJP.Tnodes();
	for(int i = 0; i < glb_t; i++)
	{
		double cost = -dclock();
		cqpropw_arg->x = (4*i)%GJP.Sites(0);
		cqpropw_arg->y = (4*i)%GJP.Sites(1);
		cqpropw_arg->z = (4*i)%GJP.Sites(2);
		cqpropw_arg->t = i;
		QPropWPointSrc cprop(*lat, cqpropw_arg, &carg);

		cost+=dclock();
		sprintf(fname,"charm point Prop(),t=%d",i);
		print_time(cname,fname,cost);

#pragma omp parallel for
		for(int ind = 0; ind < GJP.VolNodeSites(); ind++)
		{
			loop[ind]=(*lpoint[i])[ind]-cprop[ind];
		}
		if(lmcpoint[i] != NULL) delete lmcpoint[i];
		lmcpoint[i] = new QPropWS(loop, 0);
	}
	sprintf(fname,"charm point Prop(), Total");
	total_cost += dclock();
	print_time(cname, fname, total_cost);
}

void AlgKlKs::LightPointProp()
{
	char fname [500];
	double total_cost=-dclock();
	CommonArg carg("prop", "");
	char filename[1024];
	sprintf(filename, "%s_prop.txt", common_arg->filename);
	carg.set_filename(filename);

//	Lattice &lat = AlgLattice();

	//wall source propagators
///	int glb_t = GJP.TnodeSites()*GJP.Tnodes();
	for(int i = 0; i < glb_t; i++)
	{
		double cost = -dclock();
		lqpropw_arg->x = (4*i)%GJP.Sites(0);
		lqpropw_arg->y = (4*i)%GJP.Sites(1);
		lqpropw_arg->z = (4*i)%GJP.Sites(2);
		lqpropw_arg->t = i;
		QPropWPointSrc lprop(*lat, lqpropw_arg, &carg);
		cost+=dclock();
		if(lpoint[i] != NULL) delete lpoint[i];
		lpoint[i] = new QPropWS(lprop, 0);
		sprintf(fname,"Light Point Prop(), t=%d", i);
		print_time(cname,fname,cost);
	}
	sprintf(fname,"light point Prop(), Total");
	total_cost += dclock();
	print_time(cname, fname, total_cost);
}
