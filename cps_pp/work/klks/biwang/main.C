//cps class
//#include <util/lattice/fbfm.h>
#include <alg/common_arg.h> 
#include <alg/do_arg.h>     
#include <alg/meas_arg.h>   
#include <alg/cg_arg.h>     
#include <alg/qpropw_arg.h>   
#include <alg/qpropw.h>   
#include <alg/alg_rnd_gauge.h>
#include <alg/eigcg_arg.h>
#include <alg/alg_fix_gauge.h>
#include <util/qioarg.h>    
#include <alg/alg_meas.h>   
#include <config.h>		   
#include <util/lattice.h>   
#include <util/gjp.h>       
#include <util/time_cps.h>       
#include <util/verbose.h>   
#include <util/error.h>     
#include <util/qcdio.h>     
#include <util/qioarg.h>    
#include <util/ReadLatticePar.h> 

//c++ classes
#include <ctime>
#include <iostream>
#include <vector>
#include <sstream>
#include <unistd.h>
#include <sys/stat.h>
#include "algklks.h"
//#include "lanc_io.h"	//Uncommented by Bigeng


//includes origin from Chulwoo, copied from Hackathon_2018 code
 #include <omp.h>
 #include <config.h>
 #include <math.h>
 #include <util/dirac_op.h>
 #include <util/wilson.h>
 #include <comms/scu.h>
 #include <alg/alg_hmd.h>
 #include <util/omp_wrapper.h>
 #include <util/command_line.h>


//namespace include
using namespace cps;

#define decode_vml(arg_name)  {                                       \
	if ( ! arg_name.Decode(#arg_name".vml", #arg_name) )          \
	ERR.General(cname, fname, "Bad " #arg_name ".vml.\n");    \
} 

void ReadGaugeField(const MeasArg &meas_arg);
void ReadRNG(const MeasArg &meas_arg);


inline int Chdir(const char* dir)
{
  const char* fname = "Chdir(char*)";

  if(chdir(dir) != 0){
    ERR.General("", fname, "Changing to directory %s failed.\n", dir);
  }

  return 0;
}

/*
void init_bfm(int *argc, char **argv[])
{
    QDP::QDP_initialize(argc, argv);
    multi1d<int> nrow(Nd);

    for(int i = 0; i< Nd; ++i)
        nrow[i] = GJP.Sites(i);

    Layout::setLattSize(nrow);
    Layout::create();

    bfmarg::Threads(64);
    bfmarg::Reproduce(0);
    bfmarg::ReproduceChecksum(0);
    bfmarg::ReproduceMasterCheck(0);
    bfmarg::Verbose(0);

    //Fbfm::use_mixed_solver = true;	//Comment this line from Ziyuan
}
*/

// double mobius_fac=2.66667; //Ziyuan's value
double mobius_fac=2.0000;
int main(int argc,char *argv[])
{
	Start(&argc, &argv);

	const char *cname=argv[0];
	const char *fname="main(int,char**)";
	//--Argument files initialize--------------
	chdir(argv[1]);
	CommonArg common_arg("klks","");
	DoArg do_arg;
	DoArgExt doext_arg;	//Added by BIgeng, to accomodate the new features on GPU
	MeasArg meas_arg;
	QPropWArg cqpropw_arg;
	QPropWArg lqpropw_arg;
	QPropWArg sqpropw_arg;
	//LancArg lanc_arg;
	//extern QudaArg QudaParam;//Added by Bigeng, for GPU running. Already defined in the library.	


	// Added for arguement IO and IO
	char *out_file = NULL;
	CommandLine::is(argc,argv);
	printf("Current folder is: %s \n", argv[1]);
	//out_file = CommandLine::arg();
	//Chdir(CommandLine::arg()); //Added by Bigeng
	GJP.ZMobius_PC_Type (ZMOB_PC_ORIG); ////Added by Bigeng, for GPU running
	decode_vml(doext_arg); //vml file copied from Hackathon_2018

	
	FixGaugeArg fix_gauge_arg;
	lqpropw_arg.y=0;
	lqpropw_arg.z=0;
	lqpropw_arg.t=0; //Will also need source times tK[i] for the kaon correlator
	lqpropw_arg.gauge_fix_src=0;
	lqpropw_arg.gauge_fix_snk=0;
	lqpropw_arg.store_midprop=0;
	lqpropw_arg.do_half_fermion=0; 
	lqpropw_arg.cg.max_num_iter=99999;
	lqpropw_arg.cg.stop_rsd=1e-8;
	lqpropw_arg.cg.true_rsd=1e-8;
	lqpropw_arg.cg.Inverter=CG;
	lqpropw_arg.mob_arg_s=0;
	//lqpropw_arg.cg.mass=0.001;
	lqpropw_arg.cg.mass=0.01;
	cqpropw_arg.y=0;
	cqpropw_arg.z=0;
	cqpropw_arg.t=0; //Will also need source times tK[i] for the kaon correlator
	cqpropw_arg.gauge_fix_src=0;
	cqpropw_arg.gauge_fix_snk=0;
	cqpropw_arg.store_midprop=0;
	cqpropw_arg.do_half_fermion=0; 
	cqpropw_arg.cg.max_num_iter=99999;
	cqpropw_arg.cg.stop_rsd=1e-8;
	cqpropw_arg.cg.true_rsd=1e-8;
	cqpropw_arg.cg.Inverter=CG;
	cqpropw_arg.mob_arg_s=0;
	cqpropw_arg.cg.mass=0.3;

	sqpropw_arg.y=0;
	sqpropw_arg.z=0;
	sqpropw_arg.t=0; //Will also need source times tK[i] for the kaon correlator
	sqpropw_arg.gauge_fix_src=0;
	sqpropw_arg.gauge_fix_snk=0;
	sqpropw_arg.store_midprop=0;
	sqpropw_arg.do_half_fermion=0; 
	sqpropw_arg.cg.max_num_iter=99999;
	sqpropw_arg.cg.stop_rsd=1e-8;
	sqpropw_arg.cg.true_rsd=1e-8;
	sqpropw_arg.mob_arg_s=0;
	sqpropw_arg.cg.mass=0.15;

	decode_vml(do_arg);
	decode_vml(meas_arg);	
//	decode_vml(lanc_arg);		//Uncommented by Bigeng
	//decode_vml(cqpropw_arg);
	//decode_vml(lqpropw_arg);
	//decode_vml(fix_gauge_arg);
	//
	//
	#ifdef USE_QUDA	      //Added by Bigeng, for GPU running
        if (!QudaParam.Decode ("quda_arg.vml", "QudaParam")) {
        printf ("Bum quda_arg\n");
        exit (-1);
        }
        #endif	

	//----Initialization end----------------------
	fix_gauge_arg.fix_gauge_kind = FIX_GAUGE_COULOMB_T;
	fix_gauge_arg.hyperplane_start = 0; 
	fix_gauge_arg.hyperplane_step = 1; 
	fix_gauge_arg.hyperplane_num = do_arg.t_sites;
	fix_gauge_arg.stop_cond = 1e-8;
	fix_gauge_arg.max_iter_num = 10000;
	//begin to work
	char dir[512];	
	chdir(meas_arg.WorkDirectory);
	sprintf(dir, "klks_ml_%1.5f_newcps_mob_fac2.0_testall",lqpropw_arg.cg.mass);
	mkdir(dir, 0755);

	GJP.Initialize(do_arg);
	GJP.InitializeExt(doext_arg); //copied from Hackathon_2018
	//init_bfm(&argc, &argv);
	int Ls =  GJP.SnodeSites();
	//set_arg_alpha(Ls,lqpropw_arg.cg.mass,lqpropw_arg.cg.mass,mobius_fac);
	//Fbfm::current_key_mass =lqpropw_arg.cg.mass;



	//main loop for each configuration 
	int low = 2000;//atoi(argv[2]);
	int high = 2000;//atoi(argv[3]);
	int dt = 2000;//atoi(argv[4]);
	for(meas_arg.TrajCur = low; meas_arg.TrajCur <= high; meas_arg.TrajCur += dt) 
	{
		VRB.Result(cname,fname,"configuration %d start :\n",meas_arg.TrajCur);
		double cost = -dclock();
		do_arg.start_seed_value=meas_arg.TrajCur;
		LRG.Initialize();

		ReadGaugeField(meas_arg);
		//ReadRNG(meas_arg);

		//Added new RN using new RNG in CPS_QUDA, code copied from IC_Hackathon
		//char rand_out[512];
		//sprintf(rand_out, "../results/rand_out/rand.%d", meas_arg.TrajCur);
		//LRG.Write(rand_out, 0);
		//std::cout << "LRG Done!" << std::endl;
		//End of RNG in&out		


		//Lattice &lat = LatticeFactory::Create(F_CLASS_BFM, G_CLASS_NONE); 
		GnoneFmobius lat; //Added by Bigeng for GPU running
		//Fbfm &lat_bfm = dynamic_cast<Fbfm &>(lat);
		//lat_bfm.SetBfmArg(Fbfm::current_key_mass );
		sprintf(dir,"./klks_ml_%1.5f_newcps_mob_fac2.0_testall/traj_%d",lqpropw_arg.cg.mass, meas_arg.TrajCur, lqpropw_arg.cg.stop_rsd);
		common_arg.set_filename(dir);
		/*Added by Bigeng
		char lanc_dir[512];
		sprintf(lanc_dir, "/gpfs/mira-fs0/projects/LatticeQCD_3/biwang/16cube/results-%d/huge-data-lanc",resultID);

		lat.init_lanczos(lanc_arg);
		double cost_lanc = -dclock();
	//	lat.run_lanczos();
	//	lanczosWriteParNode(lat, lanc_dir);
	//	lanczosReadParNode(lat, lanc_dir,4500);
		//test_lanczos(lat, 10, 200, 1);//64^3 value used in Ziyuan's calculation.
	//	test_lanczos(lat, 10, 10, 1); //For 16^3 value, we have 100 eigenvectors. used 10*10
		cost_lanc += dclock();
		VRB.Result(cname,fname,"time spend in Lanczos = %f sec. \n", cost_lanc);
		*/ 	//Added by Bigeng
		//
		AlgFixGauge fix(lat, &common_arg, &fix_gauge_arg);
		fix.run();

	//	lat.init_lanczos(lanc_arg);
	//	double cost_lanc = -dclock();
	////	lat.run_lanczos();

		int minsep = 6;
		int nhits = 5; //Use t0 be 60
		int sep=4; //pi pi seperation
		bool do_exact = true;
		AlgKlKs klks(lat, &common_arg, &lqpropw_arg, &sqpropw_arg, &cqpropw_arg, minsep, nhits,sep, do_exact);
		klks.Run1();

		fix.free();
	
		LatticeFactory::Destroy();

		cost += dclock();
		VRB.Result(cname,fname,"configuration %d end, time spend = %f sec. \n",meas_arg.TrajCur, cost);


	}

	VRB.Result(cname,fname,"The End!\n");
	End();
	return 0;
}

void ReadGaugeField(const MeasArg &meas_arg)
{
	const char *cname = "main";
	const char *fname = "ReadGaugeField";

	GnoneFnone lat;
	char lat_file[256];
	sprintf(lat_file, "%s.%d", meas_arg.GaugeStem, meas_arg.TrajCur);
	VRB.Result(cname, fname, lat_file);
	QioArg rd_arg(lat_file, 0.001);//chkprec=0.001
	rd_arg.ConcurIONumber=meas_arg.IOconcurrency;
	ReadLatticeParallel rl;
	rl.read(lat, rd_arg);
	if(!rl.good()) 
		ERR.General(cname, fname, "Failed read lattice %s", lat_file);	
}

void ReadRNG(const MeasArg &meas_arg)
{
	const char *cname = "main";
	const char *fname = "ReadRNG(MeasArg&)";

	char rng_file[256];
	sprintf(rng_file, "%s.%d", meas_arg.RNGStem, meas_arg.TrajCur);
	if (!LRG.Read(rng_file)) 
		ERR.General(cname, fname, "Failed RNG file %s", rng_file);
}


