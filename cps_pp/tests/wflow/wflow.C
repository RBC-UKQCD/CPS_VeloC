#include<cps.h>
//-------------------------------------------------------------

USING_NAMESPACE_CPS
using namespace std;

const char *cname = "";

HmcArg hmc_arg;

ActionGaugeArg gauge_arg;
ActionQuotientArg quo_arg;
//ActionQuotientArg quo_tm_arg;
ActionRationalQuotientArg rat_quo_arg;
static int Mee;

// define all integrators; vml file arguments will not change
IntABArg ab1_arg;
IntABArg ab2_arg;
IntABArg ab3_arg;

EvoArg evo_arg;
DoArg do_arg;
PbpArg pbp_arg;
NoArg no_arg;
ApeSmearArg ape_arg;

void checkpoint(int traj);

#define decode_vml(arg_name)  do{                                       \
        if ( ! arg_name.Decode(#arg_name".vml", #arg_name) )            \
            ERR.General(cname, fname, "Bad " #arg_name ".vml.\n");      \
    } while(0)

void decode_vml_all(void)
{
    const char *fname = "decode_vml_all()";

    decode_vml(do_arg);
    decode_vml(evo_arg);
//    decode_vml(hmc_arg);
//    decode_vml(gauge_arg);
// decode_vml(quo_tm_arg);
//    decode_vml(quo_arg);
//    decode_vml(rat_quo_arg);
//    decode_vml(ab1_arg);
//    decode_vml(ab2_arg);
    // decode_vml(ab3_arg);
//    decode_vml(pbp_arg);
    // decode_vml(ape_arg);
}

void truncate_it(CommonArg *common_arg, const char stem[], int traj);
void measure_plaq(CommonArg &common_arg);
void measure_wline(CommonArg &common_arg);
void measure_tc(CommonArg &common_arg, int cycle);
void measure_pbp(CommonArg &common_arg, int traj);
void run_hmc(CommonArg &common_arg, int traj, AlgIntAB &int_ab);

void setup(int argc, char *argv[])
{
    const char *fname = "setup()";

    Start(&argc, &argv);
    CommandLine::is(argc,argv);

    if(argc < 2) {
        ERR.General(cname, fname, "Must provide VML directory.\n");
    }

//    if(chdir(argv[1]) != 0) {
    char *dir = CommandLine::arg();
    if(chdir(dir) != 0) {
        ERR.General(cname, fname, "Changing directory to %s failed.\n", dir);
    }

    decode_vml_all();

    if(chdir(evo_arg.work_directory) != 0) {
        ERR.General(cname, fname, "Changing directory to %s failed.\n", evo_arg.work_directory);
    }
    VRB.Result(cname, fname, "Reading VML files successfully.\n");

    GJP.Initialize(do_arg);
    LRG.setSerial();
    LRG.Initialize();

}

void erase_file(const char* filename)
{
  VRB.Result("","erase_file()","erasing %s\n",filename);
  FILE* f = Fopen(filename, "w");
  if(!f)
  ERR.General("","erase_file()","erasing %s failed\n",filename);
  Fclose(f);
}


void measure_wflow(int traj,
                   Float dt,
                   int num_time_steps)
{
  const char* fname = "measure_wflow_scale()";
  VRB.Result(cname, fname, "\n\n\n ************ Starting Wilson flow scale measurement *************\n\n\n");
  Float timer = dclock();

  GwilsonFnone lat;

  // save the unsmeared lattice
  LatticeContainer lat_cont;
  lat_cont.Get(lat);

  double len = dt*num_time_steps;

  CommonArg common_arg,tcharge_arg;
  char filename[512];
//  sprintf(filename, "mkdir -p results-dt%0.3f", dt);
//  system(filename);
  sprintf(filename, "./wflow.%d", dt, traj);
  erase_file(filename);
  common_arg.set_filename(filename);


  NoArg no_arg;

  AlgWilsonFlow wilson_flow(lat, &common_arg, dt);
  AlgActionDensity action_density(lat, &common_arg);
  AlgPlaq plaq(lat, &common_arg, &no_arg);
  AlgTcharge tcharge(lat, &common_arg);

  int wflow_count = 0, action_count = 0, plaq_count = 0, tcharge_count=0;
  Float wflow_time = 0, action_time = 0, plaq_time = 0, tcharge_time=0;
  Float wflow_target_t = dt*(Float)num_time_steps;

  int tc_int = (int) 1.0/dt;
  int step;
   double input_dt = dt;
    double actual_dt; // the dt that's actually been run
    double accumulate_t = 0.;
    while(accumulate_t < wflow_target_t){
      action_density.run();
      tcharge.run();
      actual_dt = wilson_flow.run_adaptive(input_dt,1e-5); // This function returns the actual dt and modify input_dt to the next target value
      accumulate_t += actual_dt;
      if(UniqueID() == 0) printf("This dt = %16.12e, next dt = %16.12e, accumulate t = %16.12e\n", actual_dt, input_dt, accumulate_t);
    }

    action_density.run();
    tcharge.run();



  // restore the unsmeared lattice
  lat_cont.Set(lat);

  VRB.Result(cname, fname, "************ Finished Wilson flow scale measurement *************");
  VRB.Result(cname, fname, "Total time                              = %e seconds\n", dclock() - timer);
  VRB.Result(cname, fname, "Time for %d wilson flow steps           = %e seconds (%e seconds each)\n", wflow_count, wflow_time, wflow_time / wflow_count);
  VRB.Result(cname, fname, "Time for %d action density measurements = %e seconds (%e seconds each)\n", action_count, action_time, action_time / action_count);
  VRB.Result(cname, fname, "Time for %d plaquette measurements      = %e seconds (%e seconds each)\n", plaq_count, plaq_time, plaq_time / plaq_count);
}

int main(int argc, char *argv[])
{
    const char *fname = "main()";

    setup(argc, argv);
//    char *hdw_label=CommandLine::arg();
//    int if_ildg = CommandLine::arg_as_int();
    int wflow_len = CommandLine::arg_as_int();
    

    Float dt= 0.02; // Hard coded, as it really doesn't matter much for adaptive smearing
    int num_time_steps = wflow_len * (int)(1./dt) ;
    VRB.Result(cname,fname,"wflow_len=%f*%d\n",dt,num_time_steps);


    int traj = evo_arg.traj_start;
    for(int conf = 0; conf< evo_arg.gauge_configurations; ++conf) {

//        checkpoint(traj);
    char lat_file[256];
    char rng_file[256];

    Float time = -dclock();


    // Save this config to disk
//    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
{
    GwilsonFnone lat;

    sprintf(lat_file,"%s.%d",evo_arg.gauge_file_stem,traj);
    VRB.Result(cname,fname,"Reading %s\n",lat_file);
    ReadLatticeSerial rl(lat,lat_file);
//    ReadLatticeParallel rl(lat,lat_file);
}

    measure_wflow(traj, dt, num_time_steps);
    time += dclock();
    print_time(fname,"measure_wflow()",time);


//    LatticeFactory::Destroy();
	traj += evo_arg.gauge_unload_period ;
    } //End config loop

    End();

    VRB.Result(cname, fname, "Program ended normally.\n");
    
    return 0;
}

