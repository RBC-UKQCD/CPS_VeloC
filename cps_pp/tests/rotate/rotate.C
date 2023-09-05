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

    char *dir = CommandLine::arg();
    if(chdir(dir) != 0) {
        ERR.General(cname, fname, "Changing directory to %s failed.\n", dir);
    }

    decode_vml_all();

//    if(chdir(evo_arg.work_directory) != 0) {
//        ERR.General(cname, fname, "Changing directory to %s failed.\n", evo_arg.work_directory);
 //   }
    VRB.Result(cname, fname, "Reading VML files successfully.\n");

    GJP.Initialize(do_arg);
    LRG.setSerial();
    LRG.Initialize();

}

int main(int argc, char *argv[])
{
    const char *fname = "main()";

    setup(argc, argv);
    char *hdw_label=CommandLine::arg();
    int if_ildg = CommandLine::arg_as_int();
    {
      GnoneFnone lat;
      WriteLatticeParallel (lat,"ckpoint_lat.out");
    }

    int traj = evo_arg.traj_start;
    for(int conf = 0; conf< evo_arg.gauge_configurations; ++conf) {

//        checkpoint(traj);
    char lat_out[256];
    char lat_in[256];

    Float time = -dclock();


    // Save this config to disk
    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);

    sprintf(lat_out,"%s.%d",evo_arg.work_directory,traj);
    FILE *fp = fopen(lat_out,"r");
    VRB.Result(cname,fname,"Output file: %s fp=%p \n",lat_out,fp);

if (!fp){
    sprintf(lat_in,"%s.%d",evo_arg.gauge_file_stem,traj);
    VRB.Result(cname,fname,"reading %s\n",lat_in);
//    ReadLatticeParallel rl(lat,lat_out);
//    ReadLatticeSerial rl(lat,lat_in);
    QioArg rd_arg(lat_in,0.001);
    rd_arg.ConcurIONumber=evo_arg.io_concurrency;
    ReadLatticeParallel rl;
    rl.readXYTZ(lat,rd_arg);
//    rl.read(lat,rd_arg);
//    exit(-42);

//    sprintf(lat_out,"ckpoint_lat.%d",traj);
    QioArg wt_arg(lat_out,0.001);
    wt_arg.ConcurIONumber=evo_arg.io_concurrency;
    WriteLatticeSerial wl;
    wl.hd.creation_date = rl.hd.creation_date;
    wl.hd.setHeader(evo_arg.ensemble_id,evo_arg.ensemble_label,traj,evo_arg.creator,hdw_label);
    cout << "creation date "<<wl.hd.creation_date<<endl;
    cout << "creation hardware "<<wl.hd.creator_hardware<<endl;
    wl.write(lat,wt_arg);

    if(!wl.good())
        ERR.General(cname,fname,"Failed writing lattice %s",lat_out);
} else
    fclose (fp);

    VRB.Result(cname,fname,"reading %s\n",lat_out);
    ReadLatticeSerial rl(lat,lat_out);
    if(!rl.good())
        ERR.General(cname,fname,"Failed to read lattice %s",lat_out);

    LatticeFactory::Destroy();
	traj += evo_arg.gauge_unload_period ;
    } //End config loop

    End();

    VRB.Result(cname, fname, "Program ended normally.\n");
    
    return 0;
}

