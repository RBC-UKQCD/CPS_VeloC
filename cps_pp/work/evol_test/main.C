#include<config.h>

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

#include<util/lattice.h>
#include<util/lattice/fbfm.h>
#include<util/random.h>
#include<util/time_cps.h>

#include<alg/alg_hmc.h>
#include<alg/common_arg.h>
#include<alg/hmc_arg.h>
#include<alg/hmd_arg.h>

#include<alg/alg_int.h>
#include<alg/int_arg.h>
#include<alg/alg_wline.h>
#include<alg/no_arg.h>
#include<alg/do_arg.h>
#include<alg/alg_plaq.h>
#include<alg/alg_pbp.h>
#include<alg/alg_tcharge.h>
#include<alg/alg_smear.h>
#include<alg/ape_smear_arg.h>
#include<alg/pbp_arg.h>
#include<alg/alg_remez.h>
#include<alg/alg_w_spect.h>
#include<alg/array_arg.h>
#include<alg/alg_fix_gauge.h>

#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/qcdio.h>
#include<util/WriteLatticePar.h>
#include<util/ReadLatticePar.h>
#include<util/qioarg.h>

#include <util/lat_cont.h>

#include <omp.h>
#include <pthread.h>
//-------------------------------------------------------------

#define HAVE_DOEXT

USING_NAMESPACE_CPS
using namespace std;

const char *cname = "";

HmcArg hmc_arg;

ActionGaugeArg gauge_arg;
ActionQuotientArg quo_arg;
ActionFermionArg frm_arg;
ActionRationalQuotientArg rat_quo_arg;

// define all integrators; vml file arguments will not change
IntABArg ab1_arg;
IntABArg ab2_arg;
IntABArg ab3_arg;

EvoArg evo_arg;
DoArg do_arg;
#ifdef HAVE_DOEXT
DoArgExt doext_arg;
#endif
PbpArg pbp_arg;
NoArg no_arg;
ApeSmearArg ape_arg;

void checkpoint(int traj);

#undef encode_vml
#define encode_vml(arg_name, traj) do{                                  \
        char vml_file[256];                                             \
        sprintf(vml_file, #arg_name".%d", traj);                        \
        if( !arg_name.Encode(vml_file, #arg_name) ){                    \
            ERR.General(cname, fname, #arg_name " encoding failed.\n"); \
        }                                                               \
    }while(0)

#undef decode_vml
#define decode_vml(arg_name)  do{                                       \
        if ( ! arg_name.Decode(#arg_name".vml", #arg_name) )            \
            ERR.General(cname, fname, "Bad " #arg_name ".vml.\n");      \
    } while(0)

void decode_vml_all(void)
{
    const char *fname = "decode_vml_all()";

    decode_vml(do_arg);
#ifdef HAVE_DOEXT
    encode_vml(doext_arg,0);
    decode_vml(doext_arg);
    encode_vml(doext_arg,0);
#endif
    decode_vml(hmc_arg);
    decode_vml(evo_arg);
    decode_vml(gauge_arg);
    decode_vml(quo_arg);
    decode_vml(frm_arg);
    decode_vml(rat_quo_arg);
    decode_vml(ab1_arg);
    decode_vml(ab2_arg);
    decode_vml(pbp_arg);
    // decode_vml(ape_arg);
}

void truncate_it(CommonArg *common_arg, const char stem[], int traj);
void measure_plaq(CommonArg &common_arg);
void measure_wline(CommonArg &common_arg);
void measure_tc(CommonArg &common_arg, int cycle);
void measure_pbp(CommonArg &common_arg, int traj);
void run_hmc(CommonArg &common_arg, int traj, AlgIntAB &int_ab);

void init(int *argc, char **argv[]){
}

void setup(int argc, char *argv[])
{
    const char *fname = "setup()";

    Start(&argc, &argv);

    if(argc < 2) {
        ERR.General(cname, fname, "Must provide VML directory.\n");
    }

    if(chdir(argv[1]) != 0) {
        ERR.General(cname, fname, "Changing directory to %s failed.\n", argv[1]);
    }

#ifdef HAVE_DOEXT
    encode_vml(doext_arg,0);
#endif
#ifdef USE_QUDA
   if ( !QudaParam.Decode("quda_arg.vml","QudaParam") )
        { printf("Bum quda_arg\n"); exit(-1);}
#endif
    decode_vml_all();

    if(chdir(evo_arg.work_directory) != 0) {
        ERR.General(cname, fname, "Changing directory to %s failed.\n", evo_arg.work_directory);
    }
    VRB.Result(cname, fname, "Reading VML files successfully.\n");

    GJP.Initialize(do_arg);
#ifdef HAVE_DOEXT
    GJP.InitializeExt(doext_arg);
#endif
//    GJP.SetMobius(32./12.);
    GJP.ZMobius_PC_Type(ZMOB_PC_ORIG );
//    GJP.ZMobius_PC_Type(ZMOB_PC_SYM1 );
    LRG.setSerial();
    LRG.Initialize();
 int threads;
 char * nthr_str = getenv("OMP_NUM_THREADS");
 if(nthr_str) sscanf(nthr_str,"%d",&threads);
 if(!UniqueID()) printf("nthreads=%d\n",threads);
 GJP.SetNthreads(threads);
 printf("GJP.Nthreads()=%d\n",GJP.Nthreads());
// exit(-43);

   


    init(&argc, &argv);

     if ( do_arg.start_conf_kind != START_CONF_FILE ) {
         checkpoint(evo_arg.traj_start);
     }
}

int main(int argc, char *argv[])
{
    const char *fname = "main()";

    setup(argc, argv);

    //////////////////////////////////////////////////////////////////////
    // creating actions and integrators
    AlgMomentum mom;
    AlgActionGauge gauge(mom, gauge_arg);
    AlgActionQuotient quotient(mom, quo_arg);
    AlgActionFermion fermion(mom, frm_arg);
    AlgActionRationalQuotient rat_quo(mom, rat_quo_arg);

    IntABArg sum_arg;
    sum_arg.A_steps = 1;
    sum_arg.B_steps = 1;
    sum_arg.level = EMBEDDED_INTEGRATOR;
    AlgIntSum sum(quotient, rat_quo, sum_arg);

    //!< Construct numerical integrators
    AlgIntAB &ab1 = AlgIntAB::Create(mom, gauge, ab1_arg);
    AlgIntAB &ab2 = AlgIntAB::Create(ab1, rat_quo,   ab2_arg);
//    AlgIntAB &ab2 = AlgIntAB::Create(ab1, quotient,   ab2_arg);
//    AlgIntAB &ab2 = AlgIntAB::Create(ab1, fermion,   ab2_arg);
//    AlgIntAB &ab2 = AlgIntAB::Create(ab1, sum,   ab2_arg);
    // AlgIntAB &ab3 = AlgIntAB::Create(ab2, sum,         ab3_arg);
    //////////////////////////////////////////////////////////////////////

    int traj = evo_arg.traj_start;
    for(int conf = 0; conf< evo_arg.gauge_configurations; ++conf) {
        CommonArg common_arg_plaq;
        CommonArg common_arg_pbp;
        CommonArg common_arg_hmc;
        // CommonArg common_arg_top;
        CommonArg common_arg_wline;

        truncate_it(&common_arg_plaq , evo_arg.plaquette_stem      , traj);
        // truncate_it(&common_arg_top  , "./results/alg_top/top"    , traj);
        truncate_it(&common_arg_wline, "./wline", traj);
        truncate_it(&common_arg_pbp  , evo_arg.pbp_stem            , traj);
        truncate_it(&common_arg_hmc  , evo_arg.evo_stem            , traj);

        // Inner trajectory loop
        for(int i = 0; i < evo_arg.gauge_unload_period; ++i, ++traj) {

            measure_plaq(common_arg_plaq);
            // measure_wline(common_arg_wline);
            // measure_tc(common_arg_top, 20);
            measure_pbp(common_arg_pbp, traj);

            run_hmc(common_arg_hmc, traj, ab2);
        }//End of inter-cfg sweep

        checkpoint(traj);

    } //End config loop

    // AlgIntAB::Destroy(ab3);
    AlgIntAB::Destroy(ab2);
    AlgIntAB::Destroy(ab1);

    End();

    VRB.Result(cname, fname, "Program ended normally.\n");
    
    return 0;
}

void checkpoint(int traj)
{
    const char *fname="checkpoint()";

    char lat_file[256];
    char rng_file[256];

    Float time = -dclock();

    // Save this config to disk
    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);

    sprintf(lat_file,"%s.%d",evo_arg.gauge_file_stem,traj);
    QioArg wt_arg(lat_file,0.001);

    wt_arg.ConcurIONumber=evo_arg.io_concurrency;
    WriteLatticeParallel wl;
    wl.setHeader(evo_arg.ensemble_id,evo_arg.ensemble_label,traj);
    wl.write(lat,wt_arg);

    if(!wl.good())
        ERR.General(cname,fname,"Failed write lattice %s",lat_file);

    LatticeFactory::Destroy();

    // Save the RNG's
    sprintf(rng_file,"%s.%d",evo_arg.rng_file_stem,traj);
    if ( !LRG.Write(rng_file) )
        ERR.General(cname,fname,"Failed write RNG file %s",rng_file);

    // Update the parameter files for restart
    do_arg.start_seed_filename = rng_file;
    do_arg.start_seed_kind = START_SEED_FILE;
    do_arg.start_conf_filename = lat_file;
    do_arg.start_conf_kind = START_CONF_FILE;
    evo_arg.traj_start     = traj;

    encode_vml(hmc_arg, traj);
    encode_vml(gauge_arg, traj);
    encode_vml(quo_arg, traj);
    encode_vml(rat_quo_arg, traj);
    encode_vml(ab1_arg, traj);
    encode_vml(ab2_arg, traj);
    // encode_vml(ab3_arg, traj);
    encode_vml(pbp_arg, traj);
    encode_vml(do_arg, traj);
    encode_vml(evo_arg, traj);

    time += dclock();
    print_flops("","checkpoint()",0,time);
}

void truncate_it(CommonArg *common_arg, const char stem[], int traj)
{
    char fnbuf[1024];
    sprintf(fnbuf, "%s.%d", stem, traj);
    FILE *truncate_it = Fopen(fnbuf, "w");
    Fclose(truncate_it);
    common_arg->set_filename(fnbuf);
}

void measure_plaq(CommonArg &common_arg)
{
    const char *fname = "measure_plaq()";

    Float dtime = -dclock();

    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_WILSON);
    AlgPlaq plaq(lat, &common_arg, &no_arg);
    plaq.run();
    LatticeFactory::Destroy();

    dtime += dclock();
    print_flops("AlgPlaq", "run()", 0, dtime);	
}

void measure_wline(CommonArg &common_arg)
{
    const char *fname = "measure_wline()";

    Float dtime = -dclock();

    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_WILSON);
    AlgWline wline(lat, &common_arg, &no_arg);
    wline.run();
    LatticeFactory::Destroy();

    dtime += dclock();
    print_flops("AlgWline", "run()", 0, dtime);	
}

void measure_tc(CommonArg &common_arg, int cycle)
{
    const char *fname = "measure_tc()";

    Float dtime = -dclock();

    // calculate the topological charge. Need to copy the lattice since
    // we need to smear it first. Use Chulwoo's "lattice container"
    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_WILSON);
    LatticeContainer lat_cont;
    // copy lattice
    lat_cont.Get(lat);

    //----- mess up lattice -------------------------
    AlgApeSmear ape(lat, &common_arg, &ape_arg);
    AlgTcharge  tcharge(lat, &common_arg);
    for (int i = 0; i < cycle; ++i) {
        VRB.Result(cname,fname,"%i\n",i);
        VRB.Result(cname,fname,"   running tcharge\n"); tcharge.run();
        VRB.Result(cname,fname,"   running ape\n"); ape.run();
        VRB.Result(cname,fname,"   running ape\n"); ape.run();
        VRB.Result(cname,fname,"   running ape\n"); ape.run();
    }
    tcharge.run();
    // restore the lattice
    lat_cont.Set(lat);
    LatticeFactory::Destroy();

    dtime += dclock();
    print_flops("AlgTcharge", "run()", 0, dtime);
}

void measure_pbp(CommonArg &common_arg, int traj)
{
    const char *fname = "measure_pbp()";

    // fix pbp_arg
    pbp_arg.src_u_s = 0;
    pbp_arg.src_l_s = GJP.SnodeSites() * GJP.Snodes() - 1;
    pbp_arg.snk_u_s = GJP.SnodeSites() * GJP.Snodes() - 1;
    pbp_arg.snk_l_s = 0;

    const int g_int = evo_arg.gauge_unload_period;
    if (traj % g_int == 0 && evo_arg.measure_pbp) {
        Float dtime = -dclock();

        LRGState rng_state;
        rng_state.GetStates();

        Lattice &lat = LatticeFactory::Create(F_CLASS_BFM, G_CLASS_NONE);
        AlgPbp pbp(lat, &common_arg, &pbp_arg);

        for(int pbp_counter = 0; pbp_counter < evo_arg.measure_pbp; pbp_counter++) {
            pbp.run();
        }
        LatticeFactory::Destroy();
        rng_state.SetStates();

        dtime += dclock();
        print_flops("AlgPbp", "run()", 0, dtime);	
    }
}

void run_hmc(CommonArg &common_arg, int traj, AlgIntAB &int_ab)
{
    const char *fname = "run_hmc()";

    Float dtime = -dclock();

    if ( (evo_arg.reproduce_interval > 0) &&
         (traj % evo_arg.reproduce_interval) == 0 ) {
        VRB.Result(cname,fname,"Running traj %d with reproduction\n",traj);
        hmc_arg.reproduce = REPRODUCE_YES;
    } else {
        VRB.Result(cname,fname,"Running traj %d without reproduction\n",traj);
        hmc_arg.reproduce = REPRODUCE_NO;	
    }
    
    AlgHmc hmc(int_ab, common_arg, hmc_arg);
    hmc.run();

    dtime += dclock();
    print_flops("AlgHmc", "run()", 0, dtime);
}
