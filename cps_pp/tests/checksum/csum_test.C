#include<cps.h>

USING_NAMESPACE_CPS

int main(int argc,char *argv[])
{

  //Change working directory
  
  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;

  int num_traj;
  //if(sscanf(argv[1],"%d",&num_traj)!=1) 
  //ERR.General(argv[0],"main()","Error Scanning argv[%d]",1);

  num_traj=10;

  do_arg.x_node_sites = 4;
  do_arg.y_node_sites = 4;
  do_arg.z_node_sites = 4;
  do_arg.t_node_sites = 2;
  do_arg.s_node_sites = 8;
  
  do_arg.x_nodes = SizeX();
  do_arg.y_nodes = SizeY();
  do_arg.z_nodes = SizeZ();
  do_arg.t_nodes = SizeT();
  do_arg.s_nodes = SizeS();
  //  do_arg.s_axis = SCU_S;
  
  
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_PRD;
  do_arg.start_conf_kind = START_CONF_DISORD;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.start_seed_value = 1357;
  do_arg.beta = 5.8;
  do_arg.dwf_height = 1.65;



    VerboseLevelType verbose_level;
    verbose_level=VERBOSE_RESULT_LEVEL;
    
    
    GJP.Initialize(do_arg);
    
    VRB.Level(verbose_level);
    //VRB.ActivateLevel(VERBOSE_SMALLOC_LEVEL);
  
  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------
  CommonArg common_arg_hmd;
  common_arg_hmd.set_filename("hmc_titest.dat");

  CommonArg common_arg_plaq;
  common_arg_plaq.set_filename("plaq_titest.dat");

  NoArg plaq_arg;

  HmdArg hmd_arg;
  hmd_arg.n_frm_masses = 1;
  hmd_arg.frm_mass[0] = 0.025;
  hmd_arg.n_bsn_masses = 1;
  hmd_arg.bsn_mass[0] = 1.0;
  hmd_arg.max_num_iter[0] = 5000;
  hmd_arg.stop_rsd[0] = 1.0E-7;
  hmd_arg.step_size = 0.01;
  hmd_arg.steps_per_traj = 10;
  hmd_arg.metropolis = METROPOLIS_NO;
  hmd_arg.reunitarize = REUNITARIZE_YES;
 
  char filename[256];
  sprintf(filename,"checksum.log.%d",UniqueID());
  FILE *fp = fopen(filename,"w"); 
  CSM.Initialize(); //use default buffer length
  CSM.Activate(CSUM_ALL);
  CSM.Deactivate(CSUM_EVL_MMP);
  CSM.Deactivate(CSUM_MMP_SUM);

  

  //----------------------------------------------------------------
  // Run HMC Phi DWF
  //----------------------------------------------------------------
  
  {
    GwilsonFdwf lat;
    //Set all the sites on the lattice to have the same random number
    
    
    AlgPlaq plaq(lat,&common_arg_plaq,&plaq_arg);
    AlgHmcPhi hmc_phi(lat,&common_arg_hmd,&hmd_arg);
    
    for (int i = 0; i < num_traj; ++i) {
      printf("TRAJ NO. %d\n",i);
      plaq.run();
      hmc_phi.run();
      CSM.Print(fp);
    }
  }
  

 return(0);
}






