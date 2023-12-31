#include<cps.h>

USING_NAMESPACE_CPS

int main(int argc,char *argv[]) {

    Start(&argc,&argv);
    DoArg do_arg;

    do_arg.x_node_sites = 4;
    do_arg.y_node_sites = 4;
    do_arg.z_node_sites = 4;
    do_arg.t_node_sites = 4;
    do_arg.x_nodes = SizeX();
    do_arg.y_nodes = SizeY();
    do_arg.z_nodes = SizeZ();
    do_arg.t_nodes = SizeT();
    do_arg.x_bc = BND_CND_PRD;
    do_arg.y_bc = BND_CND_PRD;
    do_arg.z_bc = BND_CND_PRD;
    do_arg.t_bc = BND_CND_APRD;
    do_arg.start_conf_kind = START_CONF_DISORD;
    do_arg.start_seed_kind = START_SEED_FIXED_UNIFORM;

    
    
    GJP.Initialize(do_arg);

    VRB.Level(VERBOSE_FLOW_LEVEL);
    VRB.DeactivateLevel(VERBOSE_SMALLOC_LEVEL);    

    GJP.Initialize(do_arg);

    CommonArg common;
    NoArg none;

    common.set_filename("plaqdump");

    GwilsonFnone lat;
    AlgPlaq plaquette(lat, &common, &none);
    plaquette.run();
    printf(" plaquette = %f\n", lat.SumReTrPlaqNode()/(GJP.VolNodeSites()*3.0*6.0));
    End();
    return 0;

}





  


