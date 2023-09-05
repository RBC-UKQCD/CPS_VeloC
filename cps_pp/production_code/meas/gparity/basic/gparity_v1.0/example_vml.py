import generate_vml


generate_vml.Globals.set_bc(["BND_CND_PRD","BND_CND_PRD","BND_CND_PRD","BND_CND_APRD"])
props = generate_vml.JobPropagatorArgs() 

masses = [0.01,0.045]
for m in masses:
    prd = generate_vml.PropagatorArg("prop_m%g_PRD" % m, m , "BND_CND_PRD")
    prd.setWallSource(0)

    pra = generate_vml.PropagatorArg("prop_m%g_APRD" % m, m , "BND_CND_APRD")
    pra.setWallSource(0)

    prf = generate_vml.PropagatorArg("prop_m%g_F" % m, m)
    prf.setCombinationSource( prd.tag, pra.tag, "A_PLUS_B")
    
    prb = generate_vml.PropagatorArg("prop_m%g_B" % m, m)
    prb.setCombinationSource( prd.tag, pra.tag, "A_MINUS_B")
    
    props.addPropagator(prd)
    props.addPropagator(pra)
    props.addPropagator(prf)
    props.addPropagator(prb)


props.write("prop_arg.vml")



contract = generate_vml.GparityContractArg("../configurations/ckpoint_lat.%d", 300, 5, 910, 
                                           generate_vml.FixGaugeArg("FIX_GAUGE_COULOMB_T",0,1,32,1e-08,10000))

#momenta here are in units of pi consistently
#16^3 x 32 lattice
mom_zero = [0.0,0.0,0.0]
mom_pair_zero = [mom_zero,mom_zero]
coswall = 0 #dont use cosine wall sink

#LL and HH combinations
for m in masses:
    prop_f = "prop_m%g_F" % m
    prop_b = "prop_m%g_B" % m

    fn = "../results/bilinears_%g_%g_FF" % (m,m)
    contract.addMeas( generate_vml.ContractionTypeAllBilinears(prop_f, prop_f, [mom_zero], fn) )

    fn = "../results/bilinears_%g_%g_BB" % (m,m)
    contract.addMeas( generate_vml.ContractionTypeAllBilinears(prop_b, prop_b, [mom_zero], fn) )

    fn = "../results/wallsink_bilinears_%g_%g_FF" % (m,m)
    contract.addMeas( generate_vml.ContractionTypeAllWallSinkBilinearsSpecificMomentum(prop_f, prop_f, [mom_pair_zero], coswall, fn) )

    fn = "../results/wallsink_bilinears_%g_%g_BB" % (m,m)
    contract.addMeas( generate_vml.ContractionTypeAllWallSinkBilinearsSpecificMomentum(prop_b, prop_b, [mom_pair_zero], coswall, fn) )


#HL combination
prop_l_f = "prop_m%g_F" % masses[0]
prop_l_b = "prop_m%g_B" % masses[0]
prop_h_f = "prop_m%g_F" % masses[1]
prop_h_b = "prop_m%g_B" % masses[1]

fn = "../results/bilinears_%g_%g_FF" % (masses[1],masses[0])
contract.addMeas( generate_vml.ContractionTypeAllBilinears(prop_h_f, prop_l_f, [mom_zero], fn) )

fn = "../results/bilinears_%g_%g_BB" % (masses[1],masses[0])
contract.addMeas( generate_vml.ContractionTypeAllBilinears(prop_h_b, prop_l_b, [mom_zero], fn) )

fn = "../results/wallsink_bilinears_%g_%g_FF" % (masses[1],masses[0])
contract.addMeas( generate_vml.ContractionTypeAllWallSinkBilinearsSpecificMomentum(prop_h_f, prop_l_f, [mom_pair_zero], coswall, fn) )

fn = "../results/wallsink_bilinears_%g_%g_BB" % (masses[1],masses[0])
contract.addMeas( generate_vml.ContractionTypeAllWallSinkBilinearsSpecificMomentum(prop_h_b, prop_l_b, [mom_pair_zero], coswall, fn) )

#Topological charge
fn = "../results/topq"
contract.addMeas( generate_vml.ContractionTypeTopologicalCharge(20,1,1e-12,-1,0.45, fn) )


contract.write("contract_arg.vml")
                                                                      
