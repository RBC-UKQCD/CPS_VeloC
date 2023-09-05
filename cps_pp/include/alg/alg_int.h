// -*- mode:c++; c-basic-offset:4 -*-
#ifndef INCLUDED_ALG_INT_H
#define INCLUDED_ALG_INT_H

#include<vector>
#include<set>
#include<string>
//#ifdef HAVE_VELOC
//moved to util/checksum.h>
#if 0
//#include<veloc.h>
#include<veloc.hpp>
extern veloc::client_t *veloc_cps;
#endif
//------------------------------------------------------------------
/*!\file
  \brief  Definitions of the AlgInt class and derived classes.

*/
//------------------------------------------------------------------


#include<config.h>
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/alg_meas.h>
#include <alg/cg_stats.h>
#include <util/checksum.h>
#include <alg/int_arg.h>
#include <alg/remez_arg.h>
#include <util/qcdio.h>
#include <alg/alg_eig.h> 
#include<util/time_cps.h>
CPS_START_NAMESPACE

/*!<
  The super class for all intgrators.
*/
class AlgInt :public CheckPointer {

private:
    const char *cname;
//    static int veloc_id;

protected:
    //!< the current trajectory number
    IntegratorType int_type;
#if 0
    static std::set<int> md_veloc_all;
    static std::set<int> phi_veloc_all;
    static void   *veloc_p[100];
    static size_t veloc_m_size[100];
    static int veloc_d_size[100];


    static int VeloCCounter(){
       int tmp = veloc_id++; 
       return tmp;
    }
#endif
public:
    int s_ckpt;

    static int traj_num;
  
    AlgInt();
    virtual ~AlgInt();

    //!< method to do heatbath (if necessary)
    virtual void heatbath() = 0;

    // Calculate preliminary force for force gradient
    virtual void prepare_fg(Matrix * force, Float dt_ratio);
    //!< run method evolves the integrator
    virtual void evolve (Float dt, int steps) = 0;

    //!< method to calculate energy (if necessary)
    virtual Float energy() = 0;

    //!< method to return the cost of the integrator (wrt cg iterations)
    virtual void cost(CgStats*) = 0;

    //!< method to reverse the direction of evolution (i.e. flip momenta)
    virtual void reverse() = 0;

    //!< method used to reinitialise the integrator
    virtual void init() = 0;

  //!< when G-parity active, this function (which should be do something only by AlgActionGauge) copy-conjugates the gauge field
  virtual void copyConjLattice();
int  getVer(const char *_label){
   const char *fname="getVer()";
   int veloc_v=-1;
   Float dtime2 = -dclock();
#ifdef HAVE_VELOC
  {
//       veloc_v = veloc_cps->restart_test(veloc_label.c_str(),_steps+10);
    VRB.Debug(cname,fname,"VELOC_Restart_test %s\n",_label);
    veloc_v =  veloc_cps->restart_test(_label,100);
    VRB.Result(cname,fname,"Veloc label version: %s %d\n", _label,veloc_v);
  }
#else 
  for(int i=0;i<1000;i++){
  std::string str(_label);
  str +="_"+std::to_string(i);
  std::string fn= str+"/"+str;
  FILE *fp= Fopen(ADD_ID, fn.c_str(),"r");
  if(fp){ veloc_v=i; Fclose(fp); }
  else if (i>10) break;
  }

  Float i_f=veloc_v;
  glb_min(&i_f);
  if(fabs(veloc_v-i_f)>1e-8) ERR.General(cname,fname,"i=%g i_f=%g\n",veloc_v,i_f);
  VRB.Result(cname,fname,"CPScheckpoint label version: %s  %d\n", _label,veloc_v);
#endif
   dtime2 += dclock();
   return veloc_v;
}
int  getVer(const char *cname,const char *fname){
       std::string veloc_label, veloc_label2;
       veloc_label = "MD_traj"+std::to_string(traj_num);
       return getVer(veloc_label.c_str());
}


int protect (const char *cname,const char *fname, void *p, size_t m_size, int d_size){
 int veloc_id= -1;
 veloc_id= VeloCCounter();
 Float dtime = -dclock();
#ifdef HAVE_VELOC
  VRB.Result(cname,fname,"VELOC_Mem_protect \n");
// VELOC_Mem_protect (  veloc_id  , p, m_size, d_size );
 veloc_cps->mem_protect (  veloc_id  , p, m_size, d_size );
#else
 veloc_p[veloc_id]=p;
 veloc_m_size[veloc_id]=m_size;
 veloc_d_size[veloc_id]=d_size;
#endif
 VRB.Result(cname,fname,"VELOC id size=%d %d :%d\n",veloc_id,m_size,d_size);
 dtime +=dclock();
 print_flops(fname,"protect()",0,dtime);
 return veloc_id;
}

};

/*!< 
  The super class for all numerical integrators.
*/
class AlgIntAB : public AlgInt {

private:
    const char *cname;

protected:
    AlgInt *A;
    AlgInt *B;
    int A_steps, B_steps;
    int A_calls, B_calls;
  
    IntegratorLevel level; //!< Is this the top level integrator?
    unsigned long step_cnt;
    IntABArg *ab_arg;

public:
    AlgIntAB(AlgInt &A, AlgInt &B, IntABArg &);
    virtual ~AlgIntAB();

    void heatbath();

    Float energy();

    //!< evolve method evolves the integrator
    virtual void evolve(Float dt, int steps) = 0;

    void cost(CgStats*);
    void reverse();

    //!< Dummy method
    void init();

    //!< AlgIntAB factory
    static AlgIntAB& Create(AlgInt &A, AlgInt &B, IntABArg &ab_arg);
    static void Destroy(AlgIntAB&);

    int checkpoint( const char *cname, const char *fname, int &i, int _steps) {
       int veloc_v=-1;
       Float dtime = -dclock();
       std::string veloc_label, veloc_label2;
       veloc_label = "MD_traj"+std::to_string(traj_num);
       veloc_label2 = "RNG_traj"+std::to_string(traj_num);
#ifdef HAVE_VELOC
       veloc_v = veloc_cps->restart_test(veloc_label.c_str(),_steps+10);
       VRB.Result(cname,fname,"Veloc label version i steps : %s %d %d %d\n",veloc_label.c_str(),veloc_v,i,_steps);
       if(veloc_v == VELOC_FAILURE  || (veloc_v < (i+2) ) ){ // no restart, checkpoint
//        VELOC_Checkpoint_begin(veloc_label, i+2 ); //leaving room for heatbath
//        VELOC_Checkpoint_selective(VELOC_CKPT_SOME,md_veloc_all.data(),md_veloc_all.size());
//        VELOC_Checkpoint_end(1 );
        veloc_cps->checkpoint_begin(veloc_label.c_str(), i+2 ); //leaving room for heatbath
        veloc_cps->checkpoint_mem(VELOC_CKPT_SOME,md_veloc_all);
        veloc_cps->checkpoint_end(1);
       if(veloc_v == VELOC_FAILURE  || (veloc_v <2 ) ){ // no restart, checkpoint
       VRB.Result(cname,fname, "Veloc Checkpoint label2: %s .\n",veloc_label2);
        LRG.Write(veloc_label2.c_str());
        }
        } else {
       assert(veloc_v <(_steps+2));
        i = veloc_v-2;
//        VELOC_Restart_begin (veloc_label, i+2 );
//        VELOC_Recover_selective(VELOC_CKPT_SOME,md_veloc_all.data(),md_veloc_all.size());
//        VELOC_Restart_end(1 );
        veloc_cps->restart_begin(veloc_label, i+2 ); //leaving room for heatbath
        veloc_cps->recover_mem(VELOC_CKPT_SOME,md_veloc_all);
        veloc_cps->restart_end(1);
       VRB.Result(cname,fname, "Veloc Recover label2: %s .\n",veloc_label.c_str());
        LRG.Read(veloc_label2.c_str());
       }
       dtime += dclock();
       print_flops(fname,"VeloC()",0,dtime);
#else 
       veloc_v = getVer(veloc_label.c_str());
       VRB.Result(cname,fname,"CPSCheckpoint label version i steps : %s %d %d %d\n",veloc_label.c_str(),veloc_v,i,_steps);
       if(veloc_v == -1  || (veloc_v < (i+2) ) ){ // no restart, checkpoint
       std::string string_dir= veloc_label+"_"+std::to_string(i+2);
       std::string string_cmd= "mkdir -p "+string_dir;
       if(!UniqueID()) system(string_cmd.c_str());
       Barrier();
//       Float temp=0;glb_sum(&temp);
       std::string fn= string_dir+"/"+string_dir;
       FILE *fp = Fopen(ADD_ID,fn.c_str(),"w");
//       for (int j=0;j<md_veloc_all.size();j++){
//          int id = md_veloc_all[j];
       int j=0;
       for ( auto id : md_veloc_all){
          VRB.Result(cname,fname,"md_veloc_all[%d]=%d\n",j,id);j++;
          fwrite(veloc_p[id],veloc_m_size[id],veloc_d_size[id],fp);
       }
       Fclose(fp);
       if(veloc_v == -1  || (veloc_v <2 ) ){ // no restart, checkpoint
       VRB.Result(cname,fname, "Veloc Checkpoint label2: %s .\n",veloc_label2.c_str());
        LRG.Write(veloc_label2.c_str());
        }
        } else {
       assert(veloc_v <(_steps+2));
        i = veloc_v-2;
       VRB.Result(cname,fname,"i = %d\n",i);
//       std::string fn= veloc_label+"_"+std::to_string(veloc_v);
       std::string string_dir= veloc_label+"_"+std::to_string(veloc_v);
       std::string string_cmd= "mkdir -p "+string_dir;
       if(!UniqueID()) system(string_cmd.c_str());
       Barrier();
//       Float temp=0;glb_sum(&temp);
       std::string fn= string_dir+"/"+string_dir;
       FILE *fp = Fopen(ADD_ID,fn.c_str(),"r");
//       for (int j=0;j<md_veloc_all.size();j++){
//          int id = md_veloc_all[j];
       int j=0;
       for ( auto id : md_veloc_all){
          VRB.Result(cname,fname,"md_veloc_all[%d]=%d\n",j,id);j++;
          fread(veloc_p[id],veloc_m_size[id],veloc_d_size[id],fp);
       }
       Fclose(fp);
       VRB.Result(cname,fname, "CPSCheckpoint Recover label2: %s .\n",veloc_label2.c_str());
        LRG.Read(veloc_label2.c_str());
       }
#endif
       return i;
    }

};

/*!< 
  An implementation of the 2nd order leapfrog integrater 
  (with 1 force evaluation per step).
*/
class AlgIntLeap : public AlgIntAB {

private:
    const char *cname;

public:
    AlgIntLeap(AlgInt &A, AlgInt &B, IntABArg &);
    virtual ~AlgIntLeap();

    //!< evolve method evolves the integrator
    void evolve(Float dt, int steps);
};

/*!< 
  An implementation of the 2nd order Omelyan integrater 
  (with 2 force evaluations per step).
*/
class AlgIntOmelyan : public AlgIntAB {

private:
    const char *cname;
    Float lambda;

public:
    AlgIntOmelyan(AlgInt &A, AlgInt &B, IntABArg &);
    virtual ~AlgIntOmelyan();

    //!< evolve method evolves the integrator
    void evolve(Float dt, int steps);
};

/*!< 
  An implementation of the 4th order Campostrini integrator 
  (with 3 force evaluations per step).
*/
class AlgIntCampostrini : public AlgIntAB {

private:
    const char *cname;
    Float sigma;
    Float epsilon;

public:
    AlgIntCampostrini(AlgInt &A, AlgInt &B, IntABArg &);
    virtual ~AlgIntCampostrini();

    //!< evolve method evolves the integrator
    void evolve(Float dt, int steps);
};

/*!< 
  An implementation of the 4th order Omelyan integrator with 4
  force evaluations per step.
*/
class AlgIntOmelyan44 : public AlgIntAB {

private:
    const char *cname;
    Float rho;
    Float theta;
    Float lambda;

public:
    AlgIntOmelyan44(AlgInt &A, AlgInt &B, IntABArg &);
    virtual ~AlgIntOmelyan44();

    //!< evolve method evolves the integrator
    void evolve(Float dt, int steps);
};

/*!< 
  An implementation of the 4th order Omelyan integrator with 5 force
  evaluations per step.
*/
class AlgIntOmelyan45 : public AlgIntAB {

private:
    const char *cname;
    Float theta;
    Float rho;
    float lambda;
    Float mu;

public:
    AlgIntOmelyan45(AlgInt &A, AlgInt &B, IntABArg &);
    virtual ~AlgIntOmelyan45();

    //!< evolve method evolves the integrator
    void evolve(Float dt, int steps);
};

/*!< 
  A fourth order force gradient integrator.  
*/
class AlgIntForceGrad : public AlgIntAB {

private:
    const char *cname;
    Float lambda;
    Float xi;
    Float theta;
    Float chi;
    int g_size;

protected:
    void evolve_fg(AlgInt * which_int, Float fg_dt, Float dt);

public:
    AlgIntForceGrad(AlgInt &A, AlgInt &B, IntABArg &);
    virtual ~AlgIntForceGrad();

    //!< evolve method evolves the integrator
    void evolve(Float dt, int steps);
};

/*!< 
  Produces the composite integrator of A and B.  If [A,B] != [B,A]
  then obviously order matters, the first argument passed to the
  constructor is the first integrator that is called.
*/
class AlgIntSum : public AlgIntAB {

private:
    const char *cname;

public:
    AlgIntSum(AlgInt &A, AlgInt &B, IntABArg &);
    virtual ~AlgIntSum();

    // Calculate preliminary force for force gradient
    virtual void prepare_fg(Matrix * force, Float dt_ratio);
    //!< evolve method evolves the integrator
    void evolve(Float dt, int steps);
};

/*!< 
  Super class of all Hamiltonian constituents
*/
class AlgHamiltonian : public AlgInt {

protected:
    int g_size;
    Float h_init;

private:
    const char *cname;

public:
    AlgHamiltonian();
    virtual ~AlgHamiltonian();

    virtual void heatbath() = 0;
    virtual Float energy() = 0;
    virtual void evolve(Float dt, int steps) = 0;
    virtual void cost(CgStats*) = 0;




};

/*!< 
  Class describing the conjugate momentum contribution to the Hamiltonian.
*/
class AlgMomentum : public AlgHamiltonian {

private:
    const char *cname;

    //!< This class tracks the MD time since it updates the gauge field
    char *md_time_str;
    Matrix *mom;

public:
//    int veloc_id;
    AlgMomentum();
    virtual ~AlgMomentum();

    void heatbath();

    Float energy();

    //!< evolve method evolves the integrator
    void evolve(Float dt, int steps);

    void cost(CgStats*);

    Matrix* getMom();

    void reverse();

    //!< Dummy method
    void init();
    // for HMC restarting. Not working yet
    void SaveState(std::string);
    void LoadState(std::string);
};

/*!< 
  Super class of all possible actions
*/
class AlgAction : public AlgHamiltonian {

private:
    const char *cname;

protected:
    Matrix *mom;
    ForceMeasure force_measure;
    char *force_label;
    ForceArg Fdt;

public:
    AlgAction();
    AlgAction(AlgMomentum &mom, ActionArg &action_arg);
    ~AlgAction();

    virtual void heatbath() = 0;
    virtual Float energy() = 0;
    virtual void evolve(Float dt, int steps) = 0;
    virtual void cost(CgStats*) = 0;
    void reverse();

};

/*!< 
  Super class of all possible bilinear actions.
*/
class AlgActionBilinear : public AlgAction {

private:
    const char *cname;

protected:
    ActionBilinearArg *bi_arg;

    int n_masses;
    FclassType fermion;

    //!< An array which stores the values of the masses
    Float *mass;

    //!< Maximum number of cg iterations
    int *max_num_iter;

    //!< Number of lattice sites
    int f_sites;
    //!< Number of Vectors in a Vector array
    int f_vec_count;
    //!< Number of Floats in a Vector array
    size_t f_size;
    //!< Number of checkerboards
    int Ncb;

    //!< The conjugate gradient statistics
    CgStats cg_stats;
    int cg_iter;

    //!< Pseudofermion fields, one for each mass
    Vector **phi;

    int md_steps;
    bool skip_force;
//    std::set<int> phi_veloc;
    std::vector<int> phi_veloc;

public:
    AlgActionBilinear();
    AlgActionBilinear(AlgMomentum &, ActionBilinearArg &);
    virtual ~AlgActionBilinear();
    void cost(CgStats*);
    void updateCgStats(CgArg*);

    int getNmass();
    Float getMass(int);
    FclassType getFermion();

    virtual void heatbath() = 0;
    virtual Float energy() = 0;
    virtual void evolve(Float dt, int steps) = 0;

    Vector** getPhi(){ return phi; }

    void init();
    // for HMC restarting. Not working yet
    void SaveState(std::string);
    void LoadState(std::string);

};

/*!< 
  Class describing  rational bilinear action

  S = phi^dag r^2 phi

  This class can be used to simulate <4 flavours of staggered fermions
  or <2 flavours of Wilson fermions.
*/
class AlgActionRational : public AlgActionBilinear {

private:
    const char *cname;
    int **fractionSplit;
    int **splitCheck;
    ActionRationalArg *rat_arg;

    //!< This is where the rational parameters are stored
  
    // currently force gradient calculation is using the same poles as
    // used by the usual MD step.
    RemezArg *remez_arg_md;
    RemezArg *remez_arg_mc;

    //protected:
public: //CK: temporarily public for testing
    //!< Has any evolution taken place?
    int evolved;

    //!< Has the heatbath been evaluated?
    int heatbathEval;

    //!< Has the energy been evaluated?
    int energyEval;

//    Float h_init;

    //!<  frm_cg_arg_fg is specific to force gradient integrator and has no use otherwise.
    CgArg ***frm_cg_arg_fg;
    CgArg ***frm_cg_arg_md;
    CgArg ***frm_cg_arg_mc;
    //!< Pointer to an array of structures containing solver parameters.
    /*!<
      These are the parameters corresponding to each of the dynamical fermion
      masses.      
    */

    Vector** frmn;
    //!< Array of vectors
    /*!< These will hold the solutions from the solves. */
  
    Vector** frmn_d;
    //!< Array of vectors
    /*!< These will hold the solutions from the solves multiplied by
      the D-slash operator. */ 
  
    Vector** frmn_tmp;
    //!< Used for Asqtad fraction splitting

    int total_size;
    //!< The sum of the rational approximation degrees used for the force

    int max_size;
    //!< The maximum degree of rational approximation used for the force
  
    Float *all_res;
    //!< An array holding the residues - used for asqtad force

    EigArg eig_arg;
    //!< AlgEig parameters

    char eig_file[256];
    //!< Used to store eigenvalue filename

    CommonArg ca_eig;
    //!< Used to store eigenvalue filename

    Float **lambda_low;
    Float **lambda_high;
    //!< Used to store calculated eigenvalue bounds

    //!< Automatic generation of the rational approximation.
    //Note masses are only used to bound the eigenvalues if RATIONAL_BOUNDS_AUTOMATIC is switched on
    //hence no modification is required for twisted mass fermions
    void generateApprox(Float *mass, RemezArg **remez_arg_md, 
                        RemezArg **remez_arg_mc, RationalDescr *rat);
  
    //<! Free approximations
    void destroyApprox(RemezArg *remez_arg_md, RemezArg *remez_arg_mc);
  
    //!< Allocate and setup cg arguments
    //CK: added for twisted mass fermions. For non-twisted mass fermions pass either NULL or a float array of size n_masses for the epsilon parameter
    void generateCgArg(Float *mass,
		       Float *epsilon,
                       CgArg **** cg_arg_fg,
                       CgArg **** cg_arg_md, 
                       CgArg **** cg_arg_mc, const char *label, 
                       RationalDescr *rat_des); 
    //CK: passes epsilon=NULL to the above. Has a catch to prevent mistakenly using this function for twisted mass fermions 
    void generateCgArg(Float *mass,
                       CgArg **** cg_arg_fg,
                       CgArg **** cg_arg_md, 
                       CgArg **** cg_arg_mc, const char *label, 
                       RationalDescr *rat_des);


    //<! Free cg args
    void destroyCgArg(CgArg ***cg_arg_fg,
                      CgArg ***cg_arg_md,
                      CgArg ***cg_arg_mc,
                      const char *label, RemezArg *remez_arg_md,
                      RemezArg *remez_arg_mc);

    //!< Automatic generation of required EigArg
    void generateEigArg(EigenDescr eigen);

    //!< Free EigArg
    void destroyEigArg();

    void checkApprox(Float *mass, Float *epsilon, RemezArg *remez_arg, EigenDescr eigen); //For non-twisted mass fermions pass either NULL or a float array of size n_masses for the epsilon parameter
    void checkApprox(Float *mass, RemezArg *remez_arg, EigenDescr eigen); //passes epsilon=NULL to the above. Has a catch to prevent mistakenly using this function for twisted mass fermions
    //!< Check that the approximation is still valid

public:

    AlgActionRational();
    AlgActionRational(AlgMomentum &mom, ActionRationalArg &rat_arg, int traj_num=0);
    AlgActionRational(AlgMomentum &mom, ActionBilinearArg &bi_arg);
    virtual ~AlgActionRational();

    void heatbath();
    Float energy();

    // Calculate preliminary force for force gradient
    virtual void prepare_fg(Matrix * force, Float dt_ratio);
    //!< evolve method evolves the integrator
    void evolve(Float dt, int steps);
    //!< evolve method relevant to split timescales
    void evolve(Float dt, int steps, int **fractionSplit);

    //!< Compare two approximations to avoid recalculation if possible
    static int compareApprox(RemezArg &, RemezArg &);

    void init(int traj_num);

    void setSplit(int i, int j);
    //!< Set mass i, pole j as included
    void checkSplit();
    //!< Check that all the partial fractions have been accounted for
  
    const ActionRationalArg & getRationalArg() const{ return *rat_arg; }
};

/*!< 
  Derived from AlgActionRational, this class allows independent
  evolution of partial fractions.that appear in the rational action.
*/
class AlgActionRationalSplit : public AlgActionRational {

private:
    const char *cname;
    AlgActionRational *rat;
    ActionRationalSplitArg *rat_split_arg;
    int **fractionSplit;

public:
    AlgActionRationalSplit(AlgActionRational &rat, 
                           ActionRationalSplitArg &rat_split_arg); 
			 
    virtual ~AlgActionRationalSplit();

    void heatbath();
    Float energy();

    // Calculate preliminary force for force gradient
    void prepare_fg(Matrix * force, Float dt_ratio);
    //!< evolve method evolves the integrator
    void evolve(Float dt, int steps);
  
    int getNmass();
    Float getMass(int);

    void cost(CgStats*);
};

/*!< 
  Class describing bosonic bilinear action.  Action is given by
  
  S = phi^dag (M_b^dag M_b) phi
*/
class AlgActionBoson : public AlgActionBilinear {

private:
    const char *cname;

    ActionBosonArg *bsn_arg;
    CgArg **bsn_cg_arg;   //!< Pointer to an array of solver parameters.

public:

    AlgActionBoson(AlgMomentum &mom, ActionBosonArg &boson_arg);
    virtual ~AlgActionBoson();
  
    void heatbath();

    Float energy();

    // Calculate preliminary force for force gradient
    virtual void prepare_fg(Matrix * force, Float dt_ratio);

    //!< evolve method evolves the integrator
    void evolve(Float dt, int steps);
};

/*!< 
  Class describing fermionic bilinear action

  S = phi^dag (M_f^dag M_f)^{-1} phi

  Note that a chronological inverter can be used to accelerate the
  inversion process.
*/
class AlgActionFermion : public AlgActionBilinear {

private:
    const char *cname;

    ActionFermionArg *frm_arg;

    CgArg **frm_cg_arg_md;   //!< Pointer to an array of solver parameters.
    CgArg **frm_cg_arg_fg;   //!< Pointer to an array of solver parameters for force gradient step
    CgArg **frm_cg_arg_mc;   //!< Pointer to an array of solver parameters.

    int evolved;
//    Float h_init;

    //!< Stores the history of cg solutions - used by chronological inversion
    Vector ***v;
    Vector ***cg_sol_old;
    Vector *cg_sol;

    //!< Stores the orthogonalised vectors multiplied by MatPcDagMatPc
    //!< These currently live in AlgActionFermion for a future
    //!< optimisation (chronological preconditioner)
    Vector ***vm;

    int *chrono;

    // !< Status variable controls if we can use the CG solution from a
    // !< previous force gradient solve to forecast the next normal solve.
    bool fg_forecast;

public:

    AlgActionFermion(AlgMomentum &mom, ActionFermionArg &frm_arg);
    virtual ~AlgActionFermion();
  
    void heatbath();

    Float energy();

    // Calculate preliminary force for force gradient
    virtual void prepare_fg(Matrix * force, Float dt_ratio);
    //!< evolve method evolves the integrator
    void evolve(Float dt, int steps);

    void init();
};

/*!< 
  This class is designed for simulating two flavours of domain wall
  fermion, where the action we wish to include is the given by both a
  bosonic and fermionic determinant.  The action is given by

  S = phi^dag M_b (M_f^dag M_f)^{-1} M_b^dag phi

  This class can also be used for the Hasenbusch trick for other
  fermion types.
*/
class AlgActionQuotient : public AlgActionBilinear {

private:
    const char *cname;

    ActionQuotientArg *quo_arg;

    std::vector<CgArg> bsn_cg_arg;   //!< Pointer to an array of solver parameters.

    //!< Pointer to an array of solver parameters, for force gradient
    //!< step, irrevelant if using other integrators.
    std::vector<CgArg> frm_cg_arg_fg;

    std::vector<CgArg> frm_cg_arg_md;   //!< Pointer to an array of solver parameters.
    std::vector<CgArg> frm_cg_arg_mc;   //!< Pointer to an array of solver parameters.

    std::vector<Float> bsn_mass; //!< The boson mass parameter that appears in the quotient
    std::vector<Float> frm_mass; //!< The fermion mass parameter that appears in the quotient

    // ~~added for twisted mass Wilson fermions
    std::vector<Float> bsn_mass_epsilon; //!< The boson mass parameter that appears in the quotient
    std::vector<Float> frm_mass_epsilon; //!< The fermion mass parameter that appears in the quotient

    int evolved;
//    Float h_init;

    //!< Stores the history of cg solutions - used by chronological inversion
    Vector ***v;
    Vector ***cg_sol_old;
    Vector *cg_sol;
    Vector *tmp1;
    Vector *tmp2;
    Vector *tmp3;// temporarily keep RGV for debugging

    //!< Stores the orthogonalised vectors multiplied by MatPcDagMatPc
    //!< These currently live in AlgActionQuotient for a future
    //!< optimisation (chronological preconditioner)
    Vector ***vm;

    std::vector<int> chrono;

    // !< Status variable controls if we can use the CG solution from a
    // !< previous force gradient solve to forecast the next normal solve.
    bool fg_forecast;
public:

    AlgActionQuotient(AlgMomentum &mom, ActionQuotientArg &frm_arg);
    virtual ~AlgActionQuotient();
  
    void reweight(Float *rw_fac,Float *norm);
    void heatbath();

    Float energy();

    // Calculate preliminary force for force gradient
    virtual void prepare_fg(Matrix * force, Float dt_ratio);
    //!< evolve method evolves the integrator
    void evolve(Float dt, int steps);

    void init();

    void set_skip_force(bool _skip_force) { skip_force = _skip_force; }
};

/*!< 
  This class is designed for simulating a single flavour of domain
  wall fermion, the action we wish to include is the squareroot of
  both a bosonic and fermionic determinant.  The action is given by

  S = phi^dag r_b r_f^2 r_b phi

  where r_b is the bosonic rational function and r_f is the fermionic
  rational function.
*/
class AlgActionRationalQuotient : public AlgActionRational {

private:
    const char *cname;
    int **fractionSplit;
    ActionRationalQuotientArg *rat_quo_arg;


    //protected:
public:
    //CK: temporarily public, for testing

    Float *bsn_mass;  //!< The boson mass parameter that appears in the quotient
    Float *frm_mass;  //!< The fermion mass parameter that appears in the quotient

    //CK: ~~added for twisted mass Wilson fermions (for DSDR term with G-parity square-root)
    Float* bsn_mass_epsilon; //!< The boson mass parameter that appears in the quotient
    Float* frm_mass_epsilon; //!< The fermion mass parameter that appears in the quotient

    //!< This is where the rational parameters are stored
    RemezArg *frm_remez_arg_md;
    RemezArg *frm_remez_arg_mc;
    RemezArg *bsn_remez_arg_md;
    RemezArg *bsn_remez_arg_mc;

    //!<  bsn_cg_arg_fg is specific to force gradient integrator and has no use otherwise.
    CgArg ***bsn_cg_arg_fg;
    CgArg ***bsn_cg_arg_md;
    CgArg ***bsn_cg_arg_mc;
    //!< Pointer to an array of structures containing solver parameters.
    /*!<
      These are the parameters corresponding to each of the dynamical fermion
      masses.      
    */

    Vector **eta; //!< Use to accumulate solver results

public:

    AlgActionRationalQuotient();
    AlgActionRationalQuotient(AlgMomentum &mom, 
                              ActionRationalQuotientArg &rat_quo_arg, int traj_num=0);
    virtual ~AlgActionRationalQuotient();

    void reweight(Float *rw_fac,Float *norm);
    void heatbath();
    Float energy();

    // Calculate preliminary force for force gradient
    virtual void prepare_fg(Matrix * force, Float dt_ratio);
    //!< evolve method evolves the integrator
    void evolve(Float dt, int steps);
    //!< evolve method relevant to split timescales
    void evolve(Float dt, int steps, int **fractionSplit);

    bool loadPoles(void);
    bool savePoles(void);
    bool checkPolesFile(const RemezArg &md, const RemezArg &mc, const RationalDescr &r);

    void set_skip_force(bool _skip_force) { skip_force = _skip_force; }
};

// Exact One Flavor Algorithm (EOFA)
// Describes a single flavor of Domain Wall fermion as a
// quotient det(D(mf))/det(D(mb)).
class AlgActionEOFA : public AlgActionBilinear
{
  private:
    const char* cname;
    ActionEOFAArg* eofa_arg;

  protected:
    bool fg_forecast;
    bool heatbath_forecast;
    bool heatbath_test;
    int evolved;
    int heatbathEval;
    int energyEval;
    Float h_init;
    Float* num_mass;
    Float* den_mass;
    RemezArg* LH_remez_arg;
    RemezArg* RH_remez_arg;
    CgArg** LH_cg_arg_fg;
    CgArg** LH_cg_arg_md;
    CgArg** LH_cg_arg_mc;
    CgArg*** LH_cg_arg_heatbath;
    CgArg** RH_cg_arg_fg;
    CgArg** RH_cg_arg_md;
    CgArg** RH_cg_arg_mc;
    CgArg*** RH_cg_arg_heatbath;
    Vector** frmn_tmp;
    EigArg eig_arg;
    char eig_file[256];
    CommonArg ca_eig;
    Float** lambda_low;
    Float** lambda_high;

    void generateApprox(Float* mass, RemezArg** remez_arg, EOFARationalDescr* rat);
    void destroyApprox(RemezArg* remez_arg);
    static int compareApprox(RemezArg& arg1, RemezArg& arg2);
    void generateCgArg(Float* mass, CgArg**** cg_arg, const char* l, EOFARationalDescr* rat);
    void generateCgArg(Float* mass, CgArg*** cg_arg_fg, CgArg*** cg_arg_mc, CgArg*** cg_arg_md,
        const char* l, Float* stop_rsd_fg, Float* stop_rsd_mc, Float* stop_rsd_md);
    void destroyCgArg(CgArg*** cg_arg, const char* l, RemezArg* remez_arg);
    void destroyCgArg(CgArg** cg_arg_fg, CgArg** cg_arg_mc, CgArg** cg_arg_md, const char* l);
    void generateEigArg(const EigenDescr& eigen);
    void destroyEigArg();
    void checkApprox(Float* mass, RemezArg* remez_arg, EigenDescr& eigen);

  public:
    AlgActionEOFA();
    AlgActionEOFA(AlgMomentum& mom, ActionEOFAArg& arg, bool _heatbath_forecast=true, 
        bool _heatbath_test=true, int traj_num=0);
    std::vector<double> eig_range(Float m_num, Float m_den);
    double ritz(Float m_num, Float m_den, bool compute_min);
    void Meofa(Vector* out, const Vector* in, Float m_num, Float m_den, Lattice& lat);
    void reweight(Float* rw_fac, Float* norm);
    void heatbath();
    Float energy();
    virtual void prepare_fg(Matrix* force, Float dt_ratio);
    void evolve(Float dt, int steps);
    void init(int traj_num);
    bool loadPoles();
    bool savePoles();
    bool checkPolesFile(const RemezArg& ra, const EOFARationalDescr& r);
    void set_skip_force(bool _skip_force) { skip_force = _skip_force; }
    virtual ~AlgActionEOFA();
};

/*!< 
  Class describing the pure gauge action contribution to the
  Hamiltonian.
*/
class AlgActionGauge : public AlgAction {

private:
    const char *cname;
    ActionGaugeArg *gauge_arg;
    GclassType gluon;

public:

    AlgActionGauge(AlgMomentum &mom, ActionGaugeArg &gauge_arg);
    virtual ~AlgActionGauge();

    void heatbath();

    Float energy();

    //!< preparation for force gradient evolution
    virtual void prepare_fg(Matrix * force, Float dt_ratio);
    //!< evolve method evolves the integrator
    void evolve(Float dt, int steps);

    void cost(CgStats*);

    void init();

  //!< when G-parity active, this function (which should be do something only by AlgActionGauge) copy-conjugates the gauge field
  void copyConjLattice();
};


CPS_END_NAMESPACE

#endif
