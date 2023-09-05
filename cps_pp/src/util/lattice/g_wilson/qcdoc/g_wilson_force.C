#include <config.h>
#include <math.h>
#include <util/gjp.h>
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pt.h>
#include <util/time_cps.h>
CPS_START_NAMESPACE
#define PROFILE
//------------------------------------------------------------------
// EvolveMomGforce(Matrix *mom, Float dt):
// It evolves the canonical momentum mom by dt
// using the pure gauge force.
//------------------------------------------------------------------
static const Float invs3 = -1./3.;
ForceArg Gwilson::EvolveMomGforce(Matrix *mom, Float dt){
  char *fname = "EvolveMomGforce(M*,F)";
  VRB.Func(cname,fname);
  static Matrix mt0;
  static Matrix *mp0 = &mt0;

  Float L1=0.0;
  Float L2=0.0;
  Float Linf=0.0;

#ifdef PROFILE
  Float time = -dclock();
  ForceFlops=0;
  ParTrans::PTflops=0;
#endif
  static int vol = GJP.VolNodeSites();
  const int N = 4;
  Float tmp = GJP.Beta() *invs3;
  Matrix *Unit = (Matrix *) fmalloc(vol*sizeof(Matrix));
  Matrix *tmp1[N];
  Matrix *tmp2[N];
  Matrix *result[4];
  for(int i = 0;i<4;i++){
  	result[i] = (Matrix *) fmalloc(vol*sizeof(Matrix));
  }
  for(int i = 0;i<N;i++){
  	tmp1[i] = (Matrix *) fmalloc(vol*sizeof(Matrix));
  	tmp2[i] = (Matrix *) fmalloc(vol*sizeof(Matrix));
        bzero((char *)tmp2[i],vol*sizeof(Matrix));
  }
  for(int i = 0;i<vol;i++) 
	Unit[i]=1.;
  Matrix *Units[4];
  for(int i = 0;i<N;i++) Units[i] = Unit;
  int mu,nu;
  {
    int dirs_p[] = {0,2,4,6,0,2,4};
    int dirs_m[] = {1,3,5,7,1,3,5};
    ParTransGauge pt(*this);

  
      for(nu = 1;nu<4;nu++){
        pt.run(N,tmp1,Units,dirs_m+nu);
  	pt.run(N,result,tmp1,dirs_m);
	pt.run(N,tmp1,result,dirs_p+nu);
	for(int i = 0; i<N;i++)
	vaxpy3_m(tmp2[i],&tmp,tmp1[i],tmp2[i],vol*3);
	pt.run(N,tmp1,Units,dirs_p+nu);
	pt.run(N,result,tmp1,dirs_m);
	pt.run(N,tmp1,result,dirs_m+nu);
	for(int i = 0; i<N;i++)
	vaxpy3_m(tmp2[i],&tmp,tmp1[i],tmp2[i],vol*3);
        ForceFlops +=vol*12*N;
      }
      pt.run(N,result,tmp2,dirs_p);
  }

  ForceFlops += vol*60;
  
#if 1
#pragma omp parallel for default(shared) private(mu) reduction(+:L1,L2) reduction(max:Linf)
  for (int index = 0; index < 4 * vol; index++) {
    Matrix mp1;
    int i = index % vol;
    mu = index / vol;
    Matrix *mtmp = (result[mu] + i);
    mp1.Dagger ((IFloat *) mtmp);
    mtmp->TrLessAntiHermMatrix (mp1);
    IFloat *ihp = (IFloat *) (mom + i * 4 + mu);	//The gauge momentum
//    IFloat *dotp = (IFloat *)mp0;
    IFloat *dotp2 = (IFloat *) (result[mu] + (i));
// assert (mtmp->norm () >= 0.);
    if (mtmp->norm () < 0.)
      printf ("%s::%s: Node %d: Gforce[%d][%d]=%0.12e\n", mu, i, mtmp->norm ());
    if (i < 4)
      VRB.Result (cname, fname, "Gforce[%d][%d]=%0.12e\n", mu, i, mtmp->norm ());
    fTimesV1PlusV2Single (ihp, dt, dotp2, ihp, 18);	//Update the gauge momentum
    Float norm = ((Matrix *) dotp2)->norm ();
    Float tmp = sqrt (norm);
    L1 += tmp;
    L2 += norm;
    Linf = (tmp>Linf ? tmp : Linf);
  }


#endif

  ForceFlops += vol*144;

#ifdef PROFILE
  time += dclock();
  print_flops(cname,fname,ForceFlops+ParTrans::PTflops,time);
#endif

  ffree(Unit);
  for(int i = 0;i<N;i++){
  ffree(tmp1[i]);
  ffree(tmp2[i]);
  }
  for(int i = 0;i<4;i++) 
  ffree(result[i]);

  glb_sum(&L1);
  glb_sum(&L2);
  glb_max(&Linf);

  L1 /= 4.0*GJP.VolSites();
  L2 /= 4.0*GJP.VolSites();

  return ForceArg(dt*L1, dt*sqrt(L2), dt*Linf);

}
CPS_END_NAMESPACE
