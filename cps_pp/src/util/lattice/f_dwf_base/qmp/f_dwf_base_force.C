#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of FdwfBase class.

*/
//------------------------------------------------------------------
//
// f_dwf_base.C
//
// FdwfBase is derived from FwilsonTypes and is relevant to
// domain wall fermions, moved from old Fdwf class
//
//------------------------------------------------------------------

#if  0
#define OMP_DEBUG(msg,a...) do \
if (called==1) printf("[thread=%05d] %s;%d:QMP/%s(): " msg, omp_get_thread_num(),\
__FILE__,__LINE__,__FUNCTION__,##a); \
while(0);
#else
#define OMP_DEBUG(msg,a...)
#endif

CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <math.h>
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/dwf.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/vector.h>
#include <util/random.h>
#include <util/error.h>
#include <util/time_cps.h>
#include <util/omp_wrapper.h>
#include <comms/scu.h> // GRF
#include <comms/glb.h>
CPS_START_NAMESPACE

#define Printf if(!UniqueID()) printf

static inline int offset (int *size, int *pos_p, int mu = -1){
  if (mu>3) printf("FdwfBase::offset: Error!\n");
  int pos[4];
  for(int i=0;i<4;i++){
	 pos[i]=pos_p[i];
	if (pos[i]<0 || pos[i]>=size[i]) {
	printf("pos size %d= %d %d!!\n",i,pos[i],size[i]);exit(-11);
	}
  }
  if (mu>-1) pos[mu] = (pos[mu]+1)%size[mu];
  int result = 
  pos[0] + size[0] *( pos[1] + size[1] *( pos[2] + size[2] *( pos[3])));
//  if (mu>-1) pos[mu] = (pos[mu]-1+size[mu])%size[mu];
  return result;
}

#define PROFILE
// CJ: change start
//------------------------------------------------------------------
// EvolveMomFforce(Matrix *mom, Vector *chi, Float mass, 
//                 Float dt):
// It evolves the canonical momentum mom by dt
// using the fermion force.
//------------------------------------------------------------------
ForceArg FdwfBase::EvolveMomFforce(Matrix *mom, Vector *chi, 
			   Float mass, Float dt){
  const char *fname = "EvolveMomFforce(M*,V*,F,F) [QMP version]";
  VRB.Func(cname,fname);

  Matrix *gauge = GaugeField() ;

  if (Colors() != 3)
    ERR.General(cname,fname,"Wrong nbr of colors.") ;
 
  if (SpinComponents() != 4)
    ERR.General(cname,fname,"Wrong nbr of spin comp.") ;

  if(GJP.Snodes() > 1) 
    ERR.General(cname,fname,"Not written for Snodes()>1") ;
 
  if (mom == 0)
    ERR.Pointer(cname,fname,"mom") ;
 
  if (chi == 0)
    ERR.Pointer(cname,fname,"chi") ;
 
  //----------------------------------------------------------------
  // allocate space for two CANONICAL fermion fields
  //----------------------------------------------------------------

  size_t f_size = FsiteSize() * GJP.VolNodeSites() ;
  int f_site_size_4d = 2 * Colors() * SpinComponents();
  size_t f_size_4d = f_site_size_4d * GJP.VolNodeSites() ;
 
  size_t f_size_alloc = f_size *sizeof(Float);
  //CK: Need space for both d and C\bar u^T fields stacked
  //Two 4d volumes are stacked on each Ls such that Ls can be split over multiple nodes
  if(GJP.Gparity()) f_size_alloc *=2;

  const char *str_v1 = "v1" ;
  Vector *v1 = (Vector *)fmalloc(cname,fname,str_v1,f_size_alloc);

  const char *str_v2 = "v2" ;
  Vector *v2 = (Vector *)fmalloc(cname,fname,str_v2,f_size_alloc) ;

  //----------------------------------------------------------------
  // Calculate v1, v2. Both v1, v2 must be in CANONICAL order after
  // the calculation.
  //----------------------------------------------------------------  

  VRB.Clock(cname, fname, "Before calc force vecs.\n") ;

#ifdef PROFILE
  Float time = -dclock();
  Float time0=time;
  ForceFlops=0;
#endif
  {
    CgArg cg_arg ;
    cg_arg.mass = mass ;
#ifdef PROFILE
  time = -dclock();
#endif
    DiracOpDwf dwf(*this, v1, v2, &cg_arg, CNV_FRM_NO) ;
    dwf.CalcHmdForceVecs(chi) ;
#ifdef PROFILE
  time += dclock();
  print_flops(fname,"CalcHmdForceVecs()",0,time);
  time = -dclock();
#endif

    Fconvert(v1,CANONICAL,WILSON);
    Fconvert(v2,CANONICAL,WILSON);
#ifdef PROFILE
  time += dclock();
  print_flops(fname,"Fconvert()",0,time);
#endif
  }
  ForceArg farg = EvolveMomFforceInt(mom,v1,v2,mass,dt);

  ffree(cname, fname, str_v2, v2) ;
  ffree(cname, fname, str_v1, v1) ;
  return farg;

}

//#ifdef USE_TEST
#if 1
ForceArg FdwfBase::EvolveMomFforce( Matrix* mom, // momenta
                               Vector* phi, // odd pseudo-fermion field
                               Vector* eta, // very odd pseudo-fermion field
                               Float  mass, 
                               Float dt )
{
  const char *fname = "EvolveMomFforce(M*,V*,V*,F,F)";
  if(GJP.Gparity()) ERR.General(cname,fname,"Not implemented for G-parity\n");

  VRB.Func(cname,fname);
  
  Matrix *gauge = GaugeField() ;

  if (Colors() != 3)       { ERR.General(cname,fname,"Wrong nbr of colors.") ; }
  if (SpinComponents()!=4) { ERR.General(cname,fname,"Wrong nbr of spin comp.") ;}
  if (mom == 0)            { ERR.Pointer(cname,fname,"mom") ; }
  if (phi == 0)            { ERR.Pointer(cname,fname,"phi") ; }
   
  // allocate space for two CANONICAL fermion fields

  // these are all full fermion vector sizes ( i.e. *not* preconditioned )

  const size_t f_size        ( FsiteSize() * GJP.VolNodeSites() );
  const size_t f_size_cb     ( f_size/2 ) ; // f_size must be multiple of 2
  const int f_site_size_4d( 2 * Colors() * SpinComponents() );
  const size_t f_size_4d     ( f_site_size_4d * GJP.VolNodeSites()) ;
  
  char *str_v1 = "v1" ;
  Vector *v1 = (Vector *)fmalloc(cname,fname,str_v1,f_size*sizeof(Float));

  char *str_v2 = "v2" ;
  Vector *v2 = (Vector *)fmalloc(cname,fname,str_v2,f_size*sizeof(Float));

  Float L1 = 0.0;
  Float L2 = 0.0;
  Float Linf = 0.0;

#ifdef PROFILE
  Float time = -dclock();
  ForceFlops=0;
#endif

  //CK: Bosonic force vectors
  //f_out = v1 = (eta, Deo eta)
  //f_in = v2 = (-kappa^2 phi, -kappa^2 Deo^dag phi)

  
  //Note the fermionic force vectors calculated by CalcHmdForceVecs are the following:
  //f_out = ( chi, D_eo chi )
  //f_in = ( -kappa^2 Mprec chi, -kappa^2 D_eo^dag Mprec chi )    
  //where chi is the input vector. These are the same as the bosonic force vectors if we use  phi = Mprec chi  and  eta = chi as inputs


  //Calculate v1, v2. Both must be in CANONICAL order afterwards
  {
    CgArg cg_arg ;
    cg_arg.mass = mass ;
    
    DiracOpDwf dwf(*this, v1, v2, &cg_arg, CNV_FRM_YES) ;
    Float kappa( 1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight())));

    v2->CopyVec(phi,f_size_cb);

    // rescale the input field. As the second half of the this field
    // will be constructed by acting with the PC dslash on v1, this
    // rescales *one* of the full vectors - giving rise to an overall
    // rescaling of the final answer by exactly -\kappa^2
    
    v2->VecTimesEquFloat(-kappa*kappa,f_size_cb);

    // only need one factor of -\kappa^2, so don't rescale the second
    // full vector (v2)
    v1->CopyVec(eta,f_size_cb);
        
    dwf.Dslash(v2+(f_size_cb/6), v2 , CHKB_ODD, DAG_YES); //note 6 is the size of a Vector object!
    dwf.Dslash(v1+(f_size_cb/6), v1 , CHKB_ODD, DAG_NO);
    
    // v1 and v2 are now the vectors needed to contruct the force term
    // written in ( ODD, EVEN ) ordering. They will be converted back
    // into canonical ordering when the destructor is called.
    
  }  
#ifdef PROFILE
  time += dclock();
  print_flops(fname,"Dslash()",0,time);
#endif
// There is a relative minus sign compared to EvolveMomFforce(M*,V*,F,F)

  ForceArg farg = EvolveMomFforceInt(mom,v1,v2,mass,-dt);

  ffree(cname, fname, str_v2, v2) ;
  ffree(cname, fname, str_v1, v1) ;
  return farg;

}
#endif

ForceArg FdwfBase::EvolveMomFforceInt(Matrix *mom, Vector *v1, Vector *v2, 
			   Float mass, Float dt){
  const char *fname = "EvolveMomFforceInt(M*,V*,F,F) [QMP version]";
  VRB.Func(cname,fname);
  Matrix *gauge = GaugeField() ;
  Float time=0.;

  static int called=0;
  called++;
  #if TARGET == BGQ

#else
  const int MAX_THREADS=1;
#endif

  size_t f_size = FsiteSize() * GJP.VolNodeSites() ;
  int f_site_size_4d = 2 * Colors() * SpinComponents();
  size_t f_size_4d = f_site_size_4d * GJP.VolNodeSites() ;
 
  size_t f_size_alloc = f_size *sizeof(Float);
  int f_single4dsite_alloc = FsiteSize()*sizeof(Float);
  //CK: Need space for both d and C\bar u^T fields stacked
  //Two 4d volumes are stacked on each Ls such that Ls can be split over multiple nodes
  if(GJP.Gparity()){
    f_size_alloc *=2;
    f_single4dsite_alloc*=2;
  }

  int mu, x, y, z, t, s, lx, ly, lz, lt, ls ;
  int size[5],surf[4],blklen[4],stride[4],numblk[4];
 
  size[0] = lx = GJP.XnodeSites() ;
  size[1] = ly = GJP.YnodeSites() ;
  size[2] = lz = GJP.ZnodeSites() ;
  size[3] = lt = GJP.TnodeSites() ;
  size[4] = ls = GJP.SnodeSites() ;


  blklen[0] = sizeof(Float)*FsiteSize()/size[4]; //byte size of a spin-colour vector = byte length of block of data between x and x+1
  numblk[0] = (GJP.VolNodeSites()/size[0])*size[4]; //number of spin-colour vectors in one 4d slice of the 5d volume along the x-direction
  if(GJP.Gparity()) numblk[0]*=2; //2 flavours

  stride[0] = blklen[0] * (size[0]-1); //byte offset from end of spin-colour vector to start of next vector at same x
  for (int i =1;i<4;i++){
    blklen[i] = blklen[i-1] * size[i-1]; //byte length of block of data between pos[i] and pos[i]+1
    numblk[i] = numblk[i-1] / size[i]; //number of spin-colour vectors in the (4-i)d slice of the (4-i-1) slice owned by the previous direction 
    stride[i] = blklen[i] * (size[i]-1); //byte offset from end of spin-colour vector to start of next vector at same pos[i]
  }
  long vol = GJP.VolNodeSites(); //4d volume
  for (int i =0;i<4;i++)
  if ( GJP.Nodes(i) >1){
    surf[i] = GJP.VolNodeSites()/size[i]; //3d volume at surface in i-direction
  } else surf[i]=0;

  //----------------------------------------------------------------
  // allocate buffer space for two fermion fields that are assoc
  // with only one 4-D site.
  //----------------------------------------------------------------

  const char *str_site_v1 = "site_v1" ;
  const char *str_site_v2 = "site_v2" ;

  Float *v1_buf[4],*v2_buf[4];
  Float *v1_buf_pos[4];
  Float *v2_buf_pos[4];
  for(int i =0;i<4;i++) 
    if ( GJP.Nodes(i) >1){
      v1_buf[i]=(Float *)fmalloc(cname,fname,"v1_buf",surf[i]*f_single4dsite_alloc) ;
      v2_buf[i]=(Float *)fmalloc(cname,fname,"v2_buf",surf[i]*f_single4dsite_alloc) ;
    } else {
      v1_buf[i]=v2_buf[i]=NULL;
    }
  
  Matrix *gparity_1f_mombuf;
  if(GJP.Gparity1fX()){
    gparity_1f_mombuf = (Matrix *)fmalloc(cname,fname,"gparity_1f_mombuf",4 * GJP.VolNodeSites() * sizeof(Matrix) ) ;
    for(int i=0;i<4*GJP.VolNodeSites();i++) gparity_1f_mombuf[i].ZeroMatrix();
  }

  Matrix *tmp_mat = (Matrix *)fmalloc(cname,fname,"tmp_mat",sizeof(Matrix)*3*MAX_THREADS);

  size_t f_bytes = sizeof(Float)*f_site_size_4d; //byte size of spin-colour vector
  size_t st_bytes = sizeof(Float)*f_size_4d - f_bytes; //byte offset between s-slices
  if((ls*f_bytes)!=FsiteSize()*sizeof(Float))
    ERR.General(cname,fname,"ls(%d)*f_bytes(%d)!=FsiteSize()(%d)*sizeof(Float)(%d)\n",
    ls,f_bytes,FsiteSize(),sizeof(Float));

 QMP_msgmem_t Send_mem[4];
 QMP_msgmem_t Recv_mem[4];
 QMP_msghandle_t Send[4];
 QMP_msghandle_t Recv[4];
 void *addr[2];
 size_t blksize[2];
 int nblocks[2];
 ptrdiff_t strds[2];
 for (mu=0; mu<4; mu++){
   if ( GJP.Nodes(mu) >1){
     //declare the receive buffer with size = 4d surface
     addr[0]=v1_buf[mu]; addr[1]=v2_buf[mu];
     blksize[0]=blksize[1]=(size_t)(surf[mu]*ls*f_bytes);
     if(GJP.Gparity()){ blksize[0]*=2; blksize[1]*=2; } //2 flavours

     nblocks[0]=nblocks[1]=1;
     strds[0]=strds[1]=0;
     Recv_mem[mu]=QMP_declare_strided_array_msgmem(addr,blksize,nblocks,strds,2);
     Recv[mu]=QMP_declare_receive_relative(Recv_mem[mu],mu,1,0);
     
     //declare the send buffer from the v1 and v2 vectors
     addr[0]=v1; addr[1]=v2;
     blksize[0]=blksize[1]=(size_t)blklen[mu];
     nblocks[0]=nblocks[1]=numblk[mu]; //this is twice as large as usual for G-parity. Flavour index is nested within s-index.

     strds[0]=strds[1]=(ptrdiff_t)(stride[mu]+blklen[mu]); //for direction mu this is the separation between a site and the next with the same coordinate in the mu direction
     if((size_t)( (!GJP.Gparity() && surf[mu]*ls*f_bytes != blklen[mu]*numblk[mu]) || (GJP.Gparity() && surf[mu]*ls*f_bytes*2 != blklen[mu]*numblk[mu]) ) ){
       ERR.General(cname,fname,"mu=%d: receiving bytes(%d*%d*%d) does not match with sending bytes(%d*%d)\n",mu, surf[mu],ls,f_bytes, blklen[mu],numblk[mu]);
     }
     Send_mem[mu]=QMP_declare_strided_array_msgmem(addr,blksize,nblocks,strds,2);
     Send[mu]=QMP_declare_send_relative(Send_mem[mu],mu,-1,0);
   }
 }
 VRB.Clock(cname, fname, "Before loop over links.\n") ;
  
#ifdef PROFILE
  time = -dclock();
  ForceFlops=0;
#endif

  for (mu=0; mu<4; mu++){
    if ( GJP.Nodes(mu) >1){
        QMP_start(Recv[mu]);
        QMP_start(Send[mu]);
    }
  }
#if 1
  Float L1[MAX_THREADS] ;
  Float L2[MAX_THREADS];
  Float Linf[MAX_THREADS] ;
  for(int i =0;i<MAX_THREADS;i++){
    L1[i]=L2[i]=Linf[i]=0.;
  }
//reduction(+:L1,L2)
   	long i =0;
	//split work of the 4d volume * number of directions
#pragma omp parallel for default(shared) private(mu)
    for (i=0;i<vol*4;i++){
  	int pos[4];
	long rest=i;
	mu = rest%4; rest = rest/4;
	for(int j =0; j<4;j++){
		pos[j]= rest%size[j]; rest = rest/size[j];
	}
  	if((i==0) && (called%10000)==1) Printf("omp_get_num_threads=%d",omp_get_num_threads());
	int tnum = omp_get_thread_num();
        Matrix *tmp_mat1,*tmp_mat2, *tmp_mat3;
  	tmp_mat1 = tmp_mat + 3*omp_get_thread_num();
  	tmp_mat2 = tmp_mat1+1;
	tmp_mat3 = tmp_mat2+1;
#else
  int pos[4];
  for (mu=0; mu<4; mu++)
    for (pos[3]=0; pos[3]<size[3]; pos[3]++)
    for (pos[2]=0; pos[2]<size[2]; pos[2]++)
    for (pos[1]=0; pos[1]<size[1]; pos[1]++)
    for (pos[0]=0; pos[0]<size[0]; pos[0]++){
      Matrix *tmp_mat1,*tmp_mat2, *tmp_mat3;
	  tmp_mat1 = tmp_mat;
	  tmp_mat2 = tmp_mat1+1;
	  tmp_mat3 = tmp_mat2+1;
#endif
      int gauge_offset = offset(size,pos);
      int vec_offset = f_site_size_4d*gauge_offset ;
      gauge_offset = mu+4*gauge_offset ;

      Float *v1_plus_mu ;
      Float *v2_plus_mu ;

      //for G-parity
      Float *v1_plus_mu_f2 ;
      Float *v2_plus_mu_f2 ;

      int vec_plus_mu_stride ;
      int vec_plus_mu_offset = f_site_size_4d ;

      Float coeff = -2.0 * dt ;
      Float f2_coeff = coeff; //for G-parity

      vec_plus_mu_offset *= offset(size,pos,mu);
      
      if ((GJP.Nodes(mu)>1)&&((pos[mu]+1) == size[mu]) ) {
	//do nothing if on the boundary and >1 node in mu direction
      } else {
	v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
	v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;

	v1_plus_mu_f2 = (Float *)v1+vec_plus_mu_offset + f_size_4d;
	v2_plus_mu_f2 = (Float *)v2+vec_plus_mu_offset + f_size_4d;

	if (GJP.Bc(mu)==BND_CND_GPARITY && pos[mu]+1 == size[mu]){
	  //G-parity flavour twist
	  v1_plus_mu_f2 = (Float *)v1+vec_plus_mu_offset ;
	  v2_plus_mu_f2 = (Float *)v2+vec_plus_mu_offset ;

	  v1_plus_mu = (Float *)v1+vec_plus_mu_offset + f_size_4d;
	  v2_plus_mu = (Float *)v2+vec_plus_mu_offset + f_size_4d;
	}

	//s-stride
	if(GJP.Gparity()) vec_plus_mu_stride = 2*f_size_4d - f_site_size_4d ; //2 stacked fields per s-slice
	else vec_plus_mu_stride = f_size_4d - f_site_size_4d ; 

	int vec_stride = vec_plus_mu_stride; //strides are the same for local comms
	
	if ((pos[mu]+1) == size[mu]){
	  if (GJP.NodeBc(mu)==BND_CND_APRD){
	    coeff = -coeff ;
	    f2_coeff = -f2_coeff;
	  }else if(GJP.Bc(mu)==BND_CND_GPARITY){ //NOT NodeBc
	    //at global boundary (only 1 node in this direction)
	    // d <- +CubarT,  CubarT <- -d
	    f2_coeff = -f2_coeff;
	  }
	}

        OMP_DEBUG("gauge_offset=%d vec_plus_mu_offset=%d vec_plus_mu_stride=%d thread=%d\n",gauge_offset,vec_plus_mu_offset,vec_plus_mu_stride,omp_get_thread_num());
	//printf("Doing outer product for f1 P-\n"); fflush(stdout);

	//tmp_mat1_{ij} = tr_s ( P_-^mu v1(x+mu)_i v2^\dagger(x)_j )
	sproj_tr[mu]( (IFloat *)tmp_mat1,
		      (IFloat *)v1_plus_mu,
		      (IFloat *)v2+vec_offset,
		      ls, vec_plus_mu_stride, vec_stride) ;
	//printf("Doing outer product for f1 P+\n"); fflush(stdout);

	//tmp_mat2_{ij} = tr_s ( P_+^mu v2(x+mu)_i v1^\dagger(x)_j 
	sproj_tr[mu+4]( (IFloat *)tmp_mat2,
			(IFloat *)v2_plus_mu,
			(IFloat *)v1+vec_offset,
			ls, vec_plus_mu_stride, vec_stride) ;
        OMP_DEBUG("v1_plus_mu=%p v2_plus_mu=%p,thread=%d\n",v1_plus_mu,v2_plus_mu ,omp_get_thread_num());

	*tmp_mat1 += *tmp_mat2 ;
	if(GJP.Gparity()){
	  *tmp_mat1 *= coeff; //multiply by f1 coeff here as f2 coeff might be different

	  //for flavour 2 field we must flip the projection operator used for the force term
	  //and swap the coordinates of the force vectors

	  //printf("Doing outer product for f2 P+\n"); fflush(stdout);

	  //tmp_mat2_{ij} = tr_s ( P_+^mu v1(x)_i v2^\dagger(x+mu)_j 
	  sproj_tr[mu+4]( (IFloat *)tmp_mat2,
			  (IFloat *)v1+vec_offset+f_size_4d,
			  (IFloat *)v2_plus_mu_f2,
			  ls, vec_stride, vec_plus_mu_stride) ;

	  //printf("Doing outer product for f2 P-\n"); fflush(stdout);
	  //tmp_mat3_{ij} = tr_s ( P_-^mu v2(x)_i v1^\dagger(x+mu)_j 
	  sproj_tr[mu]( (IFloat *)tmp_mat3,
			(IFloat *)v2+vec_offset+f_size_4d,
			(IFloat *)v1_plus_mu_f2,
			ls, vec_stride, vec_plus_mu_stride) ;

	  *tmp_mat2 += *tmp_mat3;

	  tmp_mat3->Trans(*tmp_mat2);//must transpose outer product on colour index
	  *tmp_mat3 *= f2_coeff;//different coefficient for the two flavours
	  *tmp_mat1 += *tmp_mat3 ;
	}

      // If GJP.Snodes > 1 sum up contributions from all s nodes
//      if(GJP.Snodes() > 1) {
//	  glb_sum_multi_dir((Float *)tmp_mat1, 4, sizeof(Matrix)/sizeof(IFloat) ) ;
//      }

      tmp_mat2->DotMEqual(*(gauge+gauge_offset), *tmp_mat1) ;
//        OMP_DEBUG(" tmp_mat2->DotMEqual(*(gauge+gauge_offset), *tmp_mat1) ; thread=%d\n",omp_get_thread_num());

      tmp_mat1->Dagger(*tmp_mat2) ;

      tmp_mat2->TrLessAntiHermMatrix(*tmp_mat1) ;
      OMP_DEBUG(" tmp_mat2->TrLessAntiHermMatrix(*tmp_mat1) ; thread=%d\n",omp_get_thread_num());

      if(!GJP.Gparity()) *tmp_mat2 *= coeff ;

      *(mom+gauge_offset) += *tmp_mat2 ;

      //G-parity 2f setup
      if(GJP.Gparity()){
	//printf("Copy conjugating\n"); fflush(stdout);
	//set force for U* link
	Matrix* mom_u = mom+gauge_offset;
	Matrix* mom_ustar = mom+gauge_offset+GJP.VolNodeSites()*4;
	mom_ustar->Conj((IFloat*)mom_u);
      }

      if(GJP.Gparity1fX()) *(gparity_1f_mombuf+gauge_offset) = *tmp_mat2;

      OMP_DEBUG("*(mom+gauge_offset) += *tmp_mat2 ; thread=%d\n",omp_get_thread_num());
#if 1
{
      Float norm = tmp_mat2->norm();
      Float tmp = sqrt(norm);
//#pragma omp atomic
      L1[tnum] += tmp;
      L2[tnum] += norm;
      Linf[tnum] = (tmp>Linf[tnum] ? tmp : Linf[tnum]);
}
#endif
    }

    } // end for x,y,z,t, mu
   OMP_DEBUG("Loop ended! %d\n",mu);
  uint64_t temp_flop = (2*9*16*ls + 18+ 198+36+24)*lx*ly*lz*lt*(4.-1./lx-1/ly-1./lz-1./lt);
  if(GJP.Gparity()) temp_flop*=2;
  ForceFlops += temp_flop;
#ifdef PROFILE
  time += dclock();
  print_flops(fname,"local",temp_flop,time);
  time =- dclock();
#endif


  for (mu=0; mu<4; mu++){
    if ( GJP.Nodes(mu) >1){
      QMP_wait(Send[mu]);
      QMP_wait(Recv[mu]);
    }
  }
#ifdef PROFILE
  time += dclock();
  print_flops(fname,"QMP_wait",0,time);
  time =- dclock();
#endif
  static int surf_initted=0;
  static int *surf_table[4];
  static int *surf_mu_offset[4];

  //record the G-parity status when the initialization was performed
  //if this changes then reinitialize
  static int gparity_init_status;
  if(surf_initted && ( gparity_init_status && !GJP.Gparity() || !gparity_init_status && GJP.Gparity() ) ) surf_initted = 0;

  if (!surf_initted){
    if(GJP.Gparity()) gparity_init_status = 1;
    else gparity_init_status = 0;

    for (mu=0; mu<4; mu++)
      if ( GJP.Nodes(mu) >1){
	surf_table[mu] = (int *)smalloc(cname,fname,"surf_table[mu]=",surf[mu]*sizeof(int));
	surf_mu_offset[mu] = (int *)smalloc(cname,fname,"surf_table[mu]=",surf[mu]*sizeof(int));
      } else {
	surf_table[mu]=surf_mu_offset[mu]=NULL;
      }
    int pos[4];
    int surf_index[4];
    for (mu=0; mu<4; mu++)
      surf_index[mu]=0;
    for (mu=0; mu<4; mu++)
      for (pos[3]=0; pos[3]<size[3]; pos[3]++)
	for (pos[2]=0; pos[2]<size[2]; pos[2]++)
	  for (pos[1]=0; pos[1]<size[1]; pos[1]++)
	    for (pos[0]=0; pos[0]<size[0]; pos[0]++)
	      if ((GJP.Nodes(mu)>1)&&((pos[mu]+1) == size[mu]) ) {
		int surf_offset=0;
		if (mu!=3) surf_offset =pos[3];
		for(int ii=2;ii>-1;ii--)
		  if(mu!=ii){
		    surf_offset = pos[ii] + size[ii]*surf_offset;	
		  }
		if (surf_offset != surf_index[mu])
		  ERR.General(cname,fname,"surf_offset(%d) != surf_index[%d](%d)\n",
			      surf_offset,mu,surf_index[mu]);
		int gauge_offset = offset(size,pos);
		*(surf_table[mu]+surf_offset) =gauge_offset;
		*(surf_mu_offset[mu]+surf_offset) =offset(size,pos,mu);
		surf_index[mu] ++;
	  
	      }
    surf_initted=1;
#ifdef PROFILE
    time += dclock();
    print_flops(fname,"setup",0,time);
    time =- dclock();
#endif
  }

//------------------------------------------------------------------
// start by summing first over direction (mu) and then over site
// to allow SCU transfers to happen face-by-face in the outermost
// loop.
//------------------------------------------------------------------

  for(int i = 0;i<4;i++){
    v1_buf_pos[i] = v1_buf[i];
    v2_buf_pos[i] = v2_buf[i];
  }

//#ifdef USE_TEST
#if 1
    for (mu=0; mu<4; mu++)
	if (surf[mu]>0){
#pragma omp parallel for
	  for (long i=0;i<surf[mu];i++){ //split surface volume over threads
	    Matrix *tmp_mat1,*tmp_mat2,*tmp_mat3;
	    int tnum = omp_get_thread_num();
	    tmp_mat1 = tmp_mat + 3*tnum;
	    tmp_mat2 = tmp_mat1+1;
	    tmp_mat3 = tmp_mat2+1;
	    int gauge_offset = *(surf_table[mu]+i);
	    int vec_offset = f_site_size_4d*gauge_offset ;
	    gauge_offset = mu+4*gauge_offset ;
	    Float *v1_plus_mu ;
	    Float *v2_plus_mu ;

	    //for G-parity
	    Float *v1_plus_mu_f2 ;
	    Float *v2_plus_mu_f2 ;
	
            v1_plus_mu = v1_buf[mu] + f_site_size_4d*i;
            v2_plus_mu = v2_buf[mu] + f_site_size_4d*i;

	    v1_plus_mu_f2 = v1_plus_mu + f_site_size_4d*surf[mu];
	    v2_plus_mu_f2 = v2_plus_mu + f_site_size_4d*surf[mu];

	    if(GJP.Bc(mu) == BND_CND_GPARITY && GJP.NodeCoor(mu) == GJP.Nodes(mu)-1){
	      //G-parity flavour twist
	      v1_plus_mu_f2 = v1_buf[mu] + f_site_size_4d*i;
	      v2_plus_mu_f2 = v2_buf[mu] + f_site_size_4d*i;

	      v1_plus_mu = v1_plus_mu_f2 + f_site_size_4d*surf[mu];
	      v2_plus_mu = v2_plus_mu_f2 + f_site_size_4d*surf[mu];
	    }

	    //jump from end of vector to next s
            int vec_plus_mu_stride = (surf[mu] -1)*f_site_size_4d ;
	    int vec_stride = f_size_4d-f_site_size_4d;
	    if(GJP.Gparity()){
	      vec_stride = 2*f_size_4d-f_site_size_4d;
	      vec_plus_mu_stride = (2*surf[mu] -1)*f_site_size_4d ;
	    }

	    //boundary coefficients
	    Float coeff = -2.0 * dt ;
	    Float f2_coeff = coeff;
	    
	    if (GJP.NodeBc(mu)==BND_CND_APRD){
	      coeff = -coeff ;
	      f2_coeff = -f2_coeff;
	    }else if(GJP.Bc(mu)==BND_CND_GPARITY && GJP.NodeCoor(mu) == GJP.Nodes(mu)-1){
	      // d <- +CubarT,  CubarT <- -d
	      f2_coeff = -f2_coeff;
	    }
	    //printf("Doing outer product for f1 P+\n"); fflush(stdout);

	    //tmp_mat1_{ij} = tr_s ( P_+^mu v1(x+mu)_i v2^\dagger(x)_j )
	    sproj_tr[mu]( (IFloat *)tmp_mat1,
			  (IFloat *)v1_plus_mu,
			  (IFloat *)v2+vec_offset,
			  ls, vec_plus_mu_stride, vec_stride) ;

	    //printf("Doing outer product for f1 P-\n"); fflush(stdout);

	    //tmp_mat2_{ij} = tr_s ( P_-^mu v2(x+mu)_i v1^\dagger(x)_j 
	    sproj_tr[mu+4]( (IFloat *)tmp_mat2,
			    (IFloat *)v2_plus_mu,
			    (IFloat *)v1+vec_offset,
			    ls, vec_plus_mu_stride, vec_stride) ;

	    *tmp_mat1 += *tmp_mat2 ;

	    if(GJP.Gparity()){
	      *tmp_mat1 *= coeff; //multiply by f1 coeff here as f2 coeff might be different

	      //for flavour 2 field we must flip the projection operator used for the force term
	      //and swap the coordinates of the force vectors

	      //printf("Doing outer product for f2 P-\n"); fflush(stdout);

	      //tmp_mat2_{ij} = tr_s ( P_-^mu v1(x)_i v2^\dagger(x+mu)_j 
	      sproj_tr[mu+4]( (IFloat *)tmp_mat2,
			      (IFloat *)v1+vec_offset+f_size_4d,
			      (IFloat *)v2_plus_mu_f2,
			      ls, vec_stride, vec_plus_mu_stride) ;

	      //printf("Doing outer product for f2 P+\n"); fflush(stdout);

	      //tmp_mat3_{ij} = tr_s ( P_+^mu v2(x)_i v1^\dagger(x+mu)_j 
	      sproj_tr[mu]( (IFloat *)tmp_mat3,
			    (IFloat *)v2+vec_offset+f_size_4d,
			    (IFloat *)v1_plus_mu_f2,
			    ls, vec_stride, vec_plus_mu_stride) ;

	      *tmp_mat2 += *tmp_mat3;

	      tmp_mat3->Trans(*tmp_mat2);//must transpose outer product on colour index
	      *tmp_mat3 *= f2_coeff;//different coefficient for the two flavours
	      *tmp_mat1 += *tmp_mat3 ;
	    }
	    
	    // If GJP.Snodes > 1 sum up contributions from all s nodes
	    //      if(GJP.Snodes() >1 ) {
	    //	  glb_sum_multi_dir((Float *)tmp_mat1, 4, sizeof(Matrix)/sizeof(IFloat) ) ;
	    //      }

	    tmp_mat2->DotMEqual(*(gauge+gauge_offset), *tmp_mat1) ;

	    tmp_mat1->Dagger(*tmp_mat2) ;

	    tmp_mat2->TrLessAntiHermMatrix(*tmp_mat1) ;

	    if(!GJP.Gparity()) *tmp_mat2 *= coeff ;

	    *(mom+gauge_offset) += *tmp_mat2 ;

	    if(GJP.Gparity()){
	      //set force for U* link
	      Matrix* mom_u = mom+gauge_offset;
	      Matrix* mom_ustar = mom+gauge_offset+GJP.VolNodeSites()*4;
	      mom_ustar->Conj((IFloat*)mom_u);
	    }
	    
	    if(GJP.Gparity1fX()) *(gparity_1f_mombuf+gauge_offset) = *tmp_mat2;

	    Float norm = tmp_mat2->norm();
	    Float tmp = sqrt(norm);
	    //#pragma omp atomic
	    {
	      L1[tnum] += tmp;
	      L2[tnum] += norm;
	      Linf[tnum] = (tmp>Linf[tnum] ? tmp : Linf[tnum]);
	    }
	  }
	}// mu

#else
    if(GJP.Gparity() || Gparity1fX()) ERR.General(cname,fname," alternative code for surface sites not implemented for G-parity");
# if 1
// reduction(+:L1,L2)
//#pragma omp parallel for default(shared) private(mu)
    for (long i=0;i<vol*4;i++){
  	int pos[4];
	long rest=i;
	for(int j =0; j<4;j++){
		pos[j]= rest%size[j]; rest = rest/size[j];
	}
	mu = rest;
        Matrix *tmp_mat1,*tmp_mat2;
  	tmp_mat1 = tmp_mat + 2*omp_get_thread_num();
  	tmp_mat2 = tmp_mat1+1;
	int tnum = omp_get_thread_num();
        if((called%10000)==1) Printf("%d %d %d %d %d thread=%d\n",mu,pos[0],pos[1],pos[2],pos[3],tnum);
#else
    int pos[4];
    for (mu=0; mu<4; mu++)
    for (pos[3]=0; pos[3]<size[3]; pos[3]++)
    for (pos[2]=0; pos[2]<size[2]; pos[2]++)
    for (pos[1]=0; pos[1]<size[1]; pos[1]++)
    for (pos[0]=0; pos[0]<size[0]; pos[0]++){
  	Matrix *tmp_mat1,*tmp_mat2;
	  tmp_mat1 = tmp_mat;
	  tmp_mat2 = tmp_mat1+1;
	int tnum=0;
#endif
      int gauge_offset = offset(size,pos);
      int vec_offset = f_site_size_4d*gauge_offset ;
      gauge_offset = mu+4*gauge_offset ;

      Float *v1_plus_mu ;
      Float *v2_plus_mu ;
      int vec_plus_mu_stride ;
      int vec_plus_mu_offset = f_site_size_4d ;

      Float coeff = -2.0 * dt ;

//      Printf("%d %d %d %d %d\n",mu,pos[0],pos[1],pos[2],pos[3]);
          vec_plus_mu_offset *= offset(size,pos,mu);
      if ((GJP.Nodes(mu)>1)&&((pos[mu]+1) == size[mu]) ) {
		int surf_offset=0;
		if (mu!=3) surf_offset =pos[3];
		for(int ii=2;ii>-1;ii--)
		if(mu!=ii){
			surf_offset = pos[ii] + size[ii]*surf_offset;	
		}
	
            v1_plus_mu = v1_buf_pos[mu] ;
            v2_plus_mu = v2_buf_pos[mu] ;
			if ((v1_buf_pos[mu]-v1_buf[mu]) != (surf_offset*f_site_size_4d))
				ERR.General(cname,fname, "(v1_buf_pos[%d]-v1_buf[%d])(%d) != (surf_offset(%d)*f_site_size_4d(%d))(%d)\n",mu,mu,(v1_buf_pos[mu]-v1_buf[mu]),surf_offset,f_site_size_4d,(surf_offset*f_site_size_4d));
			if(v1_buf_pos[mu]!=(v1_buf[mu]+surf_offset*f_site_size_4d))
			ERR.General(cname,fname,"v1_buf_pos[%d](%p)!=(v1_buf[%d](%p)+surf_offset(%d)*f_site_size_4d(%d))\n",mu,v1_buf_pos[mu],mu,v1_buf[mu],surf_offset,f_size_4d);
            vec_plus_mu_stride = (surf[mu] -1)*f_site_size_4d ;
 			v1_buf_pos[mu] += f_site_size_4d;
 			v2_buf_pos[mu] += f_site_size_4d;
            if (GJP.NodeBc(mu)==BND_CND_APRD) coeff = -coeff ;

      sproj_tr[mu]( (IFloat *)tmp_mat1,
                    (IFloat *)v1_plus_mu,
                    (IFloat *)v2+vec_offset,
                    ls, vec_plus_mu_stride, f_size_4d-f_site_size_4d) ;

      sproj_tr[mu+4]( (IFloat *)tmp_mat2,
                      (IFloat *)v2_plus_mu,
                      (IFloat *)v1+vec_offset,
                      ls, vec_plus_mu_stride, f_size_4d-f_site_size_4d) ;

      *tmp_mat1 += *tmp_mat2 ;

      // If GJP.Snodes > 1 sum up contributions from all s nodes
      if(GJP.Snodes() >1 ) {
	  glb_sum_multi_dir((Float *)tmp_mat1, 4, sizeof(Matrix)/sizeof(IFloat) ) ;
      }

      tmp_mat2->DotMEqual(*(gauge+gauge_offset), *tmp_mat1) ;

      tmp_mat1->Dagger(*tmp_mat2) ;

      tmp_mat2->TrLessAntiHermMatrix(*tmp_mat1) ;

      *tmp_mat2 *= coeff ;

      *(mom+gauge_offset) += *tmp_mat2 ;
      Float norm = tmp_mat2->norm();
      Float tmp = sqrt(norm);
//#pragma omp atomic
{
      L1[tnum] += tmp;
      L2[tnum] += norm;
      Linf[tnum] = (tmp>Linf[tnum] ? tmp : Linf[tnum]);
}
    }

    } // end for x,y,z,t,mu

#endif

    //G-parity 1f add \delta p from other 'flavour'
    //Providing you set up the 1f prop source correctly (minus sign on UR quadrant) this code also works for
    //1f G-parity in both X and Y directions
    if(GJP.Gparity1fX()){
      int momsz = GsiteSize() * GJP.VolNodeSites();
      Matrix *buf2 = (Matrix *)fmalloc(cname,fname,"buf2",momsz * sizeof(Float) ) ;

      //Communicate \delta p from first half onto second half and vice versa
      Matrix *data_buf = gparity_1f_mombuf;
      Matrix *send_buf = data_buf;
      Matrix *recv_buf = buf2;

      if(GJP.Xnodes()>1){
	//pass between nodes
	for(int i=0;i<GJP.Xnodes()/2;i++){
	  getMinusData((Float *)recv_buf, (Float *)send_buf, momsz , 0);
	  data_buf = recv_buf;
	  recv_buf = send_buf;
	  send_buf = data_buf;
	}
      }else{
	//shift field by xsites/2
#ifdef USE_OMP
#pragma omp parallel for default(shared) private(mu)
	for(long i=0;i<vol*4;i++){
#else
	for(long i=0;i<vol*4;i++){
#endif
	  //i = mu + 4*(x + Lx*(y+Ly*(z+Lz*t) ) )
	  int mu = i%4;
	  int x = (i/4) % GJP.XnodeSites();
	  int pos_rem = i/4/GJP.XnodeSites(); //(y+Ly*(z+Lz*t)

	  int x_from = (x + GJP.XnodeSites()/2) % GJP.XnodeSites();
	  int i_from = mu + 4*(x_from + GJP.XnodeSites()*pos_rem);
	  buf2[i] = gparity_1f_mombuf[i_from];
	}
	data_buf = buf2;
      }



#ifdef USE_OMP
#pragma omp parallel for default(shared) private(mu)
    for (long i=0;i<vol*4;i++){
  	int pos[4];
	long rest=i;
	mu = rest%4; rest = rest/4;
	for(int j =0; j<4;j++){
		pos[j]= rest%size[j]; rest = rest/size[j];
	}
#else
  int pos[4];
  for (mu=0; mu<4; mu++)
    for (pos[3]=0; pos[3]<size[3]; pos[3]++)
    for (pos[2]=0; pos[2]<size[2]; pos[2]++)
    for (pos[1]=0; pos[1]<size[1]; pos[1]++)
    for (pos[0]=0; pos[0]<size[0]; pos[0]++){
#endif
      int gauge_offset = offset(size,pos);
      gauge_offset = mu+4*gauge_offset ;
      
      //complex conjugate the \delta p from the other flavour
      Float *m = (Float*) &data_buf[gauge_offset];
      for(int c=1;c<18;c+=2) m[c] *= -1;
      
      //add it to the momentum at this site
      *(mom+gauge_offset) += data_buf[gauge_offset];
    }
   }

  
  if(GJP.Gparity1fX()) ffree(gparity_1f_mombuf,cname,fname,"gparity_1f_mombuf");


  temp_flop = (2*9*16*ls + 18+ 198+36+24)*lx*ly*lz*lt*(1./lx+1/ly+1./lz+1./lt);
  ForceFlops += temp_flop;
#ifdef PROFILE
  time += dclock();
  print_flops(fname,"non-local",temp_flop,time);
#endif
 
//------------------------------------------------------------------
// deallocate smalloc'd space
//------------------------------------------------------------------

  for (mu=0; mu<4; mu++)
  if ( GJP.Nodes(mu) >1){
    QMP_free_msghandle(Send[mu]);
    QMP_free_msghandle(Recv[mu]);
    QMP_free_msgmem(Send_mem[mu]);
    QMP_free_msgmem(Recv_mem[mu]);
  }

  for(int i =0;i<4;i++) {
    ffree(v1_buf[i],cname,fname,"v1_buf");
    ffree(v2_buf[i],cname,fname,"v2_buf");
  }
 
  VRB.Sfree(cname, fname, "tmp_mat", tmp_mat) ;

  for(int i =1;i<MAX_THREADS;i++){
    L1[0] += L1[i];
    L2[0] += L2[i];
    Linf[0] = (Linf[i]>Linf[0] ? Linf[i] : Linf[0]);
  }
 
  glb_sum(&L1[0]);
  glb_sum(&L2[0]);
  glb_max(&Linf[0]);

  L1[0] /= 4.0*GJP.VolSites();
  L2[0] /= 4.0*GJP.VolSites();

  return ForceArg(L1[0], sqrt(L2[0]), Linf[0]);
}
// CJ: change end

CPS_END_NAMESPACE
