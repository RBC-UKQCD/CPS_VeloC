#include <config.h>
#include <math.h>
#include <util/ReadLatticePar.h>
#include <util/time_cps.h>
#include <assert.h>

CPS_START_NAMESPACE
using namespace std;

#ifdef USE_SERIAL_IO
int ReadLatticeParallel::UseParIO=0;
#else
int ReadLatticeParallel::UseParIO=1;
#endif

#define PROFILE

void ReadLatticeParallel::read(Lattice & lat, const QioArg & rd_arg)
{
  const char * fname = "read()";
  VRB.Func(cname,fname);
  
  char loginfo[strlen(rd_arg.FileName) + 10];
  sprintf(loginfo,"Load %s",rd_arg.FileName);
  startLogging(loginfo);

#ifdef PROFILE
  struct timeval start,end;
  gettimeofday(&start,NULL);
#endif

  io_good = false;
  int error = 0;

  if (isRoot()) { // commander, analyze file header
    
    ifstream input(rd_arg.FileName);
    if ( !input.good() )
      {
	if(!UniqueID()){ printf("%s:%s Could not open file ptr %p for input.\n",cname,fname,rd_arg.FileName); fflush(stdout); }
	if(!UniqueID()){ printf("%s:%s Filename %s\n",cname,fname,rd_arg.FileName); fflush(stdout); }
	error = 1;
      }

    if(!error) {
      hd.read(input);
      input.close();
    }
  }
  if(synchronize(error) != 0)
    ERR.FileR(cname, fname, rd_arg.FileName);
  log();

  Broadcast(&hd.data_start, sizeof(streamoff));
  broadcastInt(&hd.recon_row_3);
  //  cout << "recon_row_3 = " << hd.recon_row_3 << endl;



  // check all conditions between FILE and GJP
  int nx = rd_arg.Xnodes() * rd_arg.XnodeSites();
  int ny = rd_arg.Ynodes() * rd_arg.YnodeSites();
  int nz = rd_arg.Znodes() * rd_arg.ZnodeSites();
  int nt = rd_arg.Tnodes() * rd_arg.TnodeSites();

  if(isRoot()) {
    if(hd.dimension[0] != nx || hd.dimension[1] != ny || hd.dimension[2] != nz || hd.dimension[3] != nt) {
      VRB.Flow(cname,fname,"Dimensions in file DISAGREE with GlobalJobParameter!\n");
      VRB.Flow(cname,fname,"In File: %d x %d x %d x %d\n",
	       hd.dimension[0],hd.dimension[1], hd.dimension[2], hd.dimension[3]);
      VRB.Flow(cname,fname,"In GJP:  %d x %d x %d x %d\n",nx, ny, nz, nt);
      error = 1;
    }

// Turned off the boundary check, as it is inconsistent  with NERSC convention.
// 04/03/05, CJ
  }

  if(synchronize(error) != 0)  
    ERR.General(cname, fname, "Wrong Parameters Specified\n");

  // see if file Floating Points is acceptable
  if(isRoot()) {
    fpconv.setFileFormat(hd.floating_point);
  }
  VRB.Flow(cname,fname,"FileFormat=%d",hd.floating_point);
  broadcastInt((int*)&fpconv.fileFormat);
  if(fpconv.fileFormat == FP_UNKNOWN) {
    ERR.General(cname,fname, "Data file Floating Point Format UNKNOWN\n");
  }
  
  VRB.Flow(cname,fname,"A copy of header info from file:\n");
  if(isRoot())  hd.show();

  int data_per_site = hd.recon_row_3 ? 4*12 : 4*18;

  // read lattice data, using parallel style or serial (all on node 0) style
  unsigned int csum;

  //CK: removed hardcoded IO style. Replaced with default IO style in constructor
// #if TARGET != QCDOC
//   setSerial();
// #endif

  log();

  VRB.Result(cname,fname,"Reading configuation to address: %p parIO=%d \n", rd_arg.StartConfLoadAddr,parIO());
  if(parIO()) {
    ParallelIO pario(rd_arg);
    if(!doGparityReconstructUstarField() ){
      if(!UniqueID()) printf("ReadLatticeParallel is disabling reconstruction bit on the ParallelIO object\n");
      pario.disableGparityReconstructUstarField();
    }

    if(! pario.load((char*)rd_arg.StartConfLoadAddr, data_per_site, sizeof(Matrix)*4,
		    hd, fpconv, 4, &csum))  
      ERR.General(cname,fname,"Load Failed\n");  // failed to load
  }
  else {
    SerialIO serio(rd_arg);
    if(!doGparityReconstructUstarField() ) serio.disableGparityReconstructUstarField();

    if(! serio.load((char*)rd_arg.StartConfLoadAddr, data_per_site, sizeof(Matrix)*4,
		    hd, fpconv, 4, &csum))
      ERR.General(cname,fname,"Load Failed\n");  // failed to load
  }

  log();

//  printf("Node %d: lattice read csum=%x\n",UniqueID(),csum);
  //  cout << "loader finish, csum = " << hex << csum << dec << endl << endl;
  //  cout << "loader done" << endl << endl;


  // After reading...
  // STEP 1: checksum
  if(rd_arg.Scoor() == 0)
    csum = globalSumUint(csum);
  else
    globalSumUint(0);

  if(isRoot()) {
    if( hd.checksum != csum ) {
      VRB.Flow(cname,fname, "CheckSUM error !! Header: %x  Host calc: %x\n",hd.checksum,csum);
      
      printf("Node %d: CheckSUM error !! Header: %x  Host calc: %x\n",UniqueID(),hd.checksum,csum);
      error = 1;
    }
    else
      VRB.Flow(cname,fname,"CheckSUM is ok\n");
  }

  if(synchronize(error) != 0) 
    ERR.General(cname, fname, "Checksum error\n");


  // STEP 2: reconstruct row 3
  int size_matrices( rd_arg.XnodeSites() * rd_arg.YnodeSites() 
		     * rd_arg.ZnodeSites() * rd_arg.TnodeSites() * 4); 

  Matrix * lpoint = rd_arg.StartConfLoadAddr;

  if(hd.recon_row_3) {
    VRB.Flow(cname,fname,"Reconstructing row 3\n");
    int nstacked = 1;
    
    //CK: rather than saving/loading the U* links, we save just the U links and reconstruct the U* links at load-time
    if(GJP.Gparity() && !doGparityReconstructUstarField() ) nstacked = 2;     //if this option is turned off, load and check both the U and U* links

    for(int stk=0;stk<nstacked;stk++){
    for(int mat=0; mat<size_matrices; mat++) {
	Float * rec = (Float*)&lpoint[mat + stk*size_matrices];
      // reconstruct the 3rd row
      rec[12] =  rec[2] * rec[10] - rec[3] * rec[11] - rec[4] * rec[8] + rec[5] * rec[9];
      rec[13] = -rec[2] * rec[11] - rec[3] * rec[10] + rec[4] * rec[9] + rec[5] * rec[8];
      rec[14] = -rec[0] * rec[10] + rec[1] * rec[11] + rec[4] * rec[6] - rec[5] * rec[7];
      rec[15] =  rec[0] * rec[11] + rec[1] * rec[10] - rec[4] * rec[7] - rec[5] * rec[6];
      rec[16] =  rec[0] * rec[ 8] - rec[1] * rec[ 9] - rec[2] * rec[6] + rec[3] * rec[7];
      rec[17] = -rec[0] * rec[ 9] - rec[1] * rec[ 8] + rec[2] * rec[7] + rec[3] * rec[6];
    }
    }
  }

  if(GJP.Gparity() && doGparityReconstructUstarField() ){
    //regenerate U* links
    for(int mat=0;mat<size_matrices; mat++) {
      lpoint[mat + size_matrices].Conj(lpoint[mat]);
    }
  }

  // STEP 3: check plaq and linktrace
  if(lat.GaugeField() != lpoint) lat.GaugeField(lpoint);
  Float plaq,linktrace;
  if(! CheckPlaqLinktrace(lat,rd_arg, hd.plaquette, hd.link_trace, plaq,linktrace))  
    ERR.General(cname,fname,"Plaquette (%e %e) or Link trace (%e %e) check failed\n",hd.plaquette,plaq,  hd.link_trace,linktrace);


#ifdef PROFILE
  gettimeofday(&end,NULL);
  print_flops(cname,"read",0,&start,&end);
#endif

  io_good = true;

  log();
  finishLogging();

  lat.ClearSmeared();
  VRB.FuncEnd(cname,fname);
}

void ReadLatticeParallel::readXYTZ(Lattice & lat, const QioArg & rd_arg)
{
  const char * fname = "read()";
  VRB.Func(cname,fname);
  
  char loginfo[strlen(rd_arg.FileName) + 10];
  sprintf(loginfo,"Load %s",rd_arg.FileName);
  startLogging(loginfo);

#ifdef PROFILE
  struct timeval start,end;
  gettimeofday(&start,NULL);
#endif

  io_good = false;
  int error = 0;

  if (isRoot()) { // commander, analyze file header
    
    ifstream input(rd_arg.FileName);
    if ( !input.good() )
      {
	if(!UniqueID()){ printf("%s:%s Could not open file ptr %p for input.\n",cname,fname,rd_arg.FileName); fflush(stdout); }
	if(!UniqueID()){ printf("%s:%s Filename %s\n",cname,fname,rd_arg.FileName); fflush(stdout); }
	error = 1;
      }

    if(!error) {
      hd.read(input);
      input.close();
    }
  }
  if(synchronize(error) != 0)
    ERR.FileR(cname, fname, rd_arg.FileName);
  log();

  Broadcast(&hd.data_start, sizeof(streamoff));
  broadcastInt(&hd.recon_row_3);



  // check all conditions between FILE and GJP
  int nx = rd_arg.Xnodes() * rd_arg.XnodeSites();
  int ny = rd_arg.Ynodes() * rd_arg.YnodeSites();
  int nz = rd_arg.Znodes() * rd_arg.ZnodeSites();
  int nt = rd_arg.Tnodes() * rd_arg.TnodeSites();

  if(isRoot()) {
    if(hd.dimension[0] != nx || hd.dimension[1] != ny || hd.dimension[3] != nz || hd.dimension[2] != nt) {
      VRB.Result(cname,fname,"Rotated Dimensions in file DISAGREE with GlobalJobParameter!\n");
      VRB.Result(cname,fname,"In File: %d x %d x %d x %d\n",
	       hd.dimension[0],hd.dimension[1], hd.dimension[2], hd.dimension[3]);
      VRB.Result(cname,fname,"In GJP:  %d x %d x %d x %d\n",nx, ny, nz, nt);
      error = 1;
    }

  }

  if(synchronize(error) != 0)  
    ERR.General(cname, fname, "Wrong Parameters Specified\n");

  // see if file Floating Points is acceptable
  if(isRoot()) {
    fpconv.setFileFormat(hd.floating_point);
  }
  VRB.Flow(cname,fname,"FileFormat=%d",hd.floating_point);
  broadcastInt((int*)&fpconv.fileFormat);
  if(fpconv.fileFormat == FP_UNKNOWN) {
    ERR.General(cname,fname, "Data file Floating Point Format UNKNOWN\n");
  }
  
  VRB.Flow(cname,fname,"A copy of header info from file:\n");
  if(isRoot())  hd.show();

  int data_per_site = hd.recon_row_3 ? 4*12 : 4*18;

  // read lattice data, using parallel style or serial (all on node 0) style
  unsigned int csum;

  //CK: removed hardcoded IO style. Replaced with default IO style in constructor
// #if TARGET != QCDOC
//   setSerial();
// #endif

  log();

  VRB.Flow(cname,fname,"Reading configuation to address: %p\n", rd_arg.StartConfLoadAddr);
  if(parIO()) {
    ParallelIO pario(rd_arg);
    if(!doGparityReconstructUstarField() ){
      if(!UniqueID()) printf("ReadLatticeParallel is disabling reconstruction bit on the ParallelIO object\n");
      pario.disableGparityReconstructUstarField();
    }
// Should be fixed
    if(! pario.loadXYTZ((char*)rd_arg.StartConfLoadAddr, data_per_site, sizeof(Matrix)*4,
//    if(! pario.load((char*)rd_arg.StartConfLoadAddr, data_per_site, sizeof(Matrix)*4,
		    hd, fpconv, 4, &csum))  
      ERR.General(cname,fname,"Load Failed\n");  // failed to load
  }
  else {
   ERR.General(cname,fname,"Not implemeted for Serial IO");
  }

  log();

//  printf("Node %d: lattice read csum=%x\n",UniqueID(),csum);


  // After reading...
  // STEP 1: checksum
  if(rd_arg.Scoor() == 0)
    csum = globalSumUint(csum);
  else
    globalSumUint(0);

  if(isRoot()) {
    if( hd.checksum != csum ) {
      VRB.Flow(cname,fname, "CheckSUM error !! Header: %x  Host calc: %x\n",hd.checksum,csum);
      
      printf("Node %d: CheckSUM error !! Header: %x  Host calc: %x\n",UniqueID(),hd.checksum,csum);
      error = 1;
    }
    else
      VRB.Flow(cname,fname,"CheckSUM is ok\n");
  }

  if(synchronize(error) != 0) 
    ERR.General(cname, fname, "Checksum error\n");


  // STEP 2: reconstruct row 3
  int size_matrices( rd_arg.XnodeSites() * rd_arg.YnodeSites() 
		     * rd_arg.ZnodeSites() * rd_arg.TnodeSites() * 4); 

  Matrix * lpoint = rd_arg.StartConfLoadAddr;

  if(hd.recon_row_3) {
    VRB.Flow(cname,fname,"Reconstructing row 3\n");
    int nstacked = 1;
    
    //CK: rather than saving/loading the U* links, we save just the U links and reconstruct the U* links at load-time
    if(GJP.Gparity() && !doGparityReconstructUstarField() ) nstacked = 2;     //if this option is turned off, load and check both the U and U* links

    for(int stk=0;stk<nstacked;stk++){
    for(int mat=0; mat<size_matrices; mat++) {
	Float * rec = (Float*)&lpoint[mat + stk*size_matrices];
      // reconstruct the 3rd row
      rec[12] =  rec[2] * rec[10] - rec[3] * rec[11] - rec[4] * rec[8] + rec[5] * rec[9];
      rec[13] = -rec[2] * rec[11] - rec[3] * rec[10] + rec[4] * rec[9] + rec[5] * rec[8];
      rec[14] = -rec[0] * rec[10] + rec[1] * rec[11] + rec[4] * rec[6] - rec[5] * rec[7];
      rec[15] =  rec[0] * rec[11] + rec[1] * rec[10] - rec[4] * rec[7] - rec[5] * rec[6];
      rec[16] =  rec[0] * rec[ 8] - rec[1] * rec[ 9] - rec[2] * rec[6] + rec[3] * rec[7];
      rec[17] = -rec[0] * rec[ 9] - rec[1] * rec[ 8] + rec[2] * rec[7] + rec[3] * rec[6];
    }
    }
  }

  if(GJP.Gparity() && doGparityReconstructUstarField() ){
    //regenerate U* links
    for(int mat=0;mat<size_matrices; mat++) {
      lpoint[mat + size_matrices].Conj(lpoint[mat]);
    }
  }

  // STEP 3: check plaq and linktrace
  if(lat.GaugeField() != lpoint) lat.GaugeField(lpoint);
  Float plaq,linktrace;
  if(! CheckPlaqLinktrace(lat,rd_arg, hd.plaquette, hd.link_trace, plaq,linktrace))  
    ERR.General(cname,fname,"Plaquette (%e %e) or Link trace (%e %e) check failed\n",hd.plaquette,plaq,  hd.link_trace,linktrace);
//  lat.Print();

#ifdef PROFILE
  gettimeofday(&end,NULL);
  print_flops(cname,"read",0,&start,&end);
#endif

  io_good = true;

  log();
  finishLogging();

  lat.ClearSmeared();
  VRB.FuncEnd(cname,fname);
};
bool ReadLatticeParallel::CheckPlaqLinktrace(Lattice &lat, const QioArg & rd_arg,
					     const Float plaq_inheader, const Float linktrace_inheader,
					     Float &plaq, Float &linktrace) 
{
  const char * fname = "CheckPlaqLinktrace()";
  int error = 0;

  plaq = lat.SumReTrPlaq() / 18.0 / rd_arg.VolSites() ;

  //some old lattices that did not use the reconstruct option also did not divide the plaquette by 2
  if(GJP.Gparity() && doGparityReconstructUstarField() ) plaq/=2;

  Float devplaq(0.0);
  if(isRoot()) {
    devplaq = fabs(  (plaq - plaq_inheader) / plaq ) ;
    printf("%s::%s: plaquette::\n  calc: %0.8e  header: %0.8e   rel.dev.: %0.8e\n",
	     cname,fname,plaq, plaq_inheader, devplaq);
  }

  //CK: G-parity; ReTr should be the same for U and U*
  //int nstacked = 1;
  //if(GJP.Gparity()) nstacked = 2;

  //for(int stk=0;stk<nstacked;stk++){ //after 1 iteration, m will point to second stacked field
  linktrace=0.;
  int is;
  Matrix *m =  lat.GaugeField(); 
  for(is=0;is< rd_arg.VolNodeSites()*4; is++){
    linktrace += m->ReTr();
    m++;
  }
  if(rd_arg.Scoor() == 0) 
    linktrace = globalSumFloat(linktrace) / (rd_arg.VolSites()*12.0);
  else
    globalSumFloat(0.0);

  if(isRoot()) {
    Float devlinktrace =   
      fabs(  (linktrace - linktrace_inheader) / linktrace );

    printf("%s::%s: linktrace::\n  calc: %0.8e  header: %0.8e   rel.dev.: %0.8e\n",
	     cname,fname,linktrace, linktrace_inheader, devlinktrace);
  
    Float chkprec = rd_arg.CheckPrecision;

//  CJ:  turning off the link trace test, as some lattices have a very small link trace by accident.
//    if(devplaq > chkprec || devlinktrace > chkprec) {
//      VRB.Flow(cname,fname, "Plaquette and/or Link trace different from header\n");
    if(devplaq > chkprec) {
      VRB.Flow(cname,fname, "Plaquette different from header\n");
      error = 1;
      }
    }
  //}

  if(synchronize(error) != 0) return false;

  return true;
}

// added EES for getting header info consistent on every node
int ReadLatticeParallel::getSequenceNumber(){

  int tmp=hd.sequence_number;
  broadcastInt(&tmp);
  return tmp;

}




CPS_END_NAMESPACE
