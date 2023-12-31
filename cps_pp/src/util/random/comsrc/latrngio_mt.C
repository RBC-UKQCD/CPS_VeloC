#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <util/qcdio.h>
#include <util/latrngio.h>
#include <util/qioarg.h>
#include <util/iostyle.h>
#include <util/intconv.h>
#include <util/random.h>
#include <iostream>
#include <assert.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <util/time_cps.h>
#include <unistd.h>

#ifdef USE_C11_RNG
CPS_START_NAMESPACE
using namespace std;


///////////////////////////////////////////////////////////////
// Lat Random Generator Read  functions 

LatRngRead::LatRngRead() 
  : LatRngIO(), cname("LatRngRead")
{ }

LatRngRead::~LatRngRead() 
{ }

void LatRngRead::read(RNGSTATE *mtran_dump,
		      const QioArg & rd_arg) { 
  const char * fname = "read()";
  VRB.Func(cname,fname);

  intconv.testHostFormat(sizeof(RNGSTATE));

  char loginfo[100];
  sprintf(loginfo,"Load %s",rd_arg.FileName);
  startLogging(loginfo);

  io_good = false;
  int error = 0;

  if(isRoot()) {
    // all open file and check error
    ifstream input(rd_arg.FileName);
    if ( !input.good() ) {
      //      VRB.Flow(cname,fname,"Could not open file: [%s] for input.\n",rd_arg.FileName);
      error = 1;
    }

    if(!error) {
      hd.read(input);
      input.close();
    }
  }
  // first sync point
  // executed by all, sync and share error status information
  if(synchronize(error) != 0)   
    ERR.FileR(cname, fname, rd_arg.FileName);
  log();
  int temp_start=(int)hd.data_start;
  broadcastInt(&temp_start);
  hd.data_start=temp_start;


  if(isRoot()) {
#ifdef USE_C11_MT
    if(hd.datatype != "LATTICE_RNG_C11_MT19937"){ 
#elif (defined USE_C11_RANLUX)
    if(hd.datatype != "LATTICE_RNG_C11_RANLUX48"){
#else
    if(hd.datatype != "LATTICE_RNG_C11_SITMO"){
#endif
      VRB.Flow(cname,fname,"Invalid RNG type: %s\n",hd.datatype.c_str());
      error = 1;
    }
  }
  if(synchronize(error) != 0) 
    ERR.General(cname,fname,"Invalid RNG type\n");

  // check dimensions, etc
  int nx = rd_arg.Xnodes() * rd_arg.XnodeSites();
  int ny = rd_arg.Ynodes() * rd_arg.YnodeSites();
  int nz = rd_arg.Znodes() * rd_arg.ZnodeSites();
  int nt = rd_arg.Tnodes() * rd_arg.TnodeSites();
  int ns = rd_arg.Snodes() * rd_arg.SnodeSites();

  if(isRoot()) {
    //    cout << "Lattice Dimensions in File: " << hd.dimension[0] <<"x" << hd.dimension[1] <<"x" << hd.dimension[2] <<"x" << hd.dimension[3] <<"x" << hd.dimension[4]<< endl;
    if(hd.dimension[0] != nx || hd.dimension[1] != ny || hd.dimension[2] != nz || hd.dimension[3] != nt 
       || hd.dimension[4]!=ns) {
      VRB.Flow(cname,fname,"Dimensions in file DISAGREE with GlobalJobParameter!\n");
      VRB.Flow(cname,fname,"In File: %d x %d x %d x %d x %d\n",
	       hd.dimension[0],hd.dimension[1], hd.dimension[2], hd.dimension[3], hd.dimension[4]);
      VRB.Flow(cname,fname,"In GJP:  %d x %d x %d x %d x %d\n",nx, ny, nz, nt, ns);
      error = 1;
    }
  }

  if(synchronize(error) != 0) 
    ERR.General(cname, fname, "Wrong parameters specified\n");

  // check floating point format
  if(isRoot()) {
    intconv.setFileFormat(hd.int_format);
  }
  broadcastInt((int*)&intconv.fileFormat);
  if(intconv.fileFormat == INT_UNKNOWN) {
    ERR.General(cname,fname,"Data file Integer format UNKNOWN\n");
  }

  if(isRoot())  hd.show();

  // the UGrandomGenerator has multiple type members, need a universal
  // format to make it cross-platform
  int size_rng_ints = LRG.StateSize();
  int size_rng_chars = size_rng_ints * intconv.fileIntSize();

  QioArg rng_arg(rd_arg);
  rng_arg.cutHalf();

  unsigned int csum[2] = {0}, pos_dep_csum[2] = {0};
  Float RandSum[2]={0.0}, Rand2Sum[2]={0.0};

  log();


  VRB.Result(cname,fname,"parIO()=%d\n",parIO());

  if(parIO()) {
    ParallelIO pario(rng_arg);
    VRB.Flow(cname,fname, "Start Loading  RNGs\n");
    if(! pario.load((char*)mtran_dump, size_rng_ints, size_rng_chars,
		    hd, intconv, 4,
		    &csum[0], &pos_dep_csum[0], &RandSum[0], &Rand2Sum[0]))
      ERR.General(cname, fname, "Loading failed\n");

    VRB.Flow(cname,fname,"Node %d - 4D: csum=%x, order_csum=%x\n",
	       UniqueID(),csum[0],pos_dep_csum[0]);
//    printf("Node %d - 4D: csum=%x, order_csum=%x\n",
//	       UniqueID(),csum[1],pos_dep_csum[1]);

  }
#if 1
  else {
    SerialIO serio(rng_arg);
    VRB.Result(cname,fname, "Start Loading 4-D RNGs\n");

    if(! serio.load((char*)mtran_dump, size_rng_ints, size_rng_chars,
		    hd, intconv, 4,
		    &csum[0], &pos_dep_csum[0], &RandSum[0], &Rand2Sum[0]))
      ERR.General(cname, fname, "Loading failed\n");
  }
#endif

  if (isRoot())
  cout << "Starting logging" << endl;

  log();

  if (isRoot())
  cout << "Starting globalSumUint()" << endl;

  // Step 2.1: verify checksum
//  csum[0] += csum[1];
  csum[0] = globalSumUint(csum[0]);

  if (isRoot())
  cout << "Finish globalSumUint()" << endl;

  if(isRoot()) {
    if( hd.checksum != csum[0] ) {
      VRB.Result(cname,fname, "CheckSUM error!! Header:%x  Host calc:%x\n",hd.checksum,csum[0]);
      error = 1;
    }
    else
      VRB.Result(cname,fname,"CheckSUM is ok\n");
  }

  if(synchronize(error) != 0) 
 //     ERR.General(cname,fname, "CheckSUM error\n");
  { 
    VRB.Warn(cname,fname, "CheckSUM error\n");
    log(); finishLogging();
    return;
  }

  // Step 2.2: verify position-dep. Checksum
//  pos_dep_csum[0] += pos_dep_csum[1];
  pos_dep_csum[0] = globalSumUint(pos_dep_csum[0]);
  if(isRoot()) {
    // pos_dep_csum could be absent
    if( hd.pos_dep_csum > 0 && hd.pos_dep_csum != pos_dep_csum[0] ) { 
      VRB.Result(cname,fname, "Position Dependent CheckSUM error!! Header:%x  Host_calc:%x\n",
	       hd.pos_dep_csum, pos_dep_csum[0]);
      error = 1;
    }
    else
      VRB.Result(cname,fname, "Position Dependent CheckSUM is ok\n");
  }

  if(synchronize(error) != 0) 
//    ERR.General(cname, fname, "Position-Dependent Checksum error\n");
  { 
    VRB.Warn(cname, fname, "Position-Dependent Checksum error\n");
    log(); finishLogging();
    return;
  }


  // STEP 3: Verify Rand Average and Variance
  RandSum[0] += RandSum[1];
  Rand2Sum[0] += Rand2Sum[1];
  uint64_t total_rngs_4d = rng_arg.VolSites();
  Float RandAvg = globalSumFloat(RandSum[0]) / (total_rngs_4d);
  Float RandVar = globalSumFloat(Rand2Sum[0]) / (total_rngs_4d)
                  - RandAvg * RandAvg;
  if(isRoot()) {
  VRB.Flow(cname,fname, "Average::  calc: %lf  header: %lf  rel.dev.: %lf\n",
	   RandAvg, hd.average, fabs((RandAvg-hd.average)/RandAvg));
  VRB.Flow(cname,fname, "Variance:: calc: %lf  header: %lf  rel.dev.: %lf\n",
	   RandVar, hd.variance,fabs((RandVar-hd.variance)/RandVar));
  }

  io_good = true;

  log();
  finishLogging();

  VRB.FuncEnd(cname,fname);

}


///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// LatRngWrite functions

LatRngWrite::LatRngWrite() 
  : LatRngIO(), cname("LatRngWrite")
{ }

LatRngWrite::~LatRngWrite()
{ }


#define PROFILE
void LatRngWrite::write(RNGSTATE *mtran_dump,
			const QioArg & wt_arg) {
  const char * fname = "write()";
  VRB.Func(cname,fname);

  char loginfo[100];
  sprintf(loginfo, "Unload %s", wt_arg.FileName);
  startLogging(loginfo);
#ifdef PROFILE 
  Float dtime = -dclock();
#endif

  io_good = false;
  int error = 0;


  VRB.Result(cname,fname,"wt_arg.FileIntFormat= %s\n",intconv.name(wt_arg.FileIntFormat));
  // some cross-platform issued to be discussed
  intconv.testHostFormat(sizeof(RNGSTATE));
  if(sizeof(RNGSTATE) != intconv.size(wt_arg.FileIntFormat))
  ERR.General(cname,fname,"sizeof(RNGSTATE)(%d) != intconv.size(wt_arg.FileIntFormat(%s))(%d)\n",
	sizeof(RNGSTATE) , intconv.name(wt_arg.FileIntFormat), intconv.size(wt_arg.FileIntFormat) );
  intconv.setFileFormat(wt_arg.FileIntFormat);
  //  VRB.Flow(cname,fname, "Set File INT format : %s\n", intconv.name(intconv.fileFormat));

  // size information
  int nx = wt_arg.Xnodes() * wt_arg.XnodeSites();
  int ny = wt_arg.Ynodes() * wt_arg.YnodeSites();
  int nz = wt_arg.Znodes() * wt_arg.ZnodeSites();
  int nt = wt_arg.Tnodes() * wt_arg.TnodeSites();
  int ns = wt_arg.Snodes() * wt_arg.SnodeSites();

  //  cout << "Lattice size = " << nx <<"x"<<ny<<"x"<<nz<<"x"<<nt<<"x"<<ns<<endl;

  // num of integers per RNG (NOTE: not all data in RNG are saved)
  int size_rng_ints = LRG.StateSize();
  //  cout << "size_rng_ints = " << size_rng_ints << endl;
  int size_rng_chars = size_rng_ints * intconv.fileIntSize();

//  should be set before writing
//  setParallel();
//  setSerial();
  VRB.Result(cname,fname,"parIO()=%d ConcurIONumber=%d\n",parIO(),wt_arg.ConcurIONumber);

  log();

  // start writing file
  fstream output;
//  string name_4d(wt_arg.FileName);
//  name_4d += ".4d";

  if(parIO()) {
    FILE *fp = Fopen(wt_arg.FileName,"w");
    Fclose(fp);
	Float temp;glb_sum_five(&temp);
    output.open(wt_arg.FileName);
    if(!output.good())  {
          VRB.Flow(cname,fname,"Could not open file: [%s] for output.\n",wt_arg.FileName);
      error = 1;
    }
  }
#if 1
  else {
    if(isRoot()) {
      FILE *fp = Fopen(wt_arg.FileName,"w");
      Fclose(fp);
//      fp = Fopen(name_4d.c_str(),"w");
//      Fclose(fp);
	  
      output.open(wt_arg.FileName);
//    output_4d.open(name_4d.c_str());
      if(!output.good()) {
	//	VRB.Flow(cname,fname,"Could not open file: [%s] for output.\n",wt_arg.FileName);
      printf("Node %d:Could not open file: [%s] for output.\n",UniqueID(),wt_arg.FileName);
	error = 1;
      }
    }    
  }  
#endif
//  if(synchronize(error) > 0)  
//    ERR.General(cname,fname,"Could not open file: [%s] for output.\n",wt_arg.FileName);
  if (error)
    ERR.General(cname,fname,"Could not open file: [%s] for output.\n",wt_arg.FileName);
//    printf("Node %d: says opening %s failed\n",UniqueID(),wt_arg.FileName);

  // write header
  if(isRoot()) {
    hd.init(wt_arg, intconv.fileFormat);
	VRB.Result(cname,fname,"Writing header for %s\n",wt_arg.FileName);
    hd.write(output);
//	VRB.Result(cname,fname,"Writing header for %s\n",name_4d.c_str());
//    hd.write(output_4d);
	VRB.Result(cname,fname,"Done\n");
  }
  int temp_start=(int)hd.data_start;
  broadcastInt(&temp_start); // from 0 to all
  hd.data_start=temp_start;
  log();

  unsigned int csum[2]={0}, pos_dep_csum[2] = {0};
  Float RandSum[2]={0.0}, Rand2Sum[2] = {0.0};

  QioArg rng_arg(wt_arg);
  rng_arg.cutHalf();

    VRB.Flow(cname,fname,"Start Unloading MT19937 RNGs\n");
  if(parIO()) {
    ParallelIO pario(rng_arg);
    
    if(! pario.store(output, (char*)mtran_dump, size_rng_ints, 
		     size_rng_chars, hd, intconv, 4,
		     &csum[0], &pos_dep_csum[0], &RandSum[0], &Rand2Sum[0]))
        ERR.General(cname, fname, "Unloading Failed\n");
    
    VRB.Flow(cname,fname,"Node %d - 4D: csum=%x, order_csum=%x\n",
             UniqueID(),csum[1],pos_dep_csum[1]);
    //    printf("Node %d - 4D: csum=%x, order_csum=%x\n",
    //	       UniqueID(),csum[1],pos_dep_csum[1]);
  } else {

    SerialIO serio(rng_arg);
    if(! serio.store(output, (char*)mtran_dump, size_rng_ints, 
		     size_rng_chars, hd, intconv, 4,
		     &csum[0], &pos_dep_csum[0], &RandSum[0], &Rand2Sum[0]))
      ERR.General(cname, fname, "Unloading Failed\n");
  }

  log();

  // fill in verification information
//  csum[0] += csum[1];
  csum[0] = globalSumUint(csum[0]);
//  pos_dep_csum[0] += pos_dep_csum[1];
  pos_dep_csum[0] = globalSumUint(pos_dep_csum[0]);

//  RandSum[0] += RandSum[1];
//  Rand2Sum[0] += Rand2Sum[1];
  int total_rngs_4d = rng_arg.VolSites();
//  int total_rngs_5d = total_rngs_4d * rng_arg.Snodes() * rng_arg.SnodeSites(); 
  Float RandAvg = globalSumFloat(RandSum[0]) / (total_rngs_4d);
  Float RandVar = globalSumFloat(Rand2Sum[0]) / (total_rngs_4d)
                  - RandAvg * RandAvg;

  if(isRoot()) {
    hd.fillInCheckInfo(output, csum[0], pos_dep_csum[0], RandAvg, RandVar);
    if ( !output.good() )  error = 1; 
  }

  if(synchronize(error) != 0)
    ERR.General(cname, fname, "Writing checksum and other info failed\n");
  log();

  if(parIO()) {
    output.close();
//	output_4d.close();
  } else
    if(isRoot()){
	  output.close();
//	  output_4d.close();
	}
  
  io_good = true;

  finishLogging();
#ifdef PROFILE 
  dtime +=dclock();
  print_time(cname,fname,dtime);
#endif

  VRB.FuncEnd(cname,fname);
}



CPS_END_NAMESPACE
#endif 
