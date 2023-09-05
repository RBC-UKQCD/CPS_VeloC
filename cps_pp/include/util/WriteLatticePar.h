#ifndef _WRITELATTICEPAR_H_
#define _WRITELATTICEPAR_H_
// Write the format {load,unload}_lattice 

#include <stdlib.h>	// exit()
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <util/gjp.h>
#include <util/vector.h>
#include <util/lattice.h>
#include <util/verbose.h>
#include <util/error.h>

#include <util/qioarg.h>
#include <util/fpconv.h>
#include <util/iostyle.h>
#include <util/latheader.h>


CPS_START_NAMESPACE

//----------------------------------------
// WriteLatticeParallel class
// A modification to "WriteLattice" class to enable parallel writing

class WriteLatticeParallel : public QioControl
{
  // IoStyle provides a function  IoStyle::store() 
  // which determines Parallel or Serial storing

 private:
  //    FPConv fpconv;
    int csum_pos;
    bool recon_row_3;
    const char *cname;
    int UseParIO;

 public:
    LatticeHeader hd;

 public:
    // ctor for 2-step unloading

    WriteLatticeParallel()  
      : QioControl(), cname("WriteLatticeParallel"){
      setDefault();
    }

    // ctor containing unloading behavior
    WriteLatticeParallel(Lattice & lat, const char * filename,
			 const FP_FORMAT dataFormat = FP_AUTOMATIC, const int recon_row_3 = 1)
      : QioControl(), cname("WriteLatticeParallel"){
      setDefault();
      QioArg  wt_arg(filename, dataFormat, recon_row_3);
      write(lat, wt_arg);
    }

    // ctor containing unloading behavior
    WriteLatticeParallel(Lattice & lat, const QioArg & wt_arg)
      : QioControl(), cname("WriteLatticeParallel"){
      setDefault();
      write(lat, wt_arg);
    }

    ~WriteLatticeParallel() {}

    void setHeader(const char * EnsembleId, const char * EnsembleLabel,
		   const int SequenceNumber) {
      hd.setHeader(EnsembleId, EnsembleLabel, SequenceNumber);
    }

    void write(Lattice & lat, const char * filename,
	       const FP_FORMAT dataFormat = FP_AUTOMATIC, const int recon_row_3 = 1) {
      QioArg  wt_arg(filename, dataFormat, recon_row_3);
      write(lat, wt_arg);
    }
    void write(Lattice & lat, const QioArg & wt_arg);

 public:
    inline void setParallel() { UseParIO = 1; }
    inline void setSerial() { UseParIO = 0; }
    inline int parIO() const { return UseParIO; }
    void setDefault(){
#ifdef USE_SERIAL_IO
      setSerial();
#else
      setParallel();
#endif
    }

};


class WriteLatticeSerial : public WriteLatticeParallel {
 private:
  const char * cname;

 public:
    // ctor for 2-step unloading
    WriteLatticeSerial()
      : WriteLatticeParallel(), cname("WriteLatticeParallel") {
      setDefault();
    }

    // ctor containing unloading behavior
    WriteLatticeSerial(Lattice & lat, const char * filename,
		       const FP_FORMAT dataFormat = FP_AUTOMATIC, const int recon_row_3 = 1)
      : WriteLatticeParallel(), cname("WriteLatticeParallel") {
      setDefault();
      write(lat, filename, dataFormat, recon_row_3);
    }

    // ctor containing unloading behavior
    WriteLatticeSerial(Lattice & lat, const QioArg & wt_arg)
      : WriteLatticeParallel(), cname("WriteLatticeParallel"){
      setDefault();
      write(lat, wt_arg);
    }

    ~WriteLatticeSerial() {}
};


CPS_END_NAMESPACE
#endif
