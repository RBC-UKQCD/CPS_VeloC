#ifndef __LAT_HEADER__
#define __LAT_HEADER__

#include <iostream>
#include <string>
#include <map>

#include <util/lattice.h>
#include <util/gjp.h>
#include <util/qioarg.h>
#include <util/fpconv.h>
#include <util/intconv.h>


CPS_START_NAMESPACE

// GCFheaderPar class
// header parser for parallel IO
// removed "exit()"'s, others same as class GCFheader
typedef std::map<std::string, std::string> GCFHMapParT;

class GCFheaderPar
{
private:

  GCFHMapParT headerMap;
  bool  prevFound;

public:
  
  inline bool found() { return prevFound; }

    bool    add( std::string key_eq_value );

    int     asInt   ( std::string key );
    unsigned int     asHex   ( std::string key );
    Float   asFloat ( std::string key );
    std::string  asString( std::string key );

  void Show() const;
};


// base class for all types of headers
class LatHeaderBase {
 public:
  enum HEADER_TYPES {
    LATTICE_HEADER,
    LATRNG_HEADER,
    LATTICE_NERSC
  };
  virtual enum HEADER_TYPES headerType() const = 0;
  virtual std::streamoff dataStart() const = 0;
  virtual void read(std::istream & fin) = 0;
  virtual void write(std::ostream & fout) = 0;
  virtual void fillInCheckInfo(std::ostream & fout, unsigned int cs, unsigned int pdcs,
			       const Float calc1, const Float calc2) const = 0;
  virtual void show() const = 0;

  LatHeaderBase() { }
  virtual ~LatHeaderBase() { }

  std::streamoff data_start;

 protected:
  GCFheaderPar hd;
};

// header specification/interpretation
class LatNERSCHeader : public LatHeaderBase {
 public:
     const char *cname;
  // header strings
    std::string hdr_version;
//  int recon_row_3; // determines DATATYPE = 4D_SU3_GAUGE or 4D_SU3_GAUGE_3X3
    std::string storage_format;
    int data_per_site;

  int dimension[4];
//  Float link_trace;
//  Float plaquette;

  BndCndType boundary[4]; 
  unsigned int checksum;

//    std::string ensemble_id ;
//    std::string ensemble_label ;
//  int sequence_number ;
    std::string creator ;
    std::string creator_hardware ;
    std::string creation_date ;
    std::string archive_date ;

  FP_FORMAT floating_point;

  LatNERSCHeader():LatHeaderBase(),cname("LatNERSCHeader") {
//    ensemble_id = "unspecified";
//    ensemble_label = "unspecified";
//    sequence_number = 0;
  }

  void init(const QioArg & qio_arg, FP_FORMAT FileFormat);
  void setHeader(int ndata, const char *CreatorName = NULL, const char *CreatorHardware = NULL );
    virtual void writeContent(std::ostream & fout);
    void write(std::ostream & fout,std::vector<std::string> key, std::vector<std::string> value){
    fout.seekp(0,std::ios::beg);
    fout << "BEGIN_HEADER" << std::endl;
     writeContent(fout);
	if(key.size()!=value.size())
        ERR.General(cname,"write()","key(%d) and value(%d) length mismatch\n",key.size(),value.size());
    for(size_t i =0;i<key.size();i++){
    std::cout << key[i] <<" = " << value[i] << std::endl;
    fout << key[i] <<" = " << value[i] << std::endl;
    }


    fout << "END_HEADER" << std::endl;
    data_start = fout.tellp();
    VRB.Result(cname,"write()","data_start=%d csum_pos=%d\n",data_start,csum_pos);
    }

    void write(std::ostream &fout){
	std::vector <std::string> key;
	std::vector <std::string> value;
	this->write(fout,key,value);
    }

    void fillInChecksum(std::ostream & fout, unsigned int checksum) const;
  
    virtual void read(std::istream & fin);

  void show() const { hd.Show(); }

  enum HEADER_TYPES headerType() const { return LATTICE_NERSC; }
  std::streamoff dataStart() const { return data_start; }
  void fillInCheckInfo(std::ostream & fout, unsigned int cs, unsigned int pdcs,
		       const Float calc1, const Float calc2) const {
    fillInChecksum(fout, cs);
  }

 private:
  int csum_pos;
};

template<FixGaugeType gfix_type> 
class LatGfixHeader : public LatNERSCHeader {
 public:
    const char *cname;
    Float stp_cnd; //stopping condition
//    FixGaugeType gfix_type;
    
    LatGfixHeader():LatNERSCHeader(),cname("LatGfixHeader") {}
    void writeContent(std::ostream & fout);
    void read(std::istream & fin);
//    void setHeader(FixGaugeType fix, Float stp, const char *CreatorName = NULL, const char *CreatorHardware = NULL ){
//       gfix_type=fix;	
//       stp_cnd = stp;
//       LatNERSCHeader::setHeader(18,CreatorName,CreatorHardware);
//    }

};


// header specification/interpretation
class LatticeHeader : public LatHeaderBase {
 public:
  // header strings
    std::string hdr_version;
  int recon_row_3; // determines DATATYPE = 4D_SU3_GAUGE or 4D_SU3_GAUGE_3X3
    std::string storage_format;

  int dimension[4];
  Float link_trace;
  Float plaquette;

  BndCndType boundary[4]; 
  unsigned int checksum;

    std::string ensemble_id ;
    std::string ensemble_label ;
  int sequence_number ;
    std::string creator ;
    std::string creator_hardware ;
    std::string creation_date ;
    std::string archive_date ;

  FP_FORMAT floating_point;

  LatticeHeader() {
    ensemble_id = "unspecified";
    ensemble_label = "unspecified";
    sequence_number = 0;
  }

  void init(const QioArg & qio_arg, FP_FORMAT FileFormat, Float LinkTrace, Float Plaq);
  void setHeader(const char * EnsembleId, const char * EnsembleLabel, const int SequenceNumber,const char *CreatorName = NULL, const char *CreatorHardware = NULL );
    void write(std::ostream & fout);
    void fillInChecksum(std::ostream & fout, unsigned int checksum) const;
  
    void read(std::istream & fin);

  void show() const { hd.Show(); }

  enum HEADER_TYPES headerType() const { return LATTICE_HEADER; }
  std::streamoff dataStart() const { return data_start; }
  void fillInCheckInfo(std::ostream & fout, unsigned int cs, unsigned int pdcs,
		       const Float calc1, const Float calc2) const {
    fillInChecksum(fout, cs);
  }

 private:
  int csum_pos;
};

  /* Lattice Header e.g:
BEGIN_HEADER
HDR_VERSION = 1.0
DATATYPE = 4D_SU3_GAUGE
STORAGE_FORMAT = 1.0
DIMENSION_1 = 4
DIMENSION_2 = 4
DIMENSION_3 = 4
DIMENSION_4 = 4
LINK_TRACE = -0.005858163079
PLAQUETTE  = 0.02623645351
BOUNDARY_1 = PERIODIC
BOUNDARY_2 = PERIODIC
BOUNDARY_3 = PERIODIC
BOUNDARY_4 = ANTIPERIODIC
CHECKSUM = a7db3dee
ENSEMBLE_ID = unspecified
ENSEMBLE_LABEL = unspecified
SEQUENCE_NUMBER = 0
CREATOR = unspecified
CREATOR_HARDWARE = unspecified
CREATION_DATE = unspecified
ARCHIVE_DATE = Thu Jan  1 00:18:49 1970
FLOATING_POINT = IEEE64BIG
END_HEADER
  */



// header specification/interpretation
class LatRngHeader : public LatHeaderBase {
 public:
  // header strings
  std::string hdr_version;
  std::string datatype;
  std::string storage_format;
  static const char *RNGString;


  int dimension[5];
  unsigned int checksum;
  unsigned int pos_dep_csum;

  Float average;
  Float variance;

  std::string creator ;
  std::string creator_hardware ;
  std::string creation_date ;
  std::string archive_date ;

  INT_FORMAT int_format;


  LatRngHeader():datatype(RNGString) {  }

  void init(const QioArg & qio_arg, INT_FORMAT FileFormat);
    void write(std::ostream & fout);
    void fillInCheckInfo(std::ostream & fout, unsigned int cs, unsigned int pdcs, Float avg, Float var) const;
  
    void read(std::istream & fin);

  std::streamoff dataStart() const { return data_start; }
  
  void show() const { hd.Show(); }

  enum HEADER_TYPES headerType() const { return LATRNG_HEADER; }
  

 private:
  int csum_pos;
  int pdcs_pos;
  int avg_pos;
  int var_pos;

};


/*  LatRng Header e.g:
BEGIN_HEADER
HDR_VERSION = 1.0
DATATYPE = LATTICE_RNG_5D_4D
STORAGE_FORMAT = 1.0
DIMENSION_1 = 8
DIMENSION_2 = 8
DIMENSION_3 = 8
DIMENSION_4 = 8
DIMENSION_5 = 1
CHECKSUM = 5298cc92
POS_DEP_CSUM = 38908578
INT_FORMAT = INT32BIG
AVERAGE = -0.1410818313       
VARIANCE = 0.9787710161        
END_HEADER
*/





CPS_END_NAMESPACE

#endif
