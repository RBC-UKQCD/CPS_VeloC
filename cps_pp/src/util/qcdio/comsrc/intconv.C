#include <config.h>
#include <string.h>
#include <iostream>
#include <util/data_types.h>
#include <util/intconv.h>
#include <util/qioarg.h>

CPS_START_NAMESPACE
using namespace std;

const char * INT_FORMAT_NAME[] = { 
  "n/a",
  "AUTOMATIC",
  "INT32BIG",
  "INT32LITTLE",
  "INT64BIG",
  "INT64LITTLE"
};

const int INT_FORMAT_ENTRIES = sizeof(INT_FORMAT_NAME)/sizeof(INT_FORMAT_NAME[0]);


IntConv::IntConv() 
  : fileFormat(INT_UNKNOWN), hostFormat(INT_UNKNOWN), cname("IntConv") { 
  testHostFormat();
}

IntConv::IntConv(int size) 
  : fileFormat(INT_UNKNOWN), hostFormat(INT_UNKNOWN), cname("IntConv") { 
  testHostFormat(size);
}

IntConv::~IntConv() {
}


const char * IntConv::name(const enum INT_FORMAT format) {
  return INT_FORMAT_NAME[int(format)];
}

char * IntConv::file2host(char * hbuf, const char * fdat, const int fdat_len) const {
  // trivial case
  if(hostFormat == fileFormat) {
    memcpy(hbuf,fdat,fdat_len*size(hostFormat));
    return hbuf;
  }
  assert(size(hostFormat) == size(fileFormat)) ;
  
  // if needs conversion
  // adjust endian  (currently the only conversion needed...)
  if(big_endian(hostFormat) != big_endian(fileFormat)) {
    if(size(fileFormat)==8){
    byterevn64((type64*)fdat,fdat_len);
    copy64((type64*)hbuf, (type64*)fdat, fdat_len);
    } else {
    byterevn((type32*)fdat,fdat_len);
    copy32((type32*)hbuf, (type32*)fdat, fdat_len);
    }
  }

  return hbuf;
}

char * IntConv::host2file(char *fbuf, const char * hdat, const int hdat_len) const {
  // trivial case
  if(hostFormat == fileFormat) {
    memcpy(fbuf,hdat,hdat_len*size(fileFormat));
    return fbuf;
  }

  assert(size(hostFormat) == size(fileFormat)) ;
  
  // if needs conversion
  // adjust endian (currently the only conversion needed)
  if(big_endian(hostFormat) != big_endian(fileFormat)) {
    if(size(fileFormat)==8){
    copy64((type64*)fbuf, (type64*)hdat, hdat_len);
    byterevn64((type64*)fbuf,hdat_len);
    } else {
    copy32((type32*)fbuf, (type32*)hdat, hdat_len);
    byterevn((type32*)fbuf,hdat_len);
    }
  }
//  if(big_endian(hostFormat) != big_endian(fileFormat)) {  
//    copy32((type32*)fbuf, (type32*)hdat, hdat_len);
//    byterevn((type32*)fbuf, hdat_len);
//  }

  return fbuf;
}


void IntConv::byterevn(type32 w[], int n) const {
  /*  char * buf = (char*)w;
  cout << "First 16 bytes: ";
  for(int i=0;i<16;i++) cout << hex << (unsigned int)buf[i] << " ";
  cout << dec << endl;
  */
  //  cout << "Byte reverse 32 bits" << endl;

  type32 oldv, newv;

  for(int i=0;i<n;i++) {
    oldv = w[i];
    newv = 0;
    for(int j=0;j<4;j++) {
      newv = (newv << 8) | (oldv & 0xff);
      oldv >>= 8;
    }
    w[i] = newv;
  }
  /*
  cout << "First 16 bytes: ";
  for(int i=0;i<16;i++) cout << hex << (unsigned int)buf[i] << " ";
  cout << dec << endl;
  */
}

void IntConv::byterevn64(type64 w[], int n) const{
  /*
  char * buf = (char*)w;
  cout << "First 16 bytes: ";
  for(int i=0;i<16;i++) cout << hex << (unsigned int)buf[i] << " ";
  cout << dec << endl;
  */
  //  cout << "Byte reverse 64 bits" << endl;

  type64 oldv, newv;

  for(int i=0;i<n;i++) {
    oldv = w[i];
    newv = 0;
    for(int j=0;j<8;j++) {
      newv = (newv << 8) | (oldv & 0xff);
      oldv >>= 8;
    }
    w[i] = newv;
  }
  /*
  cout << "First 16 bytes: ";
  for(int i=0;i<16;i++) cout << hex << (unsigned int)buf[i] << " ";
  cout << dec << endl;
  */
}

void IntConv::copy64(type64 tgt[], type64 src[], int n) const{
  double *s = (double*)src;
  double *t = (double*)tgt;
  for(int i=0;i<n;i++)  *t++ = *s++;
}

void IntConv::copy32(type32 tgt[], type32 src[], int n) const {
  float *s = (float*)src;
  float *t = (float*)tgt;
  for(int i=0;i<n;i++)  *t++ = *s++;
}

enum INT_FORMAT  IntConv::testHostFormat(int size ) { // test the type of CPS::Float
  const char * fname = "testHostFormat()";
  // 1. endian
  char end_check[4] = {1,0,0,0};
  uint32_t *lp = (uint32_t *)end_check;
  int host_big;

  if ( *lp == 0x1 ) { 
    //    cout << "Host is little-endian\n";
    host_big = 0;
  } else {
    //    cout << "Host is big-endian\n";
    host_big = 1;
  }

  // 2. pi test
  if (size <0 ) size = sizeof(int);
  if(size == 4) {
    if(host_big)  hostFormat = INT_32BIG;
    else          hostFormat = INT_32LITTLE;
  }
  else 
  if(size == 8) {
    if(host_big)  hostFormat = INT_64BIG;
    else          hostFormat = INT_64LITTLE;
  }
  else {
    ERR.NotImplemented(cname,fname,"IntConv::testHostFormat() : non-32 or 64 bit int not supported\n");
  }
  VRB.Result(cname,fname,"size=%d hostFormat=%s\n",size,name(hostFormat));

  return hostFormat;
}

enum INT_FORMAT  IntConv::setFileFormat(const enum INT_FORMAT dataFormat) {
  const char * fname = "setFileFormat()";

  fileFormat = dataFormat;
  if(dataFormat == INT_AUTOMATIC) 
    fileFormat = hostFormat;
  if(fileFormat == INT_UNKNOWN) {
    VRB.Flow(cname,fname,"Floating point format cannot be INT_UNKNOWN!\n");
  }
  VRB.Result(cname,fname,"fileFormat=%s\n",name(fileFormat));
//  int *crash=NULL; *crash=0;
  return fileFormat;
}

enum INT_FORMAT  IntConv::setFileFormat(const char * desc) {
  const char * fname = "setFileFormat()";
  fileFormat = INT_UNKNOWN;
  for(int i=1;i<INT_FORMAT_ENTRIES;i++) {
    if(!strcmp(INT_FORMAT_NAME[i],desc)) {
      fileFormat = INT_FORMAT(i);
      break;
    }
  }
  if(fileFormat == INT_UNKNOWN) {
    VRB.Flow(cname,fname,"Floating point format \"%s\" not recognized!\n",desc);
  }
  else if(fileFormat == INT_AUTOMATIC)  {
    fileFormat = hostFormat;
  }

  return fileFormat;
}

unsigned int IntConv::checksum(char * data, const int data_len,
			       const enum INT_FORMAT dataFormat) const{
  const char * fname = "checksum()";
  // checksum always done on 32-bits

  enum INT_FORMAT chkFormat = dataFormat;
  if(dataFormat == INT_AUTOMATIC)  chkFormat = fileFormat;

  if(chkFormat == INT_UNKNOWN) {
    VRB.Flow(cname,fname,"checksum data format not recognized!\n");
    return 0;
  }

  int csumcnt = data_len;

  if(big_endian(hostFormat) != big_endian(chkFormat))
    byterevn((type32*)data, csumcnt);

  unsigned int *buf = (unsigned int*)data;
  unsigned int s = 0;
  for(int i=0;i<csumcnt;i++)  {
    s += *buf;
    buf++;
  }

  if(big_endian(hostFormat) != big_endian(chkFormat))
    byterevn((type32*)data, csumcnt);

  return s;
}

unsigned int IntConv::posDepCsum(char * data, const int data_len, 
				 const int dimension, const QioArg & qio_arg, 
				 const int siteid, const int global_id,
				 const enum INT_FORMAT dataFormat) const {
  const char * fname = "posDepCsum()";
  // checksum always done on 32-bits

  enum INT_FORMAT chkFormat = dataFormat;
  if(dataFormat == INT_AUTOMATIC)  chkFormat = fileFormat;

  if(chkFormat == INT_UNKNOWN) {
    VRB.Flow(cname,fname,"checksum data format not recognized!\n");
    return 0;
  }

  int csumcnt = data_len;

  if(big_endian(hostFormat) != big_endian(chkFormat))
    byterevn((type32*)data, csumcnt);

  // s = pdcsum within a site
  unsigned int *buf = (unsigned int*)data;
  unsigned int s = 0;
  for(int i=0;i<csumcnt;i++)  {
    s += *buf * (i+1); // position-dep. checksum
    buf++;
  }

  if(big_endian(hostFormat) != big_endian(chkFormat))
    byterevn((type32*)data, csumcnt);

  // s * (uniqueSiteId+1) = global pdcsum
  int sid = siteid;
  int gid = global_id;

  int stk =0; //CK: for G-parity

  if(sid >= 0) { // calculate global_id via siteid
    int loc[5];
    for(int i=0;i<5;i++) {
      loc[i] = sid % qio_arg.NodeSites(i);
      sid /= qio_arg.NodeSites(i);
      loc[i] += qio_arg.NodeSites(i) * qio_arg.Coor(i);
      if(GJP.Gparity() && i==3){
	stk = sid % 2;
	sid /= 2;
      }
    }
    
    int nstacked = 1;
    if(GJP.Gparity()) nstacked =2;

    //CK: x + nx*y + nx*ny*z + nx*ny*nz*t + nx*ny*nz*nt*stk + nx*ny*nz*nt*nstacked*s
    gid = 0;
    if(dimension==5) gid = loc[4];

    gid = gid * nstacked + stk;

    for(int i=3;i>=0;i--) 
      gid = gid * qio_arg.Nodes(i) * qio_arg.NodeSites(i) + loc[i];
  }

  return s * (gid+1);
}


int IntConv::size(const enum INT_FORMAT datatype) const {
  switch(datatype) {
  case INT_32BIG:
  case INT_32LITTLE:
    return 4;
  case INT_64BIG:
  case INT_64LITTLE:
    return 8;
  default:
    return 0;
  }
}

bool IntConv::big_endian(const enum INT_FORMAT datatype) const {
  switch(datatype) {
  case INT_32BIG:
  case INT_64BIG:
    return true;
  case INT_32LITTLE:
  case INT_64LITTLE:
    return false;
  default:
    return false;
  }
}


CPS_END_NAMESPACE
