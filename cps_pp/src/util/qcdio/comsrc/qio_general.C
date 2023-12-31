#ifdef USE_QIO
#include <config.h>
#include <util/qio_general.h>
#include <util/time_cps.h>


CPS_START_NAMESPACE
using namespace std;


static int qio_sparse_num=1;
void setQioSparseNum(int n) { qio_sparse_num = n;}
int getQioSparseNum() { return qio_sparse_num;}



// the helper functions are outside the class...

int qio_node_number( const int x[])
{
  // node number of node at global-lattice-point x[]=x,y,z,t  

  // return the node, which is writing/reading at x[]


  #ifdef DEBUG_NodeNumber
  printf("UID: %i, called qio_node_number with x: [%i, %i, %i, %i], calculating...\n", UniqueID(),x[0],x[1],x[2],x[3]);
  #endif // DEBUG_NodeNumber


  //convert x in xGrid
  

  int xGrid[4]={0,0,0,0};

  int xGrid_5 = 0;


  for(int ii(0); ii < 4; ++ii)
    {
      int div = x[ii];
      while ( div >= GJP.NodeSites(ii) )
	{
	  ++xGrid[ii];
	  div -= GJP.NodeSites(ii);
	}
    }


  // case 5th dim par.
  if(GJP.Snodes() > 1 )
    {
      int div = x[3] - xGrid[3] * GJP.NodeSites(3);
      int subtr =  GJP.NodeSites(3)/GJP.Snodes();
      while ( div >= subtr )
	{
	  ++xGrid_5;
	  div -= subtr;
	}
    }

  //now get nodeNumber from xGrid
  // assuming that QMP-node numbering is the same as CPS... !!!

  int tmp(0);
  int vol(1);
  

  for( int ii(0); ii <4; ++ii)
    {
      tmp += vol*xGrid[ii];
      vol *= GJP.Nodes(ii);
    }
  

  if(GJP.Snodes() > 1)
      tmp += vol*xGrid_5;



  #ifdef DEBUG_NodeNumber
  printf("UID: %i, called qio_node_number with x: [%i, %i, %i, %i], returned: %i\n", UniqueID(),x[0],x[1],x[2],x[3],tmp);
  #endif // DEBUG_NodeNumber
  
  return tmp;
  
}


int  qio_node_index( const int x[])
{
  
  // returns local index of global coordinates


  #ifdef DEBUG_NodeIndex
  printf("UID: %i, called qio_node_index with x: [%i, %i, %i, %i] calculating...\n",UniqueID(),x[0],x[1],x[2],x[3]);
  #endif // DEBUG_NodeIndex

  
  int xLoc[4];



  for(int ii(0); ii < 4; ++ii)
    xLoc[ii] = x[ii] - GJP.NodeCoor(ii)*GJP.NodeSites(ii);
  
  // now convert local coordinate into local index

 
  
  int tmp(0);
  int vol(1);

  for( int ii(0); ii < 4; ++ii)
    {
      tmp += vol * xLoc[ii];
      vol *= GJP.NodeSites(ii);
    }



  if( tmp >= GJP.VolNodeSites() )
    { 
      printf("ERROR QIO: qio_node_index got index >= local volume");
  printf("UID: %i, called qio_node_index with x: [%i, %i, %i, %i] \n",UniqueID(),x[0],x[1],x[2],x[3]);
  printf("UID: %i, called qio_node_index with NodeCoor: [%i, %i, %i, %i] \n",UniqueID(), GJP.NodeCoor(0), GJP.NodeCoor(1), GJP.NodeCoor(2), GJP.NodeCoor(3));
  printf("UID: %i, called qio_node_index with NodeSites: [%i, %i, %i, %i] \n",UniqueID(), GJP.NodeSites(0), GJP.NodeSites(1), GJP.NodeSites(2), GJP.NodeSites(3));
      exit(-1);
    }


  #ifdef DEBUG_NodeIndex
  printf("UID: %i, called qio_node_index with x: [%i, %i, %i, %i] returned index %i\n",UniqueID(),x[0],x[1],x[2],x[3],tmp);
  #endif // DEBUG_NodeIndex

  return tmp;

}


#ifndef QIO_SPARSE_PARTFILE
// this is the original version, only works for  node == QMP_get_node_number()
void qio_get_coords( int x[], int node, int index)
{

  #ifdef DEBUG_GetCoords
  printf("UID: %i, called qio_get_coords with node: %i, index %i; calc...\n", UniqueID(), node, index);
  #endif // DEBUG_GetCoords

  // returns x[]: global(?) coordinates on this node

  //if( node != UniqueID() )
  if( node != QMP_get_node_number() )
    { 
      printf("ERROR QIO: node number mismatch %d vs %d [%d]\n",
	     node, QMP_get_node_number(),UniqueID());
 
      exit(-1); 
    }

  x[0] = GJP.XnodeCoor()*GJP.NodeSites(0);
  x[1] = GJP.YnodeCoor()*GJP.NodeSites(1);
  x[2] = GJP.ZnodeCoor()*GJP.NodeSites(2);
  x[3] = GJP.TnodeCoor()*GJP.NodeSites(3);


  int localIndex_t;

  if( GJP.Snodes() > 1 )
    {
      localIndex_t = GJP.NodeSites(3)/GJP.Snodes();
      x[3] += GJP.SnodeCoor()*GJP.NodeSites(3)/GJP.Snodes() ; 

    }
  else
    localIndex_t = GJP.NodeSites(3);

#undef USE_SLOW_INDEX
#ifdef USE_SLOW_INDEX
  int xLoc[4]={0,0,0,0};

  for( int ii(0); ii < index; ++ii)
    { 
      if( ++xLoc[0] == GJP.NodeSites(0) ) 
	{ 
	 xLoc[0]=0; 
	 if( ++xLoc[1] == GJP.NodeSites(1) )
	   {
	     xLoc[1]=0;
	     if( ++xLoc[2] == GJP.NodeSites(2) )
	       {
		 xLoc[2]=0;
		 //if( ++xLoc[3] == GJP.NodeSites(3) )
		 if( ++xLoc[3] == localIndex_t )
		   {
		     printf("ERROR in QIO: with index/global-coor-conver\n");
		     exit(-1);
		   }
	       }
	   }
	}
    }
#endif
  int xLoc2[4]={0,0,0,0};
  int rest = index;
  for(int ii=0;ii<4;ii++){
	xLoc2[ii] = rest % GJP.NodeSites(ii);
	rest = rest/GJP.NodeSites(ii);
#ifdef USE_SLOW_INDEX
	if (xLoc2[ii] != xLoc[ii])
	ERR.General("",fname,"index=%d xLoc2[%d](%d) != xLoc[%d](%d)\n",
        index,ii,xLoc2[ii],ii,xLoc[ii]);
#endif
  }
  if( xLoc2[3] >= localIndex_t ) {
     printf("ERROR in QIO: with index/global-coor-conver\n");
     exit(-1);
   }
	
  for(int ii(0); ii < 4; ++ii)
    x[ii] += xLoc2[ii];
          
  //#define DEBUG_GetCorrds
  
  #ifdef DEBUG_GetCoords
  printf("UID: %i, called qio_get_coords with node: %i, index %i; returns x: [%i, %i, %i, %i]\n", UniqueID(), node, index, x[0], x[1], x[2], x[3]);
  #endif // DEBUG_GetCoords
}

#else
// To supporthe QIO_PARTFILE format for totnode != num io node
void qio_get_coords( int x[], int node, int index)
{
  const char *fname="qio_get_coords(i[],i,i)";
  #ifdef DEBUG_GetCoords
  printf("UID: %i, called qio_get_coords with node: %i, index %i; calc...\n", UniqueID(), node, index);
  #endif // DEBUG_GetCoords

  int node_coor[4];
  int nodes[5];
  nodes[0] = GJP.Xnodes();
  nodes[1] = GJP.Ynodes();
  nodes[2] = GJP.Znodes();
  nodes[3] = GJP.Tnodes();
  nodes[4] = GJP.Snodes();
  

  int nd=node;
  node_coor[0]= nd % nodes[0];  nd /= nodes[0];
  node_coor[1]= nd % nodes[1];  nd /= nodes[1];  
  node_coor[2]= nd % nodes[2];  nd /= nodes[2];    
  node_coor[3]= nd % nodes[3];  nd /= nodes[3];    
  
  // check the assumption, that the above code is actually 
  // calculating corrdinate correctly for the local node

  if( node == QMP_get_node_number() )
    { 
      if( GJP.XnodeCoor() != node_coor[0] ||
	  GJP.YnodeCoor() != node_coor[1] ||
	  GJP.ZnodeCoor() != node_coor[2] ||
	  GJP.TnodeCoor() != node_coor[3]) {
	printf("ERROR QIO: node number mismatch\n");
	exit(-1);
      }
    }
  
     
  x[0] = node_coor[0]*GJP.NodeSites(0);
  x[1] = node_coor[1]*GJP.NodeSites(1);
  x[2] = node_coor[2]*GJP.NodeSites(2);
  x[3] = node_coor[3]*GJP.NodeSites(3);

  int localIndex_t;

  if( GJP.Snodes() > 1 )
    {
      ERR.NotImplemented("","qio_get_coords(...)","many Ls is not supported for now, as I don't understand the QIO's logic fully, sorry");

      localIndex_t = GJP.NodeSites(3)/GJP.Snodes();
      x[3] += GJP.SnodeCoor()*GJP.NodeSites(3)/GJP.Snodes() ; 

    }
  else
    localIndex_t = GJP.NodeSites(3);

#undef USE_SLOW_INDEX
#ifdef USE_SLOW_INDEX
  int xLoc[4]={0,0,0,0};

  for( int ii(0); ii < index; ++ii)
    { 
      if( ++xLoc[0] == GJP.NodeSites(0) ) 
	{ 
	 xLoc[0]=0; 
	 if( ++xLoc[1] == GJP.NodeSites(1) )
	   {
	     xLoc[1]=0;
	     if( ++xLoc[2] == GJP.NodeSites(2) )
	       {
		 xLoc[2]=0;
		 //if( ++xLoc[3] == GJP.NodeSites(3) )
		 if( ++xLoc[3] == localIndex_t )
		   {
		     printf("ERROR in QIO: with index/global-coor-conver\n");
		     exit(-1);
		   }
	       }
	   }
	}
    }
#endif
  int xLoc2[4]={0,0,0,0};
  int rest = index;
  for(int ii=0;ii<4;ii++){
	xLoc2[ii] = rest % GJP.NodeSites(ii);
	rest = rest/GJP.NodeSites(ii);
#ifdef USE_SLOW_INDEX
	if (xLoc2[ii] != xLoc[ii])
	ERR.General("",fname,"index=%d xLoc2[%d](%d) != xLoc[%d](%d)\n",
        index,ii,xLoc2[ii],ii,xLoc[ii]);
#endif
  }
  if( xLoc2[3] >= localIndex_t ) {
     printf("ERROR in QIO: with index/global-coor-conver\n");
     exit(-1);
   }
	
  for(int ii(0); ii < 4; ++ii)
    x[ii] += xLoc2[ii];
          
  //#define DEBUG_GetCorrds
  
  #ifdef DEBUG_GetCoords
  printf("UID: %i, called qio_get_coords with node: %i, index %i; returns x: [%i, %i, %i, %i]\n", UniqueID(), node, index, x[0], x[1], x[2], x[3]);
  #endif // DEBUG_GetCoords
      }
#endif

int qio_num_sites(int node)
{
  // number of sites on node
  // this is the same on all nodes !!!

  if(GJP.Snodes() > 1 )
    {
      int tmp(1);
      for(int ii(0); ii < 4; ++ii)
	tmp *= GJP.NodeSites(ii);
      
      tmp /= GJP.Snodes();
      return tmp;

    }
  else
    return GJP.VolNodeSites();
}

int qio_sites_on_node()
{
  if (GJP.Snodes() > 1)
    {
      return ( ( GJP.NodeSites(0) * GJP.NodeSites(1) * GJP.NodeSites(2) * GJP.NodeSites(3) ) / GJP.Snodes() );
    }
  else
    return GJP.VolNodeSites();
}


#ifdef QIO_SPARSE_PARTFILE
static int qio_get_node(int node)
{
  //printf("qio_get_node: %d %d\n",node,QIO_SPARSE_PARTFILE_NODES);
  int QIO_SPARSE_PARTFILE_NODES=getQioSparseNum();
  if (QIO_SPARSE_PARTFILE_NODES==0) return node;
  return (node / QIO_SPARSE_PARTFILE_NODES)*QIO_SPARSE_PARTFILE_NODES;
}
static int qio_master_io_node(){ return 0; }
#endif

// now start the qio_init-class functions


//void qio_init::qio_setLayout( QIO_Layout *layout)
void qio_init::qio_setLayout()
{
  Float time = -dclock();
  static int lattice_size[QIO_RW_DIMENSION];

  
  lattice_size[0] = GJP.Xnodes()*GJP.XnodeSites();
  lattice_size[1] = GJP.Ynodes()*GJP.YnodeSites();
  lattice_size[2] = GJP.Znodes()*GJP.ZnodeSites();
  lattice_size[3] = GJP.Tnodes()*GJP.TnodeSites();
  

  // put the helper outside class...

  layout.node_number     = qio_node_number;
  layout.node_index      = qio_node_index;
  layout.get_coords      = qio_get_coords;
  layout.num_sites       = qio_num_sites;
  layout.latsize         = lattice_size;         // lattice size x,y,z,t
  layout.latdim          = QIO_RW_DIMENSION;     // lattice dimension for QIO
  //layout.volume          = GJP.VolSites();       // lattice volume
  layout.volume          = GJP.Sites(0) * GJP.Sites(1) * GJP.Sites(2) * GJP.Sites(3);
  //layout.sites_on_node   = GJP.VolNodeSites();   // local volume
  // SLOW ??? layout.sites_on_node   = ( GJP.NodeSites(0) * GJP.NodeSites(1) * GJP.NodeSites(2) * GJP.NodeSites(3) ) / GJP.Snodes();
  layout.sites_on_node   = qio_sites_on_node();
  layout.this_node       = QMP_get_node_number();// number of this node
  layout.number_of_nodes = NumNodes();           // total number of nodes

}

void qio_init::qio_setFilesystem()
{

#ifdef QIO_SPARSE_PARTFILE
  const char *fname="qio_setFilesystem()";
  VRB.Result("qio_init",fname,"using sparse partfile with sparse_num=%d\n",getQioSparseNum());
  fs. number_io_nodes = NumNodes() / getQioSparseNum();
  fs. type = QIO_SINGLE_PATH;// or QIO_MULTIPLE_PATH
  fs. my_io_node =  qio_get_node;
  fs. master_io_node = qio_master_io_node;

  fs. io_node = NULL;
  fs. node_path=NULL; // Provide array for path for each io-node, for QIO_MULTIPLE_PATH

  pointer_fs = &fs;

#else
  VRB.Result("qio_init",fname,"using deault filesystem");
  pointer_fs = NULL; //indicating the default filesystem setting
#endif

  
}




// general functions for reading global data


void qio_putGlobal(char *buf, size_t index, int count, void *arg)
{

#ifdef DEBUG_GlobalRead
  printf("PutGlobal: UID: %i: called with index: %i, count: %i\n",UniqueID(),index, count);
#endif // DEBUG_GlobalRead
    
  Float *data = (Float *) buf;
  Float *outp = (Float *) arg;

  for(int ii(0); ii < count; ++ ii)
    *(outp + ii) = *(data + ii);

#ifdef DEBUG_GlobalRead
  printf("PutGlobal: UID: %i: wrote %f %f to GlobalData\n",UniqueID(),*(outp), *(outp+1));
#endif // DEBUG_GlobalRead

}

void qio_putGlobalSingle(char *buf, size_t index, int count, void *arg)
{

#ifdef DEBUG_GlobalRead
  printf("PutGlobalSingle: UID: %i: called with index: %i, count: %i\n",UniqueID(),index, count);
#endif // DEBUG_GlobalRead
    
  float *data = (float *) buf;
  Float *outp = (Float *) arg;

  for(int ii(0); ii < count; ++ ii)
    *(outp + ii) = *(data + ii);

#ifdef DEBUG_GlobalRead
  printf("PutGlobal: UID: %i: wrote %f %f to GlobalData\n",UniqueID(),*(outp), *(outp+1));
#endif // DEBUG_GlobalRead

}

// writing global data


void qio_getGlobal(char *buf, size_t index, int count, void *arg)
{

#ifdef DEBUG_GlobalWrite
  printf("GetGlobal: UID: %i: called with index: %i, count: %i\n",UniqueID(),index, count);
#endif // DEBUG_GlobalWrite
    

  Float *data = (Float *) arg;
  Float *outp = (Float *) buf;

  for(int ii(0); ii < count; ++ ii)
    *(outp + ii) = *(data + ii);

#ifdef DEBUG_GlobalWrite
  printf("GetGlobal: UID: %i: wrote %f %f to GlobalData\n",UniqueID(),*(outp), *(outp+1));
#endif // DEBUG_GlobalWrite

}



void qio_getGlobalSingle(char *buf, size_t index, int count, void *arg)
{

#ifdef DEBUG_GlobalWrite
  printf("GetGlobalSingle: UID: %i: called with index: %i, count: %i\n",UniqueID(),index, count);
#endif // DEBUG_GlobalWrite
    

  Float *data = (Float *) arg;
  float *outp = (float *) buf;

  for(int ii(0); ii < count; ++ ii)
    *(outp + ii) = *(data + ii);

#ifdef DEBUG_GlobalWrite
  printf("GetGlobalSingle: UID: %i: wrote %f %f to GlobalData\n",UniqueID(),*(outp), *(outp+1));
#endif // DEBUG_GlobalWrite

}




 
CPS_END_NAMESPACE
#endif
