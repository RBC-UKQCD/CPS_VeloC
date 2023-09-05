/*!\file
  \brief   Global options for compiling the Colombia code:

  Generated automatically from config.h.in by configure procedure.

  $Id: config.h.in,v 1.29.28.1 2012-11-15 18:17:08 ckelly Exp $
*/
/* Global options for compiling the Columbia code:  
 * config.h.  Generated from config.h.in by configure.
 * 
 */

#ifndef INCLUDED_CONFIG_H_
#define INCLUDED_CONFIG_H_                  //!< Prevent multiple inclusion 

#include <conf.h>

#define NOARCH 0
#define QCDOC  1
#define QCDSP  2
#define BGL    3
#define BGP    4
#define BGQ    5


//#include <version.h>
#define VERSION_MAJOR 5
#define VERSION_MINOR 3
#define VERSION_SUB 0
#define VERSION_STR "CPS_V5.3.0"

#define TARGET NOARCH
#define PARALLEL 1

#undef HAVE_BFM

// The configure procedure should make this unnecessary, but just in case...
#ifndef TARGET
#define TARGET NOARCH
#endif

#define CPS_FLOAT_ALIGN

#define CWDPREFIX(A) A



/*! Explicit casting away of the const-ness  */
#define CAST_AWAY_CONST(x) ( const_cast<char*>(x) )

/*!  Precision in the global sum (undefined gives QCDSP behaviour). */
#define GLOBALSUM_TYPE double

#define CPS_END_NAMESPACE    }  
#define CPS_START_NAMESPACE  namespace cps {
#define USING_NAMESPACE_CPS  using namespace cps;
#define CPS_NAMESPACE	     cps


#define COMMS_DATASIZE (sizeof(double))

#undef UNIFORM_SEED_TESTING
#undef UNIFORM_SEED_NO_COMMS

/* ------------------------------------------------------------------*/

#endif /* INCLUDED_CONFIG_H_ */





