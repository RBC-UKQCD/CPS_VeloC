/*!\file
  \brief Definition of EqStateArg structure.

  $Id: eq_state_arg.h,v 1.3 2004/09/02 17:00:08 zs Exp $
*/

#include<config.h>
CPS_START_NAMESPACE

#ifndef INCLUDED_EQ_STATE_ARG_H
#define INCLUDED_EQ_STATE_ARG_H        //!< Prevent multiple inclusion.


//! Container of parameters for AlgEqState. 

struct EqStateArg {
  int dir;      /*!< The direction  determining what the hyperplanes
		  in which the plaquette is computed: The plaquette is
		  computed in the hyperplanes perpendicular and parallel
		  to this direction.	
		 */
};

#endif /* !INCLUDED_EQ_STATE_ARG_H */

CPS_END_NAMESPACE
