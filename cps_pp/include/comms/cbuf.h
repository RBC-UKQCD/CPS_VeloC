#include<config.h>

/*!\file
  \brief  Functions used by the data layout conversion functions.

  $Id: common.C,v 1.5 2004/09/02 16:59:26 zs Exp $
*/
#include <util/data_types.h>

CPS_START_NAMESPACE


//! Negate a floating point array.
/*!
  \param link_size The length of the array.
  \param link The array
  \post All elements of \a link are negated.
*/
  
inline void negate_link(unsigned link_size, IFloat *link)
{
	unsigned idx ;

	for (idx=0; idx<link_size; idx++) *(link+idx) = - *(link+idx) ;
}

//! Copy an array.
/*!
  \param src The array to copy from.
  \param dst The array to copy to.
  \param site_size The length of the array.
*/

inline void site2cram(IFloat *src, IFloat *dst, unsigned site_size)
{
	unsigned offset ;

	for (offset=0; offset<site_size; offset++)
		*(dst+offset) = *(src+offset) ;
}

//! Copy and rearrange an array.
/*!
  The copy is of the form <em>\n
  dst[link_tbl[i]] = src[i]
  \n</em>
  \param src The array to copy from.
  \param dst The array to copy to.
  \param link_tbl The look-up table for the \a dst indices.
  \param site_size The length of the \a src array.
*/

inline void site2dram(IFloat *src, IFloat *dst, unsigned *link_tbl, unsigned site_size)
{
	unsigned offset ;

	for (offset=0; offset<site_size; offset++)
		*(dst+*(link_tbl+offset)) = *(src+offset) ;
}





CPS_END_NAMESPACE
