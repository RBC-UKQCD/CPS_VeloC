/*
 * Please do not edit this file.
 * It was generated using PAB's VML system.
 */

#ifndef _HYP_SMEAR_ARG_H_RPCGEN
#define _HYP_SMEAR_ARG_H_RPCGEN

#include <config.h>
#include <util/vml/types.h>
#include <util/vml/vml.h>
#include <util/enum.h>
#include <util/defines.h>
CPS_START_NAMESPACE

struct VML;
class HypSmearArg {
public:
	 bool Encode(char *filename,char *instance);
	 bool Decode(char *filename,char *instance);
	 bool Vml(VML *vmls,char *instance);
	Float tolerance;
	int orthog;
	Float c1;
	Float c2;
	Float c3;
	   HypSmearArg (  ) ;
};

/* the xdr functions */

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__STDC__) || defined(__cplusplus)
extern  bool_t vml_HypSmearArg (VML *, char *instance, HypSmearArg*);

#else /* K&R C */
extern  bool_t vml_HypSmearArg (VML *, char *instance, HypSmearArg*);

#endif /* K&R C */

#ifdef __cplusplus
}
#endif
CPS_END_NAMESPACE

#endif /* !_HYP_SMEAR_ARG_H_RPCGEN */
