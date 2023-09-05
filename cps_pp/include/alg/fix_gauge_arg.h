/*
 * Please do not edit this file.
 * It was generated using PAB's VML system.
 */

#ifndef _FIX_GAUGE_ARG_H_RPCGEN
#define _FIX_GAUGE_ARG_H_RPCGEN

#include <config.h>
#include <util/vml/types.h>
#include <util/vml/vml.h>
#include <util/enum.h>
#include <util/defines.h>
CPS_START_NAMESPACE

struct VML;
class FixGaugeArg {
public:
	 bool Encode(char *filename,char *instance);
	 bool Decode(char *filename,char *instance);
	 bool Vml(VML *vmls,char *instance);
	FixGaugeType fix_gauge_kind;
	int hyperplane_start;
	int hyperplane_step;
	int hyperplane_num;
	Float stop_cond;
	int max_iter_num;
};

/* the xdr functions */

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__STDC__) || defined(__cplusplus)
extern  bool_t vml_FixGaugeArg (VML *, char *instance, FixGaugeArg*);

#else /* K&R C */
extern  bool_t vml_FixGaugeArg (VML *, char *instance, FixGaugeArg*);

#endif /* K&R C */

#ifdef __cplusplus
}
#endif
CPS_END_NAMESPACE

#endif /* !_FIX_GAUGE_ARG_H_RPCGEN */
