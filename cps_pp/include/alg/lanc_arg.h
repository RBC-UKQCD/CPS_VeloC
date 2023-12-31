/*
 * Please do not edit this file.
 * It was generated using PAB's VML system.
 */

#ifndef _LANC_ARG_H_RPCGEN
#define _LANC_ARG_H_RPCGEN

#include <config.h>
#include <util/vml/types.h>
#include <util/vml/vml.h>
#include <util/enum.h>
#include <util/defines.h>
CPS_START_NAMESPACE

enum EigenType {
	D = 0,
	DDAG = 1,
	G5D = 2,
	DDAGD = 3,
};
typedef enum EigenType EigenType;
extern struct vml_enum_map EigenType_map[];


#include <util/vml/vml_templates.h>
class VML;
class LancArg {
public:
	 bool Encode(char *filename,char *instance);
	 bool Decode(char *filename,char *instance);
	 bool Vml(VML *vmls,char *instance);
	Float mass;
	Float stop_rsd;
	Float qr_rsd;
	enum EigenType EigenOper;
	bool_t precon;
	int N_get;
	int N_use;
	int N_true_get;
	int ch_ord;
	Float ch_alpha;
	Float ch_beta;
	bool_t ch_sh;
	Float ch_mu;
	bool_t lock;
	int maxits;
	char *fname;
	   void deep_copy(const LancArg &rhs);
	   void print(const std::string &prefix ="");
};
template<> struct rpc_deepcopy<LancArg>{
	static void doit(LancArg &into, LancArg const &from);
};

#ifndef _USE_STDLIB
#error "Cannot generate rpc_print commands without the standard library"
#endif
template<> struct rpc_print<LancArg>{
	static void doit(LancArg const &what, const std::string &prefix="" );
};


/* the xdr functions */

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__STDC__) || defined(__cplusplus)
extern  bool_t vml_EigenType (VML *, char *instance, EigenType*);
extern  bool_t vml_LancArg (VML *, char *instance, LancArg*);

#else /* K&R C */
extern  bool_t vml_EigenType (VML *, char *instance, EigenType*);
extern  bool_t vml_LancArg (VML *, char *instance, LancArg*);

#endif /* K&R C */

#ifdef __cplusplus
}
#endif
CPS_END_NAMESPACE

#endif /* !_LANC_ARG_H_RPCGEN */
