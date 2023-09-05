/*
 * Please do not edit this file.
 * It was generated using PAB's VML system.
 */

#ifndef _ARRAY_ARG_H_RPCGEN
#define _ARRAY_ARG_H_RPCGEN

#include <config.h>
#include <util/vml/types.h>
#include <util/vml/vml.h>
#include <util/enum.h>
#include <util/defines.h>
CPS_START_NAMESPACE

struct VML;
class FloatArray {
public:
	 bool Encode(char *filename,char *instance);
	 bool Decode(char *filename,char *instance);
	 bool Vml(VML *vmls,char *instance);
	   void resize (  int n_floats ) ;
	   FloatArray (  ) ;
	struct {
		u_int Floats_len;
		Float *Floats_val;
	} Floats;
};

struct VML;
class ParamArg {
public:
	 bool Encode(char *filename,char *instance);
	 bool Decode(char *filename,char *instance);
	 bool Vml(VML *vmls,char *instance);
	char *name;
	Float val;
};

struct VML;
class ParamArray {
public:
	 bool Encode(char *filename,char *instance);
	 bool Decode(char *filename,char *instance);
	 bool Vml(VML *vmls,char *instance);
	   void resize (  u_int num ) ;
	   ParamArray (  ) ;
	struct {
		u_int params_len;
		ParamArg *params_val;
	} params;
};

struct VML;
class IntArray {
public:
	 bool Encode(char *filename,char *instance);
	 bool Decode(char *filename,char *instance);
	 bool Vml(VML *vmls,char *instance);
	struct {
		u_int v_len;
		int *v_val;
	} v;
};

/* the xdr functions */

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__STDC__) || defined(__cplusplus)
extern  bool_t vml_FloatArray (VML *, char *instance, FloatArray*);
extern  bool_t vml_ParamArg (VML *, char *instance, ParamArg*);
extern  bool_t vml_ParamArray (VML *, char *instance, ParamArray*);
extern  bool_t vml_IntArray (VML *, char *instance, IntArray*);

#else /* K&R C */
extern  bool_t vml_FloatArray (VML *, char *instance, FloatArray*);
extern  bool_t vml_ParamArg (VML *, char *instance, ParamArg*);
extern  bool_t vml_ParamArray (VML *, char *instance, ParamArray*);
extern  bool_t vml_IntArray (VML *, char *instance, IntArray*);

#endif /* K&R C */

#ifdef __cplusplus
}
#endif
CPS_END_NAMESPACE

#endif /* !_ARRAY_ARG_H_RPCGEN */