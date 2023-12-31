/*
 * Please do not edit this file.
 * It was generated using PAB's VML system.
 */

#ifndef _LANCZOS_ARG_H_RPCGEN
#define _LANCZOS_ARG_H_RPCGEN

#include <config.h>
#include <util/vml/types.h>
#include <util/vml/vml.h>
#include <util/enum.h>
#include <util/defines.h>
CPS_START_NAMESPACE

struct VML;
class MatrixPolynomialArg {
public:
	 bool Encode(char *filename,char *instance);
	 bool Decode(char *filename,char *instance);
	 bool Vml(VML *vmls,char *instance);
	int Npol;
	struct {
		u_int params_len;
		Float *params_val;
	} params;
	Float *tmp1;
	Float *tmp2;
};

struct VML;
class LanczosArg {
public:
	 bool Encode(char *filename,char *instance);
	 bool Decode(char *filename,char *instance);
	 bool Vml(VML *vmls,char *instance);
	enum RitzMatType RitzMat_lanczos;
	enum RitzMatType RitzMat_convcheck;
	Float mass;
	int nk_lanczos_vectors;
	int nt_lanczos_vectors;
	int np_lanczos_vectors;
	Float eigen_shift;
	Float stop_residual;
	int maxiters;
	int n_single;
	int save;
	int conv_check;
	PrecType precision;
	int mem_save;
	char *results;
	char *file;
	MatrixPolynomialArg matpoly_arg;
};

/* the xdr functions */

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__STDC__) || defined(__cplusplus)
extern  bool_t vml_MatrixPolynomialArg (VML *, char *instance, MatrixPolynomialArg*);
extern  bool_t vml_LanczosArg (VML *, char *instance, LanczosArg*);

#else /* K&R C */
extern  bool_t vml_MatrixPolynomialArg (VML *, char *instance, MatrixPolynomialArg*);
extern  bool_t vml_LanczosArg (VML *, char *instance, LanczosArg*);

#endif /* K&R C */

#ifdef __cplusplus
}
#endif
CPS_END_NAMESPACE

#endif /* !_LANCZOS_ARG_H_RPCGEN */
