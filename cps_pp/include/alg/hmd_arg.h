/*
 * Please do not edit this file.
 * It was generated using PAB's VML system.
 */

#ifndef _HMD_ARG_H_RPCGEN
#define _HMD_ARG_H_RPCGEN

#include <config.h>
#include <util/vml/types.h>
#include <util/vml/vml.h>
#include <util/enum.h>
#include <util/defines.h>
CPS_START_NAMESPACE

typedef Float FRatVec[MAX_RAT_DEGREE];

typedef int IMassVec[MAX_HMD_MASSES];

typedef Float FMassVec[MAX_HMD_MASSES];

typedef FRatVec FRatMassVec[MAX_HMD_MASSES];

struct VML;
class HmdArg {
public:
	 bool Encode(char *filename,char *instance);
	 bool Decode(char *filename,char *instance);
	 bool Vml(VML *vmls,char *instance);
	int n_frm_masses;
	int n_bsn_masses;
	int steps_per_traj;
	Float step_size;
	MetropolisType metropolis;
	ReunitarizeType reunitarize;
	RatApproxType approx_type;
	Float spread;
	long precision;
	int isz;
	int sw;
	int chrono;
	int reproduce;
	int reproduce_attempt_limit;
	FMassVec frm_mass;
	IMassVec frm_flavors;
	FMassVec bsn_mass;
	IMassVec max_num_iter;
	FMassVec stop_rsd;
	FMassVec stop_rsd_md;
	FMassVec stop_rsd_mc;
	FieldType field_type[MAX_HMD_MASSES];
	IMassVec valid_approx;
	FMassVec lambda_low;
	FMassVec lambda_high;
	FMassVec lambda_min;
	FMassVec lambda_max;
	RhmcPolesAction rhmc_poles_action;
	char *rhmc_poles_file;
	IMassVec frm_power_num;
	IMassVec frm_power_den;
	IMassVec FRatDeg;
	IMassVec FRatDegNew;
	IMassVec SRatDegNew;
	FMassVec FRatError;
	FMassVec FRatNorm;
	FRatMassVec FRatRes;
	FRatMassVec FRatPole;
	IMassVec SRatDeg;
	FMassVec SRatError;
	FMassVec SRatNorm;
	FRatMassVec SRatRes;
	FRatMassVec SRatPole;
	FMassVec SIRatNorm;
	FRatMassVec SIRatRes;
	FRatMassVec SIRatPole;
};

struct VML;
class EvoArg {
public:
	 bool Encode(char *filename,char *instance);
	 bool Decode(char *filename,char *instance);
	 bool Vml(VML *vmls,char *instance);
	int traj_start;
	int gauge_unload_period;
	int gauge_configurations;
	int io_concurrency;
	int hdw_xcsum;
	int hdw_rcsum;
	int reproduce_interval;
	char *ensemble_id;
	char *ensemble_label;
	char *creator;
	char *gauge_file_stem;
	char *rng_file_stem;
	char *plaquette_stem;
	char *pbp_stem;
	char *evo_stem;
	char *w_spect_directory;
	char *work_directory;
	int measure_pbp;
	int measure_w_spect_interval;
};

struct VML;
class RhmcPolesState {
public:
	 bool Encode(char *filename,char *instance);
	 bool Decode(char *filename,char *instance);
	 bool Vml(VML *vmls,char *instance);
	IMassVec frm_power_num;
	IMassVec frm_power_den;
	IMassVec FRatDeg;
	IMassVec FRatDegNew;
	IMassVec SRatDegNew;
	FMassVec FRatError;
	FMassVec FRatNorm;
	FRatMassVec FRatRes;
	FRatMassVec FRatPole;
	IMassVec SRatDeg;
	FMassVec SRatError;
	FMassVec SRatNorm;
	FRatMassVec SRatRes;
	FRatMassVec SRatPole;
	FMassVec SIRatNorm;
	FRatMassVec SIRatRes;
	FRatMassVec SIRatPole;
};

/* the xdr functions */

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__STDC__) || defined(__cplusplus)
extern  bool_t vml_FRatVec (VML *, char *instance, FRatVec);
extern  bool_t vml_IMassVec (VML *, char *instance, IMassVec);
extern  bool_t vml_FMassVec (VML *, char *instance, FMassVec);
extern  bool_t vml_FRatMassVec (VML *, char *instance, FRatMassVec);
extern  bool_t vml_HmdArg (VML *, char *instance, HmdArg*);
extern  bool_t vml_EvoArg (VML *, char *instance, EvoArg*);
extern  bool_t vml_RhmcPolesState (VML *, char *instance, RhmcPolesState*);

#else /* K&R C */
extern  bool_t vml_FRatVec (VML *, char *instance, FRatVec);
extern  bool_t vml_IMassVec (VML *, char *instance, IMassVec);
extern  bool_t vml_FMassVec (VML *, char *instance, FMassVec);
extern  bool_t vml_FRatMassVec (VML *, char *instance, FRatMassVec);
extern  bool_t vml_HmdArg (VML *, char *instance, HmdArg*);
extern  bool_t vml_EvoArg (VML *, char *instance, EvoArg*);
extern  bool_t vml_RhmcPolesState (VML *, char *instance, RhmcPolesState*);

#endif /* K&R C */

#ifdef __cplusplus
}
#endif
CPS_END_NAMESPACE

#endif /* !_HMD_ARG_H_RPCGEN */
