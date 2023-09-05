#include "algklks.h"

void AlgKlKs::RunPion(const char* label) const
{
	const char* fname = "RunPion()";	
	//static Timer time(fname);
	//time.start();
	double cost = -dclock();
	
	WilsonMatrix prop1, prop2;
	
	//pion propagators
	Rcomplex * pion = new Rcomplex [glb_t];
	for(int i = 0; i < glb_t; i++) pion[i] = 0.0;
	for(int tsrc = 0; tsrc < glb_t; tsrc++)
		for(int tsnk = 0; tsnk < glb_t; tsnk++)
		{
			prop1 = lwwprop[tsrc][tsnk];
			prop2.hconj(prop1);
			int delta = (tsnk+glb_t-tsrc)%glb_t;
			pion[delta] += Trace(prop1, prop2) * (1.0/ glb_t);
		}

	//output result
	char filename[512];
	sprintf(filename, "%s_pion_%s.txt", common_arg->filename, label);
	SaveCorr(pion, filename, glb_t);

	delete [] pion;

	cost += dclock();
	print_time(cname, fname, cost);
//	time.stop();
}

void AlgKlKs::Runpipi(const char* label) const {
	const char *fname="Runpipi";
//	static Timer time(fname);
//	time.start();
	double cost=-dclock();
	WilsonMatrix Mtemp1,Mtemp2,Mtemp3;
	Rcomplex pipi[4][glb_t];
	Rcomplex ptov[glb_t];
	Rcomplex disc[2][glb_t];//disconnected piece 
	memset(pipi,0,sizeof(Rcomplex)*4*glb_t);
	for(int tsrc=0;tsrc<glb_t;tsrc++){
		////////////////pion to vaccum
		int tsrc2=(tsrc+glb_t-sep)%glb_t;
		ptov[tsrc]=Trace(lwwprop[tsrc][tsrc2],lwwprop[tsrc][tsrc2].conj_cp());
		disc[0][tsrc]=lwwprop[tsrc][tsrc].Trace();
		disc[1][tsrc]=(lwwprop[tsrc][tsrc].gl(-5)).Trace();
		for(int tsink=0;tsink<glb_t;tsink++)
		{
			////////////////////D term
			int tsink2=(tsink+glb_t+sep)%glb_t;
			int diff=(tsink-tsrc+glb_t)%glb_t;
			Mtemp1=lwwprop[tsrc][tsink];
			Mtemp1.hconj();
			Rcomplex Rtemp1=Trace(Mtemp1, lwwprop[tsrc][tsink]);

			Mtemp2=lwwprop[tsrc2][tsink2];
			Mtemp2.hconj();
			Rcomplex Rtemp2=Trace(Mtemp2, lwwprop[tsrc2][tsink2]);

			pipi[0][diff]+=0.5*Rtemp1*Rtemp2;//D term

			Mtemp1=lwwprop[tsrc][tsink2];
			Mtemp1.hconj();
			Rtemp1=Trace(Mtemp1, lwwprop[tsrc][tsink2]);

			Mtemp2=lwwprop[tsrc2][tsink];
			Mtemp2.hconj();
			Rtemp2=Trace(Mtemp2, lwwprop[tsrc2][tsink]);

			pipi[0][diff]+= 0.5*Rtemp1*Rtemp2;
			///////////////////C term 
			Mtemp1=lwwprop[tsrc][tsink];
			Mtemp1.hconj();
			Mtemp2=lwwprop[tsrc2][tsink2];
			Mtemp2.hconj();
			pipi[1][diff]=+Trace(lwwprop[tsrc][tsink2]*Mtemp1, lwwprop[tsrc2][tsink]*Mtemp2);

			/////////////////R term

			Mtemp1=lwwprop[tsrc][tsrc2];
			Mtemp1.hconj();
			Mtemp2=lwwprop[tsink2][tsink];
			Mtemp2.hconj();
			pipi[2][diff]+=0.5*Trace(lwwprop[tsrc][tsink]*Mtemp1, lwwprop[tsink2][tsrc2]*Mtemp2);
			Mtemp1=lwwprop[tsrc][tsrc2];
			Mtemp1.hconj();
			Mtemp2=lwwprop[tsink][tsink2];
			Mtemp2.hconj();
			pipi[2][diff] += 0.5*Trace(lwwprop[tsrc][tsink2]*Mtemp1, lwwprop[tsink][tsrc2]*Mtemp2);
			/////////////////V term

			Mtemp1=lwwprop[tsrc][tsrc2];
			Mtemp1.hconj();
			Mtemp2=lwwprop[tsink][tsink2];
			Mtemp2.hconj();
			pipi[3][diff]+=Trace(lwwprop[tsrc][tsrc2],Mtemp1)*Trace(lwwprop[tsink][tsink2],Mtemp2);
		}//end t sink
	}//end t source

	for(int i=0;i<4;i++)
		for(int j=0;j<glb_t;j++)
			pipi[i][j]=pipi[i][j] * (1.0/glb_t);

	char filename[512];
	sprintf(filename, "%s_pipi_%s.txt", common_arg->filename, label);
	SaveCorr(pipi[0], filename, glb_t*4);

	char filename2[512];
	sprintf(filename2, "%s_ptov_%s.txt",common_arg->filename, label);
	SaveCorr(ptov, filename2, glb_t);

	char filename3[512];
	sprintf(filename3,"%s_disc_%s.txt",common_arg->filename, label);
	SaveCorr(disc[0], filename3, glb_t);

	cost += dclock();
	print_time(cname, fname, cost);
//	time.stop();
}



void AlgKlKs::RunKaon(const char* label) const
{
	const char *fname = "RunKaon()";
//	static Timer time(fname);
//	time.start();
	double cost = -dclock();

	WilsonMatrix prop1, prop2, prop3;

	//kaon propagators
	Rcomplex * kaon = new Rcomplex [glb_t];
	Rcomplex * ssbar = new Rcomplex [glb_t];
	Rcomplex * dll = new Rcomplex [glb_t];
	Rcomplex * dss = new Rcomplex [glb_t];
	Rcomplex * dls = new Rcomplex [glb_t];
	for(int t = 0; t < glb_t; t++)
		kaon[t] = ssbar[t] = dll[t] = dss[t] = dls[t] = 0.0;
	for(int tsrc = 0; tsrc < glb_t; tsrc++)
		for(int tsnk = 0; tsnk < glb_t; tsnk++)
		{
			prop1 = lwwprop[tsrc][tsnk];
			prop2 = swwprop[tsrc][tsnk];
			prop3.hconj(prop2);
			int delta = (tsnk+glb_t-tsrc)%glb_t;
			kaon[delta] += Trace(prop1, prop3) * (1.0/ glb_t);
			ssbar[delta] += Trace(prop2, prop3) *(1.0/ glb_t);

			Rcomplex lsrc, lsnk, ssrc, ssnk;
			prop1.glV(lwwprop[tsrc][tsrc],-5);
			lsrc = prop1.Trace();

			prop1.glV(lwwprop[tsnk][tsnk],-5);
			lsnk = prop1.Trace();

			prop1.glV(swwprop[tsrc][tsrc],-5);
			ssrc = prop1.Trace();

			prop1.glV(swwprop[tsnk][tsnk],-5);
			ssnk = prop1.Trace();

			dll[delta] += lsrc*lsnk * (1.0/ glb_t);
			dss[delta] += ssrc*ssnk * (1.0/ glb_t);
			dls[delta] += (lsrc*ssnk+ssrc*lsnk) / (glb_t * 2.0);
		}

	//output result
	char filename[512];
	sprintf(filename, "%s_kaon_%s.txt", common_arg->filename, label);
	SaveCorr(kaon, filename, glb_t);

	sprintf(filename, "%s_ssbar_%s.txt", common_arg->filename, label);
	SaveCorr(ssbar, filename, glb_t);

	sprintf(filename, "%s_dll_%s.txt", common_arg->filename, label);
	SaveCorr(dll, filename, glb_t);

	sprintf(filename, "%s_dss_%s.txt", common_arg->filename, label);
	SaveCorr(dss, filename, glb_t);

	sprintf(filename, "%s_dls_%s.txt", common_arg->filename, label);
	SaveCorr(dls, filename, glb_t);

	delete [] kaon;
	delete [] ssbar;
	delete [] dll;
	delete [] dss;
	delete [] dls;

	cost += dclock();
	print_time(cname, fname, cost);
//	time.stop();
}
