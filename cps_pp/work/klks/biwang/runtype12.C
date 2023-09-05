#include "algklks.h"

void AlgKlKs::RunType12(const char* label) const
{
	const char *fname = "RunType12()";
	Float cost = -dclock();

	//const int N_threads_old = bfmarg::default_threads;
	int N_threads = 64;
	omp_set_num_threads(64);

	//wall source point sink propagator, 
	//sink point located at (x,y,z), t from ta to tb
	WilsonMatrix ltip[glb_t][T], stip[glb_t][T], ltfp[glb_t][T], stfp[glb_t][T];

	for(int tk = 0; tk < glb_t; tk++)
	{
		int ta = (tk+minsep)%glb_t;
		int tb = (tk+minsep+boxsize)%glb_t;
		int tkb = (tk+2*minsep+boxsize)%glb_t; 
		WilsonMatrix tmp;
		for(int t = 0; t < GJP.TnodeSites(); t++)
		{ 
			int physT = t + toffset;
			if(tb>=ta && (physT<ta || physT>tb)) continue;
			if(tb<ta && physT<ta && physT>tb) continue;

			int xsrc_t = (4*physT)%GJP.Sites(0);
			int ysrc_t = (4*physT)%GJP.Sites(1);
			int zsrc_t = (4*physT)%GJP.Sites(2);

			int nx = xsrc_t / GJP.XnodeSites();
			int ny = ysrc_t / GJP.YnodeSites();
			int nz = zsrc_t / GJP.ZnodeSites();
			int x = xsrc_t % GJP.XnodeSites();
			int y = ysrc_t % GJP.YnodeSites();
			int z = zsrc_t % GJP.ZnodeSites();

			Site s(x,y,z,t);
			int ind = s.Index();
			if(GJP.XnodeCoor() == nx && GJP.YnodeCoor() == ny && GJP.ZnodeCoor() == nz)
			{
				int tind = (physT+glb_t-ta)%glb_t;
				ltip[tk][tind]  = (*lwall[tk])[ind];
				ltfp[tk][tind]  = (*lwall[tkb])[ind];

				tmp.hconj((*swall[tk])[ind]);
				stip[tk][tind].grV(tmp, -5);

				tmp.hconj((*swall[tkb])[ind]);
				stfp[tk][tind].grV(tmp,-5);
			}
		}
	}

#ifdef PARALLEL
	QMP_sum_double_array((double*)ltip, 288*T*glb_t);
	QMP_sum_double_array((double*)ltfp, 288*T*glb_t);
	QMP_sum_double_array((double*)stip, 288*T*glb_t);
	QMP_sum_double_array((double*)stfp, 288*T*glb_t);
#endif

	Rcomplex dspVtrace[glb_t][4][T], dspAtrace[glb_t][4][T];   
	Rcomplex dspVtrace2[glb_t][4][T], dspAtrace2[glb_t][4][T]; 
#pragma omp parallel for 
	for(int ind=0; ind < 4*T; ind++)
	{
		int mu = ind % 4;
		int time = ind / 4;
		WilsonMatrix tmp;
		for(int tk = 0; tk < glb_t; tk++)
		{
			dspVtrace[tk][mu][time] = Trace(tmp.glV(ltip[tk][time], mu), stip[tk][time]);
			dspAtrace[tk][mu][time] = Trace(tmp.glA(ltip[tk][time], mu), stip[tk][time]);

			dspVtrace2[tk][mu][time] = Trace(tmp.glV(ltfp[tk][time], mu), stfp[tk][time]);
			dspAtrace2[tk][mu][time] = Trace(tmp.glA(ltfp[tk][time], mu), stfp[tk][time]);
		}
	}

	Rcomplex dsVtrace[glb_t][4][size_4d], dsAtrace[glb_t][4][size_4d];	
	Rcomplex dsVtrace2[glb_t][4][size_4d], dsAtrace2[glb_t][4][size_4d];	
#pragma omp parallel for
	for(int ind = 0; ind < 4*size_4d; ind++)
	{
		int mu = ind % 4;
		int site = ind / 4;
		WilsonMatrix tmp1, tmp2;
		for(int tk = 0; tk < glb_t; tk++)
		{
			int tkb = (tk+2*minsep+boxsize)%glb_t; 
			tmp1.hconj((*swall[tkb])[site]);
			tmp2.grV(tmp1, -5);
			dsVtrace[tk][mu][site] = Trace(tmp1.glV((*lwall[tkb])[site], mu), tmp2); 
			dsAtrace[tk][mu][site] = Trace(tmp1.glA((*lwall[tkb])[site], mu), tmp2); 

			tmp1.hconj((*swall[tk])[site]);
			tmp2.grV(tmp1, -5);
			dsVtrace2[tk][mu][site] = Trace(tmp1.glV((*lwall[tk])[site], mu), tmp2); 
			dsAtrace2[tk][mu][site] = Trace(tmp1.glA((*lwall[tk])[site], mu), tmp2);
		}
	}

	WilsonMatrix dsp[glb_t][T], dsp2[glb_t][T];
#pragma omp parallel for
	for(int tk = 0; tk < glb_t; tk++)
		for(int time=0; time<T; time++) 
		{
			eq_mult(dsp[tk][time], ltip[tk][time], stip[tk][time]);
			eq_mult(dsp2[tk][time], ltfp[tk][time], stfp[tk][time]);
		}

	Rcomplex _pc[8][T*T][N_threads], _pv[8][T*T][N_threads];
#pragma omp parallel for schedule(dynamic)
	for(int ind = 0; ind < T * size_4d; ind++)
	{
		int tx = ind / size_4d;
		int site = ind % size_4d;
		int physT = site / size_3d + toffset;
		int tid = omp_get_thread_num();
		for(int tk = 0; tk < glb_t; tk++)
		{
			int ta = (tk+minsep)%glb_t;
			int tb = (tk+minsep+boxsize)%glb_t;
			int tkb = (tk+2*minsep+boxsize)%glb_t; 
			if(tb>=ta && (physT<ta || physT>tb)) continue; 
			if(tb<ta && physT<ta && physT>tb) continue;
			int ty = (physT-ta+glb_t)%glb_t;

			WilsonMatrix aux1, aux2;
			WilsonMatrix tmp1, tmp2;
			WilsonMatrix uxy, uyx;
			WilsonMatrix lti, sti, ltf, stf;

			//all elements in the contractions
			lti = (*lwall[tk])[site];
			ltf = (*lwall[tkb])[site];
			tmp1.hconj((*swall[tk])[site]);
			sti.grV(tmp1, -5);			
			tmp1.hconj((*swall[tkb])[site]);
			stf.grV(tmp1, -5);
			uxy = (*lmcpoint[(tx+ta)%glb_t])[site];
			tmp1.hconj(uxy);
			tmp2.glV(tmp1, -5);
			uyx.grV(tmp2, -5);

			//all products of two internal quark propagators
			WilsonMatrix uV1[4], uA1[4], uV2[4], uA2[4]; 
			WilsonMatrix uuV1[4], uuA1[4], uuV2[4], uuA2[4];  
			for(int mu = 0; mu < 4; mu++)
			{
				//up-up quark product
				uV1[mu].glV(uxy, mu);
				eq_mult(uuV1[mu], uyx, uV1[mu]);	
				uA1[mu].glA(uxy, mu);
				eq_mult(uuA1[mu], uyx, uA1[mu]);	
				uV2[mu].glV(uyx, mu);
				eq_mult(uuV2[mu], uxy, uV2[mu]);	
				uA2[mu].glA(uyx, mu);
				eq_mult(uuA2[mu], uxy, uA2[mu]);	
			}

			Rcomplex uuVV[4][4], uuVA[4][4], uuAV[4][4], uuAA[4][4];
			for(int mu=0; mu < 4; mu++)
				for(int nu=0; nu<4; nu++)
				{
					uuVV[mu][nu] = tmp1.glV(uuV2[mu], nu).Trace();
					uuVA[mu][nu] = tmp1.glA(uuV2[mu], nu).Trace();
					uuAV[mu][nu] = tmp1.glV(uuA2[mu], nu).Trace();
					uuAA[mu][nu] = tmp1.glA(uuA2[mu], nu).Trace();
				}

			//fix \vec{x} = 0, sum all possible \vec{y}
			int index = tx*T + ty;
			//contraction 1
			for(int mu=0; mu < 4; mu++) 
				for(int nu=0; nu<4; nu++)
				{
					_pc[0][index][tid] += 
						dspVtrace[tk][mu][tx] * uuVV[mu][nu] * dsVtrace[tk][nu][site] + 
						dspVtrace[tk][mu][tx] * uuVA[mu][nu] * dsAtrace[tk][nu][site] +             
						dspAtrace[tk][mu][tx] * uuAV[mu][nu] * dsVtrace[tk][nu][site] +             
						dspAtrace[tk][mu][tx] * uuAA[mu][nu] * dsAtrace[tk][nu][site]; 

					_pv[0][index][tid] += 
						dspVtrace[tk][mu][tx] * uuAV[mu][nu] * dsAtrace[tk][nu][site] + 
						dspVtrace[tk][mu][tx] * uuAA[mu][nu] * dsVtrace[tk][nu][site] +             
						dspAtrace[tk][mu][tx] * uuVV[mu][nu] * dsAtrace[tk][nu][site] +             
						dspAtrace[tk][mu][tx] * uuVA[mu][nu] * dsVtrace[tk][nu][site]; 
				}

			//contraction 2
			WilsonMatrix ds;
			eq_mult(ds, ltf, stf);
			WilsonMatrix auxpCon(0.0), auxpVio(0.0), auxCon(0.0), auxVio(0.0);

			for(int mu = 0; mu < 4; mu++) 
			{
				tmp1.glV(dsp[tk][tx], mu);
				tmp2.grV(tmp1, mu);
				auxpCon += tmp2;
				tmp2.grA(tmp1, mu);
				auxpVio += tmp2;

				tmp1.glA(dsp[tk][tx], mu);
				tmp2.grA(tmp1, mu);
				auxpCon += tmp2;
				tmp2.grV(tmp1, mu);
				auxpVio += tmp2;

				tmp1.glV(ds, mu);
				tmp2.grV(tmp1, mu);
				auxCon += tmp2;
				tmp2.grA(tmp1, mu);
				auxVio += tmp2;

				tmp1.glA(ds, mu);
				tmp2.grA(tmp1, mu);
				auxCon += tmp2;
				tmp2.grV(tmp1, mu);
				auxVio += tmp2;
			}
			eq_mult(aux1, auxpCon, uyx); 
			eq_mult(aux2, auxCon, uxy);
			_pc[1][index][tid] += Trace(aux1, aux2);

			eq_mult(aux1, auxpVio, uyx);	
			eq_mult(aux2, auxVio, uxy);
			_pv[1][index][tid] += Trace(aux1, aux2);

			//contraction 3 and 4
			for(int mu = 0; mu < 4; mu++)
			{
				_pc[2][index][tid] += 
					Trace(auxpCon, uuV1[mu]) * dsVtrace[tk][mu][site] + 
					Trace(auxpCon, uuA1[mu]) * dsAtrace[tk][mu][site];

				_pv[2][index][tid] += 
					Trace(auxpVio, uuV1[mu]) * dsAtrace[tk][mu][site] +
					Trace(auxpVio, uuA1[mu]) * dsVtrace[tk][mu][site];

				_pc[3][index][tid] +=
					Trace(auxCon, uuV2[mu]) * dspVtrace[tk][mu][tx] +
					Trace(auxCon, uuA2[mu]) * dspAtrace[tk][mu][tx];

				_pv[3][index][tid] += 
					Trace(auxVio, uuV2[mu]) * dspAtrace[tk][mu][tx] +
					Trace(auxVio, uuA2[mu]) * dspVtrace[tk][mu][tx];
			}

			//contraction 5
			eq_mult(aux1, ltip[tk][tx], sti); 
			eq_mult(aux2, ltf, stfp[tk][tx]);

			WilsonMatrix aux1V[4], aux1A[4], aux2V[4], aux2A[4];
			for(int mu = 0; mu < 4; mu++)
			{
				aux1V[mu].glV(aux1, mu); 
				aux1A[mu].glA(aux1, mu); 
				aux2V[mu].glV(aux2, mu); 
				aux2A[mu].glA(aux2, mu); 
			}	
			Rcomplex dsdsVV, dsdsVA, dsdsAV, dsdsAA;

			for(int mu = 0; mu < 4; mu++) 
				for( int nu = 0; nu < 4; nu++)
				{
					dsdsVV = Trace(aux1V[mu], aux2V[nu]);
					dsdsVA = Trace(aux1V[mu], aux2A[nu]);
					dsdsAV = Trace(aux1A[mu], aux2V[nu]);
					dsdsAA = Trace(aux1A[mu], aux2A[nu]);

					_pc[4][index][tid] += 
						dsdsVV * uuVV[mu][nu] +
						dsdsVA * uuVA[mu][nu] +
						dsdsAV * uuAV[mu][nu] + 
						dsdsAA * uuAA[mu][nu];
					_pv[4][index][tid] +=
						dsdsVV * uuAA[mu][nu] +
						dsdsVA * uuAV[mu][nu] +
						dsdsAV * uuVA[mu][nu] +
						dsdsAA * uuVV[mu][nu];	
				}

			//contraction 6
			Rcomplex dsu1VV, dsu1VA, dsu1AV, dsu1AA, dsu2VV, dsu2VA, dsu2AV, dsu2AA;
			for(int mu = 0; mu < 4; mu++) 
				for(int nu =0; nu < 4; nu++)
				{
					dsu1VV = Trace(aux1V[mu], uV1[nu]);
					dsu1VA = Trace(aux1V[mu], uA1[nu]);
					dsu1AV = Trace(aux1A[mu], uV1[nu]);
					dsu1AA = Trace(aux1A[mu], uA1[nu]);

					dsu2VV = Trace(uV2[mu], aux2V[nu]);
					dsu2VA = Trace(uV2[mu], aux2A[nu]);
					dsu2AV = Trace(uA2[mu], aux2V[nu]);
					dsu2AA = Trace(uA2[mu], aux2A[nu]);

					_pc[5][index][tid] +=
						dsu1VV * dsu2VV + dsu1AA * dsu2AA +
						dsu1VA * dsu2VA + dsu1AV * dsu2AV;
					_pv[5][index][tid] +=
						dsu1VV * dsu2AA + dsu1VA * dsu2AV +
						dsu1AV * dsu2VA + dsu1AA * dsu2VV;				
				}

			//contraction 7 and 8
			WilsonMatrix dsds1V[4], dsds1A[4], dsds2V[4], dsds2A[4];
			for(int mu = 0; mu < 4; mu++)
			{
				eq_mult(dsds1V[mu], aux1, aux2V[mu]);
				eq_mult(dsds1A[mu], aux1, aux2A[mu]);
				eq_mult(dsds2V[mu], aux2, aux1V[mu]);
				eq_mult(dsds2A[mu], aux2, aux1A[mu]);
			}

			WilsonMatrix dsds1ConV[4], dsds1VioV[4], dsds1ConA[4], dsds1VioA[4];
			WilsonMatrix dsds2ConV[4], dsds2VioV[4], dsds2ConA[4], dsds2VioA[4];	
			for(int mu = 0; mu < 4; mu++) 
			{
				dsds1ConV[mu] = dsds1VioV[mu] = dsds1ConA[mu] = dsds1VioA[mu] = 0;
				dsds2ConV[mu] = dsds2VioV[mu] = dsds2ConA[mu] = dsds2VioA[mu] = 0;
			}
			for(int mu = 0; mu < 4; mu++) 
				for(int nu = 0; nu < 4; nu++)
				{
					tmp1.glV(dsds1V[mu], nu);
					tmp2.grV(tmp1, nu);
					dsds1ConV[mu] += tmp2;
					tmp2.grA(tmp1, nu);
					dsds1VioV[mu] += tmp2;
					tmp1.glA(dsds1V[mu], nu);
					tmp2.grA(tmp1, nu);
					dsds1ConV[mu] += tmp2;
					tmp2.grV(tmp1, nu);
					dsds1VioV[mu] += tmp2;

					tmp1.glV(dsds1A[mu], nu);
					tmp2.grV(tmp1, nu);
					dsds1ConA[mu] += tmp2;
					tmp2.grA(tmp1, nu);
					dsds1VioA[mu] += tmp2;
					tmp1.glA(dsds1A[mu], nu);
					tmp2.grA(tmp1, nu);
					dsds1ConA[mu] += tmp2;
					tmp2.grV(tmp1, nu);
					dsds1VioA[mu] += tmp2;

					tmp1.glV(dsds2V[mu], nu);
					tmp2.grV(tmp1, nu);
					dsds2ConV[mu] += tmp2;
					tmp2.grA(tmp1, nu);
					dsds2VioV[mu] += tmp2;
					tmp1.glA(dsds2V[mu], nu);
					tmp2.grA(tmp1, nu);
					dsds2ConV[mu] += tmp2;
					tmp2.grV(tmp1, nu);
					dsds2VioV[mu] += tmp2;

					tmp1.glV(dsds2A[mu], nu);
					tmp2.grV(tmp1, nu);
					dsds2ConA[mu] += tmp2;
					tmp2.grA(tmp1, nu);
					dsds2VioA[mu] += tmp2;
					tmp1.glA(dsds2A[mu], nu);
					tmp2.grA(tmp1, nu);
					dsds2ConA[mu] += tmp2;
					tmp2.grV(tmp1, nu);
					dsds2VioA[mu] += tmp2;
				}

			for(int mu = 0; mu < 4; mu++) 
			{
				_pc[6][index][tid] += 
					Trace(dsds1ConV[mu], uuV1[mu]) +
					Trace(dsds1ConA[mu], uuA1[mu]);
				_pv[6][index][tid] +=
					Trace(dsds1VioV[mu], uuA1[mu]) +
					Trace(dsds1VioA[mu], uuV1[mu]);
				_pc[7][index][tid] +=
					Trace(dsds2ConV[mu], uuV2[mu]) +
					Trace(dsds2ConA[mu], uuA2[mu]);
				_pv[7][index][tid] +=
					Trace(dsds2VioV[mu], uuA2[mu]) +
					Trace(dsds2VioA[mu], uuV2[mu]);

			}

			/*********************************************************************************************************/			
			//fix \vec{y} = 0, sum all possible \vec{x}
			index = ty*T + tx;
			//contraction 1
			for(int mu=0; mu < 4; mu++) 
				for(int nu=0; nu<4; nu++)
				{
					_pc[0][index][tid] += dspVtrace2[tk][mu][tx] * uuVV[mu][nu] * dsVtrace2[tk][nu][site] 
						+dspVtrace2[tk][mu][tx] * uuVA[mu][nu] * dsAtrace2[tk][nu][site]               
						+dspAtrace2[tk][mu][tx] * uuAV[mu][nu] * dsVtrace2[tk][nu][site]               
						+dspAtrace2[tk][mu][tx] * uuAA[mu][nu] * dsAtrace2[tk][nu][site]; 

					_pv[0][index][tid] += dspVtrace2[tk][mu][tx] * uuAV[mu][nu] * dsAtrace2[tk][nu][site] 
						+dspVtrace2[tk][mu][tx] * uuAA[mu][nu] * dsVtrace2[tk][nu][site]               
						+dspAtrace2[tk][mu][tx] * uuVV[mu][nu] * dsAtrace2[tk][nu][site]               
						+dspAtrace2[tk][mu][tx] * uuVA[mu][nu] * dsVtrace2[tk][nu][site]; 
				}

			//contraction 2
			eq_mult(ds, lti, sti);
			auxpCon = auxpVio = auxCon = auxVio = 0.0;
			for(int mu = 0; mu < 4; mu++) 
			{
				tmp1.glV(dsp2[tk][tx], mu);
				tmp2.grV(tmp1, mu);
				auxpCon += tmp2;
				tmp2.grA(tmp1, mu);
				auxpVio += tmp2;

				tmp1.glA(dsp2[tk][tx], mu);
				tmp2.grA(tmp1, mu);
				auxpCon += tmp2;
				tmp2.grV(tmp1, mu);
				auxpVio += tmp2;

				tmp1.glV(ds, mu);
				tmp2.grV(tmp1, mu);
				auxCon += tmp2;
				tmp2.grA(tmp1, mu);
				auxVio += tmp2;

				tmp1.glA(ds, mu);
				tmp2.grA(tmp1, mu);
				auxCon += tmp2;
				tmp2.grV(tmp1, mu);
				auxVio += tmp2;
			}
			eq_mult(aux1, auxpCon, uyx); 
			eq_mult(aux2, auxCon, uxy);
			_pc[1][index][tid] += Trace(aux1, aux2);

			eq_mult(aux1, auxpVio, uyx);	
			eq_mult(aux2, auxVio, uxy);
			_pv[1][index][tid] += Trace(aux1, aux2);

			//contraction 3 and 4
			for(int mu = 0; mu < 4; mu++)
			{
				_pc[3][index][tid] +=
					Trace(auxpCon, uuV1[mu]) * dsVtrace2[tk][mu][site] +
					Trace(auxpCon, uuA1[mu]) * dsAtrace2[tk][mu][site];
				_pv[3][index][tid] += 
					Trace(auxpVio, uuV1[mu]) * dsAtrace2[tk][mu][site] + 
					Trace(auxpVio, uuA1[mu]) * dsVtrace2[tk][mu][site];

				_pc[2][index][tid] +=
					Trace(auxCon, uuV2[mu]) * dspVtrace2[tk][mu][tx] +
					Trace(auxCon, uuA2[mu]) * dspAtrace2[tk][mu][tx];
				_pv[2][index][tid] +=
					Trace(auxVio, uuV2[mu]) * dspAtrace2[tk][mu][tx] +
					Trace(auxVio, uuA2[mu]) * dspVtrace2[tk][mu][tx];	
			}

			//contraction 5
			eq_mult(aux1, ltfp[tk][tx], stf); 
			eq_mult(aux2, lti, stip[tk][tx]);

			for(int mu = 0; mu < 4; mu++)
			{
				aux1V[mu].glV(aux1, mu); 
				aux1A[mu].glA(aux1, mu); 
				aux2V[mu].glV(aux2, mu); 
				aux2A[mu].glA(aux2, mu); 
			}	

			for(int mu = 0; mu < 4; mu++) 
				for( int nu = 0; nu < 4; nu++)
				{
					dsdsVV = Trace(aux1V[mu], aux2V[nu]);
					dsdsVA = Trace(aux1V[mu], aux2A[nu]);
					dsdsAV = Trace(aux1A[mu], aux2V[nu]);
					dsdsAA = Trace(aux1A[mu], aux2A[nu]);

					_pc[4][index][tid] += 
						dsdsVV * uuVV[mu][nu] +
						dsdsVA * uuVA[mu][nu] +
						dsdsAV * uuAV[mu][nu] +
						dsdsAA * uuAA[mu][nu];
					_pv[4][index][tid] +=
						dsdsVV * uuAA[mu][nu] +
						dsdsVA * uuAV[mu][nu] +
						dsdsAV * uuVA[mu][nu] +
						dsdsAA * uuVV[mu][nu];
				}

			//contraction 6
			for(int mu = 0; mu < 4; mu++) 
				for(int nu =0; nu < 4; nu++)
				{
					dsu1VV = Trace(aux1V[mu], uV1[nu]);
					dsu1VA = Trace(aux1V[mu], uA1[nu]);
					dsu1AV = Trace(aux1A[mu], uV1[nu]);
					dsu1AA = Trace(aux1A[mu], uA1[nu]);

					dsu2VV = Trace(uV2[mu], aux2V[nu]);
					dsu2VA = Trace(uV2[mu], aux2A[nu]);
					dsu2AV = Trace(uA2[mu], aux2V[nu]);
					dsu2AA = Trace(uA2[mu], aux2A[nu]);

					_pc[5][index][tid] +=
						dsu1VV * dsu2VV + dsu1AA * dsu2AA +
						dsu1VA * dsu2VA + dsu1AV * dsu2AV;
					_pv[5][index][tid] +=
						dsu1VV * dsu2AA + dsu1VA * dsu2AV +
						dsu1AV * dsu2VA + dsu1AA * dsu2VV;
				}

			//contraction 7 and 8
			for(int mu = 0; mu < 4; mu++)
			{
				eq_mult(dsds1V[mu], aux1, aux2V[mu]);
				eq_mult(dsds1A[mu], aux1, aux2A[mu]);
				eq_mult(dsds2V[mu], aux2, aux1V[mu]);
				eq_mult(dsds2A[mu], aux2, aux1A[mu]);
			}

			for(int mu = 0; mu < 4; mu++) 
			{
				dsds1ConV[mu] = dsds1VioV[mu] = dsds1ConA[mu] = dsds1VioA[mu] = 0;
				dsds2ConV[mu] = dsds2VioV[mu] = dsds2ConA[mu] = dsds2VioA[mu] = 0;
			}
			for(int mu = 0; mu < 4; mu++) 
				for(int nu = 0; nu < 4; nu++)
				{
					tmp1.glV(dsds1V[mu], nu);
					tmp2.grV(tmp1, nu);
					dsds1ConV[mu] += tmp2;
					tmp2.grA(tmp1, nu);
					dsds1VioV[mu] += tmp2;
					tmp1.glA(dsds1V[mu], nu);
					tmp2.grA(tmp1, nu);
					dsds1ConV[mu] += tmp2;
					tmp2.grV(tmp1, nu);
					dsds1VioV[mu] += tmp2;

					tmp1.glV(dsds1A[mu], nu);
					tmp2.grV(tmp1, nu);
					dsds1ConA[mu] += tmp2;
					tmp2.grA(tmp1, nu);
					dsds1VioA[mu] += tmp2;
					tmp1.glA(dsds1A[mu], nu);
					tmp2.grA(tmp1, nu);
					dsds1ConA[mu] += tmp2;
					tmp2.grV(tmp1, nu);
					dsds1VioA[mu] += tmp2;

					tmp1.glV(dsds2V[mu], nu);
					tmp2.grV(tmp1, nu);
					dsds2ConV[mu] += tmp2;
					tmp2.grA(tmp1, nu);
					dsds2VioV[mu] += tmp2;
					tmp1.glA(dsds2V[mu], nu);
					tmp2.grA(tmp1, nu);
					dsds2ConV[mu] += tmp2;
					tmp2.grV(tmp1, nu);
					dsds2VioV[mu] += tmp2;

					tmp1.glV(dsds2A[mu], nu);
					tmp2.grV(tmp1, nu);
					dsds2ConA[mu] += tmp2;
					tmp2.grA(tmp1, nu);
					dsds2VioA[mu] += tmp2;
					tmp1.glA(dsds2A[mu], nu);
					tmp2.grA(tmp1, nu);
					dsds2ConA[mu] += tmp2;
					tmp2.grV(tmp1, nu);
					dsds2VioA[mu] += tmp2;
				}

			for(int mu = 0; mu < 4; mu++) 
			{
				_pc[7][index][tid] +=
					Trace(dsds1ConV[mu], uuV1[mu]) +
					Trace(dsds1ConA[mu], uuA1[mu]);
				_pv[7][index][tid] +=
					Trace(dsds1VioV[mu], uuA1[mu]) +
					Trace(dsds1VioA[mu], uuV1[mu]);
				_pc[6][index][tid] +=
					Trace(dsds2ConV[mu], uuV2[mu]) +
					Trace(dsds2ConA[mu], uuA2[mu]);
				_pv[6][index][tid] +=
					Trace(dsds2VioV[mu], uuA2[mu]) +
					Trace(dsds2VioA[mu], uuV2[mu]);
			}	
		}
	}

	//reduction
	Rcomplex pc[8*T*T], pv[8*T*T];
	for(int i = 0; i < 8; i++)
		for(int j = 0; j < T*T; j++)
		{
			int index = i*T*T + j;
			pc[index] = pv[index] = 0.0;
			for(int tid = 0; tid < N_threads; tid++)
			{
				pc[index] += _pc[i][j][tid];
				pv[index] += _pv[i][j][tid];
			}
		}

#ifdef PARALLEL
	QMP_sum_double_array((double*)pc, 2*8*T*T);
	QMP_sum_double_array((double*)pv, 2*8*T*T);
#endif

	//normalization
	for(int i=0; i<8*T*T; i++) 
	{
		pc[i] *= 0.5 * size_3d_glb / glb_t;
		pv[i] *= 0.5 * size_3d_glb / glb_t;
	}

	//output result
	char filename[512];
	//	double sum_0=0;

	sprintf(filename, "%s_type12_pc_sep%d_%s.txt", common_arg->filename, boxsize, label);
	SaveCorr(pc, filename, 8*T*T);
	sprintf(filename, "%s_type12_pv_sep%d_%s.txt", common_arg->filename, boxsize, label);
	SaveCorr(pv, filename, 8*T*T);

	//omp_set_num_threads(N_threads_old);

	cost += dclock();
	print_time(cname, fname, cost);
}	
