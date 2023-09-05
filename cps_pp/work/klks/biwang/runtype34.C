#include "algklks.h"
void AlgKlKs::RunType34(const char* label)
{
	const char *fname = "RunType34()";
//	static Timer time(fname);
//	time.start();
	Float cost = -dclock();

//	const int N_threads_old=bfmarg::default_threads;
	const int N_threads = 64;
	int hit_sub;//subtracting the contribution from all hits are expensive, use all mode averaging
	hit_sub=(nhits>6)?6:nhits;
	omp_set_num_threads(N_threads);

	Rcomplex pc[18*T*T];
	Rcomplex pv[18*T*T];
	memset(pc,0,sizeof(Rcomplex)*18*T*T);
	memset(pv,0,sizeof(Rcomplex)*18*T*T);
	//allocate memeory
	const int mem_size = 3*T;

	rand_uppc = new WilsonMatrix* [1+hit_sub];//last one store the average
	rand_uppv = new WilsonMatrix* [1+hit_sub];
	rand_downpc = new WilsonMatrix* [1+hit_sub];
	rand_downpv = new WilsonMatrix* [1+hit_sub];

	rand_leftpc = new Rcomplex* [1+hit_sub];
	rand_leftpv = new Rcomplex* [1+hit_sub];
	rand_rightpc = new Rcomplex* [1+hit_sub];
	rand_rightpv = new Rcomplex* [1+hit_sub];

	for(int i = 0; i < hit_sub+1; i++)
	{
		rand_uppc[i] = new WilsonMatrix [mem_size];
		rand_uppv[i] = new WilsonMatrix [mem_size];
		rand_downpc[i] = new WilsonMatrix [mem_size];
		rand_downpv[i] = new WilsonMatrix [mem_size];

		rand_leftpc[i] = new Rcomplex [mem_size];
		rand_leftpv[i] = new Rcomplex [mem_size];
		rand_rightpc[i] = new Rcomplex [mem_size];
		rand_rightpv[i] = new Rcomplex [mem_size];
	}


	for(int tkaon = 0; tkaon < glb_t; tkaon++)
	{
		for(int hit = 0; hit < hit_sub; hit++)
		{
			const int ta = (tkaon+minsep)%glb_t;
			const int tb = (tkaon+minsep+boxsize)%glb_t;
			const int tkb = (tkaon+2*minsep+boxsize)%glb_t; 

			//for reduction
			WilsonMatrix _uppc[3][T][N_threads], _downpc[3][T][N_threads];
			WilsonMatrix _uppv[3][T][N_threads], _downpv[3][T][N_threads];
			Rcomplex _leftpc[3][T][N_threads], _rightpc[3][T][N_threads];
			Rcomplex _leftpv[3][T][N_threads], _rightpv[3][T][N_threads];

			memset(_uppc, 0, sizeof(WilsonMatrix)*3*T*N_threads);
			memset(_uppv, 0, sizeof(WilsonMatrix)*3*T*N_threads);
			memset(_downpc, 0, sizeof(WilsonMatrix)*3*T*N_threads);
			memset(_downpv, 0, sizeof(WilsonMatrix)*3*T*N_threads);
			memset(_leftpc, 0, sizeof(Rcomplex)*3*T*N_threads);
			memset(_leftpv, 0, sizeof(Rcomplex)*3*T*N_threads);
			memset(_rightpc, 0, sizeof(Rcomplex)*3*T*N_threads);
			memset(_rightpv, 0, sizeof(Rcomplex)*3*T*N_threads);

#pragma omp parallel for schedule(dynamic)
			for(int site = 0; site < GJP.VolNodeSites(); site++)
			{
				Site s(site);
				if(tb>=ta && (s.physT()>tb || s.physT()<ta)) continue;
				if(tb<ta && s.physT()>tb && s.physT()<ta) continue;
				int physT = s.physT();
				int tid = omp_get_thread_num();
				int diff = (physT+glb_t-ta)%glb_t;
				WilsonMatrix dwp1, swp1, dwp2, swp2, umcpp;
				WilsonMatrix tmp1, tmp2, tmp3, tmp4;
				dwp1 = (*lwall[tkaon])[site];

				tmp1.hconj((*swall[tkaon])[site]);
				swp1.grV(tmp1, -5);

				dwp2 = (*lwall[tkb])[site];

				tmp1.hconj((*swall[tkb])[site]);
				swp2.grV(tmp1, -5);

				umcpp = (*lmcloop[hit])[site];

				for(int mu = 0; mu < 4; mu++)
				{
					tmp1.glV(dwp1, mu);
					eq_mult(tmp2, swp2, tmp1);
					tmp3.glV(umcpp, mu);
					_uppc[0][diff][tid] += tmp2 * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp1);
					eq_mult(tmp3, swp2, tmp4);
					_uppc[1][diff][tid] += tmp3;
					tmp3.glA(umcpp, mu);
					_uppv[0][diff][tid] += -1.0*tmp2 * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp1);
					eq_mult(tmp3, swp2, tmp4);
					_uppv[1][diff][tid] += -1.0*tmp3;

					tmp1.glA(dwp1, mu);
					eq_mult(tmp2, swp2, tmp1);
					tmp3.glA(umcpp, mu);
					_uppc[0][diff][tid] += tmp2 * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp1);
					eq_mult(tmp3, swp2, tmp4);
					_uppc[1][diff][tid] += tmp3;
					tmp3.glV(umcpp, mu);
					_uppv[0][diff][tid] += -1.0*tmp2 * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp1);
					eq_mult(tmp3, swp2, tmp4);
					_uppv[1][diff][tid] += -1.0*tmp3;

					tmp1.glV(dwp2, mu);
					eq_mult(tmp2, swp1, tmp1);
					tmp3.glV(umcpp, mu);
					_downpc[0][diff][tid] += tmp2 * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp1);
					eq_mult(tmp3, swp1, tmp4);
					_downpc[1][diff][tid] += tmp3;
					tmp3.glA(umcpp, mu);
					_downpv[0][diff][tid] += -1.0*tmp2 * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp1);
					eq_mult(tmp3, swp1, tmp4);
					_downpv[1][diff][tid] += -1.0*tmp3;

					tmp1.glA(dwp2, mu);
					eq_mult(tmp2, swp1, tmp1);
					tmp3.glA(umcpp, mu);
					_downpc[0][diff][tid] += tmp2 * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp1);
					eq_mult(tmp3, swp1, tmp4);
					_downpc[1][diff][tid] += tmp3;
					tmp3.glV(umcpp, mu);
					_downpv[0][diff][tid] += -1.0*tmp2 * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp1);
					eq_mult(tmp3, swp1, tmp4);
					_downpv[1][diff][tid] += -1.0*tmp3;

					tmp1.glV(dwp1, mu);
					eq_mult(tmp2, tmp1, swp1);
					tmp3.glV(umcpp, mu);
					_leftpc[0][diff][tid] += tmp2.Trace() * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp2);
					_leftpc[1][diff][tid] += tmp4.Trace();
					tmp3.glA(umcpp, mu);
					_leftpv[0][diff][tid] += -1.0*tmp2.Trace() * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp2);
					_leftpv[1][diff][tid] += -1.0*tmp4.Trace();

					tmp1.glA(dwp1, mu);
					eq_mult(tmp2, tmp1, swp1);
					tmp3.glA(umcpp, mu);
					_leftpc[0][diff][tid] += tmp2.Trace() * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp2);
					_leftpc[1][diff][tid] += tmp4.Trace();
					tmp3.glV(umcpp, mu);
					_leftpv[0][diff][tid] += -1.0*tmp2.Trace() * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp2);
					_leftpv[1][diff][tid] += -1.0*tmp4.Trace();

					tmp1.glV(dwp2, mu);
					eq_mult(tmp2, tmp1, swp2);
					tmp3.glV(umcpp, mu);
					_rightpc[0][diff][tid] += tmp2.Trace() * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp2);
					_rightpc[1][diff][tid] += tmp4.Trace();
					tmp3.glA(umcpp, mu);
					_rightpv[0][diff][tid] += -1.0*tmp2.Trace() * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp2);
					_rightpv[1][diff][tid] += -1.0*tmp4.Trace();

					tmp1.glA(dwp2, mu);
					eq_mult(tmp2, tmp1,swp2);
					tmp3.glA(umcpp, mu);
					_rightpc[0][diff][tid] += tmp2.Trace() * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp2);
					_rightpc[1][diff][tid] += tmp4.Trace();
					tmp3.glV(umcpp, mu);
					_rightpv[0][diff][tid] += -1.0*tmp2.Trace() * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp2);
					_rightpv[1][diff][tid] += -1.0*tmp4.Trace();
				}

				eq_mult(tmp1, swp2, dwp1);
				_uppc[2][diff][tid] += tmp1;
				tmp1.grV(swp2, -5);
				eq_mult(tmp2, tmp1, dwp1);
				_uppv[2][diff][tid] += tmp2;

				eq_mult(tmp1, swp1, dwp2);
				_downpc[2][diff][tid] += tmp1;
				tmp1.grV(swp1, -5);
				eq_mult(tmp2, tmp1, dwp2);
				_downpv[2][diff][tid] += tmp2;

				_leftpc[2][diff][tid] += Trace(swp1, dwp1);
				tmp1.grV(swp1, -5);
				_leftpv[2][diff][tid] += Trace(tmp1, dwp1);

				_rightpc[2][diff][tid] += Trace(swp2, dwp2);
				tmp1.grV(swp2, -5);
				_rightpv[2][diff][tid] += Trace(tmp1, dwp2);
			}

			for(int i = 0; i < 3; i++)
				for(int j = 0; j < T; j++)
				{
					int ind = j + i*T ;
					rand_uppc[hit][ind] = rand_downpc[hit][ind] = rand_uppv[hit][ind] = rand_downpv[hit][ind] = 0.0;
					rand_leftpc[hit][ind] = rand_leftpv[hit][ind] = rand_rightpc[hit][ind] = rand_rightpv[hit][ind] = 0.0;

					for(int tid = 0; tid < N_threads; tid++)
					{
						rand_uppc[hit][ind] += _uppc[i][j][tid];
						rand_uppv[hit][ind] += _uppv[i][j][tid];
						rand_downpc[hit][ind] += _downpc[i][j][tid];
						rand_downpv[hit][ind] += _downpv[i][j][tid];
						rand_leftpc[hit][ind] += _leftpc[i][j][tid];
						rand_leftpv[hit][ind] += _leftpv[i][j][tid];
						rand_rightpc[hit][ind] += _rightpc[i][j][tid];
						rand_rightpv[hit][ind] += _rightpv[i][j][tid];
					}
				}
		}//end hit loop

		memset(rand_uppc[hit_sub],0,sizeof(WilsonMatrix)*3*T);
		memset(rand_uppv[hit_sub],0,sizeof(WilsonMatrix)*3*T);
		memset(rand_downpc[hit_sub],0,sizeof(WilsonMatrix)*3*T);
		memset(rand_downpv[hit_sub],0,sizeof(WilsonMatrix)*3*T);
		memset(rand_leftpc[hit_sub],0,sizeof(Rcomplex)*3*T);
		memset(rand_leftpv[hit_sub],0,sizeof(Rcomplex)*3*T);
		memset(rand_rightpc[hit_sub],0,sizeof(Rcomplex)*3*T);
		memset(rand_rightpv[hit_sub],0,sizeof(Rcomplex)*3*T);
		{
			const int ta = (tkaon+minsep)%glb_t;
			const int tb = (tkaon+minsep+boxsize)%glb_t;
			const int tkb = (tkaon+2*minsep+boxsize)%glb_t; 

			//for reduction
			WilsonMatrix _uppc[3][T][N_threads], _downpc[3][T][N_threads];
			WilsonMatrix _uppv[3][T][N_threads], _downpv[3][T][N_threads];
			Rcomplex _leftpc[3][T][N_threads], _rightpc[3][T][N_threads];
			Rcomplex _leftpv[3][T][N_threads], _rightpv[3][T][N_threads];

			memset(_uppc, 0, sizeof(WilsonMatrix)*3*T*N_threads);
			memset(_uppv, 0, sizeof(WilsonMatrix)*3*T*N_threads);
			memset(_downpc, 0, sizeof(WilsonMatrix)*3*T*N_threads);
			memset(_downpv, 0, sizeof(WilsonMatrix)*3*T*N_threads);
			memset(_leftpc, 0, sizeof(Rcomplex)*3*T*N_threads);
			memset(_leftpv, 0, sizeof(Rcomplex)*3*T*N_threads);
			memset(_rightpc, 0, sizeof(Rcomplex)*3*T*N_threads);
			memset(_rightpv, 0, sizeof(Rcomplex)*3*T*N_threads);

#pragma omp parallel for schedule(dynamic)
			for(int site = 0; site < GJP.VolNodeSites(); site++)
			{
				Site s(site);
				if(tb>=ta && (s.physT()>tb || s.physT()<ta)) continue;
				if(tb<ta && s.physT()>tb && s.physT()<ta) continue;
				int physT = s.physT();
				int tid = omp_get_thread_num();
				int diff = (physT+glb_t-ta)%glb_t;
				WilsonMatrix dwp1, swp1, dwp2, swp2, umcpp;
				WilsonMatrix tmp1, tmp2, tmp3, tmp4;
				dwp1 = (*lwall[tkaon])[site];

				tmp1.hconj((*swall[tkaon])[site]);
				swp1.grV(tmp1, -5);

				dwp2 = (*lwall[tkb])[site];

				tmp1.hconj((*swall[tkb])[site]);
				swp2.grV(tmp1, -5);

				umcpp = (*lmcloop_avg)[site];

				for(int mu = 0; mu < 4; mu++)
				{
					tmp1.glV(dwp1, mu);
					eq_mult(tmp2, swp2, tmp1);
					tmp3.glV(umcpp, mu);
					_uppc[0][diff][tid] += tmp2 * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp1);
					eq_mult(tmp3, swp2, tmp4);
					_uppc[1][diff][tid] += tmp3;
					tmp3.glA(umcpp, mu);
					_uppv[0][diff][tid] += -1.0*tmp2 * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp1);
					eq_mult(tmp3, swp2, tmp4);
					_uppv[1][diff][tid] += -1.0*tmp3;

					tmp1.glA(dwp1, mu);
					eq_mult(tmp2, swp2, tmp1);
					tmp3.glA(umcpp, mu);
					_uppc[0][diff][tid] += tmp2 * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp1);
					eq_mult(tmp3, swp2, tmp4);
					_uppc[1][diff][tid] += tmp3;
					tmp3.glV(umcpp, mu);
					_uppv[0][diff][tid] += -1.0*tmp2 * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp1);
					eq_mult(tmp3, swp2, tmp4);
					_uppv[1][diff][tid] += -1.0*tmp3;

					tmp1.glV(dwp2, mu);
					eq_mult(tmp2, swp1, tmp1);
					tmp3.glV(umcpp, mu);
					_downpc[0][diff][tid] += tmp2 * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp1);
					eq_mult(tmp3, swp1, tmp4);
					_downpc[1][diff][tid] += tmp3;
					tmp3.glA(umcpp, mu);
					_downpv[0][diff][tid] += -1.0*tmp2 * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp1);
					eq_mult(tmp3, swp1, tmp4);
					_downpv[1][diff][tid] += -1.0*tmp3;

					tmp1.glA(dwp2, mu);
					eq_mult(tmp2, swp1, tmp1);
					tmp3.glA(umcpp, mu);
					_downpc[0][diff][tid] += tmp2 * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp1);
					eq_mult(tmp3, swp1, tmp4);
					_downpc[1][diff][tid] += tmp3;
					tmp3.glV(umcpp, mu);
					_downpv[0][diff][tid] += -1.0*tmp2 * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp1);
					eq_mult(tmp3, swp1, tmp4);
					_downpv[1][diff][tid] += -1.0*tmp3;

					tmp1.glV(dwp1, mu);
					eq_mult(tmp2, tmp1, swp1);
					tmp3.glV(umcpp, mu);
					_leftpc[0][diff][tid] += tmp2.Trace() * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp2);
					_leftpc[1][diff][tid] += tmp4.Trace();
					tmp3.glA(umcpp, mu);
					_leftpv[0][diff][tid] += -1.0*tmp2.Trace() * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp2);
					_leftpv[1][diff][tid] += -1.0*tmp4.Trace();

					tmp1.glA(dwp1, mu);
					eq_mult(tmp2, tmp1, swp1);
					tmp3.glA(umcpp, mu);
					_leftpc[0][diff][tid] += tmp2.Trace() * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp2);
					_leftpc[1][diff][tid] += tmp4.Trace();
					tmp3.glV(umcpp, mu);
					_leftpv[0][diff][tid] += -1.0*tmp2.Trace() * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp2);
					_leftpv[1][diff][tid] += -1.0*tmp4.Trace();

					tmp1.glV(dwp2, mu);
					eq_mult(tmp2, tmp1, swp2);
					tmp3.glV(umcpp, mu);
					_rightpc[0][diff][tid] += tmp2.Trace() * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp2);
					_rightpc[1][diff][tid] += tmp4.Trace();
					tmp3.glA(umcpp, mu);
					_rightpv[0][diff][tid] += -1.0*tmp2.Trace() * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp2);
					_rightpv[1][diff][tid] += -1.0*tmp4.Trace();

					tmp1.glA(dwp2, mu);
					eq_mult(tmp2, tmp1,swp2);
					tmp3.glA(umcpp, mu);
					_rightpc[0][diff][tid] += tmp2.Trace() * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp2);
					_rightpc[1][diff][tid] += tmp4.Trace();
					tmp3.glV(umcpp, mu);
					_rightpv[0][diff][tid] += -1.0*tmp2.Trace() * tmp3.Trace();
					eq_mult(tmp4, tmp3, tmp2);
					_rightpv[1][diff][tid] += -1.0*tmp4.Trace();
				}

				eq_mult(tmp1, swp2, dwp1);
				_uppc[2][diff][tid] += tmp1;
				tmp1.grV(swp2, -5);
				eq_mult(tmp2, tmp1, dwp1);
				_uppv[2][diff][tid] += tmp2;

				eq_mult(tmp1, swp1, dwp2);
				_downpc[2][diff][tid] += tmp1;
				tmp1.grV(swp1, -5);
				eq_mult(tmp2, tmp1, dwp2);
				_downpv[2][diff][tid] += tmp2;

				_leftpc[2][diff][tid] += Trace(swp1, dwp1);
				tmp1.grV(swp1, -5);
				_leftpv[2][diff][tid] += Trace(tmp1, dwp1);

				_rightpc[2][diff][tid] += Trace(swp2, dwp2);
				tmp1.grV(swp2, -5);
				_rightpv[2][diff][tid] += Trace(tmp1, dwp2);
			}

			for(int i = 0; i < 3; i++){
				for(int j = 0; j < T; j++)
				{
					int ind = j + i*T ;
					for(int tid = 0; tid < N_threads; tid++)
					{
						rand_uppc[hit_sub][ind] += _uppc[i][j][tid];
						rand_uppv[hit_sub][ind] += _uppv[i][j][tid];
						rand_downpc[hit_sub][ind] += _downpc[i][j][tid];
						rand_downpv[hit_sub][ind] += _downpv[i][j][tid];
						rand_leftpc[hit_sub][ind] += _leftpc[i][j][tid];
						rand_leftpv[hit_sub][ind] += _leftpv[i][j][tid];
						rand_rightpc[hit_sub][ind] += _rightpc[i][j][tid];
						rand_rightpv[hit_sub][ind] += _rightpv[i][j][tid];
					}
				}
			}
		}



#ifdef PARALLEL
		for(int i = 0; i < hit_sub+1; i++)
		{
			QMP_sum_double_array((double*)(rand_uppc[i]), 288*3*T);
			QMP_sum_double_array((double*)(rand_uppv[i]), 288*3*T);
			QMP_sum_double_array((double*)(rand_downpc[i]), 288*3*T);
			QMP_sum_double_array((double*)(rand_downpv[i]), 288*3*T);
			QMP_sum_double_array((double*)(rand_leftpc[i]), 2*3*T);
			QMP_sum_double_array((double*)(rand_leftpv[i]), 2*3*T);
			QMP_sum_double_array((double*)(rand_rightpc[i]), 2*3*T);
			QMP_sum_double_array((double*)(rand_rightpv[i]), 2*3*T);
		}
#endif

		const int ta = (tkaon+minsep)%glb_t;
		const int tb = (tkaon+minsep+boxsize)%glb_t;
		const int tkb = (tkaon+2*minsep+boxsize)%glb_t; 

		//make a copy, for my convinience
		WilsonMatrix uppc[hit_sub+1][3][T], downpc[hit_sub+1][3][T];
		WilsonMatrix uppv[hit_sub+1][3][T], downpv[hit_sub+1][3][T];
		Rcomplex leftpc[hit_sub+1][3][T], rightpc[hit_sub+1][3][T];
		Rcomplex leftpv[hit_sub+1][3][T], rightpv[hit_sub+1][3][T];
		for(int hit = 0; hit < hit_sub+1; hit++)
			for(int i = 0; i < 3; i++)
				for(int j = 0; j < T; j++)
				{
					int ind = j + i*T ;
					uppc[hit][i][j] = rand_uppc[hit][ind];
					uppv[hit][i][j] = rand_uppv[hit][ind];
					downpc[hit][i][j] = rand_downpc[hit][ind];
					downpv[hit][i][j] = rand_downpv[hit][ind];
					leftpc[hit][i][j] = rand_leftpc[hit][ind];
					leftpv[hit][i][j] = rand_leftpv[hit][ind];
					rightpc[hit][i][j] = rand_rightpc[hit][ind];
					rightpv[hit][i][j] = rand_rightpv[hit][ind];
				}
		// save leftpc/pv, rightpc/pv to file
		char filename0[500];
		sprintf(filename0, "%s_leftpc_tk%d_%s.txt", common_arg->filename, tkaon, label);
		SaveCorr(leftpc[0][0], filename0, (hit_sub+1)*3*T);
		sprintf(filename0, "%s_leftpv_tk%d_%s.txt", common_arg->filename, tkaon, label);
		SaveCorr(leftpv[0][0], filename0, (hit_sub+1)*3*T);
		sprintf(filename0, "%s_rightpc_tkb%d_%s.txt", common_arg->filename, tkb, label);
		SaveCorr(rightpc[0][0], filename0, (hit_sub+1)*3*T);
		sprintf(filename0, "%s_rightpv_tkb%d_%s.txt", common_arg->filename, tkb, label);
		SaveCorr(rightpv[0][0], filename0, (hit_sub+1)*3*T);

#pragma omp parallel for
		for(int ind = 0; ind < T*T; ind++){
			int tx = ind / T;
			int ty = ind % T; 
			//double norm = (tx==ty)?(1.0*(nhits*nhits)/(nhits*(nhits-1))):1.0;
			double norm = (1.0*nhits)/(nhits-1);
			double factor=1.0/(-1+nhits)/hit_sub;
			pc[0*T*T+ind] += Trace(uppc[hit_sub][0][tx], downpc[hit_sub][0][ty]) *norm;
			pc[1*T*T+ind] += Trace(uppc[hit_sub][1][tx], downpc[hit_sub][1][ty]) *norm;
			pc[2*T*T+ind] += Trace(uppc[hit_sub][1][tx], downpc[hit_sub][0][ty]) *norm;
			pc[3*T*T+ind] += Trace(uppc[hit_sub][0][tx], downpc[hit_sub][1][ty]) *norm;
			pc[4*T*T+ind] += leftpc[hit_sub][0][tx] * rightpc[hit_sub][0][ty] *norm;
			pc[5*T*T+ind] += leftpc[hit_sub][1][tx] * rightpc[hit_sub][1][ty] *norm;
			pc[6*T*T+ind] += leftpc[hit_sub][1][tx] * rightpc[hit_sub][0][ty] *norm;
			pc[7*T*T+ind] += leftpc[hit_sub][0][tx] * rightpc[hit_sub][1][ty] *norm;
			pc[8*T*T+ind] += Trace(uppc[hit_sub][2][tx], downpc[hit_sub][0][ty]) ;
			pc[9*T*T+ind] += Trace(uppc[hit_sub][2][tx], downpc[hit_sub][1][ty]) ;
			pc[10*T*T+ind] += Trace(uppc[hit_sub][0][tx], downpc[hit_sub][2][ty]);
			pc[11*T*T+ind] += Trace(uppc[hit_sub][1][tx], downpc[hit_sub][2][ty]);
			pc[12*T*T+ind] += Trace(uppc[hit_sub][2][tx], downpc[hit_sub][2][ty]);
			pc[13*T*T+ind] += leftpc[hit_sub][2][tx] * rightpc[hit_sub][0][ty] ;
			pc[14*T*T+ind] += leftpc[hit_sub][2][tx] * rightpc[hit_sub][1][ty] ;
			pc[15*T*T+ind] += leftpc[hit_sub][0][tx] * rightpc[hit_sub][2][ty] ;
			pc[16*T*T+ind] += leftpc[hit_sub][1][tx] * rightpc[hit_sub][2][ty] ;
			pc[17*T*T+ind] += leftpc[hit_sub][2][tx] * rightpc[hit_sub][2][ty] ;

			pv[0*T*T+ind] += Trace(uppv[hit_sub][0][tx], downpv[hit_sub][0][ty]) *norm;
			pv[1*T*T+ind] += Trace(uppv[hit_sub][1][tx], downpv[hit_sub][1][ty]) *norm;
			pv[2*T*T+ind] += Trace(uppv[hit_sub][1][tx], downpv[hit_sub][0][ty]) *norm;
			pv[3*T*T+ind] += Trace(uppv[hit_sub][0][tx], downpv[hit_sub][1][ty])*norm;
			pv[4*T*T+ind] += leftpv[hit_sub][0][tx] * rightpv[hit_sub][0][ty] *norm;
			pv[5*T*T+ind] += leftpv[hit_sub][1][tx] * rightpv[hit_sub][1][ty] *norm;
			pv[6*T*T+ind] += leftpv[hit_sub][1][tx] * rightpv[hit_sub][0][ty]*norm;
			pv[7*T*T+ind] += leftpv[hit_sub][0][tx] * rightpv[hit_sub][1][ty] *norm;
			pv[8*T*T+ind] += Trace(uppv[hit_sub][2][tx], downpv[hit_sub][0][ty]);
			pv[9*T*T+ind] += Trace(uppv[hit_sub][2][tx], downpv[hit_sub][1][ty]);
			pv[10*T*T+ind] += Trace(uppv[hit_sub][0][tx], downpv[hit_sub][2][ty]);
			pv[11*T*T+ind] += Trace(uppv[hit_sub][1][tx], downpv[hit_sub][2][ty]);
			pv[12*T*T+ind] += Trace(uppv[hit_sub][2][tx], downpv[hit_sub][2][ty]);
			pv[13*T*T+ind] += leftpv[hit_sub][2][tx] * rightpv[hit_sub][0][ty];
			pv[14*T*T+ind] += leftpv[hit_sub][2][tx] * rightpv[hit_sub][1][ty];
			pv[15*T*T+ind] += leftpv[hit_sub][0][tx] * rightpv[hit_sub][2][ty] ;
			pv[16*T*T+ind] += leftpv[hit_sub][1][tx] * rightpv[hit_sub][2][ty] ;
			pv[17*T*T+ind] += leftpv[hit_sub][2][tx] * rightpv[hit_sub][2][ty] ;
			for(int hit=0;hit<hit_sub;hit++){
				pc[0*T*T+ind] -= Trace(uppc[hit][0][tx], downpc[hit][0][ty])*factor;
				pc[1*T*T+ind] -= Trace(uppc[hit][1][tx], downpc[hit][1][ty])*factor ;
				pc[2*T*T+ind] -= Trace(uppc[hit][1][tx], downpc[hit][0][ty]) *factor;
				pc[3*T*T+ind] -= Trace(uppc[hit][0][tx], downpc[hit][1][ty]) *factor;
				pc[4*T*T+ind] -= leftpc[hit][0][tx] * rightpc[hit][0][ty] *factor;
				pc[5*T*T+ind] -= leftpc[hit][1][tx] * rightpc[hit][1][ty] *factor;
				pc[6*T*T+ind] -= leftpc[hit][1][tx] * rightpc[hit][0][ty] *factor;
				pc[7*T*T+ind] -= leftpc[hit][0][tx] * rightpc[hit][1][ty] *factor;

				pv[0*T*T+ind] -= Trace(uppv[hit][0][tx], downpv[hit][0][ty]) *factor;
				pv[1*T*T+ind] -= Trace(uppv[hit][1][tx], downpv[hit][1][ty]) *factor;
				pv[2*T*T+ind] -= Trace(uppv[hit][1][tx], downpv[hit][0][ty]) *factor;
				pv[3*T*T+ind] -= Trace(uppv[hit][0][tx], downpv[hit][1][ty])*factor;
				pv[4*T*T+ind] -= leftpv[hit][0][tx] * rightpv[hit][0][ty] *factor;
				pv[5*T*T+ind] -= leftpv[hit][1][tx] * rightpv[hit][1][ty] *factor;
				pv[6*T*T+ind] -= leftpv[hit][1][tx] * rightpv[hit][0][ty]*factor;
				pv[7*T*T+ind] -= leftpv[hit][0][tx] * rightpv[hit][1][ty] *factor;
			}//end hit
		}//end ind
	}//end tk

	//normalization
	for(int i=0; i < 18*T*T; i++) 
	{
		pc[i] *= 1.0 / glb_t;
		pv[i] *= 1.0 / glb_t;
	}

	char filename[500];
	sprintf(filename, "%s_type34_pc_sep%d_%s.txt", common_arg->filename, boxsize, label);
	SaveCorr(pc, filename, 18*T*T);

	sprintf(filename, "%s_type34_pv_sep%d_%s.txt", common_arg->filename, boxsize, label);
	SaveCorr(pv, filename, 18*T*T);

	cost += dclock();
	print_time(cname, fname, cost);
//	omp_set_num_threads(N_threads_old);
	for(int i = 0; i < hit_sub+1; i++)
	{
		delete [] rand_uppc[i];
		delete [] rand_uppv[i];
		delete [] rand_downpc[i];
		delete [] rand_downpv[i];
		delete [] rand_leftpc[i];
		delete [] rand_leftpv[i];
		delete [] rand_rightpc[i];
		delete [] rand_rightpv[i];
	}
	delete [] rand_uppc;
	delete [] rand_uppv;
	delete [] rand_downpc;
	delete [] rand_downpv;
	delete [] rand_leftpc;
	delete [] rand_leftpv;
	delete [] rand_rightpc;
	delete [] rand_rightpv;
//	time.stop();
	syncNode();
}


