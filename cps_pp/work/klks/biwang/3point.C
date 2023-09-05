#include "algklks.h"
void AlgKlKs::RunThreePointCorr(const char* label) const
{
	const char * fname = "RunThreePointCorr()";
	//	static Timer time(fname);
	//	time.start();
	double cost = -dclock();

//	const int N_threads_old = bfmarg::default_threads;
	int N_threads = 64;
	omp_set_num_threads(64);

	int N_sep = (max_3pt_sep - min_3pt_sep)/sep_step + 1;

	Rcomplex diagram[N_sep][19][glb_t];
	memset(diagram, 0, 19*N_sep*glb_t*sizeof(Rcomplex));

	//time translation

	//	for(int delta = 0; delta < glb_t; delta++)
	for (int delta = min_3pt_sep; delta <=max_3pt_sep; delta += sep_step){
		int delta_ind = (delta - min_3pt_sep)/sep_step;
		for(int tk = 0; tk < glb_t; tk++)
		{
			int tpi = (tk+delta)%glb_t;

			WilsonMatrix dww, sww, tmp;
			dww = lwwprop[tk][tpi];
			sww.glV(swwprop[tpi][tk], -5);

			Rcomplex uwwtr, swwtr;
			tmp.glV(lwwprop[tpi][tpi], -5);
			uwwtr = tmp.Trace();
			tmp.glV(swwprop[tpi][tpi], -5);
			swwtr = tmp.Trace();

			//define this array to do reduction
			Rcomplex _diagram[19][glb_t][N_threads];
			memset(_diagram, 0, 19*glb_t*N_threads*sizeof(Rcomplex));
#pragma omp parallel for
			for(int site = 0; site < GJP.VolNodeSites(); site++)
			{
				Site s(site);
				int physT= s.physT();
				int diff = (physT + glb_t - tk)%glb_t; 
				int tid = omp_get_thread_num();

				WilsonMatrix tmp1, tmp2, tmp3;

				WilsonMatrix lwp1, lwp2, swp1, swp2, umcpp(0.0);
				lwp1 = (*lwall[tk])[site];
				lwp2 = (*lwall[tpi])[site];
				tmp1.hconj((*swall[tk])[site]);
				swp1.grV(tmp1, -5);
				tmp1.hconj((*swall[tpi])[site]);
				swp2.grV(tmp1, -5);

				umcpp = (*lmcloop_avg)[site];

				WilsonMatrix aux1, aux2, aux3, aux4, aux5;

				eq_mult(aux1, lwp1, swp1);
				tmp1.hconj(lwp2);
				tmp2.grV(tmp1, -5);
				eq_mult(aux2, lwp2, tmp2);

				eq_mult(tmp1, dww, swp1); 
				tmp2.glV(tmp1, -5);
				eq_mult(aux3, lwp2, tmp2);

				eq_mult(tmp1, lwp1, sww);
				eq_mult(aux4, tmp1, swp2);

				eq_mult(aux5, lwp2, swp2);

				Rcomplex tr1, tr2;
				for(int mu = 0; mu < 4; mu++)
				{
					//diagram 0 and 1
					tmp1.glV(aux1, mu);
					tmp2.glV(aux2, mu);
					_diagram[0][diff][tid] += tmp1.Trace() * tmp2.Trace();
					_diagram[1][diff][tid] += Trace(tmp1, tmp2);	
					//diagram 17 and 18
					tmp3.glV(aux5, mu);
					_diagram[17][diff][tid] += tmp1.Trace() * tmp3.Trace();
					_diagram[18][diff][tid] += Trace(tmp1, tmp3);

					//diagram 0 and 1
					tmp1.glA(aux1, mu);
					tmp2.glA(aux2, mu);
					_diagram[0][diff][tid] += tmp1.Trace() * tmp2.Trace();
					_diagram[1][diff][tid] += Trace(tmp1, tmp2);
					//diagram 17 and 18
					tmp3.glA(aux5, mu);
					_diagram[17][diff][tid] += tmp1.Trace() * tmp3.Trace();
					_diagram[18][diff][tid] += Trace(tmp1, tmp3);

					//diagram 2, 3, 4 ,5
					tmp1.glV(aux3, mu);
					tmp2.glV(umcpp, mu);
					tmp3.glV(aux4, mu);
					_diagram[2][diff][tid] += tmp1.Trace() * tmp2.Trace();
					_diagram[3][diff][tid] += Trace(tmp1, tmp2);
					_diagram[4][diff][tid] += tmp3.Trace() * tmp2.Trace();
					_diagram[5][diff][tid] += Trace(tmp3, tmp2);

					tmp1.glA(aux3, mu);
					tmp2.glA(umcpp, mu);
					tmp3.glA(aux4, mu);
					_diagram[2][diff][tid] += tmp1.Trace() * tmp2.Trace();
					_diagram[3][diff][tid] += Trace(tmp1, tmp2);
					_diagram[4][diff][tid] += tmp3.Trace() * tmp2.Trace();
					_diagram[5][diff][tid] += Trace(tmp3, tmp2);


					//diagram 6-11
					tmp1.glV(aux1, mu);
					tmp2.glV(umcpp, mu);
					tmp3.glA(umcpp, mu);
					tr1 = tmp1.Trace() * tmp2.Trace();
					tr2 = Trace(tmp1, tmp2);
					_diagram[6][diff][tid] += tr1 * uwwtr;
					_diagram[7][diff][tid] += tr2 * uwwtr;
					_diagram[8][diff][tid] += tr1 * swwtr;
					_diagram[9][diff][tid] += tr2 * swwtr;
					_diagram[10][diff][tid] += -1.0 * tmp1.Trace() * tmp3.Trace();
					_diagram[11][diff][tid] += -1.0 * Trace(tmp1, tmp3);

					tmp1.glA(aux1, mu);
					tmp2.glA(umcpp, mu);
					tmp3.glV(umcpp, mu);
					tr1 = tmp1.Trace() * tmp2.Trace();
					tr2 = Trace(tmp1, tmp2);
					_diagram[6][diff][tid] += tr1 * uwwtr;
					_diagram[7][diff][tid] += tr2 * uwwtr;
					_diagram[8][diff][tid] += tr1 * swwtr;
					_diagram[9][diff][tid] += tr2 * swwtr;
					_diagram[10][diff][tid] += -1.0 * tmp1.Trace() * tmp3.Trace();
					_diagram[11][diff][tid] += -1.0 * Trace(tmp1, tmp3);

				}

				//diagram 12-16
				_diagram[12][diff][tid] += aux3.Trace();
				_diagram[13][diff][tid] += aux4.Trace();
				_diagram[14][diff][tid] += aux1.Trace() * uwwtr;
				_diagram[15][diff][tid] += aux1.Trace() * swwtr;
				tmp1.grV(aux1, -5);
				_diagram[16][diff][tid] += tmp1.Trace();
			}

			//reduction
			for(int d = 0; d < 19; d++)
				for(int diff = 0; diff < glb_t; diff++)
					for(int tid = 0; tid < N_threads; tid++)
						diagram[delta_ind][d][diff] += _diagram[d][diff][tid] *(1.0/ glb_t);
		}
	}

#ifdef PARALLEL
	QMP_sum_double_array((double*)diagram, 2*19*N_sep*glb_t);
#endif

	char filename[512];
	sprintf(filename, "%s_ktox_%s.txt", common_arg->filename, label);
	SaveCorr((Rcomplex*)diagram, filename, 19*N_sep*glb_t);

//	omp_set_num_threads(N_threads_old);

	cost += dclock();
	print_time(cname, fname, cost);
	syncNode();
}

