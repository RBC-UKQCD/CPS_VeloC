#include"algklks.h"
#define Rcomplex std::complex<double>
///////we do not include the minus sign coming from the odd number of loops, 
//the minus sign appears in the contractions come from the V-A structure.
void AlgKlKs::Runk2pipi(const char* label) const
{
	const char * fname = "Runk2pipi()";
	//	static Timer time(fname);
	//	time.start();
	double cost = -dclock();
	int num_diag=12;

//	const int N_threads_old = bfmarg::default_threads;
	int N_threads = 64;
	omp_set_num_threads(64);
	int N_sep = (max_3pt_sep - min_3pt_sep)/sep_step + 1;

	Rcomplex diagram[N_sep][num_diag][glb_t];//diagram 1 5 9 17 33 mix 3 min4
	memset(diagram, 0, num_diag*N_sep*glb_t*sizeof(Rcomplex));

	//time translation

	//for(int delta = 0; delta < glb_t; delta++)//delta time between tk and tpipi
	for (int delta = min_3pt_sep; delta <=max_3pt_sep; delta += sep_step) {
		int delta_ind = (delta - min_3pt_sep) / sep_step;
		for(int tk = 0; tk < glb_t; tk++)
		{
			int tpi = (tk+delta)%glb_t;
			int tpi2=(tpi+sep)%glb_t;

			//define this array to do reduction
			Rcomplex _diagram[num_diag][glb_t][N_threads];
			memset(_diagram, 0, num_diag*glb_t*N_threads*sizeof(Rcomplex));
			Rcomplex disVterm;
			disVterm=Trace(lwwprop[tpi][tpi2].conj_cp(),lwwprop[tpi][tpi2])+Trace(lwwprop[tpi2][tpi].conj_cp(),lwwprop[tpi2][tpi]);
#pragma omp parallel for
			for(int site = 0; site < GJP.VolNodeSites(); site++)
			{
				WilsonMatrix Mtemp1, Mtemp2, Mtemp3, Mtemp4, Mtemp5, Wpart1_aux, Wpart2_aux, Wpart2s_aux,Wpart1, Wpart2, Wpart2s;
				Rcomplex Rpart1, Rpart2, Rpart2s, Raux;
				cps::SpinMatrix Spart1, Spart2, Spart2s;
				cps::Matrix Mpart1, Mpart2, Mpart2s;
				Site s(site);
				int physT= s.physT();
				int diff = (physT + glb_t - tk)%glb_t; 
				int tid = omp_get_thread_num();
				////////////////////////////////////////////////type 1
				Mtemp1=(*lwall[tpi])[s.Index()];
				Mtemp2=(*lwall[tpi])[s.Index()];
				Mtemp2.hconj();
				Mtemp3=(*lwall[tpi2])[s.Index()];
				Mtemp4=lwwprop[tk][tpi2];
				Mtemp4.gl(-5);
				Mtemp5=(*swall[tk])[s.Index()];
				Mtemp5.hconj();

				eq_mult(Wpart1_aux, Mtemp4, Mtemp5);//use Wpart1_aux as temperary storage
				eq_mult(Wpart2_aux, Mtemp3, Wpart1_aux);
				eq_mult(Wpart1_aux, Mtemp1, Mtemp2);


				for(int i=0;i<4;i++)
				{
					//VA term
					Wpart1.glV(Wpart1_aux,i);
					Wpart2.glA(Wpart2_aux,i);
					Rpart1=Wpart1.Trace();
					Rpart2=Wpart2.Trace();
					Raux=Rpart1*Rpart2;
					_diagram[0][diff][tid]-=Raux;//diagram 1
					Mpart1=SpinTrace(Wpart1);
					Mpart2=SpinTrace(Wpart2);
					Raux=Tr(Mpart1,Mpart2);
					_diagram[1][diff][tid]-=Raux;//diamgram 2
					Raux=Trace(Wpart1, Wpart2);
					_diagram[2][diff][tid]-=Raux;//diagram3
					Spart1=ColorTrace(Wpart1);
					Spart2=ColorTrace(Wpart2);
					Raux=Tr(Spart1,Spart2);
					_diagram[3][diff][tid]-=Raux;//diagram 4
					//AV term
					Wpart1.glA(Wpart1_aux,i);
					Wpart2.glV(Wpart2_aux,i);
					Rpart1=Wpart1.Trace();
					Rpart2=Wpart2.Trace();
					Raux=Rpart1*Rpart2;
					_diagram[0][diff][tid]-=Raux;//diagram 1
					Mpart1=SpinTrace(Wpart1);
					Mpart2=SpinTrace(Wpart2);
					Raux=Tr(Mpart1,Mpart2);
					_diagram[1][diff][tid]-=Raux;//diagram 2
					Raux=Trace(Wpart1, Wpart2);
					_diagram[2][diff][tid]-=Raux;//diagram 3
					Spart1=ColorTrace(Wpart1);
					Spart2=ColorTrace(Wpart2);
					Raux=Tr(Spart1,Spart2);
					_diagram[3][diff][tid]-=Raux;//diagram 4
				}

				//the second term
				Mtemp1=(*lwall[tpi2])[s.Index()];
				Mtemp2=(*lwall[tpi2])[s.Index()];
				Mtemp2.hconj();
				Mtemp3=(*lwall[tpi])[s.Index()];
				Mtemp4=lwwprop[tk][tpi];
				Mtemp4.gl(-5);
				Mtemp5=(*swall[tk])[s.Index()];
				Mtemp5.hconj();

				eq_mult(Wpart1_aux, Mtemp4, Mtemp5);//use Wpart1_aux as temperary storage
				eq_mult(Wpart2_aux, Mtemp3, Wpart1_aux);
				eq_mult(Wpart1_aux, Mtemp1, Mtemp2);

				for(int i=0;i<4;i++)
				{
					//VA term
					Wpart1.glV(Wpart1_aux,i);
					Wpart2.glA(Wpart2_aux,i);
					Rpart1=Wpart1.Trace();
					Rpart2=Wpart2.Trace();
					Raux=Rpart1*Rpart2;
					_diagram[0][diff][tid]-=Raux;//diagram 1 
					Mpart1=SpinTrace(Wpart1);
					Mpart2=SpinTrace(Wpart2);
					Raux=Tr(Mpart1,Mpart2);
					_diagram[1][diff][tid]-=Raux;
					Raux=Trace(Wpart1, Wpart2);
					_diagram[2][diff][tid]-=Raux;
					Spart1=ColorTrace(Wpart1);
					Spart2=ColorTrace(Wpart2);
					Raux=Tr(Spart1,Spart2);
					_diagram[3][diff][tid]-=Raux;
					//AV term
					Wpart1.glA(Wpart1_aux,i);
					Wpart2.glV(Wpart2_aux,i);
					Rpart1=Wpart1.Trace();
					Rpart2=Wpart2.Trace();
					Raux=Rpart1*Rpart2;
					_diagram[0][diff][tid]-=Raux;
					Mpart1=SpinTrace(Wpart1);
					Mpart2=SpinTrace(Wpart2);
					Raux=Tr(Mpart1,Mpart2);
					_diagram[1][diff][tid]-=Raux;
					Raux=Trace(Wpart1, Wpart2);
					_diagram[2][diff][tid]-=Raux;
					Spart1=ColorTrace(Wpart1);
					Spart2=ColorTrace(Wpart2);
					Raux=Tr(Spart1,Spart2);
					_diagram[3][diff][tid]-=Raux;
				}
				///////////////////////////////////////////////////////////Type 2
				Mtemp1=(*lwall[tpi])[s.Index()];
				Mtemp2=lwwprop[tpi2][tpi];
				Mtemp2.gl(-5);
				Mtemp3=(*lwall[tpi2])[s.Index()];
				Mtemp3.hconj();
				eq_mult(Wpart1_aux, Mtemp2, Mtemp3);//temporary storage
				eq_mult(Wpart2_aux, Mtemp1, Wpart1_aux);

				Mtemp1=(*lwall[tpi2])[s.Index()];
				Mtemp2=lwwprop[tpi][tpi2];
				Mtemp2.gl(-5);
				Mtemp3=(*lwall[tpi])[s.Index()];
				Mtemp3.hconj();
				eq_mult(Wpart1_aux, Mtemp2, Mtemp3);//temporary storage
				eq_mult(Wpart2s_aux, Mtemp1, Wpart1_aux);//used wpart2s as a temperary matrix
				Wpart2_aux+=Wpart2s_aux;

				Mtemp4=(*lwall[tk])[s.Index()];
				Mtemp5=(*swall[tk])[s.Index()];
				Mtemp5.hconj();
				eq_mult(Wpart1_aux, Mtemp4, Mtemp5);

				for(int i=0;i<4;i++)
				{
					//VA term
					Wpart1.glV(Wpart1_aux,i);
					Wpart2.glA(Wpart2_aux,i);
					Rpart1=Wpart1.Trace();
					Rpart2=Wpart2.Trace();
					Raux=Rpart1*Rpart2;
					_diagram[4][diff][tid]-=Raux;//diagram 5
					Mpart1=SpinTrace(Wpart1);
					Mpart2=SpinTrace(Wpart2);
					Raux=Tr(Mpart1,Mpart2);
					_diagram[5][diff][tid]-=Raux;//diagram 6
					//AV term
					Wpart1.glA(Wpart1_aux,i);
					Wpart2.glV(Wpart2_aux,i);
					Rpart1=Wpart1.Trace();
					Rpart2=Wpart2.Trace();
					Raux=Rpart1*Rpart2;
					_diagram[4][diff][tid]-=Raux;
					Mpart1=SpinTrace(Wpart1);
					Mpart2=SpinTrace(Wpart2);
					Raux=Tr(Mpart1,Mpart2);
					_diagram[5][diff][tid]-=Raux;
				}
				////////////////////////////////////////////////////////type 3
				Mtemp1=(*lwall[tpi])[s.Index()];
				Mtemp2=lwwprop[tpi][tpi2];
				Mtemp2.hconj();
				Mtemp3=lwwprop[tk][tpi2];
				Mtemp4=(*swall[tk])[s.Index()];
				Mtemp4.hconj();
				eq_mult(Wpart2_aux, Mtemp1, Mtemp2);//temporary
				eq_mult(Wpart2s_aux, Mtemp3, Mtemp4);//temporary
				eq_mult(Wpart1_aux, Wpart2_aux, Wpart2s_aux);
				Mtemp1=(*lwall[tpi2])[s.Index()];
				Mtemp2=lwwprop[tpi2][tpi];
				Mtemp2.hconj();
				Mtemp3=lwwprop[tk][tpi];
				eq_mult(Wpart2_aux, Mtemp1, Mtemp2);//temporary
				eq_mult(Wpart2s_aux, Mtemp3, Mtemp4);//temporary
				eq_mult(Mtemp1, Wpart2_aux, Wpart2s_aux);
				Wpart1_aux+=Mtemp1;

				Wpart2_aux=(*lmcloop_avg)[s.Index()];
				for(int i=0;i<4;i++)
				{
					//AA term(from VA term, gamma(5)*V -> A)
					Wpart1.glA(Wpart1_aux,i);
					Wpart2.glA(Wpart2_aux,i);
					Wpart2s.glA(Wpart2s_aux,i);
					Rpart1=Wpart1.Trace();
					Rpart2=Wpart2.Trace();
					Raux=Rpart1*Rpart2;
					_diagram[6][diff][tid]+=Raux; //diagram 7
					Mpart1=SpinTrace(Wpart1);
					Mpart2=SpinTrace(Wpart2);
					Raux=Tr(Mpart1,Mpart2);
					_diagram[7][diff][tid]+=Raux;//diagram 8
					//VV term (from AV term gamma5*A->V)
					Wpart1.glV(Wpart1_aux,i);
					Wpart2.glV(Wpart2_aux,i);
					Wpart2s.glV(Wpart2s_aux,i);
					Rpart1=Wpart1.Trace();
					Rpart2=Wpart2.Trace();
					Raux=Rpart1*Rpart2;
					_diagram[6][diff][tid]+=Raux;
					Mpart1=SpinTrace(Wpart1);
					Mpart2=SpinTrace(Wpart2);
					Raux=Tr(Mpart1,Mpart2);
					_diagram[7][diff][tid]+=Raux;
				}
				///////////////////////////////////////////////////////////////type 4
				Mtemp2=(*swall[tk])[s.Index()];
				Mtemp2.hconj();
				Mtemp1=(*lwall[tk])[s.Index()];

				eq_mult(Wpart1_aux, Mtemp1, Mtemp2);

				//light quark bubble
				Wpart2_aux=(*lmcloop_avg)[s.Index()];
				for(int i=0;i<4;i++)
				{
					//AA term(from VA term gamma5*V->A)
					Wpart1.glA(Wpart1_aux,i);
					Wpart2.glA(Wpart2_aux,i);
					Rpart1=Wpart1.Trace();
					Rpart2=Wpart2.Trace();
					Raux=Rpart1*Rpart2;
					_diagram[8][diff][tid]+=Raux*disVterm;//diagram 9
					Mpart1=SpinTrace(Wpart1);
					Mpart2=SpinTrace(Wpart2);
					Raux=Tr(Mpart1,Mpart2);
					_diagram[9][diff][tid]+=Raux*disVterm;//diagram 10
					//VV term(from AV term gamma5*A->V)
					Wpart1.glV(Wpart1_aux,i);
					Wpart2.glV(Wpart2_aux,i);
					Rpart1=Wpart1.Trace();
					Rpart2=Wpart2.Trace();
					Raux=Rpart1*Rpart2;
					_diagram[8][diff][tid]+=Raux*disVterm;
					Mpart1=SpinTrace(Wpart1);
					Mpart2=SpinTrace(Wpart2);
					Raux=Tr(Mpart1,Mpart2);
					_diagram[9][diff][tid]+=Raux*disVterm;
				}
				/////////////////////////////////////////////////////////mix diagram
				Mtemp1=(*lwall[tpi])[s.Index()];
				Mtemp2=(*swall[tk])[s.Index()];
				Mtemp2.hconj();
				Mtemp3=lwwprop[tpi][tpi2];
				Mtemp3.hconj();
				Mtemp4=lwwprop[tk][tpi2];
				Wpart1_aux=Mtemp2*Mtemp1*Mtemp3*Mtemp4;
				_diagram[10][diff][tid]+=Wpart1_aux.Trace();//mix 1
				Mtemp1=(*lwall[tk])[s.Index()];
				_diagram[11][diff][tid]-=Trace(Mtemp2,Mtemp1)*disVterm;//mix 2
				////second term
				Mtemp1=(*lwall[tpi2])[s.Index()];
				Mtemp3=lwwprop[tpi2][tpi];
				Mtemp3.hconj();
				Mtemp4=lwwprop[tk][tpi];
				Wpart1_aux=Mtemp2*Mtemp1*Mtemp3*Mtemp4;
				_diagram[10][diff][tid]+=Wpart1_aux.Trace();//mix 3


			}//end s loop

			//reduction
			for(int d = 0; d < num_diag; d++)
				for(int diff = 0; diff < glb_t; diff++)
					for(int tid = 0; tid < N_threads; tid++)
						diagram[delta_ind][d][diff] += _diagram[d][diff][tid] *0.5 *(1.0/ glb_t);
		}
	}

#ifdef PARALLEL
	QMP_sum_double_array((double*)diagram, 2*num_diag*N_sep*glb_t);
#endif

	char filename[512];
	sprintf(filename, "%s_ktopipi_%s.txt", common_arg->filename, label);
	SaveCorr(diagram[0][0], filename, num_diag*N_sep*glb_t);

//	omp_set_num_threads(N_threads_old);

	cost += dclock();
	print_time(cname, fname, cost);
}

