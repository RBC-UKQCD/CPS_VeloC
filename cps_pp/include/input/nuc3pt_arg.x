class Nuc3ptArg {
	string cname<>;
	int num_masses;
	Float mass[10];
	CgArg cg;
	int t_source;
	int num_src;
	int source_inc;
	int t_sink;
	int BoxStart;
	int BoxEnd;
	int gauss_N;
	Float gauss_W;
	int x[3];
	Float theta;
	int load_u1_lat;
	Float *u1_gauge_ptr;
	SourceType src_type;
	int DoUnPolarized;
	int DoUnPolarizedMom;
	int DoPolarized;
	int DoPolarizedMom;
	int DoHalfFermion;
	int DoPerPlusAper;
	int MaxMom2;
	int DoSS2ptF;
	int DoGa1Proj;
	int DoConserved;
	int num_mult;
	int mt[5];
	GaussianKernelLinkSmearType gauss_link_smear_type;
	int gauss_link_smear_N;
	Float gauss_link_smear_coeff;
	CalcQpropType calc_QProp;
	CalcSeqType calc_seqQ;
	string prop_file<>;
	string u_seq_prop_file<>;
	string d_seq_prop_file<>;
	string ptsrc_prop_file<>;
	string ensemble_label<>;
	int ensemble_id;
	int StartSrcSpin;
	int EndSrcSpin;
	int StartSrcColor;
	int EndSrcColor;
	int DoDisconnected;

	memfun   Nuc3ptArg (  ) ;
	memfun   void check_args (  ) ;
	memfun   int NumMasses (  ) ;
	memfun   void NumMasses (  int n ) ;
	memfun   Float Mass (  int m ) ;
};