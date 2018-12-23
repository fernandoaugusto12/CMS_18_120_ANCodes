
TString fform = "2";
void eval(TString sr, TString fs, TString frmap, TString form, double *val, double *valerr){
  TFile *file = TFile::Open("zxstatisticsOS/"+sr+"_"+fs+"_"+frmap+".root");
  TH1D *hform = (TH1D*)file->Get("form"+form);
  *val    = hform->GetBinContent(1);
  *valerr = hform->GetBinError(1);
  
  return;
}


///---------------------------------------------------------------------------
void computeFinalZX(void){
  std::cout<<"Using Formula "<<fform<<std::endl;
  
  //=================== SM Higgs ======================
  double zx_smhiggs_ss_2p2p_4mu = 24.00;
  double zx_smhiggs_ss_2p2p_4mu_err = 4.90;
  double zx_smhiggs_ss_2p2p_4e = 36.00;
  double zx_smhiggs_ss_2p2p_4e_err = 6.00;
  double zx_smhiggs_ss_2p2p_2e2mu = 8.00;
  double zx_smhiggs_ss_2p2p_2e2mu_err = 2.83;
  double zx_smhiggs_ss_2p2p_2mu2e = 56.00;
  double zx_smhiggs_ss_2p2p_2mu2e_err = 7.48;

  double zx_smhiggs_os_4mu, zx_smhiggs_os_4mu_err; eval("smhiggs","4mu","nominal",fform, &zx_smhiggs_os_4mu, &zx_smhiggs_os_4mu_err);
  double zx_smhiggs_os_av_4mu, zx_smhiggs_os_av_4mu_err; eval("smhiggs","4mu","average",fform, &zx_smhiggs_os_av_4mu, &zx_smhiggs_os_av_4mu_err);
  double zx_smhiggs_os_rw_4mu, zx_smhiggs_os_rw_4mu_err; eval("smhiggs","4mu","reweight",fform, &zx_smhiggs_os_rw_4mu, &zx_smhiggs_os_rw_4mu_err);

  double zx_smhiggs_os_4e, zx_smhiggs_os_4e_err; eval("smhiggs","4e","nominal",fform, &zx_smhiggs_os_4e, &zx_smhiggs_os_4e_err);
  double zx_smhiggs_os_av_4e, zx_smhiggs_os_av_4e_err; eval("smhiggs","4e","average",fform, &zx_smhiggs_os_av_4e, &zx_smhiggs_os_av_4e_err);
  double zx_smhiggs_os_rw_4e, zx_smhiggs_os_rw_4e_err; eval("smhiggs","4e","reweight",fform, &zx_smhiggs_os_rw_4e, &zx_smhiggs_os_rw_4e_err);

  double zx_smhiggs_os_2e2mu, zx_smhiggs_os_2e2mu_err; eval("smhiggs","2e2mu","nominal",fform, &zx_smhiggs_os_2e2mu, &zx_smhiggs_os_2e2mu_err);
  double zx_smhiggs_os_av_2e2mu, zx_smhiggs_os_av_2e2mu_err; eval("smhiggs","2e2mu","average",fform, &zx_smhiggs_os_av_2e2mu, &zx_smhiggs_os_av_2e2mu_err);
  double zx_smhiggs_os_rw_2e2mu, zx_smhiggs_os_rw_2e2mu_err; eval("smhiggs","2e2mu","reweight",fform, &zx_smhiggs_os_rw_2e2mu, &zx_smhiggs_os_rw_2e2mu_err);
  
  double zx_smhiggs_os_2mu2e, zx_smhiggs_os_2mu2e_err; eval("smhiggs","2mu2e","nominal",fform, &zx_smhiggs_os_2mu2e, &zx_smhiggs_os_2mu2e_err);
  double zx_smhiggs_os_av_2mu2e, zx_smhiggs_os_av_2mu2e_err; eval("smhiggs","2mu2e","average",fform, &zx_smhiggs_os_av_2mu2e, &zx_smhiggs_os_av_2mu2e_err);
  double zx_smhiggs_os_rw_2mu2e, zx_smhiggs_os_rw_2mu2e_err; eval("smhiggs","2mu2e","reweight",fform, &zx_smhiggs_os_rw_2mu2e, &zx_smhiggs_os_rw_2mu2e_err);


  //=================== VBF-SR ======================
  double zx_vbf_ss_2p2p_4mu = 1.00;
  double zx_vbf_ss_2p2p_4mu_err = 1.00;
  double zx_vbf_ss_2p2p_4e = 1.00;
  double zx_vbf_ss_2p2p_4e_err = 1.00;
  double zx_vbf_ss_2p2p_2e2mu = 0.00;
  double zx_vbf_ss_2p2p_2e2mu_err = 0.00;
  double zx_vbf_ss_2p2p_2mu2e = 0.00;
  double zx_vbf_ss_2p2p_2mu2e_err = 0.00;

  double zx_vbf_os_4mu, zx_vbf_os_4mu_err; eval("vbf","4mu","nominal",fform, &zx_vbf_os_4mu, &zx_vbf_os_4mu_err);
  double zx_vbf_os_av_4mu, zx_vbf_os_av_4mu_err; eval("vbf","4mu","average",fform, &zx_vbf_os_av_4mu, &zx_vbf_os_av_4mu_err);
  double zx_vbf_os_rw_4mu, zx_vbf_os_rw_4mu_err; eval("vbf","4mu","reweight",fform, &zx_vbf_os_rw_4mu, &zx_vbf_os_rw_4mu_err);

  double zx_vbf_os_4e, zx_vbf_os_4e_err; eval("vbf","4e","nominal",fform, &zx_vbf_os_4e, &zx_vbf_os_4e_err);
  double zx_vbf_os_av_4e, zx_vbf_os_av_4e_err; eval("vbf","4e","average",fform, &zx_vbf_os_av_4e, &zx_vbf_os_av_4e_err);
  double zx_vbf_os_rw_4e, zx_vbf_os_rw_4e_err; eval("vbf","4e","reweight",fform, &zx_vbf_os_rw_4e, &zx_vbf_os_rw_4e_err);

  double zx_vbf_os_2e2mu, zx_vbf_os_2e2mu_err; eval("vbf","2e2mu","nominal",fform, &zx_vbf_os_2e2mu, &zx_vbf_os_2e2mu_err);
  double zx_vbf_os_av_2e2mu, zx_vbf_os_av_2e2mu_err; eval("vbf","2e2mu","average",fform, &zx_vbf_os_av_2e2mu, &zx_vbf_os_av_2e2mu_err);
  double zx_vbf_os_rw_2e2mu, zx_vbf_os_rw_2e2mu_err; eval("vbf","2e2mu","reweight",fform, &zx_vbf_os_rw_2e2mu, &zx_vbf_os_rw_2e2mu_err);
  
  double zx_vbf_os_2mu2e, zx_vbf_os_2mu2e_err; eval("vbf","2mu2e","nominal",fform, &zx_vbf_os_2mu2e, &zx_vbf_os_2mu2e_err);
  double zx_vbf_os_av_2mu2e, zx_vbf_os_av_2mu2e_err; eval("vbf","2mu2e","average",fform, &zx_vbf_os_av_2mu2e, &zx_vbf_os_av_2mu2e_err);
  double zx_vbf_os_rw_2mu2e, zx_vbf_os_rw_2mu2e_err; eval("vbf","2mu2e","reweight",fform, &zx_vbf_os_rw_2mu2e, &zx_vbf_os_rw_2mu2e_err);
  
  
  
  
  
  std::cout<<"==========================             SM Higgs             ==================================================="<<std::endl;
  std::cout<<"========================== Final numbers for Z+X estimation ==================================================="<<std::endl;
  std::cout<<Form(" Source     | 4mu\t\t| 4e\t\t| 2e2mu\t\t| 2mu2e\t\t| 2l2l\t\t| 4l")<<std::endl;
  std::cout<<"---------------------------------------------------------------------------------------------------------------"<<std::endl;
  double zx_smhiggs_ss_2p2p_2l2l = zx_smhiggs_ss_2p2p_2e2mu + zx_smhiggs_ss_2p2p_2mu2e;
  double zx_smhiggs_ss_2p2p_2l2l_err = sqrt( pow(zx_smhiggs_ss_2p2p_2e2mu_err,2) + pow(zx_smhiggs_ss_2p2p_2mu2e_err,2) );
  double zx_smhiggs_ss_2p2p_4l = zx_smhiggs_ss_2p2p_4mu + zx_smhiggs_ss_2p2p_4e + zx_smhiggs_ss_2p2p_2e2mu + zx_smhiggs_ss_2p2p_2mu2e;
  double zx_smhiggs_ss_2p2p_4l_err = sqrt( pow(zx_smhiggs_ss_2p2p_4mu_err,2) + pow(zx_smhiggs_ss_2p2p_4e_err,2) + pow(zx_smhiggs_ss_2p2p_2e2mu_err,2) + pow(zx_smhiggs_ss_2p2p_2mu2e_err,2) );
  std::cout<<Form(" 2P2P       | %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f",
		  zx_smhiggs_ss_2p2p_4mu,zx_smhiggs_ss_2p2p_4mu_err,zx_smhiggs_ss_2p2p_4e,zx_smhiggs_ss_2p2p_4e_err,zx_smhiggs_ss_2p2p_2e2mu,zx_smhiggs_ss_2p2p_2e2mu_err,
		  zx_smhiggs_ss_2p2p_2mu2e,zx_smhiggs_ss_2p2p_2mu2e_err,zx_smhiggs_ss_2p2p_2l2l,zx_smhiggs_ss_2p2p_2l2l_err,zx_smhiggs_ss_2p2p_4l,zx_smhiggs_ss_2p2p_4l_err)<<std::endl;
  std::cout<<"---------------------------------------------------------------------------------------------------------------"<<std::endl;		  
  
  double zx_smhiggs_os_2l2l = zx_smhiggs_os_2e2mu + zx_smhiggs_os_2mu2e;
  double zx_smhiggs_os_2l2l_err = sqrt( pow(zx_smhiggs_os_2e2mu_err,2) + pow(zx_smhiggs_os_2mu2e_err,2) );
  double zx_smhiggs_os_4l = zx_smhiggs_os_4mu + zx_smhiggs_os_4e + zx_smhiggs_os_2e2mu + zx_smhiggs_os_2mu2e;
  double zx_smhiggs_os_4l_err = sqrt( pow(zx_smhiggs_os_4mu_err,2) + pow(zx_smhiggs_os_4e_err,2) + pow(zx_smhiggs_os_2e2mu_err,2) + pow(zx_smhiggs_os_2mu2e_err,2) );
  std::cout<<Form(" OS         | %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f",
		  zx_smhiggs_os_4mu,zx_smhiggs_os_4mu_err,zx_smhiggs_os_4e,zx_smhiggs_os_4e_err,zx_smhiggs_os_2e2mu,zx_smhiggs_os_2e2mu_err,
		  zx_smhiggs_os_2mu2e,zx_smhiggs_os_2mu2e_err,zx_smhiggs_os_2l2l,zx_smhiggs_os_2l2l_err,zx_smhiggs_os_4l,zx_smhiggs_os_4l_err)<<std::endl;
  std::cout<<"---------------------------------------------------------------------------------------------------------------"<<std::endl;
		  
  double zx_smhiggs_os_av_2l2l = zx_smhiggs_os_av_2e2mu + zx_smhiggs_os_av_2mu2e;
  double zx_smhiggs_os_av_2l2l_err = sqrt( pow(zx_smhiggs_os_av_2e2mu_err,2) + pow(zx_smhiggs_os_av_2mu2e_err,2) );
  double zx_smhiggs_os_av_4l = zx_smhiggs_os_av_4mu + zx_smhiggs_os_av_4e + zx_smhiggs_os_av_2e2mu + zx_smhiggs_os_av_2mu2e;
  double zx_smhiggs_os_av_4l_err = sqrt( pow(zx_smhiggs_os_av_4mu_err,2) + pow(zx_smhiggs_os_av_4e_err,2) + pow(zx_smhiggs_os_av_2e2mu_err,2) + pow(zx_smhiggs_os_av_2mu2e_err,2) );
  std::cout<<Form(" Weighted   | %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f",
		  zx_smhiggs_os_av_4mu,zx_smhiggs_os_av_4mu_err,zx_smhiggs_os_av_4e,zx_smhiggs_os_av_4e_err,zx_smhiggs_os_av_2e2mu,zx_smhiggs_os_av_2e2mu_err,
		  zx_smhiggs_os_av_2mu2e,zx_smhiggs_os_av_2mu2e_err,zx_smhiggs_os_av_2l2l,zx_smhiggs_os_av_2l2l_err,zx_smhiggs_os_av_4l,zx_smhiggs_os_av_4l_err)<<std::endl;
  std::cout<<"---------------------------------------------------------------------------------------------------------------"<<std::endl;
  
  double zx_smhiggs_os_rw_2l2l = zx_smhiggs_os_rw_2e2mu + zx_smhiggs_os_rw_2mu2e;
  double zx_smhiggs_os_rw_2l2l_err = sqrt( pow(zx_smhiggs_os_rw_2e2mu_err,2) + pow(zx_smhiggs_os_rw_2mu2e_err,2) );
  double zx_smhiggs_os_rw_4l = zx_smhiggs_os_rw_4mu + zx_smhiggs_os_rw_4e + zx_smhiggs_os_rw_2e2mu + zx_smhiggs_os_rw_2mu2e;
  double zx_smhiggs_os_rw_4l_err = sqrt( pow(zx_smhiggs_os_rw_4mu_err,2) + pow(zx_smhiggs_os_rw_4e_err,2) + pow(zx_smhiggs_os_rw_2e2mu_err,2) + pow(zx_smhiggs_os_rw_2mu2e_err,2) );
  std::cout<<Form(" Reweighted | %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f",
		  zx_smhiggs_os_rw_4mu,zx_smhiggs_os_rw_4mu_err,zx_smhiggs_os_rw_4e,zx_smhiggs_os_rw_4e_err,zx_smhiggs_os_rw_2e2mu,zx_smhiggs_os_rw_2e2mu_err,
		  zx_smhiggs_os_rw_2mu2e,zx_smhiggs_os_rw_2mu2e_err,zx_smhiggs_os_rw_2l2l,zx_smhiggs_os_rw_2l2l_err,zx_smhiggs_os_rw_4l,zx_smhiggs_os_rw_4l_err)<<std::endl;
  std::cout<<"---------------------------------------------------------------------------------------------------------------"<<std::endl;
  
  double zx_smhiggs_sys_4mu = fabs(zx_smhiggs_os_av_4mu-zx_smhiggs_os_rw_4mu);
  double zx_smhiggs_sys_4mu_err = sqrt( pow(zx_smhiggs_os_av_4mu_err,2) + pow(zx_smhiggs_os_rw_4mu_err,2) );
  double zx_smhiggs_sys_4e = fabs(zx_smhiggs_os_av_4e-zx_smhiggs_os_rw_4e);
  double zx_smhiggs_sys_4e_err = sqrt( pow(zx_smhiggs_os_av_4e_err,2) + pow(zx_smhiggs_os_rw_4e_err,2) );
  double zx_smhiggs_sys_2e2mu = fabs(zx_smhiggs_os_av_2e2mu-zx_smhiggs_os_rw_2e2mu);
  double zx_smhiggs_sys_2e2mu_err = sqrt( pow(zx_smhiggs_os_av_2e2mu_err,2) + pow(zx_smhiggs_os_rw_2e2mu_err,2) );
  double zx_smhiggs_sys_2mu2e = fabs(zx_smhiggs_os_av_2mu2e-zx_smhiggs_os_rw_2mu2e);
  double zx_smhiggs_sys_2mu2e_err = sqrt( pow(zx_smhiggs_os_av_2mu2e_err,2) + pow(zx_smhiggs_os_rw_2mu2e_err,2) );
  double zx_smhiggs_sys_2l2l = fabs(zx_smhiggs_os_av_2l2l-zx_smhiggs_os_rw_2l2l);
  double zx_smhiggs_sys_2l2l_err = sqrt( pow(zx_smhiggs_os_av_2l2l_err,2) + pow(zx_smhiggs_os_rw_2l2l_err,2) );
  double zx_smhiggs_sys_4l = sqrt(pow(zx_smhiggs_sys_4mu,2) + pow(zx_smhiggs_sys_4e,2) + pow(zx_smhiggs_sys_2l2l,2));
  double zx_smhiggs_sys_4l_err = sqrt( pow(zx_smhiggs_sys_4mu_err,2) + pow(zx_smhiggs_sys_4e_err,2) + pow(zx_smhiggs_sys_2l2l_err,2) );
  std::cout<<Form(" Syst       | %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f",
		  zx_smhiggs_sys_4mu,zx_smhiggs_sys_4mu_err,zx_smhiggs_sys_4e,zx_smhiggs_sys_4e_err,zx_smhiggs_sys_2e2mu,zx_smhiggs_sys_2e2mu_err,
		  zx_smhiggs_sys_2mu2e,zx_smhiggs_sys_2mu2e_err,zx_smhiggs_sys_2l2l,zx_smhiggs_sys_2l2l_err,zx_smhiggs_sys_4l,zx_smhiggs_sys_4l_err)<<std::endl;
  std::cout<<"==============================================================================================================="<<std::endl;
		  
  
  std::cout<<"\n\n===================================== Final numbers for Z+X estimation (stat/syst/total uncertainty) =================================================="<<std::endl;
  std::cout<<Form(" Source | 4mu\t\t\t| 4e\t\t\t| 2e2mu\t\t\t| 2mu2e\t\t\t| 2l2l\t\t\t| 4l")<<std::endl;
  std::cout<<"-------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
  std::cout<<Form(" OS     | %.2f+/-(%.2f/%.2f/%.2f)\t| %.2f+/-(%.2f/%.2f/%.2f)\t| %.2f+/-(%.2f/%.2f/%.2f)\t| %.2f+/-(%.2f/%.2f/%.2f)\t| %.2f+/-(%.2f/%.2f/%.2f)\t| %.2f+/-(%.2f/%.2f/%.2f)",
		  zx_smhiggs_os_4mu,zx_smhiggs_os_4mu_err,zx_smhiggs_sys_4mu,sqrt(pow(zx_smhiggs_os_4mu_err,2)+pow(zx_smhiggs_sys_4mu,2)),
		  zx_smhiggs_os_4e,zx_smhiggs_os_4e_err,zx_smhiggs_sys_4e,sqrt(pow(zx_smhiggs_os_4e_err,2)+pow(zx_smhiggs_sys_4e,2)),
		  zx_smhiggs_os_2e2mu,zx_smhiggs_os_2e2mu_err,zx_smhiggs_sys_2e2mu,sqrt(pow(zx_smhiggs_os_2e2mu_err,2)+pow(zx_smhiggs_sys_2e2mu,2)),
		  zx_smhiggs_os_2mu2e,zx_smhiggs_os_2mu2e_err,zx_smhiggs_sys_2mu2e,sqrt(pow(zx_smhiggs_os_2mu2e_err,2)+pow(zx_smhiggs_sys_2mu2e,2)),
		  zx_smhiggs_os_2l2l,zx_smhiggs_os_2l2l_err,zx_smhiggs_sys_2l2l,sqrt(pow(zx_smhiggs_os_2l2l_err,2)+pow(zx_smhiggs_sys_2l2l,2)),
		  zx_smhiggs_os_4l,zx_smhiggs_os_4l_err,zx_smhiggs_sys_4l,sqrt(pow(zx_smhiggs_os_4l_err,2)+pow(zx_smhiggs_sys_4l,2)))<<std::endl;  
  std::cout<<"======================================================================================================================================================="<<std::endl;


  std::cout<<"\n\n"<<std::endl;
  std::cout<<"==========================              VBF SR              ==================================================="<<std::endl;
  std::cout<<"========================== Final numbers for Z+X estimation ==================================================="<<std::endl;
  std::cout<<Form(" Source     | 4mu\t\t| 4e\t\t| 2e2mu\t\t| 2mu2e\t\t| 2l2l\t\t| 4l")<<std::endl;
  std::cout<<"---------------------------------------------------------------------------------------------------------------"<<std::endl;
  double zx_vbf_ss_2p2p_2l2l = zx_vbf_ss_2p2p_2e2mu + zx_vbf_ss_2p2p_2mu2e;
  double zx_vbf_ss_2p2p_2l2l_err = sqrt( pow(zx_vbf_ss_2p2p_2e2mu_err,2) + pow(zx_vbf_ss_2p2p_2mu2e_err,2) );
  double zx_vbf_ss_2p2p_4l = zx_vbf_ss_2p2p_4mu + zx_vbf_ss_2p2p_4e + zx_vbf_ss_2p2p_2e2mu + zx_vbf_ss_2p2p_2mu2e;
  double zx_vbf_ss_2p2p_4l_err = sqrt( pow(zx_vbf_ss_2p2p_4mu_err,2) + pow(zx_vbf_ss_2p2p_4e_err,2) + pow(zx_vbf_ss_2p2p_2e2mu_err,2) + pow(zx_vbf_ss_2p2p_2mu2e_err,2) );
  std::cout<<Form(" 2P2P       | %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f",
		  zx_vbf_ss_2p2p_4mu,zx_vbf_ss_2p2p_4mu_err,zx_vbf_ss_2p2p_4e,zx_vbf_ss_2p2p_4e_err,zx_vbf_ss_2p2p_2e2mu,zx_vbf_ss_2p2p_2e2mu_err,
		  zx_vbf_ss_2p2p_2mu2e,zx_vbf_ss_2p2p_2mu2e_err,zx_vbf_ss_2p2p_2l2l,zx_vbf_ss_2p2p_2l2l_err,zx_vbf_ss_2p2p_4l,zx_vbf_ss_2p2p_4l_err)<<std::endl;
  std::cout<<"---------------------------------------------------------------------------------------------------------------"<<std::endl;		  
  
  double zx_vbf_os_2l2l = zx_vbf_os_2e2mu + zx_vbf_os_2mu2e;
  double zx_vbf_os_2l2l_err = sqrt( pow(zx_vbf_os_2e2mu_err,2) + pow(zx_vbf_os_2mu2e_err,2) );
  double zx_vbf_os_4l = zx_vbf_os_4mu + zx_vbf_os_4e + zx_vbf_os_2e2mu + zx_vbf_os_2mu2e;
  double zx_vbf_os_4l_err = sqrt( pow(zx_vbf_os_4mu_err,2) + pow(zx_vbf_os_4e_err,2) + pow(zx_vbf_os_2e2mu_err,2) + pow(zx_vbf_os_2mu2e_err,2) );
  std::cout<<Form(" OS         | %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f",
		  zx_vbf_os_4mu,zx_vbf_os_4mu_err,zx_vbf_os_4e,zx_vbf_os_4e_err,zx_vbf_os_2e2mu,zx_vbf_os_2e2mu_err,
		  zx_vbf_os_2mu2e,zx_vbf_os_2mu2e_err,zx_vbf_os_2l2l,zx_vbf_os_2l2l_err,zx_vbf_os_4l,zx_vbf_os_4l_err)<<std::endl;
  std::cout<<"---------------------------------------------------------------------------------------------------------------"<<std::endl;
		  
  double zx_vbf_os_av_2l2l = zx_vbf_os_av_2e2mu + zx_vbf_os_av_2mu2e;
  double zx_vbf_os_av_2l2l_err = sqrt( pow(zx_vbf_os_av_2e2mu_err,2) + pow(zx_vbf_os_av_2mu2e_err,2) );
  double zx_vbf_os_av_4l = zx_vbf_os_av_4mu + zx_vbf_os_av_4e + zx_vbf_os_av_2e2mu + zx_vbf_os_av_2mu2e;
  double zx_vbf_os_av_4l_err = sqrt( pow(zx_vbf_os_av_4mu_err,2) + pow(zx_vbf_os_av_4e_err,2) + pow(zx_vbf_os_av_2e2mu_err,2) + pow(zx_vbf_os_av_2mu2e_err,2) );
  std::cout<<Form(" Weighted   | %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f",
		  zx_vbf_os_av_4mu,zx_vbf_os_av_4mu_err,zx_vbf_os_av_4e,zx_vbf_os_av_4e_err,zx_vbf_os_av_2e2mu,zx_vbf_os_av_2e2mu_err,
		  zx_vbf_os_av_2mu2e,zx_vbf_os_av_2mu2e_err,zx_vbf_os_av_2l2l,zx_vbf_os_av_2l2l_err,zx_vbf_os_av_4l,zx_vbf_os_av_4l_err)<<std::endl;
  std::cout<<"---------------------------------------------------------------------------------------------------------------"<<std::endl;
  
  double zx_vbf_os_rw_2l2l = zx_vbf_os_rw_2e2mu + zx_vbf_os_rw_2mu2e;
  double zx_vbf_os_rw_2l2l_err = sqrt( pow(zx_vbf_os_rw_2e2mu_err,2) + pow(zx_vbf_os_rw_2mu2e_err,2) );
  double zx_vbf_os_rw_4l = zx_vbf_os_rw_4mu + zx_vbf_os_rw_4e + zx_vbf_os_rw_2e2mu + zx_vbf_os_rw_2mu2e;
  double zx_vbf_os_rw_4l_err = sqrt( pow(zx_vbf_os_rw_4mu_err,2) + pow(zx_vbf_os_rw_4e_err,2) + pow(zx_vbf_os_rw_2e2mu_err,2) + pow(zx_vbf_os_rw_2mu2e_err,2) );
  std::cout<<Form(" Reweighted | %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f",
		  zx_vbf_os_rw_4mu,zx_vbf_os_rw_4mu_err,zx_vbf_os_rw_4e,zx_vbf_os_rw_4e_err,zx_vbf_os_rw_2e2mu,zx_vbf_os_rw_2e2mu_err,
		  zx_vbf_os_rw_2mu2e,zx_vbf_os_rw_2mu2e_err,zx_vbf_os_rw_2l2l,zx_vbf_os_rw_2l2l_err,zx_vbf_os_rw_4l,zx_vbf_os_rw_4l_err)<<std::endl;
  std::cout<<"---------------------------------------------------------------------------------------------------------------"<<std::endl;
  
  double zx_vbf_sys_4mu = fabs(zx_vbf_os_av_4mu-zx_vbf_os_rw_4mu);
  double zx_vbf_sys_4mu_err = sqrt( pow(zx_vbf_os_av_4mu_err,2) + pow(zx_vbf_os_rw_4mu_err,2) );
  double zx_vbf_sys_4e = fabs(zx_vbf_os_av_4e-zx_vbf_os_rw_4e);
  double zx_vbf_sys_4e_err = sqrt( pow(zx_vbf_os_av_4e_err,2) + pow(zx_vbf_os_rw_4e_err,2) );
  double zx_vbf_sys_2e2mu = fabs(zx_vbf_os_av_2e2mu-zx_vbf_os_rw_2e2mu);
  double zx_vbf_sys_2e2mu_err = sqrt( pow(zx_vbf_os_av_2e2mu_err,2) + pow(zx_vbf_os_rw_2e2mu_err,2) );
  double zx_vbf_sys_2mu2e = fabs(zx_vbf_os_av_2mu2e-zx_vbf_os_rw_2mu2e);
  double zx_vbf_sys_2mu2e_err = sqrt( pow(zx_vbf_os_av_2mu2e_err,2) + pow(zx_vbf_os_rw_2mu2e_err,2) );
  double zx_vbf_sys_2l2l = fabs(zx_vbf_os_av_2l2l-zx_vbf_os_rw_2l2l);
  double zx_vbf_sys_2l2l_err = sqrt( pow(zx_vbf_os_av_2l2l_err,2) + pow(zx_vbf_os_rw_2l2l_err,2) );
  double zx_vbf_sys_4l = sqrt(pow(zx_vbf_sys_4mu,2) + pow(zx_vbf_sys_4e,2) + pow(zx_vbf_sys_2l2l,2));
  double zx_vbf_sys_4l_err = sqrt( pow(zx_vbf_sys_4mu_err,2) + pow(zx_vbf_sys_4e_err,2) + pow(zx_vbf_sys_2l2l_err,2) );
  std::cout<<Form(" Syst       | %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f\t| %.2f+/-%.2f",
		  zx_vbf_sys_4mu,zx_vbf_sys_4mu_err,zx_vbf_sys_4e,zx_vbf_sys_4e_err,zx_vbf_sys_2e2mu,zx_vbf_sys_2e2mu_err,
		  zx_vbf_sys_2mu2e,zx_vbf_sys_2mu2e_err,zx_vbf_sys_2l2l,zx_vbf_sys_2l2l_err,zx_vbf_sys_4l,zx_vbf_sys_4l_err)<<std::endl;
  std::cout<<"==============================================================================================================="<<std::endl;
		  
  
  std::cout<<"\n\n===================================== Final numbers for Z+X estimation (stat/syst/total uncertainty) =================================================="<<std::endl;
  std::cout<<Form(" Source | 4mu\t\t\t| 4e\t\t\t| 2e2mu\t\t\t| 2mu2e\t\t\t| 2l2l\t\t\t| 4l")<<std::endl;
  std::cout<<"-------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
  std::cout<<Form(" OS     | %.2f+/-(%.2f/%.2f/%.2f)\t| %.2f+/-(%.2f/%.2f/%.2f)\t| %.2f+/-(%.2f/%.2f/%.2f)\t| %.2f+/-(%.2f/%.2f/%.2f)\t| %.2f+/-(%.2f/%.2f/%.2f)\t| %.2f+/-(%.2f/%.2f/%.2f)",
		  zx_vbf_os_4mu,zx_vbf_os_4mu_err,zx_vbf_sys_4mu,sqrt(pow(zx_vbf_os_4mu_err,2)+pow(zx_vbf_sys_4mu,2)),
		  zx_vbf_os_4e,zx_vbf_os_4e_err,zx_vbf_sys_4e,sqrt(pow(zx_vbf_os_4e_err,2)+pow(zx_vbf_sys_4e,2)),
		  zx_vbf_os_2e2mu,zx_vbf_os_2e2mu_err,zx_vbf_sys_2e2mu,sqrt(pow(zx_vbf_os_2e2mu_err,2)+pow(zx_vbf_sys_2e2mu,2)),
		  zx_vbf_os_2mu2e,zx_vbf_os_2mu2e_err,zx_vbf_sys_2mu2e,sqrt(pow(zx_vbf_os_2mu2e_err,2)+pow(zx_vbf_sys_2mu2e,2)),
		  zx_vbf_os_2l2l,zx_vbf_os_2l2l_err,zx_vbf_sys_2l2l,sqrt(pow(zx_vbf_os_2l2l_err,2)+pow(zx_vbf_sys_2l2l,2)),
		  zx_vbf_os_4l,zx_vbf_os_4l_err,zx_vbf_sys_4l,sqrt(pow(zx_vbf_os_4l_err,2)+pow(zx_vbf_sys_4l,2)))<<std::endl;  
  std::cout<<"======================================================================================================================================================="<<std::endl;

  
  
  return;
}
