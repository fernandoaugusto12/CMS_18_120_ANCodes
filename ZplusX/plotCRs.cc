#define nShifts 26 //This corresponds to how many uncertainty shifts have been done to the NNs

#define elesys 0.050
#define musys  0.135
double FRsys(double fr, double lpdgid){
  double syserr;
  if(abs(lpdgid) == 11) syserr = elesys;
  else if(abs(lpdgid) == 13) syserr = musys;
  else{
    std::cout<<"!!ERROR: PDGID "<<lpdgid<<" NOT EXPECTED!!"<<std::endl;
    return -1;
  }

  return (syserr*fr*sqrt(1 - 2*fr*(1-fr)))/pow(1-fr,2);
}

double proderr(double a, double ae, double b, double be){
  return a*b*sqrt( pow(ae/a,2) + pow(be/b,2) );
}

//Code to extract FR and plot it
void analyzeCRs(TString SR, TString fstate, TString frmap){
  std::cout<<"Running for "<<SR<<", "<<fstate<<", "<<frmap<<std::endl;
  
  gROOT->Reset();
  gROOT->SetBatch();
  TH1::SetDefaultSumw2();
  
  //Picks the FR map
  TFile *FRfile;
  TH2D *ElectronPtEtaMapFR;
  TH2D *MuonPtEtaMapFR;

  if(frmap == "nominal"){
    FRfile = TFile::Open("FakeRatePtEtaMap.root");
    ElectronPtEtaMapFR = (TH2D*)FRfile->Get("ElectronPtEtaMapFR");
    MuonPtEtaMapFR  = (TH2D*)FRfile->Get("MuonPtEtaMapFR"); 
  }
  if(frmap == "average"){
    FRfile = TFile::Open("FakeRatePtEtaMapAverage.root");
    ElectronPtEtaMapFR = (TH2D*)FRfile->Get("ElectronPtEtaMapFR");
    MuonPtEtaMapFR = (TH2D*)FRfile->Get("MuonPtEtaMapFR");
  }
  if(frmap == "reweight"){
    FRfile = TFile::Open("FakeRatePtEtaMapReweight.root");
    ElectronPtEtaMapFR = (TH2D*)FRfile->Get("ElectronPtEtaMapFR");
    MuonPtEtaMapFR = (TH2D*)FRfile->Get("MuonPtEtaMapFR");
  }
  
  //List of files
  std::vector<TString> files = {"output_ZplusX_DataRunIISummer16_03Feb2017_withKerasDiscriminants.root",
				"output_ZplusX_GluGluToZZTo4LplusZZJJTo4LEWKplusZZTo4LplusZZTo2L2Nu_13TeV_madgraph_powheg_pythia8.root",
				"output_ZplusX_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root",
				"output_ZplusX_TTTo2L2Nu_noSC_TuneCUETP8M2T4_13TeV-powheg-pythia8.root",
				"output_ZplusX_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"
  };

  int nbins=1000;
  double xmin=0, xmax=1000;
  if(SR == "vbf"){
    nbins = 14;
    xmin = 117;
    xmax = 131;
  }
  
  TString final_state = "";
  if(fstate == "4mu") final_state = "4#mu";
  else if(fstate == "2mu2e") final_state = "2#mu2e";
  else if(fstate == "2e2mu") final_state = "2e2#mu";
  else final_state = fstate;
  
  TH1D *hm4l_data_2p2f = new TH1D("hm4l_data_2p2f","",nbins,xmin,xmax);
  TH1D *hm4l_data_3p1f = new TH1D("hm4l_data_3p1f","",nbins,xmin,xmax);
  
  TH1D *hm4l_zz_2p2f = new TH1D("hm4l_zz_2p2f","",nbins,xmin,xmax);
  TH1D *hm4l_zz_3p1f = new TH1D("hm4l_zz_3p1f","",nbins,xmin,xmax);
  hm4l_zz_2p2f->SetFillColor(kAzure-4);
  hm4l_zz_3p1f->SetFillColor(kAzure-4);
  TH1D *hm4l_zz_3p1f_eq = new TH1D("hm4l_zz_3p1f_eq","",nbins,xmin,xmax);
  
  TH1D *hm4l_wz_2p2f = new TH1D("hm4l_wz_2p2f","",nbins,xmin,xmax);
  TH1D *hm4l_wz_3p1f = new TH1D("hm4l_wz_3p1f","",nbins,xmin,xmax);
  hm4l_wz_2p2f->SetFillColor(kMagenta+1);
  hm4l_wz_3p1f->SetFillColor(kMagenta+1);
  
  TH1D *hm4l_ttjets_2p2f = new TH1D("hm4l_ttjets_2p2f","",nbins,xmin,xmax);
  TH1D *hm4l_ttjets_3p1f = new TH1D("hm4l_ttjets_3p1f","",nbins,xmin,xmax);
  hm4l_ttjets_2p2f->SetFillColor(kBlue-2);
  hm4l_ttjets_3p1f->SetFillColor(kBlue-2);
  
  TH1D *hm4l_zjets_2p2f = new TH1D("hm4l_zjets_2p2f","",nbins,xmin,xmax);
  TH1D *hm4l_zjets_3p1f = new TH1D("hm4l_zjets_3p1f","",nbins,xmin,xmax);
  hm4l_zjets_2p2f->SetFillColor(kTeal-8);
  hm4l_zjets_3p1f->SetFillColor(kTeal-8);
  
  THStack *hm4lStack_2p2f = new THStack();
  THStack *hm4lStack_3p1f = new THStack();
  
  TH1D *hm4l_3p1f_from_2p2f = new TH1D("hm4l_3p1f_from_2p2f","",nbins,xmin,xmax);
  hm4l_3p1f_from_2p2f->SetLineColor(kRed);
  hm4l_3p1f_from_2p2f->SetFillStyle(0);
  
  TH1D *hm4l_3p1f_from_2p2f_eq = new TH1D("hm4l_3p1f_from_2p2f_eq","",nbins,xmin,xmax);
  TH1D *hm4l_3p1f_from_2p2f_eq_sys = new TH1D("hm4l_3p1f_from_2p2f_eq_sys","",nbins,xmin,xmax);
  TH1D *hm4l_3p1f = new TH1D("hm4l_3p1f","",nbins,xmin,xmax);  
  TH1D *hm4l_3p1f_sys = new TH1D("hm4l_3p1f_sys","",nbins,xmin,xmax);  
  
  //------------------------ For NN shape ----------------------------------------
  std::vector<TString> shift_name, shifts1, shifts2, shifts3; 
  shifts1 = {"","_l1ptDown","_l1ptUp","_l2ptDown","_l2ptUp","_l3ptDown","_l3ptUp","_l4ptDown","_l4ptUp",
	     "_j1ptDown","_j1ptUp","_j2ptDown","_j2ptUp","_j3ptDown","_j3ptUp",
	     "_metEleEnDown","_metEleEnUp","_metJetEnDown","_metJetEnUp","_metJetResDown","_metJetResUp",
	     "_metMuEnDown","_metMuEnUp","_metPhoEnDown","_metPhoEnUp","_metUncEnDown","_metUncEnUp"};
  shifts2 = {"_FRUp"};
  shifts3 = {"_FRDown"};
  if(frmap == "nominal") shift_name = shifts1;
  if(frmap == "average") shift_name = shifts2;
  if(frmap == "reweight") shift_name = shifts3;
  const int nshifts = shift_name.size();

  int nbinsx=120;
  double nnxmin=0, nnxmax=1;

  TH1D *hk57nj2_data_3p1f[nshifts];
  TH1D *hk57nj2_2p2f_3p1f[nshifts];
  
  TH1D *hk24nj3_data_3p1f[nshifts];
  TH1D *hk24nj3_2p2f_3p1f[nshifts];

  TH2D *hk57nj2_deltajj_data_3p1f[nshifts];
  TH2D *hk57nj2_deltajj_2p2f_3p1f[nshifts];
  
  TH2D *hk24nj3_deltajj_data_3p1f[nshifts];
  TH2D *hk24nj3_deltajj_2p2f_3p1f[nshifts];
  
  if(nshifts != shift_name.size()) std::cout<<"NUMBER OF SHIFTS != NUMBER OF NAMES"<<std::endl;
  for(int k=0; k<nshifts; ++k){
    hk57nj2_data_3p1f[k] = new TH1D("k57nj2_"+fstate+shift_name[k],"",nbinsx,nnxmin,nnxmax);
    hk57nj2_2p2f_3p1f[k] = new TH1D("k57nj2__"+fstate+shift_name[k],"",nbinsx,nnxmin,nnxmax);

    hk24nj3_data_3p1f[k] = new TH1D("k24nj3_"+fstate+shift_name[k],"",nbinsx,nnxmin,nnxmax);
    hk24nj3_2p2f_3p1f[k] = new TH1D("k24nj3__"+fstate+shift_name[k],"",nbinsx,nnxmin,nnxmax);
    
    hk57nj2_deltajj_data_3p1f[k] = new TH2D("k57nj2_deltajj_"+fstate+shift_name[k],"",12,nnxmin,nnxmax,10,0,9.4);
    hk57nj2_deltajj_2p2f_3p1f[k] = new TH2D("k57nj2_deltajj__"+fstate+shift_name[k],"",12,nnxmin,nnxmax,10,0,9.4);

    hk24nj3_deltajj_data_3p1f[k] = new TH2D("k24nj3_deltajj_"+fstate+shift_name[k],"",12,nnxmin,nnxmax,10,0,9.4);
    hk24nj3_deltajj_2p2f_3p1f[k] = new TH2D("k24nj3_deltajj__"+fstate+shift_name[k],"",12,nnxmin,nnxmax,10,0,9.4);
  }
  
  float sum1_2p2f=0, sum2_2p2f=0;
  Float_t f_weight=-999, f_mass4l=-999;
  Int_t f_2p2f=-999, f_3p1f=-999, f_2p2p=-999, f_Nbjets=-999, f_njets_pass=-999;
  Float_t f_lept1_pt=-999, f_lept2_pt=-999, f_lept3_pt=-999, f_lept4_pt=-999;
  Float_t f_lept1_eta=-999, f_lept2_eta=-999, f_lept3_eta=-999, f_lept4_eta=-999;
  Int_t f_lept1_pass=-999, f_lept2_pass=-999, f_lept3_pass=-999, f_lept4_pass=-999;
  Int_t f_lept1_pdgid=-999, f_lept2_pdgid=-999, f_lept3_pdgid=-999, f_lept4_pdgid=-999;
  Float_t f_k57nj2[nShifts], f_k24nj3[nShifts];
  Float_t f_jets_highpt_eta[4];
  
  double N3p1f_FRsys_sumUp=0, N3p1f_2p2f_FRsys_sumUp=0;
  double N3p1f_FRsys_sumDown=0, N3p1f_2p2f_FRsys_sumDown=0;
  double N3p1f_FRsys2_sum=0, N3p1f_2p2f_FRsys2_sum=0;
  
  for(unsigned int ifile=0; ifile<(unsigned int)files.size(); ++ifile){
    TFile *CRfile = TFile::Open("../Merged/"+files[ifile]);
    TTree *CRtree = (TTree*)CRfile->Get("ControlRegions");
    unsigned int nCR = CRtree->GetEntries();
    
    CRtree->SetBranchAddress("f_weight", &f_weight);
    CRtree->SetBranchAddress("f_mass4l", &f_mass4l);
    CRtree->SetBranchAddress("f_2p2f", &f_2p2f);
    CRtree->SetBranchAddress("f_3p1f", &f_3p1f);
    CRtree->SetBranchAddress("f_2p2p", &f_2p2p);
    CRtree->SetBranchAddress("f_lept1_pt", &f_lept1_pt);
    CRtree->SetBranchAddress("f_lept2_pt", &f_lept2_pt);
    CRtree->SetBranchAddress("f_lept3_pt", &f_lept3_pt);
    CRtree->SetBranchAddress("f_lept4_pt", &f_lept4_pt);
    CRtree->SetBranchAddress("f_lept1_eta", &f_lept1_eta);
    CRtree->SetBranchAddress("f_lept2_eta", &f_lept2_eta);
    CRtree->SetBranchAddress("f_lept3_eta", &f_lept3_eta);
    CRtree->SetBranchAddress("f_lept4_eta", &f_lept4_eta);
    CRtree->SetBranchAddress("f_lept1_pass", &f_lept1_pass);
    CRtree->SetBranchAddress("f_lept2_pass", &f_lept2_pass);
    CRtree->SetBranchAddress("f_lept3_pass", &f_lept3_pass);
    CRtree->SetBranchAddress("f_lept4_pass", &f_lept4_pass);
    CRtree->SetBranchAddress("f_lept1_pdgid", &f_lept1_pdgid);
    CRtree->SetBranchAddress("f_lept2_pdgid", &f_lept2_pdgid);
    CRtree->SetBranchAddress("f_lept3_pdgid", &f_lept3_pdgid);
    CRtree->SetBranchAddress("f_lept4_pdgid", &f_lept4_pdgid);
    CRtree->SetBranchAddress("f_njets_pass", &f_njets_pass);
    CRtree->SetBranchAddress("f_Nbjets", &f_Nbjets);
    CRtree->SetBranchAddress("f_jets_highpt_eta", &f_jets_highpt_eta);
    
    if(ifile == 0){
      CRtree->SetBranchAddress("f_k57nj2", &f_k57nj2);
      CRtree->SetBranchAddress("f_k24nj3", &f_k24nj3);
    }
       
    for(unsigned int iev=0; iev<nCR; ++iev){
      CRtree->GetEntry(iev);

      unsigned int nmu_fail = 0;
      unsigned int nele_fail = 0;
      unsigned int nele = 0;
      unsigned int nmu = 0;
      
      if( ifile == 0 && ((f_njets_pass == 2 && f_k57nj2[0] <= 0.85) || (f_njets_pass >= 3 && f_k24nj3[0] <= 0.73)) ) continue;

  
      //check final state
      if(fabs(f_lept1_pdgid) == 11){
	++nele;
	if(!f_lept1_pass) ++nele_fail;
      }
      if(fabs(f_lept1_pdgid) == 13){
	++nmu;
	if(!f_lept1_pass) ++nmu_fail;
      }//----------------------------
      if(fabs(f_lept2_pdgid) == 11){
	++nele;
	if(!f_lept2_pass) ++nele_fail;
      }
      if(fabs(f_lept2_pdgid) == 13){
	++nmu;
	if(!f_lept2_pass) ++nmu_fail;
      }//---------------------------
      if(fabs(f_lept3_pdgid) == 11){
	++nele;
	if(!f_lept3_pass) ++nele_fail;
      }
      if(fabs(f_lept3_pdgid) == 13){
	++nmu;
	if(!f_lept3_pass) ++nmu_fail;
      }//---------------------------
      if(fabs(f_lept4_pdgid) == 11){
	++nele;
	if(!f_lept4_pass) ++nele_fail;
      }
      if(fabs(f_lept4_pdgid) == 13){
	++nmu;
	if(!f_lept4_pass) ++nmu_fail;
      }
      

      if( (fstate == "4mu" && nmu == 4) || 
	  (fstate == "4e" && nele == 4) || 
	  (fstate == "2e2mu" && nmu == 2 && nele == 2 && nmu_fail > 0) ||
	  (fstate == "2mu2e" && nmu == 2 && nele == 2 && nele_fail > 0) ||
	  fstate == "4l"
      ){
	if( SR == "vbf" &&
	    !( (((f_njets_pass == 2 || f_njets_pass == 3) && f_Nbjets <= 1) || (f_njets_pass > 3 && f_Nbjets == 0)) && f_mass4l >= 118 && f_mass4l <= 130 )
	  ) continue;

	//CRs
	if(files[ifile].Contains("Data")){
	  if(f_2p2f) hm4l_data_2p2f->Fill(f_mass4l, f_weight);
	  if(f_3p1f) hm4l_data_3p1f->Fill(f_mass4l, f_weight);
	  
	  //Computing 3P1F from 2P2F (for plot)
	  if(f_2p2f){
	    float f1=0, f2=0, f3=0, f4=0;
	    if(fabs(f_lept1_pdgid)==11) f1 = ElectronPtEtaMapFR->GetBinContent( ElectronPtEtaMapFR->FindBin(f_lept1_pt,fabs(f_lept1_eta)) );
	    else f1 = MuonPtEtaMapFR->GetBinContent( MuonPtEtaMapFR->FindBin(f_lept1_pt,fabs(f_lept1_eta)) );
	    if(fabs(f_lept2_pdgid)==11) f2 = ElectronPtEtaMapFR->GetBinContent( ElectronPtEtaMapFR->FindBin(f_lept2_pt,fabs(f_lept2_eta)) );
	    else f2 = MuonPtEtaMapFR->GetBinContent( MuonPtEtaMapFR->FindBin(f_lept2_pt,fabs(f_lept2_eta)) );
	    if(fabs(f_lept3_pdgid)==11) f3 = ElectronPtEtaMapFR->GetBinContent( ElectronPtEtaMapFR->FindBin(f_lept3_pt,fabs(f_lept3_eta)) );
	    else f3 = MuonPtEtaMapFR->GetBinContent( MuonPtEtaMapFR->FindBin(f_lept3_pt,fabs(f_lept3_eta)) );
	    if(fabs(f_lept4_pdgid)==11) f4 = ElectronPtEtaMapFR->GetBinContent( ElectronPtEtaMapFR->FindBin(f_lept4_pt,fabs(f_lept4_eta)) );
	    else f4 = MuonPtEtaMapFR->GetBinContent( MuonPtEtaMapFR->FindBin(f_lept4_pt,fabs(f_lept4_eta)) );
	    
	    float weight_2p2f_for_3p1f = 0;
	    if(!f_lept1_pass) weight_2p2f_for_3p1f += f1/(1-f1);
	    if(!f_lept2_pass) weight_2p2f_for_3p1f += f2/(1-f2);
	    if(!f_lept3_pass) weight_2p2f_for_3p1f += f3/(1-f3);
	    if(!f_lept4_pass) weight_2p2f_for_3p1f += f4/(1-f4);
	    hm4l_3p1f_from_2p2f->Fill(f_mass4l, weight_2p2f_for_3p1f);
	    
	    if(!f_lept1_pass && !f_lept2_pass){sum1_2p2f+=f1/(1-f1); sum2_2p2f+=f2/(1-f2);};
	    if(!f_lept1_pass && !f_lept3_pass){sum1_2p2f+=f1/(1-f1); sum2_2p2f+=f3/(1-f3);};
	    if(!f_lept1_pass && !f_lept4_pass){sum1_2p2f+=f1/(1-f1); sum2_2p2f+=f4/(1-f4);};
	    if(!f_lept2_pass && !f_lept3_pass){sum1_2p2f+=f2/(1-f2); sum2_2p2f+=f3/(1-f3);};
	    if(!f_lept2_pass && !f_lept4_pass){sum1_2p2f+=f2/(1-f2); sum2_2p2f+=f4/(1-f4);};
	    if(!f_lept3_pass && !f_lept4_pass){sum1_2p2f+=f3/(1-f3); sum2_2p2f+=f4/(1-f4);};
	  }
	  
	  //Computing 3P1F from 2P2F (for N_bkg_SR - final estimation of Z+X)
	  if(f_2p2f){
	    float f1=0, f2=0, f3=0, f4=0;
	    if(fabs(f_lept1_pdgid)==11) f1 = ElectronPtEtaMapFR->GetBinContent( ElectronPtEtaMapFR->FindBin(f_lept1_pt,fabs(f_lept1_eta)) );
	    else f1 = MuonPtEtaMapFR->GetBinContent( MuonPtEtaMapFR->FindBin(f_lept1_pt,fabs(f_lept1_eta)) );
	    if(fabs(f_lept2_pdgid)==11) f2 = ElectronPtEtaMapFR->GetBinContent( ElectronPtEtaMapFR->FindBin(f_lept2_pt,fabs(f_lept2_eta)) );
	    else f2 = MuonPtEtaMapFR->GetBinContent( MuonPtEtaMapFR->FindBin(f_lept2_pt,fabs(f_lept2_eta)) );
	    if(fabs(f_lept3_pdgid)==11) f3 = ElectronPtEtaMapFR->GetBinContent( ElectronPtEtaMapFR->FindBin(f_lept3_pt,fabs(f_lept3_eta)) );
	    else f3 = MuonPtEtaMapFR->GetBinContent( MuonPtEtaMapFR->FindBin(f_lept3_pt,fabs(f_lept3_eta)) );
	    if(fabs(f_lept4_pdgid)==11) f4 = ElectronPtEtaMapFR->GetBinContent( ElectronPtEtaMapFR->FindBin(f_lept4_pt,fabs(f_lept4_eta)) );
	    else f4 = MuonPtEtaMapFR->GetBinContent( MuonPtEtaMapFR->FindBin(f_lept4_pt,fabs(f_lept4_eta)) );
	    
	    double weight_2p2f_for_3p1f[2] = {-1, -1};
	    double weight_2p2f_for_3p1f_syserr[2] = {-1, -1};
	    if(!f_lept1_pass){ 
	      if(weight_2p2f_for_3p1f[0] == -1){
		weight_2p2f_for_3p1f[0] = f1/(1-f1);
		weight_2p2f_for_3p1f_syserr[0] = FRsys(f1, f_lept1_pdgid);
	      }
	      else if(weight_2p2f_for_3p1f[1] == -1){
		weight_2p2f_for_3p1f[1] = f1/(1-f1);
		weight_2p2f_for_3p1f_syserr[1] = FRsys(f1, f_lept1_pdgid);
	      }
	      else std::cout<<"!!ERROR: More than 2 leptons failing at 2p2f!!"<<std::endl;
	    }
	    if(!f_lept2_pass){ 
	      if(weight_2p2f_for_3p1f[0] == -1){
		weight_2p2f_for_3p1f[0] = f2/(1-f2);
		weight_2p2f_for_3p1f_syserr[0] = FRsys(f2, f_lept2_pdgid);
	      }
	      else if(weight_2p2f_for_3p1f[1] == -1){
		weight_2p2f_for_3p1f[1] = f2/(1-f2);
		weight_2p2f_for_3p1f_syserr[1] = FRsys(f2, f_lept2_pdgid);
	      }
	      else std::cout<<"!!ERROR: More than 2 leptons failing at 2p2f!!"<<std::endl;
	    }
	    if(!f_lept3_pass){ 
	      if(weight_2p2f_for_3p1f[0] == -1){
		weight_2p2f_for_3p1f[0] = f3/(1-f3);
		weight_2p2f_for_3p1f_syserr[0] = FRsys(f3, f_lept3_pdgid);
	      }
	      else if(weight_2p2f_for_3p1f[1] == -1){
		weight_2p2f_for_3p1f[1] = f3/(1-f3);
		weight_2p2f_for_3p1f_syserr[1] = FRsys(f3, f_lept3_pdgid);
	      }
	      else std::cout<<"!!ERROR: More than 2 leptons failing at 2p2f!!"<<std::endl;
	    }
	    if(!f_lept4_pass){ 
	      if(weight_2p2f_for_3p1f[0] == -1){
		weight_2p2f_for_3p1f[0] = f4/(1-f4);
		weight_2p2f_for_3p1f_syserr[0] = FRsys(f4, f_lept4_pdgid);
	      }
	      else if(weight_2p2f_for_3p1f[1] == -1){
		weight_2p2f_for_3p1f[1] = f4/(1-f4);
		weight_2p2f_for_3p1f_syserr[1] = FRsys(f4, f_lept4_pdgid);
	      }
	      else std::cout<<"!!ERROR: More than 2 leptons failing at 2p2f!!"<<std::endl;
	    }
	    
	    hm4l_3p1f_from_2p2f_eq->Fill(f_mass4l, weight_2p2f_for_3p1f[0]*weight_2p2f_for_3p1f[1]);
	    double sigma2p2f_wi = proderr(weight_2p2f_for_3p1f[0], weight_2p2f_for_3p1f_syserr[0], weight_2p2f_for_3p1f[1], weight_2p2f_for_3p1f_syserr[1]);
	    N3p1f_2p2f_FRsys2_sum += (sigma2p2f_wi*sigma2p2f_wi);
	    N3p1f_2p2f_FRsys_sumUp += (weight_2p2f_for_3p1f[0]+weight_2p2f_for_3p1f_syserr[0])*(weight_2p2f_for_3p1f[1]+weight_2p2f_for_3p1f_syserr[1]);
	    N3p1f_2p2f_FRsys_sumDown += (weight_2p2f_for_3p1f[0]-weight_2p2f_for_3p1f_syserr[0])*(weight_2p2f_for_3p1f[1]-weight_2p2f_for_3p1f_syserr[1]);

	    for(unsigned int ishift=0; ishift<nshifts; ++ishift){
	      if(f_njets_pass == 2){
		hk57nj2_2p2f_3p1f[ishift]->Fill(f_k57nj2[ishift], weight_2p2f_for_3p1f[0]*weight_2p2f_for_3p1f[1]);
		hk57nj2_deltajj_2p2f_3p1f[ishift]->Fill(f_k57nj2[ishift], fabs(f_jets_highpt_eta[0]-f_jets_highpt_eta[1]), weight_2p2f_for_3p1f[0]*weight_2p2f_for_3p1f[1]);
	      }
	      if(f_njets_pass >= 3){
		hk24nj3_2p2f_3p1f[ishift]->Fill(f_k24nj3[ishift], weight_2p2f_for_3p1f[0]*weight_2p2f_for_3p1f[1]);
		hk24nj3_deltajj_2p2f_3p1f[ishift]->Fill(f_k24nj3[ishift], fabs(f_jets_highpt_eta[0]-f_jets_highpt_eta[1]), weight_2p2f_for_3p1f[0]*weight_2p2f_for_3p1f[1]);
	      }
	    }
	  }

	  //Weights 3P1F for final estimation of Z+X
	  if(f_3p1f){
	    float f1=0, f2=0, f3=0, f4=0;
	    if(fabs(f_lept1_pdgid)==11) f1 = ElectronPtEtaMapFR->GetBinContent( ElectronPtEtaMapFR->FindBin(f_lept1_pt,fabs(f_lept1_eta)) );
	    else f1 = MuonPtEtaMapFR->GetBinContent( MuonPtEtaMapFR->FindBin(f_lept1_pt,fabs(f_lept1_eta)) );
	    if(fabs(f_lept2_pdgid)==11) f2 = ElectronPtEtaMapFR->GetBinContent( ElectronPtEtaMapFR->FindBin(f_lept2_pt,fabs(f_lept2_eta)) );
	    else f2 = MuonPtEtaMapFR->GetBinContent( MuonPtEtaMapFR->FindBin(f_lept2_pt,fabs(f_lept2_eta)) );
	    if(fabs(f_lept3_pdgid)==11) f3 = ElectronPtEtaMapFR->GetBinContent( ElectronPtEtaMapFR->FindBin(f_lept3_pt,fabs(f_lept3_eta)) );
	    else f3 = MuonPtEtaMapFR->GetBinContent( MuonPtEtaMapFR->FindBin(f_lept3_pt,fabs(f_lept3_eta)) );
	    if(fabs(f_lept4_pdgid)==11) f4 = ElectronPtEtaMapFR->GetBinContent( ElectronPtEtaMapFR->FindBin(f_lept4_pt,fabs(f_lept4_eta)) );
	    else f4 = MuonPtEtaMapFR->GetBinContent( MuonPtEtaMapFR->FindBin(f_lept4_pt,fabs(f_lept4_eta)) );

	    float weight_3p1f = 0;
	    double weight_3p1f_syserr = 0;
	    if(!f_lept1_pass){
	      weight_3p1f = f1/(1-f1);
	      weight_3p1f_syserr = FRsys(f1, f_lept1_pdgid);
	    }
	    if(!f_lept2_pass){
	      weight_3p1f = f2/(1-f2);
	      weight_3p1f_syserr = FRsys(f2, f_lept2_pdgid);
	    }
	    if(!f_lept3_pass){
	      weight_3p1f = f3/(1-f3);
	      weight_3p1f_syserr = FRsys(f3, f_lept3_pdgid);
	    }
	    if(!f_lept4_pass){
	      weight_3p1f = f4/(1-f4);
	      weight_3p1f_syserr = FRsys(f4, f_lept4_pdgid);
	    }
	    
	    hm4l_3p1f->Fill(f_mass4l, weight_3p1f);
	    N3p1f_FRsys2_sum += (weight_3p1f_syserr*weight_3p1f_syserr);
	    N3p1f_FRsys_sumUp += (weight_3p1f+weight_3p1f_syserr);
	    N3p1f_FRsys_sumDown += (weight_3p1f-weight_3p1f_syserr);
	    

	    //NNs
	    for(unsigned int ishift=0; ishift<nshifts; ++ishift){
	      if(f_njets_pass == 2){
		hk57nj2_data_3p1f[ishift]->Fill(f_k57nj2[ishift], weight_3p1f);
		hk57nj2_deltajj_data_3p1f[ishift]->Fill(f_k57nj2[ishift], fabs(f_jets_highpt_eta[0]-f_jets_highpt_eta[1]), weight_3p1f);
	      }
	      if(f_njets_pass >= 3){
		hk24nj3_data_3p1f[ishift]->Fill(f_k24nj3[ishift], weight_3p1f);
		hk24nj3_deltajj_data_3p1f[ishift]->Fill(f_k24nj3[ishift], fabs(f_jets_highpt_eta[0]-f_jets_highpt_eta[1]), weight_3p1f);
	      }
	    }
	  }
	  
	}
	if(files[ifile].Contains("ZZ")){
	  if(f_2p2f) hm4l_zz_2p2f->Fill(f_mass4l, f_weight);
	  if(f_3p1f) hm4l_zz_3p1f->Fill(f_mass4l, f_weight);
	  if(f_3p1f){
	    float f1=0, f2=0, f3=0, f4=0;
	    if(fabs(f_lept1_pdgid)==11) f1 = ElectronPtEtaMapFR->GetBinContent( ElectronPtEtaMapFR->FindBin(f_lept1_pt,fabs(f_lept1_eta)) );
	    else f1 = MuonPtEtaMapFR->GetBinContent( MuonPtEtaMapFR->FindBin(f_lept1_pt,fabs(f_lept1_eta)) );
	    if(fabs(f_lept2_pdgid)==11) f2 = ElectronPtEtaMapFR->GetBinContent( ElectronPtEtaMapFR->FindBin(f_lept2_pt,fabs(f_lept2_eta)) );
	    else f2 = MuonPtEtaMapFR->GetBinContent( MuonPtEtaMapFR->FindBin(f_lept2_pt,fabs(f_lept2_eta)) );
	    if(fabs(f_lept3_pdgid)==11) f3 = ElectronPtEtaMapFR->GetBinContent( ElectronPtEtaMapFR->FindBin(f_lept3_pt,fabs(f_lept3_eta)) );
	    else f3 = MuonPtEtaMapFR->GetBinContent( MuonPtEtaMapFR->FindBin(f_lept3_pt,fabs(f_lept3_eta)) );
	    if(fabs(f_lept4_pdgid)==11) f4 = ElectronPtEtaMapFR->GetBinContent( ElectronPtEtaMapFR->FindBin(f_lept4_pt,fabs(f_lept4_eta)) );
	    else f4 = MuonPtEtaMapFR->GetBinContent( MuonPtEtaMapFR->FindBin(f_lept4_pt,fabs(f_lept4_eta)) );

	    float weight_3p1f = 0;
	    if(!f_lept1_pass) weight_3p1f = f1/(1-f1);
	    if(!f_lept2_pass) weight_3p1f = f2/(1-f2);
	    if(!f_lept3_pass) weight_3p1f = f3/(1-f3);
	    if(!f_lept4_pass) weight_3p1f = f4/(1-f4);
	    hm4l_zz_3p1f_eq->Fill(f_mass4l, f_weight*weight_3p1f);
	  }
	}
	if(files[ifile].Contains("WZ")){
	  if(f_2p2f) hm4l_wz_2p2f->Fill(f_mass4l, f_weight);
	  if(f_3p1f) hm4l_wz_3p1f->Fill(f_mass4l, f_weight);
	}
	if(files[ifile].Contains("TT")){
	  if(f_2p2f) hm4l_ttjets_2p2f->Fill(f_mass4l, f_weight);
	  if(f_3p1f) hm4l_ttjets_3p1f->Fill(f_mass4l, f_weight);
	}
	if(files[ifile].Contains("DY")){
	  if(f_2p2f) hm4l_zjets_2p2f->Fill(f_mass4l, f_weight);
	  if(f_3p1f) hm4l_zjets_3p1f->Fill(f_mass4l, f_weight);
	}
	
       }//finals state if
    }//End of events loop
          
    CRfile->Close();      
  }
  
  

  //Numbers from 2P2F and 3P1F region
  double nwz_2p2f=0, nwz_2p2f_err=0;
  nwz_2p2f = hm4l_wz_2p2f->IntegralAndError(-1,-1,nwz_2p2f_err);
  double nwz_3p1f=0, nwz_3p1f_err=0;
  nwz_3p1f = hm4l_wz_3p1f->IntegralAndError(-1,-1,nwz_3p1f_err);

  double nzz_2p2f=0, nzz_2p2f_err=0;
  nzz_2p2f = hm4l_zz_2p2f->IntegralAndError(-1,-1,nzz_2p2f_err);
  double nzz_3p1f=0, nzz_3p1f_err=0;
  nzz_3p1f = hm4l_zz_3p1f->IntegralAndError(-1,-1,nzz_3p1f_err);

  double ntt_2p2f=0, ntt_2p2f_err=0;
  ntt_2p2f = hm4l_ttjets_2p2f->IntegralAndError(-1,-1,ntt_2p2f_err);
  double ntt_3p1f=0, ntt_3p1f_err=0;
  ntt_3p1f = hm4l_ttjets_3p1f->IntegralAndError(-1,-1,ntt_3p1f_err);
  
  double ndy_2p2f=0, ndy_2p2f_err=0;
  ndy_2p2f = hm4l_zjets_2p2f->IntegralAndError(-1,-1,ndy_2p2f_err);
  double ndy_3p1f=0, ndy_3p1f_err=0;
  ndy_3p1f = hm4l_zjets_3p1f->IntegralAndError(-1,-1,ndy_3p1f_err);

  std::cout<<"\n==== "<<SR<<" = "<<fstate<<" = "<<frmap<<" ===="<<std::endl;
  std::cout<<     " Yields | 2P2F           | 3P1F"<<std::endl;
  std::cout<<Form(" WZ     | %.2f+/-%.2f    | %.2f+/-%.2f",nwz_2p2f,nwz_2p2f_err,nwz_3p1f,nwz_3p1f_err)<<std::endl;
  std::cout<<Form(" ZZ     | %.2f+/-%.2f    | %.2f+/-%.2f",nzz_2p2f,nzz_2p2f_err,nzz_3p1f,nzz_3p1f_err)<<std::endl;
  std::cout<<Form(" TT     | %.2f+/-%.2f  | %.2f+/-%.2f",ntt_2p2f,ntt_2p2f_err,ntt_3p1f,ntt_3p1f_err)<<std::endl;
  std::cout<<Form(" DY     | %.2f+/-%.2f | %.2f+/-%.2f",ndy_2p2f,ndy_2p2f_err,ndy_3p1f,ndy_3p1f_err)<<std::endl;
  std::cout<<"==========================================="<<std::endl;
  
  //Computes Z+X
  double N_ZZ_3p1f_err=0;
  float N_ZZ_3p1f = hm4l_zz_3p1f->IntegralAndError(-1,-1,N_ZZ_3p1f_err);
  double N_ZZ_3p1f_eq_err=0;
  float N_ZZ_3p1f_eq = hm4l_zz_3p1f_eq->IntegralAndError(-1,-1,N_ZZ_3p1f_eq_err);
  double N_bkg_3p1f_err=0;
  float N_bkg_3p1f = hm4l_3p1f_from_2p2f->IntegralAndError(-1,-1,N_bkg_3p1f_err);
  float N_3p1f_entries = hm4l_3p1f->GetEntries();
  double N_3p1f_err=0;
  float N_3p1f = hm4l_3p1f->IntegralAndError(-1,-1,N_3p1f_err);
  double N_3p1f_from_2p2f_eq_err=0;
  float N_3p1f_from_2p2f_eq = hm4l_3p1f_from_2p2f_eq->IntegralAndError(-1,-1,N_3p1f_from_2p2f_eq_err);

  std::cout<<"=========== Final estimation for Z+X ==============="<<std::endl;
  //std::cout<<Form("N_3p1f_entries:      %.1f +/- %.1f",N_3p1f_entries,sqrt(N_3p1f_entries))<<std::endl;
  //std::cout<<Form("N_3p1f:              %.1f +/- %.1f",N_3p1f,N_3p1f_err)<<std::endl;
  //std::cout<<Form("N_bkg_3p1f:          %.1f +/- %.1f",N_bkg_3p1f,N_bkg_3p1f_err)<<std::endl;
  //std::cout<<Form("N_ZZ_3p1f:           %.1f +/- %.1f",N_ZZ_3p1f,N_ZZ_3p1f_err)<<std::endl;
  //std::cout<<Form("N_ZZ_3p1f_eq:        %.1f +/- %.1f",N_ZZ_3p1f_eq,N_ZZ_3p1f_eq_err)<<std::endl;
  //std::cout<<Form("N_3p1f_from_2p2f_eq: %.1f +/- %.1f (This is Z+X if N3p1f-Nbkg3p1f-Nzz3p1f is negative - said to me)",N_3p1f_from_2p2f_eq,N_3p1f_from_2p2f_eq_err)<<std::endl;

  //std::cout<< "------- Using 1st formula --------"<<std::endl;
  //float nzx1 = (N_3p1f - N_ZZ_3p1f_eq - 2*N_3p1f_from_2p2f_eq) + N_3p1f_from_2p2f_eq;
  //error propagation
  //float nzx1_err = sqrt(pow(N_3p1f_err,2)+pow(N_ZZ_3p1f_eq,2)+pow(2*N_3p1f_from_2p2f_eq_err,2)+pow(N_3p1f_from_2p2f_eq_err,2));
  //std::cout<<Form("%.2f +/- %.2f",nzx1,nzx1_err)<<std::endl;
  
  //std::cout<< "------- Using 2nd formula --------"<<std::endl;
  //float nzx = (1-N_ZZ_3p1f/N_3p1f_entries)*N_3p1f - N_3p1f_from_2p2f_eq;
  //error propagation
  //float R1 = N_ZZ_3p1f/N_3p1f_entries;
  //float R1e = R1*sqrt(pow(N_ZZ_3p1f_err/N_ZZ_3p1f,2)+pow(sqrt(N_3p1f_entries)/N_3p1f_entries,2));
  //float C2e = (1-R1)*N_3p1f*sqrt(pow(R1e/(1-R1),2)+pow(N_3p1f_err/N_3p1f,2));
  //float nzx_err = sqrt(pow(C2e,2)+pow(N_3p1f_from_2p2f_eq_err,2));
  //std::cout<<Form("%.2f +/- %.2f",nzx,nzx_err)<<std::endl;
  
  //Operating through histograms Z+X shape and yields
  TH1D *h3p1f_fr1 = (TH1D*)hm4l_3p1f->Clone();
  h3p1f_fr1->Add(hm4l_zz_3p1f_eq, -1);
  h3p1f_fr1->Add(hm4l_3p1f_from_2p2f_eq, -1);
  double hInt1Err=0;
  float hInt1 = h3p1f_fr1->IntegralAndError(-1,-1,hInt1Err);
  std::cout<<Form("Z+X (1st formula, histogram-based) = %.2f +/- %.2f",hInt1,hInt1Err)<<std::endl;  
  
  TH1D *h3p1f_fr2 = (TH1D*)hm4l_3p1f->Clone();
  float sc_factor = (1-N_ZZ_3p1f/N_3p1f_entries);
  (*h3p1f_fr2) = (*h3p1f_fr2)*sc_factor;
  h3p1f_fr2->Add(hm4l_3p1f_from_2p2f_eq,-1);//This = fi/(1-fi)*fj/(1-fj)
  double hInt2Err=0;
  float hInt2 = h3p1f_fr2->IntegralAndError(-1,-1,hInt2Err);
  std::cout<<Form("Z+X (2nd formula, histogram-based) = %.2f +/- %.2f",hInt2,hInt2Err)<<std::endl;  
  double n3p1f_2p2f_err=0;
  float n3p1f_2p2f = hm4l_3p1f_from_2p2f_eq->IntegralAndError(-1,-1,n3p1f_2p2f_err);
  if(fstate == "4mu" || fstate == "2e2mu") std::cout<<Form("N3p1f from 2p2f = %.2f +/- %.2f",n3p1f_2p2f,n3p1f_2p2f_err)<<std::endl;
  std::cout<<"=================================================="<<std::endl;
  std::cout<<Form("N3p1f_FRsys2_sum = %.3f, sqrt() = %.3f, N3p1f_2p2f_FRsys2_sum = %.3f, sqrt() = %.3f",N3p1f_FRsys2_sum,sqrt(N3p1f_FRsys2_sum),N3p1f_2p2f_FRsys2_sum,sqrt(N3p1f_2p2f_FRsys2_sum))<<std::endl;
  std::cout<<Form("Z+X sys (loose mu)   = +/-%.2f",sqrt(N3p1f_2p2f_FRsys2_sum))<<std::endl;
  double zxerr = sqrt( pow(proderr(sc_factor, 0, N_3p1f, sqrt(N3p1f_FRsys2_sum)),2) + N3p1f_2p2f_FRsys2_sum );
  std::cout<<Form("Z+X sys (noloose mu) = +/-%.2f",zxerr)<<std::endl;
  std::cout<<Form("N3p1f_FRsys_sumUp = %.3f, N3p1f_2p2f_FRsys_sumUp = %.3f",N3p1f_FRsys_sumUp,N3p1f_2p2f_FRsys_sumUp)<<std::endl;
  std::cout<<Form("N3p1f_FRsys_sumDown = %.3f, N3p1f_2p2f_FRsys_sumDown = %.3f",N3p1f_FRsys_sumDown,N3p1f_2p2f_FRsys_sumDown)<<std::endl;
  std::cout<<Form("Z+X sys (loose mu up)   = +/-%.2f",N3p1f_2p2f_FRsys_sumUp)<<std::endl;
  std::cout<<Form("Z+X sys (loose mu down)   = +/-%.2f",N3p1f_2p2f_FRsys_sumDown)<<std::endl;
  double zxerr2 = sc_factor*N3p1f_FRsys_sumUp - N3p1f_2p2f_FRsys_sumUp;
  std::cout<<Form("Z+X sys (noloose mu up up) = +/-%.2f",zxerr2)<<std::endl;
  double zxerr3 = sc_factor*N3p1f_FRsys_sumUp - N3p1f_2p2f_FRsys_sumDown;
  std::cout<<Form("Z+X sys (noloose mu up down) = +/-%.2f",zxerr3)<<std::endl;
  double zxerr4 = sc_factor*N3p1f_FRsys_sumDown - N3p1f_2p2f_FRsys_sumDown;
  std::cout<<Form("Z+X sys (noloose mu down down) = +/-%.2f",zxerr4)<<std::endl;
  double zxerr5 = sc_factor*N3p1f_FRsys_sumDown - N3p1f_2p2f_FRsys_sumUp;
  std::cout<<Form("Z+X sys (noloose mu down up) = +/-%.2f",zxerr5)<<std::endl;
  
  TFile *zxstatistics = new TFile("zxstatisticsOS/"+SR+"_"+fstate+"_"+frmap+".root","recreate");
  TH1D *form1 = new TH1D("form1","",1,0,1);
  TH1D *form2 = new TH1D("form2","",1,0,1);
  if(fstate != "4mu" && fstate != "2e2mu"){
    h3p1f_fr2->SetName("m4lZX");
    h3p1f_fr2->Write();
    form1->SetBinContent(1,hInt1);
    form1->SetBinError(1,hInt1Err);
    form1->Write();
    form2->SetBinContent(1,hInt2);
    form2->SetBinError(1,hInt2Err);
    form2->Write();
  }else{
    hm4l_3p1f_from_2p2f_eq->SetName("m4lZX");
    hm4l_3p1f_from_2p2f_eq->Write();
    form1->SetBinContent(1,n3p1f_2p2f);
    form1->SetBinError(1,n3p1f_2p2f_err);
    form1->Write();
    form2->SetBinContent(1,n3p1f_2p2f);
    form2->SetBinError(1,n3p1f_2p2f_err);
    form2->Write();
  }
  zxstatistics->Close();
  
  
  TString ffstate = fstate;
  if(fstate == "2e2mu" || fstate == "2mu2e") ffstate = "2e2mu";
    
  TFile *k57nj2_shape = new TFile("DistributionsUncertainties_k57nj2_proc_zx_"+fstate+"_"+SR+"_"+frmap+".root","recreate");
  h3p1f_fr2->Write();
  k57nj2_shape->mkdir("ch"+ffstate);
  k57nj2_shape->cd("ch"+ffstate);
  for(unsigned int ishift=0; ishift<nshifts; ++ishift){
    if(fstate != "4mu" && fstate != "2e2mu"){
      (*hk57nj2_data_3p1f[ishift]) = (*hk57nj2_data_3p1f[ishift])*sc_factor;
      hk57nj2_data_3p1f[ishift]->Add(hk57nj2_2p2f_3p1f[ishift], -1);
      hk57nj2_data_3p1f[ishift]->SetName("zjets"+shift_name[ishift]);
      hk57nj2_data_3p1f[ishift]->Write();
    }else{
      hk57nj2_2p2f_3p1f[ishift]->SetName("zjets"+shift_name[ishift]);
      hk57nj2_2p2f_3p1f[ishift]->Write();      
    }
  }
  k57nj2_shape->Close();
    
  TFile *k24nj3_shape = new TFile("DistributionsUncertainties_k24nj3_proc_zx_"+fstate+"_"+SR+"_"+frmap+".root","recreate");
  h3p1f_fr2->Write();
  k24nj3_shape->mkdir("ch"+ffstate);
  k24nj3_shape->cd("ch"+ffstate);
  for(unsigned int ishift=0; ishift<nshifts; ++ishift){
    if(fstate != "4mu" && fstate != "2e2mu"){
      (*hk24nj3_data_3p1f[ishift]) = (*hk24nj3_data_3p1f[ishift])*sc_factor;
      hk24nj3_data_3p1f[ishift]->Add(hk24nj3_2p2f_3p1f[ishift], -1);
      hk24nj3_data_3p1f[ishift]->SetName("zjets"+shift_name[ishift]);
      hk24nj3_data_3p1f[ishift]->Write();
    }else{
      hk24nj3_2p2f_3p1f[ishift]->SetName("zjets"+shift_name[ishift]);
      hk24nj3_2p2f_3p1f[ishift]->Write();      
    }
  }
  k24nj3_shape->Close();

  TFile *k57nj2_deltajj_shape = new TFile("DistributionsUncertainties_k57nj2_deltajj_proc_zx_"+fstate+"_"+SR+"_"+frmap+".root","recreate");
  h3p1f_fr2->Write();
  k57nj2_deltajj_shape->mkdir("ch"+ffstate);
  k57nj2_deltajj_shape->cd("ch"+ffstate);
  for(unsigned int ishift=0; ishift<nshifts; ++ishift){
    if(fstate != "4mu" && fstate != "2e2mu"){
      (*hk57nj2_deltajj_data_3p1f[ishift]) = (*hk57nj2_deltajj_data_3p1f[ishift])*sc_factor;
      hk57nj2_deltajj_data_3p1f[ishift]->Add(hk57nj2_deltajj_2p2f_3p1f[ishift], -1);
      hk57nj2_deltajj_data_3p1f[ishift]->SetName("zjets2D"+shift_name[ishift]);
      hk57nj2_deltajj_data_3p1f[ishift]->Write();
    }else{
      hk57nj2_deltajj_2p2f_3p1f[ishift]->SetName("zjets2D"+shift_name[ishift]);
      hk57nj2_deltajj_2p2f_3p1f[ishift]->Write();      
    }
  }
  k57nj2_deltajj_shape->Close();
    
  TFile *k24nj3_deltajj_shape = new TFile("DistributionsUncertainties_k24nj3_deltajj_proc_zx_"+fstate+"_"+SR+"_"+frmap+".root","recreate");
  h3p1f_fr2->Write();
  k24nj3_deltajj_shape->mkdir("ch"+ffstate);
  k24nj3_deltajj_shape->cd("ch"+ffstate);
  for(unsigned int ishift=0; ishift<nshifts; ++ishift){
    if(fstate != "4mu" && fstate != "2e2mu"){
      (*hk24nj3_deltajj_data_3p1f[ishift]) = (*hk24nj3_deltajj_data_3p1f[ishift])*sc_factor;
      hk24nj3_deltajj_data_3p1f[ishift]->Add(hk24nj3_deltajj_2p2f_3p1f[ishift], -1);
      hk24nj3_deltajj_data_3p1f[ishift]->SetName("zjets2D"+shift_name[ishift]);
      hk24nj3_deltajj_data_3p1f[ishift]->Write();
    }else{
      hk24nj3_deltajj_2p2f_3p1f[ishift]->SetName("zjets2D"+shift_name[ishift]);
      hk24nj3_deltajj_2p2f_3p1f[ishift]->Write();      
    }
  }
  k24nj3_deltajj_shape->Close();
  

  //=================================== plots ================================
  //Stack histograms for 2p2f
  hm4l_zz_2p2f->Rebin(20);
  hm4l_wz_2p2f->Rebin(20);
  hm4l_ttjets_2p2f->Rebin(20);
  hm4l_zjets_2p2f->Rebin(20);
  hm4l_data_2p2f->Rebin(20);
  
  hm4lStack_2p2f->Add(hm4l_zz_2p2f);
  hm4lStack_2p2f->Add(hm4l_wz_2p2f);
  hm4lStack_2p2f->Add(hm4l_ttjets_2p2f);
  hm4lStack_2p2f->Add(hm4l_zjets_2p2f);
    
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cv","",10,10,700,700);
  hm4l_data_2p2f->SetLineColor(kBlack);
  hm4l_data_2p2f->SetMarkerColor(kBlack);
  hm4l_data_2p2f->SetMarkerSize(0.8);
  hm4l_data_2p2f->SetLineWidth(1);
  hm4l_data_2p2f->GetXaxis()->SetTitle("m_{4l}(GeV)");
  hm4l_data_2p2f->GetYaxis()->SetTitle(Form("Events/%0.fGeV",hm4l_data_2p2f->GetBinWidth(1)));
  hm4l_data_2p2f->Draw("pe");
  hm4lStack_2p2f->Draw("hist,same");
  hm4l_data_2p2f->Draw("pe,same");
  hm4l_data_2p2f->GetXaxis()->SetRangeUser(64,800);
  gPad->RedrawAxis();
    
  TLegend *leg = new TLegend(0.45,0.6,0.85,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetHeader("2P2F - "+final_state);
  leg->AddEntry(hm4l_data_2p2f,"Data","pl");
  leg->AddEntry(hm4l_zjets_2p2f,"Z+jets","f");
  leg->AddEntry(hm4l_ttjets_2p2f,"t#bar{t}+jets","f");
  leg->AddEntry(hm4l_wz_2p2f,"WZ","f");
  leg->AddEntry(hm4l_zz_2p2f,"ZZ, Z#gamma*","f");
  leg->Draw();
    
  cv->Update();
  if(frmap == "nominal") cv->Print("m4l_2p2f_"+fstate+"_"+SR+".png");
  cv->Close();
  
  //Stack histograms for 3p1f
  hm4l_zz_3p1f->Rebin(20);
  hm4l_wz_3p1f->Rebin(20);
  hm4l_ttjets_3p1f->Rebin(20);
  hm4l_zjets_3p1f->Rebin(20);
  hm4l_data_3p1f->Rebin(20);
  hm4l_3p1f_from_2p2f->Rebin(20);
  
  hm4lStack_3p1f->Add(hm4l_zz_3p1f);
  hm4lStack_3p1f->Add(hm4l_wz_3p1f);
  hm4lStack_3p1f->Add(hm4l_ttjets_3p1f);
  hm4lStack_3p1f->Add(hm4l_zjets_3p1f);

  gStyle->SetOptStat(0);
  TCanvas *cv2 = new TCanvas("cv2","",50,50,700,700);
  hm4l_data_3p1f->SetLineColor(kBlack);
  hm4l_data_3p1f->SetMarkerColor(kBlack);
  hm4l_data_3p1f->SetMarkerSize(0.8);
  hm4l_data_3p1f->SetLineWidth(1);
  hm4l_data_3p1f->GetXaxis()->SetTitle("m_{4l}(GeV)");
  hm4l_data_3p1f->GetYaxis()->SetTitle(Form("Events/%0.fGeV",hm4l_data_3p1f->GetBinWidth(1)));
  hm4l_data_3p1f->Draw("pe");
  hm4lStack_3p1f->Draw("hist,same");
  hm4l_3p1f_from_2p2f->Draw("hist,same");
  hm4l_data_3p1f->Draw("pe,same");
  hm4l_data_3p1f->GetXaxis()->SetRangeUser(64,800);
  gPad->RedrawAxis();
    
  TLegend *leg2 = new TLegend(0.45,0.6,0.85,0.9);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetHeader("3P1F - "+final_state);
  leg2->AddEntry(hm4l_data_3p1f,"Data","pl");
  leg2->AddEntry(hm4l_zjets_3p1f,"Z+jets","f");
  leg2->AddEntry(hm4l_ttjets_3p1f,"t#bar{t}+jets","f");        
  leg2->AddEntry(hm4l_wz_3p1f,"WZ","f");        
  leg2->AddEntry(hm4l_zz_3p1f,"ZZ, Z#gamma*","f");
  leg2->AddEntry(hm4l_3p1f_from_2p2f,"2P2F extr.","f");        
  leg2->Draw();

  cv2->Update();
  if(frmap == "nominal") cv2->Print("m4l_3p1f_"+fstate+"_"+SR+".png");
  cv2->Close();
  
  
  TCanvas *fcv = new TCanvas("fcv","",100,100,700,700);
  h3p1f_fr2->SetTitle("Z+X");
  h3p1f_fr2->GetXaxis()->SetTitle("m_{4l}(GeV)");
  h3p1f_fr2->SetLineColor(kBlack);
  h3p1f_fr2->SetMarkerColor(kBlack);
  h3p1f_fr2->SetMarkerSize(0.8);
  h3p1f_fr2->SetLineWidth(1);
  h3p1f_fr2->Draw("pe");
  if(SR == "smhiggs"){
    h3p1f_fr2->Rebin(20);
    h3p1f_fr2->GetXaxis()->SetRangeUser(64,800);
  }
  if(SR == "vbf"){
    h3p1f_fr2->Rebin(4);
    h3p1f_fr2->GetXaxis()->SetRangeUser(109,141);
  }
  h3p1f_fr2->GetYaxis()->SetTitle(Form("Events/%0.fGeV",h3p1f_fr2->GetBinWidth(1)));
  h3p1f_fr2->Draw("pe,same");
  
  TLegend *leg3 = new TLegend(0.3,0.8,0.9,0.9);
  leg3->SetFillColor(0);
  leg3->SetFillStyle(0);
  leg3->SetHeader("N^{Z+X}_{SR} = (1 - #frac{N^{ZZ}_{3P1F}}{N^{Data}_{3P1F}}) #sum^{N^{Data}_{3P1F}}_{i}#frac{f^{1}_{i}}{(1-f^{1}_{i})} - #sum^{N^{Data}_{2P2F}}_{j}#frac{f^{1}_{j}}{(1-f^{1}_{j})} #frac{f^{2}_{j}}{(1-f^{2}_{j})}");
  leg3->Draw();
  
  fcv->Update();  
  if(frmap == "nominal") fcv->Print("m4l_zx_"+fstate+"_"+SR+".png");
  fcv->Close();
    
  //END
} 



void plotCRs(void){
  std::vector<TString> sr = {"vbf"};
  std::vector<TString> fs = {"4mu","4e","2e2mu", "2mu2e"};
  std::vector<TString> fr = {"nominal"};//,"average","reweight"};
  
  for(unsigned int isr=0; isr<(unsigned int)sr.size(); ++isr)
    for(unsigned int ifs=0; ifs<(unsigned int)fs.size(); ++ifs)
      for(unsigned int ifr=0; ifr<(unsigned int)fr.size(); ++ifr)
	analyzeCRs(sr.at(isr), fs.at(ifs), fr.at(ifr));
      
}
