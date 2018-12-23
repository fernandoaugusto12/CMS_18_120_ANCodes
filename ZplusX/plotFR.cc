#define ele_sys 0.249
#define mu_sys 0.133

//Code to extract FR and plot it
void plotFR(TString datafile_name, TString mcfile_name){
  gROOT->SetBatch();
  TH1::SetDefaultSumw2();
  
  
  Float_t f_weight=-999;
  Float_t f_zee_mu_pt_loose=-999, f_zee_e_pt_loose=-999, f_zmumu_mu_pt_loose=-999, f_zmumu_e_pt_loose=-999;
  Float_t f_zee_mu_eta_loose=-999, f_zee_e_eta_loose=-999, f_zmumu_mu_eta_loose=-999, f_zmumu_e_eta_loose=-999;
  Float_t f_zee_mu_charge_loose=-999, f_zee_e_charge_loose=-999, f_zmumu_mu_charge_loose=-999, f_zmumu_e_charge_loose=-999;
  Float_t f_zee_mu_pt_tight=-999, f_zee_e_pt_tight=-999, f_zmumu_mu_pt_tight=-999, f_zmumu_e_pt_tight=-999;
  Float_t f_zee_mu_eta_tight=-999, f_zee_e_eta_tight=-999, f_zmumu_mu_eta_tight=-999, f_zmumu_e_eta_tight=-999;
  Float_t f_zee_mu_charge_tight=-999, f_zee_e_charge_tight=-999, f_zmumu_mu_charge_tight=-999, f_zmumu_e_charge_tight=-999;  

  
  TFile *datafile = TFile::Open(datafile_name);
  TTree *dataFR = (TTree*)datafile->Get("FakeRate");
  unsigned int ndata = dataFR->GetEntries();
  dataFR->SetBranchAddress("f_zee_mu_pt_loose", &f_zee_mu_pt_loose);
  dataFR->SetBranchAddress("f_zee_mu_eta_loose", &f_zee_mu_eta_loose);
  dataFR->SetBranchAddress("f_zee_e_pt_loose", &f_zee_e_pt_loose);
  dataFR->SetBranchAddress("f_zee_e_eta_loose", &f_zee_e_eta_loose);
  dataFR->SetBranchAddress("f_zmumu_mu_pt_loose", &f_zmumu_mu_pt_loose);
  dataFR->SetBranchAddress("f_zmumu_mu_eta_loose", &f_zmumu_mu_eta_loose);
  dataFR->SetBranchAddress("f_zmumu_e_pt_loose", &f_zmumu_e_pt_loose);
  dataFR->SetBranchAddress("f_zmumu_e_eta_loose", &f_zmumu_e_eta_loose);
  dataFR->SetBranchAddress("f_zee_mu_pt_tight", &f_zee_mu_pt_tight);
  dataFR->SetBranchAddress("f_zee_mu_eta_tight", &f_zee_mu_eta_tight);
  dataFR->SetBranchAddress("f_zee_e_pt_tight", &f_zee_e_pt_tight);
  dataFR->SetBranchAddress("f_zee_e_eta_tight", &f_zee_e_eta_tight);
  dataFR->SetBranchAddress("f_zmumu_mu_pt_tight", &f_zmumu_mu_pt_tight);
  dataFR->SetBranchAddress("f_zmumu_mu_eta_tight", &f_zmumu_mu_eta_tight);
  dataFR->SetBranchAddress("f_zmumu_e_pt_tight", &f_zmumu_e_pt_tight);
  dataFR->SetBranchAddress("f_zmumu_e_eta_tight", &f_zmumu_e_eta_tight);
  dataFR->SetBranchAddress("f_zee_mu_charge_loose", &f_zee_mu_charge_loose);
  dataFR->SetBranchAddress("f_zee_e_charge_loose", &f_zee_e_charge_loose);
  dataFR->SetBranchAddress("f_zmumu_mu_charge_loose", &f_zmumu_mu_charge_loose);
  dataFR->SetBranchAddress("f_zmumu_e_charge_loose", &f_zmumu_e_charge_loose);
  dataFR->SetBranchAddress("f_zee_mu_charge_tight", &f_zee_mu_charge_tight);
  dataFR->SetBranchAddress("f_zee_e_charge_tight", &f_zee_e_charge_tight);
  dataFR->SetBranchAddress("f_zmumu_mu_charge_tight", &f_zmumu_mu_charge_tight);
  dataFR->SetBranchAddress("f_zmumu_e_charge_tight", &f_zmumu_e_charge_tight);

  TFile *mcfile = TFile::Open(mcfile_name);
  TTree *mcFR = (TTree*)mcfile->Get("FakeRate");
  unsigned int nmc = mcFR->GetEntries();
  mcFR->SetBranchAddress("f_weight", &f_weight);
  mcFR->SetBranchAddress("f_zee_mu_pt_loose", &f_zee_mu_pt_loose);
  mcFR->SetBranchAddress("f_zee_mu_eta_loose", &f_zee_mu_eta_loose);
  mcFR->SetBranchAddress("f_zee_e_pt_loose", &f_zee_e_pt_loose);
  mcFR->SetBranchAddress("f_zee_e_eta_loose", &f_zee_e_eta_loose);
  mcFR->SetBranchAddress("f_zmumu_mu_pt_loose", &f_zmumu_mu_pt_loose);
  mcFR->SetBranchAddress("f_zmumu_mu_eta_loose", &f_zmumu_mu_eta_loose);
  mcFR->SetBranchAddress("f_zmumu_e_pt_loose", &f_zmumu_e_pt_loose);
  mcFR->SetBranchAddress("f_zmumu_e_eta_loose", &f_zmumu_e_eta_loose);
  mcFR->SetBranchAddress("f_zee_mu_pt_tight", &f_zee_mu_pt_tight);
  mcFR->SetBranchAddress("f_zee_mu_eta_tight", &f_zee_mu_eta_tight);
  mcFR->SetBranchAddress("f_zee_e_pt_tight", &f_zee_e_pt_tight);
  mcFR->SetBranchAddress("f_zee_e_eta_tight", &f_zee_e_eta_tight);
  mcFR->SetBranchAddress("f_zmumu_mu_pt_tight", &f_zmumu_mu_pt_tight);
  mcFR->SetBranchAddress("f_zmumu_mu_eta_tight", &f_zmumu_mu_eta_tight);
  mcFR->SetBranchAddress("f_zmumu_e_pt_tight", &f_zmumu_e_pt_tight);
  mcFR->SetBranchAddress("f_zmumu_e_eta_tight", &f_zmumu_e_eta_tight);
  mcFR->SetBranchAddress("f_zee_mu_charge_loose", &f_zee_mu_charge_loose);
  mcFR->SetBranchAddress("f_zee_e_charge_loose", &f_zee_e_charge_loose);
  mcFR->SetBranchAddress("f_zmumu_mu_charge_loose", &f_zmumu_mu_charge_loose);
  mcFR->SetBranchAddress("f_zmumu_e_charge_loose", &f_zmumu_e_charge_loose);
  mcFR->SetBranchAddress("f_zee_mu_charge_tight", &f_zee_mu_charge_tight);
  mcFR->SetBranchAddress("f_zee_e_charge_tight", &f_zee_e_charge_tight);
  mcFR->SetBranchAddress("f_zmumu_mu_charge_tight", &f_zmumu_mu_charge_tight);
  mcFR->SetBranchAddress("f_zmumu_e_charge_tight", &f_zmumu_e_charge_tight);
  
  //pt data 
  const int Nbins = 6;
  Double_t Bins[Nbins+1] = {0,5,10,20,30,50,1000};  
  TH1D *h_mu_pt_loose_barrel = new TH1D("h_mu_pt_loose_barrel","",Nbins,Bins);
  TH1D *h_mu_pt_tight_barrel = new TH1D("h_mu_pt_tight_barrel","",Nbins,Bins);
  TH1D *h_mu_pt_loose_endcap = new TH1D("h_mu_pt_loose_endcap","",Nbins,Bins);
  TH1D *h_mu_pt_tight_endcap = new TH1D("h_mu_pt_tight_endcap","",Nbins,Bins);

  TH1D *h_e_pt_loose_barrel = new TH1D("h_e_pt_loose_barrel","",Nbins,Bins);
  TH1D *h_e_pt_tight_barrel = new TH1D("h_e_pt_tight_barrel","",Nbins,Bins);
  TH1D *h_e_pt_loose_endcap = new TH1D("h_e_pt_loose_endcap","",Nbins,Bins);
  TH1D *h_e_pt_tight_endcap = new TH1D("h_e_pt_tight_endcap","",Nbins,Bins);

  TH1D *h_mu_pt_loose_barrel2 = new TH1D("h_mu_pt_loose_barrel2","",Nbins,Bins);
  TH1D *h_mu_pt_tight_barrel2 = new TH1D("h_mu_pt_tight_barrel2","",Nbins,Bins);
  TH1D *h_mu_pt_loose_endcap2 = new TH1D("h_mu_pt_loose_endcap2","",Nbins,Bins);
  TH1D *h_mu_pt_tight_endcap2 = new TH1D("h_mu_pt_tight_endcap2","",Nbins,Bins);

  TH1D *h_e_pt_loose_barrel2 = new TH1D("h_e_pt_loose_barrel2","",Nbins,Bins);
  TH1D *h_e_pt_tight_barrel2 = new TH1D("h_e_pt_tight_barrel2","",Nbins,Bins);
  TH1D *h_e_pt_loose_endcap2 = new TH1D("h_e_pt_loose_endcap2","",Nbins,Bins);
  TH1D *h_e_pt_tight_endcap2 = new TH1D("h_e_pt_tight_endcap2","",Nbins,Bins);
  
  //eta
  const int Nbins2 = 8;
  Double_t Bins2[Nbins2+1] = {0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.5};  
  TH1D *h_mu_eta_loose_barrel = new TH1D("h_mu_eta_loose_barrel","",Nbins2,Bins2);
  TH1D *h_mu_eta_tight_barrel = new TH1D("h_mu_eta_tight_barrel","",Nbins2,Bins2);
  TH1D *h_mu_eta_loose_endcap = new TH1D("h_mu_eta_loose_endcap","",Nbins2,Bins2);
  TH1D *h_mu_eta_tight_endcap = new TH1D("h_mu_eta_tight_endcap","",Nbins2,Bins2);

  TH1D *h_e_eta_loose_barrel = new TH1D("h_e_eta_loose_barrel","",Nbins2,Bins2);
  TH1D *h_e_eta_tight_barrel = new TH1D("h_e_eta_tight_barrel","",Nbins2,Bins2);
  TH1D *h_e_eta_loose_endcap = new TH1D("h_e_eta_loose_endcap","",Nbins2,Bins2);
  TH1D *h_e_eta_tight_endcap = new TH1D("h_e_eta_tight_endcap","",Nbins2,Bins2);

  TH1D *h_mu_eta_loose_barrel2 = new TH1D("h_mu_eta_loose_barrel2","",Nbins2,Bins2);
  TH1D *h_mu_eta_tight_barrel2 = new TH1D("h_mu_eta_tight_barrel2","",Nbins2,Bins2);
  TH1D *h_mu_eta_loose_endcap2 = new TH1D("h_mu_eta_loose_endcap2","",Nbins2,Bins2);
  TH1D *h_mu_eta_tight_endcap2 = new TH1D("h_mu_eta_tight_endcap2","",Nbins2,Bins2);

  TH1D *h_e_eta_loose_barrel2 = new TH1D("h_e_eta_loose_barrel2","",Nbins2,Bins2);
  TH1D *h_e_eta_tight_barrel2 = new TH1D("h_e_eta_tight_barrel2","",Nbins2,Bins2);
  TH1D *h_e_eta_loose_endcap2 = new TH1D("h_e_eta_loose_endcap2","",Nbins2,Bins2);
  TH1D *h_e_eta_tight_endcap2 = new TH1D("h_e_eta_tight_endcap2","",Nbins2,Bins2);
  
  //2D  
  TH2D *h_e_pt_eta_loose = new TH2D("h_e_pt_eta_loose","",Nbins,Bins,Nbins2,Bins2);
  TH2D *h_e_pt_eta_tight = new TH2D("h_e_pt_eta_tight","",Nbins,Bins,Nbins2,Bins2);

  TH2D *h_mu_pt_eta_loose = new TH2D("h_mu_pt_eta_loose","",Nbins,Bins,Nbins2,Bins2);
  TH2D *h_mu_pt_eta_tight = new TH2D("h_mu_pt_eta_tight","",Nbins,Bins,Nbins2,Bins2);

  TH2D *h_e_pt_eta_loose2 = new TH2D("h_e_pt_eta_loose2","",Nbins,Bins,Nbins2,Bins2);
  TH2D *h_e_pt_eta_tight2 = new TH2D("h_e_pt_eta_tight2","",Nbins,Bins,Nbins2,Bins2);

  TH2D *h_mu_pt_eta_loose2 = new TH2D("h_mu_pt_eta_loose2","",Nbins,Bins,Nbins2,Bins2);
  TH2D *h_mu_pt_eta_tight2 = new TH2D("h_mu_pt_eta_tight2","",Nbins,Bins,Nbins2,Bins2);
  
  float e_eta = 1.479;
  float mu_eta = 1.2;
  for(unsigned int iev=0; iev<ndata; ++iev){
    dataFR->GetEntry(iev);
    
    //pt
    if(fabs(f_zee_mu_eta_loose) <= mu_eta) h_mu_pt_loose_barrel->Fill(f_zee_mu_pt_loose);
    if(fabs(f_zee_mu_eta_loose) > mu_eta) h_mu_pt_loose_endcap->Fill(f_zee_mu_pt_loose);
    if(fabs(f_zee_mu_eta_tight) <= mu_eta) h_mu_pt_tight_barrel->Fill(f_zee_mu_pt_tight);
    if(fabs(f_zee_mu_eta_tight) > mu_eta) h_mu_pt_tight_endcap->Fill(f_zee_mu_pt_tight);
    
    if(fabs(f_zee_e_eta_loose) <= e_eta) h_e_pt_loose_barrel->Fill(f_zee_e_pt_loose);
    if(fabs(f_zee_e_eta_loose) > e_eta) h_e_pt_loose_endcap->Fill(f_zee_e_pt_loose);
    if(fabs(f_zee_e_eta_tight) <= e_eta) h_e_pt_tight_barrel->Fill(f_zee_e_pt_tight);
    if(fabs(f_zee_e_eta_tight) > e_eta) h_e_pt_tight_endcap->Fill(f_zee_e_pt_tight);
    
    if(fabs(f_zmumu_mu_eta_loose) <= mu_eta) h_mu_pt_loose_barrel->Fill(f_zmumu_mu_pt_loose);
    if(fabs(f_zmumu_mu_eta_loose) > mu_eta) h_mu_pt_loose_endcap->Fill(f_zmumu_mu_pt_loose);
    if(fabs(f_zmumu_mu_eta_tight) <= mu_eta) h_mu_pt_tight_barrel->Fill(f_zmumu_mu_pt_tight);
    if(fabs(f_zmumu_mu_eta_tight) > mu_eta) h_mu_pt_tight_endcap->Fill(f_zmumu_mu_pt_tight);
    
    if(fabs(f_zmumu_e_eta_loose) <= e_eta) h_e_pt_loose_barrel->Fill(f_zmumu_e_pt_loose);
    if(fabs(f_zmumu_e_eta_loose) > e_eta) h_e_pt_loose_endcap->Fill(f_zmumu_e_pt_loose);
    if(fabs(f_zmumu_e_eta_tight) <= e_eta) h_e_pt_tight_barrel->Fill(f_zmumu_e_pt_tight);
    if(fabs(f_zmumu_e_eta_tight) > e_eta) h_e_pt_tight_endcap->Fill(f_zmumu_e_pt_tight);
    
    //eta
    h_mu_eta_loose_barrel->Fill(f_zee_mu_eta_loose);
    h_mu_eta_loose_endcap->Fill(f_zee_mu_eta_loose);
    h_mu_eta_tight_barrel->Fill(f_zee_mu_eta_tight);
    h_mu_eta_tight_endcap->Fill(f_zee_mu_eta_tight);
    
    h_e_eta_loose_barrel->Fill(f_zee_e_eta_loose);
    h_e_eta_loose_endcap->Fill(f_zee_e_eta_loose);
    h_e_eta_tight_barrel->Fill(f_zee_e_eta_tight);
    h_e_eta_tight_endcap->Fill(f_zee_e_eta_tight);
    
    h_mu_eta_loose_barrel->Fill(f_zmumu_mu_eta_loose);
    h_mu_eta_loose_endcap->Fill(f_zmumu_mu_eta_loose);
    h_mu_eta_tight_barrel->Fill(f_zmumu_mu_eta_tight);
    h_mu_eta_tight_endcap->Fill(f_zmumu_mu_eta_tight);
    
    h_e_eta_loose_barrel->Fill(f_zmumu_e_eta_loose);
    h_e_eta_loose_endcap->Fill(f_zmumu_e_eta_loose);
    h_e_eta_tight_barrel->Fill(f_zmumu_e_eta_tight);
    h_e_eta_tight_endcap->Fill(f_zmumu_e_eta_tight);
    
    //2D
    h_e_pt_eta_loose->Fill(f_zee_e_pt_loose, fabs(f_zee_e_eta_loose));
    h_e_pt_eta_loose->Fill(f_zmumu_e_pt_loose, fabs(f_zmumu_e_eta_loose));
    h_e_pt_eta_tight->Fill(f_zee_e_pt_tight, fabs(f_zee_e_eta_tight));
    h_e_pt_eta_tight->Fill(f_zmumu_e_pt_tight, fabs(f_zmumu_e_eta_tight));

    h_mu_pt_eta_loose->Fill(f_zee_mu_pt_loose, fabs(f_zee_mu_eta_loose));
    h_mu_pt_eta_loose->Fill(f_zmumu_mu_pt_loose, fabs(f_zmumu_mu_eta_loose));
    h_mu_pt_eta_tight->Fill(f_zee_mu_pt_tight, fabs(f_zee_mu_eta_tight));
    h_mu_pt_eta_tight->Fill(f_zmumu_mu_pt_tight, fabs(f_zmumu_mu_eta_tight));
  }
  std::cout<<"=============== Data ==================="<<std::endl;
  std::cout<<"Loose electrons = "<<h_e_eta_loose_barrel->Integral()+h_e_eta_loose_endcap->Integral()<<std::endl;
  std::cout<<"Tight electrons = "<<h_e_eta_tight_barrel->Integral()+h_e_eta_tight_endcap->Integral()<<std::endl;
  std::cout<<"Loose muons = "<<h_mu_eta_loose_barrel->Integral()+h_mu_eta_loose_endcap->Integral()<<std::endl;
  std::cout<<"Tight muons = "<<h_mu_eta_tight_barrel->Integral()+h_mu_eta_tight_endcap->Integral()<<std::endl;

  //MC
  for(unsigned int iev=0; iev<nmc; ++iev){
    mcFR->GetEntry(iev);
    
    //pt
    if(fabs(f_zee_mu_eta_loose) <= mu_eta) h_mu_pt_loose_barrel2->Fill(f_zee_mu_pt_loose,f_weight);
    if(fabs(f_zee_mu_eta_loose) > mu_eta) h_mu_pt_loose_endcap2->Fill(f_zee_mu_pt_loose,f_weight);
    if(fabs(f_zee_mu_eta_tight) <= mu_eta) h_mu_pt_tight_barrel2->Fill(f_zee_mu_pt_tight,f_weight);
    if(fabs(f_zee_mu_eta_tight) > mu_eta) h_mu_pt_tight_endcap2->Fill(f_zee_mu_pt_tight,f_weight);
    
    if(fabs(f_zee_e_eta_loose) <= e_eta) h_e_pt_loose_barrel2->Fill(f_zee_e_pt_loose,f_weight);
    if(fabs(f_zee_e_eta_loose) > e_eta) h_e_pt_loose_endcap2->Fill(f_zee_e_pt_loose,f_weight);
    if(fabs(f_zee_e_eta_tight) <= e_eta) h_e_pt_tight_barrel2->Fill(f_zee_e_pt_tight,f_weight);
    if(fabs(f_zee_e_eta_tight) > e_eta) h_e_pt_tight_endcap2->Fill(f_zee_e_pt_tight,f_weight);
    
    if(fabs(f_zmumu_mu_eta_loose) <= mu_eta) h_mu_pt_loose_barrel2->Fill(f_zmumu_mu_pt_loose,f_weight);
    if(fabs(f_zmumu_mu_eta_loose) > mu_eta) h_mu_pt_loose_endcap2->Fill(f_zmumu_mu_pt_loose,f_weight);
    if(fabs(f_zmumu_mu_eta_tight) <= mu_eta) h_mu_pt_tight_barrel2->Fill(f_zmumu_mu_pt_tight,f_weight);
    if(fabs(f_zmumu_mu_eta_tight) > mu_eta) h_mu_pt_tight_endcap2->Fill(f_zmumu_mu_pt_tight,f_weight);
    
    if(fabs(f_zmumu_e_eta_loose) <= e_eta) h_e_pt_loose_barrel2->Fill(f_zmumu_e_pt_loose,f_weight);
    if(fabs(f_zmumu_e_eta_loose) > e_eta) h_e_pt_loose_endcap2->Fill(f_zmumu_e_pt_loose,f_weight);
    if(fabs(f_zmumu_e_eta_tight) <= e_eta) h_e_pt_tight_barrel2->Fill(f_zmumu_e_pt_tight,f_weight);
    if(fabs(f_zmumu_e_eta_tight) > e_eta) h_e_pt_tight_endcap2->Fill(f_zmumu_e_pt_tight,f_weight);
    
    //eta
    h_mu_eta_loose_barrel2->Fill(f_zee_mu_eta_loose,f_weight);
    h_mu_eta_loose_endcap2->Fill(f_zee_mu_eta_loose,f_weight);
    h_mu_eta_tight_barrel2->Fill(f_zee_mu_eta_tight,f_weight);
    h_mu_eta_tight_endcap2->Fill(f_zee_mu_eta_tight,f_weight);
    
    h_e_eta_loose_barrel2->Fill(f_zee_e_eta_loose,f_weight);
    h_e_eta_loose_endcap2->Fill(f_zee_e_eta_loose,f_weight);
    h_e_eta_tight_barrel2->Fill(f_zee_e_eta_tight,f_weight);
    h_e_eta_tight_endcap2->Fill(f_zee_e_eta_tight,f_weight);
    
    h_mu_eta_loose_barrel2->Fill(f_zmumu_mu_eta_loose,f_weight);
    h_mu_eta_loose_endcap2->Fill(f_zmumu_mu_eta_loose,f_weight);
    h_mu_eta_tight_barrel2->Fill(f_zmumu_mu_eta_tight,f_weight);
    h_mu_eta_tight_endcap2->Fill(f_zmumu_mu_eta_tight,f_weight);
    
    h_e_eta_loose_barrel2->Fill(f_zmumu_e_eta_loose,f_weight);
    h_e_eta_loose_endcap2->Fill(f_zmumu_e_eta_loose,f_weight);
    h_e_eta_tight_barrel2->Fill(f_zmumu_e_eta_tight,f_weight);
    h_e_eta_tight_endcap2->Fill(f_zmumu_e_eta_tight,f_weight);
    
    //2D
    h_e_pt_eta_loose2->Fill(f_zee_e_pt_loose, fabs(f_zee_e_eta_loose),f_weight);
    h_e_pt_eta_loose2->Fill(f_zmumu_e_pt_loose, fabs(f_zmumu_e_eta_loose),f_weight);
    h_e_pt_eta_tight2->Fill(f_zee_e_pt_tight, fabs(f_zee_e_eta_tight),f_weight);
    h_e_pt_eta_tight2->Fill(f_zmumu_e_pt_tight, fabs(f_zmumu_e_eta_tight),f_weight);

    h_mu_pt_eta_loose2->Fill(f_zee_mu_pt_loose, fabs(f_zee_mu_eta_loose),f_weight);
    h_mu_pt_eta_loose2->Fill(f_zmumu_mu_pt_loose, fabs(f_zmumu_mu_eta_loose),f_weight);
    h_mu_pt_eta_tight2->Fill(f_zee_mu_pt_tight, fabs(f_zee_mu_eta_tight),f_weight);
    h_mu_pt_eta_tight2->Fill(f_zmumu_mu_pt_tight, fabs(f_zmumu_mu_eta_tight),f_weight);
  }
  std::cout<<"=============== WZ + ttbar ==================="<<std::endl;
  std::cout<<"Loose electrons = "<<h_e_eta_loose_barrel2->Integral()+h_e_eta_loose_endcap2->Integral()<<std::endl;
  std::cout<<"Tight electrons = "<<h_e_eta_tight_barrel2->Integral()+h_e_eta_tight_endcap2->Integral()<<std::endl;
  std::cout<<"Loose muons = "<<h_mu_eta_loose_barrel2->Integral()+h_mu_eta_loose_endcap2->Integral()<<std::endl;
  std::cout<<"Tight muons = "<<h_mu_eta_tight_barrel2->Integral()+h_mu_eta_tight_endcap2->Integral()<<std::endl;

  
  //pt -- uncorrected fake rate
  TH1D *h_e_pt_fr_uncorr_barrel = (TH1D*)h_e_pt_tight_barrel->Clone("h_e_pt_fr_uncorr_barrel");
  h_e_pt_fr_uncorr_barrel->Divide(h_e_pt_loose_barrel);
  TH1D *h_e_pt_fr_uncorr_endcap = (TH1D*)h_e_pt_tight_endcap->Clone("h_e_pt_fr_uncorr_endcap");
  h_e_pt_fr_uncorr_endcap->Divide(h_e_pt_loose_endcap);
  TH1D *h_mu_pt_fr_uncorr_barrel = (TH1D*)h_mu_pt_tight_barrel->Clone("h_mu_pt_fr_uncorr_barrel");
  h_mu_pt_fr_uncorr_barrel->Divide(h_mu_pt_loose_barrel);
  TH1D *h_mu_pt_fr_uncorr_endcap = (TH1D*)h_mu_pt_tight_endcap->Clone("h_mu_pt_fr_uncorr_endcap");
  h_mu_pt_fr_uncorr_endcap->Divide(h_mu_pt_loose_endcap);
    
  //pt -- corrected fake rate
  TH1D *h_e_pt_fr_corr_barrel = (TH1D*)h_e_pt_tight_barrel->Clone("h_e_pt_fr_corr_barrel");
  h_e_pt_fr_corr_barrel->Add(h_e_pt_tight_barrel2,-1);
  h_e_pt_loose_barrel->Add(h_e_pt_loose_barrel2,-1);
  h_e_pt_fr_corr_barrel->Divide(h_e_pt_loose_barrel);
  
  TH1D *h_e_pt_fr_corr_endcap = (TH1D*)h_e_pt_tight_endcap->Clone("h_e_pt_fr_corr_endcap");
  h_e_pt_fr_corr_endcap->Add(h_e_pt_tight_endcap2,-1);
  h_e_pt_loose_endcap->Add(h_e_pt_loose_endcap2,-1);
  h_e_pt_fr_corr_endcap->Divide(h_e_pt_loose_endcap);

  TH1D *h_mu_pt_fr_corr_barrel = (TH1D*)h_mu_pt_tight_barrel->Clone("h_mu_pt_fr_corr_barrel");
  h_mu_pt_fr_corr_barrel->Add(h_mu_pt_tight_barrel2,-1);
  h_mu_pt_loose_barrel->Add(h_mu_pt_loose_barrel2,-1);
  h_mu_pt_fr_corr_barrel->Divide(h_mu_pt_loose_barrel);

  TH1D *h_mu_pt_fr_corr_endcap = (TH1D*)h_mu_pt_tight_endcap->Clone("h_mu_pt_fr_corr_endcap");
  h_mu_pt_fr_corr_endcap->Add(h_mu_pt_tight_endcap2,-1);
  h_mu_pt_loose_endcap->Add(h_mu_pt_loose_endcap2,-1);
  h_mu_pt_fr_corr_endcap->Divide(h_mu_pt_loose_endcap);

  //eta -- uncorrected fake rate
  TH1D *h_e_eta_fr_uncorr_barrel = (TH1D*)h_e_eta_tight_barrel->Clone("h_e_eta_fr_uncorr_barrel");
  h_e_eta_fr_uncorr_barrel->Divide(h_e_eta_loose_barrel);
  TH1D *h_e_eta_fr_uncorr_endcap = (TH1D*)h_e_eta_tight_endcap->Clone("h_e_eta_fr_uncorr_endcap");
  h_e_eta_fr_uncorr_endcap->Divide(h_e_eta_loose_endcap);
  TH1D *h_mu_eta_fr_uncorr_barrel = (TH1D*)h_mu_eta_tight_barrel->Clone("h_mu_eta_fr_uncorr_barrel");
  h_mu_eta_fr_uncorr_barrel->Divide(h_mu_eta_loose_barrel);
  TH1D *h_mu_eta_fr_uncorr_endcap = (TH1D*)h_mu_eta_tight_endcap->Clone("h_mu_eta_fr_uncorr_endcap");
  h_mu_eta_fr_uncorr_endcap->Divide(h_mu_eta_loose_endcap);
    
  //eta -- corrected fake rate
  TH1D *h_e_eta_fr_corr_barrel = (TH1D*)h_e_eta_tight_barrel->Clone("h_e_eta_fr_corr_barrel");
  h_e_eta_fr_corr_barrel->Add(h_e_eta_tight_barrel2,-1);
  h_e_eta_loose_barrel->Add(h_e_eta_loose_barrel2,-1);
  h_e_eta_fr_corr_barrel->Divide(h_e_eta_loose_barrel);
  
  TH1D *h_e_eta_fr_corr_endcap = (TH1D*)h_e_eta_tight_endcap->Clone("h_e_eta_fr_corr_endcap");
  h_e_eta_fr_corr_endcap->Add(h_e_eta_tight_endcap2,-1);
  h_e_eta_loose_endcap->Add(h_e_eta_loose_endcap2,-1);
  h_e_eta_fr_corr_endcap->Divide(h_e_eta_loose_endcap);

  TH1D *h_mu_eta_fr_corr_barrel = (TH1D*)h_mu_eta_tight_barrel->Clone("h_mu_eta_fr_corr_barrel");
  h_mu_eta_fr_corr_barrel->Add(h_mu_eta_tight_barrel2,-1);
  h_mu_eta_loose_barrel->Add(h_mu_eta_loose_barrel2,-1);
  h_mu_eta_fr_corr_barrel->Divide(h_mu_eta_loose_barrel);

  TH1D *h_mu_eta_fr_corr_endcap = (TH1D*)h_mu_eta_tight_endcap->Clone("h_mu_eta_fr_corr_endcap");
  h_mu_eta_fr_corr_endcap->Add(h_mu_eta_tight_endcap2,-1);
  h_mu_eta_loose_endcap->Add(h_mu_eta_loose_endcap2,-1);
  h_mu_eta_fr_corr_endcap->Divide(h_mu_eta_loose_endcap);

  
  gStyle->SetOptStat(0);
  
  ///------------------------ 1D plot -------------------------- 
  
  TCanvas *cv1 = new TCanvas("cv1","",100,100,850,750);
  cv1->Divide(2,2);
  cv1->cd(1);
  h_e_pt_fr_uncorr_barrel->SetMarkerColor(kBlue);
  h_e_pt_fr_uncorr_barrel->SetLineColor(kBlue);
  h_e_pt_fr_uncorr_barrel->GetYaxis()->SetRangeUser(0,0.4);
  h_e_pt_fr_uncorr_barrel->GetXaxis()->SetTitle("p_{T}(GeV)");
  h_e_pt_fr_uncorr_barrel->GetYaxis()->SetTitle("Fake Rate");
  h_e_pt_fr_uncorr_barrel->SetTitle("Electrons");  
  h_e_pt_fr_uncorr_barrel->Draw("pe");
  h_e_pt_fr_uncorr_endcap->SetMarkerColor(kRed);
  h_e_pt_fr_uncorr_endcap->SetLineColor(kRed);
  h_e_pt_fr_uncorr_endcap->Draw("pe,same");
  h_e_pt_fr_corr_barrel->SetMarkerColor(kBlue);
  h_e_pt_fr_corr_barrel->SetLineColor(kBlue);
  h_e_pt_fr_corr_barrel->SetLineStyle(2);
  h_e_pt_fr_corr_barrel->Draw("pe,same");
  h_e_pt_fr_corr_endcap->SetMarkerColor(kRed);
  h_e_pt_fr_corr_endcap->SetLineColor(kRed);
  h_e_pt_fr_corr_endcap->SetLineStyle(2);
  h_e_pt_fr_corr_endcap->Draw("pe,same");
  gPad->SetLogx();
  h_e_pt_fr_uncorr_barrel->GetXaxis()->SetMoreLogLabels();

  TLegend *leg1 = new TLegend(0.2,0.7,0.6,0.9);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(0);
  leg1->AddEntry(h_e_pt_fr_uncorr_barrel,"Barrel uncorrected","pl");
  leg1->AddEntry(h_e_pt_fr_corr_barrel,"Barrel corrected","pl");
  leg1->AddEntry(h_e_pt_fr_uncorr_endcap,"Endcap uncorrected","pl");
  leg1->AddEntry(h_e_pt_fr_corr_endcap,"Endcap corrected","pl");
  leg1->Draw();
  
  cv1->cd(2);  
  h_mu_pt_fr_uncorr_barrel->SetMarkerColor(kBlue);
  h_mu_pt_fr_uncorr_barrel->SetLineColor(kBlue);
  h_mu_pt_fr_uncorr_barrel->GetYaxis()->SetRangeUser(0,0.4);
  h_mu_pt_fr_uncorr_barrel->GetXaxis()->SetTitle("p_{T}(GeV)");
  h_mu_pt_fr_uncorr_barrel->GetYaxis()->SetTitle("Fake Rate");
  h_mu_pt_fr_uncorr_barrel->SetTitle("Muons");
  h_mu_pt_fr_uncorr_barrel->Draw("pe");
  h_mu_pt_fr_uncorr_endcap->SetMarkerColor(kRed);
  h_mu_pt_fr_uncorr_endcap->SetLineColor(kRed);
  h_mu_pt_fr_uncorr_endcap->Draw("pe,same");
  h_mu_pt_fr_corr_barrel->SetMarkerColor(kBlue);
  h_mu_pt_fr_corr_barrel->SetLineColor(kBlue);
  h_mu_pt_fr_corr_barrel->SetLineStyle(2);
  h_mu_pt_fr_corr_barrel->Draw("pe,same");
  h_mu_pt_fr_corr_endcap->SetMarkerColor(kRed);
  h_mu_pt_fr_corr_endcap->SetLineColor(kRed);
  h_mu_pt_fr_corr_endcap->SetLineStyle(2);
  h_mu_pt_fr_corr_endcap->Draw("pe,same");
  gPad->SetLogx();
  h_mu_pt_fr_uncorr_barrel->GetXaxis()->SetMoreLogLabels();
  leg1->Draw();

  cv1->cd(3);
  h_e_eta_fr_uncorr_barrel->SetMarkerColor(kBlue);
  h_e_eta_fr_uncorr_barrel->SetLineColor(kBlue);
  h_e_eta_fr_uncorr_barrel->GetYaxis()->SetRangeUser(0,0.2);
  h_e_eta_fr_uncorr_barrel->GetXaxis()->SetTitle("|#eta|");
  h_e_eta_fr_uncorr_barrel->GetYaxis()->SetTitle("Fake Rate");
  h_e_eta_fr_uncorr_barrel->SetTitle("Electrons");  
  h_e_eta_fr_uncorr_barrel->Draw("pe");
  h_e_eta_fr_corr_barrel->SetMarkerColor(kBlue);
  h_e_eta_fr_corr_barrel->SetLineColor(kBlue);
  h_e_eta_fr_corr_barrel->SetLineStyle(2);
  h_e_eta_fr_corr_barrel->Draw("pe,same");

  TLegend *leg2 = new TLegend(0.2,0.77,0.6,0.9);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->AddEntry(h_e_eta_fr_uncorr_barrel,"Uncorrected","pl");
  leg2->AddEntry(h_e_eta_fr_corr_barrel,"Corrected","pl");
  leg2->Draw();

  cv1->cd(4);
  h_mu_eta_fr_uncorr_barrel->SetMarkerColor(kBlue);
  h_mu_eta_fr_uncorr_barrel->SetLineColor(kBlue);
  h_mu_eta_fr_uncorr_barrel->GetYaxis()->SetRangeUser(0,0.2);
  h_mu_eta_fr_uncorr_barrel->GetXaxis()->SetTitle("|#eta|");
  h_mu_eta_fr_uncorr_barrel->GetYaxis()->SetTitle("Fake Rate");
  h_mu_eta_fr_uncorr_barrel->SetTitle("Muons");  
  h_mu_eta_fr_uncorr_barrel->Draw("pe");
  h_mu_eta_fr_corr_barrel->SetMarkerColor(kBlue);
  h_mu_eta_fr_corr_barrel->SetLineColor(kBlue);
  h_mu_eta_fr_corr_barrel->SetLineStyle(2);
  h_mu_eta_fr_corr_barrel->Draw("pe,same");
  leg2->Draw();
  
  cv1->Update();
  cv1->Print("fake_rate_1D_pt_eta.png");
  
  
  ///------------------------ 2D plot -------------------------- 
  gStyle->SetPaintTextFormat(".2f");
  TCanvas *cv2 = new TCanvas("cv2","",10,10,1250,550);
  cv2->Divide(2,1);
  cv2->cd(1);
  h_e_pt_eta_tight->Add(h_e_pt_eta_tight2,-1);
  h_e_pt_eta_loose->Add(h_e_pt_eta_loose2,-1);
  TH2D *h_e_pt_eta_fr_corr = (TH2D*)h_e_pt_eta_tight->Clone();
  h_e_pt_eta_fr_corr->Divide(h_e_pt_eta_loose);
  h_e_pt_eta_fr_corr->Draw("Colz,text,e");
  h_e_pt_eta_fr_corr->SetTitle("Electrons");
  h_e_pt_eta_fr_corr->GetXaxis()->SetTitle("p_{T}(GeV)");
  h_e_pt_eta_fr_corr->GetYaxis()->SetTitle("#eta");
  gPad->SetLogx();
  gPad->SetRightMargin(0.15);
  h_e_pt_eta_fr_corr->GetXaxis()->SetMoreLogLabels();
  h_e_pt_eta_fr_corr->SetMarkerColor(kBlack);
  h_e_pt_eta_fr_corr->SetMarkerSize(1.4);
  
  cv2->cd(2);
  h_mu_pt_eta_tight->Add(h_mu_pt_eta_tight2,-1);
  h_mu_pt_eta_loose->Add(h_mu_pt_eta_loose2,-1);
  TH2D *h_mu_pt_eta_fr_corr = (TH2D*)h_mu_pt_eta_tight->Clone();
  h_mu_pt_eta_fr_corr->Divide(h_mu_pt_eta_loose);
  h_mu_pt_eta_fr_corr->Draw("Colz,text,e");
  h_mu_pt_eta_fr_corr->SetTitle("Muons");
  h_mu_pt_eta_fr_corr->GetXaxis()->SetTitle("p_{T}(GeV)");
  h_mu_pt_eta_fr_corr->GetYaxis()->SetTitle("#eta");
  gPad->SetLogx();
  gPad->SetRightMargin(0.15);
  h_mu_pt_eta_fr_corr->GetXaxis()->SetMoreLogLabels();
  h_mu_pt_eta_fr_corr->SetMarkerColor(kBlack);
  h_mu_pt_eta_fr_corr->SetMarkerSize(1.4);
  
  cv2->Update();
  cv2->Print("fake_rate_2D_maps_corrected.png");
  
  
  TFile *frmap = new TFile("FakeRatePtEtaMap.root","recreate");
  h_e_pt_eta_fr_corr->SetName("ElectronPtEtaMapFR");
  h_e_pt_eta_fr_corr->SetTitle("Electron Pt vs. Eta Fake Rate Map");
  h_e_pt_eta_fr_corr->Write();
  h_mu_pt_eta_fr_corr->SetName("MuonPtEtaMapFR");
  h_mu_pt_eta_fr_corr->SetTitle("Muon Pt vs. Eta Fake Rate Map");
  h_mu_pt_eta_fr_corr->Write();
  frmap->Close();

  TFile *frmap2 = new TFile("FakeRatePtEtaMapSysUncertainty.root","recreate");
  (*h_e_pt_eta_fr_corr) = (*h_e_pt_eta_fr_corr)*(1.0+ele_sys);
  h_e_pt_eta_fr_corr->Write();
  (*h_mu_pt_eta_fr_corr) = (*h_mu_pt_eta_fr_corr)*(1.0+mu_sys);
  h_mu_pt_eta_fr_corr->Write();
  frmap2->Close();
  
  //END
}
