#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iostream>
#include <fstream>
#include <vector>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TF1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TArrow.h>
#include <TLine.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TPaveText.h>
#include <TSystem.h>
#include <TRandom3.h>
#include "ZZMatrixElement/MEMCalculators/interface/MEMCalculators.h"

//NN models
#include "/lustre/cms/store/user/mmelodea/Keras/modelsNN/k57nj2.h"
#include "/lustre/cms/store/user/mmelodea/Keras/modelsNN/k24nj3.h"



#endif

TString samples_path =  "/lustrehome/mmelodea/MonoHiggs/MonoHiggs/80X/";

using namespace MEMNames;


int loopSamples(std::vector<TString> tags, std::vector< std::vector<TString> > samples, TString tag){
  TH1D::SetDefaultSumw2();

  // Declare MEM class
  MEMs combinedMEM(13,125,"CTEQ6L");  
  
  const int maxJets = 4;
  const int maxJetsComponents = 60;
  const int maxPartons = 8;
  bool f_lhe_parton_clear[maxPartons];
  bool f_outlier;
  int f_run, f_lumi, f_event, f_lept1_pdgid, f_lept2_pdgid, f_lept3_pdgid, f_lept4_pdgid, f_Nbjets, f_lhe_npartons, f_lhe_parton_pdgid[maxPartons];
  int nbjets, f_jets_highpt_charged_hadron_multiplicity[maxJets], f_jets_highpt_neutral_hadron_multiplicity[maxJets], f_jets_highpt_photon_multiplicity[maxJets], f_jets_highpt_electron_multiplicity[maxJets];
  int f_jets_highpt_muon_multiplicity[maxJets], f_jets_highpt_hf_hadron_multiplicity[maxJets], f_jets_highpt_hf_em_multiplicity[maxJets], f_jets_highpt_charged_multiplicity[maxJets], f_jets_highpt_neutral_multiplicity[maxJets];
  int f_jets_highpt_ncomponents[maxJets], f_jets_highpt_component_pdgid[maxJets][maxJetsComponents];
  float f_lept1_pt, f_lept1_eta, f_lept1_phi, f_lept1_pfx, f_lept1_sip, f_lept2_pt, f_lept2_eta, f_lept2_phi, f_lept2_pfx, f_lept2_sip;
  float f_lept3_pt, f_lept3_eta, f_lept3_phi, f_lept3_pfx, f_lept3_sip, f_lept4_pt, f_lept4_eta, f_lept4_phi, f_lept4_pfx, f_lept4_sip;
  float f_weight, f_Djet_VAJHU, f_iso_max, f_sip_max, f_Z1mass, f_Z2mass, f_angle_costhetastar, f_angle_costheta1, f_angle_costheta2, f_angle_phi, f_angle_phistar1;
  float f_pt4l, f_eta4l, f_mass4l, f_deltajj, f_massjj, f_pfmet, f_pfmet_theta, f_pfmet_phi, f_mT, f_dphi;
  float f_jets_highpt_pt[maxJets], f_jets_highpt_eta[maxJets], f_jets_highpt_phi[maxJets], f_jets_highpt_et[maxJets], f_jets_highpt_area[maxJets], f_jets_highpt_ptd[maxJets], f_jets_highpt_charged_hadron_energy[maxJets];
  float f_jets_highpt_neutral_hadron_energy[maxJets], f_jets_highpt_photon_energy[maxJets], f_jets_highpt_electron_energy[maxJets], f_jets_highpt_muon_energy[maxJets], f_jets_highpt_hf_hadron_energy[maxJets];
  float f_jets_highpt_hf_em_energy[maxJets], f_jets_highpt_charged_em_energy[maxJets], f_jets_highpt_charged_mu_energy[maxJets], f_jets_highpt_neutral_em_energy[maxJets], f_jets_highpt_component_mt[maxJets][maxJetsComponents];
  float f_jets_highpt_component_pt[maxJets][maxJetsComponents], f_jets_highpt_component_eta[maxJets][maxJetsComponents], f_jets_highpt_component_phi[maxJets][maxJetsComponents], f_jets_highpt_component_energy[maxJets][maxJetsComponents], f_jets_highpt_component_xvertex[maxJets][maxJetsComponents];
  float f_jets_highpt_component_yvertex[maxJets][maxJetsComponents], f_jets_highpt_component_zvertex[maxJets][maxJetsComponents], f_jets_highpt_component_vertex_chi2[maxJets][maxJetsComponents];
  float f_int_weight, f_pu_weight, f_eff_weight, f_mass4lErr, f_Djet_VAJHU_UncDn, f_Djet_VAJHU_UncUp;
  float f_lept1_charge, f_lept2_charge, f_lept3_charge, f_lept4_charge, f_njets_pass, f_jets_highpt_component_charge[maxJets][maxJetsComponents];
  float f_lept1_pt_error, f_lept2_pt_error, f_lept3_pt_error, f_lept4_pt_error, f_jets_highpt_pt_error[maxJets], f_jets_highpt_btagger[maxJets];
  float f_pfmet_JetEnUp, f_pfmet_JetEnDn, f_pfmet_ElectronEnUp, f_pfmet_ElectronEnDn, f_pfmet_MuonEnUp, f_pfmet_MuonEnDn,f_pfmet_JetResUp, f_pfmet_JetResDn, f_pfmet_UnclusteredEnUp, f_pfmet_UnclusteredEnDn, f_pfmet_PhotonEnUp, f_pfmet_PhotonEnDn;
  float f_lhe_parton_pt[maxPartons], f_lhe_parton_eta[maxPartons], f_lhe_parton_phi[maxPartons], f_lhe_parton_e[maxPartons];
  float f_lept1_pt_shifts[500], f_lept2_pt_shifts[500], f_lept3_pt_shifts[500], f_lept4_pt_shifts[500], f_jets_highpt_pt_shifts[maxJets][500], f_pfmet_shifts[500];

  std::vector<TString> hfs = {"ch4mu","ch4e","ch2e2mu"};
  unsigned int nfs = hfs.size();
  std::vector<TString> hmodels = {"MELA","k57nj2","k24nj3"};
  unsigned int nmodels = hmodels.size();
  
  std::vector<TString> shift_name = {"","_l1ptDown","_l1ptUp","_l2ptDown","_l2ptUp","_l3ptDown","_l3ptUp","_l4ptDown","_l4ptUp",
				     "_j1ptDown","_j1ptUp","_j2ptDown","_j2ptUp","_j3ptDown","_j3ptUp",
				     "_metEleEnDown","_metEleEnUp","_metJetEnDown","_metJetEnUp","_metJetResDown","_metJetResUp",
				     "_metMuEnDown","_metMuEnUp","_metPhoEnDown","_metPhoEnUp","_metUncEnDown","_metUncEnUp"};
  const int nshifts = shift_name.size();
				  
  if(nshifts != shift_name.size()) std::cout<<"NUMBER OF SHIFTS != NUMBER OF NAMES"<<std::endl;
  TH1D* distributions[nmodels][nfs][nshifts];
  for(unsigned int i=0; i<nmodels; ++i){
    for(unsigned int j=0; j<nfs; ++j){
      for(int k=0; k<nshifts; ++k){
	distributions[i][j][k] = new TH1D(hmodels[i]+hfs[j]+tag+shift_name[k],"",120,0,1);
      }
    }
  }
  
  std::cout<<"Path: "<<samples_path<<std::endl;
  unsigned int ntags = tags.size();
  for(unsigned int itag=0; itag<ntags; ++itag){
    if(tags[itag] != tag) continue;
      
    unsigned int nsamples = samples[itag].size();
    for(unsigned int is=0; is<nsamples; ++is){
      TString ifile_name_root = samples_path+samples[itag][is];
      if(gSystem->AccessPathName(ifile_name_root)){
	cout<<"File "<<ifile_name_root<<" doesn't exist!"<<endl;
	continue;
      }
      else cout<<"Loading file "<<ifile_name_root<<endl;
      unsigned int ifs = 0;
      if(ifile_name_root.Contains("histos4mu_25ns"))   ifs = 0;
      if(ifile_name_root.Contains("histos4e_25ns"))    ifs = 1;
      if(ifile_name_root.Contains("histos2e2mu_25ns")) ifs = 2;
      //std::cout<<">>>>>>>>>> Final state = "<<hfs[ifs]<<std::endl;
  
      TFile *ofile = TFile::Open(ifile_name_root);
      TTree *otree = (TTree*)ofile->Get("HZZ4LeptonsAnalysisReduced");
      otree->SetBranchAddress("f_run", &f_run);
      otree->SetBranchAddress("f_lumi", &f_lumi);    
      otree->SetBranchAddress("f_event", &f_event);    
      otree->SetBranchAddress("f_weight", &f_weight);
      otree->SetBranchAddress("f_int_weight", &f_int_weight);
      otree->SetBranchAddress("f_pu_weight", &f_pu_weight);
      otree->SetBranchAddress("f_eff_weight", &f_eff_weight);
      otree->SetBranchAddress("f_lept1_pt", &f_lept1_pt);
      otree->SetBranchAddress("f_lept1_pt_error", &f_lept1_pt_error);
      otree->SetBranchAddress("f_lept1_eta", &f_lept1_eta);
      otree->SetBranchAddress("f_lept1_phi", &f_lept1_phi);
      otree->SetBranchAddress("f_lept1_charge", &f_lept1_charge);
      otree->SetBranchAddress("f_lept1_pfx", &f_lept1_pfx);
      otree->SetBranchAddress("f_lept1_sip", &f_lept1_sip);
      otree->SetBranchAddress("f_lept1_pdgid", &f_lept1_pdgid);
      otree->SetBranchAddress("f_lept2_pt", &f_lept2_pt);
      otree->SetBranchAddress("f_lept2_pt_error", &f_lept2_pt_error);
      otree->SetBranchAddress("f_lept2_eta", &f_lept2_eta);
      otree->SetBranchAddress("f_lept2_phi", &f_lept2_phi);
      otree->SetBranchAddress("f_lept2_charge", &f_lept2_charge);
      otree->SetBranchAddress("f_lept2_pfx", &f_lept2_pfx);
      otree->SetBranchAddress("f_lept2_sip", &f_lept2_sip);
      otree->SetBranchAddress("f_lept2_pdgid", &f_lept2_pdgid);
      otree->SetBranchAddress("f_lept3_pt", &f_lept3_pt);
      otree->SetBranchAddress("f_lept3_pt_error", &f_lept3_pt_error);
      otree->SetBranchAddress("f_lept3_eta", &f_lept3_eta);
      otree->SetBranchAddress("f_lept3_phi", &f_lept3_phi);
      otree->SetBranchAddress("f_lept3_charge", &f_lept3_charge);
      otree->SetBranchAddress("f_lept3_pfx", &f_lept3_pfx);
      otree->SetBranchAddress("f_lept3_sip", &f_lept3_sip);
      otree->SetBranchAddress("f_lept3_pdgid", &f_lept3_pdgid);
      otree->SetBranchAddress("f_lept4_pt", &f_lept4_pt);
      otree->SetBranchAddress("f_lept4_pt_error", &f_lept4_pt_error);
      otree->SetBranchAddress("f_lept4_eta", &f_lept4_eta);
      otree->SetBranchAddress("f_lept4_phi", &f_lept4_phi);
      otree->SetBranchAddress("f_lept4_charge", &f_lept4_charge);
      otree->SetBranchAddress("f_lept4_pfx", &f_lept4_pfx);
      otree->SetBranchAddress("f_lept4_sip", &f_lept4_sip);
      otree->SetBranchAddress("f_lept4_pdgid", &f_lept4_pdgid);
      otree->SetBranchAddress("f_iso_max", &f_iso_max);
      otree->SetBranchAddress("f_sip_max", &f_sip_max);
      otree->SetBranchAddress("f_Z1mass", &f_Z1mass);
      otree->SetBranchAddress("f_Z2mass", &f_Z2mass);
      otree->SetBranchAddress("f_angle_costhetastar", &f_angle_costhetastar);
      otree->SetBranchAddress("f_angle_costheta1", &f_angle_costheta1);
      otree->SetBranchAddress("f_angle_costheta2", &f_angle_costheta2);
      otree->SetBranchAddress("f_angle_phi", &f_angle_phi);
      otree->SetBranchAddress("f_angle_phistar1", &f_angle_phistar1);
      otree->SetBranchAddress("f_pt4l", &f_pt4l);
      otree->SetBranchAddress("f_eta4l", &f_eta4l);
      otree->SetBranchAddress("f_mass4l", &f_mass4l);
      otree->SetBranchAddress("f_mass4lErr", &f_mass4lErr);
      otree->SetBranchAddress("f_njets_pass", &f_njets_pass);
      otree->SetBranchAddress("f_deltajj", &f_deltajj);
      otree->SetBranchAddress("f_massjj", &f_massjj);
      otree->SetBranchAddress("f_Djet_VAJHU", &f_Djet_VAJHU);
      otree->SetBranchAddress("f_Djet_VAJHU_UncUp", &f_Djet_VAJHU_UncUp);
      otree->SetBranchAddress("f_Djet_VAJHU_UncDn", &f_Djet_VAJHU_UncDn);
      otree->SetBranchAddress("f_pfmet", &f_pfmet);
      otree->SetBranchAddress("f_pfmet_theta", &f_pfmet_theta);
      otree->SetBranchAddress("f_pfmet_phi", &f_pfmet_phi);
      otree->SetBranchAddress("f_pfmet_JetEnUp", &f_pfmet_JetEnUp);
      otree->SetBranchAddress("f_pfmet_JetEnDn", &f_pfmet_JetEnDn);
      otree->SetBranchAddress("f_pfmet_ElectronEnUp", &f_pfmet_ElectronEnUp);
      otree->SetBranchAddress("f_pfmet_ElectronEnDn", &f_pfmet_ElectronEnDn);
      otree->SetBranchAddress("f_pfmet_MuonEnUp", &f_pfmet_MuonEnUp);
      otree->SetBranchAddress("f_pfmet_MuonEnDn", &f_pfmet_MuonEnDn);
      otree->SetBranchAddress("f_pfmet_JetResUp", &f_pfmet_JetResUp);
      otree->SetBranchAddress("f_pfmet_JetResDn", &f_pfmet_JetResDn);
      otree->SetBranchAddress("f_pfmet_UnclusteredEnUp", &f_pfmet_UnclusteredEnUp);
      otree->SetBranchAddress("f_pfmet_UnclusteredEnDn", &f_pfmet_UnclusteredEnDn);
      otree->SetBranchAddress("f_pfmet_PhotonEnUp", &f_pfmet_PhotonEnUp);
      otree->SetBranchAddress("f_pfmet_PhotonEnDn", &f_pfmet_PhotonEnDn);
      otree->SetBranchAddress("f_mT", &f_mT);
      otree->SetBranchAddress("f_dphi", &f_dphi);
      otree->SetBranchAddress("f_Nbjets", &f_Nbjets);
      otree->SetBranchAddress("f_jets_highpt_btagger", &f_jets_highpt_btagger);
      otree->SetBranchAddress("f_jets_highpt_pt", &f_jets_highpt_pt);
      otree->SetBranchAddress("f_jets_highpt_pt_error", &f_jets_highpt_pt_error);
      otree->SetBranchAddress("f_jets_highpt_eta", &f_jets_highpt_eta);
      otree->SetBranchAddress("f_jets_highpt_phi", &f_jets_highpt_phi);
      otree->SetBranchAddress("f_jets_highpt_et", &f_jets_highpt_et);
      otree->SetBranchAddress("f_jets_highpt_area", &f_jets_highpt_area);
      otree->SetBranchAddress("f_jets_highpt_ptd", &f_jets_highpt_ptd);
      otree->SetBranchAddress("f_jets_highpt_charged_hadron_energy", &f_jets_highpt_charged_hadron_energy);
      otree->SetBranchAddress("f_jets_highpt_neutral_hadron_energy", &f_jets_highpt_neutral_hadron_energy);
      otree->SetBranchAddress("f_jets_highpt_photon_energy", &f_jets_highpt_photon_energy);
      otree->SetBranchAddress("f_jets_highpt_electron_energy", &f_jets_highpt_electron_energy);
      otree->SetBranchAddress("f_jets_highpt_muon_energy", &f_jets_highpt_muon_energy);
      otree->SetBranchAddress("f_jets_highpt_hf_hadron_energy", &f_jets_highpt_hf_hadron_energy);
      otree->SetBranchAddress("f_jets_highpt_hf_em_energy", &f_jets_highpt_hf_em_energy);
      otree->SetBranchAddress("f_jets_highpt_charged_em_energy", &f_jets_highpt_charged_em_energy);
      otree->SetBranchAddress("f_jets_highpt_charged_mu_energy", &f_jets_highpt_charged_mu_energy);
      otree->SetBranchAddress("f_jets_highpt_neutral_em_energy", &f_jets_highpt_neutral_em_energy);
      otree->SetBranchAddress("f_jets_highpt_charged_hadron_multiplicity", &f_jets_highpt_charged_hadron_multiplicity);
      otree->SetBranchAddress("f_jets_highpt_neutral_hadron_multiplicity", &f_jets_highpt_neutral_hadron_multiplicity);
      otree->SetBranchAddress("f_jets_highpt_photon_multiplicity", &f_jets_highpt_photon_multiplicity);
      otree->SetBranchAddress("f_jets_highpt_electron_multiplicity", &f_jets_highpt_electron_multiplicity);
      otree->SetBranchAddress("f_jets_highpt_muon_multiplicity", &f_jets_highpt_muon_multiplicity);
      otree->SetBranchAddress("f_jets_highpt_hf_hadron_multiplicity", &f_jets_highpt_hf_hadron_multiplicity);
      otree->SetBranchAddress("f_jets_highpt_hf_em_multiplicity", &f_jets_highpt_hf_em_multiplicity);
      otree->SetBranchAddress("f_jets_highpt_charged_multiplicity", &f_jets_highpt_charged_multiplicity);
      otree->SetBranchAddress("f_jets_highpt_neutral_multiplicity", &f_jets_highpt_neutral_multiplicity);
      otree->SetBranchAddress("f_jets_highpt_ncomponents", &f_jets_highpt_ncomponents);
      otree->SetBranchAddress("f_jets_highpt_component_pdgid", &f_jets_highpt_component_pdgid);
      otree->SetBranchAddress("f_jets_highpt_component_pt", &f_jets_highpt_component_pt);
      otree->SetBranchAddress("f_jets_highpt_component_eta", &f_jets_highpt_component_eta);
      otree->SetBranchAddress("f_jets_highpt_component_phi", &f_jets_highpt_component_phi);
      otree->SetBranchAddress("f_jets_highpt_component_energy", &f_jets_highpt_component_energy);
      otree->SetBranchAddress("f_jets_highpt_component_charge", &f_jets_highpt_component_charge);
      otree->SetBranchAddress("f_jets_highpt_component_mt", &f_jets_highpt_component_mt);
      otree->SetBranchAddress("f_jets_highpt_component_xvertex", &f_jets_highpt_component_xvertex);
      otree->SetBranchAddress("f_jets_highpt_component_yvertex", &f_jets_highpt_component_yvertex);
      otree->SetBranchAddress("f_jets_highpt_component_zvertex", &f_jets_highpt_component_zvertex);
      otree->SetBranchAddress("f_jets_highpt_component_vertex_chi2", &f_jets_highpt_component_vertex_chi2);
      otree->SetBranchAddress("f_lhe_npartons", &f_lhe_npartons);
      otree->SetBranchAddress("f_lhe_parton_clear", &f_lhe_parton_clear);
      otree->SetBranchAddress("f_lhe_parton_pdgid", &f_lhe_parton_pdgid);
      otree->SetBranchAddress("f_lhe_parton_pt", &f_lhe_parton_pt);
      otree->SetBranchAddress("f_lhe_parton_eta", &f_lhe_parton_eta);
      otree->SetBranchAddress("f_lhe_parton_phi", &f_lhe_parton_phi);
      otree->SetBranchAddress("f_lhe_parton_e", &f_lhe_parton_e);
      
      Int_t Nentries = otree->GetEntries();

      for(int ievent=0; ievent<Nentries; ievent++){
	otree->GetEntry(ievent);
      
	if( (((f_njets_pass == 2 || f_njets_pass == 3) && f_Nbjets <= 1) || (f_njets_pass > 3 && f_Nbjets == 0)) && f_mass4l >= 118 && f_mass4l <= 130 ){
	  
	  //if(f_outlier) continue;
	  float l1pt  = f_lept1_pt;
	  float l1pt_err  = f_lept1_pt_error;
	  float l1eta = f_lept1_eta;
	  float l1phi = f_lept1_phi;
	  float l2pt  = f_lept2_pt;
	  float l2pt_err  = f_lept2_pt_error;
	  float l2eta = f_lept2_eta;
	  float l2phi = f_lept2_phi;
	  float l3pt  = f_lept3_pt;
	  float l3pt_err  = f_lept3_pt_error;
	  float l3eta = f_lept3_eta;
	  float l3phi = f_lept3_phi;
	  float l4pt  = f_lept4_pt;
	  float l4pt_err  = f_lept4_pt_error;
	  float l4eta = f_lept4_eta;
	  float l4phi = f_lept4_phi;
	
	  float j1pt  = f_jets_highpt_pt[0];
	  float j1pt_err  = f_jets_highpt_pt_error[0];
	  float j1eta = f_jets_highpt_eta[0];
	  float j1phi = f_jets_highpt_phi[0];
	  float j1e   = f_jets_highpt_et[0]*TMath::CosH(f_jets_highpt_eta[0]);

	  float j2pt  = f_jets_highpt_pt[1];
	  float j2pt_err  = f_jets_highpt_pt_error[1];
	  float j2eta = f_jets_highpt_eta[1];
	  float j2phi = f_jets_highpt_phi[1];
	  float j2e   = f_jets_highpt_et[1]*TMath::CosH(f_jets_highpt_eta[1]);

	  float j3pt  = (f_jets_highpt_pt[2] != -999)? f_jets_highpt_pt[2]: 0;
	  float j3pt_err  = (f_jets_highpt_pt[2] != -999)? f_jets_highpt_pt_error[2]: 0;
	  float j3eta = (f_jets_highpt_pt[2] != -999)? f_jets_highpt_eta[2]: 0;
	  float j3phi = (f_jets_highpt_pt[2] != -999)? f_jets_highpt_phi[2]: 0;
	  float j3e   = f_jets_highpt_et[2]*TMath::CosH(f_jets_highpt_eta[2]);
	  
	  float met     = f_pfmet;
	  float met_eleEnDn = f_pfmet_ElectronEnDn;
	  float met_eleEnUp = f_pfmet_ElectronEnUp;
	  float met_muEnDn = f_pfmet_MuonEnDn;
	  float met_muEnUp = f_pfmet_MuonEnUp;
	  float met_jetEnDn = f_pfmet_JetEnDn;
	  float met_jetEnUp = f_pfmet_JetEnUp;
	  float met_jetResDn = f_pfmet_JetResDn;
	  float met_jetResUp = f_pfmet_JetResUp;
	  float met_uncEnDn = f_pfmet_UnclusteredEnDn;
	  float met_uncEnUp = f_pfmet_UnclusteredEnUp;
	  float met_phoEnDn = f_pfmet_PhotonEnDn;
	  float met_phoEnUp = f_pfmet_PhotonEnUp;
	  float njets  = f_njets_pass;
	  float nbjets = f_Nbjets;
	  
	  
	  std::vector<float> vl1pt = {l1pt,l1pt-l1pt_err,l1pt+l1pt_err};
	  std::vector<float> vl2pt = {l2pt,l2pt-l2pt_err,l2pt+l2pt_err};
	  std::vector<float> vl3pt = {l3pt,l3pt-l3pt_err,l3pt+l3pt_err};
	  std::vector<float> vl4pt = {l4pt,l4pt-l4pt_err,l4pt+l4pt_err};
	  std::vector<float> vj1pt = {j1pt,j1pt-j1pt_err,j1pt+j1pt_err};
	  std::vector<float> vj2pt = {j2pt,j2pt-j2pt_err,j2pt+j2pt_err};
	  std::vector<float> vj3pt = {j3pt,j3pt-j3pt_err,j3pt+j3pt_err};
	  std::vector<float> vmet = {met,met_eleEnDn,met_eleEnUp,met_jetEnDn,met_jetEnUp,met_jetResDn,met_jetResUp,
				     met_muEnDn,met_muEnUp,met_phoEnDn,met_phoEnUp,met_uncEnDn,met_uncEnUp};
				     
				     
	//Loop over +/- sigma on each input
	std::vector< std::vector<float> > shifts = {vl1pt,vl2pt,vl3pt,vl4pt,vj1pt,vj2pt,vj3pt,vmet};
	unsigned int ishift = 0;
	for(unsigned int i=0; i<(unsigned int)shifts.size(); ++i){
	  for(unsigned int j=0; j<(unsigned int)shifts[i].size(); ++j){
	    if(i>0 && j==0) continue; //Prevents nominal distributions to be saved more than once
	    //Reset the possible shifting inputs to their nominal value
	    l1pt  = f_lept1_pt;
	    l2pt  = f_lept2_pt;
	    l3pt  = f_lept3_pt;
	    l4pt  = f_lept4_pt;
	    j1pt  = f_jets_highpt_pt[0];
	    j2pt  = f_jets_highpt_pt[1];
	    j3pt  = (f_jets_highpt_pt[2] != -999)? f_jets_highpt_pt[2]: 0;
	    met   = f_pfmet;
	    
	    if(i==0) l1pt = shifts[i][j];
	    if(i==1) l2pt = shifts[i][j];
	    if(i==2) l3pt = shifts[i][j];
	    if(i==3) l4pt = shifts[i][j];
	    if(i==4) j1pt = shifts[i][j];
	    if(i==5) j2pt = shifts[i][j];
	    if(i==6) j3pt = shifts[i][j];
	    if(i==7) met  = shifts[i][j];
	    
	//MELA------------------------------------------------------------------------------    
	    TLorentzVector partP;
	    vector<TLorentzVector> vpartsP;
	    vector<int> partsId;
       
	    partsId.push_back(f_lept1_pdgid);
	    partsId.push_back(f_lept2_pdgid);
	    partsId.push_back(f_lept3_pdgid);
	    partsId.push_back(f_lept4_pdgid);
	    partsId.push_back(0);
	    partsId.push_back(0);
	 
	    float mass = 0;
	    mass = fabs(f_lept1_pdgid) == 13 ? 0.105 : 0.000511;
	    partP.SetPtEtaPhiM(l1pt, l1eta, l1phi, mass); vpartsP.push_back(partP);
	    mass = fabs(f_lept2_pdgid) == 13 ? 0.105 : 0.000511;
	    partP.SetPtEtaPhiM(l2pt, l2eta, l2phi, mass); vpartsP.push_back(partP);
	    mass = fabs(f_lept3_pdgid) == 13 ? 0.105 : 0.000511;
	    partP.SetPtEtaPhiM(l3pt, l3eta, l3phi, mass); vpartsP.push_back(partP);
	    mass = fabs(f_lept4_pdgid) == 13 ? 0.105 : 0.000511;
	    partP.SetPtEtaPhiM(l4pt, l4eta, l4phi, mass); vpartsP.push_back(partP);

	    partP.SetPtEtaPhiE(j1pt, j1eta, j1phi, j1e); vpartsP.push_back(partP);
	    partP.SetPtEtaPhiE(j2pt, j2eta, j2phi, j2e); vpartsP.push_back(partP);

	    double phjj_VAJHU, pvbf_VAJHU;
	    int f=combinedMEM.computeME(MEMNames::kJJ_SMHiggs_GG, MEMNames::kJHUGen, vpartsP, partsId, phjj_VAJHU); // SM gg->H+2j
	    int g=combinedMEM.computeME(MEMNames::kJJ_SMHiggs_VBF, MEMNames::kJHUGen, vpartsP, partsId, pvbf_VAJHU);  // SM VBF->H	   
	    float Djet_VAJHU = pvbf_VAJHU / ( pvbf_VAJHU + phjj_VAJHU ); // D^VBF_HJJ
            //std::cout<<"Event = "<<ievent<<" --- MELA = "<<Djet_VAJHU<<std::endl;
	    distributions[0][ifs][ishift]->Fill(Djet_VAJHU, f_weight);
            
            //------------------- Neural Networks -----------------------------------------  
            
	    std::vector<std::vector<double> > inputs = {
	      //k57nj2
	      {l1pt,l1eta,l1phi,l2pt,l2eta,l2phi,l3pt,l3eta,l3phi,l4pt,l4eta,l4phi,j1pt,j1eta,j1phi,j2pt,j2eta,j2phi,met},
	      //k24nj3
	      {l1pt,l1eta,l1phi,l2pt,l2eta,l2phi,l3pt,l3eta,l3phi,l4pt,l4eta,l4phi,j1pt,j1eta,j1phi,j2pt,j2eta,j2phi,j3pt,j3eta,j3phi,met}
	    };
            
	    if(f_njets_pass == 2){
	      distributions[1][ifs][ishift]->Fill(model_k57nj2(inputs[0]), f_weight);
	    }
	    if(f_njets_pass >= 3){
	      distributions[2][ifs][ishift]->Fill(model_k24nj3(inputs[1]), f_weight);
	    }
	    
	    ++ishift;
	  //---------------------------------------------------------------------------------//
	  }
	}//End of loop over shifts 
	}//End of "if" for VBF selection
      }            
      ofile->Close();
      //std::cout<<"Processed "<<ifile_name_root<<std::endl;
    }//end samples loop
  }//end tags loop
    
  std::cout<<"Saving histograms..."<<std::endl;
  for(unsigned int i=0; i<nmodels; ++i){
    TFile *dfile = new TFile("DistributionsUncertainties_"+hmodels[i]+"_proc_"+tag+".root","recreate");
    for(unsigned int j=0; j<nfs; ++j){
      dfile->cd();
      dfile->mkdir(hfs[j]);
      dfile->cd(hfs[j]);
      std::cout<<"Saving "<<hmodels[i]<<", "<<hfs[j]<<", "<<tag<<", Integral = "<<distributions[i][j][0]->Integral()<<std::endl;
      for(int k=0; k<nshifts; ++k){
	distributions[i][j][k]->SetName(tag+shift_name[k]);
	distributions[i][j][k]->Write();
      }
    }
    dfile->Close(); 
  }
  
  return 0;
}//End plot function


int main(int argc, char ** argv){
  
  std::vector<TString> data_obs;
  std::vector<TString> qqH;
  std::vector<TString> qqH3J;
  std::vector<TString> ggH;
  std::vector<TString> ggZZ;
  std::vector<TString> qqZZ;
  std::vector<TString> WH;
  std::vector<TString> ZH;
  std::vector<TString> ttH;
  std::vector<TString> HWW;
  std::vector<TString> TTV;
  std::vector<TString> TT4l;
  std::vector<TString> Wjets;
  std::vector<TString> VV;
  std::vector<TString> VVV;
    
  data_obs.push_back("histos4mu_25ns/output_data_obs.root");
  //DY.push_back("histos4mu_25ns/output_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  ggH.push_back("histos4mu_25ns/output_GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8.root");
  //ggH.push_back("histos4mu_25ns/output_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root");
  HWW.push_back("histos4mu_25ns/output_GluGluZH_HToWWTo2L2Nu_ZTo2L_M125_13TeV_powheg_pythia8.root");
  HWW.push_back("histos4mu_25ns/output_HWminusJ_HToWWTo2L2Nu_WToLNu_M125_13TeV_powheg_pythia8.root");
  HWW.push_back("histos4mu_25ns/output_HWplusJ_HToWWTo2L2Nu_WToLNu_M125_13TeV_powheg_pythia8.root");
  HWW.push_back("histos4mu_25ns/output_HZJ_HToWWTo2L2Nu_ZTo2L_M125_13TeV_powheg_pythia8.root");
  TT4l.push_back("histos4mu_25ns/output_TTTo2L2Nu_noSC_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
  TTV.push_back("histos4mu_25ns/output_TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root");
  TTV.push_back("histos4mu_25ns/output_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  HWW.push_back("histos4mu_25ns/output_VBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8.root");
  qqH.push_back("histos4mu_25ns/output_VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  qqH3J.push_back("histos4mu_25ns/output_VBF_HJJJ_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  Wjets.push_back("histos4mu_25ns/output_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  VV.push_back("histos4mu_25ns/output_WWTo2L2Nu_13TeV-powheg.root");
  VVV.push_back("histos4mu_25ns/output_WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  VV.push_back("histos4mu_25ns/output_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root");
  VVV.push_back("histos4mu_25ns/output_WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  WH.push_back("histos4mu_25ns/output_WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root");
  WH.push_back("histos4mu_25ns/output_WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root");
  ZH.push_back("histos4mu_25ns/output_ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8.root");
  qqZZ.push_back("histos4mu_25ns/output_ZZJJTo4L_EWK_13TeV-madgraph-pythia8.root");
  qqZZ.push_back("histos4mu_25ns/output_ZZTo2L2Nu_13TeV_powheg_pythia8.root");
  qqZZ.push_back("histos4mu_25ns/output_ZZTo4L_13TeV_powheg_pythia8.root");
  VVV.push_back("histos4mu_25ns/output_ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  HWW.push_back("histos4mu_25ns/output_bbHToWWTo2L2Nu_M-125_4FS_yb2_13TeV_amcatnlo.root");
  HWW.push_back("histos4mu_25ns/output_bbHToWWTo2L2Nu_M-125_4FS_ybyt_13TeV_amcatnlo.root");
  ttH.push_back("histos4mu_25ns/output_ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUgenV6_pythia8.root");

  data_obs.push_back("histos4e_25ns/output_data_obs.root");
  //DY.push_back("histos4e_25ns/output_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  ggH.push_back("histos4e_25ns/output_GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8.root");
  //ggH.push_back("histos4e_25ns/output_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root");
  HWW.push_back("histos4e_25ns/output_GluGluZH_HToWWTo2L2Nu_ZTo2L_M125_13TeV_powheg_pythia8.root");
  HWW.push_back("histos4e_25ns/output_HWminusJ_HToWWTo2L2Nu_WToLNu_M125_13TeV_powheg_pythia8.root");
  HWW.push_back("histos4e_25ns/output_HWplusJ_HToWWTo2L2Nu_WToLNu_M125_13TeV_powheg_pythia8.root");
  HWW.push_back("histos4e_25ns/output_HZJ_HToWWTo2L2Nu_ZTo2L_M125_13TeV_powheg_pythia8.root");
  TT4l.push_back("histos4e_25ns/output_TTTo2L2Nu_noSC_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
  TTV.push_back("histos4e_25ns/output_TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root");
  TTV.push_back("histos4e_25ns/output_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  HWW.push_back("histos4e_25ns/output_VBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8.root");
  qqH.push_back("histos4e_25ns/output_VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  qqH3J.push_back("histos4e_25ns/output_VBF_HJJJ_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  Wjets.push_back("histos4e_25ns/output_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  VV.push_back("histos4e_25ns/output_WWTo2L2Nu_13TeV-powheg.root");
  VVV.push_back("histos4e_25ns/output_WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  VV.push_back("histos4e_25ns/output_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root");
  VVV.push_back("histos4e_25ns/output_WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  WH.push_back("histos4e_25ns/output_WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root");
  WH.push_back("histos4e_25ns/output_WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root");
  ZH.push_back("histos4e_25ns/output_ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8.root");
  qqZZ.push_back("histos4e_25ns/output_ZZJJTo4L_EWK_13TeV-madgraph-pythia8.root");
  qqZZ.push_back("histos4e_25ns/output_ZZTo2L2Nu_13TeV_powheg_pythia8.root");
  qqZZ.push_back("histos4e_25ns/output_ZZTo4L_13TeV_powheg_pythia8.root");
  VVV.push_back("histos4e_25ns/output_ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  HWW.push_back("histos4e_25ns/output_bbHToWWTo2L2Nu_M-125_4FS_yb2_13TeV_amcatnlo.root");
  HWW.push_back("histos4e_25ns/output_bbHToWWTo2L2Nu_M-125_4FS_ybyt_13TeV_amcatnlo.root");
  ttH.push_back("histos4e_25ns/output_ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUgenV6_pythia8.root");

  data_obs.push_back("histos2e2mu_25ns/output_data_obs.root");
  //DY.push_back("histos2e2mu_25ns/output_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  ggH.push_back("histos2e2mu_25ns/output_GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8.root");
  //ggH.push_back("histos2e2mu_25ns/output_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root");
  HWW.push_back("histos2e2mu_25ns/output_GluGluZH_HToWWTo2L2Nu_ZTo2L_M125_13TeV_powheg_pythia8.root");
  HWW.push_back("histos2e2mu_25ns/output_HWminusJ_HToWWTo2L2Nu_WToLNu_M125_13TeV_powheg_pythia8.root");
  HWW.push_back("histos2e2mu_25ns/output_HWplusJ_HToWWTo2L2Nu_WToLNu_M125_13TeV_powheg_pythia8.root");
  HWW.push_back("histos2e2mu_25ns/output_HZJ_HToWWTo2L2Nu_ZTo2L_M125_13TeV_powheg_pythia8.root");
  TT4l.push_back("histos2e2mu_25ns/output_TTTo2L2Nu_noSC_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
  TTV.push_back("histos2e2mu_25ns/output_TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root");
  TTV.push_back("histos2e2mu_25ns/output_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  HWW.push_back("histos2e2mu_25ns/output_VBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8.root");
  qqH.push_back("histos2e2mu_25ns/output_VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  qqH3J.push_back("histos2e2mu_25ns/output_VBF_HJJJ_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  Wjets.push_back("histos2e2mu_25ns/output_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  VV.push_back("histos2e2mu_25ns/output_WWTo2L2Nu_13TeV-powheg.root");
  VVV.push_back("histos2e2mu_25ns/output_WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  VV.push_back("histos2e2mu_25ns/output_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root");
  VVV.push_back("histos2e2mu_25ns/output_WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  WH.push_back("histos2e2mu_25ns/output_WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root");
  WH.push_back("histos2e2mu_25ns/output_WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root");
  ZH.push_back("histos2e2mu_25ns/output_ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8.root");
  qqZZ.push_back("histos2e2mu_25ns/output_ZZJJTo4L_EWK_13TeV-madgraph-pythia8.root");
  qqZZ.push_back("histos2e2mu_25ns/output_ZZTo2L2Nu_13TeV_powheg_pythia8.root");
  qqZZ.push_back("histos2e2mu_25ns/output_ZZTo4L_13TeV_powheg_pythia8.root");
  VVV.push_back("histos2e2mu_25ns/output_ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  HWW.push_back("histos2e2mu_25ns/output_bbHToWWTo2L2Nu_M-125_4FS_yb2_13TeV_amcatnlo.root");
  HWW.push_back("histos2e2mu_25ns/output_bbHToWWTo2L2Nu_M-125_4FS_ybyt_13TeV_amcatnlo.root");
  ttH.push_back("histos2e2mu_25ns/output_ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUgenV6_pythia8.root");

  std::vector< std::vector<TString> > samples;
  std::vector<TString> tags;

  samples.push_back(data_obs);	tags.push_back("data_obs");
  samples.push_back(qqH);	tags.push_back("qqH");
  samples.push_back(ggH);	tags.push_back("ggH");
  samples.push_back(ggZZ);	tags.push_back("ggZZ");
  samples.push_back(qqZZ);	tags.push_back("qqZZ");
  samples.push_back(WH);	tags.push_back("WH");
  samples.push_back(ZH);	tags.push_back("ZH");
  samples.push_back(ttH);	tags.push_back("ttH");
  samples.push_back(HWW);	tags.push_back("HWW");
  samples.push_back(TTV);	tags.push_back("TTV");
  samples.push_back(VVV);	tags.push_back("VVV");
  samples.push_back(VV);	tags.push_back("VV");
  samples.push_back(Wjets);	tags.push_back("Wjets");
  samples.push_back(TT4l);	tags.push_back("TT4l");
  samples.push_back(qqH3J);	tags.push_back("qqH3J");

  TString tag = argv[1];
  std::cout<<">>>>>>>>>>>>> Looping over files for "<<tag<<" <<<<<<<<<<<<<<<<<"<<std::endl;
  int code = loopSamples(tags, samples, tag);
  std::cout<<">>>>>>>>>>>>> Finished <<<<<<<<<<<<<<<<<"<<std::endl;
  std::cout<<"CODE: "<<code<<std::endl;
     
  //Ends program
  return 0;
}
