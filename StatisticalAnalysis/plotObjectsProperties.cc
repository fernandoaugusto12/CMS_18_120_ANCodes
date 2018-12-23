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


//NN models
#include "/lustre/cms/store/user/mmelodea/Keras/modelsNN/k57nj2.h"
#include "/lustre/cms/store/user/mmelodea/Keras/modelsNN/k24nj3.h"


///Defining DR function with checking for deltaPhi < 2
#define pi_rad 3.141592654
float DR(float eta1, float eta2, float phi1, float phi2){
  float deltaPhi = fabs(phi1-phi2);
  if(deltaPhi > pi_rad)
    return sqrt(pow(eta1-eta2,2)+pow(2*pi_rad-deltaPhi,2));
  else 
    return sqrt(pow(eta1-eta2,2)+pow(phi1-phi2,2));
}

TString samples_path = "/lustrehome/mmelodea/MonoHiggs/MonoHiggs/80X/";


void loopSamples(std::vector<TString> tags, std::vector< std::vector<TString> > samples){
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
    
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

  
  
  std::vector<TString> histos_name;
  std::vector <int> nbins;
  std::vector <float> xmin, xmax;
  /*
  histos_name.push_back("lept1_pt");		nbins.push_back(25);	xmin.push_back(0);	xmax.push_back(250);
  histos_name.push_back("lept1_pt_error");	nbins.push_back(30);	xmin.push_back(0);	xmax.push_back(6);
  histos_name.push_back("lept1_eta");		nbins.push_back(20);	xmin.push_back(-2.5);	xmax.push_back(2.5);
  histos_name.push_back("lept1_phi");		nbins.push_back(20);	xmin.push_back(-3.2);	xmax.push_back(3.2);
  histos_name.push_back("lept1_pfx");		nbins.push_back(20);	xmin.push_back(0);	xmax.push_back(0.4);
  histos_name.push_back("lept1_sip");		nbins.push_back(25);	xmin.push_back(-5);	xmax.push_back(5);
  histos_name.push_back("lept2_pt");		nbins.push_back(25);	xmin.push_back(0);	xmax.push_back(250);
  histos_name.push_back("lept2_pt_error");	nbins.push_back(30);	xmin.push_back(0);	xmax.push_back(6);
  histos_name.push_back("lept2_eta");		nbins.push_back(20);	xmin.push_back(-2.5);	xmax.push_back(2.5);
  histos_name.push_back("lept2_phi");		nbins.push_back(20);	xmin.push_back(-3.2);	xmax.push_back(3.2);
  histos_name.push_back("lept2_pfx");		nbins.push_back(20);	xmin.push_back(0);	xmax.push_back(0.4);
  histos_name.push_back("lept2_sip");		nbins.push_back(25);	xmin.push_back(-5);	xmax.push_back(5);
  histos_name.push_back("lept3_pt");		nbins.push_back(25);	xmin.push_back(0);	xmax.push_back(250);
  histos_name.push_back("lept3_pt_error");	nbins.push_back(30);	xmin.push_back(0);	xmax.push_back(6);
  histos_name.push_back("lept3_eta");		nbins.push_back(20);	xmin.push_back(-2.5);	xmax.push_back(2.5);
  histos_name.push_back("lept3_phi");		nbins.push_back(20);	xmin.push_back(-3.2);	xmax.push_back(3.2);
  histos_name.push_back("lept3_pfx");		nbins.push_back(20);	xmin.push_back(0);	xmax.push_back(0.4);
  histos_name.push_back("lept3_sip");		nbins.push_back(25);	xmin.push_back(-5);	xmax.push_back(5);
  histos_name.push_back("lept4_pt");		nbins.push_back(25);	xmin.push_back(0);	xmax.push_back(250);
  histos_name.push_back("lept4_pt_error");	nbins.push_back(30);	xmin.push_back(0);	xmax.push_back(6);
  histos_name.push_back("lept4_eta");		nbins.push_back(20);	xmin.push_back(-2.5);	xmax.push_back(2.5);
  histos_name.push_back("lept4_phi");		nbins.push_back(20);	xmin.push_back(-3.2);	xmax.push_back(3.2);
  histos_name.push_back("lept4_pfx");		nbins.push_back(20);	xmin.push_back(0);	xmax.push_back(0.4);
  histos_name.push_back("lept4_sip");		nbins.push_back(25);	xmin.push_back(-5);	xmax.push_back(5);
  histos_name.push_back("Z1mass");		nbins.push_back(36);	xmin.push_back(12);	xmax.push_back(120);
  histos_name.push_back("Z2mass");		nbins.push_back(36);	xmin.push_back(12);	xmax.push_back(120);
  histos_name.push_back("costhetastar");	nbins.push_back(15);	xmin.push_back(-1.05);	xmax.push_back(1.05);
  histos_name.push_back("costheta1");		nbins.push_back(15);	xmin.push_back(-1.05);	xmax.push_back(1.05);
  histos_name.push_back("costheta2");		nbins.push_back(15);	xmin.push_back(-1.05);	xmax.push_back(1.05);
  histos_name.push_back("phi");			nbins.push_back(20);	xmin.push_back(-3.16);	xmax.push_back(3.16);
  histos_name.push_back("phistar1");		nbins.push_back(20);	xmin.push_back(-3.16);	xmax.push_back(3.16);
  histos_name.push_back("pt4l");		nbins.push_back(25);	xmin.push_back(0);	xmax.push_back(250);
  histos_name.push_back("eta4l");		nbins.push_back(25);	xmin.push_back(-5);	xmax.push_back(5);
  */
  histos_name.push_back("mass4l_smhiggs");	nbins.push_back(1000);	xmin.push_back(0);	xmax.push_back(1000);
  histos_name.push_back("mass4l_vbf");		nbins.push_back(14);	xmin.push_back(117);	xmax.push_back(131);
  //histos_name.push_back("mass4l_zoom");		nbins.push_back(4);	xmin.push_back(116);	xmax.push_back(132);
  //histos_name.push_back("njets_pass");		nbins.push_back(10);	xmin.push_back(0);	xmax.push_back(10);
  //histos_name.push_back("deltajj");		nbins.push_back(47);	xmin.push_back(0);	xmax.push_back(9.4);
  /*
  histos_name.push_back("massjj");		nbins.push_back(5);	xmin.push_back(0);	xmax.push_back(2000);
  histos_name.push_back("Djet_VAJHU");		nbins.push_back(10);	xmin.push_back(-0.05);	xmax.push_back(1.05);
  histos_name.push_back("Djet_VAJHU_error");	nbins.push_back(50);	xmin.push_back(-10);	xmax.push_back(10);
  histos_name.push_back("met");			nbins.push_back(10);	xmin.push_back(0);	xmax.push_back(200);
  histos_name.push_back("met_eta");		nbins.push_back(30);	xmin.push_back(-6);	xmax.push_back(6);
  histos_name.push_back("met_phi");		nbins.push_back(20);	xmin.push_back(-3.2);	xmax.push_back(3.2);
  histos_name.push_back("met_jetEnDn");		nbins.push_back(60);	xmin.push_back(0);	xmax.push_back(20);
  histos_name.push_back("met_jetEnUp");		nbins.push_back(60);	xmin.push_back(0);	xmax.push_back(20);
  histos_name.push_back("met_jetResDn");	nbins.push_back(60);	xmin.push_back(0);	xmax.push_back(20);
  histos_name.push_back("met_jetResUp");	nbins.push_back(60);	xmin.push_back(0);	xmax.push_back(20);
  histos_name.push_back("met_eleEnDn");		nbins.push_back(60);	xmin.push_back(0);	xmax.push_back(20);
  histos_name.push_back("met_eleEnUp");		nbins.push_back(60);	xmin.push_back(0);	xmax.push_back(20);
  histos_name.push_back("met_muEnDn");		nbins.push_back(60);	xmin.push_back(0);	xmax.push_back(20);
  histos_name.push_back("met_muEnUp");		nbins.push_back(60);	xmin.push_back(0);	xmax.push_back(20);
  histos_name.push_back("met_photEnDn");	nbins.push_back(60);	xmin.push_back(0);	xmax.push_back(20);
  histos_name.push_back("met_photEnUp");	nbins.push_back(60);	xmin.push_back(0);	xmax.push_back(20);
  histos_name.push_back("met_uncEnDn");		nbins.push_back(60);	xmin.push_back(0);	xmax.push_back(20);
  histos_name.push_back("met_uncEnUp");		nbins.push_back(60);	xmin.push_back(0);	xmax.push_back(20);
  histos_name.push_back("nbjets");		nbins.push_back(10);	xmin.push_back(0);	xmax.push_back(10);
  histos_name.push_back("dphi");		nbins.push_back(25);	xmin.push_back(-5);	xmax.push_back(5);
  histos_name.push_back("mT");			nbins.push_back(25);	xmin.push_back(0);	xmax.push_back(250);
  histos_name.push_back("lhe_npartons");	nbins.push_back(10);	xmin.push_back(0);	xmax.push_back(10);
  histos_name.push_back("lhe_parton_clear");	nbins.push_back(2);	xmin.push_back(0);	xmax.push_back(2);
  histos_name.push_back("lhe_parton_pdgid");	nbins.push_back(22);	xmin.push_back(0);	xmax.push_back(22);
  histos_name.push_back("lhe_parton_pt");	nbins.push_back(25);	xmin.push_back(0);	xmax.push_back(100);
  histos_name.push_back("lhe_parton_eta");	nbins.push_back(25);	xmin.push_back(-5);	xmax.push_back(5);
  histos_name.push_back("lhe_parton_phi");	nbins.push_back(20);	xmin.push_back(-3.2);	xmax.push_back(3.2);
  histos_name.push_back("lhe_parton_e");	nbins.push_back(25);	xmin.push_back(0);	xmax.push_back(100);
  //histos_name.push_back("weight");		nbins.push_back(10000);	xmin.push_back(0);	xmax.push_back(1);
  
  unsigned int njh = 0;
  for(unsigned int ijet=0; ijet<maxJets; ++ijet){
    histos_name.push_back(Form("jets_highpt_btagger[%i]",ijet));			nbins.push_back(50);	xmin.push_back(0);	xmax.push_back(1);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_pt[%i]",ijet));				nbins.push_back(10);	xmin.push_back(0);	xmax.push_back(200);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_pt_error[%i]",ijet));			nbins.push_back(1000);	xmin.push_back(0);	xmax.push_back(10);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_eta[%i]",ijet));				nbins.push_back(25);	xmin.push_back(-5);	xmax.push_back(5);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_phi[%i]",ijet));				nbins.push_back(20);	xmin.push_back(-3.2);	xmax.push_back(3.2);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_et[%i]",ijet));				nbins.push_back(25);	xmin.push_back(0);	xmax.push_back(300);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_area[%i]",ijet));				nbins.push_back(100);	xmin.push_back(0);	xmax.push_back(1);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_ptd[%i]",ijet));				nbins.push_back(50);	xmin.push_back(-0.05);	xmax.push_back(1.05);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_charged_hadron_energy[%i]",ijet));		nbins.push_back(50);	xmin.push_back(0);	xmax.push_back(500);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_neutral_hadron_energy[%i]",ijet));		nbins.push_back(50);	xmin.push_back(0);	xmax.push_back(500);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_photon_energy[%i]",ijet));			nbins.push_back(50);	xmin.push_back(0);	xmax.push_back(500);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_electron_energy[%i]",ijet));		nbins.push_back(50);	xmin.push_back(0);	xmax.push_back(500);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_muon_energy[%i]",ijet));			nbins.push_back(25);	xmin.push_back(0);	xmax.push_back(100);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_hf_hadron_energy[%i]",ijet));		nbins.push_back(50);	xmin.push_back(0);	xmax.push_back(500);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_hf_em_energy[%i]",ijet));			nbins.push_back(50);	xmin.push_back(0);	xmax.push_back(150);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_charged_em_energy[%i]",ijet));		nbins.push_back(50);	xmin.push_back(0);	xmax.push_back(500);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_charged_mu_energy[%i]",ijet));		nbins.push_back(25);	xmin.push_back(0);	xmax.push_back(100);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_neutral_em_energy[%i]",ijet));		nbins.push_back(50);	xmin.push_back(0);	xmax.push_back(500);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_charged_hadron_multiplicity[%i]",ijet));	nbins.push_back(20);	xmin.push_back(0);	xmax.push_back(20);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_neutral_hadron_multiplicity[%i]",ijet));	nbins.push_back(10);	xmin.push_back(0);	xmax.push_back(10);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_photon_multiplicity[%i]",ijet));		nbins.push_back(15);	xmin.push_back(0);	xmax.push_back(15);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_electron_multiplicity[%i]",ijet));		nbins.push_back(10);	xmin.push_back(0);	xmax.push_back(10);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_muon_multiplicity[%i]",ijet));		nbins.push_back(10);	xmin.push_back(0);	xmax.push_back(10);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_hf_hadron_multiplicity[%i]",ijet));		nbins.push_back(20);	xmin.push_back(0);	xmax.push_back(20);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_hf_em_multiplicity[%i]",ijet));		nbins.push_back(10);	xmin.push_back(0);	xmax.push_back(10);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_charged_multiplicity[%i]",ijet));		nbins.push_back(20);	xmin.push_back(0);	xmax.push_back(20);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_neutral_multiplicity[%i]",ijet));		nbins.push_back(20);	xmin.push_back(0);	xmax.push_back(20);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_ncomponents[%i]",ijet));			nbins.push_back(100);	xmin.push_back(0);	xmax.push_back(100);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_component_pdgid[%i]",ijet));		nbins.push_back(250);	xmin.push_back(0);	xmax.push_back(250);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_component_pt[%i]",ijet));			nbins.push_back(25);	xmin.push_back(0);	xmax.push_back(100);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_component_eta[%i]",ijet));			nbins.push_back(25);	xmin.push_back(-5);	xmax.push_back(5);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_component_phi[%i]",ijet));			nbins.push_back(20);	xmin.push_back(-3.2);	xmax.push_back(3.2);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_component_energy[%i]",ijet));		nbins.push_back(50);	xmin.push_back(0);	xmax.push_back(250);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_component_charge[%i]",ijet));		nbins.push_back(5);	xmin.push_back(-2);	xmax.push_back(3);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_component_mt[%i]",ijet));			nbins.push_back(25);	xmin.push_back(0);	xmax.push_back(50);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_component_xvertex[%i]",ijet));		nbins.push_back(40);	xmin.push_back(-20);	xmax.push_back(20);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_component_yvertex[%i]",ijet));		nbins.push_back(40);	xmin.push_back(-20);	xmax.push_back(20);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_component_zvertex[%i]",ijet));		nbins.push_back(60);	xmin.push_back(-30);	xmax.push_back(30);	if(ijet==0) ++njh;
    histos_name.push_back(Form("jets_highpt_component_vertex_chi2[%i]",ijet));		nbins.push_back(100);	xmin.push_back(0);	xmax.push_back(0.05);	if(ijet==0) ++njh;
  }
  */
  const unsigned int nhists = histos_name.size();
  std::cout<<"Generating "<<nhists<<" plots......."<<std::endl;
  unsigned int count_hists = 0;
  bool data_isPresent = true;

  
  std::vector<std::vector<double> > sig_inputs;
  std::vector<std::vector<double> > bkg_inputs;
  float dnn_serr = 0, dnn_berr = 0, ansig = 0, anbkg = 0;
  float fsig_err = 0, fbkg_err = 0;
  float fnsig = 0, fnbkg = 0, fndata = 0;
  float fs_nevents = 0, fs_yields = 0, fs_yields_err = 0;
  float sig_4mu = 0, sig_4mu_err = 0, sig_4e = 0, sig_4e_err = 0, sig_2e2mu = 0, sig_2e2mu_err = 0;
  float bkg_4mu = 0, bkg_4mu_err = 0, bkg_4e = 0, bkg_4e_err = 0, bkg_2e2mu = 0, bkg_2e2mu_err = 0;
  const unsigned int ntags = tags.size();
  TH1D *histos[ntags][4][nhists];

  TH1D *k57nj2[ntags][4];
  TH1D *k24nj3[ntags][4];
  TH2D *k57nj2_deltajj[ntags][4];
  TH2D *k24nj3_deltajj[ntags][4];

  for(unsigned int itag=0; itag<ntags; ++itag){
    
    for(unsigned int ifs=0; ifs<4; ++ifs){
      k57nj2[itag][ifs] = new TH1D("","",120,0,1);
      k24nj3[itag][ifs] = new TH1D("","",120,0,1);
      
      k57nj2_deltajj[itag][ifs] = new TH2D("","",12,0,1,10,0,9.4);
      k24nj3_deltajj[itag][ifs] = new TH2D("","",12,0,1,10,0,9.4);
    }
    
    if(tags[itag] == "Data") data_isPresent = true;
    float events = 0, yields = 0, yields_err = 0;
    for(unsigned int yhist=0; yhist<nhists; ++yhist)
      for(unsigned int ifs=0; ifs<4; ++ifs) 
	histos[itag][ifs][yhist] = new TH1D("","",nbins[yhist],xmin[yhist],xmax[yhist]);
    unsigned int nsamples = samples[itag].size();
    for(unsigned int is=0; is<nsamples; ++is){
      TString ifile_name_root = samples_path+samples[itag][is];
      if(gSystem->AccessPathName(ifile_name_root)){
	cout<<"File "<<ifile_name_root<<" doesn't exist!"<<endl;
	continue;
      }
      //else cout<<"Loading file "<<ifile_name_root<<endl;
  
      TFile *ofile = TFile::Open(ifile_name_root);
      TTree *otree = (TTree*)ofile->Get("HZZ4LeptonsAnalysisReduced");
      //otree->SetBranchAddress("f_outlier", &f_outlier);
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

      int pass_events = 0;
      for(int ievent=0; ievent<Nentries; ievent++){
	otree->GetEntry(ievent);
	//if(ievent > 10) break;
        
	  float l1pt   = f_lept1_pt;
	  float l1eta  = f_lept1_eta;
	  float l1phi  = f_lept1_phi;
	  float l2pt   = f_lept2_pt;
	  float l2eta  = f_lept2_eta;
	  float l2phi  = f_lept2_phi;
	  float l3pt   = f_lept3_pt;
	  float l3eta  = f_lept3_eta;
	  float l3phi  = f_lept3_phi;
	  float l4pt   = f_lept4_pt;
	  float l4eta  = f_lept4_eta;
	  float l4phi  = f_lept4_phi;
	  float j1pt   = f_jets_highpt_pt[0];
	  float j1eta  = f_jets_highpt_eta[0];
	  float j1phi  = f_jets_highpt_phi[0];
	  float j2pt   = f_jets_highpt_pt[1];
	  float j2eta  = f_jets_highpt_eta[1];
	  float j2phi  = f_jets_highpt_phi[1];
	  float j3pt   = (f_jets_highpt_pt[2] != -999)? f_jets_highpt_pt[2]: 0;
	  float j3eta  = (f_jets_highpt_pt[2] != -999)? f_jets_highpt_eta[2]: 0;
	  float j3phi  = (f_jets_highpt_pt[2] != -999)? f_jets_highpt_phi[2]: 0;
	  float met    = f_pfmet;
	  float njets  = f_njets_pass;
	  float nbjets = f_Nbjets;
	  
	  std::vector<std::vector<double> > inputs = {
	    //k57nj2
	    {l1pt,l1eta,l1phi,l2pt,l2eta,l2phi,l3pt,l3eta,l3phi,l4pt,l4eta,l4phi,j1pt,j1eta,j1phi,j2pt,j2eta,j2phi,met},
	    //k24nj3
	    {l1pt,l1eta,l1phi,l2pt,l2eta,l2phi,l3pt,l3eta,l3phi,l4pt,l4eta,l4phi,j1pt,j1eta,j1phi,j2pt,j2eta,j2phi,j3pt,j3eta,j3phi,met}
	  };
	  
	  double vk57nj2 = model_k57nj2(inputs[0]);
	  double vk24nj3 = model_k24nj3(inputs[1]);
	  
	  if( (f_njets_pass == 2 && vk57nj2 <= 0.85) || (f_njets_pass >= 3 && vk24nj3 <= 0.73) ) continue;
	
	int nmu=0, nele=0;      
	if(fabs(f_lept1_pdgid)==11) ++nele;
	else ++nmu;
	if(fabs(f_lept2_pdgid)==11) ++nele;
	else ++nmu;
	if(fabs(f_lept3_pdgid)==11) ++nele;
	else ++nmu;
	if(fabs(f_lept4_pdgid)==11) ++nele;
	else ++nmu;
	  
	//if(Nentries > 10 && ievent % (Nentries/10) == 0)
	  //cout<<"iEvent -------------- "<<ievent<<endl;
	histos[itag][0][0]->Fill(f_mass4l, f_weight);
	if(nmu == 4 && nele == 0) histos[itag][1][0]->Fill(f_mass4l, f_weight);
	else if(nele == 4 && nmu == 0) histos[itag][2][0]->Fill(f_mass4l, f_weight);
	else if(nmu == 2 && nele == 2) histos[itag][3][0]->Fill(f_mass4l, f_weight);
	else{
	  std::cout<<"!!Something is bad! Different number of letpons than expected!!"<<std::endl;
	}
	
	if( (((f_njets_pass == 2 || f_njets_pass == 3) && f_Nbjets <= 1) || (f_njets_pass > 3 && f_Nbjets == 0)) && f_mass4l >= 118 && f_mass4l <= 130 ){
	  histos[itag][0][1]->Fill(f_mass4l, f_weight);
	  if(nmu == 4 && nele == 0) histos[itag][1][1]->Fill(f_mass4l, f_weight);
	  else if(nele == 4 && nmu == 0) histos[itag][2][1]->Fill(f_mass4l, f_weight);
	  else if(nmu == 2 && nele == 2) histos[itag][3][1]->Fill(f_mass4l, f_weight);
	  else{
	    std::cout<<"!!Something is bad! Different number of letpons than expected!!"<<std::endl;
	  }
	  	  
	  if(f_njets_pass == 2){
	    k57nj2[itag][0]->Fill( vk57nj2, f_weight );
	    k57nj2_deltajj[itag][0]->Fill( vk57nj2, fabs(j1eta-j2eta), f_weight );
	    if(nmu == 4){
	      k57nj2[itag][1]->Fill( vk57nj2, f_weight );
	      k57nj2_deltajj[itag][1]->Fill( vk57nj2, fabs(j1eta-j2eta), f_weight );
	    }
	    if(nele == 4){
	      k57nj2[itag][2]->Fill( vk57nj2, f_weight );
	      k57nj2_deltajj[itag][2]->Fill( vk57nj2, fabs(j1eta-j2eta), f_weight );
	    }
	    if(nmu == 2 && nele == 2){
	      k57nj2[itag][3]->Fill( vk57nj2, f_weight );
	      k57nj2_deltajj[itag][3]->Fill( vk57nj2, fabs(j1eta-j2eta), f_weight );
	    }
	  }
	  if(f_njets_pass >= 3){
	    k24nj3[itag][0]->Fill( vk24nj3, f_weight );
	    k24nj3_deltajj[itag][0]->Fill( vk24nj3, fabs(j1eta-j2eta), f_weight );
	    if(nmu == 4){
	      k24nj3[itag][1]->Fill( vk24nj3, f_weight );
	      k24nj3_deltajj[itag][1]->Fill( vk24nj3, fabs(j1eta-j2eta), f_weight );
	    }
	    if(nele == 4){
	      k24nj3[itag][2]->Fill( vk24nj3, f_weight );
	      k24nj3_deltajj[itag][2]->Fill( vk24nj3, fabs(j1eta-j2eta), f_weight );
	    }
	    if(nmu == 2 && nele == 2){
	      k24nj3[itag][3]->Fill( vk24nj3, f_weight );
	      k24nj3_deltajj[itag][3]->Fill( vk24nj3, fabs(j1eta-j2eta), f_weight );
	    }
	  }
	//}
	
        //{
	  //cout<<"Pass events: "<<pass_events<<endl;
	  
	  ++events;
	  yields += f_weight;
	  yields_err += f_weight*f_weight;
	  ++pass_events;
	  if(tags[itag] == "Data"){
	    fndata += f_weight;
	  }
	  if(tags[itag] == "qqH"){
	    fnsig += f_weight;
	    fsig_err += f_weight*f_weight;
	    if(nmu == 4){sig_4mu += f_weight; sig_4mu_err += f_weight*f_weight;}
	    if(nele == 4){sig_4e += f_weight; sig_4e_err += f_weight*f_weight;}
	    if(nmu == 2 && nele == 2){sig_2e2mu += f_weight; sig_2e2mu_err += f_weight*f_weight;}
	  }
	  if(tags[itag] != "Data" && tags[itag] != "qqH" && tags[itag] != "qqH3J"){
	    fnbkg += f_weight;
	    fbkg_err += f_weight*f_weight;
	    if(nmu == 4){bkg_4mu += f_weight; bkg_4mu_err += f_weight*f_weight;}
	    if(nele == 4){bkg_4e += f_weight; bkg_4e_err += f_weight*f_weight;}
	    if(nmu == 2 && nele == 2){bkg_2e2mu += f_weight; bkg_2e2mu_err += f_weight*f_weight;}
	  }
	  //Numbers per final state
	  ++fs_nevents;
	  fs_yields += f_weight;
	  fs_yields_err += f_weight*f_weight;
	  
	  //if(ievent > 100) break;
	  //if(f_outlier) continue;

	  /*
	  histos[itag][0]->Fill(f_lept1_pt, f_weight);
	  if(f_lept1_pt != 0)
	    histos[itag][1]->Fill(f_lept1_pt_error, f_weight);
	  histos[itag][2]->Fill(f_lept1_eta, f_weight);
	  histos[itag][3]->Fill(f_lept1_phi, f_weight);
	  histos[itag][4]->Fill(f_lept1_pfx, f_weight);
	  histos[itag][5]->Fill(f_lept1_sip, f_weight);
	  histos[itag][6]->Fill(f_lept2_pt, f_weight);
	  if(f_lept2_pt != 0)
	    histos[itag][7]->Fill(f_lept2_pt_error, f_weight);
	  histos[itag][8]->Fill(f_lept2_eta, f_weight);
	  histos[itag][9]->Fill(f_lept2_phi, f_weight);
	  histos[itag][10]->Fill(f_lept2_pfx, f_weight);
	  histos[itag][11]->Fill(f_lept2_sip, f_weight);
	  histos[itag][12]->Fill(f_lept3_pt, f_weight);
	  if(f_lept3_pt != 0)
	    histos[itag][13]->Fill(f_lept3_pt_error, f_weight);
	  histos[itag][14]->Fill(f_lept3_eta, f_weight);
	  histos[itag][15]->Fill(f_lept3_phi, f_weight);
	  histos[itag][16]->Fill(f_lept3_pfx, f_weight);
	  histos[itag][17]->Fill(f_lept3_sip, f_weight);
	  histos[itag][18]->Fill(f_lept4_pt, f_weight);
	  if(f_lept4_pt != 0)
	    histos[itag][19]->Fill(f_lept4_pt_error, f_weight);
	  histos[itag][20]->Fill(f_lept4_eta, f_weight);
	  histos[itag][21]->Fill(f_lept4_phi, f_weight);
	  histos[itag][22]->Fill(f_lept4_pfx, f_weight);
	  histos[itag][23]->Fill(f_lept4_sip, f_weight);
	  histos[itag][24]->Fill(f_Z1mass, f_weight);
	  histos[itag][25]->Fill(f_Z2mass, f_weight);
	  histos[itag][26]->Fill(f_angle_costhetastar, f_weight);
	  histos[itag][27]->Fill(f_angle_costheta1, f_weight);
	  histos[itag][28]->Fill(f_angle_costheta2, f_weight);
	  histos[itag][29]->Fill(f_angle_phi, f_weight);
	  histos[itag][30]->Fill(f_angle_phistar1, f_weight);
	  histos[itag][31]->Fill(f_pt4l, f_weight);
	  histos[itag][32]->Fill(f_eta4l, f_weight);
	  if(f_mass4l >= 118 && f_mass4l <=130){
	    histos[itag][34]->Fill(f_mass4l, f_weight);
	  }
	  histos[itag][35]->Fill(f_njets_pass, f_weight);
	  histos[itag][36]->Fill(f_deltajj, f_weight);
	  histos[itag][37]->Fill(f_massjj, f_weight);
	  histos[itag][38]->Fill(f_Djet_VAJHU, f_weight);
	  if(f_Djet_VAJHU != 0)
	    histos[itag][39]->Fill(fabs(f_Djet_VAJHU_UncDn-f_Djet_VAJHU), f_weight);

	  histos[itag][40]->Fill(f_pfmet, f_weight);
	  histos[itag][41]->Fill(-TMath::Log(TMath::Tan(f_pfmet_theta/2.)), f_weight);
	  histos[itag][42]->Fill(f_pfmet_phi, f_weight);
	  histos[itag][43]->Fill(f_pfmet_JetEnDn, f_weight);
	  histos[itag][44]->Fill(f_pfmet_JetEnUp, f_weight);
	  histos[itag][45]->Fill(f_pfmet_JetResDn, f_weight);
	  histos[itag][46]->Fill(f_pfmet_JetResUp, f_weight);
	  histos[itag][47]->Fill(f_pfmet_ElectronEnDn, f_weight);
	  histos[itag][48]->Fill(f_pfmet_ElectronEnUp, f_weight);
	  histos[itag][49]->Fill(f_pfmet_MuonEnDn, f_weight);
	  histos[itag][50]->Fill(f_pfmet_MuonEnUp, f_weight);
	  histos[itag][51]->Fill(f_pfmet_PhotonEnDn, f_weight);
	  histos[itag][52]->Fill(f_pfmet_PhotonEnUp, f_weight);
	  histos[itag][53]->Fill(f_pfmet_UnclusteredEnDn, f_weight);
	  histos[itag][54]->Fill(f_pfmet_UnclusteredEnUp, f_weight);

	  histos[itag][55]->Fill(f_Nbjets, f_weight);
	  histos[itag][56]->Fill(f_dphi, f_weight);
	  histos[itag][57]->Fill(f_mT, f_weight);
	  
	  histos[itag][58]->Fill(f_lhe_npartons, f_weight);
	  for(unsigned int iparton=0; iparton<maxPartons; ++iparton){
	    if(iparton < (unsigned int)f_lhe_npartons){
	      histos[itag][59]->Fill(f_lhe_parton_clear[iparton], f_weight);
	      histos[itag][60]->Fill(f_lhe_parton_pdgid[iparton], f_weight);
	      histos[itag][61]->Fill(f_lhe_parton_pt[iparton], f_weight);
	      histos[itag][62]->Fill(f_lhe_parton_eta[iparton], f_weight);
	      histos[itag][63]->Fill(f_lhe_parton_phi[iparton], f_weight);
	      histos[itag][64]->Fill(f_lhe_parton_e[iparton], f_weight);
	    }
	  }
	  //histos[itag][65]->Fill(f_weight, 1);
	  
	  unsigned int hindex = 64;
	  for(unsigned int jjet=0; jjet<maxJets; ++jjet){
	    unsigned int jump_factor = hindex;//njh*jjet + hindex;
	    if(jjet < f_njets_pass){
	      histos[itag][1+jump_factor]->Fill(f_jets_highpt_btagger[jjet], f_weight);
	      histos[itag][2+jump_factor]->Fill(f_jets_highpt_pt[jjet], f_weight);
	      if(f_jets_highpt_pt[jjet] != 0)
		histos[itag][3+jump_factor]->Fill(+f_jets_highpt_pt_error[jjet], f_weight);
	      histos[itag][4+jump_factor]->Fill(f_jets_highpt_eta[jjet], f_weight);
	      histos[itag][5+jump_factor]->Fill(f_jets_highpt_phi[jjet], f_weight);
	      histos[itag][6+jump_factor]->Fill(f_jets_highpt_et[jjet], f_weight);
	      histos[itag][7+jump_factor]->Fill(f_jets_highpt_area[jjet], f_weight);
	      histos[itag][8+jump_factor]->Fill(f_jets_highpt_ptd[jjet], f_weight);
	      histos[itag][9+jump_factor]->Fill(f_jets_highpt_charged_hadron_energy[jjet], f_weight);
	      histos[itag][10+jump_factor]->Fill(f_jets_highpt_neutral_hadron_energy[jjet], f_weight);
	      histos[itag][11+jump_factor]->Fill(f_jets_highpt_photon_energy[jjet], f_weight);
	      histos[itag][12+jump_factor]->Fill(f_jets_highpt_electron_energy[jjet], f_weight);
	      histos[itag][13+jump_factor]->Fill(f_jets_highpt_muon_energy[jjet], f_weight);
	      histos[itag][14+jump_factor]->Fill(f_jets_highpt_hf_hadron_energy[jjet], f_weight);
	      histos[itag][15+jump_factor]->Fill(f_jets_highpt_hf_em_energy[jjet], f_weight);
	      histos[itag][16+jump_factor]->Fill(f_jets_highpt_charged_em_energy[jjet], f_weight);
	      histos[itag][17+jump_factor]->Fill(f_jets_highpt_charged_mu_energy[jjet], f_weight);
	      histos[itag][18+jump_factor]->Fill(f_jets_highpt_neutral_em_energy[jjet], f_weight);
	      histos[itag][19+jump_factor]->Fill(f_jets_highpt_charged_hadron_multiplicity[jjet], f_weight);
	      histos[itag][20+jump_factor]->Fill(f_jets_highpt_neutral_hadron_multiplicity[jjet], f_weight);
	      histos[itag][21+jump_factor]->Fill(f_jets_highpt_photon_multiplicity[jjet], f_weight);
	      histos[itag][22+jump_factor]->Fill(f_jets_highpt_electron_multiplicity[jjet], f_weight);
	      histos[itag][23+jump_factor]->Fill(f_jets_highpt_muon_multiplicity[jjet], f_weight);
	      histos[itag][24+jump_factor]->Fill(f_jets_highpt_hf_hadron_multiplicity[jjet], f_weight);
	      histos[itag][25+jump_factor]->Fill(f_jets_highpt_hf_em_multiplicity[jjet], f_weight);
	      histos[itag][26+jump_factor]->Fill(f_jets_highpt_charged_multiplicity[jjet], f_weight);
	      histos[itag][27+jump_factor]->Fill(f_jets_highpt_neutral_multiplicity[jjet], f_weight);
	      histos[itag][28+jump_factor]->Fill(f_jets_highpt_ncomponents[jjet], f_weight);
	      
	      for(int jsjet=0; jsjet<maxJetsComponents; ++jsjet){
		if(jsjet < f_jets_highpt_ncomponents[jjet]){
		  histos[itag][29+jump_factor]->Fill(fabs(f_jets_highpt_component_pdgid[jjet][jsjet]), f_weight);
		  histos[itag][30+jump_factor]->Fill(f_jets_highpt_component_pt[jjet][jsjet], f_weight);
		  histos[itag][31+jump_factor]->Fill(f_jets_highpt_component_eta[jjet][jsjet], f_weight);
		  histos[itag][32+jump_factor]->Fill(f_jets_highpt_component_phi[jjet][jsjet], f_weight);
		  histos[itag][33+jump_factor]->Fill(f_jets_highpt_component_energy[jjet][jsjet], f_weight);
		  histos[itag][34+jump_factor]->Fill(f_jets_highpt_component_charge[jjet][jsjet], f_weight);
		  histos[itag][35+jump_factor]->Fill(f_jets_highpt_component_mt[jjet][jsjet], f_weight);
		  histos[itag][36+jump_factor]->Fill(f_jets_highpt_component_xvertex[jjet][jsjet], f_weight);
		  histos[itag][37+jump_factor]->Fill(f_jets_highpt_component_yvertex[jjet][jsjet], f_weight);
		  histos[itag][38+jump_factor]->Fill(f_jets_highpt_component_zvertex[jjet][jsjet], f_weight);
		  histos[itag][39+jump_factor]->Fill(f_jets_highpt_component_vertex_chi2[jjet][jsjet], f_weight);
		}
	      }
	    }
	  }
	 */
	}//End of "if" for VBF selection
      }            
      ofile->Close();
      //cout<<ifile_name_root<<" = "<<pass_events<<endl;
      if( (is+1)%(nsamples/3) == 0){
	TString fstate = "";
	if( (is+1)/(nsamples/3) == 1 ) fstate = "4mu";
	if( (is+1)/(nsamples/3) == 2 ) fstate = "4e";
	if( (is+1)/(nsamples/3) == 3 ) fstate = "2e2mu";
	std::cout<<tags[itag]<<" -- "<<fstate<<Form(" -- events= %.0f -- yields= %.2f +/- %.2f",fs_nevents,fs_yields,sqrt(fs_yields_err))<<std::endl;
	fs_nevents = 0;
	fs_yields = 0;
	fs_yields_err = 0;
      }
    }//end samples loop
    std::cout<<tags[itag]<<Form(" -- 4l -- events= %.0f -- yields= %.2f +/- %.2f",events,yields,sqrt(yields_err))<<std::endl;
    std::cout<<"-----------------------------------------------------------"<<std::endl;
  }//end tags loop
    
  gROOT->SetBatch();
  gStyle->SetOptStat(0);
    
  
  //Gets m4l from Z+X
  //TFile *fzx = TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/zxstatisticsOS/smhiggs_4l_nominal.root");
  //TH1D *hzx = (TH1D*)fzx->Get("m4lZX");
  
  /*
  ////Plot to do
  //std::cout<<"Creating CMS tag...."<<std::endl;
  TPaveText *cms_tag = new TPaveText(.4,.95,.99,.96,"NDC");
  cms_tag->AddText("CMS #bf{Preliminary  #sqrt{s} = 13 TeV, L = 35.9fb^{-1}}");
  cms_tag->SetFillStyle(0);
  cms_tag->SetBorderSize(0);
  cms_tag->SetTextSize(0.04);
  
  //std::cout<<"Creating TLegend...."<<std::endl;
  TLegend *leg = new TLegend(0.4,0.7,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetEntrySeparation(0.1);
  leg->SetNColumns(3);
  
  TPaveText *count_tag = new TPaveText(.2,.6,4,.85,"NDC");
  if(data_isPresent) count_tag->AddText(Form("Data: %.0f",fndata));
  count_tag->AddText(Form("S: %.2f#pm%.2f",fnsig,sqrt(fsig_err)));
  count_tag->AddText(Form("B: %.2f#pm%.2f",fnbkg,sqrt(fbkg_err)));
  count_tag->SetFillStyle(0);
  count_tag->SetBorderSize(0);
  count_tag->SetTextSize(0.05);  

  //std::cout<<"Creating TCanvas...."<<std::endl;
  TCanvas *cv[nhists];
  TPad *pad1[nhists];
  TPad *pad2[nhists];

  //std::cout<<"Stacking histograms and computing statistics...."<<std::endl;
  int vbf3j_tag = 0, vbf2j_tag = 0;
  THStack *shistos[nhists];
  TGraph *hdata[nhists];
  TGraphAsymmErrors *fdata[nhists];
  TH1D *hratio[nhists];
  TH1D *mc_sum[nhists];
  TF1 *f1 = new TF1("f1","[0]",-10000,10000);
  TFile *histograms_plots = new TFile("Histograms.root","recreate");
  for(unsigned int ihist=0; ihist<nhists; ++ihist){
    //std::cout<<"Creating canvas..."<<std::endl;
    cv[ihist] = new TCanvas(histos_name[ihist],"",10,10,700,700);
    pad1[ihist] = new TPad("","",0.05,0.05,0.95,0.30);
    pad2[ihist] = new TPad("","",0.05,0.30,0.95,0.97);
    pad1[ihist]->SetTopMargin(0);
    pad1[ihist]->SetBottomMargin(0.3);
    pad2[ihist]->SetBottomMargin(0.02);
    pad1[ihist]->Draw();
    pad2[ihist]->Draw();
    
    shistos[ihist] = new THStack();
    hdata[ihist] = new TGraph();
    fdata[ihist] = new TGraphAsymmErrors();
    //std::cout<<"Going to get stack and data info..."<<std::endl;
    
    //Adds Z+X
    if(histos_name[ihist] == "mass4l"){
      hzx->Rebin(4);
      hzx->SetLineColor(kSpring-6);
      hzx->SetFillColor(kSpring-6);
      shistos[ihist]->Add( hzx );
    }
    
    float maxy = 0, maxx = 0;
    for(unsigned int itag=0; itag<ntags; ++itag){      
      if(tags[itag] == "qqH") vbf2j_tag = itag;
      if(tags[itag] == "qqH3J") vbf3j_tag = itag;
      
      if(histos[itag][ihist]->Integral() == 0) continue;
      //stack MCs and Data      
      if(tags[itag] == "Data"){
	if(histos_name[ihist] == "mass4l") histos[itag][ihist]->Rebin(4);
	int ip = 0;
	for(unsigned int ib=0; ib<(unsigned int)histos[itag][ihist]->GetNbinsX(); ++ib){
	  float x = histos[itag][ihist]->GetBinCenter(ib+1);
	  float y = histos[itag][ihist]->GetBinContent(ib+1); 
	  hdata[ihist]->SetPoint(ib,x,y);
	  if(y!=0 && x>maxx) maxx = x;
	  if(y > 0){
	    float emy = 0, epy = 0;
	    emy = -0.5 + sqrt(y+0.25);
	    epy = +0.5 + sqrt(y+0.25);
	    fdata[ihist]->SetPoint(ip,x,y);
	    fdata[ihist]->SetPointError(ip,0,0,emy,epy);
	    ++ip;
	    if(y+epy > maxy) maxy = y+epy;
	  }
	}
	hratio[ihist] = (TH1D*)histos[itag][ihist]->Clone("");
	hratio[ihist]->SetName("hratio_"+histos_name[ihist]);
	//std::cout<<"hratio in lokus Int = "<<hratio[ihist]->Integral()<<std::endl;
      }
      else{
	if(histos_name[ihist] == "mass4l") histos[itag][ihist]->Rebin(4);
	if(tags[itag] != "qqH3J"){
	  if(tags[itag] == "qqH"){
	    histos[itag][ihist]->SetLineColor(kPink-8);
	    histos[itag][ihist]->SetFillColor(kPink-9);
	  }
	  if(tags[itag] == "ggH"){
	    histos[itag][ihist]->SetLineColor(kPink-1);
	    histos[itag][ihist]->SetFillColor(kPink);
	  }
	  if(tags[itag] == "VH"){
	    histos[itag][ihist]->SetLineColor(kOrange-3);
	    histos[itag][ihist]->SetFillColor(kOrange-2);
	  }
	  if(tags[itag] == "ttH"){
	    histos[itag][ihist]->SetLineColor(kYellow+1);
	    histos[itag][ihist]->SetFillColor(kYellow);
	  }
	  if(tags[itag] == "qqZZ"){
	    histos[itag][ihist]->SetLineColor(kAzure-8);
	    histos[itag][ihist]->SetFillColor(kAzure-9);
	  }
	  if(tags[itag] == "ggZZ"){
	    histos[itag][ihist]->SetLineColor(kAzure-3);
	    histos[itag][ihist]->SetFillColor(kAzure-1);
	  }
	  shistos[ihist]->Add( histos[itag][ihist] );
	}
      }
    }
    
    
    if(ihist == 0){
      leg->AddEntry(fdata[ihist],"Data","pel");
      leg->AddEntry(histos[5][ihist],"qqH","f");
      leg->AddEntry(histos[3][ihist],"ggH","f");    
      leg->AddEntry(histos[4][ihist],"VH","f");    
      leg->AddEntry(histos[6][ihist],"ttH","f");    
      leg->AddEntry(histos[2][ihist],"qqZZ","f");    
      leg->AddEntry(histos[1][ihist],"ggZZ","f");    
    }
    if(histos_name[ihist] == "mass4l") leg->AddEntry(hzx,"Z+X","f");


    //std::cout<<"Getting stack from THStack... "<<shistos[ihist]->GetNhists()<<std::endl;
    if(shistos[ihist]->GetNhists() != 0){
      mc_sum[ihist] = ((TH1D*)(shistos[ihist]->GetStack()->Last()));
      if(ihist == 0) leg->AddEntry(mc_sum[ihist],"Stat. Uncertainty","f");
      float binning = mc_sum[ihist]->GetBinWidth(1);
      float xmin = hratio[ihist]->GetXaxis()->GetXmin();
      float xmax = hratio[ihist]->GetXaxis()->GetXmax();
      if(histos_name[ihist] == "mass4l"){
	xmin = 60;
	xmax = 300;
      }

      
      //cout<<"Drawing histograms..."<<endl;
      pad1[ihist]->cd();
      //Computes data-mc ratio
      //std::cout<<"hratio bins = "<<hratio[ihist]->GetNbinsX()<<std::endl;
      //std::cout<<"mc_sum bins = "<<mc_sum[ihist]->GetNbinsX()<<std::endl;
      hratio[ihist]->Divide(mc_sum[ihist]);
      std::cout<<"Fiting "<<histos_name[ihist]<<Form(", range: [%.2f, %.2f]",xmin,xmax)<<std::endl;
      hratio[ihist]->Fit(f1,"N","",xmin,xmax);
      double fit_value = f1->GetParameter(0);
      double fit_error = f1->GetParError(0);
      
      TGraphErrors *gfit = new TGraphErrors();
      gfit->SetPoint(0,xmin,fit_value);
      gfit->SetPointError(0,0.0,fit_error);
      gfit->SetPoint(1,xmax,fit_value);
      gfit->SetPointError(1,0.0,fit_error);
      gfit->SetMarkerSize(0);
      gfit->SetLineStyle(7);
      gfit->SetFillColor(kGray+2);
      gfit->SetFillStyle(3001);
  
      TLegend *rleg = new TLegend(0.8,0.75,0.95,0.95);
      rleg->SetFillColor(0);
      rleg->SetFillStyle(0);
      rleg->AddEntry(gfit,"Fit","fl");

      TLine *lfit = new TLine(xmin,1.0,xmax,1.0);
      lfit->SetLineColor(kRed);

      gfit->Draw("ale3");
      gfit->GetXaxis()->SetLabelSize(0.14);
      gfit->GetXaxis()->SetTitle(histos_name[ihist]);
      gfit->GetXaxis()->SetTitleSize(0.16);
      gfit->GetXaxis()->SetTitleOffset(1.0);
      gfit->GetYaxis()->SetLabelSize(0.13);
      gfit->GetYaxis()->SetTitle("Data/MC");
      gfit->GetYaxis()->SetTitleSize(0.14);
      gfit->GetYaxis()->SetTitleOffset(0.35);
      gfit->GetXaxis()->SetTickLength(0.06);
      gfit->GetXaxis()->SetTickLength(0.06);
      gfit->GetXaxis()->SetLimits(xmin,xmax);
      gfit->GetYaxis()->SetLimits(0,2);
      lfit->Draw();
      hratio[ihist]->SetMarkerColor(kBlack);
      hratio[ihist]->SetLineColor(kBlack);
      hratio[ihist]->Draw("pe,same");
      rleg->Draw();    
      
      
      float ymin = ((TH1D*)shistos[ihist]->GetStack()->First())->GetMaximum();
      float ymax = ((TH1D*)shistos[ihist]->GetStack()->Last())->GetMaximum();
      if(data_isPresent && fdata[ihist]->GetMaximum() > ymax) ymax = fdata[ihist]->GetMaximum();
  
      //std::cout<<"Plotting sum of MCs..."<<std::endl;
      mc_sum[ihist]->SetFillColor(kBlack);
      mc_sum[ihist]->SetLineColor(kBlack);
      mc_sum[ihist]->SetLineWidth(0);
      mc_sum[ihist]->SetFillStyle(3001);//3018);
      mc_sum[ihist]->SetMarkerStyle(0);
            
      pad2[ihist]->cd();
      //shistos[ihist]->GetXaxis()->SetTitle("m_{4l}");
      shistos[ihist]->GetYaxis()->SetTitle(Form("Events/%.3f",binning));
      shistos[ihist]->GetYaxis()->SetTitleOffset(1.0);
      shistos[ihist]->GetXaxis()->SetLimits(xmin,xmax);
      shistos[ihist]->Draw("hist");
      mc_sum[ihist]->Draw("e2,same");
      fdata[ihist]->SetMarkerColor(kBlack);
      fdata[ihist]->SetMarkerStyle(20);
      fdata[ihist]->Draw("ep,same");
      leg->Draw();
      cms_tag->Draw();
      count_tag->Draw();
      pad2[ihist]->RedrawAxis();
      fdata[ihist]->GetYaxis()->SetRangeUser(0,1.2*maxy);
      cv[ihist]->Update();
      cv[ihist]->Write();
      cv[ihist]->Print("Histograms/histo_"+histos_name[ihist]+".png");

      ++count_hists;
      std::cout<<"Histogram #= "<<count_hists<<std::endl;
    }
  }
  */
  
  
  TFile *histograms_plots = new TFile("Histograms.root","recreate");
  std::vector<TString> srs = {"smhiggs","vbf"};
  for(unsigned int isr=1; isr<(unsigned int)srs.size(); ++isr){
    TH1D *hzx_4mu = (TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/zxstatisticsOS/"+srs[isr]+"_4mu_nominal.root"))->Get("m4lZX");
    TH1D *hzx_4e = (TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/zxstatisticsOS/"+srs[isr]+"_4e_nominal.root"))->Get("m4lZX");
    TH1D *hzx_2e2mu = (TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/zxstatisticsOS/"+srs[isr]+"_2e2mu_nominal.root"))->Get("m4lZX");
    TH1D *hzx_2mu2e = (TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/zxstatisticsOS/"+srs[isr]+"_2mu2e_nominal.root"))->Get("m4lZX");
    
    histograms_plots->mkdir(srs[isr]);
    histograms_plots->mkdir(srs[isr]+"/4l");
    histograms_plots->cd(srs[isr]+"/4l");
    for(unsigned int itag=0; itag<ntags; ++itag){
      histos[itag][0][isr]->SetName(tags[itag]);
      histos[itag][0][isr]->Write();
      std::cout<<"Integral at 4l from "<<tags[itag]<<Form(" = %.3f",histos[itag][0][isr]->Integral())<<std::endl;
    }
    TH1D *hzx4l = (TH1D*)hzx_4mu->Clone();
    hzx4l->Add(hzx_4e);
    hzx4l->Add(hzx_2e2mu);
    hzx4l->Add(hzx_2mu2e);
    hzx4l->Write();
    
    if(srs[isr] == "vbf"){
      for(unsigned int itag=0; itag<ntags; ++itag){
	k57nj2[itag][0]->SetName("k57nj2_"+tags[itag]);
	k57nj2[itag][0]->Write();
	k24nj3[itag][0]->SetName("k24nj3_"+tags[itag]);
	k24nj3[itag][0]->Write();

	k57nj2_deltajj[itag][0]->SetName("k57nj2_deltajj_"+tags[itag]);
	k57nj2_deltajj[itag][0]->Write();
	k24nj3_deltajj[itag][0]->SetName("k24nj3_deltajj_"+tags[itag]);
	k24nj3_deltajj[itag][0]->Write();	
      }
      TH1D *k57nj2_zx = (TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/vbf/DistributionsUncertainties_k57nj2_proc_zx_4l_vbf_nominal.root"))->Get("ch4mu/zjets");
      k57nj2_zx->Add((TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/vbf/DistributionsUncertainties_k57nj2_proc_zx_4l_vbf_nominal.root"))->Get("ch4e/zjets"));
      k57nj2_zx->Add((TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/vbf/DistributionsUncertainties_k57nj2_proc_zx_4l_vbf_nominal.root"))->Get("ch2e2mu/zjets"));
      k57nj2_zx->SetName("k57nj2_zjets");
      histograms_plots->cd(srs[isr]+"/4l");
      k57nj2_zx->Write();
      
      TH1D *k24nj3_zx = (TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/vbf/DistributionsUncertainties_k24nj3_proc_zx_4l_vbf_nominal.root"))->Get("ch4mu/zjets");
      k24nj3_zx->Add((TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/vbf/DistributionsUncertainties_k24nj3_proc_zx_4l_vbf_nominal.root"))->Get("ch4e/zjets"));
      k24nj3_zx->Add((TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/vbf/DistributionsUncertainties_k24nj3_proc_zx_4l_vbf_nominal.root"))->Get("ch2e2mu/zjets"));
      k24nj3_zx->SetName("k24nj3_zjets");
      histograms_plots->cd(srs[isr]+"/4l");
      k24nj3_zx->Write();
      
      TH1D *k57nj2_deltajj_zx = (TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/vbf/DistributionsUncertainties_k57nj2_deltajj_proc_zx_4l_vbf_nominal.root"))->Get("ch4mu/zjets2D");
      k57nj2_deltajj_zx->Add((TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/vbf/DistributionsUncertainties_k57nj2_deltajj_proc_zx_4l_vbf_nominal.root"))->Get("ch4e/zjets2D"));
      k57nj2_deltajj_zx->Add((TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/vbf/DistributionsUncertainties_k57nj2_deltajj_proc_zx_4l_vbf_nominal.root"))->Get("ch2e2mu/zjets2D"));
      k57nj2_deltajj_zx->SetName("k57nj2_deltajj_zjets");
      histograms_plots->cd(srs[isr]+"/4l");
      k57nj2_deltajj_zx->Write();
      
      TH1D *k24nj3_deltajj_zx = (TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/vbf/DistributionsUncertainties_k24nj3_deltajj_proc_zx_4l_vbf_nominal.root"))->Get("ch4mu/zjets2D");
      k24nj3_deltajj_zx->Add((TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/vbf/DistributionsUncertainties_k24nj3_deltajj_proc_zx_4l_vbf_nominal.root"))->Get("ch4e/zjets2D"));
      k24nj3_deltajj_zx->Add((TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/vbf/DistributionsUncertainties_k24nj3_deltajj_proc_zx_4l_vbf_nominal.root"))->Get("ch2e2mu/zjets2D"));
      k24nj3_deltajj_zx->SetName("k24nj3_deltajj_zjets");
      histograms_plots->cd(srs[isr]+"/4l");
      k24nj3_deltajj_zx->Write();
    }

    histograms_plots->cd();
    histograms_plots->mkdir(srs[isr]+"/4mu");
    histograms_plots->cd(srs[isr]+"/4mu");
    for(unsigned int itag=0; itag<ntags; ++itag){
      histos[itag][1][isr]->SetName(tags[itag]);
      histos[itag][1][isr]->Write();
    }
    hzx_4mu->Write();

    if(srs[isr] == "vbf"){
      for(unsigned int itag=0; itag<ntags; ++itag){
	std::cout<<tags[itag]<<Form(", m4l = %.3f, k57nj2+k24nj3 = %.3f",histos[itag][1][isr]->Integral(),k57nj2[itag][1]->Integral()+k24nj3[itag][1]->Integral())<<std::endl;
	k57nj2[itag][1]->SetName("k57nj2_"+tags[itag]);
	k57nj2[itag][1]->Write();
	k24nj3[itag][1]->SetName("k24nj3_"+tags[itag]);
	k24nj3[itag][1]->Write();

	k57nj2_deltajj[itag][1]->SetName("k57nj2_deltajj_"+tags[itag]);
	k57nj2_deltajj[itag][1]->Write();
	k24nj3_deltajj[itag][1]->SetName("k24nj3_deltajj_"+tags[itag]);
	k24nj3_deltajj[itag][1]->Write();	
      }
      TH1D *k57nj2_zx = (TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/vbf/DistributionsUncertainties_k57nj2_proc_zx_4l_vbf_nominal.root"))->Get("ch4mu/zjets");
      k57nj2_zx->SetName("k57nj2_zjets");
      histograms_plots->cd(srs[isr]+"/4mu");
      k57nj2_zx->Write();
      
      TH1D *k24nj3_zx = (TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/vbf/DistributionsUncertainties_k24nj3_proc_zx_4l_vbf_nominal.root"))->Get("ch4mu/zjets");
      k24nj3_zx->SetName("k24nj3_zjets");
      histograms_plots->cd(srs[isr]+"/4mu");
      k24nj3_zx->Write();
    }
    
    histograms_plots->cd();
    histograms_plots->mkdir(srs[isr]+"/4e");
    histograms_plots->cd(srs[isr]+"/4e");
    for(unsigned int itag=0; itag<ntags; ++itag){
      histos[itag][2][isr]->SetName(tags[itag]);
      histos[itag][2][isr]->Write();
    }
    hzx_4mu->Write();

    if(srs[isr] == "vbf"){
      for(unsigned int itag=0; itag<ntags; ++itag){
	std::cout<<tags[itag]<<Form(", m4l = %.3f, k57nj2+k24nj3 = %.3f",histos[itag][2][isr]->Integral(),k57nj2[itag][2]->Integral()+k24nj3[itag][2]->Integral())<<std::endl;
	k57nj2[itag][2]->SetName("k57nj2_"+tags[itag]);
	k57nj2[itag][2]->Write();
	k24nj3[itag][2]->SetName("k24nj3_"+tags[itag]);
	k24nj3[itag][2]->Write();
	
	k57nj2_deltajj[itag][2]->SetName("k57nj2_deltajj_"+tags[itag]);
	k57nj2_deltajj[itag][2]->Write();
	k24nj3_deltajj[itag][2]->SetName("k24nj3_deltajj_"+tags[itag]);
	k24nj3_deltajj[itag][2]->Write();	
      }
      TH1D *k57nj2_zx = (TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/vbf/DistributionsUncertainties_k57nj2_proc_zx_4l_vbf_nominal.root"))->Get("ch4e/zjets");
      k57nj2_zx->SetName("k57nj2_zjets");
      histograms_plots->cd(srs[isr]+"/4e");
      k57nj2_zx->Write();
      
      TH1D *k24nj3_zx = (TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/vbf/DistributionsUncertainties_k24nj3_proc_zx_4l_vbf_nominal.root"))->Get("ch4e/zjets");
      k24nj3_zx->SetName("k24nj3_zjets");
      histograms_plots->cd(srs[isr]+"/4e");
      k24nj3_zx->Write();
    }
    
    histograms_plots->cd();
    histograms_plots->mkdir(srs[isr]+"/2e2mu");
    histograms_plots->cd(srs[isr]+"/2e2mu");
    for(unsigned int itag=0; itag<ntags; ++itag){
      histos[itag][3][isr]->SetName(tags[itag]);
      histos[itag][3][isr]->Write();
    }
    hzx_4mu->Write();

    if(srs[isr] == "vbf"){
      for(unsigned int itag=0; itag<ntags; ++itag){
	std::cout<<tags[itag]<<Form(", m4l = %.3f, k57nj2+k24nj3 = %.3f",histos[itag][3][isr]->Integral(),k57nj2[itag][3]->Integral()+k24nj3[itag][3]->Integral())<<std::endl;
	k57nj2[itag][3]->SetName("k57nj2_"+tags[itag]);
	k57nj2[itag][3]->Write();
	k24nj3[itag][3]->SetName("k24nj3_"+tags[itag]);
	k24nj3[itag][3]->Write();

	k57nj2_deltajj[itag][3]->SetName("k57nj2_deltajj_"+tags[itag]);
	k57nj2_deltajj[itag][3]->Write();
	k24nj3_deltajj[itag][3]->SetName("k24nj3_deltajj_"+tags[itag]);
	k24nj3_deltajj[itag][3]->Write();
      }
      TH1D *k57nj2_zx = (TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/vbf/DistributionsUncertainties_k57nj2_proc_zx_4l_vbf_nominal.root"))->Get("ch2e2mu/zjets");
      k57nj2_zx->SetName("k57nj2_zjets");
      histograms_plots->cd(srs[isr]+"/2e2mu");
      k57nj2_zx->Write();
      
      TH1D *k24nj3_zx = (TH1D*)(TFile::Open("/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/vbf/DistributionsUncertainties_k24nj3_proc_zx_4l_vbf_nominal.root"))->Get("ch2e2mu/zjets");
      k24nj3_zx->SetName("k24nj3_zjets");
      histograms_plots->cd(srs[isr]+"/2e2mu");
      k24nj3_zx->Write();
    }    
  }
  histograms_plots->Close();

  std::cout<<"------------------------------------------------------"<<std::endl;
  std::cout<<Form("S(4mu):   %.2f +/- %.2f, B(4mu):   %.2f +/- %.2f",sig_4mu,sqrt(sig_4mu_err),bkg_4mu,sqrt(bkg_4mu_err))<<std::endl;
  std::cout<<Form("S(4e):    %.2f +/- %.2f, B(4e):    %.2f +/- %.2f",sig_4e,sqrt(sig_4e_err),bkg_4e,sqrt(bkg_4e_err))<<std::endl;
  std::cout<<Form("S(2e2mu): %.2f +/- %.2f, B(2e2mu): %.2f +/- %.2f",sig_2e2mu,sqrt(sig_2e2mu_err),bkg_2e2mu,sqrt(bkg_2e2mu_err))<<std::endl;
  std::cout<<Form("S(4l):    %.2f +/- %.2f, B(4l):    %.2f +/- %.2f",fnsig,sqrt(fsig_err),fnbkg,sqrt(fbkg_err))<<std::endl;
  
  
  return;
}//End plot function


void plotObjectsProperties(void){
  std::cout<<"Plots of variables will be stored in the file: Histograms.root"<<std::endl;
  gSystem->Exec("mkdir -p Histograms");

  
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
  std::vector<TString> VV;
  std::vector<TString> VVV;
  std::vector<TString> Wjets;
    
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

  samples.push_back(data_obs);	tags.push_back("Data");
  samples.push_back(ggZZ);	tags.push_back("ggZZ");
  samples.push_back(qqZZ);	tags.push_back("qqZZ");
  samples.push_back(ggH);	tags.push_back("ggH");
  samples.push_back(WH);	tags.push_back("WH");
  samples.push_back(ZH);	tags.push_back("ZH");
  samples.push_back(ttH);	tags.push_back("ttH");
  samples.push_back(qqH);	tags.push_back("qqH");
  samples.push_back(qqH3J);	tags.push_back("qqH3J");
  
  //samples.push_back(HWW);	tags.push_back("HWW");
  //samples.push_back(TTV);	tags.push_back("TTV");
  //samples.push_back(VVV);	tags.push_back("VVV");
  //samples.push_back(VV);	tags.push_back("VV");
  //samples.push_back(Wjets);	tags.push_back("Wjets");
  //samples.push_back(TT4l);	tags.push_back("TT4l");
  
  loopSamples(tags, samples);
     
  //Ends program
  return;
}
