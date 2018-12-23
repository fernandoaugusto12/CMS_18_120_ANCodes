//NOTE: NN computed at the stage of Z+X selections
//They will be stored in a matrix [NN][nominal+shifts]
#include "/lustre/cms/store/user/mmelodea/Keras/modelsNN/k57nj2.h"
#include "/lustre/cms/store/user/mmelodea/Keras/modelsNN/k24nj3.h"

const int nNNs = 8;
const int nShifts = 26;//nominal+shifts


void IncludeKerasDiscs(void){
  Float_t f_lept1_pt, f_lept1_pt_error, f_lept1_eta, f_lept1_phi;
  Float_t f_lept2_pt, f_lept2_pt_error, f_lept2_eta, f_lept2_phi;
  Float_t f_lept3_pt, f_lept3_pt_error, f_lept3_eta, f_lept3_phi;
  Float_t f_lept4_pt, f_lept4_pt_error, f_lept4_eta, f_lept4_phi;
  Float_t f_pfmet, f_pfmet_ElectronEnDn, f_pfmet_ElectronEnUp, f_pfmet_JetEnDn, f_pfmet_JetEnUp, f_pfmet_JetResDn, f_pfmet_JetResUp;
  Float_t f_pfmet_MuonEnDn, f_pfmet_MuonEnUp, f_pfmet_PhotonEnDn, f_pfmet_PhotonEnUp, f_pfmet_UnclusteredEnDn, f_pfmet_UnclusteredEnUp;
  Float_t f_jets_highpt_pt[4], f_jets_highpt_pt_error[4], f_jets_highpt_eta[4], f_jets_highpt_phi[4];
  Float_t f_k57nj2[nShifts], f_k24nj3[nShifts];
  
  
  TFile *ofile = TFile::Open("../Merged/output_ZplusX_DataRunIISummer16_03Feb2017.root");
  TTree *otreeFR = (TTree*)ofile->Get("FakeRate");
  TTree *otreeCR = (TTree*)ofile->Get("ControlRegions");

  TFile *ffile = new TFile("../Merged/output_ZplusX_DataRunIISummer16_03Feb2017_withKerasDiscriminants.root","recreate");
  TTree *ftreeFR = otreeFR->CloneTree();
  TTree *ftreeCR = otreeCR->CloneTree();
  ftreeCR->SetBranchAddress("f_lept1_pt", &f_lept1_pt);
  ftreeCR->SetBranchAddress("f_lept1_pt_error", &f_lept1_pt_error);
  ftreeCR->SetBranchAddress("f_lept1_eta", &f_lept1_eta);
  ftreeCR->SetBranchAddress("f_lept1_phi", &f_lept1_phi);
  ftreeCR->SetBranchAddress("f_lept2_pt", &f_lept2_pt);
  ftreeCR->SetBranchAddress("f_lept2_pt_error", &f_lept2_pt_error);
  ftreeCR->SetBranchAddress("f_lept2_eta", &f_lept2_eta);
  ftreeCR->SetBranchAddress("f_lept2_phi", &f_lept2_phi);
  ftreeCR->SetBranchAddress("f_lept3_pt", &f_lept3_pt);
  ftreeCR->SetBranchAddress("f_lept3_pt_error", &f_lept3_pt_error);
  ftreeCR->SetBranchAddress("f_lept3_eta", &f_lept3_eta);
  ftreeCR->SetBranchAddress("f_lept3_phi", &f_lept3_phi);
  ftreeCR->SetBranchAddress("f_lept4_pt", &f_lept4_pt);
  ftreeCR->SetBranchAddress("f_lept4_pt_error", &f_lept4_pt_error);
  ftreeCR->SetBranchAddress("f_lept4_eta", &f_lept4_eta);
  ftreeCR->SetBranchAddress("f_lept4_phi", &f_lept4_phi);
  ftreeCR->SetBranchAddress("f_pfmet", &f_pfmet);
  ftreeCR->SetBranchAddress("f_pfmet_JetEnUp", &f_pfmet_JetEnUp);
  ftreeCR->SetBranchAddress("f_pfmet_JetEnDn", &f_pfmet_JetEnDn);
  ftreeCR->SetBranchAddress("f_pfmet_ElectronEnUp", &f_pfmet_ElectronEnUp);
  ftreeCR->SetBranchAddress("f_pfmet_ElectronEnDn", &f_pfmet_ElectronEnDn);
  ftreeCR->SetBranchAddress("f_pfmet_MuonEnUp", &f_pfmet_MuonEnUp);
  ftreeCR->SetBranchAddress("f_pfmet_MuonEnDn", &f_pfmet_MuonEnDn);
  ftreeCR->SetBranchAddress("f_pfmet_JetResUp", &f_pfmet_JetResUp);
  ftreeCR->SetBranchAddress("f_pfmet_JetResDn", &f_pfmet_JetResDn);
  ftreeCR->SetBranchAddress("f_pfmet_UnclusteredEnUp", &f_pfmet_UnclusteredEnUp);
  ftreeCR->SetBranchAddress("f_pfmet_UnclusteredEnDn", &f_pfmet_UnclusteredEnDn);
  ftreeCR->SetBranchAddress("f_pfmet_PhotonEnUp", &f_pfmet_PhotonEnUp);
  ftreeCR->SetBranchAddress("f_pfmet_PhotonEnDn", &f_pfmet_PhotonEnDn);
  ftreeCR->SetBranchAddress("f_jets_highpt_pt", &f_jets_highpt_pt);
  ftreeCR->SetBranchAddress("f_jets_highpt_pt_error", &f_jets_highpt_pt_error);
  ftreeCR->SetBranchAddress("f_jets_highpt_eta", &f_jets_highpt_eta);
  ftreeCR->SetBranchAddress("f_jets_highpt_phi", &f_jets_highpt_phi);
  
  
  TBranch *bk57nj2 = ftreeCR->Branch("f_k57nj2", &f_k57nj2, Form("f_k57nj2[%i]/F",nShifts));
  TBranch *bk24nj3 = ftreeCR->Branch("f_k24nj3", &f_k24nj3, Form("f_k24nj3[%i]/F",nShifts));
  
  unsigned int nCR = ftreeCR->GetEntries();
  for(unsigned int iev=0; iev<nCR; ++iev){
    ftreeCR->GetEntry(iev);
    
     //========= Gets the NNs ============
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

     float j2pt  = f_jets_highpt_pt[1];
     float j2pt_err  = f_jets_highpt_pt_error[1];
     float j2eta = f_jets_highpt_eta[1];
     float j2phi = f_jets_highpt_phi[1];

     float j3pt  = (f_jets_highpt_pt[2] != -999)? f_jets_highpt_pt[2]: 0;
     float j3pt_err  = (f_jets_highpt_pt[2] != -999)? f_jets_highpt_pt_error[2]: 0;
     float j3eta = (f_jets_highpt_pt[2] != -999)? f_jets_highpt_eta[2]: 0;
     float j3phi = (f_jets_highpt_pt[2] != -999)? f_jets_highpt_phi[2]: 0;
	
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
     
	   std::vector<std::vector<double> > inputs = {
	     //k57nj2
	     {l1pt,l1eta,l1phi,l2pt,l2eta,l2phi,l3pt,l3eta,l3phi,l4pt,l4eta,l4phi,j1pt,j1eta,j1phi,j2pt,j2eta,j2phi,met},
	     //k24nj3
	     {l1pt,l1eta,l1phi,l2pt,l2eta,l2phi,l3pt,l3eta,l3phi,l4pt,l4eta,l4phi,j1pt,j1eta,j1phi,j2pt,j2eta,j2phi,j3pt,j3eta,j3phi,met},
	   };
	   
	   f_k57nj2[ishift] = model_k57nj2(inputs[0]);
	   f_k24nj3[ishift] = model_k24nj3(inputs[1]);
	   
	   ++ishift;
       }
     }
     
     bk57nj2->Fill();
     bk24nj3->Fill();
  }
  
  ffile->cd();
  ftreeFR->Write();
  ftreeCR->Write();
  ffile->Close();
  
  return;
} 
