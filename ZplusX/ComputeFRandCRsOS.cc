#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TSpline.h>
#include <TFile.h>
#include <TTree.h>

#include <TMath.h>
#include <math.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>

#include <stdlib.h>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdlib>
#include <stdio.h>
//#include <libgen.h>

#include <utility>

#include "ZZMatrixElement/MELA/src/computeAngles.h" //me comment
#include "ZZMatrixElement/MELA/src/computeAngles.cc" //me comment
#include "ZZMatrixElement/MEMCalculators/interface/MEMCalculators.h" //me comment
//#include "PhysicsTools/Utilities/interface/LumiReWeighting.h" //me comment


//NOTE: NN computed at the stage of Z+X selections
//They will be stored in a matrix [NN][nominal+shifts]
/*
#include "/lustrehome/mmelodea/Keras/modelsNN/k41nj2.h"
#include "/lustrehome/mmelodea/Keras/modelsNN/k1177nj2.h"
#include "/lustrehome/mmelodea/Keras/modelsNN/k2nj3.h"
#include "/lustrehome/mmelodea/Keras/modelsNN/k4nj3.h"
#include "/lustrehome/mmelodea/Keras/modelsNN/k16nj2e3.h"
#include "/lustrehome/mmelodea/Keras/modelsNN/k53nj2e3.h"
#include "/lustrehome/mmelodea/Keras/modelsNN/k196nj2e3.h"
#include "/lustrehome/mmelodea/Keras/modelsNN/k543nj2e3.h"
const int nNNs = 8;
const int nShifts = 29;//nominal+shifts
*/

using namespace std;
using std::map;
using std::pair;
// using namespace RooFit;
// using namespace meMCFM;
//using namespace MEMNames; //me comment 

TString auxiliar_files = "/lustrehome/mmelodea/CMSSW_8_0_13/src/HiggsAnalysis/HiggsToZZ4Leptons/test/macros/";
const double Zmass = 91.188; // nominal Z boson mass
const double mPI = 3.141592654;

//Define functions
//=================
double DELTAPHI( double phi1, double phi2 ){
  
  if( phi1 > mPI || phi1 < -mPI || phi2 > mPI || phi2 < -mPI)
    {
      //        // cout << "Angles out of range!!! " << endl;
      //        // cout << " phi1 " << phi1 << endl;
      //        // cout << " phi2 " << phi2 << endl;
      return -999;
    }
  float dp=std::abs(phi1-phi2);
  if (dp>mPI) dp-=float(2*mPI);
  //return dp;
  return  min( fabs( phi1 - phi2 ) , 2*mPI - fabs( phi1 - phi2 ) ) ;
  
}

//   double dphi = phi1 - phi2;
//   while (dphi > TMath::Pi()) dphi -= 2 * TMath::Pi();
//   while (dphi <= -TMath::Pi()) dphi += 2 * TMath::Pi();
//   return dphi;
// }

double EAele(double recoele_scl_eta,bool use2011EA){
  /////Reham   
  double EffectiveArea=0.;
  if (use2011EA){
    if (fabs(recoele_scl_eta) >= 0.0   && fabs(recoele_scl_eta) < 1.0 )   EffectiveArea = 0.18;
    if (fabs(recoele_scl_eta) >= 1.0   && fabs(recoele_scl_eta) < 1.479 ) EffectiveArea = 0.20;
    if (fabs(recoele_scl_eta) >= 1.479 && fabs(recoele_scl_eta) < 2.0 )   EffectiveArea = 0.15;
    if (fabs(recoele_scl_eta) >= 2.0   && fabs(recoele_scl_eta) < 2.2 )   EffectiveArea = 0.19;
    if (fabs(recoele_scl_eta) >= 2.2   && fabs(recoele_scl_eta) < 2.3 )   EffectiveArea = 0.21;
    if (fabs(recoele_scl_eta) >= 2.3   && fabs(recoele_scl_eta) < 2.4 )   EffectiveArea = 0.22;
    if (fabs(recoele_scl_eta) >= 2.4 )                                       EffectiveArea = 0.29;
  }
  else {
    // if (fabs(RECOELE_ETA[index]) >= 0.0   && fabs(RECOELE_ETA[index]) < 0.8 )   EffectiveArea = 0.1830;
    // if (fabs(RECOELE_ETA[index]) >= 0.8   && fabs(RECOELE_ETA[index]) < 1.3 )   EffectiveArea = 0.1734;
    // if (fabs(RECOELE_ETA[index]) >= 1.3   && fabs(RECOELE_ETA[index]) < 2.0 )   EffectiveArea = 0.1077;
    // if (fabs(RECOELE_ETA[index]) >= 2.0   && fabs(RECOELE_ETA[index]) < 2.2 )   EffectiveArea = 0.1565;
    // if (fabs(RECOELE_ETA[index]) >= 2.2 )                                       EffectiveArea = 0.2680;
    
    if (fabs(recoele_scl_eta) >= 0.0   && fabs(recoele_scl_eta) < 1.0 )   EffectiveArea = 0.1752;
    if (fabs(recoele_scl_eta) >= 1.0   && fabs(recoele_scl_eta) < 1.479 ) EffectiveArea = 0.1862;
    if (fabs(recoele_scl_eta) >= 1.479 && fabs(recoele_scl_eta) < 2.0 )   EffectiveArea = 0.1411;
    if (fabs(recoele_scl_eta) >= 2.0   && fabs(recoele_scl_eta) < 2.2 )   EffectiveArea = 0.1534;
    if (fabs(recoele_scl_eta) >= 2.2   && fabs(recoele_scl_eta) < 2.3 )   EffectiveArea = 0.1903;
    if (fabs(recoele_scl_eta) >= 2.3   && fabs(recoele_scl_eta) < 2.4 )   EffectiveArea = 0.2243;
    if (fabs(recoele_scl_eta) >= 2.4   && fabs(recoele_scl_eta) < 5.0 )   EffectiveArea = 0.2687;
    
  }
  
  return EffectiveArea;
  
}




float kfactor_qqZZ_qcd_dPhi(float GENdPhiZZ, int finalState)
{
  
  // finalState=1 : 4e/4mu/4tau
  // finalState=2 : 2e2mu/2mutau/2e2tau
  
  float k=0.0;
  
  if (finalState==1) {        
    k+=1.515838921760*(abs(GENdPhiZZ)>0.0&&abs(GENdPhiZZ)<=0.1);
    k+=1.496256665410*(abs(GENdPhiZZ)>0.1&&abs(GENdPhiZZ)<=0.2);
    k+=1.495522061910*(abs(GENdPhiZZ)>0.2&&abs(GENdPhiZZ)<=0.3);
    k+=1.483273154250*(abs(GENdPhiZZ)>0.3&&abs(GENdPhiZZ)<=0.4);
    k+=1.465589701130*(abs(GENdPhiZZ)>0.4&&abs(GENdPhiZZ)<=0.5);
    k+=1.491500887510*(abs(GENdPhiZZ)>0.5&&abs(GENdPhiZZ)<=0.6);
    k+=1.441183580450*(abs(GENdPhiZZ)>0.6&&abs(GENdPhiZZ)<=0.7);
    k+=1.440830603990*(abs(GENdPhiZZ)>0.7&&abs(GENdPhiZZ)<=0.8);
    k+=1.414339019120*(abs(GENdPhiZZ)>0.8&&abs(GENdPhiZZ)<=0.9);
    k+=1.422534218560*(abs(GENdPhiZZ)>0.9&&abs(GENdPhiZZ)<=1.0);
    k+=1.401037066000*(abs(GENdPhiZZ)>1.0&&abs(GENdPhiZZ)<=1.1);
    k+=1.408539428810*(abs(GENdPhiZZ)>1.1&&abs(GENdPhiZZ)<=1.2);
    k+=1.381247744080*(abs(GENdPhiZZ)>1.2&&abs(GENdPhiZZ)<=1.3);
    k+=1.370553357430*(abs(GENdPhiZZ)>1.3&&abs(GENdPhiZZ)<=1.4);
    k+=1.347323316000*(abs(GENdPhiZZ)>1.4&&abs(GENdPhiZZ)<=1.5);
    k+=1.340113437450*(abs(GENdPhiZZ)>1.5&&abs(GENdPhiZZ)<=1.6);
    k+=1.312661036510*(abs(GENdPhiZZ)>1.6&&abs(GENdPhiZZ)<=1.7);
    k+=1.290055062010*(abs(GENdPhiZZ)>1.7&&abs(GENdPhiZZ)<=1.8);
    k+=1.255322614790*(abs(GENdPhiZZ)>1.8&&abs(GENdPhiZZ)<=1.9);
    k+=1.254455642450*(abs(GENdPhiZZ)>1.9&&abs(GENdPhiZZ)<=2.0);
    k+=1.224047664420*(abs(GENdPhiZZ)>2.0&&abs(GENdPhiZZ)<=2.1);
    k+=1.178816782670*(abs(GENdPhiZZ)>2.1&&abs(GENdPhiZZ)<=2.2);
    k+=1.162624827140*(abs(GENdPhiZZ)>2.2&&abs(GENdPhiZZ)<=2.3);
    k+=1.105401140940*(abs(GENdPhiZZ)>2.3&&abs(GENdPhiZZ)<=2.4);
    k+=1.074749265690*(abs(GENdPhiZZ)>2.4&&abs(GENdPhiZZ)<=2.5);
    k+=1.021864599380*(abs(GENdPhiZZ)>2.5&&abs(GENdPhiZZ)<=2.6);
    k+=0.946334793286*(abs(GENdPhiZZ)>2.6&&abs(GENdPhiZZ)<=2.7);
    k+=0.857458082628*(abs(GENdPhiZZ)>2.7&&abs(GENdPhiZZ)<=2.8);
    k+=0.716607670482*(abs(GENdPhiZZ)>2.8&&abs(GENdPhiZZ)<=2.9);
    k+=1.132841784840*(abs(GENdPhiZZ)>2.9&&abs(GENdPhiZZ)<=3.1416);
  }
  
  if (finalState==2) {
    k+=1.513834489150*(abs(GENdPhiZZ)>0.0&&abs(GENdPhiZZ)<=0.1);
    k+=1.541738780180*(abs(GENdPhiZZ)>0.1&&abs(GENdPhiZZ)<=0.2);
    k+=1.497829632510*(abs(GENdPhiZZ)>0.2&&abs(GENdPhiZZ)<=0.3);
    k+=1.534956782920*(abs(GENdPhiZZ)>0.3&&abs(GENdPhiZZ)<=0.4);
    k+=1.478217033060*(abs(GENdPhiZZ)>0.4&&abs(GENdPhiZZ)<=0.5);
    k+=1.504330859290*(abs(GENdPhiZZ)>0.5&&abs(GENdPhiZZ)<=0.6);
    k+=1.520626246850*(abs(GENdPhiZZ)>0.6&&abs(GENdPhiZZ)<=0.7);
    k+=1.507013090030*(abs(GENdPhiZZ)>0.7&&abs(GENdPhiZZ)<=0.8);
    k+=1.494243156250*(abs(GENdPhiZZ)>0.8&&abs(GENdPhiZZ)<=0.9);
    k+=1.450536096150*(abs(GENdPhiZZ)>0.9&&abs(GENdPhiZZ)<=1.0);
    k+=1.460812521660*(abs(GENdPhiZZ)>1.0&&abs(GENdPhiZZ)<=1.1);
    k+=1.471603622200*(abs(GENdPhiZZ)>1.1&&abs(GENdPhiZZ)<=1.2);
    k+=1.467700038200*(abs(GENdPhiZZ)>1.2&&abs(GENdPhiZZ)<=1.3);
    k+=1.422408690640*(abs(GENdPhiZZ)>1.3&&abs(GENdPhiZZ)<=1.4);
    k+=1.397184022730*(abs(GENdPhiZZ)>1.4&&abs(GENdPhiZZ)<=1.5);
    k+=1.375593447520*(abs(GENdPhiZZ)>1.5&&abs(GENdPhiZZ)<=1.6);
    k+=1.391901318370*(abs(GENdPhiZZ)>1.6&&abs(GENdPhiZZ)<=1.7);
    k+=1.368564350560*(abs(GENdPhiZZ)>1.7&&abs(GENdPhiZZ)<=1.8);
    k+=1.317884804290*(abs(GENdPhiZZ)>1.8&&abs(GENdPhiZZ)<=1.9);
    k+=1.314019950800*(abs(GENdPhiZZ)>1.9&&abs(GENdPhiZZ)<=2.0);
    k+=1.274641749910*(abs(GENdPhiZZ)>2.0&&abs(GENdPhiZZ)<=2.1);
    k+=1.242346606820*(abs(GENdPhiZZ)>2.1&&abs(GENdPhiZZ)<=2.2);
    k+=1.244727403840*(abs(GENdPhiZZ)>2.2&&abs(GENdPhiZZ)<=2.3);
    k+=1.146259351670*(abs(GENdPhiZZ)>2.3&&abs(GENdPhiZZ)<=2.4);
    k+=1.107804993520*(abs(GENdPhiZZ)>2.4&&abs(GENdPhiZZ)<=2.5);
    k+=1.042053646740*(abs(GENdPhiZZ)>2.5&&abs(GENdPhiZZ)<=2.6);
    k+=0.973608545141*(abs(GENdPhiZZ)>2.6&&abs(GENdPhiZZ)<=2.7);
    k+=0.872169942668*(abs(GENdPhiZZ)>2.7&&abs(GENdPhiZZ)<=2.8);
    k+=0.734505279177*(abs(GENdPhiZZ)>2.8&&abs(GENdPhiZZ)<=2.9);
    k+=1.163152837230*(abs(GENdPhiZZ)>2.9&&abs(GENdPhiZZ)<=3.1416);       
  }
  if (k==0.0) return 1.1; // if something goes wrong return inclusive k-factor
  else return k;
  
}

float kfactor_qqZZ_qcd_M(float GENmassZZ, int finalState)
{
  
  // finalState=1 : 4e/4mu/4tau
  // finalState=2 : 2e2mu/2mutau/2e2tau
  
  float k=0.0;
  
  if (finalState==1) {
    k+=1.23613311013*(abs(GENmassZZ)>0.0&&abs(GENmassZZ)<=25.0);
    k+=1.17550314639*(abs(GENmassZZ)>25.0&&abs(GENmassZZ)<=50.0);
    k+=1.17044565911*(abs(GENmassZZ)>50.0&&abs(GENmassZZ)<=75.0);
    k+=1.03141209689*(abs(GENmassZZ)>75.0&&abs(GENmassZZ)<=100.0);
    k+=1.05285574912*(abs(GENmassZZ)>100.0&&abs(GENmassZZ)<=125.0);
    k+=1.11287217794*(abs(GENmassZZ)>125.0&&abs(GENmassZZ)<=150.0);
    k+=1.13361441158*(abs(GENmassZZ)>150.0&&abs(GENmassZZ)<=175.0);
    k+=1.10355603327*(abs(GENmassZZ)>175.0&&abs(GENmassZZ)<=200.0);
    k+=1.10053981637*(abs(GENmassZZ)>200.0&&abs(GENmassZZ)<=225.0);
    k+=1.10972676811*(abs(GENmassZZ)>225.0&&abs(GENmassZZ)<=250.0);
    k+=1.12069120525*(abs(GENmassZZ)>250.0&&abs(GENmassZZ)<=275.0);
    k+=1.11589101635*(abs(GENmassZZ)>275.0&&abs(GENmassZZ)<=300.0);
    k+=1.13906170314*(abs(GENmassZZ)>300.0&&abs(GENmassZZ)<=325.0);
    k+=1.14854594271*(abs(GENmassZZ)>325.0&&abs(GENmassZZ)<=350.0);
    k+=1.14616229031*(abs(GENmassZZ)>350.0&&abs(GENmassZZ)<=375.0);
    k+=1.14573157789*(abs(GENmassZZ)>375.0&&abs(GENmassZZ)<=400.0);
    k+=1.13829430515*(abs(GENmassZZ)>400.0&&abs(GENmassZZ)<=425.0);
    k+=1.15521193686*(abs(GENmassZZ)>425.0&&abs(GENmassZZ)<=450.0);
    k+=1.13679822698*(abs(GENmassZZ)>450.0&&abs(GENmassZZ)<=475.0);
    k+=1.13223956942*(abs(GENmassZZ)>475.0);
  }
  
  if (finalState==2) {
    k+=1.25094466582*(abs(GENmassZZ)>0.0&&abs(GENmassZZ)<=25.0);
    k+=1.22459455362*(abs(GENmassZZ)>25.0&&abs(GENmassZZ)<=50.0);
    k+=1.19287368979*(abs(GENmassZZ)>50.0&&abs(GENmassZZ)<=75.0);
    k+=1.04597506451*(abs(GENmassZZ)>75.0&&abs(GENmassZZ)<=100.0);
    k+=1.08323413771*(abs(GENmassZZ)>100.0&&abs(GENmassZZ)<=125.0);
    k+=1.09994968030*(abs(GENmassZZ)>125.0&&abs(GENmassZZ)<=150.0);
    k+=1.16698455800*(abs(GENmassZZ)>150.0&&abs(GENmassZZ)<=175.0);
    k+=1.10399053155*(abs(GENmassZZ)>175.0&&abs(GENmassZZ)<=200.0);
    k+=1.10592664340*(abs(GENmassZZ)>200.0&&abs(GENmassZZ)<=225.0);
    k+=1.10690381480*(abs(GENmassZZ)>225.0&&abs(GENmassZZ)<=250.0);
    k+=1.11194928918*(abs(GENmassZZ)>250.0&&abs(GENmassZZ)<=275.0);
    k+=1.13522586553*(abs(GENmassZZ)>275.0&&abs(GENmassZZ)<=300.0);
    k+=1.11895090244*(abs(GENmassZZ)>300.0&&abs(GENmassZZ)<=325.0);
    k+=1.13898508615*(abs(GENmassZZ)>325.0&&abs(GENmassZZ)<=350.0);
    k+=1.15463977506*(abs(GENmassZZ)>350.0&&abs(GENmassZZ)<=375.0);
    k+=1.17341664594*(abs(GENmassZZ)>375.0&&abs(GENmassZZ)<=400.0);
    k+=1.20093349763*(abs(GENmassZZ)>400.0&&abs(GENmassZZ)<=425.0);
    k+=1.18915554919*(abs(GENmassZZ)>425.0&&abs(GENmassZZ)<=450.0);
    k+=1.18546007375*(abs(GENmassZZ)>450.0&&abs(GENmassZZ)<=475.0);
    k+=1.12864505708*(abs(GENmassZZ)>475.0);
    }
  
  if (k==0.0) return 1.1;
  else return k; // if something goes wrong return inclusive k-factor
  
}

float kfactor_qqZZ_qcd_Pt(float GENpTZZ, int finalState)
{
  
  // finalState=1 : 4e/4mu/4tau
  // finalState=2 : 2e2mu/2mutau/2e2tau
  
  float k=0.0;
  
  if (finalState==1) {
    k+=0.64155491983*(abs(GENpTZZ)>0.0&&abs(GENpTZZ)<=5.0);
    k+=1.09985240531*(abs(GENpTZZ)>5.0&&abs(GENpTZZ)<=10.0);
    k+=1.29390628654*(abs(GENpTZZ)>10.0&&abs(GENpTZZ)<=15.0);
    k+=1.37859998571*(abs(GENpTZZ)>15.0&&abs(GENpTZZ)<=20.0);
    k+=1.42430263312*(abs(GENpTZZ)>20.0&&abs(GENpTZZ)<=25.0);
    k+=1.45038493266*(abs(GENpTZZ)>25.0&&abs(GENpTZZ)<=30.0);
    k+=1.47015377651*(abs(GENpTZZ)>30.0&&abs(GENpTZZ)<=35.0);
    k+=1.48828685748*(abs(GENpTZZ)>35.0&&abs(GENpTZZ)<=40.0);
    k+=1.50573440448*(abs(GENpTZZ)>40.0&&abs(GENpTZZ)<=45.0);
    k+=1.50211655928*(abs(GENpTZZ)>45.0&&abs(GENpTZZ)<=50.0);
    k+=1.50918720827*(abs(GENpTZZ)>50.0&&abs(GENpTZZ)<=55.0);
    k+=1.52463089491*(abs(GENpTZZ)>55.0&&abs(GENpTZZ)<=60.0);
    k+=1.52400838378*(abs(GENpTZZ)>60.0&&abs(GENpTZZ)<=65.0);
    k+=1.52418067701*(abs(GENpTZZ)>65.0&&abs(GENpTZZ)<=70.0);
    k+=1.55424382578*(abs(GENpTZZ)>70.0&&abs(GENpTZZ)<=75.0);
    k+=1.52544284222*(abs(GENpTZZ)>75.0&&abs(GENpTZZ)<=80.0);
    k+=1.57896384602*(abs(GENpTZZ)>80.0&&abs(GENpTZZ)<=85.0);
    k+=1.53034682567*(abs(GENpTZZ)>85.0&&abs(GENpTZZ)<=90.0);
    k+=1.56147329708*(abs(GENpTZZ)>90.0&&abs(GENpTZZ)<=95.0);
    k+=1.54468169268*(abs(GENpTZZ)>95.0&&abs(GENpTZZ)<=100.0);
    k+=1.57222952415*(abs(GENpTZZ)>100.0);
  }
  
  if (finalState==2) {
    k+=0.743602533303*(abs(GENpTZZ)>0.0&&abs(GENpTZZ)<=5.0);
    k+=1.14789453219*(abs(GENpTZZ)>5.0&&abs(GENpTZZ)<=10.0);
    k+=1.33815867892*(abs(GENpTZZ)>10.0&&abs(GENpTZZ)<=15.0);
    k+=1.41420044104*(abs(GENpTZZ)>15.0&&abs(GENpTZZ)<=20.0);
    k+=1.45511318916*(abs(GENpTZZ)>20.0&&abs(GENpTZZ)<=25.0);
    k+=1.47569225244*(abs(GENpTZZ)>25.0&&abs(GENpTZZ)<=30.0);
    k+=1.49053003693*(abs(GENpTZZ)>30.0&&abs(GENpTZZ)<=35.0);
    k+=1.50622827695*(abs(GENpTZZ)>35.0&&abs(GENpTZZ)<=40.0);
    k+=1.50328889799*(abs(GENpTZZ)>40.0&&abs(GENpTZZ)<=45.0);
    k+=1.52186945281*(abs(GENpTZZ)>45.0&&abs(GENpTZZ)<=50.0);
    k+=1.52043468754*(abs(GENpTZZ)>50.0&&abs(GENpTZZ)<=55.0);
    k+=1.53977869986*(abs(GENpTZZ)>55.0&&abs(GENpTZZ)<=60.0);
    k+=1.53491994434*(abs(GENpTZZ)>60.0&&abs(GENpTZZ)<=65.0);
    k+=1.51772882172*(abs(GENpTZZ)>65.0&&abs(GENpTZZ)<=70.0);
    k+=1.54494489131*(abs(GENpTZZ)>70.0&&abs(GENpTZZ)<=75.0);
    k+=1.57762411697*(abs(GENpTZZ)>75.0&&abs(GENpTZZ)<=80.0);
    k+=1.55078339014*(abs(GENpTZZ)>80.0&&abs(GENpTZZ)<=85.0);
    k+=1.57078191891*(abs(GENpTZZ)>85.0&&abs(GENpTZZ)<=90.0);
    k+=1.56162666568*(abs(GENpTZZ)>90.0&&abs(GENpTZZ)<=95.0);
    k+=1.54183774627*(abs(GENpTZZ)>95.0&&abs(GENpTZZ)<=100.0);
    k+=1.58485762205*(abs(GENpTZZ)>100.0);
  }
  
  if (k==0.0) return 1.1;
  else return k; // if something goes wrong return inclusive k-factor
  
}


//Variables for input TTree
   UInt_t          Run;
   UInt_t          Event;
   UInt_t          LumiSection;
   Float_t         Avginstlumi;
   Int_t           num_PU_vertices;
   Int_t           PU_BunchCrossing;
   Float_t         MC_weighting;
   Int_t           RECO_nMuHLTMatch;
   Float_t         RECOMU_PT_MuHLTMatch[100];  
   Char_t          HLTPathsFired[20000];
   Float_t         MC_E[7];
   Float_t         MC_PT[7];
   Float_t         MC_ETA[7];
   Float_t         MC_THETA[7];
   Float_t         MC_PHI[7];
   Float_t         MC_MASS[7];
   Float_t         MC_PDGID[7];
   Float_t         MC_LEPT_PT[4];
   Float_t         MC_LEPT_ETA[4];
   Float_t         MC_LEPT_PHI[4];
   Float_t         MC_LEPT_THETA[4];
   Float_t         MC_LEPT_PDGID[4];
   Float_t         MC_Z_PT[2][5];
   Float_t         MC_Z_ETA[2][5];
   Float_t         MC_Z_PHI[2][5];
   Float_t         MC_Z_THETA[2][5];
   Float_t         MC_Z_MASS[2][5];
   Float_t         MC_Z_PDGID[2][5];
   Float_t         MC_fourl_MASS[50][5];
   Float_t         MC_fourl_PT[50][5];
   Float_t         MC_fourl_PDGID[50][5];
   Float_t         MC_ZZ_MASS[4][7];
   Float_t         MC_ZZ_PT[4][7];
   Float_t         MC_ZZ_ETA[4][7];
   Float_t         MC_ZZ_PHI[4][7];
   Float_t         MC_ZZ_THETA[4][7];
   Float_t         MC_ZZ_PDGID[4][7];
   Float_t         MC_GENJET_PT[4];
   Float_t         MC_GENJET_ETA[4];
   Float_t         MC_GENJET_PHI[4];
   Float_t         MC_GENMET;
   Double_t        RECORF_2e2mu_cosTheta1_spin[100];
   Double_t        RECORF_2e2mu_cosTheta2_spin[100];
   Double_t        RECORF_2e2mu_cosThetaStar_spin[100];
   Double_t        RECORF_2e2mu_Phi_spin[100];
   Double_t        RECORF_2e2mu_Phi1_spin[100];
   Double_t        RECORF_2e2mu_Phi2_spin[100];
   Double_t        RECORF_2e2mu_phi1RF_spin[100];
   Double_t        RECORF_2e2mu_phi2RF_spin[100];
   Double_t        RECORF_2e2mu_MELA[100];
   Double_t        RECORF_4e_cosTheta1_spin[100];
   Double_t        RECORF_4e_cosTheta2_spin[100];
   Double_t        RECORF_4e_cosThetaStar_spin[100];
   Double_t        RECORF_4e_Phi_spin[100];
   Double_t        RECORF_4e_Phi1_spin[100];
   Double_t        RECORF_4e_Phi2_spin[100];
   Double_t        RECORF_4e_phi1RF_spin[100];
   Double_t        RECORF_4e_phi2RF_spin[100];
   Double_t        RECORF_4e_MELA[100];
   Double_t        RECORF_4mu_cosTheta1_spin[100];
   Double_t        RECORF_4mu_cosTheta2_spin[100];
   Double_t        RECORF_4mu_cosThetaStar_spin[100];
   Double_t        RECORF_4mu_Phi_spin[100];
   Double_t        RECORF_4mu_Phi1_spin[100];
   Double_t        RECORF_4mu_Phi2_spin[100];
   Double_t        RECORF_4mu_phi1RF_spin[100];
   Double_t        RECORF_4mu_phi2RF_spin[100];
   Double_t        RECORF_4mu_MELA[100];
   Float_t         RECO_ZMM_MASS[50];
   Float_t         RECO_ZEE_MASS[50];
   Float_t         RECO_DiLep_MASS[50];
   Float_t         RECO_ZMM_PT[3][50];
   Float_t         RECO_ZEE_PT[3][50];
   Float_t         RECO_DiLep_PT[3][50];
   Float_t         RECO_ZMM_ETA[3][50];
   Float_t         RECO_ZEE_ETA[3][50];
   Float_t         RECO_DiLep_ETA[3][50];
   Float_t         RECO_ZMM_PHI[3][50];
   Float_t         RECO_ZEE_PHI[3][50];
   Float_t         RECO_DiLep_PHI[3][50];
   Float_t         RECO_ZMMss_MASS[50];
   Float_t         RECO_ZEEss_MASS[50];
   Float_t         RECO_ZEM_MASS[50];
   Float_t         RECO_ZMMss_PT[3][50];
   Float_t         RECO_ZEEss_PT[3][50];
   Float_t         RECO_ZEM_PT[3][50];
   Float_t         RECO_ZMMss_ETA[3][50];
   Float_t         RECO_ZEEss_ETA[3][50];
   Float_t         RECO_ZEM_ETA[3][50];
   Float_t         RECO_ZMMss_PHI[3][50];
   Float_t         RECO_ZEEss_PHI[3][50];
   Float_t         RECO_ZEM_PHI[3][50];
   Float_t         RECO_MMMM_MASS[7][100];
   Float_t         RECO_MMMM_PT[7][100];
   Float_t         RECO_MMMM_ETA[7][100];
   Float_t         RECO_MMMM_PHI[7][100];
   Float_t         RECO_MMMM_MASS_REFIT[100];
   Float_t         RECO_EEEE_MASS[7][100];
   Float_t         RECO_EEEE_PT[7][100];
   Float_t         RECO_EEEE_ETA[7][100];
   Float_t         RECO_EEEE_PHI[7][100];
   Float_t         RECO_EEEE_MASS_REFIT[100];
   Float_t         RECO_EEMM_MASS[7][100];
   Float_t         RECO_EEMM_PT[7][100];
   Float_t         RECO_EEMM_ETA[7][100];
   Float_t         RECO_EEMM_PHI[7][100];
   Float_t         RECO_EEMM_MASS_REFIT[100];
   Float_t         RECO_LLL0_MASS[50];
   Float_t         RECO_LLL1_MASS[50];
   Float_t         RECO_LLL2_MASS[50];
   Float_t         RECO_LLL3_MASS[50];
   Float_t         RECO_LLL0_PT[4][50];
   Float_t         RECO_LLL1_PT[4][50];
   Float_t         RECO_LLL2_PT[4][50];
   Float_t         RECO_LLL3_PT[4][50];
   Float_t         RECO_LLLl0_MASS[20];
   Float_t         RECO_LLLl1_MASS[20];
   Float_t         RECO_LLLl0_PT[5][20];
   Float_t         RECO_LLLl1_PT[5][20];
   Float_t         RECO_LLLL0ss_MASS[20];
   Float_t         RECO_LLLL0ss_PT[5][20];
   Float_t         RECO_LLLL1ss_MASS[20];
   Float_t         RECO_LLLL1ss_PT[5][20];
   Float_t         RECO_LLLL2ss_MASS[20];
   Float_t         RECO_LLLL2ss_PT[5][20];
   Float_t         RECO_LLLL_MASS[7][100];
   Float_t         RECO_LLLL_PT[7][100];
   Float_t         RECO_LLLL_ETA[7][100];
   Float_t         RECO_LLLL_PHI[7][100];
   Float_t         RECOELE_E[100];
   Float_t         RECOELE_PT[100];
   Float_t         RECOELE_PTError[100];
   Float_t         RECOELE_P[100];
   Float_t         RECOELE_ETA[100];
   Float_t         RECOELE_THETA[100];
   Float_t         RECOELE_PHI[100];
   Float_t         RECOELE_MASS[100];
   Float_t         RECOELE_CHARGE[100];
   UChar_t         RECOELE_isEcalDriven[100];
   UChar_t         RECOELE_isTrackerDriven[100];
   Float_t         RECOELE_gsftrack_NPixHits[100];
   Float_t         RECOELE_gsftrack_NStripHits[100];
   Float_t         RECOELE_gsftrack_chi2[100];
   Float_t         RECOELE_gsftrack_dxyB[100];
   Float_t         RECOELE_gsftrack_dxy[100];
   Float_t         RECOELE_gsftrack_dxyError[100];
   Float_t         RECOELE_gsftrack_dzB[100];
   Float_t         RECOELE_gsftrack_dz[100];
   Float_t         RECOELE_gsftrack_dzError[100];
   Int_t           RECOELE_gsftrack_losthits[100];
   Int_t           RECOELE_gsftrack_validhits[100];
   Int_t           RECOELE_gsftrack_expected_inner_hits[100];
   Float_t         RECOELE_scl_E[100];
   Float_t         RECOELE_scl_Et[100];
   Float_t         RECOELE_scl_Eta[100];
   Float_t         RECOELE_scl_Phi[100];
   Float_t         RECOELE_ep[100];
   Float_t         RECOELE_eSeedp[100];
   Float_t         RECOELE_eSeedpout[100];
   Float_t         RECOELE_eElepout[100];
   Float_t         RECOELE_deltaEtaIn[100];
   Float_t         RECOELE_deltaEtaSeed[100];
   Float_t         RECOELE_deltaEtaEle[100];
   Float_t         RECOELE_deltaPhiIn[100];
   Float_t         RECOELE_deltaPhiSeed[100];
   Float_t         RECOELE_deltaPhiEle[100];
   Int_t           RECOELE_isbarrel[100];
   Int_t           RECOELE_isendcap[100];
   Int_t           RECOELE_isGap[100];
   Int_t           RECOELE_isEBetaGap[100];
   Int_t           RECOELE_isEBphiGap[100];
   Int_t           RECOELE_isEEdeeGap[100];
   Int_t           RECOELE_isEEringGap[100];
   Float_t         RECOELE_sigmaIetaIeta[100];
   Float_t         RECOELE_sigmaEtaEta[100];
   Float_t         RECOELE_e15[100];
   Float_t         RECOELE_e25max[100];
   Float_t         RECOELE_e55[100];
   Float_t         RECOELE_he[100];
   Float_t         RECOELE_r9[100];
   Float_t         RECOELE_mva[100];
   Float_t         RECOELE_fbrem[100];
   Int_t           RECOELE_nbrems[100];
   Int_t           RECOELE_Class[100];
   Double_t        RECOELE_fbrem_mode[100];
   Double_t        RECOELE_fbrem_mean[100];
   Float_t         RECOELE_EGMTRACKISO[100];
   Float_t         RECOELE_EGMHCALISO[100];
   Float_t         RECOELE_EGMECALISO[100];
   Float_t         RECOELE_EGMX[100];
   Double_t        RECOELE_PFchAllPart[100];
   Double_t        RECOELE_PFchHad[100];
   Double_t        RECOELE_PFneuHad[100];
   Double_t        RECOELE_PFphoton[100];
   Double_t        RECOELE_PFPUchAllPart[100];
   Double_t        RECOELE_PFX_dB[100];
   Double_t        RECOELE_PFX_rho[100];
   Double_t        RECOELE_PFX_rho_new[100];
   Double_t        RECOELE_regEnergy[100];
   Double_t        RECOELE_regEnergyError[100];
   Float_t         RECOELE_SIP[100];
   Float_t         RECOELE_IP[100];
   Float_t         RECOELE_IPERROR[100];
   Float_t         RECOELE_SIP_KF[100];
   Float_t         RECOELE_IP_KF[100];
   Float_t         RECOELE_IPERROR_KF[100];
   Float_t         RECOELE_SIP_GD[100];
   Float_t         RECOELE_SIP_GDEEEE[100];
   Float_t         RECOELE_SIP_Std[100];
   Float_t         RECOELE_SIP_StdEEEE[100];
   Float_t         RECOELE_SIP_Kin[100];
   Float_t         RECOELE_SIP_KinEEEE[100];
   Float_t         RECOELE_STIP[100];
   Float_t         RECOELE_SLIP[100];
   Float_t         RECOELE_TIP[100];
   Float_t         RECOELE_LIP[100];
   Float_t         RECOELE_TIPERROR[100];
   Float_t         RECOELE_LIPERROR[100];
   Double_t        RECOELE_sclRawE[100];
   Double_t        RECOELE_sclX[100];
   Double_t        RECOELE_sclY[100];
   Double_t        RECOELE_sclZ[100];
   Double_t        RECOELE_seedSubdet1[100];
   Double_t        RECOELE_seedDphi1[100];
   Double_t        RECOELE_seedDrz1[100];
   Double_t        RECOELE_seedSubdet2[100];
   Double_t        RECOELE_seedDphi2[100];
   Double_t        RECOELE_seedDrz2[100];
   Double_t        RECOELE_eidVeryLoose[100];
   Double_t        RECOELE_eidLoose[100];
   Double_t        RECOELE_eidMedium[100];
   Double_t        RECOELE_eidTight[100];
   Double_t        RECOELE_eidHZZVeryLoose[100];
   Double_t        RECOELE_eidHZZLoose[100];
   Double_t        RECOELE_eidHZZMedium[100];
   Double_t        RECOELE_eidHZZTight[100];
   Double_t        RECOELE_mvaTrigV0[100];
   Double_t        RECOELE_mvaNonTrigV0[100];
   Double_t        RECOELE_COV[100][3][3];
   UChar_t         RECOMU_isPFMu[100];
   UChar_t         RECOMU_isGlobalMu[100];
   UChar_t         RECOMU_isStandAloneMu[100];
   UChar_t         RECOMU_isTrackerMu[100];
   UChar_t         RECOMU_isCaloMu[100];
   UChar_t         RECOMU_isTrackerHighPtMu[100];
   Float_t         RECOMU_E[100];
   Float_t         RECOMU_PT[100];
   Float_t RECOMU_mubesttrkPTError[100];
   Float_t         RECOMU_P[100];
   Float_t         RECOMU_ETA[100];
   Float_t         RECOMU_THETA[100];
   Float_t         RECOMU_PHI[100];
   Float_t         RECOMU_MASS[100];
   Float_t         RECOMU_CHARGE[100];
   Double_t        RECOMU_COV[100][3][3];
   Float_t         RECOMU_TRACKISO[100];
   Float_t         RECOMU_TRACKISO_SUMPT[100];
   Float_t         RECOMU_HCALISO[100];
   Float_t         RECOMU_ECALISO[100];
   Float_t         RECOMU_X[100];
   Double_t        RECOMU_PFchHad[100];
   Double_t        RECOMU_PFneuHad[100];
   Double_t        RECOMU_PFphoton[100];
   Double_t        RECOMU_PFPUchAllPart[100];
   Double_t        RECOMU_PFX_dB[100];
   Double_t        RECOMU_PFX_dB_new[100];
   Double_t        RECOMU_PFX_rho[100];
   Double_t        RECOPFPHOT_PFchHad[20];
   Double_t        RECOPFPHOT_PFneuHad[20];
   Double_t        RECOPFPHOT_PFphoton[20];
   Double_t        RECOPFPHOT_PFPUchAllPart[20];
   Double_t        RECOPFPHOT_PFX_rho[20];
   Float_t         RECOMU_SIP[100];
   Float_t         RECOMU_IP[100];
   Float_t         RECOMU_IPERROR[100];
   Float_t         RECOMU_SIP_KF[100];
   Float_t         RECOMU_IP_KF[100];
   Float_t         RECOMU_IPERROR_KF[100];
   Float_t         RECOMU_SIP_GD[100];
   Float_t         RECOMU_SIP_GDMMMM[100];
   Float_t         RECOMU_SIP_Std[100];
   Float_t         RECOMU_SIP_StdMMMM[100];
   Float_t         RECOMU_SIP_Kin[100];
   Float_t         RECOMU_SIP_KinMMMM[100];
   Float_t         RECOMU_STIP[100];
   Float_t         RECOMU_SLIP[100];
   Float_t         RECOMU_TIP[100];
   Float_t         RECOMU_LIP[100];
   Float_t         RECOMU_TIPERROR[100];
   Float_t         RECOMU_LIPERROR[100];
   Float_t         RECOMU_caloCompatibility[100];
   Float_t         RECOMU_segmentCompatibility[100];
   UInt_t          RECOMU_numberOfMatches[100];
   UInt_t          RECOMU_numberOfMatchedStations[100];
   UChar_t         RECOMU_glbmuPromptTight[100];
   UChar_t         RECOMU_trkmuArbitration[100];
   UChar_t         RECOMU_trkmu2DCompatibilityLoose[100];
   UChar_t         RECOMU_trkmu2DCompatibilityTight[100];
   UChar_t         RECOMU_trkmuOneStationLoose[100];
   UChar_t         RECOMU_trkmuOneStationTight[100];
   UChar_t         RECOMU_trkmuLastStationLoose[100];
   UChar_t         RECOMU_trkmuLastStationTight[100];
   UChar_t         RECOMU_trkmuOneStationAngLoose[100];
   UChar_t         RECOMU_trkmuOneStationAngTight[100];
   UChar_t         RECOMU_trkmuLastStationAngLoose[100];
   UChar_t         RECOMU_trkmuLastStationAngTight[100];
   UChar_t         RECOMU_trkmuLastStationOptimizedLowPtLoose[100];
   UChar_t         RECOMU_trkmuLastStationOptimizedLowPtTight[100];
   Float_t         RECOMU_mutrkPT[100];
   Float_t         RECOMU_mutrkPTError[100];
   Float_t         RECOMU_mutrkDxy[100];
   Float_t         RECOMU_mutrkDxyError[100];
   Float_t         RECOMU_mutrkDxyB[100];
   Float_t         RECOMU_mutrkDz[100];
   Float_t         RECOMU_mutrkDzError[100];
   Float_t         RECOMU_mutrkDzB[100];
   Float_t         RECOMU_mutrkChi2PerNdof[100];
   Float_t         RECOMU_mutrkCharge[100];
   Float_t         RECOMU_mutrkNHits[100];
   Float_t         RECOMU_mutrkNStripHits[100];
   Float_t         RECOMU_mutrkNPixHits[100];
   Float_t         RECOMU_mutrkNMuonHits[100];
   Float_t         RECOMU_mutrktrackerLayersWithMeasurement[100];
   Float_t         RECOMU_muInnertrkDxy[100];
   Float_t         RECOMU_muInnertrkDxyError[100];
   Float_t         RECOMU_muInnertrkDxyB[100];
   Float_t         RECOMU_muInnertrkDz[100];
   Float_t         RECOMU_muInnertrkDzError[100];
   Float_t         RECOMU_muInnertrkDzB[100];
   Float_t         RECOMU_muInnertrkChi2PerNdof[100];
   Float_t         RECOMU_muInnertrktrackerLayersWithMeasurement[100];
   Float_t         RECOMU_muInnertrkPT[100];
   Float_t         RECOMU_muInnertrkPTError[100];
   Float_t         RECOMU_muInnertrkCharge[100];
   Float_t         RECOMU_muInnertrkNHits[100];
   Float_t         RECOMU_muInnertrkNStripHits[100];
   Float_t         RECOMU_muInnertrkNPixHits[100];
   Int_t           RECOMU_mubesttrkType[100];
   Float_t         RECOMU_mubesttrkDxy[100];
   Float_t         RECOMU_mubesttrkDxyError[100];
   Float_t         RECOMU_mubesttrkDz[100];
   Float_t         RECOMU_mubesttrkDzError[100];
   Double_t        ftsigma[100];
   Double_t        gdX[100];
   Double_t        gdY[100];
   Double_t        gdZ[100];
   Double_t        ftsigmalag[100];
   Double_t        gdlagX[100];
   Double_t        gdlagY[100];
   Double_t        gdlagZ[100];
   Double_t        gdlagProb[100];
   Double_t        gdlagNdof[100];
   Double_t        ftsigmaMMMM[100];
   Double_t        gdXMMMM[100];
   Double_t        gdYMMMM[100];
   Double_t        gdZMMMM[100];
   Double_t        ftsigmalagMMMM[100];
   Double_t        gdlagXMMMM[100];
   Double_t        gdlagYMMMM[100];
   Double_t        gdlagZMMMM[100];
   Double_t        gdlagProbMMMM[100];
   Double_t        gdlagNdofMMMM[100];
   Double_t        ftsigmaEEEE[100];
   Double_t        gdXEEEE[100];
   Double_t        gdYEEEE[100];
   Double_t        gdZEEEE[100];
   Double_t        ftsigmalagEEEE[100];
   Double_t        gdlagXEEEE[100];
   Double_t        gdlagYEEEE[100];
   Double_t        gdlagZEEEE[100];
   Double_t        gdlagProbEEEE[100];
   Double_t        gdlagNdofEEEE[100];
   Double_t        StdFitVertexX[100];
   Double_t        StdFitVertexY[100];
   Double_t        StdFitVertexZ[100];
   Double_t        StdFitVertexChi2r[100];
   Double_t        StdFitVertexProb[100];
   Float_t         StdFitVertexTrack_PT[4][100];
   Float_t         StdFitVertexTrack_ETA[4][100];
   Float_t         StdFitVertexTrack_PHI[4][100];
   Double_t        KinFitVertexX[100];
   Double_t        KinFitVertexY[100];
   Double_t        KinFitVertexZ[100];
   Double_t        KinFitVertexChi2r[100];
   Double_t        KinFitVertexProb[100];
   Double_t        StdFitVertexXMMMM[100];
   Double_t        StdFitVertexYMMMM[100];
   Double_t        StdFitVertexZMMMM[100];
   Double_t        StdFitVertexChi2rMMMM[100];
   Double_t        StdFitVertexProbMMMM[100];
   Float_t         StdFitVertexTrackMMMM_PT[4][100];
   Float_t         StdFitVertexTrackMMMM_ETA[4][100];
   Float_t         StdFitVertexTrackMMMM_PHI[4][100];
   Double_t        KinFitVertexXMMMM[100];
   Double_t        KinFitVertexYMMMM[100];
   Double_t        KinFitVertexZMMMM[100];
   Double_t        KinFitVertexChi2rMMMM[100];
   Double_t        KinFitVertexProbMMMM[100];
   Double_t        StdFitVertexXEEEE[100];
   Double_t        StdFitVertexYEEEE[100];
   Double_t        StdFitVertexZEEEE[100];
   Double_t        StdFitVertexChi2rEEEE[100];
   Double_t        StdFitVertexProbEEEE[100];
   Float_t         StdFitVertexTrackEEEE_PT[4][100];
   Float_t         StdFitVertexTrackEEEE_ETA[4][100];
   Float_t         StdFitVertexTrackEEEE_PHI[4][100];
   Double_t        KinFitVertexXEEEE[100];
   Double_t        KinFitVertexYEEEE[100];
   Double_t        KinFitVertexZEEEE[100];
   Double_t        KinFitVertexChi2rEEEE[100];
   Double_t        KinFitVertexProbEEEE[100];
   Double_t        StdFitVertexChi2rMMM[50];
   Double_t        StdFitVertexProbMMM[50];
   Double_t        StdFitVertexChi2rMME[50];
   Double_t        StdFitVertexProbMME[50];
   Double_t        StdFitVertexChi2rEEE[50];
   Double_t        StdFitVertexProbEEE[50];
   Double_t        StdFitVertexChi2rMEE[50];
   Double_t        StdFitVertexProbMEE[50];
   Double_t        StdFitVertexChi2rDiLep[40];
   Double_t        StdFitVertexProbDiLep[40];
   Float_t         ConvMapDist[100];
   Float_t         ConvMapDcot[100];
   UChar_t         RECOMU_MatchingMCTruth[100];
   Float_t         RECOMU_MatchingMCpT[100];
   Float_t         RECOMU_MatchingMCEta[100];
   Float_t         RECOMU_MatchingMCPhi[100];
   UChar_t         RECOELE_MatchingMCTruth[100];
   Float_t         RECOELE_MatchingMCpT[100];
   Float_t         RECOELE_MatchingMCEta[100];
   Float_t         RECOELE_MatchingMCPhi[100];
   UChar_t         RECOPHOT_MatchingMCTruth[50];
   Float_t         RECOPHOT_MatchingMCpT[50];
   Float_t         RECOPHOT_MatchingMCEta[50];
   Float_t         RECOPHOT_MatchingMCPhi[50];
   UChar_t         RECOzMuMu_MatchingMCTruth[50];
   Float_t         RECOzMuMu_MatchingMCpT[50];
   Float_t         RECOzMuMu_MatchingMCmass[50];
   Float_t         RECOzMuMu_MatchingMCEta[50];
   Float_t         RECOzMuMu_MatchingMCPhi[50];
   UChar_t         RECOzEE_MatchingMCTruth[50];
   Float_t         RECOzEE_MatchingMCpT[50];
   Float_t         RECOzEE_MatchingMCmass[50];
   Float_t         RECOzEE_MatchingMCEta[50];
   Float_t         RECOzEE_MatchingMCPhi[50];
   UChar_t         RECOHzzEEEE_MatchingMCTruth[100];
   Float_t         RECOHzzEEEE_MatchingMCpT[100];
   Float_t         RECOHzzEEEE_MatchingMCmass[100];
   Float_t         RECOHzzEEEE_MatchingMCEta[100];
   Float_t         RECOHzzEEEE_MatchingMCPhi[100];
   UChar_t         RECOHzzEEMM_MatchingMCTruth[100];
   Float_t         RECOHzzEEMM_MatchingMCpT[100];
   Float_t         RECOHzzEEMM_MatchingMCmass[100];
   Float_t         RECOHzzEEMM_MatchingMCEta[100];
   Float_t         RECOHzzEEMM_MatchingMCPhi[100];
   UChar_t         RECOHzzMMMM_MatchingMCTruth[100];
   Float_t         RECOHzzMMMM_MatchingMCpT[100];
   Float_t         RECOHzzMMMM_MatchingMCmass[100];
   Float_t         RECOHzzMMMM_MatchingMCEta[100];
   Float_t         RECOHzzMMMM_MatchingMCPhi[100];
   Int_t           RECO_NMU;
   Int_t           RECO_NELE;
   Int_t           RECO_NTRACK;
   Float_t         RECO_TRACK_PT[200];
   Float_t         RECO_TRACK_ETA[200];
   Float_t         RECO_TRACK_PHI[200];
   Float_t         RECO_TRACK_CHI2[200];
   Float_t         RECO_TRACK_CHI2RED[200];
   Float_t         RECO_TRACK_CHI2PROB[200];
   Int_t           RECO_TRACK_NHITS[200];
   Float_t         RECO_TRACK_DXY[200];
   Float_t         RECO_TRACK_DXYERR[200];
   Float_t         RECO_TRACK_DZ[200];
   Float_t         RECO_TRACK_DZERR[200];
   Int_t           RECO_NPHOT;
   Float_t         RECOPHOT_PT[20];
   Float_t         RECOPHOT_ETA[20];
   Float_t         RECOPHOT_PHI[20];
   Float_t         RECOPHOT_THETA[20];
   Int_t           RECO_NPFPHOT;
   Float_t         RECOPFPHOT_PT[20];
   Float_t         RECOPFPHOT_PTError[20];
   Float_t         RECOPFPHOT_ETA[20];
   Float_t         RECOPFPHOT_PHI[20];
   Float_t         RECOPFPHOT_THETA[20];
   Double_t        BeamSpot_X;
   Double_t        BeamSpot_Y;
   Double_t        BeamSpot_Z;
   Int_t           RECO_NVTX;
   Float_t         RECO_VERTEX_x[15];
   Float_t         RECO_VERTEX_y[15];
   Float_t         RECO_VERTEX_z[15];
   Float_t         RECO_VERTEX_ndof[15];
   Float_t         RECO_VERTEX_chi2[15];
   Int_t           RECO_VERTEX_ntracks[15];
   Float_t         RECO_VERTEXPROB[15];
   UChar_t         RECO_VERTEX_isValid[15];
   Float_t         RECO_VERTEX_TRACK_PT[15][100];
   Int_t           RECO_PFJET_N;
   Int_t           RECO_PFJET_CHARGE[100];
   Float_t         RECO_PFJET_ET[100];
   Float_t         RECO_PFJET_PT[100];
   Float_t         RECO_PFJET_ETA[100];
   Float_t         RECO_PFJET_PHI[100];
   Int_t           RECO_PFJET_PUID[100];
   Float_t         RECO_PFJET_PUID_MVA[100];
   Double_t        RHO_ele;
   Double_t        RHO_mu;
   Float_t         RECO_CALOMET;
   Float_t         RECO_PFMET;
   Float_t         RECO_PFMET_X;
   Float_t         RECO_PFMET_Y;
   Float_t         RECO_PFMET_PHI;
   Float_t         RECO_PFMET_THETA;
   Float_t         RECO_PFMETT1;
   Float_t         RECO_PFMETT1_X;
   Float_t         RECO_PFMETT1_Y;
   Float_t         RECO_PFMETT1_PHI;
   Float_t         RECO_PFMETT1_THETA;
   Float_t         RECO_TCMET;
   Float_t         RECO_CORMETMUONS;
   Float_t         tCHighEff_BTagJet_PT[50];
   Float_t         tCHighPur_BTagJet_PT[50];
   Float_t         cSV_BTagJet_PT[50];
   Float_t         tCHighEff_BTagJet_ETA[50];
   Float_t         tCHighPur_BTagJet_ETA[50];
   Float_t         cSV_BTagJet_ETA[50];
   Float_t         tCHighEff_BTagJet_PHI[50];
   Float_t         tCHighPur_BTagJet_PHI[50];
   Float_t         cSV_BTagJet_PHI[50];
   Float_t         tCHighEff_BTagJet_DISCR[50];
   Float_t         tCHighPur_BTagJet_DISCR[50];
   Float_t         cSV_BTagJet_DISCR[50];
   Float_t         cSV_BTagJet_ET[50];
   Float_t         RECO_PFJET_PT_UncDn[200];
   Float_t         RECO_PFMET_JetEnUp;
   Float_t         RECO_PFMET_JetEnDn;
   Float_t         RECO_PFMET_ElectronEnUp;
   Float_t         RECO_PFMET_ElectronEnDn;
   Float_t         RECO_PFMET_MuonEnUp;
   Float_t         RECO_PFMET_MuonEnDn;
   Float_t         RECO_PFMET_JetResUp;
   Float_t         RECO_PFMET_JetResDn;
   Float_t         RECO_PFMET_UnclusteredEnUp;
   Float_t         RECO_PFMET_UnclusteredEnDn;
   Float_t         RECO_PFMET_PhotonEnUp;
   Float_t         RECO_PFMET_PhotonEnDn;



///===============================================================================================================
///Note: When running on data - Xsection=1, neventsPreHLT=1, neventsPostHLT=1, DATA_type=2016, MC_type=NO
///Note: When running on data - Xsection=<>pb, neventsPreHLT=<>, neventsPostHLT=<>, DATA_type=NO, MC_type=Spring16
///===============================================================================================================   
void ComputeFRandCRsOS(TString dataset_path, TString datasetName, TString store_path, double Xsection, int neventsPreHLT,
		     int neventsPostHLT, TString DATA_type, TString MC_type, int Nevents){
  
  cout<< "\nDatasetPath: "<< dataset_path
      << "\nDatasetName: "<< datasetName
      << "\nStorePath: "<< store_path
      << "\nXsection: "<< Xsection
      << "\nneventsPreHLT: " <<neventsPreHLT
      << "\nneventsPostHLT: " <<neventsPostHLT
      << "\nDATA_type: " << DATA_type
      << "\nMC_type: " << MC_type
      << "\nNevents: " << Nevents
  << endl;
  
  TString infile_name = dataset_path + "roottree_leptons_crab_" + datasetName + ".root";
  TFile *root_file = new TFile(infile_name);
  TTree *fChain = (TTree*)root_file->Get("HZZ4LeptonsAnalysis");
   fChain->SetBranchAddress("Run", &Run);
   fChain->SetBranchAddress("Event", &Event);
   fChain->SetBranchAddress("LumiSection", &LumiSection);
   fChain->SetBranchAddress("Avginstlumi", &Avginstlumi);
   fChain->SetBranchAddress("num_PU_vertices", &num_PU_vertices);
   fChain->SetBranchAddress("PU_BunchCrossing", &PU_BunchCrossing);
   fChain->SetBranchAddress("MC_weighting", &MC_weighting);
   fChain->SetBranchAddress("RECO_nMuHLTMatch", &RECO_nMuHLTMatch);
   fChain->SetBranchAddress("RECOMU_PT_MuHLTMatch", &RECOMU_PT_MuHLTMatch);
   fChain->SetBranchAddress("HLTPathsFired", &HLTPathsFired);
   fChain->SetBranchAddress("MC_E", &MC_E);
   fChain->SetBranchAddress("MC_PT", &MC_PT);
   fChain->SetBranchAddress("MC_ETA", &MC_ETA);
   fChain->SetBranchAddress("MC_THETA", &MC_THETA);
   fChain->SetBranchAddress("MC_PHI", &MC_PHI);
   fChain->SetBranchAddress("MC_MASS", &MC_MASS);
   fChain->SetBranchAddress("MC_PDGID", &MC_PDGID);
   fChain->SetBranchAddress("MC_LEPT_PT", &MC_LEPT_PT);
   fChain->SetBranchAddress("MC_LEPT_ETA", &MC_LEPT_ETA);
   fChain->SetBranchAddress("MC_LEPT_PHI", &MC_LEPT_PHI);
   fChain->SetBranchAddress("MC_LEPT_THETA", &MC_LEPT_THETA);
   fChain->SetBranchAddress("MC_LEPT_PDGID", &MC_LEPT_PDGID);
   fChain->SetBranchAddress("MC_Z_PT", &MC_Z_PT);
   fChain->SetBranchAddress("MC_Z_ETA", &MC_Z_ETA);
   fChain->SetBranchAddress("MC_Z_PHI", &MC_Z_PHI);
   fChain->SetBranchAddress("MC_Z_THETA", &MC_Z_THETA);
   fChain->SetBranchAddress("MC_Z_MASS", &MC_Z_MASS);
   fChain->SetBranchAddress("MC_Z_PDGID", &MC_Z_PDGID);
   fChain->SetBranchAddress("MC_fourl_MASS", &MC_fourl_MASS);
   fChain->SetBranchAddress("MC_fourl_PT", &MC_fourl_PT);
   fChain->SetBranchAddress("MC_fourl_PDGID", &MC_fourl_PDGID);
   fChain->SetBranchAddress("MC_ZZ_MASS", &MC_ZZ_MASS);
   fChain->SetBranchAddress("MC_ZZ_PT", &MC_ZZ_PT);
   fChain->SetBranchAddress("MC_ZZ_ETA", &MC_ZZ_ETA);
   fChain->SetBranchAddress("MC_ZZ_PHI", &MC_ZZ_PHI);
   fChain->SetBranchAddress("MC_ZZ_THETA", &MC_ZZ_THETA);
   fChain->SetBranchAddress("MC_ZZ_PDGID", &MC_ZZ_PDGID);
   fChain->SetBranchAddress("MC_GENJET_PT", &MC_GENJET_PT);
   fChain->SetBranchAddress("MC_GENJET_ETA", &MC_GENJET_ETA);
   fChain->SetBranchAddress("MC_GENJET_PHI", &MC_GENJET_PHI);       
   fChain->SetBranchAddress("MC_GENMET", &MC_GENMET);
   fChain->SetBranchAddress("RECORF_2e2mu_cosTheta1_spin", &RECORF_2e2mu_cosTheta1_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_cosTheta2_spin", &RECORF_2e2mu_cosTheta2_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_cosThetaStar_spin", &RECORF_2e2mu_cosThetaStar_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_Phi_spin", &RECORF_2e2mu_Phi_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_Phi1_spin", &RECORF_2e2mu_Phi1_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_Phi2_spin", &RECORF_2e2mu_Phi2_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_phi1RF_spin", &RECORF_2e2mu_phi1RF_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_phi2RF_spin", &RECORF_2e2mu_phi2RF_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_MELA", &RECORF_2e2mu_MELA);
   fChain->SetBranchAddress("RECORF_4e_cosTheta1_spin", &RECORF_4e_cosTheta1_spin);
   fChain->SetBranchAddress("RECORF_4e_cosTheta2_spin", &RECORF_4e_cosTheta2_spin);
   fChain->SetBranchAddress("RECORF_4e_cosThetaStar_spin", &RECORF_4e_cosThetaStar_spin);
   fChain->SetBranchAddress("RECORF_4e_Phi_spin", &RECORF_4e_Phi_spin);
   fChain->SetBranchAddress("RECORF_4e_Phi1_spin", &RECORF_4e_Phi1_spin);
   fChain->SetBranchAddress("RECORF_4e_Phi2_spin", &RECORF_4e_Phi2_spin);
   fChain->SetBranchAddress("RECORF_4e_phi1RF_spin", &RECORF_4e_phi1RF_spin);
   fChain->SetBranchAddress("RECORF_4e_phi2RF_spin", &RECORF_4e_phi2RF_spin);
   fChain->SetBranchAddress("RECORF_4e_MELA", &RECORF_4e_MELA);
   fChain->SetBranchAddress("RECORF_4mu_cosTheta1_spin", &RECORF_4mu_cosTheta1_spin);
   fChain->SetBranchAddress("RECORF_4mu_cosTheta2_spin", &RECORF_4mu_cosTheta2_spin);
   fChain->SetBranchAddress("RECORF_4mu_cosThetaStar_spin", &RECORF_4mu_cosThetaStar_spin);
   fChain->SetBranchAddress("RECORF_4mu_Phi_spin", &RECORF_4mu_Phi_spin);
   fChain->SetBranchAddress("RECORF_4mu_Phi1_spin", &RECORF_4mu_Phi1_spin);
   fChain->SetBranchAddress("RECORF_4mu_Phi2_spin", &RECORF_4mu_Phi2_spin);
   fChain->SetBranchAddress("RECORF_4mu_phi1RF_spin", &RECORF_4mu_phi1RF_spin);
   fChain->SetBranchAddress("RECORF_4mu_phi2RF_spin", &RECORF_4mu_phi2RF_spin);
   fChain->SetBranchAddress("RECORF_4mu_MELA", &RECORF_4mu_MELA);
   fChain->SetBranchAddress("RECO_ZMM_MASS", &RECO_ZMM_MASS);
   fChain->SetBranchAddress("RECO_ZEE_MASS", &RECO_ZEE_MASS);
   fChain->SetBranchAddress("RECO_DiLep_MASS", &RECO_DiLep_MASS);
   fChain->SetBranchAddress("RECO_ZMM_PT", &RECO_ZMM_PT);
   fChain->SetBranchAddress("RECO_ZEE_PT", &RECO_ZEE_PT);
   fChain->SetBranchAddress("RECO_DiLep_PT", &RECO_DiLep_PT);
   fChain->SetBranchAddress("RECO_ZMM_ETA", &RECO_ZMM_ETA);
   fChain->SetBranchAddress("RECO_ZEE_ETA", &RECO_ZEE_ETA);
   fChain->SetBranchAddress("RECO_DiLep_ETA", &RECO_DiLep_ETA);
   fChain->SetBranchAddress("RECO_ZMM_PHI", &RECO_ZMM_PHI);
   fChain->SetBranchAddress("RECO_ZEE_PHI", &RECO_ZEE_PHI);
   fChain->SetBranchAddress("RECO_DiLep_PHI", &RECO_DiLep_PHI);
   fChain->SetBranchAddress("RECO_ZMMss_MASS", &RECO_ZMMss_MASS);
   fChain->SetBranchAddress("RECO_ZEEss_MASS", &RECO_ZEEss_MASS);
   fChain->SetBranchAddress("RECO_ZEM_MASS", &RECO_ZEM_MASS);
   fChain->SetBranchAddress("RECO_ZMMss_PT", &RECO_ZMMss_PT);
   fChain->SetBranchAddress("RECO_ZEEss_PT", &RECO_ZEEss_PT);
   fChain->SetBranchAddress("RECO_ZEM_PT", &RECO_ZEM_PT);
   fChain->SetBranchAddress("RECO_ZMMss_ETA", &RECO_ZMMss_ETA);
   fChain->SetBranchAddress("RECO_ZEEss_ETA", &RECO_ZEEss_ETA);
   fChain->SetBranchAddress("RECO_ZEM_ETA", &RECO_ZEM_ETA);
   fChain->SetBranchAddress("RECO_ZMMss_PHI", &RECO_ZMMss_PHI);
   fChain->SetBranchAddress("RECO_ZEEss_PHI", &RECO_ZEEss_PHI);
   fChain->SetBranchAddress("RECO_ZEM_PHI", &RECO_ZEM_PHI);
   fChain->SetBranchAddress("RECO_MMMM_MASS", &RECO_MMMM_MASS);
   fChain->SetBranchAddress("RECO_MMMM_PT", &RECO_MMMM_PT);
   fChain->SetBranchAddress("RECO_MMMM_ETA", &RECO_MMMM_ETA);
   fChain->SetBranchAddress("RECO_MMMM_PHI", &RECO_MMMM_PHI);
   fChain->SetBranchAddress("RECO_MMMM_MASS_REFIT", &RECO_MMMM_MASS_REFIT);
   fChain->SetBranchAddress("RECO_EEEE_MASS", &RECO_EEEE_MASS);
   fChain->SetBranchAddress("RECO_EEEE_PT", &RECO_EEEE_PT);
   fChain->SetBranchAddress("RECO_EEEE_ETA", &RECO_EEEE_ETA);
   fChain->SetBranchAddress("RECO_EEEE_PHI", &RECO_EEEE_PHI);
   fChain->SetBranchAddress("RECO_EEEE_MASS_REFIT", &RECO_EEEE_MASS_REFIT);
   fChain->SetBranchAddress("RECO_EEMM_MASS", &RECO_EEMM_MASS);
   fChain->SetBranchAddress("RECO_EEMM_PT", &RECO_EEMM_PT);
   fChain->SetBranchAddress("RECO_EEMM_ETA", &RECO_EEMM_ETA);
   fChain->SetBranchAddress("RECO_EEMM_PHI", &RECO_EEMM_PHI);
   fChain->SetBranchAddress("RECO_EEMM_MASS_REFIT", &RECO_EEMM_MASS_REFIT);
   fChain->SetBranchAddress("RECO_LLL0_MASS", &RECO_LLL0_MASS);
   fChain->SetBranchAddress("RECO_LLL1_MASS", &RECO_LLL1_MASS);
   fChain->SetBranchAddress("RECO_LLL2_MASS", &RECO_LLL2_MASS);
   fChain->SetBranchAddress("RECO_LLL3_MASS", &RECO_LLL3_MASS);
   fChain->SetBranchAddress("RECO_LLL0_PT", &RECO_LLL0_PT);
   fChain->SetBranchAddress("RECO_LLL1_PT", &RECO_LLL1_PT);
   fChain->SetBranchAddress("RECO_LLL2_PT", &RECO_LLL2_PT);
   fChain->SetBranchAddress("RECO_LLL3_PT", &RECO_LLL3_PT);
   fChain->SetBranchAddress("RECO_LLLl0_MASS", &RECO_LLLl0_MASS);
   fChain->SetBranchAddress("RECO_LLLl1_MASS", &RECO_LLLl1_MASS);
   fChain->SetBranchAddress("RECO_LLLl0_PT", &RECO_LLLl0_PT);
   fChain->SetBranchAddress("RECO_LLLl1_PT", &RECO_LLLl1_PT);
   fChain->SetBranchAddress("RECO_LLLL0ss_MASS", &RECO_LLLL0ss_MASS);
   fChain->SetBranchAddress("RECO_LLLL0ss_PT", &RECO_LLLL0ss_PT);
   fChain->SetBranchAddress("RECO_LLLL1ss_MASS", &RECO_LLLL1ss_MASS);
   fChain->SetBranchAddress("RECO_LLLL1ss_PT", &RECO_LLLL1ss_PT);
   fChain->SetBranchAddress("RECO_LLLL2ss_MASS", &RECO_LLLL2ss_MASS);
   fChain->SetBranchAddress("RECO_LLLL2ss_PT", &RECO_LLLL2ss_PT);
   fChain->SetBranchAddress("RECO_LLLL_MASS", &RECO_LLLL_MASS);
   fChain->SetBranchAddress("RECO_LLLL_PT", &RECO_LLLL_PT);
   fChain->SetBranchAddress("RECO_LLLL_ETA", &RECO_LLLL_ETA);
   fChain->SetBranchAddress("RECO_LLLL_PHI", &RECO_LLLL_PHI);
   fChain->SetBranchAddress("RECOELE_E", &RECOELE_E);
   fChain->SetBranchAddress("RECOELE_PT", &RECOELE_PT);
   fChain->SetBranchAddress("RECOELE_PTError", &RECOELE_PTError);
   fChain->SetBranchAddress("RECOELE_P", &RECOELE_P);
   fChain->SetBranchAddress("RECOELE_ETA", &RECOELE_ETA);
   fChain->SetBranchAddress("RECOELE_THETA", &RECOELE_THETA);
   fChain->SetBranchAddress("RECOELE_PHI", &RECOELE_PHI);
   fChain->SetBranchAddress("RECOELE_MASS", &RECOELE_MASS);
   fChain->SetBranchAddress("RECOELE_CHARGE", &RECOELE_CHARGE);
   fChain->SetBranchAddress("RECOELE_isEcalDriven", &RECOELE_isEcalDriven);
   fChain->SetBranchAddress("RECOELE_isTrackerDriven", &RECOELE_isTrackerDriven);
   fChain->SetBranchAddress("RECOELE_gsftrack_NPixHits", &RECOELE_gsftrack_NPixHits);
   fChain->SetBranchAddress("RECOELE_gsftrack_NStripHits", &RECOELE_gsftrack_NStripHits);
   fChain->SetBranchAddress("RECOELE_gsftrack_chi2", &RECOELE_gsftrack_chi2);
   fChain->SetBranchAddress("RECOELE_gsftrack_dxyB", &RECOELE_gsftrack_dxyB);
   fChain->SetBranchAddress("RECOELE_gsftrack_dxy", &RECOELE_gsftrack_dxy);
   fChain->SetBranchAddress("RECOELE_gsftrack_dxyError", &RECOELE_gsftrack_dxyError);
   fChain->SetBranchAddress("RECOELE_gsftrack_dzB", &RECOELE_gsftrack_dzB);
   fChain->SetBranchAddress("RECOELE_gsftrack_dz", &RECOELE_gsftrack_dz);
   fChain->SetBranchAddress("RECOELE_gsftrack_dzError", &RECOELE_gsftrack_dzError);
   fChain->SetBranchAddress("RECOELE_gsftrack_losthits", &RECOELE_gsftrack_losthits);
   fChain->SetBranchAddress("RECOELE_gsftrack_validhits", &RECOELE_gsftrack_validhits);
   fChain->SetBranchAddress("RECOELE_gsftrack_expected_inner_hits", &RECOELE_gsftrack_expected_inner_hits);
   fChain->SetBranchAddress("RECOELE_scl_E", &RECOELE_scl_E);
   fChain->SetBranchAddress("RECOELE_scl_Et", &RECOELE_scl_Et);
   fChain->SetBranchAddress("RECOELE_scl_Eta", &RECOELE_scl_Eta);
   fChain->SetBranchAddress("RECOELE_scl_Phi", &RECOELE_scl_Phi);
   fChain->SetBranchAddress("RECOELE_ep", &RECOELE_ep);
   fChain->SetBranchAddress("RECOELE_eSeedp", &RECOELE_eSeedp);
   fChain->SetBranchAddress("RECOELE_eSeedpout", &RECOELE_eSeedpout);
   fChain->SetBranchAddress("RECOELE_eElepout", &RECOELE_eElepout);
   fChain->SetBranchAddress("RECOELE_deltaEtaIn", &RECOELE_deltaEtaIn);
   fChain->SetBranchAddress("RECOELE_deltaEtaSeed", &RECOELE_deltaEtaSeed);
   fChain->SetBranchAddress("RECOELE_deltaEtaEle", &RECOELE_deltaEtaEle);
   fChain->SetBranchAddress("RECOELE_deltaPhiIn", &RECOELE_deltaPhiIn);
   fChain->SetBranchAddress("RECOELE_deltaPhiSeed", &RECOELE_deltaPhiSeed);
   fChain->SetBranchAddress("RECOELE_deltaPhiEle", &RECOELE_deltaPhiEle);
   fChain->SetBranchAddress("RECOELE_isbarrel", &RECOELE_isbarrel);
   fChain->SetBranchAddress("RECOELE_isendcap", &RECOELE_isendcap);
   fChain->SetBranchAddress("RECOELE_isGap", &RECOELE_isGap);
   fChain->SetBranchAddress("RECOELE_isEBetaGap", &RECOELE_isEBetaGap);
   fChain->SetBranchAddress("RECOELE_isEBphiGap", &RECOELE_isEBphiGap);
   fChain->SetBranchAddress("RECOELE_isEEdeeGap", &RECOELE_isEEdeeGap);
   fChain->SetBranchAddress("RECOELE_isEEringGap", &RECOELE_isEEringGap);
   fChain->SetBranchAddress("RECOELE_sigmaIetaIeta", &RECOELE_sigmaIetaIeta);
   fChain->SetBranchAddress("RECOELE_sigmaEtaEta", &RECOELE_sigmaEtaEta);
   fChain->SetBranchAddress("RECOELE_e15", &RECOELE_e15);
   fChain->SetBranchAddress("RECOELE_e25max", &RECOELE_e25max);
   fChain->SetBranchAddress("RECOELE_e55", &RECOELE_e55);
   fChain->SetBranchAddress("RECOELE_he", &RECOELE_he);
   fChain->SetBranchAddress("RECOELE_r9", &RECOELE_r9);
   fChain->SetBranchAddress("RECOELE_mva", &RECOELE_mva);
   fChain->SetBranchAddress("RECOELE_fbrem", &RECOELE_fbrem);
   fChain->SetBranchAddress("RECOELE_nbrems", &RECOELE_nbrems);
   fChain->SetBranchAddress("RECOELE_Class", &RECOELE_Class);
   fChain->SetBranchAddress("RECOELE_fbrem_mode", &RECOELE_fbrem_mode);
   fChain->SetBranchAddress("RECOELE_fbrem_mean", &RECOELE_fbrem_mean);
   fChain->SetBranchAddress("RECOELE_EGMTRACKISO", &RECOELE_EGMTRACKISO);
   fChain->SetBranchAddress("RECOELE_EGMHCALISO", &RECOELE_EGMHCALISO);
   fChain->SetBranchAddress("RECOELE_EGMECALISO", &RECOELE_EGMECALISO);
   fChain->SetBranchAddress("RECOELE_EGMX", &RECOELE_EGMX);
   fChain->SetBranchAddress("RECOELE_PFchAllPart", &RECOELE_PFchAllPart);
   fChain->SetBranchAddress("RECOELE_PFchHad", &RECOELE_PFchHad);
   fChain->SetBranchAddress("RECOELE_PFneuHad", &RECOELE_PFneuHad);
   fChain->SetBranchAddress("RECOELE_PFphoton", &RECOELE_PFphoton);
   fChain->SetBranchAddress("RECOELE_PFPUchAllPart", &RECOELE_PFPUchAllPart);
   fChain->SetBranchAddress("RECOELE_PFX_dB", &RECOELE_PFX_dB);
   fChain->SetBranchAddress("RECOELE_PFX_rho", &RECOELE_PFX_rho);
   fChain->SetBranchAddress("RECOELE_regEnergy", &RECOELE_regEnergy);
   fChain->SetBranchAddress("RECOELE_regEnergyError", &RECOELE_regEnergyError);
   fChain->SetBranchAddress("RECOELE_SIP", &RECOELE_SIP);
   fChain->SetBranchAddress("RECOELE_IP", &RECOELE_IP);
   fChain->SetBranchAddress("RECOELE_IPERROR", &RECOELE_IPERROR);
   fChain->SetBranchAddress("RECOELE_SIP_KF", &RECOELE_SIP_KF);
   fChain->SetBranchAddress("RECOELE_IP_KF", &RECOELE_IP_KF);
   fChain->SetBranchAddress("RECOELE_IPERROR_KF", &RECOELE_IPERROR_KF);
   fChain->SetBranchAddress("RECOELE_SIP_GD", &RECOELE_SIP_GD);
   fChain->SetBranchAddress("RECOELE_SIP_GDEEEE", &RECOELE_SIP_GDEEEE);
   fChain->SetBranchAddress("RECOELE_SIP_Std", &RECOELE_SIP_Std);
   fChain->SetBranchAddress("RECOELE_SIP_StdEEEE", &RECOELE_SIP_StdEEEE);
   fChain->SetBranchAddress("RECOELE_SIP_Kin", &RECOELE_SIP_Kin);
   fChain->SetBranchAddress("RECOELE_SIP_KinEEEE", &RECOELE_SIP_KinEEEE);
   fChain->SetBranchAddress("RECOELE_STIP", &RECOELE_STIP);
   fChain->SetBranchAddress("RECOELE_SLIP", &RECOELE_SLIP);
   fChain->SetBranchAddress("RECOELE_TIP", &RECOELE_TIP);
   fChain->SetBranchAddress("RECOELE_LIP", &RECOELE_LIP);
   fChain->SetBranchAddress("RECOELE_TIPERROR", &RECOELE_TIPERROR);
   fChain->SetBranchAddress("RECOELE_LIPERROR", &RECOELE_LIPERROR);
   fChain->SetBranchAddress("RECOELE_sclRawE", &RECOELE_sclRawE);
   fChain->SetBranchAddress("RECOELE_sclX", &RECOELE_sclX);
   fChain->SetBranchAddress("RECOELE_sclY", &RECOELE_sclY);
   fChain->SetBranchAddress("RECOELE_sclZ", &RECOELE_sclZ);
   fChain->SetBranchAddress("RECOELE_seedSubdet1", &RECOELE_seedSubdet1);
   fChain->SetBranchAddress("RECOELE_seedDphi1", &RECOELE_seedDphi1);
   fChain->SetBranchAddress("RECOELE_seedDrz1", &RECOELE_seedDrz1);
   fChain->SetBranchAddress("RECOELE_seedSubdet2", &RECOELE_seedSubdet2);
   fChain->SetBranchAddress("RECOELE_seedDphi2", &RECOELE_seedDphi2);
   fChain->SetBranchAddress("RECOELE_seedDrz2", &RECOELE_seedDrz2);
   fChain->SetBranchAddress("RECOELE_eidVeryLoose", &RECOELE_eidVeryLoose);
   fChain->SetBranchAddress("RECOELE_eidLoose", &RECOELE_eidLoose);
   fChain->SetBranchAddress("RECOELE_eidMedium", &RECOELE_eidMedium);
   fChain->SetBranchAddress("RECOELE_eidTight", &RECOELE_eidTight);
   fChain->SetBranchAddress("RECOELE_eidHZZVeryLoose", &RECOELE_eidHZZVeryLoose);
   fChain->SetBranchAddress("RECOELE_eidHZZLoose", &RECOELE_eidHZZLoose);
   fChain->SetBranchAddress("RECOELE_eidHZZMedium", &RECOELE_eidHZZMedium);
   fChain->SetBranchAddress("RECOELE_eidHZZTight", &RECOELE_eidHZZTight);
   fChain->SetBranchAddress("RECOELE_mvaTrigV0", &RECOELE_mvaTrigV0);
   fChain->SetBranchAddress("RECOELE_mvaNonTrigV0", &RECOELE_mvaNonTrigV0);
   fChain->SetBranchAddress("RECOELE_COV", &RECOELE_COV);
   fChain->SetBranchAddress("RECOMU_isPFMu", &RECOMU_isPFMu);
   fChain->SetBranchAddress("RECOMU_isGlobalMu", &RECOMU_isGlobalMu);
   fChain->SetBranchAddress("RECOMU_isStandAloneMu", &RECOMU_isStandAloneMu);
   fChain->SetBranchAddress("RECOMU_isTrackerMu", &RECOMU_isTrackerMu);
   fChain->SetBranchAddress("RECOMU_isCaloMu", &RECOMU_isCaloMu);
   fChain->SetBranchAddress("RECOMU_isTrackerHighPtMu", &RECOMU_isTrackerHighPtMu);
   fChain->SetBranchAddress("RECOMU_E", &RECOMU_E);
   fChain->SetBranchAddress("RECOMU_PT", &RECOMU_PT);
   fChain->SetBranchAddress("RECOMU_mubesttrkPTError", &RECOMU_mubesttrkPTError);
   fChain->SetBranchAddress("RECOMU_P", &RECOMU_P);
   fChain->SetBranchAddress("RECOMU_ETA", &RECOMU_ETA);
   fChain->SetBranchAddress("RECOMU_THETA", &RECOMU_THETA);
   fChain->SetBranchAddress("RECOMU_PHI", &RECOMU_PHI);
   fChain->SetBranchAddress("RECOMU_MASS", &RECOMU_MASS);
   fChain->SetBranchAddress("RECOMU_CHARGE", &RECOMU_CHARGE);
   fChain->SetBranchAddress("RECOMU_COV", &RECOMU_COV);
   fChain->SetBranchAddress("RECOMU_TRACKISO", &RECOMU_TRACKISO);
   fChain->SetBranchAddress("RECOMU_TRACKISO_SUMPT", &RECOMU_TRACKISO_SUMPT);
   fChain->SetBranchAddress("RECOMU_HCALISO", &RECOMU_HCALISO);
   fChain->SetBranchAddress("RECOMU_ECALISO", &RECOMU_ECALISO);
   fChain->SetBranchAddress("RECOMU_X", &RECOMU_X);
   fChain->SetBranchAddress("RECOMU_PFchHad", &RECOMU_PFchHad);
   fChain->SetBranchAddress("RECOMU_PFneuHad", &RECOMU_PFneuHad);
   fChain->SetBranchAddress("RECOMU_PFphoton", &RECOMU_PFphoton);
   fChain->SetBranchAddress("RECOMU_PFPUchAllPart", &RECOMU_PFPUchAllPart);
   fChain->SetBranchAddress("RECOMU_PFX_dB", &RECOMU_PFX_dB);
   fChain->SetBranchAddress("RECOMU_PFX_rho", &RECOMU_PFX_rho);
   fChain->SetBranchAddress("RECOPFPHOT_PFchHad", &RECOPFPHOT_PFchHad);
   fChain->SetBranchAddress("RECOPFPHOT_PFneuHad", &RECOPFPHOT_PFneuHad);
   fChain->SetBranchAddress("RECOPFPHOT_PFphoton", &RECOPFPHOT_PFphoton);
   fChain->SetBranchAddress("RECOPFPHOT_PFPUchAllPart", &RECOPFPHOT_PFPUchAllPart);
   fChain->SetBranchAddress("RECOPFPHOT_PFX_rho", &RECOPFPHOT_PFX_rho);
   fChain->SetBranchAddress("RECOMU_SIP", &RECOMU_SIP);
   fChain->SetBranchAddress("RECOMU_IP", &RECOMU_IP);
   fChain->SetBranchAddress("RECOMU_IPERROR", &RECOMU_IPERROR);
   fChain->SetBranchAddress("RECOMU_SIP_KF", &RECOMU_SIP_KF);
   fChain->SetBranchAddress("RECOMU_IP_KF", &RECOMU_IP_KF);
   fChain->SetBranchAddress("RECOMU_IPERROR_KF", &RECOMU_IPERROR_KF);
   fChain->SetBranchAddress("RECOMU_SIP_GD", &RECOMU_SIP_GD);
   fChain->SetBranchAddress("RECOMU_SIP_GDMMMM", &RECOMU_SIP_GDMMMM);
   fChain->SetBranchAddress("RECOMU_SIP_Std", &RECOMU_SIP_Std);
   fChain->SetBranchAddress("RECOMU_SIP_StdMMMM", &RECOMU_SIP_StdMMMM);
   fChain->SetBranchAddress("RECOMU_SIP_Kin", &RECOMU_SIP_Kin);
   fChain->SetBranchAddress("RECOMU_SIP_KinMMMM", &RECOMU_SIP_KinMMMM);
   fChain->SetBranchAddress("RECOMU_STIP", &RECOMU_STIP);
   fChain->SetBranchAddress("RECOMU_SLIP", &RECOMU_SLIP);
   fChain->SetBranchAddress("RECOMU_TIP", &RECOMU_TIP);
   fChain->SetBranchAddress("RECOMU_LIP", &RECOMU_LIP);
   fChain->SetBranchAddress("RECOMU_TIPERROR", &RECOMU_TIPERROR);
   fChain->SetBranchAddress("RECOMU_LIPERROR", &RECOMU_LIPERROR);
   fChain->SetBranchAddress("RECOMU_caloCompatibility", &RECOMU_caloCompatibility);
   fChain->SetBranchAddress("RECOMU_segmentCompatibility", &RECOMU_segmentCompatibility);
   fChain->SetBranchAddress("RECOMU_numberOfMatches", &RECOMU_numberOfMatches);
   fChain->SetBranchAddress("RECOMU_numberOfMatchedStations", &RECOMU_numberOfMatchedStations);
   fChain->SetBranchAddress("RECOMU_glbmuPromptTight", &RECOMU_glbmuPromptTight);
   fChain->SetBranchAddress("RECOMU_trkmuArbitration", &RECOMU_trkmuArbitration);
   fChain->SetBranchAddress("RECOMU_trkmu2DCompatibilityLoose", &RECOMU_trkmu2DCompatibilityLoose);
   fChain->SetBranchAddress("RECOMU_trkmu2DCompatibilityTight", &RECOMU_trkmu2DCompatibilityTight);
   fChain->SetBranchAddress("RECOMU_trkmuOneStationLoose", &RECOMU_trkmuOneStationLoose);
   fChain->SetBranchAddress("RECOMU_trkmuOneStationTight", &RECOMU_trkmuOneStationTight);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationLoose", &RECOMU_trkmuLastStationLoose);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationTight", &RECOMU_trkmuLastStationTight);
   fChain->SetBranchAddress("RECOMU_trkmuOneStationAngLoose", &RECOMU_trkmuOneStationAngLoose);
   fChain->SetBranchAddress("RECOMU_trkmuOneStationAngTight", &RECOMU_trkmuOneStationAngTight);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationAngLoose", &RECOMU_trkmuLastStationAngLoose);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationAngTight", &RECOMU_trkmuLastStationAngTight);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationOptimizedLowPtLoose", &RECOMU_trkmuLastStationOptimizedLowPtLoose);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationOptimizedLowPtTight", &RECOMU_trkmuLastStationOptimizedLowPtTight);
   fChain->SetBranchAddress("RECOMU_mutrkPT", &RECOMU_mutrkPT);
   fChain->SetBranchAddress("RECOMU_mutrkPTError", &RECOMU_mutrkPTError);
   fChain->SetBranchAddress("RECOMU_mutrkDxy", &RECOMU_mutrkDxy);
   fChain->SetBranchAddress("RECOMU_mutrkDxyError", &RECOMU_mutrkDxyError);
   fChain->SetBranchAddress("RECOMU_mutrkDxyB", &RECOMU_mutrkDxyB);
   fChain->SetBranchAddress("RECOMU_mutrkDz", &RECOMU_mutrkDz);
   fChain->SetBranchAddress("RECOMU_mutrkDzError", &RECOMU_mutrkDzError);
   fChain->SetBranchAddress("RECOMU_mutrkDzB", &RECOMU_mutrkDzB);
   fChain->SetBranchAddress("RECOMU_mutrkChi2PerNdof", &RECOMU_mutrkChi2PerNdof);
   fChain->SetBranchAddress("RECOMU_mutrkCharge", &RECOMU_mutrkCharge);
   fChain->SetBranchAddress("RECOMU_mutrkNHits", &RECOMU_mutrkNHits);
   fChain->SetBranchAddress("RECOMU_mutrkNStripHits", &RECOMU_mutrkNStripHits);
   fChain->SetBranchAddress("RECOMU_mutrkNPixHits", &RECOMU_mutrkNPixHits);
   fChain->SetBranchAddress("RECOMU_mutrkNMuonHits", &RECOMU_mutrkNMuonHits);
   fChain->SetBranchAddress("RECOMU_mutrktrackerLayersWithMeasurement", &RECOMU_mutrktrackerLayersWithMeasurement);
   fChain->SetBranchAddress("RECOMU_muInnertrkDxy", &RECOMU_muInnertrkDxy);
   fChain->SetBranchAddress("RECOMU_muInnertrkDxyError", &RECOMU_muInnertrkDxyError);
   fChain->SetBranchAddress("RECOMU_muInnertrkDxyB", &RECOMU_muInnertrkDxyB);
   fChain->SetBranchAddress("RECOMU_muInnertrkDz", &RECOMU_muInnertrkDz);
   fChain->SetBranchAddress("RECOMU_muInnertrkDzError", &RECOMU_muInnertrkDzError);
   fChain->SetBranchAddress("RECOMU_muInnertrkDzB", &RECOMU_muInnertrkDzB);
   fChain->SetBranchAddress("RECOMU_muInnertrkChi2PerNdof", &RECOMU_muInnertrkChi2PerNdof);
   fChain->SetBranchAddress("RECOMU_muInnertrktrackerLayersWithMeasurement", &RECOMU_muInnertrktrackerLayersWithMeasurement);
   fChain->SetBranchAddress("RECOMU_muInnertrkPT", &RECOMU_muInnertrkPT);
   fChain->SetBranchAddress("RECOMU_muInnertrkPTError", &RECOMU_muInnertrkPTError);
   fChain->SetBranchAddress("RECOMU_muInnertrkCharge", &RECOMU_muInnertrkCharge);
   fChain->SetBranchAddress("RECOMU_muInnertrkNHits", &RECOMU_muInnertrkNHits);
   fChain->SetBranchAddress("RECOMU_muInnertrkNStripHits", &RECOMU_muInnertrkNStripHits);
   fChain->SetBranchAddress("RECOMU_muInnertrkNPixHits", &RECOMU_muInnertrkNPixHits);
   fChain->SetBranchAddress("RECOMU_mubesttrkType", &RECOMU_mubesttrkType);
   fChain->SetBranchAddress("RECOMU_mubesttrkDxy", &RECOMU_mubesttrkDxy);
   fChain->SetBranchAddress("RECOMU_mubesttrkDxyError", &RECOMU_mubesttrkDxyError);
   fChain->SetBranchAddress("RECOMU_mubesttrkDz", &RECOMU_mubesttrkDz);
   fChain->SetBranchAddress("RECOMU_mubesttrkDzError", &RECOMU_mubesttrkDzError);
   fChain->SetBranchAddress("ftsigma", &ftsigma);
   fChain->SetBranchAddress("gdX", &gdX);
   fChain->SetBranchAddress("gdY", &gdY);
   fChain->SetBranchAddress("gdZ", &gdZ);
   fChain->SetBranchAddress("ftsigmalag", &ftsigmalag);
   fChain->SetBranchAddress("gdlagX", &gdlagX);
   fChain->SetBranchAddress("gdlagY", &gdlagY);
   fChain->SetBranchAddress("gdlagZ", &gdlagZ);
   fChain->SetBranchAddress("gdlagProb", &gdlagProb);
   fChain->SetBranchAddress("gdlagNdof", &gdlagNdof);
   fChain->SetBranchAddress("ftsigmaMMMM", &ftsigmaMMMM);
   fChain->SetBranchAddress("gdXMMMM", &gdXMMMM);
   fChain->SetBranchAddress("gdYMMMM", &gdYMMMM);
   fChain->SetBranchAddress("gdZMMMM", &gdZMMMM);
   fChain->SetBranchAddress("ftsigmalagMMMM", &ftsigmalagMMMM);
   fChain->SetBranchAddress("gdlagXMMMM", &gdlagXMMMM);
   fChain->SetBranchAddress("gdlagYMMMM", &gdlagYMMMM);
   fChain->SetBranchAddress("gdlagZMMMM", &gdlagZMMMM);
   fChain->SetBranchAddress("gdlagProbMMMM", &gdlagProbMMMM);
   fChain->SetBranchAddress("gdlagNdofMMMM", &gdlagNdofMMMM);
   fChain->SetBranchAddress("ftsigmaEEEE", &ftsigmaEEEE);
   fChain->SetBranchAddress("gdXEEEE", &gdXEEEE);
   fChain->SetBranchAddress("gdYEEEE", &gdYEEEE);
   fChain->SetBranchAddress("gdZEEEE", &gdZEEEE);
   fChain->SetBranchAddress("ftsigmalagEEEE", &ftsigmalagEEEE);
   fChain->SetBranchAddress("gdlagXEEEE", &gdlagXEEEE);
   fChain->SetBranchAddress("gdlagYEEEE", &gdlagYEEEE);
   fChain->SetBranchAddress("gdlagZEEEE", &gdlagZEEEE);
   fChain->SetBranchAddress("gdlagProbEEEE", &gdlagProbEEEE);
   fChain->SetBranchAddress("gdlagNdofEEEE", &gdlagNdofEEEE);
   fChain->SetBranchAddress("StdFitVertexX", &StdFitVertexX);
   fChain->SetBranchAddress("StdFitVertexY", &StdFitVertexY);
   fChain->SetBranchAddress("StdFitVertexZ", &StdFitVertexZ);
   fChain->SetBranchAddress("StdFitVertexChi2r", &StdFitVertexChi2r);
   fChain->SetBranchAddress("StdFitVertexProb", &StdFitVertexProb);
   fChain->SetBranchAddress("StdFitVertexTrack_PT", &StdFitVertexTrack_PT);
   fChain->SetBranchAddress("StdFitVertexTrack_ETA", &StdFitVertexTrack_ETA);
   fChain->SetBranchAddress("StdFitVertexTrack_PHI", &StdFitVertexTrack_PHI);
   fChain->SetBranchAddress("KinFitVertexX", &KinFitVertexX);
   fChain->SetBranchAddress("KinFitVertexY", &KinFitVertexY);
   fChain->SetBranchAddress("KinFitVertexZ", &KinFitVertexZ);
   fChain->SetBranchAddress("KinFitVertexChi2r", &KinFitVertexChi2r);
   fChain->SetBranchAddress("KinFitVertexProb", &KinFitVertexProb);
   fChain->SetBranchAddress("StdFitVertexXMMMM", &StdFitVertexXMMMM);
   fChain->SetBranchAddress("StdFitVertexYMMMM", &StdFitVertexYMMMM);
   fChain->SetBranchAddress("StdFitVertexZMMMM", &StdFitVertexZMMMM);
   fChain->SetBranchAddress("StdFitVertexChi2rMMMM", &StdFitVertexChi2rMMMM);
   fChain->SetBranchAddress("StdFitVertexProbMMMM", &StdFitVertexProbMMMM);
   fChain->SetBranchAddress("StdFitVertexTrackMMMM_PT", &StdFitVertexTrackMMMM_PT);
   fChain->SetBranchAddress("StdFitVertexTrackMMMM_ETA", &StdFitVertexTrackMMMM_ETA);
   fChain->SetBranchAddress("StdFitVertexTrackMMMM_PHI", &StdFitVertexTrackMMMM_PHI);
   fChain->SetBranchAddress("KinFitVertexXMMMM", &KinFitVertexXMMMM);
   fChain->SetBranchAddress("KinFitVertexYMMMM", &KinFitVertexYMMMM);
   fChain->SetBranchAddress("KinFitVertexZMMMM", &KinFitVertexZMMMM);
   fChain->SetBranchAddress("KinFitVertexChi2rMMMM", &KinFitVertexChi2rMMMM);
   fChain->SetBranchAddress("KinFitVertexProbMMMM", &KinFitVertexProbMMMM);
   fChain->SetBranchAddress("StdFitVertexXEEEE", &StdFitVertexXEEEE);
   fChain->SetBranchAddress("StdFitVertexYEEEE", &StdFitVertexYEEEE);
   fChain->SetBranchAddress("StdFitVertexZEEEE", &StdFitVertexZEEEE);
   fChain->SetBranchAddress("StdFitVertexChi2rEEEE", &StdFitVertexChi2rEEEE);
   fChain->SetBranchAddress("StdFitVertexProbEEEE", &StdFitVertexProbEEEE);
   fChain->SetBranchAddress("StdFitVertexTrackEEEE_PT", &StdFitVertexTrackEEEE_PT);
   fChain->SetBranchAddress("StdFitVertexTrackEEEE_ETA", &StdFitVertexTrackEEEE_ETA);
   fChain->SetBranchAddress("StdFitVertexTrackEEEE_PHI", &StdFitVertexTrackEEEE_PHI);
   fChain->SetBranchAddress("KinFitVertexXEEEE", &KinFitVertexXEEEE);
   fChain->SetBranchAddress("KinFitVertexYEEEE", &KinFitVertexYEEEE);
   fChain->SetBranchAddress("KinFitVertexZEEEE", &KinFitVertexZEEEE);
   fChain->SetBranchAddress("KinFitVertexChi2rEEEE", &KinFitVertexChi2rEEEE);
   fChain->SetBranchAddress("KinFitVertexProbEEEE", &KinFitVertexProbEEEE);
   fChain->SetBranchAddress("StdFitVertexChi2rMMM", &StdFitVertexChi2rMMM);
   fChain->SetBranchAddress("StdFitVertexProbMMM", &StdFitVertexProbMMM);
   fChain->SetBranchAddress("StdFitVertexChi2rMME", &StdFitVertexChi2rMME);
   fChain->SetBranchAddress("StdFitVertexProbMME", &StdFitVertexProbMME);
   fChain->SetBranchAddress("StdFitVertexChi2rEEE", &StdFitVertexChi2rEEE);
   fChain->SetBranchAddress("StdFitVertexProbEEE", &StdFitVertexProbEEE);
   fChain->SetBranchAddress("StdFitVertexChi2rMEE", &StdFitVertexChi2rMEE);
   fChain->SetBranchAddress("StdFitVertexProbMEE", &StdFitVertexProbMEE);
   fChain->SetBranchAddress("StdFitVertexChi2rDiLep", &StdFitVertexChi2rDiLep);
   fChain->SetBranchAddress("StdFitVertexProbDiLep", &StdFitVertexProbDiLep);
   fChain->SetBranchAddress("ConvMapDist", &ConvMapDist);
   fChain->SetBranchAddress("ConvMapDcot", &ConvMapDcot);
   fChain->SetBranchAddress("RECOMU_MatchingMCTruth", &RECOMU_MatchingMCTruth);
   fChain->SetBranchAddress("RECOMU_MatchingMCpT", &RECOMU_MatchingMCpT);
   fChain->SetBranchAddress("RECOMU_MatchingMCEta", &RECOMU_MatchingMCEta);
   fChain->SetBranchAddress("RECOMU_MatchingMCPhi", &RECOMU_MatchingMCPhi);
   fChain->SetBranchAddress("RECOELE_MatchingMCTruth", &RECOELE_MatchingMCTruth);
   fChain->SetBranchAddress("RECOELE_MatchingMCpT", &RECOELE_MatchingMCpT);
   fChain->SetBranchAddress("RECOELE_MatchingMCEta", &RECOELE_MatchingMCEta);
   fChain->SetBranchAddress("RECOELE_MatchingMCPhi", &RECOELE_MatchingMCPhi);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCTruth", &RECOPHOT_MatchingMCTruth);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCpT", &RECOPHOT_MatchingMCpT);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCEta", &RECOPHOT_MatchingMCEta);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCPhi", &RECOPHOT_MatchingMCPhi);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCTruth", &RECOzMuMu_MatchingMCTruth);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCpT", &RECOzMuMu_MatchingMCpT);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCmass", &RECOzMuMu_MatchingMCmass);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCEta", &RECOzMuMu_MatchingMCEta);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCPhi", &RECOzMuMu_MatchingMCPhi);
   fChain->SetBranchAddress("RECOzEE_MatchingMCTruth", &RECOzEE_MatchingMCTruth);
   fChain->SetBranchAddress("RECOzEE_MatchingMCpT", &RECOzEE_MatchingMCpT);
   fChain->SetBranchAddress("RECOzEE_MatchingMCmass", &RECOzEE_MatchingMCmass);
   fChain->SetBranchAddress("RECOzEE_MatchingMCEta", &RECOzEE_MatchingMCEta);
   fChain->SetBranchAddress("RECOzEE_MatchingMCPhi", &RECOzEE_MatchingMCPhi);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCTruth", &RECOHzzEEEE_MatchingMCTruth);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCpT", &RECOHzzEEEE_MatchingMCpT);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCmass", &RECOHzzEEEE_MatchingMCmass);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCEta", &RECOHzzEEEE_MatchingMCEta);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCPhi", &RECOHzzEEEE_MatchingMCPhi);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCTruth", &RECOHzzEEMM_MatchingMCTruth);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCpT", &RECOHzzEEMM_MatchingMCpT);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCmass", &RECOHzzEEMM_MatchingMCmass);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCEta", &RECOHzzEEMM_MatchingMCEta);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCPhi", &RECOHzzEEMM_MatchingMCPhi);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCTruth", &RECOHzzMMMM_MatchingMCTruth);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCpT", &RECOHzzMMMM_MatchingMCpT);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCmass", &RECOHzzMMMM_MatchingMCmass);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCEta", &RECOHzzMMMM_MatchingMCEta);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCPhi", &RECOHzzMMMM_MatchingMCPhi);
   fChain->SetBranchAddress("RECO_NMU", &RECO_NMU);
   fChain->SetBranchAddress("RECO_NELE", &RECO_NELE);
   fChain->SetBranchAddress("RECO_NTRACK", &RECO_NTRACK);
   fChain->SetBranchAddress("RECO_TRACK_PT", &RECO_TRACK_PT);
   fChain->SetBranchAddress("RECO_TRACK_ETA", &RECO_TRACK_ETA);
   fChain->SetBranchAddress("RECO_TRACK_PHI", &RECO_TRACK_PHI);
   fChain->SetBranchAddress("RECO_TRACK_CHI2", &RECO_TRACK_CHI2);
   fChain->SetBranchAddress("RECO_TRACK_CHI2RED", &RECO_TRACK_CHI2RED);
   fChain->SetBranchAddress("RECO_TRACK_CHI2PROB", &RECO_TRACK_CHI2PROB);
   fChain->SetBranchAddress("RECO_TRACK_NHITS", &RECO_TRACK_NHITS);
   fChain->SetBranchAddress("RECO_TRACK_DXY", &RECO_TRACK_DXY);
   fChain->SetBranchAddress("RECO_TRACK_DXYERR", &RECO_TRACK_DXYERR);
   fChain->SetBranchAddress("RECO_TRACK_DZ", &RECO_TRACK_DZ);
   fChain->SetBranchAddress("RECO_TRACK_DZERR", &RECO_TRACK_DZERR);
   fChain->SetBranchAddress("RECO_NPHOT", &RECO_NPHOT);
   fChain->SetBranchAddress("RECOPHOT_PT", &RECOPHOT_PT);
   fChain->SetBranchAddress("RECOPHOT_ETA", &RECOPHOT_ETA);
   fChain->SetBranchAddress("RECOPHOT_PHI", &RECOPHOT_PHI);
   fChain->SetBranchAddress("RECOPHOT_THETA", &RECOPHOT_THETA);
   fChain->SetBranchAddress("RECO_NPFPHOT", &RECO_NPFPHOT);
   fChain->SetBranchAddress("RECOPFPHOT_PT", &RECOPFPHOT_PT);
   fChain->SetBranchAddress("RECOPFPHOT_PTError", &RECOPFPHOT_PTError);
   fChain->SetBranchAddress("RECOPFPHOT_ETA", &RECOPFPHOT_ETA);
   fChain->SetBranchAddress("RECOPFPHOT_PHI", &RECOPFPHOT_PHI);
   fChain->SetBranchAddress("RECOPFPHOT_THETA", &RECOPFPHOT_THETA);
   fChain->SetBranchAddress("BeamSpot_X", &BeamSpot_X);
   fChain->SetBranchAddress("BeamSpot_Y", &BeamSpot_Y);
   fChain->SetBranchAddress("BeamSpot_Z", &BeamSpot_Z);
   fChain->SetBranchAddress("RECO_NVTX", &RECO_NVTX);
   fChain->SetBranchAddress("RECO_VERTEX_x", &RECO_VERTEX_x);
   fChain->SetBranchAddress("RECO_VERTEX_y", &RECO_VERTEX_y);
   fChain->SetBranchAddress("RECO_VERTEX_z", &RECO_VERTEX_z);
   fChain->SetBranchAddress("RECO_VERTEX_ndof", &RECO_VERTEX_ndof);
   fChain->SetBranchAddress("RECO_VERTEX_chi2", &RECO_VERTEX_chi2);
   fChain->SetBranchAddress("RECO_VERTEX_ntracks", &RECO_VERTEX_ntracks);
   fChain->SetBranchAddress("RECO_VERTEXPROB", &RECO_VERTEXPROB);
   fChain->SetBranchAddress("RECO_VERTEX_isValid", &RECO_VERTEX_isValid);
   fChain->SetBranchAddress("RECO_VERTEX_TRACK_PT", &RECO_VERTEX_TRACK_PT);
   fChain->SetBranchAddress("RECO_PFJET_N", &RECO_PFJET_N);
   fChain->SetBranchAddress("RECO_PFJET_CHARGE", &RECO_PFJET_CHARGE);
   fChain->SetBranchAddress("RECO_PFJET_ET", &RECO_PFJET_ET);
   fChain->SetBranchAddress("RECO_PFJET_PT", &RECO_PFJET_PT);
   fChain->SetBranchAddress("RECO_PFJET_ETA", &RECO_PFJET_ETA);
   fChain->SetBranchAddress("RECO_PFJET_PHI", &RECO_PFJET_PHI);
   fChain->SetBranchAddress("RECO_PFJET_PUID", &RECO_PFJET_PUID);
   fChain->SetBranchAddress("RECO_PFJET_PUID_MVA", &RECO_PFJET_PUID_MVA);
   fChain->SetBranchAddress("RHO_ele", &RHO_ele);
   fChain->SetBranchAddress("RHO_mu", &RHO_mu);
   fChain->SetBranchAddress("RECO_CALOMET", &RECO_CALOMET);
   fChain->SetBranchAddress("RECO_PFMET", &RECO_PFMET);
   fChain->SetBranchAddress("RECO_PFMET_X", &RECO_PFMET_X);
   fChain->SetBranchAddress("RECO_PFMET_Y", &RECO_PFMET_Y);
   fChain->SetBranchAddress("RECO_PFMET_PHI", &RECO_PFMET_PHI);
   fChain->SetBranchAddress("RECO_PFMET_THETA", &RECO_PFMET_THETA);
   //fChain->SetBranchAddress("RECO_PFMETT1", &RECO_PFMETT1);
   //fChain->SetBranchAddress("RECO_PFMETT1_X", &RECO_PFMETT1_X);
   //fChain->SetBranchAddress("RECO_PFMETT1_Y", &RECO_PFMETT1_Y);
   //fChain->SetBranchAddress("RECO_PFMETT1_PHI", &RECO_PFMETT1_PHI);
   //fChain->SetBranchAddress("RECO_PFMETT1_THETA", &RECO_PFMETT1_THETA);
   fChain->SetBranchAddress("RECO_TCMET", &RECO_TCMET);
   fChain->SetBranchAddress("RECO_CORMETMUONS", &RECO_CORMETMUONS);
   fChain->SetBranchAddress("tCHighEff_BTagJet_PT", &tCHighEff_BTagJet_PT);
   fChain->SetBranchAddress("tCHighPur_BTagJet_PT", &tCHighPur_BTagJet_PT);
   fChain->SetBranchAddress("cSV_BTagJet_PT", &cSV_BTagJet_PT);
   fChain->SetBranchAddress("tCHighEff_BTagJet_ETA", &tCHighEff_BTagJet_ETA);
   fChain->SetBranchAddress("tCHighPur_BTagJet_ETA", &tCHighPur_BTagJet_ETA);
   fChain->SetBranchAddress("cSV_BTagJet_ETA", &cSV_BTagJet_ETA);
   fChain->SetBranchAddress("tCHighEff_BTagJet_PHI", &tCHighEff_BTagJet_PHI);
   fChain->SetBranchAddress("tCHighPur_BTagJet_PHI", &tCHighPur_BTagJet_PHI);
   fChain->SetBranchAddress("cSV_BTagJet_PHI", &cSV_BTagJet_PHI);
   fChain->SetBranchAddress("tCHighEff_BTagJet_DISCR", &tCHighEff_BTagJet_DISCR);
   fChain->SetBranchAddress("tCHighPur_BTagJet_DISCR", &tCHighPur_BTagJet_DISCR);
   fChain->SetBranchAddress("cSV_BTagJet_DISCR", &cSV_BTagJet_DISCR);
   fChain->SetBranchAddress("cSV_BTagJet_ET", &cSV_BTagJet_ET);
   fChain->SetBranchAddress("RECO_PFJET_PT_UncDn", &RECO_PFJET_PT_UncDn);
   fChain->SetBranchAddress("RECO_PFMET_JetEnUp", &RECO_PFMET_JetEnUp);
   fChain->SetBranchAddress("RECO_PFMET_JetEnDn", &RECO_PFMET_JetEnDn);
   fChain->SetBranchAddress("RECO_PFMET_ElectronEnUp", &RECO_PFMET_ElectronEnUp);
   fChain->SetBranchAddress("RECO_PFMET_ElectronEnDn", &RECO_PFMET_ElectronEnDn);
   fChain->SetBranchAddress("RECO_PFMET_MuonEnUp", &RECO_PFMET_MuonEnUp);
   fChain->SetBranchAddress("RECO_PFMET_MuonEnDn", &RECO_PFMET_MuonEnDn);
   fChain->SetBranchAddress("RECO_PFMET_JetResUp", &RECO_PFMET_JetResUp);
   fChain->SetBranchAddress("RECO_PFMET_JetResDn", &RECO_PFMET_JetResDn);
   fChain->SetBranchAddress("RECO_PFMET_UnclusteredEnUp", &RECO_PFMET_UnclusteredEnUp);
   fChain->SetBranchAddress("RECO_PFMET_UnclusteredEnDn", &RECO_PFMET_UnclusteredEnDn);
   fChain->SetBranchAddress("RECO_PFMET_PhotonEnUp", &RECO_PFMET_PhotonEnUp);
   fChain->SetBranchAddress("RECO_PFMET_PhotonEnDn", &RECO_PFMET_PhotonEnDn);
  
  
  if (fChain == 0) return;  
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  
  
  //Defines a new TTree to store samples for Z+X
  //It might be needed to train NNs with such events
  TString ofile_name = store_path + "output_ZplusX_" + datasetName + ".root";
  TFile *theFile = new TFile(ofile_name,"RECREATE");

  Float_t f_weight=-999;
  Int_t f_run=-999, f_lumi=-999, f_event=-999;
  Float_t f_zee_mu_pt_loose=-999, f_zee_e_pt_loose=-999, f_zmumu_mu_pt_loose=-999, f_zmumu_e_pt_loose=-999;
  Float_t f_zee_mu_eta_loose=-999, f_zee_e_eta_loose=-999, f_zmumu_mu_eta_loose=-999, f_zmumu_e_eta_loose=-999;
  Float_t f_zee_mu_charge_loose=-999, f_zee_e_charge_loose=-999, f_zmumu_mu_charge_loose=-999, f_zmumu_e_charge_loose=-999;
  Float_t f_zee_mu_pt_tight=-999, f_zee_e_pt_tight=-999, f_zmumu_mu_pt_tight=-999, f_zmumu_e_pt_tight=-999;
  Float_t f_zee_mu_eta_tight=-999, f_zee_e_eta_tight=-999, f_zmumu_mu_eta_tight=-999, f_zmumu_e_eta_tight=-999;
  Float_t f_zee_mu_charge_tight=-999, f_zee_e_charge_tight=-999, f_zmumu_mu_charge_tight=-999, f_zmumu_e_charge_tight=-999;
  TTree *treeFR = new TTree("FakeRate","Fake rate");
  treeFR->Branch("f_run", &f_run, "f_run/I");
  treeFR->Branch("f_lumi", &f_lumi, "f_lumi/I");
  treeFR->Branch("f_event", &f_event, "f_event/I");
  treeFR->Branch("f_weight", &f_weight, "f_weight/F");
  treeFR->Branch("f_zee_mu_pt_loose", &f_zee_mu_pt_loose, "f_zee_mu_pt_loose/F");
  treeFR->Branch("f_zee_mu_eta_loose", &f_zee_mu_eta_loose, "f_zee_mu_eta_loose/F");
  treeFR->Branch("f_zee_e_pt_loose", &f_zee_e_pt_loose, "f_zee_e_pt_loose/F");
  treeFR->Branch("f_zee_e_eta_loose", &f_zee_e_eta_loose, "f_zee_e_eta_loose/F");
  treeFR->Branch("f_zmumu_mu_pt_loose", &f_zmumu_mu_pt_loose, "f_zmumu_mu_pt_loose/F");
  treeFR->Branch("f_zmumu_mu_eta_loose", &f_zmumu_mu_eta_loose, "f_zmumu_mu_eta_loose/F");
  treeFR->Branch("f_zmumu_e_pt_loose", &f_zmumu_e_pt_loose, "f_zmumu_e_pt_loose/F");
  treeFR->Branch("f_zmumu_e_eta_loose", &f_zmumu_e_eta_loose, "f_zmumu_e_eta_loose/F");
  treeFR->Branch("f_zee_mu_pt_tight", &f_zee_mu_pt_tight, "f_zee_mu_pt_tight/F");
  treeFR->Branch("f_zee_mu_eta_tight", &f_zee_mu_eta_tight, "f_zee_mu_eta_tight/F");
  treeFR->Branch("f_zee_e_pt_tight", &f_zee_e_pt_tight, "f_zee_e_pt_tight/F");
  treeFR->Branch("f_zee_e_eta_tight", &f_zee_e_eta_tight, "f_zee_e_eta_tight/F");
  treeFR->Branch("f_zmumu_mu_pt_tight", &f_zmumu_mu_pt_tight, "f_zmumu_mu_pt_tight/F");
  treeFR->Branch("f_zmumu_mu_eta_tight", &f_zmumu_mu_eta_tight, "f_zmumu_mu_eta_tight/F");
  treeFR->Branch("f_zmumu_e_pt_tight", &f_zmumu_e_pt_tight, "f_zmumu_e_pt_tight/F");
  treeFR->Branch("f_zmumu_e_eta_tight", &f_zmumu_e_eta_tight, "f_zmumu_e_eta_tight/F");
  treeFR->Branch("f_zee_mu_charge_loose", &f_zee_mu_charge_loose, "f_zee_mu_charge_loose/F");
  treeFR->Branch("f_zee_e_charge_loose", &f_zee_e_charge_loose, "f_zee_e_charge_loose/F");
  treeFR->Branch("f_zmumu_mu_charge_loose", &f_zmumu_mu_charge_loose, "f_zmumu_mu_charge_loose/F");
  treeFR->Branch("f_zmumu_e_charge_loose", &f_zmumu_e_charge_loose, "f_zmumu_e_charge_loose/F");
  treeFR->Branch("f_zee_mu_charge_tight", &f_zee_mu_charge_tight, "f_zee_mu_charge_tight/F");
  treeFR->Branch("f_zee_e_charge_tight", &f_zee_e_charge_tight, "f_zee_e_charge_tight/F");
  treeFR->Branch("f_zmumu_mu_charge_tight", &f_zmumu_mu_charge_tight, "f_zmumu_mu_charge_tight/F");
  treeFR->Branch("f_zmumu_e_charge_tight", &f_zmumu_e_charge_tight, "f_zmumu_e_charge_tight/F");
  

  Float_t f_weight2=-999, f_Djet_VAJHU=-999;
  Float_t f_mass4l=-999, f_Z1mass=-999, f_Z2mass=-999;
  Int_t f_run2=-999, f_lumi2=-999, f_event2=-999, f_2p2f=-999, f_3p1f=-999, f_2p2p=-999, f_njets_pass=-999, f_Nbjets=-999;
  Int_t f_lept1_pdgid=-999, f_lept2_pdgid=-999, f_lept3_pdgid=-999, f_lept4_pdgid=-999;
  Int_t f_lept1_pass=-999, f_lept2_pass=-999, f_lept3_pass=-999, f_lept4_pass=-999; 
  Float_t f_lept1_pt=-999, f_lept1_pt_error=-999, f_lept1_eta=-999, f_lept1_phi=-999;
  Float_t f_lept2_pt=-999, f_lept2_pt_error=-999, f_lept2_eta=-999, f_lept2_phi=-999;
  Float_t f_lept3_pt=-999, f_lept3_pt_error=-999, f_lept3_eta=-999, f_lept3_phi=-999;
  Float_t f_lept4_pt=-999, f_lept4_pt_error=-999, f_lept4_eta=-999, f_lept4_phi=-999;
  Float_t f_jets_highpt_pt[3], f_jets_highpt_pt_error[3], f_jets_highpt_eta[3], f_jets_highpt_phi[3], f_jets_highpt_et[3];
  Float_t f_pfmet=-999, f_pfmet_JetEnUp=-999, f_pfmet_JetEnDn=-999, f_pfmet_ElectronEnUp=-999, f_pfmet_ElectronEnDn=-999, f_pfmet_MuonEnUp=-999, f_pfmet_MuonEnDn=-999, f_pfmet_JetResUp=-999, f_pfmet_JetResDn=-999, f_pfmet_UnclusteredEnUp=-999, f_pfmet_UnclusteredEnDn=-999, f_pfmet_PhotonEnUp=-999, f_pfmet_PhotonEnDn=-999;    
  //Float_t f_k41nj2[nShifts], f_k1177nj2[nShifts], f_k2nj3[nShifts], f_k4nj3[nShifts], f_k16nj2e3[nShifts], f_k53nj2e3[nShifts], f_k196nj2e3[nShifts], f_k543nj2e3[nShifts];
  TTree *treeCR = new TTree("ControlRegions","Events in the CRs 2P2F 3P1F");
  treeCR->Branch("f_run", &f_run2, "f_run/I");
  treeCR->Branch("f_lumi", &f_lumi2, "f_lumi/I");    
  treeCR->Branch("f_event", &f_event2, "f_event/I");    
  treeCR->Branch("f_weight", &f_weight2, "f_weight/F");
  treeCR->Branch("f_2p2f", &f_2p2f, "f_2p2f/I");
  treeCR->Branch("f_3p1f", &f_3p1f, "f_3p1f/I");
  treeCR->Branch("f_2p2p", &f_2p2p, "f_2p2p/I");
  treeCR->Branch("f_lept1_pt", &f_lept1_pt, "f_lept1_pt/F");
  treeCR->Branch("f_lept1_pt_error", &f_lept1_pt_error, "f_lept1_pt_error/F");
  treeCR->Branch("f_lept1_eta", &f_lept1_eta, "f_lept1_eta/F");
  treeCR->Branch("f_lept1_phi", &f_lept1_phi, "f_lept1_phi/F");
  treeCR->Branch("f_lept1_pdgid", &f_lept1_pdgid, "f_lept1_pdgid/I");
  treeCR->Branch("f_lept1_pass", &f_lept1_pass, "f_lept1_pass/I");
  treeCR->Branch("f_lept2_pt", &f_lept2_pt, "f_lept2_pt/F");
  treeCR->Branch("f_lept2_pt_error", &f_lept2_pt_error, "f_lept2_pt_error/F");
  treeCR->Branch("f_lept2_eta", &f_lept2_eta, "f_lept2_eta/F");
  treeCR->Branch("f_lept2_phi", &f_lept2_phi, "f_lept2_phi/F");
  treeCR->Branch("f_lept2_pdgid", &f_lept2_pdgid, "f_lept2_pdgid/I");
  treeCR->Branch("f_lept2_pass", &f_lept2_pass, "f_lept2_pass/I");
  treeCR->Branch("f_lept3_pt", &f_lept3_pt, "f_lept3_pt/F");
  treeCR->Branch("f_lept3_pt_error", &f_lept3_pt_error, "f_lept3_pt_error/F");
  treeCR->Branch("f_lept3_eta", &f_lept3_eta, "f_lept3_eta/F");
  treeCR->Branch("f_lept3_phi", &f_lept3_phi, "f_lept3_phi/F");
  treeCR->Branch("f_lept3_pdgid", &f_lept3_pdgid, "f_lept3_pdgid/I");
  treeCR->Branch("f_lept3_pass", &f_lept3_pass, "f_lept3_pass/I");
  treeCR->Branch("f_lept4_pt", &f_lept4_pt, "f_lept4_pt/F");
  treeCR->Branch("f_lept4_pt_error", &f_lept4_pt_error, "f_lept4_pt_error/F");
  treeCR->Branch("f_lept4_eta", &f_lept4_eta, "f_lept4_eta/F");
  treeCR->Branch("f_lept4_phi", &f_lept4_phi, "f_lept4_phi/F");
  treeCR->Branch("f_lept4_pdgid", &f_lept4_pdgid, "f_lept4_pdgid/I");
  treeCR->Branch("f_lept4_pass", &f_lept4_pass, "f_lept4_pass/I");
  treeCR->Branch("f_mass4l", &f_mass4l, "f_mass4l/F");
  treeCR->Branch("f_Z1mass", &f_Z1mass, "f_Z1mass/F");
  treeCR->Branch("f_Z2mass", &f_Z2mass, "f_Z2mass/F");
  treeCR->Branch("f_njets_pass", &f_njets_pass, "f_njets_pass/I");
  treeCR->Branch("f_pfmet", &f_pfmet, "f_pfmet/F");
  treeCR->Branch("f_pfmet_JetEnUp", &f_pfmet_JetEnUp, "f_pfmet_JetEnUp/F");
  treeCR->Branch("f_pfmet_JetEnDn", &f_pfmet_JetEnDn, "f_pfmet_JetEnDn/F");
  treeCR->Branch("f_pfmet_ElectronEnUp", &f_pfmet_ElectronEnUp, "f_pfmet_ElectronEnUp/F");
  treeCR->Branch("f_pfmet_ElectronEnDn", &f_pfmet_ElectronEnDn, "f_pfmet_ElectronEnDn/F");
  treeCR->Branch("f_pfmet_MuonEnUp", &f_pfmet_MuonEnUp, "f_pfmet_MuonEnUp/F");
  treeCR->Branch("f_pfmet_MuonEnDn", &f_pfmet_MuonEnDn, "f_pfmet_MuonEnDn/F");
  treeCR->Branch("f_pfmet_JetResUp", &f_pfmet_JetResUp, "f_pfmet_JetResUp/F");
  treeCR->Branch("f_pfmet_JetResDn", &f_pfmet_JetResDn, "f_pfmet_JetResDn/F");
  treeCR->Branch("f_pfmet_UnclusteredEnUp", &f_pfmet_UnclusteredEnUp, "f_pfmet_UnclusteredEnUp/F");
  treeCR->Branch("f_pfmet_UnclusteredEnDn", &f_pfmet_UnclusteredEnDn, "f_pfmet_UnclusteredEnDn/F");
  treeCR->Branch("f_pfmet_PhotonEnUp", &f_pfmet_PhotonEnUp, "f_pfmet_PhotonEnUp/F");
  treeCR->Branch("f_pfmet_PhotonEnDn", &f_pfmet_PhotonEnDn, "f_pfmet_PhotonEnDn/F");
  treeCR->Branch("f_Nbjets", &f_Nbjets, "f_Nbjets/I");
  treeCR->Branch("f_jets_highpt_pt", &f_jets_highpt_pt, "f_jets_highpt_pt[3]/F");
  treeCR->Branch("f_jets_highpt_pt_error", &f_jets_highpt_pt_error, "f_jets_highpt_pt_error[3]/F");
  treeCR->Branch("f_jets_highpt_eta", &f_jets_highpt_eta, "f_jets_highpt_eta[3]/F");
  treeCR->Branch("f_jets_highpt_phi", &f_jets_highpt_phi, "f_jets_highpt_phi[3]/F");
  treeCR->Branch("f_jets_highpt_et", &f_jets_highpt_et, "f_jets_highpt_et[3]/F");
  treeCR->Branch("f_Djet_VAJHU", &f_Djet_VAJHU, "f_Djet_VAJHU/F");
  /*
  treeCR->Branch("f_k41nj2", &f_k41nj2, Form("f_k41nj2[%i]/F",nShifts));
  treeCR->Branch("f_k1177nj2", &f_k1177nj2, Form("f_k1177nj2[%i]/F",nShifts));
  treeCR->Branch("f_k2nj3", &f_k2nj3, Form("f_k2nj3[%i]/F",nShifts));
  treeCR->Branch("f_k4nj3", &f_k4nj3, Form("f_k4nj3[%i]/F",nShifts));
  treeCR->Branch("f_k16nj2e3", &f_k16nj2e3, Form("f_k16nj2e3[%i]/F",nShifts));
  treeCR->Branch("f_k53nj2e3", &f_k53nj2e3, Form("f_k53nj2e3[%i]/F",nShifts));
  treeCR->Branch("f_k196nj2e3", &f_k196nj2e3, Form("f_k196nj2e3[%i]/F",nShifts));
  treeCR->Branch("f_k543nj2e3", &f_k543nj2e3, Form("f_k543nj2e3[%i]/F",nShifts));
  */

  
  
  //======================================================================
  //  Initializing variables and histograms and call functions goes here //
  //========================================================================
  // Declare MEM class //////////////////////////check //////////////////////////////////////
  //  MEMs combinedMEM(13,125,"CTEQ6L");
  
  // MuonCalibrator
  //  KalmanMuonCalibrator calibrator("DATA_76X_13TeV");
  
  //Define functions
  //-----------------
  
  // Declare MEM class
  MEMs combinedMEM(13,125,"CTEQ6L");  

  
  double DELTAPHI( double , double ) ; //call the function 
  //double EAele(int ,bool ); //must be defined in header file
  
  /////////// Pileup reweighting 2015 data vs Fall15 MC in 76x ///////////////////////
  TFile *_filePU;
  _filePU= TFile::Open(auxiliar_files+"/puWeightsMoriond17_v2.root");
  TH1D *puweight = (TH1D*)_filePU->Get("weights");
  //_filePU= TFile::Open(auxiliar_files+"/pileup_MC_80x_271036-276811_69200.root");
  //TH1D *puweight = (TH1D*)_filePU->Get("puweight");
  
  
  
   //////////////////Lepton Efficiency Scale Factrons//////////////////////////////
  
   // Load histograms
   //=================

   //electron scale factor 2015 data 

   // TFile *ele_scale_factors2015 = new TFile("IdIsoSip.root");
   // TH2F *ele_scale_2015 = (TH2F*)gDirectory->Get("hScaleFactors_IdIsoSip");
  
  TFile *ele_scale_factors_v3 = new TFile(auxiliar_files+"/ele_scale_factors_v3.root");
  TH2F *ele_scale_factors2016 = (TH2F*)gDirectory->Get("ele_scale_factors"); 
  // TFile *ele_scale_IdIsoSipfactors2015_Crack = new TFile("IdIsoSip_Cracks.root");
  TH2F *ele_scale_factors_gap2016 = (TH2F*)gDirectory->Get("ele_scale_factors_gap"); 
  
  //muon scale factor 2015 data
  
  // TFile *mu_scale_factors = new TFile("final_HZZ.root");
  // TH2F *mu_scale_2015 = (TH2F*)gDirectory->Get("FINAL");
  
  TFile *mu_scale_factors = new TFile(auxiliar_files+"/final_HZZSF_pt0_200.root");
  TH2F *mu_scale_2016 = (TH2F*)gDirectory->Get("FINAL"); 
  
  

  ///////////////////// ZZ K factor //////////////////////////////////////////
  
  // kfactor_ggZZ(float GENmassZZ, int finalState)     
  TString strSystTitle[9] ={
    "Nominal",
    "PDFScaleDn",
    "PDFScaleUp",
    "QCDScaleDn",
    "QCDScaleUp",
    "AsDn",
    "AsUp",
    "PDFReplicaDn",
    "PDFReplicaUp"
  };
  
  TFile* fin = TFile::Open(auxiliar_files+"/Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
  // Open the files
  TSpline3* ggZZ_kf[9]={NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  for(int f=0;f<9;f++){
    ggZZ_kf[f] = (TSpline3*)fin->Get(Form("sp_kfactor_%s", strSystTitle[f].Data()));
  }   
  fin->Close();
  
  /////////////////////////////////////////////////////////////////////////////////////
  
  // Book root file (for output):
  
  //TFile * theFile = new TFile(output,"RECREATE");
  
  //_______________________//
  //Initializing variables //
  //______________________//
  
  //counters
  //========
  //1)DATA
  //------
  
  int event=0;
  
  int N_0 = 0;
  
  
  int N_2 = 0; //number of events after good muon cut
  int N_3_FSR = 0; //number of events with FSR photons
  int N_3a = 0; //number of events passing cut number of Z candidates => 2
  int N_3b = 0; //number of events after 2 cuts a)cut on mass of all Zs to be >12 || <120 , b) if number of Zs afer mass cut < 2 escape event
  
  int  N_Z1_step = 0;
  
  int R_AI_AI_MMMM =0;
  int R_AI_I_MMMM =0;
  int R_I_AI_MMMM =0;
  int R_12_MMMM =0;
  
  int R_AI_AI_MMEE =0;
  int R_AI_I_MMEE =0;
  int R_I_AI_MMEE =0;
  int R_12_MMEE =0;
  
  int R_AI_AI_EEEE =0;
  int R_AI_I_EEEE =0;
  int R_I_AI_EEEE =0;
  int R_12_EEEE =0;
  
  int R_AI_AI_EEMM =0;
  int R_AI_I_EEMM =0;
  int R_I_AI_EEMM =0;
  int R_12_EEMM =0;
  
  
  int R_AI_AI_MMMM_w =0;
  int R_AI_I_MMMM_w =0;
  int R_I_AI_MMMM_w =0;
  int R_12_MMMM_w =0;
  
  int R_AI_AI_MMEE_w =0;
  int R_AI_I_MMEE_w =0;
  int R_I_AI_MMEE_w =0;
  int R_12_MMEE_w =0;
  
  int R_AI_AI_EEEE_w =0;
  int R_AI_I_EEEE_w =0;
  int R_I_AI_EEEE_w =0;
  int R_12_EEEE_w =0;
  
  int R_AI_AI_EEMM_w =0;
  int R_AI_I_EEMM_w =0;
  int R_I_AI_EEMM_w =0;
  int R_12_EEMM_w =0;
  
  
  //2)MC
  //-----
  int N_0_w = 0;
  
  int N_2_w = 0;
  int N_3_FSR_w = 0;
  int N_3a_w = 0;
  int N_3b_w = 0;
  
  //3)weight
  //---------
  
  float newweight_3=1;
  float newweight=1; // =1 for now untill calculate weight
  //Double_t eff_weight = 1.;
  
  
  
  bool debug=false;  //debug flag  -- default false
  
  
  //_______________________//
  //Initializing histograms
  //________________________//
  
  //Number of events
  //===================
  TH1F *histo_EVENT        = new TH1F("histo_EVENT", "histo_EVENT",600,0,600.);
  
  // Pileup reweighting
 //========================
  
  TH1F *hPUvertices             = new TH1F("hPUvertices", "hPUvertices",70,0.,70.);  
  TH1F *hPUvertices_ReWeighted  = new TH1F("hPUvertices_ReWeighted", "hPUvertices_ReWeighted",70,0.,70.);  
  
  //loose lepton identification
  //==============================
  // 1)loose muon 
  
  TH1F * hN_loose_mu = new TH1F("hN_loose_mu", "Number of muons before loose cut", 30 , 0. , 30. );//sum of loose muons in all events after loss id (each event contain 1  or 2 or 3 muons ...)
  hN_loose_mu->SetXTitle("N_loose_mu");
  TH1F * hIso_loose_mu = new TH1F("hIso_loose_mu", "Max value of Iso. for loose mu in each event (worth isolation)", 2000 , -10. , 10. ); //(each event contain 1  or 2 or 3 muons ...)
  hN_loose_mu->SetXTitle("Iso_loose_mu");
  TH1F * hSip_loose_mu = new TH1F("hSip_loose_mu", "Max value of SIP for loose mu in each event (worth SIP)", 1000 , -20. , 40. );
  hN_loose_mu->SetXTitle("Sip_loose_mu");
  TH1F * hIp_loose_mu = new TH1F("hIp_loose_mu", "Max value of IP for loose mu in each event (worth IP)", 1000 , -20. , 40. );
  hN_loose_mu->SetXTitle("Ip_loose_mu");
  
  // 2)loose ele
  
  TH1F * hN_loose_e = new TH1F("hN_loose_e", "Number of electrons before loose cut", 30 , 0. , 30. );
  hN_loose_mu->SetXTitle("N_loose_mu");
  TH1F * hIso_loose_e = new TH1F("hIso_loose_e", "Max value of Iso. for loose ele in each event(worth isolation)", 2000 , -10. , 10. ); //(each event contain 1  or 2 or 3 eles ...)
  hN_loose_e->SetXTitle("Iso_loose_e");
  TH1F * hSip_loose_e = new TH1F("hSip_loose_e", "Max value of SIP for loose electron in each event (worth SIP)", 1000 , -20. , 40. );
  hN_loose_e->SetXTitle("Sip_loose_e");
  TH1F * hIp_loose_e = new TH1F("hIp_loose_e", "Max value of IP for loose electron in each event (worth IP)", 1000 , -20. , 40. );
  hN_loose_e->SetXTitle("Ip_loose_e");
  
  // good leptons
  //=================
  //1) good muons
  TH1F * hN_good_mu = new TH1F("hN_good_mu", "N_good_mu", 30 , 0. , 30. );
  hN_good_mu->SetXTitle("N_good_mu");
  
  //2)good electrons
  TH1F * hN_good_ele = new TH1F("hN_good_ele", "N_good_ele", 30 , 0. , 30. );
  hN_good_ele->SetXTitle("N_good_ele");
  
  TH1F * hN_good_lep = new TH1F("hN_good_lep", "N_good_lep", 30 , 0. , 30. );
  hN_good_lep->SetXTitle("N_good_lep");
  
  //3)good photons
  
  TH1F * hN_good_phot = new TH1F("hN_good_phot", "N_good_phot", 30 , 0. , 30. );
  hN_good_phot->SetXTitle("N_good_photon");
  
  
  //Step 3
  //-------
  //after choosing Z1 and before applying FR
  
  //mumu
  
  TH1F * hMZ1_3_mumu = new TH1F("hMZ1_3_mumu", "Mass of Z1_mumu after selection step 3", 200 , -0.5 , 199.5 );
  hMZ1_3_mumu->SetXTitle("mass_Z1  (GeV)");
  TH1F * hPtZ1_3_mumu = new TH1F("hPtZ1_3_mumu", "Pt of Z1_mumu after selection step 3", 200 , -0.5 , 199.5 );
  hPtZ1_3_mumu->SetXTitle("pt_Z1  (GeV)");
  TH1F * hYZ1_3_mumu = new TH1F("hYZ1_3_mumu", "Y of Z1_mumu after selection step 3", 500 , -5. , 5.);
  hYZ1_3_mumu->SetXTitle("Y_Z1");
  
  //ee
  
  TH1F * hMZ1_3_ee = new TH1F("hMZ1_3_ee", "Mass of Z1_ee after selection step 3", 200 , -0.5 , 199.5 );
  hMZ1_3_ee->SetXTitle("mass_Z1  (GeV)");
  TH1F * hPtZ1_3_ee = new TH1F("hPtZ1_3_ee", "Pt of Z1_ee after selection step 3", 200 , -0.5 , 199.5 );
  hPtZ1_3_ee->SetXTitle("pt_Z1_ee  (GeV)");
  TH1F * hYZ1_3_ee = new TH1F("hYZ1_3_ee", "Y of Z1_ee after selection step 3", 500 , -5. , 5.);
  hYZ1_3_ee->SetXTitle("Y_Z1");
   
  
  
  TH1F * hMZ1_3 = new TH1F("hMZ1_3", "Mass of Z1 after selection step 3", 200 , -0.5 , 199.5 );
  hMZ1_3->SetXTitle("mass_Z1  (GeV)");
  TH1F * hPtZ1_3 = new TH1F("hPtZ1_3", "Pt of Z1 after selection step 3", 200 , -0.5 , 199.5 );
  hPtZ1_3->SetXTitle("pt_Z1  (GeV)");
  TH1F * hPFMET_3 = new TH1F("hPFMET_3", "PF MET after selection step 3", 1000 , 0., 1000.);
  hPFMET_3->SetXTitle("PFMET");
  
  //step 6
  //QCD suppresion
  //minimum mass of opposite charge same flavour leptons 
  //this filled before apply cut mll>4
  
  TH1F * hminMll_6 = new TH1F("hminMll_6", "minMll at selection step 6", 400 , 0. , 200.);
  hminMll_6->SetXTitle("minMll  (GeV)");
  
  
  //Fake rate (after MET cut)
  //===========
  //1) Muon fake rate Histos
  //-------------------------
  TH1F *ZplusM_Pt_DEN_Barrel   = new TH1F("ZplusM_Pt_DEN_Barrel", "ZplusM_Pt_DEN_Barrel",150,0.,150.);
  ZplusM_Pt_DEN_Barrel->SetXTitle("pT");
  TH1F *ZplusM_Pt_DEN_Endcaps  = new TH1F("ZplusM_Pt_DEN_Endcaps", "ZplusM_Pt_DEN_Endcaps",150,0.,150.);
  ZplusM_Pt_DEN_Endcaps->SetXTitle("pT");
  TH1F *ZplusM_Pt_NUM_ID_ISO_Barrel = new TH1F("ZplusM_Pt_NUM_ID_ISO_Barrel", "ZplusM_Pt_NUM_ID_ISO_Barrel",150,0.,150.);
  ZplusM_Pt_NUM_ID_ISO_Barrel->SetXTitle("pT");
  TH1F *ZplusM_Pt_NUM_ID_ISO_Endcaps = new TH1F("ZplusM_Pt_NUM_ID_ISO_Endcaps", "ZplusM_Pt_NUM_ID_ISO_Endcaps",150,0.,150.);
  ZplusM_Pt_NUM_ID_ISO_Endcaps->SetXTitle("pT");

  //Positive

  TH1F *ZplusM_Pt_DEN_Barrel_pos   = new TH1F("ZplusM_Pt_DEN_Barrel_pos", "ZplusM_Pt_DEN_Barrel_pos",150,0.,150.);
  ZplusM_Pt_DEN_Barrel_pos->SetXTitle("pT");
  TH1F *ZplusM_Pt_DEN_Endcaps_pos  = new TH1F("ZplusM_Pt_DEN_Endcaps_pos", "ZplusM_Pt_DEN_Endcaps_pos",150,0.,150.);
  ZplusM_Pt_DEN_Endcaps_pos->SetXTitle("pT");
  TH1F *ZplusM_Pt_NUM_ID_ISO_Barrel_pos = new TH1F("ZplusM_Pt_NUM_ID_ISO_Barrel_pos", "ZplusM_Pt_NUM_ID_ISO_Barrel_pos",150,0.,150.);
  ZplusM_Pt_NUM_ID_ISO_Barrel_pos->SetXTitle("pT");
  TH1F *ZplusM_Pt_NUM_ID_ISO_Endcaps_pos = new TH1F("ZplusM_Pt_NUM_ID_ISO_Endcaps_pos", "ZplusM_Pt_NUM_ID_ISO_Endcaps_pos",150,0.,150.);
  ZplusM_Pt_NUM_ID_ISO_Endcaps_pos->SetXTitle("pT");

  //Negative

  TH1F *ZplusM_Pt_DEN_Barrel_neg   = new TH1F("ZplusM_Pt_DEN_Barrel_neg", "ZplusM_Pt_DEN_Barrel_neg",150,0.,150.);
  ZplusM_Pt_DEN_Barrel_neg->SetXTitle("pT");
  TH1F *ZplusM_Pt_DEN_Endcaps_neg  = new TH1F("ZplusM_Pt_DEN_Endcaps_neg", "ZplusM_Pt_DEN_Endcaps_neg",150,0.,150.);
  ZplusM_Pt_DEN_Endcaps_neg->SetXTitle("pT");
  TH1F *ZplusM_Pt_NUM_ID_ISO_Barrel_neg = new TH1F("ZplusM_Pt_NUM_ID_ISO_Barrel_neg", "ZplusM_Pt_NUM_ID_ISO_Barrel_neg",150,0.,150.);
  ZplusM_Pt_NUM_ID_ISO_Barrel_neg->SetXTitle("pT");
  TH1F *ZplusM_Pt_NUM_ID_ISO_Endcaps_neg = new TH1F("ZplusM_Pt_NUM_ID_ISO_Endcaps_neg", "ZplusM_Pt_NUM_ID_ISO_Endcaps_neg",150,0.,150.);
  ZplusM_Pt_NUM_ID_ISO_Endcaps_neg->SetXTitle("pT");
  
  //1) electron fake rate Histos
  //-------------------------------
  
  TH1F *ZplusE_Pt_DEN_Barrel   = new TH1F("ZplusE_Pt_DEN_Barrel", "ZplusE_Pt_DEN_Barrel",150,0.,150.);
  ZplusE_Pt_DEN_Barrel->SetXTitle("pT");
  TH1F *ZplusE_Pt_DEN_Endcaps  = new TH1F("ZplusE_Pt_DEN_Endcaps", "ZplusE_Pt_DEN_Endcaps",150,0.,150.);
  ZplusE_Pt_DEN_Endcaps->SetXTitle("pT");
  TH1F *ZplusE_Pt_NUM_ID_ISO_Barrel = new TH1F("ZplusE_Pt_NUM_ID_ISO_Barrel", "ZplusE_Pt_NUM_ID_ISO_Barrel",150,0.,150.);
  ZplusE_Pt_NUM_ID_ISO_Barrel->SetXTitle("pT");
  TH1F *ZplusE_Pt_NUM_ID_ISO_Endcaps = new TH1F("ZplusE_Pt_NUM_ID_ISO_Endcaps", "ZplusE_Pt_NUM_ID_ISO_Endcaps",150,0.,150.);
  ZplusE_Pt_NUM_ID_ISO_Endcaps->SetXTitle("pT");

  //postive

  TH1F *ZplusE_Pt_DEN_Barrel_pos   = new TH1F("ZplusE_Pt_DEN_Barrel_pos", "ZplusE_Pt_DEN_Barrel_pos",150,0.,150.);
  ZplusE_Pt_DEN_Barrel_pos->SetXTitle("pT");
  TH1F *ZplusE_Pt_DEN_Endcaps_pos  = new TH1F("ZplusE_Pt_DEN_Endcaps_pos", "ZplusE_Pt_DEN_Endcaps_pos",150,0.,150.);
  ZplusE_Pt_DEN_Endcaps_pos->SetXTitle("pT");
  TH1F *ZplusE_Pt_NUM_ID_ISO_Barrel_pos = new TH1F("ZplusE_Pt_NUM_ID_ISO_Barrel_pos", "ZplusE_Pt_NUM_ID_ISO_Barrel_pos",150,0.,150.);
  ZplusE_Pt_NUM_ID_ISO_Barrel_pos->SetXTitle("pT");
  TH1F *ZplusE_Pt_NUM_ID_ISO_Endcaps_pos = new TH1F("ZplusE_Pt_NUM_ID_ISO_Endcaps_pos", "ZplusE_Pt_NUM_ID_ISO_Endcaps_pos",150,0.,150.);
  ZplusE_Pt_NUM_ID_ISO_Endcaps_pos->SetXTitle("pT");

  //Negative
  
  TH1F *ZplusE_Pt_DEN_Barrel_neg   = new TH1F("ZplusE_Pt_DEN_Barrel_neg", "ZplusE_Pt_DEN_Barrel_neg",150,0.,150.);
  ZplusE_Pt_DEN_Barrel_neg->SetXTitle("pT");
  TH1F *ZplusE_Pt_DEN_Endcaps_neg  = new TH1F("ZplusE_Pt_DEN_Endcaps_neg", "ZplusE_Pt_DEN_Endcaps_neg",150,0.,150.);
  ZplusE_Pt_DEN_Endcaps->SetXTitle("pT");
  TH1F *ZplusE_Pt_NUM_ID_ISO_Barrel_neg = new TH1F("ZplusE_Pt_NUM_ID_ISO_Barrel_neg", "ZplusE_Pt_NUM_ID_ISO_Barrel_neg",150,0.,150.);
  ZplusE_Pt_NUM_ID_ISO_Barrel->SetXTitle("pT");
  TH1F *ZplusE_Pt_NUM_ID_ISO_Endcaps_neg = new TH1F("ZplusE_Pt_NUM_ID_ISO_Endcaps_neg", "ZplusE_Pt_NUM_ID_ISO_Endcaps_neg",150,0.,150.);
  ZplusE_Pt_NUM_ID_ISO_Endcaps_neg->SetXTitle("pT");
  
 ////
  //muons  //just after Z cut , No MET cut
  
  TH1F * hMZ1_loose_mu_barrel = new TH1F("hMZ1_loose_mu_barrel", "Mass of Z1 after selection loose muon barrel", 200 , -0.5 , 199.5 );
  hMZ1_loose_mu_barrel->SetXTitle("mass_Z1  (GeV)");
  TH1F * hPFMET_loose_mu_barrel = new TH1F("hPFMET_loose_mu_barrel", "PF MET after selection loose muon barrel", 1000 , 0., 1000.);
  hPFMET_loose_mu_barrel->SetXTitle("PFMET");
  TH1F * hPT_loose_mu_barrel= new TH1F("hPT_loose_mu_barrel", "hPT_loose_mu_barrel",150,0.,150.); 
  hPT_loose_mu_barrel->SetXTitle("pT");

  TH1F * h_PT_loose_mu_barrel_pos= new TH1F("h_PT_loose_mu_barrel_pos", "h_PT_loose_mu_barrel_pos",150,0.,150.); 
  TH1F * h_PFMET_loose_mu_barrel_pos= new TH1F("h_PFMET_loose_mu_barrel_pos", "h_PFMET_loose_mu_barrel_pos", 1000 , 0., 1000.);
  TH1F * h_PT_loose_mu_barrel_neg= new TH1F("h_PT_loose_mu_barrel_neg", "h_PT_loose_mu_barrel_neg",150,0.,150.); 
  TH1F * h_PFMET_loose_mu_barrel_neg= new TH1F("h_PFMET_loose_mu_barrel_neg", "h_PFMET_loose_mu_barrel_neg", 1000 , 0., 1000.);
  
  
  TH1F * hMZ1_loose_mu_endcap = new TH1F("hMZ1_loose_mu_endcap", "Mass of Z1 after selection loose muon_endcap", 200 , -0.5 , 199.5 );
  hMZ1_loose_mu_endcap->SetXTitle("mass_Z1  (GeV)");
  TH1F * hPFMET_loose_mu_endcap = new TH1F("hPFMET_loose_mu_endcap", "PF MET after selection loose muon_endcap", 1000 , 0., 1000.);
  hPFMET_loose_mu_endcap->SetXTitle("PFMET");
  TH1F * hPT_loose_mu_endcap= new TH1F("hPT_loose_mu_endcap", "hPT_loose_mu_endcap",150,0.,150.); 
  hPT_loose_mu_endcap->SetXTitle("pT");

  TH1F * h_PT_loose_mu_endcap_pos= new TH1F("h_PT_loose_mu_endcap_pos", "h_PT_loose_mu_endcap_pos",150,0.,150.); 
  TH1F * h_PFMET_loose_mu_endcap_pos= new TH1F("h_PFMET_loose_mu_endcap_pos", "h_PFMET_loose_mu_endcap_pos", 1000 , 0., 1000.);
  TH1F * h_PT_loose_mu_endcap_neg= new TH1F("h_PT_loose_mu_endcap_neg", "h_PT_loose_mu_endcap_neg",150,0.,150.); 
  TH1F * h_PFMET_loose_mu_endcap_neg= new TH1F("h_PFMET_loose_mu_endcap_neg", "h_PFMET_loose_mu_endcap_neg", 1000 , 0., 1000.);

  
  TH1F * hMZ1_tight_mu_barrel = new TH1F("hMZ1_tight_mu_barrel", "Mass of Z1 after selection tight muon_barrel", 200 , -0.5 , 199.5 );
  hMZ1_tight_mu_barrel->SetXTitle("mass_Z1  (GeV)");
  TH1F * hPFMET_tight_mu_barrel = new TH1F("hPFMET_tight_mu_barrel", "PF MET after selection tight muon_barrel", 1000 , 0., 1000.);
  hPFMET_tight_mu_barrel->SetXTitle("PFMET");   
  TH1F * hPT_tight_mu_barrel= new TH1F("hPT_tight_mu_barrel", "hPT_tight_mu_barrel",150,0.,150.); 
  hPT_tight_mu_barrel->SetXTitle("pT");

  TH1F * h_PT_tight_mu_barrel_pos= new TH1F("h_PT_tight_mu_barrel_pos", "h_PT_tight_mu_barrel_pos",150,0.,150.); 
  TH1F * h_PFMET_tight_mu_barrel_pos= new TH1F("h_PFMET_tight_mu_barrel_pos", "h_PFMET_tight_mu_barrel_pos", 1000 , 0., 1000.);
  TH1F * h_PT_tight_mu_barrel_neg= new TH1F("h_PT_tight_mu_barrel_neg", "h_PT_tight_mu_barrel_neg",150,0.,150.); 
  TH1F * h_PFMET_tight_mu_barrel_neg= new TH1F("h_PFMET_tight_mu_barrel_neg", "h_PFMET_tight_mu_barrel_neg", 1000 , 0., 1000.); 
  
  TH1F * hMZ1_tight_mu_endcap = new TH1F("hMZ1_tight_mu_endcap", "Mass of Z1 after selection tight muon_endcap", 200 , -0.5 , 199.5 );
  hMZ1_tight_mu_endcap->SetXTitle("mass_Z1  (GeV)");   
  TH1F * hPFMET_tight_mu_endcap = new TH1F("hPFMET_tight_mu_endcap", "PF MET after selection tight muon_endcap", 1000 , 0., 1000.);
  hPFMET_tight_mu_endcap->SetXTitle("PFMET");
  TH1F * hPT_tight_mu_endcap= new TH1F("hPT_tight_mu_endcap", "hPT_tight_mu_endcap",150,0.,150.); 
  hPT_tight_mu_endcap->SetXTitle("pT");

  TH1F * h_PT_tight_mu_endcap_pos= new TH1F("h_PT_tight_mu_endcap_pos", "h_PT_tight_mu_endcap_pos",150,0.,150.); 
  TH1F * h_PFMET_tight_mu_endcap_pos= new TH1F("h_PFMET_tight_mu_endcap_pos", "h_PFMET_tight_mu_endcap_pos", 1000 , 0., 1000.);
  TH1F * h_PT_tight_mu_endcap_neg= new TH1F("h_PT_tight_mu_endcap_neg", "h_PT_tight_mu_endcap_neg",150,0.,150.); 
  TH1F * h_PFMET_tight_mu_endcap_neg= new TH1F("h_PFMET_tight_mu_endcap_neg", "h_tight_loose_mu_endcap_neg", 1000 , 0., 1000.);
  
  
  
   ///electrons
   
   TH1F * hMZ1_loose_ele_barrel = new TH1F("hMZ1_loose_ele_barrel", "Mass of Z1 after selection loose ele barrel", 200 , -0.5 , 199.5 );
   hMZ1_loose_ele_barrel->SetXTitle("mass_Z1  (GeV)");   
   TH1F * hPFMET_loose_ele_barrel = new TH1F("hPFMET_loose_ele_barrel", "PF MET after selection loose ele barrel", 1000 , 0., 1000.);
   hPFMET_loose_ele_barrel->SetXTitle("PFMET");
   TH1F * hPT_loose_ele_barrel= new TH1F("hPT_loose_ele_barrel", "hPT_loose_ele_barrel",150,0.,150.); 
   hPT_loose_ele_barrel->SetXTitle("pT");

   TH1F * h_PT_loose_ele_barrel_pos= new TH1F("h_PT_loose_ele_barrel_pos", "h_PT_loose_ele_barrel_pos",150,0.,150.); 
   TH1F * h_PFMET_loose_ele_barrel_pos= new TH1F("h_PFMET_loose_ele_barrel_pos", "h_PFMET_loose_ele_barrel_pos", 1000 , 0., 1000.);
   TH1F * h_PT_loose_ele_barrel_neg= new TH1F("h_PT_loose_ele_barrel_neg", "h_PT_loose_ele_barrel_neg",150,0.,150.); 
   TH1F * h_PFMET_loose_ele_barrel_neg= new TH1F("h_PFMET_loose_ele_barrel_neg", "h_PFMET_loose_ele_barrel_neg", 1000 , 0., 1000.);
  
   
   TH1F * hMZ1_loose_ele_endcap = new TH1F("hMZ1_loose_ele_endcap", "Mass of Z1 after selection loose ele endcap", 200 , -0.5 , 199.5 );
   hMZ1_loose_ele_endcap->SetXTitle("mass_Z1  (GeV)");   
   TH1F * hPFMET_loose_ele_endcap = new TH1F("hPFMET_loose_ele_endcap", "PF MET after selection loose ele endcap", 1000 , 0., 1000.);
   hPFMET_loose_ele_endcap->SetXTitle("PFMET");
   TH1F * hPT_loose_ele_endcap= new TH1F("hPT_loose_ele_endcap", "hPT_loose_ele_endcap",150,0.,150.); 
   hPT_loose_ele_endcap->SetXTitle("pT");

   TH1F * h_PT_loose_ele_endcap_pos= new TH1F("h_PT_loose_ele_endcap_pos", "h_PT_loose_ele_endcap_pos",150,0.,150.); 
   TH1F * h_PFMET_loose_ele_endcap_pos= new TH1F("h_PFMET_loose_ele_endcap_pos", "h_PFMET_loose_ele_endcap_pos", 1000 , 0., 1000.);
   TH1F * h_PT_loose_ele_endcap_neg= new TH1F("h_PT_loose_ele_endcap_neg", "h_PT_loose_ele_endcap_neg",150,0.,150.); 
   TH1F * h_PFMET_loose_ele_endcap_neg= new TH1F("h_PFMET_loose_ele_endcap_neg", "h_PFMET_loose_ele_endcap_neg", 1000 , 0., 1000.);
      
   TH1F * hMZ1_tight_ele_barrel = new TH1F("hMZ1_tight_ele_barrel", "Mass of Z1 after selection tight ele barrel", 200 , -0.5 , 199.5 );
   hMZ1_tight_ele_barrel->SetXTitle("mass_Z1  (GeV)");   
   TH1F * hPFMET_tight_ele_barrel = new TH1F("hPFMET_tight_ele_barrel", "PF MET after selection tight ele barrel", 1000 , 0., 1000.);
   hPFMET_tight_ele_barrel->SetXTitle("PFMET");
   TH1F * hPT_tight_ele_barrel= new TH1F("hPT_tight_ele_barrel", "hPT_tight_ele_barrel",150,0.,150.); 
   hPT_tight_ele_barrel->SetXTitle("pT");

   TH1F * h_PT_tight_ele_barrel_pos= new TH1F("h_PT_tight_ele_barrel_pos", "h_PT_tight_ele_barrel_pos",150,0.,150.); 
   TH1F * h_PFMET_tight_ele_barrel_pos= new TH1F("h_PFMET_tight_ele_barrel_pos", "h_PFMET_tight_ele_barrel_pos", 1000 , 0., 1000.);
   TH1F * h_PT_tight_ele_barrel_neg= new TH1F("h_PT_tight_ele_barrel_neg", "h_PT_tight_ele_barrel_neg",150,0.,150.); 
   TH1F * h_PFMET_tight_ele_barrel_neg= new TH1F("h_PFMET_tight_ele_barrel_neg", "h_PFMET_tight_ele_barrel_neg", 1000 , 0., 1000.);
   
   TH1F * hMZ1_tight_ele_endcap = new TH1F("hMZ1_tight_ele_endcap", "Mass of Z1 after selection tight ele endcap", 200 , -0.5 , 199.5 );
   hMZ1_tight_ele_endcap->SetXTitle("mass_Z1  (GeV)");  
   TH1F * hPFMET_tight_ele_endcap = new TH1F("hPFMET_tight_ele_endcap", "PF MET after selection tight ele endcap", 1000 , 0., 1000.);
   hPFMET_tight_ele_endcap->SetXTitle("PFMET");
   TH1F * hPT_tight_ele_endcap= new TH1F("hPT_tight_ele_endcap", "hPT_tight_ele_endcap",150,0.,150.); 
   hPT_tight_ele_endcap->SetXTitle("pT");

   TH1F * h_PT_tight_ele_endcap_pos= new TH1F("h_PT_tight_ele_endcap_pos", "h_PT_tight_ele_endcap_pos",150,0.,150.); 
   TH1F * h_PFMET_tight_ele_endcap_pos= new TH1F("h_PFMET_tight_ele_endcap_pos", "h_PFMET_tight_ele_endcap_pos", 1000 , 0., 1000.);
   TH1F * h_PT_tight_ele_endcap_neg= new TH1F("h_PT_tight_ele_endcap_neg", "h_PT_tight_ele_endcap_neg",150,0.,150.); 
   TH1F * h_PFMET_tight_ele_endcap_neg= new TH1F("h_PFMET_tight_ele_endcap_neg", "h_PFMET_tight_ele_endcap_neg", 1000 , 0., 1000.);
   
   
   //For check
   
   TH1F * hMZ_3_check = new TH1F("hMZ_3_check", "Mass of Z after selection step 3", 200 , -0.5 , 199.5 );
   hMZ_3_check->SetXTitle("mass_Z  (GeV)");   
   TH1F *  hMZ1_3_no_effW = new TH1F("hMZ1_3_no_effW", "Mass of Z after selection step 3", 200 , -0.5 , 199.5 ); 
   TH1F * hPFMET_3_no_effW = new TH1F("hPFMET_3_no_effW", "PF MET after selection step 3", 1000 , 0., 1000.);
   
   //2D histograms

   TH2F  *MET_PT_Ele_Barrel_Den  = new TH2F ("MET_PT_Ele_Barrel_Den" , "MET_PT_Ele_Barrel_Den", 1000,0,1000,1000,0,1000);
   TH2F  *MET_PT_Ele_Barrel_Num  = new TH2F ("MET_PT_Ele_Barrel_Num" , "MET_PT_Ele_Barrel_Num", 1000,0,1000,1000,0,1000);
   
   TH2F  *MET_PT_Ele_Endcap_Den  = new TH2F ("MET_PT_Ele_Endcap_Den" , "MET_PT_Ele_Endcap_Den", 1000,0,1000,1000,0,1000);
   TH2F  *MET_PT_Ele_Endcap_Num  = new TH2F ("MET_PT_Ele_Endcap_Num" , "MET_PT_Ele_Endcap_Num", 1000,0,1000,1000,0,1000);

   TH2F  *MET_PT_Mu_Barrel_Den  = new TH2F ("MET_PT_Mu_Barrel_Den" , "MET_PT_Mu_Barrel_Den", 1000,0,1000,1000,0,1000);
   TH2F  *MET_PT_Mu_Barrel_Num  = new TH2F ("MET_PT_Mu_Barrel_Num" , "MET_PT_Mu_Barrel_Num", 1000,0,1000,1000,0,1000);
   
   TH2F  *MET_PT_Mu_Endcap_Den  = new TH2F ("MET_PT_Mu_Endcap_Den" , "MET_PT_Mu_Endcap_Den", 1000,0,1000,1000,0,1000);
   TH2F  *MET_PT_Mu_Endcap_Num  = new TH2F ("MET_PT_Mu_Endcap_Num" , "MET_PT_Mu_Endcap_Num", 1000,0,1000,1000,0,1000);

   // //Rebined

   int nxbins =6;
   int nybins =9;
   double xbins[7]={0,20,40,60,80,100,200}; //MET
   double ybins[10]={5,7,10,20,30,40,50,80,100,150}; //PT
   //Muon
  

   TH2F  *MET_PT_Mu_Barrel_Den_Rebin  = new TH2F ("MET_PT_Mu_Barrel_Den_Rebin" , "MET_PT_Mu_Barrel_Den_Rebin", nxbins,xbins,nybins,ybins);
   TH2F  *MET_PT_Mu_Barrel_Num_Rebin  = new TH2F ("MET_PT_Mu_Barrel_Num_Rebin" , "MET_PT_Mu_Barrel_Num_Rebin", nxbins,xbins,nybins,ybins);
   
   TH2F  *MET_PT_Mu_Endcap_Den_Rebin  = new TH2F ("MET_PT_Mu_Endcap_Den_Rebin" , "MET_PT_Mu_Endcap_Den_Rebin", nxbins,xbins,nybins,ybins);
   TH2F  *MET_PT_Mu_Endcap_Num_Rebin  = new TH2F ("MET_PT_Mu_Endcap_Num_Rebin" , "MET_PT_Mu_Endcap_Num_Rebin", nxbins,xbins,nybins,ybins);

   TH2F  *h_MET_PT_loose_mu_barrel_pos = new TH2F ("h_MET_PT_loose_mu_barrel_pos" , "h_MET_PT_loose_mu_barrel_pos", nxbins,xbins,nybins,ybins);
   TH2F  *h_MET_PT_loose_mu_barrel_neg = new TH2F ("h_MET_PT_loose_mu_barrel_neg" , "h_MET_PT_loose_mu_barrel_neg", nxbins,xbins,nybins,ybins);
   
   TH2F  *h_MET_PT_loose_mu_endcap_pos = new TH2F ("h_MET_PT_loose_mu_endcap_pos" , "h_MET_PT_loose_mu_endcap_pos", nxbins,xbins,nybins,ybins);
   TH2F  *h_MET_PT_loose_mu_endcap_neg = new TH2F ("h_MET_PT_loose_mu_endcap_neg" , "h_MET_PT_loose_mu_endcap_neg", nxbins,xbins,nybins,ybins);

   TH2F  *h_MET_PT_tight_mu_barrel_pos = new TH2F ("h_MET_PT_tight_mu_barrel_pos" , "h_MET_PT_loose_mu_barrel_pos", nxbins,xbins,nybins,ybins);
   TH2F  *h_MET_PT_tight_mu_barrel_neg = new TH2F ("h_MET_PT_tight_mu_barrel_neg" , "h_MET_PT_loose_mu_barrel_neg", nxbins,xbins,nybins,ybins);
   
   TH2F  *h_MET_PT_tight_mu_endcap_pos = new TH2F ("h_MET_PT_tight_mu_endcap_pos" , "h_MET_PT_tight_mu_endcap_pos", nxbins,xbins,nybins,ybins);
   TH2F  *h_MET_PT_tight_mu_endcap_neg = new TH2F ("h_MET_PT_tight_mu_endcap_neg" , "h_MET_PT_tight_mu_endcap_neg", nxbins,xbins,nybins,ybins);

   //Electron

   TH2F  *MET_PT_Ele_Barrel_Den_Rebin  = new TH2F ("MET_PT_Ele_Barrel_Den_Rebin" , "MET_PT_Ele_Barrel_Den_Rebin", nxbins,xbins,nybins,ybins);
   TH2F  *MET_PT_Ele_Barrel_Num_Rebin  = new TH2F ("MET_PT_Ele_Barrel_Num_Rebin" , "MET_PT_Ele_Barrel_Num_Rebin", nxbins,xbins,nybins,ybins);
   
   TH2F  *MET_PT_Ele_Endcap_Den_Rebin  = new TH2F ("MET_PT_Ele_Endcap_Den_Rebin" , "MET_PT_Ele_Endcap_Den_Rebin", nxbins,xbins,nybins,ybins);
   TH2F  *MET_PT_Ele_Endcap_Num_Rebin  = new TH2F ("MET_PT_Ele_Endcap_Num_Rebin" , "MET_PT_Ele_Endcap_Num_Rebin", nxbins,xbins,nybins,ybins);

   TH2F  *h_MET_PT_loose_ele_barrel_pos = new TH2F ("h_MET_PT_loose_ele_barrel_pos" , "h_MET_PT_loose_ele_barrel_pos", nxbins,xbins,nybins,ybins);
   TH2F  *h_MET_PT_loose_ele_barrel_neg = new TH2F ("h_MET_PT_loose_ele_barrel_neg" , "h_MET_PT_loose_ele_barrel_neg", nxbins,xbins,nybins,ybins);
   
   TH2F  *h_MET_PT_loose_ele_endcap_pos = new TH2F ("h_MET_PT_loose_ele_endcap_pos" , "h_MET_PT_loose_ele_endcap_pos", nxbins,xbins,nybins,ybins);
   TH2F  *h_MET_PT_loose_ele_endcap_neg = new TH2F ("h_MET_PT_loose_ele_endcap_neg" , "h_MET_PT_loose_ele_endcap_neg", nxbins,xbins,nybins,ybins);

   TH2F  *h_MET_PT_tight_ele_barrel_pos = new TH2F ("h_MET_PT_tight_ele_barrel_pos" , "h_MET_PT_loose_ele_barrel_pos", nxbins,xbins,nybins,ybins);
   TH2F  *h_MET_PT_tight_ele_barrel_neg = new TH2F ("h_MET_PT_tight_ele_barrel_neg" , "h_MET_PT_loose_ele_barrel_neg", nxbins,xbins,nybins,ybins);
   
   TH2F  *h_MET_PT_tight_ele_endcap_pos = new TH2F ("h_MET_PT_tight_ele_endcap_pos" , "h_MET_PT_tight_ele_endcap_pos", nxbins,xbins,nybins,ybins);
   TH2F  *h_MET_PT_tight_ele_endcap_neg = new TH2F ("h_MET_PT_tight_ele_endcap_neg" , "h_MET_PT_tight_ele_endcap_neg", nxbins,xbins,nybins,ybins);
   
   
   //================//
   //Control regions//
   //===============//
   
   //MMMM
   
   TH1F * hM4l_MMMM_I_I = new TH1F("hM4l_MMMM_I_I", "hM4l_MMMM_I_I", 8000, 0., 2000. );
   TH1F * hMZ1_MMMM_I_I = new TH1F("hMZ1_MMMM_I_I", "hMZ1_MMMM_I_I", 800, 0., 200. );
   TH1F * hMZ2_MMMM_I_I = new TH1F("hMZ2_MMMM_I_I", "hMZ2_MMMM_I_I", 800, 0., 200. );
   
   
   TH1F * hM4l_MMMM_AI_AI = new TH1F("hM4l_MMMM_AI_AI", "hM4l_MMMM_AI_AI", 8000, 0., 2000. );
   TH1F * hM4l_MMMM_AI_I = new TH1F("hM4l_MMMM_AI_I", "hM4l_MMMM_AI_I", 8000, 0., 2000. );
   TH1F * hM4l_MMMM_I_AI = new TH1F("hM4l_MMMM_I_AI", "hM4l_MMMM_I_AI", 8000, 0., 2000. );
   TH1F * hM4l_MMMM_3p1f_total= new TH1F("hM4l_MMMM_3p1f_total", "hM4l_MMMM_3p1f_total", 8000, 0., 2000. );
   
   // for Simran

   TH1F * hPt_4mu_AI_AI_1 = new TH1F("hPt_4mu_AI_AI_1", "hPt_4mu_AI_AI_1", 200 , -0.5 , 199.5 );
   hPt_4mu_AI_AI_1->SetXTitle("pt_4mu  (GeV)");

   TH1F * hPt_4mu_AI_AI_2 = new TH1F("hPt_4mu_AI_AI_2", "hPt_4mu_AI_AI_2", 200 , 0.0 , 200.0 );
   hPt_4mu_AI_AI_2->SetXTitle("pt_4mu  (GeV)");

   TH1F * hPt_4mu_AI_AI_3 = new TH1F("hPt_4mu_AI_AI_3", "hPt_4mu_AI_AI_3", 1000 , 0.0 , 1000.0 );
   hPt_4mu_AI_AI_3->SetXTitle("pt_4mu  (GeV)");


   TH1F * hPFMET_MMMM_AI_AI = new TH1F("hPFMET_MMMM_AI_AI", "hPFMET_MMMM_AI_AI", 1000 , 0., 1000.);
   TH1F * hPFMET_MMMM_AI_I = new TH1F("hPFMET_MMMM_AI_I", "hPFMET_MMMM_AI_I", 1000 , 0., 1000.);
   TH1F * hPFMET_MMMM_I_AI = new TH1F("hPFMET_MMMM_I_AI", "hPFMET_MMMM_I_AI", 1000 , 0., 1000.);
   TH1F * hPFMET_MMMM_3p1f_total = new TH1F("hPFMET_MMMM_3p1f_total", "hPFMET_MMMM_3p1f_total", 1000 , 0., 1000.);
   
   //MMEE
   
   TH1F * hM4l_MMEE_I_I = new TH1F("hM4l_MMEE_I_I", "hM4l_MMEE_I_I", 8000, 0., 2000. );
   TH1F * hMZ1_MMEE_I_I = new TH1F("hMZ1_MMEE_I_I", "hMZ1_MMEE_I_I", 8000, 0., 2000. );
   TH1F * hMZ2_MMEE_I_I = new TH1F("hMZ2_MMEE_I_I", "hMZ2_MMEE_I_I", 8000, 0., 2000. );
   
   TH1F * hM4l_MMEE_AI_AI = new TH1F("hM4l_MMEE_AI_AI", "hM4l_MMEE_AI_AI", 8000, 0., 2000. );
   TH1F * hM4l_MMEE_AI_I = new TH1F("hM4l_MMEE_AI_I", "hM4l_MMEE_AI_I", 8000, 0., 2000. );
   TH1F * hM4l_MMEE_I_AI = new TH1F("hM4l_MMEE_I_AI", "hM4l_MMEE_I_AI", 8000, 0., 2000. );
   TH1F * hM4l_MMEE_3p1f_total= new TH1F("hM4l_MMEE_3p1f_total", "hM4l_MMEE_3p1f_total", 8000, 0., 2000. );

   TH1F * hPFMET_MMEE_AI_AI = new TH1F("hPFMET_MMEE_AI_AI", "hPFMET_MMEE_AI_AI", 1000 , 0., 1000.);
   TH1F * hPFMET_MMEE_AI_I = new TH1F("hPFMET_MMEE_AI_I", "hPFMET_MMEE_AI_I", 1000 , 0., 1000.);
   TH1F * hPFMET_MMEE_I_AI = new TH1F("hPFMET_MMEE_I_AI", "hPFMET_MMEE_I_AI", 1000 , 0., 1000.);
   TH1F * hPFMET_MMEE_3p1f_total = new TH1F("hPFMET_MMEE_3p1f_total", "hPFMET_MMEE_3p1f_total", 1000 , 0., 1000.);
   
   //EEMM
   
   TH1F * hM4l_EEMM_I_I = new TH1F("hM4l_EEMM_I_I", "hM4l_EEMM_I_I", 8000, 0., 2000. );
   TH1F * hMZ1_EEMM_I_I = new TH1F("hMZ1_EEMM_I_I", "hMZ1_EEMM_I_I", 8000, 0., 2000. );
   TH1F * hMZ2_EEMM_I_I = new TH1F("hMZ2_EEMM_I_I", "hMZ2_EEMM_I_I", 8000, 0., 2000. );
   
   TH1F * hM4l_EEMM_AI_AI = new TH1F("hM4l_EEMM_AI_AI", "hM4l_EEMM_AI_AI", 8000, 0., 2000. );
   TH1F * hM4l_EEMM_AI_I = new TH1F("hM4l_EEMM_AI_I", "hM4l_EEMM_AI_I", 8000, 0., 2000. );
   TH1F * hM4l_EEMM_I_AI = new TH1F("hM4l_EEMM_I_AI", "hM4l_EEMM_I_AI", 8000, 0., 2000. );
   TH1F * hM4l_EEMM_3p1f_total= new TH1F("hM4l_EEMM_3p1f_total", "hM4l_EEMM_3p1f_total", 8000, 0., 2000. );

   TH1F * hPFMET_EEMM_AI_AI = new TH1F("hPFMET_EEMM_AI_AI", "hPFMET_EEMM_AI_AI", 1000 , 0., 1000.);
   TH1F * hPFMET_EEMM_AI_I = new TH1F("hPFMET_EEMM_AI_I", "hPFMET_EEMM_AI_I", 1000 , 0., 1000.);
   TH1F * hPFMET_EEMM_I_AI = new TH1F("hPFMET_EEMM_I_AI", "hPFMET_EEMM_I_AI", 1000 , 0., 1000.);
   TH1F * hPFMET_EEMM_3p1f_total = new TH1F("hPFMET_EEMM_3p1f_total", "hPFMET_EEMM_3p1f_total", 1000 , 0., 1000.);
   
   //EEEE
   
   TH1F * hM4l_EEEE_I_I = new TH1F("hM4l_EEEE_I_I", "hM4l_EEEE_I_I", 8000, 0., 2000. );
   TH1F * hMZ1_EEEE_I_I = new TH1F("hMZ1_EEEE_I_I", "hMZ1_EEEE_I_I", 800, 0., 200. );
   TH1F * hMZ2_EEEE_I_I = new TH1F("hMZ2_EEEE_I_I", "hMZ2_EEEE_I_I", 800, 0., 200. );
   
   
   TH1F * hM4l_EEEE_AI_AI = new TH1F("hM4l_EEEE_AI_AI", "hM4l_EEEE_AI_AI", 8000, 0., 2000. );
   TH1F * hM4l_EEEE_AI_I = new TH1F("hM4l_EEEE_AI_I", "hM4l_EEEE_AI_I", 8000, 0., 2000. );
   TH1F * hM4l_EEEE_I_AI = new TH1F("hM4l_EEEE_I_AI", "hM4l_EEEE_I_AI", 8000, 0., 2000. );
   TH1F * hM4l_EEEE_3p1f_total= new TH1F("hM4l_EEEE_3p1f_total", "hM4l_EEEE_3p1f_total", 8000, 0., 2000. );

   TH1F * hPFMET_EEEE_AI_AI = new TH1F("hPFMET_EEEE_AI_AI", "hPFMET_EEEE_AI_AI", 1000 , 0., 1000.);
   TH1F * hPFMET_EEEE_AI_I = new TH1F("hPFMET_EEEE_AI_I", "hPFMET_EEEE_AI_I", 1000 , 0., 1000.);
   TH1F * hPFMET_EEEE_I_AI = new TH1F("hPFMET_EEEE_I_AI", "hPFMET_EEEE_I_AI", 1000 , 0., 1000.);
   TH1F * hPFMET_EEEE_3p1f_total = new TH1F("hPFMET_EEEE_3p1f_total", "hPFMET_EEEE_3p1f_total", 1000 , 0., 1000.);

   //around Z peak in 4l

   TH1F * hPFMET_MMMM_AI_AI_Zpeak = new TH1F("hPFMET_MMMM_AI_AI_Zpeak", "hPFMET_MMMM_AI_AI_Zpeak", 1000 , 0., 1000.);
   hPFMET_MMMM_AI_AI_Zpeak->SetXTitle("PF MET (GeV)");   
   TH1F * hPFMET_MMMM_AI_I_Zpeak = new TH1F("hPFMET_MMMM_AI_I_Zpeak", "hPFMET_MMMM_AI_I_Zpeak", 1000 , 0., 1000.);
   hPFMET_MMMM_AI_I_Zpeak->SetXTitle("PF MET (GeV)");
   TH1F * hPFMET_MMMM_I_AI_Zpeak = new TH1F("hPFMET_MMMM_I_AI_Zpeak", "hPFMET_MMMM_I_AI_Zpeak", 1000 , 0., 1000.);
   hPFMET_MMMM_I_AI_Zpeak->SetXTitle("PF MET (GeV)");
   TH1F * hPFMET_MMMM_3p1f_total_Zpeak = new TH1F("hPFMET_MMMM_3p1f_total_Zpeak", "hPFMET_MMMM_3p1f_total_Zpeak", 1000 , 0., 1000.);
   hPFMET_MMMM_3p1f_total_Zpeak->SetXTitle("PF MET (GeV)"); 
	    
   TH1F * hPFMET_EEMM_AI_AI_Zpeak = new TH1F("hPFMET_EEMM_AI_AI_Zpeak", "hPFMET_EEMM_AI_AI_Zpeak", 1000 , 0., 1000.);
   hPFMET_EEMM_AI_AI_Zpeak->SetXTitle("PF MET (GeV)");   
   TH1F * hPFMET_EEMM_AI_I_Zpeak = new TH1F("hPFMET_EEMM_AI_I_Zpeak", "hPFMET_EEMM_AI_I_Zpeak", 1000 , 0., 1000.);
   hPFMET_EEMM_AI_I_Zpeak->SetXTitle("PF MET (GeV)");
   TH1F * hPFMET_EEMM_I_AI_Zpeak = new TH1F("hPFMET_EEMM_I_AI_Zpeak", "hPFMET_EEMM_I_AI_Zpeak", 1000 , 0., 1000.);
   hPFMET_EEMM_I_AI_Zpeak->SetXTitle("PF MET (GeV)");
   TH1F * hPFMET_EEMM_3p1f_total_Zpeak = new TH1F("hPFMET_EEMM_3p1f_total_Zpeak", "hPFMET_EEMM_3p1f_total_Zpeak", 1000 , 0., 1000.);
   hPFMET_EEMM_3p1f_total_Zpeak->SetXTitle("PF MET (GeV)");


   TH1F * hPFMET_EEEE_AI_AI_Zpeak = new TH1F("hPFMET_EEEE_AI_AI_Zpeak", "hPFMET_EEEE_AI_AI_Zpeak", 1000 , 0., 1000.);
   hPFMET_EEEE_AI_AI_Zpeak->SetXTitle("PF MET (GeV)");   
   TH1F * hPFMET_EEEE_AI_I_Zpeak = new TH1F("hPFMET_EEEE_AI_I_Zpeak", "hPFMET_EEEE_AI_I_Zpeak", 1000 , 0., 1000.);
   hPFMET_EEEE_AI_I_Zpeak->SetXTitle("PF MET (GeV)");
   TH1F * hPFMET_EEEE_I_AI_Zpeak = new TH1F("hPFMET_EEEE_I_AI_Zpeak", "hPFMET_EEEE_I_AI_Zpeak", 1000 , 0., 1000.);
   hPFMET_EEEE_I_AI_Zpeak->SetXTitle("PF MET (GeV)");
   TH1F * hPFMET_EEEE_3p1f_total_Zpeak = new TH1F("hPFMET_EEEE_3p1f_total_Zpeak", "hPFMET_EEEE_3p1f_total_Zpeak", 1000 , 0., 1000.);
   hPFMET_EEEE_3p1f_total_Zpeak->SetXTitle("PF MET (GeV)");

   TH1F * hPFMET_MMEE_AI_AI_Zpeak = new TH1F("hPFMET_MMEE_AI_AI_Zpeak", "hPFMET_MMEE_AI_AI_Zpeak", 1000 , 0., 1000.);
   hPFMET_MMEE_AI_AI_Zpeak->SetXTitle("PF MET (GeV)");   
   TH1F * hPFMET_MMEE_AI_I_Zpeak = new TH1F("hPFMET_MMEE_AI_I_Zpeak", "hPFMET_MMEE_AI_I_Zpeak", 1000 , 0., 1000.);
   hPFMET_MMEE_AI_I_Zpeak->SetXTitle("PF MET (GeV)");
   TH1F * hPFMET_MMEE_I_AI_Zpeak = new TH1F("hPFMET_MMEE_I_AI_Zpeak", "hPFMET_MMEE_I_AI_Zpeak", 1000 , 0., 1000.);
   hPFMET_MMEE_I_AI_Zpeak->SetXTitle("PF MET (GeV)");
   TH1F * hPFMET_MMEE_3p1f_total_Zpeak = new TH1F("hPFMET_MMEE_3p1f_total_Zpeak", "hPFMET_MMEE_3p1f_total_Zpeak", 1000 , 0., 1000.);
   hPFMET_MMEE_3p1f_total_Zpeak->SetXTitle("PF MET (GeV)");
   
   /////

   const int NMOBINS = 5;
   const double MOMIN = 0.0, MOMAX = 1000.0;
   double loglinMbins[6]={0.,25.,50.,200., 500.,1000.};
   
   TH1F * hPFMET_MMMM_AI_AI_Zpeak_log = new TH1F("hPFMET_MMMM_AI_AI_Zpeak_log", "hPFMET_MMMM_AI_AI_Zpeak_log", NMOBINS ,loglinMbins );
   hPFMET_MMMM_AI_AI_Zpeak_log->Sumw2();   
   TH1F * hPFMET_MMMM_AI_I_Zpeak_log = new TH1F("hPFMET_MMMM_AI_I_Zpeak_log", "hPFMET_MMMM_AI_I_Zpeak_log", NMOBINS ,loglinMbins );
   hPFMET_MMMM_AI_I_Zpeak_log->Sumw2();   
   TH1F * hPFMET_MMMM_I_AI_Zpeak_log = new TH1F("hPFMET_MMMM_I_AI_Zpeak_log", "hPFMET_MMMM_I_AI_Zpeak_log", NMOBINS ,loglinMbins );
   hPFMET_MMMM_I_AI_Zpeak_log->Sumw2();
   TH1F * hPFMET_MMMM_3p1f_total_Zpeak_log = new TH1F("hPFMET_MMMM_3p1f_total_Zpeak_log", "hPFMET_MMMM_3p1f_total_Zpeak_log", NMOBINS ,loglinMbins );
   hPFMET_MMMM_3p1f_total_Zpeak_log->Sumw2();

   TH1F * hPFMET_EEMM_AI_AI_Zpeak_log = new TH1F("hPFMET_EEMM_AI_AI_Zpeak_log", "hPFMET_EEMM_AI_AI_Zpeak_log", NMOBINS ,loglinMbins );
   hPFMET_EEMM_AI_AI_Zpeak_log->Sumw2();   
   TH1F * hPFMET_EEMM_AI_I_Zpeak_log = new TH1F("hPFMET_EEMM_AI_I_Zpeak_log", "hPFMET_EEMM_AI_I_Zpeak_log", NMOBINS ,loglinMbins );
   hPFMET_EEMM_AI_I_Zpeak_log->Sumw2();   
   TH1F * hPFMET_EEMM_I_AI_Zpeak_log = new TH1F("hPFMET_EEMM_I_AI_Zpeak_log", "hPFMET_EEMM_I_AI_Zpeak_log", NMOBINS ,loglinMbins );
   hPFMET_EEMM_I_AI_Zpeak_log->Sumw2();
   TH1F * hPFMET_EEMM_3p1f_total_Zpeak_log = new TH1F("hPFMET_EEMM_3p1f_total_Zpeak_log", "hPFMET_EEMM_3p1f_total_Zpeak_log", NMOBINS ,loglinMbins );
   hPFMET_EEMM_3p1f_total_Zpeak_log->Sumw2();

   TH1F * hPFMET_EEEE_AI_AI_Zpeak_log = new TH1F("hPFMET_EEEE_AI_AI_Zpeak_log", "hPFMET_EEEE_AI_AI_Zpeak_log", NMOBINS ,loglinMbins );
   hPFMET_EEEE_AI_AI_Zpeak_log->Sumw2();   
   TH1F * hPFMET_EEEE_AI_I_Zpeak_log = new TH1F("hPFMET_EEEE_AI_I_Zpeak_log", "hPFMET_EEEE_AI_I_Zpeak_log", NMOBINS ,loglinMbins );
   hPFMET_EEEE_AI_I_Zpeak_log->Sumw2();   
   TH1F * hPFMET_EEEE_I_AI_Zpeak_log = new TH1F("hPFMET_EEEE_I_AI_Zpeak_log", "hPFMET_EEEE_I_AI_Zpeak_log", NMOBINS ,loglinMbins );
   hPFMET_EEEE_I_AI_Zpeak_log->Sumw2();
   TH1F * hPFMET_EEEE_3p1f_total_Zpeak_log = new TH1F("hPFMET_EEEE_3p1f_total_Zpeak_log", "hPFMET_EEEE_3p1f_total_Zpeak_log", NMOBINS ,loglinMbins );
   hPFMET_EEEE_3p1f_total_Zpeak_log->Sumw2();

   TH1F * hPFMET_MMEE_AI_AI_Zpeak_log = new TH1F("hPFMET_MMEE_AI_AI_Zpeak_log", "hPFMET_MMEE_AI_AI_Zpeak_log", NMOBINS ,loglinMbins );
   hPFMET_MMEE_AI_AI_Zpeak_log->Sumw2();   
   TH1F * hPFMET_MMEE_AI_I_Zpeak_log = new TH1F("hPFMET_MMEE_AI_I_Zpeak_log", "hPFMET_MMEE_AI_I_Zpeak_log", NMOBINS ,loglinMbins );
   hPFMET_MMEE_AI_I_Zpeak_log->Sumw2();   
   TH1F * hPFMET_MMEE_I_AI_Zpeak_log = new TH1F("hPFMET_MMEE_I_AI_Zpeak_log", "hPFMET_MMEE_I_AI_Zpeak_log", NMOBINS ,loglinMbins );
   hPFMET_MMEE_I_AI_Zpeak_log->Sumw2();
   TH1F * hPFMET_MMEE_3p1f_total_Zpeak_log = new TH1F("hPFMET_MMEE_3p1f_total_Zpeak_log", "hPFMET_MMEE_3p1f_total_Zpeak_log", NMOBINS ,loglinMbins );
   hPFMET_MMEE_3p1f_total_Zpeak_log->Sumw2();
   
   

   //=========================//
   //Control regions Step 2   //
   //=========================//
   
   //MMMM
   
   TH1F * hM4l_MMMM_I_I_Step2 = new TH1F("hM4l_MMMM_I_I_Step2", "hM4l_MMMM_I_I_Step2", 8000, 0., 2000. );
   TH1F * hMZ1_MMMM_I_I_Step2 = new TH1F("hMZ1_MMMM_I_I_Step2", "hMZ1_MMMM_I_I_Step2", 800, 0., 200. );
   TH1F * hMZ2_MMMM_I_I_Step2 = new TH1F("hMZ2_MMMM_I_I_Step2", "hMZ2_MMMM_I_I_Step2", 800, 0., 200. );
   
   
   TH1F * hM4l_MMMM_AI_AI_Step2 = new TH1F("hM4l_MMMM_AI_AI_Step2", "hM4l_MMMM_AI_AI_Step2", 8000, 0., 2000. );
   TH1F * hM4l_MMMM_AI_I_Step2 = new TH1F("hM4l_MMMM_AI_I_Step2", "hM4l_MMMM_AI_I_Step2", 8000, 0., 2000. );
   TH1F * hM4l_MMMM_I_AI_Step2 = new TH1F("hM4l_MMMM_I_AI_Step2", "hM4l_MMMM_I_AI_Step2", 8000, 0., 2000. );
   TH1F * hM4l_MMMM_3p1f_total_Step2= new TH1F("hM4l_MMMM_3p1f_total_Step2", "hM4l_MMMM_3p1f_total_Step2", 8000, 0., 2000. );

   
   TH1F * hPFMET_MMMM_AI_AI_Step2 = new TH1F("hPFMET_MMMM_AI_AI_Step2", "hPFMET_MMMM_AI_AI_Step2", 1000 , 0., 1000.);
   TH1F * hPFMET_MMMM_AI_I_Step2 = new TH1F("hPFMET_MMMM_AI_I_Step2", "hPFMET_MMMM_AI_I_Step2", 1000 , 0., 1000.);
   TH1F * hPFMET_MMMM_I_AI_Step2 = new TH1F("hPFMET_MMMM_I_AI_Step2", "hPFMET_MMMM_I_AI_Step2", 1000 , 0., 1000.);
   TH1F * hPFMET_MMMM_3p1f_total_Step2 = new TH1F("hPFMET_MMMM_3p1f_total_Step2", "hPFMET_MMMM_3p1f_total_Step2", 1000 , 0., 1000.);
   
   //MMEE
   
   TH1F * hM4l_MMEE_I_I_Step2 = new TH1F("hM4l_MMEE_I_I_Step2", "hM4l_MMEE_I_I_Step2", 8000, 0., 2000. );
   TH1F * hMZ1_MMEE_I_I_Step2 = new TH1F("hMZ1_MMEE_I_I_Step2", "hMZ1_MMEE_I_I_Step2", 8000, 0., 2000. );
   TH1F * hMZ2_MMEE_I_I_Step2 = new TH1F("hMZ2_MMEE_I_I_Step2", "hMZ2_MMEE_I_I_Step2", 8000, 0., 2000. );
   
   TH1F * hM4l_MMEE_AI_AI_Step2 = new TH1F("hM4l_MMEE_AI_AI_Step2", "hM4l_MMEE_AI_AI_Step2", 8000, 0., 2000. );
   TH1F * hM4l_MMEE_AI_I_Step2 = new TH1F("hM4l_MMEE_AI_I_Step2", "hM4l_MMEE_AI_I_Step2", 8000, 0., 2000. );
   TH1F * hM4l_MMEE_I_AI_Step2 = new TH1F("hM4l_MMEE_I_AI_Step2", "hM4l_MMEE_I_AI_Step2", 8000, 0., 2000. );
   TH1F * hM4l_MMEE_3p1f_total_Step2= new TH1F("hM4l_MMEE_3p1f_total_Step2", "hM4l_MMEE_3p1f_total_Step2", 8000, 0., 2000. );

   TH1F * hPFMET_MMEE_AI_AI_Step2 = new TH1F("hPFMET_MMEE_AI_AI_Step2", "hPFMET_MMEE_AI_AI_Step2", 1000 , 0., 1000.);
   TH1F * hPFMET_MMEE_AI_I_Step2 = new TH1F("hPFMET_MMEE_AI_I_Step2", "hPFMET_MMEE_AI_I_Step2", 1000 , 0., 1000.);
   TH1F * hPFMET_MMEE_I_AI_Step2 = new TH1F("hPFMET_MMEE_I_AI_Step2", "hPFMET_MMEE_I_AI_Step2", 1000 , 0., 1000.);
   TH1F * hPFMET_MMEE_3p1f_total_Step2 = new TH1F("hPFMET_MMEE_3p1f_total_Step2", "hPFMET_MMEE_3p1f_total_Step2", 1000 , 0., 1000.);
   
   //EEMM
   
   TH1F * hM4l_EEMM_I_I_Step2 = new TH1F("hM4l_EEMM_I_I_Step2", "hM4l_EEMM_I_I_Step2", 8000, 0., 2000. );
   TH1F * hMZ1_EEMM_I_I_Step2 = new TH1F("hMZ1_EEMM_I_I_Step2", "hMZ1_EEMM_I_I_Step2", 8000, 0., 2000. );
   TH1F * hMZ2_EEMM_I_I_Step2 = new TH1F("hMZ2_EEMM_I_I_Step2", "hMZ2_EEMM_I_I_Step2", 8000, 0., 2000. );
   
   TH1F * hM4l_EEMM_AI_AI_Step2 = new TH1F("hM4l_EEMM_AI_AI_Step2", "hM4l_EEMM_AI_AI_Step2", 8000, 0., 2000. );
   TH1F * hM4l_EEMM_AI_I_Step2 = new TH1F("hM4l_EEMM_AI_I_Step2", "hM4l_EEMM_AI_I_Step2", 8000, 0., 2000. );
   TH1F * hM4l_EEMM_I_AI_Step2 = new TH1F("hM4l_EEMM_I_AI_Step2", "hM4l_EEMM_I_AI_Step2", 8000, 0., 2000. );
   TH1F * hM4l_EEMM_3p1f_total_Step2= new TH1F("hM4l_EEMM_3p1f_total_Step2", "hM4l_EEMM_3p1f_total_Step2", 8000, 0., 2000. );

   TH1F * hPFMET_EEMM_AI_AI_Step2 = new TH1F("hPFMET_EEMM_AI_AI_Step2", "hPFMET_EEMM_AI_AI_Step2", 1000 , 0., 1000.);
   TH1F * hPFMET_EEMM_AI_I_Step2 = new TH1F("hPFMET_EEMM_AI_I_Step2", "hPFMET_EEMM_AI_I_Step2", 1000 , 0., 1000.);
   TH1F * hPFMET_EEMM_I_AI_Step2 = new TH1F("hPFMET_EEMM_I_AI_Step2", "hPFMET_EEMM_I_AI_Step2", 1000 , 0., 1000.);
   TH1F * hPFMET_EEMM_3p1f_total_Step2 = new TH1F("hPFMET_EEMM_3p1f_total_Step2", "hPFMET_EEMM_3p1f_total_Step2", 1000 , 0., 1000.);
   
   //EEEE
   
   TH1F * hM4l_EEEE_I_I_Step2 = new TH1F("hM4l_EEEE_I_I_Step2", "hM4l_EEEE_I_I_Step2", 8000, 0., 2000. );
   TH1F * hMZ1_EEEE_I_I_Step2 = new TH1F("hMZ1_EEEE_I_I_Step2", "hMZ1_EEEE_I_I_Step2", 800, 0., 200. );
   TH1F * hMZ2_EEEE_I_I_Step2 = new TH1F("hMZ2_EEEE_I_I_Step2", "hMZ2_EEEE_I_I_Step2", 800, 0., 200. );
   
   
   TH1F * hM4l_EEEE_AI_AI_Step2 = new TH1F("hM4l_EEEE_AI_AI_Step2", "hM4l_EEEE_AI_AI_Step2", 8000, 0., 2000. );
   TH1F * hM4l_EEEE_AI_I_Step2 = new TH1F("hM4l_EEEE_AI_I_Step2", "hM4l_EEEE_AI_I_Step2", 8000, 0., 2000. );
   TH1F * hM4l_EEEE_I_AI_Step2 = new TH1F("hM4l_EEEE_I_AI_Step2", "hM4l_EEEE_I_AI_Step2", 8000, 0., 2000. );
   TH1F * hM4l_EEEE_3p1f_total_Step2= new TH1F("hM4l_EEEE_3p1f_total_Step2", "hM4l_EEEE_3p1f_total_Step2", 8000, 0., 2000. );

   TH1F * hPFMET_EEEE_AI_AI_Step2 = new TH1F("hPFMET_EEEE_AI_AI_Step2", "hPFMET_EEEE_AI_AI_Step2", 1000 , 0., 1000.);
   TH1F * hPFMET_EEEE_AI_I_Step2 = new TH1F("hPFMET_EEEE_AI_I_Step2", "hPFMET_EEEE_AI_I_Step2", 1000 , 0., 1000.);
   TH1F * hPFMET_EEEE_I_AI_Step2 = new TH1F("hPFMET_EEEE_I_AI_Step2", "hPFMET_EEEE_I_AI_Step2", 1000 , 0., 1000.);
   TH1F * hPFMET_EEEE_3p1f_total_Step2 = new TH1F("hPFMET_EEEE_3p1f_total_Step2", "hPFMET_EEEE_3p1f_total_Step2", 1000 , 0., 1000.);

   //Revert Higgs mass window in the control regions to get PFMET systematics

   //In MonoHiggs step

   //MMMM

   TH1F * hPFMET_MMMM_AI_AI_Step2_rev = new TH1F("hPFMET_MMMM_AI_AI_Step2_rev", "hPFMET_MMMM_AI_AI_Step2_rev", 1000 , 0., 1000.);
   TH1F * hPFMET_MMMM_AI_I_Step2_rev = new TH1F("hPFMET_MMMM_AI_I_Step2_rev", "hPFMET_MMMM_AI_I_Step2_rev", 1000 , 0., 1000.);
   TH1F * hPFMET_MMMM_I_AI_Step2_rev = new TH1F("hPFMET_MMMM_I_AI_Step2_rev", "hPFMET_MMMM_I_AI_Step2_rev", 1000 , 0., 1000.);
   TH1F * hPFMET_MMMM_3p1f_total_Step2_rev = new TH1F("hPFMET_MMMM_3p1f_total_Step2_rev", "hPFMET_MMMM_3p1f_total_Step2_rev", 1000 , 0., 1000.);

   //MMEE

   TH1F * hPFMET_MMEE_AI_AI_Step2_rev = new TH1F("hPFMET_MMEE_AI_AI_Step2_rev", "hPFMET_MMEE_AI_AI_Step2_rev", 1000 , 0., 1000.);
   TH1F * hPFMET_MMEE_AI_I_Step2_rev = new TH1F("hPFMET_MMEE_AI_I_Step2_rev", "hPFMET_MMEE_AI_I_Step2_rev", 1000 , 0., 1000.);
   TH1F * hPFMET_MMEE_I_AI_Step2_rev = new TH1F("hPFMET_MMEE_I_AI_Step2_rev", "hPFMET_MMEE_I_AI_Step2_rev", 1000 , 0., 1000.);
   TH1F * hPFMET_MMEE_3p1f_total_Step2_rev = new TH1F("hPFMET_MMEE_3p1f_total_Step2_rev", "hPFMET_MMEE_3p1f_total_Step2_rev", 1000 , 0., 1000.);

   //EEMM

   TH1F * hPFMET_EEMM_AI_AI_Step2_rev = new TH1F("hPFMET_EEMM_AI_AI_Step2_rev", "hPFMET_EEMM_AI_AI_Step2_rev", 1000 , 0., 1000.);
   TH1F * hPFMET_EEMM_AI_I_Step2_rev = new TH1F("hPFMET_EEMM_AI_I_Step2_rev", "hPFMET_EEMM_AI_I_Step2_rev", 1000 , 0., 1000.);
   TH1F * hPFMET_EEMM_I_AI_Step2_rev = new TH1F("hPFMET_EEMM_I_AI_Step2_rev", "hPFMET_EEMM_I_AI_Step2_rev", 1000 , 0., 1000.);
   TH1F * hPFMET_EEMM_3p1f_total_Step2_rev = new TH1F("hPFMET_EEMM_3p1f_total_Step2_rev", "hPFMET_EEMM_3p1f_total_Step2_rev", 1000 , 0., 1000.);

   //EEEE

   TH1F * hPFMET_EEEE_AI_AI_Step2_rev = new TH1F("hPFMET_EEEE_AI_AI_Step2_rev", "hPFMET_EEEE_AI_AI_Step2_rev", 1000 , 0., 1000.);
   TH1F * hPFMET_EEEE_AI_I_Step2_rev = new TH1F("hPFMET_EEEE_AI_I_Step2_rev", "hPFMET_EEEE_AI_I_Step2_rev", 1000 , 0., 1000.);
   TH1F * hPFMET_EEEE_I_AI_Step2_rev = new TH1F("hPFMET_EEEE_I_AI_Step2_rev", "hPFMET_EEEE_I_AI_Step2_rev", 1000 , 0., 1000.);
   TH1F * hPFMET_EEEE_3p1f_total_Step2_rev = new TH1F("hPFMET_EEEE_3p1f_total_Step2_rev", "hPFMET_EEEE_3p1f_total_Step2_rev", 1000 , 0., 1000.);
   
   
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
   //this number nentries shows the number of input events from input root file
   //cout << "\n****************************"  <<endl;
   //cout << "Analyzing " << nentries << " entries"  <<endl;   
   
   for (Long64_t jentry=0; jentry<nentries;jentry++){
     
     //Reset variables for the ttree
     f_weight=-999;
     f_run=-999, f_lumi=-999, f_event=-999;
     f_zee_mu_pt_loose=-999, f_zee_e_pt_loose=-999, f_zmumu_mu_pt_loose=-999, f_zmumu_e_pt_loose=-999;
     f_zee_mu_eta_loose=-999, f_zee_e_eta_loose=-999, f_zmumu_mu_eta_loose=-999, f_zmumu_e_eta_loose=-999;
     f_zee_mu_charge_loose=-999, f_zee_e_charge_loose=-999, f_zmumu_mu_charge_loose=-999, f_zmumu_e_charge_loose=-999;
     f_zee_mu_pt_tight=-999, f_zee_e_pt_tight=-999, f_zmumu_mu_pt_tight=-999, f_zmumu_e_pt_tight=-999;
     f_zee_mu_eta_tight=-999, f_zee_e_eta_tight=-999, f_zmumu_mu_eta_tight=-999, f_zmumu_e_eta_tight=-999;
     f_zee_mu_charge_tight=-999, f_zee_e_charge_tight=-999, f_zmumu_mu_charge_tight=-999, f_zmumu_e_charge_tight=-999;
     
     f_weight2=-999, f_Djet_VAJHU=-999;
     f_run2=-999, f_lumi2=-999, f_event2=-999, f_2p2f=0, f_3p1f=0, f_2p2p=0;
     f_lept1_pdgid=-999, f_lept2_pdgid=-999, f_lept3_pdgid=-999, f_lept4_pdgid=-999;
     f_lept1_pass=0, f_lept2_pass=0, f_lept3_pass=0, f_lept4_pass=0;
     f_lept1_pt=-999, f_lept1_pt_error=-999, f_lept1_eta=-999, f_lept1_phi=-999;
     f_lept2_pt=-999, f_lept2_pt_error=-999, f_lept2_eta=-999, f_lept2_phi=-999;
     f_lept3_pt=-999, f_lept3_pt_error=-999, f_lept3_eta=-999, f_lept3_phi=-999;
     f_lept4_pt=-999, f_lept4_pt_error=-999, f_lept4_eta=-999, f_lept4_phi=-999;
     f_pfmet=-999, f_pfmet_JetEnUp=-999, f_pfmet_JetEnDn=-999, f_pfmet_ElectronEnUp=-999, f_pfmet_ElectronEnDn=-999, f_pfmet_MuonEnUp=-999, f_pfmet_MuonEnDn=-999, f_pfmet_JetResUp=-999, f_pfmet_JetResDn=-999, f_pfmet_UnclusteredEnUp=-999, f_pfmet_UnclusteredEnDn=-999, f_pfmet_PhotonEnUp=-999, f_pfmet_PhotonEnDn=-999;    
     f_njets_pass=-999, f_Nbjets=-999;
     
     for(int ijet=0; ijet<3; ++ijet){
      f_jets_highpt_pt[ijet] = -999;
      f_jets_highpt_pt_error[ijet] = -999;
      f_jets_highpt_eta[ijet] = -999;
      f_jets_highpt_phi[ijet] = -999;
      f_jets_highpt_et[ijet] = -999;
     }
     
     
     fChain->GetEntry(jentry);
     
     if(Nevents != -1 && jentry > Nevents) break;
     //cout<< "Processing event " << jentry <<endl;
     
     if( RECO_NMU > 100 ) RECO_NMU = 100;
     if( RECO_NELE > 100 ) RECO_NELE = 100;
     if( RECO_NPFPHOT > 30 ) RECO_NPFPHOT = 30;
     
     event++;
     histo_EVENT->Fill(event);
     
     
     //Weights
     //========
     double weight = 1.0;
     if(MC_type == "Spring16"){
        float lumifb = 35.867;
	weight = lumifb*(Xsection*1000.*neventsPostHLT/neventsPreHLT)/neventsPostHLT;
     }
     
     newweight_3=weight;
     newweight=weight;
     
     cout << "Starting weight= " << newweight << endl;
     
     // pileup reweighting 2016
     //-----------------------------------
     
     //cout<<"DATA_type= "<<DATA_type<<" and MC_type = "<<MC_type<<" and num_PU_vertices = "<<num_PU_vertices<<endl;
     
     if (DATA_type=="NO" && num_PU_vertices < 0) continue;  //make sure
     hPUvertices->Fill(num_PU_vertices,weight);
     
     double pu_weight=1.;
     
     //cout<<"DATA_type= "<<DATA_type<<endl;
     
     if (MC_type == "Spring16"){
       Int_t binx = puweight->GetXaxis()->FindBin(num_PU_vertices);
       cout << " bin x= " << binx << " " << puweight->GetBinContent(binx) << endl;	
       pu_weight=double(puweight->GetBinContent(binx));	
     }

     hPUvertices_ReWeighted->Fill(num_PU_vertices,weight*pu_weight);
     cout << "Pileup interations and weight is= " << num_PU_vertices << " " << " and weight= " << pu_weight << endl;
     
     // Changing the weight for pileup
     
     newweight_3=weight*pu_weight;
     newweight=weight*pu_weight;
     cout << "Starting weight + pileup = " << newweight << endl;
     
     ///////////////////////////////////////////////////////
     
     // ggZZ kfactor
     
     double ggzz_kf_wgt[9];
     float  weight_kfactor=1.;
     
     if (DATA_type=="NO"){
       // if( datasetName.Contains("GluGluHToZZ")){
       // for(int f=0;f<9;f++) ggzz_kf_wgt[f] = ggZZ_kf[f]->Eval(MC_MASS[0]); // Evaluate at the true m4l
       // weight_kfactor=ggzz_kf_wgt[0]; // Using the nominal one
       
       // newweight_3=weight*pu_weight*weight_kfactor;
       // newweight=weight*pu_weight*weight_kfactor;
       //}
       if ( datasetName.Contains("GluGluToZZ") ||  datasetName.Contains("GluGluToContinToZZ")){
	 //cout<<"K factor step GluGluToZZ"<<endl;
	 for(int f=0;f<9;f++) ggzz_kf_wgt[f] = ggZZ_kf[f]->Eval(MC_ZZ_MASS[0][0]); // Evaluate at the true m4l
	 weight_kfactor=ggzz_kf_wgt[0]; // Using the nominal one
	 //cout<<"weight_kfactor = "<<weight_kfactor<<endl;
	 
	 newweight_3=weight*pu_weight*weight_kfactor;
	 newweight=weight*pu_weight*weight_kfactor;
	 //cout<<"newweight = "<<newweight<<endl;
       }
       //weight_kfactor=2.3;
     }
     
     
     // qqZZ kfactor
     
     double qqzz_kf_wgt;
     weight_kfactor=1.;
     int finalState=-999;
     if (DATA_type=="NO"){
       if( datasetName.Contains("ZZTo4L") )  {
	 //cout<<"K factor step qqToZZ"<<endl;
	 for (int l=0;l<4;l++){
	   if (MC_ZZ_MASS[l][0]>0. &&
     		fabs(MC_ZZ_PDGID[l][3])==fabs(MC_ZZ_PDGID[l][4]) && 
	       fabs(MC_ZZ_PDGID[l][3])==fabs(MC_ZZ_PDGID[l][5]) &&
	       fabs(MC_ZZ_PDGID[l][3])==fabs(MC_ZZ_PDGID[l][6])) finalState=1; // 4e, 4mu, 4tau
	   else finalState=2;
	   weight_kfactor=kfactor_qqZZ_qcd_M(MC_ZZ_MASS[l][0],finalState);
	   //cout<<"weight_kfactor = "<<weight_kfactor<<endl;
	   newweight_3=weight*pu_weight*weight_kfactor;
	   newweight=weight*pu_weight*weight_kfactor;
	   //cout<<"newweight = "<<newweight<<endl;
	 }	
       }
     }
     
     //       // Weight for MCNLO samples
     
     // if( datasetName.Contains("amcatnlo")) {
     //   cout << "Reweighting sample of amcatnlo with weight= " << MC_weighting << endl;
     //   newweight=weight*pu_weight*weight_kfactor*MC_weighting;
     // }
     
     //cout<<"************************************************************************"<<endl;
     
     
     // ** Step 0
     //============
     
     //simply number of entries...
     //cout << "\n** Step 0: \nAnalyzing entry: " << jentry << " Run: " << Run << " Event: " << Event << " LumiSection: " << LumiSection << endl ; 
     N_0++; //number of events before any cut
     
     
     // Loose leptons identification 
     //==============================
     
     // 1) loose muon identification
     //--------------------------------
     
     int N_loose_mu = 0;
     double max_Iso_loose_mu = -1;
     double max_Sip_loose_mu = -1;
     double max_Ip_loose_mu = -1;
     
     //To define variable size array
     
     int* arraysize_mu = new int[1];
     arraysize_mu[0] = RECO_NMU;
     int iL_loose_mu[arraysize_mu[0]];
     delete [] arraysize_mu;
     
     for( int i = 0; i < RECO_NMU; i++ ){
       iL_loose_mu[i]=-999.;
     }
     
     
     for(int i = 0; i < RECO_NMU; i++){
       /*
       // if (debug)
       cout << "\n muon i="<< i <<" properties: "
	    << "\nRECOMU_isPFMu[i] " << int(RECOMU_isPFMu[i])
	    << "\nRECOMU_isGlobalMu[i] " << int(RECOMU_isGlobalMu[i])
	    << "\nRECOMU_isTrackerMu[i] " << int(RECOMU_isTrackerMu[i])
	    << "\nRECOMU_PT[i] " << RECOMU_PT[i]
	    << "\nfabs(RECOMU_ETA[i]) " << fabs(RECOMU_ETA[i])
	    << "\nfabs(RECOMU_PHI[i]) " << RECOMU_PHI[i]
	    << "\nfabs( RECOMU_muInnertrkDxy[i] ) " << fabs( RECOMU_muInnertrkDxy[i] )
	    << "\nfabs( RECOMU_muInnertrkDz[i] ) " << fabs( RECOMU_muInnertrkDz[i] )
	    <<"\nfabs( RECOMU_PFX_dB[i] ) " <<RECOMU_PFX_dB[i]
	    <<"\nfabs( RECOMU_SIP[i] ) " <<RECOMU_SIP[i]
	    <<"\nfabs( RECOMU_IP[i] ) " <<RECOMU_IP[i]
	    <<"\nRECOMU_numberOfMatches[i] "<<RECOMU_numberOfMatches[i]
	    <<"\nRECOMU_mubesttrkType[i]" <<RECOMU_mubesttrkType[i]
	    << endl ;
       */
       
       if( ( RECOMU_isGlobalMu[i] || (RECOMU_isTrackerMu[i] && RECOMU_numberOfMatches[i]>0)) //new part in 13 TeV  RECOMU_numberOfMatches[i]>0
	   && RECOMU_mubesttrkType[i]!=2 //new part for 13 TeV
	   && RECOMU_PT[i] > 5. 
	   && fabs(RECOMU_ETA[i]) < 2.4 
	   // && fabs(RECOMU_muInnertrkDxy[i]) < .5 && fabs(RECOMU_muInnertrkDz[i]) < 1.
	   && fabs(RECOMU_mubesttrkDxy[i]) < .5 && fabs(RECOMU_mubesttrkDz[i]) < 1. // use muon best track instead of muon inner track in 13 TeV
	   //&& fabs( RECOMU_SIP[i] ) < 4.
	   )
	 { 
	   iL_loose_mu[N_loose_mu]=i; 
	   N_loose_mu++ ;
	   if( RECOMU_PFX_dB[i] > max_Iso_loose_mu ) max_Iso_loose_mu = RECOMU_PFX_dB[i] ;
	   if( fabs( RECOMU_SIP[i] ) > max_Sip_loose_mu ) max_Sip_loose_mu = fabs( RECOMU_SIP[i] ) ;
	   if( fabs( RECOMU_IP[i] ) > max_Ip_loose_mu ) max_Ip_loose_mu = fabs( RECOMU_IP[i] ) ;
	 }
     } //end loop on loose muon 
     
     
     hN_loose_mu->Fill( N_loose_mu,newweight );
     hIso_loose_mu->Fill( max_Iso_loose_mu,newweight);
     hSip_loose_mu->Fill( max_Sip_loose_mu,newweight );
     hIp_loose_mu->Fill( max_Ip_loose_mu,newweight );
     
     // cout<<"=================================="<<endl;

     //cout<<"losse muons in the event = "<<N_loose_mu<<"with indices"<< iL_loose_mu[0] <<","<< iL_loose_mu[1] <<","<< iL_loose_mu[2] <<","<< iL_loose_mu[3] <<","<< iL_loose_mu[4] <<","<< iL_loose_mu[5] <<endl;
     
     // 2) loose electron identification
     //----------------------------------
     
     // Effective AREA
     bool tag_2011=false;
     // if (DATA_type=="2010" || DATA_type=="2011" || MC_type=="Fall11"){
     //tag_2011=true;
     //  }
     
     int N_loose_e =0;
     double max_Iso_loose_e = -1;
     double max_Sip_loose_e = -1;
     double max_Ip_loose_e = -1;
     
     int* arraysize_e = new int[1];
     arraysize_e[0] = RECO_NELE;
     int iL_loose_e[arraysize_e[0]];
     delete [] arraysize_e;
     
     for( int i = 0; i < RECO_NELE; i++ ){
       iL_loose_e[i]=-999.;
     }
     
     for( int i = 0; i < RECO_NELE; ++i ){
       /*
       // if (debug)
       cout << "\n electron i="<< i <<" properties: "
	    << "\nRECOELE_PT[i] " << RECOELE_PT[i]
	    << "\nfabs(RECOELE_ETA[i]) " << fabs(RECOELE_ETA[i])
	    << "\nfabs(RECOELE_PHI[i]) " <<RECOELE_PHI[i]
	    << "\nfabs( RECOELE_gsftrack_dxy[i] ) " << fabs( RECOELE_gsftrack_dxy[i] )
	    << "\nfabs( RECOELE_gsftrack_dz[i] ) " << fabs( RECOELE_gsftrack_dz[i] )
	    <<"\nfabs(  RECOELE_PFX_rho[i] ) " << fabs(  RECOELE_PFX_rho[i] )
	    <<"\nfabs( RECOELE_SIP[i] ) " << fabs( RECOELE_SIP[i] )
	    <<"\nfabs( RECOELE_IP[i] ) " << fabs( RECOELE_IP[i] )
	    <<"\nRECOELE_gsftrack_expected_inner_hits[i] "<< RECOELE_gsftrack_expected_inner_hits[i]
	    << endl ;
       */
       if( RECOELE_PT[i] > 7. 
	   && fabs(RECOELE_ETA[i]) < 2.5 
	   //	&& RECOELE_gsftrack_expected_inner_hits[i]<=1  not used in 13 TeV
	   && fabs(RECOELE_gsftrack_dxy[i]) < .5  && fabs(RECOELE_gsftrack_dz[i]) < 1. 
	   ) { 
	 iL_loose_e[N_loose_e]=i;	
	 N_loose_e++ ;
	 if( RECOELE_PFX_rho[i] > max_Iso_loose_e ) max_Iso_loose_e = RECOELE_PFX_rho[i] ;
	 if( fabs( RECOELE_SIP[i] ) > max_Sip_loose_e ) max_Sip_loose_e = fabs( RECOELE_SIP[i] ) ;
	 if( fabs( RECOELE_IP[i] ) > max_Ip_loose_e ) max_Ip_loose_e = fabs( RECOELE_IP[i] ) ;
       }
     }// end loop on electrons 
     
     hN_loose_e->Fill( N_loose_e,newweight );
     hIso_loose_e->Fill( max_Iso_loose_e,newweight);
     hSip_loose_e->Fill( max_Sip_loose_e,newweight );
     hIp_loose_e->Fill( max_Ip_loose_e,newweight );
     
     //cout<<"losse electrons in the event = "<<N_loose_e<<"with indices"<< iL_loose_e[0] <<","<< iL_loose_e[1] <<","<< iL_loose_e[2] <<","<< iL_loose_e[3] <<","<< iL_loose_e[4] <<","<< iL_loose_e[5] <<endl; 
     
     // Electron Cleaning  -- electrons separated from muons (deltaR > 0.05) 
     //=====================================================================
     //this cleaning is done for all muons and electrons in the event ------ no isolation cuts on muons or electrons
     // each electron will loop in all good muons in the event for emu cleaning -- electrons must be separated from muons with deltaR > 0.05
     // veto or delet any electron with muon in cone with deltaR <0.05 (deltaR mustbe > 0.05) 
     //cleaning done for muons passing tight (good) id and SIP<4 cut
 
     for(int e = 0; e < RECO_NELE; e++){
       for(int mu = 0; mu < RECO_NMU; mu++){
	 
	 if( (RECOMU_isPFMu[mu] || (RECOMU_isTrackerHighPtMu[mu] && RECOMU_PT[mu] > 200.)) //new than loose muons
	     && (RECOMU_isGlobalMu[mu] || (RECOMU_isTrackerMu[mu] && RECOMU_numberOfMatches[mu]>0))
	     && RECOMU_mubesttrkType[mu]!=2
	     && RECOMU_PT[mu] > 5. 
	     && fabs(RECOMU_ETA[mu]) < 2.4 
	     && fabs(RECOMU_mubesttrkDxy[mu]) < .5 && fabs(RECOMU_mubesttrkDz[mu]) < 1. 
	     && fabs(RECOMU_SIP[mu])<4.  //new than loose muons 
	     );
	 else continue;
	 
	 double deltaR = sqrt( pow( DELTAPHI( RECOMU_PHI[mu] , RECOELE_PHI[e] ),2) + pow(RECOMU_ETA[mu] - RECOELE_ETA[e],2) );
	 //if(debug) cout<<"deltaR = "<<deltaR<<endl;
	 if( deltaR <= 0.05 ){
	   
	   //if( debug ) cout << "Event not passing the HLT trigger paths" << endl;
	   
	   RECOELE_PT[e]  = -0.01;
	   RECOELE_ETA[e] = -99.;
	   RECOELE_PHI[e] = -99.;
	   RECOELE_SIP[e] = -99.;
	 }
       }
     }
     // end of electron muon cleaning 
     
     //3) good muon identification no isolation
     
     int iL[8]= { -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1};
     
     int N_good = 0 ;
     
     for(int i = 0; i < RECO_NMU; i++){
       /*
       if (debug) cout << "\n muon i="<< i <<" properties: "
		       << "\nRECOMU_isPFMu[i] " << int(RECOMU_isPFMu[i])
		       << "\nRECOMU_isGlobalMu[i] " << int(RECOMU_isGlobalMu[i])
		       << "\nRECOMU_isTrackerMu[i] " << int(RECOMU_isTrackerMu[i])
		       << "\nRECOMU_trkmuArbitration[i] " << int(RECOMU_trkmuArbitration[i])
		       << "\nRECOMU_PT[i] " << RECOMU_PT[i]
		       << "\nfabs(RECOMU_ETA[i]) " << fabs(RECOMU_ETA[i])
		       << "\nfabs(RECOMU_PHI[i]) " << RECOMU_PHI[i]
		       << "\nfabs( RECOMU_muInnertrkDxy[i] ) " << fabs( RECOMU_muInnertrkDxy[i] )
		       << "\nfabs( RECOMU_muInnertrkDz[i] ) " << fabs( RECOMU_muInnertrkDz[i] )
		       <<"\nfabs( RECOMU_PFX_dB[i] ) " <<RECOMU_PFX_dB[i]
		       <<"\nfabs( RECOMU_SIP[i] ) " <<RECOMU_SIP[i]
		       <<"\nfabs( RECOMU_IP[i] ) " <<RECOMU_IP[i]
		       <<"\nRECOMU_numberOfMatches[i] "<<RECOMU_numberOfMatches[i]
		       <<"\nRECOMU_mubesttrkType[i]" <<RECOMU_mubesttrkType[i]
		       << endl ;
       */
       //if( ( RECOMU_isTrackerMu[i] && !RECOMU_isGlobalMu[i] ) && !RECOMU_trkmuArbitration[i] ) continue;//arbitrated
       
       if( (RECOMU_isPFMu[i] || (RECOMU_isTrackerHighPtMu[i] && RECOMU_PT[i] > 200.))
	   && ( RECOMU_isGlobalMu[i] || (RECOMU_isTrackerMu[i] && RECOMU_numberOfMatches[i]>0))
	   && RECOMU_mubesttrkType[i]!=2	 
	   && RECOMU_PT[i] > 5. 
	   && fabs(RECOMU_ETA[i]) < 2.4 
	   && fabs(RECOMU_mubesttrkDxy[i]) < .5 && fabs(RECOMU_mubesttrkDz[i]) < 1.
	   )
	 { 
	   iL[N_good] = i ;
	   N_good++ ;
	 }
       
     } //end loop for good muons
     /*
     // if (debug)  
     cout << "\n good muons in the event"<<N_good<< "with indeces : " // this shows the index of each good muon in each event
	  << "\niL[0]: " << iL[0]
	  << "\niL[1]: " << iL[1]
	  << "\niL[2]: " << iL[2]
	  << "\niL[3]: " << iL[3]
	  << "\niL[4]: " << iL[4]
	  << "\niL[5]: " << iL[5]
	  << "\niL[6]: " << iL[6]
	  << "\niL[7]: " << iL[7]
	  << endl ; 
     */
     hN_good_mu->Fill( N_good,newweight );
     
     N_2++ ;  // fill counter
     N_2_w=N_2_w+newweight;
     
     //3) good electron identification no isolation
     
     //we need electron identification to be used in photon identification and cleaning
     
     int Ne_good = 0 ;
     int iLe[8]= { -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1}; //electrons
     
     for( int i = 0; i < RECO_NELE; i++ ){
       /*
       if (debug) cout << "\n Electron i="<< i <<" properties: "
		       << "\nRECOELE_PT[i] " << RECOELE_PT[i]
		       << "\nfabs(RECOELE_ETA[i]) " << fabs(RECOELE_ETA[i])
		       << "\nRECOELE_PFX_rho[i] " << RECOELE_PFX_rho[i]
		       << "\nfabs( RECOELE_SIP[i] ) " << fabs( RECOELE_SIP[i] )
		       << "\nRECOELE_mvaNonTrigV0[i] " << RECOELE_mvaNonTrigV0[i]
		       << "\nfabs( RECOELE_gsftrack_dxy[i] ) " << fabs( RECOELE_gsftrack_dxy[i] )
		       << "\nfabs( RECOELE_gsftrack_dz[i] ) " << fabs( RECOELE_gsftrack_dz[i] )
		       <<"\nRECOELE_gsftrack_expected_inner_hits[i] "<< RECOELE_gsftrack_expected_inner_hits[i]
		       <<"\n(fabs(RECOELE_scl_Eta[i]) "<< fabs(RECOELE_scl_Eta[i])
		       << endl ;
       */
       
       if( RECOELE_PT[i] > 7. && fabs(RECOELE_ETA[i]) < 2.5  /*&& RECOELE_gsftrack_expected_inner_hits[i]<=1 */ ) /* ok */ ;
       else continue ;
       
       bool BDT_ok = 0; // Spring16 with CMSSW_8_2_0
       
       if( RECOELE_PT[i] > 7. &&  RECOELE_PT[i] <= 10. ){
	 if( fabs(RECOELE_scl_Eta[i]) < .8 && RECOELE_mvaNonTrigV0[i] > -0.211 ) BDT_ok = 1 ;
	 if( ( fabs(RECOELE_scl_Eta[i]) >= .8 && fabs(RECOELE_scl_Eta[i]) < 1.479 )
	     && RECOELE_mvaNonTrigV0[i] > -0.396 ) BDT_ok = 1 ;
	 if( fabs(RECOELE_scl_Eta[i]) >= 1.479 && RECOELE_mvaNonTrigV0[i] > -0.215 ) BDT_ok = 1 ;
       }
       else { 
	 if( fabs(RECOELE_scl_Eta[i]) < .8 && RECOELE_mvaNonTrigV0[i] > -0.870 ) BDT_ok = 1 ;
	 if( ( fabs(RECOELE_scl_Eta[i]) >= .8 && fabs(RECOELE_scl_Eta[i]) <= 1.479 )
	     && RECOELE_mvaNonTrigV0[i] > -0.838 ) BDT_ok = 1 ;
	 if( fabs(RECOELE_scl_Eta[i]) > 1.479 && RECOELE_mvaNonTrigV0[i] > -0.763 ) BDT_ok = 1 ;
       }
       
       if( !BDT_ok ) continue ;
       
       if(//fabs( RECOELE_SIP[i] ) < 4. && 
	  fabs(RECOELE_gsftrack_dxy[i]) < .5 
	  && fabs(RECOELE_gsftrack_dz[i]) < 1. ) /* ok */ ;
       else continue ; 
       
       iLe[ Ne_good ] = i ;
       ++Ne_good ;
       
     }// end loop on electrons
     
     hN_good_ele->Fill( Ne_good,newweight );
     hN_good_lep->Fill( N_good + Ne_good,newweight );
     
     /*
     // if (debug) 
     cout << "\n good Electron in the event" <<Ne_good<<" with indeces : " //this will show indices of good electrons in each event 
	  << "\niLe[0]: " << iLe[0]
	  << "\niLe[1]: " << iLe[1]
	  << "\niLe[2]: " << iLe[2]
	  << "\niLe[3]: " << iLe[3]
	  << "\niLe[4]: " << iLe[4]
	  << "\niLe[5]: " << iLe[5]
	  << "\niLe[6]: " << iLe[6]
	  << "\niLe[7]: " << iLe[7]
	  << endl ;
     */
     
     //4) photon definition & cleaning
     //--------------------------------
     //photon preselection: pT > 2 GeV, abs eta <2.4, photon PF relative isolation less than 1.8
     //Supercluster veto: remove all PF photons that match with any elecrton passing loose ID and SIP cuts
     //and matching is according to (abs phi < 2 , abs eta < 0.05) OR (delta R < 0.15), with respect to the electron's supercluster.
     /////////////////// Check/////////////////must modify part add photon PF relative isolation less than 1.8
     
     int iLp[30];
     for( int i = 0 ; i < 30 ; i++ ){iLp[i] = -1;}
     
     int Nphotons = 0;
     
     for( int i = 0; i < RECO_NPFPHOT; i++ ){
       /*
       // if( debug )
       cout << "\n Photon i="<< i <<" properties: "
	    << "\n RECOPFPHOT_PT[i] " << RECOPFPHOT_PT[i]
	    << "\n fabs(RECOPFPHOT_ETA[i]) " << fabs(RECOPFPHOT_ETA[i])
	    << "\n RECOPFPHOT_PHI[i] " << RECOPFPHOT_PHI[i]
	    << "\n RECOPFPHOT_PFX_rho[i] " << RECOPFPHOT_PFX_rho[i]
	    << endl ;
       */
       if ( RECOPFPHOT_PT[i] > 2. && fabs(RECOPFPHOT_ETA[i]) < 2.4 && RECOPFPHOT_PFX_rho[i] < 1.8) {
	 
	 bool is_clean = 1;
	 
	 // cleaning
	 for(int e = 0; e < N_loose_e; e++){
	   if (fabs( RECOELE_SIP[iL_loose_e[e]])>=4.) continue;
	   double deltaPhi = DELTAPHI( RECOPFPHOT_PHI[i] , RECOELE_scl_Phi[iL_loose_e[e]] ) ;
	   double deltaEta = fabs( RECOPFPHOT_ETA[i] - RECOELE_scl_Eta[iL_loose_e[e]] );
	   double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[i] , RECOELE_scl_Phi[iL_loose_e[e]] ),2) + pow(RECOPFPHOT_ETA[i] - RECOELE_scl_Eta[iL_loose_e[e]],2) );
	   
	   if( ( fabs(deltaPhi) < 2 && fabs(deltaEta) < 0.05 ) || deltaR <= 0.15 ){		  
	     //if( debug )cout << "Photon not passing the electron cleaning" << endl;	
	     is_clean = 0;	      
	   }
	 } // end loop on eles		             	
	 
	 if( !is_clean ) continue ;
	 
	 iLp[ Nphotons ] = i ;
	 ++Nphotons ;
	 
       }
     }// end loop on photons
     
     hN_good_phot->Fill( Nphotons,newweight );
     /*
     //if (debug)
     cout << "Number of good photons in the event = "<<Nphotons<<"with indeces: "
	  << "\niLp[0]: " << iLp[0]
	  << "\niLp[1]: " << iLp[1]
	  << "\niLp[2]: " << iLp[2]
	  << "\niLp[3]: " << iLp[3]
	  << "\niLp[4]: " << iLp[4]
	  << "\niLp[5]: " << iLp[5]
	  << "\niLp[6]: " << iLp[6]
	  << "\niLp[7]: " << iLp[7]
	  << "\nNumber of good photons: " << Nphotons
	  << endl ;
     */
     // assign to each photon the closest lepton
     //this assignment done with loose muons and loose electrons
     // but in 7 TeV and 8 TeV was with good muons and good electrons
     
     int iLp_l[30]; //this array will be filled with lepton index with min deltaR with photon
     for( int i = 0 ; i < 30 ; i++ ){iLp_l[i] = -1;}
     int iLp_tagEM[30]; //this array will be filled with lepton type with min deltaR
     for( int i = 0 ; i < 30 ; i++ ){iLp_tagEM[i] = -1;} // tag  0: mu  1: ele
     
     float RECOPFPHOT_DR[30];
     for( int i = 0 ; i < 30 ; i++ ){ RECOPFPHOT_DR[i] = -999; }
     
     
     for( int i = 0; i < Nphotons; i++ ){
       
       double min_deltaR = 1000;
       //double min_deltaR_ET2=1000;
       int  l_min_deltaR = -1;
       int  tag_min_deltaR = -1;   // 0: mu  1: ele
       
       //loose Muon loop
       
       for(int l = 0; l < N_loose_mu; l++){ // loop on muons
	 if (fabs(RECOMU_SIP[iL_loose_mu[l]])>=4.) continue;
	 double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOMU_PHI[iL_loose_mu[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOMU_ETA[iL_loose_mu[l]],2) );
	 if(!(deltaR < 0.5 && deltaR/pow(RECOPFPHOT_PT[iLp[i]],2)<0.012) ) continue;
	 
	 if( deltaR<min_deltaR) { // the closest lepton
	   min_deltaR = deltaR;
	   l_min_deltaR = l;
	   tag_min_deltaR = 0;
	 }
	 
       }//end loop on muons 
       
       
       //loose electron loop 
       
       for(int l = 0; l < N_loose_e; ++l){ // loop on electrons
	 if (fabs(RECOELE_SIP[iL_loose_e[l]])>=4.) continue;
	 double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOELE_PHI[iL_loose_e[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOELE_ETA[iL_loose_e[l]],2) );
	 if(!(deltaR < 0.5 && deltaR/pow(RECOPFPHOT_PT[iLp[i]],2)<0.012) ) continue;
	 
	 if( deltaR < min_deltaR ) { //the closest lepton
	   min_deltaR = deltaR;
	   l_min_deltaR = l;
	   tag_min_deltaR = 1;
	 }  
       }//end loop on electrons  
       
       //  if (debug)  { { if(tag_min_deltaR==0) cout<<"photon with index = "<<iLp[i]<<" has min_deltaR = "<<min_deltaR<<" and min_deltaR_ET2 = "<<min_deltaR_ET2<<" with muon index "<< iL_loose_mu[l_min_deltaR]<<endl;
       //     if(tag_min_deltaR==1)cout<<"photon with index = "<<iLp[i]<<" has min_deltaR = "<<min_deltaR<<" and min_deltaR_ET2 = "<<min_deltaR_ET2<<" with electron index "<< iL_loose_e[l_min_deltaR]<<endl;}}
       
       
       if( min_deltaR < 0.5  ){
	 //iLp_l[ i ] = iL[l_min_deltaR];
	 if (tag_min_deltaR==0) iLp_l[ i ] = iL_loose_mu[l_min_deltaR];
	 if (tag_min_deltaR==1) iLp_l[ i ] = iL_loose_e[l_min_deltaR];
	 iLp_tagEM[ i ] = tag_min_deltaR;
	 RECOPFPHOT_DR[iLp[i]]=min_deltaR; 	
       }
       
     }//end loop on photon for( int i = 0; i < Nphotons; i++ )
     /*
     // if (debug)
     cout << "Indeces of loose leptons associated to photons: "
	  << "\niLp_l[0]: " << iLp_l[0]
	  << "\niLp_l[1]: " << iLp_l[1]
	  << "\niLp_l[2]: " << iLp_l[2]
	  << "\niLp_l[3]: " << iLp_l[3]
	  << "\niLp_l[4]: " << iLp_l[4]
	  << "\niLp_l[5]: " << iLp_l[5]
	  << "\niLp_l[6]: " << iLp_l[6]
	  << "\niLp_l[7]: " << iLp_l[7]
	  << "\niLp_l[7]: " << iLp_l[8]
	  << endl ;
     
     // if (debug)
     cout << "Tag of leptons associated to photons: (0: mu , 1:ele)"
	  << "\niLp_tagEM[0]: " << iLp_tagEM[0]
	  << "\niLp_tagEM[1]: " << iLp_tagEM[1]
	  << "\niLp_tagEM[2]: " << iLp_tagEM[2]
	  << "\niLp_tagEM[3]: " << iLp_tagEM[3]
	  << "\niLp_tagEM[4]: " << iLp_tagEM[4]
	  << "\niLp_tagEM[5]: " << iLp_tagEM[5]
	  << "\niLp_tagEM[6]: " << iLp_tagEM[6]
	  << "\niLp_tagEM[7]: " << iLp_tagEM[7]
	  << "\niLp_tagEM[7]: " << iLp_tagEM[8]
	  << endl ;
     
     if (debug)    { { for(int i=0.;i<Nphotons;i++) {
	   if (iLp_l[i]!=-1 && iLp_tagEM[i]==0) cout << "There is photon with pT= " << RECOPFPHOT_PT[iLp[i]] << " attached to a muon with pT= " << RECOMU_PT[iLp_l[i]] << endl;
	   if (iLp_l[i]!=-1 && iLp_tagEM[i]==1) cout << "There is photon with pT= " << RECOPFPHOT_PT[iLp[i]] << " attached to an electron with pT= " << RECOELE_PT[iLp_l[i]] << endl;
	 }}}
     */
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     //This part to choose only one photon to one lepton, 
     //if we have more than one photon attached to the same lepton we will choose the photon which make min_deltaR_ET2 to the lepton 
     
     
     //Nicola
     
     // Multiple photons associated to the same lepton: the lowest-ΔR(γ,l)/ETγ2 has to be selected.
     
     //Muon part
     
     double min_deltaR_ET2=1000;
     int p_min_deltaR_ET2=-1;
     
     for(int l = 0; l < N_loose_mu; ++l){ // loop on muons
       if (fabs(RECOMU_SIP[iL_loose_mu[l]])>=4.) continue; //loose ID + SIP cut
       min_deltaR_ET2=1000;
       p_min_deltaR_ET2=-1;
       
       for( int p = 0; p < Nphotons; ++p ){
	 if( iLp_l[ p ] == iL_loose_mu[l] && iLp_tagEM[ p ] == 0 )  {
	   // cout <<  "index muon" << iL_loose_mu[l] << endl;
	   double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[p]] , RECOMU_PHI[iL_loose_mu[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[p]] - RECOMU_ETA[iL_loose_mu[l]],2) );
	   double deltaR_ET2 = deltaR/pow(RECOPFPHOT_PT[iLp[p]],2);
	   if (deltaR_ET2<min_deltaR_ET2) {
	     min_deltaR_ET2=deltaR_ET2;
	     RECOPFPHOT_DR[iLp[p]]=deltaR;
	     p_min_deltaR_ET2=p;
	   }
	 }
       }
       
       if (p_min_deltaR_ET2!=-1){
	 for( int p = 0; p < Nphotons; ++p ){
	   if( iLp_l[ p ] == iL_loose_mu[l] && iLp_tagEM[ p ] == 0 )  {
	     if (p!=p_min_deltaR_ET2){
	       iLp_l[ p ] = -1;
	       iLp_tagEM[ p ] = -1;
	     }
	   }
	 }
       }
       
     }
     
     
     //Electron part     
     //
     min_deltaR_ET2=1000;
     p_min_deltaR_ET2=-1;
     
     for(int l = 0; l < N_loose_e; ++l){ // loop on electrons
       if (fabs(RECOELE_SIP[iL_loose_e[l]])>=4.) continue; //loose ID + SIP cut
       min_deltaR_ET2=1000;
       p_min_deltaR_ET2=-1;
       
       for( int p = 0; p < Nphotons; ++p ){
	 if( iLp_l[ p ] == iL_loose_e[l] && iLp_tagEM[ p ] == 1 )  {
	   // cout <<  "index electron" << iL_loose_e[l] << endl;
	   double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[p]] , RECOELE_PHI[iL_loose_e[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[p]] - RECOELE_ETA[iL_loose_e[l]],2));
	   double deltaR_ET2 = deltaR/pow(RECOPFPHOT_PT[iLp[p]],2);
	   // cout << " deltaR_ET2= " << deltaR_ET2 <<endl;
	   if (deltaR_ET2<min_deltaR_ET2){
	     min_deltaR_ET2=deltaR_ET2;
	     RECOPFPHOT_DR[iLp[p]]=deltaR;
	     p_min_deltaR_ET2=p;
	     //cout << " p_min_deltaR_ET2= " << p_min_deltaR_ET2 <<endl;
	   }
	 }	  
       }
       
       if (p_min_deltaR_ET2!=-1){
	 for( int p = 0; p < Nphotons; ++p ){
	   if( iLp_l[ p ] == iL_loose_e[l] && iLp_tagEM[ p ] == 1 )  {
	     if (p!=p_min_deltaR_ET2){
	       iLp_l[ p ] = -1;
	       iLp_tagEM[ p ] = -1;
	     }
	   }
	 }	  
       }
       
     }	
     
     
     /*
     //if( debug )
     cout << "Indeces of loose leptons associated to the photon with lowest DeltaR/ET2: "
	  << "\niLp_l[0]: " << iLp_l[0]
	  << "\niLp_l[1]: " << iLp_l[1]
	  << "\niLp_l[2]: " << iLp_l[2]
	  << "\niLp_l[3]: " << iLp_l[3]
	  << "\niLp_l[4]: " << iLp_l[4]
	  << "\niLp_l[5]: " << iLp_l[5]
	  << "\niLp_l[6]: " << iLp_l[6]
	  << "\niLp_l[7]: " << iLp_l[7]
	  << endl ;
     
     // if( debug )
     cout << "Tag of leptons associated to the photon with lowest DeltaR/ET2: (0: mu , 1:ele)"
	  << "\niLp_tagEM[0]: " << iLp_tagEM[0]
	  << "\niLp_tagEM[1]: " << iLp_tagEM[1]
	  << "\niLp_tagEM[2]: " << iLp_tagEM[2]
	  << "\niLp_tagEM[3]: " << iLp_tagEM[3]
	  << "\niLp_tagEM[4]: " << iLp_tagEM[4]
	  << "\niLp_tagEM[5]: " << iLp_tagEM[5]
	  << "\niLp_tagEM[6]: " << iLp_tagEM[6]
	  << "\niLp_tagEM[7]: " << iLp_tagEM[7]
	  << endl ;
     
     if(debug)  {  for(int i=0.;i<Nphotons;i++) {
	 if (iLp_l[i]!=-1 && iLp_tagEM[i]==0) cout << "There is photon with pT= " << RECOPFPHOT_PT[iLp[i]] << " attached to a muon with pT= " << RECOMU_PT[iLp_l[i]] << endl;
       }}
     */
     
     /////////////////////////////////////////////////////////////////////////////////////////////////// 
     
     // Define a new isolation array to allocate the contribution of photons
     float RECOMU_PFX_dB_new[100],RECOELE_PFX_rho_new[100];
     for (int i=0;i<100;i++){
       RECOMU_PFX_dB_new[i]=RECOMU_PFX_dB[i];
       RECOELE_PFX_rho_new[i]=RECOELE_PFX_rho[i];
     }
     
     // if(debug)
     //{cout<<"Muon isolation before FSR = "<<RECOMU_PFX_dB[0]<<","<<RECOMU_PFX_dB[1]<<","<<RECOMU_PFX_dB[2]<<","<<RECOMU_PFX_dB[3]<<","<<RECOMU_PFX_dB[4]<<","<<RECOMU_PFX_dB[5]<<","<<RECOMU_PFX_dB[6]<<","<<RECOMU_PFX_dB[7]<<","<<RECOMU_PFX_dB[8]<<","<<RECOMU_PFX_dB[9]<<","<<endl;
     //  cout<<"electron isolation before FSR = "<<RECOELE_PFX_rho[0]<<","<<RECOELE_PFX_rho[1]<<","<<RECOELE_PFX_rho[2]<<","<<RECOELE_PFX_rho[3]<<","<<RECOELE_PFX_rho[4]<<","<<RECOELE_PFX_rho[5]<<","<<RECOELE_PFX_rho[6]<<","<<RECOELE_PFX_rho[7]<<","<<RECOELE_PFX_rho[8]<<","<<RECOELE_PFX_rho[9]<<","<<endl;}
     
     
     
     // Exclude that photon from the isolation cone all leptons in the event passing loose ID + SIP cut if it was in the isolation cone and outside the isolation veto (ΔR>0.01 for muons and (ele->supercluster()->eta() < 1.479 || dR > 0.08) for electrons
     
     //if(debug)   cout << "Rho for electron pileup isolation correction is= " << RHO_ele << endl;
     double EffectiveArea=-9999.;
     
     
     for(int i=0.;i<Nphotons;i++) {
       if (iLp_l[i]==-1) continue;
       
       for(int e = 0; e < N_loose_e; ++e){
	 if (fabs( RECOELE_SIP[iL_loose_e[e]])>=4.) continue;
	 // double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOELE_scl_Phi[iL_loose_e[e]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOELE_scl_Eta[iL_loose_e[e]],2) );
	 double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOELE_PHI[iL_loose_e[e]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOELE_ETA[iL_loose_e[e]],2) );         
	 // if(debug)
	 //cout << "deltaR for photon subtraction= " << deltaR << endl;
	 //if(debug) cout<<"loose electrons deltaR = "<<deltaR<<endl;
	 
	 if( deltaR<=0.3 && (RECOELE_scl_Eta[iL_loose_e[e]]< 1.479 || deltaR>0.08) ){ // 0.4 is the isolation cone for electrons in 74x -> 0.3 in 76x              
	   // if( debug )
	   //cout << "Subtracking the photon isolation from the electron isolation value " << endl;
	    
	   EffectiveArea=EAele(RECOELE_scl_Eta[iL_loose_e[e]],tag_2011);
	   RECOELE_PFX_rho_new[iL_loose_e[e]]=
	     (RECOELE_PFchHad[iL_loose_e[e]]+
	      max(0.,RECOELE_PFneuHad[iL_loose_e[e]]+
		  (RECOELE_PFphoton[iL_loose_e[e]]-RECOPFPHOT_PT[iLp[i]] )-
		  max(RHO_ele,0.0)*(EffectiveArea)))/RECOELE_PT[iL_loose_e[e]];	    
	 }
       } // end loop on ele
       
       for(int l = 0; l < N_loose_mu; ++l){ // loop on muons
	 if (fabs(RECOMU_SIP[iL_loose_mu[l]])>=4.) continue;
	 double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOMU_PHI[iL_loose_mu[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOMU_ETA[iL_loose_mu[l]],2) );
	 // if(debug)
	 //cout<<"loose muon with index"<<l<<" has deltaR = "<<deltaR<<endl;
	 
	 if( deltaR<=0.3 && deltaR>0.01){ // 0.3 is the isolation cone for muons in 76x
	   // if( debug )
	   //cout << "Subtracting the photon isolation from the muon isolation value " << endl;
	   RECOMU_PFX_dB_new[iL_loose_mu[l]]=
	     (RECOMU_PFchHad[iL_loose_mu[l]]+
	      max(0.,RECOMU_PFneuHad[iL_loose_mu[l]]+
		  (RECOMU_PFphoton[iL_loose_mu[l]]-RECOPFPHOT_PT[iLp[i]] )-
		  0.5*RECOMU_PFPUchAllPart[iL_loose_mu[l]]))/RECOMU_PT[iL_loose_mu[l]];
	   
	 }
       } // end loop on mu		
     }	//end photon loop
     
     /*
     // if(debug){	
     for(int e = 0; e < N_loose_e; ++e){  if (fabs( RECOELE_SIP[iL_loose_e[e]]>=4.)) continue;
       cout<<" loose electrons isolation before photon sub  = "<<RECOELE_PFX_rho[iL_loose_e[e]]<<" and after sub = "<<RECOELE_PFX_rho_new[iL_loose_e[e]]<<endl;}
     for(int l = 0; l < N_loose_mu; ++l){   if (fabs(RECOMU_SIP[iL_loose_mu[l]])>=4.) continue;
       cout<<" loose muons isolation before photon sub  = "<<RECOMU_PFX_dB[iL_loose_mu[l]]<<" and after sub = "<<RECOMU_PFX_dB_new[iL_loose_mu[l]]<<endl;} 
     
     cout<<"Muon isolation after FSR = "<<RECOMU_PFX_dB_new[0]<<","<<RECOMU_PFX_dB_new[1]<<","<<RECOMU_PFX_dB_new[2]<<","<<RECOMU_PFX_dB_new[3]<<","<<RECOMU_PFX_dB_new[4]<<","<<RECOMU_PFX_dB_new[5]<<","<<RECOMU_PFX_dB_new[6]<<","<<RECOMU_PFX_dB_new[7]<<","<<RECOMU_PFX_dB_new[8]<<","<<RECOMU_PFX_dB_new[9]<<","<<endl;
     cout<<"electron isolation after FSR = "<<RECOELE_PFX_rho_new[0]<<","<<RECOELE_PFX_rho_new[1]<<","<<RECOELE_PFX_rho_new[2]<<","<<RECOELE_PFX_rho_new[3]<<","<<RECOELE_PFX_rho_new[4]<<","<<RECOELE_PFX_rho_new[5]<<","<<RECOELE_PFX_rho_new[6]<<","<<RECOELE_PFX_rho_new[7]<<","<<RECOELE_PFX_rho_new[8]<<","<<RECOELE_PFX_rho_new[9]<<","<<endl;
     //}
     */
     
     
     // *** end FSR
     
     //**** Step 3
     //============
     
     // Build first Z from 2 good muons
     //================================
     
     //initiate structure
     //this data structure (or type ) has 21 information for the Z1 candidate and its leptons 
     // we define a vector with data type of data structure candidateZ
     
     struct candidateZ {
       float massvalue;
       float massvalue_NOFSR;
       int ilept1;
       float pt1;
       float isol1;
       bool ilept1_FSR;
       float eta1;
       float phi1;
       int charge1;
       int lep1tag;
       int lep2tag;
       int charge2;
       float lep1E;	 
       float lep2E;	  
       int ilept2;
       float pt2;
       float isol2;
       bool ilept2_FSR;
       float eta2;
       float phi2;
       float pxZ;
       float pyZ;
       float pzZ;
       float EZ;
       bool withFSR;
       float ptFSR1;
       float ptFSR2;
       int iFSR1;
       int iFSR2;
       int tag;
     };
     
     vector<candidateZ> Zcandvector;
     Zcandvector.clear();

     
 //loop on good muons 
     
     // if( N_good + Ne_good < 2   ) continue ;
     //  if( N_good <2 && Ne_good < 2   ) continue ;
     
     //cout  << "\nStep 3: Number of good Muons: " << N_good << endl;
     //cout  << "\nStep 3: Number of loose Muons: " << N_loose_mu << endl;
     
     int Zxx_tag = 0;    // 1: Zmumu  ,  2: Zee
     int lepton1tag = 0; // 1 for mu ; 2 for ele
     int lepton2tag = 0;
     int i1 = -1; //index of the first lepton (from Z1)
     int j1 = -1; //index of the second lepton (from Z1)
     int pi1 = -1; 
     int pj1 = -1;
     
     bool has_FSR_Z1 = 0;
     TLorentzVector Lepton1,Lepton2,DiLepton,LeptonCorrection;
     
     
     for(int i =0 ; i < N_good ; i++){
       for(int j =i+1 ; j < N_good ; j++){
	 
	 
	 if (fabs(RECOMU_SIP[iL[i]])>=4.) continue;
	 if (fabs(RECOMU_SIP[iL[j]])>=4.) continue;

	 if (fabs(RECOMU_PFX_dB_new[iL[i]])>=0.35) continue; // Isolation
	 if (fabs(RECOMU_PFX_dB_new[iL[j]])>=0.35) continue;
	 
	 if(RECOMU_CHARGE[ iL[j] ] == RECOMU_CHARGE[ iL[i] ]) continue; // opposite charge
	 
	 double pxZ, pyZ, pzZ;
	 double EZ;
	 double massZ;
	 double massZ_noFSR = 0;
	 
	 double ptZ = 0;
	 double Y_Z = -9;
	 double sum_ptZ = 0.;
	 
	 int tempphotid1=-1;
	 int tempphotid2=-1;
	 int templepid=-1;
	 
	 float pTphot1=-999.;
	 float pTphot2=-999.;
	 
	 //Zmass with noFSR
	 //-------------------
	 
	 //if(debug)  cout << "\n Looking for a new pair"<< endl;
	 //if (debug) cout << "pt,eta,phi,charge= " << RECOMU_PT[ iL[i] ] << " " << RECOMU_ETA[ iL[i] ] << " " << RECOMU_PHI[ iL[i] ] << " " << RECOMU_CHARGE[ iL[i] ]<< endl;
	 //if (debug) cout << "pt,eta,phi,charge= " << RECOMU_PT[ iL[j] ] << " " << RECOMU_ETA[ iL[j] ] << " " << RECOMU_PHI[ iL[j] ] << " " << RECOMU_CHARGE[ iL[j] ]<< endl;
	 
	 //Evaluate the mass 
	 
	 Lepton1.SetPtEtaPhiM (RECOMU_PT[iL[i]],RECOMU_ETA[iL[i]], RECOMU_PHI[iL[i]], 0.105);
	 Lepton2.SetPtEtaPhiM (RECOMU_PT[iL[j]],RECOMU_ETA[iL[j]], RECOMU_PHI[iL[j]], 0.105);
	 
	 DiLepton=Lepton1+Lepton2;	  
	 massZ = DiLepton.M();
	 ptZ = DiLepton.Pt();
	 pxZ = DiLepton.Px();
	 pyZ = DiLepton.Py();
	 pzZ = DiLepton.Pz();
	 EZ = DiLepton.E();
	 Y_Z = DiLepton.Rapidity();
	 
	 sum_ptZ = RECOMU_PT[ iL[i] ] + RECOMU_PT[ iL[j] ];
	 massZ_noFSR = massZ;
	 
	 Zxx_tag=1;//Z decay to 2 good muons
	 lepton1tag = 1; // 1 mu , 2 ele
	 lepton2tag = 1;
	 
	 /*
	 if(debug)   {  cout<<"step NoFSR"<<endl;
	   cout << "Mass Z= " << massZ << " rapidity Z= " << Y_Z << endl;
	   cout<<"Lepton1,2 pt = "<<Lepton1.Pt()<<" , "<<Lepton2.Pt()<<endl;}
	 */
	 // This definition how to calculate px py pz in lorentz vector
	 
	 /* EZ = RECOMU_E[ iL[i] ] + RECOMU_E[ iL[j] ];
	    pxZ = RECOMU_PT[ iL[i] ]*cos( RECOMU_PHI[ iL[i] ] ) + RECOMU_PT[ iL[j] ]*cos( RECOMU_PHI[ iL[j] ] );
	    pyZ = RECOMU_PT[ iL[i] ]*sin( RECOMU_PHI[ iL[i] ] ) + RECOMU_PT[ iL[j] ]*sin( RECOMU_PHI[ iL[j] ] );
	    pzZ = RECOMU_P[ iL[i] ]*cos( RECOMU_THETA[ iL[i] ] ) + RECOMU_P[ iL[j] ]*cos( RECOMU_THETA[ iL[j] ] );
	    
	    ptZ = sqrt( pxZ*pxZ + pyZ*pyZ );
	    Y_Z = 0.5 * log ( (EZ + pzZ)/(EZ - pzZ) );
	    sum_ptZ = RECOMU_PT[ iL[i] ] + RECOMU_PT[ iL[j] ];*/
	 
	 // ** Association of FSR to Z
	 //------------------------------
	 
	 //if( debug ) cout  << "Step Z+FSR  " << endl;
	 
	 bool has_FSR_Z = 0;
	 int N_FSR_Z = 0;
	 double max_pt_FSR_Z = -1.;
	 int pi = -1; 
	 int pj = -1;
	 double mllp=-1;
	 
	 //here we will make loop on photons to know how many FSR in the event is it N_FSR =1 or > 1 and to know photon max pt
	 // we will loop on photons that made small deltaR with muons
	
	 for( int p = 0; p < Nphotons; p++ ){
	   
	   if( iLp_l[ p ] == iL[i] && iLp_tagEM[ p ] == 0 )  {
	     
	     // evaluate the mass
	     
	     
	     LeptonCorrection.SetPtEtaPhiM (RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]], RECOPFPHOT_PHI[iLp[p]], 0);
	     
	     Lepton1= Lepton1+ LeptonCorrection;
	     DiLepton=Lepton1+Lepton2;
	     
	     mllp=DiLepton.M();
	     ptZ = DiLepton.Pt();
	     pxZ = DiLepton.Px();
	     pyZ = DiLepton.Py();
	     pzZ = DiLepton.Pz();
	     EZ = DiLepton.E();
	     Y_Z = DiLepton.Rapidity();
	     
	     has_FSR_Z = 1; 
	     pi = p; 
	     N_FSR_Z++;
	     massZ=mllp;
	     if( RECOPFPHOT_PT[iLp[p]] > max_pt_FSR_Z ) max_pt_FSR_Z = RECOPFPHOT_PT[iLp[p]];
	     /*
	     if(debug) { cout<<"Z has FSR with index = "<<iLp[p]<<"pt = "<<RECOPFPHOT_PT[iLp[p]]<<" with first muon with index= "<<iL[i]<<"pt = "<<RECOMU_PT[iL[i]]<<endl;
	       cout<<"first muon after correction pt "<<Lepton1.Pt()<<endl;
	       cout<<"secomd muon pt = "<<Lepton2.Pt()<<endl;
	       cout<<"massZ = "<<massZ<<endl;}
	     */
	     
	   }//end of if( iLp_l[ p ] == iL[i] && iLp_tagEM[ p ] == 0 )
	   
	   if( iLp_l[ p ] == iL[j] && iLp_tagEM[ p ] == 0 )  {
	     
	     // evaluate the mass
	     
	     LeptonCorrection.SetPtEtaPhiM (RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]], RECOPFPHOT_PHI[iLp[p]], 0);
	     
	     Lepton2=Lepton2+LeptonCorrection;
	     DiLepton=Lepton1+Lepton2;
	     
	     mllp=DiLepton.M();
	     ptZ = DiLepton.Pt();
	     pxZ = DiLepton.Px();
	     pyZ = DiLepton.Py();
	     pzZ = DiLepton.Pz();
	     EZ = DiLepton.E();
	     Y_Z = DiLepton.Rapidity();
	     
	     has_FSR_Z = 1; 
	     pj = p; 
	     N_FSR_Z++;
	     massZ=mllp;
	     if( RECOPFPHOT_PT[iLp[p]] > max_pt_FSR_Z ) max_pt_FSR_Z = RECOPFPHOT_PT[iLp[p]];
	     /*
	     if(debug)  {cout<<"Z has FSR with index = "<<iLp[p]<<"pt = "<<RECOPFPHOT_PT[iLp[p]]<<" with second muon with index= "<<iL[j]<<"pt = "<<RECOMU_PT[iL[j]]<<endl;
	       cout<<"second muon after correction pt "<<Lepton2.Pt()<<endl;
	       cout<<"first muon pt "<<Lepton1.Pt();
	       cout<<"massZ = "<<massZ<<endl;}
	     */
	   }//end of if( iLp_l[ p ] == iL[j] && iLp_tagEM[ p ] == 0 )
	   
	 }// end loop on FSR photons for( int p = 0; p < Nphotons; p++ )
	 /*
	 if( debug && has_FSR_Z) {
	   cout  << " Z has FSR! " << endl;
	   cout  << "  N_FSR_Z " << N_FSR_Z << endl;
	   cout  << "  max_pt of photon FSR_Z " << max_pt_FSR_Z << endl;
	   if( pi > -1 ) cout  << "  pi " << pi << " --> index photon: " << iLp[pi] << " associated lepton: " << iLp_l[pi] << " (= "<< iL[i]<<" ? )  tag: " << iLp_tagEM[pi] << endl;
	   if( pj > -1 ) cout  << "  pj " << pj << " --> index photon: " << iLp[pj] << " associated lepton: " << iLp_l[pj] << " (= "<< iL[j]<<" ? )  tag: " << iLp_tagEM[pj] << endl;
	 }
	 else  {
	   if(debug)  cout << "No FSR photon attached" << endl;
	 }
	 */
	 
	 if(has_FSR_Z){ //there are 2 options N_FSR_Z==1 or >1
	   
	   ++N_3_FSR; // fill the counter
	   N_3_FSR_w=N_3_FSR_w+newweight;
	   
	   // do not recompute isolation here
	   //differ from 7TeV and 8TeV Runs
	   //isolation were calculated here to to substract FSR photon from  lepton isolation cone
	   /*
	   if( debug ) cout  << "Z Isolation (corrected for photon): "
			     << "\n RECOMU_PFX_dB_new[ iL[i] ] " << RECOMU_PFX_dB_new[ iL[i] ]
			     << "\n RECOMU_PFX_dB_new[ iL[j] ] " << RECOMU_PFX_dB_new[ iL[j] ]
			     << endl;
	   */
	   
	   if( pi != -1 ){
	     pTphot1=RECOPFPHOT_PT[iLp[pi]];
	     tempphotid1=iLp[pi];
	   }
	   if( pj != -1 ){
	     pTphot2=RECOPFPHOT_PT[iLp[pj]];
	     tempphotid2=iLp[pj];	      
	   }
	   
	 }//end of if(has_FSR_Z)
	 
	 
	 else { //if don't have FSR
	   
	   
	   pi=-1;
	   pj=-1;
	   /*
	   if( debug ) cout  << "Z Isolation: "  
			     << "\n RECOMU_PFX_dB_new[ iL[i] ] " << RECOMU_PFX_dB_new[ iL[i] ]
			     << "\n RECOMU_PFX_dB_new[ iL[j] ] " << RECOMU_PFX_dB_new[ iL[j] ]
			     << endl; 
	   */
	 }//end of else if don't have FSR
	 
	 // ***** end association of FSR to Z
	 
	 //if (debug) cout << "Filling a struct for Z" << endl; 
	 
	 // here we still in the 2 muons loop
	 // here we  define Z it is a pointer to data structure candidateZ
	 // In one event we have different combinations for Z as we still in the 2 muon loop
	 //the first muon with the second muon and third and fourth and so on
	 //we fill the vector Zcandvector with this 21 information as the first element of the vector for the first combinations of Z
	 //the second element of the vector is the 21 informations from second Z candidate in the same event and so on for all combinations
	 
	 candidateZ *Z = new candidateZ; //this Z is a pointer to data structure candidateZ 
	 Z->massvalue=massZ;
	 Z->massvalue_NOFSR=massZ_noFSR;
	 Z->ilept1=iL[i];
	 Z->ilept2=iL[j];
	 Z->pt1=RECOMU_PT[iL[i]];
	 Z->pt2=RECOMU_PT[iL[j]];
	 Z->eta1=RECOMU_ETA[iL[i]];
	 Z->eta2=RECOMU_ETA[iL[j]];
	 Z->phi1=RECOMU_PHI[iL[i]];
	 Z->phi2=RECOMU_PHI[iL[j]];
	 Z->charge1=RECOMU_CHARGE[iL[i]];
	 Z->charge2=RECOMU_CHARGE[iL[j]];
	 Z->lep1E=RECOMU_E[iL[i]];
	 Z->lep2E=RECOMU_E[iL[j]];
	 Z->lep1tag=lepton1tag;
	 Z->lep2tag=lepton2tag;
	 Z->isol1=RECOMU_PFX_dB_new[ iL[i] ];
	 Z->isol2=RECOMU_PFX_dB_new[ iL[j] ];
	 if( pi != -1 ) {Z->ilept1_FSR=true;
	   /*Z->ilept2_FSR=false;*/}
	 if( pj != -1 ) {/*Z->ilept1_FSR=false;*/
	   Z->ilept2_FSR=true;}
	 Z->pxZ=pxZ;
	 Z->pyZ=pyZ;
	 Z->pzZ=pzZ;
	 Z->EZ=EZ;
	 if( has_FSR_Z ) {
	   Z->withFSR=1;
	   Z->ptFSR1=pTphot1;
	   Z->ptFSR2=pTphot2;
	   Z->iFSR1=tempphotid1;
	   Z->iFSR2=tempphotid2;
	 }	      
	 else {
	   Z->withFSR=0;
	   Z->ptFSR1=0.;
	   Z->ptFSR2=0.;
	   Z->iFSR1=-1;
	   Z->iFSR2=-1;
	   Z->ilept1_FSR=false;
	   Z->ilept2_FSR=false;
	 }
	 
	 Z->tag=Zxx_tag;
	 Zcandvector.push_back(*Z);	  
	 
	 
       }//end of loop on second good muon for(int j =i+1 ; j < N_good ; j++)
     }//end of loop on first good muon for(int i =0 ; i < N_good ; i++)
     
     //end loop on muon couples 
     
     
     
     // Build first Z from 2 good electrons
     //====================================
     
     
     for(int i =0 ; i < Ne_good ; i++){
       for(int j =i+1 ; j < Ne_good ; j++){
	 
	 if (fabs(RECOELE_SIP[iLe[i]])>=4.) continue;
	 if (fabs(RECOELE_SIP[iLe[j]])>=4.) continue;
	 
	 if (fabs(RECOELE_PFX_rho_new[iLe[i]])>=0.35) continue; // Isolation cut
	 if (fabs(RECOELE_PFX_rho_new[iLe[j]])>=0.35) continue;
	 
	 if(RECOELE_CHARGE[ iLe[j] ] == RECOELE_CHARGE[ iLe[i] ]) continue; // opposite charge
	 
	 double pxZ, pyZ, pzZ;
	 double EZ;
	 double massZ;
	 double massZ_noFSR = 0;
	 
	 double ptZ = 0;
	 double Y_Z = -9;
	 double sum_ptZ = 0.;
	 
	 int tempphotid1=-1;
	 int tempphotid2=-1;
	 int templepid=-1;
	 
	 float pTphot1=-999.;
	 float pTphot2=-999.;
	 
	 //Zmass with noFSR
	 //-------------------
	 //if (debug) cout << "\n Looking for a new pair"<< endl;
	 //if (debug) cout << "pt,eta,phi,charge= " << RECOELE_PT[ iLe[i] ] << " " << RECOELE_ETA[ iLe[i] ] << " " << RECOELE_PHI[ iLe[i] ] << " " << RECOELE_CHARGE[ iLe[i] ]<< endl;
	 //if (debug) cout << "pt,eta,phi,charge= " << RECOELE_PT[ iLe[j] ] << " " << RECOELE_ETA[ iLe[j] ] << " " << RECOELE_PHI[ iLe[j] ] << " " << RECOELE_CHARGE[ iLe[j] ]<< endl;
	 
	 //Evaluate the mass 
	 
	 Lepton1.SetPtEtaPhiM (RECOELE_PT[iLe[i]],RECOELE_ETA[iLe[i]], RECOELE_PHI[iLe[i]], 0.000511);
	 Lepton2.SetPtEtaPhiM (RECOELE_PT[iLe[j]],RECOELE_ETA[iLe[j]], RECOELE_PHI[iLe[j]], 0.000511);
	 
	 DiLepton=Lepton1+Lepton2;	  
	 massZ = DiLepton.M();
	 ptZ = DiLepton.Pt();
	 pxZ = DiLepton.Px();
	 pyZ = DiLepton.Py();
	 pzZ = DiLepton.Pz();
	 EZ = DiLepton.E();
	 Y_Z = DiLepton.Rapidity();
	 
	 sum_ptZ = RECOELE_PT[ iLe[i] ] + RECOELE_PT[ iLe[j] ];
	 massZ_noFSR = massZ;
	 
	 Zxx_tag=2;//Z decay to electrons
	 lepton1tag = 2; // 1 mu , 2 ele
	 lepton2tag = 2;
	 
	 //if (debug) cout << "Mass Z= " << massZ << " rapidity Z= " << Y_Z << endl;
	 
	 
	 // This definition how to calculate px py pz in lorentz vector
	 
	 /* EZ = RECOELE_E[ iLe[i] ] + RECOELE_E[ iLe[j] ];
	    pxZ = RECOELE_PT[ iLe[i] ]*cos( RECOELE_PHI[ iLe[i] ] ) + RECOELE_PT[ iLe[j] ]*cos( RECOELE_PHI[ iLe[j] ] );
	    pyZ = RECOELE_PT[ iLe[i] ]*sin( RECOELE_PHI[ iLe[i] ] ) + RECOELE_PT[ iLe[j] ]*sin( RECOELE_PHI[ iLe[j] ] );
	    pzZ = RECOELE_P[ iLe[i] ]*cos( RECOELE_THETA[ iLe[i] ] ) + RECOELE_P[ iLe[j] ]*cos( RECOELE_THETA[ iLe[j] ] );
	    
	    ptZ = sqrt( pxZ*pxZ + pyZ*pyZ );
	    Y_Z = 0.5 * log ( (EZ + pzZ)/(EZ - pzZ) );
	    sum_ptZ = RECOELE_PT[ iLe[i] ] + RECOELE_PT[ iLe[j] ];*/
	 
	 // ** Association of FSR to Z
	 //------------------------------
	 
	 //if( debug ) cout  << "Step Z+FSR  " << endl;
	 
	 bool has_FSR_Z = 0;
	 int N_FSR_Z = 0;
	 double max_pt_FSR_Z = -1.;
	 int pi = -1; 
	 int pj = -1;
	 double mllp=-1;
	 
	 //here we will make loop on photons to know how many FSR in the event is it N_FSR =1 or > 1 and to know photon max pt
	 // we will loop on photons that made small deltaR with muons
	 
	 for( int p = 0; p < Nphotons; p++ ){
	   
	   if( iLp_l[ p ] == iLe[i] && iLp_tagEM[ p ] == 1 )  {
	     
	     // evaluate the mass
	     
	     
	     LeptonCorrection.SetPtEtaPhiM (RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]], RECOPFPHOT_PHI[iLp[p]], 0);
	     Lepton1=Lepton1+LeptonCorrection;
	     DiLepton=Lepton1+Lepton2;
	     
	     mllp = DiLepton.M();
	     ptZ = DiLepton.Pt();
	     pxZ = DiLepton.Px();
	     pyZ = DiLepton.Py();
	     pzZ = DiLepton.Pz();
	     EZ = DiLepton.E();
	     Y_Z = DiLepton.Rapidity();
	     
	     has_FSR_Z = 1; 
	     pi = p; 
	     N_FSR_Z++;
	     massZ=mllp;
	     if( RECOPFPHOT_PT[iLp[p]] > max_pt_FSR_Z ) max_pt_FSR_Z = RECOPFPHOT_PT[iLp[p]];
	     /*
	     if(debug)    {cout<<"Z has FSR with index = "<<iLp[p]<<"pt = "<<RECOPFPHOT_PT[iLp[p]]<<" with first electron with index= "<<iLe[i]<<"pt = "<<RECOELE_PT[iLe[i]]<<endl;
	       cout<<"first electron after correction pt "<<Lepton1.Pt()<<endl;
	       cout<<"secomd electron pt = "<<Lepton2.Pt()<<endl;
	       cout<<"massZ = "<<massZ<<endl;}
	     */
	   }//end of if( iLp_l[ p ] == iLe[i] && iLp_tagEM[ p ] == 1 )
	   
	   if( iLp_l[ p ] == iLe[j] && iLp_tagEM[ p ] == 1 )  {
	     
	     // evaluate the mass
	     
	     LeptonCorrection.SetPtEtaPhiM (RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]], RECOPFPHOT_PHI[iLp[p]], 0);
	     
	     Lepton2=Lepton2+LeptonCorrection;
	     DiLepton=Lepton1+Lepton2;
	     
	     mllp = DiLepton.M();	
	     ptZ = DiLepton.Pt();
	     pxZ = DiLepton.Px();
	     pyZ = DiLepton.Py();
	     pzZ = DiLepton.Pz();
	     EZ = DiLepton.E();
	     Y_Z = DiLepton.Rapidity();
	     
	     has_FSR_Z = 1; 
	     pj = p; 
	     N_FSR_Z++;
	     massZ=mllp;
	     if( RECOPFPHOT_PT[iLp[p]] > max_pt_FSR_Z ) max_pt_FSR_Z = RECOPFPHOT_PT[iLp[p]];
	     /*
	     if(debug) {cout<<"Z has FSR with index = "<<iLp[p]<<"pt = "<<RECOPFPHOT_PT[iLp[p]]<<" with second electron with index= "<<iLe[j]<<"pt = "<<RECOELE_PT[iLe[j]]<<endl;
	       cout<<"second electron after correction pt "<<Lepton2.Pt()<<endl;
	       cout<<"first electron pt "<<Lepton1.Pt();
	       cout<<"massZ = "<<massZ<<endl;}
	     */
	     
	   }//end of if( iLp_l[ p ] == iLe[j] && iLp_tagEM[ p ] == 1 )
	   
	 }// end loop on FSR photons for( int p = 0; p < Nphotons; p++ )
	 
	 /*
	 if( debug && has_FSR_Z) {
	   cout  << " Z has FSR! " << endl;
	   cout  << "  N_FSR_Z " << N_FSR_Z << endl;
	   cout  << "  max_pt of photon FSR_Z " << max_pt_FSR_Z << endl;
	   if( pi > -1 ) cout  << "  pi " << pi << " --> index photon: " << iLp[pi] << " associated lepton: " << iLp_l[pi] << " (= "<< iLe[i]<<" ? )  tag: " << iLp_tagEM[pi] << endl;
	   if( pj > -1 ) cout  << "  pj " << pj << " --> index photon: " << iLp[pj] << " associated lepton: " << iLp_l[pj] << " (= "<< iLe[j]<<" ? )  tag: " << iLp_tagEM[pj] << endl;
	 }
	 else  {
	   if(debug)   cout << "No FSR photon attached" << endl;
	 }
	 */
	 
	 
	 if(has_FSR_Z){ //there are 2 options N_FSR_Z==1 or >1
	   
	   ++N_3_FSR; // fill the counter
	   N_3_FSR_w=N_3_FSR_w+newweight;
	   
	   
	   // do not recompute isolation here
	   //differ from 7TeV and 8TeV Runs
	   //isolation were calculated here to to substract FSR photon from  lepton isolation cone
	   /*
	   if( debug ) cout  << "Z Isolation ( corrected for photon): "
			     << "\n RECOELE_PFX_rho_new[ iLe[i] ] " << RECOELE_PFX_rho_new[ iLe[i] ]
			     << "\n RECOELE_PFX_rho_new[ iLe[j] ] " << RECOELE_PFX_rho_new[ iLe[j] ]
			     << endl;
	   */
	   if( pi != -1 ){
	     pTphot1=RECOPFPHOT_PT[iLp[pi]];
	     tempphotid1=iLp[pi];	      
	   }
	   if( pj != -1 ){	      
	     pTphot2=RECOPFPHOT_PT[iLp[pj]];
	     tempphotid2=iLp[pj];	      
	   }
	 }//end of if(has_FSR_Z)
	 
	 else { //if don't have FSR
	   
	   pi=-1;
	   pj=-1;
	   
	   /*
	   if( debug ) cout  << "Z Isolation: "  
			     << "\n RECOELE_PFX_rho_new[ iLe[i] ] " << RECOELE_PFX_rho_new[ iLe[i] ]
			     << "\n RECOELE_PFX_rho_new[ iLe[j] ] " << RECOELE_PFX_rho_new[ iLe[j] ]
			     << endl;
	   */
	 }//end of else if don't have FSR
	 
	 // ***** end association of FSR to Z
	 
	 //if(debug) cout << "2e2mu: " << Zxx_tag << endl;
	 //if(debug) cout << "Filling a struct for Z" << endl; 
	 
	 
	 candidateZ *Z = new candidateZ; //this Z is a pointer to data structure candidateZ 
	 Z->massvalue=massZ;
	 Z->massvalue_NOFSR=massZ_noFSR;
	 Z->ilept1=iLe[i];
	 Z->ilept2=iLe[j];
	 Z->pt1=RECOELE_PT[iLe[i]];
	 Z->pt2=RECOELE_PT[iLe[j]];
	 Z->eta1=RECOELE_ETA[iLe[i]];
	 Z->eta2=RECOELE_ETA[iLe[j]];
	 Z->phi1=RECOELE_PHI[iLe[i]];
	 Z->phi2=RECOELE_PHI[iLe[j]];
	 Z->charge1=RECOELE_CHARGE[iLe[i]];
	 Z->charge2=RECOELE_CHARGE[iLe[j]];
	 Z->lep1E=RECOELE_E[iLe[i]];
	 Z->lep2E=RECOELE_E[iLe[j]];
	 Z->lep1tag=lepton1tag;
	 Z->lep2tag=lepton2tag;
	 Z->isol1=RECOELE_PFX_rho_new[ iLe[i] ];
	 Z->isol2=RECOELE_PFX_rho_new[ iLe[j] ];
	 if( pi != -1 ) {Z->ilept1_FSR=true;
	   /* Z->ilept2_FSR=false;*/}
	 if( pj != -1 ) {/*Z->ilept1_FSR=false;*/
	   Z->ilept2_FSR=true;}
	 Z->pxZ=pxZ;
	 Z->pyZ=pyZ;
	 Z->pzZ=pzZ;
	 Z->EZ=EZ;
	 if( has_FSR_Z ) {
	   Z->withFSR=1;
	   Z->ptFSR1=pTphot1;
	   Z->ptFSR2=pTphot2;
	   Z->iFSR1=tempphotid1;
	   Z->iFSR2=tempphotid2;
	 }	      
	 else {
	   Z->withFSR=0;
	   Z->ptFSR1=0.;
	   Z->ptFSR2=0.;
	   Z->iFSR1=-1;
	   Z->iFSR2=-1;
	   Z->ilept1_FSR=false;
	   Z->ilept2_FSR=false;
	 }
	 
	 Z->tag=Zxx_tag;	 	  
	 Zcandvector.push_back(*Z);	  
	 
	 
       }//end of loop on second good muon for(int j =i+1 ; j < N_good ; j++)
     }//end of loop on first good muon for(int i =0 ; i < N_good ; i++)
     
     //end loop on muon couples 
     /*
     if(Zcandvector.size()==0)continue;
     
     //  cout<<"========================================="<<endl;
     
     
     cout<<"number of Zs filled in vector Zcandvector  = "<<Zcandvector.size()<<endl;
     
     for(int i=0; i<Zcandvector.size(); i++){
       
       cout<<"Zcanvector properties tag = "<<Zcandvector.at(i).tag<<" ,mass =  "<<Zcandvector.at(i).massvalue<<",FSR = "<<Zcandvector.at(i).withFSR <<"photon indeces = "<<Zcandvector.at(i).iFSR1<<","<<Zcandvector.at(i).iFSR2<<" ,photons pT = "<<Zcandvector.at(i).ptFSR1<<" , "<<Zcandvector.at(i).ptFSR2<<" ,leptons FSR = "<<Zcandvector.at(i).ilept1_FSR<<" , "<<Zcandvector.at(i).ilept2_FSR<<" and leptons indices = "<<Zcandvector.at(i).ilept1<<","<<Zcandvector.at(i).ilept2<<"and pt = "<<Zcandvector.at(i).pt1<<","<<Zcandvector.at(i).pt2<<" and leptons isolation = "<<Zcandvector.at(i).isol1<<", "<<Zcandvector.at(i).isol2<<endl;
       
     }
     */
     
//====================================================================================================================================================================================================================================================================================================================//     
// Z candidates ready here - allows further cuts for FR calculation and ZZ/Higgs system building for CR's
//====================================================================================================================================================================================================================================================================================================================//
     
//=========================================================================================================================================================//
//  Start Fake rate calculation //     
//=========================================================================================================================================================//
     
     // Z1 selection
     //==============
     double pxZ1 = 0;  //Z1 kinematics
     double pyZ1 = 0;
     double pzZ1 = 0;
     double ptZ1 = 0;
     double EZ1 = 0;
     double Y_Z1 = -9;
     double massZ1 = 0;
     double massZ1_noFSR = 0;
     double sum_ptZ1 = 0.;
     
     double Z1pt1 =-999.;
     double Z1pt2 =-999.;
     double Z1eta1 = -999.;
     double Z1eta2 = -999.;
     double Z1phi1 = -999.;
     double Z1phi2 = -999.;
     double Z1charge1 = -999.;
     double Z1charge2 = -999.;
     
     int indexlep1Z1 = -1;
     int indexlep2Z1 = -1;
     int indexZ1= -1;
     int Z1tag=-999;
     int Z1lep1tag = -999;
     int Z1lep2tag = -999;
     float Z1lep1E = -999;
     float Z1lep2E = -999;
     
     int Z1withfsr = -1;
     int Z1fsr1 = -999;
     int Z1fsr2 = -999;
     double Z1ptFSR1 = -999;
     double Z1ptFSR2 = -999;
     bool Z1islep1fsr = false;
     bool Z1islep2fsr = false;
     
     
     
     // Choice of Z1 as the closest to the Z mass for FR calculation only
     //loop over all Zs in the vector "Zcandvector" to have Z closest to nominal Z mass
     //-------------------------------------------------------------------------------------
     
     
     candidateZ Z1cand = {0., -999., -999, -999., -999., false , -999., -999., -999, -999, -999, -999, -999., -999., -999, -999., -999., false, -999., -999., -999., -999., -999., -999., false, -999., -999., -999, -999 , -999};
     
     candidateZ Z1 = {0., -999., -999, -999., -999., false , -999., -999., -999, -999, -999, -999, -999., -999., -999, -999., -999., false, -999., -999., -999., -999., -999., -999., false, -999., -999., -999, -999 , -999};
     
     for (int i=0;i<Zcandvector.size();i++){
       
       //cout<<"Zcandvector massvalue = "<<Zcandvector.at(i).massvalue<<endl;
       //cout<<"Z1cand massvalue = "<<Z1cand.massvalue<<endl;
       
       if( fabs(Zcandvector.at(i).massvalue - Zmass) < fabs(Z1cand.massvalue - Zmass) ){
	 
	 Z1cand =  Zcandvector.at(i);
	 
	 //cout<<"Z1cand mass = "<<Z1cand.massvalue<<endl;
	 
       }
     } 
     
     Z1 = Z1cand;
     
     //cout<<"the Z1 has mass = "<<Z1.massvalue<<" and leptons with indices =  "<<Z1.ilept1<<","<<Z1.ilept2<<" and tag = "<<Z1.tag<<"and leptons energy "<<Z1.lep1E<<", "<<Z1.lep2E<<endl;
     
     //  R_3++;
     //  R_3_w=R_3_w+newweight;
     
     hMZ1_3_no_effW->Fill( Z1.massvalue,newweight);
     hPFMET_3_no_effW->Fill(RECO_PFMET,newweight);
     
     
     //efficiency weight
     //-------------------
     
     //cout<<"efficiency weight step after choose ZZ pair "<<endl;
     
     Double_t eff_weight_3 = 1.;
     
     //   int indexlep1Z1, indexlep2Z1;
     
     int z1lept[2] = {Z1.ilept1 , Z1.ilept2 };
     
     if (Z1.tag==1){ //Z1 is mumu
       
       //cout<<"Z1 is mumu "<<endl;
       
       for(int i = 0; i < 2; ++i){
	 
	 Double_t Pt = RECOMU_PT[ z1lept[i] ]; 
	 Double_t Eta = RECOMU_ETA[ z1lept[i] ];
	 
	 //cout<<"Z1 has muon = "<<z1lept[i]<<" with pt = "<<Pt<<" and eta = "<<Eta<<endl;
	 
	 if( (MC_type == "Spring16" || MC_type == "Moriond17" ) && DATA_type == "NO"){
	   
	   Int_t binx = mu_scale_2016->GetXaxis()->FindBin(Eta);
	   //cout<<"binx = "<<binx<<endl;
	   Int_t biny = mu_scale_2016->GetYaxis()->FindBin(Pt);
	   //cout<<"biny = "<<biny<<endl;
	   //cout<<"muon scale factor = "<<mu_scale_2016->GetBinContent(binx,biny)<<endl;
	   if (mu_scale_2016->GetBinContent(binx,biny)>0.) eff_weight_3*=mu_scale_2016->GetBinContent(binx,biny);
	   //cout<<"eff_weight_3 = "<<eff_weight_3<<endl;
	 }
       }
     }
     else if (Z1.tag==2){ //Z1 is ee
       
       //cout<<"Z1 is ee "<<endl;
       
       
       for(int i = 0; i < 2; ++i){
	 
	 Double_t Pt = RECOELE_PT[ z1lept[i] ]; 
	 Double_t Eta = RECOELE_ETA[ z1lept[i] ];
	 
	 //cout<<"Z1 has electron = "<<z1lept[i]<<" with pt = "<<Pt<<" and eta = "<<Eta<<endl;
	 
	 if( ( MC_type == "Spring16" || MC_type == "Moriond17" ) && DATA_type == "NO"){
	   
	   if(RECOELE_isGap[ z1lept[i] ]==0){
	     Int_t binx = ele_scale_factors2016->GetXaxis()->FindBin(Eta);
	     //cout<<"binx = "<<binx<<endl;
	     Int_t biny = ele_scale_factors2016->GetYaxis()->FindBin(Pt);
	     //cout<<"biny = "<<biny<<endl;
	     //cout<<"ele scale factor2016 = "<<ele_scale_factors2016->GetBinContent(binx,biny)<<endl;
	      if (ele_scale_factors2016->GetBinContent(binx,biny)>0.) eff_weight_3*=ele_scale_factors2016->GetBinContent(binx,biny);
	      //cout<<"eff_weight_3 = "<<eff_weight_3<<endl;
	      
	   }
	   else if(RECOELE_isGap[ z1lept[i] ]==1){
	     //cout<<"crack electron"<<endl;
	     Int_t binx = ele_scale_factors_gap2016->GetXaxis()->FindBin(Eta);
	     //cout<<"binx = "<<binx<<endl;
	     Int_t biny = ele_scale_factors_gap2016->GetYaxis()->FindBin(Pt);
	     //cout<<"biny = "<<biny<<endl;
	     //cout<<"ele scale factor crack2016 = "<<ele_scale_factors_gap2016->GetBinContent(binx,biny)<<endl;
	     if (ele_scale_factors_gap2016->GetBinContent(binx,biny)>0.) eff_weight_3*=ele_scale_factors_gap2016->GetBinContent(binx,biny);
	     //cout<<"eff_weight_3 = "<<eff_weight_3<<endl;
	   }
	   
	 }
       }
     }
     
     //cout<<"DATA_type = "<<DATA_type<<endl;
     
     if (DATA_type == "2016") eff_weight_3=1.; 
     
     // Changing the weight for pileup and LineShape and efficiency
     if (eff_weight_3>0.) {
       //cout<<"eff_weight_3 = "<<eff_weight_3<<endl;
       newweight_3=weight*pu_weight*eff_weight_3;}
     else newweight_3=weight*pu_weight;
     
     //cout << "Starting weight + pileup + efficiency= " << newweight << endl;
     // if(debug)
     //cout << "Efficiency Weight for the 2l: " << eff_weight_3 << " Final weight for ZZ pair = " << newweight << endl;
     
     
     //cut on Z1 mass (40 < MZ1 < 120)
     //--------------------------------
     
     
     massZ1 = Z1.massvalue;
     
     pxZ1 = Z1.pxZ;
     pyZ1 = Z1.pyZ;
     ptZ1 = sqrt( pxZ1*pxZ1 + pyZ1*pyZ1 );
     
     
     if ( Z1.massvalue > 40 && Z1.massvalue < 120 ) {
       
       if (Z1.tag ==1 ){//mumu
	 
	 hMZ1_3_mumu->Fill( massZ1,newweight_3 );
	 hPtZ1_3_mumu->Fill( ptZ1,newweight_3 );
    // hYZ1_3->Fill( Y_Z1,newweight );
       }
       else if (Z1.tag ==2 ){//ee
	 
	 hMZ1_3_ee->Fill( massZ1,newweight_3 );
	 hPtZ1_3_ee->Fill( ptZ1,newweight_3 );
	 // hYZ1_3->Fill( Y_Z1,newweight );
	 
       }
       
       hMZ1_3->Fill( massZ1,newweight_3 );
       hPtZ1_3->Fill( ptZ1,newweight_3 );
       
     }
     
     
     hPFMET_3->Fill(RECO_PFMET,newweight);
     
     
     if ( massZ1 <= 40 || massZ1 >= 120 ) continue;
     
//  if ( massZ1 <= 60 || massZ1 >= 120 ) continue; //not used it



     //Question from ARC

     //if (RECO_PFMET < 50 ) continue;
     
     ////////


     
     if ( Z1.tag ==1 ){ //Z1 is mumu
       
       //cout<<"Z1 is mumu "<<endl;
       
       //Mu FakeRate starts here
       //=======================//
       
       
       for (int i=0 ; i< N_loose_mu ; i++){
	 
	 if (N_loose_mu > 3 || N_loose_e> 0) continue;
	 
	 if (fabs( RECOMU_SIP[iL_loose_mu[i]] ) >= 4.)continue; //SIP cut for loose muon
	 
	 
	 if(iL_loose_mu[i] == Z1.ilept1 || iL_loose_mu[i]== Z1.ilept2 )continue; //exclude muons from Z1	  
	 
	 //cout<<"Muon Fake rate "<<endl;
	 
	 //cout<<"loose mu properties pt= "<<RECOMU_PT[ iL_loose_mu[i] ]<<"and energy = "<<RECOMU_E[ iL_loose_mu[i] ]<<endl;
	 
	 
	 // deltaR > 0.02 between any of 3 leptons
	 //=========================================
	 
	 
	 //cout<<"the final Z1 mass = "<<Z1.massvalue<<"with lepton indices = "<<Z1.ilept1<<","<<Z1.ilept2<<" and tag = "<<Z1.tag<<endl;
	 //cout<<"the third muon index = "<< iL_loose_mu[i]<<endl;
	 
	 //cout<<"Z1_lep1 Z1_lep2 deltaR = "<<sqrt( pow( DELTAPHI(  Z1.phi1 , Z1.phi2 ),2) + pow( Z1.eta1 - Z1.eta2 ,2) )<<endl;
	 //cout<<"Z1_lep1 third lep deltaR = "<<sqrt( pow( DELTAPHI(  Z1.phi1 , RECOMU_PHI[iL_loose_mu[i]] ),2) + pow( Z1.eta1 - RECOMU_ETA[iL_loose_mu[i]] ,2) )<<endl;
	 //cout<<"Z1_lep2 third lep deltaR = "<<sqrt( pow( DELTAPHI(  Z1.phi2 , RECOMU_PHI[iL_loose_mu[i]] ),2) + pow( Z1.eta2 - RECOMU_ETA[iL_loose_mu[i]] ,2) )<<endl;
	 
	 //Ghost cleaning
	 if ( sqrt( pow( DELTAPHI(  Z1.phi1 , Z1.phi2 ),2) + pow( Z1.eta1 - Z1.eta2 ,2) ) <= 0.02) continue; //Z1 (lep1,lep2)
	 if ( sqrt( pow( DELTAPHI(  Z1.phi1 , RECOMU_PHI[iL_loose_mu[i]] ),2) + pow( Z1.eta1 - RECOMU_ETA[iL_loose_mu[i]] ,2) ) <= 0.02) continue; //(Z1lep1 , loose mu)
	 if ( sqrt( pow( DELTAPHI(  Z1.phi2 , RECOMU_PHI[iL_loose_mu[i]] ),2) + pow( Z1.eta2 - RECOMU_ETA[iL_loose_mu[i]] ,2) ) <= 0.02) continue; //(Z1lep2 , loose mu)
	 
	 
	 ///////////////////////////////////////////////////////////////////////////////////////////////////////////
	 
	 
	 // pT >20,10 for any of three leptons
	 //=====================================
	 
	 vector<float> leptonspT;
	 leptonspT.clear();
	 
	 leptonspT.push_back(Z1.pt1);
	 leptonspT.push_back(Z1.pt2);
	 leptonspT.push_back(RECOMU_PT[ iL_loose_mu[i] ]);
	 
	 std::sort(leptonspT.rbegin(),leptonspT.rend());
	 
	 //cout<<"Z1pt1 = "<<Z1.pt1<<",Z1pt2 = "<<Z1.pt2<<" ,third muon pt  = "<<RECOMU_PT[ iL_loose_mu[i] ]<<endl;
	 
	 //cout<<"leptonspT = "<<leptonspT.at(0)<<" , "<<leptonspT.at(1)<<" , "<<leptonspT.at(2)<<endl;
	 
	 if ( leptonspT.at(0) > 20 && leptonspT.at(1) > 10 ){/*ok*/}
	 else continue;
	 
	 //////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	 
	 // QCD supression mll>4
	 //======================
	 
	 TLorentzVector Lepton1qcd,Lepton2qcd,Lepton3qcd,DiLeptonQCD;
	 
	 double ZQCD =-999;
	 double min_mass_2L =10000;
	 
	 //cout<<"Z1charge1 = "<< Z1.charge1<<" ,Z1charge2 = "<< Z1.charge2<<" ,loose mu charge = "<<RECOMU_CHARGE[iL_loose_mu[i]]<<endl;
	 
	 if (Z1.charge1 + RECOMU_CHARGE[iL_loose_mu[i]] == 0 ) //opposite charge
	   {
	     
	     Lepton1qcd.SetPtEtaPhiE(Z1.pt1, Z1.eta1, Z1.phi1, Z1.lep1E ); 
	     Lepton3qcd.SetPtEtaPhiE(RECOMU_PT[ iL_loose_mu[i] ], RECOMU_ETA[iL_loose_mu[i]], RECOMU_PHI[ iL_loose_mu[i] ], RECOMU_E[ iL_loose_mu[i] ] );
	     DiLeptonQCD=Lepton1qcd+Lepton3qcd;
	     
	     ZQCD = DiLeptonQCD.M();
	     
	   }
	 else {

	   Lepton2qcd.SetPtEtaPhiE(Z1.pt2, Z1.eta2, Z1.phi2, Z1.lep2E ); 
	   Lepton3qcd.SetPtEtaPhiE(RECOMU_PT[ iL_loose_mu[i] ], RECOMU_ETA[iL_loose_mu[i]], RECOMU_PHI[ iL_loose_mu[i] ], RECOMU_E[ iL_loose_mu[i] ] );
	   DiLeptonQCD=Lepton2qcd+Lepton3qcd;
	   
	      ZQCD = DiLeptonQCD.M();
	 }
	 
	 //cout<<"ZQCD = "<<ZQCD<<endl;
	 
	 if( Z1.massvalue < min_mass_2L ) min_mass_2L = Z1.massvalue ;
	 if( ZQCD < min_mass_2L ) min_mass_2L = ZQCD ;
	 
	 //cout<<"mmin_mass_2L = "<<min_mass_2L<<endl;
	 
	 if (min_mass_2L <=4)continue;
	 
	 
	 
	 /////////////////////////////////////////////////////////////////////////////////////////////////
	 
	 //Reduces contamination from photon conversion
	 if( fabs(massZ1 - Zmass) >= 7. ) continue;
	 
	 // Draw MZ1 and PFMET in CR Z+L 
	 //============================== 
	 
	 //Before cut (Den)
	 
	 if( fabs(RECOMU_ETA[iL_loose_mu[i]]) <= 1.2 ){
	   
	   hMZ1_loose_mu_barrel->Fill(massZ1,newweight_3);
	   hPFMET_loose_mu_barrel->Fill(RECO_PFMET,newweight);
	   hPT_loose_mu_barrel->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	   MET_PT_Mu_Barrel_Den->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);
	   MET_PT_Mu_Barrel_Den_Rebin->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);
	   //charge
	   if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==+1 ) {
	     h_PT_loose_mu_barrel_pos->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     h_PFMET_loose_mu_barrel_pos->Fill(RECO_PFMET,newweight);
	     h_MET_PT_loose_mu_barrel_pos->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight); }//positive
	   
	   else if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==-1 ){	     
	     h_PT_loose_mu_barrel_neg->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     h_PFMET_loose_mu_barrel_neg->Fill(RECO_PFMET,newweight);
	     h_MET_PT_loose_mu_barrel_neg->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);  }//negative 
	 }//Barrel
	 
	 else if( fabs(RECOMU_ETA[iL_loose_mu[i]]) > 1.2 ){
	   
	   hMZ1_loose_mu_endcap->Fill(massZ1,newweight_3);
	   hPT_loose_mu_endcap->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	   hPFMET_loose_mu_endcap->Fill(RECO_PFMET,newweight);
	   MET_PT_Mu_Endcap_Den->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);
	   MET_PT_Mu_Endcap_Den_Rebin->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);
	   //charge
	   if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==+1 ) {	     
	     h_PT_loose_mu_endcap_pos->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     h_PFMET_loose_mu_endcap_pos->Fill(RECO_PFMET,newweight);
	     h_MET_PT_loose_mu_endcap_pos->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);}//positive
	   
	   else if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==-1 ){	     
	     h_PT_loose_mu_endcap_neg->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     h_PFMET_loose_mu_endcap_neg->Fill(RECO_PFMET,newweight);
	     h_MET_PT_loose_mu_endcap_neg->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);}//negative
	 }//endcap    
	 
	 
	 //After cut (NUM)
	 
	 //if( RECOMU_isPFMu[iL_loose_mu[i]] && RECOMU_PFX_dB_new[iL_loose_mu[i]]<0.35){
	 
	 if( ( RECOMU_isPFMu[iL_loose_mu[i]] || (RECOMU_isTrackerHighPtMu[iL_loose_mu[i]] && RECOMU_PT[iL_loose_mu[i]] > 200.) ) && RECOMU_PFX_dB_new[iL_loose_mu[i]]<0.35  ){
	   
	   if( fabs(RECOMU_ETA[iL_loose_mu[i]]) <= 1.2 ) {
	     
	     hMZ1_tight_mu_barrel->Fill(massZ1,newweight_3);
	     hPT_tight_mu_barrel->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     hPFMET_tight_mu_barrel->Fill(RECO_PFMET,newweight);
	     MET_PT_Mu_Barrel_Num->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);
	     MET_PT_Mu_Barrel_Num_Rebin->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);
	     //charge
	     if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==+1 ) {
	     h_PT_tight_mu_barrel_pos->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     h_PFMET_tight_mu_barrel_pos->Fill(RECO_PFMET,newweight);
	     h_MET_PT_tight_mu_barrel_pos->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight); }//positive
	   
	   else if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==-1 ){	     
	     h_PT_tight_mu_barrel_neg->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     h_PFMET_tight_mu_barrel_neg->Fill(RECO_PFMET,newweight);
	     h_MET_PT_tight_mu_barrel_neg->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);  }//negative 
	   }//Barrel
	   
	   else if( fabs(RECOMU_ETA[iL_loose_mu[i]]) > 1.2 ){
	     
	     hMZ1_tight_mu_endcap->Fill(massZ1,newweight_3);
	     hPT_tight_mu_endcap->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     hPFMET_tight_mu_endcap->Fill(RECO_PFMET,newweight);
	     MET_PT_Mu_Endcap_Num->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);
	     MET_PT_Mu_Endcap_Num_Rebin->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);
	     //charge
	     if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==+1 ) {	     
	       h_PT_tight_mu_endcap_pos->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	       h_PFMET_tight_mu_endcap_pos->Fill(RECO_PFMET,newweight);
	       h_MET_PT_tight_mu_endcap_pos->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);}//positive
	     
	     else if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==-1 ){	     
	       h_PT_tight_mu_endcap_neg->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	       h_PFMET_tight_mu_endcap_neg->Fill(RECO_PFMET,newweight);
	       h_MET_PT_tight_mu_endcap_neg->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);}//negative
	     
	   }//endcap
	   
	 }
	 
	 
	 //////////////////////////////////////////////////////////////////////
	 
	 //use PFMET <25 or use PFMET <30
	 
	 
	 if(RECO_PFMET >= 25)continue; //this to supress prompet leptons from wz and ttbar
	 if( fabs(massZ1 - Zmass) >= 7. ) continue; //this suppress QCD effect
	 
	 //Before cut
	 //Info from the event
	 f_run   = Run;
	 f_lumi  = LumiSection;
	 f_event = Event;
	 f_weight = newweight;	 
	 f_zmumu_mu_pt_loose = RECOMU_PT[iL_loose_mu[i]];
	 f_zmumu_mu_eta_loose = RECOMU_ETA[iL_loose_mu[i]];
	 f_zmumu_mu_charge_loose = RECOMU_CHARGE[iL_loose_mu[i]];
	 
	 if(fabs(RECOMU_ETA[iL_loose_mu[i]]) <= 1.2 ){
	   ZplusM_Pt_DEN_Barrel->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	   if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==+1 ) {ZplusM_Pt_DEN_Barrel_pos->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);}
	    else if (  (RECOMU_CHARGE[iL_loose_mu[i]]) ==-1 ) {ZplusM_Pt_DEN_Barrel_neg->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);}  
	 }//Barrel
	 
	 
	 else if( fabs(RECOMU_ETA[iL_loose_mu[i]]) > 1.2 ){  
	   ZplusM_Pt_DEN_Endcaps->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	   if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==+1 ) {ZplusM_Pt_DEN_Endcaps_pos->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);}
	   else if (  (RECOMU_CHARGE[iL_loose_mu[i]]) ==-1 ) {ZplusM_Pt_DEN_Endcaps_neg->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);}
	 }//endcap
	 
	 
	 //After cut 
	 // if( RECOMU_isPFMu[iL_loose_mu[i]] && RECOMU_PFX_dB_new[iL_loose_mu[i]]<0.35  ){
	 
	 if( ( RECOMU_isPFMu[iL_loose_mu[i]] || (RECOMU_isTrackerHighPtMu[iL_loose_mu[i]] && RECOMU_PT[iL_loose_mu[i]] > 200.) ) && RECOMU_PFX_dB_new[iL_loose_mu[i]]<0.35  ){
	  //Info from the event
	  f_zmumu_mu_pt_tight = RECOMU_PT[iL_loose_mu[i]];
	  f_zmumu_mu_eta_tight = RECOMU_ETA[iL_loose_mu[i]];
	  f_zmumu_mu_charge_tight = RECOMU_CHARGE[iL_loose_mu[i]];	   
	   
	   if( fabs(RECOMU_ETA[iL_loose_mu[i]]) <= 1.2 ){
	     ZplusM_Pt_NUM_ID_ISO_Barrel->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==+1 ) {ZplusM_Pt_NUM_ID_ISO_Barrel_pos->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);}
	     else if (  (RECOMU_CHARGE[iL_loose_mu[i]]) ==-1 ) {ZplusM_Pt_NUM_ID_ISO_Barrel_neg->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);} 
	   } //Barrel 
	   else if(fabs (RECOMU_ETA[iL_loose_mu[i]]) > 1.2 ){
	     ZplusM_Pt_NUM_ID_ISO_Endcaps->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==+1 ) {ZplusM_Pt_NUM_ID_ISO_Endcaps_pos->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);}
	     else if (  (RECOMU_CHARGE[iL_loose_mu[i]]) ==-1 ) {ZplusM_Pt_NUM_ID_ISO_Endcaps_neg->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);}
	   } //endcap
	 }//tight id, iso
	 
	 //cout<<"histos filled for Mu FR"<<endl;
	 
	 //sprintf (Eventformat,"FR=%d:%d:%d",Run,LumiSection,Event);
	 
       }//end loop on loose muon  for (int i=0 ; i< N_loose ; i++)
       
       //end of Muon Fake rate

              
       //ele FakeRate starts here
       //=====================================================//
       
       for(int i = 0; i < N_loose_e ; i++){
	 
	 if (N_loose_mu > 2 || N_loose_e> 1) continue;
	 
	 if (fabs(RECOELE_SIP[iL_loose_e[i]])>=4.)continue; //SIP cut
	 
	 //cout<<"electron fake rate"<<endl;
	 

	 ///////////////////////////////////////////////////////////////////////////////////////////////////////////
	 
	 //deltaR > 0.02 between any of 3 leptons
	 //========================================
	 
	 
	 //cout<<"the final Z1 mass = "<<Z1.massvalue<<"with lepton indices = "<<Z1.ilept1<<","<<Z1.ilept2<<" and tag = "<<Z1.tag<<endl;
	 //cout<<"the third muon index = "<< iL_loose_e[i]<<endl;
	 
	 //cout<<"Z1_lep1 Z1_lep2 deltaR = "<<sqrt( pow( DELTAPHI(  Z1.phi1 , Z1.phi2 ),2) + pow( Z1.eta1 - Z1.eta2 ,2) )<<endl;
	 //cout<<"Z1_lep1 third lep deltaR = "<<sqrt( pow( DELTAPHI(  Z1.phi1 , RECOELE_PHI[iL_loose_e[i]] ),2) + pow( Z1.eta1 - RECOELE_ETA[iL_loose_e[i]] ,2) )<<endl;
	 //cout<<"Z1_lep2 third lep deltaR = "<<sqrt( pow( DELTAPHI(  Z1.phi2 , RECOELE_PHI[iL_loose_e[i]] ),2) + pow( Z1.eta2 - RECOELE_ETA[iL_loose_e[i]] ,2) )<<endl;
	 
	 if ( sqrt( pow( DELTAPHI(  Z1.phi1 , Z1.phi2 ),2) + pow( Z1.eta1 - Z1.eta2 ,2) ) <= 0.02) continue; //Z1 (lep1,lep2)
	 if ( sqrt( pow( DELTAPHI(  Z1.phi1 , RECOELE_PHI[iL_loose_e[i]] ),2) + pow( Z1.eta1 - RECOELE_ETA[iL_loose_e[i]] ,2) ) <= 0.02) continue; //(Z1lep1 , loose mu)
	 if ( sqrt( pow( DELTAPHI(  Z1.phi2 , RECOELE_PHI[iL_loose_e[i]] ),2) + pow( Z1.eta2 - RECOELE_ETA[iL_loose_e[i]] ,2) ) <= 0.02) continue; //(Z1lep2 , loose mu)
	 

	 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 
	 // pT >20,10 for any of three leptons
	 //=====================================
	 
	 vector<float> leptonspT;
	 leptonspT.clear();
	 
	 leptonspT.push_back(Z1.pt1);
	 leptonspT.push_back(Z1.pt2);
	 leptonspT.push_back(RECOELE_PT[ iL_loose_e[i] ]);
	 
	 std::sort(leptonspT.rbegin(),leptonspT.rend());
	 
	 //cout<<"Z1pt1 = "<<Z1.pt1<<",Z1pt2 = "<<Z1.pt2<<" ,third ele pt  = "<<RECOELE_PT[ iL_loose_e[i] ]<<endl;
	 
	 //cout<<"leptonspT = "<<leptonspT.at(0)<<" , "<<leptonspT.at(1)<<" , "<<leptonspT.at(2)<<endl;
	 
	 if ( leptonspT.at(0) > 20 && leptonspT.at(1) > 10 ){/*ok*/}
	 else continue;
	 
	 /////////////////////////////////////////////////////////////////////////////////////////////////////////////

	 
	 // QCD supression mll>4
	 //======================
	 
	 TLorentzVector Lepton1qcd,Lepton2qcd,Lepton3qcd,DiLeptonQCD;
	 
	 double ZQCD =-999;
	 double min_mass_2L =10000;
	 
	 //cout<<"Z1charge1 = "<< Z1.charge1<<" ,Z1charge2 = "<< Z1.charge2<<" ,loose ele charge = "<<RECOELE_CHARGE[iL_loose_e[i]]<<endl;
	 
	 if (Z1.charge1 + RECOELE_CHARGE[iL_loose_e[i]] == 0 ) //opposite charge
	   {
	     
	     Lepton1qcd.SetPtEtaPhiE(Z1.pt1, Z1.eta1, Z1.phi1, Z1.lep1E ); 
	     Lepton3qcd.SetPtEtaPhiE(RECOELE_PT[ iL_loose_e[i] ], RECOELE_ETA[iL_loose_e[i]], RECOELE_PHI[ iL_loose_e[i] ], RECOELE_E[ iL_loose_e[i] ] );
	     DiLeptonQCD=Lepton1qcd+Lepton3qcd;
	     
	     ZQCD = DiLeptonQCD.M();
	     
	   }
	 else {
	   
	   Lepton2qcd.SetPtEtaPhiE(Z1.pt2, Z1.eta2, Z1.phi2, Z1.lep2E ); 
	   Lepton3qcd.SetPtEtaPhiE(RECOELE_PT[ iL_loose_e[i] ], RECOELE_ETA[iL_loose_e[i]], RECOELE_PHI[ iL_loose_e[i] ], RECOELE_E[ iL_loose_e[i] ] );
	   DiLeptonQCD=Lepton2qcd+Lepton3qcd;
	   
	   ZQCD = DiLeptonQCD.M();
	 }
	 
	 
	 if( Z1.massvalue < min_mass_2L ) min_mass_2L = Z1.massvalue ;
	 if( ZQCD < min_mass_2L ) min_mass_2L = ZQCD ;
	 
	 //cout<<"mmin_mass_2L = "<<min_mass_2L<<endl;
	 
	 if (min_mass_2L <=4)continue;
	 
	 ///////////////////////////////////////////////////////////////////////////////////////////////
	 //good electron id
	 
	 bool BDT_ok = 0; // Fall15 with CMSSW_7_6_x
	 
	 if( RECOELE_PT[i] > 7. &&  RECOELE_PT[i] <= 10. ){
	   if( fabs(RECOELE_scl_Eta[i]) < .8 && RECOELE_mvaNonTrigV0[i] > -0.211 ) BDT_ok = 1 ;
	   if( ( fabs(RECOELE_scl_Eta[i]) >= .8 && fabs(RECOELE_scl_Eta[i]) < 1.479 )
	       && RECOELE_mvaNonTrigV0[i] > -0.396 ) BDT_ok = 1 ;
	   if( fabs(RECOELE_scl_Eta[i]) >= 1.479 && RECOELE_mvaNonTrigV0[i] > -0.215 ) BDT_ok = 1 ;
	 }
	 else { 
	   if( fabs(RECOELE_scl_Eta[i]) < .8 && RECOELE_mvaNonTrigV0[i] > -0.870 ) BDT_ok = 1 ;
	   if( ( fabs(RECOELE_scl_Eta[i]) >= .8 && fabs(RECOELE_scl_Eta[i]) <= 1.479 )
	       && RECOELE_mvaNonTrigV0[i] > -0.838 ) BDT_ok = 1 ;
	   if( fabs(RECOELE_scl_Eta[i]) > 1.479 && RECOELE_mvaNonTrigV0[i] > -0.763 ) BDT_ok = 1 ;
	 }
	 
	 /////////////////////////////////////////////////////////////////////////////////////////

	 if( massZ1 <= (Zmass-7.) || massZ1 >= (Zmass+7.) ) continue; // |MZ1 - Zmass| <10
	 
	 //Draw MZ1 and PFMET in CR Z1+L
	 //=============================
	 
	 //Before cut (DEN)
	 
	 if( fabs(RECOELE_ETA[iL_loose_e[i]]) <= 1.479 ) {
	   
	   hMZ1_loose_ele_barrel->Fill(massZ1,newweight_3);
	   hPFMET_loose_ele_barrel->Fill(RECO_PFMET,newweight);
	   hPT_loose_ele_barrel->Fill(RECOELE_PT[iL_loose_e[i]],newweight); 
	   MET_PT_Ele_Barrel_Den->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);
	   MET_PT_Ele_Barrel_Den_Rebin->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);
	   //charge
	   if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==+1 ) {
	     h_PT_loose_ele_barrel_pos->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	     h_PFMET_loose_ele_barrel_pos->Fill(RECO_PFMET,newweight);
	     h_MET_PT_loose_ele_barrel_pos->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight); }//positive
	   
	   else if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==-1 ){	     
	     h_PT_loose_ele_barrel_neg->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	     h_PFMET_loose_ele_barrel_neg->Fill(RECO_PFMET,newweight);
	     h_MET_PT_loose_ele_barrel_neg->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);  }//negative 
	 }//Barrel
	 
	 
	 else if( fabs(RECOELE_ETA[iL_loose_e[i]]) > 1.479 )  {

	   hMZ1_loose_ele_endcap->Fill(massZ1,newweight_3);
	   hPFMET_loose_ele_endcap->Fill(RECO_PFMET,newweight);
	   hPT_loose_ele_endcap->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	   MET_PT_Ele_Endcap_Den->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);
	   MET_PT_Ele_Endcap_Den_Rebin->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);
	   //charge
	   if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==+1 ) {
	     h_PT_loose_ele_endcap_pos->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	     h_PFMET_loose_ele_endcap_pos->Fill(RECO_PFMET,newweight);
	     h_MET_PT_loose_ele_endcap_pos->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight); }//positive
	   
	   else if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==-1 ){	     
	     h_PT_loose_ele_endcap_neg->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	     h_PFMET_loose_ele_endcap_neg->Fill(RECO_PFMET,newweight);
	     h_MET_PT_loose_ele_endcap_neg->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);  }//negative	   
	 }//endcap    
	 
	 
	 //After cut (NUM)
	 
	 if( BDT_ok && RECOELE_PFX_rho_new[iL_loose_e[i]]<0.35 ){
	   
	   //Barrel
	   
	   if( fabs(RECOELE_ETA[iL_loose_e[i]]) <= 1.479 ){
	     
	     hMZ1_tight_ele_barrel->Fill(massZ1,newweight_3); 
	     hPFMET_tight_ele_barrel->Fill(RECO_PFMET,newweight);
	     hPT_tight_ele_barrel->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	     MET_PT_Ele_Barrel_Num->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);
	     MET_PT_Ele_Barrel_Num_Rebin->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);
	   //charge
	   if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==+1 ) {
	     h_PT_tight_ele_barrel_pos->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	     h_PFMET_tight_ele_barrel_pos->Fill(RECO_PFMET,newweight);
	     h_MET_PT_tight_ele_barrel_pos->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight); }//positive

	   else if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==-1 ){	     
	     h_PT_tight_ele_barrel_neg->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	     h_PFMET_tight_ele_barrel_neg->Fill(RECO_PFMET,newweight);
	     h_MET_PT_tight_ele_barrel_neg->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);  }//negative 	   
	   }//Barrel
	   
	   //Endcap
	   
	   else if( fabs(RECOELE_ETA[iL_loose_e[i]]) > 1.479 ) {
	     
	     hMZ1_tight_ele_endcap->Fill(massZ1,newweight_3); 	     
	     hPFMET_tight_ele_endcap->Fill(RECO_PFMET,newweight);
	     hPT_tight_ele_endcap->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	     MET_PT_Ele_Endcap_Num->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);
	     MET_PT_Ele_Endcap_Num_Rebin->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);
	     	   //charge
	   if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==+1 ) {
	     h_PT_tight_ele_endcap_pos->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	     h_PFMET_tight_ele_endcap_pos->Fill(RECO_PFMET,newweight);
	     h_MET_PT_tight_ele_endcap_pos->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight); }//positive

	   else if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==-1 ){	     
	     h_PT_tight_ele_endcap_neg->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	     h_PFMET_tight_ele_endcap_neg->Fill(RECO_PFMET,newweight);
	     h_MET_PT_tight_ele_endcap_neg->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);  }//negative 
	   }//endcap
	   
	 }//tight ID, iso
	 
	 //////////////////////////////////////////////////////////////////////////////////////
	 
	 //use PFMET <25 or use PFMET <30
	 
	 if (RECO_PFMET >= 25)continue; //this to supress prompet leptons from wz and ttbar
	 if( fabs(massZ1 - Zmass) >= 7. ) continue; //this suppress QCD effect	 
	 
	  //Info from the event
	  f_run   = Run;
	  f_lumi  = LumiSection;
	  f_event = Event;
	  f_weight = newweight;	   
	  f_zmumu_e_pt_loose = RECOELE_PT[iL_loose_e[i]];
	  f_zmumu_e_eta_loose = RECOELE_ETA[iL_loose_e[i]];
	  f_zmumu_e_charge_loose = RECOELE_CHARGE[iL_loose_e[i]];
	 
	 if( fabs(RECOELE_ETA[iL_loose_e[i]]) <= 1.479 ){
	   ZplusE_Pt_DEN_Barrel->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	   if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==+1 ) {ZplusE_Pt_DEN_Barrel_pos->Fill(RECOELE_PT[iL_loose_e[i]],newweight);}
	   else if (  (RECOELE_CHARGE[iL_loose_e[i]]) ==-1 ) {ZplusE_Pt_DEN_Barrel_neg->Fill(RECOELE_PT[iL_loose_e[i]],newweight);}	 
	 }//Barrel
	 
	 else if( fabs(RECOELE_ETA[iL_loose_e[i]]) > 1.479 ){
	   ZplusE_Pt_DEN_Endcaps->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	   if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==+1 ) {ZplusE_Pt_DEN_Endcaps_pos->Fill(RECOELE_PT[iL_loose_e[i]],newweight);}
	   else if (  (RECOELE_CHARGE[iL_loose_e[i]]) ==-1 ) {ZplusE_Pt_DEN_Endcaps_neg->Fill(RECOELE_PT[iL_loose_e[i]],newweight);}  
	   
	 }//endcap
	 
	 
	 //After cut 
	 
	 if( BDT_ok && RECOELE_PFX_rho_new[iL_loose_e[i]]<0.35 ){
	  f_zmumu_e_pt_tight = RECOELE_PT[iL_loose_e[i]];
	  f_zmumu_e_eta_tight = RECOELE_ETA[iL_loose_e[i]];
	  f_zmumu_e_charge_tight = RECOELE_CHARGE[iL_loose_e[i]];
	   
	   //Barrel
	   
	   if( fabs(RECOELE_ETA[iL_loose_e[i]]) <= 1.479 ){
	     ZplusE_Pt_NUM_ID_ISO_Barrel->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	     if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==+1 ) {ZplusE_Pt_NUM_ID_ISO_Barrel_pos->Fill(RECOELE_PT[iL_loose_e[i]],newweight);}
	     else if (  (RECOELE_CHARGE[iL_loose_e[i]]) ==-1 ) {ZplusE_Pt_NUM_ID_ISO_Barrel_neg->Fill(RECOELE_PT[iL_loose_e[i]],newweight);}	 
	   }//Barrel
	   
	   //Endcap
	   
	   else if( fabs(RECOELE_ETA[iL_loose_e[i]]) > 1.479 ){
	     ZplusE_Pt_NUM_ID_ISO_Endcaps->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	     if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==+1 ) {ZplusE_Pt_NUM_ID_ISO_Endcaps_pos->Fill(RECOELE_PT[iL_loose_e[i]],newweight);}
	     else if (  (RECOELE_CHARGE[iL_loose_e[i]]) ==-1 ) {ZplusE_Pt_NUM_ID_ISO_Endcaps_neg->Fill(RECOELE_PT[iL_loose_e[i]],newweight);}	     
	   }	   
	 }//Tight id, iso
	 
	 //sprintf (Eventformat,"FR=%d:%d:%d",Run,LumiSection,Event); 
  
       }//end loop on loose electrons for(int i = 0; i < N_loose_e ; i++)
       
       //end of electron Fake rate
       
       // R_3_b++;
       //R_3_b_w=R_3_b_w+newweight;
       
       
     }//end if ( Z1tag ==1 )
     
     else if ( Z1.tag ==2 ){ //Z1 is ee
       
       //cout<<"Z1 is ee "<<endl;

       //Mu FakeRate starts here
       //==========================//
       
       for (int i=0 ; i< N_loose_mu ; i++){
	 
	 if(N_loose_mu > 1 || N_loose_e > 2)continue;
	 
	 
	 if (fabs( RECOMU_SIP[iL_loose_mu[i]] ) >= 4.)continue; //SIP cut for loose muon
	 
	 //cout<<"Muon Fake rate "<<endl;
	 
	 //cout<<"loose mu properties pt= "<<RECOMU_PT[ iL_loose_mu[i] ]<<"and energy = "<<RECOMU_E[ iL_loose_mu[i] ]<<endl;

	 
	 //////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 
	 // deltaR > 0.02 between any of 3 leptons
	 //=========================================
	 
	 
	 //cout<<"the final Z1 mass = "<<Z1.massvalue<<"with lepton indices = "<<Z1.ilept1<<","<<Z1.ilept2<<" and tag = "<<Z1.tag<<endl;
	 //cout<<"the third muon index = "<< iL_loose_mu[i]<<endl;
	 
	 //cout<<"Z1_lep1 Z1_lep2 deltaR = "<<sqrt( pow( DELTAPHI(  Z1.phi1 , Z1.phi2 ),2) + pow( Z1.eta1 - Z1.eta2 ,2) )<<endl;
	 //cout<<"Z1_lep1 third lep deltaR = "<<sqrt( pow( DELTAPHI(  Z1.phi1 , RECOMU_PHI[iL_loose_mu[i]] ),2) + pow( Z1.eta1 - RECOMU_ETA[iL_loose_mu[i]] ,2) )<<endl;
	 //cout<<"Z1_lep2 third lep deltaR = "<<sqrt( pow( DELTAPHI(  Z1.phi2 , RECOMU_PHI[iL_loose_mu[i]] ),2) + pow( Z1.eta2 - RECOMU_ETA[iL_loose_mu[i]] ,2) )<<endl;
	 
	 if ( sqrt( pow( DELTAPHI(  Z1.phi1 , Z1.phi2 ),2) + pow( Z1.eta1 - Z1.eta2 ,2) ) <= 0.02) continue; //Z1 (lep1,lep2)
	 if ( sqrt( pow( DELTAPHI(  Z1.phi1 , RECOMU_PHI[iL_loose_mu[i]] ),2) + pow( Z1.eta1 - RECOMU_ETA[iL_loose_mu[i]] ,2) ) <= 0.02) continue; //(Z1lep1 , loose mu)
	 if ( sqrt( pow( DELTAPHI(  Z1.phi2 , RECOMU_PHI[iL_loose_mu[i]] ),2) + pow( Z1.eta2 - RECOMU_ETA[iL_loose_mu[i]] ,2) ) <= 0.02) continue; //(Z1lep2 , loose mu)
	 
	 
	 ///////////////////////////////////////////////////////////////////////////////////////////////////////////
	 
	 
	 // pT >20,10 for any of three leptons
	 //=====================================
	 
	 vector<float> leptonspT;
	 leptonspT.clear();
	 
	 leptonspT.push_back(Z1.pt1);
	 leptonspT.push_back(Z1.pt2);
	 leptonspT.push_back(RECOMU_PT[ iL_loose_mu[i] ]);
	 
	 std::sort(leptonspT.rbegin(),leptonspT.rend());
	 
	 //cout<<"Z1pt1 = "<<Z1.pt1<<",Z1pt2 = "<<Z1.pt2<<" ,third muon pt  = "<<RECOMU_PT[ iL_loose_mu[i] ]<<endl;
	 
	 //cout<<"leptonspT = "<<leptonspT.at(0)<<" , "<<leptonspT.at(1)<<" , "<<leptonspT.at(2)<<endl;
	 
	 if ( leptonspT.at(0) > 20 && leptonspT.at(1) > 10 ){/*ok*/}
	 else continue;
	 

	 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 
	 // QCD supression mll>4
	 //======================
	 
	 TLorentzVector Lepton1qcd,Lepton2qcd,Lepton3qcd,DiLeptonQCD;
	 
	 double ZQCD =-999;
	 double min_mass_2L =10000;
	 
	 //cout<<"Z1charge1 = "<< Z1.charge1<<" ,Z1charge2 = "<< Z1.charge2<<" ,loose mu charge = "<<RECOMU_CHARGE[iL_loose_mu[i]]<<endl;
	 
	 if (Z1.charge1 + RECOMU_CHARGE[iL_loose_mu[i]] == 0 ) //opposite charge
	   {
	     
	     Lepton1qcd.SetPtEtaPhiE(Z1.pt1, Z1.eta1, Z1.phi1, Z1.lep1E ); 
	     Lepton3qcd.SetPtEtaPhiE(RECOMU_PT[ iL_loose_mu[i] ], RECOMU_ETA[iL_loose_mu[i]], RECOMU_PHI[ iL_loose_mu[i] ], RECOMU_E[ iL_loose_mu[i] ] );
	     DiLeptonQCD=Lepton1qcd+Lepton3qcd;
	     
	     ZQCD = DiLeptonQCD.M();
	     
	   }
	 else {
	   
	   Lepton2qcd.SetPtEtaPhiE(Z1.pt2, Z1.eta2, Z1.phi2, Z1.lep2E ); 
	   Lepton3qcd.SetPtEtaPhiE(RECOMU_PT[ iL_loose_mu[i] ], RECOMU_ETA[iL_loose_mu[i]], RECOMU_PHI[ iL_loose_mu[i] ], RECOMU_E[ iL_loose_mu[i] ] );
	   DiLeptonQCD=Lepton2qcd+Lepton3qcd;
	   
	   ZQCD = DiLeptonQCD.M();
	 }
	 
	 //cout<<"ZQCD = "<<ZQCD<<endl;
	 
	 if( Z1.massvalue < min_mass_2L ) min_mass_2L = Z1.massvalue ;
	 if( ZQCD < min_mass_2L ) min_mass_2L = ZQCD ;
	 
	 //cout<<"mmin_mass_2L = "<<min_mass_2L<<endl;
	 
	 if (min_mass_2L <=4)continue;
	 
	 
	 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	 if( massZ1 <= (Zmass-7.) || massZ1 >= (Zmass+7.) ) continue; 
	 
	 // Draw MZ1 and PFMET in CR Z+L 
	 //==============================
	 
	 //Before cut (Den)
	 
	 if( fabs(RECOMU_ETA[iL_loose_mu[i]]) <= 1.2 )  {
	   hMZ1_loose_mu_barrel->Fill(massZ1,newweight_3);
	   hPFMET_loose_mu_barrel->Fill(RECO_PFMET,newweight);
	   hPT_loose_mu_barrel->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	   MET_PT_Mu_Barrel_Den->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);
	   MET_PT_Mu_Barrel_Den_Rebin->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);
	   //charge
	   if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==+1 ) {
	     h_PT_loose_mu_barrel_pos->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     h_PFMET_loose_mu_barrel_pos->Fill(RECO_PFMET,newweight);
	     h_MET_PT_loose_mu_barrel_pos->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight); }//positive
	   
	   else if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==-1 ){	     
	     h_PT_loose_mu_barrel_neg->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     h_PFMET_loose_mu_barrel_neg->Fill(RECO_PFMET,newweight);
	     h_MET_PT_loose_mu_barrel_neg->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);  }//negative
	 }//Barrel
      
	 else if( fabs(RECOMU_ETA[iL_loose_mu[i]]) > 1.2 )  {
	   hMZ1_loose_mu_endcap->Fill(massZ1,newweight_3);
	   hPFMET_loose_mu_endcap->Fill(RECO_PFMET,newweight);
	   hPT_loose_mu_endcap->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	   MET_PT_Mu_Endcap_Den->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);
	   MET_PT_Mu_Endcap_Den_Rebin->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);
	   //charge
	   if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==+1 ) {	     
	     h_PT_loose_mu_endcap_pos->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     h_PFMET_loose_mu_endcap_pos->Fill(RECO_PFMET,newweight);
	     h_MET_PT_loose_mu_endcap_pos->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);}//positive
	   
	   else if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==-1 ){	     
	     h_PT_loose_mu_endcap_neg->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     h_PFMET_loose_mu_endcap_neg->Fill(RECO_PFMET,newweight);
	     h_MET_PT_loose_mu_endcap_neg->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);}//negative
	 }//endcap    
	 
	 
	 //After cut (NUM)
	 
	 // if( RECOMU_isPFMu[iL_loose_mu[i]] && RECOMU_PFX_dB_new[iL_loose_mu[i]]<0.35 ){
	 
	 if( ( RECOMU_isPFMu[iL_loose_mu[i]] || (RECOMU_isTrackerHighPtMu[iL_loose_mu[i]] && RECOMU_PT[iL_loose_mu[i]] > 200.) ) && RECOMU_PFX_dB_new[iL_loose_mu[i]]<0.35  ){
	   
	   if( fabs(RECOMU_ETA[iL_loose_mu[i]]) <= 1.2 ) {
	     hMZ1_tight_mu_barrel->Fill(massZ1,newweight_3); 
	     hPFMET_tight_mu_barrel->Fill(RECO_PFMET,newweight);
	     hPT_tight_mu_barrel->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     MET_PT_Mu_Barrel_Num->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);
	     MET_PT_Mu_Barrel_Num_Rebin->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);
	   //charge
	   if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==+1 ) {
	     h_PT_tight_mu_barrel_pos->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     h_PFMET_tight_mu_barrel_pos->Fill(RECO_PFMET,newweight);
	     h_MET_PT_tight_mu_barrel_pos->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight); }//positive
	   
	   else if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==-1 ){	     
	     h_PT_tight_mu_barrel_neg->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     h_PFMET_tight_mu_barrel_neg->Fill(RECO_PFMET,newweight);
	     h_MET_PT_tight_mu_barrel_neg->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);  }//negative
	   }//Barrel
	   
	   else if( fabs(RECOMU_ETA[iL_loose_mu[i]]) > 1.2 ){
	     hMZ1_tight_mu_endcap->Fill(massZ1,newweight_3); 
	     hPFMET_tight_mu_endcap->Fill(RECO_PFMET,newweight);
	     hPT_tight_mu_endcap->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     MET_PT_Mu_Endcap_Num->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);
	     MET_PT_Mu_Endcap_Num_Rebin->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight);  
	     //charge
	     if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==+1 ) {
	       h_PT_tight_mu_endcap_pos->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	       h_PFMET_tight_mu_endcap_pos->Fill(RECO_PFMET,newweight);
	       h_MET_PT_tight_mu_endcap_pos->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight); }//positive
	     
	     else if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==-1 ){	     
	       h_PT_tight_mu_endcap_neg->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	       h_PFMET_tight_mu_endcap_neg->Fill(RECO_PFMET,newweight);
	       h_MET_PT_tight_mu_endcap_neg->Fill(RECO_PFMET,RECOMU_PT[iL_loose_mu[i]],newweight); }//negative	     
	   }//endcap
	   
	 }//tight ID , iso
	 
	 
	 
	 /////////////////////////////////////////////////////////////////////////////////////////
	 
	 //use PFMET <25 or use PFMET <30
	 
	 if(RECO_PFMET>=25) continue ;
	 if( fabs(massZ1 - Zmass) >= 7. ) continue; //this suppress QCD effect
	 
	  //Info from the event
	  f_run   = Run;
	  f_lumi  = LumiSection;
	  f_event = Event;
	  f_weight = newweight;	   
	  f_zee_mu_pt_loose = RECOMU_PT[iL_loose_mu[i]];
	  f_zee_mu_eta_loose = RECOMU_ETA[iL_loose_mu[i]];
	  f_zee_mu_charge_loose = RECOMU_CHARGE[iL_loose_mu[i]];
	 
	 if( fabs(RECOMU_ETA[iL_loose_mu[i]]) <= 1.2 ){
	   ZplusM_Pt_DEN_Barrel->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	   if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==+1 ) {ZplusM_Pt_DEN_Barrel_pos->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);}
	   else if (  (RECOMU_CHARGE[iL_loose_mu[i]]) ==-1 ) {ZplusM_Pt_DEN_Barrel_neg->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);} 
	 }//Barrel
	 else if( fabs(RECOMU_ETA[iL_loose_mu[i]]) > 1.2 ){
	   ZplusM_Pt_DEN_Endcaps->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	   if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==+1 ) {ZplusM_Pt_DEN_Endcaps_pos->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);}
	   else if (  (RECOMU_CHARGE[iL_loose_mu[i]]) ==-1 ) {ZplusM_Pt_DEN_Endcaps_neg->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);} 
	 }//endcap
	 
	 
	 // if( RECOMU_isPFMu[iL_loose_mu[i]] && RECOMU_PFX_dB_new[iL_loose_mu[i]]<0.35 ){
	 
	 if( ( RECOMU_isPFMu[iL_loose_mu[i]] || (RECOMU_isTrackerHighPtMu[iL_loose_mu[i]] && RECOMU_PT[iL_loose_mu[i]] > 200.) ) && RECOMU_PFX_dB_new[iL_loose_mu[i]]<0.35  ){
	  f_zee_mu_pt_tight = RECOMU_PT[iL_loose_mu[i]];
	  f_zee_mu_eta_tight = RECOMU_ETA[iL_loose_mu[i]];
	  f_zee_mu_charge_tight = RECOMU_CHARGE[iL_loose_mu[i]];
	   
	   if( fabs(RECOMU_ETA[iL_loose_mu[i]]) <= 1.2 ){
	     ZplusM_Pt_NUM_ID_ISO_Barrel->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==+1 ) {ZplusM_Pt_NUM_ID_ISO_Barrel_pos->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);}
	     else if (  (RECOMU_CHARGE[iL_loose_mu[i]]) ==-1 ) {ZplusM_Pt_NUM_ID_ISO_Barrel_neg->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);} 
	   }//Barrel
	   
	   else if( fabs(RECOMU_ETA[iL_loose_mu[i]]) > 1.2 ){
	     ZplusM_Pt_NUM_ID_ISO_Endcaps->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);
	     if ( (RECOMU_CHARGE[iL_loose_mu[i]]) ==+1 ) {ZplusM_Pt_NUM_ID_ISO_Endcaps_pos->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);}
	     else if (  (RECOMU_CHARGE[iL_loose_mu[i]]) ==-1 ) {ZplusM_Pt_NUM_ID_ISO_Endcaps_neg->Fill(RECOMU_PT[iL_loose_mu[i]],newweight);} 
	   }//endcap
	   
	 }//tight id, iso
	 
	 //cout<<"histos filled for Mu FR"<<endl;
	 
	 //sprintf (Eventformat,"FR=%d:%d:%d",Run,LumiSection,Event);
	 
       }//end for loop loose muons  for (int i=0 ; i< N_loose ; i++)
       //Mu FakeRate ends here
       
       
       //ele FakeRate starts here
       //============================//
       
       for(int i = 0; i < N_loose_e ; i++){
	 
	 if(N_loose_mu > 0 || N_loose_e > 3)continue;
	 
	 if (fabs(RECOELE_SIP[iL_loose_e[i]])>=4.)continue; //SIP cut
	 
	 // if( fabs(RECOELE_PT[iL_loose_e[i]]-RECOELE_PT[indexlep1Z1])<0.0001 || fabs(RECOELE_PT[iL_loose_e[i]]-RECOELE_PT[indexlep2Z1])<0.0001 ) continue ;
	 
	 if(iL_loose_e[i] == Z1.ilept1 || iL_loose_e[i]== Z1.ilept2 )continue;	 
	 
	 //cout<<"electron fake rate "<<endl;
	 

	 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 
	 //deltaR > 0.02 between any of 3 leptons
	 //========================================
	 
	 
	 //cout<<"the final Z1 mass = "<<Z1.massvalue<<"with lepton indices = "<<Z1.ilept1<<","<<Z1.ilept2<<" and tag = "<<Z1.tag<<endl;
	 //cout<<"the third muon index = "<< iL_loose_e[i]<<endl;
	 
	 //cout<<"Z1_lep1 Z1_lep2 deltaR = "<<sqrt( pow( DELTAPHI(  Z1.phi1 , Z1.phi2 ),2) + pow( Z1.eta1 - Z1.eta2 ,2) )<<endl;
	 //cout<<"Z1_lep1 third lep deltaR = "<<sqrt( pow( DELTAPHI(  Z1.phi1 , RECOELE_PHI[iL_loose_e[i]] ),2) + pow( Z1.eta1 - RECOELE_ETA[iL_loose_e[i]] ,2) )<<endl;
	 //cout<<"Z1_lep2 third lep deltaR = "<<sqrt( pow( DELTAPHI(  Z1.phi2 , RECOELE_PHI[iL_loose_e[i]] ),2) + pow( Z1.eta2 - RECOELE_ETA[iL_loose_e[i]] ,2) )<<endl;
	 
	 if ( sqrt( pow( DELTAPHI(  Z1.phi1 , Z1.phi2 ),2) + pow( Z1.eta1 - Z1.eta2 ,2) ) <= 0.02) continue; //Z1 (lep1,lep2)
	 if ( sqrt( pow( DELTAPHI(  Z1.phi1 , RECOELE_PHI[iL_loose_e[i]] ),2) + pow( Z1.eta1 - RECOELE_ETA[iL_loose_e[i]] ,2) ) <= 0.02) continue; //(Z1lep1 , loose mu)
	 if ( sqrt( pow( DELTAPHI(  Z1.phi2 , RECOELE_PHI[iL_loose_e[i]] ),2) + pow( Z1.eta2 - RECOELE_ETA[iL_loose_e[i]] ,2) ) <= 0.02) continue; //(Z1lep2 , loose mu)
	 
	 
	 /////////////////////////////////////////////////////////////////////////////////////////////
	 
	 // pT >20,10 for any of three leptons
	 //=====================================
	 
	 vector<float> leptonspT;
	 leptonspT.clear();
	 
	 leptonspT.push_back(Z1.pt1);
	 leptonspT.push_back(Z1.pt2);
	 leptonspT.push_back(RECOELE_PT[ iL_loose_e[i] ]);
	 
	 std::sort(leptonspT.rbegin(),leptonspT.rend());
	 
	 //cout<<"Z1pt1 = "<<Z1.pt1<<",Z1pt2 = "<<Z1.pt2<<" ,third ele pt  = "<<RECOELE_PT[ iL_loose_e[i] ]<<endl;
	 
	 //cout<<"leptonspT = "<<leptonspT.at(0)<<" , "<<leptonspT.at(1)<<" , "<<leptonspT.at(2)<<endl;
	 
	 if ( leptonspT.at(0) > 20 && leptonspT.at(1) > 10 ){/*ok*/}
	 else continue;
	 
	 
	 ///////////////////////////////////////////////////////////////////////////////////////////
	 
	 // QCD supression mll>4
	 //======================
	 
	 TLorentzVector Lepton1qcd,Lepton2qcd,Lepton3qcd,DiLeptonQCD;
	 
	 double ZQCD =-999;
	 double min_mass_2L =10000;
	 
	 //cout<<"Z1charge1 = "<< Z1.charge1<<" ,Z1charge2 = "<< Z1.charge2<<" ,loose ele charge = "<<RECOELE_CHARGE[iL_loose_e[i]]<<endl;
	 
	 if (Z1.charge1 + RECOELE_CHARGE[iL_loose_e[i]] == 0 ) //opposite charge
	   {
	     
	     Lepton1qcd.SetPtEtaPhiE(Z1.pt1, Z1.eta1, Z1.phi1, Z1.lep1E ); 
	     Lepton3qcd.SetPtEtaPhiE(RECOELE_PT[ iL_loose_e[i] ], RECOELE_ETA[iL_loose_e[i]], RECOELE_PHI[ iL_loose_e[i] ], RECOELE_E[ iL_loose_e[i] ] );
	     DiLeptonQCD=Lepton1qcd+Lepton3qcd;
	     
	     ZQCD = DiLeptonQCD.M();
	     
     }
	 else {
	   
	   Lepton2qcd.SetPtEtaPhiE(Z1.pt2, Z1.eta2, Z1.phi2, Z1.lep2E ); 
	   Lepton3qcd.SetPtEtaPhiE(RECOELE_PT[ iL_loose_e[i] ], RECOELE_ETA[iL_loose_e[i]], RECOELE_PHI[ iL_loose_e[i] ], RECOELE_E[ iL_loose_e[i] ] );
	   DiLeptonQCD=Lepton2qcd+Lepton3qcd;
	   
	   ZQCD = DiLeptonQCD.M();
	 }
	 
	 
	 if( Z1.massvalue < min_mass_2L ) min_mass_2L = Z1.massvalue ;
	 if( ZQCD < min_mass_2L ) min_mass_2L = ZQCD ;
	 
	 //cout<<"mmin_mass_2L = "<<min_mass_2L<<endl;
	 
	 if (min_mass_2L <=4)continue;
	 
	 ///////////////////////////////////////////////////////////////////////////////////////////////
	 
	 bool BDT_ok = 0; // Fall15 with CMSSW_7_6_x
	 
	 if( RECOELE_PT[i] > 7. &&  RECOELE_PT[i] <= 10. ){
	   if( fabs(RECOELE_scl_Eta[i]) < .8 && RECOELE_mvaNonTrigV0[i] > -0.211 ) BDT_ok = 1 ;
	   if( ( fabs(RECOELE_scl_Eta[i]) >= .8 && fabs(RECOELE_scl_Eta[i]) < 1.479 )
	       && RECOELE_mvaNonTrigV0[i] > -0.396 ) BDT_ok = 1 ;
	   if( fabs(RECOELE_scl_Eta[i]) >= 1.479 && RECOELE_mvaNonTrigV0[i] > -0.215 ) BDT_ok = 1 ;
	 }
	 else { 
	   if( fabs(RECOELE_scl_Eta[i]) < .8 && RECOELE_mvaNonTrigV0[i] > -0.870 ) BDT_ok = 1 ;
	   if( ( fabs(RECOELE_scl_Eta[i]) >= .8 && fabs(RECOELE_scl_Eta[i]) <= 1.479 )
	       && RECOELE_mvaNonTrigV0[i] > -0.838 ) BDT_ok = 1 ;
	   if( fabs(RECOELE_scl_Eta[i]) > 1.479 && RECOELE_mvaNonTrigV0[i] > -0.763 ) BDT_ok = 1 ;
	 }
	 
	 //////////////////////////////////////////////////////////////////////////////////////

	 if( massZ1 <= (Zmass-7.) || massZ1 >= (Zmass+7.) ) continue;
	 
	 //Draw MZ1 and PFMET in CR Z1+L
	 //=============================
	 
	 //Before cut (DEN)
	 
	 if( fabs(RECOELE_ETA[iL_loose_e[i]]) <= 1.479 ) {
	   hMZ1_loose_ele_barrel->Fill(massZ1,newweight_3); 
	   hPFMET_loose_ele_barrel->Fill(RECO_PFMET,newweight);
	   hPT_loose_ele_barrel->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	   MET_PT_Ele_Barrel_Den->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);
	   MET_PT_Ele_Barrel_Den_Rebin->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);
	   //charge
	   if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==+1 ) {
	     h_PT_loose_ele_barrel_pos->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	     h_PFMET_loose_ele_barrel_pos->Fill(RECO_PFMET,newweight);
	     h_MET_PT_loose_ele_barrel_pos->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight); }//positive
	   
	   else if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==-1 ){	     
	     h_PT_loose_ele_barrel_neg->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	     h_PFMET_loose_ele_barrel_neg->Fill(RECO_PFMET,newweight);
	     h_MET_PT_loose_ele_barrel_neg->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);  }//negative 
	 }//Barrel 
	 
	 else if( fabs(RECOELE_ETA[iL_loose_e[i]]) > 1.479) {
	   hMZ1_loose_ele_endcap->Fill(massZ1,newweight_3); 
	   hPFMET_loose_ele_endcap->Fill(RECO_PFMET,newweight);
	   hPT_loose_ele_endcap->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	   MET_PT_Ele_Endcap_Den->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);
	   MET_PT_Ele_Endcap_Den_Rebin->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);
	   //charge
	   if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==+1 ) {
	     h_PT_loose_ele_endcap_pos->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	     h_PFMET_loose_ele_endcap_pos->Fill(RECO_PFMET,newweight);
	     h_MET_PT_loose_ele_endcap_pos->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight); }//positive
	   
	   else if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==-1 ){	     
	     h_PT_loose_ele_endcap_neg->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	     h_PFMET_loose_ele_endcap_neg->Fill(RECO_PFMET,newweight);
	     h_MET_PT_loose_ele_endcap_neg->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);  }//negative
	 }//endcap 
	 
	 
	 //After cut (NUM)
	 
	 if( BDT_ok && RECOELE_PFX_rho_new[iL_loose_e[i]]<0.35 ) {
	   
	   if( fabs(RECOELE_ETA[iL_loose_e[i]]) <= 1.479 ) {
	     hMZ1_tight_ele_barrel->Fill(massZ1,newweight_3); 
	     hPFMET_tight_ele_barrel->Fill(RECO_PFMET,newweight);
	     hPT_tight_ele_barrel->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	     MET_PT_Ele_Barrel_Num->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);
	     MET_PT_Ele_Barrel_Num_Rebin->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);
	     //charge
	     if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==+1 ) {
	       h_PT_tight_ele_barrel_pos->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	       h_PFMET_tight_ele_barrel_pos->Fill(RECO_PFMET,newweight);
	       h_MET_PT_tight_ele_barrel_pos->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight); }//positive
	     
	     else if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==-1 ){	     
	       h_PT_tight_ele_barrel_neg->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	       h_PFMET_tight_ele_barrel_neg->Fill(RECO_PFMET,newweight);
	       h_MET_PT_tight_ele_barrel_neg->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);  }//negative 
	   }//Barrel
	   
	   else if( fabs(RECOELE_ETA[iL_loose_e[i]]) > 1.479 ) {
	     hMZ1_tight_ele_endcap->Fill(massZ1,newweight_3); 
	     hPFMET_tight_ele_endcap->Fill(RECO_PFMET,newweight);
	     hPT_tight_ele_endcap->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	     MET_PT_Ele_Endcap_Num->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);
	     MET_PT_Ele_Endcap_Num_Rebin->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);
	     //charge
	     if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==+1 ) {
	       h_PT_tight_ele_endcap_pos->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	       h_PFMET_tight_ele_endcap_pos->Fill(RECO_PFMET,newweight);
	       h_MET_PT_tight_ele_endcap_pos->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight); }//positive

	     else if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==-1 ){	     
	       h_PT_tight_ele_endcap_neg->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	       h_PFMET_tight_ele_endcap_neg->Fill(RECO_PFMET,newweight);
	       h_MET_PT_tight_ele_endcap_neg->Fill(RECO_PFMET,RECOELE_PT[iL_loose_e[i]],newweight);  }//negative
	   }//endcap
	   
	 }//Tight id ,iso 
	 
	 ///////////////////////////////////////////////////////////////////////////////////////
	 
	 //use PFMET <25 or use PFMET <30
	 
	 
	 if(RECO_PFMET>=25) continue ;
	 if( fabs(massZ1 - Zmass) >= 7. ) continue; //this suppress QCD effect
	 
	  //Info from the event
	  f_run   = Run;
	  f_lumi  = LumiSection;
	  f_event = Event;
	  f_weight = newweight;	   
	  f_zee_e_pt_loose = RECOELE_PT[iL_loose_e[i]];
	  f_zee_e_eta_loose = RECOELE_ETA[iL_loose_e[i]];
	  f_zee_e_charge_loose = RECOELE_CHARGE[iL_loose_e[i]];

	 // if( massZ1 <= (Zmass-7.) || massZ1 >= (Zmass+7.) ) continue; //the original place
	 
	    if( fabs(RECOELE_ETA[iL_loose_e[i]]) <= 1.479 ){
	      ZplusE_Pt_DEN_Barrel->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	      if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==+1 ) {ZplusE_Pt_DEN_Barrel_pos->Fill(RECOELE_PT[iL_loose_e[i]],newweight);}
	      else if (  (RECOELE_CHARGE[iL_loose_e[i]]) ==-1 ) {ZplusE_Pt_DEN_Barrel_neg->Fill(RECOELE_PT[iL_loose_e[i]],newweight);}
	    }//Barrel
	    
	    else if( fabs(RECOELE_ETA[iL_loose_e[i]]) > 1.479 ){
	      ZplusE_Pt_DEN_Endcaps->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
	      if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==+1 ) {ZplusE_Pt_DEN_Endcaps_pos->Fill(RECOELE_PT[iL_loose_e[i]],newweight);}
	      else if (  (RECOELE_CHARGE[iL_loose_e[i]]) ==-1 ) {ZplusE_Pt_DEN_Endcaps_neg->Fill(RECOELE_PT[iL_loose_e[i]],newweight);}
	    }//endcap
	    
	    
	    if( BDT_ok && RECOELE_PFX_rho_new[iL_loose_e[i]]<0.35  ){
	      f_zee_e_pt_tight = RECOELE_PT[iL_loose_e[i]];
	      f_zee_e_eta_tight = RECOELE_ETA[iL_loose_e[i]];
	      f_zee_e_charge_tight = RECOELE_CHARGE[iL_loose_e[i]];
	      
	      if( fabs(RECOELE_ETA[iL_loose_e[i]]) <= 1.479 ){
		ZplusE_Pt_NUM_ID_ISO_Barrel->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
		if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==+1 ) {ZplusE_Pt_NUM_ID_ISO_Barrel_pos->Fill(RECOELE_PT[iL_loose_e[i]],newweight);}
		else if (  (RECOELE_CHARGE[iL_loose_e[i]]) ==-1 ) {ZplusE_Pt_NUM_ID_ISO_Barrel_neg->Fill(RECOELE_PT[iL_loose_e[i]],newweight);}
	      }
	      else if( fabs(RECOELE_ETA[iL_loose_e[i]]) > 1.479 ){
		ZplusE_Pt_NUM_ID_ISO_Endcaps->Fill(RECOELE_PT[iL_loose_e[i]],newweight);
		if ( (RECOELE_CHARGE[iL_loose_e[i]]) ==+1 ) {ZplusE_Pt_NUM_ID_ISO_Endcaps_pos->Fill(RECOELE_PT[iL_loose_e[i]],newweight);}
		else if (  (RECOELE_CHARGE[iL_loose_e[i]]) ==-1 ) {ZplusE_Pt_NUM_ID_ISO_Endcaps_neg->Fill(RECOELE_PT[iL_loose_e[i]],newweight);}
	      }
	      
	    }//Tight id ,iso	 
	    
	    //sprintf (Eventformat,"FR=%d:%d:%d",Run,LumiSection,Event);
	    
       }//end loop on loose electrons for(int i = 0; i < N_loose_e ; i++)
       
       //end of electron Fake rate
       
     }//end else if ( Z1tag ==2 )
     
     
     //cout<<"finished"<<endl;
     if(f_run != -999 && f_lumi != -999 && f_event != -999){
       /*
       cout<< "***********************************************************"
       << "\nRun:    " << f_run
       << "\nLumi:   " << f_lumi
       << "\nEvent:  " << f_event
       << "\nWeight: " << f_weight
       << "\nf_zee_mu_pt_loose: " << f_zee_mu_pt_loose
       << "\nf_zee_mu_eta_loose: " << f_zee_mu_eta_loose
       << "\nf_zee_e_pt_loose: " << f_zee_e_pt_loose
       << "\nf_zee_e_eta_loose: " << f_zee_e_eta_loose
       << "\nf_zmumu_mu_pt_loose: " << f_zmumu_mu_pt_loose
       << "\nf_zmumu_mu_eta_loose: " << f_zmumu_mu_eta_loose
       << "\nf_zmumu_e_pt_loose: " << f_zmumu_e_pt_loose
       << "\nf_zmumu_e_eta_loose: " << f_zmumu_e_eta_loose
       << "\nf_zee_mu_pt_tight: " << f_zee_mu_pt_tight
       << "\nf_zee_mu_eta_tight: " << f_zee_mu_eta_tight
       << "\nf_zee_e_pt_tight: " << f_zee_e_pt_tight
       << "\nf_zee_e_eta_tight: " << f_zee_e_eta_tight
       << "\nf_zmumu_mu_pt_tight: " << f_zmumu_mu_pt_tight
       << "\nf_zmumu_mu_eta_tight: " << f_zmumu_mu_eta_tight
       << "\nf_zmumu_e_pt_tight: " << f_zmumu_e_pt_tight
       << "\nf_zmumu_e_eta_tight: " << f_zmumu_e_eta_tight
       << "\nf_zee_mu_charge_loose: " << f_zee_mu_charge_loose
       << "\nf_zee_e_charge_loose: " << f_zee_e_charge_loose
       << "\nf_zmumu_mu_charge_loose: " << f_zmumu_mu_charge_loose
       << "\nf_zmumu_e_charge_loose: " << f_zmumu_e_charge_loose
       << "\nf_zee_mu_charge_tight: " << f_zee_mu_charge_tight
       << "\nf_zee_e_charge_tight: " << f_zee_e_charge_tight
       << "\nf_zmumu_mu_charge_tight: " << f_zmumu_mu_charge_tight
       << "\nf_zmumu_e_charge_tight: " << f_zmumu_e_charge_tight
       <<endl;
       */
       treeFR->Fill();
     }
     

//===============================================================================================================================
// Builds the Control Regions (CR) 2P2F and 3P1F (Opposite Signal (OS) Method)
//===============================================================================================================================
     //Reset weight
     //Weights
     //========
     if(MC_type == "Spring16"){
        float lumifb = 35.867;
	weight = lumifb*(Xsection*1000.*neventsPostHLT/neventsPreHLT)/neventsPostHLT;
     }
     newweight=weight;
     // pileup reweighting 2016
     //-----------------------------------
     if (DATA_type=="NO" && num_PU_vertices < 0) continue;  //make sure
     pu_weight=1.;
     if (MC_type == "Spring16"){
       Int_t binx = puweight->GetXaxis()->FindBin(num_PU_vertices);
       ////  cout << " bin x= " << binx << " " << puweight->GetBinContent(binx) << endl;
       pu_weight=double(puweight->GetBinContent(binx));
     }
     newweight=weight*pu_weight;
     
     
     //Build Z2 candidates from electrons
     //====================================
     
     int Z2xx_tag=0;
     int Z2_lepton1tag = 0; // 1 for mu ; 2 for ele
     int Z2_lepton2tag = 0;
     
     // if( N_loose_e < 2 ) continue ; 
     
     //loop on loose electrons
     
     vector<candidateZ> Z2_candvector;
     Z2_candvector.clear();
     
     
     // array<candidateZ,2> Zs;
     // vector<std::array<candidateZ, 2> > Zsv;
     
     vector<std::pair<candidateZ, candidateZ> > ZZpair;
     ZZpair.clear();
     
     
     for (int index=0;index<Zcandvector.size();index++){     
       
       for(int i = 0; i < N_loose_e; i++){
	 
	 if (fabs(RECOELE_SIP[iL_loose_e[i]])>=4.)continue; //SIP cut
	 
	 //check that this third muon not 2 muons from Z1
	 
	 //if Z1 is Zee
	 
	 ////	 cout<<"iL_loose_e[i] = "<<iL_loose_e[i]<<endl;
	 //// cout<<"Zcandvector.at(index).tag = "<<Zcandvector.at(index).tag<<"and lepton indeces = "<<Zcandvector.at(index).ilept1<<", "<<Zcandvector.at(index).ilept2<<endl;
	 
	 
	 if (Zcandvector.at(index).tag == 2){if( iL_loose_e[i]==Zcandvector.at(index).ilept1 || iL_loose_e[i]==Zcandvector.at(index).ilept2 ) continue ;}
	 
	 //// cout<<"escape this 1st electron"<<endl;
	 
	 for(int j = i + 1; j < N_loose_e; j++){
	   
	   
	   if (fabs(RECOELE_SIP[iL_loose_e[j]])>=4.)continue; //SIP cut
	   
	   //check that this third muon not 2 muons from Z1
	   
	   //if Z1 is Zee
	   
	   //// cout<<"iL_loose_e[j] = "<<iL_loose_e[j]<<endl;
	   
	   if (Zcandvector.at(index).tag == 2){if( iL_loose_e[j]==Zcandvector.at(index).ilept1 || iL_loose_e[j]==Zcandvector.at(index).ilept2 ) continue ;}
	   
	   ////  cout<<"escape this 2nd electron"<<endl;
	   
	   if(RECOELE_CHARGE[iL_loose_e[i]] == RECOELE_CHARGE[iL_loose_e[j]]) continue; //opposite charge
	   
	   
	   double pxZ, pyZ, pzZ;
	   double EZ;
	   double massZ;
	   double massZ_noFSR = 0.;
	   
	   double ptZ = 0.;
	   double Y_Z = -9.;
	   double sum_ptZ = 0.;
	   
	   int tempphotid1=-1;
	   int tempphotid2=-1;
	   int templepid=-1;
	   
	   float pTphot1=-999.;
	   float pTphot2=-999.;
	   
	   //Zmass with noFSR
	   //-------------------	     
	   
	   
	   Lepton1.SetPtEtaPhiM (RECOELE_PT[iL_loose_e[i]],RECOELE_ETA[iL_loose_e[i]], RECOELE_PHI[iL_loose_e[i]], 0.000511);
	   Lepton2.SetPtEtaPhiM (RECOELE_PT[iL_loose_e[j]],RECOELE_ETA[iL_loose_e[j]], RECOELE_PHI[iL_loose_e[j]], 0.000511);
	   
	   DiLepton=Lepton1+Lepton2;
	   massZ = DiLepton.M();
	   ptZ = DiLepton.Pt();
	   pxZ = DiLepton.Px();
	   pyZ = DiLepton.Py();
	   pzZ = DiLepton.Pz();
	   EZ = DiLepton.E();
	   Y_Z = DiLepton.Rapidity();
	   
	   sum_ptZ = RECOELE_PT[ iL_loose_e[i] ] + RECOELE_PT[ iL_loose_e[j] ];
	   massZ_noFSR = massZ;
	   
	   Z2xx_tag=2;//Z decay to electrons
	   Z2_lepton1tag = 2; // 1 for mu ; 2 for ele
	   Z2_lepton2tag = 2;
	   
	   //if( debug )  cout<< " noFSR Z1 mass = "<< massZ<<"and pt= "<<ptZ<<"and px "<<pxZ<<"and py= "<<pyZ<<"and pz= "<<pzZ<<"and rapidity = "<<Y_Z<<"and sum_pTZ= "<<sum_ptZ<<endl;
	   
	   
	   // ** Association of FSR to Z
	   //=============================
	   
	   bool has_FSR_Z = 0;
	   int N_FSR_Z = 0;
	   double max_pt_FSR_Z = -1.;
	   int pi = -1; 
	   int pj = -1;
	   double mllp=-1;
	   
	   //here we will make loop on photons to know how many FSR in the event is it N_FSR =1 or > 1 and to know photon max pt
	   // we will loop on photons that made small deltaR with muons
	   
	   for (int p=0 ;p<Nphotons ; p++ ) {
	     
	     if( iLp_l[ p ] == iL_loose_e[i] && iLp_tagEM[ p ] == 1 ){
	       
	       //Evaluate the mass 
	       
	       LeptonCorrection.SetPtEtaPhiM (RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]], RECOPFPHOT_PHI[iLp[p]], 0);
	       
	       Lepton1= Lepton1+ LeptonCorrection;
	       DiLepton= Lepton1+Lepton2;
	       
	       mllp = DiLepton.M();	
	       ptZ = DiLepton.Pt();
	       pxZ = DiLepton.Px();
	       pyZ = DiLepton.Py();
	       pzZ = DiLepton.Pz();
	       EZ = DiLepton.E();
	       Y_Z = DiLepton.Rapidity();
	       
	       //cout<< " mllp = "<< mllp<<"and pt= "<<pt<<"and px "<<px<<"and py= "<<py<<"and pz= "<<pz<<"and rapidity = "<<Y<<endl;
	       
	       has_FSR_Z = 1; 
	       pi = p; 
	       N_FSR_Z++;
	       massZ=mllp;
	       if( RECOPFPHOT_PT[iLp[p]] > max_pt_FSR_Z ) max_pt_FSR_Z = RECOPFPHOT_PT[iLp[p]];
	       
	       
	     }  // end of if( iLp_l[ p ] == iL[i] && iLp_tagEM[ p ] == 0 )
	     
	     if( iLp_l[ p ] == iL_loose_e[j] && iLp_tagEM[ p ] == 1 ){
	       
	       //evaluate the mass	
	       
	       LeptonCorrection.SetPtEtaPhiM (RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]], RECOPFPHOT_PHI[iLp[p]], 0);
	       
	       Lepton2=Lepton2+LeptonCorrection;
	       DiLepton=Lepton1+Lepton2;
	       
	       mllp = DiLepton.M();	
	       ptZ = DiLepton.Pt();
	       pxZ = DiLepton.Px();
	       pyZ = DiLepton.Py();
	       pzZ = DiLepton.Pz();
	       EZ = DiLepton.E();
	       Y_Z = DiLepton.Rapidity();
	       
	       // cout<< " mllp = "<< mllp<<"and pt= "<<pt<<"and px "<<px<<"and py= "<<py<<"and pz= "<<pz<<"and rapidity = "<<Y<<endl;
	       
	       
		has_FSR_Z = 1; 
		pj = p; 
		N_FSR_Z++;
		massZ=mllp;
		if( RECOPFPHOT_PT[iLp[p]] > max_pt_FSR_Z ) max_pt_FSR_Z = RECOPFPHOT_PT[iLp[p]];
		
	     }//end of if( iLp_l[ p ] == iL[j] && iLp_tagEM[ p ] == 0 )
	     
	   }//end of for (int p=0 ;p<Nphotons ; p++ )
	   
	   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	   
	   
	   if(has_FSR_Z){ //there are 2 options N_FSR_Z==1 or >1
	     
	     
	     // isolation 
	     /*
	     if( debug ) cout  << "Z Isolation (corrected for photon): "
			       << "\n RECOELE_PFX_rho_new[iL_loose_e[i] ] " << RECOELE_PFX_rho_new[iL_loose_e[i]]
			       << "\n RECOELE_PFX_rho_new[iL_loose_e[j] ] " << RECOELE_PFX_rho_new[iL_loose_e[j]]
			       << endl;
	     */
	     if( pi != -1 ){
	       pTphot1=RECOPFPHOT_PT[iLp[pi]];
	       tempphotid1=iLp[pi];
	       // cout<<"pi = "<<pi<<"and pTphot1 = "<<pTphot1<<endl;	      
	     }
	     if( pj != -1 ){	      
	       pTphot2=RECOPFPHOT_PT[iLp[pj]];
	       tempphotid2=iLp[pj];
	       //  cout<<"pj = "<<pi<<"and pTphot2 = "<<pTphot2<<endl;      
	     }
	   }//end of if(has_FSR_Z)
	   
	   else { //if don't have FSR
	     
	     pi=-1;
	     pj=-1;
	     /*
	     if( debug ) cout  << "Z Isolation: "  
			       << "\n RECOELE_PFX_rho_new[ iL_loose_e[i] ] " << RECOELE_PFX_rho_new[ iL_loose_e[i] ]
			       << "\n RECOELE_PFX_rho_new[ iL_loose_e[j] ] " << RECOELE_PFX_rho_new[ iL_loose_e[j] ]
			       << endl;
	     */
	   }//end of else if don't have FSR
	   
	   
	   //////////////////////////////////////////////////////////////////////////////////////////
	   
	   
	   candidateZ *Z2e = new candidateZ; //this Z is a pointer to data structure candidateZ 
	   Z2e->massvalue=massZ;
	   Z2e->massvalue_NOFSR=massZ_noFSR;
	   Z2e->ilept1=iL_loose_e[i];
	   Z2e->ilept2=iL_loose_e[j];
	   Z2e->pt1=RECOELE_PT[iL_loose_e[i]];
	   Z2e->pt2=RECOELE_PT[iL_loose_e[j]];
	   Z2e->eta1=RECOELE_ETA[iL_loose_e[i]];
	   Z2e->eta2=RECOELE_ETA[iL_loose_e[j]];
	   Z2e->phi1=RECOELE_PHI[iL_loose_e[i]];
	   Z2e->phi2=RECOELE_PHI[iL_loose_e[j]];
	   Z2e->charge1=RECOELE_CHARGE[iL_loose_e[i]];
	   Z2e->charge2=RECOELE_CHARGE[iL_loose_e[j]];
	   Z2e->lep1E=RECOELE_E[iL_loose_e[i]];
	   Z2e->lep2E=RECOELE_E[iL_loose_e[j]];
	   Z2e->lep1tag=Z2_lepton1tag;
	   Z2e->lep2tag=Z2_lepton2tag;
	   Z2e->isol1=RECOELE_PFX_rho_new[iL_loose_e[i]];
	   Z2e->isol2=RECOELE_PFX_rho_new[iL_loose_e[j]];
	   if( pi != -1 ) {Z2e->ilept1_FSR=true;
	     /*Z2e->ilept2_FSR=false;*/}
	   if( pj != -1 ) {/*Z2e->ilept1_FSR=false;*/
	     Z2e->ilept2_FSR=true;}
	   Z2e->pxZ=pxZ;
	   Z2e->pyZ=pyZ;
	   Z2e->pzZ=pzZ;
	   Z2e->EZ=EZ;
	   if( has_FSR_Z ) {
	     Z2e->withFSR=1;
	     Z2e->ptFSR1=pTphot1;
	     Z2e->ptFSR2=pTphot2;
	     Z2e->iFSR1=tempphotid1;
	     Z2e->iFSR2=tempphotid2;
	   }	      
	   else {
	     Z2e->withFSR=0;
	     Z2e->ptFSR1=0.;
	     Z2e->ptFSR2=0.;
	     Z2e->iFSR1=-1;
	     Z2e->iFSR2=-1;
	     Z2e->ilept1_FSR=false;
	     Z2e->ilept2_FSR=false;
	   }
	   
	   Z2e->tag=Z2xx_tag;
	   
	   Z2_candvector.push_back(*Z2e);
	   
	   //Zs = {Zcandvector.at(index) , Z2e};
	   //Zsv.push_back(Zs);
	   
	   
	   ZZpair.push_back( {Zcandvector.at(index) , *Z2e} );
	   //ZZpair = std::make_pair (Zcandvector.at(index),Z2e);
	   
	 }//end of for(int j = i + 1; j < N_loose_e; j++)
       }//end of for(int i = 0; i < N_loose_e; i++)
       
     }//end for (int i=0;i<Zcandvector.size();i++)
     
     
     // cout<<"loose electron indices "<<iL_loose_e[0]<<","<<iL_loose_e[1]<<","<<iL_loose_e[2]<<","<<iL_loose_e[3]<<","<<iL_loose_e[4]<<","<<iL_loose_e[5]<<endl;
     
     // cout<<"Z2_candvector for 2 electrons = "<<Z2_candvector.size()<<endl;
     
     // for (int i=0; i<Z2_candvector.size();i++){
     //   cout <<"Z2 has mass= " <<Z2_candvector.at(i).massvalue<<" and tag = "<<Z2_candvector.at(i).tag<<" and leptons indices = "<<Z2_candvector.at(i).ilept1<<","<<Z2_candvector.at(i).ilept2<<"with pt "<<Z2_candvector.at(i).pt1<<","<<Z2_candvector.at(i).pt2<<" and leptons tag = "<<Z2_candvector.at(i).lep1tag<<" , "<<Z2_candvector.at(i).lep2tag<<" iso = "<<Z2_candvector.at(i).isol1<<" , "<<Z2_candvector.at(i).isol2<<endl;
     //     }
     
     //cout<<"ZZ pairs from electrons = "<<ZZpair.size()<<endl;
     
     //====================================
     // Build Z2 candidates from muons  //
     //====================================
     
     
     for (int index=0;index<Zcandvector.size();index++){
       
       for(int i = 0; i < N_loose_mu; i++){
	 
	 if (fabs( RECOMU_SIP[iL_loose_mu[i]] ) >= 4.)continue;
	 
	 //check that this third muon not 2 muons from Z1
	 
	 //Z1 if Zmumu
	 
	 ////	 cout<<"iL_loose_mu[i] = "<<iL_loose_mu[i]<<endl;
	 ////	 cout<<"Zcandvector.at(index).tag = "<<Zcandvector.at(index).tag<<"and lepton indeces = "<<Zcandvector.at(index).ilept1<<", "<<Zcandvector.at(index).ilept2<<endl;
	 
	 if (Zcandvector.at(index).tag == 1) { if( iL_loose_mu[i]==Zcandvector.at(index).ilept1 || iL_loose_mu[i]==Zcandvector.at(index).ilept2 ) continue ;}
	 
	 ////	 cout<<"escape this 1st muon"<<endl;
	 
	 for(int j = i + 1; j < N_loose_mu; j++){
	   
	   
	   if (fabs( RECOMU_SIP[iL_loose_mu[j]] ) >= 4.)continue;
	   
	   //check that this third muon not 2 muons from Z1
	   
	   //Z1 if Zmumu
	   
	   //// cout<<"iL_loose_mu[j] = "<<iL_loose_mu[j]<<endl;
	   
	   if (Zcandvector.at(index).tag == 1) { if( iL_loose_mu[j]==Zcandvector.at(index).ilept1 || iL_loose_mu[j]==Zcandvector.at(index).ilept2 ) continue ;}
	   
	   ////  cout<<"escape this 2nd muon"<<endl;
	   
	   if(RECOMU_CHARGE[iL_loose_mu[i]] == RECOMU_CHARGE[iL_loose_mu[j]]) continue; //opposite charge
	   
	   
	   double pxZ, pyZ, pzZ;
	   double EZ;
	   double massZ;
	   double massZ_noFSR = 0.;
	   
	   double ptZ = 0.;
	   double Y_Z = -9.;
	   double sum_ptZ = 0.;
	   
	   int tempphotid1=-1;
	   int tempphotid2=-1;
	   int templepid=-1;
	   
	   float pTphot1=-999.;
	   float pTphot2=-999.;
	   
	   //evaluate the mass 
	   
	   Lepton1.SetPtEtaPhiM (RECOMU_PT[iL_loose_mu[i]],RECOMU_ETA[iL_loose_mu[i]], RECOMU_PHI[iL_loose_mu[i]], 0.105);
	   Lepton2.SetPtEtaPhiM (RECOMU_PT[iL_loose_mu[j]],RECOMU_ETA[iL_loose_mu[j]], RECOMU_PHI[iL_loose_mu[j]], 0.105);
	   
	   DiLepton=Lepton1+Lepton2;
	   massZ = DiLepton.M();
	   ptZ = DiLepton.Pt();
	   pxZ = DiLepton.Px();
	   pyZ = DiLepton.Py();
	   pzZ = DiLepton.Pz();
	   EZ = DiLepton.E();
	   Y_Z = DiLepton.Rapidity();
	   
	   sum_ptZ = RECOMU_PT[ iL_loose_mu[i] ] + RECOMU_PT[ iL_loose_mu[i] ];
	   massZ_noFSR = massZ;
	   
	   Z2xx_tag=1;//Z decay to muons
	   Z2_lepton1tag = 1; // 1 for mu ; 2 for ele
	   Z2_lepton2tag = 1;
	   
	   //if( debug )  cout<< " noFSR Z1 mass = "<< massZ<<"and pt= "<<ptZ<<"and px "<<pxZ<<"and py= "<<pyZ<<"and pz= "<<pzZ<<"and rapidity = "<<Y_Z<<"and sum_pTZ= "<<sum_ptZ<<endl;
	
	   
	   // ** Association of FSR to Z
	   //=============================
	   
	   bool has_FSR_Z = 0;
	   int N_FSR_Z = 0;
	   double max_pt_FSR_Z = -1.;
	   int pi = -1; 
	   int pj = -1;
	   double mllp=-1;
	   

	   //here we will make loop on photons to know how many FSR in the event is it N_FSR =1 or > 1 and to know photon max pt
	   // we will loop on photons that made small deltaR with muons

	   for (int p=0 ;p<Nphotons ; p++ ) {
	     
	     if( iLp_l[ p ] == iL_loose_mu[i] && iLp_tagEM[ p ] == 0 ){
	       
	       //evaluate the mass 
	       
	       LeptonCorrection.SetPtEtaPhiM (RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]], RECOPFPHOT_PHI[iLp[p]], 0);
	       
	       Lepton1= Lepton1+ LeptonCorrection;
	       DiLepton= Lepton1+Lepton2;
	       
	       mllp = DiLepton.M();	
	       ptZ = DiLepton.Pt();
	       pxZ = DiLepton.Px();
	       pyZ = DiLepton.Py();
	       pzZ = DiLepton.Pz();
	       EZ = DiLepton.E();
	       Y_Z = DiLepton.Rapidity();
	       
	       //cout<< " mllp = "<< mllp<<"and pt= "<<pt<<"and px "<<px<<"and py= "<<py<<"and pz= "<<pz<<"and rapidity = "<<Y<<endl;
	       
	       
	       has_FSR_Z = 1; 
	       pi = p; 
	       N_FSR_Z++;
	       massZ=mllp;
	       if( RECOPFPHOT_PT[iLp[p]] > max_pt_FSR_Z ) max_pt_FSR_Z = RECOPFPHOT_PT[iLp[p]];
	       
	     }  // end of if( iLp_l[ p ] == iL[i] && iLp_tagEM[ p ] == 0 )
	     
	     if( iLp_l[ p ] == iL_loose_mu[j] && iLp_tagEM[ p ] == 0 ){
	       
	       //evaluate the mass 
	       
	       LeptonCorrection.SetPtEtaPhiM (RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]], RECOPFPHOT_PHI[iLp[p]], 0);

	       Lepton2=Lepton2+LeptonCorrection;
	       DiLepton=Lepton1+Lepton2;
	       
	       mllp = DiLepton.M();	
	       ptZ = DiLepton.Pt();
	       pxZ = DiLepton.Px();
	       pyZ = DiLepton.Py();
	       pzZ = DiLepton.Pz();
	       EZ = DiLepton.E();
	       Y_Z = DiLepton.Rapidity();
	       
	       // cout<< " mllp = "<< mllp<<"and pt= "<<pt<<"and px "<<px<<"and py= "<<py<<"and pz= "<<pz<<"and rapidity = "<<Y<<endl;
	       
	       has_FSR_Z = 1; 
		pj = p; 
		N_FSR_Z++;
		massZ=mllp;
		if( RECOPFPHOT_PT[iLp[p]] > max_pt_FSR_Z ) max_pt_FSR_Z = RECOPFPHOT_PT[iLp[p]];
		
	     }//end of if( iLp_l[ p ] == iL[j] && iLp_tagEM[ p ] == 0 )
	     
	   }//end of for (int p=0 ;p<Nphotons ; p++ )
	  
	   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	   // Start adding FSR 
	   //===================
	   
	   if(has_FSR_Z){ //there are 2 options N_FSR_Z==1 or >1
	     
	     // isolation 
	     /*
	     if( debug ) cout  << "Z Isolation ( corrected for photon): "
			       << "\n RECOMU_PFX_dB_new[iL_loose_mu[i] ] " << RECOMU_PFX_dB_new[iL_loose_mu[i]]
			       << "\n RECOMU_PFX_dB_new[iL_loose_mu[j] ] " << RECOMU_PFX_dB_new[iL_loose_mu[j]]
			       << endl;
	     */
	     if( pi != -1 ){
	       pTphot1=RECOPFPHOT_PT[iLp[pi]];
	       tempphotid1=iLp[pi];
	       // cout<<"pi = "<<pi<<"and pTphot1 = "<<pTphot1<<endl;	      
	     }
	     if( pj != -1 ){	      
	       pTphot2=RECOPFPHOT_PT[iLp[pj]];
	       tempphotid2=iLp[pj];
	       //  cout<<"pj = "<<pi<<"and pTphot2 = "<<pTphot2<<endl;	      
	     }
	   }//end of if(has_FSR_Z)
	   
	   else { //if don't have FSR
	     
	     pi=-1;
	     pj=-1;
	     /*
	     if( debug ) cout  << "Z Isolation: "  
			       << "\n RECOMU_PFX_dB_new[ iL_loose_mu[i] ] " << RECOMU_PFX_dB_new[ iL_loose_mu[i] ]
			       << "\n RECOMU_PFX_dB_new[ iL_loose_mu[j] ] " << RECOMU_PFX_dB_new[ iL_loose_mu[j] ]
			       << endl;
	     */
	   }//end of else if don't have FSR
	   
	   
	   //////////////////////////////////////////////////////////////////////////////////////////
	   
	   candidateZ *Z2mu = new candidateZ; //this Z is a pointer to data structure candidateZ 
	   Z2mu->massvalue=massZ;
	   Z2mu->massvalue_NOFSR=massZ_noFSR;
	   Z2mu->ilept1=iL_loose_mu[i];
	   Z2mu->ilept2=iL_loose_mu[j];
	   Z2mu->pt1=RECOMU_PT[iL_loose_mu[i]];
	   Z2mu->pt2=RECOMU_PT[iL_loose_mu[j]];
	   Z2mu->eta1=RECOMU_ETA[iL_loose_mu[i]];
	   Z2mu->eta2=RECOMU_ETA[iL_loose_mu[j]];
	   Z2mu->phi1=RECOMU_PHI[iL_loose_mu[i]];
	   Z2mu->phi2=RECOMU_PHI[iL_loose_mu[j]];
	   Z2mu->charge1=RECOMU_CHARGE[iL_loose_mu[i]];
	   Z2mu->charge2=RECOMU_CHARGE[iL_loose_mu[j]];
	   Z2mu->lep1E=RECOMU_E[iL_loose_mu[i]];
	   Z2mu->lep2E=RECOMU_E[iL_loose_mu[j]];
	   Z2mu->lep1tag=Z2_lepton1tag;
	   Z2mu->lep2tag=Z2_lepton2tag;
	   Z2mu->isol1=RECOMU_PFX_dB_new[iL_loose_mu[i]];
	   Z2mu->isol2=RECOMU_PFX_dB_new[iL_loose_mu[j]];
	   if( pi != -1 ) {Z2mu->ilept1_FSR=true;
	     /*Z2mu->ilept2_FSR=false;*/}
	   if( pj != -1 ) {/*Z2mu->ilept1_FSR=false;*/
	     Z2mu->ilept2_FSR=true;}
	   Z2mu->pxZ=pxZ;
	   Z2mu->pyZ=pyZ;
	   Z2mu->pzZ=pzZ;
	   Z2mu->EZ=EZ;
	   if( has_FSR_Z ) {
	     Z2mu->withFSR=1;
	     Z2mu->ptFSR1=pTphot1;
	     Z2mu->ptFSR2=pTphot2;
	     Z2mu->iFSR1=tempphotid1; //photon index
	     Z2mu->iFSR2=tempphotid2;
	   }	      
	   else {
	     Z2mu->withFSR=0;
	     Z2mu->ptFSR1=0.;
	     Z2mu->ptFSR2=0.;
	     Z2mu->iFSR1=-1;
	     Z2mu->iFSR2=-1;
	     Z2mu->ilept1_FSR=false;
	     Z2mu->ilept2_FSR=false;
	   }
	   
	   Z2mu->tag=Z2xx_tag;
	   Z2_candvector.push_back(*Z2mu);
	   
	   ZZpair.push_back( {Zcandvector.at(index) , *Z2mu} );
	   
	   
	 }//end of for(int j = i + 1; j < N_loose; j++)
       }//end of for(int i = 0; i < N_loose; i++)
       
     }//end for (int i=0;i<Zcandvector.size();i++)
     
     
     
     //// cout<<"Z2 cand = "<<Z2_candvector.size()<<endl;
     /*   
     for (int i=0; i<Z2_candvector.size(); i++){
       
       //cout<<"Z2_candvector = "<< Z2_candvector.size()<<" has Z"<<i<<"with tag = "<<Z2_candvector.at(i).tag<<" and mass = "<<Z2_candvector.at(i).massvalue<<" and leptons indices = "<<Z2_candvector.at(i).ilept1<<" , "<<Z2_candvector.at(i).ilept2<<" and leptons tag = "<<Z2_candvector.at(i).lep1tag<<" , "<<Z2_candvector.at(i).lep2tag<<" and pTs= "<<Z2_candvector.at(i).pt1<<" , "<<Z2_candvector.at(i).pt2<<" iso = "<<Z2_candvector.at(i).isol1<<" , "<<Z2_candvector.at(i).isol2<<endl;
       
        cout<<"Z2_canvector properties tag = "<<Z2_candvector.at(i).tag<<" ,mass =  "<<Z2_candvector.at(i).massvalue<<",WithFSR = "<<Z2_candvector.at(i).withFSR <<"photon indeces = "<<Z2_candvector.at(i).iFSR1<<","<<Z2_candvector.at(i).iFSR2<<" ,photons pT = "<<Z2_candvector.at(i).ptFSR1<<" , "<<Z2_candvector.at(i).ptFSR2<<" ,IF leptons have FSR = "<<Z2_candvector.at(i).ilept1_FSR<<" , "<<Z2_candvector.at(i).ilept2_FSR<<" and leptons indices = "<<Z2_candvector.at(i).ilept1<<","<<Z2_candvector.at(i).ilept2<<"and pt = "<<Z2_candvector.at(i).pt1<<","<<Z2_candvector.at(i).pt2<<" and leptons isolation = "<<Z2_candvector.at(i).isol1<<", "<<Z2_candvector.at(i).isol2<<" and leptons tag = "<<Z2_candvector.at(i).lep1tag<<" , "<<Z2_candvector.at(i).lep2tag<<endl;
       
     }
     
     */   
     //// cout<<"number of ZZpairs = "<<ZZpair.size()<<endl;
     
     /*  for(int i=0; i<ZZpair.size(); i++ ){
       
       cout<<"ZZ pairs masses = "<<ZZpair[i].first.massvalue<<" , "<<ZZpair[i].second.massvalue<<"and tags = "<<ZZpair[i].first.tag<<", "<<ZZpair[i].second.tag<<endl;

       
       }*/
     
     if ( ZZpair.size() ==0 )continue;
     
     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     
     
     // Mass cut on all ZZ candidates 
     //-------------------------------
     //check the mass of all ZZ candidates to be massZ > 12  && massZ <120   
     
     
     vector<std::pair<candidateZ, candidateZ> > ZZpair_masscut;
     ZZpair_masscut.clear();
     
     for (int index=0; index<ZZpair.size();index++){
       
       if (ZZpair[index].first.massvalue <= 40 || ZZpair[index].first.massvalue >= 120) continue;
       if (ZZpair[index].second.massvalue <= 12 || ZZpair[index].second.massvalue >= 120) continue;
       
       ZZpair_masscut.push_back({ZZpair[index].first , ZZpair[index].second});
     }
     
     
     //// cout<<"number of ZZpairs after mass cut = "<<ZZpair_masscut.size()<<endl;
     
     /*  for(int i=0; i<ZZpair_masscut.size(); i++ ){
       
       cout<<"ZZ pairs after mass cuts have masses = "<<ZZpair_masscut[i].first.massvalue<<" , "<<ZZpair_masscut[i].second.massvalue<<"and tags = "<<ZZpair_masscut[i].first.tag<<", "<<ZZpair_masscut[i].second.tag<<endl;
       
       
       }*/
     
     if ( ZZpair_masscut.size() ==0 )continue;
     
     //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     
     //deltaR > 0.02 for all ZZ pairs
     //-------------------------------
     
     vector<std::pair<candidateZ, candidateZ> > ZZpair_deltaRcut;
     ZZpair_deltaRcut.clear();
     
     
     for(int index=0; index<ZZpair_masscut.size(); index++){
       
       // cout<<"Z1_lep1 Z1_lep2 deltaR = "<<sqrt( pow( DELTAPHI(  ZZpair_masscut[index].first.phi1 , ZZpair_masscut[index].first.phi2 ),2) + pow( ZZpair_masscut[index].first.eta1 - ZZpair_masscut[index].first.eta2 ,2) )<<endl;
       // cout<<"Z1_lep1 Z2_lep1 deltaR = "<<sqrt( pow( DELTAPHI( ZZpair_masscut[index].first.phi1 , ZZpair_masscut[index].second.phi1  ),2) + pow( ZZpair_masscut[index].first.eta1 - ZZpair_masscut[index].second.eta1 ,2) )<<endl;
       // cout<<"Z1_lep1 Z2_lep2 deltaR = "<<sqrt( pow( DELTAPHI(  ZZpair_masscut[index].first.phi1 , ZZpair_masscut[index].second.phi2  ),2) + pow( ZZpair_masscut[index].first.eta1 - ZZpair_masscut[index].second.eta2 ,2) )<<endl;
       // cout<<"Z1_lep2 Z2lep1 deltaR = "<<sqrt( pow( DELTAPHI(  ZZpair_masscut[index].first.phi2 , ZZpair_masscut[index].second.phi1 ),2) + pow( ZZpair_masscut[index].first.eta2 - ZZpair_masscut[index].second.eta1 ,2) )<<endl;
       // cout<<"Z1_lep2 Z2_lep2 deltaR = "<<sqrt( pow( DELTAPHI(  ZZpair_masscut[index].first.phi2 , ZZpair_masscut[index].second.phi2 ),2) + pow( ZZpair_masscut[index].first.eta2 - ZZpair_masscut[index].second.eta2 ,2) )<<endl;
       // cout<<"Z2_lep1 Z2_lep2 deltaR = "<<sqrt( pow( DELTAPHI(  ZZpair_masscut[index].second.phi1 ,  ZZpair_masscut[index].second.phi2 ),2) + pow( ZZpair_masscut[index].second.eta1 - ZZpair_masscut[index].second.eta2 ,2) )<<endl;
       
       if ( sqrt( pow( DELTAPHI(  ZZpair_masscut[index].first.phi1 , ZZpair_masscut[index].first.phi2 ),2) + pow( ZZpair_masscut[index].first.eta1 - ZZpair_masscut[index].first.eta2 ,2) ) <= 0.02) continue; //Z1 (lep1,lep2)
       if ( sqrt( pow( DELTAPHI( ZZpair_masscut[index].first.phi1 , ZZpair_masscut[index].second.phi1  ),2) + pow( ZZpair_masscut[index].first.eta1 - ZZpair_masscut[index].second.eta1 ,2) ) <= 0.02)continue; // (Z1lep1,Z2lep1)
       if ( sqrt( pow( DELTAPHI(  ZZpair_masscut[index].first.phi1 , ZZpair_masscut[index].second.phi2  ),2) + pow( ZZpair_masscut[index].first.eta1 - ZZpair_masscut[index].second.eta2 ,2) ) <= 0.02)continue; // (Z1lep1,Z2lep2)
       if ( sqrt( pow( DELTAPHI(  ZZpair_masscut[index].first.phi2 , ZZpair_masscut[index].second.phi1 ),2) + pow( ZZpair_masscut[index].first.eta2 - ZZpair_masscut[index].second.eta1 ,2) ) <= 0.02)continue; // (Z1lep2,Z2lep1)
       if ( sqrt( pow( DELTAPHI(  ZZpair_masscut[index].first.phi2 , ZZpair_masscut[index].second.phi2 ),2) + pow( ZZpair_masscut[index].first.eta2 - ZZpair_masscut[index].second.eta2 ,2) ) <= 0.02)continue; // (Z1lep2,Z2lep2)
       if ( sqrt( pow( DELTAPHI(  ZZpair_masscut[index].second.phi1 ,  ZZpair_masscut[index].second.phi2 ),2) + pow( ZZpair_masscut[index].second.eta1 - ZZpair_masscut[index].second.eta2 ,2) ) <= 0.02)continue; //Z2 (lep1,lep2)
       
       
       ZZpair_deltaRcut.push_back( {ZZpair_masscut[index].first , ZZpair_masscut[index].second} );
       
     }
     
     if( ZZpair_deltaRcut.size()==0) continue;
     
     //// cout<<" ZZpair_deltaRcut = "<< ZZpair_deltaRcut.size()<<endl;
     
     // for(int i=0; i<ZZpair_deltaRcut.size(); i++ ){
       
     //   cout<<"ZZ pairs after deltaR cuts have masses = "<<ZZpair_deltaRcut[i].first.massvalue<<" , "<<ZZpair_deltaRcut[i].second.massvalue<<"and tags = "<<ZZpair_deltaRcut[i].first.tag<<", "<<ZZpair_deltaRcut[i].second.tag<<endl;
       
       
     // }
     
     
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     
 // pT >20 ,10 for any two leptons from ZZ pair candidate
 //=======================================================
     
     vector<std::pair<candidateZ, candidateZ> > ZZpair_pTcut;
     ZZpair_pTcut.clear();
     
     vector<float> leptonspT;
     
     
     for (int index=0; index<ZZpair_deltaRcut.size();index++){
       
       leptonspT.clear();
       
       leptonspT.push_back(ZZpair_deltaRcut[index].first.pt1);
       leptonspT.push_back(ZZpair_deltaRcut[index].first.pt2);
       leptonspT.push_back(ZZpair_deltaRcut[index].second.pt1);
       leptonspT.push_back(ZZpair_deltaRcut[index].second.pt2);
       
       std::sort(leptonspT.rbegin(),leptonspT.rend());
       
       //// cout<<"Z1pt1 = "<<ZZpair_deltaRcut[index].first.pt1<<",Z1pt2 = "<<ZZpair_deltaRcut[index].first.pt2<<" ,Z2pt1 = "<<ZZpair_deltaRcut[index].second.pt1<<", Z2pt2 = "<<ZZpair_deltaRcut[index].second.pt2<<endl;
       
       //// cout<<"leptonspT = "<<leptonspT.at(0)<<" , "<<leptonspT.at(1)<<" , "<<leptonspT.at(2)<<" , "<<leptonspT.at(3)<<" , "<<endl;
       
       if ( leptonspT.at(0) > 20 && leptonspT.at(1) > 10 ) {
	 
	 ZZpair_pTcut.push_back( {ZZpair_deltaRcut[index].first , ZZpair_deltaRcut[index].second} );
       }  
     }
     
     //// cout<<"ZZpair_pTcut = "<<ZZpair_pTcut.size()<<endl;
     
     // for(int i=0; i<ZZpair_pTcut.size(); i++ ){
       
     //   cout<<"ZZ pairs after pT cut have masses = "<<ZZpair_pTcut[i].first.massvalue<<" , "<<ZZpair_pTcut[i].second.massvalue<<"and tags = "<<ZZpair_pTcut[i].first.tag<<", "<<ZZpair_pTcut[i].second.tag<<endl;
       
       
     // }
     
     if ( ZZpair_pTcut.size()==0 )continue;
     
     
     
     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     
     // QCD supression (mll > 4) for ZZ pair
     //======================================
     
     
     double min_mass_2L = 10000.;
     
     
     for (int index=0; index<ZZpair_pTcut.size();index++){
       
       TLorentzVector Lepton1qcd,Lepton2qcd,Lepton3qcd,Lepton4qcd,DiLeptonaQCD,DiLeptonbQCD;
       
       double ZaQCD =-999;
       double ZbQCD =-999;
       
       
       //// cout<<"Z1charge1 = "<<ZZpair_pTcut[index].first.charge1<<" ,Z1charge2 = "<<ZZpair_pTcut[index].first.charge2<<" ,Z2charge1 = "<<ZZpair_pTcut[index].second.charge1<<" ,Z2charge2 = "<<ZZpair_pTcut[index].second.charge2<<endl;
       
       if (ZZpair_pTcut[index].first.charge1 + ZZpair_pTcut[index].second.charge1 == 0 ) //opposite charge
	 {
	   
	   Lepton1qcd.SetPtEtaPhiE(ZZpair_pTcut[index].first.pt1, ZZpair_pTcut[index].first.eta1, ZZpair_pTcut[index].first.phi1, ZZpair_pTcut[index].first.lep1E ); 
	   Lepton3qcd.SetPtEtaPhiE(ZZpair_pTcut[index].second.pt1, ZZpair_pTcut[index].second.eta1, ZZpair_pTcut[index].second.phi1, ZZpair_pTcut[index].second.lep1E );
	   DiLeptonaQCD=Lepton1qcd+Lepton3qcd;
	   
	   Lepton2qcd.SetPtEtaPhiE(ZZpair_pTcut[index].first.pt2, ZZpair_pTcut[index].first.eta2, ZZpair_pTcut[index].first.phi2, ZZpair_pTcut[index].first.lep2E ); 
	   Lepton4qcd.SetPtEtaPhiE(ZZpair_pTcut[index].second.pt2, ZZpair_pTcut[index].second.eta2, ZZpair_pTcut[index].second.phi2, ZZpair_pTcut[index].second.lep2E );
	   DiLeptonbQCD=Lepton2qcd+Lepton4qcd;
	   
	   ZaQCD = DiLeptonaQCD.M();
	   ZbQCD = DiLeptonbQCD.M();
	   
	 }
       
       else {
	 Lepton1qcd.SetPtEtaPhiE(ZZpair_pTcut[index].first.pt1, ZZpair_pTcut[index].first.eta1, ZZpair_pTcut[index].first.phi1, ZZpair_pTcut[index].first.lep1E );
	 Lepton4qcd.SetPtEtaPhiE(ZZpair_pTcut[index].second.pt2, ZZpair_pTcut[index].second.eta2, ZZpair_pTcut[index].second.phi2, ZZpair_pTcut[index].second.lep2E );
	 DiLeptonaQCD=Lepton1qcd+Lepton4qcd;
	 
	 Lepton2qcd.SetPtEtaPhiE(ZZpair_pTcut[index].first.pt2, ZZpair_pTcut[index].first.eta2, ZZpair_pTcut[index].first.phi2, ZZpair_pTcut[index].first.lep2E );
	 Lepton3qcd.SetPtEtaPhiE(ZZpair_pTcut[index].second.pt1, ZZpair_pTcut[index].second.eta1, ZZpair_pTcut[index].second.phi1, ZZpair_pTcut[index].second.lep1E );
	 DiLeptonbQCD=Lepton2qcd+Lepton3qcd;
	 
	 ZaQCD = DiLeptonaQCD.M();
	 ZbQCD = DiLeptonbQCD.M();
	 
   }
       
       //// cout<<"ZaQCD = "<<ZaQCD<<endl;
       //// cout<<"ZbQCD = "<<ZbQCD<<endl;
	  
       if( ZZpair_pTcut[index].first.massvalue < min_mass_2L ) min_mass_2L = ZZpair_pTcut[index].first.massvalue ;
       if( ZZpair_pTcut[index].second.massvalue < min_mass_2L ) min_mass_2L = ZZpair_pTcut[index].second.massvalue ;
       if( ZaQCD < min_mass_2L ) min_mass_2L = ZaQCD ;
       if( ZbQCD < min_mass_2L ) min_mass_2L = ZbQCD ;
       
     }
     
     //// cout<<"min_mass_2L = "<<min_mass_2L<<endl;
     
     if (min_mass_2L <=4)continue;
     
     
     //////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     //smart cut for ZZ pair
     //----------------------
     
     // ZZ objects with alternative pairing (smart cut)
     
     //// cout<<"Smart cut start here"<<endl;
     
     
     vector<std::pair<candidateZ, candidateZ> > ZZpair_smartcut;
     ZZpair_smartcut.clear();
     
     
     for (int index=0; index<ZZpair_pTcut.size();index++){
       
       if ( ZZpair_pTcut[index].first.tag == ZZpair_pTcut[index].second.tag ){ // same ZZpair tag
	 
	 TLorentzVector Lepton1,Lepton2,Lepton3,Lepton4,Lepton1Correction,Lepton2Correction,Lepton3Correction,Lepton4Correction,DiLeptonZa,DiLeptonZb;
	 
	 float massZa=-999.,massZb=-999.;
	 
	 float lepton1ch=-999., lepton2ch=-999.;
	 float lepton3ch=-999., lepton4ch=-999.;
	 
	 lepton1ch= ZZpair_pTcut[index].first.charge1 ;
	 lepton2ch= ZZpair_pTcut[index].first.charge2 ;
	 
	 lepton3ch= ZZpair_pTcut[index].second.charge1 ;
	 lepton4ch= ZZpair_pTcut[index].second.charge2 ;
	 
	 ////cout<<"lepton1ch = "<<lepton1ch<<",lepton2ch = "<<lepton2ch<<" ,lepton3ch"<<lepton3ch<<" ,lepton4ch"<<lepton4ch<<endl;
	 
	 Lepton1.SetPtEtaPhiE(ZZpair_pTcut[index].first.pt1, ZZpair_pTcut[index].first.eta1, ZZpair_pTcut[index].first.phi1, ZZpair_pTcut[index].first.lep1E ); 
	 
	 
	 if (ZZpair_pTcut[index].first.ilept1_FSR){//Z1 lep1 has fsr

	   //// cout<<" Z1 lep1 has fsr "<<endl;
	   //// cout<<"Z1 lep1 fsr PT= "<<RECOPFPHOT_PT[ZZpair_pTcut[index].first.iFSR1]<<" ,and eta = "<<RECOPFPHOT_ETA[ZZpair_pTcut[index].first.iFSR1]<<" ,phi = "<<RECOPFPHOT_PHI[ZZpair_pTcut[index].first.iFSR1]<<endl;
	   
	   Lepton1Correction.SetPtEtaPhiM(RECOPFPHOT_PT[ZZpair_pTcut[index].first.iFSR1],RECOPFPHOT_ETA[ZZpair_pTcut[index].first.iFSR1],RECOPFPHOT_PHI[ZZpair_pTcut[index].first.iFSR1],0);
	   Lepton1+=Lepton1Correction;
	 }
	 
	 
	 Lepton2.SetPtEtaPhiE(ZZpair_pTcut[index].first.pt2, ZZpair_pTcut[index].first.eta2, ZZpair_pTcut[index].first.phi2, ZZpair_pTcut[index].first.lep2E );
	 
	 if (ZZpair_pTcut[index].first.ilept2_FSR){//Z1 lep2 has fsr
	   
	   //// cout<<" Z1 lep2 has fsr = "<<endl;
	   //// cout<<"Z1 lep2 fsr PT= "<<RECOPFPHOT_PT[ZZpair_pTcut[index].first.iFSR2]<<" ,and eta = "<<RECOPFPHOT_ETA[ZZpair_pTcut[index].first.iFSR2]<<" ,phi = "<<RECOPFPHOT_PHI[ZZpair_pTcut[index].first.iFSR2]<<endl;
	   
	   Lepton2Correction.SetPtEtaPhiM(RECOPFPHOT_PT[ZZpair_pTcut[index].first.iFSR2],RECOPFPHOT_ETA[ZZpair_pTcut[index].first.iFSR2],RECOPFPHOT_PHI[ZZpair_pTcut[index].first.iFSR2],0);
	   Lepton2+=Lepton2Correction;
	 }
	 
	 
	 
	 
	 Lepton3.SetPtEtaPhiE(ZZpair_pTcut[index].second.pt1, ZZpair_pTcut[index].second.eta1, ZZpair_pTcut[index].second.phi1, ZZpair_pTcut[index].second.lep1E );
	 
	 if (ZZpair_pTcut[index].second.ilept1_FSR){ //Z2 lep1 has fsr
	   
	   //// cout<<"Z2 lep1 fsr = "<<endl;
	   //// cout<<"Z2 lep1 fsr PT= "<<RECOPFPHOT_PT[ZZpair_pTcut[index].second.iFSR1]<<" ,and eta = "<<RECOPFPHOT_ETA[ZZpair_pTcut[index].second.iFSR1]<<" ,phi = "<<RECOPFPHOT_PHI[ZZpair_pTcut[index].second.iFSR1]<<endl;
	   
	   Lepton3Correction.SetPtEtaPhiM(RECOPFPHOT_PT[ZZpair_pTcut[index].second.iFSR1],RECOPFPHOT_ETA[ZZpair_pTcut[index].second.iFSR1],RECOPFPHOT_PHI[ZZpair_pTcut[index].second.iFSR1],0);
	   Lepton3+=Lepton3Correction;
	 }
	 
	 
	 
	 
	 Lepton4.SetPtEtaPhiE(ZZpair_pTcut[index].second.pt2, ZZpair_pTcut[index].second.eta2, ZZpair_pTcut[index].second.phi2, ZZpair_pTcut[index].second.lep2E );
	 
	 if (ZZpair_pTcut[index].second.ilept2_FSR){ //Z2 lep1 has fsr
	   
	   //// cout<<"Z2 lep2 fsr "<<endl;
	   /////cout<<"Z2 lep2 fsr PT= "<<RECOPFPHOT_PT[ZZpair_pTcut[index].second.iFSR2]<<" ,and eta = "<<RECOPFPHOT_ETA[ZZpair_pTcut[index].second.iFSR2]<<" ,phi = "<<RECOPFPHOT_PHI[ZZpair_pTcut[index].second.iFSR2]<<endl;
	   
	   Lepton4Correction.SetPtEtaPhiM(RECOPFPHOT_PT[ZZpair_pTcut[index].second.iFSR2],RECOPFPHOT_ETA[ZZpair_pTcut[index].second.iFSR2],RECOPFPHOT_PHI[ZZpair_pTcut[index].second.iFSR2],0);
	   Lepton4+=Lepton4Correction;
	 }
	 
	 
	 
	 
	 if ( lepton1ch + lepton3ch == 0 ) //opposite charge
	   
	   {
	     //// cout<<"step lep1 and lep3"<<endl;
	     DiLeptonZa=Lepton1+Lepton3;
	     DiLeptonZb=Lepton2+Lepton4;
	   }
	 else {
	   //// cout<<"step lep1 and lep4"<<endl;
	   DiLeptonZa=Lepton1+Lepton4;
	   DiLeptonZb=Lepton2+Lepton3;
	 }
	 
	 
	 if (fabs(DiLeptonZa.M()-Zmass) < fabs(DiLeptonZb.M()-Zmass)) {
	   massZa=DiLeptonZa.M();
	   massZb=DiLeptonZb.M();
	 }
	 else {
	   massZa=DiLeptonZb.M();
	   massZb=DiLeptonZa.M();
	 }
	 
	 //// cout<<"massZa = "<<massZa<<" ,massZb"<<massZb<<endl;
	 
	 if ( fabs(massZa-Zmass) < fabs( ZZpair_pTcut[index].first.massvalue -Zmass) && massZb<12) continue; // exclude those pairs
	 ////cout << "mass Za and b= " << massZa << " " << massZb << endl;
	 
	 ZZpair_smartcut.push_back( {ZZpair_pTcut[index].first , ZZpair_pTcut[index].second} );
	 
       }
       else { //not same ZZ pair tag
	 
	 ////cout<<"Not applying smart cut on ZZ pair with different tag "<<endl;
	 
	 ZZpair_smartcut.push_back( {ZZpair_pTcut[index].first , ZZpair_pTcut[index].second} );
	 
       }
       
     }
     
     
     //// cout<<" ZZpair_smartcut = "<< ZZpair_smartcut.size()<<endl;
     
     // for(int i=0; i< ZZpair_smartcut.size(); i++ ){
       
     //   cout<<"ZZ pairs after smart cut have masses = "<< ZZpair_smartcut[i].first.massvalue<<" , "<< ZZpair_smartcut[i].second.massvalue<<"and tags = "<< ZZpair_smartcut[i].first.tag<<", "<< ZZpair_smartcut[i].second.tag<<endl;
       
       
     // }
     
     
     if ( ZZpair_smartcut.size()==0 )continue;
     
     
     //////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     
     // M4l > 70 for ZZ pair 
     //=====================
     
     vector<std::pair<candidateZ, candidateZ> > ZZpair_mass4lcut;
     ZZpair_mass4lcut.clear();
     
     for (int index=0; index<ZZpair_smartcut.size();index++){
       
       double mass4l =-999;
       
       mass4l = sqrt( pow( ZZpair_smartcut[index].first.EZ + ZZpair_smartcut[index].second.EZ,2) - pow( ZZpair_smartcut[index].first.pxZ + ZZpair_smartcut[index].second.pxZ ,2) - pow( ZZpair_smartcut[index].first.pyZ + ZZpair_smartcut[index].second.pyZ ,2) - pow( ZZpair_smartcut[index].first.pzZ + ZZpair_smartcut[index].second.pzZ,2) );
       
       //// cout<<"mass4l = "<< mass4l<<endl;
       
       if (std::isnan(mass4l)) {
	 ////cout << "mass4l nan " << endl;
	 continue; 
       }
       
       if (mass4l<=70)continue;
       
       ZZpair_mass4lcut.push_back( {ZZpair_smartcut[index].first , ZZpair_smartcut[index].second} );
       
     }
     
     //// cout<<"ZZpair_mass4lcut = "<<ZZpair_mass4lcut.size()<<endl;
     
     // for(int i=0; i< ZZpair_mass4lcut.size(); i++ ){
       
     //   cout<<"ZZ pairs after mass4l cut have masses = "<< ZZpair_mass4lcut[i].first.massvalue<<" , "<< ZZpair_mass4lcut[i].second.massvalue<<"and tags = "<< ZZpair_mass4lcut[i].first.tag<<", "<< ZZpair_mass4lcut[i].second.tag<<endl;
       
       
     // }
     
     if (ZZpair_mass4lcut.size()==0 )continue;
     
     
     
     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     //choose the best ZZ pair
     //=======================
     
     vector<std::pair<candidateZ, candidateZ> > ZZpair_good;
     ZZpair_good.clear();
     
     Z1cand = {0., -999., -999, -999., -999., false , -999., -999., -999, -999, -999, -999, -999., -999., -999, -999., -999., false, -999., -999., -999., -999., -999., -999., false, -999., -999., -999, -999 , -999};
     
     candidateZ Z2cand = {0., -999., -999, -999., -999., false , -999., -999., -999, -999, -999, -999, -999., -999., -999, -999., -999., false, -999., -999., -999., -999., -999., -999., false, -999., -999., -999, -999 , -999};
     
     Z1 = {0., -999., -999, -999., -999., false , -999., -999., -999, -999, -999, -999, -999., -999., -999, -999., -999., false, -999., -999., -999., -999., -999., -999., false, -999., -999., -999, -999 , -999};
     
     candidateZ Z2 = {0., -999., -999, -999., -999., false , -999., -999., -999, -999, -999, -999, -999., -999., -999, -999., -999., false, -999., -999., -999., -999., -999., -999., false, -999., -999., -999, -999 , -999};
     
     for (int index=0; index< ZZpair_mass4lcut.size();index++){
       
       if( fabs(ZZpair_mass4lcut[index].first.massvalue - Zmass) < fabs(Z1cand.massvalue - Zmass) ){
	 
	 Z1cand =  ZZpair_mass4lcut[index].first;
	 Z2cand =  ZZpair_mass4lcut[index].second;
	 
	 //// cout<<"Z1cand mass = "<<Z1cand.massvalue<<" and Z2cand mass = "<<Z2cand.massvalue<<endl;
	 //// cout<<"SumpT of 2lep from Z2cand = "<< Z2cand.pt1 + Z2cand.pt2 <<endl;
	 
       }
       else if (fabs(ZZpair_mass4lcut[index].first.massvalue - Zmass) == fabs(Z1cand.massvalue - Zmass))
	 
	 {
	   
	   if ( (ZZpair_mass4lcut[index].second.pt1 + ZZpair_mass4lcut[index].second.pt2) > (Z2cand.pt1 + Z2cand.pt2) )
	     
	     Z2cand = ZZpair_mass4lcut[index].second;
	   
	 }	      
     }
     
     //// cout<<" Final Z1cand and Z2 cand: Z1cand mass = "<<Z1cand.massvalue<<" and Z2cand mass = "<<Z2cand.massvalue<<" and Z1cand tag = "<<Z1cand.tag<<" Z2cand tag = "<<Z2cand.tag<<" and Z1 lepton indices = "<<Z1cand.ilept1<<" , "<<Z1cand.ilept2<<" and Z2 lepton indices = "<<Z2cand.ilept1<<" , "<<Z2cand.ilept2<<endl;
     
     
     
     Z1 = Z1cand;
     Z2 = Z2cand;
     
     //// cout<<"Final ZZ pair: Z1 mass = "<<Z1.massvalue<<" and tag "<<Z1.tag<<" and leptons indices = "<<Z1.ilept1<<","<<Z1.ilept2<<" Z2 mass = "<<Z2.massvalue<<" and tag "<<Z2.tag<<" and leptons indices = "<<Z2.ilept1<<","<<Z2.ilept2<<endl;
     
     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     //efficiency weight
     //-------------------
     
     //// cout<<"efficiency weight step after choose ZZ pair "<<endl;
     
     Double_t eff_weight = 1.;
     
     
     z1lept[0] = Z1.ilept1;
     z1lept[1] = Z1.ilept2;
     int z2lept[2] = {Z2.ilept1, Z2.ilept2};
     
     
     if (Z1.tag==1){ //Z1 is mumu
       
       //// cout<<"Z1 is mumu "<<endl;
       
       for(int i = 0; i < 2; i++){
	 
	 Double_t Pt = RECOMU_PT[ z1lept[i] ]; 
	 Double_t Eta = RECOMU_ETA[ z1lept[i] ];
	 
	 //// cout<<"Z1 has muon = "<<z1lept[i]<<" with pt = "<<Pt<<" and eta = "<<Eta<<endl;
	  
	 if( (MC_type == "Spring16" || MC_type == "Moriond17" ) && DATA_type == "NO"){
	   
	   Int_t binx = mu_scale_2016->GetXaxis()->FindBin(Eta);
	   //// cout<<"binx = "<<binx<<endl;
	   Int_t biny = mu_scale_2016->GetYaxis()->FindBin(Pt);
	   ///// cout<<"biny = "<<biny<<endl;
	   //// cout<<"muon scale factor = "<<mu_scale_2016->GetBinContent(binx,biny)<<endl;
	   if (mu_scale_2016->GetBinContent(binx,biny)>0.) eff_weight*=mu_scale_2016->GetBinContent(binx,biny);
	   ///// cout<<"eff_weight = "<<eff_weight<<endl;
	 }
       }
     }
     else if (Z1.tag==2){ //Z1 is ee
       
       //// cout<<"Z1 is ee "<<endl;
       
       
       for(int i = 0; i < 2; i++){
	 
	 Double_t Pt = RECOELE_PT[ z1lept[i] ]; 
	 Double_t Eta = RECOELE_ETA[ z1lept[i] ];
	 
	 //// cout<<"Z1 has electron = "<<z1lept[i]<<" with pt = "<<Pt<<" and eta = "<<Eta<<endl;
	 
	 if( ( MC_type == "Spring16" || MC_type == "Moriond17" ) && DATA_type == "NO"){
	   
	   if(RECOELE_isGap[ z1lept[i] ]==0){
	     Int_t binx = ele_scale_factors2016->GetXaxis()->FindBin(Eta);
	     //// cout<<"binx = "<<binx<<endl;
	     Int_t biny = ele_scale_factors2016->GetYaxis()->FindBin(Pt);
	     ////cout<<"biny = "<<biny<<endl;
	     //// cout<<"ele scale factor2016 = "<<ele_scale_factors2016->GetBinContent(binx,biny)<<endl;
	     if (ele_scale_factors2016->GetBinContent(binx,biny)>0.) eff_weight*=ele_scale_factors2016->GetBinContent(binx,biny);
	     //// cout<<"eff_weight = "<<eff_weight<<endl;
	     
	   }
	      else if(RECOELE_isGap[ z1lept[i] ]==1){
		////cout<<"crack electron"<<endl;
		Int_t binx = ele_scale_factors_gap2016->GetXaxis()->FindBin(Eta);
		////cout<<"binx = "<<binx<<endl;
		Int_t biny = ele_scale_factors_gap2016->GetYaxis()->FindBin(Pt);
		////cout<<"biny = "<<biny<<endl;
		////cout<<"ele scale factor crack2016 = "<<ele_scale_factors_gap2016->GetBinContent(binx,biny)<<endl;
		if (ele_scale_factors_gap2016->GetBinContent(binx,biny)>0.) eff_weight*=ele_scale_factors_gap2016->GetBinContent(binx,biny);
		////cout<<"eff_weight = "<<eff_weight<<endl;
	      }
	   
	 }
       }
     }
     
     ///////////////////////////////
     
     if (Z2.tag==1){ //Z2 is mumu
       
       ////cout<<"Z2 is mumu "<<endl;
       
       for(int i = 0; i < 2; i++){
	 
	 Double_t Pt = RECOMU_PT[ z2lept[i] ]; 
	 Double_t Eta = RECOMU_ETA[ z2lept[i] ];
	 
	 ////cout<<"Z2 has muon = "<<z2lept[i]<<" with pt = "<<Pt<<" and eta = "<<Eta<<endl;
	 
	 if( ( MC_type == "Spring16" || MC_type == "Moriond17" ) && DATA_type == "NO"){
	   
	   Int_t binx = mu_scale_2016->GetXaxis()->FindBin(Eta);
	   ////cout<<"binx = "<<binx<<endl;
	   Int_t biny = mu_scale_2016->GetYaxis()->FindBin(Pt);
	   ////cout<<"biny = "<<biny<<endl;
	   ////cout<<"muon scale factor = "<<mu_scale_2016->GetBinContent(binx,biny)<<endl;
	   if (mu_scale_2016->GetBinContent(binx,biny)>0.) eff_weight*=mu_scale_2016->GetBinContent(binx,biny);
	   ////cout<<"eff_weight = "<<eff_weight<<endl;
	 }
       }
     }
     else if (Z2.tag==2){
       
       ////cout<<"Z2 is ee "<<endl;
       
       for(int i = 0; i < 2; i++){
	 Double_t Pt = RECOELE_PT[ z2lept[i] ]; 
	 Double_t Eta = RECOELE_ETA[ z2lept[i] ];
	 
	 ////cout<<"Z2 has electron = "<<z2lept[i]<<" with pt = "<<Pt<<" and eta = "<<Eta<<endl;
	 
	 if( ( MC_type == "Spring16" || MC_type == "Moriond17" ) && DATA_type == "NO"){
	   
	   if(RECOELE_isGap[ z1lept[i] ]==0){
	     Int_t binx = ele_scale_factors2016->GetXaxis()->FindBin(Eta);
	     ////cout<<"binx = "<<binx<<endl;
	     Int_t biny = ele_scale_factors2016->GetYaxis()->FindBin(Pt);
	     ////cout<<"biny = "<<biny<<endl;
	     ////cout<<"ele scale factor2016 = "<<ele_scale_factors2016->GetBinContent(binx,biny)<<endl;
	     if (ele_scale_factors2016->GetBinContent(binx,biny)>0.) eff_weight*=ele_scale_factors2016->GetBinContent(binx,biny);
	     ////cout<<"eff_weight = "<<eff_weight<<endl;
	   }
	   else if(RECOELE_isGap[ z1lept[i] ]==1){
	     ////cout<<"crack electron"<<endl;
	     Int_t binx = ele_scale_factors_gap2016->GetXaxis()->FindBin(Eta);
	     ////cout<<"binx = "<<binx<<endl;
	     Int_t biny = ele_scale_factors_gap2016->GetYaxis()->FindBin(Pt);
	     ////cout<<"biny = "<<biny<<endl;
	     ////cout<<"ele scale factor crack= "<<ele_scale_factors_gap2016->GetBinContent(binx,biny)<<endl;
	     if (ele_scale_factors_gap2016->GetBinContent(binx,biny)>0.) eff_weight*=ele_scale_factors_gap2016->GetBinContent(binx,biny);
	     ////cout<<"eff_weight = "<<eff_weight<<endl;
	   }
	   
	 }
       }
     }
     
     ////cout<<"DATA_type = "<<DATA_type<<endl;
     
     if (DATA_type == "2016") eff_weight=1.; 
     
     // Changing the weight for pileup and LineShape and efficiency
     if (eff_weight>0.) {
       ////cout<<"eff_weight = "<<eff_weight<<endl;
       newweight=weight*pu_weight*eff_weight;}
     else newweight=weight*pu_weight;
     
     ////cout << "Starting weight + pileup + efficiency= " << newweight << endl;
     // if(debug)
     ////cout << "Efficiency Weight for the 4l: " << eff_weight << " Final weight for ZZ pair = " << newweight << endl;
     
     
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     

     //Question from ARC

     //if (RECO_PFMET < 50 ) continue;
     
     ////////

//======================================================================================================================================================================     
//Insert the extra requirements for VBF used in the analysis
//======================================================================================================================================================================
     
     //Basic cuts to jets AND delta R section
     int njets_pass=0;
     TLorentzVector JET1,JET2,JET3;
     int jet1=-999,jet2=-999,jet3=-999;
     std::vector<int> mu_indexes;
     std::vector<int> e_indexes;
     if(Z1.tag==1){
       mu_indexes.push_back(Z1.ilept1);
       mu_indexes.push_back(Z1.ilept2);
     }else{
       e_indexes.push_back(Z1.ilept1);
       e_indexes.push_back(Z1.ilept2);
     }
     if(Z2.tag==1){
       mu_indexes.push_back(Z2.ilept1);
       mu_indexes.push_back(Z2.ilept2);
     }else{
       e_indexes.push_back(Z2.ilept1);
       e_indexes.push_back(Z2.ilept2);
     }
     std::vector<int> g_indexes;
     if(Z1.withFSR){
       g_indexes.push_back(Z1.iFSR1);
       g_indexes.push_back(Z1.iFSR2);
     }
     if(Z2.withFSR){
       g_indexes.push_back(Z2.iFSR1);
       g_indexes.push_back(Z2.iFSR2);
     }

     //std::cout<<"Jet selection ---- leptons so far"<<std::endl;
     //std::cout<<Form("l1: %.2f, l2: %.2f, l3: %.2f, l4: %.2f",Z1.pt1,Z1.pt2,Z2.pt1,Z2.pt2)<<std::endl;
     for(int i=0;i<RECO_PFJET_N;i++){
       bool jetfail = false;
       if(RECO_PFJET_PT[i]>30. && fabs(RECO_PFJET_ETA[i])<4.7){
      	 //Check that jet has deltaR>0.4 away from any tight lepton corrected for FSR
	 for(int mu = 0; mu < (int)mu_indexes.size(); ++mu){
	   if (fabs(RECOMU_SIP[mu_indexes[mu]])>=4.) continue;
      	   if (RECOMU_PFX_dB_new[mu_indexes[mu]]>=0.35) continue;
	   double deltaR = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOMU_PHI[mu_indexes[mu]]),2) + pow(RECO_PFJET_ETA[i]-RECOMU_ETA[mu_indexes[mu]],2) );
	   //cout << "1st lepton muon: " << " pT=" << RECOMU_PT[mu_indexes[mu]] <<" deltaR "<< deltaR <<endl;
	   if (deltaR<0.4){
	     jetfail = true;
     	     //cout << " jetfail " << jetfail <<endl;
	     break;
     	   }
     	 }
	 
      	 for(int ele = 0; ele < (int)e_indexes.size(); ++ele){
      	   if (fabs(RECOELE_SIP[e_indexes[ele]])>=4.) continue;
	   if (RECOELE_PFX_rho_new[e_indexes[ele]]>=0.35) continue;
      	   double deltaR = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOELE_PHI[e_indexes[ele]]),2) + pow(RECO_PFJET_ETA[i] - RECOELE_ETA[e_indexes[ele]],2));
     	   //cout << "1st lepton electron: " << " pT=" << RECOELE_PT[e_indexes[ele]] <<" deltaR "<< deltaR <<endl;
	   if (deltaR<0.4){
     	     jetfail = true;
     	     //cout << " jetfail " << jetfail <<endl;
	     break;
     	   }
     	 }

	 // cleaning w.r.t FSR photons attached to leptons
	 for(int j=0.;j<(int)g_indexes.size();j++){
	    if(g_indexes[j]!=-1){
	      double deltaR = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOPFPHOT_PHI[g_indexes[j]]),2) + pow(RECO_PFJET_ETA[i] - RECOPFPHOT_ETA[g_indexes[j]],2));
	      if (deltaR<0.4){
		jetfail = true;
		//cout << " jetfail " << jetfail <<endl;
		break;
	      }
	    }
	 }
	 
	 if (!jetfail){
	   //cout<< " PASS jet " <<i<<" PT= "<<RECO_PFJET_PT[i]<<" ETA= "<<RECO_PFJET_ETA[i]<<" PUID= "<<RECO_PFJET_PUID[i]<<endl;
	   njets_pass++;
	   if (njets_pass==1){
	     jet1=i;
	     JET1.SetPtEtaPhiE(RECO_PFJET_PT[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
	     f_jets_highpt_pt[0] = RECO_PFJET_PT[i];
	     f_jets_highpt_pt_error[0] = RECO_PFJET_PT[i]-RECO_PFJET_PT_UncDn[i];
	     f_jets_highpt_eta[0] = RECO_PFJET_ETA[i];
	     f_jets_highpt_phi[0] = RECO_PFJET_PHI[i];
	     f_jets_highpt_et[0] = RECO_PFJET_ET[i];
	   }
	   if (njets_pass==2){
	     jet2=i;
	     JET2.SetPtEtaPhiE(RECO_PFJET_PT[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
	     f_jets_highpt_pt[1] = RECO_PFJET_PT[i];
	     f_jets_highpt_pt_error[1] = RECO_PFJET_PT[i]-RECO_PFJET_PT_UncDn[i];
	     f_jets_highpt_eta[1] = RECO_PFJET_ETA[i];
	     f_jets_highpt_phi[1] = RECO_PFJET_PHI[i];
	     f_jets_highpt_et[1] = RECO_PFJET_ET[i];
	   }
	   if (njets_pass==3){
	     jet3=i;
	     JET3.SetPtEtaPhiE(RECO_PFJET_PT[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
	     f_jets_highpt_pt[2] = RECO_PFJET_PT[i];
	     f_jets_highpt_pt_error[2] = RECO_PFJET_PT[i]-RECO_PFJET_PT_UncDn[i];
	     f_jets_highpt_eta[2] = RECO_PFJET_ETA[i];
	     f_jets_highpt_phi[2] = RECO_PFJET_PHI[i];
	     f_jets_highpt_et[2] = RECO_PFJET_ET[i];
	   }
        }
        
       }
       else{
      	 jetfail = true;
       }
     }
     //Number of jets and mJJ,delta eta cuts // categories
     f_njets_pass = njets_pass;
     //cout << "Number of jets passing the selection is = " << njets_pass << endl;
     

     // b-tagged jets - ordered in pT
     int n_bjets=0;
     int index_bjets[2]={-999,-999};
     
     for (int i=0;i<50;i++){
       //if (cSV_BTagJet_DISCR[i] > 0.460){ // Loose
       if (cSV_BTagJet_DISCR[i] > 0.8484){ // Medium
	 //if(cSV_BTagJet_PT[i]>30. && fabs(cSV_BTagJet_ETA[i])<4.7 ) cout << "Found a bjet (pT>30 and |eta|<2.4) with pT= " << cSV_BTagJet_PT[i] << endl;	 
	 n_bjets++;
	 if (n_bjets==1) index_bjets[0]=i; 
	 if (n_bjets==2) index_bjets[1]=i;	   
       }
     }
     f_Nbjets = n_bjets;
     //cout << "Number of b-jets passing the selection= " << n_bjets << endl;
     
     
     ///-------------- Fills MET properties ---------------------
     f_pfmet=RECO_PFMET;
     f_pfmet_JetEnUp = RECO_PFMET_JetEnUp;
     f_pfmet_JetEnDn = RECO_PFMET_JetEnDn;
     f_pfmet_ElectronEnUp = RECO_PFMET_ElectronEnUp;
     f_pfmet_ElectronEnDn = RECO_PFMET_ElectronEnDn;
     f_pfmet_MuonEnUp = RECO_PFMET_MuonEnUp;
     f_pfmet_MuonEnDn = RECO_PFMET_MuonEnDn;
     f_pfmet_JetResUp = RECO_PFMET_JetResUp;
     f_pfmet_JetResDn = RECO_PFMET_JetResDn;
     f_pfmet_UnclusteredEnUp = RECO_PFMET_UnclusteredEnUp;
     f_pfmet_UnclusteredEnDn = RECO_PFMET_UnclusteredEnDn;
     f_pfmet_PhotonEnUp = RECO_PFMET_PhotonEnUp;
     f_pfmet_PhotonEnDn = RECO_PFMET_PhotonEnDn;
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------     


     //=================
     // Build CR Regions //
     //=================
     
     //Order the leptons by decreasing pT
     float z1pt1err=-999, z1pt2err=-999, z2pt1err=-999, z2pt2err=-999;
     float z1pdg1=-999, z1pdg2=-999, z2pdg1=-999, z2pdg2=-999;
     float z1Iso1=-999, z1Iso2=-999, z2Iso1=-999, z2Iso2=-999;
     float z1PF1=-999, z1PF2=-999, z2PF1=-999, z2PF2=-999;
     TLorentzVector L11P4, L12P4, L21P4, L22P4;
     int L11PID, L12PID, L21PID, L22PID;
     
     if(Z1.tag == 2){//double electron
       L11PID = (RECOELE_CHARGE[Z1.ilept1] == -1)? 11:-11;
       L12PID = (RECOELE_CHARGE[Z1.ilept2] == -1)? 11:-11;
       L11P4.SetPtEtaPhiM(RECOELE_PT[Z1.ilept1],RECOELE_ETA[Z1.ilept1],RECOELE_PHI[Z1.ilept1],0.000511);
       L12P4.SetPtEtaPhiM(RECOELE_PT[Z1.ilept2],RECOELE_ETA[Z1.ilept2],RECOELE_PHI[Z1.ilept2],0.000511);
       
       z1pt1err = RECOELE_PTError[Z1.ilept1];
       z1pt2err = RECOELE_PTError[Z1.ilept2];
       z1pdg1 = (Z1.charge1 < 0)? 11:-11;
       z1pdg2 = (Z1.charge2 < 0)? 11:-11;
       z1Iso1 = Z1.isol1;
       z1Iso2 = Z1.isol2;
       
       bool BDT_ok_1 = 0; // Spring16 with CMSSW_8_2_0 
       if( RECOELE_PT[Z1.ilept1] > 7. &&  RECOELE_PT[Z1.ilept1] <= 10. ){
	 if( fabs(RECOELE_scl_Eta[Z1.ilept1]) < .8 && RECOELE_mvaNonTrigV0[Z1.ilept1] > -0.211 ) BDT_ok_1 = 1 ;
	 if( ( fabs(RECOELE_scl_Eta[Z1.ilept1]) >= .8 && fabs(RECOELE_scl_Eta[Z1.ilept1]) < 1.479 )
	    && RECOELE_mvaNonTrigV0[Z1.ilept1] > -0.396 ) BDT_ok_1 = 1 ;
	 if( fabs(RECOELE_scl_Eta[Z1.ilept1]) >= 1.479 && RECOELE_mvaNonTrigV0[Z1.ilept1] > -0.215 ) BDT_ok_1 = 1 ;
       }
       else { 
	 if( fabs(RECOELE_scl_Eta[Z1.ilept1]) < .8 && RECOELE_mvaNonTrigV0[Z1.ilept1] > -0.870 ) BDT_ok_1 = 1 ;
	 if( ( fabs(RECOELE_scl_Eta[Z1.ilept1]) >= .8 && fabs(RECOELE_scl_Eta[Z1.ilept1]) <= 1.479 )
	    && RECOELE_mvaNonTrigV0[Z1.ilept1] > -0.838 ) BDT_ok_1 = 1 ;
	 if( fabs(RECOELE_scl_Eta[Z1.ilept1]) > 1.479 && RECOELE_mvaNonTrigV0[Z1.ilept1] > -0.763 ) BDT_ok_1 = 1 ;
       } 
       bool BDT_ok_2 = 0; // Spring16 with CMSSW_8_2_0
       if( RECOELE_PT[Z1.ilept2] > 7. &&  RECOELE_PT[Z1.ilept2] <= 10. ){
	 if( fabs(RECOELE_scl_Eta[Z1.ilept2]) < .8 && RECOELE_mvaNonTrigV0[Z1.ilept2] > -0.211 ) BDT_ok_2 = 1 ;
	 if( ( fabs(RECOELE_scl_Eta[Z1.ilept2]) >= .8 && fabs(RECOELE_scl_Eta[Z1.ilept2]) < 1.479 )
	    && RECOELE_mvaNonTrigV0[Z1.ilept2] > -0.396 ) BDT_ok_2 = 1 ;
	 if( fabs(RECOELE_scl_Eta[Z1.ilept2]) >= 1.479 && RECOELE_mvaNonTrigV0[Z1.ilept2] > -0.215 ) BDT_ok_2 = 1 ;
       }
       else {
	 if( fabs(RECOELE_scl_Eta[Z1.ilept2]) < .8 && RECOELE_mvaNonTrigV0[Z1.ilept2] > -0.870 ) BDT_ok_2 = 1 ;
	 if( ( fabs(RECOELE_scl_Eta[Z1.ilept2]) >= .8 && fabs(RECOELE_scl_Eta[Z1.ilept2]) <= 1.479 )
	    && RECOELE_mvaNonTrigV0[Z1.ilept2] > -0.838 ) BDT_ok_2 = 1 ;
	 if( fabs(RECOELE_scl_Eta[Z1.ilept2]) > 1.479 && RECOELE_mvaNonTrigV0[Z1.ilept2] > -0.763 ) BDT_ok_2 = 1 ;
       }
       z1PF1= (BDT_ok_1)? 1:0;
       z1PF2= (BDT_ok_2)? 1:0;       
     }else{//double muon
       L11PID = (RECOMU_CHARGE[Z1.ilept1] == -1)? 13:-13;
       L12PID = (RECOMU_CHARGE[Z1.ilept2] == -1)? 13:-13;
       L11P4.SetPtEtaPhiM(RECOMU_PT[Z1.ilept1],RECOMU_ETA[Z1.ilept1],RECOMU_PHI[Z1.ilept1],0.000511);
       L12P4.SetPtEtaPhiM(RECOMU_PT[Z1.ilept2],RECOMU_ETA[Z1.ilept2],RECOMU_PHI[Z1.ilept2],0.000511);
       
       z1pt1err = RECOMU_mubesttrkPTError[Z1.ilept1];
       z1pt2err = RECOMU_mubesttrkPTError[Z1.ilept2];
       z1pdg1 = (Z1.charge1 < 0)? 13:-13;
       z1pdg2 = (Z1.charge2 < 0)? 13:-13;       
       z1Iso1 = Z1.isol1;
       z1PF1  = RECOMU_isPFMu[Z1.ilept1];
       z1Iso2 = Z1.isol2;
       z1PF2  = RECOMU_isPFMu[Z1.ilept2];
     }
     if(Z2.tag == 2){//double elecrton
       L21PID = (RECOELE_CHARGE[Z2.ilept1] == -1)? 11:-11;
       L22PID = (RECOELE_CHARGE[Z2.ilept2] == -1)? 11:-11;
       L21P4.SetPtEtaPhiM(RECOELE_PT[Z2.ilept1],RECOELE_ETA[Z2.ilept1],RECOELE_PHI[Z2.ilept1],0.000511);
       L22P4.SetPtEtaPhiM(RECOELE_PT[Z2.ilept2],RECOELE_ETA[Z2.ilept2],RECOELE_PHI[Z2.ilept2],0.000511);

       z2pt1err = RECOELE_PTError[Z2.ilept1];
       z2pt2err = RECOELE_PTError[Z2.ilept2];       
       z2pdg1 = (Z2.charge1 < 0)? 11:-11;
       z2pdg2 = (Z2.charge2 < 0)? 11:-11;
       z2Iso1 = Z2.isol1;
       z2Iso2 = Z2.isol2;
       
       bool BDT_ok_1 = 0; // Spring16 with CMSSW_8_2_0 
       if( RECOELE_PT[Z2.ilept1] > 7. &&  RECOELE_PT[Z2.ilept1] <= 10. ){
	 if( fabs(RECOELE_scl_Eta[Z2.ilept1]) < .8 && RECOELE_mvaNonTrigV0[Z2.ilept1] > -0.211 ) BDT_ok_1 = 1 ;
	 if( ( fabs(RECOELE_scl_Eta[Z2.ilept1]) >= .8 && fabs(RECOELE_scl_Eta[Z2.ilept1]) < 1.479 )
	    && RECOELE_mvaNonTrigV0[Z2.ilept1] > -0.396 ) BDT_ok_1 = 1 ;
	 if( fabs(RECOELE_scl_Eta[Z2.ilept1]) >= 1.479 && RECOELE_mvaNonTrigV0[Z2.ilept1] > -0.215 ) BDT_ok_1 = 1 ;
       }
       else { 
	 if( fabs(RECOELE_scl_Eta[Z2.ilept1]) < .8 && RECOELE_mvaNonTrigV0[Z2.ilept1] > -0.870 ) BDT_ok_1 = 1 ;
	 if( ( fabs(RECOELE_scl_Eta[Z2.ilept1]) >= .8 && fabs(RECOELE_scl_Eta[Z2.ilept1]) <= 1.479 )
	    && RECOELE_mvaNonTrigV0[Z2.ilept1] > -0.838 ) BDT_ok_1 = 1 ;
	 if( fabs(RECOELE_scl_Eta[Z2.ilept1]) > 1.479 && RECOELE_mvaNonTrigV0[Z2.ilept1] > -0.763 ) BDT_ok_1 = 1 ;
       }  
       bool BDT_ok_2 = 0; // Spring16 with CMSSW_8_2_0
       if( RECOELE_PT[Z2.ilept2] > 7. &&  RECOELE_PT[Z2.ilept2] <= 10. ){
	 if( fabs(RECOELE_scl_Eta[Z2.ilept2]) < .8 && RECOELE_mvaNonTrigV0[Z2.ilept2] > -0.211 ) BDT_ok_2 = 1 ;
	 if( ( fabs(RECOELE_scl_Eta[Z2.ilept2]) >= .8 && fabs(RECOELE_scl_Eta[Z2.ilept2]) < 1.479 )
	    && RECOELE_mvaNonTrigV0[Z2.ilept2] > -0.396 ) BDT_ok_2 = 1 ;
	 if( fabs(RECOELE_scl_Eta[Z2.ilept2]) >= 1.479 && RECOELE_mvaNonTrigV0[Z2.ilept2] > -0.215 ) BDT_ok_2 = 1 ;
       }
       else {
	 if( fabs(RECOELE_scl_Eta[Z2.ilept2]) < .8 && RECOELE_mvaNonTrigV0[Z2.ilept2] > -0.870 ) BDT_ok_2 = 1 ;
	 if( ( fabs(RECOELE_scl_Eta[Z2.ilept2]) >= .8 && fabs(RECOELE_scl_Eta[Z2.ilept2]) <= 1.479 )
	    && RECOELE_mvaNonTrigV0[Z2.ilept2] > -0.838 ) BDT_ok_2 = 1 ;
	 if( fabs(RECOELE_scl_Eta[Z2.ilept2]) > 1.479 && RECOELE_mvaNonTrigV0[Z2.ilept2] > -0.763 ) BDT_ok_2 = 1 ;
       }
       z2PF1= (BDT_ok_1)? 1:0;
       z2PF2= (BDT_ok_2)? 1:0;       
     }else{//double muon
       L21PID = (RECOMU_CHARGE[Z2.ilept1] == -1)? 13:-13;
       L22PID = (RECOMU_CHARGE[Z2.ilept2] == -1)? 13:-13;
       L21P4.SetPtEtaPhiM(RECOMU_PT[Z2.ilept1],RECOMU_ETA[Z2.ilept1],RECOMU_PHI[Z2.ilept1],0.105);
       L22P4.SetPtEtaPhiM(RECOMU_PT[Z2.ilept2],RECOMU_ETA[Z2.ilept2],RECOMU_PHI[Z2.ilept2],0.105);
       
       z2pt1err = RECOMU_mubesttrkPTError[Z2.ilept1];
       z2pt2err = RECOMU_mubesttrkPTError[Z2.ilept2];       
       z2pdg1 = (Z2.charge1 < 0)? 13:-13;
       z2pdg2 = (Z2.charge2 < 0)? 13:-13;       
       z2Iso1 = Z2.isol1;
       z2PF1  = RECOMU_isPFMu[Z2.ilept1];
       z2Iso2 = Z2.isol2;
       z2PF2  = RECOMU_isPFMu[Z2.ilept2];
     }
     
     float lepts[4] = {Z1.pt1,Z1.pt2,Z2.pt1,Z2.pt2};
     std::vector< std::vector<float> > vlepts;
     vlepts.push_back({Z1.pt1,Z1.eta1,Z1.phi1,z1pdg1,z1Iso1,z1PF1,z1pt1err});
     vlepts.push_back({Z1.pt2,Z1.eta2,Z1.phi2,z1pdg2,z1Iso2,z1PF2,z1pt2err});
     vlepts.push_back({Z2.pt1,Z2.eta1,Z2.phi1,z2pdg1,z2Iso1,z2PF1,z2pt1err});
     vlepts.push_back({Z2.pt2,Z2.eta2,Z2.phi2,z2pdg2,z2Iso2,z2PF2,z2pt2err});
     int sort_indexes[4] = {-1,-1,-1,-1};
     for(int l1=0; l1<4; ++l1){
       float maxpt = 0;
       for(int l2=0; l2<4; ++l2){
	 if(lepts[l2] > maxpt){
	   maxpt = lepts[l2];
	   sort_indexes[l1] = l2;
	 }
       }
       lepts[sort_indexes[l1]] = -1;
     }
     //cout<<"Sorted leptons: "<<Form("%i, %i, %i, %i",sort_indexes[0],sort_indexes[1],sort_indexes[2],sort_indexes[3])<<endl;
     //cout<<"Sorted leptons: "<<Form("%.2f, %.2f, %.2f, %.2f",vlepts[sort_indexes[0]][0],vlepts[sort_indexes[1]][0],vlepts[sort_indexes[2]][0],vlepts[sort_indexes[3]][0])<<endl;
     
     int PF_1=100, PF_2=100;
     double Iso_1=10000., Iso_2=10000.;

     int PF_1_Step2=100, PF_2_Step2=100;
     double Iso_1_Step2=10000., Iso_2_Step2=10000.;     

     
     double ISO_CUT = 0.35;
     
     
     // ZZ pair at final
     TLorentzVector hP4,Z1P4,Z2P4;
     
     Z1P4.SetPxPyPzE(Z1.pxZ,Z1.pyZ,Z1.pzZ,Z1.EZ);
     Z2P4.SetPxPyPzE(Z2.pxZ,Z2.pyZ,Z2.pzZ,Z2.EZ);
     
     hP4 = Z1P4 + Z2P4;
     
     double final_mass4l = hP4.M();
     double final_pT4l = hP4.Pt();
     
     double massofhiggs=hP4.M();
     double massofZ1=Z1P4.M();
     double massofZ2=Z2P4.M();
          
     //// cout <<"Z2 has 2 leptons with iso = "<<Z2.isol1<<" , "<<Z2.isol2<<endl;
     
     Iso_1 = Z2.isol1;
     Iso_2 = Z2.isol2;
     
     //cout<<"Iso_1 = "<<Iso_1<<" and Iso_2 = "<<Iso_2<<endl;
     
     if ( Z1.tag==1 && Z2.tag==1 ) //Z1->mumu & Z2->mumu
       
       {
	 //cout<<"Z1 --> u u && Z2 --> u u"<<endl;
	 
	 PF_1=RECOMU_isPFMu[Z2.ilept1];
	 PF_2=RECOMU_isPFMu[Z2.ilept2];
	 
	 //cout<<"Event MMMM has Z2 with 2 muons with properties "<< "PF_1 = "<<PF_1 <<" and Iso_1 = "<<Iso_1<<"and PF_2 = "<<PF_2<<"and Iso_2 = "<<Iso_2 <<endl;
	 
	 if((Iso_1 >= 0.35 || PF_1==0) && (Iso_2 >= 0.35 || PF_2==0)){ //2p2f
	   //cout<<">>>> 4mu 2P2F <<<<"<<endl;
	   //---- Set tree entries ----
	   f_run2 = Run;
	   f_lumi2 = LumiSection;
	   f_event2 = Event;
	   f_weight2 = newweight;
	   f_2p2f = 1;
	   f_lept1_pt = vlepts[sort_indexes[0]][0];
	   f_lept1_eta = vlepts[sort_indexes[0]][1];
	   f_lept1_phi = vlepts[sort_indexes[0]][2];
	   f_lept1_pdgid = vlepts[sort_indexes[0]][3];
	   f_lept2_pt = vlepts[sort_indexes[1]][0];
	   f_lept2_eta = vlepts[sort_indexes[1]][1];
	   f_lept2_phi = vlepts[sort_indexes[1]][2];
	   f_lept2_pdgid = vlepts[sort_indexes[1]][3];
	   f_lept3_pt = vlepts[sort_indexes[2]][0];
	   f_lept3_eta = vlepts[sort_indexes[2]][1];
	   f_lept3_phi = vlepts[sort_indexes[2]][2];
	   f_lept3_pdgid = vlepts[sort_indexes[2]][3];
	   f_lept4_pt = vlepts[sort_indexes[3]][0];
	   f_lept4_eta = vlepts[sort_indexes[3]][1];
	   f_lept4_phi = vlepts[sort_indexes[3]][2];
	   f_lept4_pdgid = vlepts[sort_indexes[3]][3];
	   f_lept1_pass = (vlepts[sort_indexes[0]][4]>=0.35 || vlepts[sort_indexes[0]][5]==0)? 0:1;
	   f_lept2_pass = (vlepts[sort_indexes[1]][4]>=0.35 || vlepts[sort_indexes[1]][5]==0)? 0:1;
	   f_lept3_pass = (vlepts[sort_indexes[2]][4]>=0.35 || vlepts[sort_indexes[2]][5]==0)? 0:1;
	   f_lept4_pass = (vlepts[sort_indexes[3]][4]>=0.35 || vlepts[sort_indexes[3]][5]==0)? 0:1;
	   f_lept1_pt_error = vlepts[sort_indexes[0]][6];
	   f_lept2_pt_error = vlepts[sort_indexes[1]][6];
	   f_lept3_pt_error = vlepts[sort_indexes[2]][6];
	   f_lept4_pt_error = vlepts[sort_indexes[3]][6];

	   
	   R_AI_AI_MMMM++;
	   R_AI_AI_MMMM_w=R_AI_AI_MMMM_w+newweight;
	   
	   hM4l_MMMM_AI_AI->Fill(final_mass4l,newweight);
	   hPFMET_MMMM_AI_AI->Fill(RECO_PFMET,newweight);


           //for Simran
	   
           hPt_4mu_AI_AI_1->Fill(final_pT4l,newweight);
           hPt_4mu_AI_AI_2->Fill(final_pT4l,newweight);
	   hPt_4mu_AI_AI_3->Fill(final_pT4l,newweight);
	   
	   //cout<<"MMMM_AI_AI=== "<<final_mass4l<<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<<final_pT4l<<" "<< weight <<" 0.35"<<endl;
	   //cout<<"MMMM_AI_AI=== "<<final_mass4l<<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   // sprintf (outformat_MMMM,"MMMM_AI_AI=%.2f:%.2f:%.2f:%.2f:%.2f:%.2f",final_mass4l,Z2.pt1,Z2.eta1,Z2.pt2,Z2.eta2,ISO_CUT);
	   
	   ////cout<<" PF_1 = "<< PF_1<<" ,  PF_2 = "<< PF_2<<endl;
	   ////cout<<" Iso_1 = "<< Iso_1<<" ,  Iso_2 = "<< Iso_2<<endl;
	   
	   //sprintf (Eventformat,"MMMM_AI_AI=%d:%d:%d:%.2f:%.2f:%.2f",Run,LumiSection,Event,massofhiggs,massofZ1,massofZ2);

	   if ( abs(final_mass4l-Zmass)<=10. ){

	     hPFMET_MMMM_AI_AI_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_MMMM_AI_AI_Zpeak_log->Fill(RECO_PFMET,newweight);

	     //cout<<"MMMM_AI_AI_Zpeak=== "<<final_mass4l<<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   }
	   
	   
	 }//end AI_AI
	 
	 if((Iso_1 >= 0.35 || PF_1==0) && (Iso_2 < 0.35 && PF_2==1)){ //3p1f
	   //cout<<">>>> 4mu 3P1F - 1st lepton fails <<<<"<<endl;
	   //---- Set tree entries ----
	   f_run2 = Run;
	   f_lumi2 = LumiSection;
	   f_event2 = Event;
	   f_weight2 = newweight;
	   f_3p1f = 1;
	   f_lept1_pt = vlepts[sort_indexes[0]][0];
	   f_lept1_eta = vlepts[sort_indexes[0]][1];
	   f_lept1_phi = vlepts[sort_indexes[0]][2];
	   f_lept1_pdgid = vlepts[sort_indexes[0]][3];
	   f_lept2_pt = vlepts[sort_indexes[1]][0];
	   f_lept2_eta = vlepts[sort_indexes[1]][1];
	   f_lept2_phi = vlepts[sort_indexes[1]][2];
	   f_lept2_pdgid = vlepts[sort_indexes[1]][3];
	   f_lept3_pt = vlepts[sort_indexes[2]][0];
	   f_lept3_eta = vlepts[sort_indexes[2]][1];
	   f_lept3_phi = vlepts[sort_indexes[2]][2];
	   f_lept3_pdgid = vlepts[sort_indexes[2]][3];
	   f_lept4_pt = vlepts[sort_indexes[3]][0];
	   f_lept4_eta = vlepts[sort_indexes[3]][1];
	   f_lept4_phi = vlepts[sort_indexes[3]][2];
	   f_lept4_pdgid = vlepts[sort_indexes[3]][3];
	   f_lept1_pass = (vlepts[sort_indexes[0]][4]>=0.35 || vlepts[sort_indexes[0]][5]==0)? 0:1;
	   f_lept2_pass = (vlepts[sort_indexes[1]][4]>=0.35 || vlepts[sort_indexes[1]][5]==0)? 0:1;
	   f_lept3_pass = (vlepts[sort_indexes[2]][4]>=0.35 || vlepts[sort_indexes[2]][5]==0)? 0:1;
	   f_lept4_pass = (vlepts[sort_indexes[3]][4]>=0.35 || vlepts[sort_indexes[3]][5]==0)? 0:1;
	   f_lept1_pt_error = vlepts[sort_indexes[0]][6];
	   f_lept2_pt_error = vlepts[sort_indexes[1]][6];
	   f_lept3_pt_error = vlepts[sort_indexes[2]][6];
	   f_lept4_pt_error = vlepts[sort_indexes[3]][6];
	   
	   
	   R_AI_I_MMMM++;
	   R_AI_I_MMMM_w=R_AI_I_MMMM_w+newweight;
	   
	   hM4l_MMMM_AI_I->Fill(final_mass4l,newweight);	   
	   hM4l_MMMM_3p1f_total->Fill(final_mass4l,newweight);
	   
	   hPFMET_MMMM_AI_I->Fill(RECO_PFMET,newweight);
	   hPFMET_MMMM_3p1f_total->Fill(RECO_PFMET,newweight);
	   
	   //cout<<"MMMM_AI_I=== "<< final_mass4l <<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<<final_pT4l<<" "<< weight <<" 0.35"<<endl;
	   //cout<<"MMMM_AI_I=== "<< final_mass4l <<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   //sprintf (outformat_MMMM,"MMMM_AI_I=%.2f:%.2f:%.2f:%.2f:%.2f:%.2f",final_mass4l,Z2.pt1,Z2.eta1,Z2.pt2,Z2.eta2,ISO_CUT);
	   
	   //sprintf (Eventformat,"MMMM_AI_I=%d:%d:%d:%.2f:%.2f:%.2f",Run,LumiSection,Event,massofhiggs,massofZ1,massofZ2);

	   if ( abs(final_mass4l-Zmass)<=10.){
	     
	     hPFMET_MMMM_AI_I_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_MMMM_AI_I_Zpeak_log->Fill(RECO_PFMET,newweight);

	     hPFMET_MMMM_3p1f_total_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_MMMM_3p1f_total_Zpeak_log->Fill(RECO_PFMET,newweight);
	     
	     ////cout<<"MMMM_AI_I_Zpeak=== "<<final_mass4l<<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   }
	   
	 }//end AI_I
	 
	 if((Iso_1 < 0.35 && PF_1==1) && (Iso_2 >= 0.35 || PF_2==0)){//3p1f
	   //cout<<">>>> 4mu 3P1F - 2nd lepton fails <<<<"<<endl;
	   //---- Set tree entries ----
	   f_run2 = Run;
	   f_lumi2 = LumiSection;
	   f_event2 = Event;
	   f_weight2 = newweight;
	   f_3p1f = 1;
	   f_lept1_pt = vlepts[sort_indexes[0]][0];
	   f_lept1_eta = vlepts[sort_indexes[0]][1];
	   f_lept1_phi = vlepts[sort_indexes[0]][2];
	   f_lept1_pdgid = vlepts[sort_indexes[0]][3];
	   f_lept2_pt = vlepts[sort_indexes[1]][0];
	   f_lept2_eta = vlepts[sort_indexes[1]][1];
	   f_lept2_phi = vlepts[sort_indexes[1]][2];
	   f_lept2_pdgid = vlepts[sort_indexes[1]][3];
	   f_lept3_pt = vlepts[sort_indexes[2]][0];
	   f_lept3_eta = vlepts[sort_indexes[2]][1];
	   f_lept3_phi = vlepts[sort_indexes[2]][2];
	   f_lept3_pdgid = vlepts[sort_indexes[2]][3];
	   f_lept4_pt = vlepts[sort_indexes[3]][0];
	   f_lept4_eta = vlepts[sort_indexes[3]][1];
	   f_lept4_phi = vlepts[sort_indexes[3]][2];
	   f_lept4_pdgid = vlepts[sort_indexes[3]][3];
	   f_lept1_pass = (vlepts[sort_indexes[0]][4]>=0.35 || vlepts[sort_indexes[0]][5]==0)? 0:1;
	   f_lept2_pass = (vlepts[sort_indexes[1]][4]>=0.35 || vlepts[sort_indexes[1]][5]==0)? 0:1;
	   f_lept3_pass = (vlepts[sort_indexes[2]][4]>=0.35 || vlepts[sort_indexes[2]][5]==0)? 0:1;
	   f_lept4_pass = (vlepts[sort_indexes[3]][4]>=0.35 || vlepts[sort_indexes[3]][5]==0)? 0:1;
	   f_lept1_pt_error = vlepts[sort_indexes[0]][6];
	   f_lept2_pt_error = vlepts[sort_indexes[1]][6];
	   f_lept3_pt_error = vlepts[sort_indexes[2]][6];
	   f_lept4_pt_error = vlepts[sort_indexes[3]][6];
	   
	   
	   R_I_AI_MMMM++;
	   R_I_AI_MMMM_w=R_I_AI_MMMM_w+newweight;
	   
	   hM4l_MMMM_I_AI->Fill(final_mass4l,newweight);	   
	   hM4l_MMMM_3p1f_total->Fill(final_mass4l,newweight);

	   hPFMET_MMMM_I_AI->Fill(RECO_PFMET,newweight);
	   hPFMET_MMMM_3p1f_total->Fill(RECO_PFMET,newweight);
	   
	   //cout<<"MMMM_I_AI=== "<< final_mass4l <<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<<final_pT4l<<" "<< weight <<" 0.35"<<endl;
	   //cout<<"MMMM_I_AI=== "<< final_mass4l <<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   //sprintf (outformat_MMMM,"MMMM_I_AI=%.2f:%.2f:%.2f:%.2f:%.2f:%.2f",final_mass4l,Z2.pt1,Z2.eta1,Z2.pt2,Z2.eta2,ISO_CUT);
	   
	   ////sprintf (Eventformat,"MMMM_I_AI=%d:%d:%d:%.2f:%.2f:%.2f",Run,LumiSection,Event,massofhiggs,massofZ1,massofZ2);
	   
	   if ( abs(final_mass4l-Zmass)<=10.){
	     
	     hPFMET_MMMM_I_AI_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_MMMM_I_AI_Zpeak_log->Fill(RECO_PFMET,newweight);
	     
	     hPFMET_MMMM_3p1f_total_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_MMMM_3p1f_total_Zpeak_log->Fill(RECO_PFMET,newweight);
	     
	     //cout<<"MMMM_I_AI_Zpeak=== "<<final_mass4l<<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   }
	   
	 }//end I_AI
	 
	 if((Iso_1 < 0.35 && PF_1==1) && (Iso_2 < 0.35 && PF_2==1)){//2p2p
	   //cout<<">>>> 4mu 2P2P <<<<"<<endl;
	   //---- Set tree entries ----
	   f_run2 = Run;
	   f_lumi2 = LumiSection;
	   f_event2 = Event;
	   f_weight2 = newweight;
	   f_2p2p = 1;
	   f_lept1_pt = vlepts[sort_indexes[0]][0];
	   f_lept1_eta = vlepts[sort_indexes[0]][1];
	   f_lept1_phi = vlepts[sort_indexes[0]][2];
	   f_lept1_pdgid = vlepts[sort_indexes[0]][3];
	   f_lept2_pt = vlepts[sort_indexes[1]][0];
	   f_lept2_eta = vlepts[sort_indexes[1]][1];
	   f_lept2_phi = vlepts[sort_indexes[1]][2];
	   f_lept2_pdgid = vlepts[sort_indexes[1]][3];
	   f_lept3_pt = vlepts[sort_indexes[2]][0];
	   f_lept3_eta = vlepts[sort_indexes[2]][1];
	   f_lept3_phi = vlepts[sort_indexes[2]][2];
	   f_lept3_pdgid = vlepts[sort_indexes[2]][3];
	   f_lept4_pt = vlepts[sort_indexes[3]][0];
	   f_lept4_eta = vlepts[sort_indexes[3]][1];
	   f_lept4_phi = vlepts[sort_indexes[3]][2];
	   f_lept4_pdgid = vlepts[sort_indexes[3]][3];
	   f_lept1_pass = (vlepts[sort_indexes[0]][4]>=0.35 || vlepts[sort_indexes[0]][5]==0)? 0:1;
	   f_lept2_pass = (vlepts[sort_indexes[1]][4]>=0.35 || vlepts[sort_indexes[1]][5]==0)? 0:1;
	   f_lept3_pass = (vlepts[sort_indexes[2]][4]>=0.35 || vlepts[sort_indexes[2]][5]==0)? 0:1;
	   f_lept4_pass = (vlepts[sort_indexes[3]][4]>=0.35 || vlepts[sort_indexes[3]][5]==0)? 0:1;
	   f_lept1_pt_error = vlepts[sort_indexes[0]][6];
	   f_lept2_pt_error = vlepts[sort_indexes[1]][6];
	   f_lept3_pt_error = vlepts[sort_indexes[2]][6];
	   f_lept4_pt_error = vlepts[sort_indexes[3]][6];

	   
	   R_12_MMMM++;
	   R_12_MMMM_w=R_12_MMMM_w+newweight;
	   
	   //cout<<"MMMM_I_I=== "<< final_mass4l <<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<<final_pT4l<<" "<< weight <<" 0.35"<<endl;
	   ////sprintf (outformat_MMMM,"MMMM_I_I=%.2f:%.2f:%.2f:%.2f:%.2f:%.2f",final_mass4l,Z2.pt1,Z2.eta1,Z2.pt2,Z2.eta2,ISO_CUT);
	   
	   //sprintf (Eventformat,"MMMM_I_I=%d:%d:%d:%.2f:%.2f:%.2f",Run,LumiSection,Event,massofhiggs,massofZ1,massofZ2);
	   
	   hM4l_MMMM_I_I->Fill(massofhiggs,newweight);
	   hMZ1_MMMM_I_I->Fill(massofZ1,newweight);
	   hMZ2_MMMM_I_I->Fill(massofZ2,newweight);
	   
	 }//end I_I
	 
       }// end of Z1(MuMu) & Z2(mumu)
     
     else if ( Z1.tag==1 && Z2.tag==2 ) //Z1->mumu & Z2->ee
       
       {
	 //cout<<"Z1 --> u u && Z2 --> e e"<<endl;
	 
	 bool BDT_ok_1 = 0; // Spring16 with CMSSW_8_2_0
	 
	 if( RECOELE_PT[Z2.ilept1] > 7. &&  RECOELE_PT[Z2.ilept1] <= 10. ){
	   if( fabs(RECOELE_scl_Eta[Z2.ilept1]) < .8 && RECOELE_mvaNonTrigV0[Z2.ilept1] > -0.211 ) BDT_ok_1 = 1 ;
	   if( ( fabs(RECOELE_scl_Eta[Z2.ilept1]) >= .8 && fabs(RECOELE_scl_Eta[Z2.ilept1]) < 1.479 )
	       && RECOELE_mvaNonTrigV0[Z2.ilept1] > -0.396 ) BDT_ok_1 = 1 ;
	   if( fabs(RECOELE_scl_Eta[Z2.ilept1]) >= 1.479 && RECOELE_mvaNonTrigV0[Z2.ilept1] > -0.215 ) BDT_ok_1 = 1 ;
	 }
	 else { 
	   if( fabs(RECOELE_scl_Eta[Z2.ilept1]) < .8 && RECOELE_mvaNonTrigV0[Z2.ilept1] > -0.870 ) BDT_ok_1 = 1 ;
	   if( ( fabs(RECOELE_scl_Eta[Z2.ilept1]) >= .8 && fabs(RECOELE_scl_Eta[Z2.ilept1]) <= 1.479 )
	       && RECOELE_mvaNonTrigV0[Z2.ilept1] > -0.838 ) BDT_ok_1 = 1 ;
	   if( fabs(RECOELE_scl_Eta[Z2.ilept1]) > 1.479 && RECOELE_mvaNonTrigV0[Z2.ilept1] > -0.763 ) BDT_ok_1 = 1 ;
	 }

	 
	 /////////////////////////
	 
	 bool BDT_ok_2 = 0; // Spring16 with CMSSW_8_2_0
	 
	 if( RECOELE_PT[Z2.ilept2] > 7. &&  RECOELE_PT[Z2.ilept2] <= 10. ){
	   if( fabs(RECOELE_scl_Eta[Z2.ilept2]) < .8 && RECOELE_mvaNonTrigV0[Z2.ilept2] > -0.211 ) BDT_ok_2 = 1 ;
	   if( ( fabs(RECOELE_scl_Eta[Z2.ilept2]) >= .8 && fabs(RECOELE_scl_Eta[Z2.ilept2]) < 1.479 )
	       && RECOELE_mvaNonTrigV0[Z2.ilept2] > -0.396 ) BDT_ok_2 = 1 ;
	   if( fabs(RECOELE_scl_Eta[Z2.ilept2]) >= 1.479 && RECOELE_mvaNonTrigV0[Z2.ilept2] > -0.215 ) BDT_ok_2 = 1 ;
	 }
	 else {
	   if( fabs(RECOELE_scl_Eta[Z2.ilept2]) < .8 && RECOELE_mvaNonTrigV0[Z2.ilept2] > -0.870 ) BDT_ok_2 = 1 ;
	   if( ( fabs(RECOELE_scl_Eta[Z2.ilept2]) >= .8 && fabs(RECOELE_scl_Eta[Z2.ilept2]) <= 1.479 )
	       && RECOELE_mvaNonTrigV0[Z2.ilept2] > -0.838 ) BDT_ok_2 = 1 ;
	   if( fabs(RECOELE_scl_Eta[Z2.ilept2]) > 1.479 && RECOELE_mvaNonTrigV0[Z2.ilept2] > -0.763 ) BDT_ok_2 = 1 ;
	 }
	 
	 
	 /////////////////////////////
	 
	 PF_1=BDT_ok_1;
	 PF_2=BDT_ok_2;
	 
	 //////cout<<"Event MMEE has Z2 with 2 electrons with properties "<< "PF_1 = "<<PF_1 <<" and Iso_1 = "<<Iso_1<<"and PF_2 = "<<PF_2<<"and Iso_2 = "<<Iso_2 <<endl;
	 
	 if((Iso_1 >= 0.35 || PF_1==0) && (Iso_2 >= 0.35 || PF_2==0)){//2p2f
	   //cout<<">>>> 2mu2e 2P2F <<<<"<<endl;
	   //---- Set tree entries ----
	   f_run2 = Run;
	   f_lumi2 = LumiSection;
	   f_event2 = Event;
	   f_weight2 = newweight;
	   f_2p2f = 1;
	   f_lept1_pt = vlepts[sort_indexes[0]][0];
	   f_lept1_eta = vlepts[sort_indexes[0]][1];
	   f_lept1_phi = vlepts[sort_indexes[0]][2];
	   f_lept1_pdgid = vlepts[sort_indexes[0]][3];
	   f_lept2_pt = vlepts[sort_indexes[1]][0];
	   f_lept2_eta = vlepts[sort_indexes[1]][1];
	   f_lept2_phi = vlepts[sort_indexes[1]][2];
	   f_lept2_pdgid = vlepts[sort_indexes[1]][3];
	   f_lept3_pt = vlepts[sort_indexes[2]][0];
	   f_lept3_eta = vlepts[sort_indexes[2]][1];
	   f_lept3_phi = vlepts[sort_indexes[2]][2];
	   f_lept3_pdgid = vlepts[sort_indexes[2]][3];
	   f_lept4_pt = vlepts[sort_indexes[3]][0];
	   f_lept4_eta = vlepts[sort_indexes[3]][1];
	   f_lept4_phi = vlepts[sort_indexes[3]][2];
	   f_lept4_pdgid = vlepts[sort_indexes[3]][3];
	   f_lept1_pass = (vlepts[sort_indexes[0]][4]>=0.35 || vlepts[sort_indexes[0]][5]==0)? 0:1;
	   f_lept2_pass = (vlepts[sort_indexes[1]][4]>=0.35 || vlepts[sort_indexes[1]][5]==0)? 0:1;
	   f_lept3_pass = (vlepts[sort_indexes[2]][4]>=0.35 || vlepts[sort_indexes[2]][5]==0)? 0:1;
	   f_lept4_pass = (vlepts[sort_indexes[3]][4]>=0.35 || vlepts[sort_indexes[3]][5]==0)? 0:1;
	   f_lept1_pt_error = vlepts[sort_indexes[0]][6];
	   f_lept2_pt_error = vlepts[sort_indexes[1]][6];
	   f_lept3_pt_error = vlepts[sort_indexes[2]][6];
	   f_lept4_pt_error = vlepts[sort_indexes[3]][6];
	   
	   
	   R_AI_AI_MMEE++;
	   R_AI_AI_MMEE_w=R_AI_AI_MMEE_w+newweight;
	   
	   hM4l_MMEE_AI_AI->Fill(final_mass4l,newweight);
	   hPFMET_MMEE_AI_AI->Fill(RECO_PFMET,newweight);
	   
	   //cout<<"MMEE_AI_AI=== "<< final_mass4l <<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   ////sprintf (outformat_MMEE,"MMEE_AI_AI=%.2f:%.2f:%.2f:%.2f:%.2f:%.2f",final_mass4l,Z2.pt1,Z2.eta1,Z2.pt2,Z2.eta2,ISO_CUT);
	   
	   //sprintf (Eventformat,"MMEE_AI_AI=%d:%d:%d:%.2f:%.2f:%.2f",Run,LumiSection,Event,massofhiggs,massofZ1,massofZ2);
	   
	   if ( abs(final_mass4l-Zmass)<=10.){
	     
	     hPFMET_MMEE_AI_AI_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_MMEE_AI_AI_Zpeak_log->Fill(RECO_PFMET,newweight);
	     
	     //cout<<"MMEE_AI_AI_Zpeak=== "<<final_mass4l<<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   }
	   
	 }//AI_AI
	 
	 if((Iso_1 >= 0.35 || PF_1==0) && (Iso_2 < 0.35 && PF_2==1)){//3p1f
	   //cout<<">>>> 2mu2e 3P1F - 1st lepton fails <<<<"<<endl;
	   //---- Set tree entries ----
	   f_run2 = Run;
	   f_lumi2 = LumiSection;
	   f_event2 = Event;
	   f_weight2 = newweight;
	   f_3p1f = 1;
	   f_lept1_pt = vlepts[sort_indexes[0]][0];
	   f_lept1_eta = vlepts[sort_indexes[0]][1];
	   f_lept1_phi = vlepts[sort_indexes[0]][2];
	   f_lept1_pdgid = vlepts[sort_indexes[0]][3];
	   f_lept2_pt = vlepts[sort_indexes[1]][0];
	   f_lept2_eta = vlepts[sort_indexes[1]][1];
	   f_lept2_phi = vlepts[sort_indexes[1]][2];
	   f_lept2_pdgid = vlepts[sort_indexes[1]][3];
	   f_lept3_pt = vlepts[sort_indexes[2]][0];
	   f_lept3_eta = vlepts[sort_indexes[2]][1];
	   f_lept3_phi = vlepts[sort_indexes[2]][2];
	   f_lept3_pdgid = vlepts[sort_indexes[2]][3];
	   f_lept4_pt = vlepts[sort_indexes[3]][0];
	   f_lept4_eta = vlepts[sort_indexes[3]][1];
	   f_lept4_phi = vlepts[sort_indexes[3]][2];
	   f_lept4_pdgid = vlepts[sort_indexes[3]][3];
	   f_lept1_pass = (vlepts[sort_indexes[0]][4]>=0.35 || vlepts[sort_indexes[0]][5]==0)? 0:1;
	   f_lept2_pass = (vlepts[sort_indexes[1]][4]>=0.35 || vlepts[sort_indexes[1]][5]==0)? 0:1;
	   f_lept3_pass = (vlepts[sort_indexes[2]][4]>=0.35 || vlepts[sort_indexes[2]][5]==0)? 0:1;
	   f_lept4_pass = (vlepts[sort_indexes[3]][4]>=0.35 || vlepts[sort_indexes[3]][5]==0)? 0:1;
	   f_lept1_pt_error = vlepts[sort_indexes[0]][6];
	   f_lept2_pt_error = vlepts[sort_indexes[1]][6];
	   f_lept3_pt_error = vlepts[sort_indexes[2]][6];
	   f_lept4_pt_error = vlepts[sort_indexes[3]][6];
	   
	   
	   R_AI_I_MMEE++;
	   R_AI_I_MMEE_w=R_AI_I_MMEE_w+newweight;
	   
	   hM4l_MMEE_AI_I->Fill(final_mass4l,newweight);
	   hM4l_MMEE_3p1f_total->Fill(final_mass4l,newweight);

	   hPFMET_MMEE_AI_I->Fill(RECO_PFMET,newweight);
	   hPFMET_MMEE_3p1f_total->Fill(RECO_PFMET,newweight);
	   
	   //cout<<"MMEE_AI_I=== "<< final_mass4l <<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   ////sprintf (outformat_MMEE,"MMEE_AI_I=%.2f:%.2f:%.2f:%.2f:%.2f:%.2f",final_mass4l,Z2.pt1,Z2.eta1,Z2.pt2,Z2.eta2,ISO_CUT);
	   
	   //sprintf (Eventformat,"MMEE_AI_I=%d:%d:%d:%.2f:%.2f:%.2f",Run,LumiSection,Event,massofhiggs,massofZ1,massofZ2);

	   if ( abs(final_mass4l-Zmass)<=10.){
	     
	     hPFMET_MMEE_AI_I_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_MMEE_AI_I_Zpeak_log->Fill(RECO_PFMET,newweight);
	     
	     hPFMET_MMEE_3p1f_total_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_MMEE_3p1f_total_Zpeak_log->Fill(RECO_PFMET,newweight);
	     
	     //cout<<"MMEE_AI_I_Zpeak=== "<<final_mass4l<<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   }
	   
	 }//AI_I
	 
	 if((Iso_1 < 0.35 && PF_1==1) && (Iso_2 >= 0.35 || PF_2==0)){//3p1f
	   //cout<<">>>> 2mu2e 3P1F - 2nd lepton fails <<<<"<<endl;
	   //---- Set tree entries ----
	   f_run2 = Run;
	   f_lumi2 = LumiSection;
	   f_event2 = Event;
	   f_weight2 = newweight;
	   f_3p1f = 1;
	   f_lept1_pt = vlepts[sort_indexes[0]][0];
	   f_lept1_eta = vlepts[sort_indexes[0]][1];
	   f_lept1_phi = vlepts[sort_indexes[0]][2];
	   f_lept1_pdgid = vlepts[sort_indexes[0]][3];
	   f_lept2_pt = vlepts[sort_indexes[1]][0];
	   f_lept2_eta = vlepts[sort_indexes[1]][1];
	   f_lept2_phi = vlepts[sort_indexes[1]][2];
	   f_lept2_pdgid = vlepts[sort_indexes[1]][3];
	   f_lept3_pt = vlepts[sort_indexes[2]][0];
	   f_lept3_eta = vlepts[sort_indexes[2]][1];
	   f_lept3_phi = vlepts[sort_indexes[2]][2];
	   f_lept3_pdgid = vlepts[sort_indexes[2]][3];
	   f_lept4_pt = vlepts[sort_indexes[3]][0];
	   f_lept4_eta = vlepts[sort_indexes[3]][1];
	   f_lept4_phi = vlepts[sort_indexes[3]][2];
	   f_lept4_pdgid = vlepts[sort_indexes[3]][3];
	   f_lept1_pass = (vlepts[sort_indexes[0]][4]>=0.35 || vlepts[sort_indexes[0]][5]==0)? 0:1;
	   f_lept2_pass = (vlepts[sort_indexes[1]][4]>=0.35 || vlepts[sort_indexes[1]][5]==0)? 0:1;
	   f_lept3_pass = (vlepts[sort_indexes[2]][4]>=0.35 || vlepts[sort_indexes[2]][5]==0)? 0:1;
	   f_lept4_pass = (vlepts[sort_indexes[3]][4]>=0.35 || vlepts[sort_indexes[3]][5]==0)? 0:1;
	   f_lept1_pt_error = vlepts[sort_indexes[0]][6];
	   f_lept2_pt_error = vlepts[sort_indexes[1]][6];
	   f_lept3_pt_error = vlepts[sort_indexes[2]][6];
	   f_lept4_pt_error = vlepts[sort_indexes[3]][6];
	   
	   
	   R_I_AI_MMEE++;
	   R_I_AI_MMEE_w=R_I_AI_MMEE_w+newweight;
	   
	   hM4l_MMEE_I_AI->Fill(final_mass4l,newweight);  
	   hM4l_MMEE_3p1f_total->Fill(final_mass4l,newweight);

	   hPFMET_MMEE_I_AI->Fill(RECO_PFMET,newweight);
	   hPFMET_MMEE_3p1f_total->Fill(RECO_PFMET,newweight);
	   
	   //cout<<"MMEE_I_AI=== "<< final_mass4l <<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   ////sprintf (outformat_MMEE,"MMEE_I_AI=%.2f:%.2f:%.2f:%.2f:%.2f:%.2f",final_mass4l,Z2.pt1,Z2.eta1,Z2.pt2,Z2.eta2,ISO_CUT);
	   
	   //sprintf (Eventformat,"MMEE_I_AI=%d:%d:%d:%.2f:%.2f:%.2f",Run,LumiSection,Event,massofhiggs,massofZ1,massofZ2);

	   if ( abs(final_mass4l-Zmass)<=10.){
	     
	     hPFMET_MMEE_I_AI_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_MMEE_I_AI_Zpeak_log->Fill(RECO_PFMET,newweight);
	     
	     hPFMET_MMEE_3p1f_total_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_MMEE_3p1f_total_Zpeak_log->Fill(RECO_PFMET,newweight);
	     
	     //cout<<"MMEE_I_AI_Zpeak=== "<<final_mass4l<<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   }
	   
	 }//I_AI
	 
	 if((Iso_1 < 0.35 && PF_1==1) && (Iso_2 < 0.35 && PF_2==1)){//2p2p
	   //cout<<">>>> 2mu2e 2P2P <<<<"<<endl;
	   //---- Set tree entries ----
	   f_run2 = Run;
	   f_lumi2 = LumiSection;
	   f_event2 = Event;
	   f_weight2 = newweight;
	   f_2p2p = 1;
	   f_lept1_pt = vlepts[sort_indexes[0]][0];
	   f_lept1_eta = vlepts[sort_indexes[0]][1];
	   f_lept1_phi = vlepts[sort_indexes[0]][2];
	   f_lept1_pdgid = vlepts[sort_indexes[0]][3];
	   f_lept2_pt = vlepts[sort_indexes[1]][0];
	   f_lept2_eta = vlepts[sort_indexes[1]][1];
	   f_lept2_phi = vlepts[sort_indexes[1]][2];
	   f_lept2_pdgid = vlepts[sort_indexes[1]][3];
	   f_lept3_pt = vlepts[sort_indexes[2]][0];
	   f_lept3_eta = vlepts[sort_indexes[2]][1];
	   f_lept3_phi = vlepts[sort_indexes[2]][2];
	   f_lept3_pdgid = vlepts[sort_indexes[2]][3];
	   f_lept4_pt = vlepts[sort_indexes[3]][0];
	   f_lept4_eta = vlepts[sort_indexes[3]][1];
	   f_lept4_phi = vlepts[sort_indexes[3]][2];
	   f_lept4_pdgid = vlepts[sort_indexes[3]][3];
	   f_lept1_pass = (vlepts[sort_indexes[0]][4]>=0.35 || vlepts[sort_indexes[0]][5]==0)? 0:1;
	   f_lept2_pass = (vlepts[sort_indexes[1]][4]>=0.35 || vlepts[sort_indexes[1]][5]==0)? 0:1;
	   f_lept3_pass = (vlepts[sort_indexes[2]][4]>=0.35 || vlepts[sort_indexes[2]][5]==0)? 0:1;
	   f_lept4_pass = (vlepts[sort_indexes[3]][4]>=0.35 || vlepts[sort_indexes[3]][5]==0)? 0:1;
	   f_lept1_pt_error = vlepts[sort_indexes[0]][6];
	   f_lept2_pt_error = vlepts[sort_indexes[1]][6];
	   f_lept3_pt_error = vlepts[sort_indexes[2]][6];
	   f_lept4_pt_error = vlepts[sort_indexes[3]][6];
	   
	   
	   R_12_MMEE++;
	   R_12_MMEE_w=R_12_MMEE_w+newweight;
	   
	   //cout<<"MMEE_I_I=== "<< final_mass4l <<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   ////sprintf (outformat_MMEE,"MMEE_I_I=%.2f:%.2f:%.2f:%.2f:%.2f:%.2f",final_mass4l,Z2.pt1,Z2.eta1,Z2.pt2,Z2.eta2,ISO_CUT);
	   
	   //sprintf (Eventformat,"MMEE_I_I=%d:%d:%d:%.2f:%.2f:%.2f",Run,LumiSection,Event,massofhiggs,massofZ1,massofZ2);
	   
	   hM4l_MMEE_I_I->Fill(massofhiggs,newweight);
	   hMZ1_MMEE_I_I->Fill(massofZ1,newweight);
	   hMZ2_MMEE_I_I->Fill(massofZ2,newweight);
	   
	 }
	 
	 
       } //end of ( Z1tag==1 && Z2tag==2 )
     
     else if ( Z1.tag==2 && Z2.tag==2 ) //Z1->ee & Z2->ee
       
       {
	//cout<<"Z1 --> e e && Z2 --> e e"<<endl; 
	 bool BDT_ok_1 = 0; // Spring16 with CMSSW_8_2_0
	 
	 if( RECOELE_PT[Z2.ilept1] > 7. &&  RECOELE_PT[Z2.ilept1] <= 10. ){
	   if( fabs(RECOELE_scl_Eta[Z2.ilept1]) < .8 && RECOELE_mvaNonTrigV0[Z2.ilept1] > -0.211 ) BDT_ok_1 = 1 ;
	   if( ( fabs(RECOELE_scl_Eta[Z2.ilept1]) >= .8 && fabs(RECOELE_scl_Eta[Z2.ilept1]) < 1.479 )
	       && RECOELE_mvaNonTrigV0[Z2.ilept1] > -0.396 ) BDT_ok_1 = 1 ;
	   if( fabs(RECOELE_scl_Eta[Z2.ilept1]) >= 1.479 && RECOELE_mvaNonTrigV0[Z2.ilept1] > -0.215 ) BDT_ok_1 = 1 ;
	 }
	 else { 
	   if( fabs(RECOELE_scl_Eta[Z2.ilept1]) < .8 && RECOELE_mvaNonTrigV0[Z2.ilept1] > -0.870 ) BDT_ok_1 = 1 ;
	   if( ( fabs(RECOELE_scl_Eta[Z2.ilept1]) >= .8 && fabs(RECOELE_scl_Eta[Z2.ilept1]) <= 1.479 )
	       && RECOELE_mvaNonTrigV0[Z2.ilept1] > -0.838 ) BDT_ok_1 = 1 ;
	   if( fabs(RECOELE_scl_Eta[Z2.ilept1]) > 1.479 && RECOELE_mvaNonTrigV0[Z2.ilept1] > -0.763 ) BDT_ok_1 = 1 ;
	 }
	 
	 /////////////////////////
	 
	 bool BDT_ok_2 = 0; // Spring16 with CMSSW_8_2_0
	 
	 if( RECOELE_PT[Z2.ilept2] > 7. &&  RECOELE_PT[Z2.ilept2] <= 10. ){
	   if( fabs(RECOELE_scl_Eta[Z2.ilept2]) < .8 && RECOELE_mvaNonTrigV0[Z2.ilept2] > -0.211 ) BDT_ok_2 = 1 ;
	   if( ( fabs(RECOELE_scl_Eta[Z2.ilept2]) >= .8 && fabs(RECOELE_scl_Eta[Z2.ilept2]) < 1.479 )
	       && RECOELE_mvaNonTrigV0[Z2.ilept2] > -0.396 ) BDT_ok_2 = 1 ;
	   if( fabs(RECOELE_scl_Eta[Z2.ilept2]) >= 1.479 && RECOELE_mvaNonTrigV0[Z2.ilept2] > -0.215 ) BDT_ok_2 = 1 ;
	 }
	 else {
	   if( fabs(RECOELE_scl_Eta[Z2.ilept2]) < .8 && RECOELE_mvaNonTrigV0[Z2.ilept2] > -0.870 ) BDT_ok_2 = 1 ;
	   if( ( fabs(RECOELE_scl_Eta[Z2.ilept2]) >= .8 && fabs(RECOELE_scl_Eta[Z2.ilept2]) <= 1.479 )
	       && RECOELE_mvaNonTrigV0[Z2.ilept2] > -0.838 ) BDT_ok_2 = 1 ;
	   if( fabs(RECOELE_scl_Eta[Z2.ilept2]) > 1.479 && RECOELE_mvaNonTrigV0[Z2.ilept2] > -0.763 ) BDT_ok_2 = 1 ;
	 }
	 
	 /////////////////////////////
	 
	 PF_1=BDT_ok_1;
	 PF_2=BDT_ok_2;
	 
	 ////cout<<"Event EEEE has Z2 with 2 electrons with properties "<< "PF_1 = "<<PF_1 <<" and Iso_1 = "<<Iso_1<<"and PF_2 = "<<PF_2<<"and Iso_2 = "<<Iso_2 <<endl;
	 
	 if((Iso_1 >= 0.35 || PF_1==0) && (Iso_2 >= 0.35 || PF_2==0)){//2p2f
	   //cout<<">>>> 4e 2P2F <<<<"<<endl;
	   //---- Set tree entries ----
	   f_run2 = Run;
	   f_lumi2 = LumiSection;
	   f_event2 = Event;
	   f_weight2 = newweight;
	   f_2p2f = 1;	   
	   f_lept1_pt = vlepts[sort_indexes[0]][0];
	   f_lept1_eta = vlepts[sort_indexes[0]][1];
	   f_lept1_phi = vlepts[sort_indexes[0]][2];
	   f_lept1_pdgid = vlepts[sort_indexes[0]][3];
	   f_lept2_pt = vlepts[sort_indexes[1]][0];
	   f_lept2_eta = vlepts[sort_indexes[1]][1];
	   f_lept2_phi = vlepts[sort_indexes[1]][2];
	   f_lept2_pdgid = vlepts[sort_indexes[1]][3];
	   f_lept3_pt = vlepts[sort_indexes[2]][0];
	   f_lept3_eta = vlepts[sort_indexes[2]][1];
	   f_lept3_phi = vlepts[sort_indexes[2]][2];
	   f_lept3_pdgid = vlepts[sort_indexes[2]][3];
	   f_lept4_pt = vlepts[sort_indexes[3]][0];
	   f_lept4_eta = vlepts[sort_indexes[3]][1];
	   f_lept4_phi = vlepts[sort_indexes[3]][2];
	   f_lept4_pdgid = vlepts[sort_indexes[3]][3];
	   f_lept1_pass = (vlepts[sort_indexes[0]][4]>=0.35 || vlepts[sort_indexes[0]][5]==0)? 0:1;
	   f_lept2_pass = (vlepts[sort_indexes[1]][4]>=0.35 || vlepts[sort_indexes[1]][5]==0)? 0:1;
	   f_lept3_pass = (vlepts[sort_indexes[2]][4]>=0.35 || vlepts[sort_indexes[2]][5]==0)? 0:1;
	   f_lept4_pass = (vlepts[sort_indexes[3]][4]>=0.35 || vlepts[sort_indexes[3]][5]==0)? 0:1;
	   f_lept1_pt_error = vlepts[sort_indexes[0]][6];
	   f_lept2_pt_error = vlepts[sort_indexes[1]][6];
	   f_lept3_pt_error = vlepts[sort_indexes[2]][6];
	   f_lept4_pt_error = vlepts[sort_indexes[3]][6];
	   
	   
	   R_AI_AI_EEEE++;
	   R_AI_AI_EEEE_w=R_AI_AI_EEEE_w+newweight;
	   
	   hM4l_EEEE_AI_AI->Fill(final_mass4l,newweight);
	   hPFMET_EEEE_AI_AI->Fill(RECO_PFMET,newweight);
	   
	   //cout<<"EEEE_AI_AI=== "<< final_mass4l <<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   ////sprintf (outformat_EEEE,"EEEE_AI_AI=%.2f:%.2f:%.2f:%.2f:%.2f:%.2f",final_mass4l,Z2.pt1,Z2.eta1,Z2.pt2,Z2.eta2,ISO_CUT);
	   
	   //sprintf (Eventformat,"EEEE_AI_AI=%d:%d:%d:%.2f:%.2f:%.2f",Run,LumiSection,Event,massofhiggs,massofZ1,massofZ2);

	   	   
	   if ( abs(final_mass4l-Zmass)<=10.){
	     
	     hPFMET_EEEE_AI_AI_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_EEEE_AI_AI_Zpeak_log->Fill(RECO_PFMET,newweight);
	     
	     //cout<<"EEEE_AI_AI_Zpeak=== "<<final_mass4l<<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   }

	   
	 }//AI_AI
	 
	 if((Iso_1 >= 0.35 || PF_1==0) && (Iso_2 < 0.35 && PF_2==1)){//3p1f
	   //cout<<">>>> 4e 3P1F - 1st lepton fails <<<<"<<endl;
	   //---- Set tree entries ----
	   f_run2 = Run;
	   f_lumi2 = LumiSection;
	   f_event2 = Event;
	   f_weight2 = newweight;
	   f_3p1f = 1;	   
	   f_lept1_pt = vlepts[sort_indexes[0]][0];
	   f_lept1_eta = vlepts[sort_indexes[0]][1];
	   f_lept1_phi = vlepts[sort_indexes[0]][2];
	   f_lept1_pdgid = vlepts[sort_indexes[0]][3];
	   f_lept2_pt = vlepts[sort_indexes[1]][0];
	   f_lept2_eta = vlepts[sort_indexes[1]][1];
	   f_lept2_phi = vlepts[sort_indexes[1]][2];
	   f_lept2_pdgid = vlepts[sort_indexes[1]][3];
	   f_lept3_pt = vlepts[sort_indexes[2]][0];
	   f_lept3_eta = vlepts[sort_indexes[2]][1];
	   f_lept3_phi = vlepts[sort_indexes[2]][2];
	   f_lept3_pdgid = vlepts[sort_indexes[2]][3];
	   f_lept4_pt = vlepts[sort_indexes[3]][0];
	   f_lept4_eta = vlepts[sort_indexes[3]][1];
	   f_lept4_phi = vlepts[sort_indexes[3]][2];
	   f_lept4_pdgid = vlepts[sort_indexes[3]][3];
	   f_lept1_pass = (vlepts[sort_indexes[0]][4]>=0.35 || vlepts[sort_indexes[0]][5]==0)? 0:1;
	   f_lept2_pass = (vlepts[sort_indexes[1]][4]>=0.35 || vlepts[sort_indexes[1]][5]==0)? 0:1;
	   f_lept3_pass = (vlepts[sort_indexes[2]][4]>=0.35 || vlepts[sort_indexes[2]][5]==0)? 0:1;
	   f_lept4_pass = (vlepts[sort_indexes[3]][4]>=0.35 || vlepts[sort_indexes[3]][5]==0)? 0:1;
	   f_lept1_pt_error = vlepts[sort_indexes[0]][6];
	   f_lept2_pt_error = vlepts[sort_indexes[1]][6];
	   f_lept3_pt_error = vlepts[sort_indexes[2]][6];
	   f_lept4_pt_error = vlepts[sort_indexes[3]][6];
	   
	   
	   R_AI_I_EEEE++;
	   R_AI_I_EEEE_w=R_AI_I_EEEE_w+newweight;
	   
	   hM4l_EEEE_AI_I->Fill(final_mass4l,newweight);   
	   hM4l_EEEE_3p1f_total->Fill(final_mass4l,newweight);

	   hPFMET_EEEE_AI_I->Fill(RECO_PFMET,newweight);
	   hPFMET_EEEE_3p1f_total->Fill(RECO_PFMET,newweight);
	   
	   //cout<<"EEEE_AI_I=== "<< final_mass4l <<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   ////sprintf (outformat_EEEE,"EEEE_AI_I=%.2f:%.2f:%.2f:%.2f:%.2f:%.2f",final_mass4l,Z2.pt1,Z2.eta1,Z2.pt2,Z2.eta2,ISO_CUT);
	   
	   //sprintf (Eventformat,"EEEE_AI_I=%d:%d:%d:%.2f:%.2f:%.2f",Run,LumiSection,Event,massofhiggs,massofZ1,massofZ2);

	   if ( abs(final_mass4l-Zmass)<=10.){
	     
	     hPFMET_EEEE_AI_I_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_EEEE_AI_I_Zpeak_log->Fill(RECO_PFMET,newweight);
	     
	     hPFMET_EEEE_3p1f_total_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_EEEE_3p1f_total_Zpeak_log->Fill(RECO_PFMET,newweight);
	     
	     //cout<<"EEEE_AI_I_Zpeak=== "<<final_mass4l<<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   }
	   
	 }//AI_I
	 
	 if((Iso_1 < 0.35 && PF_1==1) && (Iso_2 >= 0.35 || PF_2==0)){//3p1f
	   //cout<<">>>> 4e 3P1F - 2nd lepton fails <<<<"<<endl;
	   //---- Set tree entries ----
	   f_run2 = Run;
	   f_lumi2 = LumiSection;
	   f_event2 = Event;
	   f_weight2 = newweight;
	   f_3p1f = 1;
	   f_lept1_pt = vlepts[sort_indexes[0]][0];
	   f_lept1_eta = vlepts[sort_indexes[0]][1];
	   f_lept1_phi = vlepts[sort_indexes[0]][2];
	   f_lept1_pdgid = vlepts[sort_indexes[0]][3];
	   f_lept2_pt = vlepts[sort_indexes[1]][0];
	   f_lept2_eta = vlepts[sort_indexes[1]][1];
	   f_lept2_phi = vlepts[sort_indexes[1]][2];
	   f_lept2_pdgid = vlepts[sort_indexes[1]][3];
	   f_lept3_pt = vlepts[sort_indexes[2]][0];
	   f_lept3_eta = vlepts[sort_indexes[2]][1];
	   f_lept3_phi = vlepts[sort_indexes[2]][2];
	   f_lept3_pdgid = vlepts[sort_indexes[2]][3];
	   f_lept4_pt = vlepts[sort_indexes[3]][0];
	   f_lept4_eta = vlepts[sort_indexes[3]][1];
	   f_lept4_phi = vlepts[sort_indexes[3]][2];
	   f_lept4_pdgid = vlepts[sort_indexes[3]][3];
	   f_lept1_pass = (vlepts[sort_indexes[0]][4]>=0.35 || vlepts[sort_indexes[0]][5]==0)? 0:1;
	   f_lept2_pass = (vlepts[sort_indexes[1]][4]>=0.35 || vlepts[sort_indexes[1]][5]==0)? 0:1;
	   f_lept3_pass = (vlepts[sort_indexes[2]][4]>=0.35 || vlepts[sort_indexes[2]][5]==0)? 0:1;
	   f_lept4_pass = (vlepts[sort_indexes[3]][4]>=0.35 || vlepts[sort_indexes[3]][5]==0)? 0:1;
	   f_lept1_pt_error = vlepts[sort_indexes[0]][6];
	   f_lept2_pt_error = vlepts[sort_indexes[1]][6];
	   f_lept3_pt_error = vlepts[sort_indexes[2]][6];
	   f_lept4_pt_error = vlepts[sort_indexes[3]][6];
	   
	   
	   R_I_AI_EEEE++;
	   R_I_AI_EEEE_w=R_I_AI_EEEE_w+newweight;
	   
	   hM4l_EEEE_I_AI->Fill(final_mass4l,newweight);  
	   hM4l_EEEE_3p1f_total->Fill(final_mass4l,newweight);

	   hPFMET_EEEE_I_AI->Fill(RECO_PFMET,newweight);
	   hPFMET_EEEE_3p1f_total->Fill(RECO_PFMET,newweight);
	   
	   //cout<<"EEEE_I_AI=== "<< final_mass4l <<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   ////sprintf (outformat_EEEE,"EEEE_I_AI=%.2f:%.2f:%.2f:%.2f:%.2f:%.2f",final_mass4l,Z2.pt1,Z2.eta1,Z2.pt2,Z2.eta2,ISO_CUT);
	   
	   //sprintf (Eventformat,"EEEE_I_AI=%d:%d:%d:%.2f:%.2f:%.2f",Run,LumiSection,Event,massofhiggs,massofZ1,massofZ2);

	   if ( abs(final_mass4l-Zmass)<=10.){
	     
	     hPFMET_EEEE_I_AI_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_EEEE_I_AI_Zpeak_log->Fill(RECO_PFMET,newweight);
	     
	     hPFMET_EEEE_3p1f_total_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_EEEE_3p1f_total_Zpeak_log->Fill(RECO_PFMET,newweight);
	     
	     //cout<<"EEEE_I_AI_Zpeak=== "<<final_mass4l<<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   }
	   
	 }//I_AI
	 
	 if((Iso_1 < 0.35 && PF_1==1) && (Iso_2 < 0.35 && PF_2==1)){//2p2p
	   //cout<<">>>> 4e 2P2P <<<<"<<endl;
	   //---- Set tree entries ----
	   f_run2 = Run;
	   f_lumi2 = LumiSection;
	   f_event2 = Event;
	   f_weight2 = newweight;
	   f_2p2p = 1;
	   f_lept1_pt = vlepts[sort_indexes[0]][0];
	   f_lept1_eta = vlepts[sort_indexes[0]][1];
	   f_lept1_phi = vlepts[sort_indexes[0]][2];
	   f_lept1_pdgid = vlepts[sort_indexes[0]][3];
	   f_lept2_pt = vlepts[sort_indexes[1]][0];
	   f_lept2_eta = vlepts[sort_indexes[1]][1];
	   f_lept2_phi = vlepts[sort_indexes[1]][2];
	   f_lept2_pdgid = vlepts[sort_indexes[1]][3];
	   f_lept3_pt = vlepts[sort_indexes[2]][0];
	   f_lept3_eta = vlepts[sort_indexes[2]][1];
	   f_lept3_phi = vlepts[sort_indexes[2]][2];
	   f_lept3_pdgid = vlepts[sort_indexes[2]][3];
	   f_lept4_pt = vlepts[sort_indexes[3]][0];
	   f_lept4_eta = vlepts[sort_indexes[3]][1];
	   f_lept4_phi = vlepts[sort_indexes[3]][2];
	   f_lept4_pdgid = vlepts[sort_indexes[3]][3];
	   f_lept1_pass = (vlepts[sort_indexes[0]][4]>=0.35 || vlepts[sort_indexes[0]][5]==0)? 0:1;
	   f_lept2_pass = (vlepts[sort_indexes[1]][4]>=0.35 || vlepts[sort_indexes[1]][5]==0)? 0:1;
	   f_lept3_pass = (vlepts[sort_indexes[2]][4]>=0.35 || vlepts[sort_indexes[2]][5]==0)? 0:1;
	   f_lept4_pass = (vlepts[sort_indexes[3]][4]>=0.35 || vlepts[sort_indexes[3]][5]==0)? 0:1;
	   f_lept1_pt_error = vlepts[sort_indexes[0]][6];
	   f_lept2_pt_error = vlepts[sort_indexes[1]][6];
	   f_lept3_pt_error = vlepts[sort_indexes[2]][6];
	   f_lept4_pt_error = vlepts[sort_indexes[3]][6];
	   
	   
	   R_12_EEEE++;
	   R_12_EEEE_w=R_12_EEEE_w+newweight;
	   
	   //cout<<"EEEE_I_I=== "<< final_mass4l <<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   ////sprintf (outformat_EEEE,"EEEE_I_I=%.2f:%.2f:%.2f:%.2f:%.2f:%.2f",final_mass4l,Z2.pt1,Z2.eta1,Z2.pt2,Z2.eta2,ISO_CUT);
	   
	   //sprintf (Eventformat,"EEEE_I_I=%d:%d:%d:%.2f:%.2f:%.2f",Run,LumiSection,Event,massofhiggs,massofZ1,massofZ2);
	   
	   hM4l_EEEE_I_I->Fill(massofhiggs,newweight);
	   hMZ1_EEEE_I_I->Fill(massofZ1,newweight);
	   hMZ2_EEEE_I_I->Fill(massofZ2,newweight);
	   
	 }
	 
       }// end of ( Z1tag==2 && Z2tag==2 )
     
     else if ( Z1.tag==2 && Z2.tag==1 ) //Z1->ee & Z2->mumu
       
       {
	 //cout<<"Z1 --> e e && Z2 --> u u"<<endl;
	 PF_1=RECOMU_isPFMu[Z2.ilept1];
	 PF_2=RECOMU_isPFMu[Z2.ilept2];
	 
	 ////cout<<"Event EEMM has Z2 with 2 muons with properties "<< "PF_1 = "<<PF_1 <<" and Iso_1 = "<<Iso_1<<"and PF_2 = "<<PF_2<<"and Iso_2 = "<<Iso_2 <<endl;
	 
	 if((Iso_1 >= 0.35 || PF_1==0) && (Iso_2 >= 0.35 || PF_2==0)){ //2p2f
	   //cout<<">>>> 2e2mu 2P2F <<<<"<<endl;
	   //---- Set tree entries ----
	   f_run2 = Run;
	   f_lumi2 = LumiSection;
	   f_event2 = Event;
	   f_weight2 = newweight;
	   f_2p2f = 1;	   
	   f_lept1_pt = vlepts[sort_indexes[0]][0];
	   f_lept1_eta = vlepts[sort_indexes[0]][1];
	   f_lept1_phi = vlepts[sort_indexes[0]][2];
	   f_lept1_pdgid = vlepts[sort_indexes[0]][3];
	   f_lept2_pt = vlepts[sort_indexes[1]][0];
	   f_lept2_eta = vlepts[sort_indexes[1]][1];
	   f_lept2_phi = vlepts[sort_indexes[1]][2];
	   f_lept2_pdgid = vlepts[sort_indexes[1]][3];
	   f_lept3_pt = vlepts[sort_indexes[2]][0];
	   f_lept3_eta = vlepts[sort_indexes[2]][1];
	   f_lept3_phi = vlepts[sort_indexes[2]][2];
	   f_lept3_pdgid = vlepts[sort_indexes[2]][3];
	   f_lept4_pt = vlepts[sort_indexes[3]][0];
	   f_lept4_eta = vlepts[sort_indexes[3]][1];
	   f_lept4_phi = vlepts[sort_indexes[3]][2];
	   f_lept4_pdgid = vlepts[sort_indexes[3]][3];
	   f_lept1_pass = (vlepts[sort_indexes[0]][4]>=0.35 || vlepts[sort_indexes[0]][5]==0)? 0:1;
	   f_lept2_pass = (vlepts[sort_indexes[1]][4]>=0.35 || vlepts[sort_indexes[1]][5]==0)? 0:1;
	   f_lept3_pass = (vlepts[sort_indexes[2]][4]>=0.35 || vlepts[sort_indexes[2]][5]==0)? 0:1;
	   f_lept4_pass = (vlepts[sort_indexes[3]][4]>=0.35 || vlepts[sort_indexes[3]][5]==0)? 0:1;
	   f_lept1_pt_error = vlepts[sort_indexes[0]][6];
	   f_lept2_pt_error = vlepts[sort_indexes[1]][6];
	   f_lept3_pt_error = vlepts[sort_indexes[2]][6];
	   f_lept4_pt_error = vlepts[sort_indexes[3]][6];
	   
	   
	   R_AI_AI_EEMM++;
	   R_AI_AI_EEMM_w=R_AI_AI_EEMM_w+newweight;
	   
	   hM4l_EEMM_AI_AI->Fill(final_mass4l,newweight);
	   hPFMET_EEMM_AI_AI->Fill(RECO_PFMET,newweight);
	   
	   //cout<<"EEMM_AI_AI=== "<<final_mass4l<<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   ////sprintf (outformat_EEMM,"EEMM_AI_AI=%.2f:%.2f:%.2f:%.2f:%.2f:%.2f",final_mass4l,Z2.pt1,Z2.eta1,Z2.pt2,Z2.eta2,ISO_CUT);
	   
	   ////cout<<" PF_1 = "<< PF_1<<" ,  PF_2 = "<< PF_2<<endl;
	   ////cout<<" Iso_1 = "<< Iso_1<<" ,  Iso_2 = "<< Iso_2<<endl;
	   
	   //sprintf (Eventformat,"EEMM_AI_AI=%d:%d:%d:%.2f:%.2f:%.2f",Run,LumiSection,Event,massofhiggs,massofZ1,massofZ2);

	   if ( abs(final_mass4l-Zmass)<=10.){
	     
	     hPFMET_EEMM_AI_AI_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_EEMM_AI_AI_Zpeak_log->Fill(RECO_PFMET,newweight);
	     
	     //cout<<"EEMM_AI_AI_Zpeak=== "<<final_mass4l<<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   }
	   
	 }//end AI_AI
	 
	 if((Iso_1 >= 0.35 || PF_1==0) && (Iso_2 < 0.35 && PF_2==1)){ //3p1f
	   //cout<<">>>> 2e2mu 3P1F - 1st lepton fails <<<<"<<endl;
	   //---- Set tree entries ----
	   f_run2 = Run;
	   f_lumi2 = LumiSection;
	   f_event2 = Event;
	   f_weight2 = newweight;
	   f_3p1f = 1;	   
	   f_lept1_pt = vlepts[sort_indexes[0]][0];
	   f_lept1_eta = vlepts[sort_indexes[0]][1];
	   f_lept1_phi = vlepts[sort_indexes[0]][2];
	   f_lept1_pdgid = vlepts[sort_indexes[0]][3];
	   f_lept2_pt = vlepts[sort_indexes[1]][0];
	   f_lept2_eta = vlepts[sort_indexes[1]][1];
	   f_lept2_phi = vlepts[sort_indexes[1]][2];
	   f_lept2_pdgid = vlepts[sort_indexes[1]][3];
	   f_lept3_pt = vlepts[sort_indexes[2]][0];
	   f_lept3_eta = vlepts[sort_indexes[2]][1];
	   f_lept3_phi = vlepts[sort_indexes[2]][2];
	   f_lept3_pdgid = vlepts[sort_indexes[2]][3];
	   f_lept4_pt = vlepts[sort_indexes[3]][0];
	   f_lept4_eta = vlepts[sort_indexes[3]][1];
	   f_lept4_phi = vlepts[sort_indexes[3]][2];
	   f_lept4_pdgid = vlepts[sort_indexes[3]][3];
	   f_lept1_pass = (vlepts[sort_indexes[0]][4]>=0.35 || vlepts[sort_indexes[0]][5]==0)? 0:1;
	   f_lept2_pass = (vlepts[sort_indexes[1]][4]>=0.35 || vlepts[sort_indexes[1]][5]==0)? 0:1;
	   f_lept3_pass = (vlepts[sort_indexes[2]][4]>=0.35 || vlepts[sort_indexes[2]][5]==0)? 0:1;
	   f_lept4_pass = (vlepts[sort_indexes[3]][4]>=0.35 || vlepts[sort_indexes[3]][5]==0)? 0:1;
	   f_lept1_pt_error = vlepts[sort_indexes[0]][6];
	   f_lept2_pt_error = vlepts[sort_indexes[1]][6];
	   f_lept3_pt_error = vlepts[sort_indexes[2]][6];
	   f_lept4_pt_error = vlepts[sort_indexes[3]][6];
	   
	   
	   R_AI_I_EEMM++;
	   R_AI_I_EEMM_w=R_AI_I_EEMM_w+newweight;
	   
	   hM4l_EEMM_AI_I->Fill(final_mass4l,newweight);  
	   hM4l_EEMM_3p1f_total->Fill(final_mass4l,newweight);

	   hPFMET_EEMM_AI_I->Fill(RECO_PFMET,newweight);
	   hPFMET_EEMM_3p1f_total->Fill(RECO_PFMET,newweight);
	   
	   //cout<<"EEMM_AI_I=== "<< final_mass4l <<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   ////sprintf (outformat_EEMM,"EEMM_AI_I=%.2f:%.2f:%.2f:%.2f:%.2f:%.2f",final_mass4l,Z2.pt1,Z2.eta1,Z2.pt2,Z2.eta2,ISO_CUT);
	   
	   //sprintf (Eventformat,"EEMM_AI_I=%d:%d:%d:%.2f:%.2f:%.2f",Run,LumiSection,Event,massofhiggs,massofZ1,massofZ2);

	   if ( abs(final_mass4l-Zmass)<=10.){
	     
	     hPFMET_EEMM_AI_I_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_EEMM_AI_I_Zpeak_log->Fill(RECO_PFMET,newweight);
	     
	     hPFMET_EEMM_3p1f_total_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_EEMM_3p1f_total_Zpeak_log->Fill(RECO_PFMET,newweight);
	     
	     //cout<<"EEMM_AI_I_Zpeak=== "<<final_mass4l<<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   }
	   
	 }//end AI_I
	 
	 if((Iso_1 < 0.35 && PF_1==1) && (Iso_2 >= 0.35 || PF_2==0)){//3p1f
	   //cout<<">>>> 2e2mu 3P1F - 2nd lepton fails <<<<"<<endl;
	   //---- Set tree entries ----
	   f_run2 = Run;
	   f_lumi2 = LumiSection;
	   f_event2 = Event;
	   f_weight2 = newweight;
	   f_3p1f = 1;	   
	   f_lept1_pt = vlepts[sort_indexes[0]][0];
	   f_lept1_eta = vlepts[sort_indexes[0]][1];
	   f_lept1_phi = vlepts[sort_indexes[0]][2];
	   f_lept1_pdgid = vlepts[sort_indexes[0]][3];
	   f_lept2_pt = vlepts[sort_indexes[1]][0];
	   f_lept2_eta = vlepts[sort_indexes[1]][1];
	   f_lept2_phi = vlepts[sort_indexes[1]][2];
	   f_lept2_pdgid = vlepts[sort_indexes[1]][3];
	   f_lept3_pt = vlepts[sort_indexes[2]][0];
	   f_lept3_eta = vlepts[sort_indexes[2]][1];
	   f_lept3_phi = vlepts[sort_indexes[2]][2];
	   f_lept3_pdgid = vlepts[sort_indexes[2]][3];
	   f_lept4_pt = vlepts[sort_indexes[3]][0];
	   f_lept4_eta = vlepts[sort_indexes[3]][1];
	   f_lept4_phi = vlepts[sort_indexes[3]][2];
	   f_lept4_pdgid = vlepts[sort_indexes[3]][3];
	   f_lept1_pass = (vlepts[sort_indexes[0]][4]>=0.35 || vlepts[sort_indexes[0]][5]==0)? 0:1;
	   f_lept2_pass = (vlepts[sort_indexes[1]][4]>=0.35 || vlepts[sort_indexes[1]][5]==0)? 0:1;
	   f_lept3_pass = (vlepts[sort_indexes[2]][4]>=0.35 || vlepts[sort_indexes[2]][5]==0)? 0:1;
	   f_lept4_pass = (vlepts[sort_indexes[3]][4]>=0.35 || vlepts[sort_indexes[3]][5]==0)? 0:1;
	   f_lept1_pt_error = vlepts[sort_indexes[0]][6];
	   f_lept2_pt_error = vlepts[sort_indexes[1]][6];
	   f_lept3_pt_error = vlepts[sort_indexes[2]][6];
	   f_lept4_pt_error = vlepts[sort_indexes[3]][6];
	   
	   
	   R_I_AI_EEMM++;
	   R_I_AI_EEMM_w=R_I_AI_EEMM_w+newweight;
	   
	   hM4l_EEMM_I_AI->Fill(final_mass4l,newweight);   
	   hM4l_EEMM_3p1f_total->Fill(final_mass4l,newweight);

	   hPFMET_EEMM_I_AI->Fill(RECO_PFMET,newweight);
	   hPFMET_EEMM_3p1f_total->Fill(RECO_PFMET,newweight);
	   
           //cout<<"EEMM_I_AI=== "<< final_mass4l <<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   ////sprintf (outformat_EEMM,"EEMM_I_AI=%.2f:%.2f:%.2f:%.2f:%.2f:%.2f",final_mass4l,Z2.pt1,Z2.eta1,Z2.pt2,Z2.eta2,ISO_CUT);
	   
	   //sprintf (Eventformat,"EEMM_I_AI=%d:%d:%d:%.2f:%.2f:%.2f",Run,LumiSection,Event,massofhiggs,massofZ1,massofZ2);

	   if ( abs(final_mass4l-Zmass)<=10.){
	     
	     hPFMET_EEMM_I_AI_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_EEMM_I_AI_Zpeak_log->Fill(RECO_PFMET,newweight);
	     
	     hPFMET_EEMM_3p1f_total_Zpeak->Fill(RECO_PFMET,newweight);
	     hPFMET_EEMM_3p1f_total_Zpeak_log->Fill(RECO_PFMET,newweight);
	     
	     //cout<<"EEMM_I_AI_Zpeak=== "<<final_mass4l<<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   }
	   
	 }//end I_AI
	 
	 if((Iso_1 < 0.35 && PF_1==1) && (Iso_2 < 0.35 && PF_2==1)){//2p2p
	   //cout<<">>>> 2e2mu 2P2P <<<<"<<endl;
	   //---- Set tree entries ----
	   f_run2 = Run;
	   f_lumi2 = LumiSection;
	   f_event2 = Event;
	   f_weight2 = newweight;
	   f_2p2p = 1;	   
	   f_lept1_pt = vlepts[sort_indexes[0]][0];
	   f_lept1_eta = vlepts[sort_indexes[0]][1];
	   f_lept1_phi = vlepts[sort_indexes[0]][2];
	   f_lept1_pdgid = vlepts[sort_indexes[0]][3];
	   f_lept2_pt = vlepts[sort_indexes[1]][0];
	   f_lept2_eta = vlepts[sort_indexes[1]][1];
	   f_lept2_phi = vlepts[sort_indexes[1]][2];
	   f_lept2_pdgid = vlepts[sort_indexes[1]][3];
	   f_lept3_pt = vlepts[sort_indexes[2]][0];
	   f_lept3_eta = vlepts[sort_indexes[2]][1];
	   f_lept3_phi = vlepts[sort_indexes[2]][2];
	   f_lept3_pdgid = vlepts[sort_indexes[2]][3];
	   f_lept4_pt = vlepts[sort_indexes[3]][0];
	   f_lept4_eta = vlepts[sort_indexes[3]][1];
	   f_lept4_phi = vlepts[sort_indexes[3]][2];
	   f_lept4_pdgid = vlepts[sort_indexes[3]][3];
	   f_lept1_pass = (vlepts[sort_indexes[0]][4]>=0.35 || vlepts[sort_indexes[0]][5]==0)? 0:1;
	   f_lept2_pass = (vlepts[sort_indexes[1]][4]>=0.35 || vlepts[sort_indexes[1]][5]==0)? 0:1;
	   f_lept3_pass = (vlepts[sort_indexes[2]][4]>=0.35 || vlepts[sort_indexes[2]][5]==0)? 0:1;
	   f_lept4_pass = (vlepts[sort_indexes[3]][4]>=0.35 || vlepts[sort_indexes[3]][5]==0)? 0:1;
	   f_lept1_pt_error = vlepts[sort_indexes[0]][6];
	   f_lept2_pt_error = vlepts[sort_indexes[1]][6];
	   f_lept3_pt_error = vlepts[sort_indexes[2]][6];
	   f_lept4_pt_error = vlepts[sort_indexes[3]][6];
	   
	   
	   R_12_EEMM++;
	   R_12_EEMM_w=R_12_EEMM_w+newweight;
	   
	   //cout<<"EEMM_I_I=== "<< final_mass4l <<" "<< Z2.pt1 <<" "<< Z2.eta1 <<" "<< Z2.pt2 <<" "<< Z2.eta2 <<" "<< RECO_PFMET <<" "<< weight <<" 0.35"<<endl;
	   ////sprintf (outformat_EEMM,"EEMM_I_I=%.2f:%.2f:%.2f:%.2f:%.2f:%.2f",final_mass4l,Z2.pt1,Z2.eta1,Z2.pt2,Z2.eta2,ISO_CUT);
	   
	   ////sprintf (Eventformat,"EEMM_I_I=%d:%d:%d:%.2f:%.2f:%.2f",Run,LumiSection,Event,massofhiggs,massofZ1,massofZ2);
	   
	   hM4l_EEMM_I_I->Fill(massofhiggs,newweight);
	   hMZ1_EEMM_I_I->Fill(massofZ1,newweight);
	   hMZ2_EEMM_I_I->Fill(massofZ2,newweight);
	   
	 }
	 
       } //end of ( Z1tag==2 && Z2tag==1 )



     //cout<<"finished"<<endl;
     //cout<<"final massofZ1= "<<massofZ1<<" with tag = "<<Z1.tag<<endl;
     //cout<<"final massofZ2= "<<massofZ2<<" with tag = "<<Z2.tag<<endl;
     //cout<<"final mass4l= "<<final_mass4l<<endl;

     f_mass4l = massofhiggs;
     f_Z1mass = massofZ1;
     f_Z2mass = massofZ2;
     
	
	
     if(f_run2 != -999 && f_lumi2 != -999 && f_event2 != -999){
     //========= Gets the NNs ============
     /*
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
     
	   std::vector<std::vector<double> > inputs = {
	     //k41nj2
	     {l1pt,l1eta,l1phi,l2pt,l2eta,l2phi,l3pt,l3eta,l3phi,l4pt,l4eta,l4phi,j1pt,j1eta,j1phi,j2pt,j2eta,j2phi,met},
	     //k1177nj2
	     {l1pt,l1eta,l1phi,l2pt,l2eta,l2phi,l3pt,l3eta,l3phi,l4pt,l4eta,l4phi,j1pt,j1eta,j1phi,j2pt,j2eta,j2phi,met},
	     //k2nj3
	     {l1pt,l1eta,l1phi,l2pt,l2eta,l2phi,l3pt,l3eta,l3phi,l4pt,l4eta,l4phi,j1pt,j1eta,j1phi,j2pt,j2eta,j2phi,j3pt,j3eta,j3phi},
	     //k4nj3
	     {l1pt,l1eta,l1phi,l2pt,l2eta,l2phi,l3pt,l3eta,l3phi,l4pt,l4eta,l4phi,j1pt,j1eta,j1phi,j2pt,j2eta,j2phi,j3pt,j3eta,j3phi},
	     //k16nj2e3
	     {l1pt,l1eta,l1phi,l2pt,l2eta,l2phi,l3pt,l3eta,l3phi,l4pt,l4eta,l4phi,j1pt,j1eta,j1phi,j2pt,j2eta,j2phi,j3pt,j3eta,j3phi},
	     //k53nj2e3
	     {l1pt,l1eta,l1phi,l2pt,l2eta,l2phi,l3pt,l3eta,l3phi,l4pt,l4eta,l4phi,j1pt,j1eta,j1phi,j2pt,j2eta,j2phi,j3pt,j3eta,j3phi,met},
	     //k196nj2e3
	     {l1pt,l1eta,l1phi,l2pt,l2eta,l2phi,l3pt,l3eta,l3phi,l4pt,l4eta,l4phi,j1pt,j1eta,j1phi,j2pt,j2eta,j2phi,j3pt,j3eta,j3phi},
	     //k543nj2e3
	     {l1pt,l1eta,l1phi,l2pt,l2eta,l2phi,l3pt,l3eta,l3phi,l4pt,l4eta,l4phi,j1pt,j1eta,j1phi,j2pt,j2eta,j2phi,j3pt,j3eta,j3phi,met,njets,nbjets}
	   };
	   
	   f_k41nj2[ishift] = model_k41nj2(inputs[0]);
	   f_k1177nj2[ishift] = model_k1177nj2(inputs[1]);
	   f_k2nj3[ishift] = model_k2nj3(inputs[2]);
	   f_k4nj3[ishift] = model_k4nj3(inputs[3]);
	   f_k16nj2e3[ishift] = model_k16nj2e3(inputs[4]);
	   f_k53nj2e3[ishift] = model_k53nj2e3(inputs[5]);
	   f_k196nj2e3[ishift] = model_k196nj2e3(inputs[6]);
	   f_k543nj2e3[ishift] = model_k543nj2e3(inputs[7]);
	   
	   ++ishift;
         }
       }
       */
     
       //MELA
       vector<TLorentzVector> partP;
       partP.push_back(L11P4);
       partP.push_back(L12P4);
       partP.push_back(L21P4);
       partP.push_back(L22P4);
       
       vector<int> partId;
       partId.push_back(L11PID);
       partId.push_back(L12PID);
       partId.push_back(L21PID);
       partId.push_back(L22PID);
       
       
       vector<TLorentzVector> partPprod;
       vector<int> partIdprod;
       partPprod.push_back(L11P4);
       partPprod.push_back(L12P4);
       partPprod.push_back(L21P4);
       partPprod.push_back(L22P4);
       if (JET1.Pt()>0.) partPprod.push_back(JET1);
       if (JET2.Pt()>0.) partPprod.push_back(JET2); // Can also use partPprod.push_back(nullFourVector) instead for integrated VBF MEs
       
       partIdprod.push_back(L11PID);
       partIdprod.push_back(L12PID);
       partIdprod.push_back(L21PID);
       partIdprod.push_back(L22PID);
       if (JET1.Pt()>0.) partIdprod.push_back(0);
       if (JET2.Pt()>0.) partIdprod.push_back(0); // For leptonic ZH in the future, this could actually be the opposite lepton flavor
       

       double p0plus_VAJHU, bkg_VAMCFM, p0plus_m4l, bkg_m4l;
       double p0minus_VAJHU, Dgg10_VAMCFM;
       double phjj_VAJHU, pvbf_VAJHU;
       
       int a=combinedMEM.computeME(MEMNames::kSMHiggs, MEMNames::kJHUGen, partP, partId, p0plus_VAJHU); // Calculation of SM gg->H->4l JHUGen ME      
       cout << "a= "  << p0plus_VAJHU << endl;
       int b=combinedMEM.computeME(MEMNames::k0minus, MEMNames::kJHUGen, partP, partId, p0minus_VAJHU); // Calculation of PS (0-, fa3=1) gg->H->4l JHUGen ME 
       int c=combinedMEM.computeME(MEMNames::kggHZZ_10, MEMNames::kMCFM, partP, partId, Dgg10_VAMCFM); // Direct calculation of Dgg (D^kin for off-shell) from MCFM MEs
       int d=combinedMEM.computeME(MEMNames::kqqZZ, MEMNames::kMCFM, partP, partId, bkg_VAMCFM); // qq->4l background calculation from MCFM
       combinedMEM.computePm4l(partP,partId, MEMNames::kNone, p0plus_m4l, bkg_m4l); // m4l probabilities for signal and background, nominal resolution
       if (njets_pass>=2){
	 int f=combinedMEM.computeME(MEMNames::kJJ_SMHiggs_GG, MEMNames::kJHUGen, partPprod, partIdprod, phjj_VAJHU); // SM gg->H+2j
	 int g=combinedMEM.computeME(MEMNames::kJJ_SMHiggs_VBF, MEMNames::kJHUGen, partPprod, partIdprod, pvbf_VAJHU);  // SM VBF->H
       }

       //f_D_bkg_kin = p0plus_VAJHU / ( p0plus_VAJHU + bkg_VAMCFM ); // D^kin_bkg
       //f_D_bkg = p0plus_VAJHU * p0plus_m4l / ( p0plus_VAJHU * p0plus_m4l + bkg_VAMCFM * bkg_m4l ); // D^kin including superMELA
       //f_D_gg = Dgg10_VAMCFM; // D_gg
       //f_D_g4 = p0plus_VAJHU / ( p0plus_VAJHU + p0minus_VAJHU ); // D_0-
       if (njets_pass>=2) f_Djet_VAJHU = pvbf_VAJHU / ( pvbf_VAJHU + phjj_VAJHU ); // D^VBF_HJJ
       else f_Djet_VAJHU=-1;
              
       
              
       
       /*
       cout<< "***********************************************************"
       << "\nRun:    " << f_run2
       << "\nLumi:   " << f_lumi2
       << "\nEvent:  " << f_event2
       << "\nWeight: " << f_weight2
       <<endl;
       */
       treeCR->Fill();
     }
     
     
     //cout<<" end of event analysis...."<<endl;
          
   }//end Loop for jentry the main loop of the program for (Long64_t jentry=0; jentry<nentries;jentry++)
   //cout<<" end of loop over events..."<<endl;
   
   
   //========================================================
   //      Write histograms and output files goes here      //
   //==========================================================
   
   //___________________________//
   //cout<<" write on output files "<<theFile->GetName()<<endl;
   //_________________________//
   
   
   theFile->cd();
   
   //Write TTree
   treeFR->Write();
   treeCR->Write();
   
   
   //___________________________//
   //     write histograms     //
   //_________________________//
   
   histo_EVENT->Write();
   
   //loose Lepton Identification
   //============================
   
   //Muons
   //------
   hN_loose_mu->Write();
   hIso_loose_mu->Write();
   hSip_loose_mu->Write();
   hIp_loose_mu->Write();
   
   //Electrons
   //----------
   
   hN_loose_e->Write();
   hIso_loose_e->Write();
   hSip_loose_e->Write();
   hIp_loose_e->Write();
   
   //Good Lepton Identification
   //============================
   
   hN_good_mu->Write();
   hN_good_ele->Write();
   hN_good_lep->Write();
   hN_good_phot->Write();
   
   //Step3 
   //=======
   
   
   //Step 3
   //-------
   //after choosing Z1 and before applying FR
   
   //mumu
   hMZ1_3_mumu->Write();
   hPtZ1_3_mumu->Write();
   hYZ1_3_mumu->Write();
   
   //ee
   hMZ1_3_ee->Write();
   hPtZ1_3_ee->Write();
   hYZ1_3_ee->Write();
   
   hMZ1_3->Write();
   hPtZ1_3->Write();
   hPFMET_3->Write();
   
   //hMZ_3_check->Write();
   
   hMZ1_3_no_effW->Write();   
   hPFMET_3_no_effW->Write();

   //No MET cut
   //////////////

   //Muon inclusive
   
   hMZ1_loose_mu_barrel->Write(); 
   hPFMET_loose_mu_barrel->Write();
   hPT_loose_mu_barrel->Write();
  
   hMZ1_loose_mu_endcap->Write(); 
   hPFMET_loose_mu_endcap->Write();
   hPT_loose_mu_endcap->Write();
   
   hMZ1_tight_mu_barrel->Write(); 
   hPFMET_tight_mu_barrel->Write();
   hPT_tight_mu_barrel->Write();
   
   hMZ1_tight_mu_endcap->Write(); 
   hPFMET_tight_mu_endcap->Write();
   hPT_tight_mu_endcap->Write();

   //positive , negative

   h_PT_loose_mu_barrel_pos->Write();
   h_PFMET_loose_mu_barrel_pos->Write();
   h_PT_loose_mu_barrel_neg->Write();
   h_PFMET_loose_mu_barrel_neg->Write();
   
   h_PT_loose_mu_endcap_pos->Write(); 
   h_PFMET_loose_mu_endcap_pos->Write();
   h_PT_loose_mu_endcap_neg->Write();
   h_PFMET_loose_mu_endcap_neg->Write();
   
   h_PT_tight_mu_barrel_pos->Write();
   h_PFMET_tight_mu_barrel_pos->Write();
   h_PT_tight_mu_barrel_neg->Write();
   h_PFMET_tight_mu_barrel_neg->Write(); 
     
   h_PT_tight_mu_endcap_pos->Write();
   h_PFMET_tight_mu_endcap_pos->Write();
   h_PT_tight_mu_endcap_neg->Write();
   h_PFMET_tight_mu_endcap_neg->Write();
   
   
   //Electrons inclusive
   
   hMZ1_loose_ele_barrel->Write(); 
   hPFMET_loose_ele_barrel->Write();
   hPT_loose_ele_barrel->Write();
   
   hMZ1_loose_ele_endcap->Write(); 
   hPFMET_loose_ele_endcap->Write();
   hPT_loose_ele_endcap->Write();
   
   hMZ1_tight_ele_barrel->Write(); 
   hPFMET_tight_ele_barrel->Write();
   hPT_tight_ele_barrel->Write();
   
   hMZ1_tight_ele_endcap->Write(); 
   hPFMET_tight_ele_endcap->Write();
   hPT_tight_ele_endcap->Write();

   //Positive ,negative

   h_PT_loose_ele_barrel_pos->Write();
   h_PFMET_loose_ele_barrel_pos->Write();
   h_PT_loose_ele_barrel_neg->Write();
   h_PFMET_loose_ele_barrel_neg->Write();
   
   h_PT_loose_ele_endcap_pos->Write(); 
   h_PFMET_loose_ele_endcap_pos->Write();
   h_PT_loose_ele_endcap_neg->Write();
   h_PFMET_loose_ele_endcap_neg->Write();
   
   h_PT_tight_ele_barrel_pos->Write();
   h_PFMET_tight_ele_barrel_pos->Write();
   h_PT_tight_ele_barrel_neg->Write();
   h_PFMET_tight_ele_barrel_neg->Write(); 
     
  h_PT_tight_ele_endcap_pos->Write();
  h_PFMET_tight_ele_endcap_pos->Write();
  h_PT_tight_ele_endcap_neg->Write();
  h_PFMET_tight_ele_endcap_neg->Write();
  

   
   //2D
  ///////

   MET_PT_Ele_Barrel_Den->Write();
   MET_PT_Ele_Barrel_Num->Write();
   MET_PT_Ele_Endcap_Den->Write();
   MET_PT_Ele_Endcap_Num->Write();

   MET_PT_Mu_Barrel_Den->Write();
   MET_PT_Mu_Barrel_Num->Write();
   MET_PT_Mu_Endcap_Den->Write();
   MET_PT_Mu_Endcap_Num->Write();
   
   //Rebin
   
   MET_PT_Ele_Barrel_Den_Rebin->Write();
   MET_PT_Ele_Barrel_Num_Rebin->Write();
   MET_PT_Ele_Endcap_Den_Rebin->Write();
   MET_PT_Ele_Endcap_Num_Rebin->Write();

   MET_PT_Mu_Barrel_Den_Rebin->Write();
   MET_PT_Mu_Barrel_Num_Rebin->Write();
   MET_PT_Mu_Endcap_Den_Rebin->Write();
   MET_PT_Mu_Endcap_Num_Rebin->Write();

   //Positive , Negative

   h_MET_PT_loose_mu_barrel_pos->Write();
   h_MET_PT_loose_mu_barrel_neg->Write();
   
   h_MET_PT_loose_mu_endcap_pos->Write();
   h_MET_PT_loose_mu_endcap_neg->Write();
   
   h_MET_PT_tight_mu_barrel_pos->Write();
   h_MET_PT_tight_mu_barrel_neg->Write();
   
   h_MET_PT_tight_mu_endcap_pos->Write();
   h_MET_PT_tight_mu_endcap_neg->Write();
   //  
   h_MET_PT_loose_ele_barrel_pos->Write();
   h_MET_PT_loose_ele_barrel_neg->Write();
   
   h_MET_PT_loose_ele_endcap_pos->Write();
   h_MET_PT_loose_ele_endcap_neg->Write();
   
   h_MET_PT_tight_ele_barrel_pos->Write();
   h_MET_PT_tight_ele_barrel_neg->Write();
   
   h_MET_PT_tight_ele_endcap_pos->Write();
   h_MET_PT_tight_ele_endcap_neg->Write();
   
   
   
   
   //Fake rate after MET cut
   //-------------
   //1) Muon Fake rate 
   //==================
   ZplusM_Pt_DEN_Barrel->Write();
   ZplusM_Pt_DEN_Endcaps->Write();
   ZplusM_Pt_NUM_ID_ISO_Barrel->Write();
   ZplusM_Pt_NUM_ID_ISO_Endcaps->Write();
   
  //Positive
   
   ZplusM_Pt_DEN_Barrel_pos->Write();      
   ZplusM_Pt_DEN_Endcaps_pos->Write();  
   ZplusM_Pt_NUM_ID_ISO_Barrel_pos->Write();  
   ZplusM_Pt_NUM_ID_ISO_Endcaps_pos->Write();
   

  //Negative

   ZplusM_Pt_DEN_Barrel_neg->Write();
   ZplusM_Pt_DEN_Endcaps_neg->Write();
   ZplusM_Pt_NUM_ID_ISO_Barrel_neg->Write();
   ZplusM_Pt_NUM_ID_ISO_Endcaps_neg->Write();
   
   
   //2) Electron Fake rate
   //======================
   ZplusE_Pt_DEN_Barrel->Write();
   ZplusE_Pt_DEN_Endcaps->Write();
   ZplusE_Pt_NUM_ID_ISO_Barrel->Write();
   ZplusE_Pt_NUM_ID_ISO_Endcaps->Write();

  //Positive
   
   ZplusE_Pt_DEN_Barrel_pos->Write();      
   ZplusE_Pt_DEN_Endcaps_pos->Write();  
   ZplusE_Pt_NUM_ID_ISO_Barrel_pos->Write();  
   ZplusE_Pt_NUM_ID_ISO_Endcaps_pos->Write();
   

  //Negative

   ZplusE_Pt_DEN_Barrel_neg->Write();
   ZplusE_Pt_DEN_Endcaps_neg->Write();
   ZplusE_Pt_NUM_ID_ISO_Barrel_neg->Write();
   ZplusE_Pt_NUM_ID_ISO_Endcaps_neg->Write();
   
   //control regions
   
   //MMMM
   
   hM4l_MMMM_I_I ->Write();
   hMZ1_MMMM_I_I ->Write();
   hMZ2_MMMM_I_I ->Write();
   
   
   hM4l_MMMM_AI_AI ->Write();
   hM4l_MMMM_AI_I ->Write();
   hM4l_MMMM_I_AI ->Write();
   hM4l_MMMM_3p1f_total->Write();

   hPFMET_MMMM_AI_AI->Write();
   hPFMET_MMMM_AI_I->Write();
   hPFMET_MMMM_I_AI->Write();
   hPFMET_MMMM_3p1f_total->Write();

   //for Simran

   hPt_4mu_AI_AI_1->Write();
   hPt_4mu_AI_AI_2->Write();
   hPt_4mu_AI_AI_3->Write();
   
   //MMEE
   
   hM4l_MMEE_AI_AI->Write();
   hM4l_MMEE_AI_I->Write();
   hM4l_MMEE_I_AI->Write();
   hM4l_MMEE_3p1f_total->Write(); 
   
   hM4l_MMEE_I_I->Write();
   hMZ1_MMEE_I_I->Write();
   hMZ2_MMEE_I_I->Write();

   hPFMET_MMEE_AI_AI->Write();
   hPFMET_MMEE_AI_I->Write();
   hPFMET_MMEE_I_AI->Write();
   hPFMET_MMEE_3p1f_total->Write();
   
   //EEMM
   
   hM4l_EEMM_AI_AI->Write();
   hM4l_EEMM_AI_I->Write();
   hM4l_EEMM_I_AI->Write();
   hM4l_EEMM_3p1f_total->Write();
   
   hM4l_EEMM_I_I->Write();
   hMZ1_EEMM_I_I->Write();
   hMZ2_EEMM_I_I->Write();

   hPFMET_EEMM_AI_AI->Write();
   hPFMET_EEMM_AI_I->Write();
   hPFMET_EEMM_I_AI->Write();
   hPFMET_EEMM_3p1f_total->Write();
   
   //EEEE
   
   hM4l_EEEE_I_I->Write();
   hMZ1_EEEE_I_I->Write();
   hMZ2_EEEE_I_I->Write();
   hM4l_EEEE_3p1f_total->Write();  
   
   hM4l_EEEE_AI_AI->Write();
   hM4l_EEEE_AI_I->Write();
   hM4l_EEEE_I_AI->Write();
   
   hPFMET_EEEE_AI_AI->Write();
   hPFMET_EEEE_AI_I->Write();
   hPFMET_EEEE_I_AI->Write();
   hPFMET_EEEE_3p1f_total->Write();


     //around Z peak in 4l

   hPFMET_MMMM_AI_AI_Zpeak->Write(); 
   hPFMET_MMMM_AI_I_Zpeak->Write();
   hPFMET_MMMM_I_AI_Zpeak->Write(); 
   hPFMET_MMMM_3p1f_total_Zpeak->Write();
 
	    
   hPFMET_EEMM_AI_AI_Zpeak->Write();  
   hPFMET_EEMM_AI_I_Zpeak->Write(); 
   hPFMET_EEMM_I_AI_Zpeak->Write();  
   hPFMET_EEMM_3p1f_total_Zpeak->Write();
   

   hPFMET_EEEE_AI_AI_Zpeak->Write();
   hPFMET_EEEE_AI_I_Zpeak->Write(); 
   hPFMET_EEEE_I_AI_Zpeak->Write();   
   hPFMET_EEEE_3p1f_total_Zpeak->Write();
   
   
   hPFMET_MMEE_AI_AI_Zpeak->Write();   
   hPFMET_MMEE_AI_I_Zpeak->Write();   
   hPFMET_MMEE_I_AI_Zpeak->Write();  
   hPFMET_MMEE_3p1f_total_Zpeak->Write();
   
   
   /////
   
   hPFMET_MMMM_AI_AI_Zpeak_log->Write();      
   hPFMET_MMMM_AI_I_Zpeak_log->Write();
   hPFMET_MMMM_I_AI_Zpeak_log->Write();   
   hPFMET_MMMM_3p1f_total_Zpeak_log->Write();
   
   
   hPFMET_EEMM_AI_AI_Zpeak_log->Write();      
   hPFMET_EEMM_AI_I_Zpeak_log->Write();    
   hPFMET_EEMM_I_AI_Zpeak_log->Write();  
   hPFMET_EEMM_3p1f_total_Zpeak_log->Write();
   
   
   hPFMET_EEEE_AI_AI_Zpeak_log->Write();     
   hPFMET_EEEE_AI_I_Zpeak_log->Write();      
   hPFMET_EEEE_I_AI_Zpeak_log->Write();   
   hPFMET_EEEE_3p1f_total_Zpeak_log->Write();
   

   hPFMET_MMEE_AI_AI_Zpeak_log->Write();      
   hPFMET_MMEE_AI_I_Zpeak_log->Write();      
   hPFMET_MMEE_I_AI_Zpeak_log->Write();   
   hPFMET_MMEE_3p1f_total_Zpeak_log->Write();

   
   /////////////////////////////////////////////

      //control regions Step 2
   
   //MMMM
   
   hM4l_MMMM_I_I_Step2 ->Write();
   hMZ1_MMMM_I_I_Step2 ->Write();
   hMZ2_MMMM_I_I_Step2 ->Write();
   
   
   hM4l_MMMM_AI_AI_Step2 ->Write();
   hM4l_MMMM_AI_I_Step2 ->Write();
   hM4l_MMMM_I_AI_Step2 ->Write();
   hM4l_MMMM_3p1f_total_Step2->Write();

   hPFMET_MMMM_AI_AI_Step2->Write();
   hPFMET_MMMM_AI_I_Step2 ->Write();
   hPFMET_MMMM_I_AI_Step2 ->Write();
   hPFMET_MMMM_3p1f_total_Step2->Write();
   
   //MMEE
   
   hM4l_MMEE_AI_AI_Step2->Write();
   hM4l_MMEE_AI_I_Step2->Write();
   hM4l_MMEE_I_AI_Step2->Write();
   hM4l_MMEE_3p1f_total_Step2->Write(); 
   
   hM4l_MMEE_I_I_Step2->Write();
   hMZ1_MMEE_I_I_Step2->Write();
   hMZ2_MMEE_I_I_Step2->Write();

   hPFMET_MMEE_AI_AI_Step2->Write();
   hPFMET_MMEE_AI_I_Step2 ->Write();
   hPFMET_MMEE_I_AI_Step2 ->Write();
   hPFMET_MMEE_3p1f_total_Step2->Write();
   
   //EEMM
   
   hM4l_EEMM_AI_AI_Step2->Write();
   hM4l_EEMM_AI_I_Step2->Write();
   hM4l_EEMM_I_AI_Step2->Write();
   hM4l_EEMM_3p1f_total_Step2->Write();
   
   hM4l_EEMM_I_I_Step2->Write();
   hMZ1_EEMM_I_I_Step2->Write();
   hMZ2_EEMM_I_I_Step2->Write();

   hPFMET_EEMM_AI_AI_Step2->Write();
   hPFMET_EEMM_AI_I_Step2 ->Write();
   hPFMET_EEMM_I_AI_Step2 ->Write();
   hPFMET_EEMM_3p1f_total_Step2->Write();
   
   //EEEE
   
   hM4l_EEEE_I_I_Step2->Write();
   hMZ1_EEEE_I_I_Step2->Write();
   hMZ2_EEEE_I_I_Step2->Write();
   hM4l_EEEE_3p1f_total_Step2->Write();  
   
   hM4l_EEEE_AI_AI_Step2->Write();
   hM4l_EEEE_AI_I_Step2->Write();
   hM4l_EEEE_I_AI_Step2->Write();

   hPFMET_EEEE_AI_AI_Step2->Write();
   hPFMET_EEEE_AI_I_Step2 ->Write();
   hPFMET_EEEE_I_AI_Step2 ->Write();
   hPFMET_EEEE_3p1f_total_Step2->Write();


   //In MonoHiggs step for systematic after revert M4l cut
   //MMMM

    hPFMET_MMMM_AI_AI_Step2_rev->Write();
    hPFMET_MMMM_AI_I_Step2_rev->Write();
    hPFMET_MMMM_I_AI_Step2_rev->Write();
    hPFMET_MMMM_3p1f_total_Step2_rev->Write();

   //MMEE

   hPFMET_MMEE_AI_AI_Step2_rev->Write(); 
   hPFMET_MMEE_AI_I_Step2_rev->Write(); 
   hPFMET_MMEE_I_AI_Step2_rev->Write(); 
   hPFMET_MMEE_3p1f_total_Step2_rev->Write(); 

   //EEMM

   hPFMET_EEMM_AI_AI_Step2_rev->Write();
   hPFMET_EEMM_AI_I_Step2_rev->Write(); 
   hPFMET_EEMM_I_AI_Step2_rev->Write(); 
   hPFMET_EEMM_3p1f_total_Step2_rev->Write();

   //EEEE

   hPFMET_EEEE_AI_AI_Step2_rev->Write();
   hPFMET_EEEE_AI_I_Step2_rev->Write(); 
   hPFMET_EEEE_I_AI_Step2_rev->Write();
   hPFMET_EEEE_3p1f_total_Step2_rev->Write();
   
   theFile->Close();
   
   //cout<<"finished"<<endl;

   return;
}//end of main file



//Main object
int main(int argc, char ** argv){
  TString dataset_path = argv[1];
  TString datasetName = argv[2];
  TString store_path = argv[3];
  TString sXsection = argv[4];
  double Xsection = sXsection.Atof();
  TString sneventsPreHLT = argv[5];
  int neventsPreHLT = sneventsPreHLT.Atoi();
  TString sneventsPostHLT = argv[6];
  int neventsPostHLT = sneventsPostHLT.Atoi();
  TString DATA_type = argv[7];
  TString MC_type = argv[8];
  TString sNevents = argv[9];
  int Nevents = sNevents.Atoi();
  std::cout<<">>>>>>>>>>>>> Looking at file "<<datasetName<<" <<<<<<<<<<<<<<<<<"<<std::endl;
  ComputeFRandCRsOS(dataset_path,datasetName,store_path,Xsection,neventsPreHLT,neventsPostHLT,DATA_type,MC_type,Nevents);
  std::cout<<">>>>>>>>>>>>>>>> Finished <<<<<<<<<<<<<<<<<"<<std::endl;
  
  
  return 0;
}
