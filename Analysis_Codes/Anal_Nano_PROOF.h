//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr 15 13:52:58 2019 by ROOT version 5.34/38
// from TTree Events/Events
// found on file: root://se01.indiacms.res.in//store/user/chatterj/NanoPost/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/NanoTestPost_2017TTBar/190408_151813/0000/tree_6.root
//////////////////////////////////////////////////////////

#ifndef Anal_Nano_PROOF_h
#define Anal_Nano_PROOF_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <math.h>
#include "TLorentzVector.h"
#include <TProofOutputFile.h>
#include <TProofServ.h>

#include <string>
#include <cmath>
#include <cassert>
#include <sstream>
#include <string>
#include "TFileCollection.h"
#include "THashList.h"
#include "TBenchmark.h"

#include <iostream>
#include <fstream>

#include "TSystem.h"
#include "boost/config.hpp"
#include "boost/lexical_cast.hpp"

#include "TMatrixDBase.h"
#include "TMatrixD.h"

// Header file for the classes stored in the TTree if any.

//const double M_PI = acos(-1);
/*
int* dec2bin(int dec, int length)
{
//const int length = nHLTmx;
int input = dec;
int istep=0;
int div[length];
while(input>0){
  input = input/2;
  div[istep] = input%2;
  istep++;
 }
*/

int* decToBinary(int n)
{
	const int length = 15;
    static int binaryNum[length];
    int i = 0;
    while (n > 0) {
        binaryNum[i] = n % 2;
        n = n / 2;
        i++;
    }
    int sum=0;
    
    return binaryNum;
    /*
    if(i>10) { i = 10; }
    // printing binary array in reverse order
    for (int j = i - 1; j >= 0; j--) { sum += binaryNum[j]*pow(10,j); }
    return sum;
    */ 
} 

int getbinid(double val, int nbmx, float* array) {
  if (val<array[0]) return -2;
  for (int ix=0; ix<=nbmx; ix++) {
    if (val < array[ix]) return ix-1;
  }
  return -3;
}

int getbinid(double val, int nbmx, double* array) {
  if (val<array[0]) return -2;
  for (int ix=0; ix<=nbmx; ix++) {
    if (val < array[ix]) return ix-1;
  }
  return -3;
}

double theta_to_eta(double theta) { return -log(tan(theta/2.)); }

double eta_to_theta(double eta){
return(2*atan(exp(-2*eta)));
}

double PhiInRange(const double& phi) {
  double phiout = phi;

  if( phiout > 2*M_PI || phiout < -2*M_PI) {
    phiout = fmod( phiout, 2*M_PI);
  }
  if (phiout <= -M_PI) phiout += 2*M_PI;
  else if (phiout >  M_PI) phiout -= 2*M_PI;

  return phiout;
}

double delta2R(double eta1, double phi1, double eta2, double phi2) {
  return sqrt(pow(eta1 - eta2,2) +pow(PhiInRange(phi1 - phi2),2));
}

double bmasscut_fun(double pt){
//return (0.5*pt-200);	
return 60;
}

double EW_toppt_cor(double pt){
return (exp(-1.08872-(pt*0.011998)) + 0.895139);
}

Double_t pol4(Double_t* x, Double_t* par){
return (par[4]*pow(x[0],4)+par[3]*pow(x[0],3)+par[2]*pow(x[0],2)+par[1]*pow(x[0],1)+par[0]);
}

Double_t pol3(Double_t* x, Double_t* par){
return (par[3]*pow(x[0],3)+par[2]*pow(x[0],2)+par[1]*pow(x[0],1)+par[0]);
}

Double_t pol2(Double_t* x, Double_t* par){
return (par[2]*pow(x[0],2)+par[1]*pow(x[0],1)+par[0]);
}

Double_t pol1(Double_t* x, Double_t* par){
return (par[1]*pow(x[0],1)+par[0]);
}

double Pol0(double* x, double* par){
return (par[0]);
}

double Pol1(double* x, double* par){
return (par[0]+par[1]*x[0]);
}

double MikkoFunc1(double *x, double *par){
return (par[3]+par[0]*pow(x[0],par[1])+par[2]*log(x[0])/x[0]);
}

double Parabol(double* x, double* par){
double xx = x[0]-par[0] ;
if(xx<0) {return(par[1]+(par[2]*pow(xx,2)));}
else { return(par[1]+(par[3]*pow(xx,2))) ; }
}

double BiFun(double* x, double* par){
double xx = x[0]-par[0] ;
if(xx<0) {return(par[1]+(par[2]*xx)+(par[3]*pow(xx,2)));}
else { return(par[1]+(par[2]*xx)+(par[4]*pow(xx,2))) ; }
}

double bdF0(double* x, double* par){
double xx = x[0]-par[0] ;
if(xx<0) {return(-par[2]-2*par[3]*xx);}
else {return(-par[2]-2*par[4]*xx);}
}
double bdF1(double* x, double* par){
return 1;
}
double bdF2(double* x, double* par){
return (x[0]-par[0]);
}
double bdF3(double* x, double* par){
return ((x[0]-par[0])*(x[0]-par[0]));
}
double bdF4(double* x, double* par){
return ((x[0]-par[0])*(x[0]-par[0]));
}


// Fixed size dimensions of array or collections stored in the TTree if any.

double BTag_MCEfficiency(string algo, int flavor, double pt, double eta){

if ((abs(flavor) ==5)){
	if (algo=="deepcsvm") {
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  20) && (pt < 30))) return 0.547121;
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  20) && (pt < 30))) return 0.550118;
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  20) && (pt < 30))) return 0.462942;
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  30) && (pt < 50))) return 0.670671;
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  30) && (pt < 50))) return 0.674576;
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  30) && (pt < 50))) return 0.61123;
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  50) && (pt < 70))) return 0.719526;
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  50) && (pt < 70))) return 0.720965;
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  50) && (pt < 70))) return 0.660123;
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  70) && (pt < 100))) return 0.73676;
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  70) && (pt < 100))) return 0.73467;
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  70) && (pt < 100))) return 0.671906;
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  100) && (pt < 140))) return 0.747847;
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  100) && (pt < 140))) return 0.740672;
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  100) && (pt < 140))) return 0.667903;
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  140) && (pt < 200))) return 0.744575;
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  140) && (pt < 200))) return 0.731165;
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  140) && (pt < 200))) return 0.645684;
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  200) && (pt < 300))) return 0.702543;
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  200) && (pt < 300))) return 0.681865;
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  200) && (pt < 300))) return 0.585971;
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  300) && (pt < 600))) return 0.589448;
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  300) && (pt < 600))) return 0.563962;
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  300) && (pt < 600))) return 0.472939;
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  600) && (pt < 10000))) return 0.400163;
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  600) && (pt < 10000))) return 0.376463;
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  600) && (pt < 10000))) return 0.317601;
		}
	}

if ((abs(flavor)!= 5) && (abs(flavor)!= 4)){
	if (algo=="deepcsvm") {
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  20) && (pt < 30))) return 0.0245662;
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  20) && (pt < 30))) return 0.02644;
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  20) && (pt < 30))) return 0.0232145;
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  30) && (pt < 50))) return 0.0240801;
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  30) && (pt < 50))) return 0.02682855;
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  30) && (pt < 50))) return 0.0256761;
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  50) && (pt < 70))) return 0.0230801;
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  50) && (pt < 70))) return 0.0251746;
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  50) && (pt < 70))) return 0.0257851;
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  70) && (pt < 100))) return 0.0235181;
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  70) && (pt < 100))) return 0.0260925;
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  70) && (pt < 100))) return 0.0263857;
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  100) && (pt < 140))) return 0.0237109;
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  100) && (pt < 140))) return 0.0262221;
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  100) && (pt < 140))) return 0.0287912;
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  140) && (pt < 200))) return 0.0229224;
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  140) && (pt < 200))) return 0.0250336;
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  140) && (pt < 200))) return 0.0246814;
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  200) && (pt < 300))) return 0.0191271;
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  200) && (pt < 300))) return 0.0201025;
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  200) && (pt < 300))) return 0.0200517;
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  300) && (pt < 600))) return 0.0125901;
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  300) && (pt < 600))) return 0.0118787;
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  300) && (pt < 600))) return 0.0136332;
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  600) && (pt < 10000))) return 0.00800836;
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  600) && (pt < 10000))) return 0.0063758;
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  600) && (pt < 10000))) return 0.00793709;
      }	
	}

return 1.0;

}

class Anal_Nano_PROOF : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  
   static const int njetmax = 25;
   static const int njesmax = 37;
   static const int nbtagmax = 9;
   static const int npartonmax = 200;
   static const int nSVmax = 100;
   static const int nTrigObjmax = 100;
   static const int npdfmax = 33;
   static const int nscalemax = 9;
   static const int npsmax = 4;
   static const int npartmx = 100;

   float weight;
   float weight_t;
   float SF_toppt;
   double weightev, weightpass, weightpass_btag;
   float btag_weight;
   float btag_weight_unc_up[nbtagmax], btag_weight_unc_dn[nbtagmax];
   int nevent_total;
   bool isMC;
   bool isSignal;
   bool isQCD;
   bool isTTBar;
   bool isTTHad, isTTSemiLep, isTTDiLep;
   bool isST;
   bool FakeAnalysis;
   
   char name[100];
   char title[100];
   
   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   Float_t         btagWeight_CSVV2;
   Float_t         btagWeight_DeepCSVB;
   Float_t         CaloMET_phi;
   Float_t         CaloMET_pt;
   Float_t         CaloMET_sumEt;
   Float_t         ChsMET_phi;
   Float_t         ChsMET_pt;
   Float_t         ChsMET_sumEt;
   UInt_t          nElectron;
   Float_t         Electron_deltaEtaSC[njetmax];   //[nElectron]
   Float_t         Electron_dr03EcalRecHitSumEt[njetmax];   //[nElectron]
   Float_t         Electron_dr03HcalDepth1TowerSumEt[njetmax];   //[nElectron]
   Float_t         Electron_dr03TkSumPt[njetmax];   //[nElectron]
   Float_t         Electron_dr03TkSumPtHEEP[njetmax];   //[nElectron]
   Float_t         Electron_dxy[njetmax];   //[nElectron]
   Float_t         Electron_dxyErr[njetmax];   //[nElectron]
   Float_t         Electron_dz[njetmax];   //[nElectron]
   Float_t         Electron_dzErr[njetmax];   //[nElectron]
   Float_t         Electron_eInvMinusPInv[njetmax];   //[nElectron]
   Float_t         Electron_energyErr[njetmax];   //[nElectron]
   Float_t         Electron_eta[njetmax];   //[nElectron]
   Float_t         Electron_hoe[njetmax];   //[nElectron]
   Float_t         Electron_ip3d[njetmax];   //[nElectron]
   Float_t         Electron_jetRelIso[njetmax];   //[nElectron]
   Float_t         Electron_mass[njetmax];   //[nElectron]
   Float_t         Electron_miniPFRelIso_all[njetmax];   //[nElectron]
   Float_t         Electron_miniPFRelIso_chg[njetmax];   //[nElectron]
   Float_t         Electron_mvaFall17V1Iso[njetmax];   //[nElectron]
   Float_t         Electron_mvaFall17V1noIso[njetmax];   //[nElectron]
   Float_t         Electron_mvaFall17V2Iso[njetmax];   //[nElectron]
   Float_t         Electron_mvaFall17V2noIso[njetmax];   //[nElectron]
   Float_t         Electron_pfRelIso03_all[njetmax];   //[nElectron]
   Float_t         Electron_pfRelIso03_chg[njetmax];   //[nElectron]
   Float_t         Electron_phi[njetmax];   //[nElectron]
   Float_t         Electron_pt[njetmax];   //[nElectron]
   Float_t         Electron_r9[njetmax];   //[nElectron]
   Float_t         Electron_sieie[njetmax];   //[nElectron]
   Float_t         Electron_sip3d[njetmax];   //[nElectron]
   Float_t         Electron_mvaTTH[njetmax];   //[nElectron]
   Float_t         Electron_eCorr[njetmax];   //[nElectron]
   Int_t           Electron_charge[njetmax];   //[nElectron]
   Int_t           Electron_cutBased[njetmax];   //[nElectron]
   Int_t           Electron_cutBased_Fall17_V1[njetmax];   //[nElectron]
   Int_t           Electron_jetIdx[njetmax];   //[nElectron]
   Int_t           Electron_pdgId[njetmax];   //[nElectron]
   Int_t           Electron_photonIdx[njetmax];   //[nElectron]
   Int_t           Electron_tightCharge[njetmax];   //[nElectron]
   Int_t           Electron_vidNestedWPBitmap[njetmax];   //[nElectron]
   Bool_t          Electron_convVeto[njetmax];   //[nElectron]
   Bool_t          Electron_cutBased_HEEP[njetmax];   //[nElectron]
   Bool_t          Electron_isPFcand[njetmax];   //[nElectron]
   UChar_t         Electron_lostHits[njetmax];   //[nElectron]
   Bool_t          Electron_mvaFall17V1Iso_WP80[njetmax];   //[nElectron]
   Bool_t          Electron_mvaFall17V1Iso_WP90[njetmax];   //[nElectron]
   Bool_t          Electron_mvaFall17V1Iso_WPL[njetmax];   //[nElectron]
   Bool_t          Electron_mvaFall17V1noIso_WP80[njetmax];   //[nElectron]
   Bool_t          Electron_mvaFall17V1noIso_WP90[njetmax];   //[nElectron]
   Bool_t          Electron_mvaFall17V1noIso_WPL[njetmax];   //[nElectron]
   Bool_t          Electron_mvaFall17V2Iso_WP80[njetmax];   //[nElectron]
   Bool_t          Electron_mvaFall17V2Iso_WP90[njetmax];   //[nElectron]
   Bool_t          Electron_mvaFall17V2Iso_WPL[njetmax];   //[nElectron]
   Bool_t          Electron_mvaFall17V2noIso_WP80[njetmax];   //[nElectron]
   Bool_t          Electron_mvaFall17V2noIso_WP90[njetmax];   //[nElectron]
   Bool_t          Electron_mvaFall17V2noIso_WPL[njetmax];   //[nElectron]
   UChar_t         Flag_ecalBadCalibFilterV2;
   UInt_t          nFatJet;
   Float_t         FatJet_area[njetmax];   //[nFatJet]
   Float_t         FatJet_btagCMVA[njetmax];   //[nFatJet]
   Float_t         FatJet_btagCSVV2[njetmax];   //[nFatJet]
   Float_t         FatJet_btagDeepB[njetmax];   //[nFatJet]
   Float_t         FatJet_btagHbb[njetmax];   //[nFatJet]
   Float_t         FatJet_deepTagMD_H4qvsQCD[njetmax];   //[nFatJet]
   Float_t         FatJet_deepTagMD_HbbvsQCD[njetmax];   //[nFatJet]
   Float_t         FatJet_deepTagMD_TvsQCD[njetmax];   //[nFatJet]
   Float_t         FatJet_deepTagMD_WvsQCD[njetmax];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZHbbvsQCD[njetmax];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZHccvsQCD[njetmax];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZbbvsQCD[njetmax];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZvsQCD[njetmax];   //[nFatJet]
   Float_t         FatJet_deepTagMD_bbvsLight[njetmax];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ccvsLight[njetmax];   //[nFatJet]
   Float_t         FatJet_deepTag_H[njetmax];   //[nFatJet]
   Float_t         FatJet_deepTag_QCD[njetmax];   //[nFatJet]
   Float_t         FatJet_deepTag_QCDothers[njetmax];   //[nFatJet]
   Float_t         FatJet_deepTag_TvsQCD[njetmax];   //[nFatJet]
   Float_t         FatJet_deepTag_WvsQCD[njetmax];   //[nFatJet]
   Float_t         FatJet_deepTag_ZvsQCD[njetmax];   //[nFatJet]
   Float_t         FatJet_eta[njetmax];   //[nFatJet]
   Float_t         FatJet_mass[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop[njetmax];   //[nFatJet]
   Float_t         FatJet_n2b1[njetmax];   //[nFatJet]
   Float_t         FatJet_n3b1[njetmax];   //[nFatJet]
   Float_t         FatJet_phi[njetmax];   //[nFatJet]
   Float_t         FatJet_pt[njetmax];   //[nFatJet]
   Float_t         FatJet_rawFactor[njetmax];   //[nFatJet]
   Float_t         FatJet_tau1[njetmax];   //[nFatJet]
   Float_t         FatJet_tau2[njetmax];   //[nFatJet]
   Float_t         FatJet_tau3[njetmax];   //[nFatJet]
   Float_t         FatJet_tau4[njetmax];   //[nFatJet]
   Float_t		   FatJet_subbtagCSVV2[njetmax];
   Float_t		   FatJet_subbtagDeepB[njetmax];
   Int_t           FatJet_jetId[njetmax];   //[nFatJet]
   Int_t           FatJet_subJetIdx1[njetmax];   //[nFatJet]
   Int_t           FatJet_subJetIdx2[njetmax];   //[nFatJet]
   Int_t 		   FatJet_subbtagId_CSVv2[njetmax];  
   Int_t 		   FatJet_subbtagId_DeepCSV[njetmax];  
   Float_t 		   FatJet_pt_jes_up[njesmax][njetmax];  
   Float_t 		   FatJet_mass_jes_up[njesmax][njetmax];  
   Float_t 		   FatJet_msoftdrop_jes_up[njesmax][njetmax];  
   Float_t 		   FatJet_pt_jer_up[njetmax];  
   Float_t 		   FatJet_mass_jer_up[njetmax];  
   Float_t		   FatJet_msoftdrop_jer_up[njetmax];  
   Float_t 		   FatJet_mass_jms_up[njetmax];  
   Float_t		   FatJet_msoftdrop_jms_up[njetmax];
   Float_t 		   FatJet_mass_jmr_up[njetmax];  
   Float_t		   FatJet_msoftdrop_jmr_up[njetmax];
   Float_t 		   FatJet_pt_jes_dn[njesmax][njetmax];  
   Float_t 		   FatJet_mass_jes_dn[njesmax][njetmax];  
   Float_t 		   FatJet_msoftdrop_jes_dn[njesmax][njetmax];  
   Float_t 		   FatJet_pt_jer_dn[njetmax];  
   Float_t 		   FatJet_mass_jer_dn[njetmax]; 
   Float_t		   FatJet_msoftdrop_jer_dn[njetmax];   
   Float_t 		   FatJet_mass_jms_dn[njetmax];  
   Float_t		   FatJet_msoftdrop_jms_dn[njetmax];
   Float_t 		   FatJet_mass_jmr_dn[njetmax];  
   Float_t		   FatJet_msoftdrop_jmr_dn[njetmax];
   Int_t 		   FatJet_MatchJet[njetmax];
   Int_t		   FatJet_MatchGenJet[njetmax];
   Float_t 		   FatJet_GenJetpt[njetmax];
   Float_t 		   FatJet_GenJeteta[njetmax];
   Float_t 		   FatJet_GenJetphi[njetmax];
   Float_t 		   FatJet_GenJetmass[njetmax];
   Bool_t		   FatJet_hashtop[njetmax];
   UInt_t          nGenJetAK8;
   Float_t         GenJetAK8_eta[njetmax];   //[nGenJetAK8]
   Float_t         GenJetAK8_mass[njetmax];   //[nGenJetAK8]
   Float_t         GenJetAK8_phi[njetmax];   //[nGenJetAK8]
   Float_t         GenJetAK8_pt[njetmax];   //[nGenJetAK8]
   UInt_t          nGenJet;
   Float_t         GenJet_eta[njetmax];   //[nGenJet]
   Float_t         GenJet_mass[njetmax];   //[nGenJet]
   Float_t         GenJet_phi[njetmax];   //[nGenJet]
   Float_t         GenJet_pt[njetmax];   //[nGenJet]
   UInt_t          nGenPart;
   Float_t         GenPart_eta[npartonmax];   //[nGenPart]
   Float_t         GenPart_mass[npartonmax];   //[nGenPart]
   Float_t         GenPart_phi[npartonmax];   //[nGenPart]
   Float_t         GenPart_pt[npartonmax];   //[nGenPart]
   Int_t           GenPart_genPartIdxMother[npartonmax];   //[nGenPart]
   Int_t           GenPart_pdgId[npartonmax];   //[nGenPart]
   Int_t           GenPart_status[npartonmax];   //[nGenPart]
   Int_t           GenPart_statusFlags[npartonmax];   //[nGenPart]
   UInt_t          nSubGenJetAK8;
   Float_t         SubGenJetAK8_eta[njetmax];   //[nSubGenJetAK8]
   Float_t         SubGenJetAK8_mass[njetmax];   //[nSubGenJetAK8]
   Float_t         SubGenJetAK8_phi[njetmax];   //[nSubGenJetAK8]
   Float_t         SubGenJetAK8_pt[njetmax];   //[nSubGenJetAK8]
   Float_t         Generator_binvar;
   Float_t         Generator_scalePDF;
   Float_t         Generator_weight;
   Float_t         Generator_x1;
   Float_t         Generator_x2;
   Float_t         Generator_xpdf1;
   Float_t         Generator_xpdf2;
   Int_t           Generator_id1;
   Int_t           Generator_id2;
   UInt_t          nGenVisTau;
   Float_t         GenVisTau_eta[njetmax];   //[nGenVisTau]
   Float_t         GenVisTau_mass[njetmax];   //[nGenVisTau]
   Float_t         GenVisTau_phi[njetmax];   //[nGenVisTau]
   Float_t         GenVisTau_pt[njetmax];   //[nGenVisTau]
   Int_t           GenVisTau_charge[njetmax];   //[nGenVisTau]
   Int_t           GenVisTau_genPartIdxMother[njetmax];   //[nGenVisTau]
   Int_t           GenVisTau_status[njetmax];   //[nGenVisTau]
   Float_t         genWeight;
   Float_t         LHEWeight_originalXWGTUP;
   UInt_t          nLHEPdfWeight;
   Float_t         LHEPdfWeight[npdfmax];   //[nLHEPdfWeight]
   UInt_t          nLHEScaleWeight;
   Float_t         LHEScaleWeight[nscalemax];   //[nLHEScaleWeight]
   UInt_t          nPSWeight;
   Float_t         PSWeight[npsmax];   //[nPSWeight]
   UInt_t          nIsoTrack;
   Float_t         IsoTrack_dxy[njetmax];   //[nIsoTrack]
   Float_t         IsoTrack_dz[njetmax];   //[nIsoTrack]
   Float_t         IsoTrack_eta[njetmax];   //[nIsoTrack]
   Float_t         IsoTrack_pfRelIso03_all[njetmax];   //[nIsoTrack]
   Float_t         IsoTrack_pfRelIso03_chg[njetmax];   //[nIsoTrack]
   Float_t         IsoTrack_phi[njetmax];   //[nIsoTrack]
   Float_t         IsoTrack_pt[njetmax];   //[nIsoTrack]
   Float_t         IsoTrack_miniPFRelIso_all[njetmax];   //[nIsoTrack]
   Float_t         IsoTrack_miniPFRelIso_chg[njetmax];   //[nIsoTrack]
   Int_t           IsoTrack_fromPV[njetmax];   //[nIsoTrack]
   Int_t           IsoTrack_pdgId[njetmax];   //[nIsoTrack]
   Bool_t          IsoTrack_isHighPurityTrack[njetmax];   //[nIsoTrack]
   Bool_t          IsoTrack_isPFcand[njetmax];   //[nIsoTrack]
   Bool_t          IsoTrack_isFromLostTrack[njetmax];   //[nIsoTrack]
   UInt_t          nJet;
   Float_t         Jet_area[njetmax];   //[nJet]
   Float_t         Jet_btagCMVA[njetmax];   //[nJet]
   Float_t         Jet_btagCSVV2[njetmax];   //[nJet]
   Float_t         Jet_btagDeepB[njetmax];   //[nJet]
   Float_t         Jet_btagDeepC[njetmax];   //[nJet]
   Float_t         Jet_btagDeepFlavB[njetmax];   //[nJet]
   Float_t         Jet_chEmEF[njetmax];   //[nJet]
   Float_t         Jet_chHEF[njetmax];   //[nJet]
   Float_t         Jet_eta[njetmax];   //[nJet]
   Float_t         Jet_mass[njetmax];   //[nJet]
   Float_t         Jet_muEF[njetmax];   //[nJet]
   Float_t         Jet_neEmEF[njetmax];   //[nJet]
   Float_t         Jet_neHEF[njetmax];   //[nJet]
   Float_t         Jet_phi[njetmax];   //[nJet]
   Float_t         Jet_pt[njetmax];   //[nJet]
   Float_t         Jet_qgl[njetmax];   //[nJet]
   Float_t         Jet_rawFactor[njetmax];   //[nJet]
   Float_t         Jet_bRegCorr[njetmax];   //[nJet]
   Float_t         Jet_bRegRes[njetmax];   //[nJet]
   Int_t           Jet_electronIdx1[njetmax];   //[nJet]
   Int_t           Jet_electronIdx2[njetmax];   //[nJet]
   Int_t           Jet_jetId[njetmax];   //[nJet]
   Int_t 		   Jet_MatchFatJet[njetmax];
   Int_t           Jet_muonIdx1[njetmax];   //[nJet]
   Int_t           Jet_muonIdx2[njetmax];   //[nJet]
   Int_t           Jet_nConstituents[njetmax];   //[nJet]
   Int_t           Jet_nElectrons[njetmax];   //[nJet]
   Int_t           Jet_nMuons[njetmax];   //[nJet]
   Int_t           Jet_puId[njetmax];   //[nJet]
   Float_t 		   Jet_CSVv2MCeff[njetmax];
   Float_t 		   Jet_DeepCSVMCeff[njetmax];
   Float_t 		   Jet_pt_jes_up[njesmax][njetmax];  
   Float_t 		   Jet_mass_jes_up[njesmax][njetmax];  
   Float_t 		   Jet_pt_jer_up[njetmax];  
   Float_t 		   Jet_mass_jer_up[njetmax];  
   Float_t 		   Jet_mass_jms_up[njetmax]; 
   Float_t 		   Jet_mass_jmr_up[njetmax];  
   Float_t 		   Jet_pt_jes_dn[njesmax][njetmax];  
   Float_t 		   Jet_mass_jes_dn[njesmax][njetmax];  
   Float_t 		   Jet_pt_jer_dn[njetmax];  
   Float_t 		   Jet_mass_jer_dn[njetmax]; 
   Float_t 		   Jet_mass_jms_dn[njetmax]; 
   Float_t 		   Jet_mass_jmr_dn[njetmax]; 
   Float_t 		   Jet_deepflavbtagSF_unc_up[nbtagmax][njetmax]; 
   Float_t 		   Jet_deepflavbtagSF_unc_dn[nbtagmax][njetmax];  
   Float_t         LHE_HT;
   Float_t         LHE_HTIncoming;
   Float_t         LHE_Vpt;
   UChar_t         LHE_Njets;
   UChar_t         LHE_Nb;
   UChar_t         LHE_Nc;
   UChar_t         LHE_Nuds;
   UChar_t         LHE_Nglu;
   UChar_t         LHE_NpNLO;
   UChar_t         LHE_NpLO;
   UInt_t          nLHEPart;
   Float_t         LHEPart_pt[npartmx];   //[nLHEPart]
   Float_t         LHEPart_eta[npartmx];   //[nLHEPart]
   Float_t         LHEPart_phi[npartmx];   //[nLHEPart]
   Float_t         LHEPart_mass[npartmx];   //[nLHEPart]
   Int_t           LHEPart_pdgId[npartmx];   //[nLHEPart]
   Float_t         METFixEE2017_MetUnclustEnUpDeltaX;
   Float_t         METFixEE2017_MetUnclustEnUpDeltaY;
   Float_t         METFixEE2017_covXX;
   Float_t         METFixEE2017_covXY;
   Float_t         METFixEE2017_covYY;
   Float_t         METFixEE2017_phi;
   Float_t         METFixEE2017_pt;
   Float_t         METFixEE2017_significance;
   Float_t         METFixEE2017_sumEt;
   Float_t         GenMET_phi;
   Float_t         GenMET_pt;
   Float_t         MET_MetUnclustEnUpDeltaX;
   Float_t         MET_MetUnclustEnUpDeltaY;
   Float_t         MET_covXX;
   Float_t         MET_covXY;
   Float_t         MET_covYY;
   Float_t         MET_phi;
   Float_t         MET_pt;
   Float_t         MET_significance;
   Float_t         MET_sumEt;
   UInt_t          nMuon;
   Float_t         Muon_dxy[njetmax];   //[nMuon]
   Float_t         Muon_dxyErr[njetmax];   //[nMuon]
   Float_t         Muon_dz[njetmax];   //[nMuon]
   Float_t         Muon_dzErr[njetmax];   //[nMuon]
   Float_t         Muon_eta[njetmax];   //[nMuon]
   Float_t         Muon_ip3d[njetmax];   //[nMuon]
   Float_t         Muon_jetPtRelv2[njetmax];   //[nMuon]
   Float_t         Muon_jetRelIso[njetmax];   //[nMuon]
   Float_t         Muon_mass[njetmax];   //[nMuon]
   Float_t         Muon_miniPFRelIso_all[njetmax];   //[nMuon]
   Float_t         Muon_miniPFRelIso_chg[njetmax];   //[nMuon]
   Float_t         Muon_pfRelIso03_all[njetmax];   //[nMuon]
   Float_t         Muon_pfRelIso03_chg[njetmax];   //[nMuon]
   Float_t         Muon_pfRelIso04_all[njetmax];   //[nMuon]
   Float_t         Muon_phi[njetmax];   //[nMuon]
   Float_t         Muon_pt[njetmax];   //[nMuon]
   Float_t         Muon_ptErr[njetmax];   //[nMuon]
   Float_t         Muon_segmentComp[njetmax];   //[nMuon]
   Float_t         Muon_sip3d[njetmax];   //[nMuon]
   Float_t         Muon_mvaTTH[njetmax];   //[nMuon]
   Int_t           Muon_charge[njetmax];   //[nMuon]
   Int_t           Muon_jetIdx[njetmax];   //[nMuon]
   Int_t           Muon_nStations[njetmax];   //[nMuon]
   Int_t           Muon_nTrackerLayers[njetmax];   //[nMuon]
   Int_t           Muon_pdgId[njetmax];   //[nMuon]
   Int_t           Muon_tightCharge[njetmax];   //[nMuon]
   UChar_t         Muon_highPtId[njetmax];   //[nMuon]
   Bool_t          Muon_inTimeMuon[njetmax];   //[nMuon]
   Bool_t          Muon_isGlobal[njetmax];   //[nMuon]
   Bool_t          Muon_isPFcand[njetmax];   //[nMuon]
   Bool_t          Muon_isTracker[njetmax];   //[nMuon]
   Bool_t          Muon_mediumId[njetmax];   //[nMuon]
   Bool_t          Muon_mediumPromptId[njetmax];   //[nMuon]
   UChar_t         Muon_miniIsoId[njetmax];   //[nMuon]
   UChar_t         Muon_multiIsoId[njetmax];   //[nMuon]
   UChar_t         Muon_mvaId[njetmax];   //[nMuon]
   UChar_t         Muon_pfIsoId[njetmax];   //[nMuon]
   Bool_t          Muon_softId[njetmax];   //[nMuon]
   Bool_t          Muon_softMvaId[njetmax];   //[nMuon]
   Bool_t          Muon_tightId[njetmax];   //[nMuon]
   UChar_t         Muon_tkIsoId[njetmax];   //[nMuon]
   Bool_t          Muon_triggerIdLoose[njetmax];   //[nMuon]
   UInt_t          nPhoton;
   Float_t         Photon_energyErr[njetmax];   //[nPhoton]
   Float_t         Photon_eta[njetmax];   //[nPhoton]
   Float_t         Photon_hoe[njetmax];   //[nPhoton]
   Float_t         Photon_mass[njetmax];   //[nPhoton]
   Float_t         Photon_mvaID[njetmax];   //[nPhoton]
   Float_t         Photon_mvaIDV1[njetmax];   //[nPhoton]
   Float_t         Photon_pfRelIso03_all[njetmax];   //[nPhoton]
   Float_t         Photon_pfRelIso03_chg[njetmax];   //[nPhoton]
   Float_t         Photon_phi[njetmax];   //[nPhoton]
   Float_t         Photon_pt[njetmax];   //[nPhoton]
   Float_t         Photon_r9[njetmax];   //[nPhoton]
   Float_t         Photon_sieie[njetmax];   //[nPhoton]
   Float_t 		   Photon_eCorr[njetmax];   //[nPhoton]
   Int_t           Photon_charge[njetmax];   //[nPhoton]
   Int_t           Photon_cutBasedBitmap[njetmax];   //[nPhoton]
   Int_t           Photon_cutBasedV1Bitmap[njetmax];   //[nPhoton]
   Int_t           Photon_electronIdx[njetmax];   //[nPhoton]
   Int_t           Photon_jetIdx[njetmax];   //[nPhoton]
   Int_t           Photon_pdgId[njetmax];   //[nPhoton]
   Int_t           Photon_vidNestedWPBitmap[njetmax];   //[nPhoton]
   Bool_t          Photon_electronVeto[njetmax];   //[nPhoton]
   Bool_t          Photon_isScEtaEB[njetmax];   //[nPhoton]
   Bool_t          Photon_isScEtaEE[njetmax];   //[nPhoton]
   Bool_t          Photon_mvaID_WP80[njetmax];   //[nPhoton]
   Bool_t          Photon_mvaID_WP90[njetmax];   //[nPhoton]
   Bool_t          Photon_pixelSeed[njetmax];   //[nPhoton]
   Float_t         Pileup_nTrueInt;
   Int_t           Pileup_nPU;
   Int_t           Pileup_sumEOOT;
   Int_t           Pileup_sumLOOT;
   Float_t         PuppiMET_phi;
   Float_t         PuppiMET_pt;
   Float_t         PuppiMET_sumEt;
   Float_t         RawMET_phi;
   Float_t         RawMET_pt;
   Float_t         RawMET_sumEt;
   Float_t         fixedGridRhoFastjetAll;
   Float_t         fixedGridRhoFastjetCentralCalo;
   Float_t         fixedGridRhoFastjetCentralNeutral;
   UInt_t          nGenDressedLepton;
   Float_t         GenDressedLepton_eta[njetmax];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_mass[njetmax];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_phi[njetmax];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_pt[njetmax];   //[nGenDressedLepton]
   Int_t           GenDressedLepton_pdgId[njetmax];   //[nGenDressedLepton]
   UInt_t          nSoftActivityJet;
   Float_t         SoftActivityJet_eta[njetmax];   //[nSoftActivityJet]
   Float_t         SoftActivityJet_phi[njetmax];   //[nSoftActivityJet]
   Float_t         SoftActivityJet_pt[njetmax];   //[nSoftActivityJet]
   Float_t         SoftActivityJetHT;
   Float_t         SoftActivityJetHT10;
   Float_t         SoftActivityJetHT2;
   Float_t         SoftActivityJetHT5;
   Int_t           SoftActivityJetNjets10;
   Int_t           SoftActivityJetNjets2;
   Int_t           SoftActivityJetNjets5;
   UInt_t          nSubJet;
   Float_t         SubJet_btagCMVA[njetmax];   //[nSubJet]
   Float_t         SubJet_btagCSVV2[njetmax];   //[nSubJet]
   Float_t         SubJet_btagDeepB[njetmax];   //[nSubJet]
   Float_t         SubJet_eta[njetmax];   //[nSubJet]
   Float_t         SubJet_mass[njetmax];   //[nSubJet]
   Float_t         SubJet_n2b1[njetmax];   //[nSubJet]
   Float_t         SubJet_n3b1[njetmax];   //[nSubJet]
   Float_t         SubJet_phi[njetmax];   //[nSubJet]
   Float_t         SubJet_pt[njetmax];   //[nSubJet]
   Float_t         SubJet_rawFactor[njetmax];   //[nSubJet]
   Float_t         SubJet_tau1[njetmax];   //[nSubJet]
   Float_t         SubJet_tau2[njetmax];   //[nSubJet]
   Float_t         SubJet_tau3[njetmax];   //[nSubJet]
   Float_t         SubJet_tau4[njetmax];   //[nSubJet]
   UInt_t          nTau;
   Float_t         Tau_chargedIso[njetmax];   //[nTau]
   Float_t         Tau_dxy[njetmax];   //[nTau]
   Float_t         Tau_dz[njetmax];   //[nTau]
   Float_t         Tau_eta[njetmax];   //[nTau]
   Float_t         Tau_leadTkDeltaEta[njetmax];   //[nTau]
   Float_t         Tau_leadTkDeltaPhi[njetmax];   //[nTau]
   Float_t         Tau_leadTkPtOverTauPt[njetmax];   //[nTau]
   Float_t         Tau_mass[njetmax];   //[nTau]
   Float_t         Tau_neutralIso[njetmax];   //[nTau]
   Float_t         Tau_phi[njetmax];   //[nTau]
   Float_t         Tau_photonsOutsideSignalCone[njetmax];   //[nTau]
   Float_t         Tau_pt[njetmax];   //[nTau]
   Float_t         Tau_puCorr[njetmax];   //[nTau]
   Float_t         Tau_rawAntiEle[njetmax];   //[nTau]
   Float_t         Tau_rawIso[njetmax];   //[nTau]
   Float_t         Tau_rawIsodR03[njetmax];   //[nTau]
   Float_t         Tau_rawMVAnewDM2017v2[njetmax];   //[nTau]
   Float_t         Tau_rawMVAoldDM[njetmax];   //[nTau]
   Float_t         Tau_rawMVAoldDM2017v1[njetmax];   //[nTau]
   Float_t         Tau_rawMVAoldDM2017v2[njetmax];   //[nTau]
   Float_t         Tau_rawMVAoldDMdR032017v2[njetmax];   //[nTau]
   Int_t           Tau_charge[njetmax];   //[nTau]
   Int_t           Tau_decayMode[njetmax];   //[nTau]
   Int_t           Tau_jetIdx[njetmax];   //[nTau]
   Int_t           Tau_rawAntiEleCat[njetmax];   //[nTau]
   UChar_t         Tau_idAntiEle[njetmax];   //[nTau]
   UChar_t         Tau_idAntiMu[njetmax];   //[nTau]
   Bool_t          Tau_idDecayMode[njetmax];   //[nTau]
   Bool_t          Tau_idDecayModeNewDMs[njetmax];   //[nTau]
   UChar_t         Tau_idMVAnewDM2017v2[njetmax];   //[nTau]
   UChar_t         Tau_idMVAoldDM[njetmax];   //[nTau]
   UChar_t         Tau_idMVAoldDM2017v1[njetmax];   //[nTau]
   UChar_t         Tau_idMVAoldDM2017v2[njetmax];   //[nTau]
   UChar_t         Tau_idMVAoldDMdR032017v2[njetmax];   //[nTau]
   Float_t         TkMET_phi;
   Float_t         TkMET_pt;
   Float_t         TkMET_sumEt;
   UInt_t          nTrigObj;
   Float_t         TrigObj_pt[nTrigObjmax];   //[nTrigObj]
   Float_t         TrigObj_eta[nTrigObjmax];   //[nTrigObj]
   Float_t         TrigObj_phi[nTrigObjmax];   //[nTrigObj]
   Float_t         TrigObj_l1pt[nTrigObjmax];   //[nTrigObj]
   Float_t         TrigObj_l1pt_2[nTrigObjmax];   //[nTrigObj]
   Float_t         TrigObj_l2pt[nTrigObjmax];   //[nTrigObj]
   Int_t           TrigObj_id[nTrigObjmax];   //[nTrigObj]
   Int_t           TrigObj_l1iso[nTrigObjmax];   //[nTrigObj]
   Int_t           TrigObj_l1charge[nTrigObjmax];   //[nTrigObj]
   Int_t           TrigObj_filterBits[nTrigObjmax];   //[nTrigObj]
   Int_t           genTtbarId;
   UInt_t          nOtherPV;
   Float_t         OtherPV_z[nSVmax];   //[nOtherPV]
   Float_t         PV_ndof;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_chi2;
   Float_t         PV_score;
   Int_t           PV_npvs;
   Int_t           PV_npvsGood;
   UInt_t          nSV;
   Float_t         SV_dlen[nSVmax];   //[nSV]
   Float_t         SV_dlenSig[nSVmax];   //[nSV]
   Float_t         SV_pAngle[nSVmax];   //[nSV]
   Int_t           Electron_genPartIdx[njetmax];   //[nElectron]
   UChar_t         Electron_genPartFlav[njetmax];   //[nElectron]
   Int_t           GenJetAK8_partonFlavour[njetmax];   //[nGenJetAK8]
   UChar_t         GenJetAK8_hadronFlavour[njetmax];   //[nGenJetAK8]
   Int_t           GenJet_partonFlavour[njetmax];   //[nGenJet]
   UChar_t         GenJet_hadronFlavour[njetmax];   //[nGenJet]
   Int_t           Jet_genJetIdx[njetmax];   //[nJet]
   Int_t           Jet_hadronFlavour[njetmax];   //[nJet]
   Int_t           Jet_partonFlavour[njetmax];   //[nJet]
   Int_t           Muon_genPartIdx[njetmax];   //[nMuon]
   UChar_t         Muon_genPartFlav[njetmax];   //[nMuon]
   Int_t           Photon_genPartIdx[njetmax];   //[nPhoton]
   UChar_t         Photon_genPartFlav[njetmax];   //[nPhoton]
   Float_t         MET_fiducialGenPhi;
   Float_t         MET_fiducialGenPt;
   UChar_t         Electron_cleanmask[njetmax];   //[nElectron]
   UChar_t         Jet_cleanmask[njetmax];   //[nJet]
   UChar_t         Muon_cleanmask[njetmax];   //[nMuon]
   UChar_t         Photon_cleanmask[njetmax];   //[nPhoton]
   UChar_t         Tau_cleanmask[njetmax];   //[nTau]
   Float_t         SV_chi2[nSVmax];   //[nSV]
   Float_t         SV_eta[nSVmax];   //[nSV]
   Float_t         SV_mass[nSVmax];   //[nSV]
   Float_t         SV_ndof[nSVmax];   //[nSV]
   Float_t         SV_phi[nSVmax];   //[nSV]
   Float_t         SV_pt[nSVmax];   //[nSV]
   Float_t         SV_x[nSVmax];   //[nSV]
   Float_t         SV_y[nSVmax];   //[nSV]
   Float_t         SV_z[nSVmax];   //[nSV]
   Int_t           Tau_genPartIdx[njetmax];   //[nTau]
   UChar_t         Tau_genPartFlav[njetmax];   //[nTau]
   Bool_t          L1simulation_step;
   Bool_t          HLTriggerFirstPath;

   Bool_t          HLT_AK8PFJet360_TrimMass30;
   Bool_t          HLT_AK8PFJet380_TrimMass30;
   Bool_t          HLT_AK8PFJet400_TrimMass30;
   Bool_t          HLT_AK8PFJet420_TrimMass30;
   Bool_t		   HLT_AK8PFHT700_TrimR0p1PT0p03Mass50;
   Bool_t          HLT_AK8PFHT750_TrimMass50;
   Bool_t          HLT_AK8PFHT800_TrimMass50;
   Bool_t          HLT_AK8PFHT850_TrimMass50;
   Bool_t          HLT_AK8PFHT900_TrimMass50;
   
   Bool_t          HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20;
   Bool_t 		   HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20;
   Bool_t          HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087;
   Bool_t          HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087;
   Bool_t          HLT_AK8DiPFJet300_200_TrimMass30;
   Bool_t          HLT_AK8DiPFJet280_200_TrimMass30;
   
   Bool_t          HLT_CaloJet500_NoJetID;
   Bool_t          HLT_CaloJet550_NoJetID;
   Bool_t          HLT_Trimuon5_3p5_2_Upsilon_Muon;
   Bool_t          HLT_DoubleEle25_CaloIdL_MW;
   Bool_t          HLT_DoubleEle27_CaloIdL_MW;
   Bool_t          HLT_DoubleEle33_CaloIdL_MW;
   Bool_t          HLT_DoubleEle24_eta2p1_WPTight_Gsf;
   Bool_t          HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350;
   Bool_t          HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350;
   Bool_t          HLT_Ele27_Ele37_CaloIdL_MW;
   Bool_t          HLT_Mu27_Ele37_CaloIdL_MW;
   Bool_t          HLT_Mu37_Ele27_CaloIdL_MW;
   Bool_t          HLT_Mu37_TkMu27;
   Bool_t          HLT_DoubleMu4_3_Bs;
   Bool_t          HLT_DoubleMu4_3_Jpsi_Displaced;
   Bool_t          HLT_DoubleMu4_JpsiTrk_Displaced;
   Bool_t          HLT_DoubleMu4_LowMassNonResonantTrk_Displaced;
   Bool_t          HLT_DoubleMu3_Trk_Tau3mu;
   Bool_t          HLT_DoubleMu4_PsiPrimeTrk_Displaced;
   Bool_t          HLT_DoubleMu4_Mass8_DZ_PFHT350;
   Bool_t          HLT_DoubleMu8_Mass8_PFHT350;
   Bool_t          HLT_Mu3_PFJet40;
   Bool_t          HLT_Mu7p5_L2Mu2_Jpsi;
   Bool_t          HLT_Mu7p5_L2Mu2_Upsilon;
   Bool_t          HLT_Mu7p5_Track2_Jpsi;
   Bool_t          HLT_Mu7p5_Track3p5_Jpsi;
   Bool_t          HLT_Mu7p5_Track7_Jpsi;
   Bool_t          HLT_Mu7p5_Track2_Upsilon;
   Bool_t          HLT_Mu7p5_Track3p5_Upsilon;
   Bool_t          HLT_Mu7p5_Track7_Upsilon;
   Bool_t          HLT_DoublePhoton33_CaloIdL;
   Bool_t          HLT_DoublePhoton70;
   Bool_t          HLT_DoublePhoton85;
   Bool_t          HLT_Ele20_WPTight_Gsf;
   Bool_t          HLT_Ele20_WPLoose_Gsf;
   Bool_t          HLT_Ele20_eta2p1_WPLoose_Gsf;
   Bool_t          HLT_DiEle27_WPTightCaloOnly_L1DoubleEG;
   Bool_t          HLT_Ele27_WPTight_Gsf;
   Bool_t          HLT_Ele32_WPTight_Gsf;
   Bool_t          HLT_Ele35_WPTight_Gsf;
   Bool_t          HLT_Ele35_WPTight_Gsf_L1EGMT;
   Bool_t          HLT_Ele38_WPTight_Gsf;
   Bool_t          HLT_Ele40_WPTight_Gsf;
   Bool_t          HLT_Ele32_WPTight_Gsf_L1DoubleEG;
   Bool_t          HLT_HT450_Beamspot;
   Bool_t          HLT_HT300_Beamspot;
   Bool_t          HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1;
   Bool_t          HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1;
   Bool_t          HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1;
   Bool_t          HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1;
   Bool_t          HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1;
   Bool_t          HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1;
   Bool_t          HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1;
   Bool_t          HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1;
   Bool_t          HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1;
   Bool_t          HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_TightID_SingleL1;
   Bool_t          HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_TightID_SingleL1;
   Bool_t          HLT_IsoMu20;
   Bool_t          HLT_IsoMu24;
   Bool_t          HLT_IsoMu24_eta2p1;
   Bool_t          HLT_IsoMu27;
   Bool_t          HLT_IsoMu30;
   Bool_t          HLT_UncorrectedJetE30_NoBPTX;
   Bool_t          HLT_UncorrectedJetE30_NoBPTX3BX;
   Bool_t          HLT_UncorrectedJetE60_NoBPTX3BX;
   Bool_t          HLT_UncorrectedJetE70_NoBPTX3BX;
   Bool_t          HLT_L1SingleMu18;
   Bool_t          HLT_L1SingleMu25;
   Bool_t          HLT_L2Mu10;
   Bool_t          HLT_L2Mu10_NoVertex_NoBPTX3BX;
   Bool_t          HLT_L2Mu10_NoVertex_NoBPTX;
   Bool_t          HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX;
   Bool_t          HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX;
   Bool_t          HLT_L2Mu50;
   Bool_t          HLT_DoubleL2Mu50;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;
   Bool_t          HLT_Mu25_TkMu0_Onia;
   Bool_t          HLT_Mu30_TkMu0_Onia;
   Bool_t          HLT_Mu20_TkMu0_Phi;
   Bool_t          HLT_Mu25_TkMu0_Phi;
   Bool_t          HLT_Mu20;
   Bool_t          HLT_Mu27;
   Bool_t          HLT_Mu50;
   Bool_t          HLT_Mu55;
   Bool_t          HLT_OldMu100;
   Bool_t          HLT_TkMu100;
   Bool_t          HLT_DiPFJet15_NoCaloMatched;
   Bool_t          HLT_DiPFJet25_NoCaloMatched;
   Bool_t          HLT_DiPFJet15_FBEta3_NoCaloMatched;
   Bool_t          HLT_DiPFJet25_FBEta3_NoCaloMatched;
   Bool_t          HLT_DiPFJetAve40;
   Bool_t          HLT_DiPFJetAve60;
   Bool_t          HLT_DiPFJetAve80;
   Bool_t          HLT_DiPFJetAve140;
   Bool_t          HLT_DiPFJetAve200;
   Bool_t          HLT_DiPFJetAve260;
   Bool_t          HLT_DiPFJetAve320;
   Bool_t          HLT_DiPFJetAve400;
   Bool_t          HLT_DiPFJetAve500;
   Bool_t          HLT_DiPFJetAve15_HFJEC;
   Bool_t          HLT_DiPFJetAve25_HFJEC;
   Bool_t          HLT_DiPFJetAve35_HFJEC;
   Bool_t          HLT_DiPFJetAve60_HFJEC;
   Bool_t          HLT_DiPFJetAve80_HFJEC;
   Bool_t          HLT_DiPFJetAve100_HFJEC;
   Bool_t          HLT_DiPFJetAve160_HFJEC;
   Bool_t          HLT_DiPFJetAve220_HFJEC;
   Bool_t          HLT_DiPFJetAve300_HFJEC;
   Bool_t          HLT_AK8PFJet40;
   Bool_t          HLT_AK8PFJet60;
   Bool_t          HLT_AK8PFJet80;
   Bool_t          HLT_AK8PFJet140;
   Bool_t          HLT_AK8PFJet200;
   Bool_t          HLT_AK8PFJet260;
   Bool_t          HLT_AK8PFJet320;
   Bool_t          HLT_AK8PFJet400;
   Bool_t          HLT_AK8PFJet450;
   Bool_t          HLT_AK8PFJet500;
   Bool_t          HLT_AK8PFJet550;
   Bool_t          HLT_PFJet40;
   Bool_t          HLT_PFJet60;
   Bool_t          HLT_PFJet80;
   Bool_t          HLT_PFJet140;
   Bool_t          HLT_PFJet200;
   Bool_t          HLT_PFJet260;
   Bool_t          HLT_PFJet320;
   Bool_t          HLT_PFJet400;
   Bool_t          HLT_PFJet450;
   Bool_t          HLT_PFJet500;
   Bool_t          HLT_PFJet550;
   Bool_t          HLT_PFJetFwd40;
   Bool_t          HLT_PFJetFwd60;
   Bool_t          HLT_PFJetFwd80;
   Bool_t          HLT_PFJetFwd140;
   Bool_t          HLT_PFJetFwd200;
   Bool_t          HLT_PFJetFwd260;
   Bool_t          HLT_PFJetFwd320;
   Bool_t          HLT_PFJetFwd400;
   Bool_t          HLT_PFJetFwd450;
   Bool_t          HLT_PFJetFwd500;
   Bool_t          HLT_AK8PFJetFwd40;
   Bool_t          HLT_AK8PFJetFwd60;
   Bool_t          HLT_AK8PFJetFwd80;
   Bool_t          HLT_AK8PFJetFwd140;
   Bool_t          HLT_AK8PFJetFwd200;
   Bool_t          HLT_AK8PFJetFwd260;
   Bool_t          HLT_AK8PFJetFwd320;
   Bool_t          HLT_AK8PFJetFwd400;
   Bool_t          HLT_AK8PFJetFwd450;
   Bool_t          HLT_AK8PFJetFwd500;
   Bool_t          HLT_PFHT180;
   Bool_t          HLT_PFHT250;
   Bool_t          HLT_PFHT370;
   Bool_t          HLT_PFHT430;
   Bool_t          HLT_PFHT510;
   Bool_t          HLT_PFHT590;
   Bool_t          HLT_PFHT680;
   Bool_t          HLT_PFHT780;
   Bool_t          HLT_PFHT890;
   Bool_t          HLT_PFHT800;
   Bool_t          HLT_PFHT900;
   Bool_t          HLT_PFHT1050;
   Bool_t          HLT_PFHT500_PFMET100_PFMHT100_IDTight;
   Bool_t          HLT_PFHT500_PFMET110_PFMHT110_IDTight;
   Bool_t          HLT_PFHT700_PFMET85_PFMHT85_IDTight;
   Bool_t          HLT_PFHT700_PFMET95_PFMHT95_IDTight;
   Bool_t          HLT_PFHT800_PFMET75_PFMHT75_IDTight;
   Bool_t          HLT_PFHT800_PFMET85_PFMHT85_IDTight;
   Bool_t          HLT_PFMET110_PFMHT110_IDTight;
   Bool_t          HLT_PFMET120_PFMHT120_IDTight;
   Bool_t          HLT_PFMET130_PFMHT130_IDTight;
   Bool_t          HLT_PFMET140_PFMHT140_IDTight;
   Bool_t          HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1;
   Bool_t          HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1;
   Bool_t          HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1;
   Bool_t          HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1;
   Bool_t          HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1;
   Bool_t          HLT_PFMET120_PFMHT120_IDTight_PFHT60;
   Bool_t          HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
   Bool_t          HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60;
   Bool_t          HLT_PFMETTypeOne110_PFMHT110_IDTight;
   Bool_t          HLT_PFMETTypeOne120_PFMHT120_IDTight;
   Bool_t          HLT_PFMETTypeOne130_PFMHT130_IDTight;
   Bool_t          HLT_PFMETTypeOne140_PFMHT140_IDTight;
   Bool_t          HLT_PFMETNoMu110_PFMHTNoMu110_IDTight;
   Bool_t          HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;
   Bool_t          HLT_PFMETNoMu130_PFMHTNoMu130_IDTight;
   Bool_t          HLT_PFMETNoMu140_PFMHTNoMu140_IDTight;
   Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight;
   Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight;
   Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight;
   Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight;
   Bool_t          HLT_L1ETMHadSeeds;
   Bool_t          HLT_CaloMHT90;
   Bool_t          HLT_CaloMET80_NotCleaned;
   Bool_t          HLT_CaloMET90_NotCleaned;
   Bool_t          HLT_CaloMET100_NotCleaned;
   Bool_t          HLT_CaloMET110_NotCleaned;
   Bool_t          HLT_CaloMET250_NotCleaned;
   Bool_t          HLT_CaloMET70_HBHECleaned;
   Bool_t          HLT_CaloMET80_HBHECleaned;
   Bool_t          HLT_CaloMET90_HBHECleaned;
   Bool_t          HLT_CaloMET100_HBHECleaned;
   Bool_t          HLT_CaloMET250_HBHECleaned;
   Bool_t          HLT_CaloMET300_HBHECleaned;
   Bool_t          HLT_CaloMET350_HBHECleaned;
   Bool_t          HLT_PFMET200_NotCleaned;
   Bool_t          HLT_PFMET200_HBHECleaned;
   Bool_t          HLT_PFMET250_HBHECleaned;
   Bool_t          HLT_PFMET300_HBHECleaned;
   Bool_t          HLT_PFMET200_HBHE_BeamHaloCleaned;
   Bool_t          HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned;
   Bool_t          HLT_MET105_IsoTrk50;
   Bool_t          HLT_MET120_IsoTrk50;
   Bool_t          HLT_SingleJet30_Mu12_SinglePFJet40;
   Bool_t          HLT_Mu12_DoublePFJets40_CaloBTagCSV_p33;
   Bool_t          HLT_Mu12_DoublePFJets100_CaloBTagCSV_p33;
   Bool_t          HLT_Mu12_DoublePFJets200_CaloBTagCSV_p33;
   Bool_t          HLT_Mu12_DoublePFJets350_CaloBTagCSV_p33;
   Bool_t          HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33;
   Bool_t          HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagCSV_p33;
   Bool_t          HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagCSV_p33;
   Bool_t          HLT_DoublePFJets40_CaloBTagCSV_p33;
   Bool_t          HLT_DoublePFJets100_CaloBTagCSV_p33;
   Bool_t          HLT_DoublePFJets200_CaloBTagCSV_p33;
   Bool_t          HLT_DoublePFJets350_CaloBTagCSV_p33;
   Bool_t          HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33;
   Bool_t          HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33;
   Bool_t          HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33;
   Bool_t          HLT_Photon300_NoHE;
   Bool_t          HLT_Mu8_TrkIsoVVL;
   Bool_t          HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ;
   Bool_t          HLT_Mu8_DiEle12_CaloIdL_TrackIdL;
   Bool_t          HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ;
   Bool_t          HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu17_TrkIsoVVL;
   Bool_t          HLT_Mu19_TrkIsoVVL;
   Bool_t          HLT_BTagMu_AK4DiJet20_Mu5;
   Bool_t          HLT_BTagMu_AK4DiJet40_Mu5;
   Bool_t          HLT_BTagMu_AK4DiJet70_Mu5;
   Bool_t          HLT_BTagMu_AK4DiJet110_Mu5;
   Bool_t          HLT_BTagMu_AK4DiJet170_Mu5;
   Bool_t          HLT_BTagMu_AK4Jet300_Mu5;
   Bool_t          HLT_BTagMu_AK8DiJet170_Mu5;
   Bool_t          HLT_BTagMu_AK8Jet300_Mu5;
   Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Mu12_DoublePhoton20;
   Bool_t          HLT_TriplePhoton_20_20_20_CaloIdLV2;
   Bool_t          HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL;
   Bool_t          HLT_TriplePhoton_30_30_10_CaloIdLV2;
   Bool_t          HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL;
   Bool_t          HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL;
   Bool_t          HLT_Photon25;
   Bool_t          HLT_Photon33;
   Bool_t          HLT_Photon50;
   Bool_t          HLT_Photon75;
   Bool_t          HLT_Photon90;
   Bool_t          HLT_Photon120;
   Bool_t          HLT_Photon150;
   Bool_t          HLT_Photon175;
   Bool_t          HLT_Photon200;
   Bool_t          HLT_Photon50_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon75_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon90_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon120_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon165_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon90_CaloIdL_PFHT700;
   Bool_t          HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90;
   Bool_t          HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95;
   Bool_t          HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55;
   Bool_t          HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55;
   Bool_t          HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55;
   Bool_t          HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55;
   Bool_t          HLT_Dimuon0_Jpsi_L1_NoOS;
   Bool_t          HLT_Dimuon0_Jpsi_NoVertexing_NoOS;
   Bool_t          HLT_Dimuon0_Jpsi;
   Bool_t          HLT_Dimuon0_Jpsi_NoVertexing;
   Bool_t          HLT_Dimuon0_Jpsi_L1_4R_0er1p5R;
   Bool_t          HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R;
   Bool_t          HLT_Dimuon0_Jpsi3p5_Muon2;
   Bool_t          HLT_Dimuon0_Upsilon_L1_4p5;
   Bool_t          HLT_Dimuon0_Upsilon_L1_5;
   Bool_t          HLT_Dimuon0_Upsilon_L1_4p5NoOS;
   Bool_t          HLT_Dimuon0_Upsilon_L1_4p5er2p0;
   Bool_t          HLT_Dimuon0_Upsilon_L1_4p5er2p0M;
   Bool_t          HLT_Dimuon0_Upsilon_NoVertexing;
   Bool_t          HLT_Dimuon0_Upsilon_L1_5M;
   Bool_t          HLT_Dimuon0_LowMass_L1_0er1p5R;
   Bool_t          HLT_Dimuon0_LowMass_L1_0er1p5;
   Bool_t          HLT_Dimuon0_LowMass;
   Bool_t          HLT_Dimuon0_LowMass_L1_4;
   Bool_t          HLT_Dimuon0_LowMass_L1_4R;
   Bool_t          HLT_Dimuon0_LowMass_L1_TM530;
   Bool_t          HLT_Dimuon0_Upsilon_Muon_L1_TM0;
   Bool_t          HLT_Dimuon0_Upsilon_Muon_NoL1Mass;
   Bool_t          HLT_TripleMu_5_3_3_Mass3p8to60_DZ;
   Bool_t          HLT_TripleMu_10_5_5_DZ;
   Bool_t          HLT_TripleMu_12_10_5;
   Bool_t          HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15;
   Bool_t          HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1;
   Bool_t          HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15;
   Bool_t          HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1;
   Bool_t          HLT_DoubleMu3_DZ_PFMET50_PFMHT60;
   Bool_t          HLT_DoubleMu3_DZ_PFMET70_PFMHT70;
   Bool_t          HLT_DoubleMu3_DZ_PFMET90_PFMHT90;
   Bool_t          HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass;
   Bool_t          HLT_DoubleMu4_Jpsi_Displaced;
   Bool_t          HLT_DoubleMu4_Jpsi_NoVertexing;
   Bool_t          HLT_DoubleMu4_JpsiTrkTrk_Displaced;
   Bool_t          HLT_DoubleMu43NoFiltersNoVtx;
   Bool_t          HLT_DoubleMu48NoFiltersNoVtx;
   Bool_t          HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL;
   Bool_t          HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL;
   Bool_t          HLT_DoubleMu20_7_Mass0to30_L1_DM4;
   Bool_t          HLT_DoubleMu20_7_Mass0to30_L1_DM4EG;
   Bool_t          HLT_HT425;
   Bool_t          HLT_HT430_DisplacedDijet40_DisplacedTrack;
   Bool_t          HLT_HT430_DisplacedDijet60_DisplacedTrack;
   Bool_t          HLT_HT430_DisplacedDijet80_DisplacedTrack;
   Bool_t          HLT_HT400_DisplacedDijet40_DisplacedTrack;
   Bool_t          HLT_HT650_DisplacedDijet60_Inclusive;
   Bool_t          HLT_HT550_DisplacedDijet80_Inclusive;
   Bool_t          HLT_HT550_DisplacedDijet60_Inclusive;
   Bool_t          HLT_HT650_DisplacedDijet80_Inclusive;
   Bool_t          HLT_HT750_DisplacedDijet80_Inclusive;
   Bool_t          HLT_DiJet110_35_Mjj650_PFMET110;
   Bool_t          HLT_DiJet110_35_Mjj650_PFMET120;
   Bool_t          HLT_DiJet110_35_Mjj650_PFMET130;
   Bool_t          HLT_TripleJet110_35_35_Mjj650_PFMET110;
   Bool_t          HLT_TripleJet110_35_35_Mjj650_PFMET120;
   Bool_t          HLT_TripleJet110_35_35_Mjj650_PFMET130;
   Bool_t          HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg;
   Bool_t          HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg;
   Bool_t          HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg;
   Bool_t          HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned;
   Bool_t          HLT_Ele28_eta2p1_WPTight_Gsf_HT150;
   Bool_t          HLT_Ele28_HighEta_SC20_Mass55;
   Bool_t          HLT_DoubleMu20_7_Mass0to30_Photon23;
   Bool_t          HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5;
   Bool_t          HLT_Ele15_IsoVVVL_PFHT450_PFMET50;
   Bool_t          HLT_Ele15_IsoVVVL_PFHT450;
   Bool_t          HLT_Ele50_IsoVVVL_PFHT450;
   Bool_t          HLT_Ele15_IsoVVVL_PFHT600;
   Bool_t          HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60;
   Bool_t          HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60;
   Bool_t          HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5;
   Bool_t          HLT_Mu15_IsoVVVL_PFHT450_PFMET50;
   Bool_t          HLT_Mu15_IsoVVVL_PFHT450;
   Bool_t          HLT_Mu50_IsoVVVL_PFHT450;
   Bool_t          HLT_Mu15_IsoVVVL_PFHT600;
   Bool_t          HLT_Dimuon10_PsiPrime_Barrel_Seagulls;
   Bool_t          HLT_Dimuon20_Jpsi_Barrel_Seagulls;
   Bool_t          HLT_Dimuon10_Upsilon_Barrel_Seagulls;
   Bool_t          HLT_Dimuon12_Upsilon_eta1p5;
   Bool_t          HLT_Dimuon14_Phi_Barrel_Seagulls;
   Bool_t          HLT_Dimuon18_PsiPrime;
   Bool_t          HLT_Dimuon25_Jpsi;
   Bool_t          HLT_Dimuon18_PsiPrime_noCorrL1;
   Bool_t          HLT_Dimuon24_Upsilon_noCorrL1;
   Bool_t          HLT_Dimuon24_Phi_noCorrL1;
   Bool_t          HLT_Dimuon25_Jpsi_noCorrL1;
   Bool_t          HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;
   Bool_t          HLT_DiMu9_Ele9_CaloIdL_TrackIdL;
   Bool_t          HLT_DoubleIsoMu20_eta2p1;
   Bool_t          HLT_DoubleIsoMu24_eta2p1;
   Bool_t          HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx;
   Bool_t          HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx;
   Bool_t          HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx;
   Bool_t          HLT_Mu8;
   Bool_t          HLT_Mu17;
   Bool_t          HLT_Mu19;
   Bool_t          HLT_Mu17_Photon30_IsoCaloId;
   Bool_t          HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          HLT_Ele8_CaloIdM_TrackIdM_PFJet30;
   Bool_t          HLT_Ele17_CaloIdM_TrackIdM_PFJet30;
   Bool_t          HLT_Ele23_CaloIdM_TrackIdM_PFJet30;
   Bool_t          HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;
   Bool_t          HLT_Ele115_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_Ele135_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_Ele145_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_Ele200_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_Ele250_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_Ele300_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_PFHT300PT30_QuadPFJet_75_60_45_40;
   Bool_t          HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0;
   Bool_t          HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2;
   Bool_t          HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2;
   Bool_t          HLT_PFHT380_SixPFJet32;
   Bool_t          HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5;
   Bool_t          HLT_PFHT430_SixPFJet40;
   Bool_t          HLT_PFHT350;
   Bool_t          HLT_PFHT350MinPFJet15;
   Bool_t          HLT_Photon60_R9Id90_CaloIdL_IsoL;
   Bool_t          HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL;
   Bool_t          HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15;
   Bool_t          HLT_FullTrack_Multiplicity85;
   Bool_t          HLT_FullTrack_Multiplicity100;
   Bool_t          HLT_FullTrack_Multiplicity130;
   Bool_t          HLT_FullTrack_Multiplicity155;
   Bool_t          HLT_ECALHT800;
   Bool_t          HLT_DiSC30_18_EIso_AND_HE_Mass70;
   Bool_t          HLT_Physics;
   Bool_t          HLT_Physics_part0;
   Bool_t          HLT_Physics_part1;
   Bool_t          HLT_Physics_part2;
   Bool_t          HLT_Physics_part3;
   Bool_t          HLT_Physics_part4;
   Bool_t          HLT_Physics_part5;
   Bool_t          HLT_Physics_part6;
   Bool_t          HLT_Physics_part7;
   Bool_t          HLT_Random;
   Bool_t          HLT_ZeroBias;
   Bool_t          HLT_ZeroBias_part0;
   Bool_t          HLT_ZeroBias_part1;
   Bool_t          HLT_ZeroBias_part2;
   Bool_t          HLT_ZeroBias_part3;
   Bool_t          HLT_ZeroBias_part4;
   Bool_t          HLT_ZeroBias_part5;
   Bool_t          HLT_ZeroBias_part6;
   Bool_t          HLT_ZeroBias_part7;
   Bool_t          HLT_AK4CaloJet30;
   Bool_t          HLT_AK4CaloJet40;
   Bool_t          HLT_AK4CaloJet50;
   Bool_t          HLT_AK4CaloJet80;
   Bool_t          HLT_AK4CaloJet100;
   Bool_t          HLT_AK4CaloJet120;
   Bool_t          HLT_AK4PFJet30;
   Bool_t          HLT_AK4PFJet50;
   Bool_t          HLT_AK4PFJet80;
   Bool_t          HLT_AK4PFJet100;
   Bool_t          HLT_AK4PFJet120;
   Bool_t          HLT_HISinglePhoton10_Eta3p1ForPPRef;
   Bool_t          HLT_HISinglePhoton20_Eta3p1ForPPRef;
   Bool_t          HLT_HISinglePhoton30_Eta3p1ForPPRef;
   Bool_t          HLT_HISinglePhoton40_Eta3p1ForPPRef;
   Bool_t          HLT_HISinglePhoton50_Eta3p1ForPPRef;
   Bool_t          HLT_HISinglePhoton60_Eta3p1ForPPRef;
   Bool_t          HLT_Photon20_HoverELoose;
   Bool_t          HLT_Photon30_HoverELoose;
   Bool_t          HLT_Photon40_HoverELoose;
   Bool_t          HLT_Photon50_HoverELoose;
   Bool_t          HLT_Photon60_HoverELoose;
   Bool_t          HLT_EcalCalibration;
   Bool_t          HLT_HcalCalibration;
   Bool_t          HLT_L1UnpairedBunchBptxMinus;
   Bool_t          HLT_L1UnpairedBunchBptxPlus;
   Bool_t          HLT_L1NotBptxOR;
   Bool_t          HLT_L1MinimumBiasHF_OR;
   Bool_t          HLT_L1MinimumBiasHF0OR;
   Bool_t          HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;
   Bool_t          HLT_HcalNZS;
   Bool_t          HLT_HcalPhiSym;
   Bool_t          HLT_HcalIsolatedbunch;
   Bool_t          HLT_IsoTrackHB;
   Bool_t          HLT_IsoTrackHE;
   Bool_t          HLT_ZeroBias_FirstCollisionAfterAbortGap;
   Bool_t          HLT_ZeroBias_IsolatedBunches;
   Bool_t          HLT_ZeroBias_FirstCollisionInTrain;
   Bool_t          HLT_ZeroBias_LastCollisionInTrain;
   Bool_t          HLT_ZeroBias_FirstBXAfterTrain;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1;
   Bool_t          HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg;
   Bool_t          HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg;
   Bool_t          HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg;
   Bool_t          HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg;
   Bool_t          HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg;
   Bool_t          HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg;
   Bool_t          HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg;
   Bool_t          HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg;
   Bool_t          HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg;
   Bool_t          HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg;
   Bool_t          HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg;
   Bool_t          HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg;
   Bool_t          HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr;
   Bool_t          HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90;
   Bool_t          HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100;
   Bool_t          HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110;
   Bool_t          HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120;
   Bool_t          HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130;
   Bool_t          HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr;
   Bool_t          HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr;
   Bool_t          HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1;
   Bool_t          HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1;
   Bool_t          HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;
   Bool_t          HLT_Rsq0p35;
   Bool_t          HLT_Rsq0p40;
   Bool_t          HLT_RsqMR300_Rsq0p09_MR200;
   Bool_t          HLT_RsqMR320_Rsq0p09_MR200;
   Bool_t          HLT_RsqMR300_Rsq0p09_MR200_4jet;
   Bool_t          HLT_RsqMR320_Rsq0p09_MR200_4jet;
   Bool_t          HLT_L1_DoubleJet30_Mass_Min400_Mu10;
   Bool_t          HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1;
   Bool_t          HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1;
   Bool_t          HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1;
   Bool_t          HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50;
   Bool_t          HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3;
   Bool_t          HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3;
   Bool_t          HLT_PFMET100_PFMHT100_IDTight_PFHT60;
   Bool_t          HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60;
   Bool_t          HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60;
   Bool_t          HLT_Mu18_Mu9_SameSign;
   Bool_t          HLT_Mu18_Mu9_SameSign_DZ;
   Bool_t          HLT_Mu18_Mu9;
   Bool_t          HLT_Mu18_Mu9_DZ;
   Bool_t          HLT_Mu20_Mu10_SameSign;
   Bool_t          HLT_Mu20_Mu10_SameSign_DZ;
   Bool_t          HLT_Mu20_Mu10;
   Bool_t          HLT_Mu20_Mu10_DZ;
   Bool_t          HLT_Mu23_Mu12_SameSign;
   Bool_t          HLT_Mu23_Mu12_SameSign_DZ;
   Bool_t          HLT_Mu23_Mu12;
   Bool_t          HLT_Mu23_Mu12_DZ;
   Bool_t          HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi;
   Bool_t          HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi;
   Bool_t          HLT_DoubleMu3_DCA_PFMET50_PFMHT60;
   Bool_t          HLT_TripleMu_5_3_3_Mass3p8to60_DCA;
   Bool_t          HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1;
   Bool_t          HLT_QuadPFJet103_88_75_15_DoubleBTagCSV_p013_p08_VBF1;
   Bool_t          HLT_QuadPFJet105_90_76_15_DoubleBTagCSV_p013_p08_VBF1;
   Bool_t          HLT_QuadPFJet111_90_80_15_DoubleBTagCSV_p013_p08_VBF1;
   Bool_t          HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2;
   Bool_t          HLT_QuadPFJet103_88_75_15_BTagCSV_p013_VBF2;
   Bool_t          HLT_QuadPFJet105_88_76_15_BTagCSV_p013_VBF2;
   Bool_t          HLT_QuadPFJet111_90_80_15_BTagCSV_p013_VBF2;
   Bool_t          HLT_QuadPFJet98_83_71_15;
   Bool_t          HLT_QuadPFJet103_88_75_15;
   Bool_t          HLT_QuadPFJet105_88_76_15;
   Bool_t          HLT_QuadPFJet111_90_80_15;
   Bool_t          HLT_AK8PFJet330_PFAK8BTagCSV_p17;
   Bool_t          HLT_AK8PFJet330_PFAK8BTagCSV_p1;
   Bool_t          HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55;
   Bool_t          HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55;
   Bool_t          HLTriggerFinalPath;
   Bool_t          Flag_HBHENoiseFilter;
   Bool_t          Flag_HBHENoiseIsoFilter;
   Bool_t          Flag_CSCTightHaloFilter;
   Bool_t          Flag_CSCTightHaloTrkMuUnvetoFilter;
   Bool_t          Flag_CSCTightHalo2015Filter;
   Bool_t          Flag_globalTightHalo2016Filter;
   Bool_t          Flag_globalSuperTightHalo2016Filter;
   Bool_t          Flag_HcalStripHaloFilter;
   Bool_t          Flag_hcalLaserEventFilter;
   Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          Flag_EcalDeadCellBoundaryEnergyFilter;
   Bool_t          Flag_ecalBadCalibFilter;
   Bool_t          Flag_goodVertices;
   Bool_t          Flag_eeBadScFilter;
   Bool_t          Flag_ecalLaserCorrFilter;
   Bool_t          Flag_trkPOGFilters;
   Bool_t          Flag_chargedHadronTrackResolutionFilter;
   Bool_t          Flag_muonBadTrackFilter;
   Bool_t          Flag_BadChargedCandidateFilter;
   Bool_t          Flag_BadPFMuonFilter;
   Bool_t          Flag_BadChargedCandidateSummer16Filter;
   Bool_t          Flag_BadPFMuonSummer16Filter;
   Bool_t          Flag_trkPOG_manystripclus53X;
   Bool_t          Flag_trkPOG_toomanystripclus53X;
   Bool_t          Flag_trkPOG_logErrorTooManyClusters;
   Bool_t          Flag_METFilters;
   Bool_t          L1Reco_step;
   Float_t         Jet_pt_raw[njetmax];   //[nJet]
   Float_t         Jet_pt_nom[njetmax];   //[nJet]
   Float_t         Jet_corr_JEC[njetmax];   //[nJet]
   Float_t         Jet_corr_JER[njetmax];   //[nJet]
   Float_t         Jet_mass_raw[njetmax];   //[nJet]
   Float_t         Jet_mass_nom[njetmax];   //[nJet]
   Float_t         MET_pt_nom;
   Float_t         MET_phi_nom;
   Float_t         Jet_pt_jerUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jerUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jmrUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jmsUp[njetmax];   //[nJet]
   Float_t         MET_pt_jerUp;
   Float_t         MET_phi_jerUp;
   Float_t         Jet_pt_jesAbsoluteStatUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteStatUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesAbsoluteStatUp;
   Float_t         MET_phi_jesAbsoluteStatUp;
   Float_t         Jet_pt_jesAbsoluteScaleUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteScaleUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesAbsoluteScaleUp;
   Float_t         MET_phi_jesAbsoluteScaleUp;
   Float_t         Jet_pt_jesAbsoluteFlavMapUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteFlavMapUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesAbsoluteFlavMapUp;
   Float_t         MET_phi_jesAbsoluteFlavMapUp;
   Float_t         Jet_pt_jesAbsoluteMPFBiasUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteMPFBiasUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesAbsoluteMPFBiasUp;
   Float_t         MET_phi_jesAbsoluteMPFBiasUp;
   Float_t         Jet_pt_jesFragmentationUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesFragmentationUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesFragmentationUp;
   Float_t         MET_phi_jesFragmentationUp;
   Float_t         Jet_pt_jesSinglePionECALUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesSinglePionECALUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesSinglePionECALUp;
   Float_t         MET_phi_jesSinglePionECALUp;
   Float_t         Jet_pt_jesSinglePionHCALUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesSinglePionHCALUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesSinglePionHCALUp;
   Float_t         MET_phi_jesSinglePionHCALUp;
   Float_t         Jet_pt_jesFlavorQCDUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesFlavorQCDUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesFlavorQCDUp;
   Float_t         MET_phi_jesFlavorQCDUp;
   Float_t         Jet_pt_jesTimePtEtaUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesTimePtEtaUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesTimePtEtaUp;
   Float_t         MET_phi_jesTimePtEtaUp;
   Float_t         Jet_pt_jesRelativeJEREC1Up[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativeJEREC1Up[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativeJEREC1Up;
   Float_t         MET_phi_jesRelativeJEREC1Up;
   Float_t         Jet_pt_jesRelativeJEREC2Up[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativeJEREC2Up[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativeJEREC2Up;
   Float_t         MET_phi_jesRelativeJEREC2Up;
   Float_t         Jet_pt_jesRelativeJERHFUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativeJERHFUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativeJERHFUp;
   Float_t         MET_phi_jesRelativeJERHFUp;
   Float_t         Jet_pt_jesRelativePtBBUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativePtBBUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativePtBBUp;
   Float_t         MET_phi_jesRelativePtBBUp;
   Float_t         Jet_pt_jesRelativePtEC1Up[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativePtEC1Up[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativePtEC1Up;
   Float_t         MET_phi_jesRelativePtEC1Up;
   Float_t         Jet_pt_jesRelativePtEC2Up[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativePtEC2Up[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativePtEC2Up;
   Float_t         MET_phi_jesRelativePtEC2Up;
   Float_t         Jet_pt_jesRelativePtHFUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativePtHFUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativePtHFUp;
   Float_t         MET_phi_jesRelativePtHFUp;
   Float_t         Jet_pt_jesRelativeBalUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativeBalUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativeBalUp;
   Float_t         MET_phi_jesRelativeBalUp;
   Float_t         Jet_pt_jesRelativeSampleUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativeSampleUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativeSampleUp;
   Float_t         MET_phi_jesRelativeSampleUp;
   Float_t         Jet_pt_jesRelativeFSRUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativeFSRUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativeFSRUp;
   Float_t         MET_phi_jesRelativeFSRUp;
   Float_t         Jet_pt_jesRelativeStatFSRUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativeStatFSRUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativeStatFSRUp;
   Float_t         MET_phi_jesRelativeStatFSRUp;
   Float_t         Jet_pt_jesRelativeStatECUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativeStatECUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativeStatECUp;
   Float_t         MET_phi_jesRelativeStatECUp;
   Float_t         Jet_pt_jesRelativeStatHFUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativeStatHFUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativeStatHFUp;
   Float_t         MET_phi_jesRelativeStatHFUp;
   Float_t         Jet_pt_jesPileUpDataMCUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesPileUpDataMCUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesPileUpDataMCUp;
   Float_t         MET_phi_jesPileUpDataMCUp;
   Float_t         Jet_pt_jesPileUpPtRefUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtRefUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesPileUpPtRefUp;
   Float_t         MET_phi_jesPileUpPtRefUp;
   Float_t         Jet_pt_jesPileUpPtBBUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtBBUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesPileUpPtBBUp;
   Float_t         MET_phi_jesPileUpPtBBUp;
   Float_t         Jet_pt_jesPileUpPtEC1Up[njetmax];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtEC1Up[njetmax];   //[nJet]
   Float_t         MET_pt_jesPileUpPtEC1Up;
   Float_t         MET_phi_jesPileUpPtEC1Up;
   Float_t         Jet_pt_jesPileUpPtEC2Up[njetmax];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtEC2Up[njetmax];   //[nJet]
   Float_t         MET_pt_jesPileUpPtEC2Up;
   Float_t         MET_phi_jesPileUpPtEC2Up;
   Float_t         Jet_pt_jesPileUpPtHFUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtHFUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesPileUpPtHFUp;
   Float_t         MET_phi_jesPileUpPtHFUp;
   Float_t         Jet_pt_jesPileUpMuZeroUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesPileUpMuZeroUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesPileUpMuZeroUp;
   Float_t         MET_phi_jesPileUpMuZeroUp;
   Float_t         Jet_pt_jesPileUpEnvelopeUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesPileUpEnvelopeUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesPileUpEnvelopeUp;
   Float_t         MET_phi_jesPileUpEnvelopeUp;
   Float_t         Jet_pt_jesSubTotalPileUpUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesSubTotalPileUpUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesSubTotalPileUpUp;
   Float_t         MET_phi_jesSubTotalPileUpUp;
   Float_t         Jet_pt_jesSubTotalRelativeUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesSubTotalRelativeUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesSubTotalRelativeUp;
   Float_t         MET_phi_jesSubTotalRelativeUp;
   Float_t         Jet_pt_jesSubTotalPtUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesSubTotalPtUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesSubTotalPtUp;
   Float_t         MET_phi_jesSubTotalPtUp;
   Float_t         Jet_pt_jesSubTotalScaleUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesSubTotalScaleUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesSubTotalScaleUp;
   Float_t         MET_phi_jesSubTotalScaleUp;
   Float_t         Jet_pt_jesSubTotalAbsoluteUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesSubTotalAbsoluteUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesSubTotalAbsoluteUp;
   Float_t         MET_phi_jesSubTotalAbsoluteUp;
   Float_t         Jet_pt_jesSubTotalMCUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesSubTotalMCUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesSubTotalMCUp;
   Float_t         MET_phi_jesSubTotalMCUp;
   Float_t         Jet_pt_jesTotalUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesTotalUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesTotalUp;
   Float_t         MET_phi_jesTotalUp;
   Float_t         Jet_pt_jesTotalNoFlavorUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesTotalNoFlavorUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesTotalNoFlavorUp;
   Float_t         MET_phi_jesTotalNoFlavorUp;
   Float_t         Jet_pt_jesTotalNoTimeUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesTotalNoTimeUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesTotalNoTimeUp;
   Float_t         MET_phi_jesTotalNoTimeUp;
   Float_t         Jet_pt_jesTotalNoFlavorNoTimeUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesTotalNoFlavorNoTimeUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesTotalNoFlavorNoTimeUp;
   Float_t         MET_phi_jesTotalNoFlavorNoTimeUp;
   Float_t         Jet_pt_jesFlavorZJetUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesFlavorZJetUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesFlavorZJetUp;
   Float_t         MET_phi_jesFlavorZJetUp;
   Float_t         Jet_pt_jesFlavorPhotonJetUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesFlavorPhotonJetUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesFlavorPhotonJetUp;
   Float_t         MET_phi_jesFlavorPhotonJetUp;
   Float_t         Jet_pt_jesFlavorPureGluonUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureGluonUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesFlavorPureGluonUp;
   Float_t         MET_phi_jesFlavorPureGluonUp;
   Float_t         Jet_pt_jesFlavorPureQuarkUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureQuarkUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesFlavorPureQuarkUp;
   Float_t         MET_phi_jesFlavorPureQuarkUp;
   Float_t         Jet_pt_jesFlavorPureCharmUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureCharmUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesFlavorPureCharmUp;
   Float_t         MET_phi_jesFlavorPureCharmUp;
   Float_t         Jet_pt_jesFlavorPureBottomUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureBottomUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesFlavorPureBottomUp;
   Float_t         MET_phi_jesFlavorPureBottomUp;
   Float_t         Jet_pt_jesTimeRunBUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesTimeRunBUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesTimeRunBUp;
   Float_t         MET_phi_jesTimeRunBUp;
   Float_t         Jet_pt_jesTimeRunCUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesTimeRunCUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesTimeRunCUp;
   Float_t         MET_phi_jesTimeRunCUp;
   Float_t         Jet_pt_jesTimeRunDEUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesTimeRunDEUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesTimeRunDEUp;
   Float_t         MET_phi_jesTimeRunDEUp;
   Float_t         Jet_pt_jesTimeRunFUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesTimeRunFUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesTimeRunFUp;
   Float_t         MET_phi_jesTimeRunFUp;
   Float_t         Jet_pt_jesCorrelationGroupMPFInSituUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupMPFInSituUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupMPFInSituUp;
   Float_t         MET_phi_jesCorrelationGroupMPFInSituUp;
   Float_t         Jet_pt_jesCorrelationGroupIntercalibrationUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupIntercalibrationUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupIntercalibrationUp;
   Float_t         MET_phi_jesCorrelationGroupIntercalibrationUp;
   Float_t         Jet_pt_jesCorrelationGroupbJESUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupbJESUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupbJESUp;
   Float_t         MET_phi_jesCorrelationGroupbJESUp;
   Float_t         Jet_pt_jesCorrelationGroupFlavorUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupFlavorUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupFlavorUp;
   Float_t         MET_phi_jesCorrelationGroupFlavorUp;
   Float_t         Jet_pt_jesCorrelationGroupUncorrelatedUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupUncorrelatedUp[njetmax];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupUncorrelatedUp;
   Float_t         MET_phi_jesCorrelationGroupUncorrelatedUp;
   Float_t         MET_pt_unclustEnUp;
   Float_t         MET_phi_unclustEnUp;
   Float_t         Jet_pt_jerDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jerDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jmrDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jmsDown[njetmax];   //[nJet]
   Float_t         MET_pt_jerDown;
   Float_t         MET_phi_jerDown;
   Float_t         Jet_pt_jesAbsoluteStatDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteStatDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesAbsoluteStatDown;
   Float_t         MET_phi_jesAbsoluteStatDown;
   Float_t         Jet_pt_jesAbsoluteScaleDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteScaleDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesAbsoluteScaleDown;
   Float_t         MET_phi_jesAbsoluteScaleDown;
   Float_t         Jet_pt_jesAbsoluteFlavMapDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteFlavMapDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesAbsoluteFlavMapDown;
   Float_t         MET_phi_jesAbsoluteFlavMapDown;
   Float_t         Jet_pt_jesAbsoluteMPFBiasDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteMPFBiasDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesAbsoluteMPFBiasDown;
   Float_t         MET_phi_jesAbsoluteMPFBiasDown;
   Float_t         Jet_pt_jesFragmentationDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesFragmentationDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesFragmentationDown;
   Float_t         MET_phi_jesFragmentationDown;
   Float_t         Jet_pt_jesSinglePionECALDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesSinglePionECALDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesSinglePionECALDown;
   Float_t         MET_phi_jesSinglePionECALDown;
   Float_t         Jet_pt_jesSinglePionHCALDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesSinglePionHCALDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesSinglePionHCALDown;
   Float_t         MET_phi_jesSinglePionHCALDown;
   Float_t         Jet_pt_jesFlavorQCDDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesFlavorQCDDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesFlavorQCDDown;
   Float_t         MET_phi_jesFlavorQCDDown;
   Float_t         Jet_pt_jesTimePtEtaDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesTimePtEtaDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesTimePtEtaDown;
   Float_t         MET_phi_jesTimePtEtaDown;
   Float_t         Jet_pt_jesRelativeJEREC1Down[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativeJEREC1Down[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativeJEREC1Down;
   Float_t         MET_phi_jesRelativeJEREC1Down;
   Float_t         Jet_pt_jesRelativeJEREC2Down[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativeJEREC2Down[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativeJEREC2Down;
   Float_t         MET_phi_jesRelativeJEREC2Down;
   Float_t         Jet_pt_jesRelativeJERHFDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativeJERHFDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativeJERHFDown;
   Float_t         MET_phi_jesRelativeJERHFDown;
   Float_t         Jet_pt_jesRelativePtBBDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativePtBBDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativePtBBDown;
   Float_t         MET_phi_jesRelativePtBBDown;
   Float_t         Jet_pt_jesRelativePtEC1Down[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativePtEC1Down[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativePtEC1Down;
   Float_t         MET_phi_jesRelativePtEC1Down;
   Float_t         Jet_pt_jesRelativePtEC2Down[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativePtEC2Down[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativePtEC2Down;
   Float_t         MET_phi_jesRelativePtEC2Down;
   Float_t         Jet_pt_jesRelativePtHFDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativePtHFDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativePtHFDown;
   Float_t         MET_phi_jesRelativePtHFDown;
   Float_t         Jet_pt_jesRelativeBalDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativeBalDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativeBalDown;
   Float_t         MET_phi_jesRelativeBalDown;
   Float_t         Jet_pt_jesRelativeSampleDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativeSampleDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativeSampleDown;
   Float_t         MET_phi_jesRelativeSampleDown;
   Float_t         Jet_pt_jesRelativeFSRDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativeFSRDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativeFSRDown;
   Float_t         MET_phi_jesRelativeFSRDown;
   Float_t         Jet_pt_jesRelativeStatFSRDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativeStatFSRDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativeStatFSRDown;
   Float_t         MET_phi_jesRelativeStatFSRDown;
   Float_t         Jet_pt_jesRelativeStatECDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativeStatECDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativeStatECDown;
   Float_t         MET_phi_jesRelativeStatECDown;
   Float_t         Jet_pt_jesRelativeStatHFDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesRelativeStatHFDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesRelativeStatHFDown;
   Float_t         MET_phi_jesRelativeStatHFDown;
   Float_t         Jet_pt_jesPileUpDataMCDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesPileUpDataMCDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesPileUpDataMCDown;
   Float_t         MET_phi_jesPileUpDataMCDown;
   Float_t         Jet_pt_jesPileUpPtRefDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtRefDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesPileUpPtRefDown;
   Float_t         MET_phi_jesPileUpPtRefDown;
   Float_t         Jet_pt_jesPileUpPtBBDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtBBDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesPileUpPtBBDown;
   Float_t         MET_phi_jesPileUpPtBBDown;
   Float_t         Jet_pt_jesPileUpPtEC1Down[njetmax];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtEC1Down[njetmax];   //[nJet]
   Float_t         MET_pt_jesPileUpPtEC1Down;
   Float_t         MET_phi_jesPileUpPtEC1Down;
   Float_t         Jet_pt_jesPileUpPtEC2Down[njetmax];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtEC2Down[njetmax];   //[nJet]
   Float_t         MET_pt_jesPileUpPtEC2Down;
   Float_t         MET_phi_jesPileUpPtEC2Down;
   Float_t         Jet_pt_jesPileUpPtHFDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtHFDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesPileUpPtHFDown;
   Float_t         MET_phi_jesPileUpPtHFDown;
   Float_t         Jet_pt_jesPileUpMuZeroDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesPileUpMuZeroDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesPileUpMuZeroDown;
   Float_t         MET_phi_jesPileUpMuZeroDown;
   Float_t         Jet_pt_jesPileUpEnvelopeDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesPileUpEnvelopeDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesPileUpEnvelopeDown;
   Float_t         MET_phi_jesPileUpEnvelopeDown;
   Float_t         Jet_pt_jesSubTotalPileUpDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesSubTotalPileUpDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesSubTotalPileUpDown;
   Float_t         MET_phi_jesSubTotalPileUpDown;
   Float_t         Jet_pt_jesSubTotalRelativeDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesSubTotalRelativeDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesSubTotalRelativeDown;
   Float_t         MET_phi_jesSubTotalRelativeDown;
   Float_t         Jet_pt_jesSubTotalPtDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesSubTotalPtDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesSubTotalPtDown;
   Float_t         MET_phi_jesSubTotalPtDown;
   Float_t         Jet_pt_jesSubTotalScaleDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesSubTotalScaleDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesSubTotalScaleDown;
   Float_t         MET_phi_jesSubTotalScaleDown;
   Float_t         Jet_pt_jesSubTotalAbsoluteDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesSubTotalAbsoluteDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesSubTotalAbsoluteDown;
   Float_t         MET_phi_jesSubTotalAbsoluteDown;
   Float_t         Jet_pt_jesSubTotalMCDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesSubTotalMCDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesSubTotalMCDown;
   Float_t         MET_phi_jesSubTotalMCDown;
   Float_t         Jet_pt_jesTotalDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesTotalDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesTotalDown;
   Float_t         MET_phi_jesTotalDown;
   Float_t         Jet_pt_jesTotalNoFlavorDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesTotalNoFlavorDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesTotalNoFlavorDown;
   Float_t         MET_phi_jesTotalNoFlavorDown;
   Float_t         Jet_pt_jesTotalNoTimeDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesTotalNoTimeDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesTotalNoTimeDown;
   Float_t         MET_phi_jesTotalNoTimeDown;
   Float_t         Jet_pt_jesTotalNoFlavorNoTimeDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesTotalNoFlavorNoTimeDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesTotalNoFlavorNoTimeDown;
   Float_t         MET_phi_jesTotalNoFlavorNoTimeDown;
   Float_t         Jet_pt_jesFlavorZJetDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesFlavorZJetDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesFlavorZJetDown;
   Float_t         MET_phi_jesFlavorZJetDown;
   Float_t         Jet_pt_jesFlavorPhotonJetDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesFlavorPhotonJetDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesFlavorPhotonJetDown;
   Float_t         MET_phi_jesFlavorPhotonJetDown;
   Float_t         Jet_pt_jesFlavorPureGluonDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureGluonDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesFlavorPureGluonDown;
   Float_t         MET_phi_jesFlavorPureGluonDown;
   Float_t         Jet_pt_jesFlavorPureQuarkDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureQuarkDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesFlavorPureQuarkDown;
   Float_t         MET_phi_jesFlavorPureQuarkDown;
   Float_t         Jet_pt_jesFlavorPureCharmDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureCharmDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesFlavorPureCharmDown;
   Float_t         MET_phi_jesFlavorPureCharmDown;
   Float_t         Jet_pt_jesFlavorPureBottomDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureBottomDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesFlavorPureBottomDown;
   Float_t         MET_phi_jesFlavorPureBottomDown;
   Float_t         Jet_pt_jesTimeRunBDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesTimeRunBDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesTimeRunBDown;
   Float_t         MET_phi_jesTimeRunBDown;
   Float_t         Jet_pt_jesTimeRunCDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesTimeRunCDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesTimeRunCDown;
   Float_t         MET_phi_jesTimeRunCDown;
   Float_t         Jet_pt_jesTimeRunDEDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesTimeRunDEDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesTimeRunDEDown;
   Float_t         MET_phi_jesTimeRunDEDown;
   Float_t         Jet_pt_jesTimeRunFDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesTimeRunFDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesTimeRunFDown;
   Float_t         MET_phi_jesTimeRunFDown;
   Float_t         Jet_pt_jesCorrelationGroupMPFInSituDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupMPFInSituDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupMPFInSituDown;
   Float_t         MET_phi_jesCorrelationGroupMPFInSituDown;
   Float_t         Jet_pt_jesCorrelationGroupIntercalibrationDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupIntercalibrationDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupIntercalibrationDown;
   Float_t         MET_phi_jesCorrelationGroupIntercalibrationDown;
   Float_t         Jet_pt_jesCorrelationGroupbJESDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupbJESDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupbJESDown;
   Float_t         MET_phi_jesCorrelationGroupbJESDown;
   Float_t         Jet_pt_jesCorrelationGroupFlavorDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupFlavorDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupFlavorDown;
   Float_t         MET_phi_jesCorrelationGroupFlavorDown;
   Float_t         Jet_pt_jesCorrelationGroupUncorrelatedDown[njetmax];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupUncorrelatedDown[njetmax];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupUncorrelatedDown;
   Float_t         MET_phi_jesCorrelationGroupUncorrelatedDown;
   Float_t         MET_pt_unclustEnDown;
   Float_t         MET_phi_unclustEnDown;
   
   Float_t         FatJet_pt_raw[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_nom[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_raw[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_nom[njetmax];   //[nFatJet]
   Float_t         FatJet_corr_JEC[njetmax];   //[nFatJet]
   Float_t         FatJet_corr_JER[njetmax];   //[nFatJet]
   Float_t         FatJet_corr_JMS[njetmax];   //[nFatJet]
   Float_t         FatJet_corr_JMR[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_raw[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_nom[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_corr_JMR[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_corr_JMS[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_corr_PUPPI[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_tau21DDT_nom[njetmax];   //[nFatJet]

   Float_t         FatJet_pt_jerUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jerUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jmrUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jmsUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jerUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jmrUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jmsUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesAbsoluteStatUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesAbsoluteStatUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesAbsoluteStatUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesAbsoluteScaleUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesAbsoluteScaleUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesAbsoluteScaleUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesAbsoluteFlavMapUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesAbsoluteFlavMapUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesAbsoluteFlavMapUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesAbsoluteMPFBiasUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesAbsoluteMPFBiasUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesAbsoluteMPFBiasUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesFragmentationUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesFragmentationUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFragmentationUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesSinglePionECALUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesSinglePionECALUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSinglePionECALUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesSinglePionHCALUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesSinglePionHCALUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSinglePionHCALUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorQCDUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorQCDUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorQCDUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesTimePtEtaUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesTimePtEtaUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimePtEtaUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeJEREC1Up[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeJEREC1Up[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeJEREC1Up[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeJEREC2Up[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeJEREC2Up[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeJEREC2Up[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeJERHFUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeJERHFUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeJERHFUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativePtBBUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativePtBBUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativePtBBUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativePtEC1Up[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativePtEC1Up[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativePtEC1Up[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativePtEC2Up[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativePtEC2Up[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativePtEC2Up[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativePtHFUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativePtHFUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativePtHFUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeBalUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeBalUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeBalUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeSampleUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeSampleUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeSampleUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeFSRUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeFSRUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeFSRUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeStatFSRUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeStatFSRUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeStatFSRUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeStatECUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeStatECUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeStatECUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeStatHFUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeStatHFUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeStatHFUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpDataMCUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpDataMCUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpDataMCUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtRefUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtRefUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtRefUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtBBUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtBBUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtBBUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtEC1Up[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtEC1Up[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtEC1Up[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtEC2Up[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtEC2Up[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtEC2Up[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtHFUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtHFUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtHFUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpMuZeroUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpMuZeroUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpMuZeroUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpEnvelopeUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpEnvelopeUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpEnvelopeUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalPileUpUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalPileUpUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalPileUpUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalRelativeUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalRelativeUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalRelativeUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalPtUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalPtUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalPtUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalScaleUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalScaleUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalScaleUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalAbsoluteUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalAbsoluteUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalAbsoluteUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalMCUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalMCUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalMCUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesTotalUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesTotalUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTotalUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesTotalNoFlavorUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesTotalNoFlavorUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTotalNoFlavorUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesTotalNoTimeUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesTotalNoTimeUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTotalNoTimeUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesTotalNoFlavorNoTimeUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesTotalNoFlavorNoTimeUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTotalNoFlavorNoTimeUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorZJetUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorZJetUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorZJetUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPhotonJetUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPhotonJetUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPhotonJetUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPureGluonUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPureGluonUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPureGluonUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPureQuarkUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPureQuarkUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPureQuarkUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPureCharmUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPureCharmUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPureCharmUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPureBottomUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPureBottomUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPureBottomUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesTimeRunBUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesTimeRunBUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimeRunBUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesTimeRunCUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesTimeRunCUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimeRunCUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesTimeRunDEUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesTimeRunDEUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimeRunDEUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesTimeRunFUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesTimeRunFUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimeRunFUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupMPFInSituUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupMPFInSituUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupMPFInSituUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupIntercalibrationUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupIntercalibrationUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupIntercalibrationUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupbJESUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupbJESUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupbJESUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupFlavorUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupFlavorUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupFlavorUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupUncorrelatedUp[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupUncorrelatedUp[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupUncorrelatedUp[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jerDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jerDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jmrDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jmsDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jerDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jmrDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jmsDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesAbsoluteStatDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesAbsoluteStatDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesAbsoluteStatDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesAbsoluteScaleDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesAbsoluteScaleDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesAbsoluteScaleDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesAbsoluteFlavMapDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesAbsoluteFlavMapDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesAbsoluteFlavMapDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesAbsoluteMPFBiasDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesAbsoluteMPFBiasDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesAbsoluteMPFBiasDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesFragmentationDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesFragmentationDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFragmentationDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesSinglePionECALDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesSinglePionECALDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSinglePionECALDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesSinglePionHCALDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesSinglePionHCALDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSinglePionHCALDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorQCDDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorQCDDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorQCDDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesTimePtEtaDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesTimePtEtaDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimePtEtaDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeJEREC1Down[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeJEREC1Down[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeJEREC1Down[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeJEREC2Down[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeJEREC2Down[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeJEREC2Down[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeJERHFDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeJERHFDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeJERHFDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativePtBBDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativePtBBDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativePtBBDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativePtEC1Down[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativePtEC1Down[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativePtEC1Down[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativePtEC2Down[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativePtEC2Down[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativePtEC2Down[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativePtHFDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativePtHFDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativePtHFDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeBalDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeBalDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeBalDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeSampleDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeSampleDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeSampleDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeFSRDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeFSRDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeFSRDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeStatFSRDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeStatFSRDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeStatFSRDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeStatECDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeStatECDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeStatECDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeStatHFDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeStatHFDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeStatHFDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpDataMCDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpDataMCDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpDataMCDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtRefDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtRefDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtRefDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtBBDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtBBDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtBBDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtEC1Down[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtEC1Down[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtEC1Down[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtEC2Down[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtEC2Down[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtEC2Down[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtHFDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtHFDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtHFDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpMuZeroDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpMuZeroDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpMuZeroDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpEnvelopeDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpEnvelopeDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpEnvelopeDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalPileUpDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalPileUpDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalPileUpDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalRelativeDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalRelativeDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalRelativeDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalPtDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalPtDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalPtDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalScaleDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalScaleDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalScaleDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalAbsoluteDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalAbsoluteDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalAbsoluteDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalMCDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalMCDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalMCDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesTotalDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesTotalDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTotalDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesTotalNoFlavorDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesTotalNoFlavorDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTotalNoFlavorDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesTotalNoTimeDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesTotalNoTimeDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTotalNoTimeDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesTotalNoFlavorNoTimeDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesTotalNoFlavorNoTimeDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTotalNoFlavorNoTimeDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorZJetDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorZJetDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorZJetDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPhotonJetDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPhotonJetDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPhotonJetDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPureGluonDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPureGluonDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPureGluonDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPureQuarkDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPureQuarkDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPureQuarkDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPureCharmDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPureCharmDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPureCharmDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPureBottomDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPureBottomDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPureBottomDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesTimeRunBDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesTimeRunBDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimeRunBDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesTimeRunCDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesTimeRunCDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimeRunCDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesTimeRunDEDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesTimeRunDEDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimeRunDEDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesTimeRunFDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesTimeRunFDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimeRunFDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupMPFInSituDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupMPFInSituDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupMPFInSituDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupIntercalibrationDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupIntercalibrationDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupIntercalibrationDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupbJESDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupbJESDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupbJESDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupFlavorDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupFlavorDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupFlavorDown[njetmax];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupUncorrelatedDown[njetmax];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupUncorrelatedDown[njetmax];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupUncorrelatedDown[njetmax];   //[nFatJet]
   Float_t         Jet_CSVbtagSF[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_up[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_down[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_shape[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_shape_up_jes[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_shape_down_jes[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_shape_up_lf[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_shape_down_lf[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_shape_up_hf[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_shape_down_hf[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_shape_up_hfstats1[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_shape_down_hfstats1[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_shape_up_hfstats2[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_shape_down_hfstats2[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_shape_up_lfstats1[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_shape_down_lfstats1[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_shape_up_lfstats2[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_shape_down_lfstats2[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_shape_up_cferr1[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_shape_down_cferr1[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_shape_up_cferr2[njetmax];   //[nJet]
   Float_t         Jet_CSVbtagSF_shape_down_cferr2[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_up[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_down[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_shape[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_shape_up_jes[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_shape_down_jes[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_shape_up_lf[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_shape_down_lf[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_shape_up_hf[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_shape_down_hf[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_shape_up_hfstats1[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_shape_down_hfstats1[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_shape_up_hfstats2[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_shape_down_hfstats2[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_shape_up_lfstats1[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_shape_down_lfstats1[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_shape_up_lfstats2[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_shape_down_lfstats2[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_shape_up_cferr1[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_shape_down_cferr1[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_shape_up_cferr2[njetmax];   //[nJet]
   Float_t         Jet_deepCSVbtagSF_shape_down_cferr2[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_up[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_down[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_shape[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_shape_up_jes[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_shape_down_jes[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_shape_up_lf[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_shape_down_lf[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_shape_up_hf[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_shape_down_hf[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_shape_up_hfstats1[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_shape_down_hfstats1[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_shape_up_hfstats2[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_shape_down_hfstats2[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_shape_up_lfstats1[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_shape_down_lfstats1[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_shape_up_lfstats2[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_shape_down_lfstats2[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_shape_up_cferr1[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_shape_down_cferr1[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_shape_up_cferr2[njetmax];   //[nJet]
   Float_t         Jet_deepflavbtagSF_shape_down_cferr2[njetmax];   //[nJet]
   Float_t         puWeight;
   Float_t         puWeightUp;
   Float_t         puWeightDown;
   Float_t         PrefireWeight;
   Float_t         PrefireWeight_Up;
   Float_t         PrefireWeight_Down;
   
   Float_t		   Event_Ht;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_btagWeight_CSVV2;   //!
   TBranch        *b_btagWeight_DeepCSVB;   //!
   TBranch        *b_CaloMET_phi;   //!
   TBranch        *b_CaloMET_pt;   //!
   TBranch        *b_CaloMET_sumEt;   //!
   TBranch        *b_ChsMET_phi;   //!
   TBranch        *b_ChsMET_pt;   //!
   TBranch        *b_ChsMET_sumEt;   //!
   TBranch        *b_nElectron;   //!
   TBranch        *b_Electron_deltaEtaSC;   //!
   TBranch        *b_Electron_dr03EcalRecHitSumEt;   //!
   TBranch        *b_Electron_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_Electron_dr03TkSumPt;   //!
   TBranch        *b_Electron_dr03TkSumPtHEEP;   //!
   TBranch        *b_Electron_dxy;   //!
   TBranch        *b_Electron_dxyErr;   //!
   TBranch        *b_Electron_dz;   //!
   TBranch        *b_Electron_dzErr;   //!
   TBranch        *b_Electron_eCorr;   //!
   TBranch        *b_Electron_eInvMinusPInv;   //!
   TBranch        *b_Electron_energyErr;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_hoe;   //!
   TBranch        *b_Electron_ip3d;   //!
   TBranch        *b_Electron_jetRelIso;   //!
   TBranch        *b_Electron_mass;   //!
   TBranch        *b_Electron_miniPFRelIso_all;   //!
   TBranch        *b_Electron_miniPFRelIso_chg;   //!
   TBranch        *b_Electron_mvaFall17V1Iso;   //!
   TBranch        *b_Electron_mvaFall17V1noIso;   //!
   TBranch        *b_Electron_mvaFall17V2Iso;   //!
   TBranch        *b_Electron_mvaFall17V2noIso;   //!
   TBranch        *b_Electron_pfRelIso03_all;   //!
   TBranch        *b_Electron_pfRelIso03_chg;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_r9;   //!
   TBranch        *b_Electron_sieie;   //!
   TBranch        *b_Electron_sip3d;   //!
   TBranch        *b_Electron_mvaTTH;   //!
   TBranch        *b_Electron_charge;   //!
   TBranch        *b_Electron_cutBased;   //!
   TBranch        *b_Electron_cutBased_Fall17_V1;   //!
   TBranch        *b_Electron_jetIdx;   //!
   TBranch        *b_Electron_pdgId;   //!
   TBranch        *b_Electron_photonIdx;   //!
   TBranch        *b_Electron_tightCharge;   //!
   TBranch        *b_Electron_vidNestedWPBitmap;   //!
   TBranch        *b_Electron_convVeto;   //!
   TBranch        *b_Electron_cutBased_HEEP;   //!
   TBranch        *b_Electron_isPFcand;   //!
   TBranch        *b_Electron_lostHits;   //!
   TBranch        *b_Electron_mvaFall17V1Iso_WP80;   //!
   TBranch        *b_Electron_mvaFall17V1Iso_WP90;   //!
   TBranch        *b_Electron_mvaFall17V1Iso_WPL;   //!
   TBranch        *b_Electron_mvaFall17V1noIso_WP80;   //!
   TBranch        *b_Electron_mvaFall17V1noIso_WP90;   //!
   TBranch        *b_Electron_mvaFall17V1noIso_WPL;   //!
   TBranch        *b_Electron_mvaFall17V2Iso_WP80;   //!
   TBranch        *b_Electron_mvaFall17V2Iso_WP90;   //!
   TBranch        *b_Electron_mvaFall17V2Iso_WPL;   //!
   TBranch        *b_Electron_mvaFall17V2noIso_WP80;   //!
   TBranch        *b_Electron_mvaFall17V2noIso_WP90;   //!
   TBranch        *b_Electron_mvaFall17V2noIso_WPL;   //!
   TBranch        *b_Flag_ecalBadCalibFilterV2;   //!
   TBranch        *b_nFatJet;   //!
   TBranch        *b_FatJet_area;   //!
   TBranch        *b_FatJet_btagCMVA;   //!
   TBranch        *b_FatJet_btagCSVV2;   //!
   TBranch        *b_FatJet_btagDeepB;   //!
   TBranch        *b_FatJet_btagHbb;   //!
   TBranch        *b_FatJet_deepTagMD_H4qvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_HbbvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_TvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_WvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_ZHbbvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_ZHccvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_ZbbvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_ZvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_bbvsLight;   //!
   TBranch        *b_FatJet_deepTagMD_ccvsLight;   //!
   TBranch        *b_FatJet_deepTag_H;   //!
   TBranch        *b_FatJet_deepTag_QCD;   //!
   TBranch        *b_FatJet_deepTag_QCDothers;   //!
   TBranch        *b_FatJet_deepTag_TvsQCD;   //!
   TBranch        *b_FatJet_deepTag_WvsQCD;   //!
   TBranch        *b_FatJet_deepTag_ZvsQCD;   //!
   TBranch        *b_FatJet_eta;   //!
   TBranch        *b_FatJet_mass;   //!
   TBranch        *b_FatJet_msoftdrop;   //!
   TBranch        *b_FatJet_n2b1;   //!
   TBranch        *b_FatJet_n3b1;   //!
   TBranch        *b_FatJet_phi;   //!
   TBranch        *b_FatJet_pt;   //!
   TBranch        *b_FatJet_rawFactor;   //!
   TBranch        *b_FatJet_tau1;   //!
   TBranch        *b_FatJet_tau2;   //!
   TBranch        *b_FatJet_tau3;   //!
   TBranch        *b_FatJet_tau4;   //!
   TBranch        *b_FatJet_jetId;   //!
   TBranch        *b_FatJet_subJetIdx1;   //!
   TBranch        *b_FatJet_subJetIdx2;   //!
   TBranch        *b_nGenJetAK8;   //!
   TBranch        *b_GenJetAK8_eta;   //!
   TBranch        *b_GenJetAK8_mass;   //!
   TBranch        *b_GenJetAK8_phi;   //!
   TBranch        *b_GenJetAK8_pt;   //!
   TBranch        *b_nGenJet;   //!
   TBranch        *b_GenJet_eta;   //!
   TBranch        *b_GenJet_mass;   //!
   TBranch        *b_GenJet_phi;   //!
   TBranch        *b_GenJet_pt;   //!
   TBranch        *b_nGenPart;   //!
   TBranch        *b_GenPart_eta;   //!
   TBranch        *b_GenPart_mass;   //!
   TBranch        *b_GenPart_phi;   //!
   TBranch        *b_GenPart_pt;   //!
   TBranch        *b_GenPart_genPartIdxMother;   //!
   TBranch        *b_GenPart_pdgId;   //!
   TBranch        *b_GenPart_status;   //!
   TBranch        *b_GenPart_statusFlags;   //!
   TBranch        *b_nSubGenJetAK8;   //!
   TBranch        *b_SubGenJetAK8_eta;   //!
   TBranch        *b_SubGenJetAK8_mass;   //!
   TBranch        *b_SubGenJetAK8_phi;   //!
   TBranch        *b_SubGenJetAK8_pt;   //!
   TBranch        *b_Generator_binvar;   //!
   TBranch        *b_Generator_scalePDF;   //!
   TBranch        *b_Generator_weight;   //!
   TBranch        *b_Generator_x1;   //!
   TBranch        *b_Generator_x2;   //!
   TBranch        *b_Generator_xpdf1;   //!
   TBranch        *b_Generator_xpdf2;   //!
   TBranch        *b_Generator_id1;   //!
   TBranch        *b_Generator_id2;   //!
   TBranch        *b_nGenVisTau;   //!
   TBranch        *b_GenVisTau_eta;   //!
   TBranch        *b_GenVisTau_mass;   //!
   TBranch        *b_GenVisTau_phi;   //!
   TBranch        *b_GenVisTau_pt;   //!
   TBranch        *b_GenVisTau_charge;   //!
   TBranch        *b_GenVisTau_genPartIdxMother;   //!
   TBranch        *b_GenVisTau_status;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_LHEWeight_originalXWGTUP;   //!
   TBranch        *b_nLHEPdfWeight;   //!
   TBranch        *b_LHEPdfWeight;   //!
   TBranch        *b_nLHEScaleWeight;   //!
   TBranch        *b_LHEScaleWeight;   //!
   TBranch        *b_nPSWeight;   //!
   TBranch        *b_PSWeight;   //!
   TBranch        *b_nIsoTrack;   //!
   TBranch        *b_IsoTrack_dxy;   //!
   TBranch        *b_IsoTrack_dz;   //!
   TBranch        *b_IsoTrack_eta;   //!
   TBranch        *b_IsoTrack_pfRelIso03_all;   //!
   TBranch        *b_IsoTrack_pfRelIso03_chg;   //!
   TBranch        *b_IsoTrack_phi;   //!
   TBranch        *b_IsoTrack_pt;   //!
   TBranch        *b_IsoTrack_miniPFRelIso_all;   //!
   TBranch        *b_IsoTrack_miniPFRelIso_chg;   //!
   TBranch        *b_IsoTrack_fromPV;   //!
   TBranch        *b_IsoTrack_pdgId;   //!
   TBranch        *b_IsoTrack_isHighPurityTrack;   //!
   TBranch        *b_IsoTrack_isPFcand;   //!
   TBranch        *b_IsoTrack_isFromLostTrack;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_area;   //!
   TBranch        *b_Jet_btagCMVA;   //!
   TBranch        *b_Jet_btagCSVV2;   //!
   TBranch        *b_Jet_btagDeepB;   //!
   TBranch        *b_Jet_btagDeepC;   //!
   TBranch        *b_Jet_btagDeepFlavB;   //!
   TBranch        *b_Jet_chEmEF;   //!
   TBranch        *b_Jet_chHEF;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_muEF;   //!
   TBranch        *b_Jet_neEmEF;   //!
   TBranch        *b_Jet_neHEF;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_qgl;   //!
   TBranch        *b_Jet_rawFactor;   //!
   TBranch        *b_Jet_bRegCorr;   //!
   TBranch        *b_Jet_bRegRes;   //!
   TBranch        *b_Jet_electronIdx1;   //!
   TBranch        *b_Jet_electronIdx2;   //!
   TBranch        *b_Jet_jetId;   //!
   TBranch        *b_Jet_muonIdx1;   //!
   TBranch        *b_Jet_muonIdx2;   //!
   TBranch        *b_Jet_nConstituents;   //!
   TBranch        *b_Jet_nElectrons;   //!
   TBranch        *b_Jet_nMuons;   //!
   TBranch        *b_Jet_puId;   //!
   TBranch        *b_LHE_HT;   //!
   TBranch        *b_LHE_HTIncoming;   //!
   TBranch        *b_LHE_Vpt;   //!
   TBranch        *b_LHE_Njets;   //!
   TBranch        *b_LHE_Nb;   //!
   TBranch        *b_LHE_Nc;   //!
   TBranch        *b_LHE_Nuds;   //!
   TBranch        *b_LHE_Nglu;   //!
   TBranch        *b_LHE_NpNLO;   //!
   TBranch        *b_LHE_NpLO;   //!
   TBranch        *b_nLHEPart;   //!
   TBranch        *b_LHEPart_pt;   //!
   TBranch        *b_LHEPart_eta;   //!
   TBranch        *b_LHEPart_phi;   //!
   TBranch        *b_LHEPart_mass;   //!
   TBranch        *b_LHEPart_pdgId;   //!
   TBranch        *b_METFixEE2017_MetUnclustEnUpDeltaX;   //!
   TBranch        *b_METFixEE2017_MetUnclustEnUpDeltaY;   //!
   TBranch        *b_METFixEE2017_covXX;   //!
   TBranch        *b_METFixEE2017_covXY;   //!
   TBranch        *b_METFixEE2017_covYY;   //!
   TBranch        *b_METFixEE2017_phi;   //!
   TBranch        *b_METFixEE2017_pt;   //!
   TBranch        *b_METFixEE2017_significance;   
   TBranch        *b_METFixEE2017_sumEt;   //!
   TBranch        *b_GenMET_phi;   //!
   TBranch        *b_GenMET_pt;   //!
   TBranch        *b_MET_MetUnclustEnUpDeltaX;   //!
   TBranch        *b_MET_MetUnclustEnUpDeltaY;   //!
   TBranch        *b_MET_covXX;   //!
   TBranch        *b_MET_covXY;   //!
   TBranch        *b_MET_covYY;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_MET_pt;   //!
   TBranch        *b_MET_significance;   //!
   TBranch        *b_MET_sumEt;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_Muon_dxy;   //!
   TBranch        *b_Muon_dxyErr;   //!
   TBranch        *b_Muon_dz;   //!
   TBranch        *b_Muon_dzErr;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_ip3d;   //!
   TBranch        *b_Muon_jetPtRelv2;   //!
   TBranch        *b_Muon_jetRelIso;   //!
   TBranch        *b_Muon_mass;   //!
   TBranch        *b_Muon_miniPFRelIso_all;   //!
   TBranch        *b_Muon_miniPFRelIso_chg;   //!
   TBranch        *b_Muon_pfRelIso03_all;   //!
   TBranch        *b_Muon_pfRelIso03_chg;   //!
   TBranch        *b_Muon_pfRelIso04_all;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_ptErr;   //!
   TBranch        *b_Muon_segmentComp;   //!
   TBranch        *b_Muon_sip3d;   //!
   TBranch        *b_Muon_mvaTTH;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_jetIdx;   //!
   TBranch        *b_Muon_nStations;   //!
   TBranch        *b_Muon_nTrackerLayers;   //!
   TBranch        *b_Muon_pdgId;   //!
   TBranch        *b_Muon_tightCharge;   //!
   TBranch        *b_Muon_highPtId;   //!
   TBranch        *b_Muon_inTimeMuon;   //!
   TBranch        *b_Muon_isGlobal;   //!
   TBranch        *b_Muon_isPFcand;   //!
   TBranch        *b_Muon_isTracker;   //!
   TBranch        *b_Muon_mediumId;   //!
   TBranch        *b_Muon_mediumPromptId;   //!
   TBranch        *b_Muon_miniIsoId;   //!
   TBranch        *b_Muon_multiIsoId;   //!
   TBranch        *b_Muon_mvaId;   //!
   TBranch        *b_Muon_pfIsoId;   //!
   TBranch        *b_Muon_softId;   //!
   TBranch        *b_Muon_softMvaId;   //!
   TBranch        *b_Muon_tightId;   //!
   TBranch        *b_Muon_tkIsoId;   //!
   TBranch        *b_Muon_triggerIdLoose;   //!
   TBranch        *b_nPhoton;   //!
   TBranch        *b_Photon_eCorr;   //!
   TBranch        *b_Photon_energyErr;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_hoe;   //!
   TBranch        *b_Photon_mass;   //!
   TBranch        *b_Photon_mvaID;   //!
   TBranch        *b_Photon_mvaIDV1;   //!
   TBranch        *b_Photon_pfRelIso03_all;   //!
   TBranch        *b_Photon_pfRelIso03_chg;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_r9;   //!
   TBranch        *b_Photon_sieie;   //!
   TBranch        *b_Photon_charge;   //!
   TBranch        *b_Photon_cutBasedBitmap;   //!
   TBranch        *b_Photon_cutBasedV1Bitmap;   //!
   TBranch        *b_Photon_electronIdx;   //!
   TBranch        *b_Photon_jetIdx;   //!
   TBranch        *b_Photon_pdgId;   //!
   TBranch        *b_Photon_vidNestedWPBitmap;   //!
   TBranch        *b_Photon_electronVeto;   //!
   TBranch        *b_Photon_isScEtaEB;   //!
   TBranch        *b_Photon_isScEtaEE;   //!
   TBranch        *b_Photon_mvaID_WP80;   //!
   TBranch        *b_Photon_mvaID_WP90;   //!
   TBranch        *b_Photon_pixelSeed;   //!
   TBranch        *b_Pileup_nTrueInt;   //!
   TBranch        *b_Pileup_nPU;   //!
   TBranch        *b_Pileup_sumEOOT;   //!
   TBranch        *b_Pileup_sumLOOT;   //!
   TBranch        *b_PuppiMET_phi;   //!
   TBranch        *b_PuppiMET_pt;   //!
   TBranch        *b_PuppiMET_sumEt;   //!
   TBranch        *b_RawMET_phi;   //!
   TBranch        *b_RawMET_pt;   //!
   TBranch        *b_RawMET_sumEt;   //!
   TBranch        *b_fixedGridRhoFastjetAll;   //!
   TBranch        *b_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_fixedGridRhoFastjetCentralNeutral;   //!
   TBranch        *b_nGenDressedLepton;   //!
   TBranch        *b_GenDressedLepton_eta;   //!
   TBranch        *b_GenDressedLepton_mass;   //!
   TBranch        *b_GenDressedLepton_phi;   //!
   TBranch        *b_GenDressedLepton_pt;   //!
   TBranch        *b_GenDressedLepton_pdgId;   //!
   TBranch        *b_nSoftActivityJet;   //!
   TBranch        *b_SoftActivityJet_eta;   //!
   TBranch        *b_SoftActivityJet_phi;   //!
   TBranch        *b_SoftActivityJet_pt;   //!
   TBranch        *b_SoftActivityJetHT;   //!
   TBranch        *b_SoftActivityJetHT10;   //!
   TBranch        *b_SoftActivityJetHT2;   //!
   TBranch        *b_SoftActivityJetHT5;   //!
   TBranch        *b_SoftActivityJetNjets10;   //!
   TBranch        *b_SoftActivityJetNjets2;   //!
   TBranch        *b_SoftActivityJetNjets5;   //!
   TBranch        *b_nSubJet;   //!
   TBranch        *b_SubJet_btagCMVA;   //!
   TBranch        *b_SubJet_btagCSVV2;   //!
   TBranch        *b_SubJet_btagDeepB;   //!
   TBranch        *b_SubJet_eta;   //!
   TBranch        *b_SubJet_mass;   //!
   TBranch        *b_SubJet_n2b1;   //!
   TBranch        *b_SubJet_n3b1;   //!
   TBranch        *b_SubJet_phi;   //!
   TBranch        *b_SubJet_pt;   //!
   TBranch        *b_SubJet_rawFactor;   //!
   TBranch        *b_SubJet_tau1;   //!
   TBranch        *b_SubJet_tau2;   //!
   TBranch        *b_SubJet_tau3;   //!
   TBranch        *b_SubJet_tau4;   //!
   TBranch        *b_nTau;   //!
   TBranch        *b_Tau_chargedIso;   //!
   TBranch        *b_Tau_dxy;   //!
   TBranch        *b_Tau_dz;   //!
   TBranch        *b_Tau_eta;   //!
   TBranch        *b_Tau_leadTkDeltaEta;   //!
   TBranch        *b_Tau_leadTkDeltaPhi;   //!
   TBranch        *b_Tau_leadTkPtOverTauPt;   //!
   TBranch        *b_Tau_mass;   //!
   TBranch        *b_Tau_neutralIso;   //!
   TBranch        *b_Tau_phi;   //!
   TBranch        *b_Tau_photonsOutsideSignalCone;   //!
   TBranch        *b_Tau_pt;   //!
   TBranch        *b_Tau_puCorr;   //!
   TBranch        *b_Tau_rawAntiEle;   //!
   TBranch        *b_Tau_rawIso;   //!
   TBranch        *b_Tau_rawIsodR03;   //!
   TBranch        *b_Tau_rawMVAnewDM2017v2;   //!
   TBranch        *b_Tau_rawMVAoldDM;   //!
   TBranch        *b_Tau_rawMVAoldDM2017v1;   //!
   TBranch        *b_Tau_rawMVAoldDM2017v2;   //!
   TBranch        *b_Tau_rawMVAoldDMdR032017v2;   //!
   TBranch        *b_Tau_charge;   //!
   TBranch        *b_Tau_decayMode;   //!
   TBranch        *b_Tau_jetIdx;   //!
   TBranch        *b_Tau_rawAntiEleCat;   //!
   TBranch        *b_Tau_idAntiEle;   //!
   TBranch        *b_Tau_idAntiMu;   //!
   TBranch        *b_Tau_idDecayMode;   //!
   TBranch        *b_Tau_idDecayModeNewDMs;   //!
   TBranch        *b_Tau_idMVAnewDM2017v2;   //!
   TBranch        *b_Tau_idMVAoldDM;   //!
   TBranch        *b_Tau_idMVAoldDM2017v1;   //!
   TBranch        *b_Tau_idMVAoldDM2017v2;   //!
   TBranch        *b_Tau_idMVAoldDMdR032017v2;   //!
   TBranch        *b_TkMET_phi;   //!
   TBranch        *b_TkMET_pt;   //!
   TBranch        *b_TkMET_sumEt;   //!
   TBranch        *b_nTrigObj;   //!
   TBranch        *b_TrigObj_pt;   //!
   TBranch        *b_TrigObj_eta;   //!
   TBranch        *b_TrigObj_phi;   //!
   TBranch        *b_TrigObj_l1pt;   //!
   TBranch        *b_TrigObj_l1pt_2;   //!
   TBranch        *b_TrigObj_l2pt;   //!
   TBranch        *b_TrigObj_id;   //!
   TBranch        *b_TrigObj_l1iso;   //!
   TBranch        *b_TrigObj_l1charge;   //!
   TBranch        *b_TrigObj_filterBits;   //!
   TBranch        *b_genTtbarId;   //!
   TBranch        *b_nOtherPV;   //!
   TBranch        *b_OtherPV_z;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_score;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_PV_npvsGood;   //!
   TBranch        *b_nSV;   //!
   TBranch        *b_SV_dlen;   //!
   TBranch        *b_SV_dlenSig;   //!
   TBranch        *b_SV_pAngle;   //!
   TBranch        *b_Electron_genPartIdx;   //!
   TBranch        *b_Electron_genPartFlav;   //!
   TBranch        *b_GenJetAK8_partonFlavour;   //!
   TBranch        *b_GenJetAK8_hadronFlavour;   //!
   TBranch        *b_GenJet_partonFlavour;   //!
   TBranch        *b_GenJet_hadronFlavour;   //!
   TBranch        *b_Jet_genJetIdx;   //!
   TBranch        *b_Jet_hadronFlavour;   //!
   TBranch        *b_Jet_partonFlavour;   //!
   TBranch        *b_Muon_genPartIdx;   //!
   TBranch        *b_Muon_genPartFlav;   //!
   TBranch        *b_Photon_genPartIdx;   //!
   TBranch        *b_Photon_genPartFlav;   //!
   TBranch        *b_MET_fiducialGenPhi;   //!
   TBranch        *b_MET_fiducialGenPt;   //!
   TBranch        *b_Electron_cleanmask;   //!
   TBranch        *b_Jet_cleanmask;   //!
   TBranch        *b_Muon_cleanmask;   //!
   TBranch        *b_Photon_cleanmask;   //!
   TBranch        *b_Tau_cleanmask;   //!
   TBranch        *b_SV_chi2;   //!
   TBranch        *b_SV_eta;   //!
   TBranch        *b_SV_mass;   //!
   TBranch        *b_SV_ndof;   //!
   TBranch        *b_SV_phi;   //!
   TBranch        *b_SV_pt;   //!
   TBranch        *b_SV_x;   //!
   TBranch        *b_SV_y;   //!
   TBranch        *b_SV_z;   //!
   TBranch        *b_Tau_genPartIdx;   //!
   TBranch        *b_Tau_genPartFlav;   //!
   TBranch        *b_L1simulation_step;   //!
   TBranch        *b_HLTriggerFirstPath;   //!
   TBranch        *b_HLT_AK8PFJet360_TrimMass30;   //!
   TBranch        *b_HLT_AK8PFJet380_TrimMass30;   //!
   TBranch        *b_HLT_AK8PFJet400_TrimMass30;   //!
   TBranch        *b_HLT_AK8PFJet420_TrimMass30;   //!
   TBranch        *b_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50;   //!
   TBranch        *b_HLT_AK8PFHT750_TrimMass50;   //!
   TBranch        *b_HLT_AK8PFHT800_TrimMass50;   //!
   TBranch        *b_HLT_AK8PFHT850_TrimMass50;   //!
   TBranch        *b_HLT_AK8PFHT900_TrimMass50;   //!
   TBranch        *b_HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20;   //!
   TBranch        *b_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20;   //!
   TBranch        *b_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087;   //!
   TBranch        *b_HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087;   //!
   TBranch        *b_HLT_AK8DiPFJet300_200_TrimMass30;   //!
   TBranch        *b_HLT_AK8DiPFJet280_200_TrimMass30;   //!
   TBranch        *b_HLT_CaloJet500_NoJetID;   //!
   TBranch        *b_HLT_CaloJet550_NoJetID;   //!
   TBranch        *b_HLT_Trimuon5_3p5_2_Upsilon_Muon;   //!
   TBranch        *b_HLT_DoubleEle25_CaloIdL_MW;   //!
   TBranch        *b_HLT_DoubleEle27_CaloIdL_MW;   //!
   TBranch        *b_HLT_DoubleEle33_CaloIdL_MW;   //!
   TBranch        *b_HLT_DoubleEle24_eta2p1_WPTight_Gsf;   //!
   TBranch        *b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350;   //!
   TBranch        *b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350;   //!
   TBranch        *b_HLT_Ele27_Ele37_CaloIdL_MW;   //!
   TBranch        *b_HLT_Mu27_Ele37_CaloIdL_MW;   //!
   TBranch        *b_HLT_Mu37_Ele27_CaloIdL_MW;   //!
   TBranch        *b_HLT_Mu37_TkMu27;   //!
   TBranch        *b_HLT_DoubleMu4_3_Bs;   //!
   TBranch        *b_HLT_DoubleMu4_3_Jpsi_Displaced;   //!
   TBranch        *b_HLT_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced;   //!
   TBranch        *b_HLT_DoubleMu3_Trk_Tau3mu;   //!
   TBranch        *b_HLT_DoubleMu4_PsiPrimeTrk_Displaced;   //!
   TBranch        *b_HLT_DoubleMu4_Mass8_DZ_PFHT350;   //!
   TBranch        *b_HLT_DoubleMu8_Mass8_PFHT350;   //!
   TBranch        *b_HLT_Mu3_PFJet40;   //!
   TBranch        *b_HLT_Mu7p5_L2Mu2_Jpsi;   //!
   TBranch        *b_HLT_Mu7p5_L2Mu2_Upsilon;   //!
   TBranch        *b_HLT_Mu7p5_Track2_Jpsi;   //!
   TBranch        *b_HLT_Mu7p5_Track3p5_Jpsi;   //!
   TBranch        *b_HLT_Mu7p5_Track7_Jpsi;   //!
   TBranch        *b_HLT_Mu7p5_Track2_Upsilon;   //!
   TBranch        *b_HLT_Mu7p5_Track3p5_Upsilon;   //!
   TBranch        *b_HLT_Mu7p5_Track7_Upsilon;   //!
   TBranch        *b_HLT_DoublePhoton33_CaloIdL;   //!
   TBranch        *b_HLT_DoublePhoton70;   //!
   TBranch        *b_HLT_DoublePhoton85;   //!
   TBranch        *b_HLT_Ele20_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele20_WPLoose_Gsf;   //!
   TBranch        *b_HLT_Ele20_eta2p1_WPLoose_Gsf;   //!
   TBranch        *b_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG;   //!
   TBranch        *b_HLT_Ele27_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele32_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele35_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele35_WPTight_Gsf_L1EGMT;   //!
   TBranch        *b_HLT_Ele38_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele40_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele32_WPTight_Gsf_L1DoubleEG;   //!
   TBranch        *b_HLT_HT450_Beamspot;   //!
   TBranch        *b_HLT_HT300_Beamspot;   //!
   TBranch        *b_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1;   //!
   TBranch        *b_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1;   //!
   TBranch        *b_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1;   //!
   TBranch        *b_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1;   //!
   TBranch        *b_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1;   //!
   TBranch        *b_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_TightID_SingleL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_TightID_SingleL1;   //!
   TBranch        *b_HLT_IsoMu20;   //!
   TBranch        *b_HLT_IsoMu24;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1;   //!
   TBranch        *b_HLT_IsoMu27;   //!
   TBranch        *b_HLT_IsoMu30;   //!
   TBranch        *b_HLT_UncorrectedJetE30_NoBPTX;   //!
   TBranch        *b_HLT_UncorrectedJetE30_NoBPTX3BX;   //!
   TBranch        *b_HLT_UncorrectedJetE60_NoBPTX3BX;   //!
   TBranch        *b_HLT_UncorrectedJetE70_NoBPTX3BX;   //!
   TBranch        *b_HLT_L1SingleMu18;   //!
   TBranch        *b_HLT_L1SingleMu25;   //!
   TBranch        *b_HLT_L2Mu10;   //!
   TBranch        *b_HLT_L2Mu10_NoVertex_NoBPTX3BX;   //!
   TBranch        *b_HLT_L2Mu10_NoVertex_NoBPTX;   //!
   TBranch        *b_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX;   //!
   TBranch        *b_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX;   //!
   TBranch        *b_HLT_L2Mu50;   //!
   TBranch        *b_HLT_DoubleL2Mu50;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;   //!
   TBranch        *b_HLT_Mu25_TkMu0_Onia;   //!
   TBranch        *b_HLT_Mu30_TkMu0_Onia;   //!
   TBranch        *b_HLT_Mu20_TkMu0_Phi;   //!
   TBranch        *b_HLT_Mu25_TkMu0_Phi;   //!
   TBranch        *b_HLT_Mu20;   //!
   TBranch        *b_HLT_Mu27;   //!
   TBranch        *b_HLT_Mu50;   //!
   TBranch        *b_HLT_Mu55;   //!
   TBranch        *b_HLT_OldMu100;   //!
   TBranch        *b_HLT_TkMu100;   //!
   TBranch        *b_HLT_DiPFJet15_NoCaloMatched;   //!
   TBranch        *b_HLT_DiPFJet25_NoCaloMatched;   //!
   TBranch        *b_HLT_DiPFJet15_FBEta3_NoCaloMatched;   //!
   TBranch        *b_HLT_DiPFJet25_FBEta3_NoCaloMatched;   //!
   TBranch        *b_HLT_DiPFJetAve40;   //!
   TBranch        *b_HLT_DiPFJetAve60;   //!
   TBranch        *b_HLT_DiPFJetAve80;   //!
   TBranch        *b_HLT_DiPFJetAve140;   //!
   TBranch        *b_HLT_DiPFJetAve200;   //!
   TBranch        *b_HLT_DiPFJetAve260;   //!
   TBranch        *b_HLT_DiPFJetAve320;   //!
   TBranch        *b_HLT_DiPFJetAve400;   //!
   TBranch        *b_HLT_DiPFJetAve500;   //!
   TBranch        *b_HLT_DiPFJetAve15_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve25_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve35_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve60_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve80_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve100_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve160_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve220_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve300_HFJEC;   //!
   TBranch        *b_HLT_AK8PFJet40;   //!
   TBranch        *b_HLT_AK8PFJet60;   //!
   TBranch        *b_HLT_AK8PFJet80;   //!
   TBranch        *b_HLT_AK8PFJet140;   //!
   TBranch        *b_HLT_AK8PFJet200;   //!
   TBranch        *b_HLT_AK8PFJet260;   //!
   TBranch        *b_HLT_AK8PFJet320;   //!
   TBranch        *b_HLT_AK8PFJet400;   //!
   TBranch        *b_HLT_AK8PFJet450;   //!
   TBranch        *b_HLT_AK8PFJet500;   //!
   TBranch        *b_HLT_AK8PFJet550;   //!
   TBranch        *b_HLT_PFJet40;   //!
   TBranch        *b_HLT_PFJet60;   //!
   TBranch        *b_HLT_PFJet80;   //!
   TBranch        *b_HLT_PFJet140;   //!
   TBranch        *b_HLT_PFJet200;   //!
   TBranch        *b_HLT_PFJet260;   //!
   TBranch        *b_HLT_PFJet320;   //!
   TBranch        *b_HLT_PFJet400;   //!
   TBranch        *b_HLT_PFJet450;   //!
   TBranch        *b_HLT_PFJet500;   //!
   TBranch        *b_HLT_PFJet550;   //!
   TBranch        *b_HLT_PFJetFwd40;   //!
   TBranch        *b_HLT_PFJetFwd60;   //!
   TBranch        *b_HLT_PFJetFwd80;   //!
   TBranch        *b_HLT_PFJetFwd140;   //!
   TBranch        *b_HLT_PFJetFwd200;   //!
   TBranch        *b_HLT_PFJetFwd260;   //!
   TBranch        *b_HLT_PFJetFwd320;   //!
   TBranch        *b_HLT_PFJetFwd400;   //!
   TBranch        *b_HLT_PFJetFwd450;   //!
   TBranch        *b_HLT_PFJetFwd500;   //!
   TBranch        *b_HLT_AK8PFJetFwd40;   //!
   TBranch        *b_HLT_AK8PFJetFwd60;   //!
   TBranch        *b_HLT_AK8PFJetFwd80;   //!
   TBranch        *b_HLT_AK8PFJetFwd140;   //!
   TBranch        *b_HLT_AK8PFJetFwd200;   //!
   TBranch        *b_HLT_AK8PFJetFwd260;   //!
   TBranch        *b_HLT_AK8PFJetFwd320;   //!
   TBranch        *b_HLT_AK8PFJetFwd400;   //!
   TBranch        *b_HLT_AK8PFJetFwd450;   //!
   TBranch        *b_HLT_AK8PFJetFwd500;   //!
   TBranch        *b_HLT_PFHT180;   //!
   TBranch        *b_HLT_PFHT250;   //!
   TBranch        *b_HLT_PFHT370;   //!
   TBranch        *b_HLT_PFHT430;   //!
   TBranch        *b_HLT_PFHT510;   //!
   TBranch        *b_HLT_PFHT590;   //!
   TBranch        *b_HLT_PFHT680;   //!
   TBranch        *b_HLT_PFHT780;   //!
   TBranch        *b_HLT_PFHT890;   //!
   TBranch        *b_HLT_PFHT800;   //!
   TBranch        *b_HLT_PFHT900;   //!
   TBranch        *b_HLT_PFHT1050;   //!
   TBranch        *b_HLT_PFHT500_PFMET100_PFMHT100_IDTight;   //!
   TBranch        *b_HLT_PFHT500_PFMET110_PFMHT110_IDTight;   //!
   TBranch        *b_HLT_PFHT700_PFMET85_PFMHT85_IDTight;   //!
   TBranch        *b_HLT_PFHT700_PFMET95_PFMHT95_IDTight;   //!
   TBranch        *b_HLT_PFHT800_PFMET75_PFMHT75_IDTight;   //!
   TBranch        *b_HLT_PFHT800_PFMET85_PFMHT85_IDTight;   //!
   TBranch        *b_HLT_PFMET110_PFMHT110_IDTight;   //!
   TBranch        *b_HLT_PFMET120_PFMHT120_IDTight;   //!
   TBranch        *b_HLT_PFMET130_PFMHT130_IDTight;   //!
   TBranch        *b_HLT_PFMET140_PFMHT140_IDTight;   //!
   TBranch        *b_HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1;   //!
   TBranch        *b_HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1;   //!
   TBranch        *b_HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1;   //!
   TBranch        *b_HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1;   //!
   TBranch        *b_HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1;   //!
   TBranch        *b_HLT_PFMET120_PFMHT120_IDTight_PFHT60;   //!
   TBranch        *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;   //!
   TBranch        *b_HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60;   //!
   TBranch        *b_HLT_PFMETTypeOne110_PFMHT110_IDTight;   //!
   TBranch        *b_HLT_PFMETTypeOne120_PFMHT120_IDTight;   //!
   TBranch        *b_HLT_PFMETTypeOne130_PFMHT130_IDTight;   //!
   TBranch        *b_HLT_PFMETTypeOne140_PFMHT140_IDTight;   //!
   TBranch        *b_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight;   //!
   TBranch        *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;   //!
   TBranch        *b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight;   //!
   TBranch        *b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight;   //!
   TBranch        *b_HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight;   //!
   TBranch        *b_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight;   //!
   TBranch        *b_HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight;   //!
   TBranch        *b_HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight;   //!
   TBranch        *b_HLT_L1ETMHadSeeds;   //!
   TBranch        *b_HLT_CaloMHT90;   //!
   TBranch        *b_HLT_CaloMET80_NotCleaned;   //!
   TBranch        *b_HLT_CaloMET90_NotCleaned;   //!
   TBranch        *b_HLT_CaloMET100_NotCleaned;   //!
   TBranch        *b_HLT_CaloMET110_NotCleaned;   //!
   TBranch        *b_HLT_CaloMET250_NotCleaned;   //!
   TBranch        *b_HLT_CaloMET70_HBHECleaned;   //!
   TBranch        *b_HLT_CaloMET80_HBHECleaned;   //!
   TBranch        *b_HLT_CaloMET90_HBHECleaned;   //!
   TBranch        *b_HLT_CaloMET100_HBHECleaned;   //!
   TBranch        *b_HLT_CaloMET250_HBHECleaned;   //!
   TBranch        *b_HLT_CaloMET300_HBHECleaned;   //!
   TBranch        *b_HLT_CaloMET350_HBHECleaned;   //!
   TBranch        *b_HLT_PFMET200_NotCleaned;   //!
   TBranch        *b_HLT_PFMET200_HBHECleaned;   //!
   TBranch        *b_HLT_PFMET250_HBHECleaned;   //!
   TBranch        *b_HLT_PFMET300_HBHECleaned;   //!
   TBranch        *b_HLT_PFMET200_HBHE_BeamHaloCleaned;   //!
   TBranch        *b_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned;   //!
   TBranch        *b_HLT_MET105_IsoTrk50;   //!
   TBranch        *b_HLT_MET120_IsoTrk50;   //!
   TBranch        *b_HLT_SingleJet30_Mu12_SinglePFJet40;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets40_CaloBTagCSV_p33;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets100_CaloBTagCSV_p33;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets200_CaloBTagCSV_p33;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets350_CaloBTagCSV_p33;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagCSV_p33;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagCSV_p33;   //!
   TBranch        *b_HLT_DoublePFJets40_CaloBTagCSV_p33;   //!
   TBranch        *b_HLT_DoublePFJets100_CaloBTagCSV_p33;   //!
   TBranch        *b_HLT_DoublePFJets200_CaloBTagCSV_p33;   //!
   TBranch        *b_HLT_DoublePFJets350_CaloBTagCSV_p33;   //!
   TBranch        *b_HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33;   //!
   TBranch        *b_HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33;   //!
   TBranch        *b_HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33;   //!
   TBranch        *b_HLT_Photon300_NoHE;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ;   //!
   TBranch        *b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL;   //!
   TBranch        *b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ;   //!
   TBranch        *b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL;   //!
   TBranch        *b_HLT_BTagMu_AK4DiJet20_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK4DiJet40_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK4DiJet70_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK4DiJet110_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK4DiJet170_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK4Jet300_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK8DiJet170_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK8Jet300_Mu5;   //!
   TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_HLT_Mu12_DoublePhoton20;   //!
   TBranch        *b_HLT_TriplePhoton_20_20_20_CaloIdLV2;   //!
   TBranch        *b_HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL;   //!
   TBranch        *b_HLT_TriplePhoton_30_30_10_CaloIdLV2;   //!
   TBranch        *b_HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL;   //!
   TBranch        *b_HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL;   //!
   TBranch        *b_HLT_Photon25;   //!
   TBranch        *b_HLT_Photon33;   //!
   TBranch        *b_HLT_Photon50;   //!
   TBranch        *b_HLT_Photon75;   //!
   TBranch        *b_HLT_Photon90;   //!
   TBranch        *b_HLT_Photon120;   //!
   TBranch        *b_HLT_Photon150;   //!
   TBranch        *b_HLT_Photon175;   //!
   TBranch        *b_HLT_Photon200;   //!
   TBranch        *b_HLT_Photon50_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon75_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon90_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon120_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon165_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon90_CaloIdL_PFHT700;   //!
   TBranch        *b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90;   //!
   TBranch        *b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95;   //!
   TBranch        *b_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55;   //!
   TBranch        *b_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55;   //!
   TBranch        *b_HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55;   //!
   TBranch        *b_HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_L1_NoOS;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_NoVertexing_NoOS;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_NoVertexing;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_L1_4R_0er1p5R;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi3p5_Muon2;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_L1_4p5;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_L1_5;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_L1_4p5NoOS;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_L1_4p5er2p0;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_L1_4p5er2p0M;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_NoVertexing;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_L1_5M;   //!
   TBranch        *b_HLT_Dimuon0_LowMass_L1_0er1p5R;   //!
   TBranch        *b_HLT_Dimuon0_LowMass_L1_0er1p5;   //!
   TBranch        *b_HLT_Dimuon0_LowMass;   //!
   TBranch        *b_HLT_Dimuon0_LowMass_L1_4;   //!
   TBranch        *b_HLT_Dimuon0_LowMass_L1_4R;   //!
   TBranch        *b_HLT_Dimuon0_LowMass_L1_TM530;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_Muon_L1_TM0;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_Muon_NoL1Mass;   //!
   TBranch        *b_HLT_TripleMu_5_3_3_Mass3p8to60_DZ;   //!
   TBranch        *b_HLT_TripleMu_10_5_5_DZ;   //!
   TBranch        *b_HLT_TripleMu_12_10_5;   //!
   TBranch        *b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15;   //!
   TBranch        *b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1;   //!
   TBranch        *b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15;   //!
   TBranch        *b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1;   //!
   TBranch        *b_HLT_DoubleMu3_DZ_PFMET50_PFMHT60;   //!
   TBranch        *b_HLT_DoubleMu3_DZ_PFMET70_PFMHT70;   //!
   TBranch        *b_HLT_DoubleMu3_DZ_PFMET90_PFMHT90;   //!
   TBranch        *b_HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass;   //!
   TBranch        *b_HLT_DoubleMu4_Jpsi_Displaced;   //!
   TBranch        *b_HLT_DoubleMu4_Jpsi_NoVertexing;   //!
   TBranch        *b_HLT_DoubleMu4_JpsiTrkTrk_Displaced;   //!
   TBranch        *b_HLT_DoubleMu43NoFiltersNoVtx;   //!
   TBranch        *b_HLT_DoubleMu48NoFiltersNoVtx;   //!
   TBranch        *b_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL;   //!
   TBranch        *b_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL;   //!
   TBranch        *b_HLT_DoubleMu20_7_Mass0to30_L1_DM4;   //!
   TBranch        *b_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG;   //!
   TBranch        *b_HLT_HT425;   //!
   TBranch        *b_HLT_HT430_DisplacedDijet40_DisplacedTrack;   //!
   TBranch        *b_HLT_HT430_DisplacedDijet60_DisplacedTrack;   //!
   TBranch        *b_HLT_HT430_DisplacedDijet80_DisplacedTrack;   //!
   TBranch        *b_HLT_HT400_DisplacedDijet40_DisplacedTrack;   //!
   TBranch        *b_HLT_HT650_DisplacedDijet60_Inclusive;   //!
   TBranch        *b_HLT_HT550_DisplacedDijet80_Inclusive;   //!
   TBranch        *b_HLT_HT550_DisplacedDijet60_Inclusive;   //!
   TBranch        *b_HLT_HT650_DisplacedDijet80_Inclusive;   //!
   TBranch        *b_HLT_HT750_DisplacedDijet80_Inclusive;   //!
   TBranch        *b_HLT_DiJet110_35_Mjj650_PFMET110;   //!
   TBranch        *b_HLT_DiJet110_35_Mjj650_PFMET120;   //!
   TBranch        *b_HLT_DiJet110_35_Mjj650_PFMET130;   //!
   TBranch        *b_HLT_TripleJet110_35_35_Mjj650_PFMET110;   //!
   TBranch        *b_HLT_TripleJet110_35_35_Mjj650_PFMET120;   //!
   TBranch        *b_HLT_TripleJet110_35_35_Mjj650_PFMET130;   //!
   TBranch        *b_HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg;   //!
   TBranch        *b_HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg;   //!
   TBranch        *b_HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg;   //!
   TBranch        *b_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned;   //!
   TBranch        *b_HLT_Ele28_eta2p1_WPTight_Gsf_HT150;   //!
   TBranch        *b_HLT_Ele28_HighEta_SC20_Mass55;   //!
   TBranch        *b_HLT_DoubleMu20_7_Mass0to30_Photon23;   //!
   TBranch        *b_HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5;   //!
   TBranch        *b_HLT_Ele15_IsoVVVL_PFHT450_PFMET50;   //!
   TBranch        *b_HLT_Ele15_IsoVVVL_PFHT450;   //!
   TBranch        *b_HLT_Ele50_IsoVVVL_PFHT450;   //!
   TBranch        *b_HLT_Ele15_IsoVVVL_PFHT600;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60;   //!
   TBranch        *b_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60;   //!
   TBranch        *b_HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5;   //!
   TBranch        *b_HLT_Mu15_IsoVVVL_PFHT450_PFMET50;   //!
   TBranch        *b_HLT_Mu15_IsoVVVL_PFHT450;   //!
   TBranch        *b_HLT_Mu50_IsoVVVL_PFHT450;   //!
   TBranch        *b_HLT_Mu15_IsoVVVL_PFHT600;   //!
   TBranch        *b_HLT_Dimuon10_PsiPrime_Barrel_Seagulls;   //!
   TBranch        *b_HLT_Dimuon20_Jpsi_Barrel_Seagulls;   //!
   TBranch        *b_HLT_Dimuon10_Upsilon_Barrel_Seagulls;   //!
   TBranch        *b_HLT_Dimuon12_Upsilon_eta1p5;   //!
   TBranch        *b_HLT_Dimuon14_Phi_Barrel_Seagulls;   //!
   TBranch        *b_HLT_Dimuon18_PsiPrime;   //!
   TBranch        *b_HLT_Dimuon25_Jpsi;   //!
   TBranch        *b_HLT_Dimuon18_PsiPrime_noCorrL1;   //!
   TBranch        *b_HLT_Dimuon24_Upsilon_noCorrL1;   //!
   TBranch        *b_HLT_Dimuon24_Phi_noCorrL1;   //!
   TBranch        *b_HLT_Dimuon25_Jpsi_noCorrL1;   //!
   TBranch        *b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;   //!
   TBranch        *b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL;   //!
   TBranch        *b_HLT_DoubleIsoMu20_eta2p1;   //!
   TBranch        *b_HLT_DoubleIsoMu24_eta2p1;   //!
   TBranch        *b_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx;   //!
   TBranch        *b_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx;   //!
   TBranch        *b_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx;   //!
   TBranch        *b_HLT_Mu8;   //!
   TBranch        *b_HLT_Mu17;   //!
   TBranch        *b_HLT_Mu19;   //!
   TBranch        *b_HLT_Mu17_Photon30_IsoCaloId;   //!
   TBranch        *b_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
   TBranch        *b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
   TBranch        *b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
   TBranch        *b_HLT_Ele8_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_HLT_Ele17_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_HLT_Ele23_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;   //!
   TBranch        *b_HLT_Ele115_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_HLT_Ele135_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_HLT_Ele145_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_HLT_Ele200_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_HLT_Ele250_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_HLT_Ele300_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_HLT_PFHT300PT30_QuadPFJet_75_60_45_40;   //!
   TBranch        *b_HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0;   //!
   TBranch        *b_HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2;   //!
   TBranch        *b_HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2;   //!
   TBranch        *b_HLT_PFHT380_SixPFJet32;   //!
   TBranch        *b_HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5;   //!
   TBranch        *b_HLT_PFHT430_SixPFJet40;   //!
   TBranch        *b_HLT_PFHT350;   //!
   TBranch        *b_HLT_PFHT350MinPFJet15;   //!
   TBranch        *b_HLT_Photon60_R9Id90_CaloIdL_IsoL;   //!
   TBranch        *b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL;   //!
   TBranch        *b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15;   //!
   TBranch        *b_HLT_FullTrack_Multiplicity85;   //!
   TBranch        *b_HLT_FullTrack_Multiplicity100;   //!
   TBranch        *b_HLT_FullTrack_Multiplicity130;   //!
   TBranch        *b_HLT_FullTrack_Multiplicity155;   //!
   TBranch        *b_HLT_ECALHT800;   //!
   TBranch        *b_HLT_DiSC30_18_EIso_AND_HE_Mass70;   //!
   TBranch        *b_HLT_Physics;   //!
   TBranch        *b_HLT_Physics_part0;   //!
   TBranch        *b_HLT_Physics_part1;   //!
   TBranch        *b_HLT_Physics_part2;   //!
   TBranch        *b_HLT_Physics_part3;   //!
   TBranch        *b_HLT_Physics_part4;   //!
   TBranch        *b_HLT_Physics_part5;   //!
   TBranch        *b_HLT_Physics_part6;   //!
   TBranch        *b_HLT_Physics_part7;   //!
   TBranch        *b_HLT_Random;   //!
   TBranch        *b_HLT_ZeroBias;   //!
   TBranch        *b_HLT_ZeroBias_part0;   //!
   TBranch        *b_HLT_ZeroBias_part1;   //!
   TBranch        *b_HLT_ZeroBias_part2;   //!
   TBranch        *b_HLT_ZeroBias_part3;   //!
   TBranch        *b_HLT_ZeroBias_part4;   //!
   TBranch        *b_HLT_ZeroBias_part5;   //!
   TBranch        *b_HLT_ZeroBias_part6;   //!
   TBranch        *b_HLT_ZeroBias_part7;   //!
   TBranch        *b_HLT_AK4CaloJet30;   //!
   TBranch        *b_HLT_AK4CaloJet40;   //!
   TBranch        *b_HLT_AK4CaloJet50;   //!
   TBranch        *b_HLT_AK4CaloJet80;   //!
   TBranch        *b_HLT_AK4CaloJet100;   //!
   TBranch        *b_HLT_AK4CaloJet120;   //!
   TBranch        *b_HLT_AK4PFJet30;   //!
   TBranch        *b_HLT_AK4PFJet50;   //!
   TBranch        *b_HLT_AK4PFJet80;   //!
   TBranch        *b_HLT_AK4PFJet100;   //!
   TBranch        *b_HLT_AK4PFJet120;   //!
   TBranch        *b_HLT_HISinglePhoton10_Eta3p1ForPPRef;   //!
   TBranch        *b_HLT_HISinglePhoton20_Eta3p1ForPPRef;   //!
   TBranch        *b_HLT_HISinglePhoton30_Eta3p1ForPPRef;   //!
   TBranch        *b_HLT_HISinglePhoton40_Eta3p1ForPPRef;   //!
   TBranch        *b_HLT_HISinglePhoton50_Eta3p1ForPPRef;   //!
   TBranch        *b_HLT_HISinglePhoton60_Eta3p1ForPPRef;   //!
   TBranch        *b_HLT_Photon20_HoverELoose;   //!
   TBranch        *b_HLT_Photon30_HoverELoose;   //!
   TBranch        *b_HLT_Photon40_HoverELoose;   //!
   TBranch        *b_HLT_Photon50_HoverELoose;   //!
   TBranch        *b_HLT_Photon60_HoverELoose;   //!
   TBranch        *b_HLT_EcalCalibration;   //!
   TBranch        *b_HLT_HcalCalibration;   //!
   TBranch        *b_HLT_L1UnpairedBunchBptxMinus;   //!
   TBranch        *b_HLT_L1UnpairedBunchBptxPlus;   //!
   TBranch        *b_HLT_L1NotBptxOR;   //!
   TBranch        *b_HLT_L1MinimumBiasHF_OR;   //!
   TBranch        *b_HLT_L1MinimumBiasHF0OR;   //!
   TBranch        *b_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;   //!
   TBranch        *b_HLT_HcalNZS;   //!
   TBranch        *b_HLT_HcalPhiSym;   //!
   TBranch        *b_HLT_HcalIsolatedbunch;   //!
   TBranch        *b_HLT_IsoTrackHB;   //!
   TBranch        *b_HLT_IsoTrackHE;   //!
   TBranch        *b_HLT_ZeroBias_FirstCollisionAfterAbortGap;   //!
   TBranch        *b_HLT_ZeroBias_IsolatedBunches;   //!
   TBranch        *b_HLT_ZeroBias_FirstCollisionInTrain;   //!
   TBranch        *b_HLT_ZeroBias_LastCollisionInTrain;   //!
   TBranch        *b_HLT_ZeroBias_FirstBXAfterTrain;   //!
   TBranch        *b_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1;   //!
   TBranch        *b_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1;   //!
   TBranch        *b_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1;   //!
   TBranch        *b_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1;   //!
   TBranch        *b_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1;   //!
   TBranch        *b_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1;   //!
   TBranch        *b_HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr;   //!
   TBranch        *b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90;   //!
   TBranch        *b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100;   //!
   TBranch        *b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110;   //!
   TBranch        *b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120;   //!
   TBranch        *b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130;   //!
   TBranch        *b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr;   //!
   TBranch        *b_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr;   //!
   TBranch        *b_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;   //!
   TBranch        *b_HLT_Rsq0p35;   //!
   TBranch        *b_HLT_Rsq0p40;   //!
   TBranch        *b_HLT_RsqMR300_Rsq0p09_MR200;   //!
   TBranch        *b_HLT_RsqMR320_Rsq0p09_MR200;   //!
   TBranch        *b_HLT_RsqMR300_Rsq0p09_MR200_4jet;   //!
   TBranch        *b_HLT_RsqMR320_Rsq0p09_MR200_4jet;   //!
   TBranch        *b_HLT_L1_DoubleJet30_Mass_Min400_Mu10;   //!
   TBranch        *b_HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1;   //!
   TBranch        *b_HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1;   //!
   TBranch        *b_HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1;   //!
   TBranch        *b_HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50;   //!
   TBranch        *b_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3;   //!
   TBranch        *b_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3;   //!
   TBranch        *b_HLT_PFMET100_PFMHT100_IDTight_PFHT60;   //!
   TBranch        *b_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60;   //!
   TBranch        *b_HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60;   //!
   TBranch        *b_HLT_Mu18_Mu9_SameSign;   //!
   TBranch        *b_HLT_Mu18_Mu9_SameSign_DZ;   //!
   TBranch        *b_HLT_Mu18_Mu9;   //!
   TBranch        *b_HLT_Mu18_Mu9_DZ;   //!
   TBranch        *b_HLT_Mu20_Mu10_SameSign;   //!
   TBranch        *b_HLT_Mu20_Mu10_SameSign_DZ;   //!
   TBranch        *b_HLT_Mu20_Mu10;   //!
   TBranch        *b_HLT_Mu20_Mu10_DZ;   //!
   TBranch        *b_HLT_Mu23_Mu12_SameSign;   //!
   TBranch        *b_HLT_Mu23_Mu12_SameSign_DZ;   //!
   TBranch        *b_HLT_Mu23_Mu12;   //!
   TBranch        *b_HLT_Mu23_Mu12_DZ;   //!
   TBranch        *b_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi;   //!
   TBranch        *b_HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi;   //!
   TBranch        *b_HLT_DoubleMu3_DCA_PFMET50_PFMHT60;   //!
   TBranch        *b_HLT_TripleMu_5_3_3_Mass3p8to60_DCA;   //!
   TBranch        *b_HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1;   //!
   TBranch        *b_HLT_QuadPFJet103_88_75_15_DoubleBTagCSV_p013_p08_VBF1;   //!
   TBranch        *b_HLT_QuadPFJet105_90_76_15_DoubleBTagCSV_p013_p08_VBF1;   //!
   TBranch        *b_HLT_QuadPFJet111_90_80_15_DoubleBTagCSV_p013_p08_VBF1;   //!
   TBranch        *b_HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2;   //!
   TBranch        *b_HLT_QuadPFJet103_88_75_15_BTagCSV_p013_VBF2;   //!
   TBranch        *b_HLT_QuadPFJet105_88_76_15_BTagCSV_p013_VBF2;   //!
   TBranch        *b_HLT_QuadPFJet111_90_80_15_BTagCSV_p013_VBF2;   //!
   TBranch        *b_HLT_QuadPFJet98_83_71_15;   //!
   TBranch        *b_HLT_QuadPFJet103_88_75_15;   //!
   TBranch        *b_HLT_QuadPFJet105_88_76_15;   //!
   TBranch        *b_HLT_QuadPFJet111_90_80_15;   //!
   TBranch        *b_HLT_AK8PFJet330_PFAK8BTagCSV_p17;   //!
   TBranch        *b_HLT_AK8PFJet330_PFAK8BTagCSV_p1;   //!
   TBranch        *b_HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55;   //!
   TBranch        *b_HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55;   //!
   TBranch        *b_HLTriggerFinalPath;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_Flag_CSCTightHaloFilter;   //!
   TBranch        *b_Flag_CSCTightHaloTrkMuUnvetoFilter;   //!
   TBranch        *b_Flag_CSCTightHalo2015Filter;   //!
   TBranch        *b_Flag_globalTightHalo2016Filter;   //!
   TBranch        *b_Flag_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_Flag_HcalStripHaloFilter;   //!
   TBranch        *b_Flag_hcalLaserEventFilter;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_EcalDeadCellBoundaryEnergyFilter;   //!
   TBranch        *b_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_ecalLaserCorrFilter;   //!
   TBranch        *b_Flag_trkPOGFilters;   //!
   TBranch        *b_Flag_chargedHadronTrackResolutionFilter;   //!
   TBranch        *b_Flag_muonBadTrackFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag_BadPFMuonFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateSummer16Filter;   //!
   TBranch        *b_Flag_BadPFMuonSummer16Filter;   //!
   TBranch        *b_Flag_trkPOG_manystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_toomanystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_logErrorTooManyClusters;   //!
   TBranch        *b_Flag_METFilters;   //!
   TBranch        *b_Jet_pt_raw;   //!
   TBranch        *b_Jet_pt_nom;   //!
   TBranch        *b_Jet_corr_JEC;   //!
   TBranch        *b_Jet_corr_JER;   //!
   TBranch        *b_Jet_mass_raw;   //!
   TBranch        *b_Jet_mass_nom;   //!
   TBranch        *b_MET_pt_nom;   //!
   TBranch        *b_MET_phi_nom;   //!
   TBranch        *b_Jet_pt_jerUp;   //!
   TBranch        *b_Jet_mass_jerUp;   //!
   TBranch        *b_Jet_mass_jmrUp;   //!
   TBranch        *b_Jet_mass_jmsUp;   //!
   TBranch        *b_MET_pt_jerUp;   //!
   TBranch        *b_MET_phi_jerUp;   //!
   TBranch        *b_Jet_pt_jesAbsoluteStatUp;   //!
   TBranch        *b_Jet_mass_jesAbsoluteStatUp;   //!
   TBranch        *b_MET_pt_jesAbsoluteStatUp;   //!
   TBranch        *b_MET_phi_jesAbsoluteStatUp;   //!
   TBranch        *b_Jet_pt_jesAbsoluteScaleUp;   //!
   TBranch        *b_Jet_mass_jesAbsoluteScaleUp;   //!
   TBranch        *b_MET_pt_jesAbsoluteScaleUp;   //!
   TBranch        *b_MET_phi_jesAbsoluteScaleUp;   //!
   TBranch        *b_Jet_pt_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_Jet_mass_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_MET_pt_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_MET_phi_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_Jet_pt_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_Jet_mass_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_MET_pt_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_MET_phi_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_Jet_pt_jesFragmentationUp;   //!
   TBranch        *b_Jet_mass_jesFragmentationUp;   //!
   TBranch        *b_MET_pt_jesFragmentationUp;   //!
   TBranch        *b_MET_phi_jesFragmentationUp;   //!
   TBranch        *b_Jet_pt_jesSinglePionECALUp;   //!
   TBranch        *b_Jet_mass_jesSinglePionECALUp;   //!
   TBranch        *b_MET_pt_jesSinglePionECALUp;   //!
   TBranch        *b_MET_phi_jesSinglePionECALUp;   //!
   TBranch        *b_Jet_pt_jesSinglePionHCALUp;   //!
   TBranch        *b_Jet_mass_jesSinglePionHCALUp;   //!
   TBranch        *b_MET_pt_jesSinglePionHCALUp;   //!
   TBranch        *b_MET_phi_jesSinglePionHCALUp;   //!
   TBranch        *b_Jet_pt_jesFlavorQCDUp;   //!
   TBranch        *b_Jet_mass_jesFlavorQCDUp;   //!
   TBranch        *b_MET_pt_jesFlavorQCDUp;   //!
   TBranch        *b_MET_phi_jesFlavorQCDUp;   //!
   TBranch        *b_Jet_pt_jesTimePtEtaUp;   //!
   TBranch        *b_Jet_mass_jesTimePtEtaUp;   //!
   TBranch        *b_MET_pt_jesTimePtEtaUp;   //!
   TBranch        *b_MET_phi_jesTimePtEtaUp;   //!
   TBranch        *b_Jet_pt_jesRelativeJEREC1Up;   //!
   TBranch        *b_Jet_mass_jesRelativeJEREC1Up;   //!
   TBranch        *b_MET_pt_jesRelativeJEREC1Up;   //!
   TBranch        *b_MET_phi_jesRelativeJEREC1Up;   //!
   TBranch        *b_Jet_pt_jesRelativeJEREC2Up;   //!
   TBranch        *b_Jet_mass_jesRelativeJEREC2Up;   //!
   TBranch        *b_MET_pt_jesRelativeJEREC2Up;   //!
   TBranch        *b_MET_phi_jesRelativeJEREC2Up;   //!
   TBranch        *b_Jet_pt_jesRelativeJERHFUp;   //!
   TBranch        *b_Jet_mass_jesRelativeJERHFUp;   //!
   TBranch        *b_MET_pt_jesRelativeJERHFUp;   //!
   TBranch        *b_MET_phi_jesRelativeJERHFUp;   //!
   TBranch        *b_Jet_pt_jesRelativePtBBUp;   //!
   TBranch        *b_Jet_mass_jesRelativePtBBUp;   //!
   TBranch        *b_MET_pt_jesRelativePtBBUp;   //!
   TBranch        *b_MET_phi_jesRelativePtBBUp;   //!
   TBranch        *b_Jet_pt_jesRelativePtEC1Up;   //!
   TBranch        *b_Jet_mass_jesRelativePtEC1Up;   //!
   TBranch        *b_MET_pt_jesRelativePtEC1Up;   //!
   TBranch        *b_MET_phi_jesRelativePtEC1Up;   //!
   TBranch        *b_Jet_pt_jesRelativePtEC2Up;   //!
   TBranch        *b_Jet_mass_jesRelativePtEC2Up;   //!
   TBranch        *b_MET_pt_jesRelativePtEC2Up;   //!
   TBranch        *b_MET_phi_jesRelativePtEC2Up;   //!
   TBranch        *b_Jet_pt_jesRelativePtHFUp;   //!
   TBranch        *b_Jet_mass_jesRelativePtHFUp;   //!
   TBranch        *b_MET_pt_jesRelativePtHFUp;   //!
   TBranch        *b_MET_phi_jesRelativePtHFUp;   //!
   TBranch        *b_Jet_pt_jesRelativeBalUp;   //!
   TBranch        *b_Jet_mass_jesRelativeBalUp;   //!
   TBranch        *b_MET_pt_jesRelativeBalUp;   //!
   TBranch        *b_MET_phi_jesRelativeBalUp;   //!
   TBranch        *b_Jet_pt_jesRelativeSampleUp;   //!
   TBranch        *b_Jet_mass_jesRelativeSampleUp;   //!
   TBranch        *b_MET_pt_jesRelativeSampleUp;   //!
   TBranch        *b_MET_phi_jesRelativeSampleUp;   //!
   TBranch        *b_Jet_pt_jesRelativeFSRUp;   //!
   TBranch        *b_Jet_mass_jesRelativeFSRUp;   //!
   TBranch        *b_MET_pt_jesRelativeFSRUp;   //!
   TBranch        *b_MET_phi_jesRelativeFSRUp;   //!
   TBranch        *b_Jet_pt_jesRelativeStatFSRUp;   //!
   TBranch        *b_Jet_mass_jesRelativeStatFSRUp;   //!
   TBranch        *b_MET_pt_jesRelativeStatFSRUp;   //!
   TBranch        *b_MET_phi_jesRelativeStatFSRUp;   //!
   TBranch        *b_Jet_pt_jesRelativeStatECUp;   //!
   TBranch        *b_Jet_mass_jesRelativeStatECUp;   //!
   TBranch        *b_MET_pt_jesRelativeStatECUp;   //!
   TBranch        *b_MET_phi_jesRelativeStatECUp;   //!
   TBranch        *b_Jet_pt_jesRelativeStatHFUp;   //!
   TBranch        *b_Jet_mass_jesRelativeStatHFUp;   //!
   TBranch        *b_MET_pt_jesRelativeStatHFUp;   //!
   TBranch        *b_MET_phi_jesRelativeStatHFUp;   //!
   TBranch        *b_Jet_pt_jesPileUpDataMCUp;   //!
   TBranch        *b_Jet_mass_jesPileUpDataMCUp;   //!
   TBranch        *b_MET_pt_jesPileUpDataMCUp;   //!
   TBranch        *b_MET_phi_jesPileUpDataMCUp;   //!
   TBranch        *b_Jet_pt_jesPileUpPtRefUp;   //!
   TBranch        *b_Jet_mass_jesPileUpPtRefUp;   //!
   TBranch        *b_MET_pt_jesPileUpPtRefUp;   //!
   TBranch        *b_MET_phi_jesPileUpPtRefUp;   //!
   TBranch        *b_Jet_pt_jesPileUpPtBBUp;   //!
   TBranch        *b_Jet_mass_jesPileUpPtBBUp;   //!
   TBranch        *b_MET_pt_jesPileUpPtBBUp;   //!
   TBranch        *b_MET_phi_jesPileUpPtBBUp;   //!
   TBranch        *b_Jet_pt_jesPileUpPtEC1Up;   //!
   TBranch        *b_Jet_mass_jesPileUpPtEC1Up;   //!
   TBranch        *b_MET_pt_jesPileUpPtEC1Up;   //!
   TBranch        *b_MET_phi_jesPileUpPtEC1Up;   //!
   TBranch        *b_Jet_pt_jesPileUpPtEC2Up;   //!
   TBranch        *b_Jet_mass_jesPileUpPtEC2Up;   //!
   TBranch        *b_MET_pt_jesPileUpPtEC2Up;   //!
   TBranch        *b_MET_phi_jesPileUpPtEC2Up;   //!
   TBranch        *b_Jet_pt_jesPileUpPtHFUp;   //!
   TBranch        *b_Jet_mass_jesPileUpPtHFUp;   //!
   TBranch        *b_MET_pt_jesPileUpPtHFUp;   //!
   TBranch        *b_MET_phi_jesPileUpPtHFUp;   //!
   TBranch        *b_Jet_pt_jesPileUpMuZeroUp;   //!
   TBranch        *b_Jet_mass_jesPileUpMuZeroUp;   //!
   TBranch        *b_MET_pt_jesPileUpMuZeroUp;   //!
   TBranch        *b_MET_phi_jesPileUpMuZeroUp;   //!
   TBranch        *b_Jet_pt_jesPileUpEnvelopeUp;   //!
   TBranch        *b_Jet_mass_jesPileUpEnvelopeUp;   //!
   TBranch        *b_MET_pt_jesPileUpEnvelopeUp;   //!
   TBranch        *b_MET_phi_jesPileUpEnvelopeUp;   //!
   TBranch        *b_Jet_pt_jesSubTotalPileUpUp;   //!
   TBranch        *b_Jet_mass_jesSubTotalPileUpUp;   //!
   TBranch        *b_MET_pt_jesSubTotalPileUpUp;   //!
   TBranch        *b_MET_phi_jesSubTotalPileUpUp;   //!
   TBranch        *b_Jet_pt_jesSubTotalRelativeUp;   //!
   TBranch        *b_Jet_mass_jesSubTotalRelativeUp;   //!
   TBranch        *b_MET_pt_jesSubTotalRelativeUp;   //!
   TBranch        *b_MET_phi_jesSubTotalRelativeUp;   //!
   TBranch        *b_Jet_pt_jesSubTotalPtUp;   //!
   TBranch        *b_Jet_mass_jesSubTotalPtUp;   //!
   TBranch        *b_MET_pt_jesSubTotalPtUp;   //!
   TBranch        *b_MET_phi_jesSubTotalPtUp;   //!
   TBranch        *b_Jet_pt_jesSubTotalScaleUp;   //!
   TBranch        *b_Jet_mass_jesSubTotalScaleUp;   //!
   TBranch        *b_MET_pt_jesSubTotalScaleUp;   //!
   TBranch        *b_MET_phi_jesSubTotalScaleUp;   //!
   TBranch        *b_Jet_pt_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_Jet_mass_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_MET_pt_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_MET_phi_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_Jet_pt_jesSubTotalMCUp;   //!
   TBranch        *b_Jet_mass_jesSubTotalMCUp;   //!
   TBranch        *b_MET_pt_jesSubTotalMCUp;   //!
   TBranch        *b_MET_phi_jesSubTotalMCUp;   //!
   TBranch        *b_Jet_pt_jesTotalUp;   //!
   TBranch        *b_Jet_mass_jesTotalUp;   //!
   TBranch        *b_MET_pt_jesTotalUp;   //!
   TBranch        *b_MET_phi_jesTotalUp;   //!
   TBranch        *b_Jet_pt_jesTotalNoFlavorUp;   //!
   TBranch        *b_Jet_mass_jesTotalNoFlavorUp;   //!
   TBranch        *b_MET_pt_jesTotalNoFlavorUp;   //!
   TBranch        *b_MET_phi_jesTotalNoFlavorUp;   //!
   TBranch        *b_Jet_pt_jesTotalNoTimeUp;   //!
   TBranch        *b_Jet_mass_jesTotalNoTimeUp;   //!
   TBranch        *b_MET_pt_jesTotalNoTimeUp;   //!
   TBranch        *b_MET_phi_jesTotalNoTimeUp;   //!
   TBranch        *b_Jet_pt_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_Jet_mass_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_MET_pt_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_MET_phi_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_Jet_pt_jesFlavorZJetUp;   //!
   TBranch        *b_Jet_mass_jesFlavorZJetUp;   //!
   TBranch        *b_MET_pt_jesFlavorZJetUp;   //!
   TBranch        *b_MET_phi_jesFlavorZJetUp;   //!
   TBranch        *b_Jet_pt_jesFlavorPhotonJetUp;   //!
   TBranch        *b_Jet_mass_jesFlavorPhotonJetUp;   //!
   TBranch        *b_MET_pt_jesFlavorPhotonJetUp;   //!
   TBranch        *b_MET_phi_jesFlavorPhotonJetUp;   //!
   TBranch        *b_Jet_pt_jesFlavorPureGluonUp;   //!
   TBranch        *b_Jet_mass_jesFlavorPureGluonUp;   //!
   TBranch        *b_MET_pt_jesFlavorPureGluonUp;   //!
   TBranch        *b_MET_phi_jesFlavorPureGluonUp;   //!
   TBranch        *b_Jet_pt_jesFlavorPureQuarkUp;   //!
   TBranch        *b_Jet_mass_jesFlavorPureQuarkUp;   //!
   TBranch        *b_MET_pt_jesFlavorPureQuarkUp;   //!
   TBranch        *b_MET_phi_jesFlavorPureQuarkUp;   //!
   TBranch        *b_Jet_pt_jesFlavorPureCharmUp;   //!
   TBranch        *b_Jet_mass_jesFlavorPureCharmUp;   //!
   TBranch        *b_MET_pt_jesFlavorPureCharmUp;   //!
   TBranch        *b_MET_phi_jesFlavorPureCharmUp;   //!
   TBranch        *b_Jet_pt_jesFlavorPureBottomUp;   //!
   TBranch        *b_Jet_mass_jesFlavorPureBottomUp;   //!
   TBranch        *b_MET_pt_jesFlavorPureBottomUp;   //!
   TBranch        *b_MET_phi_jesFlavorPureBottomUp;   //!
   TBranch        *b_Jet_pt_jesTimeRunBUp;   //!
   TBranch        *b_Jet_mass_jesTimeRunBUp;   //!
   TBranch        *b_MET_pt_jesTimeRunBUp;   //!
   TBranch        *b_MET_phi_jesTimeRunBUp;   //!
   TBranch        *b_Jet_pt_jesTimeRunCUp;   //!
   TBranch        *b_Jet_mass_jesTimeRunCUp;   //!
   TBranch        *b_MET_pt_jesTimeRunCUp;   //!
   TBranch        *b_MET_phi_jesTimeRunCUp;   //!
   TBranch        *b_Jet_pt_jesTimeRunDEUp;   //!
   TBranch        *b_Jet_mass_jesTimeRunDEUp;   //!
   TBranch        *b_MET_pt_jesTimeRunDEUp;   //!
   TBranch        *b_MET_phi_jesTimeRunDEUp;   //!
   TBranch        *b_Jet_pt_jesTimeRunFUp;   //!
   TBranch        *b_Jet_mass_jesTimeRunFUp;   //!
   TBranch        *b_MET_pt_jesTimeRunFUp;   //!
   TBranch        *b_MET_phi_jesTimeRunFUp;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_MET_pt_unclustEnUp;   //!
   TBranch        *b_MET_phi_unclustEnUp;   //!
   TBranch        *b_Jet_pt_jerDown;   //!
   TBranch        *b_Jet_mass_jerDown;   //!
   TBranch        *b_Jet_mass_jmrDown;   //!
   TBranch        *b_Jet_mass_jmsDown;   //!
   TBranch        *b_MET_pt_jerDown;   //!
   TBranch        *b_MET_phi_jerDown;   //!
   TBranch        *b_Jet_pt_jesAbsoluteStatDown;   //!
   TBranch        *b_Jet_mass_jesAbsoluteStatDown;   //!
   TBranch        *b_MET_pt_jesAbsoluteStatDown;   //!
   TBranch        *b_MET_phi_jesAbsoluteStatDown;   //!
   TBranch        *b_Jet_pt_jesAbsoluteScaleDown;   //!
   TBranch        *b_Jet_mass_jesAbsoluteScaleDown;   //!
   TBranch        *b_MET_pt_jesAbsoluteScaleDown;   //!
   TBranch        *b_MET_phi_jesAbsoluteScaleDown;   //!
   TBranch        *b_Jet_pt_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_Jet_mass_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_MET_pt_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_MET_phi_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_Jet_pt_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_Jet_mass_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_MET_pt_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_MET_phi_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_Jet_pt_jesFragmentationDown;   //!
   TBranch        *b_Jet_mass_jesFragmentationDown;   //!
   TBranch        *b_MET_pt_jesFragmentationDown;   //!
   TBranch        *b_MET_phi_jesFragmentationDown;   //!
   TBranch        *b_Jet_pt_jesSinglePionECALDown;   //!
   TBranch        *b_Jet_mass_jesSinglePionECALDown;   //!
   TBranch        *b_MET_pt_jesSinglePionECALDown;   //!
   TBranch        *b_MET_phi_jesSinglePionECALDown;   //!
   TBranch        *b_Jet_pt_jesSinglePionHCALDown;   //!
   TBranch        *b_Jet_mass_jesSinglePionHCALDown;   //!
   TBranch        *b_MET_pt_jesSinglePionHCALDown;   //!
   TBranch        *b_MET_phi_jesSinglePionHCALDown;   //!
   TBranch        *b_Jet_pt_jesFlavorQCDDown;   //!
   TBranch        *b_Jet_mass_jesFlavorQCDDown;   //!
   TBranch        *b_MET_pt_jesFlavorQCDDown;   //!
   TBranch        *b_MET_phi_jesFlavorQCDDown;   //!
   TBranch        *b_Jet_pt_jesTimePtEtaDown;   //!
   TBranch        *b_Jet_mass_jesTimePtEtaDown;   //!
   TBranch        *b_MET_pt_jesTimePtEtaDown;   //!
   TBranch        *b_MET_phi_jesTimePtEtaDown;   //!
   TBranch        *b_Jet_pt_jesRelativeJEREC1Down;   //!
   TBranch        *b_Jet_mass_jesRelativeJEREC1Down;   //!
   TBranch        *b_MET_pt_jesRelativeJEREC1Down;   //!
   TBranch        *b_MET_phi_jesRelativeJEREC1Down;   //!
   TBranch        *b_Jet_pt_jesRelativeJEREC2Down;   //!
   TBranch        *b_Jet_mass_jesRelativeJEREC2Down;   //!
   TBranch        *b_MET_pt_jesRelativeJEREC2Down;   //!
   TBranch        *b_MET_phi_jesRelativeJEREC2Down;   //!
   TBranch        *b_Jet_pt_jesRelativeJERHFDown;   //!
   TBranch        *b_Jet_mass_jesRelativeJERHFDown;   //!
   TBranch        *b_MET_pt_jesRelativeJERHFDown;   //!
   TBranch        *b_MET_phi_jesRelativeJERHFDown;   //!
   TBranch        *b_Jet_pt_jesRelativePtBBDown;   //!
   TBranch        *b_Jet_mass_jesRelativePtBBDown;   //!
   TBranch        *b_MET_pt_jesRelativePtBBDown;   //!
   TBranch        *b_MET_phi_jesRelativePtBBDown;   //!
   TBranch        *b_Jet_pt_jesRelativePtEC1Down;   //!
   TBranch        *b_Jet_mass_jesRelativePtEC1Down;   //!
   TBranch        *b_MET_pt_jesRelativePtEC1Down;   //!
   TBranch        *b_MET_phi_jesRelativePtEC1Down;   //!
   TBranch        *b_Jet_pt_jesRelativePtEC2Down;   //!
   TBranch        *b_Jet_mass_jesRelativePtEC2Down;   //!
   TBranch        *b_MET_pt_jesRelativePtEC2Down;   //!
   TBranch        *b_MET_phi_jesRelativePtEC2Down;   //!
   TBranch        *b_Jet_pt_jesRelativePtHFDown;   //!
   TBranch        *b_Jet_mass_jesRelativePtHFDown;   //!
   TBranch        *b_MET_pt_jesRelativePtHFDown;   //!
   TBranch        *b_MET_phi_jesRelativePtHFDown;   //!
   TBranch        *b_Jet_pt_jesRelativeBalDown;   //!
   TBranch        *b_Jet_mass_jesRelativeBalDown;   //!
   TBranch        *b_MET_pt_jesRelativeBalDown;   //!
   TBranch        *b_MET_phi_jesRelativeBalDown;   //!
   TBranch        *b_Jet_pt_jesRelativeSampleDown;   //!
   TBranch        *b_Jet_mass_jesRelativeSampleDown;   //!
   TBranch        *b_MET_pt_jesRelativeSampleDown;   //!
   TBranch        *b_MET_phi_jesRelativeSampleDown;   //!
   TBranch        *b_Jet_pt_jesRelativeFSRDown;   //!
   TBranch        *b_Jet_mass_jesRelativeFSRDown;   //!
   TBranch        *b_MET_pt_jesRelativeFSRDown;   //!
   TBranch        *b_MET_phi_jesRelativeFSRDown;   //!
   TBranch        *b_Jet_pt_jesRelativeStatFSRDown;   //!
   TBranch        *b_Jet_mass_jesRelativeStatFSRDown;   //!
   TBranch        *b_MET_pt_jesRelativeStatFSRDown;   //!
   TBranch        *b_MET_phi_jesRelativeStatFSRDown;   //!
   TBranch        *b_Jet_pt_jesRelativeStatECDown;   //!
   TBranch        *b_Jet_mass_jesRelativeStatECDown;   //!
   TBranch        *b_MET_pt_jesRelativeStatECDown;   //!
   TBranch        *b_MET_phi_jesRelativeStatECDown;   //!
   TBranch        *b_Jet_pt_jesRelativeStatHFDown;   //!
   TBranch        *b_Jet_mass_jesRelativeStatHFDown;   //!
   TBranch        *b_MET_pt_jesRelativeStatHFDown;   //!
   TBranch        *b_MET_phi_jesRelativeStatHFDown;   //!
   TBranch        *b_Jet_pt_jesPileUpDataMCDown;   //!
   TBranch        *b_Jet_mass_jesPileUpDataMCDown;   //!
   TBranch        *b_MET_pt_jesPileUpDataMCDown;   //!
   TBranch        *b_MET_phi_jesPileUpDataMCDown;   //!
   TBranch        *b_Jet_pt_jesPileUpPtRefDown;   //!
   TBranch        *b_Jet_mass_jesPileUpPtRefDown;   //!
   TBranch        *b_MET_pt_jesPileUpPtRefDown;   //!
   TBranch        *b_MET_phi_jesPileUpPtRefDown;   //!
   TBranch        *b_Jet_pt_jesPileUpPtBBDown;   //!
   TBranch        *b_Jet_mass_jesPileUpPtBBDown;   //!
   TBranch        *b_MET_pt_jesPileUpPtBBDown;   //!
   TBranch        *b_MET_phi_jesPileUpPtBBDown;   //!
   TBranch        *b_Jet_pt_jesPileUpPtEC1Down;   //!
   TBranch        *b_Jet_mass_jesPileUpPtEC1Down;   //!
   TBranch        *b_MET_pt_jesPileUpPtEC1Down;   //!
   TBranch        *b_MET_phi_jesPileUpPtEC1Down;   //!
   TBranch        *b_Jet_pt_jesPileUpPtEC2Down;   //!
   TBranch        *b_Jet_mass_jesPileUpPtEC2Down;   //!
   TBranch        *b_MET_pt_jesPileUpPtEC2Down;   //!
   TBranch        *b_MET_phi_jesPileUpPtEC2Down;   //!
   TBranch        *b_Jet_pt_jesPileUpPtHFDown;   //!
   TBranch        *b_Jet_mass_jesPileUpPtHFDown;   //!
   TBranch        *b_MET_pt_jesPileUpPtHFDown;   //!
   TBranch        *b_MET_phi_jesPileUpPtHFDown;   //!
   TBranch        *b_Jet_pt_jesPileUpMuZeroDown;   //!
   TBranch        *b_Jet_mass_jesPileUpMuZeroDown;   //!
   TBranch        *b_MET_pt_jesPileUpMuZeroDown;   //!
   TBranch        *b_MET_phi_jesPileUpMuZeroDown;   //!
   TBranch        *b_Jet_pt_jesPileUpEnvelopeDown;   //!
   TBranch        *b_Jet_mass_jesPileUpEnvelopeDown;   //!
   TBranch        *b_MET_pt_jesPileUpEnvelopeDown;   //!
   TBranch        *b_MET_phi_jesPileUpEnvelopeDown;   //!
   TBranch        *b_Jet_pt_jesSubTotalPileUpDown;   //!
   TBranch        *b_Jet_mass_jesSubTotalPileUpDown;   //!
   TBranch        *b_MET_pt_jesSubTotalPileUpDown;   //!
   TBranch        *b_MET_phi_jesSubTotalPileUpDown;   //!
   TBranch        *b_Jet_pt_jesSubTotalRelativeDown;   //!
   TBranch        *b_Jet_mass_jesSubTotalRelativeDown;   //!
   TBranch        *b_MET_pt_jesSubTotalRelativeDown;   //!
   TBranch        *b_MET_phi_jesSubTotalRelativeDown;   //!
   TBranch        *b_Jet_pt_jesSubTotalPtDown;   //!
   TBranch        *b_Jet_mass_jesSubTotalPtDown;   //!
   TBranch        *b_MET_pt_jesSubTotalPtDown;   //!
   TBranch        *b_MET_phi_jesSubTotalPtDown;   //!
   TBranch        *b_Jet_pt_jesSubTotalScaleDown;   //!
   TBranch        *b_Jet_mass_jesSubTotalScaleDown;   //!
   TBranch        *b_MET_pt_jesSubTotalScaleDown;   //!
   TBranch        *b_MET_phi_jesSubTotalScaleDown;   //!
   TBranch        *b_Jet_pt_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_Jet_mass_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_MET_pt_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_MET_phi_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_Jet_pt_jesSubTotalMCDown;   //!
   TBranch        *b_Jet_mass_jesSubTotalMCDown;   //!
   TBranch        *b_MET_pt_jesSubTotalMCDown;   //!
   TBranch        *b_MET_phi_jesSubTotalMCDown;   //!
   TBranch        *b_Jet_pt_jesTotalDown;   //!
   TBranch        *b_Jet_mass_jesTotalDown;   //!
   TBranch        *b_MET_pt_jesTotalDown;   //!
   TBranch        *b_MET_phi_jesTotalDown;   //!
   TBranch        *b_Jet_pt_jesTotalNoFlavorDown;   //!
   TBranch        *b_Jet_mass_jesTotalNoFlavorDown;   //!
   TBranch        *b_MET_pt_jesTotalNoFlavorDown;   //!
   TBranch        *b_MET_phi_jesTotalNoFlavorDown;   //!
   TBranch        *b_Jet_pt_jesTotalNoTimeDown;   //!
   TBranch        *b_Jet_mass_jesTotalNoTimeDown;   //!
   TBranch        *b_MET_pt_jesTotalNoTimeDown;   //!
   TBranch        *b_MET_phi_jesTotalNoTimeDown;   //!
   TBranch        *b_Jet_pt_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_Jet_mass_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_MET_pt_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_MET_phi_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_Jet_pt_jesFlavorZJetDown;   //!
   TBranch        *b_Jet_mass_jesFlavorZJetDown;   //!
   TBranch        *b_MET_pt_jesFlavorZJetDown;   //!
   TBranch        *b_MET_phi_jesFlavorZJetDown;   //!
   TBranch        *b_Jet_pt_jesFlavorPhotonJetDown;   //!
   TBranch        *b_Jet_mass_jesFlavorPhotonJetDown;   //!
   TBranch        *b_MET_pt_jesFlavorPhotonJetDown;   //!
   TBranch        *b_MET_phi_jesFlavorPhotonJetDown;   //!
   TBranch        *b_Jet_pt_jesFlavorPureGluonDown;   //!
   TBranch        *b_Jet_mass_jesFlavorPureGluonDown;   //!
   TBranch        *b_MET_pt_jesFlavorPureGluonDown;   //!
   TBranch        *b_MET_phi_jesFlavorPureGluonDown;   //!
   TBranch        *b_Jet_pt_jesFlavorPureQuarkDown;   //!
   TBranch        *b_Jet_mass_jesFlavorPureQuarkDown;   //!
   TBranch        *b_MET_pt_jesFlavorPureQuarkDown;   //!
   TBranch        *b_MET_phi_jesFlavorPureQuarkDown;   //!
   TBranch        *b_Jet_pt_jesFlavorPureCharmDown;   //!
   TBranch        *b_Jet_mass_jesFlavorPureCharmDown;   //!
   TBranch        *b_MET_pt_jesFlavorPureCharmDown;   //!
   TBranch        *b_MET_phi_jesFlavorPureCharmDown;   //!
   TBranch        *b_Jet_pt_jesFlavorPureBottomDown;   //!
   TBranch        *b_Jet_mass_jesFlavorPureBottomDown;   //!
   TBranch        *b_MET_pt_jesFlavorPureBottomDown;   //!
   TBranch        *b_MET_phi_jesFlavorPureBottomDown;   //!
   TBranch        *b_Jet_pt_jesTimeRunBDown;   //!
   TBranch        *b_Jet_mass_jesTimeRunBDown;   //!
   TBranch        *b_MET_pt_jesTimeRunBDown;   //!
   TBranch        *b_MET_phi_jesTimeRunBDown;   //!
   TBranch        *b_Jet_pt_jesTimeRunCDown;   //!
   TBranch        *b_Jet_mass_jesTimeRunCDown;   //!
   TBranch        *b_MET_pt_jesTimeRunCDown;   //!
   TBranch        *b_MET_phi_jesTimeRunCDown;   //!
   TBranch        *b_Jet_pt_jesTimeRunDEDown;   //!
   TBranch        *b_Jet_mass_jesTimeRunDEDown;   //!
   TBranch        *b_MET_pt_jesTimeRunDEDown;   //!
   TBranch        *b_MET_phi_jesTimeRunDEDown;   //!
   TBranch        *b_Jet_pt_jesTimeRunFDown;   //!
   TBranch        *b_Jet_mass_jesTimeRunFDown;   //!
   TBranch        *b_MET_pt_jesTimeRunFDown;   //!
   TBranch        *b_MET_phi_jesTimeRunFDown;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupUncorrelatedDown;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupUncorrelatedDown;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupUncorrelatedDown;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupUncorrelatedDown;   //!
   TBranch        *b_MET_pt_unclustEnDown;   //!
   TBranch        *b_MET_phi_unclustEnDown;   //!
   
   TBranch        *b_FatJet_pt_raw;   //!
   TBranch        *b_FatJet_pt_nom;   //!
   TBranch        *b_FatJet_mass_raw;   //!
   TBranch        *b_FatJet_mass_nom;   //!
   TBranch        *b_FatJet_corr_JEC;   //!
   TBranch        *b_FatJet_corr_JER;   //!
   TBranch        *b_FatJet_corr_JMS;   //!
   TBranch        *b_FatJet_corr_JMR;   //!
   TBranch        *b_FatJet_msoftdrop_raw;   //!
   TBranch        *b_FatJet_msoftdrop_nom;   //!
   TBranch        *b_FatJet_msoftdrop_corr_JMR;   //!
   TBranch        *b_FatJet_msoftdrop_corr_JMS;   //!
   TBranch        *b_FatJet_msoftdrop_corr_PUPPI;   //!
   TBranch        *b_FatJet_msoftdrop_tau21DDT_nom;   //!

   TBranch        *b_FatJet_pt_jerUp;   //!
   TBranch        *b_FatJet_mass_jerUp;   //!
   TBranch        *b_FatJet_mass_jmrUp;   //!
   TBranch        *b_FatJet_mass_jmsUp;   //!
   TBranch        *b_FatJet_msoftdrop_jerUp;   //!
   TBranch        *b_FatJet_msoftdrop_jmrUp;   //!
   TBranch        *b_FatJet_msoftdrop_jmsUp;   //!
   TBranch        *b_FatJet_pt_jesAbsoluteStatUp;   //!
   TBranch        *b_FatJet_mass_jesAbsoluteStatUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesAbsoluteStatUp;   //!
   TBranch        *b_FatJet_pt_jesAbsoluteScaleUp;   //!
   TBranch        *b_FatJet_mass_jesAbsoluteScaleUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesAbsoluteScaleUp;   //!
   TBranch        *b_FatJet_pt_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_FatJet_mass_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_FatJet_pt_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_FatJet_mass_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_FatJet_pt_jesFragmentationUp;   //!
   TBranch        *b_FatJet_mass_jesFragmentationUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesFragmentationUp;   //!
   TBranch        *b_FatJet_pt_jesSinglePionECALUp;   //!
   TBranch        *b_FatJet_mass_jesSinglePionECALUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesSinglePionECALUp;   //!
   TBranch        *b_FatJet_pt_jesSinglePionHCALUp;   //!
   TBranch        *b_FatJet_mass_jesSinglePionHCALUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesSinglePionHCALUp;   //!
   TBranch        *b_FatJet_pt_jesFlavorQCDUp;   //!
   TBranch        *b_FatJet_mass_jesFlavorQCDUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorQCDUp;   //!
   TBranch        *b_FatJet_pt_jesTimePtEtaUp;   //!
   TBranch        *b_FatJet_mass_jesTimePtEtaUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimePtEtaUp;   //!
   TBranch        *b_FatJet_pt_jesRelativeJEREC1Up;   //!
   TBranch        *b_FatJet_mass_jesRelativeJEREC1Up;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeJEREC1Up;   //!
   TBranch        *b_FatJet_pt_jesRelativeJEREC2Up;   //!
   TBranch        *b_FatJet_mass_jesRelativeJEREC2Up;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeJEREC2Up;   //!
   TBranch        *b_FatJet_pt_jesRelativeJERHFUp;   //!
   TBranch        *b_FatJet_mass_jesRelativeJERHFUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeJERHFUp;   //!
   TBranch        *b_FatJet_pt_jesRelativePtBBUp;   //!
   TBranch        *b_FatJet_mass_jesRelativePtBBUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativePtBBUp;   //!
   TBranch        *b_FatJet_pt_jesRelativePtEC1Up;   //!
   TBranch        *b_FatJet_mass_jesRelativePtEC1Up;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativePtEC1Up;   //!
   TBranch        *b_FatJet_pt_jesRelativePtEC2Up;   //!
   TBranch        *b_FatJet_mass_jesRelativePtEC2Up;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativePtEC2Up;   //!
   TBranch        *b_FatJet_pt_jesRelativePtHFUp;   //!
   TBranch        *b_FatJet_mass_jesRelativePtHFUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativePtHFUp;   //!
   TBranch        *b_FatJet_pt_jesRelativeBalUp;   //!
   TBranch        *b_FatJet_mass_jesRelativeBalUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeBalUp;   //!
   TBranch        *b_FatJet_pt_jesRelativeSampleUp;   //!
   TBranch        *b_FatJet_mass_jesRelativeSampleUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeSampleUp;   //!
   TBranch        *b_FatJet_pt_jesRelativeFSRUp;   //!
   TBranch        *b_FatJet_mass_jesRelativeFSRUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeFSRUp;   //!
   TBranch        *b_FatJet_pt_jesRelativeStatFSRUp;   //!
   TBranch        *b_FatJet_mass_jesRelativeStatFSRUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeStatFSRUp;   //!
   TBranch        *b_FatJet_pt_jesRelativeStatECUp;   //!
   TBranch        *b_FatJet_mass_jesRelativeStatECUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeStatECUp;   //!
   TBranch        *b_FatJet_pt_jesRelativeStatHFUp;   //!
   TBranch        *b_FatJet_mass_jesRelativeStatHFUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeStatHFUp;   //!
   TBranch        *b_FatJet_pt_jesPileUpDataMCUp;   //!
   TBranch        *b_FatJet_mass_jesPileUpDataMCUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpDataMCUp;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtRefUp;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtRefUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtRefUp;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtBBUp;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtBBUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtBBUp;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtEC1Up;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtEC1Up;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtEC1Up;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtEC2Up;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtEC2Up;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtEC2Up;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtHFUp;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtHFUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtHFUp;   //!
   TBranch        *b_FatJet_pt_jesPileUpMuZeroUp;   //!
   TBranch        *b_FatJet_mass_jesPileUpMuZeroUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpMuZeroUp;   //!
   TBranch        *b_FatJet_pt_jesPileUpEnvelopeUp;   //!
   TBranch        *b_FatJet_mass_jesPileUpEnvelopeUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpEnvelopeUp;   //!
   TBranch        *b_FatJet_pt_jesSubTotalPileUpUp;   //!
   TBranch        *b_FatJet_mass_jesSubTotalPileUpUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalPileUpUp;   //!
   TBranch        *b_FatJet_pt_jesSubTotalRelativeUp;   //!
   TBranch        *b_FatJet_mass_jesSubTotalRelativeUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalRelativeUp;   //!
   TBranch        *b_FatJet_pt_jesSubTotalPtUp;   //!
   TBranch        *b_FatJet_mass_jesSubTotalPtUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalPtUp;   //!
   TBranch        *b_FatJet_pt_jesSubTotalScaleUp;   //!
   TBranch        *b_FatJet_mass_jesSubTotalScaleUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalScaleUp;   //!
   TBranch        *b_FatJet_pt_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_FatJet_mass_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_FatJet_pt_jesSubTotalMCUp;   //!
   TBranch        *b_FatJet_mass_jesSubTotalMCUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalMCUp;   //!
   TBranch        *b_FatJet_pt_jesTotalUp;   //!
   TBranch        *b_FatJet_mass_jesTotalUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesTotalUp;   //!
   TBranch        *b_FatJet_pt_jesTotalNoFlavorUp;   //!
   TBranch        *b_FatJet_mass_jesTotalNoFlavorUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesTotalNoFlavorUp;   //!
   TBranch        *b_FatJet_pt_jesTotalNoTimeUp;   //!
   TBranch        *b_FatJet_mass_jesTotalNoTimeUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesTotalNoTimeUp;   //!
   TBranch        *b_FatJet_pt_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_FatJet_mass_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_FatJet_pt_jesFlavorZJetUp;   //!
   TBranch        *b_FatJet_mass_jesFlavorZJetUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorZJetUp;   //!
   TBranch        *b_FatJet_pt_jesFlavorPhotonJetUp;   //!
   TBranch        *b_FatJet_mass_jesFlavorPhotonJetUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPhotonJetUp;   //!
   TBranch        *b_FatJet_pt_jesFlavorPureGluonUp;   //!
   TBranch        *b_FatJet_mass_jesFlavorPureGluonUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPureGluonUp;   //!
   TBranch        *b_FatJet_pt_jesFlavorPureQuarkUp;   //!
   TBranch        *b_FatJet_mass_jesFlavorPureQuarkUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPureQuarkUp;   //!
   TBranch        *b_FatJet_pt_jesFlavorPureCharmUp;   //!
   TBranch        *b_FatJet_mass_jesFlavorPureCharmUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPureCharmUp;   //!
   TBranch        *b_FatJet_pt_jesFlavorPureBottomUp;   //!
   TBranch        *b_FatJet_mass_jesFlavorPureBottomUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPureBottomUp;   //!
   TBranch        *b_FatJet_pt_jesTimeRunBUp;   //!
   TBranch        *b_FatJet_mass_jesTimeRunBUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimeRunBUp;   //!
   TBranch        *b_FatJet_pt_jesTimeRunCUp;   //!
   TBranch        *b_FatJet_mass_jesTimeRunCUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimeRunCUp;   //!
   TBranch        *b_FatJet_pt_jesTimeRunDEUp;   //!
   TBranch        *b_FatJet_mass_jesTimeRunDEUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimeRunDEUp;   //!
   TBranch        *b_FatJet_pt_jesTimeRunFUp;   //!
   TBranch        *b_FatJet_mass_jesTimeRunFUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimeRunFUp;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_FatJet_pt_jerDown;   //!
   TBranch        *b_FatJet_mass_jerDown;   //!
   TBranch        *b_FatJet_mass_jmrDown;   //!
   TBranch        *b_FatJet_mass_jmsDown;   //!
   TBranch        *b_FatJet_msoftdrop_jerDown;   //!
   TBranch        *b_FatJet_msoftdrop_jmrDown;   //!
   TBranch        *b_FatJet_msoftdrop_jmsDown;   //!
   TBranch        *b_FatJet_pt_jesAbsoluteStatDown;   //!
   TBranch        *b_FatJet_mass_jesAbsoluteStatDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesAbsoluteStatDown;   //!
   TBranch        *b_FatJet_pt_jesAbsoluteScaleDown;   //!
   TBranch        *b_FatJet_mass_jesAbsoluteScaleDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesAbsoluteScaleDown;   //!
   TBranch        *b_FatJet_pt_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_FatJet_mass_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_FatJet_pt_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_FatJet_mass_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_FatJet_pt_jesFragmentationDown;   //!
   TBranch        *b_FatJet_mass_jesFragmentationDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesFragmentationDown;   //!
   TBranch        *b_FatJet_pt_jesSinglePionECALDown;   //!
   TBranch        *b_FatJet_mass_jesSinglePionECALDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesSinglePionECALDown;   //!
   TBranch        *b_FatJet_pt_jesSinglePionHCALDown;   //!
   TBranch        *b_FatJet_mass_jesSinglePionHCALDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesSinglePionHCALDown;   //!
   TBranch        *b_FatJet_pt_jesFlavorQCDDown;   //!
   TBranch        *b_FatJet_mass_jesFlavorQCDDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorQCDDown;   //!
   TBranch        *b_FatJet_pt_jesTimePtEtaDown;   //!
   TBranch        *b_FatJet_mass_jesTimePtEtaDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimePtEtaDown;   //!
   TBranch        *b_FatJet_pt_jesRelativeJEREC1Down;   //!
   TBranch        *b_FatJet_mass_jesRelativeJEREC1Down;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeJEREC1Down;   //!
   TBranch        *b_FatJet_pt_jesRelativeJEREC2Down;   //!
   TBranch        *b_FatJet_mass_jesRelativeJEREC2Down;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeJEREC2Down;   //!
   TBranch        *b_FatJet_pt_jesRelativeJERHFDown;   //!
   TBranch        *b_FatJet_mass_jesRelativeJERHFDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeJERHFDown;   //!
   TBranch        *b_FatJet_pt_jesRelativePtBBDown;   //!
   TBranch        *b_FatJet_mass_jesRelativePtBBDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativePtBBDown;   //!
   TBranch        *b_FatJet_pt_jesRelativePtEC1Down;   //!
   TBranch        *b_FatJet_mass_jesRelativePtEC1Down;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativePtEC1Down;   //!
   TBranch        *b_FatJet_pt_jesRelativePtEC2Down;   //!
   TBranch        *b_FatJet_mass_jesRelativePtEC2Down;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativePtEC2Down;   //!
   TBranch        *b_FatJet_pt_jesRelativePtHFDown;   //!
   TBranch        *b_FatJet_mass_jesRelativePtHFDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativePtHFDown;   //!
   TBranch        *b_FatJet_pt_jesRelativeBalDown;   //!
   TBranch        *b_FatJet_mass_jesRelativeBalDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeBalDown;   //!
   TBranch        *b_FatJet_pt_jesRelativeSampleDown;   //!
   TBranch        *b_FatJet_mass_jesRelativeSampleDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeSampleDown;   //!
   TBranch        *b_FatJet_pt_jesRelativeFSRDown;   //!
   TBranch        *b_FatJet_mass_jesRelativeFSRDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeFSRDown;   //!
   TBranch        *b_FatJet_pt_jesRelativeStatFSRDown;   //!
   TBranch        *b_FatJet_mass_jesRelativeStatFSRDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeStatFSRDown;   //!
   TBranch        *b_FatJet_pt_jesRelativeStatECDown;   //!
   TBranch        *b_FatJet_mass_jesRelativeStatECDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeStatECDown;   //!
   TBranch        *b_FatJet_pt_jesRelativeStatHFDown;   //!
   TBranch        *b_FatJet_mass_jesRelativeStatHFDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeStatHFDown;   //!
   TBranch        *b_FatJet_pt_jesPileUpDataMCDown;   //!
   TBranch        *b_FatJet_mass_jesPileUpDataMCDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpDataMCDown;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtRefDown;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtRefDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtRefDown;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtBBDown;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtBBDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtBBDown;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtEC1Down;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtEC1Down;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtEC1Down;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtEC2Down;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtEC2Down;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtEC2Down;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtHFDown;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtHFDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtHFDown;   //!
   TBranch        *b_FatJet_pt_jesPileUpMuZeroDown;   //!
   TBranch        *b_FatJet_mass_jesPileUpMuZeroDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpMuZeroDown;   //!
   TBranch        *b_FatJet_pt_jesPileUpEnvelopeDown;   //!
   TBranch        *b_FatJet_mass_jesPileUpEnvelopeDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpEnvelopeDown;   //!
   TBranch        *b_FatJet_pt_jesSubTotalPileUpDown;   //!
   TBranch        *b_FatJet_mass_jesSubTotalPileUpDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalPileUpDown;   //!
   TBranch        *b_FatJet_pt_jesSubTotalRelativeDown;   //!
   TBranch        *b_FatJet_mass_jesSubTotalRelativeDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalRelativeDown;   //!
   TBranch        *b_FatJet_pt_jesSubTotalPtDown;   //!
   TBranch        *b_FatJet_mass_jesSubTotalPtDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalPtDown;   //!
   TBranch        *b_FatJet_pt_jesSubTotalScaleDown;   //!
   TBranch        *b_FatJet_mass_jesSubTotalScaleDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalScaleDown;   //!
   TBranch        *b_FatJet_pt_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_FatJet_mass_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_FatJet_pt_jesSubTotalMCDown;   //!
   TBranch        *b_FatJet_mass_jesSubTotalMCDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalMCDown;   //!
   TBranch        *b_FatJet_pt_jesTotalDown;   //!
   TBranch        *b_FatJet_mass_jesTotalDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesTotalDown;   //!
   TBranch        *b_FatJet_pt_jesTotalNoFlavorDown;   //!
   TBranch        *b_FatJet_mass_jesTotalNoFlavorDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesTotalNoFlavorDown;   //!
   TBranch        *b_FatJet_pt_jesTotalNoTimeDown;   //!
   TBranch        *b_FatJet_mass_jesTotalNoTimeDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesTotalNoTimeDown;   //!
   TBranch        *b_FatJet_pt_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_FatJet_mass_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_FatJet_pt_jesFlavorZJetDown;   //!
   TBranch        *b_FatJet_mass_jesFlavorZJetDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorZJetDown;   //!
   TBranch        *b_FatJet_pt_jesFlavorPhotonJetDown;   //!
   TBranch        *b_FatJet_mass_jesFlavorPhotonJetDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPhotonJetDown;   //!
   TBranch        *b_FatJet_pt_jesFlavorPureGluonDown;   //!
   TBranch        *b_FatJet_mass_jesFlavorPureGluonDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPureGluonDown;   //!
   TBranch        *b_FatJet_pt_jesFlavorPureQuarkDown;   //!
   TBranch        *b_FatJet_mass_jesFlavorPureQuarkDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPureQuarkDown;   //!
   TBranch        *b_FatJet_pt_jesFlavorPureCharmDown;   //!
   TBranch        *b_FatJet_mass_jesFlavorPureCharmDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPureCharmDown;   //!
   TBranch        *b_FatJet_pt_jesFlavorPureBottomDown;   //!
   TBranch        *b_FatJet_mass_jesFlavorPureBottomDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPureBottomDown;   //!
   TBranch        *b_FatJet_pt_jesTimeRunBDown;   //!
   TBranch        *b_FatJet_mass_jesTimeRunBDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimeRunBDown;   //!
   TBranch        *b_FatJet_pt_jesTimeRunCDown;   //!
   TBranch        *b_FatJet_mass_jesTimeRunCDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimeRunCDown;   //!
   TBranch        *b_FatJet_pt_jesTimeRunDEDown;   //!
   TBranch        *b_FatJet_mass_jesTimeRunDEDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimeRunDEDown;   //!
   TBranch        *b_FatJet_pt_jesTimeRunFDown;   //!
   TBranch        *b_FatJet_mass_jesTimeRunFDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimeRunFDown;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupUncorrelatedDown;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupUncorrelatedDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupUncorrelatedDown;   //!
   TBranch        *b_Jet_CSVbtagSF;   //!
   TBranch        *b_Jet_CSVbtagSF_up;   //!
   TBranch        *b_Jet_CSVbtagSF_down;   //!
   TBranch        *b_Jet_CSVbtagSF_shape;   //!
   TBranch        *b_Jet_CSVbtagSF_shape_up_jes;   //!
   TBranch        *b_Jet_CSVbtagSF_shape_down_jes;   //!
   TBranch        *b_Jet_CSVbtagSF_shape_up_lf;   //!
   TBranch        *b_Jet_CSVbtagSF_shape_down_lf;   //!
   TBranch        *b_Jet_CSVbtagSF_shape_up_hf;   //!
   TBranch        *b_Jet_CSVbtagSF_shape_down_hf;   //!
   TBranch        *b_Jet_CSVbtagSF_shape_up_hfstats1;   //!
   TBranch        *b_Jet_CSVbtagSF_shape_down_hfstats1;   //!
   TBranch        *b_Jet_CSVbtagSF_shape_up_hfstats2;   //!
   TBranch        *b_Jet_CSVbtagSF_shape_down_hfstats2;   //!
   TBranch        *b_Jet_CSVbtagSF_shape_up_lfstats1;   //!
   TBranch        *b_Jet_CSVbtagSF_shape_down_lfstats1;   //!
   TBranch        *b_Jet_CSVbtagSF_shape_up_lfstats2;   //!
   TBranch        *b_Jet_CSVbtagSF_shape_down_lfstats2;   //!
   TBranch        *b_Jet_CSVbtagSF_shape_up_cferr1;   //!
   TBranch        *b_Jet_CSVbtagSF_shape_down_cferr1;   //!
   TBranch        *b_Jet_CSVbtagSF_shape_up_cferr2;   //!
   TBranch        *b_Jet_CSVbtagSF_shape_down_cferr2;   //!
   TBranch        *b_Jet_deepCSVbtagSF;   //!
   TBranch        *b_Jet_deepCSVbtagSF_up;   //!
   TBranch        *b_Jet_deepCSVbtagSF_down;   //!
   TBranch        *b_Jet_deepCSVbtagSF_shape;   //!
   TBranch        *b_Jet_deepCSVbtagSF_shape_up_jes;   //!
   TBranch        *b_Jet_deepCSVbtagSF_shape_down_jes;   //!
   TBranch        *b_Jet_deepCSVbtagSF_shape_up_lf;   //!
   TBranch        *b_Jet_deepCSVbtagSF_shape_down_lf;   //!
   TBranch        *b_Jet_deepCSVbtagSF_shape_up_hf;   //!
   TBranch        *b_Jet_deepCSVbtagSF_shape_down_hf;   //!
   TBranch        *b_Jet_deepCSVbtagSF_shape_up_hfstats1;   //!
   TBranch        *b_Jet_deepCSVbtagSF_shape_down_hfstats1;   //!
   TBranch        *b_Jet_deepCSVbtagSF_shape_up_hfstats2;   //!
   TBranch        *b_Jet_deepCSVbtagSF_shape_down_hfstats2;   //!
   TBranch        *b_Jet_deepCSVbtagSF_shape_up_lfstats1;   //!
   TBranch        *b_Jet_deepCSVbtagSF_shape_down_lfstats1;   //!
   TBranch        *b_Jet_deepCSVbtagSF_shape_up_lfstats2;   //!
   TBranch        *b_Jet_deepCSVbtagSF_shape_down_lfstats2;   //!
   TBranch        *b_Jet_deepCSVbtagSF_shape_up_cferr1;   //!
   TBranch        *b_Jet_deepCSVbtagSF_shape_down_cferr1;   //!
   TBranch        *b_Jet_deepCSVbtagSF_shape_up_cferr2;   //!
   TBranch        *b_Jet_deepCSVbtagSF_shape_down_cferr2;   //!
   TBranch        *b_Jet_deepflavbtagSF;   //!
   TBranch        *b_Jet_deepflavbtagSF_up;   //!
   TBranch        *b_Jet_deepflavbtagSF_down;   //!
   TBranch        *b_Jet_deepflavbtagSF_shape;   //!
   TBranch        *b_Jet_deepflavbtagSF_shape_up_jes;   //!
   TBranch        *b_Jet_deepflavbtagSF_shape_down_jes;   //!
   TBranch        *b_Jet_deepflavbtagSF_shape_up_lf;   //!
   TBranch        *b_Jet_deepflavbtagSF_shape_down_lf;   //!
   TBranch        *b_Jet_deepflavbtagSF_shape_up_hf;   //!
   TBranch        *b_Jet_deepflavbtagSF_shape_down_hf;   //!
   TBranch        *b_Jet_deepflavbtagSF_shape_up_hfstats1;   //!
   TBranch        *b_Jet_deepflavbtagSF_shape_down_hfstats1;   //!
   TBranch        *b_Jet_deepflavbtagSF_shape_up_hfstats2;   //!
   TBranch        *b_Jet_deepflavbtagSF_shape_down_hfstats2;   //!
   TBranch        *b_Jet_deepflavbtagSF_shape_up_lfstats1;   //!
   TBranch        *b_Jet_deepflavbtagSF_shape_down_lfstats1;   //!
   TBranch        *b_Jet_deepflavbtagSF_shape_up_lfstats2;   //!
   TBranch        *b_Jet_deepflavbtagSF_shape_down_lfstats2;   //!
   TBranch        *b_Jet_deepflavbtagSF_shape_up_cferr1;   //!
   TBranch        *b_Jet_deepflavbtagSF_shape_down_cferr1;   //!
   TBranch        *b_Jet_deepflavbtagSF_shape_up_cferr2;   //!
   TBranch        *b_Jet_deepflavbtagSF_shape_down_cferr2;   //!
   TBranch        *b_puWeight;   //!
   TBranch        *b_puWeightUp;   //!
   TBranch        *b_puWeightDown;   //!
   TBranch        *b_PrefireWeight;   //!
   TBranch        *b_PrefireWeight_Up;   //!
   TBranch        *b_PrefireWeight_Down;   //!
   
   TProofOutputFile *OutFile;
   TFile *fileOut;
   
   TTree *Tout ;
   TTree *Tout1 ;

   static const int netarange = 3;

   TH1D *hist_btAK4_pt[netarange];
   TH1D *hist_bAK4_pt[netarange];
   TH1D *hist_dbAK4_pt[netarange];
   TH1D *hist_dfbAK4_pt[netarange];
   TH1D *hist_qtAK4_pt[netarange];
   TH1D *hist_qAK4_pt[netarange];
   TH1D *hist_dqAK4_pt[netarange];
   TH1D *hist_dfqAK4_pt[netarange];
   
   static const int noperf_ptbins = 4;
   
   TH1D *hist_btag_csv_bhadron[noperf_ptbins+1];
   TH1D *hist_btag_deepcsv_bhadron[noperf_ptbins+1];
   TH1D *hist_btag_deepflav_bhadron[noperf_ptbins+1];
   
   TH1D *hist_btag_csv_qhadron[noperf_ptbins+1];
   TH1D *hist_btag_deepcsv_qhadron[noperf_ptbins+1];
   TH1D *hist_btag_deepflav_qhadron[noperf_ptbins+1];

   TH1D *hist_metfilter_pass;
   TH1D *hist_mutrig_pass;
   TH1D *hist_etrig_pass;
   TH1D *hist_photrig_pass;
   TH1D *hist_hadtrig_pass;

   TH1D *met_corpt;
   TH1D *met_corphi;
   TH1D *met_bysumEt;
   
   TH1D *hist_npv;
   TH1D *hist_npv_final;
   TH1D *hist_npv_nopuwt;
   TH1D *hist_npv_final_nopuwt;
   TH1D *hist_npu;
   TH1D *hist_npu_nopuwt;

   TH1D *hist_nmuons;
   TH1D *hist_nelectrons;
   TH1D *hist_nphotons;
   
   TH1D *hist_mass_hadtopq;
   TH1D *hist_pid_hadtopq;
   TH1D *hist_matchjet_hadtopq;
   
   TH1D *hist_jetpt_hastop;
   TH1D *hist_jetpt_hastop_DAK8pass;
   
   TH1D *hist_jetsdmass_top[noperf_ptbins];
   TH1D *hist_jettau32_top[noperf_ptbins];
   TH1D *hist_jetsubbtag_CSVv2_top[noperf_ptbins];
   TH1D *hist_jetsubbtag_DeepCSV_top[noperf_ptbins];
   TH1D *hist_jetDeepAK8_MD_top[noperf_ptbins];
   TH1D *hist_jetDeepAK8_MD_wm_top[noperf_ptbins];
   TH1D *hist_jetDeepAK8_top[noperf_ptbins];
   TH1D *hist_jetDeepAK8_wm_top[noperf_ptbins];
   TH1D *hist_jetgenmass_top[noperf_ptbins];
   TH1D *hist_jetDeepAK8_top_tausub_11[noperf_ptbins];
   TH1D *hist_jetDeepAK8_top_tausub_12[noperf_ptbins];
   TH1D *hist_jetDeepAK8_top_tausub_21[noperf_ptbins];
   TH1D *hist_jetDeepAK8_top_tausub_22[noperf_ptbins];
   
   TH1D *hist_mass_qg;
   TH1D *hist_pid_qg;
   
   TH1D *hist_jetsdmass_tbkg[noperf_ptbins];
   TH1D *hist_jettau32_tbkg[noperf_ptbins];
   TH1D *hist_jetsubbtag_CSVv2_tbkg[noperf_ptbins];
   TH1D *hist_jetsubbtag_DeepCSV_tbkg[noperf_ptbins];
   TH1D *hist_jetDeepAK8_MD_tbkg[noperf_ptbins];
   TH1D *hist_jetDeepAK8_MD_wm_tbkg[noperf_ptbins];
   TH1D *hist_jetDeepAK8_tbkg[noperf_ptbins];
   TH1D *hist_jetDeepAK8_wm_tbkg[noperf_ptbins];
   TH1D *hist_jetgenmass_tbkg[noperf_ptbins];
   TH1D *hist_jetDeepAK8_tbkg_tausub_11[noperf_ptbins];
   TH1D *hist_jetDeepAK8_tbkg_tausub_12[noperf_ptbins];
   TH1D *hist_jetDeepAK8_tbkg_tausub_21[noperf_ptbins];
   TH1D *hist_jetDeepAK8_tbkg_tausub_22[noperf_ptbins];

   TH1D *hist_jetpt_top;
   TH1D *hist_jetpt_msd_top;
   TH1D *hist_jetpt_top_md_deepak8_pass;
   TH1D *hist_jetpt_msd_top_md_deepak8_pass;
   
   TH1D *hist_jetpt_all;
   TH1D *hist_jetpt_all_DAK8pass;
   TH1D *hist_jetpt_topmsd_all;
   TH1D *hist_jetpt_topmsd_all_DAK8pass;

   TH1D *hist_njetAK8;
   TH1D *hist_topcand_AK8;
   TH1D *hist_jetptAK8;
   TH1D *hist_jetrapAK8;
   TH1D *hist_jetmassAK8;
   TH1D *hist_jetsdmassAK8;
   TH1D *hist_jetrhoAK8;
   TH1D *hist_jettau32AK8;
   TH1D *hist_jettau32AK8_mpass;
   TH1D *hist_jettau32AK8_mfail;
   TH1D *hist_jettau21AK8;
   TH1D *hist_njetsubjets_AK8_mpass;
   TH1D *hist_njetsubjets_AK8_lpass_mfail;
   TH1D *hist_jetsubbtagAK8;
   TH1D *hist_jetsubmassAK8;
   TH1D *hist_jetsubthetAK8;
   TH1D *hist_topjetsdmass;
   TH1D *hist_topjettau32AK8;
   TH1D *hist_topjetdeeptopscore;
   TH1D *hist_topjetmddeeptopscore;
   TH1D *hist_topjetgenmass;
   TH1D *hist_topjetgenpid;
   TH1D *hist_topjetpt;
   TH1D *hist_topjetptmatch;
   TH1D *hist_toppartonpt;
   TH1D *hist_topparton_matchedtopjet_pt;
   TH2D *hist_topjetptmsd_2d;
   TH2D *hist_2D_topjet_subjetmass12;
   TH2D *hist_2D_topjetsdmass_subjetmass1;
   TH2D *hist_2D_topjetsdmass_subjetmass2;
   
   TH1D *hist_njetsubjets_AK8_BMPpass_mpass;
   TH1D *hist_njetsubjets_AK8_BLPpass_BMPfail_mpass;
   TH1D *hist_nbsubjetAK8_mpass;
   TH1D *hist_jetsubmassAK8_mpass;
   TH1D *hist_jetsubptAK8_mpass;
   TH1D *hist_jetsubbtagAK8_mpass;
   TH1D *hist_jetsubbtagAK8_deepCSV_mpass;
   TH1D *hist_jetsubptAK8_btagpass_mpass;
   TH2D *hist_2D_tau32_subbtag_AK8_mpass;
   
   TH1D *hist_njetsubjets_AK8_BMPpass_mfail;
   TH1D *hist_njetsubjets_AK8_BLPpass_BMPfail_mfail;
   TH1D *hist_nbsubjetAK8_mfail;
   TH1D *hist_jetsubmassAK8_mfail;
   TH1D *hist_jetsubptAK8_mfail;
   TH1D *hist_jetsubbtagAK8_mfail;
   TH1D *hist_jetsubbtagAK8_deepCSV_mfail;
   TH1D *hist_jetsubptAK8_btagpass_mfail;
   TH2D *hist_2D_tau32_subbtag_AK8_mfail;
   
   TH2D *hist_2D_sdmass_tau32_AK8;
   TH2D *hist_2D_sdmass_subbtag_AK8;
   TH2D *hist_2D_tau32_subbtag_AK8;
   
   TH2D *hist_2D_sdmass_deeptopscore;
   TH2D *hist_2D_mass_deeptopscore;
   TH2D *hist_2D_tau32_deeptopscore;
   TH2D *hist_2D_subjetbtag_deeptopscore;
   
   TH2D *hist_2D_sdmass_deepmdtopscore;
   TH2D *hist_2D_mass_deepmdtopscore;
   TH2D *hist_2D_tau32_deepmdtopscore;
   TH2D *hist_2D_subjetbtag_deepmdtopscore;
   
   static const int ntautag = 2;
   static const int ntoptag = 8;
   static const int ntoptag_DAK8 = 8;
   static const int nbtag = 2;
   static const int nsbtag = 2;
   
   static const int ntopptbins = 3;
   float topptbins[ntopptbins+1] = {500,700,1000,7000};
   
   static const int ntopsdmassbins = 3;
   float topsdmassbins[ntopsdmassbins+1] = {0,105,210,400};
   
   static const int netabins = 3;
   float etabins[netabins+1] = {0,0.5,1.4,2.4};
  
   TH2D *hist_2D_sdmass_subbtag_CSVv2_AK8[nbtag][ntautag];
   TH2D *hist_2D_sdmass_subbtag_DeepCSV_AK8[nbtag][ntautag];
   TH2D *hist_2D_pt_btag_AK4[nsbtag][ntoptag];
   TH2D *hist_2D_sdmass_tau32_AK8_reg[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_AK8[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8[nbtag][nsbtag];
   
   TH2D *hist_2D_sdmass_subbtag_CSVv2_AK8_1[nbtag][ntautag];
   TH2D *hist_2D_sdmass_subbtag_DeepCSV_AK8_1[nbtag][ntautag];
   TH2D *hist_2D_sdmass_tbmass_1[nbtag][ntautag];
   TH2D *hist_2D_pt_btag_AK4_1[nsbtag][ntoptag];
   TH2D *hist_2D_DeepAK8_tbmass_1[nbtag];
   TH2D *hist_2D_MDDeepAK8_tbmass_1[nbtag];
   TH2D *hist_2D_sdmass_tau32_AK8_reg_1[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_AK8_1[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_1[nbtag][nsbtag];
   
   TH2D *hist_2D_sdmass_subbtag_CSVv2_AK8_2[nbtag][ntautag];
   TH2D *hist_2D_sdmass_subbtag_DeepCSV_AK8_2[nbtag][ntautag];
   TH2D *hist_2D_sdmass_tbmass_2[nbtag][ntautag];
   TH2D *hist_2D_pt_btag_AK4_2[nsbtag][ntoptag];
   TH2D *hist_2D_DeepAK8_tbmass_2[nbtag];
   TH2D *hist_2D_MDDeepAK8_tbmass_2[nbtag];
   TH2D *hist_2D_sdmass_tau32_AK8_reg_2[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_AK8_2[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_2[nbtag][nsbtag];
   
   TH2D *hist_2D_sdmass_subbtag_CSVv2_AK8_3[nbtag][ntautag];
   TH2D *hist_2D_sdmass_subbtag_DeepCSV_AK8_3[nbtag][ntautag];
   TH2D *hist_2D_sdmass_tbmass_3[nbtag][ntautag];
   TH2D *hist_2D_DeepAK8_tbmass_3[nbtag];
   TH2D *hist_2D_MDDeepAK8_tbmass_3[nbtag];
   TH2D *hist_2D_pt_btag_AK4_3[nsbtag][ntoptag];
   TH2D *hist_2D_sdmass_tau32_AK8_reg_3[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_AK8_3[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_3[nbtag][nsbtag];
   
   TH2D *hist_2D_sdmass_tbmass_deepAK8[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_tbmass_mddeepAK8[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_tbmass_mddeepAK8_3[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_tbmass_mddeepAK8_tt[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_tbmass_mddeepAK8_tt_3[nbtag][nsbtag];
   
   TH2D *hist_2D_sdmass_subbtag_CSVv2_AK8_DeepAK8[nbtag][ntautag];
   TH2D *hist_2D_sdmass_subbtag_DeepCSV_AK8_DeepAK8[nbtag][ntautag];
   TH2D *hist_2D_sdmass_tbmass[nbtag][ntautag];
   TH2D *hist_2D_DeepAK8_tbmass[nbtag];
   TH2D *hist_2D_MDDeepAK8_tbmass[nbtag];
  
   TH2D *hist_2D_pt_btag_AK4_DeepAK8[nsbtag][ntoptag];
   TH2D *hist_2D_sdmass_deeptag_AK8_reg_DeepAK8[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_AK8_DeepAK8[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_DeepAK8[nbtag][nsbtag];
   
   TH2D *hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8[nbtag][ntautag];
   TH2D *hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8[nbtag][ntautag];
   TH2D *hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8_3[nbtag][ntautag];
   TH2D *hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8_3[nbtag][ntautag];
   TH2D *hist_2D_pt_btag_AK4_md_DeepAK8[nsbtag][ntoptag];
   TH2D *hist_2D_pt_btag_AK4_md_DeepAK8_3[nsbtag][ntoptag];
   TH2D *hist_2D_pt_btag_AK4_md_DeepAK8_beta[nsbtag][ntoptag][netabins];
   TH2D *hist_2D_pt_btag_AK4_md_DeepAK8_beta_3[nsbtag][ntoptag][netabins];
   TH2D *hist_2D_sdmass_deeptag_AK8_reg_md_DeepAK8[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptag_AK8_reg_md_DeepAK8_3[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8_3[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_3[nbtag][nsbtag];
   
   TH2D *hist_2D_pt_btag_AK4_md_DeepAK8_Btag_up[nsbtag][ntoptag][nbtagmax];
   TH2D *hist_2D_pt_btag_AK4_md_DeepAK8_Btag_dn[nsbtag][ntoptag][nbtagmax];
   
   TH2D *hist_2D_pt_btag_AK4_md_DeepAK8_beta_Btag_up[nsbtag][ntoptag][netabins][nbtagmax];
   TH2D *hist_2D_pt_btag_AK4_md_DeepAK8_beta_Btag_dn[nsbtag][ntoptag][netabins][nbtagmax];
   
   TH2D *hist_2D_pt_btag_AK4_md_DeepAK8_beta_JES_up[nsbtag][ntoptag][netabins][njesmax];
   TH2D *hist_2D_pt_btag_AK4_md_DeepAK8_beta_JES_dn[nsbtag][ntoptag][netabins][njesmax];
   
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JES_up[nbtag][nsbtag][njesmax];
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JER_up[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JMS_up[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JMR_up[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_Btag_up[nbtag][nsbtag][nbtagmax];
   
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JES_dn[nbtag][nsbtag][njesmax];
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JER_dn[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JMS_dn[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JMR_dn[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_Btag_dn[nbtag][nsbtag][nbtagmax];
   
   TH2D *hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8_tt[nbtag][ntautag];
   TH2D *hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8_tt[nbtag][ntautag];
   TH2D *hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8_tt_3[nbtag][ntautag];
   TH2D *hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8_tt_3[nbtag][ntautag];
   TH2D *hist_2D_pt_btag_AK4_md_DeepAK8_tt[nsbtag][ntoptag];
   TH2D *hist_2D_pt_btag_AK4_md_DeepAK8_tt_3[nsbtag][ntoptag];
   TH2D *hist_2D_pt_btag_AK4_md_DeepAK8_beta_tt[nsbtag][ntoptag][netabins];
   TH2D *hist_2D_pt_btag_AK4_md_DeepAK8_beta_tt_3[nsbtag][ntoptag][netabins];
   TH2D *hist_2D_sdmass_deeptag_AK8_reg_md_DeepAK8_tt[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptag_AK8_reg_md_DeepAK8_tt_3[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8_tt[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8_tt_3[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_tt[nbtag][nsbtag];
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_tt_3[nbtag][nsbtag];
   
   
   TH1D *hist_delR_AK4_toptagAK8;
   TH1D *hist_deleta_AK4_toptagAK8;
   TH1D *hist_delphi_AK4_toptagAK8;
   TH1D *hist_delphi_AK4_toptagAK8_dRpass;
   TH1D *hist_deleta_AK4_toptagAK8_dRpass;
   TH1D *hist_deleta_AK4_toptagAK8_dRpass_mtb2000;
   
   TH1D *hist_delR_btag_toptag;
   TH1D *hist_delphi_btag_toptag;
   
   TH1D *hist_njetAK4;
   TH1D *hist_nbjetAK4;
   TH1D *hist_jetptAK4;
   TH1D *hist_jetrapAK4;
   TH1D *hist_jetmassAK4;
   TH1D *hist_jetsdmassAK4;
   TH1D *hist_jetbtagAK4;
   TH1D *hist_jetbtagdeepCSVAK4;
   TH1D *hist_jetbtagdeepflavAK4;
   TH1D *hist_jetpartonflavAK4;
   TH1D *hist_jetpartonflavAK4_btagged;
   TH1D *hist_bjetptAK4;
   TH2D *hist_tjetptAK8_bjetptAK4;
   
   TH1D *hist_jetmassAK4_1;
   TH1D *hist_jetbtagAK4_1;
   TH1D *hist_jetbtagdeepCSVAK4_1;
   TH1D *hist_jetbtagdeepflavAK4_1;
   TH1D *hist_jetpartonflavAK4_1;
   TH1D *hist_jetpartonflavAK4_btagged_1;
   
   TH1D *hist_biso_mass;
   TH1D *hist_biso_tbmassdiff;
   TH1D *hist_biso_isomass;
   TH1D *hist_biso_isopt;
   TH1D *hist_biso_TopAK8score;
   TH1D *hist_biso_TopAK8score_MD;
   TH1D *hist_biso_WAK8score;
   TH1D *hist_biso_WAK8score_MD;
   
   TH2D *hist_2d_biso_AK4mass_pt[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8sdmass_pt[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8TopScore_pt[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8TopScore_MD_pt[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8WScore_pt[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8WScore_MD_pt[ntoptag][nsbtag];
   
   TH2D *hist_2d_biso_AK4mass_pt_dAK8[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8sdmass_pt_dAK8[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8TopScore_pt_dAK8[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8TopScore_MD_pt_dAK8[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8WScore_pt_dAK8[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8WScore_MD_pt_dAK8[ntoptag][nsbtag];
   
   TH2D *hist_2d_biso_AK4mass_pt_mdAK8[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8sdmass_pt_mdAK8[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8TopScore_pt_mdAK8[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8TopScore_MD_pt_mdAK8[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8WScore_pt_mdAK8[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8WScore_MD_pt_mdAK8[ntoptag][nsbtag];
   
   // high b pt //
   
   TH1D *hist_biso_mass_1;
   TH1D *hist_biso_tbmassdiff_1;
   TH1D *hist_biso_isomass_1;
   TH1D *hist_biso_isopt_1;
   TH1D *hist_biso_TopAK8score_1;
   TH1D *hist_biso_TopAK8score_MD_1;
   TH1D *hist_biso_WAK8score_1;
   TH1D *hist_biso_WAK8score_MD_1;
   
   TH1D *hist_biso_mass_btag_1[nbtag];
   TH1D *hist_biso_TopAK8score_MD_btag_1[nbtag];
   TH1D *hist_biso_WAK8score_MD_btag_1[nbtag];
   
   TH2D *hist_2d_biso_AK4mass_pt_1[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8sdmass_pt_1[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8TopScore_pt_1[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8TopScore_MD_pt_1[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8WScore_pt_1[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8WScore_MD_pt_1[ntoptag][nsbtag];
   
   TH2D *hist_2d_biso_AK4mass_pt_dAK8_1[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8sdmass_pt_dAK8_1[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8TopScore_pt_dAK8_1[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8TopScore_MD_pt_dAK8_1[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8WScore_pt_dAK8_1[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8WScore_MD_pt_dAK8_1[ntoptag][nsbtag];
   
   TH2D *hist_2d_biso_AK4mass_pt_mdAK8_1[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8sdmass_pt_mdAK8_1[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8TopScore_pt_mdAK8_1[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8TopScore_MD_pt_mdAK8_1[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8WScore_pt_mdAK8_1[ntoptag][nsbtag];
   TH2D *hist_2d_biso_AK8WScore_MD_pt_mdAK8_1[ntoptag][nsbtag];
   
   TH2D *hist_2D_biso_AK8jet_subjetmass12;
   TH2D *hist_2D_biso_AK8jetsdmass_subjetmass1;
   TH2D *hist_2D_biso_AK8jetsdmass_subjetmass2;
   
   TH1D *hist_bjetAK8genmass;
   TH1D *hist_bjetAK8genmass_1;
   TH2D *hist_genmass_tbcor;
   
   // high b pt end //
   
   TH1D *hist_tbmass[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbpt[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbrap[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbpt_frac[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbtheta_diff[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbphi_diff[ntoptag][nbtag][nsbtag];
   TH1D *hist_partonflav_AK4[ntoptag][nbtag][nsbtag];
   
   TH1D *hist_topjetpt_sel[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetsdmass_sel[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetdeeptopscore_sel[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetdeepmdtopscore_sel[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetsbtag_sel[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetsdeepCSV_sel[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjettau32_sel[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetpt_sel[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetmass_sel[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetbtag_sel[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetmatch_AK8Topscore_sel[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetmatch_AK8Wscore_sel[ntoptag][nbtag][nsbtag];
   
   TH1D *hist_ht_tau32[ntoptag][nbtag][nsbtag];
   
   TH3D *hist_3D_topjetsdmass_topjettau32_tbmass[ntoptag][nbtag][nsbtag];
   
   // low b pt //
   
   TH1D *hist_tbmass_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbpt_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbrap_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbpt_frac_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbtheta_diff_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbphi_diff_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_partonflav_AK4_1[ntoptag][nbtag][nsbtag];
   
   TH1D *hist_topjetpt_sel_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetsdmass_sel_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetdeeptopscore_sel_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetdeepmdtopscore_sel_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetsbtag_sel_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetsdeepCSV_sel_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjettau32_sel_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetpt_sel_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetmass_sel_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetbtag_sel_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetmatch_AK8Topscore_sel_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetmatch_AK8Wscore_sel_1[ntoptag][nbtag][nsbtag];
   
   // low b pt end //
   
   // no bmass cut //
    
   TH1D *hist_tbmass_2[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbpt_2[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbrap_2[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbpt_frac_2[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbtheta_diff_2[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbphi_diff_2[ntoptag][nbtag][nsbtag];
   TH1D *hist_partonflav_AK4_2[ntoptag][nbtag][nsbtag];
   
   TH1D *hist_topjetpt_sel_2[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetsdmass_sel_2[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetdeeptopscore_sel_2[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetdeepmdtopscore_sel_2[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetsbtag_sel_2[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetsdeepCSV_sel_2[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjettau32_sel_2[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetpt_sel_2[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetmass_sel_2[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetbtag_sel_2[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetmatch_AK8Topscore_sel_2[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetmatch_AK8Wscore_sel_2[ntoptag][nbtag][nsbtag];
    
    // no bmass cut end //
    
   // high ht //
    
   TH1D *hist_tbmass_3[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbpt_3[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbrap_3[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbpt_frac_3[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbtheta_diff_3[ntoptag][nbtag][nsbtag];
   TH1D *hist_tbphi_diff_3[ntoptag][nbtag][nsbtag];
   TH1D *hist_partonflav_AK4_3[ntoptag][nbtag][nsbtag];
   
   TH1D *hist_topjetpt_sel_3[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetsdmass_sel_3[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetdeeptopscore_sel_3[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetdeepmdtopscore_sel_3[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetsbtag_sel_3[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjetsdeepCSV_sel_3[ntoptag][nbtag][nsbtag];
   TH1D *hist_topjettau32_sel_3[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetpt_sel_3[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetmass_sel_3[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetbtag_sel_3[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetmatch_AK8Topscore_sel_3[ntoptag][nbtag][nsbtag];
   TH1D *hist_bjetmatch_AK8Wscore_sel_3[ntoptag][nbtag][nsbtag];
   
   TH1D *hist_ht_tau32_3[ntoptag][nbtag][nsbtag];
    
    
    // high ht end //
    
    
    // deepAK8 based SR //
    
   TH1D *hist_tbmass_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbpt_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbrap_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_partonflav_AK4_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_topjetpt_sel_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetpt_msdbin_sel_DeepAK8[ntoptag_DAK8][nbtag][nsbtag][ntopsdmassbins];
   TH1D *hist_topjetsdmass_sel_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetdeeptopscore_sel_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetdeepmdtopscore_sel_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsbtag_sel_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsdeepCSV_sel_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjettau32_sel_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetpt_sel_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetmass_sel_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetbtag_sel_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetmatch_AK8Topscore_sel_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetmatch_AK8Wscore_sel_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_ht_DeepAK8[ntoptag_DAK8][nbtag][nsbtag]; 
   
   // Mass decorrelated
   
   TH1D *hist_tbmass_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbpt_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbrap_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_partonflav_AK4_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_dY[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_pfup[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_pfdn[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_dAK8up[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_dAK8dn[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_puup[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_pudn[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_prefireup[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_prefiredn[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_bcorup[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_bcordn[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_noptw[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_ptbin[ntoptag_DAK8][nbtag][nsbtag][ntopptbins];
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_ptbin[nbtag][nsbtag][ntopptbins];
   
   TH1D *hist_tbmass_md_DeepAK8_ptbin_tt[ntoptag_DAK8][nbtag][nsbtag][ntopptbins];
   TH2D *hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_ptbin_tt[nbtag][nsbtag][ntopptbins];
   
   TH1D *hist_topjetpt_sel_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetpt_msdbin_sel_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag][ntopsdmassbins];
   TH1D *hist_topjetsdmass_sel_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetdeeptopscore_sel_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetdeepmdtopscore_sel_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsbtag_sel_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsdeepCSV_sel_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjettau32_sel_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_partonflav_AK8_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetpt_sel_md_DeepAK8_beta[ntoptag_DAK8][nbtag][nsbtag][netabins];
   TH1D *hist_bjetpt_sel_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetmass_sel_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetbtag_sel_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetmatch_AK8Topscore_sel_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetmatch_AK8Wscore_sel_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetmatch_sdmass_sel_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag]; 
   
   TH1D *hist_topjetdeepmdtopscore_sel_md_DeepAK8up[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetdeepmdtopscore_sel_md_DeepAK8dn[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_bjetpt_sel_md_DeepAK8_Btag_up[ntoptag_DAK8][nbtag][nsbtag][nbtagmax];
   TH1D *hist_bjetpt_sel_md_DeepAK8_Btag_dn[ntoptag_DAK8][nbtag][nsbtag][nbtagmax];
   
   TH1D *hist_bjetpt_sel_md_DeepAK8_JES_up[ntoptag_DAK8][nbtag][nsbtag][njesmax];
   TH1D *hist_bjetpt_sel_md_DeepAK8_JES_dn[ntoptag_DAK8][nbtag][nsbtag][njesmax];
   
   TH2D *hist_tbmass_topjetpt_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH2D *hist_tbmass_topjeteta_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH2D *hist_tbmass_bjetpt_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   TH2D *hist_tbmass_bjeteta_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_ht_md_DeepAK8[ntoptag_DAK8][nbtag][nsbtag]; 
   
   TH1D *hist_tbmass_md_DeepAK8_PDF[ntoptag_DAK8][nbtag][nsbtag][npdfmax];
   TH1D *hist_topjetsdmass_sel_md_DeepAK8_PDF[ntoptag_DAK8][nbtag][nsbtag][npdfmax];
   TH1D *hist_tbmass_md_DeepAK8_Scale[ntoptag_DAK8][nbtag][nsbtag][nscalemax];
   TH1D *hist_topjetsdmass_sel_md_DeepAK8_Scale[ntoptag_DAK8][nbtag][nsbtag][nscalemax];
   TH1D *hist_tbmass_md_DeepAK8_PS[ntoptag_DAK8][nbtag][nsbtag][npsmax];
   TH1D *hist_topjetsdmass_sel_md_DeepAK8_PS[ntoptag_DAK8][nbtag][nsbtag][npsmax];
   
   TH1D *hist_tbmass_md_DeepAK8_Btag_up[ntoptag_DAK8][nbtag][nsbtag][nbtagmax];
   TH1D *hist_topjetsdmass_sel_md_DeepAK8_Btag_up[ntoptag_DAK8][nbtag][nsbtag][nbtagmax];
   TH1D *hist_tbmass_md_DeepAK8_Btag_dn[ntoptag_DAK8][nbtag][nsbtag][nbtagmax];
   TH1D *hist_topjetsdmass_sel_md_DeepAK8_Btag_dn[ntoptag_DAK8][nbtag][nsbtag][nbtagmax];
   
   TH1D *hist_tbmass_md_DeepAK8_JER_up[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsdmass_sel_md_DeepAK8_JER_up[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_JER_dn[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsdmass_sel_md_DeepAK8_JER_dn[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_JMR_up[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsdmass_sel_md_DeepAK8_JMR_up[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_JMR_dn[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsdmass_sel_md_DeepAK8_JMR_dn[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_topjetsdmass_sel_md_DeepAK8up[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsdmass_sel_md_DeepAK8dn[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_JES_up[ntoptag_DAK8][nbtag][nsbtag][njesmax];
   TH1D *hist_topjetsdmass_sel_md_DeepAK8_JES_up[ntoptag_DAK8][nbtag][nsbtag][njesmax];
   TH1D *hist_tbmass_md_DeepAK8_JES_dn[ntoptag_DAK8][nbtag][nsbtag][njesmax];
   TH1D *hist_topjetsdmass_sel_md_DeepAK8_JES_dn[ntoptag_DAK8][nbtag][nsbtag][njesmax];
   
   TH1D *hist_topjetpt_sel_md_DeepAK8_JES_up[ntoptag_DAK8][nbtag][nsbtag][njesmax];
   TH1D *hist_topjetpt_sel_md_DeepAK8_JES_dn[ntoptag_DAK8][nbtag][nsbtag][njesmax];
   
   TH1D *hist_tbmass_md_DeepAK8_JMS_up[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsdmass_sel_md_DeepAK8_JMS_up[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_JMS_dn[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsdmass_sel_md_DeepAK8_JMS_dn[ntoptag_DAK8][nbtag][nsbtag];
    
   TH1D *hist_tbmass_md_DeepAK8_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbpt_md_DeepAK8_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbrap_md_DeepAK8_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_partonflav_AK4_md_DeepAK8_3[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_topjetpt_sel_md_DeepAK8_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsdmass_sel_md_DeepAK8_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetdeeptopscore_sel_md_DeepAK8_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetdeepmdtopscore_sel_md_DeepAK8_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsbtag_sel_md_DeepAK8_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsdeepCSV_sel_md_DeepAK8_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjettau32_sel_md_DeepAK8_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetpt_sel_md_DeepAK8_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetmass_sel_md_DeepAK8_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetbtag_sel_md_DeepAK8_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetmatch_AK8Topscore_sel_md_DeepAK8_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetmatch_AK8Wscore_sel_md_DeepAK8_3[ntoptag_DAK8][nbtag][nsbtag]; 
   TH1D *hist_bjetmatch_sdmass_sel_md_DeepAK8_3[ntoptag_DAK8][nbtag][nsbtag]; 
   
   TH1D *hist_ht_md_DeepAK8_3[ntoptag_DAK8][nbtag][nsbtag]; 
   
   // for t-tbar //
   
   TH1D *hist_tbmass_md_DeepAK8_tt[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbpt_md_DeepAK8_tt[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbrap_md_DeepAK8_tt[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_partonflav_AK4_md_DeepAK8_tt[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_tt_deepak8pass[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_topjetpt_sel_md_DeepAK8_tt[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetpt_msdbin_sel_md_DeepAK8_tt[ntoptag_DAK8][nbtag][nsbtag][ntopsdmassbins];
   TH1D *hist_topjetsdmass_sel_md_DeepAK8_tt[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetdeeptopscore_sel_md_DeepAK8_tt[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetdeepmdtopscore_sel_md_DeepAK8_tt[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsbtag_sel_md_DeepAK8_tt[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsdeepCSV_sel_md_DeepAK8_tt[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjettau32_sel_md_DeepAK8_tt[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetpt_sel_md_DeepAK8_tt[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetmass_sel_md_DeepAK8_tt[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetbtag_sel_md_DeepAK8_tt[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetmatch_AK8Topscore_sel_md_DeepAK8_tt[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetmatch_AK8Wscore_sel_md_DeepAK8_tt[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetmatch_sdmass_sel_md_DeepAK8_tt[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_topjetpt_sel_md_DeepAK8_tt_deepak8pass[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsdmass_sel_md_DeepAK8_tt_deepak8pass[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_ht_md_DeepAK8_tt[ntoptag_DAK8][nbtag][nsbtag]; 
    
   TH1D *hist_tbmass_md_DeepAK8_tt_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbpt_md_DeepAK8_tt_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbrap_md_DeepAK8_tt_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_partonflav_AK4_md_DeepAK8_tt_3[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_topjetpt_sel_md_DeepAK8_tt_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsdmass_sel_md_DeepAK8_tt_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetdeeptopscore_sel_md_DeepAK8_tt_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetdeepmdtopscore_sel_md_DeepAK8_tt_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsbtag_sel_md_DeepAK8_tt_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjetsdeepCSV_sel_md_DeepAK8_tt_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_topjettau32_sel_md_DeepAK8_tt_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetpt_sel_md_DeepAK8_tt_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetmass_sel_md_DeepAK8_tt_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetbtag_sel_md_DeepAK8_tt_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetmatch_AK8Topscore_sel_md_DeepAK8_tt_3[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_bjetmatch_AK8Wscore_sel_md_DeepAK8_tt_3[ntoptag_DAK8][nbtag][nsbtag]; 
   TH1D *hist_bjetmatch_sdmass_sel_md_DeepAK8_tt_3[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_ht_md_DeepAK8_tt_3[ntoptag_DAK8][nbtag][nsbtag]; 
    
    // DeepAK8 SR ends //
     
    
 // with   AK8 b cand
 
   TH1D *hist_tbmass_AK8[ntoptag][nbtag][nsbtag];
   TH1D *hist_AK8bCand_DeeptagScore[ntoptag][nbtag][nsbtag];
   TH1D *hist_AK8bCand_DeeptagScore_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_AK8bCand_MD_DeeptagScore[ntoptag][nbtag][nsbtag];
   
   TH1D *hist_tbmass_AK8_DeepAK8[ntoptag][nbtag][nsbtag];
   TH1D *hist_AK8bCand_DeeptagScore_DeepAK8[ntoptag][nbtag][nsbtag];
   TH1D *hist_AK8bCand_DeeptagScore_DeepAK8_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_AK8bCand_MD_DeeptagScore_DeepAK8[ntoptag][nbtag][nsbtag];
   
   TH1D *hist_tbmass_AK8_MDDeepAK8[ntoptag][nbtag][nsbtag];
   TH1D *hist_AK8bCand_DeeptagScore_MDDeepAK8[ntoptag][nbtag][nsbtag];
   TH1D *hist_AK8bCand_DeeptagScore_MDDeepAK8_1[ntoptag][nbtag][nsbtag];
   TH1D *hist_AK8bCand_MD_DeeptagScore_MDDeepAK8[ntoptag][nbtag][nsbtag];
   
    
    // trigger histograms //
    
   static const int trig_msdbins = 3;
   double trig_msd_vals[trig_msdbins+1] = {0,105,210,4000};
    
   TH1D *hist_ht_mutrig[trig_msdbins];
   TH1D *hist_pt_mutrig[trig_msdbins];
   TH1D *hist_ht_emutrig[trig_msdbins];
   TH1D *hist_pt_emutrig[trig_msdbins];
   TH1D *hist_bpt_emutrig[trig_msdbins];
   TH1D *hist_ptsum_emutrig[trig_msdbins];
   TH1D *hist_mtb_emutrig[trig_msdbins];
   
   TH1D *hist_ht_HT1050_wemutrig[trig_msdbins];
   TH1D *hist_ht_AK4Pt500_wemutrig[trig_msdbins];
   TH1D *hist_ht_AK8Pt500_wemutrig[trig_msdbins];
   TH1D *hist_ht_AK8PFJet420_TrimMass30_wemutrig[trig_msdbins];
   TH1D *hist_ht_AK8DiPFJet300_200_TrimMass30_wemutrig[trig_msdbins]; 
   TH1D *hist_ht_AK8PFHT900_TrimMass50_wemutrig[trig_msdbins]; 
   TH1D *hist_ht_all_wemutrig[trig_msdbins];
   
   TH1D *hist_pt_HT1050_wemutrig[trig_msdbins];
   TH1D *hist_pt_AK4Pt500_wemutrig[trig_msdbins];
   TH1D *hist_pt_AK8Pt500_wemutrig[trig_msdbins];
   TH1D *hist_pt_AK8PFJet420_TrimMass30_wemutrig[trig_msdbins];
   TH1D *hist_pt_AK8DiPFJet300_200_TrimMass30_wemutrig[trig_msdbins]; 
   TH1D *hist_pt_AK8PFHT900_TrimMass50_wemutrig[trig_msdbins]; 
   TH1D *hist_pt_all_wemutrig[trig_msdbins];
   
   TH1D *hist_bpt_HT1050_wemutrig[trig_msdbins];
   TH1D *hist_bpt_AK4Pt500_wemutrig[trig_msdbins];
   TH1D *hist_bpt_AK8Pt500_wemutrig[trig_msdbins];
   TH1D *hist_bpt_AK8PFJet420_TrimMass30_wemutrig[trig_msdbins];
   TH1D *hist_bpt_AK8DiPFJet300_200_TrimMass30_wemutrig[trig_msdbins]; 
   TH1D *hist_bpt_AK8PFHT900_TrimMass50_wemutrig[trig_msdbins]; 
   TH1D *hist_bpt_all_wemutrig[trig_msdbins];
   
   TH1D *hist_ptsum_HT1050_wemutrig[trig_msdbins];
   TH1D *hist_ptsum_AK4Pt500_wemutrig[trig_msdbins];
   TH1D *hist_ptsum_AK8Pt500_wemutrig[trig_msdbins];
   TH1D *hist_ptsum_AK8PFJet420_TrimMass30_wemutrig[trig_msdbins];
   TH1D *hist_ptsum_AK8DiPFJet300_200_TrimMass30_wemutrig[trig_msdbins]; 
   TH1D *hist_ptsum_AK8PFHT900_TrimMass50_wemutrig[trig_msdbins]; 
   TH1D *hist_ptsum_all_wemutrig[trig_msdbins];
   
   TH1D *hist_mtb_HT1050_wemutrig[trig_msdbins];
   TH1D *hist_mtb_AK4Pt500_wemutrig[trig_msdbins];
   TH1D *hist_mtb_AK8Pt500_wemutrig[trig_msdbins];
   TH1D *hist_mtb_AK8PFJet420_TrimMass30_wemutrig[trig_msdbins];
   TH1D *hist_mtb_AK8DiPFJet300_200_TrimMass30_wemutrig[trig_msdbins]; 
   TH1D *hist_mtb_AK8PFHT900_TrimMass50_wemutrig[trig_msdbins]; 
   TH1D *hist_mtb_all_wemutrig[trig_msdbins];
   
   TF1 *bpfrat;
   TF1 *bpfrat_beta[netabins];
   TF1 *bpfrat_val_beta[netabins];
   TF1 *tagpfrat_tau32;
   TF1 *tagpfrat_MDAK8;
   
   TMatrixD err1;
   TMatrixD err;
   
   TF1 *fun_bdF[5][netabins];
   TMatrixD cov_bpfrat_beta[netabins];
   TF1 *fun_val_bdF[5][netabins];
   TMatrixD cov_bpfrat_val_beta[netabins];
   TMatrixD rvar;//(1,5);
   TMatrixD cvar;//(5,1);

float COV_sig_data[netabins][5][5] = {{{0}}};
float COV_sig_mc[netabins][5][5] = {{{0}}};
float COV_val_data[netabins][5][5] = {{{0}}};
float COV_val_mc[netabins][5][5] = {{{0}}};

//2016//

float COV_sig_data_16[netabins][5][5] = {
{
{1.90365e-12,-5.60411e-10,1.0794e-11,4.08474e-14,-9.19802e-15},
{-5.60411e-10,2.35804e-06,5.37058e-09,-8.32844e-12,-7.99846e-12},
{1.0794e-11,5.37058e-09,1.70353e-10,4.96898e-13,-1.66966e-13},
{4.08474e-14,-8.32844e-12,4.96898e-13,1.81698e-15,-4.56977e-16},
{-9.19802e-15,-7.99846e-12,-1.66966e-13,-4.56977e-16,2.77832e-16}
},
{
{1.56149e-12,-3.27474e-10,1.07444e-11,3.78193e-14,-1.33438e-14},
{-3.27474e-10,2.41218e-06,9.08781e-09,3.40452e-12,-1.75531e-11},
{1.07444e-11,9.08781e-09,2.08209e-10,5.7667e-13,-2.98237e-13},
{3.78193e-14,3.40452e-12,5.7667e-13,1.91535e-15,-7.76927e-16},
{-1.33438e-14,-1.75531e-11,-2.98237e-13,-7.76927e-16,6.58914e-16}
},
{
{2.61464e-12,-1.48967e-09,1.50529e-11,6.01234e-14,-2.11089e-14},
{-1.48967e-09,1.98124e-05,1.45848e-07,2.58697e-10,-4.93429e-10},
{1.50529e-11,1.45848e-07,1.84226e-09,4.32359e-12,-5.48712e-12},
{6.01234e-14,2.58697e-10,4.32359e-12,1.11654e-14,-1.228e-14},
{-2.11089e-14,-4.93429e-10,-5.48712e-12,-1.228e-14,2.96372e-14}
}
};


float COV_sig_mc_16[netabins][5][5] = {
{
{3.12376e-12,-7.09427e-10,5.18341e-12,2.45967e-14,-3.5707e-15},
{-7.09427e-10,1.20988e-06,-2.93397e-09,-2.46259e-11,1.3248e-12},
{5.18341e-12,-2.93397e-09,5.14641e-11,2.09027e-13,-4.37786e-14},
{2.45967e-14,-2.46259e-11,2.09027e-13,1.3638e-15,-1.63539e-16},
{-3.5707e-15,1.3248e-12,-4.37786e-14,-1.63539e-16,5.26334e-17}
},
{
{2.71499e-12,-7.07661e-10,6.55525e-12,3.06755e-14,-5.5373e-15},
{-7.07661e-10,9.87472e-07,-2.40658e-09,-2.01196e-11,9.70884e-13},
{6.55525e-12,-2.40658e-09,4.97013e-11,1.96924e-13,-4.92664e-14},
{3.06755e-14,-2.01196e-11,1.96924e-13,1.15679e-15,-1.78175e-16},
{-5.5373e-15,9.70884e-13,-4.92664e-14,-1.78175e-16,7.47354e-17}
},
{
{1.29935e-11,-2.3161e-09,6.14118e-11,2.37031e-13,-1.05931e-13},
{-2.3161e-09,4.46875e-06,8.24841e-10,-5.01571e-11,-2.24452e-11},
{6.14118e-11,8.24841e-10,5.34399e-10,1.79239e-12,-1.03058e-12},
{2.37031e-13,-5.01571e-11,1.79239e-12,7.28751e-15,-3.20959e-15},
{-1.05931e-13,-2.24452e-11,-1.03058e-12,-3.20959e-15,2.95643e-15}
}
};
  
float COV_val_data_16[netabins][5][5] = {
{
{5.53301e-14,-2.3789e-11,3.86824e-13,1.42897e-15,-3.4933e-16},
{-2.3789e-11,2.26073e-07,1.21124e-09,1.52138e-12,-1.55655e-12},
{3.86824e-13,1.21124e-09,1.8149e-11,4.50393e-14,-2.01995e-14},
{1.42897e-15,1.52138e-12,4.50393e-14,1.30241e-16,-4.79939e-17},
{-3.4933e-16,-1.55655e-12,-2.01995e-14,-4.79939e-17,4.25063e-17}
},
{
{4.8839e-14,-2.18316e-11,3.68685e-13,1.34499e-15,-4.27524e-16},
{-2.18316e-11,2.73797e-07,1.69049e-09,2.55602e-12,-2.9102e-12},
{3.68685e-13,1.69049e-09,2.28975e-11,5.4986e-14,-3.44792e-14},
{1.34499e-15,2.55602e-12,5.4986e-14,1.50084e-16,-7.91229e-17},
{-4.27524e-16,-2.9102e-12,-3.44792e-14,-7.91229e-17,1.01428e-16}
},
{
{1.71824e-13,-4.33441e-10,-2.70951e-12,-4.31189e-15,9.39154e-15},
{-4.33441e-10,2.28226e-06,1.73783e-08,3.15555e-11,-5.3492e-11},
{-2.70951e-12,1.73783e-08,1.83054e-10,3.96809e-13,-5.18427e-13},
{-4.31189e-15,3.15555e-11,3.96809e-13,9.27418e-16,-1.08316e-15},
{9.39154e-15,-5.3492e-11,-5.18427e-13,-1.08316e-15,3.99051e-15}
}
};

float COV_val_mc_16[netabins][5][5] = {
{
{1.68928e-13,-4.94081e-11,5.93937e-13,2.5815e-15,-5.16047e-16},
{-4.94081e-11,7.6614e-08,-1.04168e-10,-1.23977e-12,-6.79999e-15},
{5.93937e-13,-1.04168e-10,4.78222e-12,1.71832e-14,-4.72082e-15},
{2.5815e-15,-1.23977e-12,1.71832e-14,8.241e-17,-1.5687e-17},
{-5.16047e-16,-6.79999e-15,-4.72082e-15,-1.5687e-17,7.36793e-18}
},
{
{1.55051e-13,-4.55088e-11,6.23514e-13,2.62356e-15,-6.33146e-16},
{-4.55088e-11,7.23527e-08,-7.94273e-11,-1.09505e-12,-4.56924e-14},
{6.23514e-13,-7.94273e-11,5.19521e-12,1.82044e-14,-5.98585e-15},
{2.62356e-15,-1.09505e-12,1.82044e-14,8.14946e-17,-1.93977e-17},
{-6.33146e-16,-4.56924e-14,-5.98585e-15,-1.93977e-17,1.08286e-17}
},
{
{6.50713e-13,-6.86082e-11,4.24102e-12,1.47513e-14,-8.16368e-15},
{-6.86082e-11,3.95269e-07,1.27163e-09,-2.32367e-13,-4.77155e-12},
{4.24102e-12,1.27163e-09,5.00069e-11,1.50541e-13,-1.09172e-13},
{1.47513e-14,-2.32367e-13,1.50541e-13,5.24996e-16,-3.06482e-16},
{-8.16368e-15,-4.77155e-12,-1.09172e-13,-3.06482e-16,4.02461e-16}
}
};


//2017//
  
float COV_sig_data_17[netabins][5][5] = {
{
{14949.8,-0.147481,0.000377558,-6.57856e-06,-1.18034e-06},
{-0.147481,2.18693e-06,-2.47772e-09,6.07337e-11,8.49713e-12},
{0.000377558,-2.47772e-09,6.87188e-11,1.52528e-14,-1.09812e-13},
{-6.57856e-06,6.07337e-11,1.52528e-14,3.59065e-15,2.90738e-16},
{-1.18034e-06,8.49713e-12,-1.09812e-13,2.90738e-16,2.42781e-16}
},
{
{12140.6,-0.0725891,0.00014075,-4.15582e-06,-8.95973e-07},
{-0.0725891,8.7184e-07,3.88324e-11,2.28191e-11,3.35877e-12},
{0.00014075,3.88324e-11,3.7447e-11,6.00268e-14,-5.78547e-14},
{-4.15582e-06,2.28191e-11,6.00268e-14,1.82924e-15,1.72655e-16},
{-8.95973e-07,3.35877e-12,-5.78547e-14,1.72655e-16,1.5534e-16}
},
{
{2.81391e-07,6.35327e-316,6.35327e-316,1.63042e-322,6.90916e-310},
{2.122e-314,6.35592e-316,1.63042e-322,6.34837e-316,2000},
{6.35326e-316,3.21143e-322,6.90916e-310,6.35605e-316,6.35605e-316},
{6.35605e-316,6.356e-316,6.356e-316,6.356e-316,2.42092e-322},
{6.90917e-310,6.35592e-316,6.3465e-316,6.35521e-316,9.99276e-307}
}
};


float COV_sig_mc_17[netabins][5][5] = {
{
{92835.4,-1.8746,0.00181286,-2.50205e-06,-4.94541e-07},
{-1.8746,3.8253e-05,-3.78392e-08,4.1528e-11,1.066e-11},
{0.00181286,-3.78392e-08,4.89547e-11,1.18255e-14,-2.15493e-14},
{-2.50205e-06,4.1528e-11,1.18255e-14,5.12684e-16,-3.49912e-17},
{-4.94541e-07,1.066e-11,-2.15493e-14,-3.49912e-17,1.78466e-17}
},
{
{36030,-0.368063,0.000211904,-7.74262e-07,-6.94626e-07},
{-0.368063,4.38866e-06,3.98222e-12,9.6651e-12,4.94692e-12},
{0.000211904,3.98222e-12,1.66563e-11,1.3943e-14,-1.78408e-14},
{-7.74262e-07,9.6651e-12,1.3943e-14,4.24869e-17,-1.02544e-18},
{-6.94626e-07,4.94692e-12,-1.78408e-14,-1.02544e-18,2.70258e-17}
},
{
{36030,-0.368063,0.000211904,-7.74262e-07,-6.94626e-07},
{-0.368063,4.38866e-06,3.98222e-12,9.6651e-12,4.94692e-12},
{0.000211904,3.98222e-12,1.66563e-11,1.3943e-14,-1.78408e-14},
{-7.74262e-07,9.6651e-12,1.3943e-14,4.24869e-17,-1.02544e-18},
{-6.94626e-07,4.94692e-12,-1.78408e-14,-1.02544e-18,2.70258e-17}
}
};
  
float COV_val_data_17[netabins][5][5] = {
{
{4781.01,0.0346566,2.21963e-05,-7.71602e-07,-1.53584e-07},
{0.0346566,2.99385e-07,3.39735e-10,-5.53191e-12,-1.37885e-12},
{2.21963e-05,3.39735e-10,3.92655e-12,6.85265e-15,-5.10509e-15},
{-7.71602e-07,-5.53191e-12,6.85265e-15,1.59053e-16,1.34449e-17},
{-1.53584e-07,-1.37885e-12,-5.10509e-15,1.34449e-17,1.29224e-17}
},
{
{13543.8,0.169079,4.84842e-05,-2.57118e-06,-1.03509e-06},
{0.169079,2.14745e-06,8.07788e-10,-3.18365e-11,-1.33417e-11},
{4.84842e-05,8.07788e-10,3.51209e-12,-6.28621e-16,-9.43776e-15},
{-2.57118e-06,-3.18365e-11,-6.28621e-16,5.13683e-16,1.82542e-16},
{-1.03509e-06,-1.33417e-11,-9.43776e-15,1.82542e-16,9.49192e-17}
},
{
{18351.7,0.344317,6.76894e-05,-2.20693e-06,-2.67943e-06},
{0.344317,6.69884e-06,3.06377e-09,-3.81792e-11,-5.61477e-11},
{6.76894e-05,3.06377e-09,2.1321e-11,3.99337e-14,-7.17205e-14},
{-2.20693e-06,-3.81792e-11,3.99337e-14,3.85487e-16,1.87291e-16},
{-2.67943e-06,-5.61477e-11,-7.17205e-14,1.87291e-16,7.80599e-16}
}
};

float COV_val_mc_17[netabins][5][5] = {
{
{9873.22,0.0486991,-7.60813e-05,-1.48657e-06,-7.58036e-08},
{0.0486991,2.59885e-07,-4.15282e-10,-7.69307e-12,-3.58267e-13},
{-7.60813e-05,-4.15282e-10,1.57583e-12,1.52064e-14,-3.00497e-16},
{-1.48657e-06,-7.69307e-12,1.52064e-14,2.43901e-16,8.29308e-18},
{-7.58036e-08,-3.58267e-13,-3.00497e-16,8.29308e-18,1.64155e-18}
},
{
{12369.7,0.104873,-0.000220728,-2.31703e-06,-9.72895e-08},
{0.104873,9.02626e-07,-1.89892e-09,-1.9895e-11,-8.16865e-13},
{-0.000220728,-1.89892e-09,4.69949e-12,4.42189e-14,9.11004e-16},
{-2.31703e-06,-1.9895e-11,4.42189e-14,4.48744e-16,1.53648e-17},
{-9.72895e-08,-8.16865e-13,9.11004e-16,1.53648e-17,2.07866e-18}
},
{
{20952.4,0.218814,-0.000635207,-1.99722e-06,-9.24062e-08},
{0.218814,2.34568e-06,-6.55143e-09,-2.12728e-11,-1.41145e-12},
{-0.000635207,-6.55143e-09,2.47981e-11,7.78951e-14,-7.97034e-15},
{-1.99722e-06,-2.12728e-11,7.78951e-14,2.57641e-16,-2.16657e-17},
{-9.24062e-08,-1.41145e-12,-7.97034e-15,-2.16657e-17,4.03803e-17}
}
};  

//2018//

float COV_sig_data_18[netabins][5][5] = {
{
{21046.7,-0.567981,0.00122036,-2.123e-06,-1.1622e-06},
{-0.567981,1.5846e-05,-3.24532e-08,5.29434e-11,2.95139e-11},
{0.00122036,-3.24532e-08,1.11583e-10,7.13209e-15,-1.2428e-13},
{-2.123e-06,5.29434e-11,7.13209e-15,7.38734e-16,-5.07114e-17},
{-1.1622e-06,2.95139e-11,-1.2428e-13,-5.07114e-17,1.75927e-16}
},
{
{7974.5,-0.121119,0.000426532,-2.47602e-06,-5.76721e-07},
{-0.121119,2.15723e-06,-5.97009e-09,3.56441e-11,7.36305e-12},
{0.000426532,-5.97009e-09,4.876e-11,-5.19437e-14,-6.64741e-14},
{-2.47602e-06,3.56441e-11,-5.19437e-14,1.08287e-15,7.65771e-17},
{-5.76721e-07,7.36305e-12,-6.64741e-14,7.65771e-17,1.15481e-16}
},
{
{7974.5,-0.121119,0.000426532,-2.47602e-06,-5.76721e-07},
{-0.121119,2.15723e-06,-5.97009e-09,3.56441e-11,7.36305e-12},
{0.000426532,-5.97009e-09,4.876e-11,-5.19437e-14,-6.64741e-14},
{-2.47602e-06,3.56441e-11,-5.19437e-14,1.08287e-15,7.65771e-17},
{-5.76721e-07,7.36305e-12,-6.64741e-14,7.65771e-17,1.15481e-16}
}
};


float COV_sig_mc_18[netabins][5][5] = {
{
{158686,-3.46647,0.00224647,-4.29139e-05,-1.0743e-06},
{-3.46647,7.61768e-05,-5.02749e-08,9.28259e-10,2.40336e-11},
{0.00224647,-5.02749e-08,4.43018e-11,-5.52529e-13,-2.32129e-14},
{-4.29139e-05,9.28259e-10,-5.52529e-13,1.20318e-14,2.57669e-16},
{-1.0743e-06,2.40336e-11,-2.32129e-14,2.57669e-16,1.36656e-17}
},
{
{85657.8,-1.13208,-0.000683181,-1.91658e-05,-7.29642e-07},
{-1.13208,1.52541e-05,8.27301e-09,2.4732e-10,1.00509e-11},
{-0.000683181,8.27301e-09,1.66132e-11,1.9889e-13,-4.84301e-15},
{-1.91658e-05,2.4732e-10,1.9889e-13,4.6023e-15,1.23037e-16},
{-7.29642e-07,1.00509e-11,-4.84301e-15,1.23037e-16,1.91897e-17}
},
{
{85657.8,-1.13208,-0.000683181,-1.91658e-05,-7.29642e-07},
{-1.13208,1.52541e-05,8.27301e-09,2.4732e-10,1.00509e-11},
{-0.000683181,8.27301e-09,1.66132e-11,1.9889e-13,-4.84301e-15},
{-1.91658e-05,2.4732e-10,1.9889e-13,4.6023e-15,1.23037e-16},
{-7.29642e-07,1.00509e-11,-4.84301e-15,1.23037e-16,1.91897e-17}
}
};
  
float COV_val_data_18[netabins][5][5] = {
{
{5011,0.0287852,7.82044e-05,-7.70497e-07,-2.10777e-07},
{0.0287852,1.95696e-07,5.56129e-10,-4.40681e-12,-1.39061e-12},
{7.82044e-05,5.56129e-10,3.40344e-12,-6.13763e-15,-6.04556e-15},
{-7.70497e-07,-4.40681e-12,-6.13763e-15,1.38231e-16,2.54572e-17},
{-2.10777e-07,-1.39061e-12,-6.04556e-15,2.54572e-17,1.47877e-17}
},
{
{8234.19,0.0887495,2.22517e-05,-1.15258e-06,-4.03581e-07},
{0.0887495,9.79461e-07,3.55947e-10,-1.22925e-11,-4.56287e-12},
{2.22517e-05,3.55947e-10,2.08212e-12,2.1518e-15,-4.11973e-15},
{-1.15258e-06,-1.22925e-11,2.1518e-15,1.7742e-16,4.9016e-17},
{-4.03581e-07,-4.56287e-12,-4.11973e-15,4.9016e-17,2.74148e-17}
},
{
{6858.15,0.0978116,-4.36554e-05,-7.09131e-07,-8.28716e-07},
{0.0978116,1.5033e-06,-2.11122e-11,-9.31989e-12,-1.44245e-11},
{-4.36554e-05,-2.11122e-11,6.93024e-12,1.92369e-14,-1.81527e-14},
{-7.09131e-07,-9.31989e-12,1.92369e-14,1.12388e-16,3.86932e-17},
{-8.28716e-07,-1.44245e-11,-1.81527e-14,3.86932e-17,3.65731e-16}
}
};

float COV_val_mc_18[netabins][5][5] = {
{
{13662.4,0.0482629,-0.000127833,-1.23966e-06,-6.99785e-08},
{0.0482629,1.98797e-07,-5.24059e-10,-4.94001e-12,-2.10863e-13},
{-0.000127833,-5.24059e-10,2.35252e-12,1.62452e-14,-2.67409e-16},
{-1.23966e-06,-4.94001e-12,1.62452e-14,1.38394e-16,2.88805e-18},
{-6.99785e-08,-2.10863e-13,-2.67409e-16,2.88805e-18,1.29988e-18}
},
{
{13522.5,0.11481,-0.000357679,-1.6736e-06,-4.93379e-08},
{0.11481,9.94824e-07,-3.07045e-09,-1.45474e-11,-4.1283e-13},
{-0.000357679,-3.07045e-09,1.0413e-11,4.77812e-14,2.99155e-16},
{-1.6736e-06,-1.45474e-11,4.77812e-14,2.25146e-16,2.72915e-18},
{-4.93379e-08,-4.1283e-13,2.99155e-16,2.72915e-18,1.68166e-18}
},
{
{38603.3,0.584903,-0.00029867,-6.76321e-06,-1.82864e-06},
{0.584903,8.95583e-06,-4.30315e-09,-1.02804e-10,-2.85483e-11},
{-0.00029867,-4.30315e-09,1.36856e-11,8.7574e-14,-7.86966e-15},
{-6.76321e-06,-1.02804e-10,8.7574e-14,1.31329e-15,2.56839e-16},
{-1.82864e-06,-2.85483e-11,-7.86966e-15,2.56839e-16,1.59538e-16}
}
};

   TF1 *bpfrat_val;
//   TF1 *bpfrat_val_beta[netabins];
   TF1 *tagpfrat_tau32_val;
   TF1 *tagpfrat_MDAK8_val;
   TF1 *tagpfrat_MDAK8_val2;
   
   float lepptcut = 30;
    
   float HTcut = 900;
   float AK8ptcut = 507.1;//548.1;//507.1;
   float AK8ptcut_in = 200;
   float AK4ptcut_in = 40;
   float AK4ptcut_fi = 548.1;//548.1;
   float AK8subptcut = 30;
   float AK4masscut = 100;
   float jeteta_cut = 2.4;
   float tau32_cut = 0.54;
   float tau32_cut_loose = 10;//0.85;
   
   float rapidity_cut = 1.8;

   // 2016 //
//   float btagvalue = 0.8484; //0.9535; //tight 0.8484; // medium //csvv2 
//   float btagvalue_deepCSV = 0.6321; //0.8953  ;//tight  0.6321; //medium // deepcsv 
//   float btagvalue_deepFlavB = 0.45; //0.7221; // 0.3093;  // medium	//deepflav
	
   // 2017 //
//   float btagvalue = 0.8838;//0.9693;//tight 0.8838; //csvv2 medium
//   float btagvalue_deepCSV =  0.4941;//0.65;//myvalue 0.8001;// tight 0.4941; // deepcsv medium
//   float btagvalue_deepFlavB = 0.6; //0.6; // my final value //0.3033; //0.55;// my value 0.7489; // tight 0.3033; //deepflav medium

   // 2018 //
//   float btagvalue = 0.8838; 
//   float btagvalue_deepCSV = 0.4184; //0.7527  ;//tight  0.4184; //medium // deepcsv 
//   float btagvalue_deepFlavB = 0.6;  //0.7264; ;//tight   0.2770;  // medium	//deepflav

   float btagvalue ;
   float btagvalue_deepCSV;
   float btagvalue_deepFlavB;
   
   float dR_cut = 1.2;
   float dphi_cut = M_PI/2.;
   
   float btagvalue_deepCSV_m = 0.4941;
   float btagvalue_deepCSV_t = 0.8001;

   float tau32_cut_1 = 0.65;
   float tau32_cut_2 = 0.4;

//   float deepak8_cut = 0.895; // 2017 1% mistag rate //0.843;//0.391;
//   float deepak8_cut_md = 0.578; // 2017 0.5% mistag rate //0.843;//0.391;
//   float deepak8_cut = 0.937; // 2016 0.5% mistag rate
//   float deepak8_cut_md = 0.621; // 2016 0.5% mistag rate
//   float deepak8_cut = 0.898; // 2018 0.5% mistag rate
//   float deepak8_cut_md = 0.559; // 2018 0.5% mistag rate

  float deepak8_cut;
  float deepak8_cut_md;

   float deepak8_Wcut = 0.981;
   float deepak8_Wcut_md = 0.802;
   
   static const int noAK8ptbins = 4;
   float AK8ptbins[noAK8ptbins+1] = {300,400,480,600,10000};
   
   float minmasscut = 0;
   
   float topmasslow = 105;
   float topmasshigh = 210;
   
   float wmasslow = 60;
   float wmasshigh = 100;
   
   int nwpmbin = 60;
   float nwpmlow = 1000;
   float nwpmhigh = 4000;
   
   static const int nvwpmbin = 42;
   double wpmbins[nvwpmbin+1] = {1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,1650,1700,1750,1800,1850,1900,1950,2000,
					2050,2100,2150,2200,2300,2400,2500,2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,
					3800,4000,4500,5000};
   
   static const int nvwpmbin_trig = 45;
   double wpmbins_trig[nvwpmbin_trig+1] = {850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,1650,1700,1750,1800,1850,1900,1950,2000,
					2050,2100,2150,2200,2300,2400,2500,2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,
					3800,4000,4500,5000};
   
   static const int noptbins = 52 ;
   static const int nobptbins = 32;//9;
 
   static const int norhobins = 51;
   
   int nomassbins = 60;
   double mass_low = 0, mass_high = 360; 
   
   static const int nsdmassbins = 20;
   
   double ptbins[noptbins+1] = {30, 37, 43, 49, 56, 64, 74, 84,
     97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
     507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
     1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103};
	 // 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,
    // 4037, 4252, 4477, 5000} ;
   
//   double bptbins[nobptbins+1] = {20,30,50,70,100,140,200,300,600,1000}; 
   double bptbins[nobptbins+1] = {20, 30, 50, 74, 84,
     97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
     507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1172, 1497, 2000};
    
   float perf_ptbins[noperf_ptbins+1] = {300,500,700,1000,3000}; 
   
   float qcd_md_deepak8_eff[noperf_ptbins] = {0.027,0.027,0.035,0.045};
    
   float betarange[netarange+1] = {0,0.6,1.2,2.4};
   
   double pu_data[80] = {2.15411e-05,5.6992e-05,0.000127209,0.000284562,0.000628243,0.00130454,0.00248809,0.00434675,0.00699958,0.0104813,0.014724,0.0195614,0.0247524,0.0300163,0.0350706,0.0396636,0.0435963,0.0467326,0.0489994,0.0503793,0.0509,0.0506223,0.0496295,0.0480181,0.0458916,0.0433545,0.0405093,0.0374531,0.0342757,0.0310578,0.02787,0.0247719,0.0218122,0.0190285,0.0164478,0.0140877,0.011957,0.010057,0.00838287,0.0069247,0.00566889,0.00459925,0.00369801,0.00294673,0.00232704,0.0018212,0.00141253,0.00108574,0.000827064,0.000624371,0.000467131,0.000346366,0.000254531,0.000185382,0.000133825,9.57567e-05,6.79199e-05,4.77594e-05,3.3297e-05,2.30198e-05,1.57847e-05,1.07379e-05,7.2494e-06,4.85952e-06,3.23649e-06,2.14353e-06,1.41347e-06,9.29544e-07,6.11011e-07,4.02633e-07,2.66985e-07,1.78963e-07,1.21892e-07,8.47946e-08,6.05156e-08,4.4429e-08,3.35702e-08,2.6054e-08,2.06895e-08,1.67288e-08};

   double pu_MC[80] = {1.76003e-05,2.63752e-05,5.17764e-05,8.93703e-05,0.00010729,0.000139737,0.000241458,0.000725216,0.00129534,0.0024431,0.00503314,0.0092062,0.0146853,0.0204755,0.0267466,0.0337681,0.040181,0.0449554,0.0491183,0.0524606,0.0548204,0.0560606,0.0554619,0.0536665,0.0513695,0.0476493,0.0435073,0.0393303,0.0350598,0.0306642,0.0272167,0.0236864,0.0208034,0.0182563,0.0160933,0.0142631,0.012782,0.0115606,0.0105507,0.00957065,0.00885585,0.00826288,0.00758751,0.0069752,0.00622485,0.00547198,0.00483291,0.0040488,0.00337946,0.0027023,0.00212516,0.00160003,0.00117206,0.00085644,0.00056259,0.000367971,0.000246554,0.000160166,0.000101268,6.68301e-05,3.96657e-05,2.67149e-05,2.02967e-05,2.01132e-05,1.41779e-05,1.33955e-05,1.33216e-05,1.33078e-05,1.34937e-05,1.51125e-05,1.49794e-05,1.42329e-05,1.28154e-05,1.34488e-05,1.35528e-05,0,0,0,0,0};

//   double pu_rat17[100] = {0.000349139,0.0494244,0.0550583,0.0868762,0.0890973,0.127646,0.157756,0.197357,0.142307,0.450653,0.573772,0.724296,0.748964,0.753009,0.800193,0.849205,0.934488,1.02547,1.10455,1.16029,1.21521,1.26439,1.30349,1.32947,1.33878,1.34613,1.35251,1.35274,1.35446,1.34071,1.30643,1.25088,1.18394,1.11849,1.03882,0.973938,0.93495,0.89719,0.841632,0.804349,0.807898,0.848681,0.91417,1.02053,1.16536,1.33074,1.50073,1.57316,1.61508,1.52999,1.3599,1.18733,0.984448,0.769307,0.580557,0.420071,0.295849,0.203954,0.141305,0.0991562,0.0717071,0.0515937,0.0393459,0.0295183,0.0232182,0.0170967,0.0125657,0.0105301,0.00923461,0.00843191,0.00734179,0.00520827,0.00554619,0.00531438,0.0030063,0.00378597,0.00252698,0.00213592,0.000840044,0.000668234,0.000238454,0.000301504,0.000168228,8.43442e-05,0.000128126,0.000105654,3.70585e-05,2.7205e-05,2.49273e-05,8.16048e-06,4.77775e-05,4.69319e-05,8.48453e-05,5.11914e-05,8.38762e-06,1.4865e-06,2.75894e-05,5.48744e-06,2.14826e-06,3.68344e-07};
   double pu_rat17[100] = {0.184787,3.79003,3.41318,2.56169,1.62063,1.51684,1.27726,1.26113,0.614194,1.45583,1.49278,1.48366,1.33344,1.16929,1.07775,1.05425,1.08018,1.12791,1.16472,1.18895,1.21256,1.2385,1.26002,1.27057,1.27251,1.27173,1.27087,1.26708,1.27459,1.25156,1.22231,1.16957,1.11026,1.03788,0.968477,0.910741,0.866555,0.835536,0.788241,0.750281,0.758567,0.79348,0.858783,0.959118,1.09448,1.25686,1.41898,1.49494,1.53091,1.46185,1.33675,1.15523,0.950739,0.750143,0.569498,0.410964,0.289908,0.198831,0.137468,0.0967062,0.0692654,0.0509306,0.0384321,0.0299769,0.0240993,0.0170652,0.0124912,0.0107738,0.00962441,0.00883465,0.00829001,0.00803282,0.0078523,0.00788393,0.00629992,0.00534141,0.00546467,0.00547334,0.00592796,0.00592861,0.00624358,0.00643277,0.00651252,0.00487991,0.00427507,0.00466825,0.00400731,0.00487422,0.0048154,0.00462831,0.00371936,0.00384877,0.00172295,0.00205589,0.00618057,0.003259,0,0,0,0};
   double pu_rat17_up[100] = {0.178992,3.27102,2.79307,2.57516,1.39487,1.34319,1.23368,1.15366,0.494409,0.942598,1.06276,1.03109,1.01116,0.904216,0.825826,0.796697,0.811205,0.867756,0.932998,0.985497,1.03408,1.08505,1.12857,1.15501,1.1651,1.16695,1.17204,1.18246,1.20999,1.21206,1.20809,1.17863,1.14084,1.08891,1.03752,0.993071,0.95677,0.929451,0.879146,0.832972,0.828535,0.839762,0.868305,0.920017,1.00058,1.11339,1.24912,1.34572,1.44829,1.48763,1.49045,1.42965,1.31614,1.166,0.993993,0.802787,0.630059,0.477004,0.360771,0.274989,0.211373,0.165255,0.131414,0.107083,0.0891493,0.064786,0.0482185,0.0419048,0.0374053,0.0340798,0.0316022,0.0302138,0.0291769,0.0290428,0.023128,0.01966,0.020291,0.0206234,0.0227871,0.0233585,0.0253166,0.0269401,0.0282591,0.0220026,0.0200819,0.022903,0.0205816,0.0262654,0.0272823,0.0276255,0.0234329,0.0256411,0.0121589,0.0153934,0.0491761,0.027596,0,0,0,0};
   double pu_rat17_dn[100] = {0.193294,4.19484,4.45873,2.47598,1.86911,1.69223,1.33271,1.41359,0.888937,2.16776,2.18378,2.06454,1.74543,1.53951,1.43826,1.42375,1.43283,1.43299,1.42332,1.40779,1.39383,1.38957,1.39141,1.39303,1.39175,1.3814,1.36055,1.32982,1.30823,1.25667,1.20167,1.1248,1.04328,0.954728,0.87675,0.816267,0.773147,0.747541,0.715911,0.703578,0.746748,0.827213,0.943586,1.08963,1.24865,1.39347,1.48264,1.43477,1.32354,1.1237,0.906995,0.690058,0.50081,0.350321,0.237698,0.154823,0.099652,0.0630517,0.0406514,0.0269448,0.0183673,0.0129834,0.00951646,0.00728689,0.00581121,0.00412079,0.00304304,0.00265998,0.00241094,0.00224079,0.002119,0.00205616,0.00199868,0.00198159,0.00155356,0.00128494,0.00127604,0.00123525,0.00128814,0.00123622,0.00124541,0.00122396,0.00117878,0.000838107,0.00069499,0.000716697,0.000579737,0.0006631,0.000614823,0.000553572,0.000415994,0.000401861,0.000167673,0.000186188,0.000520112,0.000254473,0,0,0,0};
   
   double pu_rat16[100] = {0.332136,0.962193,1.20337,0.953956,1.09428,1.19764,0.789015,0.492152,0.747769,0.883546,0.966658,1.07117,1.12602,1.18204,1.20012,1.20987,1.1994,1.18141,1.1441,1.09757,1.06572,1.05093,1.05129,1.05066,1.04841,1.05855,1.07068,1.08196,1.09531,1.11184,1.0946,1.08417,1.04056,0.98214,0.912153,0.821718,0.718354,0.607989,0.500971,0.405178,0.308788,0.228402,0.163573,0.113394,0.0776009,0.051062,0.031858,0.0199973,0.0123033,0.00743405,0.00433781,0.00262109,0.00156357,0.000968549,0.000734598,0.000668329,0.00072343,0.000927116,0.00132969,0.0018984,0.00324026,0.00383716,0.00442283,0.00523648,0.00556154,0.00496787,0.00426864,0.00463808,0.00410515,0.00343628,0.00275464,0.00291527,0.00236412,0.00193081,0.00181666,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   double pu_rat16_up[100] = {0.323632,0.757628,1.13844,0.838277,0.995368,1.0786,0.719714,0.345182,0.504351,0.60627,0.633794,0.731941,0.828284,0.917174,0.958354,0.990304,1.02376,1.05169,1.05114,1.02817,1.0057,0.997764,1.01463,1.03755,1.05665,1.08588,1.11913,1.15394,1.19208,1.23585,1.24624,1.26966,1.25819,1.22921,1.18324,1.10585,1.00459,0.886452,0.765789,0.654404,0.531798,0.423485,0.329565,0.250409,0.189287,0.138533,0.0967281,0.0683142,0.047497,0.0325214,0.0214946,0.0146031,0.0095695,0.00614978,0.0043091,0.00305056,0.00221997,0.00187861,0.00198466,0.0024027,0.00383849,0.00447581,0.00519436,0.00625013,0.00677324,0.00618513,0.005439,0.00605281,0.00549063,0.00471326,0.00387696,0.00421262,0.00350947,0.00294617,0.00285097,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   double pu_rat16_dn[100] = {0.344112,1.22832,1.26582,1.08892,1.2203,1.31689,0.912527,0.762007,1.10079,1.34454,1.48998,1.52638,1.4987,1.50861,1.49484,1.44637,1.36695,1.29754,1.22745,1.16669,1.12555,1.09034,1.06374,1.03998,1.01783,1.00639,0.995607,0.98391,0.972505,0.959983,0.91466,0.873487,0.806578,0.731577,0.651992,0.561942,0.467416,0.373316,0.287326,0.214697,0.149552,0.100106,0.0643017,0.0396644,0.023984,0.0138603,0.00755841,0.00413608,0.00222283,0.00118804,0.00063654,0.000386056,0.000272434,0.000242664,0.000291798,0.000393216,0.000538729,0.000764926,0.00112959,0.0016076,0.00269961,0.00312853,0.00352046,0.00406384,0.00420435,0.00365555,0.00305527,0.00322686,0.00277438,0.0022544,0.00175318,0.00179875,0.0014132,0.00111745,0.00101725,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   
   double pu_rat18[100] =  {0,13.0204,57.5634,19.082,12.1118,8.99654,6.66583,4.91429,3.59483,2.74519,2.23511,1.90123,1.70431,1.57807,1.50263,1.4648,1.45006,1.45373,1.46061,1.46339,1.45963,1.43638,1.4021,1.35249,1.301,1.24853,1.20344,1.16808,1.13747,1.1148,1.09898,1.08944,1.08188,1.077,1.07364,1.069,1.06282,1.05298,1.03554,1.01205,0.979753,0.936484,0.885825,0.827701,0.760395,0.692411,0.620801,0.550926,0.48254,0.419756,0.361541,0.310761,0.266044,0.227951,0.194841,0.167998,0.145255,0.125278,0.109776,0.0961932,0.0837946,0.0738915,0.0643919,0.0565295,0.0493426,0.0422786,0.0363522,0.0314853,0.026336,0.022446,0.0189216,0.0164908,0.0137806,0.0114497,0.00946543,0.00846572,0.0066371,0.00564392,0.00468337,0.00390401,0.00355194,0.00233117,0.00205709,0.0015375,0.00122804,0.00095344,0.000598666,0.000539365,0.000235921,0.000101888,6.25912e-05,6.00391e-05,1.95961e-05,2.71883e-05,1.75359e-05,6.18118e-06,2.54377e-06,1.48936e-06,2.02377e-06,2.56097e-07};
   double pu_rat18_up[100] = {0,11.3701,49.1593,16.3978,10.4484,7.79227,5.70396,4.15872,3.02768,2.28549,1.82582,1.52983,1.3595,1.2554,1.19605,1.1684,1.16115,1.17185,1.18964,1.20936,1.22873,1.23491,1.23159,1.21107,1.18259,1.14644,1.11133,1.08136,1.05384,1.03331,1.01987,1.01367,1.01107,1.01298,1.01865,1.02593,1.03512,1.0447,1.05099,1.0554,1.05447,1.04466,1.02824,1.00332,0.965566,0.923431,0.871249,0.814665,0.752156,0.689408,0.624858,0.564,0.505617,0.452167,0.402,0.359344,0.321227,0.285921,0.258403,0.233682,0.210464,0.192413,0.174424,0.159861,0.146181,0.131623,0.119227,0.10899,0.0963316,0.086803,0.0773651,0.0712667,0.0629173,0.0552031,0.0481823,0.0455058,0.0376989,0.0339163,0.0298286,0.0264131,0.0255965,0.0179475,0.0169746,0.0136435,0.0117583,0.00988318,0.00674005,0.00661599,0.00316237,0.00149674,0.0010104,0.00106782,0.000384941,0.000591271,0.000423128,0.000165822,7.60044e-05,4.96232e-05,7.51979e-05,1.05862e-05};
   double pu_rat18_dn[100] = {0,15.0557,67.8751,22.3278,14.1211,10.4821,7.88069,5.86513,4.31762,3.35551,2.78627,2.40097,2.16428,2.00485,1.9056,1.85092,1.82051,1.80608,1.78719,1.75544,1.71117,1.64481,1.57234,1.49261,1.42092,1.35612,1.3043,1.26517,1.23118,1.20443,1.18302,1.16596,1.14834,1.13047,1.11055,1.08517,1.05388,1.01479,0.96502,0.907499,0.841466,0.767187,0.68971,0.610695,0.530471,0.45611,0.385995,0.32355,0.268127,0.221267,0.181416,0.149012,0.122387,0.100955,0.0832931,0.0694147,0.0579993,0.0482614,0.0406839,0.0341693,0.0284128,0.0238208,0.0196651,0.0163071,0.0134164,0.0108213,0.00875349,0.00713274,0.00561523,0.00450669,0.00357902,0.00293888,0.00231295,0.00180802,0.00140385,0.00117654,0.000861839,0.000682485,0.000525487,0.000404909,0.00033922,0.000204219,0.000164688,0.000112084,8.12391e-05,5.70485e-05,3.2298e-05,2.61592e-05,1.02574e-05,3.96059e-06,2.16985e-06,1.85204e-06,5.36884e-07,6.60936e-07,3.78607e-07,1.19189e-07,4.4536e-08,2.4673e-08,3.47283e-08,5.35281e-09};
   
   double logrhobins[norhobins+1] = {-0.088059,0.0942625,0.276584,0.458906,0.641227,0.823549,1.00587,1.18819,1.37051,1.55283,1.73516,1.91748,2.0998,2.28212,2.46444,2.64676,2.82909,3.01141,3.19373,3.37605,3.55837,3.74069,3.92302,4.10534,4.28766,4.46998,4.6523,4.83462,5.01694,5.19927,5.38159,5.56391,5.74623,5.92855,6.11087,6.2932,6.47552,6.65784,6.84016,7.02248,7.2048,7.38712,7.56945,7.75177,7.93409,8.11641,8.29873,8.48105,8.66338,8.8457,9.02802,9.21034};
    
   float sdmassbins[nsdmassbins+1] = {0,10,21,33,45,57,69,81,93,105,117,129,141,156,172,190,210,250,290,340,400};
   
//    double tau32SF[3][4] = {{0.85161,0.969492,0.834109,0.93509},{1.21039,0.994979,2.30017,1.11669},{1.41352,0.796524,0.960451,0.766512}}; 
   double tau32SF_nom[3][4] = {{0.859549,0.882645,0.863011,0.977778},{0.629492,1.33655,1.62628,0.586884},{1.51485,0.855309,0.747111,1.50843}}; 
   double tau32SF_tight[3][4] = {{0.795027,0.80322,0.741334,0.738192},{0.606943,1.58914,1.01724,1.11726},{1.32298,0.842836,1.95535,1.2135}};
   double tau32SF_loose[3][4] = {{0.91259,0.899487,0.901451,0.991445},{0.805296,1.18465,0.807588,0.97255},{1.49036,1.17164,1.41092,1.37033}};
  							
   double DeepAK8_SF_16[noAK8ptbins] = {1.01,0.88,0.94,1.01};  // 2016  for 0.5% WP
   double DeepAK8_MD_SF_16[noAK8ptbins] = {0.89,1.02,0.93,1.0};  // 2016  for 0.5% WP
   double DeepAK8_MD_SF_16_up[noAK8ptbins] = {0.97,1.07,0.97,1.05};	
   double DeepAK8_MD_SF_16_dn[noAK8ptbins] = {0.81,0.97,0.89,0.95};
   
   double DeepAK8_SF_17[noAK8ptbins] = {0.97,0.97,0.95,1.03};   // 2017  for 0.5% WP
   double DeepAK8_MD_SF_17[noAK8ptbins] = {0.95,1.00,0.98,0.98};  // 2017  for 0.5% WP
   double DeepAK8_MD_SF_17_up[noAK8ptbins] = {1.01,1.04,1.02,1.02};	
   double DeepAK8_MD_SF_17_dn[noAK8ptbins] = {0.89,0.96,0.94,0.94};	   
   
   double DeepAK8_SF_18[noAK8ptbins] = {0.89,1.06,1.00,0.99};  // 2018  for 0.5% WP
   double DeepAK8_MD_SF_18[noAK8ptbins] = {0.90,0.97,0.98,0.95};  // 2018  for 0.5% WP
   double DeepAK8_MD_SF_18_up[noAK8ptbins] = {0.95,1.00,1.01,0.98};	
   double DeepAK8_MD_SF_18_dn[noAK8ptbins] = {0.85,0.94,0.95,0.92};	
   
   double DeepAK8_SF[noAK8ptbins]={0};
   double DeepAK8_MD_SF[noAK8ptbins]={0};
   double DeepAK8_MD_SF_up[noAK8ptbins]={0};
   double DeepAK8_MD_SF_dn[noAK8ptbins]={0};
   
   double cortt_sdmass[nsdmassbins] = {0.704431,0.797327,0.957664,0.718151,0.691868,0.893276,0.510895,0.725355,0.789391,0.847687,0.911654,0.892455,0.892299,0.81116,0.713741,1.12043,0.9146,0.806993,1.34569,0.943199};

   double cor_b_sdmass_mc_18[2][nsdmassbins] = {
    {0.923822,1.04154,1.0317,1.12539,1.28814,1.12933,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{1.00665,0.950334,0.92359,1.02488,1.12131,1.22944,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
   };  //2018 QCD
   
   double cor_b_sdmass_data_18[2][nsdmassbins] = {
    {0.933685,0.955174,1.02156,1.19531,1.33257,1.38726,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0.970126,0.900778,0.991004,1.1816,1.37485,1.41539,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
   };  //2018 Data
   
 
   double cor_b_sdmass_mc_16[2][nsdmassbins] = {
    {0.981351,1.09841,0.932487,0.985458,0.982843,1.00175,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{1.10254,1.01459,0.817372,0.751854,0.77791,0.856823,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
   };  //2016 QCD
   
   double cor_b_sdmass_data_16[2][nsdmassbins] = {
    {1.00394,1.0785,0.963856,0.899902,0.906215,0.871693,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{1.06831,1.0291,0.876079,0.786874,0.772602,0.74453,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
   };  //2016 Data
    
     
   double cor_b_sdmass_mc_17[2][nsdmassbins] = {
    {0.886976,0.968629,1.16829,1.2444,1.41973,1.51195,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0.941006,0.974473,1.03761,1.13087,1.34041,1.36105,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
   };  //2017 QCD
   
   double cor_b_sdmass_data_17[2][nsdmassbins] = {
    {0.929639,0.981946,1.03214,1.18439,1.32073,1.29275,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0.96212,0.915459,1.01051,1.18366,1.35503,1.43479,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
   };  //2017 Data
 
   double cor_b_sdmass_mc[2][nsdmassbins] = {{0}};
   double cor_b_sdmass_data[2][nsdmassbins] = {{0}};
   
   static const int nohtbins = 30;
     
   double htbins[nohtbins+1] = {300, 330, 362, 395, 430, 468,
     507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
     1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1784, 2000,
     2116, 2500, 2941, 3637, 5000} ; 
   
   float sfwt_tau32;
   float sfwt_deepak8;
   float  sfwt_deepak8_md, sfwt_deepak8_md_up, sfwt_deepak8_md_dn; 

   Anal_Nano_PROOF(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~Anal_Nano_PROOF() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(Anal_Nano_PROOF,0);
};

#endif

#ifdef Anal_Nano_PROOF_cxx
void Anal_Nano_PROOF::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("btagWeight_CSVV2", &btagWeight_CSVV2, &b_btagWeight_CSVV2);
   fChain->SetBranchAddress("btagWeight_DeepCSVB", &btagWeight_DeepCSVB, &b_btagWeight_DeepCSVB);
   fChain->SetBranchAddress("CaloMET_phi", &CaloMET_phi, &b_CaloMET_phi);
   fChain->SetBranchAddress("CaloMET_pt", &CaloMET_pt, &b_CaloMET_pt);
   fChain->SetBranchAddress("CaloMET_sumEt", &CaloMET_sumEt, &b_CaloMET_sumEt);
   fChain->SetBranchAddress("ChsMET_phi", &ChsMET_phi, &b_ChsMET_phi);
   fChain->SetBranchAddress("ChsMET_pt", &ChsMET_pt, &b_ChsMET_pt);
   fChain->SetBranchAddress("ChsMET_sumEt", &ChsMET_sumEt, &b_ChsMET_sumEt);
   fChain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
   fChain->SetBranchAddress("Electron_deltaEtaSC", Electron_deltaEtaSC, &b_Electron_deltaEtaSC);
   fChain->SetBranchAddress("Electron_dr03EcalRecHitSumEt", Electron_dr03EcalRecHitSumEt, &b_Electron_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEt", Electron_dr03HcalDepth1TowerSumEt, &b_Electron_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("Electron_dr03TkSumPt", Electron_dr03TkSumPt, &b_Electron_dr03TkSumPt);
   fChain->SetBranchAddress("Electron_dr03TkSumPtHEEP", Electron_dr03TkSumPtHEEP, &b_Electron_dr03TkSumPtHEEP);
   fChain->SetBranchAddress("Electron_dxy", Electron_dxy, &b_Electron_dxy);
   fChain->SetBranchAddress("Electron_dxyErr", Electron_dxyErr, &b_Electron_dxyErr);
   fChain->SetBranchAddress("Electron_dz", Electron_dz, &b_Electron_dz);
   fChain->SetBranchAddress("Electron_dzErr", Electron_dzErr, &b_Electron_dzErr);
   fChain->SetBranchAddress("Electron_eCorr", Electron_eCorr, &b_Electron_eCorr);
   fChain->SetBranchAddress("Electron_eInvMinusPInv", Electron_eInvMinusPInv, &b_Electron_eInvMinusPInv);
   fChain->SetBranchAddress("Electron_energyErr", Electron_energyErr, &b_Electron_energyErr);
   fChain->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron_hoe", Electron_hoe, &b_Electron_hoe);
   fChain->SetBranchAddress("Electron_ip3d", Electron_ip3d, &b_Electron_ip3d);
   fChain->SetBranchAddress("Electron_jetRelIso", Electron_jetRelIso, &b_Electron_jetRelIso);
   fChain->SetBranchAddress("Electron_mass", Electron_mass, &b_Electron_mass);
   fChain->SetBranchAddress("Electron_miniPFRelIso_all", Electron_miniPFRelIso_all, &b_Electron_miniPFRelIso_all);
   fChain->SetBranchAddress("Electron_miniPFRelIso_chg", Electron_miniPFRelIso_chg, &b_Electron_miniPFRelIso_chg);
   fChain->SetBranchAddress("Electron_mvaFall17V1Iso", Electron_mvaFall17V1Iso, &b_Electron_mvaFall17V1Iso);
   fChain->SetBranchAddress("Electron_mvaFall17V1noIso", Electron_mvaFall17V1noIso, &b_Electron_mvaFall17V1noIso);
   fChain->SetBranchAddress("Electron_mvaFall17V2Iso", Electron_mvaFall17V2Iso, &b_Electron_mvaFall17V2Iso);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso", Electron_mvaFall17V2noIso, &b_Electron_mvaFall17V2noIso);
   fChain->SetBranchAddress("Electron_pfRelIso03_all", Electron_pfRelIso03_all, &b_Electron_pfRelIso03_all);
   fChain->SetBranchAddress("Electron_pfRelIso03_chg", Electron_pfRelIso03_chg, &b_Electron_pfRelIso03_chg);
   fChain->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_r9", Electron_r9, &b_Electron_r9);
   fChain->SetBranchAddress("Electron_sieie", Electron_sieie, &b_Electron_sieie);
   fChain->SetBranchAddress("Electron_sip3d", Electron_sip3d, &b_Electron_sip3d);
   fChain->SetBranchAddress("Electron_mvaTTH", Electron_mvaTTH, &b_Electron_mvaTTH);
   fChain->SetBranchAddress("Electron_charge", Electron_charge, &b_Electron_charge);
   fChain->SetBranchAddress("Electron_cutBased", Electron_cutBased, &b_Electron_cutBased);
   fChain->SetBranchAddress("Electron_cutBased_Fall17_V1", Electron_cutBased_Fall17_V1, &b_Electron_cutBased_Fall17_V1);
   fChain->SetBranchAddress("Electron_jetIdx", Electron_jetIdx, &b_Electron_jetIdx);
   fChain->SetBranchAddress("Electron_pdgId", Electron_pdgId, &b_Electron_pdgId);
   fChain->SetBranchAddress("Electron_photonIdx", Electron_photonIdx, &b_Electron_photonIdx);
   fChain->SetBranchAddress("Electron_tightCharge", Electron_tightCharge, &b_Electron_tightCharge);
   fChain->SetBranchAddress("Electron_vidNestedWPBitmap", Electron_vidNestedWPBitmap, &b_Electron_vidNestedWPBitmap);
   fChain->SetBranchAddress("Electron_convVeto", Electron_convVeto, &b_Electron_convVeto);
   fChain->SetBranchAddress("Electron_cutBased_HEEP", Electron_cutBased_HEEP, &b_Electron_cutBased_HEEP);
   fChain->SetBranchAddress("Electron_isPFcand", Electron_isPFcand, &b_Electron_isPFcand);
   fChain->SetBranchAddress("Electron_lostHits", Electron_lostHits, &b_Electron_lostHits);
   fChain->SetBranchAddress("Electron_mvaFall17V1Iso_WP80", Electron_mvaFall17V1Iso_WP80, &b_Electron_mvaFall17V1Iso_WP80);
   fChain->SetBranchAddress("Electron_mvaFall17V1Iso_WP90", Electron_mvaFall17V1Iso_WP90, &b_Electron_mvaFall17V1Iso_WP90);
   fChain->SetBranchAddress("Electron_mvaFall17V1Iso_WPL", Electron_mvaFall17V1Iso_WPL, &b_Electron_mvaFall17V1Iso_WPL);
   fChain->SetBranchAddress("Electron_mvaFall17V1noIso_WP80", Electron_mvaFall17V1noIso_WP80, &b_Electron_mvaFall17V1noIso_WP80);
   fChain->SetBranchAddress("Electron_mvaFall17V1noIso_WP90", Electron_mvaFall17V1noIso_WP90, &b_Electron_mvaFall17V1noIso_WP90);
   fChain->SetBranchAddress("Electron_mvaFall17V1noIso_WPL", Electron_mvaFall17V1noIso_WPL, &b_Electron_mvaFall17V1noIso_WPL);
   fChain->SetBranchAddress("Electron_mvaFall17V2Iso_WP80", Electron_mvaFall17V2Iso_WP80, &b_Electron_mvaFall17V2Iso_WP80);
   fChain->SetBranchAddress("Electron_mvaFall17V2Iso_WP90", Electron_mvaFall17V2Iso_WP90, &b_Electron_mvaFall17V2Iso_WP90);
   fChain->SetBranchAddress("Electron_mvaFall17V2Iso_WPL", Electron_mvaFall17V2Iso_WPL, &b_Electron_mvaFall17V2Iso_WPL);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso_WP80", Electron_mvaFall17V2noIso_WP80, &b_Electron_mvaFall17V2noIso_WP80);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso_WP90", Electron_mvaFall17V2noIso_WP90, &b_Electron_mvaFall17V2noIso_WP90);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso_WPL", Electron_mvaFall17V2noIso_WPL, &b_Electron_mvaFall17V2noIso_WPL);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilterV2", &Flag_ecalBadCalibFilterV2, &b_Flag_ecalBadCalibFilterV2);
   fChain->SetBranchAddress("nFatJet", &nFatJet, &b_nFatJet);
   fChain->SetBranchAddress("FatJet_area", FatJet_area, &b_FatJet_area);
   fChain->SetBranchAddress("FatJet_btagCMVA", FatJet_btagCMVA, &b_FatJet_btagCMVA);
   fChain->SetBranchAddress("FatJet_btagCSVV2", FatJet_btagCSVV2, &b_FatJet_btagCSVV2);
   fChain->SetBranchAddress("FatJet_btagDeepB", FatJet_btagDeepB, &b_FatJet_btagDeepB);
   fChain->SetBranchAddress("FatJet_btagHbb", FatJet_btagHbb, &b_FatJet_btagHbb);
   fChain->SetBranchAddress("FatJet_deepTagMD_H4qvsQCD", FatJet_deepTagMD_H4qvsQCD, &b_FatJet_deepTagMD_H4qvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_HbbvsQCD", FatJet_deepTagMD_HbbvsQCD, &b_FatJet_deepTagMD_HbbvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_TvsQCD", FatJet_deepTagMD_TvsQCD, &b_FatJet_deepTagMD_TvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_WvsQCD", FatJet_deepTagMD_WvsQCD, &b_FatJet_deepTagMD_WvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_ZHbbvsQCD", FatJet_deepTagMD_ZHbbvsQCD, &b_FatJet_deepTagMD_ZHbbvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_ZHccvsQCD", FatJet_deepTagMD_ZHccvsQCD, &b_FatJet_deepTagMD_ZHccvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_ZbbvsQCD", FatJet_deepTagMD_ZbbvsQCD, &b_FatJet_deepTagMD_ZbbvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_ZvsQCD", FatJet_deepTagMD_ZvsQCD, &b_FatJet_deepTagMD_ZvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_bbvsLight", FatJet_deepTagMD_bbvsLight, &b_FatJet_deepTagMD_bbvsLight);
   fChain->SetBranchAddress("FatJet_deepTagMD_ccvsLight", FatJet_deepTagMD_ccvsLight, &b_FatJet_deepTagMD_ccvsLight);
   fChain->SetBranchAddress("FatJet_deepTag_H", FatJet_deepTag_H, &b_FatJet_deepTag_H);
   fChain->SetBranchAddress("FatJet_deepTag_QCD", FatJet_deepTag_QCD, &b_FatJet_deepTag_QCD);
   fChain->SetBranchAddress("FatJet_deepTag_QCDothers", FatJet_deepTag_QCDothers, &b_FatJet_deepTag_QCDothers);
   fChain->SetBranchAddress("FatJet_deepTag_TvsQCD", FatJet_deepTag_TvsQCD, &b_FatJet_deepTag_TvsQCD);
   fChain->SetBranchAddress("FatJet_deepTag_WvsQCD", FatJet_deepTag_WvsQCD, &b_FatJet_deepTag_WvsQCD);
   fChain->SetBranchAddress("FatJet_deepTag_ZvsQCD", FatJet_deepTag_ZvsQCD, &b_FatJet_deepTag_ZvsQCD);
   fChain->SetBranchAddress("FatJet_eta", FatJet_eta, &b_FatJet_eta);
   fChain->SetBranchAddress("FatJet_mass", FatJet_mass, &b_FatJet_mass);
   fChain->SetBranchAddress("FatJet_msoftdrop", FatJet_msoftdrop, &b_FatJet_msoftdrop);
   fChain->SetBranchAddress("FatJet_n2b1", FatJet_n2b1, &b_FatJet_n2b1);
   fChain->SetBranchAddress("FatJet_n3b1", FatJet_n3b1, &b_FatJet_n3b1);
   fChain->SetBranchAddress("FatJet_phi", FatJet_phi, &b_FatJet_phi);
   fChain->SetBranchAddress("FatJet_pt", FatJet_pt, &b_FatJet_pt);
   fChain->SetBranchAddress("FatJet_rawFactor", FatJet_rawFactor, &b_FatJet_rawFactor);
   fChain->SetBranchAddress("FatJet_tau1", FatJet_tau1, &b_FatJet_tau1);
   fChain->SetBranchAddress("FatJet_tau2", FatJet_tau2, &b_FatJet_tau2);
   fChain->SetBranchAddress("FatJet_tau3", FatJet_tau3, &b_FatJet_tau3);
   fChain->SetBranchAddress("FatJet_tau4", FatJet_tau4, &b_FatJet_tau4);
   fChain->SetBranchAddress("FatJet_jetId", FatJet_jetId, &b_FatJet_jetId);
   fChain->SetBranchAddress("FatJet_subJetIdx1", FatJet_subJetIdx1, &b_FatJet_subJetIdx1);
   fChain->SetBranchAddress("FatJet_subJetIdx2", FatJet_subJetIdx2, &b_FatJet_subJetIdx2);
   fChain->SetBranchAddress("nGenJetAK8", &nGenJetAK8, &b_nGenJetAK8);
   fChain->SetBranchAddress("GenJetAK8_eta", GenJetAK8_eta, &b_GenJetAK8_eta);
   fChain->SetBranchAddress("GenJetAK8_mass", GenJetAK8_mass, &b_GenJetAK8_mass);
   fChain->SetBranchAddress("GenJetAK8_phi", GenJetAK8_phi, &b_GenJetAK8_phi);
   fChain->SetBranchAddress("GenJetAK8_pt", GenJetAK8_pt, &b_GenJetAK8_pt);
   fChain->SetBranchAddress("nGenJet", &nGenJet, &b_nGenJet);
   fChain->SetBranchAddress("GenJet_eta", GenJet_eta, &b_GenJet_eta);
   fChain->SetBranchAddress("GenJet_mass", GenJet_mass, &b_GenJet_mass);
   fChain->SetBranchAddress("GenJet_phi", GenJet_phi, &b_GenJet_phi);
   fChain->SetBranchAddress("GenJet_pt", GenJet_pt, &b_GenJet_pt);
   fChain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
   fChain->SetBranchAddress("GenPart_eta", GenPart_eta, &b_GenPart_eta);
   fChain->SetBranchAddress("GenPart_mass", GenPart_mass, &b_GenPart_mass);
   fChain->SetBranchAddress("GenPart_phi", GenPart_phi, &b_GenPart_phi);
   fChain->SetBranchAddress("GenPart_pt", GenPart_pt, &b_GenPart_pt);
   fChain->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother, &b_GenPart_genPartIdxMother);
   fChain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId, &b_GenPart_pdgId);
   fChain->SetBranchAddress("GenPart_status", GenPart_status, &b_GenPart_status);
   fChain->SetBranchAddress("GenPart_statusFlags", GenPart_statusFlags, &b_GenPart_statusFlags);
   fChain->SetBranchAddress("nSubGenJetAK8", &nSubGenJetAK8, &b_nSubGenJetAK8);
   fChain->SetBranchAddress("SubGenJetAK8_eta", SubGenJetAK8_eta, &b_SubGenJetAK8_eta);
   fChain->SetBranchAddress("SubGenJetAK8_mass", SubGenJetAK8_mass, &b_SubGenJetAK8_mass);
   fChain->SetBranchAddress("SubGenJetAK8_phi", SubGenJetAK8_phi, &b_SubGenJetAK8_phi);
   fChain->SetBranchAddress("SubGenJetAK8_pt", SubGenJetAK8_pt, &b_SubGenJetAK8_pt);
   fChain->SetBranchAddress("Generator_binvar", &Generator_binvar, &b_Generator_binvar);
   fChain->SetBranchAddress("Generator_scalePDF", &Generator_scalePDF, &b_Generator_scalePDF);
   fChain->SetBranchAddress("Generator_weight", &Generator_weight, &b_Generator_weight);
   fChain->SetBranchAddress("Generator_x1", &Generator_x1, &b_Generator_x1);
   fChain->SetBranchAddress("Generator_x2", &Generator_x2, &b_Generator_x2);
   fChain->SetBranchAddress("Generator_xpdf1", &Generator_xpdf1, &b_Generator_xpdf1);
   fChain->SetBranchAddress("Generator_xpdf2", &Generator_xpdf2, &b_Generator_xpdf2);
   fChain->SetBranchAddress("Generator_id1", &Generator_id1, &b_Generator_id1);
   fChain->SetBranchAddress("Generator_id2", &Generator_id2, &b_Generator_id2);
   fChain->SetBranchAddress("nGenVisTau", &nGenVisTau, &b_nGenVisTau);
   fChain->SetBranchAddress("GenVisTau_eta", GenVisTau_eta, &b_GenVisTau_eta);
   fChain->SetBranchAddress("GenVisTau_mass", GenVisTau_mass, &b_GenVisTau_mass);
   fChain->SetBranchAddress("GenVisTau_phi", GenVisTau_phi, &b_GenVisTau_phi);
   fChain->SetBranchAddress("GenVisTau_pt", GenVisTau_pt, &b_GenVisTau_pt);
   fChain->SetBranchAddress("GenVisTau_charge", GenVisTau_charge, &b_GenVisTau_charge);
   fChain->SetBranchAddress("GenVisTau_genPartIdxMother", GenVisTau_genPartIdxMother, &b_GenVisTau_genPartIdxMother);
   fChain->SetBranchAddress("GenVisTau_status", GenVisTau_status, &b_GenVisTau_status);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("LHEWeight_originalXWGTUP", &LHEWeight_originalXWGTUP, &b_LHEWeight_originalXWGTUP);
   fChain->SetBranchAddress("nLHEPdfWeight", &nLHEPdfWeight, &b_nLHEPdfWeight);
   fChain->SetBranchAddress("LHEPdfWeight", LHEPdfWeight, &b_LHEPdfWeight);
   fChain->SetBranchAddress("nLHEScaleWeight", &nLHEScaleWeight, &b_nLHEScaleWeight);
   fChain->SetBranchAddress("LHEScaleWeight", LHEScaleWeight, &b_LHEScaleWeight);
   fChain->SetBranchAddress("nPSWeight", &nPSWeight, &b_nPSWeight);
   fChain->SetBranchAddress("PSWeight", PSWeight, &b_PSWeight);
   fChain->SetBranchAddress("nIsoTrack", &nIsoTrack, &b_nIsoTrack);
   fChain->SetBranchAddress("IsoTrack_dxy", IsoTrack_dxy, &b_IsoTrack_dxy);
   fChain->SetBranchAddress("IsoTrack_dz", IsoTrack_dz, &b_IsoTrack_dz);
   fChain->SetBranchAddress("IsoTrack_eta", IsoTrack_eta, &b_IsoTrack_eta);
   fChain->SetBranchAddress("IsoTrack_pfRelIso03_all", IsoTrack_pfRelIso03_all, &b_IsoTrack_pfRelIso03_all);
   fChain->SetBranchAddress("IsoTrack_pfRelIso03_chg", IsoTrack_pfRelIso03_chg, &b_IsoTrack_pfRelIso03_chg);
   fChain->SetBranchAddress("IsoTrack_phi", IsoTrack_phi, &b_IsoTrack_phi);
   fChain->SetBranchAddress("IsoTrack_pt", IsoTrack_pt, &b_IsoTrack_pt);
   fChain->SetBranchAddress("IsoTrack_miniPFRelIso_all", IsoTrack_miniPFRelIso_all, &b_IsoTrack_miniPFRelIso_all);
   fChain->SetBranchAddress("IsoTrack_miniPFRelIso_chg", IsoTrack_miniPFRelIso_chg, &b_IsoTrack_miniPFRelIso_chg);
   fChain->SetBranchAddress("IsoTrack_fromPV", IsoTrack_fromPV, &b_IsoTrack_fromPV);
   fChain->SetBranchAddress("IsoTrack_pdgId", IsoTrack_pdgId, &b_IsoTrack_pdgId);
   fChain->SetBranchAddress("IsoTrack_isHighPurityTrack", IsoTrack_isHighPurityTrack, &b_IsoTrack_isHighPurityTrack);
   fChain->SetBranchAddress("IsoTrack_isPFcand", IsoTrack_isPFcand, &b_IsoTrack_isPFcand);
   fChain->SetBranchAddress("IsoTrack_isFromLostTrack", IsoTrack_isFromLostTrack, &b_IsoTrack_isFromLostTrack);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_area", Jet_area, &b_Jet_area);
   fChain->SetBranchAddress("Jet_btagCMVA", Jet_btagCMVA, &b_Jet_btagCMVA);
   fChain->SetBranchAddress("Jet_btagCSVV2", Jet_btagCSVV2, &b_Jet_btagCSVV2);
   fChain->SetBranchAddress("Jet_btagDeepB", Jet_btagDeepB, &b_Jet_btagDeepB);
   fChain->SetBranchAddress("Jet_btagDeepC", Jet_btagDeepC, &b_Jet_btagDeepC);
   fChain->SetBranchAddress("Jet_btagDeepFlavB", Jet_btagDeepFlavB, &b_Jet_btagDeepFlavB);
   fChain->SetBranchAddress("Jet_chEmEF", Jet_chEmEF, &b_Jet_chEmEF);
   fChain->SetBranchAddress("Jet_chHEF", Jet_chHEF, &b_Jet_chHEF);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_muEF", Jet_muEF, &b_Jet_muEF);
   fChain->SetBranchAddress("Jet_neEmEF", Jet_neEmEF, &b_Jet_neEmEF);
   fChain->SetBranchAddress("Jet_neHEF", Jet_neHEF, &b_Jet_neHEF);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_qgl", Jet_qgl, &b_Jet_qgl);
   fChain->SetBranchAddress("Jet_rawFactor", Jet_rawFactor, &b_Jet_rawFactor);
   fChain->SetBranchAddress("Jet_bRegCorr", Jet_bRegCorr, &b_Jet_bRegCorr);
   fChain->SetBranchAddress("Jet_bRegRes", Jet_bRegRes, &b_Jet_bRegRes);
   fChain->SetBranchAddress("Jet_electronIdx1", Jet_electronIdx1, &b_Jet_electronIdx1);
   fChain->SetBranchAddress("Jet_electronIdx2", Jet_electronIdx2, &b_Jet_electronIdx2);
   fChain->SetBranchAddress("Jet_jetId", Jet_jetId, &b_Jet_jetId);
   fChain->SetBranchAddress("Jet_muonIdx1", Jet_muonIdx1, &b_Jet_muonIdx1);
   fChain->SetBranchAddress("Jet_muonIdx2", Jet_muonIdx2, &b_Jet_muonIdx2);
   fChain->SetBranchAddress("Jet_nConstituents", Jet_nConstituents, &b_Jet_nConstituents);
   fChain->SetBranchAddress("Jet_nElectrons", Jet_nElectrons, &b_Jet_nElectrons);
   fChain->SetBranchAddress("Jet_nMuons", Jet_nMuons, &b_Jet_nMuons);
   fChain->SetBranchAddress("Jet_puId", Jet_puId, &b_Jet_puId);
   fChain->SetBranchAddress("LHE_HT", &LHE_HT, &b_LHE_HT);
   fChain->SetBranchAddress("LHE_HTIncoming", &LHE_HTIncoming, &b_LHE_HTIncoming);
   fChain->SetBranchAddress("LHE_Vpt", &LHE_Vpt, &b_LHE_Vpt);
   fChain->SetBranchAddress("LHE_Njets", &LHE_Njets, &b_LHE_Njets);
   fChain->SetBranchAddress("LHE_Nb", &LHE_Nb, &b_LHE_Nb);
   fChain->SetBranchAddress("LHE_Nc", &LHE_Nc, &b_LHE_Nc);
   fChain->SetBranchAddress("LHE_Nuds", &LHE_Nuds, &b_LHE_Nuds);
   fChain->SetBranchAddress("LHE_Nglu", &LHE_Nglu, &b_LHE_Nglu);
   fChain->SetBranchAddress("LHE_NpNLO", &LHE_NpNLO, &b_LHE_NpNLO);
   fChain->SetBranchAddress("LHE_NpLO", &LHE_NpLO, &b_LHE_NpLO);
   fChain->SetBranchAddress("nLHEPart", &nLHEPart, &b_nLHEPart);
   fChain->SetBranchAddress("LHEPart_pt", LHEPart_pt, &b_LHEPart_pt);
   fChain->SetBranchAddress("LHEPart_eta", LHEPart_eta, &b_LHEPart_eta);
   fChain->SetBranchAddress("LHEPart_phi", LHEPart_phi, &b_LHEPart_phi);
   fChain->SetBranchAddress("LHEPart_mass", LHEPart_mass, &b_LHEPart_mass);
   fChain->SetBranchAddress("LHEPart_pdgId", LHEPart_pdgId, &b_LHEPart_pdgId);
   /*
   fChain->SetBranchAddress("METFixEE2017_MetUnclustEnUpDeltaX", &METFixEE2017_MetUnclustEnUpDeltaX, &b_METFixEE2017_MetUnclustEnUpDeltaX);
   fChain->SetBranchAddress("METFixEE2017_MetUnclustEnUpDeltaY", &METFixEE2017_MetUnclustEnUpDeltaY, &b_METFixEE2017_MetUnclustEnUpDeltaY);
   fChain->SetBranchAddress("METFixEE2017_phi", &METFixEE2017_phi, &b_METFixEE2017_phi);
   fChain->SetBranchAddress("METFixEE2017_pt", &METFixEE2017_pt, &b_METFixEE2017_pt);
   fChain->SetBranchAddress("METFixEE2017_sumEt", &METFixEE2017_sumEt, &b_METFixEE2017_sumEt);
   */ 
   
   fChain->SetBranchAddress("METFixEE2017_MetUnclustEnUpDeltaX", &METFixEE2017_MetUnclustEnUpDeltaX, &b_METFixEE2017_MetUnclustEnUpDeltaX);
   fChain->SetBranchAddress("METFixEE2017_MetUnclustEnUpDeltaY", &METFixEE2017_MetUnclustEnUpDeltaY, &b_METFixEE2017_MetUnclustEnUpDeltaY);
   fChain->SetBranchAddress("METFixEE2017_covXX", &METFixEE2017_covXX, &b_METFixEE2017_covXX);
   fChain->SetBranchAddress("METFixEE2017_covXY", &METFixEE2017_covXY, &b_METFixEE2017_covXY);
   fChain->SetBranchAddress("METFixEE2017_covYY", &METFixEE2017_covYY, &b_METFixEE2017_covYY);
   fChain->SetBranchAddress("METFixEE2017_phi", &METFixEE2017_phi, &b_METFixEE2017_phi);
   fChain->SetBranchAddress("METFixEE2017_pt", &METFixEE2017_pt, &b_METFixEE2017_pt);
   fChain->SetBranchAddress("METFixEE2017_significance", &METFixEE2017_significance, &b_METFixEE2017_significance);
   fChain->SetBranchAddress("METFixEE2017_sumEt", &METFixEE2017_sumEt, &b_METFixEE2017_sumEt);
   
   fChain->SetBranchAddress("GenMET_phi", &GenMET_phi, &b_GenMET_phi);
   fChain->SetBranchAddress("GenMET_pt", &GenMET_pt, &b_GenMET_pt);
   fChain->SetBranchAddress("MET_MetUnclustEnUpDeltaX", &MET_MetUnclustEnUpDeltaX, &b_MET_MetUnclustEnUpDeltaX);
   fChain->SetBranchAddress("MET_MetUnclustEnUpDeltaY", &MET_MetUnclustEnUpDeltaY, &b_MET_MetUnclustEnUpDeltaY);
   fChain->SetBranchAddress("MET_covXX", &MET_covXX, &b_MET_covXX);
   fChain->SetBranchAddress("MET_covXY", &MET_covXY, &b_MET_covXY);
   fChain->SetBranchAddress("MET_covYY", &MET_covYY, &b_MET_covYY);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("MET_pt", &MET_pt, &b_MET_pt);
   fChain->SetBranchAddress("MET_significance", &MET_significance, &b_MET_significance);
   fChain->SetBranchAddress("MET_sumEt", &MET_sumEt, &b_MET_sumEt);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("Muon_dxy", Muon_dxy, &b_Muon_dxy);
   fChain->SetBranchAddress("Muon_dxyErr", Muon_dxyErr, &b_Muon_dxyErr);
   fChain->SetBranchAddress("Muon_dz", Muon_dz, &b_Muon_dz);
   fChain->SetBranchAddress("Muon_dzErr", Muon_dzErr, &b_Muon_dzErr);
   fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_ip3d", Muon_ip3d, &b_Muon_ip3d);
   fChain->SetBranchAddress("Muon_jetPtRelv2", Muon_jetPtRelv2, &b_Muon_jetPtRelv2);
   fChain->SetBranchAddress("Muon_jetRelIso", Muon_jetRelIso, &b_Muon_jetRelIso);
   fChain->SetBranchAddress("Muon_mass", Muon_mass, &b_Muon_mass);
   fChain->SetBranchAddress("Muon_miniPFRelIso_all", Muon_miniPFRelIso_all, &b_Muon_miniPFRelIso_all);
   fChain->SetBranchAddress("Muon_miniPFRelIso_chg", Muon_miniPFRelIso_chg, &b_Muon_miniPFRelIso_chg);
   fChain->SetBranchAddress("Muon_pfRelIso03_all", Muon_pfRelIso03_all, &b_Muon_pfRelIso03_all);
   fChain->SetBranchAddress("Muon_pfRelIso03_chg", Muon_pfRelIso03_chg, &b_Muon_pfRelIso03_chg);
   fChain->SetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all, &b_Muon_pfRelIso04_all);
   fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_ptErr", Muon_ptErr, &b_Muon_ptErr);
   fChain->SetBranchAddress("Muon_segmentComp", Muon_segmentComp, &b_Muon_segmentComp);
   fChain->SetBranchAddress("Muon_sip3d", Muon_sip3d, &b_Muon_sip3d);
   fChain->SetBranchAddress("Muon_mvaTTH", Muon_mvaTTH, &b_Muon_mvaTTH);
   fChain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_jetIdx", Muon_jetIdx, &b_Muon_jetIdx);
   fChain->SetBranchAddress("Muon_nStations", Muon_nStations, &b_Muon_nStations);
   fChain->SetBranchAddress("Muon_nTrackerLayers", Muon_nTrackerLayers, &b_Muon_nTrackerLayers);
   fChain->SetBranchAddress("Muon_pdgId", Muon_pdgId, &b_Muon_pdgId);
   fChain->SetBranchAddress("Muon_tightCharge", Muon_tightCharge, &b_Muon_tightCharge);
   fChain->SetBranchAddress("Muon_highPtId", Muon_highPtId, &b_Muon_highPtId);
   fChain->SetBranchAddress("Muon_inTimeMuon", Muon_inTimeMuon, &b_Muon_inTimeMuon);
   fChain->SetBranchAddress("Muon_isGlobal", Muon_isGlobal, &b_Muon_isGlobal);
   fChain->SetBranchAddress("Muon_isPFcand", Muon_isPFcand, &b_Muon_isPFcand);
   fChain->SetBranchAddress("Muon_isTracker", Muon_isTracker, &b_Muon_isTracker);
   fChain->SetBranchAddress("Muon_mediumId", Muon_mediumId, &b_Muon_mediumId);
   fChain->SetBranchAddress("Muon_mediumPromptId", Muon_mediumPromptId, &b_Muon_mediumPromptId);
   fChain->SetBranchAddress("Muon_miniIsoId", Muon_miniIsoId, &b_Muon_miniIsoId);
   fChain->SetBranchAddress("Muon_multiIsoId", Muon_multiIsoId, &b_Muon_multiIsoId);
   fChain->SetBranchAddress("Muon_mvaId", Muon_mvaId, &b_Muon_mvaId);
   fChain->SetBranchAddress("Muon_pfIsoId", Muon_pfIsoId, &b_Muon_pfIsoId);
   fChain->SetBranchAddress("Muon_softId", Muon_softId, &b_Muon_softId);
   fChain->SetBranchAddress("Muon_softMvaId", Muon_softMvaId, &b_Muon_softMvaId);
   fChain->SetBranchAddress("Muon_tightId", Muon_tightId, &b_Muon_tightId);
   fChain->SetBranchAddress("Muon_tkIsoId", Muon_tkIsoId, &b_Muon_tkIsoId);
   fChain->SetBranchAddress("Muon_triggerIdLoose", Muon_triggerIdLoose, &b_Muon_triggerIdLoose);
   fChain->SetBranchAddress("nPhoton", &nPhoton, &b_nPhoton);
   fChain->SetBranchAddress("Photon_eCorr", Photon_eCorr, &b_Photon_eCorr);
   fChain->SetBranchAddress("Photon_energyErr", Photon_energyErr, &b_Photon_energyErr);
   fChain->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_hoe", Photon_hoe, &b_Photon_hoe);
   fChain->SetBranchAddress("Photon_mass", Photon_mass, &b_Photon_mass);
   fChain->SetBranchAddress("Photon_mvaID", Photon_mvaID, &b_Photon_mvaID);
   fChain->SetBranchAddress("Photon_mvaIDV1", Photon_mvaIDV1, &b_Photon_mvaIDV1);
   fChain->SetBranchAddress("Photon_pfRelIso03_all", Photon_pfRelIso03_all, &b_Photon_pfRelIso03_all);
   fChain->SetBranchAddress("Photon_pfRelIso03_chg", Photon_pfRelIso03_chg, &b_Photon_pfRelIso03_chg);
   fChain->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon_pt", Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon_r9", Photon_r9, &b_Photon_r9);
   fChain->SetBranchAddress("Photon_sieie", Photon_sieie, &b_Photon_sieie);
   fChain->SetBranchAddress("Photon_charge", Photon_charge, &b_Photon_charge);
   fChain->SetBranchAddress("Photon_cutBasedBitmap", Photon_cutBasedBitmap, &b_Photon_cutBasedBitmap);
   fChain->SetBranchAddress("Photon_cutBasedV1Bitmap", Photon_cutBasedV1Bitmap, &b_Photon_cutBasedV1Bitmap);
   fChain->SetBranchAddress("Photon_electronIdx", Photon_electronIdx, &b_Photon_electronIdx);
   fChain->SetBranchAddress("Photon_jetIdx", Photon_jetIdx, &b_Photon_jetIdx);
   fChain->SetBranchAddress("Photon_pdgId", Photon_pdgId, &b_Photon_pdgId);
   fChain->SetBranchAddress("Photon_vidNestedWPBitmap", Photon_vidNestedWPBitmap, &b_Photon_vidNestedWPBitmap);
   fChain->SetBranchAddress("Photon_electronVeto", Photon_electronVeto, &b_Photon_electronVeto);
   fChain->SetBranchAddress("Photon_isScEtaEB", Photon_isScEtaEB, &b_Photon_isScEtaEB);
   fChain->SetBranchAddress("Photon_isScEtaEE", Photon_isScEtaEE, &b_Photon_isScEtaEE);
   fChain->SetBranchAddress("Photon_mvaID_WP80", Photon_mvaID_WP80, &b_Photon_mvaID_WP80);
   fChain->SetBranchAddress("Photon_mvaID_WP90", Photon_mvaID_WP90, &b_Photon_mvaID_WP90);
   fChain->SetBranchAddress("Photon_pixelSeed", Photon_pixelSeed, &b_Photon_pixelSeed);
   fChain->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt, &b_Pileup_nTrueInt);
   fChain->SetBranchAddress("Pileup_nPU", &Pileup_nPU, &b_Pileup_nPU);
   fChain->SetBranchAddress("Pileup_sumEOOT", &Pileup_sumEOOT, &b_Pileup_sumEOOT);
   fChain->SetBranchAddress("Pileup_sumLOOT", &Pileup_sumLOOT, &b_Pileup_sumLOOT);
   fChain->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi, &b_PuppiMET_phi);
   fChain->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt, &b_PuppiMET_pt);
   fChain->SetBranchAddress("PuppiMET_sumEt", &PuppiMET_sumEt, &b_PuppiMET_sumEt);
   fChain->SetBranchAddress("RawMET_phi", &RawMET_phi, &b_RawMET_phi);
   fChain->SetBranchAddress("RawMET_pt", &RawMET_pt, &b_RawMET_pt);
   fChain->SetBranchAddress("RawMET_sumEt", &RawMET_sumEt, &b_RawMET_sumEt);
   fChain->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, &b_fixedGridRhoFastjetCentralCalo);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, &b_fixedGridRhoFastjetCentralNeutral);
   fChain->SetBranchAddress("nGenDressedLepton", &nGenDressedLepton, &b_nGenDressedLepton);
   fChain->SetBranchAddress("GenDressedLepton_eta", GenDressedLepton_eta, &b_GenDressedLepton_eta);
   fChain->SetBranchAddress("GenDressedLepton_mass", GenDressedLepton_mass, &b_GenDressedLepton_mass);
   fChain->SetBranchAddress("GenDressedLepton_phi", GenDressedLepton_phi, &b_GenDressedLepton_phi);
   fChain->SetBranchAddress("GenDressedLepton_pt", GenDressedLepton_pt, &b_GenDressedLepton_pt);
   fChain->SetBranchAddress("GenDressedLepton_pdgId", GenDressedLepton_pdgId, &b_GenDressedLepton_pdgId);
   fChain->SetBranchAddress("nSoftActivityJet", &nSoftActivityJet, &b_nSoftActivityJet);
   fChain->SetBranchAddress("SoftActivityJet_eta", SoftActivityJet_eta, &b_SoftActivityJet_eta);
   fChain->SetBranchAddress("SoftActivityJet_phi", SoftActivityJet_phi, &b_SoftActivityJet_phi);
   fChain->SetBranchAddress("SoftActivityJet_pt", SoftActivityJet_pt, &b_SoftActivityJet_pt);
   fChain->SetBranchAddress("SoftActivityJetHT", &SoftActivityJetHT, &b_SoftActivityJetHT);
   fChain->SetBranchAddress("SoftActivityJetHT10", &SoftActivityJetHT10, &b_SoftActivityJetHT10);
   fChain->SetBranchAddress("SoftActivityJetHT2", &SoftActivityJetHT2, &b_SoftActivityJetHT2);
   fChain->SetBranchAddress("SoftActivityJetHT5", &SoftActivityJetHT5, &b_SoftActivityJetHT5);
   fChain->SetBranchAddress("SoftActivityJetNjets10", &SoftActivityJetNjets10, &b_SoftActivityJetNjets10);
   fChain->SetBranchAddress("SoftActivityJetNjets2", &SoftActivityJetNjets2, &b_SoftActivityJetNjets2);
   fChain->SetBranchAddress("SoftActivityJetNjets5", &SoftActivityJetNjets5, &b_SoftActivityJetNjets5);
   fChain->SetBranchAddress("nSubJet", &nSubJet, &b_nSubJet);
   fChain->SetBranchAddress("SubJet_btagCMVA", SubJet_btagCMVA, &b_SubJet_btagCMVA);
   fChain->SetBranchAddress("SubJet_btagCSVV2", SubJet_btagCSVV2, &b_SubJet_btagCSVV2);
   fChain->SetBranchAddress("SubJet_btagDeepB", SubJet_btagDeepB, &b_SubJet_btagDeepB);
   fChain->SetBranchAddress("SubJet_eta", SubJet_eta, &b_SubJet_eta);
   fChain->SetBranchAddress("SubJet_mass", SubJet_mass, &b_SubJet_mass);
   fChain->SetBranchAddress("SubJet_n2b1", SubJet_n2b1, &b_SubJet_n2b1);
   fChain->SetBranchAddress("SubJet_n3b1", SubJet_n3b1, &b_SubJet_n3b1);
   fChain->SetBranchAddress("SubJet_phi", SubJet_phi, &b_SubJet_phi);
   fChain->SetBranchAddress("SubJet_pt", SubJet_pt, &b_SubJet_pt);
   fChain->SetBranchAddress("SubJet_rawFactor", SubJet_rawFactor, &b_SubJet_rawFactor);
   fChain->SetBranchAddress("SubJet_tau1", SubJet_tau1, &b_SubJet_tau1);
   fChain->SetBranchAddress("SubJet_tau2", SubJet_tau2, &b_SubJet_tau2);
   fChain->SetBranchAddress("SubJet_tau3", SubJet_tau3, &b_SubJet_tau3);
   fChain->SetBranchAddress("SubJet_tau4", SubJet_tau4, &b_SubJet_tau4);
   fChain->SetBranchAddress("nTau", &nTau, &b_nTau);
   fChain->SetBranchAddress("Tau_chargedIso", Tau_chargedIso, &b_Tau_chargedIso);
   fChain->SetBranchAddress("Tau_dxy", Tau_dxy, &b_Tau_dxy);
   fChain->SetBranchAddress("Tau_dz", Tau_dz, &b_Tau_dz);
   fChain->SetBranchAddress("Tau_eta", Tau_eta, &b_Tau_eta);
   fChain->SetBranchAddress("Tau_leadTkDeltaEta", Tau_leadTkDeltaEta, &b_Tau_leadTkDeltaEta);
   fChain->SetBranchAddress("Tau_leadTkDeltaPhi", Tau_leadTkDeltaPhi, &b_Tau_leadTkDeltaPhi);
   fChain->SetBranchAddress("Tau_leadTkPtOverTauPt", Tau_leadTkPtOverTauPt, &b_Tau_leadTkPtOverTauPt);
   fChain->SetBranchAddress("Tau_mass", Tau_mass, &b_Tau_mass);
   fChain->SetBranchAddress("Tau_neutralIso", Tau_neutralIso, &b_Tau_neutralIso);
   fChain->SetBranchAddress("Tau_phi", Tau_phi, &b_Tau_phi);
   fChain->SetBranchAddress("Tau_photonsOutsideSignalCone", Tau_photonsOutsideSignalCone, &b_Tau_photonsOutsideSignalCone);
   fChain->SetBranchAddress("Tau_pt", Tau_pt, &b_Tau_pt);
   fChain->SetBranchAddress("Tau_puCorr", Tau_puCorr, &b_Tau_puCorr);
   fChain->SetBranchAddress("Tau_rawAntiEle", Tau_rawAntiEle, &b_Tau_rawAntiEle);
   fChain->SetBranchAddress("Tau_rawIso", Tau_rawIso, &b_Tau_rawIso);
   fChain->SetBranchAddress("Tau_rawIsodR03", Tau_rawIsodR03, &b_Tau_rawIsodR03);
   fChain->SetBranchAddress("Tau_rawMVAnewDM2017v2", Tau_rawMVAnewDM2017v2, &b_Tau_rawMVAnewDM2017v2);
   fChain->SetBranchAddress("Tau_rawMVAoldDM", Tau_rawMVAoldDM, &b_Tau_rawMVAoldDM);
   fChain->SetBranchAddress("Tau_rawMVAoldDM2017v1", Tau_rawMVAoldDM2017v1, &b_Tau_rawMVAoldDM2017v1);
   fChain->SetBranchAddress("Tau_rawMVAoldDM2017v2", Tau_rawMVAoldDM2017v2, &b_Tau_rawMVAoldDM2017v2);
   fChain->SetBranchAddress("Tau_rawMVAoldDMdR032017v2", Tau_rawMVAoldDMdR032017v2, &b_Tau_rawMVAoldDMdR032017v2);
   fChain->SetBranchAddress("Tau_charge", Tau_charge, &b_Tau_charge);
   fChain->SetBranchAddress("Tau_decayMode", Tau_decayMode, &b_Tau_decayMode);
   fChain->SetBranchAddress("Tau_jetIdx", Tau_jetIdx, &b_Tau_jetIdx);
   fChain->SetBranchAddress("Tau_rawAntiEleCat", Tau_rawAntiEleCat, &b_Tau_rawAntiEleCat);
   fChain->SetBranchAddress("Tau_idAntiEle", Tau_idAntiEle, &b_Tau_idAntiEle);
   fChain->SetBranchAddress("Tau_idAntiMu", Tau_idAntiMu, &b_Tau_idAntiMu);
   fChain->SetBranchAddress("Tau_idDecayMode", Tau_idDecayMode, &b_Tau_idDecayMode);
   fChain->SetBranchAddress("Tau_idDecayModeNewDMs", Tau_idDecayModeNewDMs, &b_Tau_idDecayModeNewDMs);
   fChain->SetBranchAddress("Tau_idMVAnewDM2017v2", Tau_idMVAnewDM2017v2, &b_Tau_idMVAnewDM2017v2);
   fChain->SetBranchAddress("Tau_idMVAoldDM", Tau_idMVAoldDM, &b_Tau_idMVAoldDM);
   fChain->SetBranchAddress("Tau_idMVAoldDM2017v1", Tau_idMVAoldDM2017v1, &b_Tau_idMVAoldDM2017v1);
   fChain->SetBranchAddress("Tau_idMVAoldDM2017v2", Tau_idMVAoldDM2017v2, &b_Tau_idMVAoldDM2017v2);
   fChain->SetBranchAddress("Tau_idMVAoldDMdR032017v2", Tau_idMVAoldDMdR032017v2, &b_Tau_idMVAoldDMdR032017v2);
   fChain->SetBranchAddress("TkMET_phi", &TkMET_phi, &b_TkMET_phi);
   fChain->SetBranchAddress("TkMET_pt", &TkMET_pt, &b_TkMET_pt);
   fChain->SetBranchAddress("TkMET_sumEt", &TkMET_sumEt, &b_TkMET_sumEt);
   fChain->SetBranchAddress("nTrigObj", &nTrigObj, &b_nTrigObj);
   fChain->SetBranchAddress("TrigObj_pt", TrigObj_pt, &b_TrigObj_pt);
   fChain->SetBranchAddress("TrigObj_eta", TrigObj_eta, &b_TrigObj_eta);
   fChain->SetBranchAddress("TrigObj_phi", TrigObj_phi, &b_TrigObj_phi);
   fChain->SetBranchAddress("TrigObj_l1pt", TrigObj_l1pt, &b_TrigObj_l1pt);
   fChain->SetBranchAddress("TrigObj_l1pt_2", TrigObj_l1pt_2, &b_TrigObj_l1pt_2);
   fChain->SetBranchAddress("TrigObj_l2pt", TrigObj_l2pt, &b_TrigObj_l2pt);
   fChain->SetBranchAddress("TrigObj_id", TrigObj_id, &b_TrigObj_id);
   fChain->SetBranchAddress("TrigObj_l1iso", TrigObj_l1iso, &b_TrigObj_l1iso);
   fChain->SetBranchAddress("TrigObj_l1charge", TrigObj_l1charge, &b_TrigObj_l1charge);
   fChain->SetBranchAddress("TrigObj_filterBits", TrigObj_filterBits, &b_TrigObj_filterBits);
   fChain->SetBranchAddress("genTtbarId", &genTtbarId, &b_genTtbarId);
   fChain->SetBranchAddress("nOtherPV", &nOtherPV, &b_nOtherPV);
   fChain->SetBranchAddress("OtherPV_z", OtherPV_z, &b_OtherPV_z);
   fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   fChain->SetBranchAddress("nSV", &nSV, &b_nSV);
   fChain->SetBranchAddress("SV_dlen", SV_dlen, &b_SV_dlen);
   fChain->SetBranchAddress("SV_dlenSig", SV_dlenSig, &b_SV_dlenSig);
   fChain->SetBranchAddress("SV_pAngle", SV_pAngle, &b_SV_pAngle);
   fChain->SetBranchAddress("Electron_genPartIdx", Electron_genPartIdx, &b_Electron_genPartIdx);
   fChain->SetBranchAddress("Electron_genPartFlav", Electron_genPartFlav, &b_Electron_genPartFlav);
   fChain->SetBranchAddress("GenJetAK8_partonFlavour", GenJetAK8_partonFlavour, &b_GenJetAK8_partonFlavour);
   fChain->SetBranchAddress("GenJetAK8_hadronFlavour", GenJetAK8_hadronFlavour, &b_GenJetAK8_hadronFlavour);
   fChain->SetBranchAddress("GenJet_partonFlavour", GenJet_partonFlavour, &b_GenJet_partonFlavour);
   fChain->SetBranchAddress("GenJet_hadronFlavour", GenJet_hadronFlavour, &b_GenJet_hadronFlavour);
   fChain->SetBranchAddress("Jet_genJetIdx", Jet_genJetIdx, &b_Jet_genJetIdx);
   fChain->SetBranchAddress("Jet_hadronFlavour", Jet_hadronFlavour, &b_Jet_hadronFlavour);
   fChain->SetBranchAddress("Jet_partonFlavour", Jet_partonFlavour, &b_Jet_partonFlavour);
   fChain->SetBranchAddress("Muon_genPartIdx", Muon_genPartIdx, &b_Muon_genPartIdx);
   fChain->SetBranchAddress("Muon_genPartFlav", Muon_genPartFlav, &b_Muon_genPartFlav);
   fChain->SetBranchAddress("Photon_genPartIdx", Photon_genPartIdx, &b_Photon_genPartIdx);
   fChain->SetBranchAddress("Photon_genPartFlav", Photon_genPartFlav, &b_Photon_genPartFlav);
   fChain->SetBranchAddress("MET_fiducialGenPhi", &MET_fiducialGenPhi, &b_MET_fiducialGenPhi);
   fChain->SetBranchAddress("MET_fiducialGenPt", &MET_fiducialGenPt, &b_MET_fiducialGenPt);
   fChain->SetBranchAddress("Electron_cleanmask", Electron_cleanmask, &b_Electron_cleanmask);
   fChain->SetBranchAddress("Jet_cleanmask", Jet_cleanmask, &b_Jet_cleanmask);
   fChain->SetBranchAddress("Muon_cleanmask", Muon_cleanmask, &b_Muon_cleanmask);
   fChain->SetBranchAddress("Photon_cleanmask", Photon_cleanmask, &b_Photon_cleanmask);
   fChain->SetBranchAddress("Tau_cleanmask", Tau_cleanmask, &b_Tau_cleanmask);
   fChain->SetBranchAddress("SV_chi2", SV_chi2, &b_SV_chi2);
   fChain->SetBranchAddress("SV_eta", SV_eta, &b_SV_eta);
   fChain->SetBranchAddress("SV_mass", SV_mass, &b_SV_mass);
   fChain->SetBranchAddress("SV_ndof", SV_ndof, &b_SV_ndof);
   fChain->SetBranchAddress("SV_phi", SV_phi, &b_SV_phi);
   fChain->SetBranchAddress("SV_pt", SV_pt, &b_SV_pt);
   fChain->SetBranchAddress("SV_x", SV_x, &b_SV_x);
   fChain->SetBranchAddress("SV_y", SV_y, &b_SV_y);
   fChain->SetBranchAddress("SV_z", SV_z, &b_SV_z);
   fChain->SetBranchAddress("Tau_genPartIdx", Tau_genPartIdx, &b_Tau_genPartIdx);
   fChain->SetBranchAddress("Tau_genPartFlav", Tau_genPartFlav, &b_Tau_genPartFlav);
   fChain->SetBranchAddress("L1simulation_step", &L1simulation_step, &b_L1simulation_step);
   fChain->SetBranchAddress("HLTriggerFirstPath", &HLTriggerFirstPath, &b_HLTriggerFirstPath);
   fChain->SetBranchAddress("HLT_AK8PFJet360_TrimMass30", &HLT_AK8PFJet360_TrimMass30, &b_HLT_AK8PFJet360_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8PFJet380_TrimMass30", &HLT_AK8PFJet380_TrimMass30, &b_HLT_AK8PFJet380_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8PFJet400_TrimMass30", &HLT_AK8PFJet400_TrimMass30, &b_HLT_AK8PFJet400_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8PFJet420_TrimMass30", &HLT_AK8PFJet420_TrimMass30, &b_HLT_AK8PFJet420_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50", &HLT_AK8PFHT700_TrimR0p1PT0p03Mass50, &b_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50);
   fChain->SetBranchAddress("HLT_AK8PFHT750_TrimMass50", &HLT_AK8PFHT750_TrimMass50, &b_HLT_AK8PFHT750_TrimMass50);
   fChain->SetBranchAddress("HLT_AK8PFHT800_TrimMass50", &HLT_AK8PFHT800_TrimMass50, &b_HLT_AK8PFHT800_TrimMass50);
   fChain->SetBranchAddress("HLT_AK8PFHT850_TrimMass50", &HLT_AK8PFHT850_TrimMass50, &b_HLT_AK8PFHT850_TrimMass50);
   fChain->SetBranchAddress("HLT_AK8PFHT900_TrimMass50", &HLT_AK8PFHT900_TrimMass50, &b_HLT_AK8PFHT900_TrimMass50);
   fChain->SetBranchAddress("HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20", &HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20, &b_HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20);
   fChain->SetBranchAddress("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20", &HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20, &b_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20);
   fChain->SetBranchAddress("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087", &HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087, &b_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087);
   fChain->SetBranchAddress("HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087", &HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087, &b_HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087);
   fChain->SetBranchAddress("HLT_AK8DiPFJet300_200_TrimMass30", &HLT_AK8DiPFJet300_200_TrimMass30, &b_HLT_AK8DiPFJet300_200_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8DiPFJet280_200_TrimMass30", &HLT_AK8DiPFJet280_200_TrimMass30, &b_HLT_AK8DiPFJet280_200_TrimMass30);
   fChain->SetBranchAddress("HLT_CaloJet500_NoJetID", &HLT_CaloJet500_NoJetID, &b_HLT_CaloJet500_NoJetID);
   fChain->SetBranchAddress("HLT_CaloJet550_NoJetID", &HLT_CaloJet550_NoJetID, &b_HLT_CaloJet550_NoJetID);
   fChain->SetBranchAddress("HLT_Trimuon5_3p5_2_Upsilon_Muon", &HLT_Trimuon5_3p5_2_Upsilon_Muon, &b_HLT_Trimuon5_3p5_2_Upsilon_Muon);
   fChain->SetBranchAddress("HLT_DoubleEle25_CaloIdL_MW", &HLT_DoubleEle25_CaloIdL_MW, &b_HLT_DoubleEle25_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_DoubleEle27_CaloIdL_MW", &HLT_DoubleEle27_CaloIdL_MW, &b_HLT_DoubleEle27_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_DoubleEle33_CaloIdL_MW", &HLT_DoubleEle33_CaloIdL_MW, &b_HLT_DoubleEle33_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_DoubleEle24_eta2p1_WPTight_Gsf", &HLT_DoubleEle24_eta2p1_WPTight_Gsf, &b_HLT_DoubleEle24_eta2p1_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350, &b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350);
   fChain->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350, &b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350);
   fChain->SetBranchAddress("HLT_Ele27_Ele37_CaloIdL_MW", &HLT_Ele27_Ele37_CaloIdL_MW, &b_HLT_Ele27_Ele37_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_Mu27_Ele37_CaloIdL_MW", &HLT_Mu27_Ele37_CaloIdL_MW, &b_HLT_Mu27_Ele37_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_Mu37_Ele27_CaloIdL_MW", &HLT_Mu37_Ele27_CaloIdL_MW, &b_HLT_Mu37_Ele27_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_Mu37_TkMu27", &HLT_Mu37_TkMu27, &b_HLT_Mu37_TkMu27);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Bs", &HLT_DoubleMu4_3_Bs, &b_HLT_DoubleMu4_3_Bs);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Jpsi_Displaced", &HLT_DoubleMu4_3_Jpsi_Displaced, &b_HLT_DoubleMu4_3_Jpsi_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu4_JpsiTrk_Displaced", &HLT_DoubleMu4_JpsiTrk_Displaced, &b_HLT_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced", &HLT_DoubleMu4_LowMassNonResonantTrk_Displaced, &b_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu3_Trk_Tau3mu", &HLT_DoubleMu3_Trk_Tau3mu, &b_HLT_DoubleMu3_Trk_Tau3mu);
   fChain->SetBranchAddress("HLT_DoubleMu4_PsiPrimeTrk_Displaced", &HLT_DoubleMu4_PsiPrimeTrk_Displaced, &b_HLT_DoubleMu4_PsiPrimeTrk_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu4_Mass8_DZ_PFHT350", &HLT_DoubleMu4_Mass8_DZ_PFHT350, &b_HLT_DoubleMu4_Mass8_DZ_PFHT350);
   fChain->SetBranchAddress("HLT_DoubleMu8_Mass8_PFHT350", &HLT_DoubleMu8_Mass8_PFHT350, &b_HLT_DoubleMu8_Mass8_PFHT350);
   fChain->SetBranchAddress("HLT_Mu3_PFJet40", &HLT_Mu3_PFJet40, &b_HLT_Mu3_PFJet40);
   fChain->SetBranchAddress("HLT_Mu7p5_L2Mu2_Jpsi", &HLT_Mu7p5_L2Mu2_Jpsi, &b_HLT_Mu7p5_L2Mu2_Jpsi);
   fChain->SetBranchAddress("HLT_Mu7p5_L2Mu2_Upsilon", &HLT_Mu7p5_L2Mu2_Upsilon, &b_HLT_Mu7p5_L2Mu2_Upsilon);
   fChain->SetBranchAddress("HLT_Mu7p5_Track2_Jpsi", &HLT_Mu7p5_Track2_Jpsi, &b_HLT_Mu7p5_Track2_Jpsi);
   fChain->SetBranchAddress("HLT_Mu7p5_Track3p5_Jpsi", &HLT_Mu7p5_Track3p5_Jpsi, &b_HLT_Mu7p5_Track3p5_Jpsi);
   fChain->SetBranchAddress("HLT_Mu7p5_Track7_Jpsi", &HLT_Mu7p5_Track7_Jpsi, &b_HLT_Mu7p5_Track7_Jpsi);
   fChain->SetBranchAddress("HLT_Mu7p5_Track2_Upsilon", &HLT_Mu7p5_Track2_Upsilon, &b_HLT_Mu7p5_Track2_Upsilon);
   fChain->SetBranchAddress("HLT_Mu7p5_Track3p5_Upsilon", &HLT_Mu7p5_Track3p5_Upsilon, &b_HLT_Mu7p5_Track3p5_Upsilon);
   fChain->SetBranchAddress("HLT_Mu7p5_Track7_Upsilon", &HLT_Mu7p5_Track7_Upsilon, &b_HLT_Mu7p5_Track7_Upsilon);
   fChain->SetBranchAddress("HLT_DoublePhoton33_CaloIdL", &HLT_DoublePhoton33_CaloIdL, &b_HLT_DoublePhoton33_CaloIdL);
   fChain->SetBranchAddress("HLT_DoublePhoton70", &HLT_DoublePhoton70, &b_HLT_DoublePhoton70);
   fChain->SetBranchAddress("HLT_DoublePhoton85", &HLT_DoublePhoton85, &b_HLT_DoublePhoton85);
   fChain->SetBranchAddress("HLT_Ele20_WPTight_Gsf", &HLT_Ele20_WPTight_Gsf, &b_HLT_Ele20_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele20_WPLoose_Gsf", &HLT_Ele20_WPLoose_Gsf, &b_HLT_Ele20_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_Ele20_eta2p1_WPLoose_Gsf", &HLT_Ele20_eta2p1_WPLoose_Gsf, &b_HLT_Ele20_eta2p1_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_DiEle27_WPTightCaloOnly_L1DoubleEG", &HLT_DiEle27_WPTightCaloOnly_L1DoubleEG, &b_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG);
   fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_Ele27_WPTight_Gsf, &b_HLT_Ele27_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf, &b_HLT_Ele32_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf", &HLT_Ele35_WPTight_Gsf, &b_HLT_Ele35_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf_L1EGMT", &HLT_Ele35_WPTight_Gsf_L1EGMT, &b_HLT_Ele35_WPTight_Gsf_L1EGMT);
   fChain->SetBranchAddress("HLT_Ele38_WPTight_Gsf", &HLT_Ele38_WPTight_Gsf, &b_HLT_Ele38_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele40_WPTight_Gsf", &HLT_Ele40_WPTight_Gsf, &b_HLT_Ele40_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf_L1DoubleEG", &HLT_Ele32_WPTight_Gsf_L1DoubleEG, &b_HLT_Ele32_WPTight_Gsf_L1DoubleEG);
   fChain->SetBranchAddress("HLT_HT450_Beamspot", &HLT_HT450_Beamspot, &b_HLT_HT450_Beamspot);
   fChain->SetBranchAddress("HLT_HT300_Beamspot", &HLT_HT300_Beamspot, &b_HLT_HT300_Beamspot);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1", &HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1, &b_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1", &HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1, &b_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1", &HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1, &b_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1", &HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1, &b_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1", &HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1, &b_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1", &HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1, &b_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1); 
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1", &HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1, &b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1, &b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1", &HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1, &b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_TightID_SingleL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_TightID_SingleL1, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_TightID_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_TightID_SingleL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_TightID_SingleL1, &b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_TightID_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu20", &HLT_IsoMu20, &b_HLT_IsoMu20);
   fChain->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24, &b_HLT_IsoMu24);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1", &HLT_IsoMu24_eta2p1, &b_HLT_IsoMu24_eta2p1);
   fChain->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27, &b_HLT_IsoMu27);
   fChain->SetBranchAddress("HLT_IsoMu30", &HLT_IsoMu30, &b_HLT_IsoMu30);
   fChain->SetBranchAddress("HLT_UncorrectedJetE30_NoBPTX", &HLT_UncorrectedJetE30_NoBPTX, &b_HLT_UncorrectedJetE30_NoBPTX);
   fChain->SetBranchAddress("HLT_UncorrectedJetE30_NoBPTX3BX", &HLT_UncorrectedJetE30_NoBPTX3BX, &b_HLT_UncorrectedJetE30_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_UncorrectedJetE60_NoBPTX3BX", &HLT_UncorrectedJetE60_NoBPTX3BX, &b_HLT_UncorrectedJetE60_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_UncorrectedJetE70_NoBPTX3BX", &HLT_UncorrectedJetE70_NoBPTX3BX, &b_HLT_UncorrectedJetE70_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_L1SingleMu18", &HLT_L1SingleMu18, &b_HLT_L1SingleMu18);
   fChain->SetBranchAddress("HLT_L1SingleMu25", &HLT_L1SingleMu25, &b_HLT_L1SingleMu25);
   fChain->SetBranchAddress("HLT_L2Mu10", &HLT_L2Mu10, &b_HLT_L2Mu10);
   fChain->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX3BX", &HLT_L2Mu10_NoVertex_NoBPTX3BX, &b_HLT_L2Mu10_NoVertex_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX", &HLT_L2Mu10_NoVertex_NoBPTX, &b_HLT_L2Mu10_NoVertex_NoBPTX);
   fChain->SetBranchAddress("HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX, &b_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX, &b_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_L2Mu50", &HLT_L2Mu50, &b_HLT_L2Mu50);
   fChain->SetBranchAddress("HLT_DoubleL2Mu50", &HLT_DoubleL2Mu50, &b_HLT_DoubleL2Mu50);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8);
   fChain->SetBranchAddress("HLT_Mu25_TkMu0_Onia", &HLT_Mu25_TkMu0_Onia, &b_HLT_Mu25_TkMu0_Onia);
   fChain->SetBranchAddress("HLT_Mu30_TkMu0_Onia", &HLT_Mu30_TkMu0_Onia, &b_HLT_Mu30_TkMu0_Onia);
   fChain->SetBranchAddress("HLT_Mu20_TkMu0_Phi", &HLT_Mu20_TkMu0_Phi, &b_HLT_Mu20_TkMu0_Phi);
   fChain->SetBranchAddress("HLT_Mu25_TkMu0_Phi", &HLT_Mu25_TkMu0_Phi, &b_HLT_Mu25_TkMu0_Phi);
   fChain->SetBranchAddress("HLT_Mu20", &HLT_Mu20, &b_HLT_Mu20);
   fChain->SetBranchAddress("HLT_Mu27", &HLT_Mu27, &b_HLT_Mu27);
   fChain->SetBranchAddress("HLT_Mu50", &HLT_Mu50, &b_HLT_Mu50);
   fChain->SetBranchAddress("HLT_Mu55", &HLT_Mu55, &b_HLT_Mu55);
   fChain->SetBranchAddress("HLT_OldMu100", &HLT_OldMu100, &b_HLT_OldMu100);
   fChain->SetBranchAddress("HLT_TkMu100", &HLT_TkMu100, &b_HLT_TkMu100);
   fChain->SetBranchAddress("HLT_DiPFJet15_NoCaloMatched", &HLT_DiPFJet15_NoCaloMatched, &b_HLT_DiPFJet15_NoCaloMatched);
   fChain->SetBranchAddress("HLT_DiPFJet25_NoCaloMatched", &HLT_DiPFJet25_NoCaloMatched, &b_HLT_DiPFJet25_NoCaloMatched);
   fChain->SetBranchAddress("HLT_DiPFJet15_FBEta3_NoCaloMatched", &HLT_DiPFJet15_FBEta3_NoCaloMatched, &b_HLT_DiPFJet15_FBEta3_NoCaloMatched);
   fChain->SetBranchAddress("HLT_DiPFJet25_FBEta3_NoCaloMatched", &HLT_DiPFJet25_FBEta3_NoCaloMatched, &b_HLT_DiPFJet25_FBEta3_NoCaloMatched);
   fChain->SetBranchAddress("HLT_DiPFJetAve40", &HLT_DiPFJetAve40, &b_HLT_DiPFJetAve40);
   fChain->SetBranchAddress("HLT_DiPFJetAve60", &HLT_DiPFJetAve60, &b_HLT_DiPFJetAve60);
   fChain->SetBranchAddress("HLT_DiPFJetAve80", &HLT_DiPFJetAve80, &b_HLT_DiPFJetAve80);
   fChain->SetBranchAddress("HLT_DiPFJetAve140", &HLT_DiPFJetAve140, &b_HLT_DiPFJetAve140);
   fChain->SetBranchAddress("HLT_DiPFJetAve200", &HLT_DiPFJetAve200, &b_HLT_DiPFJetAve200);
   fChain->SetBranchAddress("HLT_DiPFJetAve260", &HLT_DiPFJetAve260, &b_HLT_DiPFJetAve260);
   fChain->SetBranchAddress("HLT_DiPFJetAve320", &HLT_DiPFJetAve320, &b_HLT_DiPFJetAve320);
   fChain->SetBranchAddress("HLT_DiPFJetAve400", &HLT_DiPFJetAve400, &b_HLT_DiPFJetAve400);
   fChain->SetBranchAddress("HLT_DiPFJetAve500", &HLT_DiPFJetAve500, &b_HLT_DiPFJetAve500);
   fChain->SetBranchAddress("HLT_DiPFJetAve15_HFJEC", &HLT_DiPFJetAve15_HFJEC, &b_HLT_DiPFJetAve15_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve25_HFJEC", &HLT_DiPFJetAve25_HFJEC, &b_HLT_DiPFJetAve25_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve35_HFJEC", &HLT_DiPFJetAve35_HFJEC, &b_HLT_DiPFJetAve35_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve60_HFJEC", &HLT_DiPFJetAve60_HFJEC, &b_HLT_DiPFJetAve60_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve80_HFJEC", &HLT_DiPFJetAve80_HFJEC, &b_HLT_DiPFJetAve80_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve100_HFJEC", &HLT_DiPFJetAve100_HFJEC, &b_HLT_DiPFJetAve100_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve160_HFJEC", &HLT_DiPFJetAve160_HFJEC, &b_HLT_DiPFJetAve160_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve220_HFJEC", &HLT_DiPFJetAve220_HFJEC, &b_HLT_DiPFJetAve220_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve300_HFJEC", &HLT_DiPFJetAve300_HFJEC, &b_HLT_DiPFJetAve300_HFJEC);
   fChain->SetBranchAddress("HLT_AK8PFJet40", &HLT_AK8PFJet40, &b_HLT_AK8PFJet40);
   fChain->SetBranchAddress("HLT_AK8PFJet60", &HLT_AK8PFJet60, &b_HLT_AK8PFJet60);
   fChain->SetBranchAddress("HLT_AK8PFJet80", &HLT_AK8PFJet80, &b_HLT_AK8PFJet80);
   fChain->SetBranchAddress("HLT_AK8PFJet140", &HLT_AK8PFJet140, &b_HLT_AK8PFJet140);
   fChain->SetBranchAddress("HLT_AK8PFJet200", &HLT_AK8PFJet200, &b_HLT_AK8PFJet200);
   fChain->SetBranchAddress("HLT_AK8PFJet260", &HLT_AK8PFJet260, &b_HLT_AK8PFJet260);
   fChain->SetBranchAddress("HLT_AK8PFJet320", &HLT_AK8PFJet320, &b_HLT_AK8PFJet320);
   fChain->SetBranchAddress("HLT_AK8PFJet400", &HLT_AK8PFJet400, &b_HLT_AK8PFJet400);
   fChain->SetBranchAddress("HLT_AK8PFJet450", &HLT_AK8PFJet450, &b_HLT_AK8PFJet450);
   fChain->SetBranchAddress("HLT_AK8PFJet500", &HLT_AK8PFJet500, &b_HLT_AK8PFJet500);
   fChain->SetBranchAddress("HLT_AK8PFJet550", &HLT_AK8PFJet550, &b_HLT_AK8PFJet550);
   fChain->SetBranchAddress("HLT_PFJet40", &HLT_PFJet40, &b_HLT_PFJet40);
   fChain->SetBranchAddress("HLT_PFJet60", &HLT_PFJet60, &b_HLT_PFJet60);
   fChain->SetBranchAddress("HLT_PFJet80", &HLT_PFJet80, &b_HLT_PFJet80);
   fChain->SetBranchAddress("HLT_PFJet140", &HLT_PFJet140, &b_HLT_PFJet140);
   fChain->SetBranchAddress("HLT_PFJet200", &HLT_PFJet200, &b_HLT_PFJet200);
   fChain->SetBranchAddress("HLT_PFJet260", &HLT_PFJet260, &b_HLT_PFJet260);
   fChain->SetBranchAddress("HLT_PFJet320", &HLT_PFJet320, &b_HLT_PFJet320);
   fChain->SetBranchAddress("HLT_PFJet400", &HLT_PFJet400, &b_HLT_PFJet400);
   fChain->SetBranchAddress("HLT_PFJet450", &HLT_PFJet450, &b_HLT_PFJet450);
   fChain->SetBranchAddress("HLT_PFJet500", &HLT_PFJet500, &b_HLT_PFJet500);
   fChain->SetBranchAddress("HLT_PFJet550", &HLT_PFJet550, &b_HLT_PFJet550);
   fChain->SetBranchAddress("HLT_PFJetFwd40", &HLT_PFJetFwd40, &b_HLT_PFJetFwd40);
   fChain->SetBranchAddress("HLT_PFJetFwd60", &HLT_PFJetFwd60, &b_HLT_PFJetFwd60);
   fChain->SetBranchAddress("HLT_PFJetFwd80", &HLT_PFJetFwd80, &b_HLT_PFJetFwd80);
   fChain->SetBranchAddress("HLT_PFJetFwd140", &HLT_PFJetFwd140, &b_HLT_PFJetFwd140);
   fChain->SetBranchAddress("HLT_PFJetFwd200", &HLT_PFJetFwd200, &b_HLT_PFJetFwd200);
   fChain->SetBranchAddress("HLT_PFJetFwd260", &HLT_PFJetFwd260, &b_HLT_PFJetFwd260);
   fChain->SetBranchAddress("HLT_PFJetFwd320", &HLT_PFJetFwd320, &b_HLT_PFJetFwd320);
   fChain->SetBranchAddress("HLT_PFJetFwd400", &HLT_PFJetFwd400, &b_HLT_PFJetFwd400);
   fChain->SetBranchAddress("HLT_PFJetFwd450", &HLT_PFJetFwd450, &b_HLT_PFJetFwd450);
   fChain->SetBranchAddress("HLT_PFJetFwd500", &HLT_PFJetFwd500, &b_HLT_PFJetFwd500);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd40", &HLT_AK8PFJetFwd40, &b_HLT_AK8PFJetFwd40);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd60", &HLT_AK8PFJetFwd60, &b_HLT_AK8PFJetFwd60);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd80", &HLT_AK8PFJetFwd80, &b_HLT_AK8PFJetFwd80);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd140", &HLT_AK8PFJetFwd140, &b_HLT_AK8PFJetFwd140);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd200", &HLT_AK8PFJetFwd200, &b_HLT_AK8PFJetFwd200);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd260", &HLT_AK8PFJetFwd260, &b_HLT_AK8PFJetFwd260);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd320", &HLT_AK8PFJetFwd320, &b_HLT_AK8PFJetFwd320);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd400", &HLT_AK8PFJetFwd400, &b_HLT_AK8PFJetFwd400);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd450", &HLT_AK8PFJetFwd450, &b_HLT_AK8PFJetFwd450);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd500", &HLT_AK8PFJetFwd500, &b_HLT_AK8PFJetFwd500);
   fChain->SetBranchAddress("HLT_PFHT180", &HLT_PFHT180, &b_HLT_PFHT180);
   fChain->SetBranchAddress("HLT_PFHT250", &HLT_PFHT250, &b_HLT_PFHT250);
   fChain->SetBranchAddress("HLT_PFHT370", &HLT_PFHT370, &b_HLT_PFHT370);
   fChain->SetBranchAddress("HLT_PFHT430", &HLT_PFHT430, &b_HLT_PFHT430);
   fChain->SetBranchAddress("HLT_PFHT510", &HLT_PFHT510, &b_HLT_PFHT510);
   fChain->SetBranchAddress("HLT_PFHT590", &HLT_PFHT590, &b_HLT_PFHT590);
   fChain->SetBranchAddress("HLT_PFHT680", &HLT_PFHT680, &b_HLT_PFHT680);
   fChain->SetBranchAddress("HLT_PFHT780", &HLT_PFHT780, &b_HLT_PFHT780);
   fChain->SetBranchAddress("HLT_PFHT890", &HLT_PFHT890, &b_HLT_PFHT890);
   fChain->SetBranchAddress("HLT_PFHT800", &HLT_PFHT800, &b_HLT_PFHT800);
   fChain->SetBranchAddress("HLT_PFHT900", &HLT_PFHT900, &b_HLT_PFHT900);
   fChain->SetBranchAddress("HLT_PFHT1050", &HLT_PFHT1050, &b_HLT_PFHT1050);
   fChain->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight", &HLT_PFHT500_PFMET100_PFMHT100_IDTight, &b_HLT_PFHT500_PFMET100_PFMHT100_IDTight);
   fChain->SetBranchAddress("HLT_PFHT500_PFMET110_PFMHT110_IDTight", &HLT_PFHT500_PFMET110_PFMHT110_IDTight, &b_HLT_PFHT500_PFMET110_PFMHT110_IDTight);
   fChain->SetBranchAddress("HLT_PFHT700_PFMET85_PFMHT85_IDTight", &HLT_PFHT700_PFMET85_PFMHT85_IDTight, &b_HLT_PFHT700_PFMET85_PFMHT85_IDTight);
   fChain->SetBranchAddress("HLT_PFHT700_PFMET95_PFMHT95_IDTight", &HLT_PFHT700_PFMET95_PFMHT95_IDTight, &b_HLT_PFHT700_PFMET95_PFMHT95_IDTight);
   fChain->SetBranchAddress("HLT_PFHT800_PFMET75_PFMHT75_IDTight", &HLT_PFHT800_PFMET75_PFMHT75_IDTight, &b_HLT_PFHT800_PFMET75_PFMHT75_IDTight);
   fChain->SetBranchAddress("HLT_PFHT800_PFMET85_PFMHT85_IDTight", &HLT_PFHT800_PFMET85_PFMHT85_IDTight, &b_HLT_PFHT800_PFMET85_PFMHT85_IDTight);
   fChain->SetBranchAddress("HLT_PFMET110_PFMHT110_IDTight", &HLT_PFMET110_PFMHT110_IDTight, &b_HLT_PFMET110_PFMHT110_IDTight);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight", &HLT_PFMET120_PFMHT120_IDTight, &b_HLT_PFMET120_PFMHT120_IDTight);
   fChain->SetBranchAddress("HLT_PFMET130_PFMHT130_IDTight", &HLT_PFMET130_PFMHT130_IDTight, &b_HLT_PFMET130_PFMHT130_IDTight);
   fChain->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight", &HLT_PFMET140_PFMHT140_IDTight, &b_HLT_PFMET140_PFMHT140_IDTight);
   fChain->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1", &HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1, &b_HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1);
   fChain->SetBranchAddress("HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1", &HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1, &b_HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1", &HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1, &b_HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1);
   fChain->SetBranchAddress("HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1", &HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1, &b_HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1);
   fChain->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1", &HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1, &b_HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_PFHT60", &HLT_PFMET120_PFMHT120_IDTight_PFHT60, &b_HLT_PFMET120_PFMHT120_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60", &HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60, &b_HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_PFMETTypeOne110_PFMHT110_IDTight", &HLT_PFMETTypeOne110_PFMHT110_IDTight, &b_HLT_PFMETTypeOne110_PFMHT110_IDTight);
   fChain->SetBranchAddress("HLT_PFMETTypeOne120_PFMHT120_IDTight", &HLT_PFMETTypeOne120_PFMHT120_IDTight, &b_HLT_PFMETTypeOne120_PFMHT120_IDTight);
   fChain->SetBranchAddress("HLT_PFMETTypeOne130_PFMHT130_IDTight", &HLT_PFMETTypeOne130_PFMHT130_IDTight, &b_HLT_PFMETTypeOne130_PFMHT130_IDTight);
   fChain->SetBranchAddress("HLT_PFMETTypeOne140_PFMHT140_IDTight", &HLT_PFMETTypeOne140_PFMHT140_IDTight, &b_HLT_PFMETTypeOne140_PFMHT140_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_PFMETNoMu110_PFMHTNoMu110_IDTight, &b_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight", &HLT_PFMETNoMu130_PFMHTNoMu130_IDTight, &b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight, &b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight);
   fChain->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight);
   fChain->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight);
   fChain->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight);
   fChain->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight);
   fChain->SetBranchAddress("HLT_L1ETMHadSeeds", &HLT_L1ETMHadSeeds, &b_HLT_L1ETMHadSeeds);
   fChain->SetBranchAddress("HLT_CaloMHT90", &HLT_CaloMHT90, &b_HLT_CaloMHT90);
   fChain->SetBranchAddress("HLT_CaloMET80_NotCleaned", &HLT_CaloMET80_NotCleaned, &b_HLT_CaloMET80_NotCleaned);
   fChain->SetBranchAddress("HLT_CaloMET90_NotCleaned", &HLT_CaloMET90_NotCleaned, &b_HLT_CaloMET90_NotCleaned);
   fChain->SetBranchAddress("HLT_CaloMET100_NotCleaned", &HLT_CaloMET100_NotCleaned, &b_HLT_CaloMET100_NotCleaned);
   fChain->SetBranchAddress("HLT_CaloMET110_NotCleaned", &HLT_CaloMET110_NotCleaned, &b_HLT_CaloMET110_NotCleaned);
   fChain->SetBranchAddress("HLT_CaloMET250_NotCleaned", &HLT_CaloMET250_NotCleaned, &b_HLT_CaloMET250_NotCleaned);
   fChain->SetBranchAddress("HLT_CaloMET70_HBHECleaned", &HLT_CaloMET70_HBHECleaned, &b_HLT_CaloMET70_HBHECleaned);
   fChain->SetBranchAddress("HLT_CaloMET80_HBHECleaned", &HLT_CaloMET80_HBHECleaned, &b_HLT_CaloMET80_HBHECleaned);
   fChain->SetBranchAddress("HLT_CaloMET90_HBHECleaned", &HLT_CaloMET90_HBHECleaned, &b_HLT_CaloMET90_HBHECleaned);
   fChain->SetBranchAddress("HLT_CaloMET100_HBHECleaned", &HLT_CaloMET100_HBHECleaned, &b_HLT_CaloMET100_HBHECleaned);
   fChain->SetBranchAddress("HLT_CaloMET250_HBHECleaned", &HLT_CaloMET250_HBHECleaned, &b_HLT_CaloMET250_HBHECleaned);
   fChain->SetBranchAddress("HLT_CaloMET300_HBHECleaned", &HLT_CaloMET300_HBHECleaned, &b_HLT_CaloMET300_HBHECleaned);
   fChain->SetBranchAddress("HLT_CaloMET350_HBHECleaned", &HLT_CaloMET350_HBHECleaned, &b_HLT_CaloMET350_HBHECleaned);
   fChain->SetBranchAddress("HLT_PFMET200_NotCleaned", &HLT_PFMET200_NotCleaned, &b_HLT_PFMET200_NotCleaned);
   fChain->SetBranchAddress("HLT_PFMET200_HBHECleaned", &HLT_PFMET200_HBHECleaned, &b_HLT_PFMET200_HBHECleaned);
   fChain->SetBranchAddress("HLT_PFMET250_HBHECleaned", &HLT_PFMET250_HBHECleaned, &b_HLT_PFMET250_HBHECleaned);
   fChain->SetBranchAddress("HLT_PFMET300_HBHECleaned", &HLT_PFMET300_HBHECleaned, &b_HLT_PFMET300_HBHECleaned);
   fChain->SetBranchAddress("HLT_PFMET200_HBHE_BeamHaloCleaned", &HLT_PFMET200_HBHE_BeamHaloCleaned, &b_HLT_PFMET200_HBHE_BeamHaloCleaned);
   fChain->SetBranchAddress("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned", &HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned, &b_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned);
   fChain->SetBranchAddress("HLT_MET105_IsoTrk50", &HLT_MET105_IsoTrk50, &b_HLT_MET105_IsoTrk50);
   fChain->SetBranchAddress("HLT_MET120_IsoTrk50", &HLT_MET120_IsoTrk50, &b_HLT_MET120_IsoTrk50);
   fChain->SetBranchAddress("HLT_SingleJet30_Mu12_SinglePFJet40", &HLT_SingleJet30_Mu12_SinglePFJet40, &b_HLT_SingleJet30_Mu12_SinglePFJet40);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets40_CaloBTagCSV_p33", &HLT_Mu12_DoublePFJets40_CaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets40_CaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets100_CaloBTagCSV_p33", &HLT_Mu12_DoublePFJets100_CaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets100_CaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets200_CaloBTagCSV_p33", &HLT_Mu12_DoublePFJets200_CaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets200_CaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets350_CaloBTagCSV_p33", &HLT_Mu12_DoublePFJets350_CaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets350_CaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33", &HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagCSV_p33", &HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagCSV_p33", &HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_DoublePFJets40_CaloBTagCSV_p33", &HLT_DoublePFJets40_CaloBTagCSV_p33, &b_HLT_DoublePFJets40_CaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_DoublePFJets100_CaloBTagCSV_p33", &HLT_DoublePFJets100_CaloBTagCSV_p33, &b_HLT_DoublePFJets100_CaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_DoublePFJets200_CaloBTagCSV_p33", &HLT_DoublePFJets200_CaloBTagCSV_p33, &b_HLT_DoublePFJets200_CaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_DoublePFJets350_CaloBTagCSV_p33", &HLT_DoublePFJets350_CaloBTagCSV_p33, &b_HLT_DoublePFJets350_CaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33", &HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33, &b_HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33", &HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33, &b_HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33", &HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33, &b_HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_Photon300_NoHE", &HLT_Photon300_NoHE, &b_HLT_Photon300_NoHE);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL", &HLT_Mu8_TrkIsoVVL, &b_HLT_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ, &b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ);
   fChain->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL, &b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ, &b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ);
   fChain->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350, &b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL", &HLT_Mu19_TrkIsoVVL, &b_HLT_Mu19_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet20_Mu5", &HLT_BTagMu_AK4DiJet20_Mu5, &b_HLT_BTagMu_AK4DiJet20_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet40_Mu5", &HLT_BTagMu_AK4DiJet40_Mu5, &b_HLT_BTagMu_AK4DiJet40_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet70_Mu5", &HLT_BTagMu_AK4DiJet70_Mu5, &b_HLT_BTagMu_AK4DiJet70_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet110_Mu5", &HLT_BTagMu_AK4DiJet110_Mu5, &b_HLT_BTagMu_AK4DiJet110_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet170_Mu5", &HLT_BTagMu_AK4DiJet170_Mu5, &b_HLT_BTagMu_AK4DiJet170_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4Jet300_Mu5", &HLT_BTagMu_AK4Jet300_Mu5, &b_HLT_BTagMu_AK4Jet300_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK8DiJet170_Mu5", &HLT_BTagMu_AK8DiJet170_Mu5, &b_HLT_BTagMu_AK8DiJet170_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK8Jet300_Mu5", &HLT_BTagMu_AK8Jet300_Mu5, &b_HLT_BTagMu_AK8Jet300_Mu5);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Mu12_DoublePhoton20", &HLT_Mu12_DoublePhoton20, &b_HLT_Mu12_DoublePhoton20);
   fChain->SetBranchAddress("HLT_TriplePhoton_20_20_20_CaloIdLV2", &HLT_TriplePhoton_20_20_20_CaloIdLV2, &b_HLT_TriplePhoton_20_20_20_CaloIdLV2);
   fChain->SetBranchAddress("HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL", &HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL, &b_HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL);
   fChain->SetBranchAddress("HLT_TriplePhoton_30_30_10_CaloIdLV2", &HLT_TriplePhoton_30_30_10_CaloIdLV2, &b_HLT_TriplePhoton_30_30_10_CaloIdLV2);
   fChain->SetBranchAddress("HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL", &HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL, &b_HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL);
   fChain->SetBranchAddress("HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL", &HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL, &b_HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL);
//   fChain->SetBranchAddress("HLT_Photon25", &HLT_Photon25, &b_HLT_Photon25);
   fChain->SetBranchAddress("HLT_Photon33", &HLT_Photon33, &b_HLT_Photon33);
   fChain->SetBranchAddress("HLT_Photon50", &HLT_Photon50, &b_HLT_Photon50);
   fChain->SetBranchAddress("HLT_Photon75", &HLT_Photon75, &b_HLT_Photon75);
   fChain->SetBranchAddress("HLT_Photon90", &HLT_Photon90, &b_HLT_Photon90);
   fChain->SetBranchAddress("HLT_Photon120", &HLT_Photon120, &b_HLT_Photon120);
   fChain->SetBranchAddress("HLT_Photon150", &HLT_Photon150, &b_HLT_Photon150);
   fChain->SetBranchAddress("HLT_Photon175", &HLT_Photon175, &b_HLT_Photon175);
   fChain->SetBranchAddress("HLT_Photon200", &HLT_Photon200, &b_HLT_Photon200);
   fChain->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM", &HLT_Photon50_R9Id90_HE10_IsoM, &b_HLT_Photon50_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM", &HLT_Photon75_R9Id90_HE10_IsoM, &b_HLT_Photon75_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon90_R9Id90_HE10_IsoM", &HLT_Photon90_R9Id90_HE10_IsoM, &b_HLT_Photon90_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon120_R9Id90_HE10_IsoM", &HLT_Photon120_R9Id90_HE10_IsoM, &b_HLT_Photon120_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon165_R9Id90_HE10_IsoM", &HLT_Photon165_R9Id90_HE10_IsoM, &b_HLT_Photon165_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon90_CaloIdL_PFHT700", &HLT_Photon90_CaloIdL_PFHT700, &b_HLT_Photon90_CaloIdL_PFHT700);
   fChain->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", &HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90, &b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90);
   fChain->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95", &HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95, &b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95);
   fChain->SetBranchAddress("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55", &HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55, &b_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55);
   fChain->SetBranchAddress("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55", &HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55, &b_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55);
   fChain->SetBranchAddress("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55", &HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55, &b_HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55);
   fChain->SetBranchAddress("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55", &HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55, &b_HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_L1_NoOS", &HLT_Dimuon0_Jpsi_L1_NoOS, &b_HLT_Dimuon0_Jpsi_L1_NoOS);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing_NoOS", &HLT_Dimuon0_Jpsi_NoVertexing_NoOS, &b_HLT_Dimuon0_Jpsi_NoVertexing_NoOS);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi", &HLT_Dimuon0_Jpsi, &b_HLT_Dimuon0_Jpsi);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing", &HLT_Dimuon0_Jpsi_NoVertexing, &b_HLT_Dimuon0_Jpsi_NoVertexing);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R", &HLT_Dimuon0_Jpsi_L1_4R_0er1p5R, &b_HLT_Dimuon0_Jpsi_L1_4R_0er1p5R);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R", &HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R, &b_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi3p5_Muon2", &HLT_Dimuon0_Jpsi3p5_Muon2, &b_HLT_Dimuon0_Jpsi3p5_Muon2);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5", &HLT_Dimuon0_Upsilon_L1_4p5, &b_HLT_Dimuon0_Upsilon_L1_4p5);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_5", &HLT_Dimuon0_Upsilon_L1_5, &b_HLT_Dimuon0_Upsilon_L1_5);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5NoOS", &HLT_Dimuon0_Upsilon_L1_4p5NoOS, &b_HLT_Dimuon0_Upsilon_L1_4p5NoOS);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5er2p0", &HLT_Dimuon0_Upsilon_L1_4p5er2p0, &b_HLT_Dimuon0_Upsilon_L1_4p5er2p0);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5er2p0M", &HLT_Dimuon0_Upsilon_L1_4p5er2p0M, &b_HLT_Dimuon0_Upsilon_L1_4p5er2p0M);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_NoVertexing", &HLT_Dimuon0_Upsilon_NoVertexing, &b_HLT_Dimuon0_Upsilon_NoVertexing);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_5M", &HLT_Dimuon0_Upsilon_L1_5M, &b_HLT_Dimuon0_Upsilon_L1_5M);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass_L1_0er1p5R", &HLT_Dimuon0_LowMass_L1_0er1p5R, &b_HLT_Dimuon0_LowMass_L1_0er1p5R);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass_L1_0er1p5", &HLT_Dimuon0_LowMass_L1_0er1p5, &b_HLT_Dimuon0_LowMass_L1_0er1p5);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass", &HLT_Dimuon0_LowMass, &b_HLT_Dimuon0_LowMass);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass_L1_4", &HLT_Dimuon0_LowMass_L1_4, &b_HLT_Dimuon0_LowMass_L1_4);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass_L1_4R", &HLT_Dimuon0_LowMass_L1_4R, &b_HLT_Dimuon0_LowMass_L1_4R);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass_L1_TM530", &HLT_Dimuon0_LowMass_L1_TM530, &b_HLT_Dimuon0_LowMass_L1_TM530);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_L1_TM0", &HLT_Dimuon0_Upsilon_Muon_L1_TM0, &b_HLT_Dimuon0_Upsilon_Muon_L1_TM0);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_NoL1Mass", &HLT_Dimuon0_Upsilon_Muon_NoL1Mass, &b_HLT_Dimuon0_Upsilon_Muon_NoL1Mass);
   fChain->SetBranchAddress("HLT_TripleMu_5_3_3_Mass3p8to60_DZ", &HLT_TripleMu_5_3_3_Mass3p8to60_DZ, &b_HLT_TripleMu_5_3_3_Mass3p8to60_DZ);
   fChain->SetBranchAddress("HLT_TripleMu_10_5_5_DZ", &HLT_TripleMu_10_5_5_DZ, &b_HLT_TripleMu_10_5_5_DZ);
   fChain->SetBranchAddress("HLT_TripleMu_12_10_5", &HLT_TripleMu_12_10_5, &b_HLT_TripleMu_12_10_5);
   fChain->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15);
   fChain->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1);
   fChain->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15);
   fChain->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1);
   fChain->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET50_PFMHT60", &HLT_DoubleMu3_DZ_PFMET50_PFMHT60, &b_HLT_DoubleMu3_DZ_PFMET50_PFMHT60);
   fChain->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET70_PFMHT70", &HLT_DoubleMu3_DZ_PFMET70_PFMHT70, &b_HLT_DoubleMu3_DZ_PFMET70_PFMHT70);
   fChain->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET90_PFMHT90", &HLT_DoubleMu3_DZ_PFMET90_PFMHT90, &b_HLT_DoubleMu3_DZ_PFMET90_PFMHT90);
   fChain->SetBranchAddress("HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass", &HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass, &b_HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass);
   fChain->SetBranchAddress("HLT_DoubleMu4_Jpsi_Displaced", &HLT_DoubleMu4_Jpsi_Displaced, &b_HLT_DoubleMu4_Jpsi_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu4_Jpsi_NoVertexing", &HLT_DoubleMu4_Jpsi_NoVertexing, &b_HLT_DoubleMu4_Jpsi_NoVertexing);
   fChain->SetBranchAddress("HLT_DoubleMu4_JpsiTrkTrk_Displaced", &HLT_DoubleMu4_JpsiTrkTrk_Displaced, &b_HLT_DoubleMu4_JpsiTrkTrk_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu43NoFiltersNoVtx", &HLT_DoubleMu43NoFiltersNoVtx, &b_HLT_DoubleMu43NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_DoubleMu48NoFiltersNoVtx", &HLT_DoubleMu48NoFiltersNoVtx, &b_HLT_DoubleMu48NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL", &HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL, &b_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL);
   fChain->SetBranchAddress("HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL", &HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL, &b_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL);
   fChain->SetBranchAddress("HLT_DoubleMu20_7_Mass0to30_L1_DM4", &HLT_DoubleMu20_7_Mass0to30_L1_DM4, &b_HLT_DoubleMu20_7_Mass0to30_L1_DM4);
   fChain->SetBranchAddress("HLT_DoubleMu20_7_Mass0to30_L1_DM4EG", &HLT_DoubleMu20_7_Mass0to30_L1_DM4EG, &b_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG);
   fChain->SetBranchAddress("HLT_HT425", &HLT_HT425, &b_HLT_HT425);
   fChain->SetBranchAddress("HLT_HT430_DisplacedDijet40_DisplacedTrack", &HLT_HT430_DisplacedDijet40_DisplacedTrack, &b_HLT_HT430_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT430_DisplacedDijet60_DisplacedTrack", &HLT_HT430_DisplacedDijet60_DisplacedTrack, &b_HLT_HT430_DisplacedDijet60_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT430_DisplacedDijet80_DisplacedTrack", &HLT_HT430_DisplacedDijet80_DisplacedTrack, &b_HLT_HT430_DisplacedDijet80_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT400_DisplacedDijet40_DisplacedTrack", &HLT_HT400_DisplacedDijet40_DisplacedTrack, &b_HLT_HT400_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT650_DisplacedDijet60_Inclusive", &HLT_HT650_DisplacedDijet60_Inclusive, &b_HLT_HT650_DisplacedDijet60_Inclusive);
   fChain->SetBranchAddress("HLT_HT550_DisplacedDijet80_Inclusive", &HLT_HT550_DisplacedDijet80_Inclusive, &b_HLT_HT550_DisplacedDijet80_Inclusive);
   fChain->SetBranchAddress("HLT_HT550_DisplacedDijet60_Inclusive", &HLT_HT550_DisplacedDijet60_Inclusive, &b_HLT_HT550_DisplacedDijet60_Inclusive);
   fChain->SetBranchAddress("HLT_HT650_DisplacedDijet80_Inclusive", &HLT_HT650_DisplacedDijet80_Inclusive, &b_HLT_HT650_DisplacedDijet80_Inclusive);
   fChain->SetBranchAddress("HLT_HT750_DisplacedDijet80_Inclusive", &HLT_HT750_DisplacedDijet80_Inclusive, &b_HLT_HT750_DisplacedDijet80_Inclusive);
   fChain->SetBranchAddress("HLT_DiJet110_35_Mjj650_PFMET110", &HLT_DiJet110_35_Mjj650_PFMET110, &b_HLT_DiJet110_35_Mjj650_PFMET110);
   fChain->SetBranchAddress("HLT_DiJet110_35_Mjj650_PFMET120", &HLT_DiJet110_35_Mjj650_PFMET120, &b_HLT_DiJet110_35_Mjj650_PFMET120);
   fChain->SetBranchAddress("HLT_DiJet110_35_Mjj650_PFMET130", &HLT_DiJet110_35_Mjj650_PFMET130, &b_HLT_DiJet110_35_Mjj650_PFMET130);
   fChain->SetBranchAddress("HLT_TripleJet110_35_35_Mjj650_PFMET110", &HLT_TripleJet110_35_35_Mjj650_PFMET110, &b_HLT_TripleJet110_35_35_Mjj650_PFMET110);
   fChain->SetBranchAddress("HLT_TripleJet110_35_35_Mjj650_PFMET120", &HLT_TripleJet110_35_35_Mjj650_PFMET120, &b_HLT_TripleJet110_35_35_Mjj650_PFMET120);
   fChain->SetBranchAddress("HLT_TripleJet110_35_35_Mjj650_PFMET130", &HLT_TripleJet110_35_35_Mjj650_PFMET130, &b_HLT_TripleJet110_35_35_Mjj650_PFMET130);
   fChain->SetBranchAddress("HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg", &HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg, &b_HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg", &HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg, &b_HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg", &HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg, &b_HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned", &HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned, &b_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned);
   fChain->SetBranchAddress("HLT_Ele28_eta2p1_WPTight_Gsf_HT150", &HLT_Ele28_eta2p1_WPTight_Gsf_HT150, &b_HLT_Ele28_eta2p1_WPTight_Gsf_HT150);
   fChain->SetBranchAddress("HLT_Ele28_HighEta_SC20_Mass55", &HLT_Ele28_HighEta_SC20_Mass55, &b_HLT_Ele28_HighEta_SC20_Mass55);
   fChain->SetBranchAddress("HLT_DoubleMu20_7_Mass0to30_Photon23", &HLT_DoubleMu20_7_Mass0to30_Photon23, &b_HLT_DoubleMu20_7_Mass0to30_Photon23);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5", &HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5, &b_HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450_PFMET50", &HLT_Ele15_IsoVVVL_PFHT450_PFMET50, &b_HLT_Ele15_IsoVVVL_PFHT450_PFMET50);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450", &HLT_Ele15_IsoVVVL_PFHT450, &b_HLT_Ele15_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("HLT_Ele50_IsoVVVL_PFHT450", &HLT_Ele50_IsoVVVL_PFHT450, &b_HLT_Ele50_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT600", &HLT_Ele15_IsoVVVL_PFHT600, &b_HLT_Ele15_IsoVVVL_PFHT600);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60", &HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60, &b_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60);
   fChain->SetBranchAddress("HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60", &HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60, &b_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5", &HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5, &b_HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450_PFMET50", &HLT_Mu15_IsoVVVL_PFHT450_PFMET50, &b_HLT_Mu15_IsoVVVL_PFHT450_PFMET50);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450", &HLT_Mu15_IsoVVVL_PFHT450, &b_HLT_Mu15_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("HLT_Mu50_IsoVVVL_PFHT450", &HLT_Mu50_IsoVVVL_PFHT450, &b_HLT_Mu50_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT600", &HLT_Mu15_IsoVVVL_PFHT600, &b_HLT_Mu15_IsoVVVL_PFHT600);
   fChain->SetBranchAddress("HLT_Dimuon10_PsiPrime_Barrel_Seagulls", &HLT_Dimuon10_PsiPrime_Barrel_Seagulls, &b_HLT_Dimuon10_PsiPrime_Barrel_Seagulls);
   fChain->SetBranchAddress("HLT_Dimuon20_Jpsi_Barrel_Seagulls", &HLT_Dimuon20_Jpsi_Barrel_Seagulls, &b_HLT_Dimuon20_Jpsi_Barrel_Seagulls);
   fChain->SetBranchAddress("HLT_Dimuon10_Upsilon_Barrel_Seagulls", &HLT_Dimuon10_Upsilon_Barrel_Seagulls, &b_HLT_Dimuon10_Upsilon_Barrel_Seagulls);
   fChain->SetBranchAddress("HLT_Dimuon12_Upsilon_eta1p5", &HLT_Dimuon12_Upsilon_eta1p5, &b_HLT_Dimuon12_Upsilon_eta1p5);
   fChain->SetBranchAddress("HLT_Dimuon14_Phi_Barrel_Seagulls", &HLT_Dimuon14_Phi_Barrel_Seagulls, &b_HLT_Dimuon14_Phi_Barrel_Seagulls);
   fChain->SetBranchAddress("HLT_Dimuon18_PsiPrime", &HLT_Dimuon18_PsiPrime, &b_HLT_Dimuon18_PsiPrime);
   fChain->SetBranchAddress("HLT_Dimuon25_Jpsi", &HLT_Dimuon25_Jpsi, &b_HLT_Dimuon25_Jpsi);
   fChain->SetBranchAddress("HLT_Dimuon18_PsiPrime_noCorrL1", &HLT_Dimuon18_PsiPrime_noCorrL1, &b_HLT_Dimuon18_PsiPrime_noCorrL1);
   fChain->SetBranchAddress("HLT_Dimuon24_Upsilon_noCorrL1", &HLT_Dimuon24_Upsilon_noCorrL1, &b_HLT_Dimuon24_Upsilon_noCorrL1);
   fChain->SetBranchAddress("HLT_Dimuon24_Phi_noCorrL1", &HLT_Dimuon24_Phi_noCorrL1, &b_HLT_Dimuon24_Phi_noCorrL1);
   fChain->SetBranchAddress("HLT_Dimuon25_Jpsi_noCorrL1", &HLT_Dimuon25_Jpsi_noCorrL1, &b_HLT_Dimuon25_Jpsi_noCorrL1);
   fChain->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ, &b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ);
   fChain->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL, &b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT_DoubleIsoMu20_eta2p1", &HLT_DoubleIsoMu20_eta2p1, &b_HLT_DoubleIsoMu20_eta2p1);
   fChain->SetBranchAddress("HLT_DoubleIsoMu24_eta2p1", &HLT_DoubleIsoMu24_eta2p1, &b_HLT_DoubleIsoMu24_eta2p1);
   fChain->SetBranchAddress("HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx", &HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx, &b_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx", &HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx, &b_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx", &HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx, &b_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_Mu8", &HLT_Mu8, &b_HLT_Mu8);
   fChain->SetBranchAddress("HLT_Mu17", &HLT_Mu17, &b_HLT_Mu17);
   fChain->SetBranchAddress("HLT_Mu19", &HLT_Mu19, &b_HLT_Mu19);
   fChain->SetBranchAddress("HLT_Mu17_Photon30_IsoCaloId", &HLT_Mu17_Photon30_IsoCaloId, &b_HLT_Mu17_Photon30_IsoCaloId);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdM_TrackIdM_PFJet30", &HLT_Ele8_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele8_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdM_TrackIdM_PFJet30", &HLT_Ele17_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele17_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele23_CaloIdM_TrackIdM_PFJet30", &HLT_Ele23_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele23_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", &HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165, &b_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
   fChain->SetBranchAddress("HLT_Ele115_CaloIdVT_GsfTrkIdT", &HLT_Ele115_CaloIdVT_GsfTrkIdT, &b_HLT_Ele115_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_Ele135_CaloIdVT_GsfTrkIdT", &HLT_Ele135_CaloIdVT_GsfTrkIdT, &b_HLT_Ele135_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_Ele145_CaloIdVT_GsfTrkIdT", &HLT_Ele145_CaloIdVT_GsfTrkIdT, &b_HLT_Ele145_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_Ele200_CaloIdVT_GsfTrkIdT", &HLT_Ele200_CaloIdVT_GsfTrkIdT, &b_HLT_Ele200_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_Ele250_CaloIdVT_GsfTrkIdT", &HLT_Ele250_CaloIdVT_GsfTrkIdT, &b_HLT_Ele250_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_Ele300_CaloIdVT_GsfTrkIdT", &HLT_Ele300_CaloIdVT_GsfTrkIdT, &b_HLT_Ele300_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_PFHT300PT30_QuadPFJet_75_60_45_40", &HLT_PFHT300PT30_QuadPFJet_75_60_45_40, &b_HLT_PFHT300PT30_QuadPFJet_75_60_45_40);
   fChain->SetBranchAddress("HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0", &HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0, &b_HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0);
   fChain->SetBranchAddress("HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2", &HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2, &b_HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2);
   fChain->SetBranchAddress("HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2", &HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2, &b_HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2);
   fChain->SetBranchAddress("HLT_PFHT380_SixPFJet32", &HLT_PFHT380_SixPFJet32, &b_HLT_PFHT380_SixPFJet32);
   fChain->SetBranchAddress("HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5", &HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5, &b_HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5);
   fChain->SetBranchAddress("HLT_PFHT430_SixPFJet40", &HLT_PFHT430_SixPFJet40, &b_HLT_PFHT430_SixPFJet40);
   fChain->SetBranchAddress("HLT_PFHT350", &HLT_PFHT350, &b_HLT_PFHT350);
   fChain->SetBranchAddress("HLT_PFHT350MinPFJet15", &HLT_PFHT350MinPFJet15, &b_HLT_PFHT350MinPFJet15);
   fChain->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL", &HLT_Photon60_R9Id90_CaloIdL_IsoL, &b_HLT_Photon60_R9Id90_CaloIdL_IsoL);
   fChain->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL, &b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL);
   fChain->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15, &b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15);
   fChain->SetBranchAddress("HLT_FullTrack_Multiplicity85", &HLT_FullTrack_Multiplicity85, &b_HLT_FullTrack_Multiplicity85);
   fChain->SetBranchAddress("HLT_FullTrack_Multiplicity100", &HLT_FullTrack_Multiplicity100, &b_HLT_FullTrack_Multiplicity100);
   fChain->SetBranchAddress("HLT_FullTrack_Multiplicity130", &HLT_FullTrack_Multiplicity130, &b_HLT_FullTrack_Multiplicity130);
   fChain->SetBranchAddress("HLT_FullTrack_Multiplicity155", &HLT_FullTrack_Multiplicity155, &b_HLT_FullTrack_Multiplicity155);
   fChain->SetBranchAddress("HLT_ECALHT800", &HLT_ECALHT800, &b_HLT_ECALHT800);
   fChain->SetBranchAddress("HLT_DiSC30_18_EIso_AND_HE_Mass70", &HLT_DiSC30_18_EIso_AND_HE_Mass70, &b_HLT_DiSC30_18_EIso_AND_HE_Mass70);
   fChain->SetBranchAddress("HLT_Physics", &HLT_Physics, &b_HLT_Physics);
   fChain->SetBranchAddress("HLT_Physics_part0", &HLT_Physics_part0, &b_HLT_Physics_part0);
   fChain->SetBranchAddress("HLT_Physics_part1", &HLT_Physics_part1, &b_HLT_Physics_part1);
   fChain->SetBranchAddress("HLT_Physics_part2", &HLT_Physics_part2, &b_HLT_Physics_part2);
   fChain->SetBranchAddress("HLT_Physics_part3", &HLT_Physics_part3, &b_HLT_Physics_part3);
   fChain->SetBranchAddress("HLT_Physics_part4", &HLT_Physics_part4, &b_HLT_Physics_part4);
   fChain->SetBranchAddress("HLT_Physics_part5", &HLT_Physics_part5, &b_HLT_Physics_part5);
   fChain->SetBranchAddress("HLT_Physics_part6", &HLT_Physics_part6, &b_HLT_Physics_part6);
   fChain->SetBranchAddress("HLT_Physics_part7", &HLT_Physics_part7, &b_HLT_Physics_part7);
   fChain->SetBranchAddress("HLT_Random", &HLT_Random, &b_HLT_Random);
   fChain->SetBranchAddress("HLT_ZeroBias", &HLT_ZeroBias, &b_HLT_ZeroBias);
   fChain->SetBranchAddress("HLT_ZeroBias_part0", &HLT_ZeroBias_part0, &b_HLT_ZeroBias_part0);
   fChain->SetBranchAddress("HLT_ZeroBias_part1", &HLT_ZeroBias_part1, &b_HLT_ZeroBias_part1);
   fChain->SetBranchAddress("HLT_ZeroBias_part2", &HLT_ZeroBias_part2, &b_HLT_ZeroBias_part2);
   fChain->SetBranchAddress("HLT_ZeroBias_part3", &HLT_ZeroBias_part3, &b_HLT_ZeroBias_part3);
   fChain->SetBranchAddress("HLT_ZeroBias_part4", &HLT_ZeroBias_part4, &b_HLT_ZeroBias_part4);
   fChain->SetBranchAddress("HLT_ZeroBias_part5", &HLT_ZeroBias_part5, &b_HLT_ZeroBias_part5);
   fChain->SetBranchAddress("HLT_ZeroBias_part6", &HLT_ZeroBias_part6, &b_HLT_ZeroBias_part6);
   fChain->SetBranchAddress("HLT_ZeroBias_part7", &HLT_ZeroBias_part7, &b_HLT_ZeroBias_part7);
   fChain->SetBranchAddress("HLT_AK4CaloJet30", &HLT_AK4CaloJet30, &b_HLT_AK4CaloJet30);
   fChain->SetBranchAddress("HLT_AK4CaloJet40", &HLT_AK4CaloJet40, &b_HLT_AK4CaloJet40);
   fChain->SetBranchAddress("HLT_AK4CaloJet50", &HLT_AK4CaloJet50, &b_HLT_AK4CaloJet50);
   fChain->SetBranchAddress("HLT_AK4CaloJet80", &HLT_AK4CaloJet80, &b_HLT_AK4CaloJet80);
   fChain->SetBranchAddress("HLT_AK4CaloJet100", &HLT_AK4CaloJet100, &b_HLT_AK4CaloJet100);
   fChain->SetBranchAddress("HLT_AK4CaloJet120", &HLT_AK4CaloJet120, &b_HLT_AK4CaloJet120);
   fChain->SetBranchAddress("HLT_AK4PFJet30", &HLT_AK4PFJet30, &b_HLT_AK4PFJet30);
   fChain->SetBranchAddress("HLT_AK4PFJet50", &HLT_AK4PFJet50, &b_HLT_AK4PFJet50);
   fChain->SetBranchAddress("HLT_AK4PFJet80", &HLT_AK4PFJet80, &b_HLT_AK4PFJet80);
   fChain->SetBranchAddress("HLT_AK4PFJet100", &HLT_AK4PFJet100, &b_HLT_AK4PFJet100);
   fChain->SetBranchAddress("HLT_AK4PFJet120", &HLT_AK4PFJet120, &b_HLT_AK4PFJet120);
   fChain->SetBranchAddress("HLT_HISinglePhoton10_Eta3p1ForPPRef", &HLT_HISinglePhoton10_Eta3p1ForPPRef, &b_HLT_HISinglePhoton10_Eta3p1ForPPRef);
   fChain->SetBranchAddress("HLT_HISinglePhoton20_Eta3p1ForPPRef", &HLT_HISinglePhoton20_Eta3p1ForPPRef, &b_HLT_HISinglePhoton20_Eta3p1ForPPRef);
   fChain->SetBranchAddress("HLT_HISinglePhoton30_Eta3p1ForPPRef", &HLT_HISinglePhoton30_Eta3p1ForPPRef, &b_HLT_HISinglePhoton30_Eta3p1ForPPRef);
   fChain->SetBranchAddress("HLT_HISinglePhoton40_Eta3p1ForPPRef", &HLT_HISinglePhoton40_Eta3p1ForPPRef, &b_HLT_HISinglePhoton40_Eta3p1ForPPRef);
   fChain->SetBranchAddress("HLT_HISinglePhoton50_Eta3p1ForPPRef", &HLT_HISinglePhoton50_Eta3p1ForPPRef, &b_HLT_HISinglePhoton50_Eta3p1ForPPRef);
   fChain->SetBranchAddress("HLT_HISinglePhoton60_Eta3p1ForPPRef", &HLT_HISinglePhoton60_Eta3p1ForPPRef, &b_HLT_HISinglePhoton60_Eta3p1ForPPRef);
   fChain->SetBranchAddress("HLT_Photon20_HoverELoose", &HLT_Photon20_HoverELoose, &b_HLT_Photon20_HoverELoose);
   fChain->SetBranchAddress("HLT_Photon30_HoverELoose", &HLT_Photon30_HoverELoose, &b_HLT_Photon30_HoverELoose);
   fChain->SetBranchAddress("HLT_Photon40_HoverELoose", &HLT_Photon40_HoverELoose, &b_HLT_Photon40_HoverELoose);
   fChain->SetBranchAddress("HLT_Photon50_HoverELoose", &HLT_Photon50_HoverELoose, &b_HLT_Photon50_HoverELoose);
   fChain->SetBranchAddress("HLT_Photon60_HoverELoose", &HLT_Photon60_HoverELoose, &b_HLT_Photon60_HoverELoose);
   fChain->SetBranchAddress("HLT_EcalCalibration", &HLT_EcalCalibration, &b_HLT_EcalCalibration);
   fChain->SetBranchAddress("HLT_HcalCalibration", &HLT_HcalCalibration, &b_HLT_HcalCalibration);
   fChain->SetBranchAddress("HLT_L1UnpairedBunchBptxMinus", &HLT_L1UnpairedBunchBptxMinus, &b_HLT_L1UnpairedBunchBptxMinus);
   fChain->SetBranchAddress("HLT_L1UnpairedBunchBptxPlus", &HLT_L1UnpairedBunchBptxPlus, &b_HLT_L1UnpairedBunchBptxPlus);
   fChain->SetBranchAddress("HLT_L1NotBptxOR", &HLT_L1NotBptxOR, &b_HLT_L1NotBptxOR);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF_OR", &HLT_L1MinimumBiasHF_OR, &b_HLT_L1MinimumBiasHF_OR);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF0OR", &HLT_L1MinimumBiasHF0OR, &b_HLT_L1MinimumBiasHF0OR);
   fChain->SetBranchAddress("HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142, &b_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142);
   fChain->SetBranchAddress("HLT_HcalNZS", &HLT_HcalNZS, &b_HLT_HcalNZS);
   fChain->SetBranchAddress("HLT_HcalPhiSym", &HLT_HcalPhiSym, &b_HLT_HcalPhiSym);
   fChain->SetBranchAddress("HLT_HcalIsolatedbunch", &HLT_HcalIsolatedbunch, &b_HLT_HcalIsolatedbunch);
   fChain->SetBranchAddress("HLT_IsoTrackHB", &HLT_IsoTrackHB, &b_HLT_IsoTrackHB);
   fChain->SetBranchAddress("HLT_IsoTrackHE", &HLT_IsoTrackHE, &b_HLT_IsoTrackHE);
   fChain->SetBranchAddress("HLT_ZeroBias_FirstCollisionAfterAbortGap", &HLT_ZeroBias_FirstCollisionAfterAbortGap, &b_HLT_ZeroBias_FirstCollisionAfterAbortGap);
   fChain->SetBranchAddress("HLT_ZeroBias_IsolatedBunches", &HLT_ZeroBias_IsolatedBunches, &b_HLT_ZeroBias_IsolatedBunches);
   fChain->SetBranchAddress("HLT_ZeroBias_FirstCollisionInTrain", &HLT_ZeroBias_FirstCollisionInTrain, &b_HLT_ZeroBias_FirstCollisionInTrain);
   fChain->SetBranchAddress("HLT_ZeroBias_LastCollisionInTrain", &HLT_ZeroBias_LastCollisionInTrain, &b_HLT_ZeroBias_LastCollisionInTrain);
   fChain->SetBranchAddress("HLT_ZeroBias_FirstBXAfterTrain", &HLT_ZeroBias_FirstBXAfterTrain, &b_HLT_ZeroBias_FirstBXAfterTrain);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1);
   fChain->SetBranchAddress("HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg", &HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg, &b_HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg", &HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg, &b_HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg, &b_HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg, &b_HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg, &b_HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg, &b_HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg", &HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg, &b_HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg", &HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg, &b_HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg, &b_HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg, &b_HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg, &b_HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg, &b_HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr);
   fChain->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90, &b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90);
   fChain->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100, &b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100);
   fChain->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110, &b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110);
   fChain->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120, &b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120);
   fChain->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130, &b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130);
   fChain->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr, &b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr);
   fChain->SetBranchAddress("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr", &HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr, &b_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr);
   fChain->SetBranchAddress("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1", &HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1, &b_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", &HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL, &b_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT_Rsq0p35", &HLT_Rsq0p35, &b_HLT_Rsq0p35);
   fChain->SetBranchAddress("HLT_Rsq0p40", &HLT_Rsq0p40, &b_HLT_Rsq0p40);
   fChain->SetBranchAddress("HLT_RsqMR300_Rsq0p09_MR200", &HLT_RsqMR300_Rsq0p09_MR200, &b_HLT_RsqMR300_Rsq0p09_MR200);
   fChain->SetBranchAddress("HLT_RsqMR320_Rsq0p09_MR200", &HLT_RsqMR320_Rsq0p09_MR200, &b_HLT_RsqMR320_Rsq0p09_MR200);
   fChain->SetBranchAddress("HLT_RsqMR300_Rsq0p09_MR200_4jet", &HLT_RsqMR300_Rsq0p09_MR200_4jet, &b_HLT_RsqMR300_Rsq0p09_MR200_4jet);
   fChain->SetBranchAddress("HLT_RsqMR320_Rsq0p09_MR200_4jet", &HLT_RsqMR320_Rsq0p09_MR200_4jet, &b_HLT_RsqMR320_Rsq0p09_MR200_4jet);
   fChain->SetBranchAddress("HLT_L1_DoubleJet30_Mass_Min400_Mu10", &HLT_L1_DoubleJet30_Mass_Min400_Mu10, &b_HLT_L1_DoubleJet30_Mass_Min400_Mu10);
   fChain->SetBranchAddress("HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1", &HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1, &b_HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1", &HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1, &b_HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1", &HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1, &b_HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1);
   fChain->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50", &HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50, &b_HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50);
   fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3", &HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3, &b_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3);
   fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3", &HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3, &b_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3);
   fChain->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight_PFHT60", &HLT_PFMET100_PFMHT100_IDTight_PFHT60, &b_HLT_PFMET100_PFMHT100_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60", &HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60, &b_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60", &HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60, &b_HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_Mu18_Mu9_SameSign", &HLT_Mu18_Mu9_SameSign, &b_HLT_Mu18_Mu9_SameSign);
   fChain->SetBranchAddress("HLT_Mu18_Mu9_SameSign_DZ", &HLT_Mu18_Mu9_SameSign_DZ, &b_HLT_Mu18_Mu9_SameSign_DZ);
   fChain->SetBranchAddress("HLT_Mu18_Mu9", &HLT_Mu18_Mu9, &b_HLT_Mu18_Mu9);
   fChain->SetBranchAddress("HLT_Mu18_Mu9_DZ", &HLT_Mu18_Mu9_DZ, &b_HLT_Mu18_Mu9_DZ);
   fChain->SetBranchAddress("HLT_Mu20_Mu10_SameSign", &HLT_Mu20_Mu10_SameSign, &b_HLT_Mu20_Mu10_SameSign);
   fChain->SetBranchAddress("HLT_Mu20_Mu10_SameSign_DZ", &HLT_Mu20_Mu10_SameSign_DZ, &b_HLT_Mu20_Mu10_SameSign_DZ);
   fChain->SetBranchAddress("HLT_Mu20_Mu10", &HLT_Mu20_Mu10, &b_HLT_Mu20_Mu10);
   fChain->SetBranchAddress("HLT_Mu20_Mu10_DZ", &HLT_Mu20_Mu10_DZ, &b_HLT_Mu20_Mu10_DZ);
   fChain->SetBranchAddress("HLT_Mu23_Mu12_SameSign", &HLT_Mu23_Mu12_SameSign, &b_HLT_Mu23_Mu12_SameSign);
   fChain->SetBranchAddress("HLT_Mu23_Mu12_SameSign_DZ", &HLT_Mu23_Mu12_SameSign_DZ, &b_HLT_Mu23_Mu12_SameSign_DZ);
   fChain->SetBranchAddress("HLT_Mu23_Mu12", &HLT_Mu23_Mu12, &b_HLT_Mu23_Mu12);
   fChain->SetBranchAddress("HLT_Mu23_Mu12_DZ", &HLT_Mu23_Mu12_DZ, &b_HLT_Mu23_Mu12_DZ);
   fChain->SetBranchAddress("HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi", &HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi, &b_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi);
   fChain->SetBranchAddress("HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi", &HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi, &b_HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi);
   fChain->SetBranchAddress("HLT_DoubleMu3_DCA_PFMET50_PFMHT60", &HLT_DoubleMu3_DCA_PFMET50_PFMHT60, &b_HLT_DoubleMu3_DCA_PFMET50_PFMHT60);
   fChain->SetBranchAddress("HLT_TripleMu_5_3_3_Mass3p8to60_DCA", &HLT_TripleMu_5_3_3_Mass3p8to60_DCA, &b_HLT_TripleMu_5_3_3_Mass3p8to60_DCA);
   fChain->SetBranchAddress("HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1", &HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1, &b_HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1);
   fChain->SetBranchAddress("HLT_QuadPFJet103_88_75_15_DoubleBTagCSV_p013_p08_VBF1", &HLT_QuadPFJet103_88_75_15_DoubleBTagCSV_p013_p08_VBF1, &b_HLT_QuadPFJet103_88_75_15_DoubleBTagCSV_p013_p08_VBF1);
   fChain->SetBranchAddress("HLT_QuadPFJet105_90_76_15_DoubleBTagCSV_p013_p08_VBF1", &HLT_QuadPFJet105_90_76_15_DoubleBTagCSV_p013_p08_VBF1, &b_HLT_QuadPFJet105_90_76_15_DoubleBTagCSV_p013_p08_VBF1);
   fChain->SetBranchAddress("HLT_QuadPFJet111_90_80_15_DoubleBTagCSV_p013_p08_VBF1", &HLT_QuadPFJet111_90_80_15_DoubleBTagCSV_p013_p08_VBF1, &b_HLT_QuadPFJet111_90_80_15_DoubleBTagCSV_p013_p08_VBF1);
   fChain->SetBranchAddress("HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2", &HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2, &b_HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2);
   fChain->SetBranchAddress("HLT_QuadPFJet103_88_75_15_BTagCSV_p013_VBF2", &HLT_QuadPFJet103_88_75_15_BTagCSV_p013_VBF2, &b_HLT_QuadPFJet103_88_75_15_BTagCSV_p013_VBF2);
   fChain->SetBranchAddress("HLT_QuadPFJet105_88_76_15_BTagCSV_p013_VBF2", &HLT_QuadPFJet105_88_76_15_BTagCSV_p013_VBF2, &b_HLT_QuadPFJet105_88_76_15_BTagCSV_p013_VBF2);
   fChain->SetBranchAddress("HLT_QuadPFJet111_90_80_15_BTagCSV_p013_VBF2", &HLT_QuadPFJet111_90_80_15_BTagCSV_p013_VBF2, &b_HLT_QuadPFJet111_90_80_15_BTagCSV_p013_VBF2);
   fChain->SetBranchAddress("HLT_QuadPFJet98_83_71_15", &HLT_QuadPFJet98_83_71_15, &b_HLT_QuadPFJet98_83_71_15);
   fChain->SetBranchAddress("HLT_QuadPFJet103_88_75_15", &HLT_QuadPFJet103_88_75_15, &b_HLT_QuadPFJet103_88_75_15);
   fChain->SetBranchAddress("HLT_QuadPFJet105_88_76_15", &HLT_QuadPFJet105_88_76_15, &b_HLT_QuadPFJet105_88_76_15);
   fChain->SetBranchAddress("HLT_QuadPFJet111_90_80_15", &HLT_QuadPFJet111_90_80_15, &b_HLT_QuadPFJet111_90_80_15);
   fChain->SetBranchAddress("HLT_AK8PFJet330_PFAK8BTagCSV_p17", &HLT_AK8PFJet330_PFAK8BTagCSV_p17, &b_HLT_AK8PFJet330_PFAK8BTagCSV_p17);
   fChain->SetBranchAddress("HLT_AK8PFJet330_PFAK8BTagCSV_p1", &HLT_AK8PFJet330_PFAK8BTagCSV_p1, &b_HLT_AK8PFJet330_PFAK8BTagCSV_p1);
   fChain->SetBranchAddress("HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55", &HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55, &b_HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55);
   fChain->SetBranchAddress("HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55", &HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55, &b_HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55);
   fChain->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath, &b_HLTriggerFinalPath);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter", &Flag_CSCTightHaloTrkMuUnvetoFilter, &b_Flag_CSCTightHaloTrkMuUnvetoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, &b_Flag_CSCTightHalo2015Filter);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_HcalStripHaloFilter", &Flag_HcalStripHaloFilter, &b_Flag_HcalStripHaloFilter);
   fChain->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, &b_Flag_EcalDeadCellBoundaryEnergyFilter);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
   fChain->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters, &b_Flag_trkPOGFilters);
   fChain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter, &b_Flag_chargedHadronTrackResolutionFilter);
   fChain->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter, &b_Flag_muonBadTrackFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateSummer16Filter", &Flag_BadChargedCandidateSummer16Filter, &b_Flag_BadChargedCandidateSummer16Filter);
   fChain->SetBranchAddress("Flag_BadPFMuonSummer16Filter", &Flag_BadPFMuonSummer16Filter, &b_Flag_BadPFMuonSummer16Filter);
   fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, &b_Flag_trkPOG_manystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, &b_Flag_trkPOG_toomanystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, &b_Flag_trkPOG_logErrorTooManyClusters);
   fChain->SetBranchAddress("Flag_METFilters", &Flag_METFilters, &b_Flag_METFilters);
   
   fChain->SetBranchAddress("Jet_pt_raw", Jet_pt_raw, &b_Jet_pt_raw);
   fChain->SetBranchAddress("Jet_pt_nom", Jet_pt_nom, &b_Jet_pt_nom);
   fChain->SetBranchAddress("Jet_corr_JEC", Jet_corr_JEC, &b_Jet_corr_JEC);
   fChain->SetBranchAddress("Jet_corr_JER", Jet_corr_JER, &b_Jet_corr_JER);
   fChain->SetBranchAddress("Jet_mass_raw", Jet_mass_raw, &b_Jet_mass_raw);
   fChain->SetBranchAddress("Jet_mass_nom", Jet_mass_nom, &b_Jet_mass_nom);
   fChain->SetBranchAddress("MET_pt_nom", &MET_pt_nom, &b_MET_pt_nom);
   fChain->SetBranchAddress("MET_phi_nom", &MET_phi_nom, &b_MET_phi_nom);
   fChain->SetBranchAddress("Jet_pt_jerUp", Jet_pt_jerUp, &b_Jet_pt_jerUp);
   fChain->SetBranchAddress("Jet_mass_jerUp", Jet_mass_jerUp, &b_Jet_mass_jerUp);
   fChain->SetBranchAddress("Jet_mass_jmrUp", Jet_mass_jmrUp, &b_Jet_mass_jmrUp);
   fChain->SetBranchAddress("Jet_mass_jmsUp", Jet_mass_jmsUp, &b_Jet_mass_jmsUp);
   fChain->SetBranchAddress("MET_pt_jerUp", &MET_pt_jerUp, &b_MET_pt_jerUp);
   fChain->SetBranchAddress("MET_phi_jerUp", &MET_phi_jerUp, &b_MET_phi_jerUp);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteStatUp", Jet_pt_jesAbsoluteStatUp, &b_Jet_pt_jesAbsoluteStatUp);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteStatUp", Jet_mass_jesAbsoluteStatUp, &b_Jet_mass_jesAbsoluteStatUp);
   fChain->SetBranchAddress("MET_pt_jesAbsoluteStatUp", &MET_pt_jesAbsoluteStatUp, &b_MET_pt_jesAbsoluteStatUp);
   fChain->SetBranchAddress("MET_phi_jesAbsoluteStatUp", &MET_phi_jesAbsoluteStatUp, &b_MET_phi_jesAbsoluteStatUp);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteScaleUp", Jet_pt_jesAbsoluteScaleUp, &b_Jet_pt_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteScaleUp", Jet_mass_jesAbsoluteScaleUp, &b_Jet_mass_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("MET_pt_jesAbsoluteScaleUp", &MET_pt_jesAbsoluteScaleUp, &b_MET_pt_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("MET_phi_jesAbsoluteScaleUp", &MET_phi_jesAbsoluteScaleUp, &b_MET_phi_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteFlavMapUp", Jet_pt_jesAbsoluteFlavMapUp, &b_Jet_pt_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteFlavMapUp", Jet_mass_jesAbsoluteFlavMapUp, &b_Jet_mass_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("MET_pt_jesAbsoluteFlavMapUp", &MET_pt_jesAbsoluteFlavMapUp, &b_MET_pt_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("MET_phi_jesAbsoluteFlavMapUp", &MET_phi_jesAbsoluteFlavMapUp, &b_MET_phi_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteMPFBiasUp", Jet_pt_jesAbsoluteMPFBiasUp, &b_Jet_pt_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteMPFBiasUp", Jet_mass_jesAbsoluteMPFBiasUp, &b_Jet_mass_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("MET_pt_jesAbsoluteMPFBiasUp", &MET_pt_jesAbsoluteMPFBiasUp, &b_MET_pt_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("MET_phi_jesAbsoluteMPFBiasUp", &MET_phi_jesAbsoluteMPFBiasUp, &b_MET_phi_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("Jet_pt_jesFragmentationUp", Jet_pt_jesFragmentationUp, &b_Jet_pt_jesFragmentationUp);
   fChain->SetBranchAddress("Jet_mass_jesFragmentationUp", Jet_mass_jesFragmentationUp, &b_Jet_mass_jesFragmentationUp);
   fChain->SetBranchAddress("MET_pt_jesFragmentationUp", &MET_pt_jesFragmentationUp, &b_MET_pt_jesFragmentationUp);
   fChain->SetBranchAddress("MET_phi_jesFragmentationUp", &MET_phi_jesFragmentationUp, &b_MET_phi_jesFragmentationUp);
   fChain->SetBranchAddress("Jet_pt_jesSinglePionECALUp", Jet_pt_jesSinglePionECALUp, &b_Jet_pt_jesSinglePionECALUp);
   fChain->SetBranchAddress("Jet_mass_jesSinglePionECALUp", Jet_mass_jesSinglePionECALUp, &b_Jet_mass_jesSinglePionECALUp);
   fChain->SetBranchAddress("MET_pt_jesSinglePionECALUp", &MET_pt_jesSinglePionECALUp, &b_MET_pt_jesSinglePionECALUp);
   fChain->SetBranchAddress("MET_phi_jesSinglePionECALUp", &MET_phi_jesSinglePionECALUp, &b_MET_phi_jesSinglePionECALUp);
   fChain->SetBranchAddress("Jet_pt_jesSinglePionHCALUp", Jet_pt_jesSinglePionHCALUp, &b_Jet_pt_jesSinglePionHCALUp);
   fChain->SetBranchAddress("Jet_mass_jesSinglePionHCALUp", Jet_mass_jesSinglePionHCALUp, &b_Jet_mass_jesSinglePionHCALUp);
   fChain->SetBranchAddress("MET_pt_jesSinglePionHCALUp", &MET_pt_jesSinglePionHCALUp, &b_MET_pt_jesSinglePionHCALUp);
   fChain->SetBranchAddress("MET_phi_jesSinglePionHCALUp", &MET_phi_jesSinglePionHCALUp, &b_MET_phi_jesSinglePionHCALUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorQCDUp", Jet_pt_jesFlavorQCDUp, &b_Jet_pt_jesFlavorQCDUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorQCDUp", Jet_mass_jesFlavorQCDUp, &b_Jet_mass_jesFlavorQCDUp);
   fChain->SetBranchAddress("MET_pt_jesFlavorQCDUp", &MET_pt_jesFlavorQCDUp, &b_MET_pt_jesFlavorQCDUp);
   fChain->SetBranchAddress("MET_phi_jesFlavorQCDUp", &MET_phi_jesFlavorQCDUp, &b_MET_phi_jesFlavorQCDUp);
   fChain->SetBranchAddress("Jet_pt_jesTimePtEtaUp", Jet_pt_jesTimePtEtaUp, &b_Jet_pt_jesTimePtEtaUp);
   fChain->SetBranchAddress("Jet_mass_jesTimePtEtaUp", Jet_mass_jesTimePtEtaUp, &b_Jet_mass_jesTimePtEtaUp);
   fChain->SetBranchAddress("MET_pt_jesTimePtEtaUp", &MET_pt_jesTimePtEtaUp, &b_MET_pt_jesTimePtEtaUp);
   fChain->SetBranchAddress("MET_phi_jesTimePtEtaUp", &MET_phi_jesTimePtEtaUp, &b_MET_phi_jesTimePtEtaUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeJEREC1Up", Jet_pt_jesRelativeJEREC1Up, &b_Jet_pt_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("Jet_mass_jesRelativeJEREC1Up", Jet_mass_jesRelativeJEREC1Up, &b_Jet_mass_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("MET_pt_jesRelativeJEREC1Up", &MET_pt_jesRelativeJEREC1Up, &b_MET_pt_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("MET_phi_jesRelativeJEREC1Up", &MET_phi_jesRelativeJEREC1Up, &b_MET_phi_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("Jet_pt_jesRelativeJEREC2Up", Jet_pt_jesRelativeJEREC2Up, &b_Jet_pt_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("Jet_mass_jesRelativeJEREC2Up", Jet_mass_jesRelativeJEREC2Up, &b_Jet_mass_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("MET_pt_jesRelativeJEREC2Up", &MET_pt_jesRelativeJEREC2Up, &b_MET_pt_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("MET_phi_jesRelativeJEREC2Up", &MET_phi_jesRelativeJEREC2Up, &b_MET_phi_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("Jet_pt_jesRelativeJERHFUp", Jet_pt_jesRelativeJERHFUp, &b_Jet_pt_jesRelativeJERHFUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeJERHFUp", Jet_mass_jesRelativeJERHFUp, &b_Jet_mass_jesRelativeJERHFUp);
   fChain->SetBranchAddress("MET_pt_jesRelativeJERHFUp", &MET_pt_jesRelativeJERHFUp, &b_MET_pt_jesRelativeJERHFUp);
   fChain->SetBranchAddress("MET_phi_jesRelativeJERHFUp", &MET_phi_jesRelativeJERHFUp, &b_MET_phi_jesRelativeJERHFUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtBBUp", Jet_pt_jesRelativePtBBUp, &b_Jet_pt_jesRelativePtBBUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtBBUp", Jet_mass_jesRelativePtBBUp, &b_Jet_mass_jesRelativePtBBUp);
   fChain->SetBranchAddress("MET_pt_jesRelativePtBBUp", &MET_pt_jesRelativePtBBUp, &b_MET_pt_jesRelativePtBBUp);
   fChain->SetBranchAddress("MET_phi_jesRelativePtBBUp", &MET_phi_jesRelativePtBBUp, &b_MET_phi_jesRelativePtBBUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtEC1Up", Jet_pt_jesRelativePtEC1Up, &b_Jet_pt_jesRelativePtEC1Up);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtEC1Up", Jet_mass_jesRelativePtEC1Up, &b_Jet_mass_jesRelativePtEC1Up);
   fChain->SetBranchAddress("MET_pt_jesRelativePtEC1Up", &MET_pt_jesRelativePtEC1Up, &b_MET_pt_jesRelativePtEC1Up);
   fChain->SetBranchAddress("MET_phi_jesRelativePtEC1Up", &MET_phi_jesRelativePtEC1Up, &b_MET_phi_jesRelativePtEC1Up);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtEC2Up", Jet_pt_jesRelativePtEC2Up, &b_Jet_pt_jesRelativePtEC2Up);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtEC2Up", Jet_mass_jesRelativePtEC2Up, &b_Jet_mass_jesRelativePtEC2Up);
   fChain->SetBranchAddress("MET_pt_jesRelativePtEC2Up", &MET_pt_jesRelativePtEC2Up, &b_MET_pt_jesRelativePtEC2Up);
   fChain->SetBranchAddress("MET_phi_jesRelativePtEC2Up", &MET_phi_jesRelativePtEC2Up, &b_MET_phi_jesRelativePtEC2Up);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtHFUp", Jet_pt_jesRelativePtHFUp, &b_Jet_pt_jesRelativePtHFUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtHFUp", Jet_mass_jesRelativePtHFUp, &b_Jet_mass_jesRelativePtHFUp);
   fChain->SetBranchAddress("MET_pt_jesRelativePtHFUp", &MET_pt_jesRelativePtHFUp, &b_MET_pt_jesRelativePtHFUp);
   fChain->SetBranchAddress("MET_phi_jesRelativePtHFUp", &MET_phi_jesRelativePtHFUp, &b_MET_phi_jesRelativePtHFUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeBalUp", Jet_pt_jesRelativeBalUp, &b_Jet_pt_jesRelativeBalUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeBalUp", Jet_mass_jesRelativeBalUp, &b_Jet_mass_jesRelativeBalUp);
   fChain->SetBranchAddress("MET_pt_jesRelativeBalUp", &MET_pt_jesRelativeBalUp, &b_MET_pt_jesRelativeBalUp);
   fChain->SetBranchAddress("MET_phi_jesRelativeBalUp", &MET_phi_jesRelativeBalUp, &b_MET_phi_jesRelativeBalUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeSampleUp", Jet_pt_jesRelativeSampleUp, &b_Jet_pt_jesRelativeSampleUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeSampleUp", Jet_mass_jesRelativeSampleUp, &b_Jet_mass_jesRelativeSampleUp);
   fChain->SetBranchAddress("MET_pt_jesRelativeSampleUp", &MET_pt_jesRelativeSampleUp, &b_MET_pt_jesRelativeSampleUp);
   fChain->SetBranchAddress("MET_phi_jesRelativeSampleUp", &MET_phi_jesRelativeSampleUp, &b_MET_phi_jesRelativeSampleUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeFSRUp", Jet_pt_jesRelativeFSRUp, &b_Jet_pt_jesRelativeFSRUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeFSRUp", Jet_mass_jesRelativeFSRUp, &b_Jet_mass_jesRelativeFSRUp);
   fChain->SetBranchAddress("MET_pt_jesRelativeFSRUp", &MET_pt_jesRelativeFSRUp, &b_MET_pt_jesRelativeFSRUp);
   fChain->SetBranchAddress("MET_phi_jesRelativeFSRUp", &MET_phi_jesRelativeFSRUp, &b_MET_phi_jesRelativeFSRUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeStatFSRUp", Jet_pt_jesRelativeStatFSRUp, &b_Jet_pt_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeStatFSRUp", Jet_mass_jesRelativeStatFSRUp, &b_Jet_mass_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("MET_pt_jesRelativeStatFSRUp", &MET_pt_jesRelativeStatFSRUp, &b_MET_pt_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("MET_phi_jesRelativeStatFSRUp", &MET_phi_jesRelativeStatFSRUp, &b_MET_phi_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeStatECUp", Jet_pt_jesRelativeStatECUp, &b_Jet_pt_jesRelativeStatECUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeStatECUp", Jet_mass_jesRelativeStatECUp, &b_Jet_mass_jesRelativeStatECUp);
   fChain->SetBranchAddress("MET_pt_jesRelativeStatECUp", &MET_pt_jesRelativeStatECUp, &b_MET_pt_jesRelativeStatECUp);
   fChain->SetBranchAddress("MET_phi_jesRelativeStatECUp", &MET_phi_jesRelativeStatECUp, &b_MET_phi_jesRelativeStatECUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeStatHFUp", Jet_pt_jesRelativeStatHFUp, &b_Jet_pt_jesRelativeStatHFUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeStatHFUp", Jet_mass_jesRelativeStatHFUp, &b_Jet_mass_jesRelativeStatHFUp);
   fChain->SetBranchAddress("MET_pt_jesRelativeStatHFUp", &MET_pt_jesRelativeStatHFUp, &b_MET_pt_jesRelativeStatHFUp);
   fChain->SetBranchAddress("MET_phi_jesRelativeStatHFUp", &MET_phi_jesRelativeStatHFUp, &b_MET_phi_jesRelativeStatHFUp);
   fChain->SetBranchAddress("Jet_pt_jesPileUpDataMCUp", Jet_pt_jesPileUpDataMCUp, &b_Jet_pt_jesPileUpDataMCUp);
   fChain->SetBranchAddress("Jet_mass_jesPileUpDataMCUp", Jet_mass_jesPileUpDataMCUp, &b_Jet_mass_jesPileUpDataMCUp);
   fChain->SetBranchAddress("MET_pt_jesPileUpDataMCUp", &MET_pt_jesPileUpDataMCUp, &b_MET_pt_jesPileUpDataMCUp);
   fChain->SetBranchAddress("MET_phi_jesPileUpDataMCUp", &MET_phi_jesPileUpDataMCUp, &b_MET_phi_jesPileUpDataMCUp);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtRefUp", Jet_pt_jesPileUpPtRefUp, &b_Jet_pt_jesPileUpPtRefUp);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtRefUp", Jet_mass_jesPileUpPtRefUp, &b_Jet_mass_jesPileUpPtRefUp);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtRefUp", &MET_pt_jesPileUpPtRefUp, &b_MET_pt_jesPileUpPtRefUp);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtRefUp", &MET_phi_jesPileUpPtRefUp, &b_MET_phi_jesPileUpPtRefUp);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtBBUp", Jet_pt_jesPileUpPtBBUp, &b_Jet_pt_jesPileUpPtBBUp);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtBBUp", Jet_mass_jesPileUpPtBBUp, &b_Jet_mass_jesPileUpPtBBUp);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtBBUp", &MET_pt_jesPileUpPtBBUp, &b_MET_pt_jesPileUpPtBBUp);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtBBUp", &MET_phi_jesPileUpPtBBUp, &b_MET_phi_jesPileUpPtBBUp);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtEC1Up", Jet_pt_jesPileUpPtEC1Up, &b_Jet_pt_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtEC1Up", Jet_mass_jesPileUpPtEC1Up, &b_Jet_mass_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtEC1Up", &MET_pt_jesPileUpPtEC1Up, &b_MET_pt_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtEC1Up", &MET_phi_jesPileUpPtEC1Up, &b_MET_phi_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtEC2Up", Jet_pt_jesPileUpPtEC2Up, &b_Jet_pt_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtEC2Up", Jet_mass_jesPileUpPtEC2Up, &b_Jet_mass_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtEC2Up", &MET_pt_jesPileUpPtEC2Up, &b_MET_pt_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtEC2Up", &MET_phi_jesPileUpPtEC2Up, &b_MET_phi_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtHFUp", Jet_pt_jesPileUpPtHFUp, &b_Jet_pt_jesPileUpPtHFUp);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtHFUp", Jet_mass_jesPileUpPtHFUp, &b_Jet_mass_jesPileUpPtHFUp);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtHFUp", &MET_pt_jesPileUpPtHFUp, &b_MET_pt_jesPileUpPtHFUp);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtHFUp", &MET_phi_jesPileUpPtHFUp, &b_MET_phi_jesPileUpPtHFUp);
   fChain->SetBranchAddress("Jet_pt_jesPileUpMuZeroUp", Jet_pt_jesPileUpMuZeroUp, &b_Jet_pt_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("Jet_mass_jesPileUpMuZeroUp", Jet_mass_jesPileUpMuZeroUp, &b_Jet_mass_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("MET_pt_jesPileUpMuZeroUp", &MET_pt_jesPileUpMuZeroUp, &b_MET_pt_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("MET_phi_jesPileUpMuZeroUp", &MET_phi_jesPileUpMuZeroUp, &b_MET_phi_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("Jet_pt_jesPileUpEnvelopeUp", Jet_pt_jesPileUpEnvelopeUp, &b_Jet_pt_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("Jet_mass_jesPileUpEnvelopeUp", Jet_mass_jesPileUpEnvelopeUp, &b_Jet_mass_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("MET_pt_jesPileUpEnvelopeUp", &MET_pt_jesPileUpEnvelopeUp, &b_MET_pt_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("MET_phi_jesPileUpEnvelopeUp", &MET_phi_jesPileUpEnvelopeUp, &b_MET_phi_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalPileUpUp", Jet_pt_jesSubTotalPileUpUp, &b_Jet_pt_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalPileUpUp", Jet_mass_jesSubTotalPileUpUp, &b_Jet_mass_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("MET_pt_jesSubTotalPileUpUp", &MET_pt_jesSubTotalPileUpUp, &b_MET_pt_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("MET_phi_jesSubTotalPileUpUp", &MET_phi_jesSubTotalPileUpUp, &b_MET_phi_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalRelativeUp", Jet_pt_jesSubTotalRelativeUp, &b_Jet_pt_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalRelativeUp", Jet_mass_jesSubTotalRelativeUp, &b_Jet_mass_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("MET_pt_jesSubTotalRelativeUp", &MET_pt_jesSubTotalRelativeUp, &b_MET_pt_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("MET_phi_jesSubTotalRelativeUp", &MET_phi_jesSubTotalRelativeUp, &b_MET_phi_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalPtUp", Jet_pt_jesSubTotalPtUp, &b_Jet_pt_jesSubTotalPtUp);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalPtUp", Jet_mass_jesSubTotalPtUp, &b_Jet_mass_jesSubTotalPtUp);
   fChain->SetBranchAddress("MET_pt_jesSubTotalPtUp", &MET_pt_jesSubTotalPtUp, &b_MET_pt_jesSubTotalPtUp);
   fChain->SetBranchAddress("MET_phi_jesSubTotalPtUp", &MET_phi_jesSubTotalPtUp, &b_MET_phi_jesSubTotalPtUp);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalScaleUp", Jet_pt_jesSubTotalScaleUp, &b_Jet_pt_jesSubTotalScaleUp);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalScaleUp", Jet_mass_jesSubTotalScaleUp, &b_Jet_mass_jesSubTotalScaleUp);
   fChain->SetBranchAddress("MET_pt_jesSubTotalScaleUp", &MET_pt_jesSubTotalScaleUp, &b_MET_pt_jesSubTotalScaleUp);
   fChain->SetBranchAddress("MET_phi_jesSubTotalScaleUp", &MET_phi_jesSubTotalScaleUp, &b_MET_phi_jesSubTotalScaleUp);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalAbsoluteUp", Jet_pt_jesSubTotalAbsoluteUp, &b_Jet_pt_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalAbsoluteUp", Jet_mass_jesSubTotalAbsoluteUp, &b_Jet_mass_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("MET_pt_jesSubTotalAbsoluteUp", &MET_pt_jesSubTotalAbsoluteUp, &b_MET_pt_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("MET_phi_jesSubTotalAbsoluteUp", &MET_phi_jesSubTotalAbsoluteUp, &b_MET_phi_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalMCUp", Jet_pt_jesSubTotalMCUp, &b_Jet_pt_jesSubTotalMCUp);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalMCUp", Jet_mass_jesSubTotalMCUp, &b_Jet_mass_jesSubTotalMCUp);
   fChain->SetBranchAddress("MET_pt_jesSubTotalMCUp", &MET_pt_jesSubTotalMCUp, &b_MET_pt_jesSubTotalMCUp);
   fChain->SetBranchAddress("MET_phi_jesSubTotalMCUp", &MET_phi_jesSubTotalMCUp, &b_MET_phi_jesSubTotalMCUp);
   fChain->SetBranchAddress("Jet_pt_jesTotalUp", Jet_pt_jesTotalUp, &b_Jet_pt_jesTotalUp);
   fChain->SetBranchAddress("Jet_mass_jesTotalUp", Jet_mass_jesTotalUp, &b_Jet_mass_jesTotalUp);
   fChain->SetBranchAddress("MET_pt_jesTotalUp", &MET_pt_jesTotalUp, &b_MET_pt_jesTotalUp);
   fChain->SetBranchAddress("MET_phi_jesTotalUp", &MET_phi_jesTotalUp, &b_MET_phi_jesTotalUp);
   fChain->SetBranchAddress("Jet_pt_jesTotalNoFlavorUp", Jet_pt_jesTotalNoFlavorUp, &b_Jet_pt_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("Jet_mass_jesTotalNoFlavorUp", Jet_mass_jesTotalNoFlavorUp, &b_Jet_mass_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("MET_pt_jesTotalNoFlavorUp", &MET_pt_jesTotalNoFlavorUp, &b_MET_pt_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("MET_phi_jesTotalNoFlavorUp", &MET_phi_jesTotalNoFlavorUp, &b_MET_phi_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("Jet_pt_jesTotalNoTimeUp", Jet_pt_jesTotalNoTimeUp, &b_Jet_pt_jesTotalNoTimeUp);
   fChain->SetBranchAddress("Jet_mass_jesTotalNoTimeUp", Jet_mass_jesTotalNoTimeUp, &b_Jet_mass_jesTotalNoTimeUp);
   fChain->SetBranchAddress("MET_pt_jesTotalNoTimeUp", &MET_pt_jesTotalNoTimeUp, &b_MET_pt_jesTotalNoTimeUp);
   fChain->SetBranchAddress("MET_phi_jesTotalNoTimeUp", &MET_phi_jesTotalNoTimeUp, &b_MET_phi_jesTotalNoTimeUp);
   fChain->SetBranchAddress("Jet_pt_jesTotalNoFlavorNoTimeUp", Jet_pt_jesTotalNoFlavorNoTimeUp, &b_Jet_pt_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("Jet_mass_jesTotalNoFlavorNoTimeUp", Jet_mass_jesTotalNoFlavorNoTimeUp, &b_Jet_mass_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("MET_pt_jesTotalNoFlavorNoTimeUp", &MET_pt_jesTotalNoFlavorNoTimeUp, &b_MET_pt_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("MET_phi_jesTotalNoFlavorNoTimeUp", &MET_phi_jesTotalNoFlavorNoTimeUp, &b_MET_phi_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorZJetUp", Jet_pt_jesFlavorZJetUp, &b_Jet_pt_jesFlavorZJetUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorZJetUp", Jet_mass_jesFlavorZJetUp, &b_Jet_mass_jesFlavorZJetUp);
   fChain->SetBranchAddress("MET_pt_jesFlavorZJetUp", &MET_pt_jesFlavorZJetUp, &b_MET_pt_jesFlavorZJetUp);
   fChain->SetBranchAddress("MET_phi_jesFlavorZJetUp", &MET_phi_jesFlavorZJetUp, &b_MET_phi_jesFlavorZJetUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPhotonJetUp", Jet_pt_jesFlavorPhotonJetUp, &b_Jet_pt_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPhotonJetUp", Jet_mass_jesFlavorPhotonJetUp, &b_Jet_mass_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("MET_pt_jesFlavorPhotonJetUp", &MET_pt_jesFlavorPhotonJetUp, &b_MET_pt_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("MET_phi_jesFlavorPhotonJetUp", &MET_phi_jesFlavorPhotonJetUp, &b_MET_phi_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureGluonUp", Jet_pt_jesFlavorPureGluonUp, &b_Jet_pt_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureGluonUp", Jet_mass_jesFlavorPureGluonUp, &b_Jet_mass_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("MET_pt_jesFlavorPureGluonUp", &MET_pt_jesFlavorPureGluonUp, &b_MET_pt_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("MET_phi_jesFlavorPureGluonUp", &MET_phi_jesFlavorPureGluonUp, &b_MET_phi_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureQuarkUp", Jet_pt_jesFlavorPureQuarkUp, &b_Jet_pt_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureQuarkUp", Jet_mass_jesFlavorPureQuarkUp, &b_Jet_mass_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("MET_pt_jesFlavorPureQuarkUp", &MET_pt_jesFlavorPureQuarkUp, &b_MET_pt_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("MET_phi_jesFlavorPureQuarkUp", &MET_phi_jesFlavorPureQuarkUp, &b_MET_phi_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureCharmUp", Jet_pt_jesFlavorPureCharmUp, &b_Jet_pt_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureCharmUp", Jet_mass_jesFlavorPureCharmUp, &b_Jet_mass_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("MET_pt_jesFlavorPureCharmUp", &MET_pt_jesFlavorPureCharmUp, &b_MET_pt_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("MET_phi_jesFlavorPureCharmUp", &MET_phi_jesFlavorPureCharmUp, &b_MET_phi_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureBottomUp", Jet_pt_jesFlavorPureBottomUp, &b_Jet_pt_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureBottomUp", Jet_mass_jesFlavorPureBottomUp, &b_Jet_mass_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("MET_pt_jesFlavorPureBottomUp", &MET_pt_jesFlavorPureBottomUp, &b_MET_pt_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("MET_phi_jesFlavorPureBottomUp", &MET_phi_jesFlavorPureBottomUp, &b_MET_phi_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunBUp", Jet_pt_jesTimeRunBUp, &b_Jet_pt_jesTimeRunBUp);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunBUp", Jet_mass_jesTimeRunBUp, &b_Jet_mass_jesTimeRunBUp);
   fChain->SetBranchAddress("MET_pt_jesTimeRunBUp", &MET_pt_jesTimeRunBUp, &b_MET_pt_jesTimeRunBUp);
   fChain->SetBranchAddress("MET_phi_jesTimeRunBUp", &MET_phi_jesTimeRunBUp, &b_MET_phi_jesTimeRunBUp);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunCUp", Jet_pt_jesTimeRunCUp, &b_Jet_pt_jesTimeRunCUp);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunCUp", Jet_mass_jesTimeRunCUp, &b_Jet_mass_jesTimeRunCUp);
   fChain->SetBranchAddress("MET_pt_jesTimeRunCUp", &MET_pt_jesTimeRunCUp, &b_MET_pt_jesTimeRunCUp);
   fChain->SetBranchAddress("MET_phi_jesTimeRunCUp", &MET_phi_jesTimeRunCUp, &b_MET_phi_jesTimeRunCUp);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunDEUp", Jet_pt_jesTimeRunDEUp, &b_Jet_pt_jesTimeRunDEUp);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunDEUp", Jet_mass_jesTimeRunDEUp, &b_Jet_mass_jesTimeRunDEUp);
   fChain->SetBranchAddress("MET_pt_jesTimeRunDEUp", &MET_pt_jesTimeRunDEUp, &b_MET_pt_jesTimeRunDEUp);
   fChain->SetBranchAddress("MET_phi_jesTimeRunDEUp", &MET_phi_jesTimeRunDEUp, &b_MET_phi_jesTimeRunDEUp);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunFUp", Jet_pt_jesTimeRunFUp, &b_Jet_pt_jesTimeRunFUp);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunFUp", Jet_mass_jesTimeRunFUp, &b_Jet_mass_jesTimeRunFUp);
   fChain->SetBranchAddress("MET_pt_jesTimeRunFUp", &MET_pt_jesTimeRunFUp, &b_MET_pt_jesTimeRunFUp);
   fChain->SetBranchAddress("MET_phi_jesTimeRunFUp", &MET_phi_jesTimeRunFUp, &b_MET_phi_jesTimeRunFUp);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupMPFInSituUp", Jet_pt_jesCorrelationGroupMPFInSituUp, &b_Jet_pt_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupMPFInSituUp", Jet_mass_jesCorrelationGroupMPFInSituUp, &b_Jet_mass_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupMPFInSituUp", &MET_pt_jesCorrelationGroupMPFInSituUp, &b_MET_pt_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupMPFInSituUp", &MET_phi_jesCorrelationGroupMPFInSituUp, &b_MET_phi_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupIntercalibrationUp", Jet_pt_jesCorrelationGroupIntercalibrationUp, &b_Jet_pt_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupIntercalibrationUp", Jet_mass_jesCorrelationGroupIntercalibrationUp, &b_Jet_mass_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupIntercalibrationUp", &MET_pt_jesCorrelationGroupIntercalibrationUp, &b_MET_pt_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupIntercalibrationUp", &MET_phi_jesCorrelationGroupIntercalibrationUp, &b_MET_phi_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupbJESUp", Jet_pt_jesCorrelationGroupbJESUp, &b_Jet_pt_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupbJESUp", Jet_mass_jesCorrelationGroupbJESUp, &b_Jet_mass_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupbJESUp", &MET_pt_jesCorrelationGroupbJESUp, &b_MET_pt_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupbJESUp", &MET_phi_jesCorrelationGroupbJESUp, &b_MET_phi_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupFlavorUp", Jet_pt_jesCorrelationGroupFlavorUp, &b_Jet_pt_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupFlavorUp", Jet_mass_jesCorrelationGroupFlavorUp, &b_Jet_mass_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupFlavorUp", &MET_pt_jesCorrelationGroupFlavorUp, &b_MET_pt_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupFlavorUp", &MET_phi_jesCorrelationGroupFlavorUp, &b_MET_phi_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupUncorrelatedUp", Jet_pt_jesCorrelationGroupUncorrelatedUp, &b_Jet_pt_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupUncorrelatedUp", Jet_mass_jesCorrelationGroupUncorrelatedUp, &b_Jet_mass_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupUncorrelatedUp", &MET_pt_jesCorrelationGroupUncorrelatedUp, &b_MET_pt_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupUncorrelatedUp", &MET_phi_jesCorrelationGroupUncorrelatedUp, &b_MET_phi_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("MET_pt_unclustEnUp", &MET_pt_unclustEnUp, &b_MET_pt_unclustEnUp);
   fChain->SetBranchAddress("MET_phi_unclustEnUp", &MET_phi_unclustEnUp, &b_MET_phi_unclustEnUp);
   fChain->SetBranchAddress("Jet_pt_jerDown", Jet_pt_jerDown, &b_Jet_pt_jerDown);
   fChain->SetBranchAddress("Jet_mass_jerDown", Jet_mass_jerDown, &b_Jet_mass_jerDown);
   fChain->SetBranchAddress("Jet_mass_jmrDown", Jet_mass_jmrDown, &b_Jet_mass_jmrDown);
   fChain->SetBranchAddress("Jet_mass_jmsDown", Jet_mass_jmsDown, &b_Jet_mass_jmsDown);
   fChain->SetBranchAddress("MET_pt_jerDown", &MET_pt_jerDown, &b_MET_pt_jerDown);
   fChain->SetBranchAddress("MET_phi_jerDown", &MET_phi_jerDown, &b_MET_phi_jerDown);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteStatDown", Jet_pt_jesAbsoluteStatDown, &b_Jet_pt_jesAbsoluteStatDown);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteStatDown", Jet_mass_jesAbsoluteStatDown, &b_Jet_mass_jesAbsoluteStatDown);
   fChain->SetBranchAddress("MET_pt_jesAbsoluteStatDown", &MET_pt_jesAbsoluteStatDown, &b_MET_pt_jesAbsoluteStatDown);
   fChain->SetBranchAddress("MET_phi_jesAbsoluteStatDown", &MET_phi_jesAbsoluteStatDown, &b_MET_phi_jesAbsoluteStatDown);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteScaleDown", Jet_pt_jesAbsoluteScaleDown, &b_Jet_pt_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteScaleDown", Jet_mass_jesAbsoluteScaleDown, &b_Jet_mass_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("MET_pt_jesAbsoluteScaleDown", &MET_pt_jesAbsoluteScaleDown, &b_MET_pt_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("MET_phi_jesAbsoluteScaleDown", &MET_phi_jesAbsoluteScaleDown, &b_MET_phi_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteFlavMapDown", Jet_pt_jesAbsoluteFlavMapDown, &b_Jet_pt_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteFlavMapDown", Jet_mass_jesAbsoluteFlavMapDown, &b_Jet_mass_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("MET_pt_jesAbsoluteFlavMapDown", &MET_pt_jesAbsoluteFlavMapDown, &b_MET_pt_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("MET_phi_jesAbsoluteFlavMapDown", &MET_phi_jesAbsoluteFlavMapDown, &b_MET_phi_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteMPFBiasDown", Jet_pt_jesAbsoluteMPFBiasDown, &b_Jet_pt_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteMPFBiasDown", Jet_mass_jesAbsoluteMPFBiasDown, &b_Jet_mass_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("MET_pt_jesAbsoluteMPFBiasDown", &MET_pt_jesAbsoluteMPFBiasDown, &b_MET_pt_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("MET_phi_jesAbsoluteMPFBiasDown", &MET_phi_jesAbsoluteMPFBiasDown, &b_MET_phi_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("Jet_pt_jesFragmentationDown", Jet_pt_jesFragmentationDown, &b_Jet_pt_jesFragmentationDown);
   fChain->SetBranchAddress("Jet_mass_jesFragmentationDown", Jet_mass_jesFragmentationDown, &b_Jet_mass_jesFragmentationDown);
   fChain->SetBranchAddress("MET_pt_jesFragmentationDown", &MET_pt_jesFragmentationDown, &b_MET_pt_jesFragmentationDown);
   fChain->SetBranchAddress("MET_phi_jesFragmentationDown", &MET_phi_jesFragmentationDown, &b_MET_phi_jesFragmentationDown);
   fChain->SetBranchAddress("Jet_pt_jesSinglePionECALDown", Jet_pt_jesSinglePionECALDown, &b_Jet_pt_jesSinglePionECALDown);
   fChain->SetBranchAddress("Jet_mass_jesSinglePionECALDown", Jet_mass_jesSinglePionECALDown, &b_Jet_mass_jesSinglePionECALDown);
   fChain->SetBranchAddress("MET_pt_jesSinglePionECALDown", &MET_pt_jesSinglePionECALDown, &b_MET_pt_jesSinglePionECALDown);
   fChain->SetBranchAddress("MET_phi_jesSinglePionECALDown", &MET_phi_jesSinglePionECALDown, &b_MET_phi_jesSinglePionECALDown);
   fChain->SetBranchAddress("Jet_pt_jesSinglePionHCALDown", Jet_pt_jesSinglePionHCALDown, &b_Jet_pt_jesSinglePionHCALDown);
   fChain->SetBranchAddress("Jet_mass_jesSinglePionHCALDown", Jet_mass_jesSinglePionHCALDown, &b_Jet_mass_jesSinglePionHCALDown);
   fChain->SetBranchAddress("MET_pt_jesSinglePionHCALDown", &MET_pt_jesSinglePionHCALDown, &b_MET_pt_jesSinglePionHCALDown);
   fChain->SetBranchAddress("MET_phi_jesSinglePionHCALDown", &MET_phi_jesSinglePionHCALDown, &b_MET_phi_jesSinglePionHCALDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorQCDDown", Jet_pt_jesFlavorQCDDown, &b_Jet_pt_jesFlavorQCDDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorQCDDown", Jet_mass_jesFlavorQCDDown, &b_Jet_mass_jesFlavorQCDDown);
   fChain->SetBranchAddress("MET_pt_jesFlavorQCDDown", &MET_pt_jesFlavorQCDDown, &b_MET_pt_jesFlavorQCDDown);
   fChain->SetBranchAddress("MET_phi_jesFlavorQCDDown", &MET_phi_jesFlavorQCDDown, &b_MET_phi_jesFlavorQCDDown);
   fChain->SetBranchAddress("Jet_pt_jesTimePtEtaDown", Jet_pt_jesTimePtEtaDown, &b_Jet_pt_jesTimePtEtaDown);
   fChain->SetBranchAddress("Jet_mass_jesTimePtEtaDown", Jet_mass_jesTimePtEtaDown, &b_Jet_mass_jesTimePtEtaDown);
   fChain->SetBranchAddress("MET_pt_jesTimePtEtaDown", &MET_pt_jesTimePtEtaDown, &b_MET_pt_jesTimePtEtaDown);
   fChain->SetBranchAddress("MET_phi_jesTimePtEtaDown", &MET_phi_jesTimePtEtaDown, &b_MET_phi_jesTimePtEtaDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeJEREC1Down", Jet_pt_jesRelativeJEREC1Down, &b_Jet_pt_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("Jet_mass_jesRelativeJEREC1Down", Jet_mass_jesRelativeJEREC1Down, &b_Jet_mass_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("MET_pt_jesRelativeJEREC1Down", &MET_pt_jesRelativeJEREC1Down, &b_MET_pt_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("MET_phi_jesRelativeJEREC1Down", &MET_phi_jesRelativeJEREC1Down, &b_MET_phi_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("Jet_pt_jesRelativeJEREC2Down", Jet_pt_jesRelativeJEREC2Down, &b_Jet_pt_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("Jet_mass_jesRelativeJEREC2Down", Jet_mass_jesRelativeJEREC2Down, &b_Jet_mass_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("MET_pt_jesRelativeJEREC2Down", &MET_pt_jesRelativeJEREC2Down, &b_MET_pt_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("MET_phi_jesRelativeJEREC2Down", &MET_phi_jesRelativeJEREC2Down, &b_MET_phi_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("Jet_pt_jesRelativeJERHFDown", Jet_pt_jesRelativeJERHFDown, &b_Jet_pt_jesRelativeJERHFDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeJERHFDown", Jet_mass_jesRelativeJERHFDown, &b_Jet_mass_jesRelativeJERHFDown);
   fChain->SetBranchAddress("MET_pt_jesRelativeJERHFDown", &MET_pt_jesRelativeJERHFDown, &b_MET_pt_jesRelativeJERHFDown);
   fChain->SetBranchAddress("MET_phi_jesRelativeJERHFDown", &MET_phi_jesRelativeJERHFDown, &b_MET_phi_jesRelativeJERHFDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtBBDown", Jet_pt_jesRelativePtBBDown, &b_Jet_pt_jesRelativePtBBDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtBBDown", Jet_mass_jesRelativePtBBDown, &b_Jet_mass_jesRelativePtBBDown);
   fChain->SetBranchAddress("MET_pt_jesRelativePtBBDown", &MET_pt_jesRelativePtBBDown, &b_MET_pt_jesRelativePtBBDown);
   fChain->SetBranchAddress("MET_phi_jesRelativePtBBDown", &MET_phi_jesRelativePtBBDown, &b_MET_phi_jesRelativePtBBDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtEC1Down", Jet_pt_jesRelativePtEC1Down, &b_Jet_pt_jesRelativePtEC1Down);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtEC1Down", Jet_mass_jesRelativePtEC1Down, &b_Jet_mass_jesRelativePtEC1Down);
   fChain->SetBranchAddress("MET_pt_jesRelativePtEC1Down", &MET_pt_jesRelativePtEC1Down, &b_MET_pt_jesRelativePtEC1Down);
   fChain->SetBranchAddress("MET_phi_jesRelativePtEC1Down", &MET_phi_jesRelativePtEC1Down, &b_MET_phi_jesRelativePtEC1Down);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtEC2Down", Jet_pt_jesRelativePtEC2Down, &b_Jet_pt_jesRelativePtEC2Down);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtEC2Down", Jet_mass_jesRelativePtEC2Down, &b_Jet_mass_jesRelativePtEC2Down);
   fChain->SetBranchAddress("MET_pt_jesRelativePtEC2Down", &MET_pt_jesRelativePtEC2Down, &b_MET_pt_jesRelativePtEC2Down);
   fChain->SetBranchAddress("MET_phi_jesRelativePtEC2Down", &MET_phi_jesRelativePtEC2Down, &b_MET_phi_jesRelativePtEC2Down);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtHFDown", Jet_pt_jesRelativePtHFDown, &b_Jet_pt_jesRelativePtHFDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtHFDown", Jet_mass_jesRelativePtHFDown, &b_Jet_mass_jesRelativePtHFDown);
   fChain->SetBranchAddress("MET_pt_jesRelativePtHFDown", &MET_pt_jesRelativePtHFDown, &b_MET_pt_jesRelativePtHFDown);
   fChain->SetBranchAddress("MET_phi_jesRelativePtHFDown", &MET_phi_jesRelativePtHFDown, &b_MET_phi_jesRelativePtHFDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeBalDown", Jet_pt_jesRelativeBalDown, &b_Jet_pt_jesRelativeBalDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeBalDown", Jet_mass_jesRelativeBalDown, &b_Jet_mass_jesRelativeBalDown);
   fChain->SetBranchAddress("MET_pt_jesRelativeBalDown", &MET_pt_jesRelativeBalDown, &b_MET_pt_jesRelativeBalDown);
   fChain->SetBranchAddress("MET_phi_jesRelativeBalDown", &MET_phi_jesRelativeBalDown, &b_MET_phi_jesRelativeBalDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeSampleDown", Jet_pt_jesRelativeSampleDown, &b_Jet_pt_jesRelativeSampleDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeSampleDown", Jet_mass_jesRelativeSampleDown, &b_Jet_mass_jesRelativeSampleDown);
   fChain->SetBranchAddress("MET_pt_jesRelativeSampleDown", &MET_pt_jesRelativeSampleDown, &b_MET_pt_jesRelativeSampleDown);
   fChain->SetBranchAddress("MET_phi_jesRelativeSampleDown", &MET_phi_jesRelativeSampleDown, &b_MET_phi_jesRelativeSampleDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeFSRDown", Jet_pt_jesRelativeFSRDown, &b_Jet_pt_jesRelativeFSRDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeFSRDown", Jet_mass_jesRelativeFSRDown, &b_Jet_mass_jesRelativeFSRDown);
   fChain->SetBranchAddress("MET_pt_jesRelativeFSRDown", &MET_pt_jesRelativeFSRDown, &b_MET_pt_jesRelativeFSRDown);
   fChain->SetBranchAddress("MET_phi_jesRelativeFSRDown", &MET_phi_jesRelativeFSRDown, &b_MET_phi_jesRelativeFSRDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeStatFSRDown", Jet_pt_jesRelativeStatFSRDown, &b_Jet_pt_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeStatFSRDown", Jet_mass_jesRelativeStatFSRDown, &b_Jet_mass_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("MET_pt_jesRelativeStatFSRDown", &MET_pt_jesRelativeStatFSRDown, &b_MET_pt_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("MET_phi_jesRelativeStatFSRDown", &MET_phi_jesRelativeStatFSRDown, &b_MET_phi_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeStatECDown", Jet_pt_jesRelativeStatECDown, &b_Jet_pt_jesRelativeStatECDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeStatECDown", Jet_mass_jesRelativeStatECDown, &b_Jet_mass_jesRelativeStatECDown);
   fChain->SetBranchAddress("MET_pt_jesRelativeStatECDown", &MET_pt_jesRelativeStatECDown, &b_MET_pt_jesRelativeStatECDown);
   fChain->SetBranchAddress("MET_phi_jesRelativeStatECDown", &MET_phi_jesRelativeStatECDown, &b_MET_phi_jesRelativeStatECDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeStatHFDown", Jet_pt_jesRelativeStatHFDown, &b_Jet_pt_jesRelativeStatHFDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeStatHFDown", Jet_mass_jesRelativeStatHFDown, &b_Jet_mass_jesRelativeStatHFDown);
   fChain->SetBranchAddress("MET_pt_jesRelativeStatHFDown", &MET_pt_jesRelativeStatHFDown, &b_MET_pt_jesRelativeStatHFDown);
   fChain->SetBranchAddress("MET_phi_jesRelativeStatHFDown", &MET_phi_jesRelativeStatHFDown, &b_MET_phi_jesRelativeStatHFDown);
   fChain->SetBranchAddress("Jet_pt_jesPileUpDataMCDown", Jet_pt_jesPileUpDataMCDown, &b_Jet_pt_jesPileUpDataMCDown);
   fChain->SetBranchAddress("Jet_mass_jesPileUpDataMCDown", Jet_mass_jesPileUpDataMCDown, &b_Jet_mass_jesPileUpDataMCDown);
   fChain->SetBranchAddress("MET_pt_jesPileUpDataMCDown", &MET_pt_jesPileUpDataMCDown, &b_MET_pt_jesPileUpDataMCDown);
   fChain->SetBranchAddress("MET_phi_jesPileUpDataMCDown", &MET_phi_jesPileUpDataMCDown, &b_MET_phi_jesPileUpDataMCDown);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtRefDown", Jet_pt_jesPileUpPtRefDown, &b_Jet_pt_jesPileUpPtRefDown);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtRefDown", Jet_mass_jesPileUpPtRefDown, &b_Jet_mass_jesPileUpPtRefDown);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtRefDown", &MET_pt_jesPileUpPtRefDown, &b_MET_pt_jesPileUpPtRefDown);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtRefDown", &MET_phi_jesPileUpPtRefDown, &b_MET_phi_jesPileUpPtRefDown);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtBBDown", Jet_pt_jesPileUpPtBBDown, &b_Jet_pt_jesPileUpPtBBDown);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtBBDown", Jet_mass_jesPileUpPtBBDown, &b_Jet_mass_jesPileUpPtBBDown);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtBBDown", &MET_pt_jesPileUpPtBBDown, &b_MET_pt_jesPileUpPtBBDown);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtBBDown", &MET_phi_jesPileUpPtBBDown, &b_MET_phi_jesPileUpPtBBDown);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtEC1Down", Jet_pt_jesPileUpPtEC1Down, &b_Jet_pt_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtEC1Down", Jet_mass_jesPileUpPtEC1Down, &b_Jet_mass_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtEC1Down", &MET_pt_jesPileUpPtEC1Down, &b_MET_pt_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtEC1Down", &MET_phi_jesPileUpPtEC1Down, &b_MET_phi_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtEC2Down", Jet_pt_jesPileUpPtEC2Down, &b_Jet_pt_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtEC2Down", Jet_mass_jesPileUpPtEC2Down, &b_Jet_mass_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtEC2Down", &MET_pt_jesPileUpPtEC2Down, &b_MET_pt_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtEC2Down", &MET_phi_jesPileUpPtEC2Down, &b_MET_phi_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtHFDown", Jet_pt_jesPileUpPtHFDown, &b_Jet_pt_jesPileUpPtHFDown);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtHFDown", Jet_mass_jesPileUpPtHFDown, &b_Jet_mass_jesPileUpPtHFDown);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtHFDown", &MET_pt_jesPileUpPtHFDown, &b_MET_pt_jesPileUpPtHFDown);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtHFDown", &MET_phi_jesPileUpPtHFDown, &b_MET_phi_jesPileUpPtHFDown);
   fChain->SetBranchAddress("Jet_pt_jesPileUpMuZeroDown", Jet_pt_jesPileUpMuZeroDown, &b_Jet_pt_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("Jet_mass_jesPileUpMuZeroDown", Jet_mass_jesPileUpMuZeroDown, &b_Jet_mass_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("MET_pt_jesPileUpMuZeroDown", &MET_pt_jesPileUpMuZeroDown, &b_MET_pt_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("MET_phi_jesPileUpMuZeroDown", &MET_phi_jesPileUpMuZeroDown, &b_MET_phi_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("Jet_pt_jesPileUpEnvelopeDown", Jet_pt_jesPileUpEnvelopeDown, &b_Jet_pt_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("Jet_mass_jesPileUpEnvelopeDown", Jet_mass_jesPileUpEnvelopeDown, &b_Jet_mass_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("MET_pt_jesPileUpEnvelopeDown", &MET_pt_jesPileUpEnvelopeDown, &b_MET_pt_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("MET_phi_jesPileUpEnvelopeDown", &MET_phi_jesPileUpEnvelopeDown, &b_MET_phi_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalPileUpDown", Jet_pt_jesSubTotalPileUpDown, &b_Jet_pt_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalPileUpDown", Jet_mass_jesSubTotalPileUpDown, &b_Jet_mass_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("MET_pt_jesSubTotalPileUpDown", &MET_pt_jesSubTotalPileUpDown, &b_MET_pt_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("MET_phi_jesSubTotalPileUpDown", &MET_phi_jesSubTotalPileUpDown, &b_MET_phi_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalRelativeDown", Jet_pt_jesSubTotalRelativeDown, &b_Jet_pt_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalRelativeDown", Jet_mass_jesSubTotalRelativeDown, &b_Jet_mass_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("MET_pt_jesSubTotalRelativeDown", &MET_pt_jesSubTotalRelativeDown, &b_MET_pt_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("MET_phi_jesSubTotalRelativeDown", &MET_phi_jesSubTotalRelativeDown, &b_MET_phi_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalPtDown", Jet_pt_jesSubTotalPtDown, &b_Jet_pt_jesSubTotalPtDown);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalPtDown", Jet_mass_jesSubTotalPtDown, &b_Jet_mass_jesSubTotalPtDown);
   fChain->SetBranchAddress("MET_pt_jesSubTotalPtDown", &MET_pt_jesSubTotalPtDown, &b_MET_pt_jesSubTotalPtDown);
   fChain->SetBranchAddress("MET_phi_jesSubTotalPtDown", &MET_phi_jesSubTotalPtDown, &b_MET_phi_jesSubTotalPtDown);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalScaleDown", Jet_pt_jesSubTotalScaleDown, &b_Jet_pt_jesSubTotalScaleDown);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalScaleDown", Jet_mass_jesSubTotalScaleDown, &b_Jet_mass_jesSubTotalScaleDown);
   fChain->SetBranchAddress("MET_pt_jesSubTotalScaleDown", &MET_pt_jesSubTotalScaleDown, &b_MET_pt_jesSubTotalScaleDown);
   fChain->SetBranchAddress("MET_phi_jesSubTotalScaleDown", &MET_phi_jesSubTotalScaleDown, &b_MET_phi_jesSubTotalScaleDown);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalAbsoluteDown", Jet_pt_jesSubTotalAbsoluteDown, &b_Jet_pt_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalAbsoluteDown", Jet_mass_jesSubTotalAbsoluteDown, &b_Jet_mass_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("MET_pt_jesSubTotalAbsoluteDown", &MET_pt_jesSubTotalAbsoluteDown, &b_MET_pt_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("MET_phi_jesSubTotalAbsoluteDown", &MET_phi_jesSubTotalAbsoluteDown, &b_MET_phi_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalMCDown", Jet_pt_jesSubTotalMCDown, &b_Jet_pt_jesSubTotalMCDown);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalMCDown", Jet_mass_jesSubTotalMCDown, &b_Jet_mass_jesSubTotalMCDown);
   fChain->SetBranchAddress("MET_pt_jesSubTotalMCDown", &MET_pt_jesSubTotalMCDown, &b_MET_pt_jesSubTotalMCDown);
   fChain->SetBranchAddress("MET_phi_jesSubTotalMCDown", &MET_phi_jesSubTotalMCDown, &b_MET_phi_jesSubTotalMCDown);
   fChain->SetBranchAddress("Jet_pt_jesTotalDown", Jet_pt_jesTotalDown, &b_Jet_pt_jesTotalDown);
   fChain->SetBranchAddress("Jet_mass_jesTotalDown", Jet_mass_jesTotalDown, &b_Jet_mass_jesTotalDown);
   fChain->SetBranchAddress("MET_pt_jesTotalDown", &MET_pt_jesTotalDown, &b_MET_pt_jesTotalDown);
   fChain->SetBranchAddress("MET_phi_jesTotalDown", &MET_phi_jesTotalDown, &b_MET_phi_jesTotalDown);
   fChain->SetBranchAddress("Jet_pt_jesTotalNoFlavorDown", Jet_pt_jesTotalNoFlavorDown, &b_Jet_pt_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("Jet_mass_jesTotalNoFlavorDown", Jet_mass_jesTotalNoFlavorDown, &b_Jet_mass_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("MET_pt_jesTotalNoFlavorDown", &MET_pt_jesTotalNoFlavorDown, &b_MET_pt_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("MET_phi_jesTotalNoFlavorDown", &MET_phi_jesTotalNoFlavorDown, &b_MET_phi_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("Jet_pt_jesTotalNoTimeDown", Jet_pt_jesTotalNoTimeDown, &b_Jet_pt_jesTotalNoTimeDown);
   fChain->SetBranchAddress("Jet_mass_jesTotalNoTimeDown", Jet_mass_jesTotalNoTimeDown, &b_Jet_mass_jesTotalNoTimeDown);
   fChain->SetBranchAddress("MET_pt_jesTotalNoTimeDown", &MET_pt_jesTotalNoTimeDown, &b_MET_pt_jesTotalNoTimeDown);
   fChain->SetBranchAddress("MET_phi_jesTotalNoTimeDown", &MET_phi_jesTotalNoTimeDown, &b_MET_phi_jesTotalNoTimeDown);
   fChain->SetBranchAddress("Jet_pt_jesTotalNoFlavorNoTimeDown", Jet_pt_jesTotalNoFlavorNoTimeDown, &b_Jet_pt_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("Jet_mass_jesTotalNoFlavorNoTimeDown", Jet_mass_jesTotalNoFlavorNoTimeDown, &b_Jet_mass_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("MET_pt_jesTotalNoFlavorNoTimeDown", &MET_pt_jesTotalNoFlavorNoTimeDown, &b_MET_pt_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("MET_phi_jesTotalNoFlavorNoTimeDown", &MET_phi_jesTotalNoFlavorNoTimeDown, &b_MET_phi_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorZJetDown", Jet_pt_jesFlavorZJetDown, &b_Jet_pt_jesFlavorZJetDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorZJetDown", Jet_mass_jesFlavorZJetDown, &b_Jet_mass_jesFlavorZJetDown);
   fChain->SetBranchAddress("MET_pt_jesFlavorZJetDown", &MET_pt_jesFlavorZJetDown, &b_MET_pt_jesFlavorZJetDown);
   fChain->SetBranchAddress("MET_phi_jesFlavorZJetDown", &MET_phi_jesFlavorZJetDown, &b_MET_phi_jesFlavorZJetDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPhotonJetDown", Jet_pt_jesFlavorPhotonJetDown, &b_Jet_pt_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPhotonJetDown", Jet_mass_jesFlavorPhotonJetDown, &b_Jet_mass_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("MET_pt_jesFlavorPhotonJetDown", &MET_pt_jesFlavorPhotonJetDown, &b_MET_pt_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("MET_phi_jesFlavorPhotonJetDown", &MET_phi_jesFlavorPhotonJetDown, &b_MET_phi_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureGluonDown", Jet_pt_jesFlavorPureGluonDown, &b_Jet_pt_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureGluonDown", Jet_mass_jesFlavorPureGluonDown, &b_Jet_mass_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("MET_pt_jesFlavorPureGluonDown", &MET_pt_jesFlavorPureGluonDown, &b_MET_pt_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("MET_phi_jesFlavorPureGluonDown", &MET_phi_jesFlavorPureGluonDown, &b_MET_phi_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureQuarkDown", Jet_pt_jesFlavorPureQuarkDown, &b_Jet_pt_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureQuarkDown", Jet_mass_jesFlavorPureQuarkDown, &b_Jet_mass_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("MET_pt_jesFlavorPureQuarkDown", &MET_pt_jesFlavorPureQuarkDown, &b_MET_pt_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("MET_phi_jesFlavorPureQuarkDown", &MET_phi_jesFlavorPureQuarkDown, &b_MET_phi_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureCharmDown", Jet_pt_jesFlavorPureCharmDown, &b_Jet_pt_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureCharmDown", Jet_mass_jesFlavorPureCharmDown, &b_Jet_mass_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("MET_pt_jesFlavorPureCharmDown", &MET_pt_jesFlavorPureCharmDown, &b_MET_pt_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("MET_phi_jesFlavorPureCharmDown", &MET_phi_jesFlavorPureCharmDown, &b_MET_phi_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureBottomDown", Jet_pt_jesFlavorPureBottomDown, &b_Jet_pt_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureBottomDown", Jet_mass_jesFlavorPureBottomDown, &b_Jet_mass_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("MET_pt_jesFlavorPureBottomDown", &MET_pt_jesFlavorPureBottomDown, &b_MET_pt_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("MET_phi_jesFlavorPureBottomDown", &MET_phi_jesFlavorPureBottomDown, &b_MET_phi_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunBDown", Jet_pt_jesTimeRunBDown, &b_Jet_pt_jesTimeRunBDown);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunBDown", Jet_mass_jesTimeRunBDown, &b_Jet_mass_jesTimeRunBDown);
   fChain->SetBranchAddress("MET_pt_jesTimeRunBDown", &MET_pt_jesTimeRunBDown, &b_MET_pt_jesTimeRunBDown);
   fChain->SetBranchAddress("MET_phi_jesTimeRunBDown", &MET_phi_jesTimeRunBDown, &b_MET_phi_jesTimeRunBDown);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunCDown", Jet_pt_jesTimeRunCDown, &b_Jet_pt_jesTimeRunCDown);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunCDown", Jet_mass_jesTimeRunCDown, &b_Jet_mass_jesTimeRunCDown);
   fChain->SetBranchAddress("MET_pt_jesTimeRunCDown", &MET_pt_jesTimeRunCDown, &b_MET_pt_jesTimeRunCDown);
   fChain->SetBranchAddress("MET_phi_jesTimeRunCDown", &MET_phi_jesTimeRunCDown, &b_MET_phi_jesTimeRunCDown);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunDEDown", Jet_pt_jesTimeRunDEDown, &b_Jet_pt_jesTimeRunDEDown);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunDEDown", Jet_mass_jesTimeRunDEDown, &b_Jet_mass_jesTimeRunDEDown);
   fChain->SetBranchAddress("MET_pt_jesTimeRunDEDown", &MET_pt_jesTimeRunDEDown, &b_MET_pt_jesTimeRunDEDown);
   fChain->SetBranchAddress("MET_phi_jesTimeRunDEDown", &MET_phi_jesTimeRunDEDown, &b_MET_phi_jesTimeRunDEDown);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunFDown", Jet_pt_jesTimeRunFDown, &b_Jet_pt_jesTimeRunFDown);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunFDown", Jet_mass_jesTimeRunFDown, &b_Jet_mass_jesTimeRunFDown);
   fChain->SetBranchAddress("MET_pt_jesTimeRunFDown", &MET_pt_jesTimeRunFDown, &b_MET_pt_jesTimeRunFDown);
   fChain->SetBranchAddress("MET_phi_jesTimeRunFDown", &MET_phi_jesTimeRunFDown, &b_MET_phi_jesTimeRunFDown);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupMPFInSituDown", Jet_pt_jesCorrelationGroupMPFInSituDown, &b_Jet_pt_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupMPFInSituDown", Jet_mass_jesCorrelationGroupMPFInSituDown, &b_Jet_mass_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupMPFInSituDown", &MET_pt_jesCorrelationGroupMPFInSituDown, &b_MET_pt_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupMPFInSituDown", &MET_phi_jesCorrelationGroupMPFInSituDown, &b_MET_phi_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupIntercalibrationDown", Jet_pt_jesCorrelationGroupIntercalibrationDown, &b_Jet_pt_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupIntercalibrationDown", Jet_mass_jesCorrelationGroupIntercalibrationDown, &b_Jet_mass_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupIntercalibrationDown", &MET_pt_jesCorrelationGroupIntercalibrationDown, &b_MET_pt_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupIntercalibrationDown", &MET_phi_jesCorrelationGroupIntercalibrationDown, &b_MET_phi_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupbJESDown", Jet_pt_jesCorrelationGroupbJESDown, &b_Jet_pt_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupbJESDown", Jet_mass_jesCorrelationGroupbJESDown, &b_Jet_mass_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupbJESDown", &MET_pt_jesCorrelationGroupbJESDown, &b_MET_pt_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupbJESDown", &MET_phi_jesCorrelationGroupbJESDown, &b_MET_phi_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupFlavorDown", Jet_pt_jesCorrelationGroupFlavorDown, &b_Jet_pt_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupFlavorDown", Jet_mass_jesCorrelationGroupFlavorDown, &b_Jet_mass_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupFlavorDown", &MET_pt_jesCorrelationGroupFlavorDown, &b_MET_pt_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupFlavorDown", &MET_phi_jesCorrelationGroupFlavorDown, &b_MET_phi_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupUncorrelatedDown", Jet_pt_jesCorrelationGroupUncorrelatedDown, &b_Jet_pt_jesCorrelationGroupUncorrelatedDown);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupUncorrelatedDown", Jet_mass_jesCorrelationGroupUncorrelatedDown, &b_Jet_mass_jesCorrelationGroupUncorrelatedDown);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupUncorrelatedDown", &MET_pt_jesCorrelationGroupUncorrelatedDown, &b_MET_pt_jesCorrelationGroupUncorrelatedDown);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupUncorrelatedDown", &MET_phi_jesCorrelationGroupUncorrelatedDown, &b_MET_phi_jesCorrelationGroupUncorrelatedDown);
   fChain->SetBranchAddress("MET_pt_unclustEnDown", &MET_pt_unclustEnDown, &b_MET_pt_unclustEnDown);
   fChain->SetBranchAddress("MET_phi_unclustEnDown", &MET_phi_unclustEnDown, &b_MET_phi_unclustEnDown);
   
   fChain->SetBranchAddress("FatJet_pt_raw", FatJet_pt_raw, &b_FatJet_pt_raw);
   fChain->SetBranchAddress("FatJet_pt_nom", FatJet_pt_nom, &b_FatJet_pt_nom);
   fChain->SetBranchAddress("FatJet_mass_raw", FatJet_mass_raw, &b_FatJet_mass_raw);
   fChain->SetBranchAddress("FatJet_mass_nom", FatJet_mass_nom, &b_FatJet_mass_nom);
   fChain->SetBranchAddress("FatJet_corr_JEC", FatJet_corr_JEC, &b_FatJet_corr_JEC);
   fChain->SetBranchAddress("FatJet_corr_JER", FatJet_corr_JER, &b_FatJet_corr_JER);
   fChain->SetBranchAddress("FatJet_corr_JMS", FatJet_corr_JMS, &b_FatJet_corr_JMS);
   fChain->SetBranchAddress("FatJet_corr_JMR", FatJet_corr_JMR, &b_FatJet_corr_JMR);
   fChain->SetBranchAddress("FatJet_msoftdrop_raw", FatJet_msoftdrop_raw, &b_FatJet_msoftdrop_raw);
   fChain->SetBranchAddress("FatJet_msoftdrop_nom", FatJet_msoftdrop_nom, &b_FatJet_msoftdrop_nom);
   fChain->SetBranchAddress("FatJet_msoftdrop_corr_JMR", FatJet_msoftdrop_corr_JMR, &b_FatJet_msoftdrop_corr_JMR);
   fChain->SetBranchAddress("FatJet_msoftdrop_corr_JMS", FatJet_msoftdrop_corr_JMS, &b_FatJet_msoftdrop_corr_JMS);
   fChain->SetBranchAddress("FatJet_msoftdrop_corr_PUPPI", FatJet_msoftdrop_corr_PUPPI, &b_FatJet_msoftdrop_corr_PUPPI);
   fChain->SetBranchAddress("FatJet_msoftdrop_tau21DDT_nom", FatJet_msoftdrop_tau21DDT_nom, &b_FatJet_msoftdrop_tau21DDT_nom);

   fChain->SetBranchAddress("FatJet_pt_jerUp", FatJet_pt_jerUp, &b_FatJet_pt_jerUp);
   fChain->SetBranchAddress("FatJet_mass_jerUp", FatJet_mass_jerUp, &b_FatJet_mass_jerUp);
   fChain->SetBranchAddress("FatJet_mass_jmrUp", FatJet_mass_jmrUp, &b_FatJet_mass_jmrUp);
   fChain->SetBranchAddress("FatJet_mass_jmsUp", FatJet_mass_jmsUp, &b_FatJet_mass_jmsUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jerUp", FatJet_msoftdrop_jerUp, &b_FatJet_msoftdrop_jerUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jmrUp", FatJet_msoftdrop_jmrUp, &b_FatJet_msoftdrop_jmrUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jmsUp", FatJet_msoftdrop_jmsUp, &b_FatJet_msoftdrop_jmsUp);
   fChain->SetBranchAddress("FatJet_pt_jesAbsoluteStatUp", FatJet_pt_jesAbsoluteStatUp, &b_FatJet_pt_jesAbsoluteStatUp);
   fChain->SetBranchAddress("FatJet_mass_jesAbsoluteStatUp", FatJet_mass_jesAbsoluteStatUp, &b_FatJet_mass_jesAbsoluteStatUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteStatUp", FatJet_msoftdrop_jesAbsoluteStatUp, &b_FatJet_msoftdrop_jesAbsoluteStatUp);
   fChain->SetBranchAddress("FatJet_pt_jesAbsoluteScaleUp", FatJet_pt_jesAbsoluteScaleUp, &b_FatJet_pt_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("FatJet_mass_jesAbsoluteScaleUp", FatJet_mass_jesAbsoluteScaleUp, &b_FatJet_mass_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteScaleUp", FatJet_msoftdrop_jesAbsoluteScaleUp, &b_FatJet_msoftdrop_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("FatJet_pt_jesAbsoluteFlavMapUp", FatJet_pt_jesAbsoluteFlavMapUp, &b_FatJet_pt_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("FatJet_mass_jesAbsoluteFlavMapUp", FatJet_mass_jesAbsoluteFlavMapUp, &b_FatJet_mass_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteFlavMapUp", FatJet_msoftdrop_jesAbsoluteFlavMapUp, &b_FatJet_msoftdrop_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("FatJet_pt_jesAbsoluteMPFBiasUp", FatJet_pt_jesAbsoluteMPFBiasUp, &b_FatJet_pt_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("FatJet_mass_jesAbsoluteMPFBiasUp", FatJet_mass_jesAbsoluteMPFBiasUp, &b_FatJet_mass_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteMPFBiasUp", FatJet_msoftdrop_jesAbsoluteMPFBiasUp, &b_FatJet_msoftdrop_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("FatJet_pt_jesFragmentationUp", FatJet_pt_jesFragmentationUp, &b_FatJet_pt_jesFragmentationUp);
   fChain->SetBranchAddress("FatJet_mass_jesFragmentationUp", FatJet_mass_jesFragmentationUp, &b_FatJet_mass_jesFragmentationUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFragmentationUp", FatJet_msoftdrop_jesFragmentationUp, &b_FatJet_msoftdrop_jesFragmentationUp);
   fChain->SetBranchAddress("FatJet_pt_jesSinglePionECALUp", FatJet_pt_jesSinglePionECALUp, &b_FatJet_pt_jesSinglePionECALUp);
   fChain->SetBranchAddress("FatJet_mass_jesSinglePionECALUp", FatJet_mass_jesSinglePionECALUp, &b_FatJet_mass_jesSinglePionECALUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSinglePionECALUp", FatJet_msoftdrop_jesSinglePionECALUp, &b_FatJet_msoftdrop_jesSinglePionECALUp);
   fChain->SetBranchAddress("FatJet_pt_jesSinglePionHCALUp", FatJet_pt_jesSinglePionHCALUp, &b_FatJet_pt_jesSinglePionHCALUp);
   fChain->SetBranchAddress("FatJet_mass_jesSinglePionHCALUp", FatJet_mass_jesSinglePionHCALUp, &b_FatJet_mass_jesSinglePionHCALUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSinglePionHCALUp", FatJet_msoftdrop_jesSinglePionHCALUp, &b_FatJet_msoftdrop_jesSinglePionHCALUp);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorQCDUp", FatJet_pt_jesFlavorQCDUp, &b_FatJet_pt_jesFlavorQCDUp);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorQCDUp", FatJet_mass_jesFlavorQCDUp, &b_FatJet_mass_jesFlavorQCDUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorQCDUp", FatJet_msoftdrop_jesFlavorQCDUp, &b_FatJet_msoftdrop_jesFlavorQCDUp);
   fChain->SetBranchAddress("FatJet_pt_jesTimePtEtaUp", FatJet_pt_jesTimePtEtaUp, &b_FatJet_pt_jesTimePtEtaUp);
   fChain->SetBranchAddress("FatJet_mass_jesTimePtEtaUp", FatJet_mass_jesTimePtEtaUp, &b_FatJet_mass_jesTimePtEtaUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimePtEtaUp", FatJet_msoftdrop_jesTimePtEtaUp, &b_FatJet_msoftdrop_jesTimePtEtaUp);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeJEREC1Up", FatJet_pt_jesRelativeJEREC1Up, &b_FatJet_pt_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeJEREC1Up", FatJet_mass_jesRelativeJEREC1Up, &b_FatJet_mass_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeJEREC1Up", FatJet_msoftdrop_jesRelativeJEREC1Up, &b_FatJet_msoftdrop_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeJEREC2Up", FatJet_pt_jesRelativeJEREC2Up, &b_FatJet_pt_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeJEREC2Up", FatJet_mass_jesRelativeJEREC2Up, &b_FatJet_mass_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeJEREC2Up", FatJet_msoftdrop_jesRelativeJEREC2Up, &b_FatJet_msoftdrop_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeJERHFUp", FatJet_pt_jesRelativeJERHFUp, &b_FatJet_pt_jesRelativeJERHFUp);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeJERHFUp", FatJet_mass_jesRelativeJERHFUp, &b_FatJet_mass_jesRelativeJERHFUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeJERHFUp", FatJet_msoftdrop_jesRelativeJERHFUp, &b_FatJet_msoftdrop_jesRelativeJERHFUp);
   fChain->SetBranchAddress("FatJet_pt_jesRelativePtBBUp", FatJet_pt_jesRelativePtBBUp, &b_FatJet_pt_jesRelativePtBBUp);
   fChain->SetBranchAddress("FatJet_mass_jesRelativePtBBUp", FatJet_mass_jesRelativePtBBUp, &b_FatJet_mass_jesRelativePtBBUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativePtBBUp", FatJet_msoftdrop_jesRelativePtBBUp, &b_FatJet_msoftdrop_jesRelativePtBBUp);
   fChain->SetBranchAddress("FatJet_pt_jesRelativePtEC1Up", FatJet_pt_jesRelativePtEC1Up, &b_FatJet_pt_jesRelativePtEC1Up);
   fChain->SetBranchAddress("FatJet_mass_jesRelativePtEC1Up", FatJet_mass_jesRelativePtEC1Up, &b_FatJet_mass_jesRelativePtEC1Up);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativePtEC1Up", FatJet_msoftdrop_jesRelativePtEC1Up, &b_FatJet_msoftdrop_jesRelativePtEC1Up);
   fChain->SetBranchAddress("FatJet_pt_jesRelativePtEC2Up", FatJet_pt_jesRelativePtEC2Up, &b_FatJet_pt_jesRelativePtEC2Up);
   fChain->SetBranchAddress("FatJet_mass_jesRelativePtEC2Up", FatJet_mass_jesRelativePtEC2Up, &b_FatJet_mass_jesRelativePtEC2Up);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativePtEC2Up", FatJet_msoftdrop_jesRelativePtEC2Up, &b_FatJet_msoftdrop_jesRelativePtEC2Up);
   fChain->SetBranchAddress("FatJet_pt_jesRelativePtHFUp", FatJet_pt_jesRelativePtHFUp, &b_FatJet_pt_jesRelativePtHFUp);
   fChain->SetBranchAddress("FatJet_mass_jesRelativePtHFUp", FatJet_mass_jesRelativePtHFUp, &b_FatJet_mass_jesRelativePtHFUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativePtHFUp", FatJet_msoftdrop_jesRelativePtHFUp, &b_FatJet_msoftdrop_jesRelativePtHFUp);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeBalUp", FatJet_pt_jesRelativeBalUp, &b_FatJet_pt_jesRelativeBalUp);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeBalUp", FatJet_mass_jesRelativeBalUp, &b_FatJet_mass_jesRelativeBalUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeBalUp", FatJet_msoftdrop_jesRelativeBalUp, &b_FatJet_msoftdrop_jesRelativeBalUp);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeSampleUp", FatJet_pt_jesRelativeSampleUp, &b_FatJet_pt_jesRelativeSampleUp);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeSampleUp", FatJet_mass_jesRelativeSampleUp, &b_FatJet_mass_jesRelativeSampleUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeSampleUp", FatJet_msoftdrop_jesRelativeSampleUp, &b_FatJet_msoftdrop_jesRelativeSampleUp);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeFSRUp", FatJet_pt_jesRelativeFSRUp, &b_FatJet_pt_jesRelativeFSRUp);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeFSRUp", FatJet_mass_jesRelativeFSRUp, &b_FatJet_mass_jesRelativeFSRUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeFSRUp", FatJet_msoftdrop_jesRelativeFSRUp, &b_FatJet_msoftdrop_jesRelativeFSRUp);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeStatFSRUp", FatJet_pt_jesRelativeStatFSRUp, &b_FatJet_pt_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeStatFSRUp", FatJet_mass_jesRelativeStatFSRUp, &b_FatJet_mass_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeStatFSRUp", FatJet_msoftdrop_jesRelativeStatFSRUp, &b_FatJet_msoftdrop_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeStatECUp", FatJet_pt_jesRelativeStatECUp, &b_FatJet_pt_jesRelativeStatECUp);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeStatECUp", FatJet_mass_jesRelativeStatECUp, &b_FatJet_mass_jesRelativeStatECUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeStatECUp", FatJet_msoftdrop_jesRelativeStatECUp, &b_FatJet_msoftdrop_jesRelativeStatECUp);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeStatHFUp", FatJet_pt_jesRelativeStatHFUp, &b_FatJet_pt_jesRelativeStatHFUp);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeStatHFUp", FatJet_mass_jesRelativeStatHFUp, &b_FatJet_mass_jesRelativeStatHFUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeStatHFUp", FatJet_msoftdrop_jesRelativeStatHFUp, &b_FatJet_msoftdrop_jesRelativeStatHFUp);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpDataMCUp", FatJet_pt_jesPileUpDataMCUp, &b_FatJet_pt_jesPileUpDataMCUp);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpDataMCUp", FatJet_mass_jesPileUpDataMCUp, &b_FatJet_mass_jesPileUpDataMCUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpDataMCUp", FatJet_msoftdrop_jesPileUpDataMCUp, &b_FatJet_msoftdrop_jesPileUpDataMCUp);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtRefUp", FatJet_pt_jesPileUpPtRefUp, &b_FatJet_pt_jesPileUpPtRefUp);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtRefUp", FatJet_mass_jesPileUpPtRefUp, &b_FatJet_mass_jesPileUpPtRefUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtRefUp", FatJet_msoftdrop_jesPileUpPtRefUp, &b_FatJet_msoftdrop_jesPileUpPtRefUp);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtBBUp", FatJet_pt_jesPileUpPtBBUp, &b_FatJet_pt_jesPileUpPtBBUp);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtBBUp", FatJet_mass_jesPileUpPtBBUp, &b_FatJet_mass_jesPileUpPtBBUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtBBUp", FatJet_msoftdrop_jesPileUpPtBBUp, &b_FatJet_msoftdrop_jesPileUpPtBBUp);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtEC1Up", FatJet_pt_jesPileUpPtEC1Up, &b_FatJet_pt_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtEC1Up", FatJet_mass_jesPileUpPtEC1Up, &b_FatJet_mass_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtEC1Up", FatJet_msoftdrop_jesPileUpPtEC1Up, &b_FatJet_msoftdrop_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtEC2Up", FatJet_pt_jesPileUpPtEC2Up, &b_FatJet_pt_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtEC2Up", FatJet_mass_jesPileUpPtEC2Up, &b_FatJet_mass_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtEC2Up", FatJet_msoftdrop_jesPileUpPtEC2Up, &b_FatJet_msoftdrop_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtHFUp", FatJet_pt_jesPileUpPtHFUp, &b_FatJet_pt_jesPileUpPtHFUp);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtHFUp", FatJet_mass_jesPileUpPtHFUp, &b_FatJet_mass_jesPileUpPtHFUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtHFUp", FatJet_msoftdrop_jesPileUpPtHFUp, &b_FatJet_msoftdrop_jesPileUpPtHFUp);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpMuZeroUp", FatJet_pt_jesPileUpMuZeroUp, &b_FatJet_pt_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpMuZeroUp", FatJet_mass_jesPileUpMuZeroUp, &b_FatJet_mass_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpMuZeroUp", FatJet_msoftdrop_jesPileUpMuZeroUp, &b_FatJet_msoftdrop_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpEnvelopeUp", FatJet_pt_jesPileUpEnvelopeUp, &b_FatJet_pt_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpEnvelopeUp", FatJet_mass_jesPileUpEnvelopeUp, &b_FatJet_mass_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpEnvelopeUp", FatJet_msoftdrop_jesPileUpEnvelopeUp, &b_FatJet_msoftdrop_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalPileUpUp", FatJet_pt_jesSubTotalPileUpUp, &b_FatJet_pt_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalPileUpUp", FatJet_mass_jesSubTotalPileUpUp, &b_FatJet_mass_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalPileUpUp", FatJet_msoftdrop_jesSubTotalPileUpUp, &b_FatJet_msoftdrop_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalRelativeUp", FatJet_pt_jesSubTotalRelativeUp, &b_FatJet_pt_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalRelativeUp", FatJet_mass_jesSubTotalRelativeUp, &b_FatJet_mass_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalRelativeUp", FatJet_msoftdrop_jesSubTotalRelativeUp, &b_FatJet_msoftdrop_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalPtUp", FatJet_pt_jesSubTotalPtUp, &b_FatJet_pt_jesSubTotalPtUp);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalPtUp", FatJet_mass_jesSubTotalPtUp, &b_FatJet_mass_jesSubTotalPtUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalPtUp", FatJet_msoftdrop_jesSubTotalPtUp, &b_FatJet_msoftdrop_jesSubTotalPtUp);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalScaleUp", FatJet_pt_jesSubTotalScaleUp, &b_FatJet_pt_jesSubTotalScaleUp);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalScaleUp", FatJet_mass_jesSubTotalScaleUp, &b_FatJet_mass_jesSubTotalScaleUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalScaleUp", FatJet_msoftdrop_jesSubTotalScaleUp, &b_FatJet_msoftdrop_jesSubTotalScaleUp);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalAbsoluteUp", FatJet_pt_jesSubTotalAbsoluteUp, &b_FatJet_pt_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalAbsoluteUp", FatJet_mass_jesSubTotalAbsoluteUp, &b_FatJet_mass_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalAbsoluteUp", FatJet_msoftdrop_jesSubTotalAbsoluteUp, &b_FatJet_msoftdrop_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalMCUp", FatJet_pt_jesSubTotalMCUp, &b_FatJet_pt_jesSubTotalMCUp);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalMCUp", FatJet_mass_jesSubTotalMCUp, &b_FatJet_mass_jesSubTotalMCUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalMCUp", FatJet_msoftdrop_jesSubTotalMCUp, &b_FatJet_msoftdrop_jesSubTotalMCUp);
   fChain->SetBranchAddress("FatJet_pt_jesTotalUp", FatJet_pt_jesTotalUp, &b_FatJet_pt_jesTotalUp);
   fChain->SetBranchAddress("FatJet_mass_jesTotalUp", FatJet_mass_jesTotalUp, &b_FatJet_mass_jesTotalUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTotalUp", FatJet_msoftdrop_jesTotalUp, &b_FatJet_msoftdrop_jesTotalUp);
   fChain->SetBranchAddress("FatJet_pt_jesTotalNoFlavorUp", FatJet_pt_jesTotalNoFlavorUp, &b_FatJet_pt_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("FatJet_mass_jesTotalNoFlavorUp", FatJet_mass_jesTotalNoFlavorUp, &b_FatJet_mass_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTotalNoFlavorUp", FatJet_msoftdrop_jesTotalNoFlavorUp, &b_FatJet_msoftdrop_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("FatJet_pt_jesTotalNoTimeUp", FatJet_pt_jesTotalNoTimeUp, &b_FatJet_pt_jesTotalNoTimeUp);
   fChain->SetBranchAddress("FatJet_mass_jesTotalNoTimeUp", FatJet_mass_jesTotalNoTimeUp, &b_FatJet_mass_jesTotalNoTimeUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTotalNoTimeUp", FatJet_msoftdrop_jesTotalNoTimeUp, &b_FatJet_msoftdrop_jesTotalNoTimeUp);
   fChain->SetBranchAddress("FatJet_pt_jesTotalNoFlavorNoTimeUp", FatJet_pt_jesTotalNoFlavorNoTimeUp, &b_FatJet_pt_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("FatJet_mass_jesTotalNoFlavorNoTimeUp", FatJet_mass_jesTotalNoFlavorNoTimeUp, &b_FatJet_mass_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTotalNoFlavorNoTimeUp", FatJet_msoftdrop_jesTotalNoFlavorNoTimeUp, &b_FatJet_msoftdrop_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorZJetUp", FatJet_pt_jesFlavorZJetUp, &b_FatJet_pt_jesFlavorZJetUp);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorZJetUp", FatJet_mass_jesFlavorZJetUp, &b_FatJet_mass_jesFlavorZJetUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorZJetUp", FatJet_msoftdrop_jesFlavorZJetUp, &b_FatJet_msoftdrop_jesFlavorZJetUp);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPhotonJetUp", FatJet_pt_jesFlavorPhotonJetUp, &b_FatJet_pt_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPhotonJetUp", FatJet_mass_jesFlavorPhotonJetUp, &b_FatJet_mass_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPhotonJetUp", FatJet_msoftdrop_jesFlavorPhotonJetUp, &b_FatJet_msoftdrop_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPureGluonUp", FatJet_pt_jesFlavorPureGluonUp, &b_FatJet_pt_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPureGluonUp", FatJet_mass_jesFlavorPureGluonUp, &b_FatJet_mass_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureGluonUp", FatJet_msoftdrop_jesFlavorPureGluonUp, &b_FatJet_msoftdrop_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPureQuarkUp", FatJet_pt_jesFlavorPureQuarkUp, &b_FatJet_pt_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPureQuarkUp", FatJet_mass_jesFlavorPureQuarkUp, &b_FatJet_mass_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureQuarkUp", FatJet_msoftdrop_jesFlavorPureQuarkUp, &b_FatJet_msoftdrop_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPureCharmUp", FatJet_pt_jesFlavorPureCharmUp, &b_FatJet_pt_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPureCharmUp", FatJet_mass_jesFlavorPureCharmUp, &b_FatJet_mass_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureCharmUp", FatJet_msoftdrop_jesFlavorPureCharmUp, &b_FatJet_msoftdrop_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPureBottomUp", FatJet_pt_jesFlavorPureBottomUp, &b_FatJet_pt_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPureBottomUp", FatJet_mass_jesFlavorPureBottomUp, &b_FatJet_mass_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureBottomUp", FatJet_msoftdrop_jesFlavorPureBottomUp, &b_FatJet_msoftdrop_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("FatJet_pt_jesTimeRunBUp", FatJet_pt_jesTimeRunBUp, &b_FatJet_pt_jesTimeRunBUp);
   fChain->SetBranchAddress("FatJet_mass_jesTimeRunBUp", FatJet_mass_jesTimeRunBUp, &b_FatJet_mass_jesTimeRunBUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimeRunBUp", FatJet_msoftdrop_jesTimeRunBUp, &b_FatJet_msoftdrop_jesTimeRunBUp);
   fChain->SetBranchAddress("FatJet_pt_jesTimeRunCUp", FatJet_pt_jesTimeRunCUp, &b_FatJet_pt_jesTimeRunCUp);
   fChain->SetBranchAddress("FatJet_mass_jesTimeRunCUp", FatJet_mass_jesTimeRunCUp, &b_FatJet_mass_jesTimeRunCUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimeRunCUp", FatJet_msoftdrop_jesTimeRunCUp, &b_FatJet_msoftdrop_jesTimeRunCUp);
   fChain->SetBranchAddress("FatJet_pt_jesTimeRunDEUp", FatJet_pt_jesTimeRunDEUp, &b_FatJet_pt_jesTimeRunDEUp);
   fChain->SetBranchAddress("FatJet_mass_jesTimeRunDEUp", FatJet_mass_jesTimeRunDEUp, &b_FatJet_mass_jesTimeRunDEUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimeRunDEUp", FatJet_msoftdrop_jesTimeRunDEUp, &b_FatJet_msoftdrop_jesTimeRunDEUp);
   fChain->SetBranchAddress("FatJet_pt_jesTimeRunFUp", FatJet_pt_jesTimeRunFUp, &b_FatJet_pt_jesTimeRunFUp);
   fChain->SetBranchAddress("FatJet_mass_jesTimeRunFUp", FatJet_mass_jesTimeRunFUp, &b_FatJet_mass_jesTimeRunFUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimeRunFUp", FatJet_msoftdrop_jesTimeRunFUp, &b_FatJet_msoftdrop_jesTimeRunFUp);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupMPFInSituUp", FatJet_pt_jesCorrelationGroupMPFInSituUp, &b_FatJet_pt_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupMPFInSituUp", FatJet_mass_jesCorrelationGroupMPFInSituUp, &b_FatJet_mass_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupMPFInSituUp", FatJet_msoftdrop_jesCorrelationGroupMPFInSituUp, &b_FatJet_msoftdrop_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupIntercalibrationUp", FatJet_pt_jesCorrelationGroupIntercalibrationUp, &b_FatJet_pt_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupIntercalibrationUp", FatJet_mass_jesCorrelationGroupIntercalibrationUp, &b_FatJet_mass_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupIntercalibrationUp", FatJet_msoftdrop_jesCorrelationGroupIntercalibrationUp, &b_FatJet_msoftdrop_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupbJESUp", FatJet_pt_jesCorrelationGroupbJESUp, &b_FatJet_pt_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupbJESUp", FatJet_mass_jesCorrelationGroupbJESUp, &b_FatJet_mass_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupbJESUp", FatJet_msoftdrop_jesCorrelationGroupbJESUp, &b_FatJet_msoftdrop_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupFlavorUp", FatJet_pt_jesCorrelationGroupFlavorUp, &b_FatJet_pt_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupFlavorUp", FatJet_mass_jesCorrelationGroupFlavorUp, &b_FatJet_mass_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupFlavorUp", FatJet_msoftdrop_jesCorrelationGroupFlavorUp, &b_FatJet_msoftdrop_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupUncorrelatedUp", FatJet_pt_jesCorrelationGroupUncorrelatedUp, &b_FatJet_pt_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupUncorrelatedUp", FatJet_mass_jesCorrelationGroupUncorrelatedUp, &b_FatJet_mass_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupUncorrelatedUp", FatJet_msoftdrop_jesCorrelationGroupUncorrelatedUp, &b_FatJet_msoftdrop_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("FatJet_pt_jerDown", FatJet_pt_jerDown, &b_FatJet_pt_jerDown);
   fChain->SetBranchAddress("FatJet_mass_jerDown", FatJet_mass_jerDown, &b_FatJet_mass_jerDown);
   fChain->SetBranchAddress("FatJet_mass_jmrDown", FatJet_mass_jmrDown, &b_FatJet_mass_jmrDown);
   fChain->SetBranchAddress("FatJet_mass_jmsDown", FatJet_mass_jmsDown, &b_FatJet_mass_jmsDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jerDown", FatJet_msoftdrop_jerDown, &b_FatJet_msoftdrop_jerDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jmrDown", FatJet_msoftdrop_jmrDown, &b_FatJet_msoftdrop_jmrDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jmsDown", FatJet_msoftdrop_jmsDown, &b_FatJet_msoftdrop_jmsDown);
   fChain->SetBranchAddress("FatJet_pt_jesAbsoluteStatDown", FatJet_pt_jesAbsoluteStatDown, &b_FatJet_pt_jesAbsoluteStatDown);
   fChain->SetBranchAddress("FatJet_mass_jesAbsoluteStatDown", FatJet_mass_jesAbsoluteStatDown, &b_FatJet_mass_jesAbsoluteStatDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteStatDown", FatJet_msoftdrop_jesAbsoluteStatDown, &b_FatJet_msoftdrop_jesAbsoluteStatDown);
   fChain->SetBranchAddress("FatJet_pt_jesAbsoluteScaleDown", FatJet_pt_jesAbsoluteScaleDown, &b_FatJet_pt_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("FatJet_mass_jesAbsoluteScaleDown", FatJet_mass_jesAbsoluteScaleDown, &b_FatJet_mass_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteScaleDown", FatJet_msoftdrop_jesAbsoluteScaleDown, &b_FatJet_msoftdrop_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("FatJet_pt_jesAbsoluteFlavMapDown", FatJet_pt_jesAbsoluteFlavMapDown, &b_FatJet_pt_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("FatJet_mass_jesAbsoluteFlavMapDown", FatJet_mass_jesAbsoluteFlavMapDown, &b_FatJet_mass_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteFlavMapDown", FatJet_msoftdrop_jesAbsoluteFlavMapDown, &b_FatJet_msoftdrop_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("FatJet_pt_jesAbsoluteMPFBiasDown", FatJet_pt_jesAbsoluteMPFBiasDown, &b_FatJet_pt_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("FatJet_mass_jesAbsoluteMPFBiasDown", FatJet_mass_jesAbsoluteMPFBiasDown, &b_FatJet_mass_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteMPFBiasDown", FatJet_msoftdrop_jesAbsoluteMPFBiasDown, &b_FatJet_msoftdrop_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("FatJet_pt_jesFragmentationDown", FatJet_pt_jesFragmentationDown, &b_FatJet_pt_jesFragmentationDown);
   fChain->SetBranchAddress("FatJet_mass_jesFragmentationDown", FatJet_mass_jesFragmentationDown, &b_FatJet_mass_jesFragmentationDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFragmentationDown", FatJet_msoftdrop_jesFragmentationDown, &b_FatJet_msoftdrop_jesFragmentationDown);
   fChain->SetBranchAddress("FatJet_pt_jesSinglePionECALDown", FatJet_pt_jesSinglePionECALDown, &b_FatJet_pt_jesSinglePionECALDown);
   fChain->SetBranchAddress("FatJet_mass_jesSinglePionECALDown", FatJet_mass_jesSinglePionECALDown, &b_FatJet_mass_jesSinglePionECALDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSinglePionECALDown", FatJet_msoftdrop_jesSinglePionECALDown, &b_FatJet_msoftdrop_jesSinglePionECALDown);
   fChain->SetBranchAddress("FatJet_pt_jesSinglePionHCALDown", FatJet_pt_jesSinglePionHCALDown, &b_FatJet_pt_jesSinglePionHCALDown);
   fChain->SetBranchAddress("FatJet_mass_jesSinglePionHCALDown", FatJet_mass_jesSinglePionHCALDown, &b_FatJet_mass_jesSinglePionHCALDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSinglePionHCALDown", FatJet_msoftdrop_jesSinglePionHCALDown, &b_FatJet_msoftdrop_jesSinglePionHCALDown);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorQCDDown", FatJet_pt_jesFlavorQCDDown, &b_FatJet_pt_jesFlavorQCDDown);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorQCDDown", FatJet_mass_jesFlavorQCDDown, &b_FatJet_mass_jesFlavorQCDDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorQCDDown", FatJet_msoftdrop_jesFlavorQCDDown, &b_FatJet_msoftdrop_jesFlavorQCDDown);
   fChain->SetBranchAddress("FatJet_pt_jesTimePtEtaDown", FatJet_pt_jesTimePtEtaDown, &b_FatJet_pt_jesTimePtEtaDown);
   fChain->SetBranchAddress("FatJet_mass_jesTimePtEtaDown", FatJet_mass_jesTimePtEtaDown, &b_FatJet_mass_jesTimePtEtaDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimePtEtaDown", FatJet_msoftdrop_jesTimePtEtaDown, &b_FatJet_msoftdrop_jesTimePtEtaDown);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeJEREC1Down", FatJet_pt_jesRelativeJEREC1Down, &b_FatJet_pt_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeJEREC1Down", FatJet_mass_jesRelativeJEREC1Down, &b_FatJet_mass_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeJEREC1Down", FatJet_msoftdrop_jesRelativeJEREC1Down, &b_FatJet_msoftdrop_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeJEREC2Down", FatJet_pt_jesRelativeJEREC2Down, &b_FatJet_pt_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeJEREC2Down", FatJet_mass_jesRelativeJEREC2Down, &b_FatJet_mass_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeJEREC2Down", FatJet_msoftdrop_jesRelativeJEREC2Down, &b_FatJet_msoftdrop_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeJERHFDown", FatJet_pt_jesRelativeJERHFDown, &b_FatJet_pt_jesRelativeJERHFDown);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeJERHFDown", FatJet_mass_jesRelativeJERHFDown, &b_FatJet_mass_jesRelativeJERHFDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeJERHFDown", FatJet_msoftdrop_jesRelativeJERHFDown, &b_FatJet_msoftdrop_jesRelativeJERHFDown);
   fChain->SetBranchAddress("FatJet_pt_jesRelativePtBBDown", FatJet_pt_jesRelativePtBBDown, &b_FatJet_pt_jesRelativePtBBDown);
   fChain->SetBranchAddress("FatJet_mass_jesRelativePtBBDown", FatJet_mass_jesRelativePtBBDown, &b_FatJet_mass_jesRelativePtBBDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativePtBBDown", FatJet_msoftdrop_jesRelativePtBBDown, &b_FatJet_msoftdrop_jesRelativePtBBDown);
   fChain->SetBranchAddress("FatJet_pt_jesRelativePtEC1Down", FatJet_pt_jesRelativePtEC1Down, &b_FatJet_pt_jesRelativePtEC1Down);
   fChain->SetBranchAddress("FatJet_mass_jesRelativePtEC1Down", FatJet_mass_jesRelativePtEC1Down, &b_FatJet_mass_jesRelativePtEC1Down);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativePtEC1Down", FatJet_msoftdrop_jesRelativePtEC1Down, &b_FatJet_msoftdrop_jesRelativePtEC1Down);
   fChain->SetBranchAddress("FatJet_pt_jesRelativePtEC2Down", FatJet_pt_jesRelativePtEC2Down, &b_FatJet_pt_jesRelativePtEC2Down);
   fChain->SetBranchAddress("FatJet_mass_jesRelativePtEC2Down", FatJet_mass_jesRelativePtEC2Down, &b_FatJet_mass_jesRelativePtEC2Down);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativePtEC2Down", FatJet_msoftdrop_jesRelativePtEC2Down, &b_FatJet_msoftdrop_jesRelativePtEC2Down);
   fChain->SetBranchAddress("FatJet_pt_jesRelativePtHFDown", FatJet_pt_jesRelativePtHFDown, &b_FatJet_pt_jesRelativePtHFDown);
   fChain->SetBranchAddress("FatJet_mass_jesRelativePtHFDown", FatJet_mass_jesRelativePtHFDown, &b_FatJet_mass_jesRelativePtHFDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativePtHFDown", FatJet_msoftdrop_jesRelativePtHFDown, &b_FatJet_msoftdrop_jesRelativePtHFDown);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeBalDown", FatJet_pt_jesRelativeBalDown, &b_FatJet_pt_jesRelativeBalDown);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeBalDown", FatJet_mass_jesRelativeBalDown, &b_FatJet_mass_jesRelativeBalDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeBalDown", FatJet_msoftdrop_jesRelativeBalDown, &b_FatJet_msoftdrop_jesRelativeBalDown);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeSampleDown", FatJet_pt_jesRelativeSampleDown, &b_FatJet_pt_jesRelativeSampleDown);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeSampleDown", FatJet_mass_jesRelativeSampleDown, &b_FatJet_mass_jesRelativeSampleDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeSampleDown", FatJet_msoftdrop_jesRelativeSampleDown, &b_FatJet_msoftdrop_jesRelativeSampleDown);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeFSRDown", FatJet_pt_jesRelativeFSRDown, &b_FatJet_pt_jesRelativeFSRDown);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeFSRDown", FatJet_mass_jesRelativeFSRDown, &b_FatJet_mass_jesRelativeFSRDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeFSRDown", FatJet_msoftdrop_jesRelativeFSRDown, &b_FatJet_msoftdrop_jesRelativeFSRDown);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeStatFSRDown", FatJet_pt_jesRelativeStatFSRDown, &b_FatJet_pt_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeStatFSRDown", FatJet_mass_jesRelativeStatFSRDown, &b_FatJet_mass_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeStatFSRDown", FatJet_msoftdrop_jesRelativeStatFSRDown, &b_FatJet_msoftdrop_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeStatECDown", FatJet_pt_jesRelativeStatECDown, &b_FatJet_pt_jesRelativeStatECDown);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeStatECDown", FatJet_mass_jesRelativeStatECDown, &b_FatJet_mass_jesRelativeStatECDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeStatECDown", FatJet_msoftdrop_jesRelativeStatECDown, &b_FatJet_msoftdrop_jesRelativeStatECDown);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeStatHFDown", FatJet_pt_jesRelativeStatHFDown, &b_FatJet_pt_jesRelativeStatHFDown);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeStatHFDown", FatJet_mass_jesRelativeStatHFDown, &b_FatJet_mass_jesRelativeStatHFDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeStatHFDown", FatJet_msoftdrop_jesRelativeStatHFDown, &b_FatJet_msoftdrop_jesRelativeStatHFDown);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpDataMCDown", FatJet_pt_jesPileUpDataMCDown, &b_FatJet_pt_jesPileUpDataMCDown);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpDataMCDown", FatJet_mass_jesPileUpDataMCDown, &b_FatJet_mass_jesPileUpDataMCDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpDataMCDown", FatJet_msoftdrop_jesPileUpDataMCDown, &b_FatJet_msoftdrop_jesPileUpDataMCDown);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtRefDown", FatJet_pt_jesPileUpPtRefDown, &b_FatJet_pt_jesPileUpPtRefDown);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtRefDown", FatJet_mass_jesPileUpPtRefDown, &b_FatJet_mass_jesPileUpPtRefDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtRefDown", FatJet_msoftdrop_jesPileUpPtRefDown, &b_FatJet_msoftdrop_jesPileUpPtRefDown);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtBBDown", FatJet_pt_jesPileUpPtBBDown, &b_FatJet_pt_jesPileUpPtBBDown);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtBBDown", FatJet_mass_jesPileUpPtBBDown, &b_FatJet_mass_jesPileUpPtBBDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtBBDown", FatJet_msoftdrop_jesPileUpPtBBDown, &b_FatJet_msoftdrop_jesPileUpPtBBDown);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtEC1Down", FatJet_pt_jesPileUpPtEC1Down, &b_FatJet_pt_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtEC1Down", FatJet_mass_jesPileUpPtEC1Down, &b_FatJet_mass_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtEC1Down", FatJet_msoftdrop_jesPileUpPtEC1Down, &b_FatJet_msoftdrop_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtEC2Down", FatJet_pt_jesPileUpPtEC2Down, &b_FatJet_pt_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtEC2Down", FatJet_mass_jesPileUpPtEC2Down, &b_FatJet_mass_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtEC2Down", FatJet_msoftdrop_jesPileUpPtEC2Down, &b_FatJet_msoftdrop_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtHFDown", FatJet_pt_jesPileUpPtHFDown, &b_FatJet_pt_jesPileUpPtHFDown);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtHFDown", FatJet_mass_jesPileUpPtHFDown, &b_FatJet_mass_jesPileUpPtHFDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtHFDown", FatJet_msoftdrop_jesPileUpPtHFDown, &b_FatJet_msoftdrop_jesPileUpPtHFDown);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpMuZeroDown", FatJet_pt_jesPileUpMuZeroDown, &b_FatJet_pt_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpMuZeroDown", FatJet_mass_jesPileUpMuZeroDown, &b_FatJet_mass_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpMuZeroDown", FatJet_msoftdrop_jesPileUpMuZeroDown, &b_FatJet_msoftdrop_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpEnvelopeDown", FatJet_pt_jesPileUpEnvelopeDown, &b_FatJet_pt_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpEnvelopeDown", FatJet_mass_jesPileUpEnvelopeDown, &b_FatJet_mass_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpEnvelopeDown", FatJet_msoftdrop_jesPileUpEnvelopeDown, &b_FatJet_msoftdrop_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalPileUpDown", FatJet_pt_jesSubTotalPileUpDown, &b_FatJet_pt_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalPileUpDown", FatJet_mass_jesSubTotalPileUpDown, &b_FatJet_mass_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalPileUpDown", FatJet_msoftdrop_jesSubTotalPileUpDown, &b_FatJet_msoftdrop_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalRelativeDown", FatJet_pt_jesSubTotalRelativeDown, &b_FatJet_pt_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalRelativeDown", FatJet_mass_jesSubTotalRelativeDown, &b_FatJet_mass_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalRelativeDown", FatJet_msoftdrop_jesSubTotalRelativeDown, &b_FatJet_msoftdrop_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalPtDown", FatJet_pt_jesSubTotalPtDown, &b_FatJet_pt_jesSubTotalPtDown);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalPtDown", FatJet_mass_jesSubTotalPtDown, &b_FatJet_mass_jesSubTotalPtDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalPtDown", FatJet_msoftdrop_jesSubTotalPtDown, &b_FatJet_msoftdrop_jesSubTotalPtDown);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalScaleDown", FatJet_pt_jesSubTotalScaleDown, &b_FatJet_pt_jesSubTotalScaleDown);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalScaleDown", FatJet_mass_jesSubTotalScaleDown, &b_FatJet_mass_jesSubTotalScaleDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalScaleDown", FatJet_msoftdrop_jesSubTotalScaleDown, &b_FatJet_msoftdrop_jesSubTotalScaleDown);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalAbsoluteDown", FatJet_pt_jesSubTotalAbsoluteDown, &b_FatJet_pt_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalAbsoluteDown", FatJet_mass_jesSubTotalAbsoluteDown, &b_FatJet_mass_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalAbsoluteDown", FatJet_msoftdrop_jesSubTotalAbsoluteDown, &b_FatJet_msoftdrop_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalMCDown", FatJet_pt_jesSubTotalMCDown, &b_FatJet_pt_jesSubTotalMCDown);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalMCDown", FatJet_mass_jesSubTotalMCDown, &b_FatJet_mass_jesSubTotalMCDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalMCDown", FatJet_msoftdrop_jesSubTotalMCDown, &b_FatJet_msoftdrop_jesSubTotalMCDown);
   fChain->SetBranchAddress("FatJet_pt_jesTotalDown", FatJet_pt_jesTotalDown, &b_FatJet_pt_jesTotalDown);
   fChain->SetBranchAddress("FatJet_mass_jesTotalDown", FatJet_mass_jesTotalDown, &b_FatJet_mass_jesTotalDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTotalDown", FatJet_msoftdrop_jesTotalDown, &b_FatJet_msoftdrop_jesTotalDown);
   fChain->SetBranchAddress("FatJet_pt_jesTotalNoFlavorDown", FatJet_pt_jesTotalNoFlavorDown, &b_FatJet_pt_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("FatJet_mass_jesTotalNoFlavorDown", FatJet_mass_jesTotalNoFlavorDown, &b_FatJet_mass_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTotalNoFlavorDown", FatJet_msoftdrop_jesTotalNoFlavorDown, &b_FatJet_msoftdrop_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("FatJet_pt_jesTotalNoTimeDown", FatJet_pt_jesTotalNoTimeDown, &b_FatJet_pt_jesTotalNoTimeDown);
   fChain->SetBranchAddress("FatJet_mass_jesTotalNoTimeDown", FatJet_mass_jesTotalNoTimeDown, &b_FatJet_mass_jesTotalNoTimeDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTotalNoTimeDown", FatJet_msoftdrop_jesTotalNoTimeDown, &b_FatJet_msoftdrop_jesTotalNoTimeDown);
   fChain->SetBranchAddress("FatJet_pt_jesTotalNoFlavorNoTimeDown", FatJet_pt_jesTotalNoFlavorNoTimeDown, &b_FatJet_pt_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("FatJet_mass_jesTotalNoFlavorNoTimeDown", FatJet_mass_jesTotalNoFlavorNoTimeDown, &b_FatJet_mass_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTotalNoFlavorNoTimeDown", FatJet_msoftdrop_jesTotalNoFlavorNoTimeDown, &b_FatJet_msoftdrop_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorZJetDown", FatJet_pt_jesFlavorZJetDown, &b_FatJet_pt_jesFlavorZJetDown);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorZJetDown", FatJet_mass_jesFlavorZJetDown, &b_FatJet_mass_jesFlavorZJetDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorZJetDown", FatJet_msoftdrop_jesFlavorZJetDown, &b_FatJet_msoftdrop_jesFlavorZJetDown);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPhotonJetDown", FatJet_pt_jesFlavorPhotonJetDown, &b_FatJet_pt_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPhotonJetDown", FatJet_mass_jesFlavorPhotonJetDown, &b_FatJet_mass_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPhotonJetDown", FatJet_msoftdrop_jesFlavorPhotonJetDown, &b_FatJet_msoftdrop_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPureGluonDown", FatJet_pt_jesFlavorPureGluonDown, &b_FatJet_pt_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPureGluonDown", FatJet_mass_jesFlavorPureGluonDown, &b_FatJet_mass_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureGluonDown", FatJet_msoftdrop_jesFlavorPureGluonDown, &b_FatJet_msoftdrop_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPureQuarkDown", FatJet_pt_jesFlavorPureQuarkDown, &b_FatJet_pt_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPureQuarkDown", FatJet_mass_jesFlavorPureQuarkDown, &b_FatJet_mass_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureQuarkDown", FatJet_msoftdrop_jesFlavorPureQuarkDown, &b_FatJet_msoftdrop_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPureCharmDown", FatJet_pt_jesFlavorPureCharmDown, &b_FatJet_pt_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPureCharmDown", FatJet_mass_jesFlavorPureCharmDown, &b_FatJet_mass_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureCharmDown", FatJet_msoftdrop_jesFlavorPureCharmDown, &b_FatJet_msoftdrop_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPureBottomDown", FatJet_pt_jesFlavorPureBottomDown, &b_FatJet_pt_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPureBottomDown", FatJet_mass_jesFlavorPureBottomDown, &b_FatJet_mass_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureBottomDown", FatJet_msoftdrop_jesFlavorPureBottomDown, &b_FatJet_msoftdrop_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("FatJet_pt_jesTimeRunBDown", FatJet_pt_jesTimeRunBDown, &b_FatJet_pt_jesTimeRunBDown);
   fChain->SetBranchAddress("FatJet_mass_jesTimeRunBDown", FatJet_mass_jesTimeRunBDown, &b_FatJet_mass_jesTimeRunBDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimeRunBDown", FatJet_msoftdrop_jesTimeRunBDown, &b_FatJet_msoftdrop_jesTimeRunBDown);
   fChain->SetBranchAddress("FatJet_pt_jesTimeRunCDown", FatJet_pt_jesTimeRunCDown, &b_FatJet_pt_jesTimeRunCDown);
   fChain->SetBranchAddress("FatJet_mass_jesTimeRunCDown", FatJet_mass_jesTimeRunCDown, &b_FatJet_mass_jesTimeRunCDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimeRunCDown", FatJet_msoftdrop_jesTimeRunCDown, &b_FatJet_msoftdrop_jesTimeRunCDown);
   fChain->SetBranchAddress("FatJet_pt_jesTimeRunDEDown", FatJet_pt_jesTimeRunDEDown, &b_FatJet_pt_jesTimeRunDEDown);
   fChain->SetBranchAddress("FatJet_mass_jesTimeRunDEDown", FatJet_mass_jesTimeRunDEDown, &b_FatJet_mass_jesTimeRunDEDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimeRunDEDown", FatJet_msoftdrop_jesTimeRunDEDown, &b_FatJet_msoftdrop_jesTimeRunDEDown);
   fChain->SetBranchAddress("FatJet_pt_jesTimeRunFDown", FatJet_pt_jesTimeRunFDown, &b_FatJet_pt_jesTimeRunFDown);
   fChain->SetBranchAddress("FatJet_mass_jesTimeRunFDown", FatJet_mass_jesTimeRunFDown, &b_FatJet_mass_jesTimeRunFDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimeRunFDown", FatJet_msoftdrop_jesTimeRunFDown, &b_FatJet_msoftdrop_jesTimeRunFDown);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupMPFInSituDown", FatJet_pt_jesCorrelationGroupMPFInSituDown, &b_FatJet_pt_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupMPFInSituDown", FatJet_mass_jesCorrelationGroupMPFInSituDown, &b_FatJet_mass_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupMPFInSituDown", FatJet_msoftdrop_jesCorrelationGroupMPFInSituDown, &b_FatJet_msoftdrop_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupIntercalibrationDown", FatJet_pt_jesCorrelationGroupIntercalibrationDown, &b_FatJet_pt_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupIntercalibrationDown", FatJet_mass_jesCorrelationGroupIntercalibrationDown, &b_FatJet_mass_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupIntercalibrationDown", FatJet_msoftdrop_jesCorrelationGroupIntercalibrationDown, &b_FatJet_msoftdrop_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupbJESDown", FatJet_pt_jesCorrelationGroupbJESDown, &b_FatJet_pt_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupbJESDown", FatJet_mass_jesCorrelationGroupbJESDown, &b_FatJet_mass_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupbJESDown", FatJet_msoftdrop_jesCorrelationGroupbJESDown, &b_FatJet_msoftdrop_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupFlavorDown", FatJet_pt_jesCorrelationGroupFlavorDown, &b_FatJet_pt_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupFlavorDown", FatJet_mass_jesCorrelationGroupFlavorDown, &b_FatJet_mass_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupFlavorDown", FatJet_msoftdrop_jesCorrelationGroupFlavorDown, &b_FatJet_msoftdrop_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupUncorrelatedDown", FatJet_pt_jesCorrelationGroupUncorrelatedDown, &b_FatJet_pt_jesCorrelationGroupUncorrelatedDown);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupUncorrelatedDown", FatJet_mass_jesCorrelationGroupUncorrelatedDown, &b_FatJet_mass_jesCorrelationGroupUncorrelatedDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupUncorrelatedDown", FatJet_msoftdrop_jesCorrelationGroupUncorrelatedDown, &b_FatJet_msoftdrop_jesCorrelationGroupUncorrelatedDown);
/*
   fChain->SetBranchAddress("Jet_CSVbtagSF", Jet_CSVbtagSF, &b_Jet_CSVbtagSF);
   fChain->SetBranchAddress("Jet_CSVbtagSF_up", Jet_CSVbtagSF_up, &b_Jet_CSVbtagSF_up);
   fChain->SetBranchAddress("Jet_CSVbtagSF_down", Jet_CSVbtagSF_down, &b_Jet_CSVbtagSF_down);
   fChain->SetBranchAddress("Jet_CSVbtagSF_shape", Jet_CSVbtagSF_shape, &b_Jet_CSVbtagSF_shape);
   fChain->SetBranchAddress("Jet_CSVbtagSF_shape_up_jes", Jet_CSVbtagSF_shape_up_jes, &b_Jet_CSVbtagSF_shape_up_jes);
   fChain->SetBranchAddress("Jet_CSVbtagSF_shape_down_jes", Jet_CSVbtagSF_shape_down_jes, &b_Jet_CSVbtagSF_shape_down_jes);
   fChain->SetBranchAddress("Jet_CSVbtagSF_shape_up_lf", Jet_CSVbtagSF_shape_up_lf, &b_Jet_CSVbtagSF_shape_up_lf);
   fChain->SetBranchAddress("Jet_CSVbtagSF_shape_down_lf", Jet_CSVbtagSF_shape_down_lf, &b_Jet_CSVbtagSF_shape_down_lf);
   fChain->SetBranchAddress("Jet_CSVbtagSF_shape_up_hf", Jet_CSVbtagSF_shape_up_hf, &b_Jet_CSVbtagSF_shape_up_hf);
   fChain->SetBranchAddress("Jet_CSVbtagSF_shape_down_hf", Jet_CSVbtagSF_shape_down_hf, &b_Jet_CSVbtagSF_shape_down_hf);
   fChain->SetBranchAddress("Jet_CSVbtagSF_shape_up_hfstats1", Jet_CSVbtagSF_shape_up_hfstats1, &b_Jet_CSVbtagSF_shape_up_hfstats1);
   fChain->SetBranchAddress("Jet_CSVbtagSF_shape_down_hfstats1", Jet_CSVbtagSF_shape_down_hfstats1, &b_Jet_CSVbtagSF_shape_down_hfstats1);
   fChain->SetBranchAddress("Jet_CSVbtagSF_shape_up_hfstats2", Jet_CSVbtagSF_shape_up_hfstats2, &b_Jet_CSVbtagSF_shape_up_hfstats2);
   fChain->SetBranchAddress("Jet_CSVbtagSF_shape_down_hfstats2", Jet_CSVbtagSF_shape_down_hfstats2, &b_Jet_CSVbtagSF_shape_down_hfstats2);
   fChain->SetBranchAddress("Jet_CSVbtagSF_shape_up_lfstats1", Jet_CSVbtagSF_shape_up_lfstats1, &b_Jet_CSVbtagSF_shape_up_lfstats1);
   fChain->SetBranchAddress("Jet_CSVbtagSF_shape_down_lfstats1", Jet_CSVbtagSF_shape_down_lfstats1, &b_Jet_CSVbtagSF_shape_down_lfstats1);
   fChain->SetBranchAddress("Jet_CSVbtagSF_shape_up_lfstats2", Jet_CSVbtagSF_shape_up_lfstats2, &b_Jet_CSVbtagSF_shape_up_lfstats2);
   fChain->SetBranchAddress("Jet_CSVbtagSF_shape_down_lfstats2", Jet_CSVbtagSF_shape_down_lfstats2, &b_Jet_CSVbtagSF_shape_down_lfstats2);
   fChain->SetBranchAddress("Jet_CSVbtagSF_shape_up_cferr1", Jet_CSVbtagSF_shape_up_cferr1, &b_Jet_CSVbtagSF_shape_up_cferr1);
   fChain->SetBranchAddress("Jet_CSVbtagSF_shape_down_cferr1", Jet_CSVbtagSF_shape_down_cferr1, &b_Jet_CSVbtagSF_shape_down_cferr1);
   fChain->SetBranchAddress("Jet_CSVbtagSF_shape_up_cferr2", Jet_CSVbtagSF_shape_up_cferr2, &b_Jet_CSVbtagSF_shape_up_cferr2);
   fChain->SetBranchAddress("Jet_CSVbtagSF_shape_down_cferr2", Jet_CSVbtagSF_shape_down_cferr2, &b_Jet_CSVbtagSF_shape_down_cferr2);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF", Jet_deepCSVbtagSF, &b_Jet_deepCSVbtagSF);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_up", Jet_deepCSVbtagSF_up, &b_Jet_deepCSVbtagSF_up);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_down", Jet_deepCSVbtagSF_down, &b_Jet_deepCSVbtagSF_down);  
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_shape", Jet_deepCSVbtagSF_shape, &b_Jet_deepCSVbtagSF_shape);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_shape_up_jes", Jet_deepCSVbtagSF_shape_up_jes, &b_Jet_deepCSVbtagSF_shape_up_jes);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_shape_down_jes", Jet_deepCSVbtagSF_shape_down_jes, &b_Jet_deepCSVbtagSF_shape_down_jes);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_shape_up_lf", Jet_deepCSVbtagSF_shape_up_lf, &b_Jet_deepCSVbtagSF_shape_up_lf);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_shape_down_lf", Jet_deepCSVbtagSF_shape_down_lf, &b_Jet_deepCSVbtagSF_shape_down_lf);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_shape_up_hf", Jet_deepCSVbtagSF_shape_up_hf, &b_Jet_deepCSVbtagSF_shape_up_hf);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_shape_down_hf", Jet_deepCSVbtagSF_shape_down_hf, &b_Jet_deepCSVbtagSF_shape_down_hf);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_shape_up_hfstats1", Jet_deepCSVbtagSF_shape_up_hfstats1, &b_Jet_deepCSVbtagSF_shape_up_hfstats1);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_shape_down_hfstats1", Jet_deepCSVbtagSF_shape_down_hfstats1, &b_Jet_deepCSVbtagSF_shape_down_hfstats1);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_shape_up_hfstats2", Jet_deepCSVbtagSF_shape_up_hfstats2, &b_Jet_deepCSVbtagSF_shape_up_hfstats2);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_shape_down_hfstats2", Jet_deepCSVbtagSF_shape_down_hfstats2, &b_Jet_deepCSVbtagSF_shape_down_hfstats2);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_shape_up_lfstats1", Jet_deepCSVbtagSF_shape_up_lfstats1, &b_Jet_deepCSVbtagSF_shape_up_lfstats1);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_shape_down_lfstats1", Jet_deepCSVbtagSF_shape_down_lfstats1, &b_Jet_deepCSVbtagSF_shape_down_lfstats1);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_shape_up_lfstats2", Jet_deepCSVbtagSF_shape_up_lfstats2, &b_Jet_deepCSVbtagSF_shape_up_lfstats2);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_shape_down_lfstats2", Jet_deepCSVbtagSF_shape_down_lfstats2, &b_Jet_deepCSVbtagSF_shape_down_lfstats2);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_shape_up_cferr1", Jet_deepCSVbtagSF_shape_up_cferr1, &b_Jet_deepCSVbtagSF_shape_up_cferr1);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_shape_down_cferr1", Jet_deepCSVbtagSF_shape_down_cferr1, &b_Jet_deepCSVbtagSF_shape_down_cferr1);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_shape_up_cferr2", Jet_deepCSVbtagSF_shape_up_cferr2, &b_Jet_deepCSVbtagSF_shape_up_cferr2);
   fChain->SetBranchAddress("Jet_deepCSVbtagSF_shape_down_cferr2", Jet_deepCSVbtagSF_shape_down_cferr2, &b_Jet_deepCSVbtagSF_shape_down_cferr2);
*/
   fChain->SetBranchAddress("Jet_deepflavbtagSF", Jet_deepflavbtagSF, &b_Jet_deepflavbtagSF);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_up", Jet_deepflavbtagSF_up, &b_Jet_deepflavbtagSF_up);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_down", Jet_deepflavbtagSF_down, &b_Jet_deepflavbtagSF_down);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_shape", Jet_deepflavbtagSF_shape, &b_Jet_deepflavbtagSF_shape);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_shape_up_jes", Jet_deepflavbtagSF_shape_up_jes, &b_Jet_deepflavbtagSF_shape_up_jes);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_shape_down_jes", Jet_deepflavbtagSF_shape_down_jes, &b_Jet_deepflavbtagSF_shape_down_jes);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_shape_up_lf", Jet_deepflavbtagSF_shape_up_lf, &b_Jet_deepflavbtagSF_shape_up_lf);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_shape_down_lf", Jet_deepflavbtagSF_shape_down_lf, &b_Jet_deepflavbtagSF_shape_down_lf);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_shape_up_hf", Jet_deepflavbtagSF_shape_up_hf, &b_Jet_deepflavbtagSF_shape_up_hf);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_shape_down_hf", Jet_deepflavbtagSF_shape_down_hf, &b_Jet_deepflavbtagSF_shape_down_hf);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_shape_up_hfstats1", Jet_deepflavbtagSF_shape_up_hfstats1, &b_Jet_deepflavbtagSF_shape_up_hfstats1);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_shape_down_hfstats1", Jet_deepflavbtagSF_shape_down_hfstats1, &b_Jet_deepflavbtagSF_shape_down_hfstats1);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_shape_up_hfstats2", Jet_deepflavbtagSF_shape_up_hfstats2, &b_Jet_deepflavbtagSF_shape_up_hfstats2);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_shape_down_hfstats2", Jet_deepflavbtagSF_shape_down_hfstats2, &b_Jet_deepflavbtagSF_shape_down_hfstats2);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_shape_up_lfstats1", Jet_deepflavbtagSF_shape_up_lfstats1, &b_Jet_deepflavbtagSF_shape_up_lfstats1);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_shape_down_lfstats1", Jet_deepflavbtagSF_shape_down_lfstats1, &b_Jet_deepflavbtagSF_shape_down_lfstats1);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_shape_up_lfstats2", Jet_deepflavbtagSF_shape_up_lfstats2, &b_Jet_deepflavbtagSF_shape_up_lfstats2);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_shape_down_lfstats2", Jet_deepflavbtagSF_shape_down_lfstats2, &b_Jet_deepflavbtagSF_shape_down_lfstats2);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_shape_up_cferr1", Jet_deepflavbtagSF_shape_up_cferr1, &b_Jet_deepflavbtagSF_shape_up_cferr1);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_shape_down_cferr1", Jet_deepflavbtagSF_shape_down_cferr1, &b_Jet_deepflavbtagSF_shape_down_cferr1);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_shape_up_cferr2", Jet_deepflavbtagSF_shape_up_cferr2, &b_Jet_deepflavbtagSF_shape_up_cferr2);
   fChain->SetBranchAddress("Jet_deepflavbtagSF_shape_down_cferr2", Jet_deepflavbtagSF_shape_down_cferr2, &b_Jet_deepflavbtagSF_shape_down_cferr2);
   fChain->SetBranchAddress("puWeight", &puWeight, &b_puWeight);
   fChain->SetBranchAddress("puWeightUp", &puWeightUp, &b_puWeightUp);
   fChain->SetBranchAddress("puWeightDown", &puWeightDown, &b_puWeightDown);
   fChain->SetBranchAddress("PrefireWeight", &PrefireWeight, &b_PrefireWeight);
   fChain->SetBranchAddress("PrefireWeight_Up", &PrefireWeight_Up, &b_PrefireWeight_Up);
   fChain->SetBranchAddress("PrefireWeight_Down", &PrefireWeight_Down, &b_PrefireWeight_Down);

}

Bool_t Anal_Nano_PROOF::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef Anal_Nano_PROOF_cxx
