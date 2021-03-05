//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Aug 17 22:35:01 2020 by ROOT version 5.34/36
// from TTree Events/Events
// found on file: root://se01.indiacms.res.in//store/user/chatterj/NanoPost/WprimeToTB_TToHad_M-3300_RH_TuneCP5_13TeV-comphep-pythia8/NanoPost_2017v6_Wp_RH_THad_M3300/200603_105457/0000/tree_1.root
//////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

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

#include "LHAPDF/LHAPDF.h"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

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

double Pol2(double* x, double* par){
return (par[0]+par[1]*x[0]+par[2]*x[0]*x[0]);
}

double bdF0_pol2(double* x, double* par){
return 1;
}
double bdF1_pol2(double* x, double* par){
return x[0];
}
double bdF2_pol2(double* x, double* par){
return (x[0]*x[0]);
}

double SF_TOP(double alpha, double beta, double pt0, double pt1)
{
	double sfwt = sqrt(exp(alpha-beta*pt0) * exp(alpha-beta*pt1));
	return sfwt;
}


   static const int njetmax = 25;
   static const int njesmax = 37;
   static const int nbtagmax = 9;
   static const int npartmax = 200;
   static const int npdfmax = 45;//33;
   static const int nscalemax = 9;
   static const int npsmax = 4;
   static const int npartmx = 100;
   
   float weight;
   float weight_t;
   float SF_toppt;
   double weightev, weightpass, weightpass_btag;
   double weight_gen_pdf[npdfmax];
   double weight_gen_scale[nscalemax];
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
   bool TopTagging;
   
   char name[1000];
   char title[1000];
 
   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   Float_t         HTXS_Higgs_pt;
   Float_t         HTXS_Higgs_y;
   Int_t           HTXS_stage1_1_cat_pTjet25GeV;
   Int_t           HTXS_stage1_1_cat_pTjet30GeV;
   Int_t           HTXS_stage1_1_fine_cat_pTjet25GeV;
   Int_t           HTXS_stage1_1_fine_cat_pTjet30GeV;
   Int_t           HTXS_stage_0;
   Int_t           HTXS_stage_1_pTjet25;
   Int_t           HTXS_stage_1_pTjet30;
   UChar_t         HTXS_njets25;
   UChar_t         HTXS_njets30;
   Float_t         btagWeight_CSVV2;
   Float_t         btagWeight_DeepCSVB;
   Float_t         CaloMET_phi;
   Float_t         CaloMET_pt;
   Float_t         CaloMET_sumEt;
   Float_t         ChsMET_phi;
   Float_t         ChsMET_pt;
   Float_t         ChsMET_sumEt;
   UInt_t          nCorrT1METJet;
   Float_t         CorrT1METJet_area[njetmax];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_eta[njetmax];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_muonSubtrFactor[njetmax];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_phi[njetmax];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_rawPt[njetmax];   //[nCorrT1METJet]
   UInt_t          nElectron;
   Float_t         Electron_deltaEtaSC[njetmax];
   Float_t         Electron_dr03EcalRecHitSumEt[njetmax];
   Float_t         Electron_dr03HcalDepth1TowerSumEt[njetmax];
   Float_t         Electron_dr03TkSumPt[njetmax];
   Float_t         Electron_dr03TkSumPtHEEP[njetmax];
   Float_t         Electron_dxy[njetmax];
   Float_t         Electron_dxyErr[njetmax];
   Float_t         Electron_dz[njetmax];
   Float_t         Electron_dzErr[njetmax];
   Float_t         Electron_eCorr[njetmax];
   Float_t         Electron_eInvMinusPInv[njetmax];
   Float_t         Electron_energyErr[njetmax];
   Float_t         Electron_eta[njetmax];
   Float_t         Electron_hoe[njetmax];
   Float_t         Electron_ip3d[njetmax];
   Float_t         Electron_jetPtRelv2[njetmax];
   Float_t         Electron_jetRelIso[njetmax];
   Float_t         Electron_mass[njetmax];
   Float_t         Electron_miniPFRelIso_all[njetmax];
   Float_t         Electron_miniPFRelIso_chg[njetmax];
   Float_t         Electron_mvaFall17V1Iso[njetmax];
   Float_t         Electron_mvaFall17V1noIso[njetmax];
   Float_t         Electron_mvaFall17V2Iso[njetmax];
   Float_t         Electron_mvaFall17V2noIso[njetmax];
   Float_t         Electron_pfRelIso03_all[njetmax];
   Float_t         Electron_pfRelIso03_chg[njetmax];
   Float_t         Electron_phi[njetmax];
   Float_t         Electron_pt[njetmax];
   Float_t         Electron_r9[njetmax];
   Float_t         Electron_sieie[njetmax];
   Float_t         Electron_sip3d[njetmax];
   Float_t         Electron_mvaTTH[njetmax];
   Int_t           Electron_charge[njetmax];
   Int_t           Electron_cutBased[njetmax];
   Int_t           Electron_cutBased_Fall17_V1[njetmax];
   Int_t           Electron_jetIdx[njetmax];
   Int_t           Electron_pdgId[njetmax];
   Int_t           Electron_photonIdx[njetmax];
   Int_t           Electron_tightCharge[njetmax];
   Int_t           Electron_vidNestedWPBitmap[njetmax];
   Int_t           Electron_vidNestedWPBitmapHEEP[njetmax];
   Bool_t          Electron_convVeto[njetmax];
   Bool_t          Electron_cutBased_HEEP[njetmax];
   Bool_t          Electron_isPFcand[njetmax];
   UChar_t         Electron_lostHits[njetmax];
   Bool_t          Electron_mvaFall17V1Iso_WP80[njetmax];
   Bool_t          Electron_mvaFall17V1Iso_WP90[njetmax];
   Bool_t          Electron_mvaFall17V1Iso_WPL[njetmax];
   Bool_t          Electron_mvaFall17V1noIso_WP80[njetmax];
   Bool_t          Electron_mvaFall17V1noIso_WP90[njetmax];
   Bool_t          Electron_mvaFall17V1noIso_WPL[njetmax];
   Bool_t          Electron_mvaFall17V2Iso_WP80[njetmax];
   Bool_t          Electron_mvaFall17V2Iso_WP90[njetmax];
   Bool_t          Electron_mvaFall17V2Iso_WPL[njetmax];
   Bool_t          Electron_mvaFall17V2noIso_WP80[njetmax];
   Bool_t          Electron_mvaFall17V2noIso_WP90[njetmax];
   Bool_t          Electron_mvaFall17V2noIso_WPL[njetmax];
   UChar_t         Electron_seedGain[njetmax];
   Bool_t          Flag_ecalBadCalibFilterV2;
   UInt_t          nFatJet;
   Float_t         FatJet_area[njetmax];
   Float_t         FatJet_btagCMVA[njetmax];
   Float_t         FatJet_btagCSVV2[njetmax];
   Float_t         FatJet_btagDDBvL[njetmax];
   Float_t         FatJet_btagDDBvL_noMD[njetmax];
   Float_t         FatJet_btagDDCvB[njetmax];
   Float_t         FatJet_btagDDCvB_noMD[njetmax];
   Float_t         FatJet_btagDDCvL[njetmax];
   Float_t         FatJet_btagDDCvL_noMD[njetmax];
   Float_t         FatJet_btagDeepB[njetmax];
   Float_t         FatJet_btagHbb[njetmax];
   Float_t         FatJet_deepTagMD_H4qvsQCD[njetmax];
   Float_t         FatJet_deepTagMD_HbbvsQCD[njetmax];
   Float_t         FatJet_deepTagMD_TvsQCD[njetmax];
   Float_t         FatJet_deepTagMD_WvsQCD[njetmax];
   Float_t         FatJet_deepTagMD_ZHbbvsQCD[njetmax];
   Float_t         FatJet_deepTagMD_ZHccvsQCD[njetmax];
   Float_t         FatJet_deepTagMD_ZbbvsQCD[njetmax];
   Float_t         FatJet_deepTagMD_ZvsQCD[njetmax];
   Float_t         FatJet_deepTagMD_bbvsLight[njetmax];
   Float_t         FatJet_deepTagMD_ccvsLight[njetmax];
   Float_t         FatJet_deepTag_H[njetmax];
   Float_t         FatJet_deepTag_QCD[njetmax];
   Float_t         FatJet_deepTag_QCDothers[njetmax];
   Float_t         FatJet_deepTag_TvsQCD[njetmax];
   Float_t         FatJet_deepTag_WvsQCD[njetmax];
   Float_t         FatJet_deepTag_ZvsQCD[njetmax];
   Float_t         FatJet_eta[njetmax];
   Float_t         FatJet_mass[njetmax];
   Float_t         FatJet_msoftdrop[njetmax];
   Float_t         FatJet_n2b1[njetmax];
   Float_t         FatJet_n3b1[njetmax];
   Float_t         FatJet_phi[njetmax];
   Float_t         FatJet_pt[njetmax];
   Float_t         FatJet_rawFactor[njetmax];
   Float_t         FatJet_tau1[njetmax];
   Float_t         FatJet_tau2[njetmax];
   Float_t         FatJet_tau3[njetmax];
   Float_t         FatJet_tau4[njetmax];
   Float_t		   FatJet_subbtagCSVV2[njetmax];
   Float_t		   FatJet_subbtagDeepB[njetmax];
   Int_t           FatJet_jetId[njetmax];
   Int_t           FatJet_subJetIdx1[njetmax];
   Int_t           FatJet_subJetIdx2[njetmax];
   Int_t 		   FatJet_subbtagId_CSVv2[njetmax];  
   Int_t 		   FatJet_subbtagId_DeepCSV[njetmax];  
   UInt_t          nFsrPhoton;
   Float_t         FsrPhoton_dROverEt2[njetmax];
   Float_t         FsrPhoton_eta[njetmax];
   Float_t         FsrPhoton_phi[njetmax];
   Float_t         FsrPhoton_pt[njetmax];
   Float_t         FsrPhoton_relIso03[njetmax];
   Int_t           FsrPhoton_muonIdx[njetmax];
   UInt_t          nGenJetAK8;
   Float_t         GenJetAK8_eta[njetmax];
   Float_t         GenJetAK8_mass[njetmax];
   Float_t         GenJetAK8_phi[njetmax];
   Float_t         GenJetAK8_pt[njetmax];
   UInt_t          nGenJet;
   Float_t         GenJet_eta[njetmax];   //[nGenJet]
   Float_t         GenJet_mass[njetmax];   //[nGenJet]
   Float_t         GenJet_phi[njetmax];   //[nGenJet]
   Float_t         GenJet_pt[njetmax];   //[nGenJet]
   UInt_t          nGenPart;
   Float_t         GenPart_eta[npartmax];
   Float_t         GenPart_mass[npartmax];
   Float_t         GenPart_phi[npartmax];
   Float_t         GenPart_pt[npartmax];
   Int_t           GenPart_genPartIdxMother[npartmax];
   Int_t           GenPart_pdgId[npartmax];
   Int_t           GenPart_status[npartmax];
   Int_t           GenPart_statusFlags[npartmax];
   UInt_t          nSubGenJetAK8;
   Float_t         SubGenJetAK8_eta[njetmax];
   Float_t         SubGenJetAK8_mass[njetmax];
   Float_t         SubGenJetAK8_phi[njetmax];
   Float_t         SubGenJetAK8_pt[njetmax];
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
   Float_t         GenVisTau_eta[npartmax];
   Float_t         GenVisTau_mass[npartmax];
   Float_t         GenVisTau_phi[npartmax];
   Float_t         GenVisTau_pt[npartmax];
   Int_t           GenVisTau_charge[npartmax];
   Int_t           GenVisTau_genPartIdxMother[npartmax];
   Int_t           GenVisTau_status[npartmax];
   Float_t         genWeight;
   Float_t         LHEWeight_originalXWGTUP;
   UInt_t          nLHEPdfWeight;
   Float_t         LHEPdfWeight[npdfmax];   //[nLHEPdfWeight]
   UInt_t          nLHEReweightingWeight;
   Float_t         LHEReweightingWeight[npdfmax];   //[nLHEReweightingWeight]
   UInt_t          nLHEScaleWeight;
   Float_t         LHEScaleWeight[nscalemax];   //[nLHEScaleWeight]
   UInt_t          nPSWeight;
   Float_t         PSWeight[npsmax];   //[nPSWeight]
   UInt_t          nIsoTrack;
   Float_t         IsoTrack_dxy[njetmax];
   Float_t         IsoTrack_dz[njetmax];
   Float_t         IsoTrack_eta[njetmax];
   Float_t         IsoTrack_pfRelIso03_all[njetmax];
   Float_t         IsoTrack_pfRelIso03_chg[njetmax];
   Float_t         IsoTrack_phi[njetmax];
   Float_t         IsoTrack_pt[njetmax];
   Float_t         IsoTrack_miniPFRelIso_all[njetmax];
   Float_t         IsoTrack_miniPFRelIso_chg[njetmax];
   Int_t           IsoTrack_fromPV[njetmax];
   Int_t           IsoTrack_pdgId[njetmax];
   Bool_t          IsoTrack_isHighPurityTrack[njetmax];
   Bool_t          IsoTrack_isPFcand[njetmax];
   Bool_t          IsoTrack_isFromLostTrack[njetmax];
   UInt_t          nJet;
   Float_t         Jet_area[njetmax];
   Float_t         Jet_btagCMVA[njetmax];
   Float_t         Jet_btagCSVV2[njetmax];
   Float_t         Jet_btagDeepB[njetmax];
   Float_t         Jet_btagDeepC[njetmax];
   Float_t         Jet_btagDeepFlavB[njetmax];
   Float_t         Jet_btagDeepFlavC[njetmax];
   Float_t         Jet_chEmEF[njetmax];
   Float_t         Jet_chHEF[njetmax];
   Float_t         Jet_eta[njetmax];
   Float_t         Jet_jercCHF[njetmax];
   Float_t         Jet_jercCHPUF[njetmax];
   Float_t         Jet_mass[njetmax];
   Float_t         Jet_muEF[njetmax];
   Float_t         Jet_muonSubtrFactor[njetmax];
   Float_t         Jet_neEmEF[njetmax];
   Float_t         Jet_neHEF[njetmax];
   Float_t         Jet_phi[njetmax];
   Float_t         Jet_pt[njetmax];
   Float_t         Jet_qgl[njetmax];
   Float_t         Jet_rawFactor[njetmax];
   Float_t         Jet_bRegCorr[njetmax];
   Float_t         Jet_bRegRes[njetmax];
   Int_t           Jet_electronIdx1[njetmax];
   Int_t           Jet_electronIdx2[njetmax];
   Int_t           Jet_jetId[njetmax];
   Int_t 		   Jet_MatchFatJet[njetmax];
   Int_t           Jet_muonIdx1[njetmax];
   Int_t           Jet_muonIdx2[njetmax];
   Int_t           Jet_nConstituents[njetmax];
   Int_t           Jet_nElectrons[njetmax];
   Int_t           Jet_nMuons[njetmax];
   Int_t           Jet_puId[njetmax];
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
   Float_t         L1PreFiringWeight_Dn;
   Float_t         L1PreFiringWeight_Nom;
   Float_t         L1PreFiringWeight_Up;
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
   Float_t         LHEPart_pt[npartmax];
   Float_t         LHEPart_eta[npartmax];
   Float_t         LHEPart_phi[npartmax];
   Float_t         LHEPart_mass[npartmax];
   Int_t           LHEPart_pdgId[npartmax];
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
   Float_t         Muon_dxy[njetmax];
   Float_t         Muon_dxyErr[njetmax];
   Float_t         Muon_dz[njetmax];
   Float_t         Muon_dzErr[njetmax];
   Float_t         Muon_eta[njetmax];
   Float_t         Muon_ip3d[njetmax];
   Float_t         Muon_jetPtRelv2[njetmax];
   Float_t         Muon_jetRelIso[njetmax];
   Float_t         Muon_mass[njetmax];
   Float_t         Muon_miniPFRelIso_all[njetmax];
   Float_t         Muon_miniPFRelIso_chg[njetmax];
   Float_t         Muon_pfRelIso03_all[njetmax];
   Float_t         Muon_pfRelIso03_chg[njetmax];
   Float_t         Muon_pfRelIso04_all[njetmax];
   Float_t         Muon_phi[njetmax];
   Float_t         Muon_pt[njetmax];
   Float_t         Muon_ptErr[njetmax];
   Float_t         Muon_segmentComp[njetmax];
   Float_t         Muon_sip3d[njetmax];
   Float_t         Muon_tkRelIso[njetmax];
   Float_t         Muon_tunepRelPt[njetmax];
   Float_t         Muon_mvaLowPt[njetmax];
   Float_t         Muon_mvaTTH[njetmax];
   Int_t           Muon_charge[njetmax];
   Int_t           Muon_jetIdx[njetmax];
   Int_t           Muon_nStations[njetmax];
   Int_t           Muon_nTrackerLayers[njetmax];
   Int_t           Muon_pdgId[njetmax];
   Int_t           Muon_tightCharge[njetmax];
   Int_t           Muon_fsrPhotonIdx[njetmax];
   UChar_t         Muon_highPtId[njetmax];
   Bool_t          Muon_inTimeMuon[njetmax];
   Bool_t          Muon_isGlobal[njetmax];
   Bool_t          Muon_isPFcand[njetmax];
   Bool_t          Muon_isTracker[njetmax];
   Bool_t          Muon_looseId[njetmax];
   Bool_t          Muon_mediumId[njetmax];
   Bool_t          Muon_mediumPromptId[njetmax];
   UChar_t         Muon_miniIsoId[njetmax];
   UChar_t         Muon_multiIsoId[njetmax];
   UChar_t         Muon_mvaId[njetmax];
   UChar_t         Muon_pfIsoId[njetmax];
   Bool_t          Muon_softId[njetmax];
   Bool_t          Muon_softMvaId[njetmax];
   Bool_t          Muon_tightId[njetmax];
   UChar_t         Muon_tkIsoId[njetmax];
   Bool_t          Muon_triggerIdLoose[njetmax];
   UInt_t          nPhoton;
   Float_t         Photon_eCorr[njetmax];
   Float_t         Photon_energyErr[njetmax];
   Float_t         Photon_eta[njetmax];
   Float_t         Photon_hoe[njetmax];
   Float_t         Photon_mass[njetmax];
   Float_t         Photon_mvaID[njetmax];
   Float_t         Photon_mvaIDV1[njetmax];
   Float_t         Photon_pfRelIso03_all[njetmax];
   Float_t         Photon_pfRelIso03_chg[njetmax];
   Float_t         Photon_phi[njetmax];
   Float_t         Photon_pt[njetmax];
   Float_t         Photon_r9[njetmax];
   Float_t         Photon_sieie[njetmax];
   Int_t           Photon_charge[njetmax];
   Int_t           Photon_cutBasedBitmap[njetmax];
   Int_t           Photon_cutBasedV1Bitmap[njetmax];
   Int_t           Photon_electronIdx[njetmax];
   Int_t           Photon_jetIdx[njetmax];
   Int_t           Photon_pdgId[njetmax];
   Int_t           Photon_vidNestedWPBitmap[njetmax];
   Bool_t          Photon_electronVeto[njetmax];
   Bool_t          Photon_isScEtaEB[njetmax];
   Bool_t          Photon_isScEtaEE[njetmax];
   Bool_t          Photon_mvaID_WP80[njetmax];
   Bool_t          Photon_mvaID_WP90[njetmax];
   Bool_t          Photon_pixelSeed[njetmax];
   UChar_t         Photon_seedGain[njetmax];
   Float_t         Pileup_nTrueInt;
   Float_t         Pileup_pudensity;
   Float_t         Pileup_gpudensity;
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
   Float_t         fixedGridRhoFastjetCentral;
   Float_t         fixedGridRhoFastjetCentralCalo;
   Float_t         fixedGridRhoFastjetCentralChargedPileUp;
   Float_t         fixedGridRhoFastjetCentralNeutral;
   UInt_t          nGenDressedLepton;
   Float_t         GenDressedLepton_eta[2];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_mass[2];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_phi[2];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_pt[2];   //[nGenDressedLepton]
   Int_t           GenDressedLepton_pdgId[2];   //[nGenDressedLepton]
   Bool_t          GenDressedLepton_hasTauAnc[2];   //[nGenDressedLepton]
   UInt_t          nSoftActivityJet;
   Float_t         SoftActivityJet_eta[6];   //[nSoftActivityJet]
   Float_t         SoftActivityJet_phi[6];   //[nSoftActivityJet]
   Float_t         SoftActivityJet_pt[6];   //[nSoftActivityJet]
   Float_t         SoftActivityJetHT;
   Float_t         SoftActivityJetHT10;
   Float_t         SoftActivityJetHT2;
   Float_t         SoftActivityJetHT5;
   Int_t           SoftActivityJetNjets10;
   Int_t           SoftActivityJetNjets2;
   Int_t           SoftActivityJetNjets5;
   UInt_t          nSubJet;
   Float_t         SubJet_btagCMVA[njetmax];
   Float_t         SubJet_btagCSVV2[njetmax];
   Float_t         SubJet_btagDeepB[njetmax];
   Float_t         SubJet_eta[njetmax];
   Float_t         SubJet_mass[njetmax];
   Float_t         SubJet_n2b1[njetmax];
   Float_t         SubJet_n3b1[njetmax];
   Float_t         SubJet_phi[njetmax];
   Float_t         SubJet_pt[njetmax];
   Float_t         SubJet_rawFactor[njetmax];
   Float_t         SubJet_tau1[njetmax];
   Float_t         SubJet_tau2[njetmax];
   Float_t         SubJet_tau3[njetmax];
   Float_t         SubJet_tau4[njetmax];
   UInt_t          nTau;
   Float_t         Tau_chargedIso[njetmax];
   Float_t         Tau_dxy[njetmax];
   Float_t         Tau_dz[njetmax];
   Float_t         Tau_eta[njetmax];
   Float_t         Tau_leadTkDeltaEta[njetmax];
   Float_t         Tau_leadTkDeltaPhi[njetmax];
   Float_t         Tau_leadTkPtOverTauPt[njetmax];
   Float_t         Tau_mass[njetmax];
   Float_t         Tau_neutralIso[njetmax];
   Float_t         Tau_phi[njetmax];
   Float_t         Tau_photonsOutsideSignalCone[njetmax];
   Float_t         Tau_pt[njetmax];
   Float_t         Tau_puCorr[njetmax];
   Float_t         Tau_rawAntiEle[njetmax];
   Float_t         Tau_rawAntiEle2018[njetmax];
   Float_t         Tau_rawDeepTau2017v2p1VSe[njetmax];
   Float_t         Tau_rawDeepTau2017v2p1VSjet[njetmax];
   Float_t         Tau_rawDeepTau2017v2p1VSmu[njetmax];
   Float_t         Tau_rawIso[njetmax];
   Float_t         Tau_rawIsodR03[njetmax];
   Float_t         Tau_rawMVAnewDM2017v2[njetmax];
   Float_t         Tau_rawMVAoldDM[njetmax];
   Float_t         Tau_rawMVAoldDM2017v1[njetmax];
   Float_t         Tau_rawMVAoldDM2017v2[njetmax];
   Float_t         Tau_rawMVAoldDMdR032017v2[njetmax];
   Int_t           Tau_charge[njetmax];
   Int_t           Tau_decayMode[njetmax];
   Int_t           Tau_jetIdx[njetmax];
   Int_t           Tau_rawAntiEleCat[njetmax];
   Int_t           Tau_rawAntiEleCat2018[njetmax];
   UChar_t         Tau_idAntiEle[njetmax];
   UChar_t         Tau_idAntiEle2018[njetmax];
   UChar_t         Tau_idAntiMu[njetmax];
   Bool_t          Tau_idDecayMode[njetmax];
   Bool_t          Tau_idDecayModeNewDMs[njetmax];
   UChar_t         Tau_idDeepTau2017v2p1VSe[njetmax];
   UChar_t         Tau_idDeepTau2017v2p1VSjet[njetmax];
   UChar_t         Tau_idDeepTau2017v2p1VSmu[njetmax];
   UChar_t         Tau_idMVAnewDM2017v2[njetmax];
   UChar_t         Tau_idMVAoldDM[njetmax];
   UChar_t         Tau_idMVAoldDM2017v1[njetmax];
   UChar_t         Tau_idMVAoldDM2017v2[njetmax];
   UChar_t         Tau_idMVAoldDMdR032017v2[njetmax];
   Float_t         TkMET_phi;
   Float_t         TkMET_pt;
   Float_t         TkMET_sumEt;
   UInt_t          nTrigObj;
   Float_t         TrigObj_pt[njetmax];
   Float_t         TrigObj_eta[njetmax];
   Float_t         TrigObj_phi[njetmax];
   Float_t         TrigObj_l1pt[njetmax];
   Float_t         TrigObj_l1pt_2[njetmax];
   Float_t         TrigObj_l2pt[njetmax];
   Int_t           TrigObj_id[njetmax];
   Int_t           TrigObj_l1iso[njetmax];
   Int_t           TrigObj_l1charge[njetmax];
   Int_t           TrigObj_filterBits[njetmax];
   Int_t           genTtbarId;
   UInt_t          nOtherPV;
   Float_t         OtherPV_z[njetmax];   //[nOtherPV]
   Float_t         PV_ndof;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_chi2;
   Float_t         PV_score;
   Int_t           PV_npvs;
   Int_t           PV_npvsGood;
   UInt_t          nSV;
   Float_t         SV_dlen[njetmax];
   Float_t         SV_dlenSig[njetmax];
   Float_t         SV_dxy[njetmax];
   Float_t         SV_dxySig[njetmax];
   Float_t         SV_pAngle[njetmax];
   Int_t           Electron_genPartIdx[njetmax];
   UChar_t         Electron_genPartFlav[njetmax];
   Int_t           GenJetAK8_partonFlavour[njetmax];
   UChar_t         GenJetAK8_hadronFlavour[njetmax];
   Int_t           GenJet_partonFlavour[21];   //[nGenJet]
   UChar_t         GenJet_hadronFlavour[21];   //[nGenJet]
   Int_t           Jet_genJetIdx[njetmax];
   Int_t           Jet_hadronFlavour[njetmax];
   Int_t           Jet_partonFlavour[njetmax];
   Int_t           Muon_genPartIdx[njetmax];
   UChar_t         Muon_genPartFlav[njetmax];
   Int_t           Photon_genPartIdx[njetmax];
   UChar_t         Photon_genPartFlav[njetmax];
   Float_t         MET_fiducialGenPhi;
   Float_t         MET_fiducialGenPt;
   UChar_t         Electron_cleanmask[njetmax];
   UChar_t         Jet_cleanmask[njetmax];
   UChar_t         Muon_cleanmask[njetmax];
   UChar_t         Photon_cleanmask[njetmax];
   UChar_t         Tau_cleanmask[njetmax];
   Float_t         SV_chi2[njetmax];
   Float_t         SV_eta[njetmax];
   Float_t         SV_mass[njetmax];
   Float_t         SV_ndof[njetmax];
   Float_t         SV_phi[njetmax];
   Float_t         SV_pt[njetmax];
   Float_t         SV_x[njetmax];
   Float_t         SV_y[njetmax];
   Float_t         SV_z[njetmax];
   Int_t           Tau_genPartIdx[njetmax];
   UChar_t         Tau_genPartFlav[njetmax];
   Bool_t          L1simulation_step;
   Bool_t          HLTriggerFirstPath;
   Bool_t          HLT_AK8PFJet360_TrimMass30;
   Bool_t          HLT_AK8PFJet380_TrimMass30;
   Bool_t          HLT_AK8PFJet400_TrimMass30;
   Bool_t          HLT_AK8PFJet420_TrimMass30;
   Bool_t          HLT_AK8PFHT750_TrimMass50;
   Bool_t          HLT_AK8PFHT800_TrimMass50;
   Bool_t          HLT_AK8PFHT850_TrimMass50;
   Bool_t          HLT_AK8PFHT900_TrimMass50;
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
   Bool_t          L1_AlwaysTrue;
   Bool_t          L1_BPTX_AND_Ref1_VME;
   Bool_t          L1_BPTX_AND_Ref3_VME;
   Bool_t          L1_BPTX_AND_Ref4_VME;
   Bool_t          L1_BPTX_BeamGas_B1_VME;
   Bool_t          L1_BPTX_BeamGas_B2_VME;
   Bool_t          L1_BPTX_BeamGas_Ref1_VME;
   Bool_t          L1_BPTX_BeamGas_Ref2_VME;
   Bool_t          L1_BPTX_NotOR_VME;
   Bool_t          L1_BPTX_OR_Ref3_VME;
   Bool_t          L1_BPTX_OR_Ref4_VME;
   Bool_t          L1_BPTX_RefAND_VME;
   Bool_t          L1_BptxMinus;
   Bool_t          L1_BptxOR;
   Bool_t          L1_BptxPlus;
   Bool_t          L1_BptxXOR;
   Bool_t          L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;
   Bool_t          L1_DoubleEG6_HTT240er;
   Bool_t          L1_DoubleEG6_HTT250er;
   Bool_t          L1_DoubleEG6_HTT255er;
   Bool_t          L1_DoubleEG6_HTT270er;
   Bool_t          L1_DoubleEG6_HTT300er;
   Bool_t          L1_DoubleEG8er2p6_HTT255er;
   Bool_t          L1_DoubleEG8er2p6_HTT270er;
   Bool_t          L1_DoubleEG8er2p6_HTT300er;
   Bool_t          L1_DoubleEG_15_10;
   Bool_t          L1_DoubleEG_18_17;
   Bool_t          L1_DoubleEG_20_18;
   Bool_t          L1_DoubleEG_22_10;
   Bool_t          L1_DoubleEG_22_12;
   Bool_t          L1_DoubleEG_22_15;
   Bool_t          L1_DoubleEG_23_10;
   Bool_t          L1_DoubleEG_24_17;
   Bool_t          L1_DoubleEG_25_12;
   Bool_t          L1_DoubleEG_25_13;
   Bool_t          L1_DoubleEG_25_14;
   Bool_t          L1_DoubleEG_LooseIso23_10;
   Bool_t          L1_DoubleEG_LooseIso24_10;
   Bool_t          L1_DoubleIsoTau28er2p1;
   Bool_t          L1_DoubleIsoTau30er2p1;
   Bool_t          L1_DoubleIsoTau32er2p1;
   Bool_t          L1_DoubleIsoTau33er2p1;
   Bool_t          L1_DoubleIsoTau34er2p1;
   Bool_t          L1_DoubleIsoTau35er2p1;
   Bool_t          L1_DoubleIsoTau36er2p1;
   Bool_t          L1_DoubleIsoTau38er2p1;
   Bool_t          L1_DoubleJet100er2p3_dEta_Max1p6;
   Bool_t          L1_DoubleJet100er2p7;
   Bool_t          L1_DoubleJet112er2p3_dEta_Max1p6;
   Bool_t          L1_DoubleJet112er2p7;
   Bool_t          L1_DoubleJet120er2p7;
   Bool_t          L1_DoubleJet150er2p7;
   Bool_t          L1_DoubleJet30_Mass_Min300_dEta_Max1p5;
   Bool_t          L1_DoubleJet30_Mass_Min320_dEta_Max1p5;
   Bool_t          L1_DoubleJet30_Mass_Min340_dEta_Max1p5;
   Bool_t          L1_DoubleJet30_Mass_Min360_dEta_Max1p5;
   Bool_t          L1_DoubleJet30_Mass_Min380_dEta_Max1p5;
   Bool_t          L1_DoubleJet30_Mass_Min400_Mu10;
   Bool_t          L1_DoubleJet30_Mass_Min400_Mu6;
   Bool_t          L1_DoubleJet30_Mass_Min400_dEta_Max1p5;
   Bool_t          L1_DoubleJet35_rmovlp_IsoTau45_Mass_Min450;
   Bool_t          L1_DoubleJet40er2p7;
   Bool_t          L1_DoubleJet50er2p7;
   Bool_t          L1_DoubleJet60er2p7;
   Bool_t          L1_DoubleJet60er2p7_ETM100;
   Bool_t          L1_DoubleJet60er2p7_ETM60;
   Bool_t          L1_DoubleJet60er2p7_ETM70;
   Bool_t          L1_DoubleJet60er2p7_ETM80;
   Bool_t          L1_DoubleJet60er2p7_ETM90;
   Bool_t          L1_DoubleJet80er2p7;
   Bool_t          L1_DoubleJet_100_30_DoubleJet30_Mass_Min620;
   Bool_t          L1_DoubleJet_100_35_DoubleJet35_Mass_Min620;
   Bool_t          L1_DoubleJet_110_35_DoubleJet35_Mass_Min620;
   Bool_t          L1_DoubleJet_110_40_DoubleJet40_Mass_Min620;
   Bool_t          L1_DoubleJet_115_35_DoubleJet35_Mass_Min620;
   Bool_t          L1_DoubleJet_115_40_DoubleJet40_Mass_Min620;
   Bool_t          L1_DoubleJet_90_30_DoubleJet30_Mass_Min620;
   Bool_t          L1_DoubleLooseIsoEG22er2p1;
   Bool_t          L1_DoubleLooseIsoEG24er2p1;
   Bool_t          L1_DoubleMu0;
   Bool_t          L1_DoubleMu0_ETM40;
   Bool_t          L1_DoubleMu0_ETM55;
   Bool_t          L1_DoubleMu0_ETM60;
   Bool_t          L1_DoubleMu0_ETM65;
   Bool_t          L1_DoubleMu0_ETM70;
   Bool_t          L1_DoubleMu0_SQ;
   Bool_t          L1_DoubleMu0_SQ_OS;
   Bool_t          L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4;
   Bool_t          L1_DoubleMu0er1p4_dEta_Max1p8_OS;
   Bool_t          L1_DoubleMu0er1p5_SQ_OS;
   Bool_t          L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4;
   Bool_t          L1_DoubleMu0er1p5_SQ_dR_Max1p4;
   Bool_t          L1_DoubleMu0er2_SQ_dR_Max1p4;
   Bool_t          L1_DoubleMu18er2p1;
   Bool_t          L1_DoubleMu22er2p1;
   Bool_t          L1_DoubleMu3_OS_DoubleEG7p5Upsilon;
   Bool_t          L1_DoubleMu3_SQ_ETMHF40_Jet60_OR_DoubleJet30;
   Bool_t          L1_DoubleMu3_SQ_ETMHF50_Jet60_OR_DoubleJet30;
   Bool_t          L1_DoubleMu3_SQ_ETMHF60_Jet60_OR_DoubleJet30;
   Bool_t          L1_DoubleMu3_SQ_ETMHF70_Jet60_OR_DoubleJet30;
   Bool_t          L1_DoubleMu3_SQ_ETMHF80_Jet60_OR_DoubleJet30;
   Bool_t          L1_DoubleMu3_SQ_HTT100er;
   Bool_t          L1_DoubleMu3_SQ_HTT200er;
   Bool_t          L1_DoubleMu3_SQ_HTT220er;
   Bool_t          L1_DoubleMu3_SQ_HTT240er;
   Bool_t          L1_DoubleMu4_OS_EG12;
   Bool_t          L1_DoubleMu4_SQ_OS;
   Bool_t          L1_DoubleMu4_SQ_OS_dR_Max1p2;
   Bool_t          L1_DoubleMu4p5_SQ;
   Bool_t          L1_DoubleMu4p5_SQ_OS;
   Bool_t          L1_DoubleMu4p5_SQ_OS_dR_Max1p2;
   Bool_t          L1_DoubleMu4p5er2p0_SQ_OS;
   Bool_t          L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18;
   Bool_t          L1_DoubleMu5Upsilon_OS_DoubleEG3;
   Bool_t          L1_DoubleMu5_OS_EG12;
   Bool_t          L1_DoubleMu5_SQ_OS;
   Bool_t          L1_DoubleMu5_SQ_OS_Mass7to18;
   Bool_t          L1_DoubleMu6_SQ_OS;
   Bool_t          L1_DoubleMu7_EG7;
   Bool_t          L1_DoubleMu7_SQ_EG7;
   Bool_t          L1_DoubleMu8_SQ;
   Bool_t          L1_DoubleMu_10_0_dEta_Max1p8;
   Bool_t          L1_DoubleMu_11_4;
   Bool_t          L1_DoubleMu_12_5;
   Bool_t          L1_DoubleMu_12_8;
   Bool_t          L1_DoubleMu_13_6;
   Bool_t          L1_DoubleMu_15_5;
   Bool_t          L1_DoubleMu_15_5_SQ;
   Bool_t          L1_DoubleMu_15_7;
   Bool_t          L1_DoubleMu_15_7_SQ;
   Bool_t          L1_DoubleMu_15_7_SQ_Mass_Min4;
   Bool_t          L1_DoubleMu_20_2_SQ_Mass_Max20;
   Bool_t          L1_DoubleTau50er2p1;
   Bool_t          L1_DoubleTau70er2p1;
   Bool_t          L1_EG25er2p1_HTT125er;
   Bool_t          L1_EG27er2p1_HTT200er;
   Bool_t          L1_ETM100;
   Bool_t          L1_ETM100_Jet60_dPhi_Min0p4;
   Bool_t          L1_ETM105;
   Bool_t          L1_ETM110;
   Bool_t          L1_ETM110_Jet60_dPhi_Min0p4;
   Bool_t          L1_ETM115;
   Bool_t          L1_ETM120;
   Bool_t          L1_ETM150;
   Bool_t          L1_ETM30;
   Bool_t          L1_ETM40;
   Bool_t          L1_ETM50;
   Bool_t          L1_ETM60;
   Bool_t          L1_ETM70;
   Bool_t          L1_ETM75;
   Bool_t          L1_ETM75_Jet60_dPhi_Min0p4;
   Bool_t          L1_ETM80;
   Bool_t          L1_ETM80_Jet60_dPhi_Min0p4;
   Bool_t          L1_ETM85;
   Bool_t          L1_ETM90;
   Bool_t          L1_ETM90_Jet60_dPhi_Min0p4;
   Bool_t          L1_ETM95;
   Bool_t          L1_ETMHF100;
   Bool_t          L1_ETMHF100_HTT60er;
   Bool_t          L1_ETMHF100_Jet60_OR_DiJet30woTT28;
   Bool_t          L1_ETMHF100_Jet60_OR_DoubleJet30;
   Bool_t          L1_ETMHF100_Jet90_OR_DoubleJet45_OR_TripleJet30;
   Bool_t          L1_ETMHF110;
   Bool_t          L1_ETMHF110_HTT60er;
   Bool_t          L1_ETMHF110_Jet60_OR_DiJet30woTT28;
   Bool_t          L1_ETMHF110_Jet90_OR_DoubleJet45_OR_TripleJet30;
   Bool_t          L1_ETMHF120;
   Bool_t          L1_ETMHF120_HTT60er;
   Bool_t          L1_ETMHF120_Jet60_OR_DiJet30woTT28;
   Bool_t          L1_ETMHF150;
   Bool_t          L1_ETMHF70;
   Bool_t          L1_ETMHF70_Jet90_OR_DoubleJet45_OR_TripleJet30;
   Bool_t          L1_ETMHF80;
   Bool_t          L1_ETMHF80_HTT60er;
   Bool_t          L1_ETMHF80_Jet90_OR_DoubleJet45_OR_TripleJet30;
   Bool_t          L1_ETMHF90;
   Bool_t          L1_ETMHF90_HTT60er;
   Bool_t          L1_ETMHF90_Jet90_OR_DoubleJet45_OR_TripleJet30;
   Bool_t          L1_ETT100_BptxAND;
   Bool_t          L1_ETT110_BptxAND;
   Bool_t          L1_ETT40_BptxAND;
   Bool_t          L1_ETT50_BptxAND;
   Bool_t          L1_ETT60_BptxAND;
   Bool_t          L1_ETT70_BptxAND;
   Bool_t          L1_ETT75_BptxAND;
   Bool_t          L1_ETT80_BptxAND;
   Bool_t          L1_ETT85_BptxAND;
   Bool_t          L1_ETT90_BptxAND;
   Bool_t          L1_ETT95_BptxAND;
   Bool_t          L1_FirstBunchAfterTrain;
   Bool_t          L1_FirstBunchInTrain;
   Bool_t          L1_FirstCollisionInOrbit;
   Bool_t          L1_FirstCollisionInTrain;
   Bool_t          L1_HTT120er;
   Bool_t          L1_HTT160er;
   Bool_t          L1_HTT200er;
   Bool_t          L1_HTT220er;
   Bool_t          L1_HTT240er;
   Bool_t          L1_HTT250er_QuadJet_70_55_40_35_er2p5;
   Bool_t          L1_HTT255er;
   Bool_t          L1_HTT270er;
   Bool_t          L1_HTT280er;
   Bool_t          L1_HTT280er_QuadJet_70_55_40_35_er2p5;
   Bool_t          L1_HTT300er;
   Bool_t          L1_HTT300er_QuadJet_70_55_40_35_er2p5;
   Bool_t          L1_HTT320er;
   Bool_t          L1_HTT320er_QuadJet_70_55_40_40_er2p4;
   Bool_t          L1_HTT320er_QuadJet_70_55_40_40_er2p5;
   Bool_t          L1_HTT320er_QuadJet_70_55_45_45_er2p5;
   Bool_t          L1_HTT340er;
   Bool_t          L1_HTT340er_QuadJet_70_55_40_40_er2p5;
   Bool_t          L1_HTT340er_QuadJet_70_55_45_45_er2p5;
   Bool_t          L1_HTT380er;
   Bool_t          L1_HTT400er;
   Bool_t          L1_HTT450er;
   Bool_t          L1_HTT500er;
   Bool_t          L1_IsoEG33_Mt40;
   Bool_t          L1_IsoEG33_Mt44;
   Bool_t          L1_IsoEG33_Mt48;
   Bool_t          L1_IsoTau40er_ETM100;
   Bool_t          L1_IsoTau40er_ETM105;
   Bool_t          L1_IsoTau40er_ETM110;
   Bool_t          L1_IsoTau40er_ETM115;
   Bool_t          L1_IsoTau40er_ETM120;
   Bool_t          L1_IsoTau40er_ETM80;
   Bool_t          L1_IsoTau40er_ETM85;
   Bool_t          L1_IsoTau40er_ETM90;
   Bool_t          L1_IsoTau40er_ETM95;
   Bool_t          L1_IsoTau40er_ETMHF100;
   Bool_t          L1_IsoTau40er_ETMHF110;
   Bool_t          L1_IsoTau40er_ETMHF120;
   Bool_t          L1_IsoTau40er_ETMHF80;
   Bool_t          L1_IsoTau40er_ETMHF90;
   Bool_t          L1_IsolatedBunch;
   Bool_t          L1_LastCollisionInTrain;
   Bool_t          L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3;
   Bool_t          L1_LooseIsoEG24er2p1_HTT100er;
   Bool_t          L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3;
   Bool_t          L1_LooseIsoEG24er2p1_Jet26er2p7_dR_Min0p3;
   Bool_t          L1_LooseIsoEG24er2p1_TripleJet_26er2p7_26_26er2p7;
   Bool_t          L1_LooseIsoEG26er2p1_HTT100er;
   Bool_t          L1_LooseIsoEG26er2p1_Jet34er2p7_dR_Min0p3;
   Bool_t          L1_LooseIsoEG28er2p1_HTT100er;
   Bool_t          L1_LooseIsoEG28er2p1_Jet34er2p7_dR_Min0p3;
   Bool_t          L1_LooseIsoEG30er2p1_Jet34er2p7_dR_Min0p3;
   Bool_t          L1_MU20_EG15;
   Bool_t          L1_MinimumBiasHF0_AND_BptxAND;
   Bool_t          L1_MinimumBiasHF0_OR_BptxAND;
   Bool_t          L1_Mu10er2p1_ETM30;
   Bool_t          L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6;
   Bool_t          L1_Mu12_EG10;
   Bool_t          L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6;
   Bool_t          L1_Mu14er2p1_ETM30;
   Bool_t          L1_Mu15_HTT100er;
   Bool_t          L1_Mu18_HTT100er;
   Bool_t          L1_Mu18_Jet24er2p7;
   Bool_t          L1_Mu18er2p1_IsoTau26er2p1;
   Bool_t          L1_Mu18er2p1_Tau24er2p1;
   Bool_t          L1_Mu20_EG10;
   Bool_t          L1_Mu20_EG17;
   Bool_t          L1_Mu20_LooseIsoEG6;
   Bool_t          L1_Mu20er2p1_IsoTau26er2p1;
   Bool_t          L1_Mu20er2p1_IsoTau27er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau28er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau30er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau32er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau33er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau34er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau35er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau36er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau38er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau40er2p1;
   Bool_t          L1_Mu22er2p1_Tau50er2p1;
   Bool_t          L1_Mu22er2p1_Tau70er2p1;
   Bool_t          L1_Mu23_EG10;
   Bool_t          L1_Mu23_LooseIsoEG10;
   Bool_t          L1_Mu3_Jet120er2p7_dEta_Max0p4_dPhi_Max0p4;
   Bool_t          L1_Mu3_Jet16er2p7_dEta_Max0p4_dPhi_Max0p4;
   Bool_t          L1_Mu3_Jet30er2p5;
   Bool_t          L1_Mu3_Jet60er2p7_dEta_Max0p4_dPhi_Max0p4;
   Bool_t          L1_Mu5_EG15;
   Bool_t          L1_Mu5_EG20;
   Bool_t          L1_Mu5_EG23;
   Bool_t          L1_Mu5_LooseIsoEG18;
   Bool_t          L1_Mu5_LooseIsoEG20;
   Bool_t          L1_Mu6_DoubleEG10;
   Bool_t          L1_Mu6_DoubleEG17;
   Bool_t          L1_Mu6_HTT200er;
   Bool_t          L1_Mu6_HTT240er;
   Bool_t          L1_Mu6_HTT250er;
   Bool_t          L1_Mu7_EG23;
   Bool_t          L1_Mu7_LooseIsoEG20;
   Bool_t          L1_Mu7_LooseIsoEG23;
   Bool_t          L1_Mu8_HTT150er;
   Bool_t          L1_NotBptxOR;
   Bool_t          L1_QuadJet36er2p7_IsoTau52er2p1;
   Bool_t          L1_QuadJet36er2p7_Tau52;
   Bool_t          L1_QuadJet40er2p7;
   Bool_t          L1_QuadJet50er2p7;
   Bool_t          L1_QuadJet60er2p7;
   Bool_t          L1_QuadMu0;
   Bool_t          L1_SingleEG10;
   Bool_t          L1_SingleEG15;
   Bool_t          L1_SingleEG18;
   Bool_t          L1_SingleEG24;
   Bool_t          L1_SingleEG26;
   Bool_t          L1_SingleEG28;
   Bool_t          L1_SingleEG2_BptxAND;
   Bool_t          L1_SingleEG30;
   Bool_t          L1_SingleEG32;
   Bool_t          L1_SingleEG34;
   Bool_t          L1_SingleEG34er2p1;
   Bool_t          L1_SingleEG36;
   Bool_t          L1_SingleEG36er2p1;
   Bool_t          L1_SingleEG38;
   Bool_t          L1_SingleEG38er2p1;
   Bool_t          L1_SingleEG40;
   Bool_t          L1_SingleEG42;
   Bool_t          L1_SingleEG45;
   Bool_t          L1_SingleEG5;
   Bool_t          L1_SingleEG50;
   Bool_t          L1_SingleIsoEG18;
   Bool_t          L1_SingleIsoEG18er2p1;
   Bool_t          L1_SingleIsoEG20;
   Bool_t          L1_SingleIsoEG20er2p1;
   Bool_t          L1_SingleIsoEG22;
   Bool_t          L1_SingleIsoEG22er2p1;
   Bool_t          L1_SingleIsoEG24;
   Bool_t          L1_SingleIsoEG24er2p1;
   Bool_t          L1_SingleIsoEG26;
   Bool_t          L1_SingleIsoEG26er2p1;
   Bool_t          L1_SingleIsoEG28;
   Bool_t          L1_SingleIsoEG28er2p1;
   Bool_t          L1_SingleIsoEG30;
   Bool_t          L1_SingleIsoEG30er2p1;
   Bool_t          L1_SingleIsoEG32;
   Bool_t          L1_SingleIsoEG32er2p1;
   Bool_t          L1_SingleIsoEG33er2p1;
   Bool_t          L1_SingleIsoEG34;
   Bool_t          L1_SingleIsoEG34er2p1;
   Bool_t          L1_SingleIsoEG35;
   Bool_t          L1_SingleIsoEG35er2p1;
   Bool_t          L1_SingleIsoEG36;
   Bool_t          L1_SingleIsoEG36er2p1;
   Bool_t          L1_SingleIsoEG37;
   Bool_t          L1_SingleIsoEG38;
   Bool_t          L1_SingleIsoEG38er2p1;
   Bool_t          L1_SingleIsoEG40;
   Bool_t          L1_SingleIsoEG40er2p1;
   Bool_t          L1_SingleJet120;
   Bool_t          L1_SingleJet120_FWD;
   Bool_t          L1_SingleJet12_BptxAND;
   Bool_t          L1_SingleJet140;
   Bool_t          L1_SingleJet150;
   Bool_t          L1_SingleJet16;
   Bool_t          L1_SingleJet160;
   Bool_t          L1_SingleJet170;
   Bool_t          L1_SingleJet180;
   Bool_t          L1_SingleJet20;
   Bool_t          L1_SingleJet200;
   Bool_t          L1_SingleJet20er2p7_NotBptxOR;
   Bool_t          L1_SingleJet20er2p7_NotBptxOR_3BX;
   Bool_t          L1_SingleJet35;
   Bool_t          L1_SingleJet35_FWD;
   Bool_t          L1_SingleJet35_HFm;
   Bool_t          L1_SingleJet35_HFp;
   Bool_t          L1_SingleJet43er2p7_NotBptxOR_3BX;
   Bool_t          L1_SingleJet46er2p7_NotBptxOR_3BX;
   Bool_t          L1_SingleJet60;
   Bool_t          L1_SingleJet60_FWD;
   Bool_t          L1_SingleJet60_HFm;
   Bool_t          L1_SingleJet60_HFp;
   Bool_t          L1_SingleJet90;
   Bool_t          L1_SingleJet90_FWD;
   Bool_t          L1_SingleMu0_BMTF;
   Bool_t          L1_SingleMu0_EMTF;
   Bool_t          L1_SingleMu0_OMTF;
   Bool_t          L1_SingleMu10_LowQ;
   Bool_t          L1_SingleMu11_LowQ;
   Bool_t          L1_SingleMu12_LowQ_BMTF;
   Bool_t          L1_SingleMu12_LowQ_EMTF;
   Bool_t          L1_SingleMu12_LowQ_OMTF;
   Bool_t          L1_SingleMu14er2p1;
   Bool_t          L1_SingleMu16;
   Bool_t          L1_SingleMu16er2p1;
   Bool_t          L1_SingleMu18;
   Bool_t          L1_SingleMu18er2p1;
   Bool_t          L1_SingleMu20;
   Bool_t          L1_SingleMu20er2p1;
   Bool_t          L1_SingleMu22;
   Bool_t          L1_SingleMu22_BMTF;
   Bool_t          L1_SingleMu22_EMTF;
   Bool_t          L1_SingleMu22_OMTF;
   Bool_t          L1_SingleMu22er2p1;
   Bool_t          L1_SingleMu25;
   Bool_t          L1_SingleMu3;
   Bool_t          L1_SingleMu30;
   Bool_t          L1_SingleMu5;
   Bool_t          L1_SingleMu7;
   Bool_t          L1_SingleMuCosmics;
   Bool_t          L1_SingleMuCosmics_BMTF;
   Bool_t          L1_SingleMuCosmics_EMTF;
   Bool_t          L1_SingleMuCosmics_OMTF;
   Bool_t          L1_SingleMuOpen;
   Bool_t          L1_SingleMuOpen_NotBptxOR;
   Bool_t          L1_SingleMuOpen_NotBptxOR_3BX;
   Bool_t          L1_SingleTau100er2p1;
   Bool_t          L1_SingleTau120er2p1;
   Bool_t          L1_SingleTau130er2p1;
   Bool_t          L1_SingleTau140er2p1;
   Bool_t          L1_SingleTau20;
   Bool_t          L1_SingleTau80er2p1;
   Bool_t          L1_TripleEG_14_10_8;
   Bool_t          L1_TripleEG_18_17_8;
   Bool_t          L1_TripleEG_LooseIso20_10_5;
   Bool_t          L1_TripleJet_100_85_72_VBF;
   Bool_t          L1_TripleJet_105_85_76_VBF;
   Bool_t          L1_TripleJet_84_68_48_VBF;
   Bool_t          L1_TripleJet_88_72_56_VBF;
   Bool_t          L1_TripleJet_92_76_64_VBF;
   Bool_t          L1_TripleJet_98_83_71_VBF;
   Bool_t          L1_TripleMu0;
   Bool_t          L1_TripleMu0_OQ;
   Bool_t          L1_TripleMu3;
   Bool_t          L1_TripleMu3_SQ;
   Bool_t          L1_TripleMu_4_4_4;
   Bool_t          L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17;
   Bool_t          L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_8to14;
   Bool_t          L1_TripleMu_5SQ_3SQ_0OQ;
   Bool_t          L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9;
   Bool_t          L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9;
   Bool_t          L1_TripleMu_5_0_0;
   Bool_t          L1_TripleMu_5_3_3;
   Bool_t          L1_TripleMu_5_3p5_2p5;
   Bool_t          L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17;
   Bool_t          L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17;
   Bool_t          L1_TripleMu_5_5_3;
   Bool_t          L1_UnpairedBunchBptxMinus;
   Bool_t          L1_UnpairedBunchBptxPlus;
   Bool_t          L1_ZeroBias;
   Bool_t          L1_ZeroBias_copy;
   Float_t         Jet_pt_raw[njetmax];
   Float_t         Jet_pt_nom[njetmax];
   Float_t         Jet_mass_raw[njetmax];
   Float_t         Jet_mass_nom[njetmax];
   Float_t         Jet_corr_JEC[njetmax];
   Float_t         Jet_corr_JER[njetmax];
   Float_t         MET_pt_nom;
   Float_t         MET_phi_nom;
   Float_t         MET_pt_jer;
   Float_t         MET_phi_jer;
   Float_t         MET_pt_jerUp;
   Float_t         MET_phi_jerUp;
   Float_t         Jet_pt_jesAbsoluteStatUp[njetmax];
   Float_t         Jet_mass_jesAbsoluteStatUp[njetmax];
   Float_t         MET_pt_jesAbsoluteStatUp;
   Float_t         MET_phi_jesAbsoluteStatUp;
   Float_t         Jet_pt_jesAbsoluteScaleUp[njetmax];
   Float_t         Jet_mass_jesAbsoluteScaleUp[njetmax];
   Float_t         MET_pt_jesAbsoluteScaleUp;
   Float_t         MET_phi_jesAbsoluteScaleUp;
   Float_t         Jet_pt_jesAbsoluteFlavMapUp[njetmax];
   Float_t         Jet_mass_jesAbsoluteFlavMapUp[njetmax];
   Float_t         MET_pt_jesAbsoluteFlavMapUp;
   Float_t         MET_phi_jesAbsoluteFlavMapUp;
   Float_t         Jet_pt_jesAbsoluteMPFBiasUp[njetmax];
   Float_t         Jet_mass_jesAbsoluteMPFBiasUp[njetmax];
   Float_t         MET_pt_jesAbsoluteMPFBiasUp;
   Float_t         MET_phi_jesAbsoluteMPFBiasUp;
   Float_t         Jet_pt_jesFragmentationUp[njetmax];
   Float_t         Jet_mass_jesFragmentationUp[njetmax];
   Float_t         MET_pt_jesFragmentationUp;
   Float_t         MET_phi_jesFragmentationUp;
   Float_t         Jet_pt_jesSinglePionECALUp[njetmax];
   Float_t         Jet_mass_jesSinglePionECALUp[njetmax];
   Float_t         MET_pt_jesSinglePionECALUp;
   Float_t         MET_phi_jesSinglePionECALUp;
   Float_t         Jet_pt_jesSinglePionHCALUp[njetmax];
   Float_t         Jet_mass_jesSinglePionHCALUp[njetmax];
   Float_t         MET_pt_jesSinglePionHCALUp;
   Float_t         MET_phi_jesSinglePionHCALUp;
   Float_t         Jet_pt_jesFlavorQCDUp[njetmax];
   Float_t         Jet_mass_jesFlavorQCDUp[njetmax];
   Float_t         MET_pt_jesFlavorQCDUp;
   Float_t         MET_phi_jesFlavorQCDUp;
   Float_t         Jet_pt_jesTimePtEtaUp[njetmax];
   Float_t         Jet_mass_jesTimePtEtaUp[njetmax];
   Float_t         MET_pt_jesTimePtEtaUp;
   Float_t         MET_phi_jesTimePtEtaUp;
   Float_t         Jet_pt_jesRelativeJEREC1Up[njetmax];
   Float_t         Jet_mass_jesRelativeJEREC1Up[njetmax];
   Float_t         MET_pt_jesRelativeJEREC1Up;
   Float_t         MET_phi_jesRelativeJEREC1Up;
   Float_t         Jet_pt_jesRelativeJEREC2Up[njetmax];
   Float_t         Jet_mass_jesRelativeJEREC2Up[njetmax];
   Float_t         MET_pt_jesRelativeJEREC2Up;
   Float_t         MET_phi_jesRelativeJEREC2Up;
   Float_t         Jet_pt_jesRelativeJERHFUp[njetmax];
   Float_t         Jet_mass_jesRelativeJERHFUp[njetmax];
   Float_t         MET_pt_jesRelativeJERHFUp;
   Float_t         MET_phi_jesRelativeJERHFUp;
   Float_t         Jet_pt_jesRelativePtBBUp[njetmax];
   Float_t         Jet_mass_jesRelativePtBBUp[njetmax];
   Float_t         MET_pt_jesRelativePtBBUp;
   Float_t         MET_phi_jesRelativePtBBUp;
   Float_t         Jet_pt_jesRelativePtEC1Up[njetmax];
   Float_t         Jet_mass_jesRelativePtEC1Up[njetmax];
   Float_t         MET_pt_jesRelativePtEC1Up;
   Float_t         MET_phi_jesRelativePtEC1Up;
   Float_t         Jet_pt_jesRelativePtEC2Up[njetmax];
   Float_t         Jet_mass_jesRelativePtEC2Up[njetmax];
   Float_t         MET_pt_jesRelativePtEC2Up;
   Float_t         MET_phi_jesRelativePtEC2Up;
   Float_t         Jet_pt_jesRelativePtHFUp[njetmax];
   Float_t         Jet_mass_jesRelativePtHFUp[njetmax];
   Float_t         MET_pt_jesRelativePtHFUp;
   Float_t         MET_phi_jesRelativePtHFUp;
   Float_t         Jet_pt_jesRelativeBalUp[njetmax];
   Float_t         Jet_mass_jesRelativeBalUp[njetmax];
   Float_t         MET_pt_jesRelativeBalUp;
   Float_t         MET_phi_jesRelativeBalUp;
   Float_t         Jet_pt_jesRelativeSampleUp[njetmax];
   Float_t         Jet_mass_jesRelativeSampleUp[njetmax];
   Float_t         MET_pt_jesRelativeSampleUp;
   Float_t         MET_phi_jesRelativeSampleUp;
   Float_t         Jet_pt_jesRelativeFSRUp[njetmax];
   Float_t         Jet_mass_jesRelativeFSRUp[njetmax];
   Float_t         MET_pt_jesRelativeFSRUp;
   Float_t         MET_phi_jesRelativeFSRUp;
   Float_t         Jet_pt_jesRelativeStatFSRUp[njetmax];
   Float_t         Jet_mass_jesRelativeStatFSRUp[njetmax];
   Float_t         MET_pt_jesRelativeStatFSRUp;
   Float_t         MET_phi_jesRelativeStatFSRUp;
   Float_t         Jet_pt_jesRelativeStatECUp[njetmax];
   Float_t         Jet_mass_jesRelativeStatECUp[njetmax];
   Float_t         MET_pt_jesRelativeStatECUp;
   Float_t         MET_phi_jesRelativeStatECUp;
   Float_t         Jet_pt_jesRelativeStatHFUp[njetmax];
   Float_t         Jet_mass_jesRelativeStatHFUp[njetmax];
   Float_t         MET_pt_jesRelativeStatHFUp;
   Float_t         MET_phi_jesRelativeStatHFUp;
   Float_t         Jet_pt_jesPileUpDataMCUp[njetmax];
   Float_t         Jet_mass_jesPileUpDataMCUp[njetmax];
   Float_t         MET_pt_jesPileUpDataMCUp;
   Float_t         MET_phi_jesPileUpDataMCUp;
   Float_t         Jet_pt_jesPileUpPtRefUp[njetmax];
   Float_t         Jet_mass_jesPileUpPtRefUp[njetmax];
   Float_t         MET_pt_jesPileUpPtRefUp;
   Float_t         MET_phi_jesPileUpPtRefUp;
   Float_t         Jet_pt_jesPileUpPtBBUp[njetmax];
   Float_t         Jet_mass_jesPileUpPtBBUp[njetmax];
   Float_t         MET_pt_jesPileUpPtBBUp;
   Float_t         MET_phi_jesPileUpPtBBUp;
   Float_t         Jet_pt_jesPileUpPtEC1Up[njetmax];
   Float_t         Jet_mass_jesPileUpPtEC1Up[njetmax];
   Float_t         MET_pt_jesPileUpPtEC1Up;
   Float_t         MET_phi_jesPileUpPtEC1Up;
   Float_t         Jet_pt_jesPileUpPtEC2Up[njetmax];
   Float_t         Jet_mass_jesPileUpPtEC2Up[njetmax];
   Float_t         MET_pt_jesPileUpPtEC2Up;
   Float_t         MET_phi_jesPileUpPtEC2Up;
   Float_t         Jet_pt_jesPileUpPtHFUp[njetmax];
   Float_t         Jet_mass_jesPileUpPtHFUp[njetmax];
   Float_t         MET_pt_jesPileUpPtHFUp;
   Float_t         MET_phi_jesPileUpPtHFUp;
   Float_t         Jet_pt_jesPileUpMuZeroUp[njetmax];
   Float_t         Jet_mass_jesPileUpMuZeroUp[njetmax];
   Float_t         MET_pt_jesPileUpMuZeroUp;
   Float_t         MET_phi_jesPileUpMuZeroUp;
   Float_t         Jet_pt_jesPileUpEnvelopeUp[njetmax];
   Float_t         Jet_mass_jesPileUpEnvelopeUp[njetmax];
   Float_t         MET_pt_jesPileUpEnvelopeUp;
   Float_t         MET_phi_jesPileUpEnvelopeUp;
   Float_t         Jet_pt_jesSubTotalPileUpUp[njetmax];
   Float_t         Jet_mass_jesSubTotalPileUpUp[njetmax];
   Float_t         MET_pt_jesSubTotalPileUpUp;
   Float_t         MET_phi_jesSubTotalPileUpUp;
   Float_t         Jet_pt_jesSubTotalRelativeUp[njetmax];
   Float_t         Jet_mass_jesSubTotalRelativeUp[njetmax];
   Float_t         MET_pt_jesSubTotalRelativeUp;
   Float_t         MET_phi_jesSubTotalRelativeUp;
   Float_t         Jet_pt_jesSubTotalPtUp[njetmax];
   Float_t         Jet_mass_jesSubTotalPtUp[njetmax];
   Float_t         MET_pt_jesSubTotalPtUp;
   Float_t         MET_phi_jesSubTotalPtUp;
   Float_t         Jet_pt_jesSubTotalScaleUp[njetmax];
   Float_t         Jet_mass_jesSubTotalScaleUp[njetmax];
   Float_t         MET_pt_jesSubTotalScaleUp;
   Float_t         MET_phi_jesSubTotalScaleUp;
   Float_t         Jet_pt_jesSubTotalAbsoluteUp[njetmax];
   Float_t         Jet_mass_jesSubTotalAbsoluteUp[njetmax];
   Float_t         MET_pt_jesSubTotalAbsoluteUp;
   Float_t         MET_phi_jesSubTotalAbsoluteUp;
   Float_t         Jet_pt_jesSubTotalMCUp[njetmax];
   Float_t         Jet_mass_jesSubTotalMCUp[njetmax];
   Float_t         MET_pt_jesSubTotalMCUp;
   Float_t         MET_phi_jesSubTotalMCUp;
   Float_t         Jet_pt_jesTotalUp[njetmax];
   Float_t         Jet_mass_jesTotalUp[njetmax];
   Float_t         MET_pt_jesTotalUp;
   Float_t         MET_phi_jesTotalUp;
   Float_t         Jet_pt_jesTotalNoFlavorUp[njetmax];
   Float_t         Jet_mass_jesTotalNoFlavorUp[njetmax];
   Float_t         MET_pt_jesTotalNoFlavorUp;
   Float_t         MET_phi_jesTotalNoFlavorUp;
   Float_t         Jet_pt_jesTotalNoTimeUp[njetmax];
   Float_t         Jet_mass_jesTotalNoTimeUp[njetmax];
   Float_t         MET_pt_jesTotalNoTimeUp;
   Float_t         MET_phi_jesTotalNoTimeUp;
   Float_t         Jet_pt_jesTotalNoFlavorNoTimeUp[njetmax];
   Float_t         Jet_mass_jesTotalNoFlavorNoTimeUp[njetmax];
   Float_t         MET_pt_jesTotalNoFlavorNoTimeUp;
   Float_t         MET_phi_jesTotalNoFlavorNoTimeUp;
   Float_t         Jet_pt_jesFlavorZJetUp[njetmax];
   Float_t         Jet_mass_jesFlavorZJetUp[njetmax];
   Float_t         MET_pt_jesFlavorZJetUp;
   Float_t         MET_phi_jesFlavorZJetUp;
   Float_t         Jet_pt_jesFlavorPhotonJetUp[njetmax];
   Float_t         Jet_mass_jesFlavorPhotonJetUp[njetmax];
   Float_t         MET_pt_jesFlavorPhotonJetUp;
   Float_t         MET_phi_jesFlavorPhotonJetUp;
   Float_t         Jet_pt_jesFlavorPureGluonUp[njetmax];
   Float_t         Jet_mass_jesFlavorPureGluonUp[njetmax];
   Float_t         MET_pt_jesFlavorPureGluonUp;
   Float_t         MET_phi_jesFlavorPureGluonUp;
   Float_t         Jet_pt_jesFlavorPureQuarkUp[njetmax];
   Float_t         Jet_mass_jesFlavorPureQuarkUp[njetmax];
   Float_t         MET_pt_jesFlavorPureQuarkUp;
   Float_t         MET_phi_jesFlavorPureQuarkUp;
   Float_t         Jet_pt_jesFlavorPureCharmUp[njetmax];
   Float_t         Jet_mass_jesFlavorPureCharmUp[njetmax];
   Float_t         MET_pt_jesFlavorPureCharmUp;
   Float_t         MET_phi_jesFlavorPureCharmUp;
   Float_t         Jet_pt_jesFlavorPureBottomUp[njetmax];
   Float_t         Jet_mass_jesFlavorPureBottomUp[njetmax];
   Float_t         MET_pt_jesFlavorPureBottomUp;
   Float_t         MET_phi_jesFlavorPureBottomUp;
   Float_t         Jet_pt_jesTimeRunBUp[njetmax];
   Float_t         Jet_mass_jesTimeRunBUp[njetmax];
   Float_t         MET_pt_jesTimeRunBUp;
   Float_t         MET_phi_jesTimeRunBUp;
   Float_t         Jet_pt_jesTimeRunCUp[njetmax];
   Float_t         Jet_mass_jesTimeRunCUp[njetmax];
   Float_t         MET_pt_jesTimeRunCUp;
   Float_t         MET_phi_jesTimeRunCUp;
   Float_t         Jet_pt_jesTimeRunDEUp[njetmax];
   Float_t         Jet_mass_jesTimeRunDEUp[njetmax];
   Float_t         MET_pt_jesTimeRunDEUp;
   Float_t         MET_phi_jesTimeRunDEUp;
   Float_t         Jet_pt_jesTimeRunFUp[njetmax];
   Float_t         Jet_mass_jesTimeRunFUp[njetmax];
   Float_t         MET_pt_jesTimeRunFUp;
   Float_t         MET_phi_jesTimeRunFUp;
   Float_t         Jet_pt_jesCorrelationGroupMPFInSituUp[njetmax];
   Float_t         Jet_mass_jesCorrelationGroupMPFInSituUp[njetmax];
   Float_t         MET_pt_jesCorrelationGroupMPFInSituUp;
   Float_t         MET_phi_jesCorrelationGroupMPFInSituUp;
   Float_t         Jet_pt_jesCorrelationGroupIntercalibrationUp[njetmax];
   Float_t         Jet_mass_jesCorrelationGroupIntercalibrationUp[njetmax];
   Float_t         MET_pt_jesCorrelationGroupIntercalibrationUp;
   Float_t         MET_phi_jesCorrelationGroupIntercalibrationUp;
   Float_t         Jet_pt_jesCorrelationGroupbJESUp[njetmax];
   Float_t         Jet_mass_jesCorrelationGroupbJESUp[njetmax];
   Float_t         MET_pt_jesCorrelationGroupbJESUp;
   Float_t         MET_phi_jesCorrelationGroupbJESUp;
   Float_t         Jet_pt_jesCorrelationGroupFlavorUp[njetmax];
   Float_t         Jet_mass_jesCorrelationGroupFlavorUp[njetmax];
   Float_t         MET_pt_jesCorrelationGroupFlavorUp;
   Float_t         MET_phi_jesCorrelationGroupFlavorUp;
   Float_t         Jet_pt_jesCorrelationGroupUncorrelatedUp[njetmax];
   Float_t         Jet_mass_jesCorrelationGroupUncorrelatedUp[njetmax];
   Float_t         MET_pt_jesCorrelationGroupUncorrelatedUp;
   Float_t         MET_phi_jesCorrelationGroupUncorrelatedUp;
   Float_t         MET_pt_unclustEnUp;
   Float_t         MET_phi_unclustEnUp;
   Float_t         Jet_pt_jerDown[njetmax];
   Float_t         Jet_mass_jerDown[njetmax];
   Float_t         Jet_mass_jmsDown[njetmax];   
   Float_t         Jet_mass_jmrDown[njetmax];   
   Float_t         MET_pt_jerDown;
   Float_t         MET_phi_jerDown;
   Float_t         Jet_pt_jesAbsoluteStatDown[njetmax];
   Float_t         Jet_mass_jesAbsoluteStatDown[njetmax];
   Float_t         MET_pt_jesAbsoluteStatDown;
   Float_t         MET_phi_jesAbsoluteStatDown;
   Float_t         Jet_pt_jesAbsoluteScaleDown[njetmax];
   Float_t         Jet_mass_jesAbsoluteScaleDown[njetmax];
   Float_t         MET_pt_jesAbsoluteScaleDown;
   Float_t         MET_phi_jesAbsoluteScaleDown;
   Float_t         Jet_pt_jesAbsoluteFlavMapDown[njetmax];
   Float_t         Jet_mass_jesAbsoluteFlavMapDown[njetmax];
   Float_t         MET_pt_jesAbsoluteFlavMapDown;
   Float_t         MET_phi_jesAbsoluteFlavMapDown;
   Float_t         Jet_pt_jesAbsoluteMPFBiasDown[njetmax];
   Float_t         Jet_mass_jesAbsoluteMPFBiasDown[njetmax];
   Float_t         MET_pt_jesAbsoluteMPFBiasDown;
   Float_t         MET_phi_jesAbsoluteMPFBiasDown;
   Float_t         Jet_pt_jesFragmentationDown[njetmax];
   Float_t         Jet_mass_jesFragmentationDown[njetmax];
   Float_t         MET_pt_jesFragmentationDown;
   Float_t         MET_phi_jesFragmentationDown;
   Float_t         Jet_pt_jesSinglePionECALDown[njetmax];
   Float_t         Jet_mass_jesSinglePionECALDown[njetmax];
   Float_t         MET_pt_jesSinglePionECALDown;
   Float_t         MET_phi_jesSinglePionECALDown;
   Float_t         Jet_pt_jesSinglePionHCALDown[njetmax];
   Float_t         Jet_mass_jesSinglePionHCALDown[njetmax];
   Float_t         MET_pt_jesSinglePionHCALDown;
   Float_t         MET_phi_jesSinglePionHCALDown;
   Float_t         Jet_pt_jesFlavorQCDDown[njetmax];
   Float_t         Jet_mass_jesFlavorQCDDown[njetmax];
   Float_t         MET_pt_jesFlavorQCDDown;
   Float_t         MET_phi_jesFlavorQCDDown;
   Float_t         Jet_pt_jesTimePtEtaDown[njetmax];
   Float_t         Jet_mass_jesTimePtEtaDown[njetmax];
   Float_t         MET_pt_jesTimePtEtaDown;
   Float_t         MET_phi_jesTimePtEtaDown;
   Float_t         Jet_pt_jesRelativeJEREC1Down[njetmax];
   Float_t         Jet_mass_jesRelativeJEREC1Down[njetmax];
   Float_t         MET_pt_jesRelativeJEREC1Down;
   Float_t         MET_phi_jesRelativeJEREC1Down;
   Float_t         Jet_pt_jesRelativeJEREC2Down[njetmax];
   Float_t         Jet_mass_jesRelativeJEREC2Down[njetmax];
   Float_t         MET_pt_jesRelativeJEREC2Down;
   Float_t         MET_phi_jesRelativeJEREC2Down;
   Float_t         Jet_pt_jesRelativeJERHFDown[njetmax];
   Float_t         Jet_mass_jesRelativeJERHFDown[njetmax];
   Float_t         MET_pt_jesRelativeJERHFDown;
   Float_t         MET_phi_jesRelativeJERHFDown;
   Float_t         Jet_pt_jesRelativePtBBDown[njetmax];
   Float_t         Jet_mass_jesRelativePtBBDown[njetmax];
   Float_t         MET_pt_jesRelativePtBBDown;
   Float_t         MET_phi_jesRelativePtBBDown;
   Float_t         Jet_pt_jesRelativePtEC1Down[njetmax];
   Float_t         Jet_mass_jesRelativePtEC1Down[njetmax];
   Float_t         MET_pt_jesRelativePtEC1Down;
   Float_t         MET_phi_jesRelativePtEC1Down;
   Float_t         Jet_pt_jesRelativePtEC2Down[njetmax];
   Float_t         Jet_mass_jesRelativePtEC2Down[njetmax];
   Float_t         MET_pt_jesRelativePtEC2Down;
   Float_t         MET_phi_jesRelativePtEC2Down;
   Float_t         Jet_pt_jesRelativePtHFDown[njetmax];
   Float_t         Jet_mass_jesRelativePtHFDown[njetmax];
   Float_t         MET_pt_jesRelativePtHFDown;
   Float_t         MET_phi_jesRelativePtHFDown;
   Float_t         Jet_pt_jesRelativeBalDown[njetmax];
   Float_t         Jet_mass_jesRelativeBalDown[njetmax];
   Float_t         MET_pt_jesRelativeBalDown;
   Float_t         MET_phi_jesRelativeBalDown;
   Float_t         Jet_pt_jesRelativeSampleDown[njetmax];
   Float_t         Jet_mass_jesRelativeSampleDown[njetmax];
   Float_t         MET_pt_jesRelativeSampleDown;
   Float_t         MET_phi_jesRelativeSampleDown;
   Float_t         Jet_pt_jesRelativeFSRDown[njetmax];
   Float_t         Jet_mass_jesRelativeFSRDown[njetmax];
   Float_t         MET_pt_jesRelativeFSRDown;
   Float_t         MET_phi_jesRelativeFSRDown;
   Float_t         Jet_pt_jesRelativeStatFSRDown[njetmax];
   Float_t         Jet_mass_jesRelativeStatFSRDown[njetmax];
   Float_t         MET_pt_jesRelativeStatFSRDown;
   Float_t         MET_phi_jesRelativeStatFSRDown;
   Float_t         Jet_pt_jesRelativeStatECDown[njetmax];
   Float_t         Jet_mass_jesRelativeStatECDown[njetmax];
   Float_t         MET_pt_jesRelativeStatECDown;
   Float_t         MET_phi_jesRelativeStatECDown;
   Float_t         Jet_pt_jesRelativeStatHFDown[njetmax];
   Float_t         Jet_mass_jesRelativeStatHFDown[njetmax];
   Float_t         MET_pt_jesRelativeStatHFDown;
   Float_t         MET_phi_jesRelativeStatHFDown;
   Float_t         Jet_pt_jesPileUpDataMCDown[njetmax];
   Float_t         Jet_mass_jesPileUpDataMCDown[njetmax];
   Float_t         MET_pt_jesPileUpDataMCDown;
   Float_t         MET_phi_jesPileUpDataMCDown;
   Float_t         Jet_pt_jesPileUpPtRefDown[njetmax];
   Float_t         Jet_mass_jesPileUpPtRefDown[njetmax];
   Float_t         MET_pt_jesPileUpPtRefDown;
   Float_t         MET_phi_jesPileUpPtRefDown;
   Float_t         Jet_pt_jesPileUpPtBBDown[njetmax];
   Float_t         Jet_mass_jesPileUpPtBBDown[njetmax];
   Float_t         MET_pt_jesPileUpPtBBDown;
   Float_t         MET_phi_jesPileUpPtBBDown;
   Float_t         Jet_pt_jesPileUpPtEC1Down[njetmax];
   Float_t         Jet_mass_jesPileUpPtEC1Down[njetmax];
   Float_t         MET_pt_jesPileUpPtEC1Down;
   Float_t         MET_phi_jesPileUpPtEC1Down;
   Float_t         Jet_pt_jesPileUpPtEC2Down[njetmax];
   Float_t         Jet_mass_jesPileUpPtEC2Down[njetmax];
   Float_t         MET_pt_jesPileUpPtEC2Down;
   Float_t         MET_phi_jesPileUpPtEC2Down;
   Float_t         Jet_pt_jesPileUpPtHFDown[njetmax];
   Float_t         Jet_mass_jesPileUpPtHFDown[njetmax];
   Float_t         MET_pt_jesPileUpPtHFDown;
   Float_t         MET_phi_jesPileUpPtHFDown;
   Float_t         Jet_pt_jesPileUpMuZeroDown[njetmax];
   Float_t         Jet_mass_jesPileUpMuZeroDown[njetmax];
   Float_t         MET_pt_jesPileUpMuZeroDown;
   Float_t         MET_phi_jesPileUpMuZeroDown;
   Float_t         Jet_pt_jesPileUpEnvelopeDown[njetmax];
   Float_t         Jet_mass_jesPileUpEnvelopeDown[njetmax];
   Float_t         MET_pt_jesPileUpEnvelopeDown;
   Float_t         MET_phi_jesPileUpEnvelopeDown;
   Float_t         Jet_pt_jesSubTotalPileUpDown[njetmax];
   Float_t         Jet_mass_jesSubTotalPileUpDown[njetmax];
   Float_t         MET_pt_jesSubTotalPileUpDown;
   Float_t         MET_phi_jesSubTotalPileUpDown;
   Float_t         Jet_pt_jesSubTotalRelativeDown[njetmax];
   Float_t         Jet_mass_jesSubTotalRelativeDown[njetmax];
   Float_t         MET_pt_jesSubTotalRelativeDown;
   Float_t         MET_phi_jesSubTotalRelativeDown;
   Float_t         Jet_pt_jesSubTotalPtDown[njetmax];
   Float_t         Jet_mass_jesSubTotalPtDown[njetmax];
   Float_t         MET_pt_jesSubTotalPtDown;
   Float_t         MET_phi_jesSubTotalPtDown;
   Float_t         Jet_pt_jesSubTotalScaleDown[njetmax];
   Float_t         Jet_mass_jesSubTotalScaleDown[njetmax];
   Float_t         MET_pt_jesSubTotalScaleDown;
   Float_t         MET_phi_jesSubTotalScaleDown;
   Float_t         Jet_pt_jesSubTotalAbsoluteDown[njetmax];
   Float_t         Jet_mass_jesSubTotalAbsoluteDown[njetmax];
   Float_t         MET_pt_jesSubTotalAbsoluteDown;
   Float_t         MET_phi_jesSubTotalAbsoluteDown;
   Float_t         Jet_pt_jesSubTotalMCDown[njetmax];
   Float_t         Jet_mass_jesSubTotalMCDown[njetmax];
   Float_t         MET_pt_jesSubTotalMCDown;
   Float_t         MET_phi_jesSubTotalMCDown;
   Float_t         Jet_pt_jesTotalDown[njetmax];
   Float_t         Jet_mass_jesTotalDown[njetmax];
   Float_t         MET_pt_jesTotalDown;
   Float_t         MET_phi_jesTotalDown;
   Float_t         Jet_pt_jesTotalNoFlavorDown[njetmax];
   Float_t         Jet_mass_jesTotalNoFlavorDown[njetmax];
   Float_t         MET_pt_jesTotalNoFlavorDown;
   Float_t         MET_phi_jesTotalNoFlavorDown;
   Float_t         Jet_pt_jesTotalNoTimeDown[njetmax];
   Float_t         Jet_mass_jesTotalNoTimeDown[njetmax];
   Float_t         MET_pt_jesTotalNoTimeDown;
   Float_t         MET_phi_jesTotalNoTimeDown;
   Float_t         Jet_pt_jesTotalNoFlavorNoTimeDown[njetmax];
   Float_t         Jet_mass_jesTotalNoFlavorNoTimeDown[njetmax];
   Float_t         MET_pt_jesTotalNoFlavorNoTimeDown;
   Float_t         MET_phi_jesTotalNoFlavorNoTimeDown;
   Float_t         Jet_pt_jesFlavorZJetDown[njetmax];
   Float_t         Jet_mass_jesFlavorZJetDown[njetmax];
   Float_t         MET_pt_jesFlavorZJetDown;
   Float_t         MET_phi_jesFlavorZJetDown;
   Float_t         Jet_pt_jesFlavorPhotonJetDown[njetmax];
   Float_t         Jet_mass_jesFlavorPhotonJetDown[njetmax];
   Float_t         MET_pt_jesFlavorPhotonJetDown;
   Float_t         MET_phi_jesFlavorPhotonJetDown;
   Float_t         Jet_pt_jesFlavorPureGluonDown[njetmax];
   Float_t         Jet_mass_jesFlavorPureGluonDown[njetmax];
   Float_t         MET_pt_jesFlavorPureGluonDown;
   Float_t         MET_phi_jesFlavorPureGluonDown;
   Float_t         Jet_pt_jesFlavorPureQuarkDown[njetmax];
   Float_t         Jet_mass_jesFlavorPureQuarkDown[njetmax];
   Float_t         MET_pt_jesFlavorPureQuarkDown;
   Float_t         MET_phi_jesFlavorPureQuarkDown;
   Float_t         Jet_pt_jesFlavorPureCharmDown[njetmax];
   Float_t         Jet_mass_jesFlavorPureCharmDown[njetmax];
   Float_t         MET_pt_jesFlavorPureCharmDown;
   Float_t         MET_phi_jesFlavorPureCharmDown;
   Float_t         Jet_pt_jesFlavorPureBottomDown[njetmax];
   Float_t         Jet_mass_jesFlavorPureBottomDown[njetmax];
   Float_t         MET_pt_jesFlavorPureBottomDown;
   Float_t         MET_phi_jesFlavorPureBottomDown;
   Float_t         Jet_pt_jesTimeRunBDown[njetmax];
   Float_t         Jet_mass_jesTimeRunBDown[njetmax];
   Float_t         MET_pt_jesTimeRunBDown;
   Float_t         MET_phi_jesTimeRunBDown;
   Float_t         Jet_pt_jesTimeRunCDown[njetmax];
   Float_t         Jet_mass_jesTimeRunCDown[njetmax];
   Float_t         MET_pt_jesTimeRunCDown;
   Float_t         MET_phi_jesTimeRunCDown;
   Float_t         Jet_pt_jesTimeRunDEDown[njetmax];
   Float_t         Jet_mass_jesTimeRunDEDown[njetmax];
   Float_t         MET_pt_jesTimeRunDEDown;
   Float_t         MET_phi_jesTimeRunDEDown;
   Float_t         Jet_pt_jesTimeRunFDown[njetmax];
   Float_t         Jet_mass_jesTimeRunFDown[njetmax];
   Float_t         MET_pt_jesTimeRunFDown;
   Float_t         MET_phi_jesTimeRunFDown;
   Float_t         Jet_pt_jesCorrelationGroupMPFInSituDown[njetmax];
   Float_t         Jet_mass_jesCorrelationGroupMPFInSituDown[njetmax];
   Float_t         MET_pt_jesCorrelationGroupMPFInSituDown;
   Float_t         MET_phi_jesCorrelationGroupMPFInSituDown;
   Float_t         Jet_pt_jesCorrelationGroupIntercalibrationDown[njetmax];
   Float_t         Jet_mass_jesCorrelationGroupIntercalibrationDown[njetmax];
   Float_t         MET_pt_jesCorrelationGroupIntercalibrationDown;
   Float_t         MET_phi_jesCorrelationGroupIntercalibrationDown;
   Float_t         Jet_pt_jesCorrelationGroupbJESDown[njetmax];
   Float_t         Jet_mass_jesCorrelationGroupbJESDown[njetmax];
   Float_t         MET_pt_jesCorrelationGroupbJESDown;
   Float_t         MET_phi_jesCorrelationGroupbJESDown;
   Float_t         Jet_pt_jesCorrelationGroupFlavorDown[njetmax];
   Float_t         Jet_mass_jesCorrelationGroupFlavorDown[njetmax];
   Float_t         MET_pt_jesCorrelationGroupFlavorDown;
   Float_t         MET_phi_jesCorrelationGroupFlavorDown;
   Float_t         Jet_pt_jesCorrelationGroupUncorrelatedDown[njetmax];
   Float_t         Jet_mass_jesCorrelationGroupUncorrelatedDown[njetmax];
   Float_t         MET_pt_jesCorrelationGroupUncorrelatedDown;
   Float_t         MET_phi_jesCorrelationGroupUncorrelatedDown;
   Float_t         MET_pt_unclustEnDown;
   Float_t         MET_phi_unclustEnDown;
   Float_t         Jet_pt_jerUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jerUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jmrUp[njetmax];   //[nJet]
   Float_t         Jet_mass_jmsUp[njetmax];   //[nJet]
   Float_t         FatJet_pt_raw[njetmax];
   Float_t         FatJet_pt_nom[njetmax];
   Float_t         FatJet_mass_raw[njetmax];
   Float_t         FatJet_mass_nom[njetmax];
   Float_t         FatJet_corr_JEC[njetmax];
   Float_t         FatJet_corr_JER[njetmax];
   Float_t         FatJet_corr_JMS[njetmax];
   Float_t         FatJet_corr_JMR[njetmax];
   Float_t         FatJet_msoftdrop_raw[njetmax];
   Float_t         FatJet_msoftdrop_nom[njetmax];
   Float_t         FatJet_msoftdrop_corr_JMR[njetmax];
   Float_t         FatJet_msoftdrop_corr_JMS[njetmax];
   Float_t         FatJet_msoftdrop_corr_PUPPI[njetmax];
   Float_t         FatJet_msoftdrop_tau21DDT_nom[njetmax];
   Float_t         FatJet_pt_jerUp[njetmax];
   Float_t         FatJet_mass_jerUp[njetmax];
   Float_t         FatJet_mass_jmrUp[njetmax];
   Float_t         FatJet_mass_jmsUp[njetmax];
   Float_t         FatJet_msoftdrop_jerUp[njetmax];
   Float_t         FatJet_msoftdrop_jmrUp[njetmax];
   Float_t         FatJet_msoftdrop_jmsUp[njetmax];
   Float_t         FatJet_msoftdrop_tau21DDT_jerUp[njetmax];
   Float_t         FatJet_msoftdrop_tau21DDT_jmrUp[njetmax];
   Float_t         FatJet_msoftdrop_tau21DDT_jmsUp[njetmax];
   Float_t         FatJet_pt_jesAbsoluteStatUp[njetmax];
   Float_t         FatJet_mass_jesAbsoluteStatUp[njetmax];
   Float_t         FatJet_msoftdrop_jesAbsoluteStatUp[njetmax];
   Float_t         FatJet_pt_jesAbsoluteScaleUp[njetmax];
   Float_t         FatJet_mass_jesAbsoluteScaleUp[njetmax];
   Float_t         FatJet_msoftdrop_jesAbsoluteScaleUp[njetmax];
   Float_t         FatJet_pt_jesAbsoluteFlavMapUp[njetmax];
   Float_t         FatJet_mass_jesAbsoluteFlavMapUp[njetmax];
   Float_t         FatJet_msoftdrop_jesAbsoluteFlavMapUp[njetmax];
   Float_t         FatJet_pt_jesAbsoluteMPFBiasUp[njetmax];
   Float_t         FatJet_mass_jesAbsoluteMPFBiasUp[njetmax];
   Float_t         FatJet_msoftdrop_jesAbsoluteMPFBiasUp[njetmax];
   Float_t         FatJet_pt_jesFragmentationUp[njetmax];
   Float_t         FatJet_mass_jesFragmentationUp[njetmax];
   Float_t         FatJet_msoftdrop_jesFragmentationUp[njetmax];
   Float_t         FatJet_pt_jesSinglePionECALUp[njetmax];
   Float_t         FatJet_mass_jesSinglePionECALUp[njetmax];
   Float_t         FatJet_msoftdrop_jesSinglePionECALUp[njetmax];
   Float_t         FatJet_pt_jesSinglePionHCALUp[njetmax];
   Float_t         FatJet_mass_jesSinglePionHCALUp[njetmax];
   Float_t         FatJet_msoftdrop_jesSinglePionHCALUp[njetmax];
   Float_t         FatJet_pt_jesFlavorQCDUp[njetmax];
   Float_t         FatJet_mass_jesFlavorQCDUp[njetmax];
   Float_t         FatJet_msoftdrop_jesFlavorQCDUp[njetmax];
   Float_t         FatJet_pt_jesTimePtEtaUp[njetmax];
   Float_t         FatJet_mass_jesTimePtEtaUp[njetmax];
   Float_t         FatJet_msoftdrop_jesTimePtEtaUp[njetmax];
   Float_t         FatJet_pt_jesRelativeJEREC1Up[njetmax];
   Float_t         FatJet_mass_jesRelativeJEREC1Up[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativeJEREC1Up[njetmax];
   Float_t         FatJet_pt_jesRelativeJEREC2Up[njetmax];
   Float_t         FatJet_mass_jesRelativeJEREC2Up[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativeJEREC2Up[njetmax];
   Float_t         FatJet_pt_jesRelativeJERHFUp[njetmax];
   Float_t         FatJet_mass_jesRelativeJERHFUp[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativeJERHFUp[njetmax];
   Float_t         FatJet_pt_jesRelativePtBBUp[njetmax];
   Float_t         FatJet_mass_jesRelativePtBBUp[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativePtBBUp[njetmax];
   Float_t         FatJet_pt_jesRelativePtEC1Up[njetmax];
   Float_t         FatJet_mass_jesRelativePtEC1Up[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativePtEC1Up[njetmax];
   Float_t         FatJet_pt_jesRelativePtEC2Up[njetmax];
   Float_t         FatJet_mass_jesRelativePtEC2Up[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativePtEC2Up[njetmax];
   Float_t         FatJet_pt_jesRelativePtHFUp[njetmax];
   Float_t         FatJet_mass_jesRelativePtHFUp[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativePtHFUp[njetmax];
   Float_t         FatJet_pt_jesRelativeBalUp[njetmax];
   Float_t         FatJet_mass_jesRelativeBalUp[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativeBalUp[njetmax];
   Float_t         FatJet_pt_jesRelativeSampleUp[njetmax];
   Float_t         FatJet_mass_jesRelativeSampleUp[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativeSampleUp[njetmax];
   Float_t         FatJet_pt_jesRelativeFSRUp[njetmax];
   Float_t         FatJet_mass_jesRelativeFSRUp[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativeFSRUp[njetmax];
   Float_t         FatJet_pt_jesRelativeStatFSRUp[njetmax];
   Float_t         FatJet_mass_jesRelativeStatFSRUp[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativeStatFSRUp[njetmax];
   Float_t         FatJet_pt_jesRelativeStatECUp[njetmax];
   Float_t         FatJet_mass_jesRelativeStatECUp[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativeStatECUp[njetmax];
   Float_t         FatJet_pt_jesRelativeStatHFUp[njetmax];
   Float_t         FatJet_mass_jesRelativeStatHFUp[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativeStatHFUp[njetmax];
   Float_t         FatJet_pt_jesPileUpDataMCUp[njetmax];
   Float_t         FatJet_mass_jesPileUpDataMCUp[njetmax];
   Float_t         FatJet_msoftdrop_jesPileUpDataMCUp[njetmax];
   Float_t         FatJet_pt_jesPileUpPtRefUp[njetmax];
   Float_t         FatJet_mass_jesPileUpPtRefUp[njetmax];
   Float_t         FatJet_msoftdrop_jesPileUpPtRefUp[njetmax];
   Float_t         FatJet_pt_jesPileUpPtBBUp[njetmax];
   Float_t         FatJet_mass_jesPileUpPtBBUp[njetmax];
   Float_t         FatJet_msoftdrop_jesPileUpPtBBUp[njetmax];
   Float_t         FatJet_pt_jesPileUpPtEC1Up[njetmax];
   Float_t         FatJet_mass_jesPileUpPtEC1Up[njetmax];
   Float_t         FatJet_msoftdrop_jesPileUpPtEC1Up[njetmax];
   Float_t         FatJet_pt_jesPileUpPtEC2Up[njetmax];
   Float_t         FatJet_mass_jesPileUpPtEC2Up[njetmax];
   Float_t         FatJet_msoftdrop_jesPileUpPtEC2Up[njetmax];
   Float_t         FatJet_pt_jesPileUpPtHFUp[njetmax];
   Float_t         FatJet_mass_jesPileUpPtHFUp[njetmax];
   Float_t         FatJet_msoftdrop_jesPileUpPtHFUp[njetmax];
   Float_t         FatJet_pt_jesPileUpMuZeroUp[njetmax];
   Float_t         FatJet_mass_jesPileUpMuZeroUp[njetmax];
   Float_t         FatJet_msoftdrop_jesPileUpMuZeroUp[njetmax];
   Float_t         FatJet_pt_jesPileUpEnvelopeUp[njetmax];
   Float_t         FatJet_mass_jesPileUpEnvelopeUp[njetmax];
   Float_t         FatJet_msoftdrop_jesPileUpEnvelopeUp[njetmax];
   Float_t         FatJet_pt_jesSubTotalPileUpUp[njetmax];
   Float_t         FatJet_mass_jesSubTotalPileUpUp[njetmax];
   Float_t         FatJet_msoftdrop_jesSubTotalPileUpUp[njetmax];
   Float_t         FatJet_pt_jesSubTotalRelativeUp[njetmax];
   Float_t         FatJet_mass_jesSubTotalRelativeUp[njetmax];
   Float_t         FatJet_msoftdrop_jesSubTotalRelativeUp[njetmax];
   Float_t         FatJet_pt_jesSubTotalPtUp[njetmax];
   Float_t         FatJet_mass_jesSubTotalPtUp[njetmax];
   Float_t         FatJet_msoftdrop_jesSubTotalPtUp[njetmax];
   Float_t         FatJet_pt_jesSubTotalScaleUp[njetmax];
   Float_t         FatJet_mass_jesSubTotalScaleUp[njetmax];
   Float_t         FatJet_msoftdrop_jesSubTotalScaleUp[njetmax];
   Float_t         FatJet_pt_jesSubTotalAbsoluteUp[njetmax];
   Float_t         FatJet_mass_jesSubTotalAbsoluteUp[njetmax];
   Float_t         FatJet_msoftdrop_jesSubTotalAbsoluteUp[njetmax];
   Float_t         FatJet_pt_jesSubTotalMCUp[njetmax];
   Float_t         FatJet_mass_jesSubTotalMCUp[njetmax];
   Float_t         FatJet_msoftdrop_jesSubTotalMCUp[njetmax];
   Float_t         FatJet_pt_jesTotalUp[njetmax];
   Float_t         FatJet_mass_jesTotalUp[njetmax];
   Float_t         FatJet_msoftdrop_jesTotalUp[njetmax];
   Float_t         FatJet_pt_jesTotalNoFlavorUp[njetmax];
   Float_t         FatJet_mass_jesTotalNoFlavorUp[njetmax];
   Float_t         FatJet_msoftdrop_jesTotalNoFlavorUp[njetmax];
   Float_t         FatJet_pt_jesTotalNoTimeUp[njetmax];
   Float_t         FatJet_mass_jesTotalNoTimeUp[njetmax];
   Float_t         FatJet_msoftdrop_jesTotalNoTimeUp[njetmax];
   Float_t         FatJet_pt_jesTotalNoFlavorNoTimeUp[njetmax];
   Float_t         FatJet_mass_jesTotalNoFlavorNoTimeUp[njetmax];
   Float_t         FatJet_msoftdrop_jesTotalNoFlavorNoTimeUp[njetmax];
   Float_t         FatJet_pt_jesFlavorZJetUp[njetmax];
   Float_t         FatJet_mass_jesFlavorZJetUp[njetmax];
   Float_t         FatJet_msoftdrop_jesFlavorZJetUp[njetmax];
   Float_t         FatJet_pt_jesFlavorPhotonJetUp[njetmax];
   Float_t         FatJet_mass_jesFlavorPhotonJetUp[njetmax];
   Float_t         FatJet_msoftdrop_jesFlavorPhotonJetUp[njetmax];
   Float_t         FatJet_pt_jesFlavorPureGluonUp[njetmax];
   Float_t         FatJet_mass_jesFlavorPureGluonUp[njetmax];
   Float_t         FatJet_msoftdrop_jesFlavorPureGluonUp[njetmax];
   Float_t         FatJet_pt_jesFlavorPureQuarkUp[njetmax];
   Float_t         FatJet_mass_jesFlavorPureQuarkUp[njetmax];
   Float_t         FatJet_msoftdrop_jesFlavorPureQuarkUp[njetmax];
   Float_t         FatJet_pt_jesFlavorPureCharmUp[njetmax];
   Float_t         FatJet_mass_jesFlavorPureCharmUp[njetmax];
   Float_t         FatJet_msoftdrop_jesFlavorPureCharmUp[njetmax];
   Float_t         FatJet_pt_jesFlavorPureBottomUp[njetmax];
   Float_t         FatJet_mass_jesFlavorPureBottomUp[njetmax];
   Float_t         FatJet_msoftdrop_jesFlavorPureBottomUp[njetmax];
   Float_t         FatJet_pt_jesTimeRunBUp[njetmax];
   Float_t         FatJet_mass_jesTimeRunBUp[njetmax];
   Float_t         FatJet_msoftdrop_jesTimeRunBUp[njetmax];
   Float_t         FatJet_pt_jesTimeRunCUp[njetmax];
   Float_t         FatJet_mass_jesTimeRunCUp[njetmax];
   Float_t         FatJet_msoftdrop_jesTimeRunCUp[njetmax];
   Float_t         FatJet_pt_jesTimeRunDEUp[njetmax];
   Float_t         FatJet_mass_jesTimeRunDEUp[njetmax];
   Float_t         FatJet_msoftdrop_jesTimeRunDEUp[njetmax];
   Float_t         FatJet_pt_jesTimeRunFUp[njetmax];
   Float_t         FatJet_mass_jesTimeRunFUp[njetmax];
   Float_t         FatJet_msoftdrop_jesTimeRunFUp[njetmax];
   Float_t         FatJet_pt_jesCorrelationGroupMPFInSituUp[njetmax];
   Float_t         FatJet_mass_jesCorrelationGroupMPFInSituUp[njetmax];
   Float_t         FatJet_msoftdrop_jesCorrelationGroupMPFInSituUp[njetmax];
   Float_t         FatJet_pt_jesCorrelationGroupIntercalibrationUp[njetmax];
   Float_t         FatJet_mass_jesCorrelationGroupIntercalibrationUp[njetmax];
   Float_t         FatJet_msoftdrop_jesCorrelationGroupIntercalibrationUp[njetmax];
   Float_t         FatJet_pt_jesCorrelationGroupbJESUp[njetmax];
   Float_t         FatJet_mass_jesCorrelationGroupbJESUp[njetmax];
   Float_t         FatJet_msoftdrop_jesCorrelationGroupbJESUp[njetmax];
   Float_t         FatJet_pt_jesCorrelationGroupFlavorUp[njetmax];
   Float_t         FatJet_mass_jesCorrelationGroupFlavorUp[njetmax];
   Float_t         FatJet_msoftdrop_jesCorrelationGroupFlavorUp[njetmax];
   Float_t         FatJet_pt_jesCorrelationGroupUncorrelatedUp[njetmax];
   Float_t         FatJet_mass_jesCorrelationGroupUncorrelatedUp[njetmax];
   Float_t         FatJet_msoftdrop_jesCorrelationGroupUncorrelatedUp[njetmax];
   Float_t         FatJet_pt_jerDown[njetmax];
   Float_t         FatJet_mass_jerDown[njetmax];
   Float_t         FatJet_mass_jmrDown[njetmax];
   Float_t         FatJet_mass_jmsDown[njetmax];
   Float_t         FatJet_msoftdrop_jerDown[njetmax];
   Float_t         FatJet_msoftdrop_jmrDown[njetmax];
   Float_t         FatJet_msoftdrop_jmsDown[njetmax];
   Float_t         FatJet_msoftdrop_tau21DDT_jerDown[njetmax];
   Float_t         FatJet_msoftdrop_tau21DDT_jmrDown[njetmax];
   Float_t         FatJet_msoftdrop_tau21DDT_jmsDown[njetmax];
   Float_t         FatJet_pt_jesAbsoluteStatDown[njetmax];
   Float_t         FatJet_mass_jesAbsoluteStatDown[njetmax];
   Float_t         FatJet_msoftdrop_jesAbsoluteStatDown[njetmax];
   Float_t         FatJet_pt_jesAbsoluteScaleDown[njetmax];
   Float_t         FatJet_mass_jesAbsoluteScaleDown[njetmax];
   Float_t         FatJet_msoftdrop_jesAbsoluteScaleDown[njetmax];
   Float_t         FatJet_pt_jesAbsoluteFlavMapDown[njetmax];
   Float_t         FatJet_mass_jesAbsoluteFlavMapDown[njetmax];
   Float_t         FatJet_msoftdrop_jesAbsoluteFlavMapDown[njetmax];
   Float_t         FatJet_pt_jesAbsoluteMPFBiasDown[njetmax];
   Float_t         FatJet_mass_jesAbsoluteMPFBiasDown[njetmax];
   Float_t         FatJet_msoftdrop_jesAbsoluteMPFBiasDown[njetmax];
   Float_t         FatJet_pt_jesFragmentationDown[njetmax];
   Float_t         FatJet_mass_jesFragmentationDown[njetmax];
   Float_t         FatJet_msoftdrop_jesFragmentationDown[njetmax];
   Float_t         FatJet_pt_jesSinglePionECALDown[njetmax];
   Float_t         FatJet_mass_jesSinglePionECALDown[njetmax];
   Float_t         FatJet_msoftdrop_jesSinglePionECALDown[njetmax];
   Float_t         FatJet_pt_jesSinglePionHCALDown[njetmax];
   Float_t         FatJet_mass_jesSinglePionHCALDown[njetmax];
   Float_t         FatJet_msoftdrop_jesSinglePionHCALDown[njetmax];
   Float_t         FatJet_pt_jesFlavorQCDDown[njetmax];
   Float_t         FatJet_mass_jesFlavorQCDDown[njetmax];
   Float_t         FatJet_msoftdrop_jesFlavorQCDDown[njetmax];
   Float_t         FatJet_pt_jesTimePtEtaDown[njetmax];
   Float_t         FatJet_mass_jesTimePtEtaDown[njetmax];
   Float_t         FatJet_msoftdrop_jesTimePtEtaDown[njetmax];
   Float_t         FatJet_pt_jesRelativeJEREC1Down[njetmax];
   Float_t         FatJet_mass_jesRelativeJEREC1Down[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativeJEREC1Down[njetmax];
   Float_t         FatJet_pt_jesRelativeJEREC2Down[njetmax];
   Float_t         FatJet_mass_jesRelativeJEREC2Down[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativeJEREC2Down[njetmax];
   Float_t         FatJet_pt_jesRelativeJERHFDown[njetmax];
   Float_t         FatJet_mass_jesRelativeJERHFDown[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativeJERHFDown[njetmax];
   Float_t         FatJet_pt_jesRelativePtBBDown[njetmax];
   Float_t         FatJet_mass_jesRelativePtBBDown[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativePtBBDown[njetmax];
   Float_t         FatJet_pt_jesRelativePtEC1Down[njetmax];
   Float_t         FatJet_mass_jesRelativePtEC1Down[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativePtEC1Down[njetmax];
   Float_t         FatJet_pt_jesRelativePtEC2Down[njetmax];
   Float_t         FatJet_mass_jesRelativePtEC2Down[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativePtEC2Down[njetmax];
   Float_t         FatJet_pt_jesRelativePtHFDown[njetmax];
   Float_t         FatJet_mass_jesRelativePtHFDown[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativePtHFDown[njetmax];
   Float_t         FatJet_pt_jesRelativeBalDown[njetmax];
   Float_t         FatJet_mass_jesRelativeBalDown[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativeBalDown[njetmax];
   Float_t         FatJet_pt_jesRelativeSampleDown[njetmax];
   Float_t         FatJet_mass_jesRelativeSampleDown[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativeSampleDown[njetmax];
   Float_t         FatJet_pt_jesRelativeFSRDown[njetmax];
   Float_t         FatJet_mass_jesRelativeFSRDown[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativeFSRDown[njetmax];
   Float_t         FatJet_pt_jesRelativeStatFSRDown[njetmax];
   Float_t         FatJet_mass_jesRelativeStatFSRDown[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativeStatFSRDown[njetmax];
   Float_t         FatJet_pt_jesRelativeStatECDown[njetmax];
   Float_t         FatJet_mass_jesRelativeStatECDown[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativeStatECDown[njetmax];
   Float_t         FatJet_pt_jesRelativeStatHFDown[njetmax];
   Float_t         FatJet_mass_jesRelativeStatHFDown[njetmax];
   Float_t         FatJet_msoftdrop_jesRelativeStatHFDown[njetmax];
   Float_t         FatJet_pt_jesPileUpDataMCDown[njetmax];
   Float_t         FatJet_mass_jesPileUpDataMCDown[njetmax];
   Float_t         FatJet_msoftdrop_jesPileUpDataMCDown[njetmax];
   Float_t         FatJet_pt_jesPileUpPtRefDown[njetmax];
   Float_t         FatJet_mass_jesPileUpPtRefDown[njetmax];
   Float_t         FatJet_msoftdrop_jesPileUpPtRefDown[njetmax];
   Float_t         FatJet_pt_jesPileUpPtBBDown[njetmax];
   Float_t         FatJet_mass_jesPileUpPtBBDown[njetmax];
   Float_t         FatJet_msoftdrop_jesPileUpPtBBDown[njetmax];
   Float_t         FatJet_pt_jesPileUpPtEC1Down[njetmax];
   Float_t         FatJet_mass_jesPileUpPtEC1Down[njetmax];
   Float_t         FatJet_msoftdrop_jesPileUpPtEC1Down[njetmax];
   Float_t         FatJet_pt_jesPileUpPtEC2Down[njetmax];
   Float_t         FatJet_mass_jesPileUpPtEC2Down[njetmax];
   Float_t         FatJet_msoftdrop_jesPileUpPtEC2Down[njetmax];
   Float_t         FatJet_pt_jesPileUpPtHFDown[njetmax];
   Float_t         FatJet_mass_jesPileUpPtHFDown[njetmax];
   Float_t         FatJet_msoftdrop_jesPileUpPtHFDown[njetmax];
   Float_t         FatJet_pt_jesPileUpMuZeroDown[njetmax];
   Float_t         FatJet_mass_jesPileUpMuZeroDown[njetmax];
   Float_t         FatJet_msoftdrop_jesPileUpMuZeroDown[njetmax];
   Float_t         FatJet_pt_jesPileUpEnvelopeDown[njetmax];
   Float_t         FatJet_mass_jesPileUpEnvelopeDown[njetmax];
   Float_t         FatJet_msoftdrop_jesPileUpEnvelopeDown[njetmax];
   Float_t         FatJet_pt_jesSubTotalPileUpDown[njetmax];
   Float_t         FatJet_mass_jesSubTotalPileUpDown[njetmax];
   Float_t         FatJet_msoftdrop_jesSubTotalPileUpDown[njetmax];
   Float_t         FatJet_pt_jesSubTotalRelativeDown[njetmax];
   Float_t         FatJet_mass_jesSubTotalRelativeDown[njetmax];
   Float_t         FatJet_msoftdrop_jesSubTotalRelativeDown[njetmax];
   Float_t         FatJet_pt_jesSubTotalPtDown[njetmax];
   Float_t         FatJet_mass_jesSubTotalPtDown[njetmax];
   Float_t         FatJet_msoftdrop_jesSubTotalPtDown[njetmax];
   Float_t         FatJet_pt_jesSubTotalScaleDown[njetmax];
   Float_t         FatJet_mass_jesSubTotalScaleDown[njetmax];
   Float_t         FatJet_msoftdrop_jesSubTotalScaleDown[njetmax];
   Float_t         FatJet_pt_jesSubTotalAbsoluteDown[njetmax];
   Float_t         FatJet_mass_jesSubTotalAbsoluteDown[njetmax];
   Float_t         FatJet_msoftdrop_jesSubTotalAbsoluteDown[njetmax];
   Float_t         FatJet_pt_jesSubTotalMCDown[njetmax];
   Float_t         FatJet_mass_jesSubTotalMCDown[njetmax];
   Float_t         FatJet_msoftdrop_jesSubTotalMCDown[njetmax];
   Float_t         FatJet_pt_jesTotalDown[njetmax];
   Float_t         FatJet_mass_jesTotalDown[njetmax];
   Float_t         FatJet_msoftdrop_jesTotalDown[njetmax];
   Float_t         FatJet_pt_jesTotalNoFlavorDown[njetmax];
   Float_t         FatJet_mass_jesTotalNoFlavorDown[njetmax];
   Float_t         FatJet_msoftdrop_jesTotalNoFlavorDown[njetmax];
   Float_t         FatJet_pt_jesTotalNoTimeDown[njetmax];
   Float_t         FatJet_mass_jesTotalNoTimeDown[njetmax];
   Float_t         FatJet_msoftdrop_jesTotalNoTimeDown[njetmax];
   Float_t         FatJet_pt_jesTotalNoFlavorNoTimeDown[njetmax];
   Float_t         FatJet_mass_jesTotalNoFlavorNoTimeDown[njetmax];
   Float_t         FatJet_msoftdrop_jesTotalNoFlavorNoTimeDown[njetmax];
   Float_t         FatJet_pt_jesFlavorZJetDown[njetmax];
   Float_t         FatJet_mass_jesFlavorZJetDown[njetmax];
   Float_t         FatJet_msoftdrop_jesFlavorZJetDown[njetmax];
   Float_t         FatJet_pt_jesFlavorPhotonJetDown[njetmax];
   Float_t         FatJet_mass_jesFlavorPhotonJetDown[njetmax];
   Float_t         FatJet_msoftdrop_jesFlavorPhotonJetDown[njetmax];
   Float_t         FatJet_pt_jesFlavorPureGluonDown[njetmax];
   Float_t         FatJet_mass_jesFlavorPureGluonDown[njetmax];
   Float_t         FatJet_msoftdrop_jesFlavorPureGluonDown[njetmax];
   Float_t         FatJet_pt_jesFlavorPureQuarkDown[njetmax];
   Float_t         FatJet_mass_jesFlavorPureQuarkDown[njetmax];
   Float_t         FatJet_msoftdrop_jesFlavorPureQuarkDown[njetmax];
   Float_t         FatJet_pt_jesFlavorPureCharmDown[njetmax];
   Float_t         FatJet_mass_jesFlavorPureCharmDown[njetmax];
   Float_t         FatJet_msoftdrop_jesFlavorPureCharmDown[njetmax];
   Float_t         FatJet_pt_jesFlavorPureBottomDown[njetmax];
   Float_t         FatJet_mass_jesFlavorPureBottomDown[njetmax];
   Float_t         FatJet_msoftdrop_jesFlavorPureBottomDown[njetmax];
   Float_t         FatJet_pt_jesTimeRunBDown[njetmax];
   Float_t         FatJet_mass_jesTimeRunBDown[njetmax];
   Float_t         FatJet_msoftdrop_jesTimeRunBDown[njetmax];
   Float_t         FatJet_pt_jesTimeRunCDown[njetmax];
   Float_t         FatJet_mass_jesTimeRunCDown[njetmax];
   Float_t         FatJet_msoftdrop_jesTimeRunCDown[njetmax];
   Float_t         FatJet_pt_jesTimeRunDEDown[njetmax];
   Float_t         FatJet_mass_jesTimeRunDEDown[njetmax];
   Float_t         FatJet_msoftdrop_jesTimeRunDEDown[njetmax];
   Float_t         FatJet_pt_jesTimeRunFDown[njetmax];
   Float_t         FatJet_mass_jesTimeRunFDown[njetmax];
   Float_t         FatJet_msoftdrop_jesTimeRunFDown[njetmax];
   Float_t         FatJet_pt_jesCorrelationGroupMPFInSituDown[njetmax];
   Float_t         FatJet_mass_jesCorrelationGroupMPFInSituDown[njetmax];
   Float_t         FatJet_msoftdrop_jesCorrelationGroupMPFInSituDown[njetmax];
   Float_t         FatJet_pt_jesCorrelationGroupIntercalibrationDown[njetmax];
   Float_t         FatJet_mass_jesCorrelationGroupIntercalibrationDown[njetmax];
   Float_t         FatJet_msoftdrop_jesCorrelationGroupIntercalibrationDown[njetmax];
   Float_t         FatJet_pt_jesCorrelationGroupbJESDown[njetmax];
   Float_t         FatJet_mass_jesCorrelationGroupbJESDown[njetmax];
   Float_t         FatJet_msoftdrop_jesCorrelationGroupbJESDown[njetmax];
   Float_t         FatJet_pt_jesCorrelationGroupFlavorDown[njetmax];
   Float_t         FatJet_mass_jesCorrelationGroupFlavorDown[njetmax];
   Float_t         FatJet_msoftdrop_jesCorrelationGroupFlavorDown[njetmax];
   Float_t         FatJet_pt_jesCorrelationGroupUncorrelatedDown[njetmax];
   Float_t         FatJet_mass_jesCorrelationGroupUncorrelatedDown[njetmax];
   Float_t         FatJet_msoftdrop_jesCorrelationGroupUncorrelatedDown[njetmax];
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
   Float_t         Jet_deepflavbtagSF[njetmax];
   Float_t         Jet_deepflavbtagSF_up[njetmax];
   Float_t         Jet_deepflavbtagSF_down[njetmax];
   Float_t         Jet_deepflavbtagSF_shape[njetmax];
   Float_t         Jet_deepflavbtagSF_shape_up_jes[njetmax];
   Float_t         Jet_deepflavbtagSF_shape_down_jes[njetmax];
   Float_t         Jet_deepflavbtagSF_shape_up_lf[njetmax];
   Float_t         Jet_deepflavbtagSF_shape_down_lf[njetmax];
   Float_t         Jet_deepflavbtagSF_shape_up_hf[njetmax];
   Float_t         Jet_deepflavbtagSF_shape_down_hf[njetmax];
   Float_t         Jet_deepflavbtagSF_shape_up_hfstats1[njetmax];
   Float_t         Jet_deepflavbtagSF_shape_down_hfstats1[njetmax];
   Float_t         Jet_deepflavbtagSF_shape_up_hfstats2[njetmax];
   Float_t         Jet_deepflavbtagSF_shape_down_hfstats2[njetmax];
   Float_t         Jet_deepflavbtagSF_shape_up_lfstats1[njetmax];
   Float_t         Jet_deepflavbtagSF_shape_down_lfstats1[njetmax];
   Float_t         Jet_deepflavbtagSF_shape_up_lfstats2[njetmax];
   Float_t         Jet_deepflavbtagSF_shape_down_lfstats2[njetmax];
   Float_t         Jet_deepflavbtagSF_shape_up_cferr1[njetmax];
   Float_t         Jet_deepflavbtagSF_shape_down_cferr1[njetmax];
   Float_t         Jet_deepflavbtagSF_shape_up_cferr2[njetmax];
   Float_t         Jet_deepflavbtagSF_shape_down_cferr2[njetmax];
   Float_t         puWeight;
   Float_t         puWeightUp;
   Float_t         puWeightDown;
   Float_t         PrefireWeight;
   Float_t         PrefireWeight_Up;
   Float_t         PrefireWeight_Down;
   Float_t		   Event_Ht;
   bool			   had_trig;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_HTXS_Higgs_pt;   //!
   TBranch        *b_HTXS_Higgs_y;   //!
   TBranch        *b_HTXS_stage1_1_cat_pTjet25GeV;   //!
   TBranch        *b_HTXS_stage1_1_cat_pTjet30GeV;   //!
   TBranch        *b_HTXS_stage1_1_fine_cat_pTjet25GeV;   //!
   TBranch        *b_HTXS_stage1_1_fine_cat_pTjet30GeV;   //!
   TBranch        *b_HTXS_stage_0;   //!
   TBranch        *b_HTXS_stage_1_pTjet25;   //!
   TBranch        *b_HTXS_stage_1_pTjet30;   //!
   TBranch        *b_HTXS_njets25;   //!
   TBranch        *b_HTXS_njets30;   //!
   TBranch        *b_btagWeight_CSVV2;   //!
   TBranch        *b_btagWeight_DeepCSVB;   //!
   TBranch        *b_CaloMET_phi;   //!
   TBranch        *b_CaloMET_pt;   //!
   TBranch        *b_CaloMET_sumEt;   //!
   TBranch        *b_ChsMET_phi;   //!
   TBranch        *b_ChsMET_pt;   //!
   TBranch        *b_ChsMET_sumEt;   //!
   TBranch        *b_nCorrT1METJet;   //!
   TBranch        *b_CorrT1METJet_area;   //!
   TBranch        *b_CorrT1METJet_eta;   //!
   TBranch        *b_CorrT1METJet_muonSubtrFactor;   //!
   TBranch        *b_CorrT1METJet_phi;   //!
   TBranch        *b_CorrT1METJet_rawPt;   //!
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
   TBranch        *b_Electron_jetPtRelv2;   //!
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
   TBranch        *b_Electron_vidNestedWPBitmapHEEP;   //!
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
   TBranch        *b_Electron_seedGain;   //!
   TBranch        *b_Flag_ecalBadCalibFilterV2;   //!
   TBranch        *b_nFatJet;   //!
   TBranch        *b_FatJet_area;   //!
   TBranch        *b_FatJet_btagCMVA;   //!
   TBranch        *b_FatJet_btagCSVV2;   //!
   TBranch        *b_FatJet_btagDDBvL;   //!
   TBranch        *b_FatJet_btagDDBvL_noMD;   //!
   TBranch        *b_FatJet_btagDDCvB;   //!
   TBranch        *b_FatJet_btagDDCvB_noMD;   //!
   TBranch        *b_FatJet_btagDDCvL;   //!
   TBranch        *b_FatJet_btagDDCvL_noMD;   //!
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
   TBranch        *b_nFsrPhoton;   //!
   TBranch        *b_FsrPhoton_dROverEt2;   //!
   TBranch        *b_FsrPhoton_eta;   //!
   TBranch        *b_FsrPhoton_phi;   //!
   TBranch        *b_FsrPhoton_pt;   //!
   TBranch        *b_FsrPhoton_relIso03;   //!
   TBranch        *b_FsrPhoton_muonIdx;   //!
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
   TBranch        *b_nLHEReweightingWeight;   //!
   TBranch        *b_LHEReweightingWeight;   //!
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
   TBranch        *b_Jet_btagDeepFlavC;   //!
   TBranch        *b_Jet_chEmEF;   //!
   TBranch        *b_Jet_chHEF;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_jercCHF;   //!
   TBranch        *b_Jet_jercCHPUF;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_muEF;   //!
   TBranch        *b_Jet_muonSubtrFactor;   //!
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
   TBranch        *b_L1PreFiringWeight_Dn;   //!
   TBranch        *b_L1PreFiringWeight_Nom;   //!
   TBranch        *b_L1PreFiringWeight_Up;   //!
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
   TBranch        *b_METFixEE2017_significance;   //!
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
   TBranch        *b_Muon_tkRelIso;   //!
   TBranch        *b_Muon_tunepRelPt;   //!
   TBranch        *b_Muon_mvaLowPt;   //!
   TBranch        *b_Muon_mvaTTH;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_jetIdx;   //!
   TBranch        *b_Muon_nStations;   //!
   TBranch        *b_Muon_nTrackerLayers;   //!
   TBranch        *b_Muon_pdgId;   //!
   TBranch        *b_Muon_tightCharge;   //!
   TBranch        *b_Muon_fsrPhotonIdx;   //!
   TBranch        *b_Muon_highPtId;   //!
   TBranch        *b_Muon_inTimeMuon;   //!
   TBranch        *b_Muon_isGlobal;   //!
   TBranch        *b_Muon_isPFcand;   //!
   TBranch        *b_Muon_isTracker;   //!
   TBranch        *b_Muon_looseId;   //!
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
   TBranch        *b_Photon_seedGain;   //!
   TBranch        *b_Pileup_nTrueInt;   //!
   TBranch        *b_Pileup_pudensity;   //!
   TBranch        *b_Pileup_gpudensity;   //!
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
   TBranch        *b_fixedGridRhoFastjetCentral;   //!
   TBranch        *b_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_fixedGridRhoFastjetCentralChargedPileUp;   //!
   TBranch        *b_fixedGridRhoFastjetCentralNeutral;   //!
   TBranch        *b_nGenDressedLepton;   //!
   TBranch        *b_GenDressedLepton_eta;   //!
   TBranch        *b_GenDressedLepton_mass;   //!
   TBranch        *b_GenDressedLepton_phi;   //!
   TBranch        *b_GenDressedLepton_pt;   //!
   TBranch        *b_GenDressedLepton_pdgId;   //!
   TBranch        *b_GenDressedLepton_hasTauAnc;   //!
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
   TBranch        *b_Tau_rawAntiEle2018;   //!
   TBranch        *b_Tau_rawDeepTau2017v2p1VSe;   //!
   TBranch        *b_Tau_rawDeepTau2017v2p1VSjet;   //!
   TBranch        *b_Tau_rawDeepTau2017v2p1VSmu;   //!
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
   TBranch        *b_Tau_rawAntiEleCat2018;   //!
   TBranch        *b_Tau_idAntiEle;   //!
   TBranch        *b_Tau_idAntiEle2018;   //!
   TBranch        *b_Tau_idAntiMu;   //!
   TBranch        *b_Tau_idDecayMode;   //!
   TBranch        *b_Tau_idDecayModeNewDMs;   //!
   TBranch        *b_Tau_idDeepTau2017v2p1VSe;   //!
   TBranch        *b_Tau_idDeepTau2017v2p1VSjet;   //!
   TBranch        *b_Tau_idDeepTau2017v2p1VSmu;   //!
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
   TBranch        *b_SV_dxy;   //!
   TBranch        *b_SV_dxySig;   //!
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
   TBranch        *b_HLT_AK8PFHT750_TrimMass50;   //!
   TBranch        *b_HLT_AK8PFHT800_TrimMass50;   //!
   TBranch        *b_HLT_AK8PFHT850_TrimMass50;   //!
   TBranch        *b_HLT_AK8PFHT900_TrimMass50;   //!
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
   TBranch        *b_L1_AlwaysTrue;   //!
   TBranch        *b_L1_BPTX_AND_Ref1_VME;   //!
   TBranch        *b_L1_BPTX_AND_Ref3_VME;   //!
   TBranch        *b_L1_BPTX_AND_Ref4_VME;   //!
   TBranch        *b_L1_BPTX_BeamGas_B1_VME;   //!
   TBranch        *b_L1_BPTX_BeamGas_B2_VME;   //!
   TBranch        *b_L1_BPTX_BeamGas_Ref1_VME;   //!
   TBranch        *b_L1_BPTX_BeamGas_Ref2_VME;   //!
   TBranch        *b_L1_BPTX_NotOR_VME;   //!
   TBranch        *b_L1_BPTX_OR_Ref3_VME;   //!
   TBranch        *b_L1_BPTX_OR_Ref4_VME;   //!
   TBranch        *b_L1_BPTX_RefAND_VME;   //!
   TBranch        *b_L1_BptxMinus;   //!
   TBranch        *b_L1_BptxOR;   //!
   TBranch        *b_L1_BptxPlus;   //!
   TBranch        *b_L1_BptxXOR;   //!
   TBranch        *b_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;   //!
   TBranch        *b_L1_DoubleEG6_HTT240er;   //!
   TBranch        *b_L1_DoubleEG6_HTT250er;   //!
   TBranch        *b_L1_DoubleEG6_HTT255er;   //!
   TBranch        *b_L1_DoubleEG6_HTT270er;   //!
   TBranch        *b_L1_DoubleEG6_HTT300er;   //!
   TBranch        *b_L1_DoubleEG8er2p6_HTT255er;   //!
   TBranch        *b_L1_DoubleEG8er2p6_HTT270er;   //!
   TBranch        *b_L1_DoubleEG8er2p6_HTT300er;   //!
   TBranch        *b_L1_DoubleEG_15_10;   //!
   TBranch        *b_L1_DoubleEG_18_17;   //!
   TBranch        *b_L1_DoubleEG_20_18;   //!
   TBranch        *b_L1_DoubleEG_22_10;   //!
   TBranch        *b_L1_DoubleEG_22_12;   //!
   TBranch        *b_L1_DoubleEG_22_15;   //!
   TBranch        *b_L1_DoubleEG_23_10;   //!
   TBranch        *b_L1_DoubleEG_24_17;   //!
   TBranch        *b_L1_DoubleEG_25_12;   //!
   TBranch        *b_L1_DoubleEG_25_13;   //!
   TBranch        *b_L1_DoubleEG_25_14;   //!
   TBranch        *b_L1_DoubleEG_LooseIso23_10;   //!
   TBranch        *b_L1_DoubleEG_LooseIso24_10;   //!
   TBranch        *b_L1_DoubleIsoTau28er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau30er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau32er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau33er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau34er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau35er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau36er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau38er2p1;   //!
   TBranch        *b_L1_DoubleJet100er2p3_dEta_Max1p6;   //!
   TBranch        *b_L1_DoubleJet100er2p7;   //!
   TBranch        *b_L1_DoubleJet112er2p3_dEta_Max1p6;   //!
   TBranch        *b_L1_DoubleJet112er2p7;   //!
   TBranch        *b_L1_DoubleJet120er2p7;   //!
   TBranch        *b_L1_DoubleJet150er2p7;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min300_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min320_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min340_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min360_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min380_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min400_Mu10;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min400_Mu6;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min400_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet35_rmovlp_IsoTau45_Mass_Min450;   //!
   TBranch        *b_L1_DoubleJet40er2p7;   //!
   TBranch        *b_L1_DoubleJet50er2p7;   //!
   TBranch        *b_L1_DoubleJet60er2p7;   //!
   TBranch        *b_L1_DoubleJet60er2p7_ETM100;   //!
   TBranch        *b_L1_DoubleJet60er2p7_ETM60;   //!
   TBranch        *b_L1_DoubleJet60er2p7_ETM70;   //!
   TBranch        *b_L1_DoubleJet60er2p7_ETM80;   //!
   TBranch        *b_L1_DoubleJet60er2p7_ETM90;   //!
   TBranch        *b_L1_DoubleJet80er2p7;   //!
   TBranch        *b_L1_DoubleJet_100_30_DoubleJet30_Mass_Min620;   //!
   TBranch        *b_L1_DoubleJet_100_35_DoubleJet35_Mass_Min620;   //!
   TBranch        *b_L1_DoubleJet_110_35_DoubleJet35_Mass_Min620;   //!
   TBranch        *b_L1_DoubleJet_110_40_DoubleJet40_Mass_Min620;   //!
   TBranch        *b_L1_DoubleJet_115_35_DoubleJet35_Mass_Min620;   //!
   TBranch        *b_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620;   //!
   TBranch        *b_L1_DoubleJet_90_30_DoubleJet30_Mass_Min620;   //!
   TBranch        *b_L1_DoubleLooseIsoEG22er2p1;   //!
   TBranch        *b_L1_DoubleLooseIsoEG24er2p1;   //!
   TBranch        *b_L1_DoubleMu0;   //!
   TBranch        *b_L1_DoubleMu0_ETM40;   //!
   TBranch        *b_L1_DoubleMu0_ETM55;   //!
   TBranch        *b_L1_DoubleMu0_ETM60;   //!
   TBranch        *b_L1_DoubleMu0_ETM65;   //!
   TBranch        *b_L1_DoubleMu0_ETM70;   //!
   TBranch        *b_L1_DoubleMu0_SQ;   //!
   TBranch        *b_L1_DoubleMu0_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4;   //!
   TBranch        *b_L1_DoubleMu0er1p4_dEta_Max1p8_OS;   //!
   TBranch        *b_L1_DoubleMu0er1p5_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4;   //!
   TBranch        *b_L1_DoubleMu0er1p5_SQ_dR_Max1p4;   //!
   TBranch        *b_L1_DoubleMu0er2_SQ_dR_Max1p4;   //!
   TBranch        *b_L1_DoubleMu18er2p1;   //!
   TBranch        *b_L1_DoubleMu22er2p1;   //!
   TBranch        *b_L1_DoubleMu3_OS_DoubleEG7p5Upsilon;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF40_Jet60_OR_DoubleJet30;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF50_Jet60_OR_DoubleJet30;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF60_Jet60_OR_DoubleJet30;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF70_Jet60_OR_DoubleJet30;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF80_Jet60_OR_DoubleJet30;   //!
   TBranch        *b_L1_DoubleMu3_SQ_HTT100er;   //!
   TBranch        *b_L1_DoubleMu3_SQ_HTT200er;   //!
   TBranch        *b_L1_DoubleMu3_SQ_HTT220er;   //!
   TBranch        *b_L1_DoubleMu3_SQ_HTT240er;   //!
   TBranch        *b_L1_DoubleMu4_OS_EG12;   //!
   TBranch        *b_L1_DoubleMu4_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu4_SQ_OS_dR_Max1p2;   //!
   TBranch        *b_L1_DoubleMu4p5_SQ;   //!
   TBranch        *b_L1_DoubleMu4p5_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu4p5_SQ_OS_dR_Max1p2;   //!
   TBranch        *b_L1_DoubleMu4p5er2p0_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18;   //!
   TBranch        *b_L1_DoubleMu5Upsilon_OS_DoubleEG3;   //!
   TBranch        *b_L1_DoubleMu5_OS_EG12;   //!
   TBranch        *b_L1_DoubleMu5_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu5_SQ_OS_Mass7to18;   //!
   TBranch        *b_L1_DoubleMu6_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu7_EG7;   //!
   TBranch        *b_L1_DoubleMu7_SQ_EG7;   //!
   TBranch        *b_L1_DoubleMu8_SQ;   //!
   TBranch        *b_L1_DoubleMu_10_0_dEta_Max1p8;   //!
   TBranch        *b_L1_DoubleMu_11_4;   //!
   TBranch        *b_L1_DoubleMu_12_5;   //!
   TBranch        *b_L1_DoubleMu_12_8;   //!
   TBranch        *b_L1_DoubleMu_13_6;   //!
   TBranch        *b_L1_DoubleMu_15_5;   //!
   TBranch        *b_L1_DoubleMu_15_5_SQ;   //!
   TBranch        *b_L1_DoubleMu_15_7;   //!
   TBranch        *b_L1_DoubleMu_15_7_SQ;   //!
   TBranch        *b_L1_DoubleMu_15_7_SQ_Mass_Min4;   //!
   TBranch        *b_L1_DoubleMu_20_2_SQ_Mass_Max20;   //!
   TBranch        *b_L1_DoubleTau50er2p1;   //!
   TBranch        *b_L1_DoubleTau70er2p1;   //!
   TBranch        *b_L1_EG25er2p1_HTT125er;   //!
   TBranch        *b_L1_EG27er2p1_HTT200er;   //!
   TBranch        *b_L1_ETM100;   //!
   TBranch        *b_L1_ETM100_Jet60_dPhi_Min0p4;   //!
   TBranch        *b_L1_ETM105;   //!
   TBranch        *b_L1_ETM110;   //!
   TBranch        *b_L1_ETM110_Jet60_dPhi_Min0p4;   //!
   TBranch        *b_L1_ETM115;   //!
   TBranch        *b_L1_ETM120;   //!
   TBranch        *b_L1_ETM150;   //!
   TBranch        *b_L1_ETM30;   //!
   TBranch        *b_L1_ETM40;   //!
   TBranch        *b_L1_ETM50;   //!
   TBranch        *b_L1_ETM60;   //!
   TBranch        *b_L1_ETM70;   //!
   TBranch        *b_L1_ETM75;   //!
   TBranch        *b_L1_ETM75_Jet60_dPhi_Min0p4;   //!
   TBranch        *b_L1_ETM80;   //!
   TBranch        *b_L1_ETM80_Jet60_dPhi_Min0p4;   //!
   TBranch        *b_L1_ETM85;   //!
   TBranch        *b_L1_ETM90;   //!
   TBranch        *b_L1_ETM90_Jet60_dPhi_Min0p4;   //!
   TBranch        *b_L1_ETM95;   //!
   TBranch        *b_L1_ETMHF100;   //!
   TBranch        *b_L1_ETMHF100_HTT60er;   //!
   TBranch        *b_L1_ETMHF100_Jet60_OR_DiJet30woTT28;   //!
   TBranch        *b_L1_ETMHF100_Jet60_OR_DoubleJet30;   //!
   TBranch        *b_L1_ETMHF100_Jet90_OR_DoubleJet45_OR_TripleJet30;   //!
   TBranch        *b_L1_ETMHF110;   //!
   TBranch        *b_L1_ETMHF110_HTT60er;   //!
   TBranch        *b_L1_ETMHF110_Jet60_OR_DiJet30woTT28;   //!
   TBranch        *b_L1_ETMHF110_Jet90_OR_DoubleJet45_OR_TripleJet30;   //!
   TBranch        *b_L1_ETMHF120;   //!
   TBranch        *b_L1_ETMHF120_HTT60er;   //!
   TBranch        *b_L1_ETMHF120_Jet60_OR_DiJet30woTT28;   //!
   TBranch        *b_L1_ETMHF150;   //!
   TBranch        *b_L1_ETMHF70;   //!
   TBranch        *b_L1_ETMHF70_Jet90_OR_DoubleJet45_OR_TripleJet30;   //!
   TBranch        *b_L1_ETMHF80;   //!
   TBranch        *b_L1_ETMHF80_HTT60er;   //!
   TBranch        *b_L1_ETMHF80_Jet90_OR_DoubleJet45_OR_TripleJet30;   //!
   TBranch        *b_L1_ETMHF90;   //!
   TBranch        *b_L1_ETMHF90_HTT60er;   //!
   TBranch        *b_L1_ETMHF90_Jet90_OR_DoubleJet45_OR_TripleJet30;   //!
   TBranch        *b_L1_ETT100_BptxAND;   //!
   TBranch        *b_L1_ETT110_BptxAND;   //!
   TBranch        *b_L1_ETT40_BptxAND;   //!
   TBranch        *b_L1_ETT50_BptxAND;   //!
   TBranch        *b_L1_ETT60_BptxAND;   //!
   TBranch        *b_L1_ETT70_BptxAND;   //!
   TBranch        *b_L1_ETT75_BptxAND;   //!
   TBranch        *b_L1_ETT80_BptxAND;   //!
   TBranch        *b_L1_ETT85_BptxAND;   //!
   TBranch        *b_L1_ETT90_BptxAND;   //!
   TBranch        *b_L1_ETT95_BptxAND;   //!
   TBranch        *b_L1_FirstBunchAfterTrain;   //!
   TBranch        *b_L1_FirstBunchInTrain;   //!
   TBranch        *b_L1_FirstCollisionInOrbit;   //!
   TBranch        *b_L1_FirstCollisionInTrain;   //!
   TBranch        *b_L1_HTT120er;   //!
   TBranch        *b_L1_HTT160er;   //!
   TBranch        *b_L1_HTT200er;   //!
   TBranch        *b_L1_HTT220er;   //!
   TBranch        *b_L1_HTT240er;   //!
   TBranch        *b_L1_HTT250er_QuadJet_70_55_40_35_er2p5;   //!
   TBranch        *b_L1_HTT255er;   //!
   TBranch        *b_L1_HTT270er;   //!
   TBranch        *b_L1_HTT280er;   //!
   TBranch        *b_L1_HTT280er_QuadJet_70_55_40_35_er2p5;   //!
   TBranch        *b_L1_HTT300er;   //!
   TBranch        *b_L1_HTT300er_QuadJet_70_55_40_35_er2p5;   //!
   TBranch        *b_L1_HTT320er;   //!
   TBranch        *b_L1_HTT320er_QuadJet_70_55_40_40_er2p4;   //!
   TBranch        *b_L1_HTT320er_QuadJet_70_55_40_40_er2p5;   //!
   TBranch        *b_L1_HTT320er_QuadJet_70_55_45_45_er2p5;   //!
   TBranch        *b_L1_HTT340er;   //!
   TBranch        *b_L1_HTT340er_QuadJet_70_55_40_40_er2p5;   //!
   TBranch        *b_L1_HTT340er_QuadJet_70_55_45_45_er2p5;   //!
   TBranch        *b_L1_HTT380er;   //!
   TBranch        *b_L1_HTT400er;   //!
   TBranch        *b_L1_HTT450er;   //!
   TBranch        *b_L1_HTT500er;   //!
   TBranch        *b_L1_IsoEG33_Mt40;   //!
   TBranch        *b_L1_IsoEG33_Mt44;   //!
   TBranch        *b_L1_IsoEG33_Mt48;   //!
   TBranch        *b_L1_IsoTau40er_ETM100;   //!
   TBranch        *b_L1_IsoTau40er_ETM105;   //!
   TBranch        *b_L1_IsoTau40er_ETM110;   //!
   TBranch        *b_L1_IsoTau40er_ETM115;   //!
   TBranch        *b_L1_IsoTau40er_ETM120;   //!
   TBranch        *b_L1_IsoTau40er_ETM80;   //!
   TBranch        *b_L1_IsoTau40er_ETM85;   //!
   TBranch        *b_L1_IsoTau40er_ETM90;   //!
   TBranch        *b_L1_IsoTau40er_ETM95;   //!
   TBranch        *b_L1_IsoTau40er_ETMHF100;   //!
   TBranch        *b_L1_IsoTau40er_ETMHF110;   //!
   TBranch        *b_L1_IsoTau40er_ETMHF120;   //!
   TBranch        *b_L1_IsoTau40er_ETMHF80;   //!
   TBranch        *b_L1_IsoTau40er_ETMHF90;   //!
   TBranch        *b_L1_IsolatedBunch;   //!
   TBranch        *b_L1_LastCollisionInTrain;   //!
   TBranch        *b_L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG24er2p1_HTT100er;   //!
   TBranch        *b_L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG24er2p1_Jet26er2p7_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG24er2p1_TripleJet_26er2p7_26_26er2p7;   //!
   TBranch        *b_L1_LooseIsoEG26er2p1_HTT100er;   //!
   TBranch        *b_L1_LooseIsoEG26er2p1_Jet34er2p7_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG28er2p1_HTT100er;   //!
   TBranch        *b_L1_LooseIsoEG28er2p1_Jet34er2p7_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG30er2p1_Jet34er2p7_dR_Min0p3;   //!
   TBranch        *b_L1_MU20_EG15;   //!
   TBranch        *b_L1_MinimumBiasHF0_AND_BptxAND;   //!
   TBranch        *b_L1_MinimumBiasHF0_OR_BptxAND;   //!
   TBranch        *b_L1_Mu10er2p1_ETM30;   //!
   TBranch        *b_L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6;   //!
   TBranch        *b_L1_Mu12_EG10;   //!
   TBranch        *b_L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6;   //!
   TBranch        *b_L1_Mu14er2p1_ETM30;   //!
   TBranch        *b_L1_Mu15_HTT100er;   //!
   TBranch        *b_L1_Mu18_HTT100er;   //!
   TBranch        *b_L1_Mu18_Jet24er2p7;   //!
   TBranch        *b_L1_Mu18er2p1_IsoTau26er2p1;   //!
   TBranch        *b_L1_Mu18er2p1_Tau24er2p1;   //!
   TBranch        *b_L1_Mu20_EG10;   //!
   TBranch        *b_L1_Mu20_EG17;   //!
   TBranch        *b_L1_Mu20_LooseIsoEG6;   //!
   TBranch        *b_L1_Mu20er2p1_IsoTau26er2p1;   //!
   TBranch        *b_L1_Mu20er2p1_IsoTau27er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau28er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau30er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau32er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau33er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau34er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau35er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau36er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau38er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau40er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_Tau50er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_Tau70er2p1;   //!
   TBranch        *b_L1_Mu23_EG10;   //!
   TBranch        *b_L1_Mu23_LooseIsoEG10;   //!
   TBranch        *b_L1_Mu3_Jet120er2p7_dEta_Max0p4_dPhi_Max0p4;   //!
   TBranch        *b_L1_Mu3_Jet16er2p7_dEta_Max0p4_dPhi_Max0p4;   //!
   TBranch        *b_L1_Mu3_Jet30er2p5;   //!
   TBranch        *b_L1_Mu3_Jet60er2p7_dEta_Max0p4_dPhi_Max0p4;   //!
   TBranch        *b_L1_Mu5_EG15;   //!
   TBranch        *b_L1_Mu5_EG20;   //!
   TBranch        *b_L1_Mu5_EG23;   //!
   TBranch        *b_L1_Mu5_LooseIsoEG18;   //!
   TBranch        *b_L1_Mu5_LooseIsoEG20;   //!
   TBranch        *b_L1_Mu6_DoubleEG10;   //!
   TBranch        *b_L1_Mu6_DoubleEG17;   //!
   TBranch        *b_L1_Mu6_HTT200er;   //!
   TBranch        *b_L1_Mu6_HTT240er;   //!
   TBranch        *b_L1_Mu6_HTT250er;   //!
   TBranch        *b_L1_Mu7_EG23;   //!
   TBranch        *b_L1_Mu7_LooseIsoEG20;   //!
   TBranch        *b_L1_Mu7_LooseIsoEG23;   //!
   TBranch        *b_L1_Mu8_HTT150er;   //!
   TBranch        *b_L1_NotBptxOR;   //!
   TBranch        *b_L1_QuadJet36er2p7_IsoTau52er2p1;   //!
   TBranch        *b_L1_QuadJet36er2p7_Tau52;   //!
   TBranch        *b_L1_QuadJet40er2p7;   //!
   TBranch        *b_L1_QuadJet50er2p7;   //!
   TBranch        *b_L1_QuadJet60er2p7;   //!
   TBranch        *b_L1_QuadMu0;   //!
   TBranch        *b_L1_SingleEG10;   //!
   TBranch        *b_L1_SingleEG15;   //!
   TBranch        *b_L1_SingleEG18;   //!
   TBranch        *b_L1_SingleEG24;   //!
   TBranch        *b_L1_SingleEG26;   //!
   TBranch        *b_L1_SingleEG28;   //!
   TBranch        *b_L1_SingleEG2_BptxAND;   //!
   TBranch        *b_L1_SingleEG30;   //!
   TBranch        *b_L1_SingleEG32;   //!
   TBranch        *b_L1_SingleEG34;   //!
   TBranch        *b_L1_SingleEG34er2p1;   //!
   TBranch        *b_L1_SingleEG36;   //!
   TBranch        *b_L1_SingleEG36er2p1;   //!
   TBranch        *b_L1_SingleEG38;   //!
   TBranch        *b_L1_SingleEG38er2p1;   //!
   TBranch        *b_L1_SingleEG40;   //!
   TBranch        *b_L1_SingleEG42;   //!
   TBranch        *b_L1_SingleEG45;   //!
   TBranch        *b_L1_SingleEG5;   //!
   TBranch        *b_L1_SingleEG50;   //!
   TBranch        *b_L1_SingleIsoEG18;   //!
   TBranch        *b_L1_SingleIsoEG18er2p1;   //!
   TBranch        *b_L1_SingleIsoEG20;   //!
   TBranch        *b_L1_SingleIsoEG20er2p1;   //!
   TBranch        *b_L1_SingleIsoEG22;   //!
   TBranch        *b_L1_SingleIsoEG22er2p1;   //!
   TBranch        *b_L1_SingleIsoEG24;   //!
   TBranch        *b_L1_SingleIsoEG24er2p1;   //!
   TBranch        *b_L1_SingleIsoEG26;   //!
   TBranch        *b_L1_SingleIsoEG26er2p1;   //!
   TBranch        *b_L1_SingleIsoEG28;   //!
   TBranch        *b_L1_SingleIsoEG28er2p1;   //!
   TBranch        *b_L1_SingleIsoEG30;   //!
   TBranch        *b_L1_SingleIsoEG30er2p1;   //!
   TBranch        *b_L1_SingleIsoEG32;   //!
   TBranch        *b_L1_SingleIsoEG32er2p1;   //!
   TBranch        *b_L1_SingleIsoEG33er2p1;   //!
   TBranch        *b_L1_SingleIsoEG34;   //!
   TBranch        *b_L1_SingleIsoEG34er2p1;   //!
   TBranch        *b_L1_SingleIsoEG35;   //!
   TBranch        *b_L1_SingleIsoEG35er2p1;   //!
   TBranch        *b_L1_SingleIsoEG36;   //!
   TBranch        *b_L1_SingleIsoEG36er2p1;   //!
   TBranch        *b_L1_SingleIsoEG37;   //!
   TBranch        *b_L1_SingleIsoEG38;   //!
   TBranch        *b_L1_SingleIsoEG38er2p1;   //!
   TBranch        *b_L1_SingleIsoEG40;   //!
   TBranch        *b_L1_SingleIsoEG40er2p1;   //!
   TBranch        *b_L1_SingleJet120;   //!
   TBranch        *b_L1_SingleJet120_FWD;   //!
   TBranch        *b_L1_SingleJet12_BptxAND;   //!
   TBranch        *b_L1_SingleJet140;   //!
   TBranch        *b_L1_SingleJet150;   //!
   TBranch        *b_L1_SingleJet16;   //!
   TBranch        *b_L1_SingleJet160;   //!
   TBranch        *b_L1_SingleJet170;   //!
   TBranch        *b_L1_SingleJet180;   //!
   TBranch        *b_L1_SingleJet20;   //!
   TBranch        *b_L1_SingleJet200;   //!
   TBranch        *b_L1_SingleJet20er2p7_NotBptxOR;   //!
   TBranch        *b_L1_SingleJet20er2p7_NotBptxOR_3BX;   //!
   TBranch        *b_L1_SingleJet35;   //!
   TBranch        *b_L1_SingleJet35_FWD;   //!
   TBranch        *b_L1_SingleJet35_HFm;   //!
   TBranch        *b_L1_SingleJet35_HFp;   //!
   TBranch        *b_L1_SingleJet43er2p7_NotBptxOR_3BX;   //!
   TBranch        *b_L1_SingleJet46er2p7_NotBptxOR_3BX;   //!
   TBranch        *b_L1_SingleJet60;   //!
   TBranch        *b_L1_SingleJet60_FWD;   //!
   TBranch        *b_L1_SingleJet60_HFm;   //!
   TBranch        *b_L1_SingleJet60_HFp;   //!
   TBranch        *b_L1_SingleJet90;   //!
   TBranch        *b_L1_SingleJet90_FWD;   //!
   TBranch        *b_L1_SingleMu0_BMTF;   //!
   TBranch        *b_L1_SingleMu0_EMTF;   //!
   TBranch        *b_L1_SingleMu0_OMTF;   //!
   TBranch        *b_L1_SingleMu10_LowQ;   //!
   TBranch        *b_L1_SingleMu11_LowQ;   //!
   TBranch        *b_L1_SingleMu12_LowQ_BMTF;   //!
   TBranch        *b_L1_SingleMu12_LowQ_EMTF;   //!
   TBranch        *b_L1_SingleMu12_LowQ_OMTF;   //!
   TBranch        *b_L1_SingleMu14er2p1;   //!
   TBranch        *b_L1_SingleMu16;   //!
   TBranch        *b_L1_SingleMu16er2p1;   //!
   TBranch        *b_L1_SingleMu18;   //!
   TBranch        *b_L1_SingleMu18er2p1;   //!
   TBranch        *b_L1_SingleMu20;   //!
   TBranch        *b_L1_SingleMu20er2p1;   //!
   TBranch        *b_L1_SingleMu22;   //!
   TBranch        *b_L1_SingleMu22_BMTF;   //!
   TBranch        *b_L1_SingleMu22_EMTF;   //!
   TBranch        *b_L1_SingleMu22_OMTF;   //!
   TBranch        *b_L1_SingleMu22er2p1;   //!
   TBranch        *b_L1_SingleMu25;   //!
   TBranch        *b_L1_SingleMu3;   //!
   TBranch        *b_L1_SingleMu30;   //!
   TBranch        *b_L1_SingleMu5;   //!
   TBranch        *b_L1_SingleMu7;   //!
   TBranch        *b_L1_SingleMuCosmics;   //!
   TBranch        *b_L1_SingleMuCosmics_BMTF;   //!
   TBranch        *b_L1_SingleMuCosmics_EMTF;   //!
   TBranch        *b_L1_SingleMuCosmics_OMTF;   //!
   TBranch        *b_L1_SingleMuOpen;   //!
   TBranch        *b_L1_SingleMuOpen_NotBptxOR;   //!
   TBranch        *b_L1_SingleMuOpen_NotBptxOR_3BX;   //!
   TBranch        *b_L1_SingleTau100er2p1;   //!
   TBranch        *b_L1_SingleTau120er2p1;   //!
   TBranch        *b_L1_SingleTau130er2p1;   //!
   TBranch        *b_L1_SingleTau140er2p1;   //!
   TBranch        *b_L1_SingleTau20;   //!
   TBranch        *b_L1_SingleTau80er2p1;   //!
   TBranch        *b_L1_TripleEG_14_10_8;   //!
   TBranch        *b_L1_TripleEG_18_17_8;   //!
   TBranch        *b_L1_TripleEG_LooseIso20_10_5;   //!
   TBranch        *b_L1_TripleJet_100_85_72_VBF;   //!
   TBranch        *b_L1_TripleJet_105_85_76_VBF;   //!
   TBranch        *b_L1_TripleJet_84_68_48_VBF;   //!
   TBranch        *b_L1_TripleJet_88_72_56_VBF;   //!
   TBranch        *b_L1_TripleJet_92_76_64_VBF;   //!
   TBranch        *b_L1_TripleJet_98_83_71_VBF;   //!
   TBranch        *b_L1_TripleMu0;   //!
   TBranch        *b_L1_TripleMu0_OQ;   //!
   TBranch        *b_L1_TripleMu3;   //!
   TBranch        *b_L1_TripleMu3_SQ;   //!
   TBranch        *b_L1_TripleMu_4_4_4;   //!
   TBranch        *b_L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17;   //!
   TBranch        *b_L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_8to14;   //!
   TBranch        *b_L1_TripleMu_5SQ_3SQ_0OQ;   //!
   TBranch        *b_L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9;   //!
   TBranch        *b_L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9;   //!
   TBranch        *b_L1_TripleMu_5_0_0;   //!
   TBranch        *b_L1_TripleMu_5_3_3;   //!
   TBranch        *b_L1_TripleMu_5_3p5_2p5;   //!
   TBranch        *b_L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17;   //!
   TBranch        *b_L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17;   //!
   TBranch        *b_L1_TripleMu_5_5_3;   //!
   TBranch        *b_L1_UnpairedBunchBptxMinus;   //!
   TBranch        *b_L1_UnpairedBunchBptxPlus;   //!
   TBranch        *b_L1_ZeroBias;   //!
   TBranch        *b_L1_ZeroBias_copy;   //!
   TBranch        *b_Jet_pt_raw;   //!
   TBranch        *b_Jet_pt_nom;   //!
   TBranch        *b_Jet_mass_raw;   //!
   TBranch        *b_Jet_mass_nom;   //!
   TBranch        *b_Jet_corr_JEC;   //!
   TBranch        *b_Jet_corr_JER;   //!
   TBranch        *b_MET_pt_nom;   //!
   TBranch        *b_MET_phi_nom;   //!
   TBranch        *b_MET_pt_jer;   //!
   TBranch        *b_MET_phi_jer;   //!
   TBranch        *b_Jet_pt_jerUp;   //!
   TBranch        *b_Jet_mass_jerUp;   //!
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
   TBranch        *b_FatJet_msoftdrop_tau21DDT_jerUp;   //!
   TBranch        *b_FatJet_msoftdrop_tau21DDT_jmrUp;   //!
   TBranch        *b_FatJet_msoftdrop_tau21DDT_jmsUp;   //!
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
   TBranch        *b_FatJet_msoftdrop_tau21DDT_jerDown;   //!
   TBranch        *b_FatJet_msoftdrop_tau21DDT_jmrDown;   //!
   TBranch        *b_FatJet_msoftdrop_tau21DDT_jmsDown;   //!
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
    
   TTree *Tout ;
   TTree *Tout1 ;
   
   int ncpdf;
   int ncscale;

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
   
   TH1D *hist_pt_GEN_ak8;
   TH1D *hist_pt_GEN_ak8_lead;
   TH1D *hist_eta_GEN_ak8;
   TH1D *hist_eta_GEN_ak8_lead;
   
   TH1D *hist_pt_GEN_ak4;
   TH1D *hist_pt_GEN_ak4_lead;
   TH1D *hist_eta_GEN_ak4;
   TH1D *hist_eta_GEN_ak4_lead;
   
   TH1D *hist_mtb_GEN;

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
   TH1D *hist_nmuons_trig_pass;
   
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
   TH2D *hist_2D_pt_btag_AK4_md_DeepAK8_beta2[nsbtag][ntoptag][netabins];
   TH2D *hist_2D_mtb_btag_AK4_md_DeepAK8_beta[nsbtag][ntoptag][netabins];
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
   TH1D *hist_biso_TopAK8score_masscut_1;
   TH1D *hist_biso_TopAK8score_MD_masscut_1;
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
   
   TH1D *hist_tbmass_md_DeepAK8_fine[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_fine_puup[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_fine_pudn[ntoptag_DAK8][nbtag][nsbtag];
   TH2D *h2d_mtb_npu[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_eta1[ntoptag_DAK8][nbtag][nsbtag][netabins];
   TH1D *hist_tbmass_md_DeepAK8_eta2[ntoptag_DAK8][nbtag][nsbtag][netabins];
   
   TH1D *hist_tbmass_md_DeepAK8_dY[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_trigwt[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_pfup[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_pfdn[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_2[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_pfup_2[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_pfdn_2[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_dAK8up[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_dAK8dn[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_puup[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_pudn[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_prefireup[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_prefiredn[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_bcorup[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_bcordn[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_bcorup_2[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_bcordn_2[ntoptag_DAK8][nbtag][nsbtag];
   
   TH1D *hist_tbmass_md_DeepAK8_noptw[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_alphaup[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_alphadn[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_betaup[ntoptag_DAK8][nbtag][nsbtag];
   TH1D *hist_tbmass_md_DeepAK8_betadn[ntoptag_DAK8][nbtag][nsbtag];
   
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
   
   TH1D *hist_tbmass_md_DeepAK8_JES_up_2[ntoptag_DAK8][nbtag][nsbtag][njesmax];
   TH1D *hist_tbmass_md_DeepAK8_JES_dn_2[ntoptag_DAK8][nbtag][nsbtag][njesmax];
   
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
   
   float btagvalue ;
   float btagvalue_deepCSV;
   float btagvalue_deepFlavB;
   
   float dR_cut = 1.2;
   float dphi_cut = M_PI/2.;
   
   float btagvalue_deepCSV_m = 0.4941;
   float btagvalue_deepCSV_t = 0.8001;

   float tau32_cut_1 = 0.65;
   float tau32_cut_2 = 0.4;
   
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
    
   double sdmassbins[nsdmassbins+1] = {0,10,21,33,45,57,69,81,93,105,117,129,141,156,172,190,210,250,290,340,400};
   
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
    {0.987142,1.06626,0.928712,1.00808,0.98586,1.13492,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{1.10974,0.991804,0.783162,0.767303,0.802142,0.838826,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
   };  //2016 QCD
   
   double cor_b_sdmass_data_16[2][nsdmassbins] = {
    {1.00931,1.07519,0.949899,0.883415,0.917168,0.908525,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{1.07408,1.02367,0.843254,0.784695,0.758625,0.785022,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
   };  //2016 Data
    
 
  double cor_b_tight_sdmass_mc_16[2][nsdmassbins] = {
   {0.944359,0.986263,0.968458,1.14505,1.20321,1.44714,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
   {1.01366,0.956484,0.917361,0.984545,1.10504,1.19774,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    };  //2016 QCD b tight
    
  double cor_b_tight_sdmass_data_16[2][nsdmassbins] = {
    {0.972227,1.02831,1.00535,0.995104,1.11508,1.0687,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{1.01784,0.998772,0.91469,0.980248,0.990516,1.04135,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
   };  //2016 Data b tight
      
     
   double cor_b_sdmass_mc_17[2][nsdmassbins] = {
    {0.894106,1.00168,1.08224,1.23211,1.42184,1.44549,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0.948099,0.965375,1.01162,1.14372,1.32161,1.33012,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
   };  //2017 QCD
   
   double cor_b_sdmass_data_17[2][nsdmassbins] = {
    {0.929639,0.981946,1.03214,1.18439,1.32073,1.29275,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0.96212,0.915459,1.01051,1.18366,1.35503,1.43479,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
   };  //2017 Data
 
   double cor_b_sdmass_mc[2][nsdmassbins] = {{0}};
   double cor_b_sdmass_data[2][nsdmassbins] = {{0}};
   
   float sfwt_tau32;
   float sfwt_deepak8;
   float sfwt_deepak8_md, sfwt_deepak8_md_up, sfwt_deepak8_md_dn; 
   
   static const int nohtbins = 30;
     
double htbins[nohtbins+1] = {300, 330, 362, 395, 430, 468,
     507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
     1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1784, 2000,
     2116, 2500, 2941, 3637, 5000} ; 

float trig_weight_factors[ntopsdmassbins][nohtbins];

float trig_weight_factors_2016_QCD[ntopsdmassbins][nohtbins] ={
 {1,1,1,1,1.10345,0.805232,0.182536,0.197765,0.380912,0.670897,0.834252,0.918474,0.982272,1.0125,1.00699,1.00193,1.00054,1.00024,1.00011,1,1.00011,1,1,1,1,1.0011,1,1,1,1},
 {1,1,1,1,1,1,1.12,0.808403,0.643965,0.847031,0.985109,0.995466,0.999472,1.0007,1.00054,1.00009,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
 {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}
 };
float trig_weight_factors_2016_TT[ntopsdmassbins][nohtbins] ={
 {1,1,1,1,1.0704,1.09308,0.759866,0.764611,0.957048,1.09061,1.10632,1.07927,1.0546,1.02967,1.01028,1.00242,1.00029,0.999902,0.999861,0.99987,1.00011,1,0.999135,1,1,1.0011,1,1,1,1},
 {1,1,1,1,1,1,1.04278,0.972386,0.965671,1.03579,1.02911,1.00574,1.00067,1.00018,1.0003,1.00005,0.999761,0.999835,0.999962,0.99983,0.999896,1,1,1,1,1,1,1,1,1},
 {1,1,1,1,1,1,1,1,1,1,1,0.994339,1,1,0.999469,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}
 };
float trig_weight_factors_2016_ST[ntopsdmassbins][nohtbins] ={
 {1,1,1,1,1.08177,1.05229,0.621771,0.588201,0.75973,0.956921,1.03433,1.05128,1.0467,1.0289,1.00468,1.00274,1.00054,1.00024,1.00011,1,1.00011,1,1,1,1,1.0011,1,1,1,1},
 {1,1,1,1,1,1,1.12,1.00739,0.855615,0.995367,1.0131,1.00596,1.00192,1.0007,1.00054,1.00009,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
 {1,1,1,1,1,1,1,1,1,1,1,0.849112,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}
 };

float trig_weight_factors_2017_QCD[ntopsdmassbins][nohtbins]={
{1.52547,1,1.52041,1.27411,0.954823,1.63933,1.49198,1.78913,1.63537,1.67684,1.75807,1.71015,1.65438,1.49059,1.25432,1.09214,1.02392,1.00513,1.00158,1.00038,1.00021,1.00052,1.00025,1.00011,1,1,0.992649,1.00118,1,1.03017},
{1,1,1,1,1,644855,3.99427,1.15384,1.26838,2.61535,3.64925,3.31482,2.35012,1.6682,1.28066,1.07543,1.01339,1.0023,1.00063,1.00038,1.00067,1.00151,1,1.00037,1,1,1,1,1,1},
{1,1,1,1,1,1,1,1,1,1,5.68444,1.01553,1.33748,1.29369,1.11957,1.03858,1.01762,1.00597,1.00064,1,1,1,1,1,1,1,1,1,1,1}
 };
float trig_weight_factors_2017_TT[ntopsdmassbins][nohtbins]={
{1.52547,0.869586,1.52041,1.20361,2.66645,4.91677,16.7263,17.1385,11.7532,7.81608,5.92001,4.12231,2.84473,2.02241,1.48749,1.18026,1.04225,1.00709,0.999004,1.00094,1.00056,0.999149,1.00025,0.997346,1,1,1.00034,1.00118,1,1.03017},
{1,1,1,1,1,487617,3.92897,4.19968,4.32724,3.99868,3.95255,3.17049,2.30299,1.67155,1.28043,1.07588,1.01315,1.00246,1.00053,1.00038,1.00067,1.00151,1,1.00037,1,1,1,1,1,1},
{1,1,1,1,1,1,1,1,1,1,4.26425,1.37182,1.32734,1.27765,1.1428,1.0539,1.02123,1.00597,1.00064,1,1,1,1,1,1,1,1,1,1,1}
 };
float trig_weight_factors_2017_ST[ntopsdmassbins][nohtbins]{
{1.52547,1,1.52041,1.27411,3.03556,3.83096,9.5946,11.198,8.27829,5.7801,4.76236,3.71654,2.67955,1.94486,1.46723,1.1766,1.04315,1.00835,1.00216,1.00094,1.00056,1.00075,1.00025,1.00011,1,1,1.00034,1.00118,1,1.03017},
{1,1,1,1,1,197.61,8.34002,1.73837,4.15467,2.85729,2.96534,2.54389,2.1116,1.63588,1.27316,1.07704,1.01339,1.00257,1.00063,1.00038,1.00067,1.00151,1,1.00037,1,1,1,1,1,1},
{1,1,1,1,1,1,1,1,1,1,5.68444,1.54298,1.42353,0.933951,1.14851,1.05392,1.02123,1.00597,1.00064,1,1,1,1,1,1,1,1,1,1,1}
 };

float trig_weight_factors_2018_QCD[ntopsdmassbins][nohtbins]{
{1.52547,1,1.52041,1.27411,0.954823,1.63933,1.49198,1.78913,1.63537,1.67684,1.75807,1.71015,1.65438,1.49059,1.25432,1.09214,1.02392,1.00513,1.00158,1.00038,1.00021,1.00052,1.00025,1.00011,1,1,0.992649,1.00118,1,1.03017},
{1,1,1,1,1,644855,3.99427,1.15384,1.26838,2.61535,3.64925,3.31482,2.35012,1.6682,1.28066,1.07543,1.01339,1.0023,1.00063,1.00038,1.00067,1.00151,1,1.00037,1,1,1,1,1,1},
{1,1,1,1,1,1,1,1,1,1,5.68444,1.01553,1.33748,1.29369,1.11957,1.03858,1.01762,1.00597,1.00064,1,1,1,1,1,1,1,1,1,1,1}
 };
float trig_weight_factors_2018_TT[ntopsdmassbins][nohtbins]{
{1.52547,0.869586,1.52041,1.20361,2.66645,4.91677,16.7263,17.1385,11.7532,7.81608,5.92001,4.12231,2.84473,2.02241,1.48749,1.18026,1.04225,1.00709,0.999004,1.00094,1.00056,0.999149,1.00025,0.997346,1,1,1.00034,1.00118,1,1.03017},
{1,1,1,1,1,487617,3.92897,4.19968,4.32724,3.99868,3.95255,3.17049,2.30299,1.67155,1.28043,1.07588,1.01315,1.00246,1.00053,1.00038,1.00067,1.00151,1,1.00037,1,1,1,1,1,1},
{1,1,1,1,1,1,1,1,1,1,4.26425,1.37182,1.32734,1.27765,1.1428,1.0539,1.02123,1.00597,1.00064,1,1,1,1,1,1,1,1,1,1,1}
 };
float trig_weight_factors_2018_ST[ntopsdmassbins][nohtbins]{
{1.52547,1,1.52041,1.27411,3.03556,3.83096,9.5946,11.198,8.27829,5.7801,4.76236,3.71654,2.67955,1.94486,1.46723,1.1766,1.04315,1.00835,1.00216,1.00094,1.00056,1.00075,1.00025,1.00011,1,1,1.00034,1.00118,1,1.03017},
{1,1,1,1,1,197.61,8.34002,1.73837,4.15467,2.85729,2.96534,2.54389,2.1116,1.63588,1.27316,1.07704,1.01339,1.00257,1.00063,1.00038,1.00067,1.00151,1,1.00037,1,1,1,1,1,1},
{1,1,1,1,1,1,1,1,1,1,5.68444,1.54298,1.42353,0.933951,1.14851,1.05392,1.02123,1.00597,1.00064,1,1,1,1,1,1,1,1,1,1,1}
 };

