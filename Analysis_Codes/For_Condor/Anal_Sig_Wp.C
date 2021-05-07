#include "Anal_Sig_Wp.h"
#include <TStyle.h>
#include <TCanvas.h>

#define Btagger_DeepJet
//#define Btagger_DeepCSV
//#define Btagger_CSVv2
#define Btagger_DeepJet_wt
//#define Btagger_DeepCSV_wt

//#define Anal_2017
//#define Data_2017B
//#define Anal_2016
//#define Data_2016H
#define Anal_2018
//#define HEM_Cor_ON

//#define B_Tight

//#define LHAPDFX

using namespace std;

void fillarray(bool isMC)
{
 unsigned int fmu = 0;
   
 for(unsigned int imu=0; imu < nMuon; imu++){
	   
	if(!(Muon_tightId[imu])) continue;
	if(!(Muon_pfIsoId[imu]>=2)) continue;
//	if(Muon_pfRelIso04_all[imu] >= 0.15) continue;
	
	if(fabs(Muon_pt[imu]) < lepptcut) continue;
	if(fabs(Muon_eta[imu]) > 2.4) continue;
	if(fabs(Muon_dxy[imu])>0.045 || fabs(Muon_dz[imu])>0.2) continue;
	   
	Muon_pt[fmu] = Muon_pt[imu];
	Muon_eta[fmu] = Muon_eta[imu];
	Muon_phi[fmu] = Muon_phi[imu];
	Muon_mass[fmu] = Muon_mass[imu];
	Muon_charge[fmu] = Muon_charge[imu];
	Muon_jetIdx[fmu] = Muon_jetIdx[imu];
	Muon_tightId[fmu] = Muon_tightId[imu];
	Muon_pfIsoId[fmu] = Muon_pfIsoId[imu];
	Muon_jetPtRelv2[fmu] = Muon_jetPtRelv2[imu];
	   
	fmu++; 
	if(fmu>=njetmax) break;	
   }
	   
   nMuon = fmu;
   
 unsigned int fel = 0;
   
 for(unsigned iel=0; iel < nElectron; iel++){
	   
//  if(Electron_cutBased[iel]!=4) continue;
	if(Electron_mvaFall17V2Iso_WP90[iel]==0) continue;
	if(!(Electron_convVeto[iel]) || (Electron_lostHits[iel]>1)) continue;
	
	if(fabs(Electron_pt[iel]) < lepptcut || fabs(Electron_eta[iel])>2.4) continue;
	if(fabs(Electron_dxy[iel])>0.045 || fabs(Electron_dz[iel])>0.2) continue;
	   
	Electron_pt[fel] = Electron_pt[iel];
	Electron_eta[fel] = Electron_eta[iel];
	Electron_phi[fel] = Electron_phi[iel];
	Electron_mass[fel] = Electron_mass[iel];
	Electron_jetIdx[fel] = Electron_jetIdx[iel];
	Electron_photonIdx[fel] = Electron_photonIdx[iel];
	Electron_cutBased_HEEP[fel] = Electron_cutBased_HEEP[iel];
	Electron_cutBased[fel] = Electron_cutBased[iel];
	   
	fel++;
	if(fel>=njetmax) break;
   }
	   
	nElectron = fel;
	
 unsigned int fpho = 0;
   
 for(unsigned ipho=0; ipho < nPhoton; ipho++){
	   
	if(!(Photon_mvaID_WP80[ipho])) continue;
	if(fabs(Photon_pt[ipho]) < lepptcut || fabs(Photon_eta[ipho])>jeteta_cut) continue;
	if(!(Photon_electronVeto[ipho]) || (Photon_pixelSeed[ipho])) continue;
	
	Photon_pt[fpho] = Photon_pt[ipho];
	Photon_eta[fpho] = Photon_eta[ipho];
	Photon_phi[fpho] = Photon_phi[ipho];
	Photon_mass[fpho] = Photon_mass[ipho];
	
	fpho++;
	if(fpho>=njetmax) break;
   }
	   
   nPhoton = fpho;

   int fjet = 0;
   
   for(unsigned ijet=0; ijet<nFatJet; ijet++){
	  
	   if(ijet==0 &&  (FatJet_pt[ijet]<200. || abs(FatJet_eta[ijet])>4.0)) break;
 
	   if(FatJet_jetId[ijet]==0) continue;
	   int jetid = *(decToBinary(FatJet_jetId[ijet])+1);//((decToBinary(FatJet_jetId[ijet]))/10)%10;  // tight id
	   if(jetid!=1) continue;
	   
	   if(FatJet_pt[ijet]<30.) continue;
	   if(abs(FatJet_eta[ijet])>jeteta_cut) continue;
	   
	   bool mutag = false;
	   for(unsigned imu=0; imu<nMuon; imu++){
		   if(delta2R(Muon_eta[imu],Muon_phi[imu],FatJet_eta[ijet],FatJet_phi[ijet])<0.8){
			   mutag = true;
			   break;
			   }
		   }
	   
	   bool etag = false;
	   for(unsigned iel=0; iel<nElectron; iel++){
		   if(delta2R(Electron_eta[iel],Electron_phi[iel],FatJet_eta[ijet],FatJet_phi[ijet])<0.8){
			   etag = true;
			   break;
			   }
		   }
	   
	   if(mutag||etag) continue;

	   FatJet_area[fjet] = FatJet_area[ijet];
	   FatJet_btagCMVA[fjet] = FatJet_btagCMVA[ijet];
	   FatJet_btagCSVV2[fjet] = FatJet_btagCSVV2[ijet];
	   FatJet_btagDeepB[fjet] = FatJet_btagDeepB[ijet];
	   FatJet_btagHbb[fjet] = FatJet_btagHbb[ijet];
	   FatJet_deepTagMD_H4qvsQCD[fjet] = FatJet_deepTagMD_H4qvsQCD[ijet];
	   FatJet_deepTagMD_HbbvsQCD[fjet] = FatJet_deepTagMD_HbbvsQCD[ijet];
	   FatJet_deepTagMD_TvsQCD[fjet] = FatJet_deepTagMD_TvsQCD[ijet];
	   FatJet_deepTagMD_WvsQCD[fjet] = FatJet_deepTagMD_WvsQCD[ijet];
	   FatJet_deepTagMD_ZHbbvsQCD[fjet] = FatJet_deepTagMD_ZHbbvsQCD[ijet];
	   FatJet_deepTagMD_ZHccvsQCD[fjet] = FatJet_deepTagMD_ZHccvsQCD[ijet];
	   FatJet_deepTagMD_ZbbvsQCD[fjet] = FatJet_deepTagMD_ZbbvsQCD[ijet];
	   FatJet_deepTagMD_ZvsQCD[fjet] = FatJet_deepTagMD_ZvsQCD[ijet];
	   FatJet_deepTagMD_bbvsLight[fjet] = FatJet_deepTagMD_bbvsLight[ijet];
	   FatJet_deepTagMD_ccvsLight[fjet] = FatJet_deepTagMD_ccvsLight[ijet];
	   FatJet_deepTag_TvsQCD[fjet] = FatJet_deepTag_TvsQCD[ijet];
	   FatJet_deepTag_WvsQCD[fjet] = FatJet_deepTag_WvsQCD[ijet];
	   FatJet_deepTag_ZvsQCD[fjet] = FatJet_deepTag_ZvsQCD[ijet];
	   FatJet_eta[fjet] = FatJet_eta[ijet];
	   #ifdef NoJEC
	   FatJet_mass[fjet] = FatJet_mass[ijet];
	   FatJet_msoftdrop[fjet] = FatJet_msoftdrop[ijet];
	   FatJet_pt[fjet] = FatJet_pt[ijet];
	   #else
	   if(isMC){
	   FatJet_mass[fjet] = (FatJet_mass_raw[ijet]*FatJet_corr_JEC[ijet]*FatJet_corr_JER[ijet]);	   
	   FatJet_msoftdrop[fjet] = FatJet_msoftdrop_raw[ijet]*FatJet_corr_JER[ijet];//*FatJet_msoftdrop_corr_JMS[ijet]*FatJet_msoftdrop_corr_JMR[ijet];
	   }else{
		   FatJet_mass[fjet] = (FatJet_mass_raw[ijet]*FatJet_corr_JEC[ijet]);		
		   FatJet_msoftdrop[fjet] = FatJet_msoftdrop_raw[ijet];//*FatJet_msoftdrop_corr_JMS[ijet];
		   }
	   FatJet_pt[fjet] = FatJet_pt_nom[ijet];
	   #endif
	   
	   FatJet_n2b1[fjet] = FatJet_n2b1[ijet];
	   FatJet_n3b1[fjet] = FatJet_n3b1[ijet];
	   FatJet_phi[fjet] = FatJet_phi[ijet];
	   FatJet_rawFactor[fjet] = FatJet_rawFactor[ijet];
	   FatJet_tau1[fjet] = FatJet_tau1[ijet];
	   FatJet_tau2[fjet] = FatJet_tau2[ijet];
	   FatJet_tau3[fjet] = FatJet_tau3[ijet];
	   FatJet_tau4[fjet] = FatJet_tau4[ijet];
	   FatJet_jetId[fjet] = FatJet_jetId[ijet];
	   FatJet_subJetIdx1[fjet] = FatJet_subJetIdx1[ijet];
	   FatJet_subJetIdx2[fjet] = FatJet_subJetIdx2[ijet];
	  
	   if(isMC){
		   
		   float mindR = 0.4;
		   int matchgen = -1;
		   
		   for(unsigned igen=0; igen<nGenJetAK8; igen++){
				float delR = delta2R(FatJet_eta[fjet],FatJet_phi[fjet],GenJetAK8_eta[igen],GenJetAK8_phi[igen]);
				if(delR < mindR) { mindR = delR; matchgen = igen; }
		   }
	   
		   if(matchgen>=0){
				FatJet_MatchGenJet[fjet] = matchgen;
				FatJet_GenJetpt[fjet] = GenJetAK8_pt[matchgen];
				FatJet_GenJeteta[fjet] = GenJetAK8_eta[matchgen];
				FatJet_GenJetphi[fjet] = GenJetAK8_phi[matchgen];
				FatJet_GenJetmass[fjet] = GenJetAK8_mass[matchgen];
			}else{
			    FatJet_MatchGenJet[fjet] = FatJet_GenJetpt[fjet] = FatJet_GenJeteta[fjet] = FatJet_GenJetphi[fjet] = FatJet_GenJetmass[fjet] = -100;
			    }
	   }
	   
	     
	   #if defined(Anal_2018) && defined(HEM_Cor_ON)
	   if(isMC){
		   if(FatJet_phi[fjet]>-1.57 && FatJet_phi[fjet]<-0.87){
			   if(FatJet_eta[fjet]>-2.5 && FatJet_eta[fjet]<-1.3) { 
		  //	   scale down jES by 20 % for jets with -1.57 <phi< -0.87 and -2.5<eta<-1.3
				   FatJet_mass[fjet] = 0.8*FatJet_mass[fjet];
				   FatJet_msoftdrop[fjet] = 0.8*FatJet_msoftdrop[fjet];
				   FatJet_pt[fjet] = 0.8*FatJet_pt[fjet];
				   }
			   if(FatJet_eta[fjet]>-3.0 && FatJet_eta[fjet]<-2.5) {
		//	   	   scale down jES by 35 % for jets with -1.57 <phi< -0.87 and -3.0<eta<-2.5	   
				   FatJet_mass[fjet] = 0.65*FatJet_mass[fjet];
				   FatJet_msoftdrop[fjet] = 0.65*FatJet_msoftdrop[fjet];
				   FatJet_pt[fjet] = 0.65*FatJet_pt[fjet];
				   }	   
			   }
		   }
	   #endif
	   
	   // systematics //
	   
	   FatJet_pt_jer_up[fjet] = FatJet_pt_jerUp[ijet];
	   FatJet_mass_jer_up[fjet] = FatJet_mass_jerUp[ijet];
	   FatJet_msoftdrop_jer_up[fjet] = FatJet_msoftdrop_jerUp[ijet];
	   
	   FatJet_pt_jer_dn[fjet] = FatJet_pt_jerDown[ijet];
	   FatJet_mass_jer_dn[fjet] = FatJet_mass_jerDown[ijet];
	   FatJet_msoftdrop_jer_dn[fjet] = FatJet_msoftdrop_jerDown[ijet];
	   
	   FatJet_mass_jms_up[fjet] = FatJet_mass_jmsUp[ijet];
	   FatJet_msoftdrop_jms_up[fjet] = FatJet_msoftdrop_jmsUp[ijet];
	   FatJet_mass_jmr_up[fjet] = FatJet_mass_jmrUp[ijet];
	   FatJet_msoftdrop_jmr_up[fjet] = FatJet_msoftdrop_jmrUp[ijet];
	   
	   FatJet_mass_jms_dn[fjet] = FatJet_mass_jmsDown[ijet];
	   FatJet_msoftdrop_jms_dn[fjet] = FatJet_msoftdrop_jmsDown[ijet];
	   FatJet_mass_jmr_dn[fjet] = FatJet_mass_jmrDown[ijet];
	   FatJet_msoftdrop_jmr_dn[fjet] = FatJet_msoftdrop_jmrDown[ijet];
	   	
	   // jes //
	   
	   FatJet_pt_jes_up[0][fjet] = FatJet_pt_jesAbsoluteStatUp[ijet];
	   FatJet_pt_jes_up[1][fjet] = FatJet_pt_jesAbsoluteScaleUp[ijet];
	   FatJet_pt_jes_up[2][fjet] = FatJet_pt_jesAbsoluteFlavMapUp[ijet];
	   FatJet_pt_jes_up[3][fjet] = FatJet_pt_jesAbsoluteMPFBiasUp[ijet];
	   FatJet_pt_jes_up[4][fjet] = FatJet_pt_jesFragmentationUp[ijet];
	   FatJet_pt_jes_up[5][fjet] = FatJet_pt_jesSinglePionECALUp[ijet];
	   FatJet_pt_jes_up[6][fjet] = FatJet_pt_jesSinglePionHCALUp[ijet];
	   FatJet_pt_jes_up[7][fjet] = FatJet_pt_jesFlavorQCDUp[ijet];
	   FatJet_pt_jes_up[8][fjet] = FatJet_pt_jesTimePtEtaUp[ijet];
	   FatJet_pt_jes_up[9][fjet] = FatJet_pt_jesRelativeJEREC1Up[ijet];
	   FatJet_pt_jes_up[10][fjet] = FatJet_pt_jesRelativeJEREC2Up[ijet];
	   FatJet_pt_jes_up[11][fjet] = FatJet_pt_jesRelativeJERHFUp[ijet];
	   FatJet_pt_jes_up[12][fjet] = FatJet_pt_jesRelativePtBBUp[ijet];
	   FatJet_pt_jes_up[13][fjet] = FatJet_pt_jesRelativePtEC1Up[ijet];
	   FatJet_pt_jes_up[14][fjet] = FatJet_pt_jesRelativePtEC2Up[ijet];
	   FatJet_pt_jes_up[15][fjet] = FatJet_pt_jesRelativePtHFUp[ijet];
	   FatJet_pt_jes_up[16][fjet] = FatJet_pt_jesRelativeBalUp[ijet];
	   FatJet_pt_jes_up[17][fjet] = FatJet_pt_jesRelativeSampleUp[ijet];
	   FatJet_pt_jes_up[18][fjet] = FatJet_pt_jesRelativeFSRUp[ijet];
	   FatJet_pt_jes_up[19][fjet] = FatJet_pt_jesRelativeStatFSRUp[ijet];
	   FatJet_pt_jes_up[20][fjet] = FatJet_pt_jesRelativeStatECUp[ijet];
	   FatJet_pt_jes_up[21][fjet] = FatJet_pt_jesRelativeStatHFUp[ijet];
	   FatJet_pt_jes_up[22][fjet] = FatJet_pt_jesPileUpDataMCUp[ijet];
	   FatJet_pt_jes_up[23][fjet] = FatJet_pt_jesPileUpPtRefUp[ijet];
	   FatJet_pt_jes_up[24][fjet] = FatJet_pt_jesPileUpPtBBUp[ijet];
	   FatJet_pt_jes_up[25][fjet] = FatJet_pt_jesPileUpPtEC1Up[ijet];
	   FatJet_pt_jes_up[26][fjet] = FatJet_pt_jesPileUpPtEC2Up[ijet];
	   FatJet_pt_jes_up[27][fjet] = FatJet_pt_jesPileUpPtHFUp[ijet];
	   FatJet_pt_jes_up[28][fjet] = FatJet_pt_jesPileUpMuZeroUp[ijet];
	   FatJet_pt_jes_up[29][fjet] = FatJet_pt_jesPileUpEnvelopeUp[ijet];
	   FatJet_pt_jes_up[30][fjet] = FatJet_pt_jesSubTotalPileUpUp[ijet];
	   FatJet_pt_jes_up[31][fjet] = FatJet_pt_jesSubTotalRelativeUp[ijet];
	   FatJet_pt_jes_up[32][fjet] = FatJet_pt_jesSubTotalPtUp[ijet];
	   FatJet_pt_jes_up[33][fjet] = FatJet_pt_jesSubTotalScaleUp[ijet];
	   FatJet_pt_jes_up[34][fjet] = FatJet_pt_jesSubTotalAbsoluteUp[ijet];
	   FatJet_pt_jes_up[35][fjet] = FatJet_pt_jesSubTotalMCUp[ijet];
	   FatJet_pt_jes_up[36][fjet] = FatJet_pt_jesTotalUp[ijet];
	   
	   FatJet_mass_jes_up[0][fjet] = FatJet_mass_jesAbsoluteStatUp[ijet];
	   FatJet_mass_jes_up[1][fjet] = FatJet_mass_jesAbsoluteScaleUp[ijet];
	   FatJet_mass_jes_up[2][fjet] = FatJet_mass_jesAbsoluteFlavMapUp[ijet];
	   FatJet_mass_jes_up[3][fjet] = FatJet_mass_jesAbsoluteMPFBiasUp[ijet];
	   FatJet_mass_jes_up[4][fjet] = FatJet_mass_jesFragmentationUp[ijet];
	   FatJet_mass_jes_up[5][fjet] = FatJet_mass_jesSinglePionECALUp[ijet];
	   FatJet_mass_jes_up[6][fjet] = FatJet_mass_jesSinglePionHCALUp[ijet];
	   FatJet_mass_jes_up[7][fjet] = FatJet_mass_jesFlavorQCDUp[ijet];
	   FatJet_mass_jes_up[8][fjet] = FatJet_mass_jesTimePtEtaUp[ijet];
	   FatJet_mass_jes_up[9][fjet] = FatJet_mass_jesRelativeJEREC1Up[ijet];
	   FatJet_mass_jes_up[10][fjet] = FatJet_mass_jesRelativeJEREC2Up[ijet];
	   FatJet_mass_jes_up[11][fjet] = FatJet_mass_jesRelativeJERHFUp[ijet];
	   FatJet_mass_jes_up[12][fjet] = FatJet_mass_jesRelativePtBBUp[ijet];
	   FatJet_mass_jes_up[13][fjet] = FatJet_mass_jesRelativePtEC1Up[ijet];
	   FatJet_mass_jes_up[14][fjet] = FatJet_mass_jesRelativePtEC2Up[ijet];
	   FatJet_mass_jes_up[15][fjet] = FatJet_mass_jesRelativePtHFUp[ijet];
	   FatJet_mass_jes_up[16][fjet] = FatJet_mass_jesRelativeBalUp[ijet];
	   FatJet_mass_jes_up[17][fjet] = FatJet_mass_jesRelativeSampleUp[ijet];
	   FatJet_mass_jes_up[18][fjet] = FatJet_mass_jesRelativeFSRUp[ijet];
	   FatJet_mass_jes_up[19][fjet] = FatJet_mass_jesRelativeStatFSRUp[ijet];
	   FatJet_mass_jes_up[20][fjet] = FatJet_mass_jesRelativeStatECUp[ijet];
	   FatJet_mass_jes_up[21][fjet] = FatJet_mass_jesRelativeStatHFUp[ijet];
	   FatJet_mass_jes_up[22][fjet] = FatJet_mass_jesPileUpDataMCUp[ijet];
	   FatJet_mass_jes_up[23][fjet] = FatJet_mass_jesPileUpPtRefUp[ijet];
	   FatJet_mass_jes_up[24][fjet] = FatJet_mass_jesPileUpPtBBUp[ijet];
	   FatJet_mass_jes_up[25][fjet] = FatJet_mass_jesPileUpPtEC1Up[ijet];
	   FatJet_mass_jes_up[26][fjet] = FatJet_mass_jesPileUpPtEC2Up[ijet];
	   FatJet_mass_jes_up[27][fjet] = FatJet_mass_jesPileUpPtHFUp[ijet];
	   FatJet_mass_jes_up[28][fjet] = FatJet_mass_jesPileUpMuZeroUp[ijet];
	   FatJet_mass_jes_up[29][fjet] = FatJet_mass_jesPileUpEnvelopeUp[ijet];
	   FatJet_mass_jes_up[30][fjet] = FatJet_mass_jesSubTotalPileUpUp[ijet];
	   FatJet_mass_jes_up[31][fjet] = FatJet_mass_jesSubTotalRelativeUp[ijet];
	   FatJet_mass_jes_up[32][fjet] = FatJet_mass_jesSubTotalPtUp[ijet];
	   FatJet_mass_jes_up[33][fjet] = FatJet_mass_jesSubTotalScaleUp[ijet];
	   FatJet_mass_jes_up[34][fjet] = FatJet_mass_jesSubTotalAbsoluteUp[ijet];
	   FatJet_mass_jes_up[35][fjet] = FatJet_mass_jesSubTotalMCUp[ijet];
	   FatJet_mass_jes_up[36][fjet] = FatJet_mass_jesTotalUp[ijet];
	   
	   FatJet_msoftdrop_jes_up[0][fjet] = FatJet_msoftdrop_jesAbsoluteStatUp[ijet];
	   FatJet_msoftdrop_jes_up[1][fjet] = FatJet_msoftdrop_jesAbsoluteScaleUp[ijet];
	   FatJet_msoftdrop_jes_up[2][fjet] = FatJet_msoftdrop_jesAbsoluteFlavMapUp[ijet];
	   FatJet_msoftdrop_jes_up[3][fjet] = FatJet_msoftdrop_jesAbsoluteMPFBiasUp[ijet];
	   FatJet_msoftdrop_jes_up[4][fjet] = FatJet_msoftdrop_jesFragmentationUp[ijet];
	   FatJet_msoftdrop_jes_up[5][fjet] = FatJet_msoftdrop_jesSinglePionECALUp[ijet];
	   FatJet_msoftdrop_jes_up[6][fjet] = FatJet_msoftdrop_jesSinglePionHCALUp[ijet];
	   FatJet_msoftdrop_jes_up[7][fjet] = FatJet_msoftdrop_jesFlavorQCDUp[ijet];
	   FatJet_msoftdrop_jes_up[8][fjet] = FatJet_msoftdrop_jesTimePtEtaUp[ijet];
	   FatJet_msoftdrop_jes_up[9][fjet] = FatJet_msoftdrop_jesRelativeJEREC1Up[ijet];
	   FatJet_msoftdrop_jes_up[10][fjet] = FatJet_msoftdrop_jesRelativeJEREC2Up[ijet];
	   FatJet_msoftdrop_jes_up[11][fjet] = FatJet_msoftdrop_jesRelativeJERHFUp[ijet];
	   FatJet_msoftdrop_jes_up[12][fjet] = FatJet_msoftdrop_jesRelativePtBBUp[ijet];
	   FatJet_msoftdrop_jes_up[13][fjet] = FatJet_msoftdrop_jesRelativePtEC1Up[ijet];
	   FatJet_msoftdrop_jes_up[14][fjet] = FatJet_msoftdrop_jesRelativePtEC2Up[ijet];
	   FatJet_msoftdrop_jes_up[15][fjet] = FatJet_msoftdrop_jesRelativePtHFUp[ijet];
	   FatJet_msoftdrop_jes_up[16][fjet] = FatJet_msoftdrop_jesRelativeBalUp[ijet];
	   FatJet_msoftdrop_jes_up[17][fjet] = FatJet_msoftdrop_jesRelativeSampleUp[ijet];
	   FatJet_msoftdrop_jes_up[18][fjet] = FatJet_msoftdrop_jesRelativeFSRUp[ijet];
	   FatJet_msoftdrop_jes_up[19][fjet] = FatJet_msoftdrop_jesRelativeStatFSRUp[ijet];
	   FatJet_msoftdrop_jes_up[20][fjet] = FatJet_msoftdrop_jesRelativeStatECUp[ijet];
	   FatJet_msoftdrop_jes_up[21][fjet] = FatJet_msoftdrop_jesRelativeStatHFUp[ijet];
	   FatJet_msoftdrop_jes_up[22][fjet] = FatJet_msoftdrop_jesPileUpDataMCUp[ijet];
	   FatJet_msoftdrop_jes_up[23][fjet] = FatJet_msoftdrop_jesPileUpPtRefUp[ijet];
	   FatJet_msoftdrop_jes_up[24][fjet] = FatJet_msoftdrop_jesPileUpPtBBUp[ijet];
	   FatJet_msoftdrop_jes_up[25][fjet] = FatJet_msoftdrop_jesPileUpPtEC1Up[ijet];
	   FatJet_msoftdrop_jes_up[26][fjet] = FatJet_msoftdrop_jesPileUpPtEC2Up[ijet];
	   FatJet_msoftdrop_jes_up[27][fjet] = FatJet_msoftdrop_jesPileUpPtHFUp[ijet];
	   FatJet_msoftdrop_jes_up[28][fjet] = FatJet_msoftdrop_jesPileUpMuZeroUp[ijet];
	   FatJet_msoftdrop_jes_up[29][fjet] = FatJet_msoftdrop_jesPileUpEnvelopeUp[ijet];
	   FatJet_msoftdrop_jes_up[30][fjet] = FatJet_msoftdrop_jesSubTotalPileUpUp[ijet];
	   FatJet_msoftdrop_jes_up[31][fjet] = FatJet_msoftdrop_jesSubTotalRelativeUp[ijet];
	   FatJet_msoftdrop_jes_up[32][fjet] = FatJet_msoftdrop_jesSubTotalPtUp[ijet];
	   FatJet_msoftdrop_jes_up[33][fjet] = FatJet_msoftdrop_jesSubTotalScaleUp[ijet];
	   FatJet_msoftdrop_jes_up[34][fjet] = FatJet_msoftdrop_jesSubTotalAbsoluteUp[ijet];
	   FatJet_msoftdrop_jes_up[35][fjet] = FatJet_msoftdrop_jesSubTotalMCUp[ijet];
	   FatJet_msoftdrop_jes_up[36][fjet] = FatJet_msoftdrop_jesTotalUp[ijet];
	   
	   FatJet_pt_jes_dn[0][fjet] = FatJet_pt_jesAbsoluteStatDown[ijet];
	   FatJet_pt_jes_dn[1][fjet] = FatJet_pt_jesAbsoluteScaleDown[ijet];
	   FatJet_pt_jes_dn[2][fjet] = FatJet_pt_jesAbsoluteFlavMapDown[ijet];
	   FatJet_pt_jes_dn[3][fjet] = FatJet_pt_jesAbsoluteMPFBiasDown[ijet];
	   FatJet_pt_jes_dn[4][fjet] = FatJet_pt_jesFragmentationDown[ijet];
	   FatJet_pt_jes_dn[5][fjet] = FatJet_pt_jesSinglePionECALDown[ijet];
	   FatJet_pt_jes_dn[6][fjet] = FatJet_pt_jesSinglePionHCALDown[ijet];
	   FatJet_pt_jes_dn[7][fjet] = FatJet_pt_jesFlavorQCDDown[ijet];
	   FatJet_pt_jes_dn[8][fjet] = FatJet_pt_jesTimePtEtaDown[ijet];
	   FatJet_pt_jes_dn[9][fjet] = FatJet_pt_jesRelativeJEREC1Down[ijet];
	   FatJet_pt_jes_dn[10][fjet] = FatJet_pt_jesRelativeJEREC2Down[ijet];
	   FatJet_pt_jes_dn[11][fjet] = FatJet_pt_jesRelativeJERHFDown[ijet];
	   FatJet_pt_jes_dn[12][fjet] = FatJet_pt_jesRelativePtBBDown[ijet];
	   FatJet_pt_jes_dn[13][fjet] = FatJet_pt_jesRelativePtEC1Down[ijet];
	   FatJet_pt_jes_dn[14][fjet] = FatJet_pt_jesRelativePtEC2Down[ijet];
	   FatJet_pt_jes_dn[15][fjet] = FatJet_pt_jesRelativePtHFDown[ijet];
	   FatJet_pt_jes_dn[16][fjet] = FatJet_pt_jesRelativeBalDown[ijet];
	   FatJet_pt_jes_dn[17][fjet] = FatJet_pt_jesRelativeSampleDown[ijet];
	   FatJet_pt_jes_dn[18][fjet] = FatJet_pt_jesRelativeFSRDown[ijet];
	   FatJet_pt_jes_dn[19][fjet] = FatJet_pt_jesRelativeStatFSRDown[ijet];
	   FatJet_pt_jes_dn[20][fjet] = FatJet_pt_jesRelativeStatECDown[ijet];
	   FatJet_pt_jes_dn[21][fjet] = FatJet_pt_jesRelativeStatHFDown[ijet];
	   FatJet_pt_jes_dn[22][fjet] = FatJet_pt_jesPileUpDataMCDown[ijet];
	   FatJet_pt_jes_dn[23][fjet] = FatJet_pt_jesPileUpPtRefDown[ijet];
	   FatJet_pt_jes_dn[24][fjet] = FatJet_pt_jesPileUpPtBBDown[ijet];
	   FatJet_pt_jes_dn[25][fjet] = FatJet_pt_jesPileUpPtEC1Down[ijet];
	   FatJet_pt_jes_dn[26][fjet] = FatJet_pt_jesPileUpPtEC2Down[ijet];
	   FatJet_pt_jes_dn[27][fjet] = FatJet_pt_jesPileUpPtHFDown[ijet];
	   FatJet_pt_jes_dn[28][fjet] = FatJet_pt_jesPileUpMuZeroDown[ijet];
	   FatJet_pt_jes_dn[29][fjet] = FatJet_pt_jesPileUpEnvelopeDown[ijet];
	   FatJet_pt_jes_dn[30][fjet] = FatJet_pt_jesSubTotalPileUpDown[ijet];
	   FatJet_pt_jes_dn[31][fjet] = FatJet_pt_jesSubTotalRelativeDown[ijet];
	   FatJet_pt_jes_dn[32][fjet] = FatJet_pt_jesSubTotalPtDown[ijet];
	   FatJet_pt_jes_dn[33][fjet] = FatJet_pt_jesSubTotalScaleDown[ijet];
	   FatJet_pt_jes_dn[34][fjet] = FatJet_pt_jesSubTotalAbsoluteDown[ijet];
	   FatJet_pt_jes_dn[35][fjet] = FatJet_pt_jesSubTotalMCDown[ijet];
	   FatJet_pt_jes_dn[36][fjet] = FatJet_pt_jesTotalDown[ijet];
	   
	   FatJet_mass_jes_dn[0][fjet] = FatJet_mass_jesAbsoluteStatDown[ijet];
	   FatJet_mass_jes_dn[1][fjet] = FatJet_mass_jesAbsoluteScaleDown[ijet];
	   FatJet_mass_jes_dn[2][fjet] = FatJet_mass_jesAbsoluteFlavMapDown[ijet];
	   FatJet_mass_jes_dn[3][fjet] = FatJet_mass_jesAbsoluteMPFBiasDown[ijet];
	   FatJet_mass_jes_dn[4][fjet] = FatJet_mass_jesFragmentationDown[ijet];
	   FatJet_mass_jes_dn[5][fjet] = FatJet_mass_jesSinglePionECALDown[ijet];
	   FatJet_mass_jes_dn[6][fjet] = FatJet_mass_jesSinglePionHCALDown[ijet];
	   FatJet_mass_jes_dn[7][fjet] = FatJet_mass_jesFlavorQCDDown[ijet];
	   FatJet_mass_jes_dn[8][fjet] = FatJet_mass_jesTimePtEtaDown[ijet];
	   FatJet_mass_jes_dn[9][fjet] = FatJet_mass_jesRelativeJEREC1Down[ijet];
	   FatJet_mass_jes_dn[10][fjet] = FatJet_mass_jesRelativeJEREC2Down[ijet];
	   FatJet_mass_jes_dn[11][fjet] = FatJet_mass_jesRelativeJERHFDown[ijet];
	   FatJet_mass_jes_dn[12][fjet] = FatJet_mass_jesRelativePtBBDown[ijet];
	   FatJet_mass_jes_dn[13][fjet] = FatJet_mass_jesRelativePtEC1Down[ijet];
	   FatJet_mass_jes_dn[14][fjet] = FatJet_mass_jesRelativePtEC2Down[ijet];
	   FatJet_mass_jes_dn[15][fjet] = FatJet_mass_jesRelativePtHFDown[ijet];
	   FatJet_mass_jes_dn[16][fjet] = FatJet_mass_jesRelativeBalDown[ijet];
	   FatJet_mass_jes_dn[17][fjet] = FatJet_mass_jesRelativeSampleDown[ijet];
	   FatJet_mass_jes_dn[18][fjet] = FatJet_mass_jesRelativeFSRDown[ijet];
	   FatJet_mass_jes_dn[19][fjet] = FatJet_mass_jesRelativeStatFSRDown[ijet];
	   FatJet_mass_jes_dn[20][fjet] = FatJet_mass_jesRelativeStatECDown[ijet];
	   FatJet_mass_jes_dn[21][fjet] = FatJet_mass_jesRelativeStatHFDown[ijet];
	   FatJet_mass_jes_dn[22][fjet] = FatJet_mass_jesPileUpDataMCDown[ijet];
	   FatJet_mass_jes_dn[23][fjet] = FatJet_mass_jesPileUpPtRefDown[ijet];
	   FatJet_mass_jes_dn[24][fjet] = FatJet_mass_jesPileUpPtBBDown[ijet];
	   FatJet_mass_jes_dn[25][fjet] = FatJet_mass_jesPileUpPtEC1Down[ijet];
	   FatJet_mass_jes_dn[26][fjet] = FatJet_mass_jesPileUpPtEC2Down[ijet];
	   FatJet_mass_jes_dn[27][fjet] = FatJet_mass_jesPileUpPtHFDown[ijet];
	   FatJet_mass_jes_dn[28][fjet] = FatJet_mass_jesPileUpMuZeroDown[ijet];
	   FatJet_mass_jes_dn[29][fjet] = FatJet_mass_jesPileUpEnvelopeDown[ijet];
	   FatJet_mass_jes_dn[30][fjet] = FatJet_mass_jesSubTotalPileUpDown[ijet];
	   FatJet_mass_jes_dn[31][fjet] = FatJet_mass_jesSubTotalRelativeDown[ijet];
	   FatJet_mass_jes_dn[32][fjet] = FatJet_mass_jesSubTotalPtDown[ijet];
	   FatJet_mass_jes_dn[33][fjet] = FatJet_mass_jesSubTotalScaleDown[ijet];
	   FatJet_mass_jes_dn[34][fjet] = FatJet_mass_jesSubTotalAbsoluteDown[ijet];
	   FatJet_mass_jes_dn[35][fjet] = FatJet_mass_jesSubTotalMCDown[ijet];
	   FatJet_mass_jes_dn[36][fjet] = FatJet_mass_jesTotalDown[ijet];
	   
	   
	   FatJet_msoftdrop_jes_dn[0][fjet] = FatJet_msoftdrop_jesAbsoluteStatDown[ijet];
	   FatJet_msoftdrop_jes_dn[1][fjet] = FatJet_msoftdrop_jesAbsoluteScaleDown[ijet];
	   FatJet_msoftdrop_jes_dn[2][fjet] = FatJet_msoftdrop_jesAbsoluteFlavMapDown[ijet];
	   FatJet_msoftdrop_jes_dn[3][fjet] = FatJet_msoftdrop_jesAbsoluteMPFBiasDown[ijet];
	   FatJet_msoftdrop_jes_dn[4][fjet] = FatJet_msoftdrop_jesFragmentationDown[ijet];
	   FatJet_msoftdrop_jes_dn[5][fjet] = FatJet_msoftdrop_jesSinglePionECALDown[ijet];
	   FatJet_msoftdrop_jes_dn[6][fjet] = FatJet_msoftdrop_jesSinglePionHCALDown[ijet];
	   FatJet_msoftdrop_jes_dn[7][fjet] = FatJet_msoftdrop_jesFlavorQCDDown[ijet];
	   FatJet_msoftdrop_jes_dn[8][fjet] = FatJet_msoftdrop_jesTimePtEtaDown[ijet];
	   FatJet_msoftdrop_jes_dn[9][fjet] = FatJet_msoftdrop_jesRelativeJEREC1Down[ijet];
	   FatJet_msoftdrop_jes_dn[10][fjet] = FatJet_msoftdrop_jesRelativeJEREC2Down[ijet];
	   FatJet_msoftdrop_jes_dn[11][fjet] = FatJet_msoftdrop_jesRelativeJERHFDown[ijet];
	   FatJet_msoftdrop_jes_dn[12][fjet] = FatJet_msoftdrop_jesRelativePtBBDown[ijet];
	   FatJet_msoftdrop_jes_dn[13][fjet] = FatJet_msoftdrop_jesRelativePtEC1Down[ijet];
	   FatJet_msoftdrop_jes_dn[14][fjet] = FatJet_msoftdrop_jesRelativePtEC2Down[ijet];
	   FatJet_msoftdrop_jes_dn[15][fjet] = FatJet_msoftdrop_jesRelativePtHFDown[ijet];
	   FatJet_msoftdrop_jes_dn[16][fjet] = FatJet_msoftdrop_jesRelativeBalDown[ijet];
	   FatJet_msoftdrop_jes_dn[17][fjet] = FatJet_msoftdrop_jesRelativeSampleDown[ijet];
	   FatJet_msoftdrop_jes_dn[18][fjet] = FatJet_msoftdrop_jesRelativeFSRDown[ijet];
	   FatJet_msoftdrop_jes_dn[19][fjet] = FatJet_msoftdrop_jesRelativeStatFSRDown[ijet];
	   FatJet_msoftdrop_jes_dn[20][fjet] = FatJet_msoftdrop_jesRelativeStatECDown[ijet];
	   FatJet_msoftdrop_jes_dn[21][fjet] = FatJet_msoftdrop_jesRelativeStatHFDown[ijet];
	   FatJet_msoftdrop_jes_dn[22][fjet] = FatJet_msoftdrop_jesPileUpDataMCDown[ijet];
	   FatJet_msoftdrop_jes_dn[23][fjet] = FatJet_msoftdrop_jesPileUpPtRefDown[ijet];
	   FatJet_msoftdrop_jes_dn[24][fjet] = FatJet_msoftdrop_jesPileUpPtBBDown[ijet];
	   FatJet_msoftdrop_jes_dn[25][fjet] = FatJet_msoftdrop_jesPileUpPtEC1Down[ijet];
	   FatJet_msoftdrop_jes_dn[26][fjet] = FatJet_msoftdrop_jesPileUpPtEC2Down[ijet];
	   FatJet_msoftdrop_jes_dn[27][fjet] = FatJet_msoftdrop_jesPileUpPtHFDown[ijet];
	   FatJet_msoftdrop_jes_dn[28][fjet] = FatJet_msoftdrop_jesPileUpMuZeroDown[ijet];
	   FatJet_msoftdrop_jes_dn[29][fjet] = FatJet_msoftdrop_jesPileUpEnvelopeDown[ijet];
	   FatJet_msoftdrop_jes_dn[30][fjet] = FatJet_msoftdrop_jesSubTotalPileUpDown[ijet];
	   FatJet_msoftdrop_jes_dn[31][fjet] = FatJet_msoftdrop_jesSubTotalRelativeDown[ijet];
	   FatJet_msoftdrop_jes_dn[32][fjet] = FatJet_msoftdrop_jesSubTotalPtDown[ijet];
	   FatJet_msoftdrop_jes_dn[33][fjet] = FatJet_msoftdrop_jesSubTotalScaleDown[ijet];
	   FatJet_msoftdrop_jes_dn[34][fjet] = FatJet_msoftdrop_jesSubTotalAbsoluteDown[ijet];
	   FatJet_msoftdrop_jes_dn[35][fjet] = FatJet_msoftdrop_jesSubTotalMCDown[ijet];
	   FatJet_msoftdrop_jes_dn[36][fjet] = FatJet_msoftdrop_jesTotalDown[ijet];
	   
	   // jes ends
	   
	   // systematics end //
	   
	   fjet++;
	   if(fjet>=njetmax) break;
	   }
   
   nFatJet = fjet;
   
   
   // sort by pt AK8 //
   
   for(unsigned ijet=0; ijet < nFatJet; ijet++){
	 for(unsigned kjet=(ijet+1); kjet < nFatJet; kjet++){
		   
		if(FatJet_pt[kjet]>FatJet_pt[ijet]){
			
			float tmppt;
			
			tmppt = FatJet_pt[ijet];
			FatJet_pt[ijet] = FatJet_pt[kjet];
			FatJet_pt[kjet] = tmppt;
		    
		    tmppt = FatJet_area[ijet];
			FatJet_area[ijet] = FatJet_area[kjet];
			FatJet_area[kjet] = tmppt;
			
			tmppt = FatJet_btagCMVA[ijet];
			FatJet_btagCMVA[ijet] = FatJet_btagCMVA[kjet];
			FatJet_btagCMVA[kjet] = tmppt;
			
			tmppt = FatJet_btagCSVV2[ijet];
			FatJet_btagCSVV2[ijet] = FatJet_btagCSVV2[kjet];
			FatJet_btagCSVV2[kjet] = tmppt;
			
			tmppt = FatJet_btagDeepB[ijet];
			FatJet_btagDeepB[ijet] = FatJet_btagDeepB[kjet];
			FatJet_btagDeepB[kjet] = tmppt;
			
			tmppt = FatJet_btagHbb[ijet];
			FatJet_btagHbb[ijet] = FatJet_btagHbb[kjet];
			FatJet_btagHbb[kjet] = tmppt;
			
			tmppt = FatJet_deepTagMD_H4qvsQCD[ijet];
			FatJet_deepTagMD_H4qvsQCD[ijet] = FatJet_deepTagMD_H4qvsQCD[kjet];
			FatJet_deepTagMD_H4qvsQCD[kjet] = tmppt;
			
			tmppt = FatJet_deepTagMD_HbbvsQCD[ijet];
			FatJet_deepTagMD_HbbvsQCD[ijet] = FatJet_deepTagMD_HbbvsQCD[kjet];
			FatJet_deepTagMD_HbbvsQCD[kjet] = tmppt;
			
			tmppt = FatJet_deepTagMD_TvsQCD[ijet];
			FatJet_deepTagMD_TvsQCD[ijet] = FatJet_deepTagMD_TvsQCD[kjet];
			FatJet_deepTagMD_TvsQCD[kjet] = tmppt;
			
			tmppt = FatJet_deepTagMD_WvsQCD[ijet];
			FatJet_deepTagMD_WvsQCD[ijet] = FatJet_deepTagMD_WvsQCD[kjet];
			FatJet_deepTagMD_WvsQCD[kjet] = tmppt;
			
			tmppt = FatJet_deepTagMD_ZHbbvsQCD[ijet];
			FatJet_deepTagMD_ZHbbvsQCD[ijet] = FatJet_deepTagMD_ZHbbvsQCD[kjet];
			FatJet_deepTagMD_ZHbbvsQCD[kjet] = tmppt;
			
			tmppt = FatJet_deepTagMD_ZHccvsQCD[ijet];
			FatJet_deepTagMD_ZHccvsQCD[ijet] = FatJet_deepTagMD_ZHccvsQCD[kjet];
			FatJet_deepTagMD_ZHccvsQCD[kjet] = tmppt;
			
			tmppt = FatJet_deepTagMD_ZbbvsQCD[ijet];
			FatJet_deepTagMD_ZbbvsQCD[ijet] = FatJet_deepTagMD_ZbbvsQCD[kjet];
			FatJet_deepTagMD_ZbbvsQCD[kjet] = tmppt;
			
			tmppt = FatJet_deepTagMD_ZvsQCD[ijet];
			FatJet_deepTagMD_ZvsQCD[ijet] = FatJet_deepTagMD_ZvsQCD[kjet];
			FatJet_deepTagMD_ZvsQCD[kjet] = tmppt;
			
			tmppt = FatJet_deepTagMD_bbvsLight[ijet];
			FatJet_deepTagMD_bbvsLight[ijet] = FatJet_deepTagMD_bbvsLight[kjet];
			FatJet_deepTagMD_bbvsLight[kjet] = tmppt;
			
			tmppt = FatJet_deepTagMD_ccvsLight[ijet];
			FatJet_deepTagMD_ccvsLight[ijet] = FatJet_deepTagMD_ccvsLight[kjet];
			FatJet_deepTagMD_ccvsLight[kjet] = tmppt;
			
			tmppt = FatJet_deepTag_TvsQCD[ijet];
			FatJet_deepTag_TvsQCD[ijet] = FatJet_deepTag_TvsQCD[kjet];
			FatJet_deepTag_TvsQCD[kjet] = tmppt;
			
			tmppt = FatJet_deepTag_WvsQCD[ijet];
			FatJet_deepTag_WvsQCD[ijet] = FatJet_deepTag_WvsQCD[kjet];
			FatJet_deepTag_WvsQCD[kjet] = tmppt;
			
			tmppt = FatJet_deepTag_ZvsQCD[ijet];
			FatJet_deepTag_ZvsQCD[ijet] = FatJet_deepTag_ZvsQCD[kjet];
			FatJet_deepTag_ZvsQCD[kjet] = tmppt;
			
			tmppt = FatJet_eta[ijet];
			FatJet_eta[ijet] = FatJet_eta[kjet];
			FatJet_eta[kjet] = tmppt;
			
			tmppt = FatJet_mass[ijet];
			FatJet_mass[ijet] = FatJet_mass[kjet];
			FatJet_mass[kjet] = tmppt;
			
			tmppt = FatJet_msoftdrop[ijet];
			FatJet_msoftdrop[ijet] = FatJet_msoftdrop[kjet];
			FatJet_msoftdrop[kjet] = tmppt;
			
			tmppt = FatJet_phi[ijet];
			FatJet_phi[ijet] = FatJet_phi[kjet];
			FatJet_phi[kjet] = tmppt;
			
			tmppt = FatJet_rawFactor[ijet];
			FatJet_rawFactor[ijet] = FatJet_rawFactor[kjet];
			FatJet_rawFactor[kjet] = tmppt;
			
			tmppt = FatJet_n2b1[ijet];
			FatJet_n2b1[ijet] = FatJet_n2b1[kjet];
			FatJet_n2b1[kjet] = tmppt;
			
			tmppt = FatJet_n3b1[ijet];
			FatJet_n3b1[ijet] = FatJet_n3b1[kjet];
			FatJet_n3b1[kjet] = tmppt;
			
			tmppt = FatJet_tau1[ijet];
			FatJet_tau1[ijet] = FatJet_tau1[kjet];
			FatJet_tau1[kjet] = tmppt;
			
			tmppt = FatJet_tau2[ijet];
			FatJet_tau2[ijet] = FatJet_tau2[kjet];
			FatJet_tau2[kjet] = tmppt;
			
			tmppt = FatJet_tau3[ijet];
			FatJet_tau3[ijet] = FatJet_tau3[kjet];
			FatJet_tau3[kjet] = tmppt;
			
			tmppt = FatJet_tau4[ijet];
			FatJet_tau4[ijet] = FatJet_tau4[kjet];
			FatJet_tau4[kjet] = tmppt;
			
			tmppt = FatJet_GenJetpt[ijet];
			FatJet_GenJetpt[ijet] = FatJet_GenJetpt[kjet];
			FatJet_GenJetpt[kjet] = tmppt;
			
			tmppt = FatJet_GenJeteta[ijet];
			FatJet_GenJeteta[ijet] = FatJet_GenJeteta[kjet];
			FatJet_GenJeteta[kjet] = tmppt;
			
			tmppt = FatJet_GenJetphi[ijet];
			FatJet_GenJetphi[ijet] = FatJet_GenJetphi[kjet];
			FatJet_GenJetphi[kjet] = tmppt;
			
			tmppt = FatJet_GenJetmass[ijet];
			FatJet_GenJetmass[ijet] = FatJet_GenJetphi[kjet];
			FatJet_GenJetmass[kjet] = tmppt;
			
			// jes
			
			for(int ijes=0; ijes<njesmax; ijes++){
				
				tmppt = FatJet_pt_jes_up[ijes][ijet];
				FatJet_pt_jes_up[ijes][ijet] = FatJet_pt_jes_up[ijes][kjet];
				FatJet_pt_jes_up[ijes][kjet] = tmppt;
				
				tmppt = FatJet_mass_jes_up[ijes][ijet];
				FatJet_mass_jes_up[ijes][ijet] = FatJet_mass_jes_up[ijes][kjet];
				FatJet_mass_jes_up[ijes][kjet] = tmppt;
				
				tmppt = FatJet_msoftdrop_jes_up[ijes][ijet];
				FatJet_msoftdrop_jes_up[ijes][ijet] = FatJet_msoftdrop_jes_up[ijes][kjet];
				FatJet_msoftdrop_jes_up[ijes][kjet] = tmppt;
				
				tmppt = FatJet_pt_jes_dn[ijes][ijet];
				FatJet_pt_jes_dn[ijes][ijet] = FatJet_pt_jes_dn[ijes][kjet];
				FatJet_pt_jes_dn[ijes][kjet] = tmppt;
				
				tmppt = FatJet_mass_jes_dn[ijes][ijet];
				FatJet_mass_jes_dn[ijes][ijet] = FatJet_mass_jes_dn[ijes][kjet];
				FatJet_mass_jes_dn[ijes][kjet] = tmppt;
				
				tmppt = FatJet_msoftdrop_jes_dn[ijes][ijet];
				FatJet_msoftdrop_jes_dn[ijes][ijet] = FatJet_msoftdrop_jes_dn[ijes][kjet];
				FatJet_msoftdrop_jes_dn[ijes][kjet] = tmppt;
				
			}
			
			// jer
			
			tmppt = FatJet_pt_jer_up[ijet];
			FatJet_pt_jer_up[ijet] = FatJet_pt_jer_up[kjet];
			FatJet_pt_jer_up[kjet] = tmppt;
			
			tmppt = FatJet_mass_jer_up[ijet];
			FatJet_mass_jer_up[ijet] = FatJet_mass_jer_up[kjet];
			FatJet_mass_jer_up[kjet] = tmppt;
			
			tmppt = FatJet_msoftdrop_jer_up[ijet];
			FatJet_msoftdrop_jer_up[ijet] = FatJet_msoftdrop_jer_up[kjet];
			FatJet_msoftdrop_jer_up[kjet] = tmppt;
			
			tmppt = FatJet_pt_jer_dn[ijet];
			FatJet_pt_jer_dn[ijet] = FatJet_pt_jer_dn[kjet];
			FatJet_pt_jer_dn[kjet] = tmppt;
			
			tmppt = FatJet_mass_jer_dn[ijet];
			FatJet_mass_jer_dn[ijet] = FatJet_mass_jer_dn[kjet];
			FatJet_mass_jer_dn[kjet] = tmppt;
			
			tmppt = FatJet_msoftdrop_jer_dn[ijet];
			FatJet_msoftdrop_jer_dn[ijet] = FatJet_msoftdrop_jer_dn[kjet];
			FatJet_msoftdrop_jer_dn[kjet] = tmppt;
			
			// jms
			
			tmppt = FatJet_mass_jms_up[ijet];
			FatJet_mass_jms_up[ijet] = FatJet_mass_jms_up[kjet];
			FatJet_mass_jms_up[kjet] = tmppt;
			
			tmppt = FatJet_msoftdrop_jms_up[ijet];
			FatJet_msoftdrop_jms_up[ijet] = FatJet_msoftdrop_jms_up[kjet];
			FatJet_msoftdrop_jms_up[kjet] = tmppt;
			
			tmppt = FatJet_mass_jms_dn[ijet];
			FatJet_mass_jms_dn[ijet] = FatJet_mass_jms_dn[kjet];
			FatJet_mass_jms_dn[kjet] = tmppt;
			
			tmppt = FatJet_msoftdrop_jms_dn[ijet];
			FatJet_msoftdrop_jms_dn[ijet] = FatJet_msoftdrop_jms_dn[kjet];
			FatJet_msoftdrop_jms_dn[kjet] = tmppt;
			
			tmppt = FatJet_mass_jmr_up[ijet];
			FatJet_mass_jmr_up[ijet] = FatJet_mass_jmr_up[kjet];
			FatJet_mass_jmr_up[kjet] = tmppt;
			
			tmppt = FatJet_msoftdrop_jmr_up[ijet];
			FatJet_msoftdrop_jmr_up[ijet] = FatJet_msoftdrop_jmr_up[kjet];
			FatJet_msoftdrop_jmr_up[kjet] = tmppt;
			
			tmppt = FatJet_mass_jmr_dn[ijet];
			FatJet_mass_jmr_dn[ijet] = FatJet_mass_jmr_dn[kjet];
			FatJet_mass_jmr_dn[kjet] = tmppt;
			
			tmppt = FatJet_msoftdrop_jmr_dn[ijet];
			FatJet_msoftdrop_jmr_dn[ijet] = FatJet_msoftdrop_jmr_dn[kjet];
			FatJet_msoftdrop_jmr_dn[kjet] = tmppt;
			
			int tmpid = FatJet_jetId[ijet];
			FatJet_jetId[ijet] = FatJet_jetId[kjet];
			FatJet_jetId[kjet] = tmpid;
			
			tmpid = FatJet_subJetIdx1[ijet];
			FatJet_subJetIdx1[ijet] = FatJet_subJetIdx1[kjet];
			FatJet_subJetIdx1[kjet] = tmpid;
			
			tmpid = FatJet_subJetIdx2[ijet];
			FatJet_subJetIdx2[ijet] = FatJet_subJetIdx2[kjet];
			FatJet_subJetIdx2[kjet] = tmpid;
			
			tmpid = FatJet_MatchGenJet[ijet];
			FatJet_MatchGenJet[ijet] = FatJet_MatchGenJet[kjet];
			FatJet_MatchGenJet[kjet] = tmpid;
		     
		     }
		   
		}//kjet
	}//ijet
   
   // pt sorting ends for AK8 jet
   
   if(nFatJet > 0){
   for(unsigned ijet=0; ijet<nFatJet; ijet++){

	  float maxsubjetbtagvalue_i = -100;
	  
	  for(unsigned isub=0; isub<nSubJet; isub++){
				if((int(isub)!=FatJet_subJetIdx1[ijet])&&(int(isub)!=FatJet_subJetIdx2[ijet])) continue;	
				if(SubJet_btagCSVV2[isub] > maxsubjetbtagvalue_i) {   
						maxsubjetbtagvalue_i = SubJet_btagCSVV2[isub]; 
						FatJet_subbtagId_CSVv2[ijet] = int(isub);
						}
			}
	    FatJet_subbtagCSVV2[ijet] = maxsubjetbtagvalue_i;
	  
	  float maxsubjetdeepBvalue_i = -100;
	  
	  for(unsigned isub=0; isub<nSubJet; isub++){
				if((int(isub)!=FatJet_subJetIdx1[ijet])&&(int(isub)!=FatJet_subJetIdx2[ijet])) continue;	
				if(SubJet_btagDeepB[isub] > maxsubjetdeepBvalue_i) {  
					 maxsubjetdeepBvalue_i = SubJet_btagDeepB[isub]; 
					 FatJet_subbtagId_DeepCSV[ijet] = int(isub);
					 }
			}
		FatJet_subbtagDeepB[ijet] = maxsubjetdeepBvalue_i;
	  
	  }
   }
 
   fjet = 0;
   
   for(unsigned ijet=0; ijet<nJet; ijet++){
	
	   if(ijet==0 && (Jet_pt[ijet]<150 || abs(Jet_eta[ijet])>3.5)) break;
 
	   if(Jet_jetId[ijet]==0) continue;
	
	   int jetid = *(decToBinary(Jet_jetId[ijet])+1);//((decToBinary(Jet_jetId[ijet]))/10)%10; // tight id
	   if(jetid!=1) continue;
	   
	   if(Jet_pt[ijet]<30.) continue;
//	   if(Jet_pt[ijet]< 0.9*AK4ptcut_fi) continue;
	   if(fabs(Jet_eta[ijet])>jeteta_cut) continue;
	  
	   bool mutag = false;
	   for(unsigned imu=0; imu<nMuon; imu++){
		   if(delta2R(Muon_eta[imu],Muon_phi[imu],Jet_eta[ijet],Jet_phi[ijet])<0.4){
			   mutag = true;
			   break;
			   }
		   }
	   
	   bool etag = false;
	   for(unsigned iel=0; iel<nElectron; iel++){
		   if(delta2R(Electron_eta[iel],Electron_phi[iel],Jet_eta[ijet],Jet_phi[ijet])<0.4){
			   etag = true;
			   break;
			   }
		   }
	   
	   if(mutag||etag) continue;
	   
	   Jet_area[fjet] = Jet_area[ijet];
	   Jet_btagCMVA[fjet] = Jet_btagCMVA[ijet];
	   Jet_btagCSVV2[fjet] = Jet_btagCSVV2[ijet];
	   Jet_btagDeepB[fjet] = Jet_btagDeepB[ijet];
	   Jet_btagDeepC[fjet] = Jet_btagDeepC[ijet];
	   Jet_btagDeepFlavB[fjet] = Jet_btagDeepFlavB[ijet];
	   Jet_chEmEF[fjet] = Jet_chEmEF[ijet];
	   Jet_chHEF[fjet] = Jet_chHEF[ijet];
	   Jet_eta[fjet] = Jet_eta[ijet];
	   #ifdef NoJEC
	   Jet_mass[fjet] = Jet_mass[ijet];
	   Jet_pt[fjet] = Jet_pt[ijet];
	   #else
	   Jet_mass[fjet] = Jet_mass_nom[ijet];
	   Jet_pt[fjet] = Jet_pt_nom[ijet];
	   #endif
	   Jet_muEF[fjet] = Jet_muEF[ijet];
	   Jet_neEmEF[fjet] = Jet_neEmEF[ijet];
	   Jet_neHEF[fjet] = Jet_neHEF[ijet];
	   Jet_phi[fjet] = Jet_phi[ijet];
	   Jet_qgl[fjet] = Jet_qgl[ijet];
	   Jet_rawFactor[fjet] = Jet_rawFactor[ijet];
	   Jet_bRegCorr[fjet] = Jet_bRegCorr[ijet];
	   Jet_bRegRes[fjet] = Jet_bRegRes[ijet];
	   Jet_electronIdx1[fjet] = Jet_electronIdx1[ijet];
	   Jet_electronIdx2[fjet] = Jet_electronIdx2[ijet];
	   Jet_jetId[fjet] = Jet_jetId[ijet];
	   Jet_muonIdx1[fjet] = Jet_muonIdx1[ijet];
	   Jet_muonIdx2[fjet] = Jet_muonIdx2[ijet];
	   Jet_nConstituents[fjet] = Jet_nConstituents[ijet];
	   Jet_nElectrons[fjet] = Jet_nElectrons[ijet];
	   Jet_nMuons[fjet] = Jet_nMuons[ijet];
	   Jet_puId[fjet] = Jet_puId[ijet];
	   if(isMC){
	   Jet_genJetIdx[fjet] = Jet_genJetIdx[ijet];
	   Jet_partonFlavour[fjet] = Jet_partonFlavour[ijet];
	   Jet_hadronFlavour[fjet] = Jet_hadronFlavour[ijet];
	   }
//	   Jet_DeepCSVMCeff[fjet] = 1;//BTag_MCEfficiency("deepcsvm",int(Jet_hadronFlavour[fjet]),Jet_pt[fjet],fabs(Jet_eta[fjet]));
	
	   Jet_deepflavbtagSF_shape[fjet] = Jet_deepflavbtagSF_shape[ijet];
//	   Jet_deepCSVbtagSF_shape[fjet] = Jet_deepCSVbtagSF_shape[ijet];
//	   Jet_CSVbtagSF_shape[fjet] = Jet_CSVbtagSF_shape[ijet];

	   #ifdef Anal_2018
	   if(isMC){
		   if(Jet_phi[fjet]>-1.57 && Jet_phi[fjet]<-0.87){
			   if(Jet_eta[fjet]>-2.5 && Jet_eta[fjet]<-1.3) { 
		  //	   scale down jES by 20 % for jets with -1.57 <phi< -0.87 and -2.5<eta<-1.3
				   Jet_mass[fjet] = 0.8*Jet_mass[fjet];
				   Jet_pt[fjet] = 0.8*Jet_pt[fjet];
				   }
			   if(Jet_eta[fjet]>-3.0 && Jet_eta[fjet]<-2.5) {
		//	   	   scale down jES by 35 % for jets with -1.57 <phi< -0.87 and -3.0<eta<-2.5	   
				   Jet_mass[fjet] = 0.65*Jet_mass[fjet];
				   Jet_pt[fjet] = 0.65*Jet_pt[fjet];
				   }	   
			   }
		   }
	   #endif
	    
	   // systematics //
	   // jer
	   Jet_pt_jer_up[fjet] = Jet_pt_jerUp[ijet];
	   Jet_mass_jer_up[fjet] = Jet_mass_jerUp[ijet];
	   
	   Jet_pt_jer_dn[fjet] = Jet_pt_jerDown[ijet];
	   Jet_mass_jer_dn[fjet] = Jet_mass_jerDown[ijet];
	   //jer ends
	   
	   // jms & jmr //
	   
	   Jet_mass_jms_up[fjet] = Jet_mass_jmsUp[ijet]; //Jet_mass_jmsUp[ijet];
	   Jet_mass_jmr_up[fjet] = Jet_mass_jmrUp[ijet]; //Jet_mass_jmrUp[ijet];
	   
	   Jet_mass_jms_dn[fjet] = Jet_mass_jmsDown[ijet]; //Jet_mass_jmsDown[ijet];
	   Jet_mass_jmr_dn[fjet] = Jet_mass_jmrDown[ijet]; //Jet_mass_jmrDown[ijet];
	   
	   // jms & jmr ends//
	
	   // jes //
	   Jet_pt_jes_up[0][fjet] = Jet_pt_jesAbsoluteStatUp[ijet];
	   Jet_pt_jes_up[1][fjet] = Jet_pt_jesAbsoluteScaleUp[ijet];
	   Jet_pt_jes_up[2][fjet] = Jet_pt_jesAbsoluteFlavMapUp[ijet];
	   Jet_pt_jes_up[3][fjet] = Jet_pt_jesAbsoluteMPFBiasUp[ijet];
	   Jet_pt_jes_up[4][fjet] = Jet_pt_jesFragmentationUp[ijet];
	   Jet_pt_jes_up[5][fjet] = Jet_pt_jesSinglePionECALUp[ijet];
	   Jet_pt_jes_up[6][fjet] = Jet_pt_jesSinglePionHCALUp[ijet];
	   Jet_pt_jes_up[7][fjet] = Jet_pt_jesFlavorQCDUp[ijet];
	   Jet_pt_jes_up[8][fjet] = Jet_pt_jesTimePtEtaUp[ijet];
	   Jet_pt_jes_up[9][fjet] = Jet_pt_jesRelativeJEREC1Up[ijet];
	   Jet_pt_jes_up[10][fjet] = Jet_pt_jesRelativeJEREC2Up[ijet];
	   Jet_pt_jes_up[11][fjet] = Jet_pt_jesRelativeJERHFUp[ijet];
	   Jet_pt_jes_up[12][fjet] = Jet_pt_jesRelativePtBBUp[ijet];
	   Jet_pt_jes_up[13][fjet] = Jet_pt_jesRelativePtEC1Up[ijet];
	   Jet_pt_jes_up[14][fjet] = Jet_pt_jesRelativePtEC2Up[ijet];
	   Jet_pt_jes_up[15][fjet] = Jet_pt_jesRelativePtHFUp[ijet];
	   Jet_pt_jes_up[16][fjet] = Jet_pt_jesRelativeBalUp[ijet];
	   Jet_pt_jes_up[17][fjet] = Jet_pt_jesRelativeSampleUp[ijet];
	   Jet_pt_jes_up[18][fjet] = Jet_pt_jesRelativeFSRUp[ijet];
	   Jet_pt_jes_up[19][fjet] = Jet_pt_jesRelativeStatFSRUp[ijet];
	   Jet_pt_jes_up[20][fjet] = Jet_pt_jesRelativeStatECUp[ijet];
	   Jet_pt_jes_up[21][fjet] = Jet_pt_jesRelativeStatHFUp[ijet];
	   Jet_pt_jes_up[22][fjet] = Jet_pt_jesPileUpDataMCUp[ijet];
	   Jet_pt_jes_up[23][fjet] = Jet_pt_jesPileUpPtRefUp[ijet];
	   Jet_pt_jes_up[24][fjet] = Jet_pt_jesPileUpPtBBUp[ijet];
	   Jet_pt_jes_up[25][fjet] = Jet_pt_jesPileUpPtEC1Up[ijet];
	   Jet_pt_jes_up[26][fjet] = Jet_pt_jesPileUpPtEC2Up[ijet];
	   Jet_pt_jes_up[27][fjet] = Jet_pt_jesPileUpPtHFUp[ijet];
	   Jet_pt_jes_up[28][fjet] = Jet_pt_jesPileUpMuZeroUp[ijet];
	   Jet_pt_jes_up[29][fjet] = Jet_pt_jesPileUpEnvelopeUp[ijet];
	   Jet_pt_jes_up[30][fjet] = Jet_pt_jesSubTotalPileUpUp[ijet];
	   Jet_pt_jes_up[31][fjet] = Jet_pt_jesSubTotalRelativeUp[ijet];
	   Jet_pt_jes_up[32][fjet] = Jet_pt_jesSubTotalPtUp[ijet];
	   Jet_pt_jes_up[33][fjet] = Jet_pt_jesSubTotalScaleUp[ijet];
	   Jet_pt_jes_up[34][fjet] = Jet_pt_jesSubTotalAbsoluteUp[ijet];
	   Jet_pt_jes_up[35][fjet] = Jet_pt_jesSubTotalMCUp[ijet];
	   Jet_pt_jes_up[36][fjet] = Jet_pt_jesTotalUp[ijet];
	   
	   Jet_mass_jes_up[0][fjet] = Jet_mass_jesAbsoluteStatUp[ijet];
	   Jet_mass_jes_up[1][fjet] = Jet_mass_jesAbsoluteScaleUp[ijet];
	   Jet_mass_jes_up[2][fjet] = Jet_mass_jesAbsoluteFlavMapUp[ijet];
	   Jet_mass_jes_up[3][fjet] = Jet_mass_jesAbsoluteMPFBiasUp[ijet];
	   Jet_mass_jes_up[4][fjet] = Jet_mass_jesFragmentationUp[ijet];
	   Jet_mass_jes_up[5][fjet] = Jet_mass_jesSinglePionECALUp[ijet];
	   Jet_mass_jes_up[6][fjet] = Jet_mass_jesSinglePionHCALUp[ijet];
	   Jet_mass_jes_up[7][fjet] = Jet_mass_jesFlavorQCDUp[ijet];
	   Jet_mass_jes_up[8][fjet] = Jet_mass_jesTimePtEtaUp[ijet];
	   Jet_mass_jes_up[9][fjet] = Jet_mass_jesRelativeJEREC1Up[ijet];
	   Jet_mass_jes_up[10][fjet] = Jet_mass_jesRelativeJEREC2Up[ijet];
	   Jet_mass_jes_up[11][fjet] = Jet_mass_jesRelativeJERHFUp[ijet];
	   Jet_mass_jes_up[12][fjet] = Jet_mass_jesRelativePtBBUp[ijet];
	   Jet_mass_jes_up[13][fjet] = Jet_mass_jesRelativePtEC1Up[ijet];
	   Jet_mass_jes_up[14][fjet] = Jet_mass_jesRelativePtEC2Up[ijet];
	   Jet_mass_jes_up[15][fjet] = Jet_mass_jesRelativePtHFUp[ijet];
	   Jet_mass_jes_up[16][fjet] = Jet_mass_jesRelativeBalUp[ijet];
	   Jet_mass_jes_up[17][fjet] = Jet_mass_jesRelativeSampleUp[ijet];
	   Jet_mass_jes_up[18][fjet] = Jet_mass_jesRelativeFSRUp[ijet];
	   Jet_mass_jes_up[19][fjet] = Jet_mass_jesRelativeStatFSRUp[ijet];
	   Jet_mass_jes_up[20][fjet] = Jet_mass_jesRelativeStatECUp[ijet];
	   Jet_mass_jes_up[21][fjet] = Jet_mass_jesRelativeStatHFUp[ijet];
	   Jet_mass_jes_up[22][fjet] = Jet_mass_jesPileUpDataMCUp[ijet];
	   Jet_mass_jes_up[23][fjet] = Jet_mass_jesPileUpPtRefUp[ijet];
	   Jet_mass_jes_up[24][fjet] = Jet_mass_jesPileUpPtBBUp[ijet];
	   Jet_mass_jes_up[25][fjet] = Jet_mass_jesPileUpPtEC1Up[ijet];
	   Jet_mass_jes_up[26][fjet] = Jet_mass_jesPileUpPtEC2Up[ijet];
	   Jet_mass_jes_up[27][fjet] = Jet_mass_jesPileUpPtHFUp[ijet];
	   Jet_mass_jes_up[28][fjet] = Jet_mass_jesPileUpMuZeroUp[ijet];
	   Jet_mass_jes_up[29][fjet] = Jet_mass_jesPileUpEnvelopeUp[ijet];
	   Jet_mass_jes_up[30][fjet] = Jet_mass_jesSubTotalPileUpUp[ijet];
	   Jet_mass_jes_up[31][fjet] = Jet_mass_jesSubTotalRelativeUp[ijet];
	   Jet_mass_jes_up[32][fjet] = Jet_mass_jesSubTotalPtUp[ijet];
	   Jet_mass_jes_up[33][fjet] = Jet_mass_jesSubTotalScaleUp[ijet];
	   Jet_mass_jes_up[34][fjet] = Jet_mass_jesSubTotalAbsoluteUp[ijet];
	   Jet_mass_jes_up[35][fjet] = Jet_mass_jesSubTotalMCUp[ijet];
	   Jet_mass_jes_up[36][fjet] = Jet_mass_jesTotalUp[ijet];
	   
	   Jet_pt_jes_dn[0][fjet] = Jet_pt_jesAbsoluteStatDown[ijet];
	   Jet_pt_jes_dn[1][fjet] = Jet_pt_jesAbsoluteScaleDown[ijet];
	   Jet_pt_jes_dn[2][fjet] = Jet_pt_jesAbsoluteFlavMapDown[ijet];
	   Jet_pt_jes_dn[3][fjet] = Jet_pt_jesAbsoluteMPFBiasDown[ijet];
	   Jet_pt_jes_dn[4][fjet] = Jet_pt_jesFragmentationDown[ijet];
	   Jet_pt_jes_dn[5][fjet] = Jet_pt_jesSinglePionECALDown[ijet];
	   Jet_pt_jes_dn[6][fjet] = Jet_pt_jesSinglePionHCALDown[ijet];
	   Jet_pt_jes_dn[7][fjet] = Jet_pt_jesFlavorQCDDown[ijet];
	   Jet_pt_jes_dn[8][fjet] = Jet_pt_jesTimePtEtaDown[ijet];
	   Jet_pt_jes_dn[9][fjet] = Jet_pt_jesRelativeJEREC1Down[ijet];
	   Jet_pt_jes_dn[10][fjet] = Jet_pt_jesRelativeJEREC2Down[ijet];
	   Jet_pt_jes_dn[11][fjet] = Jet_pt_jesRelativeJERHFDown[ijet];
	   Jet_pt_jes_dn[12][fjet] = Jet_pt_jesRelativePtBBDown[ijet];
	   Jet_pt_jes_dn[13][fjet] = Jet_pt_jesRelativePtEC1Down[ijet];
	   Jet_pt_jes_dn[14][fjet] = Jet_pt_jesRelativePtEC2Down[ijet];
	   Jet_pt_jes_dn[15][fjet] = Jet_pt_jesRelativePtHFDown[ijet];
	   Jet_pt_jes_dn[16][fjet] = Jet_pt_jesRelativeBalDown[ijet];
	   Jet_pt_jes_dn[17][fjet] = Jet_pt_jesRelativeSampleDown[ijet];
	   Jet_pt_jes_dn[18][fjet] = Jet_pt_jesRelativeFSRDown[ijet];
	   Jet_pt_jes_dn[19][fjet] = Jet_pt_jesRelativeStatFSRDown[ijet];
	   Jet_pt_jes_dn[20][fjet] = Jet_pt_jesRelativeStatECDown[ijet];
	   Jet_pt_jes_dn[21][fjet] = Jet_pt_jesRelativeStatHFDown[ijet];
	   Jet_pt_jes_dn[22][fjet] = Jet_pt_jesPileUpDataMCDown[ijet];
	   Jet_pt_jes_dn[23][fjet] = Jet_pt_jesPileUpPtRefDown[ijet];
	   Jet_pt_jes_dn[24][fjet] = Jet_pt_jesPileUpPtBBDown[ijet];
	   Jet_pt_jes_dn[25][fjet] = Jet_pt_jesPileUpPtEC1Down[ijet];
	   Jet_pt_jes_dn[26][fjet] = Jet_pt_jesPileUpPtEC2Down[ijet];
	   Jet_pt_jes_dn[27][fjet] = Jet_pt_jesPileUpPtHFDown[ijet];
	   Jet_pt_jes_dn[28][fjet] = Jet_pt_jesPileUpMuZeroDown[ijet];
	   Jet_pt_jes_dn[29][fjet] = Jet_pt_jesPileUpEnvelopeDown[ijet];
	   Jet_pt_jes_dn[30][fjet] = Jet_pt_jesSubTotalPileUpDown[ijet];
	   Jet_pt_jes_dn[31][fjet] = Jet_pt_jesSubTotalRelativeDown[ijet];
	   Jet_pt_jes_dn[32][fjet] = Jet_pt_jesSubTotalPtDown[ijet];
	   Jet_pt_jes_dn[33][fjet] = Jet_pt_jesSubTotalScaleDown[ijet];
	   Jet_pt_jes_dn[34][fjet] = Jet_pt_jesSubTotalAbsoluteDown[ijet];
	   Jet_pt_jes_dn[35][fjet] = Jet_pt_jesSubTotalMCDown[ijet];
	   Jet_pt_jes_dn[36][fjet] = Jet_pt_jesTotalDown[ijet];
	   
	   Jet_mass_jes_dn[0][fjet] = Jet_mass_jesAbsoluteStatDown[ijet];
	   Jet_mass_jes_dn[1][fjet] = Jet_mass_jesAbsoluteScaleDown[ijet];
	   Jet_mass_jes_dn[2][fjet] = Jet_mass_jesAbsoluteFlavMapDown[ijet];
	   Jet_mass_jes_dn[3][fjet] = Jet_mass_jesAbsoluteMPFBiasDown[ijet];
	   Jet_mass_jes_dn[4][fjet] = Jet_mass_jesFragmentationDown[ijet];
	   Jet_mass_jes_dn[5][fjet] = Jet_mass_jesSinglePionECALDown[ijet];
	   Jet_mass_jes_dn[6][fjet] = Jet_mass_jesSinglePionHCALDown[ijet];
	   Jet_mass_jes_dn[7][fjet] = Jet_mass_jesFlavorQCDDown[ijet];
	   Jet_mass_jes_dn[8][fjet] = Jet_mass_jesTimePtEtaDown[ijet];
	   Jet_mass_jes_dn[9][fjet] = Jet_mass_jesRelativeJEREC1Down[ijet];
	   Jet_mass_jes_dn[10][fjet] = Jet_mass_jesRelativeJEREC2Down[ijet];
	   Jet_mass_jes_dn[11][fjet] = Jet_mass_jesRelativeJERHFDown[ijet];
	   Jet_mass_jes_dn[12][fjet] = Jet_mass_jesRelativePtBBDown[ijet];
	   Jet_mass_jes_dn[13][fjet] = Jet_mass_jesRelativePtEC1Down[ijet];
	   Jet_mass_jes_dn[14][fjet] = Jet_mass_jesRelativePtEC2Down[ijet];
	   Jet_mass_jes_dn[15][fjet] = Jet_mass_jesRelativePtHFDown[ijet];
	   Jet_mass_jes_dn[16][fjet] = Jet_mass_jesRelativeBalDown[ijet];
	   Jet_mass_jes_dn[17][fjet] = Jet_mass_jesRelativeSampleDown[ijet];
	   Jet_mass_jes_dn[18][fjet] = Jet_mass_jesRelativeFSRDown[ijet];
	   Jet_mass_jes_dn[19][fjet] = Jet_mass_jesRelativeStatFSRDown[ijet];
	   Jet_mass_jes_dn[20][fjet] = Jet_mass_jesRelativeStatECDown[ijet];
	   Jet_mass_jes_dn[21][fjet] = Jet_mass_jesRelativeStatHFDown[ijet];
	   Jet_mass_jes_dn[22][fjet] = Jet_mass_jesPileUpDataMCDown[ijet];
	   Jet_mass_jes_dn[23][fjet] = Jet_mass_jesPileUpPtRefDown[ijet];
	   Jet_mass_jes_dn[24][fjet] = Jet_mass_jesPileUpPtBBDown[ijet];
	   Jet_mass_jes_dn[25][fjet] = Jet_mass_jesPileUpPtEC1Down[ijet];
	   Jet_mass_jes_dn[26][fjet] = Jet_mass_jesPileUpPtEC2Down[ijet];
	   Jet_mass_jes_dn[27][fjet] = Jet_mass_jesPileUpPtHFDown[ijet];
	   Jet_mass_jes_dn[28][fjet] = Jet_mass_jesPileUpMuZeroDown[ijet];
	   Jet_mass_jes_dn[29][fjet] = Jet_mass_jesPileUpEnvelopeDown[ijet];
	   Jet_mass_jes_dn[30][fjet] = Jet_mass_jesSubTotalPileUpDown[ijet];
	   Jet_mass_jes_dn[31][fjet] = Jet_mass_jesSubTotalRelativeDown[ijet];
	   Jet_mass_jes_dn[32][fjet] = Jet_mass_jesSubTotalPtDown[ijet];
	   Jet_mass_jes_dn[33][fjet] = Jet_mass_jesSubTotalScaleDown[ijet];
	   Jet_mass_jes_dn[34][fjet] = Jet_mass_jesSubTotalAbsoluteDown[ijet];
	   Jet_mass_jes_dn[35][fjet] = Jet_mass_jesSubTotalMCDown[ijet];
	   Jet_mass_jes_dn[36][fjet] = Jet_mass_jesTotalDown[ijet];
	   
	   // jes ends
	
	   Jet_deepflavbtagSF_unc_up[0][fjet] = Jet_deepflavbtagSF_shape_up_jes[ijet];
	   Jet_deepflavbtagSF_unc_up[1][fjet] = Jet_deepflavbtagSF_shape_up_lf[ijet];
	   Jet_deepflavbtagSF_unc_up[2][fjet] = Jet_deepflavbtagSF_shape_up_hf[ijet];
	   Jet_deepflavbtagSF_unc_up[3][fjet] = Jet_deepflavbtagSF_shape_up_hfstats1[ijet];
	   Jet_deepflavbtagSF_unc_up[4][fjet] = Jet_deepflavbtagSF_shape_up_hfstats2[ijet];
	   Jet_deepflavbtagSF_unc_up[5][fjet] = Jet_deepflavbtagSF_shape_up_lfstats1[ijet];
	   Jet_deepflavbtagSF_unc_up[6][fjet] = Jet_deepflavbtagSF_shape_up_lfstats2[ijet];
	   Jet_deepflavbtagSF_unc_up[7][fjet] = Jet_deepflavbtagSF_shape_up_cferr1[ijet];
	   Jet_deepflavbtagSF_unc_up[8][fjet] = Jet_deepflavbtagSF_shape_up_cferr2[ijet];
	   
	   Jet_deepflavbtagSF_unc_dn[0][fjet] = Jet_deepflavbtagSF_shape_down_jes[ijet];
	   Jet_deepflavbtagSF_unc_dn[1][fjet] = Jet_deepflavbtagSF_shape_down_lf[ijet];
	   Jet_deepflavbtagSF_unc_dn[2][fjet] = Jet_deepflavbtagSF_shape_down_hf[ijet];
	   Jet_deepflavbtagSF_unc_dn[3][fjet] = Jet_deepflavbtagSF_shape_down_hfstats1[ijet];
	   Jet_deepflavbtagSF_unc_dn[4][fjet] = Jet_deepflavbtagSF_shape_down_hfstats2[ijet];
	   Jet_deepflavbtagSF_unc_dn[5][fjet] = Jet_deepflavbtagSF_shape_down_lfstats1[ijet];
	   Jet_deepflavbtagSF_unc_dn[6][fjet] = Jet_deepflavbtagSF_shape_down_lfstats2[ijet];
	   Jet_deepflavbtagSF_unc_dn[7][fjet] = Jet_deepflavbtagSF_shape_down_cferr1[ijet];
	   Jet_deepflavbtagSF_unc_dn[8][fjet] = Jet_deepflavbtagSF_shape_down_cferr2[ijet];
	   
	   // systematics end //
	      
	   fjet++;
	   if(fjet>=njetmax) break;
	   }
   
   nJet = fjet;
   
   // sort by pt AK4 //
   
   for(unsigned ijet=0; ijet < nJet; ijet++){
	 for(unsigned kjet=(ijet+1); kjet < nJet; kjet++){
		   
		if(Jet_pt[kjet]>Jet_pt[ijet]){
			
			float tmppt;
			
			tmppt = Jet_pt[ijet];
			Jet_pt[ijet] = Jet_pt[kjet];
			Jet_pt[kjet] = tmppt;
		    
		    tmppt = Jet_area[ijet];
			Jet_area[ijet] = Jet_area[kjet];
			Jet_area[kjet] = tmppt;
			
			tmppt = Jet_btagCMVA[ijet];
			Jet_btagCMVA[ijet] = Jet_btagCMVA[kjet];
			Jet_btagCMVA[kjet] = tmppt;
			
			tmppt = Jet_btagCSVV2[ijet];
			Jet_btagCSVV2[ijet] = Jet_btagCSVV2[kjet];
			Jet_btagCSVV2[kjet] = tmppt;
			
			tmppt = Jet_btagDeepB[ijet];
			Jet_btagDeepB[ijet] = Jet_btagDeepB[kjet];
			Jet_btagDeepB[kjet] = tmppt;
			
			tmppt = Jet_btagDeepC[ijet];
			Jet_btagDeepC[ijet] = Jet_btagDeepC[kjet];
			Jet_btagDeepC[kjet] = tmppt;
			
			tmppt = Jet_btagDeepFlavB[ijet];
			Jet_btagDeepFlavB[ijet] = Jet_btagDeepFlavB[kjet];
			Jet_btagDeepFlavB[kjet] = tmppt;
			
			tmppt = Jet_chEmEF[ijet];
			Jet_chEmEF[ijet] = Jet_chEmEF[kjet];
			Jet_chEmEF[kjet] = tmppt;
			
			tmppt = Jet_chHEF[ijet];
			Jet_chHEF[ijet] = Jet_chHEF[kjet];
			Jet_chHEF[kjet] = tmppt;
			
			tmppt = Jet_eta[ijet];
			Jet_eta[ijet] = Jet_eta[kjet];
			Jet_eta[kjet] = tmppt;
			
			tmppt = Jet_mass[ijet];
			Jet_mass[ijet] = Jet_mass[kjet];
			Jet_mass[kjet] = tmppt;
			
			tmppt = Jet_muEF[ijet];
			Jet_muEF[ijet] = Jet_muEF[kjet];
			Jet_muEF[kjet] = tmppt;
			
			tmppt = Jet_neEmEF[ijet];
			Jet_neEmEF[ijet] = Jet_neEmEF[kjet];
			Jet_neEmEF[kjet] = tmppt;
			
			tmppt = Jet_neHEF[ijet];
			Jet_neHEF[ijet] = Jet_neHEF[kjet];
			Jet_neHEF[kjet] = tmppt;
			
			tmppt = Jet_phi[ijet];
			Jet_phi[ijet] = Jet_phi[kjet];
			Jet_phi[kjet] = tmppt;
			
			tmppt = Jet_qgl[ijet];
			Jet_qgl[ijet] = Jet_qgl[kjet];
			Jet_qgl[kjet] = tmppt;
			
			tmppt = Jet_rawFactor[ijet];
			Jet_rawFactor[ijet] = Jet_rawFactor[kjet];
			Jet_rawFactor[kjet] = tmppt;
			
			tmppt = Jet_bRegCorr[ijet];
			Jet_bRegCorr[ijet] = Jet_bRegCorr[kjet];
			Jet_bRegCorr[kjet] = tmppt;
			
			tmppt = Jet_bRegRes[ijet];
			Jet_bRegRes[ijet] = Jet_bRegRes[kjet];
			Jet_bRegRes[kjet] = tmppt;
			
			tmppt = Jet_deepflavbtagSF_shape[ijet];
			Jet_deepflavbtagSF_shape[ijet] = Jet_deepflavbtagSF_shape[kjet];
			Jet_deepflavbtagSF_shape[kjet] = tmppt;
			
			for(int ijes=0; ijes<njesmax; ijes++){
				
				tmppt = Jet_pt_jes_up[ijes][ijet];
				Jet_pt_jes_up[ijes][ijet] = Jet_pt_jes_up[ijes][kjet];
				Jet_pt_jes_up[ijes][kjet] = tmppt;
				
				tmppt = Jet_mass_jes_up[ijes][ijet];
				Jet_mass_jes_up[ijes][ijet] = Jet_mass_jes_up[ijes][kjet];
				Jet_mass_jes_up[ijes][kjet] = tmppt;
				
				tmppt = Jet_pt_jes_dn[ijes][ijet];
				Jet_pt_jes_dn[ijes][ijet] = Jet_pt_jes_dn[ijes][kjet];
				Jet_pt_jes_dn[ijes][kjet] = tmppt;
				
				tmppt = Jet_mass_jes_dn[ijes][ijet];
				Jet_mass_jes_dn[ijes][ijet] = Jet_mass_jes_dn[ijes][kjet];
				Jet_mass_jes_dn[ijes][kjet] = tmppt;
				
			}
			
			for(int ibtag=0; ibtag<nbtagmax; ibtag++){
			
				tmppt = Jet_deepflavbtagSF_unc_up[ibtag][ijet];
				Jet_deepflavbtagSF_unc_up[ibtag][ijet] = Jet_deepflavbtagSF_unc_up[ibtag][kjet];
				Jet_deepflavbtagSF_unc_up[ibtag][kjet] = tmppt;
				
				tmppt = Jet_deepflavbtagSF_unc_dn[ibtag][ijet];
				Jet_deepflavbtagSF_unc_dn[ibtag][ijet] = Jet_deepflavbtagSF_unc_dn[ibtag][kjet];
				Jet_deepflavbtagSF_unc_dn[ibtag][kjet] = tmppt;
			
			}
			
			tmppt = Jet_pt_jer_up[ijet];
			Jet_pt_jer_up[ijet] = Jet_pt_jer_up[kjet];
			Jet_pt_jer_up[kjet] = tmppt;
			
			tmppt = Jet_mass_jer_up[ijet];
			Jet_mass_jer_up[ijet] = Jet_mass_jer_up[kjet];
			Jet_mass_jer_up[kjet] = tmppt;
			
			tmppt = Jet_mass_jms_up[ijet];
			Jet_mass_jms_up[ijet] = Jet_mass_jms_up[kjet];
			Jet_mass_jms_up[kjet] = tmppt;
			
			tmppt = Jet_mass_jmr_up[ijet];
			Jet_mass_jmr_up[ijet] = Jet_mass_jmr_up[kjet];
			Jet_mass_jmr_up[kjet] = tmppt;
			
			tmppt = Jet_pt_jer_dn[ijet];
			Jet_pt_jer_dn[ijet] = Jet_pt_jer_dn[kjet];
			Jet_pt_jer_dn[kjet] = tmppt;
			
			tmppt = Jet_mass_jer_dn[ijet];
			Jet_mass_jer_dn[ijet] = Jet_mass_jer_dn[kjet];
			Jet_mass_jer_dn[kjet] = tmppt;
			
			tmppt = Jet_mass_jms_dn[ijet];
			Jet_mass_jms_dn[ijet] = Jet_mass_jms_dn[kjet];
			Jet_mass_jms_dn[kjet] = tmppt;
			
			tmppt = Jet_mass_jmr_dn[ijet];
			Jet_mass_jmr_dn[ijet] = Jet_mass_jmr_dn[kjet];
			Jet_mass_jmr_dn[kjet] = tmppt;
			
			
			int tmpid = Jet_electronIdx1[ijet];
			Jet_electronIdx1[ijet] = Jet_electronIdx1[kjet];
			Jet_electronIdx1[kjet] = tmpid;
			
			tmpid = Jet_electronIdx2[ijet];
			Jet_electronIdx2[ijet] = Jet_electronIdx2[kjet];
			Jet_electronIdx2[kjet] = tmpid;
			
			tmpid = Jet_MatchFatJet[ijet];
			Jet_MatchFatJet[ijet] = Jet_MatchFatJet[kjet];
			Jet_MatchFatJet[kjet] = tmpid;
			
			tmpid = Jet_muonIdx1[ijet];
			Jet_muonIdx1[ijet] = Jet_muonIdx1[kjet];
			Jet_muonIdx1[kjet] = tmpid;
			
			tmpid = Jet_muonIdx2[ijet];
			Jet_muonIdx2[ijet] = Jet_muonIdx2[kjet];
			Jet_muonIdx2[kjet] = tmpid;
			
			tmpid = Jet_nConstituents[ijet];
			Jet_nConstituents[ijet] = Jet_nConstituents[kjet];
			Jet_nConstituents[kjet] = tmpid;
			
			tmpid = Jet_nElectrons[ijet];
			Jet_nElectrons[ijet] = Jet_nElectrons[kjet];
			Jet_nElectrons[kjet] = tmpid;
			
			tmpid = Jet_puId[ijet];
			Jet_puId[ijet] = Jet_puId[kjet];
			Jet_puId[kjet] = tmpid;
		 
			if(isMC){

			tmpid = Jet_genJetIdx[ijet];
            Jet_genJetIdx[ijet] = Jet_genJetIdx[kjet];
            Jet_genJetIdx[kjet] = tmpid;

            tmpid = Jet_hadronFlavour[ijet];
            Jet_hadronFlavour[ijet] = Jet_hadronFlavour[kjet];
            Jet_hadronFlavour[kjet] = tmpid;
	
			tmpid = Jet_partonFlavour[ijet];
			Jet_partonFlavour[ijet] = Jet_partonFlavour[kjet];
			Jet_partonFlavour[kjet] = tmpid;
				}
		     }
		   
		}
	}
   
   for(unsigned ijet=0; ijet<nJet; ijet++){
      
	int AK8matchb = -1;
	double minR = 0.4;
	
	for(unsigned kjet=0; kjet<nFatJet; kjet++){
		if(delta2R(FatJet_eta[kjet],FatJet_phi[kjet],Jet_eta[ijet],Jet_phi[ijet])<minR && (Jet_pt[ijet]*1./FatJet_pt[kjet] > 0.2) && (Jet_pt[ijet]*1./FatJet_pt[kjet] < 5)){
			minR = delta2R(FatJet_eta[kjet],FatJet_phi[kjet],Jet_eta[ijet],Jet_phi[ijet]);
			AK8matchb = kjet;
			}
		}
		Jet_MatchFatJet[ijet] = AK8matchb;
    } 
   
   for(unsigned ijet=0; ijet<nFatJet; ijet++){
     
    int AK4matchAK8 = -1;
	double minR = 0.4;
	  
	for(unsigned kjet=0; kjet<nJet; kjet++){
		if(delta2R(FatJet_eta[ijet],FatJet_phi[ijet],Jet_eta[kjet],Jet_phi[kjet])<minR && (Jet_pt[kjet]*1./FatJet_pt[ijet] > 0.2) && (Jet_pt[kjet]*1./FatJet_pt[ijet] < 5)){
			minR = delta2R(FatJet_eta[ijet],FatJet_phi[ijet],Jet_eta[kjet],Jet_phi[kjet]) ;
			AK4matchAK8 = kjet;
			}
	  }
		FatJet_MatchJet[ijet] = AK4matchAK8;
    }
   
   // pt sorting ends for AK4 jet
   
  Event_Ht = 0;    
  for(unsigned ijet=0; ijet<nJet; ijet++){
	 if(Jet_pt[ijet] > 30. && abs(Jet_eta[ijet])<=jeteta_cut){
		Event_Ht += Jet_pt[ijet];
	}
  }
     
     
  int gjet = 0;
  
  for(int ijet=0; ijet<nGenJet; ijet++){
	  
	if(fabs(GenJet_eta[ijet])>5) continue;
	if(GenJet_pt[ijet] < 20) continue;
	
	GenJet_pt[gjet] = GenJet_pt[ijet];
	GenJet_eta[gjet] = GenJet_eta[ijet];
	GenJet_phi[gjet] = GenJet_phi[ijet];
	GenJet_mass[gjet] = GenJet_mass[ijet];
	
	gjet++;
	if(gjet >= njetmax) break;
  }
    
    nGenJet = gjet;
    
  gjet = 0;
  
  for(int ijet=0; ijet<nGenJetAK8; ijet++){
	  
	if(fabs(GenJetAK8_eta[ijet])>5) continue;
	if(GenJetAK8_pt[ijet] < 20) continue;
	
	GenJetAK8_pt[gjet] = GenJetAK8_pt[ijet];
	GenJetAK8_eta[gjet] = GenJetAK8_eta[ijet];
	GenJetAK8_phi[gjet] = GenJetAK8_phi[ijet];
	GenJetAK8_mass[gjet] = GenJetAK8_mass[ijet];
	
	gjet++;
	if(gjet >= njetmax) break;
  }
    
    nGenJetAK8 = gjet;  
   
   gjet = 0;
  
   for(int ijet=0; ijet<nCorrT1METJet; ijet++){
	 
	CorrT1METJet_area[gjet] = CorrT1METJet_area[ijet];
	CorrT1METJet_eta[gjet] = CorrT1METJet_eta[ijet];
	CorrT1METJet_muonSubtrFactor[gjet] = CorrT1METJet_muonSubtrFactor[ijet];
	CorrT1METJet_phi[gjet] = CorrT1METJet_phi[ijet];
	CorrT1METJet_rawPt[gjet] = CorrT1METJet_rawPt[ijet];
	
	gjet++;
	if(gjet >= njetmax) break;
  }
    
    nCorrT1METJet = gjet; 
     
}


int main(int argc, char **argv)
//int main()
{
   isMC = true;
   isQCD = false;
   FakeAnalysis = false;
   isSignal = true;
   isTTBar = false;
   isTTHad = false; isTTSemiLep = false; isTTDiLep = false;
   isST = false;
   
   TopTagging = true;
   if(!isMC){ TopTagging = false; }
   
   char rootfiles[100];

   char outfile[1000];
   char outfilx[1000];
   char infile[1000];
   char datafile[1000];

//   cout <<"Give the input file name"<<endl;
//   cin>> rootfiles;
   int input_mass = atoi(argv[1]);
   sprintf(rootfiles,"Sig_LH_Wp_CompHEP_2018_M%d.log",input_mass);
   
   int len = strlen(rootfiles);
//   strncpy(outfilx, rootfiles, len-4);
   sprintf(outfilx,"Output_Wp_LH_CompHEP_2018_M%c%c%c%c_NanoAOD.root",rootfiles[len-8],rootfiles[len-7],rootfiles[len-6],rootfiles[len-5]);
//   outfilx[len-4]='\0';
//   sprintf(outfilx,"Output_Wp_LH_CompHEP_2018_M%d_NanoAOD.root",input_mass);
   cout<<outfilx<<endl;
//   sprintf (outfile,"%s.root",outfilx);

//   sprintf(outfilx,"Output_Wp_LH_CompHEP_2018_M%d_NanoAOD.root",in_mass);

   TFile *fileout = new TFile(outfilx,"recreate");
   
   Tout = new TTree("Tout", "WeightInfo");

   Tout->Branch("nevent_total", &nevent_total, "nevent_total/I");
   Tout->Branch("weightev", &weightev, "weightev/D");
   Tout->Branch("ncpdf", &ncpdf, "ncpdf/I");
   Tout->Branch("weight_gen_pdf", weight_gen_pdf, "weight_gen_pdf[ncpdf]/D");
   Tout->Branch("ncscale", &ncscale, "ncscale/I");
   Tout->Branch("weight_gen_scale",weight_gen_scale,"weight_gen_scale[ncscale]/D");
   
   Tout1 = new TTree("Tout_Pass", "WeightInfoPass");
   Tout1->Branch("weightpass", &weightpass, "weightpass/D");
   Tout1->Branch("weightpass_btag", &weightpass_btag, "weightpass_btag/D");

   for(int ieta=0; ieta<netarange; ieta++){
	   
	   sprintf(name,"JetpT_AK4_EtaBin%i_BHadron",ieta+1);
	   hist_btAK4_pt[ieta] = new TH1D(name,name,nobptbins,bptbins);
	   hist_btAK4_pt[ieta]->Sumw2();
	   
	   sprintf(name,"BJetpT_AK4_EtaBin%i_BHadron_CSVv2",ieta+1);
	   hist_bAK4_pt[ieta] = new TH1D(name,name,nobptbins,bptbins);
	   hist_bAK4_pt[ieta]->Sumw2();
	   
	   sprintf(name,"BJetpT_AK4_EtaBin%i_BHadron_DeepCSV",ieta+1);
	   hist_dbAK4_pt[ieta] = new TH1D(name,name,nobptbins,bptbins);
	   hist_dbAK4_pt[ieta]->Sumw2();
	   
	   sprintf(name,"BJetpT_AK4_EtaBin%i_BHadron_DeepFlavB",ieta+1);
	   hist_dfbAK4_pt[ieta] = new TH1D(name,name,nobptbins,bptbins);
	   hist_dfbAK4_pt[ieta]->Sumw2();
	   
	   sprintf(name,"JetpT_AK4_EtaBin%i_QHadron",ieta+1);
	   hist_qtAK4_pt[ieta] = new TH1D(name,name,nobptbins,bptbins);
	   hist_qtAK4_pt[ieta]->Sumw2();
	   
	   sprintf(name,"QJetpT_AK4_EtaBin%i_QHadron_CSVv2",ieta+1);
	   hist_qAK4_pt[ieta] = new TH1D(name,name,nobptbins,bptbins);
	   hist_qAK4_pt[ieta]->Sumw2();
	   
	   sprintf(name,"QJetpT_AK4_EtaBin%i_QHadron_DeepCSV",ieta+1);
	   hist_dqAK4_pt[ieta] = new TH1D(name,name,nobptbins,bptbins);
	   hist_dqAK4_pt[ieta]->Sumw2();
	   
	   sprintf(name,"QJetpT_AK4_EtaBin%i_QHadron_DeepFlavB",ieta+1);
	   hist_dfqAK4_pt[ieta] = new TH1D(name,name,nobptbins,bptbins);
	   hist_dfqAK4_pt[ieta]->Sumw2();
	   
	   }
   
   for(int ipt=0; ipt < (noperf_ptbins+1); ipt++){
   
   sprintf(name,"BTagScore_CSVv2_BHadron_ptbin%i",ipt+1);
   if(ipt<noperf_ptbins){	sprintf(title,"CSVv2 Score : p_{T} in [%i,%i] GeV",int(perf_ptbins[ipt]),int(perf_ptbins[ipt+1]));}
   else { sprintf(title,"CSVv2 Score : Inclusive in p_{T}"); }
   hist_btag_csv_bhadron[ipt] = new TH1D(name,title,100,0,1.0);
   hist_btag_csv_bhadron[ipt]->Sumw2();
   
   sprintf(name,"BTagScore_DeepCSV_BHadron_ptbin%i",ipt+1);
   if(ipt<noperf_ptbins){	sprintf(title,"DeepCSV Score : p_{T} in [%i,%i] GeV",int(perf_ptbins[ipt]),int(perf_ptbins[ipt+1]));}
   else { sprintf(title,"DeepCSV Score : Inclusive in p_{T}"); }
   hist_btag_deepcsv_bhadron[ipt] = new TH1D(name,title,100,0,1.0);
   hist_btag_deepcsv_bhadron[ipt]->Sumw2();
   
   sprintf(name,"BTagScore_DeepFlavour_BHadron_ptbin%i",ipt+1);
   if(ipt<noperf_ptbins){	sprintf(title,"DeepFlavour Score : p_{T} in [%i,%i] GeV",int(perf_ptbins[ipt]),int(perf_ptbins[ipt+1]));}
   else { sprintf(title,"DeepFlavour Score : Inclusive in p_{T}"); }
   hist_btag_deepflav_bhadron[ipt] = new TH1D(name,title,100,0,1.0);
   hist_btag_deepflav_bhadron[ipt]->Sumw2();
   
   sprintf(name,"BTagScore_CSVv2_QHadron_ptbin%i",ipt+1);
   if(ipt<noperf_ptbins){	sprintf(title,"CSVv2 Score : p_{T} in [%i,%i] GeV",int(perf_ptbins[ipt]),int(perf_ptbins[ipt+1])); }
   else { sprintf(title,"CSVv2 Score : Inclusive in p_{T}"); }
   hist_btag_csv_qhadron[ipt] = new TH1D(name,title,100,0,1.0);
   hist_btag_csv_qhadron[ipt]->Sumw2();
   
   sprintf(name,"BTagScore_DeepCSV_QHadron_ptbin%i",ipt+1);
   if(ipt<noperf_ptbins){	sprintf(title,"DeepCSV Score : p_{T} in [%i,%i] GeV",int(perf_ptbins[ipt]),int(perf_ptbins[ipt+1]));}
   else { sprintf(title,"DeepCSV Score : Inclusive in p_{T}"); }
   hist_btag_deepcsv_qhadron[ipt] = new TH1D(name,title,100,0,1.0);
   hist_btag_deepcsv_qhadron[ipt]->Sumw2();
   
   sprintf(name,"BTagScore_DeepFlavour_QHadron_ptbin%i",ipt+1);
   if(ipt<noperf_ptbins){	sprintf(title,"DeepFlavour Score : p_{T} in [%i,%i] GeV",int(perf_ptbins[ipt]),int(perf_ptbins[ipt+1]));}
   else { sprintf(title,"DeepFlavour Score : Inclusive in p_{T}"); }
   hist_btag_deepflav_qhadron[ipt] = new TH1D(name,title,100,0,1.0);
   hist_btag_deepflav_qhadron[ipt]->Sumw2();
   
   }
   
   sprintf(name,"Pt_AK8Jet_HasTop");
   sprintf(title,"P_{T} (GeV) of AK8 Jets");
   hist_jetpt_hastop = new TH1D(name,title,noperf_ptbins,perf_ptbins);
   hist_jetpt_hastop->Sumw2();
   
   sprintf(name,"Pt_AK8Jet_HasTop_MDDeepAK8_Pass");
   sprintf(title,"P_{T} (GeV) of AK8 Jets");
   hist_jetpt_hastop_DAK8pass = new TH1D(name,title,noperf_ptbins,perf_ptbins);
   hist_jetpt_hastop_DAK8pass->Sumw2();
   
   sprintf(name,"Pt_AK8Jet_All");
   sprintf(title,"P_{T} (GeV) of AK8 Jets");
   hist_jetpt_all = new TH1D(name,title,noperf_ptbins,perf_ptbins);
   hist_jetpt_all->Sumw2();
   
   sprintf(name,"Pt_AK8Jet_All_MDDeepAK8_Pass");
   sprintf(title,"P_{T} (GeV) of AK8 Jets");
   hist_jetpt_all_DAK8pass = new TH1D(name,title,noperf_ptbins,perf_ptbins);
   hist_jetpt_all_DAK8pass->Sumw2();
   
   sprintf(name,"Pt_AK8Jet_All_TopSDMass");
   sprintf(title,"P_{T} (GeV) of AK8 Jets");
   hist_jetpt_topmsd_all = new TH1D(name,title,noperf_ptbins,perf_ptbins);
   hist_jetpt_topmsd_all->Sumw2();
   
   sprintf(name,"Pt_AK8Jet_All_TopSDMass_MDDeepAK8_Pass");
   sprintf(title,"P_{T} (GeV) of AK8 Jets");
   hist_jetpt_topmsd_all_DAK8pass = new TH1D(name,title,noperf_ptbins,perf_ptbins);
   hist_jetpt_topmsd_all_DAK8pass->Sumw2();
   
   sprintf(name,"N_PV_nopuwt");
   sprintf(title,"# of Primary Vertices");
   hist_npv_nopuwt = new TH1D(name,title,100,-0.1,99.9);//80,-0.1,79.9);
   hist_npv_nopuwt->Sumw2();
   
   sprintf(name,"N_PV");
   sprintf(title,"# of Primary Vertices");
   hist_npv = new TH1D(name,title,100,-0.1,99.9);//80,-0.1,79.9);
   hist_npv->Sumw2();
   
   sprintf(name,"N_PV_Final");
   sprintf(title,"# of Primary Vertices");
   hist_npv_final = new TH1D(name,title,100,-0.1,99.9);//80,-0.1,79.9);
   hist_npv_final->Sumw2();
   
   sprintf(name,"N_PV_Final_nopuwt");
   sprintf(title,"# of Primary Vertices (w/o pu weight)");
   hist_npv_final_nopuwt = new TH1D(name,title,100,-0.1,99.9);//80,-0.1,79.9);
   hist_npv_final_nopuwt->Sumw2();
   
   sprintf(name,"N_PU");
   sprintf(title,"# of Pile-up Vertices");
   hist_npu = new TH1D(name,title,100,0,100);//80,-0.1,79.9);
   hist_npu->Sumw2();
   
   sprintf(name,"N_PU_nopuwt");
   sprintf(title,"# of Pile-up Vertices");
   hist_npu_nopuwt = new TH1D(name,title,100,0,100);//80,-0.1,79.9);
   hist_npu_nopuwt->Sumw2();
   
   sprintf(name,"MET_Filter_Passed");
   sprintf(title,"MET Filter Passing Score");
   hist_metfilter_pass= new TH1D(name,title,2,-0.5,1.5);
   hist_metfilter_pass->Sumw2();
   
   sprintf(name,"Mu_Trigger_Passed");
   sprintf(title,"Muon trigger passing score");
   hist_mutrig_pass= new TH1D(name,title,2,-0.5,1.5);
   hist_mutrig_pass->Sumw2();
   
   sprintf(name,"Elec_Trigger_Passed");
   sprintf(title,"Electron trigger passing score");
   hist_etrig_pass= new TH1D(name,title,2,-0.5,1.5);
   hist_etrig_pass->Sumw2();
   
   sprintf(name,"Pho_Trigger_Passed");
   sprintf(title,"HLT_AK8CHS350_TrimMass30 trigger passing score");
   hist_photrig_pass= new TH1D(name,title,2,-0.5,1.5);
   hist_photrig_pass->Sumw2();
   
   if(isMC){
   
   sprintf(name,"PT_GEN_AK8");
   sprintf(title,"P_{T} of AK8 GEN jets");
   hist_pt_GEN_ak8 = new TH1D(name,title,noptbins,ptbins);
   hist_pt_GEN_ak8->Sumw2();
   
   sprintf(name,"PT_GEN_AK8_Lead");
   sprintf(title,"P_{T} of leading AK8 GEN jet");
   hist_pt_GEN_ak8_lead = new TH1D(name,title,noptbins,ptbins);
   hist_pt_GEN_ak8_lead->Sumw2();
   
   sprintf(name,"PT_GEN_AK4");
   sprintf(title,"P_{T} of AK4 GEN jets");
   hist_pt_GEN_ak4 = new TH1D(name,title,noptbins,ptbins);
   hist_pt_GEN_ak4->Sumw2();
   
   sprintf(name,"PT_GEN_AK4_Lead");
   sprintf(title,"P_{T} of leading AK4 GEN jet");
   hist_pt_GEN_ak4_lead = new TH1D(name,title,noptbins,ptbins);
   hist_pt_GEN_ak4_lead->Sumw2();
   
   sprintf(name,"Eta_GEN_AK8");
   sprintf(title,"Rapidity of AK8 GEN jets");
   hist_eta_GEN_ak8 = new TH1D(name,title,50,-5.,5.);
   hist_eta_GEN_ak8->Sumw2();
   
   sprintf(name,"Eta_GEN_AK8_Lead");
   sprintf(title,"Rapidity of leading AK8 GEN jet");
   hist_eta_GEN_ak8_lead = new TH1D(name,title,50,-5.,5.);
   hist_eta_GEN_ak8_lead->Sumw2();
   
   sprintf(name,"Eta_GEN_AK4");
   sprintf(title,"Rapidity of AK4 GEN jets");
   hist_eta_GEN_ak4 = new TH1D(name,title,50,-5.,5.);
   hist_eta_GEN_ak4->Sumw2();
   
   sprintf(name,"Eta_GEN_AK4_Lead");
   sprintf(title,"Rapidity of leading AK4 GEN jet");
   hist_eta_GEN_ak4_lead = new TH1D(name,title,50,-5.,5.);
   hist_eta_GEN_ak4_lead->Sumw2();
   
   sprintf(name,"M_tb_GEN");
   sprintf(title,"M_{tb} spectra at GEN level");
   hist_mtb_GEN = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_mtb_GEN->Sumw2();
   
	}
   
   sprintf(name,"Had_Trigger_Passed");
   sprintf(title,"Hadronic trigger passing score");
   hist_hadtrig_pass= new TH1D(name,title,2,-0.5,1.5);
   hist_hadtrig_pass->Sumw2();
   
   sprintf(name,"N_Muons");
   sprintf(title,"# of muons");
   hist_nmuons= new TH1D(name,title,5,-0.5,4.5);
   hist_nmuons->Sumw2();
   
   sprintf(name,"N_Electrons");
   sprintf(title,"# of electrons");
   hist_nelectrons= new TH1D(name,title,5,-0.5,4.5);
   hist_nelectrons->Sumw2();
   
   sprintf(name,"N_Photons");
   sprintf(title,"# of photons");
   hist_nphotons= new TH1D(name,title,5,-0.5,4.5);
   hist_nphotons->Sumw2();
   
   sprintf(name,"N_Muons_trigpass");
   sprintf(title,"# of muons (passing hadronic trigger)");
   hist_nmuons_trig_pass= new TH1D(name,title,5,-0.5,4.5);
   hist_nmuons_trig_pass->Sumw2();
   
   sprintf(name,"MET_corpt");
   sprintf(title,"Corrected MET p_{T}");
   met_corpt = new TH1D(name,title,100,0,200);
   met_corpt->Sumw2();
   
   sprintf(name,"MET_corphi");
   sprintf(title,"Corrected MET #phi");
   met_corphi = new TH1D(name,title,50,-M_PI,M_PI);
   met_corphi->Sumw2();
   
   sprintf(name,"MET_corpt_by_sumEt");
   sprintf(title,"Corrected MET P_{T} / SumET");
   met_bysumEt = new TH1D(name,title,200,-0.001,0.999);
   met_bysumEt->Sumw2();
   
   // top tagging histograms //
   
   sprintf(name,"Mass_Hadronic_TopParton");
   sprintf(title,"M (GeV) of top parton");
   hist_mass_hadtopq = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_mass_hadtopq->Sumw2();
   
   sprintf(name,"PID_HardParton_TopDaughter");
   sprintf(title,"PID of hard partons top daughter");
   hist_pid_hadtopq  = new TH1D(name,title,24,-0.1,23.9);
   hist_pid_hadtopq ->Sumw2();
   
   sprintf(name,"Match_AK8Jet_TopParton");
   sprintf(title,"Match AK8Jet of Top Parton");
   hist_matchjet_hadtopq = new TH1D(name,title,20,-1.1,19.9);
   hist_matchjet_hadtopq->Sumw2();
   
   for(int ipt=0; ipt<noperf_ptbins; ipt++){
   
   sprintf(name,"SDMass_TopMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"M_{SD} (GeV) of Top matched AK8 jet");
   hist_jetsdmass_top[ipt] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetsdmass_top[ipt]->Sumw2();
   
   sprintf(name,"Tau32_TopMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"#tau_{3} / #tau_{2} of Top matched AK8 jet");
   hist_jettau32_top[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jettau32_top[ipt]->Sumw2();
   
   sprintf(name,"Subjetbtag_CSVv2_TopMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"Subjet B-tag Value (CSVv2) of Top matched AK8 jet");
   hist_jetsubbtag_CSVv2_top[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetsubbtag_CSVv2_top[ipt]->Sumw2();
   
   sprintf(name,"Subjetbtag_DeepCSV_TopMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"Subjet B-tag Value (DeepCSV) of Top matched AK8 jet");
   hist_jetsubbtag_DeepCSV_top[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetsubbtag_DeepCSV_top[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_MD_TopMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of Top matched AK8 jet");
   hist_jetDeepAK8_MD_top[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetDeepAK8_MD_top[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_MD_withMasscut_TopMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of Top matched AK8 jet (m_{SD}:(105,220) GeV");
   hist_jetDeepAK8_MD_wm_top[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetDeepAK8_MD_wm_top[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_TopMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score (Mass Decorrelated) of Top matched AK8 jet");
   hist_jetDeepAK8_top[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetDeepAK8_top[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_withMasscut_TopMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score (Mass Decorrelated) of Top matched AK8 jet (m_{SD}:(105,220) GeV");
   hist_jetDeepAK8_wm_top[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetDeepAK8_wm_top[ipt]->Sumw2();
   
   sprintf(name,"GenJetMass_TopMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"Gen Jet Mass (GeV) of Top matched AK8 jet");
   hist_jetgenmass_top[ipt] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetgenmass_top[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_TopMatched_AK8Jet_tau32_1_subDeepCSV_1_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of Top matched AK8 jet : #tau_{3} / #tau_{2} < xx && subjet DeepCSV score > yy");
   hist_jetDeepAK8_top_tausub_11[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetDeepAK8_top_tausub_11[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_TopMatched_AK8Jet_tau32_1_subDeepCSV_2_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of Top matched AK8 jet : #tau_{3} / #tau_{2} < xx && subjet DeepCSV score > yy");
   hist_jetDeepAK8_top_tausub_12[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetDeepAK8_top_tausub_12[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_TopMatched_AK8Jet_tau32_2_subDeepCSV_1_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of Top matched AK8 jet : #tau_{3} / #tau_{2} < xx && subjet DeepCSV score > yy");
   hist_jetDeepAK8_top_tausub_21[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetDeepAK8_top_tausub_21[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_TopMatched_AK8Jet_tau32_2_subDeepCSV_2_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of Top matched AK8 jet : #tau_{3} / #tau_{2} < xx && subjet DeepCSV score > yy");
   hist_jetDeepAK8_top_tausub_22[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetDeepAK8_top_tausub_22[ipt]->Sumw2();
   
	} //ipt+1
   /////
   
   sprintf(name,"Mass_HardParton");
   sprintf(title,"M (GeV) of Hard Parton");
   hist_mass_qg = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_mass_qg->Sumw2();
   
   sprintf(name,"PID_HardParton");
   sprintf(title,"PID of Hard Parton");
   hist_pid_qg = new TH1D(name,title,24,-0.1,23.9);
   hist_pid_qg->Sumw2();
   
   for(int ipt=0; ipt<noperf_ptbins; ipt++){
   
   sprintf(name,"SDMass_QMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"M_{SD} (GeV) of q/g matched AK8 jet");
   hist_jetsdmass_tbkg[ipt] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetsdmass_tbkg[ipt]->Sumw2();
   
   sprintf(name,"Tau32_QMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"#tau_{3} / #tau_{2} of q/g matched AK8 jet");
   hist_jettau32_tbkg[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jettau32_tbkg[ipt]->Sumw2();
   
   sprintf(name,"Subjetbtag_CSVv2_QMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"Subjet b tag value (CSVv2) of q/g matched AK8 jet");
   hist_jetsubbtag_CSVv2_tbkg[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetsubbtag_CSVv2_tbkg[ipt]->Sumw2();
   
   sprintf(name,"Subjetbtag_DeepCSV_QMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"Subjet b tag value (DeepCSV) of q/g matched AK8 jet");
   hist_jetsubbtag_DeepCSV_tbkg[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetsubbtag_DeepCSV_tbkg[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_MD_QMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"DeepAK8 score of q/g matched AK8 jet");
   hist_jetDeepAK8_MD_tbkg[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetDeepAK8_MD_tbkg[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_MD_withMasscut_QMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"DeepAK8 score of q/g matched AK8 jet (m_{SD}:(105,220) GeV");
   hist_jetDeepAK8_MD_wm_tbkg[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetDeepAK8_MD_wm_tbkg[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_QMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"DeepAK8 score (mass decorrelated) of q/g matched AK8 jet");
   hist_jetDeepAK8_tbkg[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetDeepAK8_tbkg[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_withMasscut_QMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"DeepAK8 score of q/g matched AK8 jet (m_{SD}:(105,220) GeV");
   hist_jetDeepAK8_wm_tbkg[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetDeepAK8_wm_tbkg[ipt]->Sumw2();
   
   sprintf(name,"GenJetMass_QMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"Gen jet mass (GeV) of q/g matched AK8 jet");
   hist_jetgenmass_tbkg[ipt] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetgenmass_tbkg[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_QMatched_AK8Jet_tau32_1_subDeepCSV_1_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of q/g matched AK8 jet : #tau_{3} / #tau_{2} < xx && subjet DeepCSV score > yy");
   hist_jetDeepAK8_tbkg_tausub_11[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetDeepAK8_tbkg_tausub_11[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_QMatched_AK8Jet_tau32_1_subDeepCSV_2_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of q/g matched AK8 jet : #tau_{3} / #tau_{2} < xx && subjet DeepCSV score > yy");
   hist_jetDeepAK8_tbkg_tausub_12[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetDeepAK8_tbkg_tausub_12[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_QMatched_AK8Jet_tau32_2_subDeepCSV_1_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of q/g matched AK8 jet : #tau_{3} / #tau_{2} < xx && subjet DeepCSV score > yy");
   hist_jetDeepAK8_tbkg_tausub_21[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetDeepAK8_tbkg_tausub_21[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_QMatched_AK8Jet_tau32_2_subDeepCSV_2_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of q/g matched AK8 jet : #tau_{3} / #tau_{2} < xx && subjet DeepCSV score > yy");
   hist_jetDeepAK8_tbkg_tausub_22[ipt] = new TH1D(name,title,50,0,1.0);
   hist_jetDeepAK8_tbkg_tausub_22[ipt]->Sumw2();
   
	}
   
   // end of top tagging histograms //
   
   sprintf(name,"NJets_AK8CHS");
   sprintf(title,"# of AK8CHS jets");
   hist_njetAK8 = new TH1D(name,title,10,-0.1,9.9);
   hist_njetAK8->Sumw2();
   
   sprintf(name,"Jet_Index_TopCand_AK8CHS");
   sprintf(title,"AK8 Top Candidate Jet Index");
   hist_topcand_AK8 = new TH1D(name,title,5,-0.1,4.9);
   hist_topcand_AK8->Sumw2();
   
   sprintf(name,"Pt_AK8CHS");
   sprintf(title,"P_{T} of AK8CHS jets");
   hist_jetptAK8 = new TH1D(name,title,noptbins,ptbins);
   hist_jetptAK8->Sumw2();
   
   sprintf(name,"Mass_AK8CHS");
   sprintf(title,"Mass of AK8CHS jets");
   hist_jetmassAK8 = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetmassAK8->Sumw2();
   
   sprintf(name,"SDMass_AK8CHS");
   sprintf(title,"Soft drop mass of AK8CHS jets");
   hist_jetsdmassAK8 = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetsdmassAK8->Sumw2();
   
   sprintf(name,"y_AK8CHS");
   sprintf(title,"Rapidity of AK8CHS jets");
   hist_jetrapAK8 = new TH1D(name,title,50,-5.,5.);
   hist_jetrapAK8->Sumw2();
  
   sprintf(name,"rho_AK8CHS");
   sprintf(title,"#rho (-2*log(m_{SD}/p_{T})) of AK8CHS jets");
   hist_jetrhoAK8 = new TH1D(name,title,norhobins,logrhobins);
   hist_jetrhoAK8->Sumw2();
   
   sprintf(name,"tau32_AK8CHS");
   sprintf(title,"#tau_3 / #tau_2 of AK8CHS jets");
   hist_jettau32AK8 = new TH1D(name,title,50,0,1.0);
   hist_jettau32AK8->Sumw2();
   
   sprintf(name,"tau21_AK8CHS");
   sprintf(title,"#tau_2 / #tau_1 of AK8CHS jets");
   hist_jettau21AK8 = new TH1D(name,title,50,-0.1,0.99);
   hist_jettau21AK8->Sumw2();
   
   sprintf(name,"N_Subjets_passing_medium_WP");
   sprintf(title,"# of subjets Passing Medium WP");
   hist_njetsubjets_AK8_mpass = new TH1D(name,title,3,-0.1,2.9);
   hist_njetsubjets_AK8_mpass->Sumw2();
   
   sprintf(name,"N_Subjets_passing_loose_WP_fail_medium_WP");
   sprintf(title,"# of subjets Passing Loose WP && Failing Medium WP");
   hist_njetsubjets_AK8_lpass_mfail = new TH1D(name,title,3,-0.1,2.9);
   hist_njetsubjets_AK8_lpass_mfail->Sumw2();
   
   sprintf(name,"subjetbtag_AK8CHS");
   sprintf(title,"Suubjet CSVv2 Score of AK8CHS jets");
   hist_jetsubbtagAK8 = new TH1D(name,title,20,0,1.0);
   hist_jetsubbtagAK8->Sumw2();
   
   sprintf(name,"subjetmass_AK8CHS");
   sprintf(title,"Suubjet Mass of AK8CHS jets");
   hist_jetsubmassAK8 = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetsubmassAK8->Sumw2();
   
   sprintf(name,"subjetcostheta_AK8CHS");
   sprintf(title,"cos(#Delta#theta^{tb_{t}})");
   hist_jetsubthetAK8 = new TH1D(name,title,50,-1.,+1.);
   hist_jetsubthetAK8->Sumw2();
   
   sprintf(name,"SDMass_AK8CHS_TopCand");
   sprintf(title,"Soft drop mass of AK8CHS top candidate jet");
   hist_topjetsdmass = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass->Sumw2();
   
   sprintf(name,"Tau32_AK8CHS_TopCand");
   sprintf(title,"#tau_{32} of AK8CHS top candidate jet");
   hist_topjettau32AK8 = new TH1D(name,title,50,0,1.0);
   hist_topjettau32AK8->Sumw2();
   
   sprintf(name,"DeepTagTvsQCD_AK8CHS_TopCand");
   sprintf(title,"DeepTagTvsQCD of AK8CHS top candidate jet");
   hist_topjetdeeptopscore = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeeptopscore->Sumw2();
   
   sprintf(name,"DeepTagTvsQCD_MD_AK8CHS_TopCand");
   sprintf(title,"DeepTagTvsQCD (Mass Decorrelated) Score of AK8CHS top candidate jet");
   hist_topjetmddeeptopscore = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetmddeeptopscore->Sumw2();
   
   sprintf(name,"Pt_AK8CHS_TopCand");
   sprintf(title,"P_{T} of AK8CHS jets");
   hist_topjetpt= new TH1D(name,title,noptbins,ptbins);
   hist_topjetpt->Sumw2();
   
   sprintf(name,"Pt_Msd_AK8CHS_TopCand");
   sprintf(title,"Correlation between m_{SD} and p_{T} of top candidate jet");
   hist_topjetptmsd_2d = new TH2D(name,title,nsdmassbins,sdmassbins,noptbins,ptbins);
   hist_topjetptmsd_2d->Sumw2();
   
   sprintf(name,"Pt_AK8CHS_TopCand_PartonMatched");
   sprintf(title,"P_{T} of AK8CHS Jets (matched to t parton)");
   hist_topjetptmatch= new TH1D(name,title,noptbins,ptbins);
   hist_topjetptmatch->Sumw2();
   
   sprintf(name,"Pt_AK8CHS_Top_Parton");
   sprintf(title,"P_{T} of top parton)");
   hist_toppartonpt= new TH1D(name,title,noptbins,ptbins);
   hist_toppartonpt->Sumw2();
   
   sprintf(name,"Pt_AK8CHS_TopCand_Matched_to_Parton");
   sprintf(title,"P_{T} of AK8CHS Jets (matched to t parton)");
   hist_topparton_matchedtopjet_pt= new TH1D(name,title,noptbins,ptbins);
   hist_topparton_matchedtopjet_pt->Sumw2();
   
   sprintf(name,"GenMass_AK8CHS_TopCand");
   sprintf(title,"Mass of GEN jet corresponding to AK8CHS top candidate jet");
   hist_topjetgenmass = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetgenmass->Sumw2();
   
   sprintf(name,"GenPID_AK8CHS_TopCand");
   sprintf(title,"Parton Flavour of GEN jet corresponding to AK8CHS top candidate jet");
   hist_topjetgenpid = new TH1D(name,title,22,-0.1,21.9);
   hist_topjetgenpid->Sumw2();
   
   sprintf(name,"H2D_SDMass_Tau32_AK8CHS");
   sprintf(title,"2D Correlation of m_{SD} and #tau_{32} of AK8CHS jets");
   hist_2D_sdmass_tau32_AK8 = new TH2D(name,title,nsdmassbins,sdmassbins,50,float(0),float(1.0));
   hist_2D_sdmass_tau32_AK8->Sumw2();
   
   sprintf(name,"H2D_SDMass_subbtag_AK8CHS");
   sprintf(title,"2D Correlation of m_{SD} and subjet b-tag AK8CHS jets");
   hist_2D_sdmass_subbtag_AK8 = new TH2D(name,title,nsdmassbins,sdmassbins,50,float(0),float(1.0));
   hist_2D_sdmass_subbtag_AK8->Sumw2();
   
   sprintf(name,"H2D_Tau32_subbtag_AK8CHS");
   sprintf(title,"2D Correlation of #tau_{32} and subjet b-tag of AK8CHS jets");
   hist_2D_tau32_subbtag_AK8 = new TH2D(name,title,50,0,1.0,20,0,1.0);
   hist_2D_tau32_subbtag_AK8->Sumw2();
   
   sprintf(name,"TopCand_SubJetMass_12");
   sprintf(title,"2D Correlation of subjet Mass of AK8 jet t-Cand");
   hist_2D_topjet_subjetmass12 = new TH2D(name,title,nsdmassbins,sdmassbins,nsdmassbins,sdmassbins);
   hist_2D_topjet_subjetmass12->Sumw2();
   
   sprintf(name,"TopCand_SubJetMass1_vs_JetSDMass");
   sprintf(title,"2D Correlation of higher subjet Mass and SDMass AK8 jet t-Cand");
   hist_2D_topjetsdmass_subjetmass1 = new TH2D(name,title,nsdmassbins,sdmassbins,nsdmassbins,sdmassbins);
   hist_2D_topjetsdmass_subjetmass1->Sumw2();
   
   sprintf(name,"TopCand_SubJetMass2_vs_JetSDMass");
   sprintf(title,"2D Correlation of Lower subjet Mass and SDMass AK8 jet t-Cand");
   hist_2D_topjetsdmass_subjetmass2 = new TH2D(name,title,nsdmassbins,sdmassbins,nsdmassbins,sdmassbins);
   hist_2D_topjetsdmass_subjetmass2->Sumw2();
   
   // m-pass //
   
   sprintf(name,"tau32_mpass_AK8CHS");
   sprintf(title,"#tau_3 / #tau_2 of AK8CHS Jets within (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_jettau32AK8_mpass = new TH1D(name,title,50,0,1.0);
   hist_jettau32AK8_mpass->Sumw2();
   
   sprintf(name,"N_Subjets_passing_medium_WP_mpass");
   sprintf(title,"# of subjets Passing Medium WP within (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_njetsubjets_AK8_BMPpass_mpass = new TH1D(name,title,3,-0.1,2.9);
   hist_njetsubjets_AK8_BMPpass_mpass->Sumw2();
   
   sprintf(name,"NBtagSubJets_mpass_AK8CHS");
   sprintf(title,"# of Btagged subjets AK8CHS Jets within (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_nbsubjetAK8_mpass = new TH1D(name,title,3,-0.1,2.9);
   hist_nbsubjetAK8_mpass->Sumw2();
   
   sprintf(name,"N_Subjets_passing_loose_WP_fail_medium_WP_mpass");
   sprintf(title,"# of subjets Passing Loose WP && Failing Medium WP within (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_njetsubjets_AK8_BLPpass_BMPfail_mpass = new TH1D(name,title,3,-0.1,2.9);
   hist_njetsubjets_AK8_BLPpass_BMPfail_mpass->Sumw2();
   
   sprintf(name,"subjetmass_mpass_AK8CHS");
   sprintf(title,"Suubjet Mass of AK8CHS Jets within (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_jetsubmassAK8_mpass = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetsubmassAK8_mpass->Sumw2();
   
   sprintf(name,"subjetpt_mpass_AK8CHS");
   sprintf(title,"Suubjet P_{T} of AK8CHS Jets within (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_jetsubptAK8_mpass = new TH1D(name,title,noptbins,ptbins);
   hist_jetsubptAK8_mpass->Sumw2();
   
   sprintf(name,"subjetbtag_mpass_AK8CHS");
   sprintf(title,"Suubjet CSVv2 Score of AK8CHS Jets within (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_jetsubbtagAK8_mpass = new TH1D(name,title,20,-0.001,0.999);
   hist_jetsubbtagAK8_mpass->Sumw2();
   
   sprintf(name,"subjetbtag_deepCSV_mpass_AK8CHS");
   sprintf(title,"Suubjet DeepCSV Score of AK8CHS Jets within (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_jetsubbtagAK8_deepCSV_mpass = new TH1D(name,title,20,-0.001,0.999);
   hist_jetsubbtagAK8_deepCSV_mpass->Sumw2();
   
   sprintf(name,"subjetpt_btagpass_mpass_AK8CHS");
   sprintf(title,"B-tagged Suubjet P_{T} of AK8CHS Jets within (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_jetsubptAK8_btagpass_mpass = new TH1D(name,title,noptbins,ptbins);
   hist_jetsubptAK8_btagpass_mpass->Sumw2();
   
   sprintf(name,"H2D_Tau32_subbtag_mpass_AK8CHS");
   sprintf(title,"2D Correlation of #tau_{32} and subjet b-tag of AK8CHS Jets  within (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_2D_tau32_subbtag_AK8_mpass = new TH2D(name,title,50,0,1.0,20,0,1.0);
   hist_2D_tau32_subbtag_AK8_mpass->Sumw2();
   
   // m-fail //
   
   sprintf(name,"tau32_mfail_AK8CHS");
   sprintf(title,"#tau_3 / #tau_2 of AK8CHS Jets outside (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_jettau32AK8_mfail = new TH1D(name,title,50,0,1.0);
   hist_jettau32AK8_mfail->Sumw2();
   
   sprintf(name,"N_Subjets_passing_medium_WP_mfail");
   sprintf(title,"# of subjets Passing Medium WP outside (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_njetsubjets_AK8_BMPpass_mfail = new TH1D(name,title,3,-0.1,2.9);
   hist_njetsubjets_AK8_BMPpass_mfail->Sumw2();
   
   sprintf(name,"NBtagSubJets_mfail_AK8CHS");
   sprintf(title,"# of Btagged subjets AK8CHS Jets outside (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_nbsubjetAK8_mfail = new TH1D(name,title,3,-0.1,2.9);
   hist_nbsubjetAK8_mfail->Sumw2();
   
   sprintf(name,"N_Subjets_passing_loose_WP_fail_medium_WP_mfail");
   sprintf(title,"# of subjets Passing Loose WP && Failing Medium WP outside (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_njetsubjets_AK8_BLPpass_BMPfail_mfail = new TH1D(name,title,3,-0.1,2.9);
   hist_njetsubjets_AK8_BLPpass_BMPfail_mfail->Sumw2();
   
   sprintf(name,"subjetmass_mfail_AK8CHS");
   sprintf(title,"Suubjet Mass of AK8CHS Jets outside (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_jetsubmassAK8_mfail = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetsubmassAK8_mfail->Sumw2();
   
   sprintf(name,"subjetpt_mfail_AK8CHS");
   sprintf(title,"Suubjet P_{T} of AK8CHS Jets outside (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_jetsubptAK8_mfail = new TH1D(name,title,noptbins,ptbins);
   hist_jetsubptAK8_mfail->Sumw2();
   
   sprintf(name,"subjetbtag_mfail_AK8CHS");
   sprintf(title,"Suubjet CSVv2 Score of AK8CHS Jets outside (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_jetsubbtagAK8_mfail = new TH1D(name,title,20,0,1.0);
   hist_jetsubbtagAK8_mfail->Sumw2();
   
   sprintf(name,"subjetbtag_deepCSV_mfail_AK8CHS");
   sprintf(title,"Suubjet DeepCSV Score of AK8CHS Jets outside (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_jetsubbtagAK8_deepCSV_mfail = new TH1D(name,title,20,0,1.0);
   hist_jetsubbtagAK8_deepCSV_mfail->Sumw2();
   
   sprintf(name,"subjetpt_btagpass_mfail_AK8CHS");
   sprintf(title,"B-tagged Suubjet P_{T} of AK8CHS Jets outside (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_jetsubptAK8_btagpass_mfail = new TH1D(name,title,noptbins,ptbins);
   hist_jetsubptAK8_btagpass_mfail->Sumw2();
   
   sprintf(name,"H2D_Tau32_subbtag_mfail_AK8CHS");
   sprintf(title,"2D Correlation of #tau_{32} and subjet b-tag of AK8CHS Jets  outside (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_2D_tau32_subbtag_AK8_mfail = new TH2D(name,title,50,0,1.0,20,0,1.0);
   hist_2D_tau32_subbtag_AK8_mfail->Sumw2();
   
   sprintf(name,"H2D_SDMass_TopvsQCD_AK8CHS");
   sprintf(title,"2D Correlation of M_{SD} and DeepAK8 Toptag Score");
   hist_2D_sdmass_deeptopscore = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.001,0.999);
   hist_2D_sdmass_deeptopscore->Sumw2();
   
   sprintf(name,"H2D_Mass_TopvsQCD_AK8CHS");
   sprintf(title,"2D Correlation of Mass and DeepAK8 Toptag Score");
   hist_2D_mass_deeptopscore = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.001,0.999);
   hist_2D_mass_deeptopscore->Sumw2();
   
   sprintf(name,"H2D_Tau32_TopvsQCD_AK8CHS");
   sprintf(title,"2D Correlation of #tau_{32} and DeepAK8 Toptag Score");
   hist_2D_tau32_deeptopscore = new TH2D(name,title,50,-0.001,0.999,50,-0.001,0.999);
   hist_2D_tau32_deeptopscore->Sumw2();
   
   sprintf(name,"H2D_Subjetbtag_TopvsQCD_AK8CHS");
   sprintf(title,"2D Correlation of subjet b-tag (DeepCSV) and DeepAK8 Toptag Score");
   hist_2D_subjetbtag_deeptopscore = new TH2D(name,title,50,-0.001,0.999,50,-0.001,0.999);
   hist_2D_subjetbtag_deeptopscore->Sumw2();
   
   sprintf(name,"H2D_SDMass_MD_TopvsQCD_AK8CHS");
   sprintf(title,"2D Correlation of M_{SD} and Mass Decorrelated DeepAK8 Toptag Score");
   hist_2D_sdmass_deepmdtopscore = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.001,0.999);
   hist_2D_sdmass_deepmdtopscore->Sumw2();
   
   sprintf(name,"H2D_Mass_MD_TopvsQCD_AK8CHS");
   sprintf(title,"2D Correlation of Mass and Mass Decorrelated DeepAK8 Toptag Score");
   hist_2D_mass_deepmdtopscore = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.001,0.999);
   hist_2D_mass_deepmdtopscore->Sumw2();
   
   sprintf(name,"H2D_Tau32_MD_TopvsQCD_AK8CHS");
   sprintf(title,"2D Correlation of #tau_{32} and Mass Decorrelated DeepAK8 Toptag Score");
   hist_2D_tau32_deepmdtopscore = new TH2D(name,title,50,-0.001,0.999,50,-0.001,0.999);
   hist_2D_tau32_deepmdtopscore->Sumw2();
   
   sprintf(name,"H2D_Subjetbtag_MD_TopvsQCD_AK8CHS");
   sprintf(title,"2D Correlation of subjet b-tag (DeepCSV) and Mass Decorrelated DeepAK8 Toptag Score");
   hist_2D_subjetbtag_deepmdtopscore = new TH2D(name,title,50,-0.001,0.999,50,-0.001,0.999);
   hist_2D_subjetbtag_deepmdtopscore->Sumw2();
   
		for(int ib=0; ib<nbtag; ib++){
		for(int itau=0; itau<ntautag; itau++){
    
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_Tau32Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		if(itau==0) { sprintf(title,"2D Correlation of Soft drop mass  and subjet b-tag (CSVv2) of AK8CHS Jets | #tau_{32} < %0.2f btag %i",tau32_cut,ib); }
		if(itau==1) { sprintf(title,"2D Correlation of Soft drop mass  and subjet b-tag (CSVv2) of AK8CHS Jets | 0.5 < #tau_{32} < %0.2f btag %i",tau32_cut_loose,ib); }
		if(itau==2) { sprintf(title,"2D Correlation of Soft drop mass  and subjet b-tag (CSVv2) of AK8CHS Jets | #tau_{32} > %0.2f btag %i",tau32_cut_loose,ib); }
		
		hist_2D_sdmass_subbtag_CSVv2_AK8[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_subbtag_CSVv2_AK8[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_Tau32Tag%i_AK4btag%i_AK8CHS_lowbpt",itau+1,ib+1);
		hist_2D_sdmass_subbtag_CSVv2_AK8_1[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_subbtag_CSVv2_AK8_1[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_Tau32Tag%i_AK4btag%i_AK8CHS_nobmasscut",itau+1,ib+1);
		hist_2D_sdmass_subbtag_CSVv2_AK8_2[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_subbtag_CSVv2_AK8_2[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_Tau32Tag%i_AK4btag%i_AK8CHS_highht",itau+1,ib+1);
		hist_2D_sdmass_subbtag_CSVv2_AK8_3[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_subbtag_CSVv2_AK8_3[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_DeepCSV_Tau32Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		if(itau==0) { sprintf(title,"2D Correlation of Soft drop mass  and subjet b-tag (DeepCSV) of AK8CHS Jets | #tau_{32} < %0.2f btag %i",tau32_cut,ib); }
		if(itau==1) { sprintf(title,"2D Correlation of Soft drop mass  and subjet b-tag (DeepCSV) of AK8CHS Jets | 0.5 < #tau_{32} < %0.2f btag %i",tau32_cut_loose,ib); }
		if(itau==2) { sprintf(title,"2D Correlation of Soft drop mass  and subjet b-tag (DeepCSV) of AK8CHS Jets | #tau_{32} > %0.2f btag %i",tau32_cut_loose,ib); }
		
		hist_2D_sdmass_subbtag_DeepCSV_AK8[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_subbtag_DeepCSV_AK8[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_DeepCSV_Tau32Tag%i_AK4btag%i_AK8CHS_lowbpt",itau+1,ib+1);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_1[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_1[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_DeepCSV_Tau32Tag%i_AK4btag%i_AK8CHS_nobmasscut",itau+1,ib+1);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_2[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_2[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_DeepCSV_Tau32Tag%i_AK4btag%i_AK8CHS_highht",itau+1,ib+1);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_3[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_3[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_DeepAK8Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		hist_2D_sdmass_subbtag_CSVv2_AK8_DeepAK8[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_subbtag_CSVv2_AK8_DeepAK8[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_DeepCSV_DeepAK8Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_DeepAK8[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_DeepAK8[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_Tau32_DeepAK8Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		hist_2D_sdmass_deeptag_AK8_reg_DeepAK8[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_deeptag_AK8_reg_DeepAK8[ib][itau]->Sumw2();
		
		
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_MDDeepAK8Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_MDDeepAK8Tag%i_AK4btag%i_AK8CHS_hight",itau+1,ib+1);
		hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8_3[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8_3[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_DeepCSV_MDDeepAK8Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_DeepCSV_MDDeepAK8Tag%i_AK4btag%i_AK8CHS_hight",itau+1,ib+1);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8_3[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8_3[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_Tau32_MDDeepAK8Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		hist_2D_sdmass_deeptag_AK8_reg_md_DeepAK8[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_deeptag_AK8_reg_md_DeepAK8[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_Tau32_MDDeepAK8Tag%i_AK4btag%i_AK8CHS_hight",itau+1,ib+1);
		hist_2D_sdmass_deeptag_AK8_reg_md_DeepAK8_3[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_deeptag_AK8_reg_md_DeepAK8_3[ib][itau]->Sumw2();
		
		// t-tbar (MD) //
		
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_DeepAK8Tag%i_AK4btag%i_AK8CHS_ttbar",itau+1,ib+1);
		hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8_tt[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8_tt[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_DeepAK8Tag%i_AK4btag%i_AK8CHS_ttbar_hight",itau+1,ib+1);
		hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8_tt_3[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0.0,1.0);
		hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8_tt_3[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_DeepCSV_MDDeepAK8Tag%i_AK4btag%i_AK8CHS_ttbar",itau+1,ib+1);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8_tt[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8_tt[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_DeepCSV_MDDeepAK8Tag%i_AK4btag%i_AK8CHS_ttbar_hight",itau+1,ib+1);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8_tt_3[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8_tt_3[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_Tau32_MDDeepAK8Tag%i_AK4btag%i_AK8CHS_ttbar",itau+1,ib+1);
		hist_2D_sdmass_deeptag_AK8_reg_md_DeepAK8_tt[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_deeptag_AK8_reg_md_DeepAK8_tt[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_Tau32_MDDeepAK8Tag%i_AK4btag%i_AK8CHS_ttbar_hight",itau+1,ib+1);
		hist_2D_sdmass_deeptag_AK8_reg_md_DeepAK8_tt_3[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_deeptag_AK8_reg_md_DeepAK8_tt_3[ib][itau]->Sumw2();
		
		// t-tbar (MD) ends
		
		sprintf(name,"H2D_SDMass_tbMass_TauTag%i_AK4btag%i_AK8CHS_lowbpt",itau+1,ib+1);
		hist_2D_sdmass_tbmass_1[ib][itau] = new TH2D(name,name,nsdmassbins,sdmassbins,nvwpmbin,wpmbins);
		hist_2D_sdmass_tbmass_1[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_tbMass_TauTag%i_AK4btag%i_AK8CHS_nobmasscut",itau+1,ib+1);
		hist_2D_sdmass_tbmass_2[ib][itau] = new TH2D(name,name,nsdmassbins,sdmassbins,nvwpmbin,wpmbins);
		hist_2D_sdmass_tbmass_2[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_tbMass_TauTag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		hist_2D_sdmass_tbmass[ib][itau] = new TH2D(name,name,nsdmassbins,sdmassbins,nvwpmbin,wpmbins);
		hist_2D_sdmass_tbmass[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_tbMass_TauTag%i_AK4btag%i_AK8CHS_highht",itau+1,ib+1);
		hist_2D_sdmass_tbmass_3[ib][itau] = new TH2D(name,name,nsdmassbins,sdmassbins,nvwpmbin,wpmbins);
		hist_2D_sdmass_tbmass_3[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_tbMass_DeepAK8Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		hist_2D_sdmass_tbmass_deepAK8[ib][itau] = new TH2D(name,name,nsdmassbins,sdmassbins,nvwpmbin,wpmbins);
		hist_2D_sdmass_tbmass_deepAK8[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_tbMass_MDDeepAK8Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		hist_2D_sdmass_tbmass_mddeepAK8[ib][itau] = new TH2D(name,name,nsdmassbins,sdmassbins,nvwpmbin,wpmbins);
		hist_2D_sdmass_tbmass_mddeepAK8[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_tbMass_MDDeepAK8Tag%i_AK4btag%i_AK8CHS_highht",itau+1,ib+1);
		hist_2D_sdmass_tbmass_mddeepAK8_3[ib][itau] = new TH2D(name,name,nsdmassbins,sdmassbins,nvwpmbin,wpmbins);
		hist_2D_sdmass_tbmass_mddeepAK8_3[ib][itau]->Sumw2();
		
		
		// t-tbar
		
		sprintf(name,"H2D_SDMass_tbMass_MDDeepAK8Tag%i_AK4btag%i_AK8CHS_ttbar",itau+1,ib+1);
		hist_2D_sdmass_tbmass_mddeepAK8_tt[ib][itau] = new TH2D(name,name,nsdmassbins,sdmassbins,nvwpmbin,wpmbins);
		hist_2D_sdmass_tbmass_mddeepAK8_tt[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_tbMass_MDDeepAK8Tag%i_AK4btag%i_AK8CHS_ttbar_highht",itau+1,ib+1);
		hist_2D_sdmass_tbmass_mddeepAK8_tt_3[ib][itau] = new TH2D(name,name,nsdmassbins,sdmassbins,nvwpmbin,wpmbins);
		hist_2D_sdmass_tbmass_mddeepAK8_tt_3[ib][itau]->Sumw2();
		
		//t-tbar ends
		
			}
		
		sprintf(name,"H2D_DeepAK8_tbMass_AK4btag%i_AK8CHS_lowbpt",ib+1);
		hist_2D_DeepAK8_tbmass_1[ib] = new TH2D(name,name,20,0.0,1.0,nvwpmbin,wpmbins);
		hist_2D_DeepAK8_tbmass_1[ib]->Sumw2();
		
		sprintf(name,"H2D_DeepAK8_tbMass_AK4btag%i_AK8CHS_nobmasscut",ib+1);
		hist_2D_DeepAK8_tbmass_2[ib] = new TH2D(name,name,20,0.0,1.0,nvwpmbin,wpmbins);
		hist_2D_DeepAK8_tbmass_2[ib]->Sumw2();
		
		sprintf(name,"H2D_DeepAK8_tbMass_AK4btag%i_AK8CHS",ib+1);
		hist_2D_DeepAK8_tbmass[ib] = new TH2D(name,name,20,0.0,1.0,nvwpmbin,wpmbins);
		hist_2D_DeepAK8_tbmass[ib]->Sumw2();
		
		sprintf(name,"H2D_DeepAK8_tbMass_AK4btag%i_AK8CHS_highht",ib+1);
		hist_2D_DeepAK8_tbmass_3[ib] = new TH2D(name,name,20,0.0,1.0,nvwpmbin,wpmbins);
		hist_2D_DeepAK8_tbmass_3[ib]->Sumw2();
		
		sprintf(name,"H2D_MD_DeepAK8_tbMass_AK4btag%i_AK8CHS_lowbpt",ib+1);
		hist_2D_MDDeepAK8_tbmass_1[ib] = new TH2D(name,name,20,0.0,1.0,nvwpmbin,wpmbins);
		hist_2D_MDDeepAK8_tbmass_1[ib]->Sumw2();
		
		sprintf(name,"H2D_MD_DeepAK8_tbMass_AK4btag%i_AK8CHS_nobmasscut",ib+1);
		hist_2D_MDDeepAK8_tbmass_2[ib] = new TH2D(name,name,20,0.0,1.0,nvwpmbin,wpmbins);
		hist_2D_MDDeepAK8_tbmass_2[ib]->Sumw2();
		
		sprintf(name,"H2D_MD_DeepAK8_tbMass_AK4btag%i_AK8CHS",ib+1);
		hist_2D_MDDeepAK8_tbmass[ib] = new TH2D(name,name,20,0.0,1.0,nvwpmbin,wpmbins);
		hist_2D_MDDeepAK8_tbmass[ib]->Sumw2();
		
		sprintf(name,"H2D_MD_DeepAK8_tbMass_AK4btag%i_AK8CHS_highht",ib+1);
		hist_2D_MDDeepAK8_tbmass_3[ib] = new TH2D(name,name,20,0.0,1.0,nvwpmbin,wpmbins);
		hist_2D_MDDeepAK8_tbmass_3[ib]->Sumw2();
			
		for(int isb=0; isb<nsbtag; isb++){
		
		sprintf(name,"H2D_SDMass_tau32_subjetbtag%i_AK4btag%i_AK8CHS",isb+1,ib+1);
		sprintf(title,"2D Correlation of Soft drop mass  and  #tau_{32} | subjet b-tag %i AK4 b-tag %i",isb,ib);
		
		hist_2D_sdmass_tau32_AK8_reg[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
		hist_2D_sdmass_tau32_AK8_reg[ib][isb]->Sumw2();
		
		sprintf(name,"H2D_SDMass_tau32_subjetbtag%i_AK4btag%i_AK8CHS_lowbpt",isb+1,ib+1);
		
		hist_2D_sdmass_tau32_AK8_reg_1[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
		hist_2D_sdmass_tau32_AK8_reg_1[ib][isb]->Sumw2();
			
		sprintf(name,"H2D_SDMass_tau32_subjetbtag%i_AK4btag%i_AK8CHS_nobmasscut",isb+1,ib+1);
		
		hist_2D_sdmass_tau32_AK8_reg_2[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
		hist_2D_sdmass_tau32_AK8_reg_2[ib][isb]->Sumw2();	
			
		sprintf(name,"H2D_SDMass_tau32_subjetbtag%i_AK4btag%i_AK8CHS_highht",isb+1,ib+1);
		
		hist_2D_sdmass_tau32_AK8_reg_3[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
		hist_2D_sdmass_tau32_AK8_reg_3[ib][isb]->Sumw2();	
			
		sprintf(name,"H2D_SDMass_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS",isb+1,ib+1);	
	    sprintf(title,"H2D SDMass DeepAK8 subjetbtag%i AK4btag%i AK8CHS",isb+1,ib+1);	
	    
	    hist_2D_sdmass_deeptopscore_AK8[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_AK8[ib][isb]->Sumw2();
		
	    sprintf(name,"H2D_SDMass_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_lowbpt",isb+1,ib+1);	
	    
	    hist_2D_sdmass_deeptopscore_AK8_1[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_AK8_1[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_nobmasscut",isb+1,ib+1);	
	   
	    hist_2D_sdmass_deeptopscore_AK8_2[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_AK8_2[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_highht",isb+1,ib+1);	
	   
	    hist_2D_sdmass_deeptopscore_AK8_3[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_AK8_3[ib][isb]->Sumw2();
	
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS",isb+1,ib+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 SubJet btag%i AK4btag%i AK8CHS",isb+1,ib+1);	
	    	
	    hist_2D_sdmass_deeptopscore_md_AK8[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_md_AK8[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_lowbpt",isb+1,ib+1);
	    	
	    hist_2D_sdmass_deeptopscore_md_AK8_1[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_md_AK8_1[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_nobmasscut",isb+1,ib+1);
	    	
	    hist_2D_sdmass_deeptopscore_md_AK8_2[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_md_AK8_2[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_highht",isb+1,ib+1);
	    	
	    hist_2D_sdmass_deeptopscore_md_AK8_3[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_md_AK8_3[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD",isb+1,ib+1);
	    sprintf(title,"H2D SDMass DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR)",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_highht",isb+1,ib+1);
	    sprintf(title,"H2D SDMass DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR)",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8_3[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8_3[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD",isb+1,ib+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR)",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8[ib][isb]->Sumw2();
	    
	    for(int ipt=0; ipt<(ntopptbins); ipt++){
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_Topptbin%i",isb+1,ib+1,ipt+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 subjetbtag%i AK4btag%i ToppT%i AK8CHS (DeepTag MD SR)",isb+1,ib+1,ipt+1);	 
	    
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_ptbin[ib][isb][ipt] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_ptbin[ib][isb][ipt]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_Topptbin%i_ttbar",isb+1,ib+1,ipt+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 subjetbtag%i AK4btag%i ToppT%i AK8CHS (DeepTag MD SR)",isb+1,ib+1,ipt+1);	 
	    
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_ptbin_tt[ib][isb][ipt] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_ptbin_tt[ib][isb][ipt]->Sumw2();
		
		}
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_highht",isb+1,ib+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR)",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_3[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_3[ib][isb]->Sumw2();
	    
	    for(int ijes=0; ijes<njesmax; ijes++){
			sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_JES%i_up",isb+1,ib+1,ijes+1);
			sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR) JES%i Up",isb+1,ib+1,ijes+1);	 
	    
			hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JES_up[ib][isb][ijes] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
			hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JES_up[ib][isb][ijes]->Sumw2();
			
			sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_JES%i_dn",isb+1,ib+1,ijes+1);
			sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR) JES%i Down",isb+1,ib+1,ijes+1);	 
	    
			hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JES_dn[ib][isb][ijes] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
			hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JES_dn[ib][isb][ijes]->Sumw2();
		}
	    
	    for(int ibunc=0; ibunc<nbtagmax; ibunc++){
			sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_BTag%i_up",isb+1,ib+1,ibunc);
			sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR) BTag%i Up",isb+1,ib+1,ibunc+1);	 
	    
			hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_Btag_up[ib][isb][ibunc] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
			hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_Btag_up[ib][isb][ibunc]->Sumw2();
			
			sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_BTag%i_dn",isb+1,ib+1,ibunc);
			sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR) BTag%i Down",isb+1,ib+1,ibunc+1);	 
	    
			hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_Btag_dn[ib][isb][ibunc] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
			hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_Btag_dn[ib][isb][ibunc]->Sumw2();
		}
	    
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_JER_up",isb+1,ib+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR) JER Up",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JER_up[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JER_up[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_JER_dn",isb+1,ib+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR) JER Down",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JER_dn[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JER_dn[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_JMS_up",isb+1,ib+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR) JMS Up",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JMS_up[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JMS_up[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_JMS_dn",isb+1,ib+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR) JMS Down",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JMS_dn[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JMS_dn[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_JMR_up",isb+1,ib+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR) JMR Up",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JMR_up[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JMR_up[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_JMR_dn",isb+1,ib+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR) JMR Down",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JMR_dn[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JMR_dn[ib][isb]->Sumw2();
	    
	    // t-tbar //
	    
	    sprintf(name,"H2D_SDMass_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_ttbar",isb+1,ib+1);
	    sprintf(title,"H2D SDMass DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR) t-#bar{t}",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8_tt[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8_tt[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_ttbar_highht",isb+1,ib+1);
	    sprintf(title,"H2D SDMass DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR) t-#bar{t}",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8_tt_3[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8_tt_3[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_ttbar",isb+1,ib+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR) t-#bar{t}",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_tt[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_tt[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_ttbar_highht",isb+1,ib+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR) t-#bar{t}",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_tt_3[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_tt_3[ib][isb]->Sumw2();
	    
	    // t-tbar ends 
	    
	    sprintf(name,"H2D_SDMass_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag",isb+1,ib+1);
	    sprintf(title,"H2D SDMass DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag SR)",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_AK8_DeepAK8[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
	    hist_2D_sdmass_deeptopscore_AK8_DeepAK8[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag",isb+1,ib+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 subjetbtag%i AK4btag%i AK8CHS (DeepTag SR)",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_md_AK8_DeepAK8[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,float(0),float(1.0));
	    hist_2D_sdmass_deeptopscore_md_AK8_DeepAK8[ib][isb]->Sumw2();
	   
		
			} //isb
		}//ib
		
	for(int isb=0; isb<nsbtag; isb++){
		for(int it=0; it<ntoptag; it++){
		
		sprintf(name,"H2D_AK4_pt_btag_TopTag%i_Subjetbtag%i_AK8CHS",it+1,isb+1);
		if(it==0) { sprintf(title,"2D Correlation of P_{T}  and b-tag of AK4CHS Jets | tight top subjet b-tag %i",isb+1); }
		if(it==1) { sprintf(title,"2D Correlation of P_{T}  and b-tag of AK4CHS Jets | medium top subjet b-tag %i",isb+1); }
		if(it==2) { sprintf(title,"2D Correlation of P_{T}  and b-tag of AK4CHS Jets | failed top subjet b-tag %i",isb+1); }
		
		hist_2D_pt_btag_AK4[isb][it] = new TH2D(name,title,noptbins,ptbins,20,0,1.0);
		hist_2D_pt_btag_AK4[isb][it]->Sumw2();
		
		
        sprintf(name,"H2D_AK4_pt_btag_TopTag%i_Subjetbtag%i_AK8CHS_lowbpt",it+1,isb+1);
		
		hist_2D_pt_btag_AK4_1[isb][it] = new TH2D(name,title,noptbins,ptbins,20,0.0,1.0);
		hist_2D_pt_btag_AK4_1[isb][it]->Sumw2();
		
		sprintf(name,"H2D_AK4_pt_btag_TopTag%i_Subjetbtag%i_AK8CHS_nobmasscut",it+1,isb+1);
		
		hist_2D_pt_btag_AK4_2[isb][it] = new TH2D(name,title,noptbins,ptbins,20,0.0,1.0);
		hist_2D_pt_btag_AK4_2[isb][it]->Sumw2();
		
		sprintf(name,"H2D_AK4_pt_btag_TopTag%i_Subjetbtag%i_AK8CHS_highht",it+1,isb+1);
		
		hist_2D_pt_btag_AK4_3[isb][it] = new TH2D(name,title,noptbins,ptbins,20,0.0,1.0);
		hist_2D_pt_btag_AK4_3[isb][it]->Sumw2();
		
		sprintf(name,"H2D_AK4_pt_btag_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS",it+1,isb+1);
		hist_2D_pt_btag_AK4_DeepAK8[isb][it] = new TH2D(name,title,noptbins,ptbins,20,0.0,1.0);
		hist_2D_pt_btag_AK4_DeepAK8[isb][it]->Sumw2();
		
		sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS",it+1,isb+1);
		hist_2D_pt_btag_AK4_md_DeepAK8[isb][it] = new TH2D(name,title,noptbins,ptbins,20,0.0,1.0);
		hist_2D_pt_btag_AK4_md_DeepAK8[isb][it]->Sumw2();
		
		for(int ibunc=0; ibunc<nbtagmax; ibunc++){
			sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_BTag%i_up",it+1,isb+1,ibunc);
			hist_2D_pt_btag_AK4_md_DeepAK8_Btag_up[isb][it][ibunc] = new TH2D(name,title,noptbins,ptbins,20,0.0,1.0);
			hist_2D_pt_btag_AK4_md_DeepAK8_Btag_up[isb][it][ibunc]->Sumw2();
	   
			sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_BTag%i_dn",it+1,isb+1,ibunc);
			hist_2D_pt_btag_AK4_md_DeepAK8_Btag_dn[isb][it][ibunc] = new TH2D(name,title,noptbins,ptbins,20,0.0,1.0);
			hist_2D_pt_btag_AK4_md_DeepAK8_Btag_dn[isb][it][ibunc]->Sumw2();
		}
		
		sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_highht",it+1,isb+1);
		hist_2D_pt_btag_AK4_md_DeepAK8_3[isb][it] = new TH2D(name,title,noptbins,ptbins,20,0.0,1.0);
		hist_2D_pt_btag_AK4_md_DeepAK8_3[isb][it]->Sumw2();
		
		sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_ttbar",it+1,isb+1);
		hist_2D_pt_btag_AK4_md_DeepAK8_tt[isb][it] = new TH2D(name,title,noptbins,ptbins,20,0.0,1.0);
		hist_2D_pt_btag_AK4_md_DeepAK8_tt[isb][it]->Sumw2();
		
		sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_ttbar_highht",it+1,isb+1);
		hist_2D_pt_btag_AK4_md_DeepAK8_tt_3[isb][it] = new TH2D(name,title,noptbins,ptbins,20,0.0,1.0);
		hist_2D_pt_btag_AK4_md_DeepAK8_tt_3[isb][it]->Sumw2();
		
		for(int ieta=0; ieta<netabins; ieta++){
			
			sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_bEtaBin%i",it+1,isb+1,ieta+1);
			hist_2D_pt_btag_AK4_md_DeepAK8_beta[isb][it][ieta] = new TH2D(name,title,noptbins,ptbins,20,0.0,1.0);
			hist_2D_pt_btag_AK4_md_DeepAK8_beta[isb][it][ieta]->Sumw2();
			
			sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_bEtaBin%i_2",it+1,isb+1,ieta+1);
			hist_2D_pt_btag_AK4_md_DeepAK8_beta2[isb][it][ieta] = new TH2D(name,title,noptbins,ptbins,20,0.0,1.0);
			hist_2D_pt_btag_AK4_md_DeepAK8_beta2[isb][it][ieta]->Sumw2();
			
			sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_bEtaBin%i_ttbar",it+1,isb+1,ieta+1);
			hist_2D_pt_btag_AK4_md_DeepAK8_beta_tt[isb][it][ieta] = new TH2D(name,title,noptbins,ptbins,20,0.0,1.0);
			hist_2D_pt_btag_AK4_md_DeepAK8_beta_tt[isb][it][ieta]->Sumw2();
			
			sprintf(name,"H2D_mtbt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_bEtaBin%i",it+1,isb+1,ieta+1);
			hist_2D_mtb_btag_AK4_md_DeepAK8_beta[isb][it][ieta] = new TH2D(name,title,nvwpmbin,wpmbins,20,0.0,1.0);
			hist_2D_mtb_btag_AK4_md_DeepAK8_beta[isb][it][ieta]->Sumw2();
			
			for(int ibunc=0; ibunc<nbtagmax; ibunc++){
				sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_bEtaBin%i_BTag%i_up",it+1,isb+1,ieta+1,ibunc);
				hist_2D_pt_btag_AK4_md_DeepAK8_beta_Btag_up[isb][it][ieta][ibunc] = new TH2D(name,title,noptbins,ptbins,20,0.0,1.0);
				hist_2D_pt_btag_AK4_md_DeepAK8_beta_Btag_up[isb][it][ieta][ibunc]->Sumw2();
	   
				sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_bEtaBin%i_BTag%i_dn",it+1,isb+1,ieta+1,ibunc);
				hist_2D_pt_btag_AK4_md_DeepAK8_beta_Btag_dn[isb][it][ieta][ibunc] = new TH2D(name,title,noptbins,ptbins,20,0.0,1.0);
				hist_2D_pt_btag_AK4_md_DeepAK8_beta_Btag_dn[isb][it][ieta][ibunc]->Sumw2();
					}
					
			for(int ijes=0; ijes<njesmax; ijes++){
				sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_bEtaBin%i_JES%i_up",it+1,isb+1,ieta+1,ijes+1);
				hist_2D_pt_btag_AK4_md_DeepAK8_beta_JES_up[isb][it][ieta][ijes] = new TH2D(name,title,noptbins,ptbins,20,0.0,1.0);
				hist_2D_pt_btag_AK4_md_DeepAK8_beta_JES_up[isb][it][ieta][ijes]->Sumw2();
	   
				sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_bEtaBin%i_JES%i_dn",it+1,isb+1,ieta+1,ijes+1);
				hist_2D_pt_btag_AK4_md_DeepAK8_beta_JES_dn[isb][it][ieta][ijes] = new TH2D(name,title,noptbins,ptbins,20,0.0,1.0);
				hist_2D_pt_btag_AK4_md_DeepAK8_beta_JES_dn[isb][it][ieta][ijes]->Sumw2();
					}		
			
				}
		
			}
		}
   
   ////
    
   ////
   
   sprintf(name,"Jet_pT_Top_Matched_Jet");
   sprintf(title,"P_{T} (GeV) of top matched jet");
   hist_jetpt_top = new TH1D(name,title,noptbins,ptbins);
   hist_jetpt_top->Sumw2();
   
   sprintf(name,"Jet_pT_Top_Matched_Jet_TopSDMass");
   sprintf(title,"P_{T} (GeV) of top matched jet (105<M_{SD}<210 GeV)");
   hist_jetpt_msd_top = new TH1D(name,title,noptbins,ptbins);
   hist_jetpt_msd_top->Sumw2();
   
   sprintf(name,"Jet_pT_Top_Matched_Jet_passing_MDDeepAK8");
   sprintf(title,"P_{T} (GeV) of top matched jet (passing (MD)DeepAK8 cut)");
   hist_jetpt_top_md_deepak8_pass = new TH1D(name,title,noptbins,ptbins);
   hist_jetpt_top_md_deepak8_pass->Sumw2();
   
   sprintf(name,"Jet_pT_Top_Matched_Jet_TopSDMass_passing_MDDeepAK8");
   sprintf(title,"P_{T} (GeV) of top matched jet (passing (MD)DeepAK8 cut & 105<M_{SD}<210 GeV)");
   hist_jetpt_msd_top_md_deepak8_pass = new TH1D(name,title,noptbins,ptbins);
   hist_jetpt_msd_top_md_deepak8_pass->Sumw2();
   
   sprintf(name,"DeltaR_AK4_TopTaggedAK8");
   sprintf(title,"#Delta R (AK4 jets, top candidate AK8 jet");
   hist_delR_AK4_toptagAK8 = new TH1D(name,title,400,0,4.0);
   hist_delR_AK4_toptagAK8->Sumw2();
   
   sprintf(name,"DeltaEta_AK4_TopTaggedAK8");
   sprintf(title,"#Delta Rapidity (AK4 jets, top candidate AK8 jet");
   hist_deleta_AK4_toptagAK8 = new TH1D(name,title,100,-5.,5.);
   hist_deleta_AK4_toptagAK8->Sumw2();
   
   sprintf(name,"DeltaPhi_AK4_TopTaggedAK8");
   sprintf(title,"#Delta #Phi (AK4 jets, top candidate AK8 jet");
   hist_delphi_AK4_toptagAK8 = new TH1D(name,title,150,-M_PI,M_PI);
   hist_delphi_AK4_toptagAK8->Sumw2();
   
   sprintf(name,"DeltaPhi_AK4_TopTaggedAK8_Passing_DeltaR");
   sprintf(title,"#Delta #Phi (AK4 jets, top candidate AK8 jet passing #DeltaR"); 
   hist_delphi_AK4_toptagAK8_dRpass = new TH1D(name,title,150,-M_PI,M_PI);
   hist_delphi_AK4_toptagAK8_dRpass->Sumw2();
   
   sprintf(name,"DeltaEta_AK4_TopTaggedAK8_Passing_DeltaR");
   sprintf(title,"#Delta Rapidity (AK4 jets, top candidate AK8 jet passing #DeltaR"); 
   hist_deleta_AK4_toptagAK8_dRpass = new TH1D(name,title,100,-5,5);
   hist_deleta_AK4_toptagAK8_dRpass->Sumw2();
   
   sprintf(name,"DeltaEta_AK4_TopTaggedAK8_Passing_DeltaR_mtb2000");
   sprintf(title,"#Delta Rapidity (AK4 jets, top candidate AK8 jet passing #DeltaR (M_{tb} > 2000 GeV)"); 
   hist_deleta_AK4_toptagAK8_dRpass_mtb2000 = new TH1D(name,title,100,-5,5);
   hist_deleta_AK4_toptagAK8_dRpass_mtb2000->Sumw2();
   
   sprintf(name,"DeltaR_BTaggedAK4_TopTaggedAK8");
   sprintf(title,"#Delta R (B Candidate AK4 jet, top candidate AK8 jet");
   hist_delR_btag_toptag = new TH1D(name,title,400,0,4.0);
   hist_delR_btag_toptag->Sumw2();
   
   sprintf(name,"DeltaPhi_BTaggedAK4_TopTaggedAK8");
   sprintf(title,"#Delta #Phi (B Candidate AK4 jet, top candidate AK8 jet");
   hist_delphi_btag_toptag = new TH1D(name,title,150,-M_PI,M_PI);
   hist_delphi_btag_toptag->Sumw2();
   
   
   
   sprintf(name,"NJets_AK4CHS");
   sprintf(title,"# of AK4CHS jets");
   hist_njetAK4 = new TH1D(name,title,10,-0.1,9.9);
   hist_njetAK4->Sumw2();
   
   sprintf(name,"NBJets_AK4CHS");
   sprintf(title,"# of B-Tagged AK4CHS jets");
   hist_nbjetAK4 = new TH1D(name,title,10,-0.1,9.9);
   hist_nbjetAK4->Sumw2();
   
   sprintf(name,"Pt_AK4CHS");
   sprintf(title,"P_{T} of AK4CHS jets");
   hist_jetptAK4 = new TH1D(name,title,noptbins,ptbins);
   hist_jetptAK4->Sumw2();
   
   sprintf(name,"Pt_AK4CHS_BCand");
   sprintf(title,"P_{T} of AK4CHS b candidate jet");
   hist_bjetptAK4 = new TH1D(name,title,noptbins,ptbins);
   hist_bjetptAK4->Sumw2();
   
   sprintf(name,"Pt_2D_AK8CHS_TopCand_AK4CHS_BCand");
   sprintf(title,"P_{T} Correlation of AK8CHS top and AK4CHS b candidate jets");
   hist_tjetptAK8_bjetptAK4 = new TH2D(name,title,noptbins,ptbins,noptbins,ptbins);
   hist_tjetptAK8_bjetptAK4->Sumw2();
   
   sprintf(name,"Mass_AK4CHS");
   sprintf(title,"Mass of AK4CHS jets");
   hist_jetmassAK4 = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetmassAK4->Sumw2();
  
   sprintf(name,"y_AK4CHS");
   sprintf(title,"Rapidity of AK4CHS jets");
   hist_jetrapAK4 = new TH1D(name,title,50,-5.,5.);
   hist_jetrapAK4->Sumw2();
   
   sprintf(name,"Btag_AK4CHS");
   sprintf(title,"CSVv2 Score of AK4CHS jets");
   hist_jetbtagAK4 = new TH1D(name,title,20,0,1.0);
   hist_jetbtagAK4->Sumw2();
   
   sprintf(name,"Btag_DeepCSV_AK4CHS");
   sprintf(title,"DeepCSV Score of AK4CHS jets");
   hist_jetbtagdeepCSVAK4 = new TH1D(name,title,20,0,1.0);
   hist_jetbtagdeepCSVAK4->Sumw2();
   
   sprintf(name,"Btag_DeepJetFlavB_AK4CHS");
   sprintf(title,"DeepJet B Flavour Score of AK4CHS jets");
   hist_jetbtagdeepflavAK4 = new TH1D(name,title,20,0,1.0);
   hist_jetbtagdeepflavAK4->Sumw2();
   
   sprintf(name,"Parton_Flavour_AK4 jets");
   sprintf(title,"Parton Flavour of AK4 jets");
   hist_jetpartonflavAK4 = new TH1D(name,title,22,-0.01,21.99);
   hist_jetpartonflavAK4->Sumw2();
   
   sprintf(name,"Parton_Flavour_BTagged_AK4 jets");
   sprintf(title,"Parton Flavour of BTagged AK4 jets");
   hist_jetpartonflavAK4_btagged = new TH1D(name,title,22,-0.01,21.99);
   hist_jetpartonflavAK4_btagged->Sumw2();
   
   sprintf(name,"Mass_AK4CHS_highpt");
   sprintf(title,"Mass of AK4CHS jets");
   hist_jetmassAK4_1 = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetmassAK4_1->Sumw2();
   
   sprintf(name,"Btag_AK4CHS_highpt");
   sprintf(title,"CSVv2 Score of AK4CHS jets");
   hist_jetbtagAK4_1 = new TH1D(name,title,20,0,1.0);
   hist_jetbtagAK4_1->Sumw2();
   
   sprintf(name,"Btag_DeepCSV_AK4CHS_highpt");
   sprintf(title,"DeepCSV Score of AK4CHS jets");
   hist_jetbtagdeepCSVAK4_1 = new TH1D(name,title,20,0,1.0);
   hist_jetbtagdeepCSVAK4_1->Sumw2();
   
   sprintf(name,"Btag_DeepJetFlavB_AK4CHS_highpt");
   sprintf(title,"DeepJet B Flavour Score of AK4CHS jets");
   hist_jetbtagdeepflavAK4_1 = new TH1D(name,title,20,0,1.0);
   hist_jetbtagdeepflavAK4_1->Sumw2();
   
   sprintf(name,"Parton_Flavour_BTagged_AK4 jets_highbpt");
   sprintf(title,"Parton Flavour of BTagged AK4 jets");
   hist_jetpartonflavAK4_btagged_1 = new TH1D(name,title,22,-0.01,21.99);
   hist_jetpartonflavAK4_btagged_1->Sumw2();
   
   sprintf(name,"Parton_Flavour_AK4 jets_highpt");
   sprintf(title,"Parton Flavour of AK4 jets");
   hist_jetpartonflavAK4_1 = new TH1D(name,title,22,-0.01,21.99);
   hist_jetpartonflavAK4_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8JetMass");
   sprintf(title,"Soft drop mass of AK8 jet associated to b-Jet");
   hist_biso_mass = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_biso_mass->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8AK4JetMassDiff");
   sprintf(title,"Difference of Mass of AK8 jet associated to b-Jet and b-Jet");
   hist_biso_tbmassdiff = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_biso_tbmassdiff->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK4byAK8_Mass");
   sprintf(title,"Ratio of Soft drop mass of AK8 jet associated to b-Jet to b-Jet Mass");
   hist_biso_isomass = new TH1D(name,title,20,-0.001,0.999);
   hist_biso_isomass->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK4byAK8_Pt");
   sprintf(title,"Ratio of P_{T} of AK8 jet associated to b-Jet to b-Jet");
   hist_biso_isopt = new TH1D(name,title,20,-0.001,0.999);
   hist_biso_isopt->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_TopTagScore");
   sprintf(title,"Top-tag score of AK8 jet associated to b-Jet to b-Jet");
   hist_biso_TopAK8score = new TH1D(name,title,100,-0.001,0.999);
   hist_biso_TopAK8score->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_MDTopTagScore");
   sprintf(title,"MD Top-tag score of AK8 jet associated to b-Jet to b-Jet");
   hist_biso_TopAK8score_MD = new TH1D(name,title,100,-0.001,0.999);
   hist_biso_TopAK8score_MD->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_WTagScore");
   sprintf(title,"W-tag score of AK8 jet associated to b-Jet to b-Jet");
   hist_biso_WAK8score = new TH1D(name,title,100,-0.001,0.999);
   hist_biso_WAK8score->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_MDWTagScore");
   sprintf(title,"MD W-tag score of AK8 jet associated to b-Jet to b-Jet");
   hist_biso_WAK8score_MD = new TH1D(name,title,100,-0.001,0.999);
   hist_biso_WAK8score_MD->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_GenJetMass");
   sprintf(title,"Mass of GEN AK8 jet associated to b-Jet to b-Jet");
   hist_bjetAK8genmass = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_bjetAK8genmass->Sumw2();
   
   for(int it=0; it<ntoptag; it++){
	   for(int isb=0; isb<nsbtag; isb++){
   
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK4Mass_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK4 Jet Mass, corresponding to b-Cand");
		hist_2d_biso_AK4mass_pt[it][isb] = new TH2D(name,title,noptbins,ptbins,nsdmassbins,sdmassbins);
		hist_2d_biso_AK4mass_pt[it][isb]->Sumw2();
   
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK8SDMass_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet Soft drop mass, corresponding to b-Cand");
		hist_2d_biso_AK8sdmass_pt[it][isb] = new TH2D(name,title,noptbins,ptbins,nsdmassbins,sdmassbins);
		hist_2d_biso_AK8sdmass_pt[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK8TopScore_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet Top-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8TopScore_pt[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8TopScore_pt[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_MDAK8TopScore_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet MD Top-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8TopScore_MD_pt[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8TopScore_MD_pt[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK8WScore_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet W-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8WScore_pt[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8WScore_pt[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_MDAK8WScore_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet MD W-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8WScore_MD_pt[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8WScore_MD_pt[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK4Mass_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK4 Jet Mass, corresponding to b-Cand");
		hist_2d_biso_AK4mass_pt_dAK8[it][isb] = new TH2D(name,title,noptbins,ptbins,nsdmassbins,sdmassbins);
		hist_2d_biso_AK4mass_pt_dAK8[it][isb]->Sumw2();
   
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK8SDMass_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet Soft drop mass, corresponding to b-Cand");
		hist_2d_biso_AK8sdmass_pt_dAK8[it][isb] = new TH2D(name,title,noptbins,ptbins,nsdmassbins,sdmassbins);
		hist_2d_biso_AK8sdmass_pt_dAK8[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK8TopScore_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet Top-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8TopScore_pt_dAK8[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8TopScore_pt_dAK8[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_MDAK8TopScore_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet Top-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8TopScore_MD_pt_dAK8[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8TopScore_MD_pt_dAK8[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK8WScore_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet W-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8WScore_pt_dAK8[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8WScore_pt_dAK8[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_MDAK8WScore_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet MD W-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8WScore_MD_pt_dAK8[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8WScore_MD_pt_dAK8[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK4Mass_MD_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK4 Jet Mass, corresponding to b-Cand");
		hist_2d_biso_AK4mass_pt_mdAK8[it][isb] = new TH2D(name,title,noptbins,ptbins,nsdmassbins,sdmassbins);
		hist_2d_biso_AK4mass_pt_mdAK8[it][isb]->Sumw2();
   
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK8SDMass_MD_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet Soft drop mass, corresponding to b-Cand");
		hist_2d_biso_AK8sdmass_pt_mdAK8[it][isb] = new TH2D(name,title,noptbins,ptbins,nsdmassbins,sdmassbins);
		hist_2d_biso_AK8sdmass_pt_mdAK8[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK8TopScore_MD_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet Top-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8TopScore_pt_mdAK8[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8TopScore_pt_mdAK8[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_MDAK8TopScore_MD_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet MD Top-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8TopScore_MD_pt_mdAK8[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8TopScore_MD_pt_mdAK8[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK8WScore_MD_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet W-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8WScore_pt_mdAK8[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8WScore_pt_mdAK8[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_MDAK8WScore_MD_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet MD W-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8WScore_MD_pt_mdAK8[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8WScore_MD_pt_mdAK8[it][isb]->Sumw2();
		
		}
	}
   
   // high b pt //
   
   sprintf(name,"IsoCheck_bjet_AK8JetMass_highbpt");
   sprintf(title,"Soft drop mass of AK8 jet associated to b-Jet");
   hist_biso_mass_1 = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_biso_mass_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8AK4JetMassDiff_highbpt");
   sprintf(title,"Difference of Mass of AK8 jet associated to b-Jet and b-Jet");
   hist_biso_tbmassdiff_1 = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_biso_tbmassdiff_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK4byAK8_Mass_highbpt");
   sprintf(title,"Ratio of Soft drop mass of AK8 jet associated to b-Jet to b-Jet Mass");
   hist_biso_isomass_1 = new TH1D(name,title,20,-0.001,0.999);
   hist_biso_isomass_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK4byAK8_Pt_highbpt");
   sprintf(title,"Ratio of P_{T} of AK8 jet associated to b-Jet to b-Jet");
   hist_biso_isopt_1 = new TH1D(name,title,20,-0.001,0.999);
   hist_biso_isopt_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_TopTagScore_highbpt");
   sprintf(title,"Top-tag score of AK8 jet associated to b-Jet to b-Jet");
   hist_biso_TopAK8score_1 = new TH1D(name,title,100,-0.001,0.999);
   hist_biso_TopAK8score_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_MDTopTagScore_highbpt");
   sprintf(title,"MD Top-tag score of AK8 jet associated to b-Jet to b-Jet");
   hist_biso_TopAK8score_MD_1 = new TH1D(name,title,100,-0.001,0.999);
   hist_biso_TopAK8score_MD_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_TopTagScore_masscut_highbpt");
   sprintf(title,"Top-tag score of AK8 jet associated to b-Jet to b-Jet (M_{SD} > 60 GeV)");
   hist_biso_TopAK8score_masscut_1 = new TH1D(name,title,100,-0.001,0.999);
   hist_biso_TopAK8score_masscut_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_MDTopTagScore_masscut_highbpt");
   sprintf(title,"MD Top-tag score of AK8 jet associated to b-Jet to b-Jet (M_{SD} > 60 GeV)");
   hist_biso_TopAK8score_MD_masscut_1 = new TH1D(name,title,100,-0.001,0.999);
   hist_biso_TopAK8score_MD_masscut_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_WTagScore_highbpt");
   sprintf(title,"W-tag score of AK8 jet associated to b-Jet to b-Jet");
   hist_biso_WAK8score_1 = new TH1D(name,title,100,-0.001,0.999);
   hist_biso_WAK8score_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_MDWTagScore_highbpt");
   sprintf(title,"W-tag score of AK8 jet associated to b-Jet to b-Jet");
   hist_biso_WAK8score_MD_1 = new TH1D(name,title,100,-0.001,0.999);
   hist_biso_WAK8score_MD_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_GenJetMass_highpt");
   sprintf(title,"Mass of GEN AK8 jet associated to b-Jet to b-Jet");
   hist_bjetAK8genmass_1 = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_bjetAK8genmass_1->Sumw2();
   
   for(int ib=0; ib<nbtag; ib++){
	   sprintf(name,"IsoCheck_bjet_AK8JetMass_highbpt_btag%i",ib);
	   sprintf(title,"Soft drop mass of AK8 jet associated to b-Jet B tag %i",ib);
	   hist_biso_mass_btag_1[ib] = new TH1D(name,title,nsdmassbins,sdmassbins);
	   hist_biso_mass_btag_1[ib]->Sumw2();
	   
	   sprintf(name,"IsoCheck_bjet_AK8_MDTopTagScore_highbpt_btag%i",ib);
	   sprintf(title,"MD Top-tag score of AK8 jet associated to b-Jet to b-Jet B tag %i",ib);
	   hist_biso_TopAK8score_MD_btag_1[ib] = new TH1D(name,title,100,-0.001,0.999);
	   hist_biso_TopAK8score_MD_btag_1[ib]->Sumw2();
	   
	   sprintf(name,"IsoCheck_bjet_AK8_MDWTagScore_highbpt_btag%i",ib);
	   sprintf(title,"W-tag score of AK8 jet associated to b-Jet to b-Jet B tag %i",ib);
	   hist_biso_WAK8score_MD_btag_1[ib] = new TH1D(name,title,100,-0.001,0.999);
	   hist_biso_WAK8score_MD_btag_1[ib]->Sumw2();
	}
   
   for(int it=0; it<ntoptag; it++){
	   for(int isb=0; isb<nsbtag; isb++){
   
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK4Mass_highbpt_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK4 Jet Mass, corresponding to b-Cand");
		hist_2d_biso_AK4mass_pt_1[it][isb] = new TH2D(name,title,noptbins,ptbins,nsdmassbins,sdmassbins);
		hist_2d_biso_AK4mass_pt_1[it][isb]->Sumw2();
   
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK8SDMass_highbpt_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet Soft drop mass, corresponding to b-Cand");
		hist_2d_biso_AK8sdmass_pt_1[it][isb] = new TH2D(name,title,noptbins,ptbins,nsdmassbins,sdmassbins);
		hist_2d_biso_AK8sdmass_pt_1[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK8TopScore_highbpt_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet Top-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8TopScore_pt_1[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8TopScore_pt_1[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_MDAK8TopScore_highbpt_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet MD Top-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8TopScore_MD_pt_1[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8TopScore_MD_pt_1[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK8WScore_highbpt_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet W-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8WScore_pt_1[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8WScore_pt_1[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_MDAK8WScore_highbpt_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet MD W-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8WScore_MD_pt_1[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8WScore_MD_pt_1[it][isb]->Sumw2();
		
		//deepak8
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK4Mass_highbpt_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK4 Jet Mass, corresponding to b-Cand");
		hist_2d_biso_AK4mass_pt_dAK8_1[it][isb] = new TH2D(name,title,noptbins,ptbins,nsdmassbins,sdmassbins);
		hist_2d_biso_AK4mass_pt_dAK8_1[it][isb]->Sumw2();
   
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK8SDMass_highbpt_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet Soft drop mass, corresponding to b-Cand");
		hist_2d_biso_AK8sdmass_pt_dAK8_1[it][isb] = new TH2D(name,title,noptbins,ptbins,nsdmassbins,sdmassbins);
		hist_2d_biso_AK8sdmass_pt_dAK8_1[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK8TopScore_highbpt_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet Top-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8TopScore_pt_dAK8_1[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8TopScore_pt_dAK8_1[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_MDAK8TopScore_highbpt_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet MD Top-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8TopScore_MD_pt_dAK8_1[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8TopScore_MD_pt_dAK8_1[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK8WScore_highbpt_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet W-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8WScore_pt_dAK8_1[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8WScore_pt_dAK8_1[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_MDAK8WScore_highbpt_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet MD W-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8WScore_MD_pt_dAK8_1[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8WScore_MD_pt_dAK8_1[it][isb]->Sumw2();
		
		//deepak8 (md)
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK4Mass_highbpt_MD_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK4 Jet Mass, corresponding to b-Cand");
		hist_2d_biso_AK4mass_pt_mdAK8_1[it][isb] = new TH2D(name,title,noptbins,ptbins,nsdmassbins,sdmassbins);
		hist_2d_biso_AK4mass_pt_mdAK8_1[it][isb]->Sumw2();
   
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK8SDMass_highbpt_MD_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet Soft drop mass, corresponding to b-Cand");
		hist_2d_biso_AK8sdmass_pt_mdAK8_1[it][isb] = new TH2D(name,title,noptbins,ptbins,nsdmassbins,sdmassbins);
		hist_2d_biso_AK8sdmass_pt_mdAK8_1[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK8TopScore_highbpt_MD_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet Top-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8TopScore_pt_mdAK8_1[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8TopScore_pt_mdAK8_1[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_MDAK8TopScore_highbpt_MD_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet MD Top-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8TopScore_MD_pt_mdAK8_1[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8TopScore_MD_pt_mdAK8_1[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_AK8WScore_highbpt_MD_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet W-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8WScore_pt_mdAK8_1[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8WScore_pt_mdAK8_1[it][isb]->Sumw2();
		
		sprintf(name,"IsoCheck_bjet_AK8Pt_MDAK8WScore_highbpt_MD_DeepAK8TopTag_TopTag%i_Subbtag%i",it,isb);
		sprintf(title,"2D Correlation of AK8 jet p_{T} and AK8 jet MD W-tag Score Mass, corresponding to b-Cand");
		hist_2d_biso_AK8WScore_MD_pt_mdAK8_1[it][isb] = new TH2D(name,title,noptbins,ptbins,100,-0.001,0.999);
		hist_2d_biso_AK8WScore_MD_pt_mdAK8_1[it][isb]->Sumw2();
   
		}
	}
   
   sprintf(name,"IsoCheck_bjet_AK8_SubJetMass_12");
   sprintf(title,"2D Correlation of subjet Mass AK8 jet, corresponding to b-Cand");
   hist_2D_biso_AK8jet_subjetmass12 = new TH2D(name,title,nsdmassbins,sdmassbins,nsdmassbins,sdmassbins);
   hist_2D_biso_AK8jet_subjetmass12->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_SubJetMass1_JetSDMass");
   sprintf(title,"2D Correlation of higher subjet Mass to AK8 jet SD Mass, corresponding to b-Cand");
   hist_2D_biso_AK8jetsdmass_subjetmass1 = new TH2D(name,title,nsdmassbins,sdmassbins,nsdmassbins,sdmassbins);
   hist_2D_biso_AK8jetsdmass_subjetmass1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_SubJetMass2_JetSDMass");
   sprintf(title,"2D Correlation of lower subjet Mass to AK8 jet SD Mass, corresponding to b-Cand");
   hist_2D_biso_AK8jetsdmass_subjetmass2 = new TH2D(name,title,nsdmassbins,sdmassbins,nsdmassbins,sdmassbins);
   hist_2D_biso_AK8jetsdmass_subjetmass2->Sumw2();
   
   sprintf(name,"Correlation_GenMass_AK8top_bAK8_highpt");
   sprintf(title,"2D Correlation of top candidate gen mass and gen mass of AK8 jet corresponding to b-Cand");
   hist_genmass_tbcor = new TH2D(name,title,nsdmassbins,sdmassbins,nsdmassbins,sdmassbins);
   hist_genmass_tbcor->Sumw2();
   
   // high b pt end //
   
   for(int it=0; it<ntoptag; it++){
    for(int ib=0; ib<nbtag; ib++){
		for(int isb=0; isb<nsbtag; isb++){
   
   sprintf(name,"tb_Mass_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Pt_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Pt of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbpt[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_tbpt[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_y_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Rapidity of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbrap[it][ib][isb] = new TH1D(name,title,50,-5.,5.);
   hist_tbrap[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_z_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Min(p_{T}^{t},p_{T}^{b}) / (p_{T}^{t}+p_{T}^{b})  | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1); 
   hist_tbpt_frac[it][ib][isb] = new TH1D(name,title,50,0,0.5);
   hist_tbpt_frac[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_deltheta_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"cos(#Delta#theta^{tb}) | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbtheta_diff[it][ib][isb] = new TH1D(name,title,50,-1.,+1.);
   hist_tbtheta_diff[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_delphi_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"#Delta#phi^{tb} | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbphi_diff[it][ib][isb] = new TH1D(name,title,150,-M_PI,M_PI);
   hist_tbphi_diff[it][ib][isb]->Sumw2();
   
   sprintf(name,"BCandidate_Jet_Parton_Flavour_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"BCandidate Jet Parton Flavour | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_partonflav_AK4[it][ib][isb] = new TH1D(name,title,22,-0.01,21.99);
   hist_partonflav_AK4[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Pt_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet P_{T} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetpt_sel[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_topjetpt_sel[it][ib][isb]->Sumw2();
    
   sprintf(name,"TopJet_SDMass_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_DeepTagTvsQCD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetdeeptopscore_sel[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeeptopscore_sel[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_MD_DeepTagTvsQCD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet Mass-Decorrelated DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetdeepmdtopscore_sel[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeepmdtopscore_sel[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet SubJet CSVv2 Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsbtag_sel[it][ib][isb] = new TH1D(name,title,20,0,1.0);
   hist_topjetsbtag_sel[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepCSV_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet SubJet DeepCSV Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdeepCSV_sel[it][ib][isb] = new TH1D(name,title,20,0,1.0);
   hist_topjetsdeepCSV_sel[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Tau32_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet #tau_{32} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjettau32_sel[it][ib][isb] = new TH1D(name,title,50,0,1.0);
   hist_topjettau32_sel[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_Pt_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"BJet P_{T} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetpt_sel[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_bjetpt_sel[it][ib][isb]->Sumw2();
    
   sprintf(name,"BJet_Mass_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"BJet Mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmass_sel[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_bjetmass_sel[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_BTag_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"BJet CSVv2 Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetbtag_sel[it][ib][isb] = new TH1D(name,title,20,0,1.0);
   hist_bjetbtag_sel[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_TopTagScore_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Top-tag Score of AK8 jet corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Topscore_sel[it][ib][isb] = new TH1D(name,title,100,0,1.0);
   hist_bjetmatch_AK8Topscore_sel[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_WTagScore_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"W-tag Score of AK8 jet corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Wscore_sel[it][ib][isb] = new TH1D(name,title,100,0,1.0);
   hist_bjetmatch_AK8Wscore_sel[it][ib][isb]->Sumw2();
   
   sprintf(name,"Event_HT_DeepAK8MD_toptag%i_btag%i_subbtag%i_tau32",it,ib,isb);
   sprintf(title,"H_{T} | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_ht_tau32[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_ht_tau32[it][ib][isb]->Sumw2();

   sprintf(name,"TopJet_SDMass_Tau32_tbMass_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet Soft drop mass #tau{32} tb Inv Mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_3D_topjetsdmass_topjettau32_tbmass[it][ib][isb] = new TH3D(name,title,nomassbins,mass_low,mass_high,50,0,1.0,nwpmbin,nwpmlow,nwpmhigh);
   hist_3D_topjetsdmass_topjettau32_tbmass[it][ib][isb]->Sumw2();
   
   // AK8 b Candidate//
   
   sprintf(name,"tb_Mass_AK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_AK8[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_AK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJetAK8_DeepTagTvsQCD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"AK8 b-Jet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_AK8bCand_DeeptagScore[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   
   sprintf(name,"BJetAK8_DeepTagTvsQCD_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"AK8 b-Jet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_AK8bCand_DeeptagScore_1[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   
   sprintf(name,"BJetAK8_MD_DeepTagTvsQCD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"AK8 b-Jet Mass Decorrelated DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_AK8bCand_MD_DeeptagScore[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   
   sprintf(name,"tb_Mass_AK8_DeepAK8TopTag_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_AK8_DeepAK8[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_AK8_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJetAK8_DeepTagTvsQCD_DeepAK8TopTag_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"AK8 b-Jet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_AK8bCand_DeeptagScore_DeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_AK8bCand_DeeptagScore_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJetAK8_DeepTagTvsQCD_DeepAK8TopTag_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"AK8 b-Jet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_AK8bCand_DeeptagScore_DeepAK8_1[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_AK8bCand_DeeptagScore_DeepAK8_1[it][ib][isb]->Sumw2();
  
   sprintf(name,"BJetAK8_MD_DeepTagTvsQCD_DeepAK8TopTag_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"AK8 b-Jet Mass Decorrelated DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_AK8bCand_MD_DeeptagScore_DeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_AK8bCand_MD_DeeptagScore_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_AK8_MD_DeepAK8TopTag_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_AK8_MDDeepAK8[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_AK8_MDDeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJetAK8_DeepTagTvsQCD_MD_DeepAK8TopTag_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"AK8 b-Jet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_AK8bCand_DeeptagScore_MDDeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_AK8bCand_DeeptagScore_MDDeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJetAK8_DeepTagTvsQCD_MD_DeepAK8TopTag_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"AK8 b-Jet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_AK8bCand_DeeptagScore_MDDeepAK8_1[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_AK8bCand_DeeptagScore_MDDeepAK8_1[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJetAK8_MD_DeepTagTvsQCD_MD_DeepAK8TopTag_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"AK8 b-Jet Mass Decorrelated DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_AK8bCand_MD_DeeptagScore_MDDeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_AK8bCand_MD_DeeptagScore_MDDeepAK8[it][ib][isb]->Sumw2();
   
			}//isb
		}//ib
	}//it


for(int it=0; it<ntoptag_DAK8; it++){
  for(int ib=0; ib<nbtag; ib++){
	for(int isb=0; isb<nsbtag; isb++){
   
   sprintf(name,"tb_Mass_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_DeepAK8[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Pt_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Pt of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbpt_DeepAK8[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_tbpt_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_y_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Rapidity of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbrap_DeepAK8[it][ib][isb] = new TH1D(name,title,50,-5.,5.);
   hist_tbrap_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BCandidate_Jet_Parton_Flavour_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"BCandidate Jet Parton Flavour | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_partonflav_AK4_DeepAK8[it][ib][isb] = new TH1D(name,title,22,-0.01,21.99);
   hist_partonflav_AK4_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Pt_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet P_{T} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetpt_sel_DeepAK8[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_topjetpt_sel_DeepAK8[it][ib][isb]->Sumw2();
   
   for(int  imsd=0; imsd<ntopsdmassbins; imsd++){
   sprintf(name,"TopJet_Pt_DeepAK8_toptag%i_btag%i_subbtag%i_msdbin%i",it,ib,isb,imsd);
   sprintf(title,"TopJet P_{T} Top WP %i B WP %i SubB WP %i M_{SD} bin %i",it+1,ib+1,isb+1,imsd+1); 
   hist_topjetpt_msdbin_sel_DeepAK8[it][ib][isb][imsd] = new TH1D(name,title,noptbins,ptbins);
   hist_topjetpt_msdbin_sel_DeepAK8[it][ib][isb][imsd]->Sumw2();
   }
    
   sprintf(name,"TopJet_SDMass_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel_DeepAK8[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_DeepTagTvsQCD_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetdeeptopscore_sel_DeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeeptopscore_sel_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_MD_DeepTagTvsQCD_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet Mass-Decorrelated DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetdeepmdtopscore_sel_DeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeepmdtopscore_sel_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet SubJet CSVv2 Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsbtag_sel_DeepAK8[it][ib][isb] = new TH1D(name,title,20,0,1.0);
   hist_topjetsbtag_sel_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepCSV_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet SubJet DeepCSV Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdeepCSV_sel_DeepAK8[it][ib][isb] = new TH1D(name,title,20,0,1.0);
   hist_topjetsdeepCSV_sel_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Tau32_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet #tau_{32} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjettau32_sel_DeepAK8[it][ib][isb] = new TH1D(name,title,50,0,1.0);
   hist_topjettau32_sel_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_Pt_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"BJet P_{T} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetpt_sel_DeepAK8[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_bjetpt_sel_DeepAK8[it][ib][isb]->Sumw2();
    
   sprintf(name,"BJet_Mass_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"BJet Mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmass_sel_DeepAK8[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_bjetmass_sel_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_BTag_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"BJet CSVv2 Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetbtag_sel_DeepAK8[it][ib][isb] = new TH1D(name,title,20,0,1.0);
   hist_bjetbtag_sel_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_TopTagScore_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Top-tag Score of AK8 jet corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Topscore_sel_DeepAK8[it][ib][isb] = new TH1D(name,title,100,0,1.0);
   hist_bjetmatch_AK8Topscore_sel_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_WTagScore_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"W-tag Score of AK8 jet corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Wscore_sel_DeepAK8[it][ib][isb] = new TH1D(name,title,100,0,1.0);
   hist_bjetmatch_AK8Wscore_sel_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"Event_HT_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"H_{T} | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_ht_DeepAK8[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_ht_DeepAK8[it][ib][isb]->Sumw2();
   
   // Mass Decorrelated //
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_option2",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (optional p/f ratio)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_2[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_passing_dy_cut",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (#Delta y <1.8)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_dY[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_dY[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_trigwgt",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (trigger weight)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_trigwt[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_trigwt[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_pf_up",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (p/f up)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_pfup[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_pfup[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_pf_dn",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (p/f down)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_pfdn[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_pfdn[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_pf2_up",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (p/f up : option 2)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_pfup_2[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_pfup_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_pf2_dn",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (p/f down : option 2)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_pfdn_2[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_pfdn_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_mdAK8_up",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (DeepAK8 SF up)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_dAK8up[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_dAK8up[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_mdAK8_dn",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (DeepAK8 SF down)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_dAK8dn[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_dAK8dn[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_pu_up",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (PU up)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_puup[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_puup[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_pu_dn",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (PU down)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_pudn[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_pudn[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_prefire_up",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (prefire up)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_prefireup[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_prefireup[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_prefire_dn",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (prefire down)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_prefiredn[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_prefiredn[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_fine",it,ib,isb);
   hist_tbmass_md_DeepAK8_fine[it][ib][isb] = new TH1D(name,title,100,1000,5000);
   hist_tbmass_md_DeepAK8_fine[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_fine_puup",it,ib,isb);
   hist_tbmass_md_DeepAK8_fine_puup[it][ib][isb] = new TH1D(name,title,100,1000,5000);
   hist_tbmass_md_DeepAK8_fine_puup[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_fine_pudn",it,ib,isb);
   hist_tbmass_md_DeepAK8_fine_pudn[it][ib][isb] = new TH1D(name,title,100,1000,5000);
   hist_tbmass_md_DeepAK8_fine_pudn[it][ib][isb]->Sumw2();
   
   sprintf(name,"TwoD_tb_Mass_npv_DeepAK8MD_toptag%i_btag%i_subbtag%i_fine",it,ib,isb);
   sprintf(title,"2D correlation of M_{tb} vs n_{PU}");
   h2d_mtb_npu[it][ib][isb] = new TH2D(name,name,100,1000,5000,100,0,100);
   h2d_mtb_npu[it][ib][isb]->Sumw2();
   
   
   if(FakeAnalysis){
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_bmsdcor_up",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (b m_{sd} correction up)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_bcorup[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_bcorup[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_bmsdcor_dn",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (b m_{sd} correction down)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_bcordn[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_bcordn[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_bmsdcor_up2",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (b m_{sd} correction up)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_bcorup_2[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_bcorup_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_bmsdcor_dn2",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (b m_{sd} correction down)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_bcordn_2[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_bcordn_2[it][ib][isb]->Sumw2();
   
   }
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_noptw",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (w/o top p_{T} weight)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_noptw[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_noptw[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_alphaup",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_alphaup[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_alphaup[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_alphadn",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_alphadn[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_alphadn[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_betaup",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_betaup[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_betaup[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_betadn",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_betadn[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_betadn[it][ib][isb]->Sumw2();
   
   for(int ipt=0; ipt<ntopptbins; ipt++){
	sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_Topptbin%i",it,ib,isb,ipt);
	sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i ToppT %i",it+1,ib+1,isb+1,ipt+1);
	hist_tbmass_md_DeepAK8_ptbin[it][ib][isb][ipt] = new TH1D(name,title,nvwpmbin,wpmbins);
	hist_tbmass_md_DeepAK8_ptbin[it][ib][isb][ipt]->Sumw2();
	
	sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_Topptbin%i_ttbar",it,ib,isb,ipt);
	sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i ToppT %i",it+1,ib+1,isb+1,ipt+1);
	hist_tbmass_md_DeepAK8_ptbin_tt[it][ib][isb][ipt] = new TH1D(name,title,nvwpmbin,wpmbins);
	hist_tbmass_md_DeepAK8_ptbin_tt[it][ib][isb][ipt]->Sumw2();
   }
   
   for(int ipdf=0; ipdf<npdfmax; ipdf++){
	   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_PDF%i",it,ib,isb,ipdf);
	   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i PDF %i",it+1,ib+1,isb+1,ipdf+1);
	   hist_tbmass_md_DeepAK8_PDF[it][ib][isb][ipdf] = new TH1D(name,title,nvwpmbin,wpmbins);
	   hist_tbmass_md_DeepAK8_PDF[it][ib][isb][ipdf]->Sumw2();
	}
	   
   for(int iscale=0; iscale<nscalemax; iscale++){
	   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_Scale%i",it,ib,isb,iscale);
	   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i Scale %i",it+1,ib+1,isb+1,iscale+1);
	   hist_tbmass_md_DeepAK8_Scale[it][ib][isb][iscale] = new TH1D(name,title,nvwpmbin,wpmbins);
	   hist_tbmass_md_DeepAK8_Scale[it][ib][isb][iscale]->Sumw2();
   }
   
   for(int ips=0; ips<npsmax; ips++){
	   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_PS%i",it,ib,isb,ips);
	   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i PS %i",it+1,ib+1,isb+1,ips+1);
	   hist_tbmass_md_DeepAK8_PS[it][ib][isb][ips] = new TH1D(name,title,nvwpmbin,wpmbins);
	   hist_tbmass_md_DeepAK8_PS[it][ib][isb][ips]->Sumw2();
   }
   
   for(int ibunc=0; ibunc<nbtagmax; ibunc++){
	   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_BTag%i_up",it,ib,isb,ibunc);
	   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i BTag %i (Up)",it+1,ib+1,isb+1,ibunc+1);
	   hist_tbmass_md_DeepAK8_Btag_up[it][ib][isb][ibunc] = new TH1D(name,title,nvwpmbin,wpmbins);
	   hist_tbmass_md_DeepAK8_Btag_up[it][ib][isb][ibunc]->Sumw2();
	   
	   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_BTag%i_dn",it,ib,isb,ibunc);
	   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i BTag %i (Down)",it+1,ib+1,isb+1,ibunc+1);
	   hist_tbmass_md_DeepAK8_Btag_dn[it][ib][isb][ibunc] = new TH1D(name,title,nvwpmbin,wpmbins);
	   hist_tbmass_md_DeepAK8_Btag_dn[it][ib][isb][ibunc]->Sumw2();
   }
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_JER_up",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i JER (Up)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_JER_up[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_JER_up[it][ib][isb]->Sumw2();
	  
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_JER_dn",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i JER (Down)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_JER_dn[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_JER_dn[it][ib][isb]->Sumw2();	  
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_JMS_up",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i JMS (Up)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_JMS_up[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_JMS_up[it][ib][isb]->Sumw2();
	  
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_JMS_dn",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i JMS (Down)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_JMS_dn[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_JMS_dn[it][ib][isb]->Sumw2();	   
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_JMR_up",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i JMR (Up)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_JMR_up[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_JMR_up[it][ib][isb]->Sumw2();
	  
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_JMR_dn",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i JMR (Down)",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_JMR_dn[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_JMR_dn[it][ib][isb]->Sumw2();	  
   
   for(int ijes=0; ijes<njesmax; ijes++){
	   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_JES%i_up",it,ib,isb,ijes+1);
	   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i JES %i (Up)",it+1,ib+1,isb+1,ijes+1);
	   hist_tbmass_md_DeepAK8_JES_up[it][ib][isb][ijes] = new TH1D(name,title,nvwpmbin,wpmbins);
	   hist_tbmass_md_DeepAK8_JES_up[it][ib][isb][ijes]->Sumw2();
	  
	   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_JES%i_dn",it,ib,isb,ijes+1);
	   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i JES %i (Down)",it+1,ib+1,isb+1,ijes+1);
	   hist_tbmass_md_DeepAK8_JES_dn[it][ib][isb][ijes] = new TH1D(name,title,nvwpmbin,wpmbins);
	   hist_tbmass_md_DeepAK8_JES_dn[it][ib][isb][ijes]->Sumw2();	  
	   
	   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_JES%i_up2",it,ib,isb,ijes+1);
	   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i JES %i (Up)",it+1,ib+1,isb+1,ijes+1);
	   hist_tbmass_md_DeepAK8_JES_up_2[it][ib][isb][ijes] = new TH1D(name,title,nvwpmbin,wpmbins);
	   hist_tbmass_md_DeepAK8_JES_up_2[it][ib][isb][ijes]->Sumw2();
	  
	   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_JES%i_dn2",it,ib,isb,ijes+1);
	   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i JES %i (Down)",it+1,ib+1,isb+1,ijes+1);
	   hist_tbmass_md_DeepAK8_JES_dn_2[it][ib][isb][ijes] = new TH1D(name,title,nvwpmbin,wpmbins);
	   hist_tbmass_md_DeepAK8_JES_dn_2[it][ib][isb][ijes]->Sumw2();	  
	}
   
   sprintf(name,"tb_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Pt of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbpt_md_DeepAK8[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_tbpt_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_y_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Rapidity of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbrap_md_DeepAK8[it][ib][isb] = new TH1D(name,title,50,-5.,5.);
   hist_tbrap_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BCandidate_Jet_Parton_Flavour_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"BCandidate Jet Parton Flavour | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_partonflav_AK4_md_DeepAK8[it][ib][isb] = new TH1D(name,title,22,-0.01,21.99);
   hist_partonflav_AK4_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet P_{T} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetpt_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_topjetpt_sel_md_DeepAK8[it][ib][isb]->Sumw2();
    
   for(int  imsd=0; imsd<ntopsdmassbins; imsd++){
   sprintf(name,"TopJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_msdbin%i",it,ib,isb,imsd);
   sprintf(title,"TopJet P_{T} Top WP %i B WP %i SubB WP %i M_{SD} bin %i",it+1,ib+1,isb+1,imsd+1); 
   hist_topjetpt_msdbin_sel_md_DeepAK8[it][ib][isb][imsd] = new TH1D(name,title,noptbins,ptbins);
   hist_topjetpt_msdbin_sel_md_DeepAK8[it][ib][isb][imsd]->Sumw2();
	}
    
   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   for(int ipdf=0; ipdf<npdfmax; ipdf++){
	   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_PDF%i",it,ib,isb,ipdf);
	   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i PDF %i",it+1,ib+1,isb+1,ipdf+1);
	   hist_topjetsdmass_sel_md_DeepAK8_PDF[it][ib][isb][ipdf] = new TH1D(name,title,nsdmassbins,sdmassbins);
	   hist_topjetsdmass_sel_md_DeepAK8_PDF[it][ib][isb][ipdf]->Sumw2();
	}
	   
   for(int iscale=0; iscale<nscalemax; iscale++){
	   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_Scale%i",it,ib,isb,iscale);
	   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i Scale %i",it+1,ib+1,isb+1,iscale+1);
	   hist_topjetsdmass_sel_md_DeepAK8_Scale[it][ib][isb][iscale] = new TH1D(name,title,nsdmassbins,sdmassbins);
	   hist_topjetsdmass_sel_md_DeepAK8_Scale[it][ib][isb][iscale]->Sumw2();
   }
  
   for(int ips=0; ips<npsmax; ips++){
	   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_PS%i",it,ib,isb,ips);
	   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i PS %i",it+1,ib+1,isb+1,ips+1);
	   hist_topjetsdmass_sel_md_DeepAK8_PS[it][ib][isb][ips] = new TH1D(name,title,nsdmassbins,sdmassbins);
	   hist_topjetsdmass_sel_md_DeepAK8_PS[it][ib][isb][ips]->Sumw2();
   }
   
   for(int ibunc=0; ibunc<nbtagmax; ibunc++){
	   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_BTag%i_up",it,ib,isb,ibunc);
	   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i BTag %i (Up)",it+1,ib+1,isb+1,ibunc+1);
	   hist_topjetsdmass_sel_md_DeepAK8_Btag_up[it][ib][isb][ibunc] = new TH1D(name,title,nsdmassbins,sdmassbins);
	   hist_topjetsdmass_sel_md_DeepAK8_Btag_up[it][ib][isb][ibunc]->Sumw2();
	   
	   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_BTag%i_dn",it,ib,isb,ibunc);
	   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i BTag %i (Down)",it+1,ib+1,isb+1,ibunc+1);
	   hist_topjetsdmass_sel_md_DeepAK8_Btag_dn[it][ib][isb][ibunc] = new TH1D(name,title,nsdmassbins,sdmassbins);
	   hist_topjetsdmass_sel_md_DeepAK8_Btag_dn[it][ib][isb][ibunc]->Sumw2();
   }
   
   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_JER_up",it,ib,isb);
   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i JER (Up)",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel_md_DeepAK8_JER_up[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel_md_DeepAK8_JER_up[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_JER_dn",it,ib,isb);
   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i JER (Down)",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel_md_DeepAK8_JER_dn[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel_md_DeepAK8_JER_dn[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_JMS_up",it,ib,isb);
   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i JMS (Up)",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel_md_DeepAK8_JMS_up[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel_md_DeepAK8_JMS_up[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_JMS_dn",it,ib,isb);
   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i JMS (Down)",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel_md_DeepAK8_JMS_dn[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel_md_DeepAK8_JMS_dn[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_JMR_up",it,ib,isb);
   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i JMR (Up)",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel_md_DeepAK8_JMR_up[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel_md_DeepAK8_JMR_up[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_JMR_dn",it,ib,isb);
   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i JMR (Down)",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel_md_DeepAK8_JMR_dn[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel_md_DeepAK8_JMR_dn[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_DeepAK8SFup",it,ib,isb);
   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel_md_DeepAK8up[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel_md_DeepAK8up[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_DeepAK8SFdn",it,ib,isb);
   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel_md_DeepAK8dn[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel_md_DeepAK8dn[it][ib][isb]->Sumw2();
   
   for(int ijes=0; ijes<njesmax; ijes++){
	   
	   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_JES%i_up",it,ib,isb,ijes+1);
	   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i JES %i (Up)",it+1,ib+1,isb+1,ijes+1);
	   hist_topjetsdmass_sel_md_DeepAK8_JES_up[it][ib][isb][ijes] = new TH1D(name,title,nsdmassbins,sdmassbins);
	   hist_topjetsdmass_sel_md_DeepAK8_JES_up[it][ib][isb][ijes]->Sumw2();
   
	   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_JES%i_dn",it,ib,isb,ijes+1);
	   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i JES %i (Down)",it+1,ib+1,isb+1,ijes+1);
	   hist_topjetsdmass_sel_md_DeepAK8_JES_dn[it][ib][isb][ijes] = new TH1D(name,title,nsdmassbins,sdmassbins);
	   hist_topjetsdmass_sel_md_DeepAK8_JES_dn[it][ib][isb][ijes]->Sumw2();
	   
	   sprintf(name,"TopJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_JES%i_up",it,ib,isb,ijes+1);
	   sprintf(title,"TopJet p_{T} Top WP %i B WP %i SubB WP %i JES %i (Up)",it+1,ib+1,isb+1,ijes+1);
	   hist_topjetpt_sel_md_DeepAK8_JES_up[it][ib][isb][ijes] = new TH1D(name,title,noptbins,ptbins);
	   hist_topjetpt_sel_md_DeepAK8_JES_up[it][ib][isb][ijes]->Sumw2();
   
	   sprintf(name,"TopJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_JES%i_dn",it,ib,isb,ijes+1);
	   sprintf(title,"TopJet p_{T} Top WP %i B WP %i SubB WP %i JES %i (Down)",it+1,ib+1,isb+1,ijes+1);
	   hist_topjetpt_sel_md_DeepAK8_JES_dn[it][ib][isb][ijes] = new TH1D(name,title,noptbins,ptbins);
	   hist_topjetpt_sel_md_DeepAK8_JES_dn[it][ib][isb][ijes]->Sumw2();
	   
	}	   
   
   sprintf(name,"TopJet_DeepTagTvsQCD_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetdeeptopscore_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeeptopscore_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_MD_DeepTagTvsQCD_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet Mass-Decorrelated DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetdeepmdtopscore_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeepmdtopscore_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_MD_DeepTagTvsQCD_DeepAK8MD_toptag%i_btag%i_subbtag%i_DeepAK8SFup",it,ib,isb);
   sprintf(title,"TopJet Mass-Decorrelated DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetdeepmdtopscore_sel_md_DeepAK8up[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeepmdtopscore_sel_md_DeepAK8up[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_MD_DeepTagTvsQCD_DeepAK8MD_toptag%i_btag%i_subbtag%i_DeepAK8SFdn",it,ib,isb);
   sprintf(title,"TopJet Mass-Decorrelated DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetdeepmdtopscore_sel_md_DeepAK8dn[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeepmdtopscore_sel_md_DeepAK8dn[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet SubJet CSVv2 Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsbtag_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,20,0,1.0);
   hist_topjetsbtag_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepCSV_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet SubJet DeepCSV Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdeepCSV_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,20,0,1.0);
   hist_topjetsdeepCSV_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Tau32_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet #tau_{32} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjettau32_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,50,0,1.0);
   hist_topjettau32_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_PartonFlavor_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet parton flavor Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_partonflav_AK8_md_DeepAK8[it][ib][isb] = new TH1D(name,title,22,-0.01,21.99);
   hist_partonflav_AK8_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"BJet P_{T} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetpt_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_bjetpt_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   for(int ieta=0; ieta<netabins; ieta++){
   
	sprintf(name,"BJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_Eta%i",it,ib,isb,ieta+1);
	sprintf(title,"BJet P_{T} Top WP %i B WP %i SubB WP %i Eta %i",it+1,ib+1,isb+1,ieta+1);
	hist_bjetpt_sel_md_DeepAK8_beta[it][ib][isb][ieta] = new TH1D(name,title,noptbins,ptbins);
	hist_bjetpt_sel_md_DeepAK8_beta[it][ib][isb][ieta]->Sumw2();
    
    sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_Eta%i_1",it,ib,isb,ieta+1);
    sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
    hist_tbmass_md_DeepAK8_eta1[it][ib][isb][ieta] = new TH1D(name,title,nvwpmbin,wpmbins);
    hist_tbmass_md_DeepAK8_eta1[it][ib][isb][ieta]->Sumw2();
    
    sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_Eta%i_2",it,ib,isb,ieta+1);
    sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
    hist_tbmass_md_DeepAK8_eta2[it][ib][isb][ieta] = new TH1D(name,title,nvwpmbin,wpmbins);
    hist_tbmass_md_DeepAK8_eta2[it][ib][isb][ieta]->Sumw2();
    
	}
    
    for(int ibunc=0; ibunc<nbtagmax; ibunc++){
	   sprintf(name,"BJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_BTag%i_up",it,ib,isb,ibunc);
	   sprintf(title,"BJet P_{T} Top WP %i B WP %i SubB WP %i BTag %i (Up)",it+1,ib+1,isb+1,ibunc+1);
	   hist_bjetpt_sel_md_DeepAK8_Btag_up[it][ib][isb][ibunc] = new TH1D(name,title,noptbins,ptbins);
	   hist_bjetpt_sel_md_DeepAK8_Btag_up[it][ib][isb][ibunc]->Sumw2();
	   
	   sprintf(name,"BJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_BTag%i_dn",it,ib,isb,ibunc);
	   sprintf(title,"BJet P_{T} Top WP %i B WP %i SubB WP %i BTag %i (Down)",it+1,ib+1,isb+1,ibunc+1);
	   hist_bjetpt_sel_md_DeepAK8_Btag_dn[it][ib][isb][ibunc] = new TH1D(name,title,noptbins,ptbins);
	   hist_bjetpt_sel_md_DeepAK8_Btag_dn[it][ib][isb][ibunc]->Sumw2();
   }
    
   for(int ijes=0; ijes<njesmax; ijes++){
	   sprintf(name,"BJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_JES%i_up",it,ib,isb,ijes+1);
	   sprintf(title,"BJet P_{T} Top WP %i B WP %i SubB WP %i JES %i (Up)",it+1,ib+1,isb+1,ijes+1);
	   hist_bjetpt_sel_md_DeepAK8_JES_up[it][ib][isb][ijes] = new TH1D(name,title,noptbins,ptbins);
	   hist_bjetpt_sel_md_DeepAK8_JES_up[it][ib][isb][ijes]->Sumw2();
	   
	   sprintf(name,"BJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_JES%i_dn",it,ib,isb,ijes+1);
	   sprintf(title,"BJet P_{T} Top WP %i B WP %i SubB WP %i JES %i (Down)",it+1,ib+1,isb+1,ijes+1);
	   hist_bjetpt_sel_md_DeepAK8_JES_dn[it][ib][isb][ijes] = new TH1D(name,title,noptbins,ptbins);
	   hist_bjetpt_sel_md_DeepAK8_JES_dn[it][ib][isb][ijes]->Sumw2();
   } 
    
   sprintf(name,"BJet_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"BJet Mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmass_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_bjetmass_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_BTag_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"BJet CSVv2 Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetbtag_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,20,0,1.0);
   hist_bjetbtag_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_TopTagScore_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Top-tag Score of AK8 jet corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Topscore_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,100,0,1.0);
   hist_bjetmatch_AK8Topscore_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_WTagScore_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"W-tag Score of AK8 jet corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Wscore_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,100,0,1.0);
   hist_bjetmatch_AK8Wscore_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Soft drop mass of AK8 jet corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_sdmass_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_bjetmatch_sdmass_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"Event_HT_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"H_{T} | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_ht_md_DeepAK8[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_ht_md_DeepAK8[it][ib][isb]->Sumw2();
   
   // isFake OR high ht
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"Pt of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbpt_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_tbpt_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_y_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"Rapidity of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbrap_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,50,-5.,5.);
   hist_tbrap_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BCandidate_Jet_Parton_Flavour_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"BCandidate Jet Parton Flavour | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_partonflav_AK4_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,22,-0.01,21.99);
   hist_partonflav_AK4_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"TopJet P_{T} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetpt_sel_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_topjetpt_sel_md_DeepAK8_3[it][ib][isb]->Sumw2();
    
   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_DeepTagTvsQCD_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"TopJet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetdeeptopscore_sel_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeeptopscore_sel_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_MD_DeepTagTvsQCD_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"TopJet Mass-Decorrelated DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetdeepmdtopscore_sel_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeepmdtopscore_sel_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"TopJet SubJet CSVv2 Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsbtag_sel_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,20,0,1.0);
   hist_topjetsbtag_sel_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepCSV_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"TopJet SubJet DeepCSV Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdeepCSV_sel_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,20,0,1.0);
   hist_topjetsdeepCSV_sel_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Tau32_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"TopJet #tau_{32} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjettau32_sel_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,50,0,1.0);
   hist_topjettau32_sel_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"BJet P_{T} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetpt_sel_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_bjetpt_sel_md_DeepAK8_3[it][ib][isb]->Sumw2();
    
   sprintf(name,"BJet_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"BJet Mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmass_sel_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_bjetmass_sel_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_BTag_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"BJet CSVv2 Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetbtag_sel_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,20,0,1.0);
   hist_bjetbtag_sel_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_TopTagScore_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"Top-tag Score of AK8 jet corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Topscore_sel_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,100,0,1.0);
   hist_bjetmatch_AK8Topscore_sel_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_WTagScore_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"W-tag Score of AK8 jet corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Wscore_sel_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,100,0,1.0);
   hist_bjetmatch_AK8Wscore_sel_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"Soft drop mass of AK8 jet corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_sdmass_sel_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_bjetmatch_sdmass_sel_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"Event_HT_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"H_{T} | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_ht_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_ht_md_DeepAK8_3[it][ib][isb]->Sumw2();
 
   
   // 2D hists for mtb //
   
   sprintf(name,"Mtb_vs_TJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"M_{tb} vs TopJet P_{T} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_topjetpt_md_DeepAK8[it][ib][isb] = new TH2D(name,title,noptbins,ptbins,nvwpmbin,wpmbins);
   hist_tbmass_topjetpt_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"Mtb_vs_BJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"M_{tb} vs BJet P_{T} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_bjetpt_md_DeepAK8[it][ib][isb] = new TH2D(name,title,noptbins,ptbins,nvwpmbin,wpmbins);
   hist_tbmass_bjetpt_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"Mtb_vs_TJet_Eta_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"M_{tb} vs TopJet #eta Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_topjeteta_md_DeepAK8[it][ib][isb] = new TH2D(name,title,30,-3.,3.,nvwpmbin,wpmbins);
   hist_tbmass_topjeteta_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"Mtb_vs_BJet_Eta_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"M_{tb} vs BJet #eta Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_bjeteta_md_DeepAK8[it][ib][isb] = new TH2D(name,title,30,-3.,3.,nvwpmbin,wpmbins);
   hist_tbmass_bjeteta_md_DeepAK8[it][ib][isb]->Sumw2();
   
   //2D hists for mtb ends //
   
   // MD for t-tbar //
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_DeepAK8pass",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_tt_deepak8pass[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_tt_deepak8pass[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"Pt of tb Pair | Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_tbpt_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_tbpt_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_y_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"Rapidity of tb Pair | Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_tbrap_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,50,-5.,5.);
   hist_tbrap_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   for(int  imsd=0; imsd<ntopsdmassbins; imsd++){
	sprintf(name,"TopJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_msdbin%i_ttbar",it,ib,isb,imsd);
	sprintf(title,"TopJet P_{T} Top WP %i B WP %i SubB WP %i M_{SD} bin %i",it+1,ib+1,isb+1,imsd+1); 
	hist_topjetpt_msdbin_sel_md_DeepAK8_tt[it][ib][isb][imsd] = new TH1D(name,title,noptbins,ptbins);
	hist_topjetpt_msdbin_sel_md_DeepAK8_tt[it][ib][isb][imsd]->Sumw2();
	}
   
   sprintf(name,"BCandidate_Jet_Parton_Flavour_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"BCandidate Jet Parton Flavour | Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_partonflav_AK4_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,22,-0.01,21.99);
   hist_partonflav_AK4_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"TopJet P_{T} Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjetpt_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_topjetpt_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_DeepAK8pass",it,ib,isb);
   sprintf(title,"TopJet P_{T} Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjetpt_sel_md_DeepAK8_tt_deepak8pass[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_topjetpt_sel_md_DeepAK8_tt_deepak8pass[it][ib][isb]->Sumw2();
    
   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_DeepAK8pass",it,ib,isb);
   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel_md_DeepAK8_tt_deepak8pass[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel_md_DeepAK8_tt_deepak8pass[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_DeepTagTvsQCD_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"TopJet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjetdeeptopscore_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeeptopscore_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_MD_DeepTagTvsQCD_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"TopJet Mass-Decorrelated DeepTagTvsQCD Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjetdeepmdtopscore_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeepmdtopscore_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"TopJet SubJet CSVv2 Score Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjetsbtag_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,20,0,1.0);
   hist_topjetsbtag_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepCSV_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"TopJet SubJet DeepCSV Score Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjetsdeepCSV_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,20,0,1.0);
   hist_topjetsdeepCSV_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Tau32_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"TopJet #tau_{32} Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjettau32_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,50,0,1.0);
   hist_topjettau32_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"BJet P_{T} Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_bjetpt_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_bjetpt_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
    
   sprintf(name,"BJet_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"BJet Mass Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_bjetmass_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_bjetmass_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_BTag_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"BJet CSVv2 Score Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_bjetbtag_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,20,0,1.0);
   hist_bjetbtag_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_TopTagScore_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"Top-tag Score of AK8 jet corresponding to Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Topscore_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,100,0,1.0);
   hist_bjetmatch_AK8Topscore_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_WTagScore_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"W-tag Score of AK8 jet corresponding to Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Wscore_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,100,0,1.0);
   hist_bjetmatch_AK8Wscore_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"Soft drop mass of AK8 jet corresponding to Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_bjetmatch_sdmass_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_bjetmatch_sdmass_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"Event_HT_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"H_{T} | Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_ht_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_ht_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   
   // isFake OR high ht
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"Pt of tb Pair | Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_tbpt_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_tbpt_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_y_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"Rapidity of tb Pair | Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_tbrap_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,50,-5.,5.);
   hist_tbrap_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BCandidate_Jet_Parton_Flavour_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"BCandidate Jet Parton Flavour | Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_partonflav_AK4_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,22,-0.01,21.99);
   hist_partonflav_AK4_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"TopJet P_{T} Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjetpt_sel_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_topjetpt_sel_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
    
   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"TopJet Soft drop mass Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_DeepTagTvsQCD_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"TopJet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjetdeeptopscore_sel_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeeptopscore_sel_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_MD_DeepTagTvsQCD_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"TopJet Mass-Decorrelated DeepTagTvsQCD Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjetdeepmdtopscore_sel_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeepmdtopscore_sel_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"TopJet SubJet CSVv2 Score Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjetsbtag_sel_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,20,0,1.0);
   hist_topjetsbtag_sel_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepCSV_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"TopJet SubJet DeepCSV Score Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjetsdeepCSV_sel_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,20,0,1.0);
   hist_topjetsdeepCSV_sel_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Tau32_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"TopJet #tau_{32} Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjettau32_sel_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,50,0,1.0);
   hist_topjettau32_sel_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"BJet P_{T} Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_bjetpt_sel_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_bjetpt_sel_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
    
   sprintf(name,"BJet_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"BJet Mass Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_bjetmass_sel_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_bjetmass_sel_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_BTag_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"BJet CSVv2 Score Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_bjetbtag_sel_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,20,0,1.0);
   hist_bjetbtag_sel_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_TopTagScore_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"Top-tag Score of AK8 jet corresponding to Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Topscore_sel_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,100,0,1.0);
   hist_bjetmatch_AK8Topscore_sel_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_WTagScore_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"W-tag Score of AK8 jet corresponding to Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Wscore_sel_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,100,0,1.0);
   hist_bjetmatch_AK8Wscore_sel_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"Soft drop mass of AK8 jet corresponding to Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_bjetmatch_sdmass_sel_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_bjetmatch_sdmass_sel_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"Event_HT_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"H_{T} | Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_ht_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_ht_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   // MD for t-tbar ends
   
   }
  }
 }

 // trigger histograms //
 
 for(int imsd = 0; imsd < trig_msdbins; imsd++){
 
   sprintf(name,"HT_SingleMu_Passed_Msd%i",imsd+1);
   sprintf(title,"HT | SingleMuon trigger passed");
   hist_ht_mutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_ht_mutrig[imsd]->Sumw2();
   
   sprintf(name,"LeadJetPt_SingleMu_Passed_Msd%i",imsd+1);
   sprintf(title,"Leading jet p_{T} | SingleMuon trigger passed");
   hist_pt_mutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_pt_mutrig[imsd]->Sumw2();
   
   sprintf(name,"HT_SingleMuElec_Passed_Msd%i",imsd+1);
   sprintf(title,"HT | (SingleMuon||SingleElectron) trigger passed");
   hist_ht_emutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_ht_emutrig[imsd]->Sumw2();
   
   sprintf(name,"LeadJetPt_SingleMuElec_Passed_Msd%i",imsd+1);
   sprintf(title,"Leading jet p_{T} | (SingleMuon||SingleElectron) trigger passed");
   hist_pt_emutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_pt_emutrig[imsd]->Sumw2();
   
   sprintf(name,"LeadBJetPt_SingleMuElec_Passed_Msd%i",imsd+1);
   sprintf(title,"B jet p_{T} | (SingleMuon||SingleElectron) trigger passed");
   hist_bpt_emutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_bpt_emutrig[imsd]->Sumw2();
   
   sprintf(name,"TwoJetPtSum_SingleMuElec_Passed_Msd%i",imsd+1);
   sprintf(title,"(P_{T}^{0} + P_{T}^{1})| (SingleMuon||SingleElectron) trigger passed");
   hist_ptsum_emutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_ptsum_emutrig[imsd]->Sumw2();
   
   sprintf(name,"Mtb_SingleMuElec_Passed_Msd%i",imsd+1);
   sprintf(title,"M_{tb}| (SingleMuon||SingleElectron) trigger passed");
   hist_mtb_emutrig[imsd]= new TH1D(name,title,nvwpmbin_trig,wpmbins_trig);
   hist_mtb_emutrig[imsd]->Sumw2();
 
 /////////////// Ref : Electron + Muon Trigger //////////
   
   // ht //
   
   sprintf(name,"HT_SingleMuorElec_HT1050_Passed_Msd%i",imsd+1);
   sprintf(title,"HT |  (SingleMuon || Single Electron) && PFHT1050 trigger passed");
   hist_ht_HT1050_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_ht_HT1050_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"HT_SingleMuorElec_AK4Pt500_Passed_Msd%i",imsd+1);
   sprintf(title,"HT |  (SingleMuon || Single Electron) && AK4PFJet500 trigger passed");
   hist_ht_AK4Pt500_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_ht_AK4Pt500_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"HT_SingleMuorElec_AK8Pt500_Passed_Msd%i",imsd+1);
   sprintf(title,"HT |  (SingleMuon || Single Electron) && AK8Pt500 trigger passed");
   hist_ht_AK8Pt500_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_ht_AK8Pt500_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"HT_SingleMuorElec_AK8PFJet420_TrimMass30_Passed_Msd%i",imsd+1);
   sprintf(title,"HT |  (SingleMuon || Single Electron) && AK8PFJet420_TrimMass30 trigger passed");
   hist_ht_AK8PFJet420_TrimMass30_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_ht_AK8PFJet420_TrimMass30_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"HT_SingleMuorElec_AK8HT900_TrimMass50_Passed_Msd%i",imsd+1);
   sprintf(title,"HT |  (SingleMuon || Single Electron) && AK8PFHT900_TrimMass50 trigger passed");
   hist_ht_AK8PFHT900_TrimMass50_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_ht_AK8PFHT900_TrimMass50_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"HT_SingleMuorElec_AK8DiPFJet300_200_TrimMass30_Passed_Msd%i",imsd+1);
   sprintf(title,"HT |  (SingleMuon || Single Electron) && AK8DiPFJet300_200_TrimMass30 trigger passed");
   hist_ht_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_ht_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"HT_SingleMuorElec_AnyTrig_Passed_Msd%i",imsd+1);
   sprintf(title,"HT |  (SingleMuon || Single Electron) && any trigger passed");
   hist_ht_all_wemutrig[imsd] = new TH1D(name,title,nohtbins,htbins);
   hist_ht_all_wemutrig[imsd]->Sumw2(); 

   
   // pt //
   
   sprintf(name,"LeadJetPt_SingleMuorElec_HT1050_Passed_Msd%i",imsd+1);
   sprintf(title,"Leading jet p_{T} |  (SingleMuon || Single Electron) && PFHT1050 trigger passed");
   hist_pt_HT1050_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_pt_HT1050_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"LeadJetPt_SingleMuorElec_AK4Pt500_Passed_Msd%i",imsd+1);
   sprintf(title,"Leading jet p_{T} |  (SingleMuon || Single Electron) && AK4PFJet500 trigger passed");
   hist_pt_AK4Pt500_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_pt_AK4Pt500_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"LeadJetPt_SingleMuorElec_AK8Pt500_Passed_Msd%i",imsd+1);
   sprintf(title,"Leading jet p_{T} |  (SingleMuon || Single Electron) && AK8Pt500 trigger passed");
   hist_pt_AK8Pt500_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_pt_AK8Pt500_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"LeadJetPt_SingleMuorElec_AK8PFJet420_TrimMass30_Passed_Msd%i",imsd+1);
   sprintf(title,"Leading jet p_{T} |  (SingleMuon || Single Electron) && AK8PFJet420_TrimMass30 trigger passed");
   hist_pt_AK8PFJet420_TrimMass30_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_pt_AK8PFJet420_TrimMass30_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"LeadJetPt_SingleMuorElec_AK8HT900_TrimMass50_Passed_Msd%i",imsd+1);
   sprintf(title,"Leading jet p_{T} |  (SingleMuon || Single Electron) && AK8PFHT900_TrimMass50 trigger passed");
   hist_pt_AK8PFHT900_TrimMass50_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_pt_AK8PFHT900_TrimMass50_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"LeadJetPt_SingleMuorElec_AK8DiPFJet300_200_TrimMass30_Passed_Msd%i",imsd+1);
   sprintf(title,"Leading jet p_{T} |  (SingleMuon || Single Electron) && AK8DiPFJet300_200_TrimMass30 trigger passed");
   hist_pt_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_pt_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"LeadJetPt_SingleMuorElec_AnyTrig_Passed_Msd%i",imsd+1);
   sprintf(title,"Leading jet p_{T} |  (SingleMuon || Single Electron) && any trigger passed");
   hist_pt_all_wemutrig[imsd] = new TH1D(name,title,nohtbins,htbins);
   hist_pt_all_wemutrig[imsd]->Sumw2(); 
   
   // b pt //
   
   sprintf(name,"LeadBJetPt_SingleMuorElec_HT1050_Passed_Msd%i",imsd+1);
   sprintf(title,"B jet p_{T} |  (SingleMuon || Single Electron) && PFHT1050 trigger passed");
   hist_bpt_HT1050_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_bpt_HT1050_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"LeadBJetPt_SingleMuorElec_AK4Pt500_Passed_Msd%i",imsd+1);
   sprintf(title,"B jet p_{T} |  (SingleMuon || Single Electron) && AK4PFJet500 trigger passed");
   hist_bpt_AK4Pt500_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_bpt_AK4Pt500_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"LeadBJetPt_SingleMuorElec_AK8Pt500_Passed_Msd%i",imsd+1);
   sprintf(title,"B jet p_{T} |  (SingleMuon || Single Electron) && AK8Pt500 trigger passed");
   hist_bpt_AK8Pt500_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_bpt_AK8Pt500_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"LeadBJetPt_SingleMuorElec_AK8PFJet420_TrimMass30_Passed_Msd%i",imsd+1);
   sprintf(title,"B jet p_{T} |  (SingleMuon || Single Electron) && AK8PFJet420_TrimMass30 trigger passed");
   hist_bpt_AK8PFJet420_TrimMass30_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_bpt_AK8PFJet420_TrimMass30_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"LeadBJetPt_SingleMuorElec_AK8HT900_TrimMass50_Passed_Msd%i",imsd+1);
   sprintf(title,"B jet p_{T} |  (SingleMuon || Single Electron) && AK8PFHT900_TrimMass50 trigger passed");
   hist_bpt_AK8PFHT900_TrimMass50_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_bpt_AK8PFHT900_TrimMass50_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"LeadBJetPt_SingleMuorElec_AK8DiPFJet300_200_TrimMass30_Passed_Msd%i",imsd+1);
   sprintf(title,"B jet p_{T} |  (SingleMuon || Single Electron) && AK8DiPFJet300_200_TrimMass30 trigger passed");
   hist_bpt_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_bpt_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"LeadBJetPt_SingleMuorElec_AnyTrig_Passed_Msd%i",imsd+1);
   sprintf(title,"B jet p_{T} |  (SingleMuon || Single Electron) && any trigger passed");
   hist_bpt_all_wemutrig[imsd] = new TH1D(name,title,nohtbins,htbins);
   hist_bpt_all_wemutrig[imsd]->Sumw2(); 
   
   
   // ptsum //
   
   sprintf(name,"TwoJetPtSum_SingleMuorElec_HT1050_Passed_Msd%i",imsd+1);
   sprintf(title,"(P_{T}^{0} + P_{T}^{1}) |  (SingleMuon || Single Electron) && PFHT1050 trigger passed");
   hist_ptsum_HT1050_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_ptsum_HT1050_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"TwoJetPtSum_SingleMuorElec_AK4Pt500_Passed_Msd%i",imsd+1);
   sprintf(title,"(P_{T}^{0} + P_{T}^{1}) |  (SingleMuon || Single Electron) && AK4PFJet500 trigger passed");
   hist_ptsum_AK4Pt500_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_ptsum_AK4Pt500_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"TwoJetPtSum_SingleMuorElec_AK8Pt500_Passed_Msd%i",imsd+1);
   sprintf(title,"(P_{T}^{0} + P_{T}^{1}) |  (SingleMuon || Single Electron) && AK8Pt500 trigger passed");
   hist_ptsum_AK8Pt500_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_ptsum_AK8Pt500_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"TwoJetPtSum_SingleMuorElec_AK8PFJet420_TrimMass30_Passed_Msd%i",imsd+1);
   sprintf(title,"(P_{T}^{0} + P_{T}^{1}) |  (SingleMuon || Single Electron) && AK8PFJet420_TrimMass30 trigger passed");
   hist_ptsum_AK8PFJet420_TrimMass30_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_ptsum_AK8PFJet420_TrimMass30_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"TwoJetPtSum_SingleMuorElec_AK8HT900_TrimMass50_Passed_Msd%i",imsd+1);
   sprintf(title,"(P_{T}^{0} + P_{T}^{1}) |  (SingleMuon || Single Electron) && AK8PFHT900_TrimMass50 trigger passed");
   hist_ptsum_AK8PFHT900_TrimMass50_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_ptsum_AK8PFHT900_TrimMass50_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"TwoJetPtSum_SingleMuorElec_AK8DiPFJet300_200_TrimMass30_Passed_Msd%i",imsd+1);
   sprintf(title,"(P_{T}^{0} + P_{T}^{1}) |  (SingleMuon || Single Electron) && AK8DiPFJet300_200_TrimMass30 trigger passed");
   hist_ptsum_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]= new TH1D(name,title,nohtbins,htbins);
   hist_ptsum_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"TwoJetPtSum_SingleMuorElec_AnyTrig_Passed_Msd%i",imsd+1);
   sprintf(title,"(P_{T}^{0} + P_{T}^{1}) |  (SingleMuon || Single Electron) && any trigger passed");
   hist_ptsum_all_wemutrig[imsd] = new TH1D(name,title,nohtbins,htbins);
   hist_ptsum_all_wemutrig[imsd]->Sumw2();
   
   // mtb //
   
   sprintf(name,"Mtb_SingleMuorElec_HT1050_Passed_Msd%i",imsd+1);
   sprintf(title,"M_{tb} |  (SingleMuon || Single Electron) && PFHT1050 trigger passed");
   hist_mtb_HT1050_wemutrig[imsd]= new TH1D(name,title,nvwpmbin_trig,wpmbins_trig);
   hist_mtb_HT1050_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"Mtb_SingleMuorElec_AK4Pt500_Passed_Msd%i",imsd+1);
   sprintf(title,"M_{tb} |  (SingleMuon || Single Electron) && AK4PFJet500 trigger passed");
   hist_mtb_AK4Pt500_wemutrig[imsd]= new TH1D(name,title,nvwpmbin_trig,wpmbins_trig);
   hist_mtb_AK4Pt500_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"Mtb_SingleMuorElec_AK8Pt500_Passed_Msd%i",imsd+1);
   sprintf(title,"M_{tb} |  (SingleMuon || Single Electron) && AK8Pt500 trigger passed");
   hist_mtb_AK8Pt500_wemutrig[imsd]= new TH1D(name,title,nvwpmbin_trig,wpmbins_trig);
   hist_mtb_AK8Pt500_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"Mtb_SingleMuorElec_AK8PFJet420_TrimMass30_Passed_Msd%i",imsd+1);
   sprintf(title,"M_{tb} |  (SingleMuon || Single Electron) && AK8PFJet420_TrimMass30 trigger passed");
   hist_mtb_AK8PFJet420_TrimMass30_wemutrig[imsd]= new TH1D(name,title,nvwpmbin_trig,wpmbins_trig);
   hist_mtb_AK8PFJet420_TrimMass30_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"Mtb_SingleMuorElec_AK8HT900_TrimMass50_Passed_Msd%i",imsd+1);
   sprintf(title,"M_{tb} |  (SingleMuon || Single Electron) && AK8PFHT900_TrimMass50 trigger passed");
   hist_mtb_AK8PFHT900_TrimMass50_wemutrig[imsd]= new TH1D(name,title,nvwpmbin_trig,wpmbins_trig);
   hist_mtb_AK8PFHT900_TrimMass50_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"Mtb_SingleMuorElec_AK8DiPFJet300_200_TrimMass30_Passed_Msd%i",imsd+1);
   sprintf(title,"M_{tb} |  (SingleMuon || Single Electron) && AK8DiPFJet300_200_TrimMass30 trigger passed");
   hist_mtb_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]= new TH1D(name,title,nvwpmbin_trig,wpmbins_trig);
   hist_mtb_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]->Sumw2();
   
   sprintf(name,"Mtb_SingleMuorElec_AnyTrig_Passed_Msd%i",imsd+1);
   sprintf(title,"M_{tb} |  (SingleMuon || Single Electron) && any trigger passed");
   hist_mtb_all_wemutrig[imsd] = new TH1D(name,title,nvwpmbin_trig,wpmbins_trig);
   hist_mtb_all_wemutrig[imsd]->Sumw2();  
 
 }
 // trigger histo end //
 
   cout<<"Before read\n";
 
   ifstream file_db;
   file_db.open(rootfiles);

   while(!(file_db.eof())){
   file_db >> datafile;
   cout <<"datafile name is "<<datafile/*<<" with weight "<<weight<<" "<<weight2*/<<endl;
  if (strstr(datafile,"#")) continue;

   if(file_db.eof()) break;

   sprintf(infile, "%s", datafile);
   cout<<"infile "<<infile<<endl;
 //  TFile* fileIn = new TFile(infile, "read");
   TFile* fileIn = TFile::Open(infile);

   if ( fileIn->IsZombie() ) continue;
   if ( fileIn->Recover() == 0 ) continue;
   
   #ifdef Anal_2016
   
	btagvalue = 0.8484; //0.9535; //tight 0.8484; // medium //csvv2 
	btagvalue_deepCSV = 0.6321; //0.8953  ;//tight  0.6321; //medium // deepcsv 
//   btagvalue_deepFlavB = 0.45; //0.7221; // 0.3093;  // medium	//deepflav
	#ifdef B_Tight
	btagvalue_deepFlavB = 0.70221;
	#else
	btagvalue_deepFlavB = 0.45;
	#endif
   
	deepak8_cut = 0.929;    // 2016 0.5% mistag rate
	deepak8_cut_md = 0.632; // 2016 0.5% mistag rate
   
	for(int ipt=0; ipt<noAK8ptbins; ipt++){
		DeepAK8_SF[ipt] = DeepAK8_SF_16[ipt];
		DeepAK8_MD_SF[ipt] = DeepAK8_MD_SF_16[ipt];
		DeepAK8_MD_SF_up[ipt] = DeepAK8_MD_SF_16_up[ipt];
		DeepAK8_MD_SF_dn[ipt] = DeepAK8_MD_SF_16_dn[ipt];
	}
	
	for(int imsd=0; imsd<ntopsdmassbins; imsd++){
		for(int ht=0; ht<nohtbins; ht++){
			if(isQCD){
				trig_weight_factors[imsd][ht] = trig_weight_factors_2016_QCD[imsd][ht];
			}
			if(isTTBar){
				trig_weight_factors[imsd][ht] = trig_weight_factors_2016_TT[imsd][ht];
			}
			if(isST){
				trig_weight_factors[imsd][ht] = trig_weight_factors_2016_ST[imsd][ht];
			}
		}
	}
	 
	#endif
	
	#ifdef Anal_2017
   
	btagvalue = 0.8838;//0.9693;//tight 0.8838; //csvv2 medium
	btagvalue_deepCSV =  0.4941;//0.65;//myvalue 0.8001;// tight 0.4941; // deepcsv medium
	btagvalue_deepFlavB = 0.6; //0.6; // my final value //0.3033; //0.55;// my value 0.7489; // tight 0.3033; //deepflav medium
   
	deepak8_cut = 0.884;    // 2017 0.5% mistag rate
	deepak8_cut_md = 0.554; // 2017 0.5% mistag rate
   
	for(int ipt=0; ipt<noAK8ptbins; ipt++){
		DeepAK8_SF[ipt] = DeepAK8_SF_17[ipt];
		DeepAK8_MD_SF[ipt] = DeepAK8_MD_SF_17[ipt];
		DeepAK8_MD_SF_up[ipt] = DeepAK8_MD_SF_17_up[ipt];
		DeepAK8_MD_SF_dn[ipt] = DeepAK8_MD_SF_17_dn[ipt];
	}
	   
    for(int imsd=0; imsd<ntopsdmassbins; imsd++){
		for(int ht=0; ht<nohtbins; ht++){
			if(isQCD){
				trig_weight_factors[imsd][ht] = trig_weight_factors_2017_QCD[imsd][ht];
			}
			if(isTTBar){
				trig_weight_factors[imsd][ht] = trig_weight_factors_2017_TT[imsd][ht];
			}
			if(isST){
				trig_weight_factors[imsd][ht] = trig_weight_factors_2017_ST[imsd][ht];
			}
		}
	}
   
   	
	#endif

	#ifdef Anal_2018
	btagvalue = 0.8838; 
	btagvalue_deepCSV = 0.4184; //0.7527  ;//tight  0.4184; //medium // deepcsv 
	btagvalue_deepFlavB = 0.6;  //0.7264; ;//tight   0.2770;  // medium	//deepflav
   
	deepak8_cut = 0.922;      // 2018 0.5% mistag rate
	deepak8_cut_md = 0.685;   // 2018 0.5% mistag rate
   
	for(int ipt=0; ipt<noAK8ptbins; ipt++){
		DeepAK8_SF[ipt] = DeepAK8_SF_18[ipt];
		DeepAK8_MD_SF[ipt] = DeepAK8_MD_SF_18[ipt];
		DeepAK8_MD_SF_up[ipt] = DeepAK8_MD_SF_18_up[ipt];
		DeepAK8_MD_SF_dn[ipt] = DeepAK8_MD_SF_18_dn[ipt];
	}
   
	for(int imsd=0; imsd<ntopsdmassbins; imsd++){
		for(int ht=0; ht<nohtbins; ht++){
			if(isQCD){
				trig_weight_factors[imsd][ht] = trig_weight_factors_2018_QCD[imsd][ht];
			}
			if(isTTBar){
				trig_weight_factors[imsd][ht] = trig_weight_factors_2018_TT[imsd][ht];
			}
			if(isST){
				trig_weight_factors[imsd][ht] = trig_weight_factors_2018_ST[imsd][ht];
			}
		}
	}
   
	#endif

	#ifdef LHAPDFX
	sprintf(name,"cteq66");
	LHAPDF::initPDFSet(name, LHAPDF::LHGRID, 0);
  	#endif

   TTree *T1;
   T1 = (TTree*)fileIn->Get("Events");
	
   T1->SetBranchAddress("run", &run, &b_run);
   T1->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   T1->SetBranchAddress("event", &event, &b_event);
   T1->SetBranchAddress("HTXS_Higgs_pt", &HTXS_Higgs_pt, &b_HTXS_Higgs_pt);
   T1->SetBranchAddress("HTXS_Higgs_y", &HTXS_Higgs_y, &b_HTXS_Higgs_y);
   T1->SetBranchAddress("HTXS_stage1_1_cat_pTjet25GeV", &HTXS_stage1_1_cat_pTjet25GeV, &b_HTXS_stage1_1_cat_pTjet25GeV);
   T1->SetBranchAddress("HTXS_stage1_1_cat_pTjet30GeV", &HTXS_stage1_1_cat_pTjet30GeV, &b_HTXS_stage1_1_cat_pTjet30GeV);
   T1->SetBranchAddress("HTXS_stage1_1_fine_cat_pTjet25GeV", &HTXS_stage1_1_fine_cat_pTjet25GeV, &b_HTXS_stage1_1_fine_cat_pTjet25GeV);
   T1->SetBranchAddress("HTXS_stage1_1_fine_cat_pTjet30GeV", &HTXS_stage1_1_fine_cat_pTjet30GeV, &b_HTXS_stage1_1_fine_cat_pTjet30GeV);
   T1->SetBranchAddress("HTXS_stage_0", &HTXS_stage_0, &b_HTXS_stage_0);
   T1->SetBranchAddress("HTXS_stage_1_pTjet25", &HTXS_stage_1_pTjet25, &b_HTXS_stage_1_pTjet25);
   T1->SetBranchAddress("HTXS_stage_1_pTjet30", &HTXS_stage_1_pTjet30, &b_HTXS_stage_1_pTjet30);
   T1->SetBranchAddress("HTXS_njets25", &HTXS_njets25, &b_HTXS_njets25);
   T1->SetBranchAddress("HTXS_njets30", &HTXS_njets30, &b_HTXS_njets30);
   T1->SetBranchAddress("btagWeight_CSVV2", &btagWeight_CSVV2, &b_btagWeight_CSVV2);
   T1->SetBranchAddress("btagWeight_DeepCSVB", &btagWeight_DeepCSVB, &b_btagWeight_DeepCSVB);
   T1->SetBranchAddress("CaloMET_phi", &CaloMET_phi, &b_CaloMET_phi);
   T1->SetBranchAddress("CaloMET_pt", &CaloMET_pt, &b_CaloMET_pt);
   T1->SetBranchAddress("CaloMET_sumEt", &CaloMET_sumEt, &b_CaloMET_sumEt);
   T1->SetBranchAddress("ChsMET_phi", &ChsMET_phi, &b_ChsMET_phi);
   T1->SetBranchAddress("ChsMET_pt", &ChsMET_pt, &b_ChsMET_pt);
   T1->SetBranchAddress("ChsMET_sumEt", &ChsMET_sumEt, &b_ChsMET_sumEt);
   T1->SetBranchAddress("nCorrT1METJet", &nCorrT1METJet, &b_nCorrT1METJet);
   T1->SetBranchAddress("CorrT1METJet_area", CorrT1METJet_area, &b_CorrT1METJet_area);
   T1->SetBranchAddress("CorrT1METJet_eta", CorrT1METJet_eta, &b_CorrT1METJet_eta);
   T1->SetBranchAddress("CorrT1METJet_muonSubtrFactor", CorrT1METJet_muonSubtrFactor, &b_CorrT1METJet_muonSubtrFactor);
   T1->SetBranchAddress("CorrT1METJet_phi", CorrT1METJet_phi, &b_CorrT1METJet_phi);
   T1->SetBranchAddress("CorrT1METJet_rawPt", CorrT1METJet_rawPt, &b_CorrT1METJet_rawPt);
   T1->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
   T1->SetBranchAddress("Electron_deltaEtaSC", Electron_deltaEtaSC, &b_Electron_deltaEtaSC);
   T1->SetBranchAddress("Electron_dr03EcalRecHitSumEt", Electron_dr03EcalRecHitSumEt, &b_Electron_dr03EcalRecHitSumEt);
   T1->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEt", Electron_dr03HcalDepth1TowerSumEt, &b_Electron_dr03HcalDepth1TowerSumEt);
   T1->SetBranchAddress("Electron_dr03TkSumPt", Electron_dr03TkSumPt, &b_Electron_dr03TkSumPt);
   T1->SetBranchAddress("Electron_dr03TkSumPtHEEP", Electron_dr03TkSumPtHEEP, &b_Electron_dr03TkSumPtHEEP);
   T1->SetBranchAddress("Electron_dxy", Electron_dxy, &b_Electron_dxy);
   T1->SetBranchAddress("Electron_dxyErr", Electron_dxyErr, &b_Electron_dxyErr);
   T1->SetBranchAddress("Electron_dz", Electron_dz, &b_Electron_dz);
   T1->SetBranchAddress("Electron_dzErr", Electron_dzErr, &b_Electron_dzErr);
   T1->SetBranchAddress("Electron_eCorr", Electron_eCorr, &b_Electron_eCorr);
   T1->SetBranchAddress("Electron_eInvMinusPInv", Electron_eInvMinusPInv, &b_Electron_eInvMinusPInv);
   T1->SetBranchAddress("Electron_energyErr", Electron_energyErr, &b_Electron_energyErr);
   T1->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
   T1->SetBranchAddress("Electron_hoe", Electron_hoe, &b_Electron_hoe);
   T1->SetBranchAddress("Electron_ip3d", Electron_ip3d, &b_Electron_ip3d);
   T1->SetBranchAddress("Electron_jetPtRelv2", Electron_jetPtRelv2, &b_Electron_jetPtRelv2);
   T1->SetBranchAddress("Electron_jetRelIso", Electron_jetRelIso, &b_Electron_jetRelIso);
   T1->SetBranchAddress("Electron_mass", Electron_mass, &b_Electron_mass);
   T1->SetBranchAddress("Electron_miniPFRelIso_all", Electron_miniPFRelIso_all, &b_Electron_miniPFRelIso_all);
   T1->SetBranchAddress("Electron_miniPFRelIso_chg", Electron_miniPFRelIso_chg, &b_Electron_miniPFRelIso_chg);
   T1->SetBranchAddress("Electron_mvaFall17V1Iso", Electron_mvaFall17V1Iso, &b_Electron_mvaFall17V1Iso);
   T1->SetBranchAddress("Electron_mvaFall17V1noIso", Electron_mvaFall17V1noIso, &b_Electron_mvaFall17V1noIso);
   T1->SetBranchAddress("Electron_mvaFall17V2Iso", Electron_mvaFall17V2Iso, &b_Electron_mvaFall17V2Iso);
   T1->SetBranchAddress("Electron_mvaFall17V2noIso", Electron_mvaFall17V2noIso, &b_Electron_mvaFall17V2noIso);
   T1->SetBranchAddress("Electron_pfRelIso03_all", Electron_pfRelIso03_all, &b_Electron_pfRelIso03_all);
   T1->SetBranchAddress("Electron_pfRelIso03_chg", Electron_pfRelIso03_chg, &b_Electron_pfRelIso03_chg);
   T1->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
   T1->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
   T1->SetBranchAddress("Electron_r9", Electron_r9, &b_Electron_r9);
   T1->SetBranchAddress("Electron_sieie", Electron_sieie, &b_Electron_sieie);
   T1->SetBranchAddress("Electron_sip3d", Electron_sip3d, &b_Electron_sip3d);
   T1->SetBranchAddress("Electron_mvaTTH", Electron_mvaTTH, &b_Electron_mvaTTH);
   T1->SetBranchAddress("Electron_charge", Electron_charge, &b_Electron_charge);
   T1->SetBranchAddress("Electron_cutBased", Electron_cutBased, &b_Electron_cutBased);
   T1->SetBranchAddress("Electron_cutBased_Fall17_V1", Electron_cutBased_Fall17_V1, &b_Electron_cutBased_Fall17_V1);
   T1->SetBranchAddress("Electron_jetIdx", Electron_jetIdx, &b_Electron_jetIdx);
   T1->SetBranchAddress("Electron_pdgId", Electron_pdgId, &b_Electron_pdgId);
   T1->SetBranchAddress("Electron_photonIdx", Electron_photonIdx, &b_Electron_photonIdx);
   T1->SetBranchAddress("Electron_tightCharge", Electron_tightCharge, &b_Electron_tightCharge);
   T1->SetBranchAddress("Electron_vidNestedWPBitmap", Electron_vidNestedWPBitmap, &b_Electron_vidNestedWPBitmap);
   T1->SetBranchAddress("Electron_vidNestedWPBitmapHEEP", Electron_vidNestedWPBitmapHEEP, &b_Electron_vidNestedWPBitmapHEEP);
   T1->SetBranchAddress("Electron_convVeto", Electron_convVeto, &b_Electron_convVeto);
   T1->SetBranchAddress("Electron_cutBased_HEEP", Electron_cutBased_HEEP, &b_Electron_cutBased_HEEP);
   T1->SetBranchAddress("Electron_isPFcand", Electron_isPFcand, &b_Electron_isPFcand);
   T1->SetBranchAddress("Electron_lostHits", Electron_lostHits, &b_Electron_lostHits);
   T1->SetBranchAddress("Electron_mvaFall17V1Iso_WP80", Electron_mvaFall17V1Iso_WP80, &b_Electron_mvaFall17V1Iso_WP80);
   T1->SetBranchAddress("Electron_mvaFall17V1Iso_WP90", Electron_mvaFall17V1Iso_WP90, &b_Electron_mvaFall17V1Iso_WP90);
   T1->SetBranchAddress("Electron_mvaFall17V1Iso_WPL", Electron_mvaFall17V1Iso_WPL, &b_Electron_mvaFall17V1Iso_WPL);
   T1->SetBranchAddress("Electron_mvaFall17V1noIso_WP80", Electron_mvaFall17V1noIso_WP80, &b_Electron_mvaFall17V1noIso_WP80);
   T1->SetBranchAddress("Electron_mvaFall17V1noIso_WP90", Electron_mvaFall17V1noIso_WP90, &b_Electron_mvaFall17V1noIso_WP90);
   T1->SetBranchAddress("Electron_mvaFall17V1noIso_WPL", Electron_mvaFall17V1noIso_WPL, &b_Electron_mvaFall17V1noIso_WPL);
   T1->SetBranchAddress("Electron_mvaFall17V2Iso_WP80", Electron_mvaFall17V2Iso_WP80, &b_Electron_mvaFall17V2Iso_WP80);
   T1->SetBranchAddress("Electron_mvaFall17V2Iso_WP90", Electron_mvaFall17V2Iso_WP90, &b_Electron_mvaFall17V2Iso_WP90);
   T1->SetBranchAddress("Electron_mvaFall17V2Iso_WPL", Electron_mvaFall17V2Iso_WPL, &b_Electron_mvaFall17V2Iso_WPL);
   T1->SetBranchAddress("Electron_mvaFall17V2noIso_WP80", Electron_mvaFall17V2noIso_WP80, &b_Electron_mvaFall17V2noIso_WP80);
   T1->SetBranchAddress("Electron_mvaFall17V2noIso_WP90", Electron_mvaFall17V2noIso_WP90, &b_Electron_mvaFall17V2noIso_WP90);
   T1->SetBranchAddress("Electron_mvaFall17V2noIso_WPL", Electron_mvaFall17V2noIso_WPL, &b_Electron_mvaFall17V2noIso_WPL);
   T1->SetBranchAddress("Electron_seedGain", Electron_seedGain, &b_Electron_seedGain);
   T1->SetBranchAddress("Flag_ecalBadCalibFilterV2", &Flag_ecalBadCalibFilterV2, &b_Flag_ecalBadCalibFilterV2);
   T1->SetBranchAddress("nFatJet", &nFatJet, &b_nFatJet);
   T1->SetBranchAddress("FatJet_area", FatJet_area, &b_FatJet_area);
   T1->SetBranchAddress("FatJet_btagCMVA", FatJet_btagCMVA, &b_FatJet_btagCMVA);
   T1->SetBranchAddress("FatJet_btagCSVV2", FatJet_btagCSVV2, &b_FatJet_btagCSVV2);
   T1->SetBranchAddress("FatJet_btagDDBvL", FatJet_btagDDBvL, &b_FatJet_btagDDBvL);
   T1->SetBranchAddress("FatJet_btagDDBvL_noMD", FatJet_btagDDBvL_noMD, &b_FatJet_btagDDBvL_noMD);
   T1->SetBranchAddress("FatJet_btagDDCvB", FatJet_btagDDCvB, &b_FatJet_btagDDCvB);
   T1->SetBranchAddress("FatJet_btagDDCvB_noMD", FatJet_btagDDCvB_noMD, &b_FatJet_btagDDCvB_noMD);
   T1->SetBranchAddress("FatJet_btagDDCvL", FatJet_btagDDCvL, &b_FatJet_btagDDCvL);
   T1->SetBranchAddress("FatJet_btagDDCvL_noMD", FatJet_btagDDCvL_noMD, &b_FatJet_btagDDCvL_noMD);
   T1->SetBranchAddress("FatJet_btagDeepB", FatJet_btagDeepB, &b_FatJet_btagDeepB);
   T1->SetBranchAddress("FatJet_btagHbb", FatJet_btagHbb, &b_FatJet_btagHbb);
   T1->SetBranchAddress("FatJet_deepTagMD_H4qvsQCD", FatJet_deepTagMD_H4qvsQCD, &b_FatJet_deepTagMD_H4qvsQCD);
   T1->SetBranchAddress("FatJet_deepTagMD_HbbvsQCD", FatJet_deepTagMD_HbbvsQCD, &b_FatJet_deepTagMD_HbbvsQCD);
   T1->SetBranchAddress("FatJet_deepTagMD_TvsQCD", FatJet_deepTagMD_TvsQCD, &b_FatJet_deepTagMD_TvsQCD);
   T1->SetBranchAddress("FatJet_deepTagMD_WvsQCD", FatJet_deepTagMD_WvsQCD, &b_FatJet_deepTagMD_WvsQCD);
   T1->SetBranchAddress("FatJet_deepTagMD_ZHbbvsQCD", FatJet_deepTagMD_ZHbbvsQCD, &b_FatJet_deepTagMD_ZHbbvsQCD);
   T1->SetBranchAddress("FatJet_deepTagMD_ZHccvsQCD", FatJet_deepTagMD_ZHccvsQCD, &b_FatJet_deepTagMD_ZHccvsQCD);
   T1->SetBranchAddress("FatJet_deepTagMD_ZbbvsQCD", FatJet_deepTagMD_ZbbvsQCD, &b_FatJet_deepTagMD_ZbbvsQCD);
   T1->SetBranchAddress("FatJet_deepTagMD_ZvsQCD", FatJet_deepTagMD_ZvsQCD, &b_FatJet_deepTagMD_ZvsQCD);
   T1->SetBranchAddress("FatJet_deepTagMD_bbvsLight", FatJet_deepTagMD_bbvsLight, &b_FatJet_deepTagMD_bbvsLight);
   T1->SetBranchAddress("FatJet_deepTagMD_ccvsLight", FatJet_deepTagMD_ccvsLight, &b_FatJet_deepTagMD_ccvsLight);
   T1->SetBranchAddress("FatJet_deepTag_H", FatJet_deepTag_H, &b_FatJet_deepTag_H);
   T1->SetBranchAddress("FatJet_deepTag_QCD", FatJet_deepTag_QCD, &b_FatJet_deepTag_QCD);
   T1->SetBranchAddress("FatJet_deepTag_QCDothers", FatJet_deepTag_QCDothers, &b_FatJet_deepTag_QCDothers);
   T1->SetBranchAddress("FatJet_deepTag_TvsQCD", FatJet_deepTag_TvsQCD, &b_FatJet_deepTag_TvsQCD);
   T1->SetBranchAddress("FatJet_deepTag_WvsQCD", FatJet_deepTag_WvsQCD, &b_FatJet_deepTag_WvsQCD);
   T1->SetBranchAddress("FatJet_deepTag_ZvsQCD", FatJet_deepTag_ZvsQCD, &b_FatJet_deepTag_ZvsQCD);
   T1->SetBranchAddress("FatJet_eta", FatJet_eta, &b_FatJet_eta);
   T1->SetBranchAddress("FatJet_mass", FatJet_mass, &b_FatJet_mass);
   T1->SetBranchAddress("FatJet_msoftdrop", FatJet_msoftdrop, &b_FatJet_msoftdrop);
   T1->SetBranchAddress("FatJet_n2b1", FatJet_n2b1, &b_FatJet_n2b1);
   T1->SetBranchAddress("FatJet_n3b1", FatJet_n3b1, &b_FatJet_n3b1);
   T1->SetBranchAddress("FatJet_phi", FatJet_phi, &b_FatJet_phi);
   T1->SetBranchAddress("FatJet_pt", FatJet_pt, &b_FatJet_pt);
   T1->SetBranchAddress("FatJet_rawFactor", FatJet_rawFactor, &b_FatJet_rawFactor);
   T1->SetBranchAddress("FatJet_tau1", FatJet_tau1, &b_FatJet_tau1);
   T1->SetBranchAddress("FatJet_tau2", FatJet_tau2, &b_FatJet_tau2);
   T1->SetBranchAddress("FatJet_tau3", FatJet_tau3, &b_FatJet_tau3);
   T1->SetBranchAddress("FatJet_tau4", FatJet_tau4, &b_FatJet_tau4);
   T1->SetBranchAddress("FatJet_jetId", FatJet_jetId, &b_FatJet_jetId);
   T1->SetBranchAddress("FatJet_subJetIdx1", FatJet_subJetIdx1, &b_FatJet_subJetIdx1);
   T1->SetBranchAddress("FatJet_subJetIdx2", FatJet_subJetIdx2, &b_FatJet_subJetIdx2);
   T1->SetBranchAddress("nFsrPhoton", &nFsrPhoton, &b_nFsrPhoton);
   T1->SetBranchAddress("FsrPhoton_dROverEt2", FsrPhoton_dROverEt2, &b_FsrPhoton_dROverEt2);
   T1->SetBranchAddress("FsrPhoton_eta", FsrPhoton_eta, &b_FsrPhoton_eta);
   T1->SetBranchAddress("FsrPhoton_phi", FsrPhoton_phi, &b_FsrPhoton_phi);
   T1->SetBranchAddress("FsrPhoton_pt", FsrPhoton_pt, &b_FsrPhoton_pt);
   T1->SetBranchAddress("FsrPhoton_relIso03", FsrPhoton_relIso03, &b_FsrPhoton_relIso03);
   T1->SetBranchAddress("FsrPhoton_muonIdx", FsrPhoton_muonIdx, &b_FsrPhoton_muonIdx);
   T1->SetBranchAddress("nGenJetAK8", &nGenJetAK8, &b_nGenJetAK8);
   T1->SetBranchAddress("GenJetAK8_eta", GenJetAK8_eta, &b_GenJetAK8_eta);
   T1->SetBranchAddress("GenJetAK8_mass", GenJetAK8_mass, &b_GenJetAK8_mass);
   T1->SetBranchAddress("GenJetAK8_phi", GenJetAK8_phi, &b_GenJetAK8_phi);
   T1->SetBranchAddress("GenJetAK8_pt", GenJetAK8_pt, &b_GenJetAK8_pt);
   T1->SetBranchAddress("nGenJet", &nGenJet, &b_nGenJet);
   T1->SetBranchAddress("GenJet_eta", GenJet_eta, &b_GenJet_eta);
   T1->SetBranchAddress("GenJet_mass", GenJet_mass, &b_GenJet_mass);
   T1->SetBranchAddress("GenJet_phi", GenJet_phi, &b_GenJet_phi);
   T1->SetBranchAddress("GenJet_pt", GenJet_pt, &b_GenJet_pt);
   T1->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
   T1->SetBranchAddress("GenPart_eta", GenPart_eta, &b_GenPart_eta);
   T1->SetBranchAddress("GenPart_mass", GenPart_mass, &b_GenPart_mass);
   T1->SetBranchAddress("GenPart_phi", GenPart_phi, &b_GenPart_phi);
   T1->SetBranchAddress("GenPart_pt", GenPart_pt, &b_GenPart_pt);
   T1->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother, &b_GenPart_genPartIdxMother);
   T1->SetBranchAddress("GenPart_pdgId", GenPart_pdgId, &b_GenPart_pdgId);
   T1->SetBranchAddress("GenPart_status", GenPart_status, &b_GenPart_status);
   T1->SetBranchAddress("GenPart_statusFlags", GenPart_statusFlags, &b_GenPart_statusFlags);
   T1->SetBranchAddress("nSubGenJetAK8", &nSubGenJetAK8, &b_nSubGenJetAK8);
   T1->SetBranchAddress("SubGenJetAK8_eta", SubGenJetAK8_eta, &b_SubGenJetAK8_eta);
   T1->SetBranchAddress("SubGenJetAK8_mass", SubGenJetAK8_mass, &b_SubGenJetAK8_mass);
   T1->SetBranchAddress("SubGenJetAK8_phi", SubGenJetAK8_phi, &b_SubGenJetAK8_phi);
   T1->SetBranchAddress("SubGenJetAK8_pt", SubGenJetAK8_pt, &b_SubGenJetAK8_pt);
   T1->SetBranchAddress("Generator_binvar", &Generator_binvar, &b_Generator_binvar);
   T1->SetBranchAddress("Generator_scalePDF", &Generator_scalePDF, &b_Generator_scalePDF);
   T1->SetBranchAddress("Generator_weight", &Generator_weight, &b_Generator_weight);
   T1->SetBranchAddress("Generator_x1", &Generator_x1, &b_Generator_x1);
   T1->SetBranchAddress("Generator_x2", &Generator_x2, &b_Generator_x2);
   T1->SetBranchAddress("Generator_xpdf1", &Generator_xpdf1, &b_Generator_xpdf1);
   T1->SetBranchAddress("Generator_xpdf2", &Generator_xpdf2, &b_Generator_xpdf2);
   T1->SetBranchAddress("Generator_id1", &Generator_id1, &b_Generator_id1);
   T1->SetBranchAddress("Generator_id2", &Generator_id2, &b_Generator_id2);
   T1->SetBranchAddress("nGenVisTau", &nGenVisTau, &b_nGenVisTau);
   T1->SetBranchAddress("GenVisTau_eta", GenVisTau_eta, &b_GenVisTau_eta);
   T1->SetBranchAddress("GenVisTau_mass", GenVisTau_mass, &b_GenVisTau_mass);
   T1->SetBranchAddress("GenVisTau_phi", GenVisTau_phi, &b_GenVisTau_phi);
   T1->SetBranchAddress("GenVisTau_pt", GenVisTau_pt, &b_GenVisTau_pt);
   T1->SetBranchAddress("GenVisTau_charge", GenVisTau_charge, &b_GenVisTau_charge);
   T1->SetBranchAddress("GenVisTau_genPartIdxMother", GenVisTau_genPartIdxMother, &b_GenVisTau_genPartIdxMother);
   T1->SetBranchAddress("GenVisTau_status", GenVisTau_status, &b_GenVisTau_status);
   T1->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   T1->SetBranchAddress("LHEWeight_originalXWGTUP", &LHEWeight_originalXWGTUP, &b_LHEWeight_originalXWGTUP);
   T1->SetBranchAddress("nLHEPdfWeight", &nLHEPdfWeight, &b_nLHEPdfWeight);
   T1->SetBranchAddress("LHEPdfWeight", &LHEPdfWeight, &b_LHEPdfWeight);
   T1->SetBranchAddress("nLHEReweightingWeight", &nLHEReweightingWeight, &b_nLHEReweightingWeight);
   T1->SetBranchAddress("LHEReweightingWeight", &LHEReweightingWeight, &b_LHEReweightingWeight);
   T1->SetBranchAddress("nLHEScaleWeight", &nLHEScaleWeight, &b_nLHEScaleWeight);
   T1->SetBranchAddress("LHEScaleWeight", &LHEScaleWeight, &b_LHEScaleWeight);
   T1->SetBranchAddress("nPSWeight", &nPSWeight, &b_nPSWeight);
   T1->SetBranchAddress("PSWeight", PSWeight, &b_PSWeight);
   T1->SetBranchAddress("nIsoTrack", &nIsoTrack, &b_nIsoTrack);
   T1->SetBranchAddress("IsoTrack_dxy", IsoTrack_dxy, &b_IsoTrack_dxy);
   T1->SetBranchAddress("IsoTrack_dz", IsoTrack_dz, &b_IsoTrack_dz);
   T1->SetBranchAddress("IsoTrack_eta", IsoTrack_eta, &b_IsoTrack_eta);
   T1->SetBranchAddress("IsoTrack_pfRelIso03_all", IsoTrack_pfRelIso03_all, &b_IsoTrack_pfRelIso03_all);
   T1->SetBranchAddress("IsoTrack_pfRelIso03_chg", IsoTrack_pfRelIso03_chg, &b_IsoTrack_pfRelIso03_chg);
   T1->SetBranchAddress("IsoTrack_phi", IsoTrack_phi, &b_IsoTrack_phi);
   T1->SetBranchAddress("IsoTrack_pt", IsoTrack_pt, &b_IsoTrack_pt);
   T1->SetBranchAddress("IsoTrack_miniPFRelIso_all", IsoTrack_miniPFRelIso_all, &b_IsoTrack_miniPFRelIso_all);
   T1->SetBranchAddress("IsoTrack_miniPFRelIso_chg", IsoTrack_miniPFRelIso_chg, &b_IsoTrack_miniPFRelIso_chg);
   T1->SetBranchAddress("IsoTrack_fromPV", IsoTrack_fromPV, &b_IsoTrack_fromPV);
   T1->SetBranchAddress("IsoTrack_pdgId", IsoTrack_pdgId, &b_IsoTrack_pdgId);
   T1->SetBranchAddress("IsoTrack_isHighPurityTrack", IsoTrack_isHighPurityTrack, &b_IsoTrack_isHighPurityTrack);
   T1->SetBranchAddress("IsoTrack_isPFcand", IsoTrack_isPFcand, &b_IsoTrack_isPFcand);
   T1->SetBranchAddress("IsoTrack_isFromLostTrack", IsoTrack_isFromLostTrack, &b_IsoTrack_isFromLostTrack);
   T1->SetBranchAddress("nJet", &nJet, &b_nJet);
   T1->SetBranchAddress("Jet_area", Jet_area, &b_Jet_area);
   T1->SetBranchAddress("Jet_btagCMVA", Jet_btagCMVA, &b_Jet_btagCMVA);
   T1->SetBranchAddress("Jet_btagCSVV2", Jet_btagCSVV2, &b_Jet_btagCSVV2);
   T1->SetBranchAddress("Jet_btagDeepB", Jet_btagDeepB, &b_Jet_btagDeepB);
   T1->SetBranchAddress("Jet_btagDeepC", Jet_btagDeepC, &b_Jet_btagDeepC);
   T1->SetBranchAddress("Jet_btagDeepFlavB", Jet_btagDeepFlavB, &b_Jet_btagDeepFlavB);
   T1->SetBranchAddress("Jet_btagDeepFlavC", Jet_btagDeepFlavC, &b_Jet_btagDeepFlavC);
   T1->SetBranchAddress("Jet_chEmEF", Jet_chEmEF, &b_Jet_chEmEF);
   T1->SetBranchAddress("Jet_chHEF", Jet_chHEF, &b_Jet_chHEF);
   T1->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   T1->SetBranchAddress("Jet_jercCHF", Jet_jercCHF, &b_Jet_jercCHF);
   T1->SetBranchAddress("Jet_jercCHPUF", Jet_jercCHPUF, &b_Jet_jercCHPUF);
   T1->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   T1->SetBranchAddress("Jet_muEF", Jet_muEF, &b_Jet_muEF);
   T1->SetBranchAddress("Jet_muonSubtrFactor", Jet_muonSubtrFactor, &b_Jet_muonSubtrFactor);
   T1->SetBranchAddress("Jet_neEmEF", Jet_neEmEF, &b_Jet_neEmEF);
   T1->SetBranchAddress("Jet_neHEF", Jet_neHEF, &b_Jet_neHEF);
   T1->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   T1->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   T1->SetBranchAddress("Jet_qgl", Jet_qgl, &b_Jet_qgl);
   T1->SetBranchAddress("Jet_rawFactor", Jet_rawFactor, &b_Jet_rawFactor);
   T1->SetBranchAddress("Jet_bRegCorr", Jet_bRegCorr, &b_Jet_bRegCorr);
   T1->SetBranchAddress("Jet_bRegRes", Jet_bRegRes, &b_Jet_bRegRes);
   T1->SetBranchAddress("Jet_electronIdx1", Jet_electronIdx1, &b_Jet_electronIdx1);
   T1->SetBranchAddress("Jet_electronIdx2", Jet_electronIdx2, &b_Jet_electronIdx2);
   T1->SetBranchAddress("Jet_jetId", Jet_jetId, &b_Jet_jetId);
   T1->SetBranchAddress("Jet_muonIdx1", Jet_muonIdx1, &b_Jet_muonIdx1);
   T1->SetBranchAddress("Jet_muonIdx2", Jet_muonIdx2, &b_Jet_muonIdx2);
   T1->SetBranchAddress("Jet_nConstituents", Jet_nConstituents, &b_Jet_nConstituents);
   T1->SetBranchAddress("Jet_nElectrons", Jet_nElectrons, &b_Jet_nElectrons);
   T1->SetBranchAddress("Jet_nMuons", Jet_nMuons, &b_Jet_nMuons);
   T1->SetBranchAddress("Jet_puId", Jet_puId, &b_Jet_puId);
   T1->SetBranchAddress("L1PreFiringWeight_Dn", &L1PreFiringWeight_Dn, &b_L1PreFiringWeight_Dn);
   T1->SetBranchAddress("L1PreFiringWeight_Nom", &L1PreFiringWeight_Nom, &b_L1PreFiringWeight_Nom);
   T1->SetBranchAddress("L1PreFiringWeight_Up", &L1PreFiringWeight_Up, &b_L1PreFiringWeight_Up);
   T1->SetBranchAddress("LHE_HT", &LHE_HT, &b_LHE_HT);
   T1->SetBranchAddress("LHE_HTIncoming", &LHE_HTIncoming, &b_LHE_HTIncoming);
   T1->SetBranchAddress("LHE_Vpt", &LHE_Vpt, &b_LHE_Vpt);
   T1->SetBranchAddress("LHE_Njets", &LHE_Njets, &b_LHE_Njets);
   T1->SetBranchAddress("LHE_Nb", &LHE_Nb, &b_LHE_Nb);
   T1->SetBranchAddress("LHE_Nc", &LHE_Nc, &b_LHE_Nc);
   T1->SetBranchAddress("LHE_Nuds", &LHE_Nuds, &b_LHE_Nuds);
   T1->SetBranchAddress("LHE_Nglu", &LHE_Nglu, &b_LHE_Nglu);
   T1->SetBranchAddress("LHE_NpNLO", &LHE_NpNLO, &b_LHE_NpNLO);
   T1->SetBranchAddress("LHE_NpLO", &LHE_NpLO, &b_LHE_NpLO);
   T1->SetBranchAddress("nLHEPart", &nLHEPart, &b_nLHEPart);
   T1->SetBranchAddress("LHEPart_pt", LHEPart_pt, &b_LHEPart_pt);
   T1->SetBranchAddress("LHEPart_eta", LHEPart_eta, &b_LHEPart_eta);
   T1->SetBranchAddress("LHEPart_phi", LHEPart_phi, &b_LHEPart_phi);
   T1->SetBranchAddress("LHEPart_mass", LHEPart_mass, &b_LHEPart_mass);
   T1->SetBranchAddress("LHEPart_pdgId", LHEPart_pdgId, &b_LHEPart_pdgId);
   T1->SetBranchAddress("METFixEE2017_MetUnclustEnUpDeltaX", &METFixEE2017_MetUnclustEnUpDeltaX, &b_METFixEE2017_MetUnclustEnUpDeltaX);
   T1->SetBranchAddress("METFixEE2017_MetUnclustEnUpDeltaY", &METFixEE2017_MetUnclustEnUpDeltaY, &b_METFixEE2017_MetUnclustEnUpDeltaY);
   T1->SetBranchAddress("METFixEE2017_covXX", &METFixEE2017_covXX, &b_METFixEE2017_covXX);
   T1->SetBranchAddress("METFixEE2017_covXY", &METFixEE2017_covXY, &b_METFixEE2017_covXY);
   T1->SetBranchAddress("METFixEE2017_covYY", &METFixEE2017_covYY, &b_METFixEE2017_covYY);
   T1->SetBranchAddress("METFixEE2017_phi", &METFixEE2017_phi, &b_METFixEE2017_phi);
   T1->SetBranchAddress("METFixEE2017_pt", &METFixEE2017_pt, &b_METFixEE2017_pt);
   T1->SetBranchAddress("METFixEE2017_significance", &METFixEE2017_significance, &b_METFixEE2017_significance);
   T1->SetBranchAddress("METFixEE2017_sumEt", &METFixEE2017_sumEt, &b_METFixEE2017_sumEt);
   T1->SetBranchAddress("GenMET_phi", &GenMET_phi, &b_GenMET_phi);
   T1->SetBranchAddress("GenMET_pt", &GenMET_pt, &b_GenMET_pt);
   T1->SetBranchAddress("MET_MetUnclustEnUpDeltaX", &MET_MetUnclustEnUpDeltaX, &b_MET_MetUnclustEnUpDeltaX);
   T1->SetBranchAddress("MET_MetUnclustEnUpDeltaY", &MET_MetUnclustEnUpDeltaY, &b_MET_MetUnclustEnUpDeltaY);
   T1->SetBranchAddress("MET_covXX", &MET_covXX, &b_MET_covXX);
   T1->SetBranchAddress("MET_covXY", &MET_covXY, &b_MET_covXY);
   T1->SetBranchAddress("MET_covYY", &MET_covYY, &b_MET_covYY);
   T1->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   T1->SetBranchAddress("MET_pt", &MET_pt, &b_MET_pt);
   T1->SetBranchAddress("MET_significance", &MET_significance, &b_MET_significance);
   T1->SetBranchAddress("MET_sumEt", &MET_sumEt, &b_MET_sumEt);
   T1->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   T1->SetBranchAddress("Muon_dxy", Muon_dxy, &b_Muon_dxy);
   T1->SetBranchAddress("Muon_dxyErr", Muon_dxyErr, &b_Muon_dxyErr);
   T1->SetBranchAddress("Muon_dz", Muon_dz, &b_Muon_dz);
   T1->SetBranchAddress("Muon_dzErr", Muon_dzErr, &b_Muon_dzErr);
   T1->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   T1->SetBranchAddress("Muon_ip3d", Muon_ip3d, &b_Muon_ip3d);
   T1->SetBranchAddress("Muon_jetPtRelv2", Muon_jetPtRelv2, &b_Muon_jetPtRelv2);
   T1->SetBranchAddress("Muon_jetRelIso", Muon_jetRelIso, &b_Muon_jetRelIso);
   T1->SetBranchAddress("Muon_mass", Muon_mass, &b_Muon_mass);
   T1->SetBranchAddress("Muon_miniPFRelIso_all", Muon_miniPFRelIso_all, &b_Muon_miniPFRelIso_all);
   T1->SetBranchAddress("Muon_miniPFRelIso_chg", Muon_miniPFRelIso_chg, &b_Muon_miniPFRelIso_chg);
   T1->SetBranchAddress("Muon_pfRelIso03_all", Muon_pfRelIso03_all, &b_Muon_pfRelIso03_all);
   T1->SetBranchAddress("Muon_pfRelIso03_chg", Muon_pfRelIso03_chg, &b_Muon_pfRelIso03_chg);
   T1->SetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all, &b_Muon_pfRelIso04_all);
   T1->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   T1->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   T1->SetBranchAddress("Muon_ptErr", Muon_ptErr, &b_Muon_ptErr);
   T1->SetBranchAddress("Muon_segmentComp", Muon_segmentComp, &b_Muon_segmentComp);
   T1->SetBranchAddress("Muon_sip3d", Muon_sip3d, &b_Muon_sip3d);
   T1->SetBranchAddress("Muon_tkRelIso", Muon_tkRelIso, &b_Muon_tkRelIso);
   T1->SetBranchAddress("Muon_tunepRelPt", Muon_tunepRelPt, &b_Muon_tunepRelPt);
   T1->SetBranchAddress("Muon_mvaLowPt", Muon_mvaLowPt, &b_Muon_mvaLowPt);
   T1->SetBranchAddress("Muon_mvaTTH", Muon_mvaTTH, &b_Muon_mvaTTH);
   T1->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
   T1->SetBranchAddress("Muon_jetIdx", Muon_jetIdx, &b_Muon_jetIdx);
   T1->SetBranchAddress("Muon_nStations", Muon_nStations, &b_Muon_nStations);
   T1->SetBranchAddress("Muon_nTrackerLayers", Muon_nTrackerLayers, &b_Muon_nTrackerLayers);
   T1->SetBranchAddress("Muon_pdgId", Muon_pdgId, &b_Muon_pdgId);
   T1->SetBranchAddress("Muon_tightCharge", Muon_tightCharge, &b_Muon_tightCharge);
   T1->SetBranchAddress("Muon_fsrPhotonIdx", Muon_fsrPhotonIdx, &b_Muon_fsrPhotonIdx);
   T1->SetBranchAddress("Muon_highPtId", Muon_highPtId, &b_Muon_highPtId);
   T1->SetBranchAddress("Muon_inTimeMuon", Muon_inTimeMuon, &b_Muon_inTimeMuon);
   T1->SetBranchAddress("Muon_isGlobal", Muon_isGlobal, &b_Muon_isGlobal);
   T1->SetBranchAddress("Muon_isPFcand", Muon_isPFcand, &b_Muon_isPFcand);
   T1->SetBranchAddress("Muon_isTracker", Muon_isTracker, &b_Muon_isTracker);
   T1->SetBranchAddress("Muon_looseId", Muon_looseId, &b_Muon_looseId);
   T1->SetBranchAddress("Muon_mediumId", Muon_mediumId, &b_Muon_mediumId);
   T1->SetBranchAddress("Muon_mediumPromptId", Muon_mediumPromptId, &b_Muon_mediumPromptId);
   T1->SetBranchAddress("Muon_miniIsoId", Muon_miniIsoId, &b_Muon_miniIsoId);
   T1->SetBranchAddress("Muon_multiIsoId", Muon_multiIsoId, &b_Muon_multiIsoId);
   T1->SetBranchAddress("Muon_mvaId", Muon_mvaId, &b_Muon_mvaId);
   T1->SetBranchAddress("Muon_pfIsoId", Muon_pfIsoId, &b_Muon_pfIsoId);
   T1->SetBranchAddress("Muon_softId", Muon_softId, &b_Muon_softId);
   T1->SetBranchAddress("Muon_softMvaId", Muon_softMvaId, &b_Muon_softMvaId);
   T1->SetBranchAddress("Muon_tightId", Muon_tightId, &b_Muon_tightId);
   T1->SetBranchAddress("Muon_tkIsoId", Muon_tkIsoId, &b_Muon_tkIsoId);
   T1->SetBranchAddress("Muon_triggerIdLoose", Muon_triggerIdLoose, &b_Muon_triggerIdLoose);
   T1->SetBranchAddress("nPhoton", &nPhoton, &b_nPhoton);
   T1->SetBranchAddress("Photon_eCorr", Photon_eCorr, &b_Photon_eCorr);
   T1->SetBranchAddress("Photon_energyErr", Photon_energyErr, &b_Photon_energyErr);
   T1->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta);
   T1->SetBranchAddress("Photon_hoe", Photon_hoe, &b_Photon_hoe);
   T1->SetBranchAddress("Photon_mass", Photon_mass, &b_Photon_mass);
   T1->SetBranchAddress("Photon_mvaID", Photon_mvaID, &b_Photon_mvaID);
   T1->SetBranchAddress("Photon_mvaIDV1", Photon_mvaIDV1, &b_Photon_mvaIDV1);
   T1->SetBranchAddress("Photon_pfRelIso03_all", Photon_pfRelIso03_all, &b_Photon_pfRelIso03_all);
   T1->SetBranchAddress("Photon_pfRelIso03_chg", Photon_pfRelIso03_chg, &b_Photon_pfRelIso03_chg);
   T1->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);
   T1->SetBranchAddress("Photon_pt", Photon_pt, &b_Photon_pt);
   T1->SetBranchAddress("Photon_r9", Photon_r9, &b_Photon_r9);
   T1->SetBranchAddress("Photon_sieie", Photon_sieie, &b_Photon_sieie);
   T1->SetBranchAddress("Photon_charge", Photon_charge, &b_Photon_charge);
   T1->SetBranchAddress("Photon_cutBasedBitmap", Photon_cutBasedBitmap, &b_Photon_cutBasedBitmap);
   T1->SetBranchAddress("Photon_cutBasedV1Bitmap", Photon_cutBasedV1Bitmap, &b_Photon_cutBasedV1Bitmap);
   T1->SetBranchAddress("Photon_electronIdx", Photon_electronIdx, &b_Photon_electronIdx);
   T1->SetBranchAddress("Photon_jetIdx", Photon_jetIdx, &b_Photon_jetIdx);
   T1->SetBranchAddress("Photon_pdgId", Photon_pdgId, &b_Photon_pdgId);
   T1->SetBranchAddress("Photon_vidNestedWPBitmap", Photon_vidNestedWPBitmap, &b_Photon_vidNestedWPBitmap);
   T1->SetBranchAddress("Photon_electronVeto", Photon_electronVeto, &b_Photon_electronVeto);
   T1->SetBranchAddress("Photon_isScEtaEB", Photon_isScEtaEB, &b_Photon_isScEtaEB);
   T1->SetBranchAddress("Photon_isScEtaEE", Photon_isScEtaEE, &b_Photon_isScEtaEE);
   T1->SetBranchAddress("Photon_mvaID_WP80", Photon_mvaID_WP80, &b_Photon_mvaID_WP80);
   T1->SetBranchAddress("Photon_mvaID_WP90", Photon_mvaID_WP90, &b_Photon_mvaID_WP90);
   T1->SetBranchAddress("Photon_pixelSeed", Photon_pixelSeed, &b_Photon_pixelSeed);
   T1->SetBranchAddress("Photon_seedGain", Photon_seedGain, &b_Photon_seedGain);
   T1->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt, &b_Pileup_nTrueInt);
   T1->SetBranchAddress("Pileup_pudensity", &Pileup_pudensity, &b_Pileup_pudensity);
   T1->SetBranchAddress("Pileup_gpudensity", &Pileup_gpudensity, &b_Pileup_gpudensity);
   T1->SetBranchAddress("Pileup_nPU", &Pileup_nPU, &b_Pileup_nPU);
   T1->SetBranchAddress("Pileup_sumEOOT", &Pileup_sumEOOT, &b_Pileup_sumEOOT);
   T1->SetBranchAddress("Pileup_sumLOOT", &Pileup_sumLOOT, &b_Pileup_sumLOOT);
   T1->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi, &b_PuppiMET_phi);
   T1->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt, &b_PuppiMET_pt);
   T1->SetBranchAddress("PuppiMET_sumEt", &PuppiMET_sumEt, &b_PuppiMET_sumEt);
   T1->SetBranchAddress("RawMET_phi", &RawMET_phi, &b_RawMET_phi);
   T1->SetBranchAddress("RawMET_pt", &RawMET_pt, &b_RawMET_pt);
   T1->SetBranchAddress("RawMET_sumEt", &RawMET_sumEt, &b_RawMET_sumEt);
   T1->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
   T1->SetBranchAddress("fixedGridRhoFastjetCentral", &fixedGridRhoFastjetCentral, &b_fixedGridRhoFastjetCentral);
   T1->SetBranchAddress("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, &b_fixedGridRhoFastjetCentralCalo);
   T1->SetBranchAddress("fixedGridRhoFastjetCentralChargedPileUp", &fixedGridRhoFastjetCentralChargedPileUp, &b_fixedGridRhoFastjetCentralChargedPileUp);
   T1->SetBranchAddress("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, &b_fixedGridRhoFastjetCentralNeutral);
   T1->SetBranchAddress("nGenDressedLepton", &nGenDressedLepton, &b_nGenDressedLepton);
   T1->SetBranchAddress("GenDressedLepton_eta", GenDressedLepton_eta, &b_GenDressedLepton_eta);
   T1->SetBranchAddress("GenDressedLepton_mass", GenDressedLepton_mass, &b_GenDressedLepton_mass);
   T1->SetBranchAddress("GenDressedLepton_phi", GenDressedLepton_phi, &b_GenDressedLepton_phi);
   T1->SetBranchAddress("GenDressedLepton_pt", GenDressedLepton_pt, &b_GenDressedLepton_pt);
   T1->SetBranchAddress("GenDressedLepton_pdgId", GenDressedLepton_pdgId, &b_GenDressedLepton_pdgId);
   T1->SetBranchAddress("GenDressedLepton_hasTauAnc", GenDressedLepton_hasTauAnc, &b_GenDressedLepton_hasTauAnc);
   T1->SetBranchAddress("nSoftActivityJet", &nSoftActivityJet, &b_nSoftActivityJet);
   T1->SetBranchAddress("SoftActivityJet_eta", SoftActivityJet_eta, &b_SoftActivityJet_eta);
   T1->SetBranchAddress("SoftActivityJet_phi", SoftActivityJet_phi, &b_SoftActivityJet_phi);
   T1->SetBranchAddress("SoftActivityJet_pt", SoftActivityJet_pt, &b_SoftActivityJet_pt);
   T1->SetBranchAddress("SoftActivityJetHT", &SoftActivityJetHT, &b_SoftActivityJetHT);
   T1->SetBranchAddress("SoftActivityJetHT10", &SoftActivityJetHT10, &b_SoftActivityJetHT10);
   T1->SetBranchAddress("SoftActivityJetHT2", &SoftActivityJetHT2, &b_SoftActivityJetHT2);
   T1->SetBranchAddress("SoftActivityJetHT5", &SoftActivityJetHT5, &b_SoftActivityJetHT5);
   T1->SetBranchAddress("SoftActivityJetNjets10", &SoftActivityJetNjets10, &b_SoftActivityJetNjets10);
   T1->SetBranchAddress("SoftActivityJetNjets2", &SoftActivityJetNjets2, &b_SoftActivityJetNjets2);
   T1->SetBranchAddress("SoftActivityJetNjets5", &SoftActivityJetNjets5, &b_SoftActivityJetNjets5);
   T1->SetBranchAddress("nSubJet", &nSubJet, &b_nSubJet);
   T1->SetBranchAddress("SubJet_btagCMVA", SubJet_btagCMVA, &b_SubJet_btagCMVA);
   T1->SetBranchAddress("SubJet_btagCSVV2", SubJet_btagCSVV2, &b_SubJet_btagCSVV2);
   T1->SetBranchAddress("SubJet_btagDeepB", SubJet_btagDeepB, &b_SubJet_btagDeepB);
   T1->SetBranchAddress("SubJet_eta", SubJet_eta, &b_SubJet_eta);
   T1->SetBranchAddress("SubJet_mass", SubJet_mass, &b_SubJet_mass);
   T1->SetBranchAddress("SubJet_n2b1", SubJet_n2b1, &b_SubJet_n2b1);
   T1->SetBranchAddress("SubJet_n3b1", SubJet_n3b1, &b_SubJet_n3b1);
   T1->SetBranchAddress("SubJet_phi", SubJet_phi, &b_SubJet_phi);
   T1->SetBranchAddress("SubJet_pt", SubJet_pt, &b_SubJet_pt);
   T1->SetBranchAddress("SubJet_rawFactor", SubJet_rawFactor, &b_SubJet_rawFactor);
   T1->SetBranchAddress("SubJet_tau1", SubJet_tau1, &b_SubJet_tau1);
   T1->SetBranchAddress("SubJet_tau2", SubJet_tau2, &b_SubJet_tau2);
   T1->SetBranchAddress("SubJet_tau3", SubJet_tau3, &b_SubJet_tau3);
   T1->SetBranchAddress("SubJet_tau4", SubJet_tau4, &b_SubJet_tau4);
   T1->SetBranchAddress("nTau", &nTau, &b_nTau);
   T1->SetBranchAddress("Tau_chargedIso", Tau_chargedIso, &b_Tau_chargedIso);
   T1->SetBranchAddress("Tau_dxy", Tau_dxy, &b_Tau_dxy);
   T1->SetBranchAddress("Tau_dz", Tau_dz, &b_Tau_dz);
   T1->SetBranchAddress("Tau_eta", Tau_eta, &b_Tau_eta);
   T1->SetBranchAddress("Tau_leadTkDeltaEta", Tau_leadTkDeltaEta, &b_Tau_leadTkDeltaEta);
   T1->SetBranchAddress("Tau_leadTkDeltaPhi", Tau_leadTkDeltaPhi, &b_Tau_leadTkDeltaPhi);
   T1->SetBranchAddress("Tau_leadTkPtOverTauPt", Tau_leadTkPtOverTauPt, &b_Tau_leadTkPtOverTauPt);
   T1->SetBranchAddress("Tau_mass", Tau_mass, &b_Tau_mass);
   T1->SetBranchAddress("Tau_neutralIso", Tau_neutralIso, &b_Tau_neutralIso);
   T1->SetBranchAddress("Tau_phi", Tau_phi, &b_Tau_phi);
   T1->SetBranchAddress("Tau_photonsOutsideSignalCone", Tau_photonsOutsideSignalCone, &b_Tau_photonsOutsideSignalCone);
   T1->SetBranchAddress("Tau_pt", Tau_pt, &b_Tau_pt);
   T1->SetBranchAddress("Tau_puCorr", Tau_puCorr, &b_Tau_puCorr);
   T1->SetBranchAddress("Tau_rawAntiEle", Tau_rawAntiEle, &b_Tau_rawAntiEle);
   T1->SetBranchAddress("Tau_rawAntiEle2018", Tau_rawAntiEle2018, &b_Tau_rawAntiEle2018);
   T1->SetBranchAddress("Tau_rawDeepTau2017v2p1VSe", Tau_rawDeepTau2017v2p1VSe, &b_Tau_rawDeepTau2017v2p1VSe);
   T1->SetBranchAddress("Tau_rawDeepTau2017v2p1VSjet", Tau_rawDeepTau2017v2p1VSjet, &b_Tau_rawDeepTau2017v2p1VSjet);
   T1->SetBranchAddress("Tau_rawDeepTau2017v2p1VSmu", Tau_rawDeepTau2017v2p1VSmu, &b_Tau_rawDeepTau2017v2p1VSmu);
   T1->SetBranchAddress("Tau_rawIso", Tau_rawIso, &b_Tau_rawIso);
   T1->SetBranchAddress("Tau_rawIsodR03", Tau_rawIsodR03, &b_Tau_rawIsodR03);
   T1->SetBranchAddress("Tau_rawMVAnewDM2017v2", Tau_rawMVAnewDM2017v2, &b_Tau_rawMVAnewDM2017v2);
   T1->SetBranchAddress("Tau_rawMVAoldDM", Tau_rawMVAoldDM, &b_Tau_rawMVAoldDM);
   T1->SetBranchAddress("Tau_rawMVAoldDM2017v1", Tau_rawMVAoldDM2017v1, &b_Tau_rawMVAoldDM2017v1);
   T1->SetBranchAddress("Tau_rawMVAoldDM2017v2", Tau_rawMVAoldDM2017v2, &b_Tau_rawMVAoldDM2017v2);
   T1->SetBranchAddress("Tau_rawMVAoldDMdR032017v2", Tau_rawMVAoldDMdR032017v2, &b_Tau_rawMVAoldDMdR032017v2);
   T1->SetBranchAddress("Tau_charge", Tau_charge, &b_Tau_charge);
   T1->SetBranchAddress("Tau_decayMode", Tau_decayMode, &b_Tau_decayMode);
   T1->SetBranchAddress("Tau_jetIdx", Tau_jetIdx, &b_Tau_jetIdx);
   T1->SetBranchAddress("Tau_rawAntiEleCat", Tau_rawAntiEleCat, &b_Tau_rawAntiEleCat);
   T1->SetBranchAddress("Tau_rawAntiEleCat2018", Tau_rawAntiEleCat2018, &b_Tau_rawAntiEleCat2018);
   T1->SetBranchAddress("Tau_idAntiEle", Tau_idAntiEle, &b_Tau_idAntiEle);
   T1->SetBranchAddress("Tau_idAntiEle2018", Tau_idAntiEle2018, &b_Tau_idAntiEle2018);
   T1->SetBranchAddress("Tau_idAntiMu", Tau_idAntiMu, &b_Tau_idAntiMu);
   T1->SetBranchAddress("Tau_idDecayMode", Tau_idDecayMode, &b_Tau_idDecayMode);
   T1->SetBranchAddress("Tau_idDecayModeNewDMs", Tau_idDecayModeNewDMs, &b_Tau_idDecayModeNewDMs);
   T1->SetBranchAddress("Tau_idDeepTau2017v2p1VSe", Tau_idDeepTau2017v2p1VSe, &b_Tau_idDeepTau2017v2p1VSe);
   T1->SetBranchAddress("Tau_idDeepTau2017v2p1VSjet", Tau_idDeepTau2017v2p1VSjet, &b_Tau_idDeepTau2017v2p1VSjet);
   T1->SetBranchAddress("Tau_idDeepTau2017v2p1VSmu", Tau_idDeepTau2017v2p1VSmu, &b_Tau_idDeepTau2017v2p1VSmu);
   T1->SetBranchAddress("Tau_idMVAnewDM2017v2", Tau_idMVAnewDM2017v2, &b_Tau_idMVAnewDM2017v2);
   T1->SetBranchAddress("Tau_idMVAoldDM", Tau_idMVAoldDM, &b_Tau_idMVAoldDM);
   T1->SetBranchAddress("Tau_idMVAoldDM2017v1", Tau_idMVAoldDM2017v1, &b_Tau_idMVAoldDM2017v1);
   T1->SetBranchAddress("Tau_idMVAoldDM2017v2", Tau_idMVAoldDM2017v2, &b_Tau_idMVAoldDM2017v2);
   T1->SetBranchAddress("Tau_idMVAoldDMdR032017v2", Tau_idMVAoldDMdR032017v2, &b_Tau_idMVAoldDMdR032017v2);
   T1->SetBranchAddress("TkMET_phi", &TkMET_phi, &b_TkMET_phi);
   T1->SetBranchAddress("TkMET_pt", &TkMET_pt, &b_TkMET_pt);
   T1->SetBranchAddress("TkMET_sumEt", &TkMET_sumEt, &b_TkMET_sumEt);
  /*
   T1->SetBranchAddress("nTrigObj", &nTrigObj, &b_nTrigObj);
   T1->SetBranchAddress("TrigObj_pt", TrigObj_pt, &b_TrigObj_pt);
   T1->SetBranchAddress("TrigObj_eta", TrigObj_eta, &b_TrigObj_eta);
   T1->SetBranchAddress("TrigObj_phi", TrigObj_phi, &b_TrigObj_phi);
   T1->SetBranchAddress("TrigObj_l1pt", TrigObj_l1pt, &b_TrigObj_l1pt);
   T1->SetBranchAddress("TrigObj_l1pt_2", TrigObj_l1pt_2, &b_TrigObj_l1pt_2);
   T1->SetBranchAddress("TrigObj_l2pt", TrigObj_l2pt, &b_TrigObj_l2pt);
   T1->SetBranchAddress("TrigObj_id", TrigObj_id, &b_TrigObj_id);
   T1->SetBranchAddress("TrigObj_l1iso", TrigObj_l1iso, &b_TrigObj_l1iso);
   T1->SetBranchAddress("TrigObj_l1charge", TrigObj_l1charge, &b_TrigObj_l1charge);
   T1->SetBranchAddress("TrigObj_filterBits", TrigObj_filterBits, &b_TrigObj_filterBits);
   */ 
   T1->SetBranchAddress("genTtbarId", &genTtbarId, &b_genTtbarId);
   T1->SetBranchAddress("nOtherPV", &nOtherPV, &b_nOtherPV);
   T1->SetBranchAddress("OtherPV_z", OtherPV_z, &b_OtherPV_z);
   T1->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   T1->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   T1->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   T1->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   T1->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   T1->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
   T1->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   T1->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   T1->SetBranchAddress("nSV", &nSV, &b_nSV);
   T1->SetBranchAddress("SV_dlen", SV_dlen, &b_SV_dlen);
   T1->SetBranchAddress("SV_dlenSig", SV_dlenSig, &b_SV_dlenSig);
   T1->SetBranchAddress("SV_dxy", SV_dxy, &b_SV_dxy);
   T1->SetBranchAddress("SV_dxySig", SV_dxySig, &b_SV_dxySig);
   T1->SetBranchAddress("SV_pAngle", SV_pAngle, &b_SV_pAngle);
   T1->SetBranchAddress("Electron_genPartIdx", Electron_genPartIdx, &b_Electron_genPartIdx);
   T1->SetBranchAddress("Electron_genPartFlav", Electron_genPartFlav, &b_Electron_genPartFlav);
   T1->SetBranchAddress("GenJetAK8_partonFlavour", GenJetAK8_partonFlavour, &b_GenJetAK8_partonFlavour);
   T1->SetBranchAddress("GenJetAK8_hadronFlavour", GenJetAK8_hadronFlavour, &b_GenJetAK8_hadronFlavour);
   T1->SetBranchAddress("GenJet_partonFlavour", GenJet_partonFlavour, &b_GenJet_partonFlavour);
   T1->SetBranchAddress("GenJet_hadronFlavour", GenJet_hadronFlavour, &b_GenJet_hadronFlavour);
   T1->SetBranchAddress("Jet_genJetIdx", Jet_genJetIdx, &b_Jet_genJetIdx);
   T1->SetBranchAddress("Jet_hadronFlavour", Jet_hadronFlavour, &b_Jet_hadronFlavour);
   T1->SetBranchAddress("Jet_partonFlavour", Jet_partonFlavour, &b_Jet_partonFlavour);
   T1->SetBranchAddress("Muon_genPartIdx", Muon_genPartIdx, &b_Muon_genPartIdx);
   T1->SetBranchAddress("Muon_genPartFlav", Muon_genPartFlav, &b_Muon_genPartFlav);
   T1->SetBranchAddress("Photon_genPartIdx", Photon_genPartIdx, &b_Photon_genPartIdx);
   T1->SetBranchAddress("Photon_genPartFlav", Photon_genPartFlav, &b_Photon_genPartFlav);
   T1->SetBranchAddress("MET_fiducialGenPhi", &MET_fiducialGenPhi, &b_MET_fiducialGenPhi);
   T1->SetBranchAddress("MET_fiducialGenPt", &MET_fiducialGenPt, &b_MET_fiducialGenPt);
   T1->SetBranchAddress("Electron_cleanmask", Electron_cleanmask, &b_Electron_cleanmask);
   T1->SetBranchAddress("Jet_cleanmask", Jet_cleanmask, &b_Jet_cleanmask);
   T1->SetBranchAddress("Muon_cleanmask", Muon_cleanmask, &b_Muon_cleanmask);
   T1->SetBranchAddress("Photon_cleanmask", Photon_cleanmask, &b_Photon_cleanmask);
   T1->SetBranchAddress("Tau_cleanmask", Tau_cleanmask, &b_Tau_cleanmask);
   T1->SetBranchAddress("SV_chi2", SV_chi2, &b_SV_chi2);
   T1->SetBranchAddress("SV_eta", SV_eta, &b_SV_eta);
   T1->SetBranchAddress("SV_mass", SV_mass, &b_SV_mass);
   T1->SetBranchAddress("SV_ndof", SV_ndof, &b_SV_ndof);
   T1->SetBranchAddress("SV_phi", SV_phi, &b_SV_phi);
   T1->SetBranchAddress("SV_pt", SV_pt, &b_SV_pt);
   T1->SetBranchAddress("SV_x", SV_x, &b_SV_x);
   T1->SetBranchAddress("SV_y", SV_y, &b_SV_y);
   T1->SetBranchAddress("SV_z", SV_z, &b_SV_z);
   T1->SetBranchAddress("Tau_genPartIdx", Tau_genPartIdx, &b_Tau_genPartIdx);
   T1->SetBranchAddress("Tau_genPartFlav", Tau_genPartFlav, &b_Tau_genPartFlav);
   T1->SetBranchAddress("L1simulation_step", &L1simulation_step, &b_L1simulation_step);
   T1->SetBranchAddress("HLTriggerFirstPath", &HLTriggerFirstPath, &b_HLTriggerFirstPath);
   T1->SetBranchAddress("HLT_AK8PFJet360_TrimMass30", &HLT_AK8PFJet360_TrimMass30, &b_HLT_AK8PFJet360_TrimMass30);
   T1->SetBranchAddress("HLT_AK8PFJet380_TrimMass30", &HLT_AK8PFJet380_TrimMass30, &b_HLT_AK8PFJet380_TrimMass30);
   T1->SetBranchAddress("HLT_AK8PFJet400_TrimMass30", &HLT_AK8PFJet400_TrimMass30, &b_HLT_AK8PFJet400_TrimMass30);
   T1->SetBranchAddress("HLT_AK8PFJet420_TrimMass30", &HLT_AK8PFJet420_TrimMass30, &b_HLT_AK8PFJet420_TrimMass30);
   T1->SetBranchAddress("HLT_AK8PFHT750_TrimMass50", &HLT_AK8PFHT750_TrimMass50, &b_HLT_AK8PFHT750_TrimMass50);
   T1->SetBranchAddress("HLT_AK8PFHT800_TrimMass50", &HLT_AK8PFHT800_TrimMass50, &b_HLT_AK8PFHT800_TrimMass50);
   T1->SetBranchAddress("HLT_AK8PFHT850_TrimMass50", &HLT_AK8PFHT850_TrimMass50, &b_HLT_AK8PFHT850_TrimMass50);
   T1->SetBranchAddress("HLT_AK8PFHT900_TrimMass50", &HLT_AK8PFHT900_TrimMass50, &b_HLT_AK8PFHT900_TrimMass50);
   T1->SetBranchAddress("HLT_CaloJet500_NoJetID", &HLT_CaloJet500_NoJetID, &b_HLT_CaloJet500_NoJetID);
   T1->SetBranchAddress("HLT_CaloJet550_NoJetID", &HLT_CaloJet550_NoJetID, &b_HLT_CaloJet550_NoJetID);
   T1->SetBranchAddress("HLT_Trimuon5_3p5_2_Upsilon_Muon", &HLT_Trimuon5_3p5_2_Upsilon_Muon, &b_HLT_Trimuon5_3p5_2_Upsilon_Muon);
   T1->SetBranchAddress("HLT_DoubleEle25_CaloIdL_MW", &HLT_DoubleEle25_CaloIdL_MW, &b_HLT_DoubleEle25_CaloIdL_MW);
   T1->SetBranchAddress("HLT_DoubleEle27_CaloIdL_MW", &HLT_DoubleEle27_CaloIdL_MW, &b_HLT_DoubleEle27_CaloIdL_MW);
   T1->SetBranchAddress("HLT_DoubleEle33_CaloIdL_MW", &HLT_DoubleEle33_CaloIdL_MW, &b_HLT_DoubleEle33_CaloIdL_MW);
   T1->SetBranchAddress("HLT_DoubleEle24_eta2p1_WPTight_Gsf", &HLT_DoubleEle24_eta2p1_WPTight_Gsf, &b_HLT_DoubleEle24_eta2p1_WPTight_Gsf);
   T1->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350, &b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350);
   T1->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350, &b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350);
   T1->SetBranchAddress("HLT_Ele27_Ele37_CaloIdL_MW", &HLT_Ele27_Ele37_CaloIdL_MW, &b_HLT_Ele27_Ele37_CaloIdL_MW);
   T1->SetBranchAddress("HLT_Mu27_Ele37_CaloIdL_MW", &HLT_Mu27_Ele37_CaloIdL_MW, &b_HLT_Mu27_Ele37_CaloIdL_MW);
   T1->SetBranchAddress("HLT_Mu37_Ele27_CaloIdL_MW", &HLT_Mu37_Ele27_CaloIdL_MW, &b_HLT_Mu37_Ele27_CaloIdL_MW);
   T1->SetBranchAddress("HLT_Mu37_TkMu27", &HLT_Mu37_TkMu27, &b_HLT_Mu37_TkMu27);
   T1->SetBranchAddress("HLT_DoubleMu4_3_Bs", &HLT_DoubleMu4_3_Bs, &b_HLT_DoubleMu4_3_Bs);
   T1->SetBranchAddress("HLT_DoubleMu4_3_Jpsi_Displaced", &HLT_DoubleMu4_3_Jpsi_Displaced, &b_HLT_DoubleMu4_3_Jpsi_Displaced);
   T1->SetBranchAddress("HLT_DoubleMu4_JpsiTrk_Displaced", &HLT_DoubleMu4_JpsiTrk_Displaced, &b_HLT_DoubleMu4_JpsiTrk_Displaced);
   T1->SetBranchAddress("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced", &HLT_DoubleMu4_LowMassNonResonantTrk_Displaced, &b_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced);
   T1->SetBranchAddress("HLT_DoubleMu3_Trk_Tau3mu", &HLT_DoubleMu3_Trk_Tau3mu, &b_HLT_DoubleMu3_Trk_Tau3mu);
   T1->SetBranchAddress("HLT_DoubleMu4_PsiPrimeTrk_Displaced", &HLT_DoubleMu4_PsiPrimeTrk_Displaced, &b_HLT_DoubleMu4_PsiPrimeTrk_Displaced);
   T1->SetBranchAddress("HLT_DoubleMu4_Mass8_DZ_PFHT350", &HLT_DoubleMu4_Mass8_DZ_PFHT350, &b_HLT_DoubleMu4_Mass8_DZ_PFHT350);
   T1->SetBranchAddress("HLT_DoubleMu8_Mass8_PFHT350", &HLT_DoubleMu8_Mass8_PFHT350, &b_HLT_DoubleMu8_Mass8_PFHT350);
   T1->SetBranchAddress("HLT_Mu3_PFJet40", &HLT_Mu3_PFJet40, &b_HLT_Mu3_PFJet40);
   T1->SetBranchAddress("HLT_Mu7p5_L2Mu2_Jpsi", &HLT_Mu7p5_L2Mu2_Jpsi, &b_HLT_Mu7p5_L2Mu2_Jpsi);
   T1->SetBranchAddress("HLT_Mu7p5_L2Mu2_Upsilon", &HLT_Mu7p5_L2Mu2_Upsilon, &b_HLT_Mu7p5_L2Mu2_Upsilon);
   T1->SetBranchAddress("HLT_Mu7p5_Track2_Jpsi", &HLT_Mu7p5_Track2_Jpsi, &b_HLT_Mu7p5_Track2_Jpsi);
   T1->SetBranchAddress("HLT_Mu7p5_Track3p5_Jpsi", &HLT_Mu7p5_Track3p5_Jpsi, &b_HLT_Mu7p5_Track3p5_Jpsi);
   T1->SetBranchAddress("HLT_Mu7p5_Track7_Jpsi", &HLT_Mu7p5_Track7_Jpsi, &b_HLT_Mu7p5_Track7_Jpsi);
   T1->SetBranchAddress("HLT_Mu7p5_Track2_Upsilon", &HLT_Mu7p5_Track2_Upsilon, &b_HLT_Mu7p5_Track2_Upsilon);
   T1->SetBranchAddress("HLT_Mu7p5_Track3p5_Upsilon", &HLT_Mu7p5_Track3p5_Upsilon, &b_HLT_Mu7p5_Track3p5_Upsilon);
   T1->SetBranchAddress("HLT_Mu7p5_Track7_Upsilon", &HLT_Mu7p5_Track7_Upsilon, &b_HLT_Mu7p5_Track7_Upsilon);
   T1->SetBranchAddress("HLT_DoublePhoton33_CaloIdL", &HLT_DoublePhoton33_CaloIdL, &b_HLT_DoublePhoton33_CaloIdL);
   T1->SetBranchAddress("HLT_DoublePhoton70", &HLT_DoublePhoton70, &b_HLT_DoublePhoton70);
   T1->SetBranchAddress("HLT_DoublePhoton85", &HLT_DoublePhoton85, &b_HLT_DoublePhoton85);
   T1->SetBranchAddress("HLT_Ele20_WPTight_Gsf", &HLT_Ele20_WPTight_Gsf, &b_HLT_Ele20_WPTight_Gsf);
   T1->SetBranchAddress("HLT_Ele20_WPLoose_Gsf", &HLT_Ele20_WPLoose_Gsf, &b_HLT_Ele20_WPLoose_Gsf);
   T1->SetBranchAddress("HLT_Ele20_eta2p1_WPLoose_Gsf", &HLT_Ele20_eta2p1_WPLoose_Gsf, &b_HLT_Ele20_eta2p1_WPLoose_Gsf);
   T1->SetBranchAddress("HLT_DiEle27_WPTightCaloOnly_L1DoubleEG", &HLT_DiEle27_WPTightCaloOnly_L1DoubleEG, &b_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG);
   T1->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_Ele27_WPTight_Gsf, &b_HLT_Ele27_WPTight_Gsf);
   T1->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf, &b_HLT_Ele32_WPTight_Gsf);
   T1->SetBranchAddress("HLT_Ele35_WPTight_Gsf", &HLT_Ele35_WPTight_Gsf, &b_HLT_Ele35_WPTight_Gsf);
   T1->SetBranchAddress("HLT_Ele35_WPTight_Gsf_L1EGMT", &HLT_Ele35_WPTight_Gsf_L1EGMT, &b_HLT_Ele35_WPTight_Gsf_L1EGMT);
   T1->SetBranchAddress("HLT_Ele38_WPTight_Gsf", &HLT_Ele38_WPTight_Gsf, &b_HLT_Ele38_WPTight_Gsf);
   T1->SetBranchAddress("HLT_Ele40_WPTight_Gsf", &HLT_Ele40_WPTight_Gsf, &b_HLT_Ele40_WPTight_Gsf);
   T1->SetBranchAddress("HLT_Ele32_WPTight_Gsf_L1DoubleEG", &HLT_Ele32_WPTight_Gsf_L1DoubleEG, &b_HLT_Ele32_WPTight_Gsf_L1DoubleEG);
   T1->SetBranchAddress("HLT_HT450_Beamspot", &HLT_HT450_Beamspot, &b_HLT_HT450_Beamspot);
   T1->SetBranchAddress("HLT_HT300_Beamspot", &HLT_HT300_Beamspot, &b_HLT_HT300_Beamspot);
   T1->SetBranchAddress("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1", &HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1, &b_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1);
   T1->SetBranchAddress("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1", &HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1, &b_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1);
   T1->SetBranchAddress("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1", &HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1, &b_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1);
   T1->SetBranchAddress("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1", &HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1, &b_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1);
   T1->SetBranchAddress("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1", &HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1, &b_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1);
   T1->SetBranchAddress("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1", &HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1, &b_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1);
   T1->SetBranchAddress("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1", &HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1, &b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1);
   T1->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1);
   T1->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1, &b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1);
   T1->SetBranchAddress("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1", &HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1, &b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1);
   T1->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_TightID_SingleL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_TightID_SingleL1, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_TightID_SingleL1);
   T1->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_TightID_SingleL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_TightID_SingleL1, &b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_TightID_SingleL1);
   T1->SetBranchAddress("HLT_IsoMu20", &HLT_IsoMu20, &b_HLT_IsoMu20);
   T1->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24, &b_HLT_IsoMu24);
   T1->SetBranchAddress("HLT_IsoMu24_eta2p1", &HLT_IsoMu24_eta2p1, &b_HLT_IsoMu24_eta2p1);
   T1->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27, &b_HLT_IsoMu27);
   T1->SetBranchAddress("HLT_IsoMu30", &HLT_IsoMu30, &b_HLT_IsoMu30);
   T1->SetBranchAddress("HLT_UncorrectedJetE30_NoBPTX", &HLT_UncorrectedJetE30_NoBPTX, &b_HLT_UncorrectedJetE30_NoBPTX);
   T1->SetBranchAddress("HLT_UncorrectedJetE30_NoBPTX3BX", &HLT_UncorrectedJetE30_NoBPTX3BX, &b_HLT_UncorrectedJetE30_NoBPTX3BX);
   T1->SetBranchAddress("HLT_UncorrectedJetE60_NoBPTX3BX", &HLT_UncorrectedJetE60_NoBPTX3BX, &b_HLT_UncorrectedJetE60_NoBPTX3BX);
   T1->SetBranchAddress("HLT_UncorrectedJetE70_NoBPTX3BX", &HLT_UncorrectedJetE70_NoBPTX3BX, &b_HLT_UncorrectedJetE70_NoBPTX3BX);
   T1->SetBranchAddress("HLT_L1SingleMu18", &HLT_L1SingleMu18, &b_HLT_L1SingleMu18);
   T1->SetBranchAddress("HLT_L1SingleMu25", &HLT_L1SingleMu25, &b_HLT_L1SingleMu25);
   T1->SetBranchAddress("HLT_L2Mu10", &HLT_L2Mu10, &b_HLT_L2Mu10);
   T1->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX3BX", &HLT_L2Mu10_NoVertex_NoBPTX3BX, &b_HLT_L2Mu10_NoVertex_NoBPTX3BX);
   T1->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX", &HLT_L2Mu10_NoVertex_NoBPTX, &b_HLT_L2Mu10_NoVertex_NoBPTX);
   T1->SetBranchAddress("HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX, &b_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX);
   T1->SetBranchAddress("HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX, &b_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX);
   T1->SetBranchAddress("HLT_L2Mu50", &HLT_L2Mu50, &b_HLT_L2Mu50);
   T1->SetBranchAddress("HLT_DoubleL2Mu50", &HLT_DoubleL2Mu50, &b_HLT_DoubleL2Mu50);
   T1->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
   T1->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL);
   T1->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
   T1->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ);
   T1->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8);
   T1->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8);
   T1->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8);
   T1->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8);
   T1->SetBranchAddress("HLT_Mu25_TkMu0_Onia", &HLT_Mu25_TkMu0_Onia, &b_HLT_Mu25_TkMu0_Onia);
   T1->SetBranchAddress("HLT_Mu30_TkMu0_Onia", &HLT_Mu30_TkMu0_Onia, &b_HLT_Mu30_TkMu0_Onia);
   T1->SetBranchAddress("HLT_Mu20_TkMu0_Phi", &HLT_Mu20_TkMu0_Phi, &b_HLT_Mu20_TkMu0_Phi);
   T1->SetBranchAddress("HLT_Mu25_TkMu0_Phi", &HLT_Mu25_TkMu0_Phi, &b_HLT_Mu25_TkMu0_Phi);
   T1->SetBranchAddress("HLT_Mu20", &HLT_Mu20, &b_HLT_Mu20);
   T1->SetBranchAddress("HLT_Mu27", &HLT_Mu27, &b_HLT_Mu27);
   T1->SetBranchAddress("HLT_Mu50", &HLT_Mu50, &b_HLT_Mu50);
   T1->SetBranchAddress("HLT_Mu55", &HLT_Mu55, &b_HLT_Mu55);
   T1->SetBranchAddress("HLT_OldMu100", &HLT_OldMu100, &b_HLT_OldMu100);
   T1->SetBranchAddress("HLT_TkMu100", &HLT_TkMu100, &b_HLT_TkMu100);
   T1->SetBranchAddress("HLT_DiPFJet15_NoCaloMatched", &HLT_DiPFJet15_NoCaloMatched, &b_HLT_DiPFJet15_NoCaloMatched);
   T1->SetBranchAddress("HLT_DiPFJet25_NoCaloMatched", &HLT_DiPFJet25_NoCaloMatched, &b_HLT_DiPFJet25_NoCaloMatched);
   T1->SetBranchAddress("HLT_DiPFJet15_FBEta3_NoCaloMatched", &HLT_DiPFJet15_FBEta3_NoCaloMatched, &b_HLT_DiPFJet15_FBEta3_NoCaloMatched);
   T1->SetBranchAddress("HLT_DiPFJet25_FBEta3_NoCaloMatched", &HLT_DiPFJet25_FBEta3_NoCaloMatched, &b_HLT_DiPFJet25_FBEta3_NoCaloMatched);
   T1->SetBranchAddress("HLT_DiPFJetAve40", &HLT_DiPFJetAve40, &b_HLT_DiPFJetAve40);
   T1->SetBranchAddress("HLT_DiPFJetAve60", &HLT_DiPFJetAve60, &b_HLT_DiPFJetAve60);
   T1->SetBranchAddress("HLT_DiPFJetAve80", &HLT_DiPFJetAve80, &b_HLT_DiPFJetAve80);
   T1->SetBranchAddress("HLT_DiPFJetAve140", &HLT_DiPFJetAve140, &b_HLT_DiPFJetAve140);
   T1->SetBranchAddress("HLT_DiPFJetAve200", &HLT_DiPFJetAve200, &b_HLT_DiPFJetAve200);
   T1->SetBranchAddress("HLT_DiPFJetAve260", &HLT_DiPFJetAve260, &b_HLT_DiPFJetAve260);
   T1->SetBranchAddress("HLT_DiPFJetAve320", &HLT_DiPFJetAve320, &b_HLT_DiPFJetAve320);
   T1->SetBranchAddress("HLT_DiPFJetAve400", &HLT_DiPFJetAve400, &b_HLT_DiPFJetAve400);
   T1->SetBranchAddress("HLT_DiPFJetAve500", &HLT_DiPFJetAve500, &b_HLT_DiPFJetAve500);
   T1->SetBranchAddress("HLT_DiPFJetAve15_HFJEC", &HLT_DiPFJetAve15_HFJEC, &b_HLT_DiPFJetAve15_HFJEC);
   T1->SetBranchAddress("HLT_DiPFJetAve25_HFJEC", &HLT_DiPFJetAve25_HFJEC, &b_HLT_DiPFJetAve25_HFJEC);
   T1->SetBranchAddress("HLT_DiPFJetAve35_HFJEC", &HLT_DiPFJetAve35_HFJEC, &b_HLT_DiPFJetAve35_HFJEC);
   T1->SetBranchAddress("HLT_DiPFJetAve60_HFJEC", &HLT_DiPFJetAve60_HFJEC, &b_HLT_DiPFJetAve60_HFJEC);
   T1->SetBranchAddress("HLT_DiPFJetAve80_HFJEC", &HLT_DiPFJetAve80_HFJEC, &b_HLT_DiPFJetAve80_HFJEC);
   T1->SetBranchAddress("HLT_DiPFJetAve100_HFJEC", &HLT_DiPFJetAve100_HFJEC, &b_HLT_DiPFJetAve100_HFJEC);
   T1->SetBranchAddress("HLT_DiPFJetAve160_HFJEC", &HLT_DiPFJetAve160_HFJEC, &b_HLT_DiPFJetAve160_HFJEC);
   T1->SetBranchAddress("HLT_DiPFJetAve220_HFJEC", &HLT_DiPFJetAve220_HFJEC, &b_HLT_DiPFJetAve220_HFJEC);
   T1->SetBranchAddress("HLT_DiPFJetAve300_HFJEC", &HLT_DiPFJetAve300_HFJEC, &b_HLT_DiPFJetAve300_HFJEC);
   T1->SetBranchAddress("HLT_AK8PFJet40", &HLT_AK8PFJet40, &b_HLT_AK8PFJet40);
   T1->SetBranchAddress("HLT_AK8PFJet60", &HLT_AK8PFJet60, &b_HLT_AK8PFJet60);
   T1->SetBranchAddress("HLT_AK8PFJet80", &HLT_AK8PFJet80, &b_HLT_AK8PFJet80);
   T1->SetBranchAddress("HLT_AK8PFJet140", &HLT_AK8PFJet140, &b_HLT_AK8PFJet140);
   T1->SetBranchAddress("HLT_AK8PFJet200", &HLT_AK8PFJet200, &b_HLT_AK8PFJet200);
   T1->SetBranchAddress("HLT_AK8PFJet260", &HLT_AK8PFJet260, &b_HLT_AK8PFJet260);
   T1->SetBranchAddress("HLT_AK8PFJet320", &HLT_AK8PFJet320, &b_HLT_AK8PFJet320);
   T1->SetBranchAddress("HLT_AK8PFJet400", &HLT_AK8PFJet400, &b_HLT_AK8PFJet400);
   T1->SetBranchAddress("HLT_AK8PFJet450", &HLT_AK8PFJet450, &b_HLT_AK8PFJet450);
   T1->SetBranchAddress("HLT_AK8PFJet500", &HLT_AK8PFJet500, &b_HLT_AK8PFJet500);
   T1->SetBranchAddress("HLT_AK8PFJet550", &HLT_AK8PFJet550, &b_HLT_AK8PFJet550);
   T1->SetBranchAddress("HLT_PFJet40", &HLT_PFJet40, &b_HLT_PFJet40);
   T1->SetBranchAddress("HLT_PFJet60", &HLT_PFJet60, &b_HLT_PFJet60);
   T1->SetBranchAddress("HLT_PFJet80", &HLT_PFJet80, &b_HLT_PFJet80);
   T1->SetBranchAddress("HLT_PFJet140", &HLT_PFJet140, &b_HLT_PFJet140);
   T1->SetBranchAddress("HLT_PFJet200", &HLT_PFJet200, &b_HLT_PFJet200);
   T1->SetBranchAddress("HLT_PFJet260", &HLT_PFJet260, &b_HLT_PFJet260);
   T1->SetBranchAddress("HLT_PFJet320", &HLT_PFJet320, &b_HLT_PFJet320);
   T1->SetBranchAddress("HLT_PFJet400", &HLT_PFJet400, &b_HLT_PFJet400);
   T1->SetBranchAddress("HLT_PFJet450", &HLT_PFJet450, &b_HLT_PFJet450);
   T1->SetBranchAddress("HLT_PFJet500", &HLT_PFJet500, &b_HLT_PFJet500);
   T1->SetBranchAddress("HLT_PFJet550", &HLT_PFJet550, &b_HLT_PFJet550);
   T1->SetBranchAddress("HLT_PFJetFwd40", &HLT_PFJetFwd40, &b_HLT_PFJetFwd40);
   T1->SetBranchAddress("HLT_PFJetFwd60", &HLT_PFJetFwd60, &b_HLT_PFJetFwd60);
   T1->SetBranchAddress("HLT_PFJetFwd80", &HLT_PFJetFwd80, &b_HLT_PFJetFwd80);
   T1->SetBranchAddress("HLT_PFJetFwd140", &HLT_PFJetFwd140, &b_HLT_PFJetFwd140);
   T1->SetBranchAddress("HLT_PFJetFwd200", &HLT_PFJetFwd200, &b_HLT_PFJetFwd200);
   T1->SetBranchAddress("HLT_PFJetFwd260", &HLT_PFJetFwd260, &b_HLT_PFJetFwd260);
   T1->SetBranchAddress("HLT_PFJetFwd320", &HLT_PFJetFwd320, &b_HLT_PFJetFwd320);
   T1->SetBranchAddress("HLT_PFJetFwd400", &HLT_PFJetFwd400, &b_HLT_PFJetFwd400);
   T1->SetBranchAddress("HLT_PFJetFwd450", &HLT_PFJetFwd450, &b_HLT_PFJetFwd450);
   T1->SetBranchAddress("HLT_PFJetFwd500", &HLT_PFJetFwd500, &b_HLT_PFJetFwd500);
   T1->SetBranchAddress("HLT_AK8PFJetFwd40", &HLT_AK8PFJetFwd40, &b_HLT_AK8PFJetFwd40);
   T1->SetBranchAddress("HLT_AK8PFJetFwd60", &HLT_AK8PFJetFwd60, &b_HLT_AK8PFJetFwd60);
   T1->SetBranchAddress("HLT_AK8PFJetFwd80", &HLT_AK8PFJetFwd80, &b_HLT_AK8PFJetFwd80);
   T1->SetBranchAddress("HLT_AK8PFJetFwd140", &HLT_AK8PFJetFwd140, &b_HLT_AK8PFJetFwd140);
   T1->SetBranchAddress("HLT_AK8PFJetFwd200", &HLT_AK8PFJetFwd200, &b_HLT_AK8PFJetFwd200);
   T1->SetBranchAddress("HLT_AK8PFJetFwd260", &HLT_AK8PFJetFwd260, &b_HLT_AK8PFJetFwd260);
   T1->SetBranchAddress("HLT_AK8PFJetFwd320", &HLT_AK8PFJetFwd320, &b_HLT_AK8PFJetFwd320);
   T1->SetBranchAddress("HLT_AK8PFJetFwd400", &HLT_AK8PFJetFwd400, &b_HLT_AK8PFJetFwd400);
   T1->SetBranchAddress("HLT_AK8PFJetFwd450", &HLT_AK8PFJetFwd450, &b_HLT_AK8PFJetFwd450);
   T1->SetBranchAddress("HLT_AK8PFJetFwd500", &HLT_AK8PFJetFwd500, &b_HLT_AK8PFJetFwd500);
   T1->SetBranchAddress("HLT_PFHT180", &HLT_PFHT180, &b_HLT_PFHT180);
   T1->SetBranchAddress("HLT_PFHT250", &HLT_PFHT250, &b_HLT_PFHT250);
   T1->SetBranchAddress("HLT_PFHT370", &HLT_PFHT370, &b_HLT_PFHT370);
   T1->SetBranchAddress("HLT_PFHT430", &HLT_PFHT430, &b_HLT_PFHT430);
   T1->SetBranchAddress("HLT_PFHT510", &HLT_PFHT510, &b_HLT_PFHT510);
   T1->SetBranchAddress("HLT_PFHT590", &HLT_PFHT590, &b_HLT_PFHT590);
   T1->SetBranchAddress("HLT_PFHT680", &HLT_PFHT680, &b_HLT_PFHT680);
   T1->SetBranchAddress("HLT_PFHT780", &HLT_PFHT780, &b_HLT_PFHT780);
   T1->SetBranchAddress("HLT_PFHT890", &HLT_PFHT890, &b_HLT_PFHT890);
   T1->SetBranchAddress("HLT_PFHT1050", &HLT_PFHT1050, &b_HLT_PFHT1050);
   T1->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight", &HLT_PFHT500_PFMET100_PFMHT100_IDTight, &b_HLT_PFHT500_PFMET100_PFMHT100_IDTight);
   T1->SetBranchAddress("HLT_PFHT500_PFMET110_PFMHT110_IDTight", &HLT_PFHT500_PFMET110_PFMHT110_IDTight, &b_HLT_PFHT500_PFMET110_PFMHT110_IDTight);
   T1->SetBranchAddress("HLT_PFHT700_PFMET85_PFMHT85_IDTight", &HLT_PFHT700_PFMET85_PFMHT85_IDTight, &b_HLT_PFHT700_PFMET85_PFMHT85_IDTight);
   T1->SetBranchAddress("HLT_PFHT700_PFMET95_PFMHT95_IDTight", &HLT_PFHT700_PFMET95_PFMHT95_IDTight, &b_HLT_PFHT700_PFMET95_PFMHT95_IDTight);
   T1->SetBranchAddress("HLT_PFHT800_PFMET75_PFMHT75_IDTight", &HLT_PFHT800_PFMET75_PFMHT75_IDTight, &b_HLT_PFHT800_PFMET75_PFMHT75_IDTight);
   T1->SetBranchAddress("HLT_PFHT800_PFMET85_PFMHT85_IDTight", &HLT_PFHT800_PFMET85_PFMHT85_IDTight, &b_HLT_PFHT800_PFMET85_PFMHT85_IDTight);
   T1->SetBranchAddress("HLT_PFMET110_PFMHT110_IDTight", &HLT_PFMET110_PFMHT110_IDTight, &b_HLT_PFMET110_PFMHT110_IDTight);
   T1->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight", &HLT_PFMET120_PFMHT120_IDTight, &b_HLT_PFMET120_PFMHT120_IDTight);
   T1->SetBranchAddress("HLT_PFMET130_PFMHT130_IDTight", &HLT_PFMET130_PFMHT130_IDTight, &b_HLT_PFMET130_PFMHT130_IDTight);
   T1->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight", &HLT_PFMET140_PFMHT140_IDTight, &b_HLT_PFMET140_PFMHT140_IDTight);
   T1->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1", &HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1, &b_HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1);
   T1->SetBranchAddress("HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1", &HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1, &b_HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1);
   T1->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1", &HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1, &b_HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1);
   T1->SetBranchAddress("HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1", &HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1, &b_HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1);
   T1->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1", &HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1, &b_HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1);
   T1->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_PFHT60", &HLT_PFMET120_PFMHT120_IDTight_PFHT60, &b_HLT_PFMET120_PFMHT120_IDTight_PFHT60);
   T1->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
   T1->SetBranchAddress("HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60", &HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60, &b_HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60);
   T1->SetBranchAddress("HLT_PFMETTypeOne110_PFMHT110_IDTight", &HLT_PFMETTypeOne110_PFMHT110_IDTight, &b_HLT_PFMETTypeOne110_PFMHT110_IDTight);
   T1->SetBranchAddress("HLT_PFMETTypeOne120_PFMHT120_IDTight", &HLT_PFMETTypeOne120_PFMHT120_IDTight, &b_HLT_PFMETTypeOne120_PFMHT120_IDTight);
   T1->SetBranchAddress("HLT_PFMETTypeOne130_PFMHT130_IDTight", &HLT_PFMETTypeOne130_PFMHT130_IDTight, &b_HLT_PFMETTypeOne130_PFMHT130_IDTight);
   T1->SetBranchAddress("HLT_PFMETTypeOne140_PFMHT140_IDTight", &HLT_PFMETTypeOne140_PFMHT140_IDTight, &b_HLT_PFMETTypeOne140_PFMHT140_IDTight);
   T1->SetBranchAddress("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_PFMETNoMu110_PFMHTNoMu110_IDTight, &b_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight);
   T1->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
   T1->SetBranchAddress("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight", &HLT_PFMETNoMu130_PFMHTNoMu130_IDTight, &b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight);
   T1->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight, &b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight);
   T1->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight);
   T1->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight);
   T1->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight);
   T1->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight);
   T1->SetBranchAddress("HLT_L1ETMHadSeeds", &HLT_L1ETMHadSeeds, &b_HLT_L1ETMHadSeeds);
   T1->SetBranchAddress("HLT_CaloMHT90", &HLT_CaloMHT90, &b_HLT_CaloMHT90);
   T1->SetBranchAddress("HLT_CaloMET80_NotCleaned", &HLT_CaloMET80_NotCleaned, &b_HLT_CaloMET80_NotCleaned);
   T1->SetBranchAddress("HLT_CaloMET90_NotCleaned", &HLT_CaloMET90_NotCleaned, &b_HLT_CaloMET90_NotCleaned);
   T1->SetBranchAddress("HLT_CaloMET100_NotCleaned", &HLT_CaloMET100_NotCleaned, &b_HLT_CaloMET100_NotCleaned);
   T1->SetBranchAddress("HLT_CaloMET110_NotCleaned", &HLT_CaloMET110_NotCleaned, &b_HLT_CaloMET110_NotCleaned);
   T1->SetBranchAddress("HLT_CaloMET250_NotCleaned", &HLT_CaloMET250_NotCleaned, &b_HLT_CaloMET250_NotCleaned);
   T1->SetBranchAddress("HLT_CaloMET70_HBHECleaned", &HLT_CaloMET70_HBHECleaned, &b_HLT_CaloMET70_HBHECleaned);
   T1->SetBranchAddress("HLT_CaloMET80_HBHECleaned", &HLT_CaloMET80_HBHECleaned, &b_HLT_CaloMET80_HBHECleaned);
   T1->SetBranchAddress("HLT_CaloMET90_HBHECleaned", &HLT_CaloMET90_HBHECleaned, &b_HLT_CaloMET90_HBHECleaned);
   T1->SetBranchAddress("HLT_CaloMET100_HBHECleaned", &HLT_CaloMET100_HBHECleaned, &b_HLT_CaloMET100_HBHECleaned);
   T1->SetBranchAddress("HLT_CaloMET250_HBHECleaned", &HLT_CaloMET250_HBHECleaned, &b_HLT_CaloMET250_HBHECleaned);
   T1->SetBranchAddress("HLT_CaloMET300_HBHECleaned", &HLT_CaloMET300_HBHECleaned, &b_HLT_CaloMET300_HBHECleaned);
   T1->SetBranchAddress("HLT_CaloMET350_HBHECleaned", &HLT_CaloMET350_HBHECleaned, &b_HLT_CaloMET350_HBHECleaned);
   T1->SetBranchAddress("HLT_PFMET200_NotCleaned", &HLT_PFMET200_NotCleaned, &b_HLT_PFMET200_NotCleaned);
   T1->SetBranchAddress("HLT_PFMET200_HBHECleaned", &HLT_PFMET200_HBHECleaned, &b_HLT_PFMET200_HBHECleaned);
   T1->SetBranchAddress("HLT_PFMET250_HBHECleaned", &HLT_PFMET250_HBHECleaned, &b_HLT_PFMET250_HBHECleaned);
   T1->SetBranchAddress("HLT_PFMET300_HBHECleaned", &HLT_PFMET300_HBHECleaned, &b_HLT_PFMET300_HBHECleaned);
   T1->SetBranchAddress("HLT_PFMET200_HBHE_BeamHaloCleaned", &HLT_PFMET200_HBHE_BeamHaloCleaned, &b_HLT_PFMET200_HBHE_BeamHaloCleaned);
   T1->SetBranchAddress("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned", &HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned, &b_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned);
   T1->SetBranchAddress("HLT_MET105_IsoTrk50", &HLT_MET105_IsoTrk50, &b_HLT_MET105_IsoTrk50);
   T1->SetBranchAddress("HLT_MET120_IsoTrk50", &HLT_MET120_IsoTrk50, &b_HLT_MET120_IsoTrk50);
   T1->SetBranchAddress("HLT_SingleJet30_Mu12_SinglePFJet40", &HLT_SingleJet30_Mu12_SinglePFJet40, &b_HLT_SingleJet30_Mu12_SinglePFJet40);
   T1->SetBranchAddress("HLT_Mu12_DoublePFJets40_CaloBTagCSV_p33", &HLT_Mu12_DoublePFJets40_CaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets40_CaloBTagCSV_p33);
   T1->SetBranchAddress("HLT_Mu12_DoublePFJets100_CaloBTagCSV_p33", &HLT_Mu12_DoublePFJets100_CaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets100_CaloBTagCSV_p33);
   T1->SetBranchAddress("HLT_Mu12_DoublePFJets200_CaloBTagCSV_p33", &HLT_Mu12_DoublePFJets200_CaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets200_CaloBTagCSV_p33);
   T1->SetBranchAddress("HLT_Mu12_DoublePFJets350_CaloBTagCSV_p33", &HLT_Mu12_DoublePFJets350_CaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets350_CaloBTagCSV_p33);
   T1->SetBranchAddress("HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33", &HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33);
   T1->SetBranchAddress("HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagCSV_p33", &HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagCSV_p33);
   T1->SetBranchAddress("HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagCSV_p33", &HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagCSV_p33);
   T1->SetBranchAddress("HLT_DoublePFJets40_CaloBTagCSV_p33", &HLT_DoublePFJets40_CaloBTagCSV_p33, &b_HLT_DoublePFJets40_CaloBTagCSV_p33);
   T1->SetBranchAddress("HLT_DoublePFJets100_CaloBTagCSV_p33", &HLT_DoublePFJets100_CaloBTagCSV_p33, &b_HLT_DoublePFJets100_CaloBTagCSV_p33);
   T1->SetBranchAddress("HLT_DoublePFJets200_CaloBTagCSV_p33", &HLT_DoublePFJets200_CaloBTagCSV_p33, &b_HLT_DoublePFJets200_CaloBTagCSV_p33);
   T1->SetBranchAddress("HLT_DoublePFJets350_CaloBTagCSV_p33", &HLT_DoublePFJets350_CaloBTagCSV_p33, &b_HLT_DoublePFJets350_CaloBTagCSV_p33);
   T1->SetBranchAddress("HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33", &HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33, &b_HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33);
   T1->SetBranchAddress("HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33", &HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33, &b_HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33);
   T1->SetBranchAddress("HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33", &HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33, &b_HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33);
   T1->SetBranchAddress("HLT_Photon300_NoHE", &HLT_Photon300_NoHE, &b_HLT_Photon300_NoHE);
   T1->SetBranchAddress("HLT_Mu8_TrkIsoVVL", &HLT_Mu8_TrkIsoVVL, &b_HLT_Mu8_TrkIsoVVL);
   T1->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ, &b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ);
   T1->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL, &b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL);
   T1->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ, &b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ);
   T1->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350, &b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350);
   T1->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
   T1->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
   T1->SetBranchAddress("HLT_Mu17_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL);
   T1->SetBranchAddress("HLT_Mu19_TrkIsoVVL", &HLT_Mu19_TrkIsoVVL, &b_HLT_Mu19_TrkIsoVVL);
   T1->SetBranchAddress("HLT_BTagMu_AK4DiJet20_Mu5", &HLT_BTagMu_AK4DiJet20_Mu5, &b_HLT_BTagMu_AK4DiJet20_Mu5);
   T1->SetBranchAddress("HLT_BTagMu_AK4DiJet40_Mu5", &HLT_BTagMu_AK4DiJet40_Mu5, &b_HLT_BTagMu_AK4DiJet40_Mu5);
   T1->SetBranchAddress("HLT_BTagMu_AK4DiJet70_Mu5", &HLT_BTagMu_AK4DiJet70_Mu5, &b_HLT_BTagMu_AK4DiJet70_Mu5);
   T1->SetBranchAddress("HLT_BTagMu_AK4DiJet110_Mu5", &HLT_BTagMu_AK4DiJet110_Mu5, &b_HLT_BTagMu_AK4DiJet110_Mu5);
   T1->SetBranchAddress("HLT_BTagMu_AK4DiJet170_Mu5", &HLT_BTagMu_AK4DiJet170_Mu5, &b_HLT_BTagMu_AK4DiJet170_Mu5);
   T1->SetBranchAddress("HLT_BTagMu_AK4Jet300_Mu5", &HLT_BTagMu_AK4Jet300_Mu5, &b_HLT_BTagMu_AK4Jet300_Mu5);
   T1->SetBranchAddress("HLT_BTagMu_AK8DiJet170_Mu5", &HLT_BTagMu_AK8DiJet170_Mu5, &b_HLT_BTagMu_AK8DiJet170_Mu5);
   T1->SetBranchAddress("HLT_BTagMu_AK8Jet300_Mu5", &HLT_BTagMu_AK8Jet300_Mu5, &b_HLT_BTagMu_AK8Jet300_Mu5);
   T1->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   T1->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
   T1->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   T1->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
   T1->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
   T1->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
   T1->SetBranchAddress("HLT_Mu12_DoublePhoton20", &HLT_Mu12_DoublePhoton20, &b_HLT_Mu12_DoublePhoton20);
   T1->SetBranchAddress("HLT_TriplePhoton_20_20_20_CaloIdLV2", &HLT_TriplePhoton_20_20_20_CaloIdLV2, &b_HLT_TriplePhoton_20_20_20_CaloIdLV2);
   T1->SetBranchAddress("HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL", &HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL, &b_HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL);
   T1->SetBranchAddress("HLT_TriplePhoton_30_30_10_CaloIdLV2", &HLT_TriplePhoton_30_30_10_CaloIdLV2, &b_HLT_TriplePhoton_30_30_10_CaloIdLV2);
   T1->SetBranchAddress("HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL", &HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL, &b_HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL);
   T1->SetBranchAddress("HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL", &HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL, &b_HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL);
   T1->SetBranchAddress("HLT_Photon25", &HLT_Photon25, &b_HLT_Photon25);
   T1->SetBranchAddress("HLT_Photon33", &HLT_Photon33, &b_HLT_Photon33);
   T1->SetBranchAddress("HLT_Photon50", &HLT_Photon50, &b_HLT_Photon50);
   T1->SetBranchAddress("HLT_Photon75", &HLT_Photon75, &b_HLT_Photon75);
   T1->SetBranchAddress("HLT_Photon90", &HLT_Photon90, &b_HLT_Photon90);
   T1->SetBranchAddress("HLT_Photon120", &HLT_Photon120, &b_HLT_Photon120);
   T1->SetBranchAddress("HLT_Photon150", &HLT_Photon150, &b_HLT_Photon150);
   T1->SetBranchAddress("HLT_Photon175", &HLT_Photon175, &b_HLT_Photon175);
   T1->SetBranchAddress("HLT_Photon200", &HLT_Photon200, &b_HLT_Photon200);
   T1->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM", &HLT_Photon50_R9Id90_HE10_IsoM, &b_HLT_Photon50_R9Id90_HE10_IsoM);
   T1->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM", &HLT_Photon75_R9Id90_HE10_IsoM, &b_HLT_Photon75_R9Id90_HE10_IsoM);
   T1->SetBranchAddress("HLT_Photon90_R9Id90_HE10_IsoM", &HLT_Photon90_R9Id90_HE10_IsoM, &b_HLT_Photon90_R9Id90_HE10_IsoM);
   T1->SetBranchAddress("HLT_Photon120_R9Id90_HE10_IsoM", &HLT_Photon120_R9Id90_HE10_IsoM, &b_HLT_Photon120_R9Id90_HE10_IsoM);
   T1->SetBranchAddress("HLT_Photon165_R9Id90_HE10_IsoM", &HLT_Photon165_R9Id90_HE10_IsoM, &b_HLT_Photon165_R9Id90_HE10_IsoM);
   T1->SetBranchAddress("HLT_Photon90_CaloIdL_PFHT700", &HLT_Photon90_CaloIdL_PFHT700, &b_HLT_Photon90_CaloIdL_PFHT700);
   T1->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", &HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90, &b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90);
   T1->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95", &HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95, &b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95);
   T1->SetBranchAddress("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55", &HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55, &b_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55);
   T1->SetBranchAddress("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55", &HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55, &b_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55);
   T1->SetBranchAddress("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55", &HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55, &b_HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55);
   T1->SetBranchAddress("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55", &HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55, &b_HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55);
   T1->SetBranchAddress("HLT_Dimuon0_Jpsi_L1_NoOS", &HLT_Dimuon0_Jpsi_L1_NoOS, &b_HLT_Dimuon0_Jpsi_L1_NoOS);
   T1->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing_NoOS", &HLT_Dimuon0_Jpsi_NoVertexing_NoOS, &b_HLT_Dimuon0_Jpsi_NoVertexing_NoOS);
   T1->SetBranchAddress("HLT_Dimuon0_Jpsi", &HLT_Dimuon0_Jpsi, &b_HLT_Dimuon0_Jpsi);
   T1->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing", &HLT_Dimuon0_Jpsi_NoVertexing, &b_HLT_Dimuon0_Jpsi_NoVertexing);
   T1->SetBranchAddress("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R", &HLT_Dimuon0_Jpsi_L1_4R_0er1p5R, &b_HLT_Dimuon0_Jpsi_L1_4R_0er1p5R);
   T1->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R", &HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R, &b_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R);
   T1->SetBranchAddress("HLT_Dimuon0_Jpsi3p5_Muon2", &HLT_Dimuon0_Jpsi3p5_Muon2, &b_HLT_Dimuon0_Jpsi3p5_Muon2);
   T1->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5", &HLT_Dimuon0_Upsilon_L1_4p5, &b_HLT_Dimuon0_Upsilon_L1_4p5);
   T1->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_5", &HLT_Dimuon0_Upsilon_L1_5, &b_HLT_Dimuon0_Upsilon_L1_5);
   T1->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5NoOS", &HLT_Dimuon0_Upsilon_L1_4p5NoOS, &b_HLT_Dimuon0_Upsilon_L1_4p5NoOS);
   T1->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5er2p0", &HLT_Dimuon0_Upsilon_L1_4p5er2p0, &b_HLT_Dimuon0_Upsilon_L1_4p5er2p0);
   T1->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5er2p0M", &HLT_Dimuon0_Upsilon_L1_4p5er2p0M, &b_HLT_Dimuon0_Upsilon_L1_4p5er2p0M);
   T1->SetBranchAddress("HLT_Dimuon0_Upsilon_NoVertexing", &HLT_Dimuon0_Upsilon_NoVertexing, &b_HLT_Dimuon0_Upsilon_NoVertexing);
   T1->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_5M", &HLT_Dimuon0_Upsilon_L1_5M, &b_HLT_Dimuon0_Upsilon_L1_5M);
   T1->SetBranchAddress("HLT_Dimuon0_LowMass_L1_0er1p5R", &HLT_Dimuon0_LowMass_L1_0er1p5R, &b_HLT_Dimuon0_LowMass_L1_0er1p5R);
   T1->SetBranchAddress("HLT_Dimuon0_LowMass_L1_0er1p5", &HLT_Dimuon0_LowMass_L1_0er1p5, &b_HLT_Dimuon0_LowMass_L1_0er1p5);
   T1->SetBranchAddress("HLT_Dimuon0_LowMass", &HLT_Dimuon0_LowMass, &b_HLT_Dimuon0_LowMass);
   T1->SetBranchAddress("HLT_Dimuon0_LowMass_L1_4", &HLT_Dimuon0_LowMass_L1_4, &b_HLT_Dimuon0_LowMass_L1_4);
   T1->SetBranchAddress("HLT_Dimuon0_LowMass_L1_4R", &HLT_Dimuon0_LowMass_L1_4R, &b_HLT_Dimuon0_LowMass_L1_4R);
   T1->SetBranchAddress("HLT_Dimuon0_LowMass_L1_TM530", &HLT_Dimuon0_LowMass_L1_TM530, &b_HLT_Dimuon0_LowMass_L1_TM530);
   T1->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_L1_TM0", &HLT_Dimuon0_Upsilon_Muon_L1_TM0, &b_HLT_Dimuon0_Upsilon_Muon_L1_TM0);
   T1->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_NoL1Mass", &HLT_Dimuon0_Upsilon_Muon_NoL1Mass, &b_HLT_Dimuon0_Upsilon_Muon_NoL1Mass);
   T1->SetBranchAddress("HLT_TripleMu_5_3_3_Mass3p8to60_DZ", &HLT_TripleMu_5_3_3_Mass3p8to60_DZ, &b_HLT_TripleMu_5_3_3_Mass3p8to60_DZ);
   T1->SetBranchAddress("HLT_TripleMu_10_5_5_DZ", &HLT_TripleMu_10_5_5_DZ, &b_HLT_TripleMu_10_5_5_DZ);
   T1->SetBranchAddress("HLT_TripleMu_12_10_5", &HLT_TripleMu_12_10_5, &b_HLT_TripleMu_12_10_5);
   T1->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15);
   T1->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1);
   T1->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15);
   T1->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1);
   T1->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET50_PFMHT60", &HLT_DoubleMu3_DZ_PFMET50_PFMHT60, &b_HLT_DoubleMu3_DZ_PFMET50_PFMHT60);
   T1->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET70_PFMHT70", &HLT_DoubleMu3_DZ_PFMET70_PFMHT70, &b_HLT_DoubleMu3_DZ_PFMET70_PFMHT70);
   T1->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET90_PFMHT90", &HLT_DoubleMu3_DZ_PFMET90_PFMHT90, &b_HLT_DoubleMu3_DZ_PFMET90_PFMHT90);
   T1->SetBranchAddress("HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass", &HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass, &b_HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass);
   T1->SetBranchAddress("HLT_DoubleMu4_Jpsi_Displaced", &HLT_DoubleMu4_Jpsi_Displaced, &b_HLT_DoubleMu4_Jpsi_Displaced);
   T1->SetBranchAddress("HLT_DoubleMu4_Jpsi_NoVertexing", &HLT_DoubleMu4_Jpsi_NoVertexing, &b_HLT_DoubleMu4_Jpsi_NoVertexing);
   T1->SetBranchAddress("HLT_DoubleMu4_JpsiTrkTrk_Displaced", &HLT_DoubleMu4_JpsiTrkTrk_Displaced, &b_HLT_DoubleMu4_JpsiTrkTrk_Displaced);
   T1->SetBranchAddress("HLT_DoubleMu43NoFiltersNoVtx", &HLT_DoubleMu43NoFiltersNoVtx, &b_HLT_DoubleMu43NoFiltersNoVtx);
   T1->SetBranchAddress("HLT_DoubleMu48NoFiltersNoVtx", &HLT_DoubleMu48NoFiltersNoVtx, &b_HLT_DoubleMu48NoFiltersNoVtx);
   T1->SetBranchAddress("HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL", &HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL, &b_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL);
   T1->SetBranchAddress("HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL", &HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL, &b_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL);
   T1->SetBranchAddress("HLT_DoubleMu20_7_Mass0to30_L1_DM4", &HLT_DoubleMu20_7_Mass0to30_L1_DM4, &b_HLT_DoubleMu20_7_Mass0to30_L1_DM4);
   T1->SetBranchAddress("HLT_DoubleMu20_7_Mass0to30_L1_DM4EG", &HLT_DoubleMu20_7_Mass0to30_L1_DM4EG, &b_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG);
   T1->SetBranchAddress("HLT_HT425", &HLT_HT425, &b_HLT_HT425);
   T1->SetBranchAddress("HLT_HT430_DisplacedDijet40_DisplacedTrack", &HLT_HT430_DisplacedDijet40_DisplacedTrack, &b_HLT_HT430_DisplacedDijet40_DisplacedTrack);
   T1->SetBranchAddress("HLT_HT430_DisplacedDijet60_DisplacedTrack", &HLT_HT430_DisplacedDijet60_DisplacedTrack, &b_HLT_HT430_DisplacedDijet60_DisplacedTrack);
   T1->SetBranchAddress("HLT_HT430_DisplacedDijet80_DisplacedTrack", &HLT_HT430_DisplacedDijet80_DisplacedTrack, &b_HLT_HT430_DisplacedDijet80_DisplacedTrack);
   T1->SetBranchAddress("HLT_HT400_DisplacedDijet40_DisplacedTrack", &HLT_HT400_DisplacedDijet40_DisplacedTrack, &b_HLT_HT400_DisplacedDijet40_DisplacedTrack);
   T1->SetBranchAddress("HLT_HT650_DisplacedDijet60_Inclusive", &HLT_HT650_DisplacedDijet60_Inclusive, &b_HLT_HT650_DisplacedDijet60_Inclusive);
   T1->SetBranchAddress("HLT_HT550_DisplacedDijet80_Inclusive", &HLT_HT550_DisplacedDijet80_Inclusive, &b_HLT_HT550_DisplacedDijet80_Inclusive);
   T1->SetBranchAddress("HLT_HT550_DisplacedDijet60_Inclusive", &HLT_HT550_DisplacedDijet60_Inclusive, &b_HLT_HT550_DisplacedDijet60_Inclusive);
   T1->SetBranchAddress("HLT_HT650_DisplacedDijet80_Inclusive", &HLT_HT650_DisplacedDijet80_Inclusive, &b_HLT_HT650_DisplacedDijet80_Inclusive);
   T1->SetBranchAddress("HLT_HT750_DisplacedDijet80_Inclusive", &HLT_HT750_DisplacedDijet80_Inclusive, &b_HLT_HT750_DisplacedDijet80_Inclusive);
   T1->SetBranchAddress("HLT_DiJet110_35_Mjj650_PFMET110", &HLT_DiJet110_35_Mjj650_PFMET110, &b_HLT_DiJet110_35_Mjj650_PFMET110);
   T1->SetBranchAddress("HLT_DiJet110_35_Mjj650_PFMET120", &HLT_DiJet110_35_Mjj650_PFMET120, &b_HLT_DiJet110_35_Mjj650_PFMET120);
   T1->SetBranchAddress("HLT_DiJet110_35_Mjj650_PFMET130", &HLT_DiJet110_35_Mjj650_PFMET130, &b_HLT_DiJet110_35_Mjj650_PFMET130);
   T1->SetBranchAddress("HLT_TripleJet110_35_35_Mjj650_PFMET110", &HLT_TripleJet110_35_35_Mjj650_PFMET110, &b_HLT_TripleJet110_35_35_Mjj650_PFMET110);
   T1->SetBranchAddress("HLT_TripleJet110_35_35_Mjj650_PFMET120", &HLT_TripleJet110_35_35_Mjj650_PFMET120, &b_HLT_TripleJet110_35_35_Mjj650_PFMET120);
   T1->SetBranchAddress("HLT_TripleJet110_35_35_Mjj650_PFMET130", &HLT_TripleJet110_35_35_Mjj650_PFMET130, &b_HLT_TripleJet110_35_35_Mjj650_PFMET130);
   T1->SetBranchAddress("HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg", &HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg, &b_HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg);
   T1->SetBranchAddress("HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg", &HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg, &b_HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg);
   T1->SetBranchAddress("HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg", &HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg, &b_HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg);
   T1->SetBranchAddress("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned", &HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned, &b_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned);
   T1->SetBranchAddress("HLT_Ele28_eta2p1_WPTight_Gsf_HT150", &HLT_Ele28_eta2p1_WPTight_Gsf_HT150, &b_HLT_Ele28_eta2p1_WPTight_Gsf_HT150);
   T1->SetBranchAddress("HLT_Ele28_HighEta_SC20_Mass55", &HLT_Ele28_HighEta_SC20_Mass55, &b_HLT_Ele28_HighEta_SC20_Mass55);
   T1->SetBranchAddress("HLT_DoubleMu20_7_Mass0to30_Photon23", &HLT_DoubleMu20_7_Mass0to30_Photon23, &b_HLT_DoubleMu20_7_Mass0to30_Photon23);
   T1->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5", &HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5, &b_HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5);
   T1->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450_PFMET50", &HLT_Ele15_IsoVVVL_PFHT450_PFMET50, &b_HLT_Ele15_IsoVVVL_PFHT450_PFMET50);
   T1->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450", &HLT_Ele15_IsoVVVL_PFHT450, &b_HLT_Ele15_IsoVVVL_PFHT450);
   T1->SetBranchAddress("HLT_Ele50_IsoVVVL_PFHT450", &HLT_Ele50_IsoVVVL_PFHT450, &b_HLT_Ele50_IsoVVVL_PFHT450);
   T1->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT600", &HLT_Ele15_IsoVVVL_PFHT600, &b_HLT_Ele15_IsoVVVL_PFHT600);
   T1->SetBranchAddress("HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60", &HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60, &b_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60);
   T1->SetBranchAddress("HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60", &HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60, &b_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60);
   T1->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5", &HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5, &b_HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5);
   T1->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450_PFMET50", &HLT_Mu15_IsoVVVL_PFHT450_PFMET50, &b_HLT_Mu15_IsoVVVL_PFHT450_PFMET50);
   T1->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450", &HLT_Mu15_IsoVVVL_PFHT450, &b_HLT_Mu15_IsoVVVL_PFHT450);
   T1->SetBranchAddress("HLT_Mu50_IsoVVVL_PFHT450", &HLT_Mu50_IsoVVVL_PFHT450, &b_HLT_Mu50_IsoVVVL_PFHT450);
   T1->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT600", &HLT_Mu15_IsoVVVL_PFHT600, &b_HLT_Mu15_IsoVVVL_PFHT600);
   T1->SetBranchAddress("HLT_Dimuon10_PsiPrime_Barrel_Seagulls", &HLT_Dimuon10_PsiPrime_Barrel_Seagulls, &b_HLT_Dimuon10_PsiPrime_Barrel_Seagulls);
   T1->SetBranchAddress("HLT_Dimuon20_Jpsi_Barrel_Seagulls", &HLT_Dimuon20_Jpsi_Barrel_Seagulls, &b_HLT_Dimuon20_Jpsi_Barrel_Seagulls);
   T1->SetBranchAddress("HLT_Dimuon10_Upsilon_Barrel_Seagulls", &HLT_Dimuon10_Upsilon_Barrel_Seagulls, &b_HLT_Dimuon10_Upsilon_Barrel_Seagulls);
   T1->SetBranchAddress("HLT_Dimuon12_Upsilon_eta1p5", &HLT_Dimuon12_Upsilon_eta1p5, &b_HLT_Dimuon12_Upsilon_eta1p5);
   T1->SetBranchAddress("HLT_Dimuon14_Phi_Barrel_Seagulls", &HLT_Dimuon14_Phi_Barrel_Seagulls, &b_HLT_Dimuon14_Phi_Barrel_Seagulls);
   T1->SetBranchAddress("HLT_Dimuon18_PsiPrime", &HLT_Dimuon18_PsiPrime, &b_HLT_Dimuon18_PsiPrime);
   T1->SetBranchAddress("HLT_Dimuon25_Jpsi", &HLT_Dimuon25_Jpsi, &b_HLT_Dimuon25_Jpsi);
   T1->SetBranchAddress("HLT_Dimuon18_PsiPrime_noCorrL1", &HLT_Dimuon18_PsiPrime_noCorrL1, &b_HLT_Dimuon18_PsiPrime_noCorrL1);
   T1->SetBranchAddress("HLT_Dimuon24_Upsilon_noCorrL1", &HLT_Dimuon24_Upsilon_noCorrL1, &b_HLT_Dimuon24_Upsilon_noCorrL1);
   T1->SetBranchAddress("HLT_Dimuon24_Phi_noCorrL1", &HLT_Dimuon24_Phi_noCorrL1, &b_HLT_Dimuon24_Phi_noCorrL1);
   T1->SetBranchAddress("HLT_Dimuon25_Jpsi_noCorrL1", &HLT_Dimuon25_Jpsi_noCorrL1, &b_HLT_Dimuon25_Jpsi_noCorrL1);
   T1->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ, &b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ);
   T1->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL, &b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL);
   T1->SetBranchAddress("HLT_DoubleIsoMu20_eta2p1", &HLT_DoubleIsoMu20_eta2p1, &b_HLT_DoubleIsoMu20_eta2p1);
   T1->SetBranchAddress("HLT_DoubleIsoMu24_eta2p1", &HLT_DoubleIsoMu24_eta2p1, &b_HLT_DoubleIsoMu24_eta2p1);
   T1->SetBranchAddress("HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx", &HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx, &b_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx);
   T1->SetBranchAddress("HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx", &HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx, &b_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx);
   T1->SetBranchAddress("HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx", &HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx, &b_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx);
   T1->SetBranchAddress("HLT_Mu8", &HLT_Mu8, &b_HLT_Mu8);
   T1->SetBranchAddress("HLT_Mu17", &HLT_Mu17, &b_HLT_Mu17);
   T1->SetBranchAddress("HLT_Mu19", &HLT_Mu19, &b_HLT_Mu19);
   T1->SetBranchAddress("HLT_Mu17_Photon30_IsoCaloId", &HLT_Mu17_Photon30_IsoCaloId, &b_HLT_Mu17_Photon30_IsoCaloId);
   T1->SetBranchAddress("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30);
   T1->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30);
   T1->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30);
   T1->SetBranchAddress("HLT_Ele8_CaloIdM_TrackIdM_PFJet30", &HLT_Ele8_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele8_CaloIdM_TrackIdM_PFJet30);
   T1->SetBranchAddress("HLT_Ele17_CaloIdM_TrackIdM_PFJet30", &HLT_Ele17_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele17_CaloIdM_TrackIdM_PFJet30);
   T1->SetBranchAddress("HLT_Ele23_CaloIdM_TrackIdM_PFJet30", &HLT_Ele23_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele23_CaloIdM_TrackIdM_PFJet30);
   T1->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", &HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165, &b_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
   T1->SetBranchAddress("HLT_Ele115_CaloIdVT_GsfTrkIdT", &HLT_Ele115_CaloIdVT_GsfTrkIdT, &b_HLT_Ele115_CaloIdVT_GsfTrkIdT);
   T1->SetBranchAddress("HLT_Ele135_CaloIdVT_GsfTrkIdT", &HLT_Ele135_CaloIdVT_GsfTrkIdT, &b_HLT_Ele135_CaloIdVT_GsfTrkIdT);
   T1->SetBranchAddress("HLT_Ele145_CaloIdVT_GsfTrkIdT", &HLT_Ele145_CaloIdVT_GsfTrkIdT, &b_HLT_Ele145_CaloIdVT_GsfTrkIdT);
   T1->SetBranchAddress("HLT_Ele200_CaloIdVT_GsfTrkIdT", &HLT_Ele200_CaloIdVT_GsfTrkIdT, &b_HLT_Ele200_CaloIdVT_GsfTrkIdT);
   T1->SetBranchAddress("HLT_Ele250_CaloIdVT_GsfTrkIdT", &HLT_Ele250_CaloIdVT_GsfTrkIdT, &b_HLT_Ele250_CaloIdVT_GsfTrkIdT);
   T1->SetBranchAddress("HLT_Ele300_CaloIdVT_GsfTrkIdT", &HLT_Ele300_CaloIdVT_GsfTrkIdT, &b_HLT_Ele300_CaloIdVT_GsfTrkIdT);
   T1->SetBranchAddress("HLT_PFHT300PT30_QuadPFJet_75_60_45_40", &HLT_PFHT300PT30_QuadPFJet_75_60_45_40, &b_HLT_PFHT300PT30_QuadPFJet_75_60_45_40);
   T1->SetBranchAddress("HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0", &HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0, &b_HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0);
   T1->SetBranchAddress("HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2", &HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2, &b_HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2);
   T1->SetBranchAddress("HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2", &HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2, &b_HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2);
   T1->SetBranchAddress("HLT_PFHT380_SixPFJet32", &HLT_PFHT380_SixPFJet32, &b_HLT_PFHT380_SixPFJet32);
   T1->SetBranchAddress("HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5", &HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5, &b_HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5);
   T1->SetBranchAddress("HLT_PFHT430_SixPFJet40", &HLT_PFHT430_SixPFJet40, &b_HLT_PFHT430_SixPFJet40);
   T1->SetBranchAddress("HLT_PFHT350", &HLT_PFHT350, &b_HLT_PFHT350);
   T1->SetBranchAddress("HLT_PFHT350MinPFJet15", &HLT_PFHT350MinPFJet15, &b_HLT_PFHT350MinPFJet15);
   T1->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL", &HLT_Photon60_R9Id90_CaloIdL_IsoL, &b_HLT_Photon60_R9Id90_CaloIdL_IsoL);
   T1->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL, &b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL);
   T1->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15, &b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15);
   T1->SetBranchAddress("HLT_FullTrack_Multiplicity85", &HLT_FullTrack_Multiplicity85, &b_HLT_FullTrack_Multiplicity85);
   T1->SetBranchAddress("HLT_FullTrack_Multiplicity100", &HLT_FullTrack_Multiplicity100, &b_HLT_FullTrack_Multiplicity100);
   T1->SetBranchAddress("HLT_FullTrack_Multiplicity130", &HLT_FullTrack_Multiplicity130, &b_HLT_FullTrack_Multiplicity130);
   T1->SetBranchAddress("HLT_FullTrack_Multiplicity155", &HLT_FullTrack_Multiplicity155, &b_HLT_FullTrack_Multiplicity155);
   T1->SetBranchAddress("HLT_ECALHT800", &HLT_ECALHT800, &b_HLT_ECALHT800);
   T1->SetBranchAddress("HLT_DiSC30_18_EIso_AND_HE_Mass70", &HLT_DiSC30_18_EIso_AND_HE_Mass70, &b_HLT_DiSC30_18_EIso_AND_HE_Mass70);
   T1->SetBranchAddress("HLT_Physics", &HLT_Physics, &b_HLT_Physics);
   T1->SetBranchAddress("HLT_Physics_part0", &HLT_Physics_part0, &b_HLT_Physics_part0);
   T1->SetBranchAddress("HLT_Physics_part1", &HLT_Physics_part1, &b_HLT_Physics_part1);
   T1->SetBranchAddress("HLT_Physics_part2", &HLT_Physics_part2, &b_HLT_Physics_part2);
   T1->SetBranchAddress("HLT_Physics_part3", &HLT_Physics_part3, &b_HLT_Physics_part3);
   T1->SetBranchAddress("HLT_Physics_part4", &HLT_Physics_part4, &b_HLT_Physics_part4);
   T1->SetBranchAddress("HLT_Physics_part5", &HLT_Physics_part5, &b_HLT_Physics_part5);
   T1->SetBranchAddress("HLT_Physics_part6", &HLT_Physics_part6, &b_HLT_Physics_part6);
   T1->SetBranchAddress("HLT_Physics_part7", &HLT_Physics_part7, &b_HLT_Physics_part7);
   T1->SetBranchAddress("HLT_Random", &HLT_Random, &b_HLT_Random);
   T1->SetBranchAddress("HLT_ZeroBias", &HLT_ZeroBias, &b_HLT_ZeroBias);
   T1->SetBranchAddress("HLT_ZeroBias_part0", &HLT_ZeroBias_part0, &b_HLT_ZeroBias_part0);
   T1->SetBranchAddress("HLT_ZeroBias_part1", &HLT_ZeroBias_part1, &b_HLT_ZeroBias_part1);
   T1->SetBranchAddress("HLT_ZeroBias_part2", &HLT_ZeroBias_part2, &b_HLT_ZeroBias_part2);
   T1->SetBranchAddress("HLT_ZeroBias_part3", &HLT_ZeroBias_part3, &b_HLT_ZeroBias_part3);
   T1->SetBranchAddress("HLT_ZeroBias_part4", &HLT_ZeroBias_part4, &b_HLT_ZeroBias_part4);
   T1->SetBranchAddress("HLT_ZeroBias_part5", &HLT_ZeroBias_part5, &b_HLT_ZeroBias_part5);
   T1->SetBranchAddress("HLT_ZeroBias_part6", &HLT_ZeroBias_part6, &b_HLT_ZeroBias_part6);
   T1->SetBranchAddress("HLT_ZeroBias_part7", &HLT_ZeroBias_part7, &b_HLT_ZeroBias_part7);
   T1->SetBranchAddress("HLT_AK4CaloJet30", &HLT_AK4CaloJet30, &b_HLT_AK4CaloJet30);
   T1->SetBranchAddress("HLT_AK4CaloJet40", &HLT_AK4CaloJet40, &b_HLT_AK4CaloJet40);
   T1->SetBranchAddress("HLT_AK4CaloJet50", &HLT_AK4CaloJet50, &b_HLT_AK4CaloJet50);
   T1->SetBranchAddress("HLT_AK4CaloJet80", &HLT_AK4CaloJet80, &b_HLT_AK4CaloJet80);
   T1->SetBranchAddress("HLT_AK4CaloJet100", &HLT_AK4CaloJet100, &b_HLT_AK4CaloJet100);
   T1->SetBranchAddress("HLT_AK4CaloJet120", &HLT_AK4CaloJet120, &b_HLT_AK4CaloJet120);
   T1->SetBranchAddress("HLT_AK4PFJet30", &HLT_AK4PFJet30, &b_HLT_AK4PFJet30);
   T1->SetBranchAddress("HLT_AK4PFJet50", &HLT_AK4PFJet50, &b_HLT_AK4PFJet50);
   T1->SetBranchAddress("HLT_AK4PFJet80", &HLT_AK4PFJet80, &b_HLT_AK4PFJet80);
   T1->SetBranchAddress("HLT_AK4PFJet100", &HLT_AK4PFJet100, &b_HLT_AK4PFJet100);
   T1->SetBranchAddress("HLT_AK4PFJet120", &HLT_AK4PFJet120, &b_HLT_AK4PFJet120);
   T1->SetBranchAddress("HLT_HISinglePhoton10_Eta3p1ForPPRef", &HLT_HISinglePhoton10_Eta3p1ForPPRef, &b_HLT_HISinglePhoton10_Eta3p1ForPPRef);
   T1->SetBranchAddress("HLT_HISinglePhoton20_Eta3p1ForPPRef", &HLT_HISinglePhoton20_Eta3p1ForPPRef, &b_HLT_HISinglePhoton20_Eta3p1ForPPRef);
   T1->SetBranchAddress("HLT_HISinglePhoton30_Eta3p1ForPPRef", &HLT_HISinglePhoton30_Eta3p1ForPPRef, &b_HLT_HISinglePhoton30_Eta3p1ForPPRef);
   T1->SetBranchAddress("HLT_HISinglePhoton40_Eta3p1ForPPRef", &HLT_HISinglePhoton40_Eta3p1ForPPRef, &b_HLT_HISinglePhoton40_Eta3p1ForPPRef);
   T1->SetBranchAddress("HLT_HISinglePhoton50_Eta3p1ForPPRef", &HLT_HISinglePhoton50_Eta3p1ForPPRef, &b_HLT_HISinglePhoton50_Eta3p1ForPPRef);
   T1->SetBranchAddress("HLT_HISinglePhoton60_Eta3p1ForPPRef", &HLT_HISinglePhoton60_Eta3p1ForPPRef, &b_HLT_HISinglePhoton60_Eta3p1ForPPRef);
   T1->SetBranchAddress("HLT_Photon20_HoverELoose", &HLT_Photon20_HoverELoose, &b_HLT_Photon20_HoverELoose);
   T1->SetBranchAddress("HLT_Photon30_HoverELoose", &HLT_Photon30_HoverELoose, &b_HLT_Photon30_HoverELoose);
   T1->SetBranchAddress("HLT_Photon40_HoverELoose", &HLT_Photon40_HoverELoose, &b_HLT_Photon40_HoverELoose);
   T1->SetBranchAddress("HLT_Photon50_HoverELoose", &HLT_Photon50_HoverELoose, &b_HLT_Photon50_HoverELoose);
   T1->SetBranchAddress("HLT_Photon60_HoverELoose", &HLT_Photon60_HoverELoose, &b_HLT_Photon60_HoverELoose);
   T1->SetBranchAddress("HLT_EcalCalibration", &HLT_EcalCalibration, &b_HLT_EcalCalibration);
   T1->SetBranchAddress("HLT_HcalCalibration", &HLT_HcalCalibration, &b_HLT_HcalCalibration);
   T1->SetBranchAddress("HLT_L1UnpairedBunchBptxMinus", &HLT_L1UnpairedBunchBptxMinus, &b_HLT_L1UnpairedBunchBptxMinus);
   T1->SetBranchAddress("HLT_L1UnpairedBunchBptxPlus", &HLT_L1UnpairedBunchBptxPlus, &b_HLT_L1UnpairedBunchBptxPlus);
   T1->SetBranchAddress("HLT_L1NotBptxOR", &HLT_L1NotBptxOR, &b_HLT_L1NotBptxOR);
   T1->SetBranchAddress("HLT_L1MinimumBiasHF_OR", &HLT_L1MinimumBiasHF_OR, &b_HLT_L1MinimumBiasHF_OR);
   T1->SetBranchAddress("HLT_L1MinimumBiasHF0OR", &HLT_L1MinimumBiasHF0OR, &b_HLT_L1MinimumBiasHF0OR);
   T1->SetBranchAddress("HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142, &b_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142);
   T1->SetBranchAddress("HLT_HcalNZS", &HLT_HcalNZS, &b_HLT_HcalNZS);
   T1->SetBranchAddress("HLT_HcalPhiSym", &HLT_HcalPhiSym, &b_HLT_HcalPhiSym);
   T1->SetBranchAddress("HLT_HcalIsolatedbunch", &HLT_HcalIsolatedbunch, &b_HLT_HcalIsolatedbunch);
   T1->SetBranchAddress("HLT_IsoTrackHB", &HLT_IsoTrackHB, &b_HLT_IsoTrackHB);
   T1->SetBranchAddress("HLT_IsoTrackHE", &HLT_IsoTrackHE, &b_HLT_IsoTrackHE);
   T1->SetBranchAddress("HLT_ZeroBias_FirstCollisionAfterAbortGap", &HLT_ZeroBias_FirstCollisionAfterAbortGap, &b_HLT_ZeroBias_FirstCollisionAfterAbortGap);
   T1->SetBranchAddress("HLT_ZeroBias_IsolatedBunches", &HLT_ZeroBias_IsolatedBunches, &b_HLT_ZeroBias_IsolatedBunches);
   T1->SetBranchAddress("HLT_ZeroBias_FirstCollisionInTrain", &HLT_ZeroBias_FirstCollisionInTrain, &b_HLT_ZeroBias_FirstCollisionInTrain);
   T1->SetBranchAddress("HLT_ZeroBias_LastCollisionInTrain", &HLT_ZeroBias_LastCollisionInTrain, &b_HLT_ZeroBias_LastCollisionInTrain);
   T1->SetBranchAddress("HLT_ZeroBias_FirstBXAfterTrain", &HLT_ZeroBias_FirstBXAfterTrain, &b_HLT_ZeroBias_FirstBXAfterTrain);
   T1->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1);
   T1->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1);
   T1->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1);
   T1->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1);
   T1->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1);
   T1->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1);
   T1->SetBranchAddress("HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg", &HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg, &b_HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg);
   T1->SetBranchAddress("HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg", &HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg, &b_HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg);
   T1->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg, &b_HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg);
   T1->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg, &b_HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg);
   T1->SetBranchAddress("HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg, &b_HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg);
   T1->SetBranchAddress("HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg, &b_HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg);
   T1->SetBranchAddress("HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg", &HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg, &b_HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg);
   T1->SetBranchAddress("HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg", &HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg, &b_HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg);
   T1->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg, &b_HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg);
   T1->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg, &b_HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg);
   T1->SetBranchAddress("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg, &b_HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg);
   T1->SetBranchAddress("HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg, &b_HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg);
   T1->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr);
   T1->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90, &b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90);
   T1->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100, &b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100);
   T1->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110, &b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110);
   T1->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120, &b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120);
   T1->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130, &b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130);
   T1->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr, &b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr);
   T1->SetBranchAddress("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr", &HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr, &b_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr);
   T1->SetBranchAddress("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1", &HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1, &b_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1);
   T1->SetBranchAddress("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1);
   T1->SetBranchAddress("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1);
   T1->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1);
   T1->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1);
   T1->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1);
   T1->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1);
   T1->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1);
   T1->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1);
   T1->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1);
   T1->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1);
   T1->SetBranchAddress("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", &HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL, &b_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL);
   T1->SetBranchAddress("HLT_Rsq0p35", &HLT_Rsq0p35, &b_HLT_Rsq0p35);
   T1->SetBranchAddress("HLT_Rsq0p40", &HLT_Rsq0p40, &b_HLT_Rsq0p40);
   T1->SetBranchAddress("HLT_RsqMR300_Rsq0p09_MR200", &HLT_RsqMR300_Rsq0p09_MR200, &b_HLT_RsqMR300_Rsq0p09_MR200);
   T1->SetBranchAddress("HLT_RsqMR320_Rsq0p09_MR200", &HLT_RsqMR320_Rsq0p09_MR200, &b_HLT_RsqMR320_Rsq0p09_MR200);
   T1->SetBranchAddress("HLT_RsqMR300_Rsq0p09_MR200_4jet", &HLT_RsqMR300_Rsq0p09_MR200_4jet, &b_HLT_RsqMR300_Rsq0p09_MR200_4jet);
   T1->SetBranchAddress("HLT_RsqMR320_Rsq0p09_MR200_4jet", &HLT_RsqMR320_Rsq0p09_MR200_4jet, &b_HLT_RsqMR320_Rsq0p09_MR200_4jet);
   T1->SetBranchAddress("HLT_L1_DoubleJet30_Mass_Min400_Mu10", &HLT_L1_DoubleJet30_Mass_Min400_Mu10, &b_HLT_L1_DoubleJet30_Mass_Min400_Mu10);
   T1->SetBranchAddress("HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1", &HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1, &b_HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1);
   T1->SetBranchAddress("HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1", &HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1, &b_HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1);
   T1->SetBranchAddress("HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1", &HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1, &b_HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1);
   T1->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50", &HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50, &b_HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50);
   T1->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3", &HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3, &b_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3);
   T1->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3", &HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3, &b_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3);
   T1->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight_PFHT60", &HLT_PFMET100_PFMHT100_IDTight_PFHT60, &b_HLT_PFMET100_PFMHT100_IDTight_PFHT60);
   T1->SetBranchAddress("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60", &HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60, &b_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60);
   T1->SetBranchAddress("HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60", &HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60, &b_HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60);
   T1->SetBranchAddress("HLT_Mu18_Mu9_SameSign", &HLT_Mu18_Mu9_SameSign, &b_HLT_Mu18_Mu9_SameSign);
   T1->SetBranchAddress("HLT_Mu18_Mu9_SameSign_DZ", &HLT_Mu18_Mu9_SameSign_DZ, &b_HLT_Mu18_Mu9_SameSign_DZ);
   T1->SetBranchAddress("HLT_Mu18_Mu9", &HLT_Mu18_Mu9, &b_HLT_Mu18_Mu9);
   T1->SetBranchAddress("HLT_Mu18_Mu9_DZ", &HLT_Mu18_Mu9_DZ, &b_HLT_Mu18_Mu9_DZ);
   T1->SetBranchAddress("HLT_Mu20_Mu10_SameSign", &HLT_Mu20_Mu10_SameSign, &b_HLT_Mu20_Mu10_SameSign);
   T1->SetBranchAddress("HLT_Mu20_Mu10_SameSign_DZ", &HLT_Mu20_Mu10_SameSign_DZ, &b_HLT_Mu20_Mu10_SameSign_DZ);
   T1->SetBranchAddress("HLT_Mu20_Mu10", &HLT_Mu20_Mu10, &b_HLT_Mu20_Mu10);
   T1->SetBranchAddress("HLT_Mu20_Mu10_DZ", &HLT_Mu20_Mu10_DZ, &b_HLT_Mu20_Mu10_DZ);
   T1->SetBranchAddress("HLT_Mu23_Mu12_SameSign", &HLT_Mu23_Mu12_SameSign, &b_HLT_Mu23_Mu12_SameSign);
   T1->SetBranchAddress("HLT_Mu23_Mu12_SameSign_DZ", &HLT_Mu23_Mu12_SameSign_DZ, &b_HLT_Mu23_Mu12_SameSign_DZ);
   T1->SetBranchAddress("HLT_Mu23_Mu12", &HLT_Mu23_Mu12, &b_HLT_Mu23_Mu12);
   T1->SetBranchAddress("HLT_Mu23_Mu12_DZ", &HLT_Mu23_Mu12_DZ, &b_HLT_Mu23_Mu12_DZ);
   T1->SetBranchAddress("HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi", &HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi, &b_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi);
   T1->SetBranchAddress("HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi", &HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi, &b_HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi);
   T1->SetBranchAddress("HLT_DoubleMu3_DCA_PFMET50_PFMHT60", &HLT_DoubleMu3_DCA_PFMET50_PFMHT60, &b_HLT_DoubleMu3_DCA_PFMET50_PFMHT60);
   T1->SetBranchAddress("HLT_TripleMu_5_3_3_Mass3p8to60_DCA", &HLT_TripleMu_5_3_3_Mass3p8to60_DCA, &b_HLT_TripleMu_5_3_3_Mass3p8to60_DCA);
   T1->SetBranchAddress("HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1", &HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1, &b_HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1);
   T1->SetBranchAddress("HLT_QuadPFJet103_88_75_15_DoubleBTagCSV_p013_p08_VBF1", &HLT_QuadPFJet103_88_75_15_DoubleBTagCSV_p013_p08_VBF1, &b_HLT_QuadPFJet103_88_75_15_DoubleBTagCSV_p013_p08_VBF1);
   T1->SetBranchAddress("HLT_QuadPFJet105_90_76_15_DoubleBTagCSV_p013_p08_VBF1", &HLT_QuadPFJet105_90_76_15_DoubleBTagCSV_p013_p08_VBF1, &b_HLT_QuadPFJet105_90_76_15_DoubleBTagCSV_p013_p08_VBF1);
   T1->SetBranchAddress("HLT_QuadPFJet111_90_80_15_DoubleBTagCSV_p013_p08_VBF1", &HLT_QuadPFJet111_90_80_15_DoubleBTagCSV_p013_p08_VBF1, &b_HLT_QuadPFJet111_90_80_15_DoubleBTagCSV_p013_p08_VBF1);
   T1->SetBranchAddress("HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2", &HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2, &b_HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2);
   T1->SetBranchAddress("HLT_QuadPFJet103_88_75_15_BTagCSV_p013_VBF2", &HLT_QuadPFJet103_88_75_15_BTagCSV_p013_VBF2, &b_HLT_QuadPFJet103_88_75_15_BTagCSV_p013_VBF2);
   T1->SetBranchAddress("HLT_QuadPFJet105_88_76_15_BTagCSV_p013_VBF2", &HLT_QuadPFJet105_88_76_15_BTagCSV_p013_VBF2, &b_HLT_QuadPFJet105_88_76_15_BTagCSV_p013_VBF2);
   T1->SetBranchAddress("HLT_QuadPFJet111_90_80_15_BTagCSV_p013_VBF2", &HLT_QuadPFJet111_90_80_15_BTagCSV_p013_VBF2, &b_HLT_QuadPFJet111_90_80_15_BTagCSV_p013_VBF2);
   T1->SetBranchAddress("HLT_QuadPFJet98_83_71_15", &HLT_QuadPFJet98_83_71_15, &b_HLT_QuadPFJet98_83_71_15);
   T1->SetBranchAddress("HLT_QuadPFJet103_88_75_15", &HLT_QuadPFJet103_88_75_15, &b_HLT_QuadPFJet103_88_75_15);
   T1->SetBranchAddress("HLT_QuadPFJet105_88_76_15", &HLT_QuadPFJet105_88_76_15, &b_HLT_QuadPFJet105_88_76_15);
   T1->SetBranchAddress("HLT_QuadPFJet111_90_80_15", &HLT_QuadPFJet111_90_80_15, &b_HLT_QuadPFJet111_90_80_15);
   T1->SetBranchAddress("HLT_AK8PFJet330_PFAK8BTagCSV_p17", &HLT_AK8PFJet330_PFAK8BTagCSV_p17, &b_HLT_AK8PFJet330_PFAK8BTagCSV_p17);
   T1->SetBranchAddress("HLT_AK8PFJet330_PFAK8BTagCSV_p1", &HLT_AK8PFJet330_PFAK8BTagCSV_p1, &b_HLT_AK8PFJet330_PFAK8BTagCSV_p1);
   T1->SetBranchAddress("HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55", &HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55, &b_HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55);
   T1->SetBranchAddress("HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55", &HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55, &b_HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55);
   T1->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath, &b_HLTriggerFinalPath);
   T1->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   T1->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   T1->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   T1->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter", &Flag_CSCTightHaloTrkMuUnvetoFilter, &b_Flag_CSCTightHaloTrkMuUnvetoFilter);
   T1->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, &b_Flag_CSCTightHalo2015Filter);
   T1->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
   T1->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   T1->SetBranchAddress("Flag_HcalStripHaloFilter", &Flag_HcalStripHaloFilter, &b_Flag_HcalStripHaloFilter);
   T1->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
   T1->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   T1->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, &b_Flag_EcalDeadCellBoundaryEnergyFilter);
   T1->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   T1->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   T1->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   T1->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
   T1->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters, &b_Flag_trkPOGFilters);
   T1->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter, &b_Flag_chargedHadronTrackResolutionFilter);
   T1->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter, &b_Flag_muonBadTrackFilter);
   T1->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   T1->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   T1->SetBranchAddress("Flag_BadChargedCandidateSummer16Filter", &Flag_BadChargedCandidateSummer16Filter, &b_Flag_BadChargedCandidateSummer16Filter);
   T1->SetBranchAddress("Flag_BadPFMuonSummer16Filter", &Flag_BadPFMuonSummer16Filter, &b_Flag_BadPFMuonSummer16Filter);
   T1->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, &b_Flag_trkPOG_manystripclus53X);
   T1->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, &b_Flag_trkPOG_toomanystripclus53X);
   T1->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, &b_Flag_trkPOG_logErrorTooManyClusters);
   T1->SetBranchAddress("Flag_METFilters", &Flag_METFilters, &b_Flag_METFilters);
   T1->SetBranchAddress("L1_AlwaysTrue", &L1_AlwaysTrue, &b_L1_AlwaysTrue);
   T1->SetBranchAddress("L1_BPTX_AND_Ref1_VME", &L1_BPTX_AND_Ref1_VME, &b_L1_BPTX_AND_Ref1_VME);
   T1->SetBranchAddress("L1_BPTX_AND_Ref3_VME", &L1_BPTX_AND_Ref3_VME, &b_L1_BPTX_AND_Ref3_VME);
   T1->SetBranchAddress("L1_BPTX_AND_Ref4_VME", &L1_BPTX_AND_Ref4_VME, &b_L1_BPTX_AND_Ref4_VME);
   T1->SetBranchAddress("L1_BPTX_BeamGas_B1_VME", &L1_BPTX_BeamGas_B1_VME, &b_L1_BPTX_BeamGas_B1_VME);
   T1->SetBranchAddress("L1_BPTX_BeamGas_B2_VME", &L1_BPTX_BeamGas_B2_VME, &b_L1_BPTX_BeamGas_B2_VME);
   T1->SetBranchAddress("L1_BPTX_BeamGas_Ref1_VME", &L1_BPTX_BeamGas_Ref1_VME, &b_L1_BPTX_BeamGas_Ref1_VME);
   T1->SetBranchAddress("L1_BPTX_BeamGas_Ref2_VME", &L1_BPTX_BeamGas_Ref2_VME, &b_L1_BPTX_BeamGas_Ref2_VME);
   T1->SetBranchAddress("L1_BPTX_NotOR_VME", &L1_BPTX_NotOR_VME, &b_L1_BPTX_NotOR_VME);
   T1->SetBranchAddress("L1_BPTX_OR_Ref3_VME", &L1_BPTX_OR_Ref3_VME, &b_L1_BPTX_OR_Ref3_VME);
   T1->SetBranchAddress("L1_BPTX_OR_Ref4_VME", &L1_BPTX_OR_Ref4_VME, &b_L1_BPTX_OR_Ref4_VME);
   T1->SetBranchAddress("L1_BPTX_RefAND_VME", &L1_BPTX_RefAND_VME, &b_L1_BPTX_RefAND_VME);
   T1->SetBranchAddress("L1_BptxMinus", &L1_BptxMinus, &b_L1_BptxMinus);
   T1->SetBranchAddress("L1_BptxOR", &L1_BptxOR, &b_L1_BptxOR);
   T1->SetBranchAddress("L1_BptxPlus", &L1_BptxPlus, &b_L1_BptxPlus);
   T1->SetBranchAddress("L1_BptxXOR", &L1_BptxXOR, &b_L1_BptxXOR);
   T1->SetBranchAddress("L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142, &b_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142);
   T1->SetBranchAddress("L1_DoubleEG6_HTT240er", &L1_DoubleEG6_HTT240er, &b_L1_DoubleEG6_HTT240er);
   T1->SetBranchAddress("L1_DoubleEG6_HTT250er", &L1_DoubleEG6_HTT250er, &b_L1_DoubleEG6_HTT250er);
   T1->SetBranchAddress("L1_DoubleEG6_HTT255er", &L1_DoubleEG6_HTT255er, &b_L1_DoubleEG6_HTT255er);
   T1->SetBranchAddress("L1_DoubleEG6_HTT270er", &L1_DoubleEG6_HTT270er, &b_L1_DoubleEG6_HTT270er);
   T1->SetBranchAddress("L1_DoubleEG6_HTT300er", &L1_DoubleEG6_HTT300er, &b_L1_DoubleEG6_HTT300er);
   T1->SetBranchAddress("L1_DoubleEG8er2p6_HTT255er", &L1_DoubleEG8er2p6_HTT255er, &b_L1_DoubleEG8er2p6_HTT255er);
   T1->SetBranchAddress("L1_DoubleEG8er2p6_HTT270er", &L1_DoubleEG8er2p6_HTT270er, &b_L1_DoubleEG8er2p6_HTT270er);
   T1->SetBranchAddress("L1_DoubleEG8er2p6_HTT300er", &L1_DoubleEG8er2p6_HTT300er, &b_L1_DoubleEG8er2p6_HTT300er);
   T1->SetBranchAddress("L1_DoubleEG_15_10", &L1_DoubleEG_15_10, &b_L1_DoubleEG_15_10);
   T1->SetBranchAddress("L1_DoubleEG_18_17", &L1_DoubleEG_18_17, &b_L1_DoubleEG_18_17);
   T1->SetBranchAddress("L1_DoubleEG_20_18", &L1_DoubleEG_20_18, &b_L1_DoubleEG_20_18);
   T1->SetBranchAddress("L1_DoubleEG_22_10", &L1_DoubleEG_22_10, &b_L1_DoubleEG_22_10);
   T1->SetBranchAddress("L1_DoubleEG_22_12", &L1_DoubleEG_22_12, &b_L1_DoubleEG_22_12);
   T1->SetBranchAddress("L1_DoubleEG_22_15", &L1_DoubleEG_22_15, &b_L1_DoubleEG_22_15);
   T1->SetBranchAddress("L1_DoubleEG_23_10", &L1_DoubleEG_23_10, &b_L1_DoubleEG_23_10);
   T1->SetBranchAddress("L1_DoubleEG_24_17", &L1_DoubleEG_24_17, &b_L1_DoubleEG_24_17);
   T1->SetBranchAddress("L1_DoubleEG_25_12", &L1_DoubleEG_25_12, &b_L1_DoubleEG_25_12);
   T1->SetBranchAddress("L1_DoubleEG_25_13", &L1_DoubleEG_25_13, &b_L1_DoubleEG_25_13);
   T1->SetBranchAddress("L1_DoubleEG_25_14", &L1_DoubleEG_25_14, &b_L1_DoubleEG_25_14);
   T1->SetBranchAddress("L1_DoubleEG_LooseIso23_10", &L1_DoubleEG_LooseIso23_10, &b_L1_DoubleEG_LooseIso23_10);
   T1->SetBranchAddress("L1_DoubleEG_LooseIso24_10", &L1_DoubleEG_LooseIso24_10, &b_L1_DoubleEG_LooseIso24_10);
   T1->SetBranchAddress("L1_DoubleIsoTau28er2p1", &L1_DoubleIsoTau28er2p1, &b_L1_DoubleIsoTau28er2p1);
   T1->SetBranchAddress("L1_DoubleIsoTau30er2p1", &L1_DoubleIsoTau30er2p1, &b_L1_DoubleIsoTau30er2p1);
   T1->SetBranchAddress("L1_DoubleIsoTau32er2p1", &L1_DoubleIsoTau32er2p1, &b_L1_DoubleIsoTau32er2p1);
   T1->SetBranchAddress("L1_DoubleIsoTau33er2p1", &L1_DoubleIsoTau33er2p1, &b_L1_DoubleIsoTau33er2p1);
   T1->SetBranchAddress("L1_DoubleIsoTau34er2p1", &L1_DoubleIsoTau34er2p1, &b_L1_DoubleIsoTau34er2p1);
   T1->SetBranchAddress("L1_DoubleIsoTau35er2p1", &L1_DoubleIsoTau35er2p1, &b_L1_DoubleIsoTau35er2p1);
   T1->SetBranchAddress("L1_DoubleIsoTau36er2p1", &L1_DoubleIsoTau36er2p1, &b_L1_DoubleIsoTau36er2p1);
   T1->SetBranchAddress("L1_DoubleIsoTau38er2p1", &L1_DoubleIsoTau38er2p1, &b_L1_DoubleIsoTau38er2p1);
   T1->SetBranchAddress("L1_DoubleJet100er2p3_dEta_Max1p6", &L1_DoubleJet100er2p3_dEta_Max1p6, &b_L1_DoubleJet100er2p3_dEta_Max1p6);
   T1->SetBranchAddress("L1_DoubleJet100er2p7", &L1_DoubleJet100er2p7, &b_L1_DoubleJet100er2p7);
   T1->SetBranchAddress("L1_DoubleJet112er2p3_dEta_Max1p6", &L1_DoubleJet112er2p3_dEta_Max1p6, &b_L1_DoubleJet112er2p3_dEta_Max1p6);
   T1->SetBranchAddress("L1_DoubleJet112er2p7", &L1_DoubleJet112er2p7, &b_L1_DoubleJet112er2p7);
   T1->SetBranchAddress("L1_DoubleJet120er2p7", &L1_DoubleJet120er2p7, &b_L1_DoubleJet120er2p7);
   T1->SetBranchAddress("L1_DoubleJet150er2p7", &L1_DoubleJet150er2p7, &b_L1_DoubleJet150er2p7);
   T1->SetBranchAddress("L1_DoubleJet30_Mass_Min300_dEta_Max1p5", &L1_DoubleJet30_Mass_Min300_dEta_Max1p5, &b_L1_DoubleJet30_Mass_Min300_dEta_Max1p5);
   T1->SetBranchAddress("L1_DoubleJet30_Mass_Min320_dEta_Max1p5", &L1_DoubleJet30_Mass_Min320_dEta_Max1p5, &b_L1_DoubleJet30_Mass_Min320_dEta_Max1p5);
   T1->SetBranchAddress("L1_DoubleJet30_Mass_Min340_dEta_Max1p5", &L1_DoubleJet30_Mass_Min340_dEta_Max1p5, &b_L1_DoubleJet30_Mass_Min340_dEta_Max1p5);
   T1->SetBranchAddress("L1_DoubleJet30_Mass_Min360_dEta_Max1p5", &L1_DoubleJet30_Mass_Min360_dEta_Max1p5, &b_L1_DoubleJet30_Mass_Min360_dEta_Max1p5);
   T1->SetBranchAddress("L1_DoubleJet30_Mass_Min380_dEta_Max1p5", &L1_DoubleJet30_Mass_Min380_dEta_Max1p5, &b_L1_DoubleJet30_Mass_Min380_dEta_Max1p5);
   T1->SetBranchAddress("L1_DoubleJet30_Mass_Min400_Mu10", &L1_DoubleJet30_Mass_Min400_Mu10, &b_L1_DoubleJet30_Mass_Min400_Mu10);
   T1->SetBranchAddress("L1_DoubleJet30_Mass_Min400_Mu6", &L1_DoubleJet30_Mass_Min400_Mu6, &b_L1_DoubleJet30_Mass_Min400_Mu6);
   T1->SetBranchAddress("L1_DoubleJet30_Mass_Min400_dEta_Max1p5", &L1_DoubleJet30_Mass_Min400_dEta_Max1p5, &b_L1_DoubleJet30_Mass_Min400_dEta_Max1p5);
   T1->SetBranchAddress("L1_DoubleJet35_rmovlp_IsoTau45_Mass_Min450", &L1_DoubleJet35_rmovlp_IsoTau45_Mass_Min450, &b_L1_DoubleJet35_rmovlp_IsoTau45_Mass_Min450);
   T1->SetBranchAddress("L1_DoubleJet40er2p7", &L1_DoubleJet40er2p7, &b_L1_DoubleJet40er2p7);
   T1->SetBranchAddress("L1_DoubleJet50er2p7", &L1_DoubleJet50er2p7, &b_L1_DoubleJet50er2p7);
   T1->SetBranchAddress("L1_DoubleJet60er2p7", &L1_DoubleJet60er2p7, &b_L1_DoubleJet60er2p7);
   T1->SetBranchAddress("L1_DoubleJet60er2p7_ETM100", &L1_DoubleJet60er2p7_ETM100, &b_L1_DoubleJet60er2p7_ETM100);
   T1->SetBranchAddress("L1_DoubleJet60er2p7_ETM60", &L1_DoubleJet60er2p7_ETM60, &b_L1_DoubleJet60er2p7_ETM60);
   T1->SetBranchAddress("L1_DoubleJet60er2p7_ETM70", &L1_DoubleJet60er2p7_ETM70, &b_L1_DoubleJet60er2p7_ETM70);
   T1->SetBranchAddress("L1_DoubleJet60er2p7_ETM80", &L1_DoubleJet60er2p7_ETM80, &b_L1_DoubleJet60er2p7_ETM80);
   T1->SetBranchAddress("L1_DoubleJet60er2p7_ETM90", &L1_DoubleJet60er2p7_ETM90, &b_L1_DoubleJet60er2p7_ETM90);
   T1->SetBranchAddress("L1_DoubleJet80er2p7", &L1_DoubleJet80er2p7, &b_L1_DoubleJet80er2p7);
   T1->SetBranchAddress("L1_DoubleJet_100_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_100_30_DoubleJet30_Mass_Min620, &b_L1_DoubleJet_100_30_DoubleJet30_Mass_Min620);
   T1->SetBranchAddress("L1_DoubleJet_100_35_DoubleJet35_Mass_Min620", &L1_DoubleJet_100_35_DoubleJet35_Mass_Min620, &b_L1_DoubleJet_100_35_DoubleJet35_Mass_Min620);
   T1->SetBranchAddress("L1_DoubleJet_110_35_DoubleJet35_Mass_Min620", &L1_DoubleJet_110_35_DoubleJet35_Mass_Min620, &b_L1_DoubleJet_110_35_DoubleJet35_Mass_Min620);
   T1->SetBranchAddress("L1_DoubleJet_110_40_DoubleJet40_Mass_Min620", &L1_DoubleJet_110_40_DoubleJet40_Mass_Min620, &b_L1_DoubleJet_110_40_DoubleJet40_Mass_Min620);
   T1->SetBranchAddress("L1_DoubleJet_115_35_DoubleJet35_Mass_Min620", &L1_DoubleJet_115_35_DoubleJet35_Mass_Min620, &b_L1_DoubleJet_115_35_DoubleJet35_Mass_Min620);
   T1->SetBranchAddress("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620", &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620, &b_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620);
   T1->SetBranchAddress("L1_DoubleJet_90_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_90_30_DoubleJet30_Mass_Min620, &b_L1_DoubleJet_90_30_DoubleJet30_Mass_Min620);
   T1->SetBranchAddress("L1_DoubleLooseIsoEG22er2p1", &L1_DoubleLooseIsoEG22er2p1, &b_L1_DoubleLooseIsoEG22er2p1);
   T1->SetBranchAddress("L1_DoubleLooseIsoEG24er2p1", &L1_DoubleLooseIsoEG24er2p1, &b_L1_DoubleLooseIsoEG24er2p1);
   T1->SetBranchAddress("L1_DoubleMu0", &L1_DoubleMu0, &b_L1_DoubleMu0);
   T1->SetBranchAddress("L1_DoubleMu0_ETM40", &L1_DoubleMu0_ETM40, &b_L1_DoubleMu0_ETM40);
   T1->SetBranchAddress("L1_DoubleMu0_ETM55", &L1_DoubleMu0_ETM55, &b_L1_DoubleMu0_ETM55);
   T1->SetBranchAddress("L1_DoubleMu0_ETM60", &L1_DoubleMu0_ETM60, &b_L1_DoubleMu0_ETM60);
   T1->SetBranchAddress("L1_DoubleMu0_ETM65", &L1_DoubleMu0_ETM65, &b_L1_DoubleMu0_ETM65);
   T1->SetBranchAddress("L1_DoubleMu0_ETM70", &L1_DoubleMu0_ETM70, &b_L1_DoubleMu0_ETM70);
   T1->SetBranchAddress("L1_DoubleMu0_SQ", &L1_DoubleMu0_SQ, &b_L1_DoubleMu0_SQ);
   T1->SetBranchAddress("L1_DoubleMu0_SQ_OS", &L1_DoubleMu0_SQ_OS, &b_L1_DoubleMu0_SQ_OS);
   T1->SetBranchAddress("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4, &b_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4);
   T1->SetBranchAddress("L1_DoubleMu0er1p4_dEta_Max1p8_OS", &L1_DoubleMu0er1p4_dEta_Max1p8_OS, &b_L1_DoubleMu0er1p4_dEta_Max1p8_OS);
   T1->SetBranchAddress("L1_DoubleMu0er1p5_SQ_OS", &L1_DoubleMu0er1p5_SQ_OS, &b_L1_DoubleMu0er1p5_SQ_OS);
   T1->SetBranchAddress("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4, &b_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4);
   T1->SetBranchAddress("L1_DoubleMu0er1p5_SQ_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_dR_Max1p4, &b_L1_DoubleMu0er1p5_SQ_dR_Max1p4);
   T1->SetBranchAddress("L1_DoubleMu0er2_SQ_dR_Max1p4", &L1_DoubleMu0er2_SQ_dR_Max1p4, &b_L1_DoubleMu0er2_SQ_dR_Max1p4);
   T1->SetBranchAddress("L1_DoubleMu18er2p1", &L1_DoubleMu18er2p1, &b_L1_DoubleMu18er2p1);
   T1->SetBranchAddress("L1_DoubleMu22er2p1", &L1_DoubleMu22er2p1, &b_L1_DoubleMu22er2p1);
   T1->SetBranchAddress("L1_DoubleMu3_OS_DoubleEG7p5Upsilon", &L1_DoubleMu3_OS_DoubleEG7p5Upsilon, &b_L1_DoubleMu3_OS_DoubleEG7p5Upsilon);
   T1->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF40_Jet60_OR_DoubleJet30", &L1_DoubleMu3_SQ_ETMHF40_Jet60_OR_DoubleJet30, &b_L1_DoubleMu3_SQ_ETMHF40_Jet60_OR_DoubleJet30);
   T1->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF50_Jet60_OR_DoubleJet30", &L1_DoubleMu3_SQ_ETMHF50_Jet60_OR_DoubleJet30, &b_L1_DoubleMu3_SQ_ETMHF50_Jet60_OR_DoubleJet30);
   T1->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF60_Jet60_OR_DoubleJet30", &L1_DoubleMu3_SQ_ETMHF60_Jet60_OR_DoubleJet30, &b_L1_DoubleMu3_SQ_ETMHF60_Jet60_OR_DoubleJet30);
   T1->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF70_Jet60_OR_DoubleJet30", &L1_DoubleMu3_SQ_ETMHF70_Jet60_OR_DoubleJet30, &b_L1_DoubleMu3_SQ_ETMHF70_Jet60_OR_DoubleJet30);
   T1->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF80_Jet60_OR_DoubleJet30", &L1_DoubleMu3_SQ_ETMHF80_Jet60_OR_DoubleJet30, &b_L1_DoubleMu3_SQ_ETMHF80_Jet60_OR_DoubleJet30);
   T1->SetBranchAddress("L1_DoubleMu3_SQ_HTT100er", &L1_DoubleMu3_SQ_HTT100er, &b_L1_DoubleMu3_SQ_HTT100er);
   T1->SetBranchAddress("L1_DoubleMu3_SQ_HTT200er", &L1_DoubleMu3_SQ_HTT200er, &b_L1_DoubleMu3_SQ_HTT200er);
   T1->SetBranchAddress("L1_DoubleMu3_SQ_HTT220er", &L1_DoubleMu3_SQ_HTT220er, &b_L1_DoubleMu3_SQ_HTT220er);
   T1->SetBranchAddress("L1_DoubleMu3_SQ_HTT240er", &L1_DoubleMu3_SQ_HTT240er, &b_L1_DoubleMu3_SQ_HTT240er);
   T1->SetBranchAddress("L1_DoubleMu4_OS_EG12", &L1_DoubleMu4_OS_EG12, &b_L1_DoubleMu4_OS_EG12);
   T1->SetBranchAddress("L1_DoubleMu4_SQ_OS", &L1_DoubleMu4_SQ_OS, &b_L1_DoubleMu4_SQ_OS);
   T1->SetBranchAddress("L1_DoubleMu4_SQ_OS_dR_Max1p2", &L1_DoubleMu4_SQ_OS_dR_Max1p2, &b_L1_DoubleMu4_SQ_OS_dR_Max1p2);
   T1->SetBranchAddress("L1_DoubleMu4p5_SQ", &L1_DoubleMu4p5_SQ, &b_L1_DoubleMu4p5_SQ);
   T1->SetBranchAddress("L1_DoubleMu4p5_SQ_OS", &L1_DoubleMu4p5_SQ_OS, &b_L1_DoubleMu4p5_SQ_OS);
   T1->SetBranchAddress("L1_DoubleMu4p5_SQ_OS_dR_Max1p2", &L1_DoubleMu4p5_SQ_OS_dR_Max1p2, &b_L1_DoubleMu4p5_SQ_OS_dR_Max1p2);
   T1->SetBranchAddress("L1_DoubleMu4p5er2p0_SQ_OS", &L1_DoubleMu4p5er2p0_SQ_OS, &b_L1_DoubleMu4p5er2p0_SQ_OS);
   T1->SetBranchAddress("L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18", &L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18, &b_L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18);
   T1->SetBranchAddress("L1_DoubleMu5Upsilon_OS_DoubleEG3", &L1_DoubleMu5Upsilon_OS_DoubleEG3, &b_L1_DoubleMu5Upsilon_OS_DoubleEG3);
   T1->SetBranchAddress("L1_DoubleMu5_OS_EG12", &L1_DoubleMu5_OS_EG12, &b_L1_DoubleMu5_OS_EG12);
   T1->SetBranchAddress("L1_DoubleMu5_SQ_OS", &L1_DoubleMu5_SQ_OS, &b_L1_DoubleMu5_SQ_OS);
   T1->SetBranchAddress("L1_DoubleMu5_SQ_OS_Mass7to18", &L1_DoubleMu5_SQ_OS_Mass7to18, &b_L1_DoubleMu5_SQ_OS_Mass7to18);
   T1->SetBranchAddress("L1_DoubleMu6_SQ_OS", &L1_DoubleMu6_SQ_OS, &b_L1_DoubleMu6_SQ_OS);
   T1->SetBranchAddress("L1_DoubleMu7_EG7", &L1_DoubleMu7_EG7, &b_L1_DoubleMu7_EG7);
   T1->SetBranchAddress("L1_DoubleMu7_SQ_EG7", &L1_DoubleMu7_SQ_EG7, &b_L1_DoubleMu7_SQ_EG7);
   T1->SetBranchAddress("L1_DoubleMu8_SQ", &L1_DoubleMu8_SQ, &b_L1_DoubleMu8_SQ);
   T1->SetBranchAddress("L1_DoubleMu_10_0_dEta_Max1p8", &L1_DoubleMu_10_0_dEta_Max1p8, &b_L1_DoubleMu_10_0_dEta_Max1p8);
   T1->SetBranchAddress("L1_DoubleMu_11_4", &L1_DoubleMu_11_4, &b_L1_DoubleMu_11_4);
   T1->SetBranchAddress("L1_DoubleMu_12_5", &L1_DoubleMu_12_5, &b_L1_DoubleMu_12_5);
   T1->SetBranchAddress("L1_DoubleMu_12_8", &L1_DoubleMu_12_8, &b_L1_DoubleMu_12_8);
   T1->SetBranchAddress("L1_DoubleMu_13_6", &L1_DoubleMu_13_6, &b_L1_DoubleMu_13_6);
   T1->SetBranchAddress("L1_DoubleMu_15_5", &L1_DoubleMu_15_5, &b_L1_DoubleMu_15_5);
   T1->SetBranchAddress("L1_DoubleMu_15_5_SQ", &L1_DoubleMu_15_5_SQ, &b_L1_DoubleMu_15_5_SQ);
   T1->SetBranchAddress("L1_DoubleMu_15_7", &L1_DoubleMu_15_7, &b_L1_DoubleMu_15_7);
   T1->SetBranchAddress("L1_DoubleMu_15_7_SQ", &L1_DoubleMu_15_7_SQ, &b_L1_DoubleMu_15_7_SQ);
   T1->SetBranchAddress("L1_DoubleMu_15_7_SQ_Mass_Min4", &L1_DoubleMu_15_7_SQ_Mass_Min4, &b_L1_DoubleMu_15_7_SQ_Mass_Min4);
   T1->SetBranchAddress("L1_DoubleMu_20_2_SQ_Mass_Max20", &L1_DoubleMu_20_2_SQ_Mass_Max20, &b_L1_DoubleMu_20_2_SQ_Mass_Max20);
   T1->SetBranchAddress("L1_DoubleTau50er2p1", &L1_DoubleTau50er2p1, &b_L1_DoubleTau50er2p1);
   T1->SetBranchAddress("L1_DoubleTau70er2p1", &L1_DoubleTau70er2p1, &b_L1_DoubleTau70er2p1);
   T1->SetBranchAddress("L1_EG25er2p1_HTT125er", &L1_EG25er2p1_HTT125er, &b_L1_EG25er2p1_HTT125er);
   T1->SetBranchAddress("L1_EG27er2p1_HTT200er", &L1_EG27er2p1_HTT200er, &b_L1_EG27er2p1_HTT200er);
   T1->SetBranchAddress("L1_ETM100", &L1_ETM100, &b_L1_ETM100);
   T1->SetBranchAddress("L1_ETM100_Jet60_dPhi_Min0p4", &L1_ETM100_Jet60_dPhi_Min0p4, &b_L1_ETM100_Jet60_dPhi_Min0p4);
   T1->SetBranchAddress("L1_ETM105", &L1_ETM105, &b_L1_ETM105);
   T1->SetBranchAddress("L1_ETM110", &L1_ETM110, &b_L1_ETM110);
   T1->SetBranchAddress("L1_ETM110_Jet60_dPhi_Min0p4", &L1_ETM110_Jet60_dPhi_Min0p4, &b_L1_ETM110_Jet60_dPhi_Min0p4);
   T1->SetBranchAddress("L1_ETM115", &L1_ETM115, &b_L1_ETM115);
   T1->SetBranchAddress("L1_ETM120", &L1_ETM120, &b_L1_ETM120);
   T1->SetBranchAddress("L1_ETM150", &L1_ETM150, &b_L1_ETM150);
   T1->SetBranchAddress("L1_ETM30", &L1_ETM30, &b_L1_ETM30);
   T1->SetBranchAddress("L1_ETM40", &L1_ETM40, &b_L1_ETM40);
   T1->SetBranchAddress("L1_ETM50", &L1_ETM50, &b_L1_ETM50);
   T1->SetBranchAddress("L1_ETM60", &L1_ETM60, &b_L1_ETM60);
   T1->SetBranchAddress("L1_ETM70", &L1_ETM70, &b_L1_ETM70);
   T1->SetBranchAddress("L1_ETM75", &L1_ETM75, &b_L1_ETM75);
   T1->SetBranchAddress("L1_ETM75_Jet60_dPhi_Min0p4", &L1_ETM75_Jet60_dPhi_Min0p4, &b_L1_ETM75_Jet60_dPhi_Min0p4);
   T1->SetBranchAddress("L1_ETM80", &L1_ETM80, &b_L1_ETM80);
   T1->SetBranchAddress("L1_ETM80_Jet60_dPhi_Min0p4", &L1_ETM80_Jet60_dPhi_Min0p4, &b_L1_ETM80_Jet60_dPhi_Min0p4);
   T1->SetBranchAddress("L1_ETM85", &L1_ETM85, &b_L1_ETM85);
   T1->SetBranchAddress("L1_ETM90", &L1_ETM90, &b_L1_ETM90);
   T1->SetBranchAddress("L1_ETM90_Jet60_dPhi_Min0p4", &L1_ETM90_Jet60_dPhi_Min0p4, &b_L1_ETM90_Jet60_dPhi_Min0p4);
   T1->SetBranchAddress("L1_ETM95", &L1_ETM95, &b_L1_ETM95);
   T1->SetBranchAddress("L1_ETMHF100", &L1_ETMHF100, &b_L1_ETMHF100);
   T1->SetBranchAddress("L1_ETMHF100_HTT60er", &L1_ETMHF100_HTT60er, &b_L1_ETMHF100_HTT60er);
   T1->SetBranchAddress("L1_ETMHF100_Jet60_OR_DiJet30woTT28", &L1_ETMHF100_Jet60_OR_DiJet30woTT28, &b_L1_ETMHF100_Jet60_OR_DiJet30woTT28);
   T1->SetBranchAddress("L1_ETMHF100_Jet60_OR_DoubleJet30", &L1_ETMHF100_Jet60_OR_DoubleJet30, &b_L1_ETMHF100_Jet60_OR_DoubleJet30);
   T1->SetBranchAddress("L1_ETMHF100_Jet90_OR_DoubleJet45_OR_TripleJet30", &L1_ETMHF100_Jet90_OR_DoubleJet45_OR_TripleJet30, &b_L1_ETMHF100_Jet90_OR_DoubleJet45_OR_TripleJet30);
   T1->SetBranchAddress("L1_ETMHF110", &L1_ETMHF110, &b_L1_ETMHF110);
   T1->SetBranchAddress("L1_ETMHF110_HTT60er", &L1_ETMHF110_HTT60er, &b_L1_ETMHF110_HTT60er);
   T1->SetBranchAddress("L1_ETMHF110_Jet60_OR_DiJet30woTT28", &L1_ETMHF110_Jet60_OR_DiJet30woTT28, &b_L1_ETMHF110_Jet60_OR_DiJet30woTT28);
   T1->SetBranchAddress("L1_ETMHF110_Jet90_OR_DoubleJet45_OR_TripleJet30", &L1_ETMHF110_Jet90_OR_DoubleJet45_OR_TripleJet30, &b_L1_ETMHF110_Jet90_OR_DoubleJet45_OR_TripleJet30);
   T1->SetBranchAddress("L1_ETMHF120", &L1_ETMHF120, &b_L1_ETMHF120);
   T1->SetBranchAddress("L1_ETMHF120_HTT60er", &L1_ETMHF120_HTT60er, &b_L1_ETMHF120_HTT60er);
   T1->SetBranchAddress("L1_ETMHF120_Jet60_OR_DiJet30woTT28", &L1_ETMHF120_Jet60_OR_DiJet30woTT28, &b_L1_ETMHF120_Jet60_OR_DiJet30woTT28);
   T1->SetBranchAddress("L1_ETMHF150", &L1_ETMHF150, &b_L1_ETMHF150);
   T1->SetBranchAddress("L1_ETMHF70", &L1_ETMHF70, &b_L1_ETMHF70);
   T1->SetBranchAddress("L1_ETMHF70_Jet90_OR_DoubleJet45_OR_TripleJet30", &L1_ETMHF70_Jet90_OR_DoubleJet45_OR_TripleJet30, &b_L1_ETMHF70_Jet90_OR_DoubleJet45_OR_TripleJet30);
   T1->SetBranchAddress("L1_ETMHF80", &L1_ETMHF80, &b_L1_ETMHF80);
   T1->SetBranchAddress("L1_ETMHF80_HTT60er", &L1_ETMHF80_HTT60er, &b_L1_ETMHF80_HTT60er);
   T1->SetBranchAddress("L1_ETMHF80_Jet90_OR_DoubleJet45_OR_TripleJet30", &L1_ETMHF80_Jet90_OR_DoubleJet45_OR_TripleJet30, &b_L1_ETMHF80_Jet90_OR_DoubleJet45_OR_TripleJet30);
   T1->SetBranchAddress("L1_ETMHF90", &L1_ETMHF90, &b_L1_ETMHF90);
   T1->SetBranchAddress("L1_ETMHF90_HTT60er", &L1_ETMHF90_HTT60er, &b_L1_ETMHF90_HTT60er);
   T1->SetBranchAddress("L1_ETMHF90_Jet90_OR_DoubleJet45_OR_TripleJet30", &L1_ETMHF90_Jet90_OR_DoubleJet45_OR_TripleJet30, &b_L1_ETMHF90_Jet90_OR_DoubleJet45_OR_TripleJet30);
   T1->SetBranchAddress("L1_ETT100_BptxAND", &L1_ETT100_BptxAND, &b_L1_ETT100_BptxAND);
   T1->SetBranchAddress("L1_ETT110_BptxAND", &L1_ETT110_BptxAND, &b_L1_ETT110_BptxAND);
   T1->SetBranchAddress("L1_ETT40_BptxAND", &L1_ETT40_BptxAND, &b_L1_ETT40_BptxAND);
   T1->SetBranchAddress("L1_ETT50_BptxAND", &L1_ETT50_BptxAND, &b_L1_ETT50_BptxAND);
   T1->SetBranchAddress("L1_ETT60_BptxAND", &L1_ETT60_BptxAND, &b_L1_ETT60_BptxAND);
   T1->SetBranchAddress("L1_ETT70_BptxAND", &L1_ETT70_BptxAND, &b_L1_ETT70_BptxAND);
   T1->SetBranchAddress("L1_ETT75_BptxAND", &L1_ETT75_BptxAND, &b_L1_ETT75_BptxAND);
   T1->SetBranchAddress("L1_ETT80_BptxAND", &L1_ETT80_BptxAND, &b_L1_ETT80_BptxAND);
   T1->SetBranchAddress("L1_ETT85_BptxAND", &L1_ETT85_BptxAND, &b_L1_ETT85_BptxAND);
   T1->SetBranchAddress("L1_ETT90_BptxAND", &L1_ETT90_BptxAND, &b_L1_ETT90_BptxAND);
   T1->SetBranchAddress("L1_ETT95_BptxAND", &L1_ETT95_BptxAND, &b_L1_ETT95_BptxAND);
   T1->SetBranchAddress("L1_FirstBunchAfterTrain", &L1_FirstBunchAfterTrain, &b_L1_FirstBunchAfterTrain);
   T1->SetBranchAddress("L1_FirstBunchInTrain", &L1_FirstBunchInTrain, &b_L1_FirstBunchInTrain);
   T1->SetBranchAddress("L1_FirstCollisionInOrbit", &L1_FirstCollisionInOrbit, &b_L1_FirstCollisionInOrbit);
   T1->SetBranchAddress("L1_FirstCollisionInTrain", &L1_FirstCollisionInTrain, &b_L1_FirstCollisionInTrain);
   T1->SetBranchAddress("L1_HTT120er", &L1_HTT120er, &b_L1_HTT120er);
   T1->SetBranchAddress("L1_HTT160er", &L1_HTT160er, &b_L1_HTT160er);
   T1->SetBranchAddress("L1_HTT200er", &L1_HTT200er, &b_L1_HTT200er);
   T1->SetBranchAddress("L1_HTT220er", &L1_HTT220er, &b_L1_HTT220er);
   T1->SetBranchAddress("L1_HTT240er", &L1_HTT240er, &b_L1_HTT240er);
   T1->SetBranchAddress("L1_HTT250er_QuadJet_70_55_40_35_er2p5", &L1_HTT250er_QuadJet_70_55_40_35_er2p5, &b_L1_HTT250er_QuadJet_70_55_40_35_er2p5);
   T1->SetBranchAddress("L1_HTT255er", &L1_HTT255er, &b_L1_HTT255er);
   T1->SetBranchAddress("L1_HTT270er", &L1_HTT270er, &b_L1_HTT270er);
   T1->SetBranchAddress("L1_HTT280er", &L1_HTT280er, &b_L1_HTT280er);
   T1->SetBranchAddress("L1_HTT280er_QuadJet_70_55_40_35_er2p5", &L1_HTT280er_QuadJet_70_55_40_35_er2p5, &b_L1_HTT280er_QuadJet_70_55_40_35_er2p5);
   T1->SetBranchAddress("L1_HTT300er", &L1_HTT300er, &b_L1_HTT300er);
   T1->SetBranchAddress("L1_HTT300er_QuadJet_70_55_40_35_er2p5", &L1_HTT300er_QuadJet_70_55_40_35_er2p5, &b_L1_HTT300er_QuadJet_70_55_40_35_er2p5);
   T1->SetBranchAddress("L1_HTT320er", &L1_HTT320er, &b_L1_HTT320er);
   T1->SetBranchAddress("L1_HTT320er_QuadJet_70_55_40_40_er2p4", &L1_HTT320er_QuadJet_70_55_40_40_er2p4, &b_L1_HTT320er_QuadJet_70_55_40_40_er2p4);
   T1->SetBranchAddress("L1_HTT320er_QuadJet_70_55_40_40_er2p5", &L1_HTT320er_QuadJet_70_55_40_40_er2p5, &b_L1_HTT320er_QuadJet_70_55_40_40_er2p5);
   T1->SetBranchAddress("L1_HTT320er_QuadJet_70_55_45_45_er2p5", &L1_HTT320er_QuadJet_70_55_45_45_er2p5, &b_L1_HTT320er_QuadJet_70_55_45_45_er2p5);
   T1->SetBranchAddress("L1_HTT340er", &L1_HTT340er, &b_L1_HTT340er);
   T1->SetBranchAddress("L1_HTT340er_QuadJet_70_55_40_40_er2p5", &L1_HTT340er_QuadJet_70_55_40_40_er2p5, &b_L1_HTT340er_QuadJet_70_55_40_40_er2p5);
   T1->SetBranchAddress("L1_HTT340er_QuadJet_70_55_45_45_er2p5", &L1_HTT340er_QuadJet_70_55_45_45_er2p5, &b_L1_HTT340er_QuadJet_70_55_45_45_er2p5);
   T1->SetBranchAddress("L1_HTT380er", &L1_HTT380er, &b_L1_HTT380er);
   T1->SetBranchAddress("L1_HTT400er", &L1_HTT400er, &b_L1_HTT400er);
   T1->SetBranchAddress("L1_HTT450er", &L1_HTT450er, &b_L1_HTT450er);
   T1->SetBranchAddress("L1_HTT500er", &L1_HTT500er, &b_L1_HTT500er);
   T1->SetBranchAddress("L1_IsoEG33_Mt40", &L1_IsoEG33_Mt40, &b_L1_IsoEG33_Mt40);
   T1->SetBranchAddress("L1_IsoEG33_Mt44", &L1_IsoEG33_Mt44, &b_L1_IsoEG33_Mt44);
   T1->SetBranchAddress("L1_IsoEG33_Mt48", &L1_IsoEG33_Mt48, &b_L1_IsoEG33_Mt48);
   T1->SetBranchAddress("L1_IsoTau40er_ETM100", &L1_IsoTau40er_ETM100, &b_L1_IsoTau40er_ETM100);
   T1->SetBranchAddress("L1_IsoTau40er_ETM105", &L1_IsoTau40er_ETM105, &b_L1_IsoTau40er_ETM105);
   T1->SetBranchAddress("L1_IsoTau40er_ETM110", &L1_IsoTau40er_ETM110, &b_L1_IsoTau40er_ETM110);
   T1->SetBranchAddress("L1_IsoTau40er_ETM115", &L1_IsoTau40er_ETM115, &b_L1_IsoTau40er_ETM115);
   T1->SetBranchAddress("L1_IsoTau40er_ETM120", &L1_IsoTau40er_ETM120, &b_L1_IsoTau40er_ETM120);
   T1->SetBranchAddress("L1_IsoTau40er_ETM80", &L1_IsoTau40er_ETM80, &b_L1_IsoTau40er_ETM80);
   T1->SetBranchAddress("L1_IsoTau40er_ETM85", &L1_IsoTau40er_ETM85, &b_L1_IsoTau40er_ETM85);
   T1->SetBranchAddress("L1_IsoTau40er_ETM90", &L1_IsoTau40er_ETM90, &b_L1_IsoTau40er_ETM90);
   T1->SetBranchAddress("L1_IsoTau40er_ETM95", &L1_IsoTau40er_ETM95, &b_L1_IsoTau40er_ETM95);
   T1->SetBranchAddress("L1_IsoTau40er_ETMHF100", &L1_IsoTau40er_ETMHF100, &b_L1_IsoTau40er_ETMHF100);
   T1->SetBranchAddress("L1_IsoTau40er_ETMHF110", &L1_IsoTau40er_ETMHF110, &b_L1_IsoTau40er_ETMHF110);
   T1->SetBranchAddress("L1_IsoTau40er_ETMHF120", &L1_IsoTau40er_ETMHF120, &b_L1_IsoTau40er_ETMHF120);
   T1->SetBranchAddress("L1_IsoTau40er_ETMHF80", &L1_IsoTau40er_ETMHF80, &b_L1_IsoTau40er_ETMHF80);
   T1->SetBranchAddress("L1_IsoTau40er_ETMHF90", &L1_IsoTau40er_ETMHF90, &b_L1_IsoTau40er_ETMHF90);
   T1->SetBranchAddress("L1_IsolatedBunch", &L1_IsolatedBunch, &b_L1_IsolatedBunch);
   T1->SetBranchAddress("L1_LastCollisionInTrain", &L1_LastCollisionInTrain, &b_L1_LastCollisionInTrain);
   T1->SetBranchAddress("L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3, &b_L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3);
   T1->SetBranchAddress("L1_LooseIsoEG24er2p1_HTT100er", &L1_LooseIsoEG24er2p1_HTT100er, &b_L1_LooseIsoEG24er2p1_HTT100er);
   T1->SetBranchAddress("L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3", &L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3, &b_L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3);
   T1->SetBranchAddress("L1_LooseIsoEG24er2p1_Jet26er2p7_dR_Min0p3", &L1_LooseIsoEG24er2p1_Jet26er2p7_dR_Min0p3, &b_L1_LooseIsoEG24er2p1_Jet26er2p7_dR_Min0p3);
   T1->SetBranchAddress("L1_LooseIsoEG24er2p1_TripleJet_26er2p7_26_26er2p7", &L1_LooseIsoEG24er2p1_TripleJet_26er2p7_26_26er2p7, &b_L1_LooseIsoEG24er2p1_TripleJet_26er2p7_26_26er2p7);
   T1->SetBranchAddress("L1_LooseIsoEG26er2p1_HTT100er", &L1_LooseIsoEG26er2p1_HTT100er, &b_L1_LooseIsoEG26er2p1_HTT100er);
   T1->SetBranchAddress("L1_LooseIsoEG26er2p1_Jet34er2p7_dR_Min0p3", &L1_LooseIsoEG26er2p1_Jet34er2p7_dR_Min0p3, &b_L1_LooseIsoEG26er2p1_Jet34er2p7_dR_Min0p3);
   T1->SetBranchAddress("L1_LooseIsoEG28er2p1_HTT100er", &L1_LooseIsoEG28er2p1_HTT100er, &b_L1_LooseIsoEG28er2p1_HTT100er);
   T1->SetBranchAddress("L1_LooseIsoEG28er2p1_Jet34er2p7_dR_Min0p3", &L1_LooseIsoEG28er2p1_Jet34er2p7_dR_Min0p3, &b_L1_LooseIsoEG28er2p1_Jet34er2p7_dR_Min0p3);
   T1->SetBranchAddress("L1_LooseIsoEG30er2p1_Jet34er2p7_dR_Min0p3", &L1_LooseIsoEG30er2p1_Jet34er2p7_dR_Min0p3, &b_L1_LooseIsoEG30er2p1_Jet34er2p7_dR_Min0p3);
   T1->SetBranchAddress("L1_MU20_EG15", &L1_MU20_EG15, &b_L1_MU20_EG15);
   T1->SetBranchAddress("L1_MinimumBiasHF0_AND_BptxAND", &L1_MinimumBiasHF0_AND_BptxAND, &b_L1_MinimumBiasHF0_AND_BptxAND);
   T1->SetBranchAddress("L1_MinimumBiasHF0_OR_BptxAND", &L1_MinimumBiasHF0_OR_BptxAND, &b_L1_MinimumBiasHF0_OR_BptxAND);
   T1->SetBranchAddress("L1_Mu10er2p1_ETM30", &L1_Mu10er2p1_ETM30, &b_L1_Mu10er2p1_ETM30);
   T1->SetBranchAddress("L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6", &L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6, &b_L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6);
   T1->SetBranchAddress("L1_Mu12_EG10", &L1_Mu12_EG10, &b_L1_Mu12_EG10);
   T1->SetBranchAddress("L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6, &b_L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6);
   T1->SetBranchAddress("L1_Mu14er2p1_ETM30", &L1_Mu14er2p1_ETM30, &b_L1_Mu14er2p1_ETM30);
   T1->SetBranchAddress("L1_Mu15_HTT100er", &L1_Mu15_HTT100er, &b_L1_Mu15_HTT100er);
   T1->SetBranchAddress("L1_Mu18_HTT100er", &L1_Mu18_HTT100er, &b_L1_Mu18_HTT100er);
   T1->SetBranchAddress("L1_Mu18_Jet24er2p7", &L1_Mu18_Jet24er2p7, &b_L1_Mu18_Jet24er2p7);
   T1->SetBranchAddress("L1_Mu18er2p1_IsoTau26er2p1", &L1_Mu18er2p1_IsoTau26er2p1, &b_L1_Mu18er2p1_IsoTau26er2p1);
   T1->SetBranchAddress("L1_Mu18er2p1_Tau24er2p1", &L1_Mu18er2p1_Tau24er2p1, &b_L1_Mu18er2p1_Tau24er2p1);
   T1->SetBranchAddress("L1_Mu20_EG10", &L1_Mu20_EG10, &b_L1_Mu20_EG10);
   T1->SetBranchAddress("L1_Mu20_EG17", &L1_Mu20_EG17, &b_L1_Mu20_EG17);
   T1->SetBranchAddress("L1_Mu20_LooseIsoEG6", &L1_Mu20_LooseIsoEG6, &b_L1_Mu20_LooseIsoEG6);
   T1->SetBranchAddress("L1_Mu20er2p1_IsoTau26er2p1", &L1_Mu20er2p1_IsoTau26er2p1, &b_L1_Mu20er2p1_IsoTau26er2p1);
   T1->SetBranchAddress("L1_Mu20er2p1_IsoTau27er2p1", &L1_Mu20er2p1_IsoTau27er2p1, &b_L1_Mu20er2p1_IsoTau27er2p1);
   T1->SetBranchAddress("L1_Mu22er2p1_IsoTau28er2p1", &L1_Mu22er2p1_IsoTau28er2p1, &b_L1_Mu22er2p1_IsoTau28er2p1);
   T1->SetBranchAddress("L1_Mu22er2p1_IsoTau30er2p1", &L1_Mu22er2p1_IsoTau30er2p1, &b_L1_Mu22er2p1_IsoTau30er2p1);
   T1->SetBranchAddress("L1_Mu22er2p1_IsoTau32er2p1", &L1_Mu22er2p1_IsoTau32er2p1, &b_L1_Mu22er2p1_IsoTau32er2p1);
   T1->SetBranchAddress("L1_Mu22er2p1_IsoTau33er2p1", &L1_Mu22er2p1_IsoTau33er2p1, &b_L1_Mu22er2p1_IsoTau33er2p1);
   T1->SetBranchAddress("L1_Mu22er2p1_IsoTau34er2p1", &L1_Mu22er2p1_IsoTau34er2p1, &b_L1_Mu22er2p1_IsoTau34er2p1);
   T1->SetBranchAddress("L1_Mu22er2p1_IsoTau35er2p1", &L1_Mu22er2p1_IsoTau35er2p1, &b_L1_Mu22er2p1_IsoTau35er2p1);
   T1->SetBranchAddress("L1_Mu22er2p1_IsoTau36er2p1", &L1_Mu22er2p1_IsoTau36er2p1, &b_L1_Mu22er2p1_IsoTau36er2p1);
   T1->SetBranchAddress("L1_Mu22er2p1_IsoTau38er2p1", &L1_Mu22er2p1_IsoTau38er2p1, &b_L1_Mu22er2p1_IsoTau38er2p1);
   T1->SetBranchAddress("L1_Mu22er2p1_IsoTau40er2p1", &L1_Mu22er2p1_IsoTau40er2p1, &b_L1_Mu22er2p1_IsoTau40er2p1);
   T1->SetBranchAddress("L1_Mu22er2p1_Tau50er2p1", &L1_Mu22er2p1_Tau50er2p1, &b_L1_Mu22er2p1_Tau50er2p1);
   T1->SetBranchAddress("L1_Mu22er2p1_Tau70er2p1", &L1_Mu22er2p1_Tau70er2p1, &b_L1_Mu22er2p1_Tau70er2p1);
   T1->SetBranchAddress("L1_Mu23_EG10", &L1_Mu23_EG10, &b_L1_Mu23_EG10);
   T1->SetBranchAddress("L1_Mu23_LooseIsoEG10", &L1_Mu23_LooseIsoEG10, &b_L1_Mu23_LooseIsoEG10);
   T1->SetBranchAddress("L1_Mu3_Jet120er2p7_dEta_Max0p4_dPhi_Max0p4", &L1_Mu3_Jet120er2p7_dEta_Max0p4_dPhi_Max0p4, &b_L1_Mu3_Jet120er2p7_dEta_Max0p4_dPhi_Max0p4);
   T1->SetBranchAddress("L1_Mu3_Jet16er2p7_dEta_Max0p4_dPhi_Max0p4", &L1_Mu3_Jet16er2p7_dEta_Max0p4_dPhi_Max0p4, &b_L1_Mu3_Jet16er2p7_dEta_Max0p4_dPhi_Max0p4);
   T1->SetBranchAddress("L1_Mu3_Jet30er2p5", &L1_Mu3_Jet30er2p5, &b_L1_Mu3_Jet30er2p5);
   T1->SetBranchAddress("L1_Mu3_Jet60er2p7_dEta_Max0p4_dPhi_Max0p4", &L1_Mu3_Jet60er2p7_dEta_Max0p4_dPhi_Max0p4, &b_L1_Mu3_Jet60er2p7_dEta_Max0p4_dPhi_Max0p4);
   T1->SetBranchAddress("L1_Mu5_EG15", &L1_Mu5_EG15, &b_L1_Mu5_EG15);
   T1->SetBranchAddress("L1_Mu5_EG20", &L1_Mu5_EG20, &b_L1_Mu5_EG20);
   T1->SetBranchAddress("L1_Mu5_EG23", &L1_Mu5_EG23, &b_L1_Mu5_EG23);
   T1->SetBranchAddress("L1_Mu5_LooseIsoEG18", &L1_Mu5_LooseIsoEG18, &b_L1_Mu5_LooseIsoEG18);
   T1->SetBranchAddress("L1_Mu5_LooseIsoEG20", &L1_Mu5_LooseIsoEG20, &b_L1_Mu5_LooseIsoEG20);
   T1->SetBranchAddress("L1_Mu6_DoubleEG10", &L1_Mu6_DoubleEG10, &b_L1_Mu6_DoubleEG10);
   T1->SetBranchAddress("L1_Mu6_DoubleEG17", &L1_Mu6_DoubleEG17, &b_L1_Mu6_DoubleEG17);
   T1->SetBranchAddress("L1_Mu6_HTT200er", &L1_Mu6_HTT200er, &b_L1_Mu6_HTT200er);
   T1->SetBranchAddress("L1_Mu6_HTT240er", &L1_Mu6_HTT240er, &b_L1_Mu6_HTT240er);
   T1->SetBranchAddress("L1_Mu6_HTT250er", &L1_Mu6_HTT250er, &b_L1_Mu6_HTT250er);
   T1->SetBranchAddress("L1_Mu7_EG23", &L1_Mu7_EG23, &b_L1_Mu7_EG23);
   T1->SetBranchAddress("L1_Mu7_LooseIsoEG20", &L1_Mu7_LooseIsoEG20, &b_L1_Mu7_LooseIsoEG20);
   T1->SetBranchAddress("L1_Mu7_LooseIsoEG23", &L1_Mu7_LooseIsoEG23, &b_L1_Mu7_LooseIsoEG23);
   T1->SetBranchAddress("L1_Mu8_HTT150er", &L1_Mu8_HTT150er, &b_L1_Mu8_HTT150er);
   T1->SetBranchAddress("L1_NotBptxOR", &L1_NotBptxOR, &b_L1_NotBptxOR);
   T1->SetBranchAddress("L1_QuadJet36er2p7_IsoTau52er2p1", &L1_QuadJet36er2p7_IsoTau52er2p1, &b_L1_QuadJet36er2p7_IsoTau52er2p1);
   T1->SetBranchAddress("L1_QuadJet36er2p7_Tau52", &L1_QuadJet36er2p7_Tau52, &b_L1_QuadJet36er2p7_Tau52);
   T1->SetBranchAddress("L1_QuadJet40er2p7", &L1_QuadJet40er2p7, &b_L1_QuadJet40er2p7);
   T1->SetBranchAddress("L1_QuadJet50er2p7", &L1_QuadJet50er2p7, &b_L1_QuadJet50er2p7);
   T1->SetBranchAddress("L1_QuadJet60er2p7", &L1_QuadJet60er2p7, &b_L1_QuadJet60er2p7);
   T1->SetBranchAddress("L1_QuadMu0", &L1_QuadMu0, &b_L1_QuadMu0);
   T1->SetBranchAddress("L1_SingleEG10", &L1_SingleEG10, &b_L1_SingleEG10);
   T1->SetBranchAddress("L1_SingleEG15", &L1_SingleEG15, &b_L1_SingleEG15);
   T1->SetBranchAddress("L1_SingleEG18", &L1_SingleEG18, &b_L1_SingleEG18);
   T1->SetBranchAddress("L1_SingleEG24", &L1_SingleEG24, &b_L1_SingleEG24);
   T1->SetBranchAddress("L1_SingleEG26", &L1_SingleEG26, &b_L1_SingleEG26);
   T1->SetBranchAddress("L1_SingleEG28", &L1_SingleEG28, &b_L1_SingleEG28);
   T1->SetBranchAddress("L1_SingleEG2_BptxAND", &L1_SingleEG2_BptxAND, &b_L1_SingleEG2_BptxAND);
   T1->SetBranchAddress("L1_SingleEG30", &L1_SingleEG30, &b_L1_SingleEG30);
   T1->SetBranchAddress("L1_SingleEG32", &L1_SingleEG32, &b_L1_SingleEG32);
   T1->SetBranchAddress("L1_SingleEG34", &L1_SingleEG34, &b_L1_SingleEG34);
   T1->SetBranchAddress("L1_SingleEG34er2p1", &L1_SingleEG34er2p1, &b_L1_SingleEG34er2p1);
   T1->SetBranchAddress("L1_SingleEG36", &L1_SingleEG36, &b_L1_SingleEG36);
   T1->SetBranchAddress("L1_SingleEG36er2p1", &L1_SingleEG36er2p1, &b_L1_SingleEG36er2p1);
   T1->SetBranchAddress("L1_SingleEG38", &L1_SingleEG38, &b_L1_SingleEG38);
   T1->SetBranchAddress("L1_SingleEG38er2p1", &L1_SingleEG38er2p1, &b_L1_SingleEG38er2p1);
   T1->SetBranchAddress("L1_SingleEG40", &L1_SingleEG40, &b_L1_SingleEG40);
   T1->SetBranchAddress("L1_SingleEG42", &L1_SingleEG42, &b_L1_SingleEG42);
   T1->SetBranchAddress("L1_SingleEG45", &L1_SingleEG45, &b_L1_SingleEG45);
   T1->SetBranchAddress("L1_SingleEG5", &L1_SingleEG5, &b_L1_SingleEG5);
   T1->SetBranchAddress("L1_SingleEG50", &L1_SingleEG50, &b_L1_SingleEG50);
   T1->SetBranchAddress("L1_SingleIsoEG18", &L1_SingleIsoEG18, &b_L1_SingleIsoEG18);
   T1->SetBranchAddress("L1_SingleIsoEG18er2p1", &L1_SingleIsoEG18er2p1, &b_L1_SingleIsoEG18er2p1);
   T1->SetBranchAddress("L1_SingleIsoEG20", &L1_SingleIsoEG20, &b_L1_SingleIsoEG20);
   T1->SetBranchAddress("L1_SingleIsoEG20er2p1", &L1_SingleIsoEG20er2p1, &b_L1_SingleIsoEG20er2p1);
   T1->SetBranchAddress("L1_SingleIsoEG22", &L1_SingleIsoEG22, &b_L1_SingleIsoEG22);
   T1->SetBranchAddress("L1_SingleIsoEG22er2p1", &L1_SingleIsoEG22er2p1, &b_L1_SingleIsoEG22er2p1);
   T1->SetBranchAddress("L1_SingleIsoEG24", &L1_SingleIsoEG24, &b_L1_SingleIsoEG24);
   T1->SetBranchAddress("L1_SingleIsoEG24er2p1", &L1_SingleIsoEG24er2p1, &b_L1_SingleIsoEG24er2p1);
   T1->SetBranchAddress("L1_SingleIsoEG26", &L1_SingleIsoEG26, &b_L1_SingleIsoEG26);
   T1->SetBranchAddress("L1_SingleIsoEG26er2p1", &L1_SingleIsoEG26er2p1, &b_L1_SingleIsoEG26er2p1);
   T1->SetBranchAddress("L1_SingleIsoEG28", &L1_SingleIsoEG28, &b_L1_SingleIsoEG28);
   T1->SetBranchAddress("L1_SingleIsoEG28er2p1", &L1_SingleIsoEG28er2p1, &b_L1_SingleIsoEG28er2p1);
   T1->SetBranchAddress("L1_SingleIsoEG30", &L1_SingleIsoEG30, &b_L1_SingleIsoEG30);
   T1->SetBranchAddress("L1_SingleIsoEG30er2p1", &L1_SingleIsoEG30er2p1, &b_L1_SingleIsoEG30er2p1);
   T1->SetBranchAddress("L1_SingleIsoEG32", &L1_SingleIsoEG32, &b_L1_SingleIsoEG32);
   T1->SetBranchAddress("L1_SingleIsoEG32er2p1", &L1_SingleIsoEG32er2p1, &b_L1_SingleIsoEG32er2p1);
   T1->SetBranchAddress("L1_SingleIsoEG33er2p1", &L1_SingleIsoEG33er2p1, &b_L1_SingleIsoEG33er2p1);
   T1->SetBranchAddress("L1_SingleIsoEG34", &L1_SingleIsoEG34, &b_L1_SingleIsoEG34);
   T1->SetBranchAddress("L1_SingleIsoEG34er2p1", &L1_SingleIsoEG34er2p1, &b_L1_SingleIsoEG34er2p1);
   T1->SetBranchAddress("L1_SingleIsoEG35", &L1_SingleIsoEG35, &b_L1_SingleIsoEG35);
   T1->SetBranchAddress("L1_SingleIsoEG35er2p1", &L1_SingleIsoEG35er2p1, &b_L1_SingleIsoEG35er2p1);
   T1->SetBranchAddress("L1_SingleIsoEG36", &L1_SingleIsoEG36, &b_L1_SingleIsoEG36);
   T1->SetBranchAddress("L1_SingleIsoEG36er2p1", &L1_SingleIsoEG36er2p1, &b_L1_SingleIsoEG36er2p1);
   T1->SetBranchAddress("L1_SingleIsoEG37", &L1_SingleIsoEG37, &b_L1_SingleIsoEG37);
   T1->SetBranchAddress("L1_SingleIsoEG38", &L1_SingleIsoEG38, &b_L1_SingleIsoEG38);
   T1->SetBranchAddress("L1_SingleIsoEG38er2p1", &L1_SingleIsoEG38er2p1, &b_L1_SingleIsoEG38er2p1);
   T1->SetBranchAddress("L1_SingleIsoEG40", &L1_SingleIsoEG40, &b_L1_SingleIsoEG40);
   T1->SetBranchAddress("L1_SingleIsoEG40er2p1", &L1_SingleIsoEG40er2p1, &b_L1_SingleIsoEG40er2p1);
   T1->SetBranchAddress("L1_SingleJet120", &L1_SingleJet120, &b_L1_SingleJet120);
   T1->SetBranchAddress("L1_SingleJet120_FWD", &L1_SingleJet120_FWD, &b_L1_SingleJet120_FWD);
   T1->SetBranchAddress("L1_SingleJet12_BptxAND", &L1_SingleJet12_BptxAND, &b_L1_SingleJet12_BptxAND);
   T1->SetBranchAddress("L1_SingleJet140", &L1_SingleJet140, &b_L1_SingleJet140);
   T1->SetBranchAddress("L1_SingleJet150", &L1_SingleJet150, &b_L1_SingleJet150);
   T1->SetBranchAddress("L1_SingleJet16", &L1_SingleJet16, &b_L1_SingleJet16);
   T1->SetBranchAddress("L1_SingleJet160", &L1_SingleJet160, &b_L1_SingleJet160);
   T1->SetBranchAddress("L1_SingleJet170", &L1_SingleJet170, &b_L1_SingleJet170);
   T1->SetBranchAddress("L1_SingleJet180", &L1_SingleJet180, &b_L1_SingleJet180);
   T1->SetBranchAddress("L1_SingleJet20", &L1_SingleJet20, &b_L1_SingleJet20);
   T1->SetBranchAddress("L1_SingleJet200", &L1_SingleJet200, &b_L1_SingleJet200);
   T1->SetBranchAddress("L1_SingleJet20er2p7_NotBptxOR", &L1_SingleJet20er2p7_NotBptxOR, &b_L1_SingleJet20er2p7_NotBptxOR);
   T1->SetBranchAddress("L1_SingleJet20er2p7_NotBptxOR_3BX", &L1_SingleJet20er2p7_NotBptxOR_3BX, &b_L1_SingleJet20er2p7_NotBptxOR_3BX);
   T1->SetBranchAddress("L1_SingleJet35", &L1_SingleJet35, &b_L1_SingleJet35);
   T1->SetBranchAddress("L1_SingleJet35_FWD", &L1_SingleJet35_FWD, &b_L1_SingleJet35_FWD);
   T1->SetBranchAddress("L1_SingleJet35_HFm", &L1_SingleJet35_HFm, &b_L1_SingleJet35_HFm);
   T1->SetBranchAddress("L1_SingleJet35_HFp", &L1_SingleJet35_HFp, &b_L1_SingleJet35_HFp);
   T1->SetBranchAddress("L1_SingleJet43er2p7_NotBptxOR_3BX", &L1_SingleJet43er2p7_NotBptxOR_3BX, &b_L1_SingleJet43er2p7_NotBptxOR_3BX);
   T1->SetBranchAddress("L1_SingleJet46er2p7_NotBptxOR_3BX", &L1_SingleJet46er2p7_NotBptxOR_3BX, &b_L1_SingleJet46er2p7_NotBptxOR_3BX);
   T1->SetBranchAddress("L1_SingleJet60", &L1_SingleJet60, &b_L1_SingleJet60);
   T1->SetBranchAddress("L1_SingleJet60_FWD", &L1_SingleJet60_FWD, &b_L1_SingleJet60_FWD);
   T1->SetBranchAddress("L1_SingleJet60_HFm", &L1_SingleJet60_HFm, &b_L1_SingleJet60_HFm);
   T1->SetBranchAddress("L1_SingleJet60_HFp", &L1_SingleJet60_HFp, &b_L1_SingleJet60_HFp);
   T1->SetBranchAddress("L1_SingleJet90", &L1_SingleJet90, &b_L1_SingleJet90);
   T1->SetBranchAddress("L1_SingleJet90_FWD", &L1_SingleJet90_FWD, &b_L1_SingleJet90_FWD);
   T1->SetBranchAddress("L1_SingleMu0_BMTF", &L1_SingleMu0_BMTF, &b_L1_SingleMu0_BMTF);
   T1->SetBranchAddress("L1_SingleMu0_EMTF", &L1_SingleMu0_EMTF, &b_L1_SingleMu0_EMTF);
   T1->SetBranchAddress("L1_SingleMu0_OMTF", &L1_SingleMu0_OMTF, &b_L1_SingleMu0_OMTF);
   T1->SetBranchAddress("L1_SingleMu10_LowQ", &L1_SingleMu10_LowQ, &b_L1_SingleMu10_LowQ);
   T1->SetBranchAddress("L1_SingleMu11_LowQ", &L1_SingleMu11_LowQ, &b_L1_SingleMu11_LowQ);
   T1->SetBranchAddress("L1_SingleMu12_LowQ_BMTF", &L1_SingleMu12_LowQ_BMTF, &b_L1_SingleMu12_LowQ_BMTF);
   T1->SetBranchAddress("L1_SingleMu12_LowQ_EMTF", &L1_SingleMu12_LowQ_EMTF, &b_L1_SingleMu12_LowQ_EMTF);
   T1->SetBranchAddress("L1_SingleMu12_LowQ_OMTF", &L1_SingleMu12_LowQ_OMTF, &b_L1_SingleMu12_LowQ_OMTF);
   T1->SetBranchAddress("L1_SingleMu14er2p1", &L1_SingleMu14er2p1, &b_L1_SingleMu14er2p1);
   T1->SetBranchAddress("L1_SingleMu16", &L1_SingleMu16, &b_L1_SingleMu16);
   T1->SetBranchAddress("L1_SingleMu16er2p1", &L1_SingleMu16er2p1, &b_L1_SingleMu16er2p1);
   T1->SetBranchAddress("L1_SingleMu18", &L1_SingleMu18, &b_L1_SingleMu18);
   T1->SetBranchAddress("L1_SingleMu18er2p1", &L1_SingleMu18er2p1, &b_L1_SingleMu18er2p1);
   T1->SetBranchAddress("L1_SingleMu20", &L1_SingleMu20, &b_L1_SingleMu20);
   T1->SetBranchAddress("L1_SingleMu20er2p1", &L1_SingleMu20er2p1, &b_L1_SingleMu20er2p1);
   T1->SetBranchAddress("L1_SingleMu22", &L1_SingleMu22, &b_L1_SingleMu22);
   T1->SetBranchAddress("L1_SingleMu22_BMTF", &L1_SingleMu22_BMTF, &b_L1_SingleMu22_BMTF);
   T1->SetBranchAddress("L1_SingleMu22_EMTF", &L1_SingleMu22_EMTF, &b_L1_SingleMu22_EMTF);
   T1->SetBranchAddress("L1_SingleMu22_OMTF", &L1_SingleMu22_OMTF, &b_L1_SingleMu22_OMTF);
   T1->SetBranchAddress("L1_SingleMu22er2p1", &L1_SingleMu22er2p1, &b_L1_SingleMu22er2p1);
   T1->SetBranchAddress("L1_SingleMu25", &L1_SingleMu25, &b_L1_SingleMu25);
   T1->SetBranchAddress("L1_SingleMu3", &L1_SingleMu3, &b_L1_SingleMu3);
   T1->SetBranchAddress("L1_SingleMu30", &L1_SingleMu30, &b_L1_SingleMu30);
   T1->SetBranchAddress("L1_SingleMu5", &L1_SingleMu5, &b_L1_SingleMu5);
   T1->SetBranchAddress("L1_SingleMu7", &L1_SingleMu7, &b_L1_SingleMu7);
   T1->SetBranchAddress("L1_SingleMuCosmics", &L1_SingleMuCosmics, &b_L1_SingleMuCosmics);
   T1->SetBranchAddress("L1_SingleMuCosmics_BMTF", &L1_SingleMuCosmics_BMTF, &b_L1_SingleMuCosmics_BMTF);
   T1->SetBranchAddress("L1_SingleMuCosmics_EMTF", &L1_SingleMuCosmics_EMTF, &b_L1_SingleMuCosmics_EMTF);
   T1->SetBranchAddress("L1_SingleMuCosmics_OMTF", &L1_SingleMuCosmics_OMTF, &b_L1_SingleMuCosmics_OMTF);
   T1->SetBranchAddress("L1_SingleMuOpen", &L1_SingleMuOpen, &b_L1_SingleMuOpen);
   T1->SetBranchAddress("L1_SingleMuOpen_NotBptxOR", &L1_SingleMuOpen_NotBptxOR, &b_L1_SingleMuOpen_NotBptxOR);
   T1->SetBranchAddress("L1_SingleMuOpen_NotBptxOR_3BX", &L1_SingleMuOpen_NotBptxOR_3BX, &b_L1_SingleMuOpen_NotBptxOR_3BX);
   T1->SetBranchAddress("L1_SingleTau100er2p1", &L1_SingleTau100er2p1, &b_L1_SingleTau100er2p1);
   T1->SetBranchAddress("L1_SingleTau120er2p1", &L1_SingleTau120er2p1, &b_L1_SingleTau120er2p1);
   T1->SetBranchAddress("L1_SingleTau130er2p1", &L1_SingleTau130er2p1, &b_L1_SingleTau130er2p1);
   T1->SetBranchAddress("L1_SingleTau140er2p1", &L1_SingleTau140er2p1, &b_L1_SingleTau140er2p1);
   T1->SetBranchAddress("L1_SingleTau20", &L1_SingleTau20, &b_L1_SingleTau20);
   T1->SetBranchAddress("L1_SingleTau80er2p1", &L1_SingleTau80er2p1, &b_L1_SingleTau80er2p1);
   T1->SetBranchAddress("L1_TripleEG_14_10_8", &L1_TripleEG_14_10_8, &b_L1_TripleEG_14_10_8);
   T1->SetBranchAddress("L1_TripleEG_18_17_8", &L1_TripleEG_18_17_8, &b_L1_TripleEG_18_17_8);
   T1->SetBranchAddress("L1_TripleEG_LooseIso20_10_5", &L1_TripleEG_LooseIso20_10_5, &b_L1_TripleEG_LooseIso20_10_5);
   T1->SetBranchAddress("L1_TripleJet_100_85_72_VBF", &L1_TripleJet_100_85_72_VBF, &b_L1_TripleJet_100_85_72_VBF);
   T1->SetBranchAddress("L1_TripleJet_105_85_76_VBF", &L1_TripleJet_105_85_76_VBF, &b_L1_TripleJet_105_85_76_VBF);
   T1->SetBranchAddress("L1_TripleJet_84_68_48_VBF", &L1_TripleJet_84_68_48_VBF, &b_L1_TripleJet_84_68_48_VBF);
   T1->SetBranchAddress("L1_TripleJet_88_72_56_VBF", &L1_TripleJet_88_72_56_VBF, &b_L1_TripleJet_88_72_56_VBF);
   T1->SetBranchAddress("L1_TripleJet_92_76_64_VBF", &L1_TripleJet_92_76_64_VBF, &b_L1_TripleJet_92_76_64_VBF);
   T1->SetBranchAddress("L1_TripleJet_98_83_71_VBF", &L1_TripleJet_98_83_71_VBF, &b_L1_TripleJet_98_83_71_VBF);
   T1->SetBranchAddress("L1_TripleMu0", &L1_TripleMu0, &b_L1_TripleMu0);
   T1->SetBranchAddress("L1_TripleMu0_OQ", &L1_TripleMu0_OQ, &b_L1_TripleMu0_OQ);
   T1->SetBranchAddress("L1_TripleMu3", &L1_TripleMu3, &b_L1_TripleMu3);
   T1->SetBranchAddress("L1_TripleMu3_SQ", &L1_TripleMu3_SQ, &b_L1_TripleMu3_SQ);
   T1->SetBranchAddress("L1_TripleMu_4_4_4", &L1_TripleMu_4_4_4, &b_L1_TripleMu_4_4_4);
   T1->SetBranchAddress("L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17", &L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17, &b_L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17);
   T1->SetBranchAddress("L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_8to14", &L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_8to14, &b_L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_8to14);
   T1->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0OQ", &L1_TripleMu_5SQ_3SQ_0OQ, &b_L1_TripleMu_5SQ_3SQ_0OQ);
   T1->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9, &b_L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9);
   T1->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9, &b_L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9);
   T1->SetBranchAddress("L1_TripleMu_5_0_0", &L1_TripleMu_5_0_0, &b_L1_TripleMu_5_0_0);
   T1->SetBranchAddress("L1_TripleMu_5_3_3", &L1_TripleMu_5_3_3, &b_L1_TripleMu_5_3_3);
   T1->SetBranchAddress("L1_TripleMu_5_3p5_2p5", &L1_TripleMu_5_3p5_2p5, &b_L1_TripleMu_5_3p5_2p5);
   T1->SetBranchAddress("L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17, &b_L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17);
   T1->SetBranchAddress("L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17, &b_L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17);
   T1->SetBranchAddress("L1_TripleMu_5_5_3", &L1_TripleMu_5_5_3, &b_L1_TripleMu_5_5_3);
   T1->SetBranchAddress("L1_UnpairedBunchBptxMinus", &L1_UnpairedBunchBptxMinus, &b_L1_UnpairedBunchBptxMinus);
   T1->SetBranchAddress("L1_UnpairedBunchBptxPlus", &L1_UnpairedBunchBptxPlus, &b_L1_UnpairedBunchBptxPlus);
   T1->SetBranchAddress("L1_ZeroBias", &L1_ZeroBias, &b_L1_ZeroBias);
   T1->SetBranchAddress("L1_ZeroBias_copy", &L1_ZeroBias_copy, &b_L1_ZeroBias_copy);
   T1->SetBranchAddress("Jet_pt_raw", Jet_pt_raw, &b_Jet_pt_raw);
   T1->SetBranchAddress("Jet_pt_nom", Jet_pt_nom, &b_Jet_pt_nom);
   T1->SetBranchAddress("Jet_mass_raw", Jet_mass_raw, &b_Jet_mass_raw);
   T1->SetBranchAddress("Jet_mass_nom", Jet_mass_nom, &b_Jet_mass_nom);
   T1->SetBranchAddress("Jet_corr_JEC", Jet_corr_JEC, &b_Jet_corr_JEC);
   T1->SetBranchAddress("Jet_corr_JER", Jet_corr_JER, &b_Jet_corr_JER);
   T1->SetBranchAddress("MET_pt_nom", &MET_pt_nom, &b_MET_pt_nom);
   T1->SetBranchAddress("MET_phi_nom", &MET_phi_nom, &b_MET_phi_nom);
   T1->SetBranchAddress("MET_pt_jer", &MET_pt_jer, &b_MET_pt_jer);
   T1->SetBranchAddress("MET_phi_jer", &MET_phi_jer, &b_MET_phi_jer);
   T1->SetBranchAddress("Jet_pt_jerUp", Jet_pt_jerUp, &b_Jet_pt_jerUp);
   T1->SetBranchAddress("Jet_mass_jerUp", Jet_mass_jerUp, &b_Jet_mass_jerUp);
   T1->SetBranchAddress("MET_pt_jerUp", &MET_pt_jerUp, &b_MET_pt_jerUp);
   T1->SetBranchAddress("MET_phi_jerUp", &MET_phi_jerUp, &b_MET_phi_jerUp);
   T1->SetBranchAddress("Jet_pt_jesAbsoluteStatUp", Jet_pt_jesAbsoluteStatUp, &b_Jet_pt_jesAbsoluteStatUp);
   T1->SetBranchAddress("Jet_mass_jesAbsoluteStatUp", Jet_mass_jesAbsoluteStatUp, &b_Jet_mass_jesAbsoluteStatUp);
   T1->SetBranchAddress("MET_pt_jesAbsoluteStatUp", &MET_pt_jesAbsoluteStatUp, &b_MET_pt_jesAbsoluteStatUp);
   T1->SetBranchAddress("MET_phi_jesAbsoluteStatUp", &MET_phi_jesAbsoluteStatUp, &b_MET_phi_jesAbsoluteStatUp);
   T1->SetBranchAddress("Jet_pt_jesAbsoluteScaleUp", Jet_pt_jesAbsoluteScaleUp, &b_Jet_pt_jesAbsoluteScaleUp);
   T1->SetBranchAddress("Jet_mass_jesAbsoluteScaleUp", Jet_mass_jesAbsoluteScaleUp, &b_Jet_mass_jesAbsoluteScaleUp);
   T1->SetBranchAddress("MET_pt_jesAbsoluteScaleUp", &MET_pt_jesAbsoluteScaleUp, &b_MET_pt_jesAbsoluteScaleUp);
   T1->SetBranchAddress("MET_phi_jesAbsoluteScaleUp", &MET_phi_jesAbsoluteScaleUp, &b_MET_phi_jesAbsoluteScaleUp);
   T1->SetBranchAddress("Jet_pt_jesAbsoluteFlavMapUp", Jet_pt_jesAbsoluteFlavMapUp, &b_Jet_pt_jesAbsoluteFlavMapUp);
   T1->SetBranchAddress("Jet_mass_jesAbsoluteFlavMapUp", Jet_mass_jesAbsoluteFlavMapUp, &b_Jet_mass_jesAbsoluteFlavMapUp);
   T1->SetBranchAddress("MET_pt_jesAbsoluteFlavMapUp", &MET_pt_jesAbsoluteFlavMapUp, &b_MET_pt_jesAbsoluteFlavMapUp);
   T1->SetBranchAddress("MET_phi_jesAbsoluteFlavMapUp", &MET_phi_jesAbsoluteFlavMapUp, &b_MET_phi_jesAbsoluteFlavMapUp);
   T1->SetBranchAddress("Jet_pt_jesAbsoluteMPFBiasUp", Jet_pt_jesAbsoluteMPFBiasUp, &b_Jet_pt_jesAbsoluteMPFBiasUp);
   T1->SetBranchAddress("Jet_mass_jesAbsoluteMPFBiasUp", Jet_mass_jesAbsoluteMPFBiasUp, &b_Jet_mass_jesAbsoluteMPFBiasUp);
   T1->SetBranchAddress("MET_pt_jesAbsoluteMPFBiasUp", &MET_pt_jesAbsoluteMPFBiasUp, &b_MET_pt_jesAbsoluteMPFBiasUp);
   T1->SetBranchAddress("MET_phi_jesAbsoluteMPFBiasUp", &MET_phi_jesAbsoluteMPFBiasUp, &b_MET_phi_jesAbsoluteMPFBiasUp);
   T1->SetBranchAddress("Jet_pt_jesFragmentationUp", Jet_pt_jesFragmentationUp, &b_Jet_pt_jesFragmentationUp);
   T1->SetBranchAddress("Jet_mass_jesFragmentationUp", Jet_mass_jesFragmentationUp, &b_Jet_mass_jesFragmentationUp);
   T1->SetBranchAddress("MET_pt_jesFragmentationUp", &MET_pt_jesFragmentationUp, &b_MET_pt_jesFragmentationUp);
   T1->SetBranchAddress("MET_phi_jesFragmentationUp", &MET_phi_jesFragmentationUp, &b_MET_phi_jesFragmentationUp);
   T1->SetBranchAddress("Jet_pt_jesSinglePionECALUp", Jet_pt_jesSinglePionECALUp, &b_Jet_pt_jesSinglePionECALUp);
   T1->SetBranchAddress("Jet_mass_jesSinglePionECALUp", Jet_mass_jesSinglePionECALUp, &b_Jet_mass_jesSinglePionECALUp);
   T1->SetBranchAddress("MET_pt_jesSinglePionECALUp", &MET_pt_jesSinglePionECALUp, &b_MET_pt_jesSinglePionECALUp);
   T1->SetBranchAddress("MET_phi_jesSinglePionECALUp", &MET_phi_jesSinglePionECALUp, &b_MET_phi_jesSinglePionECALUp);
   T1->SetBranchAddress("Jet_pt_jesSinglePionHCALUp", Jet_pt_jesSinglePionHCALUp, &b_Jet_pt_jesSinglePionHCALUp);
   T1->SetBranchAddress("Jet_mass_jesSinglePionHCALUp", Jet_mass_jesSinglePionHCALUp, &b_Jet_mass_jesSinglePionHCALUp);
   T1->SetBranchAddress("MET_pt_jesSinglePionHCALUp", &MET_pt_jesSinglePionHCALUp, &b_MET_pt_jesSinglePionHCALUp);
   T1->SetBranchAddress("MET_phi_jesSinglePionHCALUp", &MET_phi_jesSinglePionHCALUp, &b_MET_phi_jesSinglePionHCALUp);
   T1->SetBranchAddress("Jet_pt_jesFlavorQCDUp", Jet_pt_jesFlavorQCDUp, &b_Jet_pt_jesFlavorQCDUp);
   T1->SetBranchAddress("Jet_mass_jesFlavorQCDUp", Jet_mass_jesFlavorQCDUp, &b_Jet_mass_jesFlavorQCDUp);
   T1->SetBranchAddress("MET_pt_jesFlavorQCDUp", &MET_pt_jesFlavorQCDUp, &b_MET_pt_jesFlavorQCDUp);
   T1->SetBranchAddress("MET_phi_jesFlavorQCDUp", &MET_phi_jesFlavorQCDUp, &b_MET_phi_jesFlavorQCDUp);
   T1->SetBranchAddress("Jet_pt_jesTimePtEtaUp", Jet_pt_jesTimePtEtaUp, &b_Jet_pt_jesTimePtEtaUp);
   T1->SetBranchAddress("Jet_mass_jesTimePtEtaUp", Jet_mass_jesTimePtEtaUp, &b_Jet_mass_jesTimePtEtaUp);
   T1->SetBranchAddress("MET_pt_jesTimePtEtaUp", &MET_pt_jesTimePtEtaUp, &b_MET_pt_jesTimePtEtaUp);
   T1->SetBranchAddress("MET_phi_jesTimePtEtaUp", &MET_phi_jesTimePtEtaUp, &b_MET_phi_jesTimePtEtaUp);
   T1->SetBranchAddress("Jet_pt_jesRelativeJEREC1Up", Jet_pt_jesRelativeJEREC1Up, &b_Jet_pt_jesRelativeJEREC1Up);
   T1->SetBranchAddress("Jet_mass_jesRelativeJEREC1Up", Jet_mass_jesRelativeJEREC1Up, &b_Jet_mass_jesRelativeJEREC1Up);
   T1->SetBranchAddress("MET_pt_jesRelativeJEREC1Up", &MET_pt_jesRelativeJEREC1Up, &b_MET_pt_jesRelativeJEREC1Up);
   T1->SetBranchAddress("MET_phi_jesRelativeJEREC1Up", &MET_phi_jesRelativeJEREC1Up, &b_MET_phi_jesRelativeJEREC1Up);
   T1->SetBranchAddress("Jet_pt_jesRelativeJEREC2Up", Jet_pt_jesRelativeJEREC2Up, &b_Jet_pt_jesRelativeJEREC2Up);
   T1->SetBranchAddress("Jet_mass_jesRelativeJEREC2Up", Jet_mass_jesRelativeJEREC2Up, &b_Jet_mass_jesRelativeJEREC2Up);
   T1->SetBranchAddress("MET_pt_jesRelativeJEREC2Up", &MET_pt_jesRelativeJEREC2Up, &b_MET_pt_jesRelativeJEREC2Up);
   T1->SetBranchAddress("MET_phi_jesRelativeJEREC2Up", &MET_phi_jesRelativeJEREC2Up, &b_MET_phi_jesRelativeJEREC2Up);
   T1->SetBranchAddress("Jet_pt_jesRelativeJERHFUp", Jet_pt_jesRelativeJERHFUp, &b_Jet_pt_jesRelativeJERHFUp);
   T1->SetBranchAddress("Jet_mass_jesRelativeJERHFUp", Jet_mass_jesRelativeJERHFUp, &b_Jet_mass_jesRelativeJERHFUp);
   T1->SetBranchAddress("MET_pt_jesRelativeJERHFUp", &MET_pt_jesRelativeJERHFUp, &b_MET_pt_jesRelativeJERHFUp);
   T1->SetBranchAddress("MET_phi_jesRelativeJERHFUp", &MET_phi_jesRelativeJERHFUp, &b_MET_phi_jesRelativeJERHFUp);
   T1->SetBranchAddress("Jet_pt_jesRelativePtBBUp", Jet_pt_jesRelativePtBBUp, &b_Jet_pt_jesRelativePtBBUp);
   T1->SetBranchAddress("Jet_mass_jesRelativePtBBUp", Jet_mass_jesRelativePtBBUp, &b_Jet_mass_jesRelativePtBBUp);
   T1->SetBranchAddress("MET_pt_jesRelativePtBBUp", &MET_pt_jesRelativePtBBUp, &b_MET_pt_jesRelativePtBBUp);
   T1->SetBranchAddress("MET_phi_jesRelativePtBBUp", &MET_phi_jesRelativePtBBUp, &b_MET_phi_jesRelativePtBBUp);
   T1->SetBranchAddress("Jet_pt_jesRelativePtEC1Up", Jet_pt_jesRelativePtEC1Up, &b_Jet_pt_jesRelativePtEC1Up);
   T1->SetBranchAddress("Jet_mass_jesRelativePtEC1Up", Jet_mass_jesRelativePtEC1Up, &b_Jet_mass_jesRelativePtEC1Up);
   T1->SetBranchAddress("MET_pt_jesRelativePtEC1Up", &MET_pt_jesRelativePtEC1Up, &b_MET_pt_jesRelativePtEC1Up);
   T1->SetBranchAddress("MET_phi_jesRelativePtEC1Up", &MET_phi_jesRelativePtEC1Up, &b_MET_phi_jesRelativePtEC1Up);
   T1->SetBranchAddress("Jet_pt_jesRelativePtEC2Up", Jet_pt_jesRelativePtEC2Up, &b_Jet_pt_jesRelativePtEC2Up);
   T1->SetBranchAddress("Jet_mass_jesRelativePtEC2Up", Jet_mass_jesRelativePtEC2Up, &b_Jet_mass_jesRelativePtEC2Up);
   T1->SetBranchAddress("MET_pt_jesRelativePtEC2Up", &MET_pt_jesRelativePtEC2Up, &b_MET_pt_jesRelativePtEC2Up);
   T1->SetBranchAddress("MET_phi_jesRelativePtEC2Up", &MET_phi_jesRelativePtEC2Up, &b_MET_phi_jesRelativePtEC2Up);
   T1->SetBranchAddress("Jet_pt_jesRelativePtHFUp", Jet_pt_jesRelativePtHFUp, &b_Jet_pt_jesRelativePtHFUp);
   T1->SetBranchAddress("Jet_mass_jesRelativePtHFUp", Jet_mass_jesRelativePtHFUp, &b_Jet_mass_jesRelativePtHFUp);
   T1->SetBranchAddress("MET_pt_jesRelativePtHFUp", &MET_pt_jesRelativePtHFUp, &b_MET_pt_jesRelativePtHFUp);
   T1->SetBranchAddress("MET_phi_jesRelativePtHFUp", &MET_phi_jesRelativePtHFUp, &b_MET_phi_jesRelativePtHFUp);
   T1->SetBranchAddress("Jet_pt_jesRelativeBalUp", Jet_pt_jesRelativeBalUp, &b_Jet_pt_jesRelativeBalUp);
   T1->SetBranchAddress("Jet_mass_jesRelativeBalUp", Jet_mass_jesRelativeBalUp, &b_Jet_mass_jesRelativeBalUp);
   T1->SetBranchAddress("MET_pt_jesRelativeBalUp", &MET_pt_jesRelativeBalUp, &b_MET_pt_jesRelativeBalUp);
   T1->SetBranchAddress("MET_phi_jesRelativeBalUp", &MET_phi_jesRelativeBalUp, &b_MET_phi_jesRelativeBalUp);
   T1->SetBranchAddress("Jet_pt_jesRelativeSampleUp", Jet_pt_jesRelativeSampleUp, &b_Jet_pt_jesRelativeSampleUp);
   T1->SetBranchAddress("Jet_mass_jesRelativeSampleUp", Jet_mass_jesRelativeSampleUp, &b_Jet_mass_jesRelativeSampleUp);
   T1->SetBranchAddress("MET_pt_jesRelativeSampleUp", &MET_pt_jesRelativeSampleUp, &b_MET_pt_jesRelativeSampleUp);
   T1->SetBranchAddress("MET_phi_jesRelativeSampleUp", &MET_phi_jesRelativeSampleUp, &b_MET_phi_jesRelativeSampleUp);
   T1->SetBranchAddress("Jet_pt_jesRelativeFSRUp", Jet_pt_jesRelativeFSRUp, &b_Jet_pt_jesRelativeFSRUp);
   T1->SetBranchAddress("Jet_mass_jesRelativeFSRUp", Jet_mass_jesRelativeFSRUp, &b_Jet_mass_jesRelativeFSRUp);
   T1->SetBranchAddress("MET_pt_jesRelativeFSRUp", &MET_pt_jesRelativeFSRUp, &b_MET_pt_jesRelativeFSRUp);
   T1->SetBranchAddress("MET_phi_jesRelativeFSRUp", &MET_phi_jesRelativeFSRUp, &b_MET_phi_jesRelativeFSRUp);
   T1->SetBranchAddress("Jet_pt_jesRelativeStatFSRUp", Jet_pt_jesRelativeStatFSRUp, &b_Jet_pt_jesRelativeStatFSRUp);
   T1->SetBranchAddress("Jet_mass_jesRelativeStatFSRUp", Jet_mass_jesRelativeStatFSRUp, &b_Jet_mass_jesRelativeStatFSRUp);
   T1->SetBranchAddress("MET_pt_jesRelativeStatFSRUp", &MET_pt_jesRelativeStatFSRUp, &b_MET_pt_jesRelativeStatFSRUp);
   T1->SetBranchAddress("MET_phi_jesRelativeStatFSRUp", &MET_phi_jesRelativeStatFSRUp, &b_MET_phi_jesRelativeStatFSRUp);
   T1->SetBranchAddress("Jet_pt_jesRelativeStatECUp", Jet_pt_jesRelativeStatECUp, &b_Jet_pt_jesRelativeStatECUp);
   T1->SetBranchAddress("Jet_mass_jesRelativeStatECUp", Jet_mass_jesRelativeStatECUp, &b_Jet_mass_jesRelativeStatECUp);
   T1->SetBranchAddress("MET_pt_jesRelativeStatECUp", &MET_pt_jesRelativeStatECUp, &b_MET_pt_jesRelativeStatECUp);
   T1->SetBranchAddress("MET_phi_jesRelativeStatECUp", &MET_phi_jesRelativeStatECUp, &b_MET_phi_jesRelativeStatECUp);
   T1->SetBranchAddress("Jet_pt_jesRelativeStatHFUp", Jet_pt_jesRelativeStatHFUp, &b_Jet_pt_jesRelativeStatHFUp);
   T1->SetBranchAddress("Jet_mass_jesRelativeStatHFUp", Jet_mass_jesRelativeStatHFUp, &b_Jet_mass_jesRelativeStatHFUp);
   T1->SetBranchAddress("MET_pt_jesRelativeStatHFUp", &MET_pt_jesRelativeStatHFUp, &b_MET_pt_jesRelativeStatHFUp);
   T1->SetBranchAddress("MET_phi_jesRelativeStatHFUp", &MET_phi_jesRelativeStatHFUp, &b_MET_phi_jesRelativeStatHFUp);
   T1->SetBranchAddress("Jet_pt_jesPileUpDataMCUp", Jet_pt_jesPileUpDataMCUp, &b_Jet_pt_jesPileUpDataMCUp);
   T1->SetBranchAddress("Jet_mass_jesPileUpDataMCUp", Jet_mass_jesPileUpDataMCUp, &b_Jet_mass_jesPileUpDataMCUp);
   T1->SetBranchAddress("MET_pt_jesPileUpDataMCUp", &MET_pt_jesPileUpDataMCUp, &b_MET_pt_jesPileUpDataMCUp);
   T1->SetBranchAddress("MET_phi_jesPileUpDataMCUp", &MET_phi_jesPileUpDataMCUp, &b_MET_phi_jesPileUpDataMCUp);
   T1->SetBranchAddress("Jet_pt_jesPileUpPtRefUp", Jet_pt_jesPileUpPtRefUp, &b_Jet_pt_jesPileUpPtRefUp);
   T1->SetBranchAddress("Jet_mass_jesPileUpPtRefUp", Jet_mass_jesPileUpPtRefUp, &b_Jet_mass_jesPileUpPtRefUp);
   T1->SetBranchAddress("MET_pt_jesPileUpPtRefUp", &MET_pt_jesPileUpPtRefUp, &b_MET_pt_jesPileUpPtRefUp);
   T1->SetBranchAddress("MET_phi_jesPileUpPtRefUp", &MET_phi_jesPileUpPtRefUp, &b_MET_phi_jesPileUpPtRefUp);
   T1->SetBranchAddress("Jet_pt_jesPileUpPtBBUp", Jet_pt_jesPileUpPtBBUp, &b_Jet_pt_jesPileUpPtBBUp);
   T1->SetBranchAddress("Jet_mass_jesPileUpPtBBUp", Jet_mass_jesPileUpPtBBUp, &b_Jet_mass_jesPileUpPtBBUp);
   T1->SetBranchAddress("MET_pt_jesPileUpPtBBUp", &MET_pt_jesPileUpPtBBUp, &b_MET_pt_jesPileUpPtBBUp);
   T1->SetBranchAddress("MET_phi_jesPileUpPtBBUp", &MET_phi_jesPileUpPtBBUp, &b_MET_phi_jesPileUpPtBBUp);
   T1->SetBranchAddress("Jet_pt_jesPileUpPtEC1Up", Jet_pt_jesPileUpPtEC1Up, &b_Jet_pt_jesPileUpPtEC1Up);
   T1->SetBranchAddress("Jet_mass_jesPileUpPtEC1Up", Jet_mass_jesPileUpPtEC1Up, &b_Jet_mass_jesPileUpPtEC1Up);
   T1->SetBranchAddress("MET_pt_jesPileUpPtEC1Up", &MET_pt_jesPileUpPtEC1Up, &b_MET_pt_jesPileUpPtEC1Up);
   T1->SetBranchAddress("MET_phi_jesPileUpPtEC1Up", &MET_phi_jesPileUpPtEC1Up, &b_MET_phi_jesPileUpPtEC1Up);
   T1->SetBranchAddress("Jet_pt_jesPileUpPtEC2Up", Jet_pt_jesPileUpPtEC2Up, &b_Jet_pt_jesPileUpPtEC2Up);
   T1->SetBranchAddress("Jet_mass_jesPileUpPtEC2Up", Jet_mass_jesPileUpPtEC2Up, &b_Jet_mass_jesPileUpPtEC2Up);
   T1->SetBranchAddress("MET_pt_jesPileUpPtEC2Up", &MET_pt_jesPileUpPtEC2Up, &b_MET_pt_jesPileUpPtEC2Up);
   T1->SetBranchAddress("MET_phi_jesPileUpPtEC2Up", &MET_phi_jesPileUpPtEC2Up, &b_MET_phi_jesPileUpPtEC2Up);
   T1->SetBranchAddress("Jet_pt_jesPileUpPtHFUp", Jet_pt_jesPileUpPtHFUp, &b_Jet_pt_jesPileUpPtHFUp);
   T1->SetBranchAddress("Jet_mass_jesPileUpPtHFUp", Jet_mass_jesPileUpPtHFUp, &b_Jet_mass_jesPileUpPtHFUp);
   T1->SetBranchAddress("MET_pt_jesPileUpPtHFUp", &MET_pt_jesPileUpPtHFUp, &b_MET_pt_jesPileUpPtHFUp);
   T1->SetBranchAddress("MET_phi_jesPileUpPtHFUp", &MET_phi_jesPileUpPtHFUp, &b_MET_phi_jesPileUpPtHFUp);
   T1->SetBranchAddress("Jet_pt_jesPileUpMuZeroUp", Jet_pt_jesPileUpMuZeroUp, &b_Jet_pt_jesPileUpMuZeroUp);
   T1->SetBranchAddress("Jet_mass_jesPileUpMuZeroUp", Jet_mass_jesPileUpMuZeroUp, &b_Jet_mass_jesPileUpMuZeroUp);
   T1->SetBranchAddress("MET_pt_jesPileUpMuZeroUp", &MET_pt_jesPileUpMuZeroUp, &b_MET_pt_jesPileUpMuZeroUp);
   T1->SetBranchAddress("MET_phi_jesPileUpMuZeroUp", &MET_phi_jesPileUpMuZeroUp, &b_MET_phi_jesPileUpMuZeroUp);
   T1->SetBranchAddress("Jet_pt_jesPileUpEnvelopeUp", Jet_pt_jesPileUpEnvelopeUp, &b_Jet_pt_jesPileUpEnvelopeUp);
   T1->SetBranchAddress("Jet_mass_jesPileUpEnvelopeUp", Jet_mass_jesPileUpEnvelopeUp, &b_Jet_mass_jesPileUpEnvelopeUp);
   T1->SetBranchAddress("MET_pt_jesPileUpEnvelopeUp", &MET_pt_jesPileUpEnvelopeUp, &b_MET_pt_jesPileUpEnvelopeUp);
   T1->SetBranchAddress("MET_phi_jesPileUpEnvelopeUp", &MET_phi_jesPileUpEnvelopeUp, &b_MET_phi_jesPileUpEnvelopeUp);
   T1->SetBranchAddress("Jet_pt_jesSubTotalPileUpUp", Jet_pt_jesSubTotalPileUpUp, &b_Jet_pt_jesSubTotalPileUpUp);
   T1->SetBranchAddress("Jet_mass_jesSubTotalPileUpUp", Jet_mass_jesSubTotalPileUpUp, &b_Jet_mass_jesSubTotalPileUpUp);
   T1->SetBranchAddress("MET_pt_jesSubTotalPileUpUp", &MET_pt_jesSubTotalPileUpUp, &b_MET_pt_jesSubTotalPileUpUp);
   T1->SetBranchAddress("MET_phi_jesSubTotalPileUpUp", &MET_phi_jesSubTotalPileUpUp, &b_MET_phi_jesSubTotalPileUpUp);
   T1->SetBranchAddress("Jet_pt_jesSubTotalRelativeUp", Jet_pt_jesSubTotalRelativeUp, &b_Jet_pt_jesSubTotalRelativeUp);
   T1->SetBranchAddress("Jet_mass_jesSubTotalRelativeUp", Jet_mass_jesSubTotalRelativeUp, &b_Jet_mass_jesSubTotalRelativeUp);
   T1->SetBranchAddress("MET_pt_jesSubTotalRelativeUp", &MET_pt_jesSubTotalRelativeUp, &b_MET_pt_jesSubTotalRelativeUp);
   T1->SetBranchAddress("MET_phi_jesSubTotalRelativeUp", &MET_phi_jesSubTotalRelativeUp, &b_MET_phi_jesSubTotalRelativeUp);
   T1->SetBranchAddress("Jet_pt_jesSubTotalPtUp", Jet_pt_jesSubTotalPtUp, &b_Jet_pt_jesSubTotalPtUp);
   T1->SetBranchAddress("Jet_mass_jesSubTotalPtUp", Jet_mass_jesSubTotalPtUp, &b_Jet_mass_jesSubTotalPtUp);
   T1->SetBranchAddress("MET_pt_jesSubTotalPtUp", &MET_pt_jesSubTotalPtUp, &b_MET_pt_jesSubTotalPtUp);
   T1->SetBranchAddress("MET_phi_jesSubTotalPtUp", &MET_phi_jesSubTotalPtUp, &b_MET_phi_jesSubTotalPtUp);
   T1->SetBranchAddress("Jet_pt_jesSubTotalScaleUp", Jet_pt_jesSubTotalScaleUp, &b_Jet_pt_jesSubTotalScaleUp);
   T1->SetBranchAddress("Jet_mass_jesSubTotalScaleUp", Jet_mass_jesSubTotalScaleUp, &b_Jet_mass_jesSubTotalScaleUp);
   T1->SetBranchAddress("MET_pt_jesSubTotalScaleUp", &MET_pt_jesSubTotalScaleUp, &b_MET_pt_jesSubTotalScaleUp);
   T1->SetBranchAddress("MET_phi_jesSubTotalScaleUp", &MET_phi_jesSubTotalScaleUp, &b_MET_phi_jesSubTotalScaleUp);
   T1->SetBranchAddress("Jet_pt_jesSubTotalAbsoluteUp", Jet_pt_jesSubTotalAbsoluteUp, &b_Jet_pt_jesSubTotalAbsoluteUp);
   T1->SetBranchAddress("Jet_mass_jesSubTotalAbsoluteUp", Jet_mass_jesSubTotalAbsoluteUp, &b_Jet_mass_jesSubTotalAbsoluteUp);
   T1->SetBranchAddress("MET_pt_jesSubTotalAbsoluteUp", &MET_pt_jesSubTotalAbsoluteUp, &b_MET_pt_jesSubTotalAbsoluteUp);
   T1->SetBranchAddress("MET_phi_jesSubTotalAbsoluteUp", &MET_phi_jesSubTotalAbsoluteUp, &b_MET_phi_jesSubTotalAbsoluteUp);
   T1->SetBranchAddress("Jet_pt_jesSubTotalMCUp", Jet_pt_jesSubTotalMCUp, &b_Jet_pt_jesSubTotalMCUp);
   T1->SetBranchAddress("Jet_mass_jesSubTotalMCUp", Jet_mass_jesSubTotalMCUp, &b_Jet_mass_jesSubTotalMCUp);
   T1->SetBranchAddress("MET_pt_jesSubTotalMCUp", &MET_pt_jesSubTotalMCUp, &b_MET_pt_jesSubTotalMCUp);
   T1->SetBranchAddress("MET_phi_jesSubTotalMCUp", &MET_phi_jesSubTotalMCUp, &b_MET_phi_jesSubTotalMCUp);
   T1->SetBranchAddress("Jet_pt_jesTotalUp", Jet_pt_jesTotalUp, &b_Jet_pt_jesTotalUp);
   T1->SetBranchAddress("Jet_mass_jesTotalUp", Jet_mass_jesTotalUp, &b_Jet_mass_jesTotalUp);
   T1->SetBranchAddress("MET_pt_jesTotalUp", &MET_pt_jesTotalUp, &b_MET_pt_jesTotalUp);
   T1->SetBranchAddress("MET_phi_jesTotalUp", &MET_phi_jesTotalUp, &b_MET_phi_jesTotalUp);
   T1->SetBranchAddress("Jet_pt_jesTotalNoFlavorUp", Jet_pt_jesTotalNoFlavorUp, &b_Jet_pt_jesTotalNoFlavorUp);
   T1->SetBranchAddress("Jet_mass_jesTotalNoFlavorUp", Jet_mass_jesTotalNoFlavorUp, &b_Jet_mass_jesTotalNoFlavorUp);
   T1->SetBranchAddress("MET_pt_jesTotalNoFlavorUp", &MET_pt_jesTotalNoFlavorUp, &b_MET_pt_jesTotalNoFlavorUp);
   T1->SetBranchAddress("MET_phi_jesTotalNoFlavorUp", &MET_phi_jesTotalNoFlavorUp, &b_MET_phi_jesTotalNoFlavorUp);
   T1->SetBranchAddress("Jet_pt_jesTotalNoTimeUp", Jet_pt_jesTotalNoTimeUp, &b_Jet_pt_jesTotalNoTimeUp);
   T1->SetBranchAddress("Jet_mass_jesTotalNoTimeUp", Jet_mass_jesTotalNoTimeUp, &b_Jet_mass_jesTotalNoTimeUp);
   T1->SetBranchAddress("MET_pt_jesTotalNoTimeUp", &MET_pt_jesTotalNoTimeUp, &b_MET_pt_jesTotalNoTimeUp);
   T1->SetBranchAddress("MET_phi_jesTotalNoTimeUp", &MET_phi_jesTotalNoTimeUp, &b_MET_phi_jesTotalNoTimeUp);
   T1->SetBranchAddress("Jet_pt_jesTotalNoFlavorNoTimeUp", Jet_pt_jesTotalNoFlavorNoTimeUp, &b_Jet_pt_jesTotalNoFlavorNoTimeUp);
   T1->SetBranchAddress("Jet_mass_jesTotalNoFlavorNoTimeUp", Jet_mass_jesTotalNoFlavorNoTimeUp, &b_Jet_mass_jesTotalNoFlavorNoTimeUp);
   T1->SetBranchAddress("MET_pt_jesTotalNoFlavorNoTimeUp", &MET_pt_jesTotalNoFlavorNoTimeUp, &b_MET_pt_jesTotalNoFlavorNoTimeUp);
   T1->SetBranchAddress("MET_phi_jesTotalNoFlavorNoTimeUp", &MET_phi_jesTotalNoFlavorNoTimeUp, &b_MET_phi_jesTotalNoFlavorNoTimeUp);
   T1->SetBranchAddress("Jet_pt_jesFlavorZJetUp", Jet_pt_jesFlavorZJetUp, &b_Jet_pt_jesFlavorZJetUp);
   T1->SetBranchAddress("Jet_mass_jesFlavorZJetUp", Jet_mass_jesFlavorZJetUp, &b_Jet_mass_jesFlavorZJetUp);
   T1->SetBranchAddress("MET_pt_jesFlavorZJetUp", &MET_pt_jesFlavorZJetUp, &b_MET_pt_jesFlavorZJetUp);
   T1->SetBranchAddress("MET_phi_jesFlavorZJetUp", &MET_phi_jesFlavorZJetUp, &b_MET_phi_jesFlavorZJetUp);
   T1->SetBranchAddress("Jet_pt_jesFlavorPhotonJetUp", Jet_pt_jesFlavorPhotonJetUp, &b_Jet_pt_jesFlavorPhotonJetUp);
   T1->SetBranchAddress("Jet_mass_jesFlavorPhotonJetUp", Jet_mass_jesFlavorPhotonJetUp, &b_Jet_mass_jesFlavorPhotonJetUp);
   T1->SetBranchAddress("MET_pt_jesFlavorPhotonJetUp", &MET_pt_jesFlavorPhotonJetUp, &b_MET_pt_jesFlavorPhotonJetUp);
   T1->SetBranchAddress("MET_phi_jesFlavorPhotonJetUp", &MET_phi_jesFlavorPhotonJetUp, &b_MET_phi_jesFlavorPhotonJetUp);
   T1->SetBranchAddress("Jet_pt_jesFlavorPureGluonUp", Jet_pt_jesFlavorPureGluonUp, &b_Jet_pt_jesFlavorPureGluonUp);
   T1->SetBranchAddress("Jet_mass_jesFlavorPureGluonUp", Jet_mass_jesFlavorPureGluonUp, &b_Jet_mass_jesFlavorPureGluonUp);
   T1->SetBranchAddress("MET_pt_jesFlavorPureGluonUp", &MET_pt_jesFlavorPureGluonUp, &b_MET_pt_jesFlavorPureGluonUp);
   T1->SetBranchAddress("MET_phi_jesFlavorPureGluonUp", &MET_phi_jesFlavorPureGluonUp, &b_MET_phi_jesFlavorPureGluonUp);
   T1->SetBranchAddress("Jet_pt_jesFlavorPureQuarkUp", Jet_pt_jesFlavorPureQuarkUp, &b_Jet_pt_jesFlavorPureQuarkUp);
   T1->SetBranchAddress("Jet_mass_jesFlavorPureQuarkUp", Jet_mass_jesFlavorPureQuarkUp, &b_Jet_mass_jesFlavorPureQuarkUp);
   T1->SetBranchAddress("MET_pt_jesFlavorPureQuarkUp", &MET_pt_jesFlavorPureQuarkUp, &b_MET_pt_jesFlavorPureQuarkUp);
   T1->SetBranchAddress("MET_phi_jesFlavorPureQuarkUp", &MET_phi_jesFlavorPureQuarkUp, &b_MET_phi_jesFlavorPureQuarkUp);
   T1->SetBranchAddress("Jet_pt_jesFlavorPureCharmUp", Jet_pt_jesFlavorPureCharmUp, &b_Jet_pt_jesFlavorPureCharmUp);
   T1->SetBranchAddress("Jet_mass_jesFlavorPureCharmUp", Jet_mass_jesFlavorPureCharmUp, &b_Jet_mass_jesFlavorPureCharmUp);
   T1->SetBranchAddress("MET_pt_jesFlavorPureCharmUp", &MET_pt_jesFlavorPureCharmUp, &b_MET_pt_jesFlavorPureCharmUp);
   T1->SetBranchAddress("MET_phi_jesFlavorPureCharmUp", &MET_phi_jesFlavorPureCharmUp, &b_MET_phi_jesFlavorPureCharmUp);
   T1->SetBranchAddress("Jet_pt_jesFlavorPureBottomUp", Jet_pt_jesFlavorPureBottomUp, &b_Jet_pt_jesFlavorPureBottomUp);
   T1->SetBranchAddress("Jet_mass_jesFlavorPureBottomUp", Jet_mass_jesFlavorPureBottomUp, &b_Jet_mass_jesFlavorPureBottomUp);
   T1->SetBranchAddress("MET_pt_jesFlavorPureBottomUp", &MET_pt_jesFlavorPureBottomUp, &b_MET_pt_jesFlavorPureBottomUp);
   T1->SetBranchAddress("MET_phi_jesFlavorPureBottomUp", &MET_phi_jesFlavorPureBottomUp, &b_MET_phi_jesFlavorPureBottomUp);
   T1->SetBranchAddress("Jet_pt_jesTimeRunBUp", Jet_pt_jesTimeRunBUp, &b_Jet_pt_jesTimeRunBUp);
   T1->SetBranchAddress("Jet_mass_jesTimeRunBUp", Jet_mass_jesTimeRunBUp, &b_Jet_mass_jesTimeRunBUp);
   T1->SetBranchAddress("MET_pt_jesTimeRunBUp", &MET_pt_jesTimeRunBUp, &b_MET_pt_jesTimeRunBUp);
   T1->SetBranchAddress("MET_phi_jesTimeRunBUp", &MET_phi_jesTimeRunBUp, &b_MET_phi_jesTimeRunBUp);
   T1->SetBranchAddress("Jet_pt_jesTimeRunCUp", Jet_pt_jesTimeRunCUp, &b_Jet_pt_jesTimeRunCUp);
   T1->SetBranchAddress("Jet_mass_jesTimeRunCUp", Jet_mass_jesTimeRunCUp, &b_Jet_mass_jesTimeRunCUp);
   T1->SetBranchAddress("MET_pt_jesTimeRunCUp", &MET_pt_jesTimeRunCUp, &b_MET_pt_jesTimeRunCUp);
   T1->SetBranchAddress("MET_phi_jesTimeRunCUp", &MET_phi_jesTimeRunCUp, &b_MET_phi_jesTimeRunCUp);
   T1->SetBranchAddress("Jet_pt_jesTimeRunDEUp", Jet_pt_jesTimeRunDEUp, &b_Jet_pt_jesTimeRunDEUp);
   T1->SetBranchAddress("Jet_mass_jesTimeRunDEUp", Jet_mass_jesTimeRunDEUp, &b_Jet_mass_jesTimeRunDEUp);
   T1->SetBranchAddress("MET_pt_jesTimeRunDEUp", &MET_pt_jesTimeRunDEUp, &b_MET_pt_jesTimeRunDEUp);
   T1->SetBranchAddress("MET_phi_jesTimeRunDEUp", &MET_phi_jesTimeRunDEUp, &b_MET_phi_jesTimeRunDEUp);
   T1->SetBranchAddress("Jet_pt_jesTimeRunFUp", Jet_pt_jesTimeRunFUp, &b_Jet_pt_jesTimeRunFUp);
   T1->SetBranchAddress("Jet_mass_jesTimeRunFUp", Jet_mass_jesTimeRunFUp, &b_Jet_mass_jesTimeRunFUp);
   T1->SetBranchAddress("MET_pt_jesTimeRunFUp", &MET_pt_jesTimeRunFUp, &b_MET_pt_jesTimeRunFUp);
   T1->SetBranchAddress("MET_phi_jesTimeRunFUp", &MET_phi_jesTimeRunFUp, &b_MET_phi_jesTimeRunFUp);
   T1->SetBranchAddress("Jet_pt_jesCorrelationGroupMPFInSituUp", Jet_pt_jesCorrelationGroupMPFInSituUp, &b_Jet_pt_jesCorrelationGroupMPFInSituUp);
   T1->SetBranchAddress("Jet_mass_jesCorrelationGroupMPFInSituUp", Jet_mass_jesCorrelationGroupMPFInSituUp, &b_Jet_mass_jesCorrelationGroupMPFInSituUp);
   T1->SetBranchAddress("MET_pt_jesCorrelationGroupMPFInSituUp", &MET_pt_jesCorrelationGroupMPFInSituUp, &b_MET_pt_jesCorrelationGroupMPFInSituUp);
   T1->SetBranchAddress("MET_phi_jesCorrelationGroupMPFInSituUp", &MET_phi_jesCorrelationGroupMPFInSituUp, &b_MET_phi_jesCorrelationGroupMPFInSituUp);
   T1->SetBranchAddress("Jet_pt_jesCorrelationGroupIntercalibrationUp", Jet_pt_jesCorrelationGroupIntercalibrationUp, &b_Jet_pt_jesCorrelationGroupIntercalibrationUp);
   T1->SetBranchAddress("Jet_mass_jesCorrelationGroupIntercalibrationUp", Jet_mass_jesCorrelationGroupIntercalibrationUp, &b_Jet_mass_jesCorrelationGroupIntercalibrationUp);
   T1->SetBranchAddress("MET_pt_jesCorrelationGroupIntercalibrationUp", &MET_pt_jesCorrelationGroupIntercalibrationUp, &b_MET_pt_jesCorrelationGroupIntercalibrationUp);
   T1->SetBranchAddress("MET_phi_jesCorrelationGroupIntercalibrationUp", &MET_phi_jesCorrelationGroupIntercalibrationUp, &b_MET_phi_jesCorrelationGroupIntercalibrationUp);
   T1->SetBranchAddress("Jet_pt_jesCorrelationGroupbJESUp", Jet_pt_jesCorrelationGroupbJESUp, &b_Jet_pt_jesCorrelationGroupbJESUp);
   T1->SetBranchAddress("Jet_mass_jesCorrelationGroupbJESUp", Jet_mass_jesCorrelationGroupbJESUp, &b_Jet_mass_jesCorrelationGroupbJESUp);
   T1->SetBranchAddress("MET_pt_jesCorrelationGroupbJESUp", &MET_pt_jesCorrelationGroupbJESUp, &b_MET_pt_jesCorrelationGroupbJESUp);
   T1->SetBranchAddress("MET_phi_jesCorrelationGroupbJESUp", &MET_phi_jesCorrelationGroupbJESUp, &b_MET_phi_jesCorrelationGroupbJESUp);
   T1->SetBranchAddress("Jet_pt_jesCorrelationGroupFlavorUp", Jet_pt_jesCorrelationGroupFlavorUp, &b_Jet_pt_jesCorrelationGroupFlavorUp);
   T1->SetBranchAddress("Jet_mass_jesCorrelationGroupFlavorUp", Jet_mass_jesCorrelationGroupFlavorUp, &b_Jet_mass_jesCorrelationGroupFlavorUp);
   T1->SetBranchAddress("MET_pt_jesCorrelationGroupFlavorUp", &MET_pt_jesCorrelationGroupFlavorUp, &b_MET_pt_jesCorrelationGroupFlavorUp);
   T1->SetBranchAddress("MET_phi_jesCorrelationGroupFlavorUp", &MET_phi_jesCorrelationGroupFlavorUp, &b_MET_phi_jesCorrelationGroupFlavorUp);
   T1->SetBranchAddress("Jet_pt_jesCorrelationGroupUncorrelatedUp", Jet_pt_jesCorrelationGroupUncorrelatedUp, &b_Jet_pt_jesCorrelationGroupUncorrelatedUp);
   T1->SetBranchAddress("Jet_mass_jesCorrelationGroupUncorrelatedUp", Jet_mass_jesCorrelationGroupUncorrelatedUp, &b_Jet_mass_jesCorrelationGroupUncorrelatedUp);
   T1->SetBranchAddress("MET_pt_jesCorrelationGroupUncorrelatedUp", &MET_pt_jesCorrelationGroupUncorrelatedUp, &b_MET_pt_jesCorrelationGroupUncorrelatedUp);
   T1->SetBranchAddress("MET_phi_jesCorrelationGroupUncorrelatedUp", &MET_phi_jesCorrelationGroupUncorrelatedUp, &b_MET_phi_jesCorrelationGroupUncorrelatedUp);
   T1->SetBranchAddress("MET_pt_unclustEnUp", &MET_pt_unclustEnUp, &b_MET_pt_unclustEnUp);
   T1->SetBranchAddress("MET_phi_unclustEnUp", &MET_phi_unclustEnUp, &b_MET_phi_unclustEnUp);
   T1->SetBranchAddress("Jet_pt_jerDown", Jet_pt_jerDown, &b_Jet_pt_jerDown);
   T1->SetBranchAddress("Jet_mass_jerDown", Jet_mass_jerDown, &b_Jet_mass_jerDown);
   T1->SetBranchAddress("MET_pt_jerDown", &MET_pt_jerDown, &b_MET_pt_jerDown);
   T1->SetBranchAddress("MET_phi_jerDown", &MET_phi_jerDown, &b_MET_phi_jerDown);
   T1->SetBranchAddress("Jet_pt_jesAbsoluteStatDown", Jet_pt_jesAbsoluteStatDown, &b_Jet_pt_jesAbsoluteStatDown);
   T1->SetBranchAddress("Jet_mass_jesAbsoluteStatDown", Jet_mass_jesAbsoluteStatDown, &b_Jet_mass_jesAbsoluteStatDown);
   T1->SetBranchAddress("MET_pt_jesAbsoluteStatDown", &MET_pt_jesAbsoluteStatDown, &b_MET_pt_jesAbsoluteStatDown);
   T1->SetBranchAddress("MET_phi_jesAbsoluteStatDown", &MET_phi_jesAbsoluteStatDown, &b_MET_phi_jesAbsoluteStatDown);
   T1->SetBranchAddress("Jet_pt_jesAbsoluteScaleDown", Jet_pt_jesAbsoluteScaleDown, &b_Jet_pt_jesAbsoluteScaleDown);
   T1->SetBranchAddress("Jet_mass_jesAbsoluteScaleDown", Jet_mass_jesAbsoluteScaleDown, &b_Jet_mass_jesAbsoluteScaleDown);
   T1->SetBranchAddress("MET_pt_jesAbsoluteScaleDown", &MET_pt_jesAbsoluteScaleDown, &b_MET_pt_jesAbsoluteScaleDown);
   T1->SetBranchAddress("MET_phi_jesAbsoluteScaleDown", &MET_phi_jesAbsoluteScaleDown, &b_MET_phi_jesAbsoluteScaleDown);
   T1->SetBranchAddress("Jet_pt_jesAbsoluteFlavMapDown", Jet_pt_jesAbsoluteFlavMapDown, &b_Jet_pt_jesAbsoluteFlavMapDown);
   T1->SetBranchAddress("Jet_mass_jesAbsoluteFlavMapDown", Jet_mass_jesAbsoluteFlavMapDown, &b_Jet_mass_jesAbsoluteFlavMapDown);
   T1->SetBranchAddress("MET_pt_jesAbsoluteFlavMapDown", &MET_pt_jesAbsoluteFlavMapDown, &b_MET_pt_jesAbsoluteFlavMapDown);
   T1->SetBranchAddress("MET_phi_jesAbsoluteFlavMapDown", &MET_phi_jesAbsoluteFlavMapDown, &b_MET_phi_jesAbsoluteFlavMapDown);
   T1->SetBranchAddress("Jet_pt_jesAbsoluteMPFBiasDown", Jet_pt_jesAbsoluteMPFBiasDown, &b_Jet_pt_jesAbsoluteMPFBiasDown);
   T1->SetBranchAddress("Jet_mass_jesAbsoluteMPFBiasDown", Jet_mass_jesAbsoluteMPFBiasDown, &b_Jet_mass_jesAbsoluteMPFBiasDown);
   T1->SetBranchAddress("MET_pt_jesAbsoluteMPFBiasDown", &MET_pt_jesAbsoluteMPFBiasDown, &b_MET_pt_jesAbsoluteMPFBiasDown);
   T1->SetBranchAddress("MET_phi_jesAbsoluteMPFBiasDown", &MET_phi_jesAbsoluteMPFBiasDown, &b_MET_phi_jesAbsoluteMPFBiasDown);
   T1->SetBranchAddress("Jet_pt_jesFragmentationDown", Jet_pt_jesFragmentationDown, &b_Jet_pt_jesFragmentationDown);
   T1->SetBranchAddress("Jet_mass_jesFragmentationDown", Jet_mass_jesFragmentationDown, &b_Jet_mass_jesFragmentationDown);
   T1->SetBranchAddress("MET_pt_jesFragmentationDown", &MET_pt_jesFragmentationDown, &b_MET_pt_jesFragmentationDown);
   T1->SetBranchAddress("MET_phi_jesFragmentationDown", &MET_phi_jesFragmentationDown, &b_MET_phi_jesFragmentationDown);
   T1->SetBranchAddress("Jet_pt_jesSinglePionECALDown", Jet_pt_jesSinglePionECALDown, &b_Jet_pt_jesSinglePionECALDown);
   T1->SetBranchAddress("Jet_mass_jesSinglePionECALDown", Jet_mass_jesSinglePionECALDown, &b_Jet_mass_jesSinglePionECALDown);
   T1->SetBranchAddress("MET_pt_jesSinglePionECALDown", &MET_pt_jesSinglePionECALDown, &b_MET_pt_jesSinglePionECALDown);
   T1->SetBranchAddress("MET_phi_jesSinglePionECALDown", &MET_phi_jesSinglePionECALDown, &b_MET_phi_jesSinglePionECALDown);
   T1->SetBranchAddress("Jet_pt_jesSinglePionHCALDown", Jet_pt_jesSinglePionHCALDown, &b_Jet_pt_jesSinglePionHCALDown);
   T1->SetBranchAddress("Jet_mass_jesSinglePionHCALDown", Jet_mass_jesSinglePionHCALDown, &b_Jet_mass_jesSinglePionHCALDown);
   T1->SetBranchAddress("MET_pt_jesSinglePionHCALDown", &MET_pt_jesSinglePionHCALDown, &b_MET_pt_jesSinglePionHCALDown);
   T1->SetBranchAddress("MET_phi_jesSinglePionHCALDown", &MET_phi_jesSinglePionHCALDown, &b_MET_phi_jesSinglePionHCALDown);
   T1->SetBranchAddress("Jet_pt_jesFlavorQCDDown", Jet_pt_jesFlavorQCDDown, &b_Jet_pt_jesFlavorQCDDown);
   T1->SetBranchAddress("Jet_mass_jesFlavorQCDDown", Jet_mass_jesFlavorQCDDown, &b_Jet_mass_jesFlavorQCDDown);
   T1->SetBranchAddress("MET_pt_jesFlavorQCDDown", &MET_pt_jesFlavorQCDDown, &b_MET_pt_jesFlavorQCDDown);
   T1->SetBranchAddress("MET_phi_jesFlavorQCDDown", &MET_phi_jesFlavorQCDDown, &b_MET_phi_jesFlavorQCDDown);
   T1->SetBranchAddress("Jet_pt_jesTimePtEtaDown", Jet_pt_jesTimePtEtaDown, &b_Jet_pt_jesTimePtEtaDown);
   T1->SetBranchAddress("Jet_mass_jesTimePtEtaDown", Jet_mass_jesTimePtEtaDown, &b_Jet_mass_jesTimePtEtaDown);
   T1->SetBranchAddress("MET_pt_jesTimePtEtaDown", &MET_pt_jesTimePtEtaDown, &b_MET_pt_jesTimePtEtaDown);
   T1->SetBranchAddress("MET_phi_jesTimePtEtaDown", &MET_phi_jesTimePtEtaDown, &b_MET_phi_jesTimePtEtaDown);
   T1->SetBranchAddress("Jet_pt_jesRelativeJEREC1Down", Jet_pt_jesRelativeJEREC1Down, &b_Jet_pt_jesRelativeJEREC1Down);
   T1->SetBranchAddress("Jet_mass_jesRelativeJEREC1Down", Jet_mass_jesRelativeJEREC1Down, &b_Jet_mass_jesRelativeJEREC1Down);
   T1->SetBranchAddress("MET_pt_jesRelativeJEREC1Down", &MET_pt_jesRelativeJEREC1Down, &b_MET_pt_jesRelativeJEREC1Down);
   T1->SetBranchAddress("MET_phi_jesRelativeJEREC1Down", &MET_phi_jesRelativeJEREC1Down, &b_MET_phi_jesRelativeJEREC1Down);
   T1->SetBranchAddress("Jet_pt_jesRelativeJEREC2Down", Jet_pt_jesRelativeJEREC2Down, &b_Jet_pt_jesRelativeJEREC2Down);
   T1->SetBranchAddress("Jet_mass_jesRelativeJEREC2Down", Jet_mass_jesRelativeJEREC2Down, &b_Jet_mass_jesRelativeJEREC2Down);
   T1->SetBranchAddress("MET_pt_jesRelativeJEREC2Down", &MET_pt_jesRelativeJEREC2Down, &b_MET_pt_jesRelativeJEREC2Down);
   T1->SetBranchAddress("MET_phi_jesRelativeJEREC2Down", &MET_phi_jesRelativeJEREC2Down, &b_MET_phi_jesRelativeJEREC2Down);
   T1->SetBranchAddress("Jet_pt_jesRelativeJERHFDown", Jet_pt_jesRelativeJERHFDown, &b_Jet_pt_jesRelativeJERHFDown);
   T1->SetBranchAddress("Jet_mass_jesRelativeJERHFDown", Jet_mass_jesRelativeJERHFDown, &b_Jet_mass_jesRelativeJERHFDown);
   T1->SetBranchAddress("MET_pt_jesRelativeJERHFDown", &MET_pt_jesRelativeJERHFDown, &b_MET_pt_jesRelativeJERHFDown);
   T1->SetBranchAddress("MET_phi_jesRelativeJERHFDown", &MET_phi_jesRelativeJERHFDown, &b_MET_phi_jesRelativeJERHFDown);
   T1->SetBranchAddress("Jet_pt_jesRelativePtBBDown", Jet_pt_jesRelativePtBBDown, &b_Jet_pt_jesRelativePtBBDown);
   T1->SetBranchAddress("Jet_mass_jesRelativePtBBDown", Jet_mass_jesRelativePtBBDown, &b_Jet_mass_jesRelativePtBBDown);
   T1->SetBranchAddress("MET_pt_jesRelativePtBBDown", &MET_pt_jesRelativePtBBDown, &b_MET_pt_jesRelativePtBBDown);
   T1->SetBranchAddress("MET_phi_jesRelativePtBBDown", &MET_phi_jesRelativePtBBDown, &b_MET_phi_jesRelativePtBBDown);
   T1->SetBranchAddress("Jet_pt_jesRelativePtEC1Down", Jet_pt_jesRelativePtEC1Down, &b_Jet_pt_jesRelativePtEC1Down);
   T1->SetBranchAddress("Jet_mass_jesRelativePtEC1Down", Jet_mass_jesRelativePtEC1Down, &b_Jet_mass_jesRelativePtEC1Down);
   T1->SetBranchAddress("MET_pt_jesRelativePtEC1Down", &MET_pt_jesRelativePtEC1Down, &b_MET_pt_jesRelativePtEC1Down);
   T1->SetBranchAddress("MET_phi_jesRelativePtEC1Down", &MET_phi_jesRelativePtEC1Down, &b_MET_phi_jesRelativePtEC1Down);
   T1->SetBranchAddress("Jet_pt_jesRelativePtEC2Down", Jet_pt_jesRelativePtEC2Down, &b_Jet_pt_jesRelativePtEC2Down);
   T1->SetBranchAddress("Jet_mass_jesRelativePtEC2Down", Jet_mass_jesRelativePtEC2Down, &b_Jet_mass_jesRelativePtEC2Down);
   T1->SetBranchAddress("MET_pt_jesRelativePtEC2Down", &MET_pt_jesRelativePtEC2Down, &b_MET_pt_jesRelativePtEC2Down);
   T1->SetBranchAddress("MET_phi_jesRelativePtEC2Down", &MET_phi_jesRelativePtEC2Down, &b_MET_phi_jesRelativePtEC2Down);
   T1->SetBranchAddress("Jet_pt_jesRelativePtHFDown", Jet_pt_jesRelativePtHFDown, &b_Jet_pt_jesRelativePtHFDown);
   T1->SetBranchAddress("Jet_mass_jesRelativePtHFDown", Jet_mass_jesRelativePtHFDown, &b_Jet_mass_jesRelativePtHFDown);
   T1->SetBranchAddress("MET_pt_jesRelativePtHFDown", &MET_pt_jesRelativePtHFDown, &b_MET_pt_jesRelativePtHFDown);
   T1->SetBranchAddress("MET_phi_jesRelativePtHFDown", &MET_phi_jesRelativePtHFDown, &b_MET_phi_jesRelativePtHFDown);
   T1->SetBranchAddress("Jet_pt_jesRelativeBalDown", Jet_pt_jesRelativeBalDown, &b_Jet_pt_jesRelativeBalDown);
   T1->SetBranchAddress("Jet_mass_jesRelativeBalDown", Jet_mass_jesRelativeBalDown, &b_Jet_mass_jesRelativeBalDown);
   T1->SetBranchAddress("MET_pt_jesRelativeBalDown", &MET_pt_jesRelativeBalDown, &b_MET_pt_jesRelativeBalDown);
   T1->SetBranchAddress("MET_phi_jesRelativeBalDown", &MET_phi_jesRelativeBalDown, &b_MET_phi_jesRelativeBalDown);
   T1->SetBranchAddress("Jet_pt_jesRelativeSampleDown", Jet_pt_jesRelativeSampleDown, &b_Jet_pt_jesRelativeSampleDown);
   T1->SetBranchAddress("Jet_mass_jesRelativeSampleDown", Jet_mass_jesRelativeSampleDown, &b_Jet_mass_jesRelativeSampleDown);
   T1->SetBranchAddress("MET_pt_jesRelativeSampleDown", &MET_pt_jesRelativeSampleDown, &b_MET_pt_jesRelativeSampleDown);
   T1->SetBranchAddress("MET_phi_jesRelativeSampleDown", &MET_phi_jesRelativeSampleDown, &b_MET_phi_jesRelativeSampleDown);
   T1->SetBranchAddress("Jet_pt_jesRelativeFSRDown", Jet_pt_jesRelativeFSRDown, &b_Jet_pt_jesRelativeFSRDown);
   T1->SetBranchAddress("Jet_mass_jesRelativeFSRDown", Jet_mass_jesRelativeFSRDown, &b_Jet_mass_jesRelativeFSRDown);
   T1->SetBranchAddress("MET_pt_jesRelativeFSRDown", &MET_pt_jesRelativeFSRDown, &b_MET_pt_jesRelativeFSRDown);
   T1->SetBranchAddress("MET_phi_jesRelativeFSRDown", &MET_phi_jesRelativeFSRDown, &b_MET_phi_jesRelativeFSRDown);
   T1->SetBranchAddress("Jet_pt_jesRelativeStatFSRDown", Jet_pt_jesRelativeStatFSRDown, &b_Jet_pt_jesRelativeStatFSRDown);
   T1->SetBranchAddress("Jet_mass_jesRelativeStatFSRDown", Jet_mass_jesRelativeStatFSRDown, &b_Jet_mass_jesRelativeStatFSRDown);
   T1->SetBranchAddress("MET_pt_jesRelativeStatFSRDown", &MET_pt_jesRelativeStatFSRDown, &b_MET_pt_jesRelativeStatFSRDown);
   T1->SetBranchAddress("MET_phi_jesRelativeStatFSRDown", &MET_phi_jesRelativeStatFSRDown, &b_MET_phi_jesRelativeStatFSRDown);
   T1->SetBranchAddress("Jet_pt_jesRelativeStatECDown", Jet_pt_jesRelativeStatECDown, &b_Jet_pt_jesRelativeStatECDown);
   T1->SetBranchAddress("Jet_mass_jesRelativeStatECDown", Jet_mass_jesRelativeStatECDown, &b_Jet_mass_jesRelativeStatECDown);
   T1->SetBranchAddress("MET_pt_jesRelativeStatECDown", &MET_pt_jesRelativeStatECDown, &b_MET_pt_jesRelativeStatECDown);
   T1->SetBranchAddress("MET_phi_jesRelativeStatECDown", &MET_phi_jesRelativeStatECDown, &b_MET_phi_jesRelativeStatECDown);
   T1->SetBranchAddress("Jet_pt_jesRelativeStatHFDown", Jet_pt_jesRelativeStatHFDown, &b_Jet_pt_jesRelativeStatHFDown);
   T1->SetBranchAddress("Jet_mass_jesRelativeStatHFDown", Jet_mass_jesRelativeStatHFDown, &b_Jet_mass_jesRelativeStatHFDown);
   T1->SetBranchAddress("MET_pt_jesRelativeStatHFDown", &MET_pt_jesRelativeStatHFDown, &b_MET_pt_jesRelativeStatHFDown);
   T1->SetBranchAddress("MET_phi_jesRelativeStatHFDown", &MET_phi_jesRelativeStatHFDown, &b_MET_phi_jesRelativeStatHFDown);
   T1->SetBranchAddress("Jet_pt_jesPileUpDataMCDown", Jet_pt_jesPileUpDataMCDown, &b_Jet_pt_jesPileUpDataMCDown);
   T1->SetBranchAddress("Jet_mass_jesPileUpDataMCDown", Jet_mass_jesPileUpDataMCDown, &b_Jet_mass_jesPileUpDataMCDown);
   T1->SetBranchAddress("MET_pt_jesPileUpDataMCDown", &MET_pt_jesPileUpDataMCDown, &b_MET_pt_jesPileUpDataMCDown);
   T1->SetBranchAddress("MET_phi_jesPileUpDataMCDown", &MET_phi_jesPileUpDataMCDown, &b_MET_phi_jesPileUpDataMCDown);
   T1->SetBranchAddress("Jet_pt_jesPileUpPtRefDown", Jet_pt_jesPileUpPtRefDown, &b_Jet_pt_jesPileUpPtRefDown);
   T1->SetBranchAddress("Jet_mass_jesPileUpPtRefDown", Jet_mass_jesPileUpPtRefDown, &b_Jet_mass_jesPileUpPtRefDown);
   T1->SetBranchAddress("MET_pt_jesPileUpPtRefDown", &MET_pt_jesPileUpPtRefDown, &b_MET_pt_jesPileUpPtRefDown);
   T1->SetBranchAddress("MET_phi_jesPileUpPtRefDown", &MET_phi_jesPileUpPtRefDown, &b_MET_phi_jesPileUpPtRefDown);
   T1->SetBranchAddress("Jet_pt_jesPileUpPtBBDown", Jet_pt_jesPileUpPtBBDown, &b_Jet_pt_jesPileUpPtBBDown);
   T1->SetBranchAddress("Jet_mass_jesPileUpPtBBDown", Jet_mass_jesPileUpPtBBDown, &b_Jet_mass_jesPileUpPtBBDown);
   T1->SetBranchAddress("MET_pt_jesPileUpPtBBDown", &MET_pt_jesPileUpPtBBDown, &b_MET_pt_jesPileUpPtBBDown);
   T1->SetBranchAddress("MET_phi_jesPileUpPtBBDown", &MET_phi_jesPileUpPtBBDown, &b_MET_phi_jesPileUpPtBBDown);
   T1->SetBranchAddress("Jet_pt_jesPileUpPtEC1Down", Jet_pt_jesPileUpPtEC1Down, &b_Jet_pt_jesPileUpPtEC1Down);
   T1->SetBranchAddress("Jet_mass_jesPileUpPtEC1Down", Jet_mass_jesPileUpPtEC1Down, &b_Jet_mass_jesPileUpPtEC1Down);
   T1->SetBranchAddress("MET_pt_jesPileUpPtEC1Down", &MET_pt_jesPileUpPtEC1Down, &b_MET_pt_jesPileUpPtEC1Down);
   T1->SetBranchAddress("MET_phi_jesPileUpPtEC1Down", &MET_phi_jesPileUpPtEC1Down, &b_MET_phi_jesPileUpPtEC1Down);
   T1->SetBranchAddress("Jet_pt_jesPileUpPtEC2Down", Jet_pt_jesPileUpPtEC2Down, &b_Jet_pt_jesPileUpPtEC2Down);
   T1->SetBranchAddress("Jet_mass_jesPileUpPtEC2Down", Jet_mass_jesPileUpPtEC2Down, &b_Jet_mass_jesPileUpPtEC2Down);
   T1->SetBranchAddress("MET_pt_jesPileUpPtEC2Down", &MET_pt_jesPileUpPtEC2Down, &b_MET_pt_jesPileUpPtEC2Down);
   T1->SetBranchAddress("MET_phi_jesPileUpPtEC2Down", &MET_phi_jesPileUpPtEC2Down, &b_MET_phi_jesPileUpPtEC2Down);
   T1->SetBranchAddress("Jet_pt_jesPileUpPtHFDown", Jet_pt_jesPileUpPtHFDown, &b_Jet_pt_jesPileUpPtHFDown);
   T1->SetBranchAddress("Jet_mass_jesPileUpPtHFDown", Jet_mass_jesPileUpPtHFDown, &b_Jet_mass_jesPileUpPtHFDown);
   T1->SetBranchAddress("MET_pt_jesPileUpPtHFDown", &MET_pt_jesPileUpPtHFDown, &b_MET_pt_jesPileUpPtHFDown);
   T1->SetBranchAddress("MET_phi_jesPileUpPtHFDown", &MET_phi_jesPileUpPtHFDown, &b_MET_phi_jesPileUpPtHFDown);
   T1->SetBranchAddress("Jet_pt_jesPileUpMuZeroDown", Jet_pt_jesPileUpMuZeroDown, &b_Jet_pt_jesPileUpMuZeroDown);
   T1->SetBranchAddress("Jet_mass_jesPileUpMuZeroDown", Jet_mass_jesPileUpMuZeroDown, &b_Jet_mass_jesPileUpMuZeroDown);
   T1->SetBranchAddress("MET_pt_jesPileUpMuZeroDown", &MET_pt_jesPileUpMuZeroDown, &b_MET_pt_jesPileUpMuZeroDown);
   T1->SetBranchAddress("MET_phi_jesPileUpMuZeroDown", &MET_phi_jesPileUpMuZeroDown, &b_MET_phi_jesPileUpMuZeroDown);
   T1->SetBranchAddress("Jet_pt_jesPileUpEnvelopeDown", Jet_pt_jesPileUpEnvelopeDown, &b_Jet_pt_jesPileUpEnvelopeDown);
   T1->SetBranchAddress("Jet_mass_jesPileUpEnvelopeDown", Jet_mass_jesPileUpEnvelopeDown, &b_Jet_mass_jesPileUpEnvelopeDown);
   T1->SetBranchAddress("MET_pt_jesPileUpEnvelopeDown", &MET_pt_jesPileUpEnvelopeDown, &b_MET_pt_jesPileUpEnvelopeDown);
   T1->SetBranchAddress("MET_phi_jesPileUpEnvelopeDown", &MET_phi_jesPileUpEnvelopeDown, &b_MET_phi_jesPileUpEnvelopeDown);
   T1->SetBranchAddress("Jet_pt_jesSubTotalPileUpDown", Jet_pt_jesSubTotalPileUpDown, &b_Jet_pt_jesSubTotalPileUpDown);
   T1->SetBranchAddress("Jet_mass_jesSubTotalPileUpDown", Jet_mass_jesSubTotalPileUpDown, &b_Jet_mass_jesSubTotalPileUpDown);
   T1->SetBranchAddress("MET_pt_jesSubTotalPileUpDown", &MET_pt_jesSubTotalPileUpDown, &b_MET_pt_jesSubTotalPileUpDown);
   T1->SetBranchAddress("MET_phi_jesSubTotalPileUpDown", &MET_phi_jesSubTotalPileUpDown, &b_MET_phi_jesSubTotalPileUpDown);
   T1->SetBranchAddress("Jet_pt_jesSubTotalRelativeDown", Jet_pt_jesSubTotalRelativeDown, &b_Jet_pt_jesSubTotalRelativeDown);
   T1->SetBranchAddress("Jet_mass_jesSubTotalRelativeDown", Jet_mass_jesSubTotalRelativeDown, &b_Jet_mass_jesSubTotalRelativeDown);
   T1->SetBranchAddress("MET_pt_jesSubTotalRelativeDown", &MET_pt_jesSubTotalRelativeDown, &b_MET_pt_jesSubTotalRelativeDown);
   T1->SetBranchAddress("MET_phi_jesSubTotalRelativeDown", &MET_phi_jesSubTotalRelativeDown, &b_MET_phi_jesSubTotalRelativeDown);
   T1->SetBranchAddress("Jet_pt_jesSubTotalPtDown", Jet_pt_jesSubTotalPtDown, &b_Jet_pt_jesSubTotalPtDown);
   T1->SetBranchAddress("Jet_mass_jesSubTotalPtDown", Jet_mass_jesSubTotalPtDown, &b_Jet_mass_jesSubTotalPtDown);
   T1->SetBranchAddress("MET_pt_jesSubTotalPtDown", &MET_pt_jesSubTotalPtDown, &b_MET_pt_jesSubTotalPtDown);
   T1->SetBranchAddress("MET_phi_jesSubTotalPtDown", &MET_phi_jesSubTotalPtDown, &b_MET_phi_jesSubTotalPtDown);
   T1->SetBranchAddress("Jet_pt_jesSubTotalScaleDown", Jet_pt_jesSubTotalScaleDown, &b_Jet_pt_jesSubTotalScaleDown);
   T1->SetBranchAddress("Jet_mass_jesSubTotalScaleDown", Jet_mass_jesSubTotalScaleDown, &b_Jet_mass_jesSubTotalScaleDown);
   T1->SetBranchAddress("MET_pt_jesSubTotalScaleDown", &MET_pt_jesSubTotalScaleDown, &b_MET_pt_jesSubTotalScaleDown);
   T1->SetBranchAddress("MET_phi_jesSubTotalScaleDown", &MET_phi_jesSubTotalScaleDown, &b_MET_phi_jesSubTotalScaleDown);
   T1->SetBranchAddress("Jet_pt_jesSubTotalAbsoluteDown", Jet_pt_jesSubTotalAbsoluteDown, &b_Jet_pt_jesSubTotalAbsoluteDown);
   T1->SetBranchAddress("Jet_mass_jesSubTotalAbsoluteDown", Jet_mass_jesSubTotalAbsoluteDown, &b_Jet_mass_jesSubTotalAbsoluteDown);
   T1->SetBranchAddress("MET_pt_jesSubTotalAbsoluteDown", &MET_pt_jesSubTotalAbsoluteDown, &b_MET_pt_jesSubTotalAbsoluteDown);
   T1->SetBranchAddress("MET_phi_jesSubTotalAbsoluteDown", &MET_phi_jesSubTotalAbsoluteDown, &b_MET_phi_jesSubTotalAbsoluteDown);
   T1->SetBranchAddress("Jet_pt_jesSubTotalMCDown", Jet_pt_jesSubTotalMCDown, &b_Jet_pt_jesSubTotalMCDown);
   T1->SetBranchAddress("Jet_mass_jesSubTotalMCDown", Jet_mass_jesSubTotalMCDown, &b_Jet_mass_jesSubTotalMCDown);
   T1->SetBranchAddress("MET_pt_jesSubTotalMCDown", &MET_pt_jesSubTotalMCDown, &b_MET_pt_jesSubTotalMCDown);
   T1->SetBranchAddress("MET_phi_jesSubTotalMCDown", &MET_phi_jesSubTotalMCDown, &b_MET_phi_jesSubTotalMCDown);
   T1->SetBranchAddress("Jet_pt_jesTotalDown", Jet_pt_jesTotalDown, &b_Jet_pt_jesTotalDown);
   T1->SetBranchAddress("Jet_mass_jesTotalDown", Jet_mass_jesTotalDown, &b_Jet_mass_jesTotalDown);
   T1->SetBranchAddress("MET_pt_jesTotalDown", &MET_pt_jesTotalDown, &b_MET_pt_jesTotalDown);
   T1->SetBranchAddress("MET_phi_jesTotalDown", &MET_phi_jesTotalDown, &b_MET_phi_jesTotalDown);
   T1->SetBranchAddress("Jet_pt_jesTotalNoFlavorDown", Jet_pt_jesTotalNoFlavorDown, &b_Jet_pt_jesTotalNoFlavorDown);
   T1->SetBranchAddress("Jet_mass_jesTotalNoFlavorDown", Jet_mass_jesTotalNoFlavorDown, &b_Jet_mass_jesTotalNoFlavorDown);
   T1->SetBranchAddress("MET_pt_jesTotalNoFlavorDown", &MET_pt_jesTotalNoFlavorDown, &b_MET_pt_jesTotalNoFlavorDown);
   T1->SetBranchAddress("MET_phi_jesTotalNoFlavorDown", &MET_phi_jesTotalNoFlavorDown, &b_MET_phi_jesTotalNoFlavorDown);
   T1->SetBranchAddress("Jet_pt_jesTotalNoTimeDown", Jet_pt_jesTotalNoTimeDown, &b_Jet_pt_jesTotalNoTimeDown);
   T1->SetBranchAddress("Jet_mass_jesTotalNoTimeDown", Jet_mass_jesTotalNoTimeDown, &b_Jet_mass_jesTotalNoTimeDown);
   T1->SetBranchAddress("MET_pt_jesTotalNoTimeDown", &MET_pt_jesTotalNoTimeDown, &b_MET_pt_jesTotalNoTimeDown);
   T1->SetBranchAddress("MET_phi_jesTotalNoTimeDown", &MET_phi_jesTotalNoTimeDown, &b_MET_phi_jesTotalNoTimeDown);
   T1->SetBranchAddress("Jet_pt_jesTotalNoFlavorNoTimeDown", Jet_pt_jesTotalNoFlavorNoTimeDown, &b_Jet_pt_jesTotalNoFlavorNoTimeDown);
   T1->SetBranchAddress("Jet_mass_jesTotalNoFlavorNoTimeDown", Jet_mass_jesTotalNoFlavorNoTimeDown, &b_Jet_mass_jesTotalNoFlavorNoTimeDown);
   T1->SetBranchAddress("MET_pt_jesTotalNoFlavorNoTimeDown", &MET_pt_jesTotalNoFlavorNoTimeDown, &b_MET_pt_jesTotalNoFlavorNoTimeDown);
   T1->SetBranchAddress("MET_phi_jesTotalNoFlavorNoTimeDown", &MET_phi_jesTotalNoFlavorNoTimeDown, &b_MET_phi_jesTotalNoFlavorNoTimeDown);
   T1->SetBranchAddress("Jet_pt_jesFlavorZJetDown", Jet_pt_jesFlavorZJetDown, &b_Jet_pt_jesFlavorZJetDown);
   T1->SetBranchAddress("Jet_mass_jesFlavorZJetDown", Jet_mass_jesFlavorZJetDown, &b_Jet_mass_jesFlavorZJetDown);
   T1->SetBranchAddress("MET_pt_jesFlavorZJetDown", &MET_pt_jesFlavorZJetDown, &b_MET_pt_jesFlavorZJetDown);
   T1->SetBranchAddress("MET_phi_jesFlavorZJetDown", &MET_phi_jesFlavorZJetDown, &b_MET_phi_jesFlavorZJetDown);
   T1->SetBranchAddress("Jet_pt_jesFlavorPhotonJetDown", Jet_pt_jesFlavorPhotonJetDown, &b_Jet_pt_jesFlavorPhotonJetDown);
   T1->SetBranchAddress("Jet_mass_jesFlavorPhotonJetDown", Jet_mass_jesFlavorPhotonJetDown, &b_Jet_mass_jesFlavorPhotonJetDown);
   T1->SetBranchAddress("MET_pt_jesFlavorPhotonJetDown", &MET_pt_jesFlavorPhotonJetDown, &b_MET_pt_jesFlavorPhotonJetDown);
   T1->SetBranchAddress("MET_phi_jesFlavorPhotonJetDown", &MET_phi_jesFlavorPhotonJetDown, &b_MET_phi_jesFlavorPhotonJetDown);
   T1->SetBranchAddress("Jet_pt_jesFlavorPureGluonDown", Jet_pt_jesFlavorPureGluonDown, &b_Jet_pt_jesFlavorPureGluonDown);
   T1->SetBranchAddress("Jet_mass_jesFlavorPureGluonDown", Jet_mass_jesFlavorPureGluonDown, &b_Jet_mass_jesFlavorPureGluonDown);
   T1->SetBranchAddress("MET_pt_jesFlavorPureGluonDown", &MET_pt_jesFlavorPureGluonDown, &b_MET_pt_jesFlavorPureGluonDown);
   T1->SetBranchAddress("MET_phi_jesFlavorPureGluonDown", &MET_phi_jesFlavorPureGluonDown, &b_MET_phi_jesFlavorPureGluonDown);
   T1->SetBranchAddress("Jet_pt_jesFlavorPureQuarkDown", Jet_pt_jesFlavorPureQuarkDown, &b_Jet_pt_jesFlavorPureQuarkDown);
   T1->SetBranchAddress("Jet_mass_jesFlavorPureQuarkDown", Jet_mass_jesFlavorPureQuarkDown, &b_Jet_mass_jesFlavorPureQuarkDown);
   T1->SetBranchAddress("MET_pt_jesFlavorPureQuarkDown", &MET_pt_jesFlavorPureQuarkDown, &b_MET_pt_jesFlavorPureQuarkDown);
   T1->SetBranchAddress("MET_phi_jesFlavorPureQuarkDown", &MET_phi_jesFlavorPureQuarkDown, &b_MET_phi_jesFlavorPureQuarkDown);
   T1->SetBranchAddress("Jet_pt_jesFlavorPureCharmDown", Jet_pt_jesFlavorPureCharmDown, &b_Jet_pt_jesFlavorPureCharmDown);
   T1->SetBranchAddress("Jet_mass_jesFlavorPureCharmDown", Jet_mass_jesFlavorPureCharmDown, &b_Jet_mass_jesFlavorPureCharmDown);
   T1->SetBranchAddress("MET_pt_jesFlavorPureCharmDown", &MET_pt_jesFlavorPureCharmDown, &b_MET_pt_jesFlavorPureCharmDown);
   T1->SetBranchAddress("MET_phi_jesFlavorPureCharmDown", &MET_phi_jesFlavorPureCharmDown, &b_MET_phi_jesFlavorPureCharmDown);
   T1->SetBranchAddress("Jet_pt_jesFlavorPureBottomDown", Jet_pt_jesFlavorPureBottomDown, &b_Jet_pt_jesFlavorPureBottomDown);
   T1->SetBranchAddress("Jet_mass_jesFlavorPureBottomDown", Jet_mass_jesFlavorPureBottomDown, &b_Jet_mass_jesFlavorPureBottomDown);
   T1->SetBranchAddress("MET_pt_jesFlavorPureBottomDown", &MET_pt_jesFlavorPureBottomDown, &b_MET_pt_jesFlavorPureBottomDown);
   T1->SetBranchAddress("MET_phi_jesFlavorPureBottomDown", &MET_phi_jesFlavorPureBottomDown, &b_MET_phi_jesFlavorPureBottomDown);
   T1->SetBranchAddress("Jet_pt_jesTimeRunBDown", Jet_pt_jesTimeRunBDown, &b_Jet_pt_jesTimeRunBDown);
   T1->SetBranchAddress("Jet_mass_jesTimeRunBDown", Jet_mass_jesTimeRunBDown, &b_Jet_mass_jesTimeRunBDown);
   T1->SetBranchAddress("MET_pt_jesTimeRunBDown", &MET_pt_jesTimeRunBDown, &b_MET_pt_jesTimeRunBDown);
   T1->SetBranchAddress("MET_phi_jesTimeRunBDown", &MET_phi_jesTimeRunBDown, &b_MET_phi_jesTimeRunBDown);
   T1->SetBranchAddress("Jet_pt_jesTimeRunCDown", Jet_pt_jesTimeRunCDown, &b_Jet_pt_jesTimeRunCDown);
   T1->SetBranchAddress("Jet_mass_jesTimeRunCDown", Jet_mass_jesTimeRunCDown, &b_Jet_mass_jesTimeRunCDown);
   T1->SetBranchAddress("MET_pt_jesTimeRunCDown", &MET_pt_jesTimeRunCDown, &b_MET_pt_jesTimeRunCDown);
   T1->SetBranchAddress("MET_phi_jesTimeRunCDown", &MET_phi_jesTimeRunCDown, &b_MET_phi_jesTimeRunCDown);
   T1->SetBranchAddress("Jet_pt_jesTimeRunDEDown", Jet_pt_jesTimeRunDEDown, &b_Jet_pt_jesTimeRunDEDown);
   T1->SetBranchAddress("Jet_mass_jesTimeRunDEDown", Jet_mass_jesTimeRunDEDown, &b_Jet_mass_jesTimeRunDEDown);
   T1->SetBranchAddress("MET_pt_jesTimeRunDEDown", &MET_pt_jesTimeRunDEDown, &b_MET_pt_jesTimeRunDEDown);
   T1->SetBranchAddress("MET_phi_jesTimeRunDEDown", &MET_phi_jesTimeRunDEDown, &b_MET_phi_jesTimeRunDEDown);
   T1->SetBranchAddress("Jet_pt_jesTimeRunFDown", Jet_pt_jesTimeRunFDown, &b_Jet_pt_jesTimeRunFDown);
   T1->SetBranchAddress("Jet_mass_jesTimeRunFDown", Jet_mass_jesTimeRunFDown, &b_Jet_mass_jesTimeRunFDown);
   T1->SetBranchAddress("MET_pt_jesTimeRunFDown", &MET_pt_jesTimeRunFDown, &b_MET_pt_jesTimeRunFDown);
   T1->SetBranchAddress("MET_phi_jesTimeRunFDown", &MET_phi_jesTimeRunFDown, &b_MET_phi_jesTimeRunFDown);
   T1->SetBranchAddress("Jet_pt_jesCorrelationGroupMPFInSituDown", Jet_pt_jesCorrelationGroupMPFInSituDown, &b_Jet_pt_jesCorrelationGroupMPFInSituDown);
   T1->SetBranchAddress("Jet_mass_jesCorrelationGroupMPFInSituDown", Jet_mass_jesCorrelationGroupMPFInSituDown, &b_Jet_mass_jesCorrelationGroupMPFInSituDown);
   T1->SetBranchAddress("MET_pt_jesCorrelationGroupMPFInSituDown", &MET_pt_jesCorrelationGroupMPFInSituDown, &b_MET_pt_jesCorrelationGroupMPFInSituDown);
   T1->SetBranchAddress("MET_phi_jesCorrelationGroupMPFInSituDown", &MET_phi_jesCorrelationGroupMPFInSituDown, &b_MET_phi_jesCorrelationGroupMPFInSituDown);
   T1->SetBranchAddress("Jet_pt_jesCorrelationGroupIntercalibrationDown", Jet_pt_jesCorrelationGroupIntercalibrationDown, &b_Jet_pt_jesCorrelationGroupIntercalibrationDown);
   T1->SetBranchAddress("Jet_mass_jesCorrelationGroupIntercalibrationDown", Jet_mass_jesCorrelationGroupIntercalibrationDown, &b_Jet_mass_jesCorrelationGroupIntercalibrationDown);
   T1->SetBranchAddress("MET_pt_jesCorrelationGroupIntercalibrationDown", &MET_pt_jesCorrelationGroupIntercalibrationDown, &b_MET_pt_jesCorrelationGroupIntercalibrationDown);
   T1->SetBranchAddress("MET_phi_jesCorrelationGroupIntercalibrationDown", &MET_phi_jesCorrelationGroupIntercalibrationDown, &b_MET_phi_jesCorrelationGroupIntercalibrationDown);
   T1->SetBranchAddress("Jet_pt_jesCorrelationGroupbJESDown", Jet_pt_jesCorrelationGroupbJESDown, &b_Jet_pt_jesCorrelationGroupbJESDown);
   T1->SetBranchAddress("Jet_mass_jesCorrelationGroupbJESDown", Jet_mass_jesCorrelationGroupbJESDown, &b_Jet_mass_jesCorrelationGroupbJESDown);
   T1->SetBranchAddress("MET_pt_jesCorrelationGroupbJESDown", &MET_pt_jesCorrelationGroupbJESDown, &b_MET_pt_jesCorrelationGroupbJESDown);
   T1->SetBranchAddress("MET_phi_jesCorrelationGroupbJESDown", &MET_phi_jesCorrelationGroupbJESDown, &b_MET_phi_jesCorrelationGroupbJESDown);
   T1->SetBranchAddress("Jet_pt_jesCorrelationGroupFlavorDown", Jet_pt_jesCorrelationGroupFlavorDown, &b_Jet_pt_jesCorrelationGroupFlavorDown);
   T1->SetBranchAddress("Jet_mass_jesCorrelationGroupFlavorDown", Jet_mass_jesCorrelationGroupFlavorDown, &b_Jet_mass_jesCorrelationGroupFlavorDown);
   T1->SetBranchAddress("MET_pt_jesCorrelationGroupFlavorDown", &MET_pt_jesCorrelationGroupFlavorDown, &b_MET_pt_jesCorrelationGroupFlavorDown);
   T1->SetBranchAddress("MET_phi_jesCorrelationGroupFlavorDown", &MET_phi_jesCorrelationGroupFlavorDown, &b_MET_phi_jesCorrelationGroupFlavorDown);
   T1->SetBranchAddress("Jet_pt_jesCorrelationGroupUncorrelatedDown", Jet_pt_jesCorrelationGroupUncorrelatedDown, &b_Jet_pt_jesCorrelationGroupUncorrelatedDown);
   T1->SetBranchAddress("Jet_mass_jesCorrelationGroupUncorrelatedDown", Jet_mass_jesCorrelationGroupUncorrelatedDown, &b_Jet_mass_jesCorrelationGroupUncorrelatedDown);
   T1->SetBranchAddress("MET_pt_jesCorrelationGroupUncorrelatedDown", &MET_pt_jesCorrelationGroupUncorrelatedDown, &b_MET_pt_jesCorrelationGroupUncorrelatedDown);
   T1->SetBranchAddress("MET_phi_jesCorrelationGroupUncorrelatedDown", &MET_phi_jesCorrelationGroupUncorrelatedDown, &b_MET_phi_jesCorrelationGroupUncorrelatedDown);
   T1->SetBranchAddress("MET_pt_unclustEnDown", &MET_pt_unclustEnDown, &b_MET_pt_unclustEnDown);
   T1->SetBranchAddress("MET_phi_unclustEnDown", &MET_phi_unclustEnDown, &b_MET_phi_unclustEnDown);
   T1->SetBranchAddress("FatJet_pt_raw", FatJet_pt_raw, &b_FatJet_pt_raw);
   T1->SetBranchAddress("FatJet_pt_nom", FatJet_pt_nom, &b_FatJet_pt_nom);
   T1->SetBranchAddress("FatJet_mass_raw", FatJet_mass_raw, &b_FatJet_mass_raw);
   T1->SetBranchAddress("FatJet_mass_nom", FatJet_mass_nom, &b_FatJet_mass_nom);
   T1->SetBranchAddress("FatJet_corr_JEC", FatJet_corr_JEC, &b_FatJet_corr_JEC);
   T1->SetBranchAddress("FatJet_corr_JER", FatJet_corr_JER, &b_FatJet_corr_JER);
   T1->SetBranchAddress("FatJet_corr_JMS", FatJet_corr_JMS, &b_FatJet_corr_JMS);
   T1->SetBranchAddress("FatJet_corr_JMR", FatJet_corr_JMR, &b_FatJet_corr_JMR);
   T1->SetBranchAddress("FatJet_msoftdrop_raw", FatJet_msoftdrop_raw, &b_FatJet_msoftdrop_raw);
   T1->SetBranchAddress("FatJet_msoftdrop_nom", FatJet_msoftdrop_nom, &b_FatJet_msoftdrop_nom);
   T1->SetBranchAddress("FatJet_msoftdrop_corr_JMR", FatJet_msoftdrop_corr_JMR, &b_FatJet_msoftdrop_corr_JMR);
   T1->SetBranchAddress("FatJet_msoftdrop_corr_JMS", FatJet_msoftdrop_corr_JMS, &b_FatJet_msoftdrop_corr_JMS);
   T1->SetBranchAddress("FatJet_msoftdrop_corr_PUPPI", FatJet_msoftdrop_corr_PUPPI, &b_FatJet_msoftdrop_corr_PUPPI);
   T1->SetBranchAddress("FatJet_msoftdrop_tau21DDT_nom", FatJet_msoftdrop_tau21DDT_nom, &b_FatJet_msoftdrop_tau21DDT_nom);
   T1->SetBranchAddress("FatJet_pt_jerUp", FatJet_pt_jerUp, &b_FatJet_pt_jerUp);
   T1->SetBranchAddress("FatJet_mass_jerUp", FatJet_mass_jerUp, &b_FatJet_mass_jerUp);
   T1->SetBranchAddress("FatJet_mass_jmrUp", FatJet_mass_jmrUp, &b_FatJet_mass_jmrUp);
   T1->SetBranchAddress("FatJet_mass_jmsUp", FatJet_mass_jmsUp, &b_FatJet_mass_jmsUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jerUp", FatJet_msoftdrop_jerUp, &b_FatJet_msoftdrop_jerUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jmrUp", FatJet_msoftdrop_jmrUp, &b_FatJet_msoftdrop_jmrUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jmsUp", FatJet_msoftdrop_jmsUp, &b_FatJet_msoftdrop_jmsUp);
   T1->SetBranchAddress("FatJet_msoftdrop_tau21DDT_jerUp", FatJet_msoftdrop_tau21DDT_jerUp, &b_FatJet_msoftdrop_tau21DDT_jerUp);
   T1->SetBranchAddress("FatJet_msoftdrop_tau21DDT_jmrUp", FatJet_msoftdrop_tau21DDT_jmrUp, &b_FatJet_msoftdrop_tau21DDT_jmrUp);
   T1->SetBranchAddress("FatJet_msoftdrop_tau21DDT_jmsUp", FatJet_msoftdrop_tau21DDT_jmsUp, &b_FatJet_msoftdrop_tau21DDT_jmsUp);
   T1->SetBranchAddress("FatJet_pt_jesAbsoluteStatUp", FatJet_pt_jesAbsoluteStatUp, &b_FatJet_pt_jesAbsoluteStatUp);
   T1->SetBranchAddress("FatJet_mass_jesAbsoluteStatUp", FatJet_mass_jesAbsoluteStatUp, &b_FatJet_mass_jesAbsoluteStatUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteStatUp", FatJet_msoftdrop_jesAbsoluteStatUp, &b_FatJet_msoftdrop_jesAbsoluteStatUp);
   T1->SetBranchAddress("FatJet_pt_jesAbsoluteScaleUp", FatJet_pt_jesAbsoluteScaleUp, &b_FatJet_pt_jesAbsoluteScaleUp);
   T1->SetBranchAddress("FatJet_mass_jesAbsoluteScaleUp", FatJet_mass_jesAbsoluteScaleUp, &b_FatJet_mass_jesAbsoluteScaleUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteScaleUp", FatJet_msoftdrop_jesAbsoluteScaleUp, &b_FatJet_msoftdrop_jesAbsoluteScaleUp);
   T1->SetBranchAddress("FatJet_pt_jesAbsoluteFlavMapUp", FatJet_pt_jesAbsoluteFlavMapUp, &b_FatJet_pt_jesAbsoluteFlavMapUp);
   T1->SetBranchAddress("FatJet_mass_jesAbsoluteFlavMapUp", FatJet_mass_jesAbsoluteFlavMapUp, &b_FatJet_mass_jesAbsoluteFlavMapUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteFlavMapUp", FatJet_msoftdrop_jesAbsoluteFlavMapUp, &b_FatJet_msoftdrop_jesAbsoluteFlavMapUp);
   T1->SetBranchAddress("FatJet_pt_jesAbsoluteMPFBiasUp", FatJet_pt_jesAbsoluteMPFBiasUp, &b_FatJet_pt_jesAbsoluteMPFBiasUp);
   T1->SetBranchAddress("FatJet_mass_jesAbsoluteMPFBiasUp", FatJet_mass_jesAbsoluteMPFBiasUp, &b_FatJet_mass_jesAbsoluteMPFBiasUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteMPFBiasUp", FatJet_msoftdrop_jesAbsoluteMPFBiasUp, &b_FatJet_msoftdrop_jesAbsoluteMPFBiasUp);
   T1->SetBranchAddress("FatJet_pt_jesFragmentationUp", FatJet_pt_jesFragmentationUp, &b_FatJet_pt_jesFragmentationUp);
   T1->SetBranchAddress("FatJet_mass_jesFragmentationUp", FatJet_mass_jesFragmentationUp, &b_FatJet_mass_jesFragmentationUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesFragmentationUp", FatJet_msoftdrop_jesFragmentationUp, &b_FatJet_msoftdrop_jesFragmentationUp);
   T1->SetBranchAddress("FatJet_pt_jesSinglePionECALUp", FatJet_pt_jesSinglePionECALUp, &b_FatJet_pt_jesSinglePionECALUp);
   T1->SetBranchAddress("FatJet_mass_jesSinglePionECALUp", FatJet_mass_jesSinglePionECALUp, &b_FatJet_mass_jesSinglePionECALUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesSinglePionECALUp", FatJet_msoftdrop_jesSinglePionECALUp, &b_FatJet_msoftdrop_jesSinglePionECALUp);
   T1->SetBranchAddress("FatJet_pt_jesSinglePionHCALUp", FatJet_pt_jesSinglePionHCALUp, &b_FatJet_pt_jesSinglePionHCALUp);
   T1->SetBranchAddress("FatJet_mass_jesSinglePionHCALUp", FatJet_mass_jesSinglePionHCALUp, &b_FatJet_mass_jesSinglePionHCALUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesSinglePionHCALUp", FatJet_msoftdrop_jesSinglePionHCALUp, &b_FatJet_msoftdrop_jesSinglePionHCALUp);
   T1->SetBranchAddress("FatJet_pt_jesFlavorQCDUp", FatJet_pt_jesFlavorQCDUp, &b_FatJet_pt_jesFlavorQCDUp);
   T1->SetBranchAddress("FatJet_mass_jesFlavorQCDUp", FatJet_mass_jesFlavorQCDUp, &b_FatJet_mass_jesFlavorQCDUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesFlavorQCDUp", FatJet_msoftdrop_jesFlavorQCDUp, &b_FatJet_msoftdrop_jesFlavorQCDUp);
   T1->SetBranchAddress("FatJet_pt_jesTimePtEtaUp", FatJet_pt_jesTimePtEtaUp, &b_FatJet_pt_jesTimePtEtaUp);
   T1->SetBranchAddress("FatJet_mass_jesTimePtEtaUp", FatJet_mass_jesTimePtEtaUp, &b_FatJet_mass_jesTimePtEtaUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesTimePtEtaUp", FatJet_msoftdrop_jesTimePtEtaUp, &b_FatJet_msoftdrop_jesTimePtEtaUp);
   T1->SetBranchAddress("FatJet_pt_jesRelativeJEREC1Up", FatJet_pt_jesRelativeJEREC1Up, &b_FatJet_pt_jesRelativeJEREC1Up);
   T1->SetBranchAddress("FatJet_mass_jesRelativeJEREC1Up", FatJet_mass_jesRelativeJEREC1Up, &b_FatJet_mass_jesRelativeJEREC1Up);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativeJEREC1Up", FatJet_msoftdrop_jesRelativeJEREC1Up, &b_FatJet_msoftdrop_jesRelativeJEREC1Up);
   T1->SetBranchAddress("FatJet_pt_jesRelativeJEREC2Up", FatJet_pt_jesRelativeJEREC2Up, &b_FatJet_pt_jesRelativeJEREC2Up);
   T1->SetBranchAddress("FatJet_mass_jesRelativeJEREC2Up", FatJet_mass_jesRelativeJEREC2Up, &b_FatJet_mass_jesRelativeJEREC2Up);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativeJEREC2Up", FatJet_msoftdrop_jesRelativeJEREC2Up, &b_FatJet_msoftdrop_jesRelativeJEREC2Up);
   T1->SetBranchAddress("FatJet_pt_jesRelativeJERHFUp", FatJet_pt_jesRelativeJERHFUp, &b_FatJet_pt_jesRelativeJERHFUp);
   T1->SetBranchAddress("FatJet_mass_jesRelativeJERHFUp", FatJet_mass_jesRelativeJERHFUp, &b_FatJet_mass_jesRelativeJERHFUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativeJERHFUp", FatJet_msoftdrop_jesRelativeJERHFUp, &b_FatJet_msoftdrop_jesRelativeJERHFUp);
   T1->SetBranchAddress("FatJet_pt_jesRelativePtBBUp", FatJet_pt_jesRelativePtBBUp, &b_FatJet_pt_jesRelativePtBBUp);
   T1->SetBranchAddress("FatJet_mass_jesRelativePtBBUp", FatJet_mass_jesRelativePtBBUp, &b_FatJet_mass_jesRelativePtBBUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativePtBBUp", FatJet_msoftdrop_jesRelativePtBBUp, &b_FatJet_msoftdrop_jesRelativePtBBUp);
   T1->SetBranchAddress("FatJet_pt_jesRelativePtEC1Up", FatJet_pt_jesRelativePtEC1Up, &b_FatJet_pt_jesRelativePtEC1Up);
   T1->SetBranchAddress("FatJet_mass_jesRelativePtEC1Up", FatJet_mass_jesRelativePtEC1Up, &b_FatJet_mass_jesRelativePtEC1Up);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativePtEC1Up", FatJet_msoftdrop_jesRelativePtEC1Up, &b_FatJet_msoftdrop_jesRelativePtEC1Up);
   T1->SetBranchAddress("FatJet_pt_jesRelativePtEC2Up", FatJet_pt_jesRelativePtEC2Up, &b_FatJet_pt_jesRelativePtEC2Up);
   T1->SetBranchAddress("FatJet_mass_jesRelativePtEC2Up", FatJet_mass_jesRelativePtEC2Up, &b_FatJet_mass_jesRelativePtEC2Up);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativePtEC2Up", FatJet_msoftdrop_jesRelativePtEC2Up, &b_FatJet_msoftdrop_jesRelativePtEC2Up);
   T1->SetBranchAddress("FatJet_pt_jesRelativePtHFUp", FatJet_pt_jesRelativePtHFUp, &b_FatJet_pt_jesRelativePtHFUp);
   T1->SetBranchAddress("FatJet_mass_jesRelativePtHFUp", FatJet_mass_jesRelativePtHFUp, &b_FatJet_mass_jesRelativePtHFUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativePtHFUp", FatJet_msoftdrop_jesRelativePtHFUp, &b_FatJet_msoftdrop_jesRelativePtHFUp);
   T1->SetBranchAddress("FatJet_pt_jesRelativeBalUp", FatJet_pt_jesRelativeBalUp, &b_FatJet_pt_jesRelativeBalUp);
   T1->SetBranchAddress("FatJet_mass_jesRelativeBalUp", FatJet_mass_jesRelativeBalUp, &b_FatJet_mass_jesRelativeBalUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativeBalUp", FatJet_msoftdrop_jesRelativeBalUp, &b_FatJet_msoftdrop_jesRelativeBalUp);
   T1->SetBranchAddress("FatJet_pt_jesRelativeSampleUp", FatJet_pt_jesRelativeSampleUp, &b_FatJet_pt_jesRelativeSampleUp);
   T1->SetBranchAddress("FatJet_mass_jesRelativeSampleUp", FatJet_mass_jesRelativeSampleUp, &b_FatJet_mass_jesRelativeSampleUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativeSampleUp", FatJet_msoftdrop_jesRelativeSampleUp, &b_FatJet_msoftdrop_jesRelativeSampleUp);
   T1->SetBranchAddress("FatJet_pt_jesRelativeFSRUp", FatJet_pt_jesRelativeFSRUp, &b_FatJet_pt_jesRelativeFSRUp);
   T1->SetBranchAddress("FatJet_mass_jesRelativeFSRUp", FatJet_mass_jesRelativeFSRUp, &b_FatJet_mass_jesRelativeFSRUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativeFSRUp", FatJet_msoftdrop_jesRelativeFSRUp, &b_FatJet_msoftdrop_jesRelativeFSRUp);
   T1->SetBranchAddress("FatJet_pt_jesRelativeStatFSRUp", FatJet_pt_jesRelativeStatFSRUp, &b_FatJet_pt_jesRelativeStatFSRUp);
   T1->SetBranchAddress("FatJet_mass_jesRelativeStatFSRUp", FatJet_mass_jesRelativeStatFSRUp, &b_FatJet_mass_jesRelativeStatFSRUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativeStatFSRUp", FatJet_msoftdrop_jesRelativeStatFSRUp, &b_FatJet_msoftdrop_jesRelativeStatFSRUp);
   T1->SetBranchAddress("FatJet_pt_jesRelativeStatECUp", FatJet_pt_jesRelativeStatECUp, &b_FatJet_pt_jesRelativeStatECUp);
   T1->SetBranchAddress("FatJet_mass_jesRelativeStatECUp", FatJet_mass_jesRelativeStatECUp, &b_FatJet_mass_jesRelativeStatECUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativeStatECUp", FatJet_msoftdrop_jesRelativeStatECUp, &b_FatJet_msoftdrop_jesRelativeStatECUp);
   T1->SetBranchAddress("FatJet_pt_jesRelativeStatHFUp", FatJet_pt_jesRelativeStatHFUp, &b_FatJet_pt_jesRelativeStatHFUp);
   T1->SetBranchAddress("FatJet_mass_jesRelativeStatHFUp", FatJet_mass_jesRelativeStatHFUp, &b_FatJet_mass_jesRelativeStatHFUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativeStatHFUp", FatJet_msoftdrop_jesRelativeStatHFUp, &b_FatJet_msoftdrop_jesRelativeStatHFUp);
   T1->SetBranchAddress("FatJet_pt_jesPileUpDataMCUp", FatJet_pt_jesPileUpDataMCUp, &b_FatJet_pt_jesPileUpDataMCUp);
   T1->SetBranchAddress("FatJet_mass_jesPileUpDataMCUp", FatJet_mass_jesPileUpDataMCUp, &b_FatJet_mass_jesPileUpDataMCUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesPileUpDataMCUp", FatJet_msoftdrop_jesPileUpDataMCUp, &b_FatJet_msoftdrop_jesPileUpDataMCUp);
   T1->SetBranchAddress("FatJet_pt_jesPileUpPtRefUp", FatJet_pt_jesPileUpPtRefUp, &b_FatJet_pt_jesPileUpPtRefUp);
   T1->SetBranchAddress("FatJet_mass_jesPileUpPtRefUp", FatJet_mass_jesPileUpPtRefUp, &b_FatJet_mass_jesPileUpPtRefUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtRefUp", FatJet_msoftdrop_jesPileUpPtRefUp, &b_FatJet_msoftdrop_jesPileUpPtRefUp);
   T1->SetBranchAddress("FatJet_pt_jesPileUpPtBBUp", FatJet_pt_jesPileUpPtBBUp, &b_FatJet_pt_jesPileUpPtBBUp);
   T1->SetBranchAddress("FatJet_mass_jesPileUpPtBBUp", FatJet_mass_jesPileUpPtBBUp, &b_FatJet_mass_jesPileUpPtBBUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtBBUp", FatJet_msoftdrop_jesPileUpPtBBUp, &b_FatJet_msoftdrop_jesPileUpPtBBUp);
   T1->SetBranchAddress("FatJet_pt_jesPileUpPtEC1Up", FatJet_pt_jesPileUpPtEC1Up, &b_FatJet_pt_jesPileUpPtEC1Up);
   T1->SetBranchAddress("FatJet_mass_jesPileUpPtEC1Up", FatJet_mass_jesPileUpPtEC1Up, &b_FatJet_mass_jesPileUpPtEC1Up);
   T1->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtEC1Up", FatJet_msoftdrop_jesPileUpPtEC1Up, &b_FatJet_msoftdrop_jesPileUpPtEC1Up);
   T1->SetBranchAddress("FatJet_pt_jesPileUpPtEC2Up", FatJet_pt_jesPileUpPtEC2Up, &b_FatJet_pt_jesPileUpPtEC2Up);
   T1->SetBranchAddress("FatJet_mass_jesPileUpPtEC2Up", FatJet_mass_jesPileUpPtEC2Up, &b_FatJet_mass_jesPileUpPtEC2Up);
   T1->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtEC2Up", FatJet_msoftdrop_jesPileUpPtEC2Up, &b_FatJet_msoftdrop_jesPileUpPtEC2Up);
   T1->SetBranchAddress("FatJet_pt_jesPileUpPtHFUp", FatJet_pt_jesPileUpPtHFUp, &b_FatJet_pt_jesPileUpPtHFUp);
   T1->SetBranchAddress("FatJet_mass_jesPileUpPtHFUp", FatJet_mass_jesPileUpPtHFUp, &b_FatJet_mass_jesPileUpPtHFUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtHFUp", FatJet_msoftdrop_jesPileUpPtHFUp, &b_FatJet_msoftdrop_jesPileUpPtHFUp);
   T1->SetBranchAddress("FatJet_pt_jesPileUpMuZeroUp", FatJet_pt_jesPileUpMuZeroUp, &b_FatJet_pt_jesPileUpMuZeroUp);
   T1->SetBranchAddress("FatJet_mass_jesPileUpMuZeroUp", FatJet_mass_jesPileUpMuZeroUp, &b_FatJet_mass_jesPileUpMuZeroUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesPileUpMuZeroUp", FatJet_msoftdrop_jesPileUpMuZeroUp, &b_FatJet_msoftdrop_jesPileUpMuZeroUp);
   T1->SetBranchAddress("FatJet_pt_jesPileUpEnvelopeUp", FatJet_pt_jesPileUpEnvelopeUp, &b_FatJet_pt_jesPileUpEnvelopeUp);
   T1->SetBranchAddress("FatJet_mass_jesPileUpEnvelopeUp", FatJet_mass_jesPileUpEnvelopeUp, &b_FatJet_mass_jesPileUpEnvelopeUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesPileUpEnvelopeUp", FatJet_msoftdrop_jesPileUpEnvelopeUp, &b_FatJet_msoftdrop_jesPileUpEnvelopeUp);
   T1->SetBranchAddress("FatJet_pt_jesSubTotalPileUpUp", FatJet_pt_jesSubTotalPileUpUp, &b_FatJet_pt_jesSubTotalPileUpUp);
   T1->SetBranchAddress("FatJet_mass_jesSubTotalPileUpUp", FatJet_mass_jesSubTotalPileUpUp, &b_FatJet_mass_jesSubTotalPileUpUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesSubTotalPileUpUp", FatJet_msoftdrop_jesSubTotalPileUpUp, &b_FatJet_msoftdrop_jesSubTotalPileUpUp);
   T1->SetBranchAddress("FatJet_pt_jesSubTotalRelativeUp", FatJet_pt_jesSubTotalRelativeUp, &b_FatJet_pt_jesSubTotalRelativeUp);
   T1->SetBranchAddress("FatJet_mass_jesSubTotalRelativeUp", FatJet_mass_jesSubTotalRelativeUp, &b_FatJet_mass_jesSubTotalRelativeUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesSubTotalRelativeUp", FatJet_msoftdrop_jesSubTotalRelativeUp, &b_FatJet_msoftdrop_jesSubTotalRelativeUp);
   T1->SetBranchAddress("FatJet_pt_jesSubTotalPtUp", FatJet_pt_jesSubTotalPtUp, &b_FatJet_pt_jesSubTotalPtUp);
   T1->SetBranchAddress("FatJet_mass_jesSubTotalPtUp", FatJet_mass_jesSubTotalPtUp, &b_FatJet_mass_jesSubTotalPtUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesSubTotalPtUp", FatJet_msoftdrop_jesSubTotalPtUp, &b_FatJet_msoftdrop_jesSubTotalPtUp);
   T1->SetBranchAddress("FatJet_pt_jesSubTotalScaleUp", FatJet_pt_jesSubTotalScaleUp, &b_FatJet_pt_jesSubTotalScaleUp);
   T1->SetBranchAddress("FatJet_mass_jesSubTotalScaleUp", FatJet_mass_jesSubTotalScaleUp, &b_FatJet_mass_jesSubTotalScaleUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesSubTotalScaleUp", FatJet_msoftdrop_jesSubTotalScaleUp, &b_FatJet_msoftdrop_jesSubTotalScaleUp);
   T1->SetBranchAddress("FatJet_pt_jesSubTotalAbsoluteUp", FatJet_pt_jesSubTotalAbsoluteUp, &b_FatJet_pt_jesSubTotalAbsoluteUp);
   T1->SetBranchAddress("FatJet_mass_jesSubTotalAbsoluteUp", FatJet_mass_jesSubTotalAbsoluteUp, &b_FatJet_mass_jesSubTotalAbsoluteUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesSubTotalAbsoluteUp", FatJet_msoftdrop_jesSubTotalAbsoluteUp, &b_FatJet_msoftdrop_jesSubTotalAbsoluteUp);
   T1->SetBranchAddress("FatJet_pt_jesSubTotalMCUp", FatJet_pt_jesSubTotalMCUp, &b_FatJet_pt_jesSubTotalMCUp);
   T1->SetBranchAddress("FatJet_mass_jesSubTotalMCUp", FatJet_mass_jesSubTotalMCUp, &b_FatJet_mass_jesSubTotalMCUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesSubTotalMCUp", FatJet_msoftdrop_jesSubTotalMCUp, &b_FatJet_msoftdrop_jesSubTotalMCUp);
   T1->SetBranchAddress("FatJet_pt_jesTotalUp", FatJet_pt_jesTotalUp, &b_FatJet_pt_jesTotalUp);
   T1->SetBranchAddress("FatJet_mass_jesTotalUp", FatJet_mass_jesTotalUp, &b_FatJet_mass_jesTotalUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesTotalUp", FatJet_msoftdrop_jesTotalUp, &b_FatJet_msoftdrop_jesTotalUp);
   T1->SetBranchAddress("FatJet_pt_jesTotalNoFlavorUp", FatJet_pt_jesTotalNoFlavorUp, &b_FatJet_pt_jesTotalNoFlavorUp);
   T1->SetBranchAddress("FatJet_mass_jesTotalNoFlavorUp", FatJet_mass_jesTotalNoFlavorUp, &b_FatJet_mass_jesTotalNoFlavorUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesTotalNoFlavorUp", FatJet_msoftdrop_jesTotalNoFlavorUp, &b_FatJet_msoftdrop_jesTotalNoFlavorUp);
   T1->SetBranchAddress("FatJet_pt_jesTotalNoTimeUp", FatJet_pt_jesTotalNoTimeUp, &b_FatJet_pt_jesTotalNoTimeUp);
   T1->SetBranchAddress("FatJet_mass_jesTotalNoTimeUp", FatJet_mass_jesTotalNoTimeUp, &b_FatJet_mass_jesTotalNoTimeUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesTotalNoTimeUp", FatJet_msoftdrop_jesTotalNoTimeUp, &b_FatJet_msoftdrop_jesTotalNoTimeUp);
   T1->SetBranchAddress("FatJet_pt_jesTotalNoFlavorNoTimeUp", FatJet_pt_jesTotalNoFlavorNoTimeUp, &b_FatJet_pt_jesTotalNoFlavorNoTimeUp);
   T1->SetBranchAddress("FatJet_mass_jesTotalNoFlavorNoTimeUp", FatJet_mass_jesTotalNoFlavorNoTimeUp, &b_FatJet_mass_jesTotalNoFlavorNoTimeUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesTotalNoFlavorNoTimeUp", FatJet_msoftdrop_jesTotalNoFlavorNoTimeUp, &b_FatJet_msoftdrop_jesTotalNoFlavorNoTimeUp);
   T1->SetBranchAddress("FatJet_pt_jesFlavorZJetUp", FatJet_pt_jesFlavorZJetUp, &b_FatJet_pt_jesFlavorZJetUp);
   T1->SetBranchAddress("FatJet_mass_jesFlavorZJetUp", FatJet_mass_jesFlavorZJetUp, &b_FatJet_mass_jesFlavorZJetUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesFlavorZJetUp", FatJet_msoftdrop_jesFlavorZJetUp, &b_FatJet_msoftdrop_jesFlavorZJetUp);
   T1->SetBranchAddress("FatJet_pt_jesFlavorPhotonJetUp", FatJet_pt_jesFlavorPhotonJetUp, &b_FatJet_pt_jesFlavorPhotonJetUp);
   T1->SetBranchAddress("FatJet_mass_jesFlavorPhotonJetUp", FatJet_mass_jesFlavorPhotonJetUp, &b_FatJet_mass_jesFlavorPhotonJetUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesFlavorPhotonJetUp", FatJet_msoftdrop_jesFlavorPhotonJetUp, &b_FatJet_msoftdrop_jesFlavorPhotonJetUp);
   T1->SetBranchAddress("FatJet_pt_jesFlavorPureGluonUp", FatJet_pt_jesFlavorPureGluonUp, &b_FatJet_pt_jesFlavorPureGluonUp);
   T1->SetBranchAddress("FatJet_mass_jesFlavorPureGluonUp", FatJet_mass_jesFlavorPureGluonUp, &b_FatJet_mass_jesFlavorPureGluonUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureGluonUp", FatJet_msoftdrop_jesFlavorPureGluonUp, &b_FatJet_msoftdrop_jesFlavorPureGluonUp);
   T1->SetBranchAddress("FatJet_pt_jesFlavorPureQuarkUp", FatJet_pt_jesFlavorPureQuarkUp, &b_FatJet_pt_jesFlavorPureQuarkUp);
   T1->SetBranchAddress("FatJet_mass_jesFlavorPureQuarkUp", FatJet_mass_jesFlavorPureQuarkUp, &b_FatJet_mass_jesFlavorPureQuarkUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureQuarkUp", FatJet_msoftdrop_jesFlavorPureQuarkUp, &b_FatJet_msoftdrop_jesFlavorPureQuarkUp);
   T1->SetBranchAddress("FatJet_pt_jesFlavorPureCharmUp", FatJet_pt_jesFlavorPureCharmUp, &b_FatJet_pt_jesFlavorPureCharmUp);
   T1->SetBranchAddress("FatJet_mass_jesFlavorPureCharmUp", FatJet_mass_jesFlavorPureCharmUp, &b_FatJet_mass_jesFlavorPureCharmUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureCharmUp", FatJet_msoftdrop_jesFlavorPureCharmUp, &b_FatJet_msoftdrop_jesFlavorPureCharmUp);
   T1->SetBranchAddress("FatJet_pt_jesFlavorPureBottomUp", FatJet_pt_jesFlavorPureBottomUp, &b_FatJet_pt_jesFlavorPureBottomUp);
   T1->SetBranchAddress("FatJet_mass_jesFlavorPureBottomUp", FatJet_mass_jesFlavorPureBottomUp, &b_FatJet_mass_jesFlavorPureBottomUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureBottomUp", FatJet_msoftdrop_jesFlavorPureBottomUp, &b_FatJet_msoftdrop_jesFlavorPureBottomUp);
   T1->SetBranchAddress("FatJet_pt_jesTimeRunBUp", FatJet_pt_jesTimeRunBUp, &b_FatJet_pt_jesTimeRunBUp);
   T1->SetBranchAddress("FatJet_mass_jesTimeRunBUp", FatJet_mass_jesTimeRunBUp, &b_FatJet_mass_jesTimeRunBUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesTimeRunBUp", FatJet_msoftdrop_jesTimeRunBUp, &b_FatJet_msoftdrop_jesTimeRunBUp);
   T1->SetBranchAddress("FatJet_pt_jesTimeRunCUp", FatJet_pt_jesTimeRunCUp, &b_FatJet_pt_jesTimeRunCUp);
   T1->SetBranchAddress("FatJet_mass_jesTimeRunCUp", FatJet_mass_jesTimeRunCUp, &b_FatJet_mass_jesTimeRunCUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesTimeRunCUp", FatJet_msoftdrop_jesTimeRunCUp, &b_FatJet_msoftdrop_jesTimeRunCUp);
   T1->SetBranchAddress("FatJet_pt_jesTimeRunDEUp", FatJet_pt_jesTimeRunDEUp, &b_FatJet_pt_jesTimeRunDEUp);
   T1->SetBranchAddress("FatJet_mass_jesTimeRunDEUp", FatJet_mass_jesTimeRunDEUp, &b_FatJet_mass_jesTimeRunDEUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesTimeRunDEUp", FatJet_msoftdrop_jesTimeRunDEUp, &b_FatJet_msoftdrop_jesTimeRunDEUp);
   T1->SetBranchAddress("FatJet_pt_jesTimeRunFUp", FatJet_pt_jesTimeRunFUp, &b_FatJet_pt_jesTimeRunFUp);
   T1->SetBranchAddress("FatJet_mass_jesTimeRunFUp", FatJet_mass_jesTimeRunFUp, &b_FatJet_mass_jesTimeRunFUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesTimeRunFUp", FatJet_msoftdrop_jesTimeRunFUp, &b_FatJet_msoftdrop_jesTimeRunFUp);
   T1->SetBranchAddress("FatJet_pt_jesCorrelationGroupMPFInSituUp", FatJet_pt_jesCorrelationGroupMPFInSituUp, &b_FatJet_pt_jesCorrelationGroupMPFInSituUp);
   T1->SetBranchAddress("FatJet_mass_jesCorrelationGroupMPFInSituUp", FatJet_mass_jesCorrelationGroupMPFInSituUp, &b_FatJet_mass_jesCorrelationGroupMPFInSituUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupMPFInSituUp", FatJet_msoftdrop_jesCorrelationGroupMPFInSituUp, &b_FatJet_msoftdrop_jesCorrelationGroupMPFInSituUp);
   T1->SetBranchAddress("FatJet_pt_jesCorrelationGroupIntercalibrationUp", FatJet_pt_jesCorrelationGroupIntercalibrationUp, &b_FatJet_pt_jesCorrelationGroupIntercalibrationUp);
   T1->SetBranchAddress("FatJet_mass_jesCorrelationGroupIntercalibrationUp", FatJet_mass_jesCorrelationGroupIntercalibrationUp, &b_FatJet_mass_jesCorrelationGroupIntercalibrationUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupIntercalibrationUp", FatJet_msoftdrop_jesCorrelationGroupIntercalibrationUp, &b_FatJet_msoftdrop_jesCorrelationGroupIntercalibrationUp);
   T1->SetBranchAddress("FatJet_pt_jesCorrelationGroupbJESUp", FatJet_pt_jesCorrelationGroupbJESUp, &b_FatJet_pt_jesCorrelationGroupbJESUp);
   T1->SetBranchAddress("FatJet_mass_jesCorrelationGroupbJESUp", FatJet_mass_jesCorrelationGroupbJESUp, &b_FatJet_mass_jesCorrelationGroupbJESUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupbJESUp", FatJet_msoftdrop_jesCorrelationGroupbJESUp, &b_FatJet_msoftdrop_jesCorrelationGroupbJESUp);
   T1->SetBranchAddress("FatJet_pt_jesCorrelationGroupFlavorUp", FatJet_pt_jesCorrelationGroupFlavorUp, &b_FatJet_pt_jesCorrelationGroupFlavorUp);
   T1->SetBranchAddress("FatJet_mass_jesCorrelationGroupFlavorUp", FatJet_mass_jesCorrelationGroupFlavorUp, &b_FatJet_mass_jesCorrelationGroupFlavorUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupFlavorUp", FatJet_msoftdrop_jesCorrelationGroupFlavorUp, &b_FatJet_msoftdrop_jesCorrelationGroupFlavorUp);
   T1->SetBranchAddress("FatJet_pt_jesCorrelationGroupUncorrelatedUp", FatJet_pt_jesCorrelationGroupUncorrelatedUp, &b_FatJet_pt_jesCorrelationGroupUncorrelatedUp);
   T1->SetBranchAddress("FatJet_mass_jesCorrelationGroupUncorrelatedUp", FatJet_mass_jesCorrelationGroupUncorrelatedUp, &b_FatJet_mass_jesCorrelationGroupUncorrelatedUp);
   T1->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupUncorrelatedUp", FatJet_msoftdrop_jesCorrelationGroupUncorrelatedUp, &b_FatJet_msoftdrop_jesCorrelationGroupUncorrelatedUp);
   T1->SetBranchAddress("FatJet_pt_jerDown", FatJet_pt_jerDown, &b_FatJet_pt_jerDown);
   T1->SetBranchAddress("FatJet_mass_jerDown", FatJet_mass_jerDown, &b_FatJet_mass_jerDown);
   T1->SetBranchAddress("FatJet_mass_jmrDown", FatJet_mass_jmrDown, &b_FatJet_mass_jmrDown);
   T1->SetBranchAddress("FatJet_mass_jmsDown", FatJet_mass_jmsDown, &b_FatJet_mass_jmsDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jerDown", FatJet_msoftdrop_jerDown, &b_FatJet_msoftdrop_jerDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jmrDown", FatJet_msoftdrop_jmrDown, &b_FatJet_msoftdrop_jmrDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jmsDown", FatJet_msoftdrop_jmsDown, &b_FatJet_msoftdrop_jmsDown);
   T1->SetBranchAddress("FatJet_msoftdrop_tau21DDT_jerDown", FatJet_msoftdrop_tau21DDT_jerDown, &b_FatJet_msoftdrop_tau21DDT_jerDown);
   T1->SetBranchAddress("FatJet_msoftdrop_tau21DDT_jmrDown", FatJet_msoftdrop_tau21DDT_jmrDown, &b_FatJet_msoftdrop_tau21DDT_jmrDown);
   T1->SetBranchAddress("FatJet_msoftdrop_tau21DDT_jmsDown", FatJet_msoftdrop_tau21DDT_jmsDown, &b_FatJet_msoftdrop_tau21DDT_jmsDown);
   T1->SetBranchAddress("FatJet_pt_jesAbsoluteStatDown", FatJet_pt_jesAbsoluteStatDown, &b_FatJet_pt_jesAbsoluteStatDown);
   T1->SetBranchAddress("FatJet_mass_jesAbsoluteStatDown", FatJet_mass_jesAbsoluteStatDown, &b_FatJet_mass_jesAbsoluteStatDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteStatDown", FatJet_msoftdrop_jesAbsoluteStatDown, &b_FatJet_msoftdrop_jesAbsoluteStatDown);
   T1->SetBranchAddress("FatJet_pt_jesAbsoluteScaleDown", FatJet_pt_jesAbsoluteScaleDown, &b_FatJet_pt_jesAbsoluteScaleDown);
   T1->SetBranchAddress("FatJet_mass_jesAbsoluteScaleDown", FatJet_mass_jesAbsoluteScaleDown, &b_FatJet_mass_jesAbsoluteScaleDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteScaleDown", FatJet_msoftdrop_jesAbsoluteScaleDown, &b_FatJet_msoftdrop_jesAbsoluteScaleDown);
   T1->SetBranchAddress("FatJet_pt_jesAbsoluteFlavMapDown", FatJet_pt_jesAbsoluteFlavMapDown, &b_FatJet_pt_jesAbsoluteFlavMapDown);
   T1->SetBranchAddress("FatJet_mass_jesAbsoluteFlavMapDown", FatJet_mass_jesAbsoluteFlavMapDown, &b_FatJet_mass_jesAbsoluteFlavMapDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteFlavMapDown", FatJet_msoftdrop_jesAbsoluteFlavMapDown, &b_FatJet_msoftdrop_jesAbsoluteFlavMapDown);
   T1->SetBranchAddress("FatJet_pt_jesAbsoluteMPFBiasDown", FatJet_pt_jesAbsoluteMPFBiasDown, &b_FatJet_pt_jesAbsoluteMPFBiasDown);
   T1->SetBranchAddress("FatJet_mass_jesAbsoluteMPFBiasDown", FatJet_mass_jesAbsoluteMPFBiasDown, &b_FatJet_mass_jesAbsoluteMPFBiasDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteMPFBiasDown", FatJet_msoftdrop_jesAbsoluteMPFBiasDown, &b_FatJet_msoftdrop_jesAbsoluteMPFBiasDown);
   T1->SetBranchAddress("FatJet_pt_jesFragmentationDown", FatJet_pt_jesFragmentationDown, &b_FatJet_pt_jesFragmentationDown);
   T1->SetBranchAddress("FatJet_mass_jesFragmentationDown", FatJet_mass_jesFragmentationDown, &b_FatJet_mass_jesFragmentationDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesFragmentationDown", FatJet_msoftdrop_jesFragmentationDown, &b_FatJet_msoftdrop_jesFragmentationDown);
   T1->SetBranchAddress("FatJet_pt_jesSinglePionECALDown", FatJet_pt_jesSinglePionECALDown, &b_FatJet_pt_jesSinglePionECALDown);
   T1->SetBranchAddress("FatJet_mass_jesSinglePionECALDown", FatJet_mass_jesSinglePionECALDown, &b_FatJet_mass_jesSinglePionECALDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesSinglePionECALDown", FatJet_msoftdrop_jesSinglePionECALDown, &b_FatJet_msoftdrop_jesSinglePionECALDown);
   T1->SetBranchAddress("FatJet_pt_jesSinglePionHCALDown", FatJet_pt_jesSinglePionHCALDown, &b_FatJet_pt_jesSinglePionHCALDown);
   T1->SetBranchAddress("FatJet_mass_jesSinglePionHCALDown", FatJet_mass_jesSinglePionHCALDown, &b_FatJet_mass_jesSinglePionHCALDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesSinglePionHCALDown", FatJet_msoftdrop_jesSinglePionHCALDown, &b_FatJet_msoftdrop_jesSinglePionHCALDown);
   T1->SetBranchAddress("FatJet_pt_jesFlavorQCDDown", FatJet_pt_jesFlavorQCDDown, &b_FatJet_pt_jesFlavorQCDDown);
   T1->SetBranchAddress("FatJet_mass_jesFlavorQCDDown", FatJet_mass_jesFlavorQCDDown, &b_FatJet_mass_jesFlavorQCDDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesFlavorQCDDown", FatJet_msoftdrop_jesFlavorQCDDown, &b_FatJet_msoftdrop_jesFlavorQCDDown);
   T1->SetBranchAddress("FatJet_pt_jesTimePtEtaDown", FatJet_pt_jesTimePtEtaDown, &b_FatJet_pt_jesTimePtEtaDown);
   T1->SetBranchAddress("FatJet_mass_jesTimePtEtaDown", FatJet_mass_jesTimePtEtaDown, &b_FatJet_mass_jesTimePtEtaDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesTimePtEtaDown", FatJet_msoftdrop_jesTimePtEtaDown, &b_FatJet_msoftdrop_jesTimePtEtaDown);
   T1->SetBranchAddress("FatJet_pt_jesRelativeJEREC1Down", FatJet_pt_jesRelativeJEREC1Down, &b_FatJet_pt_jesRelativeJEREC1Down);
   T1->SetBranchAddress("FatJet_mass_jesRelativeJEREC1Down", FatJet_mass_jesRelativeJEREC1Down, &b_FatJet_mass_jesRelativeJEREC1Down);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativeJEREC1Down", FatJet_msoftdrop_jesRelativeJEREC1Down, &b_FatJet_msoftdrop_jesRelativeJEREC1Down);
   T1->SetBranchAddress("FatJet_pt_jesRelativeJEREC2Down", FatJet_pt_jesRelativeJEREC2Down, &b_FatJet_pt_jesRelativeJEREC2Down);
   T1->SetBranchAddress("FatJet_mass_jesRelativeJEREC2Down", FatJet_mass_jesRelativeJEREC2Down, &b_FatJet_mass_jesRelativeJEREC2Down);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativeJEREC2Down", FatJet_msoftdrop_jesRelativeJEREC2Down, &b_FatJet_msoftdrop_jesRelativeJEREC2Down);
   T1->SetBranchAddress("FatJet_pt_jesRelativeJERHFDown", FatJet_pt_jesRelativeJERHFDown, &b_FatJet_pt_jesRelativeJERHFDown);
   T1->SetBranchAddress("FatJet_mass_jesRelativeJERHFDown", FatJet_mass_jesRelativeJERHFDown, &b_FatJet_mass_jesRelativeJERHFDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativeJERHFDown", FatJet_msoftdrop_jesRelativeJERHFDown, &b_FatJet_msoftdrop_jesRelativeJERHFDown);
   T1->SetBranchAddress("FatJet_pt_jesRelativePtBBDown", FatJet_pt_jesRelativePtBBDown, &b_FatJet_pt_jesRelativePtBBDown);
   T1->SetBranchAddress("FatJet_mass_jesRelativePtBBDown", FatJet_mass_jesRelativePtBBDown, &b_FatJet_mass_jesRelativePtBBDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativePtBBDown", FatJet_msoftdrop_jesRelativePtBBDown, &b_FatJet_msoftdrop_jesRelativePtBBDown);
   T1->SetBranchAddress("FatJet_pt_jesRelativePtEC1Down", FatJet_pt_jesRelativePtEC1Down, &b_FatJet_pt_jesRelativePtEC1Down);
   T1->SetBranchAddress("FatJet_mass_jesRelativePtEC1Down", FatJet_mass_jesRelativePtEC1Down, &b_FatJet_mass_jesRelativePtEC1Down);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativePtEC1Down", FatJet_msoftdrop_jesRelativePtEC1Down, &b_FatJet_msoftdrop_jesRelativePtEC1Down);
   T1->SetBranchAddress("FatJet_pt_jesRelativePtEC2Down", FatJet_pt_jesRelativePtEC2Down, &b_FatJet_pt_jesRelativePtEC2Down);
   T1->SetBranchAddress("FatJet_mass_jesRelativePtEC2Down", FatJet_mass_jesRelativePtEC2Down, &b_FatJet_mass_jesRelativePtEC2Down);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativePtEC2Down", FatJet_msoftdrop_jesRelativePtEC2Down, &b_FatJet_msoftdrop_jesRelativePtEC2Down);
   T1->SetBranchAddress("FatJet_pt_jesRelativePtHFDown", FatJet_pt_jesRelativePtHFDown, &b_FatJet_pt_jesRelativePtHFDown);
   T1->SetBranchAddress("FatJet_mass_jesRelativePtHFDown", FatJet_mass_jesRelativePtHFDown, &b_FatJet_mass_jesRelativePtHFDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativePtHFDown", FatJet_msoftdrop_jesRelativePtHFDown, &b_FatJet_msoftdrop_jesRelativePtHFDown);
   T1->SetBranchAddress("FatJet_pt_jesRelativeBalDown", FatJet_pt_jesRelativeBalDown, &b_FatJet_pt_jesRelativeBalDown);
   T1->SetBranchAddress("FatJet_mass_jesRelativeBalDown", FatJet_mass_jesRelativeBalDown, &b_FatJet_mass_jesRelativeBalDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativeBalDown", FatJet_msoftdrop_jesRelativeBalDown, &b_FatJet_msoftdrop_jesRelativeBalDown);
   T1->SetBranchAddress("FatJet_pt_jesRelativeSampleDown", FatJet_pt_jesRelativeSampleDown, &b_FatJet_pt_jesRelativeSampleDown);
   T1->SetBranchAddress("FatJet_mass_jesRelativeSampleDown", FatJet_mass_jesRelativeSampleDown, &b_FatJet_mass_jesRelativeSampleDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativeSampleDown", FatJet_msoftdrop_jesRelativeSampleDown, &b_FatJet_msoftdrop_jesRelativeSampleDown);
   T1->SetBranchAddress("FatJet_pt_jesRelativeFSRDown", FatJet_pt_jesRelativeFSRDown, &b_FatJet_pt_jesRelativeFSRDown);
   T1->SetBranchAddress("FatJet_mass_jesRelativeFSRDown", FatJet_mass_jesRelativeFSRDown, &b_FatJet_mass_jesRelativeFSRDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativeFSRDown", FatJet_msoftdrop_jesRelativeFSRDown, &b_FatJet_msoftdrop_jesRelativeFSRDown);
   T1->SetBranchAddress("FatJet_pt_jesRelativeStatFSRDown", FatJet_pt_jesRelativeStatFSRDown, &b_FatJet_pt_jesRelativeStatFSRDown);
   T1->SetBranchAddress("FatJet_mass_jesRelativeStatFSRDown", FatJet_mass_jesRelativeStatFSRDown, &b_FatJet_mass_jesRelativeStatFSRDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativeStatFSRDown", FatJet_msoftdrop_jesRelativeStatFSRDown, &b_FatJet_msoftdrop_jesRelativeStatFSRDown);
   T1->SetBranchAddress("FatJet_pt_jesRelativeStatECDown", FatJet_pt_jesRelativeStatECDown, &b_FatJet_pt_jesRelativeStatECDown);
   T1->SetBranchAddress("FatJet_mass_jesRelativeStatECDown", FatJet_mass_jesRelativeStatECDown, &b_FatJet_mass_jesRelativeStatECDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativeStatECDown", FatJet_msoftdrop_jesRelativeStatECDown, &b_FatJet_msoftdrop_jesRelativeStatECDown);
   T1->SetBranchAddress("FatJet_pt_jesRelativeStatHFDown", FatJet_pt_jesRelativeStatHFDown, &b_FatJet_pt_jesRelativeStatHFDown);
   T1->SetBranchAddress("FatJet_mass_jesRelativeStatHFDown", FatJet_mass_jesRelativeStatHFDown, &b_FatJet_mass_jesRelativeStatHFDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesRelativeStatHFDown", FatJet_msoftdrop_jesRelativeStatHFDown, &b_FatJet_msoftdrop_jesRelativeStatHFDown);
   T1->SetBranchAddress("FatJet_pt_jesPileUpDataMCDown", FatJet_pt_jesPileUpDataMCDown, &b_FatJet_pt_jesPileUpDataMCDown);
   T1->SetBranchAddress("FatJet_mass_jesPileUpDataMCDown", FatJet_mass_jesPileUpDataMCDown, &b_FatJet_mass_jesPileUpDataMCDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesPileUpDataMCDown", FatJet_msoftdrop_jesPileUpDataMCDown, &b_FatJet_msoftdrop_jesPileUpDataMCDown);
   T1->SetBranchAddress("FatJet_pt_jesPileUpPtRefDown", FatJet_pt_jesPileUpPtRefDown, &b_FatJet_pt_jesPileUpPtRefDown);
   T1->SetBranchAddress("FatJet_mass_jesPileUpPtRefDown", FatJet_mass_jesPileUpPtRefDown, &b_FatJet_mass_jesPileUpPtRefDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtRefDown", FatJet_msoftdrop_jesPileUpPtRefDown, &b_FatJet_msoftdrop_jesPileUpPtRefDown);
   T1->SetBranchAddress("FatJet_pt_jesPileUpPtBBDown", FatJet_pt_jesPileUpPtBBDown, &b_FatJet_pt_jesPileUpPtBBDown);
   T1->SetBranchAddress("FatJet_mass_jesPileUpPtBBDown", FatJet_mass_jesPileUpPtBBDown, &b_FatJet_mass_jesPileUpPtBBDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtBBDown", FatJet_msoftdrop_jesPileUpPtBBDown, &b_FatJet_msoftdrop_jesPileUpPtBBDown);
   T1->SetBranchAddress("FatJet_pt_jesPileUpPtEC1Down", FatJet_pt_jesPileUpPtEC1Down, &b_FatJet_pt_jesPileUpPtEC1Down);
   T1->SetBranchAddress("FatJet_mass_jesPileUpPtEC1Down", FatJet_mass_jesPileUpPtEC1Down, &b_FatJet_mass_jesPileUpPtEC1Down);
   T1->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtEC1Down", FatJet_msoftdrop_jesPileUpPtEC1Down, &b_FatJet_msoftdrop_jesPileUpPtEC1Down);
   T1->SetBranchAddress("FatJet_pt_jesPileUpPtEC2Down", FatJet_pt_jesPileUpPtEC2Down, &b_FatJet_pt_jesPileUpPtEC2Down);
   T1->SetBranchAddress("FatJet_mass_jesPileUpPtEC2Down", FatJet_mass_jesPileUpPtEC2Down, &b_FatJet_mass_jesPileUpPtEC2Down);
   T1->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtEC2Down", FatJet_msoftdrop_jesPileUpPtEC2Down, &b_FatJet_msoftdrop_jesPileUpPtEC2Down);
   T1->SetBranchAddress("FatJet_pt_jesPileUpPtHFDown", FatJet_pt_jesPileUpPtHFDown, &b_FatJet_pt_jesPileUpPtHFDown);
   T1->SetBranchAddress("FatJet_mass_jesPileUpPtHFDown", FatJet_mass_jesPileUpPtHFDown, &b_FatJet_mass_jesPileUpPtHFDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtHFDown", FatJet_msoftdrop_jesPileUpPtHFDown, &b_FatJet_msoftdrop_jesPileUpPtHFDown);
   T1->SetBranchAddress("FatJet_pt_jesPileUpMuZeroDown", FatJet_pt_jesPileUpMuZeroDown, &b_FatJet_pt_jesPileUpMuZeroDown);
   T1->SetBranchAddress("FatJet_mass_jesPileUpMuZeroDown", FatJet_mass_jesPileUpMuZeroDown, &b_FatJet_mass_jesPileUpMuZeroDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesPileUpMuZeroDown", FatJet_msoftdrop_jesPileUpMuZeroDown, &b_FatJet_msoftdrop_jesPileUpMuZeroDown);
   T1->SetBranchAddress("FatJet_pt_jesPileUpEnvelopeDown", FatJet_pt_jesPileUpEnvelopeDown, &b_FatJet_pt_jesPileUpEnvelopeDown);
   T1->SetBranchAddress("FatJet_mass_jesPileUpEnvelopeDown", FatJet_mass_jesPileUpEnvelopeDown, &b_FatJet_mass_jesPileUpEnvelopeDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesPileUpEnvelopeDown", FatJet_msoftdrop_jesPileUpEnvelopeDown, &b_FatJet_msoftdrop_jesPileUpEnvelopeDown);
   T1->SetBranchAddress("FatJet_pt_jesSubTotalPileUpDown", FatJet_pt_jesSubTotalPileUpDown, &b_FatJet_pt_jesSubTotalPileUpDown);
   T1->SetBranchAddress("FatJet_mass_jesSubTotalPileUpDown", FatJet_mass_jesSubTotalPileUpDown, &b_FatJet_mass_jesSubTotalPileUpDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesSubTotalPileUpDown", FatJet_msoftdrop_jesSubTotalPileUpDown, &b_FatJet_msoftdrop_jesSubTotalPileUpDown);
   T1->SetBranchAddress("FatJet_pt_jesSubTotalRelativeDown", FatJet_pt_jesSubTotalRelativeDown, &b_FatJet_pt_jesSubTotalRelativeDown);
   T1->SetBranchAddress("FatJet_mass_jesSubTotalRelativeDown", FatJet_mass_jesSubTotalRelativeDown, &b_FatJet_mass_jesSubTotalRelativeDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesSubTotalRelativeDown", FatJet_msoftdrop_jesSubTotalRelativeDown, &b_FatJet_msoftdrop_jesSubTotalRelativeDown);
   T1->SetBranchAddress("FatJet_pt_jesSubTotalPtDown", FatJet_pt_jesSubTotalPtDown, &b_FatJet_pt_jesSubTotalPtDown);
   T1->SetBranchAddress("FatJet_mass_jesSubTotalPtDown", FatJet_mass_jesSubTotalPtDown, &b_FatJet_mass_jesSubTotalPtDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesSubTotalPtDown", FatJet_msoftdrop_jesSubTotalPtDown, &b_FatJet_msoftdrop_jesSubTotalPtDown);
   T1->SetBranchAddress("FatJet_pt_jesSubTotalScaleDown", FatJet_pt_jesSubTotalScaleDown, &b_FatJet_pt_jesSubTotalScaleDown);
   T1->SetBranchAddress("FatJet_mass_jesSubTotalScaleDown", FatJet_mass_jesSubTotalScaleDown, &b_FatJet_mass_jesSubTotalScaleDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesSubTotalScaleDown", FatJet_msoftdrop_jesSubTotalScaleDown, &b_FatJet_msoftdrop_jesSubTotalScaleDown);
   T1->SetBranchAddress("FatJet_pt_jesSubTotalAbsoluteDown", FatJet_pt_jesSubTotalAbsoluteDown, &b_FatJet_pt_jesSubTotalAbsoluteDown);
   T1->SetBranchAddress("FatJet_mass_jesSubTotalAbsoluteDown", FatJet_mass_jesSubTotalAbsoluteDown, &b_FatJet_mass_jesSubTotalAbsoluteDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesSubTotalAbsoluteDown", FatJet_msoftdrop_jesSubTotalAbsoluteDown, &b_FatJet_msoftdrop_jesSubTotalAbsoluteDown);
   T1->SetBranchAddress("FatJet_pt_jesSubTotalMCDown", FatJet_pt_jesSubTotalMCDown, &b_FatJet_pt_jesSubTotalMCDown);
   T1->SetBranchAddress("FatJet_mass_jesSubTotalMCDown", FatJet_mass_jesSubTotalMCDown, &b_FatJet_mass_jesSubTotalMCDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesSubTotalMCDown", FatJet_msoftdrop_jesSubTotalMCDown, &b_FatJet_msoftdrop_jesSubTotalMCDown);
   T1->SetBranchAddress("FatJet_pt_jesTotalDown", FatJet_pt_jesTotalDown, &b_FatJet_pt_jesTotalDown);
   T1->SetBranchAddress("FatJet_mass_jesTotalDown", FatJet_mass_jesTotalDown, &b_FatJet_mass_jesTotalDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesTotalDown", FatJet_msoftdrop_jesTotalDown, &b_FatJet_msoftdrop_jesTotalDown);
   T1->SetBranchAddress("FatJet_pt_jesTotalNoFlavorDown", FatJet_pt_jesTotalNoFlavorDown, &b_FatJet_pt_jesTotalNoFlavorDown);
   T1->SetBranchAddress("FatJet_mass_jesTotalNoFlavorDown", FatJet_mass_jesTotalNoFlavorDown, &b_FatJet_mass_jesTotalNoFlavorDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesTotalNoFlavorDown", FatJet_msoftdrop_jesTotalNoFlavorDown, &b_FatJet_msoftdrop_jesTotalNoFlavorDown);
   T1->SetBranchAddress("FatJet_pt_jesTotalNoTimeDown", FatJet_pt_jesTotalNoTimeDown, &b_FatJet_pt_jesTotalNoTimeDown);
   T1->SetBranchAddress("FatJet_mass_jesTotalNoTimeDown", FatJet_mass_jesTotalNoTimeDown, &b_FatJet_mass_jesTotalNoTimeDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesTotalNoTimeDown", FatJet_msoftdrop_jesTotalNoTimeDown, &b_FatJet_msoftdrop_jesTotalNoTimeDown);
   T1->SetBranchAddress("FatJet_pt_jesTotalNoFlavorNoTimeDown", FatJet_pt_jesTotalNoFlavorNoTimeDown, &b_FatJet_pt_jesTotalNoFlavorNoTimeDown);
   T1->SetBranchAddress("FatJet_mass_jesTotalNoFlavorNoTimeDown", FatJet_mass_jesTotalNoFlavorNoTimeDown, &b_FatJet_mass_jesTotalNoFlavorNoTimeDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesTotalNoFlavorNoTimeDown", FatJet_msoftdrop_jesTotalNoFlavorNoTimeDown, &b_FatJet_msoftdrop_jesTotalNoFlavorNoTimeDown);
   T1->SetBranchAddress("FatJet_pt_jesFlavorZJetDown", FatJet_pt_jesFlavorZJetDown, &b_FatJet_pt_jesFlavorZJetDown);
   T1->SetBranchAddress("FatJet_mass_jesFlavorZJetDown", FatJet_mass_jesFlavorZJetDown, &b_FatJet_mass_jesFlavorZJetDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesFlavorZJetDown", FatJet_msoftdrop_jesFlavorZJetDown, &b_FatJet_msoftdrop_jesFlavorZJetDown);
   T1->SetBranchAddress("FatJet_pt_jesFlavorPhotonJetDown", FatJet_pt_jesFlavorPhotonJetDown, &b_FatJet_pt_jesFlavorPhotonJetDown);
   T1->SetBranchAddress("FatJet_mass_jesFlavorPhotonJetDown", FatJet_mass_jesFlavorPhotonJetDown, &b_FatJet_mass_jesFlavorPhotonJetDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesFlavorPhotonJetDown", FatJet_msoftdrop_jesFlavorPhotonJetDown, &b_FatJet_msoftdrop_jesFlavorPhotonJetDown);
   T1->SetBranchAddress("FatJet_pt_jesFlavorPureGluonDown", FatJet_pt_jesFlavorPureGluonDown, &b_FatJet_pt_jesFlavorPureGluonDown);
   T1->SetBranchAddress("FatJet_mass_jesFlavorPureGluonDown", FatJet_mass_jesFlavorPureGluonDown, &b_FatJet_mass_jesFlavorPureGluonDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureGluonDown", FatJet_msoftdrop_jesFlavorPureGluonDown, &b_FatJet_msoftdrop_jesFlavorPureGluonDown);
   T1->SetBranchAddress("FatJet_pt_jesFlavorPureQuarkDown", FatJet_pt_jesFlavorPureQuarkDown, &b_FatJet_pt_jesFlavorPureQuarkDown);
   T1->SetBranchAddress("FatJet_mass_jesFlavorPureQuarkDown", FatJet_mass_jesFlavorPureQuarkDown, &b_FatJet_mass_jesFlavorPureQuarkDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureQuarkDown", FatJet_msoftdrop_jesFlavorPureQuarkDown, &b_FatJet_msoftdrop_jesFlavorPureQuarkDown);
   T1->SetBranchAddress("FatJet_pt_jesFlavorPureCharmDown", FatJet_pt_jesFlavorPureCharmDown, &b_FatJet_pt_jesFlavorPureCharmDown);
   T1->SetBranchAddress("FatJet_mass_jesFlavorPureCharmDown", FatJet_mass_jesFlavorPureCharmDown, &b_FatJet_mass_jesFlavorPureCharmDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureCharmDown", FatJet_msoftdrop_jesFlavorPureCharmDown, &b_FatJet_msoftdrop_jesFlavorPureCharmDown);
   T1->SetBranchAddress("FatJet_pt_jesFlavorPureBottomDown", FatJet_pt_jesFlavorPureBottomDown, &b_FatJet_pt_jesFlavorPureBottomDown);
   T1->SetBranchAddress("FatJet_mass_jesFlavorPureBottomDown", FatJet_mass_jesFlavorPureBottomDown, &b_FatJet_mass_jesFlavorPureBottomDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureBottomDown", FatJet_msoftdrop_jesFlavorPureBottomDown, &b_FatJet_msoftdrop_jesFlavorPureBottomDown);
   T1->SetBranchAddress("FatJet_pt_jesTimeRunBDown", FatJet_pt_jesTimeRunBDown, &b_FatJet_pt_jesTimeRunBDown);
   T1->SetBranchAddress("FatJet_mass_jesTimeRunBDown", FatJet_mass_jesTimeRunBDown, &b_FatJet_mass_jesTimeRunBDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesTimeRunBDown", FatJet_msoftdrop_jesTimeRunBDown, &b_FatJet_msoftdrop_jesTimeRunBDown);
   T1->SetBranchAddress("FatJet_pt_jesTimeRunCDown", FatJet_pt_jesTimeRunCDown, &b_FatJet_pt_jesTimeRunCDown);
   T1->SetBranchAddress("FatJet_mass_jesTimeRunCDown", FatJet_mass_jesTimeRunCDown, &b_FatJet_mass_jesTimeRunCDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesTimeRunCDown", FatJet_msoftdrop_jesTimeRunCDown, &b_FatJet_msoftdrop_jesTimeRunCDown);
   T1->SetBranchAddress("FatJet_pt_jesTimeRunDEDown", FatJet_pt_jesTimeRunDEDown, &b_FatJet_pt_jesTimeRunDEDown);
   T1->SetBranchAddress("FatJet_mass_jesTimeRunDEDown", FatJet_mass_jesTimeRunDEDown, &b_FatJet_mass_jesTimeRunDEDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesTimeRunDEDown", FatJet_msoftdrop_jesTimeRunDEDown, &b_FatJet_msoftdrop_jesTimeRunDEDown);
   T1->SetBranchAddress("FatJet_pt_jesTimeRunFDown", FatJet_pt_jesTimeRunFDown, &b_FatJet_pt_jesTimeRunFDown);
   T1->SetBranchAddress("FatJet_mass_jesTimeRunFDown", FatJet_mass_jesTimeRunFDown, &b_FatJet_mass_jesTimeRunFDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesTimeRunFDown", FatJet_msoftdrop_jesTimeRunFDown, &b_FatJet_msoftdrop_jesTimeRunFDown);
   T1->SetBranchAddress("FatJet_pt_jesCorrelationGroupMPFInSituDown", FatJet_pt_jesCorrelationGroupMPFInSituDown, &b_FatJet_pt_jesCorrelationGroupMPFInSituDown);
   T1->SetBranchAddress("FatJet_mass_jesCorrelationGroupMPFInSituDown", FatJet_mass_jesCorrelationGroupMPFInSituDown, &b_FatJet_mass_jesCorrelationGroupMPFInSituDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupMPFInSituDown", FatJet_msoftdrop_jesCorrelationGroupMPFInSituDown, &b_FatJet_msoftdrop_jesCorrelationGroupMPFInSituDown);
   T1->SetBranchAddress("FatJet_pt_jesCorrelationGroupIntercalibrationDown", FatJet_pt_jesCorrelationGroupIntercalibrationDown, &b_FatJet_pt_jesCorrelationGroupIntercalibrationDown);
   T1->SetBranchAddress("FatJet_mass_jesCorrelationGroupIntercalibrationDown", FatJet_mass_jesCorrelationGroupIntercalibrationDown, &b_FatJet_mass_jesCorrelationGroupIntercalibrationDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupIntercalibrationDown", FatJet_msoftdrop_jesCorrelationGroupIntercalibrationDown, &b_FatJet_msoftdrop_jesCorrelationGroupIntercalibrationDown);
   T1->SetBranchAddress("FatJet_pt_jesCorrelationGroupbJESDown", FatJet_pt_jesCorrelationGroupbJESDown, &b_FatJet_pt_jesCorrelationGroupbJESDown);
   T1->SetBranchAddress("FatJet_mass_jesCorrelationGroupbJESDown", FatJet_mass_jesCorrelationGroupbJESDown, &b_FatJet_mass_jesCorrelationGroupbJESDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupbJESDown", FatJet_msoftdrop_jesCorrelationGroupbJESDown, &b_FatJet_msoftdrop_jesCorrelationGroupbJESDown);
   T1->SetBranchAddress("FatJet_pt_jesCorrelationGroupFlavorDown", FatJet_pt_jesCorrelationGroupFlavorDown, &b_FatJet_pt_jesCorrelationGroupFlavorDown);
   T1->SetBranchAddress("FatJet_mass_jesCorrelationGroupFlavorDown", FatJet_mass_jesCorrelationGroupFlavorDown, &b_FatJet_mass_jesCorrelationGroupFlavorDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupFlavorDown", FatJet_msoftdrop_jesCorrelationGroupFlavorDown, &b_FatJet_msoftdrop_jesCorrelationGroupFlavorDown);
   T1->SetBranchAddress("FatJet_pt_jesCorrelationGroupUncorrelatedDown", FatJet_pt_jesCorrelationGroupUncorrelatedDown, &b_FatJet_pt_jesCorrelationGroupUncorrelatedDown);
   T1->SetBranchAddress("FatJet_mass_jesCorrelationGroupUncorrelatedDown", FatJet_mass_jesCorrelationGroupUncorrelatedDown, &b_FatJet_mass_jesCorrelationGroupUncorrelatedDown);
   T1->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupUncorrelatedDown", FatJet_msoftdrop_jesCorrelationGroupUncorrelatedDown, &b_FatJet_msoftdrop_jesCorrelationGroupUncorrelatedDown);
   T1->SetBranchAddress("Jet_deepflavbtagSF", Jet_deepflavbtagSF, &b_Jet_deepflavbtagSF);
   T1->SetBranchAddress("Jet_deepflavbtagSF_up", Jet_deepflavbtagSF_up, &b_Jet_deepflavbtagSF_up);
   T1->SetBranchAddress("Jet_deepflavbtagSF_down", Jet_deepflavbtagSF_down, &b_Jet_deepflavbtagSF_down);
   T1->SetBranchAddress("Jet_deepflavbtagSF_shape", Jet_deepflavbtagSF_shape, &b_Jet_deepflavbtagSF_shape);
   T1->SetBranchAddress("Jet_deepflavbtagSF_shape_up_jes", Jet_deepflavbtagSF_shape_up_jes, &b_Jet_deepflavbtagSF_shape_up_jes);
   T1->SetBranchAddress("Jet_deepflavbtagSF_shape_down_jes", Jet_deepflavbtagSF_shape_down_jes, &b_Jet_deepflavbtagSF_shape_down_jes);
   T1->SetBranchAddress("Jet_deepflavbtagSF_shape_up_lf", Jet_deepflavbtagSF_shape_up_lf, &b_Jet_deepflavbtagSF_shape_up_lf);
   T1->SetBranchAddress("Jet_deepflavbtagSF_shape_down_lf", Jet_deepflavbtagSF_shape_down_lf, &b_Jet_deepflavbtagSF_shape_down_lf);
   T1->SetBranchAddress("Jet_deepflavbtagSF_shape_up_hf", Jet_deepflavbtagSF_shape_up_hf, &b_Jet_deepflavbtagSF_shape_up_hf);
   T1->SetBranchAddress("Jet_deepflavbtagSF_shape_down_hf", Jet_deepflavbtagSF_shape_down_hf, &b_Jet_deepflavbtagSF_shape_down_hf);
   T1->SetBranchAddress("Jet_deepflavbtagSF_shape_up_hfstats1", Jet_deepflavbtagSF_shape_up_hfstats1, &b_Jet_deepflavbtagSF_shape_up_hfstats1);
   T1->SetBranchAddress("Jet_deepflavbtagSF_shape_down_hfstats1", Jet_deepflavbtagSF_shape_down_hfstats1, &b_Jet_deepflavbtagSF_shape_down_hfstats1);
   T1->SetBranchAddress("Jet_deepflavbtagSF_shape_up_hfstats2", Jet_deepflavbtagSF_shape_up_hfstats2, &b_Jet_deepflavbtagSF_shape_up_hfstats2);
   T1->SetBranchAddress("Jet_deepflavbtagSF_shape_down_hfstats2", Jet_deepflavbtagSF_shape_down_hfstats2, &b_Jet_deepflavbtagSF_shape_down_hfstats2);
   T1->SetBranchAddress("Jet_deepflavbtagSF_shape_up_lfstats1", Jet_deepflavbtagSF_shape_up_lfstats1, &b_Jet_deepflavbtagSF_shape_up_lfstats1);
   T1->SetBranchAddress("Jet_deepflavbtagSF_shape_down_lfstats1", Jet_deepflavbtagSF_shape_down_lfstats1, &b_Jet_deepflavbtagSF_shape_down_lfstats1);
   T1->SetBranchAddress("Jet_deepflavbtagSF_shape_up_lfstats2", Jet_deepflavbtagSF_shape_up_lfstats2, &b_Jet_deepflavbtagSF_shape_up_lfstats2);
   T1->SetBranchAddress("Jet_deepflavbtagSF_shape_down_lfstats2", Jet_deepflavbtagSF_shape_down_lfstats2, &b_Jet_deepflavbtagSF_shape_down_lfstats2);
   T1->SetBranchAddress("Jet_deepflavbtagSF_shape_up_cferr1", Jet_deepflavbtagSF_shape_up_cferr1, &b_Jet_deepflavbtagSF_shape_up_cferr1);
   T1->SetBranchAddress("Jet_deepflavbtagSF_shape_down_cferr1", Jet_deepflavbtagSF_shape_down_cferr1, &b_Jet_deepflavbtagSF_shape_down_cferr1);
   T1->SetBranchAddress("Jet_deepflavbtagSF_shape_up_cferr2", Jet_deepflavbtagSF_shape_up_cferr2, &b_Jet_deepflavbtagSF_shape_up_cferr2);
   T1->SetBranchAddress("Jet_deepflavbtagSF_shape_down_cferr2", Jet_deepflavbtagSF_shape_down_cferr2, &b_Jet_deepflavbtagSF_shape_down_cferr2);
   T1->SetBranchAddress("puWeight", &puWeight, &b_puWeight);
   T1->SetBranchAddress("puWeightUp", &puWeightUp, &b_puWeightUp);
   T1->SetBranchAddress("puWeightDown", &puWeightDown, &b_puWeightDown);
   T1->SetBranchAddress("PrefireWeight", &PrefireWeight, &b_PrefireWeight);
   T1->SetBranchAddress("PrefireWeight_Up", &PrefireWeight_Up, &b_PrefireWeight_Up);
   T1->SetBranchAddress("PrefireWeight_Down", &PrefireWeight_Down, &b_PrefireWeight_Down);

// value assignement //

   int nentries = T1->GetEntries();
   cout<<"nentries "<<nentries<<endl;
   
  for (int ij=0; ij<nentries; ij++) {
  if(ij%100==0){ cout<<"event "<<ij+1<<endl; }
   fileIn->cd();

   T1->GetEntry(ij);

   fillarray(isMC);
 
   	
//cout<<"NPDF "<<LHAPDF::numberPDF()+1<<endl;

// event analysis //
     
int topjetAK8 = -1;
int topjetAK8_preli = -1;
int bjetAK4 = -1;
int bjetAK4_preli = -1;
bool bjetAK4_btagged = false;
int top_wp = -1;
int dtop_wp = -1;
int mdtop_wp = -1;
int bjet_wp = -1;
int subjetb_wp = -1;
         
weight = 1.0;
puWeight = 1.0;
   
if(isMC){
	   
      weight = Generator_weight;
      nevent_total += 1;
      weightev = Generator_weight;
      
      #ifdef LHAPDFX
				// npdfmax = 1 + LHAPDF::numberPDF();
		ncpdf = 0;
		
		for (int jk=0; jk< TMath::Min(npdfmax,1 + LHAPDF::numberPDF()) ;jk++) {
				
				LHAPDF::initPDF(jk);
				
				float ypdf1 = LHAPDF::xfx(Generator_x1,Generator_scalePDF, Generator_id1);
				float ypdf2 = LHAPDF::xfx(Generator_x2,Generator_scalePDF, Generator_id2);
				weight_gen_pdf[ncpdf] = weightev*ypdf1*ypdf2;
				ncpdf++;

				if(jk==0){
					ncscale = 0;
					weight_gen_scale[ncscale] = LHAPDF::xfx(Generator_x1,Generator_scalePDF, Generator_id1) * LHAPDF::xfx(Generator_x2,Generator_scalePDF, Generator_id2) * weightev;
					ncscale++;
					weight_gen_scale[ncscale] = LHAPDF::xfx(Generator_x1,0.5*Generator_scalePDF, Generator_id1) * LHAPDF::xfx(Generator_x2,Generator_scalePDF, Generator_id2) * weightev;
					ncscale++;
					weight_gen_scale[ncscale] = LHAPDF::xfx(Generator_x1,Generator_scalePDF, Generator_id1) * LHAPDF::xfx(Generator_x2,0.5*Generator_scalePDF, Generator_id2) * weightev;
					ncscale++;
					weight_gen_scale[ncscale] = LHAPDF::xfx(Generator_x1,2.0*Generator_scalePDF, Generator_id1) * LHAPDF::xfx(Generator_x2,Generator_scalePDF, Generator_id2) * weightev;
					ncscale++;
					weight_gen_scale[ncscale] = LHAPDF::xfx(Generator_x1,Generator_scalePDF, Generator_id1) * LHAPDF::xfx(Generator_x2,2.0*Generator_scalePDF, Generator_id2) * weightev;
					ncscale++;
					weight_gen_scale[ncscale] = LHAPDF::xfx(Generator_x1,0.5*Generator_scalePDF, Generator_id1) * LHAPDF::xfx(Generator_x2,0.5*Generator_scalePDF, Generator_id2) * weightev;
					ncscale++;
					weight_gen_scale[ncscale] = LHAPDF::xfx(Generator_x1,2.0*Generator_scalePDF, Generator_id1) * LHAPDF::xfx(Generator_x2,2.0*Generator_scalePDF, Generator_id2) * weightev;
					ncscale++;
					weight_gen_scale[ncscale] = LHAPDF::xfx(Generator_x1,0.5*Generator_scalePDF, Generator_id1) * LHAPDF::xfx(Generator_x2,2.0*Generator_scalePDF, Generator_id2) * weightev;
					ncscale++;
					weight_gen_scale[ncscale] = LHAPDF::xfx(Generator_x1,2.0*Generator_scalePDF, Generator_id1) * LHAPDF::xfx(Generator_x2,0.5*Generator_scalePDF, Generator_id2) * weightev;
					ncscale++;
					}

		}
	
	
	 #endif
      
      hist_npu_nopuwt->Fill(Pileup_nTrueInt,weight);
     
     #ifdef Anal_2017
	  if(Pileup_nTrueInt>=0 && Pileup_nTrueInt<100){
		  puWeight = pu_rat17[int(Pileup_nTrueInt)];
		  puWeightUp = pu_rat17_up[int(Pileup_nTrueInt)];
		  puWeightDown = pu_rat17_dn[int(Pileup_nTrueInt)];
	  }
	  #endif
	  
	  #ifdef Anal_2016
	  if(Pileup_nTrueInt>=0 && Pileup_nTrueInt<100){
		  puWeight = pu_rat16[int(Pileup_nTrueInt)];
		  puWeightUp = pu_rat16_up[int(Pileup_nTrueInt)];
		  puWeightDown = pu_rat16_dn[int(Pileup_nTrueInt)];
	  }
	  #endif
	  
	  #ifdef Anal_2018
	  if(Pileup_nTrueInt>=0 && Pileup_nTrueInt<100){
		  puWeight = pu_rat18[int(Pileup_nTrueInt)];
		  puWeightUp = pu_rat18_up[int(Pileup_nTrueInt)];
		  puWeightDown = pu_rat18_dn[int(Pileup_nTrueInt)];
	  }
	  #endif
 
	weight *= puWeight;
	hist_npu->Fill(Pileup_nTrueInt,weight);
	
}else{
	   nevent_total += 1;
       weightev = 1;
	}  
   
if(isMC){
	#if defined(Anal_2016) || defined(Anal_2017) 
	weight *= PrefireWeight;
	#endif
}
    
#ifdef Anal_2018
    
PrefireWeight = 1;
PrefireWeight_Up = 1;
PrefireWeight_Down = 1;
    
#endif
    
if(!isMC){
	PrefireWeight = 1;
	PrefireWeight_Up = 1;
	PrefireWeight_Down = 1;	
}
   
Tout->Fill();
   
if(isMC){
	   
   int igentop = -1;
   int igenb = -1;
	   
   double maxgmass = -100;
	   
   for(unsigned ijet=0; ijet<nGenJetAK8; ijet++){
		
		hist_pt_GEN_ak8->Fill(GenJetAK8_pt[ijet],Generator_weight);
		hist_eta_GEN_ak8->Fill(GenJetAK8_eta[ijet],Generator_weight);
			
		if(ijet==0){
				
			hist_pt_GEN_ak8_lead->Fill(GenJetAK8_pt[ijet],Generator_weight);
			hist_eta_GEN_ak8_lead->Fill(GenJetAK8_eta[ijet],Generator_weight);
			
			}
				
		if(GenJetAK8_mass[ijet] > maxgmass) { 	
			igentop = ijet;
			maxgmass = GenJetAK8_mass[ijet];
		}
			
	}
		
	for(unsigned ijet=0; ijet<nGenJet; ijet++){
		
		hist_pt_GEN_ak4->Fill(GenJet_pt[ijet],Generator_weight);
		hist_eta_GEN_ak4->Fill(GenJet_eta[ijet],Generator_weight);
			
		if(ijet==0){
				
			hist_pt_GEN_ak4_lead->Fill(GenJet_pt[ijet],Generator_weight);
			hist_eta_GEN_ak4_lead->Fill(GenJet_eta[ijet],Generator_weight);
				
			}
			
	}
		
	if(igentop>=0){
		
		for(unsigned ijet=0; ijet<nGenJet; ijet++){
			if(delta2R(GenJetAK8_eta[igentop],GenJetAK8_phi[igentop],GenJet_eta[ijet],GenJet_phi[ijet]) > 1.2){
				if (fabs(PhiInRange(GenJetAK8_phi[igentop] - GenJet_phi[ijet])) > 0.5*M_PI){
					igenb = ijet;
					break;
					}
			}
		}
	
		if(igenb>=0){
			
			TLorentzVector tak8_gen; 
			tak8_gen.SetPtEtaPhiM(GenJetAK8_pt[igentop],GenJetAK8_eta[igentop],GenJetAK8_phi[igentop],GenJetAK8_mass[igentop]);
			
			TLorentzVector bak4_gen; 
			bak4_gen.SetPtEtaPhiM(GenJet_pt[igenb],GenJet_eta[igenb],GenJet_phi[igenb],GenJet_mass[igenb]);
			
			if(tak8_gen.Pt()>500 && bak4_gen.Pt()>500){
				hist_mtb_GEN->Fill((tak8_gen+bak4_gen).M(),Generator_weight);
			}
			
		}
		
	  }
		
	}
   
   hist_npv->Fill(PV_npvsGood,weight);
   if(isMC && puWeight>1.e-9){
	hist_npv_nopuwt->Fill(PV_npvsGood,weight*1./puWeight);
	}
   
   if(PV_npvsGood < 1) continue;
   
   bool metFilterpassed = false;
   
   #if defined(Anal_2017) || defined(Anal_2018)
   
   if(isMC){
   metFilterpassed = 	    (	 Flag_goodVertices && Flag_HBHENoiseFilter 
							  && Flag_HBHENoiseIsoFilter
							  && Flag_EcalDeadCellTriggerPrimitiveFilter 
							  && Flag_globalSuperTightHalo2016Filter 
							  && Flag_ecalBadCalibFilterV2
							  && Flag_BadPFMuonFilter );
		//					  && Flag_BadChargedCandidateFilter);
	}else{
		metFilterpassed = 	 (	 Flag_goodVertices && Flag_HBHENoiseFilter 
							  && Flag_HBHENoiseIsoFilter
							  && Flag_EcalDeadCellTriggerPrimitiveFilter 
							  && Flag_globalSuperTightHalo2016Filter 
							  && Flag_ecalBadCalibFilterV2
							  && Flag_BadPFMuonFilter 
							  && Flag_eeBadScFilter    );
		//					  && Flag_BadChargedCandidateFilter);
		 }
	
   #endif
   
   #ifdef Anal_2016
   
   if(isMC){
   metFilterpassed = 	    (	 Flag_goodVertices && Flag_HBHENoiseFilter 
							  && Flag_HBHENoiseIsoFilter
							  && Flag_EcalDeadCellTriggerPrimitiveFilter 
							  && Flag_globalSuperTightHalo2016Filter 
							  && Flag_BadPFMuonFilter );
		//					  && Flag_BadChargedCandidateFilter);
	}else{
		metFilterpassed = 	 (	 Flag_goodVertices && Flag_HBHENoiseFilter 
							  && Flag_HBHENoiseIsoFilter
							  && Flag_EcalDeadCellTriggerPrimitiveFilter 
							  && Flag_globalSuperTightHalo2016Filter 
							  && Flag_BadPFMuonFilter 
							  && Flag_eeBadScFilter    );
		//					  && Flag_BadChargedCandidateFilter);
		 }
	
   #endif
      
   hist_metfilter_pass->Fill(metFilterpassed,weight);
      
   if(!metFilterpassed) continue;
   
   bool mu_trig = (HLT_IsoMu27 || HLT_Mu50);
  
   #if defined(Anal_2017) || defined(Anal_2018)
   had_trig = (HLT_AK8PFJet420_TrimMass30 || HLT_AK8PFHT900_TrimMass50 || HLT_PFJet500 || HLT_AK8PFJet500 || HLT_PFHT1050);
   #endif
   
   #ifdef Data_2017B
    had_trig = (HLT_PFJet500 || HLT_AK8PFJet500 || HLT_PFHT1050);
   #endif
   
   #ifdef Anal_2016
	   had_trig = (HLT_AK8DiPFJet280_200_TrimMass30 || 
                    HLT_AK8DiPFJet300_200_TrimMass30 ||
					HLT_AK8PFHT700_TrimR0p1PT0p03Mass50 || HLT_AK8PFJet360_TrimMass30 ||
					HLT_PFJet450 || HLT_AK8PFJet450 || HLT_PFHT800 || HLT_PFHT900);
 					
   #endif
   
   #ifdef Data_2016H
		had_trig = (HLT_AK8DiPFJet300_200_TrimMass30 || HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20 ||
					HLT_AK8PFHT700_TrimR0p1PT0p03Mass50 || HLT_AK8PFJet360_TrimMass30 ||
					HLT_PFJet450 || HLT_AK8PFJet450 || HLT_PFHT900);
   #endif					
   
   bool el_trig = false;
   
   #if defined(Anal_2017) || defined(Anal_2018)
   el_trig = (HLT_Ele40_WPTight_Gsf);
   #endif
   
   #ifdef Anal_2016
   el_trig = (HLT_Ele27_WPTight_Gsf);
   #endif
    
   hist_mutrig_pass->Fill(mu_trig,weight);
   hist_etrig_pass->Fill(el_trig,weight);
   hist_hadtrig_pass->Fill(had_trig,weight);
   
   bool jetcutAK8 = false;
   jetcutAK8 = (nFatJet>1 && FatJet_pt[0]>250);

   bool jetcutAK4 = false;
   jetcutAK4 = (nJet>1 && Jet_pt[0]>250);
   
   bool dijetcutAK4 = false;
   dijetcutAK4 = (nJet>1 && Jet_pt[0]>250 && Jet_pt[1]>250 && fabs(PhiInRange(Jet_phi[0]-Jet_phi[1])) > 0.5*M_PI);

   bool dijetcutAK8 = false;
   dijetcutAK8 = (nFatJet>1 && FatJet_pt[0]>250 && FatJet_pt[1]>250 && fabs(PhiInRange(FatJet_phi[0]-FatJet_phi[1])) > 0.5*M_PI);
 
   if(nFatJet>0){
      
      float maxscore_preli = -1000;
    
      for(unsigned ijet=0; ijet<nFatJet; ijet++){
          
        TLorentzVector AK8vec4;
        AK8vec4.SetPtEtaPhiM(FatJet_pt[ijet],FatJet_eta[ijet],FatJet_phi[ijet],FatJet_mass[ijet]);

		if(FatJet_pt[ijet] > AK8ptcut_in && fabs(FatJet_eta[ijet])<=jeteta_cut){
		  
			if(FatJet_deepTagMD_TvsQCD[ijet] > maxscore_preli){
				maxscore_preli = FatJet_deepTagMD_TvsQCD[ijet];
				topjetAK8_preli = ijet;
			}
		}
		  
	  }//ijet
  }//nFatJet>0


 bool b_found_preli = false;  	  
    
 if(topjetAK8_preli>=0){
      
  for(unsigned ijet=0; ijet<nJet; ijet++){
			  
	double dR = delta2R(FatJet_eta[topjetAK8_preli],FatJet_phi[topjetAK8_preli],Jet_eta[ijet],Jet_phi[ijet]);
	double dEta = (Jet_eta[ijet]-FatJet_eta[topjetAK8_preli]);
	double dPhi = PhiInRange(Jet_phi[ijet]-FatJet_phi[topjetAK8_preli]);
			  
	if(Jet_MatchFatJet[ijet] == topjetAK8_preli) continue;
			  
	if(dR>dR_cut){
		if(fabs(dPhi) > dphi_cut){	
					
			TLorentzVector AK4vec4; 
			AK4vec4.SetPtEtaPhiM(Jet_pt[ijet],Jet_eta[ijet],Jet_phi[ijet],Jet_mass[ijet]);
					
			if(fabs(Jet_eta[ijet])<=jeteta_cut && (Jet_pt[ijet]>AK4ptcut_in)){
						
			if(!b_found_preli){
						
				bjetAK4_preli = ijet;
				b_found_preli = true;
				break;
					
			  }
		    }
		}// dPhi
	 } // dR
  }//ijet		  
		  
}//topjetAK8

   // trigger plots fill //
   
   if((/*el_trig||*/mu_trig) && jetcutAK8 && jetcutAK4 && dijetcutAK4 && dijetcutAK8 && weight>=0){
	   
	   if(topjetAK8_preli>=0 && bjetAK4_preli>=0){
		   
		   TLorentzVector pt; pt.SetPtEtaPhiM(FatJet_pt[topjetAK8_preli],FatJet_eta[topjetAK8_preli],FatJet_phi[topjetAK8_preli],FatJet_mass[topjetAK8_preli]);
		   TLorentzVector pb; pb.SetPtEtaPhiM(Jet_pt[bjetAK4_preli],Jet_eta[bjetAK4_preli],Jet_phi[bjetAK4_preli],Jet_mass[bjetAK4_preli]);
		   TLorentzVector ptb; ptb = pt+pb;
		   
		   float twoptsum = (FatJet_pt[topjetAK8_preli]+Jet_pt[bjetAK4_preli]);
		   
		   int imsd = getbinid(FatJet_msoftdrop[topjetAK8_preli],trig_msdbins,trig_msd_vals);
		   
		   if(imsd>=0){
		   			
			hist_ht_emutrig[imsd]->Fill(Event_Ht,weight);
			hist_pt_emutrig[imsd]->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_emutrig[imsd]->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_emutrig[imsd]->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_emutrig[imsd]->Fill(ptb.M(),weight);
			}
		
		#if defined(Anal_2017) || defined(Anal_2018)
			
		if(HLT_PFHT1050 ){
			hist_ht_HT1050_wemutrig[imsd]->Fill(Event_Ht,weight);
			hist_pt_HT1050_wemutrig[imsd]->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_HT1050_wemutrig[imsd]->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_HT1050_wemutrig[imsd]->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_HT1050_wemutrig[imsd]->Fill(ptb.M(),weight);
				}
			}	
		if(HLT_PFJet500){
			hist_ht_AK4Pt500_wemutrig[imsd]->Fill(Event_Ht,weight);
			hist_pt_AK4Pt500_wemutrig[imsd]->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK4Pt500_wemutrig[imsd]->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK4Pt500_wemutrig[imsd]->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK4Pt500_wemutrig[imsd]->Fill(ptb.M(),weight);
				}
			}
		if(HLT_AK8PFJet500){
			hist_ht_AK8Pt500_wemutrig[imsd]->Fill(Event_Ht,weight);
			hist_pt_AK8Pt500_wemutrig[imsd]->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK8Pt500_wemutrig[imsd]->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK8Pt500_wemutrig[imsd]->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK8Pt500_wemutrig[imsd]->Fill(ptb.M(),weight);
				}
			}
		if(HLT_AK8PFJet420_TrimMass30){
			hist_ht_AK8PFJet420_TrimMass30_wemutrig[imsd]->Fill(Event_Ht,weight);
			hist_pt_AK8PFJet420_TrimMass30_wemutrig[imsd]->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK8PFJet420_TrimMass30_wemutrig[imsd]->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK8PFJet420_TrimMass30_wemutrig[imsd]->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK8PFJet420_TrimMass30_wemutrig[imsd]->Fill(ptb.M(),weight);
				}
			}
		
		if(HLT_AK8PFHT900_TrimMass50){
			hist_ht_AK8PFHT900_TrimMass50_wemutrig[imsd]->Fill(Event_Ht,weight);
			hist_pt_AK8PFHT900_TrimMass50_wemutrig[imsd]->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK8PFHT900_TrimMass50_wemutrig[imsd]->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK8PFHT900_TrimMass50_wemutrig[imsd]->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK8PFHT900_TrimMass50_wemutrig[imsd]->Fill(ptb.M(),weight);
				}
			}
		
		#endif
		
		#ifdef Anal_2016
		
		#ifdef Data_2016H
			
		if(HLT_PFHT900){
			hist_ht_HT1050_wemutrig[imsd]->Fill(Event_Ht,weight);
			hist_pt_HT1050_wemutrig[imsd]->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_HT1050_wemutrig[imsd]->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_HT1050_wemutrig[imsd]->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_HT1050_wemutrig[imsd]->Fill(ptb.M(),weight);
				}
			}	
		
		#else
		
		if(HLT_PFHT800||HLT_PFHT900){
			hist_ht_HT1050_wemutrig[imsd]->Fill(Event_Ht,weight);
			hist_pt_HT1050_wemutrig[imsd]->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_HT1050_wemutrig[imsd]->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_HT1050_wemutrig[imsd]->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_HT1050_wemutrig[imsd]->Fill(ptb.M(),weight);
				}
			}	
		
		#endif	
			
		if(HLT_PFJet450){
			hist_ht_AK4Pt500_wemutrig[imsd]->Fill(Event_Ht,weight);
			hist_pt_AK4Pt500_wemutrig[imsd]->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK4Pt500_wemutrig[imsd]->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK4Pt500_wemutrig[imsd]->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK4Pt500_wemutrig[imsd]->Fill(ptb.M(),weight);
				}
			}
		if(HLT_AK8PFJet450){
			hist_ht_AK8Pt500_wemutrig[imsd]->Fill(Event_Ht,weight);
			hist_pt_AK8Pt500_wemutrig[imsd]->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK8Pt500_wemutrig[imsd]->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK8Pt500_wemutrig[imsd]->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK8Pt500_wemutrig[imsd]->Fill(ptb.M(),weight);
				}
			}
		if(HLT_AK8PFJet360_TrimMass30){
			hist_ht_AK8PFJet420_TrimMass30_wemutrig[imsd]->Fill(Event_Ht,weight);
			hist_pt_AK8PFJet420_TrimMass30_wemutrig[imsd]->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK8PFJet420_TrimMass30_wemutrig[imsd]->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK8PFJet420_TrimMass30_wemutrig[imsd]->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK8PFJet420_TrimMass30_wemutrig[imsd]->Fill(ptb.M(),weight);
				}
			}
		
		if(HLT_AK8PFHT700_TrimR0p1PT0p03Mass50){
			hist_ht_AK8PFHT900_TrimMass50_wemutrig[imsd]->Fill(Event_Ht,weight);
			hist_pt_AK8PFHT900_TrimMass50_wemutrig[imsd]->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK8PFHT900_TrimMass50_wemutrig[imsd]->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK8PFHT900_TrimMass50_wemutrig[imsd]->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK8PFHT900_TrimMass50_wemutrig[imsd]->Fill(ptb.M(),weight);
				}
			}
		
		#ifdef Data_2016H
		
		if( HLT_AK8DiPFJet300_200_TrimMass30 || HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20) {
			hist_ht_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]->Fill(Event_Ht,weight);
			hist_pt_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]->Fill(ptb.M(),weight);
				}
			}
			
		#else	
		
		if( HLT_AK8DiPFJet280_200_TrimMass30 || HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20 ||
		    HLT_AK8DiPFJet300_200_TrimMass30 || HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20) {
			hist_ht_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]->Fill(Event_Ht,weight);
			hist_pt_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK8DiPFJet300_200_TrimMass30_wemutrig[imsd]->Fill(ptb.M(),weight);
				}
			}
		
		#endif
		
		#endif
					
		if(had_trig){
			hist_ht_all_wemutrig[imsd]->Fill(Event_Ht,weight);
			hist_pt_all_wemutrig[imsd]->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_all_wemutrig[imsd]->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_all_wemutrig[imsd]->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_all_wemutrig[imsd]->Fill(ptb.M(),weight);
				}
			}
		}	
	 }
   }
   // trigger fill end //
     
   hist_nmuons->Fill(nMuon,weight);
   hist_nelectrons->Fill(nElectron,weight);
   hist_nphotons->Fill(nPhoton,weight);
   
   if(had_trig){
	   hist_nmuons_trig_pass->Fill(nMuon,weight);
	   }
      
   if((nMuon  > 0) && (Muon_pt[0]>lepptcut)) continue;
   if((nElectron  > 0) && (Electron_pt[0]>lepptcut)) continue;
   if((nPhoton  > 0) && (Photon_pt[0]>lepptcut)) continue;
    
   if(!(had_trig)) continue; 
   
   if(nFatJet<1 || nJet<2) continue;
   
   if(isMC){
	weightpass = weight;
	Tout1->Fill();
   }

// pt fill to calculate b-tagging efficiency //
   
   if(isMC){
   
   for(unsigned ijet=0; ijet<nJet; ijet++){
	
	int ieta = getbinid(fabs(Jet_eta[ijet]),netarange,betarange);
	
	if(ieta>=0){
	
	if(Jet_genJetIdx[ijet]>=0 && Jet_genJetIdx[ijet]<njetmax){
		if (GenJet_pt[Jet_genJetIdx[ijet]]>8.){
	
	if(abs(Jet_hadronFlavour[ijet])==5){
		hist_btAK4_pt[ieta]->Fill(Jet_pt[ijet],weight);
		if(Jet_btagCSVV2[ijet]>btagvalue){
		hist_bAK4_pt[ieta]->Fill(Jet_pt[ijet],weight);
		}
		if(Jet_btagDeepB[ijet]>btagvalue_deepCSV){
		hist_dbAK4_pt[ieta]->Fill(Jet_pt[ijet],weight);
		}
		if(Jet_btagDeepFlavB[ijet]>btagvalue_deepFlavB){
		hist_dfbAK4_pt[ieta]->Fill(Jet_pt[ijet],weight);
		}
		
		int ipt = getbinid(Jet_pt[ijet],noperf_ptbins,perf_ptbins);
		if(ipt >= 0 && ipt<(noperf_ptbins)){
		hist_btag_csv_bhadron[ipt]->Fill(Jet_btagCSVV2[ijet],weight);
		hist_btag_deepcsv_bhadron[ipt]->Fill(Jet_btagDeepB[ijet],weight);
		hist_btag_deepflav_bhadron[ipt]->Fill(Jet_btagDeepFlavB[ijet],weight);
		}
		hist_btag_csv_bhadron[noperf_ptbins]->Fill(Jet_btagCSVV2[ijet],weight);
		hist_btag_deepcsv_bhadron[noperf_ptbins]->Fill(Jet_btagDeepB[ijet],weight);
		hist_btag_deepflav_bhadron[noperf_ptbins]->Fill(Jet_btagDeepFlavB[ijet],weight);
	}
    
   if(abs(Jet_hadronFlavour[ijet])==0){
		hist_qtAK4_pt[ieta]->Fill(Jet_pt[ijet],weight);
		if(Jet_btagCSVV2[ijet]>btagvalue){
		hist_qAK4_pt[ieta]->Fill(Jet_pt[ijet],weight);
		}
		if(Jet_btagDeepB[ijet]>btagvalue_deepCSV){
		hist_dqAK4_pt[ieta]->Fill(Jet_pt[ijet],weight);
		}
		if(Jet_btagDeepFlavB[ijet]>btagvalue_deepFlavB){
		hist_dfqAK4_pt[ieta]->Fill(Jet_pt[ijet],weight);
		}
		int ipt = getbinid(Jet_pt[ijet],noperf_ptbins,perf_ptbins);
		if(ipt >= 0 && ipt<(noperf_ptbins)){
		hist_btag_csv_qhadron[ipt]->Fill(Jet_btagCSVV2[ijet],weight);
		hist_btag_deepcsv_qhadron[ipt]->Fill(Jet_btagDeepB[ijet],weight);
		hist_btag_deepflav_qhadron[ipt]->Fill(Jet_btagDeepFlavB[ijet],weight);
			}
		hist_btag_csv_qhadron[noperf_ptbins]->Fill(Jet_btagCSVV2[ijet],weight);
		hist_btag_deepcsv_qhadron[noperf_ptbins]->Fill(Jet_btagDeepB[ijet],weight);
		hist_btag_deepflav_qhadron[noperf_ptbins]->Fill(Jet_btagDeepFlavB[ijet],weight);			
		}
		
	  } //gen pt	
     } //gen match
	}//ieta
   }//ijet
 }//isMC
   // end //
      
// if(Event_Ht<HTcut) return kFALSE; 
    
   if(isMC){
	   
	 // shape based discriminator //
	 
	#ifndef NoJEC
		float mcprod = 1.;
		float dataprod = 1.;
	    btag_weight = 1.;
		for(int ibtag=0; ibtag<nbtagmax; ibtag++){
			btag_weight_unc_up[ibtag] = 1;
			btag_weight_unc_dn[ibtag] = 1;
		}
	 
	 #ifdef Btagger_DeepJet_wt
	 UInt_t nbcandjetAK4 = TMath::Min(int(nJet),int(4));
	  for(unsigned ijet=0; ijet<nbcandjetAK4; ijet++){
			
			if(!(fabs(Jet_eta[ijet]) <= jeteta_cut && (Jet_pt[ijet]>AK4ptcut_fi && Jet_pt[ijet]<10000.))) continue;
			btag_weight *= Jet_deepflavbtagSF_shape[ijet];
	
			for(int ibtag=0; ibtag<nbtagmax; ibtag++){
				btag_weight_unc_up[ibtag] *= Jet_deepflavbtagSF_unc_up[ibtag][ijet];
				btag_weight_unc_dn[ibtag] *= Jet_deepflavbtagSF_unc_dn[ibtag][ijet];
			}
		}
	
		weight *= 	btag_weight;
		
	 #endif	
/*	 
	 #ifdef Btagger_DeepCSV_wt
	  for(unsigned ijet=0; ijet<nJet; ijet++){
			
			if(!(fabs(Jet_eta[ijet]) <=jeteta_cut && (Jet_pt[ijet]>30. && Jet_pt[ijet]<10000.))) continue;
			btag_weight *= Jet_deepCSVbtagSF_shape[ijet];
		}
	
		weight *= 	btag_weight;
	 #endif	
	 
	 #ifdef Btagger_CSVv2
	  for(unsigned ijet=0; ijet<nJet; ijet++){
		
			if(!(fabs(Jet_eta[ijet]) <=jeteta_cut && (Jet_pt[ijet]>30. && Jet_pt[ijet]<10000.))) continue;
			btag_weight *= Jet_CSVbtagSF_shape[ijet];
		}
	
		weight *= 	btag_weight;
	 #endif	
*/	 	
	#endif

	weightpass_btag = weight;
	Tout1->Fill();
	  }
	   
   #ifdef NoJEC   
   met_corpt->Fill(MET_pt,weight);
   met_corphi->Fill(MET_phi,weight);
   met_bysumEt->Fill(MET_pt*1./MET_sumEt,weight);
   #else
   met_corpt->Fill(MET_pt_nom,weight);
   met_corphi->Fill(MET_phi_nom,weight);
   met_bysumEt->Fill(MET_pt_nom*1./MET_sumEt,weight);
   #endif
   
      
   int nbjetsAK4 = 0;
   
   #ifdef Btagger_DeepJet
   
   for(unsigned ijet=0; ijet<nJet; ijet++){
			 
	      if(fabs(Jet_eta[ijet]) <= jeteta_cut && (Jet_pt[ijet]>30.)){
			if(Jet_btagDeepFlavB[ijet] > btagvalue_deepFlavB){
			  nbjetsAK4++;
			}
		 }
	  }
	  
  #endif	
	
  #ifdef Btagger_DeepCSV
  
  for(unsigned ijet=0; ijet<nJet; ijet++){
			 
	      if(fabs(Jet_eta[ijet]) <=jeteta_cut && (Jet_pt[ijet]>30.)){
			if(Jet_btagDeepB[ijet] > btagvalue_deepCSV){
			  nbjetsAK4++;
			}
		 }
	  }
  
  #endif	
  
  #ifdef Btagger_CSVv2
  
  for(unsigned ijet=0; ijet<nJet; ijet++){
			 
	      if(fabs(Jet_eta[ijet]) <=jeteta_cut && (Jet_pt[ijet]>30.)){
			if(Jet_btagCSVV2[ijet] > btagvalue){
			  nbjetsAK4++;
			}
		 }
	  }
  
  #endif
	 
	  hist_nbjetAK4->Fill(nbjetsAK4,weight);
		
	  hist_njetAK8->Fill(nFatJet,weight);
	  hist_njetAK4->Fill(nJet,weight);

// quarks from top //

int top_d[6]={-1,-1,-1,-1,-1,-1};
int ist=0; int isqg = 0;
int qp[4] = {-1,-1,-1,-1};
int bp[2] = {-1,-1};
int igtop[2] = {-1,-1};
int ngtop = 0; int nghtop = 0;
int ihtop[2] = {-1,-1};
int qg_d[6]={-1,-1,-1,-1,-1,-1};
TLorentzVector topparton_vec4;
TLorentzVector topdaugh_vec4[2][3];

weight_t = weight;

if(isMC){	
	for(unsigned ijet=0; ijet<nFatJet; ijet++){
		
		hist_jetpt_all->Fill(FatJet_pt[ijet],weight);
		
		if((FatJet_msoftdrop[ijet] > topmasslow) && (FatJet_msoftdrop[ijet] < topmasshigh)){
			hist_jetpt_topmsd_all->Fill(FatJet_pt[ijet],weight);
			}
		
		if(FatJet_deepTagMD_TvsQCD[ijet] > deepak8_cut_md){
			hist_jetpt_all_DAK8pass->Fill(FatJet_pt[ijet],weight);
			if((FatJet_msoftdrop[ijet] > topmasslow) && (FatJet_msoftdrop[ijet] < topmasshigh)){
				hist_jetpt_topmsd_all_DAK8pass->Fill(FatJet_pt[ijet],weight);
			}
		}
	}
}

if(isMC){

for(unsigned igen=0; igen<nGenPart; igen++){ 
				  
	if(!(GenPart_status[igen]==22)) continue;
	if(!(abs(GenPart_pdgId[igen])==6)) continue;
//	int ishard = ((decToBinary(GenPart_statusFlags[igen]))/10000000)%10;
	int ishard = *(decToBinary(GenPart_statusFlags[igen])+7);
	if(ishard != 1) continue;
			 
	igtop[ngtop] = int(igen);
	ngtop++;
	if(ngtop==2) break;		  
	
	}

}

if(isMC && TopTagging){

for(unsigned igen=0; igen<nGenPart; igen++){
	if(GenPart_status[igen]==0) continue;
			  
		int isprompt = *(decToBinary(GenPart_statusFlags[igen])+0);//((decToBinary(GenPart_statusFlags[igen]))%10);
		int ishard = *(decToBinary(GenPart_statusFlags[igen])+7);//((decToBinary(GenPart_statusFlags[igen]))/10000000)%10;
		int isFirstcopy = *(decToBinary(GenPart_statusFlags[igen])+12);
		int fromhard = *(decToBinary(GenPart_statusFlags[igen])+8);//((decToBinary(GenPart_statusFlags[igen]))/100000000);
		int isLastcopy = *(decToBinary(GenPart_statusFlags[igen])+13);
		int pid = abs(GenPart_pdgId[igen]);
		
		if(ngtop>0){  
			if((isprompt==1) && (fromhard==1) && (abs(GenPart_status[igen])==23) && (GenPart_genPartIdxMother[igen] >= 0) && (pid<6)){
				if((abs(GenPart_pdgId[igen])==5 && abs(GenPart_pdgId[GenPart_genPartIdxMother[igen]])==6)||(abs(GenPart_pdgId[igen])<5 && abs(GenPart_pdgId[GenPart_genPartIdxMother[igen]])==24)){
					top_d[ist] = igen;
					ist++;
			  }
		  }
			if(ist==6) break;
	    }
		
		if(ngtop==0 && isSignal){
			if((isprompt==1) && (fromhard==1) && (abs(GenPart_pdgId[igen])<6) && (ishard==0) &&(isLastcopy==1) /*&& (GenPart_statusFlags[igen]==8449)*/){
				top_d[ist] = igen;
				ist++;
				}
				if(ist==6) break;
			}
			
		if(ngtop==0 && !isSignal){
			if((isprompt==1) && (fromhard==1) && (ishard==1) && (isFirstcopy==1) && (abs(GenPart_status[igen])==23) && ((pid>=1 && pid<6)||pid==21) ){
				qg_d[isqg] = igen;
				isqg++;
				}
				if(isqg==6) break;
			}
	
	}


if(!isSignal && ngtop>0){
		  
	if(ist>=3){
		
		int iq = 0; int ib = 0;
		  
		for(int id=0; id<ist; id++){
			if(top_d[id]>=0 && abs(GenPart_pdgId[top_d[id]])<5)  { qp[iq] = top_d[id]; iq++; }
			if(iq>=4) break;
		}
		
		for(int id=0; id<ist; id++){
			if(top_d[id]>=0 && abs(GenPart_pdgId[top_d[id]])==5) { bp[ib] = top_d[id]; ib++; }
			if(ib>=2) break;
		}
		  
		if(iq>=3){
			 top_d[0] = qp[0];
			 if(GenPart_genPartIdxMother[qp[0]] == GenPart_genPartIdxMother[qp[1]]) { top_d[1] = qp[1];  top_d[3] = qp[2];  top_d[4] = qp[3];}
			 if(GenPart_genPartIdxMother[qp[0]] == GenPart_genPartIdxMother[qp[2]]) { top_d[1] = qp[2];  top_d[3] = qp[3];  top_d[4] = qp[1];}
			 if(GenPart_genPartIdxMother[qp[0]] == GenPart_genPartIdxMother[qp[3]]) { top_d[1] = qp[3];  top_d[3] = qp[1];  top_d[4] = qp[2];}
		 }
		
		if(ist==3){  
			top_d[2] = bp[0]; top_d[3] = -1; top_d[4] = -1; top_d[5] = -1;
			
			TLorentzVector p1; p1.SetPtEtaPhiM(GenPart_pt[top_d[0]],GenPart_eta[top_d[0]],GenPart_phi[top_d[0]],GenPart_mass[top_d[0]]);
			TLorentzVector p2; p2.SetPtEtaPhiM(GenPart_pt[top_d[1]],GenPart_eta[top_d[1]],GenPart_phi[top_d[1]],GenPart_mass[top_d[1]]);
			TLorentzVector p3; p3.SetPtEtaPhiM(GenPart_pt[top_d[2]],GenPart_eta[top_d[2]],GenPart_phi[top_d[2]],GenPart_mass[top_d[2]]);
	  
			topdaugh_vec4[0][0] = p1;
			topdaugh_vec4[0][1] = p2;
			topdaugh_vec4[0][2] = p3;
		}
		
		if(ist>=4){
		 
		if(GenPart_pdgId[GenPart_genPartIdxMother[top_d[0]]]*GenPart_pdgId[bp[0]]>0) { top_d[2] = bp[0];  top_d[5] = bp[1]; }
		  else { 
			  if(GenPart_pdgId[GenPart_genPartIdxMother[top_d[0]]]*GenPart_pdgId[bp[1]]>0) { top_d[2] = bp[1];  top_d[5] = bp[0];}
			   else { top_d[2] = -1; }
			   }
		  
		  TLorentzVector p1; p1.SetPtEtaPhiM(GenPart_pt[top_d[0]],GenPart_eta[top_d[0]],GenPart_phi[top_d[0]],GenPart_mass[top_d[0]]);
		  TLorentzVector p2; p2.SetPtEtaPhiM(GenPart_pt[top_d[1]],GenPart_eta[top_d[1]],GenPart_phi[top_d[1]],GenPart_mass[top_d[1]]);
		  TLorentzVector p3; p3.SetPtEtaPhiM(GenPart_pt[top_d[2]],GenPart_eta[top_d[2]],GenPart_phi[top_d[2]],GenPart_mass[top_d[2]]);
	  
		  TLorentzVector p4; p4.SetPtEtaPhiM(GenPart_pt[top_d[3]],GenPart_eta[top_d[3]],GenPart_phi[top_d[3]],GenPart_mass[top_d[3]]);
		  TLorentzVector p5; p5.SetPtEtaPhiM(GenPart_pt[top_d[4]],GenPart_eta[top_d[4]],GenPart_phi[top_d[4]],GenPart_mass[top_d[4]]);
		  TLorentzVector p6; p6.SetPtEtaPhiM(GenPart_pt[top_d[5]],GenPart_eta[top_d[5]],GenPart_phi[top_d[5]],GenPart_mass[top_d[5]]);
	  
	  
		  topdaugh_vec4[0][0] = p1;
		  topdaugh_vec4[0][1] = p2;
		  topdaugh_vec4[0][2] = p3;
		  
		  topdaugh_vec4[1][0] = p4;
		  topdaugh_vec4[1][1] = p5;
		  topdaugh_vec4[1][2] = p6;
		  
		  }	
	}


  if(ngtop > 1 && igtop[0]>=0 && igtop[1]>=0 && ist>=3){
  	
	if(top_d[0]>=0){
		
		if(abs(GenPart_pdgId[top_d[0]])<5) { nghtop++; ihtop[0] = igtop[0]; }
		if(top_d[3]>=0 && abs(GenPart_pdgId[top_d[3]])<5) { nghtop++;  ihtop[1] = igtop[1]; }
		
	}else { nghtop = 0; }
	
 }
 
 if(ngtop==1 && top_d[0]>=0 && abs(GenPart_pdgId[top_d[0]])<5){
	 ihtop[0] = igtop[0];
	 nghtop = 1;
	 }
	 
 if(ngtop<1){
	 nghtop = 0;
	 }	 

 }// !isignal


if(!isSignal && nghtop>0 && ihtop[0]>=0){	
	topparton_vec4.SetPtEtaPhiM(GenPart_pt[ihtop[0]],GenPart_eta[ihtop[0]],GenPart_phi[ihtop[0]],GenPart_mass[ihtop[0]]);
	}

if(isSignal && ist==4){
	
	int qs[2] = {-1,-1};
	int bs[2] = {-1,-1};
	
	int iq=0; int ib=0;
	
	for(int id=0; id<ist; id++){
			if(abs(GenPart_pdgId[top_d[id]])<5)  { qs[iq] = top_d[id]; iq++; }
			if(iq>=2) break;
		}
		
		for(int id=0; id<ist; id++){
			if(abs(GenPart_pdgId[top_d[id]])==5) { bs[ib] = top_d[id]; ib++; }
			if(ib>=2) break;
		}
	
	  top_d[0] = (abs(GenPart_pdgId[qs[0]])>abs(GenPart_pdgId[qs[1]]))?qs[0]:qs[1]; 
	  top_d[1] = (abs(GenPart_pdgId[qs[0]])<abs(GenPart_pdgId[qs[1]]))?qs[0]:qs[1]; 
	  
	  if(GenPart_pdgId[top_d[0]] > 0){
		  top_d[2] = (GenPart_pdgId[bs[0]]>0) ? bs[0] : bs[1];
		  top_d[3] = (GenPart_pdgId[bs[0]]<0) ? bs[0] : bs[1];
		  }
	  
	  if(GenPart_pdgId[top_d[0]] < 0){
		  top_d[2] = (GenPart_pdgId[bs[0]]<0) ? bs[0] : bs[1];
		  top_d[3] = (GenPart_pdgId[bs[0]]>0) ? bs[0] : bs[1];
		  }
	  
//	  str = TString::Format("p1 p2 p3  %i %i %i\n",GenPart_pdgId[top_d[0]],GenPart_pdgId[top_d[1]],GenPart_pdgId[top_d[2]]);
//    if(gProofServ) gProofServ->SendAsynMessage(str);
	  
	  TLorentzVector p1; p1.SetPtEtaPhiM(GenPart_pt[top_d[0]],GenPart_eta[top_d[0]],GenPart_phi[top_d[0]],GenPart_mass[top_d[0]]);
	  TLorentzVector p2; p2.SetPtEtaPhiM(GenPart_pt[top_d[1]],GenPart_eta[top_d[1]],GenPart_phi[top_d[1]],GenPart_mass[top_d[1]]);
	  TLorentzVector p3; p3.SetPtEtaPhiM(GenPart_pt[top_d[2]],GenPart_eta[top_d[2]],GenPart_phi[top_d[2]],GenPart_mass[top_d[2]]);
	  
	  topdaugh_vec4[0][0] = p1;
	  topdaugh_vec4[0][1] = p2;
	  topdaugh_vec4[0][2] = p3;
	  
	  topparton_vec4 = p1+p2+p3;
	}

}//isMC
// end of quarks from top //

// top tagging study //

int gnjetAK8_top = -1;
int pfjetAK8_top = -1;

if(isMC && nGenJetAK8>0){
	double mindiff = 1000;
	
	for(unsigned igen=0; igen<nGenJetAK8; igen++){
		if(fabs(GenJetAK8_mass[igen] - 173) < mindiff){
			mindiff = fabs(GenJetAK8_mass[igen] - 173);
			gnjetAK8_top = int(igen);
		}
	}
}

if(gnjetAK8_top>=0 && isMC && isSignal){
	
	double mindR = 0.8;
		
	for(unsigned ijet=0; ijet<nFatJet; ijet++){
			
		double dR_pg = delta2R(topparton_vec4.Eta(),topparton_vec4.Phi(),FatJet_eta[ijet],FatJet_phi[ijet]);
		if(dR_pg<mindR){
			mindR = dR_pg;
			pfjetAK8_top = int(ijet);
			}
		}
		hist_mass_hadtopq->Fill(topparton_vec4.M(),weight);
		hist_matchjet_hadtopq->Fill(pfjetAK8_top,weight);
	}

// top pt reweighing  // ttbar only

SF_toppt =  1;

if(ngtop==2 && isTTBar && isMC){
	SF_toppt =  sqrt(exp(0.0615-0.0005*GenPart_pt[igtop[0]]) * exp(0.0615-0.0005*GenPart_pt[igtop[1]]));
//	SF_toppt = SF_TOP(0.0615,0.0005,GenPart_pt[igtop[0]],GenPart_pt[igtop[1]]);
//  SF_toppt =  sqrt(exp(0.0615-0.0005*TMath::Min(float(400.),GenPart_pt[igtop[0]])) * exp(0.0615-0.0005*TMath::Min(float(400.),GenPart_pt[igtop[1]])));
  //*sqrt(EW_toppt_cor(GenPart_pt[igtop[0]])*EW_toppt_cor(GenPart_pt[igtop[1]]));
  weight_t *= SF_toppt;
}

// top pt reweighing ends

// Top tagging scale factors

// tau32 scale factor //

sfwt_tau32 = 1;
	  
if((nghtop>0) && isMC && TopTagging){

for(int ih=0; ih<nghtop; ih++){
	 
 if(ihtop[ih]<0) continue;		 
	 
 for(unsigned ijet=0; ijet<nFatJet; ijet++){
			
  double dR_pg = delta2R(GenPart_eta[ihtop[ih]],GenPart_phi[ihtop[ih]],FatJet_eta[ijet],FatJet_phi[ijet]);	
  
  if(dR_pg>0.6) continue; 
		  	  
  float sfwt = 1;
		  
  int mtag = -1;
		  
  if(ist>=3 && (top_d[0]>=0 && top_d[1]>=0 && top_d[2]>=0)){
			  
  if((delta2R(GenPart_eta[top_d[0]],GenPart_phi[top_d[0]],FatJet_eta[ijet],FatJet_phi[ijet])<0.8)&&(delta2R(GenPart_eta[top_d[1]],GenPart_phi[top_d[1]],FatJet_eta[ijet],FatJet_phi[ijet])<0.8)&&(delta2R(GenPart_eta[top_d[2]],GenPart_phi[top_d[2]],FatJet_eta[ijet],FatJet_phi[ijet])<0.8)){
	   mtag = 0;
	}else{
			  
		if(((delta2R(GenPart_eta[top_d[0]],GenPart_phi[top_d[0]],FatJet_eta[ijet],FatJet_phi[ijet])<0.8)&&(delta2R(GenPart_eta[top_d[1]],GenPart_phi[top_d[1]],FatJet_eta[ijet],FatJet_phi[ijet])<0.8)&&(delta2R(GenPart_eta[top_d[2]],GenPart_phi[top_d[2]],FatJet_eta[ijet],FatJet_phi[ijet])>=0.8))
		  || ((delta2R(GenPart_eta[top_d[1]],GenPart_phi[top_d[1]],FatJet_eta[ijet],FatJet_phi[ijet])<0.8)&&(delta2R(GenPart_eta[top_d[2]],GenPart_phi[top_d[2]],FatJet_eta[ijet],FatJet_phi[ijet])<0.8)&&(delta2R(GenPart_eta[top_d[0]],GenPart_phi[top_d[0]],FatJet_eta[ijet],FatJet_phi[ijet])>=0.8))
		  || ((delta2R(GenPart_eta[top_d[2]],GenPart_phi[top_d[2]],FatJet_eta[ijet],FatJet_phi[ijet])<0.8)&&(delta2R(GenPart_eta[top_d[0]],GenPart_phi[top_d[0]],FatJet_eta[ijet],FatJet_phi[ijet])<0.8)&&(delta2R(GenPart_eta[top_d[1]],GenPart_phi[top_d[1]],FatJet_eta[ijet],FatJet_phi[ijet])>=0.8)) ) {
				  mtag = 1;
		  }
		  else{
			  if(((delta2R(GenPart_eta[top_d[0]],GenPart_phi[top_d[0]],FatJet_eta[ijet],FatJet_phi[ijet])<0.8)&&(delta2R(GenPart_eta[top_d[1]],GenPart_phi[top_d[1]],FatJet_eta[ijet],FatJet_phi[ijet])>=0.8)&&(delta2R(GenPart_eta[top_d[2]],GenPart_phi[top_d[2]],FatJet_eta[ijet],FatJet_phi[ijet])>=0.8))
			   ||((delta2R(GenPart_eta[top_d[1]],GenPart_phi[top_d[1]],FatJet_eta[ijet],FatJet_phi[ijet])<0.8)&&(delta2R(GenPart_eta[top_d[2]],GenPart_phi[top_d[2]],FatJet_eta[ijet],FatJet_phi[ijet])>=0.8)&&(delta2R(GenPart_eta[top_d[0]],GenPart_phi[top_d[0]],FatJet_eta[ijet],FatJet_phi[ijet])>=0.8))
			   ||((delta2R(GenPart_eta[top_d[2]],GenPart_phi[top_d[2]],FatJet_eta[ijet],FatJet_phi[ijet])<0.8)&&(delta2R(GenPart_eta[top_d[0]],GenPart_phi[top_d[0]],FatJet_eta[ijet],FatJet_phi[ijet])>=0.8)&&(delta2R(GenPart_eta[top_d[1]],GenPart_phi[top_d[1]],FatJet_eta[ijet],FatJet_phi[ijet])>=0.8)) ){
				  mtag = 2;
		  }
	   }
	}
	
  }
		  
  TLorentzVector tjet; 
  tjet.SetPtEtaPhiM(FatJet_pt[ijet],FatJet_eta[ijet],FatJet_phi[ijet],FatJet_mass[ijet]);
		  
  if(mtag>=0){
	if(tjet.Pt()<=400.) { sfwt = tau32SF_nom[mtag][0]; }
	if(tjet.Pt()>400. && tjet.Pt()<=500.) { sfwt = tau32SF_nom[mtag][1]; }
	if(tjet.Pt()>500. && tjet.Pt()<=600.) { sfwt = tau32SF_nom[mtag][2]; }
	if(tjet.Pt()>600.) { sfwt = tau32SF_nom[mtag][3]; }
  }
		 
	 sfwt_tau32 *= sfwt ; 
	 
	 break;
	 
	  } //ijet
	}//ih	 
 } //for top scale factor 

sfwt_deepak8 = 1.;
sfwt_deepak8_md = 1.;
sfwt_deepak8_md_up = 1;
sfwt_deepak8_md_dn = 1;

if(isMC){
	
	float ttag_eff = 1;
	if(isQCD){
		ttag_eff = 0.0275;
		}else{
			if(isTTHad){
				ttag_eff = 0.6; // tt hadronic 
			}
			if(isTTSemiLep){
				ttag_eff = 0.54; // tt hadronic 
			}
			if(isTTDiLep){
				ttag_eff = 0.35; // tt hadronic 
			}
			if(isST||isSignal){
				ttag_eff = 0.55; // tt hadronic 
			}
		 }
	
	UInt_t ntcandjAK8 = TMath::Min(int(nFatJet),int(2));
	
	for(unsigned ijet=0; ijet<ntcandjAK8; ijet++){
	    
	    if(FatJet_pt[ijet] < AK8ptcut) continue;
	    if(FatJet_msoftdrop[ijet]<topmasslow || FatJet_msoftdrop[ijet]>topmasshigh) continue;
	    
	    if(isQCD){
			int iptbin_perf = getbinid(FatJet_pt[ijet],noperf_ptbins,perf_ptbins);
			if(iptbin_perf>=0 && iptbin_perf<noperf_ptbins){
			ttag_eff = qcd_md_deepak8_eff[iptbin_perf];
			}
		}
	    
		int iptbin = getbinid(FatJet_pt[ijet],noAK8ptbins,AK8ptbins);
		if(iptbin>=0){
			if(FatJet_deepTag_TvsQCD[ijet] >= deepak8_cut){
				sfwt_deepak8 *= DeepAK8_SF[iptbin];
			}else{
				  sfwt_deepak8 *= (1.-ttag_eff*DeepAK8_SF[iptbin])*1./(1.-ttag_eff);
			  }
			  
			if(FatJet_deepTagMD_TvsQCD[ijet] >= deepak8_cut_md){
				sfwt_deepak8_md *= DeepAK8_MD_SF[iptbin];
				sfwt_deepak8_md_up *= DeepAK8_MD_SF_up[iptbin];
				sfwt_deepak8_md_dn *= DeepAK8_MD_SF_dn[iptbin];
			}else{
				  sfwt_deepak8_md *= (1.-ttag_eff*DeepAK8_MD_SF[iptbin])*1./(1.-ttag_eff);
				  sfwt_deepak8_md_up *= (1.-ttag_eff*DeepAK8_MD_SF_up[iptbin])*1./(1.-ttag_eff);
				  sfwt_deepak8_md_dn *= (1.-ttag_eff*DeepAK8_MD_SF_dn[iptbin])*1./(1.-ttag_eff);
				 }	
		}
	}
	
	weight_t *= sfwt_deepak8_md ; 
}


if(isMC && !isSignal &&(nghtop>0)){
	
	for(unsigned ijet=0; ijet<nFatJet; ijet++){
		FatJet_hashtop[ijet] = false;
		for(int ih=0; ih<nghtop; ih++){
			double dR_pg = delta2R(GenPart_eta[ihtop[ih]],GenPart_phi[ihtop[ih]],FatJet_eta[ijet],FatJet_phi[ijet]);
			if(dR_pg < 0.6){
				FatJet_hashtop[ijet] = true;
				break;
			}
		}
	}
	
	for(unsigned ijet=0; ijet<nFatJet; ijet++){
		if(FatJet_hashtop[ijet]){
			hist_jetpt_hastop->Fill(FatJet_pt[ijet],weight);
			if(FatJet_deepTagMD_TvsQCD[ijet] > deepak8_cut_md){
				hist_jetpt_hastop_DAK8pass->Fill(FatJet_pt[ijet],weight);
			}
		}
	}
}

if(isMC && !isSignal &&(nghtop>0)){
	
	double mindR = 0.8;
	
	for(int ih=0; ih<nghtop; ih++){
	
	for(unsigned ijet=0; ijet<nFatJet; ijet++){
			
		double dR_pg = delta2R(GenPart_eta[ihtop[ih]],GenPart_phi[ihtop[ih]],FatJet_eta[ijet],FatJet_phi[ijet]);
		double dR_tdj = std::max({delta2R(topdaugh_vec4[ih][0].Eta(),topdaugh_vec4[ih][0].Phi(),FatJet_eta[ijet],FatJet_phi[ijet]),
								  delta2R(topdaugh_vec4[ih][1].Eta(),topdaugh_vec4[ih][1].Phi(),FatJet_eta[ijet],FatJet_phi[ijet]),
								  delta2R(topdaugh_vec4[ih][2].Eta(),topdaugh_vec4[ih][2].Phi(),FatJet_eta[ijet],FatJet_phi[ijet])});
								  
        double topmass_part = (topdaugh_vec4[ih][0]+topdaugh_vec4[ih][1]+topdaugh_vec4[ih][2]).M();						  
								  
		if(dR_pg < 0.6 && dR_tdj<0.6 && (topmass_part>165 && topmass_part<180)){
			
			hist_jetpt_top->Fill(FatJet_pt[ijet],weight);
			if(FatJet_msoftdrop[ijet] > topmasslow && FatJet_msoftdrop[ijet] < topmasshigh){
				hist_jetpt_msd_top->Fill(FatJet_pt[ijet],weight);
			}
			if(FatJet_deepTagMD_TvsQCD[ijet] > deepak8_cut_md){
				hist_jetpt_top_md_deepak8_pass->Fill(FatJet_pt[ijet],weight);
				if(FatJet_msoftdrop[ijet] > topmasslow && FatJet_msoftdrop[ijet] < topmasshigh){
					hist_jetpt_msd_top_md_deepak8_pass->Fill(FatJet_pt[ijet],weight);
					}
				}
				
			
			int ipt = getbinid(FatJet_pt[ijet],noperf_ptbins,perf_ptbins);
			if(ipt >= 0 && ipt<(noperf_ptbins)){
			
		 	if(fabs(FatJet_eta[ijet])<2.5){
				hist_jetsdmass_top[ipt]->Fill(FatJet_msoftdrop[ijet],weight_t);
				float tau32 = FatJet_tau3[ijet]*1./TMath::Max(float_t(1.e-6),FatJet_tau2[ijet]);
				hist_jettau32_top[ipt]->Fill(tau32,weight_t);
				hist_jetsubbtag_CSVv2_top[ipt]->Fill(FatJet_subbtagCSVV2[ijet],weight_t);
				hist_jetsubbtag_DeepCSV_top[ipt]->Fill(FatJet_subbtagDeepB[ijet],weight_t);
				hist_jetDeepAK8_MD_top[ipt]->Fill(FatJet_deepTagMD_TvsQCD[ijet],weight_t);
				hist_jetDeepAK8_top[ipt]->Fill(FatJet_deepTag_TvsQCD[ijet],weight_t);
				if(FatJet_msoftdrop[ijet] > topmasslow && FatJet_msoftdrop[ijet] < topmasshigh){
					hist_jetDeepAK8_MD_wm_top[ipt]->Fill(FatJet_deepTagMD_TvsQCD[ijet],weight_t);
					hist_jetDeepAK8_wm_top[ipt]->Fill(FatJet_deepTag_TvsQCD[ijet],weight_t);
					}
				
				hist_jetgenmass_top[ipt]->Fill(FatJet_GenJetmass[ijet],weight_t);
				if(FatJet_msoftdrop[ijet] > topmasslow && FatJet_msoftdrop[ijet] < topmasshigh){
					if(tau32 < tau32_cut_1 && FatJet_subbtagDeepB[ijet]>btagvalue_deepCSV_m){
						hist_jetDeepAK8_top_tausub_11[ipt]->Fill(FatJet_deepTag_TvsQCD[ijet],weight_t);
						}
					if(tau32 < tau32_cut_2 && FatJet_subbtagDeepB[ijet]>btagvalue_deepCSV_m){
						hist_jetDeepAK8_top_tausub_21[ipt]->Fill(FatJet_deepTag_TvsQCD[ijet],weight_t);
						}
					if(tau32 < tau32_cut_1 && FatJet_subbtagDeepB[ijet]>btagvalue_deepCSV_t){
						hist_jetDeepAK8_top_tausub_12[ipt]->Fill(FatJet_deepTag_TvsQCD[ijet],weight_t);
						}
					if(tau32 < tau32_cut_2 && FatJet_subbtagDeepB[ijet]>btagvalue_deepCSV_t){
						hist_jetDeepAK8_top_tausub_22[ipt]->Fill(FatJet_deepTag_TvsQCD[ijet],weight_t);
						}			
					}
				}
			}	
				hist_matchjet_hadtopq->Fill(ijet,weight);
			}
			
			break;
		}
	
		hist_mass_hadtopq->Fill(GenPart_mass[ihtop[ih]],weight);
		hist_pid_hadtopq->Fill(abs(GenPart_pdgId[ihtop[ih]]),weight);
	}	
	
}


if(isMC && isSignal){
	
	for(unsigned ijet=0; ijet<nFatJet; ijet++){
			
		double dR_pg = delta2R(topparton_vec4.Eta(),topparton_vec4.Phi(),FatJet_eta[ijet],FatJet_phi[ijet]);
		double dR_tdj = std::max({delta2R(topdaugh_vec4[0][0].Eta(),topdaugh_vec4[0][0].Phi(),FatJet_eta[ijet],FatJet_phi[ijet]),
								  delta2R(topdaugh_vec4[0][1].Eta(),topdaugh_vec4[0][1].Phi(),FatJet_eta[ijet],FatJet_phi[ijet]),
								  delta2R(topdaugh_vec4[0][2].Eta(),topdaugh_vec4[0][2].Phi(),FatJet_eta[ijet],FatJet_phi[ijet])});
		
		double topmass_part = (topdaugh_vec4[0][0]+topdaugh_vec4[0][1]+topdaugh_vec4[0][2]).M();
		
		if(dR_pg < 0.6 && dR_tdj<0.6 && topmass_part>165 && topmass_part<180){
			
			hist_jetpt_top->Fill(FatJet_pt[ijet],weight);
			if(FatJet_msoftdrop[ijet] > topmasslow && FatJet_msoftdrop[ijet] < topmasshigh){
				hist_jetpt_msd_top->Fill(FatJet_pt[ijet],weight);
			}
			if(FatJet_deepTagMD_TvsQCD[ijet] > deepak8_cut_md){
				hist_jetpt_top_md_deepak8_pass->Fill(FatJet_pt[ijet],weight);
				if(FatJet_msoftdrop[ijet] > topmasslow && FatJet_msoftdrop[ijet] < topmasshigh){
					hist_jetpt_msd_top_md_deepak8_pass->Fill(FatJet_pt[ijet],weight);
					}
				}
			
			int ipt = getbinid(FatJet_pt[ijet],noperf_ptbins,perf_ptbins);
			if(ipt >= 0 && ipt<(noperf_ptbins)){
			
		 	if(fabs(FatJet_eta[ijet])<2.5){
				hist_jetsdmass_top[ipt]->Fill(FatJet_msoftdrop[ijet],weight_t);
				float tau32 = FatJet_tau3[ijet]*1./TMath::Max(float_t(1.e-6),FatJet_tau2[ijet]);
				hist_jettau32_top[ipt]->Fill(tau32,weight_t);
				hist_jetsubbtag_CSVv2_top[ipt]->Fill(FatJet_subbtagCSVV2[ijet],weight_t);
				hist_jetsubbtag_DeepCSV_top[ipt]->Fill(FatJet_subbtagDeepB[ijet],weight_t);
				hist_jetDeepAK8_MD_top[ipt]->Fill(FatJet_deepTagMD_TvsQCD[ijet],weight_t);
				if(FatJet_msoftdrop[ijet] > topmasslow && FatJet_msoftdrop[ijet] < topmasshigh){
					hist_jetDeepAK8_MD_wm_top[ipt]->Fill(FatJet_deepTagMD_TvsQCD[ijet],weight_t);
					hist_jetDeepAK8_wm_top[ipt]->Fill(FatJet_deepTag_TvsQCD[ijet],weight_t);
					}
				hist_jetDeepAK8_top[ipt]->Fill(FatJet_deepTag_TvsQCD[ijet],weight_t);
				hist_jetgenmass_top[ipt]->Fill(FatJet_GenJetmass[ijet],weight_t);
				if(FatJet_msoftdrop[ijet] > topmasslow && FatJet_msoftdrop[ijet] < topmasshigh){
					if(tau32 < tau32_cut_1 && FatJet_subbtagDeepB[ijet]>btagvalue_deepCSV_m){
						hist_jetDeepAK8_top_tausub_11[ipt]->Fill(FatJet_deepTag_TvsQCD[ijet],weight_t);
						}
					if(tau32 < tau32_cut_2 && FatJet_subbtagDeepB[ijet]>btagvalue_deepCSV_m){
						hist_jetDeepAK8_top_tausub_21[ipt]->Fill(FatJet_deepTag_TvsQCD[ijet],weight_t);
						}
					if(tau32 < tau32_cut_1 && FatJet_subbtagDeepB[ijet]>btagvalue_deepCSV_t){
						hist_jetDeepAK8_top_tausub_12[ipt]->Fill(FatJet_deepTag_TvsQCD[ijet],weight_t);
						}
					if(tau32 < tau32_cut_2 && FatJet_subbtagDeepB[ijet]>btagvalue_deepCSV_t){
						hist_jetDeepAK8_top_tausub_22[ipt]->Fill(FatJet_deepTag_TvsQCD[ijet],weight_t);
						}			
					}
				}
			}
				hist_matchjet_hadtopq->Fill(ijet,weight);
			}
			
	//		break;
		}
	
		hist_mass_hadtopq->Fill(topparton_vec4.M(),weight);
			
}


if(isMC && !isSignal){
	
	for(int ih=0; ih<isqg; ih++){
		
	for(unsigned ijet=0; ijet<nFatJet; ijet++){
			
		double dR_pg = delta2R(GenPart_eta[qg_d[ih]],GenPart_phi[qg_d[ih]],FatJet_eta[ijet],FatJet_phi[ijet]);
		if(dR_pg < 0.6){
			
			int ipt = getbinid(FatJet_pt[ijet],noperf_ptbins,perf_ptbins);
			if(ipt >= 0 && ipt<(noperf_ptbins)){
			
		 	if(fabs(FatJet_eta[ijet])<2.5){
				hist_jetsdmass_tbkg[ipt]->Fill(FatJet_msoftdrop[ijet],weight_t);
				float tau32 = FatJet_tau3[ijet]*1./TMath::Max(float_t(1.e-6),FatJet_tau2[ijet]);
				hist_jettau32_tbkg[ipt]->Fill(tau32,weight_t);
				hist_jetsubbtag_CSVv2_tbkg[ipt]->Fill(FatJet_subbtagCSVV2[ijet],weight_t);
				hist_jetsubbtag_DeepCSV_tbkg[ipt]->Fill(FatJet_subbtagDeepB[ijet],weight_t);
				hist_jetDeepAK8_MD_tbkg[ipt]->Fill(FatJet_deepTagMD_TvsQCD[ijet],weight_t);
				if(FatJet_msoftdrop[ijet] > topmasslow && FatJet_msoftdrop[ijet] < topmasshigh){
					hist_jetDeepAK8_MD_wm_tbkg[ipt]->Fill(FatJet_deepTagMD_TvsQCD[ijet],weight_t);
					hist_jetDeepAK8_wm_tbkg[ipt]->Fill(FatJet_deepTag_TvsQCD[ijet],weight_t);
					}
				hist_jetDeepAK8_tbkg[ipt]->Fill(FatJet_deepTag_TvsQCD[ijet],weight_t);
				hist_jetgenmass_tbkg[ipt]->Fill(FatJet_GenJetmass[ijet],weight_t);
				if(FatJet_msoftdrop[ijet] > topmasslow && FatJet_msoftdrop[ijet] < topmasshigh){
					if(tau32 < tau32_cut_1 && FatJet_subbtagDeepB[ijet]>btagvalue_deepCSV_m){
						hist_jetDeepAK8_tbkg_tausub_11[ipt]->Fill(FatJet_deepTag_TvsQCD[ijet],weight_t);
						}
					if(tau32 < tau32_cut_2 && FatJet_subbtagDeepB[ijet]>btagvalue_deepCSV_m){
						hist_jetDeepAK8_tbkg_tausub_21[ipt]->Fill(FatJet_deepTag_TvsQCD[ijet],weight_t);
						}
					if(tau32 < tau32_cut_1 && FatJet_subbtagDeepB[ijet]>btagvalue_deepCSV_t){
						hist_jetDeepAK8_tbkg_tausub_12[ipt]->Fill(FatJet_deepTag_TvsQCD[ijet],weight_t);
						}
					if(tau32 < tau32_cut_2 && FatJet_subbtagDeepB[ijet]>btagvalue_deepCSV_t){
						hist_jetDeepAK8_tbkg_tausub_22[ipt]->Fill(FatJet_deepTag_TvsQCD[ijet],weight_t);
						}			
					}
				}
			}
				hist_matchjet_hadtopq->Fill(ijet,weight);
			}
		
			break;
		}
		
		hist_pid_qg->Fill(abs(GenPart_pdgId[qg_d[ih]]),weight);
	}

}


// end of top tagging study //
	  
//	  weight_t = weight;
	  
	  
if(nFatJet>0){
      
      float maxscore = -1000;
      
      for(unsigned ijet=0; ijet<nFatJet; ijet++){
     
		TLorentzVector AK8vec4; 
		AK8vec4.SetPtEtaPhiM(FatJet_pt[ijet],FatJet_eta[ijet],FatJet_phi[ijet],FatJet_mass[ijet]);	
		hist_jetrapAK8->Fill(AK8vec4.Rapidity(),weight); 
		 
		if(fabs(FatJet_eta[ijet])<=jeteta_cut){
			
			hist_jetptAK8->Fill(FatJet_pt[ijet],weight);
			
			if(FatJet_pt[ijet] > AK8ptcut){
			
				hist_jetmassAK8->Fill(FatJet_mass[ijet],weight);
				hist_jetsdmassAK8->Fill(FatJet_msoftdrop[ijet],weight);
				hist_jetrhoAK8->Fill(-2*log(FatJet_msoftdrop[ijet]/FatJet_pt[ijet]),weight);
				
				if(FatJet_deepTagMD_TvsQCD[ijet] > maxscore){
					maxscore = FatJet_deepTagMD_TvsQCD[ijet];
					topjetAK8 = ijet;
				}
			 }//pt-cut
		  }//y-cut
		   
	  }//ijet
	  
	  hist_topcand_AK8->Fill(topjetAK8,weight);
 }
	  
	  
bool b_found = false;  	  
    
if(topjetAK8>=0){
      
  TLorentzVector tcand4v;
  tcand4v.SetPtEtaPhiM(FatJet_pt[topjetAK8],FatJet_eta[topjetAK8],FatJet_phi[topjetAK8],FatJet_mass[topjetAK8]);
        
  for(unsigned ijet=0; ijet<nJet; ijet++){
	
	TLorentzVector jcand4v;
	jcand4v.SetPtEtaPhiM(Jet_pt[ijet],Jet_eta[ijet],Jet_phi[ijet],Jet_mass[ijet]);
        		  
	double dR = delta2R(FatJet_eta[topjetAK8],FatJet_phi[topjetAK8],Jet_eta[ijet],Jet_phi[ijet]);
	double dEta = (jcand4v.Rapidity()-tcand4v.Rapidity());
	double dPhi = PhiInRange(Jet_phi[ijet]-FatJet_phi[topjetAK8]);
			  
	hist_delR_AK4_toptagAK8->Fill(dR,weight);
	hist_deleta_AK4_toptagAK8->Fill(dEta,weight);
	hist_delphi_AK4_toptagAK8->Fill(dPhi,weight);
		  
	if(Jet_MatchFatJet[ijet] == topjetAK8) continue;
			  
	if(dR>dR_cut){
		if(had_trig){
			hist_delphi_AK4_toptagAK8_dRpass->Fill(dPhi,weight);
			hist_deleta_AK4_toptagAK8_dRpass->Fill(dEta,weight);
		}
		if(fabs(dPhi) > dphi_cut){	
					
			TLorentzVector AK4vec4; 
			AK4vec4.SetPtEtaPhiM(Jet_pt[ijet],Jet_eta[ijet],Jet_phi[ijet],Jet_mass[ijet]);
			    	  
			hist_jetptAK4->Fill(Jet_pt[ijet],weight);
			hist_jetrapAK4->Fill(AK4vec4.Rapidity(),weight);
					
			if(fabs(Jet_eta[ijet]) <= jeteta_cut && (Jet_pt[ijet]>AK4ptcut_in)){
						
			if(!b_found){
						
				bjetAK4 = ijet;
				b_found = true;
				
					hist_jetmassAK4->Fill(AK4vec4.M(),weight);
					hist_jetbtagAK4->Fill(Jet_btagCSVV2[bjetAK4],weight);
					hist_jetbtagdeepCSVAK4->Fill(Jet_btagDeepB[bjetAK4],weight);
					hist_jetbtagdeepflavAK4->Fill(Jet_btagDeepFlavB[bjetAK4],weight);
					
					if(Jet_pt[bjetAK4] > AK4ptcut_fi){
						hist_jetmassAK4_1->Fill(AK4vec4.M(),weight);
						hist_jetbtagAK4_1->Fill(Jet_btagCSVV2[bjetAK4],weight);
						hist_jetbtagdeepCSVAK4_1->Fill(Jet_btagDeepB[bjetAK4],weight);
						hist_jetbtagdeepflavAK4_1->Fill(Jet_btagDeepFlavB[bjetAK4],weight);
					}
			  }
		    }
		}// dPhi
	 } // dR
  }//ijet		  
		  
}//topjetAK8	  
	  
	  if(topjetAK8>=0){
		  
		  hist_topjetpt->Fill(FatJet_pt[topjetAK8],weight_t);
		  hist_topjetsdmass->Fill(FatJet_msoftdrop[topjetAK8],weight_t);
		  hist_topjetdeeptopscore->Fill(FatJet_deepTag_TvsQCD[topjetAK8],weight_t);
		  hist_topjetmddeeptopscore->Fill(FatJet_deepTagMD_TvsQCD[topjetAK8],weight_t);
		  
		  hist_topjetptmsd_2d->Fill(FatJet_msoftdrop[topjetAK8],FatJet_pt[topjetAK8],weight_t);
		  
		  float tau32AK8 = FatJet_tau3[topjetAK8]*1./TMath::Max(float(1.e-6),FatJet_tau2[topjetAK8]);
		  hist_topjettau32AK8->Fill(tau32AK8,weight_t);
		  
		  if(isMC){
			  
			float dR_tj = -1;
			
			if(ihtop[0]>=0 && ihtop[1]>=1) { dR_tj = min(delta2R(GenPart_eta[ihtop[0]],GenPart_phi[ihtop[0]],FatJet_eta[topjetAK8],FatJet_phi[topjetAK8]),delta2R(GenPart_eta[ihtop[1]],GenPart_phi[ihtop[1]],FatJet_eta[topjetAK8],FatJet_phi[topjetAK8]));}
			else{
			  if(ihtop[0]>=0) { dR_tj = delta2R(GenPart_eta[ihtop[0]],GenPart_phi[ihtop[0]],FatJet_eta[topjetAK8],FatJet_phi[topjetAK8]); }
			  if(ihtop[1]>=0) { dR_tj = delta2R(GenPart_eta[ihtop[1]],GenPart_phi[ihtop[1]],FatJet_eta[topjetAK8],FatJet_phi[topjetAK8]); }
			  }
			if(dR_tj >=0 && dR_tj<0.8){
			  hist_topjetptmatch->Fill(FatJet_pt[topjetAK8],weight_t);
			  }  
			
			if(ihtop[0]>=0){  
				int iclose = 0;
				if(ihtop[1]>=0 && (ihtop[0]!=ihtop[1])){
					if(delta2R(GenPart_eta[ihtop[0]],GenPart_phi[ihtop[0]],FatJet_eta[topjetAK8],FatJet_phi[topjetAK8]) > delta2R(GenPart_eta[ihtop[1]],GenPart_phi[ihtop[1]],FatJet_eta[topjetAK8],FatJet_phi[topjetAK8]) ) { iclose = 1; }
				}
				hist_toppartonpt->Fill(GenPart_pt[ihtop[iclose]],weight_t);
				if(delta2R(GenPart_eta[ihtop[iclose]],GenPart_phi[ihtop[iclose]],FatJet_eta[topjetAK8],FatJet_phi[topjetAK8])<0.8){
					hist_topparton_matchedtopjet_pt->Fill(GenPart_pt[ihtop[iclose]],weight);
					}
				}
			
			  
			if(FatJet_MatchGenJet[topjetAK8]>=0){
				hist_topjetgenmass->Fill(FatJet_GenJetmass[topjetAK8],weight_t);
				hist_topjetgenpid->Fill(GenJetAK8_partonFlavour[topjetAK8],weight_t);
			}
		 
		  if(isSignal){
				hist_toppartonpt->Fill(topparton_vec4.Pt(),weight);
				dR_tj = delta2R(topparton_vec4.Eta(),topparton_vec4.Phi(),FatJet_eta[topjetAK8],FatJet_phi[topjetAK8]) ;
				if(dR_tj<0.8){
					hist_topparton_matchedtopjet_pt->Fill(topparton_vec4.Pt(),weight);
					hist_topjetptmatch->Fill(FatJet_pt[topjetAK8],weight);
					}
			  }
		 
		  }
		  
		  float submass[2]={0}; float subpt[2]={0};
		  int isbcnt = 0;
		  for(unsigned isub=0; isub<nSubJet; isub++){	
			if((int(isub)!=FatJet_subJetIdx1[topjetAK8])&&(int(isub)!=FatJet_subJetIdx2[topjetAK8])) continue;
			submass[isbcnt] = SubJet_mass[isub];
			subpt[isbcnt] = SubJet_pt[isub];
			isbcnt++;
			if(isbcnt==2) break;
			}
		
		if(isbcnt==2){
		  (submass[0] > submass[1])? hist_2D_topjet_subjetmass12->Fill(submass[0],submass[1],weight_t) : hist_2D_topjet_subjetmass12->Fill(submass[1],submass[0],weight_t);
		   if(submass[0] > submass[1]){
			  hist_2D_topjetsdmass_subjetmass1->Fill(FatJet_msoftdrop[topjetAK8],submass[0],weight_t);
			  hist_2D_topjetsdmass_subjetmass2->Fill(FatJet_msoftdrop[topjetAK8],submass[1],weight_t);
			}else{
				hist_2D_topjetsdmass_subjetmass1->Fill(FatJet_msoftdrop[topjetAK8],submass[1],weight_t);
				hist_2D_topjetsdmass_subjetmass2->Fill(FatJet_msoftdrop[topjetAK8],submass[0],weight_t);
				}
			}
		  
		  TLorentzVector AK8vec4;
		  AK8vec4.SetPtEtaPhiM(FatJet_pt[topjetAK8],FatJet_eta[topjetAK8],FatJet_phi[topjetAK8],FatJet_mass[topjetAK8]);
		 
		  hist_2D_sdmass_tau32_AK8->Fill(FatJet_msoftdrop[topjetAK8],tau32AK8,weight_t);
		  hist_2D_sdmass_subbtag_AK8->Fill(FatJet_msoftdrop[topjetAK8],FatJet_subbtagDeepB[topjetAK8],weight_t);
		  
		  hist_2D_tau32_deeptopscore->Fill(tau32AK8,FatJet_deepTag_TvsQCD[topjetAK8],weight_t);
		  hist_2D_sdmass_deeptopscore->Fill(FatJet_msoftdrop[topjetAK8],FatJet_deepTag_TvsQCD[topjetAK8],weight_t);
		  hist_2D_mass_deeptopscore->Fill(FatJet_mass[topjetAK8],FatJet_deepTag_TvsQCD[topjetAK8],weight_t);
		  hist_2D_subjetbtag_deeptopscore->Fill(FatJet_subbtagDeepB[topjetAK8],FatJet_deepTag_TvsQCD[topjetAK8],weight_t);
		  
		  hist_2D_tau32_deepmdtopscore->Fill(tau32AK8,FatJet_deepTagMD_TvsQCD[topjetAK8],weight_t);
		  hist_2D_sdmass_deepmdtopscore->Fill(FatJet_msoftdrop[topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight_t);
		  hist_2D_mass_deepmdtopscore->Fill(FatJet_mass[topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight_t);
		  hist_2D_subjetbtag_deepmdtopscore->Fill(FatJet_subbtagDeepB[topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight_t);
		  
		  if(FatJet_msoftdrop[topjetAK8] > topmasslow && FatJet_msoftdrop[topjetAK8] < topmasshigh){
			  
			  hist_jettau32AK8_mpass->Fill(tau32AK8,weight);  
			  
			  float maxsubjetbtagvalue = -100;
			  int sub_b = -1;
			  
			  if(nSubJet > 0){
				  
			   int nbsubjetsAK8 = 0;
				  
			  for(unsigned isub=0; isub<nSubJet; isub++){
				
					if((int(isub)!=FatJet_subJetIdx1[topjetAK8])&&(int(isub)!=FatJet_subJetIdx2[topjetAK8])) continue;
		
					if(int(isub) == FatJet_subbtagId_DeepCSV[topjetAK8]){
						hist_jetsubmassAK8_mpass->Fill(SubJet_mass[isub],weight_t);
						}
						
					hist_jetsubptAK8_mpass->Fill(SubJet_pt[isub],weight);
					
					if(SubJet_btagDeepB[isub] > btagvalue){
						nbsubjetsAK8++;
						}
					
					}
				
				hist_jetsubbtagAK8_mpass->Fill(FatJet_subbtagCSVV2[topjetAK8],weight_t);
				hist_jetsubbtagAK8_deepCSV_mpass->Fill(FatJet_subbtagDeepB[topjetAK8],weight_t);
				
				if(FatJet_subbtagDeepB[topjetAK8] > btagvalue){
					hist_jetsubptAK8_btagpass_mpass->Fill(SubJet_pt[sub_b],weight_t);
				}
				hist_nbsubjetAK8_mpass->Fill(nbsubjetsAK8,weight);
				
				hist_2D_tau32_subbtag_AK8_mpass->Fill(tau32AK8,FatJet_subbtagDeepB[topjetAK8],weight_t);
				
				if((FatJet_subbtagDeepB[topjetAK8]>btagvalue) && (tau32AK8 <tau32_cut)){
			
				  TLorentzVector sub4bvec;
				  sub4bvec.SetPtEtaPhiM(SubJet_pt[sub_b],SubJet_eta[sub_b],SubJet_phi[sub_b],SubJet_mass[sub_b]);
				  hist_jetsubthetAK8->Fill(sub4bvec.Vect()*AK8vec4.Vect()*1./(sub4bvec.Vect().Mag()*AK8vec4.Vect().Mag()), weight_t);
				  
				  }
				}
			  } // top mass window
			  else{
				  
					hist_jettau32AK8_mfail->Fill(tau32AK8,weight_t); 
				  
					float fmaxsubjetbtagvalue= -100;
					int sub_bf = -1;
			  
					if(nSubJet > 0){
			 
				    int nbsubjetsAK8 = 0;
				  
					for(unsigned isub=0; isub<nSubJet; isub++){
				
						if((int(isub)!=FatJet_subJetIdx1[topjetAK8])&&(int(isub)!=FatJet_subJetIdx2[topjetAK8])) continue;
						
						if(FatJet_subbtagId_DeepCSV[topjetAK8]>=0){
							if(int(isub) == FatJet_subbtagId_DeepCSV[topjetAK8]){
								hist_jetsubmassAK8_mfail->Fill(SubJet_mass[isub],weight_t);
							}
						}
						
						hist_jetsubptAK8_mfail->Fill(SubJet_pt[isub],weight);
						
						if(SubJet_btagDeepB[isub] > btagvalue){
							nbsubjetsAK8++;
							}
						
						}
				
						hist_jetsubbtagAK8_mfail->Fill(FatJet_subbtagCSVV2[topjetAK8],weight_t);
						hist_jetsubbtagAK8_deepCSV_mfail->Fill(FatJet_subbtagDeepB[topjetAK8],weight_t);
						
						if(FatJet_subbtagDeepB[topjetAK8] > btagvalue){
						hist_jetsubptAK8_btagpass_mfail->Fill(SubJet_pt[sub_bf],weight_t);
						}
						hist_nbsubjetAK8_mfail->Fill(nbsubjetsAK8,weight_t);
						
						hist_2D_tau32_subbtag_AK8_mfail->Fill(tau32AK8,FatJet_subbtagDeepB[topjetAK8],weight_t);
					}
				  
				  if(FatJet_msoftdrop[topjetAK8] > 60. && FatJet_msoftdrop[topjetAK8]<100.){
							hist_jettau21AK8->Fill(FatJet_tau2[topjetAK8]*1./FatJet_tau1[topjetAK8],weight_t);
						}
				  }
		  
	     if(FatJet_msoftdrop[topjetAK8] >= topmasslow && FatJet_msoftdrop[topjetAK8] < topmasshigh){
		  
		  if(tau32AK8 < tau32_cut) { top_wp = 0; }
		  if(tau32AK8 >= tau32_cut) { top_wp = 1; }
		  
		  if(FatJet_deepTag_TvsQCD[topjetAK8] >= deepak8_cut) { dtop_wp = 0; }
		  if(FatJet_deepTag_TvsQCD[topjetAK8] < deepak8_cut) { dtop_wp = 1; }
		  
		  if(FatJet_deepTagMD_TvsQCD[topjetAK8] >= deepak8_cut_md) { mdtop_wp = 0; }
		  if(FatJet_deepTagMD_TvsQCD[topjetAK8] < deepak8_cut_md) { mdtop_wp = 1; }
				 
		  }
		  else{
			     if(tau32AK8 < tau32_cut) { top_wp = 2; }
				 if(tau32AK8 >= tau32_cut) { top_wp = 3; }
				 
				 if(FatJet_deepTag_TvsQCD[topjetAK8] >= deepak8_cut) { dtop_wp = 2; }
				 if(FatJet_deepTag_TvsQCD[topjetAK8] < deepak8_cut) { dtop_wp = 3; }
				 
				 if(FatJet_deepTagMD_TvsQCD[topjetAK8] >= deepak8_cut_md) { mdtop_wp = 2; }
				 if(FatJet_deepTagMD_TvsQCD[topjetAK8] < deepak8_cut_md) { mdtop_wp = 3; }
			  
				 if(tau32AK8 < tau32_cut && FatJet_msoftdrop[topjetAK8] >= minmasscut && FatJet_msoftdrop[topjetAK8] < topmasslow) { top_wp = 4; }
				 if(tau32AK8 >= tau32_cut && FatJet_msoftdrop[topjetAK8] >= minmasscut && FatJet_msoftdrop[topjetAK8] < topmasslow) { top_wp = 5; }
				 
				 if(FatJet_deepTag_TvsQCD[topjetAK8] >= deepak8_cut && FatJet_msoftdrop[topjetAK8] >= minmasscut && FatJet_msoftdrop[topjetAK8] < topmasslow) { dtop_wp = 4; }
				 if(FatJet_deepTag_TvsQCD[topjetAK8] < deepak8_cut && FatJet_msoftdrop[topjetAK8] >= minmasscut && FatJet_msoftdrop[topjetAK8] < topmasslow) { dtop_wp = 5; }
				 
				 if(FatJet_deepTagMD_TvsQCD[topjetAK8] >= deepak8_cut_md && FatJet_msoftdrop[topjetAK8] >= minmasscut && FatJet_msoftdrop[topjetAK8] < topmasslow) { mdtop_wp = 4; }
				 if(FatJet_deepTagMD_TvsQCD[topjetAK8] < deepak8_cut_md && FatJet_msoftdrop[topjetAK8] >= minmasscut && FatJet_msoftdrop[topjetAK8] < topmasslow) { mdtop_wp = 5; }
			  
				 if(tau32AK8 < tau32_cut && FatJet_msoftdrop[topjetAK8] >= topmasshigh) { top_wp = 6; }
				 if(tau32AK8 >= tau32_cut && FatJet_msoftdrop[topjetAK8] >= topmasshigh) { top_wp = 7; }
				 
				 if(FatJet_deepTag_TvsQCD[topjetAK8] >= deepak8_cut && FatJet_msoftdrop[topjetAK8] >= topmasshigh) { dtop_wp = 6; }
				 if(FatJet_deepTag_TvsQCD[topjetAK8] < deepak8_cut && FatJet_msoftdrop[topjetAK8] >= topmasshigh) { dtop_wp = 7; }
				 
				 if(FatJet_deepTagMD_TvsQCD[topjetAK8] >= deepak8_cut_md && FatJet_msoftdrop[topjetAK8] >= topmasshigh) { mdtop_wp = 6; }
				 if(FatJet_deepTagMD_TvsQCD[topjetAK8] < deepak8_cut_md && FatJet_msoftdrop[topjetAK8] >= topmasshigh) { mdtop_wp = 7; }
			 
			  }
	  
	     int nsubjets_btagged = 0;
	     for(unsigned isub=0; isub<nSubJet; isub++){
				
				if((int(isub)!=FatJet_subJetIdx1[topjetAK8])&&(int(isub)!=FatJet_subJetIdx2[topjetAK8])) continue;
				if(SubJet_btagDeepB[isub] > btagvalue_deepCSV) { nsubjets_btagged++; }
			}
	
		 if(FatJet_subbtagDeepB[topjetAK8] >= btagvalue_deepCSV) { subjetb_wp = 1; }
		 else { subjetb_wp = 0; }
	 	 
		   }//topjetAK8>=0
		
	if(bjetAK4>=0){	
			
		TLorentzVector AK4vec4; 
		AK4vec4.SetPtEtaPhiM(Jet_pt[bjetAK4],Jet_eta[bjetAK4],Jet_phi[bjetAK4],Jet_mass[bjetAK4]);
		
		hist_bjetptAK4->Fill(Jet_pt[bjetAK4],weight);
		hist_tjetptAK8_bjetptAK4->Fill(FatJet_pt[topjetAK8],Jet_pt[bjetAK4],weight);
					
		if(isMC){
					
			hist_jetpartonflavAK4->Fill(Jet_partonFlavour[bjetAK4],weight);
			if(Jet_pt[bjetAK4] > AK4ptcut_fi){
				hist_jetpartonflavAK4_1->Fill(Jet_partonFlavour[bjetAK4],weight);
			}
				
			#ifdef Btagger_DeepJet
				if(Jet_btagDeepFlavB[bjetAK4] > btagvalue_deepFlavB){
					hist_jetpartonflavAK4_btagged->Fill(Jet_partonFlavour[bjetAK4],weight);
					if(Jet_pt[bjetAK4] > AK4ptcut_fi){
						hist_jetpartonflavAK4_btagged_1->Fill(Jet_partonFlavour[bjetAK4],weight);
					}
				}
			#endif
				
			#ifdef Btagger_DeepCSV
				if(Jet_btagDeepB[bjetAK4] > btagvalue_deepCSV){
					hist_jetpartonflavAK4_btagged->Fill(Jet_partonFlavour[bjetAK4],weight);
					if(Jet_pt[bjetAK4] > AK4ptcut_fi){
						hist_jetpartonflavAK4_btagged_1->Fill(Jet_partonFlavour[bjetAK4],weight);
					}
				}
			#endif	
				
			#ifdef Btagger_CSVv2
				if(Jet_btagCSVV2[bjetAK4] > btagvalue){
					hist_jetpartonflavAK4_btagged->Fill(Jet_partonFlavour[bjetAK4],weight);
					if(Jet_pt[bjetAK4] > AK4ptcut_fi){
						hist_jetpartonflavAK4_btagged_1->Fill(Jet_partonFlavour[bjetAK4],weight);
					}
				}
			#endif	
			}//isMC
			
			#ifdef Btagger_DeepJet
			if((Jet_btagDeepFlavB[bjetAK4] >= btagvalue_deepFlavB)){
				bjet_wp = 1;
			}
			else { bjet_wp = 0; }
			#endif
				  
			#ifdef Btagger_DeepCSV
			if((Jet_btagDeepB[bjetAK4] >= btagvalue_deepCSV)){
				bjet_wp = 1;
			}
			else { bjet_wp = 0; }
			#endif
				  
			#ifdef Btagger_CSVv2
			if(Jet_btagCSVV2[bjetAK4] >= btagvalue){
				bjet_wp = 1;
			}
			else { bjet_wp = 0; }
			#endif	
					
			int AK8matchb = Jet_MatchFatJet[bjetAK4];
				    
			TLorentzVector matchtjet;
				    
			if(AK8matchb>=0){
				
				matchtjet.SetPtEtaPhiM(FatJet_pt[AK8matchb],FatJet_eta[AK8matchb],FatJet_phi[AK8matchb],FatJet_mass[AK8matchb]);
						
				hist_biso_mass->Fill(FatJet_msoftdrop[AK8matchb],weight);
				hist_biso_tbmassdiff->Fill((matchtjet-AK4vec4).M(),weight);
				hist_biso_isomass->Fill(AK4vec4.M()/FatJet_msoftdrop[AK8matchb],weight);
				hist_biso_isopt->Fill(Jet_pt[bjetAK4]/FatJet_pt[AK8matchb],weight);
				hist_biso_TopAK8score->Fill(FatJet_deepTag_TvsQCD[AK8matchb],weight);
				hist_biso_TopAK8score_MD->Fill(FatJet_deepTagMD_TvsQCD[AK8matchb],weight);
				hist_biso_WAK8score->Fill(FatJet_deepTag_WvsQCD[AK8matchb],weight);
				hist_biso_WAK8score_MD->Fill(FatJet_deepTagMD_WvsQCD[AK8matchb],weight);
						
						
				if(top_wp>=0 && subjetb_wp>=0){
					
					hist_2d_biso_AK4mass_pt[top_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],AK4vec4.M(),weight);
					hist_2d_biso_AK8sdmass_pt[top_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_msoftdrop[AK8matchb],weight);
					hist_2d_biso_AK8TopScore_pt[top_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTag_TvsQCD[AK8matchb],weight);
					hist_2d_biso_AK8TopScore_MD_pt[top_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTagMD_TvsQCD[AK8matchb],weight);
					hist_2d_biso_AK8WScore_pt[top_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTag_WvsQCD[AK8matchb],weight);	
					hist_2d_biso_AK8WScore_MD_pt[top_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTagMD_WvsQCD[AK8matchb],weight);	
				
					hist_2d_biso_AK4mass_pt_dAK8[dtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],AK4vec4.M(),weight);
					hist_2d_biso_AK8sdmass_pt_dAK8[dtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_msoftdrop[AK8matchb],weight);
					hist_2d_biso_AK8TopScore_pt_dAK8[dtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTag_TvsQCD[AK8matchb],weight);
					hist_2d_biso_AK8TopScore_MD_pt_dAK8[dtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTagMD_TvsQCD[AK8matchb],weight);
					hist_2d_biso_AK8WScore_pt_dAK8[dtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTag_WvsQCD[AK8matchb],weight);	
					hist_2d_biso_AK8WScore_MD_pt_dAK8[dtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTagMD_WvsQCD[AK8matchb],weight);	
				
					hist_2d_biso_AK4mass_pt_mdAK8[mdtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],AK4vec4.M(),weight);
					hist_2d_biso_AK8sdmass_pt_mdAK8[mdtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_msoftdrop[AK8matchb],weight);
					hist_2d_biso_AK8TopScore_pt_mdAK8[mdtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTag_TvsQCD[AK8matchb],weight);
					hist_2d_biso_AK8TopScore_MD_pt_mdAK8[mdtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTagMD_TvsQCD[AK8matchb],weight);
					hist_2d_biso_AK8WScore_pt_mdAK8[mdtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTag_WvsQCD[AK8matchb],weight);	
					hist_2d_biso_AK8WScore_MD_pt_mdAK8[mdtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTagMD_WvsQCD[AK8matchb],weight);	
				
				}
					
				if(isMC){
					if(FatJet_MatchGenJet[AK8matchb]>=0){
						hist_bjetAK8genmass->Fill(FatJet_GenJetmass[AK8matchb],weight);
					}
				}
					
				if(Jet_pt[bjetAK4] > AK4ptcut_fi){
						
					double dR = delta2R(FatJet_eta[topjetAK8],FatJet_phi[topjetAK8],Jet_eta[bjetAK4],Jet_phi[bjetAK4]);
					double dPhi = PhiInRange(Jet_phi[bjetAK4]-FatJet_phi[topjetAK8]);
						
					hist_delR_btag_toptag->Fill(dR,weight);
					hist_delphi_btag_toptag->Fill(dPhi,weight);
					
					hist_biso_mass_1->Fill(FatJet_msoftdrop[AK8matchb],weight);
					hist_biso_tbmassdiff_1->Fill((matchtjet-AK4vec4).M(),weight);
					hist_biso_isomass_1->Fill(AK4vec4.M()/FatJet_msoftdrop[AK8matchb],weight);
					hist_biso_isopt_1->Fill(Jet_pt[bjetAK4]/FatJet_pt[AK8matchb],weight);
					hist_biso_TopAK8score_1->Fill(FatJet_deepTag_TvsQCD[AK8matchb],weight);
					hist_biso_TopAK8score_MD_1->Fill(FatJet_deepTagMD_TvsQCD[AK8matchb],weight);
					hist_biso_WAK8score_1->Fill(FatJet_deepTag_WvsQCD[AK8matchb],weight);
					hist_biso_WAK8score_MD_1->Fill(FatJet_deepTagMD_WvsQCD[AK8matchb],weight);
					
					if(FatJet_msoftdrop[AK8matchb]>60){
						hist_biso_TopAK8score_masscut_1->Fill(FatJet_deepTag_TvsQCD[AK8matchb],weight);
						hist_biso_TopAK8score_MD_masscut_1->Fill(FatJet_deepTagMD_TvsQCD[AK8matchb],weight);
					}
					
					if(bjet_wp>=0){
						hist_biso_mass_btag_1[bjet_wp]->Fill(FatJet_msoftdrop[AK8matchb],weight);
						hist_biso_TopAK8score_MD_btag_1[bjet_wp]->Fill(FatJet_deepTagMD_TvsQCD[AK8matchb],weight);
						hist_biso_WAK8score_MD_btag_1[bjet_wp]->Fill(FatJet_deepTagMD_WvsQCD[AK8matchb],weight);
					}
					
					if(top_wp>=0 && subjetb_wp>=0){
					
						hist_2d_biso_AK4mass_pt_1[top_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],AK4vec4.M(),weight);
						hist_2d_biso_AK8sdmass_pt_1[top_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_msoftdrop[AK8matchb],weight);
						hist_2d_biso_AK8TopScore_pt_1[top_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTag_TvsQCD[AK8matchb],weight);
						hist_2d_biso_AK8TopScore_MD_pt_1[top_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTagMD_TvsQCD[AK8matchb],weight);
						hist_2d_biso_AK8WScore_pt_1[top_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTag_WvsQCD[AK8matchb],weight);
						hist_2d_biso_AK8WScore_MD_pt_1[top_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTag_WvsQCD[AK8matchb],weight);
					
						hist_2d_biso_AK4mass_pt_dAK8_1[dtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],AK4vec4.M(),weight);
						hist_2d_biso_AK8sdmass_pt_dAK8_1[dtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_msoftdrop[AK8matchb],weight);
						hist_2d_biso_AK8TopScore_pt_dAK8_1[dtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTag_TvsQCD[AK8matchb],weight);
						hist_2d_biso_AK8TopScore_MD_pt_dAK8_1[dtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTagMD_TvsQCD[AK8matchb],weight);
						hist_2d_biso_AK8WScore_pt_dAK8_1[dtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTag_WvsQCD[AK8matchb],weight);
						hist_2d_biso_AK8WScore_MD_pt_dAK8_1[dtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTagMD_WvsQCD[AK8matchb],weight);
					
						hist_2d_biso_AK4mass_pt_mdAK8_1[mdtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],AK4vec4.M(),weight);
						hist_2d_biso_AK8sdmass_pt_mdAK8_1[mdtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_msoftdrop[AK8matchb],weight);
						hist_2d_biso_AK8TopScore_pt_mdAK8_1[mdtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTag_TvsQCD[AK8matchb],weight);
						hist_2d_biso_AK8TopScore_MD_pt_mdAK8_1[mdtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTagMD_TvsQCD[AK8matchb],weight);
						hist_2d_biso_AK8WScore_pt_mdAK8_1[mdtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTag_WvsQCD[AK8matchb],weight);
						hist_2d_biso_AK8WScore_MD_pt_mdAK8_1[mdtop_wp][subjetb_wp]->Fill(FatJet_pt[AK8matchb],FatJet_deepTagMD_WvsQCD[AK8matchb],weight);
					
					}
							
				    if(isMC){
						if(FatJet_MatchGenJet[AK8matchb]>=0){
							hist_bjetAK8genmass_1->Fill(FatJet_GenJetmass[AK8matchb],weight);
							if(FatJet_MatchGenJet[topjetAK8]>=0){
								hist_genmass_tbcor->Fill(FatJet_GenJetmass[AK8matchb],FatJet_GenJetmass[topjetAK8],weight);
							}
						}
					}
							
					float submass[2]={0}; 
					int isbcnt = 0;
					for(unsigned isub=0; isub<nSubJet; isub++){	
						if((int(isub)!=int(FatJet_subJetIdx1[AK8matchb]))&&(int(isub)!=int(FatJet_subJetIdx2[AK8matchb]))) continue;
							submass[isbcnt] = SubJet_mass[isub];
							isbcnt++;
							if(isbcnt>1) break;
						}
		
					if(isbcnt==2){
						(submass[0] > submass[1])? hist_2D_biso_AK8jet_subjetmass12->Fill(submass[0],submass[1],weight) : hist_2D_biso_AK8jet_subjetmass12->Fill(submass[1],submass[0],weight);
						if(submass[0] > submass[1]){
							hist_2D_biso_AK8jetsdmass_subjetmass1->Fill(FatJet_msoftdrop[AK8matchb],submass[0],weight);
							hist_2D_biso_AK8jetsdmass_subjetmass2->Fill(FatJet_msoftdrop[AK8matchb],submass[1],weight);
							}else{
								hist_2D_biso_AK8jetsdmass_subjetmass1->Fill(FatJet_msoftdrop[AK8matchb],submass[1],weight);
								hist_2D_biso_AK8jetsdmass_subjetmass2->Fill(FatJet_msoftdrop[AK8matchb],submass[0],weight);
								}
						}
						
					}
				}

	  } // bjetAK4>=0
			
	if(isMC){					
	  weight_t *= 1./sfwt_deepak8_md;
	}
	  // b jet from AK8 //
	  
	  bool b_foundAK8 = false;
	  int bjetAK8 = -1;
	  int bjetAK8_wp = -1;
	  
	  if(topjetAK8>=0){
      
      for(unsigned ijet=0; ijet<nFatJet; ijet++){
			  
			  double dR = delta2R(FatJet_eta[topjetAK8],FatJet_phi[topjetAK8],FatJet_eta[ijet],FatJet_phi[ijet]);
			  double dEta = (FatJet_eta[ijet]-FatJet_eta[topjetAK8]);
			  double dPhi = PhiInRange(FatJet_phi[ijet]-FatJet_phi[topjetAK8]);
			  
			  if(int(ijet) == topjetAK8) continue;
			  
			  if(dR>dR_cut){
				
				if(fabs(dPhi) > dphi_cut){	
					
					TLorentzVector bAK8vec4; 
					bAK8vec4.SetPtEtaPhiM(FatJet_pt[ijet],FatJet_eta[ijet],FatJet_phi[ijet],FatJet_mass[ijet]);
			   	  
					if(fabs(bAK8vec4.Rapidity()) < 2.5 && (FatJet_pt[ijet]>AK4ptcut_in)){
						
					if(!b_foundAK8){
						
						bjetAK8 = ijet;
						b_foundAK8 = true;  
				     
						if((FatJet_btagDeepB[bjetAK8] > btagvalue_deepCSV)){
							bjetAK8_wp = 1;
						}
						else { bjetAK8_wp = 0; }
						
						}
					}
				 }// dPhi
			  } // dR
		   }//ijet
		}// if topjetAK8 
	  
	  //AK8 b end
	  
	  if(topjetAK8>=0 && bjetAK8>=0){
		  
		  TLorentzVector tjet; 
		  tjet.SetPtEtaPhiM(FatJet_pt[topjetAK8],FatJet_eta[topjetAK8],FatJet_phi[topjetAK8],FatJet_mass[topjetAK8]);
		  TLorentzVector bjet; 
		  bjet.SetPtEtaPhiM(FatJet_pt[bjetAK8],FatJet_eta[bjetAK8],FatJet_phi[bjetAK8],FatJet_mass[bjetAK8]);
	  
		  int ittag = top_wp;	int ibtag = bjetAK8_wp; int isubbtag = subjetb_wp;
		  int idttag = dtop_wp; int imdttag = mdtop_wp; 
		  
		  if(ittag>=0 && ibtag>=0 && isubbtag>=0){
			  
			  if(FatJet_pt[bjetAK8] > AK4ptcut_fi){
				  if(FatJet_msoftdrop[bjetAK8] < bmasscut_fun(FatJet_pt[bjetAK8])){
					hist_tbmass_AK8[ittag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_t*sfwt_tau32);
					hist_AK8bCand_DeeptagScore[ittag][ibtag][isubbtag]->Fill(FatJet_deepTag_TvsQCD[bjetAK8],weight_t*sfwt_tau32);
					hist_AK8bCand_MD_DeeptagScore[ittag][ibtag][isubbtag]->Fill(FatJet_deepTagMD_TvsQCD[bjetAK8],weight_t*sfwt_tau32);
				  }
					hist_AK8bCand_DeeptagScore_1[ittag][ibtag][isubbtag]->Fill(FatJet_deepTag_TvsQCD[bjetAK8],weight_t*sfwt_tau32);
			  }
		  }
		  
		  if(idttag>=0 && ibtag>=0 && isubbtag>=0){
			  
			  if(FatJet_pt[bjetAK8] > AK4ptcut_fi){
				  if(FatJet_msoftdrop[bjetAK8] < bmasscut_fun(FatJet_pt[bjetAK8])){
					hist_tbmass_AK8_DeepAK8[idttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_t*sfwt_deepak8);
					hist_AK8bCand_DeeptagScore_DeepAK8[idttag][ibtag][isubbtag]->Fill(FatJet_deepTag_TvsQCD[bjetAK8],weight_t*sfwt_deepak8);
					hist_AK8bCand_MD_DeeptagScore_DeepAK8[idttag][ibtag][isubbtag]->Fill(FatJet_deepTagMD_TvsQCD[bjetAK8],weight_t*sfwt_deepak8);
				  }
					hist_AK8bCand_DeeptagScore_DeepAK8_1[idttag][ibtag][isubbtag]->Fill(FatJet_deepTag_TvsQCD[bjetAK8],weight_t*sfwt_deepak8);
			  }
		  }
		  
		  if(imdttag>=0 && ibtag>=0 && isubbtag>=0){
			  
			  if(FatJet_pt[bjetAK8] > AK4ptcut_fi){
				  if(FatJet_msoftdrop[bjetAK8] < bmasscut_fun(FatJet_pt[bjetAK8])){
					hist_tbmass_AK8_MDDeepAK8[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_t*sfwt_deepak8_md);
					hist_AK8bCand_DeeptagScore_MDDeepAK8[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTag_TvsQCD[bjetAK8],weight_t*sfwt_deepak8_md);
					hist_AK8bCand_MD_DeeptagScore_MDDeepAK8[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTagMD_TvsQCD[bjetAK8],weight_t*sfwt_deepak8_md);
				  }
					hist_AK8bCand_DeeptagScore_MDDeepAK8_1[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTag_TvsQCD[bjetAK8],weight_t*sfwt_deepak8_md);
			  }
		  }
	  }
	  
	  if(topjetAK8>=0 && bjetAK4>=0){
		  
		  hist_npv_final->Fill(PV_npvsGood,weight);
		  if(puWeight>1.e-9){
		  hist_npv_final_nopuwt->Fill(PV_npvsGood,weight*1./puWeight);
		  }
		  
		  TLorentzVector tjet; 
		  tjet.SetPtEtaPhiM(FatJet_pt[topjetAK8],FatJet_eta[topjetAK8],FatJet_phi[topjetAK8],FatJet_mass[topjetAK8]);
		  TLorentzVector bjet; 
		  bjet.SetPtEtaPhiM(Jet_pt[bjetAK4],Jet_eta[bjetAK4],Jet_phi[bjetAK4],Jet_mass[bjetAK4]);
	  
		  TLorentzVector tjet_jer_up, tjet_jer_dn;
		  tjet_jer_up.SetPtEtaPhiM(FatJet_pt_jer_up[topjetAK8],FatJet_eta[topjetAK8],FatJet_phi[topjetAK8],FatJet_mass_jer_up[topjetAK8]);
		  tjet_jer_dn.SetPtEtaPhiM(FatJet_pt_jer_dn[topjetAK8],FatJet_eta[topjetAK8],FatJet_phi[topjetAK8],FatJet_mass_jer_dn[topjetAK8]);
		  TLorentzVector bjet_jer_up, bjet_jer_dn;
		  bjet_jer_up.SetPtEtaPhiM(Jet_pt_jer_up[bjetAK4],Jet_eta[bjetAK4],Jet_phi[bjetAK4],Jet_mass_jer_up[bjetAK4]);
		  bjet_jer_dn.SetPtEtaPhiM(Jet_pt_jer_dn[bjetAK4],Jet_eta[bjetAK4],Jet_phi[bjetAK4],Jet_mass_jer_dn[bjetAK4]);
		  
		  TLorentzVector tjet_jms_up, tjet_jms_dn;
		  tjet_jms_up.SetPtEtaPhiM(FatJet_pt[topjetAK8],FatJet_eta[topjetAK8],FatJet_phi[topjetAK8],FatJet_mass_jms_up[topjetAK8]);
		  tjet_jms_dn.SetPtEtaPhiM(FatJet_pt[topjetAK8],FatJet_eta[topjetAK8],FatJet_phi[topjetAK8],FatJet_mass_jms_dn[topjetAK8]);
		  TLorentzVector bjet_jms_up, bjet_jms_dn;
		  bjet_jms_up.SetPtEtaPhiM(Jet_pt[bjetAK4],Jet_eta[bjetAK4],Jet_phi[bjetAK4],Jet_mass_jms_up[bjetAK4]);
		  bjet_jms_dn.SetPtEtaPhiM(Jet_pt[bjetAK4],Jet_eta[bjetAK4],Jet_phi[bjetAK4],Jet_mass_jms_dn[bjetAK4]);  
		  
		  TLorentzVector tjet_jmr_up, tjet_jmr_dn;
		  tjet_jmr_up.SetPtEtaPhiM(FatJet_pt[topjetAK8],FatJet_eta[topjetAK8],FatJet_phi[topjetAK8],FatJet_mass_jmr_up[topjetAK8]);
		  tjet_jmr_dn.SetPtEtaPhiM(FatJet_pt[topjetAK8],FatJet_eta[topjetAK8],FatJet_phi[topjetAK8],FatJet_mass_jmr_dn[topjetAK8]);
		  TLorentzVector bjet_jmr_up, bjet_jmr_dn;
		  bjet_jmr_up.SetPtEtaPhiM(Jet_pt[bjetAK4],Jet_eta[bjetAK4],Jet_phi[bjetAK4],Jet_mass_jmr_up[bjetAK4]);
		  bjet_jmr_dn.SetPtEtaPhiM(Jet_pt[bjetAK4],Jet_eta[bjetAK4],Jet_phi[bjetAK4],Jet_mass_jmr_dn[bjetAK4]);
		  
		  TLorentzVector tjet_jes_up[njesmax], tjet_jes_dn[njesmax];
		  for(int ijes=0; ijes<njesmax; ijes++){
				tjet_jes_up[ijes].SetPtEtaPhiM(FatJet_pt_jes_up[ijes][topjetAK8],FatJet_eta[topjetAK8],FatJet_phi[topjetAK8],FatJet_mass_jes_up[ijes][topjetAK8]);
				tjet_jes_dn[ijes].SetPtEtaPhiM(FatJet_pt_jes_dn[ijes][topjetAK8],FatJet_eta[topjetAK8],FatJet_phi[topjetAK8],FatJet_mass_jes_dn[ijes][topjetAK8]);
			}
		  TLorentzVector bjet_jes_up[njesmax], bjet_jes_dn[njesmax];   
		  for(int ijes=0; ijes<njesmax; ijes++){
				bjet_jes_up[ijes].SetPtEtaPhiM(Jet_pt_jes_up[ijes][bjetAK4],Jet_eta[bjetAK4],Jet_phi[bjetAK4],Jet_mass_jes_up[ijes][bjetAK4]);
				bjet_jes_dn[ijes].SetPtEtaPhiM(Jet_pt_jes_dn[ijes][bjetAK4],Jet_eta[bjetAK4],Jet_phi[bjetAK4],Jet_mass_jes_dn[ijes][bjetAK4]);
			}
	  
		  int ittag = top_wp;	int ibtag = bjet_wp; int isubbtag = subjetb_wp;
		  int idttag = dtop_wp; int imdttag = mdtop_wp;
		  
		  if(ittag>=0 && ibtag>=0 && isubbtag>=0){
		  
		  float tau32 = FatJet_tau3[topjetAK8]*1./TMath::Max(float(1.e-6),FatJet_tau2[topjetAK8]);
		  int tautag = -1;
		  if(tau32 < tau32_cut) { tautag = 0; }
		  if(tau32 >= tau32_cut) { tautag = 1; }
		  
		  float exwt = 1.;
		  float exwt2 = 1.;
		 
		  float weight1, weight2;
		 
		  weight1 = weight_t * sfwt_tau32 * exwt;
		  weight2 = weight_t * sfwt_tau32 * exwt2;
		  
		  if(Jet_pt[bjetAK4] > AK4ptcut_fi){
			  
			  if((Jet_MatchFatJet[bjetAK4]>=0)&&(FatJet_msoftdrop[Jet_MatchFatJet[bjetAK4]] < bmasscut_fun(FatJet_pt[Jet_MatchFatJet[bjetAK4]]))){
			
					hist_tbmass[ittag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight1);
					hist_tbpt[ittag][ibtag][isubbtag]->Fill((tjet+bjet).Pt(),weight1);
					hist_tbrap[ittag][ibtag][isubbtag]->Fill((tjet+bjet).Rapidity(),weight1);
					hist_tbpt_frac[ittag][ibtag][isubbtag]->Fill(TMath::Min(tjet.Pt(),bjet.Pt())/(tjet.Pt()+bjet.Pt()),weight1);
		    			if(isMC){ 
						hist_partonflav_AK4[ittag][ibtag][isubbtag]->Fill(Jet_partonFlavour[bjetAK4],weight1);
					}
					hist_2D_sdmass_tbmass[ibtag][tautag]->Fill(FatJet_msoftdrop[topjetAK8],(tjet+bjet).M(),weight1);
					hist_2D_DeepAK8_tbmass[ibtag]->Fill(FatJet_deepTag_TvsQCD[topjetAK8],(tjet+bjet).M(),weight1);
					hist_2D_MDDeepAK8_tbmass[ibtag]->Fill(FatJet_deepTagMD_TvsQCD[topjetAK8],(tjet+bjet).M(),weight1);
					
					hist_2D_sdmass_subbtag_CSVv2_AK8[ibtag][tautag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_subbtagCSVV2[topjetAK8],weight1);
					hist_2D_sdmass_subbtag_DeepCSV_AK8[ibtag][tautag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_subbtagDeepB[topjetAK8],weight1);
					hist_2D_sdmass_tau32_AK8_reg[ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],tau32,weight1);
					#ifdef Btagger_DeepJet
						hist_2D_pt_btag_AK4[isubbtag][ittag]->Fill(Jet_pt[bjetAK4],Jet_btagDeepFlavB[bjetAK4],weight1);
					#endif
					#ifdef Btagger_DeepCSV
						hist_2D_pt_btag_AK4[isubbtag][ittag]->Fill(Jet_pt[bjetAK4],Jet_btagDeepB[bjetAK4],weight1);
					#endif
					#ifdef Btagger_CSVv2
						hist_2D_pt_btag_AK4[isubbtag][ittag]->Fill(Jet_pt[bjetAK4],Jet_btagCSVV2[bjetAK4],weight1);
					#endif
					
					hist_2D_sdmass_deeptopscore_AK8[ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_deepTag_TvsQCD[topjetAK8],weight1);
					hist_2D_sdmass_deeptopscore_md_AK8[ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight1);
		  
					hist_topjetpt_sel[ittag][ibtag][isubbtag]->Fill(FatJet_pt[topjetAK8],weight1);
					hist_topjetsdmass_sel[ittag][ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],weight1);
					hist_topjetdeeptopscore_sel[ittag][ibtag][isubbtag]->Fill(FatJet_deepTag_TvsQCD[topjetAK8],weight1);
					hist_topjetdeepmdtopscore_sel[ittag][ibtag][isubbtag]->Fill(FatJet_deepTagMD_TvsQCD[topjetAK8],weight1);
					hist_topjetsbtag_sel[ittag][ibtag][isubbtag]->Fill(FatJet_subbtagCSVV2[topjetAK8],weight1);
					hist_topjetsdeepCSV_sel[ittag][ibtag][isubbtag]->Fill(FatJet_subbtagDeepB[topjetAK8],weight1);
					hist_topjettau32_sel[ittag][ibtag][isubbtag]->Fill(tau32,weight1);
		  
					hist_bjetpt_sel[ittag][ibtag][isubbtag]->Fill(Jet_pt[bjetAK4],weight1);
					hist_bjetmass_sel[ittag][ibtag][isubbtag]->Fill(bjet.M(),weight1);
					#ifdef Btagger_DeepJet
						hist_bjetbtag_sel[ittag][ibtag][isubbtag]->Fill(Jet_btagDeepFlavB[bjetAK4],weight1);
					#endif
					#ifdef Btagger_DeepCSV
						hist_bjetbtag_sel[ittag][ibtag][isubbtag]->Fill(Jet_btagDeepB[bjetAK4],weight1);
					#endif
					#ifdef Btagger_CSVv2
						hist_bjetbtag_sel[ittag][ibtag][isubbtag]->Fill(Jet_btagCSVV2[bjetAK4],weight1);
					#endif
						
					if(Jet_MatchFatJet[bjetAK4] >= 0){
					hist_bjetmatch_AK8Topscore_sel[ittag][ibtag][isubbtag]->Fill(FatJet_deepTag_TvsQCD[Jet_MatchFatJet[bjetAK4]],weight1);
					hist_bjetmatch_AK8Wscore_sel[ittag][ibtag][isubbtag]->Fill(FatJet_deepTag_WvsQCD[Jet_MatchFatJet[bjetAK4]],weight1);
					}
					hist_3D_topjetsdmass_topjettau32_tbmass[ittag][ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],tau32,(tjet+bjet).M(),weight1);
					
					hist_ht_tau32[imdttag][ibtag][isubbtag]->Fill(Event_Ht,weight1);
					
					} // cut on AK8 jet matched to b candidate (mass cut)
			  
			  } // b pt cut
		   } // whether all have logical WP 
		   
		   
 if(idttag>=0 && ibtag>=0 && isubbtag>=0){
		  
		  float tau32 = FatJet_tau3[topjetAK8]*1./TMath::Max(float(1.e-6),FatJet_tau2[topjetAK8]);
		  int deeptag = -1;
		  if(FatJet_deepTag_TvsQCD[topjetAK8] < deepak8_cut) { deeptag = 0; }
		  if(FatJet_deepTag_TvsQCD[topjetAK8] >= deepak8_cut) { deeptag = 1; }
		  
		  float exwt = 1.;
		  float exwt2 = 1.;
		 
		  float weight1, weight2;
		  
		   if(FakeAnalysis){
			weight1 = weight_t * sfwt_deepak8 * exwt;
			weight2 = weight_t * sfwt_deepak8 * exwt2;
		    }else{
			    weight1 = weight_t * sfwt_deepak8;
				weight2 = weight_t * sfwt_deepak8 ;
			   }
		  
		  if(Jet_pt[bjetAK4] > AK4ptcut_fi){
		    
		  if((Jet_MatchFatJet[bjetAK4]>=0)&&(FatJet_msoftdrop[Jet_MatchFatJet[bjetAK4]] < bmasscut_fun(FatJet_pt[Jet_MatchFatJet[bjetAK4]]))){
		
				hist_tbmass_DeepAK8[idttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight1);
				hist_tbpt_DeepAK8[idttag][ibtag][isubbtag]->Fill((tjet+bjet).Pt(),weight1);
				hist_tbrap_DeepAK8[idttag][ibtag][isubbtag]->Fill((tjet+bjet).Rapidity(),weight1);
			        if(isMC){	
					hist_partonflav_AK4_DeepAK8[idttag][ibtag][isubbtag]->Fill(Jet_partonFlavour[bjetAK4],weight1);
				}
				hist_2D_sdmass_tbmass_deepAK8[ibtag][deeptag]->Fill(FatJet_msoftdrop[topjetAK8],(tjet+bjet).M(),weight1);
					
				hist_2D_sdmass_subbtag_CSVv2_AK8_DeepAK8[ibtag][deeptag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_subbtagCSVV2[topjetAK8],weight1);
				hist_2D_sdmass_subbtag_DeepCSV_AK8_DeepAK8[ibtag][deeptag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_subbtagDeepB[topjetAK8],weight1);
				
				#ifdef Btagger_DeepJet
					hist_2D_pt_btag_AK4_DeepAK8[isubbtag][idttag]->Fill(Jet_pt[bjetAK4],Jet_btagDeepFlavB[bjetAK4],weight1);
				#endif
				#ifdef Btagger_DeepCSV
					hist_2D_pt_btag_AK4_DeepAK8[isubbtag][idttag]->Fill(Jet_pt[bjetAK4],Jet_btagDeepB[bjetAK4],weight1);
				#endif
				#ifdef Btagger_CSVv2
					hist_2D_pt_btag_AK4_DeepAK8[isubbtag][idttag]->Fill(Jet_pt[bjetAK4],Jet_btagCSVV2[bjetAK4],weight1);
				#endif
				
				hist_2D_sdmass_deeptopscore_AK8_DeepAK8[ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_deepTag_TvsQCD[topjetAK8],weight1);
				hist_2D_sdmass_deeptopscore_md_AK8_DeepAK8[ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight1);
		  
				hist_topjetpt_sel_DeepAK8[idttag][ibtag][isubbtag]->Fill(FatJet_pt[topjetAK8],weight1);
				hist_topjetsdmass_sel_DeepAK8[idttag][ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],weight1);
				hist_topjetdeeptopscore_sel_DeepAK8[idttag][ibtag][isubbtag]->Fill(FatJet_deepTag_TvsQCD[topjetAK8],weight1);
				hist_topjetdeepmdtopscore_sel_DeepAK8[idttag][ibtag][isubbtag]->Fill(FatJet_deepTagMD_TvsQCD[topjetAK8],weight1);
				hist_topjetsbtag_sel_DeepAK8[idttag][ibtag][isubbtag]->Fill(FatJet_subbtagCSVV2[topjetAK8],weight1);
				hist_topjetsdeepCSV_sel_DeepAK8[idttag][ibtag][isubbtag]->Fill(FatJet_subbtagDeepB[topjetAK8],weight1);
				hist_topjettau32_sel_DeepAK8[idttag][ibtag][isubbtag]->Fill(tau32,weight1);
		  
				int imsdbin = getbinid(FatJet_msoftdrop[topjetAK8],ntopsdmassbins,topsdmassbins);
				if(imsdbin>=0){
					hist_topjetpt_msdbin_sel_DeepAK8[idttag][ibtag][isubbtag][imsdbin]->Fill(FatJet_pt[topjetAK8],weight1);
					}
		  
				hist_bjetpt_sel_DeepAK8[idttag][ibtag][isubbtag]->Fill(Jet_pt[bjetAK4],weight1);
				hist_bjetmass_sel_DeepAK8[idttag][ibtag][isubbtag]->Fill(bjet.M(),weight1);
				#ifdef Btagger_DeepJet
					hist_bjetbtag_sel_DeepAK8[idttag][ibtag][isubbtag]->Fill(Jet_btagDeepFlavB[bjetAK4],weight1);
				#endif
				#ifdef Btagger_DeepCSV
					hist_bjetbtag_sel_DeepAK8[idttag][ibtag][isubbtag]->Fill(Jet_btagDeepB[bjetAK4],weight1);
				#endif
				#ifdef Btagger_CSVv2
					hist_bjetbtag_sel_DeepAK8[idttag][ibtag][isubbtag]->Fill(Jet_btagCSVV2[bjetAK4],weight1);
				#endif
						
				if(Jet_MatchFatJet[bjetAK4] >= 0){
					hist_bjetmatch_AK8Topscore_sel_DeepAK8[idttag][ibtag][isubbtag]->Fill(FatJet_deepTag_TvsQCD[Jet_MatchFatJet[bjetAK4]],weight1);
					hist_bjetmatch_AK8Wscore_sel_DeepAK8[idttag][ibtag][isubbtag]->Fill(FatJet_deepTag_WvsQCD[Jet_MatchFatJet[bjetAK4]],weight1);
						}
		  
				hist_ht_DeepAK8[imdttag][ibtag][isubbtag]->Fill(Event_Ht,weight1);
		  
					} // mass cut on AK8 jet matched to b candidate		
			  
			 }  // b pt cut
		   } // whether all have logical WP 

 if(imdttag>=0 && ibtag>=0 && isubbtag>=0){
		  
		  float m_tb = (tjet+bjet).M();
		  
		  float tau32 = FatJet_tau3[topjetAK8]*1./TMath::Max(float(1.e-6),FatJet_tau2[topjetAK8]);
		  
		  int deeptag = -1;
		  if(FatJet_deepTagMD_TvsQCD[topjetAK8] < deepak8_cut_md)  { deeptag = 0; }
		  if(FatJet_deepTagMD_TvsQCD[topjetAK8] >= deepak8_cut_md) { deeptag = 1; }
		  
		  float exwt = 1.;
		  float exwt2 = 1.;
		  float exwt_sys = 0;
		  float exwt_sys2 = 0;
		  
		  float cor_bsd = 1;
		 
		  float weight1, weight2, weight3;
		  
		  int ieta = getbinid(fabs(Jet_eta[bjetAK4]),netabins,etabins);
		  
		    weight1 = weight_t * sfwt_deepak8_md * cor_bsd * exwt ;
		    weight2 = weight_t * sfwt_deepak8_md * exwt ; // without b sd mass correction
			
			int iptsum = getbinid(tjet.Pt()+bjet.Pt(),nohtbins,htbins);
			int imsdbin = getbinid(FatJet_msoftdrop[topjetAK8],ntopsdmassbins,topsdmassbins);
			float trig_wt = 1;
			if(iptsum>=0 && imsdbin>=0 && isMC){ trig_wt = trig_weight_factors[imsdbin][iptsum]; }
			float weight_trig = weight1 * trig_wt;
			
			float weightup = weight_t * sfwt_deepak8_md * cor_bsd * (exwt+exwt_sys) ;
			float weightdn = weight_t * sfwt_deepak8_md * cor_bsd * (exwt-exwt_sys) ;
			
			float weight_dAK8up = weight1 * sfwt_deepak8_md_up *1./sfwt_deepak8_md;
			float weight_dAK8dn = weight1 * sfwt_deepak8_md_dn *1./sfwt_deepak8_md;
			
			float weight_puup = weight1 * puWeightUp * 1./puWeight;
			float weight_pudn = weight1 * puWeightDown * 1./puWeight;
			
			float weight_prefireup = weight1 * PrefireWeight_Up *1./ PrefireWeight;
			float weight_prefiredn = weight1 * PrefireWeight_Down *1./ PrefireWeight;
			
			float weight_noptw;
			float weight_alphaup, weight_alphadn, weight_betaup, weight_betadn;
			weight_noptw = 1;
			weight_alphaup = 1; weight_alphadn = 1; weight_betaup = 1; weight_betadn = 1;
			
			if(isTTBar && isMC){
				weight_noptw = weight1*1./SF_toppt;
				if(ngtop==2){
				weight_alphaup = weight1 * SF_TOP(1.5*0.0615,0.0005,GenPart_pt[igtop[0]],GenPart_pt[igtop[1]])/SF_toppt;
				weight_alphadn = weight1 * SF_TOP(0.5*0.0615,0.0005,GenPart_pt[igtop[0]],GenPart_pt[igtop[1]])/SF_toppt;
				weight_betaup = weight1 * SF_TOP(0.0615,1.5*0.0005,GenPart_pt[igtop[0]],GenPart_pt[igtop[1]])/SF_toppt;
				weight_betadn = weight1 * SF_TOP(0.0615,0.5*0.0005,GenPart_pt[igtop[0]],GenPart_pt[igtop[1]])/SF_toppt;
				}
			}else{
					weight_noptw = 1;
					weight_alphaup = weight_alphadn = 1;
					weight_betaup = weight_betadn = 1;
				 }
				 
			float weight_bcor_up = weight1 * (cor_bsd + 0.5*fabs(1-cor_bsd)) *1./cor_bsd;
			float weight_bcor_dn = weight1 * (cor_bsd - 0.5*fabs(1-cor_bsd)) *1./cor_bsd;
				 
			
			weight3 = weight_t * sfwt_deepak8_md * cor_bsd * exwt2 ;
			
			float weightup3 = weight_t * sfwt_deepak8_md * cor_bsd * (exwt2+exwt_sys2) ;
			float weightdn3 = weight_t * sfwt_deepak8_md * cor_bsd * (exwt2-exwt_sys2) ;
			
			float weight_noptw3;
			if(isTTBar && isMC){
				weight_noptw3 = weight3*1./SF_toppt;
			}else{
					weight_noptw3 = 1;
				 }	 
			
			float weight_bcor_up3 = weight3 * (cor_bsd + 0.5*fabs(1-cor_bsd)) *1./cor_bsd;
			float weight_bcor_dn3 = weight3 * (cor_bsd - 0.5*fabs(1-cor_bsd)) *1./cor_bsd;
			
		  if(Jet_pt[bjetAK4] > AK4ptcut_fi){
		    
		  if((Jet_MatchFatJet[bjetAK4]>=0)&&(FatJet_msoftdrop[Jet_MatchFatJet[bjetAK4]] < bmasscut_fun(FatJet_pt[Jet_MatchFatJet[bjetAK4]]))){	
			
				hist_tbmass_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight1);
				hist_tbpt_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).Pt(),weight1);
				hist_tbrap_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).Rapidity(),weight1);
				
				hist_tbmass_md_DeepAK8_fine[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight1);
				hist_tbmass_md_DeepAK8_fine_puup[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_puup);
				hist_tbmass_md_DeepAK8_fine_pudn[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_pudn);
				h2d_mtb_npu[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),int(Pileup_nTrueInt),weight1);
				
				hist_tbmass_md_DeepAK8_trigwt[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_trig);
				
				hist_tbmass_md_DeepAK8_puup[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_puup);
				hist_tbmass_md_DeepAK8_pudn[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_pudn);
				
				hist_tbmass_md_DeepAK8_prefireup[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_prefireup);
				hist_tbmass_md_DeepAK8_prefiredn[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_prefiredn);
			
				hist_tbmass_md_DeepAK8_pfup[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weightup);
				hist_tbmass_md_DeepAK8_pfdn[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weightdn);
				
				hist_tbmass_md_DeepAK8_2[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight3);
				
				hist_tbmass_md_DeepAK8_pfup_2[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weightup3);
				hist_tbmass_md_DeepAK8_pfdn_2[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weightdn3);
			
				
				if(fabs(bjet.Rapidity()-tjet.Rapidity()) < rapidity_cut){
					hist_tbmass_md_DeepAK8_dY[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight1);
					}
				
				if(FakeAnalysis){
				
				hist_tbmass_md_DeepAK8_bcorup[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_bcor_up);
				hist_tbmass_md_DeepAK8_bcordn[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_bcor_dn);
				
				hist_tbmass_md_DeepAK8_bcorup_2[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_bcor_up3);
				hist_tbmass_md_DeepAK8_bcordn_2[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_bcor_dn3);
				
				}
				
				hist_2D_sdmass_tbmass_mddeepAK8[ibtag][deeptag]->Fill(FatJet_msoftdrop[topjetAK8],(tjet+bjet).M(),weight1);
					
				hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8[ibtag][deeptag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_subbtagCSVV2[topjetAK8],weight1);
				hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8[ibtag][deeptag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_subbtagDeepB[topjetAK8],weight1);
				
				hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8[ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_deepTag_TvsQCD[topjetAK8],weight1);
				hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8[ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight1);
				
				int iptbin = getbinid(FatJet_pt[topjetAK8],ntopptbins,topptbins);
				hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_ptbin[ibtag][isubbtag][iptbin]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight1);
				hist_tbmass_md_DeepAK8_ptbin[imdttag][ibtag][isubbtag][iptbin]->Fill((tjet+bjet).M(),weight1);
		 
				hist_topjetpt_sel_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(FatJet_pt[topjetAK8],weight1);
				if(imsdbin>=0){
					hist_topjetpt_msdbin_sel_md_DeepAK8[imdttag][ibtag][isubbtag][imsdbin]->Fill(FatJet_pt[topjetAK8],weight1);
					}
					
				hist_topjetsdmass_sel_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],weight1);
				hist_topjetdeeptopscore_sel_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTag_TvsQCD[topjetAK8],weight1);
				hist_topjetdeepmdtopscore_sel_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTagMD_TvsQCD[topjetAK8],weight1);
				hist_topjetsbtag_sel_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(FatJet_subbtagCSVV2[topjetAK8],weight1);
				hist_topjetsdeepCSV_sel_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(FatJet_subbtagDeepB[topjetAK8],weight1);
				hist_topjettau32_sel_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(tau32,weight1);
				if(FatJet_MatchJet[topjetAK8]>=0){
					hist_partonflav_AK8_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(Jet_partonFlavour[FatJet_MatchJet[topjetAK8]],weight1);
				}
				
				if(isMC){
				
				#ifndef LHAPDFX
				
				if(nLHEPdfWeight>0){
					
				for(unsigned ipdf=0; ipdf<nLHEPdfWeight; ipdf++){
					double weight_pdf;
					weight_pdf = weight1*LHEPdfWeight[int(ipdf)];
					hist_tbmass_md_DeepAK8_PDF[imdttag][ibtag][isubbtag][int(ipdf)]->Fill((tjet+bjet).M(),weight_pdf);
					hist_topjetsdmass_sel_md_DeepAK8_PDF[imdttag][ibtag][isubbtag][int(ipdf)]->Fill(FatJet_msoftdrop[topjetAK8],weight_pdf);
					}
				}
				
				#endif
				
				#ifdef LHAPDFX
				// npdfmax = 1 + LHAPDF::numberPDF();
				for (int jk=0; jk< TMath::Min(npdfmax,1 + LHAPDF::numberPDF()) ;jk++) {
				
					LHAPDF::initPDF(jk);
					float ypdf1 = LHAPDF::xfx(Generator_x1,Generator_scalePDF, Generator_id1);
					float ypdf2 = LHAPDF::xfx(Generator_x2,Generator_scalePDF, Generator_id2);
					double weight_pdf = weight1*ypdf1*ypdf2 ;//*1./(Generator_xpdf1*Generator_xpdf2);
					hist_tbmass_md_DeepAK8_PDF[imdttag][ibtag][isubbtag][jk]->Fill((tjet+bjet).M(),weight_pdf);
                    hist_topjetsdmass_sel_md_DeepAK8_PDF[imdttag][ibtag][isubbtag][jk]->Fill(FatJet_msoftdrop[topjetAK8],weight_pdf);

				}
				#endif
				
				if(nLHEScaleWeight>0){
		  
				for(unsigned iscale=0; iscale<nLHEScaleWeight; iscale++){
					double weight_scale;
					weight_scale = weight1*LHEScaleWeight[int(iscale)];
					hist_tbmass_md_DeepAK8_Scale[imdttag][ibtag][isubbtag][int(iscale)]->Fill((tjet+bjet).M(),weight_scale);
					hist_topjetsdmass_sel_md_DeepAK8_Scale[imdttag][ibtag][isubbtag][int(iscale)]->Fill(FatJet_msoftdrop[topjetAK8],weight_scale);
					}
					
				}
				
				if(nPSWeight>0){
		  
				for(unsigned ips=0; ips<nPSWeight; ips++){
					double weight_ps;
					weight_ps = weight1*PSWeight[int(ips)];
					hist_tbmass_md_DeepAK8_PS[imdttag][ibtag][isubbtag][int(ips)]->Fill((tjet+bjet).M(),weight_ps);
					hist_topjetsdmass_sel_md_DeepAK8_PS[imdttag][ibtag][isubbtag][int(ips)]->Fill(FatJet_msoftdrop[topjetAK8],weight_ps);
					}
				if(nPSWeight<npsmax){
					for(unsigned ips=nPSWeight; ips<npsmax; ips++){
						hist_tbmass_md_DeepAK8_PS[imdttag][ibtag][isubbtag][int(ips)]->Fill((tjet+bjet).M(),weight1);
						hist_topjetsdmass_sel_md_DeepAK8_PS[imdttag][ibtag][isubbtag][int(ips)]->Fill(FatJet_msoftdrop[topjetAK8],weight1);
						}
					}	
					
				}	
				  
				for(int ibunc=0; ibunc<nbtagmax; ibunc++){
					
					double weight_btag_up;
					weight_btag_up = weight1*btag_weight_unc_up[ibunc]*1./btag_weight;
					hist_tbmass_md_DeepAK8_Btag_up[imdttag][ibtag][isubbtag][ibunc]->Fill((tjet+bjet).M(),weight_btag_up);
					hist_topjetsdmass_sel_md_DeepAK8_Btag_up[imdttag][ibtag][isubbtag][ibunc]->Fill(FatJet_msoftdrop[topjetAK8],weight_btag_up);
					hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_Btag_up[ibtag][isubbtag][ibunc]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight_btag_up);
					hist_bjetpt_sel_md_DeepAK8_Btag_up[imdttag][ibtag][isubbtag][ibunc]->Fill(Jet_pt[bjetAK4],weight_btag_up);
					hist_2D_pt_btag_AK4_md_DeepAK8_Btag_up[isubbtag][imdttag][ibunc]->Fill(Jet_pt[bjetAK4],Jet_btagDeepFlavB[bjetAK4],weight_btag_up);
					if(ieta>=0){
						hist_2D_pt_btag_AK4_md_DeepAK8_beta_Btag_up[isubbtag][imdttag][ieta][ibunc]->Fill(Jet_pt[bjetAK4],Jet_btagDeepFlavB[bjetAK4],weight_btag_up);
						}
					
					double weight_btag_dn;
					weight_btag_dn = weight1*btag_weight_unc_dn[ibunc]*1./btag_weight;;
					hist_tbmass_md_DeepAK8_Btag_dn[imdttag][ibtag][isubbtag][ibunc]->Fill((tjet+bjet).M(),weight_btag_dn);
					hist_topjetsdmass_sel_md_DeepAK8_Btag_dn[imdttag][ibtag][isubbtag][ibunc]->Fill(FatJet_msoftdrop[topjetAK8],weight_btag_dn);
					hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_Btag_dn[ibtag][isubbtag][ibunc]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight_btag_dn);
					hist_bjetpt_sel_md_DeepAK8_Btag_dn[imdttag][ibtag][isubbtag][ibunc]->Fill(Jet_pt[bjetAK4],weight_btag_dn);
					hist_2D_pt_btag_AK4_md_DeepAK8_Btag_dn[isubbtag][imdttag][ibunc]->Fill(Jet_pt[bjetAK4],Jet_btagDeepFlavB[bjetAK4],weight_btag_dn);
					if(ieta>=0){
						hist_2D_pt_btag_AK4_md_DeepAK8_beta_Btag_dn[isubbtag][imdttag][ieta][ibunc]->Fill(Jet_pt[bjetAK4],Jet_btagDeepFlavB[bjetAK4],weight_btag_dn);
						}
					
					}  
				 
				 hist_tbmass_md_DeepAK8_JER_up[imdttag][ibtag][isubbtag]->Fill((tjet_jer_up+bjet_jer_up).M(),weight1);
				 hist_tbmass_md_DeepAK8_JER_dn[imdttag][ibtag][isubbtag]->Fill((tjet_jer_dn+bjet_jer_dn).M(),weight1);
				 hist_topjetsdmass_sel_md_DeepAK8_JER_up[imdttag][ibtag][isubbtag]->Fill(FatJet_msoftdrop_jer_up[topjetAK8],weight1);
				 hist_topjetsdmass_sel_md_DeepAK8_JER_dn[imdttag][ibtag][isubbtag]->Fill(FatJet_msoftdrop_jer_dn[topjetAK8],weight1);
				 hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JER_up[ibtag][isubbtag]->Fill(FatJet_msoftdrop_jer_up[topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight1);
				 hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JER_dn[ibtag][isubbtag]->Fill(FatJet_msoftdrop_jer_dn[topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight1);
					
				 hist_tbmass_md_DeepAK8_JMR_up[imdttag][ibtag][isubbtag]->Fill((tjet_jmr_up+bjet_jmr_up).M(),weight1);
				 hist_tbmass_md_DeepAK8_JMR_dn[imdttag][ibtag][isubbtag]->Fill((tjet_jmr_dn+bjet_jmr_dn).M(),weight1);
				 hist_topjetsdmass_sel_md_DeepAK8_JMR_up[imdttag][ibtag][isubbtag]->Fill(FatJet_msoftdrop_jmr_up[topjetAK8],weight1);
				 hist_topjetsdmass_sel_md_DeepAK8_JMR_dn[imdttag][ibtag][isubbtag]->Fill(FatJet_msoftdrop_jmr_dn[topjetAK8],weight1);
				 hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JMR_up[ibtag][isubbtag]->Fill(FatJet_msoftdrop_jmr_up[topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight1);
				 hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JMR_dn[ibtag][isubbtag]->Fill(FatJet_msoftdrop_jmr_dn[topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight1);
			     
				 hist_tbmass_md_DeepAK8_noptw[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_noptw);
				 
				 hist_tbmass_md_DeepAK8_alphaup[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_alphaup);
				 hist_tbmass_md_DeepAK8_alphadn[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_alphadn);
				 hist_tbmass_md_DeepAK8_betaup[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_betaup);
				 hist_tbmass_md_DeepAK8_betadn[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_betadn);
				 
				
				hist_tbmass_md_DeepAK8_dAK8up[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_dAK8up);
				hist_tbmass_md_DeepAK8_dAK8dn[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight_dAK8dn);	
				hist_topjetsdmass_sel_md_DeepAK8up[imdttag][ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],weight_dAK8up);
				hist_topjetsdmass_sel_md_DeepAK8dn[imdttag][ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],weight_dAK8up);
				hist_topjetdeepmdtopscore_sel_md_DeepAK8up[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTagMD_TvsQCD[topjetAK8],weight_dAK8up);
				hist_topjetdeepmdtopscore_sel_md_DeepAK8dn[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTagMD_TvsQCD[topjetAK8],weight_dAK8dn);
					
				}
				
				
				for(int ijes=0; ijes<njesmax; ijes++){
					 
					 hist_tbmass_md_DeepAK8_JES_up[imdttag][ibtag][isubbtag][ijes]->Fill((tjet_jes_up[ijes]+bjet_jes_up[ijes]).M(),weight1);
					 hist_tbmass_md_DeepAK8_JES_dn[imdttag][ibtag][isubbtag][ijes]->Fill((tjet_jes_dn[ijes]+bjet_jes_dn[ijes]).M(),weight1);
					 hist_topjetsdmass_sel_md_DeepAK8_JES_up[imdttag][ibtag][isubbtag][ijes]->Fill(FatJet_msoftdrop_jes_up[ijes][topjetAK8],weight1);
					 hist_topjetsdmass_sel_md_DeepAK8_JES_dn[imdttag][ibtag][isubbtag][ijes]->Fill(FatJet_msoftdrop_jes_dn[ijes][topjetAK8],weight1);
					 hist_topjetpt_sel_md_DeepAK8_JES_up[imdttag][ibtag][isubbtag][ijes]->Fill(FatJet_pt_jes_up[ijes][topjetAK8],weight1);
					 hist_topjetpt_sel_md_DeepAK8_JES_dn[imdttag][ibtag][isubbtag][ijes]->Fill(FatJet_pt_jes_dn[ijes][topjetAK8],weight1);
					 hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JES_up[ibtag][isubbtag][ijes]->Fill(FatJet_msoftdrop_jes_up[ijes][topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight1);
					 hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JES_dn[ibtag][isubbtag][ijes]->Fill(FatJet_msoftdrop_jes_dn[ijes][topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight1);
					 
					 hist_bjetpt_sel_md_DeepAK8_JES_up[imdttag][ibtag][isubbtag][ijes]->Fill(Jet_pt_jes_up[ijes][bjetAK4],weight1);
					 hist_bjetpt_sel_md_DeepAK8_JES_dn[imdttag][ibtag][isubbtag][ijes]->Fill(Jet_pt_jes_dn[ijes][bjetAK4],weight1);
					 if(ieta>=0){
						hist_2D_pt_btag_AK4_md_DeepAK8_beta_JES_up[isubbtag][imdttag][ieta][ijes]->Fill(Jet_pt_jes_up[ijes][bjetAK4],Jet_btagDeepFlavB[bjetAK4],weight1);
						hist_2D_pt_btag_AK4_md_DeepAK8_beta_JES_dn[isubbtag][imdttag][ieta][ijes]->Fill(Jet_pt_jes_dn[ijes][bjetAK4],Jet_btagDeepFlavB[bjetAK4],weight1);
						}
					
					hist_tbmass_md_DeepAK8_JES_up_2[imdttag][ibtag][isubbtag][ijes]->Fill((tjet_jes_up[ijes]+bjet_jes_up[ijes]).M(),weight3);
					hist_tbmass_md_DeepAK8_JES_dn_2[imdttag][ibtag][isubbtag][ijes]->Fill((tjet_jes_dn[ijes]+bjet_jes_dn[ijes]).M(),weight3);
					
					}
					 
					 hist_tbmass_md_DeepAK8_JMS_up[imdttag][ibtag][isubbtag]->Fill((tjet_jms_up+bjet_jms_up).M(),weight1);
					 hist_tbmass_md_DeepAK8_JMS_dn[imdttag][ibtag][isubbtag]->Fill((tjet_jms_dn+bjet_jms_dn).M(),weight1);
					 hist_topjetsdmass_sel_md_DeepAK8_JMS_up[imdttag][ibtag][isubbtag]->Fill(FatJet_msoftdrop_jms_up[topjetAK8],weight1);
					 hist_topjetsdmass_sel_md_DeepAK8_JMS_dn[imdttag][ibtag][isubbtag]->Fill(FatJet_msoftdrop_jms_dn[topjetAK8],weight1);
					 hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JMS_up[ibtag][isubbtag]->Fill(FatJet_msoftdrop_jms_up[topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight1);
					 hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_JMS_dn[ibtag][isubbtag]->Fill(FatJet_msoftdrop_jms_dn[topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight1);
						 
			

				hist_bjetpt_sel_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(Jet_pt[bjetAK4],weight1);
				if(ieta>=0){
				hist_bjetpt_sel_md_DeepAK8_beta[imdttag][ibtag][isubbtag][ieta]->Fill(Jet_pt[bjetAK4],weight1);
				hist_tbmass_md_DeepAK8_eta1[imdttag][ibtag][isubbtag][ieta]->Fill(m_tb,weight1);
				hist_tbmass_md_DeepAK8_eta2[imdttag][ibtag][isubbtag][ieta]->Fill(m_tb,weight3);
				}
				if(isMC){				
  					hist_partonflav_AK4_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(Jet_partonFlavour[bjetAK4],weight1);
				}
				if((tjet+bjet).M() > 2000){
					hist_deleta_AK4_toptagAK8_dRpass_mtb2000->Fill((bjet.Rapidity()-tjet.Rapidity()),weight1);
				}
				#ifdef Btagger_DeepJet
					hist_2D_pt_btag_AK4_md_DeepAK8[isubbtag][imdttag]->Fill(Jet_pt[bjetAK4],Jet_btagDeepFlavB[bjetAK4],weight1);
					if(ieta>=0){
						hist_2D_pt_btag_AK4_md_DeepAK8_beta[isubbtag][imdttag][ieta]->Fill(Jet_pt[bjetAK4],Jet_btagDeepFlavB[bjetAK4],weight1);
						hist_2D_pt_btag_AK4_md_DeepAK8_beta2[isubbtag][imdttag][ieta]->Fill(Jet_pt[bjetAK4],Jet_btagDeepFlavB[bjetAK4],weight3);
						hist_2D_mtb_btag_AK4_md_DeepAK8_beta[isubbtag][imdttag][ieta]->Fill(m_tb,Jet_btagDeepFlavB[bjetAK4],weight1);
						}
					hist_bjetbtag_sel_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(Jet_btagDeepFlavB[bjetAK4],weight1);	
				#endif
				
				#ifdef Btagger_DeepCSV
					hist_2D_pt_btag_AK4_md_DeepAK8[isubbtag][imdttag]->Fill(Jet_pt[bjetAK4],Jet_btagDeepB[bjetAK4],weight1);
					if(ieta>=0){
						hist_2D_pt_btag_AK4_md_DeepAK8_beta[isubbtag][imdttag][ieta]->Fill(Jet_pt[bjetAK4],Jet_btagDeepB[bjetAK4],weight1);
						}
					hist_bjetbtag_sel_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(Jet_btagDeepB[bjetAK4],weight1);	
				#endif
				
				#ifdef Btagger_CSVv2
					hist_2D_pt_btag_AK4_md_DeepAK8[isubbtag][imdttag]->Fill(Jet_pt[bjetAK4],Jet_btagCSVV2[bjetAK4],weight1);
					if(ieta>=0){
						hist_2D_pt_btag_AK4_md_DeepAK8_beta[isubbtag][imdttag][ieta]->Fill(Jet_pt[bjetAK4],Jet_btagCSVV2[bjetAK4],weight1);
						}
					hist_bjetbtag_sel_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(Jet_btagCSVV2[bjetAK4],weight1);
				#endif
				
				hist_bjetmass_sel_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(bjet.M(),weight1);
						
				if(Jet_MatchFatJet[bjetAK4] >= 0){
					hist_bjetmatch_AK8Topscore_sel_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTag_TvsQCD[Jet_MatchFatJet[bjetAK4]],weight1);
					hist_bjetmatch_AK8Wscore_sel_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTag_WvsQCD[Jet_MatchFatJet[bjetAK4]],weight1);
					hist_bjetmatch_sdmass_sel_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(FatJet_msoftdrop[Jet_MatchFatJet[bjetAK4]],weight1);
					}
				
				hist_ht_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(Event_Ht,weight1);
					
				hist_tbmass_topjetpt_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(FatJet_pt[topjetAK8],(tjet+bjet).M(),weight1);
				hist_tbmass_topjeteta_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(FatJet_eta[topjetAK8],(tjet+bjet).M(),weight1);
				hist_tbmass_bjetpt_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(Jet_pt[bjetAK4],(tjet+bjet).M(),weight1);
				hist_tbmass_bjeteta_md_DeepAK8[imdttag][ibtag][isubbtag]->Fill(Jet_eta[bjetAK4],(tjet+bjet).M(),weight1);	
					
				if(1>0){
					
					hist_tbmass_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight2);
					hist_tbpt_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).Pt(),weight2);
					hist_tbrap_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).Rapidity(),weight2);
					if(isMC){				
						hist_partonflav_AK4_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill(Jet_partonFlavour[bjetAK4],weight2);
					}
					hist_2D_sdmass_tbmass_mddeepAK8_3[ibtag][deeptag]->Fill(FatJet_msoftdrop[topjetAK8],(tjet+bjet).M(),weight2);
					
					hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8_3[ibtag][deeptag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_subbtagCSVV2[topjetAK8],weight2);
					hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8_3[ibtag][deeptag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_subbtagDeepB[topjetAK8],weight2);
				
					#ifdef Btagger_DeepJet
						hist_2D_pt_btag_AK4_md_DeepAK8_3[isubbtag][imdttag]->Fill(Jet_pt[bjetAK4],Jet_btagDeepFlavB[bjetAK4],weight2);
					#endif
					#ifdef Btagger_DeepCSV
						hist_2D_pt_btag_AK4_md_DeepAK8_3[isubbtag][imdttag]->Fill(Jet_pt[bjetAK4],Jet_btagDeepB[bjetAK4],weight2);
					#endif
					#ifdef Btagger_CSVv2
						hist_2D_pt_btag_AK4_md_DeepAK8_3[isubbtag][imdttag]->Fill(Jet_pt[bjetAK4],Jet_btagCSVV2[bjetAK4],weight2);
					#endif
				
					hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8_3[ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_deepTag_TvsQCD[topjetAK8],weight2);
					hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_3[ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight2);
		 
					hist_topjetpt_sel_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill(FatJet_pt[topjetAK8],weight2);
					hist_topjetsdmass_sel_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],weight2);
					hist_topjetdeeptopscore_sel_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTag_TvsQCD[topjetAK8],weight2);
					hist_topjetdeepmdtopscore_sel_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTagMD_TvsQCD[topjetAK8],weight2);
					hist_topjetsbtag_sel_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill(FatJet_subbtagCSVV2[topjetAK8],weight2);
					hist_topjetsdeepCSV_sel_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill(FatJet_subbtagDeepB[topjetAK8],weight2);
					hist_topjettau32_sel_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill(tau32,weight2);
		  
					hist_bjetpt_sel_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill(Jet_pt[bjetAK4],weight2);
					hist_bjetmass_sel_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill(bjet.M(),weight2);
					#ifdef Btagger_DeepJet
						hist_bjetbtag_sel_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill(Jet_btagDeepFlavB[bjetAK4],weight2);
					#endif
					#ifdef Btagger_DeepCSV
						hist_bjetbtag_sel_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill(Jet_btagDeepB[bjetAK4],weight2);
					#endif
					#ifdef Btagger_CSVv2
						hist_bjetbtag_sel_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill(Jet_btagCSVV2[bjetAK4],weight2);
					#endif
						
					if(Jet_MatchFatJet[bjetAK4] >= 0){
						hist_bjetmatch_AK8Topscore_sel_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTag_TvsQCD[Jet_MatchFatJet[bjetAK4]],weight2);
						hist_bjetmatch_AK8Wscore_sel_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTag_WvsQCD[Jet_MatchFatJet[bjetAK4]],weight2);
						hist_bjetmatch_sdmass_sel_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill(FatJet_msoftdrop[Jet_MatchFatJet[bjetAK4]],weight2);
					}
					
					
					hist_ht_md_DeepAK8_3[imdttag][ibtag][isubbtag]->Fill(Event_Ht,weight2);
					
				   }//isFakeAnalysis
					
				} // cut on AK8 jet matched to b candidate
				
			if((Jet_MatchFatJet[bjetAK4]>=0) && !(FatJet_msoftdrop[Jet_MatchFatJet[bjetAK4]] < bmasscut_fun(FatJet_pt[Jet_MatchFatJet[bjetAK4]]))){
			if((Jet_MatchFatJet[bjetAK4]>=0) && (FatJet_msoftdrop[Jet_MatchFatJet[bjetAK4]]>topmasslow && FatJet_msoftdrop[Jet_MatchFatJet[bjetAK4]]<topmasshigh)){
		  
				hist_tbmass_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight1);
				if(FatJet_deepTagMD_TvsQCD[Jet_MatchFatJet[bjetAK4]] > deepak8_cut_md){
					hist_tbmass_md_DeepAK8_tt_deepak8pass[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight1);
					}
				hist_tbpt_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).Pt(),weight1);
				hist_tbrap_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).Rapidity(),weight1);
				if(isMC){				
  					hist_partonflav_AK4_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill(Jet_partonFlavour[bjetAK4],weight1);
				}
				
				int imsdbin = getbinid(FatJet_msoftdrop[topjetAK8],ntopsdmassbins,topsdmassbins);
				if(imsdbin>=0){
					hist_topjetpt_msdbin_sel_md_DeepAK8_tt[imdttag][ibtag][isubbtag][imsdbin]->Fill(FatJet_pt[topjetAK8],weight1);
					}
				
				hist_2D_sdmass_tbmass_mddeepAK8_tt[ibtag][deeptag]->Fill(FatJet_msoftdrop[topjetAK8],(tjet+bjet).M(),weight1);
					
				hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8_tt[ibtag][deeptag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_subbtagCSVV2[topjetAK8],weight1);
				hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8_tt[ibtag][deeptag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_subbtagDeepB[topjetAK8],weight1);
				
				int ieta = getbinid(fabs(Jet_eta[bjetAK4]),netabins,etabins);
				
				#ifdef Btagger_DeepJet
					hist_2D_pt_btag_AK4_md_DeepAK8_tt[isubbtag][imdttag]->Fill(Jet_pt[bjetAK4],Jet_btagDeepFlavB[bjetAK4],weight1);
					if(ieta>=0){
						hist_2D_pt_btag_AK4_md_DeepAK8_beta_tt[isubbtag][imdttag][ieta]->Fill(Jet_pt[bjetAK4],Jet_btagDeepFlavB[bjetAK4],weight1);
						}
				#endif
				#ifdef Btagger_DeepCSV
					hist_2D_pt_btag_AK4_md_DeepAK8_tt[isubbtag][imdttag]->Fill(Jet_pt[bjetAK4],Jet_btagDeepB[bjetAK4],weight1);
					if(ieta>=0){
						hist_2D_pt_btag_AK4_md_DeepAK8_beta_tt[isubbtag][imdttag][ieta]->Fill(Jet_pt[bjetAK4],Jet_btagDeepB[bjetAK4],weight1);
						}
				#endif
				#ifdef Btagger_CSVv2
					hist_2D_pt_btag_AK4_md_DeepAK8_tt[isubbtag][imdttag]->Fill(Jet_pt[bjetAK4],Jet_btagCSVV2[bjetAK4],weight1);
					if(ieta>=0){
						hist_2D_pt_btag_AK4_md_DeepAK8_beta_tt[isubbtag][imdttag][ieta]->Fill(Jet_pt[bjetAK4],Jet_btagCSVV2[bjetAK4],weight1);
						}
				#endif
				
				hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8_tt[ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_deepTag_TvsQCD[topjetAK8],weight1);
				hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_tt[ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight1);
		 
				int iptbin = getbinid(FatJet_pt[topjetAK8],ntopptbins,topptbins);
				hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_ptbin_tt[ibtag][isubbtag][iptbin]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight1);
				hist_tbmass_md_DeepAK8_ptbin_tt[imdttag][ibtag][isubbtag][iptbin]->Fill((tjet+bjet).M(),weight1);
		 
				hist_topjetpt_sel_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill(FatJet_pt[topjetAK8],weight1);
				hist_topjetsdmass_sel_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],weight1);
				hist_topjetdeeptopscore_sel_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTag_TvsQCD[topjetAK8],weight1);
				hist_topjetdeepmdtopscore_sel_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTagMD_TvsQCD[topjetAK8],weight1);
				hist_topjetsbtag_sel_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill(FatJet_subbtagCSVV2[topjetAK8],weight1);
				hist_topjetsdeepCSV_sel_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill(FatJet_subbtagDeepB[topjetAK8],weight1);
				hist_topjettau32_sel_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill(tau32,weight1);
				
				if(FatJet_deepTagMD_TvsQCD[Jet_MatchFatJet[bjetAK4]] > deepak8_cut_md){
					hist_topjetpt_sel_md_DeepAK8_tt_deepak8pass[imdttag][ibtag][isubbtag]->Fill(FatJet_pt[topjetAK8],weight1);
					hist_topjetsdmass_sel_md_DeepAK8_tt_deepak8pass[imdttag][ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],weight1);
				}
		  
				hist_bjetpt_sel_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill(Jet_pt[bjetAK4],weight1);
				hist_bjetmass_sel_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill(bjet.M(),weight1);
				#ifdef Btagger_DeepJet
					hist_bjetbtag_sel_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill(Jet_btagDeepFlavB[bjetAK4],weight1);
				#endif
				#ifdef Btagger_DeepCSV
					hist_bjetbtag_sel_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill(Jet_btagDeepB[bjetAK4],weight1);
				#endif
				#ifdef Btagger_CSVv2
					hist_bjetbtag_sel_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill(Jet_btagCSVV2[bjetAK4],weight1);
				#endif
						
				if(Jet_MatchFatJet[bjetAK4] >= 0){
					hist_bjetmatch_AK8Topscore_sel_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTag_TvsQCD[Jet_MatchFatJet[bjetAK4]],weight1);
					hist_bjetmatch_AK8Wscore_sel_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTag_WvsQCD[Jet_MatchFatJet[bjetAK4]],weight1);
					hist_bjetmatch_sdmass_sel_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill(FatJet_msoftdrop[Jet_MatchFatJet[bjetAK4]],weight1);
					}
					
				hist_ht_md_DeepAK8_tt[imdttag][ibtag][isubbtag]->Fill(Event_Ht,weight1);	
					
				if(1>0){
					hist_tbmass_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).M(),weight2);
					hist_tbpt_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).Pt(),weight2);
					hist_tbrap_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill((tjet+bjet).Rapidity(),weight2);
					if(isMC){				
						hist_partonflav_AK4_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill(Jet_partonFlavour[bjetAK4],weight2);
					}
					hist_2D_sdmass_tbmass_mddeepAK8_tt_3[ibtag][deeptag]->Fill(FatJet_msoftdrop[topjetAK8],(tjet+bjet).M(),weight2);
					
					hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8_tt_3[ibtag][deeptag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_subbtagCSVV2[topjetAK8],weight2);
					hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8_tt_3[ibtag][deeptag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_subbtagDeepB[topjetAK8],weight2);
				
					#ifdef Btagger_DeepJet
						hist_2D_pt_btag_AK4_md_DeepAK8_tt_3[isubbtag][imdttag]->Fill(Jet_pt[bjetAK4],Jet_btagDeepFlavB[bjetAK4],weight2);
					#endif
					#ifdef Btagger_DeepCSV
						hist_2D_pt_btag_AK4_md_DeepAK8_tt_3[isubbtag][imdttag]->Fill(Jet_pt[bjetAK4],Jet_btagDeepB[bjetAK4],weight2);
					#endif
					#ifdef Btagger_CSVv2
						hist_2D_pt_btag_AK4_md_DeepAK8_tt_3[isubbtag][imdttag]->Fill(Jet_pt[bjetAK4],Jet_btagCSVV2[bjetAK4],weight2);
					#endif
				
					hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8_tt_3[ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_deepTag_TvsQCD[topjetAK8],weight2);
					hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_tt_3[ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],FatJet_deepTagMD_TvsQCD[topjetAK8],weight2);
		 
					hist_topjetpt_sel_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill(FatJet_pt[topjetAK8],weight2);
					hist_topjetsdmass_sel_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill(FatJet_msoftdrop[topjetAK8],weight2);
					hist_topjetdeeptopscore_sel_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTag_TvsQCD[topjetAK8],weight2);
					hist_topjetdeepmdtopscore_sel_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTagMD_TvsQCD[topjetAK8],weight2);
					hist_topjetsbtag_sel_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill(FatJet_subbtagCSVV2[topjetAK8],weight2);
					hist_topjetsdeepCSV_sel_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill(FatJet_subbtagDeepB[topjetAK8],weight2);
					hist_topjettau32_sel_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill(tau32,weight2);
		  
					hist_bjetpt_sel_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill(Jet_pt[bjetAK4],weight2);
					hist_bjetmass_sel_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill(bjet.M(),weight2);
					#ifdef Btagger_DeepJet
						hist_bjetbtag_sel_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill(Jet_btagDeepFlavB[bjetAK4],weight2);
					#endif
					#ifdef Btagger_DeepCSV
						hist_bjetbtag_sel_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill(Jet_btagDeepB[bjetAK4],weight2);
					#endif
					#ifdef Btagger_CSVv2
						hist_bjetbtag_sel_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill(Jet_btagCSVV2[bjetAK4],weight2);
					#endif
						
					if(Jet_MatchFatJet[bjetAK4] >= 0){
						hist_bjetmatch_AK8Topscore_sel_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTag_TvsQCD[Jet_MatchFatJet[bjetAK4]],weight2);
						hist_bjetmatch_AK8Wscore_sel_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill(FatJet_deepTag_WvsQCD[Jet_MatchFatJet[bjetAK4]],weight2);
						hist_bjetmatch_sdmass_sel_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill(FatJet_msoftdrop[Jet_MatchFatJet[bjetAK4]],weight2);
					}
					
					hist_ht_md_DeepAK8_tt_3[imdttag][ibtag][isubbtag]->Fill(Event_Ht,weight2);
					
						}//isFakeAnalysis
					
					}// top mass cut on AK8 jet matched to b candidate (for t-tbar)
				} // cut on AK8 jet matched to b candidate (for t-tbar)
			  
			}  // b pt cut
		  } // whether all have logical WP 
		
		   
		} // if top & b candidates are present	    
	  	  
	  	}//entries  
	  	  
		fileIn->cd();
		delete T1;
		delete fileIn;
    }
    
	file_db.close();

    fileout->cd();
    fileout->Write();
    fileout->Close();

//	fp<<"Total Number of events in "<<fileOut->GetName()<<" is "<<nevent_total<<endl;
}
