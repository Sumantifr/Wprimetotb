#define Anal_Nano_Muon_cxx
//#define NoJEC
#define Btagger_DeepJet
//#define Btagger_DeepCSV
//#define Btagger_CSVv2
#define Btagger_DeepJet_wt
//#define Btagger_DeepCSV_wt

#define Anal_2017
//#define Data_2017B

//#define Anal_2016
//#define Data_2016H

// The class definition in Anal_Nano_Muon.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("Anal_Nano_Muon.C")
// Root > T->Process("Anal_Nano_Muon.C","some options")
// Root > T->Process("Anal_Nano_Muon.C+")
//

#include "Anal_Nano_Muon.h"
#include <TH2.h>
#include <TStyle.h>
#include <TH2.h>
#include <TStyle.h>
#include <TH2.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TH2F.h>
#include <fstream>
#include <TProofOutputFile.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <TProofServ.h>


void Anal_Nano_Muon::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void Anal_Nano_Muon::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   double pu_data[80] = {2.15411e-05,5.6992e-05,0.000127209,0.000284562,0.000628243,0.00130454,0.00248809,0.00434675,0.00699958,0.0104813,0.014724,0.0195614,0.0247524,0.0300163,0.0350706,0.0396636,0.0435963,0.0467326,0.0489994,0.0503793,0.0509,0.0506223,0.0496295,0.0480181,0.0458916,0.0433545,0.0405093,0.0374531,0.0342757,0.0310578,0.02787,0.0247719,0.0218122,0.0190285,0.0164478,0.0140877,0.011957,0.010057,0.00838287,0.0069247,0.00566889,0.00459925,0.00369801,0.00294673,0.00232704,0.0018212,0.00141253,0.00108574,0.000827064,0.000624371,0.000467131,0.000346366,0.000254531,0.000185382,0.000133825,9.57567e-05,6.79199e-05,4.77594e-05,3.3297e-05,2.30198e-05,1.57847e-05,1.07379e-05,7.2494e-06,4.85952e-06,3.23649e-06,2.14353e-06,1.41347e-06,9.29544e-07,6.11011e-07,4.02633e-07,2.66985e-07,1.78963e-07,1.21892e-07,8.47946e-08,6.05156e-08,4.4429e-08,3.35702e-08,2.6054e-08,2.06895e-08,1.67288e-08};

   double pu_MC[80] = {1.76003e-05,2.63752e-05,5.17764e-05,8.93703e-05,0.00010729,0.000139737,0.000241458,0.000725216,0.00129534,0.0024431,0.00503314,0.0092062,0.0146853,0.0204755,0.0267466,0.0337681,0.040181,0.0449554,0.0491183,0.0524606,0.0548204,0.0560606,0.0554619,0.0536665,0.0513695,0.0476493,0.0435073,0.0393303,0.0350598,0.0306642,0.0272167,0.0236864,0.0208034,0.0182563,0.0160933,0.0142631,0.012782,0.0115606,0.0105507,0.00957065,0.00885585,0.00826288,0.00758751,0.0069752,0.00622485,0.00547198,0.00483291,0.0040488,0.00337946,0.0027023,0.00212516,0.00160003,0.00117206,0.00085644,0.00056259,0.000367971,0.000246554,0.000160166,0.000101268,6.68301e-05,3.96657e-05,2.67149e-05,2.02967e-05,2.01132e-05,1.41779e-05,1.33955e-05,1.33216e-05,1.33078e-05,1.34937e-05,1.51125e-05,1.49794e-05,1.42329e-05,1.28154e-05,1.34488e-05,1.35528e-05,0,0,0,0,0};

   double logrhobins[norhobins+1] = {-0.088059,0.0942625,0.276584,0.458906,0.641227,0.823549,1.00587,1.18819,1.37051,1.55283,1.73516,1.91748,2.0998,2.28212,2.46444,2.64676,2.82909,3.01141,3.19373,3.37605,3.55837,3.74069,3.92302,4.10534,4.28766,4.46998,4.6523,4.83462,5.01694,5.19927,5.38159,5.56391,5.74623,5.92855,6.11087,6.2932,6.47552,6.65784,6.84016,7.02248,7.2048,7.38712,7.56945,7.75177,7.93409,8.11641,8.29873,8.48105,8.66338,8.8457,9.02802,9.21034};
   
   double sdmassbins[nsdmassbins+1] = {0,10,22,34,46,58,70,82,94,106,118,130,142,158,176,196,220,250,290,340,400}; 
    
   double tau32SF[3][4] = {{0.744729,0.782221,0.777628,0.810452},{0.511675,1.01478,1.147072,1.15288},{1.21932,0.931346,0.986976,0.999642}}; 
      

   OutFile = new TProofOutputFile("Output_trial.root");
 // fOutput->Add(OutFile);

    fileOut = OutFile->OpenFile("RECREATE");

    if ( !(fileOut = OutFile->OpenFile("RECREATE")) )
    {
      Warning("SlaveBegin", "problems opening file: %s/%s",
              OutFile->GetDir(), OutFile->GetFileName());
    }
   
   isMC = false;
   isQCD = false;
   FakeAnalysis = false;
   isSignal = false;
   isTTBar = false;
   isST = false;
   usePrefireWeight = true;
   
   Tout = new TTree("Tout", "WeightInfo");

   Tout->Branch("nevent_total", &nevent_total, "nevent_total/I");
   Tout->Branch("weightev", &weightev, "weightev/D");
   
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
   hist_btag_csv_bhadron[ipt] = new TH1D(name,title,100,-0.01,0.99);
   hist_btag_csv_bhadron[ipt]->Sumw2();
   
   sprintf(name,"BTagScore_DeepCSV_BHadron_ptbin%i",ipt+1);
   if(ipt<noperf_ptbins){	sprintf(title,"DeepCSV Score : p_{T} in [%i,%i] GeV",int(perf_ptbins[ipt]),int(perf_ptbins[ipt+1]));}
   else { sprintf(title,"DeepCSV Score : Inclusive in p_{T}"); }
   hist_btag_deepcsv_bhadron[ipt] = new TH1D(name,title,100,-0.01,0.99);
   hist_btag_deepcsv_bhadron[ipt]->Sumw2();
   
   sprintf(name,"BTagScore_DeepFlavour_BHadron_ptbin%i",ipt+1);
   if(ipt<noperf_ptbins){	sprintf(title,"DeepFlavour Score : p_{T} in [%i,%i] GeV",int(perf_ptbins[ipt]),int(perf_ptbins[ipt+1]));}
   else { sprintf(title,"DeepFlavour Score : Inclusive in p_{T}"); }
   hist_btag_deepflav_bhadron[ipt] = new TH1D(name,title,100,-0.01,0.99);
   hist_btag_deepflav_bhadron[ipt]->Sumw2();
   
   sprintf(name,"BTagScore_CSVv2_QHadron_ptbin%i",ipt+1);
   if(ipt<noperf_ptbins){	sprintf(title,"CSVv2 Score : p_{T} in [%i,%i] GeV",int(perf_ptbins[ipt]),int(perf_ptbins[ipt+1])); }
   else { sprintf(title,"CSVv2 Score : Inclusive in p_{T}"); }
   hist_btag_csv_qhadron[ipt] = new TH1D(name,title,100,-0.01,0.99);
   hist_btag_csv_qhadron[ipt]->Sumw2();
   
   sprintf(name,"BTagScore_DeepCSV_QHadron_ptbin%i",ipt+1);
   if(ipt<noperf_ptbins){	sprintf(title,"DeepCSV Score : p_{T} in [%i,%i] GeV",int(perf_ptbins[ipt]),int(perf_ptbins[ipt+1]));}
   else { sprintf(title,"DeepCSV Score : Inclusive in p_{T}"); }
   hist_btag_deepcsv_qhadron[ipt] = new TH1D(name,title,100,-0.01,0.99);
   hist_btag_deepcsv_qhadron[ipt]->Sumw2();
   
   sprintf(name,"BTagScore_DeepFlavour_QHadron_ptbin%i",ipt+1);
   if(ipt<noperf_ptbins){	sprintf(title,"DeepFlavour Score : p_{T} in [%i,%i] GeV",int(perf_ptbins[ipt]),int(perf_ptbins[ipt+1]));}
   else { sprintf(title,"DeepFlavour Score : Inclusive in p_{T}"); }
   hist_btag_deepflav_qhadron[ipt] = new TH1D(name,title,100,-0.01,0.99);
   hist_btag_deepflav_qhadron[ipt]->Sumw2();
   
   }
   
   sprintf(name,"N_PV");
   sprintf(title,"# of Primary Vertices");
   hist_npv = new TH1D(name,title,100,-0.1,99.9);//80,-0.1,79.9);
   hist_npv->Sumw2();
   
   sprintf(name,"N_PV_Final");
   sprintf(title,"# of Primary Vertices");
   hist_npv_final = new TH1D(name,title,100,-0.1,99.9);//80,-0.1,79.9);
   hist_npv->Sumw2();
   
   sprintf(name,"N_PU");
   sprintf(title,"# of Pile-up Vertices");
   hist_npu = new TH1D(name,title,100,0,100);//80,-0.1,79.9);
   hist_npu->Sumw2();
   
   sprintf(name,"MET_Filter_Passed");
   sprintf(title,"MET Filter Passing Score");
   hist_metfilter_pass= new TH1D(name,title,2,-0.1,1.9);
   hist_metfilter_pass->Sumw2();
   
   sprintf(name,"Mu_Trigger_Passed");
   sprintf(title,"Muon Trigger Passing Score");
   hist_mutrig_pass= new TH1D(name,title,2,-0.1,1.9);
   hist_mutrig_pass->Sumw2();
   
   sprintf(name,"Elec_Trigger_Passed");
   sprintf(title,"Electron Trigger Passing Score");
   hist_etrig_pass= new TH1D(name,title,2,-0.1,1.9);
   hist_etrig_pass->Sumw2();
   
   sprintf(name,"Pho_Trigger_Passed");
   sprintf(title,"HLT_AK8CHS350_TrimMass30 Trigger Passing Score");
   hist_photrig_pass= new TH1D(name,title,2,-0.1,1.9);
   hist_photrig_pass->Sumw2();
   
   sprintf(name,"Had_Trigger_Passed");
   sprintf(title,"Hadronic Trigger Passing Score");
   hist_hadtrig_pass= new TH1D(name,title,2,-0.1,1.9);
   hist_hadtrig_pass->Sumw2();
   
  
   sprintf(name,"N_Muons");
   sprintf(title,"# of Muons");
   hist_nmuons= new TH1D(name,title,5,-0.1,4.9);
   hist_nmuons->Sumw2();
   
   sprintf(name,"N_Electrons");
   sprintf(title,"# of Electrons");
   hist_nelectrons= new TH1D(name,title,5,-0.1,4.9);
   hist_nelectrons->Sumw2();
   
   sprintf(name,"N_Photons");
   sprintf(title,"# of Photons");
   hist_nphotons= new TH1D(name,title,5,-0.1,4.9);
   hist_nphotons->Sumw2();
   
   sprintf(name,"N_Muons_Cut");
   sprintf(title,"# of Muons (N-1 cut)");
   hist_nmuons_cut= new TH1D(name,title,5,-0.1,4.9);
   hist_nmuons_cut->Sumw2();
   
   sprintf(name,"N_Electrons_Cut");
   sprintf(title,"# of Electrons (N-1 cut)");
   hist_nelectrons_cut= new TH1D(name,title,5,-0.1,4.9);
   hist_nelectrons_cut->Sumw2();
   
   sprintf(name,"N_Photons_Cut");
   sprintf(title,"# of Photons (N-1 cut)");
   hist_nphotons_cut= new TH1D(name,title,5,-0.1,4.9);
   hist_nphotons_cut->Sumw2();
   
   sprintf(name,"Mu_Trigger_Passed_Cut");
   sprintf(title,"Muon Trigger Passing Score (N-1 cut)");
   hist_mutrig_pass_cut= new TH1D(name,title,2,-0.1,1.9);
   hist_mutrig_pass_cut->Sumw2();
   
   sprintf(name,"Trig_Match_Passed");
   sprintf(title,"Trigger Matching Passing Score");
   hist_trig_match= new TH1D(name,title,2,-0.1,1.9);
   hist_trig_match->Sumw2();
   
   sprintf(name,"Trig_Match_Passed_Cut");
   sprintf(title,"Trigger Matching Passing Score (N-1 cut)");
   hist_trig_match_cut= new TH1D(name,title,2,-0.1,1.9);
   hist_trig_match_cut->Sumw2();
  
   sprintf(name,"MET_corpt");
   sprintf(title,"Corrected MET P_{T}");
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
   sprintf(title,"M (GeV) of Top Parton");
   hist_mass_hadtopq = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_mass_hadtopq->Sumw2();
   
   sprintf(name,"PID_HardParton_TopDaughter");
   sprintf(title,"PID of Hard Partons Top Daughter");
   hist_pid_hadtopq  = new TH1D(name,title,24,-0.1,23.9);
   hist_pid_hadtopq ->Sumw2();
   
   sprintf(name,"Match_AK8Jet_TopParton");
   sprintf(title,"Match AK8Jet of Top Parton");
   hist_matchjet_hadtopq = new TH1D(name,title,20,-1.1,19.9);
   hist_matchjet_hadtopq->Sumw2();
   
   for(int ipt=0; ipt<noperf_ptbins; ipt++){
   
   sprintf(name,"SDMass_TopMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"M_{SD} (GeV) of Top Matched AK8 Jet");
   hist_jetsdmass_top[ipt] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetsdmass_top[ipt]->Sumw2();
   
   sprintf(name,"Tau32_TopMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"#tau_{3} / #tau_{2} of Top Matched AK8 Jet");
   hist_jettau32_top[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jettau32_top[ipt]->Sumw2();
   
   sprintf(name,"Subjetbtag_CSVv2_TopMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"Subjet B-tag Value (CSVv2) of Top Matched AK8 Jet");
   hist_jetsubbtag_CSVv2_top[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jetsubbtag_CSVv2_top[ipt]->Sumw2();
   
   sprintf(name,"Subjetbtag_DeepCSV_TopMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"Subjet B-tag Value (DeepCSV) of Top Matched AK8 Jet");
   hist_jetsubbtag_DeepCSV_top[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jetsubbtag_DeepCSV_top[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_MD_TopMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of Top Matched AK8 Jet");
   hist_jetDeepAK8_MD_top[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jetDeepAK8_MD_top[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_MD_withMasscut_TopMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of Top Matched AK8 Jet (m_{SD}:(105,220) GeV");
   hist_jetDeepAK8_MD_wm_top[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jetDeepAK8_MD_wm_top[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_TopMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score (Mass Decorrelated) of Top Matched AK8 Jet");
   hist_jetDeepAK8_top[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jetDeepAK8_top[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_withMasscut_TopMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score (Mass Decorrelated) of Top Matched AK8 Jet (m_{SD}:(105,220) GeV");
   hist_jetDeepAK8_wm_top[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jetDeepAK8_wm_top[ipt]->Sumw2();
   
   sprintf(name,"GenJetMass_TopMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"Gen Jet Mass (GeV) of Top Matched AK8 Jet");
   hist_jetgenmass_top[ipt] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetgenmass_top[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_TopMatched_AK8Jet_tau32_1_subDeepCSV_1_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of Top Matched AK8 Jet : #tau_{3} / #tau_{2} < xx && Subjet DeepCSV score > yy");
   hist_jetDeepAK8_top_tausub_11[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jetDeepAK8_top_tausub_11[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_TopMatched_AK8Jet_tau32_1_subDeepCSV_2_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of Top Matched AK8 Jet : #tau_{3} / #tau_{2} < xx && Subjet DeepCSV score > yy");
   hist_jetDeepAK8_top_tausub_12[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jetDeepAK8_top_tausub_12[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_TopMatched_AK8Jet_tau32_2_subDeepCSV_1_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of Top Matched AK8 Jet : #tau_{3} / #tau_{2} < xx && Subjet DeepCSV score > yy");
   hist_jetDeepAK8_top_tausub_21[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jetDeepAK8_top_tausub_21[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_TopMatched_AK8Jet_tau32_2_subDeepCSV_2_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of Top Matched AK8 Jet : #tau_{3} / #tau_{2} < xx && Subjet DeepCSV score > yy");
   hist_jetDeepAK8_top_tausub_22[ipt] = new TH1D(name,title,50,-0.01,0.99);
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
   sprintf(title,"M_{SD} (GeV) of q/g Matched AK8 Jet");
   hist_jetsdmass_tbkg[ipt] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetsdmass_tbkg[ipt]->Sumw2();
   
   sprintf(name,"Tau32_QMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"#tau_{3} / #tau_{2} of q/g Matched AK8 Jet");
   hist_jettau32_tbkg[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jettau32_tbkg[ipt]->Sumw2();
   
   sprintf(name,"Subjetbtag_CSVv2_QMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"Subjet B-tag Value (CSVv2) of q/g Matched AK8 Jet");
   hist_jetsubbtag_CSVv2_tbkg[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jetsubbtag_CSVv2_tbkg[ipt]->Sumw2();
   
   sprintf(name,"Subjetbtag_DeepCSV_QMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"Subjet B-tag Value (DeepCSV) of q/g Matched AK8 Jet");
   hist_jetsubbtag_DeepCSV_tbkg[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jetsubbtag_DeepCSV_tbkg[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_MD_QMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of q/g Matched AK8 Jet");
   hist_jetDeepAK8_MD_tbkg[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jetDeepAK8_MD_tbkg[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_MD_withMasscut_QMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of q/g Matched AK8 Jet (m_{SD}:(105,220) GeV");
   hist_jetDeepAK8_MD_wm_tbkg[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jetDeepAK8_MD_wm_tbkg[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_QMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score (Mass Decorrelated) of q/g Matched AK8 Jet");
   hist_jetDeepAK8_tbkg[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jetDeepAK8_tbkg[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_withMasscut_QMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of q/g Matched AK8 Jet (m_{SD}:(105,220) GeV");
   hist_jetDeepAK8_wm_tbkg[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jetDeepAK8_wm_tbkg[ipt]->Sumw2();
   
   sprintf(name,"GenJetMass_QMatched_AK8Jet_pt%i",ipt+1);
   sprintf(title,"Gen Jet Mass (GeV) of q/g Matched AK8 Jet");
   hist_jetgenmass_tbkg[ipt] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetgenmass_tbkg[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_QMatched_AK8Jet_tau32_1_subDeepCSV_1_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of q/g Matched AK8 Jet : #tau_{3} / #tau_{2} < xx && Subjet DeepCSV score > yy");
   hist_jetDeepAK8_tbkg_tausub_11[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jetDeepAK8_tbkg_tausub_11[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_QMatched_AK8Jet_tau32_1_subDeepCSV_2_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of q/g Matched AK8 Jet : #tau_{3} / #tau_{2} < xx && Subjet DeepCSV score > yy");
   hist_jetDeepAK8_tbkg_tausub_12[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jetDeepAK8_tbkg_tausub_12[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_QMatched_AK8Jet_tau32_2_subDeepCSV_1_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of q/g Matched AK8 Jet : #tau_{3} / #tau_{2} < xx && Subjet DeepCSV score > yy");
   hist_jetDeepAK8_tbkg_tausub_21[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jetDeepAK8_tbkg_tausub_21[ipt]->Sumw2();
   
   sprintf(name,"DeepAK8Tag_QMatched_AK8Jet_tau32_2_subDeepCSV_2_pt%i",ipt+1);
   sprintf(title,"DeepAK8 Score of q/g Matched AK8 Jet : #tau_{3} / #tau_{2} < xx && Subjet DeepCSV score > yy");
   hist_jetDeepAK8_tbkg_tausub_22[ipt] = new TH1D(name,title,50,-0.01,0.99);
   hist_jetDeepAK8_tbkg_tausub_22[ipt]->Sumw2();
   
	}
   
   // end of top tagging histograms //
   
   sprintf(name,"NJets_AK8CHS");
   sprintf(title,"# of AK8CHS Jets");
   hist_njetAK8 = new TH1D(name,title,10,-0.1,9.9);
   hist_njetAK8->Sumw2();
   
   sprintf(name,"Pt_AK8CHS");
   sprintf(title,"P_{T} of AK8CHS Jets");
   hist_jetptAK8 = new TH1D(name,title,noptbins,ptbins);
   hist_jetptAK8->Sumw2();
   
   sprintf(name,"Mass_AK8CHS");
   sprintf(title,"Mass of AK8CHS Jets");
   hist_jetmassAK8 = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetmassAK8->Sumw2();
   
   sprintf(name,"SDMass_AK8CHS");
   sprintf(title,"Soft-Drop Mass of AK8CHS Jets");
   hist_jetsdmassAK8 = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetsdmassAK8->Sumw2();
   
   sprintf(name,"y_AK8CHS");
   sprintf(title,"Rapidity of AK8CHS Jets");
   hist_jetrapAK8 = new TH1D(name,title,50,-5.,5.);
   hist_jetrapAK8->Sumw2();
  
   sprintf(name,"rho_AK8CHS");
   sprintf(title,"#rho (-2*log(m_{SD}/p_{T})) of AK8CHS Jets");
   hist_jetrhoAK8 = new TH1D(name,title,norhobins,logrhobins);
   hist_jetrhoAK8->Sumw2();
   
   sprintf(name,"tau32_AK8CHS");
   sprintf(title,"#tau_3 / #tau_2 of AK8CHS Jets");
   hist_jettau32AK8 = new TH1D(name,title,50,-0.01,0.99);
   hist_jettau32AK8->Sumw2();
   
   sprintf(name,"tau21_AK8CHS");
   sprintf(title,"#tau_2 / #tau_1 of AK8CHS Jets");
   hist_jettau21AK8 = new TH1D(name,title,50,-0.1,0.99);
   hist_jettau21AK8->Sumw2();
   
   sprintf(name,"N_Subjets_passing_medium_WP");
   sprintf(title,"# of Subjets Passing Medium WP");
   hist_njetsubjets_AK8_mpass = new TH1D(name,title,3,-0.1,2.9);
   hist_njetsubjets_AK8_mpass->Sumw2();
   
   sprintf(name,"N_Subjets_passing_loose_WP_fail_medium_WP");
   sprintf(title,"# of Subjets Passing Loose WP && Failing Medium WP");
   hist_njetsubjets_AK8_lpass_mfail = new TH1D(name,title,3,-0.1,2.9);
   hist_njetsubjets_AK8_lpass_mfail->Sumw2();
   
   sprintf(name,"subjetbtag_AK8CHS");
   sprintf(title,"Suubjet CSVv2 Score of AK8CHS Jets");
   hist_jetsubbtagAK8 = new TH1D(name,title,20,-0.01,0.99);
   hist_jetsubbtagAK8->Sumw2();
   
   sprintf(name,"subjetmass_AK8CHS");
   sprintf(title,"Suubjet Mass of AK8CHS Jets");
   hist_jetsubmassAK8 = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetsubmassAK8->Sumw2();
   
   sprintf(name,"subjetcostheta_AK8CHS");
   sprintf(title,"cos(#Delta#theta^{tb_{t}})");
   hist_jetsubthetAK8 = new TH1D(name,title,50,-1.,+1.);
   hist_jetsubthetAK8->Sumw2();
   
   sprintf(name,"SDMass_AK8CHS_TopCand");
   sprintf(title,"Soft-Drop Mass of AK8CHS Top Candidate Jet");
   hist_topjetsdmass = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass->Sumw2();
   
   sprintf(name,"Tau32_AK8CHS_TopCand");
   sprintf(title,"#tau_{32} of AK8CHS Top Candidate Jet");
   hist_topjettau32AK8 = new TH1D(name,title,50,-0.01,0.99);
   hist_topjettau32AK8->Sumw2();
   
   sprintf(name,"DeepTagTvsQCD_AK8CHS_TopCand");
   sprintf(title,"DeepTagTvsQCD of AK8CHS Top Candidate Jet");
   hist_topjetdeeptopscore = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeeptopscore->Sumw2();
   
   sprintf(name,"DeepTagTvsQCD_MD_AK8CHS_TopCand");
   sprintf(title,"DeepTagTvsQCD (Mass Decorrelated) Score of AK8CHS Top Candidate Jet");
   hist_topjetmddeeptopscore = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetmddeeptopscore->Sumw2();
   
   sprintf(name,"Pt_AK8CHS_TopCand");
   sprintf(title,"P_{T} of AK8CHS Jets");
   hist_topjetpt= new TH1D(name,title,noptbins,ptbins);
   hist_topjetpt->Sumw2();
   
   sprintf(name,"Pt_AK8CHS_TopCand_PartonMatched");
   sprintf(title,"P_{T} of AK8CHS Jets (Matched to t parton)");
   hist_topjetptmatch= new TH1D(name,title,noptbins,ptbins);
   hist_topjetptmatch->Sumw2();
   
   sprintf(name,"Pt_AK8CHS_Top_Parton");
   sprintf(title,"P_{T} of top parton)");
   hist_toppartonpt= new TH1D(name,title,noptbins,ptbins);
   hist_toppartonpt->Sumw2();
   
   sprintf(name,"Pt_AK8CHS_TopCand_MatchedtoParton");
   sprintf(title,"P_{T} of AK8CHS Jets (Matched to t parton)");
   hist_topparton_matchedtopjet_pt= new TH1D(name,title,noptbins,ptbins);
   hist_topparton_matchedtopjet_pt->Sumw2();
   
   sprintf(name,"GenMass_AK8CHS_TopCand");
   sprintf(title,"Mass of GenJet Corresponding to AK8CHS Top Candidate Jet");
   hist_topjetgenmass = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetgenmass->Sumw2();
   
   sprintf(name,"GenPID_AK8CHS_TopCand");
   sprintf(title,"Parton Flavour of GenJet Corresponding to AK8CHS Top Candidate Jet");
   hist_topjetgenpid = new TH1D(name,title,22,-0.1,21.9);
   hist_topjetgenpid->Sumw2();
   
   sprintf(name,"H2D_SDMass_Tau32_AK8CHS");
   sprintf(title,"2D Correlation of m_{SD} and #tau_{32} of AK8CHS Jets");
   hist_2D_sdmass_tau32_AK8 = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
   hist_2D_sdmass_tau32_AK8->Sumw2();
   
   sprintf(name,"H2D_SDMass_subbtag_AK8CHS");
   sprintf(title,"2D Correlation of m_{SD} and Subjet b-tag AK8CHS Jets");
   hist_2D_sdmass_subbtag_AK8 = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
   hist_2D_sdmass_subbtag_AK8->Sumw2();
   
   sprintf(name,"H2D_Tau32_subbtag_AK8CHS");
   sprintf(title,"2D Correlation of #tau_{32} and Subjet b-tag of AK8CHS Jets");
   hist_2D_tau32_subbtag_AK8 = new TH2D(name,title,50,-0.01,0.99,20,-0.01,0.99);
   hist_2D_tau32_subbtag_AK8->Sumw2();
   
   sprintf(name,"TopCand_SubJetMass_12");
   sprintf(title,"2D Correlation of Subjet Mass of AK8 Jet t-Cand");
   hist_2D_topjet_subjetmass12 = new TH2D(name,title,nsdmassbins,sdmassbins,nsdmassbins,sdmassbins);
   hist_2D_topjet_subjetmass12->Sumw2();
   
   sprintf(name,"TopCand_SubJetMass1_vs_JetSDMass");
   sprintf(title,"2D Correlation of Higher Subjet Mass and SDMass AK8 Jet t-Cand");
   hist_2D_topjetsdmass_subjetmass1 = new TH2D(name,title,nsdmassbins,sdmassbins,nsdmassbins,sdmassbins);
   hist_2D_topjetsdmass_subjetmass1->Sumw2();
   
   sprintf(name,"TopCand_SubJetMass2_vs_JetSDMass");
   sprintf(title,"2D Correlation of Lower Subjet Mass and SDMass AK8 Jet t-Cand");
   hist_2D_topjetsdmass_subjetmass2 = new TH2D(name,title,nsdmassbins,sdmassbins,nsdmassbins,sdmassbins);
   hist_2D_topjetsdmass_subjetmass2->Sumw2();
   
   // m-pass //
   
   sprintf(name,"tau32_mpass_AK8CHS");
   sprintf(title,"#tau_3 / #tau_2 of AK8CHS Jets within (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_jettau32AK8_mpass = new TH1D(name,title,50,-0.01,0.99);
   hist_jettau32AK8_mpass->Sumw2();
   
   sprintf(name,"N_Subjets_passing_medium_WP_mpass");
   sprintf(title,"# of Subjets Passing Medium WP within (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_njetsubjets_AK8_BMPpass_mpass = new TH1D(name,title,3,-0.1,2.9);
   hist_njetsubjets_AK8_BMPpass_mpass->Sumw2();
   
   sprintf(name,"NBtagSubJets_mpass_AK8CHS");
   sprintf(title,"# of Btagged Subjets AK8CHS Jets within (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_nbsubjetAK8_mpass = new TH1D(name,title,3,-0.1,2.9);
   hist_nbsubjetAK8_mpass->Sumw2();
   
   sprintf(name,"N_Subjets_passing_loose_WP_fail_medium_WP_mpass");
   sprintf(title,"# of Subjets Passing Loose WP && Failing Medium WP within (%i,%i) GeV",int(topmasslow),int(topmasshigh));
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
   sprintf(title,"2D Correlation of #tau_{32} and Subjet b-tag of AK8CHS Jets  within (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_2D_tau32_subbtag_AK8_mpass = new TH2D(name,title,50,-0.01,0.99,20,-0.01,0.99);
   hist_2D_tau32_subbtag_AK8_mpass->Sumw2();
   
   // m-fail //
   
   sprintf(name,"tau32_mfail_AK8CHS");
   sprintf(title,"#tau_3 / #tau_2 of AK8CHS Jets outside (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_jettau32AK8_mfail = new TH1D(name,title,50,-0.01,0.99);
   hist_jettau32AK8_mfail->Sumw2();
   
   sprintf(name,"N_Subjets_passing_medium_WP_mfail");
   sprintf(title,"# of Subjets Passing Medium WP outside (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_njetsubjets_AK8_BMPpass_mfail = new TH1D(name,title,3,-0.1,2.9);
   hist_njetsubjets_AK8_BMPpass_mfail->Sumw2();
   
   sprintf(name,"NBtagSubJets_mfail_AK8CHS");
   sprintf(title,"# of Btagged Subjets AK8CHS Jets outside (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_nbsubjetAK8_mfail = new TH1D(name,title,3,-0.1,2.9);
   hist_nbsubjetAK8_mfail->Sumw2();
   
   sprintf(name,"N_Subjets_passing_loose_WP_fail_medium_WP_mfail");
   sprintf(title,"# of Subjets Passing Loose WP && Failing Medium WP outside (%i,%i) GeV",int(topmasslow),int(topmasshigh));
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
   hist_jetsubbtagAK8_mfail = new TH1D(name,title,20,-0.01,0.99);
   hist_jetsubbtagAK8_mfail->Sumw2();
   
   sprintf(name,"subjetbtag_deepCSV_mfail_AK8CHS");
   sprintf(title,"Suubjet DeepCSV Score of AK8CHS Jets outside (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_jetsubbtagAK8_deepCSV_mfail = new TH1D(name,title,20,-0.01,0.99);
   hist_jetsubbtagAK8_deepCSV_mfail->Sumw2();
   
   sprintf(name,"subjetpt_btagpass_mfail_AK8CHS");
   sprintf(title,"B-tagged Suubjet P_{T} of AK8CHS Jets outside (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_jetsubptAK8_btagpass_mfail = new TH1D(name,title,noptbins,ptbins);
   hist_jetsubptAK8_btagpass_mfail->Sumw2();
   
   sprintf(name,"H2D_Tau32_subbtag_mfail_AK8CHS");
   sprintf(title,"2D Correlation of #tau_{32} and Subjet b-tag of AK8CHS Jets  outside (%i,%i) GeV",int(topmasslow),int(topmasshigh));
   hist_2D_tau32_subbtag_AK8_mfail = new TH2D(name,title,50,-0.01,0.99,20,-0.01,0.99);
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
   sprintf(title,"2D Correlation of Subjet b-tag (DeepCSV) and DeepAK8 Toptag Score");
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
   sprintf(title,"2D Correlation of Subjet b-tag (DeepCSV) and Mass Decorrelated DeepAK8 Toptag Score");
   hist_2D_subjetbtag_deepmdtopscore = new TH2D(name,title,50,-0.001,0.999,50,-0.001,0.999);
   hist_2D_subjetbtag_deepmdtopscore->Sumw2();
   
		for(int ib=0; ib<nbtag; ib++){
		for(int itau=0; itau<ntautag; itau++){
    
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_Tau32Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		if(itau==0) { sprintf(title,"2D Correlation of Soft-Drop Mass  and Subjet b-tag (CSVv2) of AK8CHS Jets | #tau_{32} < %0.2f btag %i",tau32_cut,ib); }
		if(itau==1) { sprintf(title,"2D Correlation of Soft-Drop Mass  and Subjet b-tag (CSVv2) of AK8CHS Jets | 0.5 < #tau_{32} < %0.2f btag %i",tau32_cut_loose,ib); }
		if(itau==2) { sprintf(title,"2D Correlation of Soft-Drop Mass  and Subjet b-tag (CSVv2) of AK8CHS Jets | #tau_{32} > %0.2f btag %i",tau32_cut_loose,ib); }
		
		hist_2D_sdmass_subbtag_CSVv2_AK8[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_subbtag_CSVv2_AK8[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_Tau32Tag%i_AK4btag%i_AK8CHS_lowbpt",itau+1,ib+1);
		hist_2D_sdmass_subbtag_CSVv2_AK8_1[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_subbtag_CSVv2_AK8_1[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_Tau32Tag%i_AK4btag%i_AK8CHS_nobmasscut",itau+1,ib+1);
		hist_2D_sdmass_subbtag_CSVv2_AK8_2[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_subbtag_CSVv2_AK8_2[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_Tau32Tag%i_AK4btag%i_AK8CHS_highht",itau+1,ib+1);
		hist_2D_sdmass_subbtag_CSVv2_AK8_3[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_subbtag_CSVv2_AK8_3[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_DeepCSV_Tau32Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		if(itau==0) { sprintf(title,"2D Correlation of Soft-Drop Mass  and Subjet b-tag (DeepCSV) of AK8CHS Jets | #tau_{32} < %0.2f btag %i",tau32_cut,ib); }
		if(itau==1) { sprintf(title,"2D Correlation of Soft-Drop Mass  and Subjet b-tag (DeepCSV) of AK8CHS Jets | 0.5 < #tau_{32} < %0.2f btag %i",tau32_cut_loose,ib); }
		if(itau==2) { sprintf(title,"2D Correlation of Soft-Drop Mass  and Subjet b-tag (DeepCSV) of AK8CHS Jets | #tau_{32} > %0.2f btag %i",tau32_cut_loose,ib); }
		
		hist_2D_sdmass_subbtag_DeepCSV_AK8[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_subbtag_DeepCSV_AK8[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_DeepCSV_Tau32Tag%i_AK4btag%i_AK8CHS_lowbpt",itau+1,ib+1);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_1[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_1[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_DeepCSV_Tau32Tag%i_AK4btag%i_AK8CHS_nobmasscut",itau+1,ib+1);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_2[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_2[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_DeepCSV_Tau32Tag%i_AK4btag%i_AK8CHS_highht",itau+1,ib+1);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_3[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_3[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_DeepAK8Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		hist_2D_sdmass_subbtag_CSVv2_AK8_DeepAK8[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_subbtag_CSVv2_AK8_DeepAK8[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_DeepCSV_DeepAK8Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_DeepAK8[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_DeepAK8[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_Tau32_DeepAK8Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		hist_2D_sdmass_deeptag_AK8_reg_DeepAK8[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_deeptag_AK8_reg_DeepAK8[ib][itau]->Sumw2();
		
		
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_DeepAK8Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_DeepAK8Tag%i_AK4btag%i_AK8CHS_hight",itau+1,ib+1);
		hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8_3[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8_3[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_DeepCSV_MDDeepAK8Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_DeepCSV_MDDeepAK8Tag%i_AK4btag%i_AK8CHS_hight",itau+1,ib+1);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8_3[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8_3[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_Tau32_MDDeepAK8Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		hist_2D_sdmass_deeptag_AK8_reg_md_DeepAK8[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_deeptag_AK8_reg_md_DeepAK8[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_Tau32_MDDeepAK8Tag%i_AK4btag%i_AK8CHS_hight",itau+1,ib+1);
		hist_2D_sdmass_deeptag_AK8_reg_md_DeepAK8_3[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_deeptag_AK8_reg_md_DeepAK8_3[ib][itau]->Sumw2();
		
		// t-tbar (MD) //
		
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_DeepAK8Tag%i_AK4btag%i_AK8CHS_ttbar",itau+1,ib+1);
		hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8_tt[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8_tt[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_DeepAK8Tag%i_AK4btag%i_AK8CHS_ttbar_hight",itau+1,ib+1);
		hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8_tt_3[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8_tt_3[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_DeepCSV_MDDeepAK8Tag%i_AK4btag%i_AK8CHS_ttbar",itau+1,ib+1);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8_tt[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8_tt[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_DeepCSV_MDDeepAK8Tag%i_AK4btag%i_AK8CHS_ttbar_hight",itau+1,ib+1);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8_tt_3[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_subbtag_DeepCSV_AK8_md_DeepAK8_tt_3[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_Tau32_MDDeepAK8Tag%i_AK4btag%i_AK8CHS_ttbar",itau+1,ib+1);
		hist_2D_sdmass_deeptag_AK8_reg_md_DeepAK8_tt[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
		hist_2D_sdmass_deeptag_AK8_reg_md_DeepAK8_tt[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_Tau32_MDDeepAK8Tag%i_AK4btag%i_AK8CHS_ttbar_hight",itau+1,ib+1);
		hist_2D_sdmass_deeptag_AK8_reg_md_DeepAK8_tt_3[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,-0.01,0.99);
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
		
		sprintf(name,"H2D_SDMass_tbMass_DeepAK8Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		hist_2D_sdmass_tbmass_deepAK8[ib][itau] = new TH2D(name,name,nsdmassbins,sdmassbins,nvwpmbin,wpmbins);
		hist_2D_sdmass_tbmass_deepAK8[ib][itau]->Sumw2();
		
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
		hist_2D_DeepAK8_tbmass_1[ib] = new TH2D(name,name,20,-0.01,0.99,nvwpmbin,wpmbins);
		hist_2D_DeepAK8_tbmass_1[ib] = new TH2D(name,name,20,-0.01,0.99,nvwpmbin,wpmbins);
		
		sprintf(name,"H2D_DeepAK8_tbMass_AK4btag%i_AK8CHS_nobmasscut",ib+1);
		hist_2D_DeepAK8_tbmass_2[ib] = new TH2D(name,name,20,-0.01,0.99,nvwpmbin,wpmbins);
		hist_2D_DeepAK8_tbmass_2[ib] = new TH2D(name,name,20,-0.01,0.99,nvwpmbin,wpmbins);
		
		sprintf(name,"H2D_DeepAK8_tbMass_AK4btag%i_AK8CHS",ib+1);
		hist_2D_DeepAK8_tbmass[ib] = new TH2D(name,name,20,-0.01,0.99,nvwpmbin,wpmbins);
		hist_2D_DeepAK8_tbmass[ib] = new TH2D(name,name,20,-0.01,0.99,nvwpmbin,wpmbins);
		
		sprintf(name,"H2D_DeepAK8_tbMass_AK4btag%i_AK8CHS_highht",ib+1);
		hist_2D_DeepAK8_tbmass_3[ib] = new TH2D(name,name,20,-0.01,0.99,nvwpmbin,wpmbins);
		hist_2D_DeepAK8_tbmass_3[ib] = new TH2D(name,name,20,-0.01,0.99,nvwpmbin,wpmbins);
		
		sprintf(name,"H2D_MD_DeepAK8_tbMass_AK4btag%i_AK8CHS_lowbpt",ib+1);
		hist_2D_MDDeepAK8_tbmass_1[ib] = new TH2D(name,name,20,-0.01,0.99,nvwpmbin,wpmbins);
		hist_2D_MDDeepAK8_tbmass_1[ib] = new TH2D(name,name,20,-0.01,0.99,nvwpmbin,wpmbins);
		
		sprintf(name,"H2D_MD_DeepAK8_tbMass_AK4btag%i_AK8CHS_nobmasscut",ib+1);
		hist_2D_MDDeepAK8_tbmass_2[ib] = new TH2D(name,name,20,-0.01,0.99,nvwpmbin,wpmbins);
		hist_2D_MDDeepAK8_tbmass_2[ib] = new TH2D(name,name,20,-0.01,0.99,nvwpmbin,wpmbins);
		
		sprintf(name,"H2D_MD_DeepAK8_tbMass_AK4btag%i_AK8CHS",ib+1);
		hist_2D_MDDeepAK8_tbmass[ib] = new TH2D(name,name,20,-0.01,0.99,nvwpmbin,wpmbins);
		hist_2D_MDDeepAK8_tbmass[ib] = new TH2D(name,name,20,-0.01,0.99,nvwpmbin,wpmbins);
		
		sprintf(name,"H2D_MD_DeepAK8_tbMass_AK4btag%i_AK8CHS_highht",ib+1);
		hist_2D_MDDeepAK8_tbmass_3[ib] = new TH2D(name,name,20,-0.01,0.99,nvwpmbin,wpmbins);
		hist_2D_MDDeepAK8_tbmass_3[ib] = new TH2D(name,name,20,-0.01,0.99,nvwpmbin,wpmbins);
		
			
		for(int isb=0; isb<nsbtag; isb++){
		
		sprintf(name,"H2D_SDMass_tau32_subjetbtag%i_AK4btag%i_AK8CHS",isb+1,ib+1);
		sprintf(title,"2D Correlation of Soft-Drop Mass  and  #tau_{32} | subjet b-tag %i AK4 b-tag %i",isb,ib);
		
		hist_2D_sdmass_tau32_AK8_reg[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
		hist_2D_sdmass_tau32_AK8_reg[ib][isb]->Sumw2();
		
		sprintf(name,"H2D_SDMass_tau32_subjetbtag%i_AK4btag%i_AK8CHS_lowbpt",isb+1,ib+1);
		
		hist_2D_sdmass_tau32_AK8_reg_1[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
		hist_2D_sdmass_tau32_AK8_reg_1[ib][isb]->Sumw2();
			
		sprintf(name,"H2D_SDMass_tau32_subjetbtag%i_AK4btag%i_AK8CHS_nobmasscut",isb+1,ib+1);
		
		hist_2D_sdmass_tau32_AK8_reg_2[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
		hist_2D_sdmass_tau32_AK8_reg_2[ib][isb]->Sumw2();	
			
		sprintf(name,"H2D_SDMass_tau32_subjetbtag%i_AK4btag%i_AK8CHS_highht",isb+1,ib+1);
		
		hist_2D_sdmass_tau32_AK8_reg_3[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
		hist_2D_sdmass_tau32_AK8_reg_3[ib][isb]->Sumw2();	
			
		sprintf(name,"H2D_SDMass_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS",isb+1,ib+1);	
	    sprintf(title,"H2D SDMass DeepAK8 Subjetbtag%i AK4btag%i AK8CHS",isb+1,ib+1);	
	    
	    hist_2D_sdmass_deeptopscore_AK8[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
	    hist_2D_sdmass_deeptopscore_AK8[ib][isb]->Sumw2();
		
	    sprintf(name,"H2D_SDMass_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_lowbpt",isb+1,ib+1);	
	    
	    hist_2D_sdmass_deeptopscore_AK8_1[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
	    hist_2D_sdmass_deeptopscore_AK8_1[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_nobmasscut",isb+1,ib+1);	
	   
	    hist_2D_sdmass_deeptopscore_AK8_2[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
	    hist_2D_sdmass_deeptopscore_AK8_2[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_highht",isb+1,ib+1);	
	   
	    hist_2D_sdmass_deeptopscore_AK8_3[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
	    hist_2D_sdmass_deeptopscore_AK8_3[ib][isb]->Sumw2();
	
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS",isb+1,ib+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 SubJet btag%i AK4btag%i AK8CHS",isb+1,ib+1);	
	    	
	    hist_2D_sdmass_deeptopscore_md_AK8[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
	    hist_2D_sdmass_deeptopscore_md_AK8[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_lowbpt",isb+1,ib+1);
	    	
	    hist_2D_sdmass_deeptopscore_md_AK8_1[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
	    hist_2D_sdmass_deeptopscore_md_AK8_1[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_nobmasscut",isb+1,ib+1);
	    	
	    hist_2D_sdmass_deeptopscore_md_AK8_2[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
	    hist_2D_sdmass_deeptopscore_md_AK8_2[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_highht",isb+1,ib+1);
	    	
	    hist_2D_sdmass_deeptopscore_md_AK8_3[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
	    hist_2D_sdmass_deeptopscore_md_AK8_3[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD",isb+1,ib+1);
	    sprintf(title,"H2D SDMass DeepAK8 Subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR)",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
	    hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_highht",isb+1,ib+1);
	    sprintf(title,"H2D SDMass DeepAK8 Subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR)",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8_3[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
	    hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8_3[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD",isb+1,ib+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 Subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR)",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_highht",isb+1,ib+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 Subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR)",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_3[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_3[ib][isb]->Sumw2();
	    
	    // t-tbar //
	    
	    sprintf(name,"H2D_SDMass_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_ttbar",isb+1,ib+1);
	    sprintf(title,"H2D SDMass DeepAK8 Subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR) t-#bar{t}",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8_tt[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
	    hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8_tt[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_ttbar_highht",isb+1,ib+1);
	    sprintf(title,"H2D SDMass DeepAK8 Subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR) t-#bar{t}",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8_tt_3[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
	    hist_2D_sdmass_deeptopscore_AK8_md_DeepAK8_tt_3[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_ttbar",isb+1,ib+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 Subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR) t-#bar{t}",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_tt[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_tt[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag_MD_ttbar_highht",isb+1,ib+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 Subjetbtag%i AK4btag%i AK8CHS (DeepTag MD SR) t-#bar{t}",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_tt_3[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
	    hist_2D_sdmass_deeptopscore_md_AK8_md_DeepAK8_tt_3[ib][isb]->Sumw2();
	    
	    // t-tbar ends 
	    
	    sprintf(name,"H2D_SDMass_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag",isb+1,ib+1);
	    sprintf(title,"H2D SDMass DeepAK8 Subjetbtag%i AK4btag%i AK8CHS (DeepTag SR)",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_AK8_DeepAK8[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
	    hist_2D_sdmass_deeptopscore_AK8_DeepAK8[ib][isb]->Sumw2();
	    
	    sprintf(name,"H2D_SDMass_MD_DeepAK8_subjetbtag%i_AK4btag%i_AK8CHS_DeepTag",isb+1,ib+1);
	    sprintf(title,"H2D SDMass Mass Decorrelated DeepAK8 Subjetbtag%i AK4btag%i AK8CHS (DeepTag SR)",isb+1,ib+1);	 
	    
	    hist_2D_sdmass_deeptopscore_md_AK8_DeepAK8[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,-0.01,0.99);
	    hist_2D_sdmass_deeptopscore_md_AK8_DeepAK8[ib][isb]->Sumw2();
	   
		
			} //isb
		}//ib
		
	for(int isb=0; isb<nsbtag; isb++){
		for(int it=0; it<ntoptag; it++){
		
		sprintf(name,"H2D_AK4_pt_btag_TopTag%i_Subjetbtag%i_AK8CHS",it+1,isb+1);
		if(it==0) { sprintf(title,"2D Correlation of P_{T}  and b-tag of AK4CHS Jets | tight top subjet b-tag %i",isb+1); }
		if(it==1) { sprintf(title,"2D Correlation of P_{T}  and b-tag of AK4CHS Jets | medium top subjet b-tag %i",isb+1); }
		if(it==2) { sprintf(title,"2D Correlation of P_{T}  and b-tag of AK4CHS Jets | failed top subjet b-tag %i",isb+1); }
		
		hist_2D_pt_btag_AK4[isb][it] = new TH2D(name,title,noptbins,ptbins,20,-0.01,0.99);
		hist_2D_pt_btag_AK4[isb][it]->Sumw2();
		
		
        sprintf(name,"H2D_AK4_pt_btag_TopTag%i_Subjetbtag%i_AK8CHS_lowbpt",it+1,isb+1);
		
		hist_2D_pt_btag_AK4_1[isb][it] = new TH2D(name,title,noptbins,ptbins,20,-0.01,0.99);
		hist_2D_pt_btag_AK4_1[isb][it]->Sumw2();
		
		sprintf(name,"H2D_AK4_pt_btag_TopTag%i_Subjetbtag%i_AK8CHS_nobmasscut",it+1,isb+1);
		
		hist_2D_pt_btag_AK4_2[isb][it] = new TH2D(name,title,noptbins,ptbins,20,-0.01,0.99);
		hist_2D_pt_btag_AK4_2[isb][it]->Sumw2();
		
		sprintf(name,"H2D_AK4_pt_btag_TopTag%i_Subjetbtag%i_AK8CHS_highht",it+1,isb+1);
		
		hist_2D_pt_btag_AK4_3[isb][it] = new TH2D(name,title,noptbins,ptbins,20,-0.01,0.99);
		hist_2D_pt_btag_AK4_3[isb][it]->Sumw2();
		
		sprintf(name,"H2D_AK4_pt_btag_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS",it+1,isb+1);
		hist_2D_pt_btag_AK4_DeepAK8[isb][it] = new TH2D(name,title,noptbins,ptbins,20,-0.01,0.99);
		hist_2D_pt_btag_AK4_DeepAK8[isb][it]->Sumw2();
		
		sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS",it+1,isb+1);
		hist_2D_pt_btag_AK4_md_DeepAK8[isb][it] = new TH2D(name,title,noptbins,ptbins,20,-0.01,0.99);
		hist_2D_pt_btag_AK4_md_DeepAK8[isb][it]->Sumw2();
		
		sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_highht",it+1,isb+1);
		hist_2D_pt_btag_AK4_md_DeepAK8_3[isb][it] = new TH2D(name,title,noptbins,ptbins,20,-0.01,0.99);
		hist_2D_pt_btag_AK4_md_DeepAK8_3[isb][it]->Sumw2();
		
		sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_ttbar",it+1,isb+1);
		hist_2D_pt_btag_AK4_md_DeepAK8_tt[isb][it] = new TH2D(name,title,noptbins,ptbins,20,-0.01,0.99);
		hist_2D_pt_btag_AK4_md_DeepAK8_tt[isb][it]->Sumw2();
		
		sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_ttbar_highht",it+1,isb+1);
		hist_2D_pt_btag_AK4_md_DeepAK8_tt_3[isb][it] = new TH2D(name,title,noptbins,ptbins,20,-0.01,0.99);
		hist_2D_pt_btag_AK4_md_DeepAK8_tt_3[isb][it]->Sumw2();
		
			}
		}
   
   ////
    
   ////
   
   sprintf(name,"DeltaR_AK4_TopTaggedAK8");
   sprintf(title,"#Delta R (AK4 jets, Top Candidate AK8 jet");
   hist_delR_AK4_toptagAK8 = new TH1D(name,title,400,0,4.0);
   hist_delR_AK4_toptagAK8->Sumw2();
   
   sprintf(name,"DeltaEta_AK4_TopTaggedAK8");
   sprintf(title,"#Delta Rapidity (AK4 jets, Top Candidate AK8 jet");
   hist_deleta_AK4_toptagAK8 = new TH1D(name,title,100,-5.,5.);
   hist_deleta_AK4_toptagAK8->Sumw2();
   
   sprintf(name,"DeltaPhi_AK4_TopTaggedAK8");
   sprintf(title,"#Delta #Phi (AK4 jets, Top Candidate AK8 jet");
   hist_delphi_AK4_toptagAK8 = new TH1D(name,title,150,-M_PI,M_PI);
   hist_delphi_AK4_toptagAK8->Sumw2();
   
   sprintf(name,"DeltaPhi_AK4_TopTaggedAK8_Passing_DeltaR");
   sprintf(title,"#Delta #Phi (AK4 jets, Top Candidate AK8 jet passing #DeltaR"); 
   hist_delphi_AK4_toptagAK8_dRpass = new TH1D(name,title,150,-M_PI,M_PI);
   hist_delphi_AK4_toptagAK8_dRpass->Sumw2();
   
   sprintf(name,"DeltaEta_AK4_TopTaggedAK8_Passing_DeltaR");
   sprintf(title,"#Delta Rapidity (AK4 jets, Top Candidate AK8 jet passing #DeltaR"); 
   hist_deleta_AK4_toptagAK8_dRpass = new TH1D(name,title,100,-5,5);
   hist_deleta_AK4_toptagAK8_dRpass->Sumw2();
   
   sprintf(name,"DeltaR_BTaggedAK4_TopTaggedAK8");
   sprintf(title,"#Delta R (B Candidate AK4 jet, Top Candidate AK8 jet");
   hist_delR_btag_toptag = new TH1D(name,title,400,0,4.0);
   hist_delR_btag_toptag->Sumw2();
   
   sprintf(name,"DeltaPhi_BTaggedAK4_TopTaggedAK8");
   sprintf(title,"#Delta #Phi (B Candidate AK4 jet, Top Candidate AK8 jet");
   hist_delphi_btag_toptag = new TH1D(name,title,150,-M_PI,M_PI);
   hist_delphi_btag_toptag->Sumw2();
   
   sprintf(name,"NJets_AK4CHS");
   sprintf(title,"# of AK4CHS Jets");
   hist_njetAK4 = new TH1D(name,title,10,-0.1,9.9);
   hist_njetAK4->Sumw2();
   
   sprintf(name,"NBJets_AK4CHS");
   sprintf(title,"# of B-Tagged AK4CHS Jets");
   hist_nbjetAK4 = new TH1D(name,title,10,-0.1,9.9);
   hist_nbjetAK4->Sumw2();
   
   sprintf(name,"NJets_AK4CHS_Cut");
   sprintf(title,"# of AK4CHS Jets (N-1 cut)");
   hist_njetAK4_cut = new TH1D(name,title,10,-0.1,9.9);
   hist_njetAK4_cut->Sumw2();
   
   sprintf(name,"NBJets_AK4CHS_Cut");
   sprintf(title,"# of B-Tagged AK4CHS Jets (N-1 cut)");
   hist_nbjetAK4_cut = new TH1D(name,title,10,-0.1,9.9);
   hist_nbjetAK4_cut->Sumw2();
   
   sprintf(name,"Pt_AK4CHS");
   sprintf(title,"P_{T} of AK4CHS Jets");
   hist_jetptAK4 = new TH1D(name,title,noptbins,ptbins);
   hist_jetptAK4->Sumw2();
   
   sprintf(name,"Mass_AK4CHS");
   sprintf(title,"Mass of AK4CHS Jets");
   hist_jetmassAK4 = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetmassAK4->Sumw2();
  
   sprintf(name,"y_AK4CHS");
   sprintf(title,"Rapidity of AK4CHS Jets");
   hist_jetrapAK4 = new TH1D(name,title,50,-5.,5.);
   hist_jetrapAK4->Sumw2();
   
   sprintf(name,"Btag_AK4CHS");
   sprintf(title,"CSVv2 Score of AK4CHS Jets");
   hist_jetbtagAK4 = new TH1D(name,title,20,-0.01,0.99);
   hist_jetbtagAK4->Sumw2();
   
   sprintf(name,"Btag_DeepCSV_AK4CHS");
   sprintf(title,"DeepCSV Score of AK4CHS Jets");
   hist_jetbtagdeepCSVAK4 = new TH1D(name,title,20,-0.01,0.99);
   hist_jetbtagdeepCSVAK4->Sumw2();
   
   sprintf(name,"Btag_DeepJetFlavB_AK4CHS");
   sprintf(title,"DeepJet B Flavour Score of AK4CHS Jets");
   hist_jetbtagdeepflavAK4 = new TH1D(name,title,20,-0.01,0.99);
   hist_jetbtagdeepflavAK4->Sumw2();
   
   sprintf(name,"Parton_Flavour_AK4Jets");
   sprintf(title,"Parton Flavour of AK4Jets");
   hist_jetpartonflavAK4 = new TH1D(name,title,22,-0.01,21.99);
   hist_jetpartonflavAK4->Sumw2();
   
   sprintf(name,"Parton_Flavour_BTagged_AK4Jets");
   sprintf(title,"Parton Flavour of BTagged AK4Jets");
   hist_jetpartonflavAK4_btagged = new TH1D(name,title,22,-0.01,21.99);
   hist_jetpartonflavAK4_btagged->Sumw2();
   
   sprintf(name,"Mass_AK4CHS_highpt");
   sprintf(title,"Mass of AK4CHS Jets (P_{T} > 400 GeV)");
   hist_jetmassAK4_1 = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_jetmassAK4_1->Sumw2();
   
   sprintf(name,"Btag_AK4CHS_highpt");
   sprintf(title,"CSVv2 Score of AK4CHS Jets (P_{T} > 400 GeV)");
   hist_jetbtagAK4_1 = new TH1D(name,title,20,-0.01,0.99);
   hist_jetbtagAK4_1->Sumw2();
   
   sprintf(name,"Btag_DeepCSV_AK4CHS_highpt");
   sprintf(title,"DeepCSV Score of AK4CHS Jets (P_{T} > 400 GeV)");
   hist_jetbtagdeepCSVAK4_1 = new TH1D(name,title,20,-0.01,0.99);
   hist_jetbtagdeepCSVAK4_1->Sumw2();
   
   sprintf(name,"Btag_DeepJetFlavB_AK4CHS_highpt");
   sprintf(title,"DeepJet B Flavour Score of AK4CHS Jets (P_{T} > 400 GeV)");
   hist_jetbtagdeepflavAK4_1 = new TH1D(name,title,20,-0.01,0.99);
   hist_jetbtagdeepflavAK4_1->Sumw2();
   
   sprintf(name,"Parton_Flavour_BTagged_AK4Jets_highbpt");
   sprintf(title,"Parton Flavour of BTagged AK4Jets");
   hist_jetpartonflavAK4_btagged_1 = new TH1D(name,title,22,-0.01,21.99);
   hist_jetpartonflavAK4_btagged_1->Sumw2();
   
   sprintf(name,"Parton_Flavour_AK4Jets_highpt");
   sprintf(title,"Parton Flavour of AK4Jets (P_{T} > 400 GeV)");
   hist_jetpartonflavAK4_1 = new TH1D(name,title,22,-0.01,21.99);
   hist_jetpartonflavAK4_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8JetMass");
   sprintf(title,"Soft-Drop Mass of AK8 Jet Associated to b-Jet");
   hist_biso_mass = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_biso_mass->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8AK4JetMassDiff");
   sprintf(title,"Difference of Mass of AK8 Jet Associated to b-Jet and b-Jet");
   hist_biso_tbmassdiff = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_biso_tbmassdiff->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK4byAK8_Mass");
   sprintf(title,"Ratio of Soft-Drop Mass of AK8 Jet Associated to b-Jet to b-Jet Mass");
   hist_biso_isomass = new TH1D(name,title,20,-0.001,0.999);
   hist_biso_isomass->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK4byAK8_Pt");
   sprintf(title,"Ratio of P_{T} of AK8 Jet Associated to b-Jet to b-Jet");
   hist_biso_isopt = new TH1D(name,title,20,-0.001,0.999);
   hist_biso_isopt->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_TopTagScore");
   sprintf(title,"Top-tag Score of AK8 Jet Associated to b-Jet to b-Jet");
   hist_biso_TopAK8score = new TH1D(name,title,100,-0.001,0.999);
   hist_biso_TopAK8score->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_MDTopTagScore");
   sprintf(title,"MD Top-tag Score of AK8 Jet Associated to b-Jet to b-Jet");
   hist_biso_TopAK8score_MD = new TH1D(name,title,100,-0.001,0.999);
   hist_biso_TopAK8score_MD->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_WTagScore");
   sprintf(title,"W-tag Score of AK8 Jet Associated to b-Jet to b-Jet");
   hist_biso_WAK8score = new TH1D(name,title,100,-0.001,0.999);
   hist_biso_WAK8score->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_MDWTagScore");
   sprintf(title,"MD W-tag Score of AK8 Jet Associated to b-Jet to b-Jet");
   hist_biso_WAK8score_MD = new TH1D(name,title,100,-0.001,0.999);
   hist_biso_WAK8score_MD->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_GenJetMass");
   sprintf(title,"Mass of Gen AK8 Jet Associated to b-Jet to b-Jet");
   hist_bjetAK8genmass = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_bjetAK8genmass->Sumw2();
   
 
   // high b pt //
   
   sprintf(name,"IsoCheck_bjet_AK8JetMass_highbpt");
   sprintf(title,"Soft-Drop Mass of AK8 Jet Associated to b-Jet (p_{T}>400 GeV)");
   hist_biso_mass_1 = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_biso_mass_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8JetMass_DeepAK8cut_highbpt");
   sprintf(title,"Soft-Drop Mass of AK8 Jet Associated to b-Jet with DeepAK8 cut (p_{T}>400 GeV)");
   hist_biso_mass_deepak8cut_1 = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_biso_mass_deepak8cut_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8JetMass_Msdcut_highbpt");
   sprintf(title,"Soft-Drop Mass of AK8 Jet Associated to b-Jet with m_{SD} cut (p_{T}>400 GeV)");
   hist_biso_mass_msdcut_1 = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_biso_mass_msdcut_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8AK4JetMassDiff_highbpt");
   sprintf(title,"Difference of Mass of AK8 Jet Associated to b-Jet and b-Jet (p_{T}>400 GeV)");
   hist_biso_tbmassdiff_1 = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_biso_tbmassdiff_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK4byAK8_Mass_highbpt");
   sprintf(title,"Ratio of Soft-Drop Mass of AK8 Jet Associated to b-Jet to b-Jet Mass (p_{T}>400 GeV)");
   hist_biso_isomass_1 = new TH1D(name,title,20,-0.001,0.999);
   hist_biso_isomass_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK4byAK8_Pt_highbpt");
   sprintf(title,"Ratio of P_{T} of AK8 Jet Associated to b-Jet to b-Jet (p_{T}>400 GeV)");
   hist_biso_isopt_1 = new TH1D(name,title,20,-0.001,0.999);
   hist_biso_isopt_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_TopTagScore_highbpt");
   sprintf(title,"Top-tag Score of AK8 Jet Associated to b-Jet to b-Jet (p_{T}>400 GeV)");
   hist_biso_TopAK8score_1 = new TH1D(name,title,20,-0.001,0.999);
   hist_biso_TopAK8score_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_MDTopTagScore_highbpt");
   sprintf(title,"MD Top-tag Score of AK8 Jet Associated to b-Jet to b-Jet (p_{T}>400 GeV)");
   hist_biso_TopAK8score_MD_1 = new TH1D(name,title,20,-0.001,0.999);
   hist_biso_TopAK8score_MD_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_WTagScore_highbpt");
   sprintf(title,"W-tag Score of AK8 Jet Associated to b-Jet to b-Jet (p_{T}>400 GeV)");
   hist_biso_WAK8score_1 = new TH1D(name,title,20,-0.001,0.999);
   hist_biso_WAK8score_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_MDWTagScore_highbpt");
   sprintf(title,"W-tag Score of AK8 Jet Associated to b-Jet to b-Jet (p_{T}>400 GeV)");
   hist_biso_WAK8score_MD_1 = new TH1D(name,title,20,-0.001,0.999);
   hist_biso_WAK8score_MD_1->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_GenJetMass_highpt");
   sprintf(title,"Mass of Gen AK8 Jet Associated to b-Jet to b-Jet (p_{T}>400 GeV)");
   hist_bjetAK8genmass_1 = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_bjetAK8genmass_1->Sumw2();
   
   for(int ib=0; ib<2; ib++){
   
   sprintf(name,"IsoCheck_bjet_AK8JetMass_highbpt_btag%i",ib);
   sprintf(title,"Soft-Drop Mass of AK8 Jet Associated to b-Jet (p_{T}>400 GeV)");
   hist_biso_mass_wb_1[ib] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_biso_mass_wb_1[ib]->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8AK4JetMassDiff_highbpt_btag%i",ib);
   sprintf(title,"Difference of Mass of AK8 Jet Associated to b-Jet and b-Jet (p_{T}>400 GeV)");
   hist_biso_tbmassdiff_wb_1[ib] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_biso_tbmassdiff_wb_1[ib]->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK4byAK8_Mass_highbpt_btag%i",ib);
   sprintf(title,"Ratio of Soft-Drop Mass of AK8 Jet Associated to b-Jet to b-Jet Mass (p_{T}>400 GeV)");
   hist_biso_isomass_wb_1[ib] = new TH1D(name,title,20,-0.001,0.999);
   hist_biso_isomass_wb_1[ib]->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK4byAK8_Pt_highbpt_btag%i",ib);
   sprintf(title,"Ratio of P_{T} of AK8 Jet Associated to b-Jet to b-Jet (p_{T}>400 GeV)");
   hist_biso_isopt_wb_1[ib] = new TH1D(name,title,20,-0.001,0.999);
   hist_biso_isopt_wb_1[ib]->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_TopTagScore_highbpt_btag%i",ib);
   sprintf(title,"Top-tag Score of AK8 Jet Associated to b-Jet to b-Jet (p_{T}>400 GeV)");
   hist_biso_TopAK8score_wb_1[ib] = new TH1D(name,title,20,-0.001,0.999);
   hist_biso_TopAK8score_wb_1[ib]->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_MDTopTagScore_highbpt_btag%i",ib);
   sprintf(title,"MD Top-tag Score of AK8 Jet Associated to b-Jet to b-Jet (p_{T}>400 GeV)");
   hist_biso_TopAK8score_MD_wb_1[ib] = new TH1D(name,title,20,-0.001,0.999);
   hist_biso_TopAK8score_MD_wb_1[ib]->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_WTagScore_highbpt_btag%i",ib);
   sprintf(title,"W-tag Score of AK8 Jet Associated to b-Jet to b-Jet (p_{T}>400 GeV)");
   hist_biso_WAK8score_wb_1[ib] = new TH1D(name,title,20,-0.001,0.999);
   hist_biso_WAK8score_wb_1[ib]->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_MDWTagScore_highbpt_btag%i",ib);
   sprintf(title,"W-tag Score of AK8 Jet Associated to b-Jet to b-Jet (p_{T}>400 GeV)");
   hist_biso_WAK8score_MD_wb_1[ib] = new TH1D(name,title,20,-0.001,0.999);
   hist_biso_WAK8score_MD_wb_1[ib]->Sumw2();
   
   sprintf(name,"IsoCheck_bjet_AK8_GenJetMass_highpt_btag%i",ib);
   sprintf(title,"Mass of Gen AK8 Jet Associated to b-Jet to b-Jet (p_{T}>400 GeV)");
   hist_bjetAK8genmass_wb_1[ib] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_bjetAK8genmass_wb_1[ib]->Sumw2();
   
   sprintf(name,"MycutPass_highpt_btag%i",ib);
   sprintf(title,"Jets passing my cuts on b-jet (p_{T}>400 GeV)");
   hist_bjetpass_wb_1[ib] = new TH1D(name,title,2,-0.1,1.9);
   hist_bjetpass_wb_1[ib]->Sumw2();
   
   sprintf(name,"MycutFail_highpt_btag%i",ib);
   sprintf(title,"Jets failing my cuts on b-jet (p_{T}>400 GeV)");
   hist_bjetfail_wb_1[ib] = new TH1D(name,title,2,-0.1,1.9);
   hist_bjetfail_wb_1[ib]->Sumw2();
  
	}
   
   // t iso hist //
   
   sprintf(name,"IsoCheck_tjet_SDMass");
   sprintf(title,"Soft-Drop Mass of AK8 Jet opposite to #mu");
   hist_tiso_mass = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_tiso_mass->Sumw2();
   
   sprintf(name,"IsoCheck_tjet_SDMass_MDDeepAK8cut");
   sprintf(title,"Soft-Drop Mass of AK8 Jet opposite to #mu");
   hist_tiso_mass_deepak8cut = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_tiso_mass_deepak8cut->Sumw2();
   
   sprintf(name,"IsoCheck_tjet_SDMass_MDDeepAK8cut_FlatToppT");
   sprintf(title,"Soft-Drop Mass of AK8 Jet opposite to #mu");
   hist_tiso_mass_deepak8cut_flattoppt = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_tiso_mass_deepak8cut_flattoppt->Sumw2();
   
   sprintf(name,"IsoCheck_tjet_SDMass_MWcut");
   sprintf(title,"Soft-Drop Mass of AK8 Jet opposite to #mu");
   hist_tiso_mass_mwcut = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_tiso_mass_mwcut->Sumw2();
   
   // t iso end //
   
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
   sprintf(title,"TopJet Soft-Drop Mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
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
   hist_topjetsbtag_sel[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_topjetsbtag_sel[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepCSV_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet SubJet DeepCSV Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdeepCSV_sel[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_topjetsdeepCSV_sel[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Tau32_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet #tau_{32} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjettau32_sel[it][ib][isb] = new TH1D(name,title,50,-0.01,0.99);
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
   hist_bjetbtag_sel[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_bjetbtag_sel[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_TopTagScore_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Top-tag Score of AK8 Jet Corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Topscore_sel[it][ib][isb] = new TH1D(name,title,100,-0.01,0.99);
   hist_bjetmatch_AK8Topscore_sel[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_WTagScore_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"W-tag Score of AK8 Jet Corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Wscore_sel[it][ib][isb] = new TH1D(name,title,100,-0.01,0.99);
   hist_bjetmatch_AK8Wscore_sel[it][ib][isb]->Sumw2();
   
   sprintf(name,"Event_HT_DeepAK8MD_toptag%i_btag%i_subbtag%i_tau32",it,ib,isb);
   sprintf(title,"H_{T} | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_ht_tau32[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_ht_tau32[it][ib][isb]->Sumw2();

   // Low b pt //
   
   sprintf(name,"tb_Mass_toptag%i_btag%i_subbtag%i_lowbpt",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_1[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_1[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Pt_toptag%i_btag%i_subbtag%i_lowbpt",it,ib,isb);
   sprintf(title,"Pt of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbpt_1[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_tbpt_1[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_y_toptag%i_btag%i_subbtag%i_lowbpt",it,ib,isb);
   sprintf(title,"Rapidity of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbrap_1[it][ib][isb] = new TH1D(name,title,50,-5.,5.);
   hist_tbrap_1[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_z_toptag%i_btag%i_subbtag%i_lowbpt",it,ib,isb);
   sprintf(title,"Min(p_{T}^{t},p_{T}^{b}) / (p_{T}^{t}+p_{T}^{b})  | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1); 
   hist_tbpt_frac_1[it][ib][isb] = new TH1D(name,title,50,0,0.5);
   hist_tbpt_frac_1[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_deltheta_toptag%i_btag%i_subbtag%i_lowbpt",it,ib,isb);
   sprintf(title,"cos(#Delta#theta^{tb}) | Top WP %i B WP %i SubB WP %i ",it+1,ib+1,isb+1);
   hist_tbtheta_diff_1[it][ib][isb] = new TH1D(name,title,50,-1.,+1.);
   hist_tbtheta_diff_1[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_delphi_toptag%i_btag%i_subbtag%i_lowbpt",it,ib,isb);
   sprintf(title,"#Delta#phi^{tb} | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbphi_diff_1[it][ib][isb] = new TH1D(name,title,150,-M_PI,M_PI);
   hist_tbphi_diff_1[it][ib][isb]->Sumw2();
   
   sprintf(name,"BCandidate_Jet_Parton_Flavour_toptag%i_btag%i_subbtag%i_lowbpt",it,ib,isb);
   sprintf(title,"BCandidate Jet Parton Flavour | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_partonflav_AK4_1[it][ib][isb] = new TH1D(name,title,22,-0.01,21.99);
   hist_partonflav_AK4_1[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Pt_toptag%i_btag%i_subbtag%i_lowbpt",it,ib,isb);
   sprintf(title,"TopJet P_{T} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetpt_sel_1[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_topjetpt_sel_1[it][ib][isb]->Sumw2();
    
   sprintf(name,"TopJet_SDMass_toptag%i_btag%i_subbtag%i_lowbpt",it,ib,isb);
   sprintf(title,"TopJet Soft-Drop Mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel_1[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel_1[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_DeepTagTvsQCD_toptag%i_btag%i_subbtag%i_lowbpt",it,ib,isb);
   sprintf(title,"TopJet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetdeeptopscore_sel_1[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeeptopscore_sel_1[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_MD_DeepTagTvsQCD_toptag%i_btag%i_subbtag%i_lowbpt",it,ib,isb);
   sprintf(title,"TopJet Mass-Decorrelated DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetdeepmdtopscore_sel_1[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeepmdtopscore_sel_1[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_toptag%i_btag%i_subbtag%i_lowbpt",it,ib,isb);
   sprintf(title,"TopJet SubJet CSVv2 Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsbtag_sel_1[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_topjetsbtag_sel_1[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepCSV_toptag%i_btag%i_subbtag%i_lowbpt",it,ib,isb);
   sprintf(title,"TopJet SubJet DeepCSV Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdeepCSV_sel_1[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_topjetsdeepCSV_sel_1[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Tau32_toptag%i_btag%i_subbtag%i_lowbpt",it,ib,isb);
   sprintf(title,"TopJet #tau_{32} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjettau32_sel_1[it][ib][isb] = new TH1D(name,title,50,-0.01,0.99);
   hist_topjettau32_sel_1[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_Pt_toptag%i_btag%i_subbtag%i_lowbpt",it,ib,isb);
   sprintf(title,"BJet P_{T} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetpt_sel_1[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_bjetpt_sel_1[it][ib][isb]->Sumw2();
    
   sprintf(name,"BJet_Mass_toptag%i_btag%i_subbtag%i_lowbpt",it,ib,isb);
   sprintf(title,"BJet Mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmass_sel_1[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_bjetmass_sel_1[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_BTag_toptag%i_btag%i_subbtag%i_lowbpt",it,ib,isb);
   sprintf(title,"BJet CSVv2 Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetbtag_sel_1[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_bjetbtag_sel_1[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_TopTagScore_toptag%i_btag%i_subbtag%i_lowbpt",it,ib,isb);
   sprintf(title,"Top-tag Score of AK8 Jet Corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Topscore_sel_1[it][ib][isb] = new TH1D(name,title,100,-0.01,0.99);
   hist_bjetmatch_AK8Topscore_sel_1[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_WTagScore_toptag%i_btag%i_subbtag%i_lowbpt",it,ib,isb);
   sprintf(title,"W-tag Score of AK8 Jet Corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Wscore_sel_1[it][ib][isb] = new TH1D(name,title,100,-0.01,0.99);
   hist_bjetmatch_AK8Wscore_sel_1[it][ib][isb]->Sumw2();
   
   // Low b pt end 
   
   // no b mass cut //
   
   sprintf(name,"tb_Mass_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_2[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Pt_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"Pt of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbpt_2[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_tbpt_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_y_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"Rapidity of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbrap_2[it][ib][isb] = new TH1D(name,title,50,-5.,5.);
   hist_tbrap_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_z_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"Min(p_{T}^{t},p_{T}^{b}) / (p_{T}^{t}+p_{T}^{b})  | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1); 
   hist_tbpt_frac_2[it][ib][isb] = new TH1D(name,title,50,0,0.5);
   hist_tbpt_frac_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_deltheta_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"cos(#Delta#theta^{tb}) | Top WP %i B WP %i SubB WP %i (p_{T}^{b}<40 GeV)",it+1,ib+1,isb+1);
   hist_tbtheta_diff_2[it][ib][isb] = new TH1D(name,title,50,-1.,+1.);
   hist_tbtheta_diff_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_delphi_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"#Delta#phi^{tb} | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbphi_diff_2[it][ib][isb] = new TH1D(name,title,150,-M_PI,M_PI);
   hist_tbphi_diff_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"BCandidate_Jet_Parton_Flavour_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"BCandidate Jet Parton Flavour | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_partonflav_AK4_2[it][ib][isb] = new TH1D(name,title,22,-0.01,21.99);
   hist_partonflav_AK4_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Pt_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"TopJet P_{T} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetpt_sel_2[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_topjetpt_sel_2[it][ib][isb]->Sumw2();
    
   sprintf(name,"TopJet_SDMass_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"TopJet Soft-Drop Mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel_2[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_DeepTagTvsQCD_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"TopJet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetdeeptopscore_sel_2[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeeptopscore_sel_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_MD_DeepTagTvsQCD_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"TopJet Mass-Decorrelated DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetdeepmdtopscore_sel_2[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeepmdtopscore_sel_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"TopJet SubJet CSVv2 Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsbtag_sel_2[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_topjetsbtag_sel_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepCSV_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"TopJet SubJet DeepCSV Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdeepCSV_sel_2[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_topjetsdeepCSV_sel_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Tau32_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"TopJet #tau_{32} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjettau32_sel_2[it][ib][isb] = new TH1D(name,title,50,-0.01,0.99);
   hist_topjettau32_sel_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_Pt_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"BJet P_{T} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetpt_sel_2[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_bjetpt_sel_2[it][ib][isb]->Sumw2();
    
   sprintf(name,"BJet_Mass_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"BJet Mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmass_sel_2[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_bjetmass_sel_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_BTag_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"BJet CSVv2 Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetbtag_sel_2[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_bjetbtag_sel_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_TopTagScore_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"Top-tag Score of AK8 Jet Corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Topscore_sel_2[it][ib][isb] = new TH1D(name,title,100,-0.01,0.99);
   hist_bjetmatch_AK8Topscore_sel_2[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_WTagScore_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"W-tag Score of AK8 Jet Corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Wscore_sel_2[it][ib][isb] = new TH1D(name,title,100,-0.01,0.99);
   hist_bjetmatch_AK8Wscore_sel_2[it][ib][isb]->Sumw2();

   
   // no bmass cut end //
  
  // HT cut //
   
   sprintf(name,"tb_Mass_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_3[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Pt_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"Pt of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbpt_3[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_tbpt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_y_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"Rapidity of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbrap_3[it][ib][isb] = new TH1D(name,title,50,-5.,5.);
   hist_tbrap_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_z_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"Min(p_{T}^{t},p_{T}^{b}) / (p_{T}^{t}+p_{T}^{b})  | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1); 
   hist_tbpt_frac_3[it][ib][isb] = new TH1D(name,title,50,0,0.5);
   hist_tbpt_frac_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_deltheta_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"cos(#Delta#theta^{tb}) | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbtheta_diff_3[it][ib][isb] = new TH1D(name,title,50,-1.,+1.);
   hist_tbtheta_diff_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_delphi_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"#Delta#phi^{tb} | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbphi_diff_3[it][ib][isb] = new TH1D(name,title,150,-M_PI,M_PI);
   hist_tbphi_diff_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BCandidate_Jet_Parton_Flavour_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"BCandidate Jet Parton Flavour | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_partonflav_AK4_3[it][ib][isb] = new TH1D(name,title,22,-0.01,21.99);
   hist_partonflav_AK4_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Pt_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"TopJet P_{T} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetpt_sel_3[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_topjetpt_sel_3[it][ib][isb]->Sumw2();
    
   sprintf(name,"TopJet_SDMass_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"TopJet Soft-Drop Mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel_3[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_DeepTagTvsQCD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"TopJet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetdeeptopscore_sel_3[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeeptopscore_sel_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_MD_DeepTagTvsQCD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"TopJet Mass-Decorrelated DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetdeepmdtopscore_sel_3[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeepmdtopscore_sel_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"TopJet SubJet CSVv2 Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsbtag_sel_3[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_topjetsbtag_sel_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepCSV_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"TopJet SubJet DeepCSV Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdeepCSV_sel_3[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_topjetsdeepCSV_sel_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Tau32_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"TopJet #tau_{32} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjettau32_sel_3[it][ib][isb] = new TH1D(name,title,50,-0.01,0.99);
   hist_topjettau32_sel_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_Pt_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"BJet P_{T} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetpt_sel_3[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_bjetpt_sel_3[it][ib][isb]->Sumw2();
    
   sprintf(name,"BJet_Mass_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"BJet Mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmass_sel_3[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_bjetmass_sel_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_BTag_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"BJet CSVv2 Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetbtag_sel_3[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_bjetbtag_sel_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_TopTagScore_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"Top-tag Score of AK8 Jet Corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Topscore_sel_3[it][ib][isb] = new TH1D(name,title,100,-0.01,0.99);
   hist_bjetmatch_AK8Topscore_sel_3[it][ib][isb]->Sumw2();

   sprintf(name,"BJet_MatchAK8_WTagScore_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"W-tag Score of AK8 Jet Corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Wscore_sel_3[it][ib][isb] = new TH1D(name,title,100,-0.01,0.99);
   hist_bjetmatch_AK8Wscore_sel_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"Event_HT_DeepAK8MD_toptag%i_btag%i_subbtag%i_tau32_highht",it,ib,isb);
   sprintf(title,"H_{T} | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_ht_tau32_3[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_ht_tau32_3[it][ib][isb]->Sumw2();
   
   
   // HT cut end 
   
   sprintf(name,"TopJet_SDMass_Tau32_tbMass_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet Soft-Drop Mass #tau{32} tb Inv Mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_3D_topjetsdmass_topjettau32_tbmass[it][ib][isb] = new TH3D(name,title,nomassbins,mass_low,mass_high,50,-0.01,0.99,nwpmbin,nwpmlow,nwpmhigh);
   hist_3D_topjetsdmass_topjettau32_tbmass[it][ib][isb]->Sumw2();
   
   // AK8 b Candidate//
   
   sprintf(name,"tb_Mass_AK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_AK8[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_AK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJetAK8_DeepTagTvsQCD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"AK8 b-Jet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_AK8bCand_DeeptagScore[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_AK8bCand_DeeptagScore[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   
   sprintf(name,"BJetAK8_DeepTagTvsQCD_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"AK8 b-Jet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_AK8bCand_DeeptagScore_1[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_AK8bCand_DeeptagScore_1[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   
   sprintf(name,"BJetAK8_MD_DeepTagTvsQCD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"AK8 b-Jet Mass Decorrelated DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_AK8bCand_MD_DeeptagScore[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_AK8bCand_MD_DeeptagScore[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   
   sprintf(name,"tb_Mass_AK8_DeepAK8TopTag_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_AK8_DeepAK8[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_AK8_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJetAK8_DeepTagTvsQCD_DeepAK8TopTag_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"AK8 b-Jet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_AK8bCand_DeeptagScore_DeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_AK8bCand_DeeptagScore_DeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   
   sprintf(name,"BJetAK8_DeepTagTvsQCD_DeepAK8TopTag_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"AK8 b-Jet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_AK8bCand_DeeptagScore_DeepAK8_1[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_AK8bCand_DeeptagScore_DeepAK8_1[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   
   sprintf(name,"BJetAK8_MD_DeepTagTvsQCD_DeepAK8TopTag_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"AK8 b-Jet Mass Decorrelated DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_AK8bCand_MD_DeeptagScore_DeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_AK8bCand_MD_DeeptagScore_DeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   
   sprintf(name,"tb_Mass_AK8_MD_DeepAK8TopTag_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_tbmass_AK8_MDDeepAK8[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_AK8_MDDeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJetAK8_DeepTagTvsQCD_MD_DeepAK8TopTag_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"AK8 b-Jet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_AK8bCand_DeeptagScore_MDDeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_AK8bCand_DeeptagScore_MDDeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   
   sprintf(name,"BJetAK8_DeepTagTvsQCD_MD_DeepAK8TopTag_toptag%i_btag%i_subbtag%i_nobmasscut",it,ib,isb);
   sprintf(title,"AK8 b-Jet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_AK8bCand_DeeptagScore_MDDeepAK8_1[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_AK8bCand_DeeptagScore_MDDeepAK8_1[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   
   sprintf(name,"BJetAK8_MD_DeepTagTvsQCD_MD_DeepAK8TopTag_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"AK8 b-Jet Mass Decorrelated DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_AK8bCand_MD_DeeptagScore_MDDeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_AK8bCand_MD_DeeptagScore_MDDeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   
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
    
   sprintf(name,"TopJet_SDMass_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet Soft-Drop Mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
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
   hist_topjetsbtag_sel_DeepAK8[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_topjetsbtag_sel_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepCSV_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet SubJet DeepCSV Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdeepCSV_sel_DeepAK8[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_topjetsdeepCSV_sel_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Tau32_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet #tau_{32} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjettau32_sel_DeepAK8[it][ib][isb] = new TH1D(name,title,50,-0.01,0.99);
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
   hist_bjetbtag_sel_DeepAK8[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_bjetbtag_sel_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_TopTagScore_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Top-tag Score of AK8 Jet Corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Topscore_sel_DeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.01,0.99);
   hist_bjetmatch_AK8Topscore_sel_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_WTagScore_DeepAK8_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"W-tag Score of AK8 Jet Corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Wscore_sel_DeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.01,0.99);
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
    
   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet Soft-Drop Mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   for(int ipdf=0; ipdf<npdfmax; ipdf++){
	   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_PDF%i",it,ib,isb,ipdf);
	   sprintf(title,"TopJet Soft-Drop Mass Top WP %i B WP %i SubB WP %i PDF %i",it+1,ib+1,isb+1,ipdf+1);
	   hist_topjetsdmass_sel_md_DeepAK8_PDF[it][ib][isb][ipdf] = new TH1D(name,title,nvwpmbin,wpmbins);
	   hist_topjetsdmass_sel_md_DeepAK8_PDF[it][ib][isb][ipdf]->Sumw2();
	}
	   
   for(int iscale=0; iscale<nscalemax; iscale++){
	   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_Scale%i",it,ib,isb,iscale);
	   sprintf(title,"TopJet Soft-Drop Mass Top WP %i B WP %i SubB WP %i Scale %i",it+1,ib+1,isb+1,iscale+1);
	   hist_topjetsdmass_sel_md_DeepAK8_Scale[it][ib][isb][iscale] = new TH1D(name,title,nvwpmbin,wpmbins);
	   hist_topjetsdmass_sel_md_DeepAK8_Scale[it][ib][isb][iscale]->Sumw2();
   }
  
   for(int ips=0; ips<npsmax; ips++){
	   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_PS%i",it,ib,isb,ips);
	   sprintf(title,"TopJet Soft-Drop Mass Top WP %i B WP %i SubB WP %i PS %i",it+1,ib+1,isb+1,ips+1);
	   hist_topjetsdmass_sel_md_DeepAK8_PS[it][ib][isb][ips] = new TH1D(name,title,nvwpmbin,wpmbins);
	   hist_topjetsdmass_sel_md_DeepAK8_PS[it][ib][isb][ips]->Sumw2();
   }
   
   
   sprintf(name,"TopJet_DeepTagTvsQCD_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetdeeptopscore_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeeptopscore_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_MD_DeepTagTvsQCD_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet Mass-Decorrelated DeepTagTvsQCD Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetdeepmdtopscore_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.0001,0.9999);
   hist_topjetdeepmdtopscore_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet SubJet CSVv2 Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsbtag_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_topjetsbtag_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepCSV_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet SubJet DeepCSV Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdeepCSV_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_topjetsdeepCSV_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Tau32_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"TopJet #tau_{32} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjettau32_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,50,-0.01,0.99);
   hist_topjettau32_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"BJet P_{T} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetpt_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_bjetpt_sel_md_DeepAK8[it][ib][isb]->Sumw2();
    
   sprintf(name,"BJet_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"BJet Mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmass_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_bjetmass_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_BTag_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"BJet CSVv2 Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetbtag_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_bjetbtag_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_TopTagScore_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"Top-tag Score of AK8 Jet Corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Topscore_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.01,0.99);
   hist_bjetmatch_AK8Topscore_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_WTagScore_DeepAK8MD_toptag%i_btag%i_subbtag%i",it,ib,isb);
   sprintf(title,"W-tag Score of AK8 Jet Corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Wscore_sel_md_DeepAK8[it][ib][isb] = new TH1D(name,title,100,-0.01,0.99);
   hist_bjetmatch_AK8Wscore_sel_md_DeepAK8[it][ib][isb]->Sumw2();
   
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
   sprintf(title,"TopJet Soft-Drop Mass Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
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
   hist_topjetsbtag_sel_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_topjetsbtag_sel_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepCSV_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"TopJet SubJet DeepCSV Score Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjetsdeepCSV_sel_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_topjetsdeepCSV_sel_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Tau32_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"TopJet #tau_{32} Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_topjettau32_sel_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,50,-0.01,0.99);
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
   hist_bjetbtag_sel_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_bjetbtag_sel_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_TopTagScore_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"Top-tag Score of AK8 Jet Corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Topscore_sel_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,100,-0.01,0.99);
   hist_bjetmatch_AK8Topscore_sel_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_WTagScore_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"W-tag Score of AK8 Jet Corresponding to Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Wscore_sel_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,100,-0.01,0.99);
   hist_bjetmatch_AK8Wscore_sel_md_DeepAK8_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"Event_HT_DeepAK8MD_toptag%i_btag%i_subbtag%i_highht",it,ib,isb);
   sprintf(title,"H_{T} | Top WP %i B WP %i SubB WP %i",it+1,ib+1,isb+1);
   hist_ht_md_DeepAK8_3[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_ht_md_DeepAK8_3[it][ib][isb]->Sumw2();
 
   
   // MD for t-tbar //
   
   sprintf(name,"tb_Mass_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"Mass of tb Pair | Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_tbmass_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_tbmass_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"Pt of tb Pair | Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_tbpt_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_tbpt_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"tb_y_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"Rapidity of tb Pair | Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_tbrap_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,50,-5.,5.);
   hist_tbrap_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"BCandidate_Jet_Parton_Flavour_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"BCandidate Jet Parton Flavour | Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_partonflav_AK4_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,22,-0.01,21.99);
   hist_partonflav_AK4_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Pt_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"TopJet P_{T} Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjetpt_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,noptbins,ptbins);
   hist_topjetpt_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
    
   sprintf(name,"TopJet_SDMass_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"TopJet Soft-Drop Mass Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjetsdmass_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,nsdmassbins,sdmassbins);
   hist_topjetsdmass_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
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
   hist_topjetsbtag_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_topjetsbtag_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepCSV_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"TopJet SubJet DeepCSV Score Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjetsdeepCSV_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_topjetsdeepCSV_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Tau32_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"TopJet #tau_{32} Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjettau32_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,50,-0.01,0.99);
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
   hist_bjetbtag_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_bjetbtag_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_TopTagScore_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"Top-tag Score of AK8 Jet Corresponding to Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Topscore_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,100,-0.01,0.99);
   hist_bjetmatch_AK8Topscore_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_WTagScore_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar",it,ib,isb);
   sprintf(title,"W-tag Score of AK8 Jet Corresponding to Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Wscore_sel_md_DeepAK8_tt[it][ib][isb] = new TH1D(name,title,100,-0.01,0.99);
   hist_bjetmatch_AK8Wscore_sel_md_DeepAK8_tt[it][ib][isb]->Sumw2();
   
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
   sprintf(title,"TopJet Soft-Drop Mass Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
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
   hist_topjetsbtag_sel_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_topjetsbtag_sel_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_SubjetBTag_DeepCSV_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"TopJet SubJet DeepCSV Score Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjetsdeepCSV_sel_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_topjetsdeepCSV_sel_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"TopJet_Tau32_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"TopJet #tau_{32} Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_topjettau32_sel_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,50,-0.01,0.99);
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
   hist_bjetbtag_sel_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,20,-0.01,0.99);
   hist_bjetbtag_sel_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_TopTagScore_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"Top-tag Score of AK8 Jet Corresponding to Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Topscore_sel_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,100,-0.01,0.99);
   hist_bjetmatch_AK8Topscore_sel_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"BJet_MatchAK8_WTagScore_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"W-tag Score of AK8 Jet Corresponding to Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_bjetmatch_AK8Wscore_sel_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,100,-0.01,0.99);
   hist_bjetmatch_AK8Wscore_sel_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   sprintf(name,"Event_HT_DeepAK8MD_toptag%i_btag%i_subbtag%i_ttbar_highht",it,ib,isb);
   sprintf(title,"H_{T} | Top WP %i B WP %i SubB WP %i (t-#bar{t})",it+1,ib+1,isb+1);
   hist_ht_md_DeepAK8_tt_3[it][ib][isb] = new TH1D(name,title,nvwpmbin,wpmbins);
   hist_ht_md_DeepAK8_tt_3[it][ib][isb]->Sumw2();
   
   // MD for t-tbar ends
   
   }
  }
 }

 // trigger histograms //
 
   sprintf(name,"HT_SingleMu_Passed");
   sprintf(title,"HT | SingleMuon Trigger Passed");
   hist_ht_mutrig= new TH1D(name,title,nohtbins,htbins);
   hist_ht_mutrig->Sumw2();
   
   sprintf(name,"LeadJetPt_SingleMu_Passed");
   sprintf(title,"Leading Jet P_{T} | SingleMuon Trigger Passed");
   hist_pt_mutrig= new TH1D(name,title,nohtbins,htbins);
   hist_pt_mutrig->Sumw2();
   
   sprintf(name,"HT_SingleMuElec_Passed");
   sprintf(title,"HT | (SingleMuon||SingleElectron) Trigger Passed");
   hist_ht_emutrig= new TH1D(name,title,nohtbins,htbins);
   hist_ht_emutrig->Sumw2();
   
   sprintf(name,"LeadJetPt_SingleMuElec_Passed");
   sprintf(title,"Leading Jet P_{T} | (SingleMuon||SingleElectron) Trigger Passed");
   hist_pt_emutrig= new TH1D(name,title,nohtbins,htbins);
   hist_pt_emutrig->Sumw2();
   
   sprintf(name,"LeadBJetPt_SingleMuElec_Passed");
   sprintf(title,"B jet p_{T} | (SingleMuon||SingleElectron) trigger passed");
   hist_bpt_emutrig= new TH1D(name,title,nohtbins,htbins);
   hist_bpt_emutrig->Sumw2();
   
   sprintf(name,"TwoJetPtSum_SingleMuElec_Passed");
   sprintf(title,"(P_{T}^{0} + P_{T}^{1})| (SingleMuon||SingleElectron) Trigger Passed");
   hist_ptsum_emutrig= new TH1D(name,title,nohtbins,htbins);
   hist_ptsum_emutrig->Sumw2();
   
   sprintf(name,"Mtb_SingleMuElec_Passed");
   sprintf(title,"M_{tb}| (SingleMuon||SingleElectron) Trigger Passed");
   hist_mtb_emutrig= new TH1D(name,title,nvwpmbin_trig,wpmbins_trig);
   hist_mtb_emutrig->Sumw2();
 
 /////////////// Ref : Electron + Muon Trigger //////////
    
   // ht //
   
   sprintf(name,"HT_SingleMuorElec_HT1050_Passed");
   sprintf(title,"HT |  (SingleMuon || Single Electron) && PFHT1050 trigger passed");
   hist_ht_HT1050_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_ht_HT1050_wemutrig->Sumw2();
   
   sprintf(name,"HT_SingleMuorElec_AK4Pt500_Passed");
   sprintf(title,"HT |  (SingleMuon || Single Electron) && AK4PFJet500 trigger passed");
   hist_ht_AK4Pt500_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_ht_AK4Pt500_wemutrig->Sumw2();
   
   sprintf(name,"HT_SingleMuorElec_AK8Pt500_Passed");
   sprintf(title,"HT |  (SingleMuon || Single Electron) && AK8Pt500 trigger passed");
   hist_ht_AK8Pt500_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_ht_AK8Pt500_wemutrig->Sumw2();
   
   sprintf(name,"HT_SingleMuorElec_AK8PFJet420_TrimMass30_Passed");
   sprintf(title,"HT |  (SingleMuon || Single Electron) && AK8PFJet420_TrimMass30 trigger passed");
   hist_ht_AK8PFJet420_TrimMass30_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_ht_AK8PFJet420_TrimMass30_wemutrig->Sumw2();
   
   sprintf(name,"HT_SingleMuorElec_AK8HT900_TrimMass50_Passed");
   sprintf(title,"HT |  (SingleMuon || Single Electron) && AK8PFHT900_TrimMass50 trigger passed");
   hist_ht_AK8PFHT900_TrimMass50_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_ht_AK8PFHT900_TrimMass50_wemutrig->Sumw2();
   
   sprintf(name,"HT_SingleMuorElec_AK8DiPFJet300_200_TrimMass30_Passed");
   sprintf(title,"HT |  (SingleMuon || Single Electron) && AK8DiPFJet300_200_TrimMass30 trigger passed");
   hist_ht_AK8DiPFJet300_200_TrimMass30_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_ht_AK8DiPFJet300_200_TrimMass30_wemutrig->Sumw2();
   
   sprintf(name,"HT_SingleMuorElec_AnyTrig_Passed");
   sprintf(title,"HT |  (SingleMuon || Single Electron) && any trigger passed");
   hist_ht_all_wemutrig = new TH1D(name,title,nohtbins,htbins);
   hist_ht_all_wemutrig->Sumw2(); 

   
   // pt //
   
   sprintf(name,"LeadJetPt_SingleMuorElec_HT1050_Passed");
   sprintf(title,"Leading jet p_{T} |  (SingleMuon || Single Electron) && PFHT1050 trigger passed");
   hist_pt_HT1050_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_pt_HT1050_wemutrig->Sumw2();
   
   sprintf(name,"LeadJetPt_SingleMuorElec_AK4Pt500_Passed");
   sprintf(title,"Leading jet p_{T} |  (SingleMuon || Single Electron) && AK4PFJet500 trigger passed");
   hist_pt_AK4Pt500_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_pt_AK4Pt500_wemutrig->Sumw2();
   
   sprintf(name,"LeadJetPt_SingleMuorElec_AK8Pt500_Passed");
   sprintf(title,"Leading jet p_{T} |  (SingleMuon || Single Electron) && AK8Pt500 trigger passed");
   hist_pt_AK8Pt500_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_pt_AK8Pt500_wemutrig->Sumw2();
   
   sprintf(name,"LeadJetPt_SingleMuorElec_AK8PFJet420_TrimMass30_Passed");
   sprintf(title,"Leading jet p_{T} |  (SingleMuon || Single Electron) && AK8PFJet420_TrimMass30 trigger passed");
   hist_pt_AK8PFJet420_TrimMass30_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_pt_AK8PFJet420_TrimMass30_wemutrig->Sumw2();
   
   sprintf(name,"LeadJetPt_SingleMuorElec_AK8HT900_TrimMass50_Passed");
   sprintf(title,"Leading jet p_{T} |  (SingleMuon || Single Electron) && AK8PFHT900_TrimMass50 trigger passed");
   hist_pt_AK8PFHT900_TrimMass50_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_pt_AK8PFHT900_TrimMass50_wemutrig->Sumw2();
   
   sprintf(name,"LeadJetPt_SingleMuorElec_AK8DiPFJet300_200_TrimMass30_Passed");
   sprintf(title,"Leading jet p_{T} |  (SingleMuon || Single Electron) && AK8DiPFJet300_200_TrimMass30 trigger passed");
   hist_pt_AK8DiPFJet300_200_TrimMass30_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_pt_AK8DiPFJet300_200_TrimMass30_wemutrig->Sumw2();
   
   sprintf(name,"LeadJetPt_SingleMuorElec_AnyTrig_Passed");
   sprintf(title,"Leading jet p_{T} |  (SingleMuon || Single Electron) && any trigger passed");
   hist_pt_all_wemutrig = new TH1D(name,title,nohtbins,htbins);
   hist_pt_all_wemutrig->Sumw2(); 
   
   // b pt //
   
   sprintf(name,"LeadBJetPt_SingleMuorElec_HT1050_Passed");
   sprintf(title,"B jet p_{T} |  (SingleMuon || Single Electron) && PFHT1050 trigger passed");
   hist_bpt_HT1050_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_bpt_HT1050_wemutrig->Sumw2();
   
   sprintf(name,"LeadBJetPt_SingleMuorElec_AK4Pt500_Passed");
   sprintf(title,"B jet p_{T} |  (SingleMuon || Single Electron) && AK4PFJet500 trigger passed");
   hist_bpt_AK4Pt500_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_bpt_AK4Pt500_wemutrig->Sumw2();
   
   sprintf(name,"LeadBJetPt_SingleMuorElec_AK8Pt500_Passed");
   sprintf(title,"B jet p_{T} |  (SingleMuon || Single Electron) && AK8Pt500 trigger passed");
   hist_bpt_AK8Pt500_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_bpt_AK8Pt500_wemutrig->Sumw2();
   
   sprintf(name,"LeadBJetPt_SingleMuorElec_AK8PFJet420_TrimMass30_Passed");
   sprintf(title,"B jet p_{T} |  (SingleMuon || Single Electron) && AK8PFJet420_TrimMass30 trigger passed");
   hist_bpt_AK8PFJet420_TrimMass30_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_bpt_AK8PFJet420_TrimMass30_wemutrig->Sumw2();
   
   sprintf(name,"LeadBJetPt_SingleMuorElec_AK8HT900_TrimMass50_Passed");
   sprintf(title,"B jet p_{T} |  (SingleMuon || Single Electron) && AK8PFHT900_TrimMass50 trigger passed");
   hist_bpt_AK8PFHT900_TrimMass50_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_bpt_AK8PFHT900_TrimMass50_wemutrig->Sumw2();
   
   sprintf(name,"LeadBJetPt_SingleMuorElec_AK8DiPFJet300_200_TrimMass30_Passed");
   sprintf(title,"B jet p_{T} |  (SingleMuon || Single Electron) && AK8DiPFJet300_200_TrimMass30 trigger passed");
   hist_bpt_AK8DiPFJet300_200_TrimMass30_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_bpt_AK8DiPFJet300_200_TrimMass30_wemutrig->Sumw2();
   
   sprintf(name,"LeadBJetPt_SingleMuorElec_AnyTrig_Passed");
   sprintf(title,"B jet p_{T} |  (SingleMuon || Single Electron) && any trigger passed");
   hist_bpt_all_wemutrig = new TH1D(name,title,nohtbins,htbins);
   hist_bpt_all_wemutrig->Sumw2(); 
   
   
   // ptsum //
   
   sprintf(name,"TwoJetPtSum_SingleMuorElec_HT1050_Passed");
   sprintf(title,"(P_{T}^{0} + P_{T}^{1}) |  (SingleMuon || Single Electron) && PFHT1050 trigger passed");
   hist_ptsum_HT1050_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_ptsum_HT1050_wemutrig->Sumw2();
   
   sprintf(name,"TwoJetPtSum_SingleMuorElec_AK4Pt500_Passed");
   sprintf(title,"(P_{T}^{0} + P_{T}^{1}) |  (SingleMuon || Single Electron) && AK4PFJet500 trigger passed");
   hist_ptsum_AK4Pt500_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_ptsum_AK4Pt500_wemutrig->Sumw2();
   
   sprintf(name,"TwoJetPtSum_SingleMuorElec_AK8Pt500_Passed");
   sprintf(title,"(P_{T}^{0} + P_{T}^{1}) |  (SingleMuon || Single Electron) && AK8Pt500 trigger passed");
   hist_ptsum_AK8Pt500_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_ptsum_AK8Pt500_wemutrig->Sumw2();
   
   sprintf(name,"TwoJetPtSum_SingleMuorElec_AK8PFJet420_TrimMass30_Passed");
   sprintf(title,"(P_{T}^{0} + P_{T}^{1}) |  (SingleMuon || Single Electron) && AK8PFJet420_TrimMass30 trigger passed");
   hist_ptsum_AK8PFJet420_TrimMass30_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_ptsum_AK8PFJet420_TrimMass30_wemutrig->Sumw2();
   
   sprintf(name,"TwoJetPtSum_SingleMuorElec_AK8HT900_TrimMass50_Passed");
   sprintf(title,"(P_{T}^{0} + P_{T}^{1}) |  (SingleMuon || Single Electron) && AK8PFHT900_TrimMass50 trigger passed");
   hist_ptsum_AK8PFHT900_TrimMass50_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_ptsum_AK8PFHT900_TrimMass50_wemutrig->Sumw2();
   
   sprintf(name,"TwoJetPtSum_SingleMuorElec_AK8DiPFJet300_200_TrimMass30_Passed");
   sprintf(title,"(P_{T}^{0} + P_{T}^{1}) |  (SingleMuon || Single Electron) && AK8DiPFJet300_200_TrimMass30 trigger passed");
   hist_ptsum_AK8DiPFJet300_200_TrimMass30_wemutrig= new TH1D(name,title,nohtbins,htbins);
   hist_ptsum_AK8DiPFJet300_200_TrimMass30_wemutrig->Sumw2();
   
   sprintf(name,"TwoJetPtSum_SingleMuorElec_AnyTrig_Passed");
   sprintf(title,"(P_{T}^{0} + P_{T}^{1}) |  (SingleMuon || Single Electron) && any trigger passed");
   hist_ptsum_all_wemutrig = new TH1D(name,title,nohtbins,htbins);
   hist_ptsum_all_wemutrig->Sumw2();
   
   // mtb //
   
   sprintf(name,"Mtb_SingleMuorElec_HT1050_Passed");
   sprintf(title,"M_{tb} |  (SingleMuon || Single Electron) && PFHT1050 trigger passed");
   hist_mtb_HT1050_wemutrig= new TH1D(name,title,nvwpmbin_trig,wpmbins_trig);
   hist_mtb_HT1050_wemutrig->Sumw2();
   
   sprintf(name,"Mtb_SingleMuorElec_AK4Pt500_Passed");
   sprintf(title,"M_{tb} |  (SingleMuon || Single Electron) && AK4PFJet500 trigger passed");
   hist_mtb_AK4Pt500_wemutrig= new TH1D(name,title,nvwpmbin_trig,wpmbins_trig);
   hist_mtb_AK4Pt500_wemutrig->Sumw2();
   
   sprintf(name,"Mtb_SingleMuorElec_AK8Pt500_Passed");
   sprintf(title,"M_{tb} |  (SingleMuon || Single Electron) && AK8Pt500 trigger passed");
   hist_mtb_AK8Pt500_wemutrig= new TH1D(name,title,nvwpmbin_trig,wpmbins_trig);
   hist_mtb_AK8Pt500_wemutrig->Sumw2();
   
   sprintf(name,"Mtb_SingleMuorElec_AK8PFJet420_TrimMass30_Passed");
   sprintf(title,"M_{tb} |  (SingleMuon || Single Electron) && AK8PFJet420_TrimMass30 trigger passed");
   hist_mtb_AK8PFJet420_TrimMass30_wemutrig= new TH1D(name,title,nvwpmbin_trig,wpmbins_trig);
   hist_mtb_AK8PFJet420_TrimMass30_wemutrig->Sumw2();
   
   sprintf(name,"Mtb_SingleMuorElec_AK8HT900_TrimMass50_Passed");
   sprintf(title,"M_{tb} |  (SingleMuon || Single Electron) && AK8PFHT900_TrimMass50 trigger passed");
   hist_mtb_AK8PFHT900_TrimMass50_wemutrig= new TH1D(name,title,nvwpmbin_trig,wpmbins_trig);
   hist_mtb_AK8PFHT900_TrimMass50_wemutrig->Sumw2();
   
   sprintf(name,"Mtb_SingleMuorElec_AK8DiPFJet300_200_TrimMass30_Passed");
   sprintf(title,"M_{tb} |  (SingleMuon || Single Electron) && AK8DiPFJet300_200_TrimMass30 trigger passed");
   hist_mtb_AK8DiPFJet300_200_TrimMass30_wemutrig= new TH1D(name,title,nvwpmbin_trig,wpmbins_trig);
   hist_mtb_AK8DiPFJet300_200_TrimMass30_wemutrig->Sumw2();
   
   sprintf(name,"Mtb_SingleMuorElec_AnyTrig_Passed");
   sprintf(title,"M_{tb} |  (SingleMuon || Single Electron) && any trigger passed");
   hist_mtb_all_wemutrig = new TH1D(name,title,nvwpmbin_trig,wpmbins_trig);
   hist_mtb_all_wemutrig->Sumw2();  
 
 // trigger histo end //


// Fit function

if(!isMC){

bpfrat = new TF1("BTagger_PFRate_SR",MikkoFunc1,500,2000,4);
bpfrat->SetParameters(-1.85166e+01,-4.56084e-01,5.20671e+01,4.91319e-01);

bpfrat_val = new TF1("BTagger_PFRate_VR",MikkoFunc1,500,2000,4);
bpfrat_val->SetParameters(-6.12117e+00,-4.51636e-01,1.71053e+01,2.05272e-01);

tagpfrat_tau32 = new TF1("TopTagger_Tau32_PFRate",pol2,20,400,3);
tagpfrat_tau32->SetParameters(7.65157e-03,-2.05528e-04,3.84749e-06);

tagpfrat_tau32_val = new TF1("TopTagger_Tau32_PFRate_VR",pol2,20,400,3);
tagpfrat_tau32_val->SetParameters(1.70281e-02,-4.12120e-04,4.66585e-06);

tagpfrat_MDAK8 = new TF1("TopTagger_MDDeepAK8_PFRate",pol2,20,400,3);
tagpfrat_MDAK8->SetParameters(6.50814e-02,-3.46971e-04,1.99872e-06);

tagpfrat_MDAK8_val = new TF1("TopTagger_MDDeepAK8_PFRate_VR",pol2,20,400,3);
tagpfrat_MDAK8_val->SetParameters(4.46776e-02,-3.10293e-04,1.97450e-06);

tagpfrat_MDAK8_val2 = new TF1("TopTagger_MDDeepAK8_PFRate_VR",pol2,20,400,3);
tagpfrat_MDAK8_val2->SetParameters(1.22432e-02,-9.98876e-05,3.91304e-07);
}

if(isMC && isQCD){

bpfrat = new TF1("BTagger_PFRate_SR",MikkoFunc1,400,2000,4);
bpfrat->SetParameters(-1.13445e+00,-5.30009e-01,9.42228e-01,6.32387e-02);

bpfrat_val = new TF1("BTagger_PFRate_VR",MikkoFunc1,400,2000,4);
bpfrat_val->SetParameters(1.25644e+01,-6.72034e-01,-1.34758e+01,1.22587e-02);

tagpfrat_tau32 = new TF1("TopTagger_Tau32_PFRate",pol2,20,400,3);
tagpfrat_tau32->SetParameters(2.16355e-03,-8.54422e-05,4.06720e-06);

tagpfrat_tau32_val = new TF1("TopTagger_Tau32_PFRate_VR",pol2,20,400,3);
tagpfrat_tau32_val->SetParameters(1.03655e-02,-2.63081e-04,4.70500e-06);

tagpfrat_MDAK8 = new TF1("TopTagger_MDDeepAK8_PFRate",pol2,20,400,3);
tagpfrat_MDAK8->SetParameters(4.89420e-02,-3.68049e-04,2.27951e-06);

tagpfrat_MDAK8_val= new TF1("TopTagger_MDDeepAK8_PFRate_VR",pol2,20,400,3);
tagpfrat_MDAK8_val->SetParameters(3.55999e-02,-3.00214e-04,2.11683e-06);

tagpfrat_MDAK8_val2= new TF1("TopTagger_MDDeepAK8_PFRate_VR2",pol2,20,400,3);
tagpfrat_MDAK8_val2->SetParameters(1.04260e-02,-1.07329e-04,3.99638e-07);

}


// value assignement //

#ifdef Anal_2016
   
   btagvalue = 0.8484; //0.9535; //tight 0.8484; // medium //csvv2 
   btagvalue_deepCSV = 0.6321; //0.8953  ;//tight  0.6321; //medium // deepcsv 
   btagvalue_deepFlavB = 0.45; //0.7221; // 0.3093;  // medium	//deepflav
   
   deepak8_cut = 0.937;    // 2016 0.5% mistag rate
   deepak8_cut_md = 0.621; // 2016 0.5% mistag rate
   
   for(int ipt=0; ipt<noAK8ptbins; ipt++){
		DeepAK8_SF[ipt] = DeepAK8_SF_16[ipt];
		DeepAK8_MD_SF[ipt] = DeepAK8_MD_SF_16[ipt];
		DeepAK8_MD_SF_up[ipt] = DeepAK8_MD_SF_16_up[ipt];
		DeepAK8_MD_SF_dn[ipt] = DeepAK8_MD_SF_16_dn[ipt];
	}

#endif
	
#ifdef Anal_2017
   
   btagvalue = 0.8838;//0.9693;//tight 0.8838; //csvv2 medium
   btagvalue_deepCSV =  0.4941;//0.65;//myvalue 0.8001;// tight 0.4941; // deepcsv medium
   btagvalue_deepFlavB = 0.6; //0.6; // my final value //0.3033; //0.55;// my value 0.7489; // tight 0.3033; //deepflav medium
   
   deepak8_cut = 0.895;    // 2017 0.5% mistag rate
   deepak8_cut_md = 0.578; // 2017 0.5% mistag rate
   
   for(int ipt=0; ipt<noAK8ptbins; ipt++){
		DeepAK8_SF[ipt] = DeepAK8_SF_17[ipt];
		DeepAK8_MD_SF[ipt] = DeepAK8_MD_SF_17[ipt];
		DeepAK8_MD_SF_up[ipt] = DeepAK8_MD_SF_17_up[ipt];
		DeepAK8_MD_SF_dn[ipt] = DeepAK8_MD_SF_17_dn[ipt];
	}

#endif

#ifdef Anal_2018
   btagvalue = 0.8838; 
   btagvalue_deepCSV = 0.4184; //0.7527  ;//tight  0.4184; //medium // deepcsv 
   btagvalue_deepFlavB = 0.6;  //0.7264; ;//tight   0.2770;  // medium	//deepflav
   
   deepak8_cut = 0.898;      // 2018 0.5% mistag rate
   deepak8_cut_md = 0.559;   // 2018 0.5% mistag rate
   
   for(int ipt=0; ipt<noAK8ptbins; ipt++){
		DeepAK8_SF[ipt] = DeepAK8_SF_18[ipt];
		DeepAK8_MD_SF[ipt] = DeepAK8_MD_SF_18[ipt];
		DeepAK8_MD_SF_up[ipt] = DeepAK8_MD_SF_18_up[ipt];
		DeepAK8_MD_SF_dn[ipt] = DeepAK8_MD_SF_18_dn[ipt];
	}
   
#endif
 
 
}

Bool_t Anal_Nano_Muon::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either Anal_Nano_Muon::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

 GetEntry(entry);

 TString str;
 
  
 unsigned int fmu = 0;
   
 for(unsigned int imu=0; imu < nMuon; imu++){
	   
	if(!(Muon_tightId[imu])) continue;
//	if(Muon_pfIsoId[imu]<4) continue;
//  switch off isolation cut for non-isolated muons
//	if(Muon_pfRelIso04_all[imu] >= 0.15) continue;
	
	if(fabs(Muon_pt[imu]) < lepptcut) continue;
	if(fabs(Muon_eta[imu]) > 2.4) continue;
	if(fabs(Muon_dxy[imu])>0.045 || fabs(Muon_dz[imu])>0.2) continue;
	
	// 2d iso //
	float dR_min = 1000;
	int nearjet = -1;
	for(unsigned kjet=0; kjet<nJet; kjet++){
		if(Jet_jetId[kjet]==0) continue;
		if(delta2R(Jet_eta[kjet],Jet_phi[kjet],Muon_eta[imu],Muon_phi[imu]) < dR_min){
			dR_min = delta2R(Jet_eta[kjet],Jet_phi[kjet],Muon_eta[imu],Muon_phi[imu]) ;
			nearjet = kjet;
		}
	}
	 
	Muon_pt_nearjet[imu] = 0;
	
	if(nearjet>=0){
		TLorentzVector mu_mom; mu_mom.SetPtEtaPhiM(Muon_pt[imu],Muon_eta[imu],Muon_phi[imu],Muon_mass[imu]);
		TLorentzVector j_mom; j_mom.SetPtEtaPhiM(Jet_pt[nearjet],Jet_eta[nearjet],Jet_phi[nearjet],Jet_mass[nearjet]);
		Muon_pt_nearjet[imu] = ((mu_mom.Vect()).Perp(j_mom.Vect()));
	 }
	
	if(!(dR_min > 0.4 ||  Muon_pt_nearjet[imu] > 25.)) continue;
	// 2d iso ends //
	   
	Muon_pt[fmu] = Muon_pt[imu];
	Muon_eta[fmu] = Muon_eta[imu];
	Muon_phi[fmu] = Muon_phi[imu];
	Muon_mass[fmu] = Muon_mass[imu];
	Muon_charge[fmu] = Muon_charge[imu];
	Muon_jetIdx[fmu] = Muon_jetIdx[imu];
	Muon_tightId[fmu] = Muon_tightId[imu];
	Muon_pfIsoId[fmu] = Muon_pfIsoId[imu];
	Muon_pt_nearjet[fmu] = Muon_pt_nearjet[imu];
	   
	fmu++; 
	if(fmu>=njetmax) break;	
	
   }
	   
   nMuon = fmu;
    
 unsigned int fel = 0;
   
 for(unsigned iel=0; iel < nElectron; iel++){
	   
//  if(Electron_cutBased[iel]!=4) continue;
	if(Electron_mvaFall17V2Iso_WP80[iel]==0) continue;
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
	if(fabs(Photon_pt[ipho]) < 30 || fabs(Photon_eta[ipho])>2.4) continue;
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
	  
//	   if(ijet==0 &&  (FatJet_pt[ijet]<200. || abs(FatJet_eta[ijet])>4.0)) break;
 
	   if(FatJet_jetId[ijet]==0) continue;
	   int jetid = *(decToBinary(FatJet_jetId[ijet])+1);//((decToBinary(FatJet_jetId[ijet]))/10)%10;  // tight id
	   if(jetid==0) continue;
	   
	   bool mutag = false;
	   int close_mu = -1;
	   for(unsigned imu=0; imu<nMuon; imu++){
		   if(delta2R(Muon_eta[imu],Muon_phi[imu],FatJet_eta[ijet],FatJet_phi[ijet])<0.8){
			   mutag = true;
			   close_mu = imu;
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
	   
	   if(/*mutag||*/etag) continue;

	   #ifdef NoJEC
	   FatJet_mass[fjet] = FatJet_mass[ijet];
	   FatJet_msoftdrop[fjet] = FatJet_msoftdrop[ijet];
	   FatJet_pt[fjet] = FatJet_pt[ijet];
	   #else
	   FatJet_mass[fjet] = FatJet_mass_nom[ijet];
	   if(isMC){
	   FatJet_msoftdrop[fjet] = FatJet_msoftdrop_raw[ijet]*FatJet_corr_JER[ijet]*FatJet_msoftdrop_corr_JMS[ijet]*FatJet_msoftdrop_corr_JMR[ijet];
	   }else{
		    FatJet_msoftdrop[fjet] = FatJet_msoftdrop_raw[ijet]*FatJet_msoftdrop_corr_JMS[ijet];
		    }
	   FatJet_pt[fjet] = FatJet_pt_nom[ijet];
	   #endif
	   FatJet_eta[fjet] = FatJet_eta[ijet];
	   FatJet_phi[fjet] = FatJet_phi[ijet];
	   
	   if(mutag){
		   TLorentzVector jetmom; 
		   jetmom.SetPtEtaPhiM(FatJet_pt[fjet],FatJet_eta[fjet],FatJet_phi[fjet],FatJet_mass[fjet]);
		   TLorentzVector mumom; mumom.SetPtEtaPhiM(Muon_pt[close_mu],Muon_eta[close_mu],Muon_phi[close_mu],Muon_mass[close_mu]);
		   jetmom = jetmom - mumom;
		   FatJet_pt[fjet] = jetmom.Pt();
		   FatJet_mass[fjet] = jetmom.M();
		   FatJet_eta[fjet] = jetmom.Eta();
		   FatJet_phi[fjet] =	jetmom.Phi();
		   }

	   if(FatJet_pt[fjet]<30. || abs(FatJet_eta[fjet])>2.4) break;

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
	   FatJet_n2b1[fjet] = FatJet_n2b1[ijet];
	   FatJet_n3b1[fjet] = FatJet_n3b1[ijet];
	   FatJet_rawFactor[fjet] = FatJet_rawFactor[ijet];
	   FatJet_tau1[fjet] = FatJet_tau1[ijet];
	   FatJet_tau2[fjet] = FatJet_tau2[ijet];
	   FatJet_tau3[fjet] = FatJet_tau3[ijet];
	   FatJet_tau4[fjet] = FatJet_tau4[ijet];
	   FatJet_jetId[fjet] = FatJet_jetId[ijet];
	   FatJet_subJetIdx1[fjet] = FatJet_subJetIdx1[ijet];
	   FatJet_subJetIdx2[fjet] = FatJet_subJetIdx2[ijet];
	   
	   float mindR = 0.4;
	   int matchgen = -1;
	   if(isMC){
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
	   
	   fjet++;
	   if(fjet>=njetmax) break;
	   }
   
   nFatJet = fjet;
   
   
   // sort by pt AK8 //
   
   for(unsigned ijet=0; ijet < nFatJet; ijet++){
	 for(unsigned kjet=(ijet+1); kjet < nFatJet; kjet++){
		   
		if(FatJet_pt[kjet]>FatJet_pt[ijet]){
			
			float tmppt = FatJet_pt[ijet];
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
			FatJet_btagHbb[ijet] = FatJet_btagDeepB[kjet];
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
	}//ijrt
   
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
	
	   if(Jet_jetId[ijet]==0) continue;
	
	   int jetid = *(decToBinary(Jet_jetId[ijet])+1);//((decToBinary(Jet_jetId[ijet]))/10)%10; // tight id
	   if(jetid==0) continue;
	  
	   bool mutag = false; int close_mu = -1;
	   for(unsigned imu=0; imu<nMuon; imu++){
		   if(delta2R(Muon_eta[imu],Muon_phi[imu],Jet_eta[ijet],Jet_phi[ijet])<0.4){
			   mutag = true;
			   close_mu = imu;
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
	   
	   if(/*mutag||*/etag) continue;
	   
	   #ifdef NoJEC
	   Jet_mass[fjet] = Jet_mass[ijet];
	   Jet_pt[fjet] = Jet_pt[ijet];
	   #else
	   Jet_mass[fjet] = Jet_mass_nom[ijet];
	   Jet_pt[fjet] = Jet_pt_nom[ijet];
	   #endif
	   Jet_eta[fjet] = Jet_eta[ijet];
	   Jet_phi[fjet] = Jet_phi[ijet];
	   
	   if(mutag){
		   TLorentzVector jetmom; 
		   jetmom.SetPtEtaPhiM(Jet_pt[fjet],Jet_eta[fjet],Jet_phi[fjet],Jet_mass[fjet]);
		   TLorentzVector mumom; mumom.SetPtEtaPhiM(Muon_pt[close_mu],Muon_eta[close_mu],Muon_phi[close_mu],Muon_mass[close_mu]);
		   jetmom = jetmom - mumom;
		   Jet_pt[fjet] = jetmom.Pt();
		   Jet_mass[fjet] = jetmom.M();
		   Jet_eta[fjet] = jetmom.Eta();
		   Jet_phi[fjet] =	jetmom.Phi();
		}
	
	   if(Jet_pt[fjet]<30.) continue;
	   if(fabs(Jet_eta[fjet])>2.4) continue;
	    
	   Jet_area[fjet] = Jet_area[ijet];
	   Jet_btagCMVA[fjet] = Jet_btagCMVA[ijet];
	   Jet_btagCSVV2[fjet] = Jet_btagCSVV2[ijet];
	   Jet_btagDeepB[fjet] = Jet_btagDeepB[ijet];
	   Jet_btagDeepC[fjet] = Jet_btagDeepC[ijet];
	   Jet_btagDeepFlavB[fjet] = Jet_btagDeepFlavB[ijet];
	   Jet_chEmEF[fjet] = Jet_chEmEF[ijet];
	   Jet_chHEF[fjet] = Jet_chHEF[ijet];
	   Jet_muEF[fjet] = Jet_muEF[ijet];
	   Jet_neEmEF[fjet] = Jet_neEmEF[ijet];
	   Jet_neHEF[fjet] = Jet_neHEF[ijet];
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
	   Jet_DeepCSVMCeff[fjet] = BTag_MCEfficiency("deepcsvm",int(Jet_hadronFlavour[fjet]),Jet_pt[fjet],fabs(Jet_eta[fjet]));
	   
	   fjet++;
	   if(fjet>=njetmax) break;
	   }
   
   nJet = fjet;
   
   // sort by pt AK4 //
   
   for(unsigned ijet=0; ijet < nJet; ijet++){
	 for(unsigned kjet=(ijet+1); kjet < nJet; kjet++){
		   
		if(Jet_pt[kjet]>Jet_pt[ijet]){
			
			float tmppt = Jet_pt[ijet];
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
			
			tmppt = Jet_DeepCSVMCeff[ijet];
			Jet_DeepCSVMCeff[ijet] = Jet_DeepCSVMCeff[kjet];
			Jet_DeepCSVMCeff[kjet] = tmppt;
			
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
   
   
  int nbjetsAK4 = 0;
   
   #ifdef Btagger_DeepJet
   
   for(unsigned ijet=0; ijet<nJet; ijet++){
			 
	      if(fabs(Jet_eta[ijet]) < 2.4 && (Jet_pt[ijet]>30.)){
			if(Jet_btagDeepFlavB[ijet] > btagvalue_deepFlavB){
			  nbjetsAK4++;
			}
		 }
	  }
	  
  #endif	
	
  #ifdef Btagger_DeepCSV
  
  for(unsigned ijet=0; ijet<nJet; ijet++){
			 
	      if(fabs(Jet_eta[ijet]) < 2.4 && (Jet_pt[ijet]>30.)){
			if(Jet_btagDeepB[ijet] > btagvalue_deepCSV){
			  nbjetsAK4++;
			}
		 }
	  }
  
  #endif	
  
  #ifdef Btagger_CSVv2
  
  for(unsigned ijet=0; ijet<nJet; ijet++){
			 
	      if(fabs(Jet_eta[ijet]) < 2.4 && (Jet_pt[ijet]>30.)){
			if(Jet_btagCSVV2[ijet] > btagvalue){
			  nbjetsAK4++;
			}
		 }
	  }
  
  #endif 
   
  Event_Ht = 0;    
  for(unsigned ijet=0; ijet<nJet; ijet++){
	 if(Jet_pt[ijet] > 30.){
		Event_Ht += Jet_pt[ijet];
	}
  }
     
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
   if(isMC){
      weight = Generator_weight;
      nevent_total += 1;
      weightev = Generator_weight;
	  
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
   
   if(usePrefireWeight){
	#if defined(Anal_2016) || defined(Anal_2017) 
	weight *= PrefireWeight;
	#endif
    }
   
   Tout->Fill();
   
   hist_npv->Fill(PV_npvsGood,weight);
   
   if(PV_npvsGood < 1) return kFALSE;
   
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
      
   if(!metFilterpassed) return kFALSE;
   
   bool nisomu_trig = HLT_Mu50;//HLT_IsoMu27 ;
   bool mu_trig = (HLT_IsoMu27 || HLT_Mu50);
   
   #if defined(Anal_2017) || defined(Anal_2018)
   bool had_trig = (HLT_AK8PFJet420_TrimMass30 || HLT_AK8PFHT900_TrimMass50 || HLT_PFJet500 || HLT_AK8PFJet500 || HLT_PFHT1050);
   #endif
   
   #ifdef Data_2017B
    had_trig = (HLT_PFJet500 || HLT_AK8PFJet500 || HLT_PFHT1050);
   #endif
   
   #ifdef Anal_2016
   bool had_trig = (HLT_AK8DiPFJet280_200_TrimMass30 || HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20 || 
                    HLT_AK8DiPFJet300_200_TrimMass30 || HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20 ||
					HLT_AK8PFHT700_TrimR0p1PT0p03Mass50 || HLT_AK8PFJet360_TrimMass30 ||
					HLT_PFJet450 || HLT_AK8PFJet450 || HLT_PFHT800 || HLT_PFHT900);
 //  if(isMC){
 //	   had_trig = (had_trig || HLT_AK8DiPFJet300_200_TrimMass30);
 //	   }					
   #endif
   
   #ifdef Data_2016H
		had_trig = (HLT_AK8DiPFJet300_200_TrimMass30 || HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20 ||
					HLT_AK8PFHT700_TrimR0p1PT0p03Mass50 || HLT_AK8PFJet360_TrimMass30 ||
					HLT_PFJet450 || HLT_AK8PFJet450 || HLT_PFHT900);
   #endif					
   
   bool el_trig;
   
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
		   			
			hist_ht_emutrig->Fill(Event_Ht,weight);
			hist_pt_emutrig->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_emutrig->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_emutrig->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_emutrig->Fill(ptb.M(),weight);
			}
		
		#if defined(Anal_2017) || defined(Anal_2018)
			
		if(HLT_PFHT1050 ){
			hist_ht_HT1050_wemutrig->Fill(Event_Ht,weight);
			hist_pt_HT1050_wemutrig->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_HT1050_wemutrig->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_HT1050_wemutrig->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_HT1050_wemutrig->Fill(ptb.M(),weight);
				}
			}	
		if(HLT_PFJet500){
			hist_ht_AK4Pt500_wemutrig->Fill(Event_Ht,weight);
			hist_pt_AK4Pt500_wemutrig->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK4Pt500_wemutrig->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK4Pt500_wemutrig->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK4Pt500_wemutrig->Fill(ptb.M(),weight);
				}
			}
		if(HLT_AK8PFJet500){
			hist_ht_AK8Pt500_wemutrig->Fill(Event_Ht,weight);
			hist_pt_AK8Pt500_wemutrig->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK8Pt500_wemutrig->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK8Pt500_wemutrig->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK8Pt500_wemutrig->Fill(ptb.M(),weight);
				}
			}
		if(HLT_AK8PFJet420_TrimMass30){
			hist_ht_AK8PFJet420_TrimMass30_wemutrig->Fill(Event_Ht,weight);
			hist_pt_AK8PFJet420_TrimMass30_wemutrig->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK8PFJet420_TrimMass30_wemutrig->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK8PFJet420_TrimMass30_wemutrig->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK8PFJet420_TrimMass30_wemutrig->Fill(ptb.M(),weight);
				}
			}
			
		if(HLT_AK8PFHT900_TrimMass50){
			hist_ht_AK8PFHT900_TrimMass50_wemutrig->Fill(Event_Ht,weight);
			hist_pt_AK8PFHT900_TrimMass50_wemutrig->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK8PFHT900_TrimMass50_wemutrig->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK8PFHT900_TrimMass50_wemutrig->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK8PFHT900_TrimMass50_wemutrig->Fill(ptb.M(),weight);
				}
			}	
			
		#endif
		
		#ifdef Anal_2016
		
		#ifdef Data_2016H
		
		if(HLT_PFHT900){
			hist_ht_HT1050_wemutrig->Fill(Event_Ht,weight);
			hist_pt_HT1050_wemutrig->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_HT1050_wemutrig->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_HT1050_wemutrig->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_HT1050_wemutrig->Fill(ptb.M(),weight);
				}
			}	
		
		#else
		
		if(HLT_PFHT800){
			hist_ht_HT1050_wemutrig->Fill(Event_Ht,weight);
			hist_pt_HT1050_wemutrig->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_HT1050_wemutrig->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_HT1050_wemutrig->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_HT1050_wemutrig->Fill(ptb.M(),weight);
				}
			}	
		
		#endif	
			
		if(HLT_PFJet450){
			hist_ht_AK4Pt500_wemutrig->Fill(Event_Ht,weight);
			hist_pt_AK4Pt500_wemutrig->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK4Pt500_wemutrig->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK4Pt500_wemutrig->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK4Pt500_wemutrig->Fill(ptb.M(),weight);
				}
			}
		if(HLT_AK8PFJet450){
			hist_ht_AK8Pt500_wemutrig->Fill(Event_Ht,weight);
			hist_pt_AK8Pt500_wemutrig->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK8Pt500_wemutrig->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK8Pt500_wemutrig->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK8Pt500_wemutrig->Fill(ptb.M(),weight);
				}
			}
		if(HLT_AK8PFJet360_TrimMass30){
			hist_ht_AK8PFJet420_TrimMass30_wemutrig->Fill(Event_Ht,weight);
			hist_pt_AK8PFJet420_TrimMass30_wemutrig->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK8PFJet420_TrimMass30_wemutrig->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK8PFJet420_TrimMass30_wemutrig->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK8PFJet420_TrimMass30_wemutrig->Fill(ptb.M(),weight);
				}
			}
		
		if(HLT_AK8PFHT700_TrimR0p1PT0p03Mass50){
			hist_ht_AK8PFHT900_TrimMass50_wemutrig->Fill(Event_Ht,weight);
			hist_pt_AK8PFHT900_TrimMass50_wemutrig->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK8PFHT900_TrimMass50_wemutrig->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK8PFHT900_TrimMass50_wemutrig->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK8PFHT900_TrimMass50_wemutrig->Fill(ptb.M(),weight);
				}
			}
		
		#ifdef Data_2016H
		
		if( HLT_AK8DiPFJet300_200_TrimMass30 || HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20) {
			hist_ht_AK8DiPFJet300_200_TrimMass30_wemutrig->Fill(Event_Ht,weight);
			hist_pt_AK8DiPFJet300_200_TrimMass30_wemutrig->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK8DiPFJet300_200_TrimMass30_wemutrig->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK8DiPFJet300_200_TrimMass30_wemutrig->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK8DiPFJet300_200_TrimMass30_wemutrig->Fill(ptb.M(),weight);
				}
			}
			
		#else	
		
		if( HLT_AK8DiPFJet280_200_TrimMass30 || HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20) {
			hist_ht_AK8DiPFJet300_200_TrimMass30_wemutrig->Fill(Event_Ht,weight);
			hist_pt_AK8DiPFJet300_200_TrimMass30_wemutrig->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_AK8DiPFJet300_200_TrimMass30_wemutrig->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_AK8DiPFJet300_200_TrimMass30_wemutrig->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_AK8DiPFJet300_200_TrimMass30_wemutrig->Fill(ptb.M(),weight);
				}
			}
		
		#endif
		
		#endif
				
		if(had_trig){
			hist_ht_all_wemutrig->Fill(Event_Ht,weight);
			hist_pt_all_wemutrig->Fill(FatJet_pt[topjetAK8_preli],weight);
			hist_bpt_all_wemutrig->Fill(Jet_pt[bjetAK4_preli],weight);
			hist_ptsum_all_wemutrig->Fill(twoptsum,weight);
			if(FatJet_pt[topjetAK8_preli]>AK8ptcut && Jet_pt[bjetAK4_preli]>AK4ptcut_fi){
				hist_mtb_all_wemutrig->Fill(ptb.M(),weight);
				}
			}	
		}
     }
   // trigger fill end //
   
	bool trig_match = false;
    
    for(unsigned nmu=0; nmu<nMuon; nmu++){
		for(unsigned itr=0; itr < nTrigObj; itr++){
		float delR = delta2R(TrigObj_eta[itr],TrigObj_phi[itr],Muon_eta[nmu],Muon_phi[nmu]);
		if((delR < 0.4) && (abs(TrigObj_id[itr])==13) /*&& (TrigObj_l1iso[itr]>0)*/ && (TrigObj_filterBits[itr] & (1<<3)) == (1<<3)){
			trig_match = true;
			break;
			}
		}
		if(trig_match) { break; }
	 }
	 
     
   hist_nmuons->Fill(nMuon,weight);
   hist_nelectrons->Fill(nElectron,weight);
   hist_nphotons->Fill(nPhoton,weight);
   hist_njetAK4->Fill(nJet,weight);
   hist_nbjetAK4->Fill(nbjetsAK4,weight);
   hist_njetAK8->Fill(nFatJet,weight);
   hist_trig_match->Fill(trig_match,weight);
   
//   str = TString::Format("entry %lli nMuon %i nElectron %i nPhoton %i nisomu_trig %o trig_match %i nJet %i nbjetsAK4 %i\n",entry,nMuon,nElectron,nPhoton,nisomu_trig,trig_match,nJet,nbjetsAK4);
//   if(gProofServ) gProofServ->SendAsynMessage(str);
   
   if(nElectron==0&&nPhoton==0&&nisomu_trig&&trig_match&&(nJet>=2)&&(nbjetsAK4>=2)) { hist_nmuons_cut->Fill(nMuon,weight); }
   if(nMuon==1&&nPhoton==0&&nisomu_trig&&trig_match&&(nJet>=2)&&(nbjetsAK4>=2)) { hist_nelectrons_cut->Fill(nElectron,weight); }
   if(nMuon==1&&nElectron==0&&nisomu_trig&&trig_match&&(nJet>=2)&&(nbjetsAK4>=2)) { hist_nphotons_cut->Fill(nPhoton,weight); }
   if(nMuon==1&&nElectron==0&&nPhoton==0&&trig_match&&(nJet>=2)&&(nbjetsAK4>=2)) { hist_mutrig_pass_cut->Fill(nisomu_trig,weight); }
   if(nMuon==1&&nElectron==0&&nPhoton==0&&nisomu_trig&&(nJet>=2)&&(nbjetsAK4>=2)) { hist_trig_match_cut->Fill(trig_match,weight); }
   if(nMuon==1&&nElectron==0&&nPhoton==0&&nisomu_trig&&trig_match&&(nbjetsAK4>=2)) {  hist_njetAK4_cut->Fill(nJet,weight); }
   if(nMuon==1&&nElectron==0&&nPhoton==0&&nisomu_trig&&trig_match&&(nJet>=2)) { hist_nbjetAK4_cut->Fill(nbjetsAK4,weight); }
      
   if(nMuon!=1)  return kFALSE;
   if(nElectron  > 0) return kFALSE;
   if(nPhoton > 0) return kFALSE;
   
//   if(!(had_trig)) return kFALSE; 
   if(!(nisomu_trig)) return kFALSE; 
   if(!isSignal){
	if(!trig_match) return kFALSE;
   }
   
   if(nJet<2) return kFALSE;
   if(nbjetsAK4<2) return kFALSE;
   if(nFatJet<1) return kFALSE;
   
   if(isMC){
	
	// Muon SF //
	
	int imuetabin = -1;
	int imuptbin = -1;
	int imuptbin_trig = -1;
	
	imuetabin = getbinid(fabs(Muon_eta[0]),muetabins,mueta_bins);
	imuptbin = getbinid(fabs(Muon_pt[0]),muptbins,mupt_bins);
	#ifdef Anal_2017	
	imuptbin_trig = getbinid(fabs(Muon_pt[0]),muptbins_trg,mupt_bins_trig);
	#endif
	#ifdef Anal_2018	
	imuptbin_trig = getbinid(fabs(Muon_pt[0]),muptbins_trg_18,mupt_bins_trig_18);
	#endif
	
	double mu_weight = 1;
	
	if(imuetabin>=0 & imuptbin>=0 & imuptbin_trig>=0)
	{
	#ifdef Anal_2017	
	mu_weight = muid_SF[imuetabin][imuptbin] * muiso_SF[imuetabin][imuptbin] * mutrg_SF[imuetabin][imuptbin_trig];
	#elif defined(Anal_2018)
	mu_weight = muid_SF18[imuetabin][imuptbin] * muiso_SF18[imuetabin][imuptbin] * mutrg_SF18[imuetabin][imuptbin_trig];
	#endif
	}
	
	weight *= mu_weight;
	
//	str = TString::Format("entry %lli mu eta %f pt %f sf %f\n",entry,Muon_eta[0],Muon_pt[0],mu_weight);
 //   if(gProofServ) gProofServ->SendAsynMessage(str);
	
	 // shape based discriminator //
	
	weightpass = weight;
	Tout1->Fill();
	 
	#ifndef NoJEC
		float mcprod = 1.;
		float dataprod = 1.;
		float btag_weight = 1.;
	
		unsigned int njets1 = nJet;//(nJet>=4)?4:nJet; 
	 
	 #ifdef Btagger_DeepJet_wt
	  for(unsigned ijet=0; ijet<njets1; ijet++){
			
			if(!(fabs(Jet_eta[ijet]) < 2.4 && (Jet_pt[ijet]>30. && Jet_pt[ijet]<1000.))) continue;
			btag_weight *= Jet_deepflavbtagSF_shape[ijet];
		}
	
		weight *= 	btag_weight;
	 #endif	
	 
	 #ifdef Btagger_DeepCSV_wt
	  for(unsigned ijet=0; ijet<njets1; ijet++){
			
			if(!(fabs(Jet_eta[ijet]) < 2.4 && (Jet_pt[ijet]>30. && Jet_pt[ijet]<1000.))) continue;
			btag_weight *= Jet_deepCSVbtagSF_shape[ijet];
		}
	
		weight *= 	btag_weight;
	 #endif	
	 
	 #ifdef Btagger_CSVv2
	  for(unsigned ijet=0; ijet<njets1; ijet++){
			
			if(!(fabs(Jet_eta[ijet]) < 2.4 && (Jet_pt[ijet]>30. && Jet_pt[ijet]<1000.))) continue;
			btag_weight *= Jet_CSVbtagSF_shape[ijet];
		}
	
		weight *= 	btag_weight;
	 #endif	
	 	
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
	 
   if(MET_pt < 50.) return kFALSE; 
   
   hist_npv_final->Fill(PV_npvsGood,weight);
   
// quarks from top //

int top_d[6];
int ist=0; int isqg = 0;
int qp[4] = {-1,-1,-1,-1};
int bp[2] = {-1,-1};
int igtop[2] = {-1,-1};
int ngtop = 0; int nghtop = 0;
int ihtop[2] = {-1,-1};
int qg_d[6]={-1,-1,-1,-1,-1,-1};
TLorentzVector topparton_vec4;


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

SF_toppt = 1;
SF_flat_toppt = 1;

if(ngtop==2 && isTTBar && isMC){
	SF_toppt =  sqrt(exp(0.0615-0.0005*GenPart_pt[igtop[0]]) * exp(0.0615-0.0005*GenPart_pt[igtop[1]]));
	weight *= SF_toppt;
    SF_flat_toppt =  sqrt(exp(0.0615-0.0005*TMath::Min(float(400.),GenPart_pt[igtop[0]])) * exp(0.0615-0.0005*TMath::Min(float(400.),GenPart_pt[igtop[1]])))*sqrt(EW_toppt_cor(GenPart_pt[igtop[0]])*EW_toppt_cor(GenPart_pt[igtop[1]]));
//  weight *= SF_toppt;
}

weight_t = weight;

// deepak8 scale factor //

sfwt_deepak8 = 1.;
sfwt_deepak8_md = 1.;

if(isMC){
	
	float ttag_eff = 1;
	if(isTTBar||isST){
		ttag_eff = 0.65;
	}else{
		ttag_eff = 0.05;
	}
	
	for(unsigned ijet=0; ijet<nFatJet; ijet++){
		if(FatJet_pt[ijet] < AK8ptcut) continue;
		int iptbin = getbinid(FatJet_pt[ijet],noAK8ptbins,AK8ptbins);
		if(iptbin>=0){
			if(FatJet_deepTagMD_TvsQCD[ijet] >= deepak8_cut_md){
				sfwt_deepak8_md *= DeepAK8_MD_SF[iptbin];
			}else{
				  sfwt_deepak8_md *= (1.-ttag_eff*DeepAK8_MD_SF[iptbin])*1./(1.-ttag_eff);
			  }
		}
	}
	
	weight_t *= sfwt_deepak8_md ; 
	
}
	  
//	  weight_t = weight;  

int topmu = 0;
	  
bool b_found = false;  	  
    
if(topmu>=0){
      
  for(unsigned ijet=0; ijet<nJet; ijet++){
			  
	double dR = delta2R(Muon_eta[topmu],Muon_phi[topmu],Jet_eta[ijet],Jet_phi[ijet]);
	double dEta = (Jet_eta[ijet]-Muon_eta[topmu]);
	double dPhi = PhiInRange(Jet_phi[ijet]-Muon_phi[topmu]);	
			  
	if(dR>dR_cut){
		if(fabs(dPhi) > dphi_cut){	
					
			TLorentzVector AK4vec4; 
			AK4vec4.SetPtEtaPhiM(Jet_pt[ijet],Jet_eta[ijet],Jet_phi[ijet],Jet_mass[ijet]);
			    	  
			hist_jetptAK4->Fill(Jet_pt[ijet],weight);
			hist_jetrapAK4->Fill(AK4vec4.Rapidity(),weight);
					
			if(fabs(AK4vec4.Rapidity()) < 2.4 && (Jet_pt[ijet]>AK4ptcut_in)){
						
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

 bool top_found = false;

 for(unsigned ijet=0; ijet<nFatJet; ijet++){
			  
	double dR = delta2R(Muon_eta[topmu],Muon_phi[topmu],FatJet_eta[ijet],FatJet_phi[ijet]);
	double dEta = (FatJet_eta[ijet]-Muon_eta[topmu]);
	double dPhi = PhiInRange(FatJet_phi[ijet]-Muon_phi[topmu]);	
			  
	if(dR>dR_cut){
		if(fabs(dPhi) > dphi_cut){	
			
			TLorentzVector AK4vec8; 
			AK4vec8.SetPtEtaPhiM(FatJet_pt[ijet],FatJet_eta[ijet],FatJet_phi[ijet],FatJet_mass[ijet]);
			    	  
			if(fabs(AK4vec8.Rapidity()) < 2.4 && (AK4vec8.Pt()>AK8ptcut)){
						
			if(!top_found){						
				topjetAK8 = ijet;
				top_found = true;
				break;
			}
		 }
	  }  
   }
 }
		  
}//topmu	  
	
	
if(topjetAK8>=0){
		
	hist_tiso_mass->Fill(FatJet_msoftdrop[topjetAK8],weight);
	if(FatJet_deepTagMD_TvsQCD[topjetAK8] > deepak8_cut_md){
		hist_tiso_mass_deepak8cut->Fill(FatJet_msoftdrop[topjetAK8],weight_t);
		hist_tiso_mass_deepak8cut_flattoppt->Fill(FatJet_msoftdrop[topjetAK8],weight_t*SF_flat_toppt/SF_toppt);
	}
	
	bool had_w_found = false;
	for(unsigned ijet=0; ijet<nJet; ijet++){
		for(unsigned fjet=0; fjet<nJet; fjet++){
			TLorentzVector vec1; vec1.SetPtEtaPhiM(Jet_pt[ijet],Jet_eta[ijet],Jet_phi[ijet],Jet_mass[ijet]);
			TLorentzVector vec2; vec2.SetPtEtaPhiM(Jet_pt[fjet],Jet_eta[fjet],Jet_phi[fjet],Jet_mass[fjet]);
			if((vec1+vec2).M()>wmasslow && (vec1+vec2).M()<wmasshigh){
				had_w_found = true;
				break;
				}
		}
		if(had_w_found) break;
	}
	if(had_w_found){
			hist_tiso_mass_mwcut->Fill(FatJet_msoftdrop[topjetAK8],weight);
		}
}
 		
	if(bjetAK4>=0){	
			
		TLorentzVector AK4vec4; 
		AK4vec4.SetPtEtaPhiM(Jet_pt[bjetAK4],Jet_eta[bjetAK4],Jet_phi[bjetAK4],Jet_mass[bjetAK4]);
					
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
			if((Jet_btagDeepFlavB[bjetAK4] > btagvalue_deepFlavB)){
				bjet_wp = 1;
			}
			else { bjet_wp = 0; }
			#endif
				  
			#ifdef Btagger_DeepCSV
			if((Jet_btagDeepB[bjetAK4] > btagvalue_deepCSV)){
				bjet_wp = 1;
			}
			else { bjet_wp = 0; }
			#endif
				  
			#ifdef Btagger_CSVv2
			if(Jet_btagCSVV2[bjetAK4] > btagvalue){
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
				hist_biso_TopAK8score_MD->Fill(FatJet_deepTagMD_TvsQCD[AK8matchb],weight_t);
				hist_biso_WAK8score->Fill(FatJet_deepTag_WvsQCD[AK8matchb],weight);
				hist_biso_WAK8score_MD->Fill(FatJet_deepTagMD_WvsQCD[AK8matchb],weight);
						
					
				if(isMC){
					if(FatJet_MatchGenJet[AK8matchb]>=0){
						hist_bjetAK8genmass->Fill(FatJet_GenJetmass[AK8matchb],weight);
					}
				}
					
				if(Jet_pt[bjetAK4] > AK4ptcut_fi){
						
					double dR = delta2R(Muon_eta[topmu],Muon_phi[topmu],Jet_eta[bjetAK4],Jet_phi[bjetAK4]);
					double dPhi = PhiInRange(Jet_phi[bjetAK4]-Muon_phi[topmu]);
					
					hist_delR_btag_toptag->Fill(dR,weight);
					hist_delphi_btag_toptag->Fill(dPhi,weight);
					
					hist_biso_mass_1->Fill(FatJet_msoftdrop[AK8matchb],weight_t);
					if(topjetAK8>=0){
					if(FatJet_msoftdrop[topjetAK8]>topmasslow && FatJet_msoftdrop[topjetAK8]<topmasshigh){
						hist_biso_mass_msdcut_1->Fill(FatJet_msoftdrop[AK8matchb],weight_t);
						}
					}
					if(FatJet_deepTagMD_TvsQCD[AK8matchb] > deepak8_cut_md){
						hist_biso_mass_deepak8cut_1->Fill(FatJet_msoftdrop[AK8matchb],weight_t);
						}
					hist_biso_tbmassdiff_1->Fill((matchtjet-AK4vec4).M(),weight_t);
					hist_biso_isomass_1->Fill(AK4vec4.M()/FatJet_msoftdrop[AK8matchb],weight_t);
					hist_biso_isopt_1->Fill(Jet_pt[bjetAK4]/FatJet_pt[AK8matchb],weight_t);
					hist_biso_TopAK8score_1->Fill(FatJet_deepTag_TvsQCD[AK8matchb],weight_t);
					hist_biso_TopAK8score_MD_1->Fill(FatJet_deepTagMD_TvsQCD[AK8matchb],weight_t);
					hist_biso_WAK8score_1->Fill(FatJet_deepTag_WvsQCD[AK8matchb],weight_t);
					hist_biso_WAK8score_MD_1->Fill(FatJet_deepTagMD_WvsQCD[AK8matchb],weight_t);
				
					hist_biso_mass_wb_1[bjet_wp]->Fill(FatJet_msoftdrop[AK8matchb],weight_t);
					hist_biso_tbmassdiff_wb_1[bjet_wp]->Fill((matchtjet-AK4vec4).M(),weight_t);
					hist_biso_isomass_wb_1[bjet_wp]->Fill(AK4vec4.M()/FatJet_msoftdrop[AK8matchb],weight_t);
					hist_biso_isopt_wb_1[bjet_wp]->Fill(Jet_pt[bjetAK4]/FatJet_pt[AK8matchb],weight_t);
					hist_biso_TopAK8score_wb_1[bjet_wp]->Fill(FatJet_deepTag_TvsQCD[AK8matchb],weight_t);
					hist_biso_TopAK8score_MD_wb_1[bjet_wp]->Fill(FatJet_deepTagMD_TvsQCD[AK8matchb],weight_t);
					hist_biso_WAK8score_wb_1[bjet_wp]->Fill(FatJet_deepTag_WvsQCD[AK8matchb],weight_t);
					hist_biso_WAK8score_MD_wb_1[bjet_wp]->Fill(FatJet_deepTagMD_WvsQCD[AK8matchb],weight_t);
					
					if(isMC){
					if(FatJet_MatchGenJet[AK8matchb]>=0){
						hist_bjetAK8genmass_wb_1[bjet_wp]->Fill(FatJet_GenJetmass[AK8matchb],weight_t);
							}
						}
					
					if(FatJet_msoftdrop[AK8matchb] < bmasscut_fun(FatJet_pt[AK8matchb])){
							hist_bjetpass_wb_1[bjet_wp]->Fill(1.01,weight_t);
						}else{
							hist_bjetfail_wb_1[bjet_wp]->Fill(0.01,weight_t);
							 }
				
					} //pt-cut
				} // match AK8-jet >=0

	  } // bjetAK4>=0
		  
	  	  
   return kTRUE;
}

void Anal_Nano_Muon::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

	fileOut->cd();
	fileOut->Write();

	fOutput->Add(OutFile);

	fileOut->Close();

}

void Anal_Nano_Muon::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
