#define Anal_Nano_PROOF_cxx
//#define NoJEC
#define Btagger_DeepJet
//#define Btagger_DeepCSV
//#define Btagger_CSVv2
#define Btagger_DeepJet_wt
//#define Btagger_DeepCSV_wt

//#define Anal_2017
//#define Data_2017B

#define Anal_2016
//#define Data_2016H

//#define Anal_2018
//#define HEM_Cor_ON

//#define B_Tight

#define LHAPDFX

// The class definition in Anal_Nano_PROOF.h has been generated automatically
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
// Root > T->Process("Anal_Nano_PROOF.C")
// Root > T->Process("Anal_Nano_PROOF.C","some options")
// Root > T->Process("Anal_Nano_PROOF.C+")
//

#include "Anal_Nano_PROOF.h"
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


void Anal_Nano_PROOF::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void Anal_Nano_PROOF::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
   
   double  sdmassbins[nsdmassbins+1] = {0,10,21,33,45,57,69,81,93,105,117,129,141,156,172,190,210,250,290,340,400};

   TString option = GetOption();

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
   FakeAnalysis = true;
   isSignal = false;
   isTTBar = false;
   isTTHad = false; isTTSemiLep = false; isTTDiLep = false;
   isST = false;
   
   TopTagging = true;
   if(!isMC){ TopTagging = false; }
   
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
   hist_2D_sdmass_tau32_AK8 = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
   hist_2D_sdmass_tau32_AK8->Sumw2();
   
   sprintf(name,"H2D_SDMass_subbtag_AK8CHS");
   sprintf(title,"2D Correlation of m_{SD} and subjet b-tag AK8CHS jets");
   hist_2D_sdmass_subbtag_AK8 = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
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
		
		
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_DeepAK8Tag%i_AK4btag%i_AK8CHS",itau+1,ib+1);
		hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
		hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8[ib][itau]->Sumw2();
		
		sprintf(name,"H2D_SDMass_subbtag_CSVv2_DeepAK8Tag%i_AK4btag%i_AK8CHS_hight",itau+1,ib+1);
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
		hist_2D_sdmass_subbtag_CSVv2_AK8_md_DeepAK8_tt_3[ib][itau] = new TH2D(name,title,nsdmassbins,sdmassbins,20,0,1.0);
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
		hist_2D_DeepAK8_tbmass_1[ib] = new TH2D(name,name,20,0,1.0,nvwpmbin,wpmbins);
		hist_2D_DeepAK8_tbmass_1[ib] = new TH2D(name,name,20,0,1.0,nvwpmbin,wpmbins);
		
		sprintf(name,"H2D_DeepAK8_tbMass_AK4btag%i_AK8CHS_nobmasscut",ib+1);
		hist_2D_DeepAK8_tbmass_2[ib] = new TH2D(name,name,20,0,1.0,nvwpmbin,wpmbins);
		hist_2D_DeepAK8_tbmass_2[ib] = new TH2D(name,name,20,0,1.0,nvwpmbin,wpmbins);
		
		sprintf(name,"H2D_DeepAK8_tbMass_AK4btag%i_AK8CHS",ib+1);
		hist_2D_DeepAK8_tbmass[ib] = new TH2D(name,name,20,0,1.0,nvwpmbin,wpmbins);
		hist_2D_DeepAK8_tbmass[ib] = new TH2D(name,name,20,0,1.0,nvwpmbin,wpmbins);
		
		sprintf(name,"H2D_DeepAK8_tbMass_AK4btag%i_AK8CHS_highht",ib+1);
		hist_2D_DeepAK8_tbmass_3[ib] = new TH2D(name,name,20,0,1.0,nvwpmbin,wpmbins);
		hist_2D_DeepAK8_tbmass_3[ib] = new TH2D(name,name,20,0,1.0,nvwpmbin,wpmbins);
		
		sprintf(name,"H2D_MD_DeepAK8_tbMass_AK4btag%i_AK8CHS_lowbpt",ib+1);
		hist_2D_MDDeepAK8_tbmass_1[ib] = new TH2D(name,name,20,0,1.0,nvwpmbin,wpmbins);
		hist_2D_MDDeepAK8_tbmass_1[ib] = new TH2D(name,name,20,0,1.0,nvwpmbin,wpmbins);
		
		sprintf(name,"H2D_MD_DeepAK8_tbMass_AK4btag%i_AK8CHS_nobmasscut",ib+1);
		hist_2D_MDDeepAK8_tbmass_2[ib] = new TH2D(name,name,20,0,1.0,nvwpmbin,wpmbins);
		hist_2D_MDDeepAK8_tbmass_2[ib] = new TH2D(name,name,20,0,1.0,nvwpmbin,wpmbins);
		
		sprintf(name,"H2D_MD_DeepAK8_tbMass_AK4btag%i_AK8CHS",ib+1);
		hist_2D_MDDeepAK8_tbmass[ib] = new TH2D(name,name,20,0,1.0,nvwpmbin,wpmbins);
		hist_2D_MDDeepAK8_tbmass[ib] = new TH2D(name,name,20,0,1.0,nvwpmbin,wpmbins);
		
		sprintf(name,"H2D_MD_DeepAK8_tbMass_AK4btag%i_AK8CHS_highht",ib+1);
		hist_2D_MDDeepAK8_tbmass_3[ib] = new TH2D(name,name,20,0,1.0,nvwpmbin,wpmbins);
		hist_2D_MDDeepAK8_tbmass_3[ib] = new TH2D(name,name,20,0,1.0,nvwpmbin,wpmbins);
		
			
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
	    
	    hist_2D_sdmass_deeptopscore_md_AK8_DeepAK8[ib][isb] = new TH2D(name,title,nsdmassbins,sdmassbins,50,0,1.0);
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
		
		hist_2D_pt_btag_AK4_1[isb][it] = new TH2D(name,title,noptbins,ptbins,20,0,1.0);
		hist_2D_pt_btag_AK4_1[isb][it]->Sumw2();
		
		sprintf(name,"H2D_AK4_pt_btag_TopTag%i_Subjetbtag%i_AK8CHS_nobmasscut",it+1,isb+1);
		
		hist_2D_pt_btag_AK4_2[isb][it] = new TH2D(name,title,noptbins,ptbins,20,0,1.0);
		hist_2D_pt_btag_AK4_2[isb][it]->Sumw2();
		
		sprintf(name,"H2D_AK4_pt_btag_TopTag%i_Subjetbtag%i_AK8CHS_highht",it+1,isb+1);
		
		hist_2D_pt_btag_AK4_3[isb][it] = new TH2D(name,title,noptbins,ptbins,20,0,1.0);
		hist_2D_pt_btag_AK4_3[isb][it]->Sumw2();
		
		sprintf(name,"H2D_AK4_pt_btag_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS",it+1,isb+1);
		hist_2D_pt_btag_AK4_DeepAK8[isb][it] = new TH2D(name,title,noptbins,ptbins,20,0,1.0);
		hist_2D_pt_btag_AK4_DeepAK8[isb][it]->Sumw2();
		
		sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS",it+1,isb+1);
		hist_2D_pt_btag_AK4_md_DeepAK8[isb][it] = new TH2D(name,title,noptbins,ptbins,20,0,1.0);
		hist_2D_pt_btag_AK4_md_DeepAK8[isb][it]->Sumw2();
		
		for(int ibunc=0; ibunc<nbtagmax; ibunc++){
			sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_BTag%i_up",it+1,isb+1,ibunc);
			hist_2D_pt_btag_AK4_md_DeepAK8_Btag_up[isb][it][ibunc] = new TH2D(name,title,noptbins,ptbins,20,0,1.0);
			hist_2D_pt_btag_AK4_md_DeepAK8_Btag_up[isb][it][ibunc]->Sumw2();
	   
			sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_BTag%i_dn",it+1,isb+1,ibunc);
			hist_2D_pt_btag_AK4_md_DeepAK8_Btag_dn[isb][it][ibunc] = new TH2D(name,title,noptbins,ptbins,20,0,1.0);
			hist_2D_pt_btag_AK4_md_DeepAK8_Btag_dn[isb][it][ibunc]->Sumw2();
		}
		
		sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_highht",it+1,isb+1);
		hist_2D_pt_btag_AK4_md_DeepAK8_3[isb][it] = new TH2D(name,title,noptbins,ptbins,20,0,1.0);
		hist_2D_pt_btag_AK4_md_DeepAK8_3[isb][it]->Sumw2();
		
		sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_ttbar",it+1,isb+1);
		hist_2D_pt_btag_AK4_md_DeepAK8_tt[isb][it] = new TH2D(name,title,noptbins,ptbins,20,0,1.0);
		hist_2D_pt_btag_AK4_md_DeepAK8_tt[isb][it]->Sumw2();
		
		sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_ttbar_highht",it+1,isb+1);
		hist_2D_pt_btag_AK4_md_DeepAK8_tt_3[isb][it] = new TH2D(name,title,noptbins,ptbins,20,0,1.0);
		hist_2D_pt_btag_AK4_md_DeepAK8_tt_3[isb][it]->Sumw2();
		
		for(int ieta=0; ieta<netabins; ieta++){
			
			sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_bEtaBin%i",it+1,isb+1,ieta+1);
			hist_2D_pt_btag_AK4_md_DeepAK8_beta[isb][it][ieta] = new TH2D(name,title,noptbins,ptbins,20,0,1.0);
			hist_2D_pt_btag_AK4_md_DeepAK8_beta[isb][it][ieta]->Sumw2();
			
			sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_bEtaBin%i_2",it+1,isb+1,ieta+1);
			hist_2D_pt_btag_AK4_md_DeepAK8_beta2[isb][it][ieta] = new TH2D(name,title,noptbins,ptbins,20,0,1.0);
			hist_2D_pt_btag_AK4_md_DeepAK8_beta2[isb][it][ieta]->Sumw2();
			
			sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_bEtaBin%i_ttbar",it+1,isb+1,ieta+1);
			hist_2D_pt_btag_AK4_md_DeepAK8_beta_tt[isb][it][ieta] = new TH2D(name,title,noptbins,ptbins,20,0,1.0);
			hist_2D_pt_btag_AK4_md_DeepAK8_beta_tt[isb][it][ieta]->Sumw2();
			
			sprintf(name,"H2D_mtbt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_bEtaBin%i",it+1,isb+1,ieta+1);
			hist_2D_mtb_btag_AK4_md_DeepAK8_beta[isb][it][ieta] = new TH2D(name,title,nvwpmbin,wpmbins,20,0,1.0);
			hist_2D_mtb_btag_AK4_md_DeepAK8_beta[isb][it][ieta]->Sumw2();
			
			for(int ibunc=0; ibunc<nbtagmax; ibunc++){
				sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_bEtaBin%i_BTag%i_up",it+1,isb+1,ieta+1,ibunc);
				hist_2D_pt_btag_AK4_md_DeepAK8_beta_Btag_up[isb][it][ieta][ibunc] = new TH2D(name,title,noptbins,ptbins,20,0,1.0);
				hist_2D_pt_btag_AK4_md_DeepAK8_beta_Btag_up[isb][it][ieta][ibunc]->Sumw2();
	   
				sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_bEtaBin%i_BTag%i_dn",it+1,isb+1,ieta+1,ibunc);
				hist_2D_pt_btag_AK4_md_DeepAK8_beta_Btag_dn[isb][it][ieta][ibunc] = new TH2D(name,title,noptbins,ptbins,20,0,1.0);
				hist_2D_pt_btag_AK4_md_DeepAK8_beta_Btag_dn[isb][it][ieta][ibunc]->Sumw2();
					}
					
			for(int ijes=0; ijes<njesmax; ijes++){
				sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_bEtaBin%i_JES%i_up",it+1,isb+1,ieta+1,ijes+1);
				hist_2D_pt_btag_AK4_md_DeepAK8_beta_JES_up[isb][it][ieta][ijes] = new TH2D(name,title,noptbins,ptbins,20,0,1.0);
				hist_2D_pt_btag_AK4_md_DeepAK8_beta_JES_up[isb][it][ieta][ijes]->Sumw2();
	   
				sprintf(name,"H2D_AK4_pt_btag_MD_DeepAK8_TopTag%i_Subjetbtag%i_AK8CHS_bEtaBin%i_JES%i_dn",it+1,isb+1,ieta+1,ijes+1);
				hist_2D_pt_btag_AK4_md_DeepAK8_beta_JES_dn[isb][it][ieta][ijes] = new TH2D(name,title,noptbins,ptbins,20,0,1.0);
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


// Fit function
/*
err1 = new TMatrixD(5,1);
err = new TMatrixD(1,1);
for(int ieta=0; ieta<netabins; ieta++){
	cov_bpfrat_beta[ieta] = new TMatrixD(5,5);
	cov_bpfrat_val_beta[ieta] = new TMatrixD(5,5);
}
*/

err1.ResizeTo(5,1);
err.ResizeTo(1,1);
for(int ieta=0; ieta<netabins; ieta++){
	cov_bpfrat_beta[ieta].ResizeTo(5,5);
	cov_bpfrat_val_beta[ieta].ResizeTo(5,5);
}
rvar.ResizeTo(1,5);//new TMatrixD(1,5);
cvar.ResizeTo(5,1) ;//new TMatrixD(5,1);

err2.ResizeTo(3,1);
err_pol2.ResizeTo(1,1);
for(int ieta=0; ieta<netabins; ieta++){
	cov_bpfrat_pol2[ieta].ResizeTo(3,3);
	cov_bpfrat_val_pol2[ieta].ResizeTo(3,3);
}
rvar2.ResizeTo(1,3);//new TMatrixD(1,5);
cvar2.ResizeTo(3,1) ;//new TMatrixD(5,1);

// value assignement //

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
	
  for(int ij=0; ij<2; ij++){
	 
	 #ifdef B_Tight 
	  for(int imsd=0; imsd<nsdmassbins; imsd++){
		cor_b_sdmass_mc[ij][imsd] = cor_b_tight_sdmass_mc_16[ij][imsd];
		cor_b_sdmass_data[ij][imsd] = cor_b_tight_sdmass_data_16[ij][imsd];
		}
	 #else 
	 for(int imsd=0; imsd<nsdmassbins; imsd++){
		cor_b_sdmass_mc[ij][imsd] = cor_b_sdmass_mc_16[ij][imsd];
		cor_b_sdmass_data[ij][imsd] = cor_b_sdmass_data_16[ij][imsd];
		}
	 #endif	
	 }	

  for(int ieta=0; ieta<netabins; ieta++){
	  
	 for(int ivar=0; ivar<5; ivar++){
		for(int jvar=0; jvar<5; jvar++){
			
			COV_sig_data[ieta][ivar][jvar] = COV_sig_data_16[ieta][ivar][jvar] ;
			COV_sig_mc[ieta][ivar][jvar] = COV_sig_mc_16[ieta][ivar][jvar] ;
			COV_val_data[ieta][ivar][jvar] = COV_val_data_16[ieta][ivar][jvar] ;
			COV_val_mc[ieta][ivar][jvar] = COV_val_mc_16[ieta][ivar][jvar] ;
			
			}
		}
		
	 for(int ivar=0; ivar<3; ivar++){
		for(int jvar=0; jvar<3; jvar++){	
		
			#ifdef B_Tight
			
			COV_sig_data_pol2[ieta][ivar][jvar] = COV_sig_data_pol2_16_btight[ieta][ivar][jvar] ;
			COV_sig_mc_pol2[ieta][ivar][jvar] = COV_sig_mc_pol2_16_btight[ieta][ivar][jvar] ;
			COV_val_data_pol2[ieta][ivar][jvar] = COV_val_data_pol2_16_btight[ieta][ivar][jvar] ;
			COV_val_mc_pol2[ieta][ivar][jvar] = COV_val_mc_pol2_16_btight[ieta][ivar][jvar] ;
			
			#else
			
			COV_sig_data_pol2[ieta][ivar][jvar] = COV_sig_data_pol2_16[ieta][ivar][jvar] ;
			COV_sig_mc_pol2[ieta][ivar][jvar] = COV_sig_mc_pol2_16[ieta][ivar][jvar] ;
			COV_val_data_pol2[ieta][ivar][jvar] = COV_val_data_pol2_16[ieta][ivar][jvar] ;
			COV_val_mc_pol2[ieta][ivar][jvar] = COV_val_mc_pol2_16[ieta][ivar][jvar] ;
			
			#endif
			
			}
		}	
		
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
	
   for(int ij=0; ij<2; ij++){
	 for(int imsd=0; imsd<nsdmassbins; imsd++){
		cor_b_sdmass_mc[ij][imsd] = cor_b_sdmass_mc_17[ij][imsd];
		cor_b_sdmass_data[ij][imsd] = cor_b_sdmass_data_17[ij][imsd];
		}
	 }		
   
   for(int ieta=0; ieta<netabins; ieta++){
	   
	 for(int ivar=0; ivar<5; ivar++){
		for(int jvar=0; jvar<5; jvar++){
			
			COV_sig_data[ieta][ivar][jvar] = COV_sig_data_17[ieta][ivar][jvar] ;
			COV_sig_mc[ieta][ivar][jvar] = COV_sig_mc_17[ieta][ivar][jvar] ;
			COV_val_data[ieta][ivar][jvar] = COV_val_data_17[ieta][ivar][jvar] ;
			COV_val_mc[ieta][ivar][jvar] = COV_val_mc_17[ieta][ivar][jvar] ;
			
			}
		}
		
	 for(int ivar=0; ivar<3; ivar++){
		for(int jvar=0; jvar<3; jvar++){	
		
			COV_sig_data_pol2[ieta][ivar][jvar] = COV_sig_data_pol2_17[ieta][ivar][jvar] ;
			COV_sig_mc_pol2[ieta][ivar][jvar] = COV_sig_mc_pol2_17[ieta][ivar][jvar] ;
			COV_val_data_pol2[ieta][ivar][jvar] = COV_val_data_pol2_17[ieta][ivar][jvar] ;
			COV_val_mc_pol2[ieta][ivar][jvar] = COV_val_mc_pol2_17[ieta][ivar][jvar] ;
			
			}
		}
			
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
   
   for(int ij=0; ij<2; ij++){
	 for(int imsd=0; imsd<nsdmassbins; imsd++){
		cor_b_sdmass_mc[ij][imsd] = cor_b_sdmass_mc_18[ij][imsd];
		cor_b_sdmass_data[ij][imsd] = cor_b_sdmass_data_18[ij][imsd];
		}
	 }	 
   
   for(int ieta=0; ieta<netabins; ieta++){
	   
	 for(int ivar=0; ivar<5; ivar++){
		for(int jvar=0; jvar<5; jvar++){
			
			COV_sig_data[ieta][ivar][jvar] = COV_sig_data_18[ieta][ivar][jvar] ;
			COV_sig_mc[ieta][ivar][jvar] = COV_sig_mc_18[ieta][ivar][jvar] ;
			COV_val_data[ieta][ivar][jvar] = COV_val_data_18[ieta][ivar][jvar] ;
			COV_val_mc[ieta][ivar][jvar] = COV_val_mc_18[ieta][ivar][jvar] ;
			
			}
		}
		
	for(int ivar=0; ivar<3; ivar++){
		for(int jvar=0; jvar<3; jvar++){	
		
			COV_sig_data_pol2[ieta][ivar][jvar] = COV_sig_data_pol2_18[ieta][ivar][jvar] ;
			COV_sig_mc_pol2[ieta][ivar][jvar] = COV_sig_mc_pol2_18[ieta][ivar][jvar] ;
			COV_val_data_pol2[ieta][ivar][jvar] = COV_val_data_pol2_18[ieta][ivar][jvar] ;
			COV_val_mc_pol2[ieta][ivar][jvar] = COV_val_mc_pol2_18[ieta][ivar][jvar] ;
			
			}
		}
		
	 }
   
   for(int imsd=0; imsd<ntopsdmassbins; imsd++){
		for(int ht=0; ht<nohtbins; ht++){
			if(isQCD){
				trig_weight_factors[imsd][ht] = trig_weight_factors_2018_QCD[imsd][ht];
			}
			if(isTT){
				trig_weight_factors[imsd][ht] = trig_weight_factors_2018_TT[imsd][ht];
			}
			if(isST){
				trig_weight_factors[imsd][ht] = trig_weight_factors_2018_ST[imsd][ht];
			}
		}
	}
   
#endif


if(!isMC || (isMC && (isTTBar || isST))){

bpfrat = new TF1("BTagger_PFRate_SR",MikkoFunc1,500,2000,4);
bpfrat->SetParameters(-1.85166e+01,-4.56084e-01,5.20671e+01,4.91319e-01);

#ifdef Anal_2017

for(int ieta=0; ieta<netabins; ieta++){
	
	sprintf(name,"BTagger_PFRate_SR_Eta%i",ieta+1);
	bpfrat_beta[ieta]= new TF1(name,BiFun,1000,5000,5);
  	
  	if(ieta==0) { bpfrat_beta[ieta]->SetParameters(2.32358e+03,4.16906e-02,-5.86470e-06,4.04857e-09,-2.52014e-09); }
	if(ieta==1) { bpfrat_beta[ieta]->SetParameters(2.33773e+03,3.47750e-02,-3.01114e-06,4.48931e-09,6.58976e-10); }
	if(ieta==2) { bpfrat_beta[ieta]->SetParameters(2.40356e+03,2.69763e-02,-4.83897e-06,-2.48385e-09,2.76093e-09); }

	for(int ipar=0; ipar<5; ipar++){
		sprintf(name,"dFun_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_bdF[ipar][ieta]= new TF1(name,bdF0,1000,5000,5);}
		if(ipar==1){fun_bdF[ipar][ieta]= new TF1(name,bdF1,1000,5000,5);}
		if(ipar==2){fun_bdF[ipar][ieta]= new TF1(name,bdF2,1000,5000,5);}
		if(ipar==3){fun_bdF[ipar][ieta]= new TF1(name,bdF3,1000,5000,5);}
		if(ipar==4){fun_bdF[ipar][ieta]= new TF1(name,bdF4,1000,5000,5);}
	}
	
	for(int ipar=0; ipar<5; ipar++){
	
		if(ieta==0){ fun_bdF[ipar][ieta]->SetParameters(2.32358e+03,4.16906e-02,-5.86470e-06,4.04857e-09,-2.52014e-09); }
		if(ieta==1){ fun_bdF[ipar][ieta]->SetParameters(2.33773e+03,3.47750e-02,-3.01114e-06,4.48931e-09,6.58976e-10); }
		if(ieta==2){ fun_bdF[ipar][ieta]->SetParameters(2.40356e+03,2.69763e-02,-4.83897e-06,-2.48385e-09,2.76093e-09); }
	}
	
	for(int ipar=0; ipar<5; ipar++){
		for(int jpar=0; jpar<5; jpar++){
			cov_bpfrat_beta[ieta](ipar,jpar) = COV_sig_data[ieta][ipar][jpar];
		}
	}
	
	sprintf(name,"BTagger_PFRate_Pol2_SR_Eta%i",ieta+1);
	bpfrat_pol2[ieta]= new TF1(name,Pol2,1000,5000,3);
	
	if(ieta==0) { bpfrat_pol2[ieta]->SetParameters(6.74291e-02,-1.30191e-05,7.34549e-10); }
	if(ieta==1) { bpfrat_pol2[ieta]->SetParameters(5.93216e-02,-1.58721e-05,2.27318e-09); }
	if(ieta==2) { bpfrat_pol2[ieta]->SetParameters(3.56435e-02,-5.45445e-06,7.73411e-10); }
  	
	
	for(int ipar=0; ipar<3; ipar++){
		sprintf(name,"dFun_pol2_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_bdF_pol2[ipar][ieta]= new TF1(name,bdF0_pol2,1000,5000,3);}
		if(ipar==1){fun_bdF_pol2[ipar][ieta]= new TF1(name,bdF1_pol2,1000,5000,3);}
		if(ipar==2){fun_bdF_pol2[ipar][ieta]= new TF1(name,bdF2_pol2,1000,5000,3);}
	}
	
	for(int ipar=0; ipar<3; ipar++){
	
		if(ieta==0){ fun_bdF_pol2[ipar][ieta]->SetParameters(6.74291e-02,-1.30191e-05,7.34549e-10); }
		if(ieta==1){ fun_bdF_pol2[ipar][ieta]->SetParameters(5.93216e-02,-1.58721e-05,2.27318e-09); }
		if(ieta==2){ fun_bdF_pol2[ipar][ieta]->SetParameters(3.56435e-02,-5.45445e-06,7.73411e-10); }
	}
	
	for(int ipar=0; ipar<3; ipar++){
		for(int jpar=0; jpar<3; jpar++){
			cov_bpfrat_pol2[ieta](ipar,jpar) = COV_sig_data_pol2[ieta][ipar][jpar];
		}
	}
   
}

bpfrat_val = new TF1("BTagger_PFRate_VR",MikkoFunc1,500,2000,4);
bpfrat_val->SetParameters(-6.12117e+00,-4.51636e-01,1.71053e+01,2.05272e-01);

for(int ieta=0; ieta<netabins; ieta++){
	
	sprintf(name,"BTagger_PFRate_VR_Eta%i",ieta+1);
	bpfrat_val_beta[ieta]= new TF1(name,BiFun,1000,5000,5);
	
	if(ieta==0) { bpfrat_val_beta[ieta]->SetParameters(2.37490e+03,1.64886e-02,2.05504e-06,1.00313e-09,-2.40608e-09); }
	if(ieta==1) { bpfrat_val_beta[ieta]->SetParameters(1.60759e+03,1.45295e-02,7.95157e-07,-3.10981e-10,-7.35494e-11); }
	if(ieta==2) { bpfrat_val_beta[ieta]->SetParameters(1.75466e+03,1.48627e-02,1.55269e-06,3.36807e-09,-4.34977e-10); }
	
	for(int ipar=0; ipar<5; ipar++){
		sprintf(name,"dFun_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_val_bdF[ipar][ieta]= new TF1(name,bdF0,1000,5000,5);}
		if(ipar==1){fun_val_bdF[ipar][ieta]= new TF1(name,bdF1,1000,5000,5);}
		if(ipar==2){fun_val_bdF[ipar][ieta]= new TF1(name,bdF2,1000,5000,5);}
		if(ipar==3){fun_val_bdF[ipar][ieta]= new TF1(name,bdF3,1000,5000,5);}
		if(ipar==4){fun_val_bdF[ipar][ieta]= new TF1(name,bdF4,1000,5000,5);}
	}
	
	for(int ipar=0; ipar<5; ipar++){
		
		if(ieta==0){ fun_val_bdF[ipar][ieta]->SetParameters(2.37490e+03,1.64886e-02,2.05504e-06,1.00313e-09,-2.40608e-09); }
		if(ieta==1){ fun_val_bdF[ipar][ieta]->SetParameters(1.60759e+03,1.45295e-02,7.95157e-07,-3.10981e-10,-7.35494e-11); }
		if(ieta==2){ fun_val_bdF[ipar][ieta]->SetParameters(1.75466e+03,1.48627e-02,1.55269e-06,3.36807e-09,-4.34977e-10); }
	}
	
	for(int ipar=0; ipar<5; ipar++){
		for(int jpar=0; jpar<5; jpar++){
			cov_bpfrat_val_beta[ieta](ipar,jpar) = COV_val_data[ieta][ipar][jpar];
		}
	}
	
	sprintf(name,"BTagger_PFRate_Pol2_VR_Eta%i",ieta+1);
	bpfrat_val_pol2[ieta]= new TF1(name,Pol2,1000,5000,3);
	
	if(ieta==0) { bpfrat_val_pol2[ieta]->SetParameters(1.37009e-02,1.76742e-06,-3.39836e-10); }
	if(ieta==1) { bpfrat_val_pol2[ieta]->SetParameters(1.29249e-02,1.14622e-06,-9.56192e-11); }
	if(ieta==2) { bpfrat_val_pol2[ieta]->SetParameters(1.35903e-02,1.00618e-06,-8.49165e-11); }
  	
	
	for(int ipar=0; ipar<3; ipar++){
		sprintf(name,"dFun_val_pol2_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_val_bdF_pol2[ipar][ieta]= new TF1(name,bdF0_pol2,1000,5000,3);}
		if(ipar==1){fun_val_bdF_pol2[ipar][ieta]= new TF1(name,bdF1_pol2,1000,5000,3);}
		if(ipar==2){fun_val_bdF_pol2[ipar][ieta]= new TF1(name,bdF2_pol2,1000,5000,3);}
	}
	
	for(int ipar=0; ipar<3; ipar++){
	
		if(ieta==0){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(1.37009e-02,1.76742e-06,-3.39836e-10); }
		if(ieta==1){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(1.29249e-02,1.14622e-06,-9.56192e-11); }
		if(ieta==2){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(1.35903e-02,1.00618e-06,-8.49165e-11); }
	}
	
	for(int ipar=0; ipar<3; ipar++){
		for(int jpar=0; jpar<3; jpar++){
			cov_bpfrat_val_pol2[ieta](ipar,jpar) = COV_val_data_pol2[ieta][ipar][jpar];
		}
	}
	
}

#endif

#ifdef Anal_2018

for(int ieta=0; ieta<netabins; ieta++){
	
	sprintf(name,"BTagger_PFRate_SR_Eta%i",ieta+1);
	bpfrat_beta[ieta]= new TF1(name,BiFun,1000,5000,5);
		
	if(ieta==0) { bpfrat_beta[ieta]->SetParameters(2.29350e+03,3.86563e-02,-1.07188e-05,3.15084e-09,2.21607e-09); }
	if(ieta==1) { bpfrat_beta[ieta]->SetParameters(2.33510e+03,3.35633e-02,-2.33749e-06,7.60051e-09,5.62973e-11); }
	if(ieta==2) { bpfrat_beta[ieta]->SetParameters(2.40811e+03,2.63507e-02,-4.45930e-06,-3.34820e-10,2.91961e-09); }
	
	for(int ipar=0; ipar<5; ipar++){
		sprintf(name,"dFun_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_bdF[ipar][ieta]= new TF1(name,bdF0,1000,5000,5);}
		if(ipar==1){fun_bdF[ipar][ieta]= new TF1(name,bdF1,1000,5000,5);}
		if(ipar==2){fun_bdF[ipar][ieta]= new TF1(name,bdF2,1000,5000,5);}
		if(ipar==3){fun_bdF[ipar][ieta]= new TF1(name,bdF3,1000,5000,5);}
		if(ipar==4){fun_bdF[ipar][ieta]= new TF1(name,bdF4,1000,5000,5);}
	}
	for(int ipar=0; ipar<5; ipar++){
		if(ieta==0){ fun_bdF[ipar][ieta]->SetParameters(2.29350e+03,3.86563e-02,-1.07188e-05,3.15084e-09,2.21607e-09); }
		if(ieta==1){ fun_bdF[ipar][ieta]->SetParameters(2.33510e+03,3.35633e-02,-2.33749e-06,7.60051e-09,5.62973e-11); }
		if(ieta==2){ fun_bdF[ipar][ieta]->SetParameters(2.40811e+03,2.63507e-02,-4.45930e-06,-3.34820e-10,2.91961e-09); }
	}
	
	
	for(int ipar=0; ipar<5; ipar++){
		for(int jpar=0; jpar<5; jpar++){
			cov_bpfrat_beta[ieta](ipar,jpar) = COV_sig_data[ieta][ipar][jpar];
		}
	}
   
   	
	sprintf(name,"BTagger_PFRate_Pol2_SR_Eta%i",ieta+1);
	bpfrat_pol2[ieta]= new TF1(name,Pol2,1000,5000,3);
	
	if(ieta==0) { bpfrat_pol2[ieta]->SetParameters(7.84117e-02,-2.35057e-05,2.67991e-09); }
	if(ieta==1) { bpfrat_pol2[ieta]->SetParameters(6.61623e-02,-2.14266e-05,3.16112e-09); }
	if(ieta==2) { bpfrat_pol2[ieta]->SetParameters(4.24171e-02,-1.08117e-05,1.71723e-09); }
  	
	
	for(int ipar=0; ipar<3; ipar++){
		sprintf(name,"dFun_pol2_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_bdF_pol2[ipar][ieta]= new TF1(name,bdF0_pol2,1000,5000,3);}
		if(ipar==1){fun_bdF_pol2[ipar][ieta]= new TF1(name,bdF1_pol2,1000,5000,3);}
		if(ipar==2){fun_bdF_pol2[ipar][ieta]= new TF1(name,bdF2_pol2,1000,5000,3);}
	}
	
	for(int ipar=0; ipar<3; ipar++){
	
		if(ieta==0){ fun_bdF_pol2[ipar][ieta]->SetParameters(7.84117e-02,-2.35057e-05,2.67991e-09); }
		if(ieta==1){ fun_bdF_pol2[ipar][ieta]->SetParameters(6.61623e-02,-2.14266e-05,3.16112e-09); }
		if(ieta==2){ fun_bdF_pol2[ipar][ieta]->SetParameters(4.24171e-02,-1.08117e-05,1.71723e-09); }
	}
	
	for(int ipar=0; ipar<3; ipar++){
		for(int jpar=0; jpar<3; jpar++){
			cov_bpfrat_pol2[ieta](ipar,jpar) = COV_sig_data_pol2[ieta][ipar][jpar];
		}
	}
 
   
}

bpfrat_val = new TF1("BTagger_PFRate_VR",MikkoFunc1,500,2000,4);
bpfrat_val->SetParameters(-6.12117e+00,-4.51636e-01,1.71053e+01,2.05272e-01);

for(int ieta=0; ieta<netabins; ieta++){
	sprintf(name,"BTagger_PFRate_VR_Eta%i",ieta+1);
	bpfrat_val_beta[ieta]= new TF1(name,BiFun,500,2000,5);
	if(ieta==0) { bpfrat_val_beta[ieta]->SetParameters(2.37771e+03,1.57059e-02,1.52264e-06,8.91914e-10,-2.00822e-09); }
	if(ieta==1) { bpfrat_val_beta[ieta]->SetParameters(1.60470e+03,1.38313e-02,1.54416e-07,-2.40379e-09,5.25140e-11); }
	if(ieta==2) { bpfrat_val_beta[ieta]->SetParameters(1.88921e+03,1.41138e-02,6.94243e-07,-8.17460e-10,-2.58415e-11); }
	
	for(int ipar=0; ipar<5; ipar++){
		sprintf(name,"dFun_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_val_bdF[ipar][ieta]= new TF1(name,bdF0,1000,5000,5);}
		if(ipar==1){fun_val_bdF[ipar][ieta]= new TF1(name,bdF1,1000,5000,5);}
		if(ipar==2){fun_val_bdF[ipar][ieta]= new TF1(name,bdF2,1000,5000,5);}
		if(ipar==3){fun_val_bdF[ipar][ieta]= new TF1(name,bdF3,1000,5000,5);}
		if(ipar==4){fun_val_bdF[ipar][ieta]= new TF1(name,bdF4,1000,5000,5);}
	}
	for(int ipar=0; ipar<5; ipar++){
		if(ieta==0){ fun_val_bdF[ipar][ieta]->SetParameters(2.37771e+03,1.57059e-02,1.52264e-06,8.91914e-10,-2.00822e-09); }
		if(ieta==1){ fun_val_bdF[ipar][ieta]->SetParameters(1.60470e+03,1.38313e-02,1.54416e-07,-2.40379e-09,5.25140e-11); }
		if(ieta==2){ fun_val_bdF[ipar][ieta]->SetParameters(1.88921e+03,1.41138e-02,6.94243e-07,-8.17460e-10,-2.58415e-11); }
	}
	
	for(int ipar=0; ipar<5; ipar++){
		for(int jpar=0; jpar<5; jpar++){
			cov_bpfrat_val_beta[ieta](ipar,jpar) = COV_val_data[ieta][ipar][jpar];
		}
	}
   
	sprintf(name,"BTagger_PFRate_Pol2_VR_Eta%i",ieta+1);
	bpfrat_val_pol2[ieta]= new TF1(name,Pol2,1000,5000,3);
	
	if(ieta==0) { bpfrat_val_pol2[ieta]->SetParameters(1.40055e-02,1.18775e-06,-2.74695e-10); }
	if(ieta==1) { bpfrat_val_pol2[ieta]->SetParameters(1.23360e-02,1.14762e-06,-1.71305e-10); }
	if(ieta==2) { bpfrat_val_pol2[ieta]->SetParameters(1.17753e-02,1.47266e-06,-1.39420e-10); }
  	
	
	for(int ipar=0; ipar<3; ipar++){
		sprintf(name,"dFun_val_pol2_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_val_bdF_pol2[ipar][ieta]= new TF1(name,bdF0_pol2,1000,5000,3);}
		if(ipar==1){fun_val_bdF_pol2[ipar][ieta]= new TF1(name,bdF1_pol2,1000,5000,3);}
		if(ipar==2){fun_val_bdF_pol2[ipar][ieta]= new TF1(name,bdF2_pol2,1000,5000,3);}
	}
	
	for(int ipar=0; ipar<3; ipar++){
	
		if(ieta==0){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(1.40055e-02,1.18775e-06,-2.74695e-10); }
		if(ieta==1){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(1.23360e-02,1.14762e-06,-1.71305e-10); }
		if(ieta==2){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(1.17753e-02,1.47266e-06,-1.39420e-10); }
	}
	
	for(int ipar=0; ipar<3; ipar++){
		for(int jpar=0; jpar<3; jpar++){
			cov_bpfrat_val_pol2[ieta](ipar,jpar) = COV_val_data_pol2[ieta][ipar][jpar];
		}
	}
   
}

#endif

#ifdef Anal_2016

bpfrat = new TF1("BTagger_PFRate_SR",MikkoFunc1,500,2000,4);
bpfrat->SetParameters(-1.85166e+01,-4.56084e-01,5.20671e+01,4.91319e-01);

for(int ieta=0; ieta<netabins; ieta++){
	sprintf(name,"BTagger_PFRate_SR_Eta%i",ieta+1);
	bpfrat_beta[ieta]= new TF1(name,BiFun,1000,5000,5);	
	
	#ifdef B_Tight
	
	if(ieta==0) { bpfrat_beta[ieta]->SetParameters(2.29122e+03,4.77575e-02,-7.77317e-06,-6.32208e-09,-6.82567e-12); }
	if(ieta==1) { bpfrat_beta[ieta]->SetParameters(2.32450e+03,5.43731e-02,-2.00072e-06,-4.40247e-09,2.30975e-11); }
	if(ieta==2) { bpfrat_beta[ieta]->SetParameters(1.68235e+03,4.08733e-02,-1.36196e-06,-5.02670e-09,-5.07570e-10); }
	
	#else
	/*
	if(ieta==0) { bpfrat_beta[ieta]->SetParameters(1.51559e+03,9.42536e-02,2.13754e-05,-8.46763e-10,9.40950e-10); }
	if(ieta==1) { bpfrat_beta[ieta]->SetParameters(1.63928e+03,1.11094e-01,3.27249e-05,-6.38342e-10,-8.62401e-09); }
	if(ieta==2) { bpfrat_beta[ieta]->SetParameters(1.92606e+03,1.13426e-01,6.25910e-06,-3.19538e-08,-3.90406e-09); }
	*/
	if(ieta==0) { bpfrat_beta[ieta]->SetParameters(1.62883e+03,9.64815e-02,2.42160e-05,8.71399e-09,-1.46361e-09); }
	if(ieta==1) { bpfrat_beta[ieta]->SetParameters(1.75281e+03,1.14714e-01,2.97073e-05,-7.18638e-09,-7.69304e-09); }
	if(ieta==2) { bpfrat_beta[ieta]->SetParameters(1.88203e+03,1.13035e-01,5.61238e-06,-3.84968e-08,-3.41842e-09); }
	
	#endif
	
	for(int ipar=0; ipar<5; ipar++){
		sprintf(name,"dFun_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_bdF[ipar][ieta]= new TF1(name,bdF0,1000,5000,5);}
		if(ipar==1){fun_bdF[ipar][ieta]= new TF1(name,bdF1,1000,5000,5);}
		if(ipar==2){fun_bdF[ipar][ieta]= new TF1(name,bdF2,1000,5000,5);}
		if(ipar==3){fun_bdF[ipar][ieta]= new TF1(name,bdF3,1000,5000,5);}
		if(ipar==4){fun_bdF[ipar][ieta]= new TF1(name,bdF4,1000,5000,5);}
	}
	for(int ipar=0; ipar<5; ipar++){
		
		#ifdef B_Tight
		
		if(ieta==0){ fun_bdF[ipar][ieta]->SetParameters(2.29122e+03,4.77575e-02,-7.77317e-06,-6.32208e-09,-6.82567e-12); }
		if(ieta==1){ fun_bdF[ipar][ieta]->SetParameters(2.32450e+03,5.43731e-02,-2.00072e-06,-4.40247e-09,2.30975e-11); }
		if(ieta==2){ fun_bdF[ipar][ieta]->SetParameters(1.68235e+03,4.08733e-02,-1.36196e-06,-5.02670e-09,-5.07570e-10); }
		
		#else
		/*
		if(ieta==0){ fun_bdF[ipar][ieta]->SetParameters(1.51559e+03,9.42536e-02,2.13754e-05,-8.46763e-10,9.40950e-10); }
		if(ieta==1){ fun_bdF[ipar][ieta]->SetParameters(1.63928e+03,1.11094e-01,3.27249e-05,-6.38342e-10,-8.62401e-09); }
		if(ieta==2){ fun_bdF[ipar][ieta]->SetParameters(1.92606e+03,1.13426e-01,6.25910e-06,-3.19538e-08,-3.90406e-09); }
		*/
		if(ieta==0) { fun_bdF[ipar][ieta]->SetParameters(1.62883e+03,9.64815e-02,2.42160e-05,8.71399e-09,-1.46361e-09); }
		if(ieta==1) { fun_bdF[ipar][ieta]->SetParameters(1.75281e+03,1.14714e-01,2.97073e-05,-7.18638e-09,-7.69304e-09); }
		if(ieta==2) { fun_bdF[ipar][ieta]->SetParameters(1.88203e+03,1.13035e-01,5.61238e-06,-3.84968e-08,-3.41842e-09); }
		
		#endif
	}
	
	
	for(int ipar=0; ipar<5; ipar++){
		for(int jpar=0; jpar<5; jpar++){
			cov_bpfrat_beta[ieta](ipar,jpar) = COV_sig_data[ieta][ipar][jpar];
		}
	}
   
   	
	sprintf(name,"BTagger_PFRate_Pol2_SR_Eta%i",ieta+1);
	bpfrat_pol2[ieta]= new TF1(name,Pol2,1000,5000,3);
	
  	#ifdef B_Tight
	if(ieta==0) { bpfrat_pol2[ieta]->SetParameters(4.13346e-02,1.04929e-05,-3.26939e-09); }
	if(ieta==1) { bpfrat_pol2[ieta]->SetParameters(4.27580e-02,9.60189e-06,-1.94551e-09); }
	if(ieta==2) { bpfrat_pol2[ieta]->SetParameters(4.03561e-02,1.31197e-06,-6.63512e-10); }
	#else
	/*
	if(ieta==0) { bpfrat_pol2[ieta]->SetParameters(6.33366e-02,1.91396e-05,8.10967e-10); }
	if(ieta==1) { bpfrat_pol2[ieta]->SetParameters(3.86670e-02,5.74462e-05,-7.96772e-09); }
	if(ieta==2) { bpfrat_pol2[ieta]->SetParameters(5.23110e-02,4.58462e-05,-7.88366e-09); }
	*/
	if(ieta==0) { bpfrat_pol2[ieta]->SetParameters(6.02485e-02,2.26926e-05,-1.57336e-10); }
	if(ieta==1) { bpfrat_pol2[ieta]->SetParameters(3.94634e-02,5.63117e-05,-7.62595e-09); }
	if(ieta==2) { bpfrat_pol2[ieta]->SetParameters(5.33156e-02,4.48943e-05,-7.71915e-09); }
	#endif
	
	for(int ipar=0; ipar<3; ipar++){
		sprintf(name,"dFun_pol2_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_bdF_pol2[ipar][ieta]= new TF1(name,bdF0_pol2,1000,5000,3);}
		if(ipar==1){fun_bdF_pol2[ipar][ieta]= new TF1(name,bdF1_pol2,1000,5000,3);}
		if(ipar==2){fun_bdF_pol2[ipar][ieta]= new TF1(name,bdF2_pol2,1000,5000,3);}
	}
	
	for(int ipar=0; ipar<3; ipar++){
		
		#ifdef B_Tight
		if(ieta==0){ fun_bdF_pol2[ipar][ieta]->SetParameters(4.13346e-02,1.04929e-05,-3.26939e-09); }
		if(ieta==1){ fun_bdF_pol2[ipar][ieta]->SetParameters(4.27580e-02,9.60189e-06,-1.94551e-09); }
		if(ieta==2){ fun_bdF_pol2[ipar][ieta]->SetParameters(4.03561e-02,1.31197e-06,-6.63512e-10); }
		#else
		/*
		if(ieta==0){ fun_bdF_pol2[ipar][ieta]->SetParameters(6.33366e-02,1.91396e-05,8.10967e-10); }
		if(ieta==1){ fun_bdF_pol2[ipar][ieta]->SetParameters(3.86670e-02,5.74462e-05,-7.96772e-09); }
		if(ieta==2){ fun_bdF_pol2[ipar][ieta]->SetParameters(5.23110e-02,4.58462e-05,-7.88366e-09); }
		*/
		if(ieta==0){ fun_bdF_pol2[ipar][ieta]->SetParameters(6.02485e-02,2.26926e-05,-1.57336e-10); }
		if(ieta==1){ fun_bdF_pol2[ipar][ieta]->SetParameters(3.94634e-02,5.63117e-05,-7.62595e-09); }
		if(ieta==2){ fun_bdF_pol2[ipar][ieta]->SetParameters(5.33156e-02,4.48943e-05,-7.71915e-09); }
		#endif
		
	}
	
	for(int ipar=0; ipar<3; ipar++){
		for(int jpar=0; jpar<3; jpar++){
			cov_bpfrat_pol2[ieta](ipar,jpar) = COV_sig_data_pol2[ieta][ipar][jpar];
		}
	}
 
}

bpfrat_val = new TF1("BTagger_PFRate_VR",MikkoFunc1,500,2000,4);
bpfrat_val->SetParameters(-6.12117e+00,-4.51636e-01,1.71053e+01,2.05272e-01);

for(int ieta=0; ieta<netabins; ieta++){
	sprintf(name,"BTagger_PFRate_VR_Eta%i",ieta+1);
	bpfrat_val_beta[ieta]= new TF1(name,BiFun,500,2000,5);
	
	#ifdef B_Tight
	
	if(ieta==0) { bpfrat_val_beta[ieta]->SetParameters(2.38813e+03,2.00962e-02,4.44519e-06,1.57842e-10,-4.67227e-09); }
	if(ieta==1) { bpfrat_val_beta[ieta]->SetParameters(1.60168e+03,2.24037e-02,5.85461e-06,-9.20290e-09,-3.51063e-10); }
	if(ieta==2) { bpfrat_val_beta[ieta]->SetParameters(1.65633e+03,2.31584e-02,1.42368e-06,-8.58618e-09,-5.77951e-10); }
	
	#else
	/*
	if(ieta==0) { bpfrat_val_beta[ieta]->SetParameters(1.46952e+03,3.40098e-02,1.89387e-05,2.13629e-08,1.52039e-09); }
	if(ieta==1) { bpfrat_val_beta[ieta]->SetParameters(1.56631e+03,5.18521e-02,2.29592e-05,-1.47997e-08,4.33517e-11); }
	if(ieta==2) { bpfrat_val_beta[ieta]->SetParameters(1.74334e+03,6.85619e-02,7.39425e-06,-2.81162e-08,-9.57036e-10); }
	*/
	if(ieta==0) { bpfrat_val_beta[ieta]->SetParameters(1.38373e+03,3.25541e-02,1.85890e-05,3.45274e-08,1.50838e-09); }
	if(ieta==1) { bpfrat_val_beta[ieta]->SetParameters(1.48609e+03,4.99726e-02,2.31002e-05,-2.68049e-08,1.71155e-10); }
	if(ieta==2) { bpfrat_val_beta[ieta]->SetParameters(1.63503e+03,6.72515e-02,8.78198e-06,-3.98430e-08,-1.37735e-09); }
	
	#endif
	
	for(int ipar=0; ipar<5; ipar++){
		sprintf(name,"dFun_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_val_bdF[ipar][ieta]= new TF1(name,bdF0,1000,5000,5);}
		if(ipar==1){fun_val_bdF[ipar][ieta]= new TF1(name,bdF1,1000,5000,5);}
		if(ipar==2){fun_val_bdF[ipar][ieta]= new TF1(name,bdF2,1000,5000,5);}
		if(ipar==3){fun_val_bdF[ipar][ieta]= new TF1(name,bdF3,1000,5000,5);}
		if(ipar==4){fun_val_bdF[ipar][ieta]= new TF1(name,bdF4,1000,5000,5);}
	}
	for(int ipar=0; ipar<5; ipar++){
		
		#ifdef B_Tight
		
		if(ieta==0){ fun_val_bdF[ipar][ieta]->SetParameters(2.38813e+03,2.00962e-02,4.44519e-06,1.57842e-10,-4.67227e-09); }
		if(ieta==1){ fun_val_bdF[ipar][ieta]->SetParameters(1.60168e+03,2.24037e-02,5.85461e-06,-9.20290e-09,-3.51063e-10); }
		if(ieta==2){ fun_val_bdF[ipar][ieta]->SetParameters(1.65633e+03,2.31584e-02,1.42368e-06,-8.58618e-09,-5.77951e-10); }
		
		#else
		/*
		if(ieta==0){ fun_val_bdF[ipar][ieta]->SetParameters(1.46952e+03,3.40098e-02,1.89387e-05,2.13629e-08,1.52039e-09); }
		if(ieta==1){ fun_val_bdF[ipar][ieta]->SetParameters(1.56631e+03,5.18521e-02,2.29592e-05,-1.47997e-08,4.33517e-11); }
		if(ieta==2){ fun_val_bdF[ipar][ieta]->SetParameters(1.74334e+03,6.85619e-02,7.39425e-06,-2.81162e-08,-9.57036e-10); }
		*/
		if(ieta==0){ fun_val_bdF[ipar][ieta]->SetParameters(1.38373e+03,3.25541e-02,1.85890e-05,3.45274e-08,1.50838e-09); }
		if(ieta==1){ fun_val_bdF[ipar][ieta]->SetParameters(1.48609e+03,4.99726e-02,2.31002e-05,-2.68049e-08,1.71155e-10); }
		if(ieta==2){ fun_val_bdF[ipar][ieta]->SetParameters(1.63503e+03,6.72515e-02,8.78198e-06,-3.98430e-08,-1.37735e-09); }
		
		#endif
	}
	
	for(int ipar=0; ipar<5; ipar++){
		for(int jpar=0; jpar<5; jpar++){
			cov_bpfrat_val_beta[ieta](ipar,jpar) = COV_val_data[ieta][ipar][jpar];
		}
	}
	
	
	sprintf(name,"BTagger_PFRate_Pol2_VR_Eta%i",ieta+1);
	bpfrat_val_pol2[ieta]= new TF1(name,Pol2,1000,5000,3);
	
  	#ifdef B_Tight
  	if(ieta==0) { bpfrat_val_pol2[ieta]->SetParameters(4.98020e-03,1.03754e-05,-1.82367e-09); }
	if(ieta==1) { bpfrat_val_pol2[ieta]->SetParameters(6.55454e-03,1.18117e-05,-1.31684e-09); }
	if(ieta==2) { bpfrat_val_pol2[ieta]->SetParameters(1.55511e-02,6.07535e-06,-1.04245e-09); }
	#else
	/*
	if(ieta==0) { bpfrat_val_pol2[ieta]->SetParameters(1.75482e-02,6.64651e-06,3.28207e-09); }
	if(ieta==1) { bpfrat_val_pol2[ieta]->SetParameters(7.58131e-03,3.03083e-05,-1.49105e-09); }
	if(ieta==2) { bpfrat_val_pol2[ieta]->SetParameters(3.25262e-02,2.58994e-05,-3.54484e-09); }
	*/
	if(ieta==0) { bpfrat_val_pol2[ieta]->SetParameters(1.61478e-02,8.10430e-06,2.94908e-09); }
	if(ieta==1) { bpfrat_val_pol2[ieta]->SetParameters(6.66628e-03,3.10612e-05,-1.59688e-09); }
	if(ieta==2) { bpfrat_val_pol2[ieta]->SetParameters(3.18996e-02,2.64498e-05,-3.65144e-09); }
	#endif
	
	for(int ipar=0; ipar<3; ipar++){
		sprintf(name,"dFun_val_pol2_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_val_bdF_pol2[ipar][ieta]= new TF1(name,bdF0_pol2,1000,5000,3);}
		if(ipar==1){fun_val_bdF_pol2[ipar][ieta]= new TF1(name,bdF1_pol2,1000,5000,3);}
		if(ipar==2){fun_val_bdF_pol2[ipar][ieta]= new TF1(name,bdF2_pol2,1000,5000,3);}
	}
	
	#ifdef B_Tight
	for(int ipar=0; ipar<3; ipar++){
		if(ieta==0){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(4.98020e-03,1.03754e-05,-1.82367e-09); }
		if(ieta==1){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(6.55454e-03,1.18117e-05,-1.31684e-09); }
		if(ieta==2){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(1.55511e-02,6.07535e-06,-1.04245e-09); }
	}
	#else
	for(int ipar=0; ipar<3; ipar++){
		/*
		if(ieta==0){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(1.75482e-02,6.64651e-06,3.28207e-09); }
		if(ieta==1){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(7.58131e-03,3.03083e-05,-1.49105e-09); }
		if(ieta==2){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(3.25262e-02,2.58994e-05,-3.54484e-09); }
		*/
		if(ieta==0){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(1.61478e-02,8.10430e-06,2.94908e-09); }
		if(ieta==1){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(6.66628e-03,3.10612e-05,-1.59688e-09); }
		if(ieta==2){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(3.18996e-02,2.64498e-05,-3.65144e-09); } 
	}
	#endif
	
	for(int ipar=0; ipar<3; ipar++){
		for(int jpar=0; jpar<3; jpar++){
			cov_bpfrat_val_pol2[ieta](ipar,jpar) = COV_val_data_pol2[ieta][ipar][jpar];
		}
	}
   
}

#endif

tagpfrat_tau32 = new TF1("TopTagger_Tau32_PFRate",pol2,20,400,3);
tagpfrat_tau32->SetParameters(7.65157e-03,-2.05528e-04,3.84749e-06);

tagpfrat_tau32_val = new TF1("TopTagger_Tau32_PFRate_VR",pol2,20,400,3);
tagpfrat_tau32_val->SetParameters(1.70281e-02,-4.12120e-04,4.66585e-06);

tagpfrat_MDAK8 = new TF1("TopTagger_MDDeepAK8_PFRate",Parabol,50,400,4);
//tagpfrat_MDAK8->SetParameters(1.33772e+02,1.04802e-02,4.03748e-07,9.39636e-07); // pt uncorrected
tagpfrat_MDAK8->SetParameters(1.24900e+02,3.25593e-02,-9.37709e-07,1.39243e-06);// pt uncorrected

tagpfrat_MDAK8_val = new TF1("TopTagger_MDDeepAK8_PFRate_VR",Parabol,50,400,4);
//tagpfrat_MDAK8_val->SetParameters(1.30809e+02,7.44897e-03,-2.61440e-08,7.79831e-07); // pt uncorrected
tagpfrat_MDAK8_val->SetParameters(1.24645e+02,1.48181e-02,-5.01859e-07,8.25349e-07);// pt corrected

tagpfrat_MDAK8_val2 = new TF1("TopTagger_MDDeepAK8_PFRate_VR",pol2,20,400,3);
tagpfrat_MDAK8_val2->SetParameters(1.22432e-02,-9.98876e-05,3.91304e-07);

}

if(isMC && isQCD){

bpfrat = new TF1("BTagger_PFRate_SR",MikkoFunc1,400,2000,4);
bpfrat->SetParameters(-1.13445e+00,-5.30009e-01,9.42228e-01,6.32387e-02);

#ifdef Anal_2017

for(int ieta=0; ieta<netabins; ieta++){
	
	sprintf(name,"BTagger_PFRate_SR_Eta%i",ieta+1);
	bpfrat_beta[ieta]= new TF1(name,BiFun,1000,5000,5);
	
	if(ieta==0) { bpfrat_beta[ieta]->SetParameters(2.38815e+03,3.08787e-02,-7.72864e-06,1.67637e-09,4.59941e-10); }
	if(ieta==1) { bpfrat_beta[ieta]->SetParameters(2.47737e+03,2.53782e-02,-1.92963e-06,6.31278e-09,-4.67200e-10); }
	if(ieta==2) { bpfrat_beta[ieta]->SetParameters(2.52001e+03,2.31693e-02,-1.71081e-06,2.55212e-09,-1.86928e-10); }
	
	for(int ipar=0; ipar<5; ipar++){
		sprintf(name,"dFun_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_bdF[ipar][ieta]= new TF1(name,bdF0,1000,5000,5);}
		if(ipar==1){fun_bdF[ipar][ieta]= new TF1(name,bdF1,1000,5000,5);}
		if(ipar==2){fun_bdF[ipar][ieta]= new TF1(name,bdF2,1000,5000,5);}
		if(ipar==3){fun_bdF[ipar][ieta]= new TF1(name,bdF3,1000,5000,5);}
		if(ipar==4){fun_bdF[ipar][ieta]= new TF1(name,bdF4,1000,5000,5);}
	}
	
	for(int ipar=0; ipar<5; ipar++){
		if(ieta==0){ fun_bdF[ipar][ieta]->SetParameters(2.38815e+03,3.08787e-02,-7.72864e-06,1.67637e-09,4.59941e-10); }
		if(ieta==1){ fun_bdF[ipar][ieta]->SetParameters(2.47737e+03,2.53782e-02,-1.92963e-06,6.31278e-09,-4.67200e-10); }
		if(ieta==2){ fun_bdF[ipar][ieta]->SetParameters(2.52001e+03,2.31693e-02,-1.71081e-06,2.55212e-09,-1.86928e-10); } 
	}
	for(int ipar=0; ipar<5; ipar++){
		for(int jpar=0; jpar<5; jpar++){
			cov_bpfrat_beta[ieta](ipar,jpar) = COV_sig_mc[ieta][ipar][jpar];
		}
	}
	
		
	sprintf(name,"BTagger_PFRate_Pol2_SR_Eta%i",ieta+1);
	bpfrat_pol2[ieta]= new TF1(name,Pol2,1000,5000,3);
	
  	if(ieta==0) { bpfrat_pol2[ieta]->SetParameters(5.60160e-02,-1.26041e-05,8.77695e-10); }
	if(ieta==1) { bpfrat_pol2[ieta]->SetParameters(5.34254e-02,-1.65591e-05,2.10147e-09); }
	if(ieta==2) { bpfrat_pol2[ieta]->SetParameters(3.70577e-02,-7.59598e-06,8.27945e-10); }
	
	for(int ipar=0; ipar<3; ipar++){
		sprintf(name,"dFun_pol2_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_bdF_pol2[ipar][ieta]= new TF1(name,bdF0_pol2,1000,5000,3);}
		if(ipar==1){fun_bdF_pol2[ipar][ieta]= new TF1(name,bdF1_pol2,1000,5000,3);}
		if(ipar==2){fun_bdF_pol2[ipar][ieta]= new TF1(name,bdF2_pol2,1000,5000,3);}
	}
	
	for(int ipar=0; ipar<3; ipar++){
		if(ieta==0){ fun_bdF_pol2[ipar][ieta]->SetParameters(5.60160e-02,-1.26041e-05,8.77695e-10); }
		if(ieta==1){ fun_bdF_pol2[ipar][ieta]->SetParameters(5.34254e-02,-1.65591e-05,2.10147e-09); }
		if(ieta==2){ fun_bdF_pol2[ipar][ieta]->SetParameters(3.70577e-02,-7.59598e-06,8.27945e-10); }
	}
	
	for(int ipar=0; ipar<3; ipar++){
		for(int jpar=0; jpar<3; jpar++){
			cov_bpfrat_pol2[ieta](ipar,jpar) = COV_sig_mc_pol2[ieta][ipar][jpar];
		}
	}
 
}

bpfrat_val = new TF1("BTagger_PFRate_VR",MikkoFunc1,400,2000,4);
bpfrat_val->SetParameters(1.25644e+01,-6.72034e-01,-1.34758e+01,1.22587e-02);

for(int ieta=0; ieta<netabins; ieta++){
	sprintf(name,"BTagger_PFRate_VR_Eta%i",ieta+1);
	bpfrat_val_beta[ieta]= new TF1(name,BiFun,1000,5000,5);
	
	if(ieta==0) { bpfrat_val_beta[ieta]->SetParameters(1.83843e+03,1.12613e-02,1.90763e-06,1.65728e-09,-7.50511e-10); }
	if(ieta==1) { bpfrat_val_beta[ieta]->SetParameters(1.82668e+03,1.11596e-02,1.50070e-07,-3.56212e-10,3.90915e-12); }
	if(ieta==2) { bpfrat_val_beta[ieta]->SetParameters(2.41870e+03,1.32515e-02,-2.22957e-08,-3.12136e-10,2.74437e-10); }
	
	for(int ipar=0; ipar<5; ipar++){
		sprintf(name,"dFun_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_val_bdF[ipar][ieta]= new TF1(name,bdF0,1000,5000,5);}
		if(ipar==1){fun_val_bdF[ipar][ieta]= new TF1(name,bdF1,1000,5000,5);}
		if(ipar==2){fun_val_bdF[ipar][ieta]= new TF1(name,bdF2,1000,5000,5);}
		if(ipar==3){fun_val_bdF[ipar][ieta]= new TF1(name,bdF3,1000,5000,5);}
		if(ipar==4){fun_val_bdF[ipar][ieta]= new TF1(name,bdF4,1000,5000,5);}
	}
	for(int ipar=0; ipar<5; ipar++){
		if(ieta==0){ fun_val_bdF[ipar][ieta]->SetParameters(1.83843e+03,1.12613e-02,1.90763e-06,1.65728e-09,-7.50511e-10); }
		if(ieta==1){ fun_val_bdF[ipar][ieta]->SetParameters(1.82668e+03,1.11596e-02,1.50070e-07,-3.56212e-10,3.90915e-12); }
		if(ieta==2){ fun_val_bdF[ipar][ieta]->SetParameters(2.41870e+03,1.32515e-02,-2.22957e-08,-3.12136e-10,2.74437e-10); }
	}
	
	
	for(int ipar=0; ipar<5; ipar++){
		for(int jpar=0; jpar<5; jpar++){
			cov_bpfrat_val_beta[ieta](ipar,jpar) = COV_val_mc[ieta][ipar][jpar];
		}
	}
	
	
	sprintf(name,"BTagger_PFRate_Pol2_VR_Eta%i",ieta+1);
	bpfrat_val_pol2[ieta]= new TF1(name,Pol2,1000,5000,3);
	
	if(ieta==0) { bpfrat_val_pol2[ieta]->SetParameters(8.05651e-03,2.49117e-06,-3.70402e-10); }
	if(ieta==1) { bpfrat_val_pol2[ieta]->SetParameters(1.05382e-02,3.99037e-07,-3.94526e-11); }
	if(ieta==2) { bpfrat_val_pol2[ieta]->SetParameters(1.27873e-02,5.06370e-08,5.78976e-11); }
	
	for(int ipar=0; ipar<3; ipar++){
		sprintf(name,"dFun_val_pol2_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_val_bdF_pol2[ipar][ieta]= new TF1(name,bdF0_pol2,1000,5000,3);}
		if(ipar==1){fun_val_bdF_pol2[ipar][ieta]= new TF1(name,bdF1_pol2,1000,5000,3);}
		if(ipar==2){fun_val_bdF_pol2[ipar][ieta]= new TF1(name,bdF2_pol2,1000,5000,3);}
	}
	
	for(int ipar=0; ipar<3; ipar++){
		if(ieta==0){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(8.05651e-03,2.49117e-06,-3.70402e-10); }
		if(ieta==1){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(1.05382e-02,3.99037e-07,-3.94526e-11); }
		if(ieta==2){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(1.27873e-02,5.06370e-08,5.78976e-11); }
	}
	
	for(int ipar=0; ipar<3; ipar++){
		for(int jpar=0; jpar<3; jpar++){
			cov_bpfrat_val_pol2[ieta](ipar,jpar) = COV_val_mc_pol2[ieta][ipar][jpar];
		}
	}
   
}

#endif


#ifdef Anal_2018

for(int ieta=0; ieta<netabins; ieta++){
	
	sprintf(name,"BTagger_PFRate_SR_Eta%i",ieta+1);
	bpfrat_beta[ieta]= new TF1(name,BiFun,1000,5000,5);
	
	if(ieta==0) { bpfrat_beta[ieta]->SetParameters(2.38631e+03,3.01222e-02,-3.00223e-06,1.07079e-08,-2.09702e-09); }
	if(ieta==1) { bpfrat_beta[ieta]->SetParameters(2.45145e+03,2.74238e-02,-4.85018e-06,2.17370e-09,1.03522e-09); }
	if(ieta==2) { bpfrat_beta[ieta]->SetParameters(2.40609e+03,2.33899e-02,-8.84349e-06,-5.24304e-09,3.83809e-09); }
	
	for(int ipar=0; ipar<5; ipar++){
		sprintf(name,"dFun_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_bdF[ipar][ieta]= new TF1(name,bdF0,1000,5000,5);}
		if(ipar==1){fun_bdF[ipar][ieta]= new TF1(name,bdF1,1000,5000,5);}
		if(ipar==2){fun_bdF[ipar][ieta]= new TF1(name,bdF2,1000,5000,5);}
		if(ipar==3){fun_bdF[ipar][ieta]= new TF1(name,bdF3,1000,5000,5);}
		if(ipar==4){fun_bdF[ipar][ieta]= new TF1(name,bdF4,1000,5000,5);}
	}
	for(int ipar=0; ipar<5; ipar++){
		if(ieta==0){ fun_bdF[ipar][ieta]->SetParameters(2.38631e+03,3.01222e-02,-3.00223e-06,1.07079e-08,-2.09702e-09); }
		if(ieta==1){ fun_bdF[ipar][ieta]->SetParameters(2.45145e+03,2.74238e-02,-4.85018e-06,2.17370e-09,1.03522e-09); }
		if(ieta==2){ fun_bdF[ipar][ieta]->SetParameters(2.40609e+03,2.33899e-02,-8.84349e-06,-5.24304e-09,3.83809e-09); }
	}
	
	for(int ipar=0; ipar<5; ipar++){
		for(int jpar=0; jpar<5; jpar++){
			cov_bpfrat_beta[ieta](ipar,jpar) = COV_sig_mc[ieta][ipar][jpar];
		}
	}
	
		
	sprintf(name,"BTagger_PFRate_Pol2_SR_Eta%i",ieta+1);
	bpfrat_pol2[ieta]= new TF1(name,Pol2,1000,5000,3);
	
	if(ieta==0) { bpfrat_pol2[ieta]->SetParameters(6.73496e-02,-2.09486e-05,2.29190e-09); }
	if(ieta==1) { bpfrat_pol2[ieta]->SetParameters(4.97077e-02,-1.26293e-05,1.44608e-09); }
	if(ieta==2) { bpfrat_pol2[ieta]->SetParameters(3.86554e-02,-9.17069e-06,1.05423e-09); }
  	
	
	for(int ipar=0; ipar<3; ipar++){
		sprintf(name,"dFun_pol2_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_bdF_pol2[ipar][ieta]= new TF1(name,bdF0_pol2,1000,5000,3);}
		if(ipar==1){fun_bdF_pol2[ipar][ieta]= new TF1(name,bdF1_pol2,1000,5000,3);}
		if(ipar==2){fun_bdF_pol2[ipar][ieta]= new TF1(name,bdF2_pol2,1000,5000,3);}
	}
	
	for(int ipar=0; ipar<3; ipar++){
	
		if(ieta==0){ fun_bdF_pol2[ipar][ieta]->SetParameters(6.73496e-02,-2.09486e-05,2.29190e-09); }
		if(ieta==1){ fun_bdF_pol2[ipar][ieta]->SetParameters(4.97077e-02,-1.26293e-05,1.44608e-09); }
		if(ieta==2){ fun_bdF_pol2[ipar][ieta]->SetParameters(3.86554e-02,-9.17069e-06,1.05423e-09); }
	}
	
	for(int ipar=0; ipar<3; ipar++){
		for(int jpar=0; jpar<3; jpar++){
			cov_bpfrat_pol2[ieta](ipar,jpar) = COV_sig_mc_pol2[ieta][ipar][jpar];
		}
	}
 
}

bpfrat_val = new TF1("BTagger_PFRate_VR",MikkoFunc1,400,2000,4);
bpfrat_val->SetParameters(1.25644e+01,-6.72034e-01,-1.34758e+01,1.22587e-02);

for(int ieta=0; ieta<netabins; ieta++){
	
	sprintf(name,"BTagger_PFRate_VR_Eta%i",ieta+1);
	bpfrat_val_beta[ieta]= new TF1(name,BiFun,500,2000,5);
	if(ieta==0) { bpfrat_val_beta[ieta]->SetParameters(2.34737e+03,1.27890e-02,8.02288e-07,5.02738e-10,-1.18313e-09); }
	if(ieta==1) { bpfrat_val_beta[ieta]->SetParameters(1.81111e+03,1.20177e-02,9.54510e-07,6.15626e-10,1.53146e-10); }
	if(ieta==2) { bpfrat_val_beta[ieta]->SetParameters(1.79748e+03,1.39612e-02,-5.32411e-07,-2.99721e-09,4.19168e-10); }
	
	for(int ipar=0; ipar<5; ipar++){
		sprintf(name,"dFun_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_val_bdF[ipar][ieta]= new TF1(name,bdF0,1000,5000,5);}
		if(ipar==1){fun_val_bdF[ipar][ieta]= new TF1(name,bdF1,1000,5000,5);}
		if(ipar==2){fun_val_bdF[ipar][ieta]= new TF1(name,bdF2,1000,5000,5);}
		if(ipar==3){fun_val_bdF[ipar][ieta]= new TF1(name,bdF3,1000,5000,5);}
		if(ipar==4){fun_val_bdF[ipar][ieta]= new TF1(name,bdF4,1000,5000,5);}
	}
	for(int ipar=0; ipar<5; ipar++){
		if(ieta==0){ fun_val_bdF[ipar][ieta]->SetParameters(2.34737e+03,1.27890e-02,8.02288e-07,5.02738e-10,-1.18313e-09); }
		if(ieta==1){ fun_val_bdF[ipar][ieta]->SetParameters(1.81111e+03,1.20177e-02,9.54510e-07,6.15626e-10,1.53146e-10); }
		if(ieta==2){ fun_val_bdF[ipar][ieta]->SetParameters(1.79748e+03,1.39612e-02,-5.32411e-07,-2.99721e-09,4.19168e-10); }
	}
	
	
	for(int ipar=0; ipar<5; ipar++){
		for(int jpar=0; jpar<5; jpar++){
			cov_bpfrat_val_beta[ieta](ipar,jpar) = COV_val_mc[ieta][ipar][jpar];
		}
	}
	
	
	sprintf(name,"BTagger_PFRate_Pol2_VR_Eta%i",ieta+1);
	bpfrat_val_pol2[ieta]= new TF1(name,Pol2,1000,5000,3);
	
	if(ieta==0) { bpfrat_val_pol2[ieta]->SetParameters(1.03646e-02,2.23546e-06,-5.14962e-10); }
	if(ieta==1) { bpfrat_val_pol2[ieta]->SetParameters(1.12766e-02,2.76393e-08,2.18502e-10); }
	if(ieta==2) { bpfrat_val_pol2[ieta]->SetParameters(1.37085e-02,-2.10553e-07,1.20776e-10); }
  	
	
	for(int ipar=0; ipar<3; ipar++){
		sprintf(name,"dFun_val_pol2_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_val_bdF_pol2[ipar][ieta]= new TF1(name,bdF0_pol2,1000,5000,3);}
		if(ipar==1){fun_val_bdF_pol2[ipar][ieta]= new TF1(name,bdF1_pol2,1000,5000,3);}
		if(ipar==2){fun_val_bdF_pol2[ipar][ieta]= new TF1(name,bdF2_pol2,1000,5000,3);}
	}
	
	for(int ipar=0; ipar<3; ipar++){
	
		if(ieta==0){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(1.03646e-02,2.23546e-06,-5.14962e-10); }
		if(ieta==1){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(1.12766e-02,2.76393e-08,2.18502e-10); }
		if(ieta==2){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(1.37085e-02,-2.10553e-07,1.20776e-10); }
	}
	
	for(int ipar=0; ipar<3; ipar++){
		for(int jpar=0; jpar<3; jpar++){
			cov_bpfrat_val_pol2[ieta](ipar,jpar) = COV_val_mc_pol2[ieta][ipar][jpar];
		}
	}
   
}

#endif


#ifdef Anal_2016

for(int ieta=0; ieta<netabins; ieta++){
	sprintf(name,"BTagger_PFRate_SR_Eta%i",ieta+1);
	bpfrat_beta[ieta]= new TF1(name,BiFun,1000,5000,5);
	
	#ifdef B_Tight
	
	if(ieta==0) { bpfrat_beta[ieta]->SetParameters(2.38773e+03,3.26989e-02,-2.98803e-06,4.69225e-09,-2.05419e-09); }
	if(ieta==1) { bpfrat_beta[ieta]->SetParameters(2.52028e+03,3.27123e-02,-6.72983e-07,3.46661e-09,-1.35889e-09); }
	if(ieta==2) { bpfrat_beta[ieta]->SetParameters(2.52488e+03,2.75470e-02,-5.06421e-06,-2.47926e-09,9.59098e-10); }
	
	#else
		
	if(ieta==0) { bpfrat_beta[ieta]->SetParameters(1.95962e+03,6.98537e-02,1.18804e-05,2.03276e-08,-4.51380e-09); }
	if(ieta==1) { bpfrat_beta[ieta]->SetParameters(1.69468e+03,7.28344e-02,1.89960e-06,-1.07402e-08,2.61903e-09); }
	if(ieta==2) { bpfrat_beta[ieta]->SetParameters(2.29001e+03,7.84576e-02,6.11443e-06,1.81058e-09,-4.64201e-09); }
	
	#endif
	
	for(int ipar=0; ipar<5; ipar++){
		sprintf(name,"dFun_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_bdF[ipar][ieta]= new TF1(name,bdF0,1000,5000,5);}
		if(ipar==1){fun_bdF[ipar][ieta]= new TF1(name,bdF1,1000,5000,5);}
		if(ipar==2){fun_bdF[ipar][ieta]= new TF1(name,bdF2,1000,5000,5);}
		if(ipar==3){fun_bdF[ipar][ieta]= new TF1(name,bdF3,1000,5000,5);}
		if(ipar==4){fun_bdF[ipar][ieta]= new TF1(name,bdF4,1000,5000,5);}
	}
	for(int ipar=0; ipar<5; ipar++){
		
		#ifdef B_Tight
		
		if(ieta==0){ fun_bdF[ipar][ieta]->SetParameters(2.38773e+03,3.26989e-02,-2.98803e-06,4.69225e-09,-2.05419e-09); }
		if(ieta==1){ fun_bdF[ipar][ieta]->SetParameters(2.52028e+03,3.27123e-02,-6.72983e-07,3.46661e-09,-1.35889e-09); }
		if(ieta==2){ fun_bdF[ipar][ieta]->SetParameters(2.52488e+03,2.75470e-02,-5.06421e-06,-2.47926e-09,9.59098e-10); }
	
		#else
		
		if(ieta==0){ fun_bdF[ipar][ieta]->SetParameters(1.95962e+03,6.98537e-02,1.18804e-05,2.03276e-08,-4.51380e-09); }
		if(ieta==1){ fun_bdF[ipar][ieta]->SetParameters(1.69468e+03,7.28344e-02,1.89960e-06,-1.07402e-08,2.61903e-09); }
		if(ieta==2){ fun_bdF[ipar][ieta]->SetParameters(2.29001e+03,7.84576e-02,6.11443e-06,1.81058e-09,-4.64201e-09); }
	
		#endif
	}
	
	for(int ipar=0; ipar<5; ipar++){
		for(int jpar=0; jpar<5; jpar++){
			cov_bpfrat_beta[ieta](ipar,jpar) = COV_sig_mc[ieta][ipar][jpar];
		}
	}
	
		
	sprintf(name,"BTagger_PFRate_Pol2_SR_Eta%i",ieta+1);
	bpfrat_pol2[ieta]= new TF1(name,Pol2,1000,5000,3);
	
  	
  	#ifdef B_Tight
  	
  	if(ieta==0) { bpfrat_pol2[ieta]->SetParameters(5.03656e-02,-7.99279e-06,2.79968e-10); }
	if(ieta==1) { bpfrat_pol2[ieta]->SetParameters(4.54351e-02,-6.33923e-06,4.95724e-10); }
	if(ieta==2) { bpfrat_pol2[ieta]->SetParameters(3.32082e-02,-1.65095e-06,-2.52107e-10); }
	
	#else
	
	if(ieta==0) { bpfrat_pol2[ieta]->SetParameters(5.97915e-02,7.62850e-06,-8.32615e-10); }
	if(ieta==1) { bpfrat_pol2[ieta]->SetParameters(7.10220e-02,-2.45082e-06,1.84770e-09); }
	if(ieta==2) { bpfrat_pol2[ieta]->SetParameters(5.65278e-02,1.67292e-05,-3.03467e-09); }
	
	#endif
	
	for(int ipar=0; ipar<3; ipar++){
		sprintf(name,"dFun_pol2_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_bdF_pol2[ipar][ieta]= new TF1(name,bdF0_pol2,1000,5000,3);}
		if(ipar==1){fun_bdF_pol2[ipar][ieta]= new TF1(name,bdF1_pol2,1000,5000,3);}
		if(ipar==2){fun_bdF_pol2[ipar][ieta]= new TF1(name,bdF2_pol2,1000,5000,3);}
	}
	
	
	#ifdef B_Tight
	
	for(int ipar=0; ipar<3; ipar++){
		if(ieta==0){ fun_bdF_pol2[ipar][ieta]->SetParameters(5.03656e-02,-7.99279e-06,2.79968e-10); }
		if(ieta==1){ fun_bdF_pol2[ipar][ieta]->SetParameters(4.54351e-02,-6.33923e-06,4.95724e-10); }
		if(ieta==2){ fun_bdF_pol2[ipar][ieta]->SetParameters(3.32082e-02,-1.65095e-06,-2.52107e-10); }
	}

	#else
	
	for(int ipar=0; ipar<3; ipar++){
		
		if(ieta==0){ fun_bdF_pol2[ipar][ieta]->SetParameters(5.97915e-02,7.62850e-06,-8.32615e-10); }
		if(ieta==1){ fun_bdF_pol2[ipar][ieta]->SetParameters(7.10220e-02,-2.45082e-06,1.84770e-09); }
		if(ieta==2){ fun_bdF_pol2[ipar][ieta]->SetParameters(5.65278e-02,1.67292e-05,-3.03467e-09); }
		
	}
	
	#endif

	for(int ipar=0; ipar<3; ipar++){
		for(int jpar=0; jpar<3; jpar++){
			cov_bpfrat_pol2[ieta](ipar,jpar) = COV_sig_mc_pol2[ieta][ipar][jpar];
		}
	}
 
}

bpfrat_val = new TF1("BTagger_PFRate_VR",MikkoFunc1,400,2000,4);
bpfrat_val->SetParameters(1.25644e+01,-6.72034e-01,-1.34758e+01,1.22587e-02);

for(int ieta=0; ieta<netabins; ieta++){
	sprintf(name,"BTagger_PFRate_VR_Eta%i",ieta+1);
	bpfrat_val_beta[ieta]= new TF1(name,BiFun,500,2000,5);
	
	#ifdef B_Tight
	
	if(ieta==0) { bpfrat_val_beta[ieta]->SetParameters(2.34538e+03,1.30676e-02,3.29252e-06,1.27183e-09,-2.45126e-09); }
	if(ieta==1) { bpfrat_val_beta[ieta]->SetParameters(1.66768e+03,1.40119e-02,2.76997e-06,-3.59449e-09,-3.71193e-10); }
	if(ieta==2) { bpfrat_val_beta[ieta]->SetParameters(1.74832e+03,1.60425e-02,1.07590e-06,-6.46407e-09,-5.12396e-10); }
	
	#else

	if(ieta==0) { bpfrat_val_beta[ieta]->SetParameters(1.61913e+03,2.34176e-02,1.49306e-05,2.48973e-08,-7.50976e-10); }
	if(ieta==1) { bpfrat_val_beta[ieta]->SetParameters(1.70283e+03,3.40288e-02,1.12081e-05,-8.28958e-09,3.94102e-10); }
	if(ieta==2) { bpfrat_val_beta[ieta]->SetParameters(1.80500e+03,4.91733e-02,6.78692e-06,-1.42749e-08,-1.16174e-09); }
	
	#endif
	
	for(int ipar=0; ipar<5; ipar++){
		sprintf(name,"dFun_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_val_bdF[ipar][ieta]= new TF1(name,bdF0,1000,5000,5);}
		if(ipar==1){fun_val_bdF[ipar][ieta]= new TF1(name,bdF1,1000,5000,5);}
		if(ipar==2){fun_val_bdF[ipar][ieta]= new TF1(name,bdF2,1000,5000,5);}
		if(ipar==3){fun_val_bdF[ipar][ieta]= new TF1(name,bdF3,1000,5000,5);}
		if(ipar==4){fun_val_bdF[ipar][ieta]= new TF1(name,bdF4,1000,5000,5);}
	}
	for(int ipar=0; ipar<5; ipar++){
		
		#ifdef B_Tight
		
		if(ieta==0){ fun_val_bdF[ipar][ieta]->SetParameters(2.34538e+03,1.30676e-02,3.29252e-06,1.27183e-09,-2.45126e-09); }
		if(ieta==1){ fun_val_bdF[ipar][ieta]->SetParameters(1.66768e+03,1.40119e-02,2.76997e-06,-3.59449e-09,-3.71193e-10); }
		if(ieta==2){ fun_val_bdF[ipar][ieta]->SetParameters(1.74832e+03,1.60425e-02,1.07590e-06,-6.46407e-09,-5.12396e-10); }
		
		#else
		
		if(ieta==0){ fun_val_bdF[ipar][ieta]->SetParameters(1.61913e+03,2.34176e-02,1.49306e-05,2.48973e-08,-7.50976e-10); }
		if(ieta==1){ fun_val_bdF[ipar][ieta]->SetParameters(1.70283e+03,3.40288e-02,1.12081e-05,-8.28958e-09,3.94102e-10); }
		if(ieta==2){ fun_val_bdF[ipar][ieta]->SetParameters(1.80500e+03,4.91733e-02,6.78692e-06,-1.42749e-08,-1.16174e-09); }
		
		#endif
	}
	
	
	for(int ipar=0; ipar<5; ipar++){
		for(int jpar=0; jpar<5; jpar++){
			cov_bpfrat_val_beta[ieta](ipar,jpar) = COV_val_mc[ieta][ipar][jpar];
		}
	}
	
	
	sprintf(name,"BTagger_PFRate_Pol2_VR_Eta%i",ieta+1);
	bpfrat_val_pol2[ieta]= new TF1(name,Pol2,1000,5000,3);

  	#ifdef B_Tight
  	
  	if(ieta==0) { bpfrat_val_pol2[ieta]->SetParameters(5.28523e-03,5.42903e-06,-9.21253e-10); }
	if(ieta==1) { bpfrat_val_pol2[ieta]->SetParameters(6.85565e-03,5.13629e-06,-5.62590e-10); }
	if(ieta==2) { bpfrat_val_pol2[ieta]->SetParameters(9.22647e-03,5.23730e-06,-8.92990e-10); }
	
	#else
	
	if(ieta==0) { bpfrat_val_pol2[ieta]->SetParameters(1.00135e-02,7.02723e-06,1.17011e-09); }
	if(ieta==1) { bpfrat_val_pol2[ieta]->SetParameters(1.02765e-02,1.43040e-05,-3.73468e-10); }
	if(ieta==2) { bpfrat_val_pol2[ieta]->SetParameters(2.27449e-02,1.83936e-05,-2.37152e-09); }
	
	#endif
	
	for(int ipar=0; ipar<3; ipar++){
		sprintf(name,"dFun_val_pol2_%i_Eta%i",ipar,ieta+1);
		if(ipar==0){fun_val_bdF_pol2[ipar][ieta]= new TF1(name,bdF0_pol2,1000,5000,3);}
		if(ipar==1){fun_val_bdF_pol2[ipar][ieta]= new TF1(name,bdF1_pol2,1000,5000,3);}
		if(ipar==2){fun_val_bdF_pol2[ipar][ieta]= new TF1(name,bdF2_pol2,1000,5000,3);}
	}
	
	
	#ifdef B_Tight
	
	for(int ipar=0; ipar<3; ipar++){
		if(ieta==0){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(5.28523e-03,5.42903e-06,-9.21253e-10); }
		if(ieta==1){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(6.85565e-03,5.13629e-06,-5.62590e-10); }
		if(ieta==2){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(9.22647e-03,5.23730e-06,-8.92990e-10); }
	}
	
	#else
	
	for(int ipar=0; ipar<3; ipar++){
		if(ieta==0){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(1.00135e-02,7.02723e-06,1.17011e-09); }
		if(ieta==1){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(1.02765e-02,1.43040e-05,-3.73468e-10); }
		if(ieta==2){ fun_val_bdF_pol2[ipar][ieta]->SetParameters(2.27449e-02,1.83936e-05,-2.37152e-09); }
	}
	
	#endif
	
	for(int ipar=0; ipar<3; ipar++){
		for(int jpar=0; jpar<3; jpar++){
			cov_bpfrat_val_pol2[ieta](ipar,jpar) = COV_val_mc_pol2[ieta][ipar][jpar];
		}
	}
   
}

#endif

tagpfrat_tau32 = new TF1("TopTagger_Tau32_PFRate",pol2,20,400,3);
tagpfrat_tau32->SetParameters(2.16355e-03,-8.54422e-05,4.06720e-06);

tagpfrat_tau32_val = new TF1("TopTagger_Tau32_PFRate_VR",pol2,20,400,3);
tagpfrat_tau32_val->SetParameters(1.03655e-02,-2.63081e-04,4.70500e-06);

tagpfrat_MDAK8 = new TF1("TopTagger_MDDeepAK8_PFRate",Parabol,50,400,4);
//tagpfrat_MDAK8->SetParameters(1.37107e+02,6.80223e-03,2.45669e-07,8.88678e-07); // pt uncorrected
tagpfrat_MDAK8->SetParameters(1.27054e+02,2.60026e-02,-5.63565e-07,1.43812e-06);// pt corrected

tagpfrat_MDAK8_val = new TF1("TopTagger_MDDeepAK8_PFRate_VR",Parabol,50,400,4);
//tagpfrat_MDAK8_val->SetParameters(1.33008e+02,5.04849e-03,-8.75128e-08,7.01839e-07); // pt uncorrected
tagpfrat_MDAK8_val->SetParameters(1.25787e+02,9.81153e-03,-2.22362e-07,8.19089e-07);// pt corrected

tagpfrat_MDAK8_val2= new TF1("TopTagger_MDDeepAK8_PFRate_VR2",pol2,20,400,3);
tagpfrat_MDAK8_val2->SetParameters(1.04260e-02,-1.07329e-04,3.99638e-07);

}

sprintf(name,"$LHAPDF_DATA_PATH/cteq66");
LHAPDF::initPDFSet(name, LHAPDF::LHGRID, 0);
//LHAPDF::PDFSet set(name);
 
}

Bool_t Anal_Nano_PROOF::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either Anal_Nano_PROOF::GetEntry() or TBranch::GetEntry()
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
 TString str1;

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
	   Jet_DeepCSVMCeff[fjet] = 1;//BTag_MCEfficiency("deepcsvm",int(Jet_hadronFlavour[fjet]),Jet_pt[fjet],fabs(Jet_eta[fjet]));
	
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
   
//  if(!(Pileup_nTrueInt>=20 && Pileup_nTrueInt<40)) return kFALSE;
   
   weight = 1.0;
   puWeight = 1.0;
   
   if(isMC){
	   
      weight = Generator_weight;
      nevent_total += 1;
      weightev = Generator_weight;
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
   
   bool mu_trig = (HLT_IsoMu27 || HLT_Mu50);
  
   #if defined(Anal_2017) || defined(Anal_2018)
   bool had_trig = (HLT_AK8PFJet420_TrimMass30 || HLT_AK8PFHT900_TrimMass50 || HLT_PFJet500 || HLT_AK8PFJet500 || HLT_PFHT1050);
   #endif
   
   #ifdef Data_2017B
    had_trig = (HLT_PFJet500 || HLT_AK8PFJet500 || HLT_PFHT1050);
   #endif
   
   #ifdef Anal_2016
   bool had_trig = (HLT_AK8DiPFJet280_200_TrimMass30 || 
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
      
   if((nMuon  > 0) && (Muon_pt[0]>lepptcut))  return kFALSE;
   if((nElectron  > 0) && (Electron_pt[0]>lepptcut)) return kFALSE;
   if((nPhoton  > 0) && (Photon_pt[0]>lepptcut)) return kFALSE;
    
   if(!(had_trig)) return kFALSE; 
   
   if(nFatJet<1 || nJet<2) return kFALSE;
   
   if(isMC){
	weightpass = weight;
	Tout1->Fill();
   }

// pt fill to calculate b-tagging efficiency //
   
   if(isMC){
   
   for(unsigned ijet=0; ijet<nJet; ijet++){
	
	int ieta = getbinid(fabs(Jet_eta[ijet]),netarange,betarange);
	
	if(ieta>=0){
	
	if(Jet_genJetIdx[ijet]>=0){
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
		
//		hist_mass_qg->Fill(GenPart_mass[qg_d[ih]],weight);
//		str = TString::Format("entry %lli igen %i pid %i status %i pt %f\n",entry,qg_d[ih]+1,abs(GenPart_pdgId[qg_d[ih]]),GenPart_status[qg_d[ih]],GenPart_pt[qg_d[ih]]);
//		if(gProofServ) gProofServ->SendAsynMessage(str);
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
		 
		  if(FakeAnalysis && (!isMC || (isMC && isQCD))){
 
		  if(ittag==0 && isubbtag==1 && ibtag==0) { 
			  if(bjet.Pt()<500) { exwt = bpfrat->Eval(500); }
			  else{
			  if(bjet.Pt()>2000) { exwt = bpfrat->Eval(2000); }
			  else{
			  exwt = bpfrat->Eval(bjet.Pt());
				}
		      }
			}
			
		  if(ittag==1 && isubbtag==1 && ibtag==0) { exwt = bpfrat_val->Eval(bjet.Pt()); }
		  
		  if(ittag==1 && isubbtag==1 && ibtag==1) { exwt = tagpfrat_tau32->Eval(FatJet_msoftdrop[topjetAK8]); }
		  if(ittag==1 && isubbtag==1 && ibtag==0) { exwt2 = tagpfrat_tau32_val->Eval(FatJet_msoftdrop[topjetAK8]); }
		  
		  }
	  
		  weight1 = weight_t * sfwt_tau32 * exwt;
		  weight2 = weight_t * sfwt_tau32 * exwt2;
		  
		  if(Jet_pt[bjetAK4] > AK4ptcut_fi){
			  
			  if((Jet_MatchFatJet[bjetAK4]>=0)&&(FatJet_msoftdrop[Jet_MatchFatJet[bjetAK4]] < bmasscut_fun(FatJet_pt[Jet_MatchFatJet[bjetAK4]]))){
			 /*
				int isdmassbin = getbinid(FatJet_msoftdrop[Jet_MatchFatJet[bjetAK4]],nsdmassbins,sdmassbins);
				if(isdmassbin>=0 && isdmassbin<nsdmassbins){
					if(isMC && (isTTBar||isST)){
						weight1 *= cortt_sdmass[isdmassbin];
						weight2 *= cortt_sdmass[isdmassbin];
					}
				}
		  */
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
		  /*
			int isdmassbin = getbinid(FatJet_msoftdrop[Jet_MatchFatJet[bjetAK4]],nsdmassbins,sdmassbins);
			if(isdmassbin>=0 && isdmassbin<nsdmassbins){
				if(isMC && !isQCD && !isSignal){
						weight1 *= cortt_sdmass[isdmassbin];
						weight2 *= cortt_sdmass[isdmassbin];
					}
				}
		  */
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
		  
//		  int imsdbin = getbinid(FatJet_msoftdrop[topjetAK8],ntopsdmassbins,topsdmassbins);

		  if(FakeAnalysis){
		  
//		  if(imdttag==0 && isubbtag==1 && ibtag==0) { exwt = bpfrat->Eval(bjet.Pt()); }
//		  if(imdttag==1 && isubbtag==1 && ibtag==0) { exwt = bpfrat_val->Eval(bjet.Pt()); }

		  if((imdttag==0||imdttag==4) && ibtag==0) {
			   
			  exwt = bpfrat_beta[ieta]->Eval(m_tb);   // taking pf ratio from mtb
			  
				for(int ivar=0; ivar<5; ivar++){
					rvar[0][ivar] = fun_bdF[ivar][ieta]->Eval(m_tb);
					cvar[ivar][0] = fun_bdF[ivar][ieta]->Eval(m_tb);
				}

				err1 = (cov_bpfrat_beta[ieta] * cvar);
				err = (rvar * err1);
			    exwt_sys = err(0,0);
				exwt_sys = sqrt(exwt_sys);
				
				 exwt2 = bpfrat_pol2[ieta]->Eval(m_tb);   // taking pf ratio from mtb 
				
				for(int ivar=0; ivar<3; ivar++){
					rvar2[0][ivar] = fun_bdF_pol2[ivar][ieta]->Eval(m_tb);
					cvar2[ivar][0] = fun_bdF_pol2[ivar][ieta]->Eval(m_tb);
				}

				err2 = (cov_bpfrat_pol2[ieta] * cvar2);
				err_pol2 = (rvar2 * err2);
			    exwt_sys2 = err_pol2(0,0);
				exwt_sys2 = sqrt(exwt_sys2);
			 
			}
			  
		  if((imdttag==1||imdttag==5) && ibtag==0)  { 
	 
			  exwt = bpfrat_val_beta[ieta]->Eval(m_tb);                // taking pf ratio from mtb
			  
			  for(int ivar=0; ivar<5; ivar++){
//				  fun_val_bdF[ivar][ieta]->Eval(bjet.Pt());
				  rvar[0][ivar] = fun_val_bdF[ivar][ieta]->Eval(m_tb);
				  cvar[ivar][0] = fun_val_bdF[ivar][ieta]->Eval(m_tb);
			  }
			  
				err1 = (cov_bpfrat_val_beta[ieta] * cvar);
				err =  (rvar * err1);
				exwt_sys = err(0,0);
				exwt_sys = sqrt(exwt_sys);
	
			  
			  exwt2 = bpfrat_val_pol2[ieta]->Eval(m_tb);                // taking pf ratio from mtb
			 
			  for(int ivar=0; ivar<3; ivar++){
				  rvar2[0][ivar] = fun_val_bdF_pol2[ivar][ieta]->Eval(m_tb);
				  cvar2[ivar][0] = fun_val_bdF_pol2[ivar][ieta]->Eval(m_tb);
			  }
			  
				err2 = (cov_bpfrat_val_pol2[ieta] * cvar2);
				err_pol2 =  (rvar2 * err2);
				exwt_sys2 = err_pol2(0,0);
				exwt_sys2 = sqrt(exwt_sys2);
			
		  }
		  
		  if(exwt<0) { exwt = 1.e-5; }//0.001; }
		  if(exwt2<0) { exwt2 = 1.e-5; }//0.001; }
		  
	
			if(Jet_MatchFatJet[bjetAK4]>=0){
			
			int imsd_bin = getbinid(FatJet_msoftdrop[Jet_MatchFatJet[bjetAK4]],nsdmassbins,sdmassbins);
			
			if(imsd_bin>=0&&(imdttag==0||imdttag==4)&&(ibtag==0))
				{	
				cor_bsd = (isMC)?(cor_b_sdmass_mc[0][imsd_bin]):(cor_b_sdmass_data[0][imsd_bin]);
				}
			
			if(imsd_bin>=0&&(imdttag==1||imdttag==5)&&(ibtag==0))
				{
				cor_bsd = (isMC)?(cor_b_sdmass_mc[1][imsd_bin]):(cor_b_sdmass_data[1][imsd_bin]);
				}
			}

		}
		
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
					double weight_pdf = weight1*ypdf1*ypdf2*1./(Generator_xpdf1*Generator_xpdf2);
				
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
	  	  
   return kTRUE;
}

void Anal_Nano_PROOF::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

	fileOut->cd();
	fileOut->Write();

	fOutput->Add(OutFile);

	fileOut->Close();

}

void Anal_Nano_PROOF::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
