#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
//#include "CMS_lumi.C"

void LegendSettings(TLegend *leg, int ncolumn){
  leg->SetNColumns(ncolumn);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetLineColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->SetTextFont(42);
}

TGraphAsymmErrors* getMetEfficiency(std::string pathName_, TString inputFileName_, int ifile)
{
  std::cout << "coming inside  :  "  << ifile << std::endl;
  float pTbins_[] = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200.,
		     210., 220., 230., 240., 250., 260., 270., 280., 290., 300., 325., 350., 375., 400., 450., 600., 775., 1000.};

  TH1F* histMETAll = new TH1F("METAll", "METAll", 100, 0., 1000.);
  //TH1F* histMETAll = new TH1F("METAll", "METAll", 38, pTbins_);
  histMETAll->Sumw2(); histMETAll->Reset();
  TH1F* histMETPass = new TH1F("METPass", "METPass", 100, 0., 1000.);
  //TH1F* histMETPass = new TH1F("METPass", "METPass", 38, pTbins_);
  histMETPass->Sumw2(); histMETPass->Reset();
  //TFile *file = new TFile("/tmp/mdjordje/hltJetMetNtuple_Run2018D_LAST.root");
  //TTree* tree = (TTree*)file->Get("hltJetMetNtuple/tree");


  TChain *tree = new TChain("hltJetMetNtuple/tree");
  tree->Add("file/hltJetMetNtuple.root");

  tree->SetBranchStatus("*", 1);

  float metPt_, metPhi_;
  std::vector<std::string>* triggerResults_ = new std::vector<std::string>();
  std::vector<std::string>*  triggerResults_rerun_ = new std::vector<std::string>();
  std::vector<float>* muonPt_ = new std::vector<float>();
  std::vector<float>* muonPx_ = new std::vector<float>();
  std::vector<float>* muonPy_ = new std::vector<float>();
  std::vector<float>* muonEta_ = new std::vector<float>();
  std::vector<float>* muonDxy_ = new std::vector<float>();
  std::vector<float>* muonDz_ = new std::vector<float>();
  std::vector<float>* muonR04SumChargedHadronPt_ = new std::vector<float>();
  std::vector<float>* muonR04SumNeutralHadronEt_ = new std::vector<float>();
  std::vector<float>* muonR04SumPhotonEt_ = new std::vector<float>();
  std::vector<float>* muonR04SumPUPt_ = new std::vector<float>();
  std::vector<bool>* muonIsICHEPMedium_ = new std::vector<bool>();
  std::vector<bool>* muonIsGlobal_ = new std::vector<bool>();
  std::vector<bool>* muonIsTracker_ = new std::vector<bool>();
  std::vector<bool>* muonIsPF_ = new std::vector<bool>();

  std::vector<float>* elecPt_ = new std::vector<float>();
  std::vector<float>* elecPx_ = new std::vector<float>();
  std::vector<float>* elecPy_ = new std::vector<float>();
  std::vector<float>* elecEta_ = new std::vector<float>();
  std::vector<float>* elecDxy_ = new std::vector<float>();
  std::vector<float>* elecDz_ = new std::vector<float>();
  std::vector<bool>* elec_pass_conversion_ = new std::vector<bool>();
  std::vector<float>* elecR03SumChargedHadronPt_ = new std::vector<float>();
  std::vector<float>* elecR03SumNeutralHadronEt_ = new std::vector<float>();
  std::vector<float>* elecR03SumPhotonEt_ = new std::vector<float>();
  std::vector<float>* elecR03SumPUPt_ = new std::vector<float>();
  /*std::vector<bool>* elec_mva_tight_Spring16_v1_ = new std::vector<bool>();
  std::vector<bool>* elec_mva_medium_Spring16_v1_ = new std::vector<bool>();  */
  std::vector<bool>* elec_cutId_tight_Fall17_ = new std::vector<bool>();
  bool passMETFilter_;
  std::vector<bool>* elec_mva_wp80_Iso_Fall17_ = new std::vector<bool>();
  std::vector<bool>* elec_mva_wp90_Iso_Fall17_ = new std::vector<bool>();
  std::vector<bool>* elec_mva_wpLoose_Iso_Fall17_ = new std::vector<bool>();

  //std::vector<long>* run_ = new std::vector<long>();
  //std::vector<long>* event_ = new std::vector<long>();
  //std::vector<long>* lumi_ = new std::vector<long>();
  ULong64_t run_, event_, lumi_;

  tree->SetBranchAddress("metPt", &metPt_);
  tree->SetBranchAddress("metPhi", &metPhi_);
  tree->SetBranchAddress("triggerResults", &triggerResults_);
  tree->SetBranchAddress("triggerResults_rerun",&triggerResults_rerun_);
  tree->SetBranchAddress("muonPt", &muonPt_);
  tree->SetBranchAddress("muonPx", &muonPx_);
  tree->SetBranchAddress("muonPy", &muonPy_);
  tree->SetBranchAddress("muonEta", &muonEta_);
  tree->SetBranchAddress("muonDxy", &muonDxy_);
  tree->SetBranchAddress("muonDz", &muonDz_);
  tree->SetBranchAddress("muonR04SumChargedHadronPt", &muonR04SumChargedHadronPt_);
  tree->SetBranchAddress("muonR04SumNeutralHadronEt", &muonR04SumNeutralHadronEt_);
  tree->SetBranchAddress("muonR04SumPhotonEt", &muonR04SumPhotonEt_);
  tree->SetBranchAddress("muonR04SumPUPt", &muonR04SumPUPt_);
  tree->SetBranchAddress("muonIsICHEPMedium", &muonIsICHEPMedium_);
  tree->SetBranchAddress("muonIsGlobal", &muonIsGlobal_);
  tree->SetBranchAddress("muonIsTracker", &muonIsTracker_);
  tree->SetBranchAddress("muonIsPF", &muonIsPF_);
  tree->SetBranchAddress("passMETFilter", &passMETFilter_);

  tree->SetBranchAddress("elecPt", &elecPt_);
  tree->SetBranchAddress("elecPx", &elecPx_);
  tree->SetBranchAddress("elecPy", &elecPy_);
  tree->SetBranchAddress("elecEta", &elecEta_);
  tree->SetBranchAddress("elecDxy", &elecDxy_);
  tree->SetBranchAddress("elecDz", &elecDz_);
  tree->SetBranchAddress("elec_pass_conversion", &elec_pass_conversion_);
  tree->SetBranchAddress("elecR03SumChargedHadronPt", &elecR03SumChargedHadronPt_);
  tree->SetBranchAddress("elecR03SumNeutralHadronEt", &elecR03SumNeutralHadronEt_);
  tree->SetBranchAddress("elecR03SumPhotonEt", &elecR03SumPhotonEt_);
  tree->SetBranchAddress("elecR03SumPUPt", &elecR03SumPUPt_);
  /*tree->SetBranchAddress("elec_mva_tight_Spring16_v1", &elec_mva_tight_Spring16_v1_);
  tree->SetBranchAddress("elec_mva_medium_Spring16_v1", &elec_mva_medium_Spring16_v1_);*/
  tree->SetBranchAddress("elec_cutId_tight_Fall17", &elec_cutId_tight_Fall17_);
  tree->SetBranchAddress("elec_mva_wp80_Iso_Fall17", &elec_mva_wp80_Iso_Fall17_);
  tree->SetBranchAddress("elec_mva_wp90_Iso_Fall17", &elec_mva_wp90_Iso_Fall17_);
  tree->SetBranchAddress("elec_mva_wpLoose_Iso_Fall17", &elec_mva_wpLoose_Iso_Fall17_);

  tree->SetBranchAddress("run", &run_);
  tree->SetBranchAddress("event", &event_);
  tree->SetBranchAddress("lumi", &lumi_);

//L1
  unsigned int  passedL1MET_, passedHLTCaloMET_, passedHLTCaloMETClean_;
//L1

//L1
  tree->SetBranchAddress("passedL1MET", &passedL1MET_);
  tree->SetBranchAddress("passedHLTCaloMET", &passedHLTCaloMET_);
  tree->SetBranchAddress("passedHLTCaloMETClean", &passedHLTCaloMETClean_);
//L1


  int nEntries = tree->GetEntries();
  std::cout << "nEntries   :   "  << nEntries << std::endl;
//  for (int n = 0; n < nEntries	; n++){
   for (int n = 0; n < nEntries; n++){


    tree->GetEntry(n);


    //L1
  // if(passedL1MET_ == 0) continue;
  //  if(passedHLTCaloMET_ == 0) continue;
  //  if(passedHLTCaloMETClean_ == 0) continue;
    //L1

//	if(!passMETFilter_) continue;

    //check singleElectron trigger
    bool passSingleElectron = false;
    for(size_t i = 0; i < triggerResults_->size(); i++){
      //if((*triggerResults_)[i].find("HLT_Ele27_WPTight_Gsf_v") != std::string::npos ||
        //(*triggerResults_)[i].find("HLT_Ele27_eta2p1_WPTight_Gsf_v") != std::string::npos)passSingleElectron = true;
      if((*triggerResults_)[i].find("HLT_Ele32_WPTight_Gsf_v") != std::string::npos)passSingleElectron = true;
    }
    if(!passSingleElectron) continue;
    //std::cout<<" old met "<<metPt_<<std::endl;
    //select events with at least one electron
    std::vector<float>selElecPx_; selElecPx_.clear();
    std::vector<float>selElecPy_; selElecPy_.clear();
    for(size_t ie = 0; ie < elecPt_->size(); ie++){
     if((*elecPt_)[ie] > 30. && fabs((*elecEta_)[ie]) < 2.5 &&
	(!(fabs((*elecEta_)[ie]) > 1.4442 && fabs((*elecEta_)[ie]) < 1.5660)) &&
	fabs((*elecDxy_)[ie]) < 0.05 && fabs((*elecDz_)[ie]) < 0.1 && (*elec_pass_conversion_)[ie] == true) {
	float iso = ( (*elecR03SumChargedHadronPt_)[ie] + std::max(0., (*elecR03SumNeutralHadronEt_)[ie] + (*elecR03SumPhotonEt_)[ie] - (0.5 * (*elecR03SumPUPt_)[ie])) )/(*elecPt_)[ie];
	if(iso < 0.15){
	  //if((*elec_mva_tight_Spring16_v1_)[ie] == true) {
	  if((*elec_mva_wp80_Iso_Fall17_)[ie] == true) {
	  //if((*elec_cutId_tight_Fall17_)[ie] == true) {
     	   selElecPx_.push_back((*elecPx_)[ie]);
	   selElecPy_.push_back((*elecPy_)[ie]);
	  }
	 //}
	}
     }
    }
    if(selElecPx_.size() == 0) continue;

    //find selected muons
    std::vector<float>rejectMuon_; rejectMuon_.clear();
    std::vector<float>selMuonPx_; selMuonPx_.clear();
    std::vector<float>selMuonPy_; selMuonPy_.clear();
    for(size_t im = 0; im < muonPt_->size(); im++){
      if((*muonPt_)[im] > 10. && (*muonIsGlobal_)[im] == true){
	rejectMuon_.push_back((*muonPx_)[im]);
	float iso = ( (*muonR04SumChargedHadronPt_)[im] + std::max(0., (*muonR04SumNeutralHadronEt_)[im] + (*muonR04SumPhotonEt_)[im] - (0.5 * (*muonR04SumPUPt_)[im])) )/(*muonPt_)[im];
	if((*muonPt_)[im] > 30. && fabs((*muonEta_)[im]) < 2.4 && fabs((*muonDxy_)[im]) < 0.05 && fabs((*muonDz_)[im]) < 0.1 && iso < 0.15){
	  if((*muonIsICHEPMedium_)[im] == true) {
	  selMuonPx_.push_back((*muonPx_)[im]);
	  selMuonPy_.push_back((*muonPy_)[im]);
	}
       }
      }
    }

    if(rejectMuon_.size() > 0) continue;
    TLorentzVector newMetP4_;
    newMetP4_.SetPtEtaPhiM(metPt_, 0, metPhi_, 0);
    float newMetPx_ = newMetP4_.Px();
    float newMetPy_ = newMetP4_.Py();
    //Remove muons from MET to make it similar to CaloMet
    //for(size_t im = 0; im < selMuonPx_.size(); im++){
    //  newMetPx_ += selMuonPx_[im];
    //  newMetPy_ += selMuonPy_[im];
    //}
    float newMetPt_ = sqrt(newMetPx_*newMetPx_ + newMetPy_*newMetPy_);
    //std::cout << "checking met" << newMetPt_ << std::endl;
    histMETAll->Fill(newMetPt_);


    bool passTrigger = false;
    bool passTrigger_2 = false;

    if (ifile==0){
   // std::cout << "triggerResults"   << std::endl;
    for(size_t i = 0; i < triggerResults_->size(); i++){
      if((*triggerResults_)[i].find(pathName_) != std::string::npos)passTrigger = true;
    }
   }
    else {
   //    std::cout << "triggerResults_rerun"   << std::endl;
       for(size_t i = 0; i < triggerResults_rerun_->size(); i++){
          //std::cout << "triggerResults_rerun"   << std::endl;
         if((*triggerResults_rerun_)[i].find(pathName_) != std::string::npos)passTrigger = true;
       }

     }
    if(passTrigger)histMETPass->Fill(newMetPt_);

    //if(selMuonPx_.size() > 0 && metPt_ > 100) std::cout<<" old met "<<metPt_<<" new met "<<newMetPt_<<std::endl;

    if(passTrigger == false){
	if(metPt_ < 800 && metPt_ > 400){
		std::cout<< "MET = " << metPt_ << " MET phi = " << metPhi_ << " Run = " << run_  << " Event = " << event_ << " Lumi = " << 			lumi_ << std::endl;
		}
	}
  }

  //std::cout << "hist integral:   "  << histMETPass->Integral()  << std::endl;

  TGraphAsymmErrors* histMETEff = new TGraphAsymmErrors(histMETPass, histMETAll);

  return histMETEff;
}

void plotMetEfficiency(std::vector<std::string> pathNames, std::vector<TString> filenames, std::vector<TString> legends, TString outFileName_, TString label_string)
{
   gROOT->ProcessLine(".L tdrstyle.C");
   gROOT->ProcessLine("setTDRStyle()");

  //gROOT->LoadMacro("CMS_lumi.C");

  //writeExtraText = true;       // if extra text
  //extraText  = "Preliminary";  // default extra text is "Preliminary"
  //lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  //lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  //lumi_sqrtS = "2016, 13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)

  int colors[] = { 2, 3, 6, 4, 8, 8};
  int style[] = {24, 21 ,21, 22, 23, 20};

  std::vector<TGraphAsymmErrors*> histMETEffs;
  histMETEffs.resize(filenames.size());
  for(size_t ifile = 0; ifile < filenames.size(); ifile++){
    std::cout << "ifile:  "  << ifile << std::endl;
    histMETEffs[ifile] = (TGraphAsymmErrors*)getMetEfficiency(pathNames[ifile], filenames[ifile], ifile);
    histMETEffs[ifile]->SetMarkerStyle(style[ifile]);
    histMETEffs[ifile]->SetMarkerColor(colors[ifile]);
    histMETEffs[ifile]->SetMarkerSize(1.0);
    histMETEffs[ifile]->SetLineColor(colors[ifile]);
  }


  TCanvas* c1 = new TCanvas();
  c1->SetGridx();
  c1->SetGridy();
//  c1->SetLogy();
  c1->cd();

  TH1F* histogram_base = new TH1F("histogram_base", "", 100, 1., 999.);
  histogram_base->SetTitle("");
  histogram_base->SetStats(false);
  histogram_base->SetMinimum(0.0);
  histogram_base->SetMaximum(1.2);
  histogram_base->GetXaxis()->SetTitle("Offline E_{T}^{miss} (GeV)");
  histogram_base->GetYaxis()->SetTitle("Efficiency");
  histogram_base->Draw("hist");

  for(size_t ih = 0; ih < histMETEffs.size(); ih++){
  //  std::cout << "integral :  "  << histMETEffs[ih]->Integral()  << std::endl;
    histMETEffs[ih]->Draw("Psame");
  }

TPaveText *t = new TPaveText(0.88,0.97,0.96,1.0,"NDC");
t->SetTextSize(0.03);
t->SetBorderSize(0);
t->SetFillColor(0);
t->SetTextFont(42);
t->SetTextAlign(33);
//t->AddText("2018(13 TeV)");
t->AddText("2018(13 TeV)");
t->SetX1NDC(0);
t->SetX2NDC(0.1);
t->SetY1NDC(1.12);
t->SetY2NDC(1.2);
t->Draw();

  TLegend *leg1 = new TLegend(0.5,0.4,0.9,(0.4+0.05*(histMETEffs.size())));
  LegendSettings(leg1, 1);
  for(size_t ih = 0; ih < histMETEffs.size(); ih++){
    leg1->AddEntry(histMETEffs[ih], legends[ih], "lp");
  }
  leg1->Draw();
/*
  TPaveText* label = new TPaveText(0.2, 0.85, 0.65, 0.9, "brNDC");
  label->AddText(label_string);
  label->SetFillColor(10);
  label->SetBorderSize(0);
  label->SetTextColor(1);
  label->SetTextAlign(12);
  label->SetTextSize(0.03);
  label->SetTextFont(42);
  label->Draw();
*/
  //CMS_lumi( c1, iPeriod, 11 );

  c1->Update();
  c1->SaveAs("plots/Efficiency_"+TString(outFileName_)+"_sel_overlay.png");
  c1->SaveAs("plots/Efficiency_"+TString(outFileName_)+"_sel_overlay.pdf");
  //c1->SaveAs("plots/Efficiency_"+TString(outFileName_)+"_sel.C");
}

void plotMetEfficiencyAll_overlap()
{

  std::vector<std::string> pathNames_PFMETTypeOne;
//  pathNames_PFMETTypeOne.push_back("HLT_PFMET200_HBHE_BeamHaloCleaned_v");
//  pathNames_PFMETTypeOne.push_back("HLT_PFMET200_NotCleaned_v");
//  pathNames_PFMETTypeOne.push_back("HLT_PFMET200_HBHECleaned_v");
  pathNames_PFMETTypeOne.push_back("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v");
  pathNames_PFMETTypeOne.push_back("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v");
//  pathNames_PFMETTypeOne.push_back("HLT_PFMET250_HBHECleaned_v");

  //pathNames_PFMETTypeOne.push_back("HLT_PFMET120_PFMHT120_IDTight_v");
  vector<TString> filenames_data_2016G_Type1;
  filenames_data_2016G_Type1.push_back("/afs/cern.ch/work/d/dekumar/public/JETMET_Work/newStructure/CMSSW_10_1_11_patch1/src/HLTrigger/Configuration/test/drawEfficiency/dummy");
  filenames_data_2016G_Type1.push_back("/afs/cern.ch/work/d/dekumar/public/JETMET_Work/newStructure/CMSSW_10_1_11_patch1/src/HLTrigger/Configuration/test/drawEfficiency/dummy");
//  filenames_data_2016G_Type1.push_back("../test/hltMetNtuple.root");
//  filenames_data_2016G_Type1.push_back("../test/hltMetNtuple.root");
//  filenames_data_2016G_Type1.push_back("../test/hltMetNtuple.root");

  vector<TString>legends_data_Type1;
//  legends_data_Type1.push_back("HLT E_{T}^{miss} > 200 GeV, BH Cleaned");
//  legends_data_Type1.push_back("HLT E_{T}^{miss} > 200 GeV, No Cleaning");
//  legends_data_Type1.push_back("HLT E_{T}^{miss} > 200 GeV, HBHE Cleaned");
  legends_data_Type1.push_back("With HBHE Noise Filter");
  legends_data_Type1.push_back("Without HBHE Noise Filter");

  plotMetEfficiency(pathNames_PFMETTypeOne, filenames_data_2016G_Type1, legends_data_Type1, "HLT_PFMET_2018D", "");

}
