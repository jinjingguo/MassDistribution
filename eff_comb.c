#include "Utilities/Ntuple/VertexCompositeTree.h"
#include "tnp_weight_lowPt.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TEfficiency.h"
#include <iostream>
#include <string>
#include <map>


bool triggerMuonAcceptance(const double& pt, const double& eta)
{
  return ( fabs(eta) < 2.4 &&
        (    ( fabs(eta) < 1.2 && pt >= 3.3 ) ||
           (  1.2 <= fabs(eta) && fabs(eta) < 2.1 && pt >= 3.93-1.11*fabs(eta)) ||
             (  2.1 <= fabs(eta) && fabs(eta) < 2.4 && pt >= 1.3)
           )
       );
}

bool muonAcceptance(const double& pt, const double& eta)
{
  return ( fabs(eta) < 2.4 &&
        (    ( fabs(eta) < 0.8 && pt >= 3.3 ) ||
           (  0.8 <= fabs(eta) && fabs(eta) < 1.5 && pt >= 5.81-3.14*fabs(eta)) ||
             (  1.5 <= fabs(eta) && fabs(eta) < 2.4 && pt >= 0.8 && 1.89-0.526*fabs(eta) )
           )
       );
}

void GetEff(std::string inputFile1, std::string inputFile2, double pTBin[], int size, bool corrEff, bool muonTrig, double y1, double y2, const char* title, const char* outPath1, const char* outPath2, const char* outPath3, const char* outPath4)
{
  const auto& treeDir = "dimucontana_mc"; // For MC use dimucontana_mc

  // Load the input trees
  std::map<std::string, VertexCompositeTree> treeM;
  if (!treeM["pPb"].GetTree(inputFile1, treeDir)||!treeM["Pbp"].GetTree(inputFile2, treeDir)) { std::cout << "Invalid tree!" << std::endl; return; }

  // Initialize the canvas
  TCanvas c(title, title);
  c.Divide(3,2);
 
  // Initialize the efficiencies
  std::map<std::string, long int> nevents;
  std::map<std::string, std::unique_ptr<TEfficiency>> eff;
  for(auto& t : treeM)
  {
    auto& s = t.first;
    auto& tree = t.second;
    
    // Fill the efficiency histograms
    TH1D hPassed(Form("hPassed_%s",s.c_str()), s.c_str(), size, pTBin);
    TH1D hTotal(Form("hTotal_%s",s.c_str()), s.c_str(), size, pTBin);
    hPassed.Sumw2();
    hTotal.Sumw2();
    
    // Loop over the tree
    nevents[s] = tree.GetEntries();
    for(Long64_t jentry=0; jentry<nevents[s]; jentry++)
    { 
      if(!(jentry % 1000000)) std::cout<<"Processed "<<jentry<<" events out of "<<nevents[s]<<std::endl;
      if (tree.GetEntry(jentry)<0) { std::cout << "Invalid entry!" << std::endl;}
      
      // Loop over the generated candidates
      for(uint iGen=0; iGen<tree.candSize_gen(); iGen++){
        // Check we are in the rapidity bin of interest
        const auto& y_gen = (s=="Pbp" ? -1.0 : 1.0) * tree.y_gen()[iGen];
        if(fabs(y_gen)>y1 && fabs(y_gen)<=y2) continue;
      
        // Check that both generated muons are inside the single muon acceptance kinematic region
        const auto& pTD1_gen = tree.pTD1_gen()[iGen];
        const auto& etaD1_gen = (s=="Pbp" ? -1.0 : 1.0) * tree.EtaD1_gen()[iGen];
        const bool mu1InAccep_gen = (muonTrig ? triggerMuonAcceptance(pTD1_gen, etaD1_gen) : muonAcceptance(pTD1_gen, etaD1_gen));
        const auto& pTD2_gen = tree.pTD2_gen()[iGen];
        const auto& etaD2_gen = (s=="Pbp" ? -1.0 : 1.0) * tree.EtaD2_gen()[iGen];
        const bool mu2InAccep_gen = (muonTrig ? triggerMuonAcceptance(pTD2_gen, etaD2_gen) : muonAcceptance(pTD2_gen, etaD2_gen));
        if(!mu1InAccep_gen || !mu2InAccep_gen) continue;

	// Fill the total histogram
        const auto& pT_gen = tree.pT_gen()[iGen];
        hTotal.Fill(pT_gen);

	const auto& iReco = tree.RecIdx_gen()[iGen];
        if(iReco<0) continue;

	// Check that both reconstructed muons are inside the single muon acceptance kinematic region
        const auto& pTD1 = tree.pTD1()[iReco];
        const auto& etaD1 = (s=="Pbp" ? -1.0 : 1.0) * tree.EtaD1()[iReco];
        const bool mu1InAccep = (muonTrig ? triggerMuonAcceptance(pTD1, etaD1) : muonAcceptance(pTD1, etaD1));
        const auto& pTD2 = tree.pTD2()[iReco];
        const auto& etaD2 = (s=="Pbp" ? -1.0 : 1.0) * tree.EtaD2()[iReco];
        const bool mu2InAccep = (muonTrig ? triggerMuonAcceptance(pTD2, etaD2) : muonAcceptance(pTD2, etaD2));
        if(!mu1InAccep || !mu2InAccep) continue;
      
        // Check the passing condition
        const bool softCand = tree.softCand(iReco);
        const bool goodEvt = tree.evtSel()[0];
        const bool trigCand = (muonTrig ? (tree.trigHLT()[0] && tree.trigCand(0, iReco)) : true);
        const bool passed = (goodEvt && softCand && trigCand);

	// Derive the efficiency correction
	double effCorr = 1.0;
	if (corrEff)
  	{
	  const auto& trgCorr = (muonTrig ? (tnp_weight_trg_ppb(pTD1, etaD1, 0) * tnp_weight_trg_ppb(pTD2, etaD2, 0)) : 1.0);
	  const auto& trkCorr = 1.0;//tnp_weight_trk_ppb(pTD1, etaD1, 0) * tnp_weight_trk_ppb(pTD2, etaD2, 0);
	  const auto& muidCorr = tnp_weight_muid_ppb(pTD1, etaD1, -10) * tnp_weight_muid_ppb(pTD2, etaD2, -10);
	  effCorr = trgCorr*trkCorr*muidCorr;
	}
      
        // Fill the passing histogram
        if(passed) { hPassed.Fill(pT_gen, effCorr); }
      }
    }
      
    // Set the histograms to the efficiency
    eff[s].reset(new TEfficiency(Form("eff_%s",s.c_str()), Form("eff_%s",s.c_str()), size, pTBin));
    eff[s]->SetTotalHistogram(hTotal, "f");
    eff[s]->SetPassedHistogram(hPassed, "f");
    if(TEfficiency::CheckWeights(hPassed, hTotal)) { eff[s]->SetStatisticOption(TEfficiency::kBJeffrey); }
  }

  auto hPassed = *((TH1D*)eff["pPb"]->GetPassedHistogram());
  hPassed.Add(eff["pPb"]->GetPassedHistogram(), eff["Pbp"]->GetPassedHistogram(), 110.8/nevents["pPb"], 62.6/nevents["Pbp"]);

  auto hTotal = *((TH1D*)eff["pPb"]->GetTotalHistogram());
  hTotal.Add(eff["pPb"]->GetTotalHistogram(), eff["Pbp"]->GetTotalHistogram(), 110.8/nevents["pPb"], 62.6/nevents["Pbp"]);

  TEfficiency eff_comb("eff_comb", "eff_comb", size, pTBin);
  eff_comb.SetTotalHistogram(hTotal, "f");
  eff_comb.SetPassedHistogram(hPassed, "");
  if(TEfficiency::CheckWeights(hPassed, hTotal)) { eff_comb.SetStatisticOption(TEfficiency::kBJeffrey); }

  auto hDivide1 = *((TH1D*)eff["pPb"]->GetPassedHistogram());
  hDivide1 = hDivide1*hTotal;
  auto hTotal1 = *((TH1D*)eff["pPb"]->GetTotalHistogram());
  hTotal1 = hTotal1*hPassed;
  auto hDivide2 = *((TH1D*)eff["Pbp"]->GetPassedHistogram());
  hDivide2 = hDivide2*hTotal;
  auto hTotal2 = *((TH1D*)eff["Pbp"]->GetTotalHistogram());
  hTotal2 = hTotal2*hPassed;

  c.cd(1);
  eff["pPb"]->Draw();
  c.cd(2);
  eff["Pbp"]->Draw();
  c.cd(3);
  eff_comb.Draw();
  TFile canvasFile(outPath1, "RECREATE");
  canvasFile.cd();
  c.Write("canvas");
  canvasFile.Write();
  canvasFile.Close();
  //c.Print(outputFile);

  TFile outFile1(outPath2, "RECREATE");
  outFile1.cd();
  eff["pPb"]->Write("eff_pPb");
  outFile1.Write();
  outFile1.Close();
  TFile outFile2(outPath3, "RECREATE");
  outFile2.cd();
  eff["Pbp"]->Write("eff_Pbp");
  outFile2.Write();
  outFile2.Close();
  TFile outFile3(outPath4, "RECREATE");
  outFile3.cd();
  eff_comb.Write("eff_comb");
  outFile3.Write();
  outFile3.Close();

  return;
}
void eff_comb(){
  auto inputFile1 = "/storage1/users/wl33/DiMuTrees/pPb2016/Tree/VertexCompositeTree_JPsiToMuMu_pPb-Bst_pPb816Summer16_DiMuMC.root"; 
  auto inputFile2 = "/storage1/users/wl33/DiMuTrees/pPb2016/Tree/VertexCompositeTree_JPsiToMuMu_PbP-Bst_pPb816Summer16_DiMuMC.root";
  double pTBin[] = {6.5, 9.0, 50.0};
  int size = 2;
  bool muonTrig = true;
  bool corrEff = true;
  auto title = "PromptWTJPsi0_1.4";
  double y1 = 0.0;
  double y2 = 1.4;
  auto outPath1 = "./rootFile/PromptWT/JPsi/0_1.4/PromptWTJPsi0_1.4_canvas.root";
  auto outPath2 = "./rootFile/PromptWT/JPsi/0_1.4/PromptWTJPsi0_1.4_eff_pPb.root";
  auto outPath3 = "./rootFile/PromptWT/JPsi/0_1.4/PromptWTJPsi0_1.4_eff_Pbp.root";
  auto outPath4 = "./rootFile/PromptWT/JPsi/0_1.4/PromptWTJPsi0_1.4_eff_comb.root";
  GetEff(inputFile1, inputFile2, pTBin, size, corrEff, muonTrig, y1, y2, title, outPath1, outPath2, outPath3, outPath4);
}

