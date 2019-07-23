#include "Utilities/Ntuple/VertexCompositeTree.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TEfficiency.h"
#include <iostream>
#include <string>


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

void GetAcc(std::string inputFile, double pTBin[], int size, bool corrAcc, bool muonTrig, double y1, double y2, const char* title, const char* outPath1, const char* outPath2)
{
  const auto& treeDir = "dimuana_mc"; // For MC use dimucontana_mc

  // Load the input trees
  VertexCompositeTree tree;
  if (!tree.GetTree(inputFile, treeDir)) { std::cout << "Invalid tree!" << std::endl; return; }

  // Initialize the canvas
  TCanvas c(title, title);
 
  // Initialize the acceptance
  std::unique_ptr<TEfficiency> acc;
    
  // Fill the acceptance histograms
  TH1D hPassed("hPassed", "hPassed", size, pTBin);
  TH1D hTotal("hTotal", "hTotal", size, pTBin);
  hPassed.Sumw2();
  hTotal.Sumw2();
    
  // Loop over the tree
  const auto& nevents = tree.GetEntries();
  for(Long64_t jentry=0; jentry<nevents; jentry++)
  {
    if(!(jentry % 1000000)) std::cout<<"Processed "<<jentry<<" events out of "<<nevents<<std::endl;
    if (tree.GetEntry(jentry)<0) { std::cout << "Invalid entry!" << std::endl;}
      
    // Loop over the generated candidates
    for(uint iGen=0; iGen<tree.candSize_gen(); iGen++){
      // Check we are in the rapidity bin of interest
      const auto& y_gen = tree.y_gen()[iGen];
      if(fabs(y_gen)>y1 && fabs(y_gen)<=y2) continue;
      
      // Fill the total histogram
      const auto& pT_gen = tree.pT_gen()[iGen];
      hTotal.Fill(pT_gen);
      
      // Check that both generated muons are inside the single muon acceptance kinematic region
      const auto& pTD1_gen = tree.pTD1_gen()[iGen];
      const auto& etaD1_gen = tree.EtaD1_gen()[iGen];
      const bool mu1InAccep_gen = (muonTrig ? triggerMuonAcceptance(pTD1_gen, etaD1_gen) : muonAcceptance(pTD1_gen, etaD1_gen));
      const auto& pTD2_gen = tree.pTD2_gen()[iGen];
      const auto& etaD2_gen = tree.EtaD2_gen()[iGen];
      const bool mu2InAccep_gen = (muonTrig ? triggerMuonAcceptance(pTD2_gen, etaD2_gen) : muonAcceptance(pTD2_gen, etaD2_gen));
      if(!mu1InAccep_gen || !mu2InAccep_gen) continue;

      // Derive the acceptance correction
      double accCorr = (corrAcc ? 1.0 : 1.0); // Set to 1 for now
      
      // Fill the passing histogram
      hPassed.Fill(pT_gen, accCorr);
    }
      
    // Set the histograms to the acceptance
    acc.reset(new TEfficiency("acc", "acc", size, pTBin));
    acc->SetTotalHistogram(hTotal, "f");
    acc->SetPassedHistogram(hPassed, "f");
    if(TEfficiency::CheckWeights(hPassed, hTotal)) { acc->SetStatisticOption(TEfficiency::kBJeffrey); }
  }

  c.cd();
  acc->Draw();
  TFile canvasFile(outPath1, "RECREATE");
  canvasFile.cd();
  c.Write("canvas");
  canvasFile.Write();
  canvasFile.Close();
  //c.Print(outputFile);

  TFile outFile(outPath2, "RECREATE");
  outFile.cd();
  acc->Write("acc");
  outFile.Write();
  outFile.Close();

  return;
}
void acc(){
  auto inputFile = "/storage1/users/wl33/DiMuTrees/pPb2016/Tree/VertexCompositeTree_JPsiToMuMu_pPb-Bst_GENonly_pPb816Summer16_DiMuGENONLY.root";
  double pTBin[] = {6.5, 9.0, 50.0};
  int size = 2;
  bool muonTrig = true;
  bool corrAcc = true;
  auto title = "PromptWTJPsi0_1.4";
  double y1 = 0.0;
  double y2 = 1.4;
  auto outPath1 = "./rootFile/PromptWT/JPsi/0_1.4/PromptWTJPsi0_1.4_acc_canvas.root";
  auto outPath2 = "./rootFile/PromptWT/JPsi/0_1.4/PromptWTJPsi0_1.4_acc.root";
  GetAcc(inputFile, pTBin, size, corrAcc, muonTrig, y1, y2, title, outPath1, outPath2);
}

