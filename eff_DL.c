#include "Utilities/Ntuple/VertexCompositeTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TEfficiency.h"
#include <iostream>

void eff_DL()
{
  const auto& inputFile = "/storage1/users/wl33/DiMuTrees/pPb2016/Tree/VertexCompositeTree_JPsiToMuMu_pPb-Bst_pPb816Summer16_DiMuMC.root";
  const auto& treeDir = "dimucontana_mc"; // For MC use dimucontana_mc

  VertexCompositeTree tree;
  if (!tree.GetTree(inputFile, treeDir)) { std::cout << "Invalid tree!" << std::endl; return 0.0; }

  double pTBin[] = {6.5,9.0,50.0};
  TEfficiency *eff = new TEfficiency("eff1", "0_1.4(PromptJPsi)", 2, pTBin);
  TCanvas *c = new TCanvas();

  long int nevents = tree.GetEntries();
  for(Long64_t jentry=0; jentry<nevents; jentry++)
  { 
    if(!(jentry % 1000000)) std::cout<<"Processed "<<jentry<<" events out of "<<nevents<<std::endl;
    
    if (tree.GetEntry(jentry)<0) { std::cout << "Invalid entry!" << std::endl; return 0.0; }
    
    for(uint iReco=0; iReco<tree.candSize(); iReco++){
      const auto& pT = tree.pT()[iReco];
      const auto& p = pT*std::cosh(tree.eta()[iReco]);
      const auto& decayLen = (tree.V3DDecayLength()[iReco]*tree.V3DCosPointingAngle()[iReco])*(3.0969/p)*10;

      const auto& pD1 = tree.pTD1()[iReco]*std::cosh(tree.EtaD1()[iReco]);
      const auto& pD2 = tree.pTD2()[iReco]*std::cosh(tree.EtaD2()[iReco]);
      const bool inMuonAcceptance = (pD1 > 3.5 && fabs(tree.EtaD1()[iReco])<2.4) && (pD2 > 3.5 && fabs(tree.EtaD2()[iReco])<2.4);
      if (!inMuonAcceptance) continue;

      const bool softCand = tree.softCand(iReco);
      const bool goodEvt = tree.evtSel()[0];
      const bool trigCand = tree.trigHLT()[0] && tree.trigCand(0, iReco);
      if (!softCand || !goodEvt) continue;
      
      const bool yRange = (abs(tree.y()[iReco])>0 && abs(tree.y()[iReco])<=1.4);
      if (!yRange) continue;
      
      bool Passed = decayLen < (-0.04 + 0.1812/TMath::Power(pT,0.3096));
      eff->Fill(Passed, pT);
    }    
  }
  eff->Draw();
  c->Print("0_1.4_PromptJPsi.png");

  TFile outFile("0_1.4_PromptJPsi.root", "RECREATE");
  outFile.cd();
  eff->Write("JPsi_Efficiency");
  outFile.Write();
  outFile.Close();

}

