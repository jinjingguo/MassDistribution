#include "Utilities/Ntuple/VertexCompositeTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include <iostream>

void eff()
{
  const auto& inputFile = "/storage1/users/wl33/DiMuTrees/pPb2016/Tree/VertexCompositeTree_JPsiToMuMu_pPb-Bst_pPb816Summer16_DiMuMC.root";
  const auto& treeDir = "dimucontana_mc"; // For MC use dimucontana_mc

  VertexCompositeTree tree;
  if (!tree.GetTree(inputFile, treeDir)) { std::cout << "Invalid tree!" << std::endl; return; }

  double pTBin[] = {0.0,3.0,6.5,9.0,50.0};
  double rapBin[] = {-2.4,-1.4,0,1.4,2.4};

  TEfficiency *eff_NOTrigger = new TEfficiency("eff_NOTrigger", "pT vs rapidity(NT)", 4, rapBin, 4, pTBin);
  TEfficiency *eff_WithTrigger = new TEfficiency("eff_WithTrigger", "pT vs rapidity(WT)", 4, rapBin, 4, pTBin);
  TCanvas *canvas1 = new TCanvas("canvas1", "canvas1");
 
  long int nevents = tree.GetEntries();
  for(Long64_t jentry=0; jentry<nevents; jentry++)
  {
    if(!(jentry % 10000)) std::cout<<"Processed "<<jentry<<" events out of "<<nevents<<std::endl;
    
    if (tree.GetEntry(jentry)<0) { std::cout << "Invalid entry!" << std::endl; return; }  
  
    for(uint iGen=0; iGen<tree.candSize_gen(); iGen++){
      const auto& recoIdx = tree.RecIdx_gen()[iGen];
    
      bool genPassed = true;

      const auto& pD1_gen = tree.pTD1_gen()[iGen]*std::cosh(tree.EtaD1_gen()[iGen]);
      const auto& pD2_gen = tree.pTD2_gen()[iGen]*std::cosh(tree.EtaD2_gen()[iGen]);
      const bool inMuonAcceptance = (pD1_gen > 3.5 && fabs(tree.EtaD1_gen()[iGen])<2.4) && (pD2_gen > 3.5 && fabs(tree.EtaD2_gen()[iGen])<2.4);
      genPassed = genPassed && inMuonAcceptance;
      genPassed = genPassed && (recoIdx>=0);
      
      
      if (genPassed){
        const bool softCand = tree.softCand(recoIdx, "POG");
        const bool goodEvt = tree.evtSel()[0];
        genPassed = genPassed && (softCand && goodEvt);
      } 
      eff_NOTrigger->Fill(genPassed, tree.y_gen()[iGen], tree.pT_gen()[iGen]);
      
      if(genPassed){
        const bool trigCand = tree.trigHLT()[0] && tree.trigCand(0, recoIdx);
        genPassed = genPassed && trigCand;
      }
      eff_WithTrigger->Fill(genPassed, tree.y_gen()[iGen], tree.pT_gen()[iGen]);
    }
  }
 
  eff_NOTrigger->Draw("COLZ");
  TCanvas *canvas2 = new TCanvas("canvas2", "canvas2"); 
  eff_WithTrigger->Draw("COLZ");

  TFile outFile("JPsiEff.root", "RECREATE");
  outFile.cd();
  eff_NOTrigger->Write("JPsi_Efficiency");
  outFile.Write();
  outFile.Close();
}

