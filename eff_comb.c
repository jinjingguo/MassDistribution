#include "Utilities/Ntuple/VertexCompositeTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TEfficiency.h"
#include <iostream>
#include <string>

bool triggerMuonAcceptance(Double_t pt, Double_t eta)
{
    return ( fabs(eta) < 2.4 &&
            (    ( fabs(eta) < 1.2 && pt >= 3.3 ) ||
             (  1.2 <= fabs(eta) && fabs(eta) < 2.1 && pt >= 3.93-1.11*fabs(eta)) ||
             (  2.1 <= fabs(eta) && fabs(eta) < 2.4 && pt >= 1.3)
             )
            );
}

bool muonAcceptance(Double_t pt, Double_t eta)
{
    return ( fabs(eta) < 2.4 &&
            (    ( fabs(eta) < 0.8 && pt >= 3.3 ) ||
             (  0.8 <= fabs(eta) && fabs(eta) < 1.5 && pt >= 5.81-3.14*fabs(eta)) ||
             (  1.5 <= fabs(eta) && fabs(eta) < 2.4 && pt >= 0.8 && 1.89-0.526*fabs(eta) )
             )
            );
}


void GetEff(std::string inputFile1, std::string inputFile2, double pTBin[], int size, bool muonTrig, double y1, double y2, const char* title, const char* outPath1, const char* outPath2, const char* outPath3, const char* outPath4)
{
  //const auto& inputFile1 = "/storage1/users/wl33/DiMuTrees/pPb2016/Tree/VertexCompositeTree_JPsiToMuMu_pPb-Bst_pPb816Summer16_DiMuMC.root";
  //const auto& inputFile2 = "/storage1/users/wl33/DiMuTrees/pPb2016/Tree/VertexCompositeTree_JPsiToMuMu_PbP-Bst_pPb816Summer16_DiMuMC.root";
  const auto& treeDir = "dimucontana_mc"; // For MC use dimucontana_mc

  VertexCompositeTree tree1;
  VertexCompositeTree tree2;
  if (!tree1.GetTree(inputFile1, treeDir)||!tree2.GetTree(inputFile2, treeDir)) { std::cout << "Invalid tree!" << std::endl; }

  //double pTBin[] = {0,3,6.5,9.0,50.0};
  TEfficiency *eff_pPb = new TEfficiency("eff1", "pPb", size, pTBin);
  TCanvas *c = new TCanvas();
  c->SetTitle(title);
  c->Divide(3,2); 

  long int nevents = tree1.GetEntries();
  for(Long64_t jentry=0; jentry<nevents; jentry++)
  { 
    if(!(jentry % 1000000)) std::cout<<"Processed "<<jentry<<" events out of "<<nevents<<std::endl;
    
    if (tree1.GetEntry(jentry)<0) { std::cout << "Invalid entry!" << std::endl;}
    
    for(uint iReco=0; iReco<tree1.candSize(); iReco++){
      const auto& pT = tree1.pT()[iReco];
      const auto& p = pT*std::cosh(tree1.eta()[iReco]);
      const auto& decayLen = (tree1.V3DDecayLength()[iReco]*tree1.V3DCosPointingAngle()[iReco])*(3.0969/p)*10;

      if(muonTrig){
        if ( !triggerMuonAcceptance(tree1.pTD1_gen()[iReco] , tree1.EtaD1_gen()[iReco]) || !triggerMuonAcceptance(tree1.pTD2_gen()[iGenCand] , tree1.EtaD2_gen()[iReco]))
            continue;
      }
      else{
        if ( !muonAcceptance(tree1.pTD1_gen()[iReco] , tree1.EtaD1_gen()[iReco]) || !muonAcceptance(tree1.pTD2_gen()[iReco] , tree1.EtaD2_gen()[iReco]) ) continue;
      }

      const bool yRange = (abs(tree1.y()[iReco])>y1 && abs(tree1.y()[iReco])<=y2);
      if (!yRange) continue;
        
      const bool softCand = tree1.softCand(iReco);
      const bool goodEvt = tree1.evtSel()[0];
      const bool trigCand = tree1.trigHLT()[0] && tree1.trigCand(0, iReco);
      
      bool Passed = softCand && goodEvt && trigCand;
      eff_pPb->Fill(Passed, pT);
    }    
  }

  TEfficiency *eff_Pbp = new TEfficiency("eff2", "Pbp", size, pTBin);

  long int nevents2 = tree2.GetEntries();
  for(Long64_t jentry=0; jentry<nevents2; jentry++)
  {
    if(!(jentry % 1000000)) std::cout<<"Processed "<<jentry<<" events out of "<<nevents2<<std::endl;

    if (tree2.GetEntry(jentry)<0) { std::cout << "Invalid entry!" << std::endl;}

    for(uint iReco=0; iReco<tree2.candSize(); iReco++){
      const auto& pT = tree2.pT()[iReco];
      const auto& p = pT*std::cosh(-tree2.eta()[iReco]);
      const auto& decayLen = (tree2.V3DDecayLength()[iReco]*tree2.V3DCosPointingAngle()[iReco])*(3.0969/p)*10;

      if(muonTrig){
        if ( !triggerMuonAcceptance(tree2.pTD1_gen()[iReco] , tree2.EtaD1_gen()[iReco]) || !triggerMuonAcceptance(tree2.pTD2_gen()[iReco] , tree2.EtaD2_gen()[iReco]))
                continue;
      }
      else{
        if ( !muonAcceptance(tree2.pTD1_gen()[iReco] , tree2.EtaD1_gen()[iReco]) || !muonAcceptance(tree2.pTD2_gen()[iReco] , tree2.EtaD2_gen()[iReco]) ) continue;
      }
        
      const bool yRange = (abs(tree2.y()[iReco])>y1 && abs(tree2.y()[iReco])<=y2);
      if (!yRange) continue;

      const bool softCand = tree2.softCand(iReco);
      const bool goodEvt = tree2.evtSel()[0];
      const bool trigCand = tree2.trigHLT()[0] && tree2.trigCand(0, iReco);

      bool Passed = softCand && goodEvt && trigCand;
      eff_Pbp->Fill(Passed, pT);
    }
  }

  auto hPassed = *((TH1D*)eff_pPb->GetPassedHistogram());
  hPassed.Add(eff_pPb->GetPassedHistogram(), eff_Pbp->GetPassedHistogram(), 110.8/nevents, 62.6/nevents2);

  auto hTotal = *((TH1D*)eff_pPb->GetTotalHistogram());
  hTotal.Add(eff_pPb->GetTotalHistogram(), eff_Pbp->GetTotalHistogram(), 110.8/nevents, 62.6/nevents2);

  auto hDivide1 = *((TH1D*)eff_pPb->GetPassedHistogram());
  hDivide1 = hDivide1*hTotal;
  auto hTotal1 = *((TH1D*)eff_pPb->GetTotalHistogram());
  hTotal1 = hTotal1*hPassed;
  auto hDivide2 = *((TH1D*)eff_Pbp->GetPassedHistogram());
  hDivide2 = hDivide2*hTotal;
  auto hTotal2 = *((TH1D*)eff_Pbp->GetTotalHistogram());
  hTotal2 = hTotal2*hPassed;

  TEfficiency eff_comb;
  eff_comb.SetPassedHistogram(hPassed, "f");
  eff_comb.SetTotalHistogram(hTotal, "f");
  eff_comb.SetTitle("comb");


  c->cd(1);
  eff_pPb->Draw();
  c->cd(2);
  eff_Pbp->Draw();
  c->cd(3);
  eff_comb.Draw();
  TFile canvasFile(outPath1, "RECREATE");
  canvasFile.cd();
  c->Write("canvas");
  canvasFile.Write();
  canvasFile.Close();
  //c->Print(outputFile);

  TFile outFile1(outPath2, "RECREATE");
  outFile1.cd();
  eff_pPb->Write("eff_pPb");
  outFile1.Write();
  outFile1.Close();
  TFile outFile2(outPath3, "RECREATE");
  outFile2.cd();
  eff_Pbp->Write("eff_Pbp");
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
  auto title = "PromptWTJPsi0_1.4";
  double y1 = 0.0;
  double y2 = 1.4;
  auto outPath1 = "./rootFile/PromptWT/JPsi/0_1.4/PromptWTJPsi0_1.4_canvas.root";
  auto outPath2 = "./rootFile/PromptWT/JPsi/0_1.4/PromptWTJPsi0_1.4_eff_pPb.root";
  auto outPath3 = "./rootFile/PromptWT/JPsi/0_1.4/PromptWTJPsi0_1.4_eff_Pbp.root";
  auto outPath4 = "./rootFile/PromptWT/JPsi/0_1.4/PromptWTJPsi0_1.4_eff_comb.root";
  GetEff(inputFile1, inputFile2, pTBin, size, muonTrig, y1, y2, title, outPath1, outPath2, outPath3, outPath4);
}
