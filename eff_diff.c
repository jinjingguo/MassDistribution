#include "Utilities/Ntuple/VertexCompositeTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include <iostream>

void eff_diff()
{
  const auto& inputFile = "/storage1/users/wl33/DiMuTrees/pPb2016/Tree/VertexCompositeTree_JPsiToMuMu_pPb-Bst_pPb816Summer16_DiMuMC.root";
  const auto& treeDir = "dimucontana_mc"; // For MC use dimucontana_mc

  VertexCompositeTree tree;
  if (!tree.GetTree(inputFile, treeDir)) { std::cout << "Invalid tree!" << std::endl; return; }

  TEfficiency *eff_WithTrigger = new TEfficiency("eff", "pT vs rapidity(WT)", 50, -2.4, 2.4, 100, 0, 20);
  TEfficiency *eff_NOTrigger = new TEfficiency("eff", "pT vs rapidity(NT)", 50, -2.4, 2.4, 100, 0, 20);
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
  //canvas1 -> Divide(1,2);
  //canvas1 -> cd(1);
  eff_WithTrigger -> Draw("Colz");
  canvas1->Print("WTrigger.root");
  //canvas1 -> cd(2);
  TCanvas *canvas2 = new TCanvas("canvas2", "NOTrigger");
  eff_NOTrigger -> Draw("Colz");
  canvas2->Print("NTrigger.root");

  TH1* hPassed_WithTrigger =  eff_WithTrigger->GetCopyPassedHisto();
  TH1* hPassed_NOTrigger =  eff_NOTrigger->GetCopyPassedHisto(); 
  TH1* hTotal =  eff_WithTrigger->GetCopyTotalHisto();  

  //TCanvas *canvas2 = new TCanvas("canvas2", "canvas2");
  //canvas2 -> Divide(2,2);
  //canvas2 -> cd(1);
  //hPassed_WithTrigger->Draw("Colz");
  //canvas2 -> cd(2);
  //hPassed_NOTrigger->Draw("Colz");
  //canvas2 -> cd(3);
  //hTotal->Draw("Colz");
 
  //std::cout << "WT_entries: " << hPassed_WithTrigger->GetEntries() << std::endl;
  //std::cout << "NT_entries: " << hPassed_NOTrigger->GetEntries() << std::endl;
  //std::cout << "Tot_entries: " << hTotal->GetEntries() << std::endl;

  TH1* hPassed_Diff = (TH1*) hPassed_NOTrigger->Clone();
  hPassed_Diff->Scale(-93.85/173.42);
  hPassed_Diff->Add(hPassed_WithTrigger);
  hPassed_Diff->Divide(hTotal);
 
  hPassed_Diff->SetTitle("MB_Dimuon_Diff");
  hPassed_Diff->GetXaxis()->SetTitle("rapidity");
  hPassed_Diff->GetYaxis()->SetTitle("pT"); 
 
  const UInt_t Number = 5;
  Double_t Red[Number]    = { 1.00, 0.00, 1.00, 1.00, 0.00};
  Double_t Green[Number]  = { 0.00, 1.00, 1.00, 0.00, 0.00};
  Double_t Blue[Number]   = { 0.00, 0.00, 1.00, 1.00, 1.00};
  Double_t Length[Number] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
  Int_t nb=50;

  //TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
  hPassed_Diff->SetContour(nb);
  hPassed_Diff->SetMaximum(0.5);
  hPassed_Diff->SetMinimum(-0.5);

 
  //hPassed_Diff->Draw("colz");
  //canvas1->Print("HM_Dimuon_Diff2.pdf");
}

