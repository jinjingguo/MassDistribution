#include "Utilities/Ntuple/VertexCompositeTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TEfficiency.h"
#include <iostream>

double GetThre(double lowpt, double highpt, double lowy, double highy)
{
  const auto& inputFile = "/storage1/users/wl33/DiMuTrees/pPb2016/Tree/VertexCompositeTree_JPsiToMuMu_pPb-Bst_pPb816Summer16_DiMuMC.root";
  const auto& treeDir = "dimucontana_mc"; // For MC use dimucontana_mc

  VertexCompositeTree tree;
  if (!tree.GetTree(inputFile, treeDir)) { std::cout << "Invalid tree!" << std::endl; return 0.0; }

  TH1F *hist1 = new TH1F("hist1", "DecayLength", 1000, -2, 5);

  long int nevents = tree.GetEntries();
  for(Long64_t jentry=0; jentry<nevents; jentry++)
  { 
    if(!(jentry % 1000000)) std::cout<<"Processed "<<jentry<<" events out of "<<nevents<<std::endl;
    
    if (tree.GetEntry(jentry)<0) { std::cout << "Invalid entry!" << std::endl; return 0.0; }
    
    for(uint iGen=0; iGen<tree.candSize_gen(); iGen++){
      const auto& pT = tree.pT()[iGen];
      const auto& p = pT*std::cosh(tree.eta()[iGen]);
      const auto& decayLen = (tree.V3DDecayLength()[iGen]*tree.V3DCosPointingAngle()[iGen])*(3.0969/p)*10;
      
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
        const bool trigCand = tree.trigHLT()[0] && tree.trigCand(0, recoIdx);
        genPassed = genPassed && (softCand && goodEvt && trigCand);
      }

      if(genPassed){
        const bool pT = (tree.pT()[iGen]>lowpt && tree.pT()[iGen]<highpt);
        const bool y = (abs(tree.y()[iGen])>lowy && abs(tree.y()[iGen])<highy);
        genPassed = genPassed && (pT && y);
      }
      
      if(genPassed){
        hist1->Fill(decayLen);
      }
    }
  }

  Double_t norm = hist1->GetEntries();
  hist1->Scale(1/norm);

  Double_t x[1000], y[1000];
  Int_t n = 1000;
  Double_t thre = 0;
  for (Int_t i=0;i<n;i++) {
    x[i] = -2 + 0.007*i;
    y[i] = hist1->Integral(0,i);
    if(y[i]>=0.9){
      thre = x[i];
      break;
    }
  }

  return thre;
}

void decayLen(){
  Double_t x[14] = {6.5,7,8,9,10,11,12,13,15,17,19,21,25,30};
  Double_t y[13];
  Int_t n = 13;
  for(Int_t i=0; i<n; i++){
    y[i] = GetThre(x[i], x[i+1], 0, 1.4);
    std::cout<<"point: "<<i<<std::endl;
  }
  
  TCanvas *canvas1 = new TCanvas("canvas1", "canvas1");
  TGraph *graph1 = new TGraph(n,x,y);
  graph1->SetMarkerSize(1);
  graph1->SetMarkerStyle(20);
  graph1->Draw("PAC");
  graph1->GetXaxis()->SetTitle("pT");
  graph1->GetYaxis()->SetTitle("decayLen_90%_threshold");
  graph1->SetTitle("0_1.4");

  TF1 *func = new TF1("function", "[0]+[1]/x");
  func->SetLineColor(kRed);
  graph1->Fit(func);
  func->Draw("same");


  canvas1->Print("pT_dependence.png");
}
