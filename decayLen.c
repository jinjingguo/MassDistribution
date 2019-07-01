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

  TH1F *hist1 = new TH1F("hist1", "DecayLength", 2000, -2, 5);

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
      
      const bool softCand = tree.softCand(recoIdx);
      const bool goodEvt = tree.evtSel()[0];
      const bool trigCand = tree.trigHLT()[0] && tree.trigCand(0, recoIdx);
      if (!softCand || !goodEvt || !trigCand) continue;

      const bool pTRange = (tree.pT()[iReco]>lowpt && tree.pT()[iReco]<=highpt);
      const bool yRange = (abs(tree.y()[iReco])>lowy && abs(tree.y()[iReco])<=highy);
      if (!pTRange || !yRange) continue;
      
      hist1->Fill(decayLen);
    }
  }

  Double_t norm = hist1->GetEntries();
  hist1->Scale(1/norm);

  Double_t x[2000], y[2000];
  Int_t n = 2000;
  Double_t thre = 0;
  for (Int_t i=1;i<n;i++) {
    x[i] = -2 + 0.0035*i;
    y[i] = hist1->Integral(1,i);
    if(y[i]>=0.9){
      thre = x[i];
      break;
    }
  }

  return thre;
}

void decayLen(){
  Double_t xBin[16] = {6.5,7,8,9,10,11,12,13,15,17,19,21,25,30,35,40};
  Double_t x[15], y[15];
  Int_t n = 15;
  for(Int_t i=0; i<n; i++){
    x[i] = (xBin[i+1]+xBin[i])/2.0;
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
