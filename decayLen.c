#include "Utilities/Ntuple/VertexCompositeTree.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TEfficiency.h"
#include <iostream>

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

double GetThre(double lowpt, double highpt, double lowy, double highy, bool muonTrig, double thr)
{
  std::map<std::string , std::string> inputFile;
  inputFile["pPb"] = "/storage1/users/wl33/DiMuTrees/pPb2016/Tree/VertexCompositeTree_JPsiToMuMu_pPb-Bst_pPb816Summer16_DiMuMC.root";
  inputFile["Pbp"] = "/storage1/users/wl33/DiMuTrees/pPb2016/Tree/VertexCompositeTree_JPsiToMuMu_PbP-Bst_pPb816Summer16_DiMuMC.root";
  const auto& treeDir = "dimucontana_mc"; // For MC use dimucontana_mc

  TH1D hist("hist", "DecayLength", 6000, -4.0, 8.0);
  hist.Sumw2();
  for(const auto& s : inputFile)
  {
    VertexCompositeTree tree;
    if (!tree.GetTree(s.second, treeDir)) { std::cout << "Invalid tree for: " << s.second << "!" << std::endl; return 0.0; }

    const auto& nevents = tree.GetEntries();
    for(Long64_t jentry=0; jentry<nevents; jentry++)
    { 
      if(!(jentry % 1000000)) std::cout<<"Processed "<<jentry<<" events out of "<<nevents<<std::endl;
    
      if (tree.GetEntry(jentry)<0) { std::cout << "Invalid entry for: " << s.second << "!" << std::endl; return 0.0; }
    
      for(uint iReco=0; iReco<tree.candSize(); iReco++)
      {
	const auto& rap = (s.first=="pPb" ? 1.0 : -1.0)*tree.y()[iReco];
	const bool yRange = (fabs(rap)>lowy && fabs(rap)<=highy);
  	const bool pTRange = (tree.pT()[iReco]>lowpt && tree.pT()[iReco]<=highpt);
	if(!pTRange || !yRange) continue;
        
        const auto& pTD1 = tree.pTD1()[iReco];
        const auto& etaD1 = (s.first=="pPb" ? 1.0 : -1.0)*tree.EtaD1()[iReco];
        const bool mu1InAccep = (muonTrig ? triggerMuonAcceptance(pTD1, etaD1) : muonAcceptance(pTD1, etaD1));
        const auto& pTD2 = tree.pTD2()[iReco];
        const auto& etaD2 = (s.first=="pPb" ? 1.0 : -1.0)*tree.EtaD2()[iReco];
        const bool mu2InAccep = (muonTrig ? triggerMuonAcceptance(pTD2, etaD2) : muonAcceptance(pTD2, etaD2));
        if(!mu1InAccep || !mu2InAccep) continue;

        const bool softCand = tree.softCand(iReco);
        const bool goodEvt = tree.evtSel()[0];
        const bool trigCand = (muonTrig ? (tree.trigHLT()[0] && tree.trigCand(0, iReco)) : true);
        if(!softCand || !goodEvt || !trigCand) continue;

	const auto& pT = tree.pT()[iReco];
	const auto& eta = (s.first=="pPb" ? 1.0 : -1.0)*tree.eta()[iReco];
        const auto& p = pT*std::cosh(eta);
        const auto& decayLen = (tree.V3DDecayLength()[iReco]*tree.V3DCosPointingAngle()[iReco])*(3.0969/p)*10.0;

 	const auto& weight = (s.first=="pPb" ? 110.78 : 62.64)/(110.78 + 62.64);
        hist.Fill(decayLen, weight);
      }
    }
  }

  const auto& norm = hist.Integral();
  hist.Scale(1.0/norm);
  double thre = -999.0;
  for(int i=0; i<hist.GetXaxis()->GetNbins(); i++)
  {
    const auto& x = hist.GetBinCenter(i);
    const auto& y = hist.Integral(0, i);
    if(y>=thr)
    {
      thre = x;
      break;
    }
  }

  return thre;
}

void decayLen(){
  Double_t xBin[16] = {6.5,7,8,9,10,11,12,13,15,17,19,21,25,30,35,40};
  Double_t x[15], y[15];
  Int_t n = 15;
  for(Int_t i=0; i<15; i++){
    if(i==15){break;}
    x[i] = (xBin[i+1]+xBin[i])/2.0;
    y[i] = GetThre(xBin[i], xBin[i+1], 0, 1.4, false, 0.9);
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

  TF1 *func = new TF1("function", "[0]+[1]/(x**[2])");
  func->SetLineColor(kRed);
  graph1->Fit(func);
  func->Draw("same");


  canvas1->Print("0_1.4_0.9NT.png");
}
