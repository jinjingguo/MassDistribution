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


void DrawDLEff(int id, bool muonTrig, double pTBin[], int n, const char* title, double lowy, double highy, double a, double b, double c, const char* outfile1, const char* outfile2)
{
    std::map<std::string , std::string> inputFile;
    inputFile["pPb"] = "/storage1/users/wl33/DiMuTrees/pPb2016/Tree/VertexCompositeTree_JPsiToMuMu_pPb-Bst_pPb816Summer16_DiMuMC.root";
    inputFile["Pbp"] = "/storage1/users/wl33/DiMuTrees/pPb2016/Tree/VertexCompositeTree_JPsiToMuMu_PbP-Bst_pPb816Summer16_DiMuMC.root";
    const auto& treeDir = "dimucontana_mc"; // For MC use dimucontana_mc
    
    //double pTBin[] = {6.5,9.0,50.0};
    TEfficiency *eff = new TEfficiency("eff1", title, n, pTBin);
    TCanvas *canvas = new TCanvas();
    
    for(const auto& s : inputFile)
    {
        VertexCompositeTree tree;
        if (!tree.GetTree(s.second, treeDir)) { std::cout << "Invalid tree for: " << s.second << "!" << std::endl; return; }
        
        const auto& nevents = tree.GetEntries();
        for(Long64_t jentry=0; jentry<nevents; jentry++)
        {
            if(!(jentry % 1000000)) std::cout<<"Processed "<<jentry<<" events out of "<<nevents<<std::endl;
            
            if (tree.GetEntry(jentry)<0) { std::cout << "Invalid entry for: " << s.second << "!" << std::endl; return; }
            
            for(uint iReco=0; iReco<tree.candSize(); iReco++)
            {
                const auto& rap = (s.first=="pPb" ? 1.0 : -1.0)*tree.y()[iReco];
                const bool yRange = (fabs(rap)>lowy && fabs(rap)<=highy);
                if(!yRange) continue;
                if (fabs(tree.PID_gen()[iReco])!=id) continue;
                
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
                
                bool Passed = decayLen < (a + b/TMath::Power(pT,c));
                eff->Fill(Passed, pT);
            }
        }
  }
  eff->Draw();
  canvas->Print(outfile1);

  TFile outFile(outfile2, "RECREATE");
  outFile.cd();
  eff->Write("JPsi_Efficiency");
  outFile.Write();
  outFile.Close();
}

void eff_DL2(){
    int id = 443;
    bool muonTrig = false;
    double pTBin[] = {6.5, 9.0, 50.0};
    int n = 2;
    const char* title = "0_1.4_0.9NT_PromptJPsi";
    double lowy = 0.0;
    double highy = 1.4;
    double a = -0.0100138;
    double b = 0.122266;
    double c = 0.430967;
    const char* outfile1 = "0_1.4_0.9NT_PromptJPsi.png";
    const char* outfile2 = "0_1.4_0.9NT_PromptJPsi.root";
    DrawDLEff(id, muonTrig, pTBin, n, title, lowy, highy, a, b, c, outfile1, outfile2);
}
