#include "Utilities/Ntuple/VertexCompositeTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <iostream>
  
void DM_1()
{
  const auto& inputFile = "/storage1/users/wl33/DiMuTrees/pPb2016/Tree/VertexCompositeTree_PADoubleMuon_PARun2016C_DiMuMassMin2.root";
  const auto& treeDir = "dimucontana"; // For MC use dimucontana_mc
  
  // Extract the tree
  VertexCompositeTree tree;
  if (!tree.GetTree(inputFile, treeDir)) { std::cout << "Invalid tree!" << std::endl; return; }
  
  TH1F *hmass1 = new TH1F("hmass1", "0.2 < pT < 1.8 GeV/c", 100, 2.2, 4.2);
  hmass1 -> GetXaxis() -> SetTitle("Invariant Mass");
  hmass1 -> GetYaxis() -> SetTitle("Counts");
  TH1F *hmass2 = new TH1F("hmass2", "1.8 < pT < 3.0 GeV/c", 100, 2.2, 4.2);
  hmass2 -> GetXaxis() -> SetTitle("Invariant Mass");
  hmass2 -> GetYaxis() -> SetTitle("Counts");
    TH1F *hmass3 = new TH1F("hmass3", "3.0 < pT < 4.5 GeV/c", 100, 2.2, 4.2);
    hmass3 -> GetXaxis() -> SetTitle("Invariant Mass");
    hmass3 -> GetYaxis() -> SetTitle("Counts");
    TH1F *hmass4 = new TH1F("hmass4", "4.5 < pT < 6.0 GeV/c", 100, 2.2, 4.2);
    hmass4 -> GetXaxis() -> SetTitle("Invariant Mass");
    hmass4 -> GetYaxis() -> SetTitle("Counts");
    TH1F *hmass5 = new TH1F("hmass5", "6.0 < pT < 8.0 GeV/c", 100, 2.2, 4.2);
    hmass5 -> GetXaxis() -> SetTitle("Invariant Mass");
    hmass5 -> GetYaxis() -> SetTitle("Counts");
    TH1F *hmass6 = new TH1F("hmass6", "8.0 < pT < 10.0 GeV/c", 100, 2.2, 4.2);
    hmass6 -> GetXaxis() -> SetTitle("Invariant Mass");
    hmass6 -> GetYaxis() -> SetTitle("Counts");
    TH1F *hmass7 = new TH1F("hmass7", "10.0 < pT < 20.0 GeV/c", 100, 2.2, 4.2);
    hmass7 -> GetXaxis() -> SetTitle("Invariant Mass");
    hmass7 -> GetYaxis() -> SetTitle("Counts");
  
    TCanvas *canvas = new TCanvas("canvas");
    canvas -> Divide(3,3);
  // Loop over the tree
  //std::cout << tree.GetEntries() << std::endl;
  long int nevents = tree.GetEntries();
  for(Long64_t jentry=0; jentry<nevents; jentry++)
  {
    if(!(jentry % 10000)) std::cout<<"Processed "<<jentry<<" events out of "<<nevents<<std::endl;
    // Get the entry
    if (tree.GetEntry(jentry)<0) { std::cout << "Invalid entry!" << std::endl; return; }
    // Loop over the candidates
    for(uint iCand=0; iCand<tree.candSize(); iCand++)
    {
        if (tree.trigHLT()[0]==true&&tree.evtSel()[0]==true&&tree.softMuon1()[iCand]==true&&tree.softMuon2()[iCand]==true&&tree.Ntrkoffline()>185&&tree.Ntrkoffline()<250&&tree.mass()[iCand]>2.2&&tree.mass()[iCand]<4.2&&abs(tree.y()[iCand])<1.4)
            {
                if (tree.pT()[iCand]>0.2&&tree.pT()[iCand]<1.8) {
                    hmass1 -> Fill(tree.mass()[iCand]);
                    canvas -> cd(1);
                    hmass1 -> Draw();
                }
                else if (tree.pT()[iCand]>1.8&&tree.pT()[iCand]<3) {
                    hmass2 -> Fill(tree.mass()[iCand]);
                    canvas -> cd(2);
                    hmass2 -> Draw();
                }
                else if (tree.pT()[iCand]>3&&tree.pT()[iCand]<4.5) {
                    hmass3 -> Fill(tree.mass()[iCand]);
                    canvas -> cd(3);
                    hmass3 -> Draw();
                }
                else if (tree.pT()[iCand]>4.5&&tree.pT()[iCand]<6) {
                    hmass4 -> Fill(tree.mass()[iCand]);
                    canvas -> cd(4);
                    hmass4 -> Draw();
                }
                else if (tree.pT()[iCand]>6&&tree.pT()[iCand]<8) {
                    hmass5 -> Fill(tree.mass()[iCand]);
                    canvas -> cd(5);
                    hmass5 -> Draw();
                }
                else if (tree.pT()[iCand]>8&&tree.pT()[iCand]<10) {
                    hmass6 -> Fill(tree.mass()[iCand]);
                    canvas -> cd(6);
                    hmass6 -> Draw();
                }
                else if (tree.pT()[iCand]>10&&tree.pT()[iCand]<20) {
                    hmass7 -> Fill(tree.mass()[iCand]);
                    canvas -> cd(7);
                    hmass7 -> Draw();
                }
            }
        
    }
    //if (jentry > 10000000) break;
  }
    canvas -> Print("185_250.pdf");
}

