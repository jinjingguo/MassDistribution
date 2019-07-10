#include "Utilities/Ntuple/VertexCompositeTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <iostream>
  
void DM_Ntrk_2()
{
  const auto& inputFile = "/storage1/users/ags9/pPb2016/Charmonia/VertexCompositeTree_PADoubleMuon_PARun2016C_Charmonia.root";
  const auto& treeDir = "dimucontana"; // For MC use dimucontana_mc
  
  // Extract the tree
  VertexCompositeTree tree;
  if (!tree.GetTree(inputFile, treeDir)) { std::cout << "Invalid tree!" << std::endl; return; }
  
  TH1F *hmass1 = new TH1F("hmass1", "0 < Ntrk < 35 ", 200, 2.2, 4.2);
  hmass1 -> GetXaxis() -> SetTitle("Invariant Mass");
  hmass1 -> GetYaxis() -> SetTitle("Counts");
  TH1F *hmass2 = new TH1F("hmass2", "35 < Ntrk < 60", 200, 2.2, 4.2);
  hmass2 -> GetXaxis() -> SetTitle("Invariant Mass");
  hmass2 -> GetYaxis() -> SetTitle("Counts");
    TH1F *hmass3 = new TH1F("hmass3", "60 < Ntrk < 90 ", 200, 2.2, 4.2);
    hmass3 -> GetXaxis() -> SetTitle("Invariant Mass");
    hmass3 -> GetYaxis() -> SetTitle("Counts");
    TH1F *hmass4 = new TH1F("hmass4", "90 < Ntrk < 120", 200, 2.2, 4.2);
    hmass4 -> GetXaxis() -> SetTitle("Invariant Mass");
    hmass4 -> GetYaxis() -> SetTitle("Counts");
    TH1F *hmass5 = new TH1F("hmass5", "120 < Ntrk < 150", 200, 2.2, 4.2);
    hmass5 -> GetXaxis() -> SetTitle("Invariant Mass");
    hmass5 -> GetYaxis() -> SetTitle("Counts");
    TH1F *hmass6 = new TH1F("hmass6", "150 < Ntrk < 185", 200, 2.2, 4.2);
    hmass6 -> GetXaxis() -> SetTitle("Invariant Mass");
    hmass6 -> GetYaxis() -> SetTitle("Counts");
    TH1F *hmass7 = new TH1F("hmass7", "185 < Ntrk < 250", 200, 2.2, 4.2);
    hmass7 -> GetXaxis() -> SetTitle("Invariant Mass");
    hmass7 -> GetYaxis() -> SetTitle("Counts");
    TH1F *hmass8 = new TH1F("hmass8", "Ntrk > 250", 200, 2.2, 4.2);
    hmass8 -> GetXaxis() -> SetTitle("Invariant Mass");
    hmass8 -> GetYaxis() -> SetTitle("Counts"); 
 
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
        if (tree.trigHLT()[0]==true&&tree.PFMuon1()[iCand]==true&&tree.PFMuon2()[iCand]==true&&abs(tree.y()[iCand])<2.4&&tree.pT()[iCand]>1.8&&tree.pT()[iCand]<3)
            {
                if (tree.Ntrkoffline()>0&&tree.Ntrkoffline()<35) {
                    hmass1 -> Fill(tree.mass()[iCand]);
                    canvas -> cd(1);
                    hmass1 -> Draw();
                    gPad -> SetLogy();
                }
                else if (tree.Ntrkoffline()>35&&tree.Ntrkoffline()<60) {
                    hmass2 -> Fill(tree.mass()[iCand]);
                    canvas -> cd(2);
                    hmass2 -> Draw();
                    gPad -> SetLogy();
                }
                else if (tree.Ntrkoffline()>60&&tree.Ntrkoffline()<90) {
                    hmass3 -> Fill(tree.mass()[iCand]);
                    canvas -> cd(3);
                    hmass3 -> Draw();
                    gPad -> SetLogy();
                }
                else if (tree.Ntrkoffline()>90&&tree.Ntrkoffline()<120) {
                    hmass4 -> Fill(tree.mass()[iCand]);
                    canvas -> cd(4);
                    hmass4 -> Draw();
                    gPad -> SetLogy();
                }
                else if (tree.Ntrkoffline()>120&&tree.Ntrkoffline()<150) {
                    hmass5 -> Fill(tree.mass()[iCand]);
                    canvas -> cd(5);
                    hmass5 -> Draw();
                    gPad -> SetLogy();
                }
                else if (tree.Ntrkoffline()>150&&tree.Ntrkoffline()<185) {
                    hmass6 -> Fill(tree.mass()[iCand]);
                    canvas -> cd(6);
                    hmass6 -> Draw();
                    gPad -> SetLogy();
                }
                else if (tree.Ntrkoffline()>185&&tree.Ntrkoffline()<250) {
                    hmass7 -> Fill(tree.mass()[iCand]);
                    canvas -> cd(7);
                    hmass7 -> Draw();
                    gPad -> SetLogy();
                }
                else if (tree.Ntrkoffline()>250) {
                    hmass8 -> Fill(tree.mass()[iCand]);
                    canvas -> cd(8);
                    hmass8 -> Draw();
                    gPad -> SetLogy();
                }

            }
        
    }
    //if (jentry > 10000000) break;
  }
    canvas -> Print("1.8_3_2.pdf");
}

