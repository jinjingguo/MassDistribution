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

  TEfficiency effxAcc;
  TCanvas *canvas = new TCanvas("canvas");
  
  long int ncands = tree.candSize_gen();
  for(Long64_t iGen=0; iGen<ncands; iGen++)
  {
    if(!(iGen % 10000)) std::cout<<"Processed "<<iGen<<" cands out of "<<ncands<<std::endl;
    
    const auto& recoIdx = tree.RecIdx_gen()[iGen];
    
    bool genPassed = false;
    if (recoIdx>0&&tree.softMuon1()[recoIdx]==true&&tree.softMuon2()[recoIdx]==true&&((tree.pT()[recoIdx]>0&&abs(tree.y()[recoIdx])>1.9)||(tree.pT()[recoIdx]>3&&abs(tree.y()[recoIdx])>1.4)||(tree.pT()[recoIdx]>6.5&&abs(tree.y()[recoIdx])>0))) {
      genPassed = true;
    }

    effxAcc.Fill(genPassed, tree.pT_gen()[iGen], tree.y_gen()[iGen]);


  }
  
  effxAcc.Draw();

}

