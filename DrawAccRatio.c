#include "Utilities/Ntuple/VertexCompositeTree.h"
#include "TF1.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TEfficiency.h"
#include "TFile.h"
#include <iostream>

TGraph* DivideEff(TEfficiency* eff1, TEfficiency* eff2, const char* title, int n, double xBin[]){
  double x[n], y[n];

  for (int i=0; i<n; i++){
    double num1 = eff1->GetEfficiency(i+1);
    double num2 = eff2->GetEfficiency(i+1);
    y[i] = num1/num2;
    x[i] = (xBin[i+1]+xBin[i])/2.0;
  }  

  TGraph* graph = new TGraph(n, x, y);
  graph->SetMarkerSize(1);
  graph->SetMarkerStyle(20);
  //graph->Draw("PAC");
  graph->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  graph->GetYaxis()->SetTitle("Psi2S to JPsi ratio");
  graph->SetTitle(title);
  return graph;
}

void getRatio(const char* inputFile1, const char* inputFile2, const char* inputFile3, const char* inputFile4, const char* title, const char* outputFile, int n1, double x1Bin[], int n2, double x2Bin[]){
  TFile* file1 = TFile::Open(inputFile1);
  TEfficiency* eff_JPsi1;
  file1->TFile::GetObject("acc", eff_JPsi1);
 
  TFile* file2 = TFile::Open(inputFile2);
  TEfficiency* eff_JPsi2;
  file2->TFile::GetObject("acc", eff_JPsi2);
 
  TFile* file3 = TFile::Open(inputFile3);
  TEfficiency* eff_Psi1;
  file3->TFile::GetObject("acc", eff_Psi1);

  TFile* file4 = TFile::Open(inputFile3);
  TEfficiency* eff_Psi2;
  file3->TFile::GetObject("acc", eff_Psi2);

  TCanvas* canvas = new TCanvas;

  const char* title1 = "0_1.4"; 
  TGraph* graph1;
  graph1 = DivideEff(eff_Psi1, eff_JPsi1, title1, n1, x1Bin);
  //graph1->SetLineColor(kRed);
  //graph1->SetMarkerColor(kRed);

  const char* title2 = "1.4_2.4";
  TGraph* graph2;
  graph2 = DivideEff(eff_Psi2, eff_JPsi2, title2, n2, x2Bin);
  
  canvas->Divide(2,1);
  canvas->cd(1);
  graph1->Draw();
  canvas->cd(2);
  graph2->Draw();
  //gPad->SetLogy();

  canvas->SetTitle(title);
  canvas->Print(outputFile);

  return;
}

void DrawAccRatio(){
  const char* inputFile1 = "./rootFile/NonPromptWT/JPsi/0_1.4/NonPromptWTJPsi0_1.4_acc.root";
  const char* inputFile2 = "./rootFile/NonPromptWT/JPsi/1.4_2.4/NonPromptWTJPsi1.4_2.4_acc.root";
  const char* inputFile3 = "./rootFile/NonPromptWT/Psi2S/0_1.4/NonPromptWTPsi2S0_1.4_acc.root";
  const char* inputFile4 = "./rootFile/NonPromptWT/Psi2S/1.4_2.4/NonPromptWTPsi2S1.4_2.4_acc.root";
  const char* title = "NonPromptWTPsi_over_JPsi_acc";
  const char* outputFile = "NonPromptWTPsi_over_JPsi_acc.png";
  int n1 = 2;
  double x1Bin[]  = {6.5, 9.0, 50.0};
  int n2 = 3;
  double x2Bin[]  = {3.0, 6.5, 9.0, 50.0};
  getRatio(inputFile1, inputFile2, inputFile3, inputFile4, title, outputFile, n1, x1Bin, n2, x2Bin); 
}
