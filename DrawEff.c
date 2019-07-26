#include "Utilities/Ntuple/VertexCompositeTree.h"
#include "TF1.h"
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
  graph->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  //graph->GetYaxis()->SetTitle("ratio");
  graph->SetTitle(title);
  return graph;
}

void getRatio(const char* inputFile1, const char* inputFile2, const char* inputFile3, const char* title, const char* outputFile, int n,double xBin[]){
  TFile* file1 = TFile::Open(inputFile1);
  TEfficiency* eff1;
  file1->TFile::GetObject("eff_pPb", eff1);
  eff1->SetTitle("eff_pPb;p_{T} [GeV/c];Efficiency");

  TFile* file2 = TFile::Open(inputFile2);
  TEfficiency* eff2;
  file2->TFile::GetObject("eff_Pbp", eff2);
  eff2->SetLineColor(kBlue);
  eff2->SetTitle("eff_Pbp;p_{T} [GeV/c];Efficiency");

  TFile* file3 = TFile::Open(inputFile3);
  TEfficiency* eff3;
  file3->TFile::GetObject("eff_comb", eff3);
  eff3->SetLineColor(kRed);
  eff3->SetTitle("eff_comb;p_{T} [GeV/c];Efficiency");

  TCanvas* canvas = new TCanvas;
  canvas->Divide(1,2);
  canvas->cd(1);
  eff1->Draw();
  eff2->Draw("same");
  eff3->Draw("same");
  canvas->cd(1)->BuildLegend(0.7,0.1,0.9,0.3,"y_range");
  eff1->SetTitle(title);

  canvas->cd(2);
  gPad->Divide(2,1);
  const char* title1 = "pPb_over_comb";
  TGraph* graph1;
  graph1 = DivideEff(eff1, eff3, title1, n, xBin);
  //graph->GetYaxis()->SetTitle("ratio");
  canvas->cd(2)->cd(1);
  graph1->Draw("PAC");

  const char* title2 = "Pbp_over_comb";
  TGraph* graph2;
  graph2 = DivideEff(eff2, eff3, title2, n, xBin);
  canvas->cd(2)->cd(2);
  graph2->Draw("PAC");

  canvas->SetTitle(title);
  canvas->Print(outputFile);

  return;
}

void get(const char* inputFile1, const char* inputFile2, const char* inputFile3, const char* title, const char* outputFile, int n,double xBin[]){
  TFile* file1 = TFile::Open(inputFile1);
  TEfficiency* eff1;
  file1->TFile::GetObject("eff_pPb", eff1);
 
  TFile* file2 = TFile::Open(inputFile2);
  TEfficiency* eff2;
  file2->TFile::GetObject("eff_Pbp", eff2);
 
  TFile* file3 = TFile::Open(inputFile3);
  TEfficiency* eff3;
  file3->TFile::GetObject("eff_comb", eff3);

  TCanvas* canvas = new TCanvas;
  canvas->Divide(3,2);
  canvas->cd(1);
  eff1->Draw();
  canvas->cd(2);
  eff2->Draw();
  canvas->cd(3);
  eff3->Draw();  

  const char* title1 = "pPb_over_comb"; 
  TGraph* graph1;
  graph1 = DivideEff(eff1, eff3, title1, n, xBin);
  canvas->cd(4);
  graph1->Draw("PAL");

  const char* title2 = "Pbp_over_comb";
  TGraph* graph2;
  graph2 = DivideEff(eff2, eff3, title2, n, xBin);
  canvas->cd(5);
  graph2->Draw("PAL"); 

  canvas->SetTitle(title);
  canvas->Print(outputFile);

  return;
}

void DrawEff(){
  const char* inputFile1 = "/home/jg66/DiMuonAnalysis2019/rootFile/PromptNT/JPsi/0_1.4/PromptNTJPsi0_1.4_eff_pPb.root";
  const char* inputFile2 = "/home/jg66/DiMuonAnalysis2019/rootFile/PromptNT/JPsi/0_1.4/PromptNTJPsi0_1.4_eff_Pbp.root";
  const char* inputFile3 = "/home/jg66/DiMuonAnalysis2019/rootFile/PromptNT/JPsi/0_1.4/PromptNTJPsi0_1.4_eff_comb.root";
  const char* title = "PromptNTJPsi0_1.4_canvas";
  const char* outputFile = "./rootFile/PromptNT/JPsi/0_1.4/PromptNTJPsi0_1.4_canvas.png";
  int n = 2;
  double xBin[]  = {6.5, 9.0, 50.0};
  getRatio(inputFile1, inputFile2, inputFile3, title, outputFile, n, xBin);

 
}
