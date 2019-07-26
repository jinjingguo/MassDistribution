#include "TFile.h"
#include "TEfficiency.h"
#include "TCanvas.h"

void DrawA(const char* inputFile1, const char* inputFile2, const char* inputFile3, const char* inputFile4, const char* title, const char* outfile){
  TFile* file1 = TFile::Open(inputFile1);
  TEfficiency* eff1;
  file1->TFile::GetObject("acc", eff1);
  eff1->SetTitle("JPsi_0_1.4");

  TFile* file2 = TFile::Open(inputFile2);
  TEfficiency* eff2;
  file2->TFile::GetObject("acc", eff2);
  eff2->SetTitle("JPsi_1.4_2.4");
  eff2->SetLineColor(kBlue);

  TFile* file3 = TFile::Open(inputFile3);
  TEfficiency* eff3;
  file3->TFile::GetObject("acc", eff3);
  eff3->SetTitle("Psi2S_0_1.4");
  eff3->SetLineColor(kGreen);

  TFile* file4 = TFile::Open(inputFile4);
  TEfficiency* eff4;
  file4->TFile::GetObject("acc", eff4);
  eff4->SetTitle("Psi2S_1.4_2.4");
  eff4->SetLineColor(kRed);

  TCanvas* c = new TCanvas;
  c->SetTitle(title);
  //eff1->Draw();
  eff2->Draw();
  eff1->Draw("same");
  eff3->Draw("same");
  eff4->Draw("same");
  c->BuildLegend(0.7,0.1,0.9,0.3,"y_range");
  eff2->SetTitle("eff;p_{T} [GeV/c];Accpetance"); 
  eff2->SetTitle(title); 
  //eff2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  //eff2->GetYaxis()->SetTitle("Acceptance"); 
  c->Print(outfile);
  return;
  }

void DrawAcc(){
  auto inputFile1 = "./rootFile/NonPromptWT/JPsi/0_1.4/NonPromptWTJPsi0_1.4_acc.root";
  auto inputFile2 = "./rootFile/NonPromptWT/JPsi/1.4_2.4/NonPromptWTJPsi1.4_2.4_acc.root";
  auto inputFile3 = "./rootFile/NonPromptWT/Psi2S/0_1.4/NonPromptWTPsi2S0_1.4_acc.root";
  auto inputFile4 = "./rootFile/NonPromptWT/Psi2S/1.4_2.4/NonPromptWTPsi2S1.4_2.4_acc.root";
  auto title = "NonPromptWT_acc_canvas";
  auto outfile = "NonPromptWT_acc_canvas.png";
  DrawA(inputFile1, inputFile2, inputFile3, inputFile4, title, outfile);
}
