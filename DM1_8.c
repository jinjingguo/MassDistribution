void DM1_8(){
    TFile *file1 = new TFile("/storage1/users/wl33/DiMuTrees/pPb2016/Tree/VertexCompositeTree_PADoubleMuon_PARun2016C_DiMuMassMin2.root");
    TTree *tree1 = (TTree*)file1 -> Get("dimucontana/VertexCompositeNtuple");
    TCanvas *canvas = new TCanvas("canvas");
    canvas -> Divide(3,3);

    canvas -> cd(1);
    tree1 -> Draw("mass>>hmass1","mass>2.2&&mass<4.2&&pT>0.2&&pT<1.8&&y>-1.4&&y<1.4&&Ntrkoffline>250");
    TH1F *hmass1 = (TH1F*)gDirectory -> Get("hmass1");
    hmass1 -> GetXaxis() -> SetRangeUser(2.2, 4.2);

    canvas -> cd(2);
    tree1 -> Draw("mass>>hmass2","mass>2.2&&mass<4.2&&pT>1.8&&pT<3&&y>-1.4&&y<1.4&&Ntrkoffline>250");
    TH1F *hmass2 = (TH1F*)gDirectory -> Get("hmass2");
    hmass2 -> GetXaxis() -> SetRangeUser(2.2, 4.2);

    canvas -> cd(3);
    tree1 -> Draw("mass>>hmass3","mass>2.2&&mass<4.2&&pT>3&&pT<4.5&&y>-1.4&&y<1.4&&Ntrkoffline>250");
    TH1F *hmass3 = (TH1F*)gDirectory -> Get("hmass3");
    hmass3 -> GetXaxis() -> SetRangeUser(2.2, 4.2);

    canvas -> cd(4);
    tree1 -> Draw("mass>>hmass4","mass>2.2&&mass<4.2&&pT>4.5&&pT<6&&y>-1.4&&y<1.4&&Ntrkoffline>250");
    TH1F *hmass4 = (TH1F*)gDirectory -> Get("hmass4");
    hmass4 -> GetXaxis() -> SetRangeUser(2.2, 4.2);

    canvas -> cd(5);
    tree1 -> Draw("mass>>hmass5","mass>2.2&&mass<4.2&&pT>6&&pT<8&&y>-1.4&&y<1.4&&Ntrkoffline>250");
    TH1F *hmass5 = (TH1F*)gDirectory -> Get("hmass5");
    hmass5 -> GetXaxis() -> SetRangeUser(2.2, 4.2);

    canvas -> cd(6);
    tree1 -> Draw("mass>>hmass6","mass>2.2&&mass<4.2&&pT>8&&pT<10&&y>-1.4&&y<1.4&&Ntrkoffline>250");
    TH1F *hmass6 = (TH1F*)gDirectory -> Get("hmass6");
    hmass6 -> GetXaxis() -> SetRangeUser(2.2, 4.2);

    canvas -> cd(7);
    tree1 -> Draw("mass>>hmass7","mass>2.2&&mass<4.2&&pT>10&&pT<20&&y>-1.4&&y<1.4&&Ntrkoffline>250");
    TH1F *hmass7 = (TH1F*)gDirectory -> Get("hmass7");
    hmass7 -> GetXaxis() -> SetRangeUser(2.2, 4.2);

    canvas -> Print("250_.pdf");
}
