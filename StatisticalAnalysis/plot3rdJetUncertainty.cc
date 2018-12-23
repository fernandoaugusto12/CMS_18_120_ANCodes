void plot3rdJetUncertainty(void){
  TH1::SetDefaultSumw2();  
  
  TString sr = "vbf";
  TString channel = "2e2mu";
  TString net = "k24nj3";
  TString path = sr+"/"+channel+"/"+net+"_";
  int rebin = 1;
    
  TFile *fhistos = TFile::Open("../SamplesOverview/September20_2018/Histograms.root");
  TH1D *qqH = (TH1D*)fhistos->Get(path+"qqH");
  TH1D *qqH3J = (TH1D*)fhistos->Get(path+"qqH3J");  
  TH1D *OqqH = (TH1D*)qqH->Clone("OqqH");
  
  //Report yields
  double nqqHe = 0;
  float nqqH = qqH->IntegralAndError(-1,-1,nqqHe);
  double nqqH3Je = 0;
  float nqqH3J = qqH3J->IntegralAndError(-1,-1,nqqH3Je);
  std::cout<<Form("qqH   = %.2f+/-%.2f",nqqH,nqqHe)<<std::endl;
  std::cout<<Form("qqH3J = %.2f+/-%.2f",nqqH3J,nqqH3Je)<<std::endl;
  
  qqH->Rebin(rebin);
  qqH3J->Rebin(rebin);
      
  qqH->SetLineColor(kTeal+2);
  qqH->SetFillColor(kTeal+1);
  qqH->SetMarkerColor(kTeal+2);

  qqH3J->SetLineColor(kOrange-3);
  qqH3J->SetFillColor(kOrange-9);
  qqH3J->SetMarkerColor(kOrange-3);

  gStyle->SetOptStat(0);

  TPaveText *cms_tag = new TPaveText(.34,.95,.94,.96,"NDC");
  cms_tag->AddText("CMS #bf{Preliminary  #sqrt{s} = 13 TeV, L = 35.9fb^{-1}}");
  cms_tag->SetFillStyle(0);
  cms_tag->SetBorderSize(0);
  cms_tag->SetTextSize(0.04);
  
  //std::cout<<"Creating TLegend...."<<std::endl;
  TLegend *leg = new TLegend(0.7,0.75,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  
  leg->AddEntry(qqH,"VBF-H2J","f");
  leg->AddEntry(qqH3J,"VBF-H3J","f");    

  TH1D *hratio = (TH1D*)qqH->Clone("hratio");
  hratio->Add(qqH3J, -1);
  hratio->Divide(qqH);

  TF1 *fit = new TF1("fit","[0]*x+[1]",0.00,1.00);
  hratio->Fit(fit,"ne","",0.00,1.00);

  TGraphErrors *gfit = new TGraphErrors();
  gfit->SetPoint(0,0.00,fit->Eval(0.00));
  gfit->SetPoint(1,1.00,fit->Eval(1.00));
  gfit->SetMarkerSize(0);
  
  TLegend *rleg = new TLegend(0.5,0.8,0.9,0.98);
  rleg->SetFillColor(0);
  rleg->SetFillStyle(0);
  rleg->AddEntry(gfit,Form("Fit (%.3fx #pm %.3f)",fit->GetParameter(0),fit->GetParameter(1)),"l");
  
  TCanvas *cv = new TCanvas("cv","",10,10,700,700);
  TPad *pad1 = new TPad("","",0.05,0.05,0.95,0.30);
  TPad *pad2 = new TPad("","",0.05,0.30,0.95,0.97);
  pad1->SetTopMargin(0.02);
  pad1->SetBottomMargin(0.3);
  pad2->SetBottomMargin(0.02);
  pad1->Draw();
  pad2->Draw();
  
  pad1->cd();
  hratio->SetMarkerColor(kBlack);
  hratio->SetLineColor(kBlack);
  hratio->SetMarkerSize(0.9);
  hratio->GetXaxis()->SetLabelSize(0.14);
  hratio->GetXaxis()->SetTitle("Discriminant");
  hratio->GetXaxis()->SetTitleSize(0.13);
  hratio->GetXaxis()->SetTitleOffset(1.0);
  hratio->GetYaxis()->SetLabelSize(0.13);
  hratio->GetYaxis()->SetTitle("#frac{qqH-qqH_{3J}}{qqH}");
  hratio->GetYaxis()->SetTitleSize(0.13);
  hratio->GetYaxis()->SetTitleOffset(0.5);
  hratio->GetXaxis()->SetTickLength(0.06);
  hratio->GetXaxis()->SetTickLength(0.06);
  hratio->GetXaxis()->SetLimits(0.00,1.00);
  hratio->GetYaxis()->SetLimits(-2.00,2.00);
  hratio->Draw("pe");
  gfit->Draw("l,same");
  gfit->SetLineColor(kRed);
  rleg->Draw();

  /*Create a histogram to hold the confidence intervals*/
  //TH1D *hint = new TH1D("hint","Fitted gaussian with .95 conf.band", qqH->GetNbinsX(), 0, 1);
  //(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
  //Now the "hint" histogram has the fitted function values as the
  //bin contents and the confidence intervals as bin errors
  //hint->Draw("e3 same");
  //hratio->Draw("pe,same");
  

  
  pad2->cd();
  float ymin = 0;
  float ymax = (qqH->GetMaximum() > qqH3J->GetMaximum())? qqH->GetMaximum() : qqH3J->GetMaximum();
              
  qqH->SetMaximum(1.4*ymax);
  qqH->SetMinimum(0.);
  qqH->GetXaxis()->SetLabelSize(0);
  qqH->GetYaxis()->SetTitle(Form("Events/%.3f",qqH->GetBinWidth(1)));
  qqH->GetYaxis()->SetTitleOffset(1.4);
  qqH->SetMarkerSize(0.9);
  qqH3J->SetMarkerSize(0.9);
  qqH->Draw();
  qqH3J->Draw("same");
  leg->Draw();
  cms_tag->Draw();
  pad2->RedrawAxis();
  cv->Update();
  
  TH1D *qqH_Up = (TH1D*)OqqH->Clone("qqH_j3kShapeUp");
  TH1D *qqH_Down = (TH1D*)OqqH->Clone("qqH_j3kShapeDown");
  for(int ibin=1; ibin<=OqqH->GetNbinsX(); ++ibin){
    double x = OqqH->GetBinCenter(ibin);
    double y = OqqH->GetBinContent(ibin);
    double sigma = fit->Eval(x);
    //std::cout<<Form("x = %.3f, y = %.3f+/-%.3f, sigma = %.3f",x,y,y*fabs(sigma),sigma)<<std::endl;
    qqH_Up->SetBinContent( ibin, y*(1+fabs(sigma)) );
    qqH_Down->SetBinContent( ibin, y*(1-fabs(sigma)) );
  }
  
  TFile *m4lDatacard = new TFile("../SamplesOverview/September20_2018/m4lDatacard_"+sr+"_ch"+channel+"3rdJetUncertainty_"+net+".root","recreate");
  m4lDatacard->mkdir("ch"+channel);
  m4lDatacard->cd("ch"+channel);
  qqH_Up->Write();
  qqH_Down->Write();
  m4lDatacard->Close();
  
    
  return;
}
