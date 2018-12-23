//How to use it?
//root -l -q plotPrePostFit2.cc\(\"fitDiagnostics.root\",\"shapes_prefit\",false,false,true\)
//The 2 root files are the output from Higgs Combine Tool after making fitDiagnostics
//The "shapes_prefit" (other options "shapes_fit_b" = fitting only bkg or "shapes_fit_s" = fitting sig+bkg) specifies the folder within Combine files
//The flags determine when using each channel (4mu, 4e, 2e2mu). Note that the order changes accordingly to what was given in the datacards for Combine



void plotCombinedPrePostFit(void){
  /*
  TString orig_file_ch1_cat1 = "k57nj2/BestFit/datacard_Higgs125GeV_VBFStudy_k57nj2_ch4mu.root";
  TString orig_file_ch2_cat1 = "k57nj2/BestFit/datacard_Higgs125GeV_VBFStudy_k57nj2_ch4e.root";
  TString orig_file_ch3_cat1 = "k57nj2/BestFit/datacard_Higgs125GeV_VBFStudy_k57nj2_ch2e2mu.root";

  TString orig_file_ch1_cat2 = "k24nj3/BestFit/datacard_Higgs125GeV_VBFStudy_k24nj3_ch4mu.root";
  TString orig_file_ch2_cat2 = "k24nj3/BestFit/datacard_Higgs125GeV_VBFStudy_k24nj3_ch4e.root";
  TString orig_file_ch3_cat2 = "k24nj3/BestFit/datacard_Higgs125GeV_VBFStudy_k24nj3_ch2e2mu.root";
  */  
  
  TString orig_file_ch1_cat1 = "Cutting/datacard_Higgs125GeV_VBFStudy_k57nj2_ch4mu_cutbin2.root";
  TString orig_file_ch2_cat1 = "Cutting/datacard_Higgs125GeV_VBFStudy_k57nj2_ch4e_cutbin2.root";
  TString orig_file_ch3_cat1 = "Cutting/datacard_Higgs125GeV_VBFStudy_k57nj2_ch2e2mu_cutbin2.root";

  TString orig_file_ch1_cat2 = "Cutting/datacard_Higgs125GeV_VBFStudy_k24nj3_ch4mu_cutbin2.root";
  TString orig_file_ch2_cat2 = "Cutting/datacard_Higgs125GeV_VBFStudy_k24nj3_ch4e_cutbin2.root";
  TString orig_file_ch3_cat2 = "Cutting/datacard_Higgs125GeV_VBFStudy_k24nj3_ch2e2mu_cutbin2.root";
  
  TString combineFileCat1="k57nj2_k24nj3/BestFit/fitDiagnostics.root";
  
  bool ch1=true;
  bool ch2=true; 
  bool ch3=true; 

  TString fitfolder="shapes_prefit";
  TString out_name="fitDiagnostics";
  int rebin=1;
  
  std::cout<<"Plotting from "<<fitfolder<<std::endl;
  TString channels;
  TString sch1 = "4mu";
  TString sch2 = "4e";
  TString sch3 = "2e2mu";
  if(ch1 && !ch2 && !ch3) channels = sch1;
  else if(!ch1 && ch2 && !ch3) channels = sch2;
  else if(!ch1 && !ch2 && ch3) channels = sch3;
  else channels = "4l";  

  TH1::SetDefaultSumw2();
  

  //_______________________________________________________________________
  // Handles MC first 
  //_______________________________________________________________________
  std::vector<TString> processes = {
    "zjets",
    "qqZZ",
    "ggZZ",
    "ttH",
    "WH",
    "ZH",
    "ggH",
    "qqH"
  };
  
  unsigned int Nprocesses = processes.size();
  TH1F *histosch1[Nprocesses];
  TH1F *histosch2[Nprocesses];
  TH1F *histosch3[Nprocesses];
  TH1F *histos4l[Nprocesses];

  TFile *fccat1 = TFile::Open(combineFileCat1);
  TFile *fc1cat1 = TFile::Open(orig_file_ch1_cat1);
  TFile *fc1cat2 = TFile::Open(orig_file_ch1_cat2);
  TFile *fc2cat1 = TFile::Open(orig_file_ch2_cat1);
  TFile *fc2cat2 = TFile::Open(orig_file_ch2_cat2);
  TFile *fc3cat1 = TFile::Open(orig_file_ch3_cat1);
  TFile *fc3cat2 = TFile::Open(orig_file_ch3_cat2);

  TH1F *htmp = (TH1F*)fc1cat1->Get("ch4mu/qqH");
  int NbBins = htmp->GetNbinsX();
  double xlower = htmp->GetBinLowEdge(1);
  double xupper = htmp->GetBinLowEdge(NbBins)+htmp->GetBinWidth(1);
  std::cout<<Form("Bins = %i, xlower = %.3f, xupper = %.3f",NbBins,xlower,xupper)<<std::endl;
  TH1F *htemplate = new TH1F("","",NbBins,xlower,xupper);
  
  for(unsigned int iproc=0; iproc<Nprocesses; ++iproc){
    if(ch1){
      histosch1[iproc] = (TH1F*)htemplate->Clone();
      histosch1[iproc]->SetName(processes[iproc]+"_"+sch1);
    }
    if(ch2){
      histosch2[iproc] = (TH1F*)htemplate->Clone();
      histosch2[iproc]->SetName(processes[iproc]+"_"+sch2);
    }
    if(ch3){
      histosch3[iproc] = (TH1F*)htemplate->Clone();
      histosch3[iproc]->SetName(processes[iproc]+"_"+sch3);
    }
    
    TH1F *hch1tmp = (TH1F*)fc1cat1->Get("ch"+sch1+"/"+processes[iproc]);
    hch1tmp->Add((TH1F*)fc1cat2->Get("ch"+sch1+"/"+processes[iproc]));

    TH1F *hch2tmp = (TH1F*)fc2cat1->Get("ch"+sch2+"/"+processes[iproc]);
    hch2tmp->Add((TH1F*)fc2cat2->Get("ch"+sch2+"/"+processes[iproc]));

    TH1F *hch3tmp = (TH1F*)fc3cat1->Get("ch"+sch3+"/"+processes[iproc]);
    hch3tmp->Add((TH1F*)fc3cat2->Get("ch"+sch3+"/"+processes[iproc]));
    
    std::cout<<"Loaded datacards..."<<std::endl;
	
    for(unsigned int ibin=1; ibin<=NbBins; ++ibin){
      if(ch1){
	histosch1[iproc]->SetBinContent( ibin, hch1tmp->GetBinContent(ibin) );
	histosch1[iproc]->SetBinError( ibin, hch1tmp->GetBinError(ibin) );
      }
      if(ch2){
	histosch2[iproc]->SetBinContent( ibin, hch2tmp->GetBinContent(ibin) );
	histosch2[iproc]->SetBinError( ibin, hch2tmp->GetBinError(ibin) );
      }
      if(ch3){
	histosch3[iproc]->SetBinContent( ibin, hch3tmp->GetBinContent(ibin) );
	histosch3[iproc]->SetBinError( ibin, hch3tmp->GetBinError(ibin) );
      }
    }
    
    histos4l[iproc] = (TH1F*)htemplate->Clone();
    histos4l[iproc]->SetName(processes[iproc]+"_sum");
    if(ch1) histos4l[iproc]->Add(histosch1[iproc]);
    if(ch2) histos4l[iproc]->Add(histosch2[iproc]);
    if(ch3) histos4l[iproc]->Add(histosch3[iproc]);
    
    if(processes[iproc] == "qqH"){
      histos4l[iproc]->SetLineColor(kTeal+2);
      histos4l[iproc]->SetFillColor(kTeal+1);
    }
    if(processes[iproc] == "ggH"){
      histos4l[iproc]->SetLineColor(kRed-3);
      histos4l[iproc]->SetFillColor(kRed-9);
    }
    if(processes[iproc] == "ZH"){
      histos4l[iproc]->SetLineColor(kOrange-3);
      histos4l[iproc]->SetFillColor(kOrange-2);
    }
    if(processes[iproc] == "WH"){
      histos4l[iproc]->SetLineColor(kOrange-2);
      histos4l[iproc]->SetFillColor(kOrange-2);
    }
    if(processes[iproc] == "ttH"){
      histos4l[iproc]->SetLineColor(kViolet-6);
      histos4l[iproc]->SetFillColor(kViolet-4);
    }
    if(processes[iproc] == "qqZZ"){
      histos4l[iproc]->SetLineColor(kAzure+2);
      histos4l[iproc]->SetFillColor(kAzure-9);
    }
    if(processes[iproc] == "ggZZ"){
      histos4l[iproc]->SetLineColor(kAzure-6);
      histos4l[iproc]->SetFillColor(kAzure-3);
    }
    if(processes[iproc] == "zjets"){
      histos4l[iproc]->SetLineColor(kSpring-7);
      histos4l[iproc]->SetFillColor(kSpring-6);      
    }
    
    histos4l[iproc]->SetLineWidth(0);
    std::cout<<"Channels added for "<<processes[iproc]<<", Yields = "<<histos4l[iproc]->Integral()<<std::endl;
  }

  //To get Integrals
  TH1F *hvbf = (TH1F*)htemplate->Clone();
  hvbf->SetName("vbf");
  TH1F *hbkg = (TH1F*)htemplate->Clone();
  hbkg->SetName("hbkg");
  
  THStack *mcStack = new THStack();
  for(unsigned int iproc=0; iproc<Nprocesses; ++iproc){
    //Gets the yields from each process
    double yield = 0, yield_err = 0;
    yield = histos4l[iproc]->IntegralAndError(-1,-1,yield_err);
    std::cout<< processes[iproc] << Form(" = %.2f+/-%.2f",yield,yield_err)<<std::endl;
    
    if(processes[iproc] == "qqH") hvbf->Add(histos4l[iproc]);
    else hbkg->Add(histos4l[iproc]);
    
    //Rebining for plot
    if(rebin != 1) histos4l[iproc]->Rebin(rebin);
    mcStack->Add(histos4l[iproc]);
  }
  
  //_____________________________________________
  //Data from Combine tool is a TGraphAsymmErrors
  //_____________________________________________
  TH1F *hdata = (TH1F*)htemplate->Clone();
  hdata->SetName("orig_data");
  if(ch1){
    TFile *ofilech1cat1 = TFile::Open(orig_file_ch1_cat1);
    hdata->Add( (TH1F*)ofilech1cat1->Get("ch"+sch1+"/data_obs") );
    TFile *ofilech1cat2 = TFile::Open(orig_file_ch1_cat2);
    hdata->Add( (TH1F*)ofilech1cat2->Get("ch"+sch1+"/data_obs") );
  }
  if(ch2){
    TFile *ofilech2cat1 = TFile::Open(orig_file_ch2_cat1);
    hdata->Add( (TH1F*)ofilech2cat1->Get("ch"+sch2+"/data_obs") );
    TFile *ofilech2cat2 = TFile::Open(orig_file_ch2_cat2);
    hdata->Add( (TH1F*)ofilech2cat2->Get("ch"+sch2+"/data_obs") );
  }
  if(ch3){
    TFile *ofilech3cat1 = TFile::Open(orig_file_ch3_cat1);
    hdata->Add( (TH1F*)ofilech3cat1->Get("ch"+sch3+"/data_obs") );
    TFile *ofilech3cat2 = TFile::Open(orig_file_ch3_cat2);
    hdata->Add( (TH1F*)ofilech3cat2->Get("ch"+sch3+"/data_obs") );
  }
  double yield = 0, yield_err = 0;
  yield = hdata->IntegralAndError(-1,-1,yield_err);
  std::cout<< Form("Data = %.2f+/-%.2f",yield,yield_err)<<std::endl;

  //Rebining for plot
  if(rebin != 1) hdata->Rebin(rebin);

  //Builds ratio plot
  TH1F *hSumMC = (TH1F*)((TH1F*)mcStack->GetStack()->Last())->Clone();
  hSumMC->SetName("hSumMC");
  TH1F *hratio = (TH1F*)hdata->Clone();
  hratio->SetName("hratio");
  //hratio->Add(hSumMC,-1);
  hratio->Divide(hSumMC);
  
  
  //Builds data asym error graph plot
  TGraphAsymmErrors *dataAsymGraphErr = new TGraphAsymmErrors();
  unsigned int ipoint = 0;
  double maxwindowy = 0;
  for(unsigned int ibin=1; ibin<=NbBins; ++ibin){
    double dx = hdata->GetBinCenter(ibin);
    double dy = hdata->GetBinContent(ibin);
    if(dy != 0){
      double ely = -0.5 + sqrt(dy+0.25);
      double ehy = +0.5 + sqrt(dy+0.25);
      dataAsymGraphErr->SetPoint(ipoint,dx,dy);
      dataAsymGraphErr->SetPointError(ipoint,0,0,ely,ehy);
      ++ipoint;
      if(dy+ehy > maxwindowy) maxwindowy = dy+ehy;
    }
  }
  std::cout<<"Window y = "<<maxwindowy<<std::endl;
  
  
  //Build the plot
  TLegend *leg = new TLegend(0.4,0.65,0.90,0.90);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetEntrySeparation(0.1);
  leg->SetNColumns(3);
  leg->AddEntry(dataAsymGraphErr,"Data","pe");
  for(unsigned int iproc=1; iproc<=Nprocesses; ++iproc)
    if(processes[Nprocesses-iproc] == "WH") continue;
    else if(processes[Nprocesses-iproc] == "ZH")
      leg->AddEntry(histos4l[Nprocesses-iproc],"VH","f");
    else
      leg->AddEntry(histos4l[Nprocesses-iproc],processes[Nprocesses-iproc],"f");
  leg->AddEntry(hSumMC,"Stat.+Syst. Uncertainty","f");

  TPaveText *cms_tag = new TPaveText(.3,.95,.99,.96,"NDC");
  cms_tag->AddText("CMS #bf{Preliminary  #sqrt{s} = 13 TeV, L = 35.9fb^{-1}}");
  cms_tag->SetFillStyle(0);
  cms_tag->SetBorderSize(0);
  cms_tag->SetTextSize(0.04);
  
  TF1 *ratio_fit = new TF1("ratio_fit","[0]",0.0,1.0);
  hratio->Fit(ratio_fit,"Nel","",0.0,1.0);
  double fit_value = ratio_fit->GetParameter(0);
  double fit_error = ratio_fit->GetParError(0);
  TGraphErrors *gfit = new TGraphErrors();
  gfit->SetPoint(0,0.0,fit_value);
  gfit->SetPointError(0,0.0,fit_error);
  gfit->SetPoint(1,1.0,fit_value);
  gfit->SetPointError(1,0.0,fit_error);
  gfit->SetMarkerSize(0);
  gfit->SetLineStyle(7);
  gfit->SetFillColor(kGray+2);
  gfit->SetFillStyle(3001);
  
  TLegend *rleg = new TLegend(0.8,0.75,0.95,0.95);
  rleg->SetFillColor(0);
  rleg->SetFillStyle(0);
  rleg->AddEntry(gfit,"Fit","fl");

  
  TLine *lfit = new TLine(xlower,1.0,xupper,1.0);
  lfit->SetLineColor(kRed);

  //std::cout<<"Creating TCanvas...."<<std::endl;
  //gROOT->SetBatch();
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cv","",10,10,700,700);
  TPad *pad1 = new TPad("","",0.05,0.05,0.95,0.30);
  TPad *pad2 = new TPad("","",0.05,0.30,0.95,0.97);
  pad1->SetTopMargin(0.02);
  pad1->SetBottomMargin(0.3);
  pad2->SetBottomMargin(0.02);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  hratio->Draw("pe");
  gfit->Draw("le3,same");
  lfit->Draw();
  hratio->Draw("pe,same");
  rleg->Draw();
  hratio->SetMarkerStyle(20);
  hratio->SetMarkerColor(kBlack);
  hratio->SetLineColor(kBlack);
  hratio->GetYaxis()->SetTitle("Data/MC");
  hratio->GetXaxis()->SetTitle("Discriminant");
  hratio->GetYaxis()->SetLabelSize(0.1);
  hratio->GetYaxis()->SetTitleSize(0.14);
  hratio->GetYaxis()->SetTitleOffset(0.35);
  hratio->GetXaxis()->SetLabelSize(0.14);
  hratio->GetXaxis()->SetTitleOffset(1.0);
  hratio->GetXaxis()->SetTitleSize(0.14);
  gPad->RedrawAxis();

  pad2->cd();
  mcStack->Draw("hist");
  hSumMC->Draw("pe2,same");
  dataAsymGraphErr->Draw("pe,same");
  cms_tag->Draw();
  double binwidth = mcStack->GetXaxis()->GetBinWidth(1);
  mcStack->GetYaxis()->SetTitle(Form("Events/%.2f",binwidth));
  mcStack->GetXaxis()->SetTitleSize(0.0);
  mcStack->GetXaxis()->SetLabelSize(0.0);
  mcStack->GetYaxis()->SetTitleSize(0.055);
  mcStack->GetYaxis()->SetTitleOffset(0.85);
  mcStack->SetMaximum(1.5*maxwindowy);
  hSumMC->SetFillColor(kBlack);
  hSumMC->SetFillStyle(3001);
  hSumMC->SetMarkerSize(0);
  hSumMC->SetLineColor(kBlack);
  dataAsymGraphErr->SetMarkerColor(kBlack);
  leg->Draw();
  //gPad->SetLogy();

  cv->Update();  
  //cv->Print(out_name+"_"+fitfolder+".png");
  
  return;
}
