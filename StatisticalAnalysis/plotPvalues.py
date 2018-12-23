import ROOT
ROOT.gStyle.SetOptStat(0)

cats = [0]
pvaluesObs = [0.02874]
pvaluesObsErr = [0.000747182]
pvaluesExpMedian = [0.03966]
pvaluesExpMedianErr = [0.000872778]

for i in range(2,20):  
  obsok = False
  expok = False
  pvobs = 0
  pvobserr = 0
  pvexp = 0
  pvexperr = 0
  
  fobs = ROOT.TFile.Open('k57nj2_k24nj3/SigHybridScan/higgsCombineObs%i.HybridNew.mH120.root' % i)
  if(fobs.Get('limit')):
    obsok = True
    limit = fobs.Get('limit')
    limit.GetEntry(0)
    pvobs = limit.limit
    pvobserr = limit.limitErr
  fobs.Close()
  
  fexp = ROOT.TFile.Open('k57nj2_k24nj3/SigHybridScan/higgsCombineExp%i.HybridNew.mH120.quant0.500.root' % i)
  if(fexp.Get('limit')):
    expok = True
    limit = fexp.Get('limit')
    limit.GetEntry(0)
    pvexp = limit.limit
    pvexperr = limit.limitErr
  fexp.Close()

  if(obsok and expok):
    print 'Making for %.3f, Obs = %.5f, Exp = %.5f' % ((i/24.),pvobs,pvexp)
    cats.append(i/24.)
    pvaluesObs.append(pvobs)
    pvaluesObsErr.append(pvobserr)
    pvaluesExpMedian.append(pvexp)
    pvaluesExpMedianErr.append(pvexperr)

# new canvas
c7 = ROOT.TCanvas("c7","c7",800, 800)
c7.SetLogy()
c7.SetRightMargin(0.06)
c7.SetLeftMargin(0.2)
# dummy histogram, for axis labels, ranges, etc.
dummy = ROOT.TH1D("","", len(cats), cats[0],cats[len(cats)-1])
dummy.SetBinContent(1,0.0)
#dummy.GetXaxis().SetBinLabel(1,'Njets2')
#dummy.GetXaxis().SetBinLabel(2,'Njets3')
#dummy.GetXaxis().SetBinLabel(3,'Njets2+Njets3')
dummy.GetYaxis().SetTitle('Local p-value')
dummy.SetLineColor(0)
dummy.SetLineWidth(0)
dummy.SetFillColor(0)
dummy.SetMinimum(0.0001)
dummy.SetMaximum(1.0)
dummy.GetXaxis().SetTitle("Discriminant Cut")
dummy.Draw()

# Draw some lines corresponding to 1,2,3 sigma 
latexf = ROOT.TLatex()
latexf.SetTextSize(0.4*c7.GetTopMargin())
latexf.SetTextColor(2)
f1 = ROOT.TF1("f1","0.15866",0,cats[len(cats)-1])
f1.SetLineColor(2)
f1.SetLineWidth(2)
f1.Draw("lsame")
latexf.DrawLatex(0.75, 0.15866*1.1,"1#sigma")
f2 = ROOT.TF1("f1","0.02275",0,cats[len(cats)-1])
f2.SetLineColor(2)
f2.SetLineWidth(2)
f2.Draw("lsame")
latexf.DrawLatex(0.75, 0.02275*1.1,"2#sigma")
f3 = ROOT.TF1("f1","0.0013499",0,cats[len(cats)-1])
f3.SetLineColor(2)
f3.SetLineWidth(2)
f3.Draw("lsame")
latexf.DrawLatex(0.75, 0.0013499*1.1,"3#sigma")

# Draw the expected p-value graph
gr_exp = ROOT.TGraphErrors()
for i in range(len(cats)):
  gr_exp.SetPoint(i,cats[i],pvaluesExpMedian[i])
  #gr_exp.SetPointError(i,0.5,pvaluesExpMedianErr[i])

gr_exp.SetLineColor(4)
gr_exp.SetLineWidth(1)
gr_exp.SetFillStyle(3002)
gr_exp.SetMarkerColor(4)
gr_exp.SetMarkerSize(0.6)
gr_exp.Draw("lpsame")

# Draw the observed p-value graph
gr_obs = ROOT.TGraphErrors()
for i in range(len(cats)):
  gr_obs.SetPoint(i,cats[i],pvaluesObs[i])
  #gr_obs.SetPointError(i,0.5,pvaluesObsErr[i])

gr_obs.SetLineColor(1)
gr_obs.SetLineWidth(1)
gr_obs.SetFillStyle(3002)
gr_obs.SetMarkerSize(0.6)
gr_obs.Draw("lpsame")

latex2 = ROOT.TLatex()
latex2.SetNDC()
latex2.SetTextSize(0.5*c7.GetTopMargin())
latex2.SetTextFont(42)
latex2.SetTextAlign(31) # align right                                                                                                                              
latex2.DrawLatex(0.94, 0.95,"35.9 fb^{-1} (13 TeV)")
latex2.SetTextSize(0.7*c7.GetTopMargin())
latex2.SetTextFont(62)
latex2.SetTextAlign(11) # align right                                                                                                                              
latex2.DrawLatex(0.20, 0.95, "CMS")
latex2.DrawLatex(0.3, 0.95, "preliminary")
latex2.SetTextSize(0.6*c7.GetTopMargin())
latex2.SetTextFont(52)
latex2.SetTextAlign(11)

legend = ROOT.TLegend(.7,.18,.92,.3)
legend.AddEntry(gr_exp, "Expected", "lp")
legend.AddEntry(gr_obs, "Observed", "lp")
legend.SetShadowColor(0)
legend.SetFillColor(0)
legend.SetLineColor(0)
legend.Draw("same")


c7.Draw()

raw_input()
