import ROOT
ROOT.gStyle.SetOptStat(0)

cats = [0.5,1.5,2.5]
limObs = [3.85893,9.59719,3.79019]
limObsErr = [0.0171939,0.0331924,0.0051487]
limExp = [1.82086,7.03458,1.6616407]
limExp2Sdown = [1.18124,3.79856,1.0750312]
limExp1Sdown = [1.42502,5.17756,1.2753003]
limExp1Sup = [2.66946,10.9418,2.4447288]
limExp2Sup = [4.04309,10.9418,3.6983808]

#new canvas
c6 = ROOT.TCanvas("c6","c6",800,800)
c6.SetGridx()
c6.SetGridy()
c6.SetRightMargin(0.06)
c6.SetLeftMargin(0.2)

# dummy historgram for axes labels, ranges, etc.
dummy = ROOT.TH1D("","", 3, 0, 3)
dummy.SetBinContent(1,0.0)
dummy.GetXaxis().SetBinLabel(1,'Njets2')
dummy.GetXaxis().SetBinLabel(2,'Njets3')
dummy.GetXaxis().SetBinLabel(3,'Njets2+Njets3')
dummy.GetYaxis().SetTitle('95% C.L. limit on #mu_{qqH}')
dummy.SetLineColor(0)
dummy.SetLineWidth(0)
dummy.SetFillColor(0)
dummy.SetMinimum(0.0)
dummy.SetMaximum(19.0)
dummy.Draw()

gr_exp2 = ROOT.TGraphAsymmErrors()
for i in range(3):
  gr_exp2.SetPoint(i,cats[i],limExp[i])
  gr_exp2.SetPointError(i,0.49,0.49,limExp[i]-limExp2Sdown[i],limExp2Sup[i]-limExp[i])
gr_exp2.SetLineColor(ROOT.kYellow)
gr_exp2.SetFillColor(ROOT.kYellow)
#gr_exp2.SetFillStyle(3001)
gr_exp2.Draw("e2same")

gr_exp1 = ROOT.TGraphAsymmErrors()
for i in range(3):
  gr_exp1.SetPoint(i,cats[i],limExp[i])
  gr_exp1.SetPointError(i,0.49,0.49,limExp[i]-limExp1Sdown[i],limExp1Sup[i]-limExp[i])
gr_exp1.SetLineColor(ROOT.kGreen)
gr_exp1.SetFillColor(ROOT.kGreen)
#gr_exp1.SetFillStyle(3001)
gr_exp1.Draw("e2same")

gr_exp = ROOT.TGraphErrors()
for i in range(3):
  gr_exp.SetPoint(i,cats[i],limExp[i])
  gr_exp.SetPointError(i,0.49,0)
gr_exp.SetLineColor(1)
gr_exp.SetLineWidth(2)
gr_exp.SetLineStyle(2)
gr_exp.Draw("pesame")

gr_obs = ROOT.TGraphErrors()
for i in range(3):
  gr_obs.SetPoint(i,cats[i],limObs[i])
  gr_obs.SetPointError(i,0.49,0)
gr_obs.SetLineColor(1)
gr_obs.SetLineWidth(2)
gr_obs.Draw("pesame")

latex2 = ROOT.TLatex()
latex2.SetNDC()
latex2.SetTextSize(0.5*c6.GetTopMargin())
latex2.SetTextFont(42)
latex2.SetTextAlign(31) # align right                                                                                             
latex2.DrawLatex(0.94, 0.95,"35.9 fb^{-1} (13 TeV)")
latex2.SetTextSize(0.9*c6.GetTopMargin())
latex2.SetTextFont(40)
latex2.SetTextAlign(11) # align right                                                                                             
latex2.DrawLatex(0.20, 0.95, "CMS")
latex2.DrawLatex(0.33, 0.95, "preliminary")
latex2.SetTextSize(0.7*c6.GetTopMargin())
latex2.SetTextFont(52)
latex2.SetTextAlign(11)

legend = ROOT.TLegend(.60,.70,.90,.90)
legend.AddEntry(gr_obs , "Observed 95% CL", "lpe")
legend.AddEntry(gr_exp , "Expected 95% CL", "lpe")
legend.AddEntry(gr_exp1 , "#pm 1#sigma", "f")
legend.AddEntry(gr_exp2 , "#pm 2#sigma", "f")
legend.SetShadowColor(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetLineColor(0)
legend.Draw("same")

line = ROOT.TLine(0,1,3,1)
line.SetLineColor(ROOT.kRed)
line.SetLineWidth(2)
line.Draw()

#c6.Draw() 
ROOT.gPad.RedrawAxis()
dummy.Draw("same")

raw_input()