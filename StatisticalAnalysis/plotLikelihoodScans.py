from array import array
import numpy
import ROOT

ROOT.gStyle.SetOptStat(0)

# create arrays
r_exp1 = array('d',[])
nll_exp1 = array('d',[])
r_obs1 = array('d',[])
nll_obs1 = array('d',[])
zeros1 = array('d',[])
# get expected scan
f_exp1 = ROOT.TFile("k57nj2/LikelihoodScanExp/higgsCombineExp.MultiDimFit.mH120.root","READ")
t_exp1 = f_exp1.Get("limit")
for i in xrange(1,t_exp1.GetEntries()):
    t_exp1.GetEntry(i)
    r_exp1.append(t_exp1.r)
    nll_exp1.append(2.0*t_exp1.deltaNLL)
# get observed scan
f_obs1 = ROOT.TFile("k57nj2/LikelihoodScanObs/higgsCombineObs.MultiDimFit.mH120.root","READ")
t_obs1 = f_obs1.Get("limit")
for i in xrange(1,t_obs1.GetEntries()):
    t_obs1.GetEntry(i)
    r_obs1.append(t_obs1.r)
    nll_obs1.append(2.0*t_obs1.deltaNLL)
    zeros1.append(0.0)
# convert arrays to TVectorD
v_r_exp1 = ROOT.TVectorD(len(r_exp1),r_exp1)
v_r_obs1 = ROOT.TVectorD(len(r_obs1),r_obs1)
v_nll_exp1 = ROOT.TVectorD(len(nll_exp1),nll_exp1)
v_nll_obs1 = ROOT.TVectorD(len(nll_obs1),nll_obs1)
v_zeros1 = ROOT.TVectorD(len(zeros1),zeros1)

# create arrays
r_exp2 = array('d',[])
nll_exp2 = array('d',[])
r_obs2 = array('d',[])
nll_obs2 = array('d',[])
zeros2 = array('d',[])
# get expected scan
f_exp2 = ROOT.TFile("k24nj3/LikelihoodScanExp/higgsCombineExp.MultiDimFit.mH120.root","READ")
t_exp2 = f_exp2.Get("limit")
for i in xrange(1,t_exp2.GetEntries()):
    t_exp2.GetEntry(i)
    r_exp2.append(t_exp2.r)
    nll_exp2.append(2.0*t_exp2.deltaNLL)
# get observed scan
f_obs2 = ROOT.TFile("k24nj3/LikelihoodScanObs/higgsCombineObs.MultiDimFit.mH120.root","READ")
t_obs2 = f_obs2.Get("limit")
for i in xrange(1,t_obs2.GetEntries()):
    t_obs2.GetEntry(i)
    r_obs2.append(t_obs2.r)
    nll_obs2.append(2.0*t_obs2.deltaNLL)
    zeros2.append(0.0)
# convert arrays to TVectorD
v_r_exp2 = ROOT.TVectorD(len(r_exp2),r_exp2)
v_r_obs2 = ROOT.TVectorD(len(r_obs2),r_obs2)
v_nll_exp2 = ROOT.TVectorD(len(nll_exp2),nll_exp2)
v_nll_obs2 = ROOT.TVectorD(len(nll_obs2),nll_obs2)
v_zeros2 = ROOT.TVectorD(len(zeros2),zeros2)

# create arrays
r_exp3 = array('d',[])
nll_exp3 = array('d',[])
r_obs3 = array('d',[])
nll_obs3 = array('d',[])
zeros3 = array('d',[])
# get expected scan
f_exp3 = ROOT.TFile("k57nj2_k24nj3/LikelihoodScanExp/higgsCombineExp.MultiDimFit.mH120.root","READ")
t_exp3 = f_exp3.Get("limit")
for i in xrange(1,t_exp3.GetEntries()):
    t_exp3.GetEntry(i)
    r_exp3.append(t_exp3.r)
    nll_exp3.append(2.0*t_exp3.deltaNLL)
# get observed scan
f_obs3 = ROOT.TFile("k57nj2_k24nj3/LikelihoodScanObs/higgsCombineObs.MultiDimFit.mH120.root","READ")
t_obs3 = f_obs3.Get("limit")
for i in xrange(1,t_obs3.GetEntries()):
    t_obs3.GetEntry(i)
    r_obs3.append(t_obs3.r)
    nll_obs3.append(2.0*t_obs3.deltaNLL)
    zeros3.append(0.0)
# convert arrays to TVectorD
v_r_exp3 = ROOT.TVectorD(len(r_exp3),r_exp3)
v_r_obs3 = ROOT.TVectorD(len(r_obs3),r_obs3)
v_nll_exp3 = ROOT.TVectorD(len(nll_exp3),nll_exp3)
v_nll_obs3 = ROOT.TVectorD(len(nll_obs3),nll_obs3)
v_zeros3 = ROOT.TVectorD(len(zeros3),zeros3)


# new canvas
c9 = ROOT.TCanvas("c9","c9",800, 800)
c9.SetRightMargin(0.06)
c9.SetLeftMargin(0.2)
# dummy for axis labels, ranges, etc.
dummy = ROOT.TH1D("","", 1, -2.0,10.0)
dummy.SetBinContent(1,0.0)
dummy.GetXaxis().SetTitle('#mu_{qqH}')
dummy.GetYaxis().SetTitle('-2 #Delta lnL')
dummy.SetLineColor(0)
dummy.SetLineWidth(0)
dummy.SetFillColor(0)
dummy.SetMinimum(0.0)
dummy.SetMaximum(5.0)
dummy.Draw()
# Draw some lines for 68% CL and 95% CL
latexf = ROOT.TLatex()
latexf.SetTextSize(0.4*c9.GetTopMargin())
latexf.SetTextColor(2)
f1 = ROOT.TF1("f1","1.0",0.0,10.0)
f1.SetLineColor(2)
f1.SetLineWidth(2)
f1.Draw("lsame")
latexf.DrawLatex(2.5, 1.1,"68% CL")
f2 = ROOT.TF1("f1","3.84",0.0,10.0)
f2.SetLineColor(2)
f2.SetLineWidth(2)
f2.Draw("lsame")
latexf.DrawLatex(2.5, 3.94,"95% CL")

# draw expected scan
gr_exp1 = ROOT.TGraphAsymmErrors(v_r_exp1,v_nll_exp1,v_zeros1,v_zeros1,v_zeros1,v_zeros1)
gr_exp1.SetLineColor(4)
gr_exp1.SetLineWidth(2)
gr_exp1.SetLineStyle(2)
gr_exp1.Draw("Lsame")
# draw observed scan
gr_obs1 = ROOT.TGraphAsymmErrors(v_r_obs1,v_nll_obs1,v_zeros1,v_zeros1,v_zeros1,v_zeros1)
gr_obs1.SetLineColor(4)
gr_obs1.SetLineWidth(2)
gr_obs1.Draw("Lsame")

# draw expected scan
gr_exp2 = ROOT.TGraphAsymmErrors(v_r_exp2,v_nll_exp2,v_zeros2,v_zeros2,v_zeros2,v_zeros2)
gr_exp2.SetLineColor(6)
gr_exp2.SetLineWidth(2)
gr_exp2.SetLineStyle(2)
gr_exp2.Draw("Lsame")
# draw observed scan
gr_obs2 = ROOT.TGraphAsymmErrors(v_r_obs2,v_nll_obs2,v_zeros2,v_zeros2,v_zeros2,v_zeros2)
gr_obs2.SetLineColor(6)
gr_obs2.SetLineWidth(2)
gr_obs2.Draw("Lsame")

# draw expected scan
gr_exp3 = ROOT.TGraphAsymmErrors(v_r_exp3,v_nll_exp3,v_zeros3,v_zeros3,v_zeros3,v_zeros3)
gr_exp3.SetLineColor(1)
gr_exp3.SetLineWidth(2)
gr_exp3.SetLineStyle(2)
gr_exp3.Draw("Lsame")
# draw observed scan
gr_obs3 = ROOT.TGraphAsymmErrors(v_r_obs3,v_nll_obs3,v_zeros3,v_zeros3,v_zeros3,v_zeros3)
gr_obs3.SetLineColor(1)
gr_obs3.SetLineWidth(2)
gr_obs3.Draw("Lsame")


latex2 = ROOT.TLatex()
latex2.SetNDC()
latex2.SetTextSize(0.5*c9.GetTopMargin())
latex2.SetTextFont(42)
latex2.SetTextAlign(31) # align right                                                                                                                              
latex2.DrawLatex(0.94, 0.95,"35.9 fb^{-1} (13 TeV)")
latex2.SetTextSize(0.7*c9.GetTopMargin())
latex2.SetTextFont(62)
latex2.SetTextAlign(11) # align right                                                                                                                              
latex2.DrawLatex(0.20, 0.95, "CMS")
latex2.SetTextSize(0.6*c9.GetTopMargin())
latex2.SetTextFont(52)
latex2.SetTextAlign(11)
latex2.DrawLatex(0.3, 0.95, "preliminary")

legend = ROOT.TLegend(.7,.18,.92,.3)
legend.AddEntry(gr_obs1 , "Njets2 Obs.", "l")
legend.AddEntry(gr_exp1 , "Njets2 Exp.", "l")
legend.AddEntry(gr_obs2 , "Njets3 Obs.", "l")
legend.AddEntry(gr_exp2 , "Njets3 Exp.", "l")
legend.AddEntry(gr_obs3 , "Combined Obs.", "l")
legend.AddEntry(gr_exp3 , "Combined Exp.", "l")
legend.SetShadowColor(0)
legend.SetFillColor(0)
legend.SetLineColor(0)
legend.Draw("same")

c9.Draw()
ROOT.gPad.RedrawAxis()
 
raw_input()
