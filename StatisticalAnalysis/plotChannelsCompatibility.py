import ROOT

poi = "r"
infile='k57nj2_k24nj3/ChannelCompatibility/higgsCombineObs.ChannelCompatibilityCheck.mH120.root'
infile150='k57nj2_k24nj3/ChannelCompatibility/higgsCombinerMinMinus2.ChannelCompatibilityCheck.mH120.root'

filein = ROOT.TFile(infile,"READ")
fit_nominal = filein.Get("fit_nominal")
fit_alternate = filein.Get("fit_alternate")
rFit = fit_nominal.floatParsFinal().find(poi)

filein150 = ROOT.TFile(infile150,"READ")
fit_nominal150 = filein150.Get("fit_nominal")
fit_alternate150 = filein150.Get("fit_alternate")
rFit150 = fit_nominal150.floatParsFinal().find(poi)

prefix = "_ChannelCompatibilityCheck_"+poi+"_"

nChann = 0
iter = fit_alternate.floatParsFinal().createIterator()

while True:
  a = iter.Next()
  if a==None: break
  if (prefix in a.GetName()): nChann+=1

xmin = min(rFit.getMin(),rFit150.getMin())  
xmax = min(rFit.getMax(),rFit150.getMax())
frame = ROOT.TH2F("frame",";best fit #sigma/#sigma(SM);",1,xmin,xmax,nChann,0,nChann)

iter.Reset()
iChann = 0
points = ROOT.TGraphAsymmErrors(nChann)
points150 = ROOT.TGraphAsymmErrors(nChann)

while True:
  a = iter.Next()
  if a==None: break
  if (prefix in a.GetName()):
    channel = a.GetName()
    channel.replace(prefix,"")
    points.SetPoint(iChann, a.getVal(), iChann+0.6)
    points.SetPointError(iChann, -a.getAsymErrorLo(), a.getAsymErrorHi(), 0, 0)
    print "channel = %s, mu = %.2f, ErrLow = %.2f, ErrHigh = %.2f" % (channel,a.getVal(),a.getAsymErrorLo(),a.getAsymErrorHi())
    iChann+=1
    frame.GetYaxis().SetBinLabel(iChann, channel.replace(prefix,""))

iter = fit_alternate150.floatParsFinal().createIterator()
iter.Reset()
iChann = 0
while True:
  a = iter.Next()
  if a==None: break
  if (prefix in a.GetName()):
    channel = a.GetName()
    channel.replace(prefix,"")
    points150.SetPoint(iChann, a.getVal(), iChann+0.4)
    points150.SetPointError(iChann, -a.getAsymErrorLo(), a.getAsymErrorHi(), 0, 0)
    print "channel = %s, mu = %.2f, ErrLow = %.2f, ErrHigh = %.2f" % (channel,a.getVal(),a.getAsymErrorLo(),a.getAsymErrorHi())
    iChann+=1
    
c8 = ROOT.TCanvas("c8","c8",800,800)
c8.cd()
c8.SetTopMargin(0.07)
c8.SetBottomMargin(0.12)
c8.SetLeftMargin(0.2)
c8.SetRightMargin(0.05)
points.SetLineColor(ROOT.kRed)
points.SetMarkerColor(ROOT.kRed)
points.SetLineWidth(3)
points.SetMarkerStyle(20)
points.SetMarkerSize(1.2)
points150.SetLineColor(ROOT.kViolet)
points150.SetMarkerColor(ROOT.kViolet)
points150.SetLineWidth(3)
points150.SetMarkerStyle(20)
points150.SetMarkerSize(1.2)
frame.GetXaxis().SetTitleSize(0.05)
frame.GetXaxis().SetLabelSize(0.04)                                                                                                  
frame.GetYaxis().SetLabelSize(0.06)                                                                                                  
frame.Draw()  

ROOT.gStyle.SetOptStat(0)                                                                                                            
globalFitBand = ROOT.TBox(rFit.getVal()+rFit.getAsymErrorLo(), 0, rFit.getVal()+rFit.getAsymErrorHi(), nChann)                       
globalFitBand.SetFillColor(65)                                                                                                       
globalFitBand.SetLineStyle(0)                                                                                                        
globalFitBand.Draw("SAME")                                                                                                           
globalFitLine = ROOT.TLine(rFit.getVal(), 0, rFit.getVal(), nChann)                                                                  
globalFitLine.SetLineWidth(4)                                                                                                       
globalFitLine.SetLineColor(214)                                                                                                    
globalFitLine.Draw("SAME")                                                                                                           
print "channel = Combined, mu = %.2f, ErrLow = %.2f, ErrHigh = %.2f" % (rFit.getVal(),rFit.getAsymErrorLo(),rFit.getAsymErrorHi())

globalFitBand150 = ROOT.TBox(rFit150.getVal()+rFit150.getAsymErrorLo(), 0, rFit150.getVal()+rFit150.getAsymErrorHi(), nChann)                       
globalFitBand150.SetFillColor(ROOT.kGreen+2)                                                                                                       
globalFitBand150.SetLineStyle(0)                                                                                                        
globalFitBand150.Draw("SAME")                                                                                                           
globalFitLine150 = ROOT.TLine(rFit150.getVal(), 0, rFit150.getVal(), nChann)                                                                  
globalFitLine150.SetLineWidth(4)                                                                                                       
globalFitLine150.SetLineColor(ROOT.kBlack)                                                                                                    
globalFitLine150.Draw("SAME")                                                                                                           
points.Draw("PSAME")                                                                                                                
points150.Draw("PSAME")                                                                                                                
print "channel = Combined, mu = %.2f, ErrLow = %.2f, ErrHigh = %.2f" % (rFit150.getVal(),rFit150.getAsymErrorLo(),rFit150.getAsymErrorHi())
                                                                                                                                      
latex2 = ROOT.TLatex()                                                                                                                
latex2.SetNDC()                                                                                                                       
latex2.SetTextSize(0.5*c8.GetTopMargin())                                                                                              
latex2.SetTextFont(42)                                                                                                                
latex2.SetTextAlign(31) # align right                                                                                                 
latex2.DrawLatex(0.95, 0.95,"35.9 fb^{-1} (13 TeV)")                                                                                   
latex2.SetTextSize(0.7*c8.GetTopMargin())                                                                                              
latex2.SetTextFont(62)                                                                                                                
latex2.SetTextAlign(11) # align right                                                                                                 
latex2.DrawLatex(0.2, 0.95, "CMS")                                                                                                   
latex2.DrawLatex(0.3, 0.95, "preliminary")
latex2.SetTextSize(0.6*c8.GetTopMargin())                                                                                              
latex2.SetTextFont(52)                                                                                                                
latex2.SetTextAlign(11)                                                                                                               
#latex2.DrawLatex(0.27, 0.95, "Tutorial")                                                                                              

latex2.DrawLatex(0.65, 0.25, "#mu_{qqH}^{35.9fb^{-1}} = %.2f_{%.2f}^{+%.2f}" % (rFit.getVal(),rFit.getAsymErrorLo(),rFit.getAsymErrorHi()))
latex2.DrawLatex(0.65, 0.18, "#mu_{qqH}^{150fb^{-1}} = %.2f_{%.2f}^{+%.2f}" % (rFit150.getVal(),rFit150.getAsymErrorLo(),rFit150.getAsymErrorHi()))
                                                                                                                                      
                                                                                                                                      
ROOT.gPad.RedrawAxis()                                                                                                                
c8.Update()
c8.Draw()

raw_input()
