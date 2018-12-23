import sys, os
import argparse
import copy
import math
import ROOT
import numpy as np
from matplotlib import pyplot as pyp


##----- define a library of available shape uncertainties ---------
ushapes = {'':[],
	   'l1pt':['_l1ptDown','_l1ptUp'],
	   'l2pt':['_l2ptDown','_l2ptUp'],
	   'l3pt':['_l3ptDown','_l3ptUp'],
	   'l4pt':['_l4ptDown','_l4ptUp'],
	   'j1pt':['_j1ptDown','_j1ptUp'],
	   'j2pt':['_j2ptDown','_j2ptUp'],
	   'j3pt':['_j3ptDown','_j3ptUp'],
	   'j4pt':['_j4ptDown','_j4ptUp'],
	   'met' :['_metEleEnDown','_metEleEnUp','_metJetEnDown','_metJetEnUp','_metJetResDown','_metJetResUp',
		   '_metMuEnDown','_metMuEnUp','_metPhoEnDown','_metPhoEnUp','_metUncEnDown','_metUncEnUp'],
	   'FR'  :['_FRUp','_FRDown'],
	   'j3kShape':['_j3kShapeUp','_j3kShapeDown']
}

#---------------- function containing standard parts for Higgs datacards --------------------------
def cardTemplate(options, yields, root_card_name, ich):
  lines = []
    
  lines.append('#### Datacard for Higgs Combine Tool')
  lines.append('#### Source: %s' % options.infile)
  lines.append('---------------------------------------------------------------------------------')
  lines.append('imax\t\t\t\t1  number of channels') #the program expects only one signal
  lines.append('jmax\t\t\t\t%i  number of backgrounds' % (len(options.keys)-2)) #don't count data and signal
  lines.append('kmax\t\t\t\t*  number of nuisance parameters (sources of systematical uncertainties)')
  lines.append('---------------------------------------------------------------------------------')
  lines.append('shapes\t\t\t\t* * %s $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC' % root_card_name)
  lines.append('---------------------------------------------------------------------------------')
  lines.append('bin\t\t\t\t%s' % ich)
  lines.append('observation\t\t\t%i' % yields['data_obs'])
  lines.append('---------------------------------------------------------------------------------')

  line = 'bin\t\t\t'
  processes = [ik for ik in options.keys if ik != 'data_obs']
  for ikey in processes:
    line += '\t%s\t' % ich
  lines.append(line)
  
  line = 'process\t\t\t'
  for ikey in processes:
    line += '\t%s\t' % ikey
  lines.append(line)
  
  line = 'process\t\t\t'
  i = 0
  for ikey in processes:
    line += '\t%s\t' % i
    i += 1
  lines.append(line)

  line = 'rate\t\t\t'
  for ikey in processes:
    line += '\t%.4f\t' % yields[ikey]
  lines.append(line)
  
  lines.append('---------------------------------------------------------------------------------')
  line = 'lumi_13TeV\t\tlnN'
  for ikey in processes:
    if(ikey != 'zjets'):
      line += '\t%s\t' % '1.026'
    else:
      line += '\t-\t'
  lines.append(line)
  
  line = 'ID_reco_eff\t\tlnN'
  for ikey in processes:
    if(ikey != 'zjets'):
      if('4mu' in ich):
	line += '\t%s\t' % '1.025'
      if('4e' in ich):
	line += '\t%s\t' % '1.09'
      if('2e2mu' in ich):
	line += '\t%s\t' % '1.057'
    else:
      line += '\t-\t'
  lines.append(line)
  
  if('4mu' in ich):
    line = 'mu_en_scale\t\tlnN'
    for ikey in processes:
      if(ikey != 'zjets'):
	line += '\t%s\t' % '1.0004'
      else:
	line += '\t-\t'
    lines.append(line)
  
  if('4e' in ich):
    line = 'ele_en_scale\t\tlnN'
    for ikey in processes:
      if(ikey != 'zjets'):
	line += '\t%s\t' % '1.003'
      else:
	line += '\t-\t'
    lines.append(line)

  if('2e2mu' in ich):
    line = 'mu_ele_en_scale\t\tlnN'
    for ikey in processes:
      if(ikey != 'zjets'):      
	line += '\t%s\t' % '1.0017'
      else:
	line += '\t-\t'
    lines.append(line)
  
  line = 'JES\t\t\tlnN'
  JES = {'WH':'1.0086/0.9880','ZH':'1.0186/0.9840','ggH':'1.0141/0.9849','qqH':'1.0334/0.9775','ttH':'1.0312/0.9726','HWW':'1.0141/0.9849','ggZZ':'1.0180/0.9799','qqZZ':'1.0078/0.9869','VVV':'1.0078/0.9869','TTV':'1.0312/0.9726','zjets':'-'}
  for ikey in processes:
    line += '\t%s' % JES[ikey]
  lines.append(line)
  
  line = 'mass_resol\t\tlnN'
  for ikey in processes:
    if(ikey=='WH' or ikey == 'ZH' or ikey=='ggH' or ikey=='qqH' or ikey=='ttH' or ikey=='HWW'):
      line += '\t%s\t' % '1.20'
    else:
      line += '\t-\t'
  lines.append(line)

  line = 'b_tag\t\t\tlnN'
  for ikey in processes:
    if(ikey != 'zjets'):    
      line += '\t%s\t' % '1.01'
    else:
      line += '\t-\t'
  lines.append(line)
  
  line = 'CMS_trig\t\tlnN'
  for ikey in processes:
    if(ikey != 'zjets'):
      line += '\t%s\t' % '1.02'
    else:
      line += '\t-\t'
  lines.append(line)

  line = 'BRhiggs_hzz4l\t\tlnN'
  for ikey in processes:
    if(ikey=='WH' or ikey == 'ZH' or ikey=='ggH' or ikey=='qqH' or ikey=='ttH' or ikey=='HWW'):
      line += '\t%s\t' % '1.02'
    else:
      line += '\t-\t'
  lines.append(line)

  line = 'QCDscale_qqH\t\tlnN'
  for ikey in processes:
    if(ikey=='qqH'):
      line += '\t%s' % '1.004/0.997'
    else:
      line += '\t-\t'
  lines.append(line)
  
  line = 'QCDscale_ZH\t\tlnN'
  for ikey in processes:
    if(ikey=='ZH'):
      line += '\t%s' % '1.038/0.969'
    else:
      line += '\t-\t'
  lines.append(line)
  
  line = 'QCDscale_WH\t\tlnN'
  for ikey in processes:
    if(ikey=='WH'):
      line += '\t%s' % '1.005/0.993'
    else:
      line += '\t-\t'
  lines.append(line)
  
  line = 'QCDscale_ttH\t\tlnN'
  for ikey in processes:
    if(ikey=='ttH'):
      line += '\t%s' % '1.058/0.908'
    else:
      line += '\t-\t'
  lines.append(line)

  line = 'QCDscale_ggH\t\tlnN'
  for ikey in processes:
    if(ikey=='ggH' or ikey=='ggZZ'):
      line += '\t%s' % '1.039/0.961'
    else:
      line += '\t-\t'
  lines.append(line)

  line = 'QCDscale_qqZZ\t\tlnN'
  for ikey in processes:
    if(ikey=='qqZZ'):
      line += '\t%s' % '1.032/0.958'
    else:
      line += '\t-\t'
  lines.append(line)

  #line = 'QCDscale_TTV\t\tlnN'
  #for ikey in processes:
  #  if(ikey=='TTV'):
  #    line += '\t%s' % '0.846/1.146'
  #  else:
  #    line += '\t-\t'
  #lines.append(line)
  
  line = 'QCDscale_ggVV_bonly\tlnN'
  for ikey in processes:
    if(ikey=='ggZZ'):
      line += '\t%s\t' % '1.1'
    else:
      line += '\t-\t'
  lines.append(line)

  line = 'pdf_qqZZ\t\tlnN'
  for ikey in processes:
    if(ikey=='qqZZ'):
      line += '\t%s' % '1.031/0.966'
    else:
      line += '\t-\t'
  lines.append(line)

  line = 'pdf_Higgs_ttH\t\tlnN'
  for ikey in processes:
    if(ikey=='ttH'):
      line += '\t%s' % '1.036/0.964'
    else:
      line += '\t-\t'
  lines.append(line)

  line = 'pdf_Higgs_qq\t\tlnN'
  for ikey in processes:
    if(ikey=='WH'):
      line += '\t%s' % '1.019/0.981'
    elif(ikey=='ZH'):
      line += '\t%s' % '1.016/0.984'
    elif(ikey=='qqH'):
      line += '\t%s' % '1.021/0.979'
    else:
      line += '\t-\t'
  lines.append(line)

  line = 'pdf_Higgs_gg\t\tlnN'
  for ikey in processes:
    if(ikey=='ggH' or ikey=='ggZZ'):
      line += '\t%s' % '1.032/0.968'
    else:
      line += '\t-\t'
  lines.append(line)
  
  #line = 'pdf_TTV\t\t\tlnN'
  #for ikey in processes:
  #  if(ikey=='TTV'):
  #    line += '\t%s' % '0.974/1.025'
  #  else:
  #    line += '\t-\t'
  #lines.append(line)

  line = 'EWcorr_qqZZ\t\tlnN'
  for ikey in processes:
    if(ikey=='qqZZ'):
      line += '\t%s' % '1.0012/0.998798'
    else:
      line += '\t-\t'
  lines.append(line)

  line = 'FR\t\t\tlnN'
  for ikey in processes:
    if(ikey == 'zjets'):
      if('4mu' in ich):
	line += '\t%s\t' % '1.35'
      if('4e' in ich):
	line += '\t%s\t' % '1.32'
      if('2e2mu' in ich):
	line += '\t%s\t' % '1.33'
    else:
	line += '\t-\t'
  lines.append(line)
  
  #suspect line... it might be removed
  lines.append('CMS_zz4l_bkg_kdShape param 0.0 1 [-3,3]')
  
  #shape uncertainties for NNs
  card_shapes = copy.deepcopy(ushapes)
  for ishape in options.shapes:
    for iushape in range(len(ushapes[ishape])):
      if(card_shapes[ishape][iushape].find('Down') != -1):
	card_shapes[ishape][iushape] = card_shapes[ishape][iushape].replace('Down','')
	card_shapes[ishape][iushape] = card_shapes[ishape][iushape].replace('_','')
      if(card_shapes[ishape][iushape].find('Up') != -1):
	card_shapes[ishape][iushape] = card_shapes[ishape][iushape].replace('Up','')
	card_shapes[ishape][iushape] = card_shapes[ishape][iushape].replace('_','')
  for ishape in options.shapes:
    for iushape in card_shapes[ishape]:
      if(line.find(iushape) == -1):
	if((iushape.find('met') != -1 and iushape.find('Mu') == -1) or ishape == 'j3kShape' ):
	  line = '%s\t\tshape' % iushape
	else:
	  line = '%s\t\t\tshape' % iushape
	for ikey in processes:
	  if( (ikey != 'zjets' and ishape == 'FR') or (ikey != 'qqH' and ishape == 'j3kShape') or (ikey == 'zjets' and ishape != 'FR')):
	    line += '\t-\t'
	  else:
	    line += '\t%s\t' % '1.000' #1sigma uncertainty (0.5=0.5sigma, 2=2sigma, etc)
	lines.append(line)

  return lines
##-----------------------------------------------------------------------------------

  

def main(options):
  print '\n\n>>> Generating Datacard for Higgs Combine Tool <<<'
  print 'Reading from: %s' % options.infile
  print 'Channels: %s' % options.channels
  print 'Keys: %s' % options.keys
  print 'Shape uncertainties: %s' % options.shapes
  print '---------------- Outputing------------------------'
  print 'Datacard: datacard_%s.txt' % options.outfile
  root_card_name = 'datacard_%s.root' % options.outfile
  print 'Shape file: %s' % root_card_name
  print '--------------------------------------------------'
  
  #first checks if there's no empty bins
  print 'Checking quality of histograms...'  
  fin = ROOT.TFile(options.infile)
  for ich in options.channels:
    bbins = []
    dbins = []    
    for ikey in options.keys:
      if(ikey != 'data_obs') and (ikey not in options.signal):
	h = fin.Get('%s/%s' % (ich,ikey))
	h.Rebin(options.rebin)
	if(len(bbins) == 0):
	  for ibin in range(h.GetNbinsX()):
	    bbins.append(h.GetBinContent(ibin+1))
	else:
	  for ibin in range(h.GetNbinsX()):
	    bbins[ibin] += h.GetBinContent(ibin+1)
      if(ikey == 'data_obs'):
	h = fin.Get('%s/%s' % (ich,ikey))
	h.Rebin(options.rebin)
	if(len(dbins) == 0):
	  for ibin in range(h.GetNbinsX()):
	    dbins.append(h.GetBinContent(ibin+1))
	else:
	  for ibin in range(h.GetNbinsX()):
	    dbins[ibin] += h.GetBinContent(ibin+1)
	    
    if(0 in bbins):
      #for bc, dc in zip(bbins,dbins):
	#if(bc == 0 and dc != 0):
	  print '>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<'
	  print '>>>>>>> Bad histograms... empty bin exists for %s! <<<<<<' % ich
	  print '>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<'
  fin.Close()
  
  
  
  for ich in options.channels:
    print "Producing datacards for channel %s..." % ich
    #------ extract the content from each histogram ----
    #------- and combine them to build up datacard -----
    histos = {}
    nentries = {}
    yields = {}
    errors = {}
    d_yield = 0
    d_error = 0
    s_yield = 0
    s_error = 0
    b_yield = 0
    b_error = 0
    fin = ROOT.TFile(options.infile)
    for ikey in options.keys:
      histos[ikey] = {}
      
      #nominal histogram
      htmp = fin.Get('%s/%s' % (ich,ikey))
      htmp.Rebin(options.rebin)
      #to cut distribution
      if(options.cutatbin != 0):
	histos[ikey][ikey] = ROOT.TH1D('%s' % ikey,'%s' % ikey,htmp.GetNbinsX()-options.cutatbin+1,htmp.GetBinLowEdge(options.cutatbin),1.0)
      else:
	histos[ikey][ikey] = ROOT.TH1D('%s' % ikey,'%s' % ikey,htmp.GetNbinsX(),0.0,1.0)
      jbin = 1
      for ibin in range(1,htmp.GetNbinsX()+1):
	if(ibin < options.cutatbin):
	  continue
	histos[ikey][ikey].SetBinContent(jbin,htmp.GetBinContent(ibin))
	histos[ikey][ikey].SetBinError(jbin,htmp.GetBinError(ibin))
	jbin += 1

      #+/-1 sigma shifts
      for ishape in options.shapes:
	for iushape in ushapes[ishape]:
	  if( (ikey != 'zjets' and ishape == "FR") or (ikey != 'qqH' and ishape == 'j3kShape') ):
	    continue
	  htmp = fin.Get('%s/%s%s' % (ich,ikey,iushape))
	  htmp.Rebin(options.rebin)
	  #to cut distribution
	  if(options.cutatbin != 0):
	    histos[ikey][ikey+iushape] = ROOT.TH1D('%s%s' % (ikey,iushape),'%s%s' % (ikey,iushape),htmp.GetNbinsX()-options.cutatbin+1,htmp.GetBinLowEdge(options.cutatbin),1.0)
	  else:
	    histos[ikey][ikey+iushape] = ROOT.TH1D('%s%s' % (ikey,iushape),'%s%s' % (ikey,iushape),htmp.GetNbinsX(),0.0,1.0)
	  jbin = 1
	  for ibin in range(1,htmp.GetNbinsX()+1):
	    if(ibin < options.cutatbin):
	      continue	    
	    histos[ikey][ikey+iushape].SetBinContent(jbin,htmp.GetBinContent(ibin))
	    histos[ikey][ikey+iushape].SetBinError(jbin,htmp.GetBinError(ibin))
	    jbin += 1

      #get rates to display in the datacard
      nentries[ikey] = histos[ikey][ikey].GetEntries()
      yields[ikey] = histos[ikey][ikey].Integral()
      errors[ikey] = 0
      error = 0
      for ibin in range(1, histos[ikey][ikey].GetNbinsX()+1):
	error += histos[ikey][ikey].GetBinError(ibin)**2
      errors[ikey] = math.sqrt(error)
        
      print 'Process %s, nentries: %i, yields: %.4f, error: %.4f' % (ikey,nentries[ikey],yields[ikey],errors[ikey])
      if(ikey == 'data_obs'):
	d_yield += yields[ikey]
	d_error += error
      elif(ikey == options.signal):
	s_yield += yields[ikey]
	s_error += error
      else:
	b_yield += yields[ikey]
	b_error += error
      
  
    #---- print total sig/bkg yields and uncertainties------
    print '\n------------------- Summary ---------------------'
    line  = '  Keys: '
    line1 = 'Yields: '
    line2 = 'Errors: '
    for ikey in options.keys:
      if(ikey == 'data_obs'):
	line += '%s\t' % 'Data'
      else:
	line += '%s\t' % ikey
    for ikey in options.keys:
      line1 += '%.4f\t' % yields[ikey]
    for ikey in options.keys:
      line2 += '%.4f\t' % errors[ikey]
   
    print line    
    print line1
    print line2
    print '--------------------------------------------------'
    print 'Data = %.4f +/- %.4f' % (d_yield, math.sqrt(d_error))
    print 'SIG  = %.4f +/- %.4f' % (s_yield, math.sqrt(s_error))
    print 'BKG  = %.4f +/- %.4f' % (b_yield, math.sqrt(b_error))
  

    #------ generates the datacard based on template --------
    if(options.cutatbin == 0):
      cardoutname = 'datacard_%s_%s.txt' % (options.outfile, ich)
    else:
      cardoutname = 'datacard_%s_%s_cutbin%s.txt' % (options.outfile, ich, options.cutatbin)
    print 'Datacard: %s' % cardoutname
    cardout = open(cardoutname,'w')
    lines = cardTemplate(options, yields, cardoutname.replace('.txt','.root'), ich)
    for iline in range(len(lines)):
      if(iline < len(lines)-1):
	cardout.write(lines[iline]+'\n')
      else:
	cardout.write(lines[iline])
    cardout.close()
  
    #------ generates the root file with the shapes and uncertainties -------
    fout = ROOT.TFile(cardoutname.replace('.txt','.root'), 'recreate')
    print 'Root file: %s' % fout.GetName()
    fout.mkdir(ich)
    fout.cd(ich)
    for ikey in options.keys:
      histos[ikey][ikey].Write()
      if(ikey != 'data_obs'):
	for ishape in options.shapes:
          if( (ikey != 'zjets' and ishape == "FR") or (ikey != 'qqH' and ishape == 'j3kShape') ):
            continue
	  for iushape in ushapes[ishape]:
	    #print 'Writing histos[%s][%s]' % (ikey,ikey+iushape)
	    histos[ikey][ikey+iushape].Write()
    fout.Close()
    fin.Close()
  


#-----------------------------------------------------------  
if __name__ == '__main__':

 # Setup argument parser
 parser = argparse.ArgumentParser()

 # Add more arguments
 parser.add_argument("--infile", help="Name of root file containing histograms for combine")
 parser.add_argument("--outfile", help="Name to be used for the datacard")
 parser.add_argument("--channels", nargs='+', help="Channels name")
 parser.add_argument("--keys", nargs='+', help="Histogram names (all them, even data - note that the progam expects <data_obs> as data histogram name)")
 parser.add_argument("--signal", help="Signal process name")
 parser.add_argument("--shapes", nargs='+', help="Shape uncertainties")
 parser.add_argument("--rebin", type=int, default=1, help="Option to rebin histograms")
 parser.add_argument("--cutatbin", type=int, default=0, help="Option to remove bins < cutatbin")
 
 

 # Parse default arguments
 options = parser.parse_args()
 main(options)
