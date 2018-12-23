import os.path
import math
import matplotlib.pyplot as pyp
import numpy as np
import cPickle as pickle

#filter results... ROC auc diff between Train and Test must be at most 1%
max_auc_diff = 1
njobs = 1920
for njet in [2,3]:
  vresults = {}
  points = []
  highest_ep = [0,0,0]
  nconcluded = 0
  print 'Checking Njets%i...' % njet
  for imodel in range(njobs):
    fin_name = 'Njets%i/Model%i/results/SummaryOfResults.pkl' % (njet,imodel)
    if not os.path.exists(fin_name):
      continue
    
    nconcluded += 1
    fin = open(fin_name, 'r')
    results = pickle.load( fin )
    fin.close()
    #print fin_name


    ###------------ set the dictionary --------------------------
    if len(vresults) == 0:
      vresults['inputs'] = []
      vresults['topology'] = []
      vresults['neuron'] = []
      vresults['minimizer'] = []
      vresults['patience'] = []
      vresults['batchsize'] = []
      vresults['scaletrain'] = []
      vresults['nooutliers'] = []
      vresults['taccuracy'] = []
      vresults['vaccuracy'] = []
      vresults['nn_roctrain'] = []
      vresults['nn_roctest'] = []
      vresults['nn_trainep'] = []
      vresults['nn_testep'] = []
      vresults['nn_trainsb'] = []
      vresults['nn_testsb'] = []
      for ik in results['nn']:
	if 'fullroc' in ik:
	  vresults['nn_%s' % ik] = []
      vresults['mela_roctrain'] = []
      vresults['mela_roctest'] = []
      vresults['mela_trainep'] = []
      vresults['mela_testep'] = []
      vresults['mela_trainsb'] = []
      vresults['mela_testsb'] = [] 
      for ik in results['mela']:
	if 'fullroc' in ik:
	  vresults['mela_%s' % ik] = []
    ###-------------------------------------------------------
    if njet == 3:
      max_auc_diff = 2
    if 100*math.fabs(results['nn']['roctrain']-results['nn']['roctest']) > max_auc_diff:
      continue
  
    points.append( imodel )
    #### store usefull variables for plotting
    vresults['inputs'].append( results['nninputs'] )
    vresults['topology'].append( results['topology'] )
    vresults['neuron'].append( results['neuron'] )
    vresults['minimizer'].append( results['minimizer'] )
    vresults['patience'].append( results['patience'] )
    vresults['batchsize'].append( results['batchsize'] )
    vresults['scaletrain'].append( results['scaletrain'] )
    vresults['nooutliers'].append( results['nooutliers'] )
    vresults['taccuracy'].append( results['accuracy']['train'] )
    vresults['vaccuracy'].append( results['accuracy']['test'] )
  
    vresults['nn_roctrain'].append( results['nn']['roctrain'] )
    vresults['nn_roctest'].append( results['nn']['roctest'] )
    vresults['nn_trainep'].append( results['nn']['trainep'] )
    vresults['nn_testep'].append( results['nn']['testep'] )
    vresults['nn_trainsb'].append( results['nn']['trainsb'] )
    vresults['nn_testsb'].append( results['nn']['testsb'] )
    for ik in results['nn']:
      if 'fullroc' in ik:
	vresults['nn_%s' % ik].append( results['nn'][ik] )

    vresults['mela_roctrain'].append( results['mela']['roctrain'] )
    vresults['mela_roctest'].append( results['mela']['roctest'] )
    vresults['mela_trainep'].append( results['mela']['trainep'] )
    vresults['mela_testep'].append( results['mela']['testep'] )
    vresults['mela_trainsb'].append( results['mela']['trainsb'] )
    vresults['mela_testsb'].append( results['mela']['testsb'] )
    for ik in results['mela']:
      if 'fullroc' in ik:
	vresults['mela_%s' % ik].append( results['mela'][ik] )


  ####----------------------------
  print 'Jobs eff: %f' % (nconcluded/float(njobs))

  ####------ make the plots ---------------------------------
  colors = ['blue','red','gold','darkgreen','gray','magenta']

  #### eS*p, s/sqrt(B) and ROC AUC from test/train
  print '--------- best configs based on eSp ---------------'
  for parameter in ['inputs','topology','neuron','minimizer','scaletrain','nooutliers']:
    print 'Plotting for %s' % parameter
    for iset in ['test']:
      icolor = 0
      uinputs = []
      values = {}
      for irun in range(len(vresults['inputs'])):
	if(vresults[parameter][irun] not in uinputs):
	  uinputs.append(vresults[parameter][irun])
	  values['points_%s' % vresults[parameter][irun]] = []
	  values['dep_%s' % vresults[parameter][irun]] = []
	  values['dsb_%s' % vresults[parameter][irun]] = []
	deltaep = 100*(vresults['nn_%sep' % iset][irun]-vresults['mela_%sep' % iset][irun])
	deltasb = 100*(vresults['nn_%ssb' % iset][irun]-vresults['mela_%ssb' % iset][irun])
	values['points_%s' % vresults[parameter][irun]].append(points[irun])
	values['dep_%s' % vresults[parameter][irun]].append(deltaep)
	values['dsb_%s' % vresults[parameter][irun]].append(deltasb)
	if(parameter == 'inputs' and iset == 'test' and deltaep >= highest_ep[0]):
	  highest_ep[0] = deltaep
	  highest_ep[1] = vresults['nn_%sep' % iset][irun]
	  highest_ep[2] = points[irun]
	  ###------ print best configs ------
	  if(iset == 'test'):
	    print highest_ep
      
      pyp.rc("font", size=18)
      pyp.figure(figsize=(10,7))
      pyp.subplot(211)
      for cat in uinputs:
	catname = cat
	if(parameter == 'nooutliers'):
	  if(cat):
	    catname = 'No Outliers'
	  else:
	    catname = 'With Outliers'
	if(parameter != 'nooutliers') and ('f_lept1_pt' in cat):
	  if('f_pfmet' in cat):
	    catname = r'$4l,%ij(p_{T},\eta,\phi),MET$' % njet
	  else:
	    catname = r'$4l,%ij(p_{T},\eta,\phi)$' % njet
	
	pyp.plot(values['points_%s' % cat],values['dep_%s' % cat],'*',color=colors[uinputs.index(cat)],label='%s' % catname)
      pyp.ylim(ymin=0)
      #pyp.ylim(ymax=highest_ep[0]*1.2)
      pyp.grid(True)
      pyp.legend()
      pyp.title(r'$[NN-MELA_{VBF}](\%)$')
      pyp.ylabel(r'$\epsilon_S.\pi$')
      #------------------------------
      pyp.subplot(212)
      for cat in uinputs:
	pyp.plot(values['points_%s' % cat],values['dsb_%s' % cat],'*',color=colors[uinputs.index(cat)],label='%s' % catname)
      pyp.ylim(ymin=0)
      #pyp.ylim(ymax=highest_sb[0]*1.2)
      pyp.grid(True)
      pyp.ylabel(r'$S/\sqrt{B}$')
      pyp.xlabel('Configurations')
      pyp.savefig('SummaryOfResults_Metrics_Njets%i_%s_%s.png' % (njet,iset,parameter))
