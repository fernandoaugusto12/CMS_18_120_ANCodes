###------------------------------------------------------------------------------------###
###  This code runs over ROOT files and creates a dictionary for the variables         ###
###  Its dictionary is used then by EvaluateNeuralNetwork.py                           ###
###  Author: Miqueias Melo de Almeida                                                  ###
###  Date: 08/06/2018                                                                  ###
###------------------------------------------------------------------------------------###

#set adequate environment
import os
import sys
import argparse
import cPickle as pickle
import ROOT
import math



def FormatInputs(list_of_files, output_file, tags, keys, Nevents, Njets, MaxSubjets, Nsubjets):
  frmap = ROOT.TFile.Open('/lustre/cms/store/user/mmelodea/ZplusX/ZplusX/FakeRatePtEtaMap.root')
  elefrw = frmap.Get('ElectronPtEtaMapFR')
  mufrw  = frmap.Get('MuonPtEtaMapFR')
  
  file_names = open(list_of_files,'r')
  files_addresses = []
  print 'Retrieving files addresses...'
  for i in file_names.readlines():
    files_addresses.append( i.rstrip()  )


  #format inputs from TTree (Data & MC)
  default_value = 0
  formated_inputs = []
  event_index = {}
  mc_sumweight = {}
  for ik in keys:
    event_index[ik] = 0
    mc_sumweight[ik] = 0    
        
  print 'Starting to load files...'
  ntags = len(tags)
  for ifile in range(len(files_addresses)):
    file_name = files_addresses[ifile]
    
    imc = ''
    for itag in range(ntags):
      if(file_name.find(tags[itag]) != -1):
	imc = keys[itag]
      else:
	continue
    if(imc == ''):
      continue
  
    tfile = ROOT.TFile.Open(file_name)
    print 'Processing %s' % tfile.GetName()
    tree_name = 'HZZ4LeptonsAnalysisReduced'
    if(file_name.find('ZplusX') != -1):
      tree_name = 'ControlRegions'
    tree = tfile.Get(tree_name)
    tree.GetEntries()
    for ievt, evt in enumerate(tree):
      #VBF category, -MELA cut, +m4l cut
      if( (((evt.f_njets_pass == 2 or evt.f_njets_pass == 3) and evt.f_Nbjets <= 1) or (evt.f_njets_pass > 3 and evt.f_Nbjets == 0)) and evt.f_mass4l >= 118 and evt.f_mass4l <= 130 ):
	#### use Njets to split in categories
	if(evt.f_njets_pass < Njets):
	  continue
            
        ##ZX special selection
        if(file_name.find('ZplusX') != -1 and evt.f_2p2p):
	  continue
        
        #gets the sum of weights
        event_index[imc] += 1
        if(file_name.find('ZplusX') == -1):
	  mc_sumweight[imc] += evt.f_weight
	else:
	  f1=0
	  f2=0
	  f3=0
	  f4=0
	  if(abs(f_lept1_pdgid)==11):
	    f1 = elefrw.GetBinContent( elefrw.FindBin(evt.f_lept1_pt,abs(evt.f_lept1_eta)) )
	  else:
	    f1 = mufrw.GetBinContent( mufrw.FindBin(evt.f_lept1_pt,abs(evt.f_lept1_eta)) )
	  if(abs(f_lept2_pdgid)==11):
	    f2 = elefrw.GetBinContent( elefrw.FindBin(evt.f_lept2_pt,abs(evt.f_lept2_eta)) )
	  else:
	    f2 = mufrw.GetBinContent( mufrw.FindBin(evt.f_lept2_pt,abs(evt.f_lept2_eta)) )
	  if(abs(f_lept3_pdgid)==11):
	    f3 = elefrw.GetBinContent( elefrw.FindBin(evt.f_lept3_pt,abs(evt.f_lept3_eta)) )
	  else:
	    f3 = mufrw.GetBinContent( mufrw.FindBin(evt.f_lept3_pt,abs(evt.f_lept3_eta)) )
	  if(abs(f_lept4_pdgid)==11):
	    f4 = elefrw.GetBinContent( elefrw.FindBin(evt.f_lept4_pt,abs(evt.f_lept4_eta)) )
	  else:
	    f4 = mufrw.GetBinContent( mufrw.FindBin(evt.f_lept4_pt,abs(evt.f_lept4_eta)) )
	    
	  fr_weight = 1
	  if not evt.f_lept1_pass:
	    fr_weight *= f1/(1-f1)
	  if not evt.f_lept2_pass:
	    fr_weight *= f2/(1-f2)
	  if not evt.f_lept3_pass:
	    fr_weight *= f3/(1-f3);
	  if not evt.f_lept4_pass:
	    fr_weight *= f4/(1-f4);
	  
	  mc_sumweight[imc] += fr_weight
	  
	if(event_index[imc] > Nevents):
          break;

        event = {}
        event['mc'] =  imc
        if(file_name.find('histos4mu') != -1):
            event['final_state'] =  '4mu'
        if(file_name.find('histos4e') != -1):
            event['final_state'] =  '4e'
        if(file_name.find('histos2e2mu') != -1):
            event['final_state'] =  '2e2mu'

	#event['f_outlier'] =  evt.f_outlier
	event['mc_sumweight'] = 0 
	event['f_weight'] =  evt.f_weight
	#event['f_Djet_VAJHU'] =  evt.f_Djet_VAJHU
	#event['f_Djet_VAJHU_UncUp'] =  evt.f_Djet_VAJHU_UncUp
	#event['f_Djet_VAJHU_UncDn'] =  evt.f_Djet_VAJHU_UncDn
	event['f_lept1_pt'] =  evt.f_lept1_pt
	event['f_lept1_pt_error'] =  evt.f_lept1_pt_error
	event['f_lept1_eta'] =  evt.f_lept1_eta
	event['f_lept1_phi'] =  evt.f_lept1_phi
	#event['f_lept1_charge'] =  evt.f_lept1_charge
	#event['f_lept1_pfx'] =  evt.f_lept1_pfx
	#event['f_lept1_sip'] =  evt.f_lept1_sip
	#event['f_lept1_pdgid'] =  evt.f_lept1_pdgid
	event['f_lept2_pt'] =  evt.f_lept2_pt
	event['f_lept2_pt_error'] =  evt.f_lept2_pt_error
	event['f_lept2_eta'] =  evt.f_lept2_eta
	event['f_lept2_phi'] =  evt.f_lept2_phi
	#event['f_lept2_charge'] =  evt.f_lept2_charge
	#event['f_lept2_pfx'] =  evt.f_lept2_pfx
	#event['f_lept2_sip'] =  evt.f_lept2_sip
	#event['f_lept2_pdgid'] =  evt.f_lept2_pdgid
	event['f_lept3_pt'] =  evt.f_lept3_pt
	event['f_lept3_pt_error'] =  evt.f_lept3_pt_error
	event['f_lept3_eta'] =  evt.f_lept3_eta
	event['f_lept3_phi'] =  evt.f_lept3_phi
	#event['f_lept3_charge'] =  evt.f_lept3_charge
	#event['f_lept3_pfx'] =  evt.f_lept3_pfx
	#event['f_lept3_sip'] =  evt.f_lept3_sip
	#event['f_lept3_pdgid'] =  evt.f_lept3_pdgid
	event['f_lept4_pt'] =  evt.f_lept4_pt
	event['f_lept4_pt_error'] =  evt.f_lept4_pt_error
	event['f_lept4_eta'] =  evt.f_lept4_eta
	event['f_lept4_phi'] =  evt.f_lept4_phi
	#event['f_lept4_charge'] =  evt.f_lept4_charge
	#event['f_lept4_pfx'] =  evt.f_lept4_pfx
	#event['f_lept4_sip'] =  evt.f_lept4_sip
	#event['f_lept4_pdgid'] =  evt.f_lept4_pdgid
	#event['f_iso_max'] =  evt.f_iso_max
	#event['f_sip_max'] =  evt.f_sip_max
	#event['f_Z1mass'] =  evt.f_Z1mass
	#event['f_Z2mass'] =  evt.f_Z2mass
	#event['f_angle_costhetastar'] =  evt.f_angle_costhetastar
	#event['f_angle_costheta1'] =  evt.f_angle_costheta1
	#event['f_angle_costheta2'] =  evt.f_angle_costheta2
	#event['f_angle_phi'] =  evt.f_angle_phi
	#event['f_angle_phistar1'] =  evt.f_angle_phistar1
	#event['f_pt4l'] =  evt.f_pt4l
	#event['f_eta4l'] =  evt.f_eta4l
	event['f_mass4l'] =  evt.f_mass4l
	event['f_njets_pass'] =  evt.f_njets_pass
	#event['f_deltajj'] =  evt.f_deltajj
	#event['f_massjj'] =  evt.f_massjj
	event['f_pfmet'] =  evt.f_pfmet
        #event['f_pfmet_eta'] =  -ROOT.TMath.Log(ROOT.TMath.Tan(evt.f_pfmet_theta/2.))
        #event['f_pfmet_phi'] =  evt.f_pfmet_phi
	event['f_pfmet_JetEnUp'] =  evt.f_pfmet_JetEnUp
	event['f_pfmet_JetEnDn'] =  evt.f_pfmet_JetEnDn
	event['f_pfmet_JetResUp'] =  evt.f_pfmet_JetResUp
	event['f_pfmet_JetResDn'] =  evt.f_pfmet_JetResDn
	event['f_pfmet_ElectronEnUp'] =  evt.f_pfmet_ElectronEnUp
	event['f_pfmet_ElectronEnDn'] =  evt.f_pfmet_ElectronEnDn
	event['f_pfmet_MuonEnUp'] =  evt.f_pfmet_MuonEnUp
	event['f_pfmet_MuonEnDn'] =  evt.f_pfmet_MuonEnDn
	event['f_pfmet_UnclusteredEnUp'] =  evt.f_pfmet_UnclusteredEnUp
	event['f_pfmet_UnclusteredEnDn'] =  evt.f_pfmet_UnclusteredEnDn
	event['f_pfmet_PhotonEnUp'] =  evt.f_pfmet_PhotonEnUp
	event['f_pfmet_PhotonEnDn'] =  evt.f_pfmet_PhotonEnDn
	#event['f_mT'] =  evt.f_mT
	#event['f_dphi'] =  evt.f_dphi
	event['f_nbjets'] =  evt.f_Nbjets
	for ijet in range(Njets):
	  check = evt.f_jets_highpt_pt[ijet]
	  #event['f_jets_highpt_btagger[%i]' % ijet] =  (evt.f_jets_highpt_btagger[ijet]  if (check != -999) else default_value)
	  event['f_jets_highpt_pt[%i]' % ijet] =  (evt.f_jets_highpt_pt[ijet]  if (check != -999) else default_value)
	  event['f_jets_highpt_pt_error[%i]' % ijet] =  (evt.f_jets_highpt_pt_error[ijet]  if (check != -999) else default_value)
	  event['f_jets_highpt_eta[%i]' % ijet] =  (evt.f_jets_highpt_eta[ijet]  if (check != -999) else default_value)
	  event['f_jets_highpt_phi[%i]' % ijet] =  (evt.f_jets_highpt_phi[ijet]  if (check != -999) else default_value)
	  jetE = evt.f_jets_highpt_et[ijet]*ROOT.TMath.CosH(evt.f_jets_highpt_eta[ijet])
	  event['f_jets_highpt_e[%i]' % ijet] =  (jetE if (check != -999) else default_value)
	  #event['f_jets_highpt_area[%i]' % ijet] =  (evt.f_jets_highpt_area[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_ptd[%i]' % ijet] =  (evt.f_jets_highpt_ptd[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_charged_hadron_energy[%i]' % ijet] =  (evt.f_jets_highpt_charged_hadron_energy[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_neutral_hadron_energy[%i]' % ijet] =  (evt.f_jets_highpt_neutral_hadron_energy[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_photon_energy[%i]' % ijet] =  (evt.f_jets_highpt_photon_energy[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_electron_energy[%i]' % ijet] =  (evt.f_jets_highpt_electron_energy[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_muon_energy[%i]' % ijet] =  (evt.f_jets_highpt_muon_energy[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_hf_hadron_energy[%i]' % ijet] =  (evt.f_jets_highpt_hf_hadron_energy[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_hf_em_energy[%i]' % ijet] =  (evt.f_jets_highpt_hf_em_energy[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_charged_em_energy[%i]' % ijet] =  (evt.f_jets_highpt_charged_em_energy[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_charged_mu_energy[%i]' % ijet] =  (evt.f_jets_highpt_charged_mu_energy[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_neutral_em_energy[%i]' % ijet] =  (evt.f_jets_highpt_neutral_em_energy[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_charged_hadron_multiplicity[%i]' % ijet] =  (evt.f_jets_highpt_charged_hadron_multiplicity[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_neutral_hadron_multiplicity[%i]' % ijet] =  (evt.f_jets_highpt_neutral_hadron_multiplicity[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_photon_multiplicity[%i]' % ijet] =  (evt.f_jets_highpt_photon_multiplicity[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_electron_multiplicity[%i]' % ijet] =  (evt.f_jets_highpt_electron_multiplicity[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_muon_multiplicity[%i]' % ijet] =  (evt.f_jets_highpt_muon_multiplicity[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_hf_hadron_multiplicity[%i]' % ijet] =  (evt.f_jets_highpt_hf_hadron_multiplicity[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_hf_em_multiplicity[%i]' % ijet] =  (evt.f_jets_highpt_hf_em_multiplicity[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_charged_multiplicity[%i]' % ijet] =  (evt.f_jets_highpt_charged_multiplicity[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_neutral_multiplicity[%i]' % ijet] =  (evt.f_jets_highpt_neutral_multiplicity[ijet]  if (check != -999) else default_value)
	  #event['f_jets_highpt_ncomponents[%i]' % ijet] =  (evt.f_jets_highpt_ncomponents[ijet]  if (check != -999) else default_value)
	  #### each of this guys will have 100 objects
	  #vecs = [[],[],[],[],[],[],[],[],[],[],[]]
	  #count_jets = 0
	  #for isjet in range(ijet*MaxSubjets,(ijet+1)*MaxSubjets):
	  #  check = evt.f_jets_highpt_component_pt[isjet]
	  #  vecs[0]  = (evt.f_jets_highpt_component_pdgid[isjet]       if (check != -999) else default_value)
	  #  vecs[1]  = (evt.f_jets_highpt_component_pt[isjet]          if (check != -999) else default_value)
	  #  vecs[2]  = (evt.f_jets_highpt_component_eta[isjet]         if (check != -999) else default_value)
	  #  vecs[3]  = (evt.f_jets_highpt_component_phi[isjet]         if (check != -999) else default_value)
	  #  vecs[4]  = (evt.f_jets_highpt_component_energy[isjet]      if (check != -999) else default_value)
	  #  vecs[5]  = (evt.f_jets_highpt_component_charge[isjet]      if (check != -999) else default_value)
	  #  vecs[6]  = (evt.f_jets_highpt_component_mt[isjet]          if (check != -999) else default_value)
	  #  vecs[7]  = (evt.f_jets_highpt_component_xvertex[isjet]     if (check != -999) else default_value)
	  #  vecs[8]  = (evt.f_jets_highpt_component_yvertex[isjet]     if (check != -999) else default_value)
	  #  vecs[9]  = (evt.f_jets_highpt_component_zvertex[isjet]     if (check != -999) else default_value)
	  #  vecs[10] = (evt.f_jets_highpt_component_vertex_chi2[isjet] if (check != -999) else default_value)
	  #  count_jets += 1
	  #  if(count_jets > Nsubjets):
	  #    break
	     
	  #event['f_jets_highpt_component_pdgid[%i]' % ijet] =  vecs[0]
	  #event['f_jets_highpt_component_pt[%i]' % ijet] =  vecs[1]
	  #event['f_jets_highpt_component_eta[%i]' % ijet] =  vecs[2]
	  #event['f_jets_highpt_component_phi[%i]' % ijet] =  vecs[3]
	  #event['f_jets_highpt_component_energy[%i]' % ijet] =  vecs[4]
	  #event['f_jets_highpt_component_charge[%i]' % ijet] =  vecs[5]
	  #event['f_jets_highpt_component_mt[%i]' % ijet] =  vecs[6]
	  #event['f_jets_highpt_component_xvertex[%i]' % ijet] =  vecs[7]
	  #event['f_jets_highpt_component_yvertex[%i]' % ijet] =  vecs[8]
	  #event['f_jets_highpt_component_zvertex[%i]' % ijet] =  vecs[9]
	  #event['f_jets_highpt_component_vertex_chi2[%i]' % ijet] =  vecs[10]
	  
	####---------- append the event dictionary ------
	formated_inputs.append( event )

    tfile.Close()
    #print ('%i - Processed (%i)' % (ifile, nacepted)), file_name
    
  print ''
  print '------------ Final Dataset Organization ----------------'
  print '--------------------- MCs ------------------------------'
  for ik in keys:
    print ik," -- nevents = ",event_index[ik]," -- sumweight = ",mc_sumweight[ik]
  print '----------------- Variables ----------------------------'
  print formated_inputs[0].keys()


  #for randomizing final state and processes in the dataset
  print 'Shuffling events...'
  import random
  for iev in range(len(formated_inputs)):
    #inserts the correct mc weight
    ik = formated_inputs[iev]['mc']
    formated_inputs[iev]['mc_sumweight'] = mc_sumweight[ik]

  #shuffle twice to get more randomization
  #set seed for reproductibility
  random.seed(999)
  random.shuffle(formated_inputs)
  random.seed(348)
  random.shuffle(formated_inputs)
        
  #save the events
  print ''
  print 'Saving events...'
  fileout = open(output_file+'.pkl','w')
  pickle.dump( formated_inputs, fileout )
  fileout.close()





def main(options):
  print '>>>>>>>>>> Inputs <<<<<<<<<<<<<'
  print options.infile
  print options.outfile
  print options.tags
  print options.keys
  print options.nevents
  print options.njets
  print options.maxsubjets
  print options.nsubjets
  
  if(len(options.keys)==len(options.tags)):
    print '>>>>>>>>>>  Assigments <<<<<<<<<<'
    for key,tag in zip(options.keys,options.tags):
      print key,' -- ',tag
    print ''
    FormatInputs(options.infile, options.outfile, options.tags, options.keys, options.nevents, options.njets, options.maxsubjets, options.nsubjets)
  else:
    print 'ERROR: Incompatible length of keys and tags'
  
  
  
if __name__ == '__main__':

 # Setup argument parser
 parser = argparse.ArgumentParser()

 # Add more arguments
 parser.add_argument("--infile", help="Name of txt input file with addresses of root files")
 parser.add_argument("--outfile", help="Name of the output file")
 parser.add_argument("--tags", nargs='+', help="Tags to find the MCs (from the rootfiles) to be used")
 parser.add_argument("--keys", nargs='+', help="Tags for the MCs to be used inside Keras framework")
 parser.add_argument("--nevents", type=int, default=10000000, help="Number of events to get from each file")
 parser.add_argument("--njets", type=int, default=4, help="Number of jets to be picked")
 parser.add_argument("--maxsubjets", type=int, default=10, help="Number of jet components to be picked")
 parser.add_argument("--nsubjets", type=int, default=10, help="Number of jet components to be picked")

 # Parse default arguments
 options = parser.parse_args()
 main(options)
    
