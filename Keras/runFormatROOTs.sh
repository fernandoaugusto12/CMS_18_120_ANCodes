#!/bin/bash
source ~/loadKeras.sh
cd /lustrehome/mmelodea/Keras/Trainings/FormatROOTs

python FormatROOTs.py --infile root_files.txt \
--outfile /lustre/cms/store/user/mmelodea/StudyVBF/KerasNtuples/hzz4l_vbf_selection_m4l118-130GeV_shuffledFS \
--njets 4 --maxsubjets 0 --nsubjets 0 \
--keys qqH ggH ggZZ ggZZ ggZZ ggZZ ggZZ ggZZ ggZZ ggZZ qqZZ qqZZ qqZZ ttH HWW HWW HWW TTV HWW WH WH ZH qqH3J \
--tags VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8 \
GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8 \
GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8 \
GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8 \
GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8 \
GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8 \
GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8 \
GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8 \
GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8 \
GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8 \
ZZJJTo4L_EWK_13TeV-madgraph-pythia8 \
ZZTo2L2Nu_13TeV_powheg_pythia8 \
ZZTo4L_13TeV_powheg_pythia8 \
ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUgenV6_pythia8 \
HWminusJ_HToWWTo2L2Nu_WToLNu_M125_13TeV_powheg_pythia8 \
HWplusJ_HToWWTo2L2Nu_WToLNu_M125_13TeV_powheg_pythia8 \
HZJ_HToWWTo2L2Nu_ZTo2L_M125_13TeV_powheg_pythia8 \
TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8 \
GluGluZH_HToWWTo2L2Nu_ZTo2L_M125_13TeV_powheg_pythia8 \
WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8 \
WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8 \
ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8 \
VBF_HJJJ_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8 | tee format_inputs.log