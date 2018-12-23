####recommended Method from LHC Higgs Combination Group
methods=("FitDiagnostics --plots --saveOverallShapes --saveWithUncertainties"
         "AsymptoticLimits"
         "Significance"
         "Significance -t -1 --expectSignal 1"
         "HybridNew --LHCmode LHC-limits -H AsymptoticLimits"
         "HybridNew --LHCmode LHC-limits -H AsymptoticLimits -T 2000 --expectedFromGrid=0.025 --expectSignal=1"
         "HybridNew --LHCmode LHC-limits -H AsymptoticLimits -T 2000 --expectedFromGrid=0.16 --expectSignal=1"
         "HybridNew --LHCmode LHC-limits -H AsymptoticLimits -T 2000 --expectedFromGrid=0.5 --expectSignal=1"
         "HybridNew --LHCmode LHC-limits -H AsymptoticLimits -T 2000 --expectedFromGrid=0.84 --expectSignal=1"
         "HybridNew --LHCmode LHC-limits -H AsymptoticLimits -T 2000 --expectedFromGrid=0.975 --expectSignal=1"
         "HybridNew --LHCmode LHC-significance -H Significance"
         "HybridNew --LHCmode LHC-significance -H Significance -T 2000 --expectedFromGrid=0.025 --expectSignal=1"
         "HybridNew --LHCmode LHC-significance -H Significance -T 2000 --expectedFromGrid=0.16 --expectSignal=1"
         "HybridNew --LHCmode LHC-significance -H Significance -T 2000 --expectedFromGrid=0.5 --expectSignal=1"
         "HybridNew --LHCmode LHC-significance -H Significance -T 2000 --expectedFromGrid=0.84 --expectSignal=1"
         "HybridNew --LHCmode LHC-significance -H Significance -T 2000 --expectedFromGrid=0.975 --expectSignal=1"
)

folder=("BestFit"
        "AsymLimit"
        "ObsAsymSignificance"
        "ExpAsymSignificance"
        "ObsLimHybridNew"
        "Exp2SdownLimHybridNew"
        "Exp1SdownLimHybridNew"
        "ExpMedianLimHybridNew"
        "Exp1SupLimHybridNew"
        "Exp2SupLimHybridNew"
        "ObsSigHybridNew"
        "Exp2SdownSigHybridNew"
        "Exp1SdownSigHybridNew"
        "ExpMedianSigHybridNew"
        "Exp1SupSigHybridNew"
        "Exp2SupSigHybridNew"
)

models=(#"MELA" 
        "NN1" 
        "NN2" 
        "NN3" 
        "NN6" 
        "NN10" 
        "NN11" 
        "NN372_nj2_bs1" 
        "NN412_nj2_bs1" 
        "NN1098_nj2" 
        "NN1137_nj2" 
        "NN1177_nj2"
        "NN11_nj2" 
        "NN131_nj2" 
        "NN138_nj2" 
        "NN211_nj2" 
        "NN51_nj2" 
        "NN91_nj2"
)


for model in "${models[@]}"
do
  scheme=0
  for m in "${methods[@]}"
  do
    mkdir -p ${model}/${folder[${scheme}]}
    #cp ../*${model}_* ${model}/
    #creates the condor config file
    echo "universe = vanilla" > ${model}/${folder[${scheme}]}/condor_${model}_${folder[${scheme}]}.cfg
    echo "Executable = /lustrehome/mmelodea/Keras/Uncertainty/Sigma1/Datacards/CombineLimits/${model}/${folder[${scheme}]}/${model}_${folder[${scheme}]}.sh" >> ${model}/${folder[${scheme}]}/condor_${model}_${folder[${scheme}]}.cfg
    echo "Should_Transfer_Files = YES" >> ${model}/${folder[${scheme}]}/condor_${model}_${folder[${scheme}]}.cfg
    echo "WhenToTransferOutput = ON_EXIT" >> ${model}/${folder[${scheme}]}/condor_${model}_${folder[${scheme}]}.cfg
    echo 'Requirements = TARGET.OpSys == "LINUX"&& (TARGET.Arch != "DUMMY" )' >> ${model}/${folder[${scheme}]}/condor_${model}_${folder[${scheme}]}.cfg
    echo "notify_user = mmelodea@cern.ch" >> ${model}/${folder[${scheme}]}/condor_${model}_${folder[${scheme}]}.cfg
    echo "Queue 1" >> ${model}/${folder[${scheme}]}/condor_${model}_${folder[${scheme}]}.cfg

    #creates the executable shell
    echo '#!/bin/bash' > ${model}/${folder[${scheme}]}/${model}_${folder[${scheme}]}.sh
    echo "#${m}" >> ${model}/${folder[${scheme}]}/${model}_${folder[${scheme}]}.sh
    echo "source /lustrehome/mmelodea/logincmspos.sh" >> ${model}/${folder[${scheme}]}/${model}_${folder[${scheme}]}.sh
    echo "cd /lustrehome/mmelodea/CMSSW_8_1_0/src" >> ${model}/${folder[${scheme}]}/${model}_${folder[${scheme}]}.sh
    echo 'eval `scram runtime -sh`' >> ${model}/${folder[${scheme}]}/${model}_${folder[${scheme}]}.sh
    echo "cd /lustrehome/mmelodea/Keras/Uncertainty/Sigma1/Datacards/CombineLimits/${model}/${folder[${scheme}]}" >> ${model}/${folder[${scheme}]}/${model}_${folder[${scheme}]}.sh
    echo "combine /lustrehome/mmelodea/Keras/Uncertainty/Sigma1/Datacards/combinedCard_Higgs125GeV_VBFStudy_${model}_4l.txt -M ${m} > combineResults.log 2>&1" >> ${model}/${folder[${scheme}]}/${model}_${folder[${scheme}]}.sh
    echo "combine ${model} -M ${m}"
    condor_submit -name ettore ${model}/${folder[${scheme}]}/condor_${model}_${folder[${scheme}]}.cfg
    scheme=$((scheme+1))
  done

  ##################### Impact Plot ########################
  newfolder="ImpactPlot"
  mkdir -p ${model}/${newfolder}
  echo "universe = vanilla" > ${model}/${newfolder}/condor_${model}_${newfolder}.cfg
  echo "Executable = /lustrehome/mmelodea/Keras/Uncertainty/Sigma1/Datacards/CombineLimits/${model}/${newfolder}/${model}_${newfolder}.sh" >> ${model}/${newfolder}/condor_${model}_${newfolder}.cfg
  echo "Should_Transfer_Files = YES" >> ${model}/${newfolder}/condor_${model}_${newfolder}.cfg
  echo "WhenToTransferOutput = ON_EXIT" >> ${model}/${newfolder}/condor_${model}_${newfolder}.cfg
  echo 'Requirements = TARGET.OpSys == "LINUX"&& (TARGET.Arch != "DUMMY" )' >> ${model}/${newfolder}/condor_${model}_${newfolder}.cfg
  echo "notify_user = mmelodea@cern.ch" >> ${model}/${newfolder}/condor_${model}_${newfolder}.cfg
  echo "Queue 1" >> ${model}/${newfolder}/condor_${model}_${newfolder}.cfg

  echo '#!/bin/bash' > ${model}/${newfolder}/${model}_${newfolder}.sh
  echo "#${m}" >> ${model}/${newfolder}/${model}_${newfolder}.sh
  echo "source /lustrehome/mmelodea/logincmspos.sh" >> ${model}/${newfolder}/${model}_${newfolder}.sh
  echo "cd /lustrehome/mmelodea/CMSSW_8_1_0/src" >> ${model}/${newfolder}/${model}_${newfolder}.sh
  echo 'eval `scram runtime -sh`' >> ${model}/${newfolder}/${model}_${newfolder}.sh
  echo "cd /lustrehome/mmelodea/Keras/Uncertainty/Sigma1/Datacards/CombineLimits/" >> ${model}/${newfolder}/${model}_${newfolder}.sh
  echo "cp ImpactPlot/text2workspace.py ImpactPlot/combineTool.py ImpactPlot/plotImpacts.py ${model}/${newfolder}/" >> ${model}/${newfolder}/${model}_${newfolder}.sh
  echo "cd ${model}/${newfolder}" >> ${model}/${newfolder}/${model}_${newfolder}.sh
  echo "text2workspace.py /lustrehome/mmelodea/Keras/Uncertainty/Sigma1/Datacards/combinedCard_Higgs125GeV_VBFStudy_${model}_4l.txt -m 125" >> ${model}/${newfolder}/${model}_${newfolder}.sh
  echo "python combineTool.py -M Impacts -d /lustrehome/mmelodea/Keras/Uncertainty/Sigma1/Datacards/combinedCard_Higgs125GeV_VBFStudy_${model}_4l.root -m 125 --doInitialFit --robustFit 1" >> ${model}/${newfolder}/${model}_${newfolder}.sh
  echo "python combineTool.py -M Impacts -d /lustrehome/mmelodea/Keras/Uncertainty/Sigma1/Datacards/combinedCard_Higgs125GeV_VBFStudy_${model}_4l.root -m 125 --robustFit 1 --doFits" >> ${model}/${newfolder}/${model}_${newfolder}.sh
  echo "python combineTool.py -M Impacts -d /lustrehome/mmelodea/Keras/Uncertainty/Sigma1/Datacards/combinedCard_Higgs125GeV_VBFStudy_${model}_4l.root -m 125 -o impacts_obs_4l_${model}.json" >> ${model}/${newfolder}/${model}_${newfolder}.sh
  echo "python plotImpacts.py -i impacts_obs_4l_${model}.json -o impacts_obs_4l_${model}" >> ${model}/${newfolder}/${model}_${newfolder}.sh
  condor_submit -name ettore ${model}/${newfolder}/condor_${model}_${newfolder}.cfg
done
