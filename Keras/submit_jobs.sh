#### script to submit multiple NN training jobs at Bari RECAS
#### Author: Miqueias Melo de Almeida
#### Date: 06/06/2018

vnninputs2=( 
  "f_lept1_pt f_lept1_eta f_lept1_phi f_lept2_pt f_lept2_eta f_lept2_phi f_lept3_pt f_lept3_eta f_lept3_phi f_lept4_pt f_lept4_eta f_lept4_phi f_jets_highpt_pt[0] f_jets_highpt_eta[0] f_jets_highpt_phi[0] f_jets_highpt_pt[1] f_jets_highpt_eta[1] f_jets_highpt_phi[1]"
  #"f_lept1_pt f_lept1_eta f_lept1_phi f_lept2_pt f_lept2_eta f_lept2_phi f_lept3_pt f_lept3_eta f_lept3_phi f_lept4_pt f_lept4_eta f_lept4_phi f_jets_highpt_pt[0] f_jets_highpt_eta[0] f_jets_highpt_phi[0] f_jets_highpt_e[0] f_jets_highpt_pt[1] f_jets_highpt_eta[1] f_jets_highpt_phi[1] f_jets_highpt_e[1]"
  "f_lept1_pt f_lept1_eta f_lept1_phi f_lept2_pt f_lept2_eta f_lept2_phi f_lept3_pt f_lept3_eta f_lept3_phi f_lept4_pt f_lept4_eta f_lept4_phi f_jets_highpt_pt[0] f_jets_highpt_eta[0] f_jets_highpt_phi[0] f_jets_highpt_pt[1] f_jets_highpt_eta[1] f_jets_highpt_phi[1] f_pfmet"
  #"f_lept1_pt f_lept1_eta f_lept1_phi f_lept2_pt f_lept2_eta f_lept2_phi f_lept3_pt f_lept3_eta f_lept3_phi f_lept4_pt f_lept4_eta f_lept4_phi f_jets_highpt_pt[0] f_jets_highpt_eta[0] f_jets_highpt_phi[0] f_jets_highpt_e[0] f_jets_highpt_pt[1] f_jets_highpt_eta[1] f_jets_highpt_phi[1] f_jets_highpt_e[1] f_pfmet"
) 

vnninputs3=(
  "f_lept1_pt f_lept1_eta f_lept1_phi f_lept2_pt f_lept2_eta f_lept2_phi f_lept3_pt f_lept3_eta f_lept3_phi f_lept4_pt f_lept4_eta f_lept4_phi f_jets_highpt_pt[0] f_jets_highpt_eta[0] f_jets_highpt_phi[0] f_jets_highpt_pt[1] f_jets_highpt_eta[1] f_jets_highpt_phi[1] f_jets_highpt_pt[2] f_jets_highpt_eta[2] f_jets_highpt_phi[2]"
  #"f_lept1_pt f_lept1_eta f_lept1_phi f_lept2_pt f_lept2_eta f_lept2_phi f_lept3_pt f_lept3_eta f_lept3_phi f_lept4_pt f_lept4_eta f_lept4_phi f_jets_highpt_pt[0] f_jets_highpt_eta[0] f_jets_highpt_phi[0] f_jets_highpt_e[0] f_jets_highpt_pt[1] f_jets_highpt_eta[1] f_jets_highpt_phi[1] f_jets_highpt_e[1] f_jets_highpt_pt[2] f_jets_highpt_eta[2] f_jets_highpt_phi[2] f_jets_highpt_e[2]"
  "f_lept1_pt f_lept1_eta f_lept1_phi f_lept2_pt f_lept2_eta f_lept2_phi f_lept3_pt f_lept3_eta f_lept3_phi f_lept4_pt f_lept4_eta f_lept4_phi f_jets_highpt_pt[0] f_jets_highpt_eta[0] f_jets_highpt_phi[0] f_jets_highpt_pt[1] f_jets_highpt_eta[1] f_jets_highpt_phi[1] f_jets_highpt_pt[2] f_jets_highpt_eta[2] f_jets_highpt_phi[2] f_pfmet"
  #"f_lept1_pt f_lept1_eta f_lept1_phi f_lept2_pt f_lept2_eta f_lept2_phi f_lept3_pt f_lept3_eta f_lept3_phi f_lept4_pt f_lept4_eta f_lept4_phi f_jets_highpt_pt[0] f_jets_highpt_eta[0] f_jets_highpt_phi[0] f_jets_highpt_e[0] f_jets_highpt_pt[1] f_jets_highpt_eta[1] f_jets_highpt_phi[1] f_jets_highpt_e[1] f_jets_highpt_pt[2] f_jets_highpt_eta[2] f_jets_highpt_phi[2] f_jets_highpt_e[2] f_pfmet"
)

vminimizer=(
#  "sgd" 
#  "adam" 
  "adagrad" 
  "adadelta"
#  "rmsprop"
)

vneuron=(
#  "relu"
  "selu"
#  "tanh"
)

vpatiences=(
  "50" 
#  "600"
#  "30000"
)

vbatches=(
#  "5"
#  "32"
#  "62"
  "128"
  "786"
)

vtopologies=(
#  "9 7 5"
#  "11 9 7"
#  "15 10 5"
#  "21 13 8"
#  "10 10 10 10" 
#  "30"
#  "100"
#  "100 100 50"
  "10000"
)

dropout=(
  "0.1"
#  "0.3"
#  "0.5"
#  "0.7"
#  "0.9"
#  "0.95"
#  "0.99"
#  "0.3 0.4 0.2"
#  "0.5 0.25 0.1"
)

vweights=(
  "mc_weight" 
  #"event_weight"
)

voutliers=(
  "" 
  #"--nooutliers"
)


for njet in 2
do
  vnninputs=()
  if [ $njet -eq 2 ]
  then
    vnninputs=("${vnninputs2[@]}")
  else
    vnninputs=("${vnninputs3[@]}")
  fi

  scheme=0
  #for outliers in "${voutliers[@]}"
  #do
    for topology in "${vtopologies[@]}"
    do
      for drop in "${dropout[@]}"
      do
        for batch in "${vbatches[@]}"
        do
          for patience in "${vpatiences[@]}"
          do
	    for weight in "${vweights[@]}"
	    do
              for nninputs in "${vnninputs[@]}"
              do
                #echo "Inputs: ${nninputs}"
                for neuron in "${vneuron[@]}"
                do
                  for minimizer in "${vminimizer[@]}"
                  do
                    path="/lustre/cms/store/user/mmelodea/Keras/Trainings/"
                    if [ -e ${path}/Njets${njet}/Model${scheme}/SummaryOfResults.pkl ]
                    then
                      echo "Njets${njet}/Model${scheme}/results/SummaryOfResults.pkl is ready!"
                    else
                    #elif [ ! $(-e ${path}/Njets${njet}/Model${scheme}/results) ]
                    #if [ ${scheme} < 10 ]
                    #then
                      #creates a folder to hold and make training
                      mkdir -p Njets${njet}/Model${scheme}

                      #create sh script for training
                      echo '#!/bin/bash' > Njets${njet}/Model${scheme}/k${scheme}nj${njet}.sh
                      echo "cd /lustre/cms/store/user/mmelodea/Keras/Trainings/Njets${njet}/Model${scheme}" >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.sh
                      echo "cp /lustre/cms/store/user/mmelodea/Keras/Trainings/EvaluateNeuralNetwork.py EvaluateNeuralNetwork_${scheme}.py" >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.sh
                      #echo "cp /lustre/cms/store/user/mmelodea/Keras/Trainings/Samples/hzz4l_vbf_selection_m4l118-130GeV_shuffledFS_njets${njet}.pkl ." >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.sh
                      #echo "source /lustrehome/mmelodea/loadKeras.sh" >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.sh
                      echo "source /cvmfs/sft.cern.ch/lcg/views/LCG_94/x86_64-slc6-gcc62-opt/setup.sh" >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.sh
                      echo "python /lustre/cms/store/user/mmelodea/Keras/Trainings/Njets${njet}/Model${scheme}/EvaluateNeuralNetwork_${scheme}.py --infile /lustre/cms/store/user/mmelodea/Keras/Trainings/Samples/hzz4l_vbf_selection_m4l118-130GeV_shuffledFS_njets${njet}.pkl --keys qqH ggH ggZZ ZZJJ qqZZ ttH WH ZH ZX --signal qqH --nninputs ${nninputs} --topology ${topology} --batchsize ${batch} --patience ${patience} --scaletrain ${weight} --split 0.8 --nepochs 400 --neuron ${neuron} --minimizer ${minimizer} >& training.log" >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.sh
                      #echo "mv /lustre/cms/store/user/mmelodea/Keras/Trainings/Njets${njet}/Model${scheme}/training.log /lustre/cms/store/user/mmelodea/Keras/Trainings/Njets${njet}/Model${scheme}/results/" >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.sh
                      #echo "rm /lustre/cms/store/user/mmelodea/Keras/Trainings/Njets${njet}/Model${scheme}/hzz4l_vbf_selection_m4l118-130GeV_shuffledFS_njets${njet}.pkl" >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.sh
                      chmod +x /lustre/cms/store/user/mmelodea/Keras/Trainings/Njets${njet}/Model${scheme}/k${scheme}nj${njet}.sh
                      #chmod +x /lustre/cms/store/user/mmelodea/Keras/Trainings/Njets${njet}/Model${scheme}/EvaluateNeuralNetwork_${scheme}.py

                      #creates the cfg file to call the sh
                      echo "universe = vanilla" > Njets${njet}/Model${scheme}/k${scheme}nj${njet}.cfg
                      echo "Executable = /lustre/cms/store/user/mmelodea/Keras/Trainings/Njets${njet}/Model${scheme}/k${scheme}nj${njet}.sh" >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.cfg
                      #echo "Should_Transfer_Files = YES" >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.cfg
                      #echo "getenv = True" >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.cfg
                      #echo 'environment  = "PYTHONPATH=/lustrehome/mmelodea/virtualenv-1.9/keras/lib/python2.7/site-packages"' >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.cfg

                      echo "WhenToTransferOutput = ON_EXIT" >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.cfg
                      echo 'Requirements = TARGET.OpSys == "LINUX" && (TARGET.Arch != "DUMMY" )' >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.cfg
                      #echo 'Requirements = (TARGET.Machine == "slot1@wn-infn*.recas.ba.infn.it")' >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.cfg
                      #echo "notify_user = mmelodea@cern.ch" >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.cfg
                      #echo "Log = /lustre/cms/store/user/mmelodea/Keras/Trainings/Njets${njet}/Model${scheme}/k${scheme}nj${njet}.log" >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.cfg
                      #echo "Output = /lustre/cms/store/user/mmelodea/Keras/Trainings/Njets${njet}/Model${scheme}/k${scheme}nj${njet}.out" >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.cfg
                      #echo "Error = /lustre/cms/store/user/mmelodea/Keras/Trainings/Njets${njet}/Model${scheme}/k${scheme}nj${njet}.error" >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.cfg
                      echo "Queue " >> Njets${njet}/Model${scheme}/k${scheme}nj${njet}.cfg

                      ###then submit the job
                      echo "condor_submit -name ettore k${scheme}nj${njet}.cfg"
                      cd /lustre/cms/store/user/mmelodea/Keras/Trainings/Njets${njet}/Model${scheme}/
                      condor_submit -name ettore k${scheme}nj${njet}.cfg
                      cd /lustre/cms/store/user/mmelodea/Keras/Trainings/

                    #else
                    #  echo "Njets${njet}/Model${scheme} is running/pending!"
                    fi
                    scheme=$((scheme+1))
                  done
                done
              done
            done
          done
	done
      done
    done
  done
#done
