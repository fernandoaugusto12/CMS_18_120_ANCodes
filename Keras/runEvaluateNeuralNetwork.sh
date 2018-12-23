#### script to submit multiple NN training jobs at lxplus
#### Author: Miqueias Melo de Almeida
#### Date: 05/06/2018

vusevars=( 
  "f_lept1_pt f_lept1_eta f_lept1_phi f_lept2_pt f_lept2_eta f_lept2_phi f_lept3_pt f_lept3_eta f_lept3_phi f_lept4_pt f_lept4_eta f_lept4_phi f_jets_highpt_pt[0] f_jets_highpt_eta[0] f_jets_highpt_phi[0] f_jets_highpt_pt[1] f_jets_highpt_eta[1] f_jets_highpt_phi[1]"
  "f_lept1_pt f_lept1_eta f_lept1_phi f_lept2_pt f_lept2_eta f_lept2_phi f_lept3_pt f_lept3_eta f_lept3_phi f_lept4_pt f_lept4_eta f_lept4_phi f_jets_highpt_pt[0] f_jets_highpt_eta[0] f_jets_highpt_phi[0] f_jets_highpt_pt[1] f_jets_highpt_eta[1] f_jets_highpt_phi[1] f_jets_highpt_pt[2] f_jets_highpt_eta[2] f_jets_highpt_phi[2]"
  "f_lept1_pt f_lept1_eta f_lept1_phi f_lept2_pt f_lept2_eta f_lept2_phi f_lept3_pt f_lept3_eta f_lept3_phi f_lept4_pt f_lept4_eta f_lept4_phi f_jets_highpt_pt[0] f_jets_highpt_eta[0] f_jets_highpt_phi[0] f_jets_highpt_pt[1] f_jets_highpt_eta[1] f_jets_highpt_phi[1] f_jets_highpt_pt[2] f_jets_highpt_eta[2] f_jets_highpt_phi[2] f_jets_highpt_pt[3] f_jets_highpt_eta[3] f_jets_highpt_phi[3]"
  "f_lept1_pt f_lept1_eta f_lept1_phi f_lept2_pt f_lept2_eta f_lept2_phi f_lept3_pt f_lept3_eta f_lept3_phi f_lept4_pt f_lept4_eta f_lept4_phi f_jets_highpt_pt[0] f_jets_highpt_eta[0] f_jets_highpt_phi[0] f_jets_highpt_e[0] f_jets_highpt_pt[1] f_jets_highpt_eta[1] f_jets_highpt_phi[1] f_jets_highpt_e[1]"
  "f_lept1_pt f_lept1_eta f_lept1_phi f_lept2_pt f_lept2_eta f_lept2_phi f_lept3_pt f_lept3_eta f_lept3_phi f_lept4_pt f_lept4_eta f_lept4_phi f_jets_highpt_pt[0] f_jets_highpt_eta[0] f_jets_highpt_phi[0] f_jets_highpt_pt[1] f_jets_highpt_eta[1] f_jets_highpt_phi[1] f_pfmet"
  "f_lept1_pt f_lept1_eta f_lept1_phi f_lept2_pt f_lept2_eta f_lept2_phi f_lept3_pt f_lept3_eta f_lept3_phi f_lept4_pt f_lept4_eta f_lept4_phi f_jets_highpt_pt[0] f_jets_highpt_eta[0] f_jets_highpt_phi[0] f_jets_highpt_pt[1] f_jets_highpt_eta[1] f_jets_highpt_phi[1] f_jets_highpt_pt[2] f_jets_highpt_eta[2] f_jets_highpt_phi[2] f_pfmet"
) 

vwait_for=("100")
vsbatch=("32" "128")
vlayers=("36 18 9" "21 13 8" "30" "30 5 10" "100")
vweights=("mc_weight" "event_weight")
voutlier=()


scheme=0
for vars in "${vusevars[@]}"
do
  for layers in "${vlayers[@]}"
  do
    for sbatch in "${vsbatch[@]}"
    do
      for patience in "${vwait_for[@]}"
      do
	for weight in "${vweights[@]}"
	do
          if [ -e Model${scheme}/core* ]
          then
            rm Model${scheme}/core*
	    echo "cd /afs/cern.ch/work/m/mmelodea/private/KerasJobs/Trainings" >> keras_run${scheme}.job
	    echo "mkdir -p Model${scheme}" >> keras_run${scheme}.job
	    echo "cp EvaluateNeuralNetwork.py Model${scheme}/" >> keras_run${scheme}.job
	    echo "cp UserFunctions.py Model${scheme}/" >> keras_run${scheme}.job
	    echo "cd Model${scheme}/" >> keras_run${scheme}.job
	    echo 'source /afs/cern.ch/user/m/mmelodea/loadKeras.sh' >> keras_run${scheme}.job
	    echo '' >> keras_run${scheme}.job
	    echo '#Run the training' >> keras_run${scheme}.job
	    echo 'python EvaluateNeuralNetwork.py --infile /afs/cern.ch/work/m/mmelodea/private/KerasJobs/Trainings/hzz4l_vbf_selection_m4l118-130GeV_shuffledFS.pkl \' >> keras_run${scheme}.job
	    echo '--keys ttX VBF VH HV ttH ZZ HJJ \' >> keras_run${scheme}.job
	    echo "--nninputs ${vars} --layers ${layers} --batchsize ${sbatch} --patience ${patience} --scaletrain ${weight} > training.log" >> keras_run${scheme}.job

	    echo "bsub -J keras${scheme} -q 1nh < keras_run${scheme}.job"
	    bsub -J keras${scheme} -q 1nh < keras_run${scheme}.job
          else
            echo "Model${scheme} is running or ready!"
          fi
	  scheme=$((scheme+1))
	done
      done
    done
  done
done
