for i in MELA NN1 NN2 NN3 NN6 NN10 NN11 NN372_nj2_bs1 NN412_nj2_bs1 NN1098_nj2 NN1137_nj2 NN1177_nj2 NN11_nj2 NN131_nj2 NN138_nj2 NN211_nj2 NN51_nj2 NN91_nj2 #NN1244_nj3
do
  echo "Combining cards for ${i}"
  combineCards.py ch4mu=datacard_Higgs125GeV_VBFStudy_${i}_ch4mu.txt \
                  ch4e=datacard_Higgs125GeV_VBFStudy_${i}_ch4e.txt \
                  ch2e2mu=datacard_Higgs125GeV_VBFStudy_${i}_ch2e2mu.txt \
                  > combinedCard_Higgs125GeV_VBFStudy_${i}_4l.txt
done
