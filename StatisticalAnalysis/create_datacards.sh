shapes=("--shapes l1pt l2pt l3pt l4pt j1pt j2pt"
        "--shapes l1pt l2pt l3pt l4pt j1pt j2pt"
        "--shapes l1pt l2pt l3pt l4pt j1pt j2pt j3pt"
        "--shapes l1pt l2pt l3pt l4pt j1pt j2pt j3pt j4pt"
        "--shapes l1pt l2pt l3pt l4pt j1pt j2pt"
        "--shapes l1pt l2pt l3pt l4pt j1pt j2pt met"
        "--shapes l1pt l2pt l3pt l4pt j1pt j2pt j3pt met"
        "--shapes l1pt l2pt l3pt l4pt j1pt j2pt met"
        "--shapes l1pt l2pt l3pt l4pt j1pt j2pt met"
        "--shapes l1pt l2pt l3pt l4pt j1pt j2pt met"
        "--shapes l1pt l2pt l3pt l4pt j1pt j2pt met"
        "--shapes l1pt l2pt l3pt l4pt j1pt j2pt met"
        "--shapes l1pt l2pt l3pt l4pt j1pt j2pt met"
        "--shapes l1pt l2pt l3pt l4pt j1pt j2pt met"
        "--shapes l1pt l2pt l3pt l4pt j1pt j2pt met"
        "--shapes l1pt l2pt l3pt l4pt j1pt j2pt met"
        "--shapes l1pt l2pt l3pt l4pt j1pt j2pt met"
        "--shapes l1pt l2pt l3pt l4pt j1pt j2pt met"
        "--shapes l1pt l2pt l3pt l4pt j1pt j2pt j3pt"
)

ic=0
for c in MELA NN1 NN2 NN3 NN6 NN10 NN11 NN372_nj2_bs1 NN412_nj2_bs1 NN1098_nj2 NN1137_nj2 NN1177_nj2 NN11_nj2 NN131_nj2 NN138_nj2 NN211_nj2 NN51_nj2 NN91_nj2 NN1244_nj3
do
  python genCombineInputs.py --infile DistributionsUncertainties_${c}.root \
                               --outfile Higgs125GeV_VBFStudy_${c} \
                               --channels ch4mu ch4e ch2e2mu \
                               --keys data_obs qqH ggH WH ZH ttH ggZZ qqZZ HWW TTV VVV \
                               --signal qqH \
                               ${shapes[$ic]} \
                               --rebin -1
  ic=$((ic+1))
done
