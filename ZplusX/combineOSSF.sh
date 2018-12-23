ths=( "k41nj2" "k1177nj2" "k2nj3" "k4nj3" "k16nj2e3" "k53nj2e3" "k196nj2e3" "k543nj2e3" )
fmp=( "nominal" "average" "reweight" )
srs=( "smhiggs" "vbf" )
fss=( "4mu" "4e" "2e2mu" "2mu2e" )


for isr in "${srs[@]}"
do
  for ith in "${ths[@]}"
  do
    for imp in "${fmp[@]}"
    do
      hadd -f ${isr}/DistributionsUncertainties_${ith}_proc_zx_4l_${isr}_${imp}.root    \
	      ${isr}/DistributionsUncertainties_${ith}_proc_zx_4mu_${isr}_${imp}.root   \
	      ${isr}/DistributionsUncertainties_${ith}_proc_zx_4e_${isr}_${imp}.root    \
	      ${isr}/DistributionsUncertainties_${ith}_proc_zx_2e2mu_${isr}_${imp}.root \
	      ${isr}/DistributionsUncertainties_${ith}_proc_zx_2mu2e_${isr}_${imp}.root
    done
    
    hadd -f ${isr}/DistributionsUncertainties_${ith}_proc_zx_4l_${isr}.root         \
	    ${isr}/DistributionsUncertainties_${ith}_proc_zx_4l_${isr}_nominal.root \
	    ${isr}/DistributionsUncertainties_${ith}_proc_zx_4l_${isr}_average.root \
	    ${isr}/DistributionsUncertainties_${ith}_proc_zx_4l_${isr}_reweight.root
  done  
done