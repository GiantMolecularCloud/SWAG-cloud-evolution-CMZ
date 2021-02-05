#!/bin/bash

# ####################################
# # ******************************** #
# # * SF sequence analysis script  * #
# # ******************************** #
# ####################################


# Script to calculate temperature maps

# parameters:
NH3_data=(24.4 65.6 124.7 201.7 296.5 409.2)

#########################################################

# setup
rm -r maps/*.col_dens.map
rm -r maps/*.col_dens_err.map
rm -r maps/T??_rot.map
rm -r maps/T??_kin.map
rm -r maps/T??_rot_err.map
rm -r maps/T??_kin_err.map
rm -r maps/T??_rot.fits
rm -r maps/T??_kin.fits
rm -r maps/T??_rot_err.fits
rm -r maps/T??_kin_err.fits
rm -r maps/yaxis_maps
rm -r maps/slope_maps
mkdir -p maps/yaxis_maps
mkdir -p maps/slope_maps

#########################################################
# yaxis

for j in `seq 1 6`
do
	line="NH3_${j}-${j}"
	echo -e "\ncalculating yaxis for Boltzmann plot; file: ${line}\n"
	
	# load fits images into miriad map format
	fits in=maps/${line}.col_dens.fits \
		out=maps/${line}.col_dens.map \
		op=xyin
	fits in=maps/${line}.col_dens_err.fits \
		out=maps/${line}.col_dens_err.map \
		op=xyin
		
	# calculate yaxis images: y = N/(2J+1)/g
	# statistical weight is g=1 for para-NH3 and g=2 for ortho-NH3
	if [ "$j" == 3 ] || [ "$j" == 6 ]
	then
		g=2
	else
		g=1
	fi
	
	maths exp="log10(<maps/${line}.col_dens.map>/(2*${j}+1)/${g})" \
		mask="<maps/${line}.col_dens.map>.gt.0.0" \
		out=maps/yaxis_maps/${line}.yaxis.map
	maths exp="<maps/${line}.col_dens_err.map>/(<maps/${line}.col_dens.map>*log(10))" \
		mask="<maps/${line}.col_dens_err.map>.gt.0.0" \
		out=maps/yaxis_maps/${line}.yaxis_err.map
done


#########################################################
# rotational temperatures

for t in "12" "24" "45" "36"
do
	echo -e "\ncalculating temperature map: T${t}\n"
	t1=`echo ${t:0:1}`
	t2=`echo ${t:1:1}`

	# calculate temperature
	delta_T=`echo ${NH3_data[t1-1]} ${NH3_data[t2-1]} | awk '{print $2-$1}'`
	maths exp="(<maps/yaxis_maps/NH3_${t2}-${t2}.yaxis.map>-<maps/yaxis_maps/NH3_${t1}-${t1}.yaxis.map>)/${delta_T}" \
		mask="(<maps/yaxis_maps/NH3_${t2}-${t2}.yaxis.map>.gt.0.0).and.(<maps/yaxis_maps/NH3_${t2}-${t2}.yaxis.map>.gt.0.0)" \
		out="maps/slope_maps/slope${t1}${t2}.map"
	maths exp="sqrt(<maps/yaxis_maps/NH3_${t2}-${t2}.yaxis_err.map>**2+<maps/yaxis_maps/NH3_${t1}-${t1}.yaxis_err.map>**2)/${delta_T}" \
		mask="(<maps/yaxis_maps/NH3_${t2}-${t2}.yaxis.map>.gt.0.0).and.(<maps/yaxis_maps/NH3_${t2}-${t2}.yaxis.map>.gt.0.0)" \
		out="maps/slope_maps/slope${t1}${t2}_err.map"
	maths exp="-1.0*log10(exp(1))/<maps/slope_maps/slope${t1}${t2}.map>" \
		mask="<maps/slope_maps/slope${t1}${t2}.map>.ne.0.0" \
		out="maps/T${t1}${t2}_rot.map"
	maths exp="log10(exp(1))/(<maps/slope_maps/slope${t1}${t2}.map>**2)*<maps/slope_maps/slope${t1}${t2}_err.map>" \
		mask="<maps/slope_maps/slope${t1}${t2}.map>.ne.0.0" \
		out="maps/T${t1}${t2}_rot_err.map"
	
	# set unit
	puthd in=maps/T${t1}${t2}_rot.map/bunit \
		value='K'
	puthd in=maps/T${t1}${t2}_rot_err.map/bunit \
		value='K'

	# export to fits
	fits in=maps/T${t1}${t2}_rot.map \
		out=maps/T${t1}${t2}_rot.fits \
		op=xyout
	fits in=maps/T${t1}${t2}_rot_err.map \
		out=maps/T${t1}${t2}_rot_err.fits \
		op=xyout
done


#########################################################
# kinetic temperatures
# asymetric errors are estimated by the larger error!

echo -e "\ncalculating kinetic temperature maps\n"

# T12
maths exp="6.05*exp(0.06088*<maps/T12_rot.map>)" \
	mask="<maps/T12_rot.map>.lt.100" \
	out="maps/T12_kin.map"
maths exp="6.05*0.06088*exp(0.06088*<maps/T12_rot.map>)*<maps/T12_rot_err.map>" \
	mask="<maps/T12_rot.map>.lt.100" \
	out="maps/T12_kin_err.map"


# T24
# for 0 < T_kin < 100K
maths exp="1.467*<maps/T24_rot.map>-6.984" \
	mask="<maps/T24_rot.map>.lt.200" \
	out="maps/T24_kin.map"
maths exp="1.467*<maps/T24_rot_err.map>" \
	mask="<maps/T24_rot.map>.lt.200" \
	out="maps/T24_kin_err.map"
# for 100K < T_kin < 500K
#T_kin = 27.085*np.exp(T_rot*0.019)
#delta_T_kin = T_kin*0.019*delta_T_rot


# T45
# for 50K < T_kin < 100K
maths exp="21.024*exp(<maps/T24_rot.map>*0.0198)" \
	mask="<maps/T45_rot.map>.lt.300" \
	out="maps/T45_kin.map"
maths exp="21.024*0.0198*exp(0.0198*<maps/T45_rot.map>)*<maps/T45_rot_err.map>" \
	mask="<maps/T45_rot.map>.lt.300" \
	out="maps/T45_kin_err.map"
# for 0 < T_kin < 50K
#T_kin = 1.143*T_rot-1.611
#delta_T_kin = 1.143*delta_T_rot


# T36
# for T_rot > 50K
maths exp="28.87511*exp(0.015*<maps/T24_rot.map>)" \
	mask="<maps/T36_rot.map>.lt.400" \
	out="maps/T36_kin.map"
maths exp="28.87511*0.015*exp(0.015*<maps/T36_rot.map>)*<maps/T36_rot_err.map>" \
	mask="<maps/T36_rot.map>.lt.400" \
	out="maps/T36_kin_err.map"
# for T_rot < 50K
#T_kin = T_rot
#delta_T_kin = delta_T_rot


for t in "12" "24" "45" "36"
do
	fits in="maps/T${t}_kin.map" \
		out="maps/T${t}_kin.fits" \
		op=xyout
	fits in="maps/T${t}_kin_err.map" \
		out="maps/T${t}_kin_err.fits" \
		op=xyout
done
