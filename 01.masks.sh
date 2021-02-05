#!/bin/bash

##################
# MASKING SCRIPT #
##################

# Script to produce some masks of NH3 (1,1) data
# A mask containing only 3 sigma emission, maksed
# moment 0 map, mask of that map to use for later
# cloud fitting.

# parameters
noise=0.013	# Jy
edge_mask="SWAG_data/edge.100.cubemask"
mom_thres=5	# Jy/beam.km/s

##########################################################################


echo -e "\nmasking NH3 maps to contain only emission that can be fitted (no edge channels)\n"

# mask everything outside +-190km/s to get
# rid of edge channel artifacts that can screw up fitting
#
for j in `seq 1 6`
do
	immask in=SWAG_data/NH3_${j}-${j}.2kms.clean.gal.map \
		region="images(1,30)" \
		logic=and \
		flag=false
		
	immask in=SWAG_data/NH3_${j}-${j}.2kms.clean.gal.map \
		region="images(222,250)" \
		logic=and \
		flag=false

done




echo -e "\nmaking channel map mask\n"

rm -r masks/
mkdir masks/

# work with a new copy
cp -r SWAG_data/NH3_1-1.2kms.clean.gal.map masks/

# flag map edges
maths exp="<masks/NH3_1-1.2kms.clean.gal.map>*<${edge_mask}>" \
	out="masks/NH3_1-1.no_edge.map"

# mask at 3 sigma
thres=`echo "3 ${noise}" | awk '{printf "%.3f", $1 * $2}'`
maths exp="<masks/NH3_1-1.no_edge.map>*0.0+1.0" \
	mask="<masks/NH3_1-1.2kms.clean.gal.map>.gt.${thres}" \
	out="masks/NH3_1-1.3sigma.mask"

# delete header to interpret 0 as 0 Jy instead of masked
delhd in=masks/NH3_1-1.3sigma.mask/mask



echo -e "\ncalculating moments and moment mask\n"

# generate NH3 (1,1) moment 0 map
# this map contains only 3 sigma emission

# multiply map with mask because moment has no mask parameter
maths exp="<masks/NH3_1-1.3sigma.mask>*<masks/NH3_1-1.no_edge.map>" \
	out="masks/NH3_1-1.3sigma_masked.map"
	
# calculate moments
moment in="masks/NH3_1-1.3sigma_masked.map" \
	region="images(60,200)" \
	out="masks/NH3_1-1.3sigma_masked.mom0.map" \
	mom=0

# calculate moment mask
maths exp="<masks/NH3_1-1.3sigma_masked.mom0.map>*0.0+1.0" \
	mask="<masks/NH3_1-1.3sigma_masked.mom0.map>.gt.${mom_thres}" \
	out="masks/NH3_1-1.3sigma_masked.mom0.5Jy.mask"

# once again delete the mask header info
delhd in=masks/NH3_1-1.3sigma_masked.mom0.5Jy.mask/mask

# export to fits
fits in=masks/NH3_1-1.3sigma_masked.mom0.5Jy.mask \
	out=masks/NH3_1-1.3sigma_masked.mom0.5Jy.fits \
	op=xyout





