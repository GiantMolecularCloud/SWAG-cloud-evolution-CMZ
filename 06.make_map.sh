#!/bin/bash

# ####################################
# # ******************************** #
# # * SF sequence analysis script  * #
# # ******************************** #
# ####################################


# Script to make several maps from the fitted parameters

# parameters:
line=$1
pix_list_map='masks/NH3_1-1.3sigma_masked.mom0.5Jy.fits'
fit_dir='fits_NH3/'
fit_list=${fit_dir}${line}'.fit_list.txt'

###########################################################

# get empty map to match shape
cp -r ${pix_list_map} empty.fits
setpix 0 0 0 empty.fits

# start with an empty map for all maps
# and set unit correctly
mkdir maps/
rm -r maps/${line}.v.fits
rm -r maps/${line}.v_err.fits
rm -r maps/${line}.width.fits
rm -r maps/${line}.width_err.fits
rm -r maps/${line}.opacity.fits
rm -r maps/${line}.opacity.fits
rm -r maps/${line}.opacity_err.fits
rm -r maps/${line}.good_fits.fits
rm -r maps/${line}.col_dens.fits
rm -r maps/${line}.col_dens_err.fits

cp -r empty.fits maps/${line}.v.fits
cp -r empty.fits maps/${line}.v_err.fits
cp -r empty.fits maps/${line}.width.fits
cp -r empty.fits maps/${line}.width_err.fits
cp -r empty.fits maps/${line}.opacity.fits
cp -r empty.fits maps/${line}.opacity_err.fits
cp -r empty.fits maps/${line}.good_fits.fits
cp -r empty.fits maps/${line}.col_dens.fits
cp -r empty.fits maps/${line}.col_dens_err.fits

sethead -k BUNIT='km/s' maps/${line}.v.fits
sethead -k BUNIT='km/s' maps/${line}.v_err.fits
sethead -k BUNIT='km/s' maps/${line}.width.fits
sethead -k BUNIT='km/s' maps/${line}.width_err.fits
sethead -k BUNIT='' maps/${line}.opacity.fits
sethead -k BUNIT='' maps/${line}.opacity_err.fits
sethead -k BUNIT='FIT.GOODNESS' maps/${line}.good_fits.fits
sethead -k BUNIT='cm-2' maps/${line}.col_dens.fits
sethead -k BUNIT='cm-2' maps/${line}.col_dens_err.fits

# remove temp image
rm -r empty.fits

# set pixels according to fitted values
# create a map that shows good fits as +1, bad fits as -1 (cannot use 0 here!)
grep -v '^#' ${fit_list} | while read -r file_line
do
	# get values in this line
	good=`echo ${file_line} | awk '{print $3}'`

	if [ "${good}" == 1 ]
	then
		# read other values only if necessary
		x=`echo ${file_line} | awk '{print $1}'`
		y=`echo ${file_line} | awk '{print $2}'`
		v=`echo ${file_line} | awk '{print $6}'`
		v_err=`echo ${file_line} | awk '{print $7}'`
		width=`echo ${file_line} | awk '{print $8}'`
		width_err=`echo ${file_line} | awk '{print $9}'`
		tau=`echo ${file_line} | awk '{print $10}'`
		tau_err=`echo ${file_line} | awk '{print $11}'`
		coldens=`echo ${file_line} | awk '{print $19}'`
		coldens_err=`echo ${file_line} | awk '{print $20}'`

		# write to fits file		
		setpix ${x} ${y} ${v} maps/${line}.v.fits
		setpix ${x} ${y} ${v_err} maps/${line}.v_err.fits
		setpix ${x} ${y} ${width} maps/${line}.width.fits
		setpix ${x} ${y} ${width_err} maps/${line}.width_err.fits
		setpix ${x} ${y} ${tau} maps/${line}.opacity.fits
		setpix ${x} ${y} ${tau_err} maps/${line}.opacity_err.fits
		setpix ${x} ${y} ${good} maps/${line}.good_fits.fits
		setpix ${x} ${y} ${coldens} maps/${line}.col_dens.fits
		setpix ${x} ${y} ${coldens_err} maps/${line}.col_dens_err.fits
	
		# print status
		echo -e "processed pixel ${x},${y}"
	else
		setpix ${x} ${y} -1.0 maps/${line}.good_fits.fits
	fi
done
