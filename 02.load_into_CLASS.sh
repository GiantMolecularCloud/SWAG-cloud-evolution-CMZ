#!/bin/bash

# ####################################
# # ******************************** #
# # * SF sequence analysis script  * #
# # ******************************** #
# ####################################

# script to load miriad cubes into a
# format that class can read



#############################################################

rm -r CLASS_files
mkdir CLASS_files


# for all lines
for j in {1..6}
do

# get masked cube data
maths exp="<SWAG_data/NH3_${j}-${j}.2kms.clean.gal.map>*<masks/NH3_1-1.3sigma_masked.mom0.5Jy.mask>" \
	out="CLASS_files/NH3_${j}-${j}.masked.map" \
	options=grow

# get rid of masked channels or CLASS will cause problems
imsub in="CLASS_files/NH3_${j}-${j}.masked.map" \
	out="CLASS_files/NH3_${j}-${j}.masked.sub.map" \
	region="images(26,226)"

# export miriad images as fits
fits in=CLASS_files/NH3_${j}-${j}.masked.sub.map \
	out=CLASS_files/NH3_${j}-${j}.masked.fits \
	op=xyout

#############################################################

# modify FITS header to match pixel position with
# CLASS offset

# get galactic coordinates of pixel 0,0
# imhead lists 11 decimals
x_0=`xy2sky -g -n 11 CLASS_files/NH3_${j}-${j}.masked.fits 0 0 | awk '{printf $1}'`
y_0=`xy2sky -g -n 11 CLASS_files/NH3_${j}-${j}.masked.fits 0 0 | awk '{printf $2}'`

# set FITS header: reference pixel = 0, reference value = x0 or y0
cp -r CLASS_files/NH3_${j}-${j}.masked.fits CLASS_files/NH3_${j}-${j}.masked.regrid.fits
sethead -kl CRPIX1=0.00000000000E+00 CLASS_files/NH3_${j}-${j}.masked.regrid.fits
sethead -kl CRPIX2=0.00000000000E+00 CLASS_files/NH3_${j}-${j}.masked.regrid.fits
sethead -kl CRVAL1=${x_0} CLASS_files/NH3_${j}-${j}.masked.regrid.fits
sethead -kl CRVAL2=${y_0} CLASS_files/NH3_${j}-${j}.masked.regrid.fits



##############################################################

# import fits cubes into CLASS' gdf format
class << CLASSinput
file out CLASS_files/NH3_${j}-${j}.masked.regrid.atca single
lmv CLASS_files/NH3_${j}-${j}.masked.regrid.fits
CLASSinput

done
