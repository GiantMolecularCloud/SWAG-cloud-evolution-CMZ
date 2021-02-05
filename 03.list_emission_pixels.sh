#!/bin/bash

# ####################################
# # ******************************** #
# # * SF sequence analysis script  * #
# # ******************************** #
# ####################################

# script to load miriad cubes into a 
# format that class can read

# Which file was used to create the masks?
# The same file should be used to get the list
# of meaningful pixels
file='masks/NH3_1-1.3sigma_masked.mom0.5Jy.fits'

# parameters:
cell=5		# image resolution in arcsec per pixel

############################################################

# get image size
#size_x=`imsize ${file} | awk '{print $(NF-1)}' | awk -F'x' '{print $1}'`
#size_y=`imsize ${file} | awk '{print $(NF-1)}' | awk -F'x' '{print $2}'`


# create list with emission pixels
echo -e "# pixels with significant emission to be fitted by CLASS" > pixel_list.txt
echo -e "# x\ty\txoff\tyoff" >> pixel_list.txt


# loop over x and y in image
# not necessary to loop every single pixel
# image is much smaller, so a smaller range is much faster
for x in `seq 150 1200`
do
	for y in `seq 550 1000`
	do
	
		# list pixel if it is included in the mask
		# i.e. it contains significant emission
		pixval=`getpix ${file} ${x} ${y} | awk '{printf "%.5f", $1}'`
		
		if [[ ${pixval} == 0.00000  ]]
		then
			echo -e "found excluded pixel at ${x},${y}"
		else
			echo -e "found value ${pixval} at ${x},${y}"
			
			# pixel to Gildas offset conversion
			# to be used with custom set header only!
			# offsets are defined physically, i.e. x is negative from left to right!
			xoff=`echo ${x} ${cell} | awk '{printf "%.1f", -1 * $1 * $2}'`
			yoff=`echo ${y} ${cell} | awk '{printf "%.1f", $1 * $2}'`
			
			# save values in list
			echo -e "${x}\t${y}\t${xoff}\t${yoff}" >> pixel_list.txt
		fi
	done
done


