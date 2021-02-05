#!/bin/bash

# ####################################
# # ******************************** #
# # * SF sequence analysis script  * #
# # ******************************** #
# ####################################


# Script to fit the spectra of the 6 ammonia
# (J,J) inversion lines.

# parameters:
line=$1
pixel_list='pixel_list.txt'

# do not confuse (spectral) line with file_line
# which is the line of a file to be read!


############################################################

mkdir -p fits_NH3/fit_logs/${line}
mkdir -p fits_NH3/fit_plots/${line}


# loop over all pixels with emission
while read -r file_line
do
	# ignore comments
	[[ "${file_line}" =~ ^#.*$ ]] && continue

	# get pixels and offset from list
	x=`echo ${file_line} | awk '{printf $1}'`
	y=`echo ${file_line} | awk '{printf $2}'`
	xoff=`echo ${file_line} | awk '{printf $3}'`
	yoff=`echo ${file_line} | awk '{printf $4}'`

	# set method to fit spectrum
	# use included function for low-j NH3, pre-defined HFS files
	# for higher j and a gaussian otherwise
	if [ "${line}" == "NH3_1-1" ]; then
		method="nh3(1,1)"
	elif [ "${line}" == "NH3_2-2" ]; then
		method="nh3(2,2)"
	elif [ "${line}" == "NH3_3-3" ]; then
		method="nh3(3,3)"
	elif [ "${line}" == "NH3_4-4" ]; then
		method="hfs scripts/nh3_44.hyp"
	elif [ "${line}" == "NH3_5-5" ]; then
		method="hfs scripts/nh3_55.hyp"
	elif [ "${line}" == "NH3_6-6" ]; then
		method="hfs scripts/nh3_66.hyp"
	fi

	# set initial guess for fit
	#initial_v=`egrep -v "^\s*(#|$)" ${cloud_list} | sed -n "${i}p" | awk '!/^[[:space:]]*#/{ print $5 }'`
	#initial_vmin=`egrep -v "^\s*(#|$)" ${cloud_list} | sed -n "${i}p" | awk '!/^[[:space:]]*#/{ print $4 }'`
	#initial_vmax=`egrep -v "^\s*(#|$)" ${cloud_list} | sed -n "${i}p" | awk '!/^[[:space:]]*#/{ print $6 }'`
	#initial_dv=`echo "${initial_vmax}-(${initial_vmin})" | bc`
	#initial_guess='"0 1   0 '${initial_v}'   0 '${initial_dv}'   0 0.2"'

	# print a status info
	echo " "
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	echo "fitting ${line} at pixel ${x},${y} with method ${method}"
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	echo " "
	echo "xoff: ${xoff}"
	echo "yoff: ${yoff}"
	echo " "
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"


# fit spectra in class
# (here documents cannot have indentions)
class << CLASSinput
! load spectra file
file in CLASS_files/${line}.masked.regrid.atca
set unit v
set match 2
!
! set offset and show spectrum
! unit: arcsec
set off ${xoff} ${yoff}
find
get first
plot
!
! set fitting method
file out temp.${line}.fit single \over
method ${method}
!
! set initial guess
!lines 1 ${initial}
!
minimize
write
vis
file in temp.${line}.fit
find
print fit /output fits_NH3/fit_logs/${line}/${line}.pix_${x}_${y}.fit.txt
hardcopy fits_NH3/fit_plots/${line}/${line}.pix_${x}_${y}.png /device png /overwrite
CLASSinput
	
		
	# remove temp
	rm temp.${line}.fit

done < ${pixel_list}
