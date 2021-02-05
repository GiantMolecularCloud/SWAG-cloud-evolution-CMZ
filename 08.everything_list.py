import numpy as np
from astropy.io import fits


# ####################################
# # ******************************** #
# # * SF sequence analysis script  * #
# # ******************************** #
# ####################################


# Make a temperature list

# parameters:


############################################################

# load data
###########

print "\nloading data ...\n"

pixels = np.genfromtxt('pixel_list.txt',
		dtype = "i4, i4, f8, f8",
		names = ['x', 'y', 'xoff', 'yoff']
		)
stream = np.genfromtxt('Kruijssen_orbit.dat', 
		comments='#', 
		dtype='f8',
		names=('t', 'x', 'y', 'z', 'R', 'v_x','v_y','v_z','v_orb','l','b','v_los','mu_l','mu_b','mu_x','mu_y')
		)
fit1 = np.genfromtxt('fits_NH3/NH3_1-1.fit_list.txt',
	usecols = (0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21),
	dtype   = "i4, i4, i4, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, |S10, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8",
	names   = ["x", "y", "good", "area", "area_err", "v", "v_err", "width", "width_err", "tau", "tau_err", "baserms", " linerms", "|", "T_L", "T_L_err", "Tdv", "Tdv_err", "N", "N_err", "yaxis", "yaxis_err"]
	)
fit2 = np.genfromtxt('fits_NH3/NH3_2-2.fit_list.txt',
	usecols = (0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21),
	dtype   = "i4, i4, i4, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, |S10, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8",
	names   = ["x", "y", "good", "area", "area_err", "v", "v_err", "width", "width_err", "tau", "tau_err", "baserms", " linerms", "|", "T_L", "T_L_err", "Tdv", "Tdv_err", "N", "N_err", "yaxis", "yaxis_err"]
	)
fit3 = np.genfromtxt('fits_NH3/NH3_3-3.fit_list.txt',
	usecols = (0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21),
	dtype   = "i4, i4, i4, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, |S10, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8",
	names   = ["x", "y", "good", "area", "area_err", "v", "v_err", "width", "width_err", "tau", "tau_err", "baserms", " linerms", "|", "T_L", "T_L_err", "Tdv", "Tdv_err", "N", "N_err", "yaxis", "yaxis_err"]
	)
fit4 = np.genfromtxt('fits_NH3/NH3_4-4.fit_list.txt',
	usecols = (0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21),
	dtype   = "i4, i4, i4, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, |S10, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8",
	names   = ["x", "y", "good", "area", "area_err", "v", "v_err", "width", "width_err", "tau", "tau_err", "baserms", " linerms", "|", "T_L", "T_L_err", "Tdv", "Tdv_err", "N", "N_err", "yaxis", "yaxis_err"]
	)
fit5 = np.genfromtxt('fits_NH3/NH3_5-5.fit_list.txt',
	usecols = (0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21),
	dtype   = "i4, i4, i4, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, |S10, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8",
	names   = ["x", "y", "good", "area", "area_err", "v", "v_err", "width", "width_err", "tau", "tau_err", "baserms", " linerms", "|", "T_L", "T_L_err", "Tdv", "Tdv_err", "N", "N_err", "yaxis", "yaxis_err"]
	)
fit6 = np.genfromtxt('fits_NH3/NH3_6-6.fit_list.txt',
	usecols = (0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21),
	dtype   = "i4, i4, i4, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, |S10, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8",
	names   = ["x", "y", "good", "area", "area_err", "v", "v_err", "width", "width_err", "tau", "tau_err", "baserms", " linerms", "|", "T_L", "T_L_err", "Tdv", "Tdv_err", "N", "N_err", "yaxis", "yaxis_err"]
	)
T12_kin     = fits.open('maps/T12_kin.fits')[0].data[0]
T12_kin_err = fits.open('maps/T12_kin_err.fits')[0].data[0]
T24_kin     = fits.open('maps/T24_kin.fits')[0].data[0]
T24_kin_err = fits.open('maps/T24_kin_err.fits')[0].data[0]
T45_kin     = fits.open('maps/T45_kin.fits')[0].data[0]
T45_kin_err = fits.open('maps/T45_kin_err.fits')[0].data[0]
T36_kin     = fits.open('maps/T36_kin.fits')[0].data[0]
T36_kin_err = fits.open('maps/T36_kin_err.fits')[0].data[0]

T12_rot     = fits.open('maps/T12_rot.fits')[0].data[0]
T12_rot_err = fits.open('maps/T12_rot_err.fits')[0].data[0]
T24_rot     = fits.open('maps/T24_rot.fits')[0].data[0]
T24_rot_err = fits.open('maps/T24_rot_err.fits')[0].data[0]
T45_rot     = fits.open('maps/T45_rot.fits')[0].data[0]
T45_rot_err = fits.open('maps/T45_rot_err.fits')[0].data[0]
T36_rot     = fits.open('maps/T36_rot.fits')[0].data[0]
T36_rot_err = fits.open('maps/T36_rot_err.fits')[0].data[0]

dust = fits.open('HiGAL.dust_temp.fits')[0].data


############################################################

# empty fit/temperature/time lists
mega_list = []

# header coordinate information
SWAG_header  = fits.open('maps/T12_kin.fits')[0].header
SWAG_CRVAL1 = SWAG_header['CRVAL1']
SWAG_CRPIX1 = SWAG_header['CRPIX1']
SWAG_CDELT1 = SWAG_header['CDELT1']
SWAG_CRVAL2 = SWAG_header['CRVAL2']
SWAG_CRPIX2 = SWAG_header['CRPIX2']
SWAG_CDELT2 = SWAG_header['CDELT2']
HiGAL_header  = fits.open('HiGAL.dust_temp.fits')[0].header
HiGAL_CRVAL1 = HiGAL_header['CRVAL1']
HiGAL_CRPIX1 = HiGAL_header['CRPIX1']
HiGAL_CDELT1 = HiGAL_header['CDELT1']
HiGAL_CRVAL2 = HiGAL_header['CRVAL2']
HiGAL_CRPIX2 = HiGAL_header['CRPIX2']
HiGAL_CDELT2 = HiGAL_header['CDELT2']


# loop *only* over pixels that *might* contain a value
# most of them do because fit success rate is high
# ! python and kvis sort axis differently !
# ! pixel list starts at 1, but python counts from 0 !
print "\nmaking huge table ...\n"
for line in np.arange(len(pixels)):

	# get pixel position
	x = pixels['x'][line]
	y = pixels['y'][line]
	
	# calculate galactic coordinates at this position
	glon = SWAG_CRVAL1+(x-SWAG_CRPIX1)*SWAG_CDELT1
	glat = SWAG_CRVAL2+(y-SWAG_CRPIX2)*SWAG_CDELT2
	
	
	# get row in the various files for this pixel
	#############################################
	
	row_fit1 = np.where(np.logical_and(fit1['x'] == x, fit1['y'] == y))[0]
	row_fit2 = np.where(np.logical_and(fit2['x'] == x, fit2['y'] == y))[0]
	row_fit3 = np.where(np.logical_and(fit3['x'] == x, fit3['y'] == y))[0]
	row_fit4 = np.where(np.logical_and(fit4['x'] == x, fit4['y'] == y))[0]
	row_fit5 = np.where(np.logical_and(fit5['x'] == x, fit5['y'] == y))[0]
	row_fit6 = np.where(np.logical_and(fit6['x'] == x, fit6['y'] == y))[0]
	
	
	
	# gather values if they exist for this pixel and store in temporary list
	# set value to nan if it does not exist or the fit was bad
	########################################################################
	
	v_lsr_list = []
	if (row_fit1.size == 1):
		if (fit1['good'][row_fit1] == 1):
			temp_fit1 = [fit1['T_L'][row_fit1][0], fit1['T_L_err'][row_fit1][0], fit1['width'][row_fit1][0], fit1['width_err'][row_fit1][0], fit1['v'][row_fit1][0], fit1['v_err'][row_fit1][0], fit1['tau'][row_fit1][0], fit1['tau_err'][row_fit1][0], fit1['N'][row_fit1][0], fit1['N_err'][row_fit1][0]]
			v_lsr_list.append(fit1['v'][row_fit1])
		else:
			temp_fit1 = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
	else:
		temp_fit1 = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

	if (row_fit2.size == 1):
		if (fit2['good'][row_fit2] == 1):
			temp_fit2 = [fit2['T_L'][row_fit2][0], fit2['T_L_err'][row_fit2][0], fit2['width'][row_fit2][0], fit2['width_err'][row_fit2][0], fit2['v'][row_fit2][0], fit2['v_err'][row_fit2][0], fit2['tau'][row_fit2][0], fit2['tau_err'][row_fit2][0], fit2['N'][row_fit2][0], fit2['N_err'][row_fit2][0]]
			v_lsr_list.append(fit2['v'][row_fit2])
		else:
			temp_fit2 = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
	else:
		temp_fit2 = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

	if (row_fit3.size == 1):
		if (fit3['good'][row_fit3] == 1):
			temp_fit3 = [fit3['T_L'][row_fit3][0], fit3['T_L_err'][row_fit3][0], fit3['width'][row_fit3][0], fit3['width_err'][row_fit3][0], fit3['v'][row_fit3][0], fit3['v_err'][row_fit3][0], fit3['tau'][row_fit3][0], fit3['tau_err'][row_fit3][0], fit3['N'][row_fit3][0], fit3['N_err'][row_fit3][0]]
			v_lsr_list.append(fit3['v'][row_fit3])
		else:
			temp_fit3 = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
	else:
		temp_fit3 = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

	if (row_fit4.size == 1):
		if (fit4['good'][row_fit4] == 1):
			temp_fit4 = [fit4['T_L'][row_fit4][0], fit4['T_L_err'][row_fit4][0], fit4['width'][row_fit4][0], fit4['width_err'][row_fit4][0], fit4['v'][row_fit4][0], fit4['v_err'][row_fit4][0], fit4['tau'][row_fit4][0], fit4['tau_err'][row_fit4][0], fit4['N'][row_fit4][0], fit4['N_err'][row_fit4][0]]
			v_lsr_list.append(fit4['v'][row_fit4])
		else:
			temp_fit4 = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
	else:
		temp_fit4 = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

	if (row_fit5.size == 1):
		if (fit5['good'][row_fit5] == 1):
			temp_fit5 = [fit5['T_L'][row_fit5][0], fit5['T_L_err'][row_fit5][0], fit5['width'][row_fit5][0], fit5['width_err'][row_fit5][0], fit5['v'][row_fit5][0], fit5['v_err'][row_fit5][0], fit5['tau'][row_fit5][0], fit5['tau_err'][row_fit5][0], fit5['N'][row_fit5][0], fit5['N_err'][row_fit5][0]]
			v_lsr_list.append(fit5['v'][row_fit5])
		else:
			temp_fit5 = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
	else:
		temp_fit5 = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

	if (row_fit6.size == 1):
		if (fit6['good'][row_fit6] == 1):
			temp_fit6 = [fit6['T_L'][row_fit6][0], fit6['T_L_err'][row_fit6][0], fit6['width'][row_fit6][0], fit6['width_err'][row_fit6][0], fit6['v'][row_fit6][0], fit6['v_err'][row_fit6][0], fit6['tau'][row_fit6][0], fit6['tau_err'][row_fit6][0], fit6['N'][row_fit6][0], fit6['N_err'][row_fit6][0]]
			v_lsr_list.append(fit6['v'][row_fit6])
		else:
			temp_fit6 = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
	else:
		temp_fit6 = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
	
	# get temperatures if value is reasonable
	temp_T12_rot     = T12_rot[y-1][x-1] if np.logical_and(T12_rot[y-1][x-1] < 500, T12_rot[y-1][x-1] > 0) else np.nan
	temp_T12_rot_err = T12_rot_err[y-1][x-1] if np.logical_and(T12_rot[y-1][x-1] < 500, T12_rot[y-1][x-1] > 0) else np.nan
	temp_T24_rot     = T24_rot[y-1][x-1] if np.logical_and(T24_rot[y-1][x-1] < 500, T24_rot[y-1][x-1] > 0) else np.nan
	temp_T24_rot_err = T24_rot_err[y-1][x-1] if np.logical_and(T24_rot[y-1][x-1] < 500, T24_rot[y-1][x-1] > 0) else np.nan
	temp_T45_rot     = T45_rot[y-1][x-1] if np.logical_and(T45_rot[y-1][x-1] < 500, T45_rot[y-1][x-1] > 0) else np.nan
	temp_T45_rot_err = T45_rot_err[y-1][x-1] if np.logical_and(T45_rot[y-1][x-1] < 500, T45_rot[y-1][x-1] > 0) else np.nan
	temp_T36_rot     = T36_rot[y-1][x-1] if np.logical_and(T36_rot[y-1][x-1] < 500, T36_rot[y-1][x-1] > 0) else np.nan
	temp_T36_rot_err = T36_rot_err[y-1][x-1] if np.logical_and(T36_rot[y-1][x-1] < 500, T36_rot[y-1][x-1] > 0) else np.nan

	temp_T12_kin     = T12_kin[y-1][x-1] if np.logical_and(T12_kin[y-1][x-1] < 60, T12_rot[y-1][x-1] > 0) else np.nan
	temp_T12_kin_err = T12_kin_err[y-1][x-1] if np.logical_and(T12_kin[y-1][x-1] < 60, T12_rot[y-1][x-1] > 0) else np.nan
	temp_T24_kin     = T24_kin[y-1][x-1] if np.logical_and(T24_kin[y-1][x-1] < 200, T24_rot[y-1][x-1] > 0) else np.nan
	temp_T24_kin_err = T24_kin_err[y-1][x-1] if np.logical_and(T24_kin[y-1][x-1] < 200, T24_rot[y-1][x-1] > 0) else np.nan
	temp_T45_kin     = T45_kin[y-1][x-1] if np.logical_and(T45_kin[y-1][x-1] < 200, T45_rot[y-1][x-1] > 0) else np.nan
	temp_T45_kin_err = T45_kin_err[y-1][x-1] if np.logical_and(T45_kin[y-1][x-1] < 200, T45_rot[y-1][x-1] > 0) else np.nan
	temp_T36_kin     = T36_kin[y-1][x-1] if np.logical_and(T36_kin[y-1][x-1] < 300, T36_rot[y-1][x-1] > 0) else np.nan
	temp_T36_kin_err = T36_kin_err[y-1][x-1] if np.logical_and(T36_kin[y-1][x-1] < 300, T36_rot[y-1][x-1] > 0) else np.nan
	
	x_HiGAL = int(round((glon-HiGAL_CRVAL1)/HiGAL_CDELT1+HiGAL_CRPIX1-1))
	y_HiGAL = int(round((glat-HiGAL_CRVAL2)/HiGAL_CDELT2+HiGAL_CRPIX2-1))
	temp_dust = dust[y_HiGAL, x_HiGAL]

	temp_temps = [temp_T12_rot, temp_T12_rot_err, temp_T24_rot, temp_T24_rot_err, temp_T45_rot, temp_T45_rot_err, temp_T36_rot, temp_T36_rot_err,
			temp_T12_kin, temp_T12_kin_err, temp_T24_kin, temp_T24_kin_err, temp_T45_kin, temp_T45_kin_err, temp_T36_kin, temp_T36_kin_err,
			temp_dust]
	
	
	
	# calculate time according to Kruijssen model
	#############################################
	
	# use median cloud velocity to match to sequence
	# returns np.nan when empty
	v_lsr_mean = np.median(v_lsr_list)
	
	# build a list of matching velocity
	# i.e. select the stream that is matching velocitywise
	selected = stream.copy()
	count = 0
	while ( count < len(selected) ):
		# select orbit positions that do NOT match velocity
		if not ( np.absolute(selected['v_los'][count]-v_lsr_mean) < 20.0 ):
			# delete them
			selected = np.delete(selected, count, 0)
		# if the line matches the velocity, jump to the next one
		else:
			count += 1

	# select position on the orbit with the least distance if possible
	if ( selected.size == 0 ):
		time = np.nan
	else:
		distance_list = []
		for j in np.arange(len(selected)):
			delta_l = glon - selected['l'][j]
			delta_b = glat - selected['b'][j]
			distance_list.append(np.sqrt(delta_l**2 + delta_b**2))
		min_index = np.argmin(distance_list)
		time = selected['t'][min_index]
	
	
	# concat temporary lists
	########################
	
	mega_row = [x, y, glon, glat, v_lsr_mean, time] + temp_fit1 + temp_fit2 + temp_fit3 + temp_fit4 + temp_fit5 + temp_fit6 + temp_temps
	
	
	# update list: append the current row to mega_list
	##################################################
	
	mega_list.append(mega_row)


#################################################################################


# save lists to disk
####################

header1 = "list of many many things :P\n"
header2 = "x\ty       glon       glat    v_mean    time |    (1,1) peak     (1,1) width     (1,1) v_lsr     (1,1) opacity  (1,1) column density |    (2,2) peak     (2,2) width     (2,2) v_lsr     (2,2) opacity  (2,2) column density |    (3,3) peak     (3,3) width     (3,3) v_lsr     (3,3) opacity  (3,3) column density |    (4,4) peak     (4,4) width     (4,4) v_lsr     (4,4) opacity  (4,4) column density |    (5,5) peak     (5,5) width     (5,5) v_lsr     (5,5) opacity  (5,5) column density |    (6,6) peak     (6,6) width     (6,6) v_lsr     (6,6) opacity  (6,6) column density |       T12 rot         T24 rot         T45 rot         T36 rot         T12 kin         T24 kin         T45 kin         T36 kin     Tdust"
header3 = "[pix]\t[pix]  [deg]      [deg]    [km/s]   [Myr] |     [K]   error   [km/s]  error  [km/s]   error     [1]   error    [cm^-2]      error |     [K]   error   [km/s]  error  [km/s]   error     [1]   error    [cm^-2]      error |     [K]   error   [km/s]  error  [km/s]   error     [1]   error    [cm^-2]      error |     [K]   error   [km/s]  error  [km/s]   error     [1]   error    [cm^-2]      error |     [K]   error   [km/s]  error  [km/s]   error     [1]   error    [cm^-2]      error |     [K]   error   [km/s]  error  [km/s]   error     [1]   error    [cm^-2]      error |     [K]   error     [K]   error     [K]   error     [K]   error     [K]   error     [K]   error     [K]   error     [K]   error     [K]"
header4 = '-'*714

print "\nsaving table to disk ...\n"
np.savetxt('NH3.pixel_fit.all.txt',
	mega_list,
	fmt = "%i\t%i %+9.7f %+9.7f %7.2f %7.2f | %7.3f %7.3f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %+10.3e %+10.3e | %7.3f %7.3f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %+10.3e %+10.3e | %7.3f %7.3f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %+10.3e %+10.3e | %7.3f %7.3f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %+10.3e %+10.3e | %7.3f %7.3f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %+10.3e %+10.3e | %7.3f %7.3f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %+10.3e %+10.3e | %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f",
	header    = header1+"\n"+header2+"\n"+header3+"\n"+header4
	)

