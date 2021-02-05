import os
import numpy as np

# ####################################
# # ******************************** #
# # * SF sequence analysis script  * #
# # ******************************** #
# ####################################


# Script to make a joint list of all fit parameters

# parameters:
# J = 1			# set in ipython console
line = 'NH3_'+str(J)+'-'+str(J)
fit_dir = 'fits_NH3/'
fit_log_dir = fit_dir+'fit_logs/'+line+'/'

FWHM_x = 26.22	# arcsec
FWHM_y = 17.84	# arcsec
NH3_data = np.array([[24.4, 65.6, 124.7, 201.7, 296.5, 409.2],
	[23.69447, 23.72263, 23.87013, 24.13942, 24.53292, 25.05596],
	[0.500, 0.796, 0.894, 0.935, 0.956, 0.969]])

############################################################

# make a list of all fits #
###########################

fit_list = []

for f in os.listdir('./fits_NH3/fit_logs/'+line):
	# print status
	print "processing "+f
	
	# get pixel coordinates
	x = int(f.split("_")[2])
	yend = f.split("_")[3]
	y = int(yend.split(".")[0])
	
	# load fit log
	fit_log = np.genfromtxt(fit_log_dir+f, 
		comments = "!",
		names = "nspecx, nspecy, xoff, yoff, area, area_err, v, v_err, width, width_err, tau, tau_err, base_rms, line_rms",
		dtype = ('<i4', '<i4', '|S8', '|S8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8')
		)
		
	# set flag to 0 for bad fit
	if ( fit_log['v_err'] > 10.0 or fit_log['width'] > 50.0 or fit_log['width_err'] > 10.0 ):
		fit_log['nspecx'] = 0
		fit_log['nspecy'] = 0
	
	
	
	# calculate further values from fit results #
	#############################################
	
	K    = J
	freq = NH3_data[2][J-1]
	nu   = NH3_data[1][J-1]	# GHz
	mu   = 1.472		# Debye
	
	Tant_tau       = fit_log['area']
	delta_Tant_tau = fit_log['area_err']
	tau            = fit_log['tau']
	delta_tau      = fit_log['tau_err']
	delta_v        = fit_log['width']
	delta_delta_v  = fit_log['width_err']
	
	
	# calculate T_L as in Schilke (Diplomarbeit)
	T_L = Tant_tau/tau*(1.-np.exp(-1.0*tau))*1.226e6/(nu**2*FWHM_x*FWHM_y)	# conversion to K needed because map is in Jy/beam
	delta_T_L = np.sqrt(((1-np.exp(-1.*tau))/tau*delta_Tant_tau)**2+((Tant_tau*np.exp(-1.*tau)*(tau-np.exp(tau)+1))/tau**2*delta_tau)**2)	# gaussian error of Tant_tau and tau
	
	# calculate integral T dv
	int_T_tau = np.sqrt(np.pi)/2./np.sqrt(np.log(2))*delta_v*T_L*tau/(freq*(1-np.exp(-1.*tau)))	# K km/s
	delta_int_T_tau = np.sqrt(np.pi)/2./np.sqrt(np.log(2))*np.sqrt((T_L*tau/(freq*(1-np.exp(-1.*tau)))*delta_delta_v)**2+(delta_v*tau/(freq*(1-np.exp(-1.*tau)))*delta_T_L)**2+(T_L*delta_v/freq*np.exp(tau)*(np.exp(tau)-tau-1)/(np.exp(tau)-1)**2*delta_tau)**2)	# gaussian error of delta_v, T_L and tau
	
	# calculate N
	N_tau = 1.6698e14*J*(J+1)/(K**2)/(mu**2)/nu*int_T_tau
	delta_N_tau = 1.6698e14*J*(J+1)/(K**2)/(mu**2)/nu*delta_int_T_tau

	# y-axis in Boltzmann plot is N/(2J+1)/g
	# statistical weight is g=1 for para-NH3 and g=2 for ortho-NH3
	if (J==3) or (J==6):
		g = 2
	else:
		g = 1
	
	yaxis = np.log10(N_tau/(2*J+1)/g)
	delta_yaxis = delta_N_tau/(N_tau*np.log(10))
	
	
	
	# update the full list #
	########################
	fit_list.append([x, 
		y, 
		fit_log['nspecx'].tolist(), 
		fit_log['area'].tolist(), 
		fit_log['area_err'].tolist(), 
		fit_log['v'].tolist(), 
		fit_log['v_err'].tolist(), 
		fit_log['width'].tolist(), 
		fit_log['width_err'].tolist(), 
		fit_log['tau'].tolist(), 
		fit_log['tau_err'].tolist(), 
		fit_log['base_rms'].tolist(), 
		fit_log['line_rms'].tolist(),
		'|',
		"{0:.3e}".format(T_L),
		"{0:.3e}".format(delta_T_L),
		"{0:.3e}".format(int_T_tau),
		"{0:.3e}".format(delta_int_T_tau),
		"{0:.3e}".format(N_tau),
		"{0:.3e}".format(delta_N_tau),
		"{0:.3e}".format(yaxis),
		"{0:.3e}".format(delta_yaxis)]
		)
	numpy_fit_list = np.asarray(fit_list)
	
# write full list to disk
np.savetxt(fit_dir+line+'.fit_list.txt',
	numpy_fit_list,
	fmt = "%s",
	delimiter = "\t",
	header = "list of all fitted pixels, some fits are good, others are flagged as bad\nleft values are fit results, right hand values are calculated from them\n\npixel coo.\tfit\tline area\tv_lsr\t\tline width\topacity\t\tCLASSfit rms\t|\tline temperature\t\tintegral Tdv\t\t\tcolumn density\t\t\tN/(2J+1)/g (yaxis)\nx\ty\tgood?\tvalue\terror\tvalue\terror\tvalue\terror\tvalue\terror\tbaserms\tlinerms\t|\tvalue\t\terror\t\tvalue\t\terror\t\tvalue\t\terror\t\tvalue\t\terror\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	)
print "wrote "+fit_dir+line+".fit_list.txt to disk"
