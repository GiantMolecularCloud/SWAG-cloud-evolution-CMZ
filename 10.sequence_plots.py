from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as lines
from matplotlib import rc
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import LogNorm
import os
rc('text', usetex=True)
plt.ioff()

# Define function for LaTeX string formatting of scientific notation
def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
	if not exponent:
		exponent = int(np.floor(np.log10(abs(num))))
	coeff = round(num / float(10**exponent), decimal_digits)
	if not precision:
		precision = decimal_digits
	return r"{0:.{2}f}\cdot10^{{{1:d}}}".format(coeff, exponent, precision)


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

data = np.genfromtxt('NH3.pixel_fit.all.txt',
#	usecols = (0,1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22,23,24,25,26,27,29,30,31,32,33,34,35,36,37,38,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,62,63,64,65,66,67,68,69,70,71,73,74,75,76,77,78),
	dtype   = "i4, i4, f8, f8, f8, f8, |S10, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, |S10, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, |S10, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, |S10, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, |S10, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, |S10, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, |S10, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8",
	names   = ('x','y','glon','glat','v_mean','time','|1','peak1','peak1_err','width1','width1_err','v1','v1_err','tau1','tau1_err','N1','N1_err','|2','peak2','peak2_err','width2','width2_err','v2','v2_err','tau2','tau2_err','N2','N2_err','|3','peak3','peak3_err','width3','width3_err','v3','v3_err','tau3','tau3_err','N3','N3_err','|4','peak4','peak4_err','width4','width4_err','v4','v4_err','tau4','tau4_err','N4','N4_err','|5','peak5','peak5_err','width5','width5_err','v5','v5_err','tau5','tau5_err','N5','N5_err','|6','peak6','peak6_err','width6','width6_err','v6','v6_err','tau6','tau6_err','N6','N6_err','|7','T12_rot','T12_rot_err','T24_rot','T24_rot_err','T45_rot','T45_rot_err','T36_rot','T36_rot_err','T12_kin','T12_kin_err','T24_kin','T24_kin_err','T45_kin','T45_kin_err','T36_kin','T36_kin_err','Tdust')
	)

# calculate "phase" = time since last pericenter passage
# pericenter to pericenter: 2.03 Myr
phase = data['time']%2.03
phase_shifted = np.copy(phase)
for row in np.arange(len(phase_shifted)):
	if (phase_shifted[row] > 1.5):
		phase_shifted[row] += -2.03

# load Kruijssen orbit
model = np.genfromtxt('Kruijssen_orbit.dat',
	dtype = 'f8',
	names = ('t', 'x', 'y', 'z', 'R', 'vx', 'vy', 'vz', 'vorb', 'l', 'b', 'vlos', 'mu_l', 'mu_b', 'mu_x', 'mu_y')
	)

temp_info = [["12","24","45","36"],[80,250,250,360.0],[60,200,200,300]]


############################################################

#os.system('mkdir plots/')

# set colormap
cmap = cm.rainbow
vmin = -100.0
vmax = 100.0

# order plot objects
scatter_kwargs = {"zorder":8}
fit_kwargs     = {"zorder":10}
text_kwargs    = {"zorder":12}
limit_kwargs   = {"zorder":6}

props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)




# print statistics
##################

# all pixels in this list
npix = float(len(data))

# pixels that were found to fit to the Kruijssen orbit
match1 = np.count_nonzero(~np.isnan(data['N1']))
match2 = np.count_nonzero(~np.isnan(data['N2']))
match3 = np.count_nonzero(~np.isnan(data['N3']))
match4 = np.count_nonzero(~np.isnan(data['N4']))
match5 = np.count_nonzero(~np.isnan(data['N5']))
match6 = np.count_nonzero(~np.isnan(data['N6']))

# pixels that have a temperature
match_rot_T12 = np.count_nonzero(~np.isnan(data['T12_rot']))
match_rot_T24 = np.count_nonzero(~np.isnan(data['T24_rot']))
match_rot_T45 = np.count_nonzero(~np.isnan(data['T45_rot']))
match_rot_T36 = np.count_nonzero(~np.isnan(data['T36_rot']))
match_kin_T12 = np.count_nonzero(~np.isnan(data['T12_kin']))
match_kin_T24 = np.count_nonzero(~np.isnan(data['T24_kin']))
match_kin_T45 = np.count_nonzero(~np.isnan(data['T45_kin']))
match_kin_T36 = np.count_nonzero(~np.isnan(data['T36_kin']))

print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print " NH3 fit results"
print " all pixels in list: "+str('{:.1f}'.format(npix))
print " \t\t(1,1)\t(2,2)\t(3,3)\t(4,4)\t(5,5)\t(6,6)"
print "-------------------------------------------------------------"
print " match orbit:\t"+str('{:.0f}'.format(match1))+"\t"+str('{:.0f}'.format(match2))+"\t"+str('{:.0f}'.format(match3))+"\t"+str('{:.0f}'.format(match4))+"\t"+str('{:.0f}'.format(match5))+"\t"+str('{:.0f}'.format(match6))
print " % match orbit:\t"+str('{:.1f}'.format(match1/npix*100.0))+"\t"+str('{:.1f}'.format(match2/npix*100.0))+"\t"+str('{:.1f}'.format(match3/npix*100.0))+"\t"+str('{:.1f}'.format(match4/npix*100.0))+"\t"+str('{:.1f}'.format(match5/npix*100.0))+"\t"+str('{:.1f}'.format(match6/npix*100.0))
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print " Tij =\t\tT12\tT24\tT45\tT36"
print "-------------------------------------------------------------"
print " derived Trot:\t"+str('{:.0f}'.format(match_rot_T12))+"\t"+str('{:.0f}'.format(match_rot_T24))+"\t"+str('{:.0f}'.format(match_rot_T45))+"\t"+str('{:.0f}'.format(match_rot_T36))
print " derived Tkin:\t"+str('{:.0f}'.format(match_kin_T12))+"\t"+str('{:.0f}'.format(match_kin_T24))+"\t"+str('{:.0f}'.format(match_kin_T45))+"\t"+str('{:.0f}'.format(match_kin_T36))
print " % Tkin/npix:\t"+str('{:.1f}'.format(match_kin_T12/npix*100.0))+"\t"+str('{:.1f}'.format(match_kin_T24/npix*100.0))+"\t"+str('{:.1f}'.format(match_kin_T45/npix*100.0))+"\t"+str('{:.1f}'.format(match_kin_T36/npix*100.0))
print " % Tkin/matchj:\t"+str('{:.1f}'.format(match_kin_T12/match2*100.0))+"\t"+str('{:.1f}'.format(match_kin_T24/match4*100.0))+"\t"+str('{:.1f}'.format(match_kin_T45/match5*100.0))+"\t"+str('{:.1f}'.format(match_kin_T36/match6*100.0))
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"







# plot
######


############################################################

# T?? vs. time
print "plotting T(NH3) vs. time ..."
for ttype in ["kin", "rot"]:
	for temp in ["12", "24", "45", "36"]:
		plt.clf()
		plt.cla()
		fig = plt.figure()
		fig.subplots_adjust(top=0.8)
		ax1 = plt.subplot(1,1,1)
		ax2 = ax1.twiny()
		ax1.set_xlabel('time [Myr]')
		ax1.set_ylabel('NH$_3$ temperature T [K]')
		ax1.set_xlim(-2.5,2.5)
		ax1.xaxis.set_minor_locator(MultipleLocator(0.2))
		ax1.set_ylim(0,400)
		if (temp == "12") and (ttype == "kin"): 
			ax1.set_ylim(0,80)
			ax1.text(0.17,6, r'$\tau_{ff} = 0.34$\,Myr', color='r', ha='center', va='bottom', **text_kwargs)
			ax1.plot([0.0,0.34], [5,5], color='r', linewidth=2, **text_kwargs)
		if (temp == "24") and (ttype == "kin"): 
			ax1.set_ylim(0,280)
			ax1.text(0.17,12, r'$\tau_{ff} = 0.34$\,Myr', color='r', ha='center', va='bottom', **text_kwargs)
			ax1.plot([0.0,0.34], [10,10], color='r', linewidth=2, **text_kwargs)
		if (temp == "45") and (ttype == "kin"): 
			ax1.set_ylim(0,280)
			ax1.text(0.17,12, r'$\tau_{ff} = 0.34$\,Myr', color='r', ha='center', va='bottom', **text_kwargs)
			ax1.plot([0.0,0.34], [10,10], color='r', linewidth=2, **text_kwargs)
		if (temp == "36") and (ttype == "kin"): 
			ax1.set_ylim(0,400)
			ax1.text(0.17,17, r'$\tau_{ff} = 0.34$\,Myr', color='r', ha='center', va='bottom', **text_kwargs)
			ax1.plot([0.0,0.34], [15,15], color='r', linewidth=2, **text_kwargs)
		ax1.yaxis.set_minor_locator(MultipleLocator(10))
		ax2.set_xlim(-2.5,2.5)
		ax2.set_xticks([-2.0, -1.67, -1.52, -1.30, -1.0, 0.0, 1.0, 1.51, 2.0])
		ax2.set_xticklabels(['pericenter\nnear', 'Brick', 'cloud e/f', 'Sgr B2', 'apocenter\neast', 'pericenter\nfar', 'apocenter\nwest', 'Sgr C', 'pericenter\nnear'], rotation=90, ha='center', va='bottom')
		ax1.bar(-1.7, 500, 0.4, bottom=None, color="grey", alpha=0.5, zorder=3)
		ax1.text(-1.5, 450, "dust ridge", ha='center', bbox=props, **text_kwargs)
	
		# print results
		print "T"+temp+","+ttype+" slopes:"
	
		# fit sequences
		fit_ranges = [[-2.0,-1.0],[-0.3,0.05],[1.75,2.25]]
		fit_results = "slopes:\n"
		for fit_range_index in np.arange(len(fit_ranges)):
			fit_lim_lo = fit_ranges[fit_range_index][0]
			fit_lim_hi = fit_ranges[fit_range_index][1]
	
			# fit sequence
			to_fit = [[],[]]
			for row in np.arange(len(data)):
				if (data['time'][row] > fit_lim_lo) and (data['time'][row] < fit_lim_hi) and not (np.isnan(data['T'+temp+'_'+ttype][row])): 
					to_fit[0].append(data['time'][row])
					to_fit[1].append(data['T'+temp+'_'+ttype][row])
#			coeff, cov = np.polyfit(to_fit[0], to_fit[1], 1, cov=True)
#			ax1.plot([-2.0,-1.0],[x*coeff[0]+coeff[1] for x in [-2.0,-1.0]], linewidth=1.0, markersize=0.0, color='black', linestyle='--', label='fit $'+str(fit_lim_lo)+'<t<'+str(fit_lim_hi)+'$', **fit_kwargs)
			
			# fit binned sequence
			bin_fit = [[],[],[],[]]
			for time_bin in np.arange(fit_lim_lo,fit_lim_hi,0.01):
				bin_value = []
				for row in np.arange(len(to_fit[0])):
					if (to_fit[0][row] >= time_bin) and (to_fit[0][row] < time_bin+0.01):
						bin_value.append(to_fit[1][row])
				bin_value.sort()
				if (np.median(bin_value) > 5) and (len(bin_value) > 10):
					bin_fit[0].append(time_bin+0.01)
					bin_fit[1].append(np.median(bin_value))
					bin_fit[2].append(np.median(bin_value) - bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][0])	# lower "error" of 50% include range
					bin_fit[3].append(bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][-1] - np.median(bin_value))	# upper "error" of 50% include range
			coeff_med, cov_med = np.polyfit(bin_fit[0], bin_fit[1], 1, cov=True)
			fit_results += r"$\frac{dT}{dt} = "+str('{:3.1f}'.format(coeff_med[0])).rjust(5)+"$\\,K/Myr ($"+str(fit_lim_lo)+"<t [Myr] <"+str(fit_lim_hi)+"$)\n"
			
			# plot fits; only one legend entry
			if (fit_range_index == 0):
				ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', label='fit bin median', **fit_kwargs)
				ax1.scatter(bin_fit[0], bin_fit[1], label='0.01\,Myr bins median', s=10, marker="x", color='black', **text_kwargs)
#				ax1.errorbar(bin_fit[0], bin_fit[1], yerr=zip(*[bin_fit[2],bin_fit[3]]), fmt='none', linewidth=0.0, elinewidth=0.5, ecolor='black', **limit_kwargs)
			else:
				ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', **fit_kwargs)
				ax1.scatter(bin_fit[0], bin_fit[1], s=10, marker="x", color='black', **text_kwargs)
#				ax1.errorbar(bin_fit[0], bin_fit[1], yerr=zip(*[bin_fit[2],bin_fit[3]]), fmt='none', linewidth=0.0, elinewidth=0.5, ecolor='black', **limit_kwargs)
			
			print "["+str(fit_lim_lo)+","+str(fit_lim_hi)+"]\t"+str(coeff_med[0])+" +- "+str(cov_med[0][0])
			
		# fit results text box
		ax1.text(0.68, 0.975, fit_results[:-1], transform=ax1.transAxes, fontsize=10, va='top', ha='right', multialignment='left', bbox=props)
		
		# plot data
		scatter = ax1.scatter(data['time'], data['T'+temp+'_'+ttype], label='T$_{'+ttype+','+temp+'}$', s=1, lw=0, c=data['v_mean'], cmap=cmap, vmin=vmin, vmax=vmax, **scatter_kwargs)
		
		# set Tkin sensitivity limit
		if (temp == "12") and (ttype == "kin"): 
			ax1.axhline(y=60, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		if (temp == "24") and (ttype == "kin"): 
			ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		if (temp == "45") and (ttype == "kin"): 
			ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		if (temp == "36") and (ttype == "kin"): 
			ax1.axhline(y=300, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)

		# color bar
		cbaxes = fig.add_axes([0.9, 0.1, 0.02, 0.7]) 
		cbar = plt.colorbar(scatter, cax=cbaxes)
		cbar.set_label("v$_{los}$ [km/s]")
		
		leg = ax1.legend(fancybox=True, framealpha=0.8, fontsize=10)
		leg.set_zorder(12)
		
		plt.savefig("plots/time_vs_T"+temp+"_"+ttype+".png", dpi=300)


############################################################

# T?? vs. time (zoom in)
print "plotting T(NH3) vs. time ..."
for ttype in ["kin", "rot"]:
	for temp in ["12", "24", "45", "36"]:
		fit_ranges = [[-2.0,-1.0],[-0.3,0.05],[1.75,2.25]]
		for fit_range_index in np.arange(len(fit_ranges)):
			fit_lim_lo = fit_ranges[fit_range_index][0]
			fit_lim_hi = fit_ranges[fit_range_index][1]
			fit_results = "slope:\n"

			plt.clf()
			plt.cla()
			fig = plt.figure()
			fig.subplots_adjust(top=0.8)
			ax1 = plt.subplot(1,1,1)
			ax2 = ax1.twiny()
			ax1.set_xlabel('time [Myr]')
			ax1.set_ylabel('NH$_3$ temperature T [K]')
			ax1.set_xlim(fit_lim_lo,fit_lim_hi)
			ax1.xaxis.set_minor_locator(MultipleLocator(0.1))
			ax1.set_ylim(0,400)
			if (temp == "12") and (ttype == "kin"): 
				ax1.set_ylim(0,80)
			if (temp == "24") and (ttype == "kin"): 
				ax1.set_ylim(0,280)
			if (temp == "45") and (ttype == "kin"): 
				ax1.set_ylim(0,280)
			if (temp == "36") and (ttype == "kin"): 
				ax1.set_ylim(0,400)
			ax1.yaxis.set_minor_locator(MultipleLocator(10))
			ax2.set_xticks([-2.0, -1.67, -1.52, -1.30, -1.0, 0.0, 1.0, 1.51, 2.0])
			ax2.set_xticklabels(['pericenter\nnear', 'Brick', 'cloud e/f', 'Sgr B2', 'apocenter\neast', 'pericenter\nfar', 'apocenter\nwest', 'Sgr C', 'pericenter\nnear'], rotation=90, ha='center', va='bottom')
			ax2.set_xlim(fit_lim_lo,fit_lim_hi)
#			ax1.bar(-1.7, 500, 0.4, bottom=None, color="grey", alpha=0.5, zorder=3)
#			ax1.text(-1.5, 450, "dust ridge", ha='center', bbox=props, **text_kwargs)
			ax1.text(0.17,12, r'$\tau_{ff} = 0.34$\,Myr', color='r', ha='center', va='bottom', **text_kwargs)
			ax1.plot([0.0,0.34], [10,10], color='r', linewidth=2, **text_kwargs)

	
			# fit sequence
			to_fit = [[],[]]
			for row in np.arange(len(data)):
				if (data['time'][row] > fit_lim_lo) and (data['time'][row] < fit_lim_hi) and not (np.isnan(data['T'+temp+'_'+ttype][row])): 
					to_fit[0].append(data['time'][row])
					to_fit[1].append(data['T'+temp+'_'+ttype][row])
#			coeff, cov = np.polyfit(to_fit[0], to_fit[1], 1, cov=True)
#			ax1.plot([fit_lim_lo,-1.0],[x*coeff[0]+coeff[1] for x in [fit_lim_lo,-1.0]], linewidth=1.0, markersize=0.0, color='black', linestyle='--', label='fit $'+str(fit_lim_lo)+'<t<'+str(fit_lim_hi)+'$', **fit_kwargs)
			
			# fit binned sequence
			bin_fit = [[],[],[],[]]
			for time_bin in np.arange(fit_lim_lo,fit_lim_hi,0.01):
				bin_value = []
				for row in np.arange(len(to_fit[0])):
					if (to_fit[0][row] >= time_bin) and (to_fit[0][row] < time_bin+0.01):
						bin_value.append(to_fit[1][row])
				bin_value.sort()
				if (np.median(bin_value) > 5) and (len(bin_value) > 10):
					bin_fit[0].append(time_bin+0.01)
					bin_fit[1].append(np.median(bin_value))
					bin_fit[2].append(np.median(bin_value) - bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][0])	# lower "error" of 50% include range
					bin_fit[3].append(bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][-1] - np.median(bin_value))	# upper "error" of 50% include range
			coeff_med, cov_med = np.polyfit(bin_fit[0], bin_fit[1], 1, cov=True)
			fit_results += r"$\frac{dT}{dt} = "+str('{:3.1f}'.format(coeff_med[0])).rjust(5)+"$\\,K/Myr ($"+str(fit_lim_lo)+"<t [Myr] <"+str(fit_lim_hi)+"$)\n"
			
			# plot fits
			ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', label='fit bin median', **fit_kwargs)
			ax1.scatter(bin_fit[0], bin_fit[1], label='0.01\,Myr bins median', s=10, marker="x", color='black', **text_kwargs)
#			ax1.errorbar(bin_fit[0], bin_fit[1], yerr=zip(*[bin_fit[2],bin_fit[3]]), fmt='none', linewidth=0.0, elinewidth=0.5, ecolor='black', **limit_kwargs)
			
			# fit results text box
			ax1.text(0.68, 0.975, fit_results[:-1], transform=ax1.transAxes, fontsize=10, va='top', ha='right', multialignment='left', bbox=props)
	
			# plot data
			scatter = ax1.scatter(data['time'], data['T'+temp+'_'+ttype], label='T$_{'+ttype+','+temp+'}$', s=3, lw=0, c=data['v_mean'], cmap=cmap, vmin=vmin, vmax=vmax, **scatter_kwargs)
			
			# set Tkin sensitivity limit
			if (temp == "12") and (ttype == "kin"): 
				ax1.axhline(y=60, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			if (temp == "24") and (ttype == "kin"): 
				ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			if (temp == "45") and (ttype == "kin"): 
				ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			if (temp == "36") and (ttype == "kin"): 
				ax1.axhline(y=300, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		
			# color bar
			cbaxes = fig.add_axes([0.9, 0.1, 0.02, 0.7]) 
			cbar = plt.colorbar(scatter, cax=cbaxes)
			cbar.set_label("v$_{los}$ [km/s]")
	
			leg = ax1.legend(fancybox=True, framealpha=0.8, fontsize=10)
			leg.set_zorder(12)
			
			plt.savefig("plots/time_vs_T"+temp+"_"+ttype+".zoom_"+str(fit_lim_lo)+"_"+str(fit_lim_hi)+".png", dpi=300)


############################################################

# T?? vs. time (density)
print "plotting T(NH3) vs. time ..."
for ttype in ["kin", "rot"]:
	for temp in ["12", "24", "45", "36"]:
		plt.clf()
		plt.cla()
		fig = plt.figure()
		fig.subplots_adjust(top=0.8)
		ax1 = plt.subplot(1,1,1)
		ax2 = ax1.twiny()
		ax1.set_xlabel('time [Myr]')
		ax1.set_ylabel('NH$_3$ temperature T [K]')
		ax1.set_xlim(-2.5,2.5)
		ax1.xaxis.set_minor_locator(MultipleLocator(0.2))
		ax1.set_ylim(0, 400)
		if (temp == "12"):
			ax1.set_ylim(0,80)
		if (temp == "24"):
			ax1.set_ylim(0,280)
		if (temp == "45"):
			ax1.set_ylim(0,280)
		if (temp == "36"):
			ax1.set_ylim(0,400)
		ax1.yaxis.set_minor_locator(MultipleLocator(10))
		ax2.set_xlim(-2.5,2.5)
		ax2.set_xticks([-2.0, -1.67, -1.52, -1.30, -1.0, 0.0, 1.0, 1.51, 2.0])
		ax2.set_xticklabels(['pericenter\nnear', 'Brick', 'cloud e/f', 'Sgr B2', 'apocenter\neast', 'pericenter\nfar', 'apocenter\nwest', 'Sgr C', 'pericenter\nnear'], rotation=90, ha='center', va='bottom')
		ax1.bar(-1.7, 500, 0.4, bottom=None, color="grey", alpha=0.5, zorder=3)
		ax1.text(-1.5, 380, "dust ridge", ha='center', bbox=props, **text_kwargs)
		ax1.text(0.17,12, r'$\tau_{ff} = 0.34$\,Myr', color='r', ha='center', va='bottom', **text_kwargs)
		ax1.plot([0.0,0.34], [10,10], color='r', linewidth=2, **text_kwargs)

	
		# fit sequences
		fit_ranges = [[-2.0,-1.0],[-0.3,0.05],[1.75,2.25]]
		fit_results = "slopes:\n"
		for fit_range_index in np.arange(len(fit_ranges)):
			fit_lim_lo = fit_ranges[fit_range_index][0]
			fit_lim_hi = fit_ranges[fit_range_index][1]
	
			# fit sequence
			to_fit = [[],[]]
			for row in np.arange(len(data)):
				if (data['time'][row] > fit_lim_lo) and (data['time'][row] < fit_lim_hi) and not (np.isnan(data['T'+temp+'_'+ttype][row])): 
					to_fit[0].append(data['time'][row])
					to_fit[1].append(data['T'+temp+'_'+ttype][row])
#			coeff, cov = np.polyfit(to_fit[0], to_fit[1], 1, cov=True)
#			ax1.plot([-2.0,-1.0],[x*coeff[0]+coeff[1] for x in [-2.0,-1.0]], linewidth=1.0, markersize=0.0, color='black', linestyle='--', label='fit $'+str(fit_lim_lo)+'<t<'+str(fit_lim_hi)+'$', **fit_kwargs)
			
			# fit binned sequence
			bin_fit = [[],[],[],[]]
			for time_bin in np.arange(fit_lim_lo,fit_lim_hi,0.01):
				bin_value = []
				for row in np.arange(len(to_fit[0])):
					if (to_fit[0][row] >= time_bin) and (to_fit[0][row] < time_bin+0.01):
						bin_value.append(to_fit[1][row])
				bin_value.sort()
				if (np.median(bin_value) > 5) and (len(bin_value) > 10):
					bin_fit[0].append(time_bin+0.01)
					bin_fit[1].append(np.median(bin_value))
					bin_fit[2].append(np.median(bin_value) - bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][0])	# lower "error" of 50% include range
					bin_fit[3].append(bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][-1] - np.median(bin_value))	# upper "error" of 50% include range
			coeff_med, cov_med = np.polyfit(bin_fit[0], bin_fit[1], 1, cov=True)
			fit_results += r"$\frac{dT}{dt} = "+str('{:3.1f}'.format(coeff_med[0])).rjust(5)+"$\\,K/Myr ($"+str(fit_lim_lo)+"<t [Myr] <"+str(fit_lim_hi)+"$)\n"
			
			# plot fits; only one legend entry
			if (fit_range_index == 0):
				ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', label='fit bin median', **fit_kwargs)
				ax1.scatter(bin_fit[0], bin_fit[1], label='0.01\,Myr bins median', s=10, marker="x", color='black', **text_kwargs)
#				ax1.errorbar(bin_fit[0], bin_fit[1], yerr=zip(*[bin_fit[2],bin_fit[3]]), fmt='none', linewidth=0.0, elinewidth=0.5, ecolor='black', **limit_kwargs)
			else:
				ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', **fit_kwargs)
				ax1.scatter(bin_fit[0], bin_fit[1], s=10, marker="x", color='black', **text_kwargs)
#				ax1.errorbar(bin_fit[0], bin_fit[1], yerr=zip(*[bin_fit[2],bin_fit[3]]), fmt='none', linewidth=0.0, elinewidth=0.5, ecolor='black', **limit_kwargs)
			
		# fit results text box
		ax1.text(0.68, 0.975, fit_results[:-1], transform=ax1.transAxes, fontsize=10, va='top', ha='right', multialignment='left', bbox=props)
		
		# plot data
		scatter = ax1.scatter(data['time'], data['T'+temp+'_'+ttype], label='T$_{'+ttype+','+temp+'}$', s=1, lw=0, c='dimgrey', **scatter_kwargs)
		
		# contour plot
		contour_list = [[],[]]
		for row in np.arange(len(data)):
			if not (np.isnan(data['time'][row])) and not (np.isnan(data['T'+temp+'_'+ttype][row])):
				contour_list[0].append(data['time'][row])
				contour_list[1].append(data['T'+temp+'_'+ttype][row])

		H, xedges, yedges = np.histogram2d(contour_list[0], contour_list[1], bins=100, range=[[-2.5,2.5],[0,500]])
		extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
#		levels=[25,50,75,125,175]
		levels = [i*np.amax(H) for i in np.arange(0.1,1.0,0.1)]
		contours = ax1.contour(np.transpose(H),extent=extent, levels=levels, label="10% steps of max. density", **scatter_kwargs)

		# set Tkin sensitivity limit
		if (temp == "12") and (ttype == "kin"): 
			ax1.axhline(y=60, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		if (temp == "24") and (ttype == "kin"): 
			ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		if (temp == "45") and (ttype == "kin"): 
			ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		if (temp == "36") and (ttype == "kin"): 
			ax1.axhline(y=300, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)

		# legend
		contour_legend = lines.Line2D([],[], color='b', lw=1, label='10% steps of max. density')
		leg = ax1.legend(fancybox=True, framealpha=0.8, fontsize=10)
		leg.set_zorder(12)
		
		plt.savefig("plots/time_vs_T"+temp+"_"+ttype+".density.png", dpi=300)


############################################################

# T?? vs. time (zoom in, density)
print "plotting T(NH3) vs. time ..."
for ttype in ["kin", "rot"]:
	for temp in ["12", "24", "45", "36"]:
		fit_ranges = [[-2.0,-1.0],[-0.3,0.05],[1.75,2.25]]
		for fit_range_index in np.arange(len(fit_ranges)):
			fit_lim_lo = fit_ranges[fit_range_index][0]
			fit_lim_hi = fit_ranges[fit_range_index][1]
			fit_results = "slope:\n"

			plt.clf()
			plt.cla()
			fig = plt.figure()
			fig.subplots_adjust(top=0.8)
			ax1 = plt.subplot(1,1,1)
			ax2 = ax1.twiny()
			ax1.set_xlabel('time [Myr]')
			ax1.set_ylabel('NH$_3$ temperature T [K]')
			ax1.set_xlim(fit_lim_lo,fit_lim_hi)
			ax1.xaxis.set_minor_locator(MultipleLocator(0.1))
			ax1.set_ylim(0,400)
			if (temp == "12") and (ttype == "kin"): 
				ax1.set_ylim(0,80)
			if (temp == "24") and (ttype == "kin"): 
				ax1.set_ylim(0,280)
			if (temp == "45") and (ttype == "kin"): 
				ax1.set_ylim(0,280)
			if (temp == "36") and (ttype == "kin"): 
				ax1.set_ylim(0,400)
			ax1.yaxis.set_minor_locator(MultipleLocator(10))
			ax2.set_xticks([-2.0, -1.67, -1.52, -1.30, -1.0, 0.0, 1.0, 1.51, 2.0])
			ax2.set_xticklabels(['pericenter\nnear', 'Brick', 'cloud e/f', 'Sgr B2', 'apocenter\neast', 'pericenter\nfar', 'apocenter\nwest', 'Sgr C', 'pericenter\nnear'], rotation=90, ha='center', va='bottom')
			ax2.set_xlim(fit_lim_lo,fit_lim_hi)
#			ax1.bar(-1.7, 500, 0.4, bottom=None, color="grey", alpha=0.5, zorder=3)
#			ax1.text(-1.5, 450, "dust ridge", ha='center', bbox=props, **text_kwargs)
			
			# fit sequence
			to_fit = [[],[]]
			for row in np.arange(len(data)):
				if (data['time'][row] > fit_lim_lo) and (data['time'][row] < fit_lim_hi) and not (np.isnan(data['T'+temp+'_'+ttype][row])): 
					to_fit[0].append(data['time'][row])
					to_fit[1].append(data['T'+temp+'_'+ttype][row])
#			coeff, cov = np.polyfit(to_fit[0], to_fit[1], 1, cov=True)
#			ax1.plot([fit_lim_lo,-1.0],[x*coeff[0]+coeff[1] for x in [fit_lim_lo,-1.0]], linewidth=1.0, markersize=0.0, color='black', linestyle='--', label='fit $'+str(fit_lim_lo)+'<t<'+str(fit_lim_hi)+'$', **fit_kwargs)
			
			# fit binned sequence
			bin_fit = [[],[],[],[]]
			for time_bin in np.arange(fit_lim_lo,fit_lim_hi,0.01):
				bin_value = []
				for row in np.arange(len(to_fit[0])):
					if (to_fit[0][row] >= time_bin) and (to_fit[0][row] < time_bin+0.01):
						bin_value.append(to_fit[1][row])
				bin_value.sort()
				if (np.median(bin_value) > 5) and (len(bin_value) > 10):
					bin_fit[0].append(time_bin+0.01)
					bin_fit[1].append(np.median(bin_value))
					bin_fit[2].append(np.median(bin_value) - bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][0])	# lower "error" of 50% include range
					bin_fit[3].append(bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][-1] - np.median(bin_value))	# upper "error" of 50% include range
			coeff_med, cov_med = np.polyfit(bin_fit[0], bin_fit[1], 1, cov=True)
			fit_results += r"$\frac{dT}{dt} = "+str('{:3.1f}'.format(coeff_med[0])).rjust(5)+"$\\,K/Myr ($"+str(fit_lim_lo)+"<t [Myr] <"+str(fit_lim_hi)+"$)\n"
			
			# plot fits
			ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', label='fit bin median', **fit_kwargs)
			ax1.scatter(bin_fit[0], bin_fit[1], label='0.01\,Myr bins median', s=10, marker="x", color='black', **text_kwargs)
#			ax1.errorbar(bin_fit[0], bin_fit[1], yerr=zip(*[bin_fit[2],bin_fit[3]]), fmt='none', linewidth=0.0, elinewidth=0.5, ecolor='black', **limit_kwargs)
			
			# fit results text box
			ax1.text(0.68, 0.975, fit_results[:-1], transform=ax1.transAxes, fontsize=10, va='top', ha='right', multialignment='left', bbox=props)
	
			# plot data
			scatter = ax1.scatter(data['time'], data['T'+temp+'_'+ttype], label='T$_{'+ttype+','+temp+'}$', s=1, lw=0, c='dimgrey', **scatter_kwargs)
			
			# contour plot
			contour_list = [[],[]]
			for row in np.arange(len(data)):
				if not (np.isnan(data['time'][row])) and not (np.isnan(data['T'+temp+'_'+ttype][row])):
					contour_list[0].append(data['time'][row])
					contour_list[1].append(data['T'+temp+'_'+ttype][row])
	
			H, xedges, yedges = np.histogram2d(contour_list[0], contour_list[1], bins=100, range=[[-2.5,2.5],[0,500]])
			extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
#			levels=[25,50,75,125,175]
			levels = [i*np.amax(H) for i in np.arange(0.1,1.0,0.1)]
			contours = ax1.contour(np.transpose(H),extent=extent, levels=levels, label="10% steps of max. density", **scatter_kwargs)

			
			# set Tkin sensitivity limit
			if (temp == "12") and (ttype == "kin"): 
				ax1.axhline(y=60, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			if (temp == "24") and (ttype == "kin"): 
				ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			if (temp == "45") and (ttype == "kin"): 
				ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			if (temp == "36") and (ttype == "kin"): 
				ax1.axhline(y=300, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		
			# legend
			contour_legend = lines.Line2D([],[], color='b', lw=1, label='10% steps of max. density')
			leg = ax1.legend(fancybox=True, framealpha=0.8, fontsize=10)
			leg.set_zorder(12)
			
			plt.savefig("plots/time_vs_T"+temp+"_"+ttype+".density.zoom_"+str(fit_lim_lo)+"_"+str(fit_lim_hi)+".png", dpi=300)


############################################################

# T?? vs. time (density only)
print "plotting T(NH3) vs. time ..."
for ttype in ["kin", "rot"]:
	for temp in ["12", "24", "45", "36"]:
		plt.clf()
		plt.cla()
		fig = plt.figure()
		fig.subplots_adjust(top=0.8)
		ax1 = plt.subplot(1,1,1)
		ax2 = ax1.twiny()
		ax1.set_xlabel('time [Myr]')
		ax1.set_ylabel('NH$_3$ temperature T [K]')
		ax1.set_xlim(-2.5,2.5)
		ax1.xaxis.set_minor_locator(MultipleLocator(0.2))
		ax1.set_ylim(0, 500)
		if (temp == "12"):
			ax1.set_ylim(0,80)
		if (temp == "24"):
			ax1.set_ylim(0,280)
		if (temp == "45"):
			ax1.set_ylim(0,280)
		if (temp == "36"):
			ax1.set_ylim(0,400)
		ax1.yaxis.set_minor_locator(MultipleLocator(10))
		ax2.set_xlim(-2.5,2.5)
		ax2.set_xticks([-2.0, -1.67, -1.52, -1.30, -1.0, 0.0, 1.0, 1.51, 2.0])
		ax2.set_xticklabels(['pericenter\nnear', 'Brick', 'cloud e/f', 'Sgr B2', 'apocenter\neast', 'pericenter\nfar', 'apocenter\nwest', 'Sgr C', 'pericenter\nnear'], rotation=90, ha='center', va='bottom')
		ax1.bar(-1.7, 500, 0.4, bottom=None, color="grey", alpha=0.5, zorder=3)
		ax1.text(-1.5, 450, "dust ridge", ha='center', bbox=props, **text_kwargs)
		ax1.text(0.17,12, r'$\tau_{ff} = 0.34$\,Myr', color='r', ha='center', va='bottom', **text_kwargs)
		ax1.plot([0.0,0.34], [10,10], color='r', linewidth=2, **text_kwargs)
	
		# fit sequences
		fit_ranges = [[-2.0,-1.0],[-0.3,0.05],[1.75,2.25]]
		fit_results = "slopes:\n"
		for fit_range_index in np.arange(len(fit_ranges)):
			fit_lim_lo = fit_ranges[fit_range_index][0]
			fit_lim_hi = fit_ranges[fit_range_index][1]
	
			# fit sequence
			to_fit = [[],[]]
			for row in np.arange(len(data)):
				if (data['time'][row] > fit_lim_lo) and (data['time'][row] < fit_lim_hi) and not (np.isnan(data['T'+temp+'_'+ttype][row])): 
					to_fit[0].append(data['time'][row])
					to_fit[1].append(data['T'+temp+'_'+ttype][row])
#			coeff, cov = np.polyfit(to_fit[0], to_fit[1], 1, cov=True)
#			ax1.plot([-2.0,-1.0],[x*coeff[0]+coeff[1] for x in [-2.0,-1.0]], linewidth=1.0, markersize=0.0, color='black', linestyle='--', label='fit $'+str(fit_lim_lo)+'<t<'+str(fit_lim_hi)+'$', **fit_kwargs)
			
			# fit binned sequence
			bin_fit = [[],[],[],[]]
			for time_bin in np.arange(fit_lim_lo,fit_lim_hi,0.01):
				bin_value = []
				for row in np.arange(len(to_fit[0])):
					if (to_fit[0][row] >= time_bin) and (to_fit[0][row] < time_bin+0.01):
						bin_value.append(to_fit[1][row])
				bin_value.sort()
				if (np.median(bin_value) > 5) and (len(bin_value) > 10):
					bin_fit[0].append(time_bin+0.01)
					bin_fit[1].append(np.median(bin_value))
					bin_fit[2].append(np.median(bin_value) - bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][0])	# lower "error" of 50% include range
					bin_fit[3].append(bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][-1] - np.median(bin_value))	# upper "error" of 50% include range
			coeff_med, cov_med = np.polyfit(bin_fit[0], bin_fit[1], 1, cov=True)
			fit_results += r"$\frac{dT}{dt} = "+str('{:3.1f}'.format(coeff_med[0])).rjust(5)+"$\\,K/Myr ($"+str(fit_lim_lo)+"<t [Myr] <"+str(fit_lim_hi)+"$)\n"
			
			# plot fits; only one legend entry
			if (fit_range_index == 0):
				ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', label='fit bin median', **fit_kwargs)
				ax1.scatter(bin_fit[0], bin_fit[1], label='0.01\,Myr bins median', s=10, marker="x", color='black', **text_kwargs)
				ax1.errorbar(bin_fit[0], bin_fit[1], yerr=[bin_fit[2],bin_fit[3]], fmt='none', linewidth=0.0, elinewidth=0.5, capsize=2.0, ecolor='black', **text_kwargs)
			else:
				ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', **fit_kwargs)
				ax1.scatter(bin_fit[0], bin_fit[1], s=10, marker="x", color='black', **text_kwargs)
				ax1.errorbar(bin_fit[0], bin_fit[1], yerr=[bin_fit[2],bin_fit[3]], fmt='none', linewidth=0.0, elinewidth=0.5, capsize=2.0, ecolor='black', **text_kwargs)
			
		# fit results text box
		ax1.text(0.68, 0.975, fit_results[:-1], transform=ax1.transAxes, fontsize=10, va='top', ha='right', multialignment='left', bbox=props)
		
		# plot data
#		scatter = ax1.scatter(data['time'], data['T'+temp+'_'+ttype], label='T$_{'+ttype+','+temp+'}$', s=1, lw=0, c='dimgrey', **scatter_kwargs)
		
		# contour plot
		contour_list = [[],[]]
		for row in np.arange(len(data)):
			if not (np.isnan(data['time'][row])) and not (np.isnan(data['T'+temp+'_'+ttype][row])):
				contour_list[0].append(data['time'][row])
				contour_list[1].append(data['T'+temp+'_'+ttype][row])

		H, xedges, yedges = np.histogram2d(contour_list[0], contour_list[1], bins=100, range=[[-2.5,2.5],[0,500]])
		extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
#		levels=[25,50,75,125,175]
		levels = [i*np.amax(H) for i in np.arange(0.1,1.0,0.1)]
		contours = ax1.contour(np.transpose(H),extent=extent, levels=levels, label="10% steps of max. density", **scatter_kwargs)

		# set Tkin sensitivity limit
		if (temp == "12") and (ttype == "kin"): 
			ax1.axhline(y=60, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		if (temp == "24") and (ttype == "kin"): 
			ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		if (temp == "45") and (ttype == "kin"): 
			ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		if (temp == "36") and (ttype == "kin"): 
			ax1.axhline(y=300, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)

		# legend
		contour_legend = lines.Line2D([],[], color='b', lw=1, label='10% steps of max. density')
		leg = ax1.legend(fancybox=True, framealpha=0.8, fontsize=10)
		leg.set_zorder(12)
		
		plt.savefig("plots/time_vs_T"+temp+"_"+ttype+".density_only.png", dpi=300)


############################################################

# T?? vs. time (zoom in, density only)
print "plotting T(NH3) vs. time ..."
for ttype in ["kin", "rot"]:
	for temp in ["12", "24", "45", "36"]:
		fit_ranges = [[-2.0,-1.0],[-0.3,0.05],[1.75,2.25]]
		for fit_range_index in np.arange(len(fit_ranges)):
			fit_lim_lo = fit_ranges[fit_range_index][0]
			fit_lim_hi = fit_ranges[fit_range_index][1]
			fit_results = "slope:\n"

			plt.clf()
			plt.cla()
			fig = plt.figure()
			fig.subplots_adjust(top=0.8)
			ax1 = plt.subplot(1,1,1)
			ax2 = ax1.twiny()
			ax1.set_xlabel('time [Myr]')
			ax1.set_ylabel('NH$_3$ temperature T [K]')
			ax1.set_xlim(fit_lim_lo,fit_lim_hi)
			ax1.xaxis.set_minor_locator(MultipleLocator(0.1))
			ax1.set_ylim(0,400)
			if (temp == "12") and (ttype == "kin"): 
				ax1.set_ylim(0,80)
			if (temp == "24") and (ttype == "kin"): 
				ax1.set_ylim(0,280)
			if (temp == "45") and (ttype == "kin"): 
				ax1.set_ylim(0,280)
			if (temp == "36") and (ttype == "kin"): 
				ax1.set_ylim(0,400)
			ax1.yaxis.set_minor_locator(MultipleLocator(10))
			ax2.set_xticks([-2.0, -1.67, -1.52, -1.30, -1.0, 0.0, 1.0, 1.51, 2.0])
			ax2.set_xticklabels(['pericenter\nnear', 'Brick', 'cloud e/f', 'Sgr B2', 'apocenter\neast', 'pericenter\nfar', 'apocenter\nwest', 'Sgr C', 'pericenter\nnear'], rotation=90, ha='center', va='bottom')
			ax2.set_xlim(fit_lim_lo,fit_lim_hi)
#			ax1.bar(-1.7, 500, 0.4, bottom=None, color="grey", alpha=0.5, zorder=3)
#			ax1.text(-1.5, 450, "dust ridge", ha='center', bbox=props, **text_kwargs)
		
			# fit sequence
			to_fit = [[],[]]
			for row in np.arange(len(data)):
				if (data['time'][row] > fit_lim_lo) and (data['time'][row] < fit_lim_hi) and not (np.isnan(data['T'+temp+'_'+ttype][row])): 
					to_fit[0].append(data['time'][row])
					to_fit[1].append(data['T'+temp+'_'+ttype][row])
#			coeff, cov = np.polyfit(to_fit[0], to_fit[1], 1, cov=True)
#			ax1.plot([fit_lim_lo,-1.0],[x*coeff[0]+coeff[1] for x in [fit_lim_lo,-1.0]], linewidth=1.0, markersize=0.0, color='black', linestyle='--', label='fit $'+str(fit_lim_lo)+'<t<'+str(fit_lim_hi)+'$', **fit_kwargs)
			
			# fit binned sequence
			bin_fit = [[],[],[],[]]
			for time_bin in np.arange(fit_lim_lo,fit_lim_hi,0.01):
				bin_value = []
				for row in np.arange(len(to_fit[0])):
					if (to_fit[0][row] >= time_bin) and (to_fit[0][row] < time_bin+0.01):
						bin_value.append(to_fit[1][row])
				bin_value.sort()
				if (np.median(bin_value) > 5) and (len(bin_value) > 10):
					bin_fit[0].append(time_bin+0.01)
					bin_fit[1].append(np.median(bin_value))
					bin_fit[2].append(np.median(bin_value) - bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][0])	# lower "error" of 50% include range
					bin_fit[3].append(bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][-1] - np.median(bin_value))	# upper "error" of 50% include range
			coeff_med, cov_med = np.polyfit(bin_fit[0], bin_fit[1], 1, cov=True)
			fit_results += r"$\frac{dT}{dt} = "+str('{:3.1f}'.format(coeff_med[0])).rjust(5)+"$\\,K/Myr ($"+str(fit_lim_lo)+"<t [Myr] <"+str(fit_lim_hi)+"$)\n"
			
			# plot fits
			ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', label='fit bin median', **fit_kwargs)
			ax1.scatter(bin_fit[0], bin_fit[1], label='0.01\,Myr bins median', s=10, marker="x", color='black', **text_kwargs)
			ax1.errorbar(bin_fit[0], bin_fit[1], yerr=[bin_fit[2],bin_fit[3]], fmt='none', linewidth=0.0, elinewidth=0.5, capsize=2.0, ecolor='black', **text_kwargs)
			
			# fit results text box
			ax1.text(0.68, 0.975, fit_results[:-1], transform=ax1.transAxes, fontsize=10, va='top', ha='right', multialignment='left', bbox=props)
	
			# plot data
#			scatter = ax1.scatter(data['time'], data['T'+temp+'_'+ttype], label='T$_{'+ttype+','+temp+'}$', s=1, lw=0, c='dimgrey', **scatter_kwargs)
			
			# contour plot
			contour_list = [[],[]]
			for row in np.arange(len(data)):
				if not (np.isnan(data['time'][row])) and not (np.isnan(data['T'+temp+'_'+ttype][row])):
					contour_list[0].append(data['time'][row])
					contour_list[1].append(data['T'+temp+'_'+ttype][row])
	
			H, xedges, yedges = np.histogram2d(contour_list[0], contour_list[1], bins=100, range=[[-2.5,2.5],[0,500]])
			extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
#			levels=[25,50,75,125,175]
			levels = [i*np.amax(H) for i in np.arange(0.1,1.0,0.1)]
			contours = ax1.contour(np.transpose(H),extent=extent, levels=levels, label="10% steps of max. density", **scatter_kwargs)

			
			# set Tkin sensitivity limit
			if (temp == "12") and (ttype == "kin"): 
				ax1.axhline(y=60, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			if (temp == "24") and (ttype == "kin"): 
				ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			if (temp == "45") and (ttype == "kin"): 
				ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			if (temp == "36") and (ttype == "kin"): 
				ax1.axhline(y=300, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		
			# legend
			contour_legend = lines.Line2D([],[], color='b', lw=1, label='10% steps of max. density')
			leg = ax1.legend(fancybox=True, framealpha=0.8, fontsize=10)
			leg.set_zorder(12)
			
			plt.savefig("plots/time_vs_T"+temp+"_"+ttype+".density_only.zoom_"+str(fit_lim_lo)+"_"+str(fit_lim_hi)+".png", dpi=300)


############################################################

# T?? vs. Hi-GAL dust
print "plotting T(NH3) vs. T(dust) ..."
for ttype in ["kin", "rot"]:
	for temp in ["12", "24", "45", "36"]:
		plt.clf()
		plt.cla()
		fig = plt.figure()
		fig.subplots_adjust(top=0.8)
		ax1 = plt.subplot(1,1,1)
		ax1.set_xlabel('dust temperature (Hi-GAL) [K]')
		ax1.set_ylabel('NH$_3$ temperature T [K]')
		ax1.set_xlim(15,35)
		ax1.xaxis.set_minor_locator(MultipleLocator(1))
		ax1.set_ylim(0, 300)
		ax1.yaxis.set_minor_locator(MultipleLocator(10))
		if (temp == "12"):
			ax1.set_ylim(0,80)
		if (temp == "24"):
			ax1.set_ylim(0,280)
		if (temp == "45"): 
			ax1.set_ylim(0,280)
		if (temp == "36"): 
			ax1.set_ylim(0,400)

		# set Tkin sensitivity limit
		if (temp == "12") and (ttype == "kin"): 
			ax1.axhline(y=60, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			ax1.set_ylim(0,66)
		if (temp == "24") and (ttype == "kin"): 
			ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			ax1.set_ylim(0,220)
		if (temp == "45") and (ttype == "kin"): 
			ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			ax1.set_ylim(0,220)
		if (temp == "36") and (ttype == "kin"): 
			ax1.axhline(y=300, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			ax1.set_ylim(0,330)

		# plot data
		ax1.scatter(data['Tdust'], data['T'+temp+'_'+ttype], marker='.', linewidth=0.0, s=1.0, color='blue', label='T$_{'+ttype+','+temp+'}$', **scatter_kwargs)
	
		ax1.plot(np.arange(0,100,1),np.arange(0,100,1), linewidth=1.0, markersize=0.0, color='dimgrey', label='T$_{'+ttype+','+temp+'}$ = T$_{dust}$', **fit_kwargs)

		leg = ax1.legend(fancybox=True, framealpha=None, fontsize=10)
		leg.set_zorder(12)
		
		plt.savefig("plots/Tdust_vs_T"+temp+"_"+ttype+".png", dpi=300)


############################################################

# T?? vs. Hi-GAL dust (fit)
print "plotting T(NH3) vs. T(dust) ..."
for ttype in ["kin", "rot"]:
	for temp in ["12", "24", "45", "36"]:
		fit_results = "slope:\n"
		plt.clf()
		plt.cla()
		fig = plt.figure()
		ax1 = plt.subplot(1,1,1)
		ax1.set_xlabel('dust temperature (Hi-GAL) [K]')
		ax1.set_ylabel('NH$_3$ temperature T [K]')
		ax1.set_xlim(15,35)
		ax1.xaxis.set_minor_locator(MultipleLocator(1))
		ax1.set_ylim(0, 300)
		ax1.yaxis.set_minor_locator(MultipleLocator(10))
		if (temp == "12"):
			ax1.set_ylim(0,80)
		if (temp == "24"):
			ax1.set_ylim(0,280)
		if (temp == "45"): 
			ax1.set_ylim(0,280)
		if (temp == "36"): 
			ax1.set_ylim(0,400)

		if (temp == "12") and (ttype == "kin"): 
			ax1.axhline(y=60, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			ax1.set_ylim(0,66)
		if (temp == "24") and (ttype == "kin"): 
			ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			ax1.set_ylim(0,220)
		if (temp == "45") and (ttype == "kin"): 
			ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			ax1.set_ylim(0,220)
		if (temp == "36") and (ttype == "kin"): 
			ax1.axhline(y=300, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			ax1.set_ylim(0,330)

		ax1.scatter(data['Tdust'], data['T'+temp+'_'+ttype], marker='.', linewidth=0.0, s=1.0, color='blue', label='T$_{'+ttype+','+temp+'}$', **scatter_kwargs)
	
		ax1.plot(np.arange(0,100,1),np.arange(0,100,1), linewidth=1.0, markersize=0.0, color='grey', label='T$_{'+ttype+','+temp+'}$ = T$_{dust}$', **fit_kwargs)

		to_fit = [[],[]]
		for row in np.arange(len(data)):
			if not (np.isnan(data['Tdust'][row])) and not (np.isnan(data['T'+temp+'_'+ttype][row])): 
				to_fit[0].append(data['Tdust'][row])
				to_fit[1].append(data['T'+temp+'_'+ttype][row])

		bin_fit = [[],[],[],[]]
		for T_bin in np.arange(18,30,1):
			bin_value = []
			for row in np.arange(len(to_fit[0])):
				if (to_fit[0][row] >= T_bin) and (to_fit[0][row] < T_bin+1):
					bin_value.append(to_fit[1][row])
			bin_value.sort()
			if (np.median(bin_value) > 5) and (len(bin_value) > 10):
				bin_fit[0].append(T_bin)
				bin_fit[1].append(np.median(bin_value))
				bin_fit[2].append(np.median(bin_value) - bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][0])
				bin_fit[3].append(bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][-1] - np.median(bin_value))
		coeff_med, cov_med = np.polyfit(bin_fit[0], bin_fit[1], 1, cov=True)
		fit_results += r"$\frac{dT}{dT} = "+str('{:3.1f}'.format(coeff_med[0])).rjust(5)+"$\n"
		
		ax1.plot([18,30],[x*coeff_med[0]+coeff_med[1] for x in [18,30]], linewidth=1.0, markersize=0.0, linestyle='--', color='black', label='fit bin median', **fit_kwargs)
		ax1.scatter(bin_fit[0], bin_fit[1], label='1.0\,K bins median', s=10, marker="x", color='black', **text_kwargs)
		ax1.errorbar(bin_fit[0], bin_fit[1], yerr=[bin_fit[2],bin_fit[3]], fmt='none', linewidth=0.0, elinewidth=0.5, capsize=2.0, ecolor='black', **text_kwargs)
			
		ax1.text(0.98, 0.1, fit_results[:-1], transform=ax1.transAxes, fontsize=10, va='top', ha='right', multialignment='left', bbox=props)

		leg = ax1.legend(fancybox=True, framealpha=None, fontsize=10)
		leg.set_zorder(12)
		
		plt.savefig("plots/Tdust_vs_T"+temp+"_"+ttype+".fit.png", dpi=300)


############################################################

# T?? vs. Hi-GAL dust (density)
print "plotting T(NH3) vs. T(dust) ..."
for ttype in ["kin", "rot"]:
	for temp in ["12", "24", "45", "36"]:
		plt.clf()
		plt.cla()
		fig = plt.figure()
		fig.subplots_adjust(top=0.8)
		ax1 = plt.subplot(1,1,1)
		ax1.set_xlabel('dust temperature (Hi-GAL) [K]')
		ax1.set_ylabel('NH$_3$ temperature T [K]')
		ax1.set_xlim(15,35)
		ax1.xaxis.set_minor_locator(MultipleLocator(1))
		ax1.set_ylim(0, 250)
		ax1.yaxis.set_minor_locator(MultipleLocator(10))
		if (temp == "12"):
			ax1.set_ylim(0,80)
		if (temp == "24"):
			ax1.set_ylim(0,280)
		if (temp == "45"): 
			ax1.set_ylim(0,280)
		if (temp == "36"): 
			ax1.set_ylim(0,400)

		if (temp == "12") and (ttype == "kin"): 
			ax1.axhline(y=60, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			ax1.set_ylim(0,66)
		if (temp == "24") and (ttype == "kin"): 
			ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			ax1.set_ylim(0,220)
		if (temp == "45") and (ttype == "kin"): 
			ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			ax1.set_ylim(0,220)
		if (temp == "36") and (ttype == "kin"): 
			ax1.axhline(y=300, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
			ax1.set_ylim(0,330)
	
		# plot data
		ax1.scatter(data['Tdust'], data['T'+temp+'_'+ttype], marker='.', linewidth=0.0, s=1.0, c='dimgrey', label='T$_{'+ttype+','+temp+'}$', **scatter_kwargs)
	
		ax1.plot(np.arange(0,100,1),np.arange(0,100,1), linewidth=1.0, markersize=0.0, color='grey', label='T$_{'+ttype+','+temp+'}$ = T$_{dust}$', **fit_kwargs)

		# contour
		contour_list = [[],[]]
		for row in np.arange(len(data)):
			if not (np.isnan(data['Tdust'][row])) and not (np.isnan(data['T'+temp+'_'+ttype][row])):
				contour_list[0].append(data['Tdust'][row])
				contour_list[1].append(data['T'+temp+'_'+ttype][row])
	
		H, xedges, yedges = np.histogram2d(contour_list[0], contour_list[1], bins=[50,150], range=[[15,35],[0,300]])
		extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
		levels = [i*np.amax(H) for i in np.arange(0.1,1.0,0.1)]
		contours = ax1.contour(np.transpose(H),extent=extent, levels=levels, label="10% steps of max. density", **scatter_kwargs)

		# legend
		leg = ax1.legend(fancybox=True, framealpha=None, fontsize=10)
		leg.set_zorder(12)
		
		plt.savefig("plots/Tdust_vs_T"+temp+"_"+ttype+".density.png", dpi=300)


############################################################

# T?? vs. FWHM (density)
print "plotting T(NH3) vs. FWHM ..."
for ttype in ["kin"]: #["kin", "rot"]:
	for temp in ["12", "24", "45", "36"]:
		plt.clf()
		plt.cla()
		fig = plt.figure()
		ax1 = plt.subplot(1,1,1)
		ax1.set_xlabel('FWHM [km/s]')
		ax1.set_ylabel('NH$_3$ temperature T [K]')
		ax1.set_xlim(0,50)
		ax1.xaxis.set_minor_locator(MultipleLocator(2))
		ax1.set_ylim(0, 400)
		if (temp == "12"):
			ax1.set_ylim(0,80)
		if (temp == "24"):
			ax1.set_ylim(0,280)
		if (temp == "45"):
			ax1.set_ylim(0,280)
		if (temp == "36"):
			ax1.set_ylim(0,400)
		ax1.yaxis.set_minor_locator(MultipleLocator(10))
	
		# contour plot
		contour_list = [[],[]]
		for row in np.arange(len(data)):
			if (data['width'+temp[0]][row] > 3.5) and (data['T'+temp+'_'+ttype][row] > 10.0):
				contour_list[0].append(data['width'+temp[0]][row])
				contour_list[1].append(data['T'+temp+'_'+ttype][row])

		# plot data
		scatter = ax1.scatter(contour_list[0], contour_list[1], label='T$_{'+ttype+','+temp+'}$', s=1, lw=0, c='dimgrey', **scatter_kwargs)

		H, xedges, yedges = np.histogram2d(contour_list[0], contour_list[1], bins=100, range=[[0,50],[0,400]])
		extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
#		levels=[25,50,75,125,175]
		levels = [i*np.amax(H) for i in np.arange(0.1,1.0,0.1)]
		contours = ax1.contour(np.transpose(H),extent=extent, levels=levels, label="10% steps of max. density", **scatter_kwargs)

# fit without binning is a bad idea		
#		# fit
#		coeff, cov = np.polyfit(contour_list[0], contour_list[1], 1, cov=True)
#		ax1.plot([0,50],[x*coeff[0]+coeff[1] for x in [0,50]], linewidth=1.0, markersize=0.0, color='black', linestyle='--', label='linear fit', **fit_kwargs)
#		
#		fit_results = "slope:\n"
#		fit_results += r"$\frac{dT}{dFWHM} = "+str('{:3.1f}'.format(coeff[0])).rjust(5)+"$\\,K/(km/s)\n"
#		ax1.text(0.68, 0.975, fit_results[:-1], transform=ax1.transAxes, fontsize=10, va='top', ha='right', multialignment='left', bbox=props)

		# fit
		fit_results = ""
		to_fit = [[],[]]
		for row in np.arange(len(data)):
			if not (np.isnan(data['width'+temp[0]][row])) and not (np.isnan(data['T'+temp+'_'+ttype][row])): 
				to_fit[0].append(data['width'+temp[0]][row])
				to_fit[1].append(data['T'+temp+'_'+ttype][row])

		bin_fit = [[],[],[],[]]
		for T_bin in np.arange(3,40,2):
			bin_value = []
			for row in np.arange(len(to_fit[0])):
				if (to_fit[0][row] >= T_bin) and (to_fit[0][row] < T_bin+1):
					bin_value.append(to_fit[1][row])
			bin_value.sort()
			if (np.median(bin_value) > 5) and (len(bin_value) > 10):
				bin_fit[0].append(T_bin)
				bin_fit[1].append(np.median(bin_value))
				bin_fit[2].append(np.median(bin_value) - bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][0])
				bin_fit[3].append(bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][-1] - np.median(bin_value))
		coeff_med, cov_med = np.polyfit(bin_fit[0], bin_fit[1], 1, cov=True)
		fit_results += r"$\frac{dT}{dFWHM} = "+str('{:3.1f}'.format(coeff_med[0])).rjust(5)+"$\n"
		
		ax1.plot([3,40],[x*coeff_med[0]+coeff_med[1] for x in [3,40]], linewidth=1.0, markersize=0.0, linestyle='--', color='black', label='fit bin median', **fit_kwargs)
		ax1.scatter(bin_fit[0], bin_fit[1], label='2.0\,km/s bins median', s=10, marker="x", color='black', **text_kwargs)
#		ax1.errorbar(bin_fit[0], bin_fit[1], yerr=[bin_fit[2],bin_fit[3]], fmt='none', linewidth=0.0, elinewidth=0.5, capsize=2.0, ecolor='black', **text_kwargs)
			
		ax1.text(0.975, 0.025, fit_results[:-1], transform=ax1.transAxes, fontsize=10, va='bottom', ha='right', multialignment='left', bbox=props)


		# set Tkin sensitivity limit
		if (temp == "12") and (ttype == "kin"): 
			ax1.axhline(y=60, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		if (temp == "24") and (ttype == "kin"): 
			ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		if (temp == "45") and (ttype == "kin"): 
			ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		if (temp == "36") and (ttype == "kin"): 
			ax1.axhline(y=300, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)

		# legend
		contour_legend = lines.Line2D([],[], color='b', lw=1, label='10% steps of max. density')
		leg = ax1.legend(fancybox=True, framealpha=0.8, fontsize=10)
		leg.set_zorder(12)
		
		plt.savefig("plots/FWHM("+temp[0]+","+temp[0]+")_vs_T"+temp+"_"+ttype+".density.fit.png", dpi=300)


############################################################

# T?? vs. galactic longitude
print "plotting T(NH3) vs. glon ..."
for ttype in ["kin", "rot"]:
	for temp in ["12", "24", "45", "36"]:
		plt.clf()
		plt.cla()
		fig = plt.figure()
		fig.subplots_adjust(top=0.8)
		ax1 = plt.subplot(1,1,1)
		ax2 = ax1.twiny()
		ax1.set_xlabel('galactic longitude [$^\circ$]')
		ax1.set_ylabel('NH$_3$ temperature T [K]')
		ax1.set_xlim(0.8,-0.6)
		ax1.xaxis.set_minor_locator(MultipleLocator(0.1))
		ax1.set_ylim(0, 400)
		if (temp == "12"):
			ax1.set_ylim(0,80)
		if (temp == "24"):
			ax1.set_ylim(0,280)
		if (temp == "45"):
			ax1.set_ylim(0,280)
		if (temp == "36"):
			ax1.set_ylim(0,400)
		ax1.yaxis.set_minor_locator(MultipleLocator(10))
		ax2.set_xlim(1.0,-1.0)
		ax2.set_xticks([0.75, 0.612, 0.053, -0.083, -0.121, -0.788, -0.85])
		ax2.set_xticklabels(['apocenter\neast', 'Sgr B2', 'pericenter\nnear', 'far', 'near', 'Sgr C', 'apocenter\nwest'], rotation=90, ha='center', va='bottom')
		ax1.bar(0.2, 500, 0.4, bottom=None, color="grey", alpha=0.5, zorder=3)
		ax1.text(-1.5, 450, "dust ridge", ha='center', bbox=props, **text_kwargs)
	
		# plot data
		scatter = ax1.scatter(data['glon'], data['T'+temp+'_'+ttype], label='T$_{'+ttype+','+temp+'}$', s=1, lw=0, c=data['v_mean'], cmap=cmap, vmin=vmin, vmax=vmax, **scatter_kwargs)
		
		# set Tkin sensitivity limit
		if (temp == "12") and (ttype == "kin"): 
			ax1.axhline(y=60, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		if (temp == "24") and (ttype == "kin"): 
			ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		if (temp == "45") and (ttype == "kin"): 
			ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		if (temp == "36") and (ttype == "kin"): 
			ax1.axhline(y=300, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)

		# color bar
		cbaxes = fig.add_axes([0.9, 0.1, 0.02, 0.7]) 
		cbar = plt.colorbar(scatter, cax=cbaxes)
		cbar.set_label("v$_{los}$ [km/s]")
		
		leg = ax1.legend(fancybox=True, framealpha=0.8, fontsize=10)
		leg.set_zorder(12)
		
		plt.savefig("plots/glon_vs_T"+temp+"_"+ttype+".png", dpi=300)


############################################################

# v_lsr vs. time
print "plotting v_lsr vs. time ..."
for i in np.arange(1,7):
	line = str(i)
	plt.clf()
	plt.cla()
	fig = plt.figure()
	fig.subplots_adjust(top=0.8)
	ax1 = plt.subplot(1,1,1)
	ax2 = ax1.twiny()
	ax1.set_xlabel('time [Myr]')
	ax1.set_ylabel('$v_{lsr}$ [km/s]')
	ax1.set_xlim(-2.5,2.5)
	ax1.xaxis.set_minor_locator(MultipleLocator(0.2))
	ax1.set_ylim(-150, 150)
	ax1.yaxis.set_minor_locator(MultipleLocator(10))
	ax2.set_xlim(-2.5,2.5)
	ax2.set_xticks([-2.0, -1.67, -1.52, -1.30, -1.0, 0.0, 1.0, 1.51, 2.0])
	ax2.set_xticklabels(['pericenter\nnear', 'Brick', 'cloud e/f', 'Sgr B2', 'apocenter\neast', 'pericenter\nfar', 'apocenter\nwest', 'Sgr C', 'pericenter\nnear'], rotation=90, ha='center', va='bottom')
#	ax1.bar(-1.7, 50, 0.4, bottom=None, color="grey", alpha=0.5, zorder=3)
#	ax1.text(-1.5, 45, "dust ridge", ha='center', bbox=props, **text_kwargs)
#	ax1.text(0.17,-130, r'$\tau_{ff} = 0.34$\,Myr', color='r', ha='center', va='bottom', **text_kwargs)
#	ax1.plot([0.0,0.34], [-140,-140], color='r', linewidth=2, **text_kwargs)

	# plot data
	scatter = ax1.scatter(data['time'], data['v'+line], label='NH$_3$ ('+line+','+line+')', s=1, lw=0, **scatter_kwargs)
	
	leg = ax1.legend(fancybox=True, framealpha=0.8, fontsize=10, loc=4)
	leg.set_zorder(12)
	
	plt.savefig("plots/time_vs_vlsr("+line+","+line+").png", dpi=300)


############################################################

# v_lsr vs. time
print "plotting v_mean vs. time ..."
plt.clf()
plt.cla()
fig = plt.figure()
fig.subplots_adjust(top=0.8)
ax1 = plt.subplot(1,1,1)
ax2 = ax1.twiny()
ax1.set_xlabel('time [Myr]')
ax1.set_ylabel('$v_{lsr}$ [km/s]')
ax1.set_xlim(-2.5,2.5)
ax1.xaxis.set_minor_locator(MultipleLocator(0.2))
ax1.set_ylim(-150, 150)
ax1.yaxis.set_minor_locator(MultipleLocator(10))
ax2.set_xlim(-2.5,2.5)
ax2.set_xticks([-2.0, -1.67, -1.52, -1.30, -1.0, 0.0, 1.0, 1.51, 2.0])
ax2.set_xticklabels(['pericenter\nnear', 'Brick', 'cloud e/f', 'Sgr B2', 'apocenter\neast', 'pericenter\nfar', 'apocenter\nwest', 'Sgr C', 'pericenter\nnear'], rotation=90, ha='center', va='bottom')
#ax1.bar(-1.7, 50, 0.4, bottom=None, color="grey", alpha=0.5, zorder=3)
#ax1.text(-1.5, 45, "dust ridge", ha='center', bbox=props, **text_kwargs)
#ax1.text(0.17,-130, r'$\tau_{ff} = 0.34$\,Myr', color='r', ha='center', va='bottom', **text_kwargs)
#ax1.plot([0.0,0.34], [-140,-140], color='r', linewidth=2, **text_kwargs)

# plot data
scatter = ax1.scatter(data['time'], data['v_mean'], label='NH$_3$', s=1, lw=0, color='dimgrey', **scatter_kwargs)

# plot streams with color
stream1 = ax1.plot(model['t'][np.where(np.logical_and(model['t']>1.0, model['t']<2.5))], 
		model['vlos'][np.where(np.logical_and(model['t']>1.0, model['t']<2.5))], 
		label='stream 1', lw=1.5, color='r', **fit_kwargs)
stream2 = ax1.plot(model['t'][np.where(np.logical_and(model['t']>-2.5, model['t']<-1.0))], 
		model['vlos'][np.where(np.logical_and(model['t']>-2.5, model['t']<-1.0))], 
		label='stream 2', lw=1.5, color='b', **fit_kwargs)
stream3 = ax1.plot(model['t'][np.where(np.logical_and(model['t']>-1.0, model['t']<0.0))], 
		model['vlos'][np.where(np.logical_and(model['t']>-1.0, model['t']<0.0))], 
		label='stream 3', lw=1.5, color='gold', **fit_kwargs)
stream4 = ax1.plot(model['t'][np.where(np.logical_and(model['t']>0.0, model['t']<1.0))], 
		model['vlos'][np.where(np.logical_and(model['t']>0.0, model['t']<1.0))], 
		label='stream 4', lw=1.5, color='cyan', **fit_kwargs)
	
leg = ax1.legend(fancybox=True, framealpha=0.8, fontsize=10, loc=4)
leg.set_zorder(12)
	
plt.savefig("plots/time_vs_vmean.png", dpi=300)


############################################################

# line width vs. time
print "plotting FWHM vs. time ..."
for i in np.arange(1,7):
	line = str(i)
	plt.clf()
	plt.cla()
	fig = plt.figure()
	fig.subplots_adjust(top=0.8)
	ax1 = plt.subplot(1,1,1)
	ax2 = ax1.twiny()
	ax1.set_xlabel('time [Myr]')
	ax1.set_ylabel('NH$_3$ FWHM [km/s]')
	ax1.set_xlim(-2.5,2.5)
	ax1.xaxis.set_minor_locator(MultipleLocator(0.2))
	ax1.set_ylim(0, 50)
	ax1.yaxis.set_minor_locator(MultipleLocator(2))
	ax2.set_xlim(-2.5,2.5)
	ax2.set_xticks([-2.0, -1.67, -1.52, -1.30, -1.0, 0.0, 1.0, 1.51, 2.0])
	ax2.set_xticklabels(['pericenter\nnear', 'Brick', 'cloud e/f', 'Sgr B2', 'apocenter\neast', 'pericenter\nfar', 'apocenter\nwest', 'Sgr C', 'pericenter\nnear'], rotation=90, ha='center', va='bottom')
	ax1.bar(-1.7, 50, 0.4, bottom=None, color="grey", alpha=0.5, zorder=3)
	ax1.text(-1.5, 45, "dust ridge", ha='center', bbox=props, **text_kwargs)
	ax1.text(0.17,46, r'$\tau_{ff} = 0.34$\,Myr', color='r', ha='center', va='bottom', **text_kwargs)
	ax1.plot([0.0,0.34], [45,45], color='r', linewidth=2, **text_kwargs)

	# fit sequences
	fit_ranges = [[-2.0,-1.0],[-0.3,0.05],[1.75,2.25]]
	fit_results = "slopes:\n"
	for fit_range_index in np.arange(len(fit_ranges)):
		fit_lim_lo = fit_ranges[fit_range_index][0]
		fit_lim_hi = fit_ranges[fit_range_index][1]

		# fit sequence
		to_fit = [[],[]]
		for row in np.arange(len(data)):
			if (data['time'][row] > fit_lim_lo) and (data['time'][row] < fit_lim_hi) and not (np.isnan(data['width'+line][row])): 
				to_fit[0].append(data['time'][row])
				to_fit[1].append(data['width'+line][row])
#		coeff, cov = np.polyfit(to_fit[0], to_fit[1], 1, cov=True)
#		ax1.plot([-2.0,-1.0],[x*coeff[0]+coeff[1] for x in [-2.0,-1.0]], linewidth=1.0, markersize=0.0, color='black', linestyle='--', label='fit $'+str(fit_lim_lo)+'<t<'+str(fit_lim_hi)+'$', **fit_kwargs)
		
		# fit binned sequence
		bin_fit = [[],[]]
		for time_bin in np.arange(fit_lim_lo,fit_lim_hi,0.01):
			bin_value = []
			for row in np.arange(len(to_fit[0])):
				if (to_fit[0][row] >= time_bin) and (to_fit[0][row] < time_bin+0.01):
					bin_value.append([to_fit[1][row]])
			if (np.median(bin_value) > 5) and (len(bin_value) > 10):
				bin_fit[0].append(time_bin+0.01)
				bin_fit[1].append(np.median(bin_value))
		coeff_med, cov_med = np.polyfit(bin_fit[0], bin_fit[1], 1, cov=True)
		
		# plot fits; only one legend entry
		if (fit_range_index == 0):
			ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', label='fit bin median', **fit_kwargs)
			ax1.scatter(bin_fit[0], bin_fit[1], label='0.01\,Myr bins median', s=10, marker="x", color='black', **text_kwargs)
		else:
			ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', **fit_kwargs)
			ax1.scatter(bin_fit[0], bin_fit[1], s=10, marker="x", color='black', **text_kwargs)

	# plot data
	scatter = ax1.scatter(data['time'], data['width'+line], label='('+line+','+line+')', s=1, lw=0, c=data['v_mean'], cmap=cmap, vmin=vmin, vmax=vmax, **scatter_kwargs)
	
	# color bar	
	cbaxes = fig.add_axes([0.9, 0.1, 0.02, 0.7]) 
	cbar = plt.colorbar(scatter, cax=cbaxes)
	cbar.set_label("v$_{los}$ [km/s]")
		
	leg = ax1.legend(fancybox=True, framealpha=0.8, fontsize=10)
	leg.set_zorder(12)
	
	plt.savefig("plots/time_vs_FWHM("+line+","+line+").png", dpi=300)


############################################################

# line width vs. time (density)
print "plotting FWHM vs. time ..."
for i in np.arange(1,7):
	line = str(i)
	plt.clf()
	plt.cla()
	fig = plt.figure()
	fig.subplots_adjust(top=0.8)
	ax1 = plt.subplot(1,1,1)
	ax2 = ax1.twiny()
	ax1.set_xlabel('time [Myr]')
	ax1.set_ylabel('NH$_3$ FWHM [km/s]')
	ax1.set_xlim(-2.5,2.5)
	ax1.xaxis.set_minor_locator(MultipleLocator(0.2))
	ax1.set_ylim(0, 50)
	ax1.yaxis.set_minor_locator(MultipleLocator(2))
	ax2.set_xlim(-2.5,2.5)
	ax2.set_xticks([-2.0, -1.67, -1.52, -1.30, -1.0, 0.0, 1.0, 1.51, 2.0])
	ax2.set_xticklabels(['pericenter\nnear', 'Brick', 'cloud e/f', 'Sgr B2', 'apocenter\neast', 'pericenter\nfar', 'apocenter\nwest', 'Sgr C', 'pericenter\nnear'], rotation=90, ha='center', va='bottom')
	ax1.bar(-1.7, 50, 0.4, bottom=None, color="grey", alpha=0.5, zorder=3)
	ax1.text(-1.5, 45, "dust ridge", ha='center', bbox=props, **text_kwargs)
	ax1.text(0.17,46, r'$\tau_{ff} = 0.34$\,Myr', color='r', ha='center', va='bottom', **text_kwargs)
	ax1.plot([0.0,0.34], [45,45], color='r', linewidth=2, **text_kwargs)

	# fit sequences
	fit_ranges = [[-2.0,-1.0],[-0.3,0.05],[1.75,2.25]]
	fit_results = "slopes:\n"
	for fit_range_index in np.arange(len(fit_ranges)):
		fit_lim_lo = fit_ranges[fit_range_index][0]
		fit_lim_hi = fit_ranges[fit_range_index][1]

		# fit sequence
		to_fit = [[],[]]
		for row in np.arange(len(data)):
			if (data['time'][row] > fit_lim_lo) and (data['time'][row] < fit_lim_hi) and not (np.isnan(data['width'+line][row])): 
				to_fit[0].append(data['time'][row])
				to_fit[1].append(data['width'+line][row])
#		coeff, cov = np.polyfit(to_fit[0], to_fit[1], 1, cov=True)
#		ax1.plot([-2.0,-1.0],[x*coeff[0]+coeff[1] for x in [-2.0,-1.0]], linewidth=1.0, markersize=0.0, color='black', linestyle='--', label='fit $'+str(fit_lim_lo)+'<t<'+str(fit_lim_hi)+'$', **fit_kwargs)
		
		# fit binned sequence
		bin_fit = [[],[]]
		for time_bin in np.arange(fit_lim_lo,fit_lim_hi,0.01):
			bin_value = []
			for row in np.arange(len(to_fit[0])):
				if (to_fit[0][row] >= time_bin) and (to_fit[0][row] < time_bin+0.01):
					bin_value.append([to_fit[1][row]])
			if (np.median(bin_value) > 5) and (len(bin_value) > 10):
				bin_fit[0].append(time_bin+0.01)
				bin_fit[1].append(np.median(bin_value))
		coeff_med, cov_med = np.polyfit(bin_fit[0], bin_fit[1], 1, cov=True)
		
		# plot fits; only one legend entry
		if (fit_range_index == 0):
			ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', label='fit bin median', **fit_kwargs)
			ax1.scatter(bin_fit[0], bin_fit[1], label='0.01\,Myr bins median', s=10, marker="x", color='black', **text_kwargs)
		else:
			ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', **fit_kwargs)
			ax1.scatter(bin_fit[0], bin_fit[1], s=10, marker="x", color='black', **text_kwargs)

	# contour plot
	contour_list = [[],[]]
	for row in np.arange(len(data)):
		if not (np.isnan(data['time'][row])) and not (np.isnan(data['width'+line][row])):
			contour_list[0].append(data['time'][row])
			contour_list[1].append(data['width'+line][row])
	
	H, xedges, yedges = np.histogram2d(contour_list[0], contour_list[1], bins=100, range=[[-2.5,2.5],[10,45]])
	extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
	#levels=[25,50,75,125,175]
	levels = [i*np.amax(H) for i in np.arange(0.1,1.0,0.1)]
	contours = ax1.contour(np.transpose(H),extent=extent, levels=levels, label="10% steps of max. density", **scatter_kwargs)


	# plot data
	scatter = ax1.scatter(data['time'], data['width'+line], label='('+line+','+line+')', s=1, lw=0, c='dimgrey', **scatter_kwargs)
	
	# legend
	contour_legend = lines.Line2D([],[], color='b', lw=1, label='10% steps of max. density')		
	leg = ax1.legend(fancybox=True, framealpha=0.8, fontsize=10)
	leg.set_zorder(12)
	
	plt.savefig("plots/time_vs_FWHM("+line+","+line+").density.png", dpi=300)


############################################################

# column density vs. opacity (density)
print "plotting column density vs. opacity ..."
for j in np.arange(1,7):
	J = str(j)
	plt.clf()
	plt.cla()
	fig = plt.figure()
	ax1 = plt.subplot(1,1,1)
	ax1.set_title('NH$_3$ ('+J+','+J+')')
	ax1.set_xlabel('$\log_{10}$ column density [cm$^{-2}$]')
	ax1.set_ylabel('opacity')
	ax1.set_xlim(12.0,16.0)
	ax1.xaxis.set_minor_locator(MultipleLocator(0.5))
	ax1.set_ylim(0, 12.5)
	ax1.yaxis.set_minor_locator(MultipleLocator(1))

	# plot data
	scatter = ax1.scatter(np.log10(data['N'+J]), data['tau'+J], label='J = '+J, s=1, lw=0, c='dimgrey', **scatter_kwargs)
	
	# countours
	contour_list = [[],[]]
	for row in np.arange(len(data)):
		if not (np.isnan(data['N'+J][row])) and not (np.isnan(data['tau'+J][row])):
			if (data['tau'+J][row] > 11):
				plt.plot(np.log10(data['N'+J][row]), 12, marker='^', color='r', **fit_kwargs)
			else:
				contour_list[0].append(np.log10(data['N'+J][row]))
				contour_list[1].append(data['tau'+J][row])

	H, xedges, yedges = np.histogram2d(contour_list[0], contour_list[1], bins=100, range=[[12.0,16.0],[0.1,11]])
 	extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
	levels = [5,15,25,35,45,55]
	contours = ax1.contour(np.transpose(H),extent=extent, levels=levels, label="10% steps of max. density", **scatter_kwargs)

	plt.savefig("plots/col_dens_vs_opacity.("+J+","+J+").density.png", dpi=300)


############################################################

# Hi-GAL dust temperature vs. time
print "plotting T(dust) vs. time ..."
plt.clf()
plt.cla()
fig = plt.figure()
fig.subplots_adjust(top=0.8)
ax1 = plt.subplot(1,1,1)
ax2 = ax1.twiny()
ax1.set_xlabel('time [Myr]')
ax1.set_ylabel('Hi-GAL dust temperature T [K]')
ax1.set_xlim(-2.5,2.5)
ax1.xaxis.set_minor_locator(MultipleLocator(0.2))
ax1.set_ylim(10,45)
ax1.yaxis.set_minor_locator(MultipleLocator(5))
ax2.set_xlim(-2.5,2.5)
ax2.set_xticks([-2.0, -1.67, -1.52, -1.30, -1.0, 0.0, 1.0, 1.51, 2.0])
ax2.set_xticklabels(['pericenter\nnear', 'Brick', 'cloud e/f', 'Sgr B2', 'apocenter\neast', 'pericenter\nfar', 'apocenter\nwest', 'Sgr C', 'pericenter\nnear'], rotation=90, ha='center', va='bottom')
ax1.bar(-1.7, 500, 0.4, bottom=None, color="grey", alpha=0.5, zorder=3)
ax1.text(-1.5, 450, "dust ridge", ha='center', bbox=props, **text_kwargs)
ax1.text(0.17,13, r'$\tau_{ff} = 0.34$\,Myr', color='r', ha='center', va='bottom', **text_kwargs)
ax1.plot([0.0,0.34], [12,12], color='r', linewidth=2, **text_kwargs)

	
# fit sequences
fit_ranges = [[-2.0,-1.0],[-0.3,0.05],[1.75,2.25]]
fit_results = "slopes:\n"
for fit_range_index in np.arange(len(fit_ranges)):
	fit_lim_lo = fit_ranges[fit_range_index][0]
	fit_lim_hi = fit_ranges[fit_range_index][1]
	
	# fit sequence
	to_fit = [[],[]]
	for row in np.arange(len(data)):
		if (data['time'][row] > fit_lim_lo) and (data['time'][row] < fit_lim_hi) and not (np.isnan(data['Tdust'][row])): 
			to_fit[0].append(data['time'][row])
			to_fit[1].append(data['Tdust'][row])
#	coeff, cov = np.polyfit(to_fit[0], to_fit[1], 1, cov=True)
#	ax1.plot([-2.0,-1.0],[x*coeff[0]+coeff[1] for x in [-2.0,-1.0]], linewidth=1.0, markersize=0.0, color='black', linestyle='--', label='fit $'+str(fit_lim_lo)+'<t<'+str(fit_lim_hi)+'$', **fit_kwargs)
			
	# fit binned sequence
	bin_fit = [[],[],[],[]]
	for time_bin in np.arange(fit_lim_lo,fit_lim_hi,0.01):
		bin_value = []
		for row in np.arange(len(to_fit[0])):
			if (to_fit[0][row] >= time_bin) and (to_fit[0][row] < time_bin+0.01):
				bin_value.append(to_fit[1][row])
		bin_value.sort()
		if (np.median(bin_value) > 5) and (len(bin_value) > 10):
			bin_fit[0].append(time_bin+0.01)
			bin_fit[1].append(np.median(bin_value))
			bin_fit[2].append(np.median(bin_value) - bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][0])	# lower "error" of 50% include range
			bin_fit[3].append(bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][-1] - np.median(bin_value))	# upper "error" of 50% include range
	coeff_med, cov_med = np.polyfit(bin_fit[0], bin_fit[1], 1, cov=True)
	fit_results += r"$\frac{dT}{dt} = "+str('{:3.1f}'.format(coeff_med[0])).rjust(5)+"$\\,K/Myr ($"+str(fit_lim_lo)+"<t [Myr] <"+str(fit_lim_hi)+"$)\n"
			
	# plot fits; only one legend entry
	if (fit_range_index == 0):
		ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', label='fit bin median', **fit_kwargs)
		ax1.scatter(bin_fit[0], bin_fit[1], label='0.01\,Myr bins median', s=10, marker="x", color='black', **text_kwargs)
#		ax1.errorbar(bin_fit[0], bin_fit[1], yerr=zip(*[bin_fit[2],bin_fit[3]]), fmt='none', linewidth=0.0, elinewidth=0.5, ecolor='black', **limit_kwargs)
	else:
		ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', **fit_kwargs)
		ax1.scatter(bin_fit[0], bin_fit[1], s=10, marker="x", color='black', **text_kwargs)
#		ax1.errorbar(bin_fit[0], bin_fit[1], yerr=zip(*[bin_fit[2],bin_fit[3]]), fmt='none', linewidth=0.0, elinewidth=0.5, ecolor='black', **limit_kwargs)
		
# fit results text box
ax1.text(0.68, 0.975, fit_results[:-1], transform=ax1.transAxes, fontsize=10, va='top', ha='right', multialignment='left', bbox=props)
		
# plot data
scatter = ax1.scatter(data['time'], data['Tdust'], label='T$_{dust}$', s=1, lw=0, c='dimgrey', **scatter_kwargs)
		
# contour plot
contour_list = [[],[]]
for row in np.arange(len(data)):
	if not (np.isnan(data['time'][row])) and not (np.isnan(data['Tdust'][row])):
		contour_list[0].append(data['time'][row])
		contour_list[1].append(data['Tdust'][row])

H, xedges, yedges = np.histogram2d(contour_list[0], contour_list[1], bins=100, range=[[-2.5,2.5],[10,45]])
extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
#levels=[25,50,75,125,175]
levels = [i*np.amax(H) for i in np.arange(0.1,1.0,0.1)]
contours = ax1.contour(np.transpose(H),extent=extent, levels=levels, label="10% steps of max. density", **scatter_kwargs)

# legend
contour_legend = lines.Line2D([],[], color='b', lw=1, label='10% steps of max. density')
leg = ax1.legend(fancybox=True, framealpha=0.8, fontsize=10)
leg.set_zorder(12)
		
plt.savefig("plots/time_vs_Tdust.png", dpi=300)


############################################################

# T?? vs. phase (density, single fit)
print "plotting T(NH3) vs. phase ..."
for ttype in ["kin", "rot"]:
	for temp in np.arange(len(temp_info[0])):
		plt.clf()
		plt.cla()
		fig = plt.figure()
		ax1 = plt.subplot(1,1,1)
		ax1.set_xlabel('time since last pericenter passage [Myr]')
		ax1.set_ylabel('NH$_3$ temperature T [K]')
		ax1.set_xlim(-0.6,1.6)
		ax1.xaxis.set_minor_locator(MultipleLocator(0.1))
		ax1.set_ylim(0,)
		ax1.set_ylim(0,temp_info[1][temp])
		ax1.yaxis.set_minor_locator(MultipleLocator(10))
		ax1.text(0.17,0.05*temp_info[1][temp], r'$\tau_{ff} = 0.34$\,Myr', color='r', ha='center', va='bottom', **text_kwargs)
		ax1.plot([0.0,0.34], [0.03*temp_info[1][temp],0.03*temp_info[1][temp]], color='r', linewidth=2, **text_kwargs)
		ax1.bar(-0.25, 500, 0.5, bottom=None, color="grey", alpha=0.5, zorder=3)
		ax1.bar(0.25, 500, 0.75, bottom=None, color="grey", alpha=0.5, zorder=3)
		ax1.text(0.0, 0.95*temp_info[1][temp], "compression\ndominated", fontsize=10, ha='center', va='top', bbox=props, **text_kwargs)
		ax1.text(0.625, 0.95*temp_info[1][temp], "star formation\ndominated", fontsize=10, ha='center', va='top', bbox=props, **text_kwargs)

		# set Tkin sensitivity limit
		if (ttype == "kin"): 
			ax1.axhline(y=temp_info[2][temp], linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp_info[0][temp]+'}$ sensitivity limit', **limit_kwargs)

		# fit sequences
		fit_results = "slope:\n"
		fit_lim_lo = -0.25
		fit_lim_hi = 1.0
	
		# fit sequence
		to_fit = [[],[]]
		for row in np.arange(len(data)):
			if (phase_shifted[row] > fit_lim_lo) and (phase_shifted[row] < fit_lim_hi) and not (np.isnan(data['T'+temp_info[0][temp]+'_'+ttype][row])): 
				to_fit[0].append(phase_shifted[row])
				to_fit[1].append(data['T'+temp_info[0][temp]+'_'+ttype][row])
#		coeff, cov = np.polyfit(to_fit[0], to_fit[1], 1, cov=True)
#		ax1.plot([-2.0,-1.0],[x*coeff[0]+coeff[1] for x in [-2.0,-1.0]], linewidth=1.0, markersize=0.0, color='black', linestyle='--', label='fit $'+str(fit_lim_lo)+'<t<'+str(fit_lim_hi)+'$', **fit_kwargs)
			
		# fit binned sequence
		bin_fit = [[],[],[],[]]
		for time_bin in np.arange(fit_lim_lo,fit_lim_hi,0.01):
			bin_value = []
			for row in np.arange(len(to_fit[0])):
				if (to_fit[0][row] >= time_bin) and (to_fit[0][row] < time_bin+0.01):
					bin_value.append(to_fit[1][row])
			bin_value.sort()
			if (np.median(bin_value) > 5) and (len(bin_value) > 10):
				bin_fit[0].append(time_bin+0.01)
				bin_fit[1].append(np.median(bin_value))
				bin_fit[2].append(np.median(bin_value) - bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][0])	# lower "error" of 50% include range
				bin_fit[3].append(bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][-1] - np.median(bin_value))	# upper "error" of 50% include range
		coeff_med, cov_med = np.polyfit(bin_fit[0], bin_fit[1], 1, cov=True)
		fit_results += r"$\frac{dT}{dt} = "+str('{:3.1f}'.format(coeff_med[0])).rjust(5)+"$\\,K/Myr ($"+str(fit_lim_lo)+"<t [Myr] <"+str(fit_lim_hi)+"$)\n"
			
		# plot fits; only one legend entry
		ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', label='fit bin median', **fit_kwargs)
		ax1.scatter(bin_fit[0], bin_fit[1], label='0.01\,Myr bins median', s=10, marker="x", color='black', **text_kwargs)
#		ax1.errorbar(bin_fit[0], bin_fit[1], yerr=zip(*[bin_fit[2],bin_fit[3]]), fmt='none', linewidth=0.0, elinewidth=0.5, ecolor='black', **limit_kwargs)
			
		# fit results text box
		ax1.text(0.98, 0.03, fit_results[:-1], transform=ax1.transAxes, fontsize=10, va='bottom', ha='right', multialignment='left', bbox=props, **text_kwargs)
		
		# plot data
		scatter = ax1.scatter(phase_shifted, data['T'+temp_info[0][temp]+'_'+ttype], label='T$_{'+ttype+','+temp_info[0][temp]+'}$', s=1, lw=0, c='dimgrey', **scatter_kwargs)
		
		# contour plot
		contour_list = [[],[]]
		for row in np.arange(len(data)):
			if not (np.isnan(phase_shifted[row])) and not (np.isnan(data['T'+temp_info[0][temp]+'_'+ttype][row])):
				contour_list[0].append(phase_shifted[row])
				contour_list[1].append(data['T'+temp_info[0][temp]+'_'+ttype][row])

		H, xedges, yedges = np.histogram2d(contour_list[0], contour_list[1], bins=100, range=[[-2.5,2.5],[0,500]])
		extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
#		levels=[25,50,75,125,175]
		levels = [i*np.amax(H) for i in np.arange(0.1,1.0,0.1)]
		contours = ax1.contour(np.transpose(H),extent=extent, levels=levels, label="10% steps of max. density", **scatter_kwargs)

		# legend
		contour_legend = lines.Line2D([],[], color='b', lw=1, label='10% steps of max. density')
		leg = ax1.legend(fancybox=True, framealpha=0.8, fontsize=10)
		leg.set_zorder(12)
		
		plt.savefig("plots/phase_vs_T"+temp_info[0][temp]+"_"+ttype+".density.png", dpi=300)
		
############################################################

# T?? vs. phase (density, two fits)
print "plotting T(NH3) vs. phase ..."
for ttype in ["kin", "rot"]:
	for temp in np.arange(len(temp_info[0])):
		plt.clf()
		plt.cla()
		fig = plt.figure()
		ax1 = plt.subplot(1,1,1)
		ax1.set_xlabel('time since last pericenter passage [Myr]')
		ax1.set_ylabel('NH$_3$ temperature T [K]')
		ax1.set_xlim(-0.6,1.6)
		ax1.xaxis.set_minor_locator(MultipleLocator(0.1))
		ax1.set_ylim(0,temp_info[1][temp])
		ax1.yaxis.set_minor_locator(MultipleLocator(10))
		ax1.text(0.17,0.05*temp_info[1][temp], r'$\tau_{ff} = 0.34$\,Myr', color='r', ha='center', va='bottom', **text_kwargs)
		ax1.plot([0.0,0.34], [0.03*temp_info[1][temp],0.03*temp_info[1][temp]], color='r', linewidth=2, **text_kwargs)
		ax1.bar(-0.25, 500, 0.5, bottom=None, color="grey", alpha=0.5, zorder=3)
		ax1.bar(0.25, 500, 0.75, bottom=None, color="grey", alpha=0.5, zorder=3)
		ax1.text(0.0, 0.95*temp_info[1][temp], "compression\ndominated", fontsize=10, ha='center', va='top', bbox=props, **text_kwargs)
		ax1.text(0.625, 0.95*temp_info[1][temp], "star formation\ndominated", fontsize=10, ha='center', va='top', bbox=props, **text_kwargs)
	
		# set Tkin sensitivity limit
		if (ttype == "kin"): 
			ax1.axhline(y=temp_info[2][temp], linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp_info[0][temp]+'}$ sensitivity limit', **limit_kwargs)

		# fit sequences
		fit_ranges = [[-0.25,0.25],[0.25,0.9]]
		fit_results = "slopes:\n"
		for fit_range_index in np.arange(len(fit_ranges)):
			fit_lim_lo = fit_ranges[fit_range_index][0]
			fit_lim_hi = fit_ranges[fit_range_index][1]
	
			# fit sequence
			to_fit = [[],[]]
			for row in np.arange(len(data)):
				if (phase_shifted[row] > fit_lim_lo) and (phase_shifted[row] < fit_lim_hi) and not (np.isnan(data['T'+temp_info[0][temp]+'_'+ttype][row])): 
					to_fit[0].append(phase_shifted[row])
					to_fit[1].append(data['T'+temp_info[0][temp]+'_'+ttype][row])
#			coeff, cov = np.polyfit(to_fit[0], to_fit[1], 1, cov=True)
#			ax1.plot([-2.0,-1.0],[x*coeff[0]+coeff[1] for x in [-2.0,-1.0]], linewidth=1.0, markersize=0.0, color='black', linestyle='--', label='fit $'+str(fit_lim_lo)+'<t<'+str(fit_lim_hi)+'$', **fit_kwargs)
			
			# fit binned sequence
			bin_fit = [[],[],[],[]]
			for time_bin in np.arange(fit_lim_lo,fit_lim_hi,0.01):
				bin_value = []
				for row in np.arange(len(to_fit[0])):
					if (to_fit[0][row] >= time_bin) and (to_fit[0][row] < time_bin+0.01):
						bin_value.append(to_fit[1][row])
				bin_value.sort()
				if (np.median(bin_value) > 5) and (len(bin_value) > 10):
					bin_fit[0].append(time_bin+0.01)
					bin_fit[1].append(np.median(bin_value))
					bin_fit[2].append(np.median(bin_value) - bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][0])	# lower "error" of 50% include range
					bin_fit[3].append(bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][-1] - np.median(bin_value))	# upper "error" of 50% include range
			coeff_med, cov_med = np.polyfit(bin_fit[0], bin_fit[1], 1, cov=True)
			fit_results += r"$\frac{dT}{dt} = "+str('{:3.1f}'.format(coeff_med[0])).rjust(5)+"$\\,K/Myr ($"+str(fit_lim_lo)+"<t [Myr] <"+str(fit_lim_hi)+"$)\n"
			
			# plot fits; only one legend entry
			if (fit_range_index == 0):
				ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', label='fit bin median', **fit_kwargs)
				ax1.scatter(bin_fit[0], bin_fit[1], label='0.01\,Myr bins median', s=10, marker="x", color='black', **text_kwargs)
#				ax1.errorbar(bin_fit[0], bin_fit[1], yerr=zip(*[bin_fit[2],bin_fit[3]]), fmt='none', linewidth=0.0, elinewidth=0.5, ecolor='black', **limit_kwargs)
			else:
				ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', **fit_kwargs)
				ax1.scatter(bin_fit[0], bin_fit[1], s=10, marker="x", color='black', **text_kwargs)
#				ax1.errorbar(bin_fit[0], bin_fit[1], yerr=zip(*[bin_fit[2],bin_fit[3]]), fmt='none', linewidth=0.0, elinewidth=0.5, ecolor='black', **limit_kwargs)
			
		# fit results text box
		ax1.text(0.98, 0.03, fit_results[:-1], transform=ax1.transAxes, fontsize=10, va='bottom', ha='right', multialignment='left', bbox=props, **text_kwargs)
		
		# plot data
		scatter = ax1.scatter(phase_shifted, data['T'+temp_info[0][temp]+'_'+ttype], label='T$_{'+ttype+','+temp_info[0][temp]+'}$', s=1, lw=0, c='dimgrey', **scatter_kwargs)
		
		# contour plot
		contour_list = [[],[]]
		for row in np.arange(len(data)):
			if not (np.isnan(phase_shifted[row])) and not (np.isnan(data['T'+temp_info[0][temp]+'_'+ttype][row])):
				contour_list[0].append(phase_shifted[row])
				contour_list[1].append(data['T'+temp_info[0][temp]+'_'+ttype][row])

		H, xedges, yedges = np.histogram2d(contour_list[0], contour_list[1], bins=100, range=[[-2.5,2.5],[0,500]])
		extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
#		levels=[25,50,75,125,175]
		levels = [i*np.amax(H) for i in np.arange(0.1,1.0,0.1)]
		contours = ax1.contour(np.transpose(H),extent=extent, levels=levels, label="10% steps of max. density", **scatter_kwargs)

		# legend
		contour_legend = lines.Line2D([],[], color='b', lw=1, label='10% steps of max. density')
		leg = ax1.legend(fancybox=True, framealpha=0.8, fontsize=10)
		leg.set_zorder(12)
		
		plt.savefig("plots/phase_vs_T"+temp_info[0][temp]+"_"+ttype+".density.2fit.png", dpi=300)
		
		
############################################################

# opacity vs. time (density, fit)
print "plotting opacity vs. time ..."
for j in np.arange(1,7):
	J = str(j)
	plt.clf()
	plt.cla()
	fig = plt.figure()
	fig.subplots_adjust(top=0.8)
	ax1 = plt.subplot(1,1,1)
	ax2 = ax1.twiny()
	ax1.set_xlabel('time [Myr]')
	ax1.set_ylabel(r'opacity $\tau$')
	ax1.set_xlim(-2.5,2.5)
	ax1.xaxis.set_minor_locator(MultipleLocator(0.1))
	ax1.set_ylim(0, 8)
	ax1.set_yscale("linear")
	ax1.yaxis.set_minor_locator(MultipleLocator(2))
	ax2.set_xlim(-2.5,2.5)
	ax2.set_xticks([-2.0, -1.67, -1.52, -1.30, -1.0, 0.0, 1.0, 1.51, 2.0])
	ax2.set_xticklabels(['pericenter\nnear', 'Brick', 'cloud e/f', 'Sgr B2', 'apocenter\neast', 'pericenter\nfar', 'apocenter\nwest', 'Sgr C', 'pericenter\nnear'], rotation=90, ha='center', va='bottom')
	ax1.bar(-1.7, 500, 0.4, bottom=None, color="grey", alpha=0.5, zorder=3)
	props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
	ax1.text(-1.5, 28, "dust ridge", ha='center', bbox=props, **text_kwargs)
	ax1.text(0.17,12, r'$\tau_{ff} = 0.34$\,Myr', color='r', ha='center', va='bottom', **text_kwargs)
	ax1.plot([0.0,0.34], [10,10], color='r', linewidth=2, **text_kwargs)
	
	# fit sequences
	fit_ranges = [[-2.0,-1.0],[-0.3,0.05],[1.4,2.25]]
	fit_results = "slopes:\n"
	for fit_range_index in np.arange(len(fit_ranges)):
		fit_lim_lo = fit_ranges[fit_range_index][0]
		fit_lim_hi = fit_ranges[fit_range_index][1]

		# fit sequence
		to_fit = [[],[]]
		for row in np.arange(len(data)):
			if (data['time'][row] > fit_lim_lo) and (data['time'][row] < fit_lim_hi) and not (np.isnan(data['tau'+J][row])): 
				to_fit[0].append(data['time'][row])
				to_fit[1].append(data['tau'+J][row])
#		coeff, cov = np.polyfit(to_fit[0], to_fit[1], 1, cov=True)
#		ax1.plot([-2.0,-1.0],[x*coeff[0]+coeff[1] for x in [-2.0,-1.0]], linewidth=1.0, markersize=0.0, color='black', linestyle='--', label='fit $'+str(fit_lim_lo)+'<t<'+str(fit_lim_hi)+'$', **fit_kwargs)
		
		# fit binned sequence
		bin_fit = [[],[],[],[]]
		for time_bin in np.arange(fit_lim_lo,fit_lim_hi,0.01):
			bin_value = []
			for row in np.arange(len(to_fit[0])):
				if (to_fit[0][row] >= time_bin) and (to_fit[0][row] < time_bin+0.01):
					bin_value.append(to_fit[1][row])
			bin_value.sort()
			if (np.median(bin_value) < 6) and (len(bin_value) > 10):
				bin_fit[0].append(time_bin+0.01)
				bin_fit[1].append(np.median(bin_value))
				bin_fit[2].append(np.median(bin_value) - bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][0])	# lower "error" of 50% include range
				bin_fit[3].append(bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][-1] - np.median(bin_value))	# upper "error" of 50% include range
		coeff_med, cov_med = np.polyfit(bin_fit[0], bin_fit[1], 1, cov=True)
		fit_results += r"$\frac{d\tau}{dt} = "+str('{:3.1f}'.format(coeff_med[0])).rjust(5)+"$\\,Myr$^{-1}$ ($"+str(fit_lim_lo)+"<t [Myr] <"+str(fit_lim_hi)+"$)\n"
		
		# plot fits; only one legend entry
		if (fit_range_index == 0):
			ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', label='fit bin median', **fit_kwargs)
			ax1.scatter(bin_fit[0], bin_fit[1], label='0.01\,Myr bins median', s=10, marker="x", color='black', **text_kwargs)
#			ax1.errorbar(bin_fit[0], bin_fit[1], yerr=zip(*[bin_fit[2],bin_fit[3]]), fmt='none', linewidth=0.0, elinewidth=0.5, ecolor='black', **limit_kwargs)
		else:
			ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', **fit_kwargs)
			ax1.scatter(bin_fit[0], bin_fit[1], s=10, marker="x", color='black', **text_kwargs)
#			ax1.errorbar(bin_fit[0], bin_fit[1], yerr=zip(*[bin_fit[2],bin_fit[3]]), fmt='none', linewidth=0.0, elinewidth=0.5, ecolor='black', **limit_kwargs)
			
	# fit results text box
	ax1.text(0.68, 0.975, fit_results[:-1], transform=ax1.transAxes, fontsize=10, va='top', ha='right', multialignment='left', bbox=props)
		
	# plot data
	plotable = [[],[]]
	too_high = [[],[]]
	for row in np.arange(len(data)):
		if (data['tau'+J][row] < 6):
			plotable[0].append(data['time'][row])
			plotable[1].append(data['tau'+J][row])
		else:
			too_high[0].append(data['time'][row])
			too_high[1].append(data['tau'+J][row])
		
	scatter = ax1.scatter(plotable[0], plotable[1], label='NH$_3$ ('+J+','+J+r') ($\tau \le 6$)', s=1, lw=0, c='dimgrey', **scatter_kwargs)
	scatter2 = ax1.scatter(too_high[0], [6 for i in np.arange(len(too_high[0]))], label='NH$_3$ ('+J+','+J+r') ($\tau > 6$)', marker='^', s=2, lw=0, c='dimgrey', **scatter_kwargs)
		
	# contour plot
	contour_list = [[],[]]
	for row in np.arange(len(data)):
		if not (np.isnan(data['time'][row])) and not (np.isnan(data['tau'+J][row])):
			contour_list[0].append(data['time'][row])
			contour_list[1].append(data['tau'+J][row])

	H, xedges, yedges = np.histogram2d(contour_list[0], contour_list[1], bins=100, range=[[-2.5,2.5],[0,30]])
	extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
#	levels=[25,50,75,125,175]
	levels = [i*np.amax(H) for i in np.arange(0.1,1.0,0.1)]
	contours = ax1.contour(np.transpose(H),extent=extent, levels=levels, label="10% steps of max. density", **scatter_kwargs)

	# legend
	contour_legend = lines.Line2D([],[], color='b', lw=1, label='10% steps of max. density')
	leg = ax1.legend(fancybox=True, framealpha=0.8, fontsize=10)
	leg.set_zorder(12)
	
	plt.savefig("plots/time_vs_opacity"+J+".density.fit.png", dpi=300)


############################################################

# T?? vs. FWHM (density)
print "plotting T(NH3) vs. FWHM ..."
for ttype in ["kin"]: #["kin", "rot"]:
	for temp in ["12", "24", "45", "36"]:
		plt.clf()
		plt.cla()
		fig = plt.figure()
		ax1 = plt.subplot(1,1,1)
		ax1.set_xlabel('FWHM [km/s]')
		ax1.set_ylabel('NH$_3$ temperature T [K]')
		ax1.set_xlim(0,50)
		ax1.xaxis.set_minor_locator(MultipleLocator(2))
		ax1.set_ylim(0, 400)
		if (temp == "12"):
			ax1.set_ylim(0,80)
		if (temp == "24"):
			ax1.set_ylim(0,280)
		if (temp == "45"):
			ax1.set_ylim(0,280)
		if (temp == "36"):
			ax1.set_ylim(0,400)
		ax1.yaxis.set_minor_locator(MultipleLocator(10))
	
		# contour plot
		contour_list = [[],[]]
		for row in np.arange(len(data)):
			if (data['width'+temp[0]][row] > 3.5) and (data['T'+temp+'_'+ttype][row] > 10.0):
				contour_list[0].append(data['width'+temp[0]][row])
				contour_list[1].append(data['T'+temp+'_'+ttype][row])

		# plot data
		scatter = ax1.scatter(contour_list[0], contour_list[1], label='T$_{'+ttype+','+temp+'}$', s=1, lw=0, c='dimgrey', **scatter_kwargs)

		H, xedges, yedges = np.histogram2d(contour_list[0], contour_list[1], bins=100, range=[[0,50],[0,400]])
		extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
#		levels=[25,50,75,125,175]
		levels = [i*np.amax(H) for i in np.arange(0.1,1.0,0.1)]
		contours = ax1.contour(np.transpose(H),extent=extent, levels=levels, label="10% steps of max. density", **scatter_kwargs)

# fit without binning is a bad idea		
#		# fit
#		coeff, cov = np.polyfit(contour_list[0], contour_list[1], 1, cov=True)
#		ax1.plot([0,50],[x*coeff[0]+coeff[1] for x in [0,50]], linewidth=1.0, markersize=0.0, color='black', linestyle='--', label='linear fit', **fit_kwargs)
#		
#		fit_results = "slope:\n"
#		fit_results += r"$\frac{dT}{dFWHM} = "+str('{:3.1f}'.format(coeff[0])).rjust(5)+"$\\,K/(km/s)\n"
#		ax1.text(0.68, 0.975, fit_results[:-1], transform=ax1.transAxes, fontsize=10, va='top', ha='right', multialignment='left', bbox=props)

	# fit sequences
	fit_ranges = [[-2.0,-1.0],[-0.3,0.05],[1.75,2.25]]
	fit_results = "slopes:\n"
	for fit_range_index in np.arange(len(fit_ranges)):
		fit_lim_lo = fit_ranges[fit_range_index][0]
		fit_lim_hi = fit_ranges[fit_range_index][1]

		# fit sequence
		to_fit = [[],[]]
		for row in np.arange(len(data)):
			if (data['time'][row] > fit_lim_lo) and (data['time'][row] < fit_lim_hi) and not (np.isnan(data['width'+line][row])): 
				to_fit[0].append(data['time'][row])
				to_fit[1].append(data['width'+line][row])
#		coeff, cov = np.polyfit(to_fit[0], to_fit[1], 1, cov=True)
#		ax1.plot([-2.0,-1.0],[x*coeff[0]+coeff[1] for x in [-2.0,-1.0]], linewidth=1.0, markersize=0.0, color='black', linestyle='--', label='fit $'+str(fit_lim_lo)+'<t<'+str(fit_lim_hi)+'$', **fit_kwargs)
		
		# fit binned sequence
		bin_fit = [[],[]]
		for time_bin in np.arange(fit_lim_lo,fit_lim_hi,0.01):
			bin_value = []
			for row in np.arange(len(to_fit[0])):
				if (to_fit[0][row] >= time_bin) and (to_fit[0][row] < time_bin+0.01):
					bin_value.append([to_fit[1][row]])
			if (np.median(bin_value) > 5) and (len(bin_value) > 10):
				bin_fit[0].append(time_bin+0.01)
				bin_fit[1].append(np.median(bin_value))
		coeff_med, cov_med = np.polyfit(bin_fit[0], bin_fit[1], 1, cov=True)
		
		# plot fits; only one legend entry
		if (fit_range_index == 0):
			ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', label='fit bin median', **fit_kwargs)
			ax1.scatter(bin_fit[0], bin_fit[1], label='0.01\,Myr bins median', s=10, marker="x", color='black', **text_kwargs)
		else:
			ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', **fit_kwargs)
			ax1.scatter(bin_fit[0], bin_fit[1], s=10, marker="x", color='black', **text_kwargs)


		# set Tkin sensitivity limit
		if (temp == "12") and (ttype == "kin"): 
			ax1.axhline(y=60, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		if (temp == "24") and (ttype == "kin"): 
			ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		if (temp == "45") and (ttype == "kin"): 
			ax1.axhline(y=200, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)
		if (temp == "36") and (ttype == "kin"): 
			ax1.axhline(y=300, linewidth=1.0, linestyle='--', color='grey', label='T$_{'+temp+'}$ sensitivity limit', **limit_kwargs)

		# legend
		contour_legend = lines.Line2D([],[], color='b', lw=1, label='10% steps of max. density')
		leg = ax1.legend(fancybox=True, framealpha=0.8, fontsize=10)
		leg.set_zorder(12)
		
		plt.savefig("plots/FWHM("+temp[0]+","+temp[0]+")_vs_T"+temp+"_"+ttype+".density.png", dpi=300)


############################################################

# column density vs. time (density, fit)
print "plotting column density vs. time ..."
for j in np.arange(1,7):
	J = str(j)
	plt.clf()
	plt.cla()
	fig = plt.figure()
	fig.subplots_adjust(top=0.8)
	ax1 = plt.subplot(1,1,1)
	ax2 = ax1.twiny()
	ax1.set_xlabel('time [Myr]')
	ax1.set_ylabel(r'log$_{10}$ column density N [cm$^{-2}$]')
	ax1.set_yscale("linear")
	ax1.set_xlim(-2.5,2.5)
	ax1.xaxis.set_minor_locator(MultipleLocator(0.1))
	ax1.set_ylim(13, 17)
	ax2.set_xlim(-2.5,2.5)
	ax2.set_xticks([-2.0, -1.67, -1.52, -1.30, -1.0, 0.0, 1.0, 1.51, 2.0])
	ax2.set_xticklabels(['pericenter\nnear', 'Brick', 'cloud e/f', 'Sgr B2', 'apocenter\neast', 'pericenter\nfar', 'apocenter\nwest', 'Sgr C', 'pericenter\nnear'], rotation=90, ha='center', va='bottom')
	ax1.bar(-1.7, 500, 0.4, bottom=None, color="grey", alpha=0.5, zorder=3)
	props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
	ax1.text(-1.5, 28, "dust ridge", ha='center', bbox=props, **text_kwargs)
	ax1.text(0.17,12, r'$\tau_{ff} = 0.34$\,Myr', color='r', ha='center', va='bottom', **text_kwargs)
	ax1.plot([0.0,0.34], [10,10], color='r', linewidth=2, **text_kwargs)
		
	# fit sequences
	fit_ranges = [[-2.0,-1.0],[-0.3,0.05],[1.4,2.25]]
	fit_results = "slopes:\n"
	for fit_range_index in np.arange(len(fit_ranges)):
		fit_lim_lo = fit_ranges[fit_range_index][0]
		fit_lim_hi = fit_ranges[fit_range_index][1]
	
		# fit sequence
		to_fit = [[],[]]
		for row in np.arange(len(data)):
			if (data['time'][row] > fit_lim_lo) and (data['time'][row] < fit_lim_hi) and not (np.isnan(np.log10(data['N'+J][row]))): 
				to_fit[0].append(data['time'][row])
				to_fit[1].append(np.log10(data['N'+J][row]))
			
		# fit binned sequence
		bin_fit = [[],[],[],[]]
		for time_bin in np.arange(fit_lim_lo,fit_lim_hi,0.01):
			bin_value = []
			for row in np.arange(len(to_fit[0])):
				if (to_fit[0][row] >= time_bin) and (to_fit[0][row] < time_bin+0.01):
					bin_value.append(to_fit[1][row])
			bin_value.sort()
			if (len(bin_value) > 10):
				bin_fit[0].append(time_bin+0.01)
				bin_fit[1].append(np.median(bin_value))
				bin_fit[2].append(np.median(bin_value) - bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][0])
				bin_fit[3].append(bin_value[int(len(bin_value) * .25) : int(len(bin_value) * .75)][-1] - np.median(bin_value))
		coeff_med, cov_med = np.polyfit(bin_fit[0], bin_fit[1], 1, cov=True)
		fit_results += r"$\frac{dN}{dt} = "+str(sci_notation(coeff_med[0],2)).rjust(5)+"$\\,cm$^{-2}$\,Myr$^{-1}$ ($"+str(fit_lim_lo)+"<t [Myr] <"+str(fit_lim_hi)+"$)\n"
			
		# plot fits; only one legend entry
		if (fit_range_index == 0):
			ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', label='fit bin median', **fit_kwargs)
			ax1.scatter(bin_fit[0], bin_fit[1], label='0.01\,Myr bins median', s=10, marker="x", color='black', **text_kwargs)
		else:
			ax1.plot([fit_lim_lo,fit_lim_hi],[x*coeff_med[0]+coeff_med[1] for x in [fit_lim_lo,fit_lim_hi]], linewidth=1.0, markersize=0.0, color='black', **fit_kwargs)
			ax1.scatter(bin_fit[0], bin_fit[1], s=10, marker="x", color='black', **text_kwargs)
				
	# fit results text box
	ax1.text(0.68, 0.975, fit_results[:-1], transform=ax1.transAxes, fontsize=10, va='top', ha='right', multialignment='left', bbox=props)
		
	# plot data
	scatter = ax1.scatter(data['time'], np.log10(data['N'+J]), label='NH$_3$ ('+J+','+J+r')', s=1, lw=0, c='dimgrey', **scatter_kwargs)
	
	# contour plot
	contour_list = [[],[]]
	for row in np.arange(len(data)):
		if not (np.isnan(data['time'][row])) and not (np.isnan(np.log10(data['N'+J][row]))):
			contour_list[0].append(data['time'][row])
			contour_list[1].append(np.log10(data['N'+J][row]))
	
	H, xedges, yedges = np.histogram2d(contour_list[0], contour_list[1], bins=100, range=[[-2.5,2.5],[13,17]])
	extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
	levels = [i*np.amax(H) for i in np.arange(0.1,1.0,0.1)]
	contours = ax1.contour(np.transpose(H),extent=extent, levels=levels, label="10% steps of max. density", **scatter_kwargs)
	
	# legend
	contour_legend = lines.Line2D([],[], color='b', lw=1, label='10% steps of max. density')
	leg = ax1.legend(fancybox=True, framealpha=0.8, fontsize=10)
	leg.set_zorder(12)
	
	plt.savefig("plots/time_vs_col_dens"+J+".density.fit.png", dpi=300)

############################################################
