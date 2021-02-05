import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit
from scipy.misc import factorial
rc('text', usetex=True)
plt.ioff()


# ####################################
# # ******************************** #
# # * SF sequence analysis script  * #
# # ******************************** #
# ####################################


# Get a bunch of histograms for illustration

# parameters:
fit_dir = 'fits_NH3/'

############################################################

#os.system('mkdir plots/')

# load data
data1 = np.genfromtxt(fit_dir+'NH3_1-1.fit_list.txt',
	usecols = (0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21),
	names   = "x, y, good, area, area_err, v, v_err, width, width_err, tau, tau_err, baserms, linerms, T_L, T_L_err, Tdv, Tdv_err, N, N_err, yaxis, yaxis_err",
	dtype   = "i4, i4, i4, f8,f8,f8,f8,f8,f8,f8,f8,f8,f8, f8,f8,f8,f8,f8,f8,f8,f8"
	)
data2 = np.genfromtxt(fit_dir+'NH3_2-2.fit_list.txt',
	usecols = (0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21),
	names   = "x, y, good, area, area_err, v, v_err, width, width_err, tau, tau_err, baserms, linerms, T_L, T_L_err, Tdv, Tdv_err, N, N_err, yaxis, yaxis_err",
	dtype   = "i4, i4, i4, f8,f8,f8,f8,f8,f8,f8,f8,f8,f8, f8,f8,f8,f8,f8,f8,f8,f8"
	)
data3 = np.genfromtxt(fit_dir+'NH3_3-3.fit_list.txt',
	usecols = (0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21),
	names   = "x, y, good, area, area_err, v, v_err, width, width_err, tau, tau_err, baserms, linerms, T_L, T_L_err, Tdv, Tdv_err, N, N_err, yaxis, yaxis_err",
	dtype   = "i4, i4, i4, f8,f8,f8,f8,f8,f8,f8,f8,f8,f8, f8,f8,f8,f8,f8,f8,f8,f8"
	)
data4 = np.genfromtxt(fit_dir+'NH3_4-4.fit_list.txt',
	usecols = (0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21),
	names   = "x, y, good, area, area_err, v, v_err, width, width_err, tau, tau_err, baserms, linerms, T_L, T_L_err, Tdv, Tdv_err, N, N_err, yaxis, yaxis_err",
	dtype   = "i4, i4, i4, f8,f8,f8,f8,f8,f8,f8,f8,f8,f8, f8,f8,f8,f8,f8,f8,f8,f8"
	)
data5 = np.genfromtxt(fit_dir+'NH3_5-5.fit_list.txt',
	usecols = (0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21),
	names   = "x, y, good, area, area_err, v, v_err, width, width_err, tau, tau_err, baserms, linerms, T_L, T_L_err, Tdv, Tdv_err, N, N_err, yaxis, yaxis_err",
	dtype   = "i4, i4, i4, f8,f8,f8,f8,f8,f8,f8,f8,f8,f8, f8,f8,f8,f8,f8,f8,f8,f8"
	)
data6 = np.genfromtxt(fit_dir+'NH3_6-6.fit_list.txt',
	usecols = (0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21),
	names   = "x, y, good, area, area_err, v, v_err, width, width_err, tau, tau_err, baserms, linerms, T_L, T_L_err, Tdv, Tdv_err, N, N_err, yaxis, yaxis_err",
	dtype   = "i4, i4, i4, f8,f8,f8,f8,f8,f8,f8,f8,f8,f8, f8,f8,f8,f8,f8,f8,f8,f8"
	)

############################################################


# get number of good and bad fits
npix1 = float(len(data1['good']))
npix2 = float(len(data2['good']))
npix3 = float(len(data3['good']))
npix4 = float(len(data4['good']))
npix5 = float(len(data5['good']))
npix6 = float(len(data6['good']))
bad_good1 = np.bincount(data1['good'])
bad_good2 = np.bincount(data2['good'])
bad_good3 = np.bincount(data3['good'])
bad_good4 = np.bincount(data4['good'])
bad_good5 = np.bincount(data5['good'])
bad_good6 = np.bincount(data6['good'])
good_temp1 = np.count_nonzero(~np.isnan(data1['yaxis']))
good_temp2 = np.count_nonzero(~np.isnan(data2['yaxis']))
good_temp3 = np.count_nonzero(~np.isnan(data3['yaxis']))
good_temp4 = np.count_nonzero(~np.isnan(data4['yaxis']))
good_temp5 = np.count_nonzero(~np.isnan(data5['yaxis']))
good_temp6 = np.count_nonzero(~np.isnan(data6['yaxis']))
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print " NH3 fit results:"
print " \t\t(1,1)\t(2,2)\t(3,3)\t(4,4)\t(5,5)\t(6,6)"
print "-------------------------------------------------------------"
print " all fits:\t24383\t24383\t24383\t24383\t24383\t24383"
print " fitted:\t"+str(int(npix1))+"\t"+str(int(npix2))+"\t"+str(int(npix3))+"\t"+str(int(npix4))+"\t"+str(int(npix5))+"\t"+str(int(npix6))
print " % success:\t"+str('{:.1f}'.format(npix1/24383*100.0))+"\t"+str('{:.1f}'.format(npix2/24383*100.0))+"\t"+str('{:.1f}'.format(npix3/24383*100.0))+"\t"+str('{:.1f}'.format(npix4/24383*100.0))+"\t"+str('{:.1f}'.format(npix5/24383*100.0))+"\t"+str('{:.1f}'.format(npix6/24383*100.0))
print " thereof are ..."
print " bad fits:\t"+str(bad_good1[0])+"\t"+str(bad_good2[0])+"\t"+str(bad_good3[0])+"\t"+str(bad_good4[0])+"\t"+str(bad_good5[0])+"\t"+str(bad_good6[0])
print " good fits:\t"+str(bad_good1[1])+"\t"+str(bad_good2[1])+"\t"+str(bad_good3[1])+"\t"+str(bad_good4[1])+"\t"+str(bad_good5[1])+"\t"+str(bad_good6[1])
print " % success:\t"+str('{:.1f}'.format(bad_good1[1]/npix1*100.0))+"\t"+str('{:.1f}'.format(bad_good2[1]/npix2*100.0))+"\t"+str('{:.1f}'.format(bad_good3[1]/npix3*100.0))+"\t"+str('{:.1f}'.format(bad_good4[1]/npix4*100.0))+"\t"+str('{:.1f}'.format(bad_good5[1]/npix5*100.0))+"\t"+str('{:.1f}'.format(bad_good6[1]/npix6*100.0))
print " thereof could be derived ..."
print " temperature:\t"+str(good_temp1)+"\t"+str(good_temp2)+"\t"+str(good_temp3)+"\t"+str(good_temp4)+"\t"+str(good_temp5)+"\t"+str(good_temp6)
print " % success:\t"+str('{:.1f}'.format(good_temp1/npix1*100.0))+"\t"+str('{:.1f}'.format(good_temp2/npix3*100.0))+"\t"+str('{:.1f}'.format(good_temp3/npix3*100.0))+"\t"+str('{:.1f}'.format(good_temp4/npix4*100.0))+"\t"+str('{:.1f}'.format(good_temp5/npix5*100.0))+"\t"+str('{:.1f}'.format(good_temp6/npix6*100.0))
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"


############################################################

# plot histograms of velocity, line width, opacity, column density
alevel = 0.75	# alpha channel in plots


############################################################

# velocity
plt.clf()
plt.cla()
fig = plt.figure()
ax1 = plt.subplot(1,1,1)
n1, bins1, patches1 = ax1.hist(data1['v'], bins=50, label='(1,1)', facecolor='darkblue', alpha=alevel)
n2, bins2, patches2 = ax1.hist(data2['v'], bins=50, label='(2,2)', facecolor='mediumblue', alpha=alevel)
n3, bins3, patches3 = ax1.hist(data3['v'], bins=50, label='(3,3)', facecolor='blue', alpha=alevel)
n4, bins4, patches4 = ax1.hist(data4['v'], bins=50, label='(4,4)', facecolor='dodgerblue', alpha=alevel)
n5, bins5, patches5 = ax1.hist(data5['v'], bins=50, label='(5,5)', facecolor='skyblue', alpha=alevel)
n6, bins6, patches6 = ax1.hist(data6['v'], bins=50, label='(6,6)', facecolor='cyan', alpha=alevel)
ax1.set_xlabel(r'v$_{lsr}$ [km/s]')
ax1.set_ylabel('pixel count')
ax1.set_title('Velocity distribution of fitted pixels')
ax1.set_xlim(-300,300)
ax1.grid(True)
ax1.legend(loc=2)
plt.savefig("plots/hist.velocity.png", dpi=300)


############################################################

# line width
plt.clf()
plt.cla()
fig = plt.figure()
ax1 = plt.subplot(1,1,1)
n1, bins1, patches1 = ax1.hist(data1['width'], bins=50, range=(0,100), label='(1,1)', facecolor='darkblue', alpha=alevel)
n2, bins2, patches2 = ax1.hist(data2['width'], bins=50, range=(0,100), label='(2,2)', facecolor='mediumblue', alpha=alevel)
n3, bins3, patches3 = ax1.hist(data3['width'], bins=50, range=(0,100), label='(3,3)', facecolor='blue', alpha=alevel)
n4, bins4, patches4 = ax1.hist(data4['width'], bins=50, range=(0,100), label='(4,4)', facecolor='dodgerblue', alpha=alevel)
n5, bins5, patches5 = ax1.hist(data5['width'], bins=50, range=(0,100), label='(5,5)', facecolor='skyblue', alpha=alevel)
n6, bins6, patches6 = ax1.hist(data6['width'], bins=50, range=(0,100), label='(6,6)', facecolor='cyan', alpha=alevel)
ax1.set_xlabel(r'v$_{lsr}$ [km/s]')
ax1.set_ylabel('pixel count')
ax1.set_title('FWHM distribution of fitted pixels')
ax1.set_xlim(0,80)
ax1.grid(True)
ax1.legend(loc=1)
plt.savefig("plots/hist.line_width.png", dpi=300)


############################################################

# line width normalized
plt.clf()
plt.cla()
fig = plt.figure()
ax1 = plt.subplot(1,1,1)
#n1, bins1, patches1 = ax1.hist(data1['width'], bins=50, range=(0,100), normed=True, histtype='step', label='(1,1)', linewidth=2.0, color='mediumblue')
#n2, bins2, patches2 = ax1.hist(data2['width'], bins=50, range=(0,100), normed=True, histtype='step', label='(2,2)', linewidth=2.0, color='cyan')
#n3, bins3, patches3 = ax1.hist(data3['width'], bins=50, range=(0,100), normed=True, histtype='step', label='(3,3)', linewidth=2.0, color='darkgreen')
#n4, bins4, patches4 = ax1.hist(data4['width'], bins=50, range=(0,100), normed=True, histtype='step', label='(4,4)', linewidth=2.0, color='gold')
#n5, bins5, patches5 = ax1.hist(data5['width'], bins=50, range=(0,100), normed=True, histtype='step', label='(5,5)', linewidth=2.0, color='darkorange')
#n6, bins6, patches6 = ax1.hist(data6['width'], bins=50, range=(0,100), normed=True, histtype='step', label='(6,6)', linewidth=2.0, color='crimson')
entries, bin_edges, patches_all = ax1.hist([data1['width'],data2['width'],data3['width'],data4['width'],data5['width'],data6['width']], 
	bins=30, range=(0,60), normed=True, histtype='bar', label=['(1,1)','(2,2)','(3,3)','(4,4)','(5,5)','(6,6)'], linewidth=0.5, 
#	color=['darkblue','mediumblue','blue','dodgerblue','skyblue','cyan'])
	color=['blue','cyan','lime','yellow','orange','red'])

ax1.set_xlabel(r'v$_{lsr}$ [km/s]')
ax1.set_ylabel('normalised pixel count')
ax1.set_title('FWHM distribution of fitted pixels')
ax1.set_xlim(0,60)
ax1.grid(True)
ax1.legend(loc=1)
plt.savefig("plots/hist.line_width.norm.png", dpi=300)


############################################################

# opacity
plt.clf()
plt.cla()
fig = plt.figure()
ax1 = plt.subplot(1,1,1)
n1, bins1, patches1 = ax1.hist(data1['tau'], bins=40, range=(0,20), label='(1,1)', histtype='step', color='blue')
n2, bins2, patches2 = ax1.hist(data2['tau'], bins=40, range=(0,20), label='(2,2)', histtype='step', color='cyan')
n3, bins3, patches3 = ax1.hist(data3['tau'], bins=40, range=(0,20), label='(3,3)', histtype='step', color='lime')
n4, bins4, patches4 = ax1.hist(data4['tau'], bins=40, range=(0,20), label='(4,4)', histtype='step', color='yellow')
n5, bins5, patches5 = ax1.hist(data5['tau'], bins=40, range=(0,20), label='(5,5)', histtype='step', color='orange')
n6, bins6, patches6 = ax1.hist(data6['tau'], bins=40, range=(0,20), label='(6,6)', histtype='step', color='red')
ax1.set_xlabel(r'opacity')
ax1.set_ylabel('pixel count')
ax1.set_title('Opacity distribution of fitted pixels')
ax1.set_xlim(0,20)
ax1.xaxis.set_minor_locator(MultipleLocator(2))
ax1.grid(True)
ax1.legend(loc=1)
plt.savefig("plots/hist.opacity.png", dpi=300)


############################################################

# opacity
plt.clf()
plt.cla()
fig = plt.figure()
ax1 = plt.subplot(1,1,1)
n1, bins1, patches1 = ax1.hist(data1['tau'], bins=40, range=(0,20), label='(1,1)', facecolor='darkblue', alpha=alevel)
n2, bins2, patches2 = ax1.hist(data2['tau'], bins=40, range=(0,20), label='(2,2)', facecolor='mediumblue', alpha=alevel)
n3, bins3, patches3 = ax1.hist(data3['tau'], bins=40, range=(0,20), label='(3,3)', facecolor='blue', alpha=alevel)
n4, bins4, patches4 = ax1.hist(data4['tau'], bins=40, range=(0,20), label='(4,4)', facecolor='dodgerblue', alpha=alevel)
n5, bins5, patches5 = ax1.hist(data5['tau'], bins=40, range=(0,20), label='(5,5)', facecolor='skyblue', alpha=alevel)
n6, bins6, patches6 = ax1.hist(data6['tau'], bins=40, range=(0,20), label='(6,6)', facecolor='cyan', alpha=alevel)
ax1.set_xlabel(r'opacity')
ax1.set_ylabel('pixel count')
ax1.set_title('Opacity distribution of fitted pixels')
ax1.set_xlim(0,20)
ax1.xaxis.set_minor_locator(MultipleLocator(2))
ax1.grid(True)
ax1.legend(loc=1)
plt.savefig("plots/hist.opacity2.png", dpi=300)


############################################################

# opacity
plt.clf()
plt.cla()
fig = plt.figure()
ax1 = plt.subplot(1,1,1)
n_all, bins_all, patches_all = ax1.hist([data1['tau'],data2['tau'],data3['tau'],data4['tau'],data5['tau'],data6['tau']], 
	bins=40, range=(0,20), normed=False, histtype='bar', label=['(1,1)','(2,2)','(3,3)','(4,4)','(5,5)','(6,6)'], linewidth=0.5, 
#	color=['darkblue','mediumblue','blue','dodgerblue','skyblue','cyan'])
	color=['blue','cyan','lime','yellow','orange','red'])
ax1.set_xlabel(r'opacity')
ax1.set_ylabel('normalised pixel count')
ax1.set_title('Opacity distribution of fitted pixels')
ax1.set_xlim(0,20)
ax1.xaxis.set_minor_locator(MultipleLocator(2))
ax1.grid(True)
ax1.legend(loc=1)
plt.savefig("plots/hist.opacity3.png", dpi=300)


############################################################

# opacity normalized
plt.clf()
plt.cla()
fig = plt.figure()
ax1 = plt.subplot(1,1,1)
n_all, bins_all, patches_all = ax1.hist([data1['tau'],data2['tau'],data3['tau'],data4['tau'],data5['tau'],data6['tau']], 
	bins=40, range=(0,20), normed=True, histtype='bar', label=['(1,1)','(2,2)','(3,3)','(4,4)','(5,5)','(6,6)'], linewidth=0.5, 
#	color=['darkblue','mediumblue','blue','dodgerblue','skyblue','cyan'])
	color=['blue','cyan','lime','yellow','orange','red'])
ax1.set_xlabel(r'opacity')
ax1.set_ylabel('normalised pixel count')
ax1.set_title('Opacity distribution of fitted pixels')
ax1.set_xlim(0,20)
ax1.xaxis.set_minor_locator(MultipleLocator(2))
ax1.grid(True)
ax1.legend(loc=1)
plt.savefig("plots/hist.opacity.norm.png", dpi=300)


############################################################

# column density
plt.clf()
plt.cla()
fig = plt.figure()
ax1 = plt.subplot(1,1,1)
n1, bins1, patches1 = ax1.hist(data1['N']/10**16, bins=50, range=(-0.1,0.8), label='(1,1)', facecolor='darkblue', alpha=alevel)
n2, bins2, patches2 = ax1.hist(data2['N']/10**16, bins=50, range=(-0.1,0.8), label='(2,2)', facecolor='mediumblue', alpha=alevel)
n3, bins3, patches3 = ax1.hist(data3['N']/10**16, bins=50, range=(-0.1,0.8), label='(3,3)', facecolor='blue', alpha=alevel)
n4, bins4, patches4 = ax1.hist(data4['N']/10**16, bins=50, range=(-0.1,0.8), label='(4,4)', facecolor='dodgerblue', alpha=alevel)
n5, bins5, patches5 = ax1.hist(data5['N']/10**16, bins=50, range=(-0.1,0.8), label='(5,5)', facecolor='skyblue', alpha=alevel)
n6, bins6, patches6 = ax1.hist(data6['N']/10**16, bins=50, range=(-0.1,0.8), label='(6,6)', facecolor='cyan', alpha=alevel)
ax1.set_xlabel(r'N [$10^{16}$\,cm$^{-2}$]')
ax1.set_ylabel('pixel count')
ax1.set_title('Column density distribution of fitted pixels')
ax1.set_xlim(-0.1,1.8)
ax1.grid(True)
ax1.legend(loc=1)
plt.savefig("plots/hist.col_dens.png", dpi=300)

############################################################
