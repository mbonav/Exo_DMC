#!/usr/bin/env python
# coding: utf-8

__author__ = 'Mariangela Bonavita'
__version__ = 'v1.0-alpha'
__all__ = ['exodmc']


import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib import colors as mcolors
from matplotlib.ticker import ScalarFormatter
from scipy import interpolate
from numpy import random as rn
import scipy.ndimage as ndimage
import time
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)




##########################################
###              EXO-DMC               ###
### EXOplanet Detection Map Calculator ###
##########################################

"""



Class: exodmc

Monte Carlo tool for the statistical analysis of exoplanet surveys results.
Combines the information on the target stars with the instrument detection limits
to estimate the probability of detection of a given synthetic planet population,
ultimately generating detection probability maps.

Parameters:
star_ID: list, required. ID of the target(s)
star_dist: list, required. Target(s) distance, in pc

Methods:

1) set_grid
Change range or resolution of grid over companions are generated.
	- x_min: float, optional argument, lower limit for grid x axis (default = 1)
	- x_max: float, optional argument, upper limit for grid x axis (default = 1000)
	- nx: int, optional argument, number of steps in the grid x axis (default = 100)
	- xlog: boolean, optional argument. If True the x axis will be uniformly spaced in log
	- y_min: float, optional argument, lower limit for grid y axis (default = 0.5)
	- y_max: float, optional argument, upper limit for grid y axis (default = 75)
	- ny: int, optional argument, number of steps in the grid y axis (default = 100)
	- ylog: boolean, optional argument. If True the y axis will be uniformly spaced in log
	- ngen: float, optional argument, default=1000. Number of orbital elements sets to be generated for each point in the grid.
		All orbital parameters are uniformly distributed by default, except for the eccentricity.
	- e_dist: string, optional. Desired eccentricity distribution. Can be uniform ('uni') or Gaussian ('gauss', default)
	- e_mu: float, optional, mean of the gaussian eccentricity distribution, default is 0.0
	- e_sigma: float, optional, sigma of the gaussian eccentricity distribution, default is 0.3

2) DImode
	Estimates the detection probability map for DI data.
	Parameters:
	- xlim: list, required. Detection limit(s) projected separations (arcsec)
	- ylim: list, required. Minimum detectable mass (in Mjup) at each xlim.
	- lxunit: string, optional. Unit for xlim, default is arcseconds ('as').
		Can also be set to 'au' or 'mas'
	- verbose: if True (default), the code provide the runtime for each target
	- plot: if True (default) a .png file with the detection probability map is produced for each target

"""

class exodmc(object):

	def __init__(self, star_ID, star_dist):

		if len(np.shape(star_ID)) <= 1: star_ID=[star_ID]
		if len(np.shape(star_dist)) <= 1: star_dist=[star_dist]


		self.ID = star_ID
		self.dpc = star_dist # in pc
		self.set_grid()

	def set_grid(self, x_min=0.1, x_max=1000., nx=100,  logx=False, y_min=0.1, y_max=100., ny=100, logy=False, ngen=1000, e_dist='gauss', mu=0, sigma=0.3):

		self.x_min = x_min
		self.x_max = x_max
		self.x_nsteps = nx
		self.logx = logx

		self.y_min = y_min
		self.y_max = y_max
		self.y_nsteps = ny
		self.logy = logy

		self.norb = ngen
		self.e_dist = e_dist
		self.e_mu = mu
		self.e_sigma = sigma

		self.sma = np.linspace(self.x_min, self.x_max, self.x_nsteps)
		if self.logx is True: self.sma = np.logspace(np.log10(self.x_min), np.log10(self.x_max), self.x_nsteps)

		self.M2 = np.linspace(self.y_min, self.y_max, self.y_nsteps)	# range of M2 in Mjup
		if self.logy is True: self.M2 = np.logspace(np.log10(self.y_min), np.log10(self.y_max), self.y_nsteps)

		if e_dist=='gauss': self.ecc = np.abs(rn.normal(mu, sigma, self.norb))
		else: self.ecc = rn.random_sample(self.norb)
		self.exq = np.sqrt(1-self.ecc*self.ecc)
		self.Omega_Node = rn.random_sample(self.norb)*2.*np.pi # Longitude of node ranges between 0 and 2pi
		self.Omega_Peri = rn.random_sample(self.norb)*2.*np.pi # Longitude of periastron ranges between 0 and 2pi
		self.omega = self.Omega_Peri-self.Omega_Node
		cosi=2*rn.random_sample(self.norb) -1.
		self.irad=np.arccos(cosi)
		self.T0 = rn.random_sample(self.norb) # T peri in fraction of period

		# Thiele—Innes elements

		A1=(np.cos(self.Omega_Peri)*np.cos(self.Omega_Node))-(np.sin(self.Omega_Peri)*np.sin(self.Omega_Node)*np.cos(self.irad))
		B1=(np.cos(self.Omega_Peri)*np.sin(self.Omega_Node))+(np.sin(self.Omega_Peri)*np.cos(self.Omega_Node)*np.cos(self.irad))
		F1=(-1*np.sin(self.Omega_Peri)*np.cos(self.Omega_Node))-(np.cos(self.Omega_Peri)*np.sin(self.Omega_Node)*np.cos(self.irad))
		G1=(-1*np.sin(self.Omega_Peri)*np.sin(self.Omega_Node))+(np.cos(self.Omega_Peri)*np.cos(self.Omega_Node)*np.cos(self.irad))

		self.M=rn.random_sample(self.norb)*2*np.pi # mean anomaly
		E0=self.M + np.sin(self.M) * self.ecc + ((self.ecc**2)/2)*np.sin(2*self.M)
		M0=E0-self.ecc*np.sin(E0)
		self.E1=E0 + (self.M-M0)/1 - self.ecc*np.cos(E0) # eccentric anomaly
		self.nurad=2*np.arctan((np.sqrt(self.exq))*np.tan(self.E1/2.))

		x1=np.cos(self.nurad)-self.ecc[np.newaxis,:]
		y1=np.sqrt(1-self.ecc[np.newaxis,:])*np.sin(self.nurad)
		y2=(B1[np.newaxis,:]*x1 + G1[np.newaxis,:]*y1)
		x2=(A1[np.newaxis,:]*x1 + F1[np.newaxis,:]*y1)

		# radius vector and projected separation (arcsecs)
		self.rad=(np.sqrt(x2**2 + y2**2)).T
		self.rho=(self.rad[:,np.newaxis]*(self.sma[:,np.newaxis]/self.dpc)).T
		#self.rho=self.rad[:,np.newaxis]*(self.sma/self.dpc) # projected separation in AU



	def DImode(self, xlim, ylim, lxunit='as', lyunit='Mjup', verbose=True, plot=True, savefig=True):

		if len(np.shape(xlim)) is 1: xlim=[np.array(xlim)]
		if len(np.shape(ylim)) is 1: ylim=[np.array(ylim)]

		ns=np.size(self.dpc)
		detmap = []
		self.detflag = []
		for ll in range(ns):
			start = time.time()
			det=np.zeros((self.x_nsteps, self.y_nsteps, self.norb))

			if lxunit == 'au': xlim[ll] = xlim[ll]/self.dpc[ll]
			if lxunit == 'mas': xlim[ll] = xlim[ll]/1000.

			s=np.array(np.where(ylim[ll] < self.y_max))
			if np.size(s) > 1:
				max_mass=np.nanmax(ylim[ll][s])
				min_mass=np.nanmin(ylim[ll][s])

			rlim=np.interp(self.rho[ll], xlim[ll], ylim[ll], right=np.nan,left=np.nan)
			mm=np.where((self.M2 > min_mass) & (self.M2 < max_mass))
			for i in range(self.x_nsteps):
				ff=np.where((self.rho[ll,i].any() < np.min(xlim[ll])) & (self.rho[ll,i].any() > np.max(xlim[ll])))
				for j in range(self.y_nsteps):
					if (self.M2[j] > min_mass and self.M2[j] < max_mass):
						index=np.where(rlim[i]<self.M2[j])
						if np.size(index) > 1: det[i,j,index]=1
				if np.size(ff) > 1: det[i,j,ff]=0
				mm=np.where(self.M2 > max_mass)[0]
				if np.size(mm) is not 0: det[:,mm,:]=np.tile(det[:,mm[0]-1,:][:,np.newaxis,:], [np.size(mm),1])
			map=np.sum(det, axis=2)/self.norb
			detmap.append(map)
			self.detflag.append(det)
			end = time.time()
			hours, rem = divmod(end-start, 3600)
			minutes, seconds = divmod(rem, 60)
			if verbose is True: print(self.ID[ll], "time elapsed - {:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))

			if plot is True:
				fig = plt.figure(figsize=(8, 6))
				plt.rc('font', family='serif', size='20')
				plt.rc('text', usetex=True)
				plt.rc('xtick', labelsize='15')
				plt.rc('ytick', labelsize='15')
				ax = fig.add_axes([0.15, 0.15, 0.8, 0.7])
				ax.set_xscale('log')
				ax.set_yscale('log')
				ax.set_ylabel(" Mass (M$_{Jup}$)")
				ax.set_xlabel(" Semi major axis (au) ")

				levels=[10,20,50,70,90,95,99,100]
				norm = mcolors.Normalize(0, 100)
				cf0=ax.contourf(self.sma, self.M2, map.T*100, norm=norm, levels=np.arange(0,100,0.1), extend='neither', cmap='Blues', antialiased=False, zorder=0)
				contours = plt.contour(self.sma, self.M2, map.T*100, levels, cmap='bone', zorder=1, linewidths=1)
				CB = plt.colorbar(cf0,  extend='both', cmap='Blues', ticks=levels)
				CB.add_lines(contours)
				CB.set_ticks(levels)
				CB.ax.set_yticklabels(["{:.0f}".format(i)+"\%" for i in CB.get_ticks()]) # set ticks of your format

				plt.title(self.ID[ll])
				if savefig is True: plt.savefig(self.ID[ll]+'_detprob.png', dpi=300)
				#plt.show()
		return detmap
