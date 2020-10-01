
from Exo_DMC import *
from astropy.io import ascii

ID=['test01']
dist=([10])
# generate the syntetic planet population with the standard setup
# Inputs: 
# star_ID: list, required. ID of the target(s)
# star_dist: list, required. Target(s) distance, in pc

map=exodmc(ID, dist)

# the set_grid method allows to change range or resolution of the grid 
# parameters: 
#	- x_min: float, optional argument, lower limit for grid x axis (default = 1)
#	- x_max: float, optional argument, upper limit for grid x axis (default = 1000)
#	- nx: int, optional argument, number of steps in the grid x axis (default = 100)
#	- xlog: boolean, optional argument. If True the x axis will be uniformly spaced in log
#	- y_min: float, optional argument, lower limit for grid y axis (default = 0.5)
#	- y_max: float, optional argument, upper limit for grid y axis (default = 75)
#	- ny: int, optional argument, number of steps in the grid y axis (default = 100)
#	- ylog: boolean, optional argument. If True the y axis will be uniformly spaced in log
#	- ngen: float, optional argument, default=1000. Number of orbital elements sets to be generated for each point in the grid.
#		Note that all orbital parameters are uniformly distributed by default, except for the eccentricity.
#	- e_dist: string, optional. Desired eccentricity distribution. Can be uniform ('uni') or Gaussian ('gauss', default)
#	- e_mu: float, optional, mean of the gaussian eccentricity distribution, default is 0.0
#	- e_sigma: float, optional, sigma of the gaussian eccentricity distribution, default is 0.3

map.set_grid(x_min=1, x_max=100, logx=True)

#read in the detection limit: rho (arcsec) vs mass (mjup)
lim = ascii.read('detlim.dat')
xlim = [lim.field('sep')]
ylim = [lim.field('mlim')]

# the DImode method generates the detection maps and the plots 
# parameters:
#	- xlim: list, required. Detection limit(s) projected separations (arcsec)
#	- ylim: list, required. Minimum detectable mass (in Mjup) at each xlim.
#	- lxunit: string, optional. Unit for xlim, default is arcseconds ('as'), can also be set to 'au' or 'mas'
# - verbose: if True (default), the code provide the runtime for each target 
# - plot: if True (default) a .png file with the detection probability map is produced for each target

prob = map.DImode(xlim, ylim)


