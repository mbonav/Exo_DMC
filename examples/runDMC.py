
from Exo_DMC import *
from astropy.io import ascii

ID=['test01','test02']
dist=([10, 10])
# generate the syntetic planet population with the standard setup
# x_min=0.1, x_max=1000., nx=100,  logx=False, y_min=0.1, y_max=100., ny=100, logy=False, ngen=1000, e_dist='gauss', mu=0, sigma=0.3):
map=exodmc(ID, dist)

# change the grid 
map.set_grid(x_min=1, x_max=100, logx=True)

#read in the detection limit: rho (arcsec) vs mass (mjup)
lim = ascii.read('detlim.txt')
xlim = [lim.field('col1')]
ylim = [lim.field('col2')]

# get the detection maps and the plots 
prob = map.DImode(xlim, ylim, verbose=True, plot=True)
