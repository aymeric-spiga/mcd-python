#! /usr/bin/env python
import numpy as np # numerics
from mcd import mcd # MCD python wrapper
import ppplot # for plots
import ppcompute # for derivatives
import planets # for N2 computation

##############################
# this script computes       #
# GW saturation ratio        #
# from MCD field             #
# see Spiga et al. GRL 2012  #
##############################
Ls = 30.
Ls = 90.
Ls = 270. 
##############

# the constant in the saturation ratio
# ... see Spiga et al. GRL 2012
lambdah = 150e3
fzero = 7.5e-7
alpha = fzero*lambdah/(2.*np.pi)

# prepare MCD request
rq = mcd()
rq.xdate = Ls
   # xzs should not be 0 otherwise divide by u=0
rq.xzs,rq.xze = 40000.,180000. 

# request zonal mean
#rq.zonalmean()
rq.lon = 0. ; rq.secalt()

# coordinates
z = rq.ycoord
lat = rq.xcoord
# zonal wind
u = np.transpose(rq.zonwindtab)
# density
rho = np.transpose(rq.denstab)
# static stability
temp = np.transpose(rq.temptab)
dummy,dTdz = ppcompute.deriv2d(temp,lat,z)
bv = planets.Mars.N2(T0=temp,dTdz=dTdz)

# compute saturation ratio
s = np.sqrt(alpha*bv/(rho*np.abs(u)*u*u))
w = np.where(s>1) ; s[w] = 1. # S ratio is <= 1

# plot bounds -- same as Spiga et al. 2012
lb,ub = -2.,-0.5
w = np.where(np.log10(s)<lb) ; s[w] = 10.**(lb)
w = np.where(np.log10(s)>ub) ; s[w] = 10.**(ub)

# plot
pl = ppplot.plot2d()
pl.f = np.log10(s)
pl.x = lat ; pl.xlabel = "Latitude"
pl.y = z/1000. ; pl.ylabel = "Altitude (km)"
pl.c = u
pl.colorbar = "hot_r"
pl.fmt = "%.1f"
pl.units = r'$\log_{10}\mathcal{S}$'
pl.makesave(mode="png",filename="saturatio_%.0f" % (Ls))
