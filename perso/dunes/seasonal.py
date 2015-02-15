#! /usr/bin/env python
import ppplot
import numpy as np


dalon = +154.94
dalt = 6


# load data
data = np.loadtxt("u.txt")
lat = data[:,0]
lon = data[:,1]
ls = data[:,2]
u = data[:,3:15]
data = np.loadtxt("v.txt")
v = data[:,3:15]

# computations
wind = u**2 + v**2
wind = np.sqrt(wind)

# get evolution over season
ind = np.where(lon == dalon)

# make a simple plot
pl = ppplot.plot1d()
pl.x = range(0,360,30)
pl.xlabel = "L$_s$ ($^{\circ}$)"
pl.ylabel = "Wind speed (m s$^{-1}$)"
pl.out = "png"
pl.filename = "lsdune"

pl.f = u[ind][:,dalt]
pl.legend = "WE component"
pl.make()

pl.f = v[ind][:,dalt]
pl.legend = "SN component"
pl.make()

pl.f = wind[ind][:,dalt]
pl.legend = "total wind"
pl.make()

# save plot
ppplot.save(mode="png",filename="seasonal",includedate=False)
