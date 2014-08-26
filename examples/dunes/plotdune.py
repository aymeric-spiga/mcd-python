#! /usr/bin/env python
import ppplot
import numpy as np

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

# make a simple plot
pl = ppplot.plot1d()
pl.x = range(0,24,2)
pl.xlabel = "Local time (Martian hours)"
pl.ylabel = "Wind speed (m s$^{-1}$)"
pl.out = "png"
pl.filename = "plotdune"

pl.f = u[0,:]
pl.legend = "WE component"
pl.make()

pl.f = v[0,:]
pl.legend = "SN component"
pl.make()

pl.f = wind[0,:]
pl.legend = "total wind"
pl.make()

# save plot
ppplot.save(mode="png",filename="plotdune",includedate=False)
