#! /usr/bin/env python
import ppplot
import numpy as np
import matplotlib.pyplot as mpl

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

# make a simple histogram
mpl.hist(np.ravel(wind),bins=20)
mpl.xlabel("Wind speed (m s$^{-1}$)")
mpl.ylabel("Population")
ppplot.save(mode="png",filename="histodune",includedate=False)
