#! /usr/bin/env python
from mcd import mcd
import numpy as np

##122 127

ph = mcd()

ph.lat = 68.22
ph.lon = 234.25
ph.xdate = 122. 
ph.xz = 5000.
ph.loct = 3.
ph.zkey = 3

yy = []
for iii in np.arange(-50,50,1):
  ph.perturkey = 4 # large-scale + small-scale
  ph.seedin = iii  #random number generator seed (unused if perturkey=0)
  ph.update()
  yy.append(ph.temp)

print np.mean(yy)
print np.std(yy)
