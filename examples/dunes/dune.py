#! /usr/bin/env python
from mcd import mcd
import numpy as np
import time as timelib

time0 = timelib.time()

test = mcd()
test.xz = 100.

tablat,tablon = np.loadtxt("dune_lat_long.txt",unpack=True)
nnn = tablat.size

ufile = open('u.txt', 'w')
vfile = open('v.txt', 'w')

for lsls in range(0,360,30):
  for iii in range(nnn):  
    test.xdate = lsls
    test.lat, test.lon = tablat[iii], tablon[iii]
    test.diurnal()
    latlonls = "%+07.2f %+07.2f %3.0f " % (tablat[iii], tablon[iii], lsls)
    tab = test.zonwindtab
    u = "%+06.2f "*tab.size % tuple(tab)
    tab = test.merwindtab
    v = "%+06.2f "*tab.size % tuple(tab)
    ufile.write( latlonls+u+" \n" )
    vfile.write( latlonls+v+" \n" )
  print lsls,timelib.time()-time0

ufile.close()
vfile.close()
