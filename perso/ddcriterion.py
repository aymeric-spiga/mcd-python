#! /usr/bin/env python
from mcd import mcd
import numpy as np
import ppplot

dd = mcd()
dd.toversion5()
dd.fixedlt = True
dd.loct = 13.
dd.zkey = 3
dd.xz = 20.

for ls in [45,90,135,180,225,270]:

    print ls
    dd.xdate = ls
    dd.latlon(ndx=64,ndy=48)
    
    maxw = dd.extvartab[:,:,47]
    z0 = dd.extvartab[:,:,29]
    stress = dd.extvartab[:,:,51]
    rho = dd.extvartab[:,:,92]
    hgt = dd.extvartab[:,:,4]

    tsurf = dd.extvartab[:,:,15]
    temp = dd.extvartab[:,:,93]
   
 
    ustar = np.sqrt(stress/rho)
    wstar = maxw / 2.75
    
    limsmall = 0.5
    ustar[np.where(ustar<limsmall)] = limsmall
    
    ratio = wstar / ustar

    ratio = ((tsurf-temp)/temp)*(wstar/ustar)
    ratio[np.where(ratio<0.)] = 0.

    #ratio = (tsurf-temp)/temp
    ratio = tsurf-temp

    pl = ppplot.plot2d()    
    pl.f = ratio
    #pl.c = hgt
    pl.x = dd.xcoord
    pl.y = dd.ycoord
    #pl.vmin = 0.
    #pl.vmax = 3.
    pl.back = "ddmap"
    pl.proj = "cyl"
    pl.trans = 0.5
    pl.colorbar = "hot"
    pl.leftcorrect = True
    pl.make()
    ppplot.save(mode="png",filename="test"+str(ls))
    ppplot.close()
