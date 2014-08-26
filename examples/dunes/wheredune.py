#! /usr/bin/env python
import ppclass
import numpy as np

## POINTS
tablat,tablon = np.loadtxt("dune_lat_long.txt",unpack=True)
nnn = tablat.size
forplot = open('add_text.txt', 'w')
for iii in range(nnn):
    forplot.write( "%+07.2f ; %+07.2f ; * ; r \n" % (tablon[iii],tablat[iii]) )
forplot.close()

## MAP
pl = ppclass.pp()
pl.file = "/home/aymeric/Science/MODELES/MESOSCALE/LMD_MM_MARS/WPS_GEOG/surface_new.nc"
pl.var = "z0"
pl.proj = "robin"
pl.colorbar = "Greys"
pl.vmin = 1.e-2
pl.vmax = 2. 
pl.showcb = False
pl.title = ""
pl.div = 40
pl.out = "png"
pl.filename = "z0dune"
pl.includedate = False
pl.getplot()

pl2 = ppclass.pp()
pl2 << pl
pl2.trans = 0.0
pl2.back = "vishires"
pl2.filename = "wheredune"
pl2.getplot()
