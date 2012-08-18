#! /usr/bin/env python

from mcd import mcd

gale = mcd()


gale.lat = -4.6
gale.lon = 137.4
gale.loct = 15.
gale.xz = 1.

gale.xdate = 150.6

gale.update()
gale.printmcd()
gale.printallextvar()
#gale.printextvar("rho")

#gale.seasonal()
#gale.plot1d(["tsurf","u","v"])

#gale.xdate = 270.
#gale.diurnal()
#gale.plot1d(["u","v"])

#gale.latlon()
#gale.map2d("tsurf")

import matplotlib.pyplot as mpl
mpl.show()
