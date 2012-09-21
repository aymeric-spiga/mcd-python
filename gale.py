#! /usr/bin/env python

#import sys
#
#class NullDevice():
#    def write(self, s):
#        pass
#
#original_stdout = sys.stdout  # keep a reference to STDOUT
#
#sys.stdout = NullDevice()  # redirect the real STDOUT
#sys.stderr = NullDevice()  # redirect the real STDOUT

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

gale.seasonal()
gale.plot1d(["tsurf","u","v"])

import matplotlib.pyplot as mpl
mpl.savefig("temp.png",dpi=85,bbox_inches='tight',pad_inches=0.25)


#gale.xdate = 270.
#gale.diurnal()
#gale.plot1d(["u","v"])

gale.latlon()
gale.map2d("tsurf")

import matplotlib.pyplot as mpl
#mpl.show()
mpl.savefig("img/temp.png",dpi=110,bbox_inches='tight',pad_inches=0.4)
